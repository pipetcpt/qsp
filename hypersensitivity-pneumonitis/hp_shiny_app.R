## ============================================================
## Hypersensitivity Pneumonitis (HP) — Shiny QSP Dashboard
## ============================================================
## Interactive simulation dashboard for HP immune-fibrotic model
## 7 Treatment scenarios | 6 Analysis tabs | Real-time PK/PD/QSP
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)
library(DT)

## ─────────────────────────────────────────────────────────────
## MODEL CODE (inline for Shiny portability)
## ─────────────────────────────────────────────────────────────
hp_code <- '
$PARAM
AG_exposure=1.0 AG_avoidance=0.0 k_ag_clear=0.5
k_mac_act=0.8 k_mac_base=0.05 k_mac_death=0.15 k_mac_M2=0.05
k_neu_rec=0.6 k_neu_death=0.8
k_TNF_prod=0.4 k_TNF_deg=1.2 k_IL6_prod=0.35 k_IL6_deg=1.5
k_IL12_prod=0.3 k_IL12_deg=1.0 k_IFNg_prod=0.5 k_IFNg_deg=1.2
k_TGFb_prod=0.25 k_TGFb_deg=0.8 k_IL17_prod=0.3 k_IL17_deg=1.0
k_IL10_prod=0.15 k_IL10_deg=0.8
k_Th1_diff=0.4 k_Th17_diff=0.25 k_Treg_diff=0.15
k_Th1_death=0.1 k_Th17_death=0.12 k_Treg_death=0.08
T_naive_0=1.0 Treg_inhib=0.3
k_gran_form=0.2 k_gran_res=0.15 k_gran_fib=0.08
k_fib_act=0.3 k_fib_death=0.2 k_myo_diff=0.25 k_myo_death=0.15
k_col_prod=0.15 k_col_deg=0.05 k_ROS_prod=0.3 k_ROS_deg=0.5
FVC_0=90.0 DLCO_0=85.0
k_FVC_col=0.08 k_DLCO_col=0.10 k_FVC_recov=0.01 k_DLCO_recov=0.005
ka_PDN=1.5 CL_PDN=10.0 Vd_PDN=35.0 F_PDN=0.82
ka_MMF=2.0 CL_MPA=16.0 Vd_MPA=55.0 F_MMF=0.94
ka_Nint=0.8 CL_Nint=85.0 Vd_Nint=1050.0 F_Nint=0.047
EC50_PDN_inflam=150.0 EC50_MPA_prolif=0.5 EC50_Nint_fib=200.0
Emax_PDN=0.85 Emax_MPA=0.80 Emax_Nint=0.75
KL6_base=250.0 k_KL6_col=2.0 k_KL6_deg=0.3

$CMT
AG_lung M_M1 M_M2 Neutrophil
C_TNF C_IL6 C_IL12 C_IFNg C_TGFb C_IL17 C_IL10
T_Th1 T_Th17 T_Treg Granuloma
Fibroblast Myofib Collagen ROS
PDN_gut PDN_cent MPA_gut MPA_cent Nint_gut Nint_cent
FVC DLCO KL6_serum

$INIT
AG_lung=0 M_M1=0.1 M_M2=0.05 Neutrophil=0.05
C_TNF=0 C_IL6=0 C_IL12=0 C_IFNg=0 C_TGFb=0 C_IL17=0 C_IL10=0
T_Th1=0.1 T_Th17=0.05 T_Treg=0.1 Granuloma=0
Fibroblast=0.05 Myofib=0 Collagen=0 ROS=0
PDN_gut=0 PDN_cent=0 MPA_gut=0 MPA_cent=0 Nint_gut=0 Nint_cent=0
FVC=90 DLCO=85 KL6_serum=250

$ODE
double PDN_conc = PDN_cent / Vd_PDN * 1000;
double MPA_conc = MPA_cent / Vd_MPA;
double Nint_conc = Nint_cent / Vd_Nint * 1000;
double E_PDN = Emax_PDN * PDN_conc / (EC50_PDN_inflam + PDN_conc);
double E_MPA = Emax_MPA * MPA_conc / (EC50_MPA_prolif + MPA_conc);
double E_Nint = Emax_Nint * Nint_conc / (EC50_Nint_fib + Nint_conc);
double AG_in = AG_exposure * (1.0 - AG_avoidance);
dxdt_AG_lung = AG_in - k_ag_clear * AG_lung;
double mac_stim = k_mac_act * AG_lung * (1 + 0.3*C_IFNg) * (1 - 0.8*E_PDN);
dxdt_M_M1 = mac_stim + k_mac_base - k_mac_death*M_M1 - k_mac_M2*C_TGFb*M_M1;
dxdt_M_M2 = k_mac_M2*C_TGFb*M_M1 - k_mac_death*M_M2;
dxdt_Neutrophil = k_neu_rec*(C_TNF+C_IL17)*(1-0.5*E_PDN) - k_neu_death*Neutrophil;
dxdt_C_TNF  = k_TNF_prod*M_M1*(1-E_PDN) - k_TNF_deg*C_TNF;
dxdt_C_IL6  = k_IL6_prod*(M_M1+T_Th17)*(1-E_PDN) - k_IL6_deg*C_IL6;
dxdt_C_IL12 = k_IL12_prod*M_M1 - k_IL12_deg*C_IL12;
dxdt_C_IFNg = k_IFNg_prod*T_Th1 + 0.1*M_M1 - k_IFNg_deg*C_IFNg;
dxdt_C_TGFb = k_TGFb_prod*(M_M2+Myofib) - k_TGFb_deg*C_TGFb;
dxdt_C_IL17 = k_IL17_prod*T_Th17 - k_IL17_deg*C_IL17;
dxdt_C_IL10 = k_IL10_prod*(T_Treg+M_M2) - k_IL10_deg*C_IL10;
double Th1_diff = k_Th1_diff*C_IL12*T_naive_0*(1-Treg_inhib*T_Treg)*(1-E_PDN)*(1-E_MPA);
dxdt_T_Th1  = Th1_diff - k_Th1_death*T_Th1;
double Th17_diff = k_Th17_diff*C_IL6*C_TGFb*T_naive_0*(1-Treg_inhib*T_Treg)*(1-E_MPA);
dxdt_T_Th17 = Th17_diff - k_Th17_death*T_Th17;
dxdt_T_Treg = k_Treg_diff*C_TGFb*C_IL10*T_naive_0 - k_Treg_death*T_Treg;
double gran_form = k_gran_form*M_M1*T_Th1;
double gran_res  = k_gran_res*C_IL10*Granuloma;
double gran_fib  = k_gran_fib*Granuloma*(1-C_IL10)*(1-E_PDN);
dxdt_Granuloma = gran_form - gran_res - gran_fib;
dxdt_Fibroblast = k_fib_act*C_TGFb*(1+0.3*gran_fib)*(1-0.5*E_PDN) - k_fib_death*Fibroblast;
dxdt_Myofib = k_myo_diff*C_TGFb*Fibroblast*(1-E_Nint) - k_myo_death*Myofib;
dxdt_Collagen = k_col_prod*Myofib*(1-0.6*E_Nint) - k_col_deg*Collagen;
dxdt_ROS = k_ROS_prod*(M_M2+Neutrophil) - k_ROS_deg*ROS;
dxdt_PDN_gut  = -ka_PDN*PDN_gut;
dxdt_PDN_cent = ka_PDN*F_PDN*PDN_gut - (CL_PDN/Vd_PDN)*PDN_cent;
dxdt_MPA_gut  = -ka_MMF*MPA_gut;
dxdt_MPA_cent = ka_MMF*F_MMF*MPA_gut - (CL_MPA/Vd_MPA)*MPA_cent;
dxdt_Nint_gut  = -ka_Nint*Nint_gut;
dxdt_Nint_cent = ka_Nint*F_Nint*Nint_gut - (CL_Nint/Vd_Nint)*Nint_cent;
dxdt_FVC  = -k_FVC_col*Collagen*FVC/100.0 + k_FVC_recov*(FVC_0-FVC);
dxdt_DLCO = -k_DLCO_col*Collagen*DLCO/100.0 + k_DLCO_recov*(DLCO_0-DLCO);
dxdt_KL6_serum = k_KL6_col*(Collagen+ROS+M_M2) - k_KL6_deg*(KL6_serum-KL6_base);

$TABLE
double PDN_Cconc = PDN_cent/Vd_PDN*1000;
double MPA_Cconc = MPA_cent/Vd_MPA;
double Nint_Cconc = Nint_cent/Vd_Nint*1000;
double Inflam = C_TNF + C_IL6 + C_IFNg;
double Fibrosis_idx = Collagen + Myofib;

$CAPTURE
PDN_Cconc MPA_Cconc Nint_Cconc
FVC DLCO KL6_serum Inflam Fibrosis_idx
M_M1 M_M2 T_Th1 T_Th17 T_Treg
C_TNF C_IL6 C_TGFb C_IFNg C_IL10
Granuloma Collagen Myofib ROS Neutrophil
'

## ─────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(
    title = HTML("HP-QSP Model"),
    titleWidth = 250
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",  tabName = "tab_patient",  icon = icon("user-md")),
      menuItem("PK Profiles",      tabName = "tab_pk",       icon = icon("chart-line")),
      menuItem("Immune Response",  tabName = "tab_immune",   icon = icon("shield-alt")),
      menuItem("Fibrosis & QSP",   tabName = "tab_fibrosis", icon = icon("lungs")),
      menuItem("Clinical Endpoints", tabName = "tab_clinical", icon = icon("stethoscope")),
      menuItem("Scenario Comparison", tabName = "tab_compare", icon = icon("exchange-alt")),
      menuItem("Biomarkers",       tabName = "tab_biomarker",icon = icon("flask")),
      menuItem("References",       tabName = "tab_refs",     icon = icon("book"))
    ),
    hr(),
    ## ── PATIENT PROFILE ──────────────────────────────────
    h5("Patient Parameters", style = "padding-left:15px; color:#ECF0F1"),
    sliderInput("pt_fvc",  "Baseline FVC (%)",   55, 100, 82, 1),
    sliderInput("pt_dlco", "Baseline DLCO (%)",  40,  95, 72, 1),
    hr(),
    ## ── ANTIGEN EXPOSURE ──────────────────────────────────
    h5("Antigen Exposure", style = "padding-left:15px; color:#ECF0F1"),
    sliderInput("ag_exp",  "Antigen Level (0-2)", 0, 2, 1.0, 0.1),
    sliderInput("ag_avoid","Avoidance Efficacy",  0, 1, 0.0, 0.05),
    hr(),
    ## ── TREATMENT SELECTION ───────────────────────────────
    h5("Drug Therapy", style = "padding-left:15px; color:#ECF0F1"),
    checkboxInput("use_pdn",  "Prednisolone",  FALSE),
    conditionalPanel(
      condition = "input.use_pdn == true",
      sliderInput("pdn_dose", "PDN dose (mg/d)", 5, 80, 40, 5)
    ),
    checkboxInput("use_mmf",  "MMF",  FALSE),
    conditionalPanel(
      condition = "input.use_mmf == true",
      sliderInput("mmf_dose", "MMF dose (mg BID)", 500, 2000, 1500, 250)
    ),
    checkboxInput("use_nint", "Nintedanib", FALSE),
    conditionalPanel(
      condition = "input.use_nint == true",
      sliderInput("nint_dose","NINT dose (mg BID)", 50, 200, 150, 50)
    ),
    hr(),
    sliderInput("sim_days", "Simulation (days)", 90, 1095, 730, 30),
    actionButton("run_sim", "Run Simulation", class = "btn-success btn-block",
                 icon = icon("play"))
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box { border-radius: 8px; }
      .value-box { border-radius: 8px; }
    "))),

    tabItems(
      ## ─────────────────────────────────────────────────────
      ## TAB 1: PATIENT PROFILE
      ## ─────────────────────────────────────────────────────
      tabItem("tab_patient",
        fluidRow(
          valueBoxOutput("vb_fvc",   width = 3),
          valueBoxOutput("vb_dlco",  width = 3),
          valueBoxOutput("vb_kl6",   width = 3),
          valueBoxOutput("vb_inflam",width = 3)
        ),
        fluidRow(
          box(width = 12, title = "Disease Overview: Hypersensitivity Pneumonitis",
              status = "primary", solidHeader = TRUE,
              p(strong("Hypersensitivity Pneumonitis (HP)"), "is an immune-mediated interstitial
              lung disease caused by repeated inhalation of organic antigens (avian proteins,
              thermophilic bacteria, molds). The disease spectrum ranges from acute/subacute
              (reversible with antigen avoidance) to chronic/fibrotic HP (progressive, poor prognosis)."),
              p("Key pathogenic features include: Th1-dominant CD4+ lymphocytic alveolitis,
              non-caseating granuloma formation, TGF-β-driven myofibroblast differentiation,
              and progressive collagen deposition (UIP or NSIP-like pattern on HRCT)."),
              tags$ul(
                tags$li(strong("Diagnostic criteria:"), "antigen exposure history + BAL lymphocytosis >20%
                         + HRCT pattern + lung biopsy/BAL CD4:CD8 ratio"),
                tags$li(strong("Key biomarkers:"), "Serum KL-6 (>500 U/mL), SP-D, specific IgG precipitins,
                         BAL lymphocytosis, FVC decline"),
                tags$li(strong("Prognosis:"), "Progressive HP 5-year mortality 20–40%; UIP pattern = worse outcome"),
                tags$li(strong("Treatment pillars:"), "1) Antigen avoidance (most important)
                         2) Corticosteroids (acute/subacute) 3) Immunosuppressants (MMF/AZA)
                         4) Nintedanib (antifibrotic, approved 2022 for fibrotic ILD)")
              )
          )
        )
      ),

      ## ─────────────────────────────────────────────────────
      ## TAB 2: PK PROFILES
      ## ─────────────────────────────────────────────────────
      tabItem("tab_pk",
        fluidRow(
          box(width = 12, title = "Drug Concentration-Time Profiles",
              status = "info", solidHeader = TRUE,
              plotlyOutput("plot_pk", height = "450px"))
        ),
        fluidRow(
          box(width = 6, title = "Prednisolone PK Parameters",
              status = "warning", solidHeader = FALSE,
              tableOutput("tbl_pk_pdn")),
          box(width = 6, title = "Nintedanib PK Parameters",
              status = "warning", solidHeader = FALSE,
              tableOutput("tbl_pk_nint"))
        )
      ),

      ## ─────────────────────────────────────────────────────
      ## TAB 3: IMMUNE RESPONSE
      ## ─────────────────────────────────────────────────────
      tabItem("tab_immune",
        fluidRow(
          box(width = 6, title = "T Cell Dynamics (Th1 / Th17 / Treg)",
              status = "primary", solidHeader = TRUE,
              plotlyOutput("plot_tcell", height = "380px")),
          box(width = 6, title = "Macrophage Polarization (M1/M2)",
              status = "primary", solidHeader = TRUE,
              plotlyOutput("plot_mac", height = "380px"))
        ),
        fluidRow(
          box(width = 6, title = "Pro-inflammatory Cytokines",
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plot_cytokine_pro", height = "340px")),
          box(width = 6, title = "Regulatory / Fibrotic Mediators",
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_cytokine_reg", height = "340px"))
        )
      ),

      ## ─────────────────────────────────────────────────────
      ## TAB 4: FIBROSIS & QSP
      ## ─────────────────────────────────────────────────────
      tabItem("tab_fibrosis",
        fluidRow(
          box(width = 6, title = "Granuloma Dynamics",
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_gran", height = "340px")),
          box(width = 6, title = "Fibrosis Cascade",
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plot_fibrosis", height = "340px"))
        ),
        fluidRow(
          box(width = 6, title = "Myofibroblast & Collagen Deposition",
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plot_collagen", height = "340px")),
          box(width = 6, title = "Oxidative Stress (ROS)",
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_ros", height = "340px"))
        )
      ),

      ## ─────────────────────────────────────────────────────
      ## TAB 5: CLINICAL ENDPOINTS
      ## ─────────────────────────────────────────────────────
      tabItem("tab_clinical",
        fluidRow(
          box(width = 6, title = "FVC (% Predicted) — Lung Function",
              status = "success", solidHeader = TRUE,
              plotlyOutput("plot_fvc", height = "380px")),
          box(width = 6, title = "DLCO (% Predicted) — Gas Transfer",
              status = "success", solidHeader = TRUE,
              plotlyOutput("plot_dlco", height = "380px"))
        ),
        fluidRow(
          box(width = 12, title = "Key Clinical Milestones",
              status = "info", solidHeader = TRUE,
              DTOutput("tbl_milestones"))
        )
      ),

      ## ─────────────────────────────────────────────────────
      ## TAB 6: SCENARIO COMPARISON
      ## ─────────────────────────────────────────────────────
      tabItem("tab_compare",
        fluidRow(
          box(width = 12, title = "Multi-Scenario FVC Comparison (7 Treatment Arms)",
              status = "primary", solidHeader = TRUE,
              plotlyOutput("plot_compare_fvc", height = "420px"))
        ),
        fluidRow(
          box(width = 6, title = "Collagen Deposition Comparison",
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plot_compare_col", height = "340px")),
          box(width = 6, title = "Inflammation Index Comparison",
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_compare_inflam", height = "340px"))
        )
      ),

      ## ─────────────────────────────────────────────────────
      ## TAB 7: BIOMARKERS
      ## ─────────────────────────────────────────────────────
      tabItem("tab_biomarker",
        fluidRow(
          box(width = 6, title = "Serum KL-6 (Fibrosis Biomarker)",
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_kl6", height = "350px")),
          box(width = 6, title = "Composite Inflammation Index",
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plot_inflam", height = "350px"))
        ),
        fluidRow(
          box(width = 12, title = "Biomarker Summary Table",
              status = "info", solidHeader = TRUE,
              DTOutput("tbl_biomarker"))
        )
      ),

      ## ─────────────────────────────────────────────────────
      ## TAB 8: REFERENCES
      ## ─────────────────────────────────────────────────────
      tabItem("tab_refs",
        fluidRow(
          box(width = 12, title = "Key References",
              status = "info", solidHeader = TRUE,
              tags$ul(
                tags$li("Morisset J, et al. (2020) Mycophenolate or Azathioprine for HP. Lancet Respir Med."),
                tags$li("Raghu G, et al. (2021) Nintedanib in fibrotic ILD. NEJM."),
                tags$li("Giménez A, et al. (2018) Prognosis of fibrotic HP. ERJ."),
                tags$li("Walsh SLF, et al. (2014) Natural history of HP. Thorax."),
                tags$li("Fernández Pérez ER, et al. (2018) HP diagnosis/management. ATS."),
                tags$li("Vasakova M, et al. (2017) HP guidelines. Respirology.")
              )
          )
        )
      )
    ) # tabItems
  ) # dashboardBody
) # dashboardPage

## ─────────────────────────────────────────────────────────────
## SERVER
## ─────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## Build model once
  mod <- reactive({ mcode("hp_qsp_shiny", hp_code, quiet = TRUE) })

  ## ── RUN SIMULATION ────────────────────────────────────────
  sim_data <- eventReactive(input$run_sim, {
    m <- mod()
    m <- param(m,
      AG_exposure  = input$ag_exp,
      AG_avoidance = input$ag_avoid,
      FVC_0  = input$pt_fvc,
      DLCO_0 = input$pt_dlco
    )
    m <- init(m, FVC = input$pt_fvc, DLCO = input$pt_dlco)

    # Build dosing events
    dose_ev <- NULL

    if (input$use_pdn) {
      e_pdn <- ev(amt = input$pdn_dose, cmt = "PDN_gut",
                  ii = 24, addl = input$sim_days - 1, time = 0)
      dose_ev <- if (is.null(dose_ev)) e_pdn else ev(dose_ev, e_pdn)
    }

    if (input$use_mmf) {
      e_mmf <- ev(amt = input$mmf_dose, cmt = "MPA_gut",
                  ii = 12, addl = input$sim_days * 2 - 1, time = 0)
      dose_ev <- if (is.null(dose_ev)) e_mmf else ev(dose_ev, e_mmf)
    }

    if (input$use_nint) {
      e_nint <- ev(amt = input$nint_dose, cmt = "Nint_gut",
                   ii = 12, addl = input$sim_days * 2 - 1, time = 0)
      dose_ev <- if (is.null(dose_ev)) e_nint else ev(dose_ev, e_nint)
    }

    if (!is.null(dose_ev)) {
      mrgsim(m, ev = dose_ev, end = input$sim_days, delta = 1) %>% as.data.frame()
    } else {
      mrgsim(m, end = input$sim_days, delta = 1) %>% as.data.frame()
    }
  }, ignoreNULL = FALSE)

  ## ── VALUE BOXES ───────────────────────────────────────────
  output$vb_fvc <- renderValueBox({
    d <- sim_data()
    last_fvc <- round(tail(d$FVC, 1), 1)
    color <- if (last_fvc >= 70) "green" else if (last_fvc >= 55) "yellow" else "red"
    valueBox(paste0(last_fvc, "%"), "Final FVC", icon = icon("lungs"), color = color)
  })

  output$vb_dlco <- renderValueBox({
    d <- sim_data()
    last_dlco <- round(tail(d$DLCO, 1), 1)
    color <- if (last_dlco >= 60) "green" else if (last_dlco >= 45) "yellow" else "red"
    valueBox(paste0(last_dlco, "%"), "Final DLCO", icon = icon("heartbeat"), color = color)
  })

  output$vb_kl6 <- renderValueBox({
    d <- sim_data()
    last_kl6 <- round(tail(d$KL6_serum, 1), 0)
    color <- if (last_kl6 < 500) "green" else if (last_kl6 < 1000) "yellow" else "red"
    valueBox(paste0(last_kl6, " U/mL"), "Final KL-6", icon = icon("flask"), color = color)
  })

  output$vb_inflam <- renderValueBox({
    d <- sim_data()
    last_inf <- round(tail(d$Inflam, 1), 3)
    color <- if (last_inf < 0.3) "green" else if (last_inf < 0.8) "yellow" else "red"
    valueBox(round(last_inf, 2), "Inflammation\nIndex", icon = icon("fire"), color = color)
  })

  ## ── PK PLOT ───────────────────────────────────────────────
  output$plot_pk <- renderPlotly({
    d <- sim_data()
    df <- data.frame(
      time = d$time,
      Prednisolone_ngmL = d$PDN_Cconc,
      MPA_ugmL = d$MPA_Cconc,
      Nintedanib_ngmL = d$Nint_Cconc
    ) %>%
      pivot_longer(-time, names_to = "Drug", values_to = "Concentration")

    p <- ggplot(df, aes(x = time, y = Concentration, color = Drug)) +
      geom_line(linewidth = 1) +
      facet_wrap(~Drug, scales = "free_y") +
      labs(title = "Drug Concentration-Time Profiles",
           x = "Days", y = "Concentration") +
      theme_bw(base_size = 11) + theme(legend.position = "none")
    ggplotly(p) %>% layout(margin = list(l=60, r=20, b=60, t=60))
  })

  output$tbl_pk_pdn <- renderTable({
    data.frame(
      Parameter = c("Bioavailability (F)", "ka", "CL", "Vd", "t½", "EC50 (anti-inflam)"),
      Value = c("82%", "1.5 /h", "10 L/h", "35 L", "~3.4 h", "150 ng/mL"),
      Source = c("Literature", "Pop PK", "Pop PK", "Pop PK", "Calculated", "Estimated")
    )
  })

  output$tbl_pk_nint <- renderTable({
    data.frame(
      Parameter = c("Bioavailability (F)", "ka", "CL", "Vd", "t½", "EC50 (antifibrotic)"),
      Value = c("4.7%", "0.8 /h", "85 L/h", "1050 L", "~10.8 h", "200 ng/mL"),
      Source = c("FDA label", "Pop PK", "Pop PK", "Pop PK", "Calculated", "Estimated")
    )
  })

  ## ── IMMUNE: T CELL PLOT ───────────────────────────────────
  output$plot_tcell <- renderPlotly({
    d <- sim_data()
    df <- data.frame(time=d$time, Th1=d$T_Th1, Th17=d$T_Th17, Treg=d$T_Treg) %>%
      pivot_longer(-time, names_to="CellType", values_to="Level")
    p <- ggplot(df, aes(x=time, y=Level, color=CellType)) + geom_line(linewidth=1.1) +
      scale_color_manual(values=c("Th1"="#E74C3C","Th17"="#E67E22","Treg"="#27AE60")) +
      labs(title="T Cell Subsets", x="Days", y="Level (AU)") + theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plot_mac <- renderPlotly({
    d <- sim_data()
    df <- data.frame(time=d$time, M1=d$M_M1, M2=d$M_M2) %>%
      pivot_longer(-time, names_to="Type", values_to="Level")
    p <- ggplot(df, aes(x=time, y=Level, color=Type)) + geom_line(linewidth=1.1) +
      scale_color_manual(values=c("M1"="#E74C3C","M2"="#8E44AD")) +
      labs(title="Macrophage Polarization", x="Days", y="Level (AU)") + theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plot_cytokine_pro <- renderPlotly({
    d <- sim_data()
    df <- data.frame(time=d$time, `TNF-a`=d$C_TNF, `IL-6`=d$C_IL6,
                     `IFN-g`=d$C_IFNg, `IL-17`=d$C_IL17) %>%
      pivot_longer(-time, names_to="Cytokine", values_to="Level")
    p <- ggplot(df, aes(x=time, y=Level, color=Cytokine)) + geom_line(linewidth=1) +
      labs(title="Pro-inflammatory Cytokines", x="Days", y="Level (AU)") + theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plot_cytokine_reg <- renderPlotly({
    d <- sim_data()
    df <- data.frame(time=d$time, `TGF-b1`=d$C_TGFb, `IL-10`=d$C_IL10) %>%
      pivot_longer(-time, names_to="Cytokine", values_to="Level")
    p <- ggplot(df, aes(x=time, y=Level, color=Cytokine)) + geom_line(linewidth=1) +
      scale_color_manual(values=c("TGF-b1"="#8E44AD","IL-10"="#27AE60")) +
      labs(title="Regulatory / Fibrotic Mediators", x="Days", y="Level (AU)") + theme_bw(base_size=11)
    ggplotly(p)
  })

  ## ── FIBROSIS PLOTS ────────────────────────────────────────
  output$plot_gran <- renderPlotly({
    d <- sim_data()
    p <- ggplot(data.frame(time=d$time, gran=d$Granuloma), aes(x=time, y=gran)) +
      geom_line(linewidth=1.2, color="#CA6F1E") + geom_area(fill="#CA6F1E", alpha=0.3) +
      labs(title="Granuloma Burden", x="Days", y="Granuloma (AU)") + theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plot_fibrosis <- renderPlotly({
    d <- sim_data()
    df <- data.frame(time=d$time, Fibroblast=d$Fibroblast, Myofib=d$Myofib) %>%
      pivot_longer(-time, names_to="Cell", values_to="Level")
    p <- ggplot(df, aes(x=time, y=Level, color=Cell)) + geom_line(linewidth=1.1) +
      scale_color_manual(values=c("Fibroblast"="#F1948A","Myofib"="#E74C3C")) +
      labs(title="Fibroblast → Myofibroblast", x="Days", y="Level (AU)") + theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plot_collagen <- renderPlotly({
    d <- sim_data()
    p <- ggplot(data.frame(time=d$time, col=d$Collagen), aes(x=time, y=col)) +
      geom_line(linewidth=1.3, color="#922B21") + geom_area(fill="#922B21", alpha=0.3) +
      labs(title="Collagen Deposition (Fibrosis Burden)", x="Days", y="Collagen (AU)") + theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plot_ros <- renderPlotly({
    d <- sim_data()
    p <- ggplot(data.frame(time=d$time, ros=d$ROS), aes(x=time, y=ros)) +
      geom_line(linewidth=1.2, color="#F39C12") +
      labs(title="Reactive Oxygen Species", x="Days", y="ROS (AU)") + theme_bw(base_size=11)
    ggplotly(p)
  })

  ## ── CLINICAL ENDPOINT PLOTS ───────────────────────────────
  output$plot_fvc <- renderPlotly({
    d <- sim_data()
    p <- ggplot(data.frame(time=d$time, fvc=d$FVC), aes(x=time, y=fvc)) +
      geom_line(linewidth=1.5, color="#27AE60") +
      geom_hline(yintercept=70, linetype="dashed", color="red") +
      annotate("text", x=30, y=71.5, label="Transplant ref. ~70%", color="red", size=3) +
      labs(title="FVC (% Predicted)", x="Days", y="FVC (%)") +
      coord_cartesian(ylim=c(30,100)) + theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plot_dlco <- renderPlotly({
    d <- sim_data()
    p <- ggplot(data.frame(time=d$time, dlco=d$DLCO), aes(x=time, y=dlco)) +
      geom_line(linewidth=1.5, color="#2471A3") +
      geom_hline(yintercept=40, linetype="dashed", color="red") +
      labs(title="DLCO (% Predicted)", x="Days", y="DLCO (%)") +
      coord_cartesian(ylim=c(20,95)) + theme_bw(base_size=11)
    ggplotly(p)
  })

  output$tbl_milestones <- renderDT({
    d <- sim_data()
    timepoints <- c(90, 180, 365, 548, 730)
    timepoints <- timepoints[timepoints <= max(d$time)]
    milestones <- d %>%
      filter(time %in% timepoints) %>%
      select(time, FVC, DLCO, KL6_serum, Collagen, Inflam) %>%
      mutate(across(where(is.numeric), ~round(., 2))) %>%
      rename(Day=time, `FVC (%)`=FVC, `DLCO (%)`=DLCO,
             `KL-6 (U/mL)`=KL6_serum, `Collagen (AU)`=Collagen, `Inflam Index`=Inflam)
    datatable(milestones, options=list(pageLength=10), rownames=FALSE)
  })

  ## ── SCENARIO COMPARISON ───────────────────────────────────
  scenario_data <- reactive({
    m <- mod()
    scenarios <- list(
      list(label="1.No Tx", ag_av=0, pdn=0, mmf=0, nint=0),
      list(label="2.Avoidance", ag_av=0.9, pdn=0, mmf=0, nint=0),
      list(label="3.PDN 40mg", ag_av=0, pdn=40, mmf=0, nint=0),
      list(label="4.MMF 1500", ag_av=0, pdn=0, mmf=1500, nint=0),
      list(label="5.Nint 150", ag_av=0, pdn=0, mmf=0, nint=150),
      list(label="6.PDN+MMF", ag_av=0.5, pdn=40, mmf=1500, nint=0),
      list(label="7.Nint+Avoid", ag_av=0.95, pdn=0, mmf=0, nint=150)
    )
    DAYS <- min(input$sim_days, 730)
    bind_rows(lapply(scenarios, function(sc) {
      mp <- param(m, AG_exposure=1.0, AG_avoidance=sc$ag_av,
                  FVC_0=input$pt_fvc, DLCO_0=input$pt_dlco)
      mp <- init(mp, FVC=input$pt_fvc, DLCO=input$pt_dlco)
      dose_ev <- NULL
      if (sc$pdn > 0) {
        e1 <- ev(amt=sc$pdn, cmt="PDN_gut", ii=24, addl=DAYS-1)
        dose_ev <- if(is.null(dose_ev)) e1 else ev(dose_ev, e1)
      }
      if (sc$mmf > 0) {
        e2 <- ev(amt=sc$mmf, cmt="MPA_gut", ii=12, addl=DAYS*2-1)
        dose_ev <- if(is.null(dose_ev)) e2 else ev(dose_ev, e2)
      }
      if (sc$nint > 0) {
        e3 <- ev(amt=sc$nint, cmt="Nint_gut", ii=12, addl=DAYS*2-1)
        dose_ev <- if(is.null(dose_ev)) e3 else ev(dose_ev, e3)
      }
      out <- if (!is.null(dose_ev)) {
        mrgsim(mp, ev=dose_ev, end=DAYS, delta=7) %>% as.data.frame()
      } else {
        mrgsim(mp, end=DAYS, delta=7) %>% as.data.frame()
      }
      out$Scenario <- sc$label
      out
    }))
  })

  output$plot_compare_fvc <- renderPlotly({
    d <- scenario_data()
    colors7 <- c("#E74C3C","#27AE60","#3498DB","#9B59B6","#F39C12","#1ABC9C","#2C3E50")
    p <- ggplot(d, aes(x=time, y=FVC, color=Scenario)) + geom_line(linewidth=1.1) +
      scale_color_manual(values=colors7) +
      geom_hline(yintercept=70, linetype="dashed", color="red") +
      labs(title="FVC Comparison — 7 Treatment Scenarios",
           x="Days", y="FVC (% predicted)") +
      coord_cartesian(ylim=c(30,100)) + theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plot_compare_col <- renderPlotly({
    d <- scenario_data()
    colors7 <- c("#E74C3C","#27AE60","#3498DB","#9B59B6","#F39C12","#1ABC9C","#2C3E50")
    p <- ggplot(d, aes(x=time, y=Collagen, color=Scenario)) + geom_line(linewidth=1.1) +
      scale_color_manual(values=colors7) +
      labs(title="Collagen Deposition Comparison", x="Days", y="Collagen (AU)") + theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plot_compare_inflam <- renderPlotly({
    d <- scenario_data()
    colors7 <- c("#E74C3C","#27AE60","#3498DB","#9B59B6","#F39C12","#1ABC9C","#2C3E50")
    p <- ggplot(d, aes(x=time, y=Inflam, color=Scenario)) + geom_line(linewidth=1.1) +
      scale_color_manual(values=colors7) +
      labs(title="Inflammation Index Comparison", x="Days", y="Inflam Index (AU)") + theme_bw(base_size=11)
    ggplotly(p)
  })

  ## ── BIOMARKER PLOTS ───────────────────────────────────────
  output$plot_kl6 <- renderPlotly({
    d <- sim_data()
    p <- ggplot(data.frame(time=d$time, kl6=d$KL6_serum), aes(x=time, y=kl6)) +
      geom_line(linewidth=1.3, color="#F39C12") +
      geom_hline(yintercept=500, linetype="dashed", color="red") +
      annotate("text", x=30, y=520, label="Diagnostic threshold (500 U/mL)", color="red", size=3) +
      labs(title="Serum KL-6", x="Days", y="KL-6 (U/mL)") + theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plot_inflam <- renderPlotly({
    d <- sim_data()
    p <- ggplot(data.frame(time=d$time, inf=d$Inflam), aes(x=time, y=inf)) +
      geom_line(linewidth=1.3, color="#E74C3C") +
      geom_area(fill="#E74C3C", alpha=0.25) +
      labs(title="Composite Inflammation Index\n(TNF-α + IL-6 + IFN-γ)",
           x="Days", y="Index (AU)") + theme_bw(base_size=11)
    ggplotly(p)
  })

  output$tbl_biomarker <- renderDT({
    d <- sim_data()
    timepoints <- c(0, 90, 180, 365, 730)
    timepoints <- timepoints[timepoints <= max(d$time)]
    df <- d %>%
      filter(time %in% timepoints) %>%
      select(time, KL6_serum, Inflam, Collagen, Myofib, Granuloma, ROS) %>%
      mutate(across(where(is.numeric), ~round(., 3))) %>%
      rename(Day=time, `KL-6 (U/mL)`=KL6_serum, `Inflam Index`=Inflam,
             `Collagen (AU)`=Collagen, `Myofib (AU)`=Myofib,
             `Granuloma (AU)`=Granuloma, `ROS (AU)`=ROS)
    datatable(df, options=list(pageLength=10), rownames=FALSE)
  })
}

## ─────────────────────────────────────────────────────────────
## RUN
## ─────────────────────────────────────────────────────────────
shinyApp(ui, server)
