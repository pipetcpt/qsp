## =============================================================================
## ALS QSP Shiny Dashboard
## Amyotrophic Lateral Sclerosis — Interactive Simulation
## Tabs: Patient Profile · Drug PK · Biomarkers · Clinical Endpoints ·
##       Scenario Comparison · Mechanistic Pathways
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)

## ── Inline mrgsolve model ────────────────────────────────────────────────────
als_code <- '
$PARAM
k_MN_death=0.00080 k_MN_death_U=0.00060 MN_upper_0=1.0 MN_lower_0=1.0
k_SOD1_syn=0.50  k_SOD1_deg=0.10 k_SOD1_mis=0.08 k_SOD1_clr=0.015 f_SOD1_mut=0.0
k_TDP_export=0.020 k_TDP_import=0.150 k_TDP_agg=0.010 k_TDP_clr=0.008 TDP43_nuc_0=10.0
Glu_base=1.0 k_Glu_rel=2.0 k_EAAT2=1.80 k_EAAT2_ALS=0.40
k_Ca_entry=0.50 k_Ca_efflux=1.20 Ca_i_0=0.10
k_ROS_prod=0.40 k_ROS_scav=0.60 k_GSH_syn=0.50 k_GSH_cons=0.12 ROS_0=0.50 GSH_0=4.0
k_Mito_dam=0.06 k_Mito_rep=0.03 Mito_0=1.0
k_Mic_act=0.12 k_Mic_res=0.060 k_TNFa_prod=0.60 k_TNFa_deg=0.35
k_IL1b_prod=0.40 k_IL1b_deg=0.25 Mic_0=0.10
k_BDNF_prod=0.25 k_BDNF_deg=0.18 BDNF_0=1.0 BDNF_prot_wt=0.30
k_NfL_rel=0.10 k_NfL_clr=0.025 NfL_CSF_0=5.0
ALSFRS_0=48.0 FVC_0=100.0 k_FVC_dec=0.0060
F_RIL=0.60 ka_RIL=0.80 CL_RIL=28.0 V1_RIL=245.0 Q_RIL=15.0 V2_RIL=112.0
IC50_RIL=0.50 Emax_RIL=0.60
CL_EDA=18.0 V_EDA=120.0 IC50_EDA=1.20 Emax_EDA=0.70
ka_TOF=0.030 CL_TOF=0.50 V1_TOF=15.0 Q_TOF=0.30 V2_TOF=5.0
EC50_TOF=0.10 Emax_TOF=0.80
F_PB=0.85 ka_PB=1.20 CL_PB=12.0 V_PB=50.0 EC50_PB=50.0 Emax_PB=0.50 Emax_mito_PB=0.30

$CMT DEPOT_RIL C1_RIL C2_RIL IV_EDA DEPOT_TOF C1_TOF C2_TOF DEPOT_PB C_PB
     MN_upper MN_lower SOD1_wt SOD1_mis TDP43_nuc TDP43_cyto
     Glut_syn Ca_i ROS GSH Mito Mic_act TNFa BDNF
     NfL_CSF ALSFRS FVC

$MAIN
double Cp_RIL   = C1_RIL  / V1_RIL;
double Cp_EDA   = IV_EDA  / V_EDA;
double Ccsf_TOF = C2_TOF  / V2_TOF;
double Cp_PB    = (C_PB   / V_PB) * 1000.0;
double E_RIL    = Emax_RIL * Cp_RIL / (IC50_RIL + Cp_RIL + 1e-9);
double E_EDA    = Emax_EDA * Cp_EDA / (IC50_EDA + Cp_EDA + 1e-9);
double E_TOF    = Emax_TOF * Ccsf_TOF / (EC50_TOF + Ccsf_TOF + 1e-9);
double E_PB     = Emax_PB  * Cp_PB  / (EC50_PB  + Cp_PB  + 1e-9);
double E_mPB    = Emax_mito_PB * Cp_PB / (EC50_PB + Cp_PB + 1e-9);
double total_SOD1 = SOD1_wt + SOD1_mis + 1e-9;
double frac_SOD1mis = SOD1_mis / total_SOD1;
double SOD1_burden  = f_SOD1_mut * frac_SOD1mis;
double frac_TDP_cyto = TDP43_cyto / (TDP43_nuc + TDP43_cyto + 1e-9);
double Glu_excess = fmax(0.0, Glut_syn - Glu_base);
double ROS_norm   = ROS / (ROS_0 + 1e-9);
double TNFa_norm  = TNFa / 1.0;
double MN_death_rate = k_MN_death * (1.0 + 3.0*SOD1_burden + 2.0*frac_TDP_cyto +
  1.5*ROS_norm + 1.2*Glu_excess + 0.8*TNFa_norm - BDNF_prot_wt*(BDNF/BDNF_0));
MN_death_rate = fmax(0.0001, MN_death_rate);
double EAAT2_eff  = k_EAAT2 * (k_EAAT2_ALS + (1.0 - k_EAAT2_ALS)/(1.0 + Mic_act));
double SOD1_syn_eff = k_SOD1_syn * (1.0 - E_TOF);
double NfL_release  = k_NfL_rel * MN_death_rate * (MN_upper + MN_lower);
double Mic_stim  = 0.5*SOD1_burden + 0.3*frac_TDP_cyto + 0.3*ROS_norm +
  0.5*(MN_upper_0 - MN_upper) + 0.5*(MN_lower_0 - MN_lower);

$ODE
dxdt_DEPOT_RIL = -ka_RIL * DEPOT_RIL;
dxdt_C1_RIL    = F_RIL * ka_RIL * DEPOT_RIL - (CL_RIL+Q_RIL)/V1_RIL*C1_RIL + Q_RIL/V2_RIL*C2_RIL;
dxdt_C2_RIL    = Q_RIL/V1_RIL*C1_RIL - Q_RIL/V2_RIL*C2_RIL;
dxdt_IV_EDA    = -(CL_EDA/V_EDA)*IV_EDA;
dxdt_DEPOT_TOF = -ka_TOF*DEPOT_TOF;
dxdt_C1_TOF    = ka_TOF*DEPOT_TOF - (CL_TOF+Q_TOF)/V1_TOF*C1_TOF + Q_TOF/V2_TOF*C2_TOF;
dxdt_C2_TOF    = Q_TOF/V1_TOF*C1_TOF - Q_TOF/V2_TOF*C2_TOF;
dxdt_DEPOT_PB  = -ka_PB*DEPOT_PB;
dxdt_C_PB      = F_PB*ka_PB*DEPOT_PB - (CL_PB/V_PB)*C_PB;
dxdt_MN_upper  = -k_MN_death_U * MN_death_rate * MN_upper;
dxdt_MN_lower  = -k_MN_death * MN_death_rate * MN_lower;
dxdt_SOD1_wt   = SOD1_syn_eff - k_SOD1_deg*SOD1_wt - k_SOD1_mis*SOD1_wt*f_SOD1_mut;
dxdt_SOD1_mis  = k_SOD1_mis*SOD1_wt*f_SOD1_mut - k_SOD1_clr*SOD1_mis;
dxdt_TDP43_nuc  = -k_TDP_export*TDP43_nuc + k_TDP_import*TDP43_cyto;
dxdt_TDP43_cyto =  k_TDP_export*TDP43_nuc - k_TDP_import*TDP43_cyto - k_TDP_agg*TDP43_cyto;
dxdt_Glut_syn  = k_Glu_rel*(1.0-E_RIL)*MN_lower - EAAT2_eff*Glut_syn;
dxdt_Ca_i      = k_Ca_entry*Glu_excess - k_Ca_efflux*(Ca_i - Ca_i_0);
double ROS_prod = k_ROS_prod*(1.0+SOD1_burden)*Ca_i/Mito;
double ROS_scav = k_ROS_scav*GSH*ROS + E_EDA*ROS;
dxdt_ROS       = ROS_prod - ROS_scav;
dxdt_GSH       = k_GSH_syn - k_GSH_cons*GSH - k_ROS_scav*GSH*ROS;
dxdt_Mito      = k_Mito_rep*(1.0-Mito)*(1.0+E_mPB) - k_Mito_dam*ROS*Mito;
dxdt_Mic_act   = k_Mic_act*Mic_stim*(1.0-Mic_act) - k_Mic_res*Mic_act;
dxdt_TNFa      = k_TNFa_prod*Mic_act - k_TNFa_deg*TNFa;
dxdt_BDNF      = k_BDNF_prod*(1.0-0.35*TNFa_norm) - k_BDNF_deg*BDNF;
dxdt_NfL_CSF   = NfL_release - k_NfL_clr*NfL_CSF;
double MN_frac  = 0.5*(MN_upper+MN_lower)/(0.5*(MN_upper_0+MN_lower_0));
dxdt_ALSFRS    = -k_MN_death*ALSFRS*(1.8-MN_frac)*(1.0+0.5*Mic_act);
dxdt_FVC       = -k_FVC_dec*FVC*(1.0+TNFa_norm+Mic_act);

$INIT
DEPOT_RIL=0 C1_RIL=0 C2_RIL=0 IV_EDA=0
DEPOT_TOF=0 C1_TOF=0 C2_TOF=0 DEPOT_PB=0 C_PB=0
MN_upper=1.0 MN_lower=1.0 SOD1_wt=5.0 SOD1_mis=0.0
TDP43_nuc=10.0 TDP43_cyto=0.5 Glut_syn=1.0 Ca_i=0.1
ROS=0.5 GSH=4.0 Mito=1.0 Mic_act=0.1 TNFa=0.2 BDNF=1.0
NfL_CSF=5.0 ALSFRS=48.0 FVC=100.0

$TABLE
capture Cp_RIL    = C1_RIL/V1_RIL;
capture Cp_EDA    = IV_EDA/V_EDA;
capture Ccsf_TOF  = C2_TOF/V2_TOF;
capture MN_pct    = 100.0*(MN_upper+MN_lower)/2.0;
capture ROS_norm_out = ROS/ROS_0;
'

als_mod <- mcode("als_shiny", als_code, quiet = TRUE)

## ── Scenario runner ──────────────────────────────────────────────────────────
run_als <- function(sc, days, f_SOD1 = 0.0, k_mult = 1.0,
                    use_RIL = FALSE, use_EDA = FALSE,
                    use_TOF = FALSE, use_PB = FALSE) {

  p_mod <- param(als_mod, list(f_SOD1_mut = f_SOD1,
                               k_MN_death = 0.00080 * k_mult))

  ev_list <- list()
  if (use_RIL) ev_list[["ril"]] <- ev(time = seq(0, (days-0.5)*24, 12),
                                       amt = 50, cmt = "DEPOT_RIL", evid = 1)
  if (use_EDA) {
    times_eda <- c(seq(0,13)*24,
                   unlist(lapply(1:floor(days/28), function(c) seq(28+(c-1)*28,
                                                                    min(37+(c-1)*28, days)*24, 24))))
    ev_list[["eda"]] <- ev(time = times_eda, amt = 60, cmt = "IV_EDA", evid = 1)
  }
  if (use_TOF) {
    tof_d <- unique(c(0, 14, 28, seq(56, days, 28)))
    ev_list[["tof"]] <- ev(time = tof_d * 24, amt = 100, cmt = "DEPOT_TOF", evid = 1)
  }
  if (use_PB) ev_list[["pb"]] <- ev(time = seq(0, (days-0.5)*24, 12),
                                     amt = 3000, cmt = "DEPOT_PB", evid = 1)

  all_ev <- if (length(ev_list) == 0) ev(time=0, amt=0, cmt=1, evid=2) else Reduce(c, ev_list)

  p_mod %>%
    ev(all_ev) %>%
    mrgsim(end = days * 24, delta = 6, obsonly = TRUE) %>%
    as_tibble() %>%
    mutate(time_days = time / 24)
}

## ── UI ───────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = span(icon("brain"), "ALS QSP Dashboard"),
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "profile",    icon = icon("user-md")),
      menuItem("Drug PK",             tabName = "pk",         icon = icon("pills")),
      menuItem("Biomarkers",          tabName = "biomarkers", icon = icon("vial")),
      menuItem("Clinical Endpoints",  tabName = "clinical",   icon = icon("chart-line")),
      menuItem("Scenario Comparison", tabName = "scenarios",  icon = icon("exchange-alt")),
      menuItem("Mechanistic Pathways",tabName = "pathways",   icon = icon("project-diagram"))
    ),

    hr(),
    h5(strong("  Patient Parameters"), style = "color:#B3E5FC; margin-left:8px"),
    sliderInput("days", "Duration (days):", min = 90, max = 730, value = 548, step = 30),
    selectInput("als_type", "ALS Subtype:",
                choices = c("Sporadic ALS (TDP-43)" = "0",
                            "SOD1-familial ALS"      = "1"),
                selected = "0"),
    sliderInput("prog_rate", "Progression rate multiplier:", min = 0.5, max = 2.0, value = 1.0, step = 0.1),

    hr(),
    h5(strong("  Drug Selection"), style = "color:#B3E5FC; margin-left:8px"),
    checkboxInput("use_ril",  "Riluzole 50 mg BID (PO)",        value = FALSE),
    checkboxInput("use_eda",  "Edaravone 60 mg/day (IV/PO)",    value = FALSE),
    checkboxInput("use_tof",  "Tofersen 100 mg SC q4w",         value = FALSE),
    checkboxInput("use_amx",  "AMX0035 (PB 3g + TUDCA 1g) BID", value = FALSE),

    hr(),
    actionButton("run_sim", "Run Simulation", icon = icon("play"),
                 class = "btn-primary btn-block")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .skin-blue .main-header .logo { background-color: #1A237E; font-weight: bold; }
      .skin-blue .main-header .navbar { background-color: #283593; }
      .content-wrapper { background-color: #f5f5f5; }
      .box-header { font-weight: bold; }
    "))),

    tabItems(

      ## ── TAB 1: Patient Profile ─────────────────────────────────────────────
      tabItem(tabName = "profile",
        fluidRow(
          box(width = 12, title = "ALS — Disease Overview", status = "primary",
              solidHeader = TRUE,
              HTML("<b>Amyotrophic Lateral Sclerosis (ALS)</b> is a fatal progressive neurodegenerative
              disease affecting upper and lower motor neurons. Median survival from symptom onset is
              2–5 years, with respiratory failure as the leading cause of death.
              ~90% of cases are sporadic; ~10% familial (SOD1, C9orf72, TDP-43/FUS mutations).
              Key mechanisms: <b>glutamate excitotoxicity, protein aggregation (TDP-43, SOD1),
              oxidative stress, neuroinflammation, and impaired axonal transport</b>."))
        ),
        fluidRow(
          valueBoxOutput("vbox_mn",    width = 3),
          valueBoxOutput("vbox_alsfrs",width = 3),
          valueBoxOutput("vbox_fvc",   width = 3),
          valueBoxOutput("vbox_nfl",   width = 3)
        ),
        fluidRow(
          box(width = 6, title = "Drug PK/PD Summary Table", status = "info", solidHeader = TRUE,
              DTOutput("pk_table")),
          box(width = 6, title = "Clinical Milestones", status = "warning", solidHeader = TRUE,
              DTOutput("milestone_table"))
        )
      ),

      ## ── TAB 2: Drug PK ────────────────────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(width = 6, title = "Riluzole Plasma Concentration (Cp)", status = "primary",
              solidHeader = TRUE, plotlyOutput("pk_ril_plot")),
          box(width = 6, title = "Edaravone Plasma Concentration", status = "success",
              solidHeader = TRUE, plotlyOutput("pk_eda_plot"))
        ),
        fluidRow(
          box(width = 6, title = "Tofersen CSF Concentration", status = "warning",
              solidHeader = TRUE, plotlyOutput("pk_tof_plot")),
          box(width = 6, title = "Drug Effect Summary", status = "info",
              solidHeader = TRUE, plotlyOutput("drug_effect_plot"))
        )
      ),

      ## ── TAB 3: Biomarkers ─────────────────────────────────────────────────
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(width = 6, title = "CSF Neurofilament Light Chain (NfL)", status = "primary",
              solidHeader = TRUE, plotlyOutput("nfl_plot")),
          box(width = 6, title = "Microglial Activation & TNF-α", status = "danger",
              solidHeader = TRUE, plotlyOutput("inflam_plot"))
        ),
        fluidRow(
          box(width = 6, title = "Reactive Oxygen Species & Glutathione", status = "warning",
              solidHeader = TRUE, plotlyOutput("redox_plot")),
          box(width = 6, title = "Mitochondrial Function & BDNF", status = "success",
              solidHeader = TRUE, plotlyOutput("mito_bdnf_plot"))
        )
      ),

      ## ── TAB 4: Clinical Endpoints ─────────────────────────────────────────
      tabItem(tabName = "clinical",
        fluidRow(
          box(width = 6, title = "ALSFRS-R Total Score", status = "primary",
              solidHeader = TRUE, plotlyOutput("alsfrs_plot")),
          box(width = 6, title = "FVC % Predicted", status = "danger",
              solidHeader = TRUE, plotlyOutput("fvc_plot"))
        ),
        fluidRow(
          box(width = 6, title = "Motor Neuron Survival (%)", status = "warning",
              solidHeader = TRUE, plotlyOutput("mn_plot")),
          box(width = 6, title = "Synaptic Glutamate & Intracellular Ca²⁺", status = "info",
              solidHeader = TRUE, plotlyOutput("excitotox_plot"))
        )
      ),

      ## ── TAB 5: Scenario Comparison ────────────────────────────────────────
      tabItem(tabName = "scenarios",
        fluidRow(
          box(width = 12, title = "Multi-Drug Scenario Comparison — ALSFRS-R",
              status = "primary", solidHeader = TRUE, plotlyOutput("scenario_alsfrs_plot"))
        ),
        fluidRow(
          box(width = 6, title = "Motor Neuron Survival — All Scenarios", status = "warning",
              solidHeader = TRUE, plotlyOutput("scenario_mn_plot")),
          box(width = 6, title = "CSF NfL — All Scenarios", status = "info",
              solidHeader = TRUE, plotlyOutput("scenario_nfl_plot"))
        ),
        fluidRow(
          box(width = 12, title = "Endpoint Summary at 18 months",
              status = "success", solidHeader = TRUE, DTOutput("scenario_table"))
        )
      ),

      ## ── TAB 6: Mechanistic Pathways ───────────────────────────────────────
      tabItem(tabName = "pathways",
        fluidRow(
          box(width = 6, title = "TDP-43 Nuclear/Cytoplasmic Dynamics", status = "danger",
              solidHeader = TRUE, plotlyOutput("tdp43_plot")),
          box(width = 6, title = "SOD1 Misfolding Dynamics (SOD1-ALS)", status = "warning",
              solidHeader = TRUE, plotlyOutput("sod1_plot"))
        ),
        fluidRow(
          box(width = 12, title = "Mechanistic Map (Graphviz)", status = "primary",
              solidHeader = TRUE, height = "600px",
              imageOutput("mech_map_img", height = "560px"))
        )
      )
    )
  )
)

## ── Server ───────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## ── Reactive simulation ─────────────────────────────────────────────────
  sim_data <- eventReactive(input$run_sim, {
    run_als(
      sc       = "user_defined",
      days     = input$days,
      f_SOD1   = as.numeric(input$als_type),
      k_mult   = input$prog_rate,
      use_RIL  = input$use_ril,
      use_EDA  = input$use_eda,
      use_TOF  = input$use_tof,
      use_PB   = input$use_amx
    )
  }, ignoreNULL = FALSE)

  ## Scenario comparison (all 7 arms at once)
  sc_data <- eventReactive(input$run_sim, {
    f_mut <- as.numeric(input$als_type)
    k_m   <- input$prog_rate
    days  <- input$days
    sc_list <- list(
      list(n = "Untreated",       r=F, e=F, t=F, p=F),
      list(n = "Riluzole",        r=T, e=F, t=F, p=F),
      list(n = "Edaravone",       r=F, e=T, t=F, p=F),
      list(n = "Ril+Eda",         r=T, e=T, t=F, p=F),
      list(n = "Tofersen",        r=F, e=F, t=T, p=F),
      list(n = "AMX0035",         r=F, e=F, t=F, p=T),
      list(n = "All Drugs",       r=T, e=T, t=T, p=T)
    )
    bind_rows(lapply(sc_list, function(x) {
      run_als(x$n, days, f_SOD1=f_mut, k_mult=k_m,
              use_RIL=x$r, use_EDA=x$e, use_TOF=x$t, use_PB=x$p) %>%
        mutate(scenario = x$n)
    }))
  }, ignoreNULL = FALSE)

  ## ── Value boxes ────────────────────────────────────────────────────────
  output$vbox_mn <- renderValueBox({
    d <- sim_data(); last <- tail(d, 1)
    valueBox(sprintf("%.0f%%", last$MN_pct), "Motor Neurons at End",
             icon = icon("brain"), color = if (last$MN_pct > 60) "blue" else "red")
  })
  output$vbox_alsfrs <- renderValueBox({
    d <- sim_data(); last <- tail(d, 1)
    valueBox(sprintf("%.1f", last$ALSFRS), "Final ALSFRS-R",
             icon = icon("clipboard"), color = "yellow")
  })
  output$vbox_fvc <- renderValueBox({
    d <- sim_data(); last <- tail(d, 1)
    valueBox(sprintf("%.0f%%", last$FVC), "Final FVC%",
             icon = icon("lungs"), color = if (last$FVC > 50) "green" else "red")
  })
  output$vbox_nfl <- renderValueBox({
    d <- sim_data(); last <- tail(d, 1)
    valueBox(sprintf("%.1f", last$NfL_CSF), "CSF NfL (AU)",
             icon = icon("vial"), color = "purple")
  })

  ## ── PK tables ──────────────────────────────────────────────────────────
  output$pk_table <- renderDT({
    data.frame(
      Drug      = c("Riluzole", "Edaravone", "Tofersen", "PB (AMX0035)"),
      Dose      = c("50 mg BID PO","60 mg/day IV","100 mg SC q4w","3 g BID PO"),
      F_pct     = c("60%","~100% IV","NA","85%"),
      t_half    = c("12 h","4.5 h","7 days","3 h"),
      CL        = c("28 L/h","18 L/h","0.5 L/h","12 L/h"),
      Mechanism = c("↓ Glu release","ROS scavenger","SOD1 ASO","ER stress / mito protect")
    )
  }, options = list(pageLength = 5, dom = "t"), rownames = FALSE)

  output$milestone_table <- renderDT({
    d <- sim_data()
    alsfrs_6pt  <- min(d$time_days[d$ALSFRS <= (d$ALSFRS[1] - 6)], na.rm = TRUE)
    fvc_50      <- min(d$time_days[d$FVC    <= 50], na.rm = TRUE)
    data.frame(
      Milestone = c("Time to 6-pt ALSFRS-R decline",
                    "Time to FVC ≤ 50%",
                    "ALSFRS-R at 12 months",
                    "FVC at 12 months",
                    "MN survival at 18 months"),
      Value     = c(
        if (is.finite(alsfrs_6pt)) paste0(round(alsfrs_6pt), " days") else ">18 mo",
        if (is.finite(fvc_50))     paste0(round(fvc_50),     " days") else ">18 mo",
        round(d$ALSFRS[which.min(abs(d$time_days - 365))], 1),
        paste0(round(d$FVC[which.min(abs(d$time_days - 365))], 1), "%"),
        paste0(round(tail(d$MN_pct, 1), 1), "%")
      )
    )
  }, options = list(dom = "t"), rownames = FALSE)

  ## ── PK Plots ──────────────────────────────────────────────────────────
  output$pk_ril_plot <- renderPlotly({
    d <- sim_data() %>% filter(time_days <= 7)
    plot_ly(d, x = ~time_days*24, y = ~Cp_RIL, type = "scatter", mode = "lines",
            line = list(color="#1976D2", width=2)) %>%
      add_lines(y = rep(0.5, nrow(d)), line = list(color="red", dash="dash"),
                name = "IC50 = 0.5 μg/mL") %>%
      layout(xaxis = list(title="Time (h)"), yaxis = list(title="Riluzole Cp (μg/mL)"),
             title = "Riluzole PK (first 7 days)")
  })
  output$pk_eda_plot <- renderPlotly({
    d <- sim_data() %>% filter(time_days <= 30)
    plot_ly(d, x = ~time_days, y = ~Cp_EDA, type = "scatter", mode = "lines",
            line = list(color="#388E3C", width=2)) %>%
      add_lines(y = rep(1.2, nrow(d)), line = list(color="red", dash="dash"),
                name = "IC50 = 1.2 μg/mL") %>%
      layout(xaxis = list(title="Time (days)"), yaxis = list(title="Edaravone Cp (μg/mL)"),
             title = "Edaravone PK (first 30 days)")
  })
  output$pk_tof_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~Ccsf_TOF, type = "scatter", mode = "lines",
            line = list(color="#F57C00", width=2)) %>%
      add_lines(y = rep(0.1, nrow(d)), line = list(color="red", dash="dash"),
                name = "EC50 = 0.1 μg/mL") %>%
      layout(xaxis = list(title="Time (days)"), yaxis = list(title="Tofersen CSF (μg/mL)"),
             title = "Tofersen CSF Concentration")
  })
  output$drug_effect_plot <- renderPlotly({
    d <- sim_data()
    E_ril_col <- Inf; E_eda_col <- Inf; E_tof_col <- Inf
    if ("Cp_RIL" %in% names(d)) {
      d <- d %>% mutate(
        E_ril = 0.60 * Cp_RIL / (0.50 + Cp_RIL),
        E_eda = 0.70 * Cp_EDA / (1.20 + Cp_EDA),
        E_tof = 0.80 * Ccsf_TOF / (0.10 + Ccsf_TOF)
      )
    }
    plot_ly(d, x = ~time_days) %>%
      add_lines(y = ~E_ril*100, name = "Riluzole % Glu inh.", line = list(color="#1976D2")) %>%
      add_lines(y = ~E_eda*100, name = "Edaravone % ROS scav.", line = list(color="#388E3C")) %>%
      add_lines(y = ~E_tof*100, name = "Tofersen % SOD1 KD", line = list(color="#F57C00")) %>%
      layout(xaxis = list(title="Time (days)"), yaxis = list(title="Drug Effect (%)"))
  })

  ## ── Biomarker plots ────────────────────────────────────────────────────
  output$nfl_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~NfL_CSF, type="scatter", mode="lines",
            line=list(color="#6A1B9A", width=2)) %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="CSF NfL (AU)"),
             title="CSF Neurofilament Light Chain")
  })
  output$inflam_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_days) %>%
      add_lines(y=~Mic_act*100, name="Microglia Activation %", line=list(color="#F57F17")) %>%
      add_lines(y=~TNFa*10, name="TNF-α ×10 (AU)", line=list(color="#D32F2F")) %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="Value"),
             title="Neuroinflammation Dynamics")
  })
  output$redox_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_days) %>%
      add_lines(y=~ROS_norm_out, name="ROS (normalized)", line=list(color="#E65100")) %>%
      add_lines(y=~GSH/4.0,     name="GSH (normalized)", line=list(color="#2E7D32")) %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="Normalized Level"),
             title="Oxidative Stress — ROS vs GSH")
  })
  output$mito_bdnf_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_days) %>%
      add_lines(y=~Mito*100, name="Mitochondrial function %", line=list(color="#00838F")) %>%
      add_lines(y=~BDNF*100, name="BDNF (normalized ×100)", line=list(color="#1976D2")) %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="% Baseline"),
             title="Mitochondrial Function & BDNF Trophic Support")
  })

  ## ── Clinical plots ──────────────────────────────────────────────────────
  output$alsfrs_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_days, y=~ALSFRS, type="scatter", mode="lines",
            line=list(color="#1A237E", width=2.5)) %>%
      add_lines(y=rep(36, nrow(d)), line=list(color="orange", dash="dash"),
                name="Moderate disability (36)") %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="ALSFRS-R Score", range=c(0,48)),
             title="ALSFRS-R Total Score Progression")
  })
  output$fvc_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_days, y=~FVC, type="scatter", mode="lines",
            line=list(color="#D32F2F", width=2.5)) %>%
      add_lines(y=rep(50, nrow(d)), line=list(color="red", dash="dash"), name="NIV threshold") %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="FVC %", range=c(0,105)),
             title="Forced Vital Capacity (FVC) Decline")
  })
  output$mn_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_days) %>%
      add_lines(y=~MN_upper*100, name="Upper MN", line=list(color="#1976D2")) %>%
      add_lines(y=~MN_lower*100, name="Lower MN", line=list(color="#D32F2F")) %>%
      add_lines(y=~MN_pct,       name="Average",  line=list(color="#37474F", dash="dash")) %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="Survival %"),
             title="Upper & Lower Motor Neuron Survival")
  })
  output$excitotox_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_days) %>%
      add_lines(y=~Glut_syn, name="Synaptic Glu (AU)", line=list(color="#1565C0")) %>%
      add_lines(y=~Ca_i*10,  name="Ca²⁺ ×10 (μM)",   line=list(color="#880E4F")) %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="Concentration"),
             title="Excitotoxicity — Glutamate & Calcium")
  })

  ## ── Scenario comparison ─────────────────────────────────────────────────
  sc_colors <- c(
    "Untreated"="red", "Riluzole"="blue", "Edaravone"="green",
    "Ril+Eda"="purple", "Tofersen"="orange", "AMX0035"="teal", "All Drugs"="brown"
  )
  output$scenario_alsfrs_plot <- renderPlotly({
    d <- sc_data()
    p <- plot_ly(d, x=~time_days, y=~ALSFRS, color=~scenario, type="scatter", mode="lines",
                 colors=sc_colors) %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="ALSFRS-R", range=c(0,50)),
             title="ALSFRS-R — All 7 Scenarios", legend=list(orientation="h"))
    p
  })
  output$scenario_mn_plot <- renderPlotly({
    d <- sc_data()
    plot_ly(d, x=~time_days, y=~MN_pct, color=~scenario, type="scatter", mode="lines",
            colors=sc_colors) %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="MN Survival %"))
  })
  output$scenario_nfl_plot <- renderPlotly({
    d <- sc_data()
    plot_ly(d, x=~time_days, y=~NfL_CSF, color=~scenario, type="scatter", mode="lines",
            colors=sc_colors) %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="CSF NfL (AU)"))
  })
  output$scenario_table <- renderDT({
    d <- sc_data()
    target_day <- 548
    d %>%
      group_by(scenario) %>%
      summarise(
        ALSFRS_0mo  = round(ALSFRS[which.min(abs(time_days - 0))],   1),
        ALSFRS_18mo = round(ALSFRS[which.min(abs(time_days - target_day))], 1),
        FVC_18mo    = paste0(round(FVC[which.min(abs(time_days - target_day))], 1), "%"),
        MN_surv_pct = paste0(round(MN_pct[which.min(abs(time_days - target_day))], 1), "%"),
        NfL_18mo    = round(NfL_CSF[which.min(abs(time_days - target_day))], 2),
        .groups = "drop"
      ) %>%
      rename(Scenario=scenario, `ALSFRS(0mo)`=ALSFRS_0mo,
             `ALSFRS(18mo)`=ALSFRS_18mo, `FVC(18mo)`=FVC_18mo,
             `MN Survival`=MN_surv_pct, `NfL(18mo)`=NfL_18mo)
  }, options = list(dom = "t"), rownames = FALSE)

  ## ── Mechanistic pathways ──────────────────────────────────────────────
  output$tdp43_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_days) %>%
      add_lines(y=~TDP43_nuc,  name="Nuclear TDP-43", line=list(color="#0D47A1")) %>%
      add_lines(y=~TDP43_cyto, name="Cytoplasmic TDP-43", line=list(color="#D32F2F")) %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="TDP-43 (AU)"),
             title="TDP-43 Nuclear↔Cytoplasmic Shuttle")
  })
  output$sod1_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_days) %>%
      add_lines(y=~SOD1_wt,  name="WT SOD1",       line=list(color="#2E7D32")) %>%
      add_lines(y=~SOD1_mis, name="Misfolded SOD1", line=list(color="#C62828")) %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="SOD1 (AU)"),
             title="SOD1 Protein Dynamics (relevant in SOD1-ALS)")
  })
  output$mech_map_img <- renderImage({
    img_path <- file.path(getwd(), "als_qsp_model.png")
    if (!file.exists(img_path)) {
      img_path <- file.path(dirname(getwd()), "als_qsp_model.png")
    }
    list(src = img_path, contentType = "image/png",
         width = "100%", alt = "ALS Mechanistic Map")
  }, deleteFile = FALSE)
}

## ── Launch ───────────────────────────────────────────────────────────────────
shinyApp(ui, server)
