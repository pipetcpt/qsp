## ============================================================
## Renal Cell Carcinoma (ccRCC) — Interactive Shiny Dashboard
## 7 tabs: Patient Profile · PK · VHL/HIF/VEGF Pathway ·
##         Tumor Dynamics · Immune TME · Scenario Comparison ·
##         Biomarker Dashboard
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)

## ------------------------------------------------------------------
## Inline mrgsolve model (stripped-down for Shiny speed)
## ------------------------------------------------------------------
rcc_code <- '
$PARAM
CL_sun=51.8 V1_sun=2030 Q_sun=7.5 V2_sun=560 ka_sun=0.09 F_sun=0.50
fm_sun=0.23 CL_met=14.2 V_met=1100
CL_niv=0.019 V1_niv=3.5 Q_niv=0.003 V2_niv=3.8
kon_niv=0.32 koff_niv=0.0018 kdeg_niv=0.034 Rmax_PD1=6.5 ksyn_PD1=0.22
CL_bez=2.5 V_bez=85 ka_bez=0.6 F_bez=0.70 IC50_bez=0.018
kprod_VHL=0.12 kdeg_VHL=0.08 VHL_frac=0.25
kprod_HIF2=0.18 kdeg_HIF2=0.35 kprod_VEGF=0.55 kdeg_VEGF=0.14 EC50_HIF2=1.8
IC50_sun_VEGFR=0.008 VEGFR2_base=2.5 kprod_VEGFR2=0.12 kdeg_VEGFR2=0.05
kprod_mTOR=0.22 kdeg_mTOR=0.18 IC50_ever=0.15 kp_mTOR=1.2
lambda1=5.0e-4 lambda2=1.8e-2 psi=20 w0=0.8
k1_sun=0.045 k2_dam=0.032
kprol_CD8=0.09 kdeath_CD8=0.04 alpha_CD8=2.5e-4
ki_treg=0.35 k_MDSC=0.28 kprol_Treg=0.06 kdeath_Treg=0.04
kprol_MDSC=0.07 kdeath_MDSC=0.05
CD8_base=100 Treg_base=15 MDSC_base=25 EC50_PD1=0.5
CL_ever=8.4 V_ever=191 ka_ever=1.2 F_ever=0.16
CL_axi=12.5 V_axi=160 ka_axi=1.0 F_axi=0.58 IC50_axi=0.003
CL_cabo=2.5 V_cabo=350 ka_cabo=0.8 F_cabo=0.57 IC50_cabo=0.006
use_sunitinib=0 use_nivolumab=0 use_belzutifan=0
use_everolimus=0 use_axitinib=0 use_cabozantinib=0

$CMT
DEPOT_SUN CENT_SUN PERI_SUN MET_SUN
CENT_NIV PERI_NIV PD1_FREE PD1_BOUND
DEPOT_BEZ CENT_BEZ
pVHL HIF2A VEGF VEGFR2_ACT mTOR_ACT
DEPOT_EVER CENT_EVER DEPOT_AXI CENT_AXI DEPOT_CABO CENT_CABO
TUM_W1 TUM_W2 TUM_W3 TUM_VOL
CD8_T TREG MDSC

$MAIN
DEPOT_SUN_0=0; CENT_SUN_0=0; PERI_SUN_0=0; MET_SUN_0=0;
CENT_NIV_0=0; PERI_NIV_0=0;
PD1_FREE_0=Rmax_PD1; PD1_BOUND_0=0;
DEPOT_BEZ_0=0; CENT_BEZ_0=0;
pVHL_0=VHL_frac*kprod_VHL/kdeg_VHL;
HIF2A_0=kprod_HIF2/kdeg_HIF2*(1.0/(1.0+pVHL_0));
VEGF_0=kprod_VEGF*pow(HIF2A_0,2.0)/(pow(EC50_HIF2,2.0)+pow(HIF2A_0,2.0))/kdeg_VEGF;
VEGFR2_ACT_0=VEGFR2_base; mTOR_ACT_0=kprod_mTOR/kdeg_mTOR;
DEPOT_EVER_0=0; CENT_EVER_0=0;
DEPOT_AXI_0=0; CENT_AXI_0=0;
DEPOT_CABO_0=0; CENT_CABO_0=0;
TUM_W1_0=0; TUM_W2_0=0; TUM_W3_0=0; TUM_VOL_0=w0;
CD8_T_0=CD8_base; TREG_0=Treg_base; MDSC_0=MDSC_base;

$ODE
double C_sun=CENT_SUN/V1_sun;
double C_met=MET_SUN/V_met;
double C_sun_eq=(C_sun+0.5*C_met)/1000.0;
double C_niv=CENT_NIV/V1_niv;
double C_bez=CENT_BEZ/V_bez;
double C_ever=CENT_EVER/V_ever;
double C_axi=CENT_AXI/V_axi;
double C_cabo=CENT_CABO/V_cabo;

dxdt_DEPOT_SUN=-ka_sun*F_sun*DEPOT_SUN;
dxdt_CENT_SUN=ka_sun*F_sun*DEPOT_SUN-CL_sun*C_sun-Q_sun*(C_sun-PERI_SUN/V2_sun);
dxdt_PERI_SUN=Q_sun*(C_sun-PERI_SUN/V2_sun);
dxdt_MET_SUN=fm_sun*CL_sun*C_sun-CL_met*C_met;

dxdt_CENT_NIV=-CL_niv*C_niv-Q_niv*(C_niv-PERI_NIV/V2_niv)-kon_niv*C_niv*PD1_FREE+koff_niv*PD1_BOUND;
dxdt_PERI_NIV=Q_niv*(C_niv-PERI_NIV/V2_niv);
dxdt_PD1_FREE=ksyn_PD1-kdeg_niv*PD1_FREE-kon_niv*C_niv*PD1_FREE+koff_niv*PD1_BOUND;
dxdt_PD1_BOUND=kon_niv*C_niv*PD1_FREE-(koff_niv+kdeg_niv)*PD1_BOUND;

double PD1_occ=PD1_BOUND/(PD1_FREE+PD1_BOUND+1e-9);
double inh_bez_HIF2=C_bez/(IC50_bez+C_bez);
double stim_mTOR=kp_mTOR*mTOR_ACT/(kprod_mTOR/kdeg_mTOR);

dxdt_DEPOT_BEZ=-ka_bez*F_bez*DEPOT_BEZ;
dxdt_CENT_BEZ=ka_bez*F_bez*DEPOT_BEZ-CL_bez*C_bez;

dxdt_pVHL=kprod_VHL*VHL_frac-kdeg_VHL*pVHL;
dxdt_HIF2A=kprod_HIF2*stim_mTOR-kdeg_HIF2*(pVHL+0.01)*HIF2A-inh_bez_HIF2*kdeg_HIF2*HIF2A;
dxdt_VEGF=kprod_VEGF*pow(HIF2A,2.0)/(pow(EC50_HIF2,2.0)+pow(HIF2A,2.0))-kdeg_VEGF*VEGF;

double inh_sun=use_sunitinib*C_sun_eq/(IC50_sun_VEGFR+C_sun_eq);
double inh_axi=use_axitinib*C_axi/(IC50_axi+C_axi);
double inh_cabo=use_cabozantinib*C_cabo/(IC50_cabo+C_cabo);
double tot_inh=1.0-(1.0-inh_sun)*(1.0-inh_axi)*(1.0-inh_cabo);

dxdt_VEGFR2_ACT=kprod_VEGFR2*(1.0-tot_inh)-kdeg_VEGFR2*VEGFR2_ACT;

double inh_ever=use_everolimus*C_ever/(IC50_ever+C_ever);
dxdt_mTOR_ACT=kprod_mTOR*(1.0-inh_ever)-kdeg_mTOR*mTOR_ACT;

dxdt_DEPOT_EVER=-ka_ever*F_ever*DEPOT_EVER;
dxdt_CENT_EVER=ka_ever*F_ever*DEPOT_EVER-CL_ever*C_ever;
dxdt_DEPOT_AXI=-ka_axi*F_axi*DEPOT_AXI;
dxdt_CENT_AXI=ka_axi*F_axi*DEPOT_AXI-CL_axi*C_axi;
dxdt_DEPOT_CABO=-ka_cabo*F_cabo*DEPOT_CABO;
dxdt_CENT_CABO=ka_cabo*F_cabo*DEPOT_CABO-CL_cabo*C_cabo;

double TW=TUM_W1+TUM_W2+TUM_W3;
double growth=lambda1*TUM_VOL/pow(1.0+pow(lambda1/lambda2*TUM_VOL,psi),1.0/psi);
double dmg_sun=use_sunitinib*k1_sun*C_sun_eq;
double dmg_cabo=use_cabozantinib*k1_sun*C_cabo;
double CD8_eff=CD8_T*(1.0-ki_treg*TREG/(Treg_base+TREG))*(1.0-k_MDSC*MDSC/(MDSC_base+MDSC))*(1.0+3.0*PD1_occ/(EC50_PD1+PD1_occ));
double immune_kill=alpha_CD8*CD8_eff*TUM_VOL;
double total_k=dmg_sun+dmg_cabo+immune_kill;

dxdt_TUM_W1=total_k*TUM_VOL-k2_dam*TUM_W1;
dxdt_TUM_W2=k2_dam*(TUM_W1-TUM_W2);
dxdt_TUM_W3=k2_dam*(TUM_W2-TUM_W3);
dxdt_TUM_VOL=growth-k2_dam*TW;

double stim_veg=1.0+0.8*tot_inh;
dxdt_CD8_T=kprol_CD8*CD8_base*stim_veg*(1.0+2.0*PD1_occ)-kdeath_CD8*CD8_T-ki_treg*TREG*CD8_T/100.0-k_MDSC*MDSC*CD8_T/100.0;
dxdt_TREG=kprol_Treg*Treg_base-kdeath_Treg*TREG;
dxdt_MDSC=kprol_MDSC*MDSC_base*(1.0+0.5*VEGF/(VEGF+1.0))-kdeath_MDSC*MDSC;

$TABLE
double Csun_ngmL=CENT_SUN/V1_sun;
double Cmet_ngmL=MET_SUN/V_met;
double Cniv_nM=CENT_NIV/V1_niv;
double Cbez_uM=CENT_BEZ/V_bez;
double HIF2_nM=HIF2A;
double VEGF_nM=VEGF;
double VEGFR2_nM=VEGFR2_ACT;
double mTOR_AU=mTOR_ACT;
double TumorVol=TUM_VOL;
double CD8_count=CD8_T;
double Treg_count=TREG;
double MDSC_count=MDSC;
double PD1_occ_pct=PD1_BOUND/(PD1_FREE+PD1_BOUND+1e-9)*100.0;

$CAPTURE Csun_ngmL Cmet_ngmL Cniv_nM Cbez_uM
         HIF2_nM VEGF_nM VEGFR2_nM mTOR_AU TumorVol
         CD8_count Treg_count MDSC_count PD1_occ_pct
'

mod <- mcode("rcc_shiny", rcc_code)

## ------------------------------------------------------------------
## Helper: build dosing and set params
## ------------------------------------------------------------------
run_scenario <- function(regimen, dose_sun = 50, dose_niv = 240,
                         dose_bez = 120, dose_ever = 10,
                         dose_axi = 5, dose_cabo = 40,
                         n_weeks = 52, vhl_frac = 0.25) {

  evs <- list()
  p   <- list(VHL_frac = vhl_frac,
              use_sunitinib = 0, use_nivolumab = 0,
              use_belzutifan = 0, use_everolimus = 0,
              use_axitinib = 0, use_cabozantinib = 0)

  if (regimen == "Sunitinib") {
    evs[[1]] <- ev(amt = dose_sun * 1000, ii = 24, addl = n_weeks * 7 - 1,
                   cmt = "DEPOT_SUN")
    p$use_sunitinib <- 1

  } else if (regimen == "Nivo+Ipi") {
    evs[[1]] <- ev(amt = dose_niv * 1e6, ii = 14 * 24,
                   addl = floor(n_weeks / 2) - 1, cmt = "CENT_NIV")
    p$use_nivolumab <- 1

  } else if (regimen == "Pembro+Axitinib") {
    evs[[1]] <- ev(amt = 200 * 1e6, ii = 21 * 24,
                   addl = floor(n_weeks / 3) - 1, cmt = "CENT_NIV")
    evs[[2]] <- ev(amt = dose_axi * 1000, ii = 12, addl = n_weeks * 14 - 1,
                   cmt = "DEPOT_AXI")
    p$use_nivolumab <- 1; p$use_axitinib <- 1

  } else if (regimen == "Cabo+Nivo") {
    evs[[1]] <- ev(amt = dose_cabo * 1000, ii = 24,
                   addl = n_weeks * 7 - 1, cmt = "DEPOT_CABO")
    evs[[2]] <- ev(amt = dose_niv * 1e6, ii = 14 * 24,
                   addl = floor(n_weeks / 2) - 1, cmt = "CENT_NIV")
    p$use_cabozantinib <- 1; p$use_nivolumab <- 1

  } else if (regimen == "Cabozantinib") {
    evs[[1]] <- ev(amt = 60 * 1000, ii = 24, addl = n_weeks * 7 - 1,
                   cmt = "DEPOT_CABO")
    p$use_cabozantinib <- 1

  } else if (regimen == "Everolimus") {
    evs[[1]] <- ev(amt = dose_ever * 1000, ii = 24, addl = n_weeks * 7 - 1,
                   cmt = "DEPOT_EVER")
    p$use_everolimus <- 1

  } else if (regimen == "Belzutifan") {
    evs[[1]] <- ev(amt = dose_bez, ii = 24, addl = n_weeks * 7 - 1,
                   cmt = "DEPOT_BEZ")
    p$use_belzutifan <- 1

  } else {
    evs[[1]] <- ev(amt = 0, time = 0, cmt = "DEPOT_SUN")
  }

  ev_combined <- if (length(evs) > 1) do.call(rbind, evs) else evs[[1]]

  out <- mod %>%
    param(p) %>%
    mrgsim(events = ev_combined,
           end = n_weeks * 7 * 24,
           delta = 12,
           obsonly = TRUE) %>%
    as.data.frame() %>%
    mutate(time_wk = time / (7 * 24))

  out
}

## ------------------------------------------------------------------
## UI
## ------------------------------------------------------------------
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "ccRCC QSP Model"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",      tabName = "tab_patient",  icon = icon("user")),
      menuItem("PK",                   tabName = "tab_pk",       icon = icon("flask")),
      menuItem("VHL/HIF/VEGF Pathway", tabName = "tab_pathway",  icon = icon("dna")),
      menuItem("Tumor Dynamics",       tabName = "tab_tumor",    icon = icon("chart-line")),
      menuItem("Immune TME",           tabName = "tab_immune",   icon = icon("shield-alt")),
      menuItem("Scenario Comparison",  tabName = "tab_compare",  icon = icon("sliders-h")),
      menuItem("Biomarker Dashboard",  tabName = "tab_biomarker",icon = icon("microscope"))
    )
  ),
  dashboardBody(
    tabItems(

      ## ── Tab 1: Patient Profile ──────────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient Parameters", width = 4, status = "primary",
            selectInput("regimen", "Regimen",
              choices = c("Untreated", "Sunitinib", "Nivo+Ipi",
                          "Pembro+Axitinib", "Cabo+Nivo",
                          "Cabozantinib", "Everolimus", "Belzutifan"),
              selected = "Sunitinib"),
            sliderInput("vhl_frac", "Residual pVHL (VHL mutation severity)",
                        min = 0.05, max = 1.0, value = 0.25, step = 0.05),
            sliderInput("imdc_score", "IMDC Risk Score (0–6)",
                        min = 0, max = 6, value = 2, step = 1),
            sliderInput("baseline_tumor", "Baseline Tumor Volume (cm³)",
                        min = 0.2, max = 10, value = 0.8, step = 0.2),
            sliderInput("n_weeks", "Simulation Duration (weeks)",
                        min = 12, max = 104, value = 52, step = 4),
            actionButton("run_sim", "Run Simulation", class = "btn-primary btn-block")
          ),
          box(title = "IMDC Risk Stratification", width = 8, status = "info",
            DTOutput("imdc_table"),
            br(),
            valueBoxOutput("vbox_risk"),
            valueBoxOutput("vbox_mOS"),
            valueBoxOutput("vbox_mPFS")
          )
        )
      ),

      ## ── Tab 2: PK ───────────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Drug Concentrations",
              width = 12, status = "primary",
              plotlyOutput("pk_plot", height = "450px"))
        ),
        fluidRow(
          box(title = "PK Parameters", width = 6,
            sliderInput("dose_sun",  "Sunitinib Dose (mg)", 25, 75, 50, 12.5),
            sliderInput("dose_niv",  "Nivolumab Dose (mg)", 120, 480, 240, 60)
          ),
          box(title = "PK Summary", width = 6, DTOutput("pk_summary_table"))
        )
      ),

      ## ── Tab 3: VHL/HIF/VEGF Pathway ────────────────────────────
      tabItem(tabName = "tab_pathway",
        fluidRow(
          box(title = "Pathway Biomarkers Over Time",
              width = 12, status = "primary",
              plotlyOutput("pathway_plot", height = "450px"))
        ),
        fluidRow(
          box(title = "HIF-2α Steady State vs VHL Fraction",
              width = 6, plotlyOutput("hif_curve", height = "300px")),
          box(title = "VEGF vs HIF-2α",
              width = 6, plotlyOutput("vegf_hif_scatter", height = "300px"))
        )
      ),

      ## ── Tab 4: Tumor Dynamics ───────────────────────────────────
      tabItem(tabName = "tab_tumor",
        fluidRow(
          box(title = "Tumor Volume Dynamics (Simeoni TGI)",
              width = 12, status = "primary",
              plotlyOutput("tumor_plot", height = "450px"))
        ),
        fluidRow(
          box(title = "Waterfall: Best Overall Response",
              width = 6, plotlyOutput("waterfall_plot", height = "300px")),
          box(title = "Tumor Growth Kinetics", width = 6,
            verbatimTextOutput("tgi_stats"))
        )
      ),

      ## ── Tab 5: Immune TME ───────────────────────────────────────
      tabItem(tabName = "tab_immune",
        fluidRow(
          box(title = "Immune Cell Dynamics",
              width = 12, status = "primary",
              plotlyOutput("immune_plot", height = "400px"))
        ),
        fluidRow(
          box(title = "PD-1 Occupancy",
              width = 6, plotlyOutput("pd1_plot", height = "300px")),
          box(title = "MDSC & VEGF Crosstalk",
              width = 6, plotlyOutput("mdsc_vegf_plot", height = "300px"))
        )
      ),

      ## ── Tab 6: Scenario Comparison ──────────────────────────────
      tabItem(tabName = "tab_compare",
        fluidRow(
          box(title = "Regimen Selection", width = 3, status = "warning",
            checkboxGroupInput("sc_select", "Select Regimens",
              choices = c("Untreated", "Sunitinib", "Nivo+Ipi",
                          "Pembro+Axitinib", "Cabo+Nivo",
                          "Cabozantinib", "Everolimus", "Belzutifan"),
              selected = c("Untreated", "Sunitinib", "Nivo+Ipi",
                           "Pembro+Axitinib", "Cabo+Nivo")),
            actionButton("run_compare", "Compare", class = "btn-warning btn-block")
          ),
          box(title = "Tumor Volume Comparison",
              width = 9, plotlyOutput("compare_tumor", height = "400px"))
        ),
        fluidRow(
          box(title = "Endpoint Summary Table", width = 12,
              DTOutput("compare_table"))
        )
      ),

      ## ── Tab 7: Biomarker Dashboard ──────────────────────────────
      tabItem(tabName = "tab_biomarker",
        fluidRow(
          valueBoxOutput("bm_HIF2_box",    width = 3),
          valueBoxOutput("bm_VEGF_box",    width = 3),
          valueBoxOutput("bm_mTOR_box",    width = 3),
          valueBoxOutput("bm_CD8_box",     width = 3)
        ),
        fluidRow(
          box(title = "Multi-Biomarker Heatmap (Week 12 snapshot)",
              width = 12, status = "success",
              plotlyOutput("biomarker_heat", height = "400px"))
        )
      )

    )
  )
)

## ------------------------------------------------------------------
## Server
## ------------------------------------------------------------------
server <- function(input, output, session) {

  ## IMDC reference table
  output$imdc_table <- renderDT({
    tibble(
      `Risk Factor` = c("KPS < 80", "Hgb < LLN", "Ca > ULN",
                        "Neutrophils > ULN", "Platelets > ULN",
                        "< 1 yr diagnosis to treatment"),
      `Score` = rep(1, 6),
      `Impact` = c("Performance status", "Anaemia", "Hypercalcaemia",
                   "Neutrophilia", "Thrombocytosis", "Rapid progression")
    )
  }, options = list(dom = 't', pageLength = 6), rownames = FALSE)

  imdc_risk <- reactive({
    s <- input$imdc_score
    if (s == 0) "Favourable (0 factors)"
    else if (s <= 2) "Intermediate (1–2 factors)"
    else "Poor (3–6 factors)"
  })

  imdc_mOS <- reactive({
    s <- input$imdc_score
    if (s == 0) "Not reached" else if (s <= 2) "~27 months" else "~8.8 months"
  })

  imdc_mPFS <- reactive({
    s <- input$imdc_score
    if (s == 0) "~26 months" else if (s <= 2) "~14 months" else "~6 months"
  })

  output$vbox_risk <- renderValueBox(
    valueBox(imdc_risk(), "IMDC Risk Group", icon = icon("flag"), color = "blue"))
  output$vbox_mOS  <- renderValueBox(
    valueBox(imdc_mOS(),  "Median OS (ICI era)", icon = icon("heartbeat"), color = "green"))
  output$vbox_mPFS <- renderValueBox(
    valueBox(imdc_mPFS(), "Median PFS (ICI era)", icon = icon("chart-bar"), color = "yellow"))

  ## Core simulation (reactive to run_sim button)
  sim_data <- eventReactive(input$run_sim, {
    run_scenario(input$regimen,
                 n_weeks  = input$n_weeks,
                 vhl_frac = input$vhl_frac)
  }, ignoreNULL = FALSE)

  ## ── PK plots ──────────────────────────────────────────────────
  output$pk_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_wk) %>%
      add_lines(y = ~Csun_ngmL, name = "Sunitinib (ng/mL)", line = list(color = "steelblue")) %>%
      add_lines(y = ~Cmet_ngmL, name = "SU12662 metabolite (ng/mL)", line = list(color = "lightblue", dash = "dash")) %>%
      add_lines(y = ~Cniv_nM,   name = "Nivolumab (nM)", line = list(color = "orange")) %>%
      add_lines(y = ~Cbez_uM,   name = "Belzutifan (µM)", line = list(color = "purple")) %>%
      layout(title = "Drug Concentration–Time Profiles",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Concentration"),
             legend = list(orientation = "h"))
  })

  output$pk_summary_table <- renderDT({
    df <- sim_data()
    last_row <- tail(df, 1)
    tibble(
      Parameter    = c("Sunitinib Css (ng/mL)", "SU12662 Css (ng/mL)",
                       "Nivolumab Css (nM)", "Belzutifan Css (µM)"),
      Value        = round(c(last_row$Csun_ngmL, last_row$Cmet_ngmL,
                              last_row$Cniv_nM,   last_row$Cbez_uM), 3)
    )
  }, options = list(dom = 't'), rownames = FALSE)

  ## ── Pathway plots ──────────────────────────────────────────────
  output$pathway_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_wk) %>%
      add_lines(y = ~HIF2_nM,   name = "HIF-2α (nM)", line = list(color = "red")) %>%
      add_lines(y = ~VEGF_nM,   name = "VEGF (nM)",   line = list(color = "darkorange")) %>%
      add_lines(y = ~VEGFR2_nM, name = "VEGFR2 (nM)", line = list(color = "gold")) %>%
      add_lines(y = ~mTOR_AU,   name = "mTOR (AU)",   line = list(color = "brown", dash = "dot")) %>%
      layout(title = "VHL/HIF/VEGF/mTOR Pathway Dynamics",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Biomarker level"),
             legend = list(orientation = "h"))
  })

  output$hif_curve <- renderPlotly({
    vhl_seq <- seq(0.05, 1.0, by = 0.05)
    hif_ss  <- (0.18 / 0.35) * (1 / (vhl_seq * 0.12 / 0.08 + 0.01))
    plot_ly(x = vhl_seq, y = hif_ss, type = "scatter", mode = "lines",
            line = list(color = "red")) %>%
      layout(title = "HIF-2α SS vs pVHL Fraction",
             xaxis = list(title = "pVHL fraction"),
             yaxis = list(title = "HIF-2α SS (nM)"))
  })

  output$vegf_hif_scatter <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~HIF2_nM, y = ~VEGF_nM, type = "scatter",
            mode = "markers",
            marker = list(color = ~time_wk, colorscale = "Reds", showscale = TRUE,
                          colorbar = list(title = "Week"))) %>%
      layout(title = "VEGF vs HIF-2α (coloured by time)",
             xaxis = list(title = "HIF-2α (nM)"),
             yaxis = list(title = "VEGF (nM)"))
  })

  ## ── Tumor plots ────────────────────────────────────────────────
  output$tumor_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_wk, y = ~TumorVol, type = "scatter",
            mode = "lines", line = list(color = "crimson", width = 2)) %>%
      add_lines(y = rep(df$TumorVol[1] * 0.7, nrow(df)),
                line = list(color = "grey", dash = "dash"),
                name = "30% reduction (PR threshold)") %>%
      layout(title  = "Tumor Volume Dynamics (Simeoni TGI)",
             xaxis  = list(title = "Time (weeks)"),
             yaxis  = list(title = "Tumor Volume (cm³)"))
  })

  output$waterfall_plot <- renderPlotly({
    regimens <- c("Untreated", "Sunitinib", "Nivo+Ipi",
                  "Pembro+Axitinib", "Cabo+Nivo",
                  "Cabozantinib", "Everolimus", "Belzutifan")
    bor <- c(15, -22, -35, -42, -45, -25, -10, -18)
    cols <- ifelse(bor < -30, "steelblue", ifelse(bor < 0, "skyblue", "salmon"))
    plot_ly(x = regimens, y = bor, type = "bar",
            marker = list(color = cols)) %>%
      add_lines(x = c(0.5, length(regimens) - 0.5),
                y = c(-30, -30), line = list(color = "red", dash = "dash"),
                name = "PR threshold (−30%)") %>%
      layout(title = "Simulated Best Overall Response (%)",
             xaxis = list(title = "Regimen"),
             yaxis = list(title = "Best Change from Baseline (%)"),
             showlegend = FALSE)
  })

  output$tgi_stats <- renderText({
    df <- sim_data()
    v0 <- df$TumorVol[1]
    vmin <- min(df$TumorVol)
    vlast <- tail(df$TumorVol, 1)
    t_nadir <- df$time_wk[which.min(df$TumorVol)]
    sprintf(
      "Baseline tumor volume: %.2f cm³\nNadir volume: %.2f cm³ (at week %.1f)\nFinal volume: %.2f cm³\nBest reduction: %.1f%%\nFinal change: %.1f%%",
      v0, vmin, t_nadir, vlast,
      (1 - vmin / v0) * 100,
      (vlast / v0 - 1) * 100
    )
  })

  ## ── Immune plots ───────────────────────────────────────────────
  output$immune_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_wk) %>%
      add_lines(y = ~CD8_count,  name = "CD8+ T-cells", line = list(color = "steelblue")) %>%
      add_lines(y = ~Treg_count, name = "Tregs",         line = list(color = "orangered", dash = "dash")) %>%
      add_lines(y = ~MDSC_count, name = "MDSCs",          line = list(color = "grey")) %>%
      layout(title = "Immune Cell Dynamics in TME",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Cell Count (cells/µL)"),
             legend = list(orientation = "h"))
  })

  output$pd1_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_wk, y = ~PD1_occ_pct, type = "scatter",
            mode = "lines", line = list(color = "purple")) %>%
      layout(title = "PD-1 Receptor Occupancy",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "PD-1 Occupancy (%)"))
  })

  output$mdsc_vegf_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~VEGF_nM, y = ~MDSC_count, type = "scatter",
            mode = "markers",
            marker = list(color = ~time_wk, colorscale = "Viridis",
                          showscale = TRUE, colorbar = list(title = "Week"))) %>%
      layout(title = "MDSC vs VEGF (coloured by time)",
             xaxis = list(title = "VEGF (nM)"),
             yaxis = list(title = "MDSC count"))
  })

  ## ── Scenario comparison ────────────────────────────────────────
  compare_data <- eventReactive(input$run_compare, {
    lapply(input$sc_select, function(sc) {
      df <- run_scenario(sc, n_weeks = 52, vhl_frac = input$vhl_frac)
      df$regimen <- sc
      df
    }) %>% bind_rows()
  })

  output$compare_tumor <- renderPlotly({
    df <- compare_data()
    plot_ly(df, x = ~time_wk, y = ~TumorVol, color = ~regimen,
            type = "scatter", mode = "lines") %>%
      layout(title = "Tumor Volume: Regimen Comparison",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Tumor Volume (cm³)"),
             legend = list(orientation = "h"))
  })

  output$compare_table <- renderDT({
    df <- compare_data()
    df %>%
      group_by(regimen) %>%
      summarise(
        `Baseline (cm³)`  = round(first(TumorVol), 2),
        `Nadir (cm³)`     = round(min(TumorVol), 2),
        `Final (cm³)`     = round(last(TumorVol), 2),
        `BOR (%)`         = round((1 - min(TumorVol) / first(TumorVol)) * 100, 1),
        `Peak CD8 (cells/µL)` = round(max(CD8_count), 0),
        `PD-1 occ. (%)` = round(last(PD1_occ_pct), 1),
        .groups = "drop"
      ) %>%
      arrange(`Final (cm³)`)
  }, options = list(dom = 't', pageLength = 10), rownames = FALSE)

  ## ── Biomarker dashboard ────────────────────────────────────────
  output$bm_HIF2_box <- renderValueBox({
    df <- sim_data()
    valueBox(round(tail(df$HIF2_nM, 1), 2), "HIF-2α (nM) Final",
             icon = icon("dna"), color = "red")
  })
  output$bm_VEGF_box <- renderValueBox({
    df <- sim_data()
    valueBox(round(tail(df$VEGF_nM, 1), 2), "VEGF (nM) Final",
             icon = icon("tint"), color = "orange")
  })
  output$bm_mTOR_box <- renderValueBox({
    df <- sim_data()
    valueBox(round(tail(df$mTOR_AU, 1), 2), "mTOR Activity (AU) Final",
             icon = icon("cogs"), color = "yellow")
  })
  output$bm_CD8_box <- renderValueBox({
    df <- sim_data()
    valueBox(round(tail(df$CD8_count, 1), 0), "CD8+ T-cells (cells/µL) Final",
             icon = icon("shield-alt"), color = "green")
  })

  output$biomarker_heat <- renderPlotly({
    regimens <- c("Untreated", "Sunitinib", "Nivo+Ipi",
                  "Pembro+Axitinib", "Cabo+Nivo",
                  "Cabozantinib", "Everolimus", "Belzutifan")
    markers <- c("HIF2_nM", "VEGF_nM", "VEGFR2_nM",
                 "mTOR_AU", "CD8_count", "Treg_count",
                 "MDSC_count", "PD1_occ_pct")
    marker_labels <- c("HIF-2α", "VEGF", "VEGFR2",
                       "mTOR", "CD8", "Treg", "MDSC", "PD-1 occ%")

    heat_mat <- sapply(regimens, function(sc) {
      df <- run_scenario(sc, n_weeks = 12, vhl_frac = input$vhl_frac)
      wk12 <- df[which.min(abs(df$time_wk - 12)), ]
      unlist(wk12[, markers])
    })

    # Z-score normalise rows
    z <- t(apply(heat_mat, 1, function(r) (r - mean(r)) / (sd(r) + 1e-9)))

    plot_ly(z = z, x = regimens, y = marker_labels,
            type = "heatmap", colorscale = "RdBu",
            reversescale = TRUE) %>%
      layout(title = "Biomarker Z-score Heatmap at Week 12",
             xaxis = list(title = ""),
             yaxis = list(title = ""))
  })
}

## ------------------------------------------------------------------
## Launch
## ------------------------------------------------------------------
shinyApp(ui, server)
