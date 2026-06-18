##############################################################################
# Dilated Cardiomyopathy (DCM) — Interactive Shiny QSP Dashboard
# 6 Tabs: Patient Profile | PK | PD Key Metrics | Clinical Endpoints |
#         Scenario Comparison | Biomarkers
##############################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)

# ============================================================
# EMBEDDED MODEL CODE (same ODE system as dcm_mrgsolve_model.R)
# ============================================================
code_dcm <- '
$PROB DCM QSP Model (Shiny)

$PARAM
LVEF0=25, LVEDV0=250, BNP0=800, AngII0=50, Aldo0=200,
NE0=600, Fib0=0.20, TGFb0=15, IL6_0=8, GFR0=55, Vol0=7.5, SixMWT0=300,
ON_ENA=1, ON_CAR=1, ON_SPR=1, ON_SAC=0, ON_DAPA=0, ON_IVA=0,
ka_ENA=0.69, CL_ENA=5.0, Vd_ENA=35, F_ENA=0.40,
ka_CAR=0.80, CL_CAR=120, Vd_CAR=800, F_CAR=0.25,
ka_SPR=0.50, CL_SPR=8.0, Vd_SPR=50, F_SPR=0.65,
ka_SAC=1.20, CL_SAC=3.0, Vd_SAC=40, F_SAC=0.60,
ka_DAPA=1.50, CL_DAPA=8.5, Vd_DAPA=75, F_DAPA=0.78,
Emax_ACEi=0.80, EC50_ACEi=0.05, kout_AngII=0.35, ksyn_AngII=17.5,
Emax_ARB=0.85, EC50_ARB=0.08, kout_Aldo=0.20, AngII_Aldo=0.015,
Emax_MRA=0.90, EC50_MRA=0.10, kout_NE=0.40, ksyn_NE=240,
Emax_BB=0.75, EC50_BB=0.05, NE_HR_slope=0.001,
BB_LVEF=0.010, ACEi_LVEF=0.006, ARB_LVEF=0.006,
ARNI_LVEF=0.012, DAPA_LVEF=0.008, IVA_LVEF=0.004,
AngII_LVEF=0.0003, NE_LVEF=0.0002, Fib_LVEF=0.05,
LVEF_max=55, LVEDV_min=120, kout_BNP=0.04, AngII_BNP=0.01, LVEDV_BNP=0.05,
Emax_NEP=0.75, EC50_NEP=0.03, kout_Fib=0.0005, Aldo_Fib=0.0001,
TGFb_Fib=0.0005, MRA_Fib=0.003, DAPA_Fib=0.002, Fib_max=0.60,
kout_TGFb=0.10, AngII_TGFb=0.10, Aldo_TGFb=0.05, MRA_TGFb=0.20,
kout_IL6=0.15, NE_IL6=0.008, Fib_IL6=0.80, BB_IL6=0.20,
kout_GFR=0.005, AngII_GFR=0.005, DAPA_GFR=0.002, CO_GFR=0.10,
kout_Vol=0.08, MRA_Vol=0.15, DAPA_Vol=0.20, ACEi_Vol=0.10,
kout_6MWT=0.0008, LVEF_6MWT=1.8, IL6_6MWT=5.0, Vol_6MWT=20,
HR0=85, IVA_HR=0.25

$CMT ENA_GUT ENA_CENT CAR_GUT CAR_CENT SPR_GUT SPR_CENT SAC_GUT SAC_CENT DAPA_GUT DAPA_CENT
     AngII Aldo NE LVEF LVEDV BNP Fib TGFb IL6 GFR Vol SixMWT

$MAIN
ENA_GUT_0=0; ENA_CENT_0=0; CAR_GUT_0=0; CAR_CENT_0=0;
SPR_GUT_0=0; SPR_CENT_0=0; SAC_GUT_0=0; SAC_CENT_0=0;
DAPA_GUT_0=0; DAPA_CENT_0=0;
AngII_0=AngII0; Aldo_0=Aldo0; NE_0=NE0; LVEF_0=LVEF0;
LVEDV_0=LVEDV0; BNP_0=BNP0; Fib_0=Fib0; TGFb_0=TGFb0;
IL6_0=IL6_0; GFR_0=GFR0; Vol_0=Vol0; SixMWT_0=SixMWT0;

$ODE
double C_ENA=ENA_CENT/Vd_ENA, C_CAR=CAR_CENT/Vd_CAR,
       C_SPR=SPR_CENT/Vd_SPR, C_SAC=SAC_CENT/Vd_SAC, C_DAPA=DAPA_CENT/Vd_DAPA;
dxdt_ENA_GUT=-ka_ENA*ENA_GUT;
dxdt_ENA_CENT=ka_ENA*ENA_GUT*F_ENA-CL_ENA*C_ENA;
dxdt_CAR_GUT=-ka_CAR*CAR_GUT;
dxdt_CAR_CENT=ka_CAR*CAR_GUT*F_CAR-CL_CAR*C_CAR;
dxdt_SPR_GUT=-ka_SPR*SPR_GUT;
dxdt_SPR_CENT=ka_SPR*SPR_GUT*F_SPR-CL_SPR*C_SPR;
dxdt_SAC_GUT=-ka_SAC*SAC_GUT;
dxdt_SAC_CENT=ka_SAC*SAC_GUT*F_SAC-CL_SAC*C_SAC;
dxdt_DAPA_GUT=-ka_DAPA*DAPA_GUT;
dxdt_DAPA_CENT=ka_DAPA*DAPA_GUT*F_DAPA-CL_DAPA*C_DAPA;

double E_ACEi=ON_ENA*Emax_ACEi*C_ENA/(EC50_ACEi+C_ENA);
double C_ARB=ON_SAC*0.5;
double E_ARB=ON_SAC*Emax_ARB*C_ARB/(EC50_ARB+C_ARB);
double E_BB=ON_CAR*Emax_BB*C_CAR/(EC50_BB+C_CAR);
double E_MRA=ON_SPR*Emax_MRA*C_SPR/(EC50_MRA+C_SPR);
double E_NEP=ON_SAC*Emax_NEP*C_SAC/(EC50_NEP+C_SAC);
double E_DAPA=ON_DAPA*C_DAPA/(0.2+C_DAPA);
double E_IVA=ON_IVA*IVA_HR;
double ACEi_or_ARB=fmax(E_ACEi,E_ARB);
double CO_effect=fmax(0.1,1.0-LVEF/LVEF_max);

dxdt_AngII=ksyn_AngII*(1.0-ACEi_or_ARB)-kout_AngII*AngII;
dxdt_Aldo=AngII_Aldo*AngII-kout_Aldo*Aldo;
dxdt_NE=ksyn_NE*CO_effect*(1.0-E_BB*0.5)-kout_NE*NE;

double LVEF_det=AngII_LVEF*AngII+NE_LVEF*NE+Fib_LVEF*Fib+0.001*IL6;
double LVEF_imp=BB_LVEF*E_BB+ACEi_LVEF*E_ACEi+ARB_LVEF*E_ARB+ARNI_LVEF*E_NEP+DAPA_LVEF*E_DAPA+IVA_LVEF*E_IVA;
dxdt_LVEF=(LVEF_imp-LVEF_det)*LVEF*(1.0-LVEF/LVEF_max)*(LVEF>5.0?1.0:0.0);

double LVEDV_inc=AngII*0.01+Aldo*0.005+Vol*5;
double LVEDV_dec=(E_ACEi+E_ARB+E_NEP)*2.0+E_MRA*1.0+E_DAPA*1.5;
dxdt_LVEDV=(LVEDV_inc-LVEDV_dec)*0.001*(LVEDV-LVEDV_min);

double NEP_deg=kout_BNP*(1.0-E_NEP);
dxdt_BNP=(AngII_BNP*AngII+LVEDV_BNP*(LVEDV/100.0))-NEP_deg*BNP;

dxdt_Fib=TGFb_Fib*TGFb+Aldo_Fib*Aldo*(1.0-E_MRA)-(MRA_Fib*E_MRA+DAPA_Fib*E_DAPA+kout_Fib)*Fib;
if(Fib>Fib_max)dxdt_Fib=0.0;
dxdt_TGFb=AngII_TGFb*(AngII/AngII0)+Aldo_TGFb*(Aldo/Aldo0)-(kout_TGFb+MRA_TGFb*E_MRA)*TGFb;
dxdt_IL6=NE_IL6*NE+Fib_IL6*Fib-(kout_IL6+BB_IL6*E_BB*0.1)*IL6;

double CO=((LVEF/100.0)*LVEDV*HR0*(1.0-E_IVA))/60000.0;
dxdt_GFR=CO_GFR*CO-(kout_GFR+AngII_GFR*(AngII/AngII0)*0.01)*GFR+DAPA_GFR*E_DAPA*(1.0-GFR/GFR0);
dxdt_Vol=0.5-(MRA_Vol*E_MRA+DAPA_Vol*E_DAPA+ACEi_Vol*(E_ACEi+E_ARB*0.5)+kout_Vol)*Vol;

double target6=LVEF_6MWT*LVEF-IL6_6MWT*IL6-Vol_6MWT*(Vol-5.0);
dxdt_SixMWT=0.002*(target6-SixMWT);

$TABLE
capture HR=HR0*(1.0-E_IVA)*(1.0+NE_HR_slope*(NE-NE0));
capture CO_Lmin=(LVEF/100.0)*LVEDV*HR/60000.0;
capture NT_proBNP=BNP*6.5;
capture NYHA=(LVEF<15)?4.0:(LVEF<25)?3.0:(LVEF<35)?2.5:2.0;
capture Cena=ENA_CENT/Vd_ENA;
capture Ccar=CAR_CENT/Vd_CAR;
capture Cspr=SPR_CENT/Vd_SPR;
capture Csac=SAC_CENT/Vd_SAC;
capture Cdapa=DAPA_CENT/Vd_DAPA;
'

# Compile once at startup
mod_dcm <- mread_cache("DCM_Shiny", tempdir(), code_dcm)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

run_sim <- function(params_list, duration_wk = 24,
                    on_ena, on_car, on_spr, on_sac, on_dapa, on_iva) {
  dur_h <- duration_wk * 7 * 24

  ev <- ev(time = 0)
  if (on_ena && !on_sac) ev <- ev + ev(amt = 10, cmt = "ENA_GUT", ii = 12, addl = dur_h/12 - 1)
  if (on_car)             ev <- ev + ev(amt = 25, cmt = "CAR_GUT", ii = 12, addl = dur_h/12 - 1)
  if (on_spr)             ev <- ev + ev(amt = 25, cmt = "SPR_GUT", ii = 24, addl = dur_h/24 - 1)
  if (on_sac)             ev <- ev + ev(amt = 97, cmt = "SAC_GUT", ii = 12, addl = dur_h/12 - 1)
  if (on_dapa)            ev <- ev + ev(amt = 10, cmt = "DAPA_GUT", ii = 24, addl = dur_h/24 - 1)

  flags <- c(ON_ENA  = as.numeric(on_ena && !on_sac),
             ON_CAR  = as.numeric(on_car),
             ON_SPR  = as.numeric(on_spr),
             ON_SAC  = as.numeric(on_sac),
             ON_DAPA = as.numeric(on_dapa),
             ON_IVA  = as.numeric(on_iva))

  mod_dcm %>%
    param(c(params_list, flags)) %>%
    mrgsim(events = ev, end = dur_h, delta = 24, obsonly = TRUE) %>%
    as.data.frame() %>%
    mutate(time_week = time / (24 * 7))
}

plot_theme <- theme_bw(base_size = 13) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# ============================================================
# UI
# ============================================================

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "DCM QSP Dashboard", titleWidth = 280),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "profile",   icon = icon("user-md")),
      menuItem("Pharmacokinetics",   tabName = "pk",        icon = icon("pills")),
      menuItem("PD Key Metrics",     tabName = "pd",        icon = icon("heartbeat")),
      menuItem("Clinical Endpoints", tabName = "endpoints", icon = icon("chart-line")),
      menuItem("Scenario Comparison",tabName = "compare",   icon = icon("balance-scale")),
      menuItem("Biomarkers",         tabName = "biomarkers",icon = icon("flask"))
    ),

    hr(),
    h4(" Disease Parameters", style = "padding-left:15px; color:#aaa;"),
    sliderInput("LVEF0",   "Baseline LVEF (%)",    min=10, max=40, value=25, step=1),
    sliderInput("LVEDV0",  "Baseline LVEDV (mL)",  min=150, max=350, value=250, step=10),
    sliderInput("BNP0",    "Baseline BNP (pg/mL)", min=200, max=2000, value=800, step=50),
    sliderInput("dur_wk",  "Simulation (weeks)",   min=4,  max=52, value=24, step=4),

    hr(),
    h4(" Treatment Selection", style = "padding-left:15px; color:#aaa;"),
    checkboxInput("on_ena",  "Enalapril (ACEi)",     value = TRUE),
    checkboxInput("on_sac",  "Sacubitril/Val (ARNI)", value = FALSE),
    checkboxInput("on_car",  "Carvedilol (BB)",       value = TRUE),
    checkboxInput("on_spr",  "Spironolactone (MRA)",  value = TRUE),
    checkboxInput("on_dapa", "Dapagliflozin (SGLT2i)",value = FALSE),
    checkboxInput("on_iva",  "Ivabradine (HCN4i)",    value = FALSE),

    hr(),
    actionButton("run", "Run Simulation", icon = icon("play"),
                 style = "margin:10px; background-color:#3498db; color:white;")
  ),

  dashboardBody(
    tabItems(

      # ---- TAB 1: PATIENT PROFILE ----
      tabItem(tabName = "profile",
        fluidRow(
          box(title = "DCM Pathophysiology Summary", width = 12, status = "primary", solidHeader = TRUE,
            p("Dilated Cardiomyopathy (DCM) is characterized by left ventricular dilatation
               and systolic dysfunction in the absence of abnormal loading conditions or coronary artery disease.
               Prevalence: ~1:250 in the general population. Major causes: genetic (~50%), myocarditis,
               tachycardiomyopathy, alcohol, peripartum, idiopathic."),
            p("The QSP model integrates:"),
            tags$ul(
              tags$li("RAAS: Angiotensin II, Aldosterone, ACE2/Ang(1-7) axis"),
              tags$li("SNS: Norepinephrine toxicity, β-receptor down-regulation"),
              tags$li("Cardiac mechanics: LVEF, LVEDV, stroke volume, cardiac output"),
              tags$li("Fibrosis: TGF-β → myofibroblast → collagen deposition"),
              tags$li("Inflammation: IL-6, TNF-α negative inotropy"),
              tags$li("Natriuretic peptides: BNP/NT-proBNP, neprilysin axis"),
              tags$li("Renal: eGFR, volume homeostasis, SGLT2 transporter")
            )
          )
        ),
        fluidRow(
          valueBoxOutput("vb_lvef",   width = 3),
          valueBoxOutput("vb_bnp",    width = 3),
          valueBoxOutput("vb_nyha",   width = 3),
          valueBoxOutput("vb_6mwt",   width = 3)
        ),
        fluidRow(
          box(title = "Baseline Patient Status", width = 6, status = "info",
            tableOutput("tbl_baseline")
          ),
          box(title = "Clinical Trial Evidence Base", width = 6, status = "info",
            tableOutput("tbl_trials")
          )
        )
      ),

      # ---- TAB 2: PHARMACOKINETICS ----
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Drug Plasma Concentrations Over Time", width = 12, status = "primary", solidHeader = TRUE,
            plotlyOutput("plot_pk", height = "450px")
          )
        ),
        fluidRow(
          box(title = "Drug Effect (Emax Fraction) Over Time", width = 12, status = "info",
            plotlyOutput("plot_pk_emax", height = "350px")
          )
        ),
        fluidRow(
          box(title = "PK Parameter Summary", width = 12,
            DTOutput("tbl_pk_params")
          )
        )
      ),

      # ---- TAB 3: PD KEY METRICS ----
      tabItem(tabName = "pd",
        fluidRow(
          box(title = "LVEF Over Time", width = 6, status = "primary", solidHeader = TRUE,
            plotlyOutput("plot_lvef", height = "350px")),
          box(title = "LV End-Diastolic Volume", width = 6, status = "primary", solidHeader = TRUE,
            plotlyOutput("plot_lvedv", height = "350px"))
        ),
        fluidRow(
          box(title = "Angiotensin II & Aldosterone", width = 6, status = "warning", solidHeader = TRUE,
            plotlyOutput("plot_raas", height = "350px")),
          box(title = "Norepinephrine (SNS Activation)", width = 6, status = "warning", solidHeader = TRUE,
            plotlyOutput("plot_ne", height = "350px"))
        ),
        fluidRow(
          box(title = "Heart Rate & Cardiac Output", width = 6, status = "info",
            plotlyOutput("plot_hr_co", height = "350px")),
          box(title = "Plasma Volume Index", width = 6, status = "info",
            plotlyOutput("plot_vol", height = "350px"))
        )
      ),

      # ---- TAB 4: CLINICAL ENDPOINTS ----
      tabItem(tabName = "endpoints",
        fluidRow(
          box(title = "6-Minute Walk Test (Exercise Capacity)", width = 6, status = "success", solidHeader = TRUE,
            plotlyOutput("plot_6mwt", height = "350px")),
          box(title = "NYHA Functional Class", width = 6, status = "success", solidHeader = TRUE,
            plotlyOutput("plot_nyha", height = "350px"))
        ),
        fluidRow(
          box(title = "BNP & NT-proBNP", width = 6, status = "danger", solidHeader = TRUE,
            plotlyOutput("plot_bnp", height = "350px")),
          box(title = "eGFR (Renal Function)", width = 6, status = "danger", solidHeader = TRUE,
            plotlyOutput("plot_gfr", height = "350px"))
        ),
        fluidRow(
          box(title = "Summary Table at Weeks 4 / 12 / 24", width = 12,
            DTOutput("tbl_endpoints")
          )
        )
      ),

      # ---- TAB 5: SCENARIO COMPARISON ----
      tabItem(tabName = "compare",
        fluidRow(
          box(title = "Multi-Scenario LVEF Comparison (All 5 Scenarios)", width = 12,
              status = "primary", solidHeader = TRUE,
            plotlyOutput("plot_compare_lvef", height = "400px"))
        ),
        fluidRow(
          box(title = "BNP Comparison", width = 6, status = "info",
            plotlyOutput("plot_compare_bnp", height = "350px")),
          box(title = "6MWT Comparison", width = 6, status = "success",
            plotlyOutput("plot_compare_6mwt", height = "350px"))
        ),
        fluidRow(
          box(title = "Fibrosis Comparison", width = 6, status = "warning",
            plotlyOutput("plot_compare_fib", height = "350px")),
          box(title = "Scenario Summary at Week 24", width = 6,
            DTOutput("tbl_compare"))
        )
      ),

      # ---- TAB 6: BIOMARKERS ----
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(title = "TGF-β (Pro-fibrotic Signal)", width = 6, status = "warning", solidHeader = TRUE,
            plotlyOutput("plot_tgfb", height = "350px")),
          box(title = "IL-6 (Systemic Inflammation)", width = 6, status = "danger", solidHeader = TRUE,
            plotlyOutput("plot_il6", height = "350px"))
        ),
        fluidRow(
          box(title = "Myocardial Fibrosis (%)", width = 6, status = "warning",
            plotlyOutput("plot_fib", height = "350px")),
          box(title = "Biomarker Reference Ranges", width = 6,
            DTOutput("tbl_biomarkers"))
        )
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

  # Reactive: run single-scenario simulation
  sim_data <- eventReactive(input$run, ignoreNULL = FALSE, {
    params <- list(LVEF0 = input$LVEF0, LVEDV0 = input$LVEDV0, BNP0 = input$BNP0)
    run_sim(params, input$dur_wk,
            input$on_ena, input$on_car, input$on_spr,
            input$on_sac, input$on_dapa, input$on_iva)
  })

  # Reactive: run all 5 canonical scenarios
  compare_data <- eventReactive(input$run, ignoreNULL = FALSE, {
    base_params <- list(LVEF0 = input$LVEF0, LVEDV0 = input$LVEDV0, BNP0 = input$BNP0)
    scens <- list(
      list(name="1 Placebo",             ena=F, car=F, spr=F, sac=F, dapa=F, iva=F),
      list(name="2 ACEi+BB",             ena=T, car=T, spr=F, sac=F, dapa=F, iva=F),
      list(name="3 ACEi+BB+MRA",         ena=T, car=T, spr=T, sac=F, dapa=F, iva=F),
      list(name="4 ARNI+BB+MRA",         ena=F, car=T, spr=T, sac=T, dapa=F, iva=F),
      list(name="5 ARNI+BB+MRA+SGLT2i",  ena=F, car=T, spr=T, sac=T, dapa=T, iva=F)
    )
    purrr::map_dfr(scens, function(s) {
      run_sim(base_params, input$dur_wk, s$ena, s$car, s$spr, s$sac, s$dapa, s$iva) %>%
        mutate(Scenario = s$name)
    })
  })

  # ---- Value Boxes (Tab 1) ----
  output$vb_lvef <- renderValueBox({
    d <- sim_data(); last <- tail(d, 1)
    valueBox(paste0(round(last$LVEF, 1), "%"), "Final LVEF",
             icon = icon("heart"), color = if (last$LVEF < 35) "red" else "green")
  })
  output$vb_bnp <- renderValueBox({
    d <- sim_data(); last <- tail(d, 1)
    valueBox(paste0(round(last$BNP), " pg/mL"), "Final BNP",
             icon = icon("tint"), color = if (last$BNP > 400) "red" else "yellow")
  })
  output$vb_nyha <- renderValueBox({
    d <- sim_data(); last <- tail(d, 1)
    valueBox(round(last$NYHA, 1), "NYHA Class",
             icon = icon("walking"), color = "blue")
  })
  output$vb_6mwt <- renderValueBox({
    d <- sim_data(); last <- tail(d, 1)
    valueBox(paste0(round(last$SixMWT), " m"), "Final 6MWT",
             icon = icon("shoe-prints"), color = "purple")
  })

  # ---- Tab 1 tables ----
  output$tbl_baseline <- renderTable({
    tibble(
      Parameter = c("Disease", "LVEF", "LVEDV", "BNP", "NYHA Class", "6MWT", "eGFR"),
      Value = c("Dilated Cardiomyopathy",
                paste0(input$LVEF0, "%"),
                paste0(input$LVEDV0, " mL"),
                paste0(input$BNP0, " pg/mL"),
                ifelse(input$LVEF0 < 15, "IV", ifelse(input$LVEF0 < 25, "III", "II-III")),
                paste0(round(input$LVEF0 * 1.8 - 8 * 15 - 20 * 2.5), " m"),
                "55 mL/min/1.73m²")
    )
  })

  output$tbl_trials <- renderTable({
    tibble(
      Trial = c("CONSENSUS (1987)", "COPERNICUS (2001)", "RALES (1999)",
                "PARADIGM-HF (2014)", "SHIFT (2010)", "DAPA-HF (2019)"),
      Drug = c("Enalapril", "Carvedilol", "Spironolactone",
               "Sacubitril/Val", "Ivabradine", "Dapagliflozin"),
      `CV Death RRR` = c("16%", "34%", "30%", "20%", "N/A", "18%"),
      `HHF RRR` = c("26%", "35%", "36%", "21%", "26%", "25%"),
      N = c("253", "2289", "1663", "8442", "6558", "4744")
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # ---- Tab 2: PK ----
  output$plot_pk <- renderPlotly({
    d <- sim_data() %>%
      select(time_week, Cena, Ccar, Cspr, Csac, Cdapa) %>%
      pivot_longer(-time_week, names_to = "Drug", values_to = "Conc") %>%
      mutate(Drug = recode(Drug,
        Cena="Enalaprilat", Ccar="Carvedilol", Cspr="Spironolactone",
        Csac="LBQ657 (Sacubitril)", Cdapa="Dapagliflozin"))
    p <- ggplot(d, aes(x=time_week, y=Conc, color=Drug)) +
      geom_line(size=1.1) +
      labs(x="Week", y="Concentration (µg/mL)", color="Drug") +
      plot_theme
    ggplotly(p)
  })

  output$plot_pk_emax <- renderPlotly({
    d <- sim_data()
    # Approximate effect from concentration (for display)
    d2 <- d %>%
      mutate(
        `ACEi (E_max=0.80)` = 0.80 * Cena / (0.05 + Cena),
        `BB (E_max=0.75)`   = 0.75 * Ccar / (0.05 + Ccar),
        `MRA (E_max=0.90)`  = 0.90 * Cspr / (0.10 + Cspr),
        `NEP-i (E_max=0.75)`= 0.75 * Csac / (0.03 + Csac),
        `SGLT2i`            = Cdapa / (0.2 + Cdapa)
      ) %>%
      select(time_week, `ACEi (E_max=0.80)`, `BB (E_max=0.75)`,
             `MRA (E_max=0.90)`, `NEP-i (E_max=0.75)`, `SGLT2i`) %>%
      pivot_longer(-time_week)
    p <- ggplot(d2, aes(x=time_week, y=value, color=name)) +
      geom_line(size=1.1) +
      ylim(0, 1) +
      labs(x="Week", y="Drug Effect (0-1)", color="Mechanism") +
      plot_theme
    ggplotly(p)
  })

  output$tbl_pk_params <- renderDT({
    tibble(
      Drug = c("Enalapril→Enalaprilat", "Carvedilol", "Spironolactone",
               "Sacubitril (LBQ657)", "Dapagliflozin"),
      `Dose (mg)` = c("10 BID", "25 BID", "25 QD", "97 BID", "10 QD"),
      `Bioavailability` = c("40%", "25%", "65%", "60%", "78%"),
      `CL (L/h)` = c(5.0, 120, 8.0, 3.0, 8.5),
      `Vd (L)` = c(35, 800, 50, 40, 75),
      `t½ (h)` = c("11", "7", "9", "9", "12"),
      `EC50 (µg/mL)` = c(0.05, 0.05, 0.10, 0.03, "—")
    )
  }, options = list(pageLength = 5), rownames = FALSE)

  # ---- Tab 3: PD ----
  mk_plotly <- function(df, x, y, ylab, title, hline=NULL, hlab=NULL) {
    p <- ggplot(df, aes_string(x=x, y=y)) +
      geom_line(color="#2196F3", size=1.2) +
      labs(x="Week", y=ylab, title=title) + plot_theme
    if (!is.null(hline))
      p <- p + geom_hline(yintercept=hline, linetype="dashed", color="red") +
               annotate("text", x=max(df[[x]])*0.8, y=hline*1.05, label=hlab, size=3.5, color="red")
    ggplotly(p)
  }

  output$plot_lvef  <- renderPlotly({ mk_plotly(sim_data(), "time_week", "LVEF", "LVEF (%)", "LVEF", 35, "ICD threshold") })
  output$plot_lvedv <- renderPlotly({ mk_plotly(sim_data(), "time_week", "LVEDV", "LVEDV (mL)", "LV End-Diastolic Volume") })
  output$plot_ne    <- renderPlotly({ mk_plotly(sim_data(), "time_week", "NE", "NE (pg/mL)", "Norepinephrine") })
  output$plot_vol   <- renderPlotly({ mk_plotly(sim_data(), "time_week", "Vol", "Volume (L)", "Plasma Volume Index") })

  output$plot_raas <- renderPlotly({
    d <- sim_data() %>% select(time_week, AngII, Aldo) %>%
      pivot_longer(-time_week)
    p <- ggplot(d, aes(x=time_week, y=value, color=name)) +
      geom_line(size=1.2) +
      labs(x="Week", y="pg/mL", color="Variable", title="RAAS Markers") + plot_theme
    ggplotly(p)
  })

  output$plot_hr_co <- renderPlotly({
    d <- sim_data() %>% select(time_week, HR, CO_Lmin) %>%
      pivot_longer(-time_week)
    p <- ggplot(d, aes(x=time_week, y=value, color=name)) +
      geom_line(size=1.2) +
      facet_wrap(~name, scales="free_y") +
      labs(x="Week", y="Value", color="") + plot_theme
    ggplotly(p)
  })

  # ---- Tab 4: Clinical Endpoints ----
  output$plot_6mwt <- renderPlotly({ mk_plotly(sim_data(), "time_week", "SixMWT", "6MWT (m)", "6-Minute Walk Test", 300, "Severe HF threshold") })
  output$plot_nyha <- renderPlotly({ mk_plotly(sim_data(), "time_week", "NYHA", "NYHA Class", "NYHA Functional Class") })
  output$plot_bnp  <- renderPlotly({ mk_plotly(sim_data(), "time_week", "BNP", "BNP (pg/mL)", "BNP", 400, "Risk threshold") })
  output$plot_gfr  <- renderPlotly({ mk_plotly(sim_data(), "time_week", "GFR", "eGFR (mL/min/1.73m²)", "eGFR") })

  output$tbl_endpoints <- renderDT({
    sim_data() %>%
      filter(time_week %in% c(0, 4, 12, input$dur_wk)) %>%
      group_by(time_week) %>%
      summarise(
        `LVEF (%)` = round(mean(LVEF), 1),
        `BNP (pg/mL)` = round(mean(BNP)),
        `NT-proBNP` = round(mean(NT_proBNP)),
        `6MWT (m)` = round(mean(SixMWT)),
        `HR (bpm)` = round(mean(HR)),
        `NYHA` = round(mean(NYHA), 1),
        `eGFR` = round(mean(GFR)),
        `Fibrosis (%)` = round(mean(Fib) * 100, 1)
      ) %>%
      rename(Week = time_week)
  }, rownames = FALSE)

  # ---- Tab 5: Scenario Comparison ----
  mk_compare <- function(y, ylab) {
    d <- compare_data()
    p <- ggplot(d, aes_string(x="time_week", y=y, color="Scenario")) +
      geom_line(size=1.1) +
      labs(x="Week", y=ylab, color="") + plot_theme
    ggplotly(p)
  }

  output$plot_compare_lvef <- renderPlotly({
    d <- compare_data()
    p <- ggplot(d, aes(x=time_week, y=LVEF, color=Scenario)) +
      geom_line(size=1.3) +
      geom_hline(yintercept=35, linetype="dashed", color="gray50") +
      annotate("text", x=max(d$time_week)*0.7, y=36.5, label="ICD/CRT threshold 35%", size=3.5) +
      labs(x="Week", y="LVEF (%)", title="LVEF: All 5 Treatment Scenarios", color="Regimen") +
      scale_color_brewer(palette="Dark2") + plot_theme
    ggplotly(p)
  })
  output$plot_compare_bnp  <- renderPlotly({ mk_compare("BNP",    "BNP (pg/mL)") })
  output$plot_compare_6mwt <- renderPlotly({ mk_compare("SixMWT", "6MWT (m)") })
  output$plot_compare_fib  <- renderPlotly({
    d <- compare_data() %>% mutate(Fib_pct = Fib * 100)
    p <- ggplot(d, aes(x=time_week, y=Fib_pct, color=Scenario)) +
      geom_line(size=1.1) +
      labs(x="Week", y="Fibrosis (%)", color="") + plot_theme
    ggplotly(p)
  })

  output$tbl_compare <- renderDT({
    compare_data() %>%
      filter(time_week == input$dur_wk) %>%
      group_by(Scenario) %>%
      summarise(
        `LVEF (%)` = round(mean(LVEF), 1),
        `ΔLVEF` = round(mean(LVEF) - input$LVEF0, 1),
        `BNP` = round(mean(BNP)),
        `6MWT (m)` = round(mean(SixMWT)),
        `Fibrosis (%)` = round(mean(Fib) * 100, 1)
      )
  }, rownames = FALSE)

  # ---- Tab 6: Biomarkers ----
  output$plot_tgfb <- renderPlotly({ mk_plotly(sim_data(), "time_week", "TGFb", "TGF-β (pg/mL)", "TGF-β (Pro-fibrotic Driver)") })
  output$plot_il6  <- renderPlotly({ mk_plotly(sim_data(), "time_week", "IL6", "IL-6 (pg/mL)", "IL-6 (Cardiac Depressant)") })
  output$plot_fib  <- renderPlotly({
    d <- sim_data() %>% mutate(Fib_pct = Fib * 100)
    mk_plotly(d, "time_week", "Fib_pct", "Fibrosis (%)", "Myocardial Fibrosis Fraction")
  })

  output$tbl_biomarkers <- renderDT({
    tibble(
      Biomarker = c("BNP", "NT-proBNP", "TGF-β", "IL-6", "Troponin I", "Galectin-3",
                    "PICP", "Aldosterone", "Norepinephrine", "CRP"),
      Normal = c("<100", "<300", "<12", "<3", "<0.04", "<17.8",
                 "<140", "<150", "<500", "<5"),
      Unit = c("pg/mL", "pg/mL", "pg/mL", "pg/mL", "ng/mL", "ng/mL",
               "µg/L", "pg/mL", "pg/mL", "mg/L"),
      `DCM Elevated` = c(">400", ">1000", ">15", ">8", ">0.1", ">35",
                          ">200", ">200", ">600", ">10"),
      Role = c("Wall stress", "Heart failure prognosis", "Pro-fibrotic",
               "Negative inotropy", "Myocyte injury", "Fibrosis/HF predictor",
               "Collagen I synthesis", "Na retention / fibrosis", "SNS activation", "Inflammation")
    )
  }, rownames = FALSE, options = list(pageLength = 10))
}

# ============================================================
# LAUNCH
# ============================================================
shinyApp(ui, server)
