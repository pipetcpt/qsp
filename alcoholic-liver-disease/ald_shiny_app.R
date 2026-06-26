## ============================================================
## Alcoholic Liver Disease (ALD) — Interactive Shiny Dashboard
## 7 Tabs: Patient Profile · Drug PK · Disease Dynamics ·
##         Clinical Endpoints · Scenario Comparison · Biomarkers
## ============================================================
library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

## ---- Inline mrgsolve model ----
ald_code <- '
$PARAM @annotated
k_etoh_elim  : 0.15  : Ethanol elimination (/h)
k_ADH        : 0.9   : ADH metabolism (/h)
k_CYP2E1_bas : 0.05  : CYP2E1 basal activity
k_CYP2E1_ind : 0.15  : CYP2E1 induction
km_etoh_CYP  : 50.0  : Km ethanol CYP2E1
k_AA_clear   : 2.0   : Acetaldehyde clearance
k_ROS_prod   : 0.08  : ROS production
k_ROS_clear  : 0.5   : ROS clearance
GSH0         : 5.0   : Baseline GSH
k_GSH_synth  : 0.3   : GSH synthesis
k_GSH_depl   : 0.08  : GSH depletion by ROS
LPS0         : 1.0   : Baseline LPS
k_LPS_prod   : 0.02  : LPS influx rate
k_LPS_clear  : 0.15  : LPS clearance
k_perm_etoh  : 0.005 : Ethanol→permeability
KC0          : 1.0   : Baseline KC activation
k_KC_act     : 0.3   : KC activation rate
k_KC_res     : 0.05  : KC resolution
EC50_KC_LPS  : 2.0   : LPS EC50 for KC
TNF0         : 1.0   : Baseline TNF
k_TNF_prod   : 0.5   : TNF production
k_TNF_clear  : 0.8   : TNF clearance
IL1B0        : 1.0   : Baseline IL1B
k_IL1B_prod  : 0.4   : IL-1β production
k_IL1B_clear : 0.7   : IL-1β clearance
NEUT0        : 1.0   : Baseline neutrophils
k_neut_rec   : 0.2   : Neutrophil recruitment
k_neut_clear : 0.15  : Neutrophil clearance
EC50_neut    : 3.0   : CXCL8 EC50 neutrophil
k_Hdeath_TNF : 0.04  : TNF hepatocyte death
k_Hdeath_ROS : 0.03  : ROS hepatocyte death
k_Hdeath_neut: 0.02  : Neutrophil hepatocyte death
k_Hregen     : 0.005 : Hepatocyte regeneration
k_Hregen_max : 0.02  : Max GCSF-enhanced regen
ALT0         : 30.0  : Baseline ALT
k_ALT_rel    : 0.8   : ALT release rate
k_ALT_clear  : 0.02  : ALT serum clearance
BILI0        : 1.2   : Baseline bilirubin
k_bili_prod  : 0.012 : Bilirubin production
k_bili_conj  : 0.04  : Bilirubin conjugation
INR0         : 1.0   : Baseline INR
k_clot_synth : 0.03  : Clotting factor synthesis
k_fib        : 0.0003: Fibrosis progression
k_fib_regress: 0.0001: Fibrosis regression
EC50_fib_TGF : 2.0   : TGF EC50 fibrosis
Emax_pred    : 0.75  : Pred NF-kB inhibition Emax
EC50_pred    : 150.0 : Pred EC50 (ng/mL)
PRED_ss      : 0.0   : Prednisolone steady-state (ng/mL)
Emax_NAC     : 0.65  : NAC GSH Emax
EC50_NAC     : 80.0  : NAC EC50
NAC_ss       : 0.0   : NAC steady-state (ug/mL)
Emax_GCSF    : 1.0   : G-CSF neutrophil Emax fold
Emax_pento   : 0.40  : Pentoxifylline TNF Emax
PTX_ss       : 0.0   : PTX steady-state (ng/mL)
EC50_pento   : 600.0 : PTX EC50
Emax_anakin  : 0.70  : Anakinra IL-1β Emax
EC50_anakin  : 2000.0: Anakinra EC50
ANK_ss       : 0.0   : Anakinra steady-state (ng/mL)

$INIT
ETOH=0,AA=1.0,ROS=1.0,GSH=5.0,LPS=1.0,KC=1.0,TNF=1.0,IL1B=1.0,
NEUT=1.0,H=1.0,ALT=30.0,BILI=1.2,INR=1.0,F=0.5

$MAIN
double eff_pred  = Emax_pred  * PRED_ss  / (EC50_pred  + PRED_ss);
double eff_NAC   = Emax_NAC   * NAC_ss   / (EC50_NAC   + NAC_ss);
double eff_pento = Emax_pento * PTX_ss   / (EC50_pento + PTX_ss);
double eff_anakin= Emax_anakin* ANK_ss   / (EC50_anakin + ANK_ss);
double neut_fold = Emax_GCSF;
double ROS_clear = k_ROS_clear * (GSH / GSH0);
double KC_drive  = k_KC_act * LPS / (EC50_KC_LPS + LPS) * (1.0 + 0.3*AA/5.0);
double KC_inhib  = eff_pred;
double TNF_inhib = eff_pred + eff_pento - eff_pred*eff_pento;
double IL1B_inhib= eff_pred + eff_anakin - eff_pred*eff_anakin;
double neut_drive= k_neut_rec * KC / (EC50_neut/3.0 + KC);
double H_safe    = (H > 0.01) ? H : 0.01;
double Hdeath    = H_safe * (k_Hdeath_TNF*TNF + k_Hdeath_ROS*ROS + k_Hdeath_neut*NEUT);
double regen_rate= k_Hregen*(1.0-F/4.0) + k_Hregen_max*(neut_fold-1.0)/1.8;
double fib_prog  = k_fib*KC/(EC50_fib_TGF+KC)*(1.0-H_safe);
double fib_regress=k_fib_regress*H_safe;

$ODE
double AA_prod = k_ADH*ETOH + k_CYP2E1_bas*ETOH*ETOH/(km_etoh_CYP+ETOH);
dxdt_ETOH  = -k_etoh_elim*ETOH - k_ADH*ETOH;
dxdt_AA    = AA_prod - k_AA_clear*AA;
dxdt_ROS   = k_ROS_prod*(1+k_CYP2E1_ind*ETOH/50)+0.05*KC - ROS_clear*ROS;
dxdt_GSH   = k_GSH_synth*(1+eff_NAC) - k_GSH_depl*ROS*GSH - 0.02*AA*GSH;
dxdt_LPS   = LPS0*k_LPS_prod*(1+5*k_perm_etoh*ETOH) - k_LPS_clear*LPS;
dxdt_KC    = KC_drive*(1-KC_inhib) - k_KC_res*KC;
dxdt_TNF   = k_TNF_prod*KC*(1-TNF_inhib) - k_TNF_clear*(TNF-TNF0);
dxdt_IL1B  = k_IL1B_prod*KC*(1-IL1B_inhib) - k_IL1B_clear*(IL1B-IL1B0);
dxdt_NEUT  = neut_drive*neut_fold - k_neut_clear*NEUT;
dxdt_H     = -Hdeath + regen_rate*(1.0-H_safe);
dxdt_ALT   = k_ALT_rel*Hdeath*ALT0*20 - k_ALT_clear*(ALT-ALT0);
dxdt_BILI  = k_bili_prod - k_bili_conj*H_safe*BILI;
dxdt_INR   = (k_clot_synth-k_clot_synth*H_safe)/k_clot_synth*0.01 - 0.002*(INR-INR0)*H_safe;
dxdt_F     = fib_prog - fib_regress;

$TABLE
double MELD  = 3.78*log(BILI+0.01)+11.2*log(INR+0.01)+9.57*log(1.01)+6.43;
double DF    = 4.6*(INR-1.0)*14.0+BILI;
double logit = -3.5+0.18*MELD;
double d90   = 1.0/(1+exp(-logit));
capture ALT_c=ALT; capture BILI_c=BILI; capture INR_c=INR;
capture H_c=H*100; capture GSH_c=GSH; capture ROS_c=ROS;
capture KC_c=KC; capture NEUT_c=NEUT; capture F_c=F;
capture TNF_c=TNF; capture IL1B_c=IL1B; capture MELD_c=MELD;
capture DF_c=DF; capture d90_c=d90*100;
'
mod <- mcode("ALD_Shiny", ald_code, quiet = TRUE)

## ---- Helper: run one scenario ----
run_one <- function(input, scenario_preds, duration_days = 90) {
  tfinal <- duration_days * 24

  m <- param(mod,
    PRED_ss  = scenario_preds$PRED_ss,
    NAC_ss   = scenario_preds$NAC_ss,
    Emax_GCSF= scenario_preds$GCSF_fold,
    PTX_ss   = scenario_preds$PTX_ss,
    ANK_ss   = scenario_preds$ANK_ss,
    k_etoh_elim = scenario_preds$k_etoh_elim
  )

  m <- init(m,
    ETOH = input$etoh_bac,
    AA   = input$etoh_bac / 20,
    ROS  = 1.0 + input$etoh_bac / 40,
    GSH  = max(0.5, 5.0 - input$etoh_bac / 50),
    LPS  = 1.0 + input$etoh_bac / 30,
    KC   = 1.0 + input$severity / 3,
    TNF  = 1.0 + input$severity / 2.5,
    IL1B = 1.0 + input$severity / 2.0,
    NEUT = 1.0 + input$severity / 4,
    H    = max(0.2, 1.0 - input$severity / 8),
    ALT  = input$alt_base,
    BILI = input$bili_base,
    INR  = input$inr_base,
    F    = input$fibrosis_base
  )

  sim <- mrgsim(m, end = tfinal, delta = 4)
  df  <- as.data.frame(sim)
  df$time_days <- df$time / 24
  df
}

## ---- UI ----
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "ALD QSP Model"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",   tabName = "patient",   icon = icon("user")),
      menuItem("Drug PK",           tabName = "pk",        icon = icon("pills")),
      menuItem("Disease Dynamics",  tabName = "disease",   icon = icon("bacterium")),
      menuItem("Clinical Endpoints",tabName = "endpoints", icon = icon("chart-line")),
      menuItem("Scenario Comparison",tabName="scenarios",  icon = icon("balance-scale")),
      menuItem("Biomarkers & Risk", tabName = "biomarkers",icon = icon("flask"))
    ),
    hr(),
    h5("Patient Parameters", style = "padding-left:15px; color:#ccc"),
    sliderInput("alt_base",   "Baseline ALT (IU/L)",    30, 600, 180, 10),
    sliderInput("bili_base",  "Baseline Bilirubin (mg/dL)", 0.5, 30, 10, 0.5),
    sliderInput("inr_base",   "Baseline INR",            0.8, 4.0, 1.8, 0.1),
    sliderInput("fibrosis_base","Fibrosis (Laennec 0–4)", 0, 4, 1.5, 0.5),
    sliderInput("etoh_bac",   "BAC at admission (mg/dL)",  0, 200, 60, 10),
    sliderInput("severity",   "Inflammation Severity (0–8)", 0, 8, 4, 0.5),
    selectInput("scenario",   "Treatment Scenario",
      choices = list(
        "S1: Active Drinking (No Rx)" = "S1",
        "S2: Abstinence Only"         = "S2",
        "S3: Prednisolone 40mg QD"    = "S3",
        "S4: NAC IV (GET protocol)"   = "S4",
        "S5: Prednisolone + NAC"      = "S5",
        "S6: G-CSF 5μg/kg × 5d"      = "S6",
        "S7: Prednisolone + Anakinra" = "S7"
      ), selected = "S3"),
    hr(),
    actionButton("run_btn", "Run Simulation", class = "btn-danger btn-block")
  ),

  dashboardBody(
    tabItems(

      ## ---- TAB 1: Patient Profile ----
      tabItem(tabName = "patient",
        fluidRow(
          box(title = "Alcoholic Liver Disease — Disease Overview", width = 8,
              status = "danger", solidHeader = TRUE,
              p(strong("Alcoholic Liver Disease (ALD)"), "encompasses a spectrum from simple steatosis to alcoholic hepatitis (AH), cirrhosis, and hepatocellular carcinoma (HCC). Severe AH (Maddrey's DF ≥32) carries a 30–50% 90-day mortality without treatment."),
              tags$ul(
                tags$li("Global prevalence: ~3 million ALD-related deaths per year"),
                tags$li("Mechanism: Ethanol → acetaldehyde/ROS → gut dysbiosis → LPS → KC/NLRP3 → TNF-α/IL-1β → hepatocyte death + fibrosis"),
                tags$li("Key scoring: MELD ≥20, Maddrey's DF ≥32, Lille ≥0.45 (non-responder)"),
                tags$li("Current standard: Prednisolone 40mg/d × 28d (STOPAH, Thursz 2015 NEJM)")
              )
          ),
          box(title = "Computed Severity Scores", width = 4,
              status = "warning", solidHeader = TRUE,
              uiOutput("severity_box")
          )
        ),
        fluidRow(
          box(title = "ALD Disease Spectrum", width = 12, status = "danger",
              plotOutput("disease_spectrum_plot", height = "250px"))
        )
      ),

      ## ---- TAB 2: Drug PK ----
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Prednisolone PK (2-cmt, 40 mg QD)", width = 6,
              status = "primary", solidHeader = TRUE,
              plotlyOutput("pk_pred_plot", height = "300px")),
          box(title = "NAC IV Concentration (GET protocol)", width = 6,
              status = "success", solidHeader = TRUE,
              plotlyOutput("pk_nac_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "G-CSF PK (5 μg/kg SC × 5 days)", width = 6,
              status = "info", solidHeader = TRUE,
              plotlyOutput("pk_gcsf_plot", height = "300px")),
          box(title = "Drug PK Parameters", width = 6,
              status = "warning",
              DTOutput("pk_params_table"))
        )
      ),

      ## ---- TAB 3: Disease Dynamics ----
      tabItem(tabName = "disease",
        fluidRow(
          box(title = "Oxidative Stress: ROS vs. GSH", width = 6,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("ros_gsh_plot", height = "300px")),
          box(title = "Gut-Liver Axis: LPS & Kupffer Cells", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("lps_kc_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "Cytokines: TNF-α & IL-1β", width = 6,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("cytokine_plot", height = "300px")),
          box(title = "Neutrophil Infiltration", width = 6,
              status = "info", solidHeader = TRUE,
              plotlyOutput("neut_plot", height = "300px"))
        )
      ),

      ## ---- TAB 4: Clinical Endpoints ----
      tabItem(tabName = "endpoints",
        fluidRow(
          box(title = "MELD Score Trajectory", width = 6,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("meld_plot", height = "300px")),
          box(title = "Serum ALT", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("alt_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "Bilirubin & INR", width = 6,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("bili_inr_plot", height = "300px")),
          box(title = "Hepatocyte Viability & Fibrosis", width = 6,
              status = "info", solidHeader = TRUE,
              plotlyOutput("hep_fib_plot", height = "300px"))
        )
      ),

      ## ---- TAB 5: Scenario Comparison ----
      tabItem(tabName = "scenarios",
        fluidRow(
          box(title = "Run All 7 Scenarios", width = 12,
              status = "primary", solidHeader = TRUE,
              actionButton("run_all_btn", "Run All Scenarios", class = "btn-primary"),
              hr(),
              DTOutput("scenario_comparison_table")
          )
        ),
        fluidRow(
          box(title = "MELD — All Scenarios", width = 6, status = "danger",
              plotlyOutput("meld_compare_plot", height = "350px")),
          box(title = "90-day Mortality Risk", width = 6, status = "warning",
              plotlyOutput("mortality_compare_plot", height = "350px"))
        ),
        fluidRow(
          box(title = "Bilirubin — All Scenarios", width = 6, status = "info",
              plotlyOutput("bili_compare_plot", height = "300px")),
          box(title = "Hepatocyte Viability — All Scenarios", width = 6, status = "success",
              plotlyOutput("hep_compare_plot", height = "300px"))
        )
      ),

      ## ---- TAB 6: Biomarkers & Risk ----
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(title = "Lille Score Predictor", width = 6,
              status = "warning", solidHeader = TRUE,
              p("Lille score ≥0.45 at Day 7 predicts steroid non-responders (90-day mortality >70%)."),
              numericInput("bili_d0",  "Bilirubin Day 0 (mg/dL)", value = 10, min = 0.5, max = 60),
              numericInput("bili_d7",  "Bilirubin Day 7 (mg/dL)", value = 7,  min = 0.5, max = 60),
              numericInput("creat_in", "Creatinine (mg/dL)",      value = 0.9, min = 0.3, max = 10),
              numericInput("pt_base",  "PT ratio (patient/control)", value = 1.3, min = 0.8, max = 5),
              numericInput("age_in",   "Age (years)",             value = 50, min = 18, max = 85),
              numericInput("alb_in",   "Albumin (g/dL)",          value = 2.8, min = 1, max = 5),
              actionButton("calc_lille", "Compute Lille Score", class = "btn-warning"),
              br(), br(),
              uiOutput("lille_result")
          ),
          box(title = "MELD & Maddrey's DF Calculator", width = 6,
              status = "danger", solidHeader = TRUE,
              numericInput("meld_bili",  "Bilirubin (mg/dL)",  value = 10,  min = 0.5, max = 60),
              numericInput("meld_inr",   "INR",                value = 1.8, min = 0.8, max = 10),
              numericInput("meld_creat", "Creatinine (mg/dL)", value = 0.9, min = 0.3, max = 15),
              numericInput("meld_pt",    "PT (patient, sec)",  value = 18,  min = 10,  max = 60),
              numericInput("meld_pt_c",  "PT (control, sec)",  value = 12,  min = 10,  max = 20),
              actionButton("calc_meld", "Compute Scores", class = "btn-danger"),
              br(), br(),
              uiOutput("meld_result")
          )
        ),
        fluidRow(
          box(title = "Oxidative Stress Sensitivity Analysis", width = 6,
              status = "info", solidHeader = TRUE,
              sliderInput("gsh_init", "Initial GSH level (mM)", 0.5, 5, 2, 0.1),
              plotlyOutput("sensitivity_gsh_plot", height = "280px")),
          box(title = "Steroid Responder vs Non-responder", width = 6,
              status = "warning", solidHeader = TRUE,
              sliderInput("resp_thresh", "Responder cutoff (Lille)", 0, 1, 0.45, 0.05),
              plotlyOutput("responder_plot", height = "280px"))
        )
      )
    )
  )
)

## ---- Server ----
server <- function(input, output, session) {

  ## Scenario parameter maps
  sc_params <- list(
    S1 = list(PRED_ss=0,   NAC_ss=0,   GCSF_fold=1.0, PTX_ss=0,   ANK_ss=0,    k_etoh_elim=0.05),
    S2 = list(PRED_ss=0,   NAC_ss=0,   GCSF_fold=1.0, PTX_ss=0,   ANK_ss=0,    k_etoh_elim=0.30),
    S3 = list(PRED_ss=200, NAC_ss=0,   GCSF_fold=1.0, PTX_ss=0,   ANK_ss=0,    k_etoh_elim=0.20),
    S4 = list(PRED_ss=0,   NAC_ss=150, GCSF_fold=1.0, PTX_ss=0,   ANK_ss=0,    k_etoh_elim=0.20),
    S5 = list(PRED_ss=200, NAC_ss=150, GCSF_fold=1.0, PTX_ss=0,   ANK_ss=0,    k_etoh_elim=0.20),
    S6 = list(PRED_ss=0,   NAC_ss=0,   GCSF_fold=1.8, PTX_ss=0,   ANK_ss=0,    k_etoh_elim=0.20),
    S7 = list(PRED_ss=200, NAC_ss=0,   GCSF_fold=1.0, PTX_ss=0,   ANK_ss=3500, k_etoh_elim=0.20)
  )
  sc_labels <- c(S1="Active Drinking",S2="Abstinence",S3="Prednisolone",
                 S4="NAC IV",S5="Pred+NAC",S6="G-CSF",S7="Pred+Anakinra")
  sc_colors <- c(S1="#D32F2F",S2="#1976D2",S3="#7B1FA2",S4="#388E3C",
                 S5="#F57C00",S6="#00796B",S7="#5D4037")

  ## Reactive simulation
  sim_data <- eventReactive(input$run_btn, {
    sp <- sc_params[[input$scenario]]
    run_one(input, sp)
  }, ignoreNULL = FALSE)

  all_sim_data <- eventReactive(input$run_all_btn, {
    withProgress(message = "Running all 7 scenarios...", value = 0, {
      out <- list()
      for (sc in names(sc_params)) {
        incProgress(1/7, detail = paste("Scenario", sc))
        sp  <- sc_params[[sc]]
        df  <- run_one(input, sp)
        df$scenario <- sc
        df$label    <- sc_labels[sc]
        out[[sc]]   <- df
      }
    })
    bind_rows(out)
  }, ignoreNULL = TRUE)

  ## ---- Tab 1: Severity scores ----
  output$severity_box <- renderUI({
    meld_val <- 3.78 * log(input$bili_base + 0.01) +
                11.2 * log(input$inr_base  + 0.01) +
                9.57 * log(1.01) + 6.43
    df_val   <- 4.6  * (input$inr_base - 1.0) * 14 + input$bili_base
    abic_val <- 50 * 0.1 + input$bili_base * 0.08 + input$inr_base * 0.8 + 1.0 * 0.3

    meld_col <- if (meld_val >= 25) "red" else if (meld_val >= 15) "orange" else "green"
    df_col   <- if (df_val   >= 32) "red" else "green"
    tagList(
      tags$p(style = paste0("font-size:16px; color:", meld_col),
             strong("MELD: "), round(meld_val, 1)),
      tags$p(style = paste0("font-size:16px; color:", df_col),
             strong("Maddrey's DF: "), round(df_val, 1),
             if (df_val >= 32) tags$small(" (steroid threshold)") else NULL),
      tags$p(style = "font-size:14px",
             strong("ABIC: "), round(abic_val, 2)),
      hr(),
      tags$p(style = "font-size:11px; color:#888",
             "MELD ≥25: high risk. DF ≥32: prednisolone indicated. ABIC: low<6.71, high≥9")
    )
  })

  output$disease_spectrum_plot <- renderPlot({
    stages <- data.frame(
      stage = factor(c("Normal","Steatosis","Steatohepatitis","Fibrosis","Cirrhosis","HCC"),
                     levels = c("Normal","Steatosis","Steatohepatitis","Fibrosis","Cirrhosis","HCC")),
      meld   = c(6, 8, 14, 18, 24, 30),
      width  = c(1,1,1,1,1,1)
    )
    ggplot(stages, aes(stage, meld, fill = stage)) +
      geom_col() +
      geom_text(aes(label = paste0("MELD~", meld)), vjust = -0.3, size = 3.5) +
      scale_fill_brewer(palette = "Reds") +
      labs(x = "ALD Spectrum", y = "Typical MELD Score",
           title = "ALD Disease Spectrum — Typical MELD Progression") +
      theme_bw(base_size = 11) + theme(legend.position = "none")
  })

  ## ---- Tab 2: Drug PK ----
  output$pk_pred_plot <- renderPlotly({
    t_h <- seq(0, 28 * 24, by = 0.5)
    # 2-compartment oral: simplified Bateman-like for SS
    ka <- 1.5; CL <- 8.5; Vc <- 35; Q <- 15; Vp <- 50; F <- 0.82; D <- 40
    ke <- CL / Vc; k12 <- Q / Vc; k21 <- Q / Vp
    conc <- sapply(t_h %% 24, function(t) {
      D * ka * F / Vc / (ka - ke) * (exp(-ke * t) - exp(-ka * t)) * 1e6 / 312.4
    })  # rough single-dose ng/mL (MW=360 g/mol)
    df_pk <- data.frame(time_h = t_h, conc_ng = pmax(0, conc))
    p <- plot_ly(df_pk, x = ~time_h / 24, y = ~conc_ng, type = "scatter", mode = "lines",
                 line = list(color = "#7B1FA2", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Prednisolone (ng/mL)"),
             shapes = list(list(type = "line", x0 = 0, x1 = 28,
                                y0 = 150, y1 = 150,
                                line = list(dash = "dot", color = "red"))))
    p
  })

  output$pk_nac_plot <- renderPlotly({
    # NAC IV: loading 150 mg/kg over 1h, then 50 mg/kg over 4h, then 100 mg/kg over 16h
    # Plasma ~ug/mL; simplified 1-cmt
    t_h <- seq(0, 24, by = 0.1)
    Vc <- 30; CL <- 12; ke <- CL / Vc
    D1 <- 10500; r1 <- D1 / 1; D2 <- 3500; r2 <- D2 / 4; D3 <- 7000; r3 <- D3 / 16
    conc <- sapply(t_h, function(t) {
      c1 <- if (t <= 1)  r1 / CL * (1 - exp(-ke * t))
            else         r1 / CL * (1 - exp(-ke * 1)) * exp(-ke * (t - 1))
      c2 <- if (t > 1 && t <= 5) r2 / CL * (1 - exp(-ke * (t - 1)))
            else if (t > 5) r2 / CL * (1 - exp(-ke * 4)) * exp(-ke * (t - 5))
            else 0
      c3 <- if (t > 5 && t <= 21) r3 / CL * (1 - exp(-ke * (t - 5)))
            else if (t > 21) r3 / CL * (1 - exp(-ke * 16)) * exp(-ke * (t - 21))
            else 0
      (c1 + c2 + c3) / 1000  # μg/mL
    })
    df_pk <- data.frame(time_h = t_h, conc = pmax(0, conc))
    plot_ly(df_pk, x = ~time_h, y = ~conc, type = "scatter", mode = "lines",
            line = list(color = "#388E3C", width = 2)) %>%
      layout(xaxis = list(title = "Time (h)"),
             yaxis = list(title = "NAC Plasma (μg/mL)"),
             shapes = list(list(type = "line", x0 = 0, x1 = 24,
                                y0 = 80, y1 = 80,
                                line = list(dash = "dot", color = "orange"))))
  })

  output$pk_gcsf_plot <- renderPlotly({
    t_h <- seq(0, 5 * 24, by = 0.5)
    ka <- 0.5; ke <- 0.8 / 4.5
    # 5 daily doses 350 μg → ng/mL (Vc=4.5L)
    conc <- numeric(length(t_h))
    for (d in 0:4) {
      t_off <- t_h - d * 24
      conc <- conc + ifelse(t_off >= 0,
                350 / 4.5 * ka / (ka - ke) * (exp(-ke * t_off) - exp(-ka * t_off)), 0)
    }
    df_pk <- data.frame(time_h = t_h, conc = pmax(0, conc))
    plot_ly(df_pk, x = ~time_h / 24, y = ~conc, type = "scatter", mode = "lines",
            line = list(color = "#00796B", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "G-CSF Plasma (ng/mL)"))
  })

  output$pk_params_table <- renderDT({
    df <- data.frame(
      Drug  = c("Prednisolone","Prednisolone","NAC","G-CSF","Pentoxifylline","Anakinra"),
      Param = c("CL (L/h)","Vc (L)","CL (L/h)","Vc (L)","t½ (h)","t½ (h)"),
      Value = c(8.5, 35, 12, 4.5, 1.6, 4.5),
      Source= c("Bergrem 1983","Frey 1982","Borgström 1986","Kuwabara 1994",
                "Lillibridge 1995","Bresnihan 1998")
    )
    datatable(df, options = list(dom = "t"), rownames = FALSE)
  })

  ## ---- Tabs 3 & 4: Disease dynamics ----
  make_dyn_plot <- function(df, ycol, ylab, color) {
    plot_ly(df, x = ~time_days, y = df[[ycol]], type = "scatter", mode = "lines",
            line = list(color = color, width = 2)) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = ylab))
  }

  output$ros_gsh_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_trace(y = ~ROS_c, name = "ROS (rel.)", type = "scatter", mode = "lines",
                line = list(color = "#E53935", width = 2)) %>%
      add_trace(y = ~GSH_c, name = "GSH (mM)", type = "scatter", mode = "lines",
                line = list(color = "#388E3C", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Level"),
             legend = list(orientation = "h"))
  })

  output$lps_kc_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_trace(y = ~LPS,  name = "LPS (portal, rel)", type = "scatter", mode = "lines",
                line = list(color = "#F57C00", width = 2)) %>%
      add_trace(y = ~KC_c, name = "KC Activation (rel)", type = "scatter", mode = "lines",
                line = list(color = "#7B1FA2", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Level"),
             legend = list(orientation = "h"))
  })

  output$cytokine_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_trace(y = ~TNF_c,  name = "TNF-α (rel)",  type = "scatter", mode = "lines",
                line = list(color = "#D32F2F", width = 2)) %>%
      add_trace(y = ~IL1B_c, name = "IL-1β (rel)", type = "scatter", mode = "lines",
                line = list(color = "#EF6C00", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Level (rel)"),
             legend = list(orientation = "h"))
  })

  output$neut_plot <- renderPlotly({
    df <- sim_data()
    make_dyn_plot(df, "NEUT_c", "Liver Neutrophils (rel.)", "#1565C0")
  })

  output$meld_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days, y = ~MELD_c, type = "scatter", mode = "lines",
            line = list(color = "#D32F2F", width = 2.5)) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "MELD Score"),
             shapes = list(
               list(type="line",x0=0,x1=90,y0=20,y1=20,
                    line=list(dash="dot",color="grey")),
               list(type="line",x0=0,x1=90,y0=25,y1=25,
                    line=list(dash="dash",color="red"))
             ))
  })

  output$alt_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days, y = ~ALT_c, type = "scatter", mode = "lines",
            line = list(color = "#F57C00", width = 2.5)) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "ALT (IU/L)"),
             shapes = list(list(type="line",x0=0,x1=90,y0=40,y1=40,
                                line=list(dash="dot",color="grey"))))
  })

  output$bili_inr_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_trace(y = ~BILI_c, name = "Bilirubin (mg/dL)", type = "scatter", mode = "lines",
                line = list(color = "#E65100", width = 2)) %>%
      add_trace(y = ~INR_c,  name = "INR",              type = "scatter", mode = "lines",
                line = list(color = "#B71C1C", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)"),
             legend = list(orientation = "h"))
  })

  output$hep_fib_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_trace(y = ~H_c,   name = "Viable Hepatocytes (%)", type = "scatter", mode = "lines",
                line = list(color = "#2E7D32", width = 2)) %>%
      add_trace(y = ~F_c * 25, name = "Fibrosis Score × 25", type = "scatter", mode = "lines",
                line = list(color = "#BF360C", width = 2, dash = "dash")) %>%
      layout(xaxis = list(title = "Time (days)"),
             legend = list(orientation = "h"))
  })

  ## ---- Tab 5: Scenario Comparison ----
  output$scenario_comparison_table <- renderDT({
    req(input$run_all_btn)
    df_all <- all_sim_data()
    summary_tbl <- df_all %>%
      filter(time_days >= 89) %>%
      group_by(scenario, label) %>%
      slice_tail(n = 1) %>%
      select(scenario, label, MELD_c, ALT_c, BILI_c, INR_c, H_c, F_c, d90_c) %>%
      rename(Scenario = scenario, Treatment = label,
             MELD = MELD_c, `ALT (IU/L)` = ALT_c,
             `Bilirubin` = BILI_c, INR = INR_c,
             `Viable Hep %` = H_c, `Fibrosis` = F_c,
             `D90 Mortality%` = d90_c) %>%
      mutate(across(where(is.numeric), ~round(., 1)))
    datatable(summary_tbl, options = list(dom = "t"), rownames = FALSE) %>%
      formatStyle("D90 Mortality%", backgroundColor = styleInterval(c(30, 50), c("#C8E6C9","#FFF9C4","#FFCDD2")))
  })

  comp_plot <- function(ycol, ylab) {
    req(input$run_all_btn)
    df_all <- all_sim_data()
    plot_ly(df_all, x = ~time_days, y = df_all[[ycol]],
            color = ~scenario, colors = sc_colors,
            type = "scatter", mode = "lines",
            text = ~label) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = ylab),
             legend = list(orientation = "h"))
  }

  output$meld_compare_plot      <- renderPlotly(comp_plot("MELD_c", "MELD Score"))
  output$mortality_compare_plot <- renderPlotly(comp_plot("d90_c",  "90-day Mortality (%)"))
  output$bili_compare_plot      <- renderPlotly(comp_plot("BILI_c", "Bilirubin (mg/dL)"))
  output$hep_compare_plot       <- renderPlotly(comp_plot("H_c",    "Viable Hepatocytes (%)"))

  ## ---- Tab 6: Biomarkers ----
  observeEvent(input$calc_lille, {
    # Lille score formula (Kim et al. 2010)
    bili_d0 <- input$bili_d0
    bili_d7 <- input$bili_d7
    creat   <- input$creat_in
    pt_r    <- input$pt_base
    age     <- input$age_in
    alb     <- input$alb_in
    lille <- 3.19 - 0.101 * age + 0.147 * alb + 0.0165 * bili_d0 -
             0.206 * (creat > 1.3) - 0.0065 * bili_d0 - 0.0096 * bili_d7
    # Simplified Kim formula
    r <- exp(-lille) / (1 + exp(-lille))
    col <- if (r >= 0.45) "red" else "green"
    resp <- if (r >= 0.45) "Non-responder (high risk)" else "Responder (lower risk)"

    output$lille_result <- renderUI({
      tagList(
        tags$h4(style = paste0("color:", col), paste("Lille Score:", round(r, 3))),
        tags$p(style = paste0("color:", col), strong(resp)),
        tags$p(style = "font-size:11px; color:#888",
               "Lille ≥0.45: consider switching to pentoxifylline or liver transplant evaluation")
      )
    })
  })

  observeEvent(input$calc_meld, {
    meld <- 3.78 * log(input$meld_bili + 0.01) +
            11.2 * log(input$meld_inr  + 0.01) +
            9.57 * log(input$meld_creat+ 0.01) + 6.43
    df_s  <- 4.6 * (input$meld_pt - input$meld_pt_c) + input$meld_bili
    meld  <- max(6, round(meld, 1))
    df_s  <- round(df_s, 1)
    mc <- if (meld >= 25) "red" else if (meld >= 15) "orange" else "green"
    dc <- if (df_s >= 32) "red" else "green"

    output$meld_result <- renderUI({
      tagList(
        tags$h4(style = paste0("color:", mc), paste("MELD:", meld)),
        tags$h4(style = paste0("color:", dc), paste("Maddrey's DF:", df_s)),
        if (df_s >= 32) tags$p("DF ≥32: prednisolone 40mg/day × 28 days indicated (if no contraindication)"),
        if (meld >= 25) tags$p("MELD ≥25: high transplant-free mortality; consider listing")
      )
    })
  })

  output$sensitivity_gsh_plot <- renderPlotly({
    gsh_levels <- c(0.5, 1.0, 2.0, 3.0, 5.0)
    t_days <- seq(0, 30, by = 1)
    df_list <- lapply(gsh_levels, function(g) {
      m_gsh <- init(param(mod, NAC_ss = 0, PRED_ss = 0), GSH = g, H = 0.6,
                    ROS = 3, KC = 3, TNF = 3, IL1B = 4, NEUT = 2,
                    BILI = 8, INR = 1.7, ALT = 150)
      s <- as.data.frame(mrgsim(m_gsh, end = 30 * 24, delta = 24))
      data.frame(time_days = s$time / 24, MELD_c = s$MELD_c, GSH_init = g)
    })
    df_gsh <- bind_rows(df_list)
    plot_ly(df_gsh, x = ~time_days, y = ~MELD_c, color = ~factor(GSH_init),
            type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Day"),
             yaxis = list(title = "MELD Score"),
             colorway = c("#B71C1C","#E64A19","#F57C00","#FBC02D","#388E3C"),
             legend = list(title = list(text = "Init GSH (mM)")))
  })

  output$responder_plot <- renderPlotly({
    t_days <- seq(0, 90, by = 1)
    # Responder: Pred_ss=200; Non-responder: simulate as resistance (lower Emax)
    sims <- lapply(c(0.75, 0.3), function(em) {
      m_r <- param(mod, PRED_ss = 200, Emax_pred = em, k_etoh_elim = 0.20)
      m_r <- init(m_r, H = 0.55, BILI = 10, INR = 1.8, KC = 3.5,
                  TNF = 4, IL1B = 5, NEUT = 3, ALT = 180, F = 1.5)
      s <- as.data.frame(mrgsim(m_r, end = 90 * 24, delta = 24))
      data.frame(time_days = s$time / 24, MELD_c = s$MELD_c,
                 Group = if (em > 0.5) "Responder" else "Non-responder")
    })
    df_resp <- bind_rows(sims)
    plot_ly(df_resp, x = ~time_days, y = ~MELD_c, color = ~Group,
            colors = c("Responder" = "#388E3C", "Non-responder" = "#D32F2F"),
            type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Day"),
             yaxis = list(title = "MELD Score"),
             legend = list(orientation = "h"))
  })
}

shinyApp(ui, server)
