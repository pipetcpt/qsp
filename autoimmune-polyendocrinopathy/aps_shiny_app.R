## ============================================================
## APS / APECED QSP Shiny Dashboard
## Autoimmune Polyendocrinopathy Syndrome (Type 1 / APECED)
##
## Tabs:
##  1. Patient Profile & Genetics
##  2. Immunology Dashboard (AIRE, T cells, Tregs, AutoAb)
##  3. PK Profiles (Drug concentrations)
##  4. Endocrine Function (Cortisol, PTH, Ca, Glucose, T4, TSH)
##  5. Clinical Endpoints (Organ function %, scores)
##  6. Scenario Comparison
##  7. Biomarker Monitoring (screening panel)
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

## ── mrgsolve model code (embedded) ───────────────────────────
MODEL_CODE <- '
$PARAM
AIRE_mut_sev  = 0.70
k_AutoT_prod  = 0.005
k_AutoT_clear = 0.10
k_Treg_prod   = 0.020
k_Treg_clear  = 0.08
Emax_Treg     = 0.90
EC50_Treg     = 5.0
k_Ab_prod     = 0.003
k_Ab_clear    = 0.003
Adrenal_fn0   = 100.0
k_adren_dest  = 0.002
k_adren_repair= 0.0001
Cortisol_basal= 12.0
k_cort_clear  = 2.4
ACTH_drive    = 1.5
PTG_fn0       = 100.0
k_PTG_dest    = 0.0015
k_PTG_repair  = 0.0001
PTH_basal     = 40.0
k_PTH_clear   = 6.0
Ca_normal     = 9.4
k_Ca_clear    = 0.15
Ca_GI_abs_rate= 1.2
Ca_renal_frac = 0.95
Beta_mass0    = 100.0
k_beta_dest   = 0.003
k_beta_repl   = 0.001
Ins_max       = 600.0
Ins_clear     = 0.10
Glucose_basal = 90.0
G_stim_half   = 110.0
k_G_clear     = 0.05
HGO_basal     = 4.5
InsSens       = 1.0
Thyroid_fn0   = 100.0
k_thy_dest    = 0.0018
k_thy_repair  = 0.0001
FT4_normal    = 1.2
k_T4_clear    = 0.077
TSH_basal     = 2.0
k_TSH_clear   = 0.083
FT4_set       = 1.2
HC_dose       = 0.0
CsA_dose      = 0.0
Aba_dose      = 0.0
JAKi_dose     = 0.0
IC50_CsA      = 150.0
IC50_Aba      = 2.0
IC50_RTX      = 5.0
IC50_JAKi     = 8.0
ka_HC         = 3.5
F_HC          = 0.96
k_HC_clear    = 16.0
Vd_HC         = 35.0
ka_CsA        = 2.0
F_CsA         = 0.30
k_CsA_clear   = 0.20
F_Aba         = 0.79
k_Aba_clear   = 0.017
Vd_Aba        = 0.07
k_RTX_clear   = 0.014
ka_JAKi       = 6.0
F_JAKi        = 0.74
k_JAKi_clear  = 2.8
Vd_JAKi       = 87.0

$CMT
AIRE_func AutoT_pool Treg_pool
AutoAb_adren AutoAb_PTG AutoAb_beta AutoAb_thy
Adrenal_fn Cortisol_c
PTG_fn PTH_plasma Ca_serum
Beta_mass Insulin_p Glucose_p
Thyroid_fn TSH_plasma FT4_plasma
Drug_CsA Drug_Aba Drug_RTX Drug_JAKi Drug_HC

$INIT
AIRE_func=0.3 AutoT_pool=2.0 Treg_pool=15.0
AutoAb_adren=1.0 AutoAb_PTG=1.0 AutoAb_beta=1.0 AutoAb_thy=1.0
Adrenal_fn=100.0 Cortisol_c=12.0
PTG_fn=100.0 PTH_plasma=40.0 Ca_serum=9.4
Beta_mass=100.0 Insulin_p=60.0 Glucose_p=90.0
Thyroid_fn=100.0 TSH_plasma=2.0 FT4_plasma=1.2
Drug_CsA=0 Drug_Aba=0 Drug_RTX=0 Drug_JAKi=0 Drug_HC=0

$MAIN
double AIRE_activity = 1.0 - AIRE_mut_sev;
double E_CsA  = (Drug_CsA  > 0) ? Drug_CsA  / (IC50_CsA  + Drug_CsA)  : 0;
double E_Aba  = (Drug_Aba  > 0) ? Drug_Aba  / (IC50_Aba  + Drug_Aba)  : 0;
double E_RTX  = (Drug_RTX  > 0) ? Drug_RTX  / (IC50_RTX  + Drug_RTX)  : 0;
double E_JAKi = (Drug_JAKi > 0) ? Drug_JAKi / (IC50_JAKi + Drug_JAKi) : 0;
double Immuno_suppress = (1-0.85*E_CsA)*(1-0.80*E_Aba)*(1-0.80*E_RTX)*(1-0.70*E_JAKi);
double Treg_ratio  = Treg_pool / (AutoT_pool + 0.01);
double Treg_effect = Emax_Treg * Treg_ratio / (EC50_Treg + Treg_ratio);

$ODE
dxdt_AIRE_func = 0;
double AutoT_prod = k_AutoT_prod * (1 + 4*(1 - (1.0-AIRE_mut_sev)));
dxdt_AutoT_pool  = AutoT_prod - k_AutoT_clear * AutoT_pool * (1+Treg_effect) * Immuno_suppress;
dxdt_Treg_pool   = k_Treg_prod*(1-AIRE_mut_sev) - k_Treg_clear*Treg_pool;
double Ab_supp = 1 - 0.85*E_RTX;
dxdt_AutoAb_adren = k_Ab_prod*AutoT_pool*Ab_supp - k_Ab_clear*AutoAb_adren;
dxdt_AutoAb_PTG   = k_Ab_prod*AutoT_pool*Ab_supp - k_Ab_clear*AutoAb_PTG;
dxdt_AutoAb_beta  = k_Ab_prod*AutoT_pool*Ab_supp - k_Ab_clear*AutoAb_beta;
dxdt_AutoAb_thy   = k_Ab_prod*AutoT_pool*Ab_supp - k_Ab_clear*AutoAb_thy;
double dAd = -(k_adren_dest*AutoAb_adren*AutoT_pool*(1-E_Aba*0.5)*(1-E_JAKi*0.4))
             + k_adren_repair*Adrenal_fn;
if(Adrenal_fn<=0&&dAd<0) dAd=0; if(Adrenal_fn>=100&&dAd>0) dAd=0;
dxdt_Adrenal_fn = dAd;
dxdt_Cortisol_c = Cortisol_basal*(Adrenal_fn/100.0)*ACTH_drive - k_cort_clear*Cortisol_c;
double dPTG = -(k_PTG_dest*AutoAb_PTG*AutoT_pool*(1-E_Aba*0.5)) + k_PTG_repair*PTG_fn;
if(PTG_fn<=0&&dPTG<0) dPTG=0; if(PTG_fn>=100&&dPTG>0) dPTG=0;
dxdt_PTG_fn = dPTG;
double PTH_stim = 1 + 2.0*(1-Ca_serum/Ca_normal);
if(PTH_stim<0.1) PTH_stim=0.1;
dxdt_PTH_plasma = PTH_basal*(PTG_fn/100.0)*PTH_stim - k_PTH_clear*PTH_plasma;
double PTH_nr = PTH_plasma/PTH_basal;
dxdt_Ca_serum = Ca_GI_abs_rate*PTH_nr + 0.5 - k_Ca_clear*Ca_serum*(1/Ca_renal_frac)/(PTH_nr+0.1);
double dBeta = -(k_beta_dest*AutoAb_beta*AutoT_pool*(1-E_Aba*0.5)*(1-E_CsA*0.6))
               + k_beta_repl*Beta_mass*(1-Beta_mass/100.0);
if(Beta_mass<=0&&dBeta<0) dBeta=0; if(Beta_mass>=100&&dBeta>0) dBeta=0;
dxdt_Beta_mass = dBeta;
double Cortisol_eff = Cortisol_c + Drug_HC;
double InsSens_eff = InsSens*(1-0.3*(1-(Cortisol_eff/Cortisol_basal>-1?1:Cortisol_eff/Cortisol_basal)));
if(InsSens_eff<0.1) InsSens_eff=0.1;
double Ins_sec = Ins_max*(Beta_mass/100.0)*(Glucose_p/(G_stim_half+Glucose_p));
dxdt_Insulin_p = Ins_sec - Ins_clear*Insulin_p;
dxdt_Glucose_p = HGO_basal*(1-0.7*(Insulin_p/(Insulin_p+80)))
                 - k_G_clear*InsSens_eff*Insulin_p*Glucose_p/100.0;
double dThy = -(k_thy_dest*AutoAb_thy*AutoT_pool*(1-E_JAKi*0.4)) + k_thy_repair*Thyroid_fn;
if(Thyroid_fn<=0&&dThy<0) dThy=0; if(Thyroid_fn>=100&&dThy>0) dThy=0;
dxdt_Thyroid_fn = dThy;
double TSH_drive = 1 + 1.5*(TSH_plasma/TSH_basal-1); if(TSH_drive<0.1) TSH_drive=0.1;
dxdt_FT4_plasma = FT4_normal*(Thyroid_fn/100.0)*TSH_drive - k_T4_clear*FT4_plasma;
double TSH_prod = TSH_basal*(FT4_set/(FT4_plasma+0.01))*1.5; if(TSH_prod>50) TSH_prod=50;
dxdt_TSH_plasma = TSH_prod - k_TSH_clear*TSH_plasma;
dxdt_Drug_CsA  = (CsA_dose>0?ka_CsA*F_CsA*CsA_dose*1000/70.0:0) - k_CsA_clear*Drug_CsA;
dxdt_Drug_Aba  = (Aba_dose>0?F_Aba*Aba_dose/(7*70*Vd_Aba):0) - k_Aba_clear*Drug_Aba;
dxdt_Drug_RTX  = -k_RTX_clear*Drug_RTX;
dxdt_Drug_JAKi = (JAKi_dose>0?ka_JAKi*F_JAKi*JAKi_dose/Vd_JAKi:0) - k_JAKi_clear*Drug_JAKi;
dxdt_Drug_HC   = (HC_dose>0?ka_HC*F_HC*HC_dose/Vd_HC:0) - k_HC_clear*Drug_HC;

$TABLE
capture cortisol_total = Cortisol_c + Drug_HC;
capture HbA1c_est = 3.31 + 0.0237*Glucose_p;
capture APS_comps = (cortisol_total<3?1:0)+(Ca_serum<8.0?1:0)+(Glucose_p>200?1:0)+(TSH_plasma>10?1:0);
'

mod_global <- mcode("APS_SHINY", MODEL_CODE, quiet=TRUE)

## ── Helper: run simulation ────────────────────────────────────
run_sim <- function(params_list, sim_years=5) {
  # RTX bolus events (every 6 months if dose > 0)
  rtx_dose <- params_list[["RTX_dose"]]
  ev_base <- ev(time=0, cmt=1, amt=0)  # dummy
  if(!is.null(rtx_dose) && rtx_dose > 0) {
    rtx_times <- seq(0, sim_years*365, by=180)
    ev_rtx <- ev(cmt="Drug_RTX", amt=rtx_dose*1.7, time=rtx_times, rate=-2)
    ev_base <- ev_base + ev_rtx
  }

  do.call(param, c(list(mod_global), params_list)) %>%
    init(AIRE_func = 1 - params_list[["AIRE_mut_sev"]]) %>%
    mrgsim(ev_base, end=sim_years*365, delta=14, obsonly=TRUE) %>%
    as_tibble() %>%
    mutate(Year = time / 365)
}

## ── Shiny UI ─────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "APS/APECED QSP Model", titleWidth = 280),

  dashboardSidebar(
    width = 260,
    sidebarMenu(
      menuItem("Patient Profile", tabName="profile", icon=icon("user-md")),
      menuItem("Immunology", tabName="immuno", icon=icon("shield-virus")),
      menuItem("Drug PK", tabName="pk", icon=icon("pills")),
      menuItem("Endocrine Function", tabName="endocrine", icon=icon("heartbeat")),
      menuItem("Clinical Endpoints", tabName="endpoints", icon=icon("chart-line")),
      menuItem("Scenario Comparison", tabName="scenario", icon=icon("balance-scale")),
      menuItem("Biomarker Monitoring", tabName="biomarker", icon=icon("vial"))
    ),
    hr(),
    tags$div(style="padding:10px; color:#ddd; font-size:11px;",
      strong("About:"), br(),
      "QSP model for Autoimmune Polyendocrinopathy Syndrome (APS Type 1 / APECED)",
      br(), "20-ODE mrgsolve system",
      br(), "Organ systems: Adrenal, PTG, Pancreas, Thyroid"
    )
  ),

  dashboardBody(
    tabItems(
      ## ── TAB 1: Patient Profile ────────────────────────────
      tabItem("profile",
        fluidRow(
          box(title="Genetic & Disease Parameters", status="purple", solidHeader=TRUE, width=4,
            sliderInput("AIRE_mut", "AIRE Mutation Severity (0=none, 1=complete loss)",
                        min=0, max=1, value=0.70, step=0.05),
            selectInput("APS_type", "APS Type",
                        c("APS Type 1 (APECED — AIRE mutation)"="1",
                          "APS Type 2 (Schmidt — HLA-DR3/4)"="2",
                          "APS Type 3 (Thyroid + T1DM)"="3")),
            numericInput("sim_years", "Simulation Duration (years)", value=5, min=1, max=20),
            actionButton("run_btn", "▶ Run Simulation", class="btn-success btn-lg", width="100%")
          ),
          box(title="HLA & Genetics", status="info", solidHeader=TRUE, width=4,
            selectInput("HLA_type", "HLA Haplotype",
                        c("HLA-DR3/DQ2 (High risk)"="DR3DQ2",
                          "HLA-DR4/DQ8 (Moderate risk)"="DR4DQ8",
                          "HLA-DR3/DR4 (Highest risk)"="DR3DR4",
                          "Non-risk HLA"="none")),
            selectInput("AIRE_variant", "AIRE Mutation Type",
                        c("R257X (Finnish founder, severe)"="R257X",
                          "IVS13+1G>C (British founder)"="IVS13",
                          "964del13 (severe)"="del13",
                          "R139X (moderate)"="R139X",
                          "Novel missense (mild)"="novel_mild")),
            checkboxGroupInput("components_present", "Components Present at Diagnosis",
                               c("Chronic Mucocutaneous Candidiasis"="CMC",
                                 "Hypoparathyroidism"="HP",
                                 "Addison's Disease"="AD",
                                 "T1DM"="T1DM", "Hashimoto's"="HASH"))
          ),
          box(title="Patient Characteristics", status="warning", solidHeader=TRUE, width=4,
            numericInput("pt_age", "Patient Age (years)", value=12, min=1, max=70),
            selectInput("pt_sex", "Sex", c("Female"="F","Male"="M")),
            numericInput("pt_weight", "Body Weight (kg)", value=55, min=10, max=150),
            numericInput("pt_height", "Height (cm)", value=160, min=80, max=210),
            verbatimTextOutput("pt_summary")
          )
        ),
        fluidRow(
          box(title="APS1 Disease Criteria & Timeline", status="primary", solidHeader=TRUE, width=12,
            plotlyOutput("disease_timeline", height="300px")
          )
        )
      ),

      ## ── TAB 2: Immunology Dashboard ────────────────────────
      tabItem("immuno",
        fluidRow(
          box(title="Immunosuppressant Settings", status="purple", solidHeader=TRUE, width=3,
            sliderInput("CsA_dose",  "Cyclosporine A (mg/kg/day)", 0, 10, 0, step=0.5),
            sliderInput("Aba_dose",  "Abatacept SC (mg/week)", 0, 250, 0, step=25),
            sliderInput("RTX_dose",  "Rituximab IV (mg/m², q6mo)", 0, 500, 0, step=50),
            sliderInput("JAKi_dose", "Tofacitinib (mg/day)", 0, 20, 0, step=1)
          ),
          tabBox(title="Immune System Dynamics", width=9,
            tabPanel("T Cells & Tregs",
              plotlyOutput("plot_tcells", height="380px")),
            tabPanel("Autoantibodies",
              plotlyOutput("plot_autoAb", height="380px")),
            tabPanel("Treg:AutoT Ratio",
              plotlyOutput("plot_treg_ratio", height="380px"))
          )
        ),
        fluidRow(
          box(title="Drug Effect Summary", status="info", solidHeader=TRUE, width=12,
            tableOutput("drug_effect_table")
          )
        )
      ),

      ## ── TAB 3: Drug PK ─────────────────────────────────────
      tabItem("pk",
        fluidRow(
          box(title="Hormone Replacement PK", status="yellow", solidHeader=TRUE, width=3,
            sliderInput("HC_dose",  "Hydrocortisone (mg/day)", 0, 40, 20, step=2.5),
            sliderInput("FC_dose_ui", "Fludrocortisone (µg/day)", 0, 300, 100, step=25),
            sliderInput("Ca_supp",  "Calcium Supplement (mg/day)", 0, 3000, 1000, step=250),
            sliderInput("LT4_dose", "Levothyroxine (µg/day)", 0, 200, 75, step=12.5)
          ),
          tabBox(title="Plasma Concentration vs Time", width=9,
            tabPanel("HC / Cortisol", plotlyOutput("plot_pk_hc", height="380px")),
            tabPanel("Immunosuppressants", plotlyOutput("plot_pk_immuno", height="380px")),
            tabPanel("Summary PK Table", DTOutput("pk_table"))
          )
        )
      ),

      ## ── TAB 4: Endocrine Function ──────────────────────────
      tabItem("endocrine",
        fluidRow(
          box(title="Reference Ranges", status="success", solidHeader=TRUE, width=3,
            tags$table(class="table table-condensed",
              tags$thead(tags$tr(tags$th("Marker"), tags$th("Normal Range"))),
              tags$tbody(
                tags$tr(tags$td("Cortisol (AM)"), tags$td("8–20 µg/dL")),
                tags$tr(tags$td("PTH"), tags$td("15–65 pg/mL")),
                tags$tr(tags$td("Ca²⁺ (total)"), tags$td("8.5–10.5 mg/dL")),
                tags$tr(tags$td("Glucose (fasting)"), tags$td("70–100 mg/dL")),
                tags$tr(tags$td("HbA1c"), tags$td("<5.7%")),
                tags$tr(tags$td("TSH"), tags$td("0.4–4.0 mIU/L")),
                tags$tr(tags$td("Free T4"), tags$td("0.8–1.8 ng/dL"))
              )
            ),
            hr(),
            verbatimTextOutput("current_values_txt")
          ),
          tabBox(title="Endocrine Biomarkers", width=9,
            tabPanel("HPA Axis (Cortisol/ACTH)",
              plotlyOutput("plot_hpa", height="380px")),
            tabPanel("Parathyroid / Calcium",
              plotlyOutput("plot_ptg", height="380px")),
            tabPanel("Pancreas / Glucose",
              plotlyOutput("plot_pancreas", height="380px")),
            tabPanel("Thyroid (TSH / T4)",
              plotlyOutput("plot_thyroid", height="380px"))
          )
        )
      ),

      ## ── TAB 5: Clinical Endpoints ──────────────────────────
      tabItem("endpoints",
        fluidRow(
          valueBoxOutput("vbox_adrenal", width=3),
          valueBoxOutput("vbox_PTG", width=3),
          valueBoxOutput("vbox_beta", width=3),
          valueBoxOutput("vbox_thyroid", width=3)
        ),
        fluidRow(
          box(title="Organ Function Over Time (%)", status="primary", solidHeader=TRUE, width=8,
            plotlyOutput("plot_organ_fn", height="380px")),
          box(title="Disease Burden Score", status="danger", solidHeader=TRUE, width=4,
            plotlyOutput("plot_aps_score", height="380px"))
        ),
        fluidRow(
          box(title="Endpoint Summary at Simulation End", status="info", solidHeader=TRUE, width=12,
            DTOutput("endpoint_table")
          )
        )
      ),

      ## ── TAB 6: Scenario Comparison ─────────────────────────
      tabItem("scenario",
        fluidRow(
          box(title="Predefined Scenarios", status="purple", solidHeader=TRUE, width=12,
            checkboxGroupInput("scen_select", "Select Scenarios to Compare",
                               inline=TRUE,
                               choices=c(
                                 "1. Natural History (Severe)"="s1",
                                 "2. HRT Only"="s2",
                                 "3. HRT + CsA"="s3",
                                 "4. HRT + Abatacept"="s4",
                                 "5. HRT + Rituximab"="s5",
                                 "6. HRT + JAKi"="s6",
                                 "7. Early Intervention"="s7"
                               ),
                               selected=c("s1","s2","s4"))
          )
        ),
        fluidRow(
          tabBox(title="Comparison Plots", width=12,
            tabPanel("Organ Function", plotlyOutput("comp_organ", height="420px")),
            tabPanel("Autoantibodies", plotlyOutput("comp_ab", height="420px")),
            tabPanel("Metabolic",      plotlyOutput("comp_meta", height="420px")),
            tabPanel("Summary Table",  DTOutput("comp_table"))
          )
        )
      ),

      ## ── TAB 7: Biomarker Monitoring ────────────────────────
      tabItem("biomarker",
        fluidRow(
          box(title="Annual Screening Panel (APS1)", status="warning", solidHeader=TRUE, width=12,
            p("APS1 requires lifelong annual screening for development of additional disease components."),
            fluidRow(
              column(6,
                tags$table(class="table table-bordered table-hover",
                  tags$thead(tags$tr(
                    tags$th("Biomarker"), tags$th("Purpose"), tags$th("Frequency"),
                    tags$th("Action Threshold")
                  )),
                  tags$tbody(
                    tags$tr(tags$td("Anti-21-OH Ab (IgG)"), tags$td("Addison's risk"),
                            tags$td("Annual"), tags$td(">1 U/mL: monitor cortisol")),
                    tags$tr(tags$td("Anti-NALP5 Ab"), tags$td("Hypoparathyroidism risk"),
                            tags$td("Annual"), tags$td(">10 U/mL: measure Ca, PTH")),
                    tags$tr(tags$td("Anti-GAD65 Ab"), tags$td("T1DM risk"),
                            tags$td("Annual"), tags$td(">5 U/mL: OGTT")),
                    tags$tr(tags$td("Anti-TPO / TG Ab"), tags$td("Thyroid disease"),
                            tags$td("Annual"), tags$td(">34 IU/mL: TSH")),
                    tags$tr(tags$td("Anti-IFN-α Ab"), tags$td("APS1 pathognomonic"),
                            tags$td("At diagnosis"), tags$td("Present: confirms APS1")),
                    tags$tr(tags$td("Anti-Parietal cell Ab"), tags$td("Gastric atrophy"),
                            tags$td("2-yearly"), tags$td("Positive: B12 screen")),
                    tags$tr(tags$td("Anti-17α-OH Ab"), tags$td("POI/Hypogonadism"),
                            tags$td("Annual"), tags$td("Positive: LH/FSH/estradiol"))
                  )
                )
              ),
              column(6,
                plotlyOutput("biomarker_radar", height="350px")
              )
            )
          )
        ),
        fluidRow(
          box(title="Autoantibody Dynamics Simulation", status="info", solidHeader=TRUE, width=12,
            plotlyOutput("biomarker_trend", height="380px")
          )
        )
      )
    )
  )
)

## ── Shiny Server ─────────────────────────────────────────────
server <- function(input, output, session) {

  ## Reactive simulation (single scenario from sliders)
  sim_result <- eventReactive(input$run_btn, {
    withProgress(message="Running ODE simulation...", value=0.5, {
      plist <- list(
        AIRE_mut_sev = input$AIRE_mut,
        HC_dose      = input$HC_dose,
        CsA_dose     = input$CsA_dose,
        Aba_dose     = input$Aba_dose,
        RTX_dose     = input$RTX_dose,
        JAKi_dose    = input$JAKi_dose
      )
      run_sim(plist, sim_years=input$sim_years)
    })
  }, ignoreNULL=FALSE)

  ## Pre-run defaults on startup
  observe({
    if(is.null(sim_result())) {
      run_sim(list(AIRE_mut_sev=0.70, HC_dose=20, CsA_dose=0,
                   Aba_dose=0, RTX_dose=0, JAKi_dose=0), sim_years=5)
    }
  })

  ## Patient summary text
  output$pt_summary <- renderText({
    bmi <- input$pt_weight / (input$pt_height/100)^2
    paste0("Age: ", input$pt_age, " yr  |  Sex: ", input$pt_sex, "\n",
           "BMI: ", round(bmi,1), " kg/m²\n",
           "AIRE severity: ", round(input$AIRE_mut*100,0), "% function lost\n",
           "Components present: ", length(input$components_present))
  })

  ## Disease timeline plot
  output$disease_timeline <- renderPlotly({
    df <- sim_result()
    if(is.null(df)) return(NULL)
    df_long <- df %>%
      pivot_longer(c(Adrenal_fn, PTG_fn, Beta_mass, Thyroid_fn),
                   names_to="Organ", values_to="Fn") %>%
      mutate(Organ = recode(Organ,
        "Adrenal_fn"="Adrenal", "PTG_fn"="Parathyroid",
        "Beta_mass"="Beta Cells", "Thyroid_fn"="Thyroid"))
    p <- ggplot(df_long, aes(Year, Fn, color=Organ)) +
      geom_line(size=1.2) +
      geom_hline(yintercept=c(20,50), linetype=c("dashed","dotted"),
                 color=c("red","orange"), alpha=0.6) +
      scale_color_brewer(palette="Set1") +
      labs(title="Organ Function Timeline", x="Years", y="Function (%)") +
      theme_bw(base_size=12)
    ggplotly(p)
  })

  ## T cells & Tregs
  output$plot_tcells <- renderPlotly({
    df <- sim_result(); if(is.null(df)) return(NULL)
    df_long <- df %>%
      pivot_longer(c(AutoT_pool, Treg_pool), names_to="Cell", values_to="Count") %>%
      mutate(Cell = recode(Cell, "AutoT_pool"="Autoreactive T", "Treg_pool"="Regulatory T (Treg)"))
    p <- ggplot(df_long, aes(Year, Count, color=Cell)) +
      geom_line(size=1.3) +
      scale_color_manual(values=c("#E63946","#2196F3")) +
      labs(x="Years", y="Cells/µL", title="T Cell Dynamics") + theme_bw(base_size=12)
    ggplotly(p)
  })

  ## Autoantibodies
  output$plot_autoAb <- renderPlotly({
    df <- sim_result(); if(is.null(df)) return(NULL)
    df_long <- df %>%
      pivot_longer(c(AutoAb_adren,AutoAb_PTG,AutoAb_beta,AutoAb_thy),
                   names_to="Ab", values_to="Titer") %>%
      mutate(Ab = recode(Ab,
        "AutoAb_adren"="Anti-21-OH","AutoAb_PTG"="Anti-NALP5",
        "AutoAb_beta"="Anti-GAD65","AutoAb_thy"="Anti-TPO"))
    p <- ggplot(df_long, aes(Year, Titer, color=Ab)) +
      geom_line(size=1.2) +
      scale_color_manual(values=c("#FF6B35","#004E89","#1A936F","#C9184A")) +
      labs(x="Years", y="Titer (U/mL)", title="Autoantibody Dynamics") + theme_bw(base_size=12)
    ggplotly(p)
  })

  ## Treg:AutoT ratio
  output$plot_treg_ratio <- renderPlotly({
    df <- sim_result(); if(is.null(df)) return(NULL)
    df2 <- df %>% mutate(Ratio = Treg_pool / (AutoT_pool + 0.01))
    p <- ggplot(df2, aes(Year, Ratio)) +
      geom_line(size=1.5, color="#9C27B0") +
      geom_hline(yintercept=5, linetype="dashed", color="green", alpha=0.7) +
      labs(x="Years", y="Treg:AutoT Ratio",
           title="Treg:Autoreactive T Ratio (>5 protective)",
           subtitle="Dashed: protective threshold") + theme_bw(base_size=12)
    ggplotly(p)
  })

  ## Drug effect table
  output$drug_effect_table <- renderTable({
    df <- sim_result(); if(is.null(df)) return(NULL)
    last <- tail(df, 1)
    data.frame(
      Drug = c("Cyclosporine A","Abatacept","Rituximab","Tofacitinib","Hydrocortisone"),
      `Plasma Conc` = c(
        paste0(round(last$Drug_CsA,1)," ng/mL"),
        paste0(round(last$Drug_Aba,2)," µg/mL"),
        paste0(round(last$Drug_RTX,2)," µg/mL"),
        paste0(round(last$Drug_JAKi,1)," ng/mL"),
        paste0(round(last$Drug_HC,2)," µg/dL")
      ),
      `IC50/Target` = c("150 ng/mL","2 µg/mL","5 µg/mL","8 ng/mL","GR agonist"),
      `Mechanism` = c("Calcineurin inhibition","CD28/B7 blockade",
                       "CD20 B cell depletion","JAK1/3 inhibition",
                       "GR activation → anti-inflam")
    )
  })

  ## HC PK plot
  output$plot_pk_hc <- renderPlotly({
    df <- sim_result(); if(is.null(df)) return(NULL)
    df2 <- df %>% mutate(Cortisol_total = cortisol_total)
    p <- ggplot(df2, aes(Year)) +
      geom_line(aes(y=Cortisol_c, color="Endogenous Cortisol"), size=1.2) +
      geom_line(aes(y=Drug_HC,    color="HC (exogenous)"), size=1.2) +
      geom_line(aes(y=Cortisol_total, color="Total Cortisol"), size=1.5, linetype="dashed") +
      geom_hline(yintercept=c(3,8,20), linetype="dotted",
                 color=c("red","orange","green"), alpha=0.7) +
      scale_color_manual(values=c("#E63946","#2196F3","#1B4332")) +
      labs(x="Years", y="Cortisol (µg/dL)",
           title="Cortisol — Endogenous + HC Replacement",
           subtitle="Red=adrenal failure threshold; Green=normal range") + theme_bw(base_size=12)
    ggplotly(p)
  })

  ## Immunosuppressant PK
  output$plot_pk_immuno <- renderPlotly({
    df <- sim_result(); if(is.null(df)) return(NULL)
    df_long <- df %>%
      pivot_longer(c(Drug_CsA, Drug_Aba, Drug_JAKi),
                   names_to="Drug", values_to="Conc") %>%
      mutate(Drug = recode(Drug, "Drug_CsA"="CsA (ng/mL)",
                           "Drug_Aba"="Abatacept (µg/mL×10)",
                           "Drug_JAKi"="Tofacitinib (ng/mL)"))
    p <- ggplot(df_long, aes(Year, Conc, color=Drug)) +
      geom_line(size=1.2) +
      scale_color_brewer(palette="Dark2") +
      labs(x="Years", y="Concentration", title="Immunosuppressant Drug Concentrations") +
      theme_bw(base_size=12)
    ggplotly(p)
  })

  ## PK summary table
  output$pk_table <- renderDT({
    df <- sim_result(); if(is.null(df)) return(NULL)
    df %>%
      filter(time %% 182 < 14) %>%
      select(Year, Drug_HC, Drug_CsA, Drug_Aba, Drug_RTX, Drug_JAKi) %>%
      mutate(across(where(is.numeric), ~round(., 2))) %>%
      rename("Year"=Year, "HC (µg/dL)"=Drug_HC, "CsA (ng/mL)"=Drug_CsA,
             "Aba (µg/mL)"=Drug_Aba, "RTX (µg/mL)"=Drug_RTX, "JAKi (ng/mL)"=Drug_JAKi) %>%
      datatable(options=list(pageLength=10))
  })

  ## HPA axis plot
  output$plot_hpa <- renderPlotly({
    df <- sim_result(); if(is.null(df)) return(NULL)
    p <- ggplot(df, aes(Year)) +
      geom_line(aes(y=cortisol_total, color="Cortisol (µg/dL)"), size=1.3) +
      geom_ribbon(aes(ymin=8, ymax=20), fill="green", alpha=0.1) +
      labs(x="Years", y="Cortisol (µg/dL)",
           title="HPA Axis — Cortisol Dynamics",
           subtitle="Green shading: normal range (8-20 µg/dL)") + theme_bw(base_size=12)
    ggplotly(p)
  })

  ## PTG/Ca plot
  output$plot_ptg <- renderPlotly({
    df <- sim_result(); if(is.null(df)) return(NULL)
    p <- ggplot(df, aes(Year)) +
      geom_line(aes(y=PTH_plasma/10, color="PTH/10 (pg/mL)"), size=1.3) +
      geom_line(aes(y=Ca_serum, color="Ca²⁺ (mg/dL)"), size=1.3) +
      geom_hline(yintercept=c(8.5, 10.5), linetype="dashed", color=c("red","red"), alpha=0.5) +
      geom_ribbon(aes(ymin=8.5, ymax=10.5), fill="blue", alpha=0.08) +
      scale_color_manual(values=c("#E63946","#2196F3")) +
      labs(x="Years", y="Value", title="Parathyroid & Calcium",
           subtitle="Blue shading: normal Ca range (8.5-10.5 mg/dL)") + theme_bw(base_size=12)
    ggplotly(p)
  })

  ## Pancreas/glucose
  output$plot_pancreas <- renderPlotly({
    df <- sim_result(); if(is.null(df)) return(NULL)
    p <- ggplot(df, aes(Year)) +
      geom_line(aes(y=Glucose_p, color="Glucose (mg/dL)"), size=1.3) +
      geom_line(aes(y=HbA1c_est*10, color="HbA1c×10 (%)"), size=1.3, linetype="dashed") +
      geom_hline(yintercept=c(70,180), linetype="dotted", color=c("red","orange"), alpha=0.5) +
      scale_color_manual(values=c("#4CAF50","#FF9800")) +
      labs(x="Years", y="Value", title="Pancreatic Function — Glucose & HbA1c") + theme_bw(base_size=12)
    ggplotly(p)
  })

  ## Thyroid
  output$plot_thyroid <- renderPlotly({
    df <- sim_result(); if(is.null(df)) return(NULL)
    p <- ggplot(df, aes(Year)) +
      geom_line(aes(y=TSH_plasma, color="TSH (mIU/L)"), size=1.3) +
      geom_line(aes(y=FT4_plasma*10, color="FT4×10 (ng/dL)"), size=1.3) +
      geom_ribbon(aes(ymin=0.4, ymax=4.0), fill="green", alpha=0.08) +
      scale_color_manual(values=c("#9C27B0","#FF5722")) +
      labs(x="Years", y="Value", title="Thyroid — TSH & Free T4",
           subtitle="Green: normal TSH range (0.4-4.0 mIU/L)") + theme_bw(base_size=12)
    ggplotly(p)
  })

  ## Value boxes (final time point)
  make_status_box <- function(val, threshold, label, icon_name, above_bad=FALSE) {
    status_ok  <- ifelse(above_bad, val <= threshold, val >= threshold)
    valueBox(paste0(round(val,0),"%"), label, icon=icon(icon_name),
             color=ifelse(status_ok, "green", "red"))
  }

  output$vbox_adrenal  <- renderValueBox({
    df <- sim_result(); val <- ifelse(is.null(df), NA, tail(df$Adrenal_fn,1))
    valueBox(paste0(round(val,0),"%"), "Adrenal Function",
             icon=icon("capsules"), color=ifelse(val > 20,"green","red"))
  })
  output$vbox_PTG <- renderValueBox({
    df <- sim_result(); val <- ifelse(is.null(df), NA, tail(df$PTG_fn,1))
    valueBox(paste0(round(val,0),"%"), "Parathyroid Function",
             icon=icon("skull-crossbones"), color=ifelse(val > 20,"green","red"))
  })
  output$vbox_beta <- renderValueBox({
    df <- sim_result(); val <- ifelse(is.null(df), NA, tail(df$Beta_mass,1))
    valueBox(paste0(round(val,0),"%"), "Beta Cell Mass",
             icon=icon("syringe"), color=ifelse(val > 20,"green","red"))
  })
  output$vbox_thyroid <- renderValueBox({
    df <- sim_result(); val <- ifelse(is.null(df), NA, tail(df$Thyroid_fn,1))
    valueBox(paste0(round(val,0),"%"), "Thyroid Function",
             icon=icon("heartbeat"), color=ifelse(val > 20,"green","red"))
  })

  ## Organ function plot
  output$plot_organ_fn <- renderPlotly({
    df <- sim_result(); if(is.null(df)) return(NULL)
    df_long <- df %>%
      pivot_longer(c(Adrenal_fn, PTG_fn, Beta_mass, Thyroid_fn),
                   names_to="Organ", values_to="Fn") %>%
      mutate(Organ = recode(Organ,
        "Adrenal_fn"="Adrenal","PTG_fn"="Parathyroid",
        "Beta_mass"="Beta Cells","Thyroid_fn"="Thyroid"))
    p <- ggplot(df_long, aes(Year, Fn, color=Organ, fill=Organ)) +
      geom_line(size=1.4) +
      geom_hline(yintercept=20, linetype="dashed", color="red", alpha=0.6) +
      scale_color_manual(values=c("#E63946","#2196F3","#4CAF50","#FF9800")) +
      labs(x="Years", y="Function (%)", title="Target Organ Function") + theme_bw(base_size=12)
    ggplotly(p)
  })

  ## APS disease score
  output$plot_aps_score <- renderPlotly({
    df <- sim_result(); if(is.null(df)) return(NULL)
    p <- ggplot(df, aes(Year, APS_comps)) +
      geom_step(size=1.5, color="#E63946") +
      scale_y_continuous(breaks=0:4) +
      labs(x="Years", y="Active Components", title="APS Disease Burden Score",
           subtitle="0=remission, 4=all four major organs failing") + theme_bw(base_size=12)
    ggplotly(p)
  })

  ## Endpoint table
  output$endpoint_table <- renderDT({
    df <- sim_result(); if(is.null(df)) return(NULL)
    last <- tail(df, 1)
    data.frame(
      Endpoint       = c("Adrenal Function (%)","PTH (pg/mL)","Serum Ca (mg/dL)",
                          "Beta Cell Mass (%)","Fasting Glucose (mg/dL)","HbA1c (est. %)",
                          "Thyroid Function (%)","TSH (mIU/L)","Free T4 (ng/dL)",
                          "Anti-21-OH Ab (U/mL)","Disease Components"),
      Value          = c(round(last$Adrenal_fn,1), round(last$PTH_plasma,1),
                          round(last$Ca_serum,2),  round(last$Beta_mass,1),
                          round(last$Glucose_p,1), round(last$HbA1c_est,2),
                          round(last$Thyroid_fn,1),round(last$TSH_plasma,2),
                          round(last$FT4_plasma,3),round(last$AutoAb_adren,1),
                          last$APS_comps),
      `Normal Range` = c(">80%","15-65","8.5-10.5",">80%","70-100","<5.7%",
                          ">80%","0.4-4.0","0.8-1.8","<1","0")
    ) %>%
      datatable(options=list(pageLength=15, dom="t"))
  })

  ## ── Scenario comparison (pre-defined) ─────────────────────
  scenario_defs <- list(
    s1=list(AIRE_mut_sev=0.90,HC_dose=0, CsA_dose=0,Aba_dose=0,RTX_dose=0,JAKi_dose=0),
    s2=list(AIRE_mut_sev=0.90,HC_dose=20,CsA_dose=0,Aba_dose=0,RTX_dose=0,JAKi_dose=0),
    s3=list(AIRE_mut_sev=0.90,HC_dose=20,CsA_dose=3.5,Aba_dose=0,RTX_dose=0,JAKi_dose=0),
    s4=list(AIRE_mut_sev=0.90,HC_dose=20,CsA_dose=0,Aba_dose=125,RTX_dose=0,JAKi_dose=0),
    s5=list(AIRE_mut_sev=0.90,HC_dose=20,CsA_dose=0,Aba_dose=0,RTX_dose=375,JAKi_dose=0),
    s6=list(AIRE_mut_sev=0.90,HC_dose=20,CsA_dose=0,Aba_dose=0,RTX_dose=0,JAKi_dose=10),
    s7=list(AIRE_mut_sev=0.30,HC_dose=15,CsA_dose=0,Aba_dose=0,RTX_dose=0,JAKi_dose=0)
  )
  scenario_names <- c(
    s1="1. Natural History",s2="2. HRT Only",s3="3. HRT+CsA",
    s4="4. HRT+Abatacept",s5="5. HRT+Rituximab",s6="6. HRT+JAKi",s7="7. Early Interv."
  )
  pal7 <- c("#E63946","#2196F3","#4CAF50","#FF9800","#9C27B0","#00BCD4","#795548")
  names(pal7) <- names(scenario_defs)

  scenario_data <- reactive({
    req(input$scen_select)
    sims <- lapply(input$scen_select, function(s) {
      df <- run_sim(scenario_defs[[s]], sim_years=5)
      df$ScenName <- scenario_names[s]
      df$ScenID   <- s
      df
    })
    bind_rows(sims)
  })

  output$comp_organ <- renderPlotly({
    df <- scenario_data()
    df_long <- df %>%
      pivot_longer(c(Adrenal_fn,PTG_fn,Beta_mass,Thyroid_fn),
                   names_to="Organ",values_to="Fn") %>%
      mutate(Organ=recode(Organ,"Adrenal_fn"="Adrenal","PTG_fn"="Parathyroid",
                          "Beta_mass"="Beta Cells","Thyroid_fn"="Thyroid"))
    p <- ggplot(df_long, aes(Year,Fn,color=ScenName)) +
      geom_line(size=1.1) +
      facet_wrap(~Organ,nrow=2) +
      scale_color_brewer(palette="Dark2") +
      labs(x="Years",y="Function (%)",title="Organ Function — Scenario Comparison") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$comp_ab <- renderPlotly({
    df <- scenario_data()
    df_long <- df %>%
      pivot_longer(c(AutoAb_adren,AutoAb_PTG,AutoAb_beta,AutoAb_thy),
                   names_to="Ab",values_to="Titer") %>%
      mutate(Ab=recode(Ab,"AutoAb_adren"="Anti-21-OH","AutoAb_PTG"="Anti-NALP5",
                       "AutoAb_beta"="Anti-GAD65","AutoAb_thy"="Anti-TPO"))
    p <- ggplot(df_long, aes(Year,Titer,color=ScenName)) +
      geom_line(size=1.1) + facet_wrap(~Ab,nrow=2) +
      scale_color_brewer(palette="Dark2") +
      labs(x="Years",y="Titer (U/mL)",title="Autoantibody Comparison") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$comp_meta <- renderPlotly({
    df <- scenario_data()
    df_long <- df %>%
      pivot_longer(c(cortisol_total,Ca_serum,Glucose_p,FT4_plasma),
                   names_to="Marker",values_to="Val") %>%
      mutate(Marker=recode(Marker,"cortisol_total"="Cortisol (µg/dL)",
                           "Ca_serum"="Calcium (mg/dL)","Glucose_p"="Glucose (mg/dL)",
                           "FT4_plasma"="FT4 (ng/dL)"))
    p <- ggplot(df_long, aes(Year,Val,color=ScenName)) +
      geom_line(size=1.1) + facet_wrap(~Marker,scales="free_y",nrow=2) +
      scale_color_brewer(palette="Dark2") +
      labs(x="Years",y="Value",title="Metabolic Biomarkers — Comparison") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$comp_table <- renderDT({
    df <- scenario_data()
    df %>%
      filter(time == max(time)) %>%
      select(ScenName, Adrenal_fn, PTG_fn, Beta_mass, Thyroid_fn,
             cortisol_total, Ca_serum, Glucose_p, HbA1c_est,
             TSH_plasma, FT4_plasma, APS_comps) %>%
      mutate(across(where(is.numeric),~round(.,2))) %>%
      rename("Scenario"=ScenName,"Adrenal %"=Adrenal_fn,"PTG %"=PTG_fn,
             "Beta %"=Beta_mass,"Thyroid %"=Thyroid_fn,"Cortisol"=cortisol_total,
             "Ca"=Ca_serum,"Glucose"=Glucose_p,"HbA1c"=HbA1c_est,
             "TSH"=TSH_plasma,"FT4"=FT4_plasma,"Components"=APS_comps) %>%
      datatable(options=list(dom="t", pageLength=10))
  })

  ## Biomarker radar
  output$biomarker_radar <- renderPlotly({
    df <- sim_result(); if(is.null(df)) return(NULL)
    last <- tail(df, 1)
    categories <- c("Adrenal","Parathyroid","Pancreas","Thyroid","Immune control","Ca Balance")
    vals <- c(
      pmin(last$Adrenal_fn/100, 1),
      pmin(last$PTG_fn/100, 1),
      pmin(last$Beta_mass/100, 1),
      pmin(last$Thyroid_fn/100, 1),
      pmin(last$Treg_pool/(last$AutoT_pool*5+0.1), 1),
      pmin((last$Ca_serum - 7)/(10.5 - 7), 1)
    )
    plot_ly(
      type="scatterpolar",
      r=c(vals, vals[1]),
      theta=c(categories, categories[1]),
      fill="toself",
      fillcolor="rgba(150,50,200,0.3)",
      line=list(color="purple")
    ) %>%
      layout(polar=list(radialaxis=list(visible=TRUE, range=c(0,1))),
             title="Organ System Status (1.0 = normal)")
  })

  ## Biomarker trend
  output$biomarker_trend <- renderPlotly({
    df <- sim_result(); if(is.null(df)) return(NULL)
    df_long <- df %>%
      pivot_longer(c(AutoAb_adren, AutoAb_PTG, AutoAb_beta, AutoAb_thy),
                   names_to="Ab", values_to="Titer") %>%
      mutate(Ab = recode(Ab, "AutoAb_adren"="Anti-21-OH (Adrenal)",
                         "AutoAb_PTG"="Anti-NALP5 (PTG)",
                         "AutoAb_beta"="Anti-GAD65 (Pancreas)",
                         "AutoAb_thy"="Anti-TPO (Thyroid)"))
    p <- ggplot(df_long, aes(Year, Titer, color=Ab)) +
      geom_line(size=1.3) +
      geom_hline(yintercept=1, linetype="dashed", color="gray50", alpha=0.7) +
      scale_color_manual(values=c("#E63946","#2196F3","#4CAF50","#FF9800")) +
      labs(x="Years", y="Autoantibody Titer (U/mL)",
           title="Biomarker Monitoring — Annual Autoantibody Screening",
           subtitle="Dashed: normal/low baseline (1 U/mL)") + theme_bw(base_size=12)
    ggplotly(p)
  })

  ## Current values text
  output$current_values_txt <- renderText({
    df <- sim_result(); if(is.null(df)) return("Run simulation first")
    last <- tail(df, 1)
    paste0(
      "== Current Values (End of Simulation) ==\n",
      "Cortisol: ",  round(last$cortisol_total,1), " µg/dL\n",
      "PTH:      ",  round(last$PTH_plasma,0),     " pg/mL\n",
      "Ca²⁺:     ",  round(last$Ca_serum,2),       " mg/dL\n",
      "Glucose:  ",  round(last$Glucose_p,0),       " mg/dL\n",
      "HbA1c:    ",  round(last$HbA1c_est,1),       "%\n",
      "TSH:      ",  round(last$TSH_plasma,2),      " mIU/L\n",
      "FT4:      ",  round(last$FT4_plasma,3),      " ng/dL\n",
      "Components: ",last$APS_comps
    )
  })
}

shinyApp(ui, server)
