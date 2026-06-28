## ============================================================
## OCD QSP Shiny Interactive Dashboard
## Obsessive-Compulsive Disorder — QSP Model Explorer
## ============================================================
library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

## ============================================================
## MRGSOLVE MODEL (compiled once at startup)
## ============================================================
ocd_code <- '
$PARAM
CL_SSRI=28, V1_SSRI=490, V2_SSRI=910, Q_SSRI=18, KA_SSRI=0.70, F_SSRI=0.44
BBB_k12=0.30, BBB_k21=0.15
SERT_EC50=1.20, SERT_n=1.50, SERT_Emax=1.00
k_5HT_rel=0.08, k_5HT_rup=0.40, k_5HT_deg=0.05, k_5HT_base=1.00, auto_inhib=0.30
k_des=0.012, k_res=0.008
OFC_base=1.00, OFC_5HT_k=0.25, OFC_Glu_k=0.40, OFC_tau=24.0
Caud_base=1.00, Caud_OFC_k=0.35, Caud_5HT_k=0.20, Caud_tau=18.0
Dir_base=1.20, Ind_base=0.80, Dir_D2_k=0.30
Thal_base=1.00, Thal_GPi_k=0.50
BDNF_base=1.00, BDNF_5HT_k=0.15, k_BDNF_syn=0.003, k_BDNF_deg=0.002
D2_EC50=3.00, D2_n=1.20
ERP_kmax=0.40, ERP_onset=0.20
YBOCS0=28.0, YBOCS_OFC=0.35, YBOCS_Caud=0.35, YBOCS_Anx=0.20, YBOCS_tau=72.0
Anx_base=1.00, Anx_OFC_k=0.40, Anx_tau=12.0
SSRI_FLAG=0, CMI_FLAG=0, AUG_FLAG=0, ERP_FLAG=0
CL_CMI=24, V1_CMI=630, KA_CMI=0.90, F_CMI=0.36, SERT_EC50_CMI=0.40, KA_DCMI=0.04, CL_DCMI=8.0
CL_RISP=18, V1_RISP=300, KA_RISP=1.20, F_RISP=0.70
ETA_CL=1, ETA_V=1, ETA_YBOCS0=1

$CMT
AG_SSRI A1_SSRI A2_SSRI A_CNS AG_CMI A1_CMI A_DCMI AG_RISP A1_RISP
SERT_OCC HT5_SYN DES_5HT1 OFC_ACT CAUD_ACT THAL_ACT
DIR_PATH IND_PATH BDNF_LV ERP_EFF YBOCS ANXIETY D2R_OCC

$MAIN
double CL_eff=CL_SSRI*ETA_CL; double V1_eff=V1_SSRI*ETA_V;
double kel_SSRI=CL_eff/V1_eff; double k12_SSRI=Q_SSRI/V1_eff;
double k21_SSRI=Q_SSRI/V2_SSRI;
double Cp_SSRI=1000.0*A1_SSRI/V1_eff;
double Ccns_SSRI=1000.0*A_CNS/(V1_eff*0.2);
double Cp_CMI=1000.0*A1_CMI/V1_CMI; double Cp_DCMI=1000.0*A_DCMI/(V1_CMI*0.8);
double Cp_RISP=1000.0*A1_RISP/V1_RISP;
double SERT_tSSRI=SSRI_FLAG*SERT_Emax*pow(Ccns_SSRI,SERT_n)/(pow(SERT_EC50,SERT_n)+pow(Ccns_SSRI,SERT_n));
double SERT_tCMI=CMI_FLAG*SERT_Emax*(Cp_CMI+2.0*Cp_DCMI)/(SERT_EC50_CMI+(Cp_CMI+2.0*Cp_DCMI));
double SERT_target=fmin(fmax(SERT_tSSRI,SERT_tCMI),SERT_Emax);
double D2R_target=AUG_FLAG*pow(Cp_RISP,D2_n)/(pow(D2_EC50,D2_n)+pow(Cp_RISP,D2_n));
double ERP_norm=ERP_FLAG*ERP_kmax*ERP_EFF;
HT5_SYN_0=k_5HT_base; OFC_ACT_0=OFC_base; CAUD_ACT_0=Caud_base;
THAL_ACT_0=Thal_base; DIR_PATH_0=Dir_base; IND_PATH_0=Ind_base;
BDNF_LV_0=BDNF_base; YBOCS_0=YBOCS0*ETA_YBOCS0; ANXIETY_0=Anx_base;
SERT_OCC_0=0.0; DES_5HT1_0=0.0; ERP_EFF_0=0.0; D2R_OCC_0=0.0;

$ODE
dxdt_AG_SSRI=-KA_SSRI*AG_SSRI;
dxdt_A1_SSRI=SSRI_FLAG*KA_SSRI*AG_SSRI*F_SSRI-(kel_SSRI+k12_SSRI)*A1_SSRI+k21_SSRI*A2_SSRI-BBB_k12*A1_SSRI+BBB_k21*A_CNS;
dxdt_A2_SSRI=k12_SSRI*A1_SSRI-k21_SSRI*A2_SSRI;
dxdt_A_CNS=BBB_k12*A1_SSRI-BBB_k21*A_CNS;
dxdt_AG_CMI=-KA_CMI*AG_CMI;
dxdt_A1_CMI=CMI_FLAG*KA_CMI*AG_CMI*F_CMI-(CL_CMI/V1_CMI)*A1_CMI-KA_DCMI*A1_CMI;
dxdt_A_DCMI=KA_DCMI*A1_CMI-(CL_DCMI/V1_CMI)*A_DCMI;
dxdt_AG_RISP=-KA_RISP*AG_RISP;
dxdt_A1_RISP=AUG_FLAG*KA_RISP*AG_RISP*F_RISP-(CL_RISP/V1_RISP)*A1_RISP;
dxdt_SERT_OCC=0.30*(SERT_target-SERT_OCC);
dxdt_D2R_OCC=0.40*(D2R_target-D2R_OCC);
double rel_5HT=k_5HT_rel*(1.0-auto_inhib*(1.0-DES_5HT1));
double rup_5HT=k_5HT_rup*(1.0-SERT_OCC)*HT5_SYN;
double deg_5HT=k_5HT_deg*HT5_SYN;
dxdt_HT5_SYN=rel_5HT-rup_5HT-deg_5HT;
dxdt_DES_5HT1=k_des*(HT5_SYN-1.0)*(1.0-DES_5HT1)-k_res*DES_5HT1;
double OFC_target=OFC_base*(1.0+Caud_OFC_k*(CAUD_ACT-1.0))*(1.0-OFC_5HT_k*(HT5_SYN-1.0))*(1.0-ERP_norm);
dxdt_OFC_ACT=(OFC_target-OFC_ACT)/OFC_tau;
double Caud_target=Caud_base*(1.0+Caud_OFC_k*(OFC_ACT-1.0))*(1.0-Caud_5HT_k*(HT5_SYN-1.0));
dxdt_CAUD_ACT=(Caud_target-CAUD_ACT)/Caud_tau;
double Dir_target=Dir_base*(1.0-0.15*(HT5_SYN-1.0));
double Ind_target=Ind_base*(1.0+Dir_D2_k*D2R_OCC+0.10*(HT5_SYN-1.0));
dxdt_DIR_PATH=0.05*(Dir_target-DIR_PATH);
dxdt_IND_PATH=0.05*(Ind_target-IND_PATH);
double GPi_act=fmax(0.0,DIR_PATH-IND_PATH);
double Thal_target=Thal_base*(1.0-Thal_GPi_k*GPi_act*0.30);
dxdt_THAL_ACT=(Thal_target-THAL_ACT)/48.0;
dxdt_BDNF_LV=k_BDNF_syn*(HT5_SYN-1.0)*(1.0-BDNF_LV)+k_BDNF_syn*BDNF_5HT_k-k_BDNF_deg*(BDNF_LV-BDNF_base);
dxdt_ERP_EFF=ERP_FLAG*ERP_onset*(ERP_kmax-ERP_EFF)*0.05-0.002*ERP_EFF*(1.0-ERP_FLAG);
double Anx_target=Anx_base*(1.0+Anx_OFC_k*(OFC_ACT-1.0))*(1.0-0.20*(BDNF_LV-1.0));
dxdt_ANXIETY=(Anx_target-ANXIETY)/Anx_tau;
double YBOCS_target=YBOCS0*ETA_YBOCS0*(1.0+YBOCS_OFC*(OFC_ACT-1.0))*(1.0+YBOCS_Caud*(CAUD_ACT-1.0))*(1.0+YBOCS_Anx*(ANXIETY-1.0));
YBOCS_target=fmax(0.0,fmin(40.0,YBOCS_target));
dxdt_YBOCS=(YBOCS_target-YBOCS)/YBOCS_tau;

$TABLE
capture Cp_ng=1000.0*A1_SSRI/(V1_SSRI*ETA_V);
capture Ccns_ng=1000.0*A_CNS/((V1_SSRI*ETA_V)*0.2);
capture Cp_CMI_ng=1000.0*A1_CMI/V1_CMI;
capture Cp_RISP_ng=1000.0*A1_RISP/V1_RISP;
capture SERT_pct=SERT_OCC*100.0;
capture D2R_pct=D2R_OCC*100.0;
capture HT5_norm=HT5_SYN; capture OFC_norm=OFC_ACT; capture CAUD_norm=CAUD_ACT;
capture THAL_norm=THAL_ACT; capture BDNF_norm=BDNF_LV;
capture DIR_norm=DIR_PATH; capture IND_norm=IND_PATH;
capture YBOCS_score=YBOCS;
capture YBOCS_pct_chg=(YBOCS0*ETA_YBOCS0-YBOCS)/(YBOCS0*ETA_YBOCS0)*100.0;
capture Responder=(YBOCS_pct_chg>=35.0)?1.0:0.0;
capture In_remission=(YBOCS<=12.0)?1.0:0.0;
capture Anxiety_norm=ANXIETY;
capture GPi_norm=fmax(0.0,DIR_PATH-IND_PATH);
'

mod <- mcode("ocd_shiny", ocd_code)

## ============================================================
## HELPER: Run simulation given UI inputs
## ============================================================
run_sim <- function(
  YBOCS0_in = 28,
  SSRI_dose = 200, SSRI_on = TRUE,
  CMI_dose  = 0,   CMI_on  = FALSE,
  RISP_dose = 0,   AUG_on  = FALSE,
  ERP_on    = FALSE, ERP_kmax_in = 0.40,
  sim_weeks = 52,
  SERT_EC50_in = 1.20, CL_in = 28
) {
  tt <- seq(0, sim_weeks * 7 * 24, by = 24)

  ev_list <- list()
  if (SSRI_on && SSRI_dose > 0) {
    ev_list[["ssri"]] <- ev(amt = SSRI_dose, cmt = 1,
                             ii = 24, addl = sim_weeks * 7 - 1)
  }
  if (CMI_on && CMI_dose > 0) {
    ev_list[["cmi"]] <- ev(amt = CMI_dose, cmt = 5,
                            ii = 24, addl = sim_weeks * 7 - 1)
  }
  if (AUG_on && RISP_dose > 0) {
    # Risperidone starts at week 12
    start_h <- 12 * 7 * 24
    n_doses  <- max(0, sim_weeks * 7 - 12 * 7)
    if (n_doses > 0)
      ev_list[["risp"]] <- ev(amt = RISP_dose, cmt = 8,
                               time = start_h, ii = 24, addl = n_doses - 1)
  }

  ev_combined <- if (length(ev_list) == 0) NULL else {
    do.call(c, ev_list)
  }

  out <- mod %>%
    param(
      SSRI_FLAG = as.integer(SSRI_on),
      CMI_FLAG  = as.integer(CMI_on),
      AUG_FLAG  = as.integer(AUG_on),
      ERP_FLAG  = as.integer(ERP_on),
      YBOCS0    = YBOCS0_in,
      SERT_EC50 = SERT_EC50_in,
      CL_SSRI   = CL_in,
      ERP_kmax  = ERP_kmax_in
    )

  if (!is.null(ev_combined)) {
    out <- out %>% mrgsim(ev = ev_combined, end = max(tt), delta = 24)
  } else {
    out <- out %>% mrgsim(end = max(tt), delta = 24)
  }

  out %>%
    as.data.frame() %>%
    mutate(Day  = time / 24,
           Week = time / 168)
}

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "OCD QSP Model"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",       tabName = "tab_patient", icon = icon("user")),
      menuItem("PK — Drug Levels",      tabName = "tab_pk",      icon = icon("pills")),
      menuItem("PD — Neurotransmitters",tabName = "tab_pd",      icon = icon("brain")),
      menuItem("CSTC Circuit",          tabName = "tab_cstc",    icon = icon("project-diagram")),
      menuItem("Clinical Endpoints",    tabName = "tab_clinical", icon = icon("chart-line")),
      menuItem("Scenario Comparison",   tabName = "tab_scenario", icon = icon("exchange-alt")),
      menuItem("Biomarkers",            tabName = "tab_biomarker",icon = icon("flask")),
      menuItem("About",                 tabName = "tab_about",   icon = icon("info-circle"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .box-title { font-weight:bold; }
      .skin-blue .main-header .logo { background-color:#1A5276; }
      .skin-blue .main-header .navbar { background-color:#2E86C1; }
    "))),

    tabItems(
      ## --------------------------------------------------------
      ## TAB 1: Patient Profile
      ## --------------------------------------------------------
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient Characteristics", status = "primary",
              solidHeader = TRUE, width = 4,
              sliderInput("YBOCS0", "Baseline Y-BOCS Score",
                          min = 16, max = 40, value = 28, step = 1),
              selectInput("severity", "OCD Severity Category",
                          choices = c("Mild (16–23)" = "mild",
                                      "Moderate (24–31)" = "moderate",
                                      "Severe (32–40)" = "severe"),
                          selected = "moderate"),
              checkboxGroupInput("comorbid", "Comorbidities",
                                 choices = c("Major Depression" = "mdd",
                                             "Anxiety Disorder" = "anx",
                                             "Tics/Tourette" = "tic",
                                             "OC Personality" = "ocp"),
                                 selected = NULL),
              sliderInput("age", "Patient Age (years)", 20, 65, 32, 1),
              selectInput("sex", "Sex", choices = c("Female", "Male"), selected = "Female"),
              numericInput("weight", "Body Weight (kg)", value = 68, min = 40, max = 150)
          ),
          box(title = "Pharmacogenomics", status = "warning",
              solidHeader = TRUE, width = 4,
              selectInput("CYP2D6", "CYP2D6 Metabolizer Status",
                          choices = c("Poor (PM)" = "PM",
                                      "Intermediate (IM)" = "IM",
                                      "Normal (NM)" = "NM",
                                      "Ultra-rapid (UM)" = "UM"),
                          selected = "NM"),
              selectInput("SERT_poly", "SLC6A4 Genotype (5-HTTLPR)",
                          choices = c("S/S (low expression)" = "SS",
                                      "S/L (intermediate)" = "SL",
                                      "L/L (high expression)" = "LL"),
                          selected = "SL"),
              selectInput("COMT_poly", "COMT Val158Met",
                          choices = c("Val/Val (high activity)" = "VV",
                                      "Val/Met (intermediate)" = "VM",
                                      "Met/Met (low activity)" = "MM"),
                          selected = "VM"),
              helpText("Pharmacogenomics adjusts simulated PK/PD parameters.")
          ),
          box(title = "Patient Summary", status = "success",
              solidHeader = TRUE, width = 4,
              tableOutput("patient_summary"),
              hr(),
              h4("OCD Severity Assessment"),
              plotOutput("severity_gauge", height = "200px")
          )
        ),
        fluidRow(
          box(title = "Disease Background: OCD Pathophysiology", status = "info",
              solidHeader = TRUE, width = 12,
              column(6,
                h4("CSTC Circuit in OCD"),
                tags$ul(
                  tags$li(strong("OFC hyperactivity:"), " Generates intrusive thoughts (obsessions)"),
                  tags$li(strong("Caudate hypermetabolism:"), " Maintains urge to perform compulsions"),
                  tags$li(strong("Direct > Indirect pathway:"), " Thalamus released from inhibition → runaway loop"),
                  tags$li(strong("Amygdala hyperactivation:"), " Amplifies threat perception and anxiety")
                ),
                h4("Key Neurotransmitter Imbalances"),
                tags$ul(
                  tags$li("5-HT: hypoactivity at postsynaptic sites (despite SSRIs needed)"),
                  tags$li("DA: elevated striatal DA; imbalanced D1/D2 signaling"),
                  tags$li("Glu: OFC hyperglutamatergic drive onto caudate"),
                  tags$li("GABA: reduced inhibitory interneuron function in striatum")
                )
              ),
              column(6,
                h4("Treatment Targets"),
                tags$table(class = "table table-bordered",
                  tags$thead(tags$tr(tags$th("Target"), tags$th("Drug Class"), tags$th("Effect"))),
                  tags$tbody(
                    tags$tr(tags$td("SERT"), tags$td("SSRI / Clomipramine"), tags$td("↑ Synaptic 5-HT → OFC normalization")),
                    tags$tr(tags$td("D2R"), tags$td("Antipsychotics"), tags$td("Restores direct/indirect pathway balance")),
                    tags$tr(tags$td("NMDA-R"), tags$td("Memantine, Ketamine"), tags$td("Glutamate modulation in OFC/caudate")),
                    tags$tr(tags$td("NMDA Gly site"), tags$td("D-Cycloserine"), tags$td("Augments ERP fear extinction")),
                    tags$tr(tags$td("CSTC circuit"), tags$td("ERP / CBT"), tags$td("Top-down OFC normalization"))
                  )
                )
              )
          )
        )
      ),

      ## --------------------------------------------------------
      ## TAB 2: PK
      ## --------------------------------------------------------
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Pharmacokinetic Parameters", status = "primary",
              solidHeader = TRUE, width = 4,
              h4("SSRI (Sertraline)"),
              sliderInput("SSRI_dose", "Sertraline Dose (mg/day)",
                          25, 400, 200, 25),
              sliderInput("CL_SSRI", "CL (L/h)", 10, 60, 28, 1),
              sliderInput("SERT_EC50", "SERT EC50 (ng/mL)",
                          0.2, 5.0, 1.2, 0.1),
              hr(),
              h4("Clomipramine"),
              checkboxInput("CMI_on", "Use Clomipramine instead of SSRI", FALSE),
              sliderInput("CMI_dose", "Clomipramine Dose (mg/day)",
                          25, 300, 150, 25),
              hr(),
              sliderInput("sim_weeks_pk", "Simulation Duration (weeks)",
                          4, 52, 24, 4)
          ),
          box(title = "Plasma & CNS Drug Concentrations", status = "primary",
              solidHeader = TRUE, width = 8,
              plotlyOutput("pk_plot", height = "400px"),
              hr(),
              plotlyOutput("sert_plot", height = "250px")
          )
        ),
        fluidRow(
          box(title = "PK Summary Table", status = "info",
              solidHeader = TRUE, width = 6,
              DTOutput("pk_table")
          ),
          box(title = "PK Concepts", status = "info",
              solidHeader = TRUE, width = 6,
              h5("Key PK Facts — Sertraline"),
              tags$ul(
                tags$li("t½ ≈ 26 h (steady state in ~5 days)"),
                tags$li("Vd ≈ 7 L/kg (highly lipophilic)"),
                tags$li("Protein binding ≈ 98%"),
                tags$li("CYP2D6/2C19 metabolism"),
                tags$li("CNS:Plasma ratio ≈ 2–3×"),
                tags$li("≥80% SERT occupancy required for OCD efficacy"),
                tags$li("Effective plasma: ~30–200 ng/mL (therapeutic range)")
              ),
              h5("SERT Occupancy & Response"),
              tags$p("Unlike depression (where 50–60% SERT occupancy suffices),
                     OCD typically requires ≥80% SERT occupancy, correlating
                     with higher SSRI doses than those used in MDD.
                     (Zitterl et al. 2008, Neuropsychopharmacology)")
          )
        )
      ),

      ## --------------------------------------------------------
      ## TAB 3: Neurotransmitters (PD)
      ## --------------------------------------------------------
      tabItem(tabName = "tab_pd",
        fluidRow(
          box(title = "Treatment Selection", status = "primary",
              solidHeader = TRUE, width = 3,
              checkboxInput("SSRI_on_pd", "Sertraline", TRUE),
              sliderInput("SSRI_dose_pd", "SSRI Dose (mg/day)", 25, 400, 200, 25),
              checkboxInput("CMI_on_pd", "Clomipramine", FALSE),
              sliderInput("CMI_dose_pd", "CMI Dose (mg/day)", 25, 300, 150, 25),
              checkboxInput("ERP_on_pd", "ERP Therapy", FALSE),
              sliderInput("ERP_kmax_pd", "Max ERP Effect (0–1)", 0.1, 0.6, 0.4, 0.05),
              sliderInput("sim_wk_pd", "Duration (weeks)", 4, 52, 24, 4)
          ),
          box(title = "5-HT Dynamics & SERT Occupancy", status = "danger",
              solidHeader = TRUE, width = 9,
              plotlyOutput("ht5_plot", height = "350px"),
              hr(),
              plotlyOutput("des_plot", height = "200px")
          )
        ),
        fluidRow(
          box(title = "Synaptic 5-HT vs SERT Occupancy Relationship",
              status = "info", solidHeader = TRUE, width = 6,
              plotlyOutput("sert_5ht_scatter", height = "300px")
          ),
          box(title = "5-HT1A Desensitization (Critical for Delayed Response)",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("des_timeline", height = "300px"),
              helpText("The 2–4 week delay in SSRI response is partly explained by
                       progressive desensitization of the 5-HT1A somatodendritic
                       autoreceptor, which then allows sustained 5-HT elevation.")
          )
        )
      ),

      ## --------------------------------------------------------
      ## TAB 4: CSTC Circuit
      ## --------------------------------------------------------
      tabItem(tabName = "tab_cstc",
        fluidRow(
          box(title = "CSTC Circuit Dynamics", status = "primary",
              solidHeader = TRUE, width = 12,
              plotlyOutput("cstc_plot", height = "450px")
          )
        ),
        fluidRow(
          box(title = "Direct vs Indirect Pathway Balance",
              status = "info", solidHeader = TRUE, width = 6,
              plotlyOutput("dir_ind_plot", height = "300px"),
              helpText("OCD = excess direct (Go) pathway → thalamus released from inhibition → OFC loop hyperdrive.")
          ),
          box(title = "BDNF & Neuroplasticity",
              status = "success", solidHeader = TRUE, width = 6,
              plotlyOutput("bdnf_plot", height = "300px"),
              helpText("Chronic 5-HT elevation via SSRIs → ↑BDNF → enhanced neuroplasticity and synaptic remodeling.")
          )
        ),
        fluidRow(
          box(title = "CSTC Circuit Overview", status = "warning",
              solidHeader = TRUE, width = 12,
              column(6,
                h4("Direct Pathway (Go) — HYPERACTIVE in OCD"),
                tags$p("OFC → Caudate (D1 neurons) → GPi/SNr (inhibited) → MD Thalamus (disinhibited) → OFC (loop)"),
                tags$p(style="color:#C0392B; font-weight:bold;",
                       "Result: Thalamus over-activates OFC → compulsive loop")
              ),
              column(6,
                h4("Indirect Pathway (No-Go) — HYPOACTIVE in OCD"),
                tags$p("OFC → Caudate (D2 neurons) → GPe (inhibited) → STN (disinhibited) → GPi/SNr (activated) → Thalamus (inhibited)"),
                tags$p(style="color:#1E8449; font-weight:bold;",
                       "Treatment goal: restore indirect pathway to suppress compulsive loop")
              )
          )
        )
      ),

      ## --------------------------------------------------------
      ## TAB 5: Clinical Endpoints
      ## --------------------------------------------------------
      tabItem(tabName = "tab_clinical",
        fluidRow(
          box(title = "Treatment Settings", status = "primary",
              solidHeader = TRUE, width = 3,
              checkboxInput("SSRI_on_c", "Sertraline", TRUE),
              sliderInput("SSRI_dose_c", "SSRI Dose (mg/day)", 25, 400, 200, 25),
              checkboxInput("CMI_on_c", "Clomipramine", FALSE),
              sliderInput("CMI_dose_c", "CMI Dose (mg/day)", 25, 300, 150, 25),
              checkboxInput("AUG_on_c", "Risperidone Augmentation", FALSE),
              sliderInput("RISP_dose_c", "Risperidone Dose (mg/day)", 0.5, 3.0, 1.0, 0.5),
              checkboxInput("ERP_on_c", "ERP Therapy", FALSE),
              sliderInput("YBOCS0_c", "Baseline Y-BOCS", 16, 40, 28, 1),
              sliderInput("sim_wk_c", "Duration (weeks)", 8, 104, 52, 4)
          ),
          box(title = "Y-BOCS Score Over Time", status = "danger",
              solidHeader = TRUE, width = 9,
              plotlyOutput("ybocs_plot", height = "400px"),
              fluidRow(
                valueBoxOutput("vbox_ybocs_final", width = 3),
                valueBoxOutput("vbox_pct_chg", width = 3),
                valueBoxOutput("vbox_responder", width = 3),
                valueBoxOutput("vbox_remission", width = 3)
              )
          )
        ),
        fluidRow(
          box(title = "Anxiety State Over Time", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("anxiety_plot", height = "300px")
          ),
          box(title = "Response & Remission Timeline",
              status = "success", solidHeader = TRUE, width = 6,
              plotlyOutput("response_plot", height = "300px"),
              helpText("Response = ≥35% Y-BOCS reduction from baseline. Remission = Y-BOCS ≤12.")
          )
        )
      ),

      ## --------------------------------------------------------
      ## TAB 6: Scenario Comparison
      ## --------------------------------------------------------
      tabItem(tabName = "tab_scenario",
        fluidRow(
          box(title = "Multi-Scenario Comparison",
              status = "primary", solidHeader = TRUE, width = 12,
              column(3,
                h4("Simulation Settings"),
                sliderInput("sc_YBOCS0", "Baseline Y-BOCS", 16, 40, 28, 1),
                sliderInput("sc_weeks", "Duration (weeks)", 12, 104, 52, 4),
                hr(),
                h4("SSRI Dose Range"),
                sliderInput("sc_ssri_dose", "Sertraline (mg/day)", 25, 400, 200, 25),
                hr(),
                checkboxInput("sc_include_untreated", "Show Untreated", TRUE),
                checkboxInput("sc_include_ssri", "Show SSRI", TRUE),
                checkboxInput("sc_include_cmi", "Show Clomipramine", TRUE),
                checkboxInput("sc_include_aug", "Show SSRI+Augmentation", TRUE),
                checkboxInput("sc_include_combo", "Show SSRI+ERP", TRUE),
                checkboxInput("sc_include_erp", "Show ERP Alone", FALSE)
              ),
              column(9,
                plotlyOutput("scenario_ybocs", height = "350px"),
                hr(),
                plotlyOutput("scenario_sert", height = "250px")
              )
          )
        ),
        fluidRow(
          box(title = "Outcome Table at Week 12, 24, 52",
              status = "info", solidHeader = TRUE, width = 12,
              DTOutput("scenario_table")
          )
        )
      ),

      ## --------------------------------------------------------
      ## TAB 7: Biomarkers
      ## --------------------------------------------------------
      tabItem(tabName = "tab_biomarker",
        fluidRow(
          box(title = "Neuroimaging & Biochemical Biomarkers",
              status = "primary", solidHeader = TRUE, width = 12,
              tabsetPanel(
                tabPanel("Circuit Biomarkers",
                  plotlyOutput("biomarker_circuit", height = "400px"),
                  hr(),
                  helpText("OFC and caudate hypermetabolism (FDG-PET) normalizes with successful treatment.")
                ),
                tabPanel("Neurochemical Biomarkers",
                  plotlyOutput("biomarker_neuro", height = "400px"),
                  helpText("SERT occupancy (PET), 5-HIAA (CSF), BDNF (serum) as treatment response biomarkers.")
                ),
                tabPanel("Pharmacogenomics Impact",
                  fluidRow(
                    column(6, plotlyOutput("pgx_sert", height = "300px")),
                    column(6, plotlyOutput("pgx_ybocs", height = "300px"))
                  ),
                  helpText("CYP2D6 genotype significantly affects SSRI plasma levels and thus SERT occupancy and efficacy.")
                ),
                tabPanel("DBS Simulation",
                  sliderInput("dbs_effect", "DBS Effect on CSTC Circuit (0=off, 1=max)",
                              0, 1, 0, 0.1),
                  plotlyOutput("dbs_plot", height = "350px"),
                  helpText("Deep Brain Stimulation (DBS) targets VC/VS or STN, disrupting the CSTC hyperactivity loop.")
                )
              )
          )
        )
      ),

      ## --------------------------------------------------------
      ## TAB 8: About
      ## --------------------------------------------------------
      tabItem(tabName = "tab_about",
        fluidRow(
          box(title = "About This OCD QSP Model", status = "info",
              solidHeader = TRUE, width = 12,
              column(6,
                h3("Model Description"),
                tags$p("This dashboard implements a Quantitative Systems Pharmacology (QSP)
                        model for Obsessive-Compulsive Disorder (OCD), integrating:"),
                tags$ul(
                  tags$li("2-compartment PK models for sertraline, clomipramine, and risperidone"),
                  tags$li("Blood-brain barrier transport to CNS compartment"),
                  tags$li("SERT occupancy kinetics (Emax/Hill model)"),
                  tags$li("Serotonin synaptic dynamics with 5-HT1A autoreceptor desensitization"),
                  tags$li("Dopamine system and D2R occupancy by antipsychotics"),
                  tags$li("CSTC circuit ODE model (OFC, Caudate, Thalamus, Direct/Indirect pathways)"),
                  tags$li("BDNF neuroplasticity dynamics"),
                  tags$li("ERP/CBT therapy effects on OFC normalization"),
                  tags$li("Y-BOCS clinical endpoint with realistic delay")
                ),
                h3("Key Clinical Calibration"),
                tags$ul(
                  tags$li("Baseline Y-BOCS ≈ 28 (moderate-severe OCD)"),
                  tags$li("SSRI response rate ≈ 40–60% (≥35% Y-BOCS reduction)"),
                  tags$li("Clomipramine ≈ slightly superior efficacy to SSRIs"),
                  tags$li("ERP + SSRI combination > monotherapy"),
                  tags$li("≥80% SERT occupancy required for OCD efficacy (vs 60% for depression)"),
                  tags$li("Delayed response: 6–12 weeks (5-HT1A desensitization)")
                )
              ),
              column(6,
                h3("Key References"),
                tags$ol(
                  tags$li("Soomro GM et al. (2008) Cochrane Rev — SSRIs for OCD"),
                  tags$li("Foa EB et al. (2005) JAMA — ERP vs Clomipramine RCT"),
                  tags$li("Zitterl W et al. (2008) Neuropsychopharmacology — SERT occupancy"),
                  tags$li("Goodman WK et al. (1989) JAMA — Y-BOCS development"),
                  tags$li("Saxena S & Rauch SL (2004) Annu Rev Neurosci — CSTC circuit"),
                  tags$li("Bloch MH et al. (2006) Mol Psychiatry — augmentation meta-analysis"),
                  tags$li("Pittenger C et al. (2011) Biol Psychiatry — glutamate in OCD"),
                  tags$li("Milad MR & Rauch SL (2012) Trends Cogn Sci — OCD and beyond")
                ),
                h3("Model Limitations"),
                tags$ul(
                  tags$li("Simplified 2-cmt PK; ignores active metabolites for SSRIs"),
                  tags$li("CSTC model is phenomenological, not first-principles neural"),
                  tags$li("Neuroinflammation not explicitly modeled"),
                  tags$li("ERP effect is parameterized as a smoothed accumulator"),
                  tags$li("Not validated against patient-level clinical trial data"),
                  tags$li("For educational and research use only")
                )
              )
          )
        )
      )
    ) # tabItems
  ) # dashboardBody
) # dashboardPage

## ============================================================
## SERVER
## ============================================================
server <- function(input, output, session) {

  ## -- Shared reactive simulation ----------------------------
  sim_data <- reactive({
    # Tab-5 driven main simulation
    run_sim(
      YBOCS0_in  = input$YBOCS0_c,
      SSRI_dose  = input$SSRI_dose_c,
      SSRI_on    = input$SSRI_on_c,
      CMI_dose   = input$CMI_dose_c,
      CMI_on     = input$CMI_on_c,
      RISP_dose  = input$RISP_dose_c,
      AUG_on     = input$AUG_on_c,
      ERP_on     = input$ERP_on_c,
      sim_weeks  = input$sim_wk_c,
      SERT_EC50_in = 1.20,
      CL_in      = 28
    )
  })

  sim_pk <- reactive({
    run_sim(
      SSRI_dose  = input$SSRI_dose,
      SSRI_on    = TRUE,
      CMI_on     = input$CMI_on,
      CMI_dose   = input$CMI_dose,
      sim_weeks  = input$sim_weeks_pk,
      SERT_EC50_in = input$SERT_EC50,
      CL_in      = input$CL_SSRI
    )
  })

  sim_pd <- reactive({
    run_sim(
      SSRI_dose  = input$SSRI_dose_pd,
      SSRI_on    = input$SSRI_on_pd,
      CMI_on     = input$CMI_on_pd,
      CMI_dose   = input$CMI_dose_pd,
      ERP_on     = input$ERP_on_pd,
      ERP_kmax_in= input$ERP_kmax_pd,
      sim_weeks  = input$sim_wk_pd
    )
  })

  ## -- TAB 1: Patient Summary --------------------------------
  output$patient_summary <- renderTable({
    data.frame(
      Parameter = c("Age", "Weight", "Sex", "Baseline Y-BOCS",
                    "OCD Duration (est.)", "CYP2D6"),
      Value = c(paste(input$age, "years"),
                paste(input$weight, "kg"),
                input$sex,
                paste(input$YBOCS0, "/40"),
                ">2 years (typical)",
                input$CYP2D6)
    )
  }, bordered = TRUE, striped = TRUE)

  output$severity_gauge <- renderPlot({
    score <- input$YBOCS0
    df <- data.frame(x = 1, y = score)
    ggplot(df, aes(x = x, y = y)) +
      geom_bar(stat = "identity", fill = ifelse(score <= 23, "#27AE60",
                                                 ifelse(score <= 31, "#E67E22", "#E74C3C")),
               width = 0.5) +
      geom_hline(yintercept = c(16, 24, 32, 40), linetype = "dashed",
                 color = c("#27AE60","#E67E22","#E74C3C","#C0392B")) +
      annotate("text", x = 1.4, y = c(20, 28, 36),
               label = c("Mild", "Moderate", "Severe"), size = 4) +
      scale_y_continuous(limits = c(0, 40)) +
      labs(x = "", y = "Y-BOCS Score", title = "") +
      theme_minimal() + theme(axis.text.x = element_blank())
  })

  ## -- TAB 2: PK Plots ---------------------------------------
  output$pk_plot <- renderPlotly({
    d <- sim_pk()
    p <- ggplot(d, aes(x = Week)) +
      geom_line(aes(y = Cp_ng, color = "Plasma Cp"), linewidth = 1) +
      geom_line(aes(y = Ccns_ng, color = "CNS Ccns"), linewidth = 1,
                linetype = "dashed") +
      scale_color_manual(values = c("Plasma Cp" = "#2E86C1", "CNS Ccns" = "#E74C3C")) +
      geom_hline(yintercept = c(30, 200), linetype = "dotted", color = "gray") +
      annotate("text", x = max(d$Week), y = 32, label = "30 ng/mL (lower therapeutic)",
               hjust = 1, size = 3, color = "gray50") +
      annotate("text", x = max(d$Week), y = 202, label = "200 ng/mL (upper)",
               hjust = 1, size = 3, color = "gray50") +
      labs(x = "Week", y = "Concentration (ng/mL)", color = "Compartment",
           title = "SSRI Plasma and CNS Concentrations") +
      theme_bw()
    ggplotly(p)
  })

  output$sert_plot <- renderPlotly({
    d <- sim_pk()
    p <- ggplot(d, aes(x = Week, y = SERT_pct)) +
      geom_line(color = "#E74C3C", linewidth = 1.2) +
      geom_hline(yintercept = 80, linetype = "dashed", color = "navy") +
      annotate("text", x = max(d$Week), y = 82, label = "≥80% required for OCD",
               hjust = 1, size = 3, color = "navy") +
      labs(x = "Week", y = "SERT Occupancy (%)", title = "SERT Occupancy") +
      ylim(0, 100) + theme_bw()
    ggplotly(p)
  })

  output$pk_table <- renderDT({
    d <- sim_pk() %>%
      filter(Week %in% c(1, 2, 4, 8, 12, 24)) %>%
      select(Week, Cp_ng, Ccns_ng, SERT_pct) %>%
      mutate_if(is.numeric, ~ round(., 1))
    names(d) <- c("Week", "Plasma Cp (ng/mL)", "CNS Ccns (ng/mL)", "SERT Occ (%)")
    datatable(d, options = list(dom = "t", pageLength = 10))
  })

  ## -- TAB 3: PD Plots ---------------------------------------
  output$ht5_plot <- renderPlotly({
    d <- sim_pd()
    p <- ggplot(d, aes(x = Week)) +
      geom_line(aes(y = HT5_norm, color = "Synaptic 5-HT"), linewidth = 1.2) +
      geom_line(aes(y = SERT_pct / 100, color = "SERT Occupancy"), linewidth = 1) +
      geom_hline(yintercept = 1.0, linetype = "dotted", color = "black") +
      scale_color_manual(values = c("Synaptic 5-HT" = "#E74C3C",
                                     "SERT Occupancy" = "#2E86C1")) +
      labs(x = "Week", y = "Normalized Value",
           color = "Variable", title = "Synaptic 5-HT Dynamics") +
      theme_bw()
    ggplotly(p)
  })

  output$des_plot <- renderPlotly({
    d <- sim_pd()
    p <- ggplot(d, aes(x = Week, y = DES_5HT1)) +
      geom_line(color = "#8E44AD", linewidth = 1.2) +
      labs(x = "Week", y = "5-HT1A Desensitization (0–1)",
           title = "5-HT1A Autoreceptor Desensitization") +
      ylim(0, 1) + theme_bw()
    ggplotly(p)
  })

  output$sert_5ht_scatter <- renderPlotly({
    d <- sim_pd()
    p <- ggplot(d, aes(x = SERT_pct, y = HT5_norm, color = Week)) +
      geom_path(linewidth = 0.8) +
      geom_point(size = 1) +
      scale_color_viridis_c() +
      labs(x = "SERT Occupancy (%)", y = "Synaptic 5-HT (normalized)",
           title = "SERT Occupancy vs Synaptic 5-HT", color = "Week") +
      theme_bw()
    ggplotly(p)
  })

  output$des_timeline <- renderPlotly({
    d <- sim_pd()
    p <- ggplot(d, aes(x = Week)) +
      geom_area(aes(y = DES_5HT1), fill = "#8E44AD", alpha = 0.3) +
      geom_line(aes(y = DES_5HT1, color = "5-HT1A Desens."), linewidth = 1) +
      geom_line(aes(y = HT5_norm - 1, color = "ΔSynaptic 5-HT"), linewidth = 1,
                linetype = "dashed") +
      scale_color_manual(values = c("5-HT1A Desens." = "#8E44AD",
                                     "ΔSynaptic 5-HT" = "#E74C3C")) +
      labs(x = "Week", y = "Value", color = "",
           title = "Desensitization Drives Delayed SSRI Response") +
      theme_bw()
    ggplotly(p)
  })

  ## -- TAB 4: CSTC -------------------------------------------
  output$cstc_plot <- renderPlotly({
    d <- sim_data()
    df_long <- d %>%
      select(Week, OFC_norm, CAUD_norm, THAL_norm, ANXIETY_norm = Anxiety_norm) %>%
      pivot_longer(-Week, names_to = "Region", values_to = "Activity") %>%
      mutate(Region = recode(Region,
        OFC_norm = "OFC (Orbitofrontal Cortex)",
        CAUD_norm = "Caudate Nucleus",
        THAL_norm = "MD Thalamus",
        ANXIETY_norm = "Anxiety State"
      ))
    p <- ggplot(df_long, aes(x = Week, y = Activity, color = Region)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 1.0, linetype = "dotted", color = "black") +
      annotate("text", x = max(d$Week), y = 1.02, label = "Normal baseline",
               hjust = 1, size = 3) +
      scale_color_manual(values = c(
        "OFC (Orbitofrontal Cortex)" = "#E74C3C",
        "Caudate Nucleus" = "#2E86C1",
        "MD Thalamus" = "#8E44AD",
        "Anxiety State" = "#E67E22"
      )) +
      labs(x = "Week", y = "Activity (normalized, 1.0 = normal)",
           color = "Brain Region/State", title = "CSTC Circuit Activity Normalization") +
      theme_bw()
    ggplotly(p)
  })

  output$dir_ind_plot <- renderPlotly({
    d <- sim_data()
    df <- d %>%
      select(Week, DIR_norm, IND_norm) %>%
      pivot_longer(-Week, names_to = "Pathway", values_to = "Activity") %>%
      mutate(Pathway = recode(Pathway,
        DIR_norm = "Direct (Go) Pathway",
        IND_norm = "Indirect (No-Go) Pathway"
      ))
    p <- ggplot(df, aes(x = Week, y = Activity, color = Pathway)) +
      geom_line(linewidth = 1.2) +
      geom_hline(yintercept = 1.0, linetype = "dotted") +
      scale_color_manual(values = c("Direct (Go) Pathway" = "#E74C3C",
                                     "Indirect (No-Go) Pathway" = "#1E8449")) +
      labs(x = "Week", y = "Pathway Activation",
           color = "", title = "Go vs No-Go Pathway Balance") + theme_bw()
    ggplotly(p)
  })

  output$bdnf_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = Week, y = BDNF_norm)) +
      geom_line(color = "#27AE60", linewidth = 1.2) +
      geom_hline(yintercept = 1.0, linetype = "dotted") +
      labs(x = "Week", y = "BDNF Level (normalized)",
           title = "BDNF Neuroplasticity Dynamics") + theme_bw()
    ggplotly(p)
  })

  ## -- TAB 5: Clinical ----------------------------------------
  output$ybocs_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = Week, y = YBOCS_score)) +
      geom_line(color = "#C0392B", linewidth = 1.5) +
      geom_hline(yintercept = 12, linetype = "dashed", color = "navy") +
      geom_hline(yintercept = input$YBOCS0_c * 0.65,
                 linetype = "dotted", color = "darkgreen") +
      annotate("text", x = max(d$Week), y = 13.5, label = "Remission (Y-BOCS≤12)",
               hjust = 1, size = 3, color = "navy") +
      annotate("text", x = max(d$Week), y = input$YBOCS0_c * 0.65 + 1.5,
               label = "Response (35% reduction)", hjust = 1, size = 3, color = "darkgreen") +
      scale_x_continuous(breaks = seq(0, input$sim_wk_c, 4)) +
      ylim(0, 40) +
      labs(x = "Week", y = "Y-BOCS Score",
           title = paste("Y-BOCS Trajectory (Baseline =", input$YBOCS0_c, ")")) +
      theme_bw()
    ggplotly(p)
  })

  final_vals <- reactive({
    d <- sim_data()
    last <- tail(d, 1)
    list(
      ybocs  = round(last$YBOCS_score, 1),
      pct    = round(last$YBOCS_pct_chg, 1),
      resp   = last$Responder,
      remiss = last$In_remission
    )
  })

  output$vbox_ybocs_final <- renderValueBox({
    v <- final_vals()
    valueBox(v$ybocs, "Final Y-BOCS", icon = icon("chart-bar"),
             color = ifelse(v$ybocs <= 12, "green", ifelse(v$ybocs <= 20, "yellow", "red")))
  })
  output$vbox_pct_chg <- renderValueBox({
    v <- final_vals()
    valueBox(paste0(v$pct, "%"), "Y-BOCS Reduction", icon = icon("arrow-down"),
             color = ifelse(v$pct >= 35, "green", ifelse(v$pct >= 25, "yellow", "red")))
  })
  output$vbox_responder <- renderValueBox({
    v <- final_vals()
    valueBox(ifelse(v$resp == 1, "YES", "NO"), "Responder (≥35%)",
             icon = icon("check"), color = ifelse(v$resp == 1, "green", "red"))
  })
  output$vbox_remission <- renderValueBox({
    v <- final_vals()
    valueBox(ifelse(v$remiss == 1, "YES", "NO"), "Remission (≤12)",
             icon = icon("star"), color = ifelse(v$remiss == 1, "green", "yellow"))
  })

  output$anxiety_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = Week, y = Anxiety_norm)) +
      geom_line(color = "#E67E22", linewidth = 1.2) +
      geom_hline(yintercept = 1.0, linetype = "dotted") +
      labs(x = "Week", y = "Anxiety State (normalized)", title = "Anxiety State Dynamics") +
      theme_bw()
    ggplotly(p)
  })

  output$response_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = Week)) +
      geom_line(aes(y = YBOCS_pct_chg, color = "% Y-BOCS Reduction"), linewidth = 1.2) +
      geom_hline(yintercept = 35, linetype = "dashed", color = "#27AE60") +
      annotate("text", x = 2, y = 37, label = "Response threshold (35%)", size = 3) +
      scale_color_manual(values = c("% Y-BOCS Reduction" = "#2E86C1")) +
      labs(x = "Week", y = "% Y-BOCS Reduction", color = "",
           title = "Treatment Response Trajectory") + theme_bw()
    ggplotly(p)
  })

  ## -- TAB 6: Scenario Comparison ----------------------------
  sc_data <- reactive({
    weeks <- input$sc_weeks
    ybocs0 <- input$sc_YBOCS0
    dose <- input$sc_ssri_dose

    scenarios <- list()
    if (input$sc_include_untreated) {
      scenarios[["1: Untreated"]] <- run_sim(YBOCS0_in=ybocs0, SSRI_on=FALSE, sim_weeks=weeks)
    }
    if (input$sc_include_ssri) {
      scenarios[[paste0("2: Sertraline ", dose, "mg")]] <-
        run_sim(YBOCS0_in=ybocs0, SSRI_dose=dose, SSRI_on=TRUE, sim_weeks=weeks)
    }
    if (input$sc_include_cmi) {
      scenarios[["3: Clomipramine 250mg"]] <-
        run_sim(YBOCS0_in=ybocs0, CMI_on=TRUE, CMI_dose=250, SSRI_on=FALSE, sim_weeks=weeks)
    }
    if (input$sc_include_aug) {
      scenarios[["4: SSRI+Risperidone"]] <-
        run_sim(YBOCS0_in=ybocs0, SSRI_dose=dose, SSRI_on=TRUE,
                AUG_on=TRUE, RISP_dose=1.5, sim_weeks=weeks)
    }
    if (input$sc_include_combo) {
      scenarios[["5: SSRI+ERP"]] <-
        run_sim(YBOCS0_in=ybocs0, SSRI_dose=dose, SSRI_on=TRUE,
                ERP_on=TRUE, ERP_kmax_in=0.45, sim_weeks=weeks)
    }
    if (input$sc_include_erp) {
      scenarios[["6: ERP Alone"]] <-
        run_sim(YBOCS0_in=ybocs0, SSRI_on=FALSE, ERP_on=TRUE, ERP_kmax_in=0.35, sim_weeks=weeks)
    }

    bind_rows(lapply(names(scenarios), function(nm) {
      scenarios[[nm]] %>% mutate(Scenario = nm)
    }))
  })

  output$scenario_ybocs <- renderPlotly({
    d <- sc_data()
    p <- ggplot(d, aes(x = Week, y = YBOCS_score, color = Scenario)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 12, linetype = "dashed", color = "navy") +
      labs(x = "Week", y = "Y-BOCS Score", title = "Y-BOCS — Scenario Comparison") +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$scenario_sert <- renderPlotly({
    d <- sc_data()
    p <- ggplot(d, aes(x = Week, y = SERT_pct, color = Scenario)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 80, linetype = "dashed", color = "navy") +
      labs(x = "Week", y = "SERT Occupancy (%)", title = "SERT Occupancy — Scenario Comparison") +
      ylim(0, 100) + theme_bw() + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$scenario_table <- renderDT({
    d <- sc_data() %>%
      filter(round(Week) %in% c(12, 24, 52)) %>%
      group_by(Scenario, Week = round(Week)) %>%
      summarise(
        `Y-BOCS`   = round(mean(YBOCS_score), 1),
        `%Reduction`= round(mean(YBOCS_pct_chg), 1),
        `SERT Occ%` = round(mean(SERT_pct), 1),
        `OFC Act`   = round(mean(OFC_norm), 3),
        Response    = ifelse(mean(YBOCS_pct_chg) >= 35, "YES", "no"),
        Remission   = ifelse(mean(YBOCS_score) <= 12, "YES", "no"),
        .groups = "drop"
      )
    datatable(d, options = list(dom = "t", pageLength = 30))
  })

  ## -- TAB 7: Biomarkers --------------------------------------
  output$biomarker_circuit <- renderPlotly({
    d <- sim_data()
    df <- d %>%
      select(Week, OFC_norm, CAUD_norm, THAL_norm) %>%
      pivot_longer(-Week) %>%
      mutate(name = recode(name, OFC_norm="OFC (FDG-PET)",
                           CAUD_norm="Caudate (FDG-PET)",
                           THAL_norm="MD Thalamus (fMRI)"))
    p <- ggplot(df, aes(x = Week, y = value, color = name)) +
      geom_line(linewidth = 1.2) +
      geom_hline(yintercept = 1.0, linetype = "dotted") +
      labs(x = "Week", y = "Activity (normalized)",
           color = "Biomarker (Method)", title = "Neuroimaging Biomarkers") +
      theme_bw()
    ggplotly(p)
  })

  output$biomarker_neuro <- renderPlotly({
    d <- sim_data()
    df <- d %>%
      select(Week, SERT_pct, HT5_norm, BDNF_norm) %>%
      mutate(SERT_pct = SERT_pct / 100) %>%
      pivot_longer(-Week) %>%
      mutate(name = recode(name, SERT_pct="SERT Occ (PET)",
                           HT5_norm="5-HT (CSF 5-HIAA proxy)",
                           BDNF_norm="BDNF (Serum)"))
    p <- ggplot(df, aes(x = Week, y = value, color = name)) +
      geom_line(linewidth = 1.2) +
      geom_hline(yintercept = 1.0, linetype = "dotted") +
      labs(x = "Week", y = "Biomarker Value (normalized or fraction)",
           color = "Biomarker", title = "Neurochemical Biomarkers") +
      theme_bw()
    ggplotly(p)
  })

  output$pgx_sert <- renderPlotly({
    # Simulate different CYP2D6 genotypes → different CL → different SERT occ
    cyp_cl <- c(PM = 12, IM = 20, NM = 28, UM = 50)
    pgx_df <- bind_rows(lapply(names(cyp_cl), function(geno) {
      run_sim(SSRI_dose = 200, SSRI_on = TRUE, sim_weeks = 24,
              CL_in = cyp_cl[[geno]]) %>%
        mutate(Genotype = geno)
    }))
    p <- ggplot(pgx_df, aes(x = Week, y = SERT_pct, color = Genotype)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 80, linetype = "dashed", color = "navy") +
      labs(x = "Week", y = "SERT Occupancy (%)",
           title = "CYP2D6 Genotype → SERT Occupancy\n(Sertraline 200 mg)") +
      ylim(0, 100) + theme_bw()
    ggplotly(p)
  })

  output$pgx_ybocs <- renderPlotly({
    cyp_cl <- c(PM = 12, IM = 20, NM = 28, UM = 50)
    pgx_df <- bind_rows(lapply(names(cyp_cl), function(geno) {
      run_sim(SSRI_dose = 200, SSRI_on = TRUE, sim_weeks = 24,
              CL_in = cyp_cl[[geno]]) %>%
        mutate(Genotype = geno)
    }))
    p <- ggplot(pgx_df, aes(x = Week, y = YBOCS_score, color = Genotype)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 12, linetype = "dashed", color = "navy") +
      labs(x = "Week", y = "Y-BOCS Score",
           title = "CYP2D6 Genotype → Y-BOCS Response\n(Sertraline 200 mg)") +
      ylim(0, 40) + theme_bw()
    ggplotly(p)
  })

  output$dbs_plot <- renderPlotly({
    dbs <- input$dbs_effect
    d <- run_sim(SSRI_on = FALSE, ERP_on = FALSE) %>%
      mutate(OFC_dbs  = OFC_norm  * (1 - 0.30 * dbs),
             CAUD_dbs = CAUD_norm * (1 - 0.25 * dbs),
             YBOCS_dbs = pmax(0, YBOCS_score * (1 - 0.40 * dbs)))
    df <- d %>%
      select(Week, OFC_dbs, CAUD_dbs, YBOCS_dbs) %>%
      pivot_longer(-Week, names_to = "Metric", values_to = "Value") %>%
      mutate(Metric = recode(Metric, OFC_dbs = "OFC Activity",
                             CAUD_dbs = "Caudate Activity",
                             YBOCS_dbs = "Y-BOCS (scaled)"))
    p <- ggplot(df, aes(x = Week, y = Value, color = Metric)) +
      geom_line(linewidth = 1) +
      labs(x = "Week", y = "Value",
           title = paste0("DBS Effect = ", dbs, " (VC/VS or STN target)"),
           color = "") +
      theme_bw()
    ggplotly(p)
  })
}

## ============================================================
## LAUNCH
## ============================================================
shinyApp(ui = ui, server = server)
