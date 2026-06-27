##############################################################################
# Hypereosinophilic Syndrome (HES) QSP — Interactive Shiny Dashboard
# 6 tabs: Patient Profile · PK · Eosinophil Kinetics ·
#         Organ Damage · Treatment Comparison · Biomarker Correlations
##############################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(patchwork)

# ─── Embed model code ────────────────────────────────────────────────────────
hes_model_code <- '
$PARAM
AEC0=3000 IL5_0=5.0 kprod_IL5=0.05 kdeg_IL5=0.01
Emax_IL5=3.0 EC50_IL5=2.0 Emax_surv=0.7 EC50_surv=1.5
kprol_EoP=0.050 kmat_EoP=0.020 kmat_EoI=0.015 krel_BM=0.030 kdeath_EoP=0.008
kel_Eo=0.038 Vd_Eo=5.0 k_tissue=0.010
k_fibrosis=0.0002 krev_fibro=0.0001 FIBROSIS_0=0.05
k_pulm=0.0003 krev_pulm=0.0002 PULM_0=0.10
F_MEPO=0.80 ka_MEPO=0.0035 CL_MEPO=0.0072 Vc_MEPO=3.5
Q_MEPO=0.012 Vp_MEPO=3.5 kon_MEPO=10.0 koff_MEPO=0.001 kdeg_TMDD=0.05
F_BENRA=0.59 ka_BENRA=0.006 CL_BENRA=0.004 Vc_BENRA=3.0
Q_BENRA=0.015 Vp_BENRA=3.8 Emax_ADCC=0.95 EC50_ADCC=0.01 ADCC_hill=1.5
F_IMAT=0.98 ka_IMAT=0.30 CL_IMAT=12.0 Vc_IMAT=110.0 IMAT_IC50=0.10 IMAT_hill=1.0
F_PRED=0.82 ka_PRED=0.80 CL_PRED=3.5 Vc_PRED=38.0
Emax_PRED=0.80 EC50_PRED=0.05 Emax_APO=0.60 EC50_APO=0.08
USE_MEPO=0 USE_BENRA=0 USE_IMAT=0 USE_PRED=0
CLONAL_HES=0 CLONAL_FOLD=5.0

$CMT MEPO_DEPOT MEPO_C1 MEPO_C2 TMDD
BENRA_DEPOT BENRA_C1 BENRA_C2
IMAT_GUT IMAT_C PRED_GUT PRED_C
EoP EoI EoM_BM EO_BLOOD IL5
FIBROSIS PULM_SCORE

$INIT MEPO_DEPOT=0 MEPO_C1=0 MEPO_C2=0 TMDD=0
BENRA_DEPOT=0 BENRA_C1=0 BENRA_C2=0
IMAT_GUT=0 IMAT_C=0 PRED_GUT=0 PRED_C=0
EoP=100 EoI=80 EoM_BM=150 EO_BLOOD=3000 IL5=5.0
FIBROSIS=0.05 PULM_SCORE=0.10

$ODE
double IL5_free=IL5>0?IL5:0;
double MEPO_Cp=MEPO_C1/Vc_MEPO;
double PRED_Cp=PRED_C;
double pred_il5_inhib=0.0;
if(USE_PRED>0.5) pred_il5_inhib=Emax_PRED*PRED_Cp/(EC50_PRED+PRED_Cp);
double IMAT_Cp=IMAT_C;
double imat_inhib=0.0;
if(USE_IMAT>0.5&&CLONAL_HES>0.5) imat_inhib=IMAT_Cp/(IMAT_IC50+IMAT_Cp);
double clonal_factor=1.0;
if(CLONAL_HES>0.5) clonal_factor=CLONAL_FOLD*(1.0-imat_inhib);
double PROD_IL5=kprod_IL5*(1.0-pred_il5_inhib);
double DEG_IL5=kdeg_IL5*IL5;
double TMDD_on=0.0,TMDD_off=0.0;
if(USE_MEPO>0.5){TMDD_on=kon_MEPO*(MEPO_C1/Vc_MEPO)*IL5_free;TMDD_off=koff_MEPO*TMDD;}
dxdt_IL5=PROD_IL5-DEG_IL5-TMDD_on+TMDD_off;
double MEPO_ABS=(USE_MEPO>0.5)?(ka_MEPO*F_MEPO*MEPO_DEPOT):0.0;
dxdt_MEPO_DEPOT=-ka_MEPO*MEPO_DEPOT;
dxdt_MEPO_C1=MEPO_ABS-(CL_MEPO/Vc_MEPO)*MEPO_C1-Q_MEPO*(MEPO_C1/Vc_MEPO-MEPO_C2/Vp_MEPO)-TMDD_on*Vc_MEPO+TMDD_off*Vc_MEPO;
dxdt_MEPO_C2=Q_MEPO*(MEPO_C1/Vc_MEPO-MEPO_C2/Vp_MEPO);
dxdt_TMDD=TMDD_on-TMDD_off-kdeg_TMDD*TMDD;
double BENRA_ABS=(USE_BENRA>0.5)?(ka_BENRA*F_BENRA*BENRA_DEPOT):0.0;
dxdt_BENRA_DEPOT=-ka_BENRA*BENRA_DEPOT;
dxdt_BENRA_C1=BENRA_ABS-(CL_BENRA/Vc_BENRA)*BENRA_C1-Q_BENRA*(BENRA_C1/Vc_BENRA-BENRA_C2/Vp_BENRA);
dxdt_BENRA_C2=Q_BENRA*(BENRA_C1/Vc_BENRA-BENRA_C2/Vp_BENRA);
dxdt_IMAT_GUT=-(USE_IMAT>0.5?ka_IMAT:0.0)*IMAT_GUT;
dxdt_IMAT_C=(USE_IMAT>0.5?ka_IMAT*F_IMAT*IMAT_GUT:0.0)-(CL_IMAT/Vc_IMAT)*IMAT_C;
dxdt_PRED_GUT=-(USE_PRED>0.5?ka_PRED:0.0)*PRED_GUT;
dxdt_PRED_C=(USE_PRED>0.5?ka_PRED*F_PRED*PRED_GUT:0.0)-(CL_PRED/Vc_PRED)*PRED_C;
double IL5_Eprol=Emax_IL5*IL5_free/(EC50_IL5+IL5_free);
double pred_apo=0.0;
if(USE_PRED>0.5) pred_apo=Emax_APO*PRED_Cp/(EC50_APO+PRED_Cp);
double BENRA_Cp=BENRA_C1/Vc_BENRA;
double adcc_eff=0.0;
if(USE_BENRA>0.5){double Bn=pow(BENRA_Cp,ADCC_hill);double En=pow(EC50_ADCC,ADCC_hill);adcc_eff=Emax_ADCC*Bn/(En+Bn);}
double PROL_RATE=kprol_EoP*(1.0+IL5_Eprol)*clonal_factor*EoP;
double DEATH_EoP=kdeath_EoP*(1.0+pred_apo)*EoP;
dxdt_EoP=PROL_RATE-kmat_EoP*EoP-DEATH_EoP;
dxdt_EoI=kmat_EoP*EoP-kmat_EoI*EoI;
dxdt_EoM_BM=kmat_EoI*EoI-krel_BM*EoM_BM;
double IL5_surv=Emax_surv*IL5_free/(EC50_surv+IL5_free);
double kel_eff=kel_Eo*(1.0-IL5_surv);
double adcc_elim=adcc_eff*EO_BLOOD;
dxdt_EO_BLOOD=krel_BM*EoM_BM-kel_eff*EO_BLOOD-k_tissue*EO_BLOOD-adcc_elim;
double EO_norm=EO_BLOOD/500.0;
double FIBRO_IN=k_fibrosis*(EO_norm>1.0?EO_norm-1.0:0.0)*(1.0-FIBROSIS);
dxdt_FIBROSIS=FIBRO_IN-krev_fibro*FIBROSIS;
double PULM_IN=k_pulm*(EO_norm>1.0?EO_norm-1.0:0.0)*(1.0-PULM_SCORE);
dxdt_PULM_SCORE=PULM_IN-krev_pulm*PULM_SCORE;

$TABLE
double AEC_obs=EO_BLOOD;
double MEPO_Cobs=MEPO_C1/Vc_MEPO;
double BENRA_Cobs=BENRA_C1/Vc_BENRA;
double IMAT_Cobs=IMAT_C;
double PRED_Cobs=PRED_C;
double PCT_CHG=(AEC_obs-AEC0)/AEC0*100.0;
double RESP_300=(AEC_obs<300)?1.0:0.0;
double RESP_1500=(AEC_obs<1500)?1.0:0.0;

$CAPTURE AEC_obs MEPO_Cobs BENRA_Cobs IMAT_Cobs PRED_Cobs
PCT_CHG RESP_300 RESP_1500 IL5 FIBROSIS PULM_SCORE EoP EoI EoM_BM
'
mod <- mcode("HES_Shiny", hes_model_code, quiet = TRUE)

# ─── Helper: run simulation ──────────────────────────────────────────────────
run_sim <- function(params, events, end_wk = 52) {
  end_h <- end_wk * 7 * 24
  do.call(param, c(list(mod), params)) %>%
    mrgsim(events = events, end = end_h, delta = 24) %>%
    as.data.frame() %>%
    mutate(time_wk = time / (7 * 24))
}

build_events <- function(drug, dose_mepo, dose_benra, dose_imat, dose_pred) {
  evlist <- NULL
  if (drug %in% c("Mepolizumab", "Both Anti-IL5")) {
    e <- ev(cmt = "MEPO_DEPOT", amt = dose_mepo * 1000,
            ii = 4 * 7 * 24, addl = 12, time = 0)
    evlist <- if (is.null(evlist)) e else c(evlist, e)
  }
  if (drug %in% c("Benralizumab", "Both Anti-IL5")) {
    loading <- ev(cmt = "BENRA_DEPOT", amt = dose_benra * 1000,
                  ii = 4 * 7 * 24, addl = 2, time = 0)
    maint   <- ev(cmt = "BENRA_DEPOT", amt = dose_benra * 1000,
                  ii = 8 * 7 * 24, addl = 5, time = 3 * 4 * 7 * 24)
    e <- c(loading, maint)
    evlist <- if (is.null(evlist)) e else c(evlist, e)
  }
  if (drug == "Imatinib (Clonal)") {
    e <- ev(cmt = "IMAT_GUT", amt = dose_imat,
            ii = 24, addl = 52 * 7 - 1, time = 0)
    evlist <- if (is.null(evlist)) e else c(evlist, e)
  }
  if (drug == "Prednisolone") {
    e <- ev(cmt = "PRED_GUT", amt = dose_pred,
            ii = 24, addl = 52 * 7 - 1, time = 0)
    evlist <- if (is.null(evlist)) e else c(evlist, e)
  }
  if (is.null(evlist)) ev(cmt = "EO_BLOOD", amt = 0, time = 0)
  else evlist
}

make_params <- function(drug, aec0, clonal) {
  p <- list(
    AEC0 = aec0, EO_BLOOD = aec0, EoP = aec0 / 30,
    USE_MEPO  = as.integer(drug %in% c("Mepolizumab", "Both Anti-IL5")),
    USE_BENRA = as.integer(drug %in% c("Benralizumab", "Both Anti-IL5")),
    USE_IMAT  = as.integer(drug == "Imatinib (Clonal)"),
    USE_PRED  = as.integer(drug == "Prednisolone"),
    CLONAL_HES = as.integer(clonal)
  )
  p
}

# ─────────────────────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(
    title = span(icon("heartbeat"), "HES QSP Dashboard"),
    titleWidth = 320
  ),

  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("Patient Profile",        tabName = "tab_patient",  icon = icon("user-md")),
      menuItem("Pharmacokinetics",       tabName = "tab_pk",       icon = icon("flask")),
      menuItem("Eosinophil Kinetics",    tabName = "tab_eo",       icon = icon("tint")),
      menuItem("Organ Damage",           tabName = "tab_organ",    icon = icon("heart")),
      menuItem("Treatment Comparison",   tabName = "tab_compare",  icon = icon("balance-scale")),
      menuItem("Biomarker Correlations", tabName = "tab_biom",     icon = icon("chart-line"))
    ),
    hr(),
    h5("── Patient Parameters ──", style = "color:#ccc; padding-left:10px"),
    sliderInput("aec0",  "Baseline AEC (cells/µL)", 500, 15000, 3000, step = 500),
    checkboxInput("clonal", "Clonal HES (FIP1L1-PDGFRA+)", FALSE),
    sliderInput("clonal_fold", "Clonal EoP fold ×", 1, 10, 5, step = 0.5),
    hr(),
    h5("── Drug Selection ──", style = "color:#ccc; padding-left:10px"),
    selectInput("drug", "Treatment",
                choices = c("Untreated", "Prednisolone", "Mepolizumab",
                            "Benralizumab", "Imatinib (Clonal)", "Both Anti-IL5"),
                selected = "Mepolizumab"),
    conditionalPanel("input.drug == 'Mepolizumab' || input.drug == 'Both Anti-IL5'",
      sliderInput("dose_mepo",  "Mepolizumab dose (mg)", 100, 600, 300, step = 100)),
    conditionalPanel("input.drug == 'Benralizumab' || input.drug == 'Both Anti-IL5'",
      sliderInput("dose_benra", "Benralizumab dose (mg)", 10, 60, 30, step = 10)),
    conditionalPanel("input.drug == 'Imatinib (Clonal)'",
      sliderInput("dose_imat",  "Imatinib dose (mg/day)", 50, 400, 100, step = 50)),
    conditionalPanel("input.drug == 'Prednisolone'",
      sliderInput("dose_pred",  "Prednisolone dose (mg/day)", 10, 100, 70, step = 5)),
    sliderInput("sim_wk", "Simulation duration (weeks)", 12, 104, 52, step = 4),
    actionButton("run_btn", "Run Simulation", icon = icon("play"),
                 class = "btn-success", style = "width:90%;margin:10px;")
  ),

  dashboardBody(
    tags$head(tags$style(HTML(
      ".content-wrapper { background-color: #f4f6f9; }
       .box-header h3 { font-size: 16px; font-weight: bold; }
       .nav-tabs-custom .nav-tabs li a { font-size: 13px; }"
    ))),
    tabItems(

      # ── TAB 1: Patient Profile ─────────────────────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(width = 12, title = "Hypereosinophilic Syndrome — Disease Overview",
              status = "danger", solidHeader = TRUE,
              fluidRow(
                valueBoxOutput("vbox_aec",    width = 3),
                valueBoxOutput("vbox_hes_type", width = 3),
                valueBoxOutput("vbox_cardiac", width = 3),
                valueBoxOutput("vbox_pulm",   width = 3)
              )
          )
        ),
        fluidRow(
          box(width = 6, title = "Disease Summary", status = "info", solidHeader = TRUE,
              HTML("
                <b>Hypereosinophilic Syndrome (HES)</b> is defined by sustained peripheral blood
                eosinophilia ≥1500 cells/µL for ≥1 month with evidence of eosinophil-mediated
                organ damage, in the absence of secondary causes.
                <br><br>
                <b>Key Subtypes:</b>
                <ul>
                  <li><b>Lymphocytic HES (L-HES)</b>: Aberrant Th2 T-cell clones → IL-5 overproduction</li>
                  <li><b>Myeloid/Clonal HES (M-HES)</b>: FIP1L1-PDGFRA fusion → constitutive TK activation</li>
                  <li><b>Reactive HES</b>: Secondary to parasites, atopy, malignancy</li>
                  <li><b>Idiopathic HES (IHES)</b>: No identifiable cause</li>
                </ul>
                <b>Organs affected:</b> Heart (Löffler endocarditis → EMF), Lungs, GI tract,
                Skin, Nervous system
              ")
          ),
          box(width = 6, title = "Selected Patient Profile", status = "success", solidHeader = TRUE,
              DTOutput("patient_profile_table")
          )
        ),
        fluidRow(
          box(width = 12, title = "Mechanism of Action Summary", status = "warning", solidHeader = TRUE,
              HTML("
                <div style='display:flex;gap:20px;flex-wrap:wrap'>
                  <div style='flex:1;min-width:200px;background:#e8f5e9;padding:12px;border-radius:8px'>
                    <b style='color:#2e7d32'>Mepolizumab</b><br>
                    Anti-IL-5 IgG1 mAb → prevents IL-5 binding to IL-5Rα → ↓ EoP proliferation,
                    ↓ BM release, ↓ eosinophil survival
                  </div>
                  <div style='flex:1;min-width:200px;background:#e3f2fd;padding:12px;border-radius:8px'>
                    <b style='color:#1565c0'>Benralizumab</b><br>
                    Anti-IL-5Rα afucosylated IgG1 → blocks IL-5Rα + ADCC-mediated NK-cell depletion
                    → near-complete AEC depletion within weeks
                  </div>
                  <div style='flex:1;min-width:200px;background:#fff3e0;padding:12px;border-radius:8px'>
                    <b style='color:#e65100'>Imatinib</b><br>
                    BCR-ABL/PDGFRA/c-Kit TKI → FIP1L1-PDGFRA inhibition → ↓ clonal EoP
                    expansion → dramatic AEC reduction in M-HES
                  </div>
                  <div style='flex:1;min-width:200px;background:#fce4ec;padding:12px;border-radius:8px'>
                    <b style='color:#880e4f'>Prednisolone</b><br>
                    GR agonist → NF-κB inhibition → ↓ IL-5/IL-4/IL-13, promotes Eo apoptosis;
                    effective but long-term AEs limit use
                  </div>
                </div>
              ")
          )
        )
      ),

      # ── TAB 2: PK ─────────────────────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(width = 12, title = "Drug Pharmacokinetics", status = "primary",
              solidHeader = TRUE, plotlyOutput("pk_plot", height = "500px"))
        ),
        fluidRow(
          box(width = 6, title = "PK Parameters Summary", status = "info", solidHeader = TRUE,
              DTOutput("pk_table")),
          box(width = 6, title = "Drug Exposure Metrics (AUC, Cmax, Tmax)", status = "success",
              solidHeader = TRUE, DTOutput("pk_metrics"))
        )
      ),

      # ── TAB 3: Eosinophil Kinetics ────────────────────────────────────────
      tabItem(tabName = "tab_eo",
        fluidRow(
          box(width = 8, title = "Absolute Eosinophil Count (AEC) Over Time",
              status = "danger", solidHeader = TRUE,
              plotlyOutput("aec_plot", height = "380px")),
          box(width = 4, title = "AEC Response Metrics at Week 24",
              status = "warning", solidHeader = TRUE,
              valueBoxOutput("resp_300_box",  width = 12),
              valueBoxOutput("resp_1500_box", width = 12),
              hr(),
              uiOutput("aec_interpretation")
          )
        ),
        fluidRow(
          box(width = 6, title = "Bone Marrow Eosinophilopoiesis (BM Compartments)",
              status = "primary", solidHeader = TRUE,
              plotlyOutput("bm_plot", height = "320px")),
          box(width = 6, title = "IL-5 Dynamics and % Change AEC",
              status = "success", solidHeader = TRUE,
              plotlyOutput("il5_pct_plot", height = "320px"))
        )
      ),

      # ── TAB 4: Organ Damage ───────────────────────────────────────────────
      tabItem(tabName = "tab_organ",
        fluidRow(
          box(width = 6, title = "Cardiac Fibrosis Progression (Löffler → EMF → RCM)",
              status = "danger", solidHeader = TRUE,
              plotlyOutput("cardiac_plot", height = "360px")),
          box(width = 6, title = "Pulmonary Infiltration Score Over Time",
              status = "info", solidHeader = TRUE,
              plotlyOutput("pulm_plot", height = "360px"))
        ),
        fluidRow(
          box(width = 12, title = "Organ Damage Interpretation", status = "warning",
              solidHeader = TRUE,
              HTML("
                <b>Cardiac score interpretation:</b>
                <ul>
                  <li>0.0–0.1: Minimal endocardial involvement (Stage I Löffler)</li>
                  <li>0.1–0.3: Moderate fibrosis, early tricuspid/mitral regurgitation (Stage II)</li>
                  <li>0.3–0.6: Significant restrictive cardiomyopathy, wall thickening</li>
                  <li>0.6–1.0: Severe EMF, biventricular compromise, high mortality risk</li>
                </ul>
                <b>Pulmonary score interpretation:</b>
                <ul>
                  <li>0.0–0.1: No significant pulmonary infiltrates</li>
                  <li>0.1–0.3: Mild infiltrates, minimal DLCO reduction</li>
                  <li>0.3–0.6: Moderate eosinophilic pneumonia, pleural effusion possible</li>
                  <li>0.6–1.0: Severe respiratory compromise, pulmonary hypertension risk</li>
                </ul>
              ")
          )
        )
      ),

      # ── TAB 5: Treatment Comparison ───────────────────────────────────────
      tabItem(tabName = "tab_compare",
        fluidRow(
          box(width = 12, title = "5 Treatment Scenarios — Side-by-Side Comparison",
              status = "primary", solidHeader = TRUE,
              plotlyOutput("compare_plot", height = "500px"))
        ),
        fluidRow(
          box(width = 12, title = "Week 24 Comparative Summary Table",
              status = "success", solidHeader = TRUE,
              DTOutput("compare_table"))
        )
      ),

      # ── TAB 6: Biomarker Correlations ────────────────────────────────────
      tabItem(tabName = "tab_biom",
        fluidRow(
          box(width = 6, title = "AEC vs Cardiac Fibrosis Score",
              status = "danger", solidHeader = TRUE,
              plotlyOutput("biom_cardiac_plot", height = "380px")),
          box(width = 6, title = "IL-5 vs AEC Relationship",
              status = "info", solidHeader = TRUE,
              plotlyOutput("biom_il5_plot", height = "380px"))
        ),
        fluidRow(
          box(width = 6, title = "AEC vs Pulmonary Score",
              status = "primary", solidHeader = TRUE,
              plotlyOutput("biom_pulm_plot", height = "320px")),
          box(width = 6, title = "Serum Biomarker Trajectories",
              status = "success", solidHeader = TRUE,
              plotlyOutput("biom_multi_plot", height = "320px"))
        )
      )
    ) # end tabItems
  ) # end dashboardBody
)

# ─────────────────────────────────────────────────────────────────────────────
# SERVER
# ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Reactive simulation ───────────────────────────────────────────────────
  sim_result <- eventReactive(input$run_btn, {
    params <- make_params(input$drug, input$aec0, input$clonal)
    if (input$clonal) params$CLONAL_FOLD <- input$clonal_fold
    events <- build_events(input$drug, input$dose_mepo, input$dose_benra,
                           input$dose_imat, input$dose_pred)
    run_sim(params, events, end_wk = input$sim_wk)
  }, ignoreNULL = FALSE)

  # All-scenario comparison (fixed)
  compare_result <- reactive({
    end_h <- 52 * 7 * 24
    scenarios <- list(
      list(drug = "Untreated",     aec = 3000, clonal = FALSE, label = "Untreated"),
      list(drug = "Prednisolone",  aec = 3000, clonal = FALSE, label = "Prednisolone 70mg"),
      list(drug = "Mepolizumab",   aec = 3000, clonal = FALSE, label = "Mepolizumab 300mg q4w"),
      list(drug = "Benralizumab",  aec = 3000, clonal = FALSE, label = "Benralizumab 30mg q4w→q8w"),
      list(drug = "Imatinib (Clonal)", aec = 5000, clonal = TRUE, label = "Imatinib 100mg/day (Clonal)")
    )
    bind_rows(lapply(scenarios, function(sc) {
      p <- make_params(sc$drug, sc$aec, sc$clonal)
      e <- build_events(sc$drug, 300, 30, 100, 70)
      run_sim(p, e, end_wk = 52) %>% mutate(scenario = sc$label)
    }))
  })

  # ── Tab 1: Patient Profile value boxes ───────────────────────────────────
  output$vbox_aec <- renderValueBox({
    valueBox(
      format(input$aec0, big.mark = ","), "Baseline AEC (cells/µL)",
      icon = icon("tint"), color = if (input$aec0 > 5000) "red" else if (input$aec0 > 1500) "orange" else "green"
    )
  })
  output$vbox_hes_type <- renderValueBox({
    typ <- if (input$clonal) "Clonal (M-HES)" else "Reactive / Idiopathic"
    valueBox(typ, "HES Subtype", icon = icon("dna"), color = if (input$clonal) "purple" else "blue")
  })
  output$vbox_cardiac <- renderValueBox({
    d <- sim_result()
    last <- d %>% slice_tail(n = 1)
    valueBox(round(last$FIBROSIS, 3), "Cardiac Fibrosis Score", icon = icon("heart"), color = "red")
  })
  output$vbox_pulm <- renderValueBox({
    d <- sim_result()
    last <- d %>% slice_tail(n = 1)
    valueBox(round(last$PULM_SCORE, 3), "Pulmonary Score", icon = icon("lungs"), color = "navy")
  })
  output$patient_profile_table <- renderDT({
    df <- data.frame(
      Parameter  = c("Baseline AEC", "HES Subtype", "Treatment", "CLONAL_FOLD",
                     "Sim Duration (wk)"),
      Value      = c(paste(input$aec0, "cells/µL"),
                     if (input$clonal) "M-HES (Clonal)" else "L-HES / Idiopathic",
                     input$drug,
                     ifelse(input$clonal, input$clonal_fold, "N/A"),
                     input$sim_wk)
    )
    datatable(df, options = list(dom = "t", paging = FALSE), rownames = FALSE)
  })

  # ── Tab 2: PK ─────────────────────────────────────────────────────────────
  output$pk_plot <- renderPlotly({
    d <- sim_result()
    plot_ly() %>%
      add_lines(data = d, x = ~time_wk, y = ~MEPO_Cobs,  name = "Mepolizumab (µg/mL)", line = list(color = "#27AE60")) %>%
      add_lines(data = d, x = ~time_wk, y = ~BENRA_Cobs, name = "Benralizumab (µg/mL)", line = list(color = "#3498DB")) %>%
      add_lines(data = d, x = ~time_wk, y = ~IMAT_Cobs,  name = "Imatinib (µg/mL)",     line = list(color = "#E67E22")) %>%
      add_lines(data = d, x = ~time_wk, y = ~PRED_Cobs,  name = "Prednisolone (µg/mL)", line = list(color = "#D4AC0D")) %>%
      layout(xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Concentration (µg/mL)"),
             legend = list(orientation = "h"))
  })

  output$pk_table <- renderDT({
    df <- data.frame(
      Drug           = c("Mepolizumab", "Benralizumab", "Imatinib", "Prednisolone"),
      `Route`        = c("SC", "SC", "PO", "PO"),
      `t½ (days)`    = c(22, 15, 18, 0.83),
      `F (%)`        = c(80, 59, 98, 82),
      `Vc (L)`       = c(3.5, 3.0, 110, 38),
      `CL (L/h)`     = c(0.0072, 0.004, 12, 3.5),
      check.names = FALSE
    )
    datatable(df, options = list(dom = "t", paging = FALSE), rownames = FALSE)
  })
  output$pk_metrics <- renderDT({
    d <- sim_result()
    calc_metrics <- function(col, drug) {
      vals <- d[[col]]
      cmax <- max(vals, na.rm = TRUE)
      tmax <- d$time_wk[which.max(vals)]
      auc  <- sum(diff(d$time) * (vals[-1] + vals[-length(vals)]) / 2, na.rm = TRUE)
      data.frame(Drug = drug, Cmax_ugmL = round(cmax, 3),
                 Tmax_wk = round(tmax, 2), AUC_ugmL_h = round(auc, 1))
    }
    df <- bind_rows(
      calc_metrics("MEPO_Cobs",  "Mepolizumab"),
      calc_metrics("BENRA_Cobs", "Benralizumab"),
      calc_metrics("IMAT_Cobs",  "Imatinib"),
      calc_metrics("PRED_Cobs",  "Prednisolone")
    )
    datatable(df, options = list(dom = "t", paging = FALSE), rownames = FALSE)
  })

  # ── Tab 3: Eosinophil Kinetics ────────────────────────────────────────────
  output$aec_plot <- renderPlotly({
    d <- sim_result()
    plot_ly(d, x = ~time_wk) %>%
      add_lines(y = ~AEC_obs, name = "AEC (cells/µL)", line = list(color = "#E74C3C", width = 2)) %>%
      add_lines(y = rep(1500, nrow(d)), name = "HES threshold (1500)", line = list(dash = "dash", color = "orange")) %>%
      add_lines(y = rep(300,  nrow(d)), name = "Target AEC (300)", line = list(dash = "dot", color = "green")) %>%
      layout(xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "AEC (cells/µL)"),
             legend = list(orientation = "h"))
  })

  output$resp_300_box <- renderValueBox({
    d <- sim_result()
    wk24 <- d %>% filter(abs(time_wk - 24) < 0.6) %>% slice(1)
    val  <- if (wk24$RESP_300 > 0.5) "YES" else "NO"
    valueBox(val, "AEC < 300 at Wk24", icon = icon("check"), color = if (val == "YES") "green" else "red", width = 12)
  })
  output$resp_1500_box <- renderValueBox({
    d <- sim_result()
    wk24 <- d %>% filter(abs(time_wk - 24) < 0.6) %>% slice(1)
    val  <- if (wk24$RESP_1500 > 0.5) "YES" else "NO"
    valueBox(val, "AEC < 1500 at Wk24", icon = icon("check"), color = if (val == "YES") "green" else "orange", width = 12)
  })
  output$aec_interpretation <- renderUI({
    d <- sim_result()
    last <- d %>% slice_tail(n = 1)
    aec  <- last$AEC_obs
    tags$div(
      style = "font-size:13px",
      HTML(paste0(
        "<b>Final AEC:</b> ", round(aec), " cells/µL<br>",
        "<b>Status:</b> ",
        if (aec < 300) "<span style='color:green'>Complete response (< 300/µL)</span>"
        else if (aec < 1500) "<span style='color:orange'>Partial control (< 1500/µL)</span>"
        else "<span style='color:red'>HES active (≥ 1500/µL)</span>"
      ))
    )
  })

  output$bm_plot <- renderPlotly({
    d <- sim_result()
    plot_ly(d, x = ~time_wk) %>%
      add_lines(y = ~EoP,    name = "EoP (progenitors)", line = list(color = "#E74C3C")) %>%
      add_lines(y = ~EoI,    name = "EoI (immature)",    line = list(color = "#F39C12")) %>%
      add_lines(y = ~EoM_BM, name = "EoM-BM (mature)",  line = list(color = "#9B59B6")) %>%
      layout(xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "BM Eosinophil Units"),
             legend = list(orientation = "h"))
  })

  output$il5_pct_plot <- renderPlotly({
    d <- sim_result()
    ay <- list(overlaying = "y", side = "right", title = "% Change AEC")
    plot_ly(d, x = ~time_wk) %>%
      add_lines(y = ~IL5,    name = "IL-5 (pg/mL)", line = list(color = "#27AE60")) %>%
      add_lines(y = ~PCT_CHG, name = "% ΔUE (AEC)", yaxis = "y2", line = list(color = "#3498DB", dash = "dash")) %>%
      layout(yaxis = list(title = "IL-5 (pg/mL)"),
             yaxis2 = ay,
             xaxis  = list(title = "Time (weeks)"),
             legend = list(orientation = "h"))
  })

  # ── Tab 4: Organ Damage ───────────────────────────────────────────────────
  output$cardiac_plot <- renderPlotly({
    d <- sim_result()
    plot_ly(d, x = ~time_wk, y = ~FIBROSIS) %>%
      add_lines(name = "Cardiac Fibrosis", line = list(color = "#E74C3C", width = 2)) %>%
      add_lines(y = rep(0.3, nrow(d)), name = "Moderate threshold", line = list(dash = "dash", color = "orange")) %>%
      add_lines(y = rep(0.6, nrow(d)), name = "Severe threshold", line = list(dash = "dot", color = "red")) %>%
      layout(xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Fibrosis Score (0–1)", range = c(0, 1)),
             legend = list(orientation = "h"))
  })
  output$pulm_plot <- renderPlotly({
    d <- sim_result()
    plot_ly(d, x = ~time_wk, y = ~PULM_SCORE) %>%
      add_lines(name = "Pulmonary Score", line = list(color = "#2E86C1", width = 2)) %>%
      layout(xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Pulmonary Score (0–1)", range = c(0, 1)))
  })

  # ── Tab 5: Comparison ─────────────────────────────────────────────────────
  colors5 <- c(
    "Untreated"                   = "#E74C3C",
    "Prednisolone 70mg"           = "#F39C12",
    "Mepolizumab 300mg q4w"       = "#27AE60",
    "Benralizumab 30mg q4w→q8w"   = "#3498DB",
    "Imatinib 100mg/day (Clonal)" = "#9B59B6"
  )
  output$compare_plot <- renderPlotly({
    d <- compare_result()
    plt <- plot_ly()
    for (sc in unique(d$scenario)) {
      sub <- d %>% filter(scenario == sc)
      plt <- plt %>%
        add_lines(data = sub, x = ~time_wk, y = ~AEC_obs, name = sc,
                  line = list(color = colors5[[sc]], width = 2))
    }
    plt %>%
      add_lines(y = rep(1500, 100), x = seq(0, 52, length.out = 100),
                name = "HES threshold", line = list(dash = "dash", color = "orange")) %>%
      layout(xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "AEC (cells/µL)"),
             legend = list(orientation = "h"))
  })
  output$compare_table <- renderDT({
    d <- compare_result() %>%
      filter(abs(time_wk - 24) < 0.6) %>%
      group_by(scenario) %>% slice(1) %>% ungroup() %>%
      select(scenario, AEC_obs, PCT_CHG, RESP_300, RESP_1500, FIBROSIS, PULM_SCORE, IL5) %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
    datatable(d, options = list(dom = "t", paging = FALSE, scrollX = TRUE), rownames = FALSE)
  })

  # ── Tab 6: Biomarkers ─────────────────────────────────────────────────────
  output$biom_cardiac_plot <- renderPlotly({
    d <- sim_result()
    plot_ly(d, x = ~AEC_obs, y = ~FIBROSIS, color = ~time_wk,
            colors = c("#FFF176", "#E74C3C")) %>%
      add_markers(size = 4) %>%
      layout(xaxis = list(title = "AEC (cells/µL)"),
             yaxis = list(title = "Cardiac Fibrosis Score"))
  })
  output$biom_il5_plot <- renderPlotly({
    d <- sim_result()
    plot_ly(d, x = ~IL5, y = ~AEC_obs, color = ~time_wk,
            colors = c("#B3E5FC", "#01579B")) %>%
      add_markers(size = 4) %>%
      layout(xaxis = list(title = "IL-5 (pg/mL)"),
             yaxis = list(title = "AEC (cells/µL)"))
  })
  output$biom_pulm_plot <- renderPlotly({
    d <- sim_result()
    plot_ly(d, x = ~AEC_obs, y = ~PULM_SCORE, color = ~time_wk,
            colors = c("#E0F7FA", "#006064")) %>%
      add_markers(size = 4) %>%
      layout(xaxis = list(title = "AEC (cells/µL)"),
             yaxis = list(title = "Pulmonary Score"))
  })
  output$biom_multi_plot <- renderPlotly({
    d <- sim_result()
    plot_ly(d, x = ~time_wk) %>%
      add_lines(y = ~(FIBROSIS * 100), name = "Cardiac (×100)", line = list(color = "#E74C3C")) %>%
      add_lines(y = ~(PULM_SCORE * 100), name = "Pulmonary (×100)", line = list(color = "#2E86C1")) %>%
      add_lines(y = ~IL5, name = "IL-5 (pg/mL)", line = list(color = "#27AE60")) %>%
      layout(xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Score / pg/mL"),
             legend = list(orientation = "h"))
  })

} # end server

# ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
