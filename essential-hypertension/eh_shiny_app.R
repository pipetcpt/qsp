## =============================================================================
## Essential Hypertension — Interactive Shiny Dashboard
## =============================================================================
## Tabs (6+):
##   1. Patient Profile & Disease Overview
##   2. Drug Pharmacokinetics (PK)
##   3. RAAS & Biomarker Dynamics
##   4. Hemodynamic Endpoints (SBP, DBP, MAP, HR)
##   5. Treatment Scenario Comparison
##   6. Long-term Organ Protection (LVM, eGFR)
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(patchwork)

## ─── Embed mrgsolve model (shared with eh_mrgsolve_model.R) ──────────────────
code_eh <- '
$PARAM @annotated
  KA_ACEI:1.2:ACEI ka; F_ACEI:0.28:ACEI F; V1_ACEI:8.0:ACEI Vc; V2_ACEI:32.0:ACEI Vp
  CL_ACEI:6.5:ACEI CL; Q_ACEI:3.0:ACEI Q
  KA_ARB:1.1:ARB ka; F_ARB:0.33:ARB F; V1_ARB:14.0:ARB Vc; V2_ARB:45.0:ARB Vp
  CL_ARB:5.2:ARB CL; Q_ARB:3.5:ARB Q
  KA_CCB:0.25:CCB ka; F_CCB:0.64:CCB F; V1_CCB:21.0:CCB Vc; V2_CCB:400.0:CCB Vp
  CL_CCB:3.5:CCB CL; Q_CCB:8.0:CCB Q
  KA_BB:1.3:BB ka; F_BB:0.80:BB F; V1_BB:12.0:BB Vc; V2_BB:100.0:BB Vp
  CL_BB:9.0:BB CL; Q_BB:4.5:BB Q
  KA_HCTZ:1.5:HCTZ ka; F_HCTZ:0.70:HCTZ F; V1_HCTZ:4.0:HCTZ Vc; CL_HCTZ:18.0:HCTZ CL
  ANGII0:15.0:AngII baseline; KPROD_AII:0.15:AngII prod; KDEG_AII:0.15:AngII deg
  IC50_ACEI:0.005:ACEI IC50; IC50_ARB:0.02:ARB IC50
  ALDO0:180.0:Aldo baseline; KPROD_AL:0.12:Aldo prod; KDEG_AL:0.22:Aldo deg
  SNS0:1.0:SNS baseline; KRET_SNS:0.30:SNS rate; IC50_BB:0.10:BB IC50
  NO0:1.0:NO baseline; KPROD_NO:0.50:NO prod; KDEG_NO:0.50:NO deg
  TPR0:1.0:TPR baseline; KTPR_RET:0.08:TPR rate
  ALPHA_AII:0.40:AngII->TPR; ALPHA_SNS:0.30:SNS->TPR; ALPHA_NO:0.20:NO->TPR; ALPHA_CCB:0.25:CCB->TPR
  HR0:70.0:HR baseline; SV0:70.0:SV baseline; ALPHA_HR:0.20:SNS->HR; BETA_BB:0.30:BB->HR
  PV0:3.2:PV baseline; KPVRET:0.06:PV rate; HCTZ_PV:0.10:HCTZ->PV
  MAP0:100.0:MAP baseline; PP0:50.0:PP baseline
  LVM0:210.0:LVM baseline; KLVM_ON:0.002:LVH growth; KLVM_RET:0.0005:LVM regress; MAP_THRESH:93.0:MAP threshold
  EGFR0:72.0:eGFR baseline; KEGFR_DEC:0.0001:eGFR decline; KEGFR_RET:0.0003:eGFR recovery

$CMT ACEI_C ACEI_P ARB_C ARB_P CCB_C CCB_P BB_C BB_P HCTZ_C
     ANGII ALDO SNS_T NO_IDX TPR_N CO_L PV_L LVM_G EGFR_ML

$MAIN
  ANGII_0=ANGII0; ALDO_0=ALDO0; SNS_T_0=SNS0; NO_IDX_0=NO0;
  TPR_N_0=TPR0; CO_L_0=MAP0/80.0; PV_L_0=PV0; LVM_G_0=LVM0; EGFR_ML_0=EGFR0;

$ODE
  double ACEI_Cp=ACEI_C/V1_ACEI;
  dxdt_ACEI_C=-CL_ACEI/V1_ACEI*ACEI_C-Q_ACEI/V1_ACEI*ACEI_C+Q_ACEI/V2_ACEI*ACEI_P;
  dxdt_ACEI_P=Q_ACEI/V1_ACEI*ACEI_C-Q_ACEI/V2_ACEI*ACEI_P;
  double ARB_Cp=ARB_C/V1_ARB;
  dxdt_ARB_C=-CL_ARB/V1_ARB*ARB_C-Q_ARB/V1_ARB*ARB_C+Q_ARB/V2_ARB*ARB_P;
  dxdt_ARB_P=Q_ARB/V1_ARB*ARB_C-Q_ARB/V2_ARB*ARB_P;
  double CCB_Cp=CCB_C/V1_CCB;
  dxdt_CCB_C=-CL_CCB/V1_CCB*CCB_C-Q_CCB/V1_CCB*CCB_C+Q_CCB/V2_CCB*CCB_P;
  dxdt_CCB_P=Q_CCB/V1_CCB*CCB_C-Q_CCB/V2_CCB*CCB_P;
  double BB_Cp=BB_C/V1_BB;
  dxdt_BB_C=-CL_BB/V1_BB*BB_C-Q_BB/V1_BB*BB_C+Q_BB/V2_BB*BB_P;
  dxdt_BB_P=Q_BB/V1_BB*BB_C-Q_BB/V2_BB*BB_P;
  double HCTZ_Cp=HCTZ_C/V1_HCTZ;
  dxdt_HCTZ_C=-CL_HCTZ/V1_HCTZ*HCTZ_C;
  double ACE_inhib=ACEI_Cp/(ACEI_Cp+IC50_ACEI);
  double AT1R_block=ARB_Cp/(ARB_Cp+IC50_ARB);
  double BB_block=BB_Cp/(BB_Cp+IC50_BB);
  double VGCC_block=CCB_Cp/(CCB_Cp+0.003);
  double NCC_inhib=HCTZ_Cp/(HCTZ_Cp+0.02);
  double ACE_activity=1.0-ACE_inhib;
  double AngII_feedback=1.0+0.8*AT1R_block;
  double AngII_prod=KPROD_AII*ANGII0*ACE_activity*AngII_feedback;
  dxdt_ANGII=AngII_prod-KDEG_AII*ANGII;
  double Aldo_stim=(ANGII/ANGII0)*(1.0-AT1R_block);
  dxdt_ALDO=KPROD_AL*ALDO0*Aldo_stim-KDEG_AL*ALDO;
  double SNS_input=SNS0*(1.0-0.15*BB_block);
  dxdt_SNS_T=KRET_SNS*(SNS_input-SNS_T);
  double NO_target=NO0*(1.0-0.3*(ANGII/ANGII0-1.0)+0.4*ACE_inhib);
  if(NO_target<0.1) NO_target=0.1; if(NO_target>2.5) NO_target=2.5;
  dxdt_NO_IDX=KPROD_NO*NO_target-KDEG_NO*NO_IDX;
  double TPR_target=TPR0*(1.0+ALPHA_AII*(ANGII/ANGII0-1.0)+ALPHA_SNS*(SNS_T/SNS0-1.0)
    -ALPHA_NO*(NO_IDX/NO0-1.0)*NO_IDX/NO0-ALPHA_CCB*VGCC_block-AT1R_block*0.15-ACE_inhib*0.12);
  if(TPR_target<0.3) TPR_target=0.3;
  dxdt_TPR_N=KTPR_RET*(TPR_target-TPR_N);
  double HR_current=HR0*(1.0+ALPHA_HR*(SNS_T/SNS0-1.0)-BETA_BB*BB_block);
  if(HR_current<40) HR_current=40;
  double SV_current=SV0*(PV_L/PV0)*(1.0-0.12*(TPR_N-1.0));
  if(SV_current<20) SV_current=20;
  dxdt_CO_L=0.5*((HR_current*SV_current/1000.0)-CO_L);
  double PV_target=PV0*(1.0+0.08*(ALDO/ALDO0-1.0)-HCTZ_PV*NCC_inhib);
  if(PV_target<1.5) PV_target=1.5;
  dxdt_PV_L=KPVRET*(PV_target-PV_L);
  double MAP_calc=CO_L*TPR_N*80.0;
  double LVM_stimulus=(MAP_calc>MAP_THRESH)?(MAP_calc-MAP_THRESH):0.0;
  double LVM_regression=(MAP_calc<MAP_THRESH)?KLVM_RET*(LVM_G-180.0):0.0;
  dxdt_LVM_G=KLVM_ON*LVM_stimulus-LVM_regression;
  double EGFR_decline=(MAP_calc>MAP_THRESH)?KEGFR_DEC*(MAP_calc-MAP_THRESH)*EGFR_ML:0.0;
  double EGFR_recovery=(MAP_calc<=MAP_THRESH)?KEGFR_RET*(EGFR0-EGFR_ML):0.0;
  dxdt_EGFR_ML=-EGFR_decline+EGFR_recovery;

$TABLE
  double MAP_out=CO_L*TPR_N*80.0;
  double ART_proxy=1.0+0.005*(LVM_G-LVM0);
  double PP_cur=PP0*(1.0+0.5*(ART_proxy-1.0));
  capture MAP=MAP_out; capture SBP=MAP_out+PP_cur*2.0/3.0;
  capture DBP=MAP_out-PP_cur/3.0; capture PP=PP_cur;
  capture HR=HR0*(1.0+ALPHA_HR*(SNS_T/SNS0-1.0)-BETA_BB*(BB_C/V1_BB/(BB_C/V1_BB+IC50_BB)));
  capture CO=CO_L; capture TPR=TPR_N; capture PV=PV_L;
  capture AngII=ANGII; capture Aldo=ALDO; capture NO=NO_IDX;
  capture LVM=LVM_G; capture eGFR=EGFR_ML;
  capture Cp_ACEI=ACEI_C/V1_ACEI; capture Cp_ARB=ARB_C/V1_ARB;
  capture Cp_CCB=CCB_C/V1_CCB; capture Cp_BB=BB_C/V1_BB;
  capture Cp_HCTZ=HCTZ_C/V1_HCTZ;
'

mod <- mread("eh_shiny", tempdir(), code_eh)

## ─── Helper: build simulation ────────────────────────────────────────────────
run_sim <- function(dose_acei = 0, dose_arb = 0, dose_ccb = 0,
                    dose_bb = 0, dose_hctz = 0,
                    weeks = 12, delta = 6) {
  dur  <- weeks * 7 * 24
  evts <- ev(ID = 1, amt = 0, cmt = "ACEI_C", time = 0)
  add_dose <- function(ev_list, cmt, dose_mg, F_val) {
    if (dose_mg > 0)
      ev_list <- ev_list + ev(ID = 1, amt = dose_mg * F_val,
                              cmt = cmt, ii = 24,
                              addl = floor(dur / 24) - 1)
    ev_list
  }
  evts <- add_dose(evts, "ACEI_C", dose_acei, 0.28)
  evts <- add_dose(evts, "ARB_C",  dose_arb,  0.33)
  evts <- add_dose(evts, "CCB_C",  dose_ccb,  0.64)
  evts <- add_dose(evts, "BB_C",   dose_bb,   0.80)
  evts <- add_dose(evts, "HCTZ_C", dose_hctz, 0.70)

  mod %>% ev(evts) %>%
    mrgsim(end = dur, delta = delta) %>%
    as_tibble() %>%
    mutate(time_wk = time / 168)
}

## ─── UI ─────────────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "Essential Hypertension QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "tab_patient",  icon = icon("user-md")),
      menuItem("Drug PK",            tabName = "tab_pk",       icon = icon("pills")),
      menuItem("RAAS & Biomarkers",  tabName = "tab_raas",     icon = icon("vials")),
      menuItem("Hemodynamics",       tabName = "tab_hemo",     icon = icon("heartbeat")),
      menuItem("Scenario Comparison",tabName = "tab_scenario", icon = icon("chart-bar")),
      menuItem("Organ Protection",   tabName = "tab_organ",    icon = icon("kidneys"))
    ),
    hr(),
    h5("Drug Doses (mg/day)", style = "padding-left:15px; color:#bbb"),
    sliderInput("dose_acei", "Ramipril (ACEI)", 0, 20,  0, step = 2.5),
    sliderInput("dose_arb",  "Losartan (ARB)",  0, 100, 0, step = 25),
    sliderInput("dose_ccb",  "Amlodipine (CCB)",0, 10,  0, step = 2.5),
    sliderInput("dose_bb",   "Bisoprolol (BB)", 0, 20,  0, step = 2.5),
    sliderInput("dose_hctz", "HCTZ (Diuretic)", 0, 25,  0, step = 6.25),
    sliderInput("sim_weeks", "Simulation (weeks)", 4, 52, 24),
    actionButton("run_sim", "▶ Run Simulation", class = "btn-primary btn-block",
                 style = "margin: 8px 10px; width: 85%")
  ),

  dashboardBody(
    tabItems(

      ## TAB 1 — Patient Profile ───────────────────────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient Characteristics", width = 6, status = "primary",
            selectInput("pt_sex",   "Sex",   c("Male", "Female")),
            sliderInput("pt_age",   "Age (years)", 30, 80, 55),
            sliderInput("pt_bmi",   "BMI (kg/m²)", 18, 45, 28),
            sliderInput("pt_sbp0",  "Pre-treatment SBP (mmHg)", 130, 200, 158),
            sliderInput("pt_dbp0",  "Pre-treatment DBP (mmHg)",  80, 120,  96),
            checkboxGroupInput("pt_risk", "Additional Risk Factors",
              choices = c("Smoking", "Diabetes", "Dyslipidemia",
                          "Sleep Apnea", "Family History", "Chronic Kidney Disease"),
              selected = c("Dyslipidemia")
            )
          ),
          box(title = "ESH/ESC Risk Classification", width = 6, status = "warning",
            tableOutput("risk_class_tbl"),
            hr(),
            h5("Disease Overview: Pathophysiology"),
            p("Essential (primary) hypertension accounts for 90-95% of all hypertension cases.
              Key mechanistic drivers include:"),
            tags$ul(
              tags$li(strong("RAAS over-activation:"), " AngII-AT1R signaling drives vasoconstriction, aldosterone release, and vascular remodeling"),
              tags$li(strong("Sympathetic excess:"), " Norepinephrine increases HR, cardiac output, and renal renin release"),
              tags$li(strong("Endothelial dysfunction:"), " Reduced NO bioavailability and increased ET-1 impair vasodilation"),
              tags$li(strong("Renal Na+ retention:"), " Impaired pressure-natriuresis relationship sustains elevated plasma volume"),
              tags$li(strong("Vascular remodeling:"), " Structural changes (intima-media thickening, arterial stiffness) persist after acute signals")
            )
          )
        ),
        fluidRow(
          box(title = "Antihypertensive Drug Mechanisms", width = 12, status = "info",
            DT::dataTableOutput("drug_mech_tbl")
          )
        )
      ),

      ## TAB 2 — Drug PK ────────────────────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Plasma Concentration–Time Profiles (first 72 h)", width = 12,
            plotlyOutput("pk_plot", height = "450px"))
        ),
        fluidRow(
          box(title = "Drug PK Parameters Summary", width = 6,
            tableOutput("pk_params_tbl")),
          box(title = "Receptor Occupancy (% blocked at 24 h)", width = 6,
            plotlyOutput("receptor_occ_plot", height = "280px"))
        )
      ),

      ## TAB 3 — RAAS & Biomarkers ──────────────────────────────────────────────
      tabItem(tabName = "tab_raas",
        fluidRow(
          box(title = "Angiotensin II Dynamics", width = 6,
            plotlyOutput("angii_plot", height = "320px")),
          box(title = "Aldosterone Dynamics", width = 6,
            plotlyOutput("aldo_plot", height = "320px"))
        ),
        fluidRow(
          box(title = "NO Index (Endothelial Function)", width = 6,
            plotlyOutput("no_plot", height = "300px")),
          box(title = "Plasma Volume & SNS Tone", width = 6,
            plotlyOutput("pv_sns_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "Biomarker Summary at Selected Time", width = 12,
            sliderInput("bm_time", "View at week:", 0, 24, 12),
            tableOutput("biomarker_tbl"))
        )
      ),

      ## TAB 4 — Hemodynamics ───────────────────────────────────────────────────
      tabItem(tabName = "tab_hemo",
        fluidRow(
          valueBoxOutput("vbox_sbp", width = 3),
          valueBoxOutput("vbox_dbp", width = 3),
          valueBoxOutput("vbox_map", width = 3),
          valueBoxOutput("vbox_hr",  width = 3)
        ),
        fluidRow(
          box(title = "SBP / DBP Over Time", width = 8,
            plotlyOutput("bp_plot", height = "380px")),
          box(title = "Hemodynamic Targets", width = 4,
            tags$table(class = "table table-condensed table-bordered",
              tags$thead(tags$tr(tags$th("Target"), tags$th("Goal"))),
              tags$tbody(
                tags$tr(tags$td("SBP"),       tags$td("< 130 mmHg")),
                tags$tr(tags$td("DBP"),       tags$td("< 80 mmHg")),
                tags$tr(tags$td("MAP"),       tags$td("< 93 mmHg")),
                tags$tr(tags$td("HR (BB)"),   tags$td("55–70 bpm")),
                tags$tr(tags$td("LVM index"), tags$td("< 115 g/m²")),
                tags$tr(tags$td("eGFR"),      tags$td("> 60 mL/min"))
              )
            )
          )
        ),
        fluidRow(
          box(title = "Heart Rate & Cardiac Output", width = 6,
            plotlyOutput("hr_co_plot", height = "300px")),
          box(title = "TPR & Pulse Pressure", width = 6,
            plotlyOutput("tpr_pp_plot", height = "300px"))
        )
      ),

      ## TAB 5 — Scenario Comparison ────────────────────────────────────────────
      tabItem(tabName = "tab_scenario",
        fluidRow(
          box(title = "Preset Treatment Scenarios", width = 12, status = "info",
            selectInput("preset", "Select Scenario",
              choices = c(
                "Untreated"                   = "none",
                "ACEI: Ramipril 10 mg"        = "acei",
                "ARB: Losartan 100 mg"        = "arb",
                "CCB: Amlodipine 10 mg"       = "ccb",
                "BB: Bisoprolol 10 mg"        = "bb",
                "HCTZ: 25 mg"                 = "hctz",
                "ACEI+CCB+HCTZ (triple)"      = "triple",
                "ARB+CCB+HCTZ (triple)"       = "triple2"
              ),
              multiple = TRUE,
              selected = c("none", "acei", "ccb", "triple")
            ),
            actionButton("run_compare", "Compare Scenarios", class = "btn-success btn-sm")
          )
        ),
        fluidRow(
          box(title = "SBP Comparison Over 24 Weeks", width = 6,
            plotlyOutput("scen_sbp_plot", height = "360px")),
          box(title = "AngII Comparison", width = 6,
            plotlyOutput("scen_angii_plot", height = "360px"))
        ),
        fluidRow(
          box(title = "Week 12 Endpoint Summary Table", width = 12,
            DT::dataTableOutput("scen_summary_tbl"))
        )
      ),

      ## TAB 6 — Organ Protection ───────────────────────────────────────────────
      tabItem(tabName = "tab_organ",
        fluidRow(
          box(title = "LV Mass (Hypertrophy Index) — 24 weeks", width = 6,
            plotlyOutput("lvm_plot", height = "350px")),
          box(title = "eGFR (Renal Function) — 24 weeks", width = 6,
            plotlyOutput("egfr_plot", height = "350px"))
        ),
        fluidRow(
          box(title = "Cardiovascular Risk Reduction Estimate", width = 6,
            plotlyOutput("cv_risk_plot", height = "300px")),
          box(title = "Organ Damage Markers", width = 6, status = "warning",
            p(strong("LV Hypertrophy (LVH):"), " LV mass index > 115 g/m² (men) / 95 g/m² (women) signals target organ damage. Each 10 mmHg SBP reduction associated with 10–15% LVM regression."),
            hr(),
            p(strong("Chronic Kidney Disease (CKD):"), " MAP > 93 mmHg drives glomerular hypertension → proteinuria → eGFR decline. RAAS blockade (ACEI/ARB) reduces proteinuria independently of BP reduction."),
            hr(),
            p(strong("Stroke risk:"), " Every 20/10 mmHg rise in SBP/DBP doubles cardiovascular mortality (Lewington et al., Lancet 2002)."),
            hr(),
            p(strong("Retinopathy:"), " Arteriolar changes (cotton-wool spots, AV nicking) reflect chronic arterial hypertension and correlate with cerebrovascular risk.")
          )
        )
      )

    ) # tabItems
  ) # dashboardBody
) # dashboardPage


## ─── SERVER ─────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  ## Reactive: run custom simulation
  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Running QSP simulation...", {
      run_sim(
        dose_acei = input$dose_acei,
        dose_arb  = input$dose_arb,
        dose_ccb  = input$dose_ccb,
        dose_bb   = input$dose_bb,
        dose_hctz = input$dose_hctz,
        weeks     = input$sim_weeks,
        delta     = 6
      )
    })
  }, ignoreNULL = FALSE)

  ## Reactive: pre-defined scenarios
  scen_data <- eventReactive(input$run_compare, {
    preset_map <- list(
      none   = c(0,   0,   0,  0,  0),
      acei   = c(10,  0,   0,  0,  0),
      arb    = c(0,   100, 0,  0,  0),
      ccb    = c(0,   0,   10, 0,  0),
      bb     = c(0,   0,   0,  10, 0),
      hctz   = c(0,   0,   0,  0,  25),
      triple = c(5,   0,   5,  0,  12.5),
      triple2= c(0,   50,  5,  0,  12.5)
    )
    scen_labels <- c(
      none    = "Untreated",
      acei    = "ACEI (Ramipril 10)",
      arb     = "ARB (Losartan 100)",
      ccb     = "CCB (Amlodipine 10)",
      bb      = "BB (Bisoprolol 10)",
      hctz    = "HCTZ 25 mg",
      triple  = "ACEI+CCB+HCTZ",
      triple2 = "ARB+CCB+HCTZ"
    )

    chosen <- input$preset
    if (is.null(chosen) || length(chosen) == 0) chosen <- c("none", "acei", "triple")

    withProgress(message = "Running scenario comparison...", {
      bind_rows(lapply(chosen, function(s) {
        d <- preset_map[[s]]
        run_sim(d[1], d[2], d[3], d[4], d[5], weeks = 24, delta = 6) %>%
          mutate(scenario = scen_labels[[s]])
      }))
    })
  }, ignoreNULL = FALSE)

  ## ── TAB 1 outputs ──────────────────────────────────────────────────────────
  output$risk_class_tbl <- renderTable({
    sbp <- input$pt_sbp0; dbp <- input$pt_dbp0
    stage <- dplyr::case_when(
      sbp < 130 & dbp < 80 ~ "Normal",
      sbp < 130 & dbp < 80 ~ "Elevated",
      sbp < 140 | dbp < 90 ~ "Stage 1 Hypertension",
      sbp < 160 | dbp < 100~ "Stage 2 Hypertension",
      TRUE ~ "Stage 3 Hypertension"
    )
    n_risks <- length(input$pt_risk)
    cv_risk <- dplyr::case_when(
      stage == "Normal" ~ "Low",
      stage == "Elevated" & n_risks < 3 ~ "Moderate",
      stage %in% c("Stage 1 Hypertension") & n_risks < 3 ~ "Moderate",
      n_risks >= 3 ~ "High",
      TRUE ~ "Very High"
    )
    data.frame(
      Parameter  = c("BP Stage", "Risk Factors", "CV Risk Category",
                     "Target SBP", "Recommended Drug"),
      Value      = c(stage, as.character(n_risks), cv_risk,
                     "< 130 mmHg",
                     ifelse(n_risks >= 2, "ACEI/ARB + CCB ± HCTZ", "ACEI or ARB monotherapy"))
    )
  })

  output$drug_mech_tbl <- DT::renderDataTable({
    DT::datatable(data.frame(
      Drug_Class  = c("ACE Inhibitor","ARB","CCB (Dihydropyridine)","Beta-blocker","Thiazide Diuretic"),
      Example     = c("Ramipril","Losartan","Amlodipine","Bisoprolol","HCTZ"),
      Mechanism   = c(
        "Inhibits ACE → ↓AngII → ↓AT1R → vasodilation; ↑bradykinin → ↑NO",
        "Blocks AT1R → ↓AngII vasoconstriction, ↓aldosterone; reactive ↑renin",
        "Blocks L-type VGCC → ↓Ca²⁺ in VSM → vasodilation; no reflex tachy",
        "Blocks β1-AR → ↓HR, ↓CO, ↓renin; reduces sympathetic cardiac output",
        "Inhibits NCC → ↓Na⁺ reabsorption in DCT → ↓plasma volume → ↓CO"
      ),
      SBP_Reduction = c("−8 to −12","−8 to −12","−8 to −14","−6 to −10","−7 to −11"),
      Special_Benefit = c("CKD, proteinuria, HF","CKD, proteinuria, HF","ISH, angina","AF, HF, post-MI","Resistant HTN, edema")
    ), options = list(pageLength = 5, dom = "t"), rownames = FALSE)
  })

  ## ── TAB 2 — PK ─────────────────────────────────────────────────────────────
  output$pk_plot <- renderPlotly({
    d <- sim_data() %>% filter(time <= 72)
    pk_long <- d %>%
      select(time, Cp_ACEI, Cp_ARB, Cp_CCB, Cp_BB, Cp_HCTZ) %>%
      pivot_longer(-time, names_to = "Drug", values_to = "Cp") %>%
      filter(Cp > 1e-6) %>%
      mutate(Drug = recode(Drug,
        Cp_ACEI = "Ramiprilat (ACEI)",
        Cp_ARB  = "EXP3174 (ARB)",
        Cp_CCB  = "Amlodipine (CCB)",
        Cp_BB   = "Bisoprolol (BB)",
        Cp_HCTZ = "HCTZ"
      ))
    p <- ggplot(pk_long, aes(time, Cp, color = Drug)) +
      geom_line(size = 0.9) +
      labs(x = "Time (h)", y = "Plasma Conc. (mg/L)", title = "Drug PK — First 72 h") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_params_tbl <- renderTable({
    data.frame(
      Drug   = c("Ramiprilat","EXP3174","Amlodipine","Bisoprolol","HCTZ"),
      Dose   = c(paste0(input$dose_acei," mg"),paste0(input$dose_arb," mg"),
                 paste0(input$dose_ccb," mg"),paste0(input$dose_bb," mg"),
                 paste0(input$dose_hctz," mg")),
      t_half = c("9–18 h","6–9 h","30–50 h","9–12 h","5–15 h"),
      Tmax   = c("1–2 h","6–8 h","6–12 h","2–4 h","1–5 h"),
      Vd     = c("~56 L","~34 L","~21 L/kg","~3.2 L/kg","~0.8 L/kg"),
      F_oral = c("28%","33%","64%","80%","70%")
    )
  })

  output$receptor_occ_plot <- renderPlotly({
    d <- sim_data() %>% filter(time == max(time))
    vals <- c(
      ACE_inhib  = 100 * (d$Cp_ACEI / (d$Cp_ACEI + 0.005)),
      AT1R_block = 100 * (d$Cp_ARB  / (d$Cp_ARB  + 0.02 )),
      VGCC_block = 100 * (d$Cp_CCB  / (d$Cp_CCB  + 0.003)),
      BB_block   = 100 * (d$Cp_BB   / (d$Cp_BB   + 0.10 )),
      NCC_inhib  = 100 * (d$Cp_HCTZ / (d$Cp_HCTZ + 0.02 ))
    )
    df <- data.frame(Target = names(vals), Occupancy = round(vals, 1))
    p <- ggplot(df, aes(Target, Occupancy, fill = Target)) +
      geom_col() + ylim(0, 100) +
      labs(y = "Receptor occupancy (%)", title = "At End of Simulation") +
      theme_bw() + theme(legend.position = "none")
    ggplotly(p)
  })

  ## ── TAB 3 — RAAS ───────────────────────────────────────────────────────────
  output$angii_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_wk, AngII)) + geom_line(color = "#1E88E5", size=0.9) +
      geom_hline(yintercept = 15, linetype = "dashed") +
      labs(x = "Weeks", y = "AngII (pg/mL)", title = "Angiotensin II") + theme_bw()
    ggplotly(p)
  })
  output$aldo_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_wk, Aldo)) + geom_line(color = "#FB8C00", size=0.9) +
      labs(x = "Weeks", y = "Aldosterone (pmol/L)") + theme_bw()
    ggplotly(p)
  })
  output$no_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_wk, NO)) + geom_line(color = "#43A047", size=0.9) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(x = "Weeks", y = "NO index (normalized)") + theme_bw()
    ggplotly(p)
  })
  output$pv_sns_plot <- renderPlotly({
    d <- sim_data() %>%
      select(time_wk, PV, SNS_T) %>%
      pivot_longer(-time_wk)
    p <- ggplot(d, aes(time_wk, value, color = name)) + geom_line(size=0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Weeks") + theme_bw() + theme(legend.position = "none")
    ggplotly(p)
  })
  output$biomarker_tbl <- renderTable({
    target_wk <- input$bm_time
    d <- sim_data() %>%
      filter(abs(time_wk - target_wk) == min(abs(time_wk - target_wk))) %>%
      slice(1)
    data.frame(
      Biomarker   = c("AngII (pg/mL)","Aldosterone (pmol/L)","NO Index","Plasma Volume (L)","SNS Tone","TPR (norm.)"),
      Value       = round(c(d$AngII, d$Aldo, d$NO, d$PV, d$SNS_T, d$TPR), 2),
      Reference   = c("8–25","110–860","1.0 (normal)","3.0–3.5","1.0 (normal)","1.0 (normal)")
    )
  })

  ## ── TAB 4 — Hemodynamics ───────────────────────────────────────────────────
  last_vals <- reactive({
    d <- sim_data()
    d[nrow(d), ]
  })
  output$vbox_sbp <- renderValueBox({
    v <- round(last_vals()$SBP, 0)
    valueBox(paste0(v, " mmHg"), "SBP (final)",
      icon = icon("tachometer-alt"), color = ifelse(v > 130, "red", "green"))
  })
  output$vbox_dbp <- renderValueBox({
    v <- round(last_vals()$DBP, 0)
    valueBox(paste0(v, " mmHg"), "DBP (final)",
      icon = icon("tachometer-alt"), color = ifelse(v > 80, "yellow", "green"))
  })
  output$vbox_map <- renderValueBox({
    v <- round(last_vals()$MAP, 0)
    valueBox(paste0(v, " mmHg"), "MAP (final)",
      icon = icon("heart"), color = ifelse(v > 93, "red", "green"))
  })
  output$vbox_hr <- renderValueBox({
    v <- round(last_vals()$HR, 0)
    valueBox(paste0(v, " bpm"), "HR (final)",
      icon = icon("heartbeat"), color = ifelse(v > 80, "yellow", "green"))
  })
  output$bp_plot <- renderPlotly({
    d <- sim_data()
    d2 <- d %>% select(time_wk, SBP, DBP) %>%
      pivot_longer(-time_wk, names_to = "Measure", values_to = "mmHg")
    p <- ggplot(d2, aes(time_wk, mmHg, color = Measure)) + geom_line(size=0.9) +
      geom_hline(aes(yintercept = yint, linetype = lbl),
        data = data.frame(yint = c(130,80), lbl = c("SBP target","DBP target"))) +
      scale_color_manual(values = c(SBP="#E53935", DBP="#1E88E5")) +
      labs(x = "Weeks", y = "mmHg", title = "Blood Pressure Over Time") + theme_bw()
    ggplotly(p)
  })
  output$hr_co_plot <- renderPlotly({
    d <- sim_data() %>% select(time_wk, HR, CO) %>%
      pivot_longer(-time_wk, names_to = "Measure", values_to = "Value")
    p <- ggplot(d, aes(time_wk, Value, color = Measure)) + geom_line(size=0.9) +
      facet_wrap(~Measure, scales = "free_y") +
      labs(x = "Weeks") + theme_bw() + theme(legend.position = "none")
    ggplotly(p)
  })
  output$tpr_pp_plot <- renderPlotly({
    d <- sim_data() %>% select(time_wk, TPR, PP) %>%
      pivot_longer(-time_wk, names_to = "Measure", values_to = "Value")
    p <- ggplot(d, aes(time_wk, Value, color = Measure)) + geom_line(size=0.9) +
      facet_wrap(~Measure, scales = "free_y") +
      labs(x = "Weeks") + theme_bw() + theme(legend.position = "none")
    ggplotly(p)
  })

  ## ── TAB 5 — Scenarios ──────────────────────────────────────────────────────
  output$scen_sbp_plot <- renderPlotly({
    d <- scen_data()
    p <- ggplot(d, aes(time_wk, SBP, color = scenario)) + geom_line(size=0.8) +
      geom_hline(yintercept = 130, linetype = "dashed") +
      labs(x="Weeks", y="SBP (mmHg)", title="SBP by Scenario") + theme_bw()
    ggplotly(p)
  })
  output$scen_angii_plot <- renderPlotly({
    d <- scen_data()
    p <- ggplot(d, aes(time_wk, AngII, color = scenario)) + geom_line(size=0.8) +
      labs(x="Weeks", y="AngII (pg/mL)", title="Angiotensin II by Scenario") + theme_bw()
    ggplotly(p)
  })
  output$scen_summary_tbl <- DT::renderDataTable({
    d <- scen_data() %>%
      group_by(scenario) %>%
      filter(abs(time_wk - 12) == min(abs(time_wk - 12))) %>%
      slice(1) %>%
      select(scenario, SBP, DBP, MAP, HR, CO, AngII, Aldo, eGFR, LVM) %>%
      mutate(across(where(is.numeric), ~ round(., 1)))
    DT::datatable(d, rownames = FALSE,
      options = list(pageLength = 10, dom = "t"))
  })

  ## ── TAB 6 — Organ Protection ───────────────────────────────────────────────
  output$lvm_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_wk, LVM)) + geom_line(color="#8E24AA", size=0.9) +
      geom_hline(yintercept = 200, linetype="dashed", color="gray40") +
      labs(x="Weeks", y="LV Mass (g)", title="Left Ventricular Mass") + theme_bw()
    ggplotly(p)
  })
  output$egfr_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_wk, eGFR)) + geom_line(color="#00ACC1", size=0.9) +
      geom_hline(yintercept = 60, linetype="dashed", color="gray40") +
      labs(x="Weeks", y="eGFR (mL/min/1.73m²)", title="Renal Function (eGFR)") + theme_bw()
    ggplotly(p)
  })
  output$cv_risk_plot <- renderPlotly({
    d <- sim_data()
    d2 <- d %>%
      mutate(rel_risk = exp(0.05 * (SBP - 115) / 20)) %>%
      filter(time_wk %in% seq(0, 24, by = 4))
    p <- ggplot(d2, aes(time_wk, rel_risk)) +
      geom_line(color="#E53935", size=0.9) + geom_point(color="#E53935") +
      labs(x="Weeks", y="Relative CV Risk", title="Estimated CV Risk (SBP-based)") +
      theme_bw()
    ggplotly(p)
  })
}

shinyApp(ui, server)
