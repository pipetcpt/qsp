## =============================================================================
## ATTR Amyloidosis QSP Shiny Dashboard
## Interactive Simulation: TTR Stabilizers · siRNA · ASO · Natural History
## =============================================================================
library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)
library(scales)

## ─── Inline model (minimal version for Shiny) ────────────────────────────────
ATTR_CODE <- '
$PARAM
ATTR_type=1, WT=75, AGE=72,
kin_mRNA=0.01440, kout_mRNA=0.01440,
ka_TAF=0.42, CL_TAF=0.96, Vc_TAF=16.0, Vp_TAF=32.0, Q_TAF=1.20, F_TAF=0.99,
Emax_stab=0.80, EC50_stab=0.80, n_stab=1.50,
ka_VUT=0.080, CL_VUT=0.120, Vc_VUT=5.80, F_VUT=0.820,
Imax_VUT=0.830, IC50_VUT=0.450,
ka_INO=0.100, CL_INO=0.040, Vc_INO=78.0, F_INO=0.700,
Imax_INO=0.720, IC50_INO=28.0,
CL_PAT=0.180, Vc_PAT=3.30, Vp_PAT=5.60, Q_PAT=0.960,
Imax_PAT=0.800, IC50_PAT=0.600,
ksyn_TET=0.180, kout_TET=0.0098, kdis_TET=0.020,
kagg_MONO=0.100, kdeg_MONO=0.150,
kfib_OLIGO=0.050, kdeg_OLIGO=0.030, kdeg_FIB=0.0001,
f_heart_wt=0.400, f_nerve_wt=0.050,
f_heart_v=0.200, f_nerve_v=0.350,
kdet_EF=0.00015, krec_EF=0.00050, LVEF_base=62.0, LVEF_min=20.0,
kin_BNP=0.500, kout_BNP=0.0800, BNP_base=200.0,
kin_inflam=0.0100, kout_inflam=0.0500, FIB50_inf=1.50,
kin_NIS=0.00200, kout_NIS=0.000050, NIS_base=5.0,
kdet_mBMI=0.000040, NIS50_mBMI=50.0, mBMI_base=1000.0,
kdet_eGFR=0.000020, eGFR_base=72.0

$INIT
A_TAF_GUT=0, A_TAF_C=0, A_TAF_P=0,
A_VUT_SC=0, A_VUT_C=0,
A_INO_SC=0, A_INO_C=0,
A_PAT_C=0, A_PAT_P=0,
TTR_MRNA=1.0, TTR_TET=7.78, TTR_MONO=0.50, TTR_OLIGO=0.01,
FIB_HRT=0.10, FIB_NRV=0.05, FIB_SYS=0.05,
INFLAM=0.012, LVEF=62.0, NT_proBNP=200.0,
NIS=5.0, mBMI=1000.0, eGFR=72.0

$ODE
double f_heart = (ATTR_type==1) ? 0.400 : 0.200;
double f_nerve = (ATTR_type==1) ? 0.050 : 0.350;
double f_sys   = 1.0 - f_heart - f_nerve;
if(f_sys<0.05) f_sys=0.05;

dxdt_A_TAF_GUT = -ka_TAF*A_TAF_GUT;
double C_TAF   = A_TAF_C/Vc_TAF;
dxdt_A_TAF_C   = F_TAF*ka_TAF*A_TAF_GUT + Q_TAF*(A_TAF_P/Vp_TAF-C_TAF) - (CL_TAF/Vc_TAF)*A_TAF_C;
dxdt_A_TAF_P   = Q_TAF*(C_TAF - A_TAF_P/Vp_TAF);

dxdt_A_VUT_SC  = -ka_VUT*A_VUT_SC;
double C_VUT   = A_VUT_C/Vc_VUT;
dxdt_A_VUT_C   = F_VUT*ka_VUT*A_VUT_SC - (CL_VUT/Vc_VUT)*A_VUT_C;

dxdt_A_INO_SC  = -ka_INO*A_INO_SC;
double C_INO   = A_INO_C/Vc_INO;
dxdt_A_INO_C   = F_INO*ka_INO*A_INO_SC - (CL_INO/Vc_INO)*A_INO_C;

double C_PAT   = A_PAT_C/Vc_PAT;
dxdt_A_PAT_C   = R_PAT + Q_PAT*(A_PAT_P/Vp_PAT-C_PAT) - (CL_PAT/Vc_PAT)*A_PAT_C;
dxdt_A_PAT_P   = Q_PAT*(C_PAT - A_PAT_P/Vp_PAT);

double E_stab  = Emax_stab*pow(C_TAF,n_stab)/(pow(EC50_stab,n_stab)+pow(C_TAF,n_stab));
if(E_stab>0.95) E_stab=0.95;
double E_VUT   = Imax_VUT*C_VUT/(IC50_VUT+C_VUT);
double E_PAT   = Imax_PAT*C_PAT/(IC50_PAT+C_PAT);
double E_INO   = Imax_INO*C_INO/(IC50_INO+C_INO);
double E_RNA   = E_VUT+E_PAT+E_INO;
if(E_RNA>0.95) E_RNA=0.95;

dxdt_TTR_MRNA  = kin_mRNA*(1.0-E_RNA) - kout_mRNA*TTR_MRNA;
dxdt_TTR_TET   = ksyn_TET*TTR_MRNA - kdis_TET*(1.0-E_stab)*TTR_TET - kout_TET*TTR_TET;
dxdt_TTR_MONO  = 2.0*kdis_TET*(1.0-E_stab)*TTR_TET - kagg_MONO*TTR_MONO - kdeg_MONO*TTR_MONO;
dxdt_TTR_OLIGO = kagg_MONO*TTR_MONO - kfib_OLIGO*TTR_OLIGO - kdeg_OLIGO*TTR_OLIGO;
dxdt_FIB_HRT   = kfib_OLIGO*TTR_OLIGO*f_heart - kdeg_FIB*FIB_HRT;
dxdt_FIB_NRV   = kfib_OLIGO*TTR_OLIGO*f_nerve - kdeg_FIB*FIB_NRV;
dxdt_FIB_SYS   = kfib_OLIGO*TTR_OLIGO*f_sys   - kdeg_FIB*FIB_SYS;
dxdt_INFLAM    = kin_inflam*FIB_HRT/(FIB50_inf+FIB_HRT) - kout_inflam*INFLAM;
double EF_det  = kdet_EF*FIB_HRT*(1.0+INFLAM)*LVEF;
double EF_rec  = krec_EF*(LVEF_base-LVEF);
dxdt_LVEF      = (LVEF>LVEF_min) ? EF_rec-EF_det : 0;
dxdt_NT_proBNP = kin_BNP*(INFLAM+1.0/(LVEF/100.0+0.01)) - kout_BNP*NT_proBNP;
dxdt_NIS       = kin_NIS*FIB_NRV - kout_NIS*NIS;
dxdt_mBMI      = -kdet_mBMI*(NIS/(NIS+NIS50_mBMI))*mBMI;
dxdt_eGFR      = -kdet_eGFR*FIB_SYS*eGFR;

$TABLE
capture C_TAF_ug=A_TAF_C/Vc_TAF, C_VUT_ng=(A_VUT_C/Vc_VUT)*1000,
        C_INO_ng=(A_INO_C/Vc_INO)*1000, C_PAT_ug=A_PAT_C/Vc_PAT,
        E_stab_pct=E_stab*100, E_RNA_pct=E_RNA*100,
        TTR_mRNA_rel=TTR_MRNA, TTR_TET_AU=TTR_TET,
        FIB_HRT_AU=FIB_HRT, FIB_NRV_AU=FIB_NRV,
        LVEF_pct=LVEF, NTproBNP=NT_proBNP, NIS_score=NIS,
        mBMI_val=mBMI, eGFR_val=eGFR, INFLAM_lvl=INFLAM
'

mod <- mcode("ATTR_shiny", ATTR_CODE, quiet=TRUE)

## ─── Helpers ─────────────────────────────────────────────────────────────────
build_ev <- function(drug, dose_mg, attr_type, duration_yr) {
  dur_h  <- duration_yr * 8760
  switch(drug,
    "None" = ev(time=0, amt=0, cmt=1, ii=0, addl=0),
    "Tafamidis 61mg QD" = ev(time=0, amt=61, cmt=1, ii=24,
                              addl=ceiling(dur_h/24)-1),
    "Vutrisiran 25mg Q3M" = ev(time=0, amt=25, cmt=4, ii=2184,
                                addl=ceiling(dur_h/2184)-1),
    "Inotersen 300mg QW" = ev(time=0, amt=300, cmt=6, ii=168,
                               addl=ceiling(dur_h/168)-1),
    "Patisiran 0.3mg/kg Q3W" = ev(time=0, amt=0.3*75, cmt=8, rate=-2,
                                   ii=504, addl=ceiling(dur_h/504)-1)
  )
}

run_sim <- function(ev_obj, attr_type, duration_yr) {
  dur_h <- duration_yr * 8760
  mod %>%
    param(ATTR_type = attr_type) %>%
    ev(ev_obj) %>%
    mrgsim(end=dur_h, delta=24) %>%
    as_tibble() %>%
    mutate(time_yr = time/8760)
}

## ─────────────────────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(
    title = tags$span("ATTR Amyloidosis QSP", style="font-size:15px; font-weight:bold;"),
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      id = "tabs",
      menuItem("Patient Profile",      tabName="tab_patient",   icon=icon("user")),
      menuItem("Drug PK",              tabName="tab_pk",        icon=icon("pills")),
      menuItem("TTR Misfolding",       tabName="tab_ttr",       icon=icon("dna")),
      menuItem("Cardiac Outcomes",     tabName="tab_cardiac",   icon=icon("heartbeat")),
      menuItem("Neurological Outcomes",tabName="tab_neuro",     icon=icon("brain")),
      menuItem("Scenario Comparison",  tabName="tab_compare",   icon=icon("chart-bar")),
      menuItem("Biomarker Dashboard",  tabName="tab_biomarker", icon=icon("chart-line")),
      menuItem("Model Info",           tabName="tab_info",      icon=icon("info-circle"))
    ),
    tags$hr(),
    tags$div(style="padding:10px;",
      h5("Simulation Settings", style="color:#ecf0f1; font-weight:bold;"),
      selectInput("attr_type", "ATTR Phenotype:",
                  choices=c("ATTRwt (Cardiac)"="1","ATTRv (Neuropathic)"="2")),
      selectInput("drug1", "Treatment:",
                  choices=c("None","Tafamidis 61mg QD","Vutrisiran 25mg Q3M",
                            "Inotersen 300mg QW","Patisiran 0.3mg/kg Q3W")),
      sliderInput("duration_yr","Simulation (years):", min=1, max=5, value=3, step=0.5),
      sliderInput("pt_age","Patient Age:",min=50,max=85,value=72,step=1),
      sliderInput("pt_wt","Weight (kg):",min=40,max=120,value=75,step=5),
      actionButton("run_btn","▶ Run Simulation",
                   class="btn-success", style="width:100%; margin-top:8px;")
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper,.right-side{background-color:#f4f6f9;}
      .box-header{background:#2c3e50 !important; color:white !important;}
      .box-header h3{color:white !important;}
      .value-box .inner{color:white;}
    "))),

    tabItems(

      ## ── Tab 1: Patient Profile ─────────────────────────────────────────────
      tabItem(tabName="tab_patient",
        fluidRow(
          valueBoxOutput("vbox_disease", width=3),
          valueBoxOutput("vbox_drug",    width=3),
          valueBoxOutput("vbox_lvef",    width=3),
          valueBoxOutput("vbox_nis",     width=3)
        ),
        fluidRow(
          box(title="Patient & Disease Overview", width=6, status="primary",
            tags$table(class="table table-condensed",
              tags$thead(tags$tr(tags$th("Parameter"), tags$th("Value"), tags$th("Unit"))),
              tags$tbody(
                tags$tr(tags$td("ATTR Phenotype"), tags$td(textOutput("pt_phenotype",inline=TRUE)), tags$td("—")),
                tags$tr(tags$td("Age"), tags$td(textOutput("pt_age_disp",inline=TRUE)), tags$td("years")),
                tags$tr(tags$td("Weight"), tags$td(textOutput("pt_wt_disp",inline=TRUE)), tags$td("kg")),
                tags$tr(tags$td("Treatment"), tags$td(textOutput("pt_drug_disp",inline=TRUE)), tags$td("—")),
                tags$tr(tags$td("Duration"), tags$td(textOutput("pt_dur_disp",inline=TRUE)), tags$td("years")),
                tags$tr(tags$td("Sim timepoints"), tags$td(textOutput("pt_ntp",inline=TRUE)), tags$td("—"))
              )
            )
          ),
          box(title="ATTR Disease Pathophysiology", width=6, status="info",
            tags$ul(
              tags$li(tags$b("Liver TTR synthesis:"), " ~10 μg/mL plasma TTR; t½~72h"),
              tags$li(tags$b("Tetramer dissociation:"), " Rate-limiting step; ATTRv mutations accelerate"),
              tags$li(tags$b("Toxic oligomers:"), " Soluble pre-fibrillar species; NLRP3/IL-1β"),
              tags$li(tags$b("Cardiac (ATTRwt):"), " Biventricular wall thickening; HFpEF→HFrEF"),
              tags$li(tags$b("Neuropathic (ATTRv):"), " Length-dependent sensorimotor polyneuropathy"),
              tags$li(tags$b("Drug targets:"), " TTR synthesis (siRNA/ASO), tetramer stability, fibrils")
            )
          )
        ),
        fluidRow(
          box(title="Biomarker Trajectory (quick view)", width=12, status="success",
              plotlyOutput("profile_plot", height="320px"))
        )
      ),

      ## ── Tab 2: Drug PK ────────────────────────────────────────────────────
      tabItem(tabName="tab_pk",
        fluidRow(
          box(title="Tafamidis Plasma Concentration (μg/mL)", width=6, status="primary",
              plotlyOutput("pk_taf_plot", height="260px")),
          box(title="RNA-Based Therapy Concentration", width=6, status="success",
              plotlyOutput("pk_rna_plot", height="260px"))
        ),
        fluidRow(
          box(title="Drug Effect: Stabilization vs mRNA Knockdown", width=6, status="warning",
              plotlyOutput("effect_plot", height="260px")),
          box(title="TTR Tetramer Level (Treatment Effect)", width=6, status="danger",
              plotlyOutput("tet_plot", height="260px"))
        )
      ),

      ## ── Tab 3: TTR Misfolding ─────────────────────────────────────────────
      tabItem(tabName="tab_ttr",
        fluidRow(
          box(title="TTR mRNA (% baseline)", width=6, status="primary",
              plotlyOutput("mrna_plot", height="280px")),
          box(title="TTR Tetramer & Monomer Levels (AU)", width=6, status="info",
              plotlyOutput("tet_mono_plot", height="280px"))
        ),
        fluidRow(
          box(title="Toxic Oligomers & Fibril Accumulation (AU)", width=6, status="danger",
              plotlyOutput("fibril_plot", height="280px")),
          box(title="Fibril Distribution by Organ", width=6, status="warning",
              plotlyOutput("fibril_dist_plot", height="280px"))
        )
      ),

      ## ── Tab 4: Cardiac Outcomes ───────────────────────────────────────────
      tabItem(tabName="tab_cardiac",
        fluidRow(
          valueBoxOutput("vbox_lvef2",   width=3),
          valueBoxOutput("vbox_bnp",     width=3),
          valueBoxOutput("vbox_inflam",  width=3),
          valueBoxOutput("vbox_fibhrt",  width=3)
        ),
        fluidRow(
          box(title="LVEF Over Time (%)", width=6, status="danger",
              plotlyOutput("lvef_plot", height="280px")),
          box(title="NT-proBNP Over Time (pg/mL)", width=6, status="warning",
              plotlyOutput("bnp_plot", height="280px"))
        ),
        fluidRow(
          box(title="Cardiac Fibril Load & Inflammation", width=6, status="primary",
              plotlyOutput("card_fibril_inflam_plot", height="280px")),
          box(title="Projected NYHA Class Trajectory", width=6, status="info",
              plotlyOutput("nyha_plot", height="280px"))
        )
      ),

      ## ── Tab 5: Neurological Outcomes ──────────────────────────────────────
      tabItem(tabName="tab_neuro",
        fluidRow(
          valueBoxOutput("vbox_nis2",  width=3),
          valueBoxOutput("vbox_mbmi",  width=3),
          valueBoxOutput("vbox_egfr",  width=3),
          valueBoxOutput("vbox_fibnrv",width=3)
        ),
        fluidRow(
          box(title="NIS Score Over Time (0–244)", width=6, status="danger",
              plotlyOutput("nis_plot", height="280px")),
          box(title="Modified BMI Over Time (g/m²)", width=6, status="warning",
              plotlyOutput("mbmi_plot", height="280px"))
        ),
        fluidRow(
          box(title="Nerve Fibril Load Over Time (AU)", width=6, status="primary",
              plotlyOutput("nrv_fibril_plot", height="280px")),
          box(title="Renal Function (eGFR, mL/min/1.73m²)", width=6, status="info",
              plotlyOutput("egfr_plot", height="280px"))
        )
      ),

      ## ── Tab 6: Scenario Comparison ────────────────────────────────────────
      tabItem(tabName="tab_compare",
        fluidRow(
          box(title="Select Scenarios to Compare", width=12, status="primary",
            checkboxGroupInput("compare_drugs","",
              choices=c("None (Natural History)",
                        "Tafamidis 61mg QD",
                        "Vutrisiran 25mg Q3M",
                        "Inotersen 300mg QW",
                        "Patisiran 0.3mg/kg Q3W"),
              selected=c("None (Natural History)","Tafamidis 61mg QD","Vutrisiran 25mg Q3M"),
              inline=TRUE)
          )
        ),
        fluidRow(
          box(title="LVEF Comparison", width=6, status="danger",
              plotlyOutput("comp_lvef", height="260px")),
          box(title="NIS Score Comparison", width=6, status="warning",
              plotlyOutput("comp_nis", height="260px"))
        ),
        fluidRow(
          box(title="TTR mRNA % Reduction Comparison", width=6, status="success",
              plotlyOutput("comp_mrna", height="260px")),
          box(title="Comparison Summary Table (at 18 months)", width=6, status="info",
              DTOutput("compare_table"))
        )
      ),

      ## ── Tab 7: Biomarker Dashboard ────────────────────────────────────────
      tabItem(tabName="tab_biomarker",
        fluidRow(
          box(title="Cardiac Biomarker Panel", width=6, status="danger",
              plotlyOutput("bm_cardiac", height="320px")),
          box(title="Neurological Biomarker Panel", width=6, status="primary",
              plotlyOutput("bm_neuro", height="320px"))
        ),
        fluidRow(
          box(title="Biomarker Summary at Key Timepoints", width=12, status="success",
              DTOutput("bm_summary_table"))
        )
      ),

      ## ── Tab 8: Model Info ─────────────────────────────────────────────────
      tabItem(tabName="tab_info",
        fluidRow(
          box(title="Model Overview", width=6, status="primary",
            tags$ul(
              tags$li(tags$b("Disease:"), " Transthyretin Amyloidosis (ATTR)"),
              tags$li(tags$b("ODE Compartments:"), " 25 (9 PK + 16 PD)"),
              tags$li(tags$b("Drug classes:"), " TTR stabilizer, siRNA (LNP/GalNAc), 2'-MOE ASO"),
              tags$li(tags$b("Scenarios:"), " Natural history (ATTRwt/ATTRv), Tafamidis, Patisiran, Vutrisiran, Inotersen"),
              tags$li(tags$b("Clinical calibration:"), " ATTR-ACT, APOLLO, HELIOS-A, NEURO-TTR, ATTRiBUTE-CM")
            )
          ),
          box(title="Key References", width=6, status="info",
            tags$ul(
              tags$li("Maurer et al. (2018). Tafamidis in ATTRwt-CM. NEJM 379:1007-1016 (ATTR-ACT)"),
              tags$li("Adams et al. (2018). Patisiran in hATTR. NEJM 379:11-21 (APOLLO)"),
              tags$li("Gillmore et al. (2021). Vutrisiran in ATTRv. NEJM 385:1407-1416 (HELIOS-A)"),
              tags$li("Benson et al. (2018). Inotersen in hATTR. Lancet 391:1864-1873 (NEURO-TTR)"),
              tags$li("Fontana et al. (2015). Mechanisms of ATTR Amyloidosis. NEJM")
            )
          )
        )
      )
    )
  )
)

## ─────────────────────────────────────────────────────────────────────────────
## Server
## ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## ── Simulation reactive ──────────────────────────────────────────────────
  sim_data <- eventReactive(input$run_btn, {
    req(input$drug1, input$attr_type, input$duration_yr)
    withProgress(message="Running ATTR simulation...", {
      ev_obj <- build_ev(input$drug1, 0, as.numeric(input$attr_type), input$duration_yr)
      run_sim(ev_obj, as.numeric(input$attr_type), input$duration_yr)
    })
  }, ignoreNULL = FALSE)

  ## ── Patient profile text outputs ─────────────────────────────────────────
  output$pt_phenotype  <- renderText(ifelse(input$attr_type=="1","ATTRwt (Cardiac)","ATTRv (Neuropathic)"))
  output$pt_age_disp   <- renderText(input$pt_age)
  output$pt_wt_disp    <- renderText(input$pt_wt)
  output$pt_drug_disp  <- renderText(input$drug1)
  output$pt_dur_disp   <- renderText(input$duration_yr)
  output$pt_ntp        <- renderText({ d <- sim_data(); nrow(d) })

  ## ── Value boxes ──────────────────────────────────────────────────────────
  last_row <- reactive({
    d <- sim_data()
    tail(d, 1)
  })

  output$vbox_disease <- renderValueBox({
    valueBox(ifelse(input$attr_type=="1","ATTRwt","ATTRv"),"ATTR Phenotype",
             icon=icon("dna"), color="blue")
  })
  output$vbox_drug <- renderValueBox({
    valueBox(input$drug1,"Treatment", icon=icon("pills"), color="green")
  })
  output$vbox_lvef <- renderValueBox({
    lr <- last_row()
    valueBox(sprintf("%.1f%%", lr$LVEF_pct), "Final LVEF",
             icon=icon("heartbeat"), color=ifelse(lr$LVEF_pct>=50,"green","red"))
  })
  output$vbox_nis <- renderValueBox({
    lr <- last_row()
    valueBox(sprintf("%.0f", lr$NIS_score), "Final NIS Score",
             icon=icon("brain"), color=ifelse(lr$NIS_score<40,"green","orange"))
  })
  output$vbox_lvef2 <- renderValueBox({
    lr <- last_row()
    valueBox(sprintf("%.1f%%", lr$LVEF_pct),"Final LVEF",
             icon=icon("heartbeat"),color=ifelse(lr$LVEF_pct>=50,"green","red"))
  })
  output$vbox_bnp <- renderValueBox({
    lr <- last_row()
    valueBox(sprintf("%.0f pg/mL", lr$NTproBNP),"NT-proBNP",
             icon=icon("chart-line"),color=ifelse(lr$NTproBNP<900,"orange","red"))
  })
  output$vbox_inflam <- renderValueBox({
    lr <- last_row()
    valueBox(sprintf("%.3f", lr$INFLAM_lvl),"Cardiac Inflammation",
             icon=icon("fire"),color="orange")
  })
  output$vbox_fibhrt <- renderValueBox({
    lr <- last_row()
    valueBox(sprintf("%.2f AU", lr$FIB_HRT_AU),"Cardiac Fibril Load",
             icon=icon("layer-group"),color="red")
  })
  output$vbox_nis2 <- renderValueBox({
    lr <- last_row()
    valueBox(sprintf("%.1f", lr$NIS_score),"NIS Score",
             icon=icon("brain"),color=ifelse(lr$NIS_score<40,"green","orange"))
  })
  output$vbox_mbmi <- renderValueBox({
    lr <- last_row()
    valueBox(sprintf("%.0f g/m²", lr$mBMI_val),"mBMI",
             icon=icon("weight"),color=ifelse(lr$mBMI_val>820,"green","red"))
  })
  output$vbox_egfr <- renderValueBox({
    lr <- last_row()
    valueBox(sprintf("%.0f", lr$eGFR_val),"eGFR (mL/min/1.73m²)",
             icon=icon("kidneys"),color=ifelse(lr$eGFR_val>=60,"green","orange"))
  })
  output$vbox_fibnrv <- renderValueBox({
    lr <- last_row()
    valueBox(sprintf("%.2f AU", lr$FIB_NRV_AU),"Nerve Fibril Load",
             icon=icon("bolt"),color="purple")
  })

  ## ── Helper for plotly line ────────────────────────────────────────────────
  make_plotly <- function(d, x, y, ylab, color="#2980b9", title="") {
    plot_ly(d, x=~time_yr, y=as.formula(paste0("~",y)),
            type="scatter", mode="lines",
            line=list(color=color, width=2)) %>%
      layout(xaxis=list(title="Time (years)"),
             yaxis=list(title=ylab), title=title,
             hovermode="closest")
  }

  ## ── Tab 1: Profile plot ──────────────────────────────────────────────────
  output$profile_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr) %>%
      add_lines(y=~LVEF_pct, name="LVEF (%)", line=list(color="#e74c3c")) %>%
      add_lines(y=~NIS_score, name="NIS", line=list(color="#3498db")) %>%
      add_lines(y=~TTR_mRNA_rel*100, name="TTR mRNA%", line=list(color="#2ecc71")) %>%
      layout(xaxis=list(title="Time (years)"), yaxis=list(title="Value"),
             legend=list(orientation="h"), title="Key Biomarkers Overview")
  })

  ## ── Tab 2: PK plots ──────────────────────────────────────────────────────
  output$pk_taf_plot <- renderPlotly({
    d <- sim_data() %>% filter(time <= 90*24)
    make_plotly(d,"time_yr","C_TAF_ug","Tafamidis (μg/mL)","#2980b9","Tafamidis PK")
  })
  output$pk_rna_plot <- renderPlotly({
    d <- sim_data() %>% filter(time <= 90*24)
    plot_ly(d, x=~time_yr) %>%
      add_lines(y=~C_VUT_ng, name="Vutrisiran (ng/mL)", line=list(color="#9b59b6")) %>%
      add_lines(y=~C_INO_ng, name="Inotersen (ng/mL)", line=list(color="#f39c12")) %>%
      add_lines(y=~C_PAT_ug*1000, name="Patisiran (ng/mL)", line=list(color="#27ae60")) %>%
      layout(xaxis=list(title="Time (years)"),yaxis=list(title="Conc (ng/mL)"),
             legend=list(orientation="h"))
  })
  output$effect_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr) %>%
      add_lines(y=~E_stab_pct, name="TTR Stabilization%", line=list(color="#2980b9")) %>%
      add_lines(y=~E_RNA_pct,  name="mRNA Knockdown%",    line=list(color="#27ae60")) %>%
      layout(xaxis=list(title="Time (years)"),yaxis=list(title="Effect (%)"),
             legend=list(orientation="h"))
  })
  output$tet_plot <- renderPlotly({
    make_plotly(sim_data(),"time_yr","TTR_TET_AU","TTR Tetramer (AU)","#3498db","TTR Tetramer Level")
  })

  ## ── Tab 3: TTR misfolding ────────────────────────────────────────────────
  output$mrna_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr, y=~TTR_mRNA_rel*100, type="scatter", mode="lines",
            line=list(color="#27ae60", width=2)) %>%
      layout(xaxis=list(title="Time (years)"),yaxis=list(title="TTR mRNA (% baseline)"),
             shapes=list(list(type="line",x0=0,x1=max(d$time_yr),
                              y0=20,y1=20,line=list(dash="dot",color="gray"))))
  })
  output$tet_mono_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr) %>%
      add_lines(y=~TTR_TET_AU, name="Tetramer", line=list(color="#3498db")) %>%
      add_lines(y=~FIB_HRT_AU, name="Cardiac fibrils", line=list(color="#e74c3c")) %>%
      layout(xaxis=list(title="Time (years)"),yaxis=list(title="AU"),
             legend=list(orientation="h"))
  })
  output$fibril_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr) %>%
      add_lines(y=~TTR_OLIGO_AU, name="Oligomers", line=list(color="#9b59b6")) %>%
      add_lines(y=~FIB_HRT_AU,   name="Cardiac Fibrils", line=list(color="#e74c3c")) %>%
      add_lines(y=~FIB_NRV_AU,   name="Nerve Fibrils", line=list(color="#2980b9")) %>%
      layout(xaxis=list(title="Time (years)"),yaxis=list(title="AU"),
             legend=list(orientation="h"))
  })
  output$fibril_dist_plot <- renderPlotly({
    d <- sim_data()
    last_d <- tail(d, 1)
    vals <- c(last_d$FIB_HRT_AU, last_d$FIB_NRV_AU,
              last_d$FIB_HRT_AU + last_d$FIB_NRV_AU * 0.14)
    plot_ly(labels=c("Cardiac","Peripheral Nerve","Systemic"),
            values=vals, type="pie",
            marker=list(colors=c("#e74c3c","#3498db","#2ecc71"))) %>%
      layout(title="Fibril Distribution at End of Simulation")
  })

  ## ── Tab 4: Cardiac ────────────────────────────────────────────────────────
  output$lvef_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr, y=~LVEF_pct, type="scatter", mode="lines",
            line=list(color="#e74c3c", width=2)) %>%
      layout(xaxis=list(title="Time (years)"),yaxis=list(title="LVEF (%)"),
             shapes=list(list(type="line",x0=0,x1=max(d$time_yr),y0=50,y1=50,
                              line=list(dash="dot",color="orange"))))
  })
  output$bnp_plot <- renderPlotly({
    make_plotly(sim_data(),"time_yr","NTproBNP","NT-proBNP (pg/mL)","#f39c12","NT-proBNP")
  })
  output$card_fibril_inflam_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr) %>%
      add_lines(y=~FIB_HRT_AU, name="Cardiac Fibrils (AU)", line=list(color="#e74c3c")) %>%
      add_lines(y=~INFLAM_lvl, name="Inflammation Index", line=list(color="#f39c12")) %>%
      layout(xaxis=list(title="Time (years)"),yaxis=list(title="AU / Index"),
             legend=list(orientation="h"))
  })
  output$nyha_plot <- renderPlotly({
    d <- sim_data() %>%
      mutate(NYHA_est = 1 + pmin(3, pmax(0, (62 - LVEF_pct)/10)))
    make_plotly(d,"time_yr","NYHA_est","Estimated NYHA Class","#c0392b","NYHA Class Estimate")
  })

  ## ── Tab 5: Neurological ──────────────────────────────────────────────────
  output$nis_plot <- renderPlotly({
    make_plotly(sim_data(),"time_yr","NIS_score","NIS Score","#3498db","Neuropathy Impairment Score")
  })
  output$mbmi_plot <- renderPlotly({
    make_plotly(sim_data(),"time_yr","mBMI_val","mBMI (g/m²)","#27ae60","Modified BMI")
  })
  output$nrv_fibril_plot <- renderPlotly({
    make_plotly(sim_data(),"time_yr","FIB_NRV_AU","Nerve Fibril Load (AU)","#9b59b6","Nerve Fibrils")
  })
  output$egfr_plot <- renderPlotly({
    make_plotly(sim_data(),"time_yr","eGFR_val","eGFR (mL/min/1.73m²)","#1abc9c","Renal Function")
  })

  ## ── Tab 6: Scenario comparison ───────────────────────────────────────────
  comp_data <- reactive({
    req(input$compare_drugs)
    attr_t <- as.numeric(input$attr_type)
    dur    <- input$duration_yr
    drugs  <- input$compare_drugs
    map_drug <- function(d) ifelse(d == "None (Natural History)", "None", d)
    withProgress(message="Running comparisons...", value=0, {
      bind_rows(lapply(seq_along(drugs), function(i) {
        incProgress(1/length(drugs))
        ev_o <- build_ev(map_drug(drugs[i]), 0, attr_t, dur)
        run_sim(ev_o, attr_t, dur) %>% mutate(Scenario=drugs[i])
      }))
    })
  })

  comp_colors <- c(
    "None (Natural History)"   = "#e74c3c",
    "Tafamidis 61mg QD"        = "#2980b9",
    "Vutrisiran 25mg Q3M"      = "#9b59b6",
    "Inotersen 300mg QW"       = "#f39c12",
    "Patisiran 0.3mg/kg Q3W"   = "#27ae60"
  )

  output$comp_lvef <- renderPlotly({
    d <- comp_data()
    plot_ly(d, x=~time_yr, y=~LVEF_pct, color=~Scenario, type="scatter", mode="lines",
            colors=comp_colors[unique(d$Scenario)]) %>%
      layout(xaxis=list(title="Years"),yaxis=list(title="LVEF (%)"),
             legend=list(orientation="h"))
  })
  output$comp_nis <- renderPlotly({
    d <- comp_data()
    plot_ly(d, x=~time_yr, y=~NIS_score, color=~Scenario, type="scatter", mode="lines",
            colors=comp_colors[unique(d$Scenario)]) %>%
      layout(xaxis=list(title="Years"),yaxis=list(title="NIS Score"),
             legend=list(orientation="h"))
  })
  output$comp_mrna <- renderPlotly({
    d <- comp_data()
    plot_ly(d, x=~time_yr, y=~TTR_mRNA_rel*100, color=~Scenario,
            type="scatter", mode="lines",
            colors=comp_colors[unique(d$Scenario)]) %>%
      layout(xaxis=list(title="Years"),yaxis=list(title="TTR mRNA (% baseline)"),
             legend=list(orientation="h"))
  })
  output$compare_table <- renderDT({
    d <- comp_data()
    target_time <- 18 * 730.5
    d %>% group_by(Scenario) %>%
      filter(abs(time - target_time) == min(abs(time - target_time))) %>%
      summarise(
        `LVEF (%)` = round(mean(LVEF_pct), 1),
        `NT-proBNP` = round(mean(NTproBNP)),
        `NIS Score` = round(mean(NIS_score), 1),
        `mBMI` = round(mean(mBMI_val)),
        `TTR mRNA (%)` = round(mean(TTR_mRNA_rel*100), 1),
        `eGFR` = round(mean(eGFR_val), 1),
        .groups="drop"
      ) %>%
      datatable(options=list(dom='t', pageLength=10), rownames=FALSE)
  })

  ## ── Tab 7: Biomarker dashboard ───────────────────────────────────────────
  output$bm_cardiac <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr) %>%
      add_lines(y=~LVEF_pct,   name="LVEF (%)",        yaxis="y",  line=list(color="#e74c3c")) %>%
      add_lines(y=~NTproBNP/100, name="NT-proBNP/100",  yaxis="y2", line=list(color="#f39c12")) %>%
      add_lines(y=~INFLAM_lvl*100, name="Inflam×100",   yaxis="y",  line=list(color="#c0392b",dash="dot")) %>%
      layout(
        xaxis = list(title="Years"),
        yaxis = list(title="LVEF (%) / Inflam×100"),
        yaxis2 = list(title="NT-proBNP/100 (pg/mL)", overlaying="y", side="right"),
        legend = list(orientation="h")
      )
  })
  output$bm_neuro <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr) %>%
      add_lines(y=~NIS_score, name="NIS Score",  line=list(color="#3498db")) %>%
      add_lines(y=~mBMI_val/10, name="mBMI/10",  line=list(color="#27ae60")) %>%
      add_lines(y=~eGFR_val, name="eGFR",         line=list(color="#1abc9c")) %>%
      layout(xaxis=list(title="Years"),yaxis=list(title="Score / Value"),
             legend=list(orientation="h"))
  })
  output$bm_summary_table <- renderDT({
    d <- sim_data()
    timepoints <- c(0, 365*1, 365*2, 365*3)
    bind_rows(lapply(timepoints, function(tp) {
      r <- d %>% filter(abs(time - tp*24) == min(abs(time - tp*24))) %>% head(1)
      tibble(
        `Time (yr)` = tp/365,
        `LVEF (%)` = round(r$LVEF_pct, 1),
        `NT-proBNP` = round(r$NTproBNP),
        `FIB_HRT (AU)` = round(r$FIB_HRT_AU, 3),
        `NIS Score` = round(r$NIS_score, 1),
        `mBMI (g/m²)` = round(r$mBMI_val),
        `eGFR` = round(r$eGFR_val, 1),
        `TTR mRNA (%)` = round(r$TTR_mRNA_rel*100, 1)
      )
    })) %>%
      datatable(options=list(dom='t'), rownames=FALSE)
  })
}

## ─── Launch ──────────────────────────────────────────────────────────────────
shinyApp(ui, server)
