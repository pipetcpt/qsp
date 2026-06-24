## ============================================================
## Acromegaly QSP Interactive Shiny Dashboard
## GH-Secreting Pituitary Adenoma — PK/PD Simulation
##
## Tabs:
##  1. Patient Profile & Disease Characteristics
##  2. Drug PK — Plasma Concentration Profiles
##  3. GH & IGF-1 (Primary PD Endpoints)
##  4. Cardiovascular & Metabolic Outcomes
##  5. Treatment Scenario Comparison
##  6. Biomarker Dashboard & Target Achievement
##  7. Tumor Volume & Disease Progression
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(tidyr)
library(DT)

## ---- Inline mrgsolve model (minimal for Shiny) ----
mod_code_shiny <- '
$PARAM
KA_SSA=0.0060 V1_SSA=14 V2_SSA=55 Q_SSA=4.5 CL_SSA=9.0
KA_PEG=0.040  V1_PEG=7  V2_PEG=45 Q_PEG=2.5 CL_PEG=0.45
IC50_SSTR2=0.4 IMAX_SSA=0.90 IC50_PEG_GHR=5.0 IMAX_PEG=0.97
KOUT_GH=0.693 GH_BASE=15 FB_IGF1=0.001
KACT_STAT=0.50 KDEG_STAT=0.35 SOCS2_FB=0.40
KPROD_IGF1=1.2 KCLEAR_IGF1=0.010
KG_ADENOM=0.00005 KK_ADENOM=5.0
KIN_LVH=0.002 KOUT_LVH=0.0004
KINS_GH=0.08 KOUT_GLUC=0.003
KARTH=0.000010 KARTH_REM=0.000001
DRUG_TYPE=1 USE_CAB=0 PAS_GLUC=0

$CMT DEPOT_SSA CENT_SSA PERI_SSA CENT_PEG PERI_PEG
     GH_PLASMA STAT5b_ACT IGF1_PLASMA
     ADENOM_VOL LVH_IDX GLUCOSE ARTH_SCORE

$INIT
GH_PLASMA=20 STAT5b_ACT=1.0 IGF1_PLASMA=450
ADENOM_VOL=1.5 LVH_IDX=145 GLUCOSE=105 ARTH_SCORE=10

$ODE
dxdt_DEPOT_SSA = -KA_SSA*DEPOT_SSA;
dxdt_CENT_SSA  =  KA_SSA*DEPOT_SSA
                 -(CL_SSA/V1_SSA)*CENT_SSA
                 -(Q_SSA/V1_SSA)*CENT_SSA
                 +(Q_SSA/V2_SSA)*PERI_SSA;
dxdt_PERI_SSA  = (Q_SSA/V1_SSA)*CENT_SSA-(Q_SSA/V2_SSA)*PERI_SSA;
dxdt_CENT_PEG  =  KA_PEG*PERI_PEG
                 -(CL_PEG/V1_PEG)*CENT_PEG
                 -(Q_PEG/V1_PEG)*CENT_PEG
                 +(Q_PEG/V2_PEG)*PERI_PEG;
dxdt_PERI_PEG  = -(CL_PEG/V2_PEG)*PERI_PEG
                  -(Q_PEG/V2_PEG)*PERI_PEG
                  +(Q_PEG/V1_PEG)*CENT_PEG;
double SSTR2_occ = CENT_SSA/(CENT_SSA+IC50_SSTR2+1e-9);
double SSA_eff   = IMAX_SSA*SSTR2_occ;
double PEG_eff   = IMAX_PEG*CENT_PEG/(CENT_PEG+IC50_PEG_GHR+1e-9);
double GH_prod   = GH_BASE*(1-SSA_eff)/(1+FB_IGF1*IGF1_PLASMA);
double cab_red   = USE_CAB*0.20;
GH_prod = GH_prod*(1-cab_red)*(ADENOM_VOL/1.5);
double PAS_add   = (DRUG_TYPE==3)?PAS_GLUC:0;
dxdt_GH_PLASMA  = GH_prod - KOUT_GH*GH_PLASMA;
double GHR_avail = (1-PEG_eff);
dxdt_STAT5b_ACT = KACT_STAT*GH_PLASMA*GHR_avail
                  - KDEG_STAT*(1+SOCS2_FB*STAT5b_ACT)*STAT5b_ACT;
double igf_prod  = KPROD_IGF1*STAT5b_ACT;
dxdt_IGF1_PLASMA = 50*igf_prod - KCLEAR_IGF1*IGF1_PLASMA;
dxdt_ADENOM_VOL = KG_ADENOM*ADENOM_VOL*(1-ADENOM_VOL/KK_ADENOM)*(1-0.10*SSA_eff);
dxdt_LVH_IDX    = KIN_LVH*GH_PLASMA - KOUT_LVH*(LVH_IDX-90);
dxdt_GLUCOSE    = KINS_GH*GH_PLASMA + PAS_add*0.001 - KOUT_GLUC*(GLUCOSE-85);
dxdt_ARTH_SCORE = KARTH*IGF1_PLASMA - KARTH_REM*ARTH_SCORE;

$TABLE
double SSA_conc = CENT_SSA;
double PEG_conc = CENT_PEG;
double GH_val   = GH_PLASMA;
double IGF1_val = IGF1_PLASMA;
double LVmass   = LVH_IDX;
double Gluc     = GLUCOSE;
double Arth     = ARTH_SCORE;
double Tumor    = ADENOM_VOL;
double SSTR_occ_out = IMAX_SSA*CENT_SSA/(CENT_SSA+IC50_SSTR2+1e-9);
double PEG_occ_out  = IMAX_PEG*CENT_PEG/(CENT_PEG+IC50_PEG_GHR+1e-9);
double GH_ctrl  = (GH_PLASMA < 2.5) ? 1.0 : 0.0;
double IGF1_ctrl= (IGF1_PLASMA < 250) ? 1.0 : 0.0;
double Remiss   = GH_ctrl*IGF1_ctrl;

$CAPTURE SSA_conc PEG_conc GH_val IGF1_val LVmass Gluc Arth Tumor
         SSTR_occ_out PEG_occ_out GH_ctrl IGF1_ctrl Remiss
'

mod_shiny <- mcode_cache("acromegaly_shiny", mod_code_shiny)

## ---- Color palette ----
pal_scenario <- c(
  "Untreated"           = "#e74c3c",
  "Surgery"             = "#27ae60",
  "Octreotide LAR"      = "#3498db",
  "Lanreotide AG"       = "#9b59b6",
  "Pasireotide LAR"     = "#e67e22",
  "Pegvisomant"         = "#1abc9c",
  "Combo SSA+Peg"       = "#f39c12",
  "Surgery+SSA rescue"  = "#2ecc71"
)

theme_acro_shiny <- theme_bw(base_size=13) +
  theme(
    legend.position="bottom",
    panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold", size=14),
    strip.background=element_rect(fill="#2c3e50"),
    strip.text=element_text(color="white")
  )

## ====================================================================
## UI
## ====================================================================
ui <- dashboardPage(
  skin="blue",
  dashboardHeader(title="Acromegaly QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",       tabName="tab_patient",   icon=icon("user")),
      menuItem("Drug PK",               tabName="tab_pk",        icon=icon("pills")),
      menuItem("GH & IGF-1 (PD)",       tabName="tab_pd",        icon=icon("chart-line")),
      menuItem("CV & Metabolic",        tabName="tab_cv",        icon=icon("heartbeat")),
      menuItem("Scenario Comparison",   tabName="tab_compare",   icon=icon("table")),
      menuItem("Biomarker Dashboard",   tabName="tab_biomarker", icon=icon("vials")),
      menuItem("Tumor Progression",     tabName="tab_tumor",     icon=icon("microscope"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .info-box { min-height: 70px; }
      .info-box-icon { height: 70px; line-height: 70px; }
      .info-box-content { padding-top: 10px; }
    "))),

    tabItems(

      ## ================================================================
      ## TAB 1: Patient Profile
      ## ================================================================
      tabItem(tabName="tab_patient",
        fluidRow(
          box(title="Patient Characteristics", status="primary", solidHeader=TRUE,
              width=4,
              numericInput("pt_age",   "Age (years)",             value=45, min=18, max=80),
              selectInput("pt_sex",    "Sex", choices=c("Male","Female")),
              numericInput("pt_height","Height (cm)",             value=170, min=140, max=200),
              numericInput("pt_weight","Weight (kg)",             value=80, min=40, max=150)
          ),
          box(title="Disease Parameters", status="warning", solidHeader=TRUE,
              width=4,
              selectInput("adenom_type", "Adenoma Type",
                          choices=c("Microadenoma (<10mm)","Macroadenoma (≥10mm)",
                                    "Giant adenoma (≥40mm)")),
              numericInput("adenom_vol", "Adenoma Volume (cm³)", value=1.5, min=0.1, max=10, step=0.1),
              numericInput("gh_base",    "Baseline GH (ng/mL)",  value=20, min=1, max=200),
              numericInput("igf1_base",  "Baseline IGF-1 (ng/mL)", value=450, min=50, max=2000),
              selectInput("gnas_mut", "GNAS Mutation",
                          choices=c("Unknown","Present (~40%)","Absent")),
              selectInput("sstr2_expr", "SSTR2 Expression (IHC)",
                          choices=c("High","Moderate","Low","Unknown"))
          ),
          box(title="Treatment Selection", status="success", solidHeader=TRUE,
              width=4,
              selectInput("drug_sel", "Primary Treatment",
                          choices=c("Untreated","Surgery","Octreotide LAR",
                                    "Lanreotide AG","Pasireotide LAR",
                                    "Pegvisomant","Combo SSA+Peg","Surgery+SSA rescue")),
              conditionalPanel(
                condition="input.drug_sel == 'Octreotide LAR'",
                sliderInput("oct_dose", "Dose (mg)", min=10, max=30, value=20, step=10)
              ),
              conditionalPanel(
                condition="input.drug_sel == 'Lanreotide AG'",
                sliderInput("lan_dose", "Dose (mg)", min=60, max=120, value=90, step=30)
              ),
              conditionalPanel(
                condition="input.drug_sel == 'Pasireotide LAR'",
                sliderInput("pas_dose", "Dose (mg)", min=40, max=60, value=40, step=20)
              ),
              conditionalPanel(
                condition="input.drug_sel == 'Pegvisomant'",
                sliderInput("peg_dose", "Daily Dose (mg)", min=10, max=30, value=15, step=5)
              ),
              checkboxInput("use_cab", "Add Cabergoline", FALSE),
              numericInput("sim_months", "Simulation Duration (months)", value=12, min=3, max=36),
              actionButton("run_sim", "Run Simulation", class="btn-primary btn-lg")
          )
        ),
        fluidRow(
          infoBoxOutput("ib_gh",   width=3),
          infoBoxOutput("ib_igf1", width=3),
          infoBoxOutput("ib_lvmi", width=3),
          infoBoxOutput("ib_gluc", width=3)
        ),
        fluidRow(
          box(title="Disease Overview", status="info", solidHeader=TRUE, width=12,
              uiOutput("disease_overview"))
        )
      ),

      ## ================================================================
      ## TAB 2: Drug PK
      ## ================================================================
      tabItem(tabName="tab_pk",
        fluidRow(
          box(title="PK Parameter Inputs", status="primary", solidHeader=TRUE, width=3,
              helpText("SSA = Somatostatin Analog; Peg = Pegvisomant"),
              numericInput("v1_ssa",  "V1_SSA (L)",  value=14),
              numericInput("v2_ssa",  "V2_SSA (L)",  value=55),
              numericInput("cl_ssa",  "CL_SSA (L/h)",value=9.0),
              numericInput("q_ssa",   "Q_SSA (L/h)", value=4.5),
              numericInput("ka_ssa",  "Ka_SSA (h⁻¹)",value=0.006),
              hr(),
              numericInput("v1_peg",  "V1_Peg (L)",  value=7),
              numericInput("cl_peg",  "CL_Peg (L/h)",value=0.45),
              numericInput("ka_peg",  "Ka_Peg (h⁻¹)",value=0.040)
          ),
          box(title="Drug Concentration–Time Profiles", status="success", solidHeader=TRUE,
              width=9,
              plotlyOutput("pk_plot", height="450px")
          )
        ),
        fluidRow(
          box(title="PK Parameters Summary", status="info", solidHeader=TRUE, width=12,
              DTOutput("pk_table"))
        )
      ),

      ## ================================================================
      ## TAB 3: GH & IGF-1 PD
      ## ================================================================
      tabItem(tabName="tab_pd",
        fluidRow(
          box(title="PD Control Thresholds", status="primary", solidHeader=TRUE, width=3,
              numericInput("gh_ctrl_thresh",   "GH Control Threshold (ng/mL)", value=2.5),
              numericInput("igf1_ctrl_thresh", "IGF-1 Normal Upper Limit (ng/mL)", value=250),
              numericInput("gh_cure_ogtt",     "GH Cure (OGTT nadir, ng/mL)", value=1.0),
              helpText("IGF-1 normalization = primary marker of disease control."),
              helpText("GH nadir after 75g OGTT < 1 ng/mL = biochemical cure.")
          ),
          box(title="GH & IGF-1 Response Over Time", status="success", solidHeader=TRUE,
              width=9,
              plotlyOutput("pd_plot", height="450px"))
        ),
        fluidRow(
          box(title="SSTR2 & GHR Receptor Occupancy", status="warning", solidHeader=TRUE,
              width=6,
              plotlyOutput("receptor_plot", height="350px")),
          box(title="Monthly Biochemical Control (%)", status="info", solidHeader=TRUE,
              width=6,
              plotlyOutput("remission_plot", height="350px"))
        )
      ),

      ## ================================================================
      ## TAB 4: Cardiovascular & Metabolic
      ## ================================================================
      tabItem(tabName="tab_cv",
        fluidRow(
          box(title="CV & Metabolic Outcomes", status="danger", solidHeader=TRUE, width=12,
              plotlyOutput("cv_metab_plot", height="500px"))
        ),
        fluidRow(
          box(title="LV Mass Index", status="danger", solidHeader=TRUE, width=6,
              plotlyOutput("lvh_plot", height="300px")),
          box(title="Fasting Glucose & Pasireotide Hyperglycemia", status="warning",
              solidHeader=TRUE, width=6,
              plotlyOutput("glucose_plot", height="300px"))
        )
      ),

      ## ================================================================
      ## TAB 5: Scenario Comparison
      ## ================================================================
      tabItem(tabName="tab_compare",
        fluidRow(
          box(title="Multi-Scenario Comparison", status="primary", solidHeader=TRUE,
              width=12,
              checkboxGroupInput("scen_sel", "Select Scenarios:",
                choices=c("Untreated","Surgery","Octreotide LAR","Lanreotide AG",
                          "Pasireotide LAR","Pegvisomant","Combo SSA+Peg",
                          "Surgery+SSA rescue"),
                selected=c("Untreated","Octreotide LAR","Pegvisomant","Surgery"),
                inline=TRUE)
          )
        ),
        fluidRow(
          box(title="GH Response", width=6, plotlyOutput("comp_gh", height="300px")),
          box(title="IGF-1 Response", width=6, plotlyOutput("comp_igf1", height="300px"))
        ),
        fluidRow(
          box(title="12-Month Outcome Summary Table", status="info",
              solidHeader=TRUE, width=12,
              DTOutput("compare_table"))
        )
      ),

      ## ================================================================
      ## TAB 6: Biomarker Dashboard
      ## ================================================================
      tabItem(tabName="tab_biomarker",
        fluidRow(
          valueBoxOutput("vb_gh_ctrl",   width=3),
          valueBoxOutput("vb_igf1_ctrl", width=3),
          valueBoxOutput("vb_lvmi_ctrl", width=3),
          valueBoxOutput("vb_gluc_ctrl", width=3)
        ),
        fluidRow(
          box(title="Biomarker Trajectory Heatmap", status="primary",
              solidHeader=TRUE, width=8,
              plotlyOutput("biomarker_heat", height="400px")),
          box(title="Target Achievement Gauge", status="success",
              solidHeader=TRUE, width=4,
              plotlyOutput("gauge_plot", height="400px"))
        ),
        fluidRow(
          box(title="Arthropathy Progression", status="warning",
              solidHeader=TRUE, width=6,
              plotlyOutput("arth_plot", height="300px")),
          box(title="Biochemical Remission Over Time", status="success",
              solidHeader=TRUE, width=6,
              plotlyOutput("remission_time", height="300px"))
        )
      ),

      ## ================================================================
      ## TAB 7: Tumor Progression
      ## ================================================================
      tabItem(tabName="tab_tumor",
        fluidRow(
          box(title="Tumor Volume Dynamics", status="primary",
              solidHeader=TRUE, width=8,
              plotlyOutput("tumor_plot", height="450px")),
          box(title="Tumor Growth Parameters", status="warning",
              solidHeader=TRUE, width=4,
              numericInput("kg_rate", "Tumor Growth Rate (cm³/h × 10⁻⁵)",
                           value=5, min=1, max=20),
              numericInput("kk_max",  "Carrying Capacity (cm³)",
                           value=5.0, min=1, max=20),
              numericInput("ssa_antiprof", "SSA Antiproliferative Effect (%)",
                           value=15, min=0, max=40),
              helpText("SSA agents have modest antiproliferative effects (15-30% tumor shrinkage in responsive patients)."),
              hr(),
              helpText("Surgery targets: microadenoma ~85% cure; macroadenoma ~50%."),
              helpText("SRS reduces GH secretion over 3-10 years.")
          )
        ),
        fluidRow(
          box(title="MRI Tumor Size Simulation vs. Clinical Data",
              status="info", solidHeader=TRUE, width=12,
              plotlyOutput("tumor_compare", height="350px"))
        )
      )

    ) # end tabItems
  ) # end dashboardBody
)

## ====================================================================
## SERVER
## ====================================================================
server <- function(input, output, session) {

  ## ---- Reactive: run simulation ----
  sim_data <- eventReactive(input$run_sim, {
    req(input$drug_sel)

    t_end <- input$sim_months * 720  # hours
    tgrid <- seq(0, t_end, by=6)
    monthly <- seq(0, t_end, by=672)

    drug <- input$drug_sel
    gh0  <- input$gh_base
    igf0 <- input$igf1_base
    av0  <- input$adenom_vol

    # Build events based on drug selection
    build_ev <- function(drug, oct_d, lan_d, pas_d, peg_d) {
      if(drug == "Untreated") {
        return(ev(time=0, amt=0, cmt=1))
      } else if(drug == "Surgery") {
        return(ev(time=0, amt=0, cmt=1))
      } else if(drug == "Octreotide LAR") {
        evs <- lapply(monthly, function(t) ev(time=t, amt=oct_d, cmt=1))
        return(do.call(c, evs))
      } else if(drug == "Lanreotide AG") {
        evs <- lapply(monthly, function(t) ev(time=t, amt=lan_d, cmt=1))
        return(do.call(c, evs))
      } else if(drug == "Pasireotide LAR") {
        evs <- lapply(monthly, function(t) ev(time=t, amt=pas_d, cmt=1))
        return(do.call(c, evs))
      } else if(drug == "Pegvisomant") {
        daily <- seq(0, t_end, by=24)
        evs <- lapply(daily, function(t) ev(time=t, amt=peg_d, cmt=5))
        return(do.call(c, evs))
      } else if(drug == "Combo SSA+Peg") {
        daily <- seq(0, t_end, by=24)
        ev_ssa <- lapply(monthly, function(t) ev(time=t, amt=20, cmt=1))
        ev_peg <- lapply(daily,   function(t) ev(time=t, amt=10, cmt=5))
        return(do.call(c, c(ev_ssa, ev_peg)))
      } else if(drug == "Surgery+SSA rescue") {
        late_monthly <- monthly[monthly > 720]
        evs <- lapply(late_monthly, function(t) ev(time=t, amt=30, cmt=1))
        return(do.call(c, c(list(ev(time=0, amt=0, cmt=1)), evs)))
      }
    }

    events <- build_ev(drug, input$oct_dose, input$lan_dose,
                       input$pas_dose, input$peg_dose)

    # Set parameters
    surg_done <- ifelse(drug %in% c("Surgery","Surgery+SSA rescue"), 1, 0)
    surg_kill <- ifelse(drug == "Surgery", 0.70,
                        ifelse(drug == "Surgery+SSA rescue", 0.45, 0.0))
    dt_code <- switch(drug,
      "Octreotide LAR" = 1,
      "Lanreotide AG"  = 2,
      "Pasireotide LAR"= 3,
      "Pegvisomant"    = 4,
      "Combo SSA+Peg"  = 5,
      1)
    pas_gluc_add <- ifelse(drug == "Pasireotide LAR", 15, 0)
    gh_init <- ifelse(surg_done==1, gh0*(1-surg_kill), gh0)
    igf1_init<- ifelse(surg_done==1, igf0*0.8, igf0)
    av_init  <- ifelse(surg_done==1, av0*(1-surg_kill*0.7), av0)

    mod_shiny %>%
      param(DRUG_TYPE=dt_code, USE_CAB=as.numeric(input$use_cab),
            PAS_GLUC=pas_gluc_add,
            GH_BASE=gh0*max(1-surg_kill,0.3),
            KA_SSA=input$ka_ssa, V1_SSA=input$v1_ssa,
            V2_SSA=input$v2_ssa, CL_SSA=input$cl_ssa, Q_SSA=input$q_ssa,
            V1_PEG=input$v1_peg, CL_PEG=input$cl_peg, KA_PEG=input$ka_peg) %>%
      init(GH_PLASMA=gh_init, IGF1_PLASMA=igf1_init,
           ADENOM_VOL=av_init) %>%
      ev(events) %>%
      mrgsim(tgrid=tgrid) %>%
      as.data.frame() %>%
      mutate(time_mo = time/720, scenario=drug)
  }, ignoreNULL=FALSE)

  ## ---- Info boxes ----
  output$ib_gh <- renderInfoBox({
    d <- sim_data(); if(is.null(d)) return(infoBox("GH","—","red"))
    last <- tail(d, 1)
    ctrl <- last$GH_val < 2.5
    infoBox("GH (final)", paste0(round(last$GH_val,1)," ng/mL"),
            subtitle=ifelse(ctrl,"✓ Controlled","✗ Elevated"),
            icon=icon("chart-line"),
            color=ifelse(ctrl,"green","red"))
  })
  output$ib_igf1 <- renderInfoBox({
    d <- sim_data(); if(is.null(d)) return(infoBox("IGF-1","—","red"))
    last <- tail(d,1)
    ctrl <- last$IGF1_val < 250
    infoBox("IGF-1 (final)", paste0(round(last$IGF1_val,0)," ng/mL"),
            subtitle=ifelse(ctrl,"✓ Normalized","✗ Elevated"),
            icon=icon("vials"), color=ifelse(ctrl,"green","orange"))
  })
  output$ib_lvmi <- renderInfoBox({
    d <- sim_data(); if(is.null(d)) return(infoBox("LV Mass","—","yellow"))
    last <- tail(d,1)
    ctrl <- last$LVmass < 115
    infoBox("LV Mass Index", paste0(round(last$LVmass,0)," g/m²"),
            subtitle=ifelse(ctrl,"✓ Normal","✗ LVH"),
            icon=icon("heartbeat"), color=ifelse(ctrl,"green","yellow"))
  })
  output$ib_gluc <- renderInfoBox({
    d <- sim_data(); if(is.null(d)) return(infoBox("Glucose","—","yellow"))
    last <- tail(d,1)
    ctrl <- last$Gluc < 100
    infoBox("Fasting Glucose", paste0(round(last$Gluc,0)," mg/dL"),
            subtitle=ifelse(ctrl,"✓ Normal","✗ Elevated"),
            icon=icon("tint"), color=ifelse(ctrl,"green","yellow"))
  })

  ## ---- Disease overview text ----
  output$disease_overview <- renderUI({
    HTML(paste0(
      "<div style='padding:10px;'>",
      "<h4>Acromegaly — QSP Model Summary</h4>",
      "<p><b>Pathophysiology:</b> GH-secreting pituitary adenoma causes chronic GH/IGF-1 excess,
      leading to acral growth, arthropathy, cardiovascular complications (LVH, cardiomyopathy),
      metabolic disturbances (insulin resistance, T2DM), sleep apnea, and increased CV mortality (SMR 1.3-2.0).</p>",
      "<p><b>Treatment Goals:</b> Biochemical remission = random GH <2.5 ng/mL OR GH nadir (OGTT) <1 ng/mL
      PLUS age/sex-normalized IGF-1.</p>",
      "<p><b>Current Drug:</b> <b>", input$drug_sel, "</b><br>",
      "Simulation: <b>", input$sim_months, " months</b></p>",
      "<table border='1' cellpadding='5' style='border-collapse:collapse;width:100%;'>",
      "<tr style='background:#2c3e50;color:white;'>
        <th>Drug</th><th>IGF-1 Norm Rate</th><th>GH Control</th><th>Tumor Shrinkage</th><th>Key AE</th></tr>",
      "<tr><td>Octreotide LAR</td><td>25-35%</td><td>GH<2.5: 57%</td><td>20-30%</td><td>GI, cholelithiasis</td></tr>",
      "<tr><td>Lanreotide AG</td><td>27-38%</td><td>similar</td><td>25-30%</td><td>GI, cholelithiasis</td></tr>",
      "<tr><td>Pasireotide LAR</td><td>38-48%</td><td>higher</td><td>40-50%</td><td>Hyperglycemia 57%</td></tr>",
      "<tr><td>Pegvisomant</td><td>63-95%</td><td>via IGF-1</td><td>none (may grow)</td><td>LFT elevation 2-5%</td></tr>",
      "<tr><td>Surgery</td><td>85% (micro) / 50% (macro)</td><td>best</td><td>complete</td><td>hypopituitarism</td></tr>",
      "</table></div>"
    ))
  })

  ## ---- PK Plot ----
  output$pk_plot <- renderPlotly({
    d <- sim_data(); req(d)
    p <- d %>%
      pivot_longer(cols=c(SSA_conc, PEG_conc),
                   names_to="drug_comp", values_to="conc") %>%
      mutate(drug_comp = recode(drug_comp,
        SSA_conc="SSA Central (ng/mL)",
        PEG_conc="Pegvisomant Central (ng/mL)")) %>%
      filter(conc > 0.001) %>%
      ggplot(aes(x=time_mo, y=conc, color=drug_comp)) +
      geom_line(linewidth=1) +
      facet_wrap(~drug_comp, scales="free_y") +
      labs(title="Drug PK — Central Compartment Concentration",
           x="Time (months)", y="Concentration (ng/mL)") +
      theme_acro_shiny + theme(legend.position="none")
    ggplotly(p)
  })

  ## ---- PK table ----
  output$pk_table <- renderDT({
    tribble(
      ~Drug,             ~"t1/2 (h)", ~"Vd (L)",~"CL (L/h)",~"Bioavail.",~"Route",   ~"Frequency",
      "Octreotide LAR",  115,           70,       9.0,        "~100% (IM)",  "IM",      "q28d",
      "Lanreotide AG",   ~92,           62,       7.5,        "~70% (SC)",   "Deep SC", "q28d",
      "Pasireotide LAR", ~550,          200,      5.0,        "~100% (IM)",  "IM",      "q28d",
      "Pegvisomant",     138,           52,       0.45,       "~60-80% (SC)","SC",      "Daily"
    ) %>%
      datatable(options=list(pageLength=10, scrollX=TRUE), rownames=FALSE)
  })

  ## ---- GH & IGF-1 PD ----
  output$pd_plot <- renderPlotly({
    d <- sim_data(); req(d)
    p <- d %>%
      pivot_longer(cols=c(GH_val, IGF1_val),
                   names_to="marker", values_to="value") %>%
      mutate(marker=recode(marker, GH_val="GH (ng/mL)", IGF1_val="IGF-1 (ng/mL)")) %>%
      ggplot(aes(x=time_mo, y=value, color=marker)) +
      geom_line(linewidth=1.2) +
      geom_hline(data=data.frame(marker=c("GH (ng/mL)","IGF-1 (ng/mL)"),
                                 thresh=c(input$gh_ctrl_thresh, input$igf1_ctrl_thresh)),
                 aes(yintercept=thresh), linetype="dashed", color="red") +
      facet_wrap(~marker, scales="free_y") +
      labs(title="GH & IGF-1 PD Response", x="Time (months)", y="Level") +
      theme_acro_shiny + theme(legend.position="none")
    ggplotly(p)
  })

  output$receptor_plot <- renderPlotly({
    d <- sim_data(); req(d)
    p <- d %>%
      pivot_longer(cols=c(SSTR_occ_out, PEG_occ_out),
                   names_to="receptor", values_to="occupancy") %>%
      mutate(receptor=recode(receptor,
        SSTR_occ_out="SSTR2 Occupancy",
        PEG_occ_out="GHR Blockade (Pegvisomant)"),
        occupancy=occupancy*100) %>%
      ggplot(aes(x=time_mo, y=occupancy, color=receptor)) +
      geom_line(linewidth=1.1) +
      labs(title="Receptor Occupancy / GHR Blockade",
           x="Time (months)", y="Occupancy (%)", color="") +
      ylim(0,100) + theme_acro_shiny
    ggplotly(p)
  })

  output$remission_plot <- renderPlotly({
    d <- sim_data(); req(d)
    p <- d %>%
      group_by(time_mo_bin = floor(time_mo)) %>%
      summarise(remiss_pct = mean(Remiss, na.rm=TRUE)*100, .groups="drop") %>%
      ggplot(aes(x=time_mo_bin, y=remiss_pct)) +
      geom_col(fill="#27ae60", alpha=0.8) +
      labs(title="Monthly Biochemical Control (%)",
           x="Month", y="Control (%)") +
      ylim(0,100) + theme_acro_shiny
    ggplotly(p)
  })

  ## ---- CV & Metabolic ----
  output$cv_metab_plot <- renderPlotly({
    d <- sim_data(); req(d)
    p <- d %>%
      pivot_longer(cols=c(LVmass, Gluc, Arth),
                   names_to="param", values_to="value") %>%
      mutate(param=recode(param,
        LVmass="LV Mass Index (g/m²)",
        Gluc="Fasting Glucose (mg/dL)",
        Arth="Arthropathy Score (0-100)")) %>%
      ggplot(aes(x=time_mo, y=value, color=param)) +
      geom_line(linewidth=1.1) +
      facet_wrap(~param, scales="free_y") +
      labs(title="Cardiovascular & Metabolic Outcomes",
           x="Time (months)", y="Value") +
      theme_acro_shiny + theme(legend.position="none")
    ggplotly(p)
  })

  output$lvh_plot <- renderPlotly({
    d <- sim_data(); req(d)
    p <- d %>%
      ggplot(aes(x=time_mo, y=LVmass)) +
      geom_line(color="#c0392b", linewidth=1.2) +
      geom_hline(yintercept=115, linetype="dashed", color="red") +
      annotate("text", x=max(d$time_mo)*0.7, y=117, label="Normal limit (M=115)", size=3) +
      labs(title="LV Mass Index", x="Time (months)", y="LV Mass Index (g/m²)") +
      theme_acro_shiny
    ggplotly(p)
  })

  output$glucose_plot <- renderPlotly({
    d <- sim_data(); req(d)
    p <- d %>%
      ggplot(aes(x=time_mo, y=Gluc)) +
      geom_line(color="#e67e22", linewidth=1.2) +
      geom_hline(yintercept=100, linetype="dashed", color="blue") +
      geom_hline(yintercept=126, linetype="dashed", color="red") +
      annotate("text", x=max(d$time_mo)*0.6, y=102, label="Impaired fasting (100)", size=3) +
      annotate("text", x=max(d$time_mo)*0.6, y=128, label="Diabetes (126)", size=3, color="red") +
      labs(title="Fasting Glucose", x="Time (months)", y="Fasting Glucose (mg/dL)") +
      theme_acro_shiny
    ggplotly(p)
  })

  ## ---- Scenario Comparison ----
  multi_sim <- reactive({
    t_end  <- 8760
    tgrid  <- seq(0, t_end, by=24)
    monthly <- seq(0, t_end, by=672)
    daily   <- seq(0, t_end, by=24)

    scenarios <- list(
      "Untreated"    = list(ev=ev(time=0,amt=0,cmt=1), p=list(DRUG_TYPE=1)),
      "Surgery"      = list(ev=ev(time=0,amt=0,cmt=1),
                            p=list(DRUG_TYPE=1), init=list(GH_PLASMA=6,IGF1_PLASMA=300,ADENOM_VOL=0.45)),
      "Octreotide LAR" = list(
        ev=do.call(c, lapply(monthly, function(t) ev(time=t,amt=30,cmt=1))),
        p=list(DRUG_TYPE=1)),
      "Lanreotide AG" = list(
        ev=do.call(c, lapply(monthly, function(t) ev(time=t,amt=120,cmt=1))),
        p=list(DRUG_TYPE=2, KA_SSA=0.004, CL_SSA=7.5)),
      "Pasireotide LAR" = list(
        ev=do.call(c, lapply(monthly, function(t) ev(time=t,amt=60,cmt=1))),
        p=list(DRUG_TYPE=3, KA_SSA=0.0055, CL_SSA=5.0, V1_SSA=100, PAS_GLUC=15)),
      "Pegvisomant" = list(
        ev=do.call(c, lapply(daily, function(t) ev(time=t,amt=15,cmt=5))),
        p=list(DRUG_TYPE=4)),
      "Combo SSA+Peg" = list(
        ev=do.call(c, c(lapply(monthly, function(t) ev(time=t,amt=20,cmt=1)),
                        lapply(daily,   function(t) ev(time=t,amt=10,cmt=5)))),
        p=list(DRUG_TYPE=5)),
      "Surgery+SSA rescue" = list(
        ev=do.call(c, c(list(ev(time=0,amt=0,cmt=1)),
                        lapply(monthly[monthly>720], function(t) ev(time=t,amt=30,cmt=1)))),
        p=list(DRUG_TYPE=1),
        init=list(GH_PLASMA=12, IGF1_PLASMA=400, ADENOM_VOL=0.9))
    )

    sel <- input$scen_sel
    bind_rows(lapply(sel, function(nm) {
      sc  <- scenarios[[nm]]
      m   <- do.call(param, c(list(mod_shiny), sc$p))
      if(!is.null(sc$init)) m <- do.call(init, c(list(m), sc$init))
      m %>% ev(sc$ev) %>%
        mrgsim(tgrid=tgrid) %>%
        as.data.frame() %>%
        mutate(scenario=nm, time_mo=time/720)
    }))
  })

  output$comp_gh <- renderPlotly({
    d <- multi_sim(); req(d)
    p <- d %>%
      ggplot(aes(x=time_mo, y=GH_val, color=scenario)) +
      geom_line(linewidth=0.9) +
      geom_hline(yintercept=2.5, linetype="dashed", color="red") +
      labs(title="GH by Scenario", x="Month", y="GH (ng/mL)") +
      theme_acro_shiny
    ggplotly(p)
  })

  output$comp_igf1 <- renderPlotly({
    d <- multi_sim(); req(d)
    p <- d %>%
      ggplot(aes(x=time_mo, y=IGF1_val, color=scenario)) +
      geom_line(linewidth=0.9) +
      geom_hline(yintercept=250, linetype="dashed", color="red") +
      labs(title="IGF-1 by Scenario", x="Month", y="IGF-1 (ng/mL)") +
      theme_acro_shiny
    ggplotly(p)
  })

  output$compare_table <- renderDT({
    d <- multi_sim(); req(d)
    d %>%
      filter(abs(time_mo - round(time_mo)) < 0.1) %>%
      filter(time_mo %in% c(0,3,6,12)) %>%
      group_by(scenario, time_mo) %>%
      summarise(
        GH=round(mean(GH_val),1),
        IGF1=round(mean(IGF1_val),0),
        LVMI=round(mean(LVmass),0),
        Glucose=round(mean(Gluc),0),
        Tumor_cm3=round(mean(Tumor),2),
        Remission_pct=round(mean(Remiss)*100,0),
        .groups="drop") %>%
      rename(Month=time_mo, `GH (ng/mL)`=GH, `IGF-1 (ng/mL)`=IGF1,
             `LVMI (g/m²)`=LVMI, `Glucose (mg/dL)`=Glucose,
             `Tumor (cm³)`=Tumor_cm3, `Remission %`=Remission_pct) %>%
      datatable(options=list(pageLength=20, scrollX=TRUE),
                rownames=FALSE,
                filter="top")
  })

  ## ---- Biomarker Dashboard ----
  output$vb_gh_ctrl <- renderValueBox({
    d <- sim_data(); req(d)
    pct <- mean(d$GH_ctrl, na.rm=TRUE)*100
    valueBox(paste0(round(pct),"% of time"),
             "GH Controlled (<2.5 ng/mL)",
             icon=icon("check-circle"),
             color=ifelse(pct>50,"green","red"))
  })
  output$vb_igf1_ctrl <- renderValueBox({
    d <- sim_data(); req(d)
    pct <- mean(d$IGF1_ctrl, na.rm=TRUE)*100
    valueBox(paste0(round(pct),"% of time"),
             "IGF-1 Normalized",
             icon=icon("vials"),
             color=ifelse(pct>50,"green","orange"))
  })
  output$vb_lvmi_ctrl <- renderValueBox({
    d <- sim_data(); req(d)
    final_lvmi <- tail(d$LVmass,1)
    valueBox(paste0(round(final_lvmi)," g/m²"),
             "Final LV Mass Index",
             icon=icon("heartbeat"),
             color=ifelse(final_lvmi<115,"green","yellow"))
  })
  output$vb_gluc_ctrl <- renderValueBox({
    d <- sim_data(); req(d)
    final_gluc <- tail(d$Gluc,1)
    valueBox(paste0(round(final_gluc)," mg/dL"),
             "Final Fasting Glucose",
             icon=icon("tint"),
             color=ifelse(final_gluc<100,"green",
                          ifelse(final_gluc<126,"yellow","red")))
  })

  output$biomarker_heat <- renderPlotly({
    d <- sim_data(); req(d)
    d_norm <- d %>%
      filter(abs(time_mo - round(time_mo)) < 0.15) %>%
      group_by(time_mo=round(time_mo)) %>%
      summarise(
        GH_norm     = GH_val/20,
        IGF1_norm   = IGF1_val/450,
        LVMI_norm   = LVmass/145,
        Glucose_norm= Gluc/105,
        Arth_norm   = Arth/10,
        Tumor_norm  = Tumor/1.5,
        .groups="drop") %>%
      pivot_longer(cols=-time_mo, names_to="biomarker", values_to="norm_val") %>%
      mutate(biomarker=gsub("_norm","",biomarker))

    p <- d_norm %>%
      ggplot(aes(x=time_mo, y=biomarker, fill=norm_val)) +
      geom_tile(color="white") +
      scale_fill_gradient2(low="#27ae60", mid="#f39c12", high="#e74c3c",
                           midpoint=0.7, name="Normalized\nValue") +
      labs(title="Biomarker Heatmap (normalized to baseline)",
           x="Month", y="Biomarker") +
      theme_acro_shiny
    ggplotly(p)
  })

  output$gauge_plot <- renderPlotly({
    d <- sim_data(); req(d)
    final <- tail(d,1)
    targets_met <- sum(c(
      final$GH_val < 2.5,
      final$IGF1_val < 250,
      final$LVmass < 115,
      final$Gluc < 100,
      final$Remiss > 0.5
    ))
    plot_ly(type="indicator", mode="gauge+number",
            value=targets_met,
            title=list(text="Targets Met (out of 5)"),
            gauge=list(
              axis=list(range=list(0,5)),
              bar=list(color="#27ae60"),
              steps=list(
                list(range=c(0,2), color="#e74c3c"),
                list(range=c(2,4), color="#f39c12"),
                list(range=c(4,5), color="#27ae60")
              )
            ))
  })

  output$arth_plot <- renderPlotly({
    d <- sim_data(); req(d)
    p <- d %>%
      ggplot(aes(x=time_mo, y=Arth)) +
      geom_line(color="#784212", linewidth=1.2) +
      labs(title="Arthropathy Score",
           subtitle="IGF-1 driven joint damage (partially irreversible)",
           x="Time (months)", y="Score (0-100)") +
      theme_acro_shiny
    ggplotly(p)
  })

  output$remission_time <- renderPlotly({
    d <- sim_data(); req(d)
    p <- d %>%
      group_by(time_bin=floor(time_mo)) %>%
      summarise(remiss=mean(Remiss,na.rm=TRUE)*100,.groups="drop") %>%
      ggplot(aes(x=time_bin, y=remiss)) +
      geom_area(fill="#27ae60", alpha=0.5) +
      geom_line(color="#27ae60", linewidth=1) +
      labs(title="Biochemical Remission Over Time",
           x="Month", y="% in Remission") +
      ylim(0,100) + theme_acro_shiny
    ggplotly(p)
  })

  ## ---- Tumor Progression ----
  output$tumor_plot <- renderPlotly({
    d <- sim_data(); req(d)
    p <- d %>%
      ggplot(aes(x=time_mo, y=Tumor)) +
      geom_line(color="#8e44ad", linewidth=1.3) +
      geom_hline(yintercept=1.0, linetype="dashed", color="blue") +
      annotate("text", x=max(d$time_mo)*0.5, y=1.05,
               label="Microadenoma/Macroadenoma boundary (10mm)", size=3) +
      labs(title="Pituitary Adenoma Volume Over Time",
           x="Time (months)", y="Tumor Volume (cm³)") +
      theme_acro_shiny
    ggplotly(p)
  })

  output$tumor_compare <- renderPlotly({
    # Simulated vs approximate clinical SSA tumor shrinkage data
    clinical_ref <- tibble(
      month=c(3,6,12,24),
      Oct_LAR_shrink=c(10,18,22,25),
      Pas_LAR_shrink=c(15,28,35,42),
      label=c("3mo","6mo","12mo","24mo")
    ) %>% pivot_longer(cols=c(Oct_LAR_shrink, Pas_LAR_shrink),
                       names_to="drug", values_to="pct_shrink")

    p <- ggplot(clinical_ref, aes(x=month, y=pct_shrink, color=drug, linetype="Clinical")) +
      geom_point(size=3) + geom_line() +
      labs(title="Tumor Shrinkage: Simulation vs. Clinical Reference Data",
           subtitle="Reference: Colao et al. 2009, Gadelha et al. 2014",
           x="Month", y="Tumor Volume Reduction (%)", color="Drug") +
      scale_x_continuous(breaks=c(3,6,12,24)) +
      theme_acro_shiny
    ggplotly(p)
  })

}

## ---- Launch ----
shinyApp(ui, server)
