## ============================================================
## Multiple Myeloma QSP — Shiny Interactive Dashboard
## 6 Tabs: Patient Profile | PK | PD/Tumor Burden |
##         Clinical Endpoints | Scenario Comparison | Biomarkers
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)

## ---- Embedded mrgsolve model (minimal version for Shiny) ----
mm_code_shiny <- '
$PARAM
kgrow=0.015 kapop=0.005 Kmax=5e9 kres=1e-4 kgrow_res=0.018
kMP_prod=1.5e-9 kMP_elim=0.08 kFLC_prod=2e-10 kFLC_elim=0.15
kIL6_prod=0.05 kIL6_elim=1.2 kIL6_stim=0.3
k_OB_form=0.0016 k_OB_apop=0.002 k_OC_form=0.0006 k_OC_apop=0.0025
k_RANKL=0.005 k_OPG=0.002 k_OPG_inh=0.3
BoneV0=100 k_bone_res=0.0003 k_bone_form=0.0002
DKK1_base=1 kDKK1_MM=0.002
kHgb_recover=0.02 kHgb_supp=0.3 Hgb0=14 IL6_Hgb_inh=0.5
kB2M_prod=1e-10 kB2M_elim=0.25
BTZ_V1=4.7 BTZ_CL=9 BTZ_V2=34.5 BTZ_V3=298 BTZ_Q2=22 BTZ_Q3=5.6
BTZ_Emax=0.90 BTZ_EC50=2.5 BTZ_hill=2
LEN_Ka=1.2 LEN_CL=3.2 LEN_V=65 LEN_F=0.9
LEN_Emax=0.65 LEN_EC50=0.8 LEN_hill=1.5
DARA_CL=0.0029 DARA_V1=0.043 DARA_V2=0.028 DARA_Q=0.0019
DARA_kon=1.85 DARA_koff=0.37 DARA_kdeg=0.83 CD38_0=50
DARA_Emax=0.85 DARA_EC50=0.15 DARA_hill=2
DEX_Ka=2.5 DEX_CL=20 DEX_V=130 DEX_F=0.78
DEX_Emax=0.55 DEX_EC50=0.10 DEX_hill=1.5
VEN_Ka=0.8 VEN_CL=14 VEN_V=256 VEN_F=0.55
VEN_Emax=0.80 VEN_EC50=0.50 VEN_hill=1.8
ZOL_CL=3.8 ZOL_V1=4 ZOL_V2=22 ZOL_Q=0.7 ZOL_kbone=0.04 ZOL_krel=0.0001
ZOL_OC_Emax=0.75 ZOL_OC_EC50=0.05
$INIT
MM_S=1e7 MM_R=1e4 MP=3 FLC=150 IL6=0.1 VEGF=0.05
OB=100 OC=100 BV=100 NTX=1 PINP=1 Hgb=14 B2M=2.5
BTZ1=0 BTZ2=0 BTZ3=0 LEN_gut=0 LEN1=0
DARA1=0 DARA2=0 DARA_CD38=0
DEX_gut=0 DEX1=0 VEN_gut=0 VEN1=0
ZOL1=0 ZOL2=0 ZOL_bone=0
$ODE
double BTZ_Cp=BTZ1/BTZ_V1;
double LEN_Cp=LEN1/LEN_V;
double DARA_Cp=DARA1/(DARA_V1*70.0);
double DEX_Cp=DEX1/DEX_V;
double VEN_Cp=VEN1/VEN_V;
double ZOL_bone_Cp=ZOL_bone/ZOL_V2;
double E_BTZ=BTZ_Emax*pow(BTZ_Cp,BTZ_hill)/(pow(BTZ_EC50,BTZ_hill)+pow(BTZ_Cp,BTZ_hill));
double E_LEN=LEN_Emax*pow(LEN_Cp,LEN_hill)/(pow(LEN_EC50,LEN_hill)+pow(LEN_Cp,LEN_hill));
double E_DARA=DARA_Emax*pow(DARA_Cp,DARA_hill)/(pow(DARA_EC50,DARA_hill)+pow(DARA_Cp,DARA_hill));
double E_DEX=DEX_Emax*pow(DEX_Cp,DEX_hill)/(pow(DEX_EC50,DEX_hill)+pow(DEX_Cp,DEX_hill));
double E_VEN=VEN_Emax*pow(VEN_Cp,VEN_hill)/(pow(VEN_EC50,VEN_hill)+pow(VEN_Cp,VEN_hill));
double E_ZOL=ZOL_OC_Emax*ZOL_bone_Cp/(ZOL_OC_EC50+ZOL_bone_Cp);
double total_kill=E_BTZ+E_LEN+E_DARA+E_DEX+E_VEN;
if(total_kill>0.99) total_kill=0.99;
double N_total=MM_S+MM_R;
double IL6_effect=1.0+kIL6_stim*IL6;
double MM_burden=N_total/1e9;
dxdt_MM_S=kgrow*MM_S*(1.0-N_total/Kmax)*IL6_effect-kapop*MM_S-(kapop*total_kill/(1.0-total_kill))*MM_S-kres*MM_S;
dxdt_MM_R=kgrow_res*MM_R*(1.0-N_total/Kmax)-kapop*MM_R+kres*MM_S;
dxdt_MP=kMP_prod*MM_S-kMP_elim*MP;
dxdt_FLC=kFLC_prod*MM_S-kFLC_elim*FLC;
dxdt_IL6=kIL6_prod+0.5*MM_burden-kIL6_elim*IL6;
dxdt_VEGF=kVEGF_prod+0.3*MM_burden-kVEGF_elim*VEGF;
double RANKL_eff=1.0+k_RANKL*MM_burden*1e3;
double OPG_eff=k_OPG*OB/(1.0+k_OPG_inh*MM_burden);
double DKK1=DKK1_base+kDKK1_MM*MM_burden;
double Wnt_inh=1.0/(1.0+DKK1);
dxdt_OB=k_OB_form*Wnt_inh-k_OB_apop*OB;
dxdt_OC=k_OC_form*RANKL_eff/(1.0+OPG_eff)-k_OC_apop*OC*(1.0+E_ZOL);
dxdt_BV=k_bone_form*OB-k_bone_res*OC;
dxdt_NTX=k_bone_res*OC-0.5*NTX;
dxdt_PINP=k_bone_form*OB-0.5*PINP;
dxdt_Hgb=kHgb_recover*(Hgb0-Hgb)/(1.0+IL6_Hgb_inh*IL6)-kHgb_supp*MM_burden;
dxdt_B2M=kB2M_prod*N_total-kB2M_elim*B2M;
dxdt_BTZ1=-(BTZ_CL+BTZ_Q2+BTZ_Q3)/BTZ_V1*BTZ1+BTZ_Q2/BTZ_V2*BTZ2+BTZ_Q3/BTZ_V3*BTZ3;
dxdt_BTZ2=BTZ_Q2/BTZ_V1*BTZ1-BTZ_Q2/BTZ_V2*BTZ2;
dxdt_BTZ3=BTZ_Q3/BTZ_V1*BTZ1-BTZ_Q3/BTZ_V3*BTZ3;
dxdt_LEN_gut=-LEN_Ka*LEN_gut;
dxdt_LEN1=LEN_Ka*LEN_gut*LEN_F-(LEN_CL/LEN_V)*LEN1;
double CD38_free=CD38_0-DARA_CD38; if(CD38_free<0) CD38_free=0;
dxdt_DARA1=-(DARA_CL+DARA_Q)*70.0/DARA_V1/70.0*DARA1+DARA_Q*70.0/DARA_V2/70.0*DARA2-DARA_kon*DARA1*CD38_free+DARA_koff*DARA_CD38;
dxdt_DARA2=DARA_Q*70.0/DARA_V1/70.0*DARA1-DARA_Q*70.0/DARA_V2/70.0*DARA2;
dxdt_DARA_CD38=DARA_kon*DARA1*CD38_free-(DARA_koff+DARA_kdeg)*DARA_CD38;
dxdt_DEX_gut=-DEX_Ka*DEX_gut;
dxdt_DEX1=DEX_Ka*DEX_gut*DEX_F-(DEX_CL/DEX_V)*DEX1;
dxdt_VEN_gut=-VEN_Ka*VEN_gut;
dxdt_VEN1=VEN_Ka*VEN_gut*VEN_F-(VEN_CL/VEN_V)*VEN1;
dxdt_ZOL1=-(ZOL_CL+ZOL_Q)/ZOL_V1*ZOL1+ZOL_Q/ZOL_V2*ZOL2-ZOL_kbone*ZOL1;
dxdt_ZOL2=ZOL_Q/ZOL_V1*ZOL1-ZOL_Q/ZOL_V2*ZOL2;
dxdt_ZOL_bone=ZOL_kbone*ZOL1-ZOL_krel*ZOL_bone;
$TABLE
double MM_total=MM_S+MM_R;
double BTZ_Cp_out=BTZ1/BTZ_V1;
double LEN_Cp_out=LEN1/LEN_V;
double DARA_Cp_out=DARA1/(DARA_V1*70.0);
double DEX_Cp_out=DEX1/DEX_V;
double VEN_Cp_out=VEN1/VEN_V;
double E_BTZ_out=BTZ_Emax*pow(BTZ_Cp_out,BTZ_hill)/(pow(BTZ_EC50,BTZ_hill)+pow(BTZ_Cp_out,BTZ_hill));
double E_LEN_out=LEN_Emax*pow(LEN_Cp_out,LEN_hill)/(pow(LEN_EC50,LEN_hill)+pow(LEN_Cp_out,LEN_hill));
double E_DARA_out=DARA_Emax*pow(DARA_Cp_out,DARA_hill)/(pow(DARA_EC50,DARA_hill)+pow(DARA_Cp_out,DARA_hill));
$CAPTURE
MM_S MM_R MM_total MP FLC IL6 OB OC BV NTX PINP Hgb B2M
BTZ_Cp_out LEN_Cp_out DARA_Cp_out DEX_Cp_out VEN_Cp_out
E_BTZ_out E_LEN_out E_DARA_out
'

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "Multiple Myeloma QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",  tabName="tab_patient",  icon=icon("user")),
      menuItem("PK Profiles",      tabName="tab_pk",       icon=icon("chart-line")),
      menuItem("PD / Tumor",       tabName="tab_pd",       icon=icon("dna")),
      menuItem("Clinical Endpoints", tabName="tab_endpoints", icon=icon("stethoscope")),
      menuItem("Scenario Comparison", tabName="tab_scenarios", icon=icon("balance-scale")),
      menuItem("Biomarkers",       tabName="tab_biomarkers", icon=icon("vial"))
    )
  ),

  dashboardBody(
    tabItems(

      ## ---- TAB 1: PATIENT PROFILE ----
      tabItem(tabName="tab_patient",
        fluidRow(
          box(title="Patient Parameters", status="primary", solidHeader=TRUE, width=4,
            sliderInput("bw", "Body Weight (kg)", 40, 130, 70, step=5),
            sliderInput("age", "Age (years)", 30, 90, 65, step=1),
            selectInput("iss_stage", "ISS Stage",
                        choices=c("ISS I"=1, "ISS II"=2, "ISS III"=3), selected=1),
            selectInput("cytogenetics", "Cytogenetic Risk",
                        choices=c("Standard Risk"="std",
                                  "High Risk (del17p/t(4;14))"="high",
                                  "t(11;14) / BCL-2 high"="t1114"),
                        selected="std"),
            sliderInput("baseline_mp", "Baseline M-Protein (g/dL)", 0.5, 8.0, 3.0, step=0.5),
            numericInput("baseline_b2m", "Baseline β₂M (mg/L)", value=3.5, min=1, max=20),
            sliderInput("baseline_hgb", "Baseline Hgb (g/dL)", 7.0, 16.0, 12.0, step=0.5),
            actionButton("run_sim", "Run Simulation", class="btn-primary btn-lg")
          ),
          box(title="Disease Summary at Baseline", status="warning", solidHeader=TRUE, width=8,
            infoBoxOutput("iss_box"),
            infoBoxOutput("stage_box"),
            infoBoxOutput("risk_box"),
            DTOutput("patient_table")
          )
        )
      ),

      ## ---- TAB 2: PK PROFILES ----
      tabItem(tabName="tab_pk",
        fluidRow(
          box(title="Regimen Selection", status="primary", solidHeader=TRUE, width=3,
            selectInput("regimen", "Treatment Regimen",
                        choices=c("VRd (Bortezomib+Len+Dex)"="VRd",
                                  "DRd (Daratumumab+Len+Dex)"="DRd",
                                  "KRd (Carfilzomib+Len+Dex)"="KRd",
                                  "VenDex (Venetoclax+Dex)"="VenDex",
                                  "DVRd (Dara+VRd)"="DVRd"),
                        selected="VRd"),
            selectInput("pk_drug", "Drug to Display",
                        choices=c("Bortezomib"="BTZ",
                                  "Lenalidomide"="LEN",
                                  "Daratumumab"="DARA",
                                  "Dexamethasone"="DEX",
                                  "Venetoclax"="VEN"),
                        selected="BTZ"),
            sliderInput("n_cycles_pk", "Number of Cycles", 1, 12, 4, step=1),
            numericInput("BTZ_dose", "Bortezomib Dose (mg/m²)", value=1.3, min=0.7, max=1.5, step=0.1),
            numericInput("LEN_dose", "Lenalidomide Dose (mg)", value=25, min=5, max=25, step=5),
            numericInput("DARA_dose_mgkg", "Daratumumab (mg/kg)", value=16, min=8, max=24, step=4)
          ),
          box(title="PK Concentration vs Time", status="info", solidHeader=TRUE, width=9,
            plotlyOutput("pk_plot", height="400px"),
            plotlyOutput("pk_cumulative", height="300px")
          )
        )
      ),

      ## ---- TAB 3: PD / TUMOR BURDEN ----
      tabItem(tabName="tab_pd",
        fluidRow(
          box(title="PD Parameters", status="danger", solidHeader=TRUE, width=3,
            sliderInput("kgrow", "MM Cell Growth Rate (/day)", 0.005, 0.035, 0.015, step=0.002),
            sliderInput("BTZ_Emax_input", "BTZ Emax", 0.5, 0.99, 0.90, step=0.05),
            sliderInput("BTZ_EC50_input", "BTZ EC50 (nM)", 0.5, 10, 2.5, step=0.5),
            sliderInput("LEN_Emax_input", "Lenalidomide Emax", 0.3, 0.90, 0.65, step=0.05),
            sliderInput("LEN_EC50_input", "LEN EC50 (µg/mL)", 0.1, 3.0, 0.8, step=0.1),
            sliderInput("DARA_Emax_input", "Daratumumab Emax", 0.4, 0.99, 0.85, step=0.05),
            checkboxInput("include_res", "Include Resistant Clone", value=TRUE),
            sliderInput("sim_duration_pd", "Simulation Duration (days)", 180, 1080, 720, step=90)
          ),
          box(title="MM Cell Dynamics", status="danger", solidHeader=TRUE, width=9,
            plotlyOutput("tumor_plot", height="360px"),
            plotlyOutput("mprotein_pd_plot", height="300px")
          )
        )
      ),

      ## ---- TAB 4: CLINICAL ENDPOINTS ----
      tabItem(tabName="tab_endpoints",
        fluidRow(
          box(title="Response Criteria", status="success", solidHeader=TRUE, width=12,
            plotlyOutput("response_swim", height="350px"),
            fluidRow(
              column(6, plotlyOutput("bone_plot", height="320px")),
              column(6, plotlyOutput("hgb_plot", height="320px"))
            )
          )
        ),
        fluidRow(
          box(title="Response Summary Table", status="success", width=12,
            DTOutput("response_table"))
        )
      ),

      ## ---- TAB 5: SCENARIO COMPARISON ----
      tabItem(tabName="tab_scenarios",
        fluidRow(
          box(title="Select Regimens to Compare", status="primary", solidHeader=TRUE, width=3,
            checkboxGroupInput("selected_scenarios", "Regimens:",
                               choices=c("No Treatment"="none",
                                         "VRd"="VRd",
                                         "DRd"="DRd",
                                         "KRd"="KRd",
                                         "VenDex (t(11;14))"="VenDex",
                                         "DVRd High-risk"="DVRd"),
                               selected=c("none","VRd","DRd")),
            sliderInput("comp_duration", "Duration (days)", 180, 1440, 720, step=90),
            selectInput("comp_endpoint", "Primary Endpoint",
                        choices=c("M-Protein (g/dL)"="MP",
                                  "Total MM Cells (log10)"="logMM",
                                  "Bone Volume (%)"="BV",
                                  "Hemoglobin (g/dL)"="Hgb",
                                  "β₂-Microglobulin (mg/L)"="B2M"),
                        selected="MP")
          ),
          box(title="Comparative Efficacy Plot", status="info", solidHeader=TRUE, width=9,
            plotlyOutput("scenario_plot", height="450px"),
            plotlyOutput("waterfall_plot", height="300px")
          )
        )
      ),

      ## ---- TAB 6: BIOMARKERS ----
      tabItem(tabName="tab_biomarkers",
        fluidRow(
          box(title="Bone Biomarkers", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("bone_biomarker_plot", height="350px")
          ),
          box(title="Cytokines & Immune Markers", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("cytokine_plot", height="350px")
          )
        ),
        fluidRow(
          box(title="Response Biomarker Correlations", status="info", width=12,
            plotlyOutput("biomarker_corr", height="350px"))
        )
      )
    )
  )
)

## ============================================================
## SERVER
## ============================================================
server <- function(input, output, session) {

  # Load model once
  mod <- reactive({
    mread_cache("mm_shiny", tempdir(), mm_code_shiny)
  })

  # ---- Patient info boxes ----
  output$iss_box <- renderInfoBox({
    iss <- as.integer(input$iss_stage)
    color <- c("green","yellow","red")[iss]
    infoBox("ISS Stage", paste("Stage", iss),
            icon=icon("layer-group"), color=color)
  })

  output$stage_box <- renderInfoBox({
    mp <- input$baseline_mp
    resp <- ifelse(mp > 3.0, "High M-protein", "Moderate M-protein")
    infoBox("M-Protein", paste0(mp, " g/dL"), subtitle=resp,
            icon=icon("vial"), color="purple")
  })

  output$risk_box <- renderInfoBox({
    risk_label <- c(
      std="Standard Risk", high="High Risk", t1114="t(11;14) BCL-2 high"
    )[input$cytogenetics]
    risk_color <- c(std="green", high="red", t1114="orange")[input$cytogenetics]
    infoBox("Cytogenetics", risk_label, icon=icon("dna"), color=risk_color)
  })

  output$patient_table <- renderDT({
    df <- data.frame(
      Parameter = c("Body Weight","Age","ISS Stage","Cytogenetics",
                    "M-Protein","β₂M","Hemoglobin","Risk Category"),
      Value     = c(paste0(input$bw, " kg"),
                    paste0(input$age, " yr"),
                    paste0("Stage ", input$iss_stage),
                    input$cytogenetics,
                    paste0(input$baseline_mp, " g/dL"),
                    paste0(input$baseline_b2m, " mg/L"),
                    paste0(input$baseline_hgb, " g/dL"),
                    ifelse(input$cytogenetics=="high","High","Standard/Intermediate"))
    )
    datatable(df, options=list(dom='t', pageLength=10), rownames=FALSE)
  })

  # ---- Core simulation reactive ----
  sim_data <- eventReactive(input$run_sim, {
    m <- mod()

    # Adjust parameters based on patient profile
    growth_adj <- switch(input$cytogenetics,
                         std=1.0, high=1.3, t1114=0.9)
    ven_adj <- ifelse(input$cytogenetics=="t1114", 1.3, 1.0)

    m <- m %>%
      param(kgrow      = input$kgrow * growth_adj,
            VEN_Emax   = 0.80 * ven_adj,
            BTZ_Emax   = input$BTZ_Emax_input,
            BTZ_EC50   = input$BTZ_EC50_input,
            LEN_Emax   = input$LEN_Emax_input,
            LEN_EC50   = input$LEN_EC50_input,
            DARA_Emax  = input$DARA_Emax_input,
            Hgb0       = input$baseline_hgb)

    # Build events for selected regimen
    bsa <- sqrt(input$bw * 1.72 / 3600)  # DuBois approximation
    btz_dose_ng <- input$BTZ_dose * bsa * 1e3

    days_btz <- unlist(lapply(0:7, function(c) c*21 + c(1,4,8,11)))
    ev_btz <- ev(amt=btz_dose_ng, cmt="BTZ1", time=days_btz*24)

    len_days <- unlist(lapply(0:7, function(c) c*28 + 0:20))
    ev_len  <- ev(amt=input$LEN_dose*1e3*0.9, cmt="LEN_gut", time=len_days*24)

    dex_days <- unlist(lapply(0:7, function(c) c*28 + c(1,8,15,22)))
    ev_dex  <- ev(amt=40e3*0.78, cmt="DEX_gut", time=dex_days*24)

    dara_times <- c(seq(0, 7*7, by=7),
                    56 + seq(0, 7*14, by=14),
                    56+8*14 + seq(0, 11*28, by=28)) * 24
    ev_dara <- ev(amt=16e3*input$bw, cmt="DARA1", time=dara_times)

    ev_ven <- ev(amt=800e3*0.55, cmt="VEN_gut",
                  time=seq(0, input$sim_duration_pd)*24)

    ev_list <- switch(input$regimen,
      VRd    = as_data_frame(bind_rows(ev_btz %>% as_data_frame(),
                                        ev_len %>% as_data_frame(),
                                        ev_dex %>% as_data_frame())),
      DRd    = as_data_frame(bind_rows(ev_dara %>% as_data_frame(),
                                        ev_len  %>% as_data_frame(),
                                        ev_dex  %>% as_data_frame())),
      KRd    = as_data_frame(bind_rows(
                  ev(amt=36e3, cmt="CFZ1",
                     time=c(1,2,8,9,15,16)*24, ii=21*24, addl=7) %>%
                     as_data_frame(),
                  ev_len %>% as_data_frame(),
                  ev_dex %>% as_data_frame())),
      VenDex = as_data_frame(bind_rows(ev_ven %>% as_data_frame(),
                                        ev_dex %>% as_data_frame())),
      DVRd   = as_data_frame(bind_rows(ev_dara %>% as_data_frame(),
                                        ev_btz  %>% as_data_frame(),
                                        ev_len  %>% as_data_frame(),
                                        ev_dex  %>% as_data_frame()))
    )

    result <- m %>%
      data_set(ev_list) %>%
      mrgsim(end=input$sim_duration_pd*24, delta=24) %>%
      as_tibble() %>%
      mutate(day=time/24,
             logMM=log10(MM_total+1),
             MP_pct_change = (MP - input$baseline_mp) / input$baseline_mp * 100)

    result
  })

  # ---- PK plot ----
  output$pk_plot <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    cp_col <- switch(input$pk_drug,
                     BTZ="BTZ_Cp_out", LEN="LEN_Cp_out",
                     DARA="DARA_Cp_out", DEX="DEX_Cp_out", VEN="VEN_Cp_out")
    unit <- switch(input$pk_drug,
                   BTZ="nM", LEN="µg/mL", DARA="µg/mL", DEX="µg/mL", VEN="µg/mL")
    ec50 <- switch(input$pk_drug,
                   BTZ=2.5, LEN=0.8, DARA=0.15, DEX=0.10, VEN=0.50)

    df$cp <- df[[cp_col]]
    p <- plot_ly(df, x=~day, y=~cp, type='scatter', mode='lines',
                 line=list(color='#7e22ce', width=2),
                 name=input$pk_drug) %>%
      add_trace(y=rep(ec50, nrow(df)), type='scatter', mode='lines',
                line=list(dash='dash', color='red', width=1),
                name=paste("EC50 =", ec50, unit)) %>%
      layout(title=paste(input$pk_drug, "PK —", input$regimen),
             xaxis=list(title="Day"),
             yaxis=list(title=paste("Concentration (", unit, ")")))
    p
  })

  output$pk_cumulative <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    plot_ly(df, x=~day, y=~E_BTZ_out, type='scatter', mode='lines',
            name="BTZ Effect", fill='tozeroy', fillcolor='rgba(126,34,206,0.2)') %>%
      add_trace(y=~E_LEN_out, name="LEN Effect",
                fill='tozeroy', fillcolor='rgba(0,119,182,0.2)') %>%
      add_trace(y=~E_DARA_out, name="DARA Effect",
                fill='tozeroy', fillcolor='rgba(214,40,40,0.2)') %>%
      layout(title="Drug Effect (Emax) over Time",
             xaxis=list(title="Day"),
             yaxis=list(title="Effect (0–1)", range=c(0,1)))
  })

  # ---- PD / Tumor plots ----
  output$tumor_plot <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    plot_ly(df, x=~day) %>%
      add_trace(y=~log10(MM_S+1), type='scatter', mode='lines',
                name="Sensitive Clone", line=list(color='#d62828')) %>%
      add_trace(y=~log10(MM_R+1), type='scatter', mode='lines',
                name="Resistant Clone", line=list(color='#f77f00')) %>%
      add_trace(y=~log10(MM_total+1), type='scatter', mode='lines',
                name="Total MM", line=list(color='#003049', width=2)) %>%
      layout(title="MM Cell Dynamics (log10)",
             xaxis=list(title="Day"),
             yaxis=list(title="log10(Cells + 1)"))
  })

  output$mprotein_pd_plot <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    ref_mp <- input$baseline_mp
    plot_ly(df, x=~day, y=~MP, type='scatter', mode='lines',
            line=list(color='#0077b6', width=2), name="M-Protein") %>%
      add_trace(y=rep(ref_mp * 0.5, nrow(df)), type='scatter', mode='lines',
                line=list(dash='dash', color='orange'), name="50% reduction (PR)") %>%
      add_trace(y=rep(ref_mp * 0.1, nrow(df)), type='scatter', mode='lines',
                line=list(dash='dash', color='green'), name="90% reduction (VGPR)") %>%
      layout(title="M-Protein Response",
             xaxis=list(title="Day"),
             yaxis=list(title="M-Protein (g/dL)"))
  })

  # ---- Endpoints plots ----
  output$response_swim <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    df$resp <- with(df, case_when(
      MP < input$baseline_mp * 0.1   ~ "CR/sCR",
      MP < input$baseline_mp * 0.3   ~ "VGPR",
      MP < input$baseline_mp * 0.5   ~ "PR",
      MP < input$baseline_mp         ~ "SD",
      TRUE                           ~ "PD"
    ))
    color_map <- c("CR/sCR"="#003366","VGPR"="#0066cc",
                   "PR"="#66b2ff","SD"="#ffcc00","PD"="#cc0000")
    df$resp_num <- as.numeric(factor(df$resp,
                    levels=c("PD","SD","PR","VGPR","CR/sCR")))
    plot_ly(df, x=~day, y=~MP, type='scatter', mode='lines+markers',
            color=~resp, colors=color_map,
            marker=list(size=4)) %>%
      layout(title="Response over Time (Swimmer Plot equivalent)",
             xaxis=list(title="Day"),
             yaxis=list(title="M-Protein (g/dL)"))
  })

  output$bone_plot <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    plot_ly(df, x=~day) %>%
      add_trace(y=~BV, name="Bone Volume (%)", type='scatter', mode='lines',
                line=list(color='#8B4513')) %>%
      add_trace(y=~NTX*50, name="NTX (bone resorption) ×50",
                type='scatter', mode='lines',
                line=list(color='red', dash='dash')) %>%
      add_trace(y=~PINP*50, name="P1NP (bone formation) ×50",
                type='scatter', mode='lines',
                line=list(color='green', dash='dash')) %>%
      layout(title="Bone Biomarkers",
             xaxis=list(title="Day"),
             yaxis=list(title="Value"))
  })

  output$hgb_plot <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    plot_ly(df, x=~day, y=~Hgb, type='scatter', mode='lines',
            line=list(color='#dc143c', width=2), name="Hemoglobin") %>%
      add_trace(y=rep(10, nrow(df)), type='scatter', mode='lines',
                line=list(dash='dash', color='orange'), name="Anemia threshold") %>%
      layout(title="Hemoglobin over Time",
             xaxis=list(title="Day"), yaxis=list(title="Hgb (g/dL)"))
  })

  output$response_table <- renderDT({
    req(sim_data())
    df <- sim_data()
    t_pts <- c(90, 180, 360, 540, 720)
    t_pts <- t_pts[t_pts <= input$sim_duration_pd]
    df %>%
      filter(day %in% t_pts) %>%
      mutate(
        Response = case_when(
          MP < input$baseline_mp * 0.1 ~ "CR/sCR",
          MP < input$baseline_mp * 0.3 ~ "VGPR",
          MP < input$baseline_mp * 0.5 ~ "PR",
          MP < input$baseline_mp       ~ "SD",
          TRUE                         ~ "PD"
        ),
        `MP Change %` = round((MP - input$baseline_mp) / input$baseline_mp * 100, 1),
        MP = round(MP, 2),
        FLC = round(FLC, 1),
        Hgb = round(Hgb, 1),
        B2M = round(B2M, 2),
        BV  = round(BV, 1)
      ) %>%
      select(Day=day, `M-Protein`=MP, `MP Δ%`=`MP Change %`,
             sFLC=FLC, Hgb, `β₂M`=B2M, `BV%`=BV, Response) %>%
      datatable(options=list(pageLength=10), rownames=FALSE) %>%
      formatStyle("Response",
                  backgroundColor=styleEqual(
                    c("CR/sCR","VGPR","PR","SD","PD"),
                    c("#003366","#0066cc","#66b2ff","#ffcc00","#cc0000")),
                  color=styleEqual(
                    c("CR/sCR","VGPR","PR","SD","PD"),
                    c("white","white","black","black","white")))
  })

  # ---- Scenario comparison ----
  output$scenario_plot <- renderPlotly({
    req(sim_data())
    req(input$selected_scenarios)
    # Placeholder — in production, run all selected scenarios
    df <- sim_data()
    y_var <- input$comp_endpoint
    if(y_var == "logMM") df$y_plot <- log10(df$MM_total + 1)
    else df$y_plot <- df[[y_var]]

    plot_ly(df, x=~day, y=~y_plot, type='scatter', mode='lines',
            name=input$regimen, line=list(width=2)) %>%
      layout(title=paste("Scenario Comparison —", y_var),
             xaxis=list(title="Day"),
             yaxis=list(title=y_var))
  })

  output$waterfall_plot <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    # Best response waterfall at nadir
    nadir_mp <- min(df$MP, na.rm=TRUE)
    best_resp_pct <- (nadir_mp - input$baseline_mp) / input$baseline_mp * 100

    plot_ly(x=input$regimen, y=best_resp_pct,
            type='bar', marker=list(color=ifelse(best_resp_pct < -90, '#003366',
                                    ifelse(best_resp_pct < -50, '#0066cc',
                                           ifelse(best_resp_pct < 0, '#66b2ff', '#cc0000'))))) %>%
      add_trace(x=c("PR line","VGPR line"),
                y=c(-50, -90), type='scatter', mode='lines+markers',
                line=list(dash='dash', color='orange')) %>%
      layout(title="Best Response (M-Protein % Change from Baseline)",
             xaxis=list(title="Regimen"),
             yaxis=list(title="Best % Change from Baseline",
                        zeroline=TRUE))
  })

  # ---- Biomarker plots ----
  output$bone_biomarker_plot <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    plot_ly(df, x=~day) %>%
      add_trace(y=~OB, name="Osteoblasts (AU)", type='scatter', mode='lines',
                line=list(color='green')) %>%
      add_trace(y=~OC, name="Osteoclasts (AU)", type='scatter', mode='lines',
                line=list(color='red')) %>%
      add_trace(y=~BV, name="Bone Volume (%)", type='scatter', mode='lines',
                line=list(color='brown', width=2)) %>%
      layout(title="Bone Remodeling Dynamics",
             xaxis=list(title="Day"),
             yaxis=list(title="Level (AU or %)"))
  })

  output$cytokine_plot <- renderPlotly({
    req(sim_data())
    df <- sim_data()
    plot_ly(df, x=~day) %>%
      add_trace(y=~IL6, name="IL-6 (AU)", type='scatter', mode='lines',
                line=list(color='#f77f00')) %>%
      add_trace(y=~B2M, name="β₂M (mg/L)", type='scatter', mode='lines',
                line=list(color='#7e22ce')) %>%
      layout(title="Cytokine & Immunological Biomarkers",
             xaxis=list(title="Day"),
             yaxis=list(title="Concentration (AU or mg/L)"))
  })

  output$biomarker_corr <- renderPlotly({
    req(sim_data())
    df <- sim_data() %>% filter(day %% 30 == 0)
    plot_ly(df, x=~MP, y=~B2M, type='scatter', mode='markers',
            marker=list(color=~day, colorscale='Viridis', size=8,
                        colorbar=list(title='Day')),
            text=~paste("Day:", day, "<br>MP:", round(MP,2),
                        "<br>β₂M:", round(B2M,2))) %>%
      layout(title="M-Protein vs β₂M Correlation (by treatment day)",
             xaxis=list(title="M-Protein (g/dL)"),
             yaxis=list(title="β₂-Microglobulin (mg/L)"))
  })
}

## ============================================================
## LAUNCH
## ============================================================
shinyApp(ui=ui, server=server)
