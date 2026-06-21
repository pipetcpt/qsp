## ============================================================
## HIV/AIDS QSP Shiny Dashboard
## 6 Tabs: Patient Profile · Viral Kinetics · CD4/Immune ·
##         Drug PK · Scenario Comparison · Reservoir & Resistance
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(tidyr)
library(DT)

# ─────────────────────────────────────────────────────────────
# mrgsolve Model (embedded)
# ─────────────────────────────────────────────────────────────
hiv_code <- '
$PARAM
  use_TDF=0,use_FTC=0,use_DTG=0,use_BIC=0,use_EFV=0,use_DRV=0,
  ka_TDF=1.0,CL_TDF=16.0,Vc_TDF=45.0,k_IC_TDF=0.010,k_eg_TDF=0.00462,
  ka_FTC=1.5,CL_FTC=13.7,Vc_FTC=213.0,k_IC_FTC=0.020,k_eg_FTC=0.01782,
  ka_DTG=0.9,CL_DTG=1.0,Vc_DTG=17.4,
  ka_BIC=0.5,CL_BIC=0.41,Vc_BIC=25.0,
  ka_EFV=0.5,CL_EFV=4.5,Vc_EFV=270.0,
  ka_DRV=1.2,CL_DRV=1.0,Vc_DRV=88.4,
  IC50_TFV=100.0,IC50_FTC=1000.0,IC50_DTG=2.0,
  IC50_BIC=1.7,IC50_EFV=3000.0,IC50_DRV=0.5,n_H=1.5,
  s_T=10.0,d_T=0.01,beta=5e-8,delta_I=1.0,
  p_V=3000.0,c_V=23.0,k_lat=0.02,r_lat=1e-5,
  s_E=5.0,p_E=0.5,d_E=0.05,k_kill=0.004,
  k_IL6=0.01,d_IL6=0.10,IL6_0=2.0,
  k_mut=1e-6,k_rd=0.001

$CMT TDF_GUT TDF_PLASMA TFV_DP FTC_PLASMA FTC_TP
     DTG_PLASMA BIC_PLASMA EFV_PLASMA DRV_PLASMA
     T_CELL I_CELL V_FREE L_CELL E_CELL
     INFLAM RESIST VL_LOG CD4_SMOOTH

$INIT TDF_GUT=0,TDF_PLASMA=0,TFV_DP=0,FTC_PLASMA=0,FTC_TP=0,
      DTG_PLASMA=0,BIC_PLASMA=0,EFV_PLASMA=0,DRV_PLASMA=0,
      T_CELL=1000,I_CELL=0.01,V_FREE=1.0,L_CELL=0.001,E_CELL=100,
      INFLAM=2.0,RESIST=0,VL_LOG=3.0,CD4_SMOOTH=1000

$ODE
  double ke_TDF=CL_TDF/Vc_TDF;
  dxdt_TDF_GUT   =-ka_TDF*TDF_GUT;
  dxdt_TDF_PLASMA= ka_TDF*TDF_GUT/Vc_TDF - ke_TDF*TDF_PLASMA;
  double TFV_nM=TDF_PLASMA*1000.0/287.2;
  dxdt_TFV_DP    = k_IC_TDF*TFV_nM - k_eg_TDF*TFV_DP;
  double ke_FTC=CL_FTC/Vc_FTC;
  dxdt_FTC_PLASMA=-ke_FTC*FTC_PLASMA;
  double FTC_nM=FTC_PLASMA*1000.0/247.2;
  dxdt_FTC_TP    = k_IC_FTC*FTC_nM - k_eg_FTC*FTC_TP;
  double ke_DTG=CL_DTG/Vc_DTG;
  dxdt_DTG_PLASMA=-ke_DTG*DTG_PLASMA;
  double ke_BIC=CL_BIC/Vc_BIC;
  dxdt_BIC_PLASMA=-ke_BIC*BIC_PLASMA;
  double ke_EFV=CL_EFV/Vc_EFV;
  dxdt_EFV_PLASMA=-ke_EFV*EFV_PLASMA;
  double ke_DRV=CL_DRV/Vc_DRV;
  dxdt_DRV_PLASMA=-ke_DRV*DRV_PLASMA;
  double DTG_nM=DTG_PLASMA*1000.0/419.4;
  double BIC_nM=BIC_PLASMA*1000.0/449.4;
  double EFV_nM=EFV_PLASMA*1000.0/315.7;
  double DRV_nM=DRV_PLASMA*1000.0/547.7;
  double eTFV=pow(TFV_DP,n_H)/(pow(IC50_TFV,n_H)+pow(TFV_DP,n_H));
  double eFTC=pow(FTC_TP,n_H)/(pow(IC50_FTC,n_H)+pow(FTC_TP,n_H));
  double eDTG=pow(DTG_nM,n_H)/(pow(IC50_DTG,n_H)+pow(DTG_nM,n_H));
  double eBIC=pow(BIC_nM,n_H)/(pow(IC50_BIC,n_H)+pow(BIC_nM,n_H));
  double eEFV=pow(EFV_nM,n_H)/(pow(IC50_EFV,n_H)+pow(EFV_nM,n_H));
  double eDRV=pow(DRV_nM,n_H)/(pow(IC50_DRV,n_H)+pow(DRV_nM,n_H));
  double eNRTI=1.0-(1.0-eTFV*use_TDF)*(1.0-eFTC*use_FTC);
  double eRT=1.0-(1.0-eNRTI)*(1.0-eEFV*use_EFV);
  double eINSTI=(use_DTG>0.5)?eDTG:((use_BIC>0.5)?eBIC:0.0);
  double eta=1.0-(1.0-eRT)*(1.0-eINSTI)*(1.0-eDRV*use_DRV);
  if(eta>0.9999)eta=0.9999; if(eta<0)eta=0;
  double T=(T_CELL>0)?T_CELL:0.0;
  double I=(I_CELL>0)?I_CELL:0.0;
  double V=(V_FREE>0)?V_FREE:0.0;
  double L=(L_CELL>0)?L_CELL:0.0;
  double E=(E_CELL>0)?E_CELL:0.0;
  double R=(RESIST>0)?RESIST:0.0;
  double beta_eff=beta*(1.0-eta)*(1.0+2.0*R);
  dxdt_T_CELL=s_T - d_T*T - beta_eff*V*T*(1.0-k_lat) + r_lat*L;
  dxdt_I_CELL=beta_eff*V*T*(1.0-k_lat) - delta_I*I - k_kill*E*I;
  double pi_supp=1.0-eDRV*use_DRV*0.5;
  dxdt_V_FREE=p_V*I*pi_supp - c_V*V + p_V*r_lat*L*0.1;
  dxdt_L_CELL=beta_eff*V*T*k_lat - r_lat*L - k_kill*E*L*0.01;
  dxdt_E_CELL=s_E + p_E*I/(0.5+I) - d_E*E;
  dxdt_INFLAM=k_IL6*I - d_IL6*(INFLAM-IL6_0);
  dxdt_RESIST=k_mut*V*(1.0-eta) - k_rd*R;
  double VL_now=(V*1000.0>1.0)?log10(V*1000.0):0.0;
  dxdt_VL_LOG=(VL_now-VL_LOG)*0.5;
  dxdt_CD4_SMOOTH=((T-I)-CD4_SMOOTH)*0.5;

$TABLE
  double C_TFV=TDF_PLASMA, C_FTC=FTC_PLASMA,
         C_DTG=DTG_PLASMA, C_BIC=BIC_PLASMA,
         C_EFV=EFV_PLASMA, C_DRV=DRV_PLASMA;
  double TFV_DP_nM=TFV_DP, FTC_TP_nM=FTC_TP;
  double CD4_count=(T_CELL-I_CELL>0)?(T_CELL-I_CELL):0.0;
  double VL_copies=V_FREE*1000.0;
  double VL_log10=(VL_copies>1.0)?log10(VL_copies):0.0;
  double VL_supp=(VL_copies<50.0)?1.0:0.0;
  double AIDS_risk=(CD4_count<200.0)?1.0:0.0;
  double Lat_IUPM=L_CELL*1e3;
  double IL6_conc=INFLAM;
  double Res_score=RESIST;
  double CTL_count=E_CELL;
  double DTG_n=DTG_PLASMA*1000.0/419.4;
  double BIC_n=BIC_PLASMA*1000.0/449.4;
  double EFV_n=EFV_PLASMA*1000.0/315.7;
  double DRV_n=DRV_PLASMA*1000.0/547.7;
  double eTFV2=pow(TFV_DP,n_H)/(pow(IC50_TFV,n_H)+pow(TFV_DP,n_H));
  double eFTC2=pow(FTC_TP,n_H)/(pow(IC50_FTC,n_H)+pow(FTC_TP,n_H));
  double eDTG2=pow(DTG_n,n_H)/(pow(IC50_DTG,n_H)+pow(DTG_n,n_H));
  double eBIC2=pow(BIC_n,n_H)/(pow(IC50_BIC,n_H)+pow(BIC_n,n_H));
  double eEFV2=pow(EFV_n,n_H)/(pow(IC50_EFV,n_H)+pow(EFV_n,n_H));
  double eDRV2=pow(DRV_n,n_H)/(pow(IC50_DRV,n_H)+pow(DRV_n,n_H));
  double eNRTI2=1.0-(1.0-eTFV2*use_TDF)*(1.0-eFTC2*use_FTC);
  double eRT2=1.0-(1.0-eNRTI2)*(1.0-eEFV2*use_EFV);
  double eI2=(use_DTG>0.5)?eDTG2:((use_BIC>0.5)?eBIC2:0.0);
  double Eta=1.0-(1.0-eRT2)*(1.0-eI2)*(1.0-eDRV2*use_DRV);
  if(Eta>0.9999)Eta=0.9999; if(Eta<0)Eta=0.0;

$CAPTURE
  C_TFV C_FTC C_DTG C_BIC C_EFV C_DRV
  TFV_DP_nM FTC_TP_nM
  CD4_count VL_copies VL_log10 VL_supp AIDS_risk
  Lat_IUPM IL6_conc Res_score Eta CTL_count
'

mod_base <- mrgsolve::mcode("hiv_shiny", hiv_code, quiet=TRUE)

# ─────────────────────────────────────────────────────────────
# Helper: build dosing events
# ─────────────────────────────────────────────────────────────
build_ev <- function(drugs, t_start, t_end) {
  dmap <- list(
    TDF = list(cmt="TDF_GUT",    amt=300,            ii=24),
    FTC = list(cmt="FTC_PLASMA", amt=200/213.0*0.93, ii=24),
    DTG = list(cmt="DTG_PLASMA", amt=50/17.4*0.53,   ii=24),
    BIC = list(cmt="BIC_PLASMA", amt=50/25.0*0.95,   ii=24),
    EFV = list(cmt="EFV_PLASMA", amt=600/270.0*0.42, ii=24),
    DRV = list(cmt="DRV_PLASMA", amt=800/88.4*0.82,  ii=24)
  )
  evs <- lapply(drugs, function(d) {
    p <- dmap[[d]]
    addl <- max(0L, as.integer(floor((t_end - t_start) / p$ii)))
    ev(cmt=p$cmt, amt=p$amt, ii=p$ii, addl=addl, time=t_start)
  })
  do.call(c, evs)
}

# ─────────────────────────────────────────────────────────────
# Simulation function
# ─────────────────────────────────────────────────────────────
run_sim <- function(sim_yrs, regimen, cd4_init=1000, vl_init=10000,
                    art_start=0, sti_on=FALSE, prep=FALSE,
                    custom_CL_TDF=16, custom_CL_DTG=1.0) {
  sim_days <- sim_yrs * 365
  drugs <- switch(regimen,
    "TDF/FTC/DTG"  = c("TDF","FTC","DTG"),
    "TAF/FTC/BIC"  = c("TDF","FTC","BIC"),
    "TDF/FTC/EFV"  = c("TDF","FTC","EFV"),
    "DRV/r+DTG"    = c("DRV","DTG"),
    "TDF/FTC"      = c("TDF","FTC"),
    NULL
  )

  pars <- list(
    use_TDF = as.numeric("TDF" %in% drugs),
    use_FTC = as.numeric("FTC" %in% drugs),
    use_DTG = as.numeric("DTG" %in% drugs),
    use_BIC = as.numeric("BIC" %in% drugs),
    use_EFV = as.numeric("EFV" %in% drugs),
    use_DRV = as.numeric("DRV" %in% drugs),
    CL_TDF  = custom_CL_TDF,
    CL_DTG  = custom_CL_DTG
  )
  if (regimen == "TAF/FTC/BIC") {
    pars$k_IC_TDF <- 0.05
    pars$ka_TDF   <- 2.5
  }

  v0 <- max(vl_init / 1000, 0.001)
  i0 <- max(vl_init / 300000, 1e-6)

  m <- mod_base %>%
    param(pars) %>%
    init(T_CELL=cd4_init, V_FREE=v0, I_CELL=i0,
         L_CELL=if(prep) 1e-9 else 0.001,
         CD4_SMOOTH=cd4_init, VL_LOG=log10(max(vl_init,1)))

  if (is.null(drugs) || regimen == "No ART") {
    out <- mrgsim(m, end=sim_days, delta=1, obsonly=TRUE)
  } else if (sti_on) {
    ev1 <- build_ev(drugs, art_start, sim_days %/% 3)
    ev2 <- build_ev(drugs, sim_days * 2 %/% 3, sim_days)
    out <- mrgsim(m, events=c(ev1, ev2), end=sim_days, delta=1, obsonly=TRUE)
  } else {
    evts <- build_ev(drugs, art_start, sim_days)
    out  <- mrgsim(m, events=evts, end=sim_days, delta=1, obsonly=TRUE)
  }
  as_tibble(out)
}

# ─────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title="HIV/AIDS QSP Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("① Patient Profile",   tabName="tab_patient",  icon=icon("user")),
      menuItem("② Viral Kinetics",    tabName="tab_vk",       icon=icon("virus")),
      menuItem("③ CD4 & Immunity",    tabName="tab_cd4",      icon=icon("shield-halved")),
      menuItem("④ Drug PK",           tabName="tab_pk",       icon=icon("pills")),
      menuItem("⑤ Scenario Compare",  tabName="tab_compare",  icon=icon("chart-bar")),
      menuItem("⑥ Reservoir & Resist",tabName="tab_resist",   icon=icon("dna"))
    ),
    hr(),
    h5("Simulation Settings", style="padding-left:15px; color:#ccc"),
    sliderInput("sim_yrs","Simulation (years)",min=1,max=10,value=5,step=1),
    sliderInput("cd4_init","Baseline CD4 (cells/µL)",min=50,max=1200,value=500,step=50),
    sliderInput("vl_init","Baseline VL (copies/mL)",min=1000,max=1e6,value=50000,step=1000),
    selectInput("regimen","ART Regimen",
                choices=c("No ART","TDF/FTC/DTG","TAF/FTC/BIC","TDF/FTC/EFV","DRV/r+DTG","TDF/FTC"),
                selected="TDF/FTC/DTG"),
    sliderInput("art_start","ART Start (days)",min=0,max=360,value=0,step=30),
    checkboxInput("sti_on","Structured Treatment Interruption",FALSE),
    checkboxInput("prep_mode","PrEP Mode (low VL exposure)",FALSE),
    sliderInput("cl_tdf","TDF CL (L/h)",min=5,max=40,value=16,step=1),
    sliderInput("cl_dtg","DTG CL (L/h)",min=0.3,max=3,value=1.0,step=0.1),
    actionButton("run_btn","▶ Run Simulation",class="btn-danger btn-block")
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .box-header { background-color: #c0392b !important; }
      .content-wrapper, .right-side { background-color: #f9f9f9; }
    "))),
    tabItems(

      # ── Tab 1: Patient Profile ─────────────────────────────
      tabItem(tabName="tab_patient",
        fluidRow(
          valueBoxOutput("vbox_cd4",  width=3),
          valueBoxOutput("vbox_vl",   width=3),
          valueBoxOutput("vbox_risk", width=3),
          valueBoxOutput("vbox_prep", width=3)
        ),
        fluidRow(
          box(title="Patient HIV Status Summary", width=6, solidHeader=TRUE, status="danger",
            tableOutput("tbl_patient")
          ),
          box(title="ART Regimen Pharmacology", width=6, solidHeader=TRUE, status="danger",
            tableOutput("tbl_regimen")
          )
        ),
        fluidRow(
          box(title="HIV Disease Stage (CDC/WHO Classification)", width=12, solidHeader=TRUE, status="warning",
            p("Disease staging is determined by CD4 count and viral load at baseline."),
            tableOutput("tbl_staging")
          )
        )
      ),

      # ── Tab 2: Viral Kinetics ──────────────────────────────
      tabItem(tabName="tab_vk",
        fluidRow(
          box(title="Plasma Viral Load (log₁₀ copies/mL)", width=12, solidHeader=TRUE, status="danger",
            plotlyOutput("plot_vl", height="380px")
          )
        ),
        fluidRow(
          box(title="Viral Kinetics — Biphasic Decline", width=6, solidHeader=TRUE,
            plotlyOutput("plot_vl_phases", height="320px")
          ),
          box(title="Viral Load Suppression Probability", width=6, solidHeader=TRUE,
            plotlyOutput("plot_vl_supp", height="320px")
          )
        )
      ),

      # ── Tab 3: CD4 & Immunity ──────────────────────────────
      tabItem(tabName="tab_cd4",
        fluidRow(
          box(title="CD4⁺ T Cell Count (cells/µL)", width=8, solidHeader=TRUE, status="danger",
            plotlyOutput("plot_cd4", height="380px")
          ),
          box(title="Key Immunologic Thresholds", width=4, solidHeader=TRUE,
            tableOutput("tbl_cd4_thresh")
          )
        ),
        fluidRow(
          box(title="CD8⁺ CTL Dynamics", width=6, solidHeader=TRUE,
            plotlyOutput("plot_ctl", height="300px")
          ),
          box(title="IL-6 Systemic Inflammation", width=6, solidHeader=TRUE,
            plotlyOutput("plot_il6", height="300px")
          )
        )
      ),

      # ── Tab 4: Drug PK ─────────────────────────────────────
      tabItem(tabName="tab_pk",
        fluidRow(
          box(title="Plasma Drug Concentrations (µg/mL)", width=8, solidHeader=TRUE, status="danger",
            plotlyOutput("plot_pk_plasma", height="380px")
          ),
          box(title="PK Parameters", width=4, solidHeader=TRUE,
            tableOutput("tbl_pk_params")
          )
        ),
        fluidRow(
          box(title="Intracellular Active Metabolites (nM)", width=6, solidHeader=TRUE,
            plotlyOutput("plot_pk_ic", height="300px")
          ),
          box(title="Overall ART Efficacy η (%)", width=6, solidHeader=TRUE,
            plotlyOutput("plot_eta", height="300px")
          )
        )
      ),

      # ── Tab 5: Scenario Comparison ─────────────────────────
      tabItem(tabName="tab_compare",
        fluidRow(
          box(title="Multi-Scenario: Viral Load Comparison", width=12, solidHeader=TRUE, status="danger",
            plotlyOutput("plot_comp_vl", height="400px")
          )
        ),
        fluidRow(
          box(title="Multi-Scenario: CD4 Recovery", width=6, solidHeader=TRUE,
            plotlyOutput("plot_comp_cd4", height="300px")
          ),
          box(title="Scenario Summary Table (Week 48)", width=6, solidHeader=TRUE,
            DTOutput("tbl_compare")
          )
        )
      ),

      # ── Tab 6: Reservoir & Resistance ─────────────────────
      tabItem(tabName="tab_resist",
        fluidRow(
          box(title="HIV Latent Reservoir (log₁₀ IUPM)", width=6, solidHeader=TRUE, status="danger",
            plotlyOutput("plot_reservoir", height="350px")
          ),
          box(title="Drug Resistance Score (0–100%)", width=6, solidHeader=TRUE, status="danger",
            plotlyOutput("plot_resist", height="350px")
          )
        ),
        fluidRow(
          box(title="Key Resistance Mutations by Drug Class", width=6, solidHeader=TRUE,
            tableOutput("tbl_mutations")
          ),
          box(title="Time to Virologic Failure Risk", width=6, solidHeader=TRUE,
            plotlyOutput("plot_vf_risk", height="280px")
          )
        )
      )
    )
  )
)

# ─────────────────────────────────────────────────────────────
# SERVER
# ─────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Reactive simulation ───────────────────────────────────
  sim_data <- eventReactive(input$run_btn, {
    withProgress(message="Running HIV simulation...", value=0.3, {
      tryCatch(
        run_sim(
          sim_yrs    = input$sim_yrs,
          regimen    = input$regimen,
          cd4_init   = input$cd4_init,
          vl_init    = input$vl_init,
          art_start  = input$art_start,
          sti_on     = input$sti_on,
          prep       = input$prep_mode,
          custom_CL_TDF = input$cl_tdf,
          custom_CL_DTG = input$cl_dtg
        ),
        error = function(e) {
          showNotification(paste("Simulation error:", e$message), type="error")
          NULL
        }
      )
    })
  }, ignoreNULL=FALSE)

  # Auto-run on start
  observe({ input$run_btn })
  isolate({ sim_data() })

  # ── Multi-scenario (all regimens) ─────────────────────────
  all_scens <- reactive({
    input$run_btn
    withProgress(message="Running all regimens...", value=0.1, {
      regs <- c("No ART","TDF/FTC/DTG","TAF/FTC/BIC","TDF/FTC/EFV","DRV/r+DTG")
      lapply(regs, function(r) {
        tryCatch(
          run_sim(input$sim_yrs, r, input$cd4_init, input$vl_init) %>%
            mutate(regimen=r),
          error=function(e) NULL
        )
      }) %>% bind_rows()
    })
  })

  # ── Tab 1: Patient Profile ─────────────────────────────────
  output$vbox_cd4 <- renderValueBox({
    valueBox(
      value = paste0(input$cd4_init, " cells/µL"),
      subtitle = "Baseline CD4 Count",
      icon = icon("shield-halved"),
      color = if (input$cd4_init < 200) "red" else if (input$cd4_init < 350) "yellow" else "green"
    )
  })

  output$vbox_vl <- renderValueBox({
    valueBox(
      value = paste0(round(log10(input$vl_init), 1), " log₁₀"),
      subtitle = paste0("Baseline VL (", format(input$vl_init, big.mark=","), " c/mL)"),
      icon = icon("virus"),
      color = if (input$vl_init > 100000) "red" else if (input$vl_init > 10000) "yellow" else "blue"
    )
  })

  output$vbox_risk <- renderValueBox({
    stage <- if (input$cd4_init < 200) "AIDS (Stage 3)"
             else if (input$cd4_init < 350) "Symptomatic HIV (Stage 2)"
             else "Asymptomatic HIV (Stage 1)"
    color <- if (input$cd4_init < 200) "red" else if (input$cd4_init < 350) "yellow" else "blue"
    valueBox(stage, "CDC/WHO HIV Stage", icon=icon("hospital"), color=color)
  })

  output$vbox_prep <- renderValueBox({
    valueBox(
      value = input$regimen,
      subtitle = "Selected ART Regimen",
      icon = icon("pills"),
      color = "purple"
    )
  })

  output$tbl_patient <- renderTable({
    data.frame(
      Parameter = c("Baseline CD4 (cells/µL)","Baseline VL (log₁₀)","VL (copies/mL)",
                    "HIV Stage (CDC)","ART Regimen","ART Start (day)",
                    "STI Mode","Simulation Duration"),
      Value = c(
        input$cd4_init,
        round(log10(input$vl_init), 2),
        format(input$vl_init, big.mark=","),
        if (input$cd4_init < 200) "Stage 3 (AIDS)"
          else if (input$cd4_init < 350) "Stage 2"
          else "Stage 1",
        input$regimen,
        input$art_start,
        if (input$sti_on) "Yes (1/3–2/3 off)" else "No",
        paste0(input$sim_yrs, " years")
      )
    )
  })

  output$tbl_regimen <- renderTable({
    data.frame(
      Regimen         = c("TDF/FTC/DTG","TAF/FTC/BIC","TDF/FTC/EFV","DRV/r+DTG","TDF/FTC"),
      Brand           = c("Triumeq (GSK)","Biktarvy (Gilead)","Atripla (BMS/Gilead)",
                          "Salvage (PI+INSTI)","PrEP (Truvada)"),
      `Class`         = c("NRTI×2+INSTI","NRTI×2+INSTI","NRTI×2+NNRTI","PI/r+INSTI","NRTI×2"),
      `Resistance Barrier` = c("High","High","Low (K103N)","High","N/A (PrEP)"),
      `WHO Preferred` = c("Yes (2022)","Yes (2022)","LMIC option","2nd-line","PrEP")
    )
  })

  output$tbl_staging <- renderTable({
    data.frame(
      Stage = c("1 (Asymptomatic)","2 (Symptomatic)","3 (AIDS)","3 (AIDS)"),
      `CD4 (cells/µL)` = c("≥500","200–499","<200","Any"),
      VL = c("Any","Any","Any","VL >50 on ART"),
      Criteria = c("No AIDS-defining illness","Mild unexplained symptoms",
                   "AIDS-defining conditions","Virologic failure on ART"),
      `WHO Guideline` = c("Monitor","Start ART","Start ART urgently","Salvage regimen")
    )
  })

  # ── Tab 2: Viral Kinetics ──────────────────────────────────
  output$plot_vl <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x=time/365, y=VL_log10)) +
      geom_line(color="#c0392b", linewidth=1) +
      geom_hline(yintercept=log10(50), linetype="dashed", color="steelblue") +
      geom_hline(yintercept=log10(200), linetype="dotted", color="orange") +
      annotate("text", x=input$sim_yrs*0.9, y=log10(50)+0.15,
               label="VL < 50 c/mL (suppressed)", size=3, color="steelblue") +
      scale_x_continuous("Time (Years)", breaks=0:input$sim_yrs) +
      scale_y_continuous("Viral Load (log₁₀ copies/mL)", limits=c(0,8)) +
      labs(title=paste0("Viral Load — ", input$regimen)) +
      theme_bw(base_size=11)
    ggplotly(p) %>% layout(hovermode="x unified")
  })

  output$plot_vl_phases <- renderPlotly({
    df <- sim_data(); req(df)
    df14 <- df %>% filter(time <= 90)
    p <- ggplot(df14, aes(x=time, y=VL_log10)) +
      geom_line(color="#c0392b", linewidth=1) +
      scale_x_continuous("Time (Days)", breaks=seq(0,90,14)) +
      scale_y_continuous("log₁₀ VL") +
      labs(title="Biphasic VL Decline (First 90 days)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plot_vl_supp <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x=time/365, y=cumsum(VL_supp)/pmax(1,row_number())*100)) +
      geom_line(color="#27ae60", linewidth=1) +
      scale_x_continuous("Time (Years)") +
      scale_y_continuous("Cumulative % VL Suppressed", limits=c(0,100)) +
      labs(title="Cumulative Virologic Suppression %") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  # ── Tab 3: CD4 & Immunity ──────────────────────────────────
  output$plot_cd4 <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x=time/365, y=CD4_count)) +
      geom_line(color="#2980b9", linewidth=1) +
      geom_hline(yintercept=200, linetype="dashed", color="red",   alpha=0.7) +
      geom_hline(yintercept=500, linetype="dashed", color="green4",alpha=0.7) +
      geom_ribbon(aes(ymin=200, ymax=500), fill="yellow", alpha=0.08) +
      scale_x_continuous("Time (Years)", breaks=0:input$sim_yrs) +
      scale_y_continuous("CD4⁺ Count (cells/µL)", limits=c(0, NA)) +
      labs(title="CD4⁺ T Cell Count Trajectory") +
      theme_bw(base_size=11)
    ggplotly(p) %>% layout(hovermode="x unified")
  })

  output$tbl_cd4_thresh <- renderTable({
    data.frame(
      `CD4 (cells/µL)` = c("<200","200–350","350–500",">500"),
      Interpretation = c("AIDS-defining risk","Start ART (WHO)","Normal low","Normal"),
      Action = c("Urgent ART + OI prophylaxis","Start ART","ART recommended","Maintain ART")
    )
  })

  output$plot_ctl <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x=time/365, y=CTL_count)) +
      geom_line(color="#8e44ad", linewidth=1) +
      scale_x_continuous("Time (Years)") +
      scale_y_continuous("CD8⁺ CTL (cells/µL)") +
      labs(title="CD8⁺ CTL Effector Dynamics") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plot_il6 <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x=time/365, y=IL6_conc)) +
      geom_line(color="#e67e22", linewidth=1) +
      geom_hline(yintercept=2, linetype="dashed", color="grey50") +
      scale_x_continuous("Time (Years)") +
      scale_y_continuous("IL-6 (pg/mL)") +
      labs(title="Systemic Inflammation (IL-6)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  # ── Tab 4: Drug PK ─────────────────────────────────────────
  output$plot_pk_plasma <- renderPlotly({
    df <- sim_data(); req(df)
    pk30 <- df %>% filter(time <= 30) %>%
      select(time, C_TFV, C_FTC, C_DTG, C_BIC, C_EFV, C_DRV) %>%
      pivot_longer(-time, names_to="Drug", values_to="Conc")
    # Only drugs with non-zero concentrations
    pk30 <- pk30 %>% group_by(Drug) %>% filter(max(Conc,na.rm=TRUE)>0) %>% ungroup()
    drug_labs <- c(C_TFV="TFV (plasma)",C_FTC="FTC (plasma)",C_DTG="DTG (plasma)",
                   C_BIC="BIC (plasma)",C_EFV="EFV (plasma)",C_DRV="DRV (plasma)")
    pk30 <- pk30 %>% mutate(Drug=recode(Drug,!!!drug_labs))
    p <- ggplot(pk30, aes(x=time, y=Conc, color=Drug)) +
      geom_line(linewidth=0.9) +
      scale_x_continuous("Time (Days)", breaks=seq(0,30,5)) +
      scale_y_continuous("Plasma Conc. (µg/mL)") +
      scale_color_brewer(palette="Dark2") +
      labs(title="Plasma Drug Concentrations — Days 1–30", color=NULL) +
      theme_bw(base_size=11) +
      theme(legend.position="bottom")
    ggplotly(p)
  })

  output$tbl_pk_params <- renderTable({
    data.frame(
      Drug = c("TDF/TFV","FTC","DTG","BIC","EFV","DRV/r"),
      CL_Lh = c(input$cl_tdf, 13.7, input$cl_dtg, 0.41, 4.5, 1.0),
      Vc_L  = c(45, 213, 17.4, 25, 270, 88.4),
      `t1/2_h` = round(c(45/input$cl_tdf, 213/13.7, 17.4/input$cl_dtg,
                          25/0.41, 270/4.5, 88.4/1.0)*0.693, 1)
    )
  })

  output$plot_pk_ic <- renderPlotly({
    df <- sim_data(); req(df)
    ic <- df %>% filter(time <= 30) %>%
      select(time, TFV_DP_nM, FTC_TP_nM) %>%
      pivot_longer(-time, names_to="Metabolite", values_to="Conc_nM") %>%
      mutate(Metabolite=recode(Metabolite, TFV_DP_nM="TFV-DP (IC)", FTC_TP_nM="FTC-TP (IC)"))
    p <- ggplot(ic, aes(x=time, y=Conc_nM, color=Metabolite)) +
      geom_line(linewidth=0.9) +
      geom_hline(yintercept=100,  linetype="dashed", color="#E66101", alpha=0.7) +
      geom_hline(yintercept=1000, linetype="dashed", color="#5E3C99", alpha=0.7) +
      scale_color_manual(values=c("TFV-DP (IC)"="#E66101","FTC-TP (IC)"="#5E3C99")) +
      scale_x_continuous("Time (Days)", breaks=seq(0,30,5)) +
      scale_y_continuous("Intracellular Conc. (nM)") +
      labs(title="TFV-DP & FTC-TP Intracellular", color=NULL) +
      theme_bw(base_size=11) +
      theme(legend.position="bottom")
    ggplotly(p)
  })

  output$plot_eta <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x=time/365, y=Eta*100)) +
      geom_line(color="#27ae60", linewidth=1) +
      geom_hline(yintercept=95, linetype="dashed", color="grey40") +
      annotate("text", x=input$sim_yrs*0.85, y=96,
               label="95% efficacy threshold", size=3, color="grey40") +
      scale_x_continuous("Time (Years)") +
      scale_y_continuous("ART Efficacy η (%)", limits=c(0,101)) +
      labs(title="Combined ART Antiviral Efficacy") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  # ── Tab 5: Scenario Comparison ─────────────────────────────
  output$plot_comp_vl <- renderPlotly({
    df <- all_scens(); req(nrow(df) > 0)
    p <- ggplot(df, aes(x=time/365, y=VL_log10, color=regimen)) +
      geom_line(linewidth=0.8) +
      geom_hline(yintercept=log10(50), linetype="dashed", color="grey40") +
      scale_x_continuous("Time (Years)", breaks=0:input$sim_yrs) +
      scale_y_continuous("Viral Load (log₁₀ copies/mL)", limits=c(0,8)) +
      scale_color_brewer(palette="Set1") +
      labs(title="Viral Load Comparison — All Regimens", color="Regimen") +
      theme_bw(base_size=11) +
      theme(legend.position="bottom")
    ggplotly(p)
  })

  output$plot_comp_cd4 <- renderPlotly({
    df <- all_scens(); req(nrow(df) > 0)
    p <- ggplot(df, aes(x=time/365, y=CD4_count, color=regimen)) +
      geom_line(linewidth=0.8) +
      geom_hline(yintercept=200, linetype="dashed", color="red", alpha=0.6) +
      scale_x_continuous("Time (Years)") +
      scale_y_continuous("CD4⁺ (cells/µL)") +
      scale_color_brewer(palette="Set1") +
      labs(title="CD4 Recovery — All Regimens", color="Regimen") +
      theme_bw(base_size=11) +
      theme(legend.position="bottom")
    ggplotly(p)
  })

  output$tbl_compare <- renderDT({
    df <- all_scens(); req(nrow(df) > 0)
    df %>%
      group_by(regimen) %>%
      summarise(
        `VL Week 48 (log₁₀)` = round(VL_log10[which.min(abs(time-336))],2),
        `CD4 Week 48`         = round(CD4_count[which.min(abs(time-336))]),
        `% Suppressed`        = round(mean(VL_supp,na.rm=TRUE)*100,1),
        `CD4 Nadir`           = round(min(CD4_count,na.rm=TRUE)),
        `AIDS Risk`           = max(AIDS_risk,na.rm=TRUE),
        .groups="drop"
      ) %>%
      datatable(rownames=FALSE, options=list(dom="t",ordering=FALSE),
                class="compact stripe hover") %>%
      formatStyle("AIDS Risk", backgroundColor=styleInterval(0.5,c("white","#ffcccc")))
  })

  # ── Tab 6: Reservoir & Resistance ─────────────────────────
  output$plot_reservoir <- renderPlotly({
    df <- sim_data(); req(df)
    df$lat <- pmax(df$Lat_IUPM, 1e-6)
    p <- ggplot(df, aes(x=time/365, y=log10(lat))) +
      geom_line(color="#8e44ad", linewidth=1) +
      geom_hline(yintercept=log10(1), linetype="dashed", color="grey50") +
      annotate("text", x=input$sim_yrs*0.85, y=0.1,
               label="1 IUPM (functional cure threshold)", size=3, color="grey40") +
      scale_x_continuous("Time (Years)") +
      scale_y_continuous("log₁₀ IUPM") +
      labs(title="HIV Latent Reservoir Decay on ART") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plot_resist <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x=time/365, y=Res_score*100)) +
      geom_line(color="#e74c3c", linewidth=1) +
      scale_x_continuous("Time (Years)") +
      scale_y_continuous("Resistance Score (%)", limits=c(0,100)) +
      labs(title="Drug Resistance Mutation Score") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$tbl_mutations <- renderTable({
    data.frame(
      `Drug Class` = c("NRTI","NRTI","NNRTI","INSTI","INSTI","PI"),
      Drug = c("FTC/3TC","TDF/TAF","EFV/NVP","DTG","RAL/EVG","DRV"),
      `Key Mutation` = c("M184V/I","K65R","K103N","R263K","Q148H/G140S","I50V/L76V"),
      `FC (fold-change IC50)` = c(">500×",">10×",">30×","<2× (high barrier)",">50×",">10×"),
      `Barrier to Resistance` = c("Low","Medium","Low","Very High","Medium","High")
    )
  })

  output$plot_vf_risk <- renderPlotly({
    df <- sim_data(); req(df)
    # Virologic failure = VL rebounds >200 after suppression
    df <- df %>%
      mutate(VF_risk = Res_score * (1 - VL_supp))
    p <- ggplot(df, aes(x=time/365, y=VF_risk*100)) +
      geom_area(fill="#e74c3c", alpha=0.4) +
      geom_line(color="#c0392b", linewidth=1) +
      scale_x_continuous("Time (Years)") +
      scale_y_continuous("Virologic Failure Risk Index (%)", limits=c(0,100)) +
      labs(title="Virologic Failure Risk Index") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

}

shinyApp(ui=ui, server=server)
