## ============================================================
## AKI QSP Shiny App — Interactive Dashboard
## Acute Kidney Injury (급성 신손상)
## ============================================================
## Tabs: 1-Patient Profile | 2-Drug PK | 3-Biomarkers |
##       4-Clinical Endpoints | 5-Scenario Comparison | 6-Fibrosis/CKD
## ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mrgsolve)
library(DT)

## ─── Inline mrgsolve model (simplified for Shiny) ─────────────────────────
aki_code <- '
$PARAM
AKI_type=1 isch_sev=0.5 lps_amt=0 tox_dose=0
GFR0=100 TCV0=1 SCR0=0.9
k_inj=0.12 k_apop=0.07 k_nec=0.03 k_ferr=0.02
k_gfr_rec=0.004 k_gfr_inj=0.09
k_backleak=0.015 k_obstruct=0.012
k_il6_prod=0.25 k_il6_deg=0.10
k_tnfa_prod=0.18 k_tnfa_deg=0.14
k_ros_prod=0.18 k_ros_deg=0.14
k_gsh_prod=0.05 k_gsh_dep=0.09
k_atp_prod=0.20 k_atp_dep=0.14
k_ngal_prod=0.45 k_ngal_deg=0.11
k_kim1_prod=0.28 k_kim1_deg=0.07
k_cysC_prod=0.05 k_cysC_deg=0.03
k_repair=0.003 k_egf=0.09 k_wnt=0.04
k_tgfb_prod=0.018 k_tgfb_deg=0.05
k_fibrosis=0.0008 k_myo_act=0.013 k_myo_deg=0.018
Vd_fur1=8 Vd_fur2=14 CL_fur=8 Q_fur=4 F_fur=0.6
NKCC2_IC50=0.5 OAT_eff=1.0
CL_ne=250 Vd_ne=10 Emax_map=0.35 EC50_ne=0.01
CL_nac=15 Vd_nac=30 Emax_gsh=0.55 EC50_nac=10 Emax_no=0.28
CRRT_on=0 CL_crrt=3 CL_crrt_il6=0.5
EC50_ros=0.5 EC50_atp=0.4 n_hill=2

$CMT
FUR_CENT FUR_PERI FUR_GUT
NE_CENT NAC_CENT
ATP ROS GSH TCV GFR IL6 TNFa
NGAL KIM1 SCR CysC
REPAIR_CAP TGFb MYO FIBROSIS

$MAIN
double inj_stim = 0;
if(AKI_type==1) inj_stim = isch_sev;
if(AKI_type==2) inj_stim = tox_dose/20;
if(AKI_type==3) inj_stim = lps_amt/10;
double FUR_TUB = FUR_CENT * OAT_eff * TCV;
double E_NKCC2 = FUR_TUB/(FUR_TUB+NKCC2_IC50);
double E_NE = Emax_map*NE_CENT/(NE_CENT+EC50_ne);
double RBF = 1.0+E_NE;
double E_NAC_gsh = Emax_gsh*NAC_CENT/(NAC_CENT+EC50_nac);
double E_NAC_no  = Emax_no *NAC_CENT/(NAC_CENT+EC50_nac);
double ROS_frac  = ROS/(ROS+EC50_ros);
double ATP_frac  = (1-ATP)/(1-ATP+EC50_atp);
double inj_comp  = 0.5*ROS_frac+0.5*ATP_frac;
double nfkb = 0.5*(inj_stim+ROS+(1-TCV));
double GPx4 = GSH;
double ferr = k_ferr*(1-GPx4)*inj_comp;
double GFR_n = GFR/GFR0;
double CR_elim = (GFR_n*0.08+CRRT_on*CL_crrt/100.0)*SCR;
double repair_s = TCV*REPAIR_CAP*k_egf;
double diuresis = 1+3*E_NKCC2;

$ODE
dxdt_FUR_GUT  = -1.5*FUR_GUT;
dxdt_FUR_CENT = F_fur*1.5*FUR_GUT-(CL_fur/Vd_fur1)*FUR_CENT
                -(Q_fur/Vd_fur1)*FUR_CENT+(Q_fur/Vd_fur2)*FUR_PERI;
dxdt_FUR_PERI = (Q_fur/Vd_fur1)*FUR_CENT-(Q_fur/Vd_fur2)*FUR_PERI;
dxdt_NE_CENT  = -(CL_ne/Vd_ne)*NE_CENT;
dxdt_NAC_CENT = -(CL_nac/Vd_nac)*NAC_CENT;
dxdt_ATP = k_atp_prod*(1-inj_stim)*(1-FIBROSIS)-k_atp_dep*(1+inj_stim+0.5*(1-TCV));
dxdt_ROS = k_ros_prod*(1+2*inj_stim+(1-TCV)+(1-ATP))-k_ros_deg*GSH*ROS-0.2*E_NAC_no*ROS;
dxdt_GSH = k_gsh_prod*(1+E_NAC_gsh)-k_gsh_dep*ROS*GSH;
dxdt_TCV = -(k_apop*inj_comp*TCV*pow(ROS,n_hill)+k_nec*(1-ATP)*TCV+ferr*TCV)
           +repair_s*(1-TCV)*(1-FIBROSIS);
dxdt_GFR = k_gfr_rec*TCV*RBF*(GFR0-GFR)-k_gfr_inj*(1-TCV)*(1+k_backleak+k_obstruct)*GFR;
dxdt_IL6 = k_il6_prod*(1+5*nfkb+lps_amt/5)-k_il6_deg*IL6-CRRT_on*CL_crrt_il6*IL6/5;
dxdt_TNFa= k_tnfa_prod*(1+3*nfkb+lps_amt/8)-k_tnfa_deg*TNFa;
dxdt_NGAL= k_ngal_prod*(1-TCV)*(1+3*inj_stim)-k_ngal_deg*NGAL;
dxdt_KIM1= k_kim1_prod*(1-TCV)*(1+2*inj_stim)-k_kim1_deg*KIM1;
dxdt_SCR = 0.005-CR_elim;
dxdt_CysC= k_cysC_prod-k_cysC_deg*(GFR/GFR0)*CysC;
dxdt_REPAIR_CAP = k_repair*(1-FIBROSIS)*(1+k_wnt)-0.05*(IL6/(IL6+20))*REPAIR_CAP;
dxdt_TGFb= k_tgfb_prod*(1+3*(1-TCV)+nfkb)-k_tgfb_deg*TGFb;
dxdt_MYO = k_myo_act*TGFb*(1-MYO)-k_myo_deg*MYO;
dxdt_FIBROSIS = k_fibrosis*MYO*(1-FIBROSIS);

$TABLE
double AKI_stage = 0;
double sr = SCR/SCR0;
if(sr>=1.5) AKI_stage=1; if(sr>=2.0) AKI_stage=2; if(sr>=3.0) AKI_stage=3;
double eGFR = GFR;
double UO = 1.0*(1+3*E_NKCC2)*TCV*(GFR/GFR0);

$CAPTURE AKI_stage eGFR UO FUR_CENT NE_CENT NAC_CENT
ATP ROS GSH TCV IL6 TNFa NGAL KIM1 SCR CysC
TGFb MYO FIBROSIS REPAIR_CAP diuresis E_NKCC2
'

aki_mod <- mcode("AKI_Shiny", aki_code)
aki_mod <- aki_mod %>%
  init(FUR_CENT=0,FUR_PERI=0,FUR_GUT=0,NE_CENT=0,NAC_CENT=0,
       ATP=1,ROS=0.1,GSH=1,TCV=1,GFR=100,IL6=5,TNFa=3,
       NGAL=5,KIM1=0.5,SCR=0.9,CysC=0.8,
       REPAIR_CAP=1,TGFb=1,MYO=0,FIBROSIS=0)

## Helper: run one scenario
run_aki <- function(mod, params, events=NULL, end=168, delta=0.5) {
  m <- do.call(param, c(list(mod), params))
  if (is.null(events))
    out <- mrgsim(m, end=end, delta=delta)
  else
    out <- mrgsim(m, events=events, end=end, delta=delta)
  as_tibble(out)
}

## ─── UI ────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  dashboardHeader(title="AKI QSP Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("① Patient Profile",   tabName="patient",   icon=icon("user")),
      menuItem("② Drug PK",           tabName="pk",        icon=icon("pills")),
      menuItem("③ Kidney Biomarkers", tabName="biomarkers",icon=icon("flask")),
      menuItem("④ Clinical Endpoints",tabName="endpoints", icon=icon("heartbeat")),
      menuItem("⑤ Scenario Compare",  tabName="compare",   icon=icon("chart-bar")),
      menuItem("⑥ AKI-to-CKD",       tabName="fibrosis",  icon=icon("dna"))
    ),
    hr(),
    h5("Global Settings", style="padding-left:15px; color:#aaa"),
    selectInput("aki_type","AKI Subtype",
                choices=c("Ischemia-Reperfusion (IRI)"=1,
                          "Nephrotoxic (NTX)"=2,
                          "Sepsis-Associated (SA-AKI)"=3), selected=1),
    conditionalPanel("input.aki_type==1",
      sliderInput("isch_sev","Ischemia Severity",0,1,0.5,0.05)),
    conditionalPanel("input.aki_type==2",
      sliderInput("tox_dose","Cisplatin Dose (mg/kg)",0,100,70,5)),
    conditionalPanel("input.aki_type==3",
      sliderInput("lps_amt","LPS / Sepsis Severity (0-10)",0,10,6,0.5)),
    sliderInput("scr0_val","Baseline sCr (mg/dL)",0.5,2.5,0.9,0.1),
    sliderInput("gfr0_val","Baseline GFR",30,130,100,5),
    checkboxInput("crrt_on","CRRT Active",FALSE),
    hr(),
    h5("Treatment", style="padding-left:15px; color:#aaa"),
    checkboxInput("use_fur","Furosemide 40mg IV q12h",FALSE),
    checkboxInput("use_ne","Norepinephrine Infusion",FALSE),
    sliderInput("ne_rate","NE rate (mcg/kg/min)",0,0.5,0.1,0.05),
    checkboxInput("use_nac","NAC IV Prophylaxis",FALSE),
    numericInput("sim_days","Simulation (days)",7,1,30,1)
  ),
  dashboardBody(
    tabItems(
      ## ─── TAB 1: Patient Profile ─────────────────────────────
      tabItem("patient",
        fluidRow(
          box(title="Patient & Disease Summary", width=4, status="primary",
            verbatimTextOutput("pt_summary")),
          box(title="KDIGO AKI Risk Factors", width=4, status="warning",
            tags$ul(
              tags$li("CKD stage ≥ 3 (GFR < 60)"),
              tags$li("Age > 65 years"),
              tags$li("Diabetes mellitus"),
              tags$li("Heart failure / cardiogenic shock"),
              tags$li("Major surgery / CPB"),
              tags$li("Sepsis / ICU admission"),
              tags$li("Nephrotoxin exposure (aminoglycosides, contrast)")
            )
          ),
          box(title="KDIGO AKI Staging", width=4, status="danger",
            tableOutput("kdigo_tbl"))
        ),
        fluidRow(
          box(title="Simulated GFR & sCr — Overview", width=12,
            plotOutput("patient_overview_plot", height=300))
        )
      ),
      ## ─── TAB 2: Drug PK ─────────────────────────────────────
      tabItem("pk",
        fluidRow(
          box(title="Furosemide Plasma Concentration (mg/L)", width=6,
            status="info", plotOutput("fur_pk_plot", height=260)),
          box(title="Norepinephrine & NAC Plasma Levels", width=6,
            status="info", plotOutput("ne_nac_pk_plot", height=260))
        ),
        fluidRow(
          box(title="NKCC2 Inhibition (Furosemide)", width=6,
            plotOutput("nkcc2_plot", height=260)),
          box(title="Diuresis Fold-Increase", width=6,
            plotOutput("diuresis_plot", height=260))
        )
      ),
      ## ─── TAB 3: Biomarkers ──────────────────────────────────
      tabItem("biomarkers",
        fluidRow(
          box(title="NGAL (Early Tubular Injury)", width=6,
            status="warning", plotOutput("ngal_plot", height=260)),
          box(title="KIM-1 (Brush Border Shedding)", width=6,
            status="warning", plotOutput("kim1_plot", height=260))
        ),
        fluidRow(
          box(title="Cystatin C vs Creatinine", width=6,
            plotOutput("cysc_scr_plot", height=260)),
          box(title="Cellular Stress: ROS & GSH", width=6,
            plotOutput("ros_gsh_plot", height=260))
        ),
        fluidRow(
          box(title="Tubular Cell Viability", width=12,
            plotOutput("tcv_plot", height=240))
        )
      ),
      ## ─── TAB 4: Clinical Endpoints ──────────────────────────
      tabItem("endpoints",
        fluidRow(
          box(title="eGFR Trajectory", width=6,
            status="danger", plotOutput("gfr_endpoint_plot", height=260)),
          box(title="AKI Stage over Time (KDIGO)", width=6,
            status="danger", plotOutput("aki_stage_plot", height=260))
        ),
        fluidRow(
          box(title="Urine Output (mL/kg/h)", width=6,
            plotOutput("uo_plot", height=260)),
          box(title="Inflammatory Cytokines", width=6,
            plotOutput("cytokine_plot", height=260))
        )
      ),
      ## ─── TAB 5: Scenario Comparison ─────────────────────────
      tabItem("compare",
        fluidRow(
          box(title="Comparative GFR — All Scenarios", width=12, status="primary",
            plotOutput("compare_gfr_plot", height=320))
        ),
        fluidRow(
          box(title="Comparative sCr", width=6,
            plotOutput("compare_scr_plot", height=280)),
          box(title="Comparative NGAL", width=6,
            plotOutput("compare_ngal_plot", height=280))
        ),
        fluidRow(
          box(title="Clinical Summary Table", width=12,
            DTOutput("compare_tbl"))
        )
      ),
      ## ─── TAB 6: AKI-to-CKD Transition ──────────────────────
      tabItem("fibrosis",
        fluidRow(
          box(title="TGF-β1 Dynamics", width=6, status="warning",
            plotOutput("tgfb_plot", height=260)),
          box(title="Myofibroblast & Fibrosis Index", width=6, status="warning",
            plotOutput("fibrosis_plot", height=260))
        ),
        fluidRow(
          box(title="Long-term GFR Recovery Trajectory", width=8,
            plotOutput("longterm_gfr_plot", height=260)),
          box(title="AKI-to-CKD Risk Summary", width=4,
            verbatimTextOutput("ckd_risk_text"))
        )
      )
    )
  )
)

## ─── SERVER ────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  sim_result <- reactive({
    params <- list(
      AKI_type = as.numeric(input$aki_type),
      isch_sev = input$isch_sev,
      tox_dose = input$tox_dose,
      lps_amt  = input$lps_amt,
      GFR0     = input$gfr0_val,
      SCR0     = input$scr0_val,
      CRRT_on  = as.numeric(input$crrt_on)
    )
    events <- NULL
    evlist <- list()
    if (input$use_fur)
      evlist[[length(evlist)+1]] <- ev(amt=40, ii=12, addl=23, cmt="FUR_GUT", time=2)
    if (input$use_ne)
      evlist[[length(evlist)+1]] <- ev(amt=input$ne_rate*60, rate=input$ne_rate*60,
                                       cmt="NE_CENT", time=0)
    if (input$use_nac)
      evlist[[length(evlist)+1]] <- ev(amt=10000, rate=1000, cmt="NAC_CENT", time=0)
    if (length(evlist) > 0) events <- do.call("+", evlist)
    run_aki(aki_mod, params, events,
            end = input$sim_days * 24, delta = 0.5)
  })

  ## All-scenarios for comparison
  all_scenarios <- reactive({
    pbase <- list(AKI_type=1, isch_sev=0.6, GFR0=100, SCR0=0.9, CRRT_on=0)
    s1 <- run_aki(aki_mod, pbase, NULL, 168) %>% mutate(Scenario="① IRI – No Tx")
    ev_fur <- ev(amt=40,ii=12,addl=13,cmt="FUR_GUT",time=2)
    s2 <- run_aki(aki_mod, pbase, ev_fur, 168) %>% mutate(Scenario="② IRI + Furosemide")
    ev_ne <- ev(amt=6,rate=6,cmt="NE_CENT",time=0)
    s3 <- run_aki(aki_mod, pbase, ev_ne, 168)  %>% mutate(Scenario="③ IRI + NE Support")
    ev_nac <- ev(amt=10000,rate=1000,cmt="NAC_CENT",time=-2)
    s4 <- run_aki(aki_mod, pbase, ev_nac,168)  %>% mutate(Scenario="④ IRI + NAC")
    psa <- modifyList(pbase, list(AKI_type=3, lps_amt=7, isch_sev=0.3))
    s5 <- run_aki(aki_mod, psa,  NULL, 168)    %>% mutate(Scenario="⑤ SA-AKI – No Tx")
    psa2 <- modifyList(psa, list(CRRT_on=1))
    s6 <- run_aki(aki_mod, psa2, ev_fur+ev_ne, 168) %>% mutate(Scenario="⑥ SA-AKI + Combo+CRRT")
    bind_rows(s1,s2,s3,s4,s5,s6)
  })

  ## ── Tab 1 ──
  output$pt_summary <- renderPrint({
    cat("AKI Type      :", switch(input$aki_type,
         "1"="Ischemia-Reperfusion","2"="Nephrotoxic","3"="Sepsis-Associated"),"\n")
    cat("Baseline GFR  :", input$gfr0_val, "mL/min/1.73m²\n")
    cat("Baseline sCr  :", input$scr0_val, "mg/dL\n")
    cat("CRRT          :", ifelse(input$crrt_on,"ON","OFF"), "\n")
    cat("Furosemide    :", ifelse(input$use_fur,"40mg IV q12h","none"), "\n")
    cat("NE infusion   :", ifelse(input$use_ne,paste(input$ne_rate,"mcg/kg/min"),"none"),"\n")
    cat("NAC           :", ifelse(input$use_nac,"IV prophylaxis","none"),"\n")
  })
  output$kdigo_tbl <- renderTable({
    data.frame(
      Stage = c("0 (None)","1","2","3"),
      sCr_criteria = c("< 1.5× baseline",
                       "1.5–1.9× OR +0.3 mg/dL","2.0–2.9×","≥ 3× OR dialysis"),
      UO_criteria  = c("≥ 0.5 mL/kg/h","< 0.5 mL/kg/h (>6h)",
                       "< 0.5 mL/kg/h (>12h)","< 0.3 mL/kg/h (>24h)"),
      Mortality    = c("Ref","~20%","~30%","~50%")
    )
  })
  output$patient_overview_plot <- renderPlot({
    d <- sim_result()
    par(mfrow=c(1,2), mar=c(4,4,2,1))
    plot(d$time/24, d$eGFR, type="l", lwd=2, col="#1565C0",
         xlab="Days", ylab="eGFR (mL/min/1.73m²)", main="GFR Trajectory")
    abline(h=60, lty=2, col="gray")
    plot(d$time/24, d$SCR, type="l", lwd=2, col="#B71C1C",
         xlab="Days", ylab="sCr (mg/dL)", main="Serum Creatinine")
    abline(h=c(input$scr0_val*1.5, input$scr0_val*2, input$scr0_val*3),
           lty=c(3,2,1), col=c("#EF5350","#E53935","#B71C1C"))
  })

  ## ── Tab 2 ──
  output$fur_pk_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(time, FUR_CENT)) +
      geom_line(color="#1565C0", size=1.2) +
      labs(x="Time (h)", y="Furosemide (mg/L)", title="Furosemide Central Conc.") +
      theme_bw()
  })
  output$ne_nac_pk_plot <- renderPlot({
    d <- sim_result() %>% select(time, NE_CENT, NAC_CENT) %>%
      pivot_longer(-time, names_to="Drug", values_to="Conc")
    ggplot(d, aes(time, Conc, color=Drug)) +
      geom_line(size=1.1) +
      labs(x="Time (h)", y="Plasma Conc.", title="NE & NAC Concentrations") +
      theme_bw()
  })
  output$nkcc2_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(time, E_NKCC2*100)) +
      geom_line(color="#0277BD", size=1.2) +
      labs(x="Time (h)", y="NKCC2 Inhibition (%)", title="Furosemide: NKCC2 Effect") +
      theme_bw()
  })
  output$diuresis_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(time, diuresis)) +
      geom_line(color="#0288D1", size=1.2) +
      geom_hline(yintercept=1, linetype="dashed", color="gray") +
      labs(x="Time (h)", y="Diuresis Fold-Increase", title="Furosemide Diuretic Response") +
      theme_bw()
  })

  ## ── Tab 3 ──
  output$ngal_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(time, NGAL)) +
      geom_line(color="#E65100", size=1.3) +
      geom_hline(yintercept=150, linetype="dashed", color="red", alpha=0.6) +
      annotate("text", x=max(d$time)*0.7, y=158, label="Diagnostic threshold ~150 ng/mL",
               color="red", size=3) +
      labs(x="Time (h)", y="NGAL (ng/mL)", title="NGAL — Tubular Injury Marker") +
      theme_bw()
  })
  output$kim1_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(time, KIM1)) +
      geom_line(color="#BF360C", size=1.3) +
      labs(x="Time (h)", y="KIM-1 (ng/mL)", title="KIM-1 — Brush Border Shedding") +
      theme_bw()
  })
  output$cysc_scr_plot <- renderPlot({
    d <- sim_result() %>%
      select(time, SCR, CysC) %>%
      mutate(SCR_norm=SCR/max(SCR), CysC_norm=CysC/max(CysC)) %>%
      pivot_longer(-time, names_to="Marker", values_to="Value")
    ggplot(filter(d, Marker %in% c("SCR_norm","CysC_norm")),
           aes(time, Value, color=Marker)) +
      geom_line(size=1.1) +
      scale_color_manual(values=c(SCR_norm="#C62828", CysC_norm="#1565C0"),
                         labels=c("Creatinine (norm)","Cystatin C (norm)")) +
      labs(x="Time (h)", y="Normalized Level", title="Cystatin C vs Creatinine",
           color="Marker") + theme_bw()
  })
  output$ros_gsh_plot <- renderPlot({
    d <- sim_result() %>%
      select(time, ROS, GSH) %>%
      pivot_longer(-time, names_to="Variable", values_to="Value")
    ggplot(d, aes(time, Value, color=Variable)) +
      geom_line(size=1.2) +
      scale_color_manual(values=c(ROS="#D32F2F", GSH="#388E3C")) +
      labs(x="Time (h)", y="Normalized Level", title="Oxidative Stress: ROS vs GSH",
           color="Variable") + theme_bw()
  })
  output$tcv_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(time, TCV)) +
      geom_line(color="#4A148C", size=1.4) +
      geom_hline(yintercept=0.25, linetype="dashed", color="#880E4F") +
      annotate("text", x=max(d$time)*0.6, y=0.27,
               label="Tubular atrophy risk threshold", color="#880E4F", size=3.5) +
      labs(x="Time (h)", y="TCV (0–1)", title="Tubular Cell Viability") +
      scale_y_continuous(limits=c(0,1.05)) +
      theme_bw()
  })

  ## ── Tab 4 ──
  output$gfr_endpoint_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(time/24, eGFR)) +
      geom_line(color="#1565C0", size=1.4) +
      geom_hline(yintercept=c(60, 30, 15), linetype="dashed",
                 color=c("orange","red","darkred")) +
      annotate("text", x=0.1, y=c(62, 32, 17),
               label=c("CKD G3","CKD G4","CKD G5"), color=c("orange","red","darkred"),
               hjust=0, size=3) +
      labs(x="Days", y="eGFR (mL/min/1.73m²)", title="eGFR Trajectory") +
      theme_bw()
  })
  output$aki_stage_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(time/24, AKI_stage)) +
      geom_step(color="#B71C1C", size=1.3) +
      scale_y_continuous(breaks=0:3, labels=c("None","Stage 1","Stage 2","Stage 3")) +
      labs(x="Days", y="AKI Stage (KDIGO)", title="AKI Stage over Time") +
      theme_bw()
  })
  output$uo_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(time/24, UO)) +
      geom_line(color="#00796B", size=1.2) +
      geom_hline(yintercept=c(0.5, 0.3), linetype=c("dashed","dotted"),
                 color=c("orange","red")) +
      annotate("text", x=0.2, y=c(0.52, 0.32), label=c("Oliguria","Severe Oliguria"),
               hjust=0, size=3, color=c("orange","red")) +
      labs(x="Days", y="UO (mL/kg/h)", title="Urine Output") +
      theme_bw()
  })
  output$cytokine_plot <- renderPlot({
    d <- sim_result() %>%
      select(time, IL6, TNFa) %>%
      pivot_longer(-time, names_to="Cytokine", values_to="Level")
    ggplot(d, aes(time/24, Level, color=Cytokine)) +
      geom_line(size=1.2) +
      scale_color_manual(values=c(IL6="#2E7D32", TNFa="#827717")) +
      labs(x="Days", y="Cytokine (pg/mL)", title="IL-6 & TNF-α Dynamics",
           color="Cytokine") + theme_bw()
  })

  ## ── Tab 5 ──
  output$compare_gfr_plot <- renderPlot({
    d <- all_scenarios()
    ggplot(d, aes(time/24, eGFR, color=Scenario)) +
      geom_line(size=1.1) +
      geom_hline(yintercept=60, linetype="dashed", color="gray40") +
      labs(x="Days", y="eGFR (mL/min/1.73m²)",
           title="All Scenarios: eGFR Comparison") +
      theme_bw() + theme(legend.position="right") +
      guides(color=guide_legend(ncol=1))
  })
  output$compare_scr_plot <- renderPlot({
    d <- all_scenarios()
    ggplot(d, aes(time/24, SCR, color=Scenario)) +
      geom_line(size=1) +
      labs(x="Days", y="sCr (mg/dL)", title="Serum Creatinine") +
      theme_bw() + theme(legend.position="bottom") +
      guides(color=guide_legend(ncol=2, byrow=TRUE))
  })
  output$compare_ngal_plot <- renderPlot({
    d <- all_scenarios()
    ggplot(d, aes(time/24, NGAL, color=Scenario)) +
      geom_line(size=1) +
      labs(x="Days", y="NGAL (ng/mL)", title="NGAL Biomarker") +
      theme_bw() + theme(legend.position="bottom") +
      guides(color=guide_legend(ncol=2, byrow=TRUE))
  })
  output$compare_tbl <- renderDT({
    d <- all_scenarios() %>%
      group_by(Scenario) %>%
      summarise(
        `Peak sCr (mg/dL)`  = round(max(SCR),2),
        `Max AKI Stage`     = max(AKI_stage),
        `Nadir GFR`         = round(min(eGFR),1),
        `Peak NGAL (ng/mL)` = round(max(NGAL),0),
        `Peak IL-6 (pg/mL)` = round(max(IL6),0),
        `Min TCV`           = round(min(TCV),3),
        `Final GFR`         = round(last(eGFR),1),
        .groups="drop"
      )
    datatable(d, options=list(pageLength=10, scrollX=TRUE), rownames=FALSE)
  })

  ## ── Tab 6 ──
  output$tgfb_plot <- renderPlot({
    d <- sim_result()
    ggplot(d, aes(time/24, TGFb)) +
      geom_line(color="#F57F17", size=1.3) +
      labs(x="Days", y="TGF-β1 (pg/mL)", title="TGF-β1: Pro-Fibrotic Driver") +
      theme_bw()
  })
  output$fibrosis_plot <- renderPlot({
    d <- sim_result()
    ggplot(d) +
      geom_line(aes(time/24, FIBROSIS*100, color="Fibrosis Index (%)"), size=1.2) +
      geom_line(aes(time/24, MYO*100, color="Myofibroblast (%)"), size=1, linetype="dashed") +
      scale_color_manual(values=c("Fibrosis Index (%)"="#D32F2F",
                                  "Myofibroblast (%)"="#EF6C00")) +
      labs(x="Days", y="Index (%)", title="Fibrosis & Myofibroblast Burden",
           color="Variable") + theme_bw()
  })
  output$longterm_gfr_plot <- renderPlot({
    params_lt <- list(
      AKI_type = as.numeric(input$aki_type),
      isch_sev = input$isch_sev, lps_amt=input$lps_amt, tox_dose=input$tox_dose,
      GFR0=input$gfr0_val, SCR0=input$scr0_val, CRRT_on=0
    )
    d <- run_aki(aki_mod, params_lt, NULL, end=720, delta=2)  # 30 days
    ggplot(d, aes(time/24, eGFR)) +
      geom_line(color="#1565C0", size=1.3) +
      geom_hline(yintercept=c(input$gfr0_val, 60), linetype=c("dotted","dashed"),
                 color=c("gray","orange")) +
      labs(x="Days", y="eGFR (mL/min/1.73m²)",
           title="Long-term GFR Recovery (30 days)") +
      theme_bw()
  })
  output$ckd_risk_text <- renderPrint({
    d <- sim_result()
    final_gfr  <- round(last(d$eGFR), 1)
    max_stage  <- max(d$AKI_stage)
    final_fibr <- round(last(d$FIBROSIS)*100, 2)
    final_tgfb <- round(last(d$TGFb), 2)
    cat("─── AKI-to-CKD Risk Assessment ───\n\n")
    cat("Final eGFR     :", final_gfr, "mL/min/1.73m²\n")
    cat("Max AKI Stage  :", max_stage, "\n")
    cat("Fibrosis Index :", final_fibr, "%\n")
    cat("TGF-β1 (final) :", final_tgfb, "pg/mL\n\n")
    ckd_risk <- if(max_stage >= 2 && final_gfr < 60) "HIGH"
                else if(max_stage >= 1 && final_fibr > 5) "MODERATE"
                else "LOW"
    cat("CKD Progression Risk:", ckd_risk, "\n\n")
    cat("KDIGO recommendation:\n")
    cat(" ≥ Stage 2 → 3-month follow-up\n")
    cat(" ≥ Stage 3 → nephrology referral\n")
    cat(" GFR < 60 at 3m → CKD confirmed\n")
  })
}

## ─── Launch ───────────────────────────────────────────────────────────────
shinyApp(ui, server)
