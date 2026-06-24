## ============================================================
## Multiple Sclerosis QSP — Interactive Shiny Dashboard
## Author  : Claude Code Routine (CCR)
## Date    : 2026-06-16
## ============================================================
## Tabs:
##   1. Patient Profile & Disease Phenotype
##   2. Drug PK — Concentration-Time Profiles
##   3. Immunology — Peripheral & CNS Immune Dynamics
##   4. Disease Progression — MRI & BBB
##   5. Clinical Endpoints — EDSS, ARR, NfL
##   6. Scenario Comparison — Multiple Drugs
##   7. Biomarkers & Safety Monitoring
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)
library(DT)
library(scales)
library(patchwork)

## ── Embed mrgsolve model code ─────────────────────────────────────────────────

ms_model_code <- '
$PARAM @annotated
dose_ifnb:0:IFN-beta active
dose_nat:0:Natalizumab active
dose_ocre:0:Ocrelizumab active
dose_sip:0:Siponimod active
dose_dmf:0:DMF active
dose_clad:0:Cladribine active

ka_ifnb:0.30:/h  CL_ifnb:5.0:L/h  V1_ifnb:12.0:L  Q_ifnb:1.5:L/h  V2_ifnb:20.0:L  F_ifnb:0.40:
CL_nat:0.014:L/h  V1_nat:3.5:L  Q_nat:0.50:L/h  V2_nat:6.0:L
kon_nat:0.013:/h  koff_nat:0.003:/h  kint_nat:0.005:/h  Rtot_nat:30.0:nM  ksyn_R:0.15:nM/h  kdeg_R:0.005:/h
CL_ocre:0.012:L/h  V1_ocre:3.0:L  Q_ocre:0.40:L/h  V2_ocre:5.0:L  EC50_ocre:1.5:ug/mL  Emax_ocre:0.97:
ka_sip:0.50:/h  CL_sip:2.2:L/h  V1_sip:30.0:L  Q_sip:3.0:L/h  V2_sip:60.0:L  EC50_sip:8.0:ng/mL  Emax_sip:0.75:
ka_dmf:1.20:/h  CL_dmf:25.0:L/h  V1_dmf:50.0:L  EC50_dmf:0.40:ug/mL  EC50_nrf2:0.20:ug/mL
ka_clad:0.40:/h  CL_clad:15.0:L/h  V1_clad:40.0:L  EC50_clad:5.0:ng/mL  Emax_clad:0.90:

kin_Th1:0.10:/h  kout_Th1:0.10:/h
kin_Th17:0.06:/h  kout_Th17:0.06:/h
kin_Treg:0.04:/h  kout_Treg:0.04:/h  Treg_ss:1.0:  k_Treg_inh:0.30:
kin_B:0.08:/h  kout_B:0.08:/h
kin_BBB:0.02:/h  k_Th17_BBB:0.02:  k_Th1_BBB:0.01:  BBB_ss:0.85:
kin_cTh1:0.008:/h  kin_cTh17:0.005:/h  kout_cTh1:0.015:/h  kout_cTh17:0.015:/h  k_nat_inf:0.90:
kin_micro:0.004:/h  kout_micro:0.006:/h  k_IL10_micro:0.40:
kin_oligo:0.005:/h  kout_oligo:0.005:/h  k_micro_oligo:0.04:  k_ROS_oligo:0.02:
kin_OPC:0.008:/h  kout_OPC:0.003:/h  kdiff_OPC:0.002:/h  k_LINGO_inh:0.50:
kin_myelin:0.003:/h  kdmyelin:0.02:  k_remyelin:0.005:/h
k_axon_loss:0.010:  k_axon_rep:0.002:/h
kin_NfL:0.0020:/h  k_axon_NfL:0.040:  kelim_NfL:0.030:/h
kin_GFAP:0.002:/h  k_astro_GFAP:0.020:  kelim_GFAP:0.025:/h
k_EDSS_acc:0.0005:/h  k_EDSS_max:10.0:
k_relapse_base:1.5:  k_relapse_myelin:3.0:

$CMT
A1_ifnb C1_ifnb C2_ifnb
C1_nat C2_nat RC_nat R_nat
C1_ocre C2_ocre
C1_sip C2_sip
C1_dmf
C1_clad
Th1 Th17 Treg Bcell
BBB cTh1 cTh17 Micro
Oligo OPC Myelin Axon
NfL GFAP

$MAIN
double E_ifnb  = (dose_ifnb>0) ? C1_ifnb/(C1_ifnb+500.0) : 0;
double I_VCAM  = E_ifnb*0.50;
double I_TNF   = E_ifnb*0.40;
double E_Treg_ifnb = E_ifnb*0.30;
double occ_nat_val = (dose_nat>0) ? RC_nat/(RC_nat+R_nat+1e-6) : 0;
double I_BBB_nat = occ_nat_val*k_nat_inf;
double E_ocre  = (dose_ocre>0) ? Emax_ocre*C1_ocre/(C1_ocre+EC50_ocre) : 0;
double E_sip   = (dose_sip>0)  ? Emax_sip*C1_sip/(C1_sip+EC50_sip) : 0;
double I_NFKB  = (dose_dmf>0)  ? C1_dmf/(C1_dmf+EC50_dmf)*0.45 : 0;
double E_NRF2  = (dose_dmf>0)  ? C1_dmf/(C1_dmf+EC50_nrf2)*0.55 : 0;
double I_Th_dmf = I_NFKB*0.50;
double I_clad  = (dose_clad>0) ? Emax_clad*C1_clad/(C1_clad+EC50_clad) : 0;
double I_lymph = 1.0-(1.0-E_sip)*(1.0-I_clad)*(1.0-I_Th_dmf);
double I_micro_Treg = k_IL10_micro*(Treg/Treg_ss);
double ROS_prot = 1.0-E_NRF2;
double EDSS_v = k_EDSS_max*(1.0-Axon)*0.8+k_EDSS_max*(1.0-Myelin)*0.2;
EDSS_v = (EDSS_v>10.0)?10.0:(EDSS_v<0.0?0.0:EDSS_v);
double myel_lr = kdmyelin*Micro*Myelin;
double ARR_v = k_relapse_base*myel_lr/(kdmyelin*0.5)
              *(1.0-occ_nat_val*0.68)*(1.0-E_sip*0.55)
              *(1.0-E_ocre*0.47)*(1.0-I_clad*0.58)
              *(1.0-E_ifnb*0.34)*(1.0-I_NFKB*0.49);

$ODE
dxdt_A1_ifnb = -ka_ifnb*A1_ifnb;
dxdt_C1_ifnb = (dose_ifnb>0?F_ifnb*ka_ifnb*A1_ifnb:0.0)-(CL_ifnb/V1_ifnb+Q_ifnb/V1_ifnb)*C1_ifnb+(Q_ifnb/V2_ifnb)*C2_ifnb;
dxdt_C2_ifnb = (Q_ifnb/V1_ifnb)*C1_ifnb-(Q_ifnb/V2_ifnb)*C2_ifnb;
dxdt_C1_nat  = -(CL_nat/V1_nat+Q_nat/V1_nat)*C1_nat+(Q_nat/V2_nat)*C2_nat-kon_nat*C1_nat*R_nat+koff_nat*RC_nat;
dxdt_C2_nat  = (Q_nat/V1_nat)*C1_nat-(Q_nat/V2_nat)*C2_nat;
dxdt_RC_nat  = kon_nat*C1_nat*R_nat-(koff_nat+kint_nat)*RC_nat;
dxdt_R_nat   = ksyn_R-kdeg_R*R_nat-kon_nat*C1_nat*R_nat+koff_nat*RC_nat;
dxdt_C1_ocre = -(CL_ocre/V1_ocre+Q_ocre/V1_ocre)*C1_ocre+(Q_ocre/V2_ocre)*C2_ocre;
dxdt_C2_ocre = (Q_ocre/V1_ocre)*C1_ocre-(Q_ocre/V2_ocre)*C2_ocre;
dxdt_C1_sip  = ka_sip*C2_sip-(CL_sip/V1_sip+Q_sip/V1_sip)*C1_sip+(Q_sip/V2_sip)*C2_sip;
dxdt_C2_sip  = -(ka_sip+Q_sip/V2_sip)*C2_sip+(Q_sip/V1_sip)*C1_sip;
dxdt_C1_dmf  = -(CL_dmf/V1_dmf)*C1_dmf;
dxdt_C1_clad = -(CL_clad/V1_clad)*C1_clad;

dxdt_Th1  = kin_Th1*(1.0-I_lymph)*(1.0-I_TNF*0.5)+kin_Th1*E_Treg_ifnb*0.2-kout_Th1*Th1-kout_Th1*I_clad*Th1;
dxdt_Th17 = kin_Th17*(1.0-I_lymph)*(1.0/(1.0+k_Treg_inh*Treg/Treg_ss))-kout_Th17*Th17;
dxdt_Treg = kin_Treg*(1.0+E_Treg_ifnb)-kout_Treg*Treg;
dxdt_Bcell = kin_B*(1.0-E_ocre)*(1.0-I_lymph)-kout_B*Bcell;
dxdt_BBB   = kin_BBB*(1.0-BBB)-BBB*(k_Th17_BBB*Th17+k_Th1_BBB*Th1)*(1.0-I_VCAM)+kin_BBB*I_VCAM;
double BBB_op = 1.0-BBB;
dxdt_cTh1  = kin_cTh1*Th1*(1.0+BBB_op*2.0)*(1.0-I_BBB_nat)-kout_cTh1*cTh1;
dxdt_cTh17 = kin_cTh17*Th17*(1.0+BBB_op*3.0)*(1.0-I_BBB_nat)-kout_cTh17*cTh17;
dxdt_Micro = kin_micro*(cTh1+2.0*cTh17+0.5*(1.0-Bcell))*(1.0-I_micro_Treg)-kout_micro*Micro;
dxdt_Oligo = kin_oligo-kout_oligo*Oligo-Oligo*(k_micro_oligo*Micro+k_ROS_oligo*(1.0-E_NRF2)*Micro);
dxdt_OPC   = kin_OPC-kout_OPC*OPC-kdiff_OPC*OPC*(1.0-k_LINGO_inh*Micro);
dxdt_Myelin = k_remyelin*OPC*(1.0-Myelin)-kdmyelin*Micro*Myelin*ROS_prot+kin_myelin*(1.0-Myelin);
dxdt_Axon   = k_axon_rep*Myelin-k_axon_loss*(1.0-Myelin)-k_axon_loss*0.5*Micro*Axon;
dxdt_NfL    = kin_NfL+k_axon_NfL*(1.0-Axon)*(1.0+Micro)-kelim_NfL*NfL;
dxdt_GFAP   = kin_GFAP+k_astro_GFAP*(cTh17+Micro*0.5)*(1.0-BBB)-kelim_GFAP*GFAP;

$TABLE
capture EDSS=EDSS_v; capture ARR=ARR_v;
capture T2lesion=(1.0-Myelin)*15000.0;
capture occ_nat=occ_nat_val*100.0;
capture E_ocre_pct=E_ocre*100.0;
capture E_sip_pct=E_sip*100.0;

$INIT
A1_ifnb=0,C1_ifnb=0,C2_ifnb=0,
C1_nat=0,C2_nat=0,RC_nat=0,R_nat=30.0,
C1_ocre=0,C2_ocre=0,
C1_sip=0,C2_sip=0,
C1_dmf=0,C1_clad=0,
Th1=1.5,Th17=1.8,Treg=0.7,Bcell=1.2,
BBB=0.80,cTh1=0.30,cTh17=0.25,
Micro=0.50,Oligo=0.80,OPC=0.90,
Myelin=0.75,Axon=0.85,
NfL=12.0,GFAP=180.0
'

ms <- mrgsolve::mcode("ms_shiny", ms_model_code)

## ── Drug color scheme ─────────────────────────────────────────────────────────

drug_colors <- c(
  none = "#7F8C8D", ifnb = "#3498DB", nat = "#E74C3C",
  ocre = "#9B59B6", sip  = "#F39C12", dmf = "#27AE60", clad = "#E67E22"
)
drug_labels <- c(
  none = "Untreated (RRMS)", ifnb = "IFN-β (Avonex)", nat = "Natalizumab",
  ocre = "Ocrelizumab", sip  = "Siponimod", dmf  = "Dimethyl Fumarate",
  clad = "Cladribine"
)

## ── Simulation runner ─────────────────────────────────────────────────────────

run_sim <- function(drug, duration_wk = 104, weight_kg = 70,
                    th1_init = 1.5, th17_init = 1.8, myelin_init = 0.75) {
  init_vals <- c(Th1 = th1_init, Th17 = th17_init, Myelin = myelin_init,
                 Treg = 0.7, Bcell = 1.2, BBB = 0.80,
                 cTh1 = 0.30, cTh17 = 0.25, Micro = 0.50,
                 Oligo = 0.80, OPC = 0.90, Axon = 0.85,
                 NfL = 12.0, GFAP = 180.0)

  obs_t <- sort(unique(c(0, seq(0, 12*7*24, by=7*24),
                          seq(12*7*24, duration_wk*7*24, by=4*7*24))))

  dose_events <- switch(drug,
    "ifnb" = ev(cmt="A1_ifnb", amt=30, ii=7*24,
                addl=round(duration_wk*7*24/(7*24))-1, time=0),
    "nat"  = ev(cmt="C1_nat",  amt=300000/3.5, ii=4*7*24,
                addl=round(duration_wk*7*24/(4*7*24))-1, time=0),
    "ocre" = {
      df <- data.frame(
        cmt  = c("C1_ocre","C1_ocre","C1_ocre"),
        amt  = c(300000/3.0, 300000/3.0, 600000/3.0),
        time = c(0, 14*24, 24*7*24), evid = 1)
      as.ev(df)
    },
    "sip"  = ev(cmt="C2_sip",  amt=2000*0.84/30.0, ii=24,
                addl=round(duration_wk*7*24/24)-1, time=0),
    "dmf"  = ev(cmt="C1_dmf",  amt=240000/50.0*0.25, ii=12,
                addl=round(duration_wk*7*24/12)-1, time=0),
    "clad" = {
      df <- data.frame(
        cmt = "C1_clad",
        amt = rep(3.5*weight_kg/8*1000/40.0, 8),
        time = c(0,24,48,72,96,120,144,168), evid=1)
      as.ev(df)
    },
    NULL  # none
  )

  p_override <- switch(drug,
    "ifnb" = list(dose_ifnb=1), "nat"  = list(dose_nat=1),
    "ocre" = list(dose_ocre=1), "sip"  = list(dose_sip=1),
    "dmf"  = list(dose_dmf=1),  "clad" = list(dose_clad=1),
    list()
  )

  sim <- ms %>%
    init(.as_list(init_vals)) %>%
    param(.as_list(p_override))

  if (!is.null(dose_events)) {
    sim <- sim %>% ev(dose_events)
  }

  sim %>%
    mrgsim(end = duration_wk*7*24, delta = 24, obsonly = TRUE) %>%
    as_tibble() %>%
    mutate(drug = drug, time_wk = time/(7*24), time_yr = time/(365.25*24))
}

.as_list <- function(v) as.list(v)

## ═══════════════════════════════════════════════════════════════════════════════
## UI
## ═══════════════════════════════════════════════════════════════════════════════

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "MS QSP Dashboard", titleWidth = 280),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "patient",  icon = icon("user-md")),
      menuItem("Drug PK",             tabName = "pk",       icon = icon("flask")),
      menuItem("Immunology",          tabName = "immuno",   icon = icon("microscope")),
      menuItem("Disease Progression", tabName = "mri",      icon = icon("brain")),
      menuItem("Clinical Endpoints",  tabName = "clinical", icon = icon("chart-line")),
      menuItem("Scenario Comparison", tabName = "scenario", icon = icon("balance-scale")),
      menuItem("Biomarkers & Safety", tabName = "biomarker",icon = icon("heartbeat"))
    ),
    hr(),
    h5("Simulation Settings", style = "color:#BDC3C7; padding-left:15px;"),
    sliderInput("duration_wk", "Duration (weeks):", min=24, max=208, value=104, step=4),
    sliderInput("weight_kg",   "Body weight (kg):", min=40, max=140, value=70, step=5),
    hr(),
    h5("RRMS Disease Severity", style = "color:#BDC3C7; padding-left:15px;"),
    sliderInput("th1_init",    "Th1 (baseline 1.5):",   min=1.0, max=3.0, value=1.5, step=0.1),
    sliderInput("th17_init",   "Th17 (baseline 1.8):",  min=1.0, max=3.5, value=1.8, step=0.1),
    sliderInput("myelin_init", "Myelin (baseline 0.75):",min=0.5, max=1.0, value=0.75, step=0.05)
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .box-header { background-color: #2C3E50 !important; color: white; }
      .nav-tabs-custom>.tab-content { padding: 10px; }
      .info-box { min-height: 70px; }
    "))),

    tabItems(

      ## ── Tab 1: Patient Profile ───────────────────────────────────────────────
      tabItem(tabName = "patient",
        fluidRow(
          box(title = "Patient & Disease Characteristics", width = 6, status = "primary",
            radioButtons("ms_type", "MS Phenotype:",
              choices = c("RRMS (Relapsing-Remitting)" = "rrms",
                          "SPMS (Secondary Progressive)" = "spms",
                          "PPMS (Primary Progressive)"   = "ppms"),
              selected = "rrms"),
            selectInput("drug_pt", "Select Drug Monotherapy:",
              choices = c("Untreated" = "none", "IFN-β (Avonex)" = "ifnb",
                          "Natalizumab" = "nat", "Ocrelizumab" = "ocre",
                          "Siponimod" = "sip", "Dimethyl Fumarate" = "dmf",
                          "Cladribine" = "clad"),
              selected = "none"),
            numericInput("edss_init", "Baseline EDSS:", value = 2.5, min=0, max=9, step=0.5),
            numericInput("age_init",  "Patient Age (years):", value = 35, min=18, max=65),
            selectInput("sex_pt", "Sex:", choices = c("Female"="F","Male"="M"), selected="F"),
            numericInput("disease_dur", "Disease Duration (years):", value=3, min=0, max=20),
            actionButton("run_pt", "Run Simulation", class = "btn-primary btn-lg btn-block")
          ),
          box(title = "Baseline Disease State", width = 6, status = "info",
            plotOutput("pt_radar", height = "300px"),
            hr(),
            fluidRow(
              infoBoxOutput("ib_edss",   width = 4),
              infoBoxOutput("ib_arr",    width = 4),
              infoBoxOutput("ib_nfl",    width = 4)
            )
          )
        ),
        fluidRow(
          box(title = "Mechanistic Overview", width = 12, status = "warning",
            HTML('<div style="text-align:center">
              <img src="https://raw.githubusercontent.com/pipetcpt/qsp/main/multiple-sclerosis/ms_qsp.svg"
                   style="max-width:100%; height:auto" alt="MS QSP Map" />
              <p style="color:#7F8C8D; font-size:11px; margin-top:5px;">
                MS QSP Mechanistic Map — 12 subgraph clusters, 160+ nodes
              </p>
            </div>')
          )
        )
      ),

      ## ── Tab 2: Drug PK ───────────────────────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Drug Selection & PK Settings", width = 3, status = "primary",
            checkboxGroupInput("drugs_pk", "Select Drugs:",
              choices = c("IFN-β" = "ifnb", "Natalizumab" = "nat",
                          "Ocrelizumab" = "ocre", "Siponimod" = "sip",
                          "DMF" = "dmf", "Cladribine" = "clad"),
              selected = c("nat", "sip")),
            selectInput("pk_var", "PK Variable:",
              choices = c("Natalizumab [µg/mL]" = "C1_nat",
                          "Siponimod [ng/mL]"   = "C1_sip",
                          "IFN-β central [IU/mL]" = "C1_ifnb",
                          "Ocrelizumab [µg/mL]" = "C1_ocre",
                          "MMF [µg/mL]"         = "C1_dmf",
                          "Cladribine [ng/mL]"  = "C1_clad"),
              selected = "C1_nat"),
            actionButton("run_pk", "Update PK Plots", class = "btn-info btn-block")
          ),
          box(title = "Concentration-Time Profile", width = 9, status = "info",
            plotlyOutput("pk_plot", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Natalizumab TMDD — Receptor Occupancy", width = 6, status = "warning",
            plotlyOutput("tmdd_plot", height = "300px")
          ),
          box(title = "Ocrelizumab — B Cell Depletion Kinetics", width = 6, status = "warning",
            plotlyOutput("bcell_plot", height = "300px")
          )
        )
      ),

      ## ── Tab 3: Immunology ────────────────────────────────────────────────────
      tabItem(tabName = "immuno",
        fluidRow(
          box(title = "Drug for Immunology Analysis", width = 3, status = "primary",
            selectInput("drug_imm", "Drug:",
              choices = setNames(names(drug_labels), drug_labels), selected = "nat"),
            checkboxGroupInput("imm_vars", "Variables to Plot:",
              choices = c("Th1" = "Th1", "Th17" = "Th17", "Treg" = "Treg",
                          "B cells" = "Bcell", "CNS Th1" = "cTh1",
                          "CNS Th17" = "cTh17", "Microglia" = "Micro"),
              selected = c("Th1","Th17","Treg","cTh1","Micro")),
            actionButton("run_imm", "Update", class = "btn-info btn-block")
          ),
          box(title = "Peripheral Immune Cell Dynamics", width = 9, status = "info",
            plotlyOutput("imm_plot", height = "420px")
          )
        ),
        fluidRow(
          box(title = "CNS Infiltration & Microglia", width = 6, status = "warning",
            plotlyOutput("cns_plot", height = "300px")
          ),
          box(title = "BBB Integrity", width = 6, status = "warning",
            plotlyOutput("bbb_plot", height = "300px")
          )
        )
      ),

      ## ── Tab 4: Disease Progression ───────────────────────────────────────────
      tabItem(tabName = "mri",
        fluidRow(
          box(title = "Drug for MRI Analysis", width = 3, status = "primary",
            selectInput("drug_mri", "Drug:",
              choices = setNames(names(drug_labels), drug_labels), selected = "nat"),
            actionButton("run_mri", "Update", class = "btn-info btn-block"),
            hr(),
            h5("MRI Reference Ranges:"),
            tags$ul(
              tags$li("T2 volume typical RRMS: 5,000–30,000 mm³"),
              tags$li("Myelin integrity >0.9 = near-normal"),
              tags$li("Axonal integrity >0.9 = low disability")
            )
          ),
          box(title = "Demyelination & Remyelination Dynamics", width = 9, status = "info",
            plotlyOutput("myelin_plot", height = "420px")
          )
        ),
        fluidRow(
          box(title = "T2 Lesion Volume Trajectory", width = 6, status = "warning",
            plotlyOutput("t2_plot", height = "300px")
          ),
          box(title = "Axonal Integrity Over Time", width = 6, status = "warning",
            plotlyOutput("axon_plot", height = "300px")
          )
        )
      ),

      ## ── Tab 5: Clinical Endpoints ────────────────────────────────────────────
      tabItem(tabName = "clinical",
        fluidRow(
          box(title = "Drug Selection", width = 3, status = "primary",
            selectInput("drug_clin", "Drug:",
              choices = setNames(names(drug_labels), drug_labels), selected = "nat"),
            actionButton("run_clin", "Update", class = "btn-info btn-block"),
            hr(),
            h5("EDSS Scale Guide:"),
            tags$table(class="table table-condensed", style="font-size:11px",
              tags$tr(tags$th("EDSS"), tags$th("Description")),
              tags$tr(tags$td("0"),    tags$td("Normal exam")),
              tags$tr(tags$td("1-1.5"),tags$td("No disability")),
              tags$tr(tags$td("2-2.5"),tags$td("Minimal disability")),
              tags$tr(tags$td("3-4.5"),tags$td("Moderate disability")),
              tags$tr(tags$td("5-6"),  tags$td("Walking aid needed")),
              tags$tr(tags$td("6.5-8"),tags$td("Bilateral walking aid")),
              tags$tr(tags$td("9-10"), tags$td("Bedbound / death"))
            )
          ),
          box(title = "EDSS Trajectory", width = 9, status = "info",
            plotlyOutput("edss_plot", height = "420px")
          )
        ),
        fluidRow(
          box(title = "Annualized Relapse Rate", width = 6, status = "warning",
            plotlyOutput("arr_plot", height = "300px")
          ),
          box(title = "NEDA-3 Component Monitor", width = 6, status = "warning",
            plotlyOutput("neda_plot", height = "300px")
          )
        )
      ),

      ## ── Tab 6: Scenario Comparison ───────────────────────────────────────────
      tabItem(tabName = "scenario",
        fluidRow(
          box(title = "Select Drugs to Compare", width = 3, status = "primary",
            checkboxGroupInput("drugs_comp", "Drugs:",
              choices = setNames(names(drug_labels), drug_labels),
              selected = c("none","nat","ocre","sip","dmf")),
            selectInput("comp_endpoint", "Endpoint:",
              choices = c("EDSS" = "EDSS", "ARR" = "ARR",
                          "T2 Lesion (mm³)" = "T2lesion",
                          "Myelin Integrity" = "Myelin",
                          "Axonal Integrity" = "Axon",
                          "Serum NfL (pg/mL)" = "NfL"),
              selected = "EDSS"),
            actionButton("run_comp", "Run All Scenarios",
                         class = "btn-danger btn-lg btn-block")
          ),
          box(title = "Multi-Drug Comparison Plot", width = 9, status = "danger",
            plotlyOutput("comp_plot", height = "420px")
          )
        ),
        fluidRow(
          box(title = "2-Year Summary Table", width = 12, status = "info",
            DTOutput("summary_table")
          )
        )
      ),

      ## ── Tab 7: Biomarkers & Safety ───────────────────────────────────────────
      tabItem(tabName = "biomarker",
        fluidRow(
          box(title = "Drug Selection", width = 3, status = "primary",
            selectInput("drug_bm", "Drug:",
              choices = setNames(names(drug_labels), drug_labels), selected = "ocre"),
            actionButton("run_bm", "Update", class = "btn-info btn-block"),
            hr(),
            h5("Biomarker Reference Values:"),
            tags$ul(
              tags$li("NfL healthy <10 pg/mL"),
              tags$li("NfL active RRMS: 15-40 pg/mL"),
              tags$li("GFAP healthy <150 pg/mL"),
              tags$li("GFAP active MS: 200-600 pg/mL")
            ),
            hr(),
            h5("Safety Monitoring Thresholds:"),
            tags$ul(
              tags$li(HTML("Siponimod: lymphocyte <0.2×10⁹/L → hold")),
              tags$li("Natalizumab: JCV index >1.5 + >24m → high PML risk"),
              tags$li("Cladribine: CD4 <200/µL → withhold")
            )
          ),
          box(title = "NfL & GFAP Biomarker Trajectories", width = 9, status = "info",
            plotlyOutput("bm_plot", height = "420px")
          )
        ),
        fluidRow(
          box(title = "Lymphocyte Count & Safety (Siponimod)", width = 6, status = "warning",
            plotlyOutput("lymph_plot", height = "300px")
          ),
          box(title = "Receptor Occupancy / Bioavailability", width = 6, status = "warning",
            plotlyOutput("occ_plot", height = "300px")
          )
        )
      )
    )
  )
)

## ═══════════════════════════════════════════════════════════════════════════════
## SERVER
## ═══════════════════════════════════════════════════════════════════════════════

server <- function(input, output, session) {

  ## Reactive simulation runners
  sim_pt <- eventReactive(input$run_pt, {
    withProgress(message = "Running simulation...", {
      run_sim(input$drug_pt, input$duration_wk, input$weight_kg,
              input$th1_init, input$th17_init, input$myelin_init)
    })
  })

  sim_imm <- eventReactive(input$run_imm, {
    run_sim(input$drug_imm, input$duration_wk, input$weight_kg,
            input$th1_init, input$th17_init, input$myelin_init)
  })

  sim_mri <- eventReactive(input$run_mri, {
    run_sim(input$drug_mri, input$duration_wk, input$weight_kg,
            input$th1_init, input$th17_init, input$myelin_init)
  })

  sim_clin <- eventReactive(input$run_clin, {
    run_sim(input$drug_clin, input$duration_wk, input$weight_kg,
            input$th1_init, input$th17_init, input$myelin_init)
  })

  sim_bm <- eventReactive(input$run_bm, {
    run_sim(input$drug_bm, input$duration_wk, input$weight_kg,
            input$th1_init, input$th17_init, input$myelin_init)
  })

  sim_comp <- eventReactive(input$run_comp, {
    withProgress(message = "Simulating all scenarios...", {
      bind_rows(lapply(input$drugs_comp, function(d) {
        setProgress(detail = paste("Running:", drug_labels[d]))
        run_sim(d, input$duration_wk, input$weight_kg,
                input$th1_init, input$th17_init, input$myelin_init)
      }))
    })
  })

  ## Shared PK simulation
  sim_pk <- eventReactive(input$run_pk, {
    bind_rows(lapply(input$drugs_pk, function(d) {
      run_sim(d, min(input$duration_wk, 52), input$weight_kg,
              input$th1_init, input$th17_init, input$myelin_init)
    }))
  })

  ## ── InfoBoxes ────────────────────────────────────────────────────────────────
  output$ib_edss <- renderInfoBox({
    infoBox("Baseline EDSS", input$edss_init,
            icon = icon("walking"), color = "light-blue", fill = TRUE)
  })
  output$ib_arr <- renderInfoBox({
    infoBox("Expected ARR", "1.5/yr",
            subtitle = "Untreated RRMS",
            icon = icon("redo"), color = "orange", fill = TRUE)
  })
  output$ib_nfl <- renderInfoBox({
    infoBox("NfL", "12 pg/mL",
            subtitle = "At baseline",
            icon = icon("dna"), color = "red", fill = TRUE)
  })

  ## ── Patient radar placeholder ─────────────────────────────────────────────
  output$pt_radar <- renderPlot({
    vals <- c(Th1 = input$th1_init/3, Th17 = input$th17_init/3.5,
              `1-Treg` = 1-0.7, `1-Myelin` = 1-input$myelin_init,
              `1-Axon` = 0.15, Micro = 0.5)
    df <- data.frame(metric = names(vals), value = unname(vals))
    ggplot(df, aes(x = metric, y = value, fill = metric)) +
      geom_col(show.legend = FALSE) +
      scale_fill_manual(values = c("#E74C3C","#F39C12","#9B59B6",
                                    "#E67E22","#C0392B","#8E44AD")) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(title = "Disease Severity Profile",
           y = "Normalized Severity", x = "") +
      theme_bw() + coord_flip()
  })

  ## ── PK plots ─────────────────────────────────────────────────────────────────
  output$pk_plot <- renderPlotly({
    req(sim_pk())
    d <- sim_pk()
    yvar <- input$pk_var
    if (!(yvar %in% names(d))) {
      yvar <- "C1_nat"
    }
    p <- ggplot(d, aes(x = time_wk, y = .data[[yvar]],
                        color = drug, linetype = drug)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = drug_colors, labels = drug_labels,
                         name = "Drug") +
      labs(title = paste("PK:", yvar), x = "Time (weeks)", y = yvar) +
      theme_bw()
    ggplotly(p)
  })

  output$tmdd_plot <- renderPlotly({
    req(sim_pk())
    d <- sim_pk() %>% filter(drug == "nat")
    if (nrow(d) == 0) return(plotly_empty())
    p <- ggplot(d, aes(x = time_wk, y = occ_nat)) +
      geom_line(color = "#E74C3C", linewidth = 1.3) +
      geom_hline(yintercept = 80, linetype = 2, color = "grey40") +
      labs(title = "Natalizumab: α4-integrin Receptor Occupancy",
           x = "Weeks", y = "% Receptor Saturation") +
      ylim(0, 100) + theme_bw()
    ggplotly(p)
  })

  output$bcell_plot <- renderPlotly({
    req(sim_pk())
    d <- sim_pk() %>% filter(drug == "ocre")
    if (nrow(d) == 0) return(plotly_empty())
    p <- ggplot(d, aes(x = time_wk, y = Bcell)) +
      geom_line(color = "#9B59B6", linewidth = 1.3) +
      geom_hline(yintercept = 0.05, linetype = 2, color = "grey40") +
      labs(title = "Ocrelizumab: B Cell Depletion Kinetics",
           x = "Weeks", y = "B cells (normalized, 1=baseline)") +
      ylim(0, 1.5) + theme_bw()
    ggplotly(p)
  })

  ## ── Immunology plots ──────────────────────────────────────────────────────────
  output$imm_plot <- renderPlotly({
    req(sim_imm())
    d <- sim_imm()
    vars <- intersect(input$imm_vars, names(d))
    if (length(vars) == 0) vars <- c("Th1", "Th17", "Treg")
    d_long <- d %>%
      select(time_wk, drug, all_of(vars)) %>%
      pivot_longer(all_of(vars), names_to = "variable", values_to = "value")
    p <- ggplot(d_long, aes(x = time_wk, y = value,
                             color = variable, linetype = variable)) +
      geom_line(linewidth = 1.1) +
      labs(title = paste("Immune Dynamics —", drug_labels[input$drug_imm]),
           x = "Weeks", y = "Normalized cell level") +
      theme_bw()
    ggplotly(p)
  })

  output$cns_plot <- renderPlotly({
    req(sim_imm())
    d <- sim_imm() %>%
      select(time_wk, cTh1, cTh17, Micro) %>%
      pivot_longer(c(cTh1, cTh17, Micro), names_to = "cell", values_to = "val")
    p <- ggplot(d, aes(x = time_wk, y = val, color = cell)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c(cTh1="#E74C3C", cTh17="#F39C12", Micro="#8E44AD")) +
      labs(title = "CNS Immune Infiltration", x = "Weeks", y = "Level") +
      theme_bw()
    ggplotly(p)
  })

  output$bbb_plot <- renderPlotly({
    req(sim_imm())
    p <- ggplot(sim_imm(), aes(x = time_wk, y = BBB)) +
      geom_line(color = "#27AE60", linewidth = 1.3) +
      geom_hline(yintercept = 0.95, linetype = 2, color = "grey40") +
      ylim(0, 1) +
      labs(title = "BBB Integrity", x = "Weeks",
           y = "Integrity (0=disrupted, 1=intact)") +
      theme_bw()
    ggplotly(p)
  })

  ## ── MRI / Disease plots ───────────────────────────────────────────────────────
  output$myelin_plot <- renderPlotly({
    req(sim_mri())
    d <- sim_mri() %>%
      select(time_wk, Myelin, OPC, Oligo) %>%
      pivot_longer(c(Myelin, OPC, Oligo), names_to = "var", values_to = "val")
    p <- ggplot(d, aes(x = time_wk, y = val, color = var)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c(Myelin="#F39C12", OPC="#3498DB", Oligo="#27AE60")) +
      ylim(0, 1.1) +
      labs(title = "Myelination Biology", x = "Weeks",
           y = "Normalized level (1=healthy)") +
      theme_bw()
    ggplotly(p)
  })

  output$t2_plot <- renderPlotly({
    req(sim_mri())
    p <- ggplot(sim_mri(), aes(x = time_wk, y = T2lesion)) +
      geom_line(color = "#E74C3C", linewidth = 1.3) +
      labs(title = "T2 Lesion Volume", x = "Weeks", y = "Volume (mm³)") +
      theme_bw()
    ggplotly(p)
  })

  output$axon_plot <- renderPlotly({
    req(sim_mri())
    p <- ggplot(sim_mri(), aes(x = time_wk, y = Axon)) +
      geom_line(color = "#8E44AD", linewidth = 1.3) +
      ylim(0, 1) +
      labs(title = "Axonal Integrity", x = "Weeks",
           y = "Axonal integrity (1=healthy)") +
      theme_bw()
    ggplotly(p)
  })

  ## ── Clinical endpoint plots ───────────────────────────────────────────────────
  output$edss_plot <- renderPlotly({
    req(sim_clin())
    p <- ggplot(sim_clin(), aes(x = time_wk, y = EDSS)) +
      geom_line(color = drug_colors[input$drug_clin], linewidth = 1.4) +
      geom_ribbon(aes(ymin = EDSS * 0.9, ymax = EDSS * 1.1),
                  fill = drug_colors[input$drug_clin], alpha = 0.15) +
      ylim(0, 6) +
      labs(title = paste("EDSS Progression —", drug_labels[input$drug_clin]),
           x = "Weeks", y = "EDSS Score") +
      theme_bw()
    ggplotly(p)
  })

  output$arr_plot <- renderPlotly({
    req(sim_clin())
    p <- ggplot(sim_clin(), aes(x = time_wk, y = ARR)) +
      geom_line(color = drug_colors[input$drug_clin], linewidth = 1.2) +
      geom_hline(yintercept = 0.5, linetype = 2, color = "grey40") +
      labs(title = "Annualized Relapse Rate", x = "Weeks",
           y = "ARR (events/year)") +
      theme_bw()
    ggplotly(p)
  })

  output$neda_plot <- renderPlotly({
    req(sim_clin())
    d <- sim_clin() %>%
      mutate(
        no_relapse = as.integer(ARR < 0.3),
        no_mri     = as.integer(T2lesion < quantile(T2lesion, 0.25, na.rm=TRUE) * 1.2),
        no_progress = as.integer(EDSS <= first(EDSS) + 0.5),
        NEDA3 = as.integer(no_relapse & no_mri & no_progress)
      )
    p <- ggplot(d, aes(x = time_wk, y = cumsum(NEDA3)/pmax(row_number(), 1) * 100)) +
      geom_line(color = "#27AE60", linewidth = 1.3) +
      labs(title = "NEDA-3 Achievement Rate (Cumulative %)",
           x = "Weeks", y = "% Visits with NEDA-3") +
      ylim(0, 100) + theme_bw()
    ggplotly(p)
  })

  ## ── Scenario comparison ───────────────────────────────────────────────────────
  output$comp_plot <- renderPlotly({
    req(sim_comp())
    d <- sim_comp()
    yvar <- input$comp_endpoint
    if (!(yvar %in% names(d))) yvar <- "EDSS"
    p <- ggplot(d, aes(x = time_wk, y = .data[[yvar]],
                        color = drug, linetype = drug)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = drug_colors, labels = drug_labels,
                         name = "Treatment") +
      scale_linetype_manual(values = c(none="solid",ifnb="dashed",
                                       nat="solid",ocre="dotdash",
                                       sip="dashed",dmf="dotted",
                                       clad="longdash"),
                            labels = drug_labels, name = "Treatment") +
      labs(title = paste("Scenario Comparison:", yvar),
           x = "Weeks", y = yvar) +
      theme_bw()
    ggplotly(p)
  })

  output$summary_table <- renderDT({
    req(sim_comp())
    d <- sim_comp()
    baseline_arr <- d %>%
      filter(drug == "none") %>%
      summarise(base_arr = mean(ARR, na.rm = TRUE)) %>%
      pull(base_arr)

    tbl <- d %>%
      group_by(drug) %>%
      summarise(
        Treatment    = drug_labels[first(drug)],
        `Mean EDSS`  = round(mean(EDSS, na.rm = TRUE), 2),
        `Final EDSS` = round(last(EDSS), 2),
        `Mean ARR`   = round(mean(ARR, na.rm = TRUE), 3),
        `ARR Reduction (%)` = round((1 - mean(ARR, na.rm=TRUE)/baseline_arr)*100, 1),
        `Final Myelin`= round(last(Myelin), 3),
        `Final NfL`  = round(last(NfL), 1),
        `Final T2 (mm³)` = round(last(T2lesion), 0),
        .groups = "drop"
      ) %>%
      select(-drug) %>%
      arrange(`Mean ARR`)

    datatable(tbl, options = list(dom = "t", pageLength = 10),
              class = "table-striped table-hover") %>%
      formatStyle(
        "ARR Reduction (%)",
        background = styleColorBar(c(0, 80), "#27AE60"),
        backgroundSize = "100% 88%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center"
      )
  })

  ## ── Biomarkers ────────────────────────────────────────────────────────────────
  output$bm_plot <- renderPlotly({
    req(sim_bm())
    d <- sim_bm() %>%
      select(time_wk, NfL, GFAP) %>%
      pivot_longer(c(NfL, GFAP), names_to = "biomarker", values_to = "value")
    p <- ggplot(d, aes(x = time_wk, y = value, color = biomarker)) +
      geom_line(linewidth = 1.2) +
      facet_wrap(~biomarker, scales = "free_y") +
      scale_color_manual(values = c(NfL = "#E74C3C", GFAP = "#3498DB")) +
      geom_hline(data = data.frame(biomarker = c("NfL","GFAP"),
                                   yint = c(10, 150)),
                 aes(yintercept = yint), linetype = 2, color = "grey40") +
      labs(title = "Neurodegeneration Biomarkers", x = "Weeks",
           y = "Concentration (pg/mL)") +
      theme_bw()
    ggplotly(p)
  })

  output$lymph_plot <- renderPlotly({
    req(sim_bm())
    # Use siponimod for lymphocyte monitoring
    d_sip <- tryCatch(
      run_sim("sip", input$duration_wk, input$weight_kg,
              input$th1_init, input$th17_init, input$myelin_init),
      error = function(e) sim_bm()
    )
    p <- ggplot(d_sip, aes(x = time_wk, y = (Th1 + Th17 + Bcell + Treg) * 0.7)) +
      geom_line(color = "#F39C12", linewidth = 1.3) +
      geom_hline(yintercept = 0.2, linetype = 2, color = "red") +
      annotate("text", x = 2, y = 0.25, label = "Safety threshold (0.2×10⁹/L)",
               color = "red", size = 3) +
      labs(title = "Siponimod: Lymphocyte Safety Monitoring",
           x = "Weeks", y = "Lymphocyte count index") +
      theme_bw()
    ggplotly(p)
  })

  output$occ_plot <- renderPlotly({
    req(sim_bm())
    d <- sim_bm()
    if ("occ_nat" %in% names(d) && any(d$occ_nat > 0, na.rm = TRUE)) {
      p <- ggplot(d, aes(x = time_wk, y = occ_nat)) +
        geom_line(color = "#E74C3C", linewidth = 1.2) +
        geom_hline(yintercept = 80, linetype = 2, color = "grey40") +
        ylim(0, 100) +
        labs(title = "Receptor/Target Occupancy (%)",
             x = "Weeks", y = "% Occupancy") +
        theme_bw()
    } else {
      p <- ggplot(d, aes(x = time_wk, y = E_ocre_pct)) +
        geom_line(color = "#9B59B6", linewidth = 1.2) +
        labs(title = "Ocrelizumab B cell depletion (%)",
             x = "Weeks", y = "% Depletion") +
        theme_bw()
    }
    ggplotly(p)
  })
}

## ── Launch ────────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
