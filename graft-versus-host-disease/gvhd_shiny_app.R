## ============================================================
## GvHD QSP Shiny Dashboard
## Graft-versus-Host Disease Interactive Simulator
## ============================================================
## Tabs:
##  1. Patient & HSCT Profile
##  2. Drug PK Dashboard
##  3. Immune Cell Dynamics
##  4. Cytokine Network
##  5. Organ Damage & Clinical Endpoints
##  6. Scenario Comparison
##  7. Biomarkers
##  8. Mechanistic Map
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
## mrgsolve model (embedded)
## ============================================================

gvhd_code <- '
$PARAM @annotated
CsA_F:0.30:CsA bioavailability
CsA_ka:0.80:CsA absorption rate (1/h)
CsA_CL:25.0:CsA clearance (L/h)
CsA_V1:300:CsA central volume (L)
CsA_Q:50.0:CsA inter-cmpt CL (L/h)
CsA_V2:2000:CsA peripheral volume (L)
TAC_F:0.25:TAC bioavailability
TAC_ka:0.60:TAC absorption rate (1/h)
TAC_CL:2.50:TAC clearance (L/h)
TAC_V1:30.0:TAC central volume (L)
TAC_Q:3.00:TAC inter-cmpt CL (L/h)
TAC_V2:200:TAC peripheral volume (L)
PRED_F:0.99:Prednisone bioavailability
PRED_ka:1.20:PRED absorption rate (1/h)
PRED_CL:15.0:PRED clearance (L/h)
PRED_V1:50.0:PRED central volume (L)
RUX_F:0.95:Ruxolitinib bioavailability
RUX_ka:2.00:RUX absorption rate (1/h)
RUX_CL:17.7:RUX clearance (L/h)
RUX_V1:72.0:RUX central volume (L)
RUX_Q:10.0:RUX inter-cmpt CL (L/h)
RUX_V2:100:RUX peripheral volume (L)
BELU_F:0.80:Belumosudil bioavailability
BELU_ka:0.50:BELU absorption rate (1/h)
BELU_CL:6.00:BELU clearance (L/h)
BELU_V1:175:BELU central volume (L)
MMF_F:0.94:MMF/MPA bioavailability
MMF_ka:1.50:MMF absorption rate (1/h)
MMF_CL:12.0:MPA clearance (L/h)
MMF_V1:3.50:MPA central volume (L)
IC50_CsA:150.0:CsA IC50 calcineurin (ng/mL)
IC50_TAC:10.0:TAC IC50 calcineurin (ng/mL)
Imax_CN:0.85:Max calcineurin inhibition
n_CN:1.50:CNI Hill coefficient
IC50_PRED:50.0:PRED IC50 cytokine supp (ng/mL)
Imax_PRED:0.80:Max PRED effect
IC50_RUX:280.0:RUX IC50 JAK (nM)
Imax_RUX:0.90:Max JAK inhibition
IC50_BELU:100.0:BELU IC50 ROCK2 (nM)
Imax_BELU:0.85:Max ROCK2 inhibition
IC50_MPA:1.50:MPA IC50 proliferation (ug/mL)
Imax_MPA:0.85:Max MMF/MPA effect
kact_T:0.03:T cell activation rate (1/h)
kprol_T:0.05:T cell proliferation rate (1/h)
kdiff_Th1:0.02:Th1 differentiation rate
kdiff_Th17:0.015:Th17 differentiation rate
kdiff_Treg:0.010:Treg differentiation rate
kdeath_T:0.008:T cell apoptosis rate (1/h)
kprod_TNFa:0.040:TNF-a production rate
kprod_IFNg:0.035:IFN-g production rate
kprod_IL17:0.020:IL-17A production rate
kprod_IL10:0.025:IL-10 production rate
kprod_TGFb:0.015:TGF-b production rate
kprod_IL6:0.030:IL-6 production rate
kdeg_Cyt:0.30:Cytokine degradation rate (1/h)
kdam_skin:0.005:Skin damage rate
krep_skin:0.002:Skin repair rate
kdam_gut:0.008:Gut damage rate
krep_gut:0.003:Gut repair rate
kdam_liver:0.004:Liver damage rate
krep_liver:0.002:Liver repair rate
kdam_lung:0.003:Lung damage rate
krep_lung:0.001:Lung repair rate
kfibr:0.002:Fibrosis rate
Allo_stim:1.00:Alloantigen stimulus
GvL_effect:0.30:GvL preservation

$CMT CsA_gut CsA_cent CsA_peri TAC_gut TAC_cent TAC_peri
PRED_gut PRED_cent RUX_gut RUX_cent RUX_peri BELU_gut BELU_cent
MPA_gut MPA_cent
Th1 Th17 Treg CD8_eff Bcell
TNFa IFNg IL17A IL10 TGFb IL6
Skin_dmg Gut_dmg Liver_dmg Lung_dmg Fibrosis

$INIT CsA_gut=0,CsA_cent=0,CsA_peri=0,TAC_gut=0,TAC_cent=0,TAC_peri=0
PRED_gut=0,PRED_cent=0,RUX_gut=0,RUX_cent=0,RUX_peri=0,BELU_gut=0,BELU_cent=0
MPA_gut=0,MPA_cent=0
Th1=1,Th17=1,Treg=1,CD8_eff=1,Bcell=1
TNFa=0.1,IFNg=0.1,IL17A=0.1,IL10=0.5,TGFb=0.5,IL6=0.1
Skin_dmg=0,Gut_dmg=0,Liver_dmg=0,Lung_dmg=0,Fibrosis=0

$ODE
double CsA_conc = CsA_cent / CsA_V1 * 1000;
double TAC_conc = TAC_cent / TAC_V1 * 1000;
double PRED_conc= PRED_cent/ PRED_V1* 1000;
double RUX_conc_nM = (RUX_cent / RUX_V1 * 1000) / 0.306;
double BELU_conc_nM= (BELU_cent/ BELU_V1* 1000) / 0.532;
double MPA_conc = MPA_cent / MMF_V1;
double E_CsA = Imax_CN * pow(CsA_conc,n_CN)/(pow(IC50_CsA,n_CN)+pow(CsA_conc,n_CN));
double E_TAC = Imax_CN * pow(TAC_conc,n_CN)/(pow(IC50_TAC,n_CN)+pow(TAC_conc,n_CN));
double E_CNI = 1-(1-E_CsA)*(1-E_TAC);
double E_PRED= Imax_PRED*PRED_conc/(IC50_PRED+PRED_conc);
double E_RUX = Imax_RUX *RUX_conc_nM/(IC50_RUX +RUX_conc_nM);
double E_BELU= Imax_BELU*BELU_conc_nM/(IC50_BELU+BELU_conc_nM);
double E_MPA = Imax_MPA *MPA_conc/(IC50_MPA+MPA_conc);
double E_drug_inflam = 1-(1-E_CNI)*(1-E_PRED)*(1-E_RUX);
double E_drug_prolif = 1-(1-E_CNI)*(1-E_MPA);
double E_Treg_drug = 0.3*E_RUX + 0.4*E_BELU;
dxdt_CsA_gut  = -CsA_ka*CsA_gut;
dxdt_CsA_cent =  CsA_F*CsA_ka*CsA_gut-(CsA_CL+CsA_Q)/CsA_V1*CsA_cent+CsA_Q/CsA_V2*CsA_peri;
dxdt_CsA_peri =  CsA_Q/CsA_V1*CsA_cent-CsA_Q/CsA_V2*CsA_peri;
dxdt_TAC_gut  = -TAC_ka*TAC_gut;
dxdt_TAC_cent =  TAC_F*TAC_ka*TAC_gut-(TAC_CL+TAC_Q)/TAC_V1*TAC_cent+TAC_Q/TAC_V2*TAC_peri;
dxdt_TAC_peri =  TAC_Q/TAC_V1*TAC_cent-TAC_Q/TAC_V2*TAC_peri;
dxdt_PRED_gut  = -PRED_ka*PRED_gut;
dxdt_PRED_cent =  PRED_F*PRED_ka*PRED_gut-PRED_CL/PRED_V1*PRED_cent;
dxdt_RUX_gut  = -RUX_ka*RUX_gut;
dxdt_RUX_cent =  RUX_F*RUX_ka*RUX_gut-(RUX_CL+RUX_Q)/RUX_V1*RUX_cent+RUX_Q/RUX_V2*RUX_peri;
dxdt_RUX_peri =  RUX_Q/RUX_V1*RUX_cent-RUX_Q/RUX_V2*RUX_peri;
dxdt_BELU_gut  = -BELU_ka*BELU_gut;
dxdt_BELU_cent =  BELU_F*BELU_ka*BELU_gut-BELU_CL/BELU_V1*BELU_cent;
dxdt_MPA_gut  = -MMF_ka*MPA_gut;
dxdt_MPA_cent =  MMF_F*MMF_ka*MPA_gut-MMF_CL/MMF_V1*MPA_cent;
double T_activ = kact_T*Allo_stim*(1-E_drug_prolif);
dxdt_Th1    = kprol_T*Allo_stim*(1+0.5*IFNg/(IFNg+1))*(1-E_drug_prolif)*Th1 + T_activ*kdiff_Th1 - (kdeath_T+0.2*IL10/(IL10+1))*Th1;
dxdt_Th17   = kprol_T*Allo_stim*(1+0.4*IL6/(IL6+2))*(1-E_RUX)*(1-E_BELU*0.7)*(1-E_MPA)*Th17 + T_activ*kdiff_Th17 - (kdeath_T+0.3*IL10/(IL10+1))*Th17;
dxdt_Treg   = kprol_T*TGFb/(TGFb+1)*(1+E_Treg_drug)*(1-E_MPA*0.3)*Treg + T_activ*kdiff_Treg - (kdeath_T+0.15*TNFa/(TNFa+2)+0.15*IL6/(IL6+2))*Treg;
dxdt_CD8_eff= kprol_T*Allo_stim*(1+0.3*IFNg/(IFNg+1))*(1-E_drug_prolif)*CD8_eff + T_activ*0.5 - (kdeath_T+0.2*IL10/(IL10+1))*CD8_eff;
dxdt_Bcell  = kprol_T*0.4*(1-E_MPA)*Bcell - kdeath_T*0.5*Bcell;
dxdt_TNFa  = kprod_TNFa*(Th1+0.5*CD8_eff)*(1-E_drug_inflam) - kdeg_Cyt*TNFa;
dxdt_IFNg  = kprod_IFNg*(Th1+0.7*CD8_eff)*(1-E_PRED*0.6-E_RUX*0.4) - kdeg_Cyt*IFNg;
dxdt_IL17A = kprod_IL17*Th17*(1-E_RUX*0.7-E_BELU*0.5)*(1-E_PRED*0.4) - kdeg_Cyt*IL17A;
dxdt_IL10  = kprod_IL10*Treg*(1+E_Treg_drug) - kdeg_Cyt*IL10;
dxdt_TGFb  = kprod_TGFb*(Treg+0.3*Bcell) - kdeg_Cyt*TGFb;
dxdt_IL6   = kprod_IL6*(1+0.5*TNFa/(TNFa+1))*(1-E_RUX*0.5-E_PRED*0.5) - kdeg_Cyt*IL6;
double inflam_driver = (TNFa+IFNg+0.7*IL17A)/3;
double anti_inflam = IL10;
dxdt_Skin_dmg  = kdam_skin*(Th1+CD8_eff)*inflam_driver*(1-E_drug_inflam)*(1-Skin_dmg) - krep_skin*Treg*(anti_inflam+1)*Skin_dmg;
dxdt_Gut_dmg   = kdam_gut*inflam_driver*CD8_eff*(1-E_drug_inflam)*(1-Gut_dmg) - krep_gut*(1+anti_inflam)*(1-Fibrosis*0.3)*Gut_dmg;
dxdt_Liver_dmg = kdam_liver*(TNFa+IFNg)*0.5*(1-E_drug_inflam)*(1-Liver_dmg) - krep_liver*(1+anti_inflam*0.5)*Liver_dmg;
dxdt_Lung_dmg  = kdam_lung*(IL17A+0.3*TNFa)*(1-E_drug_inflam*0.3)*(1-Lung_dmg) - krep_lung*Lung_dmg;
dxdt_Fibrosis  = kfibr*TGFb*(1-E_BELU)*(1-Fibrosis) - 0.0005*Fibrosis;

$TABLE
double CsA_C0 = CsA_cent/CsA_V1*1000;
double TAC_C0 = TAC_cent/TAC_V1*1000;
double PRED_C = PRED_cent/PRED_V1*1000;
double RUX_C  = RUX_cent/RUX_V1*1000;
double BELU_C = BELU_cent/BELU_V1*1000;
double MPA_C  = MPA_cent/MMF_V1;
double E_CNI_out  = 1-(1-(Imax_CN*pow(CsA_C0,n_CN)/(pow(IC50_CsA,n_CN)+pow(CsA_C0,n_CN))))*(1-(Imax_CN*pow(TAC_C0,n_CN)/(pow(IC50_TAC,n_CN)+pow(TAC_C0,n_CN))));
double E_RUX_out  = Imax_RUX*(RUX_C/0.306)/(IC50_RUX+(RUX_C/0.306));
double E_BELU_out = Imax_BELU*(BELU_C/0.532)/(IC50_BELU+(BELU_C/0.532));
double aGvHD_Grade = (Skin_dmg+Gut_dmg+Liver_dmg)/3.0*3.0;
double cGvHD_Score = (Skin_dmg+Gut_dmg+Liver_dmg+Lung_dmg+Fibrosis)/5.0*3.0;
double ST2_bio = 10+500*Gut_dmg;
double REG3a_bio = 10+200*Gut_dmg;
double TNFR1_bio = 1+5*TNFa;

$CAPTURE CsA_C0 TAC_C0 RUX_C BELU_C PRED_C MPA_C
E_CNI_out E_RUX_out E_BELU_out
aGvHD_Grade cGvHD_Score
ST2_bio REG3a_bio TNFR1_bio
Th1 Th17 Treg CD8_eff Bcell
TNFa IFNg IL17A IL10 TGFb IL6
Skin_dmg Gut_dmg Liver_dmg Lung_dmg Fibrosis
'

mod <- mcode("gvhd_shiny", gvhd_code, quiet = TRUE)

## ============================================================
## Helper: run simulation
## ============================================================

run_sim <- function(
    weight = 70, allo_stim = 1.0,
    csa_dose = 0, csa_dur = 100,
    tac_dose = 0, tac_dur = 100,
    pred_dose = 0, pred_dur = 14, pred_start = 20,
    rux_dose = 0, rux_dur = 180, rux_start = 30,
    belu_dose = 0, belu_dur = 180, belu_start = 30,
    mmf_dose = 0, mmf_dur = 100,
    sim_days = 365) {

  evs <- ev(amt = 0, cmt = 1, time = 0)  # placeholder

  if (csa_dose > 0) {
    csa_total <- csa_dose * weight  # mg/day
    evs <- ev_seq(evs, ev(amt = csa_total/2, cmt = 1, ii = 12,
                          addl = ceiling(csa_dur * 24/12) - 1, time = 0))
  }
  if (tac_dose > 0) {
    tac_total <- tac_dose * weight  # mg/day
    evs <- ev_seq(evs, ev(amt = tac_total/2, cmt = 4, ii = 12,
                          addl = ceiling(tac_dur * 24/12) - 1, time = 0))
  }
  if (pred_dose > 0) {
    pred_total <- pred_dose * weight
    evs <- ev_seq(evs, ev(amt = pred_total, cmt = 7, ii = 24,
                          addl = pred_dur - 1, time = pred_start * 24))
  }
  if (rux_dose > 0) {
    evs <- ev_seq(evs, ev(amt = rux_dose, cmt = 9, ii = 12,
                          addl = ceiling(rux_dur * 24/12) - 1, time = rux_start * 24))
  }
  if (belu_dose > 0) {
    evs <- ev_seq(evs, ev(amt = belu_dose, cmt = 12, ii = 24,
                          addl = belu_dur - 1, time = belu_start * 24))
  }
  if (mmf_dose > 0) {
    evs <- ev_seq(evs, ev(amt = mmf_dose, cmt = 14, ii = 12,
                          addl = ceiling(mmf_dur * 24/12) - 1, time = 0))
  }

  mrgsim(
    mod %>% param(Allo_stim = allo_stim),
    ev = evs,
    end = sim_days * 24, delta = 12,
    obsaug = TRUE
  ) %>% as_tibble() %>% mutate(Day = time / 24)
}

## ============================================================
## UI
## ============================================================

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "GvHD QSP Dashboard", titleWidth = 280),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient & HSCT Profile",  tabName = "patient",   icon = icon("user")),
      menuItem("Drug PK Dashboard",        tabName = "pk",        icon = icon("pills")),
      menuItem("Immune Cell Dynamics",     tabName = "immune",    icon = icon("microscope")),
      menuItem("Cytokine Network",         tabName = "cytokine",  icon = icon("project-diagram")),
      menuItem("Organ Damage & Endpoints", tabName = "organ",     icon = icon("heartbeat")),
      menuItem("Scenario Comparison",      tabName = "scenario",  icon = icon("chart-line")),
      menuItem("Biomarkers",               tabName = "biomarker", icon = icon("vial")),
      menuItem("Mechanistic Map",          tabName = "map",       icon = icon("sitemap"))
    ),
    hr(),
    h5("Patient Parameters", style = "margin-left:15px; color:#ddd;"),
    sliderInput("weight", "Body Weight (kg)", 40, 120, 70, step = 5),
    sliderInput("allo_stim", "Alloantigen Stimulus (0-1)", 0, 1, 1, step = 0.1),
    sliderInput("sim_days", "Simulation Duration (days)", 90, 730, 365, step = 30),
    hr(),
    h5("Prophylaxis Drugs", style = "margin-left:15px; color:#ddd;"),
    sliderInput("csa_dose", "CsA (mg/kg/day)", 0, 8, 3, step = 0.5),
    sliderInput("csa_dur", "CsA Duration (days)", 0, 365, 100, step = 10),
    sliderInput("tac_dose", "TAC (mg/kg/day x100)", 0, 8, 0, step = 0.5),
    sliderInput("mmf_dose", "MMF/dose (mg BID)", 0, 3000, 1500, step = 250),
    hr(),
    h5("Treatment Drugs", style = "margin-left:15px; color:#ddd;"),
    sliderInput("pred_dose", "Prednisone (mg/kg/day)", 0, 2, 0, step = 0.1),
    sliderInput("pred_dur", "Pred Duration (days)", 7, 90, 14, step = 7),
    sliderInput("pred_start", "Pred Start (day post-HSCT)", 1, 60, 20, step = 1),
    sliderInput("rux_dose", "Ruxolitinib (mg BID)", 0, 25, 0, step = 5),
    sliderInput("rux_start", "Rux Start Day", 1, 180, 30, step = 5),
    sliderInput("belu_dose", "Belumosudil (mg QD)", 0, 400, 0, step = 100),
    sliderInput("belu_start", "Belu Start Day", 1, 360, 30, step = 5),
    actionButton("run_sim", "Run Simulation", class = "btn-primary btn-block",
                 icon = icon("play"))
  ),

  dashboardBody(
    tabItems(
      ## TAB 1: Patient & HSCT Profile
      tabItem(tabName = "patient",
        fluidRow(
          box(title = "Patient & HSCT Profile", status = "primary", solidHeader = TRUE, width = 12,
              fluidRow(
                valueBoxOutput("vb_weight", width = 3),
                valueBoxOutput("vb_allo", width = 3),
                valueBoxOutput("vb_sim_days", width = 3),
                valueBoxOutput("vb_hsct_type", width = 3)
              )
          )
        ),
        fluidRow(
          box(title = "GvHD Risk Overview", status = "info", solidHeader = TRUE, width = 6,
              plotlyOutput("risk_radar", height = "350px")),
          box(title = "Key HSCT / GvHD Facts", status = "info", solidHeader = TRUE, width = 6,
              h4("GvHD Epidemiology"),
              p("• Acute GvHD occurs in 30-50% of matched sibling, 50-70% of MUD transplants"),
              p("• Chronic GvHD develops in 40-70% of aGvHD survivors"),
              p("• GvHD is the leading cause of non-relapse mortality (NRM) post-HSCT"),
              h4("Key Drug Targets"),
              tags$ul(
                tags$li(strong("Calcineurin (CsA/TAC):"), " blocks NFAT → IL-2 transcription"),
                tags$li(strong("JAK1/2 (Ruxolitinib):"), " blocks IFN-γ, IL-6, IL-12 signaling; expands Treg"),
                tags$li(strong("ROCK2 (Belumosudil):"), " IRF4/STAT3 → Th17↓, Treg↑, fibrosis↓"),
                tags$li(strong("BTK (Ibrutinib):"), " B cell GC reaction, Th17 pathway"),
                tags$li(strong("IMPDH (MMF):"), " lymphocyte de novo purine synthesis blocked")
              ),
              h4("FDA-approved for GvHD:"),
              tags$ul(
                tags$li("Ruxolitinib (REACH1, 2, 3 trials)"),
                tags$li("Belumosudil (ROCKstar trial)"),
                tags$li("Ibrutinib (chronic GvHD 2L+)")
              )
          )
        )
      ),

      ## TAB 2: Drug PK Dashboard
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Drug Concentration-Time Profiles", status = "primary", solidHeader = TRUE, width = 12,
              fluidRow(
                column(6, plotlyOutput("pk_csa_plot", height = "280px")),
                column(6, plotlyOutput("pk_rux_plot", height = "280px"))
              ),
              fluidRow(
                column(6, plotlyOutput("pk_tac_plot", height = "280px")),
                column(6, plotlyOutput("pk_belu_plot", height = "280px"))
              )
          )
        ),
        fluidRow(
          box(title = "PK Summary Statistics", status = "info", solidHeader = TRUE, width = 6,
              DT::dataTableOutput("pk_summary_table")),
          box(title = "Drug Effect (PD) Over Time", status = "success", solidHeader = TRUE, width = 6,
              plotlyOutput("pd_effects_plot", height = "300px"))
        )
      ),

      ## TAB 3: Immune Cell Dynamics
      tabItem(tabName = "immune",
        fluidRow(
          box(title = "T Cell Subset Dynamics", status = "warning", solidHeader = TRUE, width = 8,
              plotlyOutput("tcell_plot", height = "400px")),
          box(title = "Th17/Treg Balance (cGvHD key)", status = "danger", solidHeader = TRUE, width = 4,
              plotlyOutput("th17_treg_ratio", height = "400px"))
        ),
        fluidRow(
          box(title = "B Cell Dynamics", status = "info", solidHeader = TRUE, width = 6,
              plotlyOutput("bcell_plot", height = "300px")),
          box(title = "Immunological Milieu Summary", status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("immune_heatmap", height = "300px"))
        )
      ),

      ## TAB 4: Cytokine Network
      tabItem(tabName = "cytokine",
        fluidRow(
          box(title = "Pro-inflammatory Cytokines", status = "danger", solidHeader = TRUE, width = 6,
              plotlyOutput("pro_cytokine_plot", height = "350px")),
          box(title = "Anti-inflammatory Cytokines", status = "success", solidHeader = TRUE, width = 6,
              plotlyOutput("anti_cytokine_plot", height = "350px"))
        ),
        fluidRow(
          box(title = "Cytokine Summary Table (Day 30 & 90 & 180)", status = "info", solidHeader = TRUE, width = 12,
              DT::dataTableOutput("cytokine_table"))
        )
      ),

      ## TAB 5: Organ Damage & Endpoints
      tabItem(tabName = "organ",
        fluidRow(
          valueBoxOutput("vb_aGvHD", width = 3),
          valueBoxOutput("vb_cGvHD", width = 3),
          valueBoxOutput("vb_fibrosis", width = 3),
          valueBoxOutput("vb_ffs", width = 3)
        ),
        fluidRow(
          box(title = "Organ Damage Scores Over Time", status = "danger", solidHeader = TRUE, width = 8,
              plotlyOutput("organ_damage_plot", height = "380px")),
          box(title = "Organ Damage at Key Timepoints", status = "info", solidHeader = TRUE, width = 4,
              DT::dataTableOutput("organ_table"))
        ),
        fluidRow(
          box(title = "aGvHD Grade (Glucksberg Proxy)", status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("agvhd_grade_plot", height = "300px")),
          box(title = "cGvHD Score (NIH Global Proxy)", status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("cgvhd_score_plot", height = "300px"))
        )
      ),

      ## TAB 6: Scenario Comparison
      tabItem(tabName = "scenario",
        fluidRow(
          box(title = "Preset Scenarios", status = "primary", solidHeader = TRUE, width = 12,
              p("Click 'Run Comparison' to simulate all 6 preset treatment scenarios"),
              actionButton("run_compare", "Run All Scenarios", class = "btn-success btn-lg",
                           icon = icon("chart-bar"))
          )
        ),
        fluidRow(
          box(title = "aGvHD Grade Comparison", status = "danger", solidHeader = TRUE, width = 6,
              plotlyOutput("compare_agvhd", height = "350px")),
          box(title = "cGvHD Score Comparison", status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("compare_cgvhd", height = "350px"))
        ),
        fluidRow(
          box(title = "Fibrosis Development", status = "info", solidHeader = TRUE, width = 6,
              plotlyOutput("compare_fibr", height = "300px")),
          box(title = "Summary Outcomes Table", status = "success", solidHeader = TRUE, width = 6,
              DT::dataTableOutput("compare_table"))
        )
      ),

      ## TAB 7: Biomarkers
      tabItem(tabName = "biomarker",
        fluidRow(
          box(title = "GvHD Biomarker Panel Over Time", status = "primary", solidHeader = TRUE, width = 12,
              plotlyOutput("biomarker_plot", height = "400px"))
        ),
        fluidRow(
          box(title = "Biomarker Reference Thresholds", status = "info", solidHeader = TRUE, width = 6,
              DT::dataTableOutput("biomarker_table")),
          box(title = "Biomarker Interpretation", status = "warning", solidHeader = TRUE, width = 6,
              h4("ST2 (Soluble IL-33 Receptor)"),
              p("• Cutoff: >33 ng/mL → predicts grade 3-4 GI GvHD"),
              p("• Ann Arbor algorithm day+7/14 post-HSCT predicts GvHD grade"),
              h4("REG3α (Regenerating Islet-Derived Protein 3α)"),
              p("• Paneth cell marker of intestinal damage"),
              p("• >>23 ng/mL associated with poor GI GvHD outcome"),
              h4("sTNFR1"),
              p("• Shed from TNF-α stimulated cells"),
              p("• Predicts skin and systemic GvHD severity"),
              h4("CXCL9"),
              p("• IFN-γ-inducible chemokine for T cell homing"),
              p("• Predicts aGvHD in peripheral blood")
          )
        )
      ),

      ## TAB 8: Mechanistic Map
      tabItem(tabName = "map",
        fluidRow(
          box(title = "GvHD Mechanistic Network Overview", status = "primary", solidHeader = TRUE, width = 12,
              h4("Full Mechanistic Map"),
              p("The Graphviz DOT file (gvhd_qsp_model.dot) encodes 100+ nodes across 14 pathophysiological clusters:"),
              tags$ol(
                tags$li(strong("HSCT Context & Conditioning"), " - conditioning regimen, tissue damage, DAMP/PAMP release"),
                tags$li(strong("Antigen Presentation & T Cell Priming"), " - host DCs, direct/indirect alloreactivity, co-stimulation"),
                tags$li(strong("Donor T Cell Activation & Differentiation"), " - Th1/Th17/Treg/CD8 polarization, intracellular signaling (NFAT, NF-κB, JAK-STAT)"),
                tags$li(strong("Cytokine Network"), " - TNF-α, IFN-γ, IL-6, IL-17A, IL-10, TGF-β, BAFF"),
                tags$li(strong("Target Organ - Skin"), " - lichenoid/sclerotic lesions, mLSS score, pruritus"),
                tags$li(strong("Target Organ - Gut"), " - enterocyte apoptosis, crypt damage, ST2/REG3α biomarkers"),
                tags$li(strong("Target Organ - Liver"), " - bile duct damage, cholestasis, Glucksberg grading"),
                tags$li(strong("Target Organ - Lung"), " - bronchiolitis obliterans (BOS), FEV1, CLAD score"),
                tags$li(strong("B Cell Pathology"), " - Tfh-B GC reaction, autoantibodies, Breg, BTK pathway"),
                tags$li(strong("Fibrosis Pathway"), " - TGF-β/SMAD, ROCK2, EMT, myofibroblast activation"),
                tags$li(strong("Drug PK: CNI"), " - CsA 2-compartment, TAC 2-compartment, CYP3A4/5 metabolism"),
                tags$li(strong("Drug PK/PD: Ruxolitinib"), " - JAK1/2 inhibition, STAT3/5 blockade, Treg expansion"),
                tags$li(strong("Other Drugs"), " - Corticosteroids (NF-κB), Belumosudil (ROCK2), MMF (IMPDH), Ibrutinib (BTK)"),
                tags$li(strong("Clinical Endpoints & Biomarkers"), " - Glucksberg/NIH scores, ORR, FFS, OS, NRM")
              ),
              hr(),
              tags$img(src = "gvhd_qsp_model.png", width = "100%",
                       alt = "GvHD Mechanistic Map (PNG) - Open SVG for interactive version"),
              p("Open gvhd_qsp_model.svg in a browser for the full interactive mechanistic map.")
          )
        )
      )
    )
  )
)

## ============================================================
## SERVER
## ============================================================

server <- function(input, output, session) {

  # Reactive: run simulation
  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Running GvHD simulation...", {
      run_sim(
        weight     = input$weight,
        allo_stim  = input$allo_stim,
        csa_dose   = input$csa_dose,
        csa_dur    = input$csa_dur,
        tac_dose   = input$tac_dose / 100,
        mmf_dose   = input$mmf_dose,
        pred_dose  = input$pred_dose,
        pred_dur   = input$pred_dur,
        pred_start = input$pred_start,
        rux_dose   = input$rux_dose,
        rux_start  = input$rux_start,
        belu_dose  = input$belu_dose,
        belu_start = input$belu_start,
        sim_days   = input$sim_days
      )
    })
  }, ignoreNULL = FALSE)

  # Tab 1 Value Boxes
  output$vb_weight    <- renderValueBox(valueBox(paste(input$weight, "kg"), "Body Weight", icon = icon("weight"), color = "blue"))
  output$vb_allo      <- renderValueBox(valueBox(input$allo_stim, "Alloantigen Stimulus", icon = icon("dna"), color = "red"))
  output$vb_sim_days  <- renderValueBox(valueBox(input$sim_days, "Simulation Days", icon = icon("calendar"), color = "green"))
  output$vb_hsct_type <- renderValueBox(valueBox("Allo-HSCT", "Transplant Type", icon = icon("syringe"), color = "purple"))

  # Risk radar chart
  output$risk_radar <- renderPlotly({
    dat <- sim_data()
    if (is.null(dat)) return(NULL)
    last_row <- tail(dat, 1)
    plot_ly(
      type = "scatterpolar", mode = "lines+markers",
      r = c(last_row$Skin_dmg, last_row$Gut_dmg, last_row$Liver_dmg,
            last_row$Lung_dmg, last_row$Fibrosis, last_row$Skin_dmg),
      theta = c("Skin", "Gut", "Liver", "Lung", "Fibrosis", "Skin"),
      fill = "toself", fillcolor = "rgba(231,76,60,0.3)",
      line = list(color = "#E74C3C")
    ) %>% layout(polar = list(radialaxis = list(range = c(0,1))),
                 title = paste0("End-of-Simulation Organ Damage\n(Day ", input$sim_days, ")"))
  })

  # Tab 2: PK plots
  make_pk_plot <- function(dat, col, title, ylab, target_lo = NA, target_hi = NA) {
    p <- plot_ly(dat, x = ~Day, y = ~get(col), type = "scatter", mode = "lines",
                 line = list(color = "#2980B9", width = 2)) %>%
      layout(title = title, xaxis = list(title = "Days Post-HSCT"),
             yaxis = list(title = ylab))
    if (!is.na(target_lo))
      p <- p %>% add_segments(x = 0, xend = input$sim_days, y = target_lo, yend = target_lo,
                               line = list(dash = "dash", color = "green"), name = "Target Low")
    if (!is.na(target_hi))
      p <- p %>% add_segments(x = 0, xend = input$sim_days, y = target_hi, yend = target_hi,
                               line = list(dash = "dash", color = "red"), name = "Target High")
    p
  }

  output$pk_csa_plot <- renderPlotly({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    make_pk_plot(dat, "CsA_C0", "CsA Trough Concentration", "CsA (ng/mL)", 100, 300)
  })
  output$pk_tac_plot <- renderPlotly({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    make_pk_plot(dat, "TAC_C0", "Tacrolimus Trough", "TAC (ng/mL)", 5, 15)
  })
  output$pk_rux_plot <- renderPlotly({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    make_pk_plot(dat, "RUX_C", "Ruxolitinib Concentration", "RUX (ng/mL)")
  })
  output$pk_belu_plot <- renderPlotly({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    make_pk_plot(dat, "BELU_C", "Belumosudil Concentration", "BELU (ng/mL)")
  })

  output$pk_summary_table <- DT::renderDataTable({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    tibble(
      Drug = c("CsA (ng/mL)", "TAC (ng/mL)", "Ruxolitinib (ng/mL)", "Belumosudil (ng/mL)"),
      Max  = c(max(dat$CsA_C0), max(dat$TAC_C0), max(dat$RUX_C), max(dat$BELU_C)),
      Min  = c(min(dat$CsA_C0), min(dat$TAC_C0), min(dat$RUX_C), min(dat$BELU_C)),
      Mean = c(mean(dat$CsA_C0), mean(dat$TAC_C0), mean(dat$RUX_C), mean(dat$BELU_C))
    ) %>% mutate(across(where(is.numeric), ~round(., 2)))
  }, options = list(pageLength = 5))

  output$pd_effects_plot <- renderPlotly({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    plot_ly(dat, x = ~Day) %>%
      add_lines(y = ~E_CNI_out * 100, name = "CNI Inhibition (%)", line = list(color = "#E74C3C")) %>%
      add_lines(y = ~E_RUX_out * 100, name = "JAK Inhibition (%)", line = list(color = "#3498DB")) %>%
      add_lines(y = ~E_BELU_out * 100, name = "ROCK2 Inhibition (%)", line = list(color = "#9B59B6")) %>%
      layout(title = "Drug PD Effects Over Time",
             xaxis = list(title = "Days"),
             yaxis = list(title = "% Inhibition", range = c(0,100)))
  })

  # Tab 3: Immune
  output$tcell_plot <- renderPlotly({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    plot_ly(dat, x = ~Day) %>%
      add_lines(y = ~Th1,     name = "Th1 (IFN-γ, TNF-α)", line = list(color = "#E74C3C")) %>%
      add_lines(y = ~Th17,    name = "Th17 (IL-17A)",       line = list(color = "#8E44AD")) %>%
      add_lines(y = ~Treg,    name = "Treg (FoxP3+)",       line = list(color = "#27AE60")) %>%
      add_lines(y = ~CD8_eff, name = "CD8+ CTL",            line = list(color = "#E67E22")) %>%
      layout(title = "T Cell Subset Dynamics",
             xaxis = list(title = "Days Post-HSCT"),
             yaxis = list(title = "Relative Pool Size"),
             legend = list(x = 0.7, y = 1))
  })

  output$th17_treg_ratio <- renderPlotly({
    dat <- sim_data() %>% filter(time %% 24 == 0) %>%
      mutate(Ratio = Th17 / pmax(Treg, 0.01))
    plot_ly(dat, x = ~Day, y = ~Ratio, type = "scatter", mode = "lines",
            line = list(color = "#C0392B", width = 2)) %>%
      add_segments(x = 0, xend = max(dat$Day), y = 1, yend = 1,
                   line = list(dash = "dash", color = "gray"), name = "Balanced (=1)") %>%
      layout(title = "Th17/Treg Ratio\n(>1 = cGvHD risk)",
             xaxis = list(title = "Days"), yaxis = list(title = "Th17/Treg Ratio"))
  })

  output$bcell_plot <- renderPlotly({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    plot_ly(dat, x = ~Day, y = ~Bcell, type = "scatter", mode = "lines",
            line = list(color = "#48C9B0", width = 2)) %>%
      layout(title = "B Cell Pool Dynamics", xaxis = list(title = "Days"),
             yaxis = list(title = "Relative B Cell Pool"))
  })

  output$immune_heatmap <- renderPlotly({
    dat <- sim_data() %>%
      filter(time %% (24*7) == 0) %>%
      select(Day, Th1, Th17, Treg, CD8_eff, Bcell) %>%
      pivot_longer(-Day, names_to = "Cell", values_to = "Value")
    plot_ly(dat, x = ~Day, y = ~Cell, z = ~Value, type = "heatmap",
            colorscale = "RdBu", reversescale = TRUE) %>%
      layout(title = "Immune Cell Heatmap (Weekly)", xaxis = list(title = "Days"))
  })

  # Tab 4: Cytokines
  output$pro_cytokine_plot <- renderPlotly({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    plot_ly(dat, x = ~Day) %>%
      add_lines(y = ~TNFa,  name = "TNF-α", line = list(color = "#E74C3C")) %>%
      add_lines(y = ~IFNg,  name = "IFN-γ", line = list(color = "#C0392B")) %>%
      add_lines(y = ~IL17A, name = "IL-17A",line = list(color = "#8E44AD")) %>%
      add_lines(y = ~IL6,   name = "IL-6",  line = list(color = "#E67E22")) %>%
      layout(title = "Pro-inflammatory Cytokines",
             xaxis = list(title = "Days"), yaxis = list(title = "Concentration (ng/mL)"))
  })

  output$anti_cytokine_plot <- renderPlotly({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    plot_ly(dat, x = ~Day) %>%
      add_lines(y = ~IL10, name = "IL-10 (anti-inflam)", line = list(color = "#27AE60")) %>%
      add_lines(y = ~TGFb, name = "TGF-β (Treg/fibrosis)", line = list(color = "#2980B9")) %>%
      layout(title = "Anti-inflammatory / Regulatory Cytokines",
             xaxis = list(title = "Days"), yaxis = list(title = "Concentration (ng/mL)"))
  })

  output$cytokine_table <- DT::renderDataTable({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    bind_rows(
      dat %>% filter(abs(Day - 30)  < 1) %>% head(1) %>% mutate(Timepoint = "Day 30"),
      dat %>% filter(abs(Day - 90)  < 1) %>% head(1) %>% mutate(Timepoint = "Day 90"),
      dat %>% filter(abs(Day - 180) < 1) %>% head(1) %>% mutate(Timepoint = "Day 180")
    ) %>% select(Timepoint, TNFa, IFNg, IL17A, IL10, TGFb, IL6) %>%
      mutate(across(where(is.numeric), ~round(., 3)))
  })

  # Tab 5: Organ damage
  output$vb_aGvHD <- renderValueBox({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    max_grade <- max(dat$aGvHD_Grade, na.rm = TRUE)
    valueBox(round(max_grade, 2), "Max aGvHD Grade (0-3)", icon = icon("exclamation-triangle"),
             color = ifelse(max_grade < 1, "green", ifelse(max_grade < 2, "yellow", "red")))
  })

  output$vb_cGvHD <- renderValueBox({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    last_val <- tail(dat$cGvHD_Score, 1)
    valueBox(round(last_val, 2), "Final cGvHD Score (0-3)", icon = icon("lungs"),
             color = ifelse(last_val < 1, "green", ifelse(last_val < 2, "yellow", "red")))
  })

  output$vb_fibrosis <- renderValueBox({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    last_val <- tail(dat$Fibrosis, 1)
    valueBox(paste0(round(last_val * 100), "%"), "Final Fibrosis Index", icon = icon("layer-group"),
             color = ifelse(last_val < 0.2, "green", ifelse(last_val < 0.5, "yellow", "red")))
  })

  output$vb_ffs <- renderValueBox({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    last_grade <- tail(dat$aGvHD_Grade, 1)
    last_cgvhd <- tail(dat$cGvHD_Score, 1)
    ffs <- ifelse(last_grade < 2 && last_cgvhd < 1.5, "Yes", "No")
    valueBox(ffs, "Failure-Free Survival", icon = icon("check"),
             color = ifelse(ffs == "Yes", "green", "red"))
  })

  output$organ_damage_plot <- renderPlotly({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    plot_ly(dat, x = ~Day) %>%
      add_lines(y = ~Skin_dmg  * 100, name = "Skin",   line = list(color = "#E67E22")) %>%
      add_lines(y = ~Gut_dmg   * 100, name = "Gut",    line = list(color = "#27AE60")) %>%
      add_lines(y = ~Liver_dmg * 100, name = "Liver",  line = list(color = "#F39C12")) %>%
      add_lines(y = ~Lung_dmg  * 100, name = "Lung",   line = list(color = "#8E44AD")) %>%
      add_lines(y = ~Fibrosis  * 100, name = "Fibrosis",line=list(color="#566573",dash="dash")) %>%
      layout(title = "Organ Damage Scores Over Time",
             xaxis = list(title = "Days Post-HSCT"),
             yaxis = list(title = "Damage Score (%)", range = c(0, 100)))
  })

  output$organ_table <- DT::renderDataTable({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    bind_rows(
      dat %>% filter(abs(Day - 30)  < 1) %>% head(1) %>% mutate(Day_label = "Day 30"),
      dat %>% filter(abs(Day - 90)  < 1) %>% head(1) %>% mutate(Day_label = "Day 90"),
      dat %>% filter(abs(Day - 180) < 1) %>% head(1) %>% mutate(Day_label = "Day 180"),
      tail(dat, 1) %>% mutate(Day_label = paste("Day", input$sim_days))
    ) %>% select(Day_label, Skin_dmg, Gut_dmg, Liver_dmg, Lung_dmg, Fibrosis) %>%
      mutate(across(where(is.numeric), ~round(. * 100, 1))) %>%
      rename(Day = Day_label, `Skin %` = Skin_dmg, `Gut %` = Gut_dmg,
             `Liver %` = Liver_dmg, `Lung %` = Lung_dmg, `Fibrosis %` = Fibrosis)
  }, options = list(pageLength = 4, dom = "t"))

  output$agvhd_grade_plot <- renderPlotly({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    plot_ly(dat, x = ~Day, y = ~aGvHD_Grade, type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(231,76,60,0.2)",
            line = list(color = "#E74C3C")) %>%
      add_segments(x = 0, xend = max(dat$Day), y = 2, yend = 2,
                   line = list(dash = "dash", color = "orange"), name = "Grade 2 threshold") %>%
      layout(title = "aGvHD Glucksberg Grade",
             yaxis = list(title = "Grade (0-3)", range = c(0, 3)))
  })

  output$cgvhd_score_plot <- renderPlotly({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    plot_ly(dat, x = ~Day, y = ~cGvHD_Score, type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(155,89,182,0.2)",
            line = list(color = "#8E44AD")) %>%
      layout(title = "cGvHD NIH Global Score",
             yaxis = list(title = "NIH Score (0-3)", range = c(0, 3)))
  })

  # Tab 6: Scenario comparison
  compare_data <- eventReactive(input$run_compare, {
    scenarios <- list(
      list(name = "No Prophylaxis",         csa = 0,   tac = 0,   mmf = 0,    rux = 0,  belu = 0, pred = 0),
      list(name = "CsA Monoproph",           csa = 3,   tac = 0,   mmf = 0,    rux = 0,  belu = 0, pred = 0),
      list(name = "CsA + MMF",              csa = 3,   tac = 0,   mmf = 1500, rux = 0,  belu = 0, pred = 0),
      list(name = "TAC + MMF (standard)",   csa = 0,   tac = 0.03,mmf = 1500, rux = 0,  belu = 0, pred = 0),
      list(name = "CsA → Ruxolitinib",      csa = 3,   tac = 0,   mmf = 1500, rux = 10, belu = 0, pred = 0),
      list(name = "CsA → Belumosudil",      csa = 3,   tac = 0,   mmf = 1500, rux = 0,  belu = 200, pred = 0)
    )
    withProgress(message = "Running scenario comparison...", {
      bind_rows(lapply(scenarios, function(sc) {
        run_sim(weight = 70, csa_dose = sc$csa, tac_dose = sc$tac, mmf_dose = sc$mmf,
                rux_dose = sc$rux, belu_dose = sc$belu, pred_dose = sc$pred,
                sim_days = 365) %>%
          mutate(scenario = sc$name)
      }))
    })
  })

  output$compare_agvhd <- renderPlotly({
    dat <- compare_data() %>% filter(time %% 24 == 0)
    dat %>% plot_ly(x = ~Day, y = ~aGvHD_Grade, color = ~scenario, type = "scatter", mode = "lines") %>%
      layout(title = "aGvHD Grade Comparison", yaxis = list(title = "Grade", range = c(0,3)))
  })

  output$compare_cgvhd <- renderPlotly({
    dat <- compare_data() %>% filter(time %% 24 == 0)
    dat %>% plot_ly(x = ~Day, y = ~cGvHD_Score, color = ~scenario, type = "scatter", mode = "lines") %>%
      layout(title = "cGvHD Score Comparison", yaxis = list(title = "Score", range = c(0,3)))
  })

  output$compare_fibr <- renderPlotly({
    dat <- compare_data() %>% filter(time %% 24 == 0)
    dat %>% plot_ly(x = ~Day, y = ~Fibrosis, color = ~scenario, type = "scatter", mode = "lines") %>%
      layout(title = "Fibrosis Index", yaxis = list(title = "Fibrosis (0-1)"))
  })

  output$compare_table <- DT::renderDataTable({
    dat <- compare_data() %>% filter(time %% 24 == 0)
    dat %>% group_by(scenario) %>%
      summarise(
        `Max aGvHD` = round(max(aGvHD_Grade), 2),
        `End cGvHD` = round(tail(cGvHD_Score, 1), 2),
        `End Fibrosis` = round(tail(Fibrosis, 1), 3),
        `Max Th17` = round(max(Th17), 2),
        `Min Treg` = round(min(Treg), 2),
        .groups = "drop"
      )
  }, options = list(pageLength = 6, dom = "t"))

  # Tab 7: Biomarkers
  output$biomarker_plot <- renderPlotly({
    dat <- sim_data() %>% filter(time %% 24 == 0)
    plot_ly(dat, x = ~Day) %>%
      add_lines(y = ~ST2_bio,   name = "ST2 (ng/mL)", line = list(color = "#E74C3C")) %>%
      add_lines(y = ~REG3a_bio, name = "REG3α (ng/mL)",line=list(color="#F39C12")) %>%
      add_lines(y = ~TNFR1_bio, name = "sTNFR1 (ng/mL)",line=list(color="#9B59B6")) %>%
      add_segments(x = 0, xend = max(dat$Day), y = 33, yend = 33,
                   line = list(dash = "dot", color = "#E74C3C"), name = "ST2 Cutoff (33 ng/mL)") %>%
      layout(title = "GvHD Biomarker Panel",
             xaxis = list(title = "Days Post-HSCT"),
             yaxis = list(title = "Biomarker Concentration (ng/mL)"))
  })

  output$biomarker_table <- DT::renderDataTable({
    tibble(
      Biomarker = c("ST2 (sST2)", "REG3α", "sTNFR1", "sIL2Rα (CD25)", "Elafin", "CXCL9"),
      Normal    = c("<33 ng/mL", "<10 ng/mL", "<2 ng/mL", "<500 U/mL", "<10 ng/mL", "<100 pg/mL"),
      Cutoff    = c(">33 predicts grade 3-4 GI", ">23 poor GI outcome",
                    ">2.5 ng/mL: aGvHD", ">1000 U/mL: active aGvHD",
                    "Skin GvHD marker", "T cell homing, aGvHD"),
      `Tissue Source` = c("GI epithelium", "Paneth cell", "TNF-α stimulated cells",
                          "Activated T cells", "Keratinocytes", "IFN-γ induced")
    )
  }, options = list(pageLength = 6, dom = "t"))
}

## ============================================================
## Run App
## ============================================================

shinyApp(ui = ui, server = server)
