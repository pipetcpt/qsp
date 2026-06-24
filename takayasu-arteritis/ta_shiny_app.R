## ============================================================
## Takayasu Arteritis QSP — Interactive Shiny Dashboard
## Author : QSP Disease Model Library (CCR)
## Date   : 2026-06-19
## ============================================================
## Tabs:
##   1. Patient Profile & Disease Overview
##   2. Pharmacokinetics (Drug Concentration)
##   3. PD Key Indicators (IL-6, TNF, T cells)
##   4. Clinical Endpoints (NIH, ITAS, CRP, ESR)
##   5. Scenario Comparison (multi-drug)
##   6. Biomarkers (PET-CT, MRI-VWT, Serum markers)
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
## mrgsolve model code (inline)
## ============================================================
ta_code <- '
$PARAM
  ka_PRED=1.50, Vc_PRED=28.0, Vp_PRED=56.0, CLd_PRED=40.0, CL_PRED=18.0,
  ke0_PRED=0.25, F_PRED=0.82,
  ka_TCZ=0.0083, Vc_TCZ=3.5, Vp_TCZ=2.0, CLd_TCZ=0.30, CL_TCZ=0.18,
  CLMM_TCZ=0.08, KM_TCZ=1.20, F_TCZ=0.80,
  ka_MTX=0.90, Vc_MTX=18.0, CL_MTX=5.0, kpg_MTX=0.12, kdpg_MTX=0.018, F_MTX=0.70,
  Vc_IFX=3.0, Vp_IFX=1.8, CLd_IFX=0.22, CL_IFX=0.16,
  ksyn_IL6=0.45, kdeg_IL6=0.25, kon_IL6R=0.001, koff_IL6R=0.002,
  ksyn_sR=3.0, kdeg_sR=0.015,
  ksyn_TNF=0.20, kdeg_TNF=0.30,
  ksyn_Th1=0.10, kdeg_Th1=0.02, inh_Treg_Th1=0.025,
  ksyn_Th17=0.07, kdeg_Th17=0.02, inh_Treg_Th17=0.030,
  ksyn_Treg=0.02, kdeg_Treg=0.018,
  ksyn_VWI=0.008, kdeg_VWI=0.004,
  kprog_ST=0.0005, kreg_ST=0.0001,
  ksyn_CRP=0.60, kdeg_CRP=0.025,
  ksyn_PET=0.0015, kdeg_PET=0.003,
  kVWT=0.002,
  Emax_PRED=0.85, EC50_PRED=0.10, hill_PRED=1.5,
  Emax_TCZ=0.95, EC50_TCZ=0.50, hill_TCZ=1.2,
  Emax_MTX=0.55, EC50_MTX=0.05, hill_MTX=1.0,
  Emax_IFX=0.90, EC50_IFX=0.80, hill_IFX=1.3,
  amp_IL6_Th1=0.006, amp_IL6_Th17=0.004, amp_IL6_TNF=0.05,
  amp_VWI_IL6=0.20, amp_VWI_TNF=0.15, amp_VWI_Th17=0.08,
  amp_ST_VWI=0.40, ESR_base=20.0, kESR=2.5, CRP_base=5.0

$CMT PRED_GUT PRED_C PRED_P PRED_EFF
     TCZ_SC TCZ_C TCZ_P
     MTX_GUT MTX_C MTX_PG
     IFX_C IFX_P
     IL6 sIL6R IL6_cmplx TNF TH1 TH17 TREG
     VWI ST CRP PET VWT

$INIT
  IL6=8, sIL6R=200, IL6_cmplx=0, TNF=2, TH1=50, TH17=35, TREG=8,
  VWI=6, ST=0, CRP=45, PET=3.5, VWT=4.5

$ODE
  double Inh_PRED = Emax_PRED*pow(PRED_EFF,hill_PRED)/(pow(EC50_PRED,hill_PRED)+pow(PRED_EFF,hill_PRED));
  double Occ_TCZ  = Emax_TCZ*pow(TCZ_C,hill_TCZ)/(pow(EC50_TCZ,hill_TCZ)+pow(TCZ_C,hill_TCZ));
  double Inh_MTX  = Emax_MTX*pow(MTX_PG,hill_MTX)/(pow(EC50_MTX,hill_MTX)+pow(MTX_PG,hill_MTX));
  double Inh_IFX  = Emax_IFX*pow(IFX_C,hill_IFX)/(pow(EC50_IFX,hill_IFX)+pow(IFX_C,hill_IFX));

  dxdt_PRED_GUT = -ka_PRED*PRED_GUT;
  dxdt_PRED_C   = ka_PRED*F_PRED*PRED_GUT/Vc_PRED-(CL_PRED/Vc_PRED)*PRED_C-(CLd_PRED/Vc_PRED)*PRED_C+(CLd_PRED/Vp_PRED)*PRED_P;
  dxdt_PRED_P   = (CLd_PRED/Vc_PRED)*PRED_C-(CLd_PRED/Vp_PRED)*PRED_P;
  dxdt_PRED_EFF = ke0_PRED*(PRED_C-PRED_EFF);
  dxdt_TCZ_SC   = -ka_TCZ*TCZ_SC;
  double CL_TCZ_tot=CL_TCZ+CLMM_TCZ*TCZ_C/(KM_TCZ+TCZ_C);
  dxdt_TCZ_C    = ka_TCZ*F_TCZ*TCZ_SC/Vc_TCZ-(CL_TCZ_tot/Vc_TCZ)*TCZ_C-(CLd_TCZ/Vc_TCZ)*TCZ_C+(CLd_TCZ/Vp_TCZ)*TCZ_P;
  dxdt_TCZ_P    = (CLd_TCZ/Vc_TCZ)*TCZ_C-(CLd_TCZ/Vp_TCZ)*TCZ_P;
  dxdt_MTX_GUT  = -ka_MTX*MTX_GUT;
  dxdt_MTX_C    = ka_MTX*F_MTX*MTX_GUT/Vc_MTX-(CL_MTX/Vc_MTX)*MTX_C-kpg_MTX*MTX_C;
  dxdt_MTX_PG   = kpg_MTX*MTX_C-kdpg_MTX*MTX_PG;
  dxdt_IFX_C    = -(CL_IFX/Vc_IFX)*IFX_C-(CLd_IFX/Vc_IFX)*IFX_C+(CLd_IFX/Vp_IFX)*IFX_P;
  dxdt_IFX_P    = (CLd_IFX/Vc_IFX)*IFX_C-(CLd_IFX/Vp_IFX)*IFX_P;

  double IL6_syn=ksyn_IL6+amp_IL6_Th1*TH1+amp_IL6_Th17*TH17+amp_IL6_TNF*TNF;
  dxdt_IL6      = IL6_syn*(1-Inh_PRED)-kdeg_IL6*(1+Occ_TCZ*0.2)*IL6-kon_IL6R*IL6*sIL6R+koff_IL6R*IL6_cmplx;
  dxdt_sIL6R    = ksyn_sR*(1+2.5*Occ_TCZ)-kdeg_sR*sIL6R-kon_IL6R*IL6*sIL6R+koff_IL6R*IL6_cmplx;
  dxdt_IL6_cmplx= kon_IL6R*IL6*sIL6R-koff_IL6R*IL6_cmplx-kdeg_IL6*IL6_cmplx;
  dxdt_TNF      = ksyn_TNF+0.03*TH1-kdeg_TNF*TNF*(1-Inh_IFX)*(1-0.6*Inh_PRED);
  dxdt_TH1      = ksyn_Th1*(1+0.05*IL6)-kdeg_Th1*TH1*(1-Inh_MTX)*(1-Inh_PRED*0.5)-inh_Treg_Th1*TREG*TH1;
  dxdt_TH17     = ksyn_Th17*(1+0.04*IL6)*(1-Occ_TCZ*0.8)-kdeg_Th17*TH17*(1-Inh_MTX*0.7)*(1-Inh_PRED*0.4)-inh_Treg_Th17*TREG*TH17;
  dxdt_TREG     = ksyn_Treg*(1+0.3*Inh_PRED+0.2*Occ_TCZ)-kdeg_Treg*TREG;
  double VWI_drive=amp_VWI_IL6*IL6+amp_VWI_TNF*TNF+amp_VWI_Th17*TH17;
  double Drug_inh_VWI=1-(1-Inh_PRED)*(1-Occ_TCZ*0.9)*(1-Inh_IFX*0.7)*(1-Inh_MTX*0.3);
  dxdt_VWI      = ksyn_VWI*VWI_drive*(1-Drug_inh_VWI)-kdeg_VWI*VWI;
  dxdt_ST       = kprog_ST*amp_ST_VWI*VWI-kreg_ST*ST;
  dxdt_CRP      = ksyn_CRP*IL6*(1-Occ_TCZ*0.95)-kdeg_CRP*CRP;
  dxdt_PET      = ksyn_PET*VWI-kdeg_PET*PET;
  dxdt_VWT      = kVWT*VWI-0.001*VWT;

$TABLE
  double ESR_now=ESR_base+kESR*(CRP-CRP_base);
  double NIH_SCORE=2*(CRP>20?1:CRP/20)+2*(ESR_now>40?1:ESR_now/40)+3*(VWI/10)+4*(ST/50)+3*(PET/4);
  if(NIH_SCORE>20) NIH_SCORE=20;
  double ITAS=1.5*(CRP>10?1:0)+1.5*(VWI>5?1:0)+3.0*(ST>20?1:0)+2.0*(PET>2.5?1:0);
  capture ESR_now NIH_SCORE ITAS
'
mod_shiny <- mcode("TA_shiny", ta_code, quiet = TRUE)

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "red",

  dashboardHeader(
    title = "Takayasu Arteritis QSP",
    titleWidth = 320
  ),

  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("Patient Profile", tabName = "profile",
               icon = icon("user-circle")),
      menuItem("Pharmacokinetics", tabName = "pk",
               icon = icon("flask")),
      menuItem("PD Key Indicators", tabName = "pd",
               icon = icon("dna")),
      menuItem("Clinical Endpoints", tabName = "endpoints",
               icon = icon("heartbeat")),
      menuItem("Scenario Comparison", tabName = "scenarios",
               icon = icon("balance-scale")),
      menuItem("Biomarkers", tabName = "biomarkers",
               icon = icon("microscope"))
    ),
    hr(),
    tags$div(style = "padding: 10px;",
      tags$h4("Global Parameters", style = "color:#F48FB1;"),

      sliderInput("wt", "Body Weight (kg)", 40, 100, 65, step = 5),
      sliderInput("sim_days", "Simulation Duration (days)", 90, 730, 365, step = 30),

      tags$h4("Disease Activity at Baseline", style = "color:#F48FB1;"),
      sliderInput("VWI_init", "VWI Score (0-10)", 0, 10, 6, step = 0.5),
      sliderInput("IL6_init", "Serum IL-6 (pg/mL)", 1, 50, 8, step = 1),
      sliderInput("CRP_init", "CRP (mg/L)", 5, 150, 45, step = 5),

      tags$h4("Drug Selection", style = "color:#F48FB1;"),
      checkboxInput("use_pred", "Prednisone", value = TRUE),
      conditionalPanel("input.use_pred",
        sliderInput("dose_pred", "Prednisone dose (mg/day)", 5, 80, 65, step = 5),
        sliderInput("taper_pred", "Taper start (day)", 30, 180, 60, step = 10),
        sliderInput("taper_to_pred", "Taper target (mg/day)", 2.5, 20, 10, step = 2.5)
      ),
      checkboxInput("use_tcz", "Tocilizumab 162 mg SC q2w", value = FALSE),
      checkboxInput("use_mtx", "Methotrexate 15 mg/week PO", value = FALSE),
      checkboxInput("use_ifx", "Infliximab 5 mg/kg IV", value = FALSE),
      actionButton("run_sim", "Run Simulation", icon = icon("play"),
                   class = "btn-warning", width = "100%")
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper, .right-side { background-color: #FAFAFA; }
      .box { border-top: 3px solid #E65100; }
      .value-box .inner h3 { font-size: 26px; }
    "))),

    tabItems(

      ## ---- TAB 1: Patient Profile ----
      tabItem(tabName = "profile",
        fluidRow(
          valueBoxOutput("vb_nih", width = 3),
          valueBoxOutput("vb_crp", width = 3),
          valueBoxOutput("vb_il6", width = 3),
          valueBoxOutput("vb_pet", width = 3)
        ),
        fluidRow(
          box(width = 8, title = "Disease Overview: Takayasu Arteritis",
              solidHeader = TRUE, status = "danger",
              tags$div(
                tags$h4("Classification (ACR/EULAR 2022)"),
                tags$p("Takayasu arteritis (TA) is a granulomatous large-vessel vasculitis
                primarily affecting the aorta and its major branches. Onset typically
                occurs in women <50 years. The disease progresses through two phases:
                (1) inflammatory (active vasculitis with systemic symptoms) and
                (2) fibrotic/occlusive (stenosis, aneurysm formation)."),
                tags$h4("Angiographic Pattern (Numano Classification)"),
                tags$ul(
                  tags$li("Type I: Aortic arch branches only"),
                  tags$li("Type IIa: Ascending aorta + arch"),
                  tags$li("Type IIb: Ascending + arch + descending thoracic aorta"),
                  tags$li("Type III: Descending thoracic + abdominal aorta"),
                  tags$li("Type IV: Abdominal aorta only"),
                  tags$li("Type V: Combined (most common worldwide)")
                ),
                tags$h4("Key Pathogenic Drivers Modeled"),
                tags$ul(
                  tags$li("IL-6 / STAT3 axis → Th17 polarization, CRP, vessel wall inflammation"),
                  tags$li("TNF-α → NF-κB → ICAM-1/VCAM-1 upregulation → leukocyte adhesion"),
                  tags$li("CD8+ CTL & NK cells → perforin/granzyme-B → smooth muscle injury"),
                  tags$li("Intimal hyperplasia → luminal stenosis → ischemic complications"),
                  tags$li("Vasa vasorum neovascularization → FDG-PET avidity")
                )
              )
          ),
          box(width = 4, title = "Baseline Summary",
              solidHeader = TRUE, status = "warning",
              tableOutput("baseline_table")
          )
        ),
        fluidRow(
          box(width = 12, title = "Disease Natural History — NIH Score Trajectory",
              solidHeader = TRUE, status = "danger",
              plotlyOutput("natural_history_plot", height = "320px")
          )
        )
      ),

      ## ---- TAB 2: Pharmacokinetics ----
      tabItem(tabName = "pk",
        fluidRow(
          box(width = 6, title = "Prednisone/Prednisolone Plasma Concentration",
              solidHeader = TRUE, status = "warning",
              plotlyOutput("pk_pred_plot", height = "300px")),
          box(width = 6, title = "Tocilizumab Plasma Concentration (SC q2w)",
              solidHeader = TRUE, status = "warning",
              plotlyOutput("pk_tcz_plot", height = "300px"))
        ),
        fluidRow(
          box(width = 6, title = "Methotrexate Plasma & Intracellular Polyglutamates",
              solidHeader = TRUE, status = "warning",
              plotlyOutput("pk_mtx_plot", height = "300px")),
          box(width = 6, title = "Infliximab Plasma Concentration (5 mg/kg q6w)",
              solidHeader = TRUE, status = "warning",
              plotlyOutput("pk_ifx_plot", height = "300px"))
        ),
        fluidRow(
          box(width = 12, title = "PK Summary Table (Cmax, Cmin, AUC)",
              solidHeader = TRUE, status = "primary",
              DTOutput("pk_summary_table"))
        )
      ),

      ## ---- TAB 3: PD Key Indicators ----
      tabItem(tabName = "pd",
        fluidRow(
          box(width = 6, title = "Serum IL-6 Dynamics",
              solidHeader = TRUE, status = "danger",
              plotlyOutput("il6_plot", height = "280px"),
              tags$small("Note: TCZ causes paradoxical serum IL-6 rise (blocked signaling).")
          ),
          box(width = 6, title = "TNF-α Dynamics",
              solidHeader = TRUE, status = "danger",
              plotlyOutput("tnf_plot", height = "280px"))
        ),
        fluidRow(
          box(width = 6, title = "T Cell Subsets (Th1 / Th17 / Treg)",
              solidHeader = TRUE, status = "success",
              plotlyOutput("tcell_plot", height = "280px")),
          box(width = 6, title = "Vessel Wall Inflammation Index (VWI 0-10)",
              solidHeader = TRUE, status = "danger",
              plotlyOutput("vwi_plot", height = "280px"))
        ),
        fluidRow(
          box(width = 12, title = "Drug Effect Occupancy / Inhibition",
              solidHeader = TRUE, status = "warning",
              plotlyOutput("drug_effect_plot", height = "260px"))
        )
      ),

      ## ---- TAB 4: Clinical Endpoints ----
      tabItem(tabName = "endpoints",
        fluidRow(
          box(width = 6, title = "NIH Disease Activity Score (0-20)",
              solidHeader = TRUE, status = "danger",
              plotlyOutput("nih_plot", height = "280px"),
              tags$small("Remission: NIH < 4, Low activity: 4-8, Active: >8")),
          box(width = 6, title = "ITAS 2010 Score",
              solidHeader = TRUE, status = "danger",
              plotlyOutput("itas_plot", height = "280px"))
        ),
        fluidRow(
          box(width = 6, title = "CRP (mg/L)",
              solidHeader = TRUE, status = "warning",
              plotlyOutput("crp_plot", height = "280px")),
          box(width = 6, title = "ESR (mm/hr)",
              solidHeader = TRUE, status = "warning",
              plotlyOutput("esr_plot", height = "280px"))
        ),
        fluidRow(
          box(width = 12, title = "Arterial Stenosis Index (%)",
              solidHeader = TRUE, status = "danger",
              plotlyOutput("stenosis_plot", height = "260px"))
        )
      ),

      ## ---- TAB 5: Scenario Comparison ----
      tabItem(tabName = "scenarios",
        fluidRow(
          box(width = 12, title = "Multi-Scenario Comparison: NIH Score",
              solidHeader = TRUE, status = "danger",
              plotlyOutput("scenario_nih", height = "320px"))
        ),
        fluidRow(
          box(width = 6, title = "Scenario Comparison: CRP",
              solidHeader = TRUE, status = "warning",
              plotlyOutput("scenario_crp", height = "280px")),
          box(width = 6, title = "Scenario Comparison: VWI",
              solidHeader = TRUE, status = "danger",
              plotlyOutput("scenario_vwi", height = "280px"))
        ),
        fluidRow(
          box(width = 12, title = "Year-1 Outcome Table",
              solidHeader = TRUE, status = "primary",
              DTOutput("scenario_table"))
        )
      ),

      ## ---- TAB 6: Biomarkers ----
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(width = 6, title = "PET-CT FDG SUVmax (Vascular Activity)",
              solidHeader = TRUE, status = "primary",
              plotlyOutput("pet_plot", height = "280px"),
              tags$small("PET positivity threshold: SUVmax > 2.5 (liver reference method)")),
          box(width = 6, title = "MRI Vessel Wall Thickness (mm)",
              solidHeader = TRUE, status = "primary",
              plotlyOutput("vwt_plot", height = "280px"),
              tags$small("Normal: < 2 mm; Thickened: ≥ 3 mm"))
        ),
        fluidRow(
          box(width = 6, title = "Serum IL-6 vs. CRP Correlation",
              solidHeader = TRUE, status = "danger",
              plotlyOutput("il6_crp_scatter", height = "280px")),
          box(width = 6, title = "Soluble IL-6R (ng/mL) — TCZ pharmacodynamics",
              solidHeader = TRUE, status = "info",
              plotlyOutput("sil6r_plot", height = "280px"),
              tags$small("TCZ paradox: serum IL-6 ↑ and sIL-6R ↑ after blockade (non-functional complex)"))
        ),
        fluidRow(
          box(width = 12, title = "Biomarker Summary at Key Timepoints",
              solidHeader = TRUE, status = "primary",
              DTOutput("biomarker_table"))
        )
      )
    )
  )
)

## ============================================================
## SERVER
## ============================================================
server <- function(input, output, session) {

  ## ---- Reactive: run simulation ----
  sim_data <- eventReactive(input$run_sim, {
    sim_end <- input$sim_days * 24
    dt <- 4

    ## Build dosing events
    ev_list <- ev(time = 0, amt = 0, cmt = "PRED_GUT")

    if (input$use_pred) {
      times_init  <- seq(0, input$taper_pred * 24, by = 24)
      times_taper <- seq((input$taper_pred + 7) * 24, sim_end, by = 24)
      dt_  <- max(1, length(times_taper))
      dose_taper <- seq(input$dose_pred, input$taper_to_pred, length.out = dt_)
      ev_list <- c(
        ev(time = times_init,  amt = input$dose_pred, cmt = "PRED_GUT"),
        ev(time = times_taper, amt = dose_taper,      cmt = "PRED_GUT")
      )
    }

    if (input$use_tcz) {
      times_tcz <- seq(0, sim_end, by = 14 * 24)
      ev_list <- c(ev_list, ev(time = times_tcz, amt = 162, cmt = "TCZ_SC"))
    }

    if (input$use_mtx) {
      times_mtx <- seq(0, sim_end, by = 7 * 24)
      ev_list <- c(ev_list, ev(time = times_mtx, amt = 15, cmt = "MTX_GUT"))
    }

    if (input$use_ifx) {
      dose_ifx <- 5 * input$wt
      ind_t <- c(0, 2, 6) * 7 * 24
      mnt_t <- seq(12 * 7 * 24, sim_end, by = 6 * 7 * 24)
      ev_list <- c(ev_list, ev(time = c(ind_t, mnt_t), amt = dose_ifx, cmt = "IFX_C"))
    }

    mod_shiny %>%
      init(IL6 = input$IL6_init, CRP = input$CRP_init, VWI = input$VWI_init,
           TNF = 2, TH1 = 50, TH17 = 35, TREG = 8, PET = 3.5, VWT = 4.5) %>%
      ev(ev_list) %>%
      mrgsim(end = sim_end, delta = dt) %>%
      as_tibble() %>%
      mutate(time_days = time / 24)
  }, ignoreNULL = FALSE)

  ## ---- Run all 5 scenarios ----
  all_scenarios <- eventReactive(input$run_sim, {
    sim_end <- input$sim_days * 24
    dt <- 4

    scen_list <- list(
      "No Treatment"         = ev(time = 0, amt = 0, cmt = "PRED_GUT"),
      "Prednisone"           = ev(time = seq(0, sim_end, by=24), amt = 65, cmt = "PRED_GUT"),
      "Pred + MTX"           = c(ev(time=seq(0,sim_end,by=24),amt=65,cmt="PRED_GUT"),
                                  ev(time=seq(0,sim_end,by=7*24),amt=15,cmt="MTX_GUT")),
      "Pred + TCZ"           = c(ev(time=seq(0,sim_end,by=24),amt=65,cmt="PRED_GUT"),
                                  ev(time=seq(0,sim_end,by=14*24),amt=162,cmt="TCZ_SC")),
      "Pred + Infliximab"    = c(ev(time=seq(0,sim_end,by=24),amt=65,cmt="PRED_GUT"),
                                  ev(time=c(0,14,42,84,126,168,210)*24,amt=5*65,cmt="IFX_C"))
    )

    bind_rows(lapply(names(scen_list), function(nm) {
      mod_shiny %>%
        init(IL6=input$IL6_init, CRP=input$CRP_init, VWI=input$VWI_init,
             TNF=2, TH1=50, TH17=35, TREG=8, PET=3.5, VWT=4.5) %>%
        ev(scen_list[[nm]]) %>%
        mrgsim(end=sim_end, delta=dt) %>%
        as_tibble() %>%
        mutate(time_days=time/24, Scenario=nm)
    }))
  }, ignoreNULL = FALSE)

  ## ---- Value Boxes ----
  get_last <- function(col) {
    d <- sim_data()
    if (is.null(d) || nrow(d) == 0) return(NA)
    tail(d[[col]], 1)
  }

  output$vb_nih <- renderValueBox({
    v <- round(get_last("NIH_SCORE"), 1)
    col <- if (!is.na(v) && v < 4) "green" else if (!is.na(v) && v < 8) "yellow" else "red"
    valueBox(value = v, subtitle = "NIH Activity Score", icon = icon("stethoscope"), color = col)
  })
  output$vb_crp <- renderValueBox({
    v <- round(get_last("CRP"), 1)
    col <- if (!is.na(v) && v < 10) "green" else "orange"
    valueBox(value = paste0(v, " mg/L"), subtitle = "CRP", icon = icon("vials"), color = col)
  })
  output$vb_il6 <- renderValueBox({
    v <- round(get_last("IL6"), 2)
    valueBox(value = paste0(v, " pg/mL"), subtitle = "Serum IL-6", icon = icon("atom"), color = "red")
  })
  output$vb_pet <- renderValueBox({
    v <- round(get_last("PET"), 2)
    col <- if (!is.na(v) && v < 2.5) "green" else "red"
    valueBox(value = v, subtitle = "PET-CT SUVmax", icon = icon("radiation"), color = col)
  })

  ## ---- Baseline Table ----
  output$baseline_table <- renderTable({
    data.frame(
      Parameter = c("Body Weight", "Sim. Duration", "Initial IL-6",
                    "Initial CRP", "Initial VWI", "Prednisone",
                    "Tocilizumab", "Methotrexate", "Infliximab"),
      Value     = c(paste0(input$wt, " kg"),
                    paste0(input$sim_days, " days"),
                    paste0(input$IL6_init, " pg/mL"),
                    paste0(input$CRP_init, " mg/L"),
                    as.character(input$VWI_init),
                    if(input$use_pred) paste0(input$dose_pred," mg/d") else "OFF",
                    if(input$use_tcz) "162 mg SC q2w" else "OFF",
                    if(input$use_mtx) "15 mg/week" else "OFF",
                    if(input$use_ifx) "5 mg/kg IV" else "OFF")
    )
  })

  ## ---- Natural History ----
  output$natural_history_plot <- renderPlotly({
    d <- all_scenarios() %>% filter(Scenario == "No Treatment")
    p <- ggplot(d, aes(time_days, NIH_SCORE)) +
      geom_line(color = "#D32F2F", linewidth = 1.2) +
      geom_hline(yintercept = 4, linetype = "dashed") +
      labs(x = "Time (days)", y = "NIH Score", title = "Natural History — No Treatment") +
      theme_bw()
    ggplotly(p)
  })

  ## ---- PK Plots ----
  pk_plot <- function(col, title, color, ylab = "Conc. (mg/L)") {
    renderPlotly({
      d <- sim_data()
      p <- ggplot(d, aes(time_days, .data[[col]])) +
        geom_line(color = color, linewidth = 1) +
        labs(title = title, x = "Time (days)", y = ylab) + theme_bw()
      ggplotly(p)
    })
  }
  output$pk_pred_plot <- pk_plot("PRED_C", "Prednisolone Plasma", "#FF6F00")
  output$pk_tcz_plot  <- pk_plot("TCZ_C",  "Tocilizumab Plasma",  "#1565C0")
  output$pk_mtx_plot  <- pk_plot("MTX_C",  "MTX Plasma",          "#388E3C", "Conc. (µmol/L)")
  output$pk_ifx_plot  <- pk_plot("IFX_C",  "Infliximab Plasma",   "#6A1B9A")

  output$pk_summary_table <- renderDT({
    d <- sim_data()
    drugs <- list(
      list(col="PRED_C", name="Prednisolone"),
      list(col="TCZ_C",  name="Tocilizumab"),
      list(col="MTX_C",  name="Methotrexate"),
      list(col="IFX_C",  name="Infliximab")
    )
    rows <- lapply(drugs, function(dr) {
      vals <- d[[dr$col]]
      data.frame(Drug=dr$name, Cmax=round(max(vals),3),
                 Cmin=round(min(vals[vals>0]),6),
                 AUC_approx=round(sum(diff(d$time/24)*head(vals,-1)),2))
    })
    bind_rows(rows)
  }, options = list(pageLength = 5))

  ## ---- PD Plots ----
  pd_plot <- function(col, title, color) {
    renderPlotly({
      d <- sim_data()
      p <- ggplot(d, aes(time_days, .data[[col]])) +
        geom_line(color = color, linewidth = 1) +
        labs(title = title, x = "Time (days)", y = col) + theme_bw()
      ggplotly(p)
    })
  }
  output$il6_plot  <- pd_plot("IL6",  "Serum IL-6 (pg/mL)",  "#E91E63")
  output$tnf_plot  <- pd_plot("TNF",  "TNF-α (pg/mL)",       "#C62828")
  output$vwi_plot  <- pd_plot("VWI",  "VWI Score (0-10)",    "#6A1B9A")

  output$tcell_plot <- renderPlotly({
    d <- sim_data() %>%
      select(time_days, TH1, TH17, TREG) %>%
      pivot_longer(c(TH1, TH17, TREG), names_to = "Cell", values_to = "Count")
    p <- ggplot(d, aes(time_days, Count, color = Cell)) +
      geom_line(linewidth = 0.9) +
      scale_color_manual(values = c(TH1="#D32F2F", TH17="#1565C0", TREG="#2E7D32")) +
      labs(title = "T Cell Subsets", x = "Time (days)", y = "Count (cells/µL)") + theme_bw()
    ggplotly(p)
  })

  output$drug_effect_plot <- renderPlotly({
    d <- sim_data() %>%
      mutate(
        pred_eff = Emax_PRED * PRED_EFF^1.5 / (0.10^1.5 + PRED_EFF^1.5),
        tcz_eff  = 0.95 * TCZ_C^1.2  / (0.50^1.2  + TCZ_C^1.2),
        mtx_eff  = 0.55 * MTX_PG^1.0 / (0.05^1.0  + MTX_PG^1.0),
        ifx_eff  = 0.90 * IFX_C^1.3  / (0.80^1.3  + IFX_C^1.3)
      ) %>%
      select(time_days, pred_eff, tcz_eff, mtx_eff, ifx_eff) %>%
      pivot_longer(-time_days, names_to = "Drug", values_to = "Effect")
    p <- ggplot(d, aes(time_days, Effect, color = Drug)) +
      geom_line(linewidth = 0.8) +
      scale_color_manual(values=c(pred_eff="#FF6F00",tcz_eff="#1565C0",
                                   mtx_eff="#388E3C",ifx_eff="#6A1B9A")) +
      scale_y_continuous(labels=scales::percent_format()) +
      labs(title="Drug Effect (% Emax)", x="Time (days)", y="Fractional Emax") + theme_bw()
    ggplotly(p)
  })

  ## ---- Clinical Endpoint Plots ----
  output$nih_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_days, NIH_SCORE)) +
      geom_line(color="#D32F2F", linewidth=1.2) +
      geom_hline(yintercept=4, linetype="dashed", color="gray40") +
      labs(title="NIH Disease Activity Score", x="Time (days)", y="NIH Score (0-20)") +
      theme_bw()
    ggplotly(p)
  })
  output$itas_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_days, ITAS)) +
      geom_line(color="#880E4F", linewidth=1.2) +
      labs(title="ITAS 2010 Score", x="Time (days)", y="ITAS Score") + theme_bw()
    ggplotly(p)
  })
  output$crp_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_days, CRP)) +
      geom_line(color="#E65100", linewidth=1) +
      geom_hline(yintercept=10, linetype="dashed") +
      labs(title="CRP (mg/L)", x="Time (days)", y="CRP (mg/L)") + theme_bw()
    ggplotly(p)
  })
  output$esr_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_days, ESR_now)) +
      geom_line(color="#BF360C", linewidth=1) +
      geom_hline(yintercept=40, linetype="dashed") +
      labs(title="ESR (mm/hr)", x="Time (days)", y="ESR (mm/hr)") + theme_bw()
    ggplotly(p)
  })
  output$stenosis_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_days, ST)) +
      geom_line(color="#B71C1C", linewidth=1) +
      labs(title="Arterial Stenosis Index (%)", x="Time (days)", y="Stenosis (%)") + theme_bw()
    ggplotly(p)
  })

  ## ---- Scenario Comparison ----
  scen_colors <- c("No Treatment"="#D32F2F","Prednisone"="#F57C00",
                   "Pred + MTX"="#388E3C","Pred + TCZ"="#1565C0",
                   "Pred + Infliximab"="#6A1B9A")

  scenario_plot <- function(col, title) {
    renderPlotly({
      d <- all_scenarios()
      p <- ggplot(d, aes(time_days, .data[[col]], color = Scenario)) +
        geom_line(linewidth = 0.9) +
        scale_color_manual(values = scen_colors) +
        labs(title = title, x = "Time (days)", y = col) + theme_bw() +
        theme(legend.position = "bottom")
      ggplotly(p)
    })
  }
  output$scenario_nih <- scenario_plot("NIH_SCORE", "NIH Score by Scenario")
  output$scenario_crp <- scenario_plot("CRP",        "CRP by Scenario")
  output$scenario_vwi <- scenario_plot("VWI",        "VWI by Scenario")

  output$scenario_table <- renderDT({
    all_scenarios() %>%
      filter(abs(time_days - input$sim_days) < 0.5) %>%
      group_by(Scenario) %>% slice(1) %>%
      select(Scenario, IL6, CRP, VWI, ST, PET, NIH_SCORE, ITAS) %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
  })

  ## ---- Biomarker Plots ----
  output$pet_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_days, PET)) +
      geom_line(color="#1565C0", linewidth=1) +
      geom_hline(yintercept=2.5, linetype="dashed") +
      labs(title="FDG-PET SUVmax", x="Time (days)", y="SUVmax") + theme_bw()
    ggplotly(p)
  })
  output$vwt_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_days, VWT)) +
      geom_line(color="#4A148C", linewidth=1) +
      geom_hline(yintercept=3, linetype="dashed") +
      labs(title="MRI Vessel Wall Thickness (mm)", x="Time (days)", y="VWT (mm)") + theme_bw()
    ggplotly(p)
  })
  output$il6_crp_scatter <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(IL6, CRP, color=time_days)) +
      geom_path(linewidth=0.7) +
      geom_point(size=0.6, alpha=0.5) +
      scale_color_gradient(low="#FFCCBC", high="#B71C1C") +
      labs(title="IL-6 vs. CRP (Phase Plot)", x="IL-6 (pg/mL)", y="CRP (mg/L)") + theme_bw()
    ggplotly(p)
  })
  output$sil6r_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_days, sIL6R)) +
      geom_line(color="#00838F", linewidth=1) +
      labs(title="Soluble IL-6R (ng/mL)", x="Time (days)", y="sIL-6R (ng/mL)") + theme_bw()
    ggplotly(p)
  })
  output$biomarker_table <- renderDT({
    d <- sim_data()
    key_days <- c(0, 30, 90, 180, 365)
    bind_rows(lapply(key_days, function(td) {
      row <- d %>% filter(abs(time_days - td) == min(abs(time_days - td))) %>% slice(1)
      data.frame(Day=td, IL6=round(row$IL6,2), CRP=round(row$CRP,1),
                 VWI=round(row$VWI,2), PET=round(row$PET,2),
                 VWT=round(row$VWT,2), NIH=round(row$NIH_SCORE,1))
    }))
  })
}

## ============================================================
## Run App
## ============================================================
shinyApp(ui = ui, server = server)
