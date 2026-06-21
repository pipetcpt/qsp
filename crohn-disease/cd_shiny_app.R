# =============================================================================
# Crohn's Disease QSP — Interactive Shiny Dashboard
# 7 tabs: Patient Profile · PK · Cytokines/PD · Clinical Endpoints ·
#         Scenario Comparison · Biomarker Correlation · Dosing Optimizer
# =============================================================================

library(shiny)
library(bslib)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(patchwork)

# --------------------------------------------------------------------------
# Inline model code (same as cd_mrgsolve_model.R, condensed for Shiny)
# --------------------------------------------------------------------------
cd_code <- '
$PROB Crohn Disease QSP — Shiny Interactive Model

$PARAM
  Vc_IFX=4.0, Vp_IFX=3.5, CL_IFX=0.55, Q_IFX=1.40,
  V_ADA=7.8, CL_ADA=0.35, ka_ADA=0.38, F_ADA=0.64,
  Vc_UST=3.5, Vp_UST=7.3, CL_UST=0.19, Q_UST=0.91,
  Vc_VDZ=3.3, CL_VDZ=0.49,
  V_TGN=1500, CL_TGN=10.5, ka_TGN=0.040, F_AZA=0.47,
  V_PRED=45.0, CL_PRED=7.2,
  V_UPA=85.0, CL_UPA=24.0, ka_UPA=2.5, F_UPA=0.79,
  TNF0=100, kprod_TNF=5.0, kdeg_TNF=0.05,
  Emax_IFX_TNF=0.93, EC50_IFX_TNF=1.0,
  Emax_ADA_TNF=0.90, EC50_ADA_TNF=0.8,
  IL12230=80, kprod_IL1223=4.0, kdeg_IL1223=0.05,
  Emax_UST=0.88, EC50_UST=0.8,
  IL170=60, kprod_IL17=3.0, kdeg_IL17=0.05, kIL23_IL17=0.02,
  Th170=1.0, kprol_Th17=0.06, kdeg_Th17=0.06, kIL23_Th17=0.03, kTreg_Th17=0.02,
  Th10=1.0, kprol_Th1=0.06, kdeg_Th1=0.06, kIL12_Th1=0.025,
  Treg0=1.0, kprol_Treg=0.05, kdeg_Treg=0.05,
  Neut0=1.0, kprod_Neut=0.05, kdeg_Neut=0.05, kIL17_Neut=0.03, kTNF_Neut=0.02,
  MI0=65, kprod_MI=3.25, kdeg_MI=0.05, kTNF_MI=0.015, kIL17_MI=0.010, kTreg_MI=0.010, kNeut_MI=0.008,
  CRP0=25, kprod_CRP=0.096, kdeg_CRP=0.25,
  FC0=800, kprod_FC=2.0, kdeg_FC=0.06, kNeut_FC=25.0, kMI_FC=0.55,
  BMD0=1.0, kform_BMD=0.000195, kresorb_BMD=0.000200, kTNF_BMD=0.000080, kPRED_BMD=0.000150,
  Hgb0=10.5, kprod_Hgb=0.0913, kdeg_Hgb=0.007, kMI_Hgb=0.00195,
  Emax_VDZ_MI=0.60, EC50_VDZ_MI=0.4,
  Emax_TGN_MI=0.45, EC50_TGN_MI=235,
  Emax_PRED=0.72, EC50_PRED=0.05,
  Emax_UPA=0.74, EC50_UPA=0.08

$CMT A1_IFX A2_IFX ADA_dep A_ADA A1_UST A2_UST A_VDZ A_TGN A_PRED A_UPA
     TNF IL1223 IL17 Th17 Th1 Treg Neut MI CRP FC BMD Hgb

$INIT A1_IFX=0,A2_IFX=0,ADA_dep=0,A_ADA=0,A1_UST=0,A2_UST=0,
      A_VDZ=0,A_TGN=0,A_PRED=0,A_UPA=0,
      TNF=100,IL1223=80,IL17=60,Th17=1,Th1=1,Treg=1,Neut=1,
      MI=65,CRP=25,FC=800,BMD=1.0,Hgb=10.5

$ODE
  double C1_IFX = A1_IFX/Vc_IFX;
  double C_ADA  = A_ADA/V_ADA;
  double C1_UST = A1_UST/Vc_UST;
  double C_VDZ  = A_VDZ/Vc_VDZ;
  double C_TGN_r= (A_TGN/V_TGN)*1000;
  double C_PRED = A_PRED/V_PRED;
  double C_UPA  = A_UPA/V_UPA;

  double k10i=CL_IFX/Vc_IFX, k12i=Q_IFX/Vc_IFX, k21i=Q_IFX/Vp_IFX;
  dxdt_A1_IFX = -k10i*A1_IFX - k12i*A1_IFX + k21i*A2_IFX;
  dxdt_A2_IFX =  k12i*A1_IFX - k21i*A2_IFX;
  dxdt_ADA_dep= -ka_ADA*ADA_dep;
  dxdt_A_ADA  =  ka_ADA*F_ADA*ADA_dep - (CL_ADA/V_ADA)*A_ADA;
  double k10u=CL_UST/Vc_UST, k12u=Q_UST/Vc_UST, k21u=Q_UST/Vp_UST;
  dxdt_A1_UST = -k10u*A1_UST - k12u*A1_UST + k21u*A2_UST;
  dxdt_A2_UST =  k12u*A1_UST - k21u*A2_UST;
  dxdt_A_VDZ  = -(CL_VDZ/Vc_VDZ)*A_VDZ;
  dxdt_A_TGN  = ka_TGN*F_AZA*(A_TGN>0?A_TGN:0) - (CL_TGN/V_TGN)*A_TGN;
  dxdt_A_PRED = -(CL_PRED/V_PRED)*A_PRED;
  dxdt_A_UPA  = -(CL_UPA/V_UPA)*A_UPA;

  double E_IFX = Emax_IFX_TNF*C1_IFX/(EC50_IFX_TNF+C1_IFX);
  double E_ADA = Emax_ADA_TNF*C_ADA/(EC50_ADA_TNF+C_ADA);
  double E_antiTNF = 1-(1-E_IFX)*(1-E_ADA);
  if(E_antiTNF>0.97) E_antiTNF=0.97;
  double E_UST = Emax_UST*C1_UST/(EC50_UST+C1_UST);
  double E_VDZ = Emax_VDZ_MI*C_VDZ/(EC50_VDZ_MI+C_VDZ);
  double E_TGN = Emax_TGN_MI*C_TGN_r/(EC50_TGN_MI+C_TGN_r);
  double E_PRED= Emax_PRED*C_PRED/(EC50_PRED+C_PRED);
  double E_UPA = Emax_UPA*C_UPA/(EC50_UPA+C_UPA);
  double E_comb= 1-(1-E_antiTNF)*(1-E_UST*0.55)*(1-E_VDZ*0.40)*(1-E_TGN*0.30)*(1-E_PRED)*(1-E_UPA*0.85);

  double sT=1+0.4*(MI/MI0-1); if(sT<0.1)sT=0.1;
  dxdt_TNF = kprod_TNF*sT*(1-E_antiTNF) - kdeg_TNF*TNF;
  double sI=1+0.3*(MI/MI0-1); if(sI<0.1)sI=0.1;
  dxdt_IL1223 = kprod_IL1223*sI*(1-E_UST)*(1-E_UPA*0.30) - kdeg_IL1223*IL1223;
  double dIL17=1+kIL23_IL17*(IL1223/IL12230-1); if(dIL17<0.1)dIL17=0.1;
  dxdt_IL17 = kprod_IL17*dIL17*(Th17/Th170)*(1-E_UPA*0.55) - kdeg_IL17*IL17;
  double Ts=1+kIL23_Th17*(IL1223/IL12230-1); if(Ts<0.05)Ts=0.05;
  double Tsu=kTreg_Th17*(Treg/Treg0-1);
  dxdt_Th17 = kprol_Th17*Th17*Ts*(1-Tsu)*(1-E_VDZ*0.75)*(1-E_TGN*0.45)*(1-E_UPA*0.65) - kdeg_Th17*Th17;
  double Th1s=1+kIL12_Th1*(IL1223/IL12230-1); if(Th1s<0.05)Th1s=0.05;
  dxdt_Th1 = kprol_Th1*Th1*Th1s*(1-E_VDZ*0.75)*(1-E_TGN*0.40)*(1-E_UPA*0.50) - kdeg_Th1*Th1;
  dxdt_Treg = kprol_Treg*Treg*(1+E_comb*0.15) - kdeg_Treg*Treg;
  double Nd=1+kIL17_Neut*(IL17/IL170-1)+kTNF_Neut*(TNF/TNF0-1); if(Nd<0.1)Nd=0.1;
  dxdt_Neut = kprod_Neut*Neut*Nd*(1-E_comb*0.80) - kdeg_Neut*Neut;
  double MIs=kTNF_MI*(TNF/TNF0-1)+kIL17_MI*(IL17/IL170-1)+kNeut_MI*(Neut/Neut0-1);
  double MIsupp=kTreg_MI*(Treg/Treg0-1);
  double MI_net=1+MIs-MIsupp; if(MI_net<0.05)MI_net=0.05;
  dxdt_MI = kprod_MI*MI_net*(1-E_comb) - kdeg_MI*MI;
  dxdt_CRP= kprod_CRP*MI - kdeg_CRP*CRP; if(CRP<0)CRP=0;
  dxdt_FC = kprod_FC + kNeut_FC*Neut + kMI_FC*MI*(1-E_comb*0.90) - kdeg_FC*FC; if(FC<0)FC=0;
  double BMDloss=(kresorb_BMD+kTNF_BMD*TNF/TNF0+kPRED_BMD*C_PRED)*BMD;
  dxdt_BMD= kform_BMD*BMD - BMDloss;
  double Hp=kprod_Hgb*(1-kMI_Hgb*MI); if(Hp<0)Hp=0;
  dxdt_Hgb= Hp - kdeg_Hgb*Hgb; if(Hgb<0)Hgb=0;

$TABLE
  double C1_IFX_c=A1_IFX/Vc_IFX;
  double C_ADA_c=A_ADA/V_ADA;
  double C1_UST_c=A1_UST/Vc_UST;
  double C_VDZ_c=A_VDZ/Vc_VDZ;
  double C_TGN_c=(A_TGN/V_TGN)*1000;
  double CDAI_c=150+(MI-MI0)*2.2; if(CDAI_c<0)CDAI_c=0;
  double HBI_c=5+(MI-MI0)*0.13; if(HBI_c<0)HBI_c=0;
  double SES_CD_c=(MI/100)*24; if(SES_CD_c<0)SES_CD_c=0;
  double Remission_c=(CDAI_c<150)?1:0;
  double MucHeal_c=(SES_CD_c<=2)?1:0;
  double VDZ_RO=C_VDZ_c/(0.7+C_VDZ_c);

  capture C_IFX=C1_IFX_c, C_ADA=C_ADA_c, C_UST=C1_UST_c, C_VDZ=C_VDZ_c, C_TGN=C_TGN_c;
  capture TNFo=TNF, IL1223o=IL1223, IL17o=IL17, Th17o=Th17, Treg_o=Treg, Neut_o=Neut;
  capture MI_o=MI, CRP_o=CRP, FC_o=FC, BMD_o=BMD, Hgb_o=Hgb;
  capture CDAI_o=CDAI_c, HBI_o=HBI_c, SES_o=SES_CD_c;
  capture Remis=Remission_c, MucHeal=MucHeal_c, VDZ_RO=VDZ_RO;

$CAPTURE C_IFX C_ADA C_UST C_VDZ C_TGN
         TNFo IL1223o IL17o Th17o Treg_o Neut_o
         MI_o CRP_o FC_o BMD_o Hgb_o
         CDAI_o HBI_o SES_o Remis MucHeal VDZ_RO
'

mod <- suppressMessages(mcode("cd_shiny", cd_code))

# Helper: build event table from UI inputs
build_events <- function(drug, dose_ifx, dose_ada, dose_ust, dose_vdz,
                         aza_dose, pred_dose, upa_dose, sim_weeks, wt) {
  ev_list <- list()

  if ("IFX" %in% drug && dose_ifx > 0) {
    amt <- dose_ifx * wt
    times <- unique(c(0, 14, 42, seq(98, sim_weeks*7, by=56)))
    ev_list[["IFX"]] <- data.frame(time=times, amt=amt, cmt="A1_IFX", evid=1, ii=0, addl=0)
  }
  if ("ADA" %in% drug && dose_ada > 0) {
    times <- unique(c(0, 14, 28, seq(42, sim_weeks*7, by=14)))
    amts  <- c(160, 80, dose_ada, rep(dose_ada, length(seq(42, sim_weeks*7, by=14))))
    amts  <- amts[seq_along(times)]
    ev_list[["ADA"]] <- data.frame(time=times, amt=amts, cmt="ADA_dep", evid=1, ii=0, addl=0)
  }
  if ("UST" %in% drug && dose_ust > 0) {
    iv_amt <- round(dose_ust * wt / 10) * 10
    sc_t   <- seq(56, sim_weeks*7, by=56)
    ev_list[["UST"]] <- rbind(
      data.frame(time=0,    amt=iv_amt, cmt="A1_UST", evid=1, ii=0, addl=0),
      data.frame(time=sc_t, amt=90,     cmt="A1_UST", evid=1, ii=0, addl=0)
    )
  }
  if ("VDZ" %in% drug && dose_vdz > 0) {
    times <- unique(c(0, 14, 42, seq(98, sim_weeks*7, by=56)))
    ev_list[["VDZ"]] <- data.frame(time=times, amt=dose_vdz, cmt="A_VDZ", evid=1, ii=0, addl=0)
  }
  if (aza_dose > 0) {
    ev_list[["AZA"]] <- data.frame(time=0:(sim_weeks*7), amt=aza_dose,
                                    cmt="A_TGN", evid=1, ii=0, addl=0)
  }
  if (pred_dose > 0) {
    # 4-wk induction, 4-wk taper
    amts <- c(rep(pred_dose, 28), rep(pred_dose*0.75, 14),
              rep(pred_dose*0.5, 14), rep(pred_dose*0.25, 7))
    ev_list[["PRED"]] <- data.frame(time=0:(length(amts)-1), amt=amts,
                                     cmt="A_PRED", evid=1, ii=0, addl=0)
  }
  if (upa_dose > 0) {
    ind_t   <- 0:55
    maint_t <- 56:(sim_weeks*7)
    ev_list[["UPA"]] <- rbind(
      data.frame(time=ind_t,   amt=upa_dose, cmt="A_UPA", evid=1, ii=0, addl=0),
      data.frame(time=maint_t, amt=30,       cmt="A_UPA", evid=1, ii=0, addl=0)
    )
  }
  if (length(ev_list) == 0) {
    return(ev(time=0, amt=0, cmt="A1_IFX"))
  }
  df <- bind_rows(ev_list) %>% arrange(time)
  as.ev(df)
}

# Predefined comparison scenarios
SCENARIOS <- list(
  "No Treatment"       = function(wt, sw) ev(time=0, amt=0, cmt="A1_IFX"),
  "Infliximab (5mg/kg)"= function(wt, sw) {
    amt <- 5*wt; times <- unique(c(0,14,42,seq(98,sw*7,by=56)))
    as.ev(data.frame(time=times,amt=amt,cmt="A1_IFX",evid=1,ii=0,addl=0))
  },
  "Adalimumab (160→40mg)" = function(wt, sw) {
    times <- unique(c(0,14,28,seq(42,sw*7,by=14)))
    amts  <- c(160,80,40,rep(40,length(seq(42,sw*7,by=14))))
    amts  <- amts[seq_along(times)]
    as.ev(data.frame(time=times,amt=amts,cmt="ADA_dep",evid=1,ii=0,addl=0))
  },
  "Ustekinumab (6mg/kg IV→90mg SC)" = function(wt, sw) {
    iv_amt <- round(6*wt/10)*10; sc_t <- seq(56,sw*7,by=56)
    as.ev(rbind(data.frame(time=0,amt=iv_amt,cmt="A1_UST",evid=1,ii=0,addl=0),
                data.frame(time=sc_t,amt=90,cmt="A1_UST",evid=1,ii=0,addl=0)))
  },
  "Vedolizumab (300mg IV)" = function(wt, sw) {
    times <- unique(c(0,14,42,seq(98,sw*7,by=56)))
    as.ev(data.frame(time=times,amt=300,cmt="A_VDZ",evid=1,ii=0,addl=0))
  },
  "Upadacitinib (45mg QD)" = function(wt, sw) {
    ind_t <- 0:55; mt <- 56:(sw*7)
    as.ev(rbind(data.frame(time=ind_t,amt=45,cmt="A_UPA",evid=1,ii=0,addl=0),
                data.frame(time=mt,  amt=30,cmt="A_UPA",evid=1,ii=0,addl=0)))
  },
  "IFX + AZA (combo)"     = function(wt, sw) {
    amt <- 5*wt; times <- unique(c(0,14,42,seq(98,sw*7,by=56)))
    ifx_ev <- data.frame(time=times,amt=amt,cmt="A1_IFX",evid=1,ii=0,addl=0)
    aza_ev <- data.frame(time=0:(sw*7),amt=150,cmt="A_TGN",evid=1,ii=0,addl=0)
    as.ev(bind_rows(ifx_ev, aza_ev))
  }
)

SCEN_COLORS <- c(
  "No Treatment"="gray50","Infliximab (5mg/kg)"="#E53935",
  "Adalimumab (160→40mg)"="#8E24AA","Ustekinumab (6mg/kg IV→90mg SC)"="#1E88E5",
  "Vedolizumab (300mg IV)"="#00897B","Upadacitinib (45mg QD)"="#F4511E",
  "IFX + AZA (combo)"="#3949AB"
)

# ===========================================================================
# UI
# ===========================================================================
ui <- page_navbar(
  title = "Crohn's Disease QSP Dashboard",
  theme = bs_theme(bootswatch = "flatly", primary = "#1565C0"),

  # ── Tab 1: Patient Profile ─────────────────────────────────────────────
  nav_panel("Patient Profile",
    layout_sidebar(
      sidebar = sidebar(
        title = "Patient Parameters",
        numericInput("wt",    "Body Weight (kg)",  70, 30, 150, 5),
        numericInput("age",   "Age (years)",        35, 10, 80, 1),
        selectInput("sex",    "Sex",   c("Male","Female")),
        selectInput("disease_loc", "Disease Location",
                    c("Ileal (L1)", "Colonic (L2)", "Ileocolonic (L3)",
                      "Upper GI (L4)")),
        selectInput("disease_beh", "Disease Behaviour",
                    c("Non-stricturing/penetrating (B1)",
                      "Stricturing (B2)", "Penetrating (B3)")),
        numericInput("init_cdai", "Baseline CDAI", 280, 150, 600, 10),
        numericInput("init_crp",  "Baseline CRP (mg/L)", 30, 1, 200, 1),
        numericInput("init_fc",   "Baseline FC (μg/g)", 1000, 50, 5000, 50),
        selectInput("prior_bio", "Prior Biologic Use",
                    c("Bio-naïve","One prior biologic (anti-TNF)",
                      "Two+ prior biologics","Anti-TNF failure")),
        numericInput("sim_weeks", "Simulation Duration (weeks)", 78, 12, 156, 4),
        actionButton("btn_profile", "Update Profile", class = "btn-primary w-100 mt-2")
      ),
      fluidRow(
        column(6,
          card(
            card_header("Patient Summary"),
            tableOutput("pt_summary_tbl")
          )
        ),
        column(6,
          card(
            card_header("Baseline Biomarker Severity"),
            plotOutput("baseline_radar", height = "300px")
          )
        )
      ),
      fluidRow(
        column(12,
          card(
            card_header("Disease Activity Classification"),
            verbatimTextOutput("disease_class")
          )
        )
      )
    )
  ),

  # ── Tab 2: Pharmacokinetics ────────────────────────────────────────────
  nav_panel("Pharmacokinetics",
    layout_sidebar(
      sidebar = sidebar(
        title = "Drug & Dosing",
        checkboxGroupInput("drug_pk", "Select Drug(s)",
                           choices = c("IFX","ADA","UST","VDZ"),
                           selected = c("IFX","ADA")),
        conditionalPanel("input.drug_pk.includes('IFX')",
          numericInput("dose_ifx", "IFX Dose (mg/kg)", 5, 1, 10, 0.5)
        ),
        conditionalPanel("input.drug_pk.includes('ADA')",
          numericInput("dose_ada", "ADA Maintenance (mg)", 40, 20, 80, 20)
        ),
        conditionalPanel("input.drug_pk.includes('UST')",
          numericInput("dose_ust", "UST IV induction (mg/kg)", 6, 1, 10, 1)
        ),
        conditionalPanel("input.drug_pk.includes('VDZ')",
          numericInput("dose_vdz", "VDZ dose (mg)", 300, 100, 600, 100)
        ),
        numericInput("aza_pk", "AZA (mg/day, 0=off)", 0, 0, 200, 25),
        numericInput("pred_pk", "Prednisone induction (mg/day, 0=off)", 0, 0, 60, 10),
        numericInput("upa_pk",  "Upadacitinib induction (mg/day, 0=off)", 0, 0, 60, 15),
        actionButton("btn_pk", "Run PK Simulation", class="btn-success w-100 mt-2")
      ),
      fluidRow(
        column(12,
          card(card_header("Drug Concentration–Time Profiles"),
               plotlyOutput("pk_conc_plot", height="450px"))
        )
      ),
      fluidRow(
        column(6,
          card(card_header("Receptor Occupancy (VDZ α4β7)"),
               plotlyOutput("vdz_ro_plot", height="300px"))
        ),
        column(6,
          card(card_header("6-TGN RBC Levels (AZA)"),
               plotlyOutput("tgn_plot", height="300px"))
        )
      ),
      fluidRow(
        column(12,
          card(card_header("PK Summary Table"),
               DTOutput("pk_tbl"))
        )
      )
    )
  ),

  # ── Tab 3: PD / Biomarkers ─────────────────────────────────────────────
  nav_panel("Cytokines & PD",
    layout_sidebar(
      sidebar = sidebar(
        title = "PD Settings",
        checkboxGroupInput("drug_pd", "Active Drug(s)",
                           choices = c("IFX","ADA","UST","VDZ","AZA","PRED","UPA"),
                           selected = "IFX"),
        sliderInput("ifx_dose_pd", "IFX dose (mg/kg)", 1, 10, 5, 0.5),
        numericInput("aza_pd",  "AZA (mg/day)", 150, 0, 250, 25),
        numericInput("pred_pd", "Prednisone (mg/day, induction)", 0, 0, 60, 10),
        numericInput("upa_pd",  "Upadacitinib (mg/day)", 0, 0, 60, 15),
        actionButton("btn_pd", "Run PD Simulation", class="btn-warning w-100 mt-2")
      ),
      fluidRow(
        column(6, card(card_header("TNF-α"), plotlyOutput("tnf_plot", height="270px"))),
        column(6, card(card_header("IL-17A"), plotlyOutput("il17_plot", height="270px")))
      ),
      fluidRow(
        column(6, card(card_header("IL-12/23 Composite"), plotlyOutput("il1223_plot", height="270px"))),
        column(6, card(card_header("T Cell Populations"), plotlyOutput("tcell_plot", height="270px")))
      ),
      fluidRow(
        column(6, card(card_header("Mucosal Inflammation Score"), plotlyOutput("mi_plot", height="270px"))),
        column(6, card(card_header("CRP & Fecal Calprotectin"), plotlyOutput("crp_fc_plot", height="270px")))
      )
    )
  ),

  # ── Tab 4: Clinical Endpoints ──────────────────────────────────────────
  nav_panel("Clinical Endpoints",
    layout_sidebar(
      sidebar = sidebar(
        title = "Endpoint Settings",
        selectInput("drug_clin", "Drug Regimen",
                    choices = c("IFX 5mg/kg","ADA 160/80/40mg",
                                "UST 6mg/kg IV+SC","VDZ 300mg",
                                "UPA 45→30mg","IFX+AZA"),
                    selected = "IFX 5mg/kg"),
        numericInput("wt_clin", "Weight (kg)", 70, 30, 150, 5),
        numericInput("sw_clin", "Simulation Weeks", 78, 12, 156, 4),
        actionButton("btn_clin", "Run", class="btn-info w-100 mt-2")
      ),
      fluidRow(
        column(6, card(card_header("CDAI Over Time"),
                       plotlyOutput("cdai_plot", height="320px"))),
        column(6, card(card_header("SES-CD (Endoscopic Score)"),
                       plotlyOutput("ses_plot",  height="320px")))
      ),
      fluidRow(
        column(4, card(card_header("Hemoglobin"),
                       plotlyOutput("hgb_plot", height="280px"))),
        column(4, card(card_header("Bone Mineral Density"),
                       plotlyOutput("bmd_plot", height="280px"))),
        column(4, card(card_header("Response/Remission Timeline"),
                       plotlyOutput("remis_plot", height="280px")))
      )
    )
  ),

  # ── Tab 5: Scenario Comparison ─────────────────────────────────────────
  nav_panel("Scenario Comparison",
    layout_sidebar(
      sidebar = sidebar(
        title = "Scenario Settings",
        checkboxGroupInput("scen_sel", "Select Scenarios",
                           choices = names(SCENARIOS),
                           selected = c("No Treatment","Infliximab (5mg/kg)",
                                        "Adalimumab (160→40mg)",
                                        "Ustekinumab (6mg/kg IV→90mg SC)")),
        numericInput("wt_scen", "Patient Weight (kg)", 70, 30, 150, 5),
        numericInput("sw_scen", "Simulation Weeks", 78, 12, 156, 4),
        actionButton("btn_scen", "Run Comparison", class="btn-danger w-100 mt-2")
      ),
      fluidRow(
        column(6, card(card_header("CDAI"), plotlyOutput("scen_cdai", height="300px"))),
        column(6, card(card_header("CRP"),  plotlyOutput("scen_crp",  height="300px")))
      ),
      fluidRow(
        column(6, card(card_header("Fecal Calprotectin"), plotlyOutput("scen_fc", height="300px"))),
        column(6, card(card_header("Mucosal Inflammation"), plotlyOutput("scen_mi", height="300px")))
      ),
      fluidRow(
        column(12,
          card(card_header("Week 52 Summary Table"),
               DTOutput("scen_tbl"))
        )
      )
    )
  ),

  # ── Tab 6: Biomarker Correlations ──────────────────────────────────────
  nav_panel("Biomarker Correlations",
    layout_sidebar(
      sidebar = sidebar(
        title = "Correlation Settings",
        selectInput("corr_drug", "Drug Regimen",
                    choices = names(SCENARIOS)[-1],
                    selected = "Infliximab (5mg/kg)"),
        numericInput("wt_corr", "Weight (kg)", 70, 30, 150, 5),
        numericInput("n_pop",   "Population N", 150, 50, 500, 50),
        sliderInput("cv_pk", "PK CV% (IIV)", 10, 60, 35, 5),
        actionButton("btn_corr", "Run Population Sim", class="btn-secondary w-100 mt-2")
      ),
      fluidRow(
        column(6, card(card_header("CRP vs CDAI"), plotlyOutput("corr_crp_cdai", height="300px"))),
        column(6, card(card_header("FC vs SES-CD"), plotlyOutput("corr_fc_ses",   height="300px")))
      ),
      fluidRow(
        column(6, card(card_header("Drug Conc vs CRP (Week 52)"),
                       plotlyOutput("corr_drugconc_crp", height="300px"))),
        column(6, card(card_header("TNF vs Mucosal Inflammation"),
                       plotlyOutput("corr_tnf_mi", height="300px")))
      )
    )
  ),

  # ── Tab 7: Dosing Optimizer ────────────────────────────────────────────
  nav_panel("Dosing Optimizer",
    layout_sidebar(
      sidebar = sidebar(
        title = "Optimization Settings",
        selectInput("opt_drug", "Optimize",
                    c("Infliximab","Adalimumab","Ustekinumab","Vedolizumab")),
        numericInput("wt_opt", "Weight (kg)", 70, 30, 150, 5),
        sliderInput("target_crp", "Target CRP (mg/L)", 1, 30, 5, 1),
        sliderInput("target_fc",  "Target FC (μg/g)", 50, 500, 150, 50),
        numericInput("sw_opt", "Simulation Weeks", 52, 12, 104, 4),
        actionButton("btn_opt", "Run Dose-Response", class="btn-dark w-100 mt-2")
      ),
      fluidRow(
        column(6, card(card_header("Dose–Response: CRP at Week 52"),
                       plotlyOutput("opt_crp_plot", height="350px"))),
        column(6, card(card_header("Dose–Response: CDAI at Week 52"),
                       plotlyOutput("opt_cdai_plot", height="350px")))
      ),
      fluidRow(
        column(12,
          card(card_header("Trough Concentration vs Outcome"),
               plotlyOutput("trough_outcome", height="350px"))
        )
      )
    )
  )
)

# ===========================================================================
# SERVER
# ===========================================================================
server <- function(input, output, session) {

  # ─── Helpers ─────────────────────────────────────────────────────────
  run_sim <- function(ev_obj, wt=70, sw=78) {
    mrgsim(mod, ev_obj, end=sw*7, delta=1, obsonly=TRUE) %>%
      as.data.frame() %>%
      mutate(time_wk = time/7)
  }

  my_plotly <- function(p) ggplotly(p, tooltip=c("x","y","colour")) %>%
    layout(legend=list(orientation="h", y=-0.2))

  # ─── Tab 1: Patient Profile ───────────────────────────────────────────
  output$pt_summary_tbl <- renderTable({
    data.frame(
      Parameter = c("Weight","Age","Sex","Disease Location",
                    "Disease Behaviour","Baseline CDAI",
                    "Baseline CRP","Baseline FC","Prior Biologic"),
      Value     = c(paste(input$wt,"kg"), paste(input$age,"years"),
                    input$sex, input$disease_loc, input$disease_beh,
                    input$init_cdai, paste(input$init_crp,"mg/L"),
                    paste(input$init_fc,"μg/g"), input$prior_bio)
    )
  })

  output$disease_class <- renderText({
    cdai <- input$init_cdai
    cls  <- if (cdai < 150) "Remission (CDAI < 150)"
            else if (cdai < 220) "Mild activity (CDAI 150–219)"
            else if (cdai < 450) "Moderate activity (CDAI 220–449)"
            else "Severe activity (CDAI ≥ 450)"
    paste0("Disease Class: ", cls, "\n",
           "STRIDE-II Target: Clinical remission + endoscopic remission (SES-CD ≤ 2)\n",
           "Biomarker Targets: CRP < 5 mg/L, FC < 150 μg/g\n\n",
           "EBM Treatment Recommendation:\n",
           if (input$prior_bio == "Bio-naïve")
             "  → First-line: IFX, ADA, UST, or VDZ ± immunomodulator\n  → Consider combination (IFX+AZA) for anti-drug antibody prevention"
           else
             "  → Second-line: Switch mechanism (e.g., anti-TNF → UST/VDZ/JAKi)\n  → Consider JAK inhibitor (upadacitinib) for moderate-severe refractory CD")
  })

  output$baseline_radar <- renderPlot({
    d <- data.frame(
      marker = c("CRP","FC","CDAI"),
      value  = c(
        min(input$init_crp / 50, 1),
        min(input$init_fc / 2000, 1),
        min((input$init_cdai - 150) / 450, 1)
      )
    )
    ggplot(d, aes(x=marker, y=value, fill=marker)) +
      geom_col(show.legend=FALSE) +
      scale_fill_manual(values=c("#E53935","#FB8C00","#1E88E5")) +
      ylim(0, 1) +
      labs(y="Severity (normalized)", x="", title="Baseline Severity Profile") +
      theme_minimal(base_size=13)
  })

  # ─── Tab 2: PK ────────────────────────────────────────────────────────
  pk_data <- eventReactive(input$btn_pk, {
    ev_obj <- build_events(
      drug     = input$drug_pk,
      dose_ifx = isolate(input$dose_ifx),
      dose_ada = isolate(input$dose_ada),
      dose_ust = isolate(input$dose_ust),
      dose_vdz = isolate(input$dose_vdz),
      aza_dose = isolate(input$aza_pk),
      pred_dose= isolate(input$pred_pk),
      upa_dose = isolate(input$upa_pk),
      sim_weeks= isolate(input$sim_weeks),
      wt       = isolate(input$wt)
    )
    run_sim(ev_obj, wt=input$wt, sw=input$sim_weeks)
  })

  output$pk_conc_plot <- renderPlotly({
    df <- pk_data()
    long <- df %>%
      select(time_wk, C_IFX, C_ADA, C_UST, C_VDZ) %>%
      pivot_longer(-time_wk, names_to="Drug", values_to="Conc") %>%
      mutate(Drug = recode(Drug, C_IFX="Infliximab", C_ADA="Adalimumab",
                           C_UST="Ustekinumab", C_VDZ="Vedolizumab"))
    p <- ggplot(long, aes(time_wk, Conc, color=Drug)) +
      geom_line(linewidth=0.8) +
      labs(x="Time (weeks)", y="Concentration (mcg/mL)", title="Drug PK") +
      theme_bw(base_size=12)
    my_plotly(p)
  })

  output$vdz_ro_plot <- renderPlotly({
    df <- pk_data()
    p <- ggplot(df, aes(time_wk, VDZ_RO*100)) +
      geom_line(color="#00897B", linewidth=0.9) +
      geom_hline(yintercept=95, linetype="dashed") +
      labs(x="Weeks", y="α4β7 RO (%)", title="VDZ Receptor Occupancy") +
      theme_bw(base_size=12)
    my_plotly(p)
  })

  output$tgn_plot <- renderPlotly({
    df <- pk_data()
    p <- ggplot(df, aes(time_wk, C_TGN)) +
      geom_line(color="#4DB6AC", linewidth=0.9) +
      geom_hline(yintercept=c(235,450), linetype=c("dashed","dotted"),
                 color=c("orange","red")) +
      annotate("text",x=2,y=250,label="Therapeutic lower",size=3,hjust=0) +
      annotate("text",x=2,y=465,label="Therapeutic upper",size=3,hjust=0) +
      labs(x="Weeks", y="6-TGN (pmol/8×10^8 RBC)", title="6-TGN Levels") +
      theme_bw(base_size=12)
    my_plotly(p)
  })

  output$pk_tbl <- renderDT({
    df <- pk_data()
    tbl <- df %>%
      filter(time_wk %in% c(0,2,6,14,22,30,38,52)) %>%
      select(time_wk, C_IFX, C_ADA, C_UST, C_VDZ, VDZ_RO, C_TGN) %>%
      rename(Week=time_wk, `IFX(mcg/mL)`=C_IFX, `ADA(mcg/mL)`=C_ADA,
             `UST(mcg/mL)`=C_UST, `VDZ(mcg/mL)`=C_VDZ,
             `VDZ RO(%)`=VDZ_RO, `6-TGN`=C_TGN) %>%
      mutate(across(where(is.numeric), ~round(., 2)))
    datatable(tbl, options=list(pageLength=10))
  })

  # ─── Tab 3: PD ────────────────────────────────────────────────────────
  pd_data <- eventReactive(input$btn_pd, {
    ev_obj <- build_events(
      drug     = input$drug_pd,
      dose_ifx = isolate(input$ifx_dose_pd),
      dose_ada = 40, dose_ust=6, dose_vdz=300,
      aza_dose = isolate(input$aza_pd),
      pred_dose= isolate(input$pred_pd),
      upa_dose = isolate(input$upa_pd),
      sim_weeks= isolate(input$sim_weeks),
      wt       = isolate(input$wt)
    )
    run_sim(ev_obj, wt=input$wt, sw=input$sim_weeks)
  })

  mk_pd_plot <- function(df, var, ylab, title, col, hline=NULL) {
    p <- ggplot(df, aes(time_wk, .data[[var]])) +
      geom_line(color=col, linewidth=1.0) +
      labs(x="Weeks", y=ylab, title=title) + theme_bw(base_size=11)
    if (!is.null(hline)) p <- p + geom_hline(yintercept=hline, linetype="dashed")
    p
  }

  output$tnf_plot    <- renderPlotly({ my_plotly(mk_pd_plot(pd_data(),"TNFo","TNF-α (pg/mL)","TNF-α","#E53935")) })
  output$il17_plot   <- renderPlotly({ my_plotly(mk_pd_plot(pd_data(),"IL17o","IL-17A (pg/mL)","IL-17A","#1E88E5")) })
  output$il1223_plot <- renderPlotly({ my_plotly(mk_pd_plot(pd_data(),"IL1223o","IL-12/23 (pg/mL)","IL-12/23","#7E57C2")) })
  output$mi_plot     <- renderPlotly({ my_plotly(mk_pd_plot(pd_data(),"MI_o","MI Score (0–100)","Mucosal Inflammation","#F57F17")) })
  output$tcell_plot  <- renderPlotly({
    df <- pd_data() %>%
      select(time_wk, Th17o, Treg_o) %>%
      pivot_longer(-time_wk, names_to="Cell", values_to="val") %>%
      mutate(Cell=recode(Cell, Th17o="Th17", Treg_o="Treg"))
    p <- ggplot(df, aes(time_wk, val, color=Cell)) +
      geom_line(linewidth=0.9) +
      scale_color_manual(values=c(Th17="#C2185B",Treg="#388E3C")) +
      labs(x="Weeks",y="Normalized cell count",title="T Cell Populations") +
      theme_bw(base_size=11)
    my_plotly(p)
  })
  output$crp_fc_plot <- renderPlotly({
    df <- pd_data() %>% select(time_wk, CRP_o, FC_o) %>%
      pivot_longer(-time_wk, names_to="Marker", values_to="val") %>%
      mutate(Marker=recode(Marker, CRP_o="CRP (mg/L)", FC_o="FC (μg/g)"))
    p <- ggplot(df, aes(time_wk, val, color=Marker)) +
      geom_line(linewidth=0.9) +
      facet_wrap(~Marker, scales="free_y") +
      scale_color_manual(values=c("CRP (mg/L)"="#FF7043","FC (μg/g)"="#795548")) +
      labs(x="Weeks",y="Value",title="Inflammatory Biomarkers") +
      theme_bw(base_size=11)
    my_plotly(p)
  })

  # ─── Tab 4: Clinical Endpoints ────────────────────────────────────────
  clin_data <- eventReactive(input$btn_clin, {
    wt <- isolate(input$wt_clin); sw <- isolate(input$sw_clin)
    ev_obj <- switch(input$drug_clin,
      "IFX 5mg/kg"  = {
        times <- unique(c(0,14,42,seq(98,sw*7,by=56)))
        as.ev(data.frame(time=times,amt=5*wt,cmt="A1_IFX",evid=1,ii=0,addl=0))
      },
      "ADA 160/80/40mg" = {
        times <- unique(c(0,14,28,seq(42,sw*7,by=14)))
        amts  <- c(160,80,40,rep(40,length(seq(42,sw*7,by=14))))
        as.ev(data.frame(time=times,amt=amts[seq_along(times)],cmt="ADA_dep",evid=1,ii=0,addl=0))
      },
      "UST 6mg/kg IV+SC" = {
        iv <- round(6*wt/10)*10; sc_t <- seq(56,sw*7,by=56)
        as.ev(rbind(data.frame(time=0,amt=iv,cmt="A1_UST",evid=1,ii=0,addl=0),
                    data.frame(time=sc_t,amt=90,cmt="A1_UST",evid=1,ii=0,addl=0)))
      },
      "VDZ 300mg" = {
        times <- unique(c(0,14,42,seq(98,sw*7,by=56)))
        as.ev(data.frame(time=times,amt=300,cmt="A_VDZ",evid=1,ii=0,addl=0))
      },
      "UPA 45→30mg" = {
        as.ev(rbind(data.frame(time=0:55,amt=45,cmt="A_UPA",evid=1,ii=0,addl=0),
                    data.frame(time=56:(sw*7),amt=30,cmt="A_UPA",evid=1,ii=0,addl=0)))
      },
      "IFX+AZA" = {
        times <- unique(c(0,14,42,seq(98,sw*7,by=56)))
        as.ev(rbind(data.frame(time=times,amt=5*wt,cmt="A1_IFX",evid=1,ii=0,addl=0),
                    data.frame(time=0:(sw*7),amt=150,cmt="A_TGN",evid=1,ii=0,addl=0)))
      }
    )
    run_sim(ev_obj, wt=wt, sw=sw)
  })

  output$cdai_plot <- renderPlotly({
    df <- clin_data()
    p <- ggplot(df, aes(time_wk, CDAI_o)) +
      geom_line(color="#1E88E5", linewidth=1.0) +
      geom_hline(yintercept=150, linetype="dashed", color="green4") +
      annotate("text",x=2,y=155,label="Remission",size=3,hjust=0,color="green4") +
      labs(x="Weeks",y="CDAI",title="Crohn's Disease Activity Index") +
      theme_bw(base_size=12)
    my_plotly(p)
  })
  output$ses_plot  <- renderPlotly({
    df <- clin_data()
    p <- ggplot(df, aes(time_wk, SES_o)) +
      geom_line(color="#8E24AA", linewidth=1.0) +
      geom_hline(yintercept=2.5, linetype="dashed", color="green4") +
      labs(x="Weeks",y="SES-CD",title="Simple Endoscopic Score") +
      theme_bw(base_size=12)
    my_plotly(p)
  })
  output$hgb_plot  <- renderPlotly({
    df <- clin_data()
    p <- ggplot(df, aes(time_wk, Hgb_o)) +
      geom_line(color="#E53935", linewidth=1.0) +
      geom_hline(yintercept=12, linetype="dashed") +
      labs(x="Weeks",y="Hgb (g/dL)",title="Hemoglobin") +
      theme_bw(base_size=12)
    my_plotly(p)
  })
  output$bmd_plot  <- renderPlotly({
    df <- clin_data()
    p <- ggplot(df, aes(time_wk, BMD_o)) +
      geom_line(color="#795548", linewidth=1.0) +
      labs(x="Weeks",y="BMD (norm)",title="Bone Mineral Density") +
      theme_bw(base_size=12)
    my_plotly(p)
  })
  output$remis_plot <- renderPlotly({
    df <- clin_data()
    p <- ggplot(df, aes(time_wk, Remis)) +
      geom_area(fill="#00897B", alpha=0.4) +
      geom_line(color="#00897B", linewidth=1.0) +
      ylim(0, 1.05) +
      labs(x="Weeks",y="In Remission (0/1)",title="Remission Status") +
      theme_bw(base_size=12)
    my_plotly(p)
  })

  # ─── Tab 5: Scenario Comparison ───────────────────────────────────────
  scen_data <- eventReactive(input$btn_scen, {
    wt <- isolate(input$wt_scen); sw <- isolate(input$sw_scen)
    sel <- isolate(input$scen_sel)
    lapply(sel, function(nm) {
      ev_obj <- SCENARIOS[[nm]](wt, sw)
      run_sim(ev_obj, wt=wt, sw=sw) %>% mutate(Scenario=nm)
    }) %>% bind_rows()
  })

  mk_scen_plot <- function(df, var, ylab, title) {
    p <- ggplot(df, aes(time_wk, .data[[var]], color=Scenario)) +
      geom_line(linewidth=0.8, alpha=0.9) +
      scale_color_manual(values=SCEN_COLORS, drop=FALSE) +
      labs(x="Weeks",y=ylab,title=title) +
      theme_bw(base_size=11) +
      theme(legend.position="bottom", legend.title=element_blank())
    ggplotly(p, tooltip=c("x","y","colour")) %>%
      layout(legend=list(orientation="h",y=-0.3))
  }

  output$scen_cdai <- renderPlotly({ mk_scen_plot(scen_data(),"CDAI_o","CDAI","CDAI") })
  output$scen_crp  <- renderPlotly({ mk_scen_plot(scen_data(),"CRP_o","CRP (mg/L)","CRP") })
  output$scen_fc   <- renderPlotly({ mk_scen_plot(scen_data(),"FC_o","FC (μg/g)","Fecal Calprotectin") })
  output$scen_mi   <- renderPlotly({ mk_scen_plot(scen_data(),"MI_o","MI Score","Mucosal Inflammation") })

  output$scen_tbl <- renderDT({
    df <- scen_data() %>%
      filter(abs(time - 364) <= 1) %>%
      group_by(Scenario) %>% slice(1) %>%
      select(Scenario, CDAI_o, MI_o, CRP_o, FC_o, Hgb_o, BMD_o, Remis, MucHeal) %>%
      rename(CDAI=CDAI_o, MI=MI_o, CRP=CRP_o, FC=FC_o,
             Hgb=Hgb_o, BMD=BMD_o, Remission=Remis, `Muc.Heal`=MucHeal) %>%
      mutate(across(where(is.numeric), ~round(.,2)))
    datatable(df, options=list(pageLength=10))
  })

  # ─── Tab 6: Biomarker Correlations ────────────────────────────────────
  corr_data <- eventReactive(input$btn_corr, {
    wt  <- isolate(input$wt_corr)
    sw  <- isolate(input$sw_scen)
    n   <- isolate(input$n_pop)
    cv  <- isolate(input$cv_pk) / 100
    nm  <- isolate(input$corr_drug)

    set.seed(123)
    idata <- tibble(
      ID      = 1:n,
      CL_IFX  = rlnorm(n, log(0.55), cv),
      Vc_IFX  = rlnorm(n, log(4.0),  cv*0.7),
      CL_ADA  = rlnorm(n, log(0.35), cv),
      V_ADA   = rlnorm(n, log(7.8),  cv*0.7),
      CL_UST  = rlnorm(n, log(0.19), cv),
      CL_VDZ  = rlnorm(n, log(0.49), cv)
    )
    ev_obj <- SCENARIOS[[nm]](wt, 78)
    mrgsim(mod, ev_obj, idata=idata, end=364, delta=7, obsonly=TRUE) %>%
      as.data.frame() %>%
      filter(abs(time-364) <= 1) %>%
      group_by(ID) %>% slice(1)
  })

  mk_corr <- function(df, xv, yv, xl, yl) {
    p <- ggplot(df, aes(.data[[xv]], .data[[yv]])) +
      geom_point(alpha=0.5, color="#1E88E5") +
      geom_smooth(method="lm", se=TRUE, color="#E53935") +
      labs(x=xl, y=yl) + theme_bw(base_size=12)
    ggplotly(p)
  }

  output$corr_crp_cdai      <- renderPlotly({ mk_corr(corr_data(),"CRP_o","CDAI_o","CRP (mg/L)","CDAI") })
  output$corr_fc_ses        <- renderPlotly({ mk_corr(corr_data(),"FC_o","SES_o","FC (μg/g)","SES-CD") })
  output$corr_drugconc_crp  <- renderPlotly({ mk_corr(corr_data(),"C_IFX","CRP_o","IFX Conc. (mcg/mL)","CRP (mg/L)") })
  output$corr_tnf_mi        <- renderPlotly({ mk_corr(corr_data(),"TNFo","MI_o","TNF-α (pg/mL)","MI Score") })

  # ─── Tab 7: Dosing Optimizer ──────────────────────────────────────────
  opt_data <- eventReactive(input$btn_opt, {
    wt <- isolate(input$wt_opt); sw <- isolate(input$sw_opt)
    drug <- isolate(input$opt_drug)
    doses <- switch(drug,
      "Infliximab"   = seq(1, 10, by=1),
      "Adalimumab"   = c(20, 40, 80, 160),
      "Ustekinumab"  = seq(1, 10, by=1),
      "Vedolizumab"  = c(100, 200, 300, 450, 600)
    )
    purrr::map_dfr(doses, function(d) {
      ev_obj <- switch(drug,
        "Infliximab" = {
          times <- unique(c(0,14,42,seq(98,sw*7,by=56)))
          as.ev(data.frame(time=times,amt=d*wt,cmt="A1_IFX",evid=1,ii=0,addl=0))
        },
        "Adalimumab" = {
          times <- unique(c(0,14,28,seq(42,sw*7,by=14)))
          amts  <- c(160,80,d,rep(d,length(seq(42,sw*7,by=14))))
          as.ev(data.frame(time=times,amt=amts[seq_along(times)],cmt="ADA_dep",evid=1,ii=0,addl=0))
        },
        "Ustekinumab" = {
          iv <- round(d*wt/10)*10; sc_t <- seq(56,sw*7,by=56)
          as.ev(rbind(data.frame(time=0,amt=iv,cmt="A1_UST",evid=1,ii=0,addl=0),
                      data.frame(time=sc_t,amt=90,cmt="A1_UST",evid=1,ii=0,addl=0)))
        },
        "Vedolizumab" = {
          times <- unique(c(0,14,42,seq(98,sw*7,by=56)))
          as.ev(data.frame(time=times,amt=d,cmt="A_VDZ",evid=1,ii=0,addl=0))
        }
      )
      sim <- mrgsim(mod, ev_obj, end=sw*7, delta=7, obsonly=TRUE) %>%
        as.data.frame() %>% filter(abs(time-sw*7)<=7) %>% slice_tail(n=1)
      data.frame(Dose=d, CRP_w52=sim$CRP_o, CDAI_w52=sim$CDAI_o,
                 FC_w52=sim$FC_o, C_drug=sim$C_IFX+sim$C_ADA+sim$C_UST+sim$C_VDZ)
    })
  })

  output$opt_crp_plot <- renderPlotly({
    df <- opt_data()
    tgt <- isolate(input$target_crp)
    p <- ggplot(df, aes(Dose, CRP_w52)) +
      geom_line(color="#E53935", linewidth=1.2) + geom_point(color="#E53935", size=3) +
      geom_hline(yintercept=tgt, linetype="dashed", color="green4") +
      labs(x="Dose", y="CRP at Week 52 (mg/L)", title="Dose–CRP Response") +
      theme_bw(base_size=12)
    ggplotly(p)
  })
  output$opt_cdai_plot <- renderPlotly({
    df <- opt_data()
    p <- ggplot(df, aes(Dose, CDAI_w52)) +
      geom_line(color="#1E88E5", linewidth=1.2) + geom_point(color="#1E88E5", size=3) +
      geom_hline(yintercept=150, linetype="dashed", color="green4") +
      labs(x="Dose", y="CDAI at Week 52", title="Dose–CDAI Response") +
      theme_bw(base_size=12)
    ggplotly(p)
  })
  output$trough_outcome <- renderPlotly({
    df <- opt_data()
    p <- ggplot(df, aes(C_drug, CRP_w52)) +
      geom_point(aes(size=CDAI_w52), color="#7E57C2", alpha=0.8) +
      geom_smooth(method="loess", se=FALSE, color="#E53935") +
      geom_vline(xintercept=3, linetype="dashed") +
      annotate("text",x=3.3,y=max(df$CRP_w52)*0.9,
               label="Target trough\n≥3 mcg/mL",size=3) +
      labs(x="Drug Trough (mcg/mL)", y="CRP at Week 52 (mg/L)",
           size="CDAI", title="Trough Concentration vs Outcome") +
      theme_bw(base_size=12)
    ggplotly(p)
  })
}

# ===========================================================================
# LAUNCH
# ===========================================================================
shinyApp(ui, server)
