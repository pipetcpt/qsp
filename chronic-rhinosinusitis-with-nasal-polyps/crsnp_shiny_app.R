## =============================================================================
## CRSwNP QSP Shiny App — Interactive Dashboard
## Chronic Rhinosinusitis with Nasal Polyps
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)
library(scales)

## ---- Inline mrgsolve model code ----
crsnp_model_code <- '
$PARAM @annotated
DRUG       : 0     : Drug (0=none,1=Dupilumab,2=Mepolizumab,3=Benralizumab,4=Omalizumab,5=Tezepelumab)
USE_INCS   : 1     : INCS (0/1)
USE_MLK    : 0     : Montelukast (0/1)
F_DUP      : 0.64  : Dupilumab F_sc
KA_DUP     : 0.14  : Dupilumab ka (1/d)
V1_DUP     : 3.10  : Dupilumab V1 (L)
V2_DUP     : 3.40  : Dupilumab V2 (L)
CL_DUP     : 0.50  : Dupilumab CL (L/d)
Q_DUP      : 0.80  : Dupilumab Q (L/d)
F_MEP      : 0.80  : Mepolizumab F_sc
KA_MEP     : 0.12  : Mepolizumab ka (1/d)
V1_MEP     : 3.60  : Mepolizumab V1 (L)
CL_MEP     : 0.22  : Mepolizumab CL (L/d)
F_BEN      : 0.56  : Benralizumab F_sc
KA_BEN     : 0.10  : Benralizumab ka (1/d)
V1_BEN     : 3.10  : Benralizumab V1 (L)
CL_BEN     : 0.35  : Benralizumab CL (L/d)
F_OMA      : 0.62  : Omalizumab F_sc
KA_OMA     : 0.10  : Omalizumab ka (1/d)
V1_OMA     : 4.00  : Omalizumab V1 (L)
CL_OMA     : 0.19  : Omalizumab CL (L/d)
F_TEZ      : 0.77  : Tezepelumab F_sc
KA_TEZ     : 0.09  : Tezepelumab ka (1/d)
V1_TEZ     : 4.00  : Tezepelumab V1 (L)
CL_TEZ     : 0.18  : Tezepelumab CL (L/d)
IC50_DUP   : 1.20  : Dupilumab IC50 (ug/mL)
IC50_MEP   : 0.50  : Mepolizumab IC50 (ug/mL)
IC50_BEN   : 0.30  : Benralizumab IC50 (ug/mL)
IC50_OMA   : 2.50  : Omalizumab IC50 (ug/mL)
IC50_TEZ   : 1.50  : Tezepelumab IC50 (ug/mL)
EMAX_INCS  : 0.50  : INCS Emax
EMAX_MLK   : 0.35  : Montelukast Emax
K_EPI_DMG  : 0.08  : Epithelial damage rate
K_EPI_REP  : 0.05  : Epithelial repair rate
EPI_SS     : 0.40  : Disease SS epithelial integrity
K_TSLP_P   : 0.50  : TSLP production
K_TSLP_D   : 0.50  : TSLP degradation
K_ILC2_P   : 0.30  : ILC2 activation
K_ILC2_D   : 0.25  : ILC2 resolution
K_TH2_P    : 0.20  : Th2 polarization
K_TH2_D    : 0.18  : Th2 resolution
K_IL4_P    : 0.40  : IL-4 production
K_IL4_D    : 0.80  : IL-4 degradation
K_IL5_P    : 0.30  : IL-5 production
K_IL5_D    : 0.60  : IL-5 degradation
K_IL13_P   : 0.30  : IL-13 production
K_IL13_D   : 0.50  : IL-13 degradation
K_IGE_P    : 0.06  : IgE production
K_IGE_D    : 0.030 : IgE degradation
IGE_BAS    : 300   : Baseline IgE (kU/L)
K_EOSB_P   : 0.50  : Blood Eos production
K_EOSB_D   : 0.10  : Blood Eos clearance
EOSB_BAS   : 500   : Baseline blood Eos
K_EOST_R   : 0.04  : Tissue Eos recruitment
K_EOST_D   : 0.20  : Tissue Eos clearance
EOST_BAS   : 80    : Baseline tissue Eos
K_GOBC_P   : 0.20  : Goblet cell expansion
K_GOBC_D   : 0.08  : Goblet cell resolution
K_TGFB_P   : 0.15  : TGF-b production
K_TGFB_D   : 0.25  : TGF-b degradation
K_VEGF_P   : 0.12  : VEGF production
K_VEGF_D   : 0.20  : VEGF degradation
K_NPS_G    : 0.012 : NPS growth rate
K_NPS_R    : 0.015 : NPS resolution rate
NPS_MAX    : 8.0   : Max NPS
NPS_BAS    : 5.5   : Baseline NPS
NPS_SNOT   : 3.5   : NPS-SNOT22 coefficient
EOST_OBS   : 2.5   : EosTissue-obstruction coeff
TGFB_FIB   : 1.8   : TGFb-fibrosis coeff

$CMT @annotated
D_SC : Dupilumab SC depot
D_C1 : Dupilumab central
D_P1 : Dupilumab peripheral
M_SC : Mepolizumab SC
M_C1 : Mepolizumab central
B_SC : Benralizumab SC
B_C1 : Benralizumab central
O_SC : Omalizumab SC
O_C1 : Omalizumab central
T_SC : Tezepelumab SC
T_C1 : Tezepelumab central
EPI TSLP ILC2 TH2 IL4 IL5 IL13 IGE EOSB EOST GOBC TGFB VEGF NPS

$MAIN
EPI_0  = EPI_SS;
TSLP_0 = K_TSLP_P * (1-EPI_SS) / K_TSLP_D;
ILC2_0 = K_ILC2_P * TSLP_0 / K_ILC2_D;
TH2_0  = 0.8;
double src0 = ILC2_0 + TH2_0;
IL4_0  = K_IL4_P * src0 / K_IL4_D;
IL5_0  = K_IL5_P * src0 / K_IL5_D;
IL13_0 = K_IL13_P * src0 / K_IL13_D;
IGE_0  = IGE_BAS; EOSB_0 = EOSB_BAS; EOST_0 = EOST_BAS;
GOBC_0 = K_GOBC_P * IL13_0 / K_GOBC_D;
TGFB_0 = K_TGFB_P * EOST_0 / K_TGFB_D;
VEGF_0 = K_VEGF_P * (EOST_0 + TGFB_0) / K_VEGF_D;
NPS_0  = NPS_BAS;

$ODE
dxdt_D_SC = -KA_DUP*D_SC;
dxdt_D_C1 = KA_DUP*F_DUP*D_SC/V1_DUP - (CL_DUP/V1_DUP)*D_C1 - (Q_DUP/V1_DUP)*D_C1 + (Q_DUP/V2_DUP)*D_P1;
dxdt_D_P1 = (Q_DUP/V1_DUP)*D_C1 - (Q_DUP/V2_DUP)*D_P1;
dxdt_M_SC = -KA_MEP*M_SC;
dxdt_M_C1 = KA_MEP*F_MEP*M_SC/V1_MEP - (CL_MEP/V1_MEP)*M_C1;
dxdt_B_SC = -KA_BEN*B_SC;
dxdt_B_C1 = KA_BEN*F_BEN*B_SC/V1_BEN - (CL_BEN/V1_BEN)*B_C1;
dxdt_O_SC = -KA_OMA*O_SC;
dxdt_O_C1 = KA_OMA*F_OMA*O_SC/V1_OMA - (CL_OMA/V1_OMA)*O_C1;
dxdt_T_SC = -KA_TEZ*T_SC;
dxdt_T_C1 = KA_TEZ*F_TEZ*T_SC/V1_TEZ - (CL_TEZ/V1_TEZ)*T_C1;

double E_DUP  = (DRUG==1) ? D_C1/(IC50_DUP+D_C1) : 0;
double E_MEP  = (DRUG==2) ? M_C1/(IC50_MEP+M_C1) : 0;
double E_BENB = (DRUG==3) ? B_C1/(IC50_BEN+B_C1) : 0;
double E_BENT = (DRUG==3) ? B_C1/(IC50_BEN*0.5+B_C1) : 0;
double E_OMA  = (DRUG==4) ? O_C1/(IC50_OMA+O_C1) : 0;
double E_TEZ  = (DRUG==5) ? T_C1/(IC50_TEZ+T_C1) : 0;
double E_INCS_eff = USE_INCS * EMAX_INCS;

dxdt_EPI  = -K_EPI_DMG*EOST*(1-EPI) + K_EPI_REP*(1-EPI)*(1+E_DUP*0.5+E_TEZ*0.5);
dxdt_TSLP = K_TSLP_P*(1-EPI) - K_TSLP_D*TSLP*(1+E_TEZ);
dxdt_ILC2 = K_ILC2_P*TSLP - K_ILC2_D*ILC2;
double il4_sig = IL4*(1-E_DUP);
dxdt_TH2  = K_TH2_P*il4_sig + 0.05*ILC2 - K_TH2_D*TH2;
double src = ILC2 + TH2;
dxdt_IL4  = K_IL4_P*src*(1-E_INCS_eff) - K_IL4_D*IL4;
dxdt_IL5  = K_IL5_P*src*(1-E_INCS_eff) - K_IL5_D*IL5*(1+E_MEP*3.0);
dxdt_IL13 = K_IL13_P*src*(1-E_INCS_eff) - K_IL13_D*IL13;
dxdt_IGE  = K_IGE_P*(IL4+IL13)*(1-E_DUP*0.8) - K_IGE_D*IGE - E_OMA*K_IGE_D*3.0*IGE;
dxdt_EOSB = K_EOSB_P*IL5*(1-E_MEP*0.9)*(EOSB_BAS/100.0) - K_EOSB_D*EOSB - E_BENB*K_EOSB_D*15.0*EOSB;
dxdt_EOST = K_EOST_R*EOSB*IL5*(1-E_INCS_eff*0.6) - K_EOST_D*EOST - E_BENT*K_EOST_D*8.0*EOST - E_DUP*K_EOST_D*0.8*EOST;
dxdt_GOBC = K_GOBC_P*IL13*(1-E_DUP*0.85) - K_GOBC_D*GOBC;
dxdt_TGFB = K_TGFB_P*EOST - K_TGFB_D*TGFB;
dxdt_VEGF = K_VEGF_P*(EOST+TGFB)*(1-E_INCS_eff*0.5) - K_VEGF_D*VEGF;
double nps_treat = E_DUP*0.040*NPS + E_MEP*0.015*NPS + E_BENT*0.015*NPS
                  + E_OMA*0.018*NPS + E_TEZ*0.022*NPS + E_INCS_eff*0.008*NPS;
dxdt_NPS  = K_NPS_G*EOST*VEGF*(NPS_MAX-NPS)/NPS_MAX - K_NPS_R*NPS - nps_treat;

$TABLE
double Cp_DUP = D_C1; double Cp_MEP = M_C1; double Cp_BEN = B_C1;
double Cp_OMA = O_C1; double Cp_TEZ = T_C1;
double OBS_VAS = fmin(10.0, fmax(0.0, 3.5 + EOST_OBS*EOST/EOST_BAS*2.0 + 0.8*NPS));
double OLFACT = fmax(0.0, 10.0 - 1.2*NPS);
double SNOT22 = fmin(110.0, fmax(0.0, NPS*NPS_SNOT + OBS_VAS*3.0 + (10-OLFACT)*2.0 + GOBC*4.0));
double LM_CT = fmin(24.0, fmax(0.0, NPS*2.0 + TGFB*TGFB_FIB));
double BLD_EOS = EOSB; double SERUM_IGE = IGE;
double FeNO = fmax(5.0, 25.0 - NPS*2.0);
double PERIOSTIN = IL13*8.0 + EOST*1.5;

$CAPTURE Cp_DUP Cp_MEP Cp_BEN Cp_OMA Cp_TEZ
EPI TSLP ILC2 TH2 IL4 IL5 IL13 IGE EOSB EOST GOBC TGFB VEGF NPS
OBS_VAS OLFACT SNOT22 LM_CT BLD_EOS SERUM_IGE FeNO PERIOSTIN
'

## ---- Load model ----
mod_base <- mread("crsnp_shiny", tempdir(), crsnp_model_code, quiet = TRUE)

## ---- Helper functions ----
drug_labels <- c(
  "0" = "No Drug",
  "1" = "Dupilumab (anti-IL-4Rα)",
  "2" = "Mepolizumab (anti-IL-5)",
  "3" = "Benralizumab (anti-IL-5Rα)",
  "4" = "Omalizumab (anti-IgE)",
  "5" = "Tezepelumab (anti-TSLP)"
)
drug_doses  <- c("0"=0,  "1"=300, "2"=100, "3"=30, "4"=300, "5"=210)
drug_intrvl <- c("0"=14, "1"=14,  "2"=28,  "3"=28, "4"=28,  "5"=28)
drug_cmt    <- c("0"=1,  "1"=1,   "2"=4,   "3"=6,  "4"=8,   "5"=10)
drug_colors <- c(
  "No Drug"                  = "#7F8C8D",
  "Dupilumab (anti-IL-4Rα)"  = "#E74C3C",
  "Mepolizumab (anti-IL-5)"  = "#8E44AD",
  "Benralizumab (anti-IL-5Rα)" = "#2980B9",
  "Omalizumab (anti-IgE)"    = "#27AE60",
  "Tezepelumab (anti-TSLP)"  = "#D35400"
)

make_dosing_ev <- function(drug_id, n_weeks = 52) {
  di <- as.character(drug_id)
  if (drug_id == 0) return(NULL)
  dose <- drug_doses[di]; intv <- drug_intrvl[di]; cmt <- drug_cmt[di]
  n_doses <- floor(n_weeks * 7 / intv) + 1
  if (drug_id == 3) {  # Benralizumab: q4wx3 then q8w
    e1 <- ev(time = c(0, 28, 56), amt = dose, cmt = cmt)
    starts <- seq(84, n_weeks * 7, by = 56)
    if (length(starts) > 0)
      e2 <- ev(time = starts, amt = dose, cmt = cmt)
    else e2 <- NULL
    return(c(e1, e2))
  }
  ev(time = seq(0, (n_doses-1)*intv, by = intv), amt = dose, cmt = cmt)
}

run_sim <- function(drug_id, use_incs, use_mlk, n_weeks, ige_bas, eosb_bas, nps_bas) {
  obs_t <- seq(0, n_weeks * 7, by = 7)
  mod_run <- mod_base %>%
    param(DRUG = drug_id, USE_INCS = use_incs, USE_MLK = use_mlk,
          IGE_BAS = ige_bas, EOSB_BAS = eosb_bas, NPS_BAS = nps_bas)
  dose_ev <- make_dosing_ev(drug_id, n_weeks)
  if (is.null(dose_ev)) {
    out <- mrgsim(mod_run, obsonly = TRUE, tgrid = obs_t, digits = 4)
  } else {
    out <- mrgsim(mod_run, events = dose_ev, obsonly = TRUE, tgrid = obs_t, digits = 4)
  }
  as.data.frame(out) %>% mutate(Week = time / 7)
}

## ============================================================================
## UI
## ============================================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "CRSwNP QSP Model"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "profile",    icon = icon("user")),
      menuItem("Drug PK",            tabName = "pk",         icon = icon("pills")),
      menuItem("Cytokine/Biomarkers",tabName = "biomarkers", icon = icon("flask")),
      menuItem("Disease Endpoints",  tabName = "endpoints",  icon = icon("chart-line")),
      menuItem("Scenario Comparison",tabName = "comparison", icon = icon("balance-scale")),
      menuItem("Responder Analysis", tabName = "responder",  icon = icon("dna")),
      menuItem("Long-term Outcomes", tabName = "longterm",   icon = icon("calendar"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML(
      ".box-header { font-weight: bold; }
       .info-box-icon { font-size: 28px; }
       .shiny-plot-output { width: 100%; }"
    ))),
    tabItems(

      ## ---- TAB 1: Patient Profile ----
      tabItem(tabName = "profile",
        fluidRow(
          box(title = "Patient Characteristics", width = 6, status = "primary",
            sliderInput("nps_bas", "Baseline NPS (0-8):", 0, 8, 5.5, 0.5),
            sliderInput("eosb_bas", "Baseline Blood Eosinophils (cells/μL):",
                        100, 2000, 500, 50),
            sliderInput("ige_bas", "Baseline Total IgE (kU/L):", 10, 2000, 300, 10),
            selectInput("endotype", "Disease Endotype:",
                        choices = c("Eosinophilic (Type 2 High)" = "eos",
                                    "Mixed (Intermediate)" = "mixed",
                                    "Non-eosinophilic (IL-17 driven)" = "noneos")),
            checkboxGroupInput("comorbid", "Comorbidities:",
                               choices = c("Asthma" = "asthma",
                                           "AERD (Aspirin-Exacerbated)" = "aerd",
                                           "Atopic Dermatitis" = "ad",
                                           "Allergic Rhinitis" = "ar"))
          ),
          box(title = "Treatment Selection", width = 6, status = "success",
            selectInput("drug_sel", "Biologic Drug:",
                        choices = drug_labels, selected = "1"),
            checkboxInput("use_incs", "Add INCS (Intranasal Corticosteroid)", value = TRUE),
            checkboxInput("use_mlk", "Add Montelukast (Leukotriene Antagonist)", value = FALSE),
            sliderInput("n_weeks", "Simulation Duration (weeks):", 12, 104, 52, 4),
            hr(),
            actionButton("run_sim", "Run Simulation", class = "btn-primary btn-lg",
                         icon = icon("play"))
          )
        ),
        fluidRow(
          valueBoxOutput("vb_nps"),
          valueBoxOutput("vb_snot"),
          valueBoxOutput("vb_eos")
        ),
        fluidRow(
          box(title = "CRSwNP Endotype Guide", width = 12, status = "info",
            p(strong("Eosinophilic (Type 2 High):"), " Blood Eos >300, IgE >100, IL-5 high. Responds best to dupilumab, benralizumab, mepolizumab."),
            p(strong("Non-eosinophilic:"), " Blood Eos <150, IL-17 driven, common in Asian populations. Less responsive to anti-IL-5 biologics."),
            p(strong("AERD (Samter's Triad):"), " NP + asthma + NSAID hypersensitivity. Fibrin-rich polyps. Aspirin desensitization + dupilumab."),
            p(strong("Diagnostic Cutoffs:"), " Blood Eos ≥300 predicts response to IL-5 pathway biologics; IgE ≥30 for omalizumab consideration.")
          )
        )
      ),

      ## ---- TAB 2: Drug PK ----
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Plasma Concentration–Time Profile", width = 9, status = "primary",
            plotlyOutput("pk_plot", height = "450px")
          ),
          box(title = "PK Parameters", width = 3, status = "info",
            tableOutput("pk_table")
          )
        ),
        fluidRow(
          box(title = "Drug Mechanism of Action", width = 12, status = "success",
            fluidRow(
              column(4, h4("Dupilumab (anti-IL-4Rα)"),
                p("Blocks shared IL-4Rα → inhibits IL-4 AND IL-13 signaling."),
                p("Reduces: goblet cell hyperplasia, eosinophil recruitment (CCL26↓), IgE class switching, tissue fibrosis."),
                p("300 mg SC q2w; Ctrough ~60-80 μg/mL at SS")
              ),
              column(4, h4("Mepolizumab/Benralizumab"),
                p(strong("Mepolizumab:"), " Neutralizes free IL-5. Reduces blood eos ~60%. Spares IL-5Rα."),
                p(strong("Benralizumab:"), " Blocks IL-5Rα + ADCC via NK cells → near-complete eosinophil depletion (>90%)."),
                p("q4w×3 then q8w maintenance dosing.")
              ),
              column(4, h4("Omalizumab / Tezepelumab"),
                p(strong("Omalizumab:"), " Neutralizes free IgE → prevents FcεRI loading → mast cell/basophil blockade."),
                p(strong("Tezepelumab:"), " Most upstream — blocks TSLP before ILC2, DC, mast cell, Th2 activation.")
              )
            )
          )
        )
      ),

      ## ---- TAB 3: Cytokines/Biomarkers ----
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(title = "Type 2 Cytokine Dynamics (IL-4, IL-5, IL-13)", width = 12, status = "primary",
            plotlyOutput("cytokine_plot", height = "380px")
          )
        ),
        fluidRow(
          box(title = "Blood Eosinophils Over Time", width = 6, status = "warning",
            plotlyOutput("eos_plot", height = "350px")
          ),
          box(title = "Total Serum IgE Over Time", width = 6, status = "info",
            plotlyOutput("ige_plot", height = "350px")
          )
        ),
        fluidRow(
          box(title = "TGF-β & VEGF (Tissue Remodeling Signals)", width = 6, status = "danger",
            plotlyOutput("remodel_plot", height = "300px")
          ),
          box(title = "TSLP & ILC2 (Upstream Innate Signals)", width = 6, status = "success",
            plotlyOutput("innate_plot", height = "300px")
          )
        )
      ),

      ## ---- TAB 4: Disease Endpoints ----
      tabItem(tabName = "endpoints",
        fluidRow(
          box(title = "Nasal Polyp Score (NPS, 0-8)", width = 6, status = "danger",
            plotlyOutput("nps_plot", height = "350px")
          ),
          box(title = "SNOT-22 Patient-Reported Outcome (0-110)", width = 6, status = "warning",
            plotlyOutput("snot_plot", height = "350px")
          )
        ),
        fluidRow(
          box(title = "Nasal Obstruction VAS (0-10)", width = 6, status = "primary",
            plotlyOutput("obs_plot", height = "300px")
          ),
          box(title = "Olfactory Function Score (0-10)", width = 6, status = "success",
            plotlyOutput("olf_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Endpoint Summary Table", width = 12, status = "info",
            DTOutput("endpoint_table")
          )
        )
      ),

      ## ---- TAB 5: Scenario Comparison ----
      tabItem(tabName = "comparison",
        fluidRow(
          box(title = "Select Drugs to Compare", width = 3, status = "primary",
            checkboxGroupInput("cmp_drugs", "Include in Comparison:",
                               choices = drug_labels[-1], selected = names(drug_labels[-1])),
            checkboxInput("cmp_incs", "All with INCS", value = TRUE),
            sliderInput("cmp_weeks", "Duration (weeks):", 12, 52, 52, 4),
            sliderInput("cmp_nps_bas", "Baseline NPS:", 0, 8, 5.5, 0.5),
            sliderInput("cmp_eosb_bas", "Baseline Eos (cells/μL):", 100, 2000, 500, 50),
            actionButton("run_compare", "Compare All", class = "btn-warning btn-lg",
                         icon = icon("random"))
          ),
          box(title = "NPS Comparison — All Drugs", width = 9, status = "danger",
            plotlyOutput("cmp_nps_plot", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Blood Eosinophil Comparison", width = 6, status = "warning",
            plotlyOutput("cmp_eos_plot", height = "320px")
          ),
          box(title = "SNOT-22 Comparison", width = 6, status = "info",
            plotlyOutput("cmp_snot_plot", height = "320px")
          )
        ),
        fluidRow(
          box(title = "Comparative Efficacy at Week 24 & 52", width = 12, status = "success",
            DTOutput("cmp_table")
          )
        )
      ),

      ## ---- TAB 6: Responder Analysis ----
      tabItem(tabName = "responder",
        fluidRow(
          box(title = "Biomarker-Stratified Response Analysis", width = 12, status = "primary",
            p("Explore how baseline biomarkers (blood eosinophils, IgE, NPS) predict treatment response."),
            fluidRow(
              column(4,
                sliderInput("resp_drug", "Drug:", 1, 5, 1, 1,
                            animate = FALSE),
                selectInput("resp_drug_sel", "Drug:", choices = drug_labels[-1], selected = "1"),
                sliderInput("resp_nps", "Baseline NPS:", 1, 8, 5.5, 0.5)
              ),
              column(8,
                plotlyOutput("resp_eos_sweep", height = "350px")
              )
            )
          )
        ),
        fluidRow(
          box(title = "Blood Eosinophil Sweep: NPS Response at Week 24", width = 6, status = "warning",
            plotlyOutput("sweep_eos_nps", height = "350px")
          ),
          box(title = "Serum IgE Sweep: NPS Response at Week 24 (Omalizumab)", width = 6, status = "info",
            plotlyOutput("sweep_ige_nps", height = "350px")
          )
        ),
        fluidRow(
          box(title = "Biomarker Thresholds Summary", width = 12, status = "success",
            tags$ul(
              tags$li(strong("IL-5 pathway biologics (Mepolizumab/Benralizumab):"), " Blood Eos ≥300 cells/μL predicts greater NPS reduction. Eos >500 → ~50% more response."),
              tags$li(strong("Dupilumab:"), " Effective across endotypes; response less dependent on blood eos (targets IL-4/IL-13 signaling)."),
              tags$li(strong("Omalizumab:"), " Total IgE 30-1500 kU/L + body weight determines dosing. Best for allergen-driven CRSwNP."),
              tags$li(strong("Tezepelumab:"), " Upstream TSLP blockade — effective regardless of eosinophilic status. Reduces ILC2, mast cells, Th2 simultaneously."),
              tags$li(strong("Periostin:"), " >25 ng/mL suggests IL-13 driven fibrotic phenotype → better dupilumab response.")
            )
          )
        )
      ),

      ## ---- TAB 7: Long-term Outcomes ----
      tabItem(tabName = "longterm",
        fluidRow(
          box(title = "Treatment Discontinuation & Relapse Analysis", width = 12, status = "danger",
            p("Simulate what happens when biologic is stopped after initial treatment course."),
            fluidRow(
              column(3,
                selectInput("lt_drug", "Drug:", choices = drug_labels[-1], selected = "1"),
                sliderInput("lt_treat_wk", "Treatment Duration (weeks):", 24, 104, 52, 4),
                sliderInput("lt_obs_wk", "Total Follow-up (weeks):", 52, 156, 104, 4),
                checkboxInput("lt_incs", "INCS Maintained Post-Stop", value = TRUE),
                actionButton("run_lt", "Run Long-term Simulation", class = "btn-danger btn-lg",
                             icon = icon("calendar"))
              ),
              column(9,
                plotlyOutput("lt_nps_plot", height = "400px")
              )
            )
          )
        ),
        fluidRow(
          box(title = "SNOT-22 Long-term Trajectory", width = 6, status = "warning",
            plotlyOutput("lt_snot_plot", height = "300px")
          ),
          box(title = "Blood Eosinophil Recovery Post-Stop", width = 6, status = "info",
            plotlyOutput("lt_eos_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Key Long-term Considerations", width = 12, status = "primary",
            fluidRow(
              column(4,
                h4("Relapse Risk"),
                p("CRSwNP is a chronic condition requiring long-term biologic therapy."),
                p("Relapse typically occurs within 3-6 months of stopping biologics."),
                p("INCS alone delays but does not prevent full relapse.")
              ),
              column(4,
                h4("Treatment Sequencing"),
                p("Surgery (FESS) may be combined with biologics for refractory disease."),
                p("Post-FESS biologics reduce polyp recurrence rates."),
                p("Dupilumab: LIBERTY NP real-world data shows sustained benefit at 2 years.")
              ),
              column(4,
                h4("Biomarker Monitoring"),
                p("Blood eos and serum IgE should be monitored every 3-6 months."),
                p("Periostin and ECP levels track disease activity."),
                p("FeNO may serve as a complementary marker for comorbid asthma.")
              )
            )
          )
        )
      )
    )
  )
)

## ============================================================================
## SERVER
## ============================================================================
server <- function(input, output, session) {

  ## Reactive simulation result (single drug)
  sim_data <- reactiveVal(NULL)

  observeEvent(input$run_sim, {
    withProgress(message = "Running simulation...", {
      df <- run_sim(
        drug_id  = as.integer(input$drug_sel),
        use_incs = as.integer(input$use_incs),
        use_mlk  = as.integer(input$use_mlk),
        n_weeks  = input$n_weeks,
        ige_bas  = input$ige_bas,
        eosb_bas = input$eosb_bas,
        nps_bas  = input$nps_bas
      )
      sim_data(df)
    })
  }, ignoreNULL = FALSE)

  # Auto-run on startup
  observe({
    if (is.null(sim_data())) {
      df <- run_sim(1, 1, 0, 52, 300, 500, 5.5)
      sim_data(df)
    }
  })

  ## Value boxes
  output$vb_nps <- renderValueBox({
    df <- sim_data(); req(df)
    wk24 <- df %>% filter(abs(Week - 24) < 0.6) %>% slice(1)
    delta <- round(wk24$NPS - df$NPS[1], 1)
    valueBox(round(wk24$NPS, 1), paste0("NPS at Wk 24 (Δ", delta, ")"),
             color = ifelse(delta < -1.5, "green", "yellow"), icon = icon("nose"))
  })
  output$vb_snot <- renderValueBox({
    df <- sim_data(); req(df)
    wk24 <- df %>% filter(abs(Week - 24) < 0.6) %>% slice(1)
    valueBox(round(wk24$SNOT22, 0), "SNOT-22 at Wk 24",
             color = ifelse(wk24$SNOT22 < 20, "green", "orange"), icon = icon("comment"))
  })
  output$vb_eos <- renderValueBox({
    df <- sim_data(); req(df)
    wk24 <- df %>% filter(abs(Week - 24) < 0.6) %>% slice(1)
    valueBox(round(wk24$BLD_EOS, 0), "Blood Eos at Wk 24 (cells/μL)",
             color = ifelse(wk24$BLD_EOS < 150, "green", "red"), icon = icon("tint"))
  })

  ## PK plot
  output$pk_plot <- renderPlotly({
    df <- sim_data(); req(df)
    drug_id <- as.integer(input$drug_sel)
    pk_col  <- c("1"="Cp_DUP","2"="Cp_MEP","3"="Cp_BEN","4"="Cp_OMA","5"="Cp_TEZ")[as.character(drug_id)]
    if (is.null(pk_col) || drug_id == 0) {
      p <- ggplot() + labs(title="No drug selected", x="Week", y="Cp (μg/mL)") + theme_bw()
    } else {
      p <- ggplot(df, aes_string(x = "Week", y = pk_col)) +
        geom_line(color = "#E74C3C", linewidth = 1.2) +
        labs(title = paste(drug_labels[as.character(drug_id)], "— Plasma Concentration"),
             x = "Week", y = "Cp (μg/mL)") +
        theme_bw(base_size = 12)
    }
    ggplotly(p)
  })

  output$pk_table <- renderTable({
    data.frame(
      Drug = c("Dupilumab","Mepolizumab","Benralizumab","Omalizumab","Tezepelumab"),
      Route = "SC",
      Dose = c("300 mg","100 mg","30 mg","75-600 mg","210 mg"),
      Interval = c("q2w","q4w","q4w→q8w","q2-4w","q4w"),
      F = c(0.64, 0.80, 0.56, 0.62, 0.77),
      t_half = c("21d","20d","15d","26d","26d")
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  ## Cytokine plot
  output$cytokine_plot <- renderPlotly({
    df <- sim_data(); req(df)
    df_long <- df %>%
      select(Week, IL4, IL5, IL13) %>%
      pivot_longer(cols = -Week, names_to = "Cytokine", values_to = "Level")
    p <- ggplot(df_long, aes(x = Week, y = Level, color = Cytokine)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("IL4"="#E74C3C","IL5"="#8E44AD","IL13"="#F39C12")) +
      labs(title = "Type 2 Cytokine Dynamics", x = "Week", y = "Concentration (AU)") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$eos_plot <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x = Week, y = BLD_EOS)) +
      geom_line(color = "#D35400", linewidth = 1.2) +
      geom_hline(yintercept = 300, linetype = "dashed", color = "red") +
      annotate("text", x = 2, y = 310, label = "300 cells/μL threshold", hjust = 0, size = 3) +
      labs(title = "Blood Eosinophil Count", x = "Week", y = "Eos (cells/μL)") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$ige_plot <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x = Week, y = SERUM_IGE)) +
      geom_line(color = "#2980B9", linewidth = 1.2) +
      labs(title = "Total Serum IgE", x = "Week", y = "IgE (kU/L)") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$remodel_plot <- renderPlotly({
    df <- sim_data(); req(df)
    df_long <- df %>% select(Week, TGFB, VEGF) %>%
      pivot_longer(-Week, names_to = "Factor", values_to = "Level")
    p <- ggplot(df_long, aes(x = Week, y = Level, color = Factor)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("TGFB"="#C0392B","VEGF"="#16A085")) +
      labs(title = "Remodeling Signals: TGF-β & VEGF", x = "Week", y = "AU") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$innate_plot <- renderPlotly({
    df <- sim_data(); req(df)
    df_long <- df %>% select(Week, TSLP, ILC2) %>%
      pivot_longer(-Week, names_to = "Signal", values_to = "Level")
    p <- ggplot(df_long, aes(x = Week, y = Level, color = Signal)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("TSLP"="#8E44AD","ILC2"="#27AE60")) +
      labs(title = "Innate Signals: TSLP & ILC2", x = "Week", y = "AU") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  ## Endpoint plots
  output$nps_plot <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x = Week, y = NPS)) +
      geom_line(color = "#C0392B", linewidth = 1.3) +
      scale_y_continuous(limits = c(0, 8), breaks = 0:8) +
      labs(title = "Nasal Polyp Score (NPS)", x = "Week", y = "NPS (0–8)") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$snot_plot <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x = Week, y = SNOT22)) +
      geom_line(color = "#E67E22", linewidth = 1.3) +
      geom_hline(yintercept = 20, linetype = "dashed", color = "steelblue") +
      labs(title = "SNOT-22 Score", x = "Week", y = "SNOT-22 (0–110)") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$obs_plot <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x = Week, y = OBS_VAS)) +
      geom_line(color = "#2471A3", linewidth = 1.2) +
      scale_y_continuous(limits = c(0, 10)) +
      labs(title = "Nasal Obstruction VAS", x = "Week", y = "VAS (0–10)") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$olf_plot <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x = Week, y = OLFACT)) +
      geom_line(color = "#1E8449", linewidth = 1.2) +
      scale_y_continuous(limits = c(0, 10)) +
      labs(title = "Olfactory Function Score", x = "Week", y = "Score (0–10, higher=better)") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$endpoint_table <- renderDT({
    df <- sim_data(); req(df)
    df %>%
      filter(Week %in% c(0, 12, 24, 36, 52)) %>%
      mutate(across(c(NPS, SNOT22, OBS_VAS, OLFACT, BLD_EOS, SERUM_IGE, FeNO),
                    \(x) round(x, 1))) %>%
      select(Week, NPS, SNOT22, OBS_VAS, OLFACT, BLD_EOS, SERUM_IGE, FeNO) %>%
      rename("Blood Eos" = BLD_EOS, "IgE (kU/L)" = SERUM_IGE,
             "Obstruction VAS" = OBS_VAS, "Olfactory" = OLFACT)
  }, options = list(pageLength = 10))

  ## Scenario comparison
  cmp_data <- reactiveVal(NULL)

  observeEvent(input$run_compare, {
    withProgress(message = "Comparing all drugs...", {
      drug_ids <- as.integer(input$cmp_drugs)
      use_incs <- as.integer(input$cmp_incs)
      results  <- lapply(drug_ids, function(did) {
        df <- run_sim(did, use_incs, 0, input$cmp_weeks,
                      300, input$cmp_eosb_bas, input$cmp_nps_bas)
        df$Drug <- drug_labels[as.character(did)]
        df
      })
      # Add no-drug baseline
      df0 <- run_sim(0, use_incs, 0, input$cmp_weeks, 300, input$cmp_eosb_bas, input$cmp_nps_bas)
      df0$Drug <- "INCS Only"
      cmp_data(bind_rows(results, df0))
    })
  }, ignoreNULL = FALSE)

  observe({
    if (is.null(cmp_data())) {
      results <- lapply(1:5, function(did) {
        df <- run_sim(did, 1, 0, 52, 300, 500, 5.5)
        df$Drug <- drug_labels[as.character(did)]
        df
      })
      df0 <- run_sim(0, 1, 0, 52, 300, 500, 5.5); df0$Drug <- "INCS Only"
      cmp_data(bind_rows(results, df0))
    }
  })

  output$cmp_nps_plot <- renderPlotly({
    df <- cmp_data(); req(df)
    p <- ggplot(df, aes(x = Week, y = NPS, color = Drug)) +
      geom_line(linewidth = 1.0) +
      scale_y_continuous(limits = c(0, 8), breaks = 0:8) +
      labs(title = "NPS Comparison — All Biologics vs. INCS Only",
           x = "Week", y = "NPS (0–8)") +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$cmp_eos_plot <- renderPlotly({
    df <- cmp_data(); req(df)
    p <- ggplot(df, aes(x = Week, y = BLD_EOS, color = Drug)) +
      geom_line(linewidth = 1.0) +
      labs(title = "Blood Eosinophils", x = "Week", y = "Eos (cells/μL)") +
      theme_bw(base_size = 12) +
      theme(legend.position = "none")
    ggplotly(p)
  })

  output$cmp_snot_plot <- renderPlotly({
    df <- cmp_data(); req(df)
    p <- ggplot(df, aes(x = Week, y = SNOT22, color = Drug)) +
      geom_line(linewidth = 1.0) +
      labs(title = "SNOT-22", x = "Week", y = "SNOT-22 (0–110)") +
      theme_bw(base_size = 12) +
      theme(legend.position = "none")
    ggplotly(p)
  })

  output$cmp_table <- renderDT({
    df <- cmp_data(); req(df)
    df %>%
      filter(Week %in% c(0, 24, 52)) %>%
      group_by(Drug, Week) %>%
      summarise(NPS = round(mean(NPS), 2), SNOT22 = round(mean(SNOT22), 1),
                BLD_EOS = round(mean(BLD_EOS), 0), .groups = "drop") %>%
      arrange(Week, NPS)
  }, options = list(pageLength = 18))

  ## Biomarker sweep (responder analysis)
  output$sweep_eos_nps <- renderPlotly({
    eos_vals <- seq(100, 2000, by = 100)
    drug_id <- as.integer(input$resp_drug_sel)
    res <- lapply(eos_vals, function(eos) {
      df <- run_sim(drug_id, 1, 0, 24, 300, eos, input$resp_nps)
      wk24 <- df %>% filter(abs(Week - 24) < 0.6) %>% slice(1)
      data.frame(Eos_Baseline = eos,
                 NPS_Wk24 = wk24$NPS,
                 Delta_NPS = wk24$NPS - input$resp_nps)
    })
    res_df <- bind_rows(res)
    p <- ggplot(res_df, aes(x = Eos_Baseline, y = Delta_NPS)) +
      geom_line(color = "#8E44AD", linewidth = 1.2) +
      geom_hline(yintercept = -0.9, linetype = "dashed", color = "red") +
      labs(title = paste("NPS Change at Wk 24 vs. Baseline Eos\n",
                         drug_labels[as.character(drug_id)]),
           x = "Baseline Blood Eos (cells/μL)",
           y = "ΔNPS at Week 24 (negative=improvement)") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  output$sweep_ige_nps <- renderPlotly({
    ige_vals <- seq(30, 1500, by = 50)
    res <- lapply(ige_vals, function(ige) {
      df <- run_sim(4, 1, 0, 24, ige, 500, input$resp_nps)  # omalizumab
      wk24 <- df %>% filter(abs(Week - 24) < 0.6) %>% slice(1)
      data.frame(IgE_Baseline = ige, Delta_NPS = wk24$NPS - input$resp_nps)
    })
    res_df <- bind_rows(res)
    p <- ggplot(res_df, aes(x = IgE_Baseline, y = Delta_NPS)) +
      geom_line(color = "#27AE60", linewidth = 1.2) +
      labs(title = "NPS Change at Wk 24 vs. Baseline IgE\nOmalizumab",
           x = "Baseline Serum IgE (kU/L)",
           y = "ΔNPS at Week 24") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })

  ## Long-term simulation
  lt_data <- reactiveVal(NULL)

  observeEvent(input$run_lt, {
    withProgress(message = "Running long-term simulation...", {
      drug_id  <- as.integer(input$lt_drug)
      obs_t1   <- seq(0, input$lt_treat_wk * 7, by = 7)
      mod_on   <- mod_base %>% param(DRUG = drug_id, USE_INCS = 1)
      dose_ev  <- make_dosing_ev(drug_id, input$lt_treat_wk)
      p1 <- if (is.null(dose_ev)) {
        mrgsim(mod_on, obsonly = TRUE, tgrid = obs_t1, digits = 4)
      } else {
        mrgsim(mod_on, events = dose_ev, obsonly = TRUE, tgrid = obs_t1, digits = 4)
      }
      p1df <- as.data.frame(p1) %>% mutate(Phase = "Treatment")

      # Post-stop
      last_row <- tail(p1df, 1)
      init_l <- as.list(last_row[, c("EPI","TSLP","ILC2","TH2","IL4","IL5",
                                      "IL13","IGE","EOSB","EOST","GOBC","TGFB","VEGF","NPS")])
      mod_off <- mod_base %>% param(DRUG = 0, USE_INCS = as.integer(input$lt_incs)) %>%
        init(init_l)
      obs_t2  <- seq(0, (input$lt_obs_wk - input$lt_treat_wk) * 7, by = 7)
      p2df <- as.data.frame(mrgsim(mod_off, obsonly = TRUE, tgrid = obs_t2, digits = 4)) %>%
        mutate(time = time + input$lt_treat_wk * 7, Phase = "Post-Stop")

      lt_data(bind_rows(p1df, p2df) %>% mutate(Week = time / 7))
    })
  })

  output$lt_nps_plot <- renderPlotly({
    df <- lt_data()
    if (is.null(df)) {
      return(plotly::plot_ly() %>% layout(title = "Press 'Run Long-term Simulation'"))
    }
    stop_wk <- df %>% filter(Phase == "Treatment") %>% pull(Week) %>% max()
    p <- ggplot(df, aes(x = Week, y = NPS, color = Phase)) +
      geom_line(linewidth = 1.2) +
      geom_vline(xintercept = stop_wk, linetype = "dotted", color = "red", linewidth = 1) +
      annotate("text", x = stop_wk + 1, y = 7, label = "Drug stopped", hjust = 0,
               color = "red", size = 3.5) +
      scale_y_continuous(limits = c(0, 8), breaks = 0:8) +
      scale_color_manual(values = c("Treatment" = "#2471A3", "Post-Stop" = "#E74C3C")) +
      labs(title = "NPS: Treatment & Post-Discontinuation Relapse",
           x = "Week", y = "NPS (0–8)") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$lt_snot_plot <- renderPlotly({
    df <- lt_data(); req(df)
    stop_wk <- df %>% filter(Phase == "Treatment") %>% pull(Week) %>% max()
    p <- ggplot(df, aes(x = Week, y = SNOT22, color = Phase)) +
      geom_line(linewidth = 1.1) +
      geom_vline(xintercept = stop_wk, linetype = "dotted", color = "red") +
      scale_color_manual(values = c("Treatment"="#E67E22","Post-Stop"="#E74C3C")) +
      labs(title = "SNOT-22 Long-term", x = "Week", y = "SNOT-22 (0–110)") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$lt_eos_plot <- renderPlotly({
    df <- lt_data(); req(df)
    stop_wk <- df %>% filter(Phase == "Treatment") %>% pull(Week) %>% max()
    p <- ggplot(df, aes(x = Week, y = BLD_EOS, color = Phase)) +
      geom_line(linewidth = 1.1) +
      geom_vline(xintercept = stop_wk, linetype = "dotted", color = "red") +
      scale_color_manual(values = c("Treatment"="#2980B9","Post-Stop"="#E74C3C")) +
      labs(title = "Blood Eosinophil Recovery", x = "Week", y = "Eos (cells/μL)") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })
}

## ============================================================================
## Run App
## ============================================================================
shinyApp(ui, server)
