################################################################################
# SLE QSP Shiny App
# Interactive dashboard for the Systemic Lupus Erythematosus QSP model
#
# Tabs:
#   1. Patient Profile & Disease Severity
#   2. Drug PK Simulator (HCQ / Belimumab / Anifrolumab / MMF / Voclosporin)
#   3. IFN-BAFF-B Cell Axis
#   4. Anti-dsDNA & Complement
#   5. Lupus Nephritis Monitor
#   6. SLEDAI Activity & Response
#   7. Population Variability
#   8. Mechanistic Map Viewer
################################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(DT)
library(plotly)

# ==============================================================================
# MODEL (inline — same as sle_model.R but stripped to essentials)
# ==============================================================================

sle_model_code <- '
$PROB SLE QSP Shiny v1.0

$PARAM @annotated
Ka_HCQ    : 0.020  : /day  HCQ absorption rate
F_HCQ     : 0.40   : -    HCQ bioavailability
Vd_HCQ    : 5200   : L    HCQ volume of distribution
CL_HCQ    : 100    : L/day HCQ clearance (t1/2 ~50 days)
CL_BELI   : 0.215  : L/day Belimumab CL
Q_BELI    : 0.45   : L/day Belimumab Q
V1_BELI   : 5.29   : L    Belimumab V1
V2_BELI   : 3.46   : L    Belimumab V2
Rin_BAFF  : 0.42   : ng/mL/day BAFF synthesis
kout_BAFF : 0.15   : /day BAFF elimination
kon_BELI  : 0.60   : 1/(ng/mL*day) BELI-BAFF association
koff_BELI : 0.002  : /day BELI-BAFF dissociation
kint_BELI : 0.20   : /day BELI-BAFF complex internalization
CL_ANIF   : 0.170  : L/day Anifrolumab CL
V_ANIF    : 6.8    : L    Anifrolumab V
Ka_MMF    : 2.5    : /day MMF absorption
F_MMF     : 0.93   : -    MMF/MPA bioavailability
CL_MPA    : 252    : L/day MPA clearance (Vd*CL_h)
Ka_VOC    : 1.2    : /day VOC absorption
F_VOC     : 0.29   : -    VOC bioavailability
Vd_VOC    : 2300   : L    VOC volume
CL_VOC    : 1512   : L/day VOC clearance
IC50_VOC  : 0.8    : ng/mL VOC IC50 calcineurin
Imax_VOC  : 0.95   : -    VOC Imax
CS_dose   : 0      : mg/day Prednisone dose
EC50_CS   : 5.0    : mg/day CS EC50
Emax_CS   : 0.80   : -    CS Emax
EC50_TLR_HCQ : 500 : ng/mL HCQ EC50 TLR inhibition
Imax_HCQ  : 0.85   : -    HCQ Imax
k_IFNA    : 0.20   : /day IFN-a clearance
k_pDC_IFN : 0.50  : /day pDC->IFN coupling
k_IFN_BAFF : 0.015 : /day IFN->BAFF upregulation
IFNa_SLE  : 8.0   : IU/mL Disease IFN-alpha
kin_TLR   : 0.10  : /day TLR baseline
kout_TLR  : 0.10  : /day TLR decay
k_IFN_score : 1.2 : -   IFN score coupling
kin_Bcell  : 3.0  : cells/uL/day B cell production
kout_Bcell : 0.030 : /day B cell turnover
k_BAFF_Bcell : 0.10 : /day BAFF->B survival
EC50_BAFF_B : 2.5 : ng/mL BAFF EC50
k_Blast   : 0.05  : /day B->Plasmablast
kout_Blast : 0.25 : /day Plasmablast decay
k_LLPC    : 0.02  : /day PB->LLPC
kout_LLPC : 0.008 : /day LLPC turnover
k_Ab      : 0.15  : /day LLPC->AntiDsDNA
kout_Ab   : 0.035 : /day AntiDsDNA CL
Anti_SLE  : 800   : IU/mL Disease anti-dsDNA
k_IC      : 0.002 : /day/(IU/mL) Ab->IC
kout_IC   : 0.15  : /day IC clearance
C3_normal : 1.0   : g/L Normal C3
kin_C3    : 0.12  : g/L/day C3 synthesis
kout_C3_base : 0.12 : /day C3 baseline consumption
k_IC_C3   : 0.15  : /day/(AU) IC->C3 consumption
kin_C4    : 0.030 : g/L/day C4 synthesis
kout_C4_base : 0.12 : /day C4 baseline consumption
k_IC_C4   : 0.10  : /day/(AU) IC->C4 consumption
UPCR_SLE  : 2.0   : mg/mg Disease UPCR
k_IC_prot : 0.80  : /day/(AU) IC->proteinuria
kout_prot : 0.35  : /day Proteinuria resolution
eGFR_ss   : 75    : mL/min eGFR baseline
k_prot_GFR : 0.008 : /day GFR worsening
k_GFR_recov : 0.003 : /day GFR recovery
GFR_min   : 15    : mL/min GFR floor

$CMT @annotated
HCQ_gut   : HCQ gut (mg)
HCQ_cent  : HCQ central (mg)
BELI_cent : Belimumab central
BELI_periph : Belimumab peripheral
BAFF_free : Free BAFF (ng/mL)
BELI_cmplx : Beli-BAFF complex
ANIF_cent : Anifrolumab central
MMF_gut   : MMF gut (mg)
MPA_cent  : MPA central (mg)
VOC_gut   : VOC gut (mg)
VOC_cent  : VOC central (mg)
TLR_act   : TLR signal (RU)
IFNa_conc : IFN-alpha (IU/mL)
IFNscore  : IFN gene score
Bcell     : B cell count (cells/uL)
PBlast    : Plasmablast (RU)
LLPC      : Long-lived plasma cell (RU)
AntiDsDNA : Anti-dsDNA (IU/mL)
IC_burden : Immune complexes (AU)
C3_serum  : C3 (g/L)
C4_serum  : C4 (g/L)
Proteinuria : UPCR (mg/mg)
eGFR_ode  : eGFR (mL/min)

$MAIN
if (NEWIND <= 1) {
  HCQ_gut_0    = 0;    HCQ_cent_0  = 0;
  BELI_cent_0  = 0;    BELI_periph_0 = 0;
  BAFF_free_0  = Rin_BAFF / kout_BAFF;
  BELI_cmplx_0 = 0;    ANIF_cent_0 = 0;
  MMF_gut_0    = 0;    MPA_cent_0  = 0;
  VOC_gut_0    = 0;    VOC_cent_0  = 0;
  TLR_act_0    = 1.5;  IFNa_conc_0 = IFNa_SLE;
  IFNscore_0   = 4.0;  Bcell_0     = 80.0;
  PBlast_0     = 15.0; LLPC_0      = 40.0;
  AntiDsDNA_0  = Anti_SLE;
  IC_burden_0  = 3.5;  C3_serum_0  = 0.6;
  C4_serum_0   = 0.10;
  Proteinuria_0 = UPCR_SLE;
  eGFR_ode_0   = eGFR_ss;
}

double C_HCQ  = HCQ_cent  / Vd_HCQ * 1000.0;
double C_BELI = BELI_cent / V1_BELI;
double C_ANIF = ANIF_cent / V_ANIF;
double C_MPA  = MPA_cent  / (3.6 * 70.0);
double C_VOC  = VOC_cent  / Vd_VOC * 1e6;

double E_HCQ    = Imax_HCQ * C_HCQ / (EC50_TLR_HCQ + C_HCQ);
double ANIF_ng  = C_ANIF * 1000.0;
double E_ANIF   = ANIF_ng / (150.0 + ANIF_ng);
double E_MPA    = C_MPA   / (0.25 + C_MPA);
double E_VOC    = Imax_VOC * C_VOC / (IC50_VOC + C_VOC);
double E_CS     = Emax_CS  * CS_dose / (EC50_CS + CS_dose);

$ODE
dxdt_HCQ_gut   = -Ka_HCQ * HCQ_gut;
dxdt_HCQ_cent  = Ka_HCQ * HCQ_gut - (CL_HCQ / Vd_HCQ) * HCQ_cent;

double k10b = CL_BELI / V1_BELI;
double k12b = Q_BELI  / V1_BELI;
double k21b = Q_BELI  / V2_BELI;
dxdt_BELI_cent   = -(k10b+k12b)*BELI_cent + k21b*BELI_periph
                    - kon_BELI*(BELI_cent/V1_BELI)*BAFF_free*V1_BELI
                    + koff_BELI*BELI_cmplx*V1_BELI;
dxdt_BELI_periph = k12b*BELI_cent - k21b*BELI_periph;
dxdt_BAFF_free   = Rin_BAFF - kout_BAFF*BAFF_free
                    - kon_BELI*(BELI_cent/V1_BELI)*BAFF_free
                    + koff_BELI*BELI_cmplx
                    + k_IFN_BAFF*IFNa_conc;
dxdt_BELI_cmplx  = kon_BELI*(BELI_cent/V1_BELI)*BAFF_free
                    - (koff_BELI+kint_BELI)*BELI_cmplx;

dxdt_ANIF_cent = -(CL_ANIF/V_ANIF)*ANIF_cent;
dxdt_MMF_gut   = -Ka_MMF * MMF_gut;
dxdt_MPA_cent  = Ka_MMF*MMF_gut*F_MMF - (CL_MPA/(3.6*70))*MPA_cent;
dxdt_VOC_gut   = -Ka_VOC * VOC_gut;
dxdt_VOC_cent  = Ka_VOC*VOC_gut*F_VOC - (CL_VOC/Vd_VOC)*VOC_cent;

double TLR_inh = E_HCQ + 0.3*E_CS;
TLR_inh = (TLR_inh > 0.95) ? 0.95 : TLR_inh;
dxdt_TLR_act  = kin_TLR*1.5*(1.0-TLR_inh) - kout_TLR*TLR_act;
dxdt_IFNa_conc = k_pDC_IFN*TLR_act - k_IFNA*IFNa_conc;
double IFN_eff = IFNa_conc*(1.0-E_ANIF);
dxdt_IFNscore = k_IFN_score*IFN_eff - 0.20*IFNscore;

double BAFF_eff = BAFF_free/(EC50_BAFF_B+BAFF_free);
double Bp = kin_Bcell*(1.0-0.5*E_MPA)*(1.0-0.3*E_CS);
double Bd = kout_Bcell*Bcell*(1.0-k_BAFF_Bcell*BAFF_eff);
dxdt_Bcell   = Bp - Bd;
double kb_e  = k_Blast*(1.0+0.5*IFNa_conc/IFNa_SLE)*(1.0-0.4*E_MPA)*(1.0-0.3*E_CS);
dxdt_PBlast  = kb_e*Bcell - kout_Blast*PBlast;
double kL_e  = k_LLPC*(1.0-0.5*E_VOC)*(1.0-0.2*E_CS);
dxdt_LLPC    = kL_e*PBlast - kout_LLPC*LLPC;
dxdt_AntiDsDNA = k_Ab*LLPC - kout_Ab*AntiDsDNA;
dxdt_IC_burden = k_IC*AntiDsDNA - kout_IC*IC_burden;
dxdt_C3_serum  = kin_C3 - (kout_C3_base + k_IC_C3*IC_burden)*C3_serum;
dxdt_C4_serum  = kin_C4 - (kout_C4_base + k_IC_C4*IC_burden)*C4_serum;
double prot_res = kout_prot*Proteinuria*(1.0+1.5*E_VOC+0.5*E_CS);
dxdt_Proteinuria = k_IC_prot*IC_burden - prot_res;
if (Proteinuria < 0.04) dxdt_Proteinuria = 0;
dxdt_eGFR_ode  = k_GFR_recov*(90.0-eGFR_ode) - k_prot_GFR*Proteinuria*eGFR_ode;
if (eGFR_ode <= GFR_min) dxdt_eGFR_ode = 0;

$TABLE
double C_HCQ_ng  = HCQ_cent / Vd_HCQ * 1000;
double C_BELI_ug = BELI_cent / V1_BELI;
double C_ANIF_ug = ANIF_cent / V_ANIF;
double C_MPA_ug  = MPA_cent  / (3.6 * 70);
double C_VOC_ng  = VOC_cent  / Vd_VOC * 1e6;
double dsDNA_pts = (AntiDsDNA>200) ? 4.0 : ((AntiDsDNA>100)?2.0:0.0);
double C3C4_pts  = ((C3_serum<0.8)?2.0:0.0)+((C4_serum<0.16)?2.0:0.0);
double rn_pts    = (Proteinuria>0.5)?4.0:((Proteinuria>0.2)?2.0:0.0);
double SLEDAI_c  = dsDNA_pts + C3C4_pts + rn_pts + 4.0;
SLEDAI_c = (SLEDAI_c>30) ? 30.0 : SLEDAI_c;
double CR_r = (Proteinuria < 0.5 && eGFR_ode > 60) ? 1.0 : 0.0;
double IFN_norm = IFNscore * 2.5;
double Bc_pct   = Bcell / 250.0 * 100;

capture C_HCQ_ng C_BELI_ug C_ANIF_ug C_MPA_ug C_VOC_ng
capture SLEDAI_c IFN_norm Bc_pct BAFF_free
capture AntiDsDNA C3_serum C4_serum Proteinuria eGFR_ode IC_burden CR_r

$CAPTURE
C_HCQ_ng C_BELI_ug C_ANIF_ug C_MPA_ug C_VOC_ng
SLEDAI_c IFN_norm Bc_pct BAFF_free
AntiDsDNA C3_serum C4_serum Proteinuria eGFR_ode IC_burden CR_r
'

sle_mod <- mcode("SLE_QSP_shiny", sle_model_code)

# ==============================================================================
# UI
# ==============================================================================

ui <- dashboardPage(
  skin = "purple",

  dashboardHeader(
    title = span(icon("dna"), "SLE QSP Dashboard"),
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      id = "sidebar",
      menuItem("Patient Profile",     tabName = "patient",   icon = icon("user")),
      menuItem("Drug PK Simulator",   tabName = "pk",        icon = icon("pills")),
      menuItem("IFN & BAFF Axis",     tabName = "ifn",       icon = icon("bolt")),
      menuItem("Anti-dsDNA & Complement", tabName = "abs",   icon = icon("shield-halved")),
      menuItem("Lupus Nephritis",     tabName = "nephritis", icon = icon("kidneys")),
      menuItem("SLEDAI & Response",   tabName = "sledai",    icon = icon("chart-line")),
      menuItem("Population Variability", tabName = "popvar", icon = icon("users")),
      menuItem("Mechanistic Map",     tabName = "map",       icon = icon("diagram-project")),
      menuItem("References",          tabName = "refs",      icon = icon("book"))
    ),
    hr(),
    # Global drug toggles
    h5("  Active Treatments", style = "color:#DDD; padding-left:15px;"),
    checkboxInput("use_hcq",  "Hydroxychloroquine (HCQ)", value = TRUE),
    checkboxInput("use_beli", "Belimumab",                value = FALSE),
    checkboxInput("use_anif", "Anifrolumab",              value = FALSE),
    checkboxInput("use_mmf",  "MMF",                      value = FALSE),
    checkboxInput("use_voc",  "Voclosporin",               value = FALSE),
    sliderInput("cs_dose", "Prednisone (mg/day)",
                min = 0, max = 60, value = 0, step = 2.5),
    sliderInput("sim_weeks", "Simulation duration (weeks)",
                min = 4, max = 104, value = 52, step = 4),
    actionButton("run_sim", "Run Simulation",
                 icon = icon("play"), width = "90%",
                 style = "background:#7B2D8B; color:white; margin:8px 5%")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background: #F5F5F5; }
      .box { border-radius:6px; }
      .value-box .inner h3 { font-size:22px; }
      .nav-tabs-custom .nav-tabs li.active a { font-weight:bold; }
    "))),

    tabItems(

      # ---- TAB 1: Patient Profile -------------------------------------------
      tabItem(tabName = "patient",
        fluidRow(
          box(title = "Patient Disease Profile", width = 6, status = "purple",
            solidHeader = TRUE,
            selectInput("disease_activity",
                        "Baseline Disease Activity",
                        choices = c("Mild (SLEDAI 4-6)"      = "mild",
                                    "Moderate (SLEDAI 6-12)" = "moderate",
                                    "Severe (SLEDAI >12)"    = "severe"),
                        selected = "moderate"),
            checkboxGroupInput("organ_involvement",
                               "Organ Involvement",
                               choices = c("Lupus Nephritis (Class III/IV)" = "LN",
                                           "Neuropsychiatric (NPSLE)"       = "CNS",
                                           "Haematologic (cytopenias)"      = "Heme",
                                           "Mucocutaneous (malar, discoid)" = "Skin",
                                           "Musculoskeletal (arthritis)"    = "MSK"),
                               selected = c("LN", "Skin")),
            sliderInput("baseline_upcr",  "Baseline UPCR (mg/mg)", 0.1, 5.0, 2.0, 0.1),
            sliderInput("baseline_egfr",  "Baseline eGFR (mL/min)", 15, 100, 75, 5),
            sliderInput("baseline_dsdna", "Baseline Anti-dsDNA (IU/mL)", 50, 2000, 800, 50),
            numericInput("body_weight", "Body Weight (kg)", 65, 40, 120)
          ),
          box(title = "Baseline Biomarker Status", width = 6, status = "warning",
            solidHeader = TRUE,
            tags$h5("Serology", style = "font-weight:bold"),
            fluidRow(
              valueBoxOutput("vb_dsdna",   width = 6),
              valueBoxOutput("vb_C3",      width = 6)
            ),
            fluidRow(
              valueBoxOutput("vb_C4",      width = 6),
              valueBoxOutput("vb_sledai",  width = 6)
            ),
            tags$h5("IFN status", style = "font-weight:bold; margin-top:10px"),
            radioButtons("ifn_status",
                         "IFN Signature",
                         choices = c("IFN-HIGH (>median)" = "high",
                                     "IFN-LOW"             = "low"),
                         selected = "high", inline = TRUE),
            tags$p("IFN-high patients respond better to anifrolumab (TULIP-1/2).",
                   style = "color:#666; font-size:12px"),
            tags$hr(),
            tags$h5("Renal Classification (ISN/RPS)", style = "font-weight:bold"),
            radioButtons("LN_class", "LN Class",
                         choices = c("Class I/II (mesangial)" = "I_II",
                                     "Class III/IV (proliferative)" = "III_IV",
                                     "Class V (membranous)"  = "V",
                                     "No nephritis"          = "none"),
                         selected = "III_IV", inline = FALSE)
          )
        ),
        fluidRow(
          box(title = "SLE Disease Activity Overview", width = 12, status = "primary",
            solidHeader = TRUE,
            fluidRow(
              column(6,
                tags$h5("SLEDAI-2K Domain Breakdown", style = "font-weight:bold"),
                plotOutput("sledai_radar", height = "280px")
              ),
              column(6,
                tags$h5("Baseline Lab Summary", style = "font-weight:bold"),
                DTOutput("baseline_table")
              )
            )
          )
        )
      ),

      # ---- TAB 2: Drug PK ---------------------------------------------------
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Drug Dose Settings", width = 3, status = "success",
            solidHeader = TRUE,
            conditionalPanel("input.use_hcq",
              sliderInput("hcq_dose",  "HCQ dose (mg/day)", 200, 600, 400, 100)),
            conditionalPanel("input.use_beli",
              sliderInput("beli_dose_mgkg", "Belimumab (mg/kg q4w)", 1, 15, 10, 1)),
            conditionalPanel("input.use_anif",
              selectInput("anif_dose", "Anifrolumab dose",
                          choices = c("150 mg IV q4w" = 150, "300 mg IV q4w" = 300),
                          selected = 300)),
            conditionalPanel("input.use_mmf",
              sliderInput("mmf_dose", "MMF dose (mg/day)", 500, 3000, 2000, 250)),
            conditionalPanel("input.use_voc",
              sliderInput("voc_dose", "Voclosporin (mg/day)", 23.7, 71.1, 47.4, 23.7)),
            tags$hr(),
            tags$p(icon("info-circle"), "Belimumab: 10 mg/kg q4w IV (weeks 0,2,4 then monthly); Anifrolumab: 300 mg IV q4w",
                   style = "font-size:11px; color:#666")
          ),
          box(title = "PK Concentration-Time Profiles", width = 9, status = "success",
            solidHeader = TRUE,
            tabsetPanel(
              tabPanel("HCQ",         plotlyOutput("pk_hcq_plot",  height = "380px")),
              tabPanel("Belimumab",   plotlyOutput("pk_beli_plot", height = "380px")),
              tabPanel("Anifrolumab", plotlyOutput("pk_anif_plot", height = "380px")),
              tabPanel("MMF/MPA",     plotlyOutput("pk_mpa_plot",  height = "380px")),
              tabPanel("Voclosporin", plotlyOutput("pk_voc_plot",  height = "380px")),
              tabPanel("All PK",      plotlyOutput("pk_all_plot",  height = "380px"))
            )
          )
        ),
        fluidRow(
          box(title = "PK Parameter Summary", width = 12, status = "info",
            solidHeader = TRUE,
            DTOutput("pk_param_table")
          )
        )
      ),

      # ---- TAB 3: IFN & BAFF -----------------------------------------------
      tabItem(tabName = "ifn",
        fluidRow(
          box(title = "Type I IFN Pathway", width = 6, status = "purple",
            solidHeader = TRUE,
            tags$p("IFN-α drives BAFF upregulation, B cell survival, and amplification of the autoimmune loop. Anifrolumab (anti-IFNAR1) blocks downstream IFN signaling."),
            plotlyOutput("ifn_score_plot", height = "320px"),
            plotlyOutput("ifna_plot",      height = "220px")
          ),
          box(title = "BAFF/BLyS Dynamics (TMDD)", width = 6, status = "warning",
            solidHeader = TRUE,
            tags$p("Belimumab binds free BAFF via TMDD kinetics. IFN-α → BAFF upregulation creates a feedback loop."),
            plotlyOutput("baff_free_plot",    height = "280px"),
            plotlyOutput("baff_complex_plot", height = "260px")
          )
        ),
        fluidRow(
          box(title = "IFN Score vs. Drug Exposure", width = 12, status = "info",
            solidHeader = TRUE,
            fluidRow(
              column(6, plotlyOutput("ifn_vs_anif_dr", height = "320px")),
              column(6,
                tags$h5("IFN Signature High vs Low", style = "font-weight:bold"),
                tags$p("TULIP-2 trial: Anifrolumab 300 mg reduced IFN score by 75% at week 52 in IFN-HIGH patients (vs ~55% IFN-LOW)."),
                tableOutput("ifn_response_table")
              )
            )
          )
        )
      ),

      # ---- TAB 4: Anti-dsDNA & Complement -----------------------------------
      tabItem(tabName = "abs",
        fluidRow(
          box(title = "Anti-dsDNA IgG Dynamics", width = 6, status = "danger",
            solidHeader = TRUE,
            tags$p("Anti-dsDNA titers reflect LLPC pool activity. Belimumab depletes B cells → ↓ plasmablasts → ↓ anti-dsDNA over weeks to months."),
            plotlyOutput("ab_dsdna_plot", height = "320px"),
            sliderInput("target_dsdna", "Target anti-dsDNA threshold (IU/mL)",
                        50, 400, 200, 25)
          ),
          box(title = "Complement C3 & C4", width = 6, status = "primary",
            solidHeader = TRUE,
            tags$p("Complement is consumed by circulating immune complexes. C3/C4 normalization indicates disease control."),
            plotlyOutput("C3_plot", height = "240px"),
            plotlyOutput("C4_plot", height = "240px")
          )
        ),
        fluidRow(
          box(title = "Immune Complex Burden", width = 12, status = "warning",
            solidHeader = TRUE,
            fluidRow(
              column(8, plotlyOutput("IC_plot", height = "300px")),
              column(4,
                tags$h5("Key IC-driven downstream effects:"),
                tags$ul(
                  tags$li("Complement activation → C3a/C5a anaphylatoxins"),
                  tags$li("Glomerular deposition → lupus nephritis"),
                  tags$li("Skin/joint deposition → vasculitis, arthritis"),
                  tags$li("CNS: anti-NMDAR access → neuropsychiatric SLE")
                ),
                tags$hr(),
                verbatimTextOutput("ic_summary_text")
              )
            )
          )
        )
      ),

      # ---- TAB 5: Lupus Nephritis ------------------------------------------
      tabItem(tabName = "nephritis",
        fluidRow(
          valueBoxOutput("vb_proteinuria_current", width = 3),
          valueBoxOutput("vb_egfr_current",        width = 3),
          valueBoxOutput("vb_cr_rate",              width = 3),
          valueBoxOutput("vb_ln_class",             width = 3)
        ),
        fluidRow(
          box(title = "Proteinuria (UPCR) Trajectory", width = 6, status = "warning",
            solidHeader = TRUE,
            tags$p("UPCR < 0.5 mg/mg at week 52 with stable eGFR = Complete Renal Response (AURORA, BLISS-LN criteria)."),
            plotlyOutput("prot_plot",      height = "320px"),
            checkboxInput("show_cr_line",  "Show CR threshold (0.5)", value = TRUE),
            checkboxInput("show_pr_line",  "Show PR threshold (1.0)", value = TRUE)
          ),
          box(title = "eGFR Trajectory", width = 6, status = "danger",
            solidHeader = TRUE,
            tags$p("GFR decline driven by sustained proteinuria; voclosporin + MMF combination shown to preserve GFR (AURORA-1)."),
            plotlyOutput("egfr_plot",      height = "320px"),
            sliderInput("target_egfr",    "eGFR target (mL/min)",
                        30, 90, 60, 5)
          )
        ),
        fluidRow(
          box(title = "Lupus Nephritis Treatment Response", width = 12, status = "success",
            solidHeader = TRUE,
            fluidRow(
              column(6, plotlyOutput("ln_response_barplot", height = "300px")),
              column(6, DTOutput("ln_response_table"))
            )
          )
        )
      ),

      # ---- TAB 6: SLEDAI & Response ----------------------------------------
      tabItem(tabName = "sledai",
        fluidRow(
          box(title = "SLEDAI-2K Dynamic Score", width = 8, status = "primary",
            solidHeader = TRUE,
            plotlyOutput("sledai_plot", height = "380px")
          ),
          box(title = "Response Thresholds", width = 4, status = "info",
            solidHeader = TRUE,
            tags$table(class = "table table-bordered table-sm",
              tags$thead(tags$tr(
                tags$th("Endpoint"), tags$th("Definition")
              )),
              tags$tbody(
                tags$tr(tags$td(tags$b("SRI-4")),
                        tags$td("SLEDAI-2K ≥4pt↓ + no new BILAG A + PGA no worsening")),
                tags$tr(tags$td(tags$b("BICLA")),
                        tags$td("All BILAG A/B resolved; no new A/B; SLEDAI no worse")),
                tags$tr(tags$td(tags$b("LLDAS")),
                        tags$td("SLEDAI-2K ≤4; no active major organ; pred ≤7.5 mg/d")),
                tags$tr(tags$td(tags$b("Remission")),
                        tags$td("SLEDAI-2K = 0 (clinical); pred ≤5 mg/d (DORIS)"))
              )
            ),
            tags$hr(),
            tags$h5("% Time in LLDAS (at Week 52)"),
            verbatimTextOutput("lldas_summary"),
            tags$h5("% Time in Remission"),
            verbatimTextOutput("remission_summary")
          )
        ),
        fluidRow(
          box(title = "Treatment Scenario Comparison", width = 12, status = "warning",
            solidHeader = TRUE,
            fluidRow(
              column(3,
                checkboxGroupInput("scen_compare",
                                   "Select scenarios to compare:",
                                   choices = c("No treatment"    = "none",
                                               "HCQ only"        = "hcq",
                                               "+ MMF"           = "mmf",
                                               "+ Belimumab"     = "beli",
                                               "+ Anifrolumab"   = "anif",
                                               "LN triple"       = "triple"),
                                   selected = c("none", "hcq", "beli", "triple"))
              ),
              column(9, plotlyOutput("scenario_compare_plot", height = "340px"))
            )
          )
        )
      ),

      # ---- TAB 7: Population Variability ------------------------------------
      tabItem(tabName = "popvar",
        fluidRow(
          box(title = "Population Simulation Settings", width = 3, status = "purple",
            solidHeader = TRUE,
            sliderInput("n_subj_pop", "N virtual patients", 50, 500, 200, 50),
            sliderInput("iiv_cv_cl", "Belimumab CL IIV (%CV)", 10, 50, 28, 5),
            sliderInput("iiv_cv_baff", "BAFF synthesis IIV (%CV)", 10, 60, 35, 5),
            sliderInput("iiv_cv_upcr", "Baseline UPCR IIV (%CV)", 10, 70, 40, 5),
            actionButton("run_pop", "Run Population Sim",
                         icon = icon("users"), width = "100%",
                         style = "background:#7B2D8B; color:white")
          ),
          box(title = "Population UPCR Trajectories", width = 9, status = "purple",
            solidHeader = TRUE,
            plotlyOutput("pop_prot_ribbon", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Complete Renal Response Rate (Population)", width = 6, status = "success",
            solidHeader = TRUE,
            plotlyOutput("pop_cr_rate_plot", height = "300px")
          ),
          box(title = "Anti-dsDNA Population Distribution at Week 52", width = 6, status = "warning",
            solidHeader = TRUE,
            plotlyOutput("pop_ab_hist", height = "300px")
          )
        )
      ),

      # ---- TAB 8: Mechanistic Map ------------------------------------------
      tabItem(tabName = "map",
        fluidRow(
          box(title = "SLE QSP Mechanistic Map (161 nodes, 9 pathway clusters)",
              width = 12, status = "purple", solidHeader = TRUE,
            fluidRow(
              column(3,
                checkboxGroupInput("map_clusters",
                                   "Show clusters:",
                                   choices = c("Cell Death & Antigen Release" = "celldeath",
                                               "Innate Immunity & IFN"        = "IFN",
                                               "B Cell Activation"            = "Bcell",
                                               "T Cell Dysregulation"         = "Tcell",
                                               "Complement System"            = "Complement",
                                               "Lupus Nephritis"              = "kidney",
                                               "Multi-Organ Involvement"      = "systemic",
                                               "Drug PK/PD"                   = "drugs",
                                               "Biomarkers & Endpoints"       = "endpoints"),
                                   selected = c("IFN", "Bcell", "drugs", "kidney")),
                tags$hr(),
                tags$p("Full map rendered via Graphviz fdp layout (force-directed). SVG/PNG in the model directory.",
                       style = "font-size:12px; color:#666")
              ),
              column(9,
                uiOutput("map_img"),
                tags$p("Click the image to open full-resolution SVG.",
                       style = "font-size:11px; color:#888; margin-top:5px")
              )
            )
          )
        )
      ),

      # ---- TAB 9: References -----------------------------------------------
      tabItem(tabName = "refs",
        fluidRow(
          box(title = "Selected Key References", width = 12, status = "info",
            solidHeader = TRUE,
            DTOutput("refs_table")
          )
        )
      )

    ) # end tabItems
  )   # end dashboardBody
)     # end dashboardPage

# ==============================================================================
# SERVER
# ==============================================================================

server <- function(input, output, session) {

  # --------------------------------------------------------------------------
  # REACTIVE: build simulation inputs
  # --------------------------------------------------------------------------
  sim_inputs <- reactive({
    list(
      hcq_dose   = ifelse(input$use_hcq,  input$hcq_dose,  0),
      beli_dose  = ifelse(input$use_beli, input$beli_dose_mgkg * 70, 0),
      anif_dose  = ifelse(input$use_anif, as.numeric(input$anif_dose), 0),
      mmf_dose   = ifelse(input$use_mmf,  input$mmf_dose,  0),
      voc_dose   = ifelse(input$use_voc,  input$voc_dose,  0),
      cs_dose    = input$cs_dose,
      sim_days   = input$sim_weeks * 7,
      upcr_sl    = input$baseline_upcr,
      egfr_ss    = input$baseline_egfr,
      anti_sl    = input$baseline_dsdna
    )
  })

  # --------------------------------------------------------------------------
  # REACTIVE: run simulation on button click
  # --------------------------------------------------------------------------
  sim_results <- eventReactive(input$run_sim, {
    withProgress(message = "Running SLE simulation...", {
      inp <- sim_inputs()

      # Build dose events
      ev_list <- list()
      if (inp$hcq_dose > 0) {
        ev_list[["HCQ"]] <- ev(amt = inp$hcq_dose/2, cmt = "HCQ_gut",
                               ii = 0.5, addl = inp$sim_days*2 - 1, time = 0)
      }
      if (inp$beli_dose > 0) {
        ev_list[["BELI"]] <- ev(amt = inp$beli_dose, cmt = "BELI_cent",
                                ii = 28, addl = floor(inp$sim_days/28)-1, time = 0)
      }
      if (inp$anif_dose > 0) {
        ev_list[["ANIF"]] <- ev(amt = inp$anif_dose, cmt = "ANIF_cent",
                                ii = 28, addl = floor(inp$sim_days/28)-1, time = 0)
      }
      if (inp$mmf_dose > 0) {
        ev_list[["MMF"]] <- ev(amt = inp$mmf_dose/2, cmt = "MMF_gut",
                               ii = 0.5, addl = inp$sim_days*2 - 1, time = 0)
      }
      if (inp$voc_dose > 0) {
        ev_list[["VOC"]] <- ev(amt = inp$voc_dose/2, cmt = "VOC_gut",
                               ii = 0.5, addl = inp$sim_days*2 - 1, time = 0)
      }

      dose_ev <- if (length(ev_list) > 0) Reduce(mrgsolve::c, ev_list) else NULL

      mod_upd <- param(sle_mod,
                       CS_dose  = inp$cs_dose,
                       UPCR_SLE = inp$upcr_sl,
                       eGFR_ss  = inp$egfr_ss,
                       Anti_SLE = inp$anti_sl)

      tg <- seq(0, inp$sim_days, by = 1)
      if (!is.null(dose_ev)) {
        out <- mrgsim(mod_upd, events = dose_ev, tgrid = tg) %>% as_tibble()
      } else {
        out <- mrgsim(mod_upd, tgrid = tg) %>% as_tibble()
      }
      out
    })
  })

  # --------------------------------------------------------------------------
  # Value boxes – patient profile
  # --------------------------------------------------------------------------
  output$vb_dsdna <- renderValueBox({
    v <- input$baseline_dsdna
    valueBox(paste(v, "IU/mL"), "Anti-dsDNA",
             icon = icon("vial"),
             color = if (v > 200) "red" else "green")
  })
  output$vb_C3 <- renderValueBox({
    valueBox("0.6 g/L", "C3 (low)", icon = icon("water-ladder"), color = "orange")
  })
  output$vb_C4 <- renderValueBox({
    valueBox("0.10 g/L", "C4 (low)", icon = icon("droplet"), color = "orange")
  })
  output$vb_sledai <- renderValueBox({
    act <- switch(input$disease_activity,
                  mild = 5, moderate = 10, severe = 18)
    valueBox(act, "SLEDAI-2K", icon = icon("chart-bar"),
             color = if (act >= 12) "red" else if (act >= 6) "orange" else "green")
  })

  # --------------------------------------------------------------------------
  # Baseline table
  # --------------------------------------------------------------------------
  output$baseline_table <- renderDT({
    df <- data.frame(
      Parameter   = c("Anti-dsDNA", "C3", "C4", "UPCR", "eGFR",
                       "IFN Score", "SLEDAI-2K"),
      Value       = c(paste(input$baseline_dsdna, "IU/mL"),
                       "0.6 g/L", "0.10 g/L",
                       paste(input$baseline_upcr, "mg/mg"),
                       paste(input$baseline_egfr, "mL/min"),
                       "4.0 (high)", "10"),
      Normal      = c("<100 IU/mL", "0.9-1.8 g/L", "0.16-0.47 g/L",
                       "<0.2 mg/mg", ">90 mL/min", "<2.0", "<4"),
      Status      = c("HIGH", "LOW", "LOW", "HIGH", "REDUCED", "HIGH", "MODERATE")
    )
    datatable(df, options = list(dom = "t", pageLength = 10),
              rownames = FALSE) %>%
      formatStyle("Status",
                  backgroundColor = styleEqual(
                    c("HIGH", "LOW", "REDUCED", "MODERATE", "NORMAL"),
                    c("#FFD0D0", "#FFE4B5", "#FFE4B5", "#FFF0C0", "#D0FFD0")))
  })

  # --------------------------------------------------------------------------
  # PK plots
  # --------------------------------------------------------------------------
  output$pk_hcq_plot <- renderPlotly({
    req(sim_results())
    df <- sim_results()
    p <- ggplot(df, aes(x = time, y = C_HCQ_ng)) +
      geom_line(color = "#2E8B57", linewidth = 1.1) +
      geom_hline(yintercept = c(500, 1000, 1200),
                 linetype = "dashed", color = c("orange", "green4", "red"),
                 alpha = 0.7) +
      labs(title = "HCQ Whole Blood Concentration",
           x = "Time (days)", y = "HCQ (ng/mL)") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_beli_plot <- renderPlotly({
    req(sim_results())
    df <- sim_results()
    p <- ggplot(df, aes(x = time, y = C_BELI_ug)) +
      geom_line(color = "#1565C0", linewidth = 1.1) +
      labs(title = "Belimumab Serum Concentration", x = "Time (days)", y = "µg/mL") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_anif_plot <- renderPlotly({
    req(sim_results())
    df <- sim_results()
    p <- ggplot(df, aes(x = time, y = C_ANIF_ug)) +
      geom_line(color = "#6A0DAD", linewidth = 1.1) +
      labs(title = "Anifrolumab Serum Concentration", x = "Time (days)", y = "µg/mL") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_mpa_plot <- renderPlotly({
    req(sim_results())
    df <- sim_results()
    p <- ggplot(df, aes(x = time, y = C_MPA_ug)) +
      geom_line(color = "#C62828", linewidth = 1.1) +
      geom_hline(yintercept = c(1, 3.5), linetype = "dashed",
                 color = c("orange", "green4")) +
      labs(title = "MPA (Mycophenolic Acid)", x = "Time (days)", y = "MPA µg/mL") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_voc_plot <- renderPlotly({
    req(sim_results())
    df <- sim_results()
    p <- ggplot(df, aes(x = time, y = C_VOC_ng)) +
      geom_line(color = "#E65100", linewidth = 1.1) +
      labs(title = "Voclosporin Plasma Concentration", x = "Time (days)", y = "ng/mL") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_all_plot <- renderPlotly({
    req(sim_results())
    df <- sim_results() %>%
      select(time, C_HCQ_ng, C_BELI_ug, C_ANIF_ug, C_MPA_ug, C_VOC_ng) %>%
      pivot_longer(-time, names_to = "drug", values_to = "conc") %>%
      mutate(drug = recode(drug,
        C_HCQ_ng  = "HCQ (ng/mL)",
        C_BELI_ug = "Belimumab (µg/mL)",
        C_ANIF_ug = "Anifrolumab (µg/mL)",
        C_MPA_ug  = "MPA (µg/mL)",
        C_VOC_ng  = "Voclosporin (ng/mL)"))
    p <- ggplot(df %>% filter(conc > 0.001), aes(x = time, y = conc, color = drug)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~ drug, scales = "free_y", ncol = 2) +
      labs(title = "All Drug PK Profiles", x = "Time (days)", y = "Concentration") +
      theme_bw() + theme(legend.position = "none")
    ggplotly(p)
  })

  # PK parameter table
  output$pk_param_table <- renderDT({
    df <- data.frame(
      Drug          = c("HCQ", "Belimumab", "Anifrolumab", "MMF/MPA", "Voclosporin"),
      Route         = c("Oral BID", "IV q4w", "IV q4w", "Oral BID", "Oral BID"),
      t_half        = c("~50 days", "~19.4 days", "~26 days", "~17 h (MPA)", "~30 h"),
      Vd            = c("~5200 L", "5.29 L (V1)", "6.8 L", "252 L (MPA)", "2300 L"),
      CL            = c("~100 L/day", "0.215 L/day", "0.170 L/day", "252 L/day", "1512 L/day"),
      MOA           = c("TLR7/9 block (lysosomotropic)",
                        "Anti-BAFF/BLyS (TMDD)",
                        "Anti-IFNAR1 → ↓ISG",
                        "IMPDH inhibition → ↓purines",
                        "Calcineurin inhibition → ↓NFAT")
    )
    datatable(df, options = list(dom = "t"), rownames = FALSE)
  })

  # --------------------------------------------------------------------------
  # IFN plots
  # --------------------------------------------------------------------------
  output$ifn_score_plot <- renderPlotly({
    req(sim_results())
    df <- sim_results()
    p <- ggplot(df, aes(x = time, y = IFN_norm)) +
      geom_line(color = "#6A0DAD", linewidth = 1.2) +
      geom_hline(yintercept = c(2, 5), linetype = "dashed",
                 color = c("green4", "red"), alpha = 0.7) +
      labs(title = "IFN Gene Signature Score", x = "Days", y = "IFN Score (normalized)") +
      theme_bw()
    ggplotly(p)
  })

  output$ifna_plot <- renderPlotly({
    req(sim_results())
    p <- ggplot(sim_results(), aes(x = time, y = IFNa_conc)) +
      geom_line(color = "#AB82FF", linewidth = 1) +
      labs(title = "IFN-α Concentration", x = "Days", y = "IFN-α (IU/mL, relative)") +
      theme_bw()
    ggplotly(p)
  })

  output$baff_free_plot <- renderPlotly({
    req(sim_results())
    p <- ggplot(sim_results(), aes(x = time, y = BAFF_free)) +
      geom_line(color = "#006400", linewidth = 1.2) +
      geom_hline(yintercept = 2.8, linetype = "dashed", color = "gray50") +
      labs(title = "Free BAFF/BLyS (TMDD)", x = "Days", y = "Free BAFF (ng/mL)") +
      theme_bw()
    ggplotly(p)
  })

  output$baff_complex_plot <- renderPlotly({
    req(sim_results())
    p <- ggplot(sim_results(), aes(x = time, y = BELI_cmplx)) +
      geom_line(color = "#228B22", linewidth = 1) +
      labs(title = "Belimumab-BAFF Complex (TMDD)", x = "Days", y = "Complex (AU)") +
      theme_bw()
    ggplotly(p)
  })

  output$ifn_response_table <- renderTable({
    data.frame(
      Population   = c("IFN-HIGH", "IFN-LOW"),
      `BICLA (Week 52)` = c("47.8%", "32.1%"),
      `IFN Score ↓` = c("74.6%", "54.9%"),
      `Flare Rate ↓` = c("49.5%", "30.1%"),
      check.names = FALSE
    )
  })

  # --------------------------------------------------------------------------
  # Anti-dsDNA & Complement
  # --------------------------------------------------------------------------
  output$ab_dsdna_plot <- renderPlotly({
    req(sim_results())
    df <- sim_results()
    p <- ggplot(df, aes(x = time, y = AntiDsDNA)) +
      geom_line(color = "#C62828", linewidth = 1.2) +
      geom_hline(yintercept = input$target_dsdna,
                 linetype = "dashed", color = "green4") +
      labs(title = "Anti-dsDNA IgG Titer", x = "Days", y = "Anti-dsDNA (IU/mL)") +
      theme_bw()
    ggplotly(p)
  })

  output$C3_plot <- renderPlotly({
    req(sim_results())
    p <- ggplot(sim_results(), aes(x = time, y = C3_serum)) +
      geom_line(color = "#1565C0", linewidth = 1.1) +
      geom_hline(yintercept = c(0.9, 1.8), linetype = "dashed",
                 color = c("red", "gray60")) +
      labs(title = "C3 Serum Level", x = "Days", y = "C3 (g/L)") + ylim(0, 2) +
      theme_bw()
    ggplotly(p)
  })

  output$C4_plot <- renderPlotly({
    req(sim_results())
    p <- ggplot(sim_results(), aes(x = time, y = C4_serum)) +
      geom_line(color = "#0D47A1", linewidth = 1.1) +
      geom_hline(yintercept = 0.16, linetype = "dashed", color = "red") +
      labs(title = "C4 Serum Level", x = "Days", y = "C4 (g/L)") + ylim(0, 0.5) +
      theme_bw()
    ggplotly(p)
  })

  output$IC_plot <- renderPlotly({
    req(sim_results())
    p <- ggplot(sim_results(), aes(x = time, y = IC_burden)) +
      geom_line(color = "#E65100", linewidth = 1.2) +
      labs(title = "Circulating Immune Complex Burden",
           x = "Days", y = "IC Burden (AU)") +
      theme_bw()
    ggplotly(p)
  })

  output$ic_summary_text <- renderPrint({
    req(sim_results())
    df <- sim_results()
    last <- tail(df, 1)
    cat(sprintf("Week %d Summary:\n", input$sim_weeks),
        sprintf("  IC burden:  %.2f AU\n", last$IC_burden),
        sprintf("  C3 serum:   %.2f g/L\n", last$C3_serum),
        sprintf("  C4 serum:   %.2f g/L\n", last$C4_serum),
        sprintf("  Anti-dsDNA: %.0f IU/mL\n", last$AntiDsDNA))
  })

  # --------------------------------------------------------------------------
  # Lupus Nephritis
  # --------------------------------------------------------------------------
  output$vb_proteinuria_current <- renderValueBox({
    req(sim_results())
    v <- tail(sim_results()$Proteinuria, 1)
    valueBox(sprintf("%.2f mg/mg", v), "UPCR at Week 52",
             icon = icon("flask"),
             color = if (v < 0.5) "green" else if (v < 1.5) "yellow" else "red")
  })
  output$vb_egfr_current <- renderValueBox({
    req(sim_results())
    v <- tail(sim_results()$eGFR_ode, 1)
    valueBox(sprintf("%.0f mL/min", v), "eGFR at Week 52",
             icon = icon("filter"),
             color = if (v > 60) "green" else if (v > 30) "orange" else "red")
  })
  output$vb_cr_rate <- renderValueBox({
    req(sim_results())
    v <- tail(sim_results()$CR_r, 1) * 100
    valueBox(sprintf("%.0f%%", v), "Complete Renal Response",
             icon = icon("check-circle"),
             color = if (v > 0) "green" else "red")
  })
  output$vb_ln_class <- renderValueBox({
    cls <- switch(input$LN_class,
                  I_II  = "Class I/II",
                  III_IV = "Class III/IV",
                  V     = "Class V",
                  none  = "No nephritis")
    valueBox(cls, "LN Classification", icon = icon("microscope"), color = "purple")
  })

  output$prot_plot <- renderPlotly({
    req(sim_results())
    df <- sim_results()
    p <- ggplot(df, aes(x = time, y = Proteinuria)) +
      geom_line(color = "#E65100", linewidth = 1.2) +
      {if (input$show_cr_line) geom_hline(yintercept = 0.5, linetype = "dashed",
                                           color = "green4") else NULL} +
      {if (input$show_pr_line) geom_hline(yintercept = 1.0, linetype = "dotted",
                                           color = "orange") else NULL} +
      labs(title = "Proteinuria (UPCR)", x = "Days", y = "UPCR (mg/mg)") +
      theme_bw()
    ggplotly(p)
  })

  output$egfr_plot <- renderPlotly({
    req(sim_results())
    df <- sim_results()
    p <- ggplot(df, aes(x = time, y = eGFR_ode)) +
      geom_line(color = "#1565C0", linewidth = 1.2) +
      geom_hline(yintercept = c(60, input$target_egfr),
                 linetype = c("dashed", "dotted"),
                 color = c("orange", "green4")) +
      labs(title = "eGFR Over Time", x = "Days", y = "eGFR (mL/min/1.73m²)") +
      ylim(0, 100) + theme_bw()
    ggplotly(p)
  })

  output$ln_response_table <- renderDT({
    df <- data.frame(
      Regimen  = c("SOC alone (CS+CYC)", "MMF + CS", "Voclosporin + MMF + CS",
                   "Belimumab + MMF + CS", "Triple (Voc+MMF+Beli+CS)"),
      CR_Rate  = c("22%", "26%", "41%", "43%", "~55% (est.)"),
      PR_Rate  = c("30%", "35%", "49%", "53%", "~65% (est.)"),
      Trial    = c("—", "ALMS", "AURORA-1", "BLISS-LN", "Exploratory"),
      stringsAsFactors = FALSE
    )
    datatable(df, options = list(dom = "t"), rownames = FALSE)
  })

  # --------------------------------------------------------------------------
  # SLEDAI
  # --------------------------------------------------------------------------
  output$sledai_plot <- renderPlotly({
    req(sim_results())
    df <- sim_results()
    p <- ggplot(df, aes(x = time, y = SLEDAI_c)) +
      geom_line(color = "#4A0080", linewidth = 1.3) +
      geom_hline(yintercept = c(4, 6, 12),
                 linetype = "dashed",
                 color = c("green4", "orange", "red"),
                 alpha = 0.8) +
      annotate("text", x = max(df$time)*0.85, y = 4.5,  label = "LLDAS (≤4)",
               color = "green4", size = 3.2) +
      annotate("text", x = max(df$time)*0.85, y = 6.7,  label = "Low activity",
               color = "orange", size = 3.2) +
      annotate("text", x = max(df$time)*0.85, y = 12.7, label = "High activity",
               color = "red", size = 3.2) +
      labs(title = "SLEDAI-2K Dynamic Score",
           x = "Days", y = "SLEDAI-2K") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$lldas_summary <- renderPrint({
    req(sim_results())
    df <- sim_results() %>% filter(time > input$sim_weeks * 7 * 0.25)
    pct <- mean(df$SLEDAI_c <= 4, na.rm = TRUE) * 100
    cat(sprintf("%.1f%% of time in LLDAS (post-stabilization)\n", pct))
  })

  output$remission_summary <- renderPrint({
    req(sim_results())
    df <- sim_results() %>% filter(time > input$sim_weeks * 7 * 0.25)
    pct <- mean(df$SLEDAI_c == 0, na.rm = TRUE) * 100
    cat(sprintf("%.1f%% of time in clinical remission\n", pct))
  })

  # --------------------------------------------------------------------------
  # Mechanistic map
  # --------------------------------------------------------------------------
  output$map_img <- renderUI({
    # Assume SVG file is at relative path ../../sle_qsp.svg from app.R
    svg_path <- "../sle_qsp.svg"
    png_path <- "../sle_qsp.png"
    tags$div(
      tags$a(href = svg_path, target = "_blank",
        tags$img(src = png_path,
                 style = "max-width:100%; border:1px solid #DDD; border-radius:4px;",
                 alt = "SLE QSP Mechanistic Map"))
    )
  })

  # --------------------------------------------------------------------------
  # References table
  # --------------------------------------------------------------------------
  output$refs_table <- renderDT({
    df <- data.frame(
      No  = 1:12,
      Reference = c(
        "Tsokos GC. N Engl J Med 2011",
        "Crow MK. Ann Rheum Dis 2023",
        "Baechler EC, et al. PNAS 2003 (IFN signature)",
        "Lood C, et al. Nat Med 2016 (NETs & IFN)",
        "Navarra SV (BLISS-52). Lancet 2011",
        "Furie R (BLISS-76). Arthritis Rheum 2011",
        "Morand EF (TULIP-2). N Engl J Med 2020",
        "Rovin BH (AURORA-1). Lancet 2021",
        "Furie RA (BLISS-LN). N Engl J Med 2020",
        "Tett SE. Br J Clin Pharmacol 1988 (HCQ PK)",
        "Ding J, et al. CPT PSP 2020 (Belimumab QSP)",
        "Forde EA, et al. J PK/PD 2021 (SLE QSP)"
      ),
      PubMed_Link = c(
        "https://pubmed.ncbi.nlm.nih.gov/22129253/",
        "https://pubmed.ncbi.nlm.nih.gov/37279922/",
        "https://pubmed.ncbi.nlm.nih.gov/12604793/",
        "https://pubmed.ncbi.nlm.nih.gov/26779811/",
        "https://pubmed.ncbi.nlm.nih.gov/21296403/",
        "https://pubmed.ncbi.nlm.nih.gov/22127708/",
        "https://pubmed.ncbi.nlm.nih.gov/31851795/",
        "https://pubmed.ncbi.nlm.nih.gov/34003766/",
        "https://pubmed.ncbi.nlm.nih.gov/32937045/",
        "https://pubmed.ncbi.nlm.nih.gov/3190987/",
        "https://pubmed.ncbi.nlm.nih.gov/32255277/",
        "https://pubmed.ncbi.nlm.nih.gov/33856624/"
      ),
      stringsAsFactors = FALSE
    )
    datatable(df, escape = FALSE,
              options = list(dom = "ft", pageLength = 12),
              rownames = FALSE) %>%
      formatStyle("No", width = "40px")
  })

}  # end server

# ==============================================================================
# LAUNCH
# ==============================================================================
shinyApp(ui, server)
