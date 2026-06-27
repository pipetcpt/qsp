################################################################################
# aHUS QSP Shiny App
# Atypical Hemolytic Uremic Syndrome
# Interactive Quantitative Systems Pharmacology Dashboard
################################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)

# ─── MODEL CODE (embedded) ────────────────────────────────────────────────────

code <- '
$PARAM @annotated
CL    : 0.31  : Clearance (L/day)
V1    : 5.30  : Central volume (L)
Q     : 0.54  : Inter-compartmental clearance (L/day)
V2    : 4.10  : Peripheral volume (L)
kon   : 0.864  : C5-drug association (nM^-1 day^-1)
koff  : 8.64e-3: C5-drug dissociation (day^-1)
kint  : 0.02   : TMDD internalization (day^-1)
C3_ss   : 8000  : C3 steady state (nM)
kC3syn  : 4.167 : C3 synthesis (nM/h)
kC3deg  : 5.208e-4 : C3 degradation (h^-1)
kAP     : 0.12  : AP activation rate (h^-1)
kAP_amp : 0.30  : AP amplification coefficient
kC3b_deg: 0.80  : C3b inactivation rate (h^-1)
C5_ss   : 395   : C5 steady state (nM)
kC5syn  : 11.85 : C5 synthesis (nM/day)
kC5deg  : 0.030 : C5 degradation (day^-1)
kC5conv : 0.015 : C5 convertase activity
kMACform: 0.25  : MAC formation rate (day^-1)
kMACclr : 0.50  : MAC clearance rate (day^-1)
CFH_factor: 0.30 : CFH function (0-1)
CFI_factor: 0.50 : CFI function (0-1)
CD46_factor:0.70 : CD46 function (0-1)
kEI_MAC  : 0.40 : Endothelial injury per MAC
kEI_rep  : 0.15 : Endothelial repair rate
EI_max   : 10.0 : Max endothelial injury score
kPLT_loss: 0.08 : Platelet consumption rate
kPLT_prod: 13.0 : Platelet production (10^9/L/day)
PLT_ss   : 250  : Normal platelets (10^9/L)
kHgb_loss: 0.04 : Hgb loss rate
kHgb_prod: 0.20 : Hgb production (g/dL/day)
Hgb_ss   : 14.0 : Normal Hgb (g/dL)
kLDH_rel : 15.0 : LDH release rate
kLDH_clr : 0.40 : LDH clearance (day^-1)
LDH_norm : 200  : Normal LDH ULN (U/L)
kHpg_cons: 0.35 : Haptoglobin consumption
kHpg_syn : 0.15 : Haptoglobin synthesis
Hpg_ss   : 1.20 : Normal haptoglobin (g/L)
GFR_ss    : 90.0 : Baseline GFR (mL/min/1.73m2)
kGFR_loss : 0.05 : GFR loss rate
kGFR_rep  : 0.10 : GFR recovery rate
sC5b9_base: 150  : Baseline sC5b-9 (ng/mL)
sC5b9_max : 3500 : Max sC5b-9 (ng/mL)
kSCT_form : 0.60 : Schistocyte formation
kSCT_clr  : 0.25 : Schistocyte clearance
kCRP_syn  : 12.0 : CRP synthesis
kCRP_clr  : 0.60 : CRP clearance (day^-1)
scenario  : 2    : Treatment scenario (1-5)

$CMT
Drug_C Drug_P Drug_C5 C3pool C3b_AP C5free C5conv MACflux
Endo_inj PLT Hgb LDH Hpg GFR Schist CRP sC5b9 CH50pct

$MAIN
if(NEWIND <= 1) {
  Drug_C_0=0; Drug_P_0=0; Drug_C5_0=0;
  C3pool_0=C3_ss; C3b_AP_0=200; C5free_0=C5_ss;
  C5conv_0=5.0; MACflux_0=2.5;
  Endo_inj_0=6.0; PLT_0=60; Hgb_0=7.5; LDH_0=1500;
  Hpg_0=0.05; GFR_0=25; Schist_0=4.5; CRP_0=45;
  sC5b9_0=sC5b9_max*0.8; CH50pct_0=60;
}
double reg_eff = CFH_factor * CFI_factor * CD46_factor;
double Cc_mgL = Drug_C / V1;
double Cc_nM  = Cc_mgL * 1000.0 / 148.0;
double fB_inhibition = (scenario==5) ? 0.95 : 0.0;
double fD_inhibition = (scenario==4) ? 0.80 : 0.0;
double AP_activity = kAP*(1.0-fB_inhibition)*(1.0-fD_inhibition)/reg_eff;

$ODE
double TMDD_bind = kon * Cc_nM * C5free - koff * Drug_C5;
dxdt_Drug_C = -(CL/V1)*Drug_C - (Q/V1)*Drug_C + (Q/V2)*Drug_P
              - kon*(Drug_C/V1)*C5free*V1 + koff*Drug_C5;
dxdt_Drug_P = (Q/V1)*Drug_C - (Q/V2)*Drug_P;
dxdt_Drug_C5= kon*(Drug_C/V1)*C5free*148.0*V1/1000.0
              - koff*Drug_C5 - kint*Drug_C5;

double C3_syn_rate=kC3syn*24; double kC3deg_d=kC3deg*24;
dxdt_C3pool = C3_syn_rate - kC3deg_d*C3pool
              - AP_activity*C3pool/(C3_ss/100.0);

double C3b_form=AP_activity*C3pool/C3_ss;
double C3b_inac=kC3b_deg*reg_eff*C3b_AP;
double amploop=kAP_amp*(C3b_AP/C3_ss);
dxdt_C3b_AP=C3b_form+amploop*C3pool/C3_ss-C3b_inac;

double C5_cleavage=kC5conv*C5conv*C5free;
double C5_bound_loss=kon*Cc_nM*C5free-koff*(Drug_C5/148.0*1000.0/V1);
dxdt_C5free=kC5syn-kC5deg*C5free-C5_cleavage-C5_bound_loss;

double C5conv_form=0.08*pow(C3b_AP/C3_ss,2.0);
double C5conv_clr=0.15*C5conv;
dxdt_C5conv=C5conv_form-C5conv_clr;

double MAC_form=kMACform*C5_cleavage;
double MAC_clr=kMACclr*MACflux;
double CD59_eff=0.60;
dxdt_MACflux=MAC_form*(1.0-CD59_eff)-MAC_clr;

double EI_drive=kEI_MAC*MACflux;
double EI_rep=kEI_rep*Endo_inj;
dxdt_Endo_inj=EI_drive-EI_rep;
if(Endo_inj>EI_max) dxdt_Endo_inj=0;
if(Endo_inj<0) dxdt_Endo_inj=0;

double PLT_loss=kPLT_loss*PLT*(Endo_inj/EI_max);
dxdt_PLT=kPLT_prod-PLT_loss;
if(PLT<5) dxdt_PLT=0;

double Hgb_loss_rate=kHgb_loss*Hgb*(Endo_inj/EI_max);
dxdt_Hgb=kHgb_prod-Hgb_loss_rate;
if(Hgb<4.0) dxdt_Hgb=0;

double hemolysis_rate=kHgb_loss*Hgb*(Endo_inj/EI_max);
dxdt_LDH=kLDH_rel*hemolysis_rate-kLDH_clr*LDH;
dxdt_Hpg=kHpg_syn-kHpg_cons*Hgb_loss_rate-0.10*Hpg;
if(Hpg<0) dxdt_Hpg=0;

double GFR_loss_rate=kGFR_loss*GFR*(Endo_inj/6.0);
double GFR_rep_rate=kGFR_rep*(GFR_ss-GFR)*(1.0-Endo_inj/EI_max);
dxdt_GFR=GFR_rep_rate-GFR_loss_rate;
if(GFR<5) dxdt_GFR=0;

dxdt_Schist=kSCT_form*(Endo_inj/EI_max)-kSCT_clr*Schist;
if(Schist<0) dxdt_Schist=0;
dxdt_CRP=kCRP_syn*MACflux-kCRP_clr*CRP;

double sC5b9_drive=MAC_form*(sC5b9_max/3.0);
dxdt_sC5b9=sC5b9_drive-0.20*sC5b9;

double C5_pct=C5free/C5_ss;
dxdt_CH50pct=0.50*(100.0*C5_pct-CH50pct);

$TABLE
double Cc_ugmL=Drug_C/V1;
double Cc_nM_out=Cc_ugmL*1000.0/148.0;
double LDH_xULN=LDH/LDH_norm;
double C5_block=1.0-C5free/C5_ss;
double TMA_flag=(PLT<150)?1.0:0.0;
double dialysis_risk=(GFR<15)?1.0:0.0;
capture Cc_ugmL Cc_nM_out LDH_xULN C5_block TMA_flag dialysis_risk;
capture C3pool C3b_AP C5free MACflux sC5b9 CH50pct;
capture Drug_C Drug_P Drug_C5 C3pool C3b_AP C5free C5conv MACflux;
capture Endo_inj PLT Hgb LDH Hpg GFR Schist CRP sC5b9 CH50pct;

$CAPTURE
Cc_ugmL Cc_nM_out LDH_xULN C5_block TMA_flag dialysis_risk
C3pool C3b_AP C5free MACflux sC5b9 CH50pct
Drug_C Drug_P Drug_C5 C5conv Endo_inj PLT Hgb LDH Hpg GFR Schist CRP
'

mod <- mcode("aHUS_Shiny", code, quiet=TRUE)

# Helper: build dosing event
build_ev <- function(scenario, dose_ecu=900, dose_maint=1200, start=1) {
  if(scenario %in% c(2,4)) {
    data.frame(
      time=c(start, start+7, start+14, start+21,
             seq(start+28, by=14, length.out=12)),
      amt=c(rep(dose_ecu,4), rep(dose_maint,12)),
      cmt=1, evid=1, rate=-2, ii=0, addl=0
    )
  } else if(scenario==3) {
    data.frame(
      time=c(start, start+15, seq(start+29, by=56, length.out=5)),
      amt=c(2400, 3000, rep(3300,5)),
      cmt=1, evid=1, rate=-2, ii=0, addl=0
    )
  } else {
    data.frame(time=0, amt=0, cmt=1, evid=0, rate=0, ii=0, addl=0)
  }
}

run_sim <- function(scenario, params, sim_days=365) {
  ev <- build_ev(scenario,
                 dose_ecu   = params$dose_ecu,
                 dose_maint = params$dose_maint)
  mod %>%
    param(
      scenario    = scenario,
      CFH_factor  = params$CFH_factor,
      CFI_factor  = params$CFI_factor,
      CD46_factor = params$CD46_factor,
      kAP         = params$kAP,
      kEI_MAC     = params$kEI_MAC,
      kEI_rep     = params$kEI_rep,
      kGFR_rep    = params$kGFR_rep
    ) %>%
    ev(ev) %>%
    mrgsim(end=sim_days, delta=1) %>%
    as.data.frame()
}

# ─── UI ────────────────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = "aHUS QSP Dashboard",
    titleWidth = 300
  ),

  dashboardSidebar(
    width = 300,
    sidebarMenu(
      id = "tabs",
      menuItem("Overview",         tabName="overview",   icon=icon("info-circle")),
      menuItem("Patient Profile",  tabName="patient",    icon=icon("user-md")),
      menuItem("PK Analysis",      tabName="pk",         icon=icon("pills")),
      menuItem("Complement PD",    tabName="complement", icon=icon("dna")),
      menuItem("TMA & Hematology", tabName="hema",       icon=icon("tint")),
      menuItem("Renal Outcomes",   tabName="renal",      icon=icon("filter")),
      menuItem("Scenario Compare", tabName="scenarios",  icon=icon("chart-bar")),
      menuItem("Biomarkers",       tabName="biomarkers", icon=icon("flask")),
      menuItem("References",       tabName="refs",       icon=icon("book"))
    ),
    tags$hr(),

    # Treatment scenario
    selectInput("scenario", "Treatment Scenario:",
      choices=list(
        "Natural History (no treatment)"    = 1,
        "Eculizumab Standard Dosing"        = 2,
        "Ravulizumab (q8w)"                 = 3,
        "Eculizumab + Danicopan (fDi)"      = 4,
        "Iptacopan (oral Factor B inhibitor)" = 5
      ), selected=2
    ),

    sliderInput("sim_days", "Simulation Duration (days):", 30, 730, 365, step=30),
    sliderInput("dose_ecu",   "Eculizumab Induction Dose (mg):", 300, 1800, 900, step=100),
    sliderInput("dose_maint", "Eculizumab Maintenance Dose (mg):", 600, 2400, 1200, step=100),

    tags$hr(),
    tags$b("Genetic Profile (Complement Regulation):"),
    sliderInput("CFH_factor",  "CFH Function (0=loss, 1=normal):", 0.0, 1.0, 0.30, step=0.05),
    sliderInput("CFI_factor",  "CFI Function:",                     0.0, 1.0, 0.50, step=0.05),
    sliderInput("CD46_factor", "CD46/MCP Function:",                0.0, 1.0, 0.70, step=0.05),
    sliderInput("kAP", "AP Activation Rate (dysreg level):", 0.01, 0.50, 0.12, step=0.01),

    tags$hr(),
    tags$b("Disease/Organ Parameters:"),
    sliderInput("kEI_MAC", "Endothelial Injury per MAC:", 0.1, 1.0, 0.40, step=0.05),
    sliderInput("kEI_rep", "Endothelial Repair Rate:",    0.05, 0.50, 0.15, step=0.05),
    sliderInput("kGFR_rep","Renal Recovery Rate:",        0.02, 0.50, 0.10, step=0.02),

    tags$hr(),
    actionButton("run", "Run Simulation", class="btn btn-success btn-block",
                 icon=icon("play"))
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .info-box-icon { background-color: rgba(0,0,0,0.2) !important; }
      .box { border-top-color: #2C3E50; }
    "))),

    tabItems(

      # ── Tab 1: Overview ──────────────────────────────────────────────────────
      tabItem("overview",
        fluidRow(
          box(width=12, status="primary", title="aHUS Disease Overview & QSP Framework",
            solidHeader=TRUE, collapsible=FALSE,
            fluidRow(
              column(6,
                h4("What is aHUS?"),
                p("Atypical Hemolytic Uremic Syndrome (aHUS) is a rare, life-threatening
                   thrombotic microangiopathy (TMA) caused by chronic, uncontrolled activation
                   of the complement alternative pathway (AP). Unlike typical HUS (STEC-associated),
                   aHUS is characterized by:"),
                tags$ul(
                  tags$li("Mutations in complement regulatory genes (CFH, CFI, CD46, C3, CFB, THBD)"),
                  tags$li("Anti-CFH autoantibodies (~10% of cases)"),
                  tags$li("Triad: thrombocytopenia + MAHA + acute kidney injury"),
                  tags$li("Normal ADAMTS13 activity (unlike TTP)"),
                  tags$li("High risk of ESRD (>50%), recurrence (40-80%), mortality (25%)")
                ),
                h4("Mechanistic Cascade"),
                tags$ol(
                  tags$li("AP dysregulation → persistent C3b generation"),
                  tags$li("C3b₂Bb (AP C5 convertase) → C5 cleavage"),
                  tags$li("C5a → inflammatory amplification (neutrophil/EC activation)"),
                  tags$li("C5b-9 (MAC) → endothelial injury → TMA"),
                  tags$li("Platelet consumption + RBC fragmentation + renal damage")
                )
              ),
              column(6,
                h4("Drug Targets in this Model"),
                tags$table(class="table table-bordered table-sm",
                  tags$thead(tags$tr(tags$th("Drug"), tags$th("Target"), tags$th("Route"), tags$th("Dosing"))),
                  tags$tbody(
                    tags$tr(tags$td("Eculizumab"), tags$td("C5"), tags$td("IV"), tags$td("900mg/wk×4, 1200mg/q2w")),
                    tags$tr(tags$td("Ravulizumab"), tags$td("C5 (recycling)"), tags$td("IV"), tags$td("q8w")),
                    tags$tr(tags$td("Iptacopan"), tags$td("Factor B"), tags$td("Oral"), tags$td("200mg BID")),
                    tags$tr(tags$td("Danicopan"), tags$td("Factor D"), tags$td("Oral"), tags$td("Add-on")),
                    tags$tr(tags$td("Avacopan"), tags$td("C5aR1"), tags$td("Oral"), tags$td("30mg BID"))
                  )
                ),
                h4("Model Compartments (18 ODEs)"),
                p("PK (3) + Complement cascade (5) + TMA/hematology (8) + Renal/inflammation (2)")
              )
            )
          )
        ),
        fluidRow(
          valueBoxOutput("box_PLT",   width=3),
          valueBoxOutput("box_GFR",   width=3),
          valueBoxOutput("box_LDH",   width=3),
          valueBoxOutput("box_CH50",  width=3)
        )
      ),

      # ── Tab 2: Patient Profile ───────────────────────────────────────────────
      tabItem("patient",
        fluidRow(
          box(width=6, status="warning", title="Complement Regulatory Status",
            solidHeader=TRUE,
            plotlyOutput("plot_cfh_radar", height="300px"),
            tags$hr(),
            h5("Regulatory Effectiveness Score:"),
            verbatimTextOutput("reg_score")
          ),
          box(width=6, status="danger", title="Clinical Presentation at Diagnosis",
            solidHeader=TRUE,
            plotlyOutput("plot_presentation", height="350px")
          )
        ),
        fluidRow(
          box(width=12, status="info", title="AP Dysregulation Cascade (at t=0)",
            solidHeader=TRUE,
            plotlyOutput("plot_complement_cascade", height="300px")
          )
        )
      ),

      # ── Tab 3: PK Analysis ───────────────────────────────────────────────────
      tabItem("pk",
        fluidRow(
          box(width=6, status="primary", title="Eculizumab PK (Concentration-Time)",
            solidHeader=TRUE,
            plotlyOutput("plot_pk_conc", height="350px")
          ),
          box(width=6, status="primary", title="TMDD: Free vs Drug-Bound C5",
            solidHeader=TRUE,
            plotlyOutput("plot_c5_tmdd", height="350px")
          )
        ),
        fluidRow(
          box(width=6, status="info", title="C5 Blockade (Fractional Inhibition)",
            solidHeader=TRUE,
            plotlyOutput("plot_c5block", height="300px")
          ),
          box(width=6, status="info", title="PK Parameters (2-Compartment Model)",
            solidHeader=TRUE,
            tags$table(class="table table-bordered",
              tags$thead(tags$tr(
                tags$th("Parameter"), tags$th("Value"), tags$th("Unit"), tags$th("Source")
              )),
              tags$tbody(
                tags$tr(tags$td("CL"),  tags$td("0.31"), tags$td("L/day"), tags$td("Menne 2015")),
                tags$tr(tags$td("V1"),  tags$td("5.30"), tags$td("L"),     tags$td("Menne 2015")),
                tags$tr(tags$td("Q"),   tags$td("0.54"), tags$td("L/day"), tags$td("Menne 2015")),
                tags$tr(tags$td("V2"),  tags$td("4.10"), tags$td("L"),     tags$td("Menne 2015")),
                tags$tr(tags$td("t½"),  tags$td("~11"), tags$td("days"),   tags$td("Calculated")),
                tags$tr(tags$td("Kd (C5)"), tags$td("~10"), tags$td("pM"), tags$td("Thomas 2012")),
                tags$tr(tags$td("kon"), tags$td("0.864"), tags$td("nM⁻¹d⁻¹"), tags$td("Estimated")),
                tags$tr(tags$td("koff"),tags$td("8.64e-3"), tags$td("d⁻¹"), tags$td("Estimated"))
              )
            )
          )
        )
      ),

      # ── Tab 4: Complement PD ─────────────────────────────────────────────────
      tabItem("complement",
        fluidRow(
          box(width=6, status="info", title="C3 Pool Dynamics",
            solidHeader=TRUE,
            plotlyOutput("plot_c3", height="300px")
          ),
          box(width=6, status="danger", title="MAC Flux (Disease Activity Driver)",
            solidHeader=TRUE,
            plotlyOutput("plot_mac", height="300px")
          )
        ),
        fluidRow(
          box(width=6, status="warning", title="C3b Alternative Pathway Accumulation",
            solidHeader=TRUE,
            plotlyOutput("plot_c3b", height="300px")
          ),
          box(width=6, status="success", title="C5 Convertase Activity",
            solidHeader=TRUE,
            plotlyOutput("plot_c5conv", height="300px")
          )
        )
      ),

      # ── Tab 5: TMA & Hematology ──────────────────────────────────────────────
      tabItem("hema",
        fluidRow(
          box(width=6, status="danger", title="Platelet Count",
            solidHeader=TRUE,
            plotlyOutput("plot_plt", height="300px")
          ),
          box(width=6, status="danger", title="Hemoglobin",
            solidHeader=TRUE,
            plotlyOutput("plot_hgb", height="300px")
          )
        ),
        fluidRow(
          box(width=4, status="warning", title="LDH (Hemolysis Marker)",
            solidHeader=TRUE,
            plotlyOutput("plot_ldh", height="280px")
          ),
          box(width=4, status="warning", title="Haptoglobin",
            solidHeader=TRUE,
            plotlyOutput("plot_hpg", height="280px")
          ),
          box(width=4, status="warning", title="Schistocytes (%)",
            solidHeader=TRUE,
            plotlyOutput("plot_schist", height="280px")
          )
        )
      ),

      # ── Tab 6: Renal Outcomes ────────────────────────────────────────────────
      tabItem("renal",
        fluidRow(
          box(width=6, status="primary", title="eGFR Trajectory",
            solidHeader=TRUE,
            plotlyOutput("plot_gfr", height="350px")
          ),
          box(width=6, status="danger", title="Endothelial Injury Score (TMA Driver)",
            solidHeader=TRUE,
            plotlyOutput("plot_ei", height="350px")
          )
        ),
        fluidRow(
          box(width=6, status="warning", title="Dialysis Risk Profile",
            solidHeader=TRUE,
            plotlyOutput("plot_dialysis", height="300px")
          ),
          box(width=6, status="info", title="GFR CKD Stage Classification",
            solidHeader=TRUE,
            plotlyOutput("plot_gfr_stage", height="300px")
          )
        )
      ),

      # ── Tab 7: Scenario Comparison ───────────────────────────────────────────
      tabItem("scenarios",
        fluidRow(
          box(width=12, status="primary",
            title="All Treatment Scenarios - Side-by-Side Comparison",
            solidHeader=TRUE,
            fluidRow(
              column(6, plotlyOutput("plot_scen_plt", height="300px")),
              column(6, plotlyOutput("plot_scen_gfr", height="300px"))
            ),
            fluidRow(
              column(6, plotlyOutput("plot_scen_hgb", height="300px")),
              column(6, plotlyOutput("plot_scen_ldh", height="300px"))
            )
          )
        ),
        fluidRow(
          box(width=12, status="info", title="Efficacy Summary Table at Day 90 & 180",
            solidHeader=TRUE,
            DTOutput("table_summary")
          )
        )
      ),

      # ── Tab 8: Biomarkers ────────────────────────────────────────────────────
      tabItem("biomarkers",
        fluidRow(
          box(width=6, status="warning", title="sC5b-9 (Soluble MAC Biomarker)",
            solidHeader=TRUE,
            plotlyOutput("plot_sc5b9", height="300px")
          ),
          box(width=6, status="info", title="CH50 (% Complement Activity)",
            solidHeader=TRUE,
            plotlyOutput("plot_ch50", height="300px")
          )
        ),
        fluidRow(
          box(width=6, status="danger", title="CRP (Inflammatory Marker)",
            solidHeader=TRUE,
            plotlyOutput("plot_crp", height="300px")
          ),
          box(width=6, status="success", title="Biomarker Correlation Dashboard",
            solidHeader=TRUE,
            plotlyOutput("plot_bm_corr", height="300px")
          )
        )
      ),

      # ── Tab 9: References ────────────────────────────────────────────────────
      tabItem("refs",
        box(width=12, status="primary", title="Key References",
          solidHeader=TRUE,
          h4("Pivotal Clinical Trials"),
          tags$ul(
            tags$li("Legendre CM et al. (2013). Terminal complement inhibitor eculizumab in atypical hemolytic-uremic syndrome. NEJM, 368(23):2169-2181. [PMID: 23738544]"),
            tags$li("Fakhouri F et al. (2017). Atypical hemolytic uremic syndrome. Lancet, 390(10105):1847-1860. [PMID: 28499565]"),
            tags$li("Cavero T et al. (2017). Eculizumab in adults with renal transplant and aHUS. JASN, 28(6):1913-1920."),
            tags$li("Licht C et al. (2015). Efficacy and safety of eculizumab in atypical HUS from 2-year extensions of phase 2 studies. Kidney Int, 87(5):1061-1073.")
          ),
          h4("PK/PD Modeling"),
          tags$ul(
            tags$li("Menne J et al. (2015). Validation of treatment strategies for atypical hemolytic uremic syndrome. JASN. Population PK model of eculizumab."),
            tags$li("Thomas TC et al. (2012). Inhibition of complement activity by humanized anti-C5 antibody and single chain Fv. Mol Immunol, 33:1389-1401."),
            tags$li("Rother RP et al. (2007). Discovery and development of the complement inhibitor eculizumab for the treatment of PNH. Nat Biotech, 25:1256-1264.")
          ),
          h4("Disease Mechanisms"),
          tags$ul(
            tags$li("Jokiranta TS (2017). HUS and atypical HUS. Blood, 129(21):2847-2856."),
            tags$li("Nester CM & Barbour T (2015). Atypical hemolytic uremic syndrome: a tale of two proteins. Clin J Am Soc Nephrol, 10(9):1610-1622."),
            tags$li("Kavanagh D et al. (2021). Atypical hemolytic uremic syndrome: a practical guide. Clin Kidney J.")
          )
        )
      )
    )
  )
)

# ─── SERVER ───────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # Reactive simulation
  sim_data <- eventReactive(input$run, {
    withProgress(message="Running simulation...", {
      params <- list(
        dose_ecu    = input$dose_ecu,
        dose_maint  = input$dose_maint,
        CFH_factor  = input$CFH_factor,
        CFI_factor  = input$CFI_factor,
        CD46_factor = input$CD46_factor,
        kAP         = input$kAP,
        kEI_MAC     = input$kEI_MAC,
        kEI_rep     = input$kEI_rep,
        kGFR_rep    = input$kGFR_rep
      )
      run_sim(as.numeric(input$scenario), params, input$sim_days)
    })
  }, ignoreNULL = FALSE)

  # Multi-scenario reactive
  all_scenarios <- eventReactive(input$run, {
    withProgress(message="Running all scenarios...", {
      params <- list(
        dose_ecu=900, dose_maint=1200,
        CFH_factor=input$CFH_factor, CFI_factor=input$CFI_factor,
        CD46_factor=input$CD46_factor, kAP=input$kAP,
        kEI_MAC=input$kEI_MAC, kEI_rep=input$kEI_rep, kGFR_rep=input$kGFR_rep
      )
      labs <- c("Natural History","Eculizumab Std","Ravulizumab",
                "Ecu+Danicopan","Iptacopan Oral")
      bind_rows(lapply(1:5, function(s) {
        d <- run_sim(s, params, input$sim_days)
        d$scenario_label <- labs[s]
        d
      }))
    })
  }, ignoreNULL = FALSE)

  scen_colors <- c(
    "Natural History"  = "#C0392B",
    "Eculizumab Std"   = "#2980B9",
    "Ravulizumab"      = "#27AE60",
    "Ecu+Danicopan"    = "#8E44AD",
    "Iptacopan Oral"   = "#E67E22"
  )

  # Value boxes (Day 90)
  last_val <- function(col) {
    d <- sim_data()
    d[[col]][nrow(d)]
  }

  output$box_PLT <- renderValueBox({
    v <- round(last_val("PLT"))
    valueBox(value=paste0(v, " x10⁹/L"), subtitle="Final Platelets",
             icon=icon("tint"), color=if(v<150) "red" else "green")
  })
  output$box_GFR <- renderValueBox({
    v <- round(last_val("GFR"))
    valueBox(value=paste0(v, " mL/min"), subtitle="Final eGFR",
             icon=icon("filter"), color=if(v<30) "red" else if(v<60) "yellow" else "green")
  })
  output$box_LDH <- renderValueBox({
    v <- round(last_val("LDH"))
    valueBox(value=paste0(v, " U/L"), subtitle="Final LDH",
             icon=icon("vials"), color=if(v>400) "red" else "green")
  })
  output$box_CH50 <- renderValueBox({
    v <- round(last_val("CH50pct"))
    valueBox(value=paste0(v, "%"), subtitle="CH50 Activity",
             icon=icon("shield-alt"), color=if(v>20) "red" else "green")
  })

  # Regulatory score
  output$reg_score <- renderText({
    reg <- round(input$CFH_factor * input$CFI_factor * input$CD46_factor, 3)
    paste0("Combined Regulatory Effectiveness: ", reg,
           "\n(Normal = 1.0, Severe disease = <0.10)")
  })

  # PK plots
  output$plot_pk_conc <- renderPlotly({
    d <- sim_data()
    p <- plot_ly(d, x=~time, y=~Cc_ugmL, type="scatter", mode="lines",
                 line=list(color="#2980B9", width=2),
                 name="Eculizumab (μg/mL)")
    p <- layout(p, xaxis=list(title="Day"),
                yaxis=list(title="Concentration (μg/mL)"),
                title="Eculizumab PK - Central Compartment",
                shapes=list(
                  list(type="line", x0=0, x1=max(d$time), y0=35, y1=35,
                       line=list(dash="dash", color="gray"))
                ))
    p
  })

  output$plot_c5_tmdd <- renderPlotly({
    d <- sim_data()
    plot_ly() %>%
      add_trace(data=d, x=~time, y=~C5free, type="scatter", mode="lines",
                name="Free C5 (nM)", line=list(color="#E74C3C")) %>%
      add_trace(data=d, x=~time, y=~Drug_C5/148, type="scatter", mode="lines",
                name="Drug-C5 Complex (proxy)", line=list(color="#27AE60")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Concentration (nM)"),
             title="TMDD: Free vs Bound C5")
  })

  output$plot_c5block <- renderPlotly({
    d <- sim_data()
    d$C5_block <- 1 - d$C5free / 395
    plot_ly(d, x=~time, y=~I(C5_block*100), type="scatter", mode="lines",
            line=list(color="#8E44AD", width=2)) %>%
      layout(xaxis=list(title="Day"),
             yaxis=list(title="C5 Blockade (%)", range=c(0,100)),
             title="Fractional C5 Inhibition")
  })

  # Complement plots
  output$plot_c3 <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~C3pool, type="scatter", mode="lines",
            line=list(color="#3498DB", width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="C3 (nM)"),
             title="C3 Pool Dynamics")
  })

  output$plot_mac <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~MACflux, type="scatter", mode="lines",
            line=list(color="#C0392B", width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="MAC (AU)"),
             title="MAC Flux (Disease Driver)")
  })

  output$plot_c3b <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~C3b_AP, type="scatter", mode="lines",
            line=list(color="#E67E22", width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="C3b AP (nM)"),
             title="Alternative Pathway C3b")
  })

  output$plot_c5conv <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~C5conv, type="scatter", mode="lines",
            line=list(color="#8E44AD", width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="C5 Convertase (AU)"),
             title="C5 Convertase Activity")
  })

  # Hematology
  output$plot_plt <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~PLT, type="scatter", mode="lines",
            line=list(color="#E74C3C", width=2)) %>%
      add_trace(x=c(0, max(d$time)), y=c(150,150), type="scatter", mode="lines",
                line=list(dash="dash", color="gray"), name="Threshold") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Platelets (x10⁹/L)"),
             title="Platelet Count")
  })

  output$plot_hgb <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~Hgb, type="scatter", mode="lines",
            line=list(color="#8E44AD", width=2)) %>%
      add_trace(x=c(0,max(d$time)), y=c(10,10), type="scatter", mode="lines",
                line=list(dash="dash", color="orange"), name="Target Hgb") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Hemoglobin (g/dL)"),
             title="Hemoglobin")
  })

  output$plot_ldh <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~LDH, type="scatter", mode="lines",
            line=list(color="#F39C12", width=2)) %>%
      add_trace(x=c(0,max(d$time)), y=c(200,200), type="scatter", mode="lines",
                line=list(dash="dash", color="gray"), name="ULN") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="LDH (U/L)"),
             title="LDH")
  })

  output$plot_hpg <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~Hpg, type="scatter", mode="lines",
            line=list(color="#1ABC9C", width=2)) %>%
      add_trace(x=c(0,max(d$time)), y=c(0.1,0.1), type="scatter", mode="lines",
                line=list(dash="dash", color="gray"), name="Low limit") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Haptoglobin (g/L)"),
             title="Haptoglobin")
  })

  output$plot_schist <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~Schist, type="scatter", mode="lines",
            line=list(color="#C0392B", width=2)) %>%
      add_trace(x=c(0,max(d$time)), y=c(1,1), type="scatter", mode="lines",
                line=list(dash="dash", color="orange"), name=">1% = TMA") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Schistocytes (%)"),
             title="Schistocyte Count")
  })

  # Renal
  output$plot_gfr <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~GFR, type="scatter", mode="lines",
            line=list(color="#2980B9", width=2)) %>%
      add_trace(x=c(0,max(d$time)), y=c(15,15), type="scatter", mode="lines",
                line=list(dash="dash", color="red"), name="ESRD threshold") %>%
      add_trace(x=c(0,max(d$time)), y=c(60,60), type="scatter", mode="lines",
                line=list(dash="dash", color="orange"), name="CKD G3") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="eGFR (mL/min/1.73m²)"),
             title="eGFR Trajectory")
  })

  output$plot_ei <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~Endo_inj, type="scatter", mode="lines",
            line=list(color="#E74C3C", width=2), fill="tozeroy",
            fillcolor="rgba(231,76,60,0.2)") %>%
      layout(xaxis=list(title="Day"),
             yaxis=list(title="Endothelial Injury Score (0-10)"),
             title="Endothelial Injury (TMA Driver)")
  })

  output$plot_dialysis <- renderPlotly({
    d <- sim_data()
    d$dialysis_risk <- as.integer(d$GFR < 15)
    plot_ly(d, x=~time, y=~GFR, type="scatter", mode="lines",
            name="eGFR", line=list(color="#2980B9")) %>%
      add_trace(x=~time, y=~I(dialysis_risk*15), type="scatter", mode="lines",
                fill="tozeroy", fillcolor="rgba(231,76,60,0.3)",
                line=list(color="transparent"), name="Dialysis Zone") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="eGFR"),
             title="Dialysis Risk Zone (GFR <15)")
  })

  output$plot_gfr_stage <- renderPlotly({
    d <- sim_data()
    d$ckd_stage <- cut(d$GFR,
      breaks=c(-Inf, 15, 30, 45, 60, 90, Inf),
      labels=c("G5 (ESRD)","G4","G3b","G3a","G2","G1"))
    stage_colors <- c("G5 (ESRD)"="#C0392B","G4"="#E74C3C","G3b"="#E67E22",
                      "G3a"="#F39C12","G2"="#F9E79F","G1"="#27AE60")
    plot_ly(d, x=~time, y=~GFR, type="scatter", mode="lines",
            line=list(color="#2980B9")) %>%
      layout(xaxis=list(title="Day"),
             yaxis=list(title="eGFR (mL/min/1.73m²)"),
             title="GFR with CKD Stage Boundaries",
             shapes=lapply(c(15,30,45,60,90), function(g) {
               list(type="line", x0=0, x1=max(d$time), y0=g, y1=g,
                    line=list(dash="dash", color="gray", width=1))
             }))
  })

  # Scenario comparison
  make_scen_plot <- function(var, title, ylab, threshold=NULL) {
    d <- all_scenarios()
    p <- plot_ly()
    for(s in unique(d$scenario_label)) {
      dd <- d[d$scenario_label==s,]
      p <- add_trace(p, data=dd, x=~time, y=as.formula(paste0("~",var)),
                     type="scatter", mode="lines",
                     name=s, line=list(color=scen_colors[s], width=2))
    }
    if(!is.null(threshold)) {
      p <- add_trace(p, x=c(0,max(d$time)), y=c(threshold,threshold),
                     type="scatter", mode="lines",
                     line=list(dash="dash", color="gray"), name="Threshold",
                     showlegend=FALSE)
    }
    layout(p, xaxis=list(title="Day"), yaxis=list(title=ylab), title=title,
           legend=list(x=0.01, y=0.01))
  }

  output$plot_scen_plt  <- renderPlotly(make_scen_plot("PLT","Platelet Count","x10⁹/L",150))
  output$plot_scen_gfr  <- renderPlotly(make_scen_plot("GFR","eGFR","mL/min/1.73m²",60))
  output$plot_scen_hgb  <- renderPlotly(make_scen_plot("Hgb","Hemoglobin","g/dL",10))
  output$plot_scen_ldh  <- renderPlotly(make_scen_plot("LDH","LDH","U/L",200))

  output$table_summary <- renderDT({
    d <- all_scenarios()
    d %>%
      filter(time %in% c(28,90,180,365)) %>%
      group_by(scenario_label, time) %>%
      summarise(
        PLT_mean=round(mean(PLT),0), Hgb_mean=round(mean(Hgb),1),
        GFR_mean=round(mean(GFR),0), LDH_xULN=round(mean(LDH/200),1),
        sC5b9=round(mean(sC5b9),0), CH50pct=round(mean(CH50pct),0),
        .groups="drop"
      ) %>%
      rename(
        "Scenario"=scenario_label, "Day"=time,
        "PLT (x10⁹/L)"=PLT_mean, "Hgb (g/dL)"=Hgb_mean,
        "GFR (mL/min)"=GFR_mean, "LDH (×ULN)"=LDH_xULN,
        "sC5b-9 (ng/mL)"=sC5b9, "CH50 (%)"=CH50pct
      ) %>%
      datatable(rownames=FALSE,
                options=list(pageLength=20, scrollX=TRUE),
                class="table-bordered table-sm") %>%
      formatStyle("Scenario", fontWeight="bold")
  })

  # Biomarkers
  output$plot_sc5b9 <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~sC5b9, type="scatter", mode="lines",
            line=list(color="#E67E22", width=2)) %>%
      add_trace(x=c(0,max(d$time)), y=c(244,244), type="scatter", mode="lines",
                line=list(dash="dash", color="gray"), name="Normal ULN") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="sC5b-9 (ng/mL)"),
             title="sC5b-9 (Soluble MAC Biomarker)")
  })

  output$plot_ch50 <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~CH50pct, type="scatter", mode="lines",
            line=list(color="#3498DB", width=2)) %>%
      add_trace(x=c(0,max(d$time)), y=c(10,10), type="scatter", mode="lines",
                line=list(dash="dash", color="green"), name="Target <10%") %>%
      layout(xaxis=list(title="Day"),
             yaxis=list(title="CH50 (% activity)", range=c(0,100)),
             title="CH50 Functional Complement Activity")
  })

  output$plot_crp <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~CRP, type="scatter", mode="lines",
            line=list(color="#C0392B", width=2)) %>%
      add_trace(x=c(0,max(d$time)), y=c(10,10), type="scatter", mode="lines",
                line=list(dash="dash", color="gray"), name="ULN (10 mg/L)") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="CRP (mg/L)"),
             title="CRP (Inflammatory Marker)")
  })

  output$plot_bm_corr <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~sC5b9, y=~PLT, type="scatter", mode="markers",
            marker=list(color=~time, colorscale="Viridis", showscale=TRUE,
                        colorbar=list(title="Day"), size=4),
            text=~paste0("Day: ", round(time), "<br>sC5b-9: ", round(sC5b9),
                         "<br>PLT: ", round(PLT))) %>%
      layout(xaxis=list(title="sC5b-9 (ng/mL)"),
             yaxis=list(title="Platelets (x10⁹/L)"),
             title="sC5b-9 vs PLT Correlation (color=time)")
  })

  output$plot_presentation <- renderPlotly({
    d <- sim_data()
    d0 <- d[1,]
    vars <- c("PLT","Hgb","GFR","sC5b9","CH50pct","Schist")
    normal <- c(250, 14, 90, 150, 100, 0.5)
    vals   <- c(d0$PLT, d0$Hgb, d0$GFR, d0$sC5b9, d0$CH50pct, d0$Schist)
    pct_norm <- pmin(vals/normal*100, 300)
    plot_ly(x=vars, y=pct_norm, type="bar",
            marker=list(color=ifelse(pct_norm>100, "#27AE60","#E74C3C"))) %>%
      layout(xaxis=list(title=""), yaxis=list(title="% of Normal"),
             title="Presentation Values (% of normal)")
  })

  output$plot_cfh_radar <- renderPlotly({
    theta <- c("CFH","CFI","CD46","CFH","CFI")
    r_patient <- c(input$CFH_factor, input$CFI_factor, input$CD46_factor,
                   input$CFH_factor, input$CFI_factor)
    r_normal  <- c(1.0, 1.0, 1.0, 1.0, 1.0)
    plot_ly(type="scatterpolar", fill="toself") %>%
      add_trace(r=r_normal, theta=c("CFH","CFI","CD46","CFH"),
                name="Normal", fillcolor="rgba(39,174,96,0.2)",
                line=list(color="#27AE60")) %>%
      add_trace(r=r_patient[1:4], theta=c("CFH","CFI","CD46","CFH"),
                name="Patient", fillcolor="rgba(231,76,60,0.2)",
                line=list(color="#E74C3C")) %>%
      layout(polar=list(radialaxis=list(range=c(0,1))),
             title="Complement Regulatory Function")
  })

  output$plot_complement_cascade <- renderPlotly({
    d <- sim_data()
    plot_ly() %>%
      add_trace(data=d, x=~time, y=~C3b_AP, name="C3b AP (nM/10)",
                y=~I(C3b_AP/10), type="scatter", mode="lines",
                line=list(color="#3498DB")) %>%
      add_trace(data=d, x=~time, y=~C5conv, name="C5 Convertase (AU)",
                type="scatter", mode="lines", line=list(color="#E67E22")) %>%
      add_trace(data=d, x=~time, y=~MACflux, name="MAC Flux (AU)",
                type="scatter", mode="lines", line=list(color="#C0392B")) %>%
      layout(xaxis=list(title="Day"),
             yaxis=list(title="Activity Level"),
             title="Complement Cascade Dynamics")
  })
}

# ─── Launch ────────────────────────────────────────────────────────────────────

shinyApp(ui, server)
