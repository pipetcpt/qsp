##############################################################################
#  Uveitis QSP Shiny Dashboard
#  Interactive simulation of uveitis inflammation, barrier disruption,
#  and treatment response across 7 drug scenarios
#
#  Tabs:
#   1. Patient Profile & Disease Setup
#   2. Pharmacokinetics (Plasma & Ocular Drug Levels)
#   3. Immune & Cytokine Dynamics
#   4. Blood-Ocular Barrier Integrity
#   5. Clinical Endpoints (VA, CME, IOP, SUN Grade)
#   6. Scenario Comparison (all 7 treatments)
#   7. Virtual Patient Population (VP simulation)
#   8. Biomarker Monitoring
##############################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)

# ─────────────────────────────────────────────────────────────────────────────
# Inline model code (same ODE structure as main model file)
# ─────────────────────────────────────────────────────────────────────────────
uvt_code <- '
$PROB Uveitis QSP — Shiny App Model
$PARAM @annotated
CL:0.45:L/h Clearance  Vd:12.0:L Volume  k12:0.10:1/h  k21:0.05:1/h
ka:0.80:1/h  F_oral:0.80:—  F_sc:0.64:—
kBAB:0.003:1/h  kBRB:0.001:1/h  CL_eye:0.04:1/h
k_depot:0.005:1/h
k_Tprol:0.02:1/h  k_Tdeath:0.008:1/h
k_Treg0:0.015:1/h  k_Tregd:0.006:1/h
T_eff0:1.0:AU  T_reg0:2.0:AU
k_APCact:0.03:1/h  k_APCdec:0.04:1/h  APC0:0.5:AU
k_Mac:0.02:1/h  k_Macd:0.015:1/h  Mac0:0.8:AU
ksyn_TNF:0.05:AU/h  kdeg_TNF:0.15:1/h
ksyn_IL6:0.04:AU/h  kdeg_IL6:0.10:1/h
ksyn_IL17:0.03:AU/h  kdeg_IL17:0.12:1/h
ksyn_VEGF:0.02:AU/h  kdeg_VEGF:0.08:1/h
TNF0:0.2:AU  IL60:0.15:AU  IL170:0.1:AU  VEGF0:0.08:AU
BAB_base:1.0:AU  k_BABdeg:0.04:1/h  k_BABrep:0.006:1/h
BRB_base:1.0:AU  k_BRBdeg:0.035:1/h  k_BRBrep:0.005:1/h
k_CellsAH:0.08:1/h  k_CellsCl:0.03:1/h
k_CMEform:0.005:1/h  k_CMEres:0.002:1/h
k_VAdeg:0.02:1/h  k_VArec:0.005:1/h
IOP_base:15.0:mmHg  k_IOPinf:0.5:mmHg/AU
Emax_cs:0.85:—  EC50_cs:5.0:AU  Hill_cs:1.5:—
Emax_anti_tnf:0.90:—  EC50_aTNF:2.0:AU  Hill_aTNF:2.0:—
Emax_aVEGF:0.80:—  EC50_aVEGF:1.0:AU
Ag_stim:0.5:AU  flare_on:1:0/1

$CMT Cgut Cp Cperiph C_ant C_post C_depot
     T_eff T_reg APC_act Macro TNF IL6 IL17 VEGF
     BAB_int BRB_int Cells_AH CME VA_def IOP_e GR_occ

$MAIN
double flare_mult = (flare_on>0.5)?3.0:1.0;
Cgut_0=0; Cp_0=0; Cperiph_0=0; C_ant_0=0; C_post_0=0; C_depot_0=0;
T_eff_0=T_eff0*flare_mult; T_reg_0=T_reg0/flare_mult;
APC_act_0=APC0*flare_mult; Macro_0=Mac0*flare_mult;
TNF_0=TNF0*flare_mult; IL6_0=IL60*flare_mult;
IL17_0=IL170*flare_mult; VEGF_0=VEGF0*flare_mult;
BAB_int_0=BAB_base*(1.0/flare_mult); BRB_int_0=BRB_base*(1.0/flare_mult);
Cells_AH_0=2.0*flare_mult; CME_0=0.3*(flare_mult-1.0);
VA_def_0=0.05*(flare_mult-1.0); IOP_e_0=0.0; GR_occ_0=0.0;

$ODE
dxdt_Cgut=-ka*Cgut;
dxdt_Cp=ka*F_oral*Cgut-(CL/Vd)*Cp-k12*Cp+k21*Cperiph-kBAB*Cp-kBRB*Cp;
dxdt_Cperiph=k12*Cp-k21*Cperiph;
dxdt_C_ant=kBAB*BAB_int*Cp+0.3*C_depot-CL_eye*C_ant;
dxdt_C_post=kBRB*BRB_int*Cp+k_depot*C_depot-CL_eye*C_post;
dxdt_C_depot=-k_depot*C_depot;
double drug_conc=C_ant+0.3*Cp;
double hill_cs=pow(drug_conc,Hill_cs);
double hill_cs50=pow(EC50_cs,Hill_cs);
double GR_effect=Emax_cs*hill_cs/(hill_cs50+hill_cs);
dxdt_GR_occ=GR_effect-GR_occ;
double hill_aTNF=pow(Cp,Hill_aTNF);
double hill_aTNF50=pow(EC50_aTNF,Hill_aTNF);
double aTNF_eff=Emax_anti_tnf*hill_aTNF/(hill_aTNF50+hill_aTNF);
double aVEGF_eff=Emax_aVEGF*C_post/(EC50_aVEGF+C_post);
double Teff_prolif=k_Tprol*T_eff*TNF*(1.0-GR_occ);
double Teff_reg=0.15*T_reg*T_eff;
dxdt_T_eff=Ag_stim+Teff_prolif-k_Tdeath*T_eff-Teff_reg-0.20*GR_occ*T_eff;
dxdt_T_reg=k_Treg0*(1.0+0.5*GR_occ)-k_Tregd*T_reg-0.1*TNF*T_reg;
dxdt_APC_act=k_APCact*Ag_stim-k_APCdec*APC_act-0.30*GR_occ*APC_act;
dxdt_Macro=k_Mac*T_eff*APC_act-k_Macd*Macro-0.25*GR_occ*Macro;
double TNF_prod=ksyn_TNF*(Macro+T_eff);
double TNF_neut=aTNF_eff*TNF;
dxdt_TNF=TNF_prod*(1.0-GR_occ*0.70)-kdeg_TNF*TNF-TNF_neut;
dxdt_IL6=ksyn_IL6*(Macro+0.5*TNF)*(1.0-GR_occ*0.65)-kdeg_IL6*IL6;
dxdt_IL17=ksyn_IL17*T_eff*(1.0-GR_occ*0.50)-kdeg_IL17*IL17;
dxdt_VEGF=ksyn_VEGF*(TNF+0.5*IL6)*(1.0-GR_occ*0.40)-kdeg_VEGF*VEGF-aVEGF_eff*VEGF;
double BAB_disrupt=k_BABdeg*(TNF+0.5*IL6)*BAB_int;
double BAB_repair=k_BABrep*(1.0+2.0*GR_occ)*(BAB_base-BAB_int);
dxdt_BAB_int=BAB_repair-BAB_disrupt;
double BRB_disrupt=k_BRBdeg*(VEGF+0.3*TNF)*BRB_int;
double BRB_repair=k_BRBrep*(1.0+1.5*GR_occ+2.0*aVEGF_eff)*(BRB_base-BRB_int);
dxdt_BRB_int=BRB_repair-BRB_disrupt;
double BAB_breach=fmax(0.0,BAB_base-BAB_int);
dxdt_Cells_AH=k_CellsAH*BAB_breach*T_eff-k_CellsCl*Cells_AH;
double BRB_breach=fmax(0.0,BRB_base-BRB_int);
dxdt_CME=k_CMEform*BRB_breach*VEGF-k_CMEres*(1.0+2.0*aVEGF_eff+GR_occ)*CME;
dxdt_VA_def=k_VAdeg*CME-k_VArec*(1.0+aVEGF_eff)*VA_def;
double IOP_inflam=k_IOPinf*(1.0-BAB_int)*Cells_AH;
double IOP_steroid=2.5*GR_occ*GR_occ;
dxdt_IOP_e=IOP_inflam+IOP_steroid-0.1*IOP_e;

$TABLE
capture BCVA_logMAR=VA_def;
capture IOP_mmHg=IOP_base+IOP_e;
capture CST_um=250+200*CME;
capture SUN_grade=fmin(4.0,Cells_AH/2.0);
capture TNF_level=TNF; capture IL6_level=IL6; capture VEGF_level=VEGF;
capture BAB=BAB_int; capture BRB=BRB_int;
capture Teff_cells=T_eff; capture Treg_cells=T_reg;
capture GR_occupancy=GR_occ;
'

mod <- suppressMessages(mcode("uvt_shiny", uvt_code))

# ─────────────────────────────────────────────────────────────────────────────
# Build dosing events
# ─────────────────────────────────────────────────────────────────────────────
build_dosing <- function(scen, dose_mg, freq_h, n_doses, anti_tnf_sc=FALSE) {
  if (anti_tnf_sc) {
    ev(amt = dose_mg, cmt = "Cp", rate = 0, ii = freq_h, addl = n_doses - 1, time = 0)
  } else {
    ev(amt = dose_mg, cmt = "Cgut", rate = 0, ii = freq_h, addl = n_doses - 1, time = 0)
  }
}

run_sim <- function(params, dosing_ev) {
  mrgsim(param_set(mod, params),
         ev = dosing_ev,
         end = 365 * 24, delta = 6) %>%
    as_tibble() %>%
    mutate(time_d = time / 24)
}

scenarios_info <- list(
  S1 = list(label="No Treatment", anti_tnf=FALSE, amt=0, ii=24, n=0),
  S2 = list(label="Topical Prednisolone 1% QID", anti_tnf=FALSE, amt=0.5, ii=6, n=60),
  S3 = list(label="Periocular Triamcinolone 40mg", anti_tnf=FALSE, amt=40, ii=999, n=1),
  S4 = list(label="IVT Dexamethasone Implant", anti_tnf=FALSE, depot=TRUE, amt=700, ii=999, n=1),
  S5 = list(label="Systemic Prednisone 1mg/kg/d", anti_tnf=FALSE, amt=60, ii=24, n=90),
  S6 = list(label="Adalimumab 40mg q2w", anti_tnf=TRUE, amt=40, ii=336, n=12),
  S7 = list(label="Combination (Pred+Adalimumab)", anti_tnf=FALSE, combo=TRUE, amt=60, ii=24, n=90)
)

# ─────────────────────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Uveitis QSP Dashboard", titleWidth = 320),
  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile", tabName="profile", icon=icon("user-md")),
      menuItem("Pharmacokinetics", tabName="pk", icon=icon("flask")),
      menuItem("Immune & Cytokines", tabName="immune", icon=icon("dna")),
      menuItem("Barrier Integrity", tabName="barrier", icon=icon("shield-alt")),
      menuItem("Clinical Endpoints", tabName="clinical", icon=icon("eye")),
      menuItem("Scenario Comparison", tabName="scenario", icon=icon("chart-bar")),
      menuItem("Virtual Patients", tabName="vp", icon=icon("users")),
      menuItem("Biomarker Monitor", tabName="biomarker", icon=icon("microscope"))
    ),
    hr(),
    h5("  Global Simulation Settings", style="color:#ECF0F1;padding-left:10px"),
    sliderInput("sim_dur", "Duration (days):", min=30, max=730, value=365, step=30),
    sliderInput("flare_severity", "Flare Severity (1-5):", min=1, max=5, value=3, step=1),
    selectInput("uv_type", "Uveitis Type:",
                choices=c("Anterior"="ant","Intermediate"="int","Posterior"="post","Panuveitis"="pan"),
                selected="post"),
    checkboxInput("hlab27", "HLA-B27 Positive", value=FALSE)
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background: #ECF0F1; }
      .box { border-radius:8px; }
      .small-box { border-radius:8px; }
    "))),
    tabItems(
      # ── Tab 1: Patient Profile ──────────────────────────────────────────────
      tabItem("profile",
        fluidRow(
          box(title="Patient & Disease Parameters", width=4, status="primary", solidHeader=TRUE,
            sliderInput("pt_age", "Age (years):", 18, 80, 42),
            selectInput("pt_sex", "Sex:", c("Female"="F","Male"="M"), "F"),
            selectInput("uveitis_etiology", "Etiology:",
                        choices=c("Idiopathic"="idio","HLA-B27 Associated"="hlab27",
                                  "Birdshot (HLA-A29)"="bird","Sarcoidosis"="sarc",
                                  "VKH Syndrome"="vkh","Behcet Disease"="behcet",
                                  "MS-Associated"="ms","JIA-Associated"="jia"),
                        selected="idio"),
            sliderInput("baseline_va", "Baseline BCVA (logMAR):", 0, 1.0, 0.2, step=0.05),
            sliderInput("baseline_iop", "Baseline IOP (mmHg):", 8, 25, 15),
            checkboxInput("chronic_uveitis", "Chronic/Recurrent Uveitis", TRUE),
            numericInput("flare_rate", "Annual Flare Rate:", 2, min=0, max=8, step=0.5)
          ),
          box(title="Treatment Selection", width=4, status="warning", solidHeader=TRUE,
            selectInput("treatment", "Scenario:",
                        choices=c("S1 — No Treatment"="S1",
                                  "S2 — Topical Prednisolone QID"="S2",
                                  "S3 — Periocular Triamcinolone 40mg"="S3",
                                  "S4 — IVT Dexamethasone Implant"="S4",
                                  "S5 — Systemic Prednisone 1mg/kg"="S5",
                                  "S6 — Adalimumab 40mg q2w"="S6",
                                  "S7 — Combination (Pred+Adalimumab)"="S7"),
                        selected="S6"),
            numericInput("dose_mg", "Dose (mg):", 40, min=1, max=200),
            numericInput("dose_freq_h", "Dosing Interval (hours):", 336, min=6, max=720),
            numericInput("n_doses", "Number of Doses:", 12, min=1, max=365),
            checkboxInput("anti_tnf_route", "Biologic (SC/IV route)", TRUE),
            actionButton("run_sim_btn", "Run Simulation", icon=icon("play"),
                         class="btn-success btn-block", style="margin-top:15px")
          ),
          box(title="Disease Severity Score", width=4, status="danger", solidHeader=TRUE,
            h4("SUN Grading Scale"),
            tableOutput("sun_scale"),
            hr(),
            h4("Activity Indices"),
            tableOutput("activity_score_tbl")
          )
        ),
        fluidRow(
          valueBoxOutput("vbox_va", width=3),
          valueBoxOutput("vbox_cst", width=3),
          valueBoxOutput("vbox_sun", width=3),
          valueBoxOutput("vbox_iop", width=3)
        )
      ),
      # ── Tab 2: Pharmacokinetics ─────────────────────────────────────────────
      tabItem("pk",
        fluidRow(
          box(title="Plasma Drug Concentration", width=6, status="primary", solidHeader=TRUE,
              plotOutput("pk_plasma", height=320)),
          box(title="Ocular Drug Concentration", width=6, status="info", solidHeader=TRUE,
              plotOutput("pk_ocular", height=320))
        ),
        fluidRow(
          box(title="GR Occupancy & PD Effect", width=6, status="warning", solidHeader=TRUE,
              plotOutput("pk_gr", height=280)),
          box(title="PK Parameters Summary", width=6, status="success", solidHeader=TRUE,
              tableOutput("pk_summary_tbl"))
        )
      ),
      # ── Tab 3: Immune & Cytokines ───────────────────────────────────────────
      tabItem("immune",
        fluidRow(
          box(title="T Effector / Regulatory T Cells", width=6, status="danger", solidHeader=TRUE,
              plotOutput("immune_tcell", height=300)),
          box(title="APC & Macrophage Activation", width=6, status="warning", solidHeader=TRUE,
              plotOutput("immune_apc", height=300))
        ),
        fluidRow(
          box(title="TNF-α & IL-6 Dynamics", width=6, status="primary", solidHeader=TRUE,
              plotOutput("cyto_tnf_il6", height=300)),
          box(title="IL-17A & VEGF Dynamics", width=6, status="info", solidHeader=TRUE,
              plotOutput("cyto_il17_vegf", height=300))
        )
      ),
      # ── Tab 4: Barrier Integrity ────────────────────────────────────────────
      tabItem("barrier",
        fluidRow(
          box(title="Blood-Aqueous Barrier (BAB) Integrity", width=6, status="primary", solidHeader=TRUE,
              plotOutput("bab_plot", height=300)),
          box(title="Blood-Retinal Barrier (BRB) Integrity", width=6, status="info", solidHeader=TRUE,
              plotOutput("brb_plot", height=300))
        ),
        fluidRow(
          box(title="Anterior Chamber Cells (SUN Grade)", width=6, status="danger", solidHeader=TRUE,
              plotOutput("cells_ah_plot", height=280)),
          box(title="BAB/BRB — Mechanistic Explanation", width=6, status="success", solidHeader=TRUE,
            HTML("<h5>Blood-Aqueous Barrier (BAB)</h5>
                 <p>Maintains immune privilege in the anterior segment. Disrupted by TNF-α, IL-1β,
                 and inflammatory cells → protein flare + WBC influx.</p>
                 <h5>Blood-Retinal Barrier (BRB)</h5>
                 <p>Inner BRB (retinal capillary endothelium) + Outer BRB (RPE tight junctions).
                 VEGF-A is the primary driver of BRB disruption → macular edema (CME/CST).</p>
                 <h5>Treatment Effects</h5>
                 <ul>
                   <li><b>Corticosteroids</b>: Restore tight junctions via GR-ZO1 signaling</li>
                   <li><b>Anti-TNF</b>: Reduce leukostasis and endothelial activation</li>
                   <li><b>Anti-VEGF (IVT)</b>: Directly block VEGF-driven BRB disruption</li>
                 </ul>")
          )
        )
      ),
      # ── Tab 5: Clinical Endpoints ───────────────────────────────────────────
      tabItem("clinical",
        fluidRow(
          box(title="Visual Acuity (logMAR) — lower is better", width=6, status="danger", solidHeader=TRUE,
              plotOutput("va_plot", height=300)),
          box(title="OCT Central Subfield Thickness (CME)", width=6, status="warning", solidHeader=TRUE,
              plotOutput("cst_plot", height=300))
        ),
        fluidRow(
          box(title="Intraocular Pressure (IOP)", width=6, status="primary", solidHeader=TRUE,
              plotOutput("iop_plot", height=280)),
          box(title="Clinical Milestones Summary", width=6, status="success", solidHeader=TRUE,
              DTOutput("milestones_tbl"))
        )
      ),
      # ── Tab 6: Scenario Comparison ──────────────────────────────────────────
      tabItem("scenario",
        fluidRow(
          box(title="All Scenarios — Visual Acuity", width=6, status="danger", solidHeader=TRUE,
              plotOutput("scen_va", height=300)),
          box(title="All Scenarios — CME (CST)", width=6, status="warning", solidHeader=TRUE,
              plotOutput("scen_cst", height=300))
        ),
        fluidRow(
          box(title="All Scenarios — TNF-α (0-90d)", width=6, status="primary", solidHeader=TRUE,
              plotOutput("scen_tnf", height=280)),
          box(title="All Scenarios — Summary Table (Day 90)", width=6, status="info", solidHeader=TRUE,
              DTOutput("scen_tbl"))
        )
      ),
      # ── Tab 7: Virtual Patients ─────────────────────────────────────────────
      tabItem("vp",
        fluidRow(
          box(title="VP Population Parameters", width=3, status="primary", solidHeader=TRUE,
            sliderInput("vp_n", "N Virtual Patients:", 50, 500, 100, step=50),
            selectInput("vp_scen", "Scenario:",
                        choices=c("S6 — Adalimumab"="S6","S5 — Systemic Pred"="S5",
                                  "S7 — Combination"="S7","S1 — No Treat"="S1"),
                        selected="S6"),
            sliderInput("vp_iiv", "IIV (% CV):", 10, 60, 30),
            actionButton("run_vp_btn", "Run VP Simulation", icon=icon("play"),
                         class="btn-info btn-block")
          ),
          box(title="BCVA Response by Disease Subtype", width=9, status="info", solidHeader=TRUE,
              plotOutput("vp_bcva", height=380))
        ),
        fluidRow(
          box(title="Responder Analysis (SUN Grade < 0.5)", width=6, status="success", solidHeader=TRUE,
              plotOutput("vp_responder", height=280)),
          box(title="Variability in CME Response", width=6, status="warning", solidHeader=TRUE,
              plotOutput("vp_cme_var", height=280))
        )
      ),
      # ── Tab 8: Biomarker Monitor ────────────────────────────────────────────
      tabItem("biomarker",
        fluidRow(
          box(title="Cytokine Biomarker Trajectories", width=6, status="primary", solidHeader=TRUE,
              plotOutput("bm_cytokines", height=300)),
          box(title="OCT & Clinical Imaging Metrics", width=6, status="info", solidHeader=TRUE,
              plotOutput("bm_oct", height=300))
        ),
        fluidRow(
          box(title="Drug Level & Anti-Drug Antibody Risk", width=6, status="warning", solidHeader=TRUE,
              plotOutput("bm_drug_level", height=280)),
          box(title="Monitoring Schedule Table", width=6, status="success", solidHeader=TRUE,
              tableOutput("monitor_schedule"))
        )
      )
    )
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# Server
# ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  colors7 <- c("#E74C3C","#F39C12","#27AE60","#2980B9","#8E44AD","#1ABC9C","#C0392B")

  # Build parameters from UI inputs
  make_params <- function() {
    flare_mult <- input$flare_severity
    list(
      flare_on   = 1,
      Ag_stim    = 0.3 * flare_mult,
      BAB_base   = 1.0,
      BRB_base   = 1.0,
      k_BABdeg   = if (input$uv_type %in% c("ant","pan")) 0.06 else 0.03,
      k_BRBdeg   = if (input$uv_type %in% c("post","pan")) 0.05 else 0.02,
      ksyn_TNF   = if (input$hlab27) 0.08 else 0.05,
      ksyn_IL17  = if (input$hlab27) 0.06 else 0.03
    )
  }

  # Make dosing event from user inputs
  make_dose_ev <- function(scen) {
    if (scen == "S1") return(ev())
    if (scen == "S4") {
      return(ev(amt=700, cmt="C_depot", rate=0, time=0))
    }
    cmt_use <- if (input$anti_tnf_route && scen %in% c("S6","S7")) "Cp" else "Cgut"
    amt <- input$dose_mg
    ii  <- input$dose_freq_h
    n   <- input$n_doses
    if (scen == "S7") {
      e1 <- ev(amt=60, cmt="Cgut", rate=0, ii=24, addl=89, time=0)
      e2 <- ev(amt=40, cmt="Cp",   rate=0, ii=336, addl=11, time=0)
      return(e1 + e2)
    }
    ev(amt=amt, cmt=cmt_use, rate=0, ii=ii, addl=max(0,n-1), time=0)
  }

  # Single scenario simulation (reactive)
  sim_data <- eventReactive(input$run_sim_btn, {
    params  <- make_params()
    dosing  <- make_dose_ev(input$treatment)
    end_h   <- input$sim_dur * 24
    out <- mrgsim(param_set(mod, params), ev=dosing,
                  end=end_h, delta=6) %>%
      as_tibble() %>%
      mutate(time_d = time/24)
    out
  }, ignoreNULL=FALSE)

  # All scenarios simulation (for comparison tab)
  all_scen_data <- eventReactive(input$run_sim_btn, {
    params <- make_params()
    scenarios_ev <- list(
      S1 = ev(),
      S2 = ev(amt=0.5, cmt="Cgut", ii=6, addl=59, time=0),
      S3 = ev(amt=40, cmt="Cgut", time=0),
      S4 = ev(amt=700, cmt="C_depot", time=0),
      S5 = ev(amt=60, cmt="Cgut", ii=24, addl=89, time=0),
      S6 = ev(amt=40, cmt="Cp", ii=336, addl=11, time=0),
      S7 = ev(amt=60, cmt="Cgut", ii=24, addl=89, time=0) +
           ev(amt=40, cmt="Cp", ii=336, addl=11, time=0)
    )
    scen_labels <- c("S1 No Treat","S2 Topical CS","S3 Perioc Triam",
                     "S4 IVT Implant","S5 Systemic Pred","S6 Adalimumab","S7 Combination")
    bind_rows(mapply(function(scen_ev, lbl) {
      mrgsim(param_set(mod, params), ev=scen_ev,
             end=input$sim_dur*24, delta=12) %>%
        as_tibble() %>%
        mutate(time_d=time/24, scenario=lbl)
    }, scenarios_ev, scen_labels, SIMPLIFY=FALSE))
  }, ignoreNULL=FALSE)

  # ── Value Boxes ──────────────────────────────────────────────────────────
  get_last <- function(col) {
    d <- sim_data()
    if (nrow(d)==0) return(NA)
    tail(d[[col]], 1)
  }
  output$vbox_va <- renderValueBox({
    val <- round(get_last("BCVA_logMAR"), 3)
    valueBox(val, "Final logMAR", icon=icon("eye"), color="red")
  })
  output$vbox_cst <- renderValueBox({
    val <- round(get_last("CST_um"), 0)
    valueBox(paste0(val,"μm"), "OCT CST", icon=icon("images"), color="orange")
  })
  output$vbox_sun <- renderValueBox({
    val <- round(get_last("SUN_grade"), 2)
    valueBox(val, "SUN Grade", icon=icon("microscope"), color=if(!is.na(val)&&val<0.5)"green" else "yellow")
  })
  output$vbox_iop <- renderValueBox({
    val <- round(get_last("IOP_mmHg"), 1)
    valueBox(paste0(val,"mmHg"), "IOP", icon=icon("tachometer-alt"), color=if(!is.na(val)&&val>21)"red" else "blue")
  })

  # ── SUN Scale table ──────────────────────────────────────────────────────
  output$sun_scale <- renderTable({
    data.frame(Grade=c("0","0.5+","1+","2+","3+","4+"),
               Description=c("None","1-5 cells/HPF","6-15 cells","16-25 cells","26-50 cells",">50 cells"))
  }, striped=TRUE, hover=TRUE, bordered=TRUE)
  output$activity_score_tbl <- renderTable({
    d <- sim_data()
    if (nrow(d)==0) return(NULL)
    row_d90 <- d %>% filter(abs(time_d - 90) == min(abs(time_d - 90))) %>% slice(1)
    data.frame(Metric=c("SUN Grade","CST (μm)","logMAR","IOP (mmHg)","TNF-α"),
               Value=round(c(row_d90$SUN_grade, row_d90$CST_um, row_d90$BCVA_logMAR,
                              row_d90$IOP_mmHg, row_d90$TNF_level), 3))
  }, striped=TRUE)

  # ── PK plots ─────────────────────────────────────────────────────────────
  output$pk_plasma <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x=time_d, y=Cp)) + geom_line(color="#2980B9", size=1.2) +
      labs(title="Plasma Drug Concentration", x="Time (days)", y="Cp (AU)") + theme_bw(13)
  })
  output$pk_ocular <- renderPlot({
    d <- sim_data() %>% select(time_d, C_ant, C_post) %>%
      pivot_longer(-time_d, names_to="compartment", values_to="conc") %>%
      mutate(compartment=recode(compartment, C_ant="Anterior Chamber", C_post="Posterior Segment"))
    ggplot(d, aes(x=time_d, y=conc, color=compartment)) + geom_line(size=1.2) +
      scale_color_manual(values=c("#E74C3C","#27AE60")) +
      labs(title="Ocular Drug Concentration", x="Time (days)", y="Concentration (AU)", color="") +
      theme_bw(13) + theme(legend.position="bottom")
  })
  output$pk_gr <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x=time_d, y=GR_occupancy)) + geom_line(color="#8E44AD", size=1.2) +
      geom_hline(yintercept=0.5, linetype="dashed", color="gray50") +
      scale_y_continuous(limits=c(0,1.05)) +
      labs(title="Glucocorticoid Receptor Occupancy", x="Time (days)", y="GR Occupancy (0-1)") +
      theme_bw(13)
  })
  output$pk_summary_tbl <- renderTable({
    data.frame(
      Parameter=c("CL (L/h)","Vd (L)","ka (1/h)","t½ (h)","kBAB (1/h)","kBRB (1/h)"),
      Value=c(0.45,12.0,0.80,round(log(2)*12/0.45,1),0.003,0.001),
      Description=c("Systemic clearance","Volume of distribution",
                    "Absorption rate","Effective half-life",
                    "BAB transfer","BRB transfer")
    )
  }, striped=TRUE, bordered=TRUE)

  # ── Immune plots ─────────────────────────────────────────────────────────
  output$immune_tcell <- renderPlot({
    d <- sim_data() %>% select(time_d, Teff_cells, Treg_cells) %>%
      pivot_longer(-time_d, names_to="cell", values_to="level") %>%
      mutate(cell=recode(cell, Teff_cells="T Effector (Th1/Th17)", Treg_cells="T Regulatory"))
    ggplot(d, aes(x=time_d, y=level, color=cell)) + geom_line(size=1.2) +
      scale_color_manual(values=c("#E74C3C","#27AE60")) +
      labs(title="T Cell Populations", x="Time (days)", y="Relative Level (AU)", color="") +
      theme_bw(13) + theme(legend.position="bottom")
  })
  output$immune_apc <- renderPlot({
    d <- sim_data() %>% select(time_d, APC_act, Macro) %>%
      pivot_longer(-time_d, names_to="cell", values_to="level") %>%
      mutate(cell=recode(cell, APC_act="Activated APCs/DCs", Macro="M1 Macrophages"))
    ggplot(d, aes(x=time_d, y=level, color=cell)) + geom_line(size=1.2) +
      scale_color_manual(values=c("#F39C12","#8E44AD")) +
      labs(title="APC & Macrophage Activation", x="Time (days)", y="AU", color="") +
      theme_bw(13) + theme(legend.position="bottom")
  })
  output$cyto_tnf_il6 <- renderPlot({
    d <- sim_data() %>% filter(time_d <= min(90, input$sim_dur)) %>%
      select(time_d, TNF_level, IL6_level) %>%
      pivot_longer(-time_d, names_to="cyto", values_to="level")
    ggplot(d, aes(x=time_d, y=level, color=cyto)) + geom_line(size=1.2) +
      scale_color_manual(values=c("#E74C3C","#3498DB")) +
      labs(title="TNF-α & IL-6 (0-90 days)", x="Days", y="AU", color="") +
      theme_bw(13) + theme(legend.position="bottom")
  })
  output$cyto_il17_vegf <- renderPlot({
    d <- sim_data() %>% filter(time_d <= min(90, input$sim_dur)) %>%
      select(time_d, IL17=IL17_level, VEGF=VEGF_level) %>%
      pivot_longer(-time_d, names_to="cyto", values_to="level")
    ggplot(d, aes(x=time_d, y=level, color=cyto)) + geom_line(size=1.2) +
      scale_color_manual(values=c("#E67E22","#1ABC9C")) +
      labs(title="IL-17A & VEGF (0-90 days)", x="Days", y="AU", color="") +
      theme_bw(13) + theme(legend.position="bottom")
  })

  # ── Barrier plots ─────────────────────────────────────────────────────────
  output$bab_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x=time_d, y=BAB)) + geom_line(color="#2ECC71", size=1.2) +
      geom_hline(yintercept=0.8, linetype="dashed", color="orange") +
      scale_y_continuous(limits=c(0,1.05)) +
      labs(title="BAB Integrity", x="Time (days)", y="Integrity (0=disrupted, 1=intact)") +
      theme_bw(13)
  })
  output$brb_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x=time_d, y=BRB)) + geom_line(color="#1ABC9C", size=1.2) +
      geom_hline(yintercept=0.8, linetype="dashed", color="orange") +
      scale_y_continuous(limits=c(0,1.05)) +
      labs(title="BRB Integrity", x="Time (days)", y="Integrity (0=disrupted, 1=intact)") +
      theme_bw(13)
  })
  output$cells_ah_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x=time_d, y=SUN_grade)) + geom_line(color="#E74C3C", size=1.2) +
      geom_hline(yintercept=0.5, linetype="dashed", color="green") +
      annotate("text", x=5, y=0.6, label="Quiescence threshold", size=3, color="green4") +
      scale_y_continuous(limits=c(0,4.5), breaks=0:4) +
      labs(title="SUN Anterior Chamber Cells Grade", x="Time (days)", y="SUN Grade") +
      theme_bw(13)
  })

  # ── Clinical endpoint plots ───────────────────────────────────────────────
  output$va_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x=time_d, y=BCVA_logMAR)) + geom_line(color="#E74C3C", size=1.2) +
      labs(title="Visual Acuity Deficit", x="Time (days)", y="logMAR") + theme_bw(13)
  })
  output$cst_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x=time_d, y=CST_um)) + geom_line(color="#F39C12", size=1.2) +
      geom_hline(yintercept=300, linetype="dashed", color="green4") +
      labs(title="OCT CST (CME)", x="Time (days)", y="CST (μm)") + theme_bw(13)
  })
  output$iop_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x=time_d, y=IOP_mmHg)) + geom_line(color="#2980B9", size=1.2) +
      geom_hline(yintercept=21, linetype="dashed", color="red") +
      labs(title="Intraocular Pressure", x="Time (days)", y="IOP (mmHg)") + theme_bw(13)
  })
  output$milestones_tbl <- renderDT({
    d <- sim_data()
    tps <- c(7, 30, 60, 90, 180, 365)
    tps_avail <- tps[tps <= input$sim_dur]
    rows <- bind_rows(lapply(tps_avail, function(tp) {
      r <- d %>% filter(abs(time_d - tp) == min(abs(time_d - tp))) %>% slice(1)
      data.frame(Day=tp, logMAR=round(r$BCVA_logMAR,3), CST=round(r$CST_um,0),
                 SUN=round(r$SUN_grade,2), IOP=round(r$IOP_mmHg,1),
                 TNF=round(r$TNF_level,3), BAB=round(r$BAB,3))
    }))
    datatable(rows, options=list(dom='t', pageLength=10), rownames=FALSE)
  })

  # ── Scenario comparison ───────────────────────────────────────────────────
  output$scen_va <- renderPlot({
    d <- all_scen_data()
    ggplot(d, aes(x=time_d, y=BCVA_logMAR, color=scenario)) + geom_line(size=0.9) +
      scale_color_manual(values=colors7) +
      labs(title="Visual Acuity — All Scenarios", x="Days", y="logMAR", color="") +
      theme_bw(13) + theme(legend.position="bottom") + guides(color=guide_legend(ncol=2))
  })
  output$scen_cst <- renderPlot({
    d <- all_scen_data()
    ggplot(d, aes(x=time_d, y=CST_um, color=scenario)) + geom_line(size=0.9) +
      scale_color_manual(values=colors7) +
      geom_hline(yintercept=300, linetype="dashed") +
      labs(title="OCT CST — All Scenarios", x="Days", y="CST (μm)", color="") +
      theme_bw(13) + theme(legend.position="bottom") + guides(color=guide_legend(ncol=2))
  })
  output$scen_tnf <- renderPlot({
    d <- all_scen_data() %>% filter(time_d <= 90)
    ggplot(d, aes(x=time_d, y=TNF_level, color=scenario)) + geom_line(size=0.9) +
      scale_color_manual(values=colors7) +
      labs(title="TNF-α (0-90d)", x="Days", y="TNF (AU)", color="") +
      theme_bw(13) + theme(legend.position="bottom") + guides(color=guide_legend(ncol=2))
  })
  output$scen_tbl <- renderDT({
    d <- all_scen_data() %>% filter(abs(time_d - 90) == min(abs(time_d - 90))) %>%
      group_by(scenario) %>% slice(1) %>%
      select(scenario, BCVA_logMAR, CST_um, SUN_grade, IOP_mmHg, TNF_level, BAB, BRB) %>%
      mutate(across(where(is.numeric), ~round(.x, 3)))
    datatable(as.data.frame(d), options=list(dom='t'), rownames=FALSE)
  })

  # ── Virtual Patients ──────────────────────────────────────────────────────
  vp_sim <- eventReactive(input$run_vp_btn, {
    n   <- input$vp_n
    iiv <- input$vp_iiv / 100
    set.seed(123)
    params_base <- make_params()
    vp_df <- tibble(
      ID      = 1:n,
      CL_i    = rlnorm(n, log(0.45), iiv),
      Vd_i    = rlnorm(n, log(12), iiv*0.8),
      ksyn_TNF_i = rlnorm(n, log(0.05), iiv*1.2),
      Ag_i    = rlnorm(n, log(0.5), iiv),
      subtype = sample(c("Th1-dominant","Th17-dominant","Mixed"),n,
                       replace=TRUE, prob=c(0.4,0.35,0.25))
    )
    scen_ev <- switch(input$vp_scen,
      S1 = ev(), S5 = ev(amt=60, cmt="Cgut", ii=24, addl=89, time=0),
      S6 = ev(amt=40, cmt="Cp", ii=336, addl=11, time=0),
      S7 = ev(amt=60, cmt="Cgut", ii=24, addl=89, time=0) +
           ev(amt=40, cmt="Cp", ii=336, addl=11, time=0)
    )
    bind_rows(lapply(1:n, function(i) {
      p  <- vp_df[i,]
      pl <- c(params_base, list(CL=p$CL_i, Vd=p$Vd_i,
                                ksyn_TNF=p$ksyn_TNF_i, Ag_stim=p$Ag_i))
      if (p$subtype == "Th17-dominant") pl$ksyn_IL17 <- pl$ksyn_TNF * 1.5
      out <- mrgsim(param_set(mod, pl), ev=scen_ev,
                    end=180*24, delta=24) %>% as_tibble() %>%
        mutate(time_d=time/24, VP_ID=i, subtype=p$subtype)
      out %>% filter(time_d %in% c(0,30,90,180))
    }))
  }, ignoreNULL=FALSE)

  output$vp_bcva <- renderPlot({
    d <- vp_sim() %>%
      group_by(subtype, time_d) %>%
      summarise(med=median(BCVA_logMAR), q25=quantile(BCVA_logMAR,0.25),
                q75=quantile(BCVA_logMAR,0.75), .groups="drop")
    ggplot(d, aes(x=time_d, y=med, color=subtype, fill=subtype)) +
      geom_line(size=1.3) +
      geom_ribbon(aes(ymin=q25, ymax=q75), alpha=0.2, color=NA) +
      scale_color_manual(values=c("#E74C3C","#2980B9","#27AE60")) +
      scale_fill_manual(values=c("#E74C3C","#2980B9","#27AE60")) +
      labs(title="BCVA Response by Disease Subtype (Median ± IQR)",
           x="Days", y="logMAR", color="", fill="") +
      theme_bw(14)
  })
  output$vp_responder <- renderPlot({
    d <- vp_sim() %>% filter(time_d %in% c(30,90,180)) %>%
      group_by(subtype, time_d) %>%
      summarise(resp=mean(SUN_grade < 0.5), .groups="drop") %>%
      mutate(time_d=factor(paste0("Day ",time_d)))
    ggplot(d, aes(x=time_d, y=resp*100, fill=subtype)) +
      geom_bar(stat="identity", position="dodge") +
      scale_fill_manual(values=c("#E74C3C","#2980B9","#27AE60")) +
      labs(title="Responder Rate (SUN < 0.5)", x="", y="Responders (%)", fill="") +
      theme_bw(13)
  })
  output$vp_cme_var <- renderPlot({
    d <- vp_sim() %>% filter(time_d == 90)
    ggplot(d, aes(x=CST_um, fill=subtype)) +
      geom_histogram(bins=20, alpha=0.7, position="identity") +
      scale_fill_manual(values=c("#E74C3C","#2980B9","#27AE60")) +
      geom_vline(xintercept=300, linetype="dashed") +
      labs(title="Distribution of CST at Day 90", x="CST (μm)", y="Count", fill="") +
      theme_bw(13)
  })

  # ── Biomarker monitoring ──────────────────────────────────────────────────
  output$bm_cytokines <- renderPlot({
    d <- sim_data() %>% filter(time_d <= min(180, input$sim_dur)) %>%
      select(time_d, TNF_level, IL6_level, VEGF_level) %>%
      pivot_longer(-time_d, names_to="biomarker", values_to="level") %>%
      mutate(biomarker=recode(biomarker, TNF_level="TNF-α", IL6_level="IL-6", VEGF_level="VEGF"))
    ggplot(d, aes(x=time_d, y=level, color=biomarker)) + geom_line(size=1.1) +
      scale_color_manual(values=c("#E74C3C","#3498DB","#27AE60")) +
      labs(title="Cytokine Biomarkers (0-180d)", x="Days", y="AU", color="") +
      theme_bw(13) + theme(legend.position="bottom")
  })
  output$bm_oct <- renderPlot({
    d <- sim_data() %>%
      select(time_d, CST_um, BCVA_logMAR) %>%
      pivot_longer(-time_d, names_to="measure", values_to="value") %>%
      mutate(measure=recode(measure, CST_um="CST (μm)", BCVA_logMAR="logMAR×200"))
    d$value[d$measure=="logMAR×200"] <- d$value[d$measure=="logMAR×200"] * 200
    ggplot(d, aes(x=time_d, y=value, color=measure)) + geom_line(size=1.1) +
      scale_color_manual(values=c("#E67E22","#8E44AD")) +
      labs(title="OCT & Visual Function Monitoring", x="Days", y="Value", color="") +
      theme_bw(13) + theme(legend.position="bottom")
  })
  output$bm_drug_level <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x=time_d, y=DrugCp)) + geom_line(color="#1ABC9C", size=1.2) +
      geom_hline(yintercept=5, linetype="dashed", color="orange") +
      annotate("text", x=10, y=5.5, label="Therapeutic trough target", size=3, color="orange") +
      labs(title="Drug Plasma Level (Trough Monitoring)", x="Days", y="Cp (AU)") +
      theme_bw(13)
  })
  output$monitor_schedule <- renderTable({
    data.frame(
      Timepoint=c("Baseline","Week 2","Week 4","Week 8","Week 12","Every 3mo (stable)"),
      `Visual Acuity`=rep("✓",6), `IOP`=rep("✓",6),
      `OCT-CST`=c("✓","✓","✓","✓","✓","✓"),
      `SUN Grade`=c("✓","✓","✓","✓","✓","✓"),
      `Drug Level`=c("—","—","—","✓","✓","✓"),
      `ADA Ab`=c("—","—","—","—","✓","✓")
    )
  }, striped=TRUE, bordered=TRUE)
}

shinyApp(ui, server)
