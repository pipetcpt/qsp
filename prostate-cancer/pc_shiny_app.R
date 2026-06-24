################################################################################
# Prostate Cancer QSP Model — Interactive Shiny Dashboard
# ============================================================
# Tabs:
#   1. Patient Profile & Disease Overview
#   2. Drug PK Profiles (GnRH agents, ARPI, Chemo)
#   3. HPG Axis & Hormone Levels (T, DHT, LH)
#   4. AR Signaling & PSA Dynamics
#   5. Tumor Kinetics & Resistance
#   6. Bone Metastasis & Biomarkers
#   7. Scenario Comparison & Clinical Endpoints
#   8. Sensitivity Analysis (PTEN, HRR, AR)
################################################################################

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(shinythemes)
library(DT)

# ==============================================================================
# MODEL CODE (inline)
# ==============================================================================

model_code <- '
$PARAM @annotated
kLH_prod:8.0:LH production (IU/L/day)
kLH_deg:0.96:LH degradation (/day)
GnRH_base:1.0:Baseline GnRH tone
kT_prod:0.25:T production (nmol/L/day per LH)
kT_deg:0.35:T degradation (/day)
T_adrenal:0.05:Adrenal T contribution
f5alpha:0.10:T to DHT conversion fraction
kDHT_deg:0.50:DHT degradation (/day)
T_baseline:15.0:Baseline T (nmol/L)
kAR_synth:0.05:AR synthesis
kAR_deg:0.05:AR basal degradation
kon_AR:2.0:DHT-AR on-rate
koff_AR:0.5:DHT-AR off-rate
k_nuc:1.5:AR nuclear translocation
k_nuc_off:0.3:AR nuclear export
kAR_nuc_deg:0.2:Nuclear AR degradation
PSA_kprod:0.002:PSA production rate
kPSA_deg:0.10:PSA degradation (/day)
k_prolif:0.06:Tumor proliferation (/day)
k_death:0.01:Tumor death (/day)
k_quiesce:0.02:Entry to quiescence
k_unquiesce:0.015:Exit quiescence
k_death_q:0.005:Quiescent cell death
TC_cap:1000.0:Tumor carrying capacity
AR_prolif_EC50:0.5:AR EC50 for proliferation
PTEN_loss:0.7:PTEN loss fraction (0-1)
kAKT_base:0.3:Baseline AKT activity
kAKT_max:1.0:Max AKT activity
k_AKT_AR:0.3:AKT-AR cross-talk
k_AKT_BCL2:0.2:AKT-BCL2 anti-apoptosis
kOC_form:0.05:Osteoclast formation
kOC_deg:0.15:Osteoclast degradation
kOB_form:0.04:Osteoblast formation
kOB_deg:0.12:Osteoblast degradation
kRANKL:0.8:RANKL-OC activation
kOPG:0.4:OPG inhibition
kBMD_form:0.003:BMD formation
kBMD_resorb:0.005:BMD resorption
k_bonehom:0.02:Bone homing rate
kLeup_rel:0.033:Leuprolide release (/day)
kLeup_elim:0.693:Leuprolide elimination
V_Leup:40.0:Leuprolide Vd (L)
GnRH_flare:3.0:GnRH flare multiplier
flare_decay:0.5:Flare decay (/day)
kDeg_abs:0.15:Degarelix absorption
kDeg_elim:0.023:Degarelix elimination
V_Deg:1000.0:Degarelix Vd (L)
Deg_EC50:0.001:Degarelix EC50 (ug/mL)
Deg_Emax:0.98:Degarelix Emax
kEnz_abs:1.5:Enzalutamide absorption
kEnz_elim:0.114:Enzalutamide elimination
V_Enz:110.0:Enzalutamide Vd
F_Enz:0.84:Enzalutamide bioavailability
Enz_EC50:3.0:Enzalutamide EC50 (uM)
Enz_Emax:0.95:Enzalutamide Emax
kAbi_abs:0.8:Abiraterone absorption
kAbi_elim:1.7:Abiraterone elimination
V_Abi:19669.0:Abiraterone Vd
F_Abi:0.10:Abiraterone bioavailability
Abi_EC50:0.05:Abiraterone EC50 (uM)
Abi_Emax:0.95:Abiraterone Emax
Abi_MW:391.6:Abiraterone MW
kDoc_elim1:3.94:Docetaxel alpha elimination
kDoc_elim2:0.231:Docetaxel beta elimination
kDoc_k12:1.5:Docetaxel k12
kDoc_k21:0.8:Docetaxel k21
V_Doc1:6.0:Docetaxel central Vd (L)
Doc_EC50:0.05:Docetaxel EC50 (uM)
Doc_Emax:0.90:Docetaxel Emax
kOla_abs:1.4:Olaparib absorption
kOla_elim:1.6:Olaparib elimination
V_Ola:167.0:Olaparib Vd
F_Ola:0.66:Olaparib bioavailability
Ola_EC50:0.1:Olaparib EC50 (uM)
Ola_Emax:0.80:Olaparib Emax
HRR_def:0.0:HRR deficiency (0=proficient)
kDen_abs:0.062:Denosumab absorption
kDen_elim:0.023:Denosumab elimination
V_Den:3.0:Denosumab Vd (L)
Den_EC50:0.5:Denosumab EC50 (ug/mL)
Den_Emax:0.95:Denosumab Emax
k_CRPC:0.003:CRPC acquisition rate
k_ARv7:0.002:ARv7 emergence rate

$CMT @annotated
LH:LH (IU/L)
T:Testosterone (nmol/L)
DHT:DHT (nmol/L)
AR_free:Free AR
AR_DHT:AR-DHT complex
AR_nuc:Nuclear AR
PSA:PSA (ng/mL)
TC_p:Proliferating tumor cells
TC_q:Quiescent tumor cells
CRPC_frac:CRPC fraction
ARv7_frac:ARv7 fraction
AKT_act:AKT activity
OC:Osteoclasts
OB:Osteoblasts
BMD:Bone mineral density
BoneMets:Bone metastasis burden
Leup_depot:Leuprolide depot (mg)
Leup_c:Leuprolide plasma (ng/mL)
Flare_eff:GnRH flare
Deg_sc:Degarelix SC (mg)
Deg_c:Degarelix plasma (ug/mL)
Enz_gut:Enzalutamide gut (mg)
Enz_c:Enzalutamide plasma (uM)
Abi_gut:Abiraterone gut (mg)
Abi_c:Abiraterone plasma (uM)
Doc_c:Docetaxel central (uM)
Doc_p:Docetaxel peripheral (uM)
Ola_gut:Olaparib gut (mg)
Ola_c:Olaparib plasma (uM)
Den_sc:Denosumab SC (mg)
Den_c:Denosumab plasma (ug/mL)

$INIT
LH=5.0
T=15.0
DHT=1.5
AR_free=1.0
AR_DHT=0.5
AR_nuc=0.3
PSA=4.0
TC_p=1.0
TC_q=0.2
CRPC_frac=0.01
ARv7_frac=0.0
AKT_act=0.0
OC=1.0
OB=1.0
BMD=1.0
BoneMets=0.0
Leup_depot=0.0
Leup_c=0.0
Flare_eff=0.0
Deg_sc=0.0
Deg_c=0.0
Enz_gut=0.0
Enz_c=0.0
Abi_gut=0.0
Abi_c=0.0
Doc_c=0.0
Doc_p=0.0
Ola_gut=0.0
Ola_c=0.0
Den_sc=0.0
Den_c=0.0

$MAIN
double Leup_suppress = Leup_c/(Leup_c+2.0);
double Enz_AR_inh = Enz_Emax*Enz_c/(Enz_c+Enz_EC50);
double Abi_T_inh = Abi_Emax*Abi_c/(Abi_c+Abi_EC50);
double Doc_kill = Doc_Emax*Doc_c/(Doc_c+Doc_EC50);
double Ola_kill = HRR_def*Ola_Emax*Ola_c/(Ola_c+Ola_EC50);
double Den_RANKL_inh = Den_Emax*Den_c/(Den_c+Den_EC50);
double Deg_blockade = Deg_Emax*Deg_c/(Deg_c+Deg_EC50);
double GnRH_agonist_effect = (1.0+GnRH_flare*Flare_eff)*(1.0-0.97*Leup_suppress);
double GnRH_total = GnRH_base*(1.0-Deg_blockade)*GnRH_agonist_effect;
double AKT_ss = kAKT_base+(kAKT_max-kAKT_base)*PTEN_loss*TC_p/(TC_p+0.5);
double AR_nuc_eff = AR_nuc*(1.0-Enz_AR_inh)+ARv7_frac*0.5;
double AR_nuc_norm = AR_nuc_eff/(AR_nuc_eff+0.3);
double prolif_eff = fmax(0.0,AR_nuc_norm+k_AKT_AR*AKT_act);
double k_prolif_eff = k_prolif*prolif_eff;
double k_death_eff = k_death*(1.0+Doc_kill+Ola_kill)*(1.0+k_AKT_BCL2*(1.0-AKT_act));
double RANKL_eff = kRANKL*BoneMets*(1.0-Den_RANKL_inh)/(1.0+kOPG*OB);

$ODE
dxdt_LH = kLH_prod*GnRH_total - kLH_deg*LH;
dxdt_T = kT_prod*LH*(1.0-Abi_T_inh)+T_adrenal - kT_deg*T;
dxdt_DHT = f5alpha*T - kDHT_deg*DHT;
double AR_bind = kon_AR*DHT*AR_free*(1.0-Enz_AR_inh);
double AR_unbind = koff_AR*AR_DHT;
dxdt_AR_free = kAR_synth - kAR_deg*AR_free - AR_bind + AR_unbind;
dxdt_AR_DHT = AR_bind - AR_unbind - k_nuc*AR_DHT;
dxdt_AR_nuc = k_nuc*AR_DHT - (k_nuc_off+kAR_nuc_deg)*AR_nuc;
double TC_total = TC_p+TC_q;
dxdt_PSA = PSA_kprod*AR_nuc_eff*TC_total - kPSA_deg*PSA;
dxdt_AKT_act = 5.0*(AKT_ss-AKT_act);
double logistic = 1.0-TC_total/TC_cap;
double CRPC_prolif = CRPC_frac*k_prolif*0.8;
double net_prolif = (k_prolif_eff+CRPC_prolif)*logistic;
dxdt_TC_p = net_prolif*TC_p - k_death_eff*TC_p - k_quiesce*TC_p + k_unquiesce*TC_q;
dxdt_TC_q = k_quiesce*TC_p - k_unquiesce*TC_q - k_death_q*TC_q*(1.0+Doc_kill*0.5);
dxdt_CRPC_frac = k_CRPC*(1.0-CRPC_frac)*(1.0+Enz_AR_inh*2.0+Abi_T_inh*1.5);
dxdt_ARv7_frac = CRPC_frac*k_ARv7*(1.0-ARv7_frac)*(1.0+Enz_AR_inh*3.0);
dxdt_BoneMets = k_bonehom*TC_p*(1.0-BoneMets/10.0);
dxdt_OC = kOC_form*(1.0+RANKL_eff) - kOC_deg*OC;
dxdt_OB = kOB_form*(1.0+BoneMets*0.5) - kOB_deg*OB;
dxdt_BMD = kBMD_form*OB - kBMD_resorb*OC;
dxdt_Leup_depot = -kLeup_rel*Leup_depot;
dxdt_Leup_c = kLeup_rel*Leup_depot/V_Leup*1000.0 - kLeup_elim*Leup_c;
dxdt_Flare_eff = -flare_decay*Flare_eff;
dxdt_Deg_sc = -kDeg_abs*Deg_sc;
dxdt_Deg_c = kDeg_abs*Deg_sc/V_Den - kDeg_elim*Deg_c;
dxdt_Enz_gut = -kEnz_abs*Enz_gut;
dxdt_Enz_c = kEnz_abs*F_Enz*Enz_gut/V_Enz - kEnz_elim*Enz_c;
dxdt_Abi_gut = -kAbi_abs*Abi_gut;
dxdt_Abi_c = kAbi_abs*F_Abi*Abi_gut/V_Abi*1e6/Abi_MW - kAbi_elim*Abi_c;
dxdt_Doc_c = -(kDoc_elim1+kDoc_k12)*Doc_c + kDoc_k21*Doc_p;
dxdt_Doc_p = kDoc_k12*Doc_c - (kDoc_k21+kDoc_elim2)*Doc_p;
dxdt_Ola_gut = -kOla_abs*Ola_gut;
dxdt_Ola_c = kOla_abs*F_Ola*Ola_gut/V_Ola - kOla_elim*Ola_c;
dxdt_Den_sc = -kDen_abs*Den_sc;
dxdt_Den_c = kDen_abs*Den_sc/V_Den - kDen_elim*Den_c;

$CAPTURE PSA T DHT AR_nuc AR_nuc CRPC_frac ARv7_frac
OC OB BMD BoneMets AKT_act Enz_c Abi_c Doc_c Leup_c Deg_c Ola_c TC_p TC_q
'

# Compile model
mod_base <- mcode("PC_Shiny", model_code, quiet = TRUE)

# Helper: build event table from UI inputs
build_events <- function(leup, deg, enz, abi, doc, ola, den,
                          sim_dur, body_sa = 1.8) {
  evs <- list()

  if (leup) {
    evs[[length(evs)+1]] <- data.frame(
      time = seq(0, sim_dur - 1, by = 28),
      cmt  = "Leup_depot", amt = 7.5, evid = 1
    )
  }
  if (deg) {
    evs[[length(evs)+1]] <- data.frame(
      time = 0, cmt = "Deg_sc", amt = 240, evid = 1
    )
    if (sim_dur > 28) {
      evs[[length(evs)+1]] <- data.frame(
        time = seq(28, sim_dur-1, by = 28),
        cmt  = "Deg_sc", amt = 80, evid = 1
      )
    }
  }
  if (enz) {
    evs[[length(evs)+1]] <- data.frame(
      time = seq(0, sim_dur - 1), cmt = "Enz_gut", amt = 160, evid = 1
    )
  }
  if (abi) {
    evs[[length(evs)+1]] <- data.frame(
      time = seq(0, sim_dur - 1), cmt = "Abi_gut", amt = 1000, evid = 1
    )
  }
  if (doc) {
    n_cycles <- min(6, floor(sim_dur / 21))
    doc_dose  <- 75 * body_sa  # mg
    doc_uM    <- doc_dose / 861.9 / 6.0 * 1e6
    evs[[length(evs)+1]] <- data.frame(
      time = seq(0, (n_cycles-1)*21, by = 21),
      cmt  = "Doc_c", amt = doc_uM, evid = 1
    )
  }
  if (ola) {
    evs[[length(evs)+1]] <- data.frame(
      time = sort(c(seq(0, sim_dur - 1), seq(0, sim_dur - 1) + 0.5)),
      cmt  = "Ola_gut", amt = 300, evid = 1
    )
  }
  if (den) {
    evs[[length(evs)+1]] <- data.frame(
      time = seq(0, sim_dur - 1, by = 28),
      cmt  = "Den_sc", amt = 120, evid = 1
    )
  }
  if (length(evs) == 0) return(ev())
  bind_rows(evs)
}

# ggplot theme
theme_qsp <- theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#2c3e50"),
    strip.text       = element_text(color = "white", face = "bold"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

# ==============================================================================
# UI
# ==============================================================================

ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel(
    div(
      h2("Prostate Cancer QSP Model Dashboard",
         style = "color:#2c3e50; font-weight:bold;"),
      h5("Quantitative Systems Pharmacology: HPG Axis · AR Signaling · Tumor Kinetics · Bone Metastasis",
         style = "color:#7f8c8d;")
    )
  ),

  navbarPage(
    title = "",
    id    = "navbar",

    # ==========================================
    # TAB 1: PATIENT PROFILE
    # ==========================================
    tabPanel("Patient Profile",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("Patient Characteristics", style = "color:#2c3e50;"),
          hr(),
          sliderInput("psa_init", "Baseline PSA (ng/mL)",
                      min = 0.1, max = 200, value = 4, step = 0.5),
          sliderInput("tc_init", "Initial Tumor Burden (normalized)",
                      min = 0.1, max = 10, value = 1, step = 0.1),
          selectInput("gleason", "Gleason Score / Grade Group",
                      choices = c("GG1 (≤6)"="gg1", "GG2 (3+4)"="gg2",
                                  "GG3 (4+3)"="gg3", "GG4 (8)"="gg4",
                                  "GG5 (9-10)"="gg5"),
                      selected = "gg3"),
          selectInput("disease_stage", "Disease Stage",
                      choices = c("Localised (T1-T2)"="local",
                                  "Locally Advanced (T3-T4)"="la",
                                  "Hormone-Sensitive Mets (mHSPC)"="mhspc",
                                  "Castration-Resistant (CRPC)"="crpc",
                                  "Metastatic CRPC (mCRPC)"="mcrpc"),
                      selected = "mhspc"),
          hr(),
          sliderInput("pten_loss", "PTEN Loss Fraction",
                      min = 0, max = 1, value = 0.7, step = 0.1),
          checkboxInput("brca2_mut", "BRCA2/HRR Mutation", value = FALSE),
          checkboxInput("arv7_pos",  "Baseline ARv7 Positive", value = FALSE),
          hr(),
          sliderInput("body_sa", "Body Surface Area (m²)",
                      min = 1.4, max = 2.4, value = 1.8, step = 0.05),
          sliderInput("sim_duration", "Simulation Duration (days)",
                      min = 180, max = 1825, value = 730, step = 30)
        ),
        mainPanel(
          width = 9,
          fluidRow(
            column(6,
              h4("Disease Overview: Prostate Cancer"),
              p("Prostate cancer (PCa) is the second most common malignancy in
                men globally. The androgen signaling axis—through the
                hypothalamic-pituitary-gonadal (HPG) axis—drives the majority of
                prostate cancer growth. Testosterone is converted to
                dihydrotestosterone (DHT) by 5α-reductase, which activates the
                androgen receptor (AR), driving PSA production and cell
                proliferation."),
              p("This QSP model captures the transition from hormone-sensitive
                to castration-resistant disease (CRPC), emergence of ARv7
                splice variants, and the role of PI3K/AKT pathway activation
                (PTEN loss) in treatment resistance. It integrates
                pharmacokinetics of 8 drug classes."),
              h5("Key Pathways Modeled:", style = "font-weight:bold;"),
              tags$ul(
                tags$li("HPG Axis: GnRH → LH → Testosterone → DHT"),
                tags$li("AR Signaling: DHT-AR binding → Nuclear AR → PSA"),
                tags$li("Tumor Kinetics: Proliferating + Quiescent cells"),
                tags$li("PI3K/AKT: PTEN loss → AKT → bypass AR"),
                tags$li("Bone Metastasis: RANKL/RANK/OPG vicious cycle"),
                tags$li("Resistance: CRPC fraction, ARv7 emergence")
              )
            ),
            column(6,
              h4("Patient Status Summary"),
              tableOutput("patient_summary_table"),
              br(),
              h4("Disease Stage Description"),
              uiOutput("stage_description")
            )
          ),
          hr(),
          fluidRow(
            column(12,
              h4("Mechanistic Map Preview"),
              tags$img(src = "pc_qsp_model.png",
                       width = "100%",
                       alt  = "Prostate Cancer QSP Mechanistic Map"),
              p("Full resolution: pc_qsp_model.svg",
                style = "color:gray; font-size:11px;")
            )
          )
        )
      )
    ),

    # ==========================================
    # TAB 2: DRUG PK PROFILES
    # ==========================================
    tabPanel("Drug PK",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("Treatment Selection"),
          hr(),
          h5("ADT (GnRH Axis)"),
          checkboxInput("use_leup", "Leuprolide 7.5 mg IM q28d", TRUE),
          checkboxInput("use_deg",  "Degarelix 240 mg SC (loading)", FALSE),
          hr(),
          h5("AR Pathway Inhibitors"),
          checkboxInput("use_enz", "Enzalutamide 160 mg QD", FALSE),
          checkboxInput("use_abi", "Abiraterone 1000 mg QD", FALSE),
          hr(),
          h5("Chemotherapy"),
          checkboxInput("use_doc", "Docetaxel 75 mg/m² q3w ×6", FALSE),
          hr(),
          h5("Other"),
          checkboxInput("use_ola", "Olaparib 300 mg BID", FALSE),
          checkboxInput("use_den", "Denosumab 120 mg SC q4w", FALSE),
          hr(),
          selectInput("pk_drug_view", "Show PK for:",
                      choices = c("Leuprolide (plasma)", "Degarelix (plasma)",
                                  "Enzalutamide (plasma μM)",
                                  "Abiraterone (plasma μM)",
                                  "Docetaxel (plasma μM)",
                                  "Olaparib (plasma μM)",
                                  "Denosumab (plasma μg/mL)"),
                      selected = "Leuprolide (plasma)"),
          actionButton("run_pk", "Run Simulation", class = "btn-primary")
        ),
        mainPanel(
          width = 9,
          plotOutput("pk_plot", height = "400px"),
          br(),
          plotOutput("hormone_pk_plot", height = "350px")
        )
      )
    ),

    # ==========================================
    # TAB 3: HPG AXIS & HORMONES
    # ==========================================
    tabPanel("HPG Axis & Hormones",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("ADT Regimen"),
          radioButtons("adt_type", "GnRH Drug:",
                       choices = c("None"="none",
                                   "Leuprolide (Agonist)"="leup",
                                   "Degarelix (Antagonist)"="deg"),
                       selected = "leup"),
          checkboxInput("add_abi_hor", "Add Abiraterone (adrenal T)", FALSE),
          hr(),
          sliderInput("t_init", "Baseline Testosterone (nmol/L)",
                      min = 5, max = 35, value = 15, step = 1),
          sliderInput("t_adrenal", "Adrenal T Contribution",
                      min = 0, max = 0.3, value = 0.05, step = 0.01),
          sliderInput("f5alpha", "5α-Reductase Activity",
                      min = 0.02, max = 0.3, value = 0.10, step = 0.01),
          actionButton("run_hor", "Update", class = "btn-primary")
        ),
        mainPanel(
          width = 9,
          fluidRow(
            column(6, plotOutput("lh_plot", height = "280px")),
            column(6, plotOutput("t_plot",  height = "280px"))
          ),
          fluidRow(
            column(6, plotOutput("dht_plot", height = "280px")),
            column(6, plotOutput("castrate_status", height = "280px"))
          )
        )
      )
    ),

    # ==========================================
    # TAB 4: AR SIGNALING & PSA
    # ==========================================
    tabPanel("AR Signaling & PSA",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("AR Pathway Parameters"),
          checkboxInput("use_leup_ar", "Leuprolide (ADT)", TRUE),
          checkboxInput("use_enz_ar",  "Enzalutamide (ARPI)", FALSE),
          checkboxInput("use_abi_ar",  "Abiraterone (CYP17A1)", FALSE),
          hr(),
          checkboxInput("arv7_ar", "ARv7-positive (resistance)", FALSE),
          sliderInput("arv7_frac_init", "Initial ARv7 Fraction",
                      min = 0, max = 0.5, value = 0, step = 0.05),
          hr(),
          sliderInput("enz_ec50", "Enzalutamide EC50 (μM)",
                      min = 0.5, max = 10, value = 3, step = 0.5),
          sliderInput("psa_init_ar", "Initial PSA (ng/mL)",
                      min = 0.5, max = 500, value = 10, step = 0.5),
          actionButton("run_ar", "Update AR Simulation", class = "btn-primary")
        ),
        mainPanel(
          width = 9,
          fluidRow(
            column(6, plotOutput("ar_nuc_plot", height = "280px")),
            column(6, plotOutput("psa_plot",    height = "280px"))
          ),
          fluidRow(
            column(6, plotOutput("psa_log_plot", height = "280px")),
            column(6, plotOutput("psa50_waterfall", height = "280px"))
          )
        )
      )
    ),

    # ==========================================
    # TAB 5: TUMOR KINETICS & RESISTANCE
    # ==========================================
    tabPanel("Tumor Kinetics & Resistance",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("Tumor Biology"),
          checkboxInput("use_leup_t",  "Leuprolide (ADT)",       TRUE),
          checkboxInput("use_enz_t",   "Enzalutamide (ARPI)",    FALSE),
          checkboxInput("use_doc_t",   "Docetaxel (Chemo)",      FALSE),
          checkboxInput("use_ola_t",   "Olaparib (PARPi)",       FALSE),
          hr(),
          sliderInput("pten_t", "PTEN Loss",
                      min = 0, max = 1, value = 0.7, step = 0.1),
          sliderInput("hrr_t",  "HRR Deficiency (0=no, 1=yes)",
                      min = 0, max = 1, value = 0, step = 0.1),
          sliderInput("k_prolif_t", "Tumor Proliferation Rate (/day)",
                      min = 0.02, max = 0.15, value = 0.06, step = 0.005),
          sliderInput("k_crpc_t", "CRPC Acquisition Rate",
                      min = 0.001, max = 0.01, value = 0.003, step = 0.0005),
          actionButton("run_tumor", "Update", class = "btn-primary")
        ),
        mainPanel(
          width = 9,
          fluidRow(
            column(6, plotOutput("tc_total_plot",  height = "280px")),
            column(6, plotOutput("crpc_plot",      height = "280px"))
          ),
          fluidRow(
            column(6, plotOutput("arv7_plot",      height = "280px")),
            column(6, plotOutput("akt_plot",       height = "280px"))
          )
        )
      )
    ),

    # ==========================================
    # TAB 6: BONE METASTASIS & BIOMARKERS
    # ==========================================
    tabPanel("Bone Metastasis",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("Bone Treatment"),
          checkboxInput("use_leup_b", "Leuprolide (ADT)", TRUE),
          checkboxInput("use_den_b",  "Denosumab (anti-RANKL)", FALSE),
          hr(),
          sliderInput("rankl_b", "RANKL Activity",
                      min = 0.1, max = 2.0, value = 0.8, step = 0.1),
          sliderInput("opg_b", "OPG Concentration",
                      min = 0.1, max = 1.0, value = 0.4, step = 0.05),
          sliderInput("k_bonehom_b", "Bone Homing Rate",
                      min = 0.005, max = 0.1, value = 0.02, step = 0.005),
          actionButton("run_bone", "Update Bone Simulation",
                       class = "btn-primary")
        ),
        mainPanel(
          width = 9,
          fluidRow(
            column(6, plotOutput("bmd_plot",      height = "280px")),
            column(6, plotOutput("bone_mets_plot", height = "280px"))
          ),
          fluidRow(
            column(6, plotOutput("oc_ob_plot",    height = "280px")),
            column(6, plotOutput("bone_markers",  height = "280px"))
          )
        )
      )
    ),

    # ==========================================
    # TAB 7: SCENARIO COMPARISON
    # ==========================================
    tabPanel("Scenario Comparison",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("Compare Treatment Scenarios"),
          checkboxGroupInput("scenarios_cmp",
            "Select Scenarios:",
            choices = c(
              "Untreated"              = "s0",
              "ADT Alone (Leuprolide)" = "s1",
              "ADT + Enzalutamide"     = "s2",
              "ADT + Abiraterone"      = "s3",
              "ADT + Docetaxel ×6"     = "s4",
              "ADT + Olaparib (HRR-def)"="s5",
              "Sequential ADT→ARPI→Doc"="s6"
            ),
            selected = c("s0","s1","s2","s3")
          ),
          selectInput("cmp_endpoint", "Primary Endpoint:",
                      choices = c("PSA (ng/mL)", "Testosterone (ng/dL)",
                                  "Tumor Burden", "BMD",
                                  "CRPC Fraction", "Bone Metastasis"),
                      selected = "PSA (ng/mL)"),
          hr(),
          sliderInput("cmp_pten", "PTEN Loss",
                      min = 0, max = 1, value = 0.7, step = 0.1),
          numericInput("cmp_dur", "Duration (days)", 730, 180, 1825, 30),
          actionButton("run_cmp", "Run Comparison", class = "btn-primary")
        ),
        mainPanel(
          width = 9,
          plotOutput("cmp_plot", height = "400px"),
          br(),
          h4("Clinical Endpoints Summary"),
          DT::dataTableOutput("endpoints_table")
        )
      )
    ),

    # ==========================================
    # TAB 8: SENSITIVITY ANALYSIS
    # ==========================================
    tabPanel("Sensitivity Analysis",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("Sensitivity Analysis Settings"),
          selectInput("sa_treatment", "Treatment Context:",
                      choices = c("ADT Alone"       = "adt",
                                  "ADT+Enzalutamide" = "adtenz",
                                  "Docetaxel"        = "doc"),
                      selected = "adtenz"),
          selectInput("sa_x_param", "X-axis Parameter:",
                      choices = c("PTEN Loss"           = "pten",
                                  "CRPC Acquisition Rate"= "k_crpc",
                                  "AR Proliferation EC50"= "ar_ec50",
                                  "Tumor Proliferation"  = "k_prolif",
                                  "Bone Homing Rate"     = "k_bonehom"),
                      selected = "pten"),
          selectInput("sa_color_param", "Color by:",
                      choices = c("HRR Deficiency"    = "hrr",
                                  "ARv7 Baseline"     = "arv7",
                                  "Gleason Score"     = "gleason"),
                      selected = "hrr"),
          selectInput("sa_endpoint", "Endpoint:",
                      choices = c("PSA at 3 months",
                                  "Tumor Burden at 1 year",
                                  "CRPC Fraction at 2 years",
                                  "BMD at 2 years"),
                      selected = "PSA at 3 months"),
          numericInput("sa_n", "# Parameter Values", 5, 3, 10, 1),
          actionButton("run_sa", "Run Sensitivity", class = "btn-primary")
        ),
        mainPanel(
          width = 9,
          plotOutput("sa_plot",     height = "380px"),
          br(),
          plotOutput("tornado_plot", height = "320px")
        )
      )
    )
  )  # end navbarPage
)

# ==============================================================================
# SERVER
# ==============================================================================

server <- function(input, output, session) {

  # ----- Helper: base simulation -----
  run_sim <- function(leup = FALSE, deg = FALSE, enz = FALSE, abi = FALSE,
                      doc = FALSE, ola = FALSE, den = FALSE,
                      duration = 730, pars = list()) {
    evts <- build_events(leup, deg, enz, abi, doc, ola, den, duration,
                          input$body_sa)
    sim_mod <- mod_base
    if (length(pars) > 0) sim_mod <- param(sim_mod, .x = pars)

    init_vals <- init(sim_mod)
    if (!is.null(input$psa_init)) init_vals["PSA"] <- input$psa_init
    if (!is.null(input$tc_init))  init_vals["TC_p"] <- input$tc_init

    sim_mod <- init(sim_mod, init_vals)
    mrgsim(sim_mod, events = evts, end = duration, delta = 1) %>%
      as_tibble()
  }

  # ----- TAB 1: Patient Summary -----
  output$patient_summary_table <- renderTable({
    data.frame(
      Parameter = c("Baseline PSA", "Tumor Burden", "Disease Stage",
                    "PTEN Status", "HRR Status", "Simulation Duration"),
      Value     = c(
        paste(input$psa_init, "ng/mL"),
        paste(input$tc_init, "(normalized)"),
        input$disease_stage,
        paste0(round(input$pten_loss * 100), "% loss"),
        ifelse(input$brca2_mut, "Deficient (BRCA2+)", "Proficient"),
        paste(input$sim_duration, "days")
      )
    )
  }, bordered = TRUE, striped = TRUE)

  output$stage_description <- renderUI({
    desc <- switch(input$disease_stage,
      "local"  = "Localised PCa: organ-confined, radical prostatectomy or RT curative.",
      "la"     = "Locally Advanced: extends beyond capsule (T3-T4), risk of systemic spread.",
      "mhspc"  = "mHSPC: metastatic but still testosterone-sensitive; ADT is backbone therapy.",
      "crpc"   = "CRPC: PSA/clinical progression despite castrate levels; AR still active.",
      "mcrpc"  = "mCRPC: metastatic, castration-resistant; life expectancy ~18-24 months."
    )
    tags$div(class = "alert alert-info", desc)
  })

  # ----- TAB 2: Drug PK -----
  pk_sim <- eventReactive(input$run_pk, {
    run_sim(
      leup = input$use_leup, deg = input$use_deg,
      enz  = input$use_enz,  abi = input$use_abi,
      doc  = input$use_doc,  ola = input$use_ola,
      den  = input$use_den,
      duration = input$sim_duration
    )
  }, ignoreNULL = FALSE)

  output$pk_plot <- renderPlot({
    df <- pk_sim()
    drug_var <- switch(input$pk_drug_view,
      "Leuprolide (plasma)"       = "Leup_c",
      "Degarelix (plasma)"        = "Deg_c",
      "Enzalutamide (plasma μM)"  = "Enz_c",
      "Abiraterone (plasma μM)"   = "Abi_c",
      "Docetaxel (plasma μM)"     = "Doc_c",
      "Olaparib (plasma μM)"      = "Ola_c",
      "Denosumab (plasma μg/mL)"  = "Den_c"
    )
    if (!drug_var %in% names(df)) return(NULL)

    ggplot(df, aes(x = time / 30.4, y = .data[[drug_var]])) +
      geom_line(color = "#2980b9", linewidth = 1.2) +
      labs(title = paste("PK Profile:", input$pk_drug_view),
           x = "Time (months)", y = input$pk_drug_view) +
      theme_qsp
  })

  output$hormone_pk_plot <- renderPlot({
    df <- pk_sim() %>%
      mutate(T_ngdL = T * 28.84) %>%
      pivot_longer(cols = c(LH, T_ngdL, DHT),
                   names_to = "Hormone", values_to = "Conc") %>%
      mutate(Hormone = recode(Hormone,
        LH     = "LH (IU/L)",
        T_ngdL = "Testosterone (ng/dL)",
        DHT    = "DHT (nmol/L)"))

    ggplot(df, aes(x = time / 30.4, y = Conc)) +
      geom_line(color = "#e74c3c", linewidth = 1.0) +
      geom_hline(data = data.frame(Hormone = "Testosterone (ng/dL)",
                                    yint = 50),
                 aes(yintercept = yint), linetype = "dashed", color = "gray") +
      facet_wrap(~ Hormone, scales = "free_y") +
      labs(title = "Hormone Levels Under Treatment",
           x = "Time (months)", y = "Concentration") +
      theme_qsp
  })

  # ----- TAB 3: HPG Axis -----
  hor_sim <- eventReactive(input$run_hor, {
    run_sim(
      leup = input$adt_type == "leup",
      deg  = input$adt_type == "deg",
      abi  = input$add_abi_hor,
      duration = input$sim_duration,
      pars = list(T_adrenal = input$t_adrenal,
                  f5alpha    = input$f5alpha)
    )
  }, ignoreNULL = FALSE)

  output$lh_plot <- renderPlot({
    ggplot(hor_sim(), aes(x = time / 30.4, y = LH)) +
      geom_line(color = "#3498db", linewidth = 1.1) +
      labs(title = "LH Dynamics", x = "Months", y = "LH (IU/L)") +
      theme_qsp
  })

  output$t_plot <- renderPlot({
    ggplot(hor_sim(), aes(x = time / 30.4, y = T * 28.84)) +
      geom_line(color = "#e74c3c", linewidth = 1.1) +
      geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
      annotate("text", x = 1, y = 55,
               label = "Castrate threshold", hjust = 0, size = 3, color = "red") +
      labs(title = "Testosterone", x = "Months", y = "Testosterone (ng/dL)") +
      theme_qsp
  })

  output$dht_plot <- renderPlot({
    ggplot(hor_sim(), aes(x = time / 30.4, y = DHT)) +
      geom_line(color = "#f39c12", linewidth = 1.1) +
      labs(title = "DHT", x = "Months", y = "DHT (nmol/L)") +
      theme_qsp
  })

  output$castrate_status <- renderPlot({
    df <- hor_sim() %>%
      mutate(Castrate = T * 28.84 < 50)
    ggplot(df, aes(x = time / 30.4, y = T * 28.84, color = Castrate)) +
      geom_line(linewidth = 1.2) +
      geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = c("FALSE" = "#e74c3c", "TRUE" = "#27ae60"),
                         labels = c("FALSE" = "Non-castrate", "TRUE" = "Castrate")) +
      labs(title = "Castration Status", x = "Months", y = "Testosterone (ng/dL)",
           color = "") +
      theme_qsp
  })

  # ----- TAB 4: AR Signaling & PSA -----
  ar_sim <- eventReactive(input$run_ar, {
    run_sim(
      leup = input$use_leup_ar,
      enz  = input$use_enz_ar,
      abi  = input$use_abi_ar,
      duration = input$sim_duration,
      pars = list(
        Enz_EC50   = input$enz_ec50,
        ARv7_frac  = ifelse(input$arv7_ar, input$arv7_frac_init, 0)
      )
    )
  }, ignoreNULL = FALSE)

  output$ar_nuc_plot <- renderPlot({
    ggplot(ar_sim(), aes(x = time / 30.4, y = AR_nuc)) +
      geom_line(color = "#9b59b6", linewidth = 1.1) +
      labs(title = "Nuclear AR Activity", x = "Months", y = "Nuclear AR (normalized)") +
      theme_qsp
  })

  output$psa_plot <- renderPlot({
    ggplot(ar_sim(), aes(x = time / 30.4, y = PSA)) +
      geom_line(color = "#8e44ad", linewidth = 1.1) +
      geom_hline(yintercept = 4, linetype = "dashed", color = "gray50") +
      labs(title = "PSA Dynamics (linear)", x = "Months", y = "PSA (ng/mL)") +
      theme_qsp
  })

  output$psa_log_plot <- renderPlot({
    ggplot(ar_sim(), aes(x = time / 30.4, y = PSA)) +
      geom_line(color = "#6c3483", linewidth = 1.1) +
      scale_y_log10() +
      geom_hline(yintercept = 4, linetype = "dashed", color = "gray50") +
      labs(title = "PSA Dynamics (log scale)", x = "Months", y = "PSA (ng/mL, log)") +
      theme_qsp
  })

  output$psa50_waterfall <- renderPlot({
    df <- ar_sim()
    psa_base <- df$PSA[df$time == 0][1]
    psa_chg  <- (df$PSA - psa_base) / psa_base * 100

    bar_times <- c(90, 180, 365, min(input$sim_duration, df$time))
    bar_chg   <- sapply(bar_times, function(t) {
      idx <- which.min(abs(df$time - t))
      psa_chg[idx]
    })

    bar_df <- data.frame(
      Time  = paste0("Day ", bar_times),
      Change = bar_chg
    )

    ggplot(bar_df, aes(x = Time, y = Change, fill = Change < -50)) +
      geom_col() +
      geom_hline(yintercept = -50, linetype = "dashed", color = "red") +
      scale_fill_manual(values = c("FALSE" = "#e74c3c", "TRUE" = "#27ae60"),
                        guide = "none") +
      labs(title = "PSA Change from Baseline (%)",
           x = "", y = "% Change in PSA") +
      theme_qsp
  })

  # ----- TAB 5: Tumor Kinetics -----
  tumor_sim <- eventReactive(input$run_tumor, {
    run_sim(
      leup = input$use_leup_t, enz = input$use_enz_t,
      doc  = input$use_doc_t,  ola = input$use_ola_t,
      duration = input$sim_duration,
      pars = list(
        PTEN_loss = input$pten_t,
        HRR_def   = input$hrr_t,
        k_prolif  = input$k_prolif_t,
        k_CRPC    = input$k_crpc_t
      )
    )
  }, ignoreNULL = FALSE)

  output$tc_total_plot <- renderPlot({
    df <- tumor_sim() %>% mutate(TC_total = TC_p + TC_q)
    ggplot(df, aes(x = time / 30.4, y = TC_total)) +
      geom_line(color = "#27ae60", linewidth = 1.1) +
      geom_ribbon(aes(ymin = TC_q, ymax = TC_total), alpha = 0.2, fill = "#27ae60") +
      geom_ribbon(aes(ymin = 0, ymax = TC_q), alpha = 0.3, fill = "#f39c12") +
      labs(title = "Tumor Burden (green=prolif, orange=quiescent)",
           x = "Months", y = "Tumor Cells (normalized)") +
      theme_qsp
  })

  output$crpc_plot <- renderPlot({
    ggplot(tumor_sim(), aes(x = time / 30.4, y = CRPC_frac * 100)) +
      geom_line(color = "#c0392b", linewidth = 1.1) +
      labs(title = "CRPC Subpopulation Emergence",
           x = "Months", y = "CRPC Fraction (%)") +
      theme_qsp
  })

  output$arv7_plot <- renderPlot({
    ggplot(tumor_sim(), aes(x = time / 30.4, y = ARv7_frac * 100)) +
      geom_line(color = "#8e44ad", linewidth = 1.1) +
      labs(title = "ARv7-Positive Fraction",
           x = "Months", y = "ARv7+ Cells (%)") +
      theme_qsp
  })

  output$akt_plot <- renderPlot({
    ggplot(tumor_sim(), aes(x = time / 30.4, y = AKT_act)) +
      geom_line(color = "#2980b9", linewidth = 1.1) +
      labs(title = "AKT Activity (PI3K/PTEN pathway)",
           x = "Months", y = "AKT Activity (normalized)") +
      theme_qsp
  })

  # ----- TAB 6: Bone Metastasis -----
  bone_sim <- eventReactive(input$run_bone, {
    run_sim(
      leup = input$use_leup_b,
      den  = input$use_den_b,
      duration = input$sim_duration,
      pars = list(
        kRANKL   = input$rankl_b,
        kOPG     = input$opg_b,
        k_bonehom = input$k_bonehom_b
      )
    )
  }, ignoreNULL = FALSE)

  output$bmd_plot <- renderPlot({
    ggplot(bone_sim(), aes(x = time / 30.4, y = BMD)) +
      geom_line(color = "#1abc9c", linewidth = 1.1) +
      geom_hline(yintercept = 0.9, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = 0.75, linetype = "dashed", color = "red") +
      labs(title = "Bone Mineral Density",
           subtitle = "Orange: Osteopenia (-1 T-score); Red: Osteoporosis (-2.5)",
           x = "Months", y = "BMD (normalized)") +
      theme_qsp
  })

  output$bone_mets_plot <- renderPlot({
    ggplot(bone_sim(), aes(x = time / 30.4, y = BoneMets)) +
      geom_line(color = "#e67e22", linewidth = 1.1) +
      labs(title = "Bone Metastasis Burden",
           x = "Months", y = "Bone Mets (normalized)") +
      theme_qsp
  })

  output$oc_ob_plot <- renderPlot({
    bone_sim() %>%
      select(time, OC, OB) %>%
      pivot_longer(-time, names_to = "Cell", values_to = "Count") %>%
      mutate(Cell = recode(Cell,
        OC = "Osteoclasts (bone resorption)",
        OB = "Osteoblasts (bone formation)")) %>%
      ggplot(aes(x = time / 30.4, y = Count, color = Cell)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("#e74c3c", "#3498db")) +
      labs(title = "Osteoclast/Osteoblast Balance",
           x = "Months", y = "Normalized Count", color = "") +
      theme_qsp
  })

  output$bone_markers <- renderPlot({
    df <- bone_sim() %>%
      mutate(
        ALP  = OB * 100,  # proxy for ALP
        CTx  = OC * 0.5   # proxy for CTx
      ) %>%
      select(time, ALP, CTx) %>%
      pivot_longer(-time, names_to = "Marker", values_to = "Value")

    ggplot(df, aes(x = time / 30.4, y = Value, color = Marker)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c(ALP = "#27ae60", CTx = "#e74c3c")) +
      labs(title = "Bone Turnover Markers",
           subtitle = "ALP (formation proxy) vs CTx (resorption proxy)",
           x = "Months", y = "Marker Level", color = "") +
      theme_qsp
  })

  # ----- TAB 7: Scenario Comparison -----
  cmp_results <- eventReactive(input$run_cmp, {
    scenario_map <- list(
      s0 = list(leup=F, deg=F, enz=F, abi=F, doc=F, ola=F, den=F, name="Untreated"),
      s1 = list(leup=T, deg=F, enz=F, abi=F, doc=F, ola=F, den=F, name="ADT Alone"),
      s2 = list(leup=T, deg=F, enz=T, abi=F, doc=F, ola=F, den=F, name="ADT+Enzalutamide"),
      s3 = list(leup=T, deg=F, enz=F, abi=T, doc=F, ola=F, den=F, name="ADT+Abiraterone"),
      s4 = list(leup=T, deg=F, enz=F, abi=F, doc=T, ola=F, den=F, name="ADT+Docetaxel ×6"),
      s5 = list(leup=T, deg=F, enz=F, abi=F, doc=F, ola=T, den=F, name="ADT+Olaparib(HRR)"),
      s6 = list(leup=T, deg=F, enz=T, abi=F, doc=T, ola=F, den=F, name="Seq ADT→ARPI→Doc")
    )

    selected  <- input$scenarios_cmp
    hrr_def_v <- ifelse("s5" %in% selected, 1.0, 0.0)
    pars_base <- list(PTEN_loss = input$cmp_pten)

    results <- lapply(selected, function(s) {
      sc <- scenario_map[[s]]
      pars <- if (s == "s5") c(pars_base, list(HRR_def = 1.0)) else pars_base
      run_sim(leup=sc$leup, deg=sc$deg, enz=sc$enz, abi=sc$abi,
              doc=sc$doc, ola=sc$ola, den=sc$den,
              duration=input$cmp_dur, pars=pars) %>%
        mutate(Scenario = sc$name)
    })
    bind_rows(results)
  }, ignoreNULL = FALSE)

  output$cmp_plot <- renderPlot({
    df <- cmp_results()
    pal <- c("#e74c3c","#3498db","#2ecc71","#f39c12",
             "#9b59b6","#1abc9c","#e67e22")[seq_along(unique(df$Scenario))]

    y_var <- switch(input$cmp_endpoint,
      "PSA (ng/mL)"       = "PSA",
      "Testosterone (ng/dL)" = "T",
      "Tumor Burden"      = "TC_p",
      "BMD"               = "BMD",
      "CRPC Fraction"     = "CRPC_frac",
      "Bone Metastasis"   = "BoneMets"
    )

    p <- ggplot(df, aes(x = time / 30.4, y = .data[[y_var]],
                         color = Scenario, group = Scenario)) +
      geom_line(linewidth = 1.0) +
      scale_color_manual(values = pal) +
      labs(title = paste("Scenario Comparison:", input$cmp_endpoint),
           x = "Time (months)", y = input$cmp_endpoint, color = "Scenario") +
      theme_qsp

    if (y_var == "PSA") p <- p + scale_y_log10()
    if (y_var == "T")   p <- p + geom_hline(yintercept = 0.06, linetype="dashed")
    p
  })

  output$endpoints_table <- DT::renderDataTable({
    df <- cmp_results()
    df %>%
      group_by(Scenario) %>%
      summarise(
        `PSA Nadir (ng/mL)` = round(min(PSA, na.rm=TRUE), 2),
        `T Nadir (nmol/L)`  = round(min(T, na.rm=TRUE), 3),
        `Castrate (<1.73)`  = any(T < 1.73, na.rm=TRUE),
        `PSA50 Response`    = any(PSA < PSA[time==0][1]*0.5, na.rm=TRUE),
        `TC Final`          = round((TC_p+TC_q)[time==max(time)][1], 2),
        `BMD Final`         = round(BMD[time==max(time)][1], 3),
        `CRPC %`            = round(CRPC_frac[time==max(time)][1]*100, 1),
        `ARv7 %`            = round(ARv7_frac[time==max(time)][1]*100, 1),
        .groups = "drop"
      ) %>%
      DT::datatable(options = list(scrollX = TRUE, pageLength = 10))
  })

  # ----- TAB 8: Sensitivity -----
  sa_results <- eventReactive(input$run_sa, {
    n_vals <- input$sa_n

    x_range <- switch(input$sa_x_param,
      pten     = seq(0, 1,    length.out = n_vals),
      k_crpc   = seq(0.001, 0.01, length.out = n_vals),
      ar_ec50  = seq(0.5, 10,  length.out = n_vals),
      k_prolif = seq(0.02, 0.15, length.out = n_vals),
      k_bonehom = seq(0.005, 0.1, length.out = n_vals)
    )

    color_vals <- switch(input$sa_color_param,
      hrr     = c(0, 1),
      arv7    = c(0, 0.2),
      gleason = c(0.04, 0.06, 0.10)
    )

    grid <- expand.grid(x_val = x_range, color_val = color_vals)

    sa_out <- mapply(function(xv, cv) {
      par_list <- list()
      par_list[[switch(input$sa_x_param,
        pten="PTEN_loss", k_crpc="k_CRPC", ar_ec50="Enz_EC50",
        k_prolif="k_prolif", k_bonehom="k_bonehom")]] <- xv

      par_list[[switch(input$sa_color_param,
        hrr="HRR_def", arv7="ARv7_frac", gleason="k_prolif")]] <- cv

      leup_v <- input$sa_treatment %in% c("adt","adtenz")
      enz_v  <- input$sa_treatment == "adtenz"
      doc_v  <- input$sa_treatment == "doc"

      out <- run_sim(leup=leup_v, enz=enz_v, doc=doc_v,
                     duration = 730, pars = par_list)

      endpoint_val <- switch(input$sa_endpoint,
        "PSA at 3 months"        = out$PSA[out$time == 90][1],
        "Tumor Burden at 1 year" = (out$TC_p + out$TC_q)[out$time == 365][1],
        "CRPC Fraction at 2 years"= out$CRPC_frac[out$time == 730][1],
        "BMD at 2 years"         = out$BMD[out$time == 730][1]
      )

      data.frame(x_val = xv, color_val = as.factor(round(cv, 3)),
                 endpoint = endpoint_val)
    }, grid$x_val, grid$color_val, SIMPLIFY = FALSE) %>%
      bind_rows()

    sa_out
  }, ignoreNULL = FALSE)

  output$sa_plot <- renderPlot({
    df <- sa_results()
    ggplot(df, aes(x = x_val, y = endpoint, color = color_val,
                   group = color_val)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2.5) +
      scale_color_brewer(palette = "Set1") +
      labs(title = paste("Sensitivity:", input$sa_endpoint, "vs", input$sa_x_param),
           x = input$sa_x_param, y = input$sa_endpoint,
           color = input$sa_color_param) +
      theme_qsp
  })

  output$tornado_plot <- renderPlot({
    df <- sa_results()
    tornado_df <- df %>%
      group_by(color_val) %>%
      summarise(
        min_val = min(endpoint, na.rm=TRUE),
        max_val = max(endpoint, na.rm=TRUE),
        range   = max_val - min_val,
        .groups = "drop"
      ) %>%
      arrange(desc(range))

    ggplot(tornado_df, aes(y = reorder(color_val, range))) +
      geom_segment(aes(x = min_val, xend = max_val,
                       yend = reorder(color_val, range)),
                   linewidth = 6, color = "#3498db", alpha = 0.7) +
      geom_vline(xintercept = median(sa_results()$endpoint, na.rm=TRUE),
                 color = "red", linetype = "dashed") +
      labs(title = "Tornado Plot: Parameter Impact Range",
           x = input$sa_endpoint, y = input$sa_color_param) +
      theme_qsp
  })
}

# ==============================================================================
# LAUNCH
# ==============================================================================
shinyApp(ui = ui, server = server)
