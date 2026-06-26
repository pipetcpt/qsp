## ============================================================
## FMF QSP Shiny Dashboard
## Familial Mediterranean Fever — Interactive Simulator
## Tabs: Patient Profile | Drug PK | Inflammasome | Attack Sim |
##       Clinical Endpoints | Scenario Comparison | Amyloidosis | Sensitivity
## ============================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(DT)

## ─── mrgsolve model (embedded) ─────────────────────────────────────────────

fmf_model_code <- '
$PARAM
ka_col=1.2, F_col=0.45, CL_col=18.0, V1_col=120.0,
Q_col=60.0, V2_col=480.0, k_leu_on=2.5, k_leu_off=0.08,
ka_ana=0.40, F_ana=0.95, CL_ana=1.8, V_ana=8.5,
ka_cana=0.012, F_cana=0.70, CL_cana=0.18, V1_cana=4.5,
Q_cana=0.4, V2_cana=3.0,
k_RhoA_ss=1.0, k_RhoA_deg=1.0,
k_phos_basal=0.5, k_phos_mut=0.05, k_dephos=0.3,
k_ASC_form=0.8, k_ASC_deg=0.4,
k_Casp1_act=1.2, k_Casp1_deg=0.6,
MEFV_severity=1.0,
k_IL1b_pro_prod=0.3, k_IL1b_pro_deg=0.15,
k_IL1b_mat_form=2.0, Km_IL1b_mat=0.5, k_IL1b_mat_deg=0.8, IL1b_0=5.0,
k_IL18_form=1.5, k_IL18_deg=0.3, IL18_0=200.0,
k_SAA_base=0.05, k_SAA_IL1b=1.5, k_SAA_deg=0.04, SAA_0=5.0,
k_CRP_base=0.01, k_CRP_IL1b=0.3, k_CRP_deg=0.03, CRP_0=3.0,
k_Neu_prod=50.0, k_Neu_circ_deg=0.03, k_Neu_migr=0.05, k_Neu_tis_deg=0.08,
Neu_circ_0=4000.0, Neu_tis_0=500.0,
Att_trig_thresh=2.0, k_Att_rise=0.3, k_Att_decay=0.1,
IC50_col_neu=8.0, Emax_col_neu=0.85,
IC50_col_PYRIN=5.0, Emax_col_PYRIN=0.60,
IC50_ana=50.0, Emax_ana=0.95,
IC50_cana=10.0, Emax_cana=0.98,
k_AA_dep=0.001, k_AA_deg=0.0002,
k_eGFR_dec=0.005, eGFR_0=90.0,
USE_COL=0, USE_ANA=0, USE_CANA=0

$CMT GUT_COL CENT_COL PERI_COL LEU_COL
     SC_ANA CENT_ANA
     SC_CANA CENT_CANA PERI_CANA
     RhoA Pyrin_p ASC Casp1
     IL1b_pro IL1b_mat IL18 SAA CRP
     Neu_circ Neu_tis Att_sev
     AA_dep eGFR

$INIT GUT_COL=0,CENT_COL=0,PERI_COL=0,LEU_COL=0,
      SC_ANA=0,CENT_ANA=0,
      SC_CANA=0,CENT_CANA=0,PERI_CANA=0,
      RhoA=1.0,Pyrin_p=0.5,ASC=0.2,Casp1=0.1,
      IL1b_pro=30.0,IL1b_mat=5.0,IL18=200.0,SAA=5.0,CRP=3.0,
      Neu_circ=4000.0,Neu_tis=500.0,Att_sev=0.0,
      AA_dep=0.0,eGFR=90.0

$ODE
double Cp_col=CENT_COL/V1_col;
double Cl_col_conc=LEU_COL;
dxdt_GUT_COL=-ka_col*GUT_COL;
dxdt_CENT_COL=ka_col*F_col*GUT_COL-(CL_col/V1_col)*CENT_COL-(Q_col/V1_col)*CENT_COL+(Q_col/V2_col)*PERI_COL;
dxdt_PERI_COL=(Q_col/V1_col)*CENT_COL-(Q_col/V2_col)*PERI_COL;
dxdt_LEU_COL=k_leu_on*Cp_col-k_leu_off*LEU_COL;
double Cp_ana=CENT_ANA/V_ana;
dxdt_SC_ANA=-ka_ana*SC_ANA;
dxdt_CENT_ANA=ka_ana*F_ana*SC_ANA-(CL_ana/V_ana)*CENT_ANA;
double Cp_cana=CENT_CANA/V1_cana;
dxdt_SC_CANA=-ka_cana*SC_CANA;
dxdt_CENT_CANA=ka_cana*F_cana*SC_CANA-(CL_cana/V1_cana)*CENT_CANA-(Q_cana/V1_cana)*CENT_CANA+(Q_cana/V2_cana)*PERI_CANA;
dxdt_PERI_CANA=(Q_cana/V1_cana)*CENT_CANA-(Q_cana/V2_cana)*PERI_CANA;
double E_col_neu=(USE_COL*Emax_col_neu*Cl_col_conc)/(IC50_col_neu+Cl_col_conc+1e-10);
double E_col_PYRIN=(USE_COL*Emax_col_PYRIN*Cl_col_conc)/(IC50_col_PYRIN+Cl_col_conc+1e-10);
double E_ana=(USE_ANA*Emax_ana*Cp_ana)/(IC50_ana+Cp_ana+1e-10);
double E_cana=(USE_CANA*Emax_cana*Cp_cana)/(IC50_cana+Cp_cana+1e-10);
double IL1b_block=1.0-fmax(E_ana,E_cana);
double k_phos_eff=k_phos_basal*(1.0-MEFV_severity*0.90)+k_phos_mut*MEFV_severity*0.90;
dxdt_RhoA=k_RhoA_ss-k_RhoA_deg*RhoA;
double k_phos_drug=k_phos_eff*(1.0+0.3*E_col_PYRIN);
dxdt_Pyrin_p=k_phos_drug*(1.0-Pyrin_p)*RhoA-k_dephos*Pyrin_p;
double active_pyrin=fmax(0.0,1.0-Pyrin_p);
dxdt_ASC=k_ASC_form*active_pyrin*(1.0-E_col_PYRIN)-k_ASC_deg*ASC;
dxdt_Casp1=k_Casp1_act*ASC-k_Casp1_deg*Casp1;
double rate_IL1b_mat=k_IL1b_mat_form*Casp1*IL1b_pro/(Km_IL1b_mat+IL1b_pro);
dxdt_IL1b_pro=k_IL1b_pro_prod-k_IL1b_pro_deg*IL1b_pro-rate_IL1b_mat;
dxdt_IL1b_mat=rate_IL1b_mat-k_IL1b_mat_deg*IL1b_mat*IL1b_block;
dxdt_IL18=k_IL18_form*Casp1-k_IL18_deg*IL18;
dxdt_SAA=k_SAA_base+k_SAA_IL1b*IL1b_mat*IL1b_block-k_SAA_deg*SAA;
dxdt_CRP=k_CRP_base+k_CRP_IL1b*IL1b_mat*IL1b_block-k_CRP_deg*CRP;
double Neu_migr_eff=k_Neu_migr*(1.0-E_col_neu)*(IL1b_mat/IL1b_0);
dxdt_Neu_circ=k_Neu_prod-k_Neu_circ_deg*Neu_circ-Neu_migr_eff*Neu_circ;
dxdt_Neu_tis=Neu_migr_eff*Neu_circ-k_Neu_tis_deg*Neu_tis;
double att_signal=fmax(0.0,IL1b_mat/IL1b_0-Att_trig_thresh);
dxdt_Att_sev=k_Att_rise*att_signal*(1.0-Att_sev/10.0)-k_Att_decay*Att_sev;
double SAA_excess=fmax(0.0,SAA-10.0);
dxdt_AA_dep=k_AA_dep*SAA_excess-k_AA_deg*AA_dep;
dxdt_eGFR=-k_eGFR_dec*AA_dep;

$CAPTURE Cp_col Cl_col_conc Cp_ana Cp_cana
         IL1b_mat IL18 SAA CRP Neu_circ Neu_tis
         Att_sev AA_dep eGFR active_pyrin ASC Casp1

$TABLE
double AIDAI_score=fmin(10.0,Att_sev);
double attack_flag=(Att_sev>2.0)?1.0:0.0;
'

## Load model at startup
MOD <- mcode("fmf_shiny", fmf_model_code, quiet = TRUE)

## ─── Helper functions ──────────────────────────────────────────────────────

run_sim <- function(mod, params, evs = NULL, end = 8760, delta = 4) {
  m <- mod %>% param(params)
  if (!is.null(evs)) {
    out <- m %>% ev(evs) %>% mrgsim(end = end, delta = delta)
  } else {
    out <- m %>% mrgsim(end = end, delta = delta)
  }
  as.data.frame(out) %>% mutate(time_wk = time / 168)
}

make_events <- function(cmt, amt, interval_h, end) {
  times <- seq(0, end, by = interval_h)
  ev(amt = amt, cmt = cmt, time = times)
}

theme_fmf <- function() {
  theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 13))
}

## ─── UI ─────────────────────────────────────────────────────────────────────

ui <- fluidPage(
  titlePanel(
    div(style = "background: linear-gradient(135deg,#1a3c5e,#2e86ab); color:white; padding:15px 20px; border-radius:8px; margin-bottom:15px;",
        h2("Familial Mediterranean Fever (FMF) — QSP Simulator",
           style = "margin:0; font-size:22px;"),
        p("PYRIN Inflammasome · IL-1β · Colchicine · IL-1 Inhibitors",
          style = "margin:4px 0 0 0; font-size:13px; opacity:0.85;"))
  ),

  tabsetPanel(id = "tabs",

    ## ── Tab 1: Patient Profile ──────────────────────────────────────────────
    tabPanel("Patient Profile",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Patient Characteristics"),
          selectInput("genotype", "MEFV Genotype",
            choices = c("Wild-type (no FMF)"    = "0.0",
                        "E148Q (mild)"           = "0.3",
                        "V726A (moderate)"       = "0.6",
                        "M680I (mod-severe)"     = "0.8",
                        "M694V/M694V (severe)"   = "1.0"),
            selected = "1.0"),
          sliderInput("age", "Age (years)", 5, 70, 30),
          selectInput("sex", "Sex", c("Male", "Female")),
          sliderInput("bmi", "BMI", 16, 45, 24, step = 0.5),
          selectInput("ethnicity", "Ethnicity",
            choices = c("Turkish", "Armenian", "Arabic", "Sephardic Jewish", "Other")),
          hr(),
          h4("Baseline Labs"),
          sliderInput("base_CRP", "Baseline CRP (mg/L)", 1, 50, 3),
          sliderInput("base_SAA", "Baseline SAA (mg/L)", 1, 100, 5),
          sliderInput("base_eGFR","Baseline eGFR (mL/min/1.73m²)", 30, 120, 90),
          hr(),
          h4("Treatment"),
          selectInput("drug_profile", "Current Drug",
            choices = c("None", "Colchicine 0.5 mg BID",
                        "Colchicine 1.0 mg QD", "Anakinra 100 mg QD", "Canakinumab 150 mg Q8W")),
          actionButton("run_profile", "Update Profile", class = "btn-primary btn-sm")
        ),
        mainPanel(width = 9,
          fluidRow(
            column(4,
              wellPanel(
                h4("FMF Classification (Eurofever/PRINTO 2019)"),
                uiOutput("classification_result")
              )
            ),
            column(4,
              wellPanel(
                h4("Predicted 1-year Outcomes"),
                tableOutput("profile_outcomes")
              )
            ),
            column(4,
              wellPanel(
                h4("Amyloidosis Risk"),
                uiOutput("amyloid_risk")
              )
            )
          ),
          plotOutput("profile_dynamics", height = "400px")
        )
      )
    ),

    ## ── Tab 2: Drug PK ─────────────────────────────────────────────────────
    tabPanel("Drug PK",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Colchicine"),
          selectInput("col_dose", "Dose",
            c("0.5 mg" = "500000", "1.0 mg" = "1000000", "1.5 mg" = "1500000")),
          selectInput("col_freq", "Frequency",
            c("QD (every 24h)" = "24", "BID (every 12h)" = "12",
              "TID (every 8h)" = "8")),
          numericInput("col_days", "Duration (days)", 14, 1, 90),
          hr(),
          h4("Anakinra"),
          sliderInput("ana_dose", "SC Dose (mg)", 50, 200, 100, step = 10),
          numericInput("ana_days", "Duration (days)", 30, 1, 365),
          hr(),
          h4("Canakinumab"),
          sliderInput("cana_dose", "SC Dose (mg)", 150, 300, 150, step = 50),
          numericInput("cana_weeks", "Number of doses", 3, 1, 12),
          hr(),
          actionButton("run_pk", "Simulate PK", class = "btn-primary btn-sm")
        ),
        mainPanel(width = 9,
          tabsetPanel(
            tabPanel("Colchicine PK",
              plotOutput("pk_col_plot", height = "400px"),
              br(),
              tableOutput("pk_col_table")
            ),
            tabPanel("Anakinra PK",
              plotOutput("pk_ana_plot", height = "400px")
            ),
            tabPanel("Canakinumab PK",
              plotOutput("pk_cana_plot", height = "400px")
            )
          )
        )
      )
    ),

    ## ── Tab 3: Inflammasome Dynamics ──────────────────────────────────────
    tabPanel("Inflammasome Dynamics",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("PYRIN Pathway Parameters"),
          sliderInput("mefv_sev", "MEFV Severity (0=WT, 1=M694V)", 0, 1, 1, step = 0.1),
          sliderInput("k_ASC_form", "ASC Formation Rate", 0.1, 2.0, 0.8, step = 0.1),
          sliderInput("k_Casp1_act", "Caspase-1 Activation Rate", 0.1, 3.0, 1.2, step = 0.1),
          sliderInput("k_phos_basal", "Basal Phosphorylation Rate", 0.1, 1.0, 0.5, step = 0.05),
          hr(),
          h4("Drug Effect on PYRIN"),
          checkboxInput("use_col_inflam", "Use Colchicine", TRUE),
          sliderInput("col_leu_conc", "Colchicine Leukocyte Conc (ng/mL)", 0, 50, 10),
          hr(),
          sliderInput("sim_days_inflam", "Simulation (days)", 3, 60, 14),
          actionButton("run_inflam", "Simulate", class = "btn-primary btn-sm")
        ),
        mainPanel(width = 9,
          plotOutput("inflam_plot", height = "600px"),
          br(),
          h4("Inflammasome Component Summary"),
          tableOutput("inflam_table")
        )
      )
    ),

    ## ── Tab 4: Attack Simulation ───────────────────────────────────────────
    tabPanel("Attack Simulation",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Attack Trigger"),
          sliderInput("att_thresh", "IL-1β Threshold (fold above baseline)", 1.0, 5.0, 2.0, step = 0.1),
          sliderInput("att_rise", "Attack Rise Rate", 0.05, 1.0, 0.3, step = 0.05),
          sliderInput("att_decay", "Attack Decay Rate", 0.02, 0.5, 0.1, step = 0.01),
          hr(),
          h4("Treatment"),
          checkboxInput("att_col", "Colchicine 1 mg QD", TRUE),
          checkboxInput("att_ana", "Anakinra 100 mg QD", FALSE),
          checkboxInput("att_cana","Canakinumab 150 mg Q8W", FALSE),
          hr(),
          sliderInput("att_mefv", "MEFV Severity", 0, 1, 1, step = 0.1),
          sliderInput("sim_weeks_att", "Simulation (weeks)", 4, 104, 52),
          actionButton("run_att", "Simulate Attacks", class = "btn-primary btn-sm")
        ),
        mainPanel(width = 9,
          plotOutput("att_plot", height = "500px"),
          br(),
          fluidRow(
            column(6, h4("Attack Statistics"), tableOutput("att_stats")),
            column(6, h4("IL-1β Timeline"), plotOutput("att_IL1b", height = "250px"))
          )
        )
      )
    ),

    ## ── Tab 5: Clinical Endpoints ──────────────────────────────────────────
    tabPanel("Clinical Endpoints",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Simulation Settings"),
          selectInput("ep_drug", "Treatment",
            c("No Treatment", "Colchicine 0.5 mg BID",
              "Colchicine 1.0 mg QD", "Anakinra 100 mg QD",
              "Canakinumab 150 mg Q8W")),
          sliderInput("ep_mefv", "MEFV Severity", 0, 1, 1, step = 0.1),
          sliderInput("ep_weeks", "Duration (weeks)", 12, 104, 52),
          actionButton("run_ep", "Calculate Endpoints", class = "btn-primary btn-sm"),
          hr(),
          h4("Clinical Targets"),
          p(strong("Colchicine response:"), "≥50% attack reduction"),
          p(strong("IL-1 inhibitor response:"), "Attack-free + SAA <10 mg/L"),
          p(strong("Complete response:"), "SAA <10 mg/L + no attacks"),
          p(strong("Amyloidosis prevention:"), "SAA <10 mg/L sustained")
        ),
        mainPanel(width = 9,
          fluidRow(
            column(3, wellPanel(uiOutput("ep_box1"))),
            column(3, wellPanel(uiOutput("ep_box2"))),
            column(3, wellPanel(uiOutput("ep_box3"))),
            column(3, wellPanel(uiOutput("ep_box4")))
          ),
          plotOutput("ep_plot", height = "400px"),
          br(),
          h4("AIDAI Score Over Time"),
          plotOutput("aidai_plot", height = "250px")
        )
      )
    ),

    ## ── Tab 6: Scenario Comparison ─────────────────────────────────────────
    tabPanel("Scenario Comparison",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Scenarios to Compare"),
          checkboxGroupInput("sc_drugs",
            "Select treatments:",
            choices = c("No Treatment",
                        "Colchicine 0.5 mg BID",
                        "Colchicine 1.0 mg QD",
                        "Anakinra 100 mg QD",
                        "Canakinumab 150 mg Q8W"),
            selected = c("No Treatment", "Colchicine 0.5 mg BID",
                         "Canakinumab 150 mg Q8W")),
          sliderInput("sc_mefv", "MEFV Severity", 0, 1, 1, step = 0.1),
          sliderInput("sc_weeks", "Duration (weeks)", 12, 104, 52),
          actionButton("run_sc", "Compare Scenarios", class = "btn-primary btn-sm"),
          hr(),
          selectInput("sc_outcome", "Primary Endpoint",
            c("IL-1β (pg/mL)" = "IL1b_mat",
              "SAA (mg/L)"    = "SAA",
              "CRP (mg/L)"    = "CRP",
              "Attack Severity" = "Att_sev",
              "Tissue Neutrophils" = "Neu_tis",
              "eGFR"          = "eGFR"))
        ),
        mainPanel(width = 9,
          plotOutput("sc_plot_main", height = "350px"),
          br(),
          h4("Comparative Outcomes Table"),
          DTOutput("sc_table"),
          br(),
          plotOutput("sc_radar", height = "300px")
        )
      )
    ),

    ## ── Tab 7: Amyloidosis Risk ────────────────────────────────────────────
    tabPanel("Amyloidosis Risk",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Long-term Projection"),
          sliderInput("amy_years", "Projection (years)", 1, 20, 10),
          sliderInput("amy_mefv", "MEFV Severity", 0, 1, 1, step = 0.1),
          hr(),
          h4("SAA Treatment Target"),
          sliderInput("amy_SAA_target", "Target SAA (mg/L)", 1, 30, 10),
          hr(),
          selectInput("amy_drug", "Treatment Arm",
            c("No Treatment", "Colchicine 0.5 mg BID",
              "Colchicine 1.0 mg QD", "Anakinra 100 mg QD",
              "Canakinumab 150 mg Q8W")),
          actionButton("run_amy", "Project Risk", class = "btn-primary btn-sm")
        ),
        mainPanel(width = 9,
          fluidRow(
            column(6, plotOutput("amy_SAA_plot", height = "300px")),
            column(6, plotOutput("amy_dep_plot", height = "300px"))
          ),
          plotOutput("amy_eGFR_plot", height = "300px"),
          br(),
          h4("Amyloidosis Risk Summary"),
          tableOutput("amy_table")
        )
      )
    ),

    ## ── Tab 8: Sensitivity Analysis ────────────────────────────────────────
    tabPanel("Sensitivity Analysis",
      sidebarLayout(
        sidebarPanel(width = 3,
          h4("Parameter to Vary"),
          selectInput("sens_param", "Parameter",
            c("Colchicine IC50 (Neutrophil)" = "IC50_col_neu",
              "Colchicine IC50 (PYRIN)"      = "IC50_col_PYRIN",
              "Anakinra IC50"                = "IC50_ana",
              "Canakinumab IC50"             = "IC50_cana",
              "IL-1β Maturation Rate"        = "k_IL1b_mat_form",
              "SAA Production Rate"          = "k_SAA_IL1b",
              "ASC Formation Rate"           = "k_ASC_form")),
          sliderInput("sens_min", "Min value", 0.01, 100, 1),
          sliderInput("sens_max", "Max value", 0.1, 500, 50),
          numericInput("sens_n", "Number of values", 5, 3, 15),
          hr(),
          selectInput("sens_outcome", "Endpoint",
            c("Mean IL-1β" = "IL1b_mat",
              "Mean SAA"   = "SAA",
              "Tissue Neutrophils" = "Neu_tis",
              "Attack Severity" = "Att_sev")),
          sliderInput("sens_mefv", "MEFV Severity", 0, 1, 1, step = 0.1),
          hr(),
          actionButton("run_sens", "Run Sensitivity", class = "btn-primary btn-sm")
        ),
        mainPanel(width = 9,
          plotOutput("sens_tornado", height = "400px"),
          br(),
          plotOutput("sens_time", height = "350px"),
          br(),
          h4("Sensitivity Results"),
          tableOutput("sens_table")
        )
      )
    )
  )
)

## ─── Server ─────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  ## helper: build event list from drug selections
  build_events <- function(drug, end_h, col_amt = 1000000) {
    evs <- switch(drug,
      "Colchicine 0.5 mg BID" = make_ev("GUT_COL", 500000, 12, end_h),
      "Colchicine 1.0 mg QD"  = make_ev("GUT_COL", 1000000, 24, end_h),
      "Anakinra 100 mg QD"    = make_ev("SC_ANA", 100000, 24, end_h),
      "Canakinumab 150 mg Q8W"= make_ev("SC_CANA", 150, 56*24, end_h),
      NULL
    )
    evs
  }

  make_ev <- function(cmt, amt, intv, end) {
    ev(amt = amt, cmt = cmt, time = seq(0, end, by = intv))
  }

  get_params <- function(drug, mefv = 1.0) {
    p <- list(MEFV_severity = mefv,
              USE_COL = 0, USE_ANA = 0, USE_CANA = 0)
    if (grepl("Colchicine", drug)) p$USE_COL <- 1
    if (grepl("Anakinra",   drug)) p$USE_ANA  <- 1
    if (grepl("Canakinumab",drug)) p$USE_CANA <- 1
    p
  }

  ## ── Tab 1: Patient Profile ────────────────────────────────────────────────
  profile_data <- eventReactive(input$run_profile, {
    mefv <- as.numeric(input$genotype)
    end  <- 8760
    drug <- input$drug_profile
    params <- get_params(drug, mefv)
    params$eGFR_0 <- input$base_eGFR
    params$SAA_0  <- input$base_SAA
    params$CRP_0  <- input$base_CRP
    evs <- build_events(drug, end)
    run_sim(MOD, params, evs, end, delta = 6)
  }, ignoreNULL = FALSE)

  output$profile_dynamics <- renderPlot({
    df <- profile_data()
    p1 <- ggplot(df, aes(time_wk, IL1b_mat)) + geom_line(color="#e74c3c", linewidth=0.9) +
          labs(x="Week",y="IL-1β (pg/mL)",title="Mature IL-1β") + theme_fmf()
    p2 <- ggplot(df, aes(time_wk, SAA)) + geom_line(color="#e67e22", linewidth=0.9) +
          geom_hline(yintercept=10,linetype="dashed",color="red") +
          labs(x="Week",y="SAA (mg/L)",title="Serum Amyloid A") + theme_fmf()
    p3 <- ggplot(df, aes(time_wk, Att_sev)) + geom_line(color="#8e44ad", linewidth=0.9) +
          labs(x="Week",y="Score",title="Attack Severity") + theme_fmf()
    p4 <- ggplot(df, aes(time_wk, eGFR)) + geom_line(color="#27ae60", linewidth=0.9) +
          labs(x="Week",y="eGFR",title="eGFR (mL/min/1.73m²)") + theme_fmf()
    (p1 + p2) / (p3 + p4)
  })

  output$classification_result <- renderUI({
    mefv <- as.numeric(input$genotype)
    category <- if (mefv >= 0.8) "Definite FMF" else if (mefv >= 0.5) "Probable FMF" else if (mefv > 0) "Possible FMF" else "Not FMF (WT)"
    color <- if (mefv >= 0.8) "#e74c3c" else if (mefv >= 0.5) "#e67e22" else if (mefv > 0) "#f39c12" else "#27ae60"
    tags$div(
      tags$span(category, style = paste0("font-size:16px; font-weight:bold; color:", color)),
      br(), br(),
      tags$small(paste0("MEFV severity score: ", mefv)),
      br(),
      tags$small(paste0("Ethnicity: ", input$ethnicity)),
      br(),
      tags$small(paste0("Age at onset: ", max(5, input$age - 15), "y"))
    )
  })

  output$profile_outcomes <- renderTable({
    df <- profile_data()
    data.frame(
      Metric = c("Mean IL-1β (pg/mL)", "Mean SAA (mg/L)", "Mean CRP (mg/L)",
                 "Attacks/Year", "Final eGFR"),
      Value  = c(
        round(mean(df$IL1b_mat), 1),
        round(mean(df$SAA), 1),
        round(mean(df$CRP), 1),
        round(sum(diff(c(0, df$attack_flag)) == 1), 0),
        round(last(df$eGFR), 1)
      )
    )
  })

  output$amyloid_risk <- renderUI({
    df <- profile_data()
    mean_SAA <- mean(df$SAA)
    risk_color <- if (mean_SAA > 30) "#e74c3c" else if (mean_SAA > 10) "#e67e22" else "#27ae60"
    risk_label <- if (mean_SAA > 30) "HIGH" else if (mean_SAA > 10) "MODERATE" else "LOW"
    tags$div(
      tags$span(paste0("Risk: ", risk_label), style = paste0("font-size:16px;font-weight:bold;color:", risk_color)),
      br(), br(),
      tags$small(paste0("Mean SAA: ", round(mean_SAA, 1), " mg/L")),
      br(),
      tags$small(if (mean_SAA > 10) "Target: SAA < 10 mg/L" else "SAA at target")
    )
  })

  ## ── Tab 2: Drug PK ───────────────────────────────────────────────────────
  pk_data <- eventReactive(input$run_pk, {
    end_h <- as.numeric(input$col_days) * 24
    end_ana <- as.numeric(input$ana_days) * 24
    end_cana <- as.numeric(input$cana_weeks) * 56 * 24

    ev_col  <- make_ev("GUT_COL", as.numeric(input$col_dose),
                       as.numeric(input$col_freq), end_h)
    ev_ana  <- make_ev("SC_ANA", input$ana_dose * 1000, 24, end_ana)
    ev_cana <- make_ev("SC_CANA", input$cana_dose, 56*24, end_cana)

    col_df  <- MOD %>% param(USE_COL=1,USE_ANA=0,USE_CANA=0,MEFV_severity=1) %>%
               ev(ev_col) %>% mrgsim(end=end_h, delta=0.5) %>% as.data.frame()
    ana_df  <- MOD %>% param(USE_COL=0,USE_ANA=1,USE_CANA=0,MEFV_severity=1) %>%
               ev(ev_ana) %>% mrgsim(end=end_ana, delta=0.5) %>% as.data.frame()
    cana_df <- MOD %>% param(USE_COL=0,USE_ANA=0,USE_CANA=1,MEFV_severity=1) %>%
               ev(ev_cana) %>% mrgsim(end=end_cana, delta=2) %>% as.data.frame()

    list(col = col_df, ana = ana_df, cana = cana_df)
  }, ignoreNULL = FALSE)

  output$pk_col_plot <- renderPlot({
    df <- pk_data()$col
    df %>%
      select(time, Cp_col, Cl_col_conc) %>%
      pivot_longer(-time) %>%
      mutate(name = ifelse(name=="Cp_col","Plasma","Leukocyte")) %>%
      ggplot(aes(time/24, value, color=name)) +
      geom_line(linewidth=0.9) +
      labs(x="Time (days)", y="Concentration (ng/mL)",
           title="Colchicine PK — Plasma vs Leukocyte", color="") +
      theme_fmf()
  })

  output$pk_ana_plot <- renderPlot({
    pk_data()$ana %>%
      ggplot(aes(time/24, Cp_ana)) +
      geom_line(color="#3498db", linewidth=0.9) +
      labs(x="Time (days)", y="Concentration (ng/mL)",
           title="Anakinra SC PK") + theme_fmf()
  })

  output$pk_cana_plot <- renderPlot({
    pk_data()$cana %>%
      ggplot(aes(time/24, Cp_cana)) +
      geom_line(color="#8e44ad", linewidth=0.9) +
      labs(x="Time (days)", y="Concentration (μg/mL)",
           title="Canakinumab SC PK") + theme_fmf()
  })

  output$pk_col_table <- renderTable({
    df <- pk_data()$col
    data.frame(
      Parameter = c("Cmax plasma (ng/mL)", "Tmax plasma (h)",
                    "Cmax leukocyte (ng/mL)", "Leukocyte:Plasma ratio"),
      Value = c(
        round(max(df$Cp_col), 2),
        df$time[which.max(df$Cp_col)],
        round(max(df$Cl_col_conc), 2),
        round(max(df$Cl_col_conc) / max(df$Cp_col + 0.001), 1)
      )
    )
  })

  ## ── Tab 3: Inflammasome ──────────────────────────────────────────────────
  inflam_data <- eventReactive(input$run_inflam, {
    end <- input$sim_days_inflam * 24
    params <- list(
      MEFV_severity  = input$mefv_sev,
      k_ASC_form     = input$k_ASC_form,
      k_Casp1_act    = input$k_Casp1_act,
      k_phos_basal   = input$k_phos_basal,
      USE_COL = 0, USE_ANA = 0, USE_CANA = 0
    )
    df_base <- run_sim(MOD, params, NULL, end, 1)

    if (input$use_col_inflam) {
      params$USE_COL <- 1
      # Inject a fixed leukocyte concentration via dummy event (simplified)
      ev_col <- ev(amt = input$col_leu_conc * 10 / 2.5,
                   cmt = "GUT_COL", time = seq(0, end, by = 12))
      df_drug <- run_sim(MOD, params, ev_col, end, 1)
      list(base = df_base, drug = df_drug)
    } else {
      list(base = df_base, drug = NULL)
    }
  }, ignoreNULL = FALSE)

  output$inflam_plot <- renderPlot({
    res <- inflam_data()
    df  <- res$base %>% mutate(arm = "No Drug")
    if (!is.null(res$drug)) df <- bind_rows(df, res$drug %>% mutate(arm = "Colchicine"))

    p1 <- ggplot(df, aes(time/24, active_pyrin, color=arm)) +
          geom_line(linewidth=0.9) +
          labs(x="Days",y="Level (AU)",title="Active (dephospho) PYRIN") + theme_fmf()
    p2 <- ggplot(df, aes(time/24, ASC, color=arm)) +
          geom_line(linewidth=0.9) +
          labs(x="Days",y="Level (AU)",title="ASC Speck") + theme_fmf()
    p3 <- ggplot(df, aes(time/24, Casp1, color=arm)) +
          geom_line(linewidth=0.9) +
          labs(x="Days",y="Level (AU)",title="Active Caspase-1") + theme_fmf()
    p4 <- ggplot(df, aes(time/24, IL1b_mat, color=arm)) +
          geom_line(linewidth=0.9) +
          labs(x="Days",y="pg/mL",title="Mature IL-1β") + theme_fmf()
    (p1 + p2) / (p3 + p4)
  })

  output$inflam_table <- renderTable({
    res <- inflam_data()
    df  <- res$base %>% mutate(arm = "No Drug")
    if (!is.null(res$drug)) df <- bind_rows(df, res$drug %>% mutate(arm = "Colchicine"))
    df %>% group_by(arm) %>%
      summarise(`SS active_PYRIN` = round(last(active_pyrin), 3),
                `SS ASC`          = round(last(ASC), 3),
                `SS Caspase-1`    = round(last(Casp1), 3),
                `SS IL-1β`        = round(last(IL1b_mat), 2))
  })

  ## ── Tab 4: Attack Simulation ──────────────────────────────────────────────
  att_data <- eventReactive(input$run_att, {
    end <- input$sim_weeks_att * 168
    params <- list(
      MEFV_severity  = input$att_mefv,
      Att_trig_thresh = input$att_thresh,
      k_Att_rise     = input$att_rise,
      k_Att_decay    = input$att_decay,
      USE_COL = as.integer(input$att_col),
      USE_ANA = as.integer(input$att_ana),
      USE_CANA = as.integer(input$att_cana)
    )
    evs <- NULL
    if (input$att_col)  evs <- c(evs, list(make_ev("GUT_COL", 1000000, 24, end)))
    if (input$att_ana)  evs <- c(evs, list(make_ev("SC_ANA", 100000, 24, end)))
    if (input$att_cana) evs <- c(evs, list(make_ev("SC_CANA", 150, 56*24, end)))

    all_ev <- if (!is.null(evs)) Reduce(rbind, evs) else NULL
    run_sim(MOD, params, all_ev, end, 4)
  }, ignoreNULL = FALSE)

  output$att_plot <- renderPlot({
    df <- att_data()
    ggplot(df, aes(time_wk, Att_sev)) +
      geom_area(fill="#e74c3c", alpha=0.3) +
      geom_line(color="#e74c3c", linewidth=0.8) +
      geom_hline(yintercept=2, linetype="dashed", color="darkred") +
      annotate("text", x=max(df$time_wk)*0.02, y=2.3,
               label="Attack threshold", hjust=0, size=3.5) +
      labs(x="Time (weeks)", y="Attack Severity (0-10)",
           title="FMF Attack Severity Over Time") +
      theme_fmf()
  })

  output$att_stats <- renderTable({
    df <- att_data()
    weeks <- max(df$time_wk)
    n_attacks <- sum(diff(c(0, df$attack_flag)) == 1)
    data.frame(
      Metric = c("Total attacks", "Attacks/year", "% time in attack",
                 "Mean severity", "Peak severity"),
      Value  = c(n_attacks,
                 round(n_attacks / weeks * 52, 1),
                 round(mean(df$attack_flag) * 100, 1),
                 round(mean(df$Att_sev[df$attack_flag > 0.5]), 2),
                 round(max(df$Att_sev), 2))
    )
  })

  output$att_IL1b <- renderPlot({
    df <- att_data()
    ggplot(df, aes(time_wk, IL1b_mat)) +
      geom_line(color="#e67e22") +
      labs(x="Week", y="IL-1β (pg/mL)") + theme_fmf()
  })

  ## ── Tab 5: Clinical Endpoints ─────────────────────────────────────────────
  ep_data <- eventReactive(input$run_ep, {
    end  <- input$ep_weeks * 168
    drug <- input$ep_drug
    params <- get_params(drug, input$ep_mefv)
    evs  <- build_events(drug, end)
    run_sim(MOD, params, evs, end, 4)
  }, ignoreNULL = FALSE)

  output$ep_plot <- renderPlot({
    df <- ep_data()
    p1 <- ggplot(df, aes(time_wk, SAA)) +
          geom_line(color="#e67e22",linewidth=0.8) +
          geom_hline(yintercept=10,linetype="dashed",color="red") +
          labs(x="Week",y="SAA (mg/L)",title="SAA") + theme_fmf()
    p2 <- ggplot(df, aes(time_wk, CRP)) +
          geom_line(color="#c0392b",linewidth=0.8) +
          geom_hline(yintercept=5,linetype="dashed",color="red") +
          labs(x="Week",y="CRP (mg/L)",title="CRP") + theme_fmf()
    p3 <- ggplot(df, aes(time_wk, Neu_tis)) +
          geom_line(color="#2980b9",linewidth=0.8) +
          labs(x="Week",y="AU",title="Tissue Neutrophils") + theme_fmf()
    p4 <- ggplot(df, aes(time_wk, IL18)) +
          geom_line(color="#8e44ad",linewidth=0.8) +
          labs(x="Week",y="pg/mL",title="IL-18") + theme_fmf()
    (p1 + p2) / (p3 + p4)
  })

  output$aidai_plot <- renderPlot({
    df <- ep_data()
    ggplot(df, aes(time_wk, AIDAI_score)) +
      geom_line(color="#2c3e50",linewidth=0.8) +
      geom_hline(yintercept=9,linetype="dashed",color="gray") +
      labs(x="Week",y="AIDAI Score",title="Auto-Inflammatory Disease Activity Index (AIDAI)") +
      theme_fmf()
  })

  make_ep_box <- function(label, value, unit = "", threshold = NULL, good_below = TRUE) {
    color <- if (!is.null(threshold)) {
      if (good_below) (if (value <= threshold) "#27ae60" else "#e74c3c")
      else            (if (value >= threshold) "#27ae60" else "#e74c3c")
    } else "#2c3e50"
    tags$div(
      tags$h5(label, style="color:#7f8c8d;font-size:11px;margin-bottom:4px;"),
      tags$span(paste0(round(value, 1), " ", unit),
                style = paste0("font-size:22px;font-weight:bold;color:", color))
    )
  }

  output$ep_box1 <- renderUI({
    df <- ep_data()
    make_ep_box("Mean SAA", mean(df$SAA), "mg/L", 10, TRUE)
  })
  output$ep_box2 <- renderUI({
    df <- ep_data()
    make_ep_box("Mean CRP", mean(df$CRP), "mg/L", 5, TRUE)
  })
  output$ep_box3 <- renderUI({
    df <- ep_data()
    n <- sum(diff(c(0, df$attack_flag)) == 1)
    weeks <- max(df$time_wk)
    make_ep_box("Attacks/Year", n / weeks * 52, "", 2, TRUE)
  })
  output$ep_box4 <- renderUI({
    df <- ep_data()
    make_ep_box("Final eGFR", last(df$eGFR), "mL/min", 60, FALSE)
  })

  ## ── Tab 6: Scenario Comparison ────────────────────────────────────────────
  sc_data <- eventReactive(input$run_sc, {
    end   <- input$sc_weeks * 168
    mefv  <- input$sc_mefv
    drugs <- input$sc_drugs

    lapply(drugs, function(drug) {
      params <- get_params(drug, mefv)
      evs    <- build_events(drug, end)
      run_sim(MOD, params, evs, end, 6) %>%
        mutate(drug = drug)
    }) %>% bind_rows()
  }, ignoreNULL = FALSE)

  output$sc_plot_main <- renderPlot({
    df  <- sc_data()
    var <- input$sc_outcome
    ggplot(df, aes_string("time_wk", var, color = "drug")) +
      geom_line(linewidth = 0.9) +
      labs(x = "Week", y = var,
           title = paste("Primary Endpoint:", var),
           color = "Treatment") +
      theme_fmf()
  })

  output$sc_table <- renderDT({
    df <- sc_data()
    df %>% group_by(drug) %>%
      summarise(
        `Mean IL-1β`     = round(mean(IL1b_mat), 2),
        `Mean SAA`        = round(mean(SAA), 2),
        `Mean CRP`        = round(mean(CRP), 2),
        `Attacks/yr`      = round(sum(diff(c(0, attack_flag)) == 1) / max(time_wk) * 52, 1),
        `Final eGFR`      = round(last(eGFR), 1),
        `% Attack Time`   = round(mean(attack_flag) * 100, 1)
      ) %>%
      rename(Treatment = drug) %>%
      datatable(options = list(dom = 't', pageLength = 10))
  })

  output$sc_radar <- renderPlot({
    df <- sc_data()
    summ <- df %>% group_by(drug) %>%
      summarise(IL1b = mean(IL1b_mat), SAA = mean(SAA),
                CRP = mean(CRP), AttRate = sum(diff(c(0, attack_flag))==1)/max(time_wk)*52,
                .groups="drop")
    # Normalize 0-1 (higher = worse)
    cols <- c("IL1b","SAA","CRP","AttRate")
    summ_n <- summ %>% mutate(across(all_of(cols), ~ (.x - min(.x)) / (max(.x) - min(.x) + 1e-6)))
    summ_long <- pivot_longer(summ_n, cols, names_to="metric", values_to="score")
    ggplot(summ_long, aes(metric, score, fill=drug, group=drug)) +
      geom_bar(stat="identity", position="dodge", alpha=0.75) +
      labs(title="Normalized Burden Score (lower = better)",
           x=NULL, y="Normalized Score (0-1)", fill="Treatment") +
      theme_fmf()
  })

  ## ── Tab 7: Amyloidosis ────────────────────────────────────────────────────
  amy_data <- eventReactive(input$run_amy, {
    end  <- input$amy_years * 8760
    drug <- input$amy_drug
    params <- get_params(drug, input$amy_mefv)
    evs  <- build_events(drug, end)
    run_sim(MOD, params, evs, end, 24) %>%
      mutate(time_yr = time / 8760)
  }, ignoreNULL = FALSE)

  output$amy_SAA_plot <- renderPlot({
    df <- amy_data()
    ggplot(df, aes(time_yr, SAA)) +
      geom_line(color="#e67e22",linewidth=0.9) +
      geom_hline(yintercept=input$amy_SAA_target,linetype="dashed",color="red") +
      labs(x="Year",y="SAA (mg/L)",title="SAA Over Time") + theme_fmf()
  })

  output$amy_dep_plot <- renderPlot({
    df <- amy_data()
    ggplot(df, aes(time_yr, AA_dep)) +
      geom_line(color="#8e44ad",linewidth=0.9) +
      labs(x="Year",y="AA Deposits (AU)",title="AA Amyloid Accumulation") + theme_fmf()
  })

  output$amy_eGFR_plot <- renderPlot({
    df <- amy_data()
    ggplot(df, aes(time_yr, eGFR)) +
      geom_line(color="#27ae60",linewidth=0.9) +
      geom_hline(yintercept=c(60,30),linetype="dashed",color="gray50") +
      annotate("text", x=0, y=62, label="CKD G3", hjust=0, size=3) +
      annotate("text", x=0, y=32, label="CKD G4", hjust=0, size=3) +
      labs(x="Year",y="eGFR (mL/min/1.73m²)",title="eGFR — Long-term Trajectory") +
      theme_fmf()
  })

  output$amy_table <- renderTable({
    df <- amy_data()
    data.frame(
      Metric = c("Final AA deposits (AU)", "Final eGFR",
                 "Time below eGFR 60 (%)", "Mean SAA (mg/L)",
                 "% time SAA > target"),
      Value  = c(
        round(last(df$AA_dep), 3),
        round(last(df$eGFR), 1),
        round(mean(df$eGFR < 60) * 100, 1),
        round(mean(df$SAA), 1),
        round(mean(df$SAA > input$amy_SAA_target) * 100, 1)
      )
    )
  })

  ## ── Tab 8: Sensitivity ─────────────────────────────────────────────────────
  sens_data <- eventReactive(input$run_sens, {
    param_name <- input$sens_param
    vals <- seq(input$sens_min, input$sens_max, length.out = input$sens_n)
    end  <- 8760

    lapply(vals, function(v) {
      p <- list(MEFV_severity = input$sens_mefv,
                USE_COL = 1, USE_ANA = 0, USE_CANA = 0)
      p[[param_name]] <- v
      ev_col <- make_ev("GUT_COL", 1000000, 24, end)
      run_sim(MOD, p, ev_col, end, 12) %>%
        mutate(param_val = v)
    }) %>% bind_rows()
  }, ignoreNULL = FALSE)

  output$sens_tornado <- renderPlot({
    df  <- sens_data()
    var <- input$sens_outcome
    summ <- df %>% group_by(param_val) %>%
      summarise(mean_val = mean(.data[[var]]), .groups="drop")
    ggplot(summ, aes(param_val, mean_val)) +
      geom_line(color="#2c3e50",linewidth=0.9) +
      geom_point(size=3, color="#2c3e50") +
      labs(x = paste("Parameter:", input$sens_param),
           y = paste("Mean", var),
           title = "Sensitivity: Parameter vs Outcome") +
      theme_fmf()
  })

  output$sens_time <- renderPlot({
    df  <- sens_data()
    var <- input$sens_outcome
    ggplot(df, aes_string("time_wk", var, color = "factor(round(param_val, 2))")) +
      geom_line(alpha=0.8) +
      labs(x="Week", y=var,
           title="Time Course by Parameter Value",
           color=input$sens_param) +
      theme_fmf()
  })

  output$sens_table <- renderTable({
    df  <- sens_data()
    var <- input$sens_outcome
    df %>% group_by(`Parameter Value` = round(param_val, 3)) %>%
      summarise(`Mean Outcome` = round(mean(.data[[var]]), 3),
                `Peak Outcome` = round(max(.data[[var]]), 3),
                .groups = "drop")
  })
}

## ─── Run ─────────────────────────────────────────────────────────────────────

shinyApp(ui = ui, server = server)
