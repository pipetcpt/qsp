##############################################################################
# Ankylosing Spondylitis QSP — Interactive Shiny Dashboard
# 7 tabs: Patient Profile · PK · Disease Activity · Structural Progression ·
#          Clinical Endpoints · Scenario Comparison · Biomarkers
##############################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)
library(DT)

# ─────────────────────────────────────────────────────────────────────────────
# Inline mrgsolve model (abbreviated from as_mrgsolve_model.R)
# ─────────────────────────────────────────────────────────────────────────────
as_model_code <- '
$PARAM
USE_ADA=0 USE_ETA=0 USE_SEC=0 USE_TOF=0 USE_UPA=0 USE_NSAID=0
ADA_dose=40  ADA_ka=0.019 ADA_F=0.64  ADA_CL=9.6   ADA_V=7600
ETA_dose=50  ETA_ka=0.048 ETA_F=0.76  ETA_CL=15.4  ETA_V=7600
SEC_dose=150 SEC_ka=0.015 SEC_F=0.73  SEC_CL=4.8   SEC_Vc=3600 SEC_Vp=2800 SEC_Q=1.2
TOF_dose=5   TOF_ka=1.5   TOF_F=0.74  TOF_CL=3500  TOF_V=87000
UPA_dose=15  UPA_ka=1.2   UPA_F=0.79  UPA_CL=8500  UPA_V=96000
NSAID_dose=500 NSAID_ka=0.5 NSAID_F=0.99 NSAID_CL=300 NSAID_V=9600
TNF_base=18  IL17_base=35  IL23_base=25  IL6_base=12  CRP_base=15
kprod_TNF=0.25 kdeg_TNF=0.018 kprod_IL17=0.30 kdeg_IL17=0.012
kprod_IL23=0.20 kdeg_IL23=0.015 kprod_IL6=0.35 kdeg_IL6=0.025
kprod_CRP=0.10 kprod_RANKL=1.2 kdeg_RANKL=0.012 kprod_OPG=0.8 kdeg_OPG=0.010
kform_OC=0.005 kdeg_OC=0.004 kform_BF=0.002 kdeg_BF=0.001
IC50_ADA_TNF=0.3 IC50_ETA_TNF=0.5 IC50_SEC_IL17=0.08
IC50_TOF_JAK=35 IC50_UPA_JAK=10 IC50_NSAID_COX=2500
Emax_TNFi=0.85 Emax_IL17i=0.90 Emax_JAKi=0.70 Emax_NSAID=0.50
BASDAI_ss=6.8 ASDAS_ss=3.9
WT=75 SEX=1 HLA_B27=1 HLAB27_eff=1.15

$CMT
ADA_SC ADA_C ETA_SC ETA_C SEC_SC SEC_C SEC_P TOF_C UPA_C NSAID_C
TNF IL17A IL23 IL6 CRP RANKL OPG OC BF Erosion mSASSS DiseaseAct

$INIT
ADA_SC=0 ADA_C=0 ETA_SC=0 ETA_C=0 SEC_SC=0 SEC_C=0 SEC_P=0
TOF_C=0 UPA_C=0 NSAID_C=0
TNF=18 IL17A=35 IL23=25 IL6=12 CRP=15
RANKL=100 OPG=80 OC=1.0 BF=1.0 Erosion=0 mSASSS=0 DiseaseAct=0.68

$MAIN
double ADA_ke = ADA_CL/ADA_V; double ETA_ke = ETA_CL/ETA_V;
double SEC_ke = SEC_CL/SEC_Vc; double TOF_ke = TOF_CL/TOF_V;
double UPA_ke = UPA_CL/UPA_V; double NSAID_ke = NSAID_CL/NSAID_V;
double cADA = USE_ADA  ? ADA_C /ADA_V *1e6 : 0.0;
double cETA = USE_ETA  ? ETA_C /ETA_V *1e6 : 0.0;
double cSEC = USE_SEC  ? SEC_C /SEC_Vc*1e6 : 0.0;
double cTOF = USE_TOF  ? TOF_C /TOF_V *1e9 : 0.0;
double cUPA = USE_UPA  ? UPA_C /UPA_V *1e9 : 0.0;
double cNSAID= USE_NSAID? NSAID_C/NSAID_V*1e6: 0.0;
double Inh_TNFi  = Emax_TNFi *(cADA+cETA)/(IC50_ADA_TNF+cADA+cETA+1e-10);
double Inh_IL17i = Emax_IL17i*cSEC/(IC50_SEC_IL17+cSEC+1e-10);
double Inh_JAKi  = Emax_JAKi*(cTOF/IC50_TOF_JAK+cUPA/IC50_UPA_JAK)/
                   (1+cTOF/IC50_TOF_JAK+cUPA/IC50_UPA_JAK+1e-10);
double Inh_NSAID = Emax_NSAID*cNSAID/(IC50_NSAID_COX+cNSAID+1e-10);
double HLA_factor = HLA_B27 ? HLAB27_eff : 1.0;
double IL23_drv  = IL23/(IL23_base+IL23);
double IL17_amp  = 1+0.3*IL17A/(IL17_base+IL17A);
double TNF_amp   = 1+0.5*TNF/(TNF_base+TNF);
double IL6_drv   = IL6/IL6_base;
double RANKL_OPG_ratio = (RANKL+1e-10)/(OPG+1e-10);
double TNF_RANKL_eff = 1+0.4*TNF/(TNF_base+TNF);
double TNF_OPG_eff   = 1/(1+0.3*TNF/(TNF_base+TNF));
double DA_raw = 0.4*TNF/TNF_base+0.4*IL17A/IL17_base+0.2*IL6/IL6_base;
double DA_norm = DA_raw/3.0;

$ODE
dxdt_ADA_SC = -ADA_ka*ADA_SC;
dxdt_ADA_C  =  ADA_ka*ADA_SC*ADA_F - ADA_ke*ADA_C;
dxdt_ETA_SC = -ETA_ka*ETA_SC;
dxdt_ETA_C  =  ETA_ka*ETA_SC*ETA_F - ETA_ke*ETA_C;
dxdt_SEC_SC = -SEC_ka*SEC_SC;
dxdt_SEC_C  =  SEC_ka*SEC_SC*SEC_F - SEC_ke*SEC_C-(SEC_Q/SEC_Vc)*SEC_C+(SEC_Q/SEC_Vp)*SEC_P;
dxdt_SEC_P  =  (SEC_Q/SEC_Vc)*SEC_C-(SEC_Q/SEC_Vp)*SEC_P;
dxdt_TOF_C  =  TOF_ka*TOF_dose*TOF_F - TOF_ke*TOF_C;
dxdt_UPA_C  =  UPA_ka*UPA_dose*UPA_F - UPA_ke*UPA_C;
dxdt_NSAID_C=  NSAID_ka*NSAID_dose*NSAID_F - NSAID_ke*NSAID_C;
dxdt_TNF    = kprod_TNF*TNF_base*IL17_amp*HLA_factor - kdeg_TNF*TNF*(1+Inh_TNFi+0.3*Inh_JAKi);
dxdt_IL17A  = kprod_IL17*IL17_base*(1+0.6*IL23_drv)*HLA_factor - kdeg_IL17*IL17A*(1+Inh_IL17i+0.2*Inh_JAKi);
dxdt_IL23   = kprod_IL23*IL23_base - kdeg_IL23*IL23*(1+0.4*Inh_JAKi);
dxdt_IL6    = kprod_IL6*IL6_base*TNF_amp - kdeg_IL6*IL6*(1+0.5*Inh_JAKi);
dxdt_CRP    = kprod_CRP*CRP_base*IL6_drv - 0.036*CRP;
dxdt_RANKL  = kprod_RANKL/24*RANKL_base*TNF_RANKL_eff - kdeg_RANKL*RANKL*(1+0.3*Inh_TNFi);
dxdt_OPG    = kprod_OPG/24*OPG_base*TNF_OPG_eff - kdeg_OPG*OPG;
dxdt_OC     = kform_OC*RANKL_OPG_ratio*(1+0.3*TNF/TNF_base) - kdeg_OC*OC*(1+0.4*Inh_TNFi);
dxdt_Erosion= 0.0003*OC*(1-0.7*Inh_TNFi);
double IL17_bone_drv = IL17A/IL17_base;
double TNF_DKK1_inh  = 1-0.4*Inh_TNFi;
dxdt_BF     = kform_BF*(1+0.5*IL17_bone_drv)*TNF_DKK1_inh - kdeg_BF*BF;
dxdt_mSASSS = 0.00015*BF*(1-0.5*Inh_IL17i-0.3*Inh_TNFi);
double DA_new = DA_norm*(1-0.8*Inh_TNFi)*(1-0.7*Inh_IL17i)*(1-0.6*Inh_JAKi)*(1-0.4*Inh_NSAID);
dxdt_DiseaseAct = 0.01*(DA_new-DiseaseAct);

$TABLE
double ADA_ugmL = ADA_C/ADA_V*1e6;
double ETA_ugmL = ETA_C/ETA_V*1e6;
double SEC_ugmL = SEC_C/SEC_Vc*1e6;
double TOF_ngmL = TOF_C/TOF_V*1e9;
double UPA_ngmL = UPA_C/UPA_V*1e9;
double BASDAI_sim = BASDAI_ss*DiseaseAct/0.68;
if(BASDAI_sim<0) BASDAI_sim=0; if(BASDAI_sim>10) BASDAI_sim=10;
double ASDAS_sim  = 0.75*DiseaseAct/0.68*ASDAS_ss;
if(ASDAS_sim<0) ASDAS_sim=0;
double logit20   = -2.2+3.5*(BASDAI_ss-BASDAI_sim)/BASDAI_ss;
double ASAS20_prob = 1/(1+exp(-logit20));
double logit40   = -3.5+3.5*(BASDAI_ss-BASDAI_sim)/BASDAI_ss;
double ASAS40_prob = 1/(1+exp(-logit40));
double RANKL_OPG_out = RANKL/(OPG+1e-10);
double Eff_TNFi = Inh_TNFi; double Eff_IL17i = Inh_IL17i; double Eff_JAKi = Inh_JAKi;

$CAPTURE
ADA_ugmL ETA_ugmL SEC_ugmL TOF_ngmL UPA_ngmL
TNF IL17A IL23 IL6 CRP RANKL OPG RANKL_OPG_out
OC Erosion BF mSASSS DiseaseAct BASDAI_sim ASDAS_sim
ASAS20_prob ASAS40_prob Eff_TNFi Eff_IL17i Eff_JAKi
'

as_mod <- mcode("AS_QSP_shiny", as_model_code, quiet = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# Helper: build dosing events
# ─────────────────────────────────────────────────────────────────────────────
make_ev <- function(drug, dur_wks) {
  dur_h <- dur_wks * 168
  switch(drug,
    ADA   = ev(cmt="ADA_SC", amt=40,  time=seq(0,dur_h-1,by=336)),
    ETA   = ev(cmt="ETA_SC", amt=50,  time=seq(0,dur_h-1,by=168)),
    SEC   = {
      load  <- c(0,1,2,3,4)*168
      maint <- seq(4*168, dur_h-1, by=4*168)
      ev(cmt="SEC_SC", amt=150, time=unique(c(load,maint)))
    },
    TOF   = ev(cmt="TOF_C", amt=5*0.74, time=seq(0,dur_h-1,by=12)),
    UPA   = ev(cmt="UPA_C", amt=15*0.79,time=seq(0,dur_h-1,by=24)),
    NSAID = ev(cmt="NSAID_C",amt=500*0.99,time=seq(0,dur_h-1,by=12)),
    PBO   = ev(amt=0, time=0, cmt="ADA_SC")
  )
}

run_sim <- function(drug, dur_wks=52, hla_pos=TRUE, wt=75) {
  flag_map <- c(ADA="USE_ADA", ETA="USE_ETA", SEC="USE_SEC",
                TOF="USE_TOF", UPA="USE_UPA", NSAID="USE_NSAID")
  par_set  <- list(WT=wt, HLA_B27=as.integer(hla_pos))
  if (drug != "PBO") par_set[[flag_map[drug]]] <- 1
  as_mod %>%
    param(par_set) %>%
    ev(make_ev(drug, dur_wks)) %>%
    mrgsim(end=dur_wks*168, delta=24) %>%
    as_tibble() %>%
    mutate(time_wk = time/168)
}

# ─────────────────────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────────────────────
drug_choices <- c(
  "Placebo"               = "PBO",
  "Adalimumab 40mg SC q2w"= "ADA",
  "Etanercept 50mg SC qw" = "ETA",
  "Secukinumab 150mg SC"  = "SEC",
  "Tofacitinib 5mg BID"   = "TOF",
  "Upadacitinib 15mg QD"  = "UPA",
  "NSAID (naproxen)"      = "NSAID"
)

scenario_colors <- c(
  PBO   = "#95A5A6", ADA = "#2E86C1", ETA = "#1A5276",
  SEC   = "#E74C3C", TOF = "#F39C12", UPA = "#8E44AD", NSAID = "#27AE60"
)

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(
    title = "AS QSP Model",
    titleWidth = 250
  ),

  dashboardSidebar(
    width = 250,
    sidebarMenu(
      menuItem("Patient Profile",      tabName = "tab_patient",   icon = icon("user")),
      menuItem("Drug PK",              tabName = "tab_pk",        icon = icon("pills")),
      menuItem("Disease Activity",     tabName = "tab_da",        icon = icon("chart-line")),
      menuItem("Structural Progress.", tabName = "tab_struct",    icon = icon("bone")),
      menuItem("Clinical Endpoints",   tabName = "tab_endpoints", icon = icon("stethoscope")),
      menuItem("Scenario Comparison",  tabName = "tab_compare",   icon = icon("balance-scale")),
      menuItem("Biomarker Dashboard",  tabName = "tab_biomarker", icon = icon("vial"))
    ),
    hr(),
    h5("  Global Settings", style="color:#ccc; margin-left:10px"),
    sliderInput("dur_wks", "Duration (weeks)", 12, 104, 52, step=4),
    checkboxInput("hla_pos", "HLA-B27 Positive", value=TRUE),
    sliderInput("wt", "Body Weight (kg)", 40, 130, 75),
    actionButton("run_btn", "Run Simulation", class="btn-primary",
                 style="margin:10px; width:220px")
  ),

  dashboardBody(
    tags$style(HTML(".box-header { background:#2E86C1 !important; color:white !important; }
                     .box { border-top:3px solid #2E86C1; }")),
    tabItems(

      # ── Tab 1: Patient Profile ──────────────────────────────────────────
      tabItem("tab_patient",
        fluidRow(
          box(title="AS Patient Profile", width=6, status="primary", solidHeader=TRUE,
            selectInput("drug1", "Treatment", choices=drug_choices, selected="ADA"),
            br(),
            h4("Disease Characteristics"),
            sliderInput("basdai_init", "Baseline BASDAI", 0, 10, 6.8, 0.1),
            sliderInput("crp_init",    "Baseline CRP (mg/L)", 1, 60, 15, 1),
            sliderInput("tnf_init",    "Baseline TNF-α (pg/mL)", 5, 80, 18, 1),
            sliderInput("il17_init",   "Baseline IL-17A (pg/mL)", 5, 100, 35, 1),
            checkboxInput("prior_tnfi", "Prior TNFi failure", value=FALSE),
            numericInput("msasss_init", "Baseline mSASSS", value=0, min=0, max=72)
          ),
          box(title="Patient Summary", width=6, status="info", solidHeader=TRUE,
            h4("Disease Classification"),
            verbatimTextOutput("pt_summary"),
            hr(),
            h4("Predicted Drug Class Suitability"),
            plotOutput("pt_radar", height="260px")
          )
        ),
        fluidRow(
          box(title="AS Classification Criteria (ASAS 2009)", width=12, status="warning",
            collapsible=TRUE, collapsed=TRUE,
            DTOutput("asas_criteria_tbl")
          )
        )
      ),

      # ── Tab 2: Drug PK ──────────────────────────────────────────────────
      tabItem("tab_pk",
        fluidRow(
          box(title="PK Parameters", width=4, status="primary", solidHeader=TRUE,
            selectInput("pk_drug", "Select Drug",
                        choices=c("Adalimumab"="ADA","Etanercept"="ETA",
                                  "Secukinumab"="SEC","Tofacitinib"="TOF",
                                  "Upadacitinib"="UPA"),
                        selected="ADA"),
            DTOutput("pk_params_tbl")
          ),
          box(title="Concentration–Time Profile", width=8, status="primary",
              solidHeader=TRUE,
            plotlyOutput("pk_plot", height="350px")
          )
        ),
        fluidRow(
          box(title="PK Summary (Cmax, Tmax, AUC)", width=12,
            DTOutput("pk_summary_tbl")
          )
        )
      ),

      # ── Tab 3: Disease Activity ─────────────────────────────────────────
      tabItem("tab_da",
        fluidRow(
          box(title="BASDAI & ASDAS-CRP Over Time", width=8, status="primary",
              solidHeader=TRUE,
            plotlyOutput("da_plot", height="380px")
          ),
          box(title="Response Thresholds", width=4, status="info",
            h5("BASDAI Response Levels"),
            tags$ul(
              tags$li("BASDAI < 4: Low / Inactive"),
              tags$li("BASDAI < 2: Partial Remission"),
              tags$li("BASDAI50: ≥50% improvement")
            ),
            h5("ASDAS-CRP Response"),
            tags$ul(
              tags$li("< 1.3: Inactive disease"),
              tags$li("1.3–2.1: Low activity"),
              tags$li("2.1–3.5: High activity"),
              tags$li("> 3.5: Very high activity")
            ),
            hr(),
            h5("Simulated Wk-24 Outcomes"),
            tableOutput("da_wk24_tbl")
          )
        ),
        fluidRow(
          box(title="CRP Dynamics", width=6,
            plotlyOutput("crp_plot", height="280px")
          ),
          box(title="IL-17A and TNF-α Suppression", width=6,
            plotlyOutput("cytokine_da_plot", height="280px")
          )
        )
      ),

      # ── Tab 4: Structural Progression ──────────────────────────────────
      tabItem("tab_struct",
        fluidRow(
          box(title="mSASSS Over Time (Structural Progression)", width=8,
              status="primary", solidHeader=TRUE,
            plotlyOutput("msasss_plot", height="360px")
          ),
          box(title="Bone Biology Explanation", width=4, status="warning",
            h5("Dual Bone Pathology in AS"),
            p("AS involves both bone erosion AND abnormal new bone formation (syndesmophytes)."),
            hr(),
            h5("Key Pathways"),
            tags$ul(
              tags$li("TNFi: reduces inflammation but may paradoxically allow more bone formation via DKK1 removal"),
              tags$li("IL-17Ai (secukinumab): reduces both inflammation AND syndesmophyte formation"),
              tags$li("JAKi: anti-inflammatory but bone formation data limited")
            ),
            hr(),
            verbatimTextOutput("bone_summary")
          )
        ),
        fluidRow(
          box(title="Bone Formation Index & Erosion", width=6,
            plotlyOutput("bone_plot", height="280px")
          ),
          box(title="RANKL:OPG Ratio", width=6,
            plotlyOutput("rankl_plot", height="280px")
          )
        )
      ),

      # ── Tab 5: Clinical Endpoints ───────────────────────────────────────
      tabItem("tab_endpoints",
        fluidRow(
          box(title="ASAS Response Probability Over Time", width=8,
              status="primary", solidHeader=TRUE,
            plotlyOutput("asas_plot", height="360px")
          ),
          box(title="Trial Benchmarks (Wk 16–24)", width=4, status="success",
            h5("Phase 3 RCT ASAS20 Rates"),
            DTOutput("trial_bench_tbl")
          )
        ),
        fluidRow(
          box(title="Response at Key Timepoints", width=12,
            DTOutput("response_tbl")
          )
        )
      ),

      # ── Tab 6: Scenario Comparison ──────────────────────────────────────
      tabItem("tab_compare",
        fluidRow(
          box(title="Select Scenarios to Compare", width=4, status="primary",
              solidHeader=TRUE,
            checkboxGroupInput("compare_drugs", "Treatments:",
              choices = drug_choices,
              selected = c("PBO","ADA","SEC","UPA")
            )
          ),
          box(title="BASDAI Comparison", width=8,
            plotlyOutput("compare_basdai", height="320px")
          )
        ),
        fluidRow(
          box(title="Cytokine Suppression Comparison", width=6,
            plotlyOutput("compare_cyt", height="280px")
          ),
          box(title="mSASSS Comparison", width=6,
            plotlyOutput("compare_msasss", height="280px")
          )
        ),
        fluidRow(
          box(title="Summary Table at Week 24", width=12,
            DTOutput("compare_tbl")
          )
        )
      ),

      # ── Tab 7: Biomarker Dashboard ──────────────────────────────────────
      tabItem("tab_biomarker",
        fluidRow(
          box(title="Drug & Biomarker Settings", width=3, status="primary",
              solidHeader=TRUE,
            selectInput("bm_drug", "Treatment", choices=drug_choices, selected="SEC"),
            hr(),
            h5("Key AS Biomarkers"),
            checkboxGroupInput("bm_vars", "Show biomarkers:",
              choices=c("TNF-α"="TNF","IL-17A"="IL17A","IL-23"="IL23",
                        "IL-6"="IL6","CRP"="CRP","RANKL"="RANKL","OPG"="OPG"),
              selected=c("TNF","IL17A","CRP")
            )
          ),
          box(title="Biomarker Trajectories", width=9,
            plotlyOutput("bm_plot", height="400px")
          )
        ),
        fluidRow(
          box(title="Drug Efficacy Metrics Over Time", width=6,
            plotlyOutput("eff_plot", height="280px")
          ),
          box(title="Osteoclast Activity & Bone Markers", width=6,
            plotlyOutput("oc_plot", height="280px")
          )
        )
      )
    )
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# SERVER
# ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Reactive: run primary simulation
  sim_data <- eventReactive(input$run_btn, {
    withProgress(message="Running simulation...", {
      run_sim(input$drug1, input$dur_wks, input$hla_pos, input$wt)
    })
  }, ignoreNULL=FALSE)

  # Reactive: comparison simulations
  compare_data <- eventReactive(input$run_btn, {
    withProgress(message="Running comparison...", {
      drugs <- input$compare_drugs
      bind_rows(lapply(drugs, function(d) {
        run_sim(d, input$dur_wks, input$hla_pos, input$wt) %>%
          mutate(drug = d,
                 scenario = names(drug_choices)[drug_choices == d])
      }))
    })
  }, ignoreNULL=FALSE)

  # ── Tab 1: Patient Summary ────────────────────────────────────────────
  output$pt_summary <- renderText({
    hla <- if(input$hla_pos) "HLA-B27 Positive" else "HLA-B27 Negative"
    prior <- if(input$prior_tnfi) "Prior TNFi failure: YES" else "Prior TNFi failure: No"
    class_as <- if(input$basdai_init >= 4) "Active AS (BASDAI ≥ 4)" else "Low-activity AS (BASDAI < 4)"
    asdas_est <- 0.75 * input$basdai_init/6.8 * 3.9
    drug_name <- names(drug_choices)[drug_choices == input$drug1]
    paste0(
      "Classification  : ", class_as, "\n",
      "HLA-B27 Status  : ", hla, "\n",
      "Body Weight     : ", input$wt, " kg\n",
      "Baseline BASDAI : ", input$basdai_init, "\n",
      "Estimated ASDAS : ", round(asdas_est, 1), "\n",
      "Baseline CRP    : ", input$crp_init, " mg/L\n",
      "Selected Drug   : ", drug_name, "\n",
      prior
    )
  })

  output$pt_radar <- renderPlot({
    # Simple bar showing predicted efficacy by drug class
    eff_data <- data.frame(
      drug  = c("TNFi","IL-17Ai","JAKi","NSAID"),
      ASAS20= c(58, 61, 52, 35),
      note  = c("Best for uveitis","Best structural","Versatile","Symptom relief")
    )
    eff_data$highlight <- eff_data$drug == c(
      ADA="TNFi",ETA="TNFi",SEC="IL-17Ai",TOF="JAKi",
      UPA="JAKi",NSAID="NSAID",PBO="none")[input$drug1]
    ggplot(eff_data, aes(x=reorder(drug,-ASAS20), y=ASAS20, fill=highlight)) +
      geom_bar(stat="identity") +
      geom_text(aes(label=paste0(ASAS20,"%")), vjust=-0.3, size=4) +
      scale_fill_manual(values=c(`TRUE`="#2E86C1",`FALSE`="#BDC3C7")) +
      labs(x=NULL, y="ASAS20 Rate (%)",
           title="Expected ASAS20 by Drug Class (Wk 24)") +
      theme_bw() + theme(legend.position="none")
  })

  output$asas_criteria_tbl <- renderDT({
    tibble(
      Criterion = c("Sacroiliitis on imaging + ≥1 SpA feature",
                    "HLA-B27 + ≥2 SpA features"),
      `SpA Features` = c(
        "Inflammatory back pain, Arthritis, Enthesitis (heel), Uveitis, Dactylitis, Psoriasis, Crohn's/colitis, Good response to NSAIDs, Family history SpA, HLA-B27, Elevated CRP",
        "(same list)")
    ) %>% datatable(options=list(dom="t"), rownames=FALSE)
  })

  # ── Tab 2: PK ─────────────────────────────────────────────────────────
  pk_params_df <- reactive({
    d <- input$pk_drug
    params <- list(
      ADA  = tibble(Parameter=c("ka (/h)","F","CL (mL/h)","V (mL)","t½ (d)"),
                    Value=c("0.019","0.64","9.6","7,600","14")),
      ETA  = tibble(Parameter=c("ka (/h)","F","CL (mL/h)","V (mL)","t½ (d)"),
                    Value=c("0.048","0.76","15.4","7,600","4")),
      SEC  = tibble(Parameter=c("ka (/h)","F","CL (mL/h)","Vc (mL)","Vp (mL)","Q (mL/h)","t½ (d)"),
                    Value=c("0.015","0.73","4.8","3,600","2,800","1.2","27")),
      TOF  = tibble(Parameter=c("ka (/h)","F","CL (mL/h)","V (mL)","t½ (h)"),
                    Value=c("1.5","0.74","3,500","87,000","3")),
      UPA  = tibble(Parameter=c("ka (/h)","F","CL (mL/h)","V (mL)","t½ (h)"),
                    Value=c("1.2","0.79","8,500","96,000","15"))
    )
    params[[d]]
  })

  output$pk_params_tbl <- renderDT({
    pk_params_df() %>% datatable(options=list(dom="t"), rownames=FALSE)
  })

  output$pk_plot <- renderPlotly({
    d <- sim_data()
    pk_col <- switch(input$pk_drug,
      ADA="ADA_ugmL", ETA="ETA_ugmL", SEC="SEC_ugmL",
      TOF="TOF_ngmL", UPA="UPA_ngmL")
    y_label <- if(input$pk_drug %in% c("TOF","UPA")) "Conc (ng/mL)" else "Conc (μg/mL)"
    p <- ggplot(d, aes(x=time_wk, y=.data[[pk_col]])) +
      geom_line(color="#2E86C1", linewidth=1.2) +
      labs(title=paste(names(drug_choices)[drug_choices==input$pk_drug], "PK"),
           x="Time (weeks)", y=y_label) +
      theme_bw()
    ggplotly(p)
  })

  output$pk_summary_tbl <- renderDT({
    d <- sim_data()
    tibble(
      Drug     = c("Adalimumab","Secukinumab","Tofacitinib","Upadacitinib"),
      Cmax     = c(
        round(max(d$ADA_ugmL),2),
        round(max(d$SEC_ugmL),2),
        round(max(d$TOF_ngmL),1),
        round(max(d$UPA_ngmL),1)),
      Unit     = c("μg/mL","μg/mL","ng/mL","ng/mL"),
      Tmax_wk  = c(
        round(d$time_wk[which.max(d$ADA_ugmL)],1),
        round(d$time_wk[which.max(d$SEC_ugmL)],1),
        round(d$time_wk[which.max(d$TOF_ngmL)],1),
        round(d$time_wk[which.max(d$UPA_ngmL)],1))
    ) %>% datatable(options=list(dom="t"), rownames=FALSE)
  })

  # ── Tab 3: Disease Activity ───────────────────────────────────────────
  output$da_plot <- renderPlotly({
    d <- sim_data() %>%
      select(time_wk, BASDAI_sim, ASDAS_sim) %>%
      pivot_longer(-time_wk, names_to="index", values_to="value")
    p <- ggplot(d, aes(x=time_wk, y=value, color=index)) +
      geom_line(linewidth=1.2) +
      geom_hline(data=data.frame(index=c("BASDAI_sim","ASDAS_sim"),
                                  thresh=c(4.0, 2.1)),
                 aes(yintercept=thresh, color=index), linetype="dashed") +
      scale_color_manual(values=c(BASDAI_sim="#E74C3C", ASDAS_sim="#3498DB"),
                         labels=c("BASDAI (0–10)","ASDAS-CRP")) +
      labs(title="BASDAI & ASDAS-CRP Over Time",
           x="Time (weeks)", y="Score", color="Index") +
      theme_bw()
    ggplotly(p)
  })

  output$da_wk24_tbl <- renderTable({
    d <- sim_data() %>% filter(abs(time_wk - 24) < 0.1) %>% slice(1)
    tibble(
      Outcome = c("BASDAI","ASDAS-CRP","CRP (mg/L)","ASAS20 (%)","ASAS40 (%)"),
      Value   = c(round(d$BASDAI_sim,1), round(d$ASDAS_sim,2),
                  round(d$CRP,1),
                  round(d$ASAS20_prob*100,0), round(d$ASAS40_prob*100,0))
    )
  })

  output$crp_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=time_wk, y=CRP)) +
      geom_line(color="#E74C3C", linewidth=1.2) +
      geom_hline(yintercept=5, linetype="dashed", color="gray50") +
      labs(title="CRP Over Time", x="Time (weeks)", y="CRP (mg/L)") +
      theme_bw()
    ggplotly(p)
  })

  output$cytokine_da_plot <- renderPlotly({
    d <- sim_data() %>%
      select(time_wk, TNF, IL17A) %>%
      pivot_longer(-time_wk, names_to="cytokine", values_to="conc")
    p <- ggplot(d, aes(x=time_wk, y=conc, color=cytokine)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=c(TNF="#E74C3C", IL17A="#3498DB"),
                         labels=c("IL-17A (pg/mL)","TNF-α (pg/mL)")) +
      labs(title="Key Cytokine Suppression",
           x="Time (weeks)", y="Concentration (pg/mL)", color="Cytokine") +
      theme_bw()
    ggplotly(p)
  })

  # ── Tab 4: Structural Progression ────────────────────────────────────
  output$msasss_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=time_wk, y=mSASSS)) +
      geom_line(color="#8E44AD", linewidth=1.5) +
      labs(title="mSASSS (Structural Progression)",
           subtitle="Lower = better; IL-17A inhibition > TNF inhibition for bone",
           x="Time (weeks)", y="mSASSS (AU)") +
      theme_bw()
    ggplotly(p)
  })

  output$bone_summary <- renderText({
    d <- sim_data() %>% filter(abs(time_wk - 52) < 0.1) %>% slice(1)
    paste0(
      "mSASSS at wk 52 : ", round(d$mSASSS, 3), "\n",
      "Erosion score   : ", round(d$Erosion, 4), "\n",
      "Bone Form. Index: ", round(d$BF, 3), "\n",
      "OC activity     : ", round(d$OC, 3)
    )
  })

  output$bone_plot <- renderPlotly({
    d <- sim_data() %>%
      select(time_wk, BF, Erosion) %>%
      pivot_longer(-time_wk, names_to="marker", values_to="value")
    p <- ggplot(d, aes(x=time_wk, y=value, color=marker)) +
      geom_line(linewidth=1.1) +
      scale_color_manual(values=c(BF="#27AE60", Erosion="#E74C3C"),
                         labels=c("Bone Formation Index","Bone Erosion Score")) +
      labs(title="Bone Remodeling Markers",
           x="Time (weeks)", y="Score (AU)", color="Marker") +
      theme_bw()
    ggplotly(p)
  })

  output$rankl_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=time_wk, y=RANKL_OPG_out)) +
      geom_line(color="#C0392B", linewidth=1.2) +
      geom_hline(yintercept=1.25, linetype="dashed", color="gray50") +
      labs(title="RANKL:OPG Ratio",
           subtitle="Ratio > 1 → net osteoclast activation",
           x="Time (weeks)", y="RANKL:OPG Ratio") +
      theme_bw()
    ggplotly(p)
  })

  # ── Tab 5: Clinical Endpoints ─────────────────────────────────────────
  output$asas_plot <- renderPlotly({
    d <- sim_data() %>%
      select(time_wk, ASAS20_prob, ASAS40_prob) %>%
      pivot_longer(-time_wk, names_to="endpoint", values_to="prob")
    p <- ggplot(d, aes(x=time_wk, y=prob*100, color=endpoint)) +
      geom_line(linewidth=1.2) +
      geom_hline(data=data.frame(endpoint=c("ASAS20_prob","ASAS40_prob"),
                                  thresh=c(22,11)),
                 aes(yintercept=thresh, color=endpoint), linetype="dashed") +
      scale_color_manual(values=c(ASAS20_prob="#2E86C1", ASAS40_prob="#E74C3C"),
                         labels=c("ASAS20","ASAS40")) +
      labs(title="ASAS20/40 Response Probability",
           subtitle="Dashed = ATLAS placebo rates",
           x="Time (weeks)", y="Response (%)", color="Endpoint") +
      theme_bw()
    ggplotly(p)
  })

  output$trial_bench_tbl <- renderDT({
    tibble(
      Trial          = c("ATLAS wk24 ADA","ATLAS wk24 PBO",
                         "MEASURE1 wk16 SEC","MEASURE1 wk16 PBO",
                         "COAST-V wk16 IXE","SELECT-AXIS1 wk14 UPA"),
      ASAS20         = c("59%","22%","61%","29%","52%","52%"),
      ASAS40         = c("45%","11%","36%","13%","48%","40%"),
      BASDAI50       = c("50%","13%","40%","14%","NA","NA")
    ) %>% datatable(options=list(dom="t", pageLength=10), rownames=FALSE)
  })

  output$response_tbl <- renderDT({
    d <- sim_data()
    wks <- c(4, 12, 16, 24, 52)
    map_dfr(wks, function(w) {
      r <- d %>% filter(abs(time_wk - w) < 0.2) %>% slice(1)
      tibble(
        `Week`       = w,
        `BASDAI`     = round(r$BASDAI_sim, 1),
        `ASDAS-CRP`  = round(r$ASDAS_sim, 2),
        `CRP (mg/L)` = round(r$CRP, 1),
        `ASAS20 (%)`  = round(r$ASAS20_prob * 100, 0),
        `ASAS40 (%)`  = round(r$ASAS40_prob * 100, 0),
        `mSASSS`      = round(r$mSASSS, 3)
      )
    }) %>% datatable(options=list(dom="t"), rownames=FALSE)
  })

  # ── Tab 6: Scenario Comparison ────────────────────────────────────────
  output$compare_basdai <- renderPlotly({
    d <- compare_data()
    p <- ggplot(d, aes(x=time_wk, y=BASDAI_sim, color=drug, group=scenario)) +
      geom_line(linewidth=1.1) +
      geom_hline(yintercept=4, linetype="dashed", color="gray50") +
      scale_color_manual(values=scenario_colors,
                         labels=names(drug_choices)[match(unique(d$drug), drug_choices)]) +
      labs(title="BASDAI Comparison", x="Time (weeks)", y="BASDAI", color="Drug") +
      theme_bw()
    ggplotly(p)
  })

  output$compare_cyt <- renderPlotly({
    d <- compare_data()
    p <- ggplot(d, aes(x=time_wk, y=TNF, color=drug)) +
      geom_line(linewidth=1) +
      scale_color_manual(values=scenario_colors) +
      labs(title="TNF-α Suppression", x="Time (weeks)", y="TNF-α (pg/mL)",
           color="Drug") +
      theme_bw()
    ggplotly(p)
  })

  output$compare_msasss <- renderPlotly({
    d <- compare_data()
    p <- ggplot(d, aes(x=time_wk, y=mSASSS, color=drug)) +
      geom_line(linewidth=1) +
      scale_color_manual(values=scenario_colors) +
      labs(title="mSASSS Progression", x="Time (weeks)", y="mSASSS (AU)",
           color="Drug") +
      theme_bw()
    ggplotly(p)
  })

  output$compare_tbl <- renderDT({
    d <- compare_data()
    d %>% filter(abs(time_wk - 24) < 0.2) %>%
      mutate(drug_name = names(drug_choices)[match(drug, drug_choices)]) %>%
      select(drug_name, BASDAI_sim, ASDAS_sim, CRP, TNF, IL17A,
             ASAS20_prob, ASAS40_prob, mSASSS) %>%
      mutate(across(c(ASAS20_prob, ASAS40_prob), ~round(.*100, 0)),
             across(c(BASDAI_sim, ASDAS_sim, CRP, TNF, IL17A, mSASSS),
                    ~round(., 2))) %>%
      rename(`Drug`="drug_name", `BASDAI`="BASDAI_sim", `ASDAS-CRP`="ASDAS_sim",
             `ASAS20 (%)`="ASAS20_prob", `ASAS40 (%)`="ASAS40_prob") %>%
      datatable(options=list(dom="t"), rownames=FALSE)
  })

  # ── Tab 7: Biomarker Dashboard ────────────────────────────────────────
  bm_data <- reactive({
    run_sim(input$bm_drug, input$dur_wks, input$hla_pos, input$wt)
  })

  output$bm_plot <- renderPlotly({
    d     <- bm_data()
    vars  <- input$bm_vars
    if (length(vars) == 0) return(NULL)
    dm <- d %>% select(time_wk, all_of(vars)) %>%
      pivot_longer(-time_wk, names_to="biomarker", values_to="value")
    bm_colors <- c(TNF="#E74C3C", IL17A="#3498DB", IL23="#9B59B6",
                   IL6="#F39C12", CRP="#E67E22", RANKL="#C0392B", OPG="#27AE60")
    p <- ggplot(dm, aes(x=time_wk, y=value, color=biomarker)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=bm_colors) +
      facet_wrap(~biomarker, scales="free_y") +
      labs(title="Biomarker Trajectories",
           x="Time (weeks)", y="Concentration", color="Biomarker") +
      theme_bw(base_size=11) +
      theme(legend.position="none")
    ggplotly(p)
  })

  output$eff_plot <- renderPlotly({
    d <- bm_data() %>%
      select(time_wk, Eff_TNFi, Eff_IL17i, Eff_JAKi) %>%
      pivot_longer(-time_wk, names_to="eff_type", values_to="inhibition")
    p <- ggplot(d, aes(x=time_wk, y=inhibition*100, color=eff_type)) +
      geom_line(linewidth=1.1) +
      scale_color_manual(values=c(Eff_TNFi="#E74C3C", Eff_IL17i="#3498DB",
                                   Eff_JAKi="#F39C12"),
                         labels=c("JAK inhibition","IL-17A inh.","TNF inh.")) +
      labs(title="Drug Efficacy (% Inhibition)", x="Time (weeks)",
           y="% Inhibition", color="Pathway") +
      ylim(0, 100) +
      theme_bw()
    ggplotly(p)
  })

  output$oc_plot <- renderPlotly({
    d <- bm_data() %>%
      select(time_wk, OC, BF, Erosion) %>%
      mutate(Erosion_scaled = Erosion * 10) %>%
      pivot_longer(-time_wk, names_to="marker", values_to="value") %>%
      filter(marker %in% c("OC","BF"))
    p <- ggplot(d, aes(x=time_wk, y=value, color=marker)) +
      geom_line(linewidth=1.1) +
      scale_color_manual(values=c(OC="#E74C3C", BF="#27AE60"),
                         labels=c("Bone Formation Index","Osteoclast Activity")) +
      labs(title="Osteoclast & Bone Formation",
           x="Time (weeks)", y="Activity (AU)", color="Marker") +
      theme_bw()
    ggplotly(p)
  })
}

# ─────────────────────────────────────────────────────────────────────────────
# LAUNCH
# ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
