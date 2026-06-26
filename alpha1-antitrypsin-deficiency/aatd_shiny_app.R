## ============================================================
## AATD QSP Shiny Dashboard — Alpha-1 Antitrypsin Deficiency
## ============================================================
## Tabs:
##  1. Patient Profile & Genotype
##  2. Drug PK — Serum AAT Levels
##  3. Pulmonary PD — NE, Elastin, FEV1
##  4. Clinical Endpoints — FEV1, SGRQ, CT, Exacerbations
##  5. Scenario Comparison — All 6 treatment arms
##  6. Biomarkers — NE, Desmosine, Liver (ALT, FIB-4, polymer)
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(shinycssloaders)

## ── Embedded mrgsolve Model ──────────────────────────────────
aatd_code <- '
$PROB AATD QSP Shiny Model

$PARAM
k_ZAAT_synth=0.070 k_poly=0.45 k_ERAD=0.15 k_auto=0.08
k_liver_inj=0.012 k_ERstress=0.18
k10_aug=0.052 k12_aug=0.032 k21_aug=0.022
ZAAT_sec_rate=0.35 k_AAT_elim=0.154 r_ELF_plasma=0.10
k_PMN_base=0.10 k_IL8_recruit=0.08 k_PMN_egress=0.35
k_NE_release=0.55 k_NE_elim=0.40 k_MMP12_base=0.05
k_MMP12_NE=0.03 k_MMP12_elim=0.25
Elastin0=100 k_elastin_deg=0.012 k_elastin_syn=0.002
FEV1_0=95 k_FEV1_Edecay=0.22 annual_decline_base=0.25
k_HSC_activ=0.025 k_HSC_deactiv=0.08
k_coll_synth=0.035 k_coll_deg=0.018
k_fibrosis_p=0.010 k_fibrosis_r=0.005
TGFb_max=2.5 EC50_TGFb=0.5
Emax_siRNA=0.88 EC50_siRNA=0.35 k_siRNA_off=0.007
Emax_NEi=0.75 EC50_NEi=0.18 ka_NEi=1.20 k_NEi_elim=2.40
k_gene_onset=0.005 k_gene_wane=0.00050 Emax_gene=0.90

$INIT
ZAAT_ER=5 ZAAT_Poly=2 HSC_act=0.10 Liver_coll=0.50 Liver_fib=0.10
AAT_C1=5.50 AAT_C2=2.50
PMN_lung=1.0 IL8_lung=0.80 NE_free=0.60 MMP12_lung=0.30
Elastin=92 FEV1_pct=90
AUG_C1=0 AUG_C2=0 NEi_A=0 NEi_C=0 siRNA_Eff=0 Gene_Eff=0

$ODE
dxdt_AUG_C1 = -k10_aug*AUG_C1 - k12_aug*AUG_C1 + k21_aug*AUG_C2;
dxdt_AUG_C2 = k12_aug*AUG_C1 - k21_aug*AUG_C2;
dxdt_NEi_A = -ka_NEi*NEi_A;
dxdt_NEi_C = ka_NEi*NEi_A - k_NEi_elim*NEi_C;
double NEi_eff = Emax_NEi*NEi_C/(EC50_NEi+NEi_C);
dxdt_siRNA_Eff = -k_siRNA_off*siRNA_Eff;
double siRNA_i = Emax_siRNA*siRNA_Eff/(EC50_siRNA+siRNA_Eff);
dxdt_Gene_Eff = -k_gene_wane*Gene_Eff;
double gene_f = Emax_gene*Gene_Eff/(1.0+Gene_Eff);
double ZAAT_end = ZAAT_sec_rate*(1.0-siRNA_i);
double Gene_M = gene_f*50.0;
dxdt_AAT_C1 = ZAAT_end+Gene_M-k_AAT_elim*AAT_C1-k12_aug*AAT_C1+k21_aug*AAT_C2;
dxdt_AAT_C2 = k12_aug*AAT_C1-k21_aug*AAT_C2-k_AAT_elim*AAT_C2*0.5;
double AAT_s = (AAT_C1>0)?AAT_C1:0.01;
dxdt_ZAAT_ER = 8.0*k_ZAAT_synth*(1.0-siRNA_i)-k_poly*ZAAT_ER-k_ERAD*ZAAT_ER;
dxdt_ZAAT_Poly = k_poly*ZAAT_ER-k_auto*ZAAT_Poly;
double li_rate = k_liver_inj*ZAAT_Poly;
double TGFb = TGFb_max*HSC_act/(EC50_TGFb+HSC_act);
dxdt_HSC_act = k_HSC_activ*li_rate+TGFb*0.15*(1.0-HSC_act/3.0)-k_HSC_deactiv*HSC_act;
dxdt_Liver_coll = k_coll_synth*HSC_act-k_coll_deg*Liver_coll;
double fib_r = k_fibrosis_p*Liver_coll-k_fibrosis_r*(4.0-Liver_fib);
dxdt_Liver_fib = (Liver_fib<4.0)?fib_r:0.0;
double IL8_p = k_PMN_base+0.3*NE_free;
dxdt_IL8_lung = IL8_p-0.30*IL8_lung;
dxdt_PMN_lung = k_PMN_base+k_IL8_recruit*IL8_lung-k_PMN_egress*PMN_lung;
double ELF_uM = AAT_s*r_ELF_plasma*0.19;
double AAT_inhib = ELF_uM/(11.0+ELF_uM);
double NE_inhib_t = AAT_inhib+(1.0-AAT_inhib)*NEi_eff;
dxdt_NE_free = k_NE_release*PMN_lung*(1.0-NE_inhib_t)-k_NE_elim*NE_free;
dxdt_MMP12_lung = k_MMP12_base+k_MMP12_NE*NE_free-k_MMP12_elim*MMP12_lung;
double ed_r = (k_elastin_deg*NE_free+k_elastin_deg*0.5*MMP12_lung)*(Elastin/Elastin0);
double ep_r = k_elastin_syn*(Elastin0-Elastin);
dxdt_Elastin = -ed_r+ep_r;
dxdt_FEV1_pct = -k_FEV1_Edecay*ed_r*0.8-annual_decline_base/365.0;

$TABLE
double AAT_uM = ((AAT_C1>0)?AAT_C1:0.01)*0.19;
double ELF_tab = AAT_uM*0.10;
double EmphIdx = 100.0*(1.0-Elastin/Elastin0);
double FEV1_o = FEV1_pct>10?FEV1_pct:10;
double Exacer = 0.5*exp(-0.04*(FEV1_o-30));
double NE_inhib_pct = 100.0*ELF_tab/(11.0+ELF_tab);
double SGRQ = 100.0-FEV1_o*0.65-(100.0-EmphIdx)*0.20;
capture AAT_uMol=AAT_uM; capture ELF_uMol=ELF_tab;
capture EmphIndex=EmphIdx; capture FEV1=FEV1_o;
capture Exacer_yr=Exacer; capture NE_inhib=NE_inhib_pct;
capture SGRQ_score=SGRQ; capture ZPolymer=ZAAT_Poly;
capture NEfree=NE_free; capture LiverFib=Liver_fib;
capture HSC=HSC_act; capture MMP12=MMP12_lung;

$CAPTURE AAT_uMol ELF_uMol EmphIndex FEV1 Exacer_yr NE_inhib
         SGRQ_score ZPolymer NEfree LiverFib HSC MMP12
'

mod_global <- mcode("AATD_shiny", aatd_code)

## ── Colour Palette ────────────────────────────────────────────
scen_pal <- c(
  "Untreated" = "#E53935",
  "Augmentation" = "#1976D2",
  "siRNA (Fazirsiran)" = "#00897B",
  "NE Inhibitor" = "#F57F17",
  "Gene Therapy" = "#6A1B9A",
  "Aug + NE Inh" = "#2E7D32"
)

## ── Helper: run simulation given input parameters ─────────────
run_sim <- function(mod, sim_yrs, use_aug, aug_dose_amt, use_sirna, use_nei, nei_dose_amt, use_gene) {
  sim_days <- sim_yrs * 365
  events   <- NULL

  if (use_aug) {
    t_aug <- seq(0, sim_days, by = 7)
    e_aug <- ev(time = t_aug, amt = aug_dose_amt, cmt = "AUG_C1", rate = -2)
    events <- if (is.null(events)) e_aug else events + e_aug
  }
  if (use_sirna) {
    t_si <- seq(0, sim_days, by = 84)
    e_si <- ev(time = t_si, amt = 1.0, cmt = "siRNA_Eff")
    events <- if (is.null(events)) e_si else events + e_si
  }
  if (use_nei) {
    t_ni <- seq(0, sim_days, by = 0.5)
    e_ni <- ev(time = t_ni, amt = nei_dose_amt, cmt = "NEi_A")
    events <- if (is.null(events)) e_ni else events + e_ni
  }
  if (use_gene) {
    e_ge <- ev(time = 30, amt = 2.0, cmt = "Gene_Eff")
    events <- if (is.null(events)) e_ge else events + e_ge
  }

  if (is.null(events)) {
    out <- mrgsim(mod, end = sim_days, delta = 7)
  } else {
    out <- mrgsim(mod, events = events, end = sim_days, delta = 7)
  }
  as_tibble(out) %>% mutate(Year = time / 365)
}

## ── UI ────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(
    title = tags$div(
      style = "font-size:15px; font-weight:bold;",
      "AATD QSP Dashboard"
    ), titleWidth = 250
  ),
  dashboardSidebar(
    width = 240,
    sidebarMenu(
      id = "tabs",
      menuItem("1. Patient Profile", tabName = "profile", icon = icon("user")),
      menuItem("2. Drug PK / AAT Levels", tabName = "pk", icon = icon("syringe")),
      menuItem("3. Pulmonary PD", tabName = "pd_lung", icon = icon("lungs")),
      menuItem("4. Clinical Endpoints", tabName = "endpoints", icon = icon("chart-line")),
      menuItem("5. Scenario Comparison", tabName = "scenarios", icon = icon("layer-group")),
      menuItem("6. Biomarkers", tabName = "biomarkers", icon = icon("flask"))
    ),
    hr(),
    h5("Patient Parameters", style = "padding-left:15px; color:#aaa;"),
    sliderInput("FEV1_base", "Baseline FEV1 (% pred)", 40, 95, 75, step = 5),
    sliderInput("sim_yrs",  "Simulation (years)", 1, 10, 5, step = 1),
    hr(),
    h5("Treatment", style = "padding-left:15px; color:#aaa;"),
    checkboxInput("use_aug",  "AAT Augmentation (IV)", value = FALSE),
    conditionalPanel("input.use_aug",
      sliderInput("aug_dose", "Dose (mg/dL eq.)", 30, 80, 55, step = 5)
    ),
    checkboxInput("use_sirna", "siRNA (Fazirsiran)", value = FALSE),
    checkboxInput("use_nei",   "NE Inhibitor",       value = FALSE),
    conditionalPanel("input.use_nei",
      sliderInput("nei_dose", "NE Inh Dose (rel.)", 0.1, 0.5, 0.30, step = 0.05)
    ),
    checkboxInput("use_gene", "Gene Therapy",       value = FALSE),
    hr(),
    actionButton("run_btn", "Run Simulation", class = "btn-primary", width = "90%")
  ),

  dashboardBody(
    tags$head(tags$style(HTML(
      ".content-wrapper { background-color: #f8f9fa; }
       .box { border-radius: 8px; }
       .ggplot-box { background: white; border-radius: 8px; padding: 10px; margin-bottom: 10px; }"
    ))),
    tabItems(

      ## ── Tab 1: Patient Profile ─────────────────────────────
      tabItem(tabName = "profile",
        fluidRow(
          box(width = 4, title = "Genotype & Disease Profile", status = "primary", solidHeader = TRUE,
            selectInput("genotype", "Genotype", c("ZZ (Severe)", "SZ (Moderate)", "MZ (Carrier)", "MM (Normal)")),
            selectInput("smoking",  "Smoking Status", c("Never", "Ex-Smoker", "Current Smoker")),
            numericInput("age",    "Age (years)",   45, 20, 80),
            numericInput("weight", "Weight (kg)",   70, 40, 120),
            selectInput("phenotype", "Disease Phenotype", c("Lung Dominant", "Liver Dominant", "Combined"))
          ),
          box(width = 4, title = "Clinical Measurements", status = "warning", solidHeader = TRUE,
            numericInput("AAT_init",  "Baseline AAT (mg/dL)", 10, 2, 200),
            numericInput("FEV1_meas", "Measured FEV1 (% pred)", 75, 20, 100),
            numericInput("CT_emph",   "CT Emphysema Index (%)", 15, 0, 60),
            numericInput("Metavir_in","Liver Fibrosis (Metavir)", 1, 0, 4)
          ),
          box(width = 4, title = "Risk Summary", status = "danger", solidHeader = TRUE,
            h4("Estimated Annual FEV1 Decline"),
            textOutput("fev1_decline_est"),
            h4("Exacerbation Risk (Annual)"),
            textOutput("exacer_risk_est"),
            h4("10-yr HCC Risk (ZZ Genotype)"),
            textOutput("hcc_risk_est"),
            h4("Recommendation"),
            textOutput("reco_txt")
          )
        ),
        fluidRow(
          box(width = 12, title = "AATD Disease Overview", status = "info", solidHeader = TRUE,
            fluidRow(
              column(6,
                h4("Disease Mechanism"),
                p("Alpha-1 Antitrypsin Deficiency (AATD) is caused by mutations in the SERPINA1 gene.
                  The most severe form (ZZ genotype, Glu342Lys) leads to misfolding of the AAT protein,
                  causing 85-90% retention in hepatocyte ER as pathological loop-sheet polymers."),
                p("DUAL PATHOLOGY:"),
                tags$ul(
                  tags$li("LIVER: Gain-of-function toxicity — Z-AAT polymer accumulation → ER stress →
                    hepatocyte apoptosis → fibrosis → cirrhosis → HCC risk 10-40× ↑"),
                  tags$li("LUNG: Loss-of-function — Insufficient serum AAT (<11 µM) → unopposed neutrophil
                    elastase → panacinar emphysema → severe COPD (FEV1 decline 2-3× normal rate)")
                )
              ),
              column(6,
                h4("Key Pharmacological Targets"),
                tableOutput("pharm_table")
              )
            )
          )
        )
      ),

      ## ── Tab 2: Drug PK / AAT Levels ───────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(width = 12, title = "Serum AAT Time Course", status = "primary", solidHeader = TRUE,
              withSpinner(plotlyOutput("pk_plot_serum", height = 380))
          )
        ),
        fluidRow(
          box(width = 6, title = "ELF AAT Level (Lung Compartment)", status = "success", solidHeader = TRUE,
              withSpinner(plotlyOutput("pk_plot_elf", height = 280))
          ),
          box(width = 6, title = "NE Inhibition in ELF (%)", status = "warning", solidHeader = TRUE,
              withSpinner(plotlyOutput("pk_plot_inhib", height = 280))
          )
        ),
        fluidRow(
          box(width = 12, title = "PK Model Summary", status = "info", solidHeader = TRUE,
            tableOutput("pk_params_table")
          )
        )
      ),

      ## ── Tab 3: Pulmonary PD ───────────────────────────────
      tabItem(tabName = "pd_lung",
        fluidRow(
          box(width = 6, title = "Free Neutrophil Elastase (Lung)", status = "danger", solidHeader = TRUE,
              withSpinner(plotlyOutput("pd_NE", height = 300))
          ),
          box(width = 6, title = "Lung Elastin Content (%)", status = "warning", solidHeader = TRUE,
              withSpinner(plotlyOutput("pd_Elastin", height = 300))
          )
        ),
        fluidRow(
          box(width = 6, title = "MMP-12 (Macrophage Elastase)", status = "warning", solidHeader = TRUE,
              withSpinner(plotlyOutput("pd_MMP12", height = 280))
          ),
          box(width = 6, title = "Lung PMN Burden", status = "primary", solidHeader = TRUE,
              withSpinner(plotlyOutput("pd_PMN", height = 280))
          )
        )
      ),

      ## ── Tab 4: Clinical Endpoints ─────────────────────────
      tabItem(tabName = "endpoints",
        fluidRow(
          box(width = 6, title = "FEV1 (% Predicted)", status = "danger", solidHeader = TRUE,
              withSpinner(plotlyOutput("ep_FEV1", height = 320))
          ),
          box(width = 6, title = "CT Emphysema Index (%)", status = "warning", solidHeader = TRUE,
              withSpinner(plotlyOutput("ep_CT", height = 320))
          )
        ),
        fluidRow(
          box(width = 6, title = "SGRQ Score (Quality of Life)", status = "info", solidHeader = TRUE,
              withSpinner(plotlyOutput("ep_SGRQ", height = 280))
          ),
          box(width = 6, title = "Annual Exacerbation Risk", status = "danger", solidHeader = TRUE,
              withSpinner(plotlyOutput("ep_Exacer", height = 280))
          )
        )
      ),

      ## ── Tab 5: Scenario Comparison ────────────────────────
      tabItem(tabName = "scenarios",
        fluidRow(
          box(width = 12, title = "All 6 Treatment Scenarios — Comparative Simulation", status = "primary", solidHeader = TRUE,
            p("Comparative simulation of all 6 treatment strategies vs untreated natural history."),
            withSpinner(plotlyOutput("sc_FEV1", height = 350))
          )
        ),
        fluidRow(
          box(width = 6, title = "Serum AAT by Scenario", status = "success", solidHeader = TRUE,
              withSpinner(plotlyOutput("sc_AAT", height = 280))
          ),
          box(width = 6, title = "Liver Z-Polymer by Scenario", status = "warning", solidHeader = TRUE,
              withSpinner(plotlyOutput("sc_Poly", height = 280))
          )
        ),
        fluidRow(
          box(width = 12, title = "5-Year Outcome Summary Table", status = "info", solidHeader = TRUE,
              withSpinner(DTOutput("sc_table"))
          )
        )
      ),

      ## ── Tab 6: Biomarkers ─────────────────────────────────
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(width = 6, title = "Hepatic Z-AAT Polymer (siRNA Biomarker)", status = "danger", solidHeader = TRUE,
              withSpinner(plotlyOutput("bm_poly", height = 280))
          ),
          box(width = 6, title = "Hepatic Fibrosis Stage (Metavir)", status = "warning", solidHeader = TRUE,
              withSpinner(plotlyOutput("bm_fib", height = 280))
          )
        ),
        fluidRow(
          box(width = 6, title = "Liver Stellate Cell Activation (HSC)", status = "primary", solidHeader = TRUE,
              withSpinner(plotlyOutput("bm_HSC", height = 280))
          ),
          box(width = 6, title = "Desmosine (Elastin Degradation) — Derived", status = "success", solidHeader = TRUE,
              withSpinner(plotlyOutput("bm_Desmo", height = 280))
          )
        ),
        fluidRow(
          box(width = 12, title = "Biomarker Reference Ranges", status = "info", solidHeader = TRUE,
              tableOutput("bm_ref_table")
          )
        )
      )
    )
  )
)

## ── Server ───────────────────────────────────────────────────
server <- function(input, output, session) {

  ## ── Reactive simulation ───────────────────────────────────
  sim_data <- eventReactive(input$run_btn, {
    run_sim(mod_global,
            sim_yrs    = input$sim_yrs,
            use_aug    = input$use_aug,
            aug_dose_amt = input$aug_dose,
            use_sirna  = input$use_sirna,
            use_nei    = input$use_nei,
            nei_dose_amt = input$nei_dose,
            use_gene   = input$use_gene)
  }, ignoreNULL = FALSE)

  ## ── All-scenarios reactive ────────────────────────────────
  scen_data <- reactive({
    sim_days <- input$sim_yrs * 365
    run_scenario <- function(use_aug, use_sirna, use_nei, use_gene, label) {
      run_sim(mod_global, input$sim_yrs, use_aug, 55, use_sirna, use_nei, 0.30, use_gene) %>%
        mutate(Scenario = label)
    }
    bind_rows(
      run_scenario(FALSE, FALSE, FALSE, FALSE, "Untreated"),
      run_scenario(TRUE,  FALSE, FALSE, FALSE, "Augmentation"),
      run_scenario(FALSE, TRUE,  FALSE, FALSE, "siRNA (Fazirsiran)"),
      run_scenario(FALSE, FALSE, TRUE,  FALSE, "NE Inhibitor"),
      run_scenario(FALSE, FALSE, FALSE, TRUE,  "Gene Therapy"),
      run_scenario(TRUE,  FALSE, TRUE,  FALSE, "Aug + NE Inh")
    )
  })

  ## ── Tab 1: Patient profile outputs ───────────────────────
  output$fev1_decline_est <- renderText({
    rate <- switch(input$smoking,
      "Never"          = "50-80 mL/year",
      "Ex-Smoker"      = "80-120 mL/year",
      "Current Smoker" = "120-200 mL/year"
    )
    paste0(rate, " (ZZ genotype, ", input$smoking, ")")
  })
  output$exacer_risk_est <- renderText({
    risk <- round(0.5 * exp(-0.04 * (input$FEV1_base - 30)), 2)
    paste0(risk, " exacerbations/year")
  })
  output$hcc_risk_est <- renderText({
    if (input$genotype == "ZZ (Severe)") "~28% lifetime risk (ZZ, male)" else "<1% (non-ZZ)"
  })
  output$reco_txt <- renderText({
    if (input$genotype == "ZZ (Severe)" && input$FEV1_base < 80)
      "Strongly consider AAT augmentation therapy per ERS/ATS guidelines"
    else if (input$genotype == "ZZ (Severe)")
      "Monitor FEV1 annually; consider augmentation if FEV1 decline ≥100 mL/yr"
    else
      "Carrier state; optimize modifiable risk factors"
  })
  output$pharm_table <- renderTable({
    data.frame(
      Target = c("NE (Neutrophil Elastase)", "Liver Z-AAT Polymer", "SERPINA1 mRNA",
                 "SERPINA1 Gene", "Airway Smooth Muscle", "Airway Inflammation"),
      Drug   = c("Alvelestat, Lonodelestat, Brensocatib", "VX-864, GSK3117391 (Correctors)",
                 "Fazirsiran, Belcesiran (siRNA/ASO)", "rAAV-SERPINA1, CRISPR/Cas9 (Gene Therapy)",
                 "LABA/LAMA (Bronchodilators)", "ICS, Roflumilast"),
      Status = c("Phase 2", "Phase 1-2", "Phase 2-3", "Phase 1-2", "Approved (COPD)", "Approved (COPD)")
    )
  })

  ## ── Tab 2: PK Plots ───────────────────────────────────────
  output$pk_plot_serum <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(Year, AAT_uMol)) +
      geom_line(colour = "#1976D2", linewidth = 1.2) +
      geom_hline(yintercept = 11, linetype = "dashed", colour = "red", linewidth = 0.8) +
      annotate("text", x = max(df$Year) * 0.7, y = 12.5,
               label = "Protective threshold (11 µM)", colour = "red", size = 3.5) +
      labs(title = "Serum AAT (µM)", x = "Year", y = "AAT (µM)") +
      theme_minimal(base_size = 12)
    ggplotly(p, tooltip = c("x", "y"))
  })

  output$pk_plot_elf <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(Year, ELF_uMol)) +
      geom_line(colour = "#2E7D32", linewidth = 1.1) +
      geom_hline(yintercept = 1.1, linetype = "dashed", colour = "red") +
      labs(title = "ELF AAT (µM)", x = "Year", y = "ELF AAT (µM)") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$pk_plot_inhib <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(Year, NE_inhib)) +
      geom_line(colour = "#FF8F00", linewidth = 1.1) +
      ylim(0, 100) +
      labs(title = "NE Inhibition in ELF (%)", x = "Year", y = "Inhibition (%)") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$pk_params_table <- renderTable({
    data.frame(
      Parameter = c("Elimination rate (k10)", "V1 (Central)", "k12 (Distrib.)", "k21 (Redistrib.)",
                    "ELF/Plasma ratio", "Half-life (t½)", "Protective threshold"),
      Value     = c("0.052 /day", "3.76 L/kg", "0.032 /day", "0.022 /day",
                    "0.10", "~4.5 days", ">11 µM (57 mg/dL)"),
      Source    = c("Prolastin-C label", "McElvaney 1997", "RAPID trial PK",
                    "RAPID trial PK", "Hubbard 1991", "Crystal 1990", "RAPID trial 2015")
    )
  })

  ## ── Tab 3: Pulmonary PD Plots ─────────────────────────────
  output$pd_NE <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(Year, NEfree)) +
      geom_line(colour = "#C62828", linewidth = 1.1) +
      geom_hline(yintercept = 0.15, linetype = "dashed", colour = "grey50") +
      annotate("text", x = 0.5, y = 0.17, label = "Normal NE (~0.15)", size = 3, colour = "grey50") +
      labs(title = "Free NE (Lung, normalized)", x = "Year", y = "NE Activity (rel. units)") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$pd_Elastin <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(Year, Elastin)) +
      geom_line(colour = "#4CAF50", linewidth = 1.1) +
      geom_hline(yintercept = 100, linetype = "dashed", colour = "grey50") +
      ylim(0, 105) +
      labs(title = "Lung Elastin Content (%)", x = "Year", y = "Elastin (% of normal)") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$pd_MMP12 <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(Year, MMP12)) +
      geom_line(colour = "#FF8F00", linewidth = 1) +
      labs(title = "MMP-12 Macrophage Elastase", x = "Year", y = "MMP-12 (normalized)") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$pd_PMN <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(Year, PMN_lung)) +
      geom_line(colour = "#7B1FA2", linewidth = 1) +
      labs(title = "Lung PMN Burden", x = "Year", y = "PMN (normalized)") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  ## ── Tab 4: Clinical Endpoint Plots ───────────────────────
  output$ep_FEV1 <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(Year, FEV1)) +
      geom_line(colour = "#1565C0", linewidth = 1.2) +
      geom_hline(yintercept = c(80, 50, 30), linetype = "dotted", colour = "grey60") +
      ylim(0, 100) +
      labs(title = "FEV1 (% Predicted)", x = "Year", y = "FEV1 %") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$ep_CT <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(Year, EmphIndex)) +
      geom_line(colour = "#BF360C", linewidth = 1.2) +
      labs(title = "CT Emphysema Index (%)", x = "Year", y = "Emphysema (%)") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$ep_SGRQ <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(Year, SGRQ_score)) +
      geom_line(colour = "#00695C", linewidth = 1.1) +
      ylim(0, 100) +
      labs(title = "SGRQ Score (QoL)", x = "Year", y = "SGRQ (0-100; ↑ = worse)") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$ep_Exacer <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(Year, Exacer_yr)) +
      geom_line(colour = "#AD1457", linewidth = 1.1) +
      labs(title = "Annual Exacerbation Risk", x = "Year", y = "Exacerbations/year") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  ## ── Tab 5: Scenario Comparison Plots ─────────────────────
  output$sc_FEV1 <- renderPlotly({
    df <- scen_data()
    p <- ggplot(df, aes(Year, FEV1, colour = Scenario)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = scen_pal) +
      geom_hline(yintercept = c(80, 50, 30), linetype = "dotted", colour = "grey60") +
      ylim(0, 100) +
      labs(title = "FEV1 (%) — All Scenarios", x = "Year", y = "FEV1 %") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    ggplotly(p) %>% layout(legend = list(orientation = "h", y = -0.25))
  })

  output$sc_AAT <- renderPlotly({
    df <- scen_data()
    p <- ggplot(df, aes(Year, AAT_uMol, colour = Scenario)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = scen_pal) +
      geom_hline(yintercept = 11, linetype = "dashed") +
      labs(title = "Serum AAT (µM)", x = "Year", y = "AAT µM") +
      theme_minimal(base_size = 12) + theme(legend.position = "none")
    ggplotly(p)
  })

  output$sc_Poly <- renderPlotly({
    df <- scen_data()
    p <- ggplot(df, aes(Year, ZPolymer, colour = Scenario)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = scen_pal) +
      labs(title = "Hepatic Z-Polymer", x = "Year", y = "Z-Polymer (rel.)") +
      theme_minimal(base_size = 12) + theme(legend.position = "none")
    ggplotly(p)
  })

  output$sc_table <- renderDT({
    df <- scen_data()
    tbl <- df %>%
      filter(time == max(time)) %>%
      mutate(
        `AAT (µM)` = round(AAT_uMol, 1),
        `FEV1 (%)` = round(FEV1, 1),
        `Emph. Index (%)` = round(EmphIndex, 1),
        `Z-Polymer` = round(ZPolymer, 2),
        `Metavir` = round(LiverFib, 2),
        `Exacer/yr` = round(Exacer_yr, 2),
        `SGRQ` = round(SGRQ_score, 1)
      ) %>%
      select(Scenario, `AAT (µM)`, `FEV1 (%)`, `Emph. Index (%)`,
             `Z-Polymer`, Metavir, `Exacer/yr`, SGRQ)
    datatable(tbl, options = list(dom = "t", paging = FALSE), rownames = FALSE) %>%
      formatStyle("FEV1 (%)", background = styleColorBar(c(0, 100), "#AED6F1"))
  })

  ## ── Tab 6: Biomarker Plots ────────────────────────────────
  output$bm_poly <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(Year, ZPolymer)) +
      geom_line(colour = "#E53935", linewidth = 1.1) +
      labs(title = "Hepatic Z-AAT Polymer", x = "Year", y = "Polymer (rel. units)") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$bm_fib <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(Year, LiverFib)) +
      geom_line(colour = "#795548", linewidth = 1.1) +
      scale_y_continuous(breaks = 0:4, limits = c(0, 4)) +
      labs(title = "Metavir Fibrosis Stage", x = "Year", y = "Metavir (0-4)") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$bm_HSC <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(Year, HSC)) +
      geom_line(colour = "#1976D2", linewidth = 1.1) +
      labs(title = "Activated HSC (Relative)", x = "Year", y = "HSC Activation") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$bm_Desmo <- renderPlotly({
    df <- sim_data()
    desmosine <- df %>%
      mutate(Desmosine = (100 - Elastin) * 2.5 + rnorm(nrow(df), 0, 0.5))
    p <- ggplot(desmosine, aes(Year, Desmosine)) +
      geom_line(colour = "#2E7D32", linewidth = 1.1) +
      labs(title = "Plasma Desmosine (Derived, nmol/L)", x = "Year", y = "Desmosine (nmol/L)") +
      theme_minimal(base_size = 12)
    ggplotly(p)
  })

  output$bm_ref_table <- renderTable({
    data.frame(
      Biomarker = c("Serum AAT (ZZ)", "Serum AAT (MM)", "ELF AAT (ZZ)", "ELF AAT (MM)",
                    "BAL NE Activity", "Plasma Desmosine", "FIB-4 Score",
                    "Liver Stiffness (kPa)"),
      `Reference Range` = c("6-8 mg/dL (2-7 µM)", "150-300 mg/dL (20-53 µM)",
                             "0.2-0.7 µM (ZZ)", "2-5 µM (MM)", "↑↑ in AATD (>400 µg/mL/hr)",
                             "4-7 nmol/L (normal)", "<1.45 (No significant fibrosis)",
                             "<7 kPa (normal liver)"),
      Significance = c("Key deficiency marker", "Normal range", "Lung antiprotease",
                       "Lung antiprotease reserve", "Proteolytic burden",
                       "Elastin degradation product", "Liver fibrosis screening",
                       "F3-F4 fibrosis if >12.5 kPa")
    )
  })
}

shinyApp(ui, server)
