## ============================================================
## IgG4-Related Disease (IgG4-RD) — Interactive Shiny Dashboard
## QSP Simulation: Rituximab / Prednisone / Dupilumab
## 7 Tabs: Overview · Patient Profile · PK · B Cell/Immuno ·
##         Cytokines/Fibrosis · Scenario Comparison · Biomarkers
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)

## ── mrgsolve model (inline for Shiny) ───────────────────────────────────────
model_code <- '
$PARAM
CL_RTX=0.35, V1_RTX=3.0, Q_RTX=0.8, V2_RTX=4.5
kon_RTX=0.55, koff_RTX=0.002, CD20_ss=180, ksyn_CD20=90, kdeg_CD20=0.5, kint_RTX=0.3
ka_PRED=3.0, CL_PRED=18.0, V_PRED=40.0, F_PRED=0.80
ka_DUP=0.08, CL_DUP=0.25, V1_DUP=4.5
kon_DUP=0.30, koff_DUP=0.001, IL4RA_ss=10, ksyn_IL4RA=0.5, kdeg_IL4RA=0.05, kint_DUP=0.15
kprolif_BNV=0.05, kdeath_BNV=0.05, kgc_BNV=0.02
kprolif_GCB=0.20, kdeath_GCB=0.18, kpb_GCB=0.10
kmat_PB=0.30, kdeath_PB=0.40, kdeath_PC=0.005
kprolif_TFH2=0.08, kdeath_TFH2=0.05, IL4_EC50_TFH2=0.5
kprolif_CTL4=0.04, kdeath_CTL4=0.04
kprod_IL4=0.10, kdeg_IL4=1.0
kprod_IL10=0.06, kdeg_IL10=0.8
kprod_TGFB=0.08, kdeg_TGFB=0.5
ksec_IgG4=0.15, kdeg_IgG4=0.015, IgG4_0=450
kact_MYOFIB=0.03, kinact_MYOFIB=0.02, TGFB_EC50_FIB=1.5
kcol_ECM=0.04, kdeg_ECM=0.005
Emax_RTX=0.98, EC50_RTX=50.0
Emax_PRED=0.85, EC50_PRED=80.0
w_IgG4=0.30, w_ECM=0.35, w_PC=0.20, w_TFH2=0.15

$INIT
CENT_RTX=0, PERI_RTX=0, CD20_FREE=180, RTX_CD20=0
GUT_PRED=0, CENT_PRED=0
SC_DUP=0, CENT_DUP=0, IL4RA_FREE=10, DUP_IL4RA=0
BNV=1.0, GCB=2.5, PB=3.0, PC=4.0, TFH2=3.5, CTL4=2.0
IgG4_SER=450, IL4=3.0, IL10=2.5, TGFB=3.0
MYOFIB=2.0, ECM=2.5, IRI=8.0

$ODE
double kon_r=kon_RTX*CENT_RTX*CD20_FREE;
double koff_r=koff_RTX*RTX_CD20;
double kint_r=kint_RTX*RTX_CD20;
dxdt_CENT_RTX=-CL_RTX/V1_RTX*CENT_RTX-Q_RTX/V1_RTX*CENT_RTX+Q_RTX/V2_RTX*PERI_RTX-kon_r+koff_r;
dxdt_PERI_RTX=Q_RTX/V1_RTX*CENT_RTX-Q_RTX/V2_RTX*PERI_RTX;
dxdt_CD20_FREE=ksyn_CD20-kdeg_CD20*CD20_FREE-kon_r+koff_r;
dxdt_RTX_CD20=kon_r-koff_r-kint_r;
dxdt_GUT_PRED=-ka_PRED*GUT_PRED;
dxdt_CENT_PRED=ka_PRED*GUT_PRED*F_PRED-CL_PRED/V_PRED*CENT_PRED;
double kon_d=kon_DUP*CENT_DUP*IL4RA_FREE;
double koff_d=koff_DUP*DUP_IL4RA;
double kint_d=kint_DUP*DUP_IL4RA;
dxdt_SC_DUP=-ka_DUP*SC_DUP;
dxdt_CENT_DUP=ka_DUP*SC_DUP-CL_DUP/V1_DUP*CENT_DUP-kon_d+koff_d;
dxdt_IL4RA_FREE=ksyn_IL4RA-kdeg_IL4RA*IL4RA_FREE-kon_d+koff_d;
dxdt_DUP_IL4RA=kon_d-koff_d-kint_d;
double CD20_occ=RTX_CD20/(RTX_CD20+CD20_FREE+0.001);
double Ekill_RTX=Emax_RTX*CD20_occ;
double Eimmu_PRED=Emax_PRED*CENT_PRED/(EC50_PRED+CENT_PRED);
double IL4RA_occ=DUP_IL4RA/(DUP_IL4RA+IL4RA_FREE+0.001);
double TFH2_stim=TFH2/(1.0+TFH2);
dxdt_BNV=kprolif_BNV*BNV-kdeath_BNV*BNV-kgc_BNV*BNV-Ekill_RTX*kprolif_BNV*BNV;
dxdt_GCB=kgc_BNV*BNV+kprolif_GCB*GCB*TFH2_stim*(1-IL4RA_occ*0.5)-kdeath_GCB*GCB-kpb_GCB*GCB-Ekill_RTX*kprolif_GCB*GCB-Eimmu_PRED*0.8*GCB;
dxdt_PB=kpb_GCB*GCB-kmat_PB*PB-kdeath_PB*PB-Ekill_RTX*(kmat_PB+kdeath_PB)*PB-Eimmu_PRED*0.9*PB;
double PC_RTX_kill=Ekill_RTX*0.15;
dxdt_PC=kmat_PB*PB-kdeath_PC*PC-PC_RTX_kill*PC-Eimmu_PRED*0.2*PC;
double IL4_drive=IL4/(IL4_EC50_TFH2+IL4);
dxdt_TFH2=kprolif_TFH2*TFH2*IL4_drive*(1-IL4RA_occ)-kdeath_TFH2*TFH2-Eimmu_PRED*0.85*TFH2;
dxdt_CTL4=kprolif_CTL4*CTL4-kdeath_CTL4*CTL4-Eimmu_PRED*0.7*CTL4;
dxdt_IL4=kprod_IL4*TFH2*(1-IL4RA_occ*0.6)-kdeg_IL4*IL4-Eimmu_PRED*0.6*IL4;
dxdt_IL10=kprod_IL10*(TFH2+0.5)-kdeg_IL10*IL10-Eimmu_PRED*0.3*IL10;
dxdt_TGFB=kprod_TGFB*CTL4*(1-IL4RA_occ*0.3)-kdeg_TGFB*TGFB-Eimmu_PRED*0.4*TGFB;
double IgG4_switch=IL4*IL10/(1+IL4*IL10);
dxdt_IgG4_SER=ksec_IgG4*PC*IgG4_switch-kdeg_IgG4*IgG4_SER;
double TGFB_eff=TGFB*TGFB/(TGFB_EC50_FIB*TGFB_EC50_FIB+TGFB*TGFB);
dxdt_MYOFIB=kact_MYOFIB*TGFB_eff-kinact_MYOFIB*MYOFIB-Eimmu_PRED*0.35*MYOFIB;
dxdt_ECM=kcol_ECM*MYOFIB-kdeg_ECM*ECM-Eimmu_PRED*0.2*ECM;
double IRI_new=(IgG4_SER/IgG4_0)*10*w_IgG4+(ECM/2.5)*10*w_ECM+(PC/4.0)*10*w_PC+(TFH2/3.5)*10*w_TFH2;
if(IRI_new<0) IRI_new=0; if(IRI_new>24) IRI_new=24;
dxdt_IRI=(IRI_new-IRI)*0.5;

$TABLE
double RTX_conc=CENT_RTX*148000/1e6*1000;
double Bcell_pct=BNV*100;
double CD20_occ_pct=RTX_CD20/(RTX_CD20+CD20_FREE+0.001)*100;
double IL4RA_occ_pct=DUP_IL4RA/(DUP_IL4RA+IL4RA_FREE+0.001)*100;
double CR_flag=(IRI<1.0)?1.0:0.0;
double PR_flag=(IRI>=1.0&&IRI<4.0)?1.0:0.0;

$CAPTURE RTX_conc Bcell_pct CD20_occ_pct IL4RA_occ_pct
         IgG4_SER CENT_PRED TFH2 PC TGFB IL4 IL10 ECM MYOFIB IRI CR_flag PR_flag
'

mod_global <- mcode("IgG4RD_shiny", model_code)

## ── Simulation function ────────────────────────────────────────────────────
run_sim <- function(scenario, rtx_dose_mg = 1000, rtx_n = 2,
                    pred_dose_mg = 40, pred_dur_wk = 24,
                    dup_dose_mg = 300, sim_months = 24,
                    baseline_IgG4 = 450, disease_severity = 1.0) {
  end_day <- sim_months * 30.4
  MW_RTX <- 148; MW_PRED <- 358.4 / 1e6; MW_DUP <- 146

  events <- list()

  if (scenario %in% c("RTX_1g2", "RTX_maint", "RTX_375x4")) {
    if (scenario == "RTX_1g2") {
      rtx_times <- c(0, 14)[1:min(rtx_n, 2)]
    } else if (scenario == "RTX_375x4") {
      rtx_times <- c(0, 7, 14, 21)
    } else {  # maintenance
      rtx_times <- c(0, 14, 180, 360)
    }
    rtx_nM <- (rtx_dose_mg / (MW_RTX * 1000)) * 1e9 / 3.0
    events <- c(events, list(ev(time = rtx_times, cmt = "CENT_RTX", amt = rtx_nM)))
  }

  if (scenario %in% c("Prednisone", "RTX_plus_GC")) {
    daily_nM <- (pred_dose_mg * 1e6 / MW_PRED) / 40.0
    pred_days <- seq(0, min(pred_dur_wk * 7, end_day - 1))
    events <- c(events, list(ev(time = pred_days, cmt = "GUT_PRED", amt = daily_nM)))
  }

  if (scenario == "Dupilumab") {
    dup_nM <- (dup_dose_mg / (MW_DUP * 1000)) * 1e9 / 4.5
    dup_times <- seq(0, end_day - 1, by = 14)
    events <- c(events, list(ev(time = dup_times, cmt = "SC_DUP", amt = dup_nM)))
  }

  if (length(events) == 0) events <- list(ev(time = 0, cmt = "CENT_RTX", amt = 0))

  all_ev <- do.call(c, events)

  # Adjust initial conditions for severity/IgG4 level
  mod_adj <- mod_global %>%
    param(IgG4_0 = baseline_IgG4) %>%
    init(IgG4_SER = baseline_IgG4,
         TFH2 = 3.5 * disease_severity,
         PC = 4.0 * disease_severity,
         ECM = 2.5 * disease_severity,
         IRI = 8.0 * disease_severity)

  tryCatch({
    mod_adj %>%
      mrgsim(events = all_ev, end = end_day, delta = 2) %>%
      as.data.frame() %>%
      mutate(time_months = time / 30.4,
             Scenario = scenario)
  }, error = function(e) {
    data.frame(time = 0, time_months = 0, IgG4_SER = baseline_IgG4,
               IRI = 8, Bcell_pct = 100, Scenario = scenario,
               error = TRUE)
  })
}

## ── UI ─────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "IgG4-RD QSP Model"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview",          tabName = "overview",   icon = icon("info-circle")),
      menuItem("Patient Profile",   tabName = "patient",    icon = icon("user")),
      menuItem("Pharmacokinetics",  tabName = "pk",         icon = icon("prescription-bottle")),
      menuItem("B Cell & Immunity", tabName = "bcell",      icon = icon("microscope")),
      menuItem("Cytokines & Fibrosis", tabName = "cytokine", icon = icon("dna")),
      menuItem("Scenario Comparison", tabName = "scenario",  icon = icon("chart-bar")),
      menuItem("Biomarkers",         tabName = "biomarker",  icon = icon("vial"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .main-header .logo { font-weight: bold; font-size: 18px; }
      .value-box { border-radius: 8px; }
      .box { border-radius: 8px; }
      .content-wrapper { background: #f8f9fa; }
    "))),
    tabItems(

      ## ── TAB 1: Overview ────────────────────────────────────────────────
      tabItem("overview",
        fluidRow(
          box(width = 12, status = "primary", solidHeader = TRUE,
              title = "IgG4-Related Disease (IgG4-RD) — QSP Overview",
              HTML('
                <p><b>IgG4-Related Disease (IgG4-RD)</b> is a systemic fibroinflammatory
                condition characterized by tumefactive lesions with dense lymphoplasmacytic
                infiltrate rich in IgG4+ plasma cells, storiform fibrosis, and obliterative
                phlebitis. It can affect virtually any organ.</p>
                <h4>Key Pathogenic Mechanisms</h4>
                <ul>
                  <li><b>Tfh2 cells</b> (PD-1+CXCR5+): Expanded, drive IgG4 class switching via IL-4/IL-21</li>
                  <li><b>Cytotoxic CD4+ T cells (SLAMF7+)</b>: Tissue-homing, produce TGF-β → fibrosis</li>
                  <li><b>IgG4+ plasmablasts/plasma cells</b>: Elevated; IgG4 = diagnostic biomarker</li>
                  <li><b>M2 macrophages</b>: CCL18, TGF-β → myofibroblast activation → storiform fibrosis</li>
                  <li><b>IL-4, IL-10, IL-13</b>: Drive IgG4 isotype switching and Th2 polarization</li>
                </ul>
                <h4>Drug Targets Modelled</h4>
                <ul>
                  <li><b>Rituximab (RTX)</b>: Anti-CD20 → B cell/plasmablast depletion → ↓IgG4</li>
                  <li><b>Prednisone</b>: Broad GC immunosuppression → ↓Tfh2, ↓cytokines, ↓fibroblasts</li>
                  <li><b>Dupilumab</b> (investigational): Anti-IL-4Rα → block IL-4/IL-13 → ↓IgG4 switch</li>
                </ul>
                <h4>Diagnostic Criteria</h4>
                <ul>
                  <li>Serum IgG4 > 135 mg/dL (sensitivity 90%, specificity 60%)</li>
                  <li>Histopathology: IgG4+ PC ratio > 40%, storiform fibrosis, obliterative phlebitis</li>
                  <li>Multi-organ involvement patterns</li>
                </ul>
              ')
          )
        ),
        fluidRow(
          valueBox(value = "23", subtitle = "ODE Compartments", color = "purple", icon = icon("cogs")),
          valueBox(value = "6",  subtitle = "Treatment Scenarios", color = "blue",   icon = icon("pills")),
          valueBox(value = "91%", subtitle = "RTX Response Rate (Khosroshahi 2012)", color = "green", icon = icon("chart-line")),
          valueBox(value = "135", subtitle = "IgG4 ULN (mg/dL)", color = "yellow", icon = icon("vial"))
        )
      ),

      ## ── TAB 2: Patient Profile ─────────────────────────────────────────
      tabItem("patient",
        fluidRow(
          box(width = 4, title = "Patient Parameters", status = "primary",
              sliderInput("age",        "Age (years)", 20, 80, 55, step = 1),
              sliderInput("igG4_base",  "Baseline Serum IgG4 (mg/dL)", 100, 2000, 450, step = 10),
              sliderInput("severity",   "Disease Severity (1=moderate, 2=severe)", 0.5, 2.0, 1.0, step = 0.1),
              selectInput("organ",      "Primary Organ Involved",
                          choices = c("Pancreas (AIP type 1)", "Kidney (TIN)",
                                      "Salivary Glands (Kuttner)", "Orbit",
                                      "Retroperitoneum", "Biliary (IRC)",
                                      "Lung (ILD)", "Multiple organs")),
              radioButtons("scenario", "Treatment",
                           choices = c("Untreated" = "Untreated",
                                       "Prednisone 40mg taper" = "Prednisone",
                                       "RTX 1g D1+D15" = "RTX_1g2",
                                       "RTX 375mg/m² ×4" = "RTX_375x4",
                                       "RTX + Maintenance" = "RTX_maint",
                                       "Dupilumab 300mg q2w" = "Dupilumab"),
                           selected = "RTX_1g2"),
              sliderInput("sim_months", "Simulation Duration (months)", 6, 36, 24, step = 3),
              actionButton("run_sim", "Run Simulation", class = "btn-primary btn-block")
          ),
          box(width = 8, title = "Patient Summary", status = "info",
              fluidRow(
                column(6, plotlyOutput("pt_igG4_plot", height = "250px")),
                column(6, plotlyOutput("pt_iri_plot",  height = "250px"))
              ),
              DT::dataTableOutput("pt_table")
          )
        )
      ),

      ## ── TAB 3: Pharmacokinetics ─────────────────────────────────────────
      tabItem("pk",
        fluidRow(
          box(width = 4, title = "PK Parameters", status = "warning",
              h4("Rituximab"),
              sliderInput("rtx_dose_mg", "Dose (mg)", 500, 2000, 1000, step = 100),
              numericInput("rtx_n_doses", "Number of doses", 2, min=1, max=8),
              hr(),
              h4("Prednisone"),
              sliderInput("pred_dose_mg", "Starting dose (mg/d)", 10, 60, 40, step = 5),
              sliderInput("pred_dur_wk", "Duration (weeks)", 4, 52, 24, step = 4),
              hr(),
              h4("Dupilumab (SC)"),
              sliderInput("dup_dose_mg", "Dose (mg)", 150, 600, 300, step = 50),
              actionButton("run_pk", "Compute PK", class = "btn-warning btn-block")
          ),
          box(width = 8, title = "Drug Concentration–Time Profiles", status = "warning",
              tabBox(width = 12,
                tabPanel("Rituximab", plotlyOutput("pk_rtx_plot", height = "300px"),
                         plotlyOutput("pk_cd20_plot", height = "200px")),
                tabPanel("Prednisone", plotlyOutput("pk_pred_plot", height = "300px")),
                tabPanel("Dupilumab", plotlyOutput("pk_dup_plot", height = "300px"),
                         plotlyOutput("pk_il4ra_plot", height = "200px"))
              )
          )
        )
      ),

      ## ── TAB 4: B Cell & Immunity ────────────────────────────────────────
      tabItem("bcell",
        fluidRow(
          box(width = 12, title = "B Cell & T Cell Dynamics", status = "danger",
              fluidRow(
                column(6, plotlyOutput("bcell_plot", height = "300px")),
                column(6, plotlyOutput("tcell_plot", height = "300px"))
              ),
              fluidRow(
                column(6, plotlyOutput("plasmablast_plot", height = "300px")),
                column(6, plotlyOutput("igG4_immuno_plot", height = "300px"))
              )
          )
        )
      ),

      ## ── TAB 5: Cytokines & Fibrosis ─────────────────────────────────────
      tabItem("cytokine",
        fluidRow(
          box(width = 12, title = "Cytokine Network & Fibrosis Markers", status = "success",
              fluidRow(
                column(4, plotlyOutput("il4_plot", height = "250px")),
                column(4, plotlyOutput("tgfb_plot", height = "250px")),
                column(4, plotlyOutput("il10_plot", height = "250px"))
              ),
              fluidRow(
                column(6, plotlyOutput("myofib_plot", height = "280px")),
                column(6, plotlyOutput("ecm_plot", height = "280px"))
              )
          )
        )
      ),

      ## ── TAB 6: Scenario Comparison ──────────────────────────────────────
      tabItem("scenario",
        fluidRow(
          box(width = 3, title = "Scenario Settings", status = "info",
              checkboxGroupInput("scen_select", "Select Scenarios:",
                                 choices = c("Untreated" = "Untreated",
                                             "Prednisone" = "Prednisone",
                                             "RTX 1g D1+D15" = "RTX_1g2",
                                             "RTX 375mg/m² ×4" = "RTX_375x4",
                                             "RTX + Maintenance" = "RTX_maint",
                                             "Dupilumab" = "Dupilumab"),
                                 selected = c("Untreated", "Prednisone", "RTX_1g2", "RTX_maint")),
              sliderInput("scen_months", "Duration (months)", 6, 36, 24),
              actionButton("run_scenarios", "Compare Scenarios", class = "btn-info btn-block")
          ),
          box(width = 9, title = "Multi-Scenario Comparison",
              fluidRow(
                column(6, plotlyOutput("comp_igG4", height = "280px")),
                column(6, plotlyOutput("comp_iri", height = "280px"))
              ),
              fluidRow(
                column(6, plotlyOutput("comp_bcell", height = "280px")),
                column(6, plotlyOutput("comp_ecm", height = "280px"))
              ),
              DT::dataTableOutput("scen_table")
          )
        )
      ),

      ## ── TAB 7: Biomarkers ────────────────────────────────────────────────
      tabItem("biomarker",
        fluidRow(
          box(width = 12, title = "Biomarker Dashboard & Clinical Endpoints", status = "primary",
              fluidRow(
                valueBoxOutput("vb_igG4",     width = 3),
                valueBoxOutput("vb_bcell",    width = 3),
                valueBoxOutput("vb_iri",      width = 3),
                valueBoxOutput("vb_response", width = 3)
              ),
              fluidRow(
                column(6,
                  h4("IgG4 Responder Index Over Time"),
                  plotlyOutput("bm_iri_plot", height = "300px")
                ),
                column(6,
                  h4("Serum IgG4 (mg/dL)"),
                  plotlyOutput("bm_igG4_plot", height = "300px")
                )
              ),
              fluidRow(
                column(6,
                  h4("Relapse Probability (IgG4-RD Recurrence)"),
                  plotlyOutput("bm_relapse_plot", height = "280px")
                ),
                column(6,
                  h4("Treatment Response Summary"),
                  DT::dataTableOutput("bm_summary_table")
                )
              )
          )
        )
      )
    ) # end tabItems
  ) # end dashboardBody
)

## ── Server ──────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## Reactive: run single simulation
  sim_data <- eventReactive(list(input$run_sim, input$run_pk), {
    scen <- isolate(input$scenario)
    run_sim(
      scenario      = scen,
      rtx_dose_mg   = isolate(input$rtx_dose_mg),
      rtx_n         = isolate(input$rtx_n_doses),
      pred_dose_mg  = isolate(input$pred_dose_mg),
      pred_dur_wk   = isolate(input$pred_dur_wk),
      dup_dose_mg   = isolate(input$dup_dose_mg),
      sim_months    = isolate(input$sim_months),
      baseline_IgG4 = isolate(input$igG4_base),
      disease_severity = isolate(input$severity)
    )
  }, ignoreNULL = FALSE)

  ## Reactive: run scenario comparison
  scen_data <- eventReactive(input$run_scenarios, {
    scenarios <- input$scen_select
    if (length(scenarios) == 0) scenarios <- "Untreated"
    bind_rows(lapply(scenarios, function(s) {
      run_sim(scenario = s, sim_months = input$scen_months,
              baseline_IgG4 = 450, disease_severity = 1.0) %>%
        mutate(Scenario = s)
    }))
  }, ignoreNULL = FALSE)

  ## ── Patient Profile Plots ─────────────────────────────────────────────
  output$pt_igG4_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time_months, y = IgG4_SER)) +
      geom_line(color = "#9C27B0", linewidth = 1.2) +
      geom_hline(yintercept = 135, linetype = "dashed", color = "red") +
      labs(x = "Months", y = "Serum IgG4 (mg/dL)", title = "IgG4 Response") +
      theme_bw()
    ggplotly(p)
  })

  output$pt_iri_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time_months, y = IRI)) +
      geom_line(color = "#F44336", linewidth = 1.2) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_hline(yintercept = 3, linetype = "dashed", color = "orange") +
      labs(x = "Months", y = "IRI (0-24)", title = "Disease Activity (IRI)") +
      ylim(0, 12) + theme_bw()
    ggplotly(p)
  })

  output$pt_table <- DT::renderDataTable({
    d <- sim_data()
    summary_pts <- d %>%
      filter(time_months %in% c(0, 3, 6, 12, 18, 24) |
               abs(time_months - round(time_months)) < 0.1) %>%
      distinct(time_months, .keep_all = TRUE) %>%
      filter(time_months %in% c(0, 3, 6, 12, 18, 24)) %>%
      select(Month = time_months, `IgG4 (mg/dL)` = IgG4_SER,
             `IRI Score` = IRI, `B Cells (%)` = Bcell_pct,
             `TFH2 (rel)` = TFH2, `ECM Index` = ECM) %>%
      mutate(across(where(is.numeric), ~ round(.x, 1)))
    DT::datatable(summary_pts, options = list(pageLength = 8, dom = 't'),
                  rownames = FALSE)
  })

  ## ── PK Plots ─────────────────────────────────────────────────────────
  output$pk_rtx_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time_months, y = RTX_conc)) +
      geom_line(color = "#E91E63", linewidth = 1.2) +
      labs(x = "Months", y = "Rituximab (μg/mL)", title = "Rituximab Plasma Concentration") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_cd20_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time_months, y = CD20_occ_pct)) +
      geom_line(color = "#880E4F", linewidth = 1.0) +
      labs(x = "Months", y = "CD20 Occupancy (%)", title = "CD20 Target Occupancy") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_pred_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time_months, y = CENT_PRED)) +
      geom_line(color = "#FF8F00", linewidth = 1.2) +
      labs(x = "Months", y = "Prednisolone (nM)", title = "Prednisolone Plasma Level") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_dup_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time_months, y = CENT_DUP * 146000 / 1e6 * 1000)) +
      geom_line(color = "#2E7D32", linewidth = 1.2) +
      labs(x = "Months", y = "Dupilumab (μg/mL)", title = "Dupilumab Plasma Concentration") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_il4ra_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time_months, y = IL4RA_occ_pct)) +
      geom_line(color = "#1B5E20", linewidth = 1.0) +
      labs(x = "Months", y = "IL-4Rα Occupancy (%)", title = "IL-4Rα Target Engagement") +
      theme_bw()
    ggplotly(p)
  })

  ## ── B Cell / T Cell Plots ─────────────────────────────────────────────
  output$bcell_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time_months, y = Bcell_pct)) +
      geom_line(color = "#1565C0", linewidth = 1.2) +
      geom_hline(yintercept = 5, linetype = "dashed", color = "red") +
      labs(x = "Months", y = "B Cells (% baseline)", title = "Naïve B Cell Level") +
      theme_bw()
    ggplotly(p)
  })

  output$tcell_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d) +
      geom_line(aes(x = time_months, y = TFH2, color = "Tfh2 (rel)"), linewidth = 1.2) +
      geom_line(aes(x = time_months, y = CTL4, color = "CD4 CTL (rel)"), linewidth = 1.2) +
      scale_color_manual(values = c("Tfh2 (rel)" = "#FFC107", "CD4 CTL (rel)" = "#F44336")) +
      labs(x = "Months", y = "Relative Units", title = "Pathogenic T Cell Subsets",
           color = NULL) +
      theme_bw() + theme(legend.position = "top")
    ggplotly(p)
  })

  output$plasmablast_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d) +
      geom_line(aes(x = time_months, y = PB, color = "Plasmablasts"), linewidth = 1.2) +
      geom_line(aes(x = time_months, y = PC, color = "Long-lived PC"), linewidth = 1.2) +
      scale_color_manual(values = c("Plasmablasts" = "#FF6F00",
                                     "Long-lived PC" = "#BF360C")) +
      labs(x = "Months", y = "Relative Units", title = "Plasma Cell Compartments",
           color = NULL) +
      theme_bw() + theme(legend.position = "top")
    ggplotly(p)
  })

  output$igG4_immuno_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time_months, y = IgG4_SER)) +
      geom_line(color = "#9C27B0", linewidth = 1.2) +
      geom_hline(yintercept = 135, linetype = "dashed", color = "red") +
      labs(x = "Months", y = "IgG4 (mg/dL)", title = "Serum IgG4 Dynamics") +
      theme_bw()
    ggplotly(p)
  })

  ## ── Cytokine / Fibrosis Plots ─────────────────────────────────────────
  output$il4_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time_months, y = IL4)) +
      geom_line(color = "#F9A825", linewidth = 1.2) +
      labs(x = "Months", y = "IL-4 (rel)", title = "IL-4 Dynamics") + theme_bw()
    ggplotly(p)
  })

  output$tgfb_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time_months, y = TGFB)) +
      geom_line(color = "#6A1B9A", linewidth = 1.2) +
      labs(x = "Months", y = "TGF-β1 (rel)", title = "TGF-β1 Dynamics") + theme_bw()
    ggplotly(p)
  })

  output$il10_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time_months, y = IL10)) +
      geom_line(color = "#388E3C", linewidth = 1.2) +
      labs(x = "Months", y = "IL-10 (rel)", title = "IL-10 Dynamics") + theme_bw()
    ggplotly(p)
  })

  output$myofib_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time_months, y = MYOFIB)) +
      geom_line(color = "#CE93D8", linewidth = 1.2) +
      geom_hline(yintercept = 1.0, linetype = "dashed") +
      labs(x = "Months", y = "Myofibroblast Activation (rel)", title = "Myofibroblast Activation") +
      theme_bw()
    ggplotly(p)
  })

  output$ecm_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time_months, y = ECM)) +
      geom_line(color = "#7B1FA2", linewidth = 1.2) +
      geom_hline(yintercept = 1.0, linetype = "dashed") +
      labs(x = "Months", y = "ECM Fibrosis Index (rel)", title = "Fibrosis Index") +
      theme_bw()
    ggplotly(p)
  })

  ## ── Scenario Comparison Plots ─────────────────────────────────────────
  output$comp_igG4 <- renderPlotly({
    d <- scen_data()
    p <- ggplot(d, aes(x = time_months, y = IgG4_SER, color = Scenario)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = 135, linetype = "dashed") +
      labs(x = "Months", y = "IgG4 (mg/dL)", title = "Serum IgG4 Comparison") +
      theme_bw() + theme(legend.position = "top", legend.text = element_text(size = 8))
    ggplotly(p)
  })

  output$comp_iri <- renderPlotly({
    d <- scen_data()
    p <- ggplot(d, aes(x = time_months, y = IRI, color = Scenario)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = 3, linetype = "dashed", color = "orange") +
      labs(x = "Months", y = "IRI (0-24)", title = "Disease Activity (IRI) Comparison") +
      ylim(0, 10) + theme_bw() + theme(legend.position = "top", legend.text = element_text(size = 8))
    ggplotly(p)
  })

  output$comp_bcell <- renderPlotly({
    d <- scen_data()
    p <- ggplot(d, aes(x = time_months, y = Bcell_pct, color = Scenario)) +
      geom_line(linewidth = 1.1) +
      labs(x = "Months", y = "B Cells (% baseline)", title = "B Cell Depletion Comparison") +
      theme_bw() + theme(legend.position = "top", legend.text = element_text(size = 8))
    ggplotly(p)
  })

  output$comp_ecm <- renderPlotly({
    d <- scen_data()
    p <- ggplot(d, aes(x = time_months, y = ECM, color = Scenario)) +
      geom_line(linewidth = 1.1) +
      labs(x = "Months", y = "ECM (rel)", title = "Fibrosis Index Comparison") +
      theme_bw() + theme(legend.position = "top", legend.text = element_text(size = 8))
    ggplotly(p)
  })

  output$scen_table <- DT::renderDataTable({
    d <- scen_data()
    tbl <- d %>%
      group_by(Scenario) %>%
      summarise(
        `IgG4 Baseline (mg/dL)` = round(first(IgG4_SER), 0),
        `IgG4 Wk12 (mg/dL)` = round(IgG4_SER[which.min(abs(time_months - 3))], 0),
        `IgG4 Wk24 (mg/dL)` = round(IgG4_SER[which.min(abs(time_months - 6))], 0),
        `IRI Baseline` = round(first(IRI), 1),
        `IRI Wk24` = round(IRI[which.min(abs(time_months - 6))], 1),
        `B-Cell Nadir (%)` = round(min(Bcell_pct), 1),
        `CR at 6mo (%)` = round(mean(CR_flag[abs(time_months - 6) < 0.3]) * 100, 0),
        .groups = "drop"
      )
    DT::datatable(tbl, options = list(dom = 't', pageLength = 10), rownames = FALSE)
  })

  ## ── Biomarker Value Boxes ──────────────────────────────────────────────
  output$vb_igG4 <- renderValueBox({
    d <- sim_data()
    last_igG4 <- tail(d$IgG4_SER, 1)
    valueBox(sprintf("%.0f mg/dL", last_igG4), "Serum IgG4 (final)",
             color = ifelse(last_igG4 > 135, "red", "green"), icon = icon("vial"))
  })

  output$vb_bcell <- renderValueBox({
    d <- sim_data()
    nadir <- min(d$Bcell_pct)
    valueBox(sprintf("%.1f%%", nadir), "B-Cell Nadir",
             color = ifelse(nadir < 5, "green", "yellow"), icon = icon("microscope"))
  })

  output$vb_iri <- renderValueBox({
    d <- sim_data()
    last_iri <- tail(d$IRI, 1)
    valueBox(sprintf("%.1f", last_iri), "IRI at End",
             color = ifelse(last_iri < 1, "green",
                            ifelse(last_iri < 4, "yellow", "red")),
             icon = icon("chart-line"))
  })

  output$vb_response <- renderValueBox({
    d <- sim_data()
    cr <- any(d$CR_flag > 0.5)
    pr <- any(d$PR_flag > 0.5)
    label <- if (cr) "Complete Response" else if (pr) "Partial Response" else "No Response"
    valueBox(label, "Best Response",
             color = if (cr) "green" else if (pr) "yellow" else "red",
             icon = icon("check-circle"))
  })

  output$bm_iri_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time_months, y = IRI)) +
      geom_line(color = "#F44336", linewidth = 1.3) +
      geom_hline(yintercept = c(0, 3), linetype = c("dotted", "dashed"),
                 color = c("black", "orange")) +
      annotate("text", x = 1, y = 3.3, label = "Partial Response", size = 3, color = "orange") +
      labs(x = "Months", y = "IgG4-RD Responder Index (IRI, 0-24)",
           title = "Disease Activity Trajectory") +
      ylim(0, 12) + theme_bw()
    ggplotly(p)
  })

  output$bm_igG4_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = time_months, y = IgG4_SER)) +
      geom_line(color = "#9C27B0", linewidth = 1.3) +
      geom_hline(yintercept = 135, linetype = "dashed", color = "red") +
      annotate("text", x = 0.5, y = 140, label = "ULN 135 mg/dL", size = 3, color = "red") +
      labs(x = "Months", y = "Serum IgG4 (mg/dL)", title = "IgG4 Biomarker Trajectory") +
      theme_bw()
    ggplotly(p)
  })

  output$bm_relapse_plot <- renderPlotly({
    # Relapse risk = function of residual PC and MBC (simplified)
    d <- sim_data()
    d <- d %>% mutate(
      relapse_risk = pmin(100, pmax(0, (PC * 15 + TFH2 * 10 + IgG4_SER / 45) / 3))
    )
    p <- ggplot(d, aes(x = time_months, y = relapse_risk)) +
      geom_area(fill = "#FFCDD2", alpha = 0.5) +
      geom_line(color = "#C62828", linewidth = 1.2) +
      geom_hline(yintercept = 40, linetype = "dashed", color = "orange") +
      labs(x = "Months", y = "Relapse Risk Score (0-100)",
           title = "Relapse Risk Projection") +
      theme_bw()
    ggplotly(p)
  })

  output$bm_summary_table <- DT::renderDataTable({
    d <- sim_data()
    timepoints <- c(0, 3, 6, 9, 12, 18, 24)
    tbl <- lapply(timepoints, function(mo) {
      row <- d[which.min(abs(d$time_months - mo)), ]
      data.frame(
        Month = mo,
        `IgG4 (mg/dL)` = round(row$IgG4_SER, 0),
        `IRI Score` = round(row$IRI, 1),
        `B-Cells (%)` = round(row$Bcell_pct, 1),
        `TGF-β (rel)` = round(row$TGFB, 2),
        `ECM Index` = round(row$ECM, 2),
        Response = ifelse(row$IRI < 1, "CR",
                          ifelse(row$IRI < 4, "PR", "NR")),
        check.names = FALSE
      )
    })
    tbl_df <- bind_rows(tbl)
    DT::datatable(tbl_df, options = list(dom = 't', pageLength = 10), rownames = FALSE) %>%
      DT::formatStyle("Response",
                       color = DT::styleEqual(c("CR","PR","NR"),
                                               c("green","orange","red")))
  })
}

## ── Launch ─────────────────────────────────────────────────────────────────
shinyApp(ui, server)
