## ============================================================
## Evans Syndrome (ES) – Shiny Dashboard
## Combined AIHA + ITP QSP Interactive Simulator
## ============================================================
## Author : Claude Code (CCR)  |  Date: 2026-06-20
## Tabs:
##   1. Patient Profile & Disease Parameters
##   2. Drug PK – Plasma/Blood Concentrations
##   3. Hematologic Response (Hgb, Plt, Reticulocyte)
##   4. Immunologic Markers (B cells, Ab, Treg, Complement)
##   5. Scenario Comparison (6 treatment arms)
##   6. Biomarker Dashboard (DAT, PAIgG, Treg%, B cell%)
## ============================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(DT)

## ============================================================
## mrgsolve model (embedded)
## ============================================================
ES_code <- '
$PARAM
B_RBC0=500, B_Plt0=500, Ab_RBC0=2.0, Ab_Plt0=2.0,
Treg0=80, Treg_min=20, RBC0=5.0, Hgb0=14.0,
Hgb_per_RBC=2.8, Retic0=1.5, BM_Ery0=100,
Plt0=250, MK0=100,
kprol_B=0.025, kdeg_B=0.020,
kB_Ab_RBC=0.004, kB_Ab_Plt=0.004, kclear_Ab=0.10,
kprol_Treg=0.015, kdeg_Treg=0.012,
ISuppression=0.4, EC50_Treg=60,
kform_C3b=0.05, kclear_C3b=0.30, C3b0=0.1,
kops_RBC=0.015, kphago_RBC=0.30,
kdeg_RBC=0.012, kinput_RBC=0.060,
EPO_EC50=7.0, kEPO_stim=3.0,
kprol_BM_Ery=0.050, kdeg_BM_Ery=0.040, kmat_Retic=0.50,
kops_Plt=0.018, kphago_Plt=0.35, kdeg_Plt=0.043,
kinput_Plt=10.8, kprol_MK=0.040, kdeg_MK=0.030, krel_Plt=0.108,
PRED_dose=0, MMF_dose=0, ELT_dose=0, SIRO_dose=0, SPLEN=0,
CL_Pred=15, Vc_Pred=45, F_Pred=0.82, ka_Pred=8,
CL_MPA=25, Vc_MPA=15, F_MPA=0.94, ka_MPA=6,
CL_Elt=12, Vc_Elt=20, F_Elt=0.52, ka_Elt=3,
CL_Siro=15, Vc_Siro=480, F_Siro=0.14, ka_Siro=2,
CL_IVIg=0.30, Vc_IVIg=3.5, CL_Rtx=0.33, Vc_Rtx=4, Km_Rtx=5,
Emax_CS_B=0.70, EC50_CS_B=0.08, Emax_CS_Mac=0.60, EC50_CS_Mac=0.06,
Emax_IVIg_FcR=0.75, EC50_IVIg=10, kcat_IgG=0.05,
Emax_Rtx_B=0.95, EC50_Rtx=1.0,
Emax_MPA_Lym=0.65, EC50_MPA=1.5,
Emax_Elt_MK=2.5, EC50_Elt=0.5,
Emax_Siro_Treg=2.0, EC50_Siro_Treg=1.5,
Splen_Phago_red=0.80,
RTX_on=0

$CMT
DEPOT_PRED CENT_PRED DEPOT_MPA CENT_MPA
DEPOT_ELT CENT_ELT DEPOT_SIRO CENT_SIRO
CENT_IVIG CENT_RTX
B_RBC AB_RBC B_PLT AB_PLT TREG C3B
ORBC RBC_CIRC RETIC BM_ERY
OPLT PLT_CIRC MK

$MAIN
DEPOT_PRED_0=0; CENT_PRED_0=0;
DEPOT_MPA_0=0; CENT_MPA_0=0;
DEPOT_ELT_0=0; CENT_ELT_0=0;
DEPOT_SIRO_0=0; CENT_SIRO_0=0;
CENT_IVIG_0=0; CENT_RTX_0=0;
B_RBC_0=B_RBC0; AB_RBC_0=Ab_RBC0;
B_PLT_0=B_Plt0; AB_PLT_0=Ab_Plt0;
TREG_0=Treg0; C3B_0=C3b0;
ORBC_0=0.1; RBC_CIRC_0=RBC0;
RETIC_0=Retic0; BM_ERY_0=BM_Ery0;
OPLT_0=5.0; PLT_CIRC_0=Plt0; MK_0=MK0;

$ODE
double dose_PRED_daily = PRED_dose * 70.0;
double C_Pred = CENT_PRED / Vc_Pred;
dxdt_DEPOT_PRED = F_Pred * dose_PRED_daily * ka_Pred - ka_Pred * DEPOT_PRED;
dxdt_CENT_PRED  = ka_Pred * DEPOT_PRED - CL_Pred * C_Pred;

double dose_MPA_daily = MMF_dose * 1000.0;
double C_MPA = CENT_MPA / Vc_MPA;
dxdt_DEPOT_MPA = F_MPA * dose_MPA_daily * ka_MPA - ka_MPA * DEPOT_MPA;
dxdt_CENT_MPA  = ka_MPA * DEPOT_MPA - CL_MPA * C_MPA;

double C_Elt = CENT_ELT / Vc_Elt;
dxdt_DEPOT_ELT = F_Elt * ELT_dose * ka_Elt - ka_Elt * DEPOT_ELT;
dxdt_CENT_ELT  = ka_Elt * DEPOT_ELT - CL_Elt * C_Elt;

double C_Siro = CENT_SIRO / Vc_Siro;
dxdt_DEPOT_SIRO = F_Siro * SIRO_dose * ka_Siro - ka_Siro * DEPOT_SIRO;
dxdt_CENT_SIRO  = ka_Siro * DEPOT_SIRO - CL_Siro * C_Siro;

double C_IVIg = CENT_IVIG / Vc_IVIg;
double IVIg_cat_extra = kcat_IgG * C_IVIg / (EC50_IVIg + C_IVIg);
dxdt_CENT_IVIG = -CL_IVIg * C_IVIg;

double C_Rtx = CENT_RTX / Vc_Rtx;
dxdt_CENT_RTX = -(CL_Rtx + CL_Rtx * Km_Rtx / (Km_Rtx + C_Rtx)) * C_Rtx;

double E_CS_B    = Emax_CS_B    * C_Pred / (EC50_CS_B    + C_Pred);
double E_CS_Mac  = Emax_CS_Mac  * C_Pred / (EC50_CS_Mac  + C_Pred);
double E_IVIg_FcR = Emax_IVIg_FcR * C_IVIg / (EC50_IVIg + C_IVIg);
double E_Rtx_B   = Emax_Rtx_B   * C_Rtx  / (EC50_Rtx    + C_Rtx);
double E_MPA_Lym = Emax_MPA_Lym * C_MPA  / (EC50_MPA    + C_MPA);
double E_Elt_MK  = 1.0 + (Emax_Elt_MK - 1.0) * C_Elt / (EC50_Elt + C_Elt);
double E_Siro_T  = 1.0 + (Emax_Siro_Treg - 1.0) * C_Siro / (EC50_Siro_Treg + C_Siro);
double Splen_f   = 1.0 - SPLEN * Splen_Phago_red;
double FcgR_f    = Splen_f * (1.0 - E_CS_Mac) * (1.0 - E_IVIg_FcR);
double Treg_sup  = ISuppression * TREG / (EC50_Treg + TREG);
double B_net_sup = (1.0 - Treg_sup) * (1.0 - E_Rtx_B) * (1.0 - E_MPA_Lym);

dxdt_TREG   = kprol_Treg * E_Siro_T * Treg0 - kdeg_Treg * TREG;
dxdt_B_RBC  = kprol_B * B_net_sup * B_RBC0 - kdeg_B * B_RBC * (1.0 + E_CS_B);
dxdt_B_PLT  = kprol_B * B_net_sup * B_Plt0 - kdeg_B * B_PLT * (1.0 + E_CS_B);
dxdt_AB_RBC = kB_Ab_RBC * B_RBC * (1.0 - E_CS_B) - (kclear_Ab + IVIg_cat_extra) * AB_RBC;
dxdt_AB_PLT = kB_Ab_Plt * B_PLT * (1.0 - E_CS_B) - (kclear_Ab + IVIg_cat_extra) * AB_PLT;
dxdt_C3B    = kform_C3b * AB_RBC * RBC_CIRC - kclear_C3b * C3B;

double kops_eff = kops_RBC * (AB_RBC + 0.5 * C3B);
dxdt_ORBC     = kops_eff * RBC_CIRC - kphago_RBC * FcgR_f * ORBC - kdeg_RBC * ORBC;

double Hgb_calc = RBC_CIRC * Hgb_per_RBC;
double EPO_f    = 1.0 + (kEPO_stim - 1.0) * EPO_EC50 / (EPO_EC50 + Hgb_calc);
dxdt_BM_ERY = kprol_BM_Ery * EPO_f * BM_Ery0 - kdeg_BM_Ery * BM_ERY;
dxdt_RETIC  = kdeg_BM_Ery * BM_ERY * 0.10 - kmat_Retic * RETIC;
dxdt_RBC_CIRC = kmat_Retic * RETIC * RBC0 / Retic0 - kops_eff * RBC_CIRC - kdeg_RBC * RBC_CIRC;

double kops_Plt_eff = kops_Plt * AB_PLT;
dxdt_OPLT     = kops_Plt_eff * PLT_CIRC - kphago_Plt * FcgR_f * OPLT - kdeg_Plt * OPLT;
double TPO_f  = 1.0 + 2.0 * Plt0 / (Plt0 + PLT_CIRC);
dxdt_MK       = kprol_MK * TPO_f * E_Elt_MK * MK0 - kdeg_MK * MK;
dxdt_PLT_CIRC = krel_Plt * MK - kops_Plt_eff * PLT_CIRC - kdeg_Plt * PLT_CIRC;

$TABLE
capture Hgb           = RBC_CIRC * Hgb_per_RBC;
capture Plt           = PLT_CIRC;
capture Retic_pct     = RETIC;
capture C_Pred_pl     = CENT_PRED / Vc_Pred;
capture C_MPA_pl      = CENT_MPA / Vc_MPA;
capture C_Rtx_pl      = CENT_RTX / Vc_Rtx;
capture C_Elt_pl      = CENT_ELT / Vc_Elt;
capture C_Siro_wh     = CENT_SIRO / Vc_Siro;
capture C_IVIg_pl     = CENT_IVIG / Vc_IVIg;
capture Treg_count    = TREG;
capture Ab_RBC_lvl    = AB_RBC;
capture Ab_Plt_lvl    = AB_PLT;
capture B_RBC_cnt     = B_RBC;
capture B_Plt_cnt     = B_PLT;
capture C3b_lvl       = C3B;
capture MK_pool       = MK;
'

es_mod <- mcode("ES_Shiny", ES_code, quiet=TRUE)

## ============================================================
## Simulation helper
## ============================================================
run_sim <- function(params, rtx_events=NULL, ivig_events=NULL,
                    sim_days=365) {
  mod_p <- es_mod %>% param(params)
  ev_list <- ev(time=0, amt=0, cmt=1)
  if (!is.null(rtx_events))  ev_list <- ev_list + rtx_events
  if (!is.null(ivig_events)) ev_list <- ev_list + ivig_events
  mod_p %>%
    ev(ev_list) %>%
    mrgsim(end=sim_days, delta=1, obsonly=TRUE) %>%
    as.data.frame()
}

## ============================================================
## UI
## ============================================================
ui <- navbarPage(
  title = "Evans Syndrome QSP Dashboard",
  theme = bslib::bs_theme(bootswatch="flatly"),

  ## ── Tab 1: Patient Profile ────────────────────────────────
  tabPanel("Patient Profile",
    sidebarLayout(
      sidebarPanel(width=4,
        h4("Patient Demographics"),
        numericInput("weight", "Body Weight (kg)", 70, 30, 120, 5),
        numericInput("age",    "Age (years)",       35,  5,  85, 1),
        selectInput("sex",     "Sex", choices=c("Female","Male")),
        tags$hr(),
        h4("Baseline Hematology"),
        numericInput("hgb_bl",  "Hemoglobin (g/dL)",       7.5, 3, 16, 0.5),
        numericInput("plt_bl",  "Platelet (×10⁹/L)",        25,  1, 400, 5),
        numericInput("retic_bl","Reticulocyte (%)",           6,  0,  20, 0.5),
        tags$hr(),
        h4("Immunology"),
        selectInput("dat_status",  "DAT Result",
                    c("IgG only","IgG + C3d","C3d only","Negative")),
        numericInput("b_cell_pct", "B cells (% lymphocytes)",  18, 1, 40, 1),
        numericInput("treg_pct",   "Treg (% CD4+ T cells)",     2, 0, 15, 0.5),
        tags$hr(),
        h4("Disease Severity"),
        selectInput("severity", "Current Severity",
                    c("Mild (Hgb 10-12, Plt 50-100)",
                      "Moderate (Hgb 8-10, Plt 20-50)",
                      "Severe (Hgb <8, Plt <20)")),
        selectInput("course", "Disease Course",
                    c("Chronic relapsing","Acute","Single episode")),
        h4("Prior Treatments"),
        checkboxGroupInput("prior_tx", NULL,
                           c("Steroids","IVIg","Rituximab","MMF",
                             "Danazol","Splenectomy"))
      ),
      mainPanel(width=8,
        fluidRow(
          column(6,
            h4("Disease Summary"),
            tableOutput("profile_table")
          ),
          column(6,
            h4("Pathophysiology Overview"),
            img(src=NULL, alt="Mechanistic map",
                style="max-width:100%; border:1px solid #ccc;
                       padding:4px; border-radius:4px;"),
            tags$p("Evans Syndrome = concurrent or sequential AIHA + ITP
             mediated by polyclonal autoantibodies targeting RBC and
             platelet surface antigens.", style="font-size:12px;color:#555;")
          )
        ),
        fluidRow(
          column(12,
            h4("Key Biomarker Targets"),
            DT::dataTableOutput("biomarker_targets")
          )
        )
      )
    )
  ),

  ## ── Tab 2: Drug PK ────────────────────────────────────────
  tabPanel("Drug PK",
    sidebarLayout(
      sidebarPanel(width=3,
        h4("Treatment Selection"),
        checkboxGroupInput("pk_drugs", "Drugs to simulate",
                           choices=c("Prednisone","IVIg","Rituximab",
                                     "MMF","Eltrombopag","Sirolimus"),
                           selected="Prednisone"),
        conditionalPanel("input.pk_drugs.includes('Prednisone')",
          numericInput("pred_dose_pk", "Prednisone (mg/kg/day)", 1.0, 0.1, 3, 0.1)
        ),
        conditionalPanel("input.pk_drugs.includes('IVIg')",
          numericInput("ivig_total_pk", "IVIg total dose (g/kg)", 1.0, 0.4, 2, 0.2),
          numericInput("ivig_day_pk",   "Day of infusion", 1, 1, 30, 1)
        ),
        conditionalPanel("input.pk_drugs.includes('Rituximab')",
          numericInput("rtx_bsa_pk", "BSA (m²)", 1.7, 1.2, 2.2, 0.1)
        ),
        conditionalPanel("input.pk_drugs.includes('MMF')",
          numericInput("mmf_dose_pk", "MMF dose (g/day)", 2.0, 0.5, 3.0, 0.5)
        ),
        conditionalPanel("input.pk_drugs.includes('Eltrombopag')",
          numericInput("elt_dose_pk", "Eltrombopag (mg/day)", 50, 25, 75, 25)
        ),
        conditionalPanel("input.pk_drugs.includes('Sirolimus')",
          numericInput("siro_dose_pk", "Sirolimus (mg/day)", 4, 1, 10, 1)
        ),
        numericInput("sim_days_pk", "Simulation days", 90, 14, 365, 14),
        actionButton("run_pk", "Simulate PK", class="btn-primary btn-sm")
      ),
      mainPanel(width=9,
        fluidRow(
          column(6, plotOutput("pk_pred_plot", height="260px")),
          column(6, plotOutput("pk_rtx_plot",  height="260px"))
        ),
        fluidRow(
          column(6, plotOutput("pk_mpa_plot",  height="260px")),
          column(6, plotOutput("pk_elt_plot",  height="260px"))
        ),
        fluidRow(
          column(6, plotOutput("pk_siro_plot", height="260px")),
          column(6, plotOutput("pk_ivig_plot", height="260px"))
        ),
        h5("PK Parameter Summary"),
        DT::dataTableOutput("pk_params_table")
      )
    )
  ),

  ## ── Tab 3: Hematologic Response ───────────────────────────
  tabPanel("Hematologic Response",
    sidebarLayout(
      sidebarPanel(width=3,
        h4("Treatment Arm"),
        numericInput("pred_dose_hem", "Prednisone (mg/kg/day)", 1.0, 0, 3, 0.1),
        checkboxInput("rtx_hem", "Add Rituximab (4 weekly doses)", FALSE),
        numericInput("mmf_dose_hem", "MMF dose (g/day, 0=off)", 0, 0, 3, 0.5),
        numericInput("elt_dose_hem", "Eltrombopag (mg/day, 0=off)", 0, 0, 75, 25),
        numericInput("siro_dose_hem", "Sirolimus (mg/day, 0=off)", 0, 0, 10, 1),
        checkboxInput("splen_hem", "Splenectomy (day 1)", FALSE),
        tags$hr(),
        numericInput("sim_days_hem", "Simulation (days)", 180, 30, 730, 30),
        actionButton("run_hem", "Run Simulation", class="btn-primary btn-sm"),
        tags$hr(),
        h5("Response Criteria"),
        tags$ul(style="font-size:11px;",
          tags$li("CR: Hgb≥12 g/dL AND Plt≥100×10⁹/L"),
          tags$li("PR: Improvement but not CR"),
          tags$li("NR: No response")
        )
      ),
      mainPanel(width=9,
        fluidRow(
          column(6, plotOutput("hem_hgb_plot",  height="280px")),
          column(6, plotOutput("hem_plt_plot",  height="280px"))
        ),
        fluidRow(
          column(6, plotOutput("hem_retic_plot", height="280px")),
          column(6, plotOutput("hem_mk_plot",    height="280px"))
        ),
        h4("Response Milestone Table"),
        DT::dataTableOutput("hem_milestone_table")
      )
    )
  ),

  ## ── Tab 4: Immunologic Markers ────────────────────────────
  tabPanel("Immunologic Markers",
    sidebarLayout(
      sidebarPanel(width=3,
        h4("Immunosuppressive Regimen"),
        numericInput("pred_dose_imm", "Prednisone (mg/kg/day)", 1.0, 0, 3, 0.1),
        checkboxInput("rtx_imm", "Rituximab (4 weekly doses)", FALSE),
        numericInput("mmf_dose_imm", "MMF (g/day)", 0, 0, 3, 0.5),
        numericInput("siro_dose_imm", "Sirolimus (mg/day)", 0, 0, 10, 1),
        numericInput("sim_days_imm", "Simulation (days)", 365, 30, 730, 30),
        actionButton("run_imm", "Run Simulation", class="btn-primary btn-sm"),
        tags$hr(),
        h5("Reference Ranges"),
        tags$ul(style="font-size:11px;",
          tags$li("Anti-RBC IgG: Normal <0.1 mg/L"),
          tags$li("Anti-Plt IgG: Normal <0.1 mg/L"),
          tags$li("Treg: Normal 80–120 cells/µL"),
          tags$li("B cells: Normal 150–600 cells/µL"),
          tags$li("C3b: Normal ~0.1 AU")
        )
      ),
      mainPanel(width=9,
        fluidRow(
          column(6, plotOutput("imm_ab_rbc_plot", height="270px")),
          column(6, plotOutput("imm_ab_plt_plot", height="270px"))
        ),
        fluidRow(
          column(6, plotOutput("imm_treg_plot", height="270px")),
          column(6, plotOutput("imm_bcell_plot", height="270px"))
        ),
        fluidRow(
          column(6, plotOutput("imm_c3b_plot",  height="270px")),
          column(6, plotOutput("imm_bm_ery_plot", height="270px"))
        )
      )
    )
  ),

  ## ── Tab 5: Scenario Comparison ────────────────────────────
  tabPanel("Scenario Comparison",
    sidebarLayout(
      sidebarPanel(width=3,
        h4("Simulation Parameters"),
        numericInput("sim_days_sc", "Simulation (days)", 365, 90, 730, 30),
        checkboxGroupInput("sc_select", "Scenarios to compare",
          choices=c(
            "Untreated"                      = "sc1",
            "Prednisone 1.5mg/kg"            = "sc2",
            "Pred 0.5mg/kg + Rituximab"      = "sc3",
            "Pred 0.25mg/kg + MMF 2g/d"      = "sc4",
            "Sirolimus 4mg + Eltrombopag 50mg" = "sc5",
            "Splenectomy + Pred 0.1mg/kg"    = "sc6"
          ),
          selected=c("sc1","sc2","sc3","sc4")
        ),
        actionButton("run_sc", "Compare Scenarios", class="btn-primary btn-sm"),
        tags$hr(),
        h5("Clinical Endpoints"),
        selectInput("sc_endpoint", "Primary endpoint",
                    c("Hemoglobin (g/dL)","Platelet (×10⁹/L)",
                      "Both (composite CR)"))
      ),
      mainPanel(width=9,
        plotOutput("sc_hgb_plot",  height="300px"),
        plotOutput("sc_plt_plot",  height="300px"),
        h4("Response at Day 180 & Day 365"),
        DT::dataTableOutput("sc_response_table")
      )
    )
  ),

  ## ── Tab 6: Biomarker Dashboard ────────────────────────────
  tabPanel("Biomarker Dashboard",
    sidebarLayout(
      sidebarPanel(width=3,
        h4("Biomarker Inputs"),
        numericInput("bm_hgb",  "Hemoglobin (g/dL)",  7.5, 3, 16, 0.5),
        numericInput("bm_plt",  "Platelet (×10⁹/L)",   25, 1, 400, 5),
        numericInput("bm_retic","Reticulocyte (%)",       8, 0, 25, 0.5),
        selectInput("bm_dat",   "DAT",
                    c("IgG +3","IgG +2","IgG +1","IgG+C3d","Negative")),
        numericInput("bm_paIgG","PAIgG (×10⁻¹⁵g/Plt)",  80, 5, 300, 5),
        numericInput("bm_treg_pct","Treg (% CD4)",     2.0, 0, 15, 0.5),
        numericInput("bm_bcell_pct","B cell (% lymph)", 18, 1, 40, 1),
        numericInput("bm_ldh",  "LDH (U/L)",          420, 120, 2000, 20),
        numericInput("bm_haptoglobin","Haptoglobin (g/L)", 0.1, 0, 2, 0.1),
        numericInput("bm_bili",  "Indirect bilirubin (µmol/L)", 55, 3, 200, 5),
        actionButton("run_bm", "Evaluate", class="btn-primary btn-sm")
      ),
      mainPanel(width=9,
        fluidRow(
          column(6,
            h4("Disease Activity Score"),
            plotOutput("bm_gauge_plot", height="220px")
          ),
          column(6,
            h4("Biomarker Radar"),
            plotOutput("bm_radar_plot", height="220px")
          )
        ),
        fluidRow(
          column(12,
            h4("Biomarker Reference Summary"),
            DT::dataTableOutput("bm_summary_table")
          )
        ),
        fluidRow(
          column(12,
            h4("Treatment Recommendation Logic"),
            uiOutput("bm_treatment_recommend")
          )
        )
      )
    )
  )
)

## ============================================================
## Server
## ============================================================
server <- function(input, output, session) {

  ## ----- Patient Profile -----
  output$profile_table <- renderTable({
    data.frame(
      Parameter = c("Weight","Age","Sex","Hgb","Platelet","Reticulocyte",
                    "DAT Status","Severity","Disease Course"),
      Value = c(
        paste0(input$weight," kg"), paste0(input$age," y"), input$sex,
        paste0(input$hgb_bl," g/dL"),
        paste0(input$plt_bl," ×10⁹/L"),
        paste0(input$retic_bl," %"),
        input$dat_status, input$severity, input$course
      )
    )
  })

  output$biomarker_targets <- DT::renderDataTable({
    data.frame(
      Biomarker = c("Hgb","Plt","Retic%","DAT","PAIgG","LDH","Haptoglobin",
                    "Indirect Bili","Treg %","B cell %"),
      Normal_Range = c(">12 g/dL",">150×10⁹/L","0.5–2%","Negative",
                       "<20×10⁻¹⁵g","<250 U/L","0.5–2.5 g/L","<20 µmol/L",
                       "5–10% CD4","5–15% lymph"),
      Current = c(
        paste0(input$hgb_bl," g/dL"),
        paste0(input$plt_bl," ×10⁹/L"),
        paste0(input$retic_bl,"%"),
        input$dat_status,"High","High","Low","High",
        paste0(input$treg_pct,"%"),
        paste0(input$b_cell_pct,"%")
      ),
      Status = c(
        ifelse(input$hgb_bl < 12, "⬇ Low","✓ Normal"),
        ifelse(input$plt_bl < 100,"⬇ Low","✓ Normal"),
        ifelse(input$retic_bl > 3,"⬆ High","✓ Normal"),
        ifelse(input$dat_status!="Negative","⚠ Positive","✓ Negative"),
        "⚠ Elevated","⚠ Elevated","⬇ Low","⚠ Elevated",
        ifelse(input$treg_pct < 3, "⬇ Low","✓ Normal"),
        ifelse(input$b_cell_pct > 20,"⬆ High","✓ Normal")
      )
    )
  }, options=list(pageLength=10, dom="t"))

  ## ----- PK Simulation -----
  pk_data <- eventReactive(input$run_pk, {
    p <- list(
      PRED_dose  = if("Prednisone" %in% input$pk_drugs) input$pred_dose_pk else 0,
      MMF_dose   = if("MMF"   %in% input$pk_drugs) input$mmf_dose_pk else 0,
      ELT_dose   = if("Eltrombopag" %in% input$pk_drugs) input$elt_dose_pk else 0,
      SIRO_dose  = if("Sirolimus"   %in% input$pk_drugs) input$siro_dose_pk else 0
    )
    rtx_ev <- if("Rituximab" %in% input$pk_drugs)
      ev(time=c(1,8,15,22), amt=input$rtx_bsa_pk*375, cmt="CENT_RTX") else NULL
    ivig_ev <- if("IVIg" %in% input$pk_drugs)
      ev(time=c(input$ivig_day_pk, input$ivig_day_pk+1),
         amt=input$ivig_total_pk*input$weight*500, cmt="CENT_IVIG") else NULL
    run_sim(p, rtx_ev, ivig_ev, input$sim_days_pk)
  })

  make_pk_plot <- function(data_r, col, ylab, title, color="#1E88E5") {
    df <- data_r(); if(is.null(df)) return(NULL)
    ggplot(df, aes_string("time", col)) + geom_line(color=color, size=1.1) +
      labs(title=title, x="Day", y=ylab) + theme_bw(base_size=10)
  }
  output$pk_pred_plot  <- renderPlot(make_pk_plot(pk_data,"C_Pred_pl","Prednisolone (mg/L)","Prednisolone","#FB8C00"))
  output$pk_rtx_plot   <- renderPlot(make_pk_plot(pk_data,"C_Rtx_pl","Rituximab (mg/L)","Rituximab","#1E88E5"))
  output$pk_mpa_plot   <- renderPlot(make_pk_plot(pk_data,"C_MPA_pl","MPA (mg/L)","MPA (MMF)","#43A047"))
  output$pk_elt_plot   <- renderPlot(make_pk_plot(pk_data,"C_Elt_pl","Eltrombopag (mg/L)","Eltrombopag","#8E24AA"))
  output$pk_siro_plot  <- renderPlot(make_pk_plot(pk_data,"C_Siro_wh","Sirolimus (ng/mL)","Sirolimus","#00897B"))
  output$pk_ivig_plot  <- renderPlot(make_pk_plot(pk_data,"C_IVIg_pl","IVIg (mg/L)","IVIg","#E53935"))

  output$pk_params_table <- DT::renderDataTable({
    data.frame(
      Drug = c("Prednisolone","Rituximab","MPA","Eltrombopag","Sirolimus","IVIg"),
      F_percent = c("82","100","94","52","14","100"),
      t_half = c("2.5 h","22 days","17 h","26 h","62 h","21 days"),
      Vd = c("0.97 L/kg","0.06 L/kg","0.21 L/kg","0.29 L/kg","6.9 L/kg","0.05 L/kg"),
      Target = c("GR (nuclear)","CD20 (B cell)","IMPDH","c-Mpl (TPO-R)","mTORC1","FcRn/FcγR")
    )
  }, options=list(dom="t"))

  ## ----- Hematologic Response -----
  hem_data <- eventReactive(input$run_hem, {
    p <- list(
      PRED_dose = input$pred_dose_hem,
      MMF_dose  = input$mmf_dose_hem,
      ELT_dose  = input$elt_dose_hem,
      SIRO_dose = input$siro_dose_hem,
      SPLEN     = as.integer(input$splen_hem)
    )
    rtx_ev <- if(input$rtx_hem) ev(time=c(1,8,15,22), amt=640, cmt="CENT_RTX") else NULL
    run_sim(p, rtx_ev, NULL, input$sim_days_hem)
  })

  output$hem_hgb_plot <- renderPlot({
    df <- hem_data()
    ggplot(df, aes(time, Hgb)) + geom_line(color="#C62828", size=1.2) +
      geom_hline(yintercept=c(8,10,12), linetype="dashed",
                 color=c("#EF5350","#FF9800","#43A047")) +
      annotate("text",x=max(df$time)*0.95,y=c(8.2,10.2,12.2),
               label=c("Severe","Moderate","CR"),size=3,color=c("#EF5350","#FF9800","#43A047")) +
      labs(title="Hemoglobin Response",x="Day",y="Hgb (g/dL)") +
      theme_bw(base_size=11)
  })
  output$hem_plt_plot <- renderPlot({
    df <- hem_data()
    ggplot(df, aes(time, Plt)) + geom_line(color="#0D47A1", size=1.2) +
      geom_hline(yintercept=c(20,50,100), linetype="dashed",
                 color=c("#EF5350","#FF9800","#43A047")) +
      annotate("text",x=max(df$time)*0.95,y=c(22,52,103),
               label=c("ICH risk","Moderate","CR"),size=3,color=c("#EF5350","#FF9800","#43A047")) +
      labs(title="Platelet Count Response",x="Day",y="Plt (×10⁹/L)") +
      theme_bw(base_size=11)
  })
  output$hem_retic_plot <- renderPlot({
    df <- hem_data()
    ggplot(df, aes(time, Retic_pct)) + geom_line(color="#6A1B9A", size=1.2) +
      geom_hline(yintercept=2, linetype="dashed", color="gray40") +
      labs(title="Reticulocyte %",x="Day",y="Reticulocyte (%)") +
      theme_bw(base_size=11)
  })
  output$hem_mk_plot <- renderPlot({
    df <- hem_data()
    ggplot(df, aes(time, MK_pool)) + geom_line(color="#00897B", size=1.2) +
      geom_hline(yintercept=100, linetype="dashed", color="gray40") +
      labs(title="Megakaryocyte Pool",x="Day",y="MK pool (AU)") +
      theme_bw(base_size=11)
  })
  output$hem_milestone_table <- DT::renderDataTable({
    df <- hem_data()
    milestone_days <- c(14, 28, 60, 90, 180, 365)
    milestone_days <- milestone_days[milestone_days <= max(df$time)]
    df %>% filter(time %in% milestone_days) %>%
      mutate(
        AIHA_Resp = case_when(Hgb>=12 ~ "CR",Hgb>=10 ~ "PR", TRUE ~ "NR"),
        ITP_Resp  = case_when(Plt>=100 ~ "CR",Plt>=50  ~ "PR", TRUE ~ "NR"),
        Composite = ifelse(Hgb>=12 & Plt>=100, "CR","Non-CR")
      ) %>%
      select(Day=time, Hgb, Plt, Retic_pct, AIHA_Resp, ITP_Resp, Composite) %>%
      mutate(across(c(Hgb,Plt,Retic_pct), ~round(.,1)))
  }, options=list(dom="t"))

  ## ----- Immunologic Markers -----
  imm_data <- eventReactive(input$run_imm, {
    p <- list(
      PRED_dose = input$pred_dose_imm,
      MMF_dose  = input$mmf_dose_imm,
      SIRO_dose = input$siro_dose_imm
    )
    rtx_ev <- if(input$rtx_imm) ev(time=c(1,8,15,22), amt=640, cmt="CENT_RTX") else NULL
    run_sim(p, rtx_ev, NULL, input$sim_days_imm)
  })

  output$imm_ab_rbc_plot <- renderPlot({
    df <- imm_data()
    ggplot(df, aes(time, Ab_RBC_lvl)) + geom_line(color="#C62828", size=1.1) +
      geom_hline(yintercept=0.1, linetype="dashed") +
      labs(title="Anti-RBC IgG Autoantibody",x="Day",y="mg/L") +
      theme_bw(base_size=10)
  })
  output$imm_ab_plt_plot <- renderPlot({
    df <- imm_data()
    ggplot(df, aes(time, Ab_Plt_lvl)) + geom_line(color="#0D47A1", size=1.1) +
      geom_hline(yintercept=0.1, linetype="dashed") +
      labs(title="Anti-Plt IgG Autoantibody",x="Day",y="mg/L") +
      theme_bw(base_size=10)
  })
  output$imm_treg_plot <- renderPlot({
    df <- imm_data()
    ggplot(df, aes(time, Treg_count)) + geom_line(color="#1B5E20", size=1.1) +
      geom_hline(yintercept=80, linetype="dashed") +
      labs(title="Treg Cells",x="Day",y="cells/µL") +
      theme_bw(base_size=10)
  })
  output$imm_bcell_plot <- renderPlot({
    df <- imm_data()
    ggplot(df, aes(time, B_RBC_cnt+B_Plt_cnt)) + geom_line(color="#B71C1C", size=1.1) +
      labs(title="Total Autoreactive B Cells",x="Day",y="cells/µL") +
      theme_bw(base_size=10)
  })
  output$imm_c3b_plot <- renderPlot({
    df <- imm_data()
    ggplot(df, aes(time, C3b_lvl)) + geom_line(color="#4A148C", size=1.1) +
      geom_hline(yintercept=0.1, linetype="dashed") +
      labs(title="Complement C3b (on RBC surface)",x="Day",y="AU") +
      theme_bw(base_size=10)
  })
  output$imm_bm_ery_plot <- renderPlot({
    df <- imm_data()
    ggplot(df, aes(time, Retic_pct)) + geom_line(color="#E65100", size=1.1) +
      geom_hline(yintercept=2, linetype="dashed") +
      labs(title="Reticulocyte (BM Compensatory Response)",x="Day",y="Retic %") +
      theme_bw(base_size=10)
  })

  ## ----- Scenario Comparison -----
  sc_params <- list(
    sc1 = list(PRED_dose=0,   MMF_dose=0, ELT_dose=0,  SIRO_dose=0, SPLEN=0, rtx=FALSE, label="Untreated"),
    sc2 = list(PRED_dose=1.5, MMF_dose=0, ELT_dose=0,  SIRO_dose=0, SPLEN=0, rtx=FALSE, label="Pred 1.5mg/kg"),
    sc3 = list(PRED_dose=0.5, MMF_dose=0, ELT_dose=0,  SIRO_dose=0, SPLEN=0, rtx=TRUE,  label="Pred+Rituximab"),
    sc4 = list(PRED_dose=0.25,MMF_dose=2, ELT_dose=0,  SIRO_dose=0, SPLEN=0, rtx=FALSE, label="Pred+MMF"),
    sc5 = list(PRED_dose=0.1, MMF_dose=0, ELT_dose=50, SIRO_dose=4, SPLEN=0, rtx=FALSE, label="Siro+Elt"),
    sc6 = list(PRED_dose=0.1, MMF_dose=0, ELT_dose=0,  SIRO_dose=0, SPLEN=1, rtx=FALSE, label="Splenectomy")
  )
  sc_colors <- c(sc1="#E53935",sc2="#FB8C00",sc3="#1E88E5",sc4="#43A047",sc5="#8E24AA",sc6="#00897B")

  sc_data <- eventReactive(input$run_sc, {
    sel <- input$sc_select
    bind_rows(lapply(sel, function(s) {
      p  <- sc_params[[s]]; rtx_ev <- if(p$rtx) ev(time=c(1,8,15,22), amt=640, cmt="CENT_RTX") else NULL
      pv <- p[!names(p) %in% c("rtx","label")]
      run_sim(pv, rtx_ev, NULL, input$sim_days_sc) %>% mutate(scenario=p$label, sc_id=s)
    }))
  })

  output$sc_hgb_plot <- renderPlot({
    df <- sc_data(); if(nrow(df)==0) return(NULL)
    sel_colors <- sc_colors[unique(df$sc_id)]
    names(sel_colors) <- unique(df$scenario)
    ggplot(df, aes(time, Hgb, color=scenario)) + geom_line(size=1.1) +
      geom_hline(yintercept=12, linetype="dashed") +
      scale_color_manual(values=sel_colors) +
      labs(title="Hemoglobin – Scenario Comparison", x="Day", y="Hgb (g/dL)", color="") +
      theme_bw(base_size=11) + theme(legend.position="bottom")
  })
  output$sc_plt_plot <- renderPlot({
    df <- sc_data(); if(nrow(df)==0) return(NULL)
    sel_colors <- sc_colors[unique(df$sc_id)]; names(sel_colors) <- unique(df$scenario)
    ggplot(df, aes(time, Plt, color=scenario)) + geom_line(size=1.1) +
      geom_hline(yintercept=100, linetype="dashed") +
      scale_color_manual(values=sel_colors) +
      labs(title="Platelet Count – Scenario Comparison", x="Day", y="Plt (×10⁹/L)", color="") +
      theme_bw(base_size=11) + theme(legend.position="bottom")
  })
  output$sc_response_table <- DT::renderDataTable({
    df <- sc_data(); if(nrow(df)==0) return(NULL)
    bind_rows(lapply(c(180, 365), function(d) {
      df %>% filter(time == min(abs(time - d) + time)) %>%
        slice(1) %>%
        mutate(Day=d,
               Hgb=round(Hgb,1), Plt=round(Plt),
               AIHA=case_when(Hgb>=12~"CR",Hgb>=10~"PR",TRUE~"NR"),
               ITP =case_when(Plt>=100~"CR",Plt>=50~"PR",TRUE~"NR"),
               CR  = ifelse(Hgb>=12 & Plt>=100,"YES","NO")) %>%
        select(Day, Scenario=scenario, Hgb, Plt, AIHA, ITP, Composite_CR=CR)
    })) %>% arrange(Day, desc(Hgb))
  }, options=list(dom="t"))

  ## ----- Biomarker Dashboard -----
  output$bm_gauge_plot <- renderPlot({
    # Simple bar gauge
    scores <- c(
      Hgb    = max(0, min(100, (12 - input$bm_hgb) / 9 * 100)),
      Plt    = max(0, min(100, (150 - input$bm_plt) / 149 * 100)),
      Retic  = max(0, min(100, (input$bm_retic - 2) / 18 * 100)),
      LDH    = max(0, min(100, (input$bm_ldh - 250) / 1750 * 100)),
      Haptogl= max(0, min(100, (2 - input$bm_haptoglobin) / 2 * 100)),
      Treg   = max(0, min(100, (5 - input$bm_treg_pct) / 5 * 100))
    )
    df <- data.frame(Marker=names(scores), Score=as.numeric(scores))
    ggplot(df, aes(reorder(Marker, Score), Score, fill=Score)) +
      geom_col() +
      scale_fill_gradient2(low="#43A047", mid="#FF9800", high="#E53935",
                           midpoint=50, limits=c(0,100)) +
      coord_flip() +
      labs(title="Disease Activity Scores (0=normal, 100=severe)",
           x="", y="Activity Score (0–100)") +
      theme_bw(base_size=10) + theme(legend.position="none")
  })

  output$bm_radar_plot <- renderPlot({
    # Spider/radar using bar chart (simplified)
    vals <- data.frame(
      category = c("Hemolysis","Thrombocytopenia","Autoimmunity","BM Stress","Organ Load"),
      value    = c(
        mean(c(max(0,12-input$bm_hgb)/9, max(0,input$bm_ldh-250)/1750)),
        (max(0,150-input$bm_plt)/149),
        mean(c(input$bm_treg_pct/5, input$bm_bcell_pct/40)),
        max(0,input$bm_retic-2)/18,
        max(0,55-input$bm_bili)/55
      )
    ) %>% mutate(value=pmin(1, pmax(0, value)))
    ggplot(vals, aes(category, value, fill=category)) +
      geom_col(show.legend=FALSE) +
      scale_fill_brewer(palette="Set2") +
      coord_polar() +
      ylim(0, 1) +
      labs(title="Domain Burden Radar", x="", y="") +
      theme_bw(base_size=9) +
      theme(axis.text.x=element_text(size=8))
  })

  output$bm_summary_table <- DT::renderDataTable({
    data.frame(
      Biomarker  = c("Hgb","Plt","Reticulocyte","DAT","PAIgG","LDH",
                     "Haptoglobin","Indirect Bili","Treg%","B cell%"),
      Patient    = c(input$bm_hgb, input$bm_plt, input$bm_retic,
                     input$bm_dat, input$bm_paIgG, input$bm_ldh,
                     input$bm_haptoglobin, input$bm_bili,
                     input$bm_treg_pct, input$bm_bcell_pct),
      Normal     = c(">12 g/dL",">150 ×10⁹/L","0.5–2%","Negative",
                     "<20","<250 U/L","0.5–2.5 g/L","<20 µmol/L",
                     "5–10%","5–15%"),
      Status     = c(
        ifelse(input$bm_hgb < 8,"Severe Low",ifelse(input$bm_hgb<12,"Low","Normal")),
        ifelse(input$bm_plt < 20,"Critical",ifelse(input$bm_plt<100,"Low","Normal")),
        ifelse(input$bm_retic > 5,"High (brisk)","Normal"),
        ifelse(input$bm_dat=="Negative","Normal","Positive"),
        ifelse(input$bm_paIgG > 50,"Elevated","Normal"),
        ifelse(input$bm_ldh > 500,"High","Normal"),
        ifelse(input$bm_haptoglobin < 0.2,"Depleted","Normal"),
        ifelse(input$bm_bili > 30,"Elevated","Normal"),
        ifelse(input$bm_treg_pct < 3,"Low","Normal"),
        ifelse(input$bm_bcell_pct > 20,"Elevated","Normal")
      )
    )
  }, options=list(dom="t"))

  output$bm_treatment_recommend <- renderUI({
    severity_score <- sum(c(
      input$bm_hgb < 8,
      input$bm_plt < 30,
      input$bm_ldh > 600,
      input$bm_haptoglobin < 0.1
    ))
    if (severity_score >= 3) {
      rec <- tags$div(
        tags$h5("⚠ SEVERE – Urgent Intervention Required", style="color:red"),
        tags$ul(
          tags$li("IVIg 1g/kg IV (days 1–2) for rapid Plt rescue"),
          tags$li("High-dose methylprednisolone 1g IV ×3 days"),
          tags$li("Consider packed RBC transfusion if Hgb < 7 g/dL"),
          tags$li("Consider rituximab 375mg/m² ×4 after stabilization")
        )
      )
    } else if (severity_score >= 1) {
      rec <- tags$div(
        tags$h5("⚡ MODERATE – Active Treatment Needed", style="color:orange"),
        tags$ul(
          tags$li("Prednisone 1–2 mg/kg/day with slow taper"),
          tags$li("If relapsed after steroids: Rituximab 375mg/m² ×4"),
          tags$li("If ITP predominant: Add eltrombopag 50mg/day"),
          tags$li("If Treg < 3%: Consider sirolimus (mTOR↓ → Treg↑)")
        )
      )
    } else {
      rec <- tags$div(
        tags$h5("✓ MILD / MONITORING", style="color:green"),
        tags$ul(
          tags$li("Close monitoring every 2–4 weeks"),
          tags$li("Steroid taper / maintenance MMF if on treatment"),
          tags$li("Watch for relapse signs (fatigue, bruising)")
        )
      )
    }
    rec
  })
}

## ============================================================
## Launch App
## ============================================================
shinyApp(ui=ui, server=server)
