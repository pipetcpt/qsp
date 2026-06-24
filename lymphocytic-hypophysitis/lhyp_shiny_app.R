## =============================================================================
## Lymphocytic Hypophysitis (LyH) — Interactive Shiny Dashboard
## File: lhyp_shiny_app.R
##
## Tabs:
##   1. Patient Profile & Disease Settings
##   2. Drug PK (Prednisolone / Azathioprine / Rituximab)
##   3. Pituitary Function & Mass
##   4. Hormone Axes (HPA · HPT · HPG · GH)
##   5. Immune Dynamics (T/B cells, APA)
##   6. Clinical Endpoints & Scenario Comparison
##   7. Biomarker Trajectories
##
## Launch: shiny::runApp("lhyp_shiny_app.R")
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)
library(DT)

## ---- Model code (inline) ----------------------------------------------------
lhyp_code <- '
$PARAM
  ka=1.20, Vc=45.0, Vp=120.0, CL_P=8.50, Q_P=25.0, F_pred=0.82,
  ka_aza=0.80, Vc_aza=30.0, CL_aza=12.0,
  Vc_rtx=3.1, CL_rtx=0.015,
  kprolif_Tn=0.050, kdeath_Tn=0.020, Tn_baseline=1000,
  kact_T=0.12, kdeath_Te=0.08, kinact_T=0.05,
  kprolif_Tr=0.030, kdeath_Tr=0.025, kinduct_Tr=0.015,
  kprolif_Bn=0.030, kdeath_Bn=0.015, kact_B=0.04, kdeath_Bp=0.020,
  Bn_baseline=500,
  kprod_APA=0.025, kdeg_APA=0.008,
  kpit_inflam=0.015, kpit_repair=0.004, kpit_fibs=0.002,
  kprod_ACTH=0.50, kdeg_ACTH=0.35, kprod_Cort=2.00, kdeg_Cort=0.12,
  IC50_Cort=400, n_Hill_Cort=2.0, Emax_Cort=0.90,
  kprod_TSH=0.80, kdeg_TSH=0.15, kprod_fT4=0.50, kdeg_fT4=0.06,
  IC50_fT4=15.0, Emax_fT4=0.85,
  kprod_GH=1.50, kdeg_GH=0.40, kprod_IGF1=0.30, kdeg_IGF1=0.045, IC50_IGF1=200,
  kprod_FSH=0.30, kdeg_FSH=0.10, kprod_LH=0.40, kdeg_LH=0.20,
  kprod_E2=5.0, kdeg_E2=0.08, IC50_E2=200,
  kprod_PRL=0.50, kdeg_PRL=0.12, k_stalk_PRL=0.30,
  kprod_ADH=0.40, kdeg_ADH=0.20, k_post_loss=0.80,
  EC50_pred=0.50, Emax_pred=0.85,
  EC50_aza=0.20, Emax_aza=0.70,
  EC50_rtx=0.05, Emax_rtx=0.90,
  PitFunc_init=0.85, PitInf_init=0.30

$CMT
  Pred_gut Pred_central Pred_periph
  AZA_gut AZA_plasma
  RTX_plasma
  Tn Te Tr Bn Bp APA
  PitInf PitFunc
  ACTH Cortisol TSH fT4 GH IGF1 FSH LH E2 PRL ADH

$MAIN
  double Cpred = Pred_central / Vc;
  double Caza  = AZA_plasma / Vc_aza;
  double Crtx  = RTX_plasma / Vc_rtx;
  double E_pred = Emax_pred * Cpred / (EC50_pred + Cpred);
  double E_aza  = Emax_aza  * Caza  / (EC50_aza  + Caza );
  double E_rtx  = Emax_rtx  * Crtx  / (EC50_rtx  + Crtx );
  double ImmunoSupp = 1.0 - (1.0 - E_pred) * (1.0 - E_aza * 0.5);
  double Pfunc = (PitFunc > 1.0) ? 1.0 : (PitFunc < 0.0) ? 0.0 : PitFunc;
  double Pinfl = (PitInf  < 0.0) ? 0.0 : PitInf;
  double Cort_FB = 1.0 - Emax_Cort * pow(Cortisol,n_Hill_Cort) /
                         (pow(IC50_Cort,n_Hill_Cort)+pow(Cortisol,n_Hill_Cort));
  double ACTH_prod = kprod_ACTH * Cort_FB * Pfunc * (1.0 - E_pred * 0.70);
  double fT4_FB  = 1.0 - Emax_fT4 * fT4 / (IC50_fT4 + fT4);
  double TSH_prod = kprod_TSH * fT4_FB * Pfunc;
  double IGF1_FB  = 1.0 - 0.80 * IGF1 / (IC50_IGF1 + IGF1);
  double GH_prod  = kprod_GH * IGF1_FB * Pfunc;
  double E2_FB    = 1.0 - 0.80 * E2 / (IC50_E2 + E2);
  double FSH_prod = kprod_FSH * E2_FB * Pfunc;
  double LH_prod  = kprod_LH  * E2_FB * Pfunc;
  double stalk_factor = 1.0 + k_stalk_PRL * Pinfl;
  double ADH_prod = kprod_ADH * (1.0 - k_post_loss * Pinfl / (1.0 + Pinfl));
  double Imm_act  = (Te + APA * 0.01);
  double pit_damage   = kpit_inflam * Imm_act * (1.0 - ImmunoSupp);
  double pit_recovery = kpit_repair * Pfunc;

$ODE
  dxdt_Pred_gut    = -ka * Pred_gut;
  dxdt_Pred_central = ka * Pred_gut - (CL_P/Vc)*Pred_central - (Q_P/Vc)*Pred_central + (Q_P/Vp)*Pred_periph;
  dxdt_Pred_periph  = (Q_P/Vc)*Pred_central - (Q_P/Vp)*Pred_periph;
  dxdt_AZA_gut    = -ka_aza * AZA_gut;
  dxdt_AZA_plasma =  ka_aza * AZA_gut - (CL_aza/Vc_aza)*AZA_plasma;
  dxdt_RTX_plasma = -(CL_rtx/Vc_rtx)*RTX_plasma;
  dxdt_Tn = kprolif_Tn/24.0*Tn_baseline - kdeath_Tn/24.0*Tn - kact_T/24.0*Tn*(1.0-ImmunoSupp);
  dxdt_Te = kact_T/24.0*Tn*(1.0-ImmunoSupp) - kdeath_Te/24.0*Te*(1.0+E_pred*2.5) - kinact_T/24.0*Te;
  dxdt_Tr = kinduct_Tr/24.0*Tn*(1.0+E_pred*1.5) + kprolif_Tr/24.0*Tr - kdeath_Tr/24.0*Tr;
  dxdt_Bn = kprolif_Bn/24.0*Bn_baseline - kdeath_Bn/24.0*Bn - kact_B/24.0*Bn*(1.0-ImmunoSupp)*(1.0-E_rtx);
  dxdt_Bp = kact_B/24.0*Bn*(1.0-ImmunoSupp)*(1.0-E_rtx) - kdeath_Bp/24.0*Bp*(1.0+E_pred*1.5);
  dxdt_APA = kprod_APA/24.0*Bp - kdeg_APA/24.0*APA;
  dxdt_PitInf  = pit_damage/24.0 - kpit_fibs/24.0*Pinfl - (kpit_repair+E_pred*0.02)/24.0*Pinfl;
  dxdt_PitFunc = -pit_damage/24.0 + pit_recovery/24.0;
  dxdt_ACTH    = ACTH_prod   - kdeg_ACTH * ACTH;
  dxdt_Cortisol= kprod_Cort * ACTH - kdeg_Cort * Cortisol;
  dxdt_TSH     = TSH_prod    - kdeg_TSH * TSH;
  dxdt_fT4     = kprod_fT4 * TSH - kdeg_fT4 * fT4;
  dxdt_GH      = GH_prod     - kdeg_GH  * GH;
  dxdt_IGF1    = kprod_IGF1 * GH - kdeg_IGF1 * IGF1;
  dxdt_FSH     = FSH_prod    - kdeg_FSH * FSH;
  dxdt_LH      = LH_prod     - kdeg_LH  * LH;
  dxdt_E2      = kprod_E2 * FSH - kdeg_E2 * E2;
  dxdt_PRL     = kprod_PRL * stalk_factor - kdeg_PRL * PRL;
  dxdt_ADH     = ADH_prod    - kdeg_ADH * ADH;

$INIT
  Pred_gut=0, Pred_central=0, Pred_periph=0,
  AZA_gut=0, AZA_plasma=0, RTX_plasma=0,
  Tn=1000, Te=50, Tr=200, Bn=500, Bp=20, APA=10,
  PitInf=PitInf_init, PitFunc=PitFunc_init,
  ACTH=22.0, Cortisol=500, TSH=2.5, fT4=15.0,
  GH=2.0, IGF1=180, FSH=5.0, LH=5.0, E2=200,
  PRL=15.0, ADH=3.0

$TABLE
  capture Cpred = Pred_central / Vc;
  capture Caza  = AZA_plasma / Vc_aza;
  capture Crtx  = RTX_plasma / Vc_rtx;
  capture PFS   = (PitFunc > 1) ? 1.0 : (PitFunc < 0) ? 0.0 : PitFunc;
'

mod_base <- mcode("LyH_shiny", lhyp_code, quiet = TRUE)

## ---- Dark ggplot theme ------------------------------------------------------
theme_dark_qsp <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.background  = element_rect(fill = "#0d1117", color = NA),
      panel.background = element_rect(fill = "#161b22", color = NA),
      panel.grid.major = element_line(color = "#21262d", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      text             = element_text(color = "#c9d1d9"),
      axis.text        = element_text(color = "#8b949e"),
      strip.text       = element_text(color = "#58a6ff", face = "bold"),
      legend.background = element_rect(fill = "#0d1117", color = NA),
      legend.key        = element_rect(fill = NA),
      legend.text      = element_text(color = "#c9d1d9"),
      plot.title       = element_text(color = "#58a6ff", face = "bold", size = 12),
      plot.subtitle    = element_text(color = "#8b949e", size = 9)
    )
}

## ---- Reference ranges -------------------------------------------------------
ref_ranges <- list(
  ACTH     = c(6, 50),
  Cortisol = c(200, 700),
  TSH      = c(0.4, 4.0),
  fT4      = c(12, 22),
  GH       = c(0.5, 5),
  IGF1     = c(100, 300),
  FSH      = c(2, 12),
  LH       = c(2, 12),
  E2       = c(100, 400),
  PRL      = c(5, 25),
  ADH      = c(1, 5)
)

## ============================================================================
## UI
## ============================================================================
ui <- dashboardPage(
  skin = "black",

  dashboardHeader(
    title = "LyH QSP Dashboard"
  ),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",         tabName = "tab_patient",   icon = icon("user-md")),
      menuItem("Drug PK",                 tabName = "tab_pk",        icon = icon("pills")),
      menuItem("Pituitary Function",      tabName = "tab_pituitary", icon = icon("brain")),
      menuItem("Hormone Axes",            tabName = "tab_hormones",  icon = icon("chart-line")),
      menuItem("Immune Dynamics",         tabName = "tab_immune",    icon = icon("shield-virus")),
      menuItem("Clinical Endpoints",      tabName = "tab_endpoints", icon = icon("stethoscope")),
      menuItem("Biomarker Trajectories",  tabName = "tab_biomarker", icon = icon("dna"))
    )
  ),

  dashboardBody(
    tags$style(HTML("
      .content-wrapper, .right-side { background-color: #0d1117; }
      .box { background-color: #161b22; border-top: 3px solid #58a6ff; color: #c9d1d9; }
      .box-header { background-color: #21262d; color: #58a6ff; }
      .main-header .logo, .main-header .navbar { background-color: #21262d; }
      .sidebar { background-color: #161b22; }
      .slider-label, .control-label { color: #c9d1d9 !important; }
      h4 { color: #58a6ff; }
    ")),

    tabItems(

      ## ========================== TAB 1: PATIENT PROFILE =====================
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Disease Severity at Onset", width = 6, solidHeader = TRUE,
            sliderInput("pit_func_init", "Initial Pituitary Function (%)",
                        min = 20, max = 100, value = 85, step = 5),
            sliderInput("pit_inf_init", "Initial Inflammatory Mass (relative)",
                        min = 0, max = 2, value = 0.3, step = 0.1),
            radioButtons("disease_type", "LyH Subtype",
                         choices = c("Lymphocytic Adenohypophysitis (LAH)" = "lah",
                                     "Lymphocytic Infundibuloneurohypophysitis (LINH)" = "linh",
                                     "Lymphocytic Panhypophysitis (LPH)" = "lph",
                                     "ICI-associated LyH (Checkpoint inhibitor)" = "ici"),
                         selected = "lah"),
            numericInput("sim_days", "Simulation Duration (days)", value = 365, min = 30, max = 730)
          ),
          box(title = "Treatment Selection", width = 6, solidHeader = TRUE,
            checkboxGroupInput("treatments", "Active Treatments",
                               choices = c("Prednisolone (oral)" = "pred",
                                           "Azathioprine (steroid-sparing)" = "aza",
                                           "Rituximab (anti-CD20)" = "rtx"),
                               selected = "pred"),
            conditionalPanel("input.treatments.includes('pred')",
              sliderInput("pred_dose", "Prednisolone Initial Dose (mg/day)", min = 5, max = 80, value = 60),
              sliderInput("pred_taper_wk", "Taper Duration (weeks)", min = 4, max = 26, value = 12),
              sliderInput("pred_maint", "Maintenance Dose (mg/day)", min = 2.5, max = 15, value = 5)
            ),
            conditionalPanel("input.treatments.includes('aza')",
              sliderInput("aza_dose", "Azathioprine Dose (mg/day)", min = 50, max = 200, value = 150),
              numericInput("aza_start", "Azathioprine Start (day)", value = 30, min = 0, max = 90)
            ),
            conditionalPanel("input.treatments.includes('rtx')",
              sliderInput("rtx_dose", "Rituximab Dose (mg/infusion)", min = 500, max = 1000, value = 1000),
              numericInput("rtx_n", "Number of Infusions", value = 2, min = 1, max = 4)
            ),
            actionButton("run_sim", "Run Simulation", class = "btn-primary", icon = icon("play"))
          )
        ),
        fluidRow(
          box(title = "Disease Overview: Lymphocytic Hypophysitis", width = 12, solidHeader = TRUE,
            p("Lymphocytic Hypophysitis (LyH) is a rare autoimmune condition characterized by lymphocytic infiltration
              of the pituitary gland, leading to progressive hypopituitarism. It predominantly affects women in the
              peripartum period and may also be triggered by immune checkpoint inhibitors (anti-CTLA4, anti-PD1)."),
            p("The model simulates multi-axis pituitary hormone deficiency (ACTH, TSH, FSH/LH, GH, ADH) and their
              response to immunosuppressive therapy (prednisolone, azathioprine, rituximab)."),
            tags$ul(
              tags$li("ACTH deficiency → secondary adrenal insufficiency (most common, ~70%)"),
              tags$li("TSH deficiency → central hypothyroidism (~40%)"),
              tags$li("FSH/LH deficiency → hypogonadism (~40%)"),
              tags$li("GH deficiency → adult GHD (~30%)"),
              tags$li("Stalk compression → hyperprolactinemia (~50%)"),
              tags$li("ADH deficiency → diabetes insipidus (in posterior LyH, ~25%)")
            )
          )
        )
      ),

      ## ========================== TAB 2: DRUG PK =============================
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Prednisolone PK Profile", width = 6, solidHeader = TRUE,
            plotlyOutput("plot_pred_pk", height = 350)
          ),
          box(title = "Azathioprine & Rituximab PK", width = 6, solidHeader = TRUE,
            plotlyOutput("plot_other_pk", height = 350)
          )
        ),
        fluidRow(
          box(title = "PK Parameters Summary", width = 12, solidHeader = TRUE,
            DTOutput("table_pk_params")
          )
        )
      ),

      ## ========================== TAB 3: PITUITARY ===========================
      tabItem(tabName = "tab_pituitary",
        fluidRow(
          box(title = "Pituitary Function Score (0–100%)", width = 6, solidHeader = TRUE,
            plotlyOutput("plot_pit_func", height = 350)
          ),
          box(title = "Pituitary Inflammatory Mass (MRI volume)", width = 6, solidHeader = TRUE,
            plotlyOutput("plot_pit_inf", height = 350)
          )
        ),
        fluidRow(
          box(title = "Pituitary Hormone Deficiency Progression", width = 12, solidHeader = TRUE,
            plotlyOutput("plot_hormone_deficiency", height = 400)
          )
        )
      ),

      ## ========================== TAB 4: HORMONE AXES ========================
      tabItem(tabName = "tab_hormones",
        fluidRow(
          box(title = "HPA Axis: ACTH & Cortisol", width = 6, solidHeader = TRUE,
            plotlyOutput("plot_hpa", height = 300)
          ),
          box(title = "HPT Axis: TSH & Free T4", width = 6, solidHeader = TRUE,
            plotlyOutput("plot_hpt", height = 300)
          )
        ),
        fluidRow(
          box(title = "GH / IGF-1 Axis", width = 6, solidHeader = TRUE,
            plotlyOutput("plot_gh", height = 300)
          ),
          box(title = "HPG Axis: FSH, LH, Estradiol", width = 6, solidHeader = TRUE,
            plotlyOutput("plot_hpg", height = 300)
          )
        )
      ),

      ## ========================== TAB 5: IMMUNE DYNAMICS =====================
      tabItem(tabName = "tab_immune",
        fluidRow(
          box(title = "T Cell Populations", width = 6, solidHeader = TRUE,
            plotlyOutput("plot_tcells", height = 350)
          ),
          box(title = "B Cells & Anti-Pituitary Antibodies", width = 6, solidHeader = TRUE,
            plotlyOutput("plot_bcells", height = 350)
          )
        ),
        fluidRow(
          box(title = "Immune Activation Score vs Pituitary Damage", width = 12, solidHeader = TRUE,
            plotlyOutput("plot_immune_pit", height = 350)
          )
        )
      ),

      ## ========================== TAB 6: CLINICAL ENDPOINTS ==================
      tabItem(tabName = "tab_endpoints",
        fluidRow(
          box(title = "Clinical Endpoint Radar Chart", width = 6, solidHeader = TRUE,
            plotlyOutput("plot_radar", height = 400)
          ),
          box(title = "Adrenal Crisis Risk & Cortisol Reserve", width = 6, solidHeader = TRUE,
            plotlyOutput("plot_adrenal_risk", height = 400)
          )
        ),
        fluidRow(
          box(title = "Summary Table at Key Timepoints", width = 12, solidHeader = TRUE,
            DTOutput("table_endpoints")
          )
        )
      ),

      ## ========================== TAB 7: BIOMARKER TRAJECTORIES ==============
      tabItem(tabName = "tab_biomarker",
        fluidRow(
          box(title = "Select Biomarker", width = 3, solidHeader = TRUE,
            selectInput("biomarker_select", "Biomarker",
                        choices = c("ACTH (pg/mL)" = "ACTH",
                                    "Cortisol (nmol/L)" = "Cortisol",
                                    "TSH (mIU/L)" = "TSH",
                                    "Free T4 (pmol/L)" = "fT4",
                                    "GH (ng/mL)" = "GH",
                                    "IGF-1 (ng/mL)" = "IGF1",
                                    "FSH (IU/L)" = "FSH",
                                    "LH (IU/L)" = "LH",
                                    "Estradiol (pmol/L)" = "E2",
                                    "Prolactin (ng/mL)" = "PRL",
                                    "ADH (pg/mL)" = "ADH",
                                    "Effector T cells (AU)" = "Te",
                                    "Regulatory T cells (AU)" = "Tr",
                                    "APA (AU)" = "APA",
                                    "Pituitary Function Score" = "PFS"),
                        selected = "Cortisol"),
            sliderInput("ref_alpha", "Reference Range Opacity", min = 0, max = 1, value = 0.15, step = 0.05),
            checkboxInput("show_ref", "Show Reference Range", value = TRUE),
            downloadButton("download_data", "Download CSV", class = "btn-sm")
          ),
          box(title = "Biomarker Trajectory", width = 9, solidHeader = TRUE,
            plotlyOutput("plot_biomarker", height = 500)
          )
        ),
        fluidRow(
          box(title = "Prolactin & ADH (Posterior Involvement)", width = 12, solidHeader = TRUE,
            plotlyOutput("plot_prl_adh", height = 350)
          )
        )
      )
    )  # end tabItems
  )    # end dashboardBody
)      # end dashboardPage

## ============================================================================
## SERVER
## ============================================================================
server <- function(input, output, session) {

  ## ---- Reactive: run simulation -------------------------------------------
  sim_result <- eventReactive(input$run_sim, {
    withProgress(message = "Running QSP simulation...", value = 0.3, {
      mod <- mod_base %>%
        param(PitFunc_init = input$pit_func_init / 100,
              PitInf_init  = input$pit_inf_init)

      e_list <- list()

      if ("pred" %in% input$treatments) {
        dose_init <- input$pred_dose
        n_taper   <- max(1, round(input$pred_taper_wk / 2))
        taper_seq <- seq(dose_init * 0.7, input$pred_maint, length.out = n_taper)
        e <- ev(ID = 1, time = 0, amt = dose_init * 0.82 / 2, cmt = "Pred_gut", evid = 1,
                ii = 12, addl = 55)
        start_t <- 28 * 24
        for (i in seq_along(taper_seq)) {
          e_tp <- ev(ID = 1, time = start_t + (i - 1) * 14 * 24,
                     amt = taper_seq[i] * 0.82 / 2, cmt = "Pred_gut", evid = 1,
                     ii = 12, addl = 27)
          e <- e + e_tp
        }
        maint_start <- (28 + input$pred_taper_wk * 7) * 24
        maint_dur   <- input$sim_days * 24 - maint_start
        if (maint_dur > 0) {
          e_maint <- ev(ID = 1, time = maint_start,
                        amt = input$pred_maint * 0.82 / 2, cmt = "Pred_gut", evid = 1,
                        ii = 12, addl = floor(maint_dur / 12) - 1)
          e <- e + e_maint
        }
        e_list <- c(e_list, list(e))
      }

      if ("aza" %in% input$treatments) {
        start_h <- input$aza_start * 24
        dur_h   <- (input$sim_days - input$aza_start) * 24
        if (dur_h > 0) {
          e_aza <- ev(ID = 1, time = start_h,
                      amt  = input$aza_dose * 0.5,
                      cmt  = "AZA_gut", evid = 1,
                      ii   = 24, addl = floor(dur_h / 24) - 1)
          e_list <- c(e_list, list(e_aza))
        }
      }

      if ("rtx" %in% input$treatments) {
        interval_h <- 180 * 24
        e_rtx_list <- lapply(seq_len(input$rtx_n), function(i) {
          ev(ID = 1, time = (i - 1) * interval_h,
             amt = input$rtx_dose, cmt = "RTX_plasma", evid = 1)
        })
        e_rtx <- Reduce("+", e_rtx_list)
        e_list <- c(e_list, list(e_rtx))
      }

      events_all <- if (length(e_list) > 0) Reduce("+", e_list) else NULL

      incProgress(0.4)

      out <- if (!is.null(events_all)) {
        mod %>% mrgsim(events = events_all, end = input$sim_days * 24, delta = 4)
      } else {
        mod %>% mrgsim(end = input$sim_days * 24, delta = 4)
      }

      incProgress(0.3)
      as.data.frame(out) %>% mutate(time_d = time / 24)
    })
  }, ignoreNULL = FALSE)

  ## ---- Plot helper ---------------------------------------------------------
  make_plotly <- function(gg) ggplotly(gg) %>%
    layout(paper_bgcolor = "#0d1117", plot_bgcolor = "#161b22",
           font = list(color = "#c9d1d9"),
           legend = list(bgcolor = "#161b22", font = list(color = "#c9d1d9")))

  ## ---- Tab 2: PK plots ----
  output$plot_pred_pk <- renderPlotly({
    df <- sim_result()
    p <- df %>%
      ggplot(aes(x = time_d, y = Cpred)) +
      geom_line(color = "#58a6ff", linewidth = 1) +
      labs(title = "Prednisolone Plasma Concentration",
           x = "Time (days)", y = "Cpred (mg/L)") +
      theme_dark_qsp()
    make_plotly(p)
  })

  output$plot_other_pk <- renderPlotly({
    df <- sim_result() %>%
      select(time_d, Caza, Crtx) %>%
      pivot_longer(c(Caza, Crtx), names_to = "drug", values_to = "conc") %>%
      mutate(drug = recode(drug, Caza = "Azathioprine (mg/L)", Crtx = "Rituximab (mg/L)"))
    p <- df %>%
      ggplot(aes(x = time_d, y = conc, color = drug)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c("Azathioprine (mg/L)" = "#3fb950",
                                     "Rituximab (mg/L)"     = "#ffa657")) +
      labs(title = "Azathioprine & Rituximab PK", x = "Time (days)", y = "Concentration (mg/L)") +
      theme_dark_qsp()
    make_plotly(p)
  })

  output$table_pk_params <- renderDT({
    tibble(
      Parameter = c("Prednisolone ka", "Prednisolone Vc", "Prednisolone CL", "Prednisolone Vp",
                    "Azathioprine ka", "Azathioprine Vc", "Azathioprine CL",
                    "Rituximab Vc", "Rituximab CL"),
      Value = c("1.20 h⁻¹", "45 L", "8.5 L/h", "120 L",
                "0.80 h⁻¹", "30 L", "12 L/h",
                "3.1 L", "0.015 L/h"),
      Note = c("Oral absorption", "Central volume", "Total clearance (CYP3A4)", "Tissue volume",
               "GI absorption", "Volume", "Includes TPMT",
               "Monoclonal Ab", "IgG catabolism")
    ) %>% datatable(options = list(pageLength = 9, dom = "t"),
                    class = "table-dark table-bordered table-sm")
  })

  ## ---- Tab 3: Pituitary ----
  output$plot_pit_func <- renderPlotly({
    df <- sim_result()
    p <- df %>%
      ggplot(aes(x = time_d, y = PFS * 100)) +
      geom_line(color = "#58a6ff", linewidth = 1.2) +
      geom_hline(yintercept = 70, linetype = "dashed", color = "#f85149", alpha = 0.7) +
      annotate("text", x = max(df$time_d) * 0.5, y = 68,
               label = "Hypopituitarism Threshold", color = "#f85149", size = 3) +
      labs(title = "Pituitary Function Score",
           x = "Time (days)", y = "Pituitary Function (%)") +
      ylim(0, 110) +
      theme_dark_qsp()
    make_plotly(p)
  })

  output$plot_pit_inf <- renderPlotly({
    df <- sim_result()
    p <- df %>%
      ggplot(aes(x = time_d, y = PitInf)) +
      geom_line(color = "#f85149", linewidth = 1.2) +
      geom_area(fill = "#f85149", alpha = 0.15) +
      labs(title = "Pituitary Inflammatory Mass (relative)",
           x = "Time (days)", y = "Mass Index (0 = normal)") +
      theme_dark_qsp()
    make_plotly(p)
  })

  output$plot_hormone_deficiency <- renderPlotly({
    df <- sim_result() %>%
      mutate(
        ACTH_def  = ifelse(ACTH < 10, 1, 0),
        TSH_def   = ifelse(TSH < 0.4, 1, 0),
        FSH_def   = ifelse(FSH < 2.0, 1, 0),
        GH_def    = ifelse(GH  < 0.5, 1, 0),
        PRL_hi    = ifelse(PRL > 25, 1, 0),
        ADH_def   = ifelse(ADH < 1.0, 1, 0)
      ) %>%
      select(time_d, ACTH_def, TSH_def, FSH_def, GH_def, PRL_hi, ADH_def) %>%
      pivot_longer(-time_d, names_to = "deficiency", values_to = "flag") %>%
      mutate(deficiency = recode(deficiency,
               ACTH_def = "ACTH Deficiency", TSH_def = "TSH Deficiency",
               FSH_def  = "FSH Deficiency",  GH_def  = "GH Deficiency",
               PRL_hi   = "Hyperprolactinemia", ADH_def = "ADH Deficiency"))
    p <- df %>%
      ggplot(aes(x = time_d, y = deficiency, fill = deficiency, alpha = flag)) +
      geom_tile(height = 0.8) +
      scale_alpha_continuous(range = c(0, 0.85), guide = "none") +
      scale_fill_brewer(palette = "Set1", guide = "none") +
      labs(title = "Hormone Deficiency Heatmap Over Time",
           x = "Time (days)", y = "") +
      theme_dark_qsp()
    make_plotly(p)
  })

  ## ---- Tab 4: Hormone axes ----
  output$plot_hpa <- renderPlotly({
    df <- sim_result() %>% select(time_d, ACTH, Cortisol) %>%
      pivot_longer(-time_d, names_to = "h", values_to = "v")
    p <- ggplot(df, aes(x = time_d, y = v, color = h)) +
      geom_line(linewidth = 1) +
      facet_wrap(~h, scales = "free_y") +
      scale_color_manual(values = c(ACTH = "#ffa657", Cortisol = "#ff7b72")) +
      labs(title = "HPA Axis: ACTH & Cortisol", x = "Time (days)", y = "Concentration") +
      theme_dark_qsp()
    make_plotly(p)
  })

  output$plot_hpt <- renderPlotly({
    df <- sim_result() %>% select(time_d, TSH, fT4) %>%
      pivot_longer(-time_d, names_to = "h", values_to = "v")
    p <- ggplot(df, aes(x = time_d, y = v, color = h)) +
      geom_line(linewidth = 1) +
      facet_wrap(~h, scales = "free_y") +
      scale_color_manual(values = c(TSH = "#d2a8ff", fT4 = "#79c0ff")) +
      labs(title = "HPT Axis: TSH & Free T4", x = "Time (days)", y = "Level") +
      theme_dark_qsp()
    make_plotly(p)
  })

  output$plot_gh <- renderPlotly({
    df <- sim_result() %>% select(time_d, GH, IGF1) %>%
      pivot_longer(-time_d, names_to = "h", values_to = "v")
    p <- ggplot(df, aes(x = time_d, y = v, color = h)) +
      geom_line(linewidth = 1) +
      facet_wrap(~h, scales = "free_y") +
      scale_color_manual(values = c(GH = "#3fb950", IGF1 = "#56d364")) +
      labs(title = "GH / IGF-1 Axis", x = "Time (days)", y = "Level") +
      theme_dark_qsp()
    make_plotly(p)
  })

  output$plot_hpg <- renderPlotly({
    df <- sim_result() %>% select(time_d, FSH, LH, E2) %>%
      pivot_longer(-time_d, names_to = "h", values_to = "v")
    p <- ggplot(df, aes(x = time_d, y = v, color = h)) +
      geom_line(linewidth = 1) +
      facet_wrap(~h, scales = "free_y") +
      scale_color_manual(values = c(FSH = "#ffa657", LH = "#f0883e", E2 = "#ff7b72")) +
      labs(title = "HPG Axis: FSH, LH, Estradiol", x = "Time (days)", y = "Level") +
      theme_dark_qsp()
    make_plotly(p)
  })

  ## ---- Tab 5: Immune dynamics ----
  output$plot_tcells <- renderPlotly({
    df <- sim_result() %>% select(time_d, Tn, Te, Tr) %>%
      pivot_longer(-time_d, names_to = "cell", values_to = "count") %>%
      mutate(cell = recode(cell, Tn = "Naive T", Te = "Effector T", Tr = "Regulatory T"))
    p <- ggplot(df, aes(x = time_d, y = count, color = cell)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c("Naive T" = "#58a6ff", "Effector T" = "#f85149", "Regulatory T" = "#3fb950")) +
      labs(title = "T Cell Populations", x = "Time (days)", y = "Count (AU)") +
      theme_dark_qsp()
    make_plotly(p)
  })

  output$plot_bcells <- renderPlotly({
    df <- sim_result() %>% select(time_d, Bn, Bp, APA) %>%
      pivot_longer(-time_d, names_to = "cell", values_to = "level") %>%
      mutate(cell = recode(cell, Bn = "Naive B", Bp = "Plasma Cells", APA = "Anti-Pituitary Abs"))
    p <- ggplot(df, aes(x = time_d, y = level, color = cell)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c("Naive B" = "#d2a8ff", "Plasma Cells" = "#ffa657",
                                     "Anti-Pituitary Abs" = "#ff7b72")) +
      labs(title = "B Cells & APA", x = "Time (days)", y = "Level (AU)") +
      theme_dark_qsp()
    make_plotly(p)
  })

  output$plot_immune_pit <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x = Te + APA * 0.01, y = PFS * 100, color = time_d)) +
      geom_path(linewidth = 1, arrow = arrow(length = unit(0.15, "cm"))) +
      scale_color_gradient(low = "#f85149", high = "#3fb950",
                            name = "Time (days)") +
      labs(title = "Immune Activation vs Pituitary Function (Phase Plot)",
           x = "Immune Activation Score (Te + APA×0.01)", y = "Pituitary Function (%)") +
      theme_dark_qsp()
    make_plotly(p)
  })

  ## ---- Tab 6: Clinical Endpoints ----
  output$plot_adrenal_risk <- renderPlotly({
    df <- sim_result() %>%
      mutate(adrenal_risk = ifelse(Cortisol < 200, 1, ifelse(Cortisol < 350, 0.5, 0)))
    p <- ggplot(df, aes(x = time_d, y = Cortisol)) +
      geom_line(color = "#58a6ff", linewidth = 1) +
      geom_ribbon(aes(ymin = 0, ymax = 200), fill = "#f85149", alpha = 0.2) +
      geom_ribbon(aes(ymin = 200, ymax = 350), fill = "#ffa657", alpha = 0.15) +
      geom_hline(yintercept = 200, linetype = "dashed", color = "#f85149") +
      geom_hline(yintercept = 350, linetype = "dashed", color = "#ffa657") +
      annotate("text", x = max(df$time_d) * 0.8, y = 160, label = "Adrenal Crisis Risk Zone", color = "#f85149", size = 3) +
      labs(title = "Cortisol Reserve & Adrenal Crisis Risk",
           x = "Time (days)", y = "Cortisol (nmol/L)") +
      theme_dark_qsp()
    make_plotly(p)
  })

  output$plot_radar <- renderPlotly({
    df_last <- sim_result() %>% filter(time_d == max(time_d)) %>% slice(1)
    cats <- c("HPA Axis", "HPT Axis", "GH Axis", "HPG Axis", "PRL/ADH", "Pit. Function")
    vals <- c(
      min(100, df_last$ACTH / 22 * 100),
      min(100, df_last$TSH  / 2.5 * 100),
      min(100, df_last$GH   / 2.0 * 100),
      min(100, df_last$FSH  / 5.0 * 100),
      min(100, max(0, 100 - (df_last$PRL - 15) / 0.15)),
      df_last$PFS * 100
    )
    plot_ly(type = "scatterpolar", r = vals, theta = cats, fill = "toself",
            fillcolor = "rgba(88, 166, 255, 0.3)",
            line = list(color = "#58a6ff")) %>%
      layout(polar = list(radialaxis = list(visible = TRUE, range = c(0, 120)),
                           bgcolor = "#161b22"),
             paper_bgcolor = "#0d1117",
             font = list(color = "#c9d1d9"),
             title = "Endocrine Axis Function Radar (% of Normal)")
  })

  output$table_endpoints <- renderDT({
    df <- sim_result()
    checkpoints <- c(0, 30, 90, 180, 365)
    checkpoints <- checkpoints[checkpoints <= max(df$time_d)]
    df %>%
      filter(time_d %in% checkpoints) %>%
      group_by(`Day` = time_d) %>%
      summarise(
        `Pit. Function (%)` = round(mean(PFS * 100), 1),
        `ACTH (pg/mL)`       = round(mean(ACTH), 1),
        `Cortisol (nmol/L)`  = round(mean(Cortisol), 0),
        `TSH (mIU/L)`        = round(mean(TSH), 2),
        `fT4 (pmol/L)`       = round(mean(fT4), 1),
        `GH (ng/mL)`         = round(mean(GH), 2),
        `IGF-1 (ng/mL)`      = round(mean(IGF1), 0),
        `PRL (ng/mL)`        = round(mean(PRL), 1),
        `ADH (pg/mL)`        = round(mean(ADH), 2),
        `APA (AU)`           = round(mean(APA), 1),
        .groups = "drop"
      ) %>%
      datatable(options = list(pageLength = 8, dom = "t"),
                class = "table-dark table-bordered table-sm")
  })

  ## ---- Tab 7: Biomarker ----
  output$plot_biomarker <- renderPlotly({
    df <- sim_result()
    bm <- input$biomarker_select
    if (!bm %in% names(df)) return(NULL)

    p <- ggplot(df, aes_string(x = "time_d", y = bm)) +
      geom_line(color = "#58a6ff", linewidth = 1.2)

    if (input$show_ref && bm %in% names(ref_ranges)) {
      rr <- ref_ranges[[bm]]
      p <- p +
        geom_ribbon(aes(ymin = rr[1], ymax = rr[2]),
                    fill = "#3fb950", alpha = input$ref_alpha, inherit.aes = FALSE,
                    data = df)
    }

    p <- p +
      labs(title = paste("Biomarker:", bm), x = "Time (days)", y = bm) +
      theme_dark_qsp()

    make_plotly(p)
  })

  output$plot_prl_adh <- renderPlotly({
    df <- sim_result() %>% select(time_d, PRL, ADH) %>%
      pivot_longer(-time_d, names_to = "hormone", values_to = "level") %>%
      mutate(hormone = recode(hormone, PRL = "Prolactin (ng/mL)", ADH = "ADH/Vasopressin (pg/mL)"))
    p <- ggplot(df, aes(x = time_d, y = level, color = hormone)) +
      geom_line(linewidth = 1.2) +
      facet_wrap(~hormone, scales = "free_y") +
      scale_color_manual(values = c("Prolactin (ng/mL)" = "#d2a8ff",
                                     "ADH/Vasopressin (pg/mL)" = "#79c0ff")) +
      labs(title = "Prolactin (Stalk Effect) & ADH (Posterior Involvement)",
           x = "Time (days)", y = "Level") +
      theme_dark_qsp()
    make_plotly(p)
  })

  ## ---- Download handler ----
  output$download_data <- downloadHandler(
    filename = function() paste0("lhyp_sim_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(sim_result(), file, row.names = FALSE)
  )
}

## ---- Launch ---------------------------------------------------------------
shinyApp(ui = ui, server = server)
