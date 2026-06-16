## ============================================================
## NASH QSP Shiny Application
## 비알코올성 지방간염 QSP 시뮬레이션 앱
## Version: 1.0  |  Date: 2026-06-16
##
## Dependencies: shiny, mrgsolve, dplyr, ggplot2, plotly, DT,
##               shinydashboard, shinyWidgets, patchwork, tidyr
## ============================================================

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)

## ---- Source the model compilation script (or inline here) -----
## source("nash_mrgsolve_model.R")  # provides nash_mod object

## ---- Inline model code (same as nash_mrgsolve_model.R $PROB block) ----
model_code <- '
$PARAM @annotated
DOSE_FXR  : 10  : mg
ka_FXR    : 1.2 : 1/h
CL_FXR    : 8.5 : L/h
Vc_FXR    : 30  : L
DOSE_PPAR : 800 : mg
ka_PPAR   : 0.9 : 1/h
CL_PPAR   : 12  : L/h
Vc_PPAR   : 45  : L
EC50_FXR  : 0.5 : mg/L
Emax_FXR  : 1.0 : dml
EC50_PPAR : 2.0 : mg/L
Emax_PPAR : 1.0 : dml
kFFA_in   : 0.15 : 1/h
kFFA_ox   : 0.08 : 1/h
kDNL      : 0.04 : 1/h
kTG_ester : 0.25 : 1/h
kTG_export: 0.06 : 1/h
kTG_lipo  : 0.03 : 1/h
TG_base   : 2.5  : dml
FFA_base  : 1.0  : dml
IR_base   : 0.5  : dml
kIR_DNL   : 1.8  : dml
kIR_FAO   : 0.5  : dml
k_NFkB    : 0.12 : 1/h
k_NFkB_deg: 0.18 : 1/h
TNFa_base : 1.0  : dml
kTNFa_pr  : 0.25 : 1/h
kTNFa_deg : 0.20 : 1/h
IL1b_base : 1.0  : dml
kIL1b_pr  : 0.20 : 1/h
kIL1b_deg : 0.22 : 1/h
kROS_gen  : 0.10 : 1/h
kROS_scav : 0.15 : 1/h
ROS_base  : 1.0  : dml
kASK1_act : 0.08 : 1/h
kASK1_deg : 0.12 : 1/h
kJNK_act  : 0.15 : 1/h
kJNK_deg  : 0.18 : 1/h
kHSC_act  : 0.05 : 1/h
kHSC_res  : 0.02 : 1/h
kTGFb_pr  : 0.08 : 1/h
kTGFb_deg : 0.10 : 1/h
kCOL_pr   : 0.04 : 1/h
kCOL_deg  : 0.005 : 1/h
COL_base  : 1.0   : dml
ALT_base  : 30    : U/L
kALT_inj  : 0.80  : dml
kALT_cl   : 0.04  : 1/h
Emax_FAO_PPAR : 1.5 : dml
EC50_FAO_PPAR : 1.2 : mg/L
Emax_FXR_antiinf : 0.70 : dml
EC50_FXR_antiinf : 0.4  : mg/L

$CMT @annotated
GUT_FXR  : Gut FXR agonist (mg)
CENT_FXR : Plasma FXR agonist (mg)
GUT_PPAR : Gut PPAR agonist (mg)
CENT_PPAR: Plasma PPAR agonist (mg)
FFA_L    : Hepatic FFA
TG_L     : Hepatic TG
NFkB_A   : Active NF-kB
TNFa     : TNF-alpha
IL1b     : IL-1 beta
ROS      : Reactive oxygen species
ASK1_A   : Activated ASK1
JNK_A    : Activated JNK
HSC_A    : Activated HSC (fraction)
TGFb1    : TGF-beta1
COL      : Hepatic collagen
ALT_p    : Plasma ALT

$MAIN
double Cp_FXR  = CENT_FXR  / Vc_FXR;
double Cp_PPAR = CENT_PPAR / Vc_PPAR;
double occ_FXR  = Emax_FXR  * Cp_FXR  / (EC50_FXR  + Cp_FXR);
double occ_PPAR = Emax_PPAR * Cp_PPAR / (EC50_PPAR + Cp_PPAR);
double lip_tox  = (FFA_L / FFA_base) * (TG_L / TG_base);
double FAO_fold = 1.0 + Emax_FAO_PPAR * Cp_PPAR / (EC50_FAO_PPAR + Cp_PPAR);
double FAO_rate = kFFA_ox * FAO_fold;
double FXR_DNL_supp  = 1.0 - 0.60 * occ_FXR;
double DNL_rate = kDNL * (1.0 + kIR_DNL * IR_base) * FXR_DNL_supp;
double FXR_NFkB_supp = 1.0 - Emax_FXR_antiinf * Cp_FXR / (EC50_FXR_antiinf + Cp_FXR);
double NFkB_prod = k_NFkB * lip_tox * FXR_NFkB_supp;
double HSC_drive = TGFb1 * (1.0 + 0.5 * TNFa) * (1.0 + 0.3 * JNK_A);
double PPAR_antifib = 1.0 - 0.50 * occ_PPAR;
double COL_prod = kCOL_pr * HSC_A * PPAR_antifib;
double inj_signal = lip_tox + JNK_A + 0.5 * TNFa;

$ODE
dxdt_GUT_FXR   = -ka_FXR * GUT_FXR;
dxdt_CENT_FXR  =  ka_FXR * GUT_FXR - (CL_FXR/Vc_FXR)*CENT_FXR;
dxdt_GUT_PPAR  = -ka_PPAR * GUT_PPAR;
dxdt_CENT_PPAR =  ka_PPAR * GUT_PPAR - (CL_PPAR/Vc_PPAR)*CENT_PPAR;
dxdt_FFA_L  = kFFA_in*FFA_base + DNL_rate - FAO_rate*FFA_L - kTG_ester*FFA_L;
dxdt_TG_L   = kTG_ester*FFA_L - kTG_export*TG_L - kTG_lipo*TG_L;
dxdt_NFkB_A = NFkB_prod - k_NFkB_deg*NFkB_A;
dxdt_TNFa   = kTNFa_pr*NFkB_A - kTNFa_deg*TNFa;
dxdt_IL1b   = kIL1b_pr*NFkB_A*(1.0+0.4*lip_tox) - kIL1b_deg*IL1b;
double Nrf2_fold = 1.0 + 0.3*occ_PPAR;
dxdt_ROS    = kROS_gen*lip_tox - kROS_scav*Nrf2_fold*ROS;
dxdt_ASK1_A = kASK1_act*ROS - kASK1_deg*ASK1_A;
dxdt_JNK_A  = kJNK_act*ASK1_A - kJNK_deg*JNK_A;
dxdt_HSC_A  = kHSC_act*HSC_drive*(1.0-HSC_A) - kHSC_res*HSC_A;
dxdt_TGFb1  = kTGFb_pr*HSC_A - kTGFb_deg*TGFb1;
dxdt_COL    = COL_prod - kCOL_deg*COL;
dxdt_ALT_p  = kALT_inj*inj_signal - kALT_cl*ALT_p;

$TABLE
double Steatosis_pct = 100.0*(TG_L/TG_base)/(1.0+TG_L/TG_base)*0.66;
double NAS_steat  = 3.0*Steatosis_pct/100.0;
double NAS_inflam = 3.0*(TNFa+IL1b)/2.0/(1.0+(TNFa+IL1b)/2.0);
double NAS_bal    = 2.0*JNK_A/(1.0+JNK_A);
double NAS_total  = NAS_steat + NAS_inflam + NAS_bal;
double Fibrosis_stage = 4.0*(COL-1.0)/(1.0+(COL-1.0));
double Fibrosis_stage2 = (Fibrosis_stage < 0.0) ? 0.0 : Fibrosis_stage;
double FIB4_proxy = ALT_p/30.0*(COL/COL_base);
double Cp_FXR_out  = CENT_FXR/Vc_FXR;
double Cp_PPAR_out = CENT_PPAR/Vc_PPAR;
double Eff_FXR  = Emax_FXR*Cp_FXR_out/(EC50_FXR+Cp_FXR_out);
double Eff_PPAR = Emax_PPAR*Cp_PPAR_out/(EC50_PPAR+Cp_PPAR_out);

$CAPTURE
Cp_FXR_out Cp_PPAR_out Eff_FXR Eff_PPAR
Steatosis_pct NAS_total NAS_steat NAS_inflam NAS_bal
Fibrosis_stage2 FIB4_proxy
TG_L FFA_L NFkB_A TNFa IL1b ROS ASK1_A JNK_A HSC_A TGFb1 COL ALT_p
'

nash_mod <- mcode("nash_qsp_shiny", model_code)

## ============================================================
## HELPER FUNCTION: Run a single scenario
## ============================================================
run_scenario <- function(
  dose_fxr, dose_ppar, ir_level, duration_wk,
  init_tg = 2.5, init_col = 1.5, init_hsc = 0.25, init_alt = 60
) {
  mod <- param(nash_mod,
               IR_base   = ir_level,
               DOSE_FXR  = dose_fxr,
               DOSE_PPAR = dose_ppar)

  init_v <- c(
    GUT_FXR = 0, CENT_FXR = 0, GUT_PPAR = 0, CENT_PPAR = 0,
    FFA_L = 1.0, TG_L = init_tg,
    NFkB_A = 1.2, TNFa = 1.0, IL1b = 1.0,
    ROS = 1.5, ASK1_A = 0.7, JNK_A = 0.6,
    HSC_A = init_hsc, TGFb1 = 0.3, COL = init_col,
    ALT_p = init_alt
  )

  n_doses <- duration_wk * 7
  ev_list <- list()

  if (dose_fxr > 0)
    ev_list <- c(ev_list, list(ev(amt = dose_fxr, cmt = 1, ii = 24, addl = n_doses - 1)))
  if (dose_ppar > 0)
    ev_list <- c(ev_list, list(ev(amt = dose_ppar, cmt = 3, ii = 24, addl = n_doses - 1)))

  if (length(ev_list) == 0)
    ev_final <- ev(time = 0, amt = 0, cmt = 1)
  else
    ev_final <- Reduce(c, ev_list)

  mod %>%
    init(init_v) %>%
    mrgsim(events = ev_final, end = n_doses * 24, delta = 12) %>%
    as_tibble() %>%
    mutate(time_wk = time / (7 * 24))
}

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = "NASH QSP Simulator",
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      id = "tabs",
      menuItem("Overview",          tabName = "overview",    icon = icon("info-circle")),
      menuItem("PK/PD Simulation",  tabName = "pkpd",        icon = icon("chart-line")),
      menuItem("Disease Endpoints", tabName = "endpoints",   icon = icon("stethoscope")),
      menuItem("Mechanism Detail",  tabName = "mechanism",   icon = icon("project-diagram")),
      menuItem("Data Table",        tabName = "datatable",   icon = icon("table")),
      menuItem("References",        tabName = "references",  icon = icon("book-open"))
    ),

    hr(),
    ## ---- GLOBAL INPUTS ----
    h5("  Patient/Disease Settings", style = "color:#90CAF9; margin-left:15px;"),

    sliderInput("ir_level",
                "Insulin Resistance Level",
                min = 0, max = 1, value = 0.7, step = 0.05),

    sliderInput("init_tg",
                "Initial Hepatic TG (fold of normal)",
                min = 1.0, max = 5.0, value = 2.5, step = 0.1),

    sliderInput("init_col",
                "Initial Fibrosis (Collagen Index)",
                min = 1.0, max = 3.0, value = 1.5, step = 0.1),

    sliderInput("duration_wk",
                "Simulation Duration (weeks)",
                min = 12, max = 104, value = 52, step = 4),

    hr(),
    h5("  Treatment Arms", style = "color:#A5D6A7; margin-left:15px;"),

    prettyCheckbox("use_notreat", "No Treatment",            TRUE,  icon = icon("check")),
    prettyCheckbox("use_fxr",    "FXR Agonist (OCA-like)",  FALSE, icon = icon("check")),
    prettyCheckbox("use_ppar",   "pan-PPAR (Lanifibranor)", FALSE, icon = icon("check")),
    prettyCheckbox("use_combo",  "Combination",             FALSE, icon = icon("check")),

    hr(),
    h5("  Dose Settings", style = "color:#F8BBD0; margin-left:15px;"),

    numericInput("dose_fxr",  "FXR Agonist dose (mg/day)", value = 10,  min = 0, max = 50),
    numericInput("dose_ppar", "pan-PPAR dose (mg/day)",    value = 800, min = 0, max = 1200),

    br(),
    actionBttn("run_sim", "Run Simulation",
               style = "gradient", color = "primary",
               icon = icon("play"), block = TRUE)
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .box { border-top-color: #1565C0; }
      .content-wrapper { background-color: #F5F7FA; }
      .ggplot-container { background: white; border-radius: 8px; padding: 10px; }
    "))),

    tabItems(

      ## ============================================================
      ## TAB 1: OVERVIEW
      ## ============================================================
      tabItem("overview",
        fluidRow(
          box(
            title = "NASH / NAFLD — QSP Model Overview",
            status = "primary", solidHeader = TRUE, width = 12,
            fluidRow(
              column(6,
                h4("Disease: 비알코올성 지방간염 (NASH)"),
                p("Non-Alcoholic Steatohepatitis (NASH) is a progressive liver disease
                  characterized by hepatic steatosis, lobular inflammation, hepatocyte
                  ballooning, and varying degrees of fibrosis. It represents the severe
                  end of the NAFLD spectrum and can progress to cirrhosis and HCC."),
                h5("Key Pathophysiological Modules:"),
                tags$ul(
                  tags$li("Adipose ↔ Liver lipid flux (FFA, TG, DNL, FAO)"),
                  tags$li("Insulin Resistance (IRS1/Akt/FOXO1/SREBP-1c axis)"),
                  tags$li("Oxidative Stress (ROS → ASK1 → JNK → Apoptosis)"),
                  tags$li("Innate Immunity (TLR4 → NF-κB → TNF-α, IL-1β, NLRP3)"),
                  tags$li("Fibrosis (TGF-β1 → HSC activation → Collagen deposition)"),
                  tags$li("Gut–Liver Axis (dysbiosis → LPS → Kupffer cells)"),
                  tags$li("Nuclear Receptors (FXR, PPARα/γ/δ, THRβ, LXR)")
                )
              ),
              column(6,
                h5("Modeled Drug Mechanisms:"),
                tags$table(
                  class = "table table-striped table-bordered",
                  tags$thead(tags$tr(
                    tags$th("Drug Class"), tags$th("Example"), tags$th("Primary Target"), tags$th("Phase")
                  )),
                  tags$tbody(
                    tags$tr(tags$td("FXR Agonist"),       tags$td("OCA, Cilofexor"), tags$td("FXR → ↓DNL, ↓Inflam"), tags$td("Phase 3 (OCA)")),
                    tags$tr(tags$td("pan-PPAR Agonist"),  tags$td("Lanifibranor"),   tags$td("PPARα/γ/δ → ↑FAO, ↓Fibr"), tags$td("Phase 3")),
                    tags$tr(tags$td("ACC Inhibitor"),     tags$td("Firsocostat"),    tags$td("ACC1/2 → ↓DNL"),      tags$td("Phase 2b")),
                    tags$tr(tags$td("ASK1 Inhibitor"),    tags$td("Selonsertib"),    tags$td("ASK1 → ↓JNK, apopt"), tags$td("Phase 3 (failed)")),
                    tags$tr(tags$td("CCR2/5 Antagonist"), tags$td("Cenicriviroc"),   tags$td("CCR2/5 → ↓Kupffer"), tags$td("Phase 3")),
                    tags$tr(tags$td("GLP-1 RA"),          tags$td("Semaglutide"),    tags$td("↓Weight, IR, inflam"), tags$td("Phase 3")),
                    tags$tr(tags$td("THRβ Agonist"),      tags$td("Resmetirom"),     tags$td("THRβ → ↑FAO, ↓TG"),  tags$td("FDA Approved 2024"))
                  )
                )
              )
            )
          )
        ),
        fluidRow(
          valueBoxOutput("vb_nodes",  width = 3),
          valueBoxOutput("vb_edges",  width = 3),
          valueBoxOutput("vb_drugs",  width = 3),
          valueBoxOutput("vb_refs",   width = 3)
        )
      ),

      ## ============================================================
      ## TAB 2: PK/PD SIMULATION
      ## ============================================================
      tabItem("pkpd",
        fluidRow(
          box(
            title = "Drug Concentration–Time & PD Occupancy", status = "primary",
            solidHeader = TRUE, width = 12,
            fluidRow(
              column(6, plotlyOutput("plot_pk_fxr",  height = "320px")),
              column(6, plotlyOutput("plot_pk_ppar", height = "320px"))
            ),
            fluidRow(
              column(6, plotlyOutput("plot_eff_fxr",  height = "280px")),
              column(6, plotlyOutput("plot_eff_ppar", height = "280px"))
            )
          )
        )
      ),

      ## ============================================================
      ## TAB 3: DISEASE ENDPOINTS
      ## ============================================================
      tabItem("endpoints",
        fluidRow(
          box(
            title = "NAS Score & Fibrosis Stage",
            status = "primary", solidHeader = TRUE, width = 12,
            fluidRow(
              column(6, plotlyOutput("plot_nas",      height = "350px")),
              column(6, plotlyOutput("plot_fibrosis", height = "350px"))
            )
          )
        ),
        fluidRow(
          box(
            title = "Hepatic TG / Steatosis & ALT",
            status = "warning", solidHeader = TRUE, width = 12,
            fluidRow(
              column(6, plotlyOutput("plot_tg",  height = "320px")),
              column(6, plotlyOutput("plot_alt", height = "320px"))
            )
          )
        )
      ),

      ## ============================================================
      ## TAB 4: MECHANISM DETAIL
      ## ============================================================
      tabItem("mechanism",
        fluidRow(
          box(
            title = "Inflammation Biomarkers (TNF-α, IL-1β, NF-κB)",
            status = "danger", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_inflam", height = "380px")
          ),
          box(
            title = "Oxidative Stress Cascade (ROS → ASK1 → JNK)",
            status = "danger", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_oxstress", height = "380px")
          )
        ),
        fluidRow(
          box(
            title = "Fibrosis Mediators (TGF-β1, HSC, Collagen)",
            status = "purple", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_fibmech", height = "380px")
          ),
          box(
            title = "Spider/Radar Chart — 52-week Biomarker Summary",
            status = "info", solidHeader = TRUE, width = 6,
            plotlyOutput("plot_radar", height = "380px")
          )
        )
      ),

      ## ============================================================
      ## TAB 5: DATA TABLE
      ## ============================================================
      tabItem("datatable",
        fluidRow(
          box(
            title = "Simulation Results (Weekly Summary)", status = "primary",
            solidHeader = TRUE, width = 12,
            selectInput("arm_select",
                        "Select Treatment Arm:",
                        choices  = c("No Treatment", "FXR Agonist (OCA-like)",
                                     "pan-PPAR (Lanifibranor)", "Combination"),
                        selected = "No Treatment"),
            downloadButton("dl_csv", "Download CSV"),
            br(), br(),
            DTOutput("results_table")
          )
        )
      ),

      ## ============================================================
      ## TAB 6: REFERENCES
      ## ============================================================
      tabItem("references",
        fluidRow(
          box(
            title = "Key Literature References", status = "info",
            solidHeader = TRUE, width = 12,
            DT::datatable(
              data.frame(
                `#` = 1:10,
                Reference = c(
                  "Diehl AM & Day C (2017). N Engl J Med 377:2063",
                  "Friedman SL et al. (2018). Nat Med 24:908",
                  "Neuschwander-Tetri BA et al. (2015). Lancet 385:956 [OCA FLINT]",
                  "Ratziu V et al. (2021). N Engl J Med 385:1547 [Lanifibranor NATIVE]",
                  "Harrison SA et al. (2023). N Engl J Med 388:2545 [Resmetirom MAESTRO]",
                  "Vasilyeva O et al. (2021). CPT Pharmacometrics Syst Pharmacol 10:1369",
                  "Bataller R & Brenner DA (2005). J Clin Invest 115:209",
                  "Kleiner DE et al. (2005). Hepatology 41:1313 [NAS scoring]",
                  "Loomba R et al. (2021). Cell Metab 34:531 [Firsocostat ATLAS]",
                  "Newsome PN et al. (2021). N Engl J Med 384:1113 [Semaglutide]"
                ),
                PubMed_ID = c("29166236", "29967350", "25468160", "34670043",
                              "37059523", "34714976", "15690074", "15461435",
                              "34270930", "33567104"),
                stringsAsFactors = FALSE
              ),
              options = list(pageLength = 10, dom = "Bfrtip"),
              rownames = FALSE
            )
          )
        )
      )
    )  # end tabItems
  )   # end dashboardBody
)

## ============================================================
## SERVER
## ============================================================
server <- function(input, output, session) {

  ## ---- Value boxes for overview ----
  output$vb_nodes  <- renderValueBox(valueBox(130, "Mechanistic Nodes",  icon = icon("circle"), color = "blue"))
  output$vb_edges  <- renderValueBox(valueBox(160, "Pathway Edges",      icon = icon("arrows-alt-h"), color = "green"))
  output$vb_drugs  <- renderValueBox(valueBox(10,  "Drug Mechanisms",    icon = icon("pills"), color = "orange"))
  output$vb_refs   <- renderValueBox(valueBox(60,  "Literature Refs",    icon = icon("book"), color = "purple"))

  ## ---- Reactive: collect simulations ----
  sims_all <- eventReactive(input$run_sim, {
    withProgress(message = "Running QSP simulations...", value = 0, {
      arms <- list()

      if (input$use_notreat) {
        incProgress(0.25, detail = "No Treatment")
        arms[["No Treatment"]] <- run_scenario(
          0, 0, input$ir_level, input$duration_wk,
          input$init_tg, input$init_col
        )
      }
      if (input$use_fxr) {
        incProgress(0.25, detail = "FXR Agonist")
        arms[["FXR Agonist (OCA-like)"]] <- run_scenario(
          input$dose_fxr, 0, input$ir_level, input$duration_wk,
          input$init_tg, input$init_col
        )
      }
      if (input$use_ppar) {
        incProgress(0.25, detail = "pan-PPAR")
        arms[["pan-PPAR (Lanifibranor)"]] <- run_scenario(
          0, input$dose_ppar, input$ir_level, input$duration_wk,
          input$init_tg, input$init_col
        )
      }
      if (input$use_combo) {
        incProgress(0.25, detail = "Combination")
        arms[["Combination"]] <- run_scenario(
          input$dose_fxr, input$dose_ppar, input$ir_level, input$duration_wk,
          input$init_tg, input$init_col
        )
      }

      if (length(arms) == 0) return(NULL)

      bind_rows(lapply(names(arms), function(n) arms[[n]] %>% mutate(arm = n)))
    })
  })

  arm_palette <- c(
    "No Treatment"           = "#B71C1C",
    "FXR Agonist (OCA-like)" = "#1565C0",
    "pan-PPAR (Lanifibranor)" = "#2E7D32",
    "Combination"            = "#6A1B9A"
  )

  make_plotly <- function(df, y_var, title, ylab) {
    req(df)
    p <- ggplot(df, aes(x = time_wk, y = .data[[y_var]], color = arm)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = arm_palette) +
      labs(title = title, x = "Time (weeks)", y = ylab, color = NULL) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom")
    ggplotly(p) %>% layout(legend = list(orientation = "h", y = -0.25))
  }

  ## ---- PK plots ----
  output$plot_pk_fxr  <- renderPlotly({
    make_plotly(sims_all(), "Cp_FXR_out",  "FXR Agonist — Plasma Conc.", "Cp (mg/L)")
  })
  output$plot_pk_ppar <- renderPlotly({
    make_plotly(sims_all(), "Cp_PPAR_out", "pan-PPAR Agonist — Plasma Conc.", "Cp (mg/L)")
  })
  output$plot_eff_fxr  <- renderPlotly({
    make_plotly(sims_all(), "Eff_FXR",  "FXR Receptor Occupancy", "Fractional Occupancy (0–1)")
  })
  output$plot_eff_ppar <- renderPlotly({
    make_plotly(sims_all(), "Eff_PPAR", "pan-PPAR Receptor Occupancy", "Fractional Occupancy (0–1)")
  })

  ## ---- Endpoint plots ----
  output$plot_nas      <- renderPlotly({ make_plotly(sims_all(), "NAS_total",      "NAS Score (0–8)",                "NAS") })
  output$plot_fibrosis <- renderPlotly({ make_plotly(sims_all(), "Fibrosis_stage2","Fibrosis Stage (Metavir proxy)",  "Stage") })
  output$plot_tg       <- renderPlotly({ make_plotly(sims_all(), "TG_L",           "Hepatic TG (normalized)",         "TG_L") })
  output$plot_alt      <- renderPlotly({ make_plotly(sims_all(), "ALT_p",          "Plasma ALT (U/L)",                "ALT (U/L)") })

  ## ---- Mechanism plots ----
  output$plot_inflam <- renderPlotly({
    df <- sims_all()
    req(df)
    df_long <- df %>%
      select(time_wk, arm, TNFa, IL1b, NFkB_A) %>%
      pivot_longer(c(TNFa, IL1b, NFkB_A), names_to = "marker", values_to = "value")
    p <- ggplot(df_long, aes(time_wk, value, color = arm, linetype = marker)) +
      geom_line(linewidth = 0.9) +
      scale_color_manual(values = arm_palette) +
      labs(title = "Inflammation Markers", x = "Time (wk)", y = "Relative Units") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$plot_oxstress <- renderPlotly({
    df <- sims_all()
    req(df)
    df_long <- df %>%
      select(time_wk, arm, ROS, ASK1_A, JNK_A) %>%
      pivot_longer(c(ROS, ASK1_A, JNK_A), names_to = "marker", values_to = "value")
    p <- ggplot(df_long, aes(time_wk, value, color = arm, linetype = marker)) +
      geom_line(linewidth = 0.9) +
      scale_color_manual(values = arm_palette) +
      labs(title = "Oxidative Stress Cascade", x = "Time (wk)", y = "Relative Units") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$plot_fibmech <- renderPlotly({
    df <- sims_all()
    req(df)
    df_long <- df %>%
      select(time_wk, arm, TGFb1, HSC_A, COL) %>%
      pivot_longer(c(TGFb1, HSC_A, COL), names_to = "marker", values_to = "value")
    p <- ggplot(df_long, aes(time_wk, value, color = arm, linetype = marker)) +
      geom_line(linewidth = 0.9) +
      scale_color_manual(values = arm_palette) +
      labs(title = "Fibrosis Mediators", x = "Time (wk)", y = "Relative Units") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
    ggplotly(p)
  })

  ## ---- Radar chart at 52 wk ----
  output$plot_radar <- renderPlotly({
    df <- sims_all()
    req(df)
    end_vals <- df %>%
      group_by(arm) %>%
      filter(time_wk == max(time_wk)) %>%
      slice_tail(n = 1) %>%
      ungroup() %>%
      select(arm, TG_L, TNFa, ROS, JNK_A, COL, ALT_p) %>%
      mutate(across(where(is.numeric), ~ . / max(., na.rm = TRUE)))

    fig <- plot_ly(type = "scatterpolar", fill = "toself")
    for (a in unique(end_vals$arm)) {
      row  <- end_vals %>% filter(arm == a)
      vals <- unlist(row[, -1])
      fig  <- fig %>% add_trace(
        r    = c(vals, vals[1]),
        theta = c(names(vals), names(vals)[1]),
        name  = a
      )
    }
    fig %>% layout(
      polar = list(radialaxis = list(visible = TRUE, range = c(0, 1))),
      title = "Normalized Biomarker Radar (at end of simulation)"
    )
  })

  ## ---- Data table ----
  output$results_table <- renderDT({
    df <- sims_all()
    req(df)
    df %>%
      filter(arm == input$arm_select,
             time_wk %% 4 < 0.5) %>%   # weekly summary
      select(time_wk, TG_L, Steatosis_pct, NAS_total, Fibrosis_stage2,
             TNFa, IL1b, ROS, JNK_A, HSC_A, COL, ALT_p) %>%
      mutate(across(where(is.numeric), ~ round(., 3))) %>%
      datatable(
        colnames = c("Week", "Hepatic TG", "Steatosis %", "NAS",
                     "Fibrosis Stage", "TNF-α", "IL-1β", "ROS",
                     "JNK", "HSC Act.", "Collagen", "ALT (U/L)"),
        options = list(pageLength = 15, scrollX = TRUE)
      )
  })

  ## ---- Download ----
  output$dl_csv <- downloadHandler(
    filename = function() paste0("nash_qsp_", input$arm_select, "_", Sys.Date(), ".csv"),
    content  = function(file) {
      df <- sims_all()
      req(df)
      write.csv(df %>% filter(arm == input$arm_select), file, row.names = FALSE)
    }
  )

}

## ============================================================
## LAUNCH
## ============================================================
shinyApp(ui, server)
