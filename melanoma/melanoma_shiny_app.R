## ============================================================
## Melanoma QSP Shiny App
## Interactive PK/PD dashboard: BRAF-mutant cutaneous melanoma
## Tabs: Patient Profile | Drug PK | Tumor Response |
##       Immune Dynamics | Scenario Comparison | Biomarkers
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)

## ============================================================
## Inline mrgsolve Model
## ============================================================
MEL_CODE <- '
$PARAM
ka_B=0.42, F_B=0.64, Vd_B=100.0, CL_B=13.0, DOSE_B=0.0,
ka_M=0.35, F_M=0.72, Vd_M=214.0, CL_M=4.9, DOSE_M=0.0,
CL_I=0.244, Vc_I=3.61, Vp_I=2.75, Q_I=0.448, DOSE_I=0.0, DOSE_C4=0.0,
kout_ERK=0.12,
EC50_BRAFi=0.28, Emax_BRAFi=0.90, hill_B=1.8,
EC50_MEKi=12.0,  Emax_MEKi=0.85, hill_M=1.5,
kg0=0.0055, kmax_tum=10.0, kd_ERK=0.0040, kd_immune=0.0050,
BRAF_V600=1.0,
kR_on=0.00005, Rmax=1.0,
kin_CD8=0.018, kout_CD8=0.010,
CD8_ICI_fac=1.8, CD8_CTLA4=1.4, Treg_sup=0.60,
EC50_PD1=12.0, EC50_CTLA4=20.0,
Treg0=0.15, Treg_ICI=0.40,
kin_IFNg=0.002, kout_IFNg=0.05, IFNg_CD8=0.02, kPDL1_IFNg=0.30,
LDH0=200.0, LDH_tum_fac=1.8, kout_LDH=0.03,
S100B_0=0.10, S100B_tum=2.0,
TMB_val=10.0, PDL1_pct=30.0

$CMT
GUT_BRAF CENT_BRAF GUT_MEK CENT_MEK CENT_ICI PERI_ICI
ERK_act RESIST TUMOR CD8_TIL PD1_RO CTLA4_RO
Treg_frac IFNg_TME LDH_ser S100B_ser

$MAIN
double Cb = CENT_BRAF / Vd_B;
double Cm = CENT_MEK  / Vd_M;
double Ci = CENT_ICI  / Vc_I;
double inh_B = BRAF_V600 * Emax_BRAFi * pow(Cb, hill_B) /
               (pow(EC50_BRAFi, hill_B) + pow(Cb, hill_B));
double inh_M = Emax_MEKi * pow(Cm, hill_M) /
               (pow(EC50_MEKi, hill_M) + pow(Cm, hill_M));
double ERK_inh = inh_B + inh_M - inh_B * inh_M;
double ERK_target = 1.0 - ERK_inh;
ERK_target = ERK_target + RESIST * (1.0 - ERK_target);
if (ERK_target < 0.0) ERK_target = 0.0;
if (ERK_target > 1.0) ERK_target = 1.0;
double PD1_eq   = Ci / (Ci + EC50_PD1);
double CTA4_eq  = Ci / (Ci + EC50_CTLA4);
double PDL1_ind = 1.0 + kPDL1_IFNg * IFNg_TME;
double PD1_sup  = 1.0 - PD1_RO;
double TMB_fac  = 1.0 + (TMB_val - 10.0) / 50.0;
double PDL1_fac = 1.0 + PDL1_pct / 100.0;
double CD8_eff  = CD8_TIL * TMB_fac * PDL1_fac *
                  (1.0 - Treg_frac * Treg_sup) *
                  (1.0 - 0.70 * PD1_sup / PDL1_ind);
if (CD8_eff < 0.0) CD8_eff = 0.0;
double kill_ERK    = kd_ERK    * (1.0 - ERK_act) * BRAF_V600;
double kill_immune = kd_immune * CD8_eff;
TUMOR_0    = 1.0; CD8_TIL_0 = 1.0; ERK_act_0 = 1.0;
RESIST_0   = 0.0; PD1_RO_0  = 0.0; CTLA4_RO_0 = 0.0;
Treg_frac_0 = Treg0; IFNg_TME_0 = 1.0;
LDH_ser_0  = LDH0; S100B_ser_0 = S100B_0;

$ODE
dxdt_GUT_BRAF  = -ka_B * GUT_BRAF;
dxdt_CENT_BRAF = F_B * ka_B * GUT_BRAF - (CL_B / Vd_B) * CENT_BRAF;
dxdt_GUT_MEK   = -ka_M * GUT_MEK;
dxdt_CENT_MEK  = F_M * ka_M * GUT_MEK - (CL_M / Vd_M) * CENT_MEK;
double Ci_c = CENT_ICI / Vc_I; double Cp_c = PERI_ICI / Vp_I;
dxdt_CENT_ICI  = -(CL_I / Vc_I) * CENT_ICI - Q_I * (Ci_c - Cp_c);
dxdt_PERI_ICI  =  Q_I * (Ci_c - Cp_c);
double kin_ERK = kout_ERK * ERK_target;
dxdt_ERK_act   = kin_ERK - kout_ERK * ERK_act;
double dp = (CENT_BRAF > 0 || CENT_MEK > 0) ? 1.0 : 0.0;
dxdt_RESIST    = kR_on * dp * (Rmax - RESIST);
double tg = (TUMOR > 0.01) ? kg0 * TUMOR * log(kmax_tum / TUMOR) : 0.0;
double tk = (kill_ERK + kill_immune) * TUMOR;
dxdt_TUMOR     = tg - tk;
double cd8r = kin_CD8 * (1.0 + (CD8_ICI_fac-1.0)*PD1_RO +
                               (CD8_CTLA4-1.0)*CTLA4_RO);
dxdt_CD8_TIL   = cd8r - kout_CD8 * CD8_TIL;
dxdt_PD1_RO    = 0.5 * (PD1_eq  - PD1_RO);
dxdt_CTLA4_RO  = 0.5 * (CTA4_eq - CTLA4_RO);
double Treg_eq = Treg0 * (1.0 - Treg_ICI * CTLA4_RO);
dxdt_Treg_frac = 0.1 * (Treg_eq - Treg_frac);
double ki2 = kin_IFNg + IFNg_CD8 * CD8_TIL * PD1_RO;
dxdt_IFNg_TME  = ki2 - kout_IFNg * IFNg_TME;
double LDH_eq  = LDH0 * (1.0 + (LDH_tum_fac-1.0)*(TUMOR-1.0));
if (LDH_eq < LDH0 * 0.5) LDH_eq = LDH0 * 0.5;
dxdt_LDH_ser   = kout_LDH * (LDH_eq - LDH_ser);
double S100_eq = S100B_0 * (1.0 + (S100B_tum-1.0)*(TUMOR-1.0));
if (S100_eq < S100B_0 * 0.3) S100_eq = S100B_0 * 0.3;
dxdt_S100B_ser = 0.05 * (S100_eq - S100B_ser);

$CAPTURE
double Cb_out = CENT_BRAF / Vd_B;
double Cm_out = CENT_MEK  / Vd_M;
double Ci_out = CENT_ICI  / Vc_I;
double ERK_pct   = ERK_act * 100.0;
double Tumor_pct = TUMOR   * 100.0;
double CD8_rel   = CD8_TIL;
double Resist_pct = RESIST * 100.0;
double PD1_pct   = PD1_RO  * 100.0;
double CTLA4_pct = CTLA4_RO * 100.0;
double Treg_pct  = Treg_frac * 100.0;
double IFNg_rel  = IFNg_TME;
double LDH_val   = LDH_ser;
double S100_val  = S100B_ser;
'

mel_mod <- mcode("melanoma_shiny_qsp", MEL_CODE)

## ============================================================
## Helper: Run Simulation
## ============================================================
run_sim <- function(mod, params, dose_B, tau_B, dose_M, dose_I, q_I,
                    dose_C4, use_seq, seq_switch_wk, duration_wk) {
  mod <- param(mod, params)
  sim_end <- duration_wk * 7 * 24

  ev_list <- list()
  if (dose_B > 0) {
    n_B <- floor(sim_end / tau_B)
    if (use_seq) n_B <- floor(seq_switch_wk * 7 * 24 / tau_B)
    ev_list$braf <- ev(amt = dose_B, ii = tau_B, addl = n_B - 1, cmt = 1)
  }
  if (dose_M > 0) {
    n_M <- if (use_seq) floor(seq_switch_wk * 7) else floor(duration_wk * 7)
    ev_list$mek <- ev(amt = dose_M, ii = 24, addl = n_M - 1, cmt = 3)
  }
  if (dose_I > 0) {
    start_t <- if (use_seq && (dose_B > 0 || dose_M > 0))
                 seq_switch_wk * 7 * 24 else 0
    n_I <- floor((sim_end - start_t) / (q_I * 24))
    if (n_I < 1) n_I <- 1
    ev_list$ici <- ev(amt = dose_I, ii = q_I * 24, addl = n_I - 1,
                      cmt = 5, rate = dose_I / 0.5, time = start_t)
  }
  if (dose_C4 > 0) {
    n_C4 <- min(4, floor(sim_end / (21 * 24)))
    ev_list$ctla4 <- ev(amt = dose_C4, ii = 21 * 24, addl = n_C4 - 1,
                        cmt = 5, rate = dose_C4 / 0.5)
  }

  total_ev <- if (length(ev_list) > 0) {
    Reduce(`+`, ev_list)
  } else ev()

  mrgsim(mod, ev = total_ev, end = sim_end, delta = 4) %>%
    as.data.frame() %>%
    mutate(time_wk = time / (7 * 24))
}

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "black",
  dashboardHeader(title = "Melanoma QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "tab_patient", icon = icon("user")),
      menuItem("Drug PK",            tabName = "tab_pk",      icon = icon("pills")),
      menuItem("Tumor Response",     tabName = "tab_tumor",   icon = icon("chart-line")),
      menuItem("Immune Dynamics",    tabName = "tab_immune",  icon = icon("shield-alt")),
      menuItem("Scenario Comparison",tabName = "tab_scen",    icon = icon("bar-chart")),
      menuItem("Biomarkers",         tabName = "tab_bio",     icon = icon("flask"))
    ),
    hr(),
    actionButton("run_sim", "Run Simulation",
                 icon = icon("play"),
                 style = "color:#fff; background-color:#e74c3c; margin:10px;
                          width:180px; font-weight:bold;")
  ),

  dashboardBody(
    tabItems(

      ## ---- Tab 1: Patient Profile ----
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient & Tumor Genetics", status = "danger",
              solidHeader = TRUE, width = 4,
              selectInput("braf_status", "BRAF Mutation Status:",
                          c("V600E (most common, ~50%)" = "V600E",
                            "V600K (~5%, more UV-assoc.)" = "V600K",
                            "NRAS Q61 mutation (~20%)" = "NRAS",
                            "NF1 loss (~15%)" = "NF1",
                            "WT / other (<10%)" = "WT")),
              selectInput("stage", "Stage:",
                          c("Stage IIIB (regional)" = "3B",
                            "Stage IIIC/D" = "3C",
                            "Stage IVA (M1a)" = "4A",
                            "Stage IVB (M1b - lung)" = "4B",
                            "Stage IVC (M1c/d - visceral/CNS)" = "4C")),
              numericInput("ecog", "ECOG Performance Status:", 0, 0, 2, 1),
              numericInput("ldh_base", "Baseline LDH (U/L):", 200, 80, 1500, 10)
          ),
          box(title = "Biomarker & Treatment History", status = "warning",
              solidHeader = TRUE, width = 4,
              numericInput("tmb_val", "Tumor Mutational Burden (mut/Mb):", 10, 1, 100, 1),
              numericInput("pdl1_pct", "PD-L1 Expression (% TPS):", 30, 0, 100, 5),
              selectInput("prior_tx", "Prior Therapy:",
                          c("Treatment-naive" = "naive",
                            "Prior BRAFi/MEKi" = "prior_braf",
                            "Prior anti-PD-1" = "prior_pd1",
                            "Prior ipilimumab" = "prior_ipi")),
              numericInput("duration_wk", "Simulation Duration (weeks):", 52, 12, 104, 4)
          ),
          box(title = "LDH Prognostic Assessment", status = "info",
              solidHeader = TRUE, width = 4,
              plotlyOutput("ldh_gauge", height = "220px"),
              hr(),
              htmlOutput("braf_info")
          )
        ),
        fluidRow(
          box(title = "Risk Stratification Summary", status = "primary",
              solidHeader = TRUE, width = 12,
              DTOutput("patient_summary_tbl")
          )
        )
      ),

      ## ---- Tab 2: Drug PK ----
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "BRAF Inhibitor Dosing", status = "primary",
              solidHeader = TRUE, width = 3,
              selectInput("braf_drug", "Drug:",
                          c("Vemurafenib 960mg BID" = "vmfnb",
                            "Dabrafenib 150mg BID" = "dbfnb",
                            "No BRAFi" = "none")),
              numericInput("dose_B", "BRAFi Dose (mg):", 960, 0, 1920, 10),
              selectInput("tau_B_sel", "Dosing Interval:",
                          c("BID (q12h)" = "12",
                            "QD (q24h)" = "24"))
          ),
          box(title = "MEK Inhibitor Dosing", status = "success",
              solidHeader = TRUE, width = 3,
              selectInput("mek_drug", "Drug:",
                          c("Trametinib 2mg QD" = "tram",
                            "Cobimetinib 60mg 3/1" = "cobi",
                            "No MEKi" = "none")),
              numericInput("dose_M", "MEKi Dose (mg):", 0, 0, 120, 1)
          ),
          box(title = "Immune Checkpoint Inhibitor", status = "danger",
              solidHeader = TRUE, width = 3,
              selectInput("ici_drug", "Drug:",
                          c("Pembrolizumab 200mg" = "pembro",
                            "Nivolumab 240mg" = "nivo",
                            "No ICI" = "none")),
              numericInput("dose_I", "ICI Dose (mg):", 0, 0, 400, 10),
              selectInput("q_I", "ICI Interval:",
                          c("q3w (21 days)" = "21",
                            "q6w (42 days)" = "42",
                            "q4w (28 days)" = "28")),
              numericInput("dose_C4", "Ipilimumab Dose (mg, 0=none):", 0, 0, 600, 10)
          ),
          box(title = "Sequential Strategy", status = "warning",
              solidHeader = TRUE, width = 3,
              checkboxInput("use_seq", "Sequential therapy (BRAFi → ICI)", FALSE),
              numericInput("seq_switch", "Switch to ICI at week:", 24, 4, 52, 4),
              helpText("BRAFi/MEKi given first, then ICI from week N")
          )
        ),
        fluidRow(
          box(title = "BRAF/MEK Inhibitor Plasma Concentration",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("pk_braf_plot", height = "300px")),
          box(title = "Immune Checkpoint Inhibitor Concentration",
              status = "danger", solidHeader = TRUE, width = 6,
              plotlyOutput("pk_ici_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "ERK Inhibition Profile",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("erk_plot", height = "280px")),
          box(title = "PD-1 / CTLA-4 Receptor Occupancy",
              status = "info", solidHeader = TRUE, width = 6,
              plotlyOutput("ro_plot", height = "280px"))
        )
      ),

      ## ---- Tab 3: Tumor Response ----
      tabItem(tabName = "tab_tumor",
        fluidRow(
          box(title = "Tumor Burden Over Time",
              status = "danger", solidHeader = TRUE, width = 8,
              plotlyOutput("tumor_plot", height = "380px")),
          box(title = "Response Assessment", status = "warning",
              solidHeader = TRUE, width = 4,
              plotlyOutput("response_gauge", height = "200px"),
              hr(),
              htmlOutput("response_text")
          )
        ),
        fluidRow(
          box(title = "Waterfall: Best %Change from Baseline (Wk 12)",
              status = "info", solidHeader = TRUE, width = 6,
              plotlyOutput("waterfall_plot", height = "280px")),
          box(title = "Resistance Emergence",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("resist_plot", height = "280px"))
        )
      ),

      ## ---- Tab 4: Immune Dynamics ----
      tabItem(tabName = "tab_immune",
        fluidRow(
          box(title = "CD8+ TIL Dynamics",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("cd8_plot", height = "300px")),
          box(title = "Treg Fraction in TME",
              status = "danger", solidHeader = TRUE, width = 6,
              plotlyOutput("treg_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "PD-1 Receptor Occupancy",
              status = "info", solidHeader = TRUE, width = 4,
              plotlyOutput("pd1_ro_plot", height = "260px")),
          box(title = "CTLA-4 Receptor Occupancy",
              status = "warning", solidHeader = TRUE, width = 4,
              plotlyOutput("ctla4_ro_plot", height = "260px")),
          box(title = "IFN-γ in Tumor Microenvironment",
              status = "success", solidHeader = TRUE, width = 4,
              plotlyOutput("ifng_plot", height = "260px"))
        ),
        fluidRow(
          box(title = "Immune Summary Table (at Wk 12 & 24)",
              status = "primary", solidHeader = TRUE, width = 12,
              DTOutput("immune_tbl"))
        )
      ),

      ## ---- Tab 5: Scenario Comparison ----
      tabItem(tabName = "tab_scen",
        fluidRow(
          box(title = "Treatment Scenarios to Compare",
              status = "primary", solidHeader = TRUE, width = 3,
              checkboxGroupInput("scen_sel", "Select Scenarios:",
                choices = c("1. Untreated"             = "s1",
                            "2. Vemurafenib BID"       = "s2",
                            "3. Dabrafenib+Trametinib" = "s3",
                            "4. Pembrolizumab q3w"     = "s4",
                            "5. Nivo+Ipi (CheckMate)"  = "s5",
                            "6. BRAFi/MEKi→Pembro"     = "s6"),
                selected = c("s1","s3","s4","s5","s6")),
              actionButton("run_comparison", "Run Comparison",
                           style = "background-color:#3498db; color:white;
                                    font-weight:bold; margin-top:10px;")
          ),
          box(title = "Tumor Burden Comparison",
              status = "danger", solidHeader = TRUE, width = 9,
              plotlyOutput("comp_tumor_plot", height = "360px"))
        ),
        fluidRow(
          box(title = "% LDL-C Reduction Waterfall (Wk 24)",
              status = "info", solidHeader = TRUE, width = 6,
              plotlyOutput("comp_waterfall", height = "300px")),
          box(title = "ESC-Equivalent Goal Attainment Table",
              status = "success", solidHeader = TRUE, width = 6,
              DTOutput("comp_table"))
        )
      ),

      ## ---- Tab 6: Biomarkers ----
      tabItem(tabName = "tab_bio",
        fluidRow(
          box(title = "Serum LDH Over Time (Tumor Necrosis)",
              status = "danger", solidHeader = TRUE, width = 6,
              plotlyOutput("ldh_plot", height = "300px")),
          box(title = "S100B Protein (Melanoma Burden Marker)",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("s100b_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "TMB Impact on Immunotherapy Response",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("tmb_impact_plot", height = "300px")),
          box(title = "PD-L1 Expression Impact",
              status = "info", solidHeader = TRUE, width = 6,
              plotlyOutput("pdl1_impact_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "Biomarker Summary Table",
              status = "success", solidHeader = TRUE, width = 12,
              DTOutput("bio_tbl"))
        )
      )
    )
  )
)

## ============================================================
## Server
## ============================================================
server <- function(input, output, session) {

  ## ---- Simulation ----
  sim_data <- eventReactive(input$run_sim, {
    braf_v600 <- if (input$braf_status %in% c("V600E","V600K")) 1.0 else 0.0
    params <- list(
      BRAF_V600   = braf_v600,
      LDH0        = input$ldh_base,
      TMB_val     = input$tmb_val,
      PDL1_pct    = input$pdl1_pct
    )
    dose_M_val <- if (input$dose_M > 0 && input$mek_drug != "none") input$dose_M else 0
    dose_I_val <- if (input$dose_I > 0 && input$ici_drug != "none") input$dose_I else 0

    run_sim(mel_mod, params,
            dose_B       = if (input$braf_drug != "none") input$dose_B else 0,
            tau_B        = as.numeric(input$tau_B_sel),
            dose_M       = dose_M_val,
            dose_I       = dose_I_val,
            q_I          = as.numeric(input$q_I),
            dose_C4      = input$dose_C4,
            use_seq      = input$use_seq,
            seq_switch_wk = input$seq_switch,
            duration_wk  = input$duration_wk)
  })

  ## ---- LDH Gauge (Patient Profile) ----
  output$ldh_gauge <- renderPlotly({
    ldh <- input$ldh_base
    plot_ly(
      type = "indicator", mode = "gauge+number",
      value = ldh,
      title = list(text = "Baseline LDH (U/L)"),
      gauge = list(
        axis = list(range = list(0, 1500)),
        bar = list(color = "#c0392b"),
        steps = list(
          list(range = c(0, 200), color = "#27ae60"),
          list(range = c(200, 400), color = "#f39c12"),
          list(range = c(400, 1500), color = "#c0392b")
        ),
        threshold = list(line = list(color = "red", width = 4),
                         thickness = 0.75, value = 400)
      )
    ) %>% layout(height = 200, margin = list(t = 30, b = 10))
  })

  output$braf_info <- renderUI({
    info <- switch(input$braf_status,
      V600E = HTML("<b>BRAF V600E</b><br>Most common (45-50%). Highly sensitive to vemurafenib/dabrafenib. ORR ~60-70% with BRAFi+MEKi."),
      V600K = HTML("<b>BRAF V600K</b><br>~5% of melanoma. More UV-associated. Slightly lower BRAFi sensitivity than V600E."),
      NRAS  = HTML("<b>NRAS Q61 mutation</b><br>~20% of melanoma. No approved targeted therapy. MEK inhibitor (binimetinib) active. Good ICI candidate."),
      NF1   = HTML("<b>NF1 loss</b><br>~15% of melanoma. RAS hyperactivation. No BRAFi benefit. MEKi may help. Strong ICI candidate (high TMB)."),
      WT    = HTML("<b>WT/Other</b><br>~10% of melanoma. No targeted therapy. ICI is standard of care.")
    )
    info
  })

  output$patient_summary_tbl <- renderDT({
    df <- data.frame(
      Parameter       = c("BRAF Status", "Stage", "ECOG PS",
                          "Baseline LDH", "TMB", "PD-L1 TPS",
                          "Prior Therapy",
                          "Expected 1st-line ORR", "5-year OS estimate"),
      Value           = c(input$braf_status, input$stage,
                          input$ecog, paste(input$ldh_base, "U/L"),
                          paste(input$tmb_val, "mut/Mb"),
                          paste(input$pdl1_pct, "%"),
                          input$prior_tx,
                          if (input$braf_status %in% c("V600E","V600K"))
                            "~67% (BRAFi+MEKi) / 33-58% (ICI)"
                          else "~33-43% (PD-1 mono / nivo+ipi)",
                          if (input$ldh_base > 400)
                            "~20-25% (poor prognosis)"
                          else "~35-52% (standard prognosis)"),
      Note            = c("Determines targeted therapy eligibility",
                          "Staging per AJCC 8th edition",
                          "0-1: fit; 2: borderline",
                          "ULN 250 U/L; >400 = poor prognosis",
                          ">10 mut/Mb: ICI benefit likely",
                          ">1%: PD-1 benefit",
                          "May limit retreatment options",
                          "Unselected population estimate",
                          "CheckMate 067 / COMBI-d long-term data")
    )
    datatable(df, options = list(dom = "t", pageLength = 9),
              rownames = FALSE, class = "compact")
  })

  ## ---- PK Plots ----
  output$pk_braf_plot <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    p <- plot_ly(d, x = ~time_wk) %>%
      add_lines(y = ~Cb_out, name = "BRAFi (µg/mL)",
                line = list(color = "#2471A3", width = 2)) %>%
      add_lines(y = ~Cm_out / 1000, name = "MEKi (µg/mL, /1000 ng)",
                line = list(color = "#27AE60", width = 2,
                            dash = "dash")) %>%
      layout(title = "BRAFi & MEKi Plasma Concentration",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Concentration (µg/mL)"),
             legend = list(orientation = "h", y = -0.2))
    p
  })

  output$pk_ici_plot <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly(d, x = ~time_wk) %>%
      add_lines(y = ~Ci_out, name = "ICI (µg/mL)",
                line = list(color = "#922B21", width = 2)) %>%
      add_segments(x = 12, xend = 12, y = 0, yend = max(d$Ci_out) * 1.1,
                   line = list(color = "gray", dash = "dot"),
                   name = "EC50 PD-1") %>%
      layout(title = "ICI Plasma Concentration (q3w dosing)",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Concentration (µg/mL)"),
             legend = list(orientation = "h", y = -0.2))
  })

  output$erk_plot <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly(d, x = ~time_wk) %>%
      add_lines(y = ~ERK_pct, name = "ERK Activity (%)",
                line = list(color = "#F39C12", width = 2)) %>%
      add_segments(x = 0, xend = max(d$time_wk), y = 100, yend = 100,
                   line = list(color = "gray40", dash = "dash"),
                   name = "Baseline") %>%
      layout(title = "ERK Relative Activity (%)",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "ERK Activity (%)", range = c(0, 110)),
             legend = list(orientation = "h", y = -0.2))
  })

  output$ro_plot <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly(d, x = ~time_wk) %>%
      add_lines(y = ~PD1_pct, name = "PD-1 RO (%)",
                line = list(color = "#8E44AD", width = 2)) %>%
      add_lines(y = ~CTLA4_pct, name = "CTLA-4 RO (%)",
                line = list(color = "#E74C3C", width = 2, dash = "dash")) %>%
      layout(title = "Receptor Occupancy (%)",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Receptor Occupancy (%)", range = c(0, 105)),
             legend = list(orientation = "h", y = -0.2))
  })

  ## ---- Tumor Response Plots ----
  output$tumor_plot <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly(d, x = ~time_wk) %>%
      add_lines(y = ~Tumor_pct, name = "Tumor Burden",
                line = list(color = "#C0392B", width = 2)) %>%
      add_segments(x = 0, xend = max(d$time_wk), y = 100, yend = 100,
                   line = list(color = "gray60", dash = "dash"),
                   name = "Baseline") %>%
      add_segments(x = 0, xend = max(d$time_wk), y = 70, yend = 70,
                   line = list(color = "green", dash = "dot"),
                   name = "PR threshold (−30%)") %>%
      add_segments(x = 0, xend = max(d$time_wk), y = 120, yend = 120,
                   line = list(color = "red", dash = "dot"),
                   name = "PD threshold (+20%)") %>%
      layout(title = "Tumor Burden (% of Baseline SLD)",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Tumor Burden (%)"),
             legend = list(orientation = "h", y = -0.2))
  })

  output$response_gauge <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    wk12 <- d %>% filter(time_wk >= 11.5, time_wk <= 12.5) %>%
      summarise(T = mean(Tumor_pct)) %>% pull(T)
    pct_change <- wk12 - 100

    plot_ly(type = "indicator", mode = "gauge+number+delta",
            value = pct_change,
            delta = list(reference = 0),
            title = list(text = "Best %Change (Wk 12)"),
            gauge = list(
              axis = list(range = list(-100, 50)),
              bar = list(color = "#c0392b"),
              steps = list(
                list(range = c(-100, -30), color = "#27ae60"),
                list(range = c(-30, 20),   color = "#f39c12"),
                list(range = c(20, 50),    color = "#e74c3c")
              )
            )) %>%
      layout(height = 180, margin = list(t = 30, b = 5))
  })

  output$response_text <- renderUI({
    req(sim_data())
    d <- sim_data()
    wk24 <- d %>% filter(time_wk >= 23.5, time_wk <= 24.5) %>%
      summarise(T = mean(Tumor_pct)) %>% pull(T)
    pct <- round(wk24 - 100, 1)
    resp <- if (pct <= -100) "Complete Response (CR)" else
            if (pct <= -30)  "Partial Response (PR)" else
            if (pct < 20)    "Stable Disease (SD)" else
                             "Progressive Disease (PD)"
    color <- if (grepl("CR|PR", resp)) "green" else
             if (grepl("SD", resp)) "orange" else "red"
    HTML(paste0("<b>Response at Wk 24:</b><br>",
                "<span style='font-size:16px; color:", color, ";'>",
                resp, "</span><br>",
                "(%change = ", pct, "%)"))
  })

  output$waterfall_plot <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    wk12_val <- d %>% filter(time_wk >= 11.5, time_wk <= 12.5) %>%
      summarise(chg = mean(Tumor_pct) - 100) %>% pull(chg)
    df_wf <- data.frame(
      Label = "Current\nRegimen",
      Change = wk12_val
    )
    plot_ly(df_wf, x = ~Label, y = ~Change, type = "bar",
            marker = list(color = ifelse(wk12_val < -30, "#27ae60",
                                         ifelse(wk12_val < 20, "#f39c12",
                                                "#c0392b")))) %>%
      add_segments(x = -0.5, xend = 0.5, y = -30, yend = -30,
                   line = list(color = "green", dash = "dash"),
                   name = "PR threshold") %>%
      layout(title = "Tumor Change (Week 12)",
             xaxis = list(title = ""),
             yaxis = list(title = "% Change from Baseline"),
             showlegend = FALSE)
  })

  output$resist_plot <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly(d, x = ~time_wk) %>%
      add_lines(y = ~Resist_pct, name = "Resistance Score",
                line = list(color = "#795548", width = 2)) %>%
      layout(title = "Acquired Resistance Accumulation",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Resistance Score (%)", range = c(0, 105)),
             legend = list(orientation = "h", y = -0.2))
  })

  ## ---- Immune Dynamics ----
  output$cd8_plot <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly(d, x = ~time_wk) %>%
      add_lines(y = ~CD8_rel, name = "CD8+ TIL",
                line = list(color = "#2471A3", width = 2)) %>%
      add_segments(x = 0, xend = max(d$time_wk), y = 1, yend = 1,
                   line = list(color = "gray", dash = "dash"),
                   name = "Baseline") %>%
      layout(title = "CD8+ TIL Dynamics (relative to baseline)",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "CD8+ TIL (relative)"),
             legend = list(orientation = "h", y = -0.2))
  })

  output$treg_plot <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly(d, x = ~time_wk) %>%
      add_lines(y = ~Treg_pct, name = "Treg Fraction (%)",
                line = list(color = "#8E44AD", width = 2)) %>%
      layout(title = "Regulatory T Cell Fraction in TME",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Treg Fraction (%)", range = c(0, 30)),
             legend = list(orientation = "h", y = -0.2))
  })

  output$pd1_ro_plot <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly(d, x = ~time_wk) %>%
      add_lines(y = ~PD1_pct, line = list(color = "#1ABC9C", width = 2)) %>%
      layout(title = "PD-1 Receptor Occupancy",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "RO (%)", range = c(0, 105)))
  })

  output$ctla4_ro_plot <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly(d, x = ~time_wk) %>%
      add_lines(y = ~CTLA4_pct, line = list(color = "#E74C3C", width = 2)) %>%
      layout(title = "CTLA-4 Receptor Occupancy",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "RO (%)", range = c(0, 105)))
  })

  output$ifng_plot <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly(d, x = ~time_wk) %>%
      add_lines(y = ~IFNg_rel, line = list(color = "#27AE60", width = 2)) %>%
      layout(title = "IFN-γ in TME (relative)",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "IFN-γ (relative)"))
  })

  output$immune_tbl <- renderDT({
    req(sim_data())
    d <- sim_data()
    tbl <- d %>%
      filter(time_wk %in% c(0, 4, 12, 24, 52)) %>%
      group_by(time_wk) %>%
      summarise(
        `CD8+ TIL (rel)` = round(mean(CD8_rel), 3),
        `Treg (%)` = round(mean(Treg_pct), 1),
        `PD-1 RO (%)` = round(mean(PD1_pct), 1),
        `CTLA-4 RO (%)` = round(mean(CTLA4_pct), 1),
        `IFN-γ (rel)` = round(mean(IFNg_rel), 3),
        `Tumor (%)` = round(mean(Tumor_pct), 1),
        .groups = "drop"
      )
    datatable(tbl, options = list(dom = "t"), rownames = FALSE, class = "compact")
  })

  ## ---- Scenario Comparison ----
  comp_data <- eventReactive(input$run_comparison, {
    scen_params <- list(
      s1 = list(dose_B = 0, dose_M = 0, dose_I = 0, dose_C4 = 0,
                kd_I = 0.005, cif = 1.8, ctf = 1.4, treg_i = 0.4),
      s2 = list(dose_B = 960, dose_M = 0, dose_I = 0, dose_C4 = 0,
                kd_I = 0.004, cif = 1.8, ctf = 1.4, treg_i = 0.4),
      s3 = list(dose_B = 150, dose_M = 2,  dose_I = 0, dose_C4 = 0,
                kd_I = 0.0045, cif = 1.8, ctf = 1.4, treg_i = 0.4),
      s4 = list(dose_B = 0, dose_M = 0, dose_I = 200, dose_C4 = 0,
                kd_I = 0.006, cif = 1.8, ctf = 1.4, treg_i = 0.4),
      s5 = list(dose_B = 0, dose_M = 0, dose_I = 200, dose_C4 = 300,
                kd_I = 0.008, cif = 2.2, ctf = 1.8, treg_i = 0.6),
      s6 = list(dose_B = 150, dose_M = 2, dose_I = 200, dose_C4 = 0,
                kd_I = 0.006, cif = 1.8, ctf = 1.4, treg_i = 0.4)
    )
    scen_names <- c(
      s1 = "1. Untreated",
      s2 = "2. Vemurafenib BID",
      s3 = "3. Dabrafenib+Trametinib",
      s4 = "4. Pembrolizumab q3w",
      s5 = "5. Nivo+Ipi",
      s6 = "6. BRAFi/MEKi→Pembro"
    )
    selected <- input$scen_sel
    do.call(bind_rows, lapply(selected, function(s) {
      sp <- scen_params[[s]]
      p <- list(BRAF_V600 = 1.0, LDH0 = input$ldh_base,
                TMB_val = input$tmb_val, PDL1_pct = input$pdl1_pct,
                kd_immune = sp$kd_I, CD8_ICI_fac = sp$cif,
                CD8_CTLA4 = sp$ctf, Treg_ICI = sp$treg_i)
      use_s <- (s == "s6")
      run_sim(mel_mod, p,
              dose_B = sp$dose_B, tau_B = 12,
              dose_M = sp$dose_M, dose_I = sp$dose_I,
              q_I = 21, dose_C4 = sp$dose_C4,
              use_seq = use_s, seq_switch_wk = 24,
              duration_wk = input$duration_wk) %>%
        mutate(Scenario = scen_names[[s]])
    }))
  })

  output$comp_tumor_plot <- renderPlotly({
    req(comp_data())
    d <- comp_data()
    colors6 <- c("1. Untreated"             = "#E53935",
                 "2. Vemurafenib BID"       = "#FF8F00",
                 "3. Dabrafenib+Trametinib" = "#F9A825",
                 "4. Pembrolizumab q3w"     = "#1E88E5",
                 "5. Nivo+Ipi"              = "#00897B",
                 "6. BRAFi/MEKi→Pembro"     = "#8E24AA")
    p <- plot_ly()
    for (scn in unique(d$Scenario)) {
      sub <- d %>% filter(Scenario == scn)
      p <- p %>% add_lines(data = sub, x = ~time_wk, y = ~Tumor_pct,
                            name = scn,
                            line = list(color = colors6[[scn]], width = 2))
    }
    p %>%
      add_segments(x = 0, xend = max(d$time_wk), y = 100, yend = 100,
                   line = list(color = "gray50", dash = "dash"),
                   showlegend = FALSE) %>%
      add_segments(x = 0, xend = max(d$time_wk), y = 70, yend = 70,
                   line = list(color = "green3", dash = "dot"),
                   showlegend = FALSE) %>%
      layout(title = "Tumor Burden Comparison (All Scenarios)",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Tumor Burden (%)"),
             legend = list(orientation = "h", y = -0.25))
  })

  output$comp_waterfall <- renderPlotly({
    req(comp_data())
    d <- comp_data()
    wf <- d %>%
      filter(time_wk >= 23.5, time_wk <= 24.5) %>%
      group_by(Scenario) %>%
      summarise(Change = mean(Tumor_pct) - 100, .groups = "drop") %>%
      arrange(Change)
    plot_ly(wf, y = ~reorder(Scenario, Change), x = ~Change,
            type = "bar", orientation = "h",
            marker = list(color = ifelse(wf$Change < -30, "#27ae60",
                                          ifelse(wf$Change < 20, "#f39c12",
                                                 "#c0392b")))) %>%
      add_segments(x = -30, xend = -30, y = 0, yend = nrow(wf) + 1,
                   line = list(color = "green", dash = "dash"),
                   showlegend = FALSE) %>%
      add_segments(x = 20, xend = 20, y = 0, yend = nrow(wf) + 1,
                   line = list(color = "red", dash = "dash"),
                   showlegend = FALSE) %>%
      layout(title = "Best % Change (Week 24 Waterfall)",
             xaxis = list(title = "% Change from Baseline"),
             yaxis = list(title = ""))
  })

  output$comp_table <- renderDT({
    req(comp_data())
    d <- comp_data()
    tbl <- d %>%
      filter(time_wk %in% c(12, 24, 52)) %>%
      group_by(Scenario, time_wk) %>%
      summarise(`Tumor (%)` = round(mean(Tumor_pct), 1),
                `CD8+ TIL` = round(mean(CD8_rel), 2),
                `Resist (%)` = round(mean(Resist_pct), 1),
                .groups = "drop") %>%
      pivot_wider(names_from = time_wk,
                  values_from = c(`Tumor (%)`, `CD8+ TIL`, `Resist (%)`))
    datatable(tbl, options = list(dom = "t", scrollX = TRUE),
              rownames = FALSE, class = "compact")
  })

  ## ---- Biomarkers ----
  output$ldh_plot <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly(d, x = ~time_wk) %>%
      add_lines(y = ~LDH_val, name = "LDH (U/L)",
                line = list(color = "#E74C3C", width = 2)) %>%
      add_segments(x = 0, xend = max(d$time_wk), y = 250, yend = 250,
                   line = list(color = "red", dash = "dash"),
                   name = "ULN (250 U/L)") %>%
      layout(title = "Serum LDH Over Time",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "LDH (U/L)"),
             legend = list(orientation = "h", y = -0.2))
  })

  output$s100b_plot <- renderPlotly({
    req(sim_data())
    d <- sim_data()
    plot_ly(d, x = ~time_wk) %>%
      add_lines(y = ~S100_val, name = "S100B (µg/L)",
                line = list(color = "#8E44AD", width = 2)) %>%
      layout(title = "S100B Protein (Melanoma Marker)",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "S100B (µg/L)"),
             legend = list(orientation = "h", y = -0.2))
  })

  output$tmb_impact_plot <- renderPlotly({
    tmb_vals <- c(2, 5, 10, 20, 50)
    results <- lapply(tmb_vals, function(t) {
      p <- list(BRAF_V600 = 1.0, LDH0 = 200, TMB_val = t, PDL1_pct = 30,
                kd_immune = 0.006)
      ev_i <- ev(amt = 200, ii = 21 * 24, addl = 17, cmt = 5, rate = 400)
      d <- mel_mod %>% param(p) %>%
        mrgsim(ev = ev_i, end = 24 * 7 * 24, delta = 24) %>%
        as.data.frame() %>%
        mutate(time_wk = time / (7 * 24), TMB = t)
      d
    })
    all_d <- bind_rows(results)

    p <- plot_ly()
    for (t in tmb_vals) {
      sub <- all_d %>% filter(TMB == t, time_wk <= 24)
      p <- p %>% add_lines(data = sub, x = ~time_wk, y = ~Tumor_pct,
                            name = paste("TMB", t, "mut/Mb"))
    }
    p %>%
      layout(title = "TMB Impact on ICI Response (Pembrolizumab)",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Tumor Burden (%)"),
             legend = list(orientation = "h", y = -0.2))
  })

  output$pdl1_impact_plot <- renderPlotly({
    pdl1_vals <- c(0, 1, 10, 50, 100)
    results <- lapply(pdl1_vals, function(pl) {
      p <- list(BRAF_V600 = 1.0, LDH0 = 200, TMB_val = 10, PDL1_pct = pl,
                kd_immune = 0.006)
      ev_i <- ev(amt = 200, ii = 21 * 24, addl = 17, cmt = 5, rate = 400)
      d <- mel_mod %>% param(p) %>%
        mrgsim(ev = ev_i, end = 24 * 7 * 24, delta = 24) %>%
        as.data.frame() %>%
        mutate(time_wk = time / (7 * 24), PDL1 = pl)
      d
    })
    all_d <- bind_rows(results)
    p <- plot_ly()
    for (pl in pdl1_vals) {
      sub <- all_d %>% filter(PDL1 == pl, time_wk <= 24)
      p <- p %>% add_lines(data = sub, x = ~time_wk, y = ~Tumor_pct,
                            name = paste("PD-L1", pl, "%"))
    }
    p %>%
      layout(title = "PD-L1 Expression Impact on ICI Response",
             xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Tumor Burden (%)"),
             legend = list(orientation = "h", y = -0.2))
  })

  output$bio_tbl <- renderDT({
    req(sim_data())
    d <- sim_data()
    tbl <- d %>%
      filter(time_wk %in% c(0, 4, 8, 12, 24, 36, 52)) %>%
      group_by(`Week` = time_wk) %>%
      summarise(
        `Tumor (%)` = round(mean(Tumor_pct), 1),
        `LDH (U/L)` = round(mean(LDH_val), 0),
        `S100B (µg/L)` = round(mean(S100_val), 3),
        `ERK (%)` = round(mean(ERK_pct), 1),
        `Resist (%)` = round(mean(Resist_pct), 1),
        `CD8 TIL` = round(mean(CD8_rel), 3),
        `PD-1 RO (%)` = round(mean(PD1_pct), 1),
        .groups = "drop"
      )
    datatable(tbl, options = list(dom = "t", pageLength = 10),
              rownames = FALSE, class = "compact")
  })
}

shinyApp(ui, server)
