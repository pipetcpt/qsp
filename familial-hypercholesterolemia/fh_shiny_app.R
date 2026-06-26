## ============================================================
## Familial Hypercholesterolemia (FH) — Interactive Shiny App
## QSP Simulation Dashboard
## ============================================================
## Tabs: (1) Patient Profile  (2) Drug PK  (3) Lipid Response
##       (4) Biomarkers       (5) Scenario Comparison
##       (6) Genetic Profile & Clinical Risk
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

## ---- inline mrgsolve model --------------------------------
fh_code <- '
$PARAM
ka_S=0.45, F_S=0.20, Vd_S=2.0, CL_S=0.60, CL_H=0.45,
ka_P=0.0245, F_P=0.72, Vc_P=3.1, Vp_P=2.6, CL_P=0.0060, Q_P=0.0043,
kTMDD=0.00015, MW_P=144000,
ka_EZE=0.35, F_EZE=0.35, Vd_EZE=2.5, CL_EZE=0.15,
kout_HMG=0.28, EC50_SI=0.08, Emax_SI=0.92,
kin_LR=0.025, kout_LR=0.020, EC50_LD=0.30, Emax_LD=2.5,
kin_PK9=10.0, kout_PK9=0.024, kon_PK9=0.012, koff_PK9=1e-5,
EC50_PK9_LR=200, kP9_LDLR=0.010,
kVLDL_s=0.040, kIDL_s=0.085, kLDL_cl=0.0065, kLDL_cl_LR=0.028,
kHDL_eq=0.012, fEZE_abs=0.55, fLomit=0.80, fBemp=0.18,
LDL0_het=280, VLDL0=32, IDL0=15, HDL0=46, TG0=155, PCSK9_0=400,
LDLR_fxn=0.50, BW=75.0, Vmix=4.0,
DOSE_S=0, DOSE_P=0, DOSE_EZE=0, DOSE_LOMT=0, DOSE_BEMP=0,
USE_INCLISIRAN=0

$CMT GUT_S CENT_S LIV_S SC_PCSK9I CENT_PCSK9I PERI_PCSK9I COMP_PK9
     GUT_EZE CENT_EZE HMGCR_rel LDLR_rel PCSK9_pl
     VLDL_C IDL_C LDL_C HDL_C TG_C

$INIT
GUT_S=0,CENT_S=0,LIV_S=0,SC_PCSK9I=0,CENT_PCSK9I=0,
PERI_PCSK9I=0,COMP_PK9=0,GUT_EZE=0,CENT_EZE=0,
HMGCR_rel=1.0,LDLR_rel=1.0,PCSK9_pl=400,
VLDL_C=32,IDL_C=15,LDL_C=280,HDL_C=46,TG_C=155

$ODE
double C_liv_S = LIV_S / (Vd_S * BW * 0.26);
double C_sys_S = CENT_S / (Vd_S * BW);
double PCSK9i_nM = CENT_PCSK9I / (MW_P * 1e-6 * Vc_P);
double C_EZE = CENT_EZE / (Vd_EZE * BW);

dxdt_GUT_S  = -ka_S * GUT_S;
dxdt_CENT_S =  ka_S * F_S * GUT_S - (CL_S/(Vd_S*BW))*CENT_S - (CL_H/(Vd_S*BW*0.26))*CENT_S;
dxdt_LIV_S  =  ka_S * F_S * GUT_S * 0.80 - (CL_H/(Vd_S*BW*0.26))*LIV_S;

dxdt_SC_PCSK9I   = -ka_P * SC_PCSK9I;
dxdt_CENT_PCSK9I =  ka_P * F_P * SC_PCSK9I - (CL_P/Vc_P)*CENT_PCSK9I
                    -(Q_P/Vc_P)*CENT_PCSK9I + (Q_P/Vp_P)*PERI_PCSK9I
                    - kTMDD * PCSK9_pl * CENT_PCSK9I;
dxdt_PERI_PCSK9I = (Q_P/Vc_P)*CENT_PCSK9I - (Q_P/Vp_P)*PERI_PCSK9I;
dxdt_COMP_PK9    = kon_PK9 * PCSK9i_nM * (PCSK9_pl/1000.0)
                   - koff_PK9 * COMP_PK9 - kout_PK9 * COMP_PK9;

dxdt_GUT_EZE  = -ka_EZE * GUT_EZE;
dxdt_CENT_EZE =  ka_EZE * F_EZE * GUT_EZE - (CL_EZE/(Vd_EZE*BW))*CENT_EZE;

double Inh_S   = Emax_SI * C_liv_S / (EC50_SI + C_liv_S);
double Bemp_eff = DOSE_BEMP * fBemp * 0.55;
dxdt_HMGCR_rel = kout_HMG*(1.0 - Inh_S - Bemp_eff) - kout_HMG*HMGCR_rel;

double HMGCR_inh_frac = 1.0 - HMGCR_rel;
double LDLR_up = Emax_LD * HMGCR_inh_frac / (EC50_LD + HMGCR_inh_frac);
double PCSK9_free_tmp = PCSK9_pl - COMP_PK9*1000.0;
if(PCSK9_free_tmp < 0) PCSK9_free_tmp = 0;
double PCSK9_eff = kP9_LDLR * PCSK9_free_tmp / (EC50_PK9_LR + PCSK9_free_tmp);
double kin_LR_eff = kin_LR * LDLR_fxn * (1.0 + LDLR_up);
dxdt_LDLR_rel = kin_LR_eff - (kout_LR + PCSK9_eff)*LDLR_rel;

double inclisiran_red = 1.0 - USE_INCLISIRAN * 0.80;
dxdt_PCSK9_pl = kin_PK9 * (1.0 + 0.40*HMGCR_inh_frac) * inclisiran_red
                - kout_PK9 * PCSK9_pl
                - kon_PK9 * PCSK9i_nM * PCSK9_pl
                + koff_PK9 * COMP_PK9 * 1000.0;

double EZE_eff_vldl = 1.0 - fEZE_abs * C_EZE / (0.05 + C_EZE);
double LOMT_eff = 1.0 - DOSE_LOMT * fLomit;
dxdt_VLDL_C = kVLDL_s*VLDL0*EZE_eff_vldl*LOMT_eff
              - kVLDL_s*VLDL_C - 0.20*HMGCR_inh_frac*VLDL_C;
dxdt_IDL_C  = kVLDL_s*VLDL_C - kIDL_s*IDL_C;
double kLDL_total = kLDL_cl + kLDL_cl_LR*LDLR_rel;
dxdt_LDL_C  = kIDL_s*IDL_C - kLDL_total*LDL_C - DOSE_BEMP*fBemp*kLDL_cl*LDL_C;
dxdt_HDL_C  = kHDL_eq*(HDL0*(1.0+0.06*HMGCR_inh_frac) - HDL_C);
dxdt_TG_C   = kHDL_eq*(TG0*(1.0-0.15*HMGCR_inh_frac) - TG_C);

$TABLE
double NonHDL_C = LDL_C + VLDL_C + IDL_C;
double TC = LDL_C + VLDL_C + HDL_C + IDL_C;
double LDL_red_pct = (LDL0_het > 0) ? (LDL0_het - LDL_C)/LDL0_het*100.0 : 0;
double PCSK9_free = (PCSK9_pl - COMP_PK9*1000.0 < 0) ? 0 : PCSK9_pl - COMP_PK9*1000.0;
double C_pcsk9i = CENT_PCSK9I / Vc_P;
double C_stat   = CENT_S / (Vd_S * BW);
double C_eze    = CENT_EZE / (Vd_EZE * BW);

$CAPTURE NonHDL_C, TC, LDL_red_pct, PCSK9_free,
         C_pcsk9i, C_stat, C_eze,
         HMGCR_rel, LDLR_rel, PCSK9_pl,
         LDL_C, HDL_C, TG_C, VLDL_C, IDL_C
'

mod_base <- mcode("FH_SHINY", fh_code, quiet = TRUE)

## ---- simulation helper ------------------------------------
run_sim <- function(params, ldl0, ldlr_fxn,
                    dose_s, dose_p, dose_eze,
                    dose_lomt, dose_bemp, use_incl,
                    bw, n_weeks) {
  n_days <- n_weeks * 7
  init_v <- Init(mod_base)
  init_v["LDL_C"]    <- ldl0
  init_v["LDLR_rel"] <- ldlr_fxn
  init_v["HDL_C"]    <- params$HDL0
  init_v["TG_C"]     <- params$TG0
  init_v["VLDL_C"]   <- params$VLDL0
  init_v["PCSK9_pl"] <- params$PCSK9_0

  m2 <- mod_base %>%
    param(LDLR_fxn = ldlr_fxn, LDL0_het = ldl0, BW = bw,
          HDL0 = params$HDL0, TG0 = params$TG0, VLDL0 = params$VLDL0,
          DOSE_S = dose_s, DOSE_P = dose_p, DOSE_EZE = dose_eze,
          DOSE_LOMT = dose_lomt, DOSE_BEMP = dose_bemp,
          USE_INCLISIRAN = use_incl) %>%
    init(init_v)

  evs <- list()
  if (dose_s   > 0) evs[["s"]] <- ev(amt = dose_s,   cmt = "GUT_S",     ii = 24,      addl = n_days - 1)
  if (dose_p   > 0) evs[["p"]] <- ev(amt = dose_p,   cmt = "SC_PCSK9I", ii = 28 * 24, addl = floor(n_days/28))
  if (dose_eze > 0) evs[["e"]] <- ev(amt = dose_eze,  cmt = "GUT_EZE",   ii = 24,      addl = n_days - 1)

  ev_all <- if (length(evs) > 0) Reduce(c, evs) else NULL

  out <- if (!is.null(ev_all)) {
    mrgsim(m2, ev_all, end = n_days * 24, delta = 6)
  } else {
    mrgsim(m2, end = n_days * 24, delta = 6)
  }
  as.data.frame(out) %>% mutate(time_days = time / 24, time_weeks = time / (24 * 7))
}

## ===========================================================
## UI
## ===========================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "FH QSP Dashboard"),
  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",      tabName = "tab_patient",  icon = icon("user")),
      menuItem("Drug PK",              tabName = "tab_pk",       icon = icon("pills")),
      menuItem("Lipid Response",       tabName = "tab_lipid",    icon = icon("chart-line")),
      menuItem("Biomarkers",           tabName = "tab_bio",      icon = icon("vial")),
      menuItem("Scenario Comparison",  tabName = "tab_compare",  icon = icon("balance-scale")),
      menuItem("Genetic Risk Profile", tabName = "tab_genetic",  icon = icon("dna"))
    ),
    hr(),
    h4("Patient Parameters", style = "padding-left:15px; color:#ECF0F1;"),
    selectInput("fh_type", "FH Genotype:",
                choices = c("Heterozygous FH (HetFH)" = "het",
                            "Homozygous FH (HomFH)"   = "hom",
                            "FDB (APOB R3527Q)"        = "fdb",
                            "Normal (control)"         = "normal"),
                selected = "het"),
    numericInput("bw", "Body Weight (kg):", value = 75, min = 40, max = 150),
    numericInput("age", "Age (years):", value = 45, min = 10, max = 80),
    selectInput("sex", "Sex:", choices = c("Male" = "M", "Female" = "F")),
    numericInput("hdl_base", "Baseline HDL-C (mg/dL):", value = 46, min = 20, max = 100),
    numericInput("tg_base", "Baseline TG (mg/dL):", value = 155, min = 50, max = 500),
    hr(),
    h4("Treatment", style = "padding-left:15px; color:#ECF0F1;"),
    sliderInput("dose_statin", "Statin Dose (mg/day):",
                min = 0, max = 80, value = 0, step = 5),
    checkboxInput("use_eze", "Ezetimibe 10 mg/day", value = FALSE),
    checkboxInput("use_pcsk9i", "Evolocumab 420 mg q4w", value = FALSE),
    checkboxInput("use_incl", "Inclisiran (siRNA, ~6-monthly)", value = FALSE),
    checkboxInput("use_lomt", "Lomitapide (HomFH only)", value = FALSE),
    checkboxInput("use_bemp", "Bempedoic Acid 180 mg/day", value = FALSE),
    sliderInput("n_weeks", "Simulation Duration (weeks):",
                min = 4, max = 104, value = 52, step = 4),
    actionButton("run_btn", "Run Simulation", class = "btn-success btn-block",
                 icon = icon("play"))
  ),

  dashboardBody(
    tabItems(

      ## ---- TAB 1: Patient Profile ----
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "FH Patient Overview", width = 6, status = "primary",
            uiOutput("patient_summary")
          ),
          box(title = "LDL-C Goal Attainment", width = 6, status = "warning",
            plotlyOutput("ldl_gauge", height = "250px"),
            uiOutput("goal_status")
          )
        ),
        fluidRow(
          box(title = "Baseline Lipid Panel", width = 6, status = "info",
            tableOutput("baseline_table")
          ),
          box(title = "10-year CVD Risk Estimate", width = 6, status = "danger",
            plotlyOutput("cvd_risk_plot", height = "250px")
          )
        )
      ),

      ## ---- TAB 2: Drug PK ----
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Statin — Hepatic Concentration", width = 6, status = "primary",
            plotlyOutput("statin_pk_plot", height = "300px")
          ),
          box(title = "PCSK9 Inhibitor — Plasma Concentration", width = 6, status = "success",
            plotlyOutput("pcsk9i_pk_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Ezetimibe Plasma Level", width = 6, status = "warning",
            plotlyOutput("eze_pk_plot", height = "300px")
          ),
          box(title = "HMGCR Relative Activity", width = 6, status = "info",
            plotlyOutput("hmgcr_plot", height = "300px")
          )
        )
      ),

      ## ---- TAB 3: Lipid Response ----
      tabItem(tabName = "tab_lipid",
        fluidRow(
          box(title = "LDL-C Time Course", width = 8, status = "danger",
            plotlyOutput("ldl_plot", height = "350px")
          ),
          box(title = "LDL-C Summary", width = 4, status = "danger",
            uiOutput("ldl_summary_box")
          )
        ),
        fluidRow(
          box(title = "HDL-C & TG Response", width = 6, status = "success",
            plotlyOutput("hdl_tg_plot", height = "300px")
          ),
          box(title = "Full Lipid Panel (Week 52)", width = 6, status = "info",
            plotlyOutput("lipid_panel_bar", height = "300px")
          )
        )
      ),

      ## ---- TAB 4: Biomarkers ----
      tabItem(tabName = "tab_bio",
        fluidRow(
          box(title = "Plasma PCSK9 (ng/mL)", width = 6, status = "primary",
            plotlyOutput("pcsk9_plot", height = "300px")
          ),
          box(title = "LDLR Surface Expression (%)", width = 6, status = "success",
            plotlyOutput("ldlr_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Non-HDL-C & ApoB Surrogate", width = 6, status = "warning",
            plotlyOutput("nonhdl_plot", height = "300px")
          ),
          box(title = "Biomarker Summary Table (Week 52)", width = 6, status = "info",
            DTOutput("bio_table")
          )
        )
      ),

      ## ---- TAB 5: Scenario Comparison ----
      tabItem(tabName = "tab_compare",
        fluidRow(
          box(title = "LDL-C Across Treatment Scenarios", width = 12, status = "primary",
            plotlyOutput("scenario_ldl_plot", height = "400px")
          )
        ),
        fluidRow(
          box(title = "% LDL Reduction — Waterfall Chart", width = 7, status = "danger",
            plotlyOutput("waterfall_plot", height = "350px")
          ),
          box(title = "Scenario Comparison Table", width = 5, status = "info",
            DTOutput("scenario_table")
          )
        )
      ),

      ## ---- TAB 6: Genetic Risk ----
      tabItem(tabName = "tab_genetic",
        fluidRow(
          box(title = "LDLR Mutation Class Effects", width = 6, status = "primary",
            plotlyOutput("mutation_plot", height = "350px")
          ),
          box(title = "PCSK9 Variant Impact on LDL-C", width = 6, status = "warning",
            plotlyOutput("pcsk9_variant_plot", height = "350px")
          )
        ),
        fluidRow(
          box(title = "CVD Risk Reduction by Treatment Duration", width = 12, status = "danger",
            plotlyOutput("cvd_reduction_plot", height = "300px")
          )
        )
      )

    ) # end tabItems
  )   # end dashboardBody
)     # end dashboardPage

## ===========================================================
## SERVER
## ===========================================================
server <- function(input, output, session) {

  ## ---- Reactive: genotype parameters ----
  geno_params <- reactive({
    switch(input$fh_type,
      "het"    = list(ldlr_fxn = 0.50, ldl0 = 280, label = "Heterozygous FH"),
      "hom"    = list(ldlr_fxn = 0.05, ldl0 = 680, label = "Homozygous FH"),
      "fdb"    = list(ldlr_fxn = 0.60, ldl0 = 210, label = "Familial Defective ApoB"),
      "normal" = list(ldlr_fxn = 1.00, ldl0 = 120, label = "Normal (non-FH)")
    )
  })

  ## ---- Reactive: run simulation on button click ----
  sim_out <- eventReactive(input$run_btn, {
    gp <- geno_params()
    run_sim(
      params    = list(HDL0 = input$hdl_base, TG0 = input$tg_base,
                       VLDL0 = 32, PCSK9_0 = 400),
      ldl0      = gp$ldl0,
      ldlr_fxn  = gp$ldlr_fxn,
      dose_s    = input$dose_statin,
      dose_p    = if (input$use_pcsk9i) 420 else 0,
      dose_eze  = if (input$use_eze)    10  else 0,
      dose_lomt = if (input$use_lomt)   1   else 0,
      dose_bemp = if (input$use_bemp)   1   else 0,
      use_incl  = if (input$use_incl)   1   else 0,
      bw        = input$bw,
      n_weeks   = input$n_weeks
    )
  })

  ## ---- TAB 1: Patient Summary ----
  output$patient_summary <- renderUI({
    gp <- geno_params()
    tags$div(
      tags$h4(gp$label, style = "color:#2E86C1;"),
      tags$ul(
        tags$li(sprintf("Age: %d yrs | Sex: %s | BW: %.0f kg",
                        input$age, input$sex, input$bw)),
        tags$li(sprintf("Baseline LDL-C: %.0f mg/dL", gp$ldl0)),
        tags$li(sprintf("Baseline HDL-C: %.0f mg/dL", input$hdl_base)),
        tags$li(sprintf("Baseline TG: %.0f mg/dL", input$tg_base)),
        tags$li(sprintf("LDLR Residual Function: %.0f%%", gp$ldlr_fxn * 100)),
        tags$li(sprintf("Expected ASCVD Risk: %s",
                        ifelse(gp$ldl0 > 300, "Very High (>20%/10yr)", "High (10-20%/10yr)")))
      )
    )
  })

  output$baseline_table <- renderTable({
    gp <- geno_params()
    data.frame(
      Parameter    = c("LDL-C", "HDL-C", "Triglycerides", "Non-HDL-C",
                       "Total Cholesterol", "ApoB100 (est.)", "PCSK9 (plasma)"),
      Value        = c(gp$ldl0, input$hdl_base, input$tg_base,
                       gp$ldl0 + 32 + 15, gp$ldl0 + input$hdl_base + 32 + 15,
                       round(0.8 * gp$ldl0, 0), "400 ng/mL"),
      Unit         = c("mg/dL","mg/dL","mg/dL","mg/dL","mg/dL","mg/dL","ng/mL"),
      `ESC Target` = c("< 55 (very high risk)","≥ 40","< 150","< 85","< 200","< 65","—")
    )
  })

  output$ldl_gauge <- renderPlotly({
    req(sim_out())
    df <- sim_out()
    ldl_now <- tail(df$LDL_C, 1)
    gp <- geno_params()
    plot_ly(
      type = "indicator", mode = "gauge+number+delta",
      value = ldl_now,
      delta = list(reference = gp$ldl0, decreasing = list(color = "green")),
      gauge = list(
        axis = list(range = list(0, max(gp$ldl0 * 1.1, 200))),
        bar  = list(color = ifelse(ldl_now <= 55, "green",
                                   ifelse(ldl_now <= 70, "orange", "red"))),
        steps = list(
          list(range = c(0, 55), color = "#A9DFBF"),
          list(range = c(55, 70), color = "#F9E79F"),
          list(range = c(70, 200), color = "#FADBD8")
        ),
        threshold = list(line = list(color = "darkred", width = 4), value = 70)
      ),
      title = list(text = "LDL-C (mg/dL) at Week 52")
    ) %>% layout(margin = list(t = 30))
  })

  output$goal_status <- renderUI({
    req(sim_out())
    df <- sim_out()
    ldl_now <- tail(df$LDL_C, 1)
    style55 <- if (ldl_now <= 55) "color:green;font-weight:bold;" else "color:red;"
    style70 <- if (ldl_now <= 70) "color:green;font-weight:bold;" else "color:orange;"
    tags$div(
      tags$p(sprintf("LDL-C at week %d: %.1f mg/dL", input$n_weeks, ldl_now)),
      tags$p(style = style55, sprintf("ESC <55 mg/dL goal: %s", if (ldl_now <= 55) "ACHIEVED ✓" else "NOT MET")),
      tags$p(style = style70, sprintf("ESC <70 mg/dL goal: %s", if (ldl_now <= 70) "ACHIEVED ✓" else "NOT MET"))
    )
  })

  output$cvd_risk_plot <- renderPlotly({
    gp <- geno_params()
    ldl_vals <- seq(50, 500, 10)
    cvd_risk <- 0.01 * exp(0.0035 * ldl_vals - 0.01 * input$hdl_base)
    df_risk <- data.frame(LDL = ldl_vals, Risk = pmin(cvd_risk * 100, 50))
    plot_ly(df_risk, x = ~LDL, y = ~Risk, type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(231,76,60,0.3)",
            line = list(color = "#E74C3C")) %>%
      add_segments(x = gp$ldl0, xend = gp$ldl0, y = 0,
                   yend = pmin(0.01 * exp(0.0035*gp$ldl0 - 0.01*input$hdl_base) * 100, 50),
                   line = list(color = "navy", dash = "dot")) %>%
      layout(title = "10-yr CVD Risk vs LDL-C (illustrative)",
             xaxis = list(title = "LDL-C (mg/dL)"),
             yaxis = list(title = "10-yr CVD Risk (%)"))
  })

  ## ---- TAB 2: Drug PK ----
  output$statin_pk_plot <- renderPlotly({
    req(sim_out()); df <- sim_out()
    plot_ly(df, x = ~time_weeks, y = ~C_stat, type = "scatter", mode = "lines",
            line = list(color = "#2E86C1")) %>%
      layout(title = "Statin Plasma Concentration",
             xaxis = list(title = "Weeks"),
             yaxis = list(title = "Concentration (μg/mL)"))
  })

  output$pcsk9i_pk_plot <- renderPlotly({
    req(sim_out()); df <- sim_out()
    plot_ly(df, x = ~time_weeks, y = ~C_pcsk9i, type = "scatter", mode = "lines",
            line = list(color = "#27AE60")) %>%
      layout(title = "Evolocumab Plasma (mg/L)",
             xaxis = list(title = "Weeks"),
             yaxis = list(title = "Concentration (mg/L)"))
  })

  output$eze_pk_plot <- renderPlotly({
    req(sim_out()); df <- sim_out()
    plot_ly(df, x = ~time_weeks, y = ~C_eze, type = "scatter", mode = "lines",
            line = list(color = "#8E44AD")) %>%
      layout(title = "Ezetimibe Plasma (mg/L)",
             xaxis = list(title = "Weeks"),
             yaxis = list(title = "Concentration (mg/L)"))
  })

  output$hmgcr_plot <- renderPlotly({
    req(sim_out()); df <- sim_out()
    plot_ly(df, x = ~time_weeks, y = ~I(HMGCR_rel * 100), type = "scatter", mode = "lines",
            line = list(color = "#E67E22")) %>%
      layout(title = "HMGCR Relative Activity (%)",
             xaxis = list(title = "Weeks"),
             yaxis = list(title = "HMGCR Activity (% baseline)", range = c(0, 110)))
  })

  ## ---- TAB 3: Lipid Response ----
  output$ldl_plot <- renderPlotly({
    req(sim_out()); df <- sim_out()
    plot_ly(df, x = ~time_weeks, y = ~LDL_C, type = "scatter", mode = "lines",
            line = list(color = "#E74C3C", width = 3), name = "LDL-C") %>%
      add_segments(x = 0, xend = max(df$time_weeks), y = 70, yend = 70,
                   line = list(color = "navy", dash = "dash"), name = "70 mg/dL goal") %>%
      add_segments(x = 0, xend = max(df$time_weeks), y = 55, yend = 55,
                   line = list(color = "darkred", dash = "dot"), name = "55 mg/dL goal") %>%
      layout(title = "LDL-C Time Course",
             xaxis = list(title = "Weeks"),
             yaxis = list(title = "LDL-C (mg/dL)"))
  })

  output$ldl_summary_box <- renderUI({
    req(sim_out())
    df <- sim_out()
    gp <- geno_params()
    ldl_final <- tail(df$LDL_C, 1)
    ldl_red   <- (gp$ldl0 - ldl_final) / gp$ldl0 * 100
    tags$div(
      tags$p(strong("Baseline LDL-C: "), sprintf("%.0f mg/dL", gp$ldl0)),
      tags$p(strong("Final LDL-C: "),    sprintf("%.1f mg/dL", ldl_final)),
      tags$p(strong("Reduction: "),      sprintf("%.1f%%", max(ldl_red, 0))),
      tags$hr(),
      tags$p(strong("ESC <55 mg/dL: "),
             if (ldl_final <= 55) tags$span("ACHIEVED", style = "color:green;") else
               tags$span("NOT MET", style = "color:red;")),
      tags$p(strong("ESC <70 mg/dL: "),
             if (ldl_final <= 70) tags$span("ACHIEVED", style = "color:green;") else
               tags$span("NOT MET", style = "color:orange;"))
    )
  })

  output$hdl_tg_plot <- renderPlotly({
    req(sim_out()); df <- sim_out()
    plot_ly(df, x = ~time_weeks, y = ~HDL_C, type = "scatter", mode = "lines",
            name = "HDL-C", line = list(color = "#27AE60")) %>%
      add_trace(y = ~TG_C, name = "TG", line = list(color = "#E67E22")) %>%
      layout(title = "HDL-C & Triglycerides",
             xaxis = list(title = "Weeks"),
             yaxis = list(title = "mg/dL"))
  })

  output$lipid_panel_bar <- renderPlotly({
    req(sim_out()); df <- sim_out()
    gp <- geno_params()
    final <- tail(df, 1)
    df_bar <- data.frame(
      Lipid = c("LDL-C","HDL-C","TG","VLDL-C","IDL-C"),
      Baseline = c(gp$ldl0, input$hdl_base, input$tg_base, 32, 15),
      Week52   = c(final$LDL_C, final$HDL_C, final$TG_C, final$VLDL_C, final$IDL_C)
    ) %>% pivot_longer(-Lipid, names_to = "Timepoint", values_to = "mgdL")
    plot_ly(df_bar, x = ~Lipid, y = ~mgdL, color = ~Timepoint,
            type = "bar", colors = c("Baseline" = "#AED6F1", "Week52" = "#2E86C1")) %>%
      layout(barmode = "group", title = "Lipid Panel: Baseline vs Week 52",
             yaxis = list(title = "mg/dL"))
  })

  ## ---- TAB 4: Biomarkers ----
  output$pcsk9_plot <- renderPlotly({
    req(sim_out()); df <- sim_out()
    plot_ly(df, x = ~time_weeks, y = ~PCSK9_pl, type = "scatter", mode = "lines",
            name = "Total PCSK9", line = list(color = "#E67E22")) %>%
      add_trace(y = ~PCSK9_free, name = "Free PCSK9", line = list(color = "#E74C3C", dash = "dash")) %>%
      layout(title = "Plasma PCSK9 Concentration",
             xaxis = list(title = "Weeks"),
             yaxis = list(title = "PCSK9 (ng/mL)"))
  })

  output$ldlr_plot <- renderPlotly({
    req(sim_out()); df <- sim_out()
    plot_ly(df, x = ~time_weeks, y = ~I(LDLR_rel * 100), type = "scatter", mode = "lines",
            line = list(color = "#27AE60", width = 3)) %>%
      layout(title = "Hepatic LDLR Surface Expression",
             xaxis = list(title = "Weeks"),
             yaxis = list(title = "LDLR (% of normal)"))
  })

  output$nonhdl_plot <- renderPlotly({
    req(sim_out()); df <- sim_out()
    plot_ly(df, x = ~time_weeks, y = ~NonHDL_C, type = "scatter", mode = "lines",
            name = "Non-HDL-C", line = list(color = "#8E44AD")) %>%
      add_segments(x = 0, xend = max(df$time_weeks), y = 85, yend = 85,
                   line = list(color = "darkred", dash = "dash"), name = "85 mg/dL target") %>%
      layout(title = "Non-HDL-C",
             xaxis = list(title = "Weeks"),
             yaxis = list(title = "Non-HDL-C (mg/dL)"))
  })

  output$bio_table <- renderDT({
    req(sim_out())
    df <- sim_out()
    gp <- geno_params()
    final <- tail(df, 1)
    datatable(data.frame(
      Biomarker  = c("LDL-C","HDL-C","TG","VLDL-C","Non-HDL-C","Total Chol",
                     "PCSK9 (total)","PCSK9 (free)","LDLR Surface","HMGCR Activity"),
      Baseline   = c(gp$ldl0, input$hdl_base, input$tg_base, 32,
                     gp$ldl0+32+15, gp$ldl0+input$hdl_base+32+15,
                     400, 400, round(gp$ldlr_fxn*100), 100),
      Week52     = round(c(final$LDL_C, final$HDL_C, final$TG_C, final$VLDL_C,
                           final$NonHDL_C, final$TC,
                           final$PCSK9_pl, final$PCSK9_free,
                           final$LDLR_rel*100, final$HMGCR_rel*100), 1),
      Unit       = c(rep("mg/dL", 6), rep("ng/mL", 2), rep("%", 2))
    ), options = list(pageLength = 10, dom = "t"))
  })

  ## ---- TAB 5: Scenario Comparison ----
  scenario_results <- reactive({
    gp <- geno_params()
    params_base <- list(HDL0 = input$hdl_base, TG0 = input$tg_base,
                        VLDL0 = 32, PCSK9_0 = 400)
    scens <- list(
      list(label = "No Treatment",          ds = 0, dp = 0, de = 0, dl = 0, db = 0, di = 0),
      list(label = "Statin 40 mg/d",        ds = 40,dp = 0, de = 0, dl = 0, db = 0, di = 0),
      list(label = "Statin + EZE",          ds = 40,dp = 0, de = 10,dl = 0, db = 0, di = 0),
      list(label = "PCSK9i (Evolocumab)",   ds = 0, dp = 420,de=0, dl = 0, db = 0, di = 0),
      list(label = "Statin + PCSK9i",       ds = 40,dp = 420,de=0, dl = 0, db = 0, di = 0),
      list(label = "Statin+EZE+PCSK9i",     ds = 40,dp = 420,de=10,dl = 0, db = 0, di = 0),
      list(label = "Inclisiran",             ds = 0, dp = 0, de = 0, dl = 0, db = 0, di = 1),
      list(label = "Lomitapide+Statin+PCSK9i", ds=40,dp=420,de=0, dl = 1, db = 0, di = 0)
    )
    results <- lapply(scens, function(sc) {
      tryCatch({
        out <- run_sim(params_base, gp$ldl0, gp$ldlr_fxn,
                       sc$ds, sc$dp, sc$de, sc$dl, sc$db, sc$di,
                       input$bw, input$n_weeks)
        final <- tail(out, 1)
        list(label = sc$label, sim = out, ldl_final = final$LDL_C,
             ldl_red = (gp$ldl0 - final$LDL_C) / gp$ldl0 * 100)
      }, error = function(e) NULL)
    })
    Filter(Negate(is.null), results)
  })

  output$scenario_ldl_plot <- renderPlotly({
    res <- scenario_results()
    colors_sc <- RColorBrewer::brewer.pal(min(length(res), 9), "Set1")
    p <- plot_ly()
    for (i in seq_along(res)) {
      r <- res[[i]]
      p <- add_trace(p, data = r$sim, x = ~time_weeks, y = ~LDL_C,
                     type = "scatter", mode = "lines", name = r$label,
                     line = list(color = colors_sc[((i-1) %% 9) + 1], width = 2))
    }
    p %>%
      add_segments(x = 0, xend = input$n_weeks, y = 70, yend = 70,
                   line = list(color = "navy", dash = "dash"), name = "70 mg/dL") %>%
      add_segments(x = 0, xend = input$n_weeks, y = 55, yend = 55,
                   line = list(color = "darkred", dash = "dot"), name = "55 mg/dL") %>%
      layout(title = "LDL-C: All Scenarios",
             xaxis = list(title = "Weeks"),
             yaxis = list(title = "LDL-C (mg/dL)"))
  })

  output$waterfall_plot <- renderPlotly({
    res <- scenario_results()
    df_wf <- data.frame(
      Scenario = sapply(res, `[[`, "label"),
      Red      = sapply(res, `[[`, "ldl_red")
    ) %>% arrange(desc(Red))
    plot_ly(df_wf, y = ~reorder(Scenario, Red), x = ~Red,
            type = "bar", orientation = "h",
            marker = list(color = ifelse(df_wf$Red >= 50, "#27AE60",
                                         ifelse(df_wf$Red >= 30, "#F39C12", "#E74C3C")))) %>%
      layout(title = "% LDL-C Reduction at Week 52",
             xaxis = list(title = "% Reduction from baseline"),
             yaxis = list(title = ""))
  })

  output$scenario_table <- renderDT({
    res <- scenario_results()
    gp <- geno_params()
    df_t <- data.frame(
      Scenario    = sapply(res, `[[`, "label"),
      LDL_Final   = round(sapply(res, `[[`, "ldl_final"), 1),
      `% Reduction` = round(sapply(res, `[[`, "ldl_red"), 1),
      Goal_55     = ifelse(sapply(res, function(r) r$ldl_final <= 55), "Yes", "No"),
      Goal_70     = ifelse(sapply(res, function(r) r$ldl_final <= 70), "Yes", "No")
    )
    datatable(df_t, options = list(pageLength = 10, dom = "t"),
              rownames = FALSE) %>%
      formatStyle("Goal_55", backgroundColor = styleEqual("Yes", "#A9DFBF")) %>%
      formatStyle("Goal_70", backgroundColor = styleEqual("Yes", "#D5F5E3"))
  })

  ## ---- TAB 6: Genetic ----
  output$mutation_plot <- renderPlotly({
    mutations <- data.frame(
      Class  = c("Normal", "Class 5\n(recycling)", "Class 4\n(clustering)",
                  "Class 3\n(binding)", "Class 2\n(transport)", "Class 1\n(null)"),
      LDLR   = c(100, 60, 35, 25, 10, 2),
      LDL_C  = c(120, 180, 230, 260, 300, 350)
    )
    plot_ly(mutations, x = ~LDLR, y = ~LDL_C, text = ~Class,
            type = "scatter", mode = "markers+text",
            textposition = "top right",
            marker = list(size = 15,
                          color = c("#27AE60","#F39C12","#E67E22","#E74C3C","#922B21","#7B241C"),
                          line  = list(width = 2, color = "white"))) %>%
      layout(title = "LDLR Mutation Class → LDL-C",
             xaxis = list(title = "LDLR Residual Function (%)", range = c(-5, 110)),
             yaxis = list(title = "Steady-state LDL-C (mg/dL)"))
  })

  output$pcsk9_variant_plot <- renderPlotly({
    variants <- data.frame(
      Variant = c("D374Y GoF\n(most severe)", "PCSK9 WT",
                  "R46L LoF\n(protective)", "Y142X LoF\n(very protective)"),
      PCSK9   = c(650, 400, 290, 220),
      LDL_C   = c(420, 280, 200, 140)
    )
    plot_ly(variants, x = ~PCSK9, y = ~LDL_C, text = ~Variant,
            type = "scatter", mode = "markers+text",
            textposition = "top right",
            marker = list(size = 16,
                          color = c("#E74C3C","#E67E22","#27AE60","#1E8449"))) %>%
      layout(title = "PCSK9 Variant → Plasma LDL-C",
             xaxis = list(title = "Plasma PCSK9 (ng/mL)"),
             yaxis = list(title = "LDL-C (mg/dL)"))
  })

  output$cvd_reduction_plot <- renderPlotly({
    # Mendelian randomisation-inspired: every 39 mg/dL lifetime LDL-C reduction = 55% CVD reduction
    years <- seq(0, 40, 1)
    ldl_red_per_year <- 80  # ~80 mg/dL lifetime LDL-C reduction
    cvd_red <- 1 - exp(-years * log(1.55) / 22)  # CTT meta-analysis: 22% per mmol/L per year
    plot_ly(x = ~years, y = ~(cvd_red * 100),
            type = "scatter", mode = "lines+markers",
            line = list(color = "#2E86C1", width = 3)) %>%
      layout(title = "Estimated Cumulative CVD Risk Reduction with LDL-C Treatment\n(based on CTT meta-analysis: 22% per mmol/L per year)",
             xaxis = list(title = "Years of Treatment"),
             yaxis = list(title = "CVD Risk Reduction (%)", range = c(0, 70)))
  })

}

## ---- Launch ----
shinyApp(ui, server)
