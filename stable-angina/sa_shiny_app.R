## =============================================================================
## Stable Angina (Chronic Coronary Syndrome) — Interactive QSP Shiny App
## mrgsolve-based simulation dashboard
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)
library(plotly)
library(patchwork)

# ─── mrgsolve model (same as sa_mrgsolve_model.R) ────────────────────────────
model_code <- '
$PROB Stable Angina QSP — Shiny App Version

$CMT
 GUT_BB CENTRAL_BB PERIPH_BB
 GUT_CCB CENTRAL_CCB PERIPH_CCB
 GUT_RAN CENTRAL_RAN
 GUT_IVA CENTRAL_IVA
 GUT_NIT CENTRAL_NIT NIT_TOL
 HR_STATE SBP_STATE CBF_STATE O2IMBAL ISCHEMIA ANGINA_SCORE EX_CAP PLAQUE

$PARAM
 BW=75 AGE=62 HR0=80 SBP0=145 CBF0=250 STENOSIS=70
 KA_BB=1.0  F_BB=0.85  CL_BB=12.5 V1_BB=50   Q_BB=8    V2_BB=250
 EC50_BB=30 EMAX_BB=0.25 EC50_SBP_BB=35 EMAX_SBP_BB=0.12
 KA_CCB=0.5 F_CCB=0.64 CL_CCB=7.0 V1_CCB=210 Q_CCB=5   V2_CCB=1200
 EC50_CCB=4.5 EMAX_CCB=0.18 EC50CBF_CCB=3.0 EMAX_CBF_CCB=0.30
 KA_RAN=1.2 F_RAN=0.76 CL_RAN=60  V1_RAN=600
 EC50_RAN=300 EMAX_RAN=0.70 EC50EX_RAN=250 EMAX_EX_RAN=0.15
 KA_IVA=2.5 F_IVA=0.40 CL_IVA=80  V1_IVA=230
 EC50_IVA=1.0 EMAX_IVA=0.28
 KA_NIT=2.0 F_NIT=1.0  CL_NIT=40  V1_NIT=260
 EC50_NIT=200 EMAX_NIT=0.20 EC50CBF_NIT=150 EMAX_CBF_NIT=0.20
 TOL_K1=0.05 TOL_K2=0.01
 KOUT_HR=0.3 KOUT_SBP=0.2 KOUT_CBF=0.5 KOUT_ISCH=0.1 KOUT_ANG=0.05
 RPP_THRESH=20000 WDMD=0.6 WSUP=0.4 STEN_K=0.035
 KPROG_PL=0.00014 KSTATIN_PL=0.00008
 STATIN_ON=0

$MAIN
 HR_STATE_0   = HR0;
 SBP_STATE_0  = SBP0;
 CBF_STATE_0  = CBF0 * (1 - STEN_K * STENOSIS);
 ANGINA_SCORE_0 = 5.0 * STENOSIS / 70.0;
 EX_CAP_0     = 7.0 * (1 - 0.03*(STENOSIS-50)/10);
 PLAQUE_0     = STENOSIS;

$ODE
 double KE_BB=CL_BB/V1_BB; double K12_BB=Q_BB/V1_BB; double K21_BB=Q_BB/V2_BB;
 dxdt_GUT_BB     = -KA_BB*GUT_BB;
 dxdt_CENTRAL_BB =  KA_BB*F_BB*GUT_BB-(KE_BB+K12_BB)*CENTRAL_BB+K21_BB*PERIPH_BB;
 dxdt_PERIPH_BB  =  K12_BB*CENTRAL_BB-K21_BB*PERIPH_BB;
 double Cp_BB = CENTRAL_BB/V1_BB*1000.0;

 double KE_CCB=CL_CCB/V1_CCB; double K12_CCB=Q_CCB/V1_CCB; double K21_CCB=Q_CCB/V2_CCB;
 dxdt_GUT_CCB     = -KA_CCB*GUT_CCB;
 dxdt_CENTRAL_CCB =  KA_CCB*F_CCB*GUT_CCB-(KE_CCB+K12_CCB)*CENTRAL_CCB+K21_CCB*PERIPH_CCB;
 dxdt_PERIPH_CCB  =  K12_CCB*CENTRAL_CCB-K21_CCB*PERIPH_CCB;
 double Cp_CCB = CENTRAL_CCB/V1_CCB*1000.0;

 dxdt_GUT_RAN     = -KA_RAN*GUT_RAN;
 dxdt_CENTRAL_RAN =  KA_RAN*F_RAN*GUT_RAN-(CL_RAN/V1_RAN)*CENTRAL_RAN;
 double Cp_RAN = CENTRAL_RAN/V1_RAN*1000.0;

 dxdt_GUT_IVA     = -KA_IVA*GUT_IVA;
 dxdt_CENTRAL_IVA =  KA_IVA*F_IVA*GUT_IVA-(CL_IVA/V1_IVA)*CENTRAL_IVA;
 double Cp_IVA = CENTRAL_IVA/V1_IVA*1000.0;

 dxdt_GUT_NIT     = -KA_NIT*GUT_NIT;
 dxdt_CENTRAL_NIT =  KA_NIT*F_NIT*GUT_NIT-(CL_NIT/V1_NIT)*CENTRAL_NIT;
 double Cp_NIT = CENTRAL_NIT/V1_NIT*1000.0;
 double tol_ss=(Cp_NIT>1.0)?1.0:0.0;
 dxdt_NIT_TOL=TOL_K1*(tol_ss-NIT_TOL)-TOL_K2*NIT_TOL*(1-tol_ss);
 double eff_NIT=(1-0.75*NIT_TOL);

 double E_BB_HR  =EMAX_BB*Cp_BB/(EC50_BB+Cp_BB);
 double E_BB_SBP =EMAX_SBP_BB*Cp_BB/(EC50_SBP_BB+Cp_BB);
 double E_CCB_SBP=EMAX_CCB*Cp_CCB/(EC50_CCB+Cp_CCB);
 double E_CCB_CBF=EMAX_CBF_CCB*Cp_CCB/(EC50CBF_CCB+Cp_CCB);
 double E_IVA_HR =EMAX_IVA*Cp_IVA/(EC50_IVA+Cp_IVA);
 double E_NIT_MVO=EMAX_NIT*Cp_NIT/(EC50_NIT+Cp_NIT)*eff_NIT;
 double E_NIT_CBF=EMAX_CBF_NIT*Cp_NIT/(EC50CBF_NIT+Cp_NIT)*eff_NIT;
 double E_RAN_INa=EMAX_RAN*Cp_RAN/(EC50_RAN+Cp_RAN);
 double E_RAN_EX =EMAX_EX_RAN*Cp_RAN/(EC50EX_RAN+Cp_RAN);

 double HR_tgt  = HR0*(1-E_BB_HR)*(1-E_IVA_HR);
 double SBP_tgt = SBP0*(1-E_BB_SBP)*(1-E_CCB_SBP);
 double CBF_base= CBF0*(1-STEN_K*PLAQUE);
 double CBF_tgt = CBF_base*(1+E_CCB_CBF+E_NIT_CBF);

 dxdt_HR_STATE  = KOUT_HR*(HR_tgt-HR_STATE);
 dxdt_SBP_STATE = KOUT_SBP*(SBP_tgt-SBP_STATE);
 dxdt_CBF_STATE = KOUT_CBF*(CBF_tgt-CBF_STATE);

 double MVO2_mod = (HR_STATE*SBP_STATE)*(1-E_NIT_MVO*0.5);
 double MVO2_eff = MVO2_mod*(1-E_RAN_INa*0.1);
 double O2_supply = CBF_STATE/CBF0;
 double O2_demand = MVO2_eff/(HR0*SBP0);
 dxdt_O2IMBAL = 0.3*((O2_demand*WDMD-O2_supply*WSUP)-O2IMBAL);

 double isch_drive=(MVO2_eff>RPP_THRESH)?(MVO2_eff-RPP_THRESH)/RPP_THRESH:0.0;
 double isch_tgt  =isch_drive*(1-E_RAN_INa*0.3);
 dxdt_ISCHEMIA = KOUT_ISCH*(isch_tgt-ISCHEMIA);

 double ang0    = 5.0*STENOSIS/70.0;
 double ang_tgt = ang0*(ISCHEMIA/0.15+0.001);
 ang_tgt = (ang_tgt>14)?14:ang_tgt;
 dxdt_ANGINA_SCORE = KOUT_ANG*(ang_tgt-ANGINA_SCORE);

 double ex0     = 7.0*(1-0.03*(STENOSIS-50)/10);
 double ex_tgt  = ex0*(1-0.5*isch_drive)*(1+E_RAN_EX);
 ex_tgt = (ex_tgt<1.0)?1.0:ex_tgt;
 dxdt_EX_CAP = 0.05*(ex_tgt-EX_CAP);

 double plaque_prog=KPROG_PL*PLAQUE*(1-PLAQUE/100.0);
 double plaque_regr=STATIN_ON*KSTATIN_PL*PLAQUE;
 dxdt_PLAQUE = plaque_prog-plaque_regr;

$TABLE
 capture HR_sim=HR_STATE; capture SBP_sim=SBP_STATE;
 capture RPP_sim=HR_STATE*SBP_STATE; capture CBF_sim=CBF_STATE;
 capture Isch=ISCHEMIA; capture Angina=ANGINA_SCORE; capture ExCap=EX_CAP;
 capture Plaque=PLAQUE; capture NIT_tol=NIT_TOL;
 capture Cp_BB=Cp_BB; capture Cp_CCB=Cp_CCB;
 capture Cp_RAN=Cp_RAN; capture Cp_IVA=Cp_IVA; capture Cp_NIT=Cp_NIT;
'

mod <- mcode("sa_shiny", model_code, quiet = TRUE)

# ─── Helper: run simulation ───────────────────────────────────────────────────
run_sim <- function(hr0, sbp0, stenosis, statin_on,
                    bb_mg, ccb_mg, ran_mg, iva_mg, nit_mg,
                    t_end = 168) {

  ev_list <- ev(time = 0, amt = 0)

  add_if <- function(cmt, amt, ii, t_end) {
    if (amt > 0)
      ev(cmt = cmt, amt = amt, ii = ii, addl = max(0, floor(t_end/ii) - 1))
    else NULL
  }
  all_ev <- c(
    add_if("GUT_BB",  bb_mg,  24, t_end),
    add_if("GUT_CCB", ccb_mg, 24, t_end),
    add_if("GUT_RAN", ran_mg, 12, t_end),
    add_if("GUT_IVA", iva_mg, 12, t_end),
    add_if("GUT_NIT", nit_mg, 12, t_end)
  )
  final_ev <- if (is.null(all_ev)) ev_list else do.call(c, all_ev)

  t_grid <- seq(0, t_end, by = 0.5)

  mod %>%
    param(HR0 = hr0, SBP0 = sbp0, STENOSIS = stenosis, STATIN_ON = statin_on) %>%
    mrgsim(final_ev, tgrid = t_grid) %>%
    as.data.frame()
}

# ─── CCS Classification ───────────────────────────────────────────────────────
classify_CCS <- function(angina_wk) {
  dplyr::case_when(
    angina_wk == 0       ~ "CCS I  — No ordinary activity limitation",
    angina_wk < 3        ~ "CCS II — Slight limitation of ordinary activity",
    angina_wk < 7        ~ "CCS III — Marked limitation of ordinary activity",
    TRUE                 ~ "CCS IV — Inability to carry on any activity"
  )
}

# ─── UI ───────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "red",

  dashboardHeader(title = "Stable Angina QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "tab_patient",   icon = icon("user")),
      menuItem("Drug PK",            tabName = "tab_pk",        icon = icon("pills")),
      menuItem("Hemodynamics",       tabName = "tab_hemo",      icon = icon("heartbeat")),
      menuItem("Ischemia & Angina",  tabName = "tab_ischemia",  icon = icon("exclamation-triangle")),
      menuItem("Scenario Comparison",tabName = "tab_scenario",  icon = icon("chart-bar")),
      menuItem("Biomarkers & Risk",  tabName = "tab_bio",       icon = icon("flask"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .info-box { min-height: 70px; }
      .info-box .info-box-content { padding: 5px 10px; }
    "))),
    tabItems(

      # ── Tab 1: Patient Profile ──────────────────────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient Characteristics", status = "danger", solidHeader = TRUE, width = 4,
            sliderInput("hr0",       "Baseline Heart Rate (bpm)",    min=50,  max=110, value=80, step=5),
            sliderInput("sbp0",      "Baseline SBP (mmHg)",          min=110, max=180, value=145, step=5),
            sliderInput("stenosis",  "Coronary Stenosis (%)",         min=50,  max=95,  value=70, step=5),
            sliderInput("age",       "Patient Age (years)",           min=40,  max=85,  value=62, step=1),
            checkboxInput("statin",  "Statin Therapy (Rosuvastatin)", value=FALSE)
          ),
          box(title = "Drug Doses", status = "primary", solidHeader = TRUE, width = 4,
            numericInput("bb_mg",  "Bisoprolol (mg QD)",    value=5,    min=0, max=20,   step=1.25),
            numericInput("ccb_mg", "Amlodipine (mg QD)",    value=5,    min=0, max=10,   step=2.5),
            numericInput("ran_mg", "Ranolazine ER (mg BID)",value=0,    min=0, max=1000, step=500),
            numericInput("iva_mg", "Ivabradine (mg BID)",   value=0,    min=0, max=7.5,  step=2.5),
            numericInput("nit_mg", "ISMN Nitrate (mg BID)", value=0,    min=0, max=60,   step=20),
            sliderInput("t_end",   "Simulation Duration (h)", min=24, max=336, value=168, step=24),
            actionButton("run_btn", "Run Simulation", class="btn-success btn-lg", icon=icon("play"))
          ),
          box(title = "Patient Summary", status = "warning", solidHeader = TRUE, width = 4,
            h4("Baseline Status"),
            verbatimTextOutput("patient_summary"),
            hr(),
            h4("Predicted CCS Class"),
            h2(textOutput("ccs_output"), style="color:#C62828; font-weight:bold;")
          )
        ),
        fluidRow(
          valueBoxOutput("vbox_hr",     width=3),
          valueBoxOutput("vbox_sbp",    width=3),
          valueBoxOutput("vbox_rpp",    width=3),
          valueBoxOutput("vbox_angina", width=3)
        )
      ),

      # ── Tab 2: Drug PK ──────────────────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Plasma Concentration Profiles", status = "primary",
              solidHeader = TRUE, width = 12,
            plotlyOutput("pk_plot", height = "450px")
          )
        ),
        fluidRow(
          box(title = "PK Parameters Summary", status = "info", solidHeader=TRUE, width=12,
            DTOutput("pk_table")
          )
        )
      ),

      # ── Tab 3: Hemodynamics ─────────────────────────────────────────────────
      tabItem(tabName = "tab_hemo",
        fluidRow(
          box(title = "Heart Rate & Blood Pressure", status = "danger",
              solidHeader = TRUE, width = 6,
            plotlyOutput("hr_sbp_plot", height="350px")
          ),
          box(title = "Rate-Pressure Product (RPP)", status = "danger",
              solidHeader = TRUE, width = 6,
            plotlyOutput("rpp_plot", height="350px")
          )
        ),
        fluidRow(
          box(title = "Coronary Blood Flow", status = "primary",
              solidHeader = TRUE, width = 6,
            plotlyOutput("cbf_plot", height="350px")
          ),
          box(title = "O₂ Supply-Demand Balance", status = "warning",
              solidHeader = TRUE, width = 6,
            plotlyOutput("o2_plot", height="350px")
          )
        )
      ),

      # ── Tab 4: Ischemia & Angina ────────────────────────────────────────────
      tabItem(tabName = "tab_ischemia",
        fluidRow(
          box(title = "Ischemia Burden Index", status = "danger",
              solidHeader = TRUE, width = 6,
            plotlyOutput("isch_plot", height="350px")
          ),
          box(title = "Angina Episodes (/week)", status = "warning",
              solidHeader = TRUE, width = 6,
            plotlyOutput("angina_plot", height="350px")
          )
        ),
        fluidRow(
          box(title = "Exercise Capacity (METs)", status = "success",
              solidHeader = TRUE, width = 6,
            plotlyOutput("excap_plot", height="350px")
          ),
          box(title = "Nitrate Tolerance State", status = "info",
              solidHeader = TRUE, width = 6,
            plotlyOutput("tol_plot", height="350px")
          )
        )
      ),

      # ── Tab 5: Scenario Comparison ──────────────────────────────────────────
      tabItem(tabName = "tab_scenario",
        fluidRow(
          box(title = "6 Standard Therapy Scenarios", status = "primary",
              solidHeader = TRUE, width = 12,
            actionButton("run_scenarios_btn", "Run All Scenarios", class="btn-primary", icon=icon("sync")),
            hr(),
            plotlyOutput("scenario_plot", height="500px")
          )
        ),
        fluidRow(
          box(title = "Scenario Summary Table (Week 1 Steady State)", status="info",
              solidHeader=TRUE, width=12,
            DTOutput("scenario_table")
          )
        )
      ),

      # ── Tab 6: Biomarkers & Risk ────────────────────────────────────────────
      tabItem(tabName = "tab_bio",
        fluidRow(
          box(title = "Plaque Progression (2-Year)", status = "danger",
              solidHeader = TRUE, width = 6,
            plotlyOutput("plaque_plot", height="380px")
          ),
          box(title = "Bisoprolol Dose-Response", status = "primary",
              solidHeader = TRUE, width = 6,
            plotlyOutput("dr_plot", height="380px")
          )
        ),
        fluidRow(
          box(title = "Cardiovascular Risk Indicators", status="warning",
              solidHeader=TRUE, width=12,
            DTOutput("risk_table")
          )
        )
      )

    )  # end tabItems
  )  # end dashboardBody
)  # end dashboardPage

# ─── SERVER ───────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  sim_data <- eventReactive(input$run_btn, {
    run_sim(
      hr0       = input$hr0,
      sbp0      = input$sbp0,
      stenosis  = input$stenosis,
      statin_on = as.integer(input$statin),
      bb_mg     = input$bb_mg,
      ccb_mg    = input$ccb_mg,
      ran_mg    = input$ran_mg,
      iva_mg    = input$iva_mg,
      nit_mg    = input$nit_mg,
      t_end     = input$t_end
    )
  }, ignoreNULL = FALSE)

  # Value boxes
  ss_stats <- reactive({
    df <- sim_data() %>% filter(time >= max(time)*0.85)
    list(
      hr  = round(mean(df$HR_sim,   na.rm=TRUE), 1),
      sbp = round(mean(df$SBP_sim,  na.rm=TRUE), 1),
      rpp = round(mean(df$RPP_sim,  na.rm=TRUE), 0),
      ang = round(mean(df$Angina,   na.rm=TRUE), 1)
    )
  })

  output$vbox_hr <- renderValueBox({
    valueBox(paste0(ss_stats()$hr, " bpm"), "Heart Rate (SS)",
             icon=icon("heartbeat"), color=if (ss_stats()$hr < 65) "green" else if (ss_stats()$hr < 80) "yellow" else "red")
  })
  output$vbox_sbp <- renderValueBox({
    valueBox(paste0(ss_stats()$sbp, " mmHg"), "Systolic BP (SS)",
             icon=icon("tachometer-alt"), color=if (ss_stats()$sbp < 130) "green" else if (ss_stats()$sbp < 145) "yellow" else "red")
  })
  output$vbox_rpp <- renderValueBox({
    valueBox(format(ss_stats()$rpp, big.mark=","), "RPP (SS)",
             icon=icon("chart-line"), color=if (ss_stats()$rpp < 18000) "green" else if (ss_stats()$rpp < 22000) "yellow" else "red")
  })
  output$vbox_angina <- renderValueBox({
    valueBox(paste0(ss_stats()$ang, "/wk"), "Angina Episodes (SS)",
             icon=icon("exclamation-circle"), color=if (ss_stats()$ang < 2) "green" else if (ss_stats()$ang < 5) "yellow" else "red")
  })

  output$patient_summary <- renderText({
    paste0(
      "Age: ", input$age, " years\n",
      "Baseline HR: ", input$hr0, " bpm\n",
      "Baseline SBP: ", input$sbp0, " mmHg\n",
      "Stenosis: ", input$stenosis, "%\n",
      "Baseline RPP: ", input$hr0 * input$sbp0, "\n",
      "Statin: ", ifelse(input$statin, "Yes", "No")
    )
  })

  output$ccs_output <- renderText({
    classify_CCS(ss_stats()$ang)
  })

  # PK plot
  output$pk_plot <- renderPlotly({
    df <- sim_data() %>%
      select(time, Cp_BB, Cp_CCB, Cp_RAN, Cp_IVA, Cp_NIT) %>%
      pivot_longer(-time, names_to="Drug", values_to="Cp") %>%
      mutate(Drug = recode(Drug,
        Cp_BB="Bisoprolol", Cp_CCB="Amlodipine",
        Cp_RAN="Ranolazine", Cp_IVA="Ivabradine", Cp_NIT="ISMN Nitrate"))

    p <- ggplot(df, aes(time, Cp, color=Drug, group=Drug)) +
      geom_line(linewidth=0.8) +
      labs(title="Plasma Concentration vs. Time",
           x="Time (h)", y="Cp (ng/mL)") +
      theme_bw(base_size=11) +
      scale_color_manual(values=c("#7B1FA2","#1565C0","#2E7D32","#E65100","#C62828"))
    ggplotly(p)
  })

  output$pk_table <- renderDT({
    tbl <- data.frame(
      Drug = c("Bisoprolol","Amlodipine","Ranolazine ER","Ivabradine","ISMN"),
      Dose_Regimen = c(paste0(input$bb_mg," mg QD"), paste0(input$ccb_mg," mg QD"),
                       paste0(input$ran_mg," mg BID"), paste0(input$iva_mg," mg BID"),
                       paste0(input$nit_mg," mg BID")),
      BA_pct = c("85%","64%","76%","40%","100%"),
      t_half = c("11h","35–50h","7h","2h","4–5h"),
      Vd_L   = c("~600","~1400","~600","~230","~260"),
      Metabolism = c("CYP2D6","CYP3A4","CYP3A4","CYP3A4","Hepatic"),
      Mechanism = c("β₁ blockade","L-type Ca²⁺ block","Late I_Na inhib","HCN/I_f block","NO→cGMP")
    )
    datatable(tbl, rownames=FALSE, options=list(pageLength=7, dom="t"),
              class="stripe compact")
  }, server=FALSE)

  # Hemodynamics
  output$hr_sbp_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df) +
      geom_line(aes(time, HR_sim,  color="HR (bpm)"),  linewidth=0.9) +
      geom_line(aes(time, SBP_sim, color="SBP (mmHg)"),linewidth=0.9) +
      geom_hline(yintercept=60, linetype="dashed", color="#9E9E9E") +
      scale_color_manual(values=c("HR (bpm)"="#E91E63","SBP (mmHg)"="#1565C0")) +
      labs(x="Time (h)", y="Value", color="Parameter") + theme_bw(base_size=11)
    ggplotly(p)
  })

  output$rpp_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time, RPP_sim/1000)) +
      geom_line(color="#AD1457", linewidth=0.9) +
      geom_hline(yintercept=20, linetype="dashed", color="red") +
      annotate("text", x=max(df$time)*0.6, y=20.8, label="Ischemic Threshold", size=3, color="red") +
      labs(x="Time (h)", y="RPP (×1000)", title="Rate-Pressure Product") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$cbf_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time, CBF_sim)) +
      geom_line(color="#1565C0", linewidth=0.9) +
      geom_hline(yintercept=100*(1-0.035*input$stenosis), linetype="dashed", color="grey60") +
      labs(x="Time (h)", y="CBF (mL/min)") + theme_bw(base_size=11)
    ggplotly(p)
  })

  output$o2_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time, O2IMBAL)) +
      geom_line(color="#F57F17", linewidth=0.9) +
      geom_hline(yintercept=0, linetype="solid", color="grey50") +
      annotate("text", x=max(df$time)*0.1, y=max(df$O2IMBAL)*0.9,
               label="+ = demand excess", size=3, color="#B71C1C") +
      labs(x="Time (h)", y="O₂ Imbalance (norm.)") + theme_bw(base_size=11)
    ggplotly(p)
  })

  # Ischemia & Angina
  output$isch_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time, Isch)) +
      geom_area(fill="#EF9A9A", alpha=0.4) + geom_line(color="#B71C1C", linewidth=0.9) +
      labs(x="Time (h)", y="Ischemia Index [0–1]") + theme_bw(base_size=11)
    ggplotly(p)
  })

  output$angina_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time, Angina)) +
      geom_line(color="#E65100", linewidth=0.9) +
      geom_hline(yintercept=c(0,3,7), linetype="dashed", color=c("green","orange","red")) +
      annotate("text", x=max(df$time)*0.7, y=c(0.5,3.5,7.5),
               label=c("CCS I","CCS II/III","CCS IV"), size=3,
               color=c("green","orange","red")) +
      labs(x="Time (h)", y="Angina Episodes/week") + theme_bw(base_size=11)
    ggplotly(p)
  })

  output$excap_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time, ExCap)) +
      geom_line(color="#2E7D32", linewidth=0.9) +
      geom_hline(yintercept=c(4,7), linetype="dashed", color=c("orange","green")) +
      labs(x="Time (h)", y="Exercise Capacity (METs)") + theme_bw(base_size=11)
    ggplotly(p)
  })

  output$tol_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time, NIT_tol)) +
      geom_line(color="#827717", linewidth=0.9) +
      labs(x="Time (h)", y="Nitrate Tolerance [0–1]",
           subtitle="Eccentric dosing (BID) prevents full tolerance") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  # Scenario comparison
  scenario_data <- eventReactive(input$run_scenarios_btn, {
    scenarios <- list(
      "S1: Untreated"             = list(bb=0, ccb=0, ran=0,   iva=0, nit=0,  stat=0),
      "S2: Bisoprolol 5mg"        = list(bb=5, ccb=0, ran=0,   iva=0, nit=0,  stat=0),
      "S3: BB + Amlodipine"       = list(bb=5, ccb=5, ran=0,   iva=0, nit=0,  stat=0),
      "S4: BB + CCB + Ranolazine" = list(bb=5, ccb=5, ran=1000,iva=0, nit=0,  stat=0),
      "S5: BB + Ivabradine"       = list(bb=5, ccb=0, ran=0,   iva=5, nit=0,  stat=0),
      "S6: BB + CCB + ISMN"       = list(bb=5, ccb=5, ran=0,   iva=0, nit=40, stat=1)
    )
    t_sim <- seq(0, 168, by=0.5)
    bind_rows(lapply(names(scenarios), function(nm) {
      s <- scenarios[[nm]]
      all_ev <- c(
        if (s$bb>0)  ev(cmt="GUT_BB",  amt=s$bb,  ii=24, addl=6) else NULL,
        if (s$ccb>0) ev(cmt="GUT_CCB", amt=s$ccb, ii=24, addl=6) else NULL,
        if (s$ran>0) ev(cmt="GUT_RAN", amt=s$ran, ii=12, addl=13) else NULL,
        if (s$iva>0) ev(cmt="GUT_IVA", amt=s$iva, ii=12, addl=13) else NULL,
        if (s$nit>0) ev(cmt="GUT_NIT", amt=s$nit, ii=12, addl=13) else NULL
      )
      final_ev <- if (is.null(all_ev)) ev(time=0, amt=0) else do.call(c, all_ev)
      mod %>%
        param(HR0=input$hr0, SBP0=input$sbp0, STENOSIS=input$stenosis, STATIN_ON=s$stat) %>%
        mrgsim(final_ev, tgrid=t_sim) %>%
        as.data.frame() %>%
        mutate(Scenario=nm)
    }))
  }, ignoreNULL=FALSE)

  output$scenario_plot <- renderPlotly({
    df <- scenario_data()
    cols <- c("#9E9E9E","#7B1FA2","#1565C0","#2E7D32","#E65100","#C62828")
    p <- ggplot(df, aes(time, RPP_sim/1000, color=Scenario)) +
      geom_line(linewidth=0.8) +
      scale_color_manual(values=cols) +
      geom_hline(yintercept=20, linetype="dashed", color="red") +
      labs(title="Rate-Pressure Product — All Scenarios",
           x="Time (h)", y="RPP (×1000)", color="Scenario") +
      theme_bw(base_size=11) + theme(legend.position="bottom")
    ggplotly(p)
  })

  output$scenario_table <- renderDT({
    df <- scenario_data() %>%
      filter(time >= 144) %>%
      group_by(Scenario) %>%
      summarise(
        `HR (bpm)`       = round(mean(HR_sim),  1),
        `SBP (mmHg)`     = round(mean(SBP_sim), 1),
        `RPP`            = round(mean(RPP_sim),  0),
        `CBF (mL/min)`   = round(mean(CBF_sim),  1),
        `Ischemia Index` = round(mean(Isch),      3),
        `Angina/wk`      = round(mean(Angina),    1),
        `ExCap (METs)`   = round(mean(ExCap),     1),
        .groups="drop"
      ) %>%
      mutate(CCS=classify_CCS(`Angina/wk`))
    datatable(df, rownames=FALSE, options=list(pageLength=8, dom="t"),
              class="stripe compact") %>%
      formatStyle(
        "Ischemia Index",
        backgroundColor=styleInterval(c(0.05,0.15), c("#C8E6C9","#FFF9C4","#FFCDD2"))
      )
  }, server=FALSE)

  # Biomarkers & Risk
  output$plaque_plot <- renderPlotly({
    t_long <- seq(0, 365*24*2, by=24)
    ev_base <- ev(time=0, amt=0)
    d_no  <- mod %>% param(STENOSIS=input$stenosis, STATIN_ON=0) %>%
      mrgsim(ev_base, tgrid=t_long) %>% as.data.frame() %>% mutate(Group="No Statin")
    d_yes <- mod %>% param(STENOSIS=input$stenosis, STATIN_ON=1) %>%
      mrgsim(ev_base, tgrid=t_long) %>% as.data.frame() %>% mutate(Group="Statin")
    df <- bind_rows(d_no, d_yes) %>% mutate(time_yr=time/(24*365))
    p <- ggplot(df, aes(time_yr, Plaque, color=Group)) +
      geom_line(linewidth=1) +
      scale_color_manual(values=c("#C62828","#2E7D32")) +
      labs(x="Years", y="Stenosis (%)", title="Long-term Plaque Trajectory") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$dr_plot <- renderPlotly({
    bb_doses <- c(1.25, 2.5, 5, 10, 20)
    t_grid   <- seq(0, 168, by=0.5)
    dr <- bind_rows(lapply(bb_doses, function(d) {
      ev_d <- ev(cmt="GUT_BB", amt=d, ii=24, addl=6)
      mod %>%
        param(HR0=input$hr0, SBP0=input$sbp0, STENOSIS=input$stenosis) %>%
        mrgsim(ev_d, tgrid=t_grid) %>%
        as.data.frame() %>%
        filter(time >= 144) %>%
        summarise(Dose_mg=d, HR_ss=mean(HR_sim), RPP_ss=mean(RPP_sim))
    }))
    p <- ggplot(dr, aes(Dose_mg, RPP_ss/1000)) +
      geom_line(color="#7B1FA2", linewidth=1) +
      geom_point(color="#7B1FA2", size=3) +
      geom_hline(yintercept=20, linetype="dashed", color="red") +
      annotate("text", x=8, y=20.5, label="Ischemic Threshold", size=3, color="red") +
      scale_x_continuous(breaks=bb_doses) +
      labs(title="Bisoprolol Dose-Response (RPP)", x="Dose (mg QD)", y="RPP ×1000 (SS)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$risk_table <- renderDT({
    s <- ss_stats()
    rpp_val <- s$rpp
    tbl <- data.frame(
      Parameter      = c("Heart Rate","Systolic BP","Rate-Pressure Product","Angina/week","Exercise Capacity"),
      Value          = c(paste0(s$hr," bpm"), paste0(s$sbp," mmHg"),
                         format(s$rpp,big.mark=","), paste0(s$ang,"/week"),
                         paste0(round(sim_data() %>% filter(time>=max(time)*0.85) %>% pull(ExCap) %>% mean(),1)," METs")),
      Target         = c("55–65 bpm","< 130 mmHg","< 20,000","< 2/week","> 7 METs"),
      Status         = c(
        ifelse(s$hr<65,"✓ On target","⚠ Above target"),
        ifelse(s$sbp<130,"✓ On target","⚠ Above target"),
        ifelse(rpp_val<20000,"✓ Below threshold","⚠ Above threshold"),
        ifelse(s$ang<2,"✓ Well controlled","⚠ Suboptimal control"),
        "—"
      )
    )
    datatable(tbl, rownames=FALSE, options=list(dom="t", pageLength=7),
              class="stripe compact") %>%
      formatStyle("Status",
        color=styleEqual(c("✓ On target","✓ Below threshold","✓ Well controlled","—"),
                         c("green","green","green","grey")),
        fontWeight="bold")
  }, server=FALSE)

}

shinyApp(ui, server)
