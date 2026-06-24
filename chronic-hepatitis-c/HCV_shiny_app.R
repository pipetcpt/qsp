## ============================================================
## Chronic Hepatitis C QSP — Interactive Shiny Dashboard
## Version: 1.0  |  2026-06-23
## Tabs:
##   1. Patient Profile & Disease Setting
##   2. Drug PK Profiles (DAA & PEG-IFN/RBV)
##   3. PD Biomarkers (ALT, Viral Load kinetics)
##   4. Clinical Endpoints (SVR, RVR, TND timeline)
##   5. Treatment Scenario Comparison (all regimens)
##   6. Immune Landscape (CTL exhaustion, NK, Treg)
## ============================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(shinydashboard)
library(DT)

## ── mrgsolve model (minimal version for Shiny) ─────────────
hcv_code <- '
$PARAM
ka_SOF=1.44 CL_SOF_Tp=0.05 V_SOF_Tp=1.0
EC50_SOF=0.35 n_SOF=1.6
CL_LED=0.25 EC50_LED=0.031 n_LED=1.2
CL_VEL=0.33 EC50_VEL=0.017
CL_GLE=1.2  EC50_GLE=0.45  n_GLE=1.8
CL_PIB=0.18 EC50_PIB=0.002
CL_RBV=0.1  EC50_RBV=2.0
k_RBV_RBC=0.02 kdeg_RBV_RBC=0.003
CL_PEGIFN=0.025 EC50_PEGIFN=0.8
beta_infect=1.2e-7 p_prod=100 c_clear=22
delta_I=0.08 d_T=0.001 s_T=10000
kprol_CTL=0.02 kdeath_CTL=0.004 keff_CTL=1e-4
kprol_NK=0.01  kdeath_NK=0.003  keff_NK=5e-5
kprol_Treg=0.005 kdeath_Treg=0.003 k_exhaust=1e-5
kprog_fibro=1e-5 krem_fibro=2e-6
kprog_HSC=0.001 kdeath_HSC=0.0005
kALT_prod=0.5 kALT_clear=0.008
kHCC_prog=1e-6

$CMT SOF_Tp LED_p VEL_p NS5A_i GLE_p PIB_p RBV_p RBV_RBC PEGIFN_p
     T_cell I_cell V_rna V_def
     CTL NK_cell Treg_HCV
     ALT Fibro_met HSC_act HCC_idx

$GLOBAL
double E_SOF_p,E_LED_p,E_GLE_p,E_PIB_p,E_RBV_p,E_PEGIFN_p;
double Ep,Ei;
double lytic_CTL,lytic_NK;

$MAIN
E_SOF_p=pow(SOF_Tp,n_SOF)/(pow(EC50_SOF,n_SOF)+pow(SOF_Tp,n_SOF));
E_LED_p=LED_p/(EC50_LED+LED_p);
E_GLE_p=pow(GLE_p,n_GLE)/(pow(EC50_GLE,n_GLE)+pow(GLE_p,n_GLE));
E_PIB_p=PIB_p/(EC50_PIB+PIB_p);
E_RBV_p=RBV_p/(EC50_RBV+RBV_p);
E_PEGIFN_p=PEGIFN_p/(EC50_PEGIFN+PEGIFN_p);
double E_NS5A=(NS5A_i>0)?NS5A_i/(EC50_VEL+NS5A_i):E_LED_p;
double prod_block=(1.0-E_SOF_p)*(1.0-E_NS5A)*(1.0-E_GLE_p)*(1.0-E_PIB_p*0.5);
Ep=1.0-prod_block;
Ei=E_SOF_p*0.4+E_RBV_p*0.3;
if(Ei>0.99) Ei=0.99;
lytic_CTL=keff_CTL*CTL*I_cell;
lytic_NK=keff_NK*NK_cell*I_cell;

$ODE
dxdt_SOF_Tp=-CL_SOF_Tp*SOF_Tp;
dxdt_LED_p=-CL_LED*LED_p;
dxdt_VEL_p=-CL_VEL*VEL_p;
dxdt_NS5A_i=-0.3*NS5A_i;
dxdt_GLE_p=-CL_GLE*GLE_p;
dxdt_PIB_p=-CL_PIB*PIB_p;
dxdt_RBV_p=-CL_RBV*RBV_p-k_RBV_RBC*RBV_p;
dxdt_RBV_RBC=k_RBV_RBC*RBV_p-kdeg_RBV_RBC*RBV_RBC;
dxdt_PEGIFN_p=-CL_PEGIFN*PEGIFN_p;
dxdt_T_cell=s_T-d_T*T_cell-(1.0-Ei)*beta_infect*V_rna*T_cell;
if(T_cell<0) dxdt_T_cell=0;
double inf_rate=(1.0-Ei)*beta_infect*V_rna*T_cell;
dxdt_I_cell=inf_rate-delta_I*I_cell-lytic_CTL-lytic_NK;
if(I_cell<0) dxdt_I_cell=0;
dxdt_V_rna=(1.0-Ep)*p_prod*I_cell-c_clear*V_rna;
if(V_rna<0) dxdt_V_rna=0;
dxdt_V_def=0.05*(1.0-Ep)*p_prod*I_cell-c_clear*V_def;
if(V_def<0) dxdt_V_def=0;
double CTL_s=kprol_CTL*I_cell;
dxdt_CTL=CTL_s-kdeath_CTL*CTL-(Treg_HCV/(Treg_HCV+100.0))*CTL*0.5-k_exhaust*CTL*I_cell;
if(CTL<0) dxdt_CTL=0;
dxdt_NK_cell=kprol_NK*(I_cell/(I_cell+1e4))*(1.0+E_PEGIFN_p)-kdeath_NK*NK_cell;
if(NK_cell<0) dxdt_NK_cell=0;
dxdt_Treg_HCV=kprol_Treg*(I_cell/(I_cell+5e4))-kdeath_Treg*Treg_HCV;
if(Treg_HCV<0) dxdt_Treg_HCV=0;
double cell_death=delta_I*I_cell+lytic_CTL+lytic_NK;
dxdt_ALT=kALT_prod*cell_death-kALT_clear*ALT;
if(ALT<0) dxdt_ALT=0;
dxdt_HSC_act=kprog_HSC*cell_death/(cell_death+1e3)-kdeath_HSC*HSC_act;
if(HSC_act<0) dxdt_HSC_act=0; if(HSC_act>100) dxdt_HSC_act=0;
double svr_f=(V_rna<10.0)?1.0:0.0;
dxdt_Fibro_met=kprog_fibro*HSC_act-krem_fibro*Fibro_met*svr_f;
if(Fibro_met<0) dxdt_Fibro_met=0; if(Fibro_met>4) dxdt_Fibro_met=0;
dxdt_HCC_idx=kHCC_prog*Fibro_met*(1.0-svr_f*0.75);
if(HCC_idx>1) dxdt_HCC_idx=0;

$TABLE
capture log10_V=(V_rna>0)?log10(V_rna):-1.0;
capture ALT_out=ALT;
capture Fibro_out=Fibro_met;
capture CTL_out=CTL;
capture NK_out=NK_cell;
capture Treg_out=Treg_HCV;
capture Ep_pct=Ep*100;
capture Ei_pct=Ei*100;
capture E_SOF_out=E_SOF_p*100;
capture E_LED_out=E_LED_p*100;
capture E_GLE_out=E_GLE_p*100;
capture TND=(V_rna<15)?1.0:0.0;
capture SVR12=(V_rna<15)?1.0:0.0;
capture HCC_risk=HCC_idx;

$INIT
SOF_Tp=0 LED_p=0 VEL_p=0 NS5A_i=0 GLE_p=0 PIB_p=0
RBV_p=0 RBV_RBC=0 PEGIFN_p=0
T_cell=2e7 I_cell=1e5 V_rna=1e6 V_def=1e3
CTL=1000 NK_cell=500 Treg_HCV=200
ALT=80 Fibro_met=1.0 HSC_act=20 HCC_idx=0.01
'

mod <- mcode("hcv_shiny", hcv_code, quiet = TRUE)

## ── Build dosing events ──────────────────────────────────────
build_hcv_events <- function(regimen, duration_wk) {
  dur_h <- duration_wk * 7 * 24
  sof <- ev(amt = 8.0,  cmt = "SOF_Tp", time = seq(0, dur_h - 24, by = 24))
  led <- ev(amt = 120,  cmt = "LED_p",  time = seq(0, dur_h - 24, by = 24))
  vel <- ev(amt = 100,  cmt = "VEL_p",  time = seq(0, dur_h - 24, by = 24))
  gle <- ev(amt = 18,   cmt = "GLE_p",  time = seq(0, dur_h - 24, by = 24))
  pib <- ev(amt = 0.15, cmt = "PIB_p",  time = seq(0, dur_h - 24, by = 24))
  rbv <- ev(amt = 3.0,  cmt = "RBV_p",  time = seq(0, dur_h - 12, by = 12))
  peg <- ev(amt = 150,  cmt = "PEGIFN_p", time = seq(0, dur_h - 168, by = 168))

  switch(regimen,
    "SOF/LED (Harvoni)"   = sof + led,
    "SOF/VEL (Epclusa)"   = sof + vel,
    "GLE/PIB (Mavyret)"   = gle + pib,
    "SOF/VEL+RBV (GT3)"   = sof + vel + rbv,
    "PEG-IFN + RBV"       = peg + rbv,
    ev()  # untreated
  )
}

run_hcv_sim <- function(regimen, duration_wk, fibrosis_base = 1.0, vl_base_log = 6) {
  evts <- build_hcv_events(regimen, duration_wk)
  followup_h <- duration_wk * 7 * 24 + 12 * 7 * 24  # treatment + 12-wk follow-up
  init_list <- list(Fibro_met = fibrosis_base,
                    V_rna = 10^vl_base_log,
                    I_cell = 10^vl_base_log / 100)
  mod %>%
    init(init_list) %>%
    ev(evts) %>%
    mrgsim(end = followup_h, delta = 24) %>%
    as.data.frame() %>%
    mutate(time_wk = time / (7 * 24))
}

## ── UI ──────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "Hepatitis C QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("① Patient Profile",        tabName = "patient",   icon = icon("user-circle")),
      menuItem("② Drug PK Profiles",        tabName = "pk",        icon = icon("pills")),
      menuItem("③ PD Biomarkers",           tabName = "pd",        icon = icon("vial")),
      menuItem("④ Clinical Endpoints",      tabName = "endpoints", icon = icon("chart-bar")),
      menuItem("⑤ Scenario Comparison",     tabName = "compare",   icon = icon("balance-scale")),
      menuItem("⑥ Immune Landscape",        tabName = "immune",    icon = icon("dna"))
    )
  ),
  dashboardBody(
    tabItems(

      ## ── TAB 1 ─────────────────────────────────────────────
      tabItem(tabName = "patient",
        fluidRow(
          box(title = "Patient & Treatment Settings", status = "primary",
              solidHeader = TRUE, width = 4,
            selectInput("regimen",     "DAA Regimen:",
              choices  = c("SOF/LED (Harvoni)", "SOF/VEL (Epclusa)",
                           "GLE/PIB (Mavyret)", "SOF/VEL+RBV (GT3)",
                           "PEG-IFN + RBV",    "Untreated"),
              selected = "SOF/LED (Harvoni)"),
            sliderInput("duration_wk", "Treatment Duration (weeks):",
              min = 8, max = 48, value = 12, step = 4),
            sliderInput("vl_base",     "Baseline Viral Load (log10 IU/mL):",
              min = 4, max = 8, value = 6, step = 0.5),
            sliderInput("fibro_base",  "Baseline Fibrosis (Metavir F0–F4):",
              min = 0, max = 4, value = 1, step = 0.5),
            actionButton("run_btn", "Run Simulation", class = "btn-primary btn-lg")
          ),
          box(title = "Patient Characteristics Summary", status = "warning",
              solidHeader = TRUE, width = 4,
            verbatimTextOutput("pt_summary")
          ),
          box(title = "Week-12 & EOT Outcomes", status = "success",
              solidHeader = TRUE, width = 4,
            h4("Virological Milestones"),
            tableOutput("svr_table"),
            hr(),
            h4("Liver Health"),
            tableOutput("liver_table")
          )
        )
      ),

      ## ── TAB 2 ─────────────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "DAA Intracellular / Plasma Concentrations", status = "primary",
              solidHeader = TRUE, width = 12,
            plotOutput("pk_plot", height = "380px")
          )
        ),
        fluidRow(
          box(title = "Drug Efficacy on Viral Production (εp)", status = "info",
              solidHeader = TRUE, width = 6,
            plotOutput("ep_plot", height = "280px")
          ),
          box(title = "Drug Efficacy on Infectivity (εi)", status = "warning",
              solidHeader = TRUE, width = 6,
            plotOutput("ei_plot", height = "280px")
          )
        )
      ),

      ## ── TAB 3 ─────────────────────────────────────────────
      tabItem(tabName = "pd",
        fluidRow(
          box(title = "Viral Load (log10 IU/mL)", status = "danger",
              solidHeader = TRUE, width = 12,
            plotOutput("vl_plot", height = "350px")
          )
        ),
        fluidRow(
          box(title = "ALT Level", status = "warning",
              solidHeader = TRUE, width = 6,
            plotOutput("alt_plot", height = "280px")
          ),
          box(title = "Infected Hepatocytes (log10)", status = "danger",
              solidHeader = TRUE, width = 6,
            plotOutput("icell_plot", height = "280px")
          )
        )
      ),

      ## ── TAB 4 ─────────────────────────────────────────────
      tabItem(tabName = "endpoints",
        fluidRow(
          box(title = "Virological Response Timeline", status = "primary",
              solidHeader = TRUE, width = 12,
            plotOutput("response_plot", height = "350px")
          )
        ),
        fluidRow(
          box(title = "Liver Fibrosis Score (Metavir)", status = "success",
              solidHeader = TRUE, width = 6,
            plotOutput("fibro_plot", height = "280px")
          ),
          box(title = "HCC Risk Accumulation Index", status = "danger",
              solidHeader = TRUE, width = 6,
            plotOutput("hcc_plot", height = "280px")
          )
        )
      ),

      ## ── TAB 5 ─────────────────────────────────────────────
      tabItem(tabName = "compare",
        fluidRow(
          box(title = "All Regimens — Viral Load Comparison", status = "primary",
              solidHeader = TRUE, width = 12,
            plotOutput("compare_vl_plot", height = "400px")
          )
        ),
        fluidRow(
          box(title = "End-of-Follow-up Summary Table", status = "success",
              solidHeader = TRUE, width = 12,
            DTOutput("compare_table")
          )
        )
      ),

      ## ── TAB 6 ─────────────────────────────────────────────
      tabItem(tabName = "immune",
        fluidRow(
          box(title = "HCV-specific CTL Dynamics", status = "primary",
              solidHeader = TRUE, width = 12,
            plotOutput("ctl_plot", height = "320px")
          )
        ),
        fluidRow(
          box(title = "NK Cells & Treg", status = "warning",
              solidHeader = TRUE, width = 6,
            plotOutput("nk_treg_plot", height = "280px")
          ),
          box(title = "Immune Balance (CTL/Treg Ratio)", status = "danger",
              solidHeader = TRUE, width = 6,
            plotOutput("immune_ratio_plot", height = "280px")
          )
        )
      )
    )
  )
)

## ── SERVER ──────────────────────────────────────────────────
server <- function(input, output, session) {

  sim_data <- eventReactive(input$run_btn, {
    withProgress(message = "Running HCV QSP simulation...", {
      run_hcv_sim(
        regimen     = input$regimen,
        duration_wk = input$duration_wk,
        fibrosis_base = input$fibro_base,
        vl_base_log   = input$vl_base
      )
    })
  }, ignoreNULL = FALSE)

  all_scen_data <- reactive({
    regimens <- c("SOF/LED (Harvoni)", "SOF/VEL (Epclusa)",
                  "GLE/PIB (Mavyret)", "PEG-IFN + RBV", "Untreated")
    purrr::map_dfr(regimens, function(r) {
      run_hcv_sim(r, ifelse(r == "PEG-IFN + RBV", 48, 12),
                  input$fibro_base, input$vl_base) %>%
        mutate(scenario = r)
    })
  })

  output$pt_summary <- renderText({
    paste0(
      "Regimen     : ", input$regimen,     "\n",
      "Duration    : ", input$duration_wk, " weeks\n",
      "Baseline VL : 10^", input$vl_base, " = ",
      format(round(10^input$vl_base), big.mark = ","), " IU/mL\n",
      "Fibrosis    : F", input$fibro_base, " (Metavir)\n",
      "Expected SVR: ",
      ifelse(input$regimen %in% c("SOF/LED (Harvoni)",
             "SOF/VEL (Epclusa)", "GLE/PIB (Mavyret)"),
             "~95–99%", "~50–80%")
    )
  })

  output$svr_table <- renderTable({
    d <- sim_data()
    wk4  <- d %>% filter(abs(time_wk - 4)  < 0.5) %>% slice(1)
    wk12 <- d %>% filter(abs(time_wk - 12) < 0.5) %>% slice(1)
    eot  <- d %>% filter(time_wk == input$duration_wk) %>% slice(1)
    data.frame(
      Milestone = c("Week 4 (RVR)", "Week 12 (eRVR)", "End of Treatment", "SVR12"),
      `VL (log10)` = round(c(wk4$log10_V, wk12$log10_V, eot$log10_V, tail(d$log10_V, 1)), 2),
      `TND` = c(wk4$TND, wk12$TND, eot$TND, tail(d$TND, 1))
    )
  })

  output$liver_table <- renderTable({
    d <- sim_data()
    last <- tail(d, 1)
    data.frame(
      Measure   = c("ALT (IU/L)", "Fibrosis (F)", "HCC Risk Index"),
      `End FU` = round(c(last$ALT_out, last$Fibro_out, last$HCC_risk), 3)
    )
  })

  ## Tab 2: PK
  output$pk_plot <- renderPlot({
    d <- sim_data() %>%
      select(time_wk, E_SOF_out, E_LED_out, E_GLE_out) %>%
      pivot_longer(-time_wk, names_to = "drug", values_to = "efficacy") %>%
      mutate(drug = recode(drug, E_SOF_out="SOF (NS5B)",
                           E_LED_out="LED/VEL (NS5A)", E_GLE_out="GLE (NS3/4A)"))
    ggplot(d, aes(time_wk, efficacy, colour = drug)) +
      geom_line(linewidth = 0.9) +
      ylim(0, 100) +
      labs(x="Time (weeks)", y="Enzyme Inhibition (%)", colour=NULL,
           title="DAA Target Inhibition Over Time") +
      theme_bw(base_size = 12) + theme(legend.position="bottom")
  })

  output$ep_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(time_wk, Ep_pct)) +
      geom_line(colour="#9B59B6", linewidth=1) +
      ylim(0, 100) +
      labs(x="Time (weeks)", y="εp (%)", title="Production Blockade (εp)") +
      theme_bw(base_size=12)
  })

  output$ei_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(time_wk, Ei_pct)) +
      geom_line(colour="#E67E22", linewidth=1) +
      ylim(0, 100) +
      labs(x="Time (weeks)", y="εi (%)", title="Infectivity Blockade (εi)") +
      theme_bw(base_size=12)
  })

  ## Tab 3: PD
  output$vl_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(time_wk, log10_V)) +
      geom_line(colour="#E74C3C", linewidth=1.2) +
      geom_hline(yintercept=log10(15), linetype="dashed") +
      annotate("text", x=0.5, y=log10(15)+0.15,
               label="LLOQ 15 IU/mL", hjust=0, size=3) +
      labs(x="Time (weeks)", y="HCV-RNA (log10 IU/mL)",
           title=paste("Viral Load —", input$regimen)) +
      theme_bw(base_size=12)
  })

  output$alt_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(time_wk, ALT_out)) +
      geom_line(colour="#E67E22", linewidth=1) +
      geom_hline(yintercept=40, linetype="dashed") +
      labs(x="Time (weeks)", y="ALT (IU/L)", title="Serum ALT") +
      theme_bw(base_size=12)
  })

  output$icell_plot <- renderPlot({
    d <- sim_data() %>%
      mutate(log10_I = ifelse(I_cell > 0, log10(I_cell+1), 0))
    ggplot(d, aes(time_wk, log10_I)) +
      geom_line(colour="#C0392B", linewidth=1) +
      labs(x="Time (weeks)", y="Infected cells (log10)",
           title="Infected Hepatocytes") +
      theme_bw(base_size=12)
  })

  ## Tab 4: Endpoints
  output$response_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(time_wk, log10_V)) +
      geom_line(colour="#9B59B6", linewidth=1.2) +
      geom_hline(yintercept=log10(15), linetype="dashed") +
      geom_vline(xintercept=input$duration_wk, linetype="dotted", colour="grey40") +
      annotate("text", x=input$duration_wk+0.3, y=5,
               label="EOT", hjust=0, size=3.5) +
      labs(x="Time (weeks)", y="HCV-RNA (log10 IU/mL)",
           title="Virological Response Timeline") +
      theme_bw(base_size=12)
  })

  output$fibro_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(time_wk, Fibro_out)) +
      geom_line(colour="#27AE60", linewidth=1) +
      scale_y_continuous(limits=c(0,4.2), breaks=0:4) +
      labs(x="Time (weeks)", y="Score (F0–F4)",
           title="Liver Fibrosis — Metavir Score") +
      theme_bw(base_size=12)
  })

  output$hcc_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(time_wk, HCC_risk)) +
      geom_line(colour="#C0392B", linewidth=1) +
      labs(x="Time (weeks)", y="HCC Risk Index (0–1)",
           title="HCC Risk Accumulation") +
      theme_bw(base_size=12)
  })

  ## Tab 5: Comparison
  cols5 <- c(
    "SOF/LED (Harvoni)"   = "#E74C3C",
    "SOF/VEL (Epclusa)"   = "#3498DB",
    "GLE/PIB (Mavyret)"   = "#2ECC71",
    "PEG-IFN + RBV"       = "#9B59B6",
    "Untreated"           = "#7F8C8D"
  )

  output$compare_vl_plot <- renderPlot({
    d <- all_scen_data()
    ggplot(d, aes(time_wk, log10_V, colour=scenario)) +
      geom_line(linewidth=0.9) +
      geom_hline(yintercept=log10(15), linetype="dashed") +
      scale_colour_manual(values=cols5) +
      labs(x="Time (weeks)", y="HCV-RNA (log10 IU/mL)",
           title="All Regimens — Viral Load Comparison", colour=NULL) +
      theme_bw(base_size=12) + theme(legend.position="bottom")
  })

  output$compare_table <- renderDT({
    d <- all_scen_data() %>%
      group_by(scenario) %>%
      filter(time == max(time)) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(across(c(log10_V, ALT_out, Fibro_out, HCC_risk), ~round(.x, 2))) %>%
      select(scenario, log10_V, ALT_out, Fibro_out, SVR12, HCC_risk) %>%
      rename(Scenario="scenario", `log10-VL`="log10_V",
             `ALT (IU/L)`="ALT_out", `Fibrosis (F)`="Fibro_out",
             `SVR12`="SVR12", `HCC Risk`="HCC_risk")
    datatable(d, options=list(dom="t", pageLength=10))
  })

  ## Tab 6: Immune
  output$ctl_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(time_wk, CTL_out)) +
      geom_line(colour="#1A5276", linewidth=1) +
      labs(x="Time (weeks)", y="CTL (cells/mL)",
           title="HCV-specific CTL Dynamics (with exhaustion)") +
      theme_bw(base_size=12)
  })

  output$nk_treg_plot <- renderPlot({
    d <- sim_data() %>%
      select(time_wk, NK_out, Treg_out) %>%
      pivot_longer(-time_wk, names_to="cell", values_to="count")
    ggplot(d, aes(time_wk, count, colour=cell)) +
      geom_line(linewidth=1) +
      scale_colour_manual(values=c(NK_out="#F39C12", Treg_out="#27AE60")) +
      labs(x="Time (weeks)", y="Cells/mL", colour=NULL,
           title="NK Cells & Treg") +
      theme_bw(base_size=12)
  })

  output$immune_ratio_plot <- renderPlot({
    d <- sim_data() %>%
      mutate(ctl_treg = CTL_out / (Treg_out + 1))
    ggplot(d, aes(time_wk, ctl_treg)) +
      geom_line(colour="#8E44AD", linewidth=1) +
      geom_hline(yintercept=5, linetype="dashed", colour="grey50") +
      labs(x="Time (weeks)", y="CTL / Treg Ratio",
           title="CTL/Treg Immune Balance") +
      theme_bw(base_size=12)
  })
}

shinyApp(ui, server)
