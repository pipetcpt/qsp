## ============================================================
## NMOSD QSP Shiny Dashboard
## Neuromyelitis Optica Spectrum Disorder
## Interactive simulation of AQP4-IgG / Complement / Biologics
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
## Inline mrgsolve model (compact version for Shiny)
## ============================================================

nmo_code <- '
$PARAM
kB_prod=0.02; kB_death=0.05; kB_act=0.01; kGC=0.10
kPB_diff=0.08; kPB_death=0.15; kPC_death=0.003
kAb_prod=0.5; kAb_CL=0.02; BAFF_stim=1.2; IL6_stim=1.5
C5_prod=10.0; C5_CL=0.14; C5_kon=0.005
kMAC_form=0.3; kMAC_CL=0.5
kAst_death=0.1; kAst_rep=0.005; Ast0=100; MAC_EC50=5.0; MAC_H=2.0
kLes_form=0.05; kLes_res=0.01
kEDSS_inc=0.002; kEDSS_rec=0.003; EDSS_max=10.0
kIL6_prod=0.5; kIL6_CL=0.3; kTNF_prod=0.2; kTNF_CL=0.5
kOligo_death=0.05; kOligo_rep=0.02
kNfL_rel=0.1; kNfL_CL=0.05
Ke_CL=0.022; Ke_V1=4.3; Ke_V2=3.2; Ke_Q=0.16
Ke_C5_kon=0.5; Ke_C5_koff=0.001; Ke_C5_kdeg=0.1
Ki_CL=0.019; Ki_V1=3.0; Ki_V2=2.5; Ki_Q=0.12; Ki_Bkill=0.5
Ks_ka=0.18; Ks_F=0.79; Ks_CL=0.65; Ks_V1=3.6; Ks_V2=2.0; Ks_Q=0.5
Ks_IL6R_kon=2.0; Ks_IL6R_koff=0.002
Kr_CL=0.33; Kr_V1=3.1; Kr_V2=2.8; Kr_Q=0.45; Kr_Bkill=0.4
Kp_ka=1.5; Kp_F=0.82; Kp_CL=14.4; Kp_V=47
Kp_Emax=0.8; Kp_EC50=50
Km_ka=0.9; Km_CL=12; Km_V=3.3; Km_Emax=0.7; Km_EC50=1.5

$CMT
Bnaive Bact PB PC Ab C5 EC_C5cx MAC Ast Oligo Lesion EDSS NfL IL6 TNFa
Eculi_C1 Eculi_C2 Ineb_C1 Ineb_C2 Satra_dep Satra_C1 Satra_C2 Satra_cx
Ritu_C1 Ritu_C2 Pred_gut Pred_C1 MPA_gut MPA_C1

$INIT
Bnaive=5; Bact=0.5; PB=0.2; PC=1.0; Ab=50.0; C5=70.0;
EC_C5cx=0; MAC=0.1; Ast=100; Oligo=100; Lesion=0.2;
EDSS=2.0; NfL=25; IL6=5; TNFa=3;
Eculi_C1=0; Eculi_C2=0; Ineb_C1=0; Ineb_C2=0; Satra_dep=0;
Satra_C1=0; Satra_C2=0; Satra_cx=0; Ritu_C1=0; Ritu_C2=0;
Pred_gut=0; Pred_C1=0; MPA_gut=0; MPA_C1=0;

$ODE
double Ke_k10 = (Ke_CL*24)/Ke_V1; double Ke_k12 = (Ke_Q*24)/Ke_V1; double Ke_k21 = (Ke_Q*24)/Ke_V2;
double Ki_k10 = (Ki_CL*24)/Ki_V1; double Ki_k12 = (Ki_Q*24)/Ki_V1; double Ki_k21 = (Ki_Q*24)/Ki_V2;
double Kr_k10 = Kr_CL/Kr_V1; double Kr_k12 = Kr_Q/Kr_V1; double Kr_k21 = Kr_Q/Kr_V2;
double Ks_k10 = Ks_CL/Ks_V1; double Ks_k12 = Ks_Q/Ks_V1; double Ks_k21 = Ks_Q/Ks_V2;
double Eculi_occ = Eculi_C1/(Eculi_C1+10.0);
double Satra_eff = Satra_cx/(Satra_cx+2.0);
double Pred_GR   = Kp_Emax*Pred_C1/(Kp_EC50+Pred_C1);
double MPA_inh   = Km_Emax*MPA_C1/(Km_EC50+MPA_C1);
double MAC_kill  = pow(MAC,MAC_H)/(pow(MAC_EC50,MAC_H)+pow(MAC,MAC_H));
double Bdrug_inh = 1.0-(1.0-MPA_inh)*(1.0-Pred_GR*0.5);
double IL6_boost = IL6_stim*IL6/(10.0+IL6);
double Ritu_kill = Kr_Bkill*Ritu_C1;
double Ineb_kill = Ki_Bkill*Ineb_C1;
double C5_free_eff = C5*(1.0-Eculi_occ);
double MAC_drive   = C5_kon*Ab*C5_free_eff/(Ab+20.0);
double Glu_tox = (100.0-Ast)/100.0;
double axon_dmg = (1.0-Oligo/100.0)+(100.0-Ast)/200.0;
double IL6_ast_prod = kIL6_prod*(100.0-Ast)/50.0;

dxdt_Bnaive  = kB_prod - kB_death*Bnaive - kB_act*Bnaive;
dxdt_Bact    = kB_act*Bnaive*BAFF_stim*(1-Bdrug_inh) - kGC*Bact - (Ineb_kill+Ritu_kill)*Bact;
dxdt_PB      = kGC*Bact*(1.0+IL6_boost) - kPB_diff*PB - kPB_death*PB - (Ineb_kill+Ritu_kill)*PB;
dxdt_PC      = kPB_diff*PB - kPC_death*PC - Ineb_kill*0.3*PC;
dxdt_Ab      = kAb_prod*PC - kAb_CL*Ab;
dxdt_C5      = C5_prod - C5_CL*C5 + Ke_C5_koff*EC_C5cx - Ke_C5_kon*Eculi_C1*C5*(1-Eculi_occ);
dxdt_EC_C5cx = Ke_C5_kon*Eculi_C1*C5 - Ke_C5_koff*EC_C5cx - Ke_C5_kdeg*EC_C5cx;
dxdt_MAC     = kMAC_form*MAC_drive*(1.0-Eculi_occ) - kMAC_CL*MAC;
dxdt_Ast     = kAst_rep*(100.0-Ast)*(Ast/100.0) - kAst_death*MAC_kill*Ast + Pred_GR*0.3*(100.0-Ast);
dxdt_Oligo   = kOligo_rep*(100.0-Oligo) - kOligo_death*Glu_tox*Oligo;
dxdt_Lesion  = kLes_form*(100.0-Ast)/100.0*MAC_kill - kLes_res*Lesion;
dxdt_EDSS    = kEDSS_inc*Lesion - kEDSS_rec*(EDSS>2.0?(EDSS-2.0):0.0);
dxdt_NfL     = kNfL_rel*axon_dmg*100.0 - kNfL_CL*NfL;
dxdt_IL6     = kIL6_prod + IL6_ast_prod - kIL6_CL*IL6*(1.0+Satra_eff);
dxdt_TNFa    = kTNF_prod*(100.0-Ast)/50.0 - kTNF_CL*TNFa - Pred_GR*0.4*TNFa;
dxdt_Eculi_C1 = -Ke_k10*Eculi_C1 - Ke_k12*Eculi_C1 + Ke_k21*Eculi_C2 - Ke_C5_kon*Eculi_C1*C5 + Ke_C5_koff*EC_C5cx;
dxdt_Eculi_C2 = Ke_k12*Eculi_C1 - Ke_k21*Eculi_C2;
dxdt_Ineb_C1  = -(Ki_k10+Ki_k12)*Ineb_C1 + Ki_k21*Ineb_C2;
dxdt_Ineb_C2  = Ki_k12*Ineb_C1 - Ki_k21*Ineb_C2;
dxdt_Satra_dep = -Ks_ka*Satra_dep;
dxdt_Satra_C1 = Ks_ka*Ks_F*Satra_dep/Ks_V1 - Ks_k10*Satra_C1 - Ks_k12*Satra_C1 + Ks_k21*Satra_C2 - Ks_IL6R_kon*Satra_C1*(10.0-Satra_cx) + Ks_IL6R_koff*Satra_cx;
dxdt_Satra_C2 = Ks_k12*Satra_C1 - Ks_k21*Satra_C2;
dxdt_Satra_cx = Ks_IL6R_kon*Satra_C1*(10.0-Satra_cx) - Ks_IL6R_koff*Satra_cx - 0.05*Satra_cx;
dxdt_Ritu_C1  = -(Kr_k10+Kr_k12)*Ritu_C1 + Kr_k21*Ritu_C2;
dxdt_Ritu_C2  = Kr_k12*Ritu_C1 - Kr_k21*Ritu_C2;
dxdt_Pred_gut = -(Kp_ka*24)*Pred_gut;
dxdt_Pred_C1  = (Kp_ka*24)*Kp_F*Pred_gut/Kp_V - (Kp_CL/Kp_V)*Pred_C1;
dxdt_MPA_gut  = -(Km_ka*24)*MPA_gut;
dxdt_MPA_C1   = (Km_ka*24)*MPA_gut/Km_V - (Km_CL/Km_V)*MPA_C1;

$TABLE
double GFAP_ser  = (100.0-Ast)*0.8+5.0;
double Bcell_pct = (PB+Bact)/5.8*100.0;

$CAPTURE GFAP_ser Bcell_pct
'

## Compile once at startup
mod <- mcode("NMOSD_Shiny", nmo_code)

## ============================================================
## Dosing helpers
## ============================================================
make_ev_dose <- function(drug, dose_nmol, times, rate_val = -2) {
  cmt_map <- c(
    "Eculizumab" = "Eculi_C1",
    "Inebilizumab" = "Ineb_C1",
    "Satralizumab" = "Satra_dep",
    "Rituximab"  = "Ritu_C1",
    "Prednisolone" = "Pred_gut",
    "MMF"        = "MPA_gut"
  )
  cmt_name <- cmt_map[drug]
  if (is.na(cmt_name)) return(NULL)
  is_sc <- drug %in% c("Satralizumab", "Prednisolone", "MMF")
  if (is_sc) {
    ev(cmt = cmt_name, amt = dose_nmol, time = times)
  } else {
    ev(cmt = cmt_name, amt = dose_nmol, time = times, rate = rate_val)
  }
}

build_dosing_regimen <- function(drugs_selected, n_months) {
  ev_list <- list()
  end_day <- n_months * 30

  if ("Eculizumab" %in% drugs_selected) {
    times_induct <- seq(0, 21, by = 7)
    times_maint  <- seq(28, end_day, by = 14)
    ev_list[["Eculi"]] <- make_ev_dose("Eculizumab", 1414,
                                        c(times_induct, times_maint))
  }
  if ("Inebilizumab" %in% drugs_selected) {
    times <- c(0, 14, seq(180, end_day, by = 180))
    ev_list[["Ineb"]] <- make_ev_dose("Inebilizumab", 667, times)
  }
  if ("Satralizumab" %in% drugs_selected) {
    times_induct <- c(0, 28, 56)
    times_maint  <- seq(112, end_day, by = 56)
    ev_list[["Satra"]] <- make_ev_dose("Satralizumab", 811,
                                        c(times_induct, times_maint))
  }
  if ("Rituximab" %in% drugs_selected) {
    times <- c(0, 14, seq(180, end_day, by = 180))
    ev_list[["Ritu"]] <- make_ev_dose("Rituximab", 2242, times)
  }
  if ("Prednisolone" %in% drugs_selected) {
    times <- seq(0, 4, by = 1)
    ev_list[["Pred"]] <- make_ev_dose("Prednisolone", 60e6, times)
  }
  if ("MMF" %in% drugs_selected) {
    times <- seq(0, end_day, by = 0.5)
    ev_list[["MMF"]] <- make_ev_dose("MMF", 1000, times)
  }

  if (length(ev_list) == 0) {
    return(ev(time = 0, amt = 0, cmt = 1))
  }
  Reduce(c, ev_list)
}

## ============================================================
## UI
## ============================================================

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "NMOSD QSP Simulator"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "patient",    icon = icon("user")),
      menuItem("Drug PK",            tabName = "pk",         icon = icon("pills")),
      menuItem("Disease Biomarkers", tabName = "biomarkers", icon = icon("vials")),
      menuItem("Clinical Endpoints", tabName = "endpoints",  icon = icon("chart-line")),
      menuItem("Scenario Comparison",tabName = "scenarios",  icon = icon("balance-scale")),
      menuItem("Mechanistic Map",    tabName = "map",        icon = icon("project-diagram"))
    ),
    hr(),
    h5("Treatment Selection", style = "color:white; padding-left:15px;"),
    checkboxGroupInput("drugs",
      label = NULL,
      choices = c("Eculizumab", "Inebilizumab", "Satralizumab",
                  "Rituximab", "Prednisolone", "MMF"),
      selected = "Eculizumab"),
    sliderInput("n_months", "Follow-up (months):",
                min = 3, max = 36, value = 24, step = 3),
    hr(),
    h5("Patient Parameters", style = "color:white; padding-left:15px;"),
    sliderInput("baseline_edss", "Baseline EDSS:", 0, 10, 2.0, 0.5),
    sliderInput("aqp4_titer",    "AQP4-IgG titer (nmol/L):", 10, 200, 50, 5),
    sliderInput("disease_act",   "Disease Activity (0-3):", 0, 3, 1, 0.5),
    actionButton("run_sim", "Run Simulation",
                 class = "btn-success", style = "margin:10px; width:90%;")
  ),

  dashboardBody(
    tags$head(tags$style(HTML(
      ".content-wrapper { background-color: #f4f6f9; }
       .box { border-top: 3px solid #3c8dbc; }
       .nav-tabs-custom .nav-tabs li.active { border-top-color: #3c8dbc; }"
    ))),

    tabItems(

      ## ---- Tab 1: Patient Profile ----
      tabItem(tabName = "patient",
        fluidRow(
          valueBoxOutput("vb_edss",   width = 3),
          valueBoxOutput("vb_aqp4",   width = 3),
          valueBoxOutput("vb_nfl",    width = 3),
          valueBoxOutput("vb_arr",    width = 3)
        ),
        fluidRow(
          box(title = "Disease Activity Profile", width = 6, status = "primary",
            plotlyOutput("plt_disease_profile", height = "300px")),
          box(title = "Clinical Context", width = 6, status = "info",
            h4("NMOSD Diagnostic Criteria (2015 IPND)"),
            tags$ul(
              tags$li("AQP4-IgG seropositive NMOSD: 1 core clinical characteristic"),
              tags$li("AQP4-IgG seronegative: 2+ core characteristics with MRI criteria"),
              tags$li("Core characteristics: LETM, Optic Neuritis, Area Postrema Synd,"),
              tags$li("  Acute brainstem, Symptomatic narcolepsy/diencephalic, Cerebral")
            ),
            hr(),
            h4("Treatment Goals"),
            tags$ul(
              tags$li(strong("Acute attack:"), "High-dose IV methylprednisolone ± PLEX"),
              tags$li(strong("Relapse prevention:"), "Long-term immunotherapy"),
              tags$li(strong("Target:"), "ARR < 0.2, EDSS stabilization, AQP4-IgG ↓")
            )
          )
        ),
        fluidRow(
          box(title = "Simulation Parameters Summary", width = 12, status = "warning",
            tableOutput("param_summary"))
        )
      ),

      ## ---- Tab 2: Drug PK ----
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Drug Concentration-Time Profiles", width = 12, status = "primary",
            plotlyOutput("plt_pk", height = "400px"))
        ),
        fluidRow(
          box(title = "PK Parameters Reference", width = 6, status = "info",
            DTOutput("tbl_pk_params")),
          box(title = "Drug Mechanism of Action", width = 6, status = "success",
            tableOutput("tbl_moa"))
        )
      ),

      ## ---- Tab 3: Disease Biomarkers ----
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(title = "AQP4-IgG Titer", width = 6, status = "danger",
            plotlyOutput("plt_aqp4", height = "280px")),
          box(title = "Complement MAC & C5", width = 6, status = "warning",
            plotlyOutput("plt_complement", height = "280px"))
        ),
        fluidRow(
          box(title = "B-Cell Count (% Baseline)", width = 6, status = "primary",
            plotlyOutput("plt_bcell", height = "280px")),
          box(title = "Serum IL-6 & TNF-α", width = 6, status = "info",
            plotlyOutput("plt_cytokines", height = "280px"))
        ),
        fluidRow(
          box(title = "Serum NfL (Neurofilament Light)", width = 6, status = "danger",
            plotlyOutput("plt_nfl", height = "280px")),
          box(title = "serum GFAP (Astrocyte Damage)", width = 6, status = "warning",
            plotlyOutput("plt_gfap", height = "280px"))
        )
      ),

      ## ---- Tab 4: Clinical Endpoints ----
      tabItem(tabName = "endpoints",
        fluidRow(
          box(title = "EDSS Trajectory", width = 6, status = "danger",
            plotlyOutput("plt_edss", height = "300px")),
          box(title = "Astrocyte & Oligodendrocyte Viability", width = 6, status = "info",
            plotlyOutput("plt_tissue", height = "300px"))
        ),
        fluidRow(
          box(title = "Active Lesion Burden", width = 6, status = "warning",
            plotlyOutput("plt_lesion", height = "280px")),
          box(title = "Annualized Relapse Rate Estimate", width = 6, status = "primary",
            plotlyOutput("plt_arr", height = "280px"))
        ),
        fluidRow(
          box(title = "Endpoint Summary at Final Timepoint", width = 12, status = "success",
            DTOutput("tbl_endpoints"))
        )
      ),

      ## ---- Tab 5: Scenario Comparison ----
      tabItem(tabName = "scenarios",
        fluidRow(
          box(title = "Scenario Setup", width = 12, status = "info",
            p("Comparison of 5 fixed treatment regimens vs. current selected regimen."),
            p("All scenarios use same patient baseline. Run simulation first."))
        ),
        fluidRow(
          box(title = "EDSS — All Scenarios", width = 12, status = "primary",
            plotlyOutput("plt_scen_edss", height = "380px"))
        ),
        fluidRow(
          box(title = "AQP4-IgG — All Scenarios", width = 6, status = "warning",
            plotlyOutput("plt_scen_aqp4", height = "300px")),
          box(title = "NfL — All Scenarios", width = 6, status = "danger",
            plotlyOutput("plt_scen_nfl", height = "300px"))
        ),
        fluidRow(
          box(title = "Scenario Comparison Table (Final Timepoint)", width = 12, status = "success",
            DTOutput("tbl_scenarios"))
        )
      ),

      ## ---- Tab 6: Mechanistic Map ----
      tabItem(tabName = "map",
        fluidRow(
          box(title = "NMOSD QSP Mechanistic Map", width = 12, status = "primary",
            p("Full pathway map: Genetic triggers → B-cell/AQP4-IgG → Complement → Astrocyte damage → Lesions → Disability"),
            img(src = "../nmo_qsp_model.svg", style = "max-width:100%; height:auto;"),
            hr(),
            h4("Key Pathogenic Steps"),
            tags$ol(
              tags$li(strong("AQP4-IgG production:"), "Plasma cells secrete AQP4-specific IgG1 that crosses the BBB via FcRn transcytosis"),
              tags$li(strong("Complement activation:"), "AQP4-IgG1 binds to AQP4-M23 OAP aggregates on astrocyte endfeet, activating classical → C5 → MAC"),
              tags$li(strong("Astrocyte necrosis:"), "MAC creates membrane pores → Ca²⁺ influx → necrosis; secondary GFAP/AQP4 loss"),
              tags$li(strong("Neuroinflammation:"), "C5a + CXCL1 recruit neutrophils & eosinophils; ROS + Glu release damage oligodendrocytes"),
              tags$li(strong("LETM formation:"), "Extensive spinal cord astrocyte loss → LETM (≥3 vertebral segments); area postrema → hiccups"),
              tags$li(strong("Disability:"), "Each relapse event may leave permanent EDSS increment; F:M=9:1, median ARR=1.8 untreated")
            )
          )
        )
      )

    ) # end tabItems
  ) # end dashboardBody
) # end dashboardPage

## ============================================================
## Server
## ============================================================

server <- function(input, output, session) {

  ## ---- Reactive simulation ----
  sim_result <- eventReactive(input$run_sim, {
    req(input$n_months)
    withProgress(message = "Simulating NMOSD model...", {
      evs <- build_dosing_regimen(input$drugs, input$n_months)
      out <- mrgsim(mod,
                    events = evs,
                    end    = input$n_months * 30,
                    delta  = 1,
                    init   = list(EDSS = input$baseline_edss,
                                  Ab   = input$aqp4_titer,
                                  MAC  = 0.1 * input$disease_act + 0.1,
                                  Lesion = 0.1 * input$disease_act + 0.1)) %>%
        as_tibble()
      out
    })
  }, ignoreNULL = FALSE)

  ## ---- Multi-scenario comparison ----
  scenarios_result <- eventReactive(input$run_sim, {
    withProgress(message = "Running scenario comparison...", {
      scenario_list <- list(
        "No Treatment"             = list(),
        "Eculizumab"               = list("Eculizumab"),
        "Inebilizumab"             = list("Inebilizumab"),
        "Satralizumab"             = list("Satralizumab"),
        "Rituximab + MMF"          = list("Rituximab", "MMF"),
        "Current Selection"        = as.list(input$drugs)
      )
      bind_rows(lapply(names(scenario_list), function(nm) {
        drugs <- scenario_list[[nm]]
        evs <- build_dosing_regimen(drugs, input$n_months)
        out <- mrgsim(mod,
                      events = evs,
                      end    = input$n_months * 30,
                      delta  = 1,
                      init   = list(EDSS = input$baseline_edss,
                                    Ab   = input$aqp4_titer)) %>%
          as_tibble()
        out$Scenario <- nm
        out
      }))
    })
  }, ignoreNULL = FALSE)

  ## ---- Value boxes ----
  output$vb_edss <- renderValueBox({
    df <- sim_result()
    last_edss <- round(tail(df$EDSS, 1), 2)
    valueBox(last_edss, "Final EDSS", icon = icon("wheelchair"),
             color = ifelse(last_edss < 3, "green", ifelse(last_edss < 6, "yellow", "red")))
  })

  output$vb_aqp4 <- renderValueBox({
    df <- sim_result()
    last_ab <- round(tail(df$Ab, 1), 1)
    valueBox(paste0(last_ab, " nmol/L"), "AQP4-IgG (final)",
             icon = icon("vial"), color = "blue")
  })

  output$vb_nfl <- renderValueBox({
    df <- sim_result()
    last_nfl <- round(tail(df$NfL, 1), 1)
    valueBox(paste0(last_nfl, " pg/mL"), "serum NfL (final)",
             icon = icon("brain"), color = "orange")
  })

  output$vb_arr <- renderValueBox({
    df <- sim_result()
    arr_est <- round(mean(df$Lesion * df$Ab / 50.0 * 1.8 / input$n_months, na.rm = TRUE), 2)
    valueBox(arr_est, "Est. ARR", icon = icon("bolt"),
             color = ifelse(arr_est < 0.3, "green", ifelse(arr_est < 1.0, "yellow", "red")))
  })

  ## ---- Disease profile ----
  output$plt_disease_profile <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time)) +
      geom_line(aes(y = EDSS,   color = "EDSS"),      size = 1) +
      geom_line(aes(y = Ast/10, color = "Ast (%/10)"), size = 1, linetype = "dashed") +
      geom_line(aes(y = NfL/25, color = "NfL (norm)"), size = 1, linetype = "dotted") +
      scale_color_manual(values = c("EDSS"="#D32F2F","Ast (%/10)"="#1565C0","NfL (norm)"="#F57F17")) +
      labs(x = "Days", y = "Value", color = NULL) +
      theme_bw()
    ggplotly(p)
  })

  output$param_summary <- renderTable({
    data.frame(
      Parameter = c("Baseline EDSS","AQP4-IgG (nmol/L)","Disease Activity",
                    "Follow-up (months)","Selected Drugs"),
      Value     = c(input$baseline_edss, input$aqp4_titer, input$disease_act,
                    input$n_months, paste(input$drugs, collapse = " + "))
    )
  })

  ## ---- PK plots ----
  output$plt_pk <- renderPlotly({
    df <- sim_result()
    pk_df <- df %>%
      select(time, Eculi_C1, Ineb_C1, Satra_C1, Ritu_C1, Pred_C1, MPA_C1) %>%
      pivot_longer(-time, names_to = "Drug", values_to = "Conc") %>%
      filter(Conc > 0.01)
    p <- ggplot(pk_df, aes(time, Conc, color = Drug)) +
      geom_line(size = 0.8) +
      facet_wrap(~Drug, scales = "free_y", ncol = 3) +
      labs(x = "Days", y = "Concentration") +
      theme_bw() + theme(legend.position = "none")
    ggplotly(p)
  })

  output$tbl_pk_params <- renderDT({
    data.frame(
      Drug = c("Eculizumab","Inebilizumab","Satralizumab","Rituximab","Prednisolone","MMF"),
      Route = c("IV","IV","SC","IV","PO","PO"),
      t_half = c("~11d","~16d","~30d","~22d","~3h","~17h"),
      Vd_L   = c("7.5","5.5","3.6","3.1","47","3.3"),
      CL_Lday= c("0.53","0.46","0.65","0.33","14.4","12"),
      Target = c("C5","CD19+","IL-6R","CD20+","GR (NF-κB)","IMPDH")
    )
  }, options = list(dom = "t", pageLength = 10))

  output$tbl_moa <- renderTable({
    data.frame(
      Drug = c("Eculizumab","Inebilizumab","Satralizumab",
               "Rituximab","Ublituximab","Prednisolone","MMF/MPA","AZA","IVIG","PE"),
      `Mechanism of Action` = c(
        "Terminal complement (C5) inhibition → ↓ MAC formation",
        "Anti-CD19 mAb → broad B-cell + plasmablast depletion",
        "Anti-IL-6R (recycling) → ↓ Th17 + plasmablast survival",
        "Anti-CD20 mAb → ADCC + CDC → memory B depletion",
        "Anti-CD20 (glycoengineered) → enhanced ADCC",
        "GR activation → ↓ NF-κB → ↓ TNF-α, IL-6, BBB stabilize",
        "IMPDH inhibition → ↓ B/T-cell proliferation",
        "6-TGN (purine analog) → ↓ B/T-cell DNA synthesis",
        "FcRn saturation → ↑ IgG clearance → ↓ AQP4-IgG",
        "Physical antibody removal → rapid ↓ AQP4-IgG titer"
      )
    )
  })

  ## ---- Biomarker plots ----
  plt_line <- function(df, var, title, ylab, color = "#D32F2F") {
    p <- ggplot(df, aes(time, .data[[var]])) +
      geom_line(color = color, size = 1) +
      labs(title = title, x = "Days", y = ylab) +
      theme_bw()
    ggplotly(p)
  }

  output$plt_aqp4      <- renderPlotly({ plt_line(sim_result(), "Ab",    "AQP4-IgG Titer",         "nmol/L", "#D32F2F") })
  output$plt_nfl       <- renderPlotly({ plt_line(sim_result(), "NfL",   "Serum NfL",              "pg/mL",  "#FF8F00") })
  output$plt_gfap      <- renderPlotly({ plt_line(sim_result(), "GFAP_ser", "serum GFAP",          "pg/mL",  "#6A1B9A") })

  output$plt_complement <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time)) +
      geom_line(aes(y = MAC, color = "MAC"), size = 1) +
      geom_line(aes(y = C5/10, color = "C5/10"), size = 1, linetype = "dashed") +
      scale_color_manual(values = c("MAC"="#880E4F","C5/10"="#F48FB1")) +
      labs(x = "Days", y = "nmol/L", color = NULL) + theme_bw()
    ggplotly(p)
  })

  output$plt_bcell <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time, Bcell_pct)) +
      geom_line(color = "#1B5E20", size = 1) +
      geom_hline(yintercept = 100, linetype = "dashed", color = "grey50") +
      labs(x = "Days", y = "B-cell count (% baseline)") + theme_bw()
    ggplotly(p)
  })

  output$plt_cytokines <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time)) +
      geom_line(aes(y = IL6,  color = "IL-6"),  size = 1) +
      geom_line(aes(y = TNFa, color = "TNF-α"), size = 1, linetype = "dashed") +
      scale_color_manual(values = c("IL-6"="#F57F17","TNF-α"="#D32F2F")) +
      labs(x = "Days", y = "pg/mL", color = NULL) + theme_bw()
    ggplotly(p)
  })

  ## ---- Clinical endpoint plots ----
  output$plt_edss <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time, EDSS)) +
      geom_line(color = "#D32F2F", size = 1.2) +
      geom_hline(yintercept = input$baseline_edss, linetype = "dashed", color = "grey60") +
      annotate("text", x = 10, y = input$baseline_edss + 0.15, label = "Baseline",
               color = "grey50", size = 3) +
      ylim(0, 10) +
      labs(x = "Days", y = "EDSS (0-10)") + theme_bw()
    ggplotly(p)
  })

  output$plt_tissue <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time)) +
      geom_line(aes(y = Ast,   color = "Astrocyte (%)"), size = 1) +
      geom_line(aes(y = Oligo, color = "Oligodendrocyte (%)"), size = 1, linetype = "dashed") +
      scale_color_manual(values = c("Astrocyte (%)"="#1565C0","Oligodendrocyte (%)"="#2E7D32")) +
      ylim(0, 110) +
      labs(x = "Days", y = "Cell Viability (%)", color = NULL) + theme_bw()
    ggplotly(p)
  })

  output$plt_lesion <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time, Lesion)) +
      geom_line(color = "#880E4F", size = 1) +
      labs(x = "Days", y = "Active Lesion Burden (a.u.)") + theme_bw()
    ggplotly(p)
  })

  output$plt_arr <- renderPlotly({
    df <- sim_result()
    df$ARR_roll <- df$Lesion * df$Ab / 50.0 * 1.8 / input$n_months
    p <- ggplot(df, aes(time, ARR_roll)) +
      geom_line(color = "#E65100", size = 1) +
      labs(x = "Days", y = "Estimated ARR (rolling)") + theme_bw()
    ggplotly(p)
  })

  output$tbl_endpoints <- renderDT({
    df <- sim_result()
    last <- tail(df, 1)
    data.frame(
      Endpoint         = c("EDSS","ΔEDSS from baseline","AQP4-IgG (nmol/L)",
                           "MAC (nmol/L-eq)","Astrocyte (%)","Oligodendrocyte (%)",
                           "NfL (pg/mL)","GFAP (pg/mL)","IL-6 (pg/mL)",
                           "Active Lesion Burden"),
      `Baseline`       = c(input$baseline_edss, 0, input$aqp4_titer,
                           0.1*input$disease_act+0.1, 100, 100, 25, 5, 5,
                           0.1*input$disease_act+0.1),
      `Final Value`    = round(c(last$EDSS, last$EDSS - input$baseline_edss,
                                 last$Ab, last$MAC, last$Ast, last$Oligo,
                                 last$NfL, last$GFAP_ser, last$IL6, last$Lesion), 2)
    )
  }, options = list(dom = "t", pageLength = 15))

  ## ---- Scenario comparison ----
  scen_summary <- reactive({
    df <- scenarios_result()
    df %>%
      group_by(Scenario, time) %>%
      summarise(
        EDSS_med = median(EDSS, na.rm = TRUE),
        Ab_med   = median(Ab,   na.rm = TRUE),
        NfL_med  = median(NfL,  na.rm = TRUE),
        .groups  = "drop"
      )
  })

  scen_cols <- c(
    "No Treatment"    = "#D32F2F",
    "Eculizumab"      = "#1565C0",
    "Inebilizumab"    = "#2E7D32",
    "Satralizumab"    = "#F57F17",
    "Rituximab + MMF" = "#6A1B9A",
    "Current Selection" = "#00838F"
  )

  output$plt_scen_edss <- renderPlotly({
    p <- ggplot(scen_summary(), aes(time, EDSS_med, color = Scenario)) +
      geom_line(size = 1) +
      scale_color_manual(values = scen_cols) +
      labs(x = "Days", y = "EDSS", color = "Scenario") + theme_bw()
    ggplotly(p)
  })

  output$plt_scen_aqp4 <- renderPlotly({
    p <- ggplot(scen_summary(), aes(time, Ab_med, color = Scenario)) +
      geom_line(size = 1) +
      scale_color_manual(values = scen_cols) +
      labs(x = "Days", y = "AQP4-IgG (nmol/L)") + theme_bw()
    ggplotly(p)
  })

  output$plt_scen_nfl <- renderPlotly({
    p <- ggplot(scen_summary(), aes(time, NfL_med, color = Scenario)) +
      geom_line(size = 1) +
      scale_color_manual(values = scen_cols) +
      labs(x = "Days", y = "NfL (pg/mL)") + theme_bw()
    ggplotly(p)
  })

  output$tbl_scenarios <- renderDT({
    df <- scenarios_result()
    df %>%
      filter(time == max(time)) %>%
      group_by(Scenario) %>%
      summarise(
        EDSS_final    = round(median(EDSS, na.rm = TRUE), 2),
        EDSS_change   = round(median(EDSS, na.rm = TRUE) - input$baseline_edss, 2),
        AQP4IgG_final = round(median(Ab, na.rm = TRUE), 1),
        Ast_final     = round(median(Ast, na.rm = TRUE), 1),
        NfL_final     = round(median(NfL, na.rm = TRUE), 1),
        Lesion_final  = round(median(Lesion, na.rm = TRUE), 3),
        .groups = "drop"
      ) %>%
      arrange(EDSS_final)
  }, options = list(dom = "t", pageLength = 10))

}

## ============================================================
## Run App
## ============================================================
shinyApp(ui, server)
