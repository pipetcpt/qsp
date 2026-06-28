################################################################################
# Allergic Rhinitis QSP – Interactive Shiny Dashboard
# 6 Tabs: Patient Profile · PK · Mediators/Biomarkers · Symptoms/TNSS ·
#         Scenario Comparison · Mechanistic Map
################################################################################

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mrgsolve)
library(DT)

# ── Inline mrgsolve model code (minimal, Shiny-compatible) ─────────────────

AR_MODEL_CODE <- '
$PARAM
  ALLERGEN_SS=1.0 K_ALLERGEN=0.1
  KSY_IGE=0.002 KDEG_IGE=0.005 K_BIND_MAST=0.1 K_OFF_MAST=0.001 MAST_TOTAL=1.0
  EC50_CROSS=0.5 HILL_CROSS=2.0 KDEG_HIST=2.0 KDEG_LT=1.0
  KHIST_PROD=5.0 KLT_PROD=2.0 KMAST_REC=0.05
  KSY_TH2=0.05 KDEG_TH2=0.1 KSY_IL4=1.0 KSY_IL5=0.8 KSY_IL13=1.2
  KDEG_IL4=0.5 KDEG_IL5=0.3 KDEG_IL13=0.4 TH2_BASE=0.5
  EOS_BLOOD_0=300 KEO_PROD=0.05 KEO_SURV=0.02 KEO_DEATH=0.08
  KEO_MIGRATE=0.01 KCHEMOKINE=2.0 KEO_TISSUE0=10.0 KEO_TIS_DEATH=0.05
  HIST_EC50=1.0 LT_EC50=0.5 EOS_EC50=50.0
  SNEEZE_MAX=3.0 RHINO_MAX=3.0 CONG_MAX=3.0 PRUR_MAX=3.0
  KA_CETI=0.9 CL_CETI=7.0 VD_CETI=70.0 F_CETI=0.70 H1_IC50_CETI=15.0 H1_HILL=1.0
  KA_FP=1.5 CL_FP_LOCAL=2.0 VD_FP_LOCAL=1.0 GR_IC50_FP=0.5 GR_HILL=1.2 FP_DOSE_BIOCONV=200.0
  KA_MLKT=0.5 CL_MLKT=45.0 VD_MLKT=10.0 F_MLKT=0.64 CYSLTR1_IC50=2.0
  KA_OMA=0.004 CL_OMA=0.14 VC_OMA=3140.0 VP_OMA=2360.0 Q_OMA=0.40 KON_IGE=1.0 KOFF_IGE=0.0001
$INIT
  AG=0 IGE_FREE=50 IGE_MAST=0.5 MAST_ACT=0 MAST_CHG=1.0
  HISTAMINE=0 CYS_LT=0 TH2=0.5 IL4=5 IL5=3 IL13=8
  EOS_B=300 EOS_N=10
  CETI_D=0 CETI_C=0 FP_LOC=0 MLKT_D=0 MLKT_C=0
  OMA_D=0 OMA_C=0 OMA_P=0 OMA_IGE=0
$ODE
  double Ag=AG;
  dxdt_AG=-K_ALLERGEN*AG;
  double degen_free=KDEG_IGE*IGE_FREE;
  double bind_rate=K_BIND_MAST*IGE_FREE*(MAST_TOTAL-IGE_MAST);
  double unbind_rate=K_OFF_MAST*IGE_MAST;
  double oma_cp=OMA_C/VC_OMA;
  double ige_capture=KON_IGE*oma_cp*IGE_FREE-KOFF_IGE*OMA_IGE;
  dxdt_IGE_FREE=KSY_IGE*200.0-degen_free-bind_rate+unbind_rate-ige_capture;
  dxdt_IGE_MAST=bind_rate-unbind_rate;
  dxdt_OMA_IGE=ige_capture;
  double crosslink_frac=pow(Ag,HILL_CROSS)/(pow(EC50_CROSS,HILL_CROSS)+pow(Ag,HILL_CROSS));
  double mast_trigger=crosslink_frac*IGE_MAST*MAST_CHG;
  dxdt_MAST_ACT=mast_trigger-0.5*MAST_ACT;
  dxdt_MAST_CHG=KMAST_REC*(1.0-MAST_CHG)-mast_trigger*MAST_CHG;
  dxdt_HISTAMINE=KHIST_PROD*MAST_ACT-KDEG_HIST*HISTAMINE;
  dxdt_CYS_LT=KLT_PROD*MAST_ACT-KDEG_LT*CYS_LT;
  double th2_drive=TH2_BASE+KSY_TH2*Ag*IGE_FREE/(1.0+IGE_FREE);
  dxdt_TH2=th2_drive-KDEG_TH2*TH2;
  dxdt_IL4=KSY_IL4*TH2-KDEG_IL4*IL4;
  dxdt_IL5=KSY_IL5*TH2-KDEG_IL5*IL5;
  dxdt_IL13=KSY_IL13*TH2-KDEG_IL13*IL13;
  double eos_prod=KEO_PROD*IL5+KEO_SURV*EOS_B*IL5/(1.0+IL5);
  double eos_death_b=KEO_DEATH*EOS_B;
  double chemokine=KCHEMOKINE*IL5*IL13;
  double eos_migrate=KEO_MIGRATE*chemokine*EOS_B;
  double eos_death_n=KEO_TIS_DEATH*EOS_N;
  dxdt_EOS_B=eos_prod-eos_death_b-eos_migrate;
  dxdt_EOS_N=eos_migrate-eos_death_n;
  double ke_ceti=CL_CETI/VD_CETI;
  dxdt_CETI_D=-KA_CETI*CETI_D;
  dxdt_CETI_C=KA_CETI*F_CETI*CETI_D-ke_ceti*CETI_C;
  double ceti_cp=CETI_C/VD_CETI*1000.0;
  dxdt_FP_LOC=KA_FP*FP_DOSE_BIOCONV/VD_FP_LOCAL-CL_FP_LOCAL*FP_LOC;
  double ke_mlkt=CL_MLKT/VD_MLKT;
  dxdt_MLKT_D=-KA_MLKT*MLKT_D;
  dxdt_MLKT_C=KA_MLKT*F_MLKT*MLKT_D-ke_mlkt*MLKT_C;
  double mlkt_cp=MLKT_C/VD_MLKT*1000.0;
  dxdt_OMA_D=-KA_OMA*OMA_D;
  dxdt_OMA_C=KA_OMA*OMA_D-(CL_OMA/VC_OMA+Q_OMA/VC_OMA)*OMA_C+Q_OMA/VP_OMA*OMA_P;
  dxdt_OMA_P=Q_OMA/VC_OMA*OMA_C-Q_OMA/VP_OMA*OMA_P;
  double h1ro=ceti_cp/(H1_IC50_CETI+ceti_cp);
  double gr_occ=pow(FP_LOC,GR_HILL)/(pow(GR_IC50_FP,GR_HILL)+pow(FP_LOC,GR_HILL));
  dxdt_IL4-=gr_occ*0.6*IL4;
  dxdt_IL5-=gr_occ*0.7*IL5;
  dxdt_IL13-=gr_occ*0.6*IL13;
  dxdt_EOS_N-=gr_occ*0.5*EOS_N;
$TABLE
  double CETI_CP=CETI_C/VD_CETI*1000.0;
  double FP_NM=FP_LOC;
  double MLKT_CP=MLKT_C/VD_MLKT*1000.0;
  double OMA_CP=OMA_C/(VC_OMA/1000.0);
  double H1_RO=CETI_CP/(H1_IC50_CETI+CETI_CP)*100.0;
  double GR_OCC=pow(FP_NM,GR_HILL)/(pow(GR_IC50_FP,GR_HILL)+pow(FP_NM,GR_HILL))*100.0;
  double CYSLTR1_INH=MLKT_CP/(CYSLTR1_IC50+MLKT_CP)*100.0;
  double H1_EFF=HISTAMINE*(1.0-H1_RO/100.0);
  double LT_EFF=CYS_LT*(1.0-CYSLTR1_INH/100.0);
  double SNEEZE=SNEEZE_MAX*H1_EFF/(HIST_EC50+H1_EFF);
  double RHINORRHEA=RHINO_MAX*(0.6*H1_EFF+0.4*LT_EFF)/(HIST_EC50+0.6*H1_EFF+0.4*LT_EFF);
  double EOS_CONG=EOS_N/(EOS_EC50+EOS_N);
  double CONGESTION=CONG_MAX*(0.5*LT_EFF/(LT_EC50+LT_EFF)+0.3*EOS_CONG+0.2);
  double IL13_N=IL13/(IL13+8.0);
  double PRURITUS=PRUR_MAX*(0.7*H1_EFF/(HIST_EC50+H1_EFF)+0.3*IL13_N);
  double TNSS=SNEEZE+RHINORRHEA+CONGESTION+PRURITUS;
  capture CETI_CP MLKT_CP OMA_CP FP_NM H1_RO GR_OCC CYSLTR1_INH
  capture SNEEZE RHINORRHEA CONGESTION PRURITUS TNSS
  capture IGE_FREE EOS_B EOS_N HISTAMINE CYS_LT IL4 IL5 IL13 MAST_ACT
'

AR <- suppressMessages(mcode("AR_shiny", AR_MODEL_CODE))

run_sim <- function(allergen_amt, ceti_dose, fp_dose, mlkt_dose, oma_dose,
                    sim_days = 84, challenge_day = 28) {
  ev_list <- list(ev(cmt = "AG", amt = allergen_amt, time = challenge_day * 24))
  if (ceti_dose > 0)
    ev_list[["ceti"]] <- ev(cmt = "CETI_D", amt = ceti_dose, time = 0, ii = 24, addl = sim_days - 1)
  if (fp_dose > 0)
    ev_list[["fp"]]   <- ev(cmt = "FP_LOC", amt = fp_dose,   time = 0, ii = 24, addl = sim_days - 1)
  if (mlkt_dose > 0)
    ev_list[["mlkt"]] <- ev(cmt = "MLKT_D", amt = mlkt_dose, time = 0, ii = 24, addl = sim_days - 1)
  if (oma_dose > 0)
    ev_list[["oma"]]  <- ev(cmt = "OMA_D",  amt = oma_dose,  time = 0, ii = 28 * 24, addl = 2)
  evs <- Reduce(c, ev_list)

  AR %>%
    ev(evs) %>%
    mrgsim(end = sim_days * 24, delta = 1) %>%
    as.data.frame() %>%
    mutate(time_d = time / 24)
}

# ── UI ─────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Allergic Rhinitis QSP Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",       tabName = "profile",    icon = icon("user")),
      menuItem("Drug PK",               tabName = "pk",         icon = icon("pills")),
      menuItem("Mediators & Biomarkers",tabName = "biomarkers", icon = icon("flask")),
      menuItem("Symptoms & TNSS",       tabName = "symptoms",   icon = icon("chart-line")),
      menuItem("Scenario Comparison",   tabName = "comparison", icon = icon("table")),
      menuItem("Mechanistic Map",       tabName = "map",        icon = icon("project-diagram"))
    ),

    hr(),
    h4("Patient Settings", style = "color:white; margin-left:10px"),
    sliderInput("allergen_amt",  "Allergen Challenge (AU)", 1, 10, 5, 0.5),
    checkboxInput("use_ceti",    "Cetirizine 10 mg QD",     FALSE),
    checkboxInput("use_fp",      "Fluticasone 200 μg/d",    FALSE),
    checkboxInput("use_mlkt",    "Montelukast 10 mg QD",    FALSE),
    checkboxInput("use_oma",     "Omalizumab 300 mg q4w",   FALSE),
    sliderInput("sim_days",      "Simulation (days)", 28, 168, 84, 14),
    actionButton("run_btn",      "Run Simulation", icon = icon("play"),
                 class = "btn-success", style = "margin:10px; width:200px")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .box { border-top: 3px solid #3c8dbc; }
      .value-box { min-height: 80px; }
    "))),
    tabItems(

      # ── Tab 1: Patient Profile ────────────────────────────────────────────
      tabItem(tabName = "profile",
        fluidRow(
          box(title = "Patient Atopic Profile & Biomarkers at Baseline",
              width = 12, status = "primary", solidHeader = TRUE,
              fluidRow(
                valueBoxOutput("vb_ige",   width = 3),
                valueBoxOutput("vb_eos",   width = 3),
                valueBoxOutput("vb_il5",   width = 3),
                valueBoxOutput("vb_tnss0", width = 3)
              )
          )
        ),
        fluidRow(
          box(title = "ARIA Classification Criteria", width = 6, status = "info",
              solidHeader = TRUE,
              tableOutput("aria_table")
          ),
          box(title = "AR Pathophysiology Overview", width = 6, status = "warning",
              solidHeader = TRUE,
              p(strong("Type I Hypersensitivity – IgE-mediated Allergic Response")),
              p("1. Sensitization: Allergen → DC → Th2 polarization → B cell IgE class switch"),
              p("2. Early Phase (0-60 min): Allergen cross-links mast cell IgE → degranulation"),
              p("   Mediators: Histamine, PGD2, LTC4/D4, PAF"),
              p("   Symptoms: Sneezing, rhinorrhea, pruritus (H1-driven)"),
              p("3. Late Phase (4-24h): Eosinophil/Th2 infiltration → IL-4/5/13 surge"),
              p("   Symptoms: Congestion, nasal blockage (LT/Eos-driven)"),
              p("4. Chronic: Mucosal remodeling, goblet hyperplasia, anosmia"),
              br(),
              p(strong("Key drug targets:")),
              tags$ul(
                tags$li("H1R: Cetirizine, loratadine, fexofenadine (2nd-gen AH)"),
                tags$li("GR (nasal): Fluticasone, mometasone, budesonide (INCS)"),
                tags$li("CysLT1R: Montelukast (LTRA)"),
                tags$li("IgE: Omalizumab (anti-IgE biologic)"),
                tags$li("IL-4Rα: Dupilumab | IL-5: Mepolizumab | AIT: SCIT/SLIT")
              )
          )
        )
      ),

      # ── Tab 2: Drug PK ────────────────────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Plasma Concentration – Time Profiles", width = 12,
              status = "primary", solidHeader = TRUE,
              plotOutput("pk_plot", height = "400px"))
        ),
        fluidRow(
          box(title = "H1 Receptor Occupancy (Cetirizine)", width = 6,
              status = "info", solidHeader = TRUE,
              plotOutput("h1ro_plot", height = "300px")),
          box(title = "GR Occupancy (Fluticasone) & CysLT1 Inhibition (Montelukast)",
              width = 6, status = "success", solidHeader = TRUE,
              plotOutput("pd_occ_plot", height = "300px"))
        )
      ),

      # ── Tab 3: Mediators & Biomarkers ─────────────────────────────────────
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(title = "Inflammatory Mediators (Histamine & CysLTs)", width = 6,
              status = "warning", solidHeader = TRUE,
              plotOutput("mediator_plot", height = "300px")),
          box(title = "Th2 Cytokines (IL-4, IL-5, IL-13)", width = 6,
              status = "danger", solidHeader = TRUE,
              plotOutput("cytokine_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "Eosinophil Dynamics (Blood & Nasal Tissue)", width = 6,
              status = "primary", solidHeader = TRUE,
              plotOutput("eos_plot", height = "300px")),
          box(title = "Free IgE & Mast Cell Activation", width = 6,
              status = "info", solidHeader = TRUE,
              plotOutput("ige_mast_plot", height = "300px"))
        )
      ),

      # ── Tab 4: Symptoms & TNSS ────────────────────────────────────────────
      tabItem(tabName = "symptoms",
        fluidRow(
          box(title = "Total Nasal Symptom Score (TNSS 0-12)", width = 12,
              status = "primary", solidHeader = TRUE,
              plotOutput("tnss_plot", height = "350px"))
        ),
        fluidRow(
          box(title = "Individual Symptom Scores (0-3 each)", width = 8,
              status = "warning", solidHeader = TRUE,
              plotOutput("symptom_panel", height = "350px")),
          box(title = "Peak & Trough Summary", width = 4,
              status = "success", solidHeader = TRUE,
              tableOutput("sym_summary_tbl"))
        )
      ),

      # ── Tab 5: Scenario Comparison ────────────────────────────────────────
      tabItem(tabName = "comparison",
        fluidRow(
          box(title = "Multi-Scenario TNSS Comparison", width = 12,
              status = "primary", solidHeader = TRUE,
              plotOutput("scenario_tnss", height = "400px"))
        ),
        fluidRow(
          box(title = "Endpoint Summary Table (Week 12)", width = 12,
              status = "info", solidHeader = TRUE,
              DTOutput("scenario_table"))
        )
      ),

      # ── Tab 6: Mechanistic Map ────────────────────────────────────────────
      tabItem(tabName = "map",
        fluidRow(
          box(title = "Allergic Rhinitis QSP Mechanistic Map",
              width = 12, status = "primary", solidHeader = TRUE,
              p("The SVG mechanistic map is stored in the repository:"),
              tags$code("allergic-rhinitis/ar_qsp_model.svg"),
              br(), br(),
              p(strong("Model Components (12 subgraph clusters):")),
              tags$ol(
                tags$li("Allergen Exposure & Environmental Triggers"),
                tags$li("Epithelial Barrier & Innate Alarmins (TSLP, IL-33, ILC2)"),
                tags$li("APC Sensitization & Th2 Polarization"),
                tags$li("B Cell & IgE Class Switch Recombination"),
                tags$li("Mast Cell Sensitization & Allergen-triggered Activation"),
                tags$li("Preformed & Newly Synthesized Mediators (Histamine, PGD2, CysLTs)"),
                tags$li("Receptor-mediated Nasal Physiology (H1R, CysLT1R, DP1/CRTH2)"),
                tags$li("Late Phase Reaction & Eosinophil Inflammation"),
                tags$li("Nasal Pathophysiology & Symptom Endpoints (TNSS, RQLQ)"),
                tags$li("PK/PD: Oral H1-Antihistamines (Cetirizine, Loratadine, Fexofenadine)"),
                tags$li("PK/PD: Intranasal Corticosteroids (Fluticasone, Mometasone, Budesonide)"),
                tags$li("Biologics, LTRA & Allergen Immunotherapy (Omalizumab, Dupilumab, SCIT/SLIT)")
              )
          )
        )
      )
    )
  )
)

# ── Server ─────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  sim_data <- eventReactive(input$run_btn, {
    run_sim(
      allergen_amt = input$allergen_amt,
      ceti_dose    = if (input$use_ceti)  10  else 0,
      fp_dose      = if (input$use_fp)    450 else 0,
      mlkt_dose    = if (input$use_mlkt)  10  else 0,
      oma_dose     = if (input$use_oma)   300 else 0,
      sim_days     = input$sim_days
    )
  }, ignoreNULL = FALSE)

  sim_data_init <- reactive({
    run_sim(5, 0, 0, 0, 0, 84)
  })

  # Multi-scenario reactive
  scenario_data <- eventReactive(input$run_btn, {
    sc_list <- list(
      list(name="Natural History",         ceti=0,   fp=0,   mlkt=0,  oma=0,   col="#E53935"),
      list(name="Cetirizine",              ceti=10,  fp=0,   mlkt=0,  oma=0,   col="#1E88E5"),
      list(name="Fluticasone",             ceti=0,   fp=450, mlkt=0,  oma=0,   col="#43A047"),
      list(name="Montelukast",             ceti=0,   fp=0,   mlkt=10, oma=0,   col="#FB8C00"),
      list(name="Ceti + FP",              ceti=10,  fp=450, mlkt=0,  oma=0,   col="#8E24AA"),
      list(name="Omalizumab",             ceti=0,   fp=0,   mlkt=0,  oma=300, col="#00897B"),
      list(name="Triple (Ceti+FP+MLKT)",  ceti=10,  fp=450, mlkt=10, oma=0,   col="#6D4C41")
    )
    bind_rows(lapply(sc_list, function(sc) {
      run_sim(input$allergen_amt, sc$ceti, sc$fp, sc$mlkt, sc$oma, input$sim_days) %>%
        mutate(scenario = sc$name, color = sc$col)
    }))
  }, ignoreNULL = FALSE)

  # ── Value boxes ──
  output$vb_ige <- renderValueBox({
    d <- sim_data_init()
    valueBox(round(mean(d$IGE_FREE[d$time_d < 1]), 1), "Free IgE (IU/mL)",
             icon = icon("tint"), color = "purple")
  })
  output$vb_eos <- renderValueBox({
    d <- sim_data_init()
    valueBox(round(mean(d$EOS_B[d$time_d < 1])), "Blood Eos (cells/μL)",
             icon = icon("circle"), color = "red")
  })
  output$vb_il5 <- renderValueBox({
    d <- sim_data_init()
    valueBox(round(mean(d$IL5[d$time_d < 1]), 2), "IL-5 (pg/mL)",
             icon = icon("virus"), color = "yellow")
  })
  output$vb_tnss0 <- renderValueBox({
    d <- sim_data()
    valueBox(round(mean(d$TNSS[d$time_d < 28 & d$time_d > 0]), 2),
             "Baseline TNSS", icon = icon("thermometer-half"), color = "blue")
  })

  # ── ARIA table ──
  output$aria_table <- renderTable({
    data.frame(
      Duration      = c("Intermittent", "Persistent"),
      Definition    = c("< 4 days/week OR < 4 weeks",
                        "≥ 4 days/week AND ≥ 4 weeks"),
      Mild          = c("Sleep normal, activities normal, no troublesome symptoms",
                        ""),
      Moderate_Severe = c("One or more: sleep disturbance, impaired activities, troublesome symptoms",
                          "")
    )
  }, striped = TRUE, hover = TRUE, spacing = "s")

  # ── PK plots ──
  output$pk_plot <- renderPlot({
    df <- sim_data()
    df_long <- df %>%
      select(time_d, CETI_CP, MLKT_CP, OMA_CP) %>%
      pivot_longer(-time_d, names_to = "drug", values_to = "conc") %>%
      filter(conc > 0.001)
    ggplot(df_long, aes(x = time_d, y = conc, color = drug)) +
      geom_line(linewidth = 1) +
      facet_wrap(~drug, scales = "free_y") +
      labs(title = "Drug Plasma Concentrations",
           x = "Time (days)", y = "Conc. (ng/mL or μg/mL)") +
      theme_bw(base_size = 12) +
      theme(legend.position = "none")
  })

  output$h1ro_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = time_d, y = H1_RO)) +
      geom_line(color = "#1E88E5", linewidth = 1) +
      geom_hline(yintercept = 80, linetype = "dashed", color = "red") +
      annotate("text", x = 5, y = 82, label = "80% target (Yanai 1995)", color = "red", size = 3) +
      labs(title = "H1 Receptor Occupancy (Cetirizine)",
           x = "Time (days)", y = "H1-RO (%)") +
      ylim(0, 100) + theme_bw(base_size = 12)
  })

  output$pd_occ_plot <- renderPlot({
    df <- sim_data()
    df_long <- df %>%
      select(time_d, GR_OCC, CYSLTR1_INH) %>%
      pivot_longer(-time_d, names_to = "target", values_to = "inhibition")
    ggplot(df_long, aes(x = time_d, y = inhibition, color = target)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c("GR_OCC" = "#43A047", "CYSLTR1_INH" = "#FB8C00"),
                         labels = c("GR Occupancy (FP)", "CysLT1 Inhibition (MLKT)")) +
      labs(title = "GR & CysLT1 Receptor Occupancy",
           x = "Time (days)", y = "Occupancy / Inhibition (%)", color = "Target") +
      ylim(0, 100) + theme_bw(base_size = 12)
  })

  # ── Biomarker plots ──
  output$mediator_plot <- renderPlot({
    df <- sim_data()
    df_long <- df %>%
      select(time_d, HISTAMINE, CYS_LT) %>%
      pivot_longer(-time_d, names_to = "mediator", values_to = "level")
    ggplot(df_long, aes(x = time_d, y = level, color = mediator)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = 28, linetype = "dashed", color = "grey50") +
      scale_color_manual(values = c("HISTAMINE" = "#F44336", "CYS_LT" = "#FF9800"),
                         labels = c("Histamine (AU)", "CysLTs (AU)")) +
      labs(title = "Inflammatory Mediators",
           x = "Time (days)", y = "Mediator (normalized AU)", color = "") +
      theme_bw(base_size = 12)
  })

  output$cytokine_plot <- renderPlot({
    df <- sim_data()
    df_long <- df %>%
      select(time_d, IL4, IL5, IL13) %>%
      pivot_longer(-time_d, names_to = "cytokine", values_to = "pg_mL")
    ggplot(df_long, aes(x = time_d, y = pg_mL, color = cytokine)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = 28, linetype = "dashed", color = "grey50") +
      scale_color_manual(values = c("IL4" = "#E91E63", "IL5" = "#9C27B0", "IL13" = "#F44336")) +
      labs(title = "Th2 Cytokines",
           x = "Time (days)", y = "Cytokine (pg/mL)", color = "") +
      theme_bw(base_size = 12)
  })

  output$eos_plot <- renderPlot({
    df <- sim_data()
    df_long <- df %>%
      select(time_d, EOS_B, EOS_N) %>%
      pivot_longer(-time_d, names_to = "compartment", values_to = "cells_uL")
    ggplot(df_long, aes(x = time_d, y = cells_uL, color = compartment)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = 28, linetype = "dashed", color = "grey50") +
      scale_color_manual(values = c("EOS_B" = "#1565C0", "EOS_N" = "#E53935"),
                         labels = c("Blood Eos", "Nasal Tissue Eos")) +
      labs(title = "Eosinophil Dynamics",
           x = "Time (days)", y = "Eosinophils (cells/μL)", color = "") +
      theme_bw(base_size = 12)
  })

  output$ige_mast_plot <- renderPlot({
    df <- sim_data()
    par(mfrow = c(1, 1))
    ggplot(df, aes(x = time_d)) +
      geom_line(aes(y = IGE_FREE, color = "Free IgE (IU/mL)"), linewidth = 1) +
      geom_line(aes(y = MAST_ACT * 50, color = "Mast Cell Activation × 50"), linewidth = 1) +
      geom_vline(xintercept = 28, linetype = "dashed", color = "grey50") +
      scale_color_manual(values = c("Free IgE (IU/mL)" = "#9C27B0", "Mast Cell Activation × 50" = "#FF5722")) +
      labs(title = "Free IgE & Mast Cell Activation",
           x = "Time (days)", y = "IgE (IU/mL) | Activation (scaled)", color = "") +
      theme_bw(base_size = 12)
  })

  # ── Symptom plots ──
  output$tnss_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x = time_d, y = TNSS)) +
      geom_line(color = "#D32F2F", linewidth = 1.2) +
      geom_vline(xintercept = 28, linetype = "dashed", color = "grey40") +
      annotate("text", x = 28.5, y = 11, label = "Allergen Challenge", size = 3.5, hjust = 0) +
      ylim(0, 12) +
      labs(title = "Total Nasal Symptom Score (TNSS)",
           x = "Time (days)", y = "TNSS (0-12)") +
      theme_bw(base_size = 13)
  })

  output$symptom_panel <- renderPlot({
    df <- sim_data()
    df_long <- df %>%
      select(time_d, SNEEZE, RHINORRHEA, CONGESTION, PRURITUS) %>%
      pivot_longer(-time_d, names_to = "symptom", values_to = "score")
    ggplot(df_long, aes(x = time_d, y = score, color = symptom)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~symptom, ncol = 2) +
      geom_vline(xintercept = 28, linetype = "dashed", color = "grey50", alpha = 0.7) +
      scale_color_manual(values = c("SNEEZE" = "#1E88E5", "RHINORRHEA" = "#43A047",
                                    "CONGESTION" = "#E53935", "PRURITUS" = "#8E24AA")) +
      labs(title = "Individual Symptom Scores", x = "Time (days)", y = "Score (0-3)", color = "Symptom") +
      ylim(0, 3) + theme_bw(base_size = 11) + theme(legend.position = "none")
  })

  output$sym_summary_tbl <- renderTable({
    df <- sim_data()
    data.frame(
      Metric    = c("Pre-challenge TNSS", "Peak TNSS", "Wk 4-12 TNSS", "Max TNSS change"),
      Value     = c(
        round(mean(df$TNSS[df$time_d < 28], na.rm = TRUE), 2),
        round(max(df$TNSS[df$time_d >= 28 & df$time_d <= 35], na.rm = TRUE), 2),
        round(mean(df$TNSS[df$time_d >= 56], na.rm = TRUE), 2),
        round(max(df$TNSS) - mean(df$TNSS[df$time_d < 28], na.rm = TRUE), 2)
      )
    )
  }, striped = TRUE, bordered = TRUE)

  # ── Scenario comparison ──
  output$scenario_tnss <- renderPlot({
    df <- scenario_data()
    ggplot(df, aes(x = time_d, y = TNSS, color = scenario)) +
      geom_line(linewidth = 0.9) +
      scale_color_manual(values = setNames(unique(df$color), unique(df$scenario))) +
      geom_vline(xintercept = 28, linetype = "dashed", color = "grey40") +
      annotate("text", x = 28.5, y = 11.5, label = "Challenge", size = 3.5, hjust = 0) +
      ylim(0, 12) +
      labs(title = "TNSS by Treatment Scenario",
           x = "Time (days)", y = "TNSS (0-12)", color = "Scenario") +
      theme_bw(base_size = 12) +
      theme(legend.position = "right", legend.text = element_text(size = 9))
  })

  output$scenario_table <- renderDT({
    df <- scenario_data()
    tbl <- df %>%
      group_by(scenario) %>%
      summarise(
        "Pre-challenge TNSS" = round(mean(TNSS[time_d < 28], na.rm = TRUE), 2),
        "Peak TNSS"          = round(max(TNSS[time_d >= 28 & time_d <= 35], na.rm = TRUE), 2),
        "Wk 12 TNSS"         = round(mean(TNSS[time_d >= 77], na.rm = TRUE), 2),
        "Free IgE (IU/mL)"   = round(mean(IGE_FREE[time_d >= 77], na.rm = TRUE), 1),
        "Blood Eos (μL)"     = round(mean(EOS_B[time_d >= 77], na.rm = TRUE), 0),
        "Nasal Eos (μL)"     = round(mean(EOS_N[time_d >= 77], na.rm = TRUE), 1),
        .groups = "drop"
      )
    datatable(tbl, options = list(pageLength = 10, scrollX = TRUE),
              class = "table-striped table-hover table-sm")
  })
}

shinyApp(ui, server)
