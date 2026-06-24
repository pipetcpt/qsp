## =============================================================================
## Pemphigus Vulgaris QSP Interactive Dashboard
## Shiny App — 6-tab interactive simulation
##
## Tabs:
##   1. Patient Profile & Disease Setup
##   2. Drug PK — Plasma Concentration Profiles
##   3. Immune Biomarkers — B cells, Tfh, Treg
##   4. Disease Activity — PDAI, Anti-Dsg3, Dsg3 loss
##   5. Treatment Scenario Comparison
##   6. Clinical Endpoints & Safety
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(shinydashboard)
library(DT)
library(plotly)

## ─────────────────────────────────────────────────────────────────────────────
## MODEL CODE (embedded)
## ─────────────────────────────────────────────────────────────────────────────

pv_code <- '
$PARAM
Cl_pred=5.1, Vc_pred=35, Q_pred=2.5, Vp_pred=60, ka_pred=1.2, F_pred=0.82,
Cl_RTX=0.014, Vc_RTX=3.1, Q_RTX=0.18, Vp_RTX=4.4, kTMDD=0.002,
Cl_MPA=15, Vc_MPA=50, ka_MPA=0.8, F_MPA=0.94,
Cl_IVIg=0.005, Vc_IVIg=3.5,
kBN_prod=0.015, kBN_die=0.004, kBN_act=0.0003,
kBGC_die=0.010, kBM_form=0.008, kBM_die=0.0004,
kSLPC_form=0.012, kSLPC_die=0.040,
kLLPC_form=0.002, kLLPC_die=0.0003,
kTfh_prod=0.002, kTfh_die=0.008,
kTreg_prod=0.001, kTreg_die=0.005,
kAb3_prod=0.0006, kAb3_die=0.0014,
kAb1_prod=0.0003, kAb1_die=0.0014,
Dsg3_base=100, kDsg3_int=0.0004, kDsg3_rest=0.0008,
kPDAI_up=0.018, kPDAI_dn=0.004,
kComp_on=0.0003, kComp_off=0.010,
kBone_loss=0.0001, kBone_rest=0.00005,
EC50_pred_T=5, Emax_pred_T=0.9, EC50_pred_B=8, Emax_pred_B=0.85,
EC50_RTX_B=0.05, Emax_RTX_B=0.99, Hill_RTX=2,
EC50_MPA_GC=2, Emax_MPA_GC=0.8,
EC50_IVIg_Ab=1, Emax_IVIg_Ab=0.75,
EC50_EFG_Ab=0.1, Emax_EFG_Ab=0.9,
BN0=50, BGC0=5, BM0=20, SLPC0=3, LLPC0=2, Tfh0=10, Treg0=8

$CMT DEPOT_pred Cp1_pred Cp2_pred CR1 CR2 DEPOT_MPA CMPA DEPOT_IVIg CIVG CEFG
     BN BGC BM SLPC LLPC Tfh Treg Ab3 Ab1
     Dsg3_loss PDAI_state Comp_act Cort_bone

$MAIN
double Epred_T = Emax_pred_T * Cp1_pred / (EC50_pred_T + Cp1_pred + 1e-10);
double Epred_B = Emax_pred_B * Cp1_pred / (EC50_pred_B + Cp1_pred + 1e-10);
double CR1h = pow(CR1, Hill_RTX);
double EC50h = pow(EC50_RTX_B, Hill_RTX);
double ERTX = Emax_RTX_B * CR1h / (EC50h + CR1h + 1e-20);
double EMPA = Emax_MPA_GC * CMPA / (EC50_MPA_GC + CMPA + 1e-10);
double EIVIg = Emax_IVIg_Ab * CIVG / (EC50_IVIg_Ab + CIVG + 1e-10);
double EEFG  = Emax_EFG_Ab  * CEFG / (EC50_EFG_Ab  + CEFG + 1e-10);
double E_Ab_cat = fmax(EIVIg, EEFG);
double B_total  = BN + BGC + BM;
double Dsg3_avail = fmax(0.0, Dsg3_base - Dsg3_loss);
double Treg_sup   = Treg / (Treg + Treg0 + 1e-10);
double Tfh_stim   = Tfh  / (Tfh0  + 1e-10);

$ODE
dxdt_DEPOT_pred = -ka_pred * DEPOT_pred;
dxdt_Cp1_pred = (ka_pred*DEPOT_pred*F_pred)/Vc_pred - (Cl_pred/Vc_pred)*Cp1_pred
                - (Q_pred/Vc_pred)*Cp1_pred + (Q_pred/Vp_pred)*Cp2_pred;
dxdt_Cp2_pred = (Q_pred/Vc_pred)*Cp1_pred - (Q_pred/Vp_pred)*Cp2_pred;
dxdt_CR1 = -(Cl_RTX/Vc_RTX)*CR1 - (Q_RTX/Vc_RTX)*CR1 + (Q_RTX/Vp_RTX)*CR2
           - kTMDD*CR1*B_total;
dxdt_CR2 = (Q_RTX/Vc_RTX)*CR1 - (Q_RTX/Vp_RTX)*CR2;
dxdt_DEPOT_MPA = -ka_MPA * DEPOT_MPA;
dxdt_CMPA = (ka_MPA*DEPOT_MPA*F_MPA)/Vc_MPA - (Cl_MPA/Vc_MPA)*CMPA;
dxdt_DEPOT_IVIg = -0.5*DEPOT_IVIg;
dxdt_CIVG = 0.5*DEPOT_IVIg/Vc_IVIg - (Cl_IVIg/Vc_IVIg)*CIVG*(1+3*EIVIg);
dxdt_CEFG = -0.04*CEFG;
dxdt_BN   = kBN_prod - kBN_die*BN*(1+Epred_B+ERTX)
            - kBN_act*BN*Tfh_stim*(1-Treg_sup)*(1-Epred_T);
dxdt_BGC  = kBN_act*BN*Tfh_stim*(1-Treg_sup)*(1-Epred_T)
            - kBGC_die*BGC - kBM_form*BGC*(1-EMPA)*(1-Epred_B)
            - kSLPC_form*BGC*(1-EMPA) - ERTX*kBGC_die*BGC;
dxdt_BM   = kBM_form*BGC*(1-ERTX) - kBM_die*BM*(1+ERTX) - kLLPC_form*BM;
dxdt_SLPC = kSLPC_form*BGC*(1-EMPA) - kSLPC_die*SLPC*(1+Epred_B);
dxdt_LLPC = kLLPC_form*BM - kLLPC_die*LLPC;
dxdt_Tfh  = kTfh_prod*(1-Epred_T) - kTfh_die*Tfh;
dxdt_Treg = kTreg_prod + 0.3*Epred_T - kTreg_die*Treg;
dxdt_Ab3  = kAb3_prod*(SLPC + LLPC*2) - kAb3_die*Ab3*(1+E_Ab_cat);
dxdt_Ab1  = kAb1_prod*(SLPC + LLPC)   - kAb1_die*Ab1*(1+E_Ab_cat);
dxdt_Dsg3_loss  = kDsg3_int*Ab3*Dsg3_avail - kDsg3_rest*Dsg3_loss;
dxdt_PDAI_state = kPDAI_up*(Dsg3_loss/Dsg3_base) - kPDAI_dn*PDAI_state;
dxdt_Comp_act   = kComp_on*Ab3*0.25 - kComp_off*Comp_act;
dxdt_Cort_bone  = kBone_loss*Cp1_pred - kBone_rest*Cort_bone;

$TABLE
double PDAI      = PDAI_state;
double Anti_Dsg3 = Ab3;
double Anti_Dsg1 = Ab1;
double Pred_Cp   = Cp1_pred;
double RTX_Cp    = CR1;
double MPA_Cp    = CMPA;
double BMD_loss  = fmin(15.0, Cort_bone*0.5);
double CR_off    = (PDAI_state < 2.0 && Cp1_pred < 0.1) ? 1.0 : 0.0;
double BSA_bl    = fmin(100.0, 0.3*PDAI_state);
double Btotal    = BN + BGC + BM;
double PCtotal   = SLPC + LLPC;

$INIT
DEPOT_pred=0, Cp1_pred=0, Cp2_pred=0, CR1=0, CR2=0,
DEPOT_MPA=0, CMPA=0, DEPOT_IVIg=0, CIVG=0, CEFG=0,
BN=50, BGC=5, BM=20, SLPC=3, LLPC=2,
Tfh=10, Treg=8, Ab3=180, Ab1=80,
Dsg3_loss=30, PDAI_state=20, Comp_act=5, Cort_bone=0
'

mod <- mcode("PV_shiny", pv_code, quiet = TRUE)

## ─────────────────────────────────────────────────────────────────────────────
## HELPER: build event object from UI inputs
## ─────────────────────────────────────────────────────────────────────────────

build_events <- function(bw, pred_dose, use_rtx, rtx_dose, rtx_schedule,
                         use_mmf, mmf_dose, use_ivig, ivig_dose,
                         use_efg, efg_dose, sim_weeks) {

  evs <- NULL
  dur_h <- sim_weeks * 7 * 24

  # Prednisolone daily (adjusted for body weight)
  pred_mg <- pred_dose * bw  # mg/kg × kg
  if (pred_mg > 0) {
    ev_p <- ev(time = 0, amt = pred_mg, cmt = "DEPOT_pred",
               ii = 24, addl = sim_weeks * 7 - 1, rate = 0)
    # Auto-taper: halve every 8 weeks
    weeks_taper <- seq(8, min(sim_weeks, 48), by = 8)
    ev_list_pred <- list(ev_p)
    for (wk in weeks_taper) {
      pred_mg <- pred_mg / 2
      if (pred_mg < 2) break
      ev_t <- ev(time = wk * 7 * 24, amt = pred_mg, cmt = "DEPOT_pred",
                 ii = 24, addl = (sim_weeks - wk) * 7 - 1, rate = 0)
      ev_list_pred <- c(ev_list_pred, list(ev_t))
    }
    evs <- Reduce(c, ev_list_pred)
  }

  # Rituximab
  if (use_rtx && rtx_dose > 0) {
    ev_r1 <- ev(time = 0, amt = rtx_dose, cmt = "CR1", rate = rtx_dose / 6)
    if (rtx_schedule == "Standard (0, 2wk)") {
      ev_r2 <- ev(time = 2 * 7 * 24, amt = rtx_dose, cmt = "CR1", rate = rtx_dose / 6)
      evs <- if (is.null(evs)) c(ev_r1, ev_r2) else c(evs, ev_r1, ev_r2)
    } else if (rtx_schedule == "4-weekly × 4") {
      for (k in 0:3) {
        ev_rk <- ev(time = k * 4 * 7 * 24, amt = rtx_dose / 4,
                    cmt = "CR1", rate = (rtx_dose / 4) / 4)
        evs <- if (is.null(evs)) ev_rk else c(evs, ev_rk)
      }
    } else {
      # Single dose
      evs <- if (is.null(evs)) ev_r1 else c(evs, ev_r1)
    }
  }

  # MMF (BID, 500 mg per dose)
  if (use_mmf && mmf_dose > 0) {
    mmf_per_dose <- mmf_dose * 1000 / 2  # g/day → mg per BID dose
    ev_m <- ev(time = 0, amt = mmf_per_dose, cmt = "DEPOT_MPA",
               ii = 12, addl = sim_weeks * 7 * 2 - 1)
    evs <- if (is.null(evs)) ev_m else c(evs, ev_m)
  }

  # IVIg (single course over 5 days)
  if (use_ivig && ivig_dose > 0) {
    ivig_mg <- ivig_dose * bw * 1000  # g/kg → mg
    ev_i <- ev(time = 0, amt = ivig_mg, cmt = "DEPOT_IVIg",
               rate = ivig_mg / (5 * 24))
    evs <- if (is.null(evs)) ev_i else c(evs, ev_i)
  }

  # Efgartigimod (q4w IV)
  if (use_efg && efg_dose > 0) {
    efg_mg <- efg_dose * bw
    for (k in 0:floor(sim_weeks / 4)) {
      t_k <- k * 4 * 7 * 24
      if (t_k < dur_h) {
        ev_e <- ev(time = t_k, amt = efg_mg, cmt = "CEFG", rate = efg_mg / 2)
        evs <- if (is.null(evs)) ev_e else c(evs, ev_e)
      }
    }
  }

  evs
}

## ─────────────────────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "Pemphigus Vulgaris QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",      tabName = "tab_patient",  icon = icon("user")),
      menuItem("Drug PK",              tabName = "tab_pk",       icon = icon("pills")),
      menuItem("Immune Biomarkers",    tabName = "tab_immune",   icon = icon("shield-alt")),
      menuItem("Disease Activity",     tabName = "tab_disease",  icon = icon("heartbeat")),
      menuItem("Scenario Comparison",  tabName = "tab_scenario", icon = icon("chart-bar")),
      menuItem("Safety & Endpoints",   tabName = "tab_safety",   icon = icon("exclamation-triangle"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f7fb; }
      .box { border-radius: 8px; }
      .box-header { border-radius: 8px 8px 0 0; }
    "))),

    tabItems(

      ## ── TAB 1: Patient Profile ─────────────────────────────────────────
      tabItem("tab_patient",
        fluidRow(
          box(title = "Patient Demographics", width = 4, status = "primary", solidHeader = TRUE,
            sliderInput("bw",   "Body Weight (kg)", 40, 120, 70),
            selectInput("sex",  "Sex", c("Female", "Male")),
            sliderInput("age",  "Age (years)", 18, 80, 45),
            selectInput("pv_type", "PV Subtype",
              c("Mucosal PV (Dsg3++, Dsg1 low)",
                "Mucocutaneous PV (Dsg3++ + Dsg1+)",
                "Cutaneous only (pemphigus foliaceus)")),
            selectInput("severity", "Initial Severity",
              c("Mild (PDAI <15)", "Moderate (PDAI 15-45)", "Severe (PDAI >45)"))
          ),
          box(title = "Baseline Biomarkers", width = 4, status = "warning", solidHeader = TRUE,
            sliderInput("init_Ab3",  "Baseline Anti-Dsg3 IgG (U/mL)", 20, 500, 180),
            sliderInput("init_Ab1",  "Baseline Anti-Dsg1 IgG (U/mL)",  0, 300,  80),
            sliderInput("init_PDAI", "Baseline PDAI Score",  0, 100, 20),
            numericInput("BN_init",  "Naive B cells (×10⁶/L)", 50, min = 5,  max = 200),
            numericInput("BM_init",  "Memory B cells (×10⁶/L)", 20, min = 2, max = 100)
          ),
          box(title = "Simulation Settings", width = 4, status = "info", solidHeader = TRUE,
            sliderInput("sim_weeks", "Simulation Duration (weeks)", 12, 104, 52),
            selectInput("sim_delta", "Output Frequency",
              c("Daily" = "24", "Every 3 days" = "72", "Weekly" = "168")),
            hr(),
            h5("Comorbidities / Risk Modifiers"),
            checkboxInput("hx_dm",   "Type 2 Diabetes (steroid risk ↑)", FALSE),
            checkboxInput("hx_bone", "Osteoporosis at baseline", FALSE),
            checkboxInput("hx_infect", "Frequent infections", FALSE),
            hr(),
            actionButton("run_sim", "Run Simulation",
                         class = "btn-success btn-lg btn-block",
                         icon = icon("play"))
          )
        ),
        fluidRow(
          box(title = "Treatment Regimen", width = 12, status = "success", solidHeader = TRUE,
            fluidRow(
              column(4,
                h5("Corticosteroids"),
                sliderInput("pred_dose", "Prednisone (mg/kg/day)", 0, 2, 1.0, step = 0.1),
                selectInput("steroid_form", "Formulation",
                  c("Oral prednisone", "Oral prednisolone", "IV methylprednisolone pulse"))
              ),
              column(4,
                h5("Rituximab"),
                checkboxInput("use_rtx", "Use Rituximab", FALSE),
                conditionalPanel("input.use_rtx",
                  sliderInput("rtx_dose", "RTX Dose (mg/infusion)", 100, 2000, 1000, step = 100),
                  selectInput("rtx_schedule", "Dosing Schedule",
                    c("Standard (0, 2wk)", "4-weekly × 4", "Single dose"))
                )
              ),
              column(4,
                h5("Adjuvant Agents"),
                checkboxInput("use_mmf", "Mycophenolate (MMF)", FALSE),
                conditionalPanel("input.use_mmf",
                  sliderInput("mmf_dose", "MMF (g/day)", 1, 3, 2, step = 0.5)
                ),
                checkboxInput("use_ivig", "IVIg", FALSE),
                conditionalPanel("input.use_ivig",
                  sliderInput("ivig_dose", "IVIg (g/kg)", 0.5, 3, 2, step = 0.5)
                ),
                checkboxInput("use_efg", "Efgartigimod (FcRn blocker)", FALSE),
                conditionalPanel("input.use_efg",
                  sliderInput("efg_dose", "Efgartigimod (mg/kg q4w)", 5, 25, 10, step = 5)
                )
              )
            )
          )
        )
      ),

      ## ── TAB 2: Drug PK ────────────────────────────────────────────────
      tabItem("tab_pk",
        fluidRow(
          box(title = "Prednisolone Plasma Concentration", width = 6,
            plotlyOutput("plot_pk_pred", height = "350px")),
          box(title = "Rituximab Plasma Concentration", width = 6,
            plotlyOutput("plot_pk_rtx", height = "350px"))
        ),
        fluidRow(
          box(title = "MPA Plasma Concentration (MMF)", width = 6,
            plotlyOutput("plot_pk_mpa", height = "300px")),
          box(title = "IVIg / Efgartigimod Levels", width = 6,
            plotlyOutput("plot_pk_other", height = "300px"))
        )
      ),

      ## ── TAB 3: Immune Biomarkers ───────────────────────────────────────
      tabItem("tab_immune",
        fluidRow(
          box(title = "B Cell Compartment Dynamics", width = 6,
            plotlyOutput("plot_bcell", height = "380px")),
          box(title = "Plasma Cell Dynamics", width = 6,
            plotlyOutput("plot_pc", height = "380px"))
        ),
        fluidRow(
          box(title = "Tfh & Treg Cells", width = 6,
            plotlyOutput("plot_tfh_treg", height = "320px")),
          box(title = "Anti-Dsg3 vs Anti-Dsg1 IgG", width = 6,
            plotlyOutput("plot_antibody", height = "320px"))
        )
      ),

      ## ── TAB 4: Disease Activity ────────────────────────────────────────
      tabItem("tab_disease",
        fluidRow(
          box(title = "PDAI Score Trajectory", width = 8,
            plotlyOutput("plot_pdai", height = "380px")),
          box(title = "Disease Milestones", width = 4,
            valueBoxOutput("vbox_pdai_wk12",  width = 12),
            valueBoxOutput("vbox_pdai_wk24",  width = 12),
            valueBoxOutput("vbox_cr_off",     width = 12)
          )
        ),
        fluidRow(
          box(title = "Dsg3 Protein Loss", width = 6,
            plotlyOutput("plot_dsg3", height = "300px")),
          box(title = "Blister Area (BSA%)", width = 6,
            plotlyOutput("plot_bsa", height = "300px"))
        )
      ),

      ## ── TAB 5: Scenario Comparison ────────────────────────────────────
      tabItem("tab_scenario",
        fluidRow(
          box(title = "Multi-scenario PDAI Comparison", width = 12,
            plotlyOutput("plot_scenarios", height = "420px"))
        ),
        fluidRow(
          box(title = "Scenario Summary Table", width = 12,
            DTOutput("tbl_scenarios"))
        )
      ),

      ## ── TAB 6: Safety & Clinical Endpoints ───────────────────────────
      tabItem("tab_safety",
        fluidRow(
          box(title = "Steroid-Induced Bone Loss", width = 6,
            plotlyOutput("plot_bone", height = "320px")),
          box(title = "Complement Activation", width = 6,
            plotlyOutput("plot_comp", height = "320px"))
        ),
        fluidRow(
          box(title = "Safety Summary", width = 6, status = "danger",
            tableOutput("tbl_safety")),
          box(title = "Treatment Response Summary", width = 6, status = "success",
            tableOutput("tbl_response"))
        )
      )
    )
  )
)

## ─────────────────────────────────────────────────────────────────────────────
## SERVER
## ─────────────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  ## ── Reactive simulation ──────────────────────────────────────────────────
  sim_result <- eventReactive(input$run_sim, {
    req(input$bw)

    evs <- build_events(
      bw          = input$bw,
      pred_dose   = input$pred_dose,
      use_rtx     = input$use_rtx,
      rtx_dose    = if (input$use_rtx) input$rtx_dose else 0,
      rtx_schedule= if (input$use_rtx) input$rtx_schedule else "Single dose",
      use_mmf     = input$use_mmf,
      mmf_dose    = if (input$use_mmf) input$mmf_dose else 0,
      use_ivig    = input$use_ivig,
      ivig_dose   = if (input$use_ivig) input$ivig_dose else 0,
      use_efg     = input$use_efg,
      efg_dose    = if (input$use_efg) input$efg_dose else 0,
      sim_weeks   = input$sim_weeks
    )

    if (is.null(evs)) {
      showNotification("No treatment selected — running disease natural history.",
                       type = "warning")
      evs <- ev(time = 0, amt = 0, cmt = "DEPOT_pred")
    }

    init_Ab3 <- input$init_Ab3
    init_PDAI <- input$init_PDAI
    delta_h <- as.numeric(input$sim_delta)

    out <- mod %>%
      param(kAb3_prod = 0.0006) %>%
      init(Ab3 = init_Ab3, Ab1 = input$init_Ab1,
           PDAI_state = init_PDAI, Dsg3_loss = init_PDAI * 1.5,
           BN = input$BN_init, BM = input$BM_init) %>%
      ev(evs) %>%
      mrgsim(end = input$sim_weeks * 7 * 24, delta = delta_h) %>%
      as.data.frame() %>%
      mutate(week = time / 168)

    out
  })

  ## ── Scenario comparison data ─────────────────────────────────────────────
  scenario_data <- reactive({
    bw <- 70
    time_end <- 52 * 7 * 24

    run_sc <- function(label, pred_mg_day, rtx_d, rtx_sch, mmf_g, ivig_g, efg_d) {
      evs <- build_events(
        bw = bw, pred_dose = pred_mg_day / bw,
        use_rtx = rtx_d > 0, rtx_dose = rtx_d, rtx_schedule = rtx_sch,
        use_mmf = mmf_g > 0, mmf_dose = mmf_g,
        use_ivig = ivig_g > 0, ivig_dose = ivig_g / bw,
        use_efg = efg_d > 0, efg_dose = efg_d,
        sim_weeks = 52
      )
      if (is.null(evs)) evs <- ev(time = 0, amt = 0, cmt = "DEPOT_pred")
      out <- mod %>% ev(evs) %>%
        mrgsim(end = time_end, delta = 24) %>%
        as.data.frame() %>%
        mutate(week = time / 168, scenario = label)
      out
    }

    sc_list <- list(
      run_sc("High-dose CS",         105, 0,    "Single dose",  0, 0, 0),
      run_sc("RTX + Low CS",          35, 1000, "Standard (0, 2wk)", 0, 0, 0),
      run_sc("MMF + Mod CS",          70, 0,    "Single dose",  2, 0, 0),
      run_sc("RTX + MMF + Low CS",    35, 1000, "Standard (0, 2wk)", 2, 0, 0),
      run_sc("Efgartigimod + Low CS", 35, 0,    "Single dose",  0, 0, 10)
    )

    bind_rows(sc_list)
  })

  ## ── TAB 2: PK PLOTS ──────────────────────────────────────────────────────
  output$plot_pk_pred <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(week, Pred_Cp)) + geom_line(color = "#E74C3C", size = 1.2) +
      labs(x = "Week", y = "Prednisolone (mg/L)") + theme_bw()
    ggplotly(p)
  })

  output$plot_pk_rtx <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(week, RTX_Cp)) + geom_line(color = "#2980B9", size = 1.2) +
      labs(x = "Week", y = "Rituximab (mg/L)") + theme_bw()
    ggplotly(p)
  })

  output$plot_pk_mpa <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(week, MPA_Cp)) + geom_line(color = "#27AE60", size = 1.2) +
      labs(x = "Week", y = "MPA (mg/L)") + theme_bw()
    ggplotly(p)
  })

  output$plot_pk_other <- renderPlotly({
    df <- sim_result()
    df2 <- df %>% select(week, CIVG, CEFG) %>%
      pivot_longer(-week, names_to = "drug", values_to = "conc")
    p <- ggplot(df2, aes(week, conc, color = drug)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = c("CIVG" = "#F39C12", "CEFG" = "#16A085")) +
      labs(x = "Week", y = "Concentration (mg/L or mg/mL)", color = "Drug") + theme_bw()
    ggplotly(p)
  })

  ## ── TAB 3: IMMUNE BIOMARKERS ─────────────────────────────────────────────
  output$plot_bcell <- renderPlotly({
    df <- sim_result() %>%
      select(week, BN, BGC, BM) %>%
      pivot_longer(-week, names_to = "type", values_to = "count")
    p <- ggplot(df, aes(week, count, color = type)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = c("BN" = "#3498DB", "BGC" = "#9B59B6", "BM" = "#1ABC9C"),
                         labels = c("BN" = "Naive B", "BGC" = "GC B", "BM" = "Memory B")) +
      labs(x = "Week", y = "B Cells (×10⁶/L)", color = "B Cell Type") + theme_bw()
    ggplotly(p)
  })

  output$plot_pc <- renderPlotly({
    df <- sim_result() %>%
      select(week, SLPC, LLPC) %>%
      pivot_longer(-week, names_to = "type", values_to = "count")
    p <- ggplot(df, aes(week, count, color = type)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = c("SLPC" = "#E67E22", "LLPC" = "#C0392B"),
                         labels = c("SLPC" = "Short-lived PC", "LLPC" = "Long-lived PC")) +
      labs(x = "Week", y = "Plasma Cells (×10⁶/L)", color = "PC Type") + theme_bw()
    ggplotly(p)
  })

  output$plot_tfh_treg <- renderPlotly({
    df <- sim_result() %>%
      select(week, Tfh, Treg) %>%
      pivot_longer(-week, names_to = "cell", values_to = "count")
    p <- ggplot(df, aes(week, count, color = cell)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = c("Tfh" = "#8E44AD", "Treg" = "#27AE60")) +
      labs(x = "Week", y = "Cell Count (arb)", color = "Cell Type") + theme_bw()
    ggplotly(p)
  })

  output$plot_antibody <- renderPlotly({
    df <- sim_result() %>%
      select(week, Anti_Dsg3, Anti_Dsg1) %>%
      pivot_longer(-week, names_to = "ab", values_to = "titer")
    p <- ggplot(df, aes(week, titer, color = ab)) +
      geom_line(size = 1.1) +
      geom_hline(yintercept = 20, linetype = "dashed", color = "#CC0000") +
      scale_color_manual(values = c("Anti_Dsg3" = "#E74C3C", "Anti_Dsg1" = "#F39C12"),
                         labels = c("Anti_Dsg3" = "Anti-Dsg3 IgG", "Anti_Dsg1" = "Anti-Dsg1 IgG")) +
      labs(x = "Week", y = "IgG Titer (U/mL)", color = "Antibody") + theme_bw()
    ggplotly(p)
  })

  ## ── TAB 4: DISEASE ACTIVITY ──────────────────────────────────────────────
  output$plot_pdai <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(week, PDAI)) +
      geom_line(color = "#E74C3C", size = 1.3) +
      geom_hline(yintercept = 2,  linetype = "dashed", color = "#2ECC71", size = 0.9) +
      geom_hline(yintercept = 8,  linetype = "dashed", color = "#F39C12", size = 0.9) +
      annotate("text", x = max(df$week) * 0.6, y = 2.8, label = "CR-off threshold (PDAI < 2)") +
      annotate("text", x = max(df$week) * 0.6, y = 8.8, label = "PR-min threshold (PDAI < 8)") +
      labs(title = "Pemphigus Disease Activity Index (PDAI)",
           x = "Week", y = "PDAI Score (0–250)") + theme_bw()
    ggplotly(p)
  })

  output$vbox_pdai_wk12 <- renderValueBox({
    df <- sim_result()
    val <- round(df$PDAI[df$week >= 12][1], 1)
    valueBox(paste0(val), "PDAI at Week 12", icon = icon("chart-line"),
             color = if (val < 8) "green" else if (val < 20) "yellow" else "red")
  })

  output$vbox_pdai_wk24 <- renderValueBox({
    df <- sim_result()
    val <- round(df$PDAI[df$week >= 24][1], 1)
    valueBox(paste0(val), "PDAI at Week 24", icon = icon("chart-line"),
             color = if (val < 8) "green" else if (val < 20) "yellow" else "red")
  })

  output$vbox_cr_off <- renderValueBox({
    df <- sim_result()
    cr_wk52 <- df$CR_off[df$week >= 52][1]
    if (is.na(cr_wk52)) cr_wk52 <- 0
    valueBox(if (cr_wk52 == 1) "Achieved" else "Not Achieved",
             "CR-off at Week 52",
             icon = icon(if (cr_wk52 == 1) "check-circle" else "times-circle"),
             color = if (cr_wk52 == 1) "green" else "red")
  })

  output$plot_dsg3 <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(week, Dsg3_loss)) +
      geom_line(color = "#F39C12", size = 1.2) +
      labs(x = "Week", y = "Dsg3 Loss (arb units)") + theme_bw()
    ggplotly(p)
  })

  output$plot_bsa <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(week, BSA_bl)) +
      geom_area(fill = "#FADBD8", alpha = 0.6) +
      geom_line(color = "#C0392B", size = 1.2) +
      labs(x = "Week", y = "Blister Area (%BSA, estimated)") + theme_bw()
    ggplotly(p)
  })

  ## ── TAB 5: SCENARIO COMPARISON ───────────────────────────────────────────
  output$plot_scenarios <- renderPlotly({
    df <- scenario_data()
    palette <- c(
      "High-dose CS"         = "#E74C3C",
      "RTX + Low CS"         = "#2980B9",
      "MMF + Mod CS"         = "#27AE60",
      "RTX + MMF + Low CS"   = "#8E44AD",
      "Efgartigimod + Low CS"= "#16A085"
    )
    p <- ggplot(df, aes(week, PDAI, color = scenario)) +
      geom_line(size = 1.1) +
      geom_hline(yintercept = 2, linetype = "dashed", color = "black") +
      scale_color_manual(values = palette) +
      labs(title = "PDAI Comparison Across Treatment Regimens",
           x = "Week", y = "PDAI Score", color = "Scenario") + theme_bw() +
      theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$tbl_scenarios <- renderDT({
    df <- scenario_data()
    tbl <- df %>%
      group_by(scenario) %>%
      summarize(
        `PDAI W4`  = round(PDAI[week >= 4  & week < 5][1],  1),
        `PDAI W12` = round(PDAI[week >= 12 & week < 13][1], 1),
        `PDAI W24` = round(PDAI[week >= 24 & week < 25][1], 1),
        `PDAI W52` = round(PDAI[week >= 52 & week < 53][1], 1),
        `Anti-Dsg3 W52` = round(Anti_Dsg3[week >= 52 & week < 53][1], 0),
        `CR-off W52`    = ifelse(CR_off[week >= 52 & week < 53][1] == 1, "Yes", "No"),
        `BMD Loss W52 (%)` = round(BMD_loss[week >= 52 & week < 53][1], 1),
        .groups = "drop"
      )
    datatable(tbl, options = list(dom = "t", pageLength = 10),
              rownames = FALSE) %>%
      formatStyle("CR-off W52",
                  color = styleEqual(c("Yes", "No"), c("green", "red")),
                  fontWeight = "bold")
  })

  ## ── TAB 6: SAFETY & ENDPOINTS ────────────────────────────────────────────
  output$plot_bone <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(week, BMD_loss)) +
      geom_line(color = "#8E44AD", size = 1.2) +
      geom_hline(yintercept = 10, linetype = "dashed", color = "#CC0000") +
      annotate("text", x = max(df$week) * 0.5, y = 10.5, label = "Osteoporosis risk threshold") +
      labs(x = "Week", y = "Estimated BMD Loss (%)") + theme_bw()
    ggplotly(p)
  })

  output$plot_comp <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(week, Comp_act)) +
      geom_line(color = "#E74C3C", size = 1.1) +
      labs(x = "Week", y = "Complement Activity (arb)") + theme_bw()
    ggplotly(p)
  })

  output$tbl_safety <- renderTable({
    df <- sim_result()
    data.frame(
      Metric           = c("Peak Pred Cp (mg/L)", "Cumul. Pred AUC (mg/L·wk)",
                           "Steroid BMD loss (%)", "RTX peak Cp (mg/L)",
                           "Max complement activation"),
      Value            = c(round(max(df$Pred_Cp, na.rm = TRUE), 2),
                           round(sum(df$Pred_Cp) * 7, 0),
                           round(max(df$BMD_loss, na.rm = TRUE), 2),
                           round(max(df$RTX_Cp,  na.rm = TRUE), 2),
                           round(max(df$Comp_act, na.rm = TRUE), 2))
    )
  })

  output$tbl_response <- renderTable({
    df <- sim_result()
    wk12_row <- df[df$week >= 12 & df$week < 13, ][1, ]
    wk24_row <- df[df$week >= 24 & df$week < 25, ][1, ]
    wk52_row <- df[df$week >= min(df$week[df$week >= 52]), ][1, ]
    data.frame(
      Timepoint        = c("Week 12", "Week 24", "Week 52"),
      `PDAI Score`     = c(round(wk12_row$PDAI, 1),
                           round(wk24_row$PDAI, 1),
                           round(wk52_row$PDAI, 1)),
      `Anti-Dsg3 (U/mL)` = c(round(wk12_row$Anti_Dsg3, 0),
                              round(wk24_row$Anti_Dsg3, 0),
                              round(wk52_row$Anti_Dsg3, 0)),
      `BSA Blisters (%)` = c(round(wk12_row$BSA_bl, 1),
                              round(wk24_row$BSA_bl, 1),
                              round(wk52_row$BSA_bl, 1)),
      `CR-off`         = c(ifelse(wk12_row$CR_off == 1, "Yes", "No"),
                           ifelse(wk24_row$CR_off == 1, "Yes", "No"),
                           ifelse(wk52_row$CR_off == 1, "Yes", "No")),
      check.names = FALSE
    )
  })
}

## ─────────────────────────────────────────────────────────────────────────────
## LAUNCH
## ─────────────────────────────────────────────────────────────────────────────

shinyApp(ui, server)
