## ============================================================
## Relapsing Polychondritis (RP) — Interactive Shiny QSP App
## ============================================================
## Tabs:
##   1. Patient Profile & Disease Overview
##   2. Pharmacokinetics (PK)
##   3. Immune & Cytokine Dynamics (PD)
##   4. Cartilage Integrity & RPDAI
##   5. Treatment Scenario Comparison
##   6. Biomarker Trajectories
##   7. Sensitivity Analysis
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(tidyr)
library(DT)

## ── Inline mrgsolve model ────────────────────────────────────
rpc_code <- '
$PROB RP QSP Shiny Model

$PARAM
ka_pred=1.5, CL_pred=15.0, V1_pred=30.0, V2_pred=45.0, Q_pred=8.0, F_pred=0.82,
kon_GR=0.20, koff_GR=0.02, ksyn_GR=10.0, kdeg_GR=0.05, MW_pred=360.4,
CL_tcz=0.008, V1_tcz=4.0, V2_tcz=6.5, Q_tcz=0.015,
kin_Tact=0.50, kout_Tact=0.40, kin_Th17=0.30, kout_Th17=0.35,
kin_Treg=0.25, kout_Treg=0.30,
kin_Bact=0.20, kout_Bact=0.18, kprod_Ab=0.05, kdeg_Ab=0.004,
kform_IC=0.10, kelim_IC=0.15, IC50_IC=5.0,
kact_Comp=0.40, kdeg_Comp=0.60,
kin_TNF=0.30, kout_TNF=1.20,
kin_IL6=0.25, kout_IL6=0.80,
kin_IL17=0.30, kout_IL17=0.90,
kin_IL1b=0.35, kout_IL1b=1.00,
kin_MMP=0.15, kout_MMP=0.40,
kdest_Cart=0.06, krep_Cart=0.015, Cart_max=100.0,
Emax_pred=0.85, EC50_pred=8.0, nH_pred=2.0,
Emax_tcz=0.95, EC50_tcz=0.30,
E_MTX=0.0, E_ABA=0.0, E_DAP=0.0,
h_IL6_Tact=0.30, h_IC_Tact=0.25, h_Treg_Tact=0.60,
h_IL6_Th17=0.80, h_Treg_Th17=0.70,
h_Tact_Bact=0.50, h_TNF_MMP=0.50, h_IL17_MMP=0.40, h_IL1b_MMP=0.30,
h_Comp_MMP=0.20, h_Tact_TNF=0.60, h_Tact_IL6=0.50,
h_TNF_IL6=0.40, h_IL1b_IL6=0.30, h_Th17_IL17=1.00,
h_IC_IL1b=0.40, h_Comp_IL1b=0.30

$CMT
DEPOT_pred C1_pred C2_pred GR_free GR_bound C1_tcz C2_tcz
Tact Th17 Treg Bact Ab_CII IC Comp
TNF IL6 IL17 IL1b MMP Cart

$INIT
DEPOT_pred=0, C1_pred=0, C2_pred=0,
GR_free=200, GR_bound=0, C1_tcz=0, C2_tcz=0,
Tact=1.25, Th17=0.857, Treg=0.833, Bact=1.111,
Ab_CII=50, IC=5, Comp=3.33,
TNF=0.25, IL6=0.3125, IL17=0.333, IL1b=0.35,
MMP=0.375, Cart=100

$MAIN
double Cp_pred    = C1_pred / V1_pred;
double Cp_pred_nM = Cp_pred * 1000.0 / MW_pred;
double Cp_tcz     = C1_tcz  / V1_tcz;
double Epred = Emax_pred * pow(Cp_pred, nH_pred) /
               (pow(EC50_pred, nH_pred) + pow(Cp_pred, nH_pred));
double Etcz  = Emax_tcz * Cp_tcz / (EC50_tcz + Cp_tcz);
double E_IL6_tot = fmax(Epred * 0.7, Etcz);
double GR_occ    = GR_bound / (GR_free + GR_bound + 0.001);
if(Cart < 0.0) Cart = 0.0;

$ODE
dxdt_DEPOT_pred = -ka_pred * DEPOT_pred;
dxdt_C1_pred    = F_pred*ka_pred*DEPOT_pred - (CL_pred+Q_pred)/V1_pred*C1_pred + Q_pred/V2_pred*C2_pred;
dxdt_C2_pred    = Q_pred/V1_pred*C1_pred - Q_pred/V2_pred*C2_pred;
dxdt_GR_free    = ksyn_GR - kdeg_GR*GR_free - kon_GR*Cp_pred_nM*GR_free + koff_GR*GR_bound;
dxdt_GR_bound   = kon_GR*Cp_pred_nM*GR_free - koff_GR*GR_bound - kdeg_GR*GR_bound;
dxdt_C1_tcz     = -(CL_tcz+Q_tcz)/V1_tcz*C1_tcz + Q_tcz/V2_tcz*C2_tcz;
dxdt_C2_tcz     = Q_tcz/V1_tcz*C1_tcz - Q_tcz/V2_tcz*C2_tcz;

double drive_Tact = 1.0 + h_IL6_Tact*IL6 + h_IC_Tact*IC;
double supp_Tact  = 1.0 + h_Treg_Tact*Treg + Epred*0.5 + E_ABA;
dxdt_Tact = kin_Tact*drive_Tact/supp_Tact - kout_Tact*Tact;

double drive_Th17 = 1.0 + h_IL6_Th17*IL6 + 0.3*IL1b;
double supp_Th17  = 1.0 + h_Treg_Th17*Treg + E_MTX*1.5 + Epred*0.4;
dxdt_Th17 = kin_Th17*drive_Th17/supp_Th17 - kout_Th17*Th17;

dxdt_Treg = kin_Treg*(1.0 + 0.3*GR_occ) - kout_Treg*Treg;

double drive_Bact = 1.0 + h_Tact_Bact*Tact;
dxdt_Bact = kin_Bact*drive_Bact*(1.0 - Epred*0.3 - E_MTX*0.4) - kout_Bact*Bact;
dxdt_Ab_CII = kprod_Ab*Bact - kdeg_Ab*Ab_CII;
dxdt_IC     = kform_IC*Ab_CII/(IC50_IC + Ab_CII)*Ab_CII - kelim_IC*IC;
dxdt_Comp   = kact_Comp*IC - kdeg_Comp*Comp;

double prod_TNF  = kin_TNF*(1.0 + h_Tact_TNF*Tact + 0.2*Comp);
dxdt_TNF = prod_TNF*(1.0 - Epred) - kout_TNF*TNF;
double prod_IL6  = kin_IL6*(1.0 + h_Tact_IL6*Tact + h_TNF_IL6*TNF + h_IL1b_IL6*IL1b);
dxdt_IL6 = prod_IL6*(1.0 - E_IL6_tot) - kout_IL6*IL6;
double prod_IL17 = kin_IL17*(1.0 + h_Th17_IL17*Th17);
dxdt_IL17 = prod_IL17*(1.0 - E_MTX*0.5) - kout_IL17*IL17;
double prod_IL1b = kin_IL1b*(1.0 + h_IC_IL1b*IC + h_Comp_IL1b*Comp + 0.2*TNF);
dxdt_IL1b = prod_IL1b*(1.0 - Epred*0.7) - kout_IL1b*IL1b;

double prod_MMP = kin_MMP*(1.0 + h_TNF_MMP*TNF + h_IL17_MMP*IL17 + h_IL1b_MMP*IL1b + h_Comp_MMP*Comp);
dxdt_MMP = prod_MMP*(1.0 - Epred*0.4) - kout_MMP*MMP;
dxdt_Cart = krep_Cart*(Cart_max - Cart) - kdest_Cart*MMP*(Cart/Cart_max);

$TABLE
capture Cp_pred_mgL = C1_pred/V1_pred;
capture Cp_tcz_mgL  = C1_tcz/V1_tcz;
capture GR_occ_pct  = GR_bound/(GR_free+GR_bound+0.001)*100;
capture E_pred_pct  = Epred*100;
capture E_tcz_pct   = Etcz*100;
capture CartPct     = Cart;
capture RPDAI_proxy = 10*(TNF+IL6+IL17) + 50*(1 - Cart/100);
capture CRP_proxy   = 0.5*IL6 + 0.2*TNF;
'

mod_base <- mcode("rpc_shiny", rpc_code, quiet = TRUE)

## ── Helper: run simulation ──────────────────────────────────
run_sim <- function(params_list, pred_dose = 0, pred_ii = 24,
                    tcz_dose = 0, sim_days = 365) {
  mod_s <- param(mod_base, params_list)
  events <- ev()
  if (pred_dose > 0)
    events <- c(events, ev(amt = pred_dose, ii = pred_ii,
                           addl = round(sim_days * 24 / pred_ii),
                           cmt = 1))
  if (tcz_dose > 0)
    events <- c(events, ev(amt = tcz_dose, ii = 672,
                           addl = round(sim_days / 28),
                           cmt = 6))
  out <- mrgsim(mod_s, events = events,
                end = sim_days * 24, delta = 24, digits = 4)
  as.data.frame(out) %>% mutate(time_d = time / 24)
}

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(
    title = "RP QSP Model",
    titleWidth = 250
  ),
  dashboardSidebar(
    width = 250,
    sidebarMenu(
      menuItem("Patient Profile",         tabName = "tab_profile",  icon = icon("user-md")),
      menuItem("Pharmacokinetics (PK)",   tabName = "tab_pk",       icon = icon("flask")),
      menuItem("Immune & Cytokine (PD)",  tabName = "tab_pd",       icon = icon("dna")),
      menuItem("Cartilage & RPDAI",       tabName = "tab_cart",     icon = icon("bone")),
      menuItem("Scenario Comparison",     tabName = "tab_scen",     icon = icon("chart-bar")),
      menuItem("Biomarker Trajectories",  tabName = "tab_bio",      icon = icon("vial")),
      menuItem("Sensitivity Analysis",    tabName = "tab_sens",     icon = icon("sliders-h"))
    ),
    hr(),
    h5("Simulation Settings", style = "padding-left:15px;color:#ccc;"),
    sliderInput("sim_days",  "Duration (days)",   min = 30, max = 730, value = 365, step = 30),
    sliderInput("pred_dose", "Prednisolone (mg/d)", min = 0, max = 80, value = 40, step = 5),
    sliderInput("tcz_dose",  "Tocilizumab (mg IV q4w)", min = 0, max = 800, value = 0, step = 80),
    checkboxInput("use_mtx", "Add Methotrexate (15mg/wk)", value = FALSE),
    checkboxInput("use_aba", "Add Abatacept",               value = FALSE),
    checkboxInput("use_dap", "Dapsone (100mg/d)",           value = FALSE),
    actionButton("run_sim",  "Run Simulation", icon = icon("play"),
                 class = "btn-primary btn-block",
                 style = "margin: 10px 15px; width: 220px;")
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .skin-purple .main-header .logo { background-color: #4a2c6e; }
      .skin-purple .main-header .navbar { background-color: #5c3384; }
      .content-wrapper { background-color: #f8f8f8; }
    "))),
    tabItems(

      ## ── Tab 1: Patient Profile ──────────────────────────────
      tabItem(tabName = "tab_profile",
        fluidRow(
          box(width = 12, title = "Relapsing Polychondritis — Disease Overview",
              status = "purple", solidHeader = TRUE,
            fluidRow(
              column(6,
                h4("Key Pathophysiology"),
                p("Relapsing polychondritis (RP) is a rare autoimmune disease
                  characterised by recurrent episodes of inflammation targeting
                  cartilaginous structures — particularly auricular, nasal, and
                  tracheobronchial cartilage."),
                p("Central autoantigens: Type II collagen, Matrilin-1, Type IX & XI
                  collagen, and COMP trigger anti-CII antibody production."),
                tags$ul(
                  tags$li("CD4+ Th1 / Th17 cells → cytokine storm (TNF-α, IL-6, IL-17A)"),
                  tags$li("Anti-CII IgG → immune complex deposition → complement activation"),
                  tags$li("Complement MAC + ROS → chondrocyte apoptosis"),
                  tags$li("MMPs (MMP-1/3/9/13) + ADAMTS-4/5 → collagen/PG degradation"),
                  tags$li("Saddle-nose deformity, cauliflower ear, tracheomalacia, SNHL")
                )
              ),
              column(6,
                h4("Disease Activity Index (RPDAI)"),
                tableOutput("rpdai_table"),
                hr(),
                h4("Epidemiology"),
                tags$ul(
                  tags$li("Incidence: ~3.5 per million/year"),
                  tags$li("Peak onset: 40–60 years"),
                  tags$li("M:F ratio ~1:1"),
                  tags$li("Mortality risk: airway collapse, cardiac, amyloidosis"),
                  tags$li("5-year survival ~74% without therapy")
                )
              )
            )
          )
        ),
        fluidRow(
          valueBoxOutput("vb_cart",  width = 3),
          valueBoxOutput("vb_rpdai", width = 3),
          valueBoxOutput("vb_crp",   width = 3),
          valueBoxOutput("vb_abcii", width = 3)
        ),
        fluidRow(
          box(width = 12, title = "Patient Parameter Settings",
              status = "info", solidHeader = TRUE,
            fluidRow(
              column(4,
                sliderInput("body_wt", "Body weight (kg)", 40, 120, 70),
                selectInput("hla_type", "HLA Type",
                            choices = c("HLA-DR4" = "dr4",
                                        "HLA-B27" = "b27",
                                        "Unknown"  = "unk"),
                            selected = "dr4"),
                sliderInput("disease_dur", "Disease duration (yr)", 0, 20, 2)
              ),
              column(4,
                sliderInput("init_abcii", "Baseline Anti-CII (nM)", 10, 200, 50),
                sliderInput("init_cart", "Baseline Cartilage (%)", 40, 100, 100),
                sliderInput("kdest_slide", "MMP Destructivity (×)", 0.5, 3.0, 1.0, step = 0.1)
              ),
              column(4,
                selectInput("organ_inv", "Primary Organ Involvement",
                            choices = c("Auricular/Nasal" = "ear",
                                        "Tracheobronchial" = "airway",
                                        "Systemic/Multiple" = "sys"),
                            selected = "ear"),
                radioButtons("severity", "Disease Severity",
                             choices = c("Mild" = "mild",
                                         "Moderate" = "mod",
                                         "Severe" = "sev"),
                             selected = "mod", inline = TRUE)
              )
            )
          )
        )
      ),

      ## ── Tab 2: PK ───────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(width = 12, title = "Pharmacokinetics", status = "primary", solidHeader = TRUE,
            tabsetPanel(
              tabPanel("Prednisolone PK",
                fluidRow(
                  column(3,
                    h5("Dosing Regimen"),
                    selectInput("pred_regimen", "Regimen",
                      choices = c("Single dose" = "single",
                                  "QD (once daily)" = "qd",
                                  "Induction → Taper" = "taper"),
                      selected = "qd"),
                    numericInput("pred_amt", "Dose (mg)", 40, min = 2.5, max = 80, step = 2.5),
                    hr(),
                    h5("PK Parameters"),
                    sliderInput("CL_pred_s", "CL (L/h)", 5, 40, 15),
                    sliderInput("V1_pred_s", "V1 (L)",   10, 80, 30)
                  ),
                  column(9,
                    plotlyOutput("pk_plot", height = "350px"),
                    br(),
                    plotlyOutput("gr_plot", height = "250px")
                  )
                )
              ),
              tabPanel("Tocilizumab PK",
                fluidRow(
                  column(3,
                    numericInput("tcz_amt_s", "TCZ Dose (mg)", 560, 100, 800, 80),
                    selectInput("tcz_ii", "Dosing Interval",
                                choices = c("q2w (336h)" = 336, "q4w (672h)" = 672),
                                selected = 672),
                    h5("PK Parameters"),
                    sliderInput("CL_tcz_s", "CL (L/h)",  0.001, 0.05, 0.008, step = 0.001),
                    sliderInput("V1_tcz_s", "V1 (L)",    1, 12, 4)
                  ),
                  column(9,
                    plotlyOutput("tcz_pk_plot", height = "400px")
                  )
                )
              ),
              tabPanel("PK Summary Table",
                DTOutput("pk_table")
              )
            )
          )
        )
      ),

      ## ── Tab 3: Immune & Cytokine PD ─────────────────────────
      tabItem(tabName = "tab_pd",
        fluidRow(
          box(width = 6, title = "T Cell Dynamics", status = "success", solidHeader = TRUE,
            plotlyOutput("tcell_plot", height = "350px")
          ),
          box(width = 6, title = "Cytokine Levels", status = "warning", solidHeader = TRUE,
            plotlyOutput("cytokine_plot", height = "350px")
          )
        ),
        fluidRow(
          box(width = 6, title = "B Cell & Antibody", status = "info", solidHeader = TRUE,
            plotlyOutput("bcell_plot", height = "300px")
          ),
          box(width = 6, title = "Complement Activation & IC", status = "danger", solidHeader = TRUE,
            plotlyOutput("complement_plot", height = "300px")
          )
        )
      ),

      ## ── Tab 4: Cartilage & RPDAI ────────────────────────────
      tabItem(tabName = "tab_cart",
        fluidRow(
          box(width = 8, title = "Cartilage Integrity Over Time",
              status = "purple", solidHeader = TRUE,
            plotlyOutput("cart_plot", height = "380px")
          ),
          box(width = 4, title = "RPDAI Proxy", status = "danger", solidHeader = TRUE,
            plotlyOutput("rpdai_plot", height = "380px")
          )
        ),
        fluidRow(
          box(width = 6, title = "MMP Activity", status = "warning", solidHeader = TRUE,
            plotlyOutput("mmp_plot", height = "280px")
          ),
          box(width = 6, title = "Organ Involvement Risk Score",
              status = "info", solidHeader = TRUE,
            plotlyOutput("organ_plot", height = "280px")
          )
        )
      ),

      ## ── Tab 5: Scenario Comparison ──────────────────────────
      tabItem(tabName = "tab_scen",
        fluidRow(
          box(width = 12, title = "Treatment Scenario Comparison",
              status = "primary", solidHeader = TRUE,
            fluidRow(
              column(3,
                h5("Scenarios to Compare"),
                checkboxGroupInput("scenarios_sel",
                  label = NULL,
                  choices = c("Untreated"            = "none",
                              "Prednisone 40mg/d"    = "pred",
                              "Pred + MTX"           = "pred_mtx",
                              "Tocilizumab"          = "tcz",
                              "Pred + TCZ"           = "pred_tcz",
                              "Abatacept"            = "aba",
                              "Dapsone"              = "dap"),
                  selected = c("none", "pred", "pred_mtx", "tcz")
                ),
                selectInput("scen_outcome", "Outcome to Plot",
                  choices = c("Cartilage Integrity" = "CartPct",
                              "RPDAI Proxy"         = "RPDAI_proxy",
                              "IL-6 (ng/mL)"        = "IL6",
                              "TNF-α (ng/mL)"       = "TNF",
                              "Anti-CII Ab (nM)"    = "Ab_CII",
                              "CRP Proxy"           = "CRP_proxy"),
                  selected = "CartPct"
                )
              ),
              column(9,
                plotlyOutput("scen_plot", height = "450px")
              )
            ),
            hr(),
            DTOutput("scen_table")
          )
        )
      ),

      ## ── Tab 6: Biomarker Trajectories ───────────────────────
      tabItem(tabName = "tab_bio",
        fluidRow(
          box(width = 12, title = "Biomarker Dashboard",
              status = "info", solidHeader = TRUE,
            fluidRow(
              column(3,
                checkboxGroupInput("bio_sel", "Select Biomarkers",
                  choices = c("CRP (proxy)"            = "CRP_proxy",
                              "Anti-CII Ab (nM)"       = "Ab_CII",
                              "IL-6 (ng/mL)"           = "IL6",
                              "TNF-α"                  = "TNF",
                              "IL-17A"                 = "IL17",
                              "IL-1β"                  = "IL1b",
                              "Immune Complex"         = "IC",
                              "Complement"             = "Comp",
                              "MMP Activity"           = "MMP",
                              "Cartilage Integrity"    = "CartPct"),
                  selected = c("CRP_proxy", "IL6", "TNF", "CartPct")
                )
              ),
              column(9,
                plotlyOutput("bio_plot", height = "500px")
              )
            )
          )
        ),
        fluidRow(
          box(width = 12, title = "Simulated Lab Values at Key Timepoints",
              status = "success", solidHeader = TRUE,
            DTOutput("bio_table")
          )
        )
      ),

      ## ── Tab 7: Sensitivity Analysis ─────────────────────────
      tabItem(tabName = "tab_sens",
        fluidRow(
          box(width = 12, title = "Parameter Sensitivity Analysis",
              status = "warning", solidHeader = TRUE,
            fluidRow(
              column(3,
                selectInput("sens_param", "Parameter",
                  choices = c("MMP destruction rate"   = "kdest_Cart",
                              "Cartilage repair rate"  = "krep_Cart",
                              "TNF production rate"    = "kin_TNF",
                              "IL-17 production rate"  = "kin_IL17",
                              "TNF→MMP coupling"       = "h_TNF_MMP",
                              "IL-17→MMP coupling"     = "h_IL17_MMP",
                              "Ab production rate"     = "kprod_Ab",
                              "Ab degradation rate"    = "kdeg_Ab",
                              "Tact input"             = "kin_Tact",
                              "Th17 input"             = "kin_Th17"),
                  selected = "kdest_Cart"
                ),
                sliderInput("sens_range", "Parameter Range (× baseline)",
                            0.1, 5.0, c(0.5, 2.0), step = 0.1),
                selectInput("sens_outcome", "Outcome",
                  choices = c("Cartilage Integrity" = "CartPct",
                              "RPDAI Proxy"         = "RPDAI_proxy",
                              "Anti-CII Ab"         = "Ab_CII"),
                  selected = "CartPct"
                ),
                numericInput("sens_day", "Evaluation day", 365, 30, 730, 30),
                actionButton("run_sens", "Run Sensitivity",
                             icon = icon("chart-line"), class = "btn-warning btn-block")
              ),
              column(9,
                plotlyOutput("sens_plot", height = "450px")
              )
            )
          )
        )
      )
    )
  )
)

## ============================================================
## SERVER
## ============================================================
server <- function(input, output, session) {

  ## ── Reactive: run main simulation ──────────────────────────
  sim_out <- eventReactive(input$run_sim, {
    params <- list(
      E_MTX = if (input$use_mtx) 0.60 else 0.0,
      E_ABA = if (input$use_aba) 0.65 else 0.0,
      E_DAP = if (input$use_dap) 0.55 else 0.0
    )
    run_sim(params,
            pred_dose = input$pred_dose,
            tcz_dose  = input$tcz_dose,
            sim_days  = input$sim_days)
  }, ignoreNULL = FALSE)

  ## ── Tab 1: value boxes ─────────────────────────────────────
  output$vb_cart <- renderValueBox({
    d <- sim_out()
    val <- round(tail(d$CartPct, 1), 1)
    col <- if (val > 80) "green" else if (val > 60) "yellow" else "red"
    valueBox(paste0(val, "%"), "Cartilage Integrity", icon = icon("bone"), color = col)
  })
  output$vb_rpdai <- renderValueBox({
    d <- sim_out()
    val <- round(tail(d$RPDAI_proxy, 1), 1)
    col <- if (val < 20) "green" else if (val < 50) "yellow" else "red"
    valueBox(val, "RPDAI Proxy", icon = icon("activity"), color = col)
  })
  output$vb_crp <- renderValueBox({
    d <- sim_out()
    val <- round(tail(d$CRP_proxy, 1), 2)
    col <- if (val < 0.5) "green" else if (val < 1.5) "yellow" else "red"
    valueBox(paste0(val, " AU"), "CRP Proxy", icon = icon("tint"), color = col)
  })
  output$vb_abcii <- renderValueBox({
    d <- sim_out()
    val <- round(tail(d$Ab_CII, 1), 1)
    col <- if (val < 30) "green" else if (val < 80) "yellow" else "red"
    valueBox(paste0(val, " nM"), "Anti-CII Ab", icon = icon("shield-alt"), color = col)
  })

  output$rpdai_table <- renderTable({
    data.frame(
      Domain = c("Auricular chondritis", "Nasal chondritis", "Laryngotracheal",
                 "Ocular inflammation", "Arthritis", "Vestibulo-cochlear",
                 "Cutaneous", "Cardiovascular", "Constitutional"),
      MaxScore = c(6, 6, 8, 6, 8, 8, 6, 6, 3)
    )
  })

  ## ── Tab 2: PK plots ────────────────────────────────────────
  output$pk_plot <- renderPlotly({
    params <- list(CL_pred = input$CL_pred_s, V1_pred = input$V1_pred_s)
    ii_h   <- if (input$pred_regimen == "single") 1e6 else 24
    addl_h <- if (input$pred_regimen == "single") 0   else round(14 * 24 / 24)
    ev_pk  <- ev(amt = input$pred_amt, ii = ii_h, addl = addl_h, cmt = 1)
    out_pk <- mrgsim(param(mod_base, params), events = ev_pk,
                     end = 14 * 24, delta = 1, quiet = TRUE)
    d <- as.data.frame(out_pk) %>% mutate(time_d = time / 24)
    plot_ly(d, x = ~time_d, y = ~Cp_pred_mgL, type = "scatter", mode = "lines",
            line = list(color = "steelblue", width = 2)) %>%
      layout(title = "Prednisolone Plasma Concentration",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Cp (mg/L)"))
  })

  output$gr_plot <- renderPlotly({
    params <- list(CL_pred = input$CL_pred_s, V1_pred = input$V1_pred_s)
    ev_pk  <- ev(amt = input$pred_amt, ii = 24, addl = 13, cmt = 1)
    out_pk <- mrgsim(param(mod_base, params), events = ev_pk,
                     end = 14 * 24, delta = 1, quiet = TRUE)
    d <- as.data.frame(out_pk) %>% mutate(time_d = time / 24)
    plot_ly(d, x = ~time_d, y = ~GR_occ_pct, type = "scatter", mode = "lines",
            line = list(color = "darkred", width = 2)) %>%
      layout(title = "GR Occupancy",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "GR Occupancy (%)"))
  })

  output$tcz_pk_plot <- renderPlotly({
    params <- list(CL_tcz = input$CL_tcz_s, V1_tcz = input$V1_tcz_s)
    ii_val <- as.numeric(input$tcz_ii)
    ev_tcz <- ev(amt = input$tcz_amt_s, ii = ii_val, addl = 5, cmt = 6)
    out_tk <- mrgsim(param(mod_base, params), events = ev_tcz,
                     end = round(ii_val * 6), delta = 6, quiet = TRUE)
    d <- as.data.frame(out_tk) %>% mutate(time_d = time / 24)
    plot_ly(d, x = ~time_d, y = ~Cp_tcz_mgL, type = "scatter", mode = "lines",
            line = list(color = "darkorange", width = 2)) %>%
      layout(title = "Tocilizumab Plasma Concentration",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Cp (mg/L)"))
  })

  output$pk_table <- renderDT({
    d <- sim_out() %>%
      filter(time_d %in% c(1, 7, 14, 30, 90, 180, 365)) %>%
      select(time_d, Cp_pred_mgL, Cp_tcz_mgL, GR_occ_pct, E_pred_pct, E_tcz_pct) %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
    datatable(d, options = list(pageLength = 10), rownames = FALSE,
              colnames = c("Day", "Pred (mg/L)", "TCZ (mg/L)",
                           "GR Occ (%)", "Pred Effect (%)", "TCZ Effect (%)"))
  })

  ## ── Tab 3: PD plots ────────────────────────────────────────
  output$tcell_plot <- renderPlotly({
    d <- sim_out() %>%
      select(time_d, Tact, Th17, Treg) %>%
      pivot_longer(-time_d, names_to = "Cell", values_to = "Count")
    plot_ly(d, x = ~time_d, y = ~Count, color = ~Cell, type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Days"), yaxis = list(title = "Cells (×10⁶)"),
             legend = list(orientation = "h"))
  })

  output$cytokine_plot <- renderPlotly({
    d <- sim_out() %>%
      select(time_d, TNF, IL6, IL17, IL1b) %>%
      pivot_longer(-time_d, names_to = "Cytokine", values_to = "Conc")
    plot_ly(d, x = ~time_d, y = ~Conc, color = ~Cytokine, type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Days"), yaxis = list(title = "Conc (ng/mL)"),
             legend = list(orientation = "h"))
  })

  output$bcell_plot <- renderPlotly({
    d <- sim_out() %>% select(time_d, Bact, Ab_CII)
    p1 <- plot_ly(d, x = ~time_d, y = ~Bact, type = "scatter", mode = "lines",
                  name = "B cells", line = list(color = "purple"))
    p2 <- plot_ly(d, x = ~time_d, y = ~Ab_CII, type = "scatter", mode = "lines",
                  name = "Anti-CII Ab (nM)", line = list(color = "magenta"), yaxis = "y2")
    subplot(p1, p2, nrows = 2, shareX = TRUE, titleY = TRUE) %>%
      layout(xaxis = list(title = "Days"))
  })

  output$complement_plot <- renderPlotly({
    d <- sim_out() %>% select(time_d, IC, Comp) %>%
      pivot_longer(-time_d, names_to = "Var", values_to = "Val")
    plot_ly(d, x = ~time_d, y = ~Val, color = ~Var, type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Days"), yaxis = list(title = "AU"))
  })

  ## ── Tab 4: Cartilage plots ─────────────────────────────────
  output$cart_plot <- renderPlotly({
    d <- sim_out()
    plot_ly(d, x = ~time_d, y = ~CartPct, type = "scatter", mode = "lines",
            line = list(color = "darkorchid", width = 2.5)) %>%
      layout(xaxis = list(title = "Days"), yaxis = list(title = "Cartilage Integrity (%)", range = c(0, 105)))
  })

  output$rpdai_plot <- renderPlotly({
    d <- sim_out()
    plot_ly(d, x = ~time_d, y = ~RPDAI_proxy, type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(255,80,80,0.2)",
            line = list(color = "red", width = 2)) %>%
      layout(xaxis = list(title = "Days"), yaxis = list(title = "RPDAI Proxy"))
  })

  output$mmp_plot <- renderPlotly({
    d <- sim_out()
    plot_ly(d, x = ~time_d, y = ~MMP, type = "scatter", mode = "lines",
            line = list(color = "darkorange", width = 2)) %>%
      layout(xaxis = list(title = "Days"), yaxis = list(title = "MMP Activity (AU)"))
  })

  output$organ_plot <- renderPlotly({
    d <- sim_out()
    d <- d %>% mutate(
      Auricular = 10 * (1 - CartPct / 100) * 6,
      Tracheal  = 10 * MMP * 0.8,
      Systemic  = CRP_proxy * 10
    ) %>% select(time_d, Auricular, Tracheal, Systemic) %>%
      pivot_longer(-time_d, names_to = "Site", values_to = "Score")
    plot_ly(d, x = ~time_d, y = ~Score, color = ~Site, type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Days"), yaxis = list(title = "Risk Score (AU)"))
  })

  ## ── Tab 5: Scenario comparison ─────────────────────────────
  scen_data <- eventReactive(input$run_sim, {
    scen_map <- list(
      none     = list(pred_dose = 0,  tcz_dose = 0,  extra = list()),
      pred     = list(pred_dose = 40, tcz_dose = 0,  extra = list()),
      pred_mtx = list(pred_dose = 40, tcz_dose = 0,  extra = list(E_MTX = 0.60)),
      tcz      = list(pred_dose = 0,  tcz_dose = 560, extra = list()),
      pred_tcz = list(pred_dose = 40, tcz_dose = 560, extra = list()),
      aba      = list(pred_dose = 0,  tcz_dose = 0,  extra = list(E_ABA = 0.65)),
      dap      = list(pred_dose = 0,  tcz_dose = 0,  extra = list(E_DAP = 0.55))
    )
    scen_labels <- c(none = "Untreated", pred = "Prednisone 40mg",
                     pred_mtx = "Pred + MTX", tcz = "Tocilizumab",
                     pred_tcz = "Pred + TCZ", aba = "Abatacept", dap = "Dapsone")
    bind_rows(lapply(input$scenarios_sel, function(sc) {
      cfg <- scen_map[[sc]]
      run_sim(cfg$extra, cfg$pred_dose, 24, cfg$tcz_dose, input$sim_days) %>%
        mutate(scenario = scen_labels[[sc]])
    }))
  })

  output$scen_plot <- renderPlotly({
    d <- scen_data()
    y_col <- input$scen_outcome
    plot_ly(d, x = ~time_d, y = as.formula(paste0("~", y_col)),
            color = ~scenario, type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Days"),
             yaxis = list(title = y_col),
             legend = list(orientation = "h"))
  })

  output$scen_table <- renderDT({
    d <- scen_data() %>%
      filter(abs(time_d - min(c(input$sim_days, 365))) < 1.5) %>%
      group_by(scenario) %>% slice(1) %>%
      select(scenario, CartPct, RPDAI_proxy, CRP_proxy, Ab_CII, TNF, IL6, IL17) %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
    datatable(d, rownames = FALSE,
              colnames = c("Scenario", "Cart (%)", "RPDAI", "CRP", "Anti-CII (nM)",
                           "TNF", "IL-6", "IL-17"))
  })

  ## ── Tab 6: Biomarker trajectories ──────────────────────────
  output$bio_plot <- renderPlotly({
    d <- sim_out() %>%
      select(time_d, CRP_proxy, Ab_CII, IL6, TNF, IL17, IL1b,
             IC, Comp, MMP, CartPct)
    sel <- input$bio_sel
    if (is.null(sel) || length(sel) == 0) sel <- "CRP_proxy"
    d_long <- d %>% select(time_d, all_of(sel)) %>%
      pivot_longer(-time_d, names_to = "Biomarker", values_to = "Value")
    plot_ly(d_long, x = ~time_d, y = ~Value, color = ~Biomarker,
            type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Days"),
             yaxis = list(title = "Value (units vary)"),
             legend = list(orientation = "h"))
  })

  output$bio_table <- renderDT({
    d <- sim_out() %>%
      filter(time_d %in% c(0, 30, 90, 180, 365)) %>%
      select(time_d, CartPct, RPDAI_proxy, CRP_proxy, Ab_CII,
             TNF, IL6, IL17, IL1b, MMP) %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
    datatable(d, rownames = FALSE,
              colnames = c("Day", "Cart (%)", "RPDAI", "CRP", "Anti-CII (nM)",
                           "TNF", "IL-6", "IL-17", "IL-1β", "MMP"))
  })

  ## ── Tab 7: Sensitivity ─────────────────────────────────────
  sens_out <- eventReactive(input$run_sens, {
    p_name   <- input$sens_param
    base_val <- param(mod_base)[[p_name]]
    mults    <- seq(input$sens_range[1], input$sens_range[2], length.out = 15)
    eval_day <- input$sens_day
    out_col  <- input$sens_outcome

    lapply(mults, function(m) {
      p_list <- setNames(list(base_val * m), p_name)
      d <- run_sim(p_list, sim_days = max(eval_day + 10, 40))
      val <- d %>% filter(abs(time_d - eval_day) < 1) %>%
             slice(1) %>% pull(all_of(out_col))
      data.frame(multiplier = m, value = ifelse(length(val) == 0, NA, val))
    }) %>% bind_rows()
  })

  output$sens_plot <- renderPlotly({
    d <- sens_out()
    p_name <- input$sens_param
    plot_ly(d, x = ~multiplier, y = ~value, type = "scatter", mode = "lines+markers",
            line = list(color = "darkorange", width = 2),
            marker = list(size = 6)) %>%
      layout(
        title  = paste("Sensitivity:", p_name, "→", input$sens_outcome, "at day", input$sens_day),
        xaxis  = list(title = paste("Parameter Multiplier (×", p_name, ")")),
        yaxis  = list(title = input$sens_outcome),
        shapes = list(list(type = "line", x0 = 1, x1 = 1, y0 = 0, y1 = 1,
                           yref = "paper", line = list(dash = "dot", color = "gray")))
      )
  })
}

## ── Launch ───────────────────────────────────────────────────
shinyApp(ui, server)
