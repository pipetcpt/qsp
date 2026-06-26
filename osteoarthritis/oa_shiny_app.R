## ============================================================
## Osteoarthritis QSP — Interactive Shiny Dashboard
## ============================================================
## Tabs: Patient Profile · Drug PK · Inflammatory Markers ·
##       Cartilage & Bone · Pain & Function · Scenario Compare ·
##       Biomarkers · Sensitivity Analysis
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)
library(DT)
library(scales)

## ── Embedded model code ─────────────────────────────────────
oa_code <- '
$PROB OA QSP Shiny Model

$PARAM
k_IL1b_syn=0.10 k_IL1b_deg=0.35 k_TNFa_syn=0.08 k_TNFa_deg=0.40
k_MMP13_syn=0.15 k_MMP13_deg=0.30 k_ADAM5_syn=0.12 k_ADAM5_deg=0.28
K_IL1_mmp=2.0 Emax_IL1_mmp=3.0 K_IL1_adam5=2.5
k_ColII_syn=0.008 k_ColII_deg=0.12 K_mmp_col=2.0
k_Agg_syn=0.010 k_Agg_deg=0.15 K_adam5_agg=1.5
k_Chondro_prof=0.003 k_Chondro_apt=0.006 Emax_IL1_apt=2.0 K_IL1_apt=5.0
k_Syn_syn=0.05 k_Syn_deg=0.10 Emax_ecm_syn=1.5 K_ecm_syn=30.0
k_OC_form=0.08 k_OC_deg=0.15 Emax_IL1_OC=1.8 K_IL1_OC=3.0
k_OB_form=0.06 k_OB_deg=0.10
k_JSW_loss=0.0012 k_JSW_repair=0.0001
k_PGE2_syn=0.20 k_PGE2_deg=0.60
k_Pain_syn=0.08 k_Pain_deg=0.12 Pain_baseline=10.0
Emax_PGE2_pain=60.0 K_PGE2_pain=0.5 Emax_Syn_pain=20.0 K_Syn_pain=1.0
k_CTXII_syn=0.05 k_CTXII_deg=0.25 k_COMP_rel=0.04 k_COMP_deg=0.20
ka_nsaid=1.20 F_nsaid=0.74 Vd_nsaid=455 CL_nsaid=27.7
Kp_nsaid_jt=0.40 IC50_cox2=0.042 Emax_cox2=0.95 hill_cox2=1.2
k_IACS_abs=0.008 k_IACS_el=0.25 Vd_IACS=99
EC50_GR=0.05 Emax_GR_NF=0.75
k_HA_deg=0.041 Emax_HA_pain=0.35 K_HA_pain=0.8 Emax_HA_IL1=0.25
k_Sprif_deg=0.044 Emax_FGF_col=0.40 K_FGF18=0.002
ka_tanz=0.0167 F_tanz=0.73 Vd_tanz=7.8 CL_tanz=0.0133
Emax_tanz_pain=45.0 K_tanz_NGF=0.001
KL_grade=2 Age_yr=62 k_age_mod=0.005
ColII_0=100 Aggrecan_0=100 Chondro_0=100 JSW_0=5.5
NSAID_flag=0 IACS_flag=0 HA_flag=0 Sprif_flag=0 Tanz_flag=0

$INIT
A_NSAID_gut=0 A_NSAID_plasma=0 A_NSAID_joint=0
A_IACS_joint=0 A_IACS_plasma=0
A_HA_joint=0 A_Sprif_joint=0 A_Tanz_depot=0 A_Tanz_plasma=0
IL1b=5.0 TNFa=3.0 MMP13=8.0 ADAM5=2.0
ColII=70.0 Aggrecan=65.0 Chondro=75.0 Synovitis=1.0
OC_act=1.3 OB_act=1.1 JSW=3.8 PGE2_jt=1.5
VASPain=45.0 uCTXII=0.40 COMP_s=12.0

$ODE
double age_mod  = 1.0 + k_age_mod * (Age_yr - 45.0);
double KL_mod   = 0.5 + 0.25 * KL_grade;
double C_NSAID_plasma = A_NSAID_plasma / Vd_nsaid;
double C_NSAID_joint  = A_NSAID_joint / (Vd_nsaid * 0.01);
double COX2_inh = Emax_cox2 * pow(C_NSAID_joint, hill_cox2) /
                  (pow(IC50_cox2, hill_cox2) + pow(C_NSAID_joint, hill_cox2));
double COX2_eff = 1.0 - COX2_inh * NSAID_flag;
double C_IACS = A_IACS_joint;
double GR_inh = (Emax_GR_NF * C_IACS / (EC50_GR + C_IACS)) * IACS_flag;
double NF_kB_eff = 1.0 - GR_inh;
double C_HA  = A_HA_joint;
double HA_pain_inh = Emax_HA_pain  * C_HA / (K_HA_pain + C_HA) * HA_flag;
double HA_IL1_inh  = Emax_HA_IL1   * C_HA / (K_HA_pain + C_HA) * HA_flag;
double C_Sprif = A_Sprif_joint;
double FGF18_eff = Emax_FGF_col * C_Sprif / (K_FGF18 + C_Sprif) * Sprif_flag;
double C_Tanz = A_Tanz_plasma / Vd_tanz;
double Tanz_pain_red = Emax_tanz_pain * C_Tanz / (K_tanz_NGF + C_Tanz) * Tanz_flag;
double ECM_loss = 100.0 - 0.5*(ColII + Aggrecan);
double ECM_loss_pos = (ECM_loss < 0) ? 0.0 : ECM_loss;
dxdt_A_NSAID_gut    = -ka_nsaid * A_NSAID_gut;
dxdt_A_NSAID_plasma = F_nsaid * ka_nsaid * A_NSAID_gut * NSAID_flag
                      - (CL_nsaid/Vd_nsaid)*A_NSAID_plasma
                      - Kp_nsaid_jt*(CL_nsaid/Vd_nsaid)*A_NSAID_plasma;
dxdt_A_NSAID_joint  = Kp_nsaid_jt*(CL_nsaid/Vd_nsaid)*A_NSAID_plasma
                      - (CL_nsaid/Vd_nsaid)*A_NSAID_joint;
dxdt_A_IACS_joint   = -k_IACS_abs * A_IACS_joint;
dxdt_A_IACS_plasma  = k_IACS_abs * A_IACS_joint - k_IACS_el * A_IACS_plasma;
dxdt_A_HA_joint     = -k_HA_deg * A_HA_joint;
dxdt_A_Sprif_joint  = -k_Sprif_deg * A_Sprif_joint;
dxdt_A_Tanz_depot   = -ka_tanz * A_Tanz_depot;
dxdt_A_Tanz_plasma  = F_tanz * ka_tanz * A_Tanz_depot - (CL_tanz/Vd_tanz)*A_Tanz_plasma;
double IL1b_pos = (IL1b < 0.01) ? 0.01 : IL1b;
dxdt_IL1b = k_IL1b_syn * Synovitis * age_mod * KL_mod * NF_kB_eff * (1.0-HA_IL1_inh)
            - k_IL1b_deg * IL1b;
dxdt_TNFa = k_TNFa_syn * Synovitis * NF_kB_eff - k_TNFa_deg * TNFa;
double IL1b_TNFa_stim = (IL1b/(K_IL1_mmp+IL1b))*Emax_IL1_mmp + TNFa/(1.0+TNFa);
dxdt_MMP13 = k_MMP13_syn * IL1b_TNFa_stim * NF_kB_eff - k_MMP13_deg * MMP13;
dxdt_ADAM5 = k_ADAM5_syn * (IL1b/(K_IL1_adam5+IL1b)) * NF_kB_eff - k_ADAM5_deg * ADAM5;
double ColII_pos = (ColII < 0.5) ? 0.5 : ColII;
double col_syn = k_ColII_syn * (Chondro/Chondro_0) * (1.0 + FGF18_eff);
double col_deg = k_ColII_deg * MMP13/(K_mmp_col+MMP13) * ColII_pos * age_mod;
dxdt_ColII = col_syn - col_deg;
double Agg_pos = (Aggrecan < 0.5) ? 0.5 : Aggrecan;
double agg_syn = k_Agg_syn * (Chondro/Chondro_0) * (1.0 + 0.5*FGF18_eff);
double agg_deg = k_Agg_deg * ADAM5/(K_adam5_agg+ADAM5) * Agg_pos * age_mod;
dxdt_Aggrecan = agg_syn - agg_deg;
double Chondro_pos = (Chondro < 1.0) ? 1.0 : Chondro;
dxdt_Chondro = k_Chondro_prof*Chondro_pos
               - (k_Chondro_apt + k_Chondro_apt*Emax_IL1_apt*IL1b/(K_IL1_apt+IL1b))
                 * Chondro_pos * age_mod;
double Syn_drive = k_Syn_syn*(1.0+Emax_ecm_syn*ECM_loss_pos/(K_ecm_syn+ECM_loss_pos))*age_mod;
dxdt_Synovitis = Syn_drive - k_Syn_deg * Synovitis * NF_kB_eff;
dxdt_OC_act = k_OC_form*(1.0+Emax_IL1_OC*IL1b/(K_IL1_OC+IL1b)) - k_OC_deg * OC_act;
dxdt_OB_act = k_OB_form*(1.0+0.2*OC_act)*NF_kB_eff - k_OB_deg * OB_act;
double JSW_pos = (JSW < 0.5) ? 0.5 : JSW;
dxdt_JSW = k_JSW_repair*(1.0+FGF18_eff)
           - k_JSW_loss*(MMP13/(K_mmp_col+MMP13)+ADAM5/(K_adam5_agg+ADAM5))*JSW_pos*age_mod;
dxdt_PGE2_jt = k_PGE2_syn * IL1b/(2.0+IL1b) * COX2_eff - k_PGE2_deg * PGE2_jt;
double PGE2_pain = Emax_PGE2_pain * PGE2_jt/(K_PGE2_pain+PGE2_jt);
double Syn_pain  = Emax_Syn_pain  * Synovitis/(K_Syn_pain+Synovitis);
double Struct_pain = Pain_baseline * (1.0 - JSW/JSW_0);
double Pain_target = PGE2_pain + Syn_pain + Struct_pain - HA_pain_inh*20.0 - Tanz_pain_red;
double Pain_clamped = (Pain_target<0)?0.0:(Pain_target>100)?100.0:Pain_target;
dxdt_VASPain = k_Pain_syn * (Pain_clamped - VASPain);
double CTX2_syn = k_CTXII_syn * col_deg * 10.0;
dxdt_uCTXII = CTX2_syn - k_CTXII_deg * uCTXII;
dxdt_COMP_s = k_COMP_rel*(1.0+Synovitis+MMP13/5.0) - k_COMP_deg * COMP_s;

$TABLE
double C_nsaid_pl = A_NSAID_plasma/Vd_nsaid;
double C_nsaid_jt = A_NSAID_joint/(Vd_nsaid*0.01);
double C_IACS_jt  = A_IACS_joint;
double C_HA_jt    = A_HA_joint;
double C_Sprif_jt = A_Sprif_joint;
double C_Tanz_pl  = A_Tanz_plasma/Vd_tanz;
double COX2_pct   = Emax_cox2*pow(C_nsaid_jt,hill_cox2)/(pow(IC50_cox2,hill_cox2)+pow(C_nsaid_jt,hill_cox2))*NSAID_flag*100.0;
double KOOS_est = 0.6*(100.0-VASPain*0.60) + 0.4*(50.0+(JSW/JSW_0)*30.0-(Synovitis/3.0)*10.0);
KOOS_est=(KOOS_est<0)?0:(KOOS_est>100)?100:KOOS_est;
double CartVol = 0.5*(ColII/ColII_0)+0.5*(Aggrecan/Aggrecan_0);

$CAPTURE
C_nsaid_pl C_nsaid_jt C_IACS_jt C_HA_jt C_Sprif_jt C_Tanz_pl COX2_pct KOOS_est CartVol
IL1b TNFa MMP13 ADAM5 ColII Aggrecan Chondro Synovitis OC_act OB_act
JSW PGE2_jt VASPain uCTXII COMP_s
'

run_sim <- function(input_vals, end_days = 730) {
  mod <- mcode_cache("oa_shiny", oa_code)

  # Build dosing events
  evs <- ev()

  if (input_vals$NSAID_flag == 1) {
    ev_nsaid <- ev(amt = input_vals$nsaid_dose, ii = 12,
                   addl = end_days * 2, cmt = "A_NSAID_gut")
    evs <- ev(evs, ev_nsaid)
  }
  if (input_vals$IACS_flag == 1) {
    times_iacs <- seq(0, end_days, by = input_vals$iacs_interval)
    ev_iacs <- ev(time = times_iacs, amt = input_vals$iacs_dose,
                  cmt = "A_IACS_joint")
    evs <- ev(evs, ev_iacs)
  }
  if (input_vals$HA_flag == 1) {
    times_ha <- c(0, 7, 14, 182, 189, 196)
    times_ha <- times_ha[times_ha <= end_days]
    ev_ha <- ev(time = times_ha, amt = input_vals$ha_dose, cmt = "A_HA_joint")
    evs <- ev(evs, ev_ha)
  }
  if (input_vals$Sprif_flag == 1) {
    times_s <- seq(0, end_days, by = 84)
    ev_s <- ev(time = times_s, amt = 0.030, cmt = "A_Sprif_joint")
    evs <- ev(evs, ev_s)
  }
  if (input_vals$Tanz_flag == 1) {
    times_t <- seq(0, end_days, by = 56)
    ev_t <- ev(time = times_t, amt = 2.5, cmt = "A_Tanz_depot")
    evs <- ev(evs, ev_t)
  }

  p_list <- list(
    KL_grade  = input_vals$KL_grade,
    Age_yr    = input_vals$Age_yr,
    NSAID_flag= input_vals$NSAID_flag,
    IACS_flag = input_vals$IACS_flag,
    HA_flag   = input_vals$HA_flag,
    Sprif_flag= input_vals$Sprif_flag,
    Tanz_flag = input_vals$Tanz_flag
  )

  init_vals <- list(
    IL1b    = 2.0 + input_vals$KL_grade * 1.5,
    TNFa    = 1.0 + input_vals$KL_grade * 0.8,
    MMP13   = 2.0 + input_vals$KL_grade * 2.0,
    ColII   = max(10, 100 - input_vals$KL_grade * 15),
    Aggrecan= max(10, 100 - input_vals$KL_grade * 18),
    Chondro = max(20, 100 - input_vals$KL_grade * 10),
    Synovitis= 0.3 + input_vals$KL_grade * 0.4,
    JSW     = max(1.0, 5.5 - input_vals$KL_grade * 0.8),
    VASPain = 10 + input_vals$KL_grade * 10
  )

  mod_r <- do.call(param, c(list(mod), p_list))
  mod_r <- do.call(init,  c(list(mod_r), init_vals))

  out <- mrgsim(mod_r, ev = evs, end = end_days, delta = 1) %>%
         as_tibble()
  out
}

## ── UI ──────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "OA QSP Dashboard", titleWidth = 280),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("1. Patient Profile",    tabName = "profile",   icon = icon("user")),
      menuItem("2. Drug PK",            tabName = "pk",        icon = icon("pills")),
      menuItem("3. Inflammatory Markers",tabName = "inflam",   icon = icon("fire")),
      menuItem("4. Cartilage & Bone",   tabName = "cartilage", icon = icon("bone")),
      menuItem("5. Pain & Function",    tabName = "pain",      icon = icon("chart-line")),
      menuItem("6. Scenario Comparison",tabName = "scenario",  icon = icon("flask")),
      menuItem("7. Biomarkers",         tabName = "biomarkers",icon = icon("microscope")),
      menuItem("8. Sensitivity Analysis",tabName = "sensitivity",icon = icon("sliders-h"))
    ),
    hr(),
    h5("Patient Profile", style = "padding-left:15px; color:#90CAF9;"),
    sliderInput("KL_grade", "KL Grade (1–4):", min=1, max=4, value=2, step=1),
    sliderInput("Age_yr",   "Age (years):",    min=40, max=85, value=62),
    hr(),
    h5("Treatment Selection", style = "padding-left:15px; color:#90CAF9;"),
    checkboxInput("NSAID_flag",  "NSAID (Celecoxib 200mg BID)",   FALSE),
    conditionalPanel("input.NSAID_flag",
      sliderInput("nsaid_dose", "NSAID dose (mg):", 100, 400, 200, 100)),
    checkboxInput("IACS_flag",   "IA Corticosteroid (Triamcinolone)", FALSE),
    conditionalPanel("input.IACS_flag",
      sliderInput("iacs_dose", "IACS dose (mg):", 20, 80, 40, 10),
      sliderInput("iacs_interval", "Interval (days):", 60, 180, 91, 7)),
    checkboxInput("HA_flag",     "IA Hyaluronic Acid (16mg×3)",   FALSE),
    conditionalPanel("input.HA_flag",
      sliderInput("ha_dose", "HA dose (mg):", 8, 32, 16, 4)),
    checkboxInput("Sprif_flag",  "Sprifermin (FGF-18 30μg Q12w)", FALSE),
    checkboxInput("Tanz_flag",   "Tanezumab (anti-NGF 2.5mg Q8w)",FALSE),
    hr(),
    sliderInput("sim_years", "Simulation period (years):", 1, 3, 2, 0.5),
    actionButton("run_sim", "Run Simulation", icon = icon("play"),
                 style = "background-color:#1565C0; color:white; width:100%;")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box { border-radius: 8px; }
      .info-box { border-radius: 8px; }
    "))),

    tabItems(

      ## TAB 1: Patient Profile
      tabItem(tabName = "profile",
        fluidRow(
          valueBoxOutput("vb_kl"),
          valueBoxOutput("vb_age"),
          valueBoxOutput("vb_pain")
        ),
        fluidRow(
          box(title = "Patient Disease State Summary", width = 12, status = "primary", solidHeader = TRUE,
            DTOutput("profile_table"))
        ),
        fluidRow(
          box(title = "Disease Progression Overview", width = 12, status = "info",
            plotlyOutput("profile_overview", height = "400px"))
        )
      ),

      ## TAB 2: Drug PK
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "NSAID Plasma & Joint Concentrations", width = 6, status = "primary",
            plotlyOutput("pk_nsaid", height = "300px")),
          box(title = "COX-2 Inhibition (%)", width = 6, status = "warning",
            plotlyOutput("pk_cox2", height = "300px"))
        ),
        fluidRow(
          box(title = "IA Corticosteroid Joint Kinetics", width = 6, status = "success",
            plotlyOutput("pk_iacs", height = "300px")),
          box(title = "IA HA / Sprifermin / Tanezumab", width = 6, status = "info",
            plotlyOutput("pk_other", height = "300px"))
        )
      ),

      ## TAB 3: Inflammatory Markers
      tabItem(tabName = "inflam",
        fluidRow(
          box(title = "IL-1β (pg/mL)", width = 4, status = "danger",
            plotlyOutput("inflam_il1b", height = "250px")),
          box(title = "TNF-α (pg/mL)", width = 4, status = "danger",
            plotlyOutput("inflam_tnf", height = "250px")),
          box(title = "Synovitis Score (0–3)", width = 4, status = "danger",
            plotlyOutput("inflam_syn", height = "250px"))
        ),
        fluidRow(
          box(title = "MMP-13 (ng/mL)", width = 4, status = "warning",
            plotlyOutput("inflam_mmp13", height = "250px")),
          box(title = "ADAMTS-5 Activity", width = 4, status = "warning",
            plotlyOutput("inflam_adam5", height = "250px")),
          box(title = "PGE2 Joint (ng/mL)", width = 4, status = "warning",
            plotlyOutput("inflam_pge2", height = "250px"))
        )
      ),

      ## TAB 4: Cartilage & Bone
      tabItem(tabName = "cartilage",
        fluidRow(
          box(title = "Collagen Type II (%)", width = 6, status = "primary",
            plotlyOutput("cart_colii", height = "280px")),
          box(title = "Aggrecan (%)", width = 6, status = "primary",
            plotlyOutput("cart_agg", height = "280px"))
        ),
        fluidRow(
          box(title = "Chondrocyte Pool (%)", width = 4, status = "success",
            plotlyOutput("cart_chondro", height = "260px")),
          box(title = "Joint Space Width (mm)", width = 4, status = "warning",
            plotlyOutput("cart_jsw", height = "260px")),
          box(title = "Osteoclast / Osteoblast Activity", width = 4, status = "info",
            plotlyOutput("cart_bone", height = "260px"))
        )
      ),

      ## TAB 5: Pain & Function
      tabItem(tabName = "pain",
        fluidRow(
          box(title = "VAS Pain Score (0–100)", width = 6, status = "danger",
            plotlyOutput("pain_vas", height = "300px")),
          box(title = "Estimated KOOS (0–100, higher=better)", width = 6, status = "primary",
            plotlyOutput("pain_koos", height = "300px"))
        ),
        fluidRow(
          box(title = "Pain Component Breakdown (at final timepoint)", width = 12, status = "warning",
            plotlyOutput("pain_components", height = "280px"))
        )
      ),

      ## TAB 6: Scenario Comparison
      tabItem(tabName = "scenario",
        fluidRow(
          box(title = "Multi-Scenario Settings", width = 12, status = "primary",
            p("Comparing 6 pre-defined scenarios: Natural history · Celecoxib · IA-CS · HA · Sprifermin · Tanezumab"),
            actionButton("run_all", "Run All Scenarios", icon = icon("play-circle"),
                         style = "background-color:#2E7D32; color:white;"))
        ),
        fluidRow(
          box(title = "VAS Pain — All Scenarios", width = 6, status = "danger",
            plotlyOutput("scen_pain", height = "300px")),
          box(title = "JSW — All Scenarios", width = 6, status = "warning",
            plotlyOutput("scen_jsw", height = "300px"))
        ),
        fluidRow(
          box(title = "KOOS — All Scenarios", width = 6, status = "primary",
            plotlyOutput("scen_koos", height = "300px")),
          box(title = "1-Year Outcome Table", width = 6, status = "info",
            DTOutput("scen_table"))
        )
      ),

      ## TAB 7: Biomarkers
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(title = "uCTX-II (nmol/mmolCr)", width = 6, status = "primary",
            plotlyOutput("bm_ctxii", height = "280px")),
          box(title = "Serum COMP (μg/mL)", width = 6, status = "info",
            plotlyOutput("bm_comp", height = "280px"))
        ),
        fluidRow(
          box(title = "Biomarker Clinical Reference Ranges", width = 12, status = "success",
            tableOutput("bm_reference"))
        )
      ),

      ## TAB 8: Sensitivity Analysis
      tabItem(tabName = "sensitivity",
        fluidRow(
          box(title = "Sensitivity Analysis Controls", width = 12, status = "primary",
            selectInput("sens_outcome", "Outcome Variable:",
                        choices = c("JSW (mm)" = "JSW",
                                    "VAS Pain"  = "VASPain",
                                    "KOOS"      = "KOOS_est",
                                    "ColII (%)" = "ColII",
                                    "uCTX-II"   = "uCTXII"),
                        selected = "JSW"),
            selectInput("sens_param", "Vary Parameter:",
                        choices = c("MMP-13 synthesis" = "k_MMP13_syn",
                                    "Collagen II degradation" = "k_ColII_deg",
                                    "ADAMTS-5 synthesis" = "k_ADAM5_syn",
                                    "JSW loss rate" = "k_JSW_loss",
                                    "Pain synthesis" = "k_Pain_syn",
                                    "IL-1b synthesis" = "k_IL1b_syn"),
                        selected = "k_MMP13_syn"),
            sliderInput("sens_range", "Parameter range (fold):",
                        min = 0.25, max = 3.0, value = c(0.7, 1.3), step = 0.05),
            actionButton("run_sens", "Run Sensitivity Analysis",
                         style = "background-color:#5C6BC0; color:white;"))
        ),
        fluidRow(
          box(title = "Sensitivity — Time Course", width = 8, status = "primary",
            plotlyOutput("sens_timecourse", height = "350px")),
          box(title = "Sensitivity — Tornado at 1 Year", width = 4, status = "warning",
            plotlyOutput("sens_tornado", height = "350px"))
        )
      )
    )
  )
)

## ── SERVER ───────────────────────────────────────────────────
server <- function(input, output, session) {

  # Reactive simulation result (current patient settings)
  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Running simulation...", {
      iv <- list(
        KL_grade    = input$KL_grade,
        Age_yr      = input$Age_yr,
        NSAID_flag  = as.integer(input$NSAID_flag),
        IACS_flag   = as.integer(input$IACS_flag),
        HA_flag     = as.integer(input$HA_flag),
        Sprif_flag  = as.integer(input$Sprif_flag),
        Tanz_flag   = as.integer(input$Tanz_flag),
        nsaid_dose  = if(input$NSAID_flag) input$nsaid_dose else 200,
        iacs_dose   = if(input$IACS_flag)  input$iacs_dose  else 40,
        iacs_interval=if(input$IACS_flag)  input$iacs_interval else 91,
        ha_dose     = if(input$HA_flag)    input$ha_dose    else 16
      )
      run_sim(iv, end_days = round(input$sim_years * 365))
    })
  }, ignoreNULL = FALSE)

  # Value boxes
  output$vb_kl   <- renderValueBox({
    valueBox(paste0("KL ", input$KL_grade), "KL Grade",
             icon = icon("bone"), color = if(input$KL_grade <= 2) "green" else "red")
  })
  output$vb_age  <- renderValueBox({
    valueBox(paste0(input$Age_yr, " yr"), "Patient Age",
             icon = icon("user"), color = "blue")
  })
  output$vb_pain <- renderValueBox({
    d <- sim_data()
    pain_now <- round(tail(d$VASPain, 1), 1)
    valueBox(pain_now, "Final VAS Pain",
             icon = icon("chart-line"),
             color = if(pain_now < 40) "green" else if(pain_now < 65) "yellow" else "red")
  })

  # Profile table
  output$profile_table <- renderDT({
    d <- sim_data()
    last <- tail(d, 1)
    tibble(
      Variable  = c("VAS Pain (0-100)", "KOOS (estimated)", "JSW (mm)",
                    "Collagen II (%)", "Aggrecan (%)", "Chondrocyte (%)",
                    "Synovitis Score", "IL-1β (pg/mL)", "MMP-13 (ng/mL)",
                    "uCTX-II (nmol/mmolCr)", "Serum COMP (μg/mL)"),
      `Baseline` = c(round(d$VASPain[1],1), round(d$KOOS_est[1],1), round(d$JSW[1],2),
                     round(d$ColII[1],1), round(d$Aggrecan[1],1), round(d$Chondro[1],1),
                     round(d$Synovitis[1],2), round(d$IL1b[1],2), round(d$MMP13[1],2),
                     round(d$uCTXII[1],3), round(d$COMP_s[1],2)),
      `Final`    = c(round(last$VASPain,1), round(last$KOOS_est,1), round(last$JSW,2),
                     round(last$ColII,1), round(last$Aggrecan,1), round(last$Chondro,1),
                     round(last$Synovitis,2), round(last$IL1b,2), round(last$MMP13,2),
                     round(last$uCTXII,3), round(last$COMP_s,2)),
      `Change`   = c(
        paste0(ifelse(last$VASPain-d$VASPain[1]>0,"+",""), round(last$VASPain-d$VASPain[1],1)),
        paste0(ifelse(last$KOOS_est-d$KOOS_est[1]>0,"+",""), round(last$KOOS_est-d$KOOS_est[1],1)),
        paste0(ifelse(last$JSW-d$JSW[1]>0,"+",""), round(last$JSW-d$JSW[1],2)),
        paste0(ifelse(last$ColII-d$ColII[1]>0,"+",""), round(last$ColII-d$ColII[1],1)),
        paste0(ifelse(last$Aggrecan-d$Aggrecan[1]>0,"+",""), round(last$Aggrecan-d$Aggrecan[1],1)),
        paste0(ifelse(last$Chondro-d$Chondro[1]>0,"+",""), round(last$Chondro-d$Chondro[1],1)),
        paste0(ifelse(last$Synovitis-d$Synovitis[1]>0,"+",""), round(last$Synovitis-d$Synovitis[1],2)),
        paste0(ifelse(last$IL1b-d$IL1b[1]>0,"+",""), round(last$IL1b-d$IL1b[1],2)),
        paste0(ifelse(last$MMP13-d$MMP13[1]>0,"+",""), round(last$MMP13-d$MMP13[1],2)),
        paste0(ifelse(last$uCTXII-d$uCTXII[1]>0,"+",""), round(last$uCTXII-d$uCTXII[1],3)),
        paste0(ifelse(last$COMP_s-d$COMP_s[1]>0,"+",""), round(last$COMP_s-d$COMP_s[1],2)))
    ) %>% datatable(options = list(pageLength = 15, dom = "t"), rownames = FALSE)
  })

  # Profile overview
  output$profile_overview <- renderPlotly({
    d <- sim_data()
    p <- d %>% select(time, VASPain, JSW, ColII, Aggrecan, KOOS_est) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, color = name)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", nrow = 1) +
      labs(x = "Day", y = "", color = "Variable") +
      theme_bw(base_size = 10)
    ggplotly(p)
  })

  # PK plots
  output$pk_nsaid <- renderPlotly({
    d <- sim_data()
    p <- d %>% select(time, C_nsaid_pl, C_nsaid_jt) %>%
      pivot_longer(-time, names_to = "Compartment", values_to = "Conc") %>%
      mutate(Compartment = recode(Compartment, C_nsaid_pl = "Plasma (μg/mL)",
                                               C_nsaid_jt = "Joint (μg/mL)")) %>%
      ggplot(aes(time, Conc, color = Compartment)) +
      geom_line() + labs(x = "Day", y = "NSAID (μg/mL)") + theme_bw()
    ggplotly(p)
  })

  output$pk_cox2 <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time, COX2_pct)) + geom_line(color = "#E53935") +
      geom_hline(yintercept = 80, linetype = "dashed", color = "gray50") +
      labs(x = "Day", y = "COX-2 Inhibition (%)") + theme_bw() +
      annotate("text", x = max(d$time)*0.8, y = 82, label = "80% threshold", size = 3)
    ggplotly(p)
  })

  output$pk_iacs <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time, C_IACS_jt)) + geom_line(color = "#FF8A65") +
      labs(x = "Day", y = "IACS Joint Conc (μg/mL)") + theme_bw()
    ggplotly(p)
  })

  output$pk_other <- renderPlotly({
    d <- sim_data()
    p <- d %>% select(time, C_HA_jt, C_Sprif_jt, C_Tanz_pl) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, color = name)) + geom_line() +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Day", y = "Concentration") + theme_bw()
    ggplotly(p)
  })

  # Inflammatory plots (helper)
  mk_inflam_plot <- function(d, var, ylab, col) {
    ggplot(d, aes_string("time", var)) + geom_line(color = col, linewidth = 1) +
      labs(x = "Day", y = ylab) + theme_bw()
  }
  output$inflam_il1b  <- renderPlotly({ ggplotly(mk_inflam_plot(sim_data(), "IL1b",  "IL-1β (pg/mL)", "#E91E63")) })
  output$inflam_tnf   <- renderPlotly({ ggplotly(mk_inflam_plot(sim_data(), "TNFa",  "TNF-α (pg/mL)", "#F44336")) })
  output$inflam_syn   <- renderPlotly({ ggplotly(mk_inflam_plot(sim_data(), "Synovitis", "Synovitis (0–3)", "#FF5722")) })
  output$inflam_mmp13 <- renderPlotly({ ggplotly(mk_inflam_plot(sim_data(), "MMP13", "MMP-13 (ng/mL)", "#FF9800")) })
  output$inflam_adam5 <- renderPlotly({ ggplotly(mk_inflam_plot(sim_data(), "ADAM5", "ADAMTS-5 (norm.)", "#FFC107")) })
  output$inflam_pge2  <- renderPlotly({ ggplotly(mk_inflam_plot(sim_data(), "PGE2_jt","PGE2 (ng/mL)", "#FF8F00")) })

  # Cartilage plots
  output$cart_colii   <- renderPlotly({ ggplotly(mk_inflam_plot(sim_data(), "ColII",   "Collagen II (%)", "#1565C0")) })
  output$cart_agg     <- renderPlotly({ ggplotly(mk_inflam_plot(sim_data(), "Aggrecan","Aggrecan (%)", "#283593")) })
  output$cart_chondro <- renderPlotly({ ggplotly(mk_inflam_plot(sim_data(), "Chondro", "Chondrocyte (%)", "#2E7D32")) })
  output$cart_jsw     <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time, JSW)) + geom_line(color = "#F57F17", linewidth = 1) +
      geom_hline(yintercept = 2.0, linetype = "dashed", color = "red", linewidth = 0.5) +
      annotate("text", x = max(d$time)*0.7, y = 2.1, label = "TKR risk ~2mm", size = 3, color = "red") +
      labs(x = "Day", y = "JSW (mm)") + theme_bw()
    ggplotly(p)
  })
  output$cart_bone <- renderPlotly({
    d <- sim_data()
    p <- d %>% select(time, OC_act, OB_act) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, color = name)) + geom_line() +
      labs(x = "Day", y = "Normalized Activity", color = "") + theme_bw()
    ggplotly(p)
  })

  # Pain plots
  output$pain_vas <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time, VASPain)) + geom_line(color = "#C62828", linewidth = 1) +
      geom_hline(yintercept = 40, linetype = "dashed", color = "#4CAF50") +
      annotate("text", x = max(d$time)*0.7, y = 42, label = "Clinically meaningful <40", size=3, color="#4CAF50") +
      labs(x = "Day", y = "VAS Pain (0–100)") + scale_y_continuous(limits = c(0, 100)) + theme_bw()
    ggplotly(p)
  })
  output$pain_koos <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time, KOOS_est)) + geom_line(color = "#1565C0", linewidth = 1) +
      geom_hline(yintercept = 70, linetype = "dashed", color = "#4CAF50") +
      labs(x = "Day", y = "KOOS (0–100)") + scale_y_continuous(limits = c(0, 100)) + theme_bw()
    ggplotly(p)
  })
  output$pain_components <- renderPlotly({
    d <- sim_data()
    last <- tail(d, 1)
    comps <- tibble(
      Component = c("PGE2-mediated", "Synovitis", "Structural", "HA relief", "Tanezumab"),
      Value = c(
        60 * last$PGE2_jt / (0.5 + last$PGE2_jt),
        20 * last$Synovitis / (1 + last$Synovitis),
        10 * (1 - last$JSW / 5.5),
        -0.35 * last$C_HA_jt / (0.8 + last$C_HA_jt) * 20 * as.integer(input$HA_flag),
        -45 * last$C_Tanz_pl / (0.001 + last$C_Tanz_pl) * as.integer(input$Tanz_flag))
    )
    p <- ggplot(comps, aes(Component, Value, fill = Value > 0)) +
      geom_col() + coord_flip() +
      scale_fill_manual(values = c("TRUE" = "#E53935", "FALSE" = "#4CAF50"),
                        labels = c("Pain relief", "Pain contribution")) +
      labs(x = "", y = "VAS Pain Contribution", fill = "") + theme_bw()
    ggplotly(p)
  })

  # Scenario comparison
  all_scen_data <- eventReactive(input$run_all, {
    withProgress(message = "Running all scenarios...", {
      scens <- list(
        list(label="1. Natural History",      NSAID_flag=0,IACS_flag=0,HA_flag=0,Sprif_flag=0,Tanz_flag=0),
        list(label="2. Celecoxib",            NSAID_flag=1,IACS_flag=0,HA_flag=0,Sprif_flag=0,Tanz_flag=0),
        list(label="3. IA Corticosteroid",    NSAID_flag=0,IACS_flag=1,HA_flag=0,Sprif_flag=0,Tanz_flag=0),
        list(label="4. IA Hyaluronic Acid",   NSAID_flag=0,IACS_flag=0,HA_flag=1,Sprif_flag=0,Tanz_flag=0),
        list(label="5. Sprifermin",           NSAID_flag=0,IACS_flag=0,HA_flag=0,Sprif_flag=1,Tanz_flag=0),
        list(label="6. Tanezumab",            NSAID_flag=0,IACS_flag=0,HA_flag=0,Sprif_flag=0,Tanz_flag=1)
      )
      purrr::map_dfr(scens, function(s) {
        iv <- c(s, list(KL_grade=input$KL_grade, Age_yr=input$Age_yr,
                        nsaid_dose=200, iacs_dose=40, iacs_interval=91, ha_dose=16))
        run_sim(iv, end_days=round(input$sim_years*365)) %>% mutate(Scenario=s$label)
      })
    })
  })

  colors6 <- c("1. Natural History"="#B71C1C","2. Celecoxib"="#1565C0",
               "3. IA Corticosteroid"="#2E7D32","4. IA Hyaluronic Acid"="#6A1B9A",
               "5. Sprifermin"="#E65100","6. Tanezumab"="#00838F")

  output$scen_pain <- renderPlotly({
    req(all_scen_data())
    p <- ggplot(all_scen_data(), aes(time, VASPain, color=Scenario)) + geom_line() +
      scale_color_manual(values=colors6) + labs(x="Day",y="VAS Pain") + theme_bw()
    ggplotly(p)
  })
  output$scen_jsw  <- renderPlotly({
    req(all_scen_data())
    p <- ggplot(all_scen_data(), aes(time, JSW, color=Scenario)) + geom_line() +
      scale_color_manual(values=colors6) + labs(x="Day",y="JSW (mm)") + theme_bw()
    ggplotly(p)
  })
  output$scen_koos <- renderPlotly({
    req(all_scen_data())
    p <- ggplot(all_scen_data(), aes(time, KOOS_est, color=Scenario)) + geom_line() +
      scale_color_manual(values=colors6) + labs(x="Day",y="KOOS") + theme_bw()
    ggplotly(p)
  })
  output$scen_table <- renderDT({
    req(all_scen_data())
    all_scen_data() %>%
      filter(abs(time - 365) < 2) %>%
      group_by(Scenario) %>%
      summarise(Pain=round(mean(VASPain),1), JSW=round(mean(JSW),2),
                KOOS=round(mean(KOOS_est),1), ColII=round(mean(ColII),1),
                uCTXII=round(mean(uCTXII),3), .groups="drop") %>%
      datatable(options=list(pageLength=10, dom="t"), rownames=FALSE)
  })

  # Biomarkers
  output$bm_ctxii <- renderPlotly({ ggplotly(mk_inflam_plot(sim_data(), "uCTXII", "uCTX-II (nmol/mmolCr)", "#AFB42B")) })
  output$bm_comp  <- renderPlotly({ ggplotly(mk_inflam_plot(sim_data(), "COMP_s",  "COMP (μg/mL)", "#00796B")) })
  output$bm_reference <- renderTable({
    tibble(
      Biomarker    = c("uCTX-II","Serum COMP","Serum HA","Serum MMP-3","hs-CRP"),
      `Normal Range`=c("<0.20 nmol/mmolCr","3–8 μg/mL","20–50 ng/mL","5–12 ng/mL","<1 mg/L"),
      `OA Elevated` =c(">0.40","10–30","100–300","30–100",">3"),
      `Reference`   =c("Lohmander 2003 Ann Rheum Dis",
                        "Spector 1995 Arthritis Rheum",
                        "Lindqvist 1997 Br J Rheum",
                        "Nakamura 2000 Arthritis Rheum",
                        "Stannus 2010 Arthritis Rheum")
    )
  })

  # Sensitivity
  sens_results <- eventReactive(input$run_sens, {
    withProgress(message = "Running sensitivity...", {
      mults <- seq(input$sens_range[1], input$sens_range[2], length.out = 5)
      mod_base <- mcode_cache("oa_shiny", oa_code)
      base_val <- param(mod_base)[[input$sens_param]]
      purrr::map_dfr(mults, function(m) {
        iv <- list(KL_grade=input$KL_grade, Age_yr=input$Age_yr,
                   NSAID_flag=0, IACS_flag=0, HA_flag=0, Sprif_flag=0, Tanz_flag=0,
                   nsaid_dose=200, iacs_dose=40, iacs_interval=91, ha_dose=16)
        d <- run_sim(iv, end_days=365)
        p_new <- setNames(list(base_val * m), input$sens_param)
        mod_s <- do.call(param, c(list(mod_base), p_new))
        out <- mrgsim(mod_s, ev=ev(), end=365, delta=7) %>% as_tibble()
        out$mult <- round(m, 2)
        out
      })
    })
  })

  output$sens_timecourse <- renderPlotly({
    req(sens_results())
    d <- sens_results()
    yvar <- input$sens_outcome
    p <- ggplot(d, aes_string("time", yvar, color="factor(mult)", group="factor(mult)")) +
      geom_line(linewidth = 0.9) +
      scale_color_viridis_d(name = "Parameter\nMultiplier") +
      labs(x = "Day", y = yvar) + theme_bw()
    ggplotly(p)
  })

  output$sens_tornado <- renderPlotly({
    req(sens_results())
    d <- sens_results()
    yvar <- input$sens_outcome
    df_end <- d %>% filter(time == max(time)) %>%
      group_by(mult) %>% summarise(val = mean(.data[[yvar]]), .groups="drop") %>%
      mutate(param = input$sens_param)
    p <- ggplot(df_end, aes(val, factor(mult), fill = val)) +
      geom_col() + scale_fill_gradient2(low="#4CAF50", high="#F44336", mid="#FFF176", midpoint=median(df_end$val)) +
      labs(x = yvar, y = "Parameter Multiplier", fill = "") + theme_bw() + coord_flip()
    ggplotly(p)
  })
}

## ── LAUNCH ───────────────────────────────────────────────────
shinyApp(ui, server)
