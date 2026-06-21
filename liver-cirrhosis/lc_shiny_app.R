##############################################################################
# Liver Cirrhosis QSP — Interactive Shiny Dashboard
# ============================================================
# Tabs:
#   1. Patient Profile & Disease Staging
#   2. PK Simulation (Propranolol, Spironolactone, Terlipressin, Rifaximin)
#   3. Portal Hypertension & Hemodynamics
#   4. Hepatic Synthetic Function & Scoring
#   5. Complications (Ascites, HE, HRS)
#   6. Scenario Comparison & Clinical Outcomes
#   7. Antifibrotic Therapy & Fibrosis Kinetics
#   (Bonus) Biomarker Correlation Panel
##############################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(shinycssloaders)

##############################################################################
# Embedded model code (simplified from lc_mrgsolve_model.R for Shiny)
##############################################################################

lc_model_code <- '
$PARAM
K_FIBRO=0.004 K_FIBRO_RES=0.001 K_HSC_ACT=0.015 K_HSC_DEACT=0.008
TGF_BASE=1.0 HILL_TGF=2.0 EC50_TGF=0.5
HVPG_BASE=6.0 K_HVPG_F=12.0 ET1_BASE=1.0 eNOS_BASE=1.0
K_ET1=0.3 K_eNOS=0.2
ALB_NORM=4.2 K_ALB_DEC=0.35 INR_NORM=1.0 K_INR_INC=0.8
BILI_NORM=0.8 K_BILI_INC=3.5
K_ASCITES=0.05 HVPG_THRESH=10.0 K_LYMPH=0.02 ALDO_BASE=1.0
K_ALDO_NA=0.04 OncP_NORM=1.0 K_ONCP_ALB=0.3
GFR_NORM=90.0 K_GFR_HVPG=0.8 K_SNS_GFR=0.5
CREAT_NORM=0.9 K_CREAT_GFR=80.0
NH3_NORM=30.0 K_NH3_LF=60.0 K_NH3_BYPASS=20.0
K_NH3_GUT=1.0 K_NH3_ELIM=0.1 HE_THRESH=50.0 HE_K=0.02
PROP_F=0.25 PROP_KA=0.5 PROP_CL=30.0 PROP_V1=150.0
PROP_Q=20.0 PROP_V2=300.0 PROP_EC50=40.0 PROP_EMAX=0.30
HR_BASE=85.0 PROP_HVPG_EMAX=0.25 PROP_HVPG_EC50=35.0
SPIRO_F=0.90 SPIRO_KA=0.3 SPIRO_CL=3.5 SPIRO_V=70.0
SPIRO_EC50=150.0 SPIRO_EMAX=0.85 SPIRO_NA_EMAX=0.7
TERL_CL=15.0 TERL_V=30.0 TERL_KE0=0.4
TERL_EC50=10.0 TERL_EMAX=0.50 TERL_HRS_EMAX=0.65
RIFAX_KOUT=0.15 RIFAX_EC50_NH3=200.0 RIFAX_EMAX_NH3=0.55
AF_FIBRO_EMAX=0.6 AF_HSC_EMAX=0.7 AF_EC50=1.0
MMP_ACT_NORM=1.0 TIMP_NORM=2.0 AF_DRUG_CONC=0.0

$CMT PROP_GUT PROP_C1 PROP_C2 SPIRO_GUT SPIRO_C1
     TERL_C TERL_CE RIFAX_GUT
     FIBRO HSC_ACT HVPG ALB BILIRUBIN INR_val ASCITES
     GFR_est NH3_blood HE_GRADE CREAT ALDO_ACT

$MAIN
double PROP_CONC = PROP_C1 / PROP_V1;
double SPIRO_CONC = SPIRO_C1 / SPIRO_V;
double TERL_CONC_EFF = TERL_CE;
double RIFAX_CONC = RIFAX_GUT / 50.0;
double PROP_HR_EFF = PROP_EMAX * PROP_CONC / (PROP_EC50 + PROP_CONC);
double HR_obs = HR_BASE * (1.0 - PROP_HR_EFF);
double PROP_HVPG_EFF = PROP_HVPG_EMAX * PROP_CONC / (PROP_HVPG_EC50 + PROP_CONC);
double SPIRO_ALDO_EFF = SPIRO_EMAX * SPIRO_CONC / (SPIRO_EC50 + SPIRO_CONC);
double SPIRO_NA_EFF = SPIRO_NA_EMAX * SPIRO_CONC / (SPIRO_EC50 + SPIRO_CONC);
double TERL_VASOC_EFF = TERL_EMAX * TERL_CONC_EFF / (TERL_EC50 + TERL_CONC_EFF);
double TERL_HRS_EFF = TERL_HRS_EMAX * TERL_CONC_EFF / (TERL_EC50 + TERL_CONC_EFF);
double RIFAX_NH3_EFF = RIFAX_EMAX_NH3 * RIFAX_CONC / (RIFAX_EC50_NH3 + RIFAX_CONC);
double HVPG_fibrosis = HVPG_BASE + K_HVPG_F * FIBRO;
double HVPG_ET1_eNOS = K_ET1 * ET1_BASE - K_eNOS * eNOS_BASE;
double HVPG_drug_mod = (1.0 - PROP_HVPG_EFF) * (1.0 - TERL_VASOC_EFF);
double HVPG_calc = (HVPG_fibrosis + HVPG_ET1_eNOS) * HVPG_drug_mod;
if(HVPG_calc < 1.0) HVPG_calc = 1.0;
double ALB_calc = ALB_NORM - K_ALB_DEC * FIBRO;
if(ALB_calc < 1.0) ALB_calc = 1.0;
double BILI_calc = BILI_NORM + K_BILI_INC * FIBRO;
double INR_calc = INR_NORM + K_INR_INC * FIBRO;
double CP_ALB = (ALB_calc > 3.5) ? 1.0 : ((ALB_calc > 2.8) ? 2.0 : 3.0);
double CP_BILI = (BILI_calc < 2.0) ? 1.0 : ((BILI_calc < 3.0) ? 2.0 : 3.0);
double CP_INR = (INR_calc < 1.7) ? 1.0 : ((INR_calc < 2.3) ? 2.0 : 3.0);
double CP_ASC = (ASCITES < 1.0) ? 1.0 : ((ASCITES < 5.0) ? 2.0 : 3.0);
double CP_HE = (HE_GRADE < 0.5) ? 1.0 : ((HE_GRADE < 2.5) ? 2.0 : 3.0);
double CHILD_PUGH = CP_ALB + CP_BILI + CP_INR + CP_ASC + CP_HE;
double CREAT_MELD = CREAT; if(CREAT_MELD<1.0) CREAT_MELD=1.0; if(CREAT_MELD>4.0) CREAT_MELD=4.0;
double BILI_MELD = BILIRUBIN; if(BILI_MELD<1.0) BILI_MELD=1.0;
double INR_MELD = INR_val; if(INR_MELD<1.0) INR_MELD=1.0;
double MELD_score = 3.78*log(BILI_MELD) + 11.2*log(INR_MELD) + 9.57*log(CREAT_MELD) + 6.43;
double Na_approx = 140.0 - 5.0*(ASCITES/5.0) - 3.0*(1.0 - GFR_est/90.0);
if(Na_approx<120.0) Na_approx=120.0; if(Na_approx>140.0) Na_approx=140.0;
double MELD_Na_score = MELD_score + 1.32*(140.0-Na_approx) - 0.24*MELD_score*(140.0-Na_approx)/100.0;
double VAR_BLEED_RISK = (HVPG_calc > 12.0) ? 0.02 + 0.025*(HVPG_calc-12.0) : 0.0;
double MORT_1YR = (MELD_score<9)?0.02:(MELD_score<19)?0.06:(MELD_score<29)?0.20:(MELD_score<39)?0.52:0.70;

$ODE
double prop_abs_rate = PROP_KA * PROP_GUT;
double prop_dist_forward = (PROP_Q/PROP_V1)*PROP_C1;
double prop_dist_back = (PROP_Q/PROP_V2)*PROP_C2;
double prop_elim = (PROP_CL/PROP_V1)*PROP_C1;
dxdt_PROP_GUT = -prop_abs_rate;
dxdt_PROP_C1 = prop_abs_rate - prop_elim - prop_dist_forward + prop_dist_back;
dxdt_PROP_C2 = prop_dist_forward - prop_dist_back;
double spiro_abs_rate = SPIRO_KA * SPIRO_GUT;
double spiro_elim = (SPIRO_CL/SPIRO_V)*SPIRO_C1;
dxdt_SPIRO_GUT = -spiro_abs_rate;
dxdt_SPIRO_C1 = spiro_abs_rate - spiro_elim;
dxdt_TERL_C = -(TERL_CL/TERL_V)*TERL_C;
dxdt_TERL_CE = TERL_KE0*(TERL_C/TERL_V - TERL_CE);
dxdt_RIFAX_GUT = -RIFAX_KOUT * RIFAX_GUT;
double TGF_effective = TGF_BASE * (1.0 + 2.0*FIBRO);
double SPIRO_ALDO_EFF2 = SPIRO_EMAX*SPIRO_CONC/(SPIRO_EC50+SPIRO_CONC);
double HSC_activation_in = K_HSC_ACT * TGF_effective*TGF_effective /
    (EC50_TGF*EC50_TGF + TGF_effective*TGF_effective) * (1.0-HSC_ACT);
double HSC_activation_net = HSC_activation_in*(1.0-AF_HSC_EMAX*AF_DRUG_CONC/(AF_EC50+AF_DRUG_CONC));
dxdt_HSC_ACT = HSC_activation_net - K_HSC_DEACT*HSC_ACT;
double fibro_production = K_FIBRO * HSC_ACT * (1.0-FIBRO);
double fibro_resolution = K_FIBRO_RES * MMP_ACT_NORM/(TIMP_NORM+MMP_ACT_NORM) * FIBRO;
double fibro_af_eff = AF_FIBRO_EMAX*AF_DRUG_CONC/(AF_EC50+AF_DRUG_CONC);
dxdt_FIBRO = fibro_production*(1.0-fibro_af_eff) - fibro_resolution;
if(FIBRO>0.999 && dxdt_FIBRO>0) dxdt_FIBRO=0;
if(FIBRO<0.001 && dxdt_FIBRO<0) dxdt_FIBRO=0;
dxdt_HVPG = 0.05*(HVPG_calc - HVPG);
dxdt_ALB = 0.02*(ALB_calc - ALB);
dxdt_BILIRUBIN = 0.02*(BILI_calc - BILIRUBIN);
dxdt_INR_val = 0.02*(INR_calc - INR_val);
double PROP_HVPG_EFF2 = PROP_HVPG_EMAX*PROP_CONC/(PROP_HVPG_EC50+PROP_CONC);
double TERL_VASOC_EFF2 = TERL_EMAX*TERL_CONC_EFF/(TERL_EC50+TERL_CONC_EFF);
double portal_drive = (HVPG > HVPG_THRESH) ? (HVPG-HVPG_THRESH) : 0.0;
double onco_factor = 1.0 + K_ONCP_ALB*(ALB_NORM-ALB)/ALB_NORM;
double aldo_current = ALDO_ACT*(1.0-SPIRO_ALDO_EFF2);
double ascites_form = K_ASCITES*portal_drive*onco_factor*(1.0+K_ALDO_NA*aldo_current);
double ascites_absorb = K_LYMPH*ASCITES;
double ascites_diuretic = SPIRO_NA_EFF*0.1*ASCITES;
dxdt_ASCITES = ascites_form - ascites_absorb - ascites_diuretic;
if(ASCITES<0.0 && dxdt_ASCITES<0) dxdt_ASCITES=0;
double aldo_drive = 1.0 + 0.5*FIBRO + 0.03*(HVPG-6.0);
if(aldo_drive<0.0) aldo_drive=0.0;
dxdt_ALDO_ACT = 0.01*(aldo_drive*(1.0-SPIRO_ALDO_EFF2) - ALDO_ACT);
double TERL_HRS_EFF2 = TERL_HRS_EMAX*TERL_CONC_EFF/(TERL_EC50+TERL_CONC_EFF);
double gfr_hvpg_drive = (HVPG>12.0) ? K_GFR_HVPG*(HVPG-12.0) : 0.0;
double gfr_sns_drive = K_SNS_GFR*(ALDO_ACT-1.0);
if(gfr_sns_drive<0.0) gfr_sns_drive=0.0;
double gfr_target = GFR_NORM - gfr_hvpg_drive - gfr_sns_drive;
gfr_target += TERL_HRS_EFF2*(GFR_NORM - gfr_target);
if(gfr_target<5.0) gfr_target=5.0;
dxdt_GFR_est = 0.02*(gfr_target - GFR_est);
double creat_target = CREAT_NORM*GFR_NORM/GFR_est;
if(creat_target<0.5) creat_target=0.5; if(creat_target>10.0) creat_target=10.0;
dxdt_CREAT = 0.01*(creat_target - CREAT);
double RIFAX_NH3_EFF2 = RIFAX_EMAX_NH3*RIFAX_CONC/(RIFAX_EC50_NH3+RIFAX_CONC);
double nh3_liver_fail = K_NH3_LF*FIBRO;
double nh3_shunt = K_NH3_BYPASS*(HVPG>10.0?(HVPG-10.0)/10.0:0.0);
double nh3_gut = K_NH3_GUT*NH3_NORM*(1.0-RIFAX_NH3_EFF2);
double nh3_target = NH3_NORM + nh3_liver_fail + nh3_shunt + nh3_gut;
dxdt_NH3_blood = 0.02*(nh3_target - NH3_blood) - K_NH3_ELIM*NH3_blood*0.01;
if(NH3_blood<10.0 && dxdt_NH3_blood<0) dxdt_NH3_blood=0;
double he_drive = (NH3_blood>HE_THRESH) ? HE_K*(NH3_blood-HE_THRESH)*(4.0-HE_GRADE) : 0.0;
double he_resolution = 0.05*HE_GRADE*(1.0-NH3_blood/(NH3_blood+HE_THRESH));
dxdt_HE_GRADE = he_drive - he_resolution;
if(HE_GRADE<0.0 && dxdt_HE_GRADE<0) dxdt_HE_GRADE=0;
if(HE_GRADE>4.0 && dxdt_HE_GRADE>0) dxdt_HE_GRADE=0;

$CAPTURE PROP_CONC SPIRO_CONC TERL_CONC_EFF RIFAX_CONC HR_obs
         HVPG_calc ALB_calc BILI_calc INR_calc CHILD_PUGH MELD_score
         MELD_Na_score VAR_BLEED_RISK MORT_1YR Na_approx

$INIT
FIBRO=0.55 HSC_ACT=0.30 HVPG=10.0 ALB=3.1 BILIRUBIN=2.5
INR_val=1.6 ASCITES=2.0 GFR_est=65.0 NH3_blood=65.0
HE_GRADE=0.5 CREAT=1.2 ALDO_ACT=1.8
PROP_GUT=0 PROP_C1=0 PROP_C2=0 SPIRO_GUT=0 SPIRO_C1=0
TERL_C=0 TERL_CE=0 RIFAX_GUT=0
'

mod_lc <- mcode("lc_shiny", lc_model_code, quiet = TRUE)

##############################################################################
# HELPER FUNCTIONS
##############################################################################

run_simulation <- function(fibro_init, hvpg_init, alb_init, bili_init,
                            inr_init, ascites_init, gfr_init, nh3_init,
                            he_init, creat_init,
                            prop_dose, prop_freq, spiro_dose,
                            rifax_on, af_conc, sim_days) {
    init_vals <- c(
        FIBRO = fibro_init, HSC_ACT = fibro_init * 0.4,
        HVPG = hvpg_init, ALB = alb_init, BILIRUBIN = bili_init,
        INR_val = inr_init, ASCITES = ascites_init, GFR_est = gfr_init,
        NH3_blood = nh3_init, HE_GRADE = he_init, CREAT = creat_init,
        ALDO_ACT = 1.0 + fibro_init,
        PROP_GUT=0, PROP_C1=0, PROP_C2=0, SPIRO_GUT=0, SPIRO_C1=0,
        TERL_C=0, TERL_CE=0, RIFAX_GUT=0
    )

    event_list <- list()

    if (prop_dose > 0) {
        prop_amt <- prop_dose * 1e3 * 0.25
        prop_times <- seq(0, (sim_days - 1) * 24, by = prop_freq)
        event_list[["prop"]] <- data.frame(
            time = prop_times, cmt = 1, amt = prop_amt
        )
    }
    if (spiro_dose > 0) {
        spiro_amt <- spiro_dose * 1e3 * 0.90
        spiro_times <- seq(0, (sim_days - 1) * 24, by = 24)
        event_list[["spiro"]] <- data.frame(
            time = spiro_times, cmt = 4, amt = spiro_amt
        )
    }
    if (rifax_on) {
        rifax_times <- seq(0, (sim_days - 1) * 24, by = 12)
        event_list[["rifax"]] <- data.frame(
            time = rifax_times, cmt = 8, amt = 550e3
        )
    }

    if (length(event_list) > 0) {
        ev_df <- do.call(rbind, event_list)
        events_obj <- ev(time = ev_df$time, cmt = ev_df$cmt, amt = ev_df$amt)
    } else {
        events_obj <- ev(time = 0, cmt = 1, amt = 0)
    }

    out <- mod_lc %>%
        init(init_vals) %>%
        param(AF_DRUG_CONC = af_conc) %>%
        mrgsim(events = events_obj, end = sim_days * 24, delta = 24) %>%
        as_tibble() %>%
        mutate(time_days = time / 24)

    return(out)
}

get_cp_class <- function(cp_score) {
    if (cp_score <= 6) "A (5-6)"
    else if (cp_score <= 9) "B (7-9)"
    else "C (10-15)"
}

meld_color <- function(meld) {
    if (meld < 10) "#27AE60"
    else if (meld < 20) "#F39C12"
    else if (meld < 30) "#E67E22"
    else "#C0392B"
}

##############################################################################
# UI
##############################################################################

ui <- dashboardPage(
    skin = "blue",

    dashboardHeader(
        title = "Liver Cirrhosis QSP Model",
        titleWidth = 300
    ),

    dashboardSidebar(
        width = 300,
        sidebarMenu(
            id = "sidebar_menu",
            menuItem("Patient Profile", tabName = "tab_patient",
                     icon = icon("user-md")),
            menuItem("Drug PK Profiles", tabName = "tab_pk",
                     icon = icon("pills")),
            menuItem("Portal Hypertension", tabName = "tab_portal",
                     icon = icon("heartbeat")),
            menuItem("Hepatic Function & Scoring", tabName = "tab_hepatic",
                     icon = icon("liver")),
            menuItem("Complications", tabName = "tab_complications",
                     icon = icon("exclamation-triangle")),
            menuItem("Scenario Comparison", tabName = "tab_scenario",
                     icon = icon("chart-line")),
            menuItem("Antifibrotic Therapy", tabName = "tab_antifib",
                     icon = icon("flask")),
            menuItem("Biomarker Panel", tabName = "tab_biomarker",
                     icon = icon("vials"))
        ),
        hr(),
        h4("Patient Parameters", style = "color:white; padding-left:15px;"),
        sliderInput("fibro_init", "Fibrosis Stage (F0–F4)",
                    min = 0.05, max = 0.99, value = 0.55, step = 0.05,
                    post = " [0-1]"),
        sliderInput("hvpg_init", "Initial HVPG (mmHg)",
                    min = 4, max = 25, value = 10, step = 1),
        sliderInput("alb_init", "Albumin (g/dL)",
                    min = 1.5, max = 4.5, value = 3.1, step = 0.1),
        sliderInput("bili_init", "Bilirubin (mg/dL)",
                    min = 0.5, max = 15, value = 2.5, step = 0.5),
        sliderInput("inr_init", "INR",
                    min = 1.0, max = 4.0, value = 1.6, step = 0.1),
        sliderInput("ascites_init", "Ascites Volume (L)",
                    min = 0, max = 15, value = 2, step = 0.5),
        sliderInput("gfr_init", "eGFR (mL/min/1.73m2)",
                    min = 10, max = 120, value = 65, step = 5),
        sliderInput("nh3_init", "Blood Ammonia (μmol/L)",
                    min = 10, max = 200, value = 65, step = 5),
        hr(),
        h4("Treatment", style = "color:white; padding-left:15px;"),
        sliderInput("prop_dose", "Propranolol Dose (mg)",
                    min = 0, max = 160, value = 40, step = 20),
        radioButtons("prop_freq", "Propranolol Frequency",
                     choices = c("BID (q12h)" = 12, "TID (q8h)" = 8,
                                 "QD (q24h)" = 24),
                     selected = 12, inline = TRUE),
        sliderInput("spiro_dose", "Spironolactone Dose (mg/day)",
                    min = 0, max = 400, value = 100, step = 50),
        checkboxInput("rifax_on", "Rifaximin 550 mg BID", value = FALSE),
        sliderInput("af_conc", "Antifibrotic Drug Level (0=none, 2=max)",
                    min = 0, max = 3, value = 0, step = 0.5),
        sliderInput("sim_days", "Simulation Duration (days)",
                    min = 30, max = 1460, value = 365, step = 30),
        actionButton("run_sim", "Run Simulation",
                     class = "btn-success btn-block",
                     icon = icon("play"))
    ),

    dashboardBody(
        tags$head(
            tags$style(HTML("
                .content-wrapper { background-color: #f4f6f9; }
                .box { border-radius: 8px; }
                .value-box { border-radius: 8px; }
                .info-box { border-radius: 8px; }
                .nav-tabs-custom > .nav-tabs > li.active {
                    border-top-color: #3c8dbc; }
                h4 { color: #2c3e50; font-weight: bold; }
            "))
        ),

        tabItems(

            # ============================================================
            # TAB 1: PATIENT PROFILE
            # ============================================================
            tabItem(tabName = "tab_patient",
                fluidRow(
                    valueBoxOutput("vbox_meld", width = 3),
                    valueBoxOutput("vbox_cp", width = 3),
                    valueBoxOutput("vbox_hvpg", width = 3),
                    valueBoxOutput("vbox_mort", width = 3)
                ),
                fluidRow(
                    box(title = "Disease Staging Summary", width = 6,
                        status = "primary", solidHeader = TRUE,
                        DTOutput("tbl_staging")),
                    box(title = "MELD Score Components", width = 6,
                        status = "warning", solidHeader = TRUE,
                        plotlyOutput("plot_meld_gauge", height = "300px"))
                ),
                fluidRow(
                    box(title = "Annual Risk Assessment", width = 12,
                        status = "danger", solidHeader = TRUE,
                        plotlyOutput("plot_risk_bar", height = "200px"))
                )
            ),

            # ============================================================
            # TAB 2: PK PROFILES
            # ============================================================
            tabItem(tabName = "tab_pk",
                fluidRow(
                    box(title = "Propranolol PK: Plasma Concentration",
                        width = 6, status = "info", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_prop_pk", height = "300px"))),
                    box(title = "Spironolactone PK: Plasma Concentration",
                        width = 6, status = "info", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_spiro_pk", height = "300px")))
                ),
                fluidRow(
                    box(title = "Propranolol PD: Heart Rate Reduction",
                        width = 6, status = "success", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_prop_pd", height = "300px"))),
                    box(title = "Rifaximin: Gut Concentration & NH3 Effect",
                        width = 6, status = "success", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_rifax_pk", height = "300px")))
                ),
                fluidRow(
                    box(title = "PK Parameter Summary", width = 12,
                        status = "primary", solidHeader = TRUE,
                        DTOutput("tbl_pk_summary"))
                )
            ),

            # ============================================================
            # TAB 3: PORTAL HYPERTENSION
            # ============================================================
            tabItem(tabName = "tab_portal",
                fluidRow(
                    infoBoxOutput("ibox_hvpg_class", width = 4),
                    infoBoxOutput("ibox_varices", width = 4),
                    infoBoxOutput("ibox_bleed_risk", width = 4)
                ),
                fluidRow(
                    box(title = "HVPG Trajectory", width = 8,
                        status = "primary", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_hvpg_traj", height = "350px"))),
                    box(title = "Portal Pressure Determinants", width = 4,
                        status = "info", solidHeader = TRUE,
                        plotlyOutput("plot_hvpg_pie", height = "350px"))
                ),
                fluidRow(
                    box(title = "Variceal Bleeding Risk Over Time", width = 12,
                        status = "danger", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_bleed_risk", height = "280px")))
                )
            ),

            # ============================================================
            # TAB 4: HEPATIC FUNCTION & SCORING
            # ============================================================
            tabItem(tabName = "tab_hepatic",
                fluidRow(
                    box(title = "MELD Score over Time", width = 6,
                        status = "warning", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_meld_traj", height = "300px"))),
                    box(title = "Child-Pugh Score over Time", width = 6,
                        status = "warning", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_cp_traj", height = "300px")))
                ),
                fluidRow(
                    box(title = "Albumin & Bilirubin", width = 6,
                        status = "primary", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_alb_bili", height = "300px"))),
                    box(title = "INR & Hepatic Reserve", width = 6,
                        status = "primary", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_inr", height = "300px")))
                )
            ),

            # ============================================================
            # TAB 5: COMPLICATIONS
            # ============================================================
            tabItem(tabName = "tab_complications",
                fluidRow(
                    box(title = "Ascites Volume", width = 6,
                        status = "info", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_ascites", height = "300px"))),
                    box(title = "Renal Function (eGFR & Creatinine)", width = 6,
                        status = "danger", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_renal", height = "300px")))
                ),
                fluidRow(
                    box(title = "Blood Ammonia & HE Grade", width = 6,
                        status = "warning", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_he", height = "300px"))),
                    box(title = "Complication Risk Timeline", width = 6,
                        status = "danger", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_comp_risk", height = "300px")))
                )
            ),

            # ============================================================
            # TAB 6: SCENARIO COMPARISON
            # ============================================================
            tabItem(tabName = "tab_scenario",
                fluidRow(
                    box(title = "Scenario Configuration", width = 12,
                        status = "primary", solidHeader = TRUE,
                        p("This tab compares the current treatment setting against natural history and standard-of-care combinations."),
                        actionButton("run_scenarios", "Compare All Scenarios",
                                     class = "btn-primary", icon = icon("balance-scale"))
                    )
                ),
                fluidRow(
                    box(title = "MELD Score — Scenario Comparison", width = 6,
                        status = "warning", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_scen_meld", height = "320px"))),
                    box(title = "Fibrosis Index — Scenario Comparison", width = 6,
                        status = "danger", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_scen_fibro", height = "320px")))
                ),
                fluidRow(
                    box(title = "HVPG — Scenario Comparison", width = 6,
                        status = "info", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_scen_hvpg", height = "320px"))),
                    box(title = "1-Year Mortality Risk — Scenario Comparison", width = 6,
                        status = "danger", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_scen_mort", height = "320px")))
                ),
                fluidRow(
                    box(title = "Summary Table — End of Simulation", width = 12,
                        status = "primary", solidHeader = TRUE,
                        DTOutput("tbl_scenarios"))
                )
            ),

            # ============================================================
            # TAB 7: ANTIFIBROTIC THERAPY
            # ============================================================
            tabItem(tabName = "tab_antifib",
                fluidRow(
                    box(title = "Fibrosis Kinetics (HSC Activation & Fibrosis Index)",
                        width = 8, status = "success", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_fibro_kinetics", height = "380px"))),
                    box(title = "Antifibrotic Target Pathways", width = 4,
                        status = "info", solidHeader = TRUE,
                        tags$div(
                            tags$h4("TGF-β / HSC Pathway"),
                            tags$ul(
                                tags$li("TGF-β1 → Smad2/3 → Collagen I/III"),
                                tags$li("PDGF-BB → HSC proliferation"),
                                tags$li("TIMP-1 → MMP inhibition → fibrosis ↑")
                            ),
                            tags$h4("Drug Targets (Investigational)"),
                            tags$ul(
                                tags$li(tags$b("FXR Agonist (OCA):"), " → SHP → TGF-β↓"),
                                tags$li(tags$b("CCR2/CCR5 Inhibitor:"), " → Monocyte↓"),
                                tags$li(tags$b("LOXL2 Inhibitor:"), " → Collagen crosslink↓"),
                                tags$li(tags$b("ASK1 Inhibitor:"), " → ROS/oxidative stress↓"),
                                tags$li(tags$b("THRβ Agonist (Resmetirom):"), " → NASH→fibrosis↓"),
                                tags$li(tags$b("GLP-1 Agonist (Semaglutide):"), " → Steatosis↓")
                            )
                        ))
                ),
                fluidRow(
                    box(title = "HVPG Response to Antifibrotic Treatment",
                        width = 6, status = "primary", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_antifib_hvpg", height = "280px"))),
                    box(title = "Hepatic Function Recovery Potential",
                        width = 6, status = "warning", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_antifib_func", height = "280px")))
                )
            ),

            # ============================================================
            # TAB 8: BIOMARKER PANEL
            # ============================================================
            tabItem(tabName = "tab_biomarker",
                fluidRow(
                    box(title = "Serum Biomarker Heatmap Over Time", width = 12,
                        status = "primary", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_biomarker_heat", height = "400px")))
                ),
                fluidRow(
                    box(title = "MELD vs HVPG Correlation", width = 6,
                        status = "info", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_meld_hvpg_corr", height = "300px"))),
                    box(title = "Albumin vs Fibrosis Stage", width = 6,
                        status = "warning", solidHeader = TRUE,
                        withSpinner(plotlyOutput("plot_alb_fibro", height = "300px")))
                )
            )
        )
    )
)

##############################################################################
# SERVER
##############################################################################

server <- function(input, output, session) {

    # Reactive simulation
    sim_out <- eventReactive(input$run_sim, {
        withProgress(message = "Running QSP Simulation...", value = 0, {
            incProgress(0.3)
            out <- run_simulation(
                fibro_init = input$fibro_init,
                hvpg_init  = input$hvpg_init,
                alb_init   = input$alb_init,
                bili_init  = input$bili_init,
                inr_init   = input$inr_init,
                ascites_init = input$ascites_init,
                gfr_init   = input$gfr_init,
                nh3_init   = input$nh3_init,
                he_init    = 0.5,
                creat_init = max(0.9, 0.9 * 90 / input$gfr_init),
                prop_dose  = input$prop_dose,
                prop_freq  = as.numeric(input$prop_freq),
                spiro_dose = input$spiro_dose,
                rifax_on   = input$rifax_on,
                af_conc    = input$af_conc,
                sim_days   = input$sim_days
            )
            incProgress(0.7)
            out
        })
    }, ignoreNULL = FALSE)

    # Scenario comparison reactive
    scenarios_out <- eventReactive(input$run_scenarios, {
        withProgress(message = "Running scenario comparisons...", value = 0, {
            base_args <- list(
                fibro_init = input$fibro_init, hvpg_init = input$hvpg_init,
                alb_init = input$alb_init, bili_init = input$bili_init,
                inr_init = input$inr_init, ascites_init = input$ascites_init,
                gfr_init = input$gfr_init, nh3_init = input$nh3_init,
                he_init = 0.5, creat_init = max(0.9, 0.9*90/input$gfr_init),
                sim_days = input$sim_days
            )
            incProgress(0.1)
            scen_list <- list(
                "Natural History"        = c(base_args, list(prop_dose=0, prop_freq=12, spiro_dose=0, rifax_on=FALSE, af_conc=0)),
                "Propranolol 40mg BID"   = c(base_args, list(prop_dose=40, prop_freq=12, spiro_dose=0, rifax_on=FALSE, af_conc=0)),
                "Spironolactone 100mg QD" = c(base_args, list(prop_dose=0, prop_freq=12, spiro_dose=100, rifax_on=FALSE, af_conc=0)),
                "Propranolol + Spiro"    = c(base_args, list(prop_dose=40, prop_freq=12, spiro_dose=100, rifax_on=FALSE, af_conc=0)),
                "Triple (+ Rifaximin)"   = c(base_args, list(prop_dose=40, prop_freq=12, spiro_dose=100, rifax_on=TRUE, af_conc=0)),
                "Antifibrotic (FXR Ag.)" = c(base_args, list(prop_dose=0, prop_freq=12, spiro_dose=0, rifax_on=FALSE, af_conc=1.5))
            )
            results <- lapply(names(scen_list), function(nm) {
                incProgress(0.1)
                do.call(run_simulation, scen_list[[nm]]) %>%
                    mutate(Scenario = nm)
            })
            bind_rows(results)
        })
    })

    # ---- TAB 1: VALUE BOXES ----
    output$vbox_meld <- renderValueBox({
        d <- sim_out()
        meld_now <- round(tail(d$MELD_score, 1), 1)
        col <- if(meld_now < 10) "green" else if(meld_now < 20) "yellow" else if(meld_now < 30) "orange" else "red"
        valueBox(meld_now, "MELD Score (End)", icon = icon("calculator"), color = col)
    })

    output$vbox_cp <- renderValueBox({
        d <- sim_out()
        cp <- round(tail(d$CHILD_PUGH, 1), 1)
        cp_class <- get_cp_class(cp)
        col <- if(cp <= 6) "green" else if(cp <= 9) "yellow" else "red"
        valueBox(paste0(round(cp, 0), " — ", cp_class), "Child-Pugh", icon = icon("star"), color = col)
    })

    output$vbox_hvpg <- renderValueBox({
        d <- sim_out()
        hvpg_now <- round(tail(d$HVPG_calc, 1), 1)
        col <- if(hvpg_now < 10) "green" else if(hvpg_now < 12) "yellow" else "red"
        valueBox(paste0(hvpg_now, " mmHg"), "HVPG (End)", icon = icon("tachometer-alt"), color = col)
    })

    output$vbox_mort <- renderValueBox({
        d <- sim_out()
        mort <- round(tail(d$MORT_1YR, 1) * 100, 1)
        col <- if(mort < 5) "green" else if(mort < 20) "yellow" else "red"
        valueBox(paste0(mort, "%"), "1-Year Mortality Risk", icon = icon("exclamation-circle"), color = col)
    })

    output$tbl_staging <- renderDT({
        d <- sim_out()
        last <- tail(d, 1)
        df <- data.frame(
            Parameter = c("Fibrosis Index", "HVPG (mmHg)", "Albumin (g/dL)",
                          "Bilirubin (mg/dL)", "INR", "Ascites (L)",
                          "eGFR (mL/min)", "Ammonia (μmol/L)",
                          "HE Grade", "Creatinine (mg/dL)",
                          "MELD", "MELD-Na", "Child-Pugh"),
            `Baseline` = c(input$fibro_init, input$hvpg_init, input$alb_init,
                           input$bili_init, input$inr_init, input$ascites_init,
                           input$gfr_init, input$nh3_init, 0.5,
                           max(0.9, 0.9*90/input$gfr_init),
                           round(d$MELD_score[1], 1),
                           round(d$MELD_Na_score[1], 1),
                           round(d$CHILD_PUGH[1], 1)),
            `End of Simulation` = c(
                round(last$FIBRO, 2), round(last$HVPG_calc, 1),
                round(last$ALB_calc, 2), round(last$BILI_calc, 2),
                round(last$INR_calc, 2), round(last$ASCITES, 1),
                round(last$GFR_est, 0), round(last$NH3_blood, 0),
                round(last$HE_GRADE, 1), round(last$CREAT, 2),
                round(last$MELD_score, 1), round(last$MELD_Na_score, 1),
                round(last$CHILD_PUGH, 1)
            )
        )
        datatable(df, options = list(pageLength = 13, dom = 't'),
                  rownames = FALSE) %>%
            formatStyle('End of Simulation',
                        color = styleInterval(c(0), c('black', 'black')))
    })

    output$plot_meld_gauge <- renderPlotly({
        d <- sim_out()
        meld_val <- round(tail(d$MELD_score, 1), 1)
        plot_ly(
            type = "indicator", mode = "gauge+number+delta",
            value = meld_val,
            delta = list(reference = d$MELD_score[1]),
            gauge = list(
                axis = list(range = list(6, 40)),
                bar = list(color = meld_color(meld_val)),
                steps = list(
                    list(range = c(6, 10), color = "#E8F5E9"),
                    list(range = c(10, 20), color = "#FFF9C4"),
                    list(range = c(20, 30), color = "#FFF3E0"),
                    list(range = c(30, 40), color = "#FFEBEE")
                ),
                threshold = list(line = list(color = "purple", width = 3),
                                 thickness = 0.75, value = 15)
            ),
            title = list(text = "MELD Score")
        ) %>% layout(margin = list(t = 50))
    })

    output$plot_risk_bar <- renderPlotly({
        d <- sim_out()
        last <- tail(d, 1)
        risks <- data.frame(
            Risk = c("Annual Variceal Bleeding", "1-Year Mortality",
                     "Decompensation (90d est.)", "HRS Risk (1yr est.)"),
            Probability = c(
                min(1, last$VAR_BLEED_RISK),
                min(1, last$MORT_1YR),
                min(1, last$FIBRO * 0.4 + (last$HVPG_calc > 12) * 0.1),
                min(1, pmax(0, (last$HVPG_calc - 16) * 0.05))
            )
        )
        plot_ly(risks, x = ~Probability, y = ~Risk, type = "bar",
                orientation = "h",
                marker = list(color = c("#F44336", "#9C27B0", "#FF9800", "#2196F3"))) %>%
            layout(xaxis = list(title = "Probability", tickformat = ".0%"),
                   yaxis = list(title = ""),
                   margin = list(l = 200))
    })

    # ---- TAB 2: PK PLOTS ----
    output$plot_prop_pk <- renderPlotly({
        d <- sim_out()
        plot_ly(d, x = ~time_days, y = ~PROP_CONC, type = "scatter",
                mode = "lines", line = list(color = "#1565C0", width = 2)) %>%
            layout(xaxis = list(title = "Time (days)"),
                   yaxis = list(title = "Propranolol Cplasma (ng/mL)"),
                   title = "")
    })

    output$plot_spiro_pk <- renderPlotly({
        d <- sim_out()
        plot_ly(d, x = ~time_days, y = ~SPIRO_CONC, type = "scatter",
                mode = "lines", line = list(color = "#6A1B9A", width = 2)) %>%
            layout(xaxis = list(title = "Time (days)"),
                   yaxis = list(title = "Spironolactone Cplasma (ng/mL)"),
                   title = "")
    })

    output$plot_prop_pd <- renderPlotly({
        d <- sim_out()
        plot_ly(d, x = ~time_days, y = ~HR_obs, type = "scatter",
                mode = "lines", line = list(color = "#D32F2F", width = 2)) %>%
            add_hline(y = 55, line = list(dash = "dash", color = "orange"),
                      annotation_text = "Target ≤60 bpm") %>%
            layout(xaxis = list(title = "Time (days)"),
                   yaxis = list(title = "Heart Rate (bpm)"))
    })

    output$plot_rifax_pk <- renderPlotly({
        d <- sim_out()
        p <- plot_ly(d, x = ~time_days)
        if (input$rifax_on) {
            p <- p %>% add_lines(y = ~RIFAX_CONC, name = "Gut Conc (ng/mL eq.)",
                                  line = list(color = "#4A148C", width = 2))
        }
        p %>% add_lines(y = ~NH3_blood, name = "Ammonia (μmol/L)",
                         line = list(color = "#FF9800", width = 2),
                         yaxis = "y2") %>%
            layout(
                xaxis = list(title = "Time (days)"),
                yaxis = list(title = "Rifaximin Gut Conc"),
                yaxis2 = list(title = "Ammonia (μmol/L)", overlaying = "y",
                               side = "right")
            )
    })

    output$tbl_pk_summary <- renderDT({
        df <- data.frame(
            Drug = c("Propranolol", "Spironolactone", "Terlipressin", "Rifaximin"),
            `PK Model` = c("2-compartment", "1-compartment", "Effect-compartment", "Gut-lumen"),
            `F (%)` = c("25", "90", "IV (100)", "<0.4"),
            `CL (L/h)` = c("30.0", "3.5", "15.0", "N/A"),
            `Vd (L)` = c("450", "70", "30", "Gut only"),
            `Half-life` = c("~5h (cirrhosis↑)", "~20h (active met)", "~6h", "~5h (gut)"),
            `Primary Action` = c("HR↓ / HVPG↓ (β1,β2-block)",
                                  "Aldosterone block / Ascites↓",
                                  "Splanchnic vasoconstriction / HRS",
                                  "Gut NH3↓ / HE prevention"),
            stringsAsFactors = FALSE
        )
        datatable(df, options = list(dom = 't', pageLength = 4), rownames = FALSE)
    })

    # ---- TAB 3: PORTAL HYPERTENSION ----
    output$ibox_hvpg_class <- renderInfoBox({
        d <- sim_out()
        hvpg <- round(tail(d$HVPG_calc, 1), 1)
        cls <- if(hvpg < 5) "Normal" else if(hvpg < 10) "Subclinical PHT" else if(hvpg < 12) "Clinical PHT" else "High-risk PHT"
        col <- if(hvpg < 5) "green" else if(hvpg < 10) "yellow" else if(hvpg < 12) "orange" else "red"
        infoBox("HVPG Classification", paste0(hvpg, " mmHg — ", cls),
                icon = icon("thermometer-half"), color = col)
    })

    output$ibox_varices <- renderInfoBox({
        d <- sim_out()
        hvpg <- tail(d$HVPG_calc, 1)
        msg <- if(hvpg < 10) "Varices unlikely" else if(hvpg < 12) "Small varices possible" else "Large varices — prophylaxis needed"
        col <- if(hvpg < 10) "green" else if(hvpg < 12) "yellow" else "red"
        infoBox("Variceal Status", msg, icon = icon("water"), color = col)
    })

    output$ibox_bleed_risk <- renderInfoBox({
        d <- sim_out()
        risk <- round(tail(d$VAR_BLEED_RISK, 1) * 100, 1)
        col <- if(risk < 5) "green" else if(risk < 15) "yellow" else "red"
        infoBox("Annual Bleeding Risk", paste0(risk, "%"),
                icon = icon("tint"), color = col)
    })

    output$plot_hvpg_traj <- renderPlotly({
        d <- sim_out()
        plot_ly(d, x = ~time_days, y = ~HVPG_calc, type = "scatter",
                mode = "lines", line = list(color = "#1565C0", width = 2.5),
                name = "HVPG") %>%
            add_hline(y = 10, line = list(dash = "dot", color = "#FF9800"),
                      annotation_text = "10 mmHg (ascites)") %>%
            add_hline(y = 12, line = list(dash = "dot", color = "#F44336"),
                      annotation_text = "12 mmHg (bleeding)") %>%
            layout(xaxis = list(title = "Time (days)"),
                   yaxis = list(title = "HVPG (mmHg)"))
    })

    output$plot_hvpg_pie <- renderPlotly({
        fib_contrib <- input$hvpg_init * 0.6
        et1_contrib <- input$hvpg_init * 0.25
        drug_reduce <- input$hvpg_init * 0.15
        plot_ly(labels = c("Fibrosis/Resistance", "ET-1/Vasoconstriction", "Drug Reduction"),
                values = c(fib_contrib, et1_contrib, drug_reduce),
                type = "pie",
                marker = list(colors = c("#D32F2F", "#FF9800", "#2196F3"))) %>%
            layout(title = "HVPG Determinants")
    })

    output$plot_bleed_risk <- renderPlotly({
        d <- sim_out()
        plot_ly(d, x = ~time_days, y = ~VAR_BLEED_RISK, type = "scatter",
                mode = "lines", fill = "tozeroy",
                line = list(color = "#C62828"),
                fillcolor = "rgba(198, 40, 40, 0.2)") %>%
            layout(xaxis = list(title = "Time (days)"),
                   yaxis = list(title = "Annual Variceal Bleeding Risk",
                                tickformat = ".0%"))
    })

    # ---- TAB 4: HEPATIC FUNCTION ----
    output$plot_meld_traj <- renderPlotly({
        d <- sim_out()
        plot_ly(d, x = ~time_days, y = ~MELD_score, type = "scatter",
                mode = "lines", line = list(color = "#FF6F00", width = 2),
                name = "MELD") %>%
            add_lines(y = ~MELD_Na_score, line = list(color = "#7B1FA2", dash = "dash"),
                      name = "MELD-Na") %>%
            add_hline(y = 15, line = list(dash = "dot", color = "purple"),
                      annotation_text = "MELD=15 (transplant)") %>%
            layout(xaxis = list(title = "Time (days)"),
                   yaxis = list(title = "Score"))
    })

    output$plot_cp_traj <- renderPlotly({
        d <- sim_out()
        plot_ly(d, x = ~time_days, y = ~CHILD_PUGH, type = "scatter",
                mode = "lines", line = list(color = "#1B5E20", width = 2)) %>%
            add_hline(y = 7, line = list(dash = "dot", color = "#FF9800"),
                      annotation_text = "CP Class B") %>%
            add_hline(y = 10, line = list(dash = "dot", color = "#F44336"),
                      annotation_text = "CP Class C") %>%
            layout(xaxis = list(title = "Time (days)"),
                   yaxis = list(title = "Child-Pugh Score"))
    })

    output$plot_alb_bili <- renderPlotly({
        d <- sim_out()
        plot_ly(d, x = ~time_days) %>%
            add_lines(y = ~ALB_calc, name = "Albumin (g/dL)",
                      line = list(color = "#1565C0", width = 2)) %>%
            add_lines(y = ~BILI_calc, name = "Bilirubin (mg/dL)",
                      line = list(color = "#F57F17", width = 2),
                      yaxis = "y2") %>%
            layout(
                xaxis = list(title = "Time (days)"),
                yaxis = list(title = "Albumin (g/dL)", range = c(0, 5)),
                yaxis2 = list(title = "Bilirubin (mg/dL)", overlaying = "y",
                               side = "right", range = c(0, 20))
            )
    })

    output$plot_inr <- renderPlotly({
        d <- sim_out()
        plot_ly(d, x = ~time_days, y = ~INR_calc, type = "scatter",
                mode = "lines", line = list(color = "#BF360C", width = 2)) %>%
            add_hline(y = 1.7, line = list(dash = "dot", color = "#FF9800"),
                      annotation_text = "INR=1.7") %>%
            add_hline(y = 2.3, line = list(dash = "dot", color = "#F44336"),
                      annotation_text = "INR=2.3") %>%
            layout(xaxis = list(title = "Time (days)"),
                   yaxis = list(title = "INR"))
    })

    # ---- TAB 5: COMPLICATIONS ----
    output$plot_ascites <- renderPlotly({
        d <- sim_out()
        plot_ly(d, x = ~time_days, y = ~ASCITES, type = "scatter",
                mode = "lines", fill = "tozeroy",
                line = list(color = "#0288D1"),
                fillcolor = "rgba(2, 136, 209, 0.2)") %>%
            add_hline(y = 5, line = list(dash = "dot", color = "#FF9800"),
                      annotation_text = "Moderate (5L)") %>%
            layout(xaxis = list(title = "Time (days)"),
                   yaxis = list(title = "Ascites Volume (L)"))
    })

    output$plot_renal <- renderPlotly({
        d <- sim_out()
        plot_ly(d, x = ~time_days) %>%
            add_lines(y = ~GFR_est, name = "eGFR (mL/min)",
                      line = list(color = "#1565C0", width = 2)) %>%
            add_lines(y = ~CREAT, name = "Creatinine (mg/dL)",
                      line = list(color = "#D32F2F", width = 2),
                      yaxis = "y2") %>%
            add_hline(y = 60, line = list(dash = "dot", color = "#FF9800"),
                      annotation_text = "GFR=60 (CKD G2)") %>%
            layout(
                xaxis = list(title = "Time (days)"),
                yaxis = list(title = "eGFR (mL/min)"),
                yaxis2 = list(title = "Creatinine (mg/dL)", overlaying = "y",
                               side = "right")
            )
    })

    output$plot_he <- renderPlotly({
        d <- sim_out()
        plot_ly(d, x = ~time_days) %>%
            add_lines(y = ~NH3_blood, name = "Ammonia (μmol/L)",
                      line = list(color = "#7B1FA2", width = 2)) %>%
            add_lines(y = ~HE_GRADE * 20, name = "HE Grade (×20)",
                      line = list(color = "#F57F17", dash = "dash", width = 2),
                      yaxis = "y2") %>%
            add_hline(y = 50, line = list(dash = "dot", color = "#F44336"),
                      annotation_text = "NH3=50 HE threshold") %>%
            layout(
                xaxis = list(title = "Time (days)"),
                yaxis = list(title = "Ammonia (μmol/L)"),
                yaxis2 = list(title = "HE Grade (×20 scaled)", overlaying = "y",
                               side = "right", range = c(0, 80))
            )
    })

    output$plot_comp_risk <- renderPlotly({
        d <- sim_out()
        hrs_risk <- pmax(0, (d$HVPG_calc - 16) * 0.05)
        sbp_risk <- pmax(0, d$FIBRO * 0.2 * (d$HVPG_calc > 12))
        plot_ly(d, x = ~time_days) %>%
            add_lines(y = ~VAR_BLEED_RISK, name = "Variceal Bleeding",
                      line = list(color = "#F44336", width = 2)) %>%
            add_lines(y = ~MORT_1YR, name = "1-yr Mortality",
                      line = list(color = "#9C27B0", width = 2)) %>%
            add_lines(y = hrs_risk, name = "HRS Risk",
                      line = list(color = "#FF9800", width = 2)) %>%
            add_lines(y = sbp_risk, name = "SBP Risk",
                      line = list(color = "#2196F3", width = 2)) %>%
            layout(xaxis = list(title = "Time (days)"),
                   yaxis = list(title = "Risk (probability)", tickformat = ".0%"))
    })

    # ---- TAB 6: SCENARIO COMPARISON ----
    output$plot_scen_meld <- renderPlotly({
        req(scenarios_out())
        d <- scenarios_out()
        colors <- c("#D32F2F","#1565C0","#2E7D32","#E65100","#6A1B9A","#00838F")
        plot_ly(d, x = ~time_days, y = ~MELD_score, color = ~Scenario,
                type = "scatter", mode = "lines",
                colors = colors) %>%
            add_hline(y = 15, line = list(dash = "dot", color = "black")) %>%
            layout(xaxis = list(title = "Time (days)"),
                   yaxis = list(title = "MELD Score"),
                   legend = list(orientation = "h", y = -0.2))
    })

    output$plot_scen_fibro <- renderPlotly({
        req(scenarios_out())
        d <- scenarios_out()
        colors <- c("#D32F2F","#1565C0","#2E7D32","#E65100","#6A1B9A","#00838F")
        plot_ly(d, x = ~time_days, y = ~FIBRO, color = ~Scenario,
                type = "scatter", mode = "lines", colors = colors) %>%
            layout(xaxis = list(title = "Time (days)"),
                   yaxis = list(title = "Fibrosis Index (0-1)"),
                   legend = list(orientation = "h", y = -0.2))
    })

    output$plot_scen_hvpg <- renderPlotly({
        req(scenarios_out())
        d <- scenarios_out()
        colors <- c("#D32F2F","#1565C0","#2E7D32","#E65100","#6A1B9A","#00838F")
        plot_ly(d, x = ~time_days, y = ~HVPG_calc, color = ~Scenario,
                type = "scatter", mode = "lines", colors = colors) %>%
            add_hline(y = 12, line = list(dash = "dot", color = "red")) %>%
            layout(xaxis = list(title = "Time (days)"),
                   yaxis = list(title = "HVPG (mmHg)"),
                   legend = list(orientation = "h", y = -0.2))
    })

    output$plot_scen_mort <- renderPlotly({
        req(scenarios_out())
        d <- scenarios_out()
        colors <- c("#D32F2F","#1565C0","#2E7D32","#E65100","#6A1B9A","#00838F")
        plot_ly(d, x = ~time_days, y = ~MORT_1YR, color = ~Scenario,
                type = "scatter", mode = "lines", colors = colors) %>%
            layout(xaxis = list(title = "Time (days)"),
                   yaxis = list(title = "1-Year Mortality Risk", tickformat = ".0%"),
                   legend = list(orientation = "h", y = -0.2))
    })

    output$tbl_scenarios <- renderDT({
        req(scenarios_out())
        d <- scenarios_out()
        d %>%
            group_by(Scenario) %>%
            filter(time_days == max(time_days)) %>%
            summarise(
                `MELD (end)` = round(last(MELD_score), 1),
                `HVPG (end)` = round(last(HVPG_calc), 1),
                `Fibrosis (end)` = round(last(FIBRO), 3),
                `Albumin (end)` = round(last(ALB_calc), 2),
                `Ascites (end, L)` = round(last(ASCITES), 1),
                `GFR (end)` = round(last(GFR_est), 0),
                `1yr Mortality` = scales::percent(last(MORT_1YR), accuracy = 0.1),
                .groups = "drop"
            ) %>%
            datatable(options = list(dom = 't', pageLength = 10), rownames = FALSE)
    })

    # ---- TAB 7: ANTIFIBROTIC ----
    output$plot_fibro_kinetics <- renderPlotly({
        d <- sim_out()
        plot_ly(d, x = ~time_days) %>%
            add_lines(y = ~FIBRO, name = "Fibrosis Index",
                      line = list(color = "#F57F17", width = 2.5)) %>%
            add_lines(y = ~HSC_ACT, name = "HSC Activation",
                      line = list(color = "#2E7D32", width = 2, dash = "dash")) %>%
            add_hline(y = 0.75, line = list(dash = "dot", color = "red"),
                      annotation_text = "Cirrhosis threshold") %>%
            layout(xaxis = list(title = "Time (days)"),
                   yaxis = list(title = "Index (0-1)"),
                   legend = list(orientation = "h"))
    })

    output$plot_antifib_hvpg <- renderPlotly({
        d <- sim_out()
        plot_ly(d, x = ~time_days, y = ~HVPG_calc, type = "scatter",
                mode = "lines", line = list(color = "#1565C0", width = 2),
                name = paste0("AF level = ", input$af_conc)) %>%
            add_hline(y = 12, line = list(dash = "dot", color = "red")) %>%
            layout(xaxis = list(title = "Time (days)"),
                   yaxis = list(title = "HVPG (mmHg)"))
    })

    output$plot_antifib_func <- renderPlotly({
        d <- sim_out()
        plot_ly(d, x = ~time_days) %>%
            add_lines(y = ~ALB_calc, name = "Albumin (g/dL)",
                      line = list(color = "#1565C0")) %>%
            add_lines(y = ~BILI_calc * 0.5, name = "Bilirubin/2 (mg/dL)",
                      line = list(color = "#FF9800", dash = "dash")) %>%
            add_lines(y = ~INR_calc, name = "INR",
                      line = list(color = "#D32F2F")) %>%
            layout(xaxis = list(title = "Time (days)"),
                   yaxis = list(title = "Value"))
    })

    # ---- TAB 8: BIOMARKER PANEL ----
    output$plot_biomarker_heat <- renderPlotly({
        d <- sim_out()
        t_pts <- seq(1, nrow(d), length.out = min(50, nrow(d)))
        d_sub <- d[round(t_pts), ]
        biomarkers <- c("FIBRO", "HVPG_calc", "ALB_calc", "BILI_calc",
                        "INR_calc", "ASCITES", "GFR_est", "NH3_blood",
                        "CREAT", "MELD_score")
        mat <- as.matrix(d_sub[, biomarkers])
        mat_scaled <- scale(mat)
        plot_ly(z = t(mat_scaled), x = round(d_sub$time_days, 0),
                y = biomarkers, type = "heatmap",
                colorscale = "RdBu", reversescale = TRUE) %>%
            layout(xaxis = list(title = "Time (days)"),
                   yaxis = list(title = "Biomarker"),
                   title = "Z-score normalized biomarker heatmap")
    })

    output$plot_meld_hvpg_corr <- renderPlotly({
        d <- sim_out()
        plot_ly(d, x = ~HVPG_calc, y = ~MELD_score,
                type = "scatter", mode = "markers",
                marker = list(color = ~time_days, colorscale = "Viridis",
                               showscale = TRUE, size = 6),
                text = ~paste("Day:", round(time_days))) %>%
            layout(xaxis = list(title = "HVPG (mmHg)"),
                   yaxis = list(title = "MELD Score"))
    })

    output$plot_alb_fibro <- renderPlotly({
        d <- sim_out()
        plot_ly(d, x = ~FIBRO, y = ~ALB_calc,
                type = "scatter", mode = "markers",
                marker = list(color = ~MELD_score, colorscale = "Reds",
                               showscale = TRUE, size = 6),
                text = ~paste("Day:", round(time_days))) %>%
            layout(xaxis = list(title = "Fibrosis Index"),
                   yaxis = list(title = "Albumin (g/dL)"))
    })
}

##############################################################################
# RUN APP
##############################################################################

shinyApp(ui = ui, server = server)
