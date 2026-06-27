##############################################################################
# Tuberculosis QSP — Interactive Shiny Dashboard
# Tabs:
#   1. Patient Profile & Disease Stage
#   2. Drug PK (RIPE Plasma Profiles)
#   3. Bacterial Dynamics (AR / SR / NR subpopulations)
#   4. Immune Response (Macrophage, Th1, Cytokines)
#   5. Clinical Endpoints (Culture conversion, smear, outcomes)
#   6. Scenario Comparison (All 6 scenarios)
#   7. Biomarker & PD Analysis (MIC targets, PKPD)
##############################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)

# ─── Inline mrgsolve model ──────────────────────────────────────────────────
tb_code <- '
$PARAM @annotated
kgAR:0.6931:AR growth (/d)
kgSR:0.0693:SR growth (/d)
kAR2SR:0.050:AR→SR (/d)
kSR2NR:0.020:SR→NR (/d)
kNR2AR:0.0005:NR reactivation (/d)
kdAR:0.040:AR natural death
kdSR:0.010:SR natural death
kdNR:0.002:NR natural death
Bmax:1.0e9:Max CFU
UM_ss:500:Macrophage SS
kUM_p:50:Mφ production
kUM_d:0.10:Mφ death
kInf:1.0e-4:Infection rate
kActv:0.50:Mφ activation by IFNg
kIM_d:0.15:IM death
kAM_d:0.08:AM death
kAMkill:2.0:AM kill rate
kTh1_p:5.0:Th1 production
kTh1_d:0.10:Th1 death
kTh1_half:50:Bacterial half-max for Th1
kIFNg_p:2.0:IFNg production
kIFNg_d:1.50:IFNg clearance
kTNFa_p:1.0:TNFa production
kTNFa_d:2.0:TNFa clearance
kIL10_p:0.30:IL10 production
kIL10_d:1.0:IL10 clearance
kIL10_sup:0.5:IL10 suppression of IFNg
phi_Th1:1.0:Th1 scaling
phi_Mact:1.0:Mφ activation scaling
RIF_dose_mg:600:RIF dose (mg)
RIF_ka:1.50:RIF ka (/h)
RIF_V:49.0:RIF volume (L)
RIF_CL_base:9.0:RIF CL (L/h)
RIF_F:0.68:RIF F
RIF_ind_max:0.60:RIF autoinduction
RIF_ind_ec50:4.0:RIF autoinduction EC50
RIF_lung_f:0.25:RIF lung penetration
INH_dose_mg:300:INH dose (mg)
INH_ka:1.20:INH ka (/h)
INH_V:38.0:INH volume (L)
INH_CL:5.0:INH CL (L/h)
INH_F:0.90:INH F
NAT2_factor:1.0:NAT2 factor
PZA_dose_mg:1500:PZA dose (mg)
PZA_ka:1.0:PZA ka (/h)
PZA_V:36.0:PZA volume (L)
PZA_CL:2.8:PZA CL (L/h)
PZA_F:0.98:PZA F
PZA_pH_act:1.0:PZA pH activation
EMB_dose_mg:1200:EMB dose (mg)
EMB_ka:1.1:EMB ka (/h)
EMB_V:60.0:EMB volume (L)
EMB_CL:9.4:EMB CL (L/h)
EMB_F:0.78:EMB F
BDQ_dose_mg:400:BDQ loading dose
BDQ_ka:0.50:BDQ ka
BDQ_V:164000:BDQ volume
BDQ_CL:91.0:BDQ CL
BDQ_F:0.45:BDQ F
BDQ_on:0.0:BDQ flag
RIF_Emax_AR:3.50:RIF Emax AR
RIF_EC50_AR:0.60:RIF EC50 AR
RIF_Emax_SR:1.50:RIF Emax SR
RIF_EC50_SR:1.50:RIF EC50 SR
RIF_hill:1.0:RIF Hill
INH_Emax_AR:4.00:INH Emax AR
INH_EC50_AR:0.10:INH EC50 AR
INH_Emax_SR:0.50:INH Emax SR
INH_EC50_SR:0.50:INH EC50 SR
INH_hill:1.0:INH Hill
PZA_Emax_SR:2.00:PZA Emax SR
PZA_EC50_SR:20.0:PZA EC50 SR
PZA_Emax_NR:1.00:PZA Emax NR
PZA_EC50_NR:30.0:PZA EC50 NR
PZA_hill:1.5:PZA Hill
EMB_Emax_AR:0.80:EMB Emax AR
EMB_EC50_AR:1.0:EMB EC50 AR
EMB_hill:1.0:EMB Hill
BDQ_Emax_SR:2.50:BDQ Emax SR
BDQ_EC50_SR:0.06:BDQ EC50 SR
BDQ_Emax_NR:2.00:BDQ Emax NR
BDQ_EC50_NR:0.10:BDQ EC50 NR
BDQ_hill:1.0:BDQ Hill
rpoB_mut:0.0:rpoB mutation
katG_mut:0.0:katG mutation
on_RIF:1.0:RIF on
on_INH:1.0:INH on
on_PZA:1.0:PZA on
on_EMB:1.0:EMB on
adhere:1.0:Adherence

$CMT RIF_gut RIF_c INH_c PZA_c EMB_c BDQ_c AR SR NR UM IM AM Th1 IFNg TNFa IL10

$GLOBAL
double RIF_conc, INH_conc, PZA_conc, EMB_conc, BDQ_conc;
double kRIF_AR,kRIF_SR,kINH_AR,kINH_SR,kPZA_SR,kPZA_NR,kEMB_AR_inh,kBDQ_SR,kBDQ_NR;

$MAIN
AR_0=10; SR_0=0; NR_0=0;
UM_0=UM_ss; IM_0=0; AM_0=0;
Th1_0=0; IFNg_0=0.001; TNFa_0=0.01; IL10_0=0.01;

$ODE
double rif_in  = on_RIF * adhere * RIF_dose_mg * RIF_F * RIF_ka;
double rif_ka  = RIF_ka * RIF_gut;
double rif_cl  = RIF_CL_base*(1+RIF_ind_max*(RIF_c/RIF_V)/(RIF_ind_ec50+RIF_c/RIF_V));
dxdt_RIF_gut   = rif_in - rif_ka;
dxdt_RIF_c     = rif_ka - (rif_cl/RIF_V)*RIF_c;
dxdt_INH_c     = on_INH*adhere*INH_dose_mg*INH_F*INH_ka - (INH_CL*NAT2_factor/INH_V)*INH_c;
dxdt_PZA_c     = on_PZA*adhere*PZA_dose_mg*PZA_F*PZA_ka - (PZA_CL/PZA_V)*PZA_c;
dxdt_EMB_c     = on_EMB*adhere*EMB_dose_mg*EMB_F*EMB_ka - (EMB_CL/EMB_V)*EMB_c;
dxdt_BDQ_c     = BDQ_on*adhere*BDQ_dose_mg*BDQ_F*BDQ_ka - (BDQ_CL/BDQ_V)*BDQ_c;

RIF_conc = (RIF_c/RIF_V)*RIF_lung_f;
INH_conc = INH_c/INH_V;
PZA_conc = (PZA_c/PZA_V)*PZA_pH_act;
EMB_conc = EMB_c/EMB_V;
BDQ_conc = BDQ_c/BDQ_V;

double H_rif=pow(RIF_conc,RIF_hill), H_inh=pow(INH_conc,INH_hill);
double H_pza=pow(PZA_conc,PZA_hill), H_emb=pow(EMB_conc,EMB_hill);
double H_bdq=pow(BDQ_conc,BDQ_hill);
kRIF_AR=(1-rpoB_mut)*RIF_Emax_AR*H_rif/(pow(RIF_EC50_AR,RIF_hill)+H_rif);
kRIF_SR=(1-rpoB_mut)*RIF_Emax_SR*H_rif/(pow(RIF_EC50_SR,RIF_hill)+H_rif);
kINH_AR=(1-katG_mut)*INH_Emax_AR*H_inh/(pow(INH_EC50_AR,INH_hill)+H_inh);
kINH_SR=(1-katG_mut)*INH_Emax_SR*H_inh/(pow(INH_EC50_SR,INH_hill)+H_inh);
kPZA_SR=PZA_Emax_SR*H_pza/(pow(PZA_EC50_SR,PZA_hill)+H_pza);
kPZA_NR=PZA_Emax_NR*H_pza/(pow(PZA_EC50_NR,PZA_hill)+H_pza);
kEMB_AR_inh=EMB_Emax_AR*H_emb/(pow(EMB_EC50_AR,EMB_hill)+H_emb);
kBDQ_SR=BDQ_on*BDQ_Emax_SR*H_bdq/(pow(BDQ_EC50_SR,BDQ_hill)+H_bdq);
kBDQ_NR=BDQ_on*BDQ_Emax_NR*H_bdq/(pow(BDQ_EC50_NR,BDQ_hill)+H_bdq);

double ARp=AR<0?0:AR, SRp=SR<0?0:SR, NRp=NR<0?0:NR;
double UMp=UM<0?0:UM, IMp=IM<0?0:IM, AMp=AM<0?0:AM;
double Tp =Th1<0?0:Th1, Gp=IFNg<0?0:IFNg;
double tot=ARp+SRp+NRp, tot_s=tot<1?1:tot;
dxdt_UM=kUM_p - kUM_d*UMp - kInf*UMp*ARp;
dxdt_IM=kInf*UMp*ARp - kIM_d*IMp;
double Gact=Gp/(kIL10_sup+(IL10<0?0:IL10));
dxdt_AM=kActv*phi_Mact*Gact*IMp - kAM_d*AMp;
double Tstim=tot_s/(kTh1_half+tot_s);
dxdt_Th1=phi_Th1*kTh1_p*Tstim - kTh1_d*Tp;
dxdt_IFNg=kIFNg_p*Tp - kIFNg_d*Gp;
dxdt_TNFa=kTNFa_p*IMp - kTNFa_d*(TNFa<0?0:TNFa);
dxdt_IL10=kIL10_p*Tp - kIL10_d*(IL10<0?0:IL10);
double ar_grow=(kgAR*ARp*(1-tot_s/Bmax) - kEMB_AR_inh*ARp) - (kRIF_AR+kINH_AR+kdAR)*ARp - kAMkill*AMp*ARp/(100+ARp) - kAR2SR*ARp;
dxdt_AR=ar_grow;
double sr_grow=kgSR*SRp*(1-tot_s/Bmax) - (kRIF_SR+kINH_SR+kPZA_SR+kBDQ_SR+kdSR)*SRp + kAR2SR*ARp - kSR2NR*SRp;
dxdt_SR=sr_grow;
dxdt_NR=kSR2NR*SRp - (kPZA_NR+kBDQ_NR+kdNR)*NRp - kNR2AR*NRp;

$TABLE
double tot_b=AR+SR+NR;
double log10_AR=log10(AR+0.001), log10_SR=log10(SR+0.001);
double log10_NR=log10(NR+0.001), log10_total=log10(tot_b+0.001);
double RIF_Cp=RIF_c/RIF_V, INH_Cp=INH_c/INH_V;
double PZA_Cp=PZA_c/PZA_V, EMB_Cp=EMB_c/EMB_V, BDQ_Cp=BDQ_c/BDQ_V;
double cult_conv=(tot_b<10)?1.0:0.0;
double smear_pos=(AR>1e4)?1.0:0.0;

$CAPTURE log10_AR log10_SR log10_NR log10_total
RIF_Cp INH_Cp PZA_Cp EMB_Cp BDQ_Cp
UM IM AM Th1 IFNg TNFa IL10
cult_conv smear_pos
'
mod <- mcode("tb_shiny", tb_code, quiet = TRUE)

# Helper: build dosing events for standard RIPE
make_ev <- function(rif=600, inh=300, pza=1500, emb=1200,
                    d_rife=60, d_ri=180, adhere=1.0) {
  set.seed(99)
  dr <- seq(0, d_ri-1); drife <- seq(0, d_rife-1)
  miss <- function(n) rbinom(n, 1, adhere)
  rbind(
    data.frame(time=dr*24,    amt=miss(length(dr))*rif,  cmt=1, evid=1),
    data.frame(time=dr*24,    amt=miss(length(dr))*inh,  cmt=3, evid=1),
    data.frame(time=drife*24, amt=miss(length(drife))*pza, cmt=4, evid=1),
    data.frame(time=drife*24, amt=miss(length(drife))*emb, cmt=5, evid=1)
  ) %>% arrange(time)
}

run_sim <- function(par_list, ev_df = NULL, end_day = 210, delta = 24) {
  m <- do.call(param, c(list(mod), par_list))
  if (!is.null(ev_df) && nrow(ev_df) > 0)
    m <- data_set(m, ev_df)
  out <- mrgsim(m, end = end_day * 24, delta = delta, hmax = 6)
  as_tibble(out) %>% mutate(day = time / 24)
}

# ─── UI ─────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "green",
  dashboardHeader(
    title = "TB QSP Model",
    titleWidth = 280
  ),
  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("① Patient Profile",     tabName = "profile",  icon = icon("user")),
      menuItem("② Drug PK (RIPE)",      tabName = "pk",       icon = icon("capsules")),
      menuItem("③ Bacterial Dynamics",  tabName = "bact",     icon = icon("bacterium")),
      menuItem("④ Immune Response",     tabName = "immune",   icon = icon("shield-halved")),
      menuItem("⑤ Clinical Endpoints",  tabName = "endpoint", icon = icon("stethoscope")),
      menuItem("⑥ Scenario Comparison", tabName = "compare",  icon = icon("chart-bar")),
      menuItem("⑦ PD / Biomarkers",     tabName = "pd",       icon = icon("vials"))
    ),
    hr(),
    h5("  Disease Parameters", style = "color:#a5d6a7; margin-left:15px"),
    sliderInput("phi_Th1",  "Th1 Scaling (HIV: 0.2)",  0.1, 1.0, 1.0, 0.05),
    sliderInput("phi_Mact", "Mφ Activation (DM: 0.4)", 0.1, 1.0, 1.0, 0.05),
    sliderInput("adhere",   "Adherence",                0.5, 1.0, 1.0, 0.05),
    hr(),
    h5("  Drug Selection", style = "color:#a5d6a7; margin-left:15px"),
    checkboxInput("on_RIF", "Rifampicin (RIF)",   TRUE),
    checkboxInput("on_INH", "Isoniazid (INH)",    TRUE),
    checkboxInput("on_PZA", "Pyrazinamide (PZA)", TRUE),
    checkboxInput("on_EMB", "Ethambutol (EMB)",   TRUE),
    checkboxInput("on_BDQ", "Bedaquiline (BDQ, DR-TB)", FALSE),
    hr(),
    h5("  Resistance", style = "color:#a5d6a7; margin-left:15px"),
    checkboxInput("rpoB", "rpoB mutation (RIF-R)", FALSE),
    checkboxInput("katG", "katG/inhA mutation (INH-R)", FALSE),
    hr(),
    numericInput("sim_days", "Simulation Duration (days)", 210, 30, 365, 30),
    actionButton("run_btn", "▶ Run Simulation",
                 class = "btn-success", style = "margin:10px; width:240px")
  ),

  dashboardBody(
    tags$head(tags$style(HTML(
      ".box-title{font-weight:bold}
       .content-wrapper{background:#f5f5f5}
       .main-header .logo{font-weight:bold}"
    ))),
    tabItems(

      # ── Tab 1: Patient Profile ────────────────────────────────────────────
      tabItem("profile",
        fluidRow(
          box(width=4, title="🦠 Disease Overview", status="success", solidHeader=TRUE,
            tags$div(style="font-size:13px",
              tags$h4("Tuberculosis (결핵)"),
              tags$p("Caused by ", tags$em("Mycobacterium tuberculosis"), " (Mtb),
                     an aerobic, slow-growing, acid-fast bacillus."),
              tags$p("TB remains one of the world's deadliest infectious diseases,
                     causing ~1.6 million deaths/year (WHO 2023)."),
              tags$hr(),
              tags$h5("Epidemiology"),
              tags$ul(
                tags$li("Global incidence: ~10.6 million new cases (2022)"),
                tags$li("~1/3 of world population latently infected (LTBI)"),
                tags$li("High burden: South Asia, sub-Saharan Africa"),
                tags$li("MDR-TB: ~450,000 cases/year")
              )
            )
          ),
          box(width=4, title="🔬 Pathophysiology", status="info", solidHeader=TRUE,
            tags$div(style="font-size:13px",
              tags$h5("Infection → Granuloma → Disease"),
              tags$ol(
                tags$li("Aerosol inhalation of 1–5 Mtb bacilli"),
                tags$li("Alveolar macrophage phagocytosis; phagosome escape"),
                tags$li("Innate response: PMN, NK, early IFN-γ"),
                tags$li("Adaptive: Th1 cells produce IFN-γ → M1 activation"),
                tags$li("Granuloma formation: organized vs. caseous"),
                tags$li("LTBI: NR persister state; reactivation risk ~5–10%")
              ),
              tags$h5("Bacterial States"),
              tags$ul(
                tags$li(tags$b("AR (Active Replicating):"), " killed by INH, RIF"),
                tags$li(tags$b("SR (Slow Replicating):"), " killed by RIF, PZA, BDQ"),
                tags$li(tags$b("NR (Non-Replicating):"), " killed by PZA, BDQ")
              )
            )
          ),
          box(width=4, title="💊 Treatment Regimens", status="warning", solidHeader=TRUE,
            tags$div(style="font-size:13px",
              tags$h5("DS-TB: Standard RIPE"),
              tags$table(class="table table-condensed",
                tags$tr(tags$th("Phase"), tags$th("Drugs"), tags$th("Duration")),
                tags$tr(tags$td("Intensive"), tags$td("HRZE"), tags$td("2 months")),
                tags$tr(tags$td("Continuation"), tags$td("HR"), tags$td("4 months"))
              ),
              tags$h5("MDR-TB: BDQ-Based (WHO 2022)"),
              tags$table(class="table table-condensed",
                tags$tr(tags$th("Drug"), tags$th("Target")),
                tags$tr(tags$td("Bedaquiline"), tags$td("ATP synthase")),
                tags$tr(tags$td("Pretomanid"), tags$td("F420")),
                tags$tr(tags$td("Linezolid"), tags$td("50S ribosome")),
                tags$tr(tags$td("PZA"), tags$td("Sterilizing"))
              ),
              tags$p(tags$b("Treatment success rate: "), "DS-TB: 85%; MDR-TB: 57%")
            )
          )
        ),
        fluidRow(
          box(width=12, title="📷 Mechanistic Map", status="primary", solidHeader=TRUE,
            tags$img(src="tb_qsp_model.png",
                     style="max-width:100%; border:1px solid #ccc; border-radius:4px;",
                     onerror="this.style.display='none'"),
            tags$p(style="color:#666; font-size:12px",
                   "Full mechanistic map. View tb_qsp_model.svg for interactive version.")
          )
        )
      ),

      # ── Tab 2: Drug PK ────────────────────────────────────────────────────
      tabItem("pk",
        fluidRow(
          box(width=8, title="💊 RIPE Drug PK Profiles", status="primary", solidHeader=TRUE,
            plotOutput("pk_plot", height="450px")
          ),
          box(width=4, title="PK Parameters", status="info", solidHeader=TRUE,
            sliderInput("rif_dose", "RIF dose (mg/day)", 300, 900, 600, 50),
            sliderInput("inh_dose", "INH dose (mg/day)", 100, 500, 300, 50),
            sliderInput("pza_dose", "PZA dose (mg/day)", 500, 2500, 1500, 100),
            sliderInput("emb_dose", "EMB dose (mg/day)", 400, 2000, 1200, 100),
            sliderInput("nat2",     "NAT2 Factor (Slow=0.4, Fast=1.8)", 0.3, 2.0, 1.0, 0.1),
            sliderInput("pk_days",  "Display days", 3, 30, 14, 1),
            tableOutput("pk_targets")
          )
        ),
        fluidRow(
          box(width=12, title="📊 Steady-State PK Summary", status="success", solidHeader=TRUE,
            DTOutput("pk_ss_table")
          )
        )
      ),

      # ── Tab 3: Bacterial Dynamics ─────────────────────────────────────────
      tabItem("bact",
        fluidRow(
          box(width=8, title="🦠 Bacterial Subpopulation Dynamics", status="danger", solidHeader=TRUE,
            plotOutput("bact_plot", height="400px")
          ),
          box(width=4, title="Bacterial Parameters", status="warning", solidHeader=TRUE,
            sliderInput("kgAR",   "AR Growth Rate (/day)", 0.1, 1.4, 0.6931, 0.05),
            sliderInput("kgSR",   "SR Growth Rate (/day)", 0.01, 0.20, 0.069, 0.01),
            sliderInput("kAR2SR", "AR→SR Transition",      0.01, 0.15, 0.05,  0.01),
            sliderInput("kSR2NR", "SR→NR Transition",      0.005, 0.05, 0.02, 0.005),
            sliderInput("Bmax_log10", "Max CFU (log10)",    7, 11, 9, 0.5)
          )
        ),
        fluidRow(
          box(width=6, title="Total Bacterial Burden", status="primary", solidHeader=TRUE,
            plotOutput("bact_total", height="300px")
          ),
          box(width=6, title="Bacterial Killing Rates by Drug", status="success", solidHeader=TRUE,
            plotOutput("kill_rate", height="300px")
          )
        )
      ),

      # ── Tab 4: Immune Response ────────────────────────────────────────────
      tabItem("immune",
        fluidRow(
          box(width=8, title="🛡️ Macrophage & T-Cell Dynamics", status="info", solidHeader=TRUE,
            plotOutput("immune_plot", height="450px")
          ),
          box(width=4, title="Immune Parameters", status="warning", solidHeader=TRUE,
            sliderInput("phi_Th1_tab", "Th1 Scaling",          0.1, 1.0, 1.0, 0.05),
            sliderInput("phi_Mact_tab","Mφ Activation Scaling", 0.1, 1.0, 1.0, 0.05),
            sliderInput("kActv",       "Mφ Activation Rate",    0.1, 2.0, 0.5, 0.1),
            sliderInput("kAMkill",     "AM Kill Rate",          0.5, 5.0, 2.0, 0.25),
            hr(),
            tags$div(style="font-size:12px",
              tags$h5("Clinical Context"),
              tags$ul(
                tags$li("HIV: phi_Th1 = 0.1–0.3 (CD4 <200)"),
                tags$li("Diabetes: phi_Mact = 0.3–0.5"),
                tags$li("Anti-TNF: reduces granuloma stability")
              )
            )
          )
        ),
        fluidRow(
          box(width=6, title="Cytokine Dynamics", status="success", solidHeader=TRUE,
            plotOutput("cytokine_plot", height="300px")
          ),
          box(width=6, title="Immune Cell Ratios", status="primary", solidHeader=TRUE,
            plotOutput("immune_ratio", height="300px")
          )
        )
      ),

      # ── Tab 5: Clinical Endpoints ─────────────────────────────────────────
      tabItem("endpoint",
        fluidRow(
          box(width=6, title="📊 Sputum Culture Conversion", status="success", solidHeader=TRUE,
            plotOutput("cult_plot", height="350px"),
            tags$p(style="color:#333; font-size:12px",
                   "Culture conversion = bacterial load < 10 CFU/mL.",
                   "WHO target: ≥80% of cohort converted by 2 months.")
          ),
          box(width=6, title="📊 Sputum Smear Status", status="primary", solidHeader=TRUE,
            plotOutput("smear_plot", height="350px"),
            tags$p(style="color:#333; font-size:12px",
                   "Smear positive when AR bacteria > 10⁴ CFU/mL.",
                   "Smear generally converts before culture.")
          )
        ),
        fluidRow(
          box(width=12, title="📋 Clinical Outcome Summary", status="info", solidHeader=TRUE,
            DTOutput("outcome_table")
          )
        )
      ),

      # ── Tab 6: Scenario Comparison ────────────────────────────────────────
      tabItem("compare",
        fluidRow(
          box(width=12, title="🔄 6-Scenario Comparison", status="success", solidHeader=TRUE,
            plotOutput("compare_plot", height="550px")
          )
        ),
        fluidRow(
          box(width=6, title="Bacterial Kinetics by Scenario", status="primary", solidHeader=TRUE,
            plotOutput("compare_bact", height="350px")
          ),
          box(width=6, title="Immune Response by Scenario", status="warning", solidHeader=TRUE,
            plotOutput("compare_immune", height="350px")
          )
        )
      ),

      # ── Tab 7: PD / Biomarkers ────────────────────────────────────────────
      tabItem("pd",
        fluidRow(
          box(width=6, title="⚗️ Drug Concentration-Effect (Emax)", status="info", solidHeader=TRUE,
            plotOutput("pd_curve", height="400px"),
            selectInput("pd_drug", "Drug:", choices=c("RIF","INH","PZA","EMB","BDQ"))
          ),
          box(width=6, title="🎯 PKPD Target Attainment", status="warning", solidHeader=TRUE,
            plotOutput("target_plot", height="400px"),
            sliderInput("pd_day", "Simulation Day:", 1, 180, 30, 5)
          )
        ),
        fluidRow(
          box(width=12, title="📊 Biomarker Summary at Selected Day", status="success", solidHeader=TRUE,
            DTOutput("biomarker_table")
          )
        )
      )
    )
  )
)

# ─── SERVER ──────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Reactive: run simulation on button click
  sim_data <- eventReactive(input$run_btn, {
    par_list <- list(
      phi_Th1 = input$phi_Th1,
      phi_Mact = input$phi_Mact,
      adhere   = input$adhere,
      on_RIF = as.numeric(input$on_RIF),
      on_INH = as.numeric(input$on_INH),
      on_PZA = as.numeric(input$on_PZA),
      on_EMB = as.numeric(input$on_EMB),
      BDQ_on = as.numeric(input$on_BDQ),
      rpoB_mut = as.numeric(input$rpoB),
      katG_mut = as.numeric(input$katG)
    )
    ev_df <- make_ev(adhere = input$adhere,
                     d_ri   = min(input$sim_days, 180),
                     d_rife = min(input$sim_days, 60))
    run_sim(par_list, ev_df, end_day = input$sim_days)
  }, ignoreNULL = FALSE)

  # Scenario comparison data
  sc_data <- reactive({
    sc <- list(
      list(label="1. Natural History",
           par=list(on_RIF=0,on_INH=0,on_PZA=0,on_EMB=0,BDQ_on=0,phi_Th1=1,phi_Mact=1,adhere=1,rpoB_mut=0,katG_mut=0)),
      list(label="2. Standard RIPE",
           par=list(on_RIF=1,on_INH=1,on_PZA=1,on_EMB=1,BDQ_on=0,phi_Th1=1,phi_Mact=1,adhere=1,rpoB_mut=0,katG_mut=0)),
      list(label="3. Poor Adherence",
           par=list(on_RIF=1,on_INH=1,on_PZA=1,on_EMB=1,BDQ_on=0,phi_Th1=1,phi_Mact=1,adhere=0.7,rpoB_mut=0,katG_mut=0)),
      list(label="4. MDR-TB (BDQ)",
           par=list(on_RIF=0,on_INH=0,on_PZA=1,on_EMB=0,BDQ_on=1,phi_Th1=1,phi_Mact=1,adhere=1,rpoB_mut=1,katG_mut=1)),
      list(label="5. HIV (CD4<200)",
           par=list(on_RIF=1,on_INH=1,on_PZA=1,on_EMB=1,BDQ_on=0,phi_Th1=0.2,phi_Mact=0.6,adhere=0.9,rpoB_mut=0,katG_mut=0)),
      list(label="6. Diabetic Host",
           par=list(on_RIF=1,on_INH=1,on_PZA=1,on_EMB=1,BDQ_on=0,phi_Th1=0.8,phi_Mact=0.4,adhere=1,rpoB_mut=0,katG_mut=0))
    )
    lapply(sc, function(s) {
      ev_df <- if (s$par$on_RIF == 1 || s$par$BDQ_on == 1)
        make_ev(adhere = s$par$adhere) else NULL
      run_sim(s$par, ev_df, end_day = 210) %>% mutate(scenario = s$label)
    }) %>% bind_rows()
  })

  sc_cols <- c(
    "1. Natural History" = "#B71C1C",
    "2. Standard RIPE"   = "#1B5E20",
    "3. Poor Adherence"  = "#F57F17",
    "4. MDR-TB (BDQ)"    = "#4A148C",
    "5. HIV (CD4<200)"   = "#880E4F",
    "6. Diabetic Host"   = "#1565C0"
  )

  theme_tb <- theme_bw(base_size=12) +
    theme(legend.position="bottom", panel.grid.minor=element_blank())

  # ── Tab 2: PK plots --------------------------------------------------------
  output$pk_plot <- renderPlot({
    d <- sim_data() %>%
      filter(day <= input$pk_days) %>%
      select(day, RIF_Cp, INH_Cp, PZA_Cp, EMB_Cp) %>%
      pivot_longer(-day, names_to="drug", values_to="conc") %>%
      mutate(drug = recode(drug,
        RIF_Cp="Rifampicin", INH_Cp="Isoniazid",
        PZA_Cp="Pyrazinamide", EMB_Cp="Ethambutol"
      ))
    ggplot(d, aes(day, conc, color=drug)) +
      geom_line(linewidth=1.2) +
      facet_wrap(~drug, scales="free_y", nrow=2) +
      scale_color_manual(values=c("#D32F2F","#1565C0","#2E7D32","#6D4C41")) +
      labs(title="Drug Plasma Concentration — Standard RIPE",
           x="Day", y="Concentration (mg/L)") +
      theme_tb + theme(legend.position="none")
  })

  output$pk_targets <- renderTable({
    data.frame(
      Drug = c("Rifampicin","Isoniazid","Pyrazinamide","Ethambutol"),
      `PK Target` = c("Cmax >8 mg/L, AUC/MIC >271",
                      "Cmax >3 mg/L, AUC/MIC >52",
                      "Cmax >20 mg/L, AUC/MIC >209",
                      "Cmax >2 mg/L"),
      MIC_Susceptible = c("≤0.5 mg/L","≤0.2 mg/L","≤50 mg/L","≤2 mg/L")
    )
  }, bordered=TRUE, striped=TRUE, hover=TRUE)

  output$pk_ss_table <- renderDT({
    d <- sim_data() %>% filter(day >= 7, day <= 14)
    data.frame(
      Drug = c("Rifampicin","Isoniazid","Pyrazinamide","Ethambutol"),
      Cmax_mg_L = round(c(max(d$RIF_Cp,na.rm=T), max(d$INH_Cp,na.rm=T),
                           max(d$PZA_Cp,na.rm=T), max(d$EMB_Cp,na.rm=T)), 2),
      Cmin_mg_L = round(c(min(d$RIF_Cp,na.rm=T), min(d$INH_Cp,na.rm=T),
                           min(d$PZA_Cp,na.rm=T), min(d$EMB_Cp,na.rm=T)), 2)
    )
  }, options=list(dom='t'), rownames=FALSE)

  # ── Tab 3: Bacterial Dynamics ──────────────────────────────────────────────
  output$bact_plot <- renderPlot({
    d <- sim_data() %>%
      select(day, log10_AR, log10_SR, log10_NR) %>%
      pivot_longer(-day, names_to="pop", values_to="log10_cfu") %>%
      mutate(pop = recode(pop,
        log10_AR="Active Replicating (AR)",
        log10_SR="Slow Replicating (SR)",
        log10_NR="Non-Replicating (NR)"))
    ggplot(d, aes(day, log10_cfu, color=pop)) +
      geom_line(linewidth=1.2) +
      geom_hline(yintercept=1, linetype="dashed", color="grey50") +
      annotate("text", x=3, y=1.4, label="Culture threshold", size=3.5) +
      scale_color_manual(values=c("#E65100","#1565C0","#6A1B9A")) +
      labs(title="Bacterial Subpopulation Dynamics",
           x="Day", y="log₁₀ CFU/mL", color="Population") +
      theme_tb
  })

  output$bact_total <- renderPlot({
    sim_data() %>%
      ggplot(aes(day, log10_total)) +
      geom_line(linewidth=1.2, color="#B71C1C") +
      geom_hline(yintercept=1, linetype="dashed", color="#1B5E20") +
      labs(title="Total Bacterial Burden", x="Day", y="log₁₀ Total CFU/mL") +
      theme_tb
  })

  output$kill_rate <- renderPlot({
    days <- seq(1, min(input$sim_days, 180))
    conc_rif <- sim_data()$RIF_Cp[match(days, round(sim_data()$day))]
    conc_rif[is.na(conc_rif)] <- 0
    rif_kill <- 3.5 * conc_rif^1 / (0.60^1 + conc_rif^1) * 0.25
    data.frame(day=days, RIF_kill=rif_kill) %>%
      ggplot(aes(day, RIF_kill)) +
      geom_line(linewidth=1.1, color="#D32F2F") +
      labs(title="RIF Effective Kill Rate on AR",
           x="Day", y="Kill Rate (/day)") +
      theme_tb
  })

  # ── Tab 4: Immune Response ─────────────────────────────────────────────────
  output$immune_plot <- renderPlot({
    d <- sim_data() %>%
      select(day, UM, IM, AM, Th1) %>%
      pivot_longer(-day, names_to="var", values_to="val") %>%
      mutate(var = recode(var,
        UM="Uninfected Mφ", IM="Infected Mφ", AM="Activated Mφ", Th1="Th1 Cells"))
    ggplot(d, aes(day, val, color=var)) +
      geom_line(linewidth=1.2) +
      facet_wrap(~var, scales="free_y", nrow=2) +
      scale_color_manual(values=c("#4CAF50","#F44336","#FF9800","#9C27B0")) +
      labs(title="Macrophage & Th1 Dynamics", x="Day", y="Cells (AU)") +
      theme_tb + theme(legend.position="none")
  })

  output$cytokine_plot <- renderPlot({
    d <- sim_data() %>%
      select(day, IFNg, TNFa, IL10) %>%
      pivot_longer(-day, names_to="cyt", values_to="val") %>%
      mutate(cyt = recode(cyt, IFNg="IFN-γ", TNFa="TNF-α", IL10="IL-10"))
    ggplot(d, aes(day, val, color=cyt)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=c("#2196F3","#FF5722","#9C27B0")) +
      facet_wrap(~cyt, scales="free_y") +
      labs(title="Key Cytokines", x="Day", y="Concentration (AU)") +
      theme_tb + theme(legend.position="none")
  })

  output$immune_ratio <- renderPlot({
    sim_data() %>%
      mutate(AM_UM_ratio = AM / (UM + 0.1),
             Th1_IM_ratio = Th1 / (IM + 0.1)) %>%
      select(day, AM_UM_ratio, Th1_IM_ratio) %>%
      pivot_longer(-day, names_to="ratio", values_to="val") %>%
      ggplot(aes(day, val, color=ratio)) +
      geom_line(linewidth=1.1) +
      scale_color_manual(values=c("#FF9800","#9C27B0")) +
      labs(title="Immune Cell Ratios", x="Day", y="Ratio",
           color="Ratio") +
      theme_tb
  })

  # ── Tab 5: Clinical Endpoints ──────────────────────────────────────────────
  output$cult_plot <- renderPlot({
    sim_data() %>%
      ggplot(aes(day, cult_conv)) +
      geom_step(linewidth=1.3, color="#1B5E20") +
      geom_hline(yintercept=0.5, linetype="dashed", color="grey50") +
      scale_y_continuous(limits=c(-0.05, 1.1)) +
      labs(title="Sputum Culture Conversion",
           x="Day", y="Culture Converted (1=Yes, 0=No)") +
      theme_tb
  })

  output$smear_plot <- renderPlot({
    sim_data() %>%
      ggplot(aes(day, smear_pos)) +
      geom_step(linewidth=1.3, color="#D32F2F") +
      scale_y_continuous(limits=c(-0.05, 1.1)) +
      labs(title="Sputum Smear Positivity",
           x="Day", y="Smear Positive (1=Yes, 0=No)") +
      theme_tb
  })

  output$outcome_table <- renderDT({
    d <- sim_data()
    day_conv <- suppressWarnings(min(d$day[d$cult_conv == 1], na.rm=TRUE))
    peak_bact <- max(d$log10_total, na.rm=TRUE)
    final_bact <- tail(d$log10_total, 1)
    smear_dur <- sum(d$smear_pos == 1, na.rm=TRUE)
    data.frame(
      Metric = c("Peak Bacterial Burden (log10)",
                 "Day of Culture Conversion",
                 "Final Bacterial Burden (log10)",
                 "Smear Positive Duration (days)",
                 "Treatment Success (binary)"),
      Value = c(
        round(peak_bact, 2),
        ifelse(is.infinite(day_conv), ">210 (failure)", round(day_conv, 0)),
        round(final_bact, 2),
        smear_dur,
        ifelse(final_bact < 1, "Yes (cured)", "No (active)")
      )
    )
  }, options=list(dom='t'), rownames=FALSE)

  # ── Tab 6: Scenario Comparison ─────────────────────────────────────────────
  output$compare_plot <- renderPlot({
    sc_data() %>%
      ggplot(aes(day, log10_total, color=scenario)) +
      geom_line(linewidth=1.1) +
      geom_hline(yintercept=1, linetype="dashed", color="grey40") +
      annotate("text", x=10, y=1.4, label="Culture\nthreshold", size=3) +
      scale_color_manual(values=sc_cols) +
      labs(title="Tuberculosis QSP — Total Bacterial Burden by Scenario",
           x="Day", y="log₁₀ Total Bacterial Burden (CFU/mL)", color="Scenario") +
      theme_tb + guides(color=guide_legend(ncol=2))
  })

  output$compare_bact <- renderPlot({
    sc_data() %>%
      select(day, scenario, log10_AR, log10_SR) %>%
      pivot_longer(c(log10_AR, log10_SR), names_to="pop", values_to="val") %>%
      mutate(pop=recode(pop, log10_AR="AR", log10_SR="SR")) %>%
      ggplot(aes(day, val, color=scenario, linetype=pop)) +
      geom_line() +
      scale_color_manual(values=sc_cols) +
      labs(title="AR vs SR Dynamics", x="Day", y="log₁₀ CFU/mL") +
      theme_tb + guides(color=guide_legend(ncol=1))
  })

  output$compare_immune <- renderPlot({
    sc_data() %>%
      ggplot(aes(day, Th1, color=scenario)) +
      geom_line(linewidth=1.0) +
      scale_color_manual(values=sc_cols) +
      labs(title="Th1 Cell Response by Scenario",
           x="Day", y="Th1 (AU)") +
      theme_tb + guides(color=guide_legend(ncol=1))
  })

  # ── Tab 7: PD / Biomarkers ─────────────────────────────────────────────────
  output$pd_curve <- renderPlot({
    conc_range <- seq(0, 10, length.out=200)
    drug <- input$pd_drug
    params <- switch(drug,
      RIF = list(Emax=3.5, EC50=0.60, hill=1.0, label="Rifampicin (AR)"),
      INH = list(Emax=4.0, EC50=0.10, hill=1.0, label="Isoniazid (AR)"),
      PZA = list(Emax=2.0, EC50=20.0, hill=1.5, label="Pyrazinamide (SR)"),
      EMB = list(Emax=0.8, EC50=1.00, hill=1.0, label="Ethambutol (AR)"),
      BDQ = list(Emax=2.5, EC50=0.06, hill=1.0, label="Bedaquiline (SR)")
    )
    conc_rng <- if (drug %in% c("BDQ")) seq(0,0.5,length.out=200) else
                if (drug == "PZA") seq(0,60,length.out=200) else conc_range
    kill <- params$Emax * conc_rng^params$hill /
            (params$EC50^params$hill + conc_rng^params$hill)
    ggplot(data.frame(conc=conc_rng, kill=kill), aes(conc, kill)) +
      geom_line(linewidth=1.3, color="#1B5E20") +
      geom_vline(xintercept=params$EC50, linetype="dashed", color="#D32F2F") +
      annotate("text", x=params$EC50*1.1, y=params$Emax*0.1,
               label=paste0("EC50=",params$EC50,"mg/L"), color="#D32F2F", size=3.5) +
      labs(title=paste("Emax PD Curve —", params$label),
           x="Drug Concentration (mg/L)", y="Kill Rate (/day)") +
      theme_tb
  })

  output$target_plot <- renderPlot({
    d <- sim_data() %>% filter(day == input$pd_day)
    if (nrow(d) == 0) d <- sim_data() %>% filter(day <= input$pd_day) %>% tail(1)
    drugs <- c("Rifampicin","Isoniazid","Pyrazinamide","Ethambutol")
    conc_vals <- c(d$RIF_Cp, d$INH_Cp, d$PZA_Cp, d$EMB_Cp)
    targets   <- c(2, 3, 20, 2)
    hit       <- conc_vals >= targets
    data.frame(Drug=drugs, Concentration=conc_vals, Target=targets, Hit=hit) %>%
      ggplot(aes(x=Drug, y=Concentration, fill=Hit)) +
      geom_col(alpha=0.8) +
      geom_point(aes(y=Target), shape=18, size=5, color="#D32F2F") +
      scale_fill_manual(values=c("TRUE"="#1B5E20","FALSE"="#B71C1C"),
                        labels=c("TRUE"="Target Achieved","FALSE"="Below Target")) +
      labs(title=paste("PKPD Target Attainment at Day", input$pd_day),
           x="Drug", y="Concentration (mg/L)", fill="") +
      theme_tb
  })

  output$biomarker_table <- renderDT({
    d <- sim_data() %>% filter(day <= input$pd_day) %>% tail(1)
    data.frame(
      Biomarker = c("Total Bacteria (log10 CFU/mL)",
                    "Active Replicating (log10)",
                    "Slow Replicating (log10)",
                    "Non-Replicating (log10)",
                    "Activated Macrophages",
                    "Th1 Cells",
                    "IFN-γ",
                    "TNF-α",
                    "IL-10",
                    "Culture Converted",
                    "Smear Positive"),
      Value = round(c(d$log10_total, d$log10_AR, d$log10_SR, d$log10_NR,
                      d$AM, d$Th1, d$IFNg, d$TNFa, d$IL10,
                      d$cult_conv, d$smear_pos), 3)
    )
  }, options=list(dom='t'), rownames=FALSE)
}

shinyApp(ui = ui, server = server)
