## =============================================================================
## PNH (Paroxysmal Nocturnal Hemoglobinuria) — Interactive Shiny Dashboard
## QSP Model: GPI Deficiency · Complement · Hemolysis · Drug PK/PD
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(shinydashboard)

## ─────────────────────────────────────────────────────────────────
## EMBED MODEL CODE (same as mrgsolve model file, condensed)
## ─────────────────────────────────────────────────────────────────
pnh_model_code <- '
$PARAM
f_PNH=0.70, kprod_RBC=0.0083, kmat_Ret=0.33, kd_NL_RBC=0.0083,
kd_PNH_base=0.071, EPO_stim=2.0, Hgb_normal=14.0, Hgb_nadir=7.0,
C3_ss=75.0, C5_ss=75.0, FB_ss=200.0, FD_ss=2.0,
k_tickover=0.005, k_amplif=0.15, k_C3b_dep=0.003,
k_C5conv=0.12, k_MAC_lysis=0.25, k_EVH=0.15,
kd_C3b=0.5, kd_MAC=1.0, kd_C3=0.03, kd_C5=0.03,
Hgb_per_RBC=0.028, kd_fHgb=0.7, kd_fHgb_hi=2.0,
Hp_init=120.0, kd_Hp=0.1, Hp_syn=12.0,
NO_ss=1.0, k_NO_prod=0.5, k_NO_scav=0.8, k_NO_deg=0.5,
LDH_ss=200.0, k_LDH_rel=400.0, k_LDH_clear=0.35,
ECU_V1=5.5, ECU_V2=7.3, ECU_CL=0.31, ECU_Q=0.96,
ECU_kon=0.45, ECU_koff=0.0015, ECU_kdeg=0.06,
RAV_V1=4.08, RAV_V2=3.51, RAV_CL=0.069, RAV_Q=0.28,
RAV_kon=0.50, RAV_koff=0.0012,
IPC_F=0.69, IPC_ka=2.8, IPC_V=85.0, IPC_CL=4.2, IPC_IC50=0.05, IPC_Hill=1.5,
DAN_F=0.80, DAN_ka=3.5, DAN_V=55.0, DAN_CL=8.0, DAN_IC50=0.02,
use_ECU=0, use_RAV=0, use_IPC=0, use_DAN=0

$CMT
PNH_Ret NL_Ret PNH_RBC NL_RBC
C3 C3b C5 MAC
fHgb Haptoglobin LDH NO_rel
ECU_C ECU_P C5_ECU
RAV_C RAV_P C5_RAV
IPC_gut IPC_plasma
DAN_gut DAN_plasma

$INIT
PNH_Ret=0.35, NL_Ret=0.125,
PNH_RBC=3.5, NL_RBC=1.5,
C3=75.0, C3b=0.0, C5=75.0, MAC=0.0,
fHgb=0.05, Haptoglobin=120.0, LDH=850.0, NO_rel=0.6,
ECU_C=0.0, ECU_P=0.0, C5_ECU=0.0,
RAV_C=0.0, RAV_P=0.0, C5_RAV=0.0,
IPC_gut=0.0, IPC_plasma=0.0,
DAN_gut=0.0, DAN_plasma=0.0

$ODE
double C5_total  = C5 + C5_ECU + C5_RAV;
double f_C5_free = C5 / (C5_total + 1e-6);
double E_IPC = 0.0;
if (use_IPC > 0.5) {
  E_IPC = pow(IPC_plasma, IPC_Hill)/(pow(IPC_IC50,IPC_Hill)+pow(IPC_plasma,IPC_Hill));
}
double E_DAN = 0.0;
if (use_DAN > 0.5) E_DAN = DAN_plasma/(DAN_IC50+DAN_plasma);
double f_AP_block = 1.0-(1.0-E_IPC)*(1.0-E_DAN);

double total_RBC = PNH_RBC + NL_RBC;
double Hgb_curr  = total_RBC * Hgb_normal / 5.0;
double EPO_fold  = 1.0+(EPO_stim-1.0)*pow(Hgb_nadir,2)/(pow(Hgb_nadir,2)+pow(Hgb_curr,2));
double prod_PNH  = kprod_RBC * f_PNH * 5.0 * EPO_fold;
double prod_NL   = kprod_RBC * (1.0-f_PNH) * 5.0 * EPO_fold;
double rate_IVH  = k_MAC_lysis * MAC * PNH_RBC * f_C5_free;
double rate_EVH  = k_EVH * C3b * PNH_RBC * (1.0-f_AP_block);

dxdt_PNH_Ret = prod_PNH - kmat_Ret*PNH_Ret;
dxdt_NL_Ret  = prod_NL  - kmat_Ret*NL_Ret;
dxdt_PNH_RBC = kmat_Ret*PNH_Ret - kd_PNH_base*PNH_RBC - rate_IVH - rate_EVH;
dxdt_NL_RBC  = kmat_Ret*NL_Ret  - kd_NL_RBC*NL_RBC;

double rate_C3_syn = kd_C3*C3_ss;
double rate_amp    = k_amplif*C3b*C3*(1.0-f_AP_block);
dxdt_C3  = rate_C3_syn - kd_C3*C3 - k_tickover*C3 - rate_amp;
double rate_C3b_dep = k_C3b_dep*C3*PNH_RBC*(1.0-f_AP_block);
dxdt_C3b = rate_C3b_dep - kd_C3b*C3b - k_EVH*C3b*PNH_RBC*(1.0-f_AP_block);

double rate_C5_syn  = kd_C5*C5_ss;
double rate_C5_conv = k_C5conv*C3b*C5*f_C5_free;
dxdt_C5 = rate_C5_syn - kd_C5*C5 - rate_C5_conv
          - ECU_kon*ECU_C*C5 + ECU_koff*C5_ECU
          - RAV_kon*RAV_C*C5 + RAV_koff*C5_RAV;
dxdt_MAC = rate_C5_conv*PNH_RBC - kd_MAC*MAC;

double Hgb_rel = (rate_IVH*PNH_RBC + rate_EVH*0.2*PNH_RBC)*Hgb_per_RBC;
double Hp_dep  = Haptoglobin/(Haptoglobin+10.0);
double fHgb_cl = kd_fHgb*Hp_dep*fHgb + kd_fHgb_hi*(1.0-Hp_dep)*fHgb*0.3;
dxdt_fHgb       = Hgb_rel - fHgb_cl;
dxdt_Haptoglobin= Hp_syn - kd_Hp*fHgb*Haptoglobin - 0.08*Haptoglobin;
dxdt_LDH        = k_LDH_rel*(rate_IVH*PNH_RBC+rate_EVH*0.1*PNH_RBC) - k_LDH_clear*(LDH-LDH_ss);
dxdt_NO_rel     = k_NO_prod - k_NO_scav*fHgb*NO_rel - k_NO_deg*NO_rel;

double ECU_k12 = ECU_Q/ECU_V1; double ECU_k21 = ECU_Q/ECU_V2; double ECU_kel=ECU_CL/ECU_V1;
dxdt_ECU_C = -ECU_kel*ECU_C - ECU_k12*ECU_C + ECU_k21*ECU_P - ECU_kon*ECU_C*C5 + ECU_koff*C5_ECU + ECU_kdeg*C5_ECU;
dxdt_ECU_P = ECU_k12*ECU_C - ECU_k21*ECU_P;
dxdt_C5_ECU= ECU_kon*ECU_C*C5 - ECU_koff*C5_ECU - ECU_kdeg*C5_ECU;

double RAV_k12=RAV_Q/RAV_V1; double RAV_k21=RAV_Q/RAV_V2; double RAV_kel=RAV_CL/RAV_V1;
dxdt_RAV_C = -RAV_kel*RAV_C - RAV_k12*RAV_C + RAV_k21*RAV_P - RAV_kon*RAV_C*C5 + RAV_koff*C5_RAV;
dxdt_RAV_P = RAV_k12*RAV_C - RAV_k21*RAV_P;
dxdt_C5_RAV= RAV_kon*RAV_C*C5 - RAV_koff*C5_RAV - 0.05*C5_RAV;

dxdt_IPC_gut    = -IPC_ka*IPC_gut;
dxdt_IPC_plasma = IPC_F*IPC_ka*IPC_gut/IPC_V - (IPC_CL/IPC_V)*IPC_plasma;
dxdt_DAN_gut    = -DAN_ka*DAN_gut;
dxdt_DAN_plasma = DAN_F*DAN_ka*DAN_gut/DAN_V - (DAN_CL/DAN_V)*DAN_plasma;

$TABLE
double tRBC = PNH_RBC + NL_RBC;
capture Hgb        = tRBC * Hgb_normal / 5.0;
capture LDH_ULN    = LDH / 250.0;
capture fHgb_out   = fHgb;
capture NO_pct     = NO_rel * 100.0;
capture PNH_pct    = PNH_RBC/(tRBC+0.001)*100.0;
capture FACIT      = 52.0*(0.6*Hgb/(Hgb_normal)+0.4*NO_rel);
capture Thrombo    = 1.0-NO_rel+0.3*(1.0-Hgb/Hgb_normal);
capture ECU_ug     = ECU_C*1000.0;
capture RAV_ug     = RAV_C*1000.0;
capture IPC_ug     = IPC_plasma;
capture DAN_ug     = DAN_plasma;
capture C5_free_out= C5;
capture C3b_out    = C3b;
capture MAC_out    = MAC;
capture Haptoglobin_out = Haptoglobin;
'

## ─────────────────────────────────────────────────────────────────
## HELPER FUNCTIONS
## ─────────────────────────────────────────────────────────────────
build_events <- function(use_ecu, use_rav, use_ipc, use_dan) {
  ev_list <- list()
  if (use_ecu) {
    ev_list[["ECU"]] <- c(
      ev(cmt="ECU_C", time=c(0,7,14,21), amt=600/5.5),
      ev(cmt="ECU_C", time=seq(28,364,by=14), amt=900/5.5)
    )
  }
  if (use_rav) {
    ev_list[["RAV"]] <- c(
      ev(cmt="RAV_C", time=0, amt=3000/4.08),
      ev(cmt="RAV_C", time=c(14,seq(70,365,by=56)), amt=3300/4.08)
    )
  }
  if (use_ipc) {
    ev_list[["IPC"]] <- ev(cmt="IPC_gut", time=seq(0,365,by=0.5), amt=200)
  }
  if (use_dan) {
    ev_list[["DAN"]] <- ev(cmt="DAN_gut", time=seq(0,365,by=0.333), amt=150)
  }
  if (length(ev_list) == 0) return(NULL)
  Reduce(c, ev_list)
}

run_sim <- function(mod_base, f_pnh, use_ecu, use_rav, use_ipc, use_dan,
                    end_time = 365, delta = 1) {
  m <- param(mod_base,
    f_PNH = f_pnh,
    use_ECU = as.integer(use_ecu),
    use_RAV = as.integer(use_rav),
    use_IPC = as.integer(use_ipc),
    use_DAN = as.integer(use_dan)
  )
  ev_dose <- build_events(use_ecu, use_rav, use_ipc, use_dan)
  if (is.null(ev_dose)) {
    mrgsim(m, end = end_time, delta = delta) %>% as_tibble()
  } else {
    mrgsim(m, events = ev_dose, end = end_time, delta = delta) %>% as_tibble()
  }
}

## ─────────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "red",

  dashboardHeader(
    title = "PNH QSP Dashboard",
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",   tabName = "patient",  icon = icon("user")),
      menuItem("Drug PK",           tabName = "pk",       icon = icon("pills")),
      menuItem("Complement",        tabName = "complement", icon = icon("virus")),
      menuItem("Hemolysis Markers", tabName = "hemolysis",  icon = icon("droplet")),
      menuItem("Clinical Endpoints",tabName = "endpoints",  icon = icon("chart-line")),
      menuItem("Scenario Comparison",tabName="scenario",   icon = icon("code-branch")),
      menuItem("Biomarkers",        tabName = "biomarkers", icon = icon("microscope")),
      menuItem("About",             tabName = "about",     icon = icon("info-circle"))
    ),

    hr(),
    h5("  Patient Parameters", style = "color:#ddd; padding-left:15px"),
    sliderInput("f_pnh",  "PNH Clone Fraction (%)",
                min = 10, max = 95, value = 70, step = 5,
                post = "%"),
    selectInput("disease_sev", "Disease Severity",
                choices = c("Subclinical (<20% clone)" = "sub",
                            "Moderate (20–50%)"        = "mod",
                            "Severe (>50%)"            = "sev"),
                selected = "sev"),

    hr(),
    h5("  Treatment", style = "color:#ddd; padding-left:15px"),
    checkboxInput("use_ecu", "Eculizumab 900mg q2w IV",   value = FALSE),
    checkboxInput("use_rav", "Ravulizumab 3300mg q8w IV", value = FALSE),
    checkboxInput("use_ipc", "Iptacopan 200mg BID PO",    value = FALSE),
    checkboxInput("use_dan", "Danicopan 150mg TID (add-on)", value = FALSE),

    hr(),
    sliderInput("sim_end", "Simulation Duration (days)",
                min = 90, max = 730, value = 365, step = 30),
    actionButton("run_sim", "Run Simulation", icon = icon("play"),
                 class = "btn-danger btn-block")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .box-title { font-weight: bold; }
      .small-box .icon-large { font-size: 50px; }
    "))),

    tabItems(

      ## ── TAB 1: Patient Profile ─────────────────────────────────
      tabItem(tabName = "patient",
        fluidRow(
          valueBoxOutput("vb_hgb",    width = 3),
          valueBoxOutput("vb_ldh",    width = 3),
          valueBoxOutput("vb_clone",  width = 3),
          valueBoxOutput("vb_facit",  width = 3)
        ),
        fluidRow(
          box(title = "Disease Overview: PNH Pathophysiology",
              width = 7, status = "danger",
              HTML('<table class="table table-sm">
              <tr><th>Component</th><th>PNH Biology</th><th>Clinical Impact</th></tr>
              <tr><td><b>PIGA Mutation</b></td><td>Loss of GPI anchor</td><td>No CD55/CD59 on PNH cells</td></tr>
              <tr><td><b>CD55 (DAF)</b></td><td>Absent on PNH RBCs</td><td>C3 convertase unchecked</td></tr>
              <tr><td><b>CD59 (MIRL)</b></td><td>Absent on PNH RBCs</td><td>MAC pore formation</td></tr>
              <tr><td><b>Intravascular Hemolysis</b></td><td>MAC-mediated</td><td>LDH↑, Hgb↓, hemoglobinuria</td></tr>
              <tr><td><b>Extravascular Hemolysis</b></td><td>C3b opsonization</td><td>Spleen/liver phagocytosis</td></tr>
              <tr><td><b>Free Hemoglobin</b></td><td>NO scavenging</td><td>Thrombosis, smooth muscle dystonias</td></tr>
              <tr><td><b>Thrombosis</b></td><td>NO↓, platelet activation</td><td>DVT, Budd-Chiari, CVST</td></tr>
              </table>'),
              plotOutput("plot_rbc_dynamics", height = 220)
          ),
          box(title = "Patient Laboratory Profile",
              width = 5, status = "warning",
              tableOutput("tbl_labs")
          )
        )
      ),

      ## ── TAB 2: Drug PK ─────────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Eculizumab / Ravulizumab Plasma Levels",
              width = 6, status = "info",
              plotOutput("plot_c5i_pk", height = 300),
              helpText("Target trough: >35 μg/mL (eculizumab). Trough <35 = risk of breakthrough hemolysis.")
          ),
          box(title = "Iptacopan & Danicopan Plasma Levels",
              width = 6, status = "success",
              plotOutput("plot_oral_pk", height = 300),
              helpText("Iptacopan IC50 (Factor B): 0.05 μg/mL. Danicopan IC50 (Factor D): 0.02 μg/mL.")
          )
        ),
        fluidRow(
          box(title = "PK Summary (Steady-State Characteristics)",
              width = 12, status = "primary",
              HTML('<table class="table table-bordered table-sm">
              <thead><tr><th>Drug</th><th>Target</th><th>Dose</th><th>Route</th><th>t½</th><th>Dosing Interval</th><th>IC50/EC50</th></tr></thead>
              <tbody>
              <tr><td>Eculizumab</td><td>C5</td><td>900mg</td><td>IV</td><td>~11d</td><td>Q2W</td><td>Kd~0.1 nM</td></tr>
              <tr><td>Ravulizumab</td><td>C5</td><td>3300mg</td><td>IV</td><td>~49d</td><td>Q8W</td><td>Kd~0.1 nM</td></tr>
              <tr><td>Iptacopan</td><td>Factor B</td><td>200mg</td><td>PO</td><td>~14h</td><td>BID</td><td>IC50 ~50 nM</td></tr>
              <tr><td>Danicopan</td><td>Factor D</td><td>150mg</td><td>PO</td><td>~8h</td><td>TID</td><td>IC50 ~10 nM</td></tr>
              <tr><td>Crovalimab</td><td>C5</td><td>340mg</td><td>SQ</td><td>~30d</td><td>Q4W</td><td>Kd~0.3 nM</td></tr>
              <tr><td>Pegcetacoplan</td><td>C3</td><td>1080mg</td><td>SQ</td><td>~8d</td><td>2×/week</td><td>Kd~1 nM</td></tr>
              </tbody></table>')
          )
        )
      ),

      ## ── TAB 3: Complement ──────────────────────────────────────
      tabItem(tabName = "complement",
        fluidRow(
          box(title = "Complement C3 & C5 Levels",
              width = 6, status = "primary",
              plotOutput("plot_c3c5", height = 300)
          ),
          box(title = "C3b on PNH RBC Surface (EVH Driver)",
              width = 6, status = "danger",
              plotOutput("plot_c3b_mac", height = 300)
          )
        ),
        fluidRow(
          box(title = "Complement Pathway — Key Nodes",
              width = 12, status = "info",
              HTML('<div style="font-size:11px">
              <b>Alternative Pathway (AP):</b>
              C3 → [tick-over] → C3(H₂O) → [FB + FD] → C3 Convertase (C3bBb) → C3b amplification loop<br>
              <b>PNH Pathology:</b> No CD55 → C3bBb not decayed → massive C3b deposition on PNH RBCs<br>
              <b>C5 Convertase:</b> C3bBbC3b → C5 → C5a (inflammatory) + C5b → + C6/7/8/C9×18 → MAC<br>
              <b>PNH Pathology:</b> No CD59 → C9 polymerizes freely → lytic MAC pore → intravascular hemolysis<br>
              <b>EVH:</b> C3b/iC3b on surface → splenic/hepatic phagocytosis (not blocked by C5 inhibitors!)<br>
              <b>Drug Sites:</b> Eculizumab/Ravulizumab → block C5 (stop IVH but NOT EVH); Iptacopan → block FB (stop both); Pegcetacoplan → block C3 (stop all complement)
              </div>'),
              br(),
              plotOutput("plot_complement_overview", height = 250)
          )
        )
      ),

      ## ── TAB 4: Hemolysis Markers ───────────────────────────────
      tabItem(tabName = "hemolysis",
        fluidRow(
          box(title = "LDH (Intravascular Hemolysis Marker)",
              width = 6, status = "warning",
              plotOutput("plot_ldh", height = 280),
              helpText("LDH > 1.5× ULN = active IVH. ULN = 250 U/L. Target: <1× ULN on treatment.")
          ),
          box(title = "Free Plasma Hemoglobin & Haptoglobin",
              width = 6, status = "danger",
              plotOutput("plot_fhgb_hp", height = 280)
          )
        ),
        fluidRow(
          box(title = "Nitric Oxide & Thrombosis Risk",
              width = 6, status = "info",
              plotOutput("plot_no", height = 260)
          ),
          box(title = "Hemolysis Metrics Summary",
              width = 6, status = "primary",
              tableOutput("tbl_hemolysis_summary")
          )
        )
      ),

      ## ── TAB 5: Clinical Endpoints ──────────────────────────────
      tabItem(tabName = "endpoints",
        fluidRow(
          box(title = "Hemoglobin & Transfusion Independence",
              width = 6, status = "success",
              plotOutput("plot_hgb", height = 280),
              helpText("TI threshold: Hgb ≥ 12 g/dL without transfusion (APPLY-PNH primary endpoint)")
          ),
          box(title = "FACIT-Fatigue Score",
              width = 6, status = "info",
              plotOutput("plot_facit", height = 280),
              helpText("FACIT-Fatigue: 0–52 (higher = less fatigue). Clinically meaningful change: ≥ 3 points.")
          )
        ),
        fluidRow(
          box(title = "Key PNH Clinical Trials",
              width = 12, status = "warning",
              HTML('<table class="table table-bordered table-sm">
              <thead><tr><th>Trial</th><th>Drug</th><th>N</th><th>Duration</th><th>Key Endpoint</th><th>Result</th></tr></thead>
              <tbody>
              <tr><td>TRIUMPH</td><td>Eculizumab</td><td>87</td><td>26 wk</td><td>TI, LDH norm</td><td>TI 49% vs 0%; LDH norm 86%</td></tr>
              <tr><td>ALXN1210-301</td><td>Ravulizumab vs ECU</td><td>246</td><td>26 wk</td><td>LDH non-inferiority</td><td>TI 73.6% vs 66.1%</td></tr>
              <tr><td>APPLY-PNH</td><td>Iptacopan vs ECU</td><td>97</td><td>24 wk</td><td>TI, Hgb ≥ 2g/dL increase</td><td>TI 51.1% vs 0% (p<0.001)</td></tr>
              <tr><td>PEGASUS</td><td>Pegcetacoplan vs ECU</td><td>80</td><td>16 wk</td><td>Hgb change</td><td>+3.84 g/dL vs −0.83 g/dL</td></tr>
              <tr><td>GALAXY</td><td>Danicopan add-on</td><td>75</td><td>24 wk</td><td>Hgb change</td><td>+1.4 g/dL above C5i</td></tr>
              <tr><td>COMMODORE 1/2</td><td>Crovalimab</td><td>89+108</td><td>24 wk</td><td>LDH, TI</td><td>Non-inferior to ECU</td></tr>
              </tbody></table>')
          )
        )
      ),

      ## ── TAB 6: Scenario Comparison ─────────────────────────────
      tabItem(tabName = "scenario",
        fluidRow(
          box(title = "Run & Compare Multiple Scenarios",
              width = 12, status = "primary",
              fluidRow(
                column(3,
                  h5("Scenario A"),
                  checkboxGroupInput("scen_A", NULL,
                    choices = list("Eculizumab"="ECU","Ravulizumab"="RAV",
                                   "Iptacopan"="IPC","Danicopan"="DAN"),
                    selected = character(0))
                ),
                column(3,
                  h5("Scenario B"),
                  checkboxGroupInput("scen_B", NULL,
                    choices = list("Eculizumab"="ECU","Ravulizumab"="RAV",
                                   "Iptacopan"="IPC","Danicopan"="DAN"),
                    selected = "ECU")
                ),
                column(3,
                  h5("Scenario C"),
                  checkboxGroupInput("scen_C", NULL,
                    choices = list("Eculizumab"="ECU","Ravulizumab"="RAV",
                                   "Iptacopan"="IPC","Danicopan"="DAN"),
                    selected = "IPC")
                ),
                column(3,
                  h5("Scenario D"),
                  checkboxGroupInput("scen_D", NULL,
                    choices = list("Eculizumab"="ECU","Ravulizumab"="RAV",
                                   "Iptacopan"="IPC","Danicopan"="DAN"),
                    selected = c("ECU","DAN"))
                )
              ),
              actionButton("run_compare", "Compare Scenarios",
                           icon = icon("balance-scale"), class = "btn-primary")
          )
        ),
        fluidRow(
          box(title = "Hemoglobin Comparison",
              width = 6, status = "success",
              plotOutput("plot_scen_hgb", height = 280)),
          box(title = "LDH Comparison",
              width = 6, status = "warning",
              plotOutput("plot_scen_ldh", height = 280))
        ),
        fluidRow(
          box(title = "C3b on PNH RBC (EVH)",
              width = 6, status = "danger",
              plotOutput("plot_scen_c3b", height = 250)),
          box(title = "FACIT Fatigue Comparison",
              width = 6, status = "info",
              plotOutput("plot_scen_facit", height = 250))
        )
      ),

      ## ── TAB 7: Biomarkers ─────────────────────────────────────
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(title = "PNH Clone Size Over Time",
              width = 6, status = "danger",
              plotOutput("plot_clone", height = 280),
              helpText("PNH clone size reflects underlying disease burden. Treatment reduces hemolysis but clone persists (except after HSCT).")
          ),
          box(title = "Free C5 Monitoring (C5 Inhibitor Therapy)",
              width = 6, status = "info",
              plotOutput("plot_free_c5", height = 280),
              helpText("Free C5 target < 0.5 μg/mL on eculizumab/ravulizumab. Breakthrough if inadequate trough.")
          )
        ),
        fluidRow(
          box(title = "Haptoglobin Levels",
              width = 6, status = "warning",
              plotOutput("plot_hp", height = 250)
          ),
          box(title = "Biomarker Reference Ranges",
              width = 6, status = "primary",
              HTML('<table class="table table-sm">
              <tr><th>Biomarker</th><th>Normal</th><th>PNH (untreated)</th><th>Target on Rx</th></tr>
              <tr><td>LDH</td><td>&lt;250 U/L</td><td>3–10× ULN</td><td>&lt;1.5× ULN</td></tr>
              <tr><td>Hgb</td><td>12–17 g/dL</td><td>6–10 g/dL</td><td>≥12 g/dL (TI)</td></tr>
              <tr><td>Free Hgb</td><td>&lt;5 mg/dL</td><td>↑↑</td><td>&lt;5 mg/dL</td></tr>
              <tr><td>Haptoglobin</td><td>30–200 mg/dL</td><td>&lt;10 mg/dL</td><td>Normal range</td></tr>
              <tr><td>Reticulocytes</td><td>0.5–2%</td><td>5–20%</td><td>↓ on treatment</td></tr>
              <tr><td>PNH clone (RBC)</td><td>0%</td><td>30–95%</td><td>Monitor (unchanged)</td></tr>
              <tr><td>Free C5</td><td>~75 μg/mL</td><td>75 μg/mL</td><td>&lt;0.5 μg/mL (C5i)</td></tr>
              <tr><td>ECU trough</td><td>N/A</td><td>N/A</td><td>&gt;35 μg/mL</td></tr>
              </table>')
          )
        )
      ),

      ## ── TAB 8: About ──────────────────────────────────────────
      tabItem(tabName = "about",
        fluidRow(
          box(title = "About This Model",
              width = 12, status = "info",
              HTML('
              <h4>Paroxysmal Nocturnal Hemoglobinuria (PNH) — QSP Model</h4>
              <p><b>Disease Summary:</b> PNH is a rare clonal disorder caused by somatic PIGA mutation in
              hematopoietic stem cells, leading to GPI anchor deficiency on blood cells. The absence of
              CD55 (DAF) and CD59 (MIRL) from PNH cell surfaces renders them unable to regulate the
              complement alternative pathway, resulting in chronic hemolysis, thrombosis, and cytopenias.</p>
              <h5>Model Architecture (24 ODEs)</h5>
              <ul>
              <li><b>Hematopoiesis (4):</b> PNH/normal reticulocytes, PNH/normal RBCs</li>
              <li><b>Complement (4):</b> C3, C3b surface, C5, MAC</li>
              <li><b>Hemolysis outputs (4):</b> Free Hgb, haptoglobin, LDH, NO</li>
              <li><b>Eculizumab PK (3):</b> central, peripheral, C5:ECU complex</li>
              <li><b>Ravulizumab PK (3):</b> central, peripheral, C5:RAV complex</li>
              <li><b>Iptacopan PK (2):</b> gut, plasma</li>
              <li><b>Danicopan PK (2):</b> gut, plasma</li>
              </ul>
              <h5>Treatment Scenarios</h5>
              <ol>
              <li><b>S0 – Untreated:</b> Natural disease progression</li>
              <li><b>S1 – Eculizumab:</b> 600mg loading × 4, then 900mg q2w IV</li>
              <li><b>S2 – Ravulizumab:</b> 3000mg loading, then 3300mg q8w IV</li>
              <li><b>S3 – Iptacopan:</b> 200mg PO BID (Factor B inhibitor)</li>
              <li><b>S4 – ECU + Danicopan:</b> C5 inhibitor + add-on FD inhibitor for EVH</li>
              </ol>
              <h5>Key References</h5>
              <ul>
              <li>Hillmen P et al. NEJM 2006;355:1233 (eculizumab pivotal TRIUMPH trial)</li>
              <li>Kulasekararaj AG et al. Blood 2019 (ravulizumab)</li>
              <li>Peffault de Latour R et al. NEJM 2024;390:994 (iptacopan APPLY)</li>
              <li>Risitano AM et al. Blood 2014;124:3508 (complement biology)</li>
              <li>Parker C et al. Blood 2005;106:3699 (PNH diagnosis)</li>
              </ul>
              <hr>
              <p><i>Model generated by Claude Code Routine (CCR) | QSP Disease Model Library | 2026-06-24</i></p>
              ')
          )
        )
      )
    )
  )
)

## ─────────────────────────────────────────────────────────────────
## SERVER
## ─────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Compile model once
  mod <- mcode("pnh_shiny", pnh_model_code, quiet = TRUE)

  # Reactive simulation result
  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Running PNH simulation...", {
      run_sim(mod,
              f_pnh    = input$f_pnh / 100,
              use_ecu  = input$use_ecu,
              use_rav  = input$use_rav,
              use_ipc  = input$use_ipc,
              use_dan  = input$use_dan,
              end_time = input$sim_end)
    })
  }, ignoreNULL = FALSE)

  # Scenario comparison reactive
  scenario_data <- eventReactive(input$run_compare, {
    withProgress(message = "Comparing scenarios...", {
      f_pnh <- input$f_pnh / 100
      scenarios <- list(
        list(name = "A: Untreated",
             ecu = "ECU" %in% input$scen_A, rav = "RAV" %in% input$scen_A,
             ipc = "IPC" %in% input$scen_A, dan = "DAN" %in% input$scen_A),
        list(name = "B",
             ecu = "ECU" %in% input$scen_B, rav = "RAV" %in% input$scen_B,
             ipc = "IPC" %in% input$scen_B, dan = "DAN" %in% input$scen_B),
        list(name = "C",
             ecu = "ECU" %in% input$scen_C, rav = "RAV" %in% input$scen_C,
             ipc = "IPC" %in% input$scen_C, dan = "DAN" %in% input$scen_C),
        list(name = "D",
             ecu = "ECU" %in% input$scen_D, rav = "RAV" %in% input$scen_D,
             ipc = "IPC" %in% input$scen_D, dan = "DAN" %in% input$scen_D)
      )
      # Update names with drug combination
      for (i in seq_along(scenarios)) {
        drugs <- c()
        if (scenarios[[i]]$ecu) drugs <- c(drugs, "ECU")
        if (scenarios[[i]]$rav) drugs <- c(drugs, "RAV")
        if (scenarios[[i]]$ipc) drugs <- c(drugs, "IPC")
        if (scenarios[[i]]$dan) drugs <- c(drugs, "DAN")
        if (length(drugs) == 0) drugs <- "None"
        scenarios[[i]]$name <- paste0("Scen ", LETTERS[i], ": ", paste(drugs, collapse="+"))
      }

      bind_rows(lapply(scenarios, function(s) {
        run_sim(mod, f_pnh, s$ecu, s$rav, s$ipc, s$dan, input$sim_end) %>%
          mutate(Scenario = s$name)
      }))
    })
  })

  ## Value boxes
  output$vb_hgb <- renderValueBox({
    d <- sim_data()
    last_hgb <- round(tail(d$Hgb, 1), 1)
    status <- if (last_hgb >= 12) "green" else if (last_hgb >= 8) "yellow" else "red"
    valueBox(paste0(last_hgb, " g/dL"), "Hemoglobin (final)",
             icon = icon("tint"), color = status)
  })

  output$vb_ldh <- renderValueBox({
    d <- sim_data()
    last_ldh <- round(tail(d$LDH_ULN, 1), 1)
    status <- if (last_ldh < 1.5) "green" else if (last_ldh < 3) "yellow" else "red"
    valueBox(paste0(last_ldh, "× ULN"), "LDH (final)",
             icon = icon("flask"), color = status)
  })

  output$vb_clone <- renderValueBox({
    d <- sim_data()
    last_clone <- round(tail(d$PNH_pct, 1), 0)
    valueBox(paste0(last_clone, "%"), "PNH Clone (RBC)",
             icon = icon("dna"), color = "orange")
  })

  output$vb_facit <- renderValueBox({
    d <- sim_data()
    last_facit <- round(tail(d$FACIT, 1), 0)
    status <- if (last_facit >= 40) "green" else if (last_facit >= 25) "yellow" else "red"
    valueBox(last_facit, "FACIT-Fatigue",
             icon = icon("running"), color = status)
  })

  ## RBC dynamics plot
  output$plot_rbc_dynamics <- renderPlot({
    d <- sim_data()
    d %>%
      select(time, PNH_RBC, NL_RBC) %>%
      pivot_longer(-time, names_to = "type", values_to = "count") %>%
      ggplot(aes(x = time, y = count, fill = type)) +
      geom_area(alpha = 0.7, position = "stack") +
      scale_fill_manual(values = c("PNH_RBC" = "#e53935", "NL_RBC" = "#43a047"),
                        labels = c("PNH RBCs", "Normal RBCs")) +
      labs(x = "Time (days)", y = "RBC (M/μL)", fill = "") +
      theme_bw(base_size = 11)
  })

  ## Lab table
  output$tbl_labs <- renderTable({
    d <- sim_data()
    last <- tail(d, 1)
    data.frame(
      Parameter = c("Hemoglobin (g/dL)", "LDH (× ULN)", "Free Hgb (g/dL)",
                    "Haptoglobin (mg/dL)", "NO (%normal)", "PNH Clone (%)",
                    "FACIT Score", "Thrombo Risk"),
      Value = c(round(last$Hgb, 1), round(last$LDH_ULN, 2), round(last$fHgb_out, 3),
                round(last$Haptoglobin_out, 0), round(last$NO_pct, 0), round(last$PNH_pct, 1),
                round(last$FACIT, 0), round(last$Thrombo, 2)),
      Status = c(
        ifelse(last$Hgb >= 12, "Normal", ifelse(last$Hgb >= 8, "Moderate anemia", "Severe anemia")),
        ifelse(last$LDH_ULN < 1.5, "Target achieved", "Active IVH"),
        ifelse(last$fHgb_out < 0.05, "Normal", "Elevated"),
        ifelse(last$Haptoglobin_out > 30, "Normal", "Depleted"),
        ifelse(last$NO_pct > 80, "Normal", "Reduced"),
        paste0(round(last$PNH_pct, 1), "% PNH"),
        ifelse(last$FACIT > 40, "Mild fatigue", ifelse(last$FACIT > 25, "Moderate", "Severe")),
        ifelse(last$Thrombo < 0.3, "Low", ifelse(last$Thrombo < 0.6, "Moderate", "High"))
      )
    )
  }, striped = TRUE, hover = TRUE)

  ## Drug PK plots
  output$plot_c5i_pk <- renderPlot({
    d <- sim_data()
    d %>%
      select(time, ECU_ug, RAV_ug) %>%
      pivot_longer(-time) %>%
      ggplot(aes(x = time, y = value, color = name)) +
      geom_line(size = 1.1) +
      geom_hline(yintercept = 35, linetype = "dashed", color = "red") +
      annotate("text", x = max(d$time)*0.8, y = 40, label = "Target trough 35 μg/mL",
               color = "red", size = 3) +
      scale_color_manual(values = c("ECU_ug" = "#1565c0", "RAV_ug" = "#0288d1"),
                         labels = c("Eculizumab", "Ravulizumab")) +
      labs(x = "Time (days)", y = "Concentration (μg/mL)", color = "") +
      theme_bw(base_size = 11)
  })

  output$plot_oral_pk <- renderPlot({
    d <- sim_data()
    d %>%
      select(time, IPC_ug, DAN_ug) %>%
      pivot_longer(-time) %>%
      ggplot(aes(x = time, y = value, color = name)) +
      geom_line(size = 0.8) +
      scale_color_manual(values = c("IPC_ug" = "#2e7d32", "DAN_ug" = "#827717"),
                         labels = c("Iptacopan", "Danicopan")) +
      labs(x = "Time (days)", y = "Concentration (μg/mL)", color = "") +
      theme_bw(base_size = 11)
  })

  ## Complement plots
  output$plot_c3c5 <- renderPlot({
    d <- sim_data()
    d %>%
      select(time, C3, C5_free_out) %>%
      pivot_longer(-time) %>%
      ggplot(aes(x = time, y = value, color = name)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = c("C3" = "#1976d2", "C5_free_out" = "#c62828"),
                         labels = c("C3 (μg/mL)", "Free C5 (μg/mL)")) +
      labs(x = "Time (days)", y = "Concentration (μg/mL)", color = "") +
      theme_bw(base_size = 11)
  })

  output$plot_c3b_mac <- renderPlot({
    d <- sim_data()
    d %>%
      select(time, C3b_out, MAC_out) %>%
      pivot_longer(-time) %>%
      ggplot(aes(x = time, y = value, color = name)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = c("C3b_out" = "#e65100", "MAC_out" = "#b71c1c"),
                         labels = c("C3b on PNH RBC", "MAC (C5b-9)")) +
      labs(x = "Time (days)", y = "Relative units", color = "") +
      theme_bw(base_size = 11)
  })

  output$plot_complement_overview <- renderPlot({
    d <- sim_data()
    d %>%
      ggplot(aes(x = time)) +
      geom_line(aes(y = C3b_out, color = "C3b (EVH driver)"), size = 1.1) +
      geom_line(aes(y = MAC_out * 0.5, color = "MAC/2 (IVH)"), size = 1.1) +
      scale_color_manual(values = c("C3b (EVH driver)" = "#e65100",
                                    "MAC/2 (IVH)" = "#b71c1c")) +
      labs(x = "Time (days)", y = "Complement activity (rel.)", color = "") +
      theme_bw(base_size = 11)
  })

  ## Hemolysis markers
  output$plot_ldh <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = time, y = LDH_ULN)) +
      geom_line(size = 1.2, color = "#f57f17") +
      geom_hline(yintercept = 1.5, linetype = "dashed", color = "red") +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "green4") +
      annotate("text", x = max(d$time)*0.8, y = 1.65, label = "IVH threshold", color = "red", size = 3) +
      annotate("text", x = max(d$time)*0.8, y = 1.1,  label = "Target", color = "green4", size = 3) +
      labs(x = "Time (days)", y = "LDH (× ULN)") +
      theme_bw(base_size = 12)
  })

  output$plot_fhgb_hp <- renderPlot({
    d <- sim_data()
    p <- ggplot(d, aes(x = time)) +
      geom_line(aes(y = fHgb_out, color = "Free Hgb (g/dL)"), size = 1.1) +
      geom_line(aes(y = Haptoglobin_out / 100, color = "Haptoglobin (×100 mg/dL)"), size = 1.1) +
      scale_color_manual(values = c("Free Hgb (g/dL)" = "#e53935",
                                    "Haptoglobin (×100 mg/dL)" = "#ff9800")) +
      labs(x = "Time (days)", y = "Level (relative scale)", color = "") +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom")
    p
  })

  output$plot_no <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = time, y = NO_pct)) +
      geom_line(size = 1.2, color = "#283593") +
      geom_hline(yintercept = 80, linetype = "dashed", color = "gray40") +
      ylim(0, 120) +
      labs(x = "Time (days)", y = "NO (% of normal)") +
      annotate("text", x = max(d$time)*0.7, y = 85, label = "Normal range floor", size = 3) +
      theme_bw(base_size = 12)
  })

  output$tbl_hemolysis_summary <- renderTable({
    d <- sim_data()
    ss <- filter(d, time >= max(d$time) * 0.5)
    data.frame(
      Metric = c("LDH (× ULN)", "Free Hgb (g/dL)", "Haptoglobin (mg/dL)", "NO (% normal)"),
      Baseline = c(round(d$LDH_ULN[1], 2), round(d$fHgb_out[1], 3),
                   round(d$Haptoglobin_out[1], 0), round(d$NO_pct[1], 0)),
      SS_Mean = c(round(mean(ss$LDH_ULN), 2), round(mean(ss$fHgb_out), 3),
                  round(mean(ss$Haptoglobin_out), 0), round(mean(ss$NO_pct), 0))
    )
  }, striped = TRUE)

  ## Clinical endpoints
  output$plot_hgb <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = time, y = Hgb)) +
      geom_line(size = 1.2, color = "#c62828") +
      geom_hline(yintercept = 12, linetype = "dashed", color = "green4", linewidth = 1) +
      geom_hline(yintercept = 8,  linetype = "dashed", color = "orange") +
      annotate("text", x = max(d$time)*0.75, y = 12.4, label = "TI threshold (12 g/dL)",
               color = "green4", size = 3.5) +
      annotate("text", x = max(d$time)*0.75, y = 8.4, label = "Transfusion trigger",
               color = "orange", size = 3.5) +
      ylim(4, 16) +
      labs(x = "Time (days)", y = "Hemoglobin (g/dL)") +
      theme_bw(base_size = 12)
  })

  output$plot_facit <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = time, y = FACIT)) +
      geom_line(size = 1.2, color = "#1565c0") +
      geom_hline(yintercept = 44, linetype = "dashed", color = "green4") +
      annotate("text", x = max(d$time)*0.7, y = 45.5, label = "Normal FACIT (≥44)", size = 3) +
      ylim(0, 52) +
      labs(x = "Time (days)", y = "FACIT-Fatigue Score (0–52)") +
      theme_bw(base_size = 12)
  })

  ## Scenario comparison plots
  scen_colors <- c("#d32f2f","#1565c0","#2e7d32","#6a1b9a","#e65100","#00838f")

  output$plot_scen_hgb <- renderPlot({
    d <- scenario_data()
    ggplot(d, aes(x = time, y = Hgb, color = Scenario)) +
      geom_line(size = 1.1) +
      geom_hline(yintercept = 12, linetype = "dashed", color = "gray40") +
      scale_color_brewer(palette = "Set1") +
      labs(x = "Time (days)", y = "Hgb (g/dL)", color = NULL) +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })

  output$plot_scen_ldh <- renderPlot({
    d <- scenario_data()
    ggplot(d, aes(x = time, y = LDH_ULN, color = Scenario)) +
      geom_line(size = 1.1) +
      geom_hline(yintercept = 1.5, linetype = "dashed") +
      scale_color_brewer(palette = "Set1") +
      labs(x = "Time (days)", y = "LDH (× ULN)", color = NULL) +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })

  output$plot_scen_c3b <- renderPlot({
    d <- scenario_data()
    ggplot(d, aes(x = time, y = C3b_out, color = Scenario)) +
      geom_line(size = 1.1) +
      scale_color_brewer(palette = "Set1") +
      labs(x = "Time (days)", y = "C3b on PNH RBC (rel.)", color = NULL) +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })

  output$plot_scen_facit <- renderPlot({
    d <- scenario_data()
    ggplot(d, aes(x = time, y = FACIT, color = Scenario)) +
      geom_line(size = 1.1) +
      scale_color_brewer(palette = "Set1") +
      labs(x = "Time (days)", y = "FACIT Score (0–52)", color = NULL) +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })

  ## Biomarkers
  output$plot_clone <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = time, y = PNH_pct)) +
      geom_line(size = 1.2, color = "#e53935") +
      ylim(0, 100) +
      labs(x = "Time (days)", y = "PNH RBC Clone (%)") +
      theme_bw(base_size = 12)
  })

  output$plot_free_c5 <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = time, y = C5_free_out)) +
      geom_line(size = 1.2, color = "#c62828") +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
      annotate("text", x = max(d$time)*0.7, y = 1.5, label = "Target <0.5 μg/mL", size = 3, color = "red") +
      labs(x = "Time (days)", y = "Free C5 (μg/mL)") +
      theme_bw(base_size = 12)
  })

  output$plot_hp <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = time, y = Haptoglobin_out)) +
      geom_line(size = 1.2, color = "#ff9800") +
      geom_hline(yintercept = 30, linetype = "dashed", color = "red") +
      annotate("text", x = max(d$time)*0.7, y = 40, label = "Lower normal (30 mg/dL)", size = 3) +
      labs(x = "Time (days)", y = "Haptoglobin (mg/dL)") +
      theme_bw(base_size = 12)
  })
}

## ─────────────────────────────────────────────────────────────────
## RUN APP
## ─────────────────────────────────────────────────────────────────
shinyApp(ui, server)
