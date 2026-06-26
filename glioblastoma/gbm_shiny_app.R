# ============================================================
# GBM QSP Shiny Dashboard
# Glioblastoma Multiforme — Interactive Simulation
# ============================================================
# Tabs: 7
#  1. Patient Profile & Treatment Setup
#  2. Drug PK (TMZ, Bevacizumab, Anti-PD1)
#  3. DNA Damage & MGMT Repair
#  4. Tumor Cell Dynamics
#  5. Scenario Comparison
#  6. Tumor Microenvironment & Biomarkers
#  7. Clinical Endpoints & Survival
# ============================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(shinydashboard)
library(DT)

# ---- Model code (same as mrgsolve model file) ----
gbm_code <- '
$PROB GBM QSP Shiny Model

$PARAM @annotated
ka=2.04:h-1 CL_tmz=11.4:L/h V1_tmz=22.5:L Q_tmz=3.0:L/h V2_tmz=50.0:L
f_brain=0.28:- k_hyd=0.6:h-1 k_O6=0.10:h-1 k_O6deg=0.02:h-1
DOSE_TMZ=0:mg/m2 BSA=1.8:m2
MGMT_meth=1:- IDH_status=0:-
kMGMT_high=0.40:h-1 kMGMT_low=0.05:h-1
CL_bev=0.207:L/day V1_bev=2.92:L Q_bev=0.553:L/day V2_bev=3.28:L
DOSE_BEV=0:mg/kg WT=70:kg
kon_bv=0.005:- koff_bv=0.0002:- kdeg_cx=0.005:-
VEGF0=45.0:pg/mL kprod_V=5.0:- kdeg_V=0.08:day-1
CL_pd1=0.22:L/day V_pd1=3.28:L DOSE_PD1=0:mg EC50_pd1=0.003:mg/L
kg_ts=0.003:day-1 kg_tr=0.0015:day-1 kg_gsc=0.001:day-1
K_cap=1e10:cells Ts0_init=9.5e8:cells Tr0_init=5e6:cells GSC0_init=2e5:cells
kd_tmz_ts=0.25:- kd_tmz_gsc=0.06:- k_resist=0.0003:day-1
alpha_gbm=0.30:Gy-1 beta_gbm=0.030:Gy-2 RT_dose_fx=2.0:Gy
RT_active=0:- OER_hyp=1.8:- frac_hyp=0.30:-
CD8_0=1.0:- Treg_0=0.5:- TAM_0=2.0:- NV0=1.0:-
k_prim=0.05:- k_exh=0.12:- k_kill_cd8=0.04:- k_sup_treg=0.10:- k_tam_sup=0.06:-
k_treg_exp=0.08:- k_treg_deg=0.05:- k_tam_rec=0.01:- k_tam_deg=0.03:-
kgNV=0.04:- kdNV=0.02:-
CONCURRENT_RT=0:- ADJUVANT_TMZ=0:- BEV_ON=0:- PD1_ON=0:- TTF_ON=0:-
TTF_kill_fx=0.06:day-1

$CMT Gut Cp_tmz Cp2_tmz Cbrain O6MeG BEV_Cp BEV_Cp2 VEGF_free BEV_VEGF
APD1_Cp PD1_occ Ts Tr GSC CD8_eff Treg_c TAM_M2 NV

$MAIN
double kMGMT = kMGMT_low + (kMGMT_high - kMGMT_low) * (1.0 - MGMT_meth);
double idh_gr = 1.0 - 0.12 * IDH_status;
double alpha_eff = alpha_gbm * (1.0 + 0.4 * O6MeG / (0.5 + O6MeG));
double RT_OER_factor = 1.0 - frac_hyp * (1.0 - 1.0/OER_hyp);
double RT_kill_frac = 0.0;
if (RT_active > 0.5 || CONCURRENT_RT > 0.5) {
  double SF_RT = exp(-alpha_eff * RT_dose_fx - beta_gbm * RT_dose_fx * RT_dose_fx);
  RT_kill_frac = (1.0 - SF_RT) * RT_OER_factor;
}
double TotalCells = Ts + Tr + GSC;
double TumorVol_mL = TotalCells / 1e9;
double TumorDiam_cm = 2.0 * pow(0.238732 * TumorVol_mL, 0.333333);
double Gompertz_Ts  = -kg_ts  * idh_gr * log((TotalCells + 1.0) / K_cap);
double Gompertz_Tr  = -kg_tr  * idh_gr * log((TotalCells + 1.0) / K_cap);
double Gompertz_GSC = -kg_gsc * idh_gr * log((TotalCells + 1.0) / K_cap);
double kill_Ts_TMZ  = kd_tmz_ts * O6MeG;
double kill_Tr_TMZ  = kd_tmz_ts * O6MeG * 0.1;
double kill_GSC_TMZ = kd_tmz_ts * O6MeG * kd_tmz_gsc;
double exhaust_rate = k_exh * (1.0 - PD1_occ);
double kill_CD8 = k_kill_cd8 * CD8_eff * (1.0 - 0.7 * Treg_c / (0.5 + Treg_c))
                                         * (1.0 - 0.5 * TAM_M2 / (2.0 + TAM_M2));
double kill_TTF = TTF_kill_fx * TTF_ON;
double VEGF_NV_eff = VEGF_free / (50.0 + VEGF_free);
double Ag_load = TotalCells / (Ts0_init + Tr0_init + GSC0_init);

if (NEWIND <= 1) {
  Ts_0 = Ts0_init; Tr_0 = Tr0_init; GSC_0 = GSC0_init;
  VEGF_free_0 = VEGF0; CD8_eff_0 = CD8_0;
  Treg_c_0 = Treg_0; TAM_M2_0 = TAM_0; NV_0 = NV0; PD1_occ_0 = 0.0;
}

$ODE
dxdt_Gut     = -ka * Gut;
dxdt_Cp_tmz  = ka * Gut / V1_tmz - (CL_tmz/V1_tmz)*Cp_tmz
               - (Q_tmz/V1_tmz)*Cp_tmz + (Q_tmz/V2_tmz)*Cp2_tmz;
dxdt_Cp2_tmz = (Q_tmz/V1_tmz)*Cp_tmz - (Q_tmz/V2_tmz)*Cp2_tmz;
dxdt_Cbrain   = (CL_tmz/V1_tmz)*Cp_tmz*f_brain - k_hyd*Cbrain;
dxdt_O6MeG    = k_O6*Cbrain - kMGMT*O6MeG - k_O6deg*O6MeG;
double rate_cx_on = kon_bv*BEV_Cp*VEGF_free;
double rate_cx_off= koff_bv*BEV_VEGF;
dxdt_BEV_Cp  = -(CL_bev/V1_bev)*BEV_Cp - (Q_bev/V1_bev)*BEV_Cp
               + (Q_bev/V2_bev)*BEV_Cp2 - rate_cx_on + rate_cx_off;
dxdt_BEV_Cp2 = (Q_bev/V1_bev)*BEV_Cp - (Q_bev/V2_bev)*BEV_Cp2;
double VEGF_prod = kprod_V*(TotalCells/(Ts0_init+1e-3))*(1.0+0.3*NV);
dxdt_VEGF_free = VEGF_prod - kdeg_V*VEGF_free - rate_cx_on + rate_cx_off;
dxdt_BEV_VEGF  = rate_cx_on - rate_cx_off - kdeg_cx*BEV_VEGF;
dxdt_APD1_Cp = -(CL_pd1/V_pd1)*APD1_Cp;
double occ_target = (PD1_ON > 0.5) ? APD1_Cp/(EC50_pd1 + APD1_Cp) : 0.0;
dxdt_PD1_occ = 0.15*(occ_target - PD1_occ);
dxdt_Ts = Ts*(Gompertz_Ts - kill_Ts_TMZ - RT_kill_frac*kg_ts - kill_CD8 - kill_TTF - k_resist) + 0.005*GSC;
dxdt_Tr = Tr*(Gompertz_Tr - kill_Tr_TMZ - kill_CD8*0.3 - kill_TTF*0.5) + k_resist*Ts;
dxdt_GSC = GSC*(Gompertz_GSC - kill_GSC_TMZ - kill_CD8*0.15 - kill_TTF*0.3);
dxdt_CD8_eff = k_prim*Ag_load/(1.0+Ag_load) - exhaust_rate*CD8_eff
               - k_sup_treg*Treg_c*CD8_eff - k_tam_sup*TAM_M2*CD8_eff;
dxdt_Treg_c = k_treg_exp*TAM_M2*(1.0+0.5*Ag_load) - k_treg_deg*Treg_c;
dxdt_TAM_M2 = k_tam_rec*(TotalCells/Ts0_init) - k_tam_deg*TAM_M2;
dxdt_NV = kgNV*VEGF_NV_eff - kdNV*NV;

$TABLE
capture TumorVol_mL=TumorVol_mL; capture TumorDiam_cm=TumorDiam_cm;
capture TotalCells=TotalCells; capture O6MeG_cap=O6MeG;
capture VEGF_cap=VEGF_free; capture NV_cap=NV;
capture PD1_occ_cap=PD1_occ; capture CD8_cap=CD8_eff;
capture Treg_cap=Treg_c; capture TAM_cap=TAM_M2;
capture kill_TMZ_cap=kill_Ts_TMZ; capture kill_CD8_cap=kill_CD8;
'

mod_shiny <- mcode("GBM_shiny", gbm_code)

# ---- Helper: build events ----
build_ev_shiny <- function(tmz_dose=75, adj_dose=150, adj_cycles=6,
                            adj_start=57, bev_dose=0, bev_start=0, bev_n=0,
                            pd1_dose=0, pd1_start=0, pd1_n=0,
                            conc_days=42, BSA=1.8, WT=70) {
  evts <- data.frame(time=numeric(0), cmt=integer(0), amt=numeric(0),
                     rate=numeric(0), evid=integer(0))
  # Concurrent TMZ
  if (conc_days > 0 && tmz_dose > 0) {
    mg <- tmz_dose * BSA
    times <- seq(0, (conc_days-1)*24, by=24)
    evts <- rbind(evts, data.frame(time=times, cmt=1L, amt=mg, rate=0, evid=1L))
  }
  # Adjuvant TMZ
  for (cyc in seq_len(adj_cycles)) {
    d0 <- adj_start + (cyc-1)*28
    for (dd in 0:4) {
      evts <- rbind(evts, data.frame(
        time=(d0+dd)*24, cmt=1L, amt=adj_dose*BSA, rate=0, evid=1L))
    }
  }
  # Bevacizumab
  if (bev_n > 0 && bev_dose > 0) {
    mg_bev <- bev_dose * WT
    for (ii in 0:(bev_n-1)) {
      evts <- rbind(evts, data.frame(
        time=(bev_start+ii*14)*24, cmt=6L, amt=mg_bev,
        rate=mg_bev/1.5, evid=1L))
    }
  }
  # Anti-PD1
  if (pd1_n > 0 && pd1_dose > 0) {
    for (ii in 0:(pd1_n-1)) {
      evts <- rbind(evts, data.frame(
        time=(pd1_start+ii*21)*24, cmt=10L, amt=pd1_dose,
        rate=pd1_dose/0.5, evid=1L))
    }
  }
  if (nrow(evts) > 0) evts <- evts[order(evts$time), ]
  evts
}

# ============================================================
# UI
# ============================================================
ui <- dashboardPage(
  skin="blue",
  dashboardHeader(title="GBM QSP Dashboard", titleWidth=280),
  dashboardSidebar(
    width=280,
    sidebarMenu(
      menuItem("Patient Profile",         tabName="tab_patient",   icon=icon("user")),
      menuItem("Drug PK",                 tabName="tab_pk",        icon=icon("flask")),
      menuItem("DNA Damage & MGMT",       tabName="tab_dna",       icon=icon("dna")),
      menuItem("Tumor Dynamics",          tabName="tab_tumor",     icon=icon("chart-line")),
      menuItem("Scenario Comparison",     tabName="tab_scenario",  icon=icon("layer-group")),
      menuItem("TME & Biomarkers",        tabName="tab_tme",       icon=icon("microscope")),
      menuItem("Clinical Endpoints",      tabName="tab_endpoints", icon=icon("heartbeat"))
    ),
    hr(),
    h5("Global Parameters", style="padding-left:15px;color:#ccc;"),
    sliderInput("sim_days","Simulation Duration (days)",
                min=90, max=730, value=365, step=30),
    sliderInput("BSA_val","BSA (m²)", min=1.2, max=2.5, value=1.8, step=0.1),
    sliderInput("WT_val","Body Weight (kg)", min=40, max=120, value=70, step=5),
    selectInput("MGMT_status","MGMT Methylation",
                choices=c("Methylated (~45%)"=1,"Unmethylated (~55%)"=0), selected=1),
    selectInput("IDH_status","IDH Status",
                choices=c("IDH Wild-type (~90%)"=0,"IDH Mutant (~10%)"=1), selected=0),
    actionButton("run_sim","▶  Run Simulation", class="btn-primary btn-block",
                 style="margin:10px 15px; width:calc(100%-30px);")
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .box-title { font-weight: bold; }
      .shiny-plot-output { min-height: 380px; }
    "))),
    tabItems(

      # ------- TAB 1: Patient Profile -------
      tabItem(tabName="tab_patient",
        fluidRow(
          box(title="Treatment Selection", status="primary", solidHeader=TRUE, width=6,
            selectInput("treatment_arm","Treatment Regimen",
              choices=c(
                "Untreated (Control)"="ctrl",
                "Stupp: Concurrent TMZ+RT (75 mg/m²/d)"="stupp_conc",
                "Stupp: Full (Concurrent + Adjuvant 150-200 mg/m²)"="stupp_full",
                "Stupp + Bevacizumab (AVAGLIO-like)"="stupp_bev",
                "Stupp + Tumor Treating Fields (EF-14)"="stupp_ttf",
                "Pembrolizumab + TMZ (recurrent GBM)"="pembro",
                "Bevacizumab Salvage (recurrent GBM)"="bev_salvage"
              ), selected="stupp_full"),
            hr(),
            h5("Concurrent TMZ Phase"),
            sliderInput("conc_dose","TMZ dose (mg/m²/d)", 60, 100, 75, 5),
            sliderInput("conc_days","Duration (days)", 21, 60, 42, 7),
            h5("Adjuvant TMZ Phase"),
            sliderInput("adj_dose","TMZ dose (mg/m²)", 100, 200, 150, 25),
            sliderInput("adj_cycles","Number of cycles", 1, 12, 6, 1)
          ),
          box(title="Patient Characteristics", status="info", solidHeader=TRUE, width=6,
            sliderInput("age_val","Age (years)", 20, 80, 58, 1),
            selectInput("kps_val","KPS Score",
                        choices=c("100 (Normal)"=100,"90"=90,"80"=80,"70"=70,"60"=60),
                        selected=80),
            selectInput("eor_val","Extent of Resection",
                        choices=c("Gross Total (GTR, >95%)"="GTR",
                                  "Sub-total (STR, 60-95%)"="STR",
                                  "Biopsy only (<10%)"="biopsy"),
                        selected="GTR"),
            sliderInput("init_cells","Initial Tumor Size (×10⁸ cells)",
                        min=0.1, max=50, value=9.5, step=0.5),
            h5("Additional Treatment"),
            checkboxInput("add_bev","Add Bevacizumab (10 mg/kg q14d)", FALSE),
            checkboxInput("add_pd1","Add Anti-PD1 Pembrolizumab (200 mg q3w)", FALSE),
            checkboxInput("add_ttf","Add Tumor Treating Fields (TTF)", FALSE)
          )
        ),
        fluidRow(
          box(title="Clinical Parameter Summary", width=12, status="warning",
            DTOutput("pt_summary_table"))
        )
      ),

      # ------- TAB 2: Drug PK -------
      tabItem(tabName="tab_pk",
        fluidRow(
          box(title="TMZ Plasma & Brain PK", status="primary", solidHeader=TRUE, width=8,
            plotOutput("pk_tmz_plot", height="360px")),
          box(title="PK Parameters", status="info", width=4,
            h5("TMZ PK (Ostermann 2004)"),
            HTML("<table class='table table-condensed'>
                 <tr><td>CL:</td><td>11.4 L/h</td></tr>
                 <tr><td>V₁:</td><td>22.5 L</td></tr>
                 <tr><td>V₂:</td><td>50.0 L</td></tr>
                 <tr><td>Q:</td><td>3.0 L/h</td></tr>
                 <tr><td>t½:</td><td>~1.8 h</td></tr>
                 <tr><td>Kp,brain:</td><td>0.28</td></tr>
                 </table>"),
            hr(),
            h5("Bevacizumab PK (Lu 2008)"),
            HTML("<table class='table table-condensed'>
                 <tr><td>CL:</td><td>0.207 L/day</td></tr>
                 <tr><td>V₁:</td><td>2.92 L</td></tr>
                 <tr><td>t½:</td><td>~20 days</td></tr>
                 </table>")
          )
        ),
        fluidRow(
          box(title="Bevacizumab PK & VEGF Binding", status="success",
              solidHeader=TRUE, width=6,
            plotOutput("pk_bev_plot", height="320px")),
          box(title="Anti-PD1 (Pembrolizumab) PK", status="warning",
              solidHeader=TRUE, width=6,
            plotOutput("pk_pd1_plot", height="320px"))
        )
      ),

      # ------- TAB 3: DNA Damage -------
      tabItem(tabName="tab_dna",
        fluidRow(
          box(title="O6-MeG Lesion Dynamics (TMZ Mechanism)", status="danger",
              solidHeader=TRUE, width=8,
            plotOutput("dna_lesion_plot", height="360px")),
          box(title="MGMT Mechanism", width=4, status="info",
            HTML("
            <h5>TMZ → MTIC → O6-MeG</h5>
            <p>TMZ undergoes spontaneous pH-dependent hydrolysis
            to MTIC, which methylates guanine at the O6 position.</p>
            <hr/>
            <h5>MGMT Repair</h5>
            <p><b>Unmethylated MGMT:</b> High MGMT expression directly
            reverses O6-MeG, conferring TMZ resistance.</p>
            <p><b>Methylated MGMT:</b> Promoter silencing → low MGMT →
            O6-MeG accumulates → MMR-mediated DSBs → apoptosis.</p>
            <hr/>
            <h5>MMR Pathway</h5>
            <p>MSH2/MSH6 sensor → MLH1/PMS2 → futile repair cycles
            → replication-associated DSBs → cell death.</p>
            ")
          )
        ),
        fluidRow(
          box(title="MGMT Status Comparison", status="warning", solidHeader=TRUE, width=12,
            plotOutput("mgmt_compare_plot", height="320px"))
        )
      ),

      # ------- TAB 4: Tumor Dynamics -------
      tabItem(tabName="tab_tumor",
        fluidRow(
          box(title="Tumor Cell Population (Ts / Tr / GSC)", status="danger",
              solidHeader=TRUE, width=8,
            plotOutput("tumor_cells_plot", height="380px")),
          box(title="Tumor Model Assumptions", width=4, status="info",
            HTML("
            <h5>Gompertz Growth Model</h5>
            <p>dTs/dt = Ts × [−kg × log(N/K) − kill_TMZ − kill_RT − kill_CD8 − k_resist]</p>
            <hr/>
            <h5>Three-Compartment Tumor</h5>
            <ul>
            <li><b>Ts (Sensitive):</b> ~9.5×10⁸ cells initial; killed by TMZ, RT, CD8</li>
            <li><b>Tr (Resistant):</b> ~5×10⁶ cells; MGMT-high/PTEN-loss subclone; low TMZ kill</li>
            <li><b>GSC (Stem):</b> ~2×10⁵ cells; drives recurrence; most resistant</li>
            </ul>
            <hr/>
            <h5>Resistance Mechanisms</h5>
            <ul>
            <li>Ts→Tr epigenetic transition (k=0.0003/day)</li>
            <li>MGMT-independent resistance pathways</li>
            <li>GSC repopulation after treatment</li>
            </ul>
            ")
          )
        ),
        fluidRow(
          box(title="Tumor Volume Trajectory (log scale)", status="primary",
              solidHeader=TRUE, width=12,
            plotOutput("tumor_vol_plot", height="280px"))
        )
      ),

      # ------- TAB 5: Scenario Comparison -------
      tabItem(tabName="tab_scenario",
        fluidRow(
          box(title="All-Scenario Tumor Volume Comparison", status="primary",
              solidHeader=TRUE, width=12,
            plotOutput("scenario_plot", height="420px"))
        ),
        fluidRow(
          box(title="Treatment Efficacy Summary", status="success",
              solidHeader=TRUE, width=12,
            DTOutput("scenario_table"))
        )
      ),

      # ------- TAB 6: TME & Biomarkers -------
      tabItem(tabName="tab_tme",
        fluidRow(
          box(title="Tumor Immune Microenvironment", status="success",
              solidHeader=TRUE, width=8,
            plotOutput("tme_plot", height="360px")),
          box(title="Key TME Components", width=4, status="info",
            HTML("
            <h5>CD8+ Effector T Cells</h5>
            <p>Primed by tumor Ag; exhausted via PD-1 pathway;
            suppressed by Tregs and M2-TAMs. The main anti-tumor effectors.</p>
            <hr/>
            <h5>Regulatory T Cells (Tregs)</h5>
            <p>Expanded by TGF-β from TAMs. Suppress CD8 T cell killing.
            Inversely correlated with survival in GBM.</p>
            <hr/>
            <h5>M2-TAMs (Glioma-Associated Macrophages)</h5>
            <p>~90% of tumor-infiltrating macrophages are M2-polarized.
            Secrete TGF-β, IL-10, VEGF, IDO1, ARG1 → immunosuppression.</p>
            ")
          )
        ),
        fluidRow(
          box(title="VEGF & Neovascularization", status="warning",
              solidHeader=TRUE, width=6,
            plotOutput("vegf_nv_plot", height="300px")),
          box(title="PD-1 Occupancy (Anti-PD1 Tx)", status="danger",
              solidHeader=TRUE, width=6,
            plotOutput("pd1_occ_plot", height="300px"))
        )
      ),

      # ------- TAB 7: Clinical Endpoints -------
      tabItem(tabName="tab_endpoints",
        fluidRow(
          box(title="Tumor Diameter Over Time (RANO)", status="primary",
              solidHeader=TRUE, width=8,
            plotOutput("diam_plot", height="360px")),
          box(title="RANO & Survival Reference", width=4, status="info",
            HTML("
            <h5>RANO Criteria</h5>
            <ul>
            <li><b>CR:</b> No enhancement, stable T2/FLAIR, no steroids</li>
            <li><b>PR:</b> ≥50% ↓ in SLD of enhancing lesion</li>
            <li><b>SD:</b> Neither PR nor PD</li>
            <li><b>PD:</b> ≥25% ↑ in SLD or new lesion</li>
            </ul>
            <hr/>
            <h5>Clinical Benchmarks (Stupp 2005 NEJM)</h5>
            <table class='table table-condensed table-bordered'>
            <tr><th>Subgroup</th><th>Median OS</th></tr>
            <tr><td>TMZ+RT overall</td><td>14.6 mo</td></tr>
            <tr><td>MGMT methylated</td><td>21.7 mo</td></tr>
            <tr><td>MGMT unmethylated</td><td>12.6 mo</td></tr>
            <tr><td>IDH-mutant GBM</td><td>~31 mo</td></tr>
            <tr><td>TTF+TMZ (EF-14)</td><td>20.9 mo</td></tr>
            </table>
            ")
          )
        ),
        fluidRow(
          box(title="Tumor Kill Rate Components", status="danger",
              solidHeader=TRUE, width=12,
            plotOutput("kill_components_plot", height="300px"))
        )
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  # --- Reactive simulation ---
  sim_result <- eventReactive(input$run_sim, {
    withProgress(message="Running GBM simulation...", {
      mgmt <- as.numeric(input$MGMT_status)
      idh  <- as.numeric(input$IDH_status)
      BSA  <- input$BSA_val
      WT   <- input$WT_val
      days <- input$sim_days

      add_bev <- input$add_bev
      add_pd1 <- input$add_pd1
      add_ttf <- input$add_ttf

      flags <- list(
        CONCURRENT_RT = as.numeric(input$treatment_arm %in%
                          c("stupp_conc","stupp_full","stupp_bev","stupp_ttf")),
        ADJUVANT_TMZ  = as.numeric(input$treatment_arm %in%
                          c("stupp_full","stupp_bev","stupp_ttf","pembro")),
        BEV_ON    = as.numeric(add_bev | input$treatment_arm %in% c("stupp_bev","bev_salvage")),
        PD1_ON    = as.numeric(add_pd1 | input$treatment_arm == "pembro"),
        TTF_ON    = as.numeric(add_ttf | input$treatment_arm == "stupp_ttf")
      )

      conc_days_use <- if(flags$CONCURRENT_RT > 0) input$conc_days else 0
      bev_n <- if(flags$BEV_ON > 0) 12 else 0
      pd1_n <- if(flags$PD1_ON > 0) 12 else 0
      adj_n <- if(flags$ADJUVANT_TMZ > 0) input$adj_cycles else 0

      evts <- build_ev_shiny(
        tmz_dose=input$conc_dose, adj_dose=input$adj_dose,
        adj_cycles=adj_n, adj_start=conc_days_use+15,
        bev_dose=if(flags$BEV_ON > 0) 10 else 0,
        bev_start=0, bev_n=bev_n,
        pd1_dose=if(flags$PD1_ON > 0) 200 else 0,
        pd1_start=0, pd1_n=pd1_n,
        conc_days=conc_days_use, BSA=BSA, WT=WT
      )

      p_list <- modifyList(flags, list(
        MGMT_meth=mgmt, IDH_status=idh, BSA=BSA, WT=WT,
        Ts0_init=input$init_cells * 1e8
      ))

      mod_run <- mod_shiny %>% param(p_list)

      if (nrow(evts) > 0) {
        out <- mrgsim(mod_run, events=evts, end=days*24, delta=6, obsonly=TRUE)
      } else {
        out <- mrgsim(mod_run, end=days*24, delta=6, obsonly=TRUE)
      }

      df <- as.data.frame(out)
      df$time_day <- df$time / 24
      df
    })
  }, ignoreNULL=FALSE)

  # All-scenario sim
  all_scenarios_sim <- eventReactive(input$run_sim, {
    withProgress(message="Running all scenarios...", {
      mgmt <- as.numeric(input$MGMT_status)
      idh  <- as.numeric(input$IDH_status)
      BSA  <- input$BSA_val; WT <- input$WT_val
      days <- input$sim_days

      scen_list <- list(
        "Untreated"            = list(flags=list(CONCURRENT_RT=0,ADJUVANT_TMZ=0,BEV_ON=0,PD1_ON=0,TTF_ON=0),
                                       ev=build_ev_shiny(conc_days=0,adj_cycles=0,bev_n=0,pd1_n=0,BSA=BSA,WT=WT)),
        "Stupp (MGMT+)"        = list(flags=list(CONCURRENT_RT=1,ADJUVANT_TMZ=1,BEV_ON=0,PD1_ON=0,TTF_ON=0),
                                       ev=build_ev_shiny(conc_days=42,adj_cycles=6,bev_n=0,pd1_n=0,BSA=BSA,WT=WT)),
        "Stupp (MGMT-)"        = list(flags=list(CONCURRENT_RT=1,ADJUVANT_TMZ=1,BEV_ON=0,PD1_ON=0,TTF_ON=0),
                                       ev=build_ev_shiny(conc_days=42,adj_cycles=6,bev_n=0,pd1_n=0,BSA=BSA,WT=WT)),
        "Stupp+BEV"            = list(flags=list(CONCURRENT_RT=1,ADJUVANT_TMZ=1,BEV_ON=1,PD1_ON=0,TTF_ON=0),
                                       ev=build_ev_shiny(conc_days=42,adj_cycles=6,bev_n=12,bev_start=0,pd1_n=0,BSA=BSA,WT=WT)),
        "Stupp+TTF"            = list(flags=list(CONCURRENT_RT=1,ADJUVANT_TMZ=1,BEV_ON=0,PD1_ON=0,TTF_ON=1),
                                       ev=build_ev_shiny(conc_days=42,adj_cycles=6,bev_n=0,pd1_n=0,BSA=BSA,WT=WT)),
        "Pembrolizumab+TMZ"    = list(flags=list(CONCURRENT_RT=0,ADJUVANT_TMZ=1,BEV_ON=0,PD1_ON=1,TTF_ON=0),
                                       ev=build_ev_shiny(conc_days=0,adj_cycles=8,adj_start=0,bev_n=0,
                                                         pd1_n=12,pd1_start=0,pd1_dose=200,BSA=BSA,WT=WT)),
        "BEV Salvage"          = list(flags=list(CONCURRENT_RT=0,ADJUVANT_TMZ=0,BEV_ON=1,PD1_ON=0,TTF_ON=0),
                                       ev=build_ev_shiny(conc_days=0,adj_cycles=0,bev_n=12,bev_start=0,pd1_n=0,BSA=BSA,WT=WT))
      )

      mgmt_vals <- c(1,1,0,1,1,1,1)

      res_all <- lapply(seq_along(scen_list), function(i) {
        nm   <- names(scen_list)[i]
        sc   <- scen_list[[i]]
        mgmt_i <- mgmt_vals[i]
        p_run <- modifyList(sc$flags, list(MGMT_meth=mgmt_i, IDH_status=idh,
                                            BSA=BSA, WT=WT,
                                            Ts0_init=input$init_cells*1e8))
        m_run <- mod_shiny %>% param(p_run)
        if (nrow(sc$ev) > 0) {
          out <- mrgsim(m_run, events=sc$ev, end=days*24, delta=6, obsonly=TRUE)
        } else {
          out <- mrgsim(m_run, end=days*24, delta=6, obsonly=TRUE)
        }
        df <- as.data.frame(out)
        df$time_day <- df$time / 24
        df$scenario <- nm
        df
      })
      bind_rows(res_all)
    })
  }, ignoreNULL=FALSE)

  # ---- Outputs ----
  output$pk_tmz_plot <- renderPlot({
    df <- sim_result()
    dfpk <- df %>% filter(time_day <= min(10, max(df$time_day)))
    ggplot(dfpk, aes(x=time_day)) +
      geom_line(aes(y=Cp_tmz, color="Plasma"), linewidth=1) +
      geom_line(aes(y=Cbrain, color="Brain ECF"), linewidth=1) +
      scale_color_manual(values=c("Plasma"="#2196F3","Brain ECF"="#FF9800")) +
      labs(title="TMZ PK: Plasma vs. Brain ECF (first 10 days)",
           x="Time (days)", y="Concentration (mg/L = µg/mL)", color="Compartment") +
      theme_bw(base_size=12) + theme(legend.position="bottom")
  })

  output$pk_bev_plot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x=time_day)) +
      geom_line(aes(y=BEV_Cp, color="BEV Plasma"), linewidth=1) +
      geom_line(aes(y=VEGF_cap, color="Free VEGF"), linewidth=1, linetype=2) +
      scale_color_manual(values=c("BEV Plasma"="#4FC3F7","Free VEGF"="#E91E63")) +
      labs(title="Bevacizumab PK & Free VEGF",
           x="Time (days)", y="Concentration (mg/L or pg/mL)", color="") +
      theme_bw(base_size=12) + theme(legend.position="bottom")
  })

  output$pk_pd1_plot <- renderPlot({
    df <- sim_result()
    p1 <- ggplot(df, aes(x=time_day, y=APD1_Cp)) +
      geom_line(color="#8BC34A", linewidth=1) +
      labs(title="Anti-PD1 (Pembrolizumab) Plasma",
           x="Time (days)", y="Conc. (mg/L)") + theme_bw(base_size=12)
    p2 <- ggplot(df, aes(x=time_day, y=PD1_occ_cap)) +
      geom_line(color="#FF5722", linewidth=1) + ylim(0,1) +
      labs(title="PD-1 Receptor Occupancy",
           x="Time (days)", y="Occupancy (fraction)") + theme_bw(base_size=12)
    gridExtra::grid.arrange(p1, p2, nrow=1)
  })

  output$dna_lesion_plot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x=time_day, y=O6MeG_cap)) +
      geom_line(color="#E91E63", linewidth=1.2) +
      labs(title="O6-MeG DNA Lesion Over Time",
           subtitle="MGMT repair reduces lesion accumulation (methylated = less repair)",
           x="Time (days)", y="O6-MeG Lesion (rel. units)") +
      theme_bw(base_size=12)
  })

  output$mgmt_compare_plot <- renderPlot({
    BSA <- input$BSA_val; WT <- input$WT_val; days <- input$sim_days
    idh <- as.numeric(input$IDH_status)
    ev  <- build_ev_shiny(conc_days=42,adj_cycles=6,bev_n=0,pd1_n=0,BSA=BSA,WT=WT)
    results <- lapply(c("Methylated"=1,"Unmethylated"=0), function(mg) {
      m_run <- mod_shiny %>% param(list(MGMT_meth=mg, IDH_status=idh, BSA=BSA, WT=WT,
                                         CONCURRENT_RT=1, ADJUVANT_TMZ=1,
                                         Ts0_init=input$init_cells*1e8))
      out <- mrgsim(m_run, events=ev, end=days*24, delta=6, obsonly=TRUE)
      df <- as.data.frame(out); df$time_day <- df$time/24; df
    })
    bind_rows(Map(function(df, nm) { df$MGMT <- nm; df },
                  results, names(results))) %>%
      ggplot(aes(x=time_day, y=O6MeG_cap, color=MGMT)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=c("Methylated"="#4CAF50","Unmethylated"="#E91E63")) +
      labs(title="O6-MeG Lesion: MGMT Methylation Comparison (Stupp Protocol)",
           x="Time (days)", y="O6-MeG Lesion (rel. units)", color="MGMT Status") +
      theme_bw(base_size=12)
  })

  output$tumor_cells_plot <- renderPlot({
    df <- sim_result() %>%
      select(time_day, Ts, Tr, GSC) %>%
      pivot_longer(-time_day, names_to="Population", values_to="Cells")
    ggplot(df, aes(x=time_day, y=Cells, color=Population)) +
      geom_line(linewidth=1) +
      scale_y_log10() +
      scale_color_manual(values=c("Ts"="#2196F3","Tr"="#E91E63","GSC"="#FF9800"),
                         labels=c("Ts"="Sensitive Cells","Tr"="Resistant Cells",
                                  "GSC"="Glioma Stem Cells")) +
      labs(title="Tumor Cell Population Dynamics (log scale)",
           x="Time (days)", y="Cell Count (log scale)", color="Population") +
      theme_bw(base_size=12) + theme(legend.position="bottom")
  })

  output$tumor_vol_plot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x=time_day, y=TumorVol_mL)) +
      geom_line(color="#3949AB", linewidth=1.2) +
      scale_y_log10() +
      geom_hline(yintercept=1, linetype="dashed", color="gray") +
      annotate("text", x=max(df$time_day)*0.05, y=1.5, label="1 mL ≈ 10⁹ cells",
               size=3.5, color="gray40") +
      labs(title="Tumor Volume (mL, log scale)",
           subtitle="Derived from total cell count (Ts + Tr + GSC)",
           x="Time (days)", y="Tumor Volume (mL)") +
      theme_bw(base_size=12)
  })

  output$scenario_plot <- renderPlot({
    df_all <- all_scenarios_sim()
    pal <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF")
    ggplot(df_all, aes(x=time_day, y=TumorVol_mL, color=scenario)) +
      geom_line(linewidth=0.9) +
      scale_color_manual(values=pal) +
      scale_y_log10() +
      labs(title="All-Scenario Tumor Volume Comparison",
           x="Time (days)", y="Tumor Volume (mL, log)", color="Treatment") +
      theme_bw(base_size=12) + theme(legend.position="bottom",
                                      legend.text=element_text(size=9)) +
      guides(color=guide_legend(nrow=3))
  })

  output$scenario_table <- renderDT({
    df_all <- all_scenarios_sim()
    df_all %>%
      group_by(scenario) %>%
      summarise(
        `Vol at 3mo (mL)` = round(TumorVol_mL[which.min(abs(time_day-90))], 4),
        `Vol at 6mo (mL)` = round(TumorVol_mL[which.min(abs(time_day-180))], 4),
        `Vol at 1yr (mL)` = round(TumorVol_mL[which.min(abs(time_day-365))], 4),
        `Min Vol (mL)`    = round(min(TumorVol_mL), 5),
        `Final Vol (mL)`  = round(tail(TumorVol_mL, 1), 4)
      ) %>%
      datatable(options=list(pageLength=10, dom='t'), rownames=FALSE) %>%
      formatStyle("scenario", fontWeight="bold")
  })

  output$tme_plot <- renderPlot({
    df <- sim_result() %>%
      select(time_day, CD8_cap, Treg_cap, TAM_cap) %>%
      pivot_longer(-time_day, names_to="Cell", values_to="Level")
    ggplot(df, aes(x=time_day, y=Level, color=Cell)) +
      geom_line(linewidth=1) +
      scale_color_manual(values=c("CD8_cap"="#4CAF50","Treg_cap"="#E91E63",
                                   "TAM_cap"="#9C27B0"),
                         labels=c("CD8_cap"="CD8+ CTL","Treg_cap"="Tregs",
                                  "TAM_cap"="M2-TAMs")) +
      labs(title="Tumor Immune Microenvironment Dynamics",
           x="Time (days)", y="Relative Cell Level", color="Cell Type") +
      theme_bw(base_size=12) + theme(legend.position="bottom")
  })

  output$vegf_nv_plot <- renderPlot({
    df <- sim_result()
    p1 <- ggplot(df, aes(x=time_day, y=VEGF_cap)) +
      geom_line(color="#2196F3", linewidth=1.2) +
      labs(title="Free VEGF-A", x="Time (days)", y="VEGF (pg/mL)") +
      theme_bw(base_size=11)
    p2 <- ggplot(df, aes(x=time_day, y=NV_cap)) +
      geom_line(color="#FF5722", linewidth=1.2) +
      labs(title="Neovascularization Index", x="Time (days)", y="NV (rel. units)") +
      theme_bw(base_size=11)
    gridExtra::grid.arrange(p1, p2, nrow=1)
  })

  output$pd1_occ_plot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x=time_day, y=PD1_occ_cap)) +
      geom_line(color="#9C27B0", linewidth=1.2) + ylim(0,1) +
      labs(title="PD-1 Receptor Occupancy",
           subtitle="Target: >90% at trough for pembrolizumab",
           x="Time (days)", y="PD-1 Occupancy (fraction)") +
      geom_hline(yintercept=0.9, linetype="dashed", color="red", alpha=0.7) +
      annotate("text", x=max(df$time_day)*0.5, y=0.93,
               label="90% target", color="red", size=3.5) +
      theme_bw(base_size=12)
  })

  output$diam_plot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x=time_day, y=TumorDiam_cm)) +
      geom_line(color="#D32F2F", linewidth=1.2) +
      geom_hline(yintercept=3.0, linetype="dashed", color="gray40") +
      annotate("text", x=max(df$time_day)*0.1, y=3.2,
               label="~3 cm threshold (surgical consideration)", size=3) +
      labs(title="Tumor Diameter Over Time (RANO Surrogate)",
           subtitle="Approximate sphere from tumor volume",
           x="Time (days)", y="Tumor Diameter (cm)") +
      theme_bw(base_size=12)
  })

  output$kill_components_plot <- renderPlot({
    df <- sim_result() %>%
      select(time_day, kill_TMZ_cap, kill_CD8_cap) %>%
      pivot_longer(-time_day, names_to="Kill_Mech", values_to="Rate")
    ggplot(df, aes(x=time_day, y=Rate, fill=Kill_Mech)) +
      geom_area(alpha=0.7, position="stack") +
      scale_fill_manual(values=c("kill_TMZ_cap"="#E91E63","kill_CD8_cap"="#4CAF50"),
                        labels=c("kill_CD8_cap"="CD8 T cell kill","kill_TMZ_cap"="TMZ kill")) +
      labs(title="Tumor Kill Rate Components",
           x="Time (days)", y="Kill Rate (day⁻¹)", fill="Mechanism") +
      theme_bw(base_size=12) + theme(legend.position="bottom")
  })

  output$pt_summary_table <- renderDT({
    data.frame(
      Parameter = c("Age","KPS","Extent of Resection","MGMT Status","IDH Status",
                    "Initial Tumor Size","BSA","Body Weight",
                    "Treatment Arm","Sim. Duration"),
      Value     = c(input$age_val, input$kps_val, input$eor_val,
                    ifelse(input$MGMT_status==1,"Methylated","Unmethylated"),
                    ifelse(input$IDH_status==0,"Wild-type","Mutant"),
                    paste0(input$init_cells," × 10⁸ cells"),
                    paste0(input$BSA_val," m²"), paste0(input$WT_val," kg"),
                    input$treatment_arm, paste0(input$sim_days," days"))
    ) %>%
      datatable(options=list(dom='t', pageLength=15), rownames=FALSE,
                colnames=c("Parameter","Value"))
  })
}

# ============================================================
# RUN APP
# ============================================================
if (!requireNamespace("gridExtra", quietly=TRUE)) install.packages("gridExtra")
library(gridExtra)

shinyApp(ui=ui, server=server)
