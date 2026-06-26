## ============================================================
## Sickle Cell Disease (SCD) — Interactive QSP Shiny App
## ============================================================
## Tabs:
##   1. Patient Profile & Disease Overview
##   2. Drug PK Profiles
##   3. Hematological Response (Hgb, HbF, Reticulocytes)
##   4. Vaso-occlusion & Biomarkers
##   5. Treatment Scenario Comparison
##   6. End-Organ & Biomarker Dashboard
## ============================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(shinydashboard)
library(tidyr)
library(scales)

## ─── Inline model (same as scd_mrgsolve_model.R) ─────────────
scd_model_code <- '
$PARAM
  ka_HU=2.0, CL_HU=3.5, Vd_HU=28.0, F_HU=0.80, MW_HU=76.1,
  ka_VOX=0.45, CL_VOX=2.1, Vd_VOX=97.0,
  kon_VOX=0.15, koff_VOX=0.005, RBC_capacity=50,
  CL_CRIZ=0.008, Vc_CRIZ=3.5, Vp_CRIZ=2.8,
  k12_CRIZ=0.015, k21_CRIZ=0.012,
  ka_LG=1.2, CL_LG=45.0, Vd_LG=18.0,
  kprod_CFU=0.015, kdiff_CFU=0.08, kdiff_RET=0.035,
  kdeath_S=0.0055, kdeath_N=0.00014,
  RET_0=45, RBC_S_0=1800, RBC_N_0=200, CFU_E_0=12,
  HbF0=0.07, kHbF_deg=0.002,
  EC50_HU_HbF=15.0, Emax_HU_HbF=0.18,
  Hgb_0=8.5, kprod_Hgb=0.003, kdeg_Hgb=0.00035,
  LDH_0=480, kprod_LDH=0.15, kdeg_LDH=0.08,
  Bili_0=42, kdeg_Bili=0.02,
  Hp_0=0.25, kdeg_Hp=0.03, kprod_Hp=0.008,
  kfree_Hb=0.12,
  NO_0=0.40, kNO_prod=0.05, kNO_scav=0.8,
  Psel_0=1.0, kPsel_deg=0.04,
  EC50_CRIZ_Psel=0.5, Imax_CRIZ_Psel=0.80,
  VOC_0=0.0004, kVOC_Psel=0.8, kVOC_NO=0.6, kVOC_HbF=0.5,
  NADH_0=0.55, kNADH_prod=0.04, kNADH_ox=0.06,
  Iron_0=900, kIron_acc=0.002, kIron_deg=0.001,
  TRV_0=2.65, kTRV_prod=0.0001, kTRV_deg=0.0002,
  eGFR_0=105, keGFR_dec=0.000005,
  DOSE_HU=0, BW=65, DOSE_VOX=0, DOSE_CRIZ=0, DOSE_LG=0

$CMT HU_gut HU_plasma VOX_plasma VOX_RBC CRIZ_C CRIZ_P LG
     CFU_E RET RBC_S RBC_N HbF_frac Hgb free_Hb
     Haptoglobin LDH Bilirubin NO P_selectin VOC NADH Iron TRV eGFR

$INIT
  HU_gut=0, HU_plasma=0, VOX_plasma=0, VOX_RBC=0,
  CRIZ_C=0, CRIZ_P=0, LG=0,
  CFU_E=12, RET=45, RBC_S=1800, RBC_N=200,
  HbF_frac=0.07, Hgb=8.5, free_Hb=8.0,
  Haptoglobin=0.25, LDH=480, Bilirubin=42,
  NO=0.40, P_selectin=1.0, VOC=0.0004,
  NADH=0.55, Iron=900, TRV=2.65, eGFR=105

$ODE
  double HU_mg=DOSE_HU*BW; double HU_umol=HU_mg*1000/MW_HU;
  double HU_rate=HU_umol/24.0;
  double VOX_rate=DOSE_VOX*1500.0*1000/392.4/24.0;
  double LG_rate=DOSE_LG*10000.0*1000/146.1/24.0;
  dxdt_HU_gut=F_HU*HU_rate-ka_HU*HU_gut;
  dxdt_HU_plasma=ka_HU*HU_gut/Vd_HU-(CL_HU/Vd_HU)*HU_plasma;
  double HU_uM=HU_plasma;
  double VOX_free=VOX_plasma;
  dxdt_VOX_plasma=VOX_rate/Vd_VOX-(CL_VOX/Vd_VOX)*VOX_plasma
    -kon_VOX*VOX_free*(RBC_capacity-VOX_RBC)+koff_VOX*VOX_RBC;
  dxdt_VOX_RBC=kon_VOX*VOX_free*(RBC_capacity-VOX_RBC)-koff_VOX*VOX_RBC;
  dxdt_CRIZ_C=-(CL_CRIZ/Vc_CRIZ)*CRIZ_C-k12_CRIZ*CRIZ_C+k21_CRIZ*CRIZ_P;
  dxdt_CRIZ_P=k12_CRIZ*CRIZ_C-k21_CRIZ*CRIZ_P;
  double CRIZ_ugml=CRIZ_C*148000/1e6;
  dxdt_LG=LG_rate/Vd_LG-(CL_LG/Vd_LG)*LG;
  double E_HU_HbF=Emax_HU_HbF*pow(HU_uM,1.5)/(pow(EC50_HU_HbF,1.5)+pow(HU_uM,1.5));
  double VOX_occ=VOX_RBC/(RBC_capacity+1e-6);
  double VOX_polymerInhib=0.45*VOX_occ;
  double I_CRIZ_Psel=Imax_CRIZ_Psel*CRIZ_ugml/(EC50_CRIZ_Psel+CRIZ_ugml);
  double E_LG_NADH=0.15*DOSE_LG*(LG/(LG+100));
  dxdt_HbF_frac=(HbF0+E_HU_HbF)*kHbF_deg-kHbF_deg*HbF_frac;
  double EPO_stim=pow(Hgb_0/(Hgb+0.01),1.5);
  dxdt_CFU_E=kprod_CFU*EPO_stim-kdiff_CFU*CFU_E;
  dxdt_RET=kdiff_CFU*CFU_E-kdiff_RET*RET;
  double eff_sickle=(1-HbF_frac)*(1-VOX_polymerInhib);
  dxdt_RBC_S=kdiff_RET*RET*eff_sickle-kdeath_S*RBC_S;
  dxdt_RBC_N=kdiff_RET*RET*(1-eff_sickle)-kdeath_N*RBC_N;
  double RBC_total=RBC_S+RBC_N;
  dxdt_Hgb=kprod_Hgb*(RBC_total/(RBC_S_0+RBC_N_0))-kdeg_Hgb*Hgb;
  double hemolysis_rate=kdeath_S*RBC_S*0.015;
  dxdt_free_Hb=hemolysis_rate-kfree_Hb*free_Hb-kdeg_Hp*Haptoglobin*free_Hb;
  dxdt_Haptoglobin=kprod_Hp-kdeg_Hp*free_Hb*Haptoglobin;
  double LDH_base_rate=kprod_LDH*(hemolysis_rate/(kdeath_S*RBC_S_0*0.015+1e-6))*LDH_0;
  dxdt_LDH=LDH_base_rate-kdeg_LDH*LDH;
  dxdt_Bilirubin=0.006*hemolysis_rate*Bili_0-kdeg_Bili*Bilirubin;
  dxdt_NO=kNO_prod-kNO_scav*free_Hb*NO-0.01*NO;
  double Psel_stim=1.0+0.3*(1-NO/NO_0);
  double Psel_inhib=1-I_CRIZ_Psel;
  dxdt_P_selectin=kPsel_deg*Psel_stim*Psel_inhib-kPsel_deg*P_selectin;
  double VOC_driver=P_selectin*kVOC_Psel*(1.0/(NO+0.1))*kVOC_NO*(1-HbF_frac*kVOC_HbF);
  double VOC_norm=Psel_0*kVOC_Psel*(1.0/(NO_0+0.1))*kVOC_NO*(1-HbF0*kVOC_HbF);
  dxdt_VOC=VOC_0*(VOC_driver/(VOC_norm+1e-9)+1e-9)-0.05*VOC;
  dxdt_NADH=kNADH_prod*(1+E_LG_NADH)-kNADH_ox*NADH;
  dxdt_Iron=kIron_acc*hemolysis_rate*100-kIron_deg*Iron;
  double TRV_drive=1+0.002*(1-NO/NO_0)+0.001*(LDH/LDH_0-1);
  dxdt_TRV=kTRV_prod*TRV_drive-kTRV_deg*TRV;
  dxdt_eGFR=-keGFR_dec*(1+0.5*(RBC_S/RBC_S_0-1))*eGFR;

$TABLE
  double VOC_annual=VOC*8760;
  double HbF_pct=HbF_frac*100;
  double ret_pct=RET/(RBC_S+RBC_N+RET)*100;
  double CRIZ_conc=CRIZ_C*148000/1e6;
  capture(VOC_annual, HbF_pct, ret_pct, CRIZ_conc,
          LDH, Bilirubin, NO, NADH, Iron, TRV, eGFR,
          Haptoglobin, HU_plasma, VOX_RBC, CRIZ_C)
'

mod_base <- mcode("SCD_Shiny", scd_model_code, quiet=TRUE)

## ─── Helper: run simulation ───────────────────────────────────
run_sim <- function(mod, hu_dose, vox_on, criz_on, lg_on, bw, sim_days, criz_dose) {
  mod_p <- param(mod, DOSE_HU=hu_dose, DOSE_VOX=as.numeric(vox_on),
                 DOSE_CRIZ=as.numeric(criz_on), DOSE_LG=as.numeric(lg_on), BW=bw)
  sim_times <- seq(0, sim_days*24, by=24)

  if (criz_on) {
    criz_dose_nmol <- (criz_dose * bw * 1e6) / 148000
    criz_ev_times  <- seq(0, sim_days*24, by=28*24)
    ev_criz <- ev(amt=criz_dose_nmol, cmt="CRIZ_C", time=criz_ev_times)
    out <- mrgsim(mod_p, events=ev_criz, times=sim_times, delta=1, end=sim_days*24)
  } else {
    out <- mrgsim(mod_p, times=sim_times, delta=1, end=sim_days*24)
  }
  as.data.frame(out) %>% mutate(time_days = time / 24)
}

## ─── Theme ────────────────────────────────────────────────────
theme_scd <- function() {
  theme_minimal(base_size=13) +
    theme(
      plot.title    = element_text(face="bold", size=14, colour="#2C3E50"),
      plot.subtitle = element_text(size=11, colour="#7F8C8D"),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom",
      legend.key.size  = unit(0.5,"cm")
    )
}

scd_red   <- "#C0392B"
scd_blue  <- "#2980B9"
scd_green <- "#27AE60"

## ═══════════════════════════════════════════════════════════════
## UI
## ═══════════════════════════════════════════════════════════════
ui <- dashboardPage(
  skin="red",
  dashboardHeader(title=tags$span(
    tags$img(src="", height="30px"),
    "Sickle Cell Disease QSP Dashboard"), titleWidth=400),

  dashboardSidebar(
    width=320,
    sidebarMenu(
      id="tabs",
      menuItem("Patient Profile",         tabName="tab_profile",   icon=icon("user-circle")),
      menuItem("Drug PK Profiles",        tabName="tab_pk",        icon=icon("pills")),
      menuItem("Hematological Response",  tabName="tab_heme",      icon=icon("tint")),
      menuItem("Vaso-occlusion & Biomarkers", tabName="tab_voc",   icon=icon("heartbeat")),
      menuItem("Scenario Comparison",     tabName="tab_compare",   icon=icon("chart-bar")),
      menuItem("End-Organ Dashboard",     tabName="tab_organ",     icon=icon("lungs"))
    ),
    hr(),
    h4("Patient Parameters", style="padding-left:15px; color:#ECF0F1"),
    sliderInput("bw",     "Body Weight (kg)",   min=30, max=120, value=65, step=1),
    sliderInput("sim_days","Simulation Duration (days)", min=30, max=730, value=365, step=30),
    hr(),
    h4("Treatment Selection", style="padding-left:15px; color:#ECF0F1"),
    checkboxInput("use_hu",   "Hydroxyurea (HU)",        value=FALSE),
    conditionalPanel("input.use_hu",
      sliderInput("hu_dose","HU Dose (mg/kg/d)", min=5, max=35, value=20, step=5)
    ),
    checkboxInput("use_vox",  "Voxelotor (1500 mg/d)",   value=FALSE),
    checkboxInput("use_criz", "Crizanlizumab (IV q4w)",  value=FALSE),
    conditionalPanel("input.use_criz",
      sliderInput("criz_dose","CRIZ Dose (mg/kg)", min=2.5, max=7.5, value=5.0, step=2.5)
    ),
    checkboxInput("use_lg",   "L-Glutamine (5g BID)",    value=FALSE),
    hr(),
    actionButton("run_sim", "Run Simulation", icon=icon("play"),
                 class="btn-primary btn-block")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #F8F9FA; }
      .box { border-radius:8px; }
      .info-box-icon { border-radius:6px 0 0 6px; }
      .valueBox { border-radius:6px; }
      .skin-red .main-header .logo { background-color:#7B241C; }
      .skin-red .main-header .navbar { background-color:#922B21; }
    "))),

    tabItems(

      ## ─── TAB 1: Patient Profile ────────────────────────────
      tabItem("tab_profile",
        fluidRow(
          box(width=12, title="Sickle Cell Disease — Disease Overview",
              status="danger", solidHeader=TRUE,
              tabsetPanel(
                tabPanel("Pathophysiology",
                  h4("Molecular Basis of SCD"),
                  p("SCD is caused by a single nucleotide substitution (GAG→GTG) in codon 6 of the ",
                    strong("HBB gene"), ", leading to the Glu6Val mutation (βS-globin). The resulting ",
                    strong("HbS tetramer (α₂βS₂)"), " undergoes deoxygenation-induced polymerization,",
                    " distorting RBCs into the characteristic sickle shape."),
                  tags$ul(
                    tags$li(strong("Polymerization:"), " Deoxy-HbS forms 14-strand polymer fibers at a critical concentration (MCHC >25 g/dL)."),
                    tags$li(strong("Sickling:"), " Repeated cycles of oxygenation/deoxygenation cause cumulative membrane damage, RBC dehydration, and reduced lifespan (10–20 days vs. 120 days)."),
                    tags$li(strong("Hemolysis:"), " Intravascular and extravascular hemolysis produces cell-free Hb, free heme, LDH, and bilirubin."),
                    tags$li(strong("Vaso-occlusion:"), " P-selectin-mediated adhesion of sickle RBCs, neutrophils, and platelets to activated endothelium triggers microvascular stasis."),
                    tags$li(strong("NO scavenging:"), " Cell-free Hb reacts with NO (k=6.4×10⁷ M⁻¹s⁻¹), reducing vasodilation and promoting vasoconstriction.")
                  )
                ),
                tabPanel("QSP Model Structure",
                  h4("Model Architecture: 24 ODEs"),
                  tableOutput("model_structure"),
                  p("The model integrates drug PK (HU, VOX, CRIZ, L-Glu), erythropoiesis, ",
                    "hemolysis, vascular biology, and end-organ dysfunction.")
                ),
                tabPanel("Treatment Targets",
                  h4("Approved Treatments & Mechanisms"),
                  tableOutput("treatment_moa")
                )
              )
          )
        ),
        fluidRow(
          valueBoxOutput("box_genotype"), valueBoxOutput("box_hgb_base"),
          valueBoxOutput("box_hbf_base"), valueBoxOutput("box_voc_base")
        )
      ),

      ## ─── TAB 2: Drug PK ────────────────────────────────────
      tabItem("tab_pk",
        fluidRow(
          box(width=12, title="Pharmacokinetic Profiles", status="info", solidHeader=TRUE,
              plotlyOutput("pk_plot", height="500px"))
        ),
        fluidRow(
          box(width=6, title="PK Parameters Summary", status="primary", solidHeader=TRUE,
              tableOutput("pk_params")),
          box(width=6, title="PK Notes", status="primary",
              h5("Hydroxyurea"),
              tags$ul(
                tags$li("Oral bioavailability: ~80–100%; Cmax: 25–100 μM at 15–35 mg/kg/d"),
                tags$li("t½: ~2–4 h; 80% excreted unchanged in urine"),
                tags$li("RNR inhibition → stress erythropoiesis → γ-globin induction → HbF↑")
              ),
              h5("Voxelotor (GBT440)"),
              tags$ul(
                tags$li("1500 mg/d oral; Cmax ~2 μM; t½: 16–25 h"),
                tags$li("Binds α-globin T-state, stabilizes HbS O₂ affinity"),
                tags$li("RBC-bound fraction >99%; reduces deoxy-HbS fraction")
              ),
              h5("Crizanlizumab"),
              tags$ul(
                tags$li("5 mg/kg IV q4w (IgG2 anti-P-selectin mAb)"),
                tags$li("t½: ~10 days; Vc=3.5 L; linear clearance"),
                tags$li("Target occupancy >80% achieved at 5 mg/kg")
              )
          )
        )
      ),

      ## ─── TAB 3: Hematological Response ─────────────────────
      tabItem("tab_heme",
        fluidRow(
          column(6, box(width=NULL, title="Hemoglobin (g/dL)", status="danger", solidHeader=TRUE,
              plotlyOutput("hgb_plot", height="320px"))),
          column(6, box(width=NULL, title="HbF Percentage (%)", status="success", solidHeader=TRUE,
              plotlyOutput("hbf_plot", height="320px")))
        ),
        fluidRow(
          column(6, box(width=NULL, title="Reticulocyte Count (%)", status="warning", solidHeader=TRUE,
              plotlyOutput("ret_plot", height="320px"))),
          column(6, box(width=NULL, title="LDH (Hemolysis Marker, U/L)", status="danger", solidHeader=TRUE,
              plotlyOutput("ldh_plot", height="320px")))
        ),
        fluidRow(
          box(width=12, title="Hematology Summary at Selected Timepoints", status="primary", solidHeader=TRUE,
              DTOutput("heme_table"))
        )
      ),

      ## ─── TAB 4: Vaso-occlusion & Biomarkers ─────────────────
      tabItem("tab_voc",
        fluidRow(
          column(6, box(width=NULL, title="Annual VOC Rate (crises/year)", status="danger", solidHeader=TRUE,
              plotlyOutput("voc_plot", height="320px"))),
          column(6, box(width=NULL, title="P-Selectin Expression (Relative)", status="warning", solidHeader=TRUE,
              plotlyOutput("psel_plot", height="320px")))
        ),
        fluidRow(
          column(6, box(width=NULL, title="Nitric Oxide Bioavailability", status="success", solidHeader=TRUE,
              plotlyOutput("no_plot", height="320px"))),
          column(6, box(width=NULL, title="Bilirubin (μmol/L)", status="warning", solidHeader=TRUE,
              plotlyOutput("bili_plot", height="320px")))
        ),
        fluidRow(
          box(width=12, title="Vaso-occlusion Biomarker Summary", status="info", solidHeader=TRUE,
              DTOutput("voc_table"))
        )
      ),

      ## ─── TAB 5: Scenario Comparison ────────────────────────
      tabItem("tab_compare",
        fluidRow(
          box(width=12, title="Multi-Scenario Treatment Comparison",
              status="primary", solidHeader=TRUE,
              p("Select up to 4 scenarios to compare side-by-side."),
              fluidRow(
                column(3, checkboxGroupInput("scen_check", "Scenarios to Compare:",
                  choices = c(
                    "No Treatment"          = "no_tx",
                    "Hydroxyurea 20 mg/kg"  = "hu20",
                    "Voxelotor"             = "vox",
                    "Crizanlizumab"         = "criz",
                    "L-Glutamine"           = "lg",
                    "HU + Voxelotor"        = "hu_vox",
                    "HU + VOX + CRIZ"       = "triple"
                  ),
                  selected = c("no_tx","hu20","vox","hu_vox")
                )),
                column(9,
                  selectInput("compare_endpoint", "Primary Endpoint:",
                    choices=c("Hemoglobin (g/dL)"="Hgb",
                              "HbF (%)"="HbF_pct",
                              "VOC Rate (crises/yr)"="VOC_annual",
                              "LDH (U/L)"="LDH",
                              "NO Index"="NO",
                              "P-Selectin"="P_selectin",
                              "TRV (m/s)"="TRV"),
                    selected="Hgb"),
                  plotlyOutput("compare_plot", height="380px")
                )
              )
          )
        ),
        fluidRow(
          box(width=12, title="Comparative Efficacy Table (1-Year Outcomes)",
              status="success", solidHeader=TRUE,
              DTOutput("compare_table"))
        )
      ),

      ## ─── TAB 6: End-Organ Dashboard ─────────────────────────
      tabItem("tab_organ",
        fluidRow(
          valueBoxOutput("kpi_hgb"),  valueBoxOutput("kpi_hbf"),
          valueBoxOutput("kpi_voc"),  valueBoxOutput("kpi_ldh")
        ),
        fluidRow(
          column(6, box(width=NULL, title="TRV — Pulmonary HTN Marker (m/s)",
              status="danger", solidHeader=TRUE,
              plotlyOutput("trv_plot", height="280px"))),
          column(6, box(width=NULL, title="eGFR — Renal Function (mL/min/1.73m²)",
              status="warning", solidHeader=TRUE,
              plotlyOutput("gfr_plot", height="280px")))
        ),
        fluidRow(
          column(6, box(width=NULL, title="Ferritin / Iron Stores (μg/L)",
              status="warning", solidHeader=TRUE,
              plotlyOutput("iron_plot", height="280px"))),
          column(6, box(width=NULL, title="Haptoglobin (g/L)",
              status="info", solidHeader=TRUE,
              plotlyOutput("hp_plot", height="280px")))
        ),
        fluidRow(
          box(width=12, title="End-Organ KPI Tracker", status="primary", solidHeader=TRUE,
              DTOutput("organ_table"))
        )
      )
    )
  )
)

## ═══════════════════════════════════════════════════════════════
## SERVER
## ═══════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  ## ── Simulation reactive ──────────────────────────────────────
  sim_data <- eventReactive(input$run_sim, ignoreNULL=FALSE, {
    withProgress(message="Running QSP simulation...", value=0.5, {
      run_sim(mod_base,
              hu_dose   = if(input$use_hu) input$hu_dose else 0,
              vox_on    = input$use_vox,
              criz_on   = input$use_criz,
              lg_on     = input$use_lg,
              bw        = input$bw,
              sim_days  = input$sim_days,
              criz_dose = if(input$use_criz) input$criz_dose else 0)
    })
  })

  ## ── Model structure table ────────────────────────────────────
  output$model_structure <- renderTable({
    data.frame(
      Module = c("Drug PK","Drug PK","Drug PK","Drug PK",
                 "Erythropoiesis","RBC Dynamics","Hemoglobin",
                 "Hemolysis","Vascular","Vaso-occlusion",
                 "Oxidative Stress","Iron/Organ"),
      `ODEs` = c("HU gut, HU plasma","VOX plasma, VOX-RBC",
                 "CRIZ central, CRIZ peripheral","L-Glutamine plasma",
                 "CFU-E, Reticulocytes","RBC sickle, RBC normal, HbF fraction",
                 "Hgb","free Hb, Haptoglobin, LDH, Bilirubin",
                 "NO, P-Selectin","VOC rate",
                 "NADH","Iron, TRV, eGFR"),
      `Count` = c(2,2,2,1,2,3,1,4,2,1,1,3),
      check.names=FALSE
    )
  })

  ## ── Treatment MoA table ──────────────────────────────────────
  output$treatment_moa <- renderTable({
    data.frame(
      Drug           = c("Hydroxyurea","Voxelotor","Crizanlizumab","L-Glutamine"),
      `Mechanism`    = c("RNR inhibitor → stress erythropoiesis → HbF↑",
                         "Anti-sickling: binds Hb α-chain, ↑O₂ affinity, prevents deoxy-HbS polymerization",
                         "Anti-P-selectin mAb: blocks RBC/WBC adhesion to endothelium",
                         "Antioxidant: NAD⁺ precursor, reduces RBC oxidative damage"),
      `Approval`     = c("FDA 1998 (adults), 2017 (children ≥2 yr)","FDA 2019","FDA 2019","FDA 2017"),
      `Key Trial`    = c("MSH (Charache 1995 NEJM)","HOPE (Vichinsky 2019 NEJM)",
                         "SUSTAIN (Ataga 2017 NEJM)","Phase III (Niihara 2018 NEJM)"),
      `Key Endpoint` = c("↓VOC 44%, ↓mortality","↑Hgb 1.0 g/dL, ↓sickling",
                         "↓VOC 45%","↓crises 25%"),
      check.names=FALSE
    )
  })

  ## ── ValueBoxes (Profile tab) ─────────────────────────────────
  output$box_genotype <- renderValueBox({
    valueBox("HbSS", "Genotype (Most Severe)", icon=icon("dna"), color="red")
  })
  output$box_hgb_base <- renderValueBox({
    valueBox("8.5 g/dL", "Baseline Hemoglobin", icon=icon("tint"), color="orange")
  })
  output$box_hbf_base <- renderValueBox({
    valueBox("7%", "Baseline HbF", icon=icon("shield-alt"), color="yellow")
  })
  output$box_voc_base <- renderValueBox({
    valueBox("3.5/yr", "Annual VOC Rate (untreated)", icon=icon("hospital"), color="red")
  })

  ## ── PK plot ──────────────────────────────────────────────────
  output$pk_plot <- renderPlotly({
    df <- sim_data()
    pk_long <- df %>%
      select(time_days, HU_plasma, VOX_RBC, CRIZ_C) %>%
      pivot_longer(-time_days, names_to="Drug", values_to="Conc") %>%
      mutate(Drug = recode(Drug,
        HU_plasma="Hydroxyurea (μM plasma)",
        VOX_RBC="Voxelotor (μM, RBC-bound)",
        CRIZ_C="Crizanlizumab (nM, central)"))

    p <- ggplot(pk_long %>% filter(time_days <= 28), aes(time_days, Conc, colour=Drug)) +
      geom_line(linewidth=1) +
      scale_colour_manual(values=c(scd_blue, scd_green, scd_red)) +
      facet_wrap(~Drug, scales="free_y") +
      labs(title="Drug PK Profiles (First 28 Days)",
           x="Time (days)", y="Concentration") +
      theme_scd() + theme(legend.position="none")
    ggplotly(p)
  })

  output$pk_params <- renderTable({
    data.frame(
      Parameter     = c("HU t½","HU Cmax (20 mg/kg)","VOX t½","VOX Cmax","CRIZ t½","CRIZ Cmax"),
      Value         = c("2–4 h","~50 μM","16–25 h","~2 μM","~10 days","~5 μg/mL"),
      Route         = c("Oral","Oral","Oral","Oral","IV","IV q4w"),
      check.names=FALSE
    )
  })

  ## ── Hematology plots ─────────────────────────────────────────
  make_plotly <- function(df, yvar, ytitle, ref_y=NULL, ref_label=NULL, color=scd_blue) {
    p <- ggplot(df, aes_string("time_days", yvar)) +
      geom_line(colour=color, linewidth=1.1) +
      labs(x="Time (days)", y=ytitle) +
      theme_scd()
    if (!is.null(ref_y)) {
      p <- p + geom_hline(yintercept=ref_y, linetype="dashed", colour="grey50") +
        annotate("text", x=5, y=ref_y+abs(ref_y)*0.02,
                 label=ref_label, size=3, colour="grey50")
    }
    ggplotly(p)
  }

  output$hgb_plot  <- renderPlotly({ make_plotly(sim_data(), "Hgb", "Hemoglobin (g/dL)", 9, "Target ≥9", scd_red) })
  output$hbf_plot  <- renderPlotly({ make_plotly(sim_data(), "HbF_pct", "HbF (%)", 20, "Target ≥20%", scd_green) })
  output$ret_plot  <- renderPlotly({ make_plotly(sim_data(), "ret_pct", "Reticulocytes (%)", color="#E67E22") })
  output$ldh_plot  <- renderPlotly({ make_plotly(sim_data(), "LDH", "LDH (U/L)", 250, "ULN 250 U/L", scd_red) })

  output$heme_table <- renderDT({
    df <- sim_data()
    snap <- df %>% filter(time_days %in% c(30, 90, 180, 365, 730)) %>%
      mutate(across(where(is.numeric), ~round(., 2))) %>%
      select(Day=time_days, Hgb, HbF_pct, ret_pct, LDH, Bilirubin)
    datatable(snap, options=list(dom="t", pageLength=10), rownames=FALSE) %>%
      formatStyle("Hgb", backgroundColor=styleInterval(c(8, 10), c("#FADBD8","#FEF9E7","#D5F5E3")))
  })

  ## ── VOC / Biomarker plots ──────────────────────────────────
  output$voc_plot  <- renderPlotly({ make_plotly(sim_data(), "VOC_annual", "Annual VOC Rate", color=scd_red) })
  output$psel_plot <- renderPlotly({ make_plotly(sim_data(), "P_selectin", "P-Selectin (relative)", color="#E67E22") })
  output$no_plot   <- renderPlotly({ make_plotly(sim_data(), "NO", "NO Bioavailability (a.u.)", color=scd_green) })
  output$bili_plot <- renderPlotly({ make_plotly(sim_data(), "Bilirubin", "Bilirubin (μmol/L)", 17, "ULN 17 μmol/L", "#F39C12") })

  output$voc_table <- renderDT({
    df <- sim_data()
    snap <- df %>% filter(time_days %in% c(30, 90, 180, 365)) %>%
      mutate(across(where(is.numeric), ~round(., 3))) %>%
      select(Day=time_days, VOC_annual, P_selectin, NO, Bilirubin, NADH, Haptoglobin)
    datatable(snap, options=list(dom="t"), rownames=FALSE)
  })

  ## ── Scenario comparison ───────────────────────────────────
  scen_map <- list(
    no_tx  = list(DOSE_HU=0,  DOSE_VOX=0, DOSE_CRIZ=0, DOSE_LG=0, criz=FALSE),
    hu20   = list(DOSE_HU=20, DOSE_VOX=0, DOSE_CRIZ=0, DOSE_LG=0, criz=FALSE),
    vox    = list(DOSE_HU=0,  DOSE_VOX=1, DOSE_CRIZ=0, DOSE_LG=0, criz=FALSE),
    criz   = list(DOSE_HU=0,  DOSE_VOX=0, DOSE_CRIZ=1, DOSE_LG=0, criz=TRUE),
    lg     = list(DOSE_HU=0,  DOSE_VOX=0, DOSE_CRIZ=0, DOSE_LG=1, criz=FALSE),
    hu_vox = list(DOSE_HU=20, DOSE_VOX=1, DOSE_CRIZ=0, DOSE_LG=0, criz=FALSE),
    triple = list(DOSE_HU=20, DOSE_VOX=1, DOSE_CRIZ=1, DOSE_LG=0, criz=TRUE)
  )
  scen_labels <- c(
    no_tx="No Treatment", hu20="Hydroxyurea 20 mg/kg",
    vox="Voxelotor", criz="Crizanlizumab",
    lg="L-Glutamine", hu_vox="HU + Voxelotor", triple="HU + VOX + CRIZ"
  )
  scen_colors <- c(
    no_tx="#7F8C8D", hu20="#2980B9", vox="#27AE60",
    criz="#E67E22", lg="#9B59B6", hu_vox="#C0392B", triple="#1ABC9C"
  )

  compare_data <- eventReactive(input$run_sim, ignoreNULL=FALSE, {
    withProgress(message="Running scenario comparisons...", value=0, {
      selected <- if(length(input$scen_check)==0) c("no_tx","hu20","vox","hu_vox") else input$scen_check
      lapply(selected, function(s) {
        setProgress(value=which(selected==s)/length(selected), message=paste("Simulating:", scen_labels[s]))
        p <- scen_map[[s]]
        df <- run_sim(mod_base, p$DOSE_HU, p$DOSE_VOX==1, p$criz, p$DOSE_LG==1,
                      input$bw, input$sim_days, if(p$criz) 5 else 0)
        df$scenario <- scen_labels[s]
        df$scen_id  <- s
        df
      }) %>% bind_rows()
    })
  })

  output$compare_plot <- renderPlotly({
    df <- compare_data()
    yvar <- input$compare_endpoint
    ylab <- switch(yvar, Hgb="Hemoglobin (g/dL)", HbF_pct="HbF (%)",
                   VOC_annual="VOC (crises/yr)", LDH="LDH (U/L)",
                   NO="NO Index", P_selectin="P-Selectin", TRV="TRV (m/s)")
    p <- ggplot(df, aes_string("time_days", yvar, colour="scenario")) +
      geom_line(linewidth=0.9) +
      scale_colour_manual(values=setNames(scen_colors[df$scen_id %>% unique()],
                                          scen_labels[df$scen_id %>% unique()])) +
      labs(x="Time (days)", y=ylab, colour="") + theme_scd()
    ggplotly(p)
  })

  output$compare_table <- renderDT({
    df <- compare_data()
    tbl <- df %>% filter(time_days >= input$sim_days - 1) %>%
      group_by(Scenario=scenario) %>%
      summarise(
        `Hgb (g/dL)`       = round(mean(Hgb), 2),
        `HbF (%)`          = round(mean(HbF_pct), 1),
        `VOC (crises/yr)`  = round(mean(VOC_annual), 2),
        `LDH (U/L)`        = round(mean(LDH), 0),
        `NO Index`         = round(mean(NO), 3),
        `TRV (m/s)`        = round(mean(TRV), 3),
        `eGFR`             = round(mean(eGFR), 1),
        .groups="drop")
    datatable(tbl, options=list(dom="t", pageLength=10), rownames=FALSE) %>%
      formatStyle("Hgb (g/dL)",
                  backgroundColor=styleInterval(c(8.5, 9.5), c("#FADBD8","#FEF9E7","#D5F5E3")))
  })

  ## ── End-Organ KPIs ───────────────────────────────────────
  kpi_data <- reactive({
    df <- sim_data()
    df %>% filter(time_days == max(time_days))
  })

  output$kpi_hgb <- renderValueBox({
    v <- round(kpi_data()$Hgb[1], 1)
    col <- if(v >= 9) "green" else if(v >= 7) "yellow" else "red"
    valueBox(paste0(v, " g/dL"), "Hemoglobin", icon=icon("tint"), color=col)
  })
  output$kpi_hbf <- renderValueBox({
    v <- round(kpi_data()$HbF_pct[1], 1)
    col <- if(v >= 20) "green" else if(v >= 10) "yellow" else "red"
    valueBox(paste0(v, "%"), "HbF", icon=icon("shield-alt"), color=col)
  })
  output$kpi_voc <- renderValueBox({
    v <- round(kpi_data()$VOC_annual[1], 1)
    col <- if(v <= 1) "green" else if(v <= 3) "yellow" else "red"
    valueBox(v, "VOC Rate (crises/yr)", icon=icon("hospital"), color=col)
  })
  output$kpi_ldh <- renderValueBox({
    v <- round(kpi_data()$LDH[1], 0)
    col <- if(v <= 300) "green" else if(v <= 500) "yellow" else "red"
    valueBox(paste0(v, " U/L"), "LDH", icon=icon("flask"), color=col)
  })

  output$trv_plot  <- renderPlotly({ make_plotly(sim_data(), "TRV", "TRV (m/s)", 2.5, "PH threshold", scd_red) })
  output$gfr_plot  <- renderPlotly({ make_plotly(sim_data(), "eGFR", "eGFR (mL/min/1.73m²)", color="#E67E22") })
  output$iron_plot <- renderPlotly({ make_plotly(sim_data(), "Iron", "Ferritin (μg/L)", 500, "Target <500", "#F39C12") })
  output$hp_plot   <- renderPlotly({ make_plotly(sim_data(), "Haptoglobin", "Haptoglobin (g/L)", 0.5, "Lower limit", scd_blue) })

  output$organ_table <- renderDT({
    df <- sim_data()
    snap <- df %>% filter(time_days %in% c(30, 90, 180, 365, 730)) %>%
      select(Day=time_days, TRV, eGFR, Iron, Haptoglobin,
             NO, P_selectin, NADH, Bilirubin) %>%
      mutate(across(where(is.numeric), ~round(., 3)))
    datatable(snap, options=list(dom="t"), rownames=FALSE) %>%
      formatStyle("TRV", backgroundColor=styleInterval(c(2.5, 3.0),
                                                         c("#D5F5E3","#FEF9E7","#FADBD8")))
  })
}

shinyApp(ui, server)
