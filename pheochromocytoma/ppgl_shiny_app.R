## ============================================================
##  Pheochromocytoma/Paraganglioma (PPGL) QSP — Shiny Dashboard
##  6 Tabs: Patient Profile · PK · Catecholamines/Biomarkers ·
##          Cardiovascular · Tumor (Malignant) · Scenario Compare
## ============================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(shinydashboard)

# ── Embedded mrgsolve model code ─────────────────────────────
ppgl_code <- '
$PARAM
ka_PHE=0.693 F_PHE=0.27 CL_PHE=5.8 V1_PHE=45 V2_PHE=180 Q_PHE=10
ka_DOX=0.21  F_DOX=0.65 V_DOX=80   t50_DOX=22
ka_MET=0.70  F_MET=0.85 CL_MET=4.5 V_MET=40 IC50_MET=85
ka_BB=0.80   F_BB=0.36  CL_BB=50   V_BB=260  IC50_BB=0.022
ka_SUNIT=0.28 F_SUNIT=0.60 CL_SUNIT=34 V1_SUNIT=2230 V2_SUNIT=1900 Q_SUNIT=6.5
ksynth_NE=0.35 ksynth_EPI=0.12 kstor=0.80 krel=0.045 krelease=1.0
kdeg_NE=0.62 kdeg_EPI=0.55 NE0=1.8 EPI0=0.28 NMN_factor=0.82
TH_base=1.0 kTH_deg=0.05 kTH_synth=0.05
TUM0=30.0 kgrowth_h=0.000125 VEGF0=120 kVEGF_synth=0.10 kVEGF_deg=0.20
IC50_SUNIT_VEGFR=0.018
SBP0=155 DBP0=98 HR0=92
E_NE_SBP=28.0 E_NE_HR=12.0 E_EPI_SBP=18.0 E_EPI_HR=22.0
kSBP=0.15 kHR=0.25 alpha_block_max=0.92
GLU0=5.5 kGLU=0.15 E_NE_GLU=0.8
FFA0=0.55 kFFA=0.12 E_EPI_FFA=0.25
CgA0=180 kCgA_synth=0.12 kCgA_deg=0.05
DOSE_PHE=0 DOSE_DOX=0 DOSE_MET=0 DOSE_BB=0 DOSE_SUNIT=0 SURGERY=0

$CMT PHE_gut PHE_C PHE_P DOX_C MET_C BB_C SUNIT_C SUNIT_P
     TH_act NE_store NE_plasma EPI_plasma TUMvol VEGF_tum
     SBP DBP HR GLU FFA CgA_plasma

$MAIN
PHE_gut_0=0; PHE_C_0=0; PHE_P_0=0; DOX_C_0=0; MET_C_0=0; BB_C_0=0;
SUNIT_C_0=0; SUNIT_P_0=0;
TH_act_0=TH_base;
NE_store_0=NE0/krel*kstor;
NE_plasma_0=NE0; EPI_plasma_0=EPI0;
TUMvol_0=TUM0; VEGF_tum_0=VEGF0;
SBP_0=SBP0; DBP_0=DBP0; HR_0=HR0;
GLU_0=GLU0; FFA_0=FFA0;
CgA_plasma_0=CgA0;

$ODE
dxdt_PHE_gut= -ka_PHE*PHE_gut;
dxdt_PHE_C = ka_PHE*F_PHE*PHE_gut/V1_PHE-(CL_PHE/V1_PHE)*PHE_C-(Q_PHE/V1_PHE)*PHE_C+(Q_PHE/V2_PHE)*PHE_P;
dxdt_PHE_P = (Q_PHE/V1_PHE)*PHE_C-(Q_PHE/V2_PHE)*PHE_P;
double CL_DOX_c=0.693/t50_DOX*V_DOX;
dxdt_DOX_C = ka_DOX*F_DOX*DOSE_DOX/V_DOX-(CL_DOX_c/V_DOX)*DOX_C;
dxdt_MET_C = ka_MET*F_MET*DOSE_MET/V_MET-(CL_MET/V_MET)*MET_C;
dxdt_BB_C  = ka_BB*F_BB*DOSE_BB/V_BB-(CL_BB/V_BB)*BB_C;
dxdt_SUNIT_C=ka_SUNIT*F_SUNIT*DOSE_SUNIT/V1_SUNIT-(CL_SUNIT/V1_SUNIT)*SUNIT_C-(Q_SUNIT/V1_SUNIT)*SUNIT_C+(Q_SUNIT/V2_SUNIT)*SUNIT_P;
dxdt_SUNIT_P=(Q_SUNIT/V1_SUNIT)*SUNIT_C-(Q_SUNIT/V2_SUNIT)*SUNIT_P;
double MET_inhib=MET_C/(MET_C+IC50_MET);
double TH_target=TH_base*(1-0.80*MET_inhib);
dxdt_TH_act=kTH_synth*TH_target-kTH_deg*TH_act;
double tumFactor=TUMvol/TUM0;
double synth_NE=ksynth_NE*TH_act*tumFactor;
double synth_EPI=ksynth_EPI*TH_act*tumFactor;
double krel_eff=krel*(1+krelease*(tumFactor-1));
dxdt_NE_store=kstor*synth_NE-krel_eff*NE_store;
dxdt_NE_plasma=krel_eff*NE_store-kdeg_NE*NE_plasma;
dxdt_EPI_plasma=synth_EPI-kdeg_EPI*EPI_plasma;
double SUNIT_inh=SUNIT_C/(SUNIT_C+IC50_SUNIT_VEGFR);
double kgrowth_eff=kgrowth_h*(1-0.65*SUNIT_inh);
double surgery_factor=(SURGERY>0)?0:1;
dxdt_TUMvol=kgrowth_eff*TUMvol*surgery_factor;
dxdt_VEGF_tum=kVEGF_synth*TUMvol-kVEGF_deg*VEGF_tum;
double PHE_eff=PHE_C/(PHE_C+0.012);
double DOX_eff=DOX_C/(DOX_C+0.002);
double alpha_block=1-(1-PHE_eff)*(1-DOX_eff);
double NE_excess=NE_plasma-NE0;
double EPI_excess=EPI_plasma-EPI0;
double SBP_target=SBP0+E_NE_SBP*(NE_excess>0?NE_excess:0)+E_EPI_SBP*(EPI_excess>0?EPI_excess:0);
SBP_target=SBP_target*(1-alpha_block_max*alpha_block);
double DBP_target=DBP0+0.65*(SBP_target-SBP0);
dxdt_SBP=kSBP*(SBP_target-SBP);
dxdt_DBP=kSBP*(DBP_target-DBP);
double beta_block=BB_C/(BB_C+IC50_BB);
double HR_target=HR0+E_NE_HR*(NE_excess>0?NE_excess:0)+E_EPI_HR*(EPI_excess>0?EPI_excess:0);
HR_target=HR_target*(1-0.40*beta_block);
dxdt_HR=kHR*(HR_target-HR);
double GLU_target=GLU0+E_NE_GLU*(NE_excess>0?NE_excess:0)+0.5*(EPI_excess>0?EPI_excess:0);
dxdt_GLU=kGLU*(GLU_target-GLU);
double FFA_target=FFA0+E_EPI_FFA*(EPI_excess>0?EPI_excess:0);
dxdt_FFA=kFFA*(FFA_target-FFA);
dxdt_CgA_plasma=kCgA_synth*TUMvol-kCgA_deg*CgA_plasma;

$TABLE
capture SBP_mmHg=SBP; capture DBP_mmHg=DBP; capture MAP_mmHg=DBP+(SBP-DBP)/3.0;
capture HR_bpm=HR; capture NE_pl=NE_plasma; capture EPI_pl=EPI_plasma;
capture NMN_pl=NE_plasma*NMN_factor; capture MN_pl=EPI_plasma*1.15;
capture TH_activity=TH_act; capture Tumor_mL=TUMvol; capture VEGF_pg=VEGF_tum;
capture Glucose_mM=GLU; capture FFA_mM=FFA; capture CgA_ng=CgA_plasma;
capture PHE_Conc=PHE_C; capture DOX_Conc=DOX_C; capture MET_Conc=MET_C;
capture BB_Conc=BB_C; capture SUNIT_Conc=SUNIT_C;
double alpha_blockFrac=1-(1-PHE_C/(PHE_C+0.012))*(1-DOX_C/(DOX_C+0.002));
capture AlphaBlock=alpha_blockFrac;
'

mod_ppgl <- mcode("PPGL_shiny", ppgl_code)

# ── Helper: create dosing events ─────────────────────────────
make_dose <- function(amt, cmt, ii, addl) {
  ev(amt=amt, cmt=cmt, ii=ii, addl=addl, time=0)
}

# ── UI ────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "PPGL QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile", tabName="profile",  icon=icon("user")),
      menuItem("Drug PK",         tabName="pk",       icon=icon("flask")),
      menuItem("Catecholamines",  tabName="catechol", icon=icon("vial")),
      menuItem("Cardiovascular",  tabName="cardio",   icon=icon("heartbeat")),
      menuItem("Tumor (Malignant)", tabName="tumor",  icon=icon("dna")),
      menuItem("Scenario Compare",  tabName="compare",icon=icon("chart-bar"))
    ),
    br(),
    div(style="padding:10px;",
      h5("Global Controls", style="color:white;"),
      sliderInput("sim_days", "Simulation duration (days)", 7, 90, 30, step=7),
      actionButton("simulate", "Run Simulation", class="btn-warning btn-block",
                   icon=icon("play"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f4f4; }
      .box-title { font-weight: bold; }
    "))),

    tabItems(

      # ── TAB 1: Patient Profile ──────────────────────────────
      tabItem(tabName="profile",
        fluidRow(
          box(title="Patient & Tumor Characteristics", width=6, status="purple",
            selectInput("ppgl_type", "PPGL Type",
                        choices=c("Sporadic Pheochromocytoma (adrenal)",
                                  "SDHB-mutant PPGL (malignant risk)",
                                  "VHL-associated Pheochromocytoma",
                                  "MEN2 / RET-mutant Pheochromocytoma",
                                  "NF1 Pheochromocytoma",
                                  "Metastatic / Malignant PPGL")),
            numericInput("tum_vol", "Baseline Tumor Volume (mL)", 30, 5, 200),
            numericInput("ne0", "Baseline Plasma NE (nmol/L)", 1.8, 0.5, 20),
            numericInput("epi0", "Baseline Plasma EPI (nmol/L)", 0.28, 0.05, 10),
            numericInput("sbp0", "Baseline SBP (mmHg)", 155, 100, 250),
            numericInput("hr0", "Baseline HR (bpm)", 92, 50, 130)
          ),
          box(title="Biochemical Profile", width=6, status="warning",
            numericInput("cga0", "Plasma CgA (ng/mL)", 180, 20, 2000),
            sliderInput("nmn_mult", "NMN multiplier (×ULN)", 1, 50, 8),
            selectInput("mutation", "Germline Mutation Status",
                        choices=c("None / Sporadic","SDHB","SDHD","SDHC",
                                  "VHL","RET (MEN2)","NF1","MAX","TMEM127")),
            checkboxInput("bilateral", "Bilateral PPGL", FALSE),
            checkboxInput("extraadrenal", "Extra-adrenal (Paraganglioma)", FALSE),
            actionButton("set_profile", "Apply Profile", class="btn-primary btn-block")
          )
        ),
        fluidRow(
          box(title="Disease Summary", width=12, status="info",
            tableOutput("profile_table")
          )
        )
      ),

      # ── TAB 2: Drug PK ─────────────────────────────────────
      tabItem(tabName="pk",
        fluidRow(
          box(title="Alpha-Blocker Selection", width=4, status="blue",
            selectInput("alpha_drug", "Preoperative Alpha-Blocker",
                        choices=c("None","Phenoxybenzamine (irreversible)",
                                  "Doxazosin (selective α₁)")),
            conditionalPanel("input.alpha_drug != 'None'",
              sliderInput("phe_dose", "Phenoxybenzamine total dose (mg/d)", 20, 120, 60),
              sliderInput("dox_dose", "Doxazosin dose (mg/d)", 2, 32, 16),
              sliderInput("preop_days", "Preop preparation (days)", 7, 28, 14)
            )
          ),
          box(title="Adjunct Drugs", width=4, status="orange",
            checkboxInput("use_met", "Metyrosine (TH inhibitor)", FALSE),
            conditionalPanel("input.use_met",
              sliderInput("met_dose", "Metyrosine dose (mg/d)", 500, 4000, 2000)),
            checkboxInput("use_bb", "Beta-Blocker (after alpha)", FALSE),
            conditionalPanel("input.use_bb",
              sliderInput("bb_dose", "Propranolol dose (mg/d)", 40, 240, 120))
          ),
          box(title="Malignant PPGL Systemic Therapy", width=4, status="red",
            checkboxInput("use_sunit", "Sunitinib (malignant PPGL)", FALSE),
            conditionalPanel("input.use_sunit",
              sliderInput("sunit_dose", "Sunitinib dose (mg/d)", 25, 50, 37.5, step=12.5)),
            checkboxInput("surgery", "Surgery at end of preop", FALSE)
          )
        ),
        fluidRow(
          box(title="PK Profiles", width=12, status="primary",
            plotOutput("pk_plot", height="400px"))
        )
      ),

      # ── TAB 3: Catecholamines & Biomarkers ─────────────────
      tabItem(tabName="catechol",
        fluidRow(
          box(title="Plasma Catecholamines & Fractionated Metanephrines",
              width=12, status="success",
            fluidRow(
              column(6, plotOutput("ne_epi_plot", height="300px")),
              column(6, plotOutput("nmn_mn_plot", height="300px"))
            )
          )
        ),
        fluidRow(
          box(title="TH Activity & Synthesis", width=6, status="warning",
            plotOutput("th_plot", height="280px")),
          box(title="Chromogranin-A (Tumor Marker)", width=6, status="info",
            plotOutput("cga_plot", height="280px"))
        ),
        fluidRow(
          box(title="Biochemical Reference Ranges", width=12, status="primary",
            tableOutput("biomarker_table"))
        )
      ),

      # ── TAB 4: Cardiovascular ──────────────────────────────
      tabItem(tabName="cardio",
        fluidRow(
          box(title="Blood Pressure Control", width=6, status="danger",
            plotOutput("bp_plot", height="320px")),
          box(title="Heart Rate Response", width=6, status="warning",
            plotOutput("hr_plot", height="320px"))
        ),
        fluidRow(
          box(title="Alpha-Blockade Fraction & Cardiovascular Risk",
              width=6, status="purple",
            plotOutput("ab_frac_plot", height="280px")),
          box(title="Metabolic Effects (Glucose & FFA)", width=6, status="orange",
            plotOutput("metabolic_plot", height="280px"))
        ),
        fluidRow(
          box(title="Hypertensive Crisis Risk Assessment", width=12, status="danger",
            infoBoxOutput("crisis_risk", width=4),
            infoBoxOutput("sbp_control", width=4),
            infoBoxOutput("hr_control", width=4))
        )
      ),

      # ── TAB 5: Tumor Dynamics (Malignant) ──────────────────
      tabItem(tabName="tumor",
        fluidRow(
          box(title="Tumor Volume Over Time", width=8, status="danger",
            plotOutput("tumor_vol_plot", height="360px")),
          box(title="Tumor Parameters", width=4, status="warning",
            sliderInput("tum_growth", "Net growth rate (% per month)", 0.5, 15, 4, step=0.5),
            sliderInput("sunit_sens", "Sunitinib sensitivity (VEGFR IC50 shift ×)", 0.1, 5, 1, step=0.1),
            hr(),
            h5("RECIST Assessment (Day 90)"),
            uiOutput("recist_ui")
          )
        ),
        fluidRow(
          box(title="VEGF & Angiogenesis Inhibition", width=6, status="info",
            plotOutput("vegf_plot", height="280px")),
          box(title="Malignant PPGL — Treatment Options Overview", width=6, status="primary",
            tableOutput("treatment_options_table"))
        )
      ),

      # ── TAB 6: Scenario Comparison ─────────────────────────
      tabItem(tabName="compare",
        fluidRow(
          box(title="Scenario Selection", width=4, status="primary",
            checkboxGroupInput("scenarios_compare",
              "Compare Scenarios:",
              choices=c("No treatment (S0)"=0,
                        "PHE preop → Surgery (S1)"=1,
                        "Doxazosin preop → Surgery (S2)"=2,
                        "PHE + MET + BB triple (S3)"=3,
                        "Sunitinib malignant (S4)"=4,
                        "Metyrosine mono (S5)"=5),
              selected=c(0,1,2,3)),
            selectInput("compare_var", "Outcome Variable",
                        choices=c("SBP (mmHg)"="SBP_mmHg",
                                  "DBP (mmHg)"="DBP_mmHg",
                                  "HR (bpm)"="HR_bpm",
                                  "NMN (nmol/L)"="NMN_pl",
                                  "MN (nmol/L)"="MN_pl",
                                  "TH Activity"="TH_activity",
                                  "CgA (ng/mL)"="CgA_ng",
                                  "Tumor Volume (mL)"="Tumor_mL",
                                  "Glucose (mmol/L)"="Glucose_mM")),
            actionButton("run_compare", "Compare Now", class="btn-success btn-block",
                         icon=icon("chart-line"))
          ),
          box(title="Comparison Plot", width=8, status="success",
            plotOutput("compare_plot", height="400px"))
        ),
        fluidRow(
          box(title="Summary Table — All Scenarios at Day 14 & Day 30", width=12, status="info",
            tableOutput("compare_table"))
        )
      )
    )
  )
)

# ── SERVER ────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Reactive simulation ─────────────────────────────────────
  sim_data <- eventReactive(list(input$simulate, input$run_compare), {
    end_h <- input$sim_days * 24
    dose_ii_phe <- 12
    dose_ii_dox <- 24
    dose_ii_met <- 6
    dose_ii_bb  <- 8
    dose_ii_su  <- 24
    dur_h <- end_h

    # Build events per scenario
    run_one <- function(scen_id) {
      evs <- NULL
      if (scen_id == 1) {
        evs <- ev(amt=input$phe_dose/dose_ii_phe/303.8, cmt="PHE_gut",
                  ii=dose_ii_phe, addl=round(dur_h/dose_ii_phe)-1)
      } else if (scen_id == 2) {
        evs <- ev(amt=input$dox_dose/dose_ii_dox/451.5, cmt="DOX_C",
                  ii=dose_ii_dox, addl=round(dur_h/dose_ii_dox)-1)
      } else if (scen_id == 3) {
        evs <- c(ev(amt=input$phe_dose/dose_ii_phe/303.8, cmt="PHE_gut",
                    ii=dose_ii_phe, addl=round(dur_h/dose_ii_phe)-1),
                 ev(amt=input$met_dose/dose_ii_met/195.2, cmt="MET_C",
                    ii=dose_ii_met, addl=round(dur_h/dose_ii_met)-1),
                 ev(amt=input$bb_dose/dose_ii_bb/259.3, cmt="BB_C",
                    ii=dose_ii_bb, addl=round(dur_h/dose_ii_bb)-1))
      } else if (scen_id == 4) {
        evs <- ev(amt=input$sunit_dose/dose_ii_su/532.6, cmt="SUNIT_C",
                  ii=dose_ii_su, addl=round(dur_h/dose_ii_su)-1)
      } else if (scen_id == 5) {
        evs <- ev(amt=input$met_dose/dose_ii_met/195.2, cmt="MET_C",
                  ii=dose_ii_met, addl=round(dur_h/dose_ii_met)-1)
      }

      p_override <- list(NE0=input$ne0, EPI0=input$epi0,
                         SBP0=input$sbp0, HR0=input$hr0,
                         TUM0=input$tum_vol, CgA0=input$cga0)

      tryCatch({
        mod2 <- param(mod_ppgl, p_override)
        if (is.null(evs)) {
          out <- mrgsim(mod2, end=end_h, delta=1) %>% as.data.frame()
        } else {
          out <- mrgsim(mod2, events=evs, end=end_h, delta=1) %>% as.data.frame()
        }
        out$scenario <- scen_id
        out$scenario_label <- c("0"="No treatment","1"="PHE preop","2"="DOX preop",
                                 "3"="PHE+MET+BB","4"="Sunitinib","5"="Metyrosine")[as.character(scen_id)]
        out
      }, error=function(e) NULL)
    }

    scens <- if (!is.null(input$scenarios_compare) && length(input$scenarios_compare)>0) {
      as.integer(input$scenarios_compare)
    } else { 0:3 }

    bind_rows(lapply(scens, run_one))
  }, ignoreNULL=FALSE)

  # ── Profile Table ───────────────────────────────────────────
  output$profile_table <- renderTable({
    data.frame(
      Parameter = c("PPGL Type","Tumor Volume","Baseline NE","Baseline EPI",
                    "Baseline SBP","Baseline HR","Chromogranin-A","Mutation"),
      Value = c(input$ppgl_type, paste0(input$tum_vol," mL"),
                paste0(input$ne0," nmol/L"), paste0(input$epi0," nmol/L"),
                paste0(input$sbp0," mmHg"), paste0(input$hr0," bpm"),
                paste0(input$cga0," ng/mL"), input$mutation),
      stringsAsFactors=FALSE)
  })

  # ── PK Plot ─────────────────────────────────────────────────
  output$pk_plot <- renderPlot({
    df <- sim_data()
    req(nrow(df)>0)
    pk_long <- df %>%
      mutate(day=time/24) %>%
      pivot_longer(c(PHE_Conc, DOX_Conc, MET_Conc, BB_Conc, SUNIT_Conc),
                   names_to="Drug", values_to="Conc") %>%
      filter(Conc > 0)
    if (nrow(pk_long)==0) { plot.new(); text(0.5,0.5,"No active drug doses"); return() }
    ggplot(pk_long, aes(x=day, y=Conc, color=Drug, linetype=scenario_label)) +
      geom_line(size=1) + scale_y_log10() +
      scale_color_brewer(palette="Set1") +
      labs(title="Drug PK — Plasma Concentrations", x="Day", y="Conc (µmol/L, log)", color="Drug", linetype="Scenario") +
      theme_bw(base_size=13)
  })

  # ── NE/EPI Plot ─────────────────────────────────────────────
  output$ne_epi_plot <- renderPlot({
    df <- sim_data()
    req(nrow(df)>0)
    ggplot(df, aes(x=time/24, color=scenario_label)) +
      geom_line(aes(y=NE_pl, linetype="NE"), size=1) +
      geom_hline(yintercept=input$ne0, linetype="dotted", color="gray40") +
      scale_color_brewer(palette="Set2") +
      labs(title="Plasma NE", x="Day", y="NE (nmol/L)", color="Scenario", linetype="") +
      theme_bw(base_size=12)
  })

  output$nmn_mn_plot <- renderPlot({
    df <- sim_data()
    req(nrow(df)>0)
    ggplot(df, aes(x=time/24)) +
      geom_line(aes(y=NMN_pl, color=scenario_label), size=1) +
      geom_hline(yintercept=0.9, linetype="dashed", color="red", linewidth=0.8) +
      annotate("text",x=0.5,y=0.92,label="3× ULN",size=3,color="red",hjust=0) +
      scale_color_brewer(palette="Set2") +
      labs(title="Plasma Normetanephrine (NMN)", x="Day", y="NMN (nmol/L)", color="Scenario") +
      theme_bw(base_size=12)
  })

  output$th_plot <- renderPlot({
    df <- sim_data()
    req(nrow(df)>0)
    ggplot(df, aes(x=time/24, y=TH_activity, color=scenario_label)) +
      geom_line(size=1) +
      geom_hline(yintercept=0.5, linetype="dashed", color="orange") +
      scale_color_brewer(palette="Set2") +
      labs(title="TH Enzyme Activity", x="Day", y="Relative TH Activity", color="Scenario") +
      theme_bw(base_size=12) + ylim(0, 1.1)
  })

  output$cga_plot <- renderPlot({
    df <- sim_data()
    req(nrow(df)>0)
    ggplot(df, aes(x=time/24, y=CgA_ng, color=scenario_label)) +
      geom_line(size=1) +
      geom_hline(yintercept=101, linetype="dashed", color="gray50") +
      annotate("text",x=0.5,y=103,label="ULN 101 ng/mL",size=3,hjust=0,color="gray40") +
      scale_color_brewer(palette="Set2") +
      labs(title="Plasma Chromogranin-A", x="Day", y="CgA (ng/mL)", color="Scenario") +
      theme_bw(base_size=12)
  })

  output$biomarker_table <- renderTable({
    data.frame(
      Biomarker = c("Plasma NMN","Plasma MN","Plasma NE","Plasma EPI",
                    "24h Urine NMN","24h Urine MN","Plasma CgA","NSE"),
      Normal = c("<0.29 nmol/L","<0.28 nmol/L","<0.6 nmol/L","<0.28 nmol/L",
                 "<0.9 µmol/d","<0.35 µmol/d","<101 ng/mL","<12 ng/mL"),
      `PPGL Cutoff (3×ULN)` = c("≥0.87","≥0.85","≥1.8","≥0.85",
                                  "≥2.7 µmol/d","≥1.05 µmol/d","≥300","≥36"),
      Sensitivity = c("97%","99%","85%","90%","87%","82%","70%","55%"),
      stringsAsFactors=FALSE, check.names=FALSE
    )
  })

  # ── BP & HR Plots ───────────────────────────────────────────
  output$bp_plot <- renderPlot({
    df <- sim_data()
    req(nrow(df)>0)
    ggplot(df, aes(x=time/24, color=scenario_label)) +
      geom_line(aes(y=SBP_mmHg), size=1) +
      geom_line(aes(y=DBP_mmHg), size=0.7, linetype="dashed") +
      geom_hline(yintercept=140, linetype="dotted", color="red") +
      geom_hline(yintercept=90, linetype="dotted", color="orange") +
      scale_color_brewer(palette="Set1") +
      labs(title="Blood Pressure (solid=SBP, dash=DBP)", x="Day", y="BP (mmHg)", color="Scenario") +
      theme_bw(base_size=12)
  })

  output$hr_plot <- renderPlot({
    df <- sim_data()
    req(nrow(df)>0)
    ggplot(df, aes(x=time/24, y=HR_bpm, color=scenario_label)) +
      geom_line(size=1) +
      geom_hline(yintercept=c(60,100), linetype="dashed", color=c("blue","red")) +
      scale_color_brewer(palette="Set1") +
      labs(title="Heart Rate", x="Day", y="HR (bpm)", color="Scenario") +
      theme_bw(base_size=12)
  })

  output$ab_frac_plot <- renderPlot({
    df <- sim_data()
    req(nrow(df)>0)
    ggplot(df, aes(x=time/24, y=AlphaBlock*100, color=scenario_label)) +
      geom_line(size=1) +
      geom_hline(yintercept=80, linetype="dashed", color="green4") +
      annotate("text",x=0.5,y=82,label="Target >80%",size=3,color="green4",hjust=0) +
      scale_color_brewer(palette="Set1") +
      labs(title="Alpha-Receptor Blockade (%)", x="Day", y="Blockade (%)", color="Scenario") +
      theme_bw(base_size=12) + ylim(0,100)
  })

  output$metabolic_plot <- renderPlot({
    df <- sim_data()
    req(nrow(df)>0)
    ggplot(df, aes(x=time/24, color=scenario_label)) +
      geom_line(aes(y=Glucose_mM), size=1) +
      geom_hline(yintercept=11.1, linetype="dashed", color="red") +
      annotate("text",x=0.5,y=11.3,label="Hyperglycemia",size=3,color="red",hjust=0) +
      scale_color_brewer(palette="Set1") +
      labs(title="Plasma Glucose", x="Day", y="Glucose (mmol/L)", color="Scenario") +
      theme_bw(base_size=12)
  })

  output$crisis_risk <- renderInfoBox({
    df <- sim_data()
    if (is.null(df) || nrow(df)==0) return(infoBox("Crisis Risk","N/A",icon=icon("exclamation")))
    max_sbp <- max(df$SBP_mmHg, na.rm=TRUE)
    risk <- if (max_sbp >= 180) "HIGH" else if (max_sbp >= 160) "MODERATE" else "LOW"
    col  <- if (max_sbp >= 180) "red" else if (max_sbp >= 160) "orange" else "green"
    infoBox("Crisis Risk", paste0(risk, "\n(Max SBP=",round(max_sbp),"mmHg)"),
            icon=icon("bolt"), color=col)
  })

  output$sbp_control <- renderInfoBox({
    df <- sim_data()
    if (is.null(df) || nrow(df)==0) return(infoBox("SBP at Day 14","N/A",icon=icon("tachometer-alt")))
    d14 <- df %>% filter(abs(time - 336) < 2)
    sbp14 <- if (nrow(d14)>0) round(mean(d14$SBP_mmHg),0) else NA
    col <- if (!is.na(sbp14) && sbp14 < 140) "green" else "red"
    infoBox("SBP at Day 14", paste0(sbp14," mmHg (target <140)"),
            icon=icon("heart"), color=col)
  })

  output$hr_control <- renderInfoBox({
    df <- sim_data()
    if (is.null(df) || nrow(df)==0) return(infoBox("HR at Day 14","N/A",icon=icon("heartbeat")))
    d14 <- df %>% filter(abs(time - 336) < 2)
    hr14 <- if (nrow(d14)>0) round(mean(d14$HR_bpm),0) else NA
    col <- if (!is.na(hr14) && hr14 <= 80) "green" else "orange"
    infoBox("HR at Day 14", paste0(hr14," bpm (target 60-80)"),
            icon=icon("heartbeat"), color=col)
  })

  # ── Tumor Plots ─────────────────────────────────────────────
  output$tumor_vol_plot <- renderPlot({
    df <- sim_data()
    req(nrow(df)>0)
    ggplot(df %>% filter(scenario %in% c(0,4)), aes(x=time/24, y=Tumor_mL, color=scenario_label)) +
      geom_line(size=1.2) +
      geom_hline(yintercept=input$tum_vol*0.70, linetype="dashed", color="green4") +
      annotate("text", x=1, y=input$tum_vol*0.72, label="PR threshold (−30%)", size=3, hjust=0, color="green4") +
      scale_color_manual(values=c("No treatment"="#E53935","Sunitinib"="#8E24AA")) +
      labs(title="Tumor Volume Dynamics", x="Day", y="Tumor Volume (mL)", color="Scenario") +
      theme_bw(base_size=13)
  })

  output$recist_ui <- renderUI({
    df <- sim_data()
    if (is.null(df) || nrow(df)==0) return(p("Run simulation first."))
    d90 <- df %>% filter(abs(time - min(max(df$time), 2160)) < 2, scenario==4)
    if (nrow(d90)==0) d90 <- df %>% filter(scenario==4) %>% tail(1)
    if (nrow(d90)==0) return(p("No sunitinib data."))
    vol <- round(mean(d90$Tumor_mL),1)
    pct <- round((vol - input$tum_vol)/input$tum_vol*100, 1)
    resp <- if (pct <= -30) "Partial Response (PR)" else if (pct <= 20) "Stable Disease (SD)" else "Progressive Disease (PD)"
    col  <- if (pct <= -30) "green" else if (pct <= 20) "orange" else "red"
    tags$div(style=paste0("background:",col,";color:white;padding:8px;border-radius:4px;"),
             strong(resp), br(), paste0("Volume: ",vol," mL (Δ=",pct,"%)"))
  })

  output$vegf_plot <- renderPlot({
    df <- sim_data()
    req(nrow(df)>0)
    ggplot(df, aes(x=time/24, y=VEGF_pg, color=scenario_label)) +
      geom_line(size=1) +
      scale_color_brewer(palette="Set2") +
      labs(title="Plasma VEGF Dynamics", x="Day", y="VEGF-A (pg/mL)", color="Scenario") +
      theme_bw(base_size=12)
  })

  output$treatment_options_table <- renderTable({
    data.frame(
      Treatment = c("Sunitinib","Cyclophosphamide/Vincristine/Dacarbazine (CVD)",
                    "¹³¹I-MIBG","¹⁷⁷Lu-DOTATATE","Cabozantinib","Pembrolizumab"),
      `Objective RR` = c("25-40%","37%","13-22%","25%","19%","25%"),
      `Median PFS` = c("13 mo","18 mo","12 mo","NR","5 mo","6 mo"),
      Evidence = c("Phase II","Phase II (Averbuch)","Phase II","Phase II","Phase II","Phase II"),
      stringsAsFactors=FALSE, check.names=FALSE
    )
  })

  # ── Compare Plots ────────────────────────────────────────────
  output$compare_plot <- renderPlot({
    df <- sim_data()
    req(nrow(df)>0)
    var <- input$compare_var
    if (!(var %in% names(df))) return()
    ggplot(df, aes_string(x="time/24", y=var, color="scenario_label")) +
      geom_line(size=1.1) +
      scale_color_brewer(palette="Set1") +
      labs(title=paste("Scenario Comparison —", var), x="Day", y=var, color="Scenario") +
      theme_bw(base_size=13)
  })

  output$compare_table <- renderTable({
    df <- sim_data()
    req(nrow(df)>0)
    df %>%
      filter(time %in% c(0, 168, 336, 504, 720)) %>%
      group_by(scenario_label, time) %>%
      summarise(Day=unique(time)/24,
                SBP=round(mean(SBP_mmHg),0),
                DBP=round(mean(DBP_mmHg),0),
                HR=round(mean(HR_bpm),0),
                NMN=round(mean(NMN_pl),2),
                CgA=round(mean(CgA_ng),0),
                TumVol=round(mean(Tumor_mL),1),
                .groups="drop") %>%
      rename(Scenario=scenario_label)
  })
}

shinyApp(ui=ui, server=server)
