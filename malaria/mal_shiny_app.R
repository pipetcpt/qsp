## ============================================================
## Malaria QSP Shiny Dashboard
## 7 Tabs: Patient/Infection Profile · PK (Artemisinin/Partner) ·
##         Parasite Dynamics · Clinical Endpoints ·
##         Scenario Comparison · Biomarkers · Resistance (K13)
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(tidyr)
library(DT)

# ─────────────────────────────────────────────────────────────
# mrgsolve Model (embedded; mirrors mal_mrgsolve_model.R)
# ─────────────────────────────────────────────────────────────
mal_code <- '
$PARAM
  use_AS=0,use_LUM=0,use_PPQ=0,use_AQ=0,use_PQ=0,K13_RES=0,
  ka_AS=4.0,ke_AS=1.386,Vc_AS=30.0,ke_DHA=0.90,Vc_DHA=70.0,frac_AS_DHA=0.90,
  ka_LUM=0.30,CL_LUM=0.14,Vc_LUM=21.0,Q_LUM=0.05,Vp_LUM=500.0,
  ka_PPQ=0.10,CL_PPQ=0.90,Vc_PPQ=150.0,Q_PPQ=0.60,Vp_PPQ=5000.0,
  ka_AQ=0.50,ke_DEAQ=0.0032,Vc_DEAQ=1500.0,
  ka_PQ=0.60,ke_PQ=0.11,Vc_PQ=220.0,
  Emax_AS=0.42,EC50_AS=5.0,hill_AS=1.8,Kres_shift=6.0,
  Emax_LUM=0.10,EC50_LUM=280.0,hill_LUM=1.0,
  Emax_PPQ=0.10,EC50_PPQ=30.0,hill_PPQ=1.0,
  Emax_AQ=0.10,EC50_AQ=40.0,hill_AQ=1.0,
  Emax_PQ_liver=0.06,EC50_PQ=15.0,hill_PQ=1.0,Emax_PQ_gam=0.30,
  G6PD_deficient=0,k_PQ_hemolysis=0.0006,
  k_RT=0.0556,k_TS=0.0833,k_SR=0.0556,burst=16.0,
  inv_eff0=0.45,immune_block=0.85,gam_commit=0.015,k_gam_mat=0.0052,
  k_liver_egress=0.0060,relapse_on=0,k_relapse=0.00028,
  RBC_base=5.0e6,k_RBC_sen=3.47e-4,k_ery_comp=0.0020,
  HB_base=13.5,hb_per_rupture=2.5e-6,k_dysery=3.0e-6,
  Tbase=37.0,Tmax_rise=3.2,K_fever=2.0e4,k_temp_decay=0.25,
  k_immune_gain=3.0e-7,k_immune_decay=0.0010,K_immune=5.0e4

$INIT
  AS_GUT=0,AS_PLASMA=0,DHA_PLASMA=0,LUM_GUT=0,LUM_CENTRAL=0,LUM_PERIPH=0,
  PPQ_GUT=0,PPQ_CENTRAL=0,PPQ_PERIPH=0,AQ_GUT=0,DEAQ_PLASMA=0,PQ_GUT=0,PQ_PLASMA=0,
  RBC_U=5.0e6,PRBC_RING=10.0,PRBC_TROPH=6.0,PRBC_SCHIZONT=4.0,
  LIVER_PARASITE=1.0e4,GAMETOCYTE=0.1,HB=13.5,TEMP=37.0,IMMUNITY=0.05

$ODE
  dxdt_AS_GUT    = -ka_AS*AS_GUT;
  dxdt_AS_PLASMA =  ka_AS*AS_GUT/Vc_AS - ke_AS*AS_PLASMA;
  dxdt_DHA_PLASMA=  frac_AS_DHA*ke_AS*AS_PLASMA*Vc_AS/Vc_DHA - ke_DHA*DHA_PLASMA;

  double keLUM=CL_LUM/Vc_LUM;
  double kcpLUM=Q_LUM/Vc_LUM;
  double kpcLUM=Q_LUM/Vp_LUM;
  dxdt_LUM_GUT    = -ka_LUM*LUM_GUT;
  dxdt_LUM_CENTRAL=  ka_LUM*LUM_GUT/Vc_LUM - keLUM*LUM_CENTRAL - kcpLUM*(LUM_CENTRAL-LUM_PERIPH);
  dxdt_LUM_PERIPH =  kpcLUM*(LUM_CENTRAL-LUM_PERIPH);

  double kePPQ=CL_PPQ/Vc_PPQ;
  double kcpPPQ=Q_PPQ/Vc_PPQ;
  double kpcPPQ=Q_PPQ/Vp_PPQ;
  dxdt_PPQ_GUT    = -ka_PPQ*PPQ_GUT;
  dxdt_PPQ_CENTRAL=  ka_PPQ*PPQ_GUT/Vc_PPQ - kePPQ*PPQ_CENTRAL - kcpPPQ*(PPQ_CENTRAL-PPQ_PERIPH);
  dxdt_PPQ_PERIPH =  kpcPPQ*(PPQ_CENTRAL-PPQ_PERIPH);

  dxdt_AQ_GUT      = -ka_AQ*AQ_GUT;
  dxdt_DEAQ_PLASMA =  ka_AQ*AQ_GUT/Vc_DEAQ - ke_DEAQ*DEAQ_PLASMA;

  dxdt_PQ_GUT   = -ka_PQ*PQ_GUT;
  dxdt_PQ_PLASMA=  ka_PQ*PQ_GUT/Vc_PQ - ke_PQ*PQ_PLASMA;

  double C_DHA=(DHA_PLASMA>0)?DHA_PLASMA*1000.0:0.0;
  double Clum_ng=(LUM_CENTRAL>0)?LUM_CENTRAL*1000.0:0.0;
  double Cppq_ng=(PPQ_CENTRAL>0)?PPQ_CENTRAL*1000.0:0.0;
  double Cdeaq_ng=(DEAQ_PLASMA>0)?DEAQ_PLASMA*1000.0:0.0;
  double Cpq_ng=(PQ_PLASMA>0)?PQ_PLASMA*1000.0:0.0;

  double EC50_AS_eff=EC50_AS*(1.0+(Kres_shift-1.0)*K13_RES);
  double E_AS =use_AS *Emax_AS *pow(C_DHA,hill_AS)  /(pow(EC50_AS_eff,hill_AS) +pow(C_DHA,hill_AS) +1e-9);
  double E_LUM=use_LUM*Emax_LUM*pow(Clum_ng,hill_LUM)/(pow(EC50_LUM,hill_LUM)+pow(Clum_ng,hill_LUM)+1e-9);
  double E_PPQ=use_PPQ*Emax_PPQ*pow(Cppq_ng,hill_PPQ)/(pow(EC50_PPQ,hill_PPQ)+pow(Cppq_ng,hill_PPQ)+1e-9);
  double E_AQ =use_AQ *Emax_AQ *pow(Cdeaq_ng,hill_AQ)/(pow(EC50_AQ,hill_AQ) +pow(Cdeaq_ng,hill_AQ)+1e-9);
  double E_partner=1.0-(1.0-E_LUM)*(1.0-E_PPQ)*(1.0-E_AQ);
  double E_PQ_liver=use_PQ*Emax_PQ_liver*pow(Cpq_ng,hill_PQ)/(pow(EC50_PQ,hill_PQ)+pow(Cpq_ng,hill_PQ)+1e-9);
  double E_PQ_gam  =use_PQ*Emax_PQ_gam  *pow(Cpq_ng,hill_PQ)/(pow(EC50_PQ,hill_PQ)+pow(Cpq_ng,hill_PQ)+1e-9);

  double RU=(RBC_U>0)?RBC_U:0.0;
  double R=(PRBC_RING>0)?PRBC_RING:0.0;
  double T=(PRBC_TROPH>0)?PRBC_TROPH:0.0;
  double S=(PRBC_SCHIZONT>0)?PRBC_SCHIZONT:0.0;
  double LP=(LIVER_PARASITE>0)?LIVER_PARASITE:0.0;
  double IMM=(IMMUNITY>0)?IMMUNITY:0.0;

  double total_para=R+T+S;
  double inv_eff=inv_eff0*(1.0-immune_block*IMM);
  double liver_flux=k_liver_egress*LP;
  double relapse_source=relapse_on*k_relapse*(LIVER_PARASITE>1.0?0.0:1.0)*5.0e3;
  double rupture_flux=k_SR*S;
  double new_infections=(liver_flux+rupture_flux*burst)*inv_eff*(RU/RBC_base);

  double kill_ring=E_AS;
  double kill_troph=E_AS*0.6+E_partner;
  double kill_schiz=E_partner*0.8;

  dxdt_LIVER_PARASITE=relapse_source-liver_flux-E_PQ_liver*LP;
  dxdt_RBC_U = k_RBC_sen*(RBC_base-RU) + k_ery_comp*(HB_base-HB>0?(HB_base-HB):0.0)*1.0e5/(1.0+k_dysery*total_para) - new_infections;
  dxdt_PRBC_RING    = new_infections - k_RT*R - kill_ring*R;
  dxdt_PRBC_TROPH   = k_RT*R - k_TS*T - kill_troph*T - gam_commit*k_RT*R;
  dxdt_PRBC_SCHIZONT= k_TS*T - k_SR*S - kill_schiz*S;
  dxdt_GAMETOCYTE   = gam_commit*k_RT*R - k_gam_mat*GAMETOCYTE - E_PQ_gam*GAMETOCYTE;

  double hemolysis=hb_per_rupture*rupture_flux + (G6PD_deficient>0.5? k_PQ_hemolysis*Cpq_ng*(HB/HB_base):0.0);
  dxdt_HB = -hemolysis + 0.0015*(HB_base-HB)/(1.0+k_dysery*total_para);

  double fever_drive=Tmax_rise*rupture_flux*burst/(rupture_flux*burst+K_fever+1e-6);
  dxdt_TEMP=(Tbase+fever_drive-TEMP)*k_temp_decay;
  dxdt_IMMUNITY=k_immune_gain*total_para/(total_para+K_immune)*(1.0-IMM) - k_immune_decay*IMM;

$TABLE
  double C_DHA_out =(DHA_PLASMA>0)?DHA_PLASMA*1000.0:0.0;
  double C_LUM_out =(LUM_CENTRAL>0)?LUM_CENTRAL*1000.0:0.0;
  double C_PPQ_out =(PPQ_CENTRAL>0)?PPQ_CENTRAL*1000.0:0.0;
  double C_DEAQ_out=(DEAQ_PLASMA>0)?DEAQ_PLASMA*1000.0:0.0;
  double C_PQ_out  =(PQ_PLASMA>0)?PQ_PLASMA*1000.0:0.0;
  double PERIPH_PARASITEMIA=PRBC_RING+0.10*PRBC_TROPH;
  double TOTAL_PARASITEMIA=PRBC_RING+PRBC_TROPH+PRBC_SCHIZONT;
  double log10_para=(PERIPH_PARASITEMIA>0.05)?log10(PERIPH_PARASITEMIA):-1.301;
  double Hb_gdl=HB; double Temp_C=TEMP; double Gam_uL=GAMETOCYTE;
  double Immunity_pct=IMMUNITY*100.0; double Liver_burden=LIVER_PARASITE;

$CAPTURE
  C_DHA_out C_LUM_out C_PPQ_out C_DEAQ_out C_PQ_out
  PERIPH_PARASITEMIA TOTAL_PARASITEMIA log10_para
  Hb_gdl Temp_C Gam_uL Immunity_pct Liver_burden
'

mod_base <- mrgsolve::mcode("malaria_shiny", mal_code, quiet = TRUE)

# ─────────────────────────────────────────────────────────────
# Dosing helper
# ─────────────────────────────────────────────────────────────
build_ev <- function(regimen, t_start) {
  switch(regimen,
    "AL (Artemether-Lumefantrine)" = c(
      ev(cmt="AS_GUT",  amt=100, ii=12, addl=5, time=t_start),
      ev(cmt="LUM_GUT", amt=480, ii=12, addl=5, time=t_start)),
    "ASAQ (Artesunate-Amodiaquine)" = c(
      ev(cmt="AS_GUT", amt=200, ii=24, addl=2, time=t_start),
      ev(cmt="AQ_GUT", amt=600, ii=24, addl=2, time=t_start)),
    "DP (Dihydroartemisinin-Piperaquine)" = c(
      ev(cmt="AS_GUT",  amt=240, ii=24, addl=2, time=t_start),
      ev(cmt="PPQ_GUT", amt=960, ii=24, addl=2, time=t_start)),
    "IV Artesunate (Severe Malaria)" = c(
      ev(cmt="AS_PLASMA", amt=2.4*70/30, ii=12, addl=2, time=t_start),
      ev(cmt="AS_GUT",  amt=100, ii=12, addl=5, time=t_start+72),
      ev(cmt="LUM_GUT", amt=480, ii=12, addl=5, time=t_start+72)),
    "No Treatment" = NULL
  )
}

use_flags <- function(regimen) {
  list(
    use_AS  = as.numeric(regimen %in% c("AL (Artemether-Lumefantrine)","ASAQ (Artesunate-Amodiaquine)",
                                         "DP (Dihydroartemisinin-Piperaquine)","IV Artesunate (Severe Malaria)")),
    use_LUM = as.numeric(regimen %in% c("AL (Artemether-Lumefantrine)","IV Artesunate (Severe Malaria)")),
    use_PPQ = as.numeric(regimen == "DP (Dihydroartemisinin-Piperaquine)"),
    use_AQ  = as.numeric(regimen == "ASAQ (Artesunate-Amodiaquine)")
  )
}

run_sim <- function(sim_days, regimen, species="P. falciparum", immunity0=0.05,
                     baseline_para=10, k13_res=FALSE, radical_cure=FALSE,
                     g6pd_def=FALSE, severe_init=FALSE) {

  pars <- use_flags(regimen)
  pars$K13_RES <- as.numeric(k13_res)
  pars$use_PQ  <- as.numeric(radical_cure)
  pars$G6PD_deficient <- as.numeric(g6pd_def)
  pars$relapse_on <- as.numeric(species == "P. vivax")

  m <- mod_base %>% param(pars) %>%
    init(IMMUNITY = immunity0,
         PRBC_RING = baseline_para * 0.5, PRBC_TROPH = baseline_para * 0.3,
         PRBC_SCHIZONT = baseline_para * 0.2,
         HB = if (severe_init) 9.0 else 13.5)

  evts <- build_ev(regimen, 0)
  if (radical_cure) evts <- c(evts, ev(cmt="PQ_GUT", amt=30, ii=24, addl=13, time=0))

  out <- if (is.null(evts)) mrgsim(m, end=sim_days*24, delta=1, obsonly=TRUE)
         else mrgsim(m, events=evts, end=sim_days*24, delta=1, obsonly=TRUE)
  as_tibble(out) %>% mutate(day = time/24)
}

# ─────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "green",
  dashboardHeader(title = "Malaria QSP Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("① Patient/Infection", tabName = "tab_patient", icon = icon("user-injured")),
      menuItem("② Drug PK",           tabName = "tab_pk",      icon = icon("pills")),
      menuItem("③ Parasite Dynamics", tabName = "tab_para",    icon = icon("bacteria")),
      menuItem("④ Clinical Endpoints",tabName = "tab_clin",    icon = icon("notes-medical")),
      menuItem("⑤ Scenario Compare",  tabName = "tab_compare", icon = icon("chart-bar")),
      menuItem("⑥ Biomarkers",        tabName = "tab_bio",     icon = icon("droplet")),
      menuItem("⑦ Resistance (K13)",  tabName = "tab_resist",  icon = icon("dna"))
    ),
    hr(),
    h5("Simulation Settings", style="padding-left:15px;color:#ccc"),
    selectInput("species", "Plasmodium Species", choices = c("P. falciparum","P. vivax"), selected="P. falciparum"),
    selectInput("regimen", "Treatment Regimen",
                choices = c("No Treatment","AL (Artemether-Lumefantrine)","ASAQ (Artesunate-Amodiaquine)",
                            "DP (Dihydroartemisinin-Piperaquine)","IV Artesunate (Severe Malaria)"),
                selected = "AL (Artemether-Lumefantrine)"),
    sliderInput("sim_days", "Simulation Duration (days)", min=7, max=200, value=40, step=1),
    sliderInput("immunity0", "Baseline Naturally-Acquired Immunity", min=0, max=0.9, value=0.05, step=0.05),
    sliderInput("baseline_para", "Baseline Parasitemia (/µL)", min=1, max=2000, value=10, step=1),
    checkboxInput("k13_res", "K13 Propeller Mutant (Artemisinin Resistance)", FALSE),
    checkboxInput("radical_cure", "Add 14-day Primaquine Radical Cure", FALSE),
    checkboxInput("g6pd_def", "G6PD Deficient (Hemolysis Risk)", FALSE),
    checkboxInput("severe_init", "Severe Malaria at Baseline (Hb 9 g/dL)", FALSE),
    actionButton("run_btn", "▶ Run Simulation", class = "btn-success btn-block")
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .box-header { background-color: #1b6b3a !important; }
      .content-wrapper, .right-side { background-color: #f7fbf8; }
    "))),
    tabItems(

      # ── Tab 1: Patient / Infection Profile ─────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          valueBoxOutput("vbox_species", width=3),
          valueBoxOutput("vbox_para0",   width=3),
          valueBoxOutput("vbox_immune",  width=3),
          valueBoxOutput("vbox_regimen", width=3)
        ),
        fluidRow(
          box(title="Infection Summary", width=6, solidHeader=TRUE, status="success",
              tableOutput("tbl_patient")),
          box(title="Mechanistic Basis (from Mechanistic Map)", width=6, solidHeader=TRUE, status="success",
              p("Sporozoite inoculation → hepatocyte invasion → exo-erythrocytic schizogony →",
                "merozoite release → 48h erythrocytic cycle (ring→trophozoite→schizont) →",
                "sequestration (PfEMP1-ICAM1/CD36/EPCR) → cytokine-driven fever & anemia."),
              p(strong("P. vivax / P. ovale:"), " form dormant hepatic hypnozoites causing relapse",
                " weeks-to-months after apparent cure unless 8-aminoquinoline radical cure is given."))
        ),
        fluidRow(
          box(title="WHO Severe Malaria Criteria Reference", width=12, solidHeader=TRUE, status="warning",
              tableOutput("tbl_severe_criteria"))
        )
      ),

      # ── Tab 2: Drug PK ──────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title="Dihydroartemisinin (Active Metabolite) Plasma Concentration", width=6,
              solidHeader=TRUE, status="success", plotlyOutput("plot_dha", height="320px")),
          box(title="Partner Drug Plasma Concentration (Post-Treatment Prophylaxis)", width=6,
              solidHeader=TRUE, status="success", plotlyOutput("plot_partner", height="320px"))
        ),
        fluidRow(
          box(title="Primaquine Plasma Concentration (if Radical Cure Selected)", width=12,
              solidHeader=TRUE, status="success", plotlyOutput("plot_pq", height="300px"))
        )
      ),

      # ── Tab 3: Parasite Dynamics ────────────────────────────
      tabItem(tabName = "tab_para",
        fluidRow(
          box(title="Peripheral Parasitemia (log10 parasites/µL)", width=12, solidHeader=TRUE,
              status="success", plotlyOutput("plot_parasitemia", height="380px"))
        ),
        fluidRow(
          box(title="Parasite Stage Composition", width=6, solidHeader=TRUE,
              plotlyOutput("plot_stages", height="320px")),
          box(title="Gametocyte Carriage (Transmission Potential)", width=6, solidHeader=TRUE,
              plotlyOutput("plot_gam", height="320px"))
        )
      ),

      # ── Tab 4: Clinical Endpoints ───────────────────────────
      tabItem(tabName = "tab_clin",
        fluidRow(
          valueBoxOutput("vbox_pct",  width=4),
          valueBoxOutput("vbox_fct",  width=4),
          valueBoxOutput("vbox_acpr", width=4)
        ),
        fluidRow(
          box(title="WWARN-Style Clinical/Parasitological Endpoints", width=12, solidHeader=TRUE,
              status="success", DTOutput("tbl_endpoints"))
        )
      ),

      # ── Tab 5: Scenario Comparison ──────────────────────────
      tabItem(tabName = "tab_compare",
        fluidRow(
          box(title="Regimen Comparison — Parasite Clearance", width=12, solidHeader=TRUE,
              status="success", plotlyOutput("plot_compare", height="420px"))
        ),
        fluidRow(
          box(title="Summary Table (All Regimens, Same Baseline)", width=12, solidHeader=TRUE,
              DTOutput("tbl_compare"))
        )
      ),

      # ── Tab 6: Biomarkers ───────────────────────────────────
      tabItem(tabName = "tab_bio",
        fluidRow(
          box(title="Hemoglobin Trajectory (Malarial Anemia)", width=6, solidHeader=TRUE,
              status="success", plotlyOutput("plot_hb", height="320px")),
          box(title="Body Temperature (Paroxysmal Fever)", width=6, solidHeader=TRUE,
              status="success", plotlyOutput("plot_temp", height="320px"))
        ),
        fluidRow(
          box(title="Naturally-Acquired Immunity Index", width=6, solidHeader=TRUE,
              plotlyOutput("plot_immunity", height="300px")),
          box(title="Liver-Stage / Hypnozoite Burden (P. vivax)", width=6, solidHeader=TRUE,
              plotlyOutput("plot_liver", height="300px"))
        )
      ),

      # ── Tab 7: Resistance ───────────────────────────────────
      tabItem(tabName = "tab_resist",
        fluidRow(
          box(title="K13 Wild-Type vs Mutant — Parasite Clearance", width=12, solidHeader=TRUE,
              status="success", plotlyOutput("plot_resist", height="400px"))
        ),
        fluidRow(
          box(title="Resistance Mechanism Reference", width=12, solidHeader=TRUE, status="warning",
              tableOutput("tbl_resist_ref"))
        )
      )
    )
  )
)

# ─────────────────────────────────────────────────────────────
# Server
# ─────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  sim_data <- eventReactive(input$run_btn, {
    run_sim(input$sim_days, input$regimen, input$species, input$immunity0,
            input$baseline_para, input$k13_res, input$radical_cure,
            input$g6pd_def, input$severe_init)
  }, ignoreNULL = FALSE)

  compare_data <- eventReactive(input$run_btn, {
    regimens <- c("No Treatment","AL (Artemether-Lumefantrine)","ASAQ (Artesunate-Amodiaquine)",
                  "DP (Dihydroartemisinin-Piperaquine)")
    bind_rows(lapply(regimens, function(r) {
      run_sim(min(input$sim_days,40), r, input$species, input$immunity0, input$baseline_para,
              input$k13_res, FALSE, FALSE, input$severe_init) %>% mutate(regimen = r)
    }))
  }, ignoreNULL = FALSE)

  resist_data <- eventReactive(input$run_btn, {
    bind_rows(
      run_sim(14, "AL (Artemether-Lumefantrine)", input$species, input$immunity0, input$baseline_para,
              FALSE, FALSE, FALSE, input$severe_init) %>% mutate(K13 = "Wild-Type"),
      run_sim(14, "AL (Artemether-Lumefantrine)", input$species, input$immunity0, input$baseline_para,
              TRUE, FALSE, FALSE, input$severe_init) %>% mutate(K13 = "C580Y Mutant")
    )
  }, ignoreNULL = FALSE)

  # ── Tab 1 outputs ──
  output$vbox_species <- renderValueBox(valueBox(input$species, "Species", icon=icon("microscope"), color="green"))
  output$vbox_para0   <- renderValueBox(valueBox(input$baseline_para, "Baseline Parasitemia (/µL)", icon=icon("bacteria"), color="orange"))
  output$vbox_immune  <- renderValueBox(valueBox(paste0(round(input$immunity0*100),"%"), "Acquired Immunity", icon=icon("shield-halved"), color="blue"))
  output$vbox_regimen <- renderValueBox(valueBox(strsplit(input$regimen," ")[[1]][1], "Regimen", icon=icon("prescription-bottle"), color="purple"))

  output$tbl_patient <- renderTable({
    data.frame(Parameter = c("Species","Regimen","Baseline Parasitemia","Baseline Immunity","G6PD Status","Severe at Baseline"),
               Value = c(input$species, input$regimen, paste(input$baseline_para,"/µL"),
                         paste0(round(input$immunity0*100),"%"),
                         ifelse(input$g6pd_def,"Deficient","Normal"),
                         ifelse(input$severe_init,"Yes (Hb 9 g/dL)","No")))
  })

  output$tbl_severe_criteria <- renderTable({
    data.frame(
      Criterion = c("Cerebral malaria","Severe anemia","Acute kidney injury","Pulmonary edema/ARDS",
                    "Hypoglycemia","Hyperparasitemia","Acidosis","Shock"),
      Threshold = c("Coma (Blantyre <=2)","Hb <5 g/dL / Hct <15%","Creatinine >3 mg/dL","SpO2 <92%",
                    "Glucose <40 mg/dL","Parasitemia >10%","Lactate >5 mmol/L / pH<7.25","SBP <80 mmHg")
    )
  })

  # ── Tab 2: PK ──
  output$plot_dha <- renderPlotly({
    d <- sim_data() %>% filter(time <= 96)
    p <- ggplot(d, aes(x=time, y=C_DHA_out)) + geom_line(color="#1b6b3a", linewidth=0.9) +
      labs(x="Time (h)", y="DHA (ng/mL)") + theme_bw()
    ggplotly(p)
  })
  output$plot_partner <- renderPlotly({
    d <- sim_data()
    d$Cpartner <- pmax(d$C_LUM_out, d$C_PPQ_out, d$C_DEAQ_out)
    p <- ggplot(d, aes(x=day, y=Cpartner)) + geom_line(color="#984EA3", linewidth=0.9) +
      labs(x="Time (days)", y="Partner Drug (ng/mL)") + theme_bw()
    ggplotly(p)
  })
  output$plot_pq <- renderPlotly({
    d <- sim_data() %>% filter(day <= 20)
    p <- ggplot(d, aes(x=day, y=C_PQ_out)) + geom_line(color="#008888", linewidth=0.9) +
      labs(x="Time (days)", y="Primaquine (ng/mL)") + theme_bw()
    ggplotly(p)
  })

  # ── Tab 3: Parasite Dynamics ──
  output$plot_parasitemia <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=day, y=log10_para)) + geom_line(color="#E41A1C", linewidth=1.0) +
      geom_hline(yintercept=log10(10), linetype="dashed", color="grey40") +
      labs(x="Time (days)", y="log10 Parasitemia (/µL)") + theme_bw()
    ggplotly(p)
  })
  output$plot_stages <- renderPlotly({
    d <- sim_data() %>% select(day, PRBC_RING, PRBC_TROPH, PRBC_SCHIZONT) %>%
      pivot_longer(-day, names_to="Stage", values_to="Count")
    p <- ggplot(d %>% filter(day<=10), aes(x=day, y=Count, color=Stage)) + geom_line(linewidth=0.8) +
      labs(x="Time (days)", y="Parasites/µL") + theme_bw() + theme(legend.position="bottom")
    ggplotly(p)
  })
  output$plot_gam <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=day, y=Gam_uL)) + geom_line(color="#FF7F00", linewidth=0.9) +
      labs(x="Time (days)", y="Mature Gametocytes (/µL)") + theme_bw()
    ggplotly(p)
  })

  # ── Tab 4: Clinical Endpoints ──
  endpoint_calc <- reactive({
    d <- sim_data()
    below <- d$time[d$PERIPH_PARASITEMIA < 10 & d$time > 12]
    pct <- if (length(below) > 0) round(min(below)) else NA
    fct <- { ft <- d$time[d$Temp_C < 37.5 & d$time > 6]; if (length(ft)>0) round(min(ft)) else NA }
    d28 <- any(d$day > 28 & d$PERIPH_PARASITEMIA >= 10, na.rm=TRUE)
    list(pct=pct, fct=fct, acpr=!isTRUE(d28))
  })
  output$vbox_pct  <- renderValueBox(valueBox(paste(endpoint_calc()$pct,"h"), "Parasite Clearance Time", icon=icon("clock"), color="green"))
  output$vbox_fct  <- renderValueBox(valueBox(paste(endpoint_calc()$fct,"h"), "Fever Clearance Time", icon=icon("temperature-low"), color="orange"))
  output$vbox_acpr <- renderValueBox(valueBox(ifelse(endpoint_calc()$acpr,"ACPR","Failure"), "Day-28 Outcome", icon=icon("check"), color=ifelse(endpoint_calc()$acpr,"green","red")))

  output$tbl_endpoints <- renderDT({
    d <- sim_data()
    datatable(data.frame(
      Endpoint = c("Peak Parasitemia (/µL)","Parasite Clearance Time (h)","Fever Clearance Time (h)",
                   "Hemoglobin Nadir (g/dL)","Peak Temperature (C)","Day-3 Microscopy Positive"),
      Value = c(round(max(d$PERIPH_PARASITEMIA)), endpoint_calc()$pct, endpoint_calc()$fct,
                round(min(d$Hb_gdl),1), round(max(d$Temp_C),1),
                ifelse(any(abs(d$time-72)<1 & d$PERIPH_PARASITEMIA>=10),"Yes","No"))
    ), options=list(dom='t'))
  })

  # ── Tab 5: Scenario Comparison ──
  output$plot_compare <- renderPlotly({
    d <- compare_data()
    p <- ggplot(d, aes(x=day, y=log10_para, color=regimen)) + geom_line(linewidth=0.9) +
      geom_hline(yintercept=log10(10), linetype="dashed", color="grey40") +
      labs(x="Time (days)", y="log10 Parasitemia (/µL)", color=NULL) + theme_bw() +
      theme(legend.position="bottom")
    ggplotly(p)
  })
  output$tbl_compare <- renderDT({
    d <- compare_data() %>% group_by(regimen) %>%
      summarise(Peak_Parasitemia = round(max(PERIPH_PARASITEMIA)),
                Hb_Nadir = round(min(Hb_gdl),1),
                Day10_Parasitemia = round(PERIPH_PARASITEMIA[which.min(abs(day-10))]), .groups="drop")
    datatable(d, options=list(dom='t'))
  })

  # ── Tab 6: Biomarkers ──
  output$plot_hb <- renderPlotly({
    p <- ggplot(sim_data(), aes(x=day, y=Hb_gdl)) + geom_line(color="#CC0000", linewidth=0.9) +
      labs(x="Time (days)", y="Hemoglobin (g/dL)") + theme_bw()
    ggplotly(p)
  })
  output$plot_temp <- renderPlotly({
    p <- ggplot(sim_data(), aes(x=day, y=Temp_C)) + geom_line(color="#CC6600", linewidth=0.9) +
      geom_hline(yintercept=37.5, linetype="dashed") +
      labs(x="Time (days)", y="Temperature (C)") + theme_bw()
    ggplotly(p)
  })
  output$plot_immunity <- renderPlotly({
    p <- ggplot(sim_data(), aes(x=day, y=Immunity_pct)) + geom_line(color="#0055AA", linewidth=0.9) +
      labs(x="Time (days)", y="Immunity Index (%)") + theme_bw()
    ggplotly(p)
  })
  output$plot_liver <- renderPlotly({
    p <- ggplot(sim_data(), aes(x=day, y=Liver_burden)) + geom_line(color="#008888", linewidth=0.9) +
      labs(x="Time (days)", y="Liver/Hypnozoite Burden") + theme_bw()
    ggplotly(p)
  })

  # ── Tab 7: Resistance ──
  output$plot_resist <- renderPlotly({
    d <- resist_data()
    p <- ggplot(d, aes(x=day, y=log10_para, color=K13)) + geom_line(linewidth=1.0) +
      geom_hline(yintercept=log10(10), linetype="dashed", color="grey40") +
      labs(x="Time (days)", y="log10 Parasitemia (/µL)", color="K13 Genotype") + theme_bw() +
      theme(legend.position="bottom")
    ggplotly(p)
  })
  output$tbl_resist_ref <- renderTable({
    data.frame(
      Marker = c("K13 C580Y / R561H","pfmdr1 amplification","pfcrt K76T","plasmepsin 2/3 amplification"),
      Association = c("Delayed artemisinin ring-stage clearance (Ashley 2014 NEJM)",
                       "Reduced lumefantrine/mefloquine susceptibility",
                       "Chloroquine/amodiaquine resistance",
                       "Piperaquine treatment failure (Amato 2017 Lancet ID)")
    )
  })
}

shinyApp(ui, server)
