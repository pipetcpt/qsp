## ============================================================
## Pseudogout (CPPD Crystal Deposition Disease) — Shiny Dashboard
## 7 Tabs: Patient Profile · Drug PK · Inflammation · Crystal
##         Dynamics · Clinical Endpoints · Scenario Comparison ·
##         Biomarkers & Sensitivity
## ============================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(shinydashboard)
library(plotly)

# ---- Embed mrgsolve model code (same as cppd_mrgsolve_model.R) ----
MODEL_CODE <- '
$PROB CPPD QSP Model (Shiny)
$PARAM
DOSE_COLCH=0, DOSE_INDO=0, DOSE_PRED=0, DOSE_ANA=0,
F_COLCH=0.45, Ka_COLCH=1.2, CL_COLCH=31.0, Vd_COLCH=1470,
F_INDO=0.98,  Ka_INDO=1.8,  CL_INDO=12.0,  Vd_INDO=24.0,
F_PRED=0.82,  Ka_PRED=1.5,  CL_PRED=6.6,   Vd_PRED=33.0,
EC50_GR=5.0,  Emax_GR=0.95,
F_ANA=0.95,   Ka_ANA=0.4,   CL_ANA=6.3,    Vd_ANA=7.0,
PPi_base=5.0, kENPP1=0.15,  kANKH=0.08,    kTNAP=0.20,
Ca_art=2.2,   Mg_art=0.4,
PPi_th=6.0,   kNuc=0.002,   kGrow=0.05,    kShed_base=0.003,
kDissolve=0.002, CrystMax=100,
kNLRP3_on=0.5, EC50_cryst_nlrp3=5.0, kNLRP3_off=0.3,
kIL1b_syn=2.0, kIL1b_deg=0.4, IL1Ra_endo=0.6,
kNeut_rec=8.0, EC50_neut=20, kNeut_deg=0.1, Neut_base=200,
kIL6_syn=1.5,  kIL6_deg=0.5,
kCRP_syn=3.0,  kCRP_deg=0.04, CRP_base=2.0,
kPGE2_syn=0.8, kPGE2_deg=1.2,
PainVAS_base=1.0, kPain_resp=0.005,
CI_init=95.0,  kCart_dam=0.0002, kCart_rep=0.00005, CI_min=20.0,
kLipox_syn=0.05, kLipox_deg=0.3,
EC50_COLCH_NLRP3=1.5, EC50_COLCH_NEUT=2.0, Emax_COLCH=0.85,
EC50_INDO_COX=0.3, Emax_INDO=0.90,
EC50_PRED_NFKB=5.0, Emax_PRED=0.88,
KD_ANA=0.5, Emax_ANA=0.92

$CMT
COLCH_depot COLCH_central INDO_depot INDO_central
PRED_depot  PRED_central  ANA_depot  ANA_central
PPi_ext Cryst_cart Cryst_SF
NLRP3_act IL1b Neutrophil IL6 PGE2 CRP PainVAS CartInteg LipoxA4

$MAIN
double Cp_COLCH = COLCH_central / Vd_COLCH * 1000;
double Cp_INDO  = INDO_central  / Vd_INDO;
double Cp_PRED  = PRED_central  / Vd_PRED * 1000;
double Cp_ANA   = ANA_central   / Vd_ANA;
double GR_occ   = Emax_GR * Cp_PRED / (EC50_GR + Cp_PRED);
double E_COLCH_NLRP3 = Emax_COLCH * Cp_COLCH / (EC50_COLCH_NLRP3 + Cp_COLCH);
double E_COLCH_NEUT  = Emax_COLCH * Cp_COLCH / (EC50_COLCH_NEUT  + Cp_COLCH);
double E_INDO_COX    = Emax_INDO  * Cp_INDO  / (EC50_INDO_COX + Cp_INDO);
double E_PRED_NFKB   = Emax_PRED  * GR_occ;
double E_ANA         = Emax_ANA   * Cp_ANA   / (KD_ANA + Cp_ANA);
double E_IL1b_total  = 1 - (1-E_COLCH_NLRP3)*(1-E_PRED_NFKB)*(1-E_ANA*0.7);
double kShed_eff     = kShed_base * (1 + 0.01*Cryst_cart);
double Neut_norm     = fmax(Neutrophil, 0);
double cryst_SF_eff  = fmax(Cryst_SF, 0);
double IL8_proxy     = fmax(IL1b,0)*0.5 + fmax(IL6,0)*0.2;
double IL8_recruit   = kNeut_rec * IL8_proxy / (EC50_neut + IL8_proxy);
double LipoxA4_eff   = fmax(LipoxA4, 0);
double resolve_factor= 1 + 0.5*LipoxA4_eff/(50+LipoxA4_eff);
double PainVAS_target= PainVAS_base + kPain_resp*fmax(IL1b,0) + 0.0008*fmax(PGE2,0);
double kDam_eff      = kCart_dam*(1+0.01*fmax(IL1b,0)+0.005*fmax(Cryst_cart,0));

$ODE
dxdt_COLCH_depot   = -Ka_COLCH*COLCH_depot;
dxdt_COLCH_central =  Ka_COLCH*F_COLCH*COLCH_depot - (CL_COLCH/Vd_COLCH)*COLCH_central;
dxdt_INDO_depot    = -Ka_INDO*INDO_depot;
dxdt_INDO_central  =  Ka_INDO*F_INDO*INDO_depot - (CL_INDO/Vd_INDO)*INDO_central;
dxdt_PRED_depot    = -Ka_PRED*PRED_depot;
dxdt_PRED_central  =  Ka_PRED*F_PRED*PRED_depot - (CL_PRED/Vd_PRED)*PRED_central;
dxdt_ANA_depot     = -Ka_ANA*ANA_depot;
dxdt_ANA_central   =  Ka_ANA*F_ANA*ANA_depot - (CL_ANA/Vd_ANA)*ANA_central;
dxdt_PPi_ext  = kENPP1*(PPi_base+3.0) + kANKH*2.0 - kTNAP*PPi_ext;
double PPi_excess = fmax(PPi_ext - PPi_th, 0);
double kNuc_eff   = kNuc*PPi_excess*Ca_art/(1+Mg_art/0.5);
dxdt_Cryst_cart= kNuc_eff*(CrystMax-Cryst_cart)+kGrow*Cryst_cart*PPi_excess/(PPi_excess+2.0)-kShed_eff*Cryst_cart;
dxdt_Cryst_SF  = kShed_eff*Cryst_cart - kDissolve*Cryst_SF - 0.02*Neut_norm*Cryst_SF/(50+Cryst_SF);
double nlrp3_stim = kNLRP3_on*cryst_SF_eff/(EC50_cryst_nlrp3+cryst_SF_eff);
double nlrp3_inhib= (E_COLCH_NLRP3+GR_occ*0.3)*NLRP3_act;
dxdt_NLRP3_act = nlrp3_stim*(1-NLRP3_act) - kNLRP3_off*NLRP3_act - nlrp3_inhib;
double IL1b_prod = kIL1b_syn*NLRP3_act*(1-IL1Ra_endo)*E_IL1b_total;
dxdt_IL1b     = IL1b_prod - kIL1b_deg*IL1b;
double neut_ingress = IL8_recruit*(1-E_COLCH_NEUT)/resolve_factor;
dxdt_Neutrophil = neut_ingress - kNeut_deg*Neut_norm;
double IL6_prod = kIL6_syn*(1+0.05*fmax(IL1b,0))*(1-E_PRED_NFKB*0.7);
dxdt_IL6      = IL6_prod - kIL6_deg*IL6;
double pge2_prod = kPGE2_syn*(1+0.02*fmax(IL1b,0))*(1-E_INDO_COX)*(1-E_PRED_NFKB*0.5);
dxdt_PGE2     = pge2_prod - kPGE2_deg*PGE2;
double crp_prod = kCRP_syn*fmax(IL6,0)/(5+fmax(IL6,0));
dxdt_CRP      = crp_prod + 0.04*CRP_base - kCRP_deg*CRP;
dxdt_PainVAS  = 0.3*(PainVAS_target - PainVAS);
dxdt_CartInteg= kCart_rep*(100-CartInteg) - kDam_eff*CartInteg;
double lipox_stim = kLipox_syn*fmax(Neut_norm,0)/(200+fmax(Neut_norm,0));
dxdt_LipoxA4  = lipox_stim*100 - kLipox_deg*LipoxA4;

$CAPTURE
Cp_COLCH Cp_INDO Cp_PRED Cp_ANA GR_occ
E_COLCH_NLRP3 E_INDO_COX E_PRED_NFKB E_ANA
'

# ---- Compile model once at startup ----
mod <- mcode("CPPD_QSP_SHINY", MODEL_CODE, quiet = TRUE)

# ---- Default initial conditions ----
INIT_DEF <- list(
  COLCH_depot=0, COLCH_central=0, INDO_depot=0, INDO_central=0,
  PRED_depot=0,  PRED_central=0,  ANA_depot=0,  ANA_central=0,
  PPi_ext=8.0, Cryst_cart=25.0, Cryst_SF=0.0, NLRP3_act=0.0,
  IL1b=0.5, Neutrophil=200, IL6=3.0, PGE2=50, CRP=2.0,
  PainVAS=1.0, CartInteg=75.0, LipoxA4=20.0
)

# ---- Simulation helper ----
run_sim <- function(drug, colch_dose, indo_dose, pred_dose, ana_dose,
                    crystal_shed, sim_days, ca_art, mg_art, ppith) {
  tend <- sim_days * 24
  tg <- seq(0, tend, by = 2)

  ev_list <- list()
  if ("Colchicine" %in% drug && colch_dose > 0)
    ev_list[["colch"]] <- ev(time=seq(0, tend, by=12), cmt="COLCH_depot", amt=colch_dose)
  if ("Indomethacin" %in% drug && indo_dose > 0)
    ev_list[["indo"]]  <- ev(time=seq(0, tend, by=8),  cmt="INDO_depot",  amt=indo_dose)
  if ("Prednisolone" %in% drug && pred_dose > 0)
    ev_list[["pred"]]  <- ev(time=seq(0, tend, by=24), cmt="PRED_depot",  amt=pred_dose)
  if ("Anakinra" %in% drug && ana_dose > 0)
    ev_list[["ana"]]   <- ev(time=seq(0, tend, by=24), cmt="ANA_depot",   amt=ana_dose)

  evs <- ev(time=0, cmt="Cryst_SF", amt=crystal_shed)
  for (e in ev_list) evs <- evs + e

  out <- mod %>%
    init(INIT_DEF) %>%
    param(Ca_art=ca_art, Mg_art=mg_art, PPi_th=ppith) %>%
    mrgsim_e(evs, tgrid=tg, output="df") %>%
    mutate(time_d = time / 24)
  out
}

scenario_sim <- function(ca_art, mg_art, ppith) {
  tend <- 10 * 24
  tg   <- seq(0, tend, by=2)
  shed <- ev(time=0, cmt="Cryst_SF", amt=15)

  run_one <- function(add_evs, label) {
    evs <- shed
    if (!is.null(add_evs)) evs <- evs + add_evs
    mod %>% init(INIT_DEF) %>%
      param(Ca_art=ca_art, Mg_art=mg_art, PPi_th=ppith) %>%
      mrgsim_e(evs, tgrid=tg, output="df") %>%
      mutate(time_d=time/24, scenario=label)
  }

  bind_rows(
    run_one(NULL, "1. Untreated"),
    run_one(ev(time=seq(0,tend,by=12), cmt="COLCH_depot", amt=0.5), "2. Colchicine 0.5mg BID"),
    run_one(ev(time=seq(0,tend,by=8),  cmt="INDO_depot",  amt=50),  "3. Indomethacin 50mg TID"),
    run_one(ev(time=seq(0,tend,by=24), cmt="PRED_depot",  amt=30),  "4. Prednisolone 30mg/d"),
    run_one(ev(time=seq(0,tend,by=24), cmt="ANA_depot",   amt=100), "5. Anakinra 100mg/d"),
    run_one(ev(time=seq(0,tend,by=12), cmt="COLCH_depot", amt=0.5) +
            ev(time=seq(0,tend,by=8),  cmt="INDO_depot",  amt=50),  "6. Colch+Indo Combo")
  )
}

scen_pal <- c(
  "1. Untreated"           = "#e74c3c",
  "2. Colchicine 0.5mg BID"= "#2980b9",
  "3. Indomethacin 50mg TID"="#27ae60",
  "4. Prednisolone 30mg/d" = "#f39c12",
  "5. Anakinra 100mg/d"    = "#8e44ad",
  "6. Colch+Indo Combo"    = "#16a085"
)

# ---- UI ----
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "CPPD / Pseudogout QSP Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient Profile",   tabName="tab_profile",  icon=icon("user-md")),
      menuItem("2. Drug PK",           tabName="tab_pk",       icon=icon("pills")),
      menuItem("3. Inflammation",      tabName="tab_inflam",   icon=icon("fire")),
      menuItem("4. Crystal Dynamics",  tabName="tab_crystal",  icon=icon("gem")),
      menuItem("5. Clinical Endpoints",tabName="tab_clinical", icon=icon("chart-line")),
      menuItem("6. Scenario Comparison",tabName="tab_scenario",icon=icon("balance-scale")),
      menuItem("7. Biomarkers & Sensitivity",tabName="tab_bio",icon=icon("flask"))
    ),
    hr(),
    h5("Simulation Controls", style="color:white; padding-left:15px"),
    checkboxGroupInput("drug", "Active Drugs:",
                       choices=c("Colchicine","Indomethacin","Prednisolone","Anakinra"),
                       selected="Colchicine"),
    conditionalPanel("input.drug.indexOf('Colchicine')>=0",
      numericInput("colch_dose","Colchicine/dose (mg):", 0.5, 0.5, 1.2, 0.5)),
    conditionalPanel("input.drug.indexOf('Indomethacin')>=0",
      numericInput("indo_dose","Indomethacin/dose (mg):", 50, 25, 100, 25)),
    conditionalPanel("input.drug.indexOf('Prednisolone')>=0",
      numericInput("pred_dose","Prednisolone/d (mg):", 30, 5, 60, 5)),
    conditionalPanel("input.drug.indexOf('Anakinra')>=0",
      numericInput("ana_dose","Anakinra/d (mg):", 100, 50, 200, 50)),
    sliderInput("crystal_shed","Crystal Shed at t=0 (µg/mL):", 5, 40, 15, 1),
    sliderInput("sim_days","Simulation Duration (days):", 3, 30, 10, 1),
    sliderInput("ca_art","Articular Ca²⁺ (mM):", 1.5, 3.5, 2.2, 0.1),
    sliderInput("mg_art","Articular Mg²⁺ (mM):", 0.1, 1.0, 0.4, 0.05),
    sliderInput("ppith","PPi Nucleation Threshold (µM):", 4, 10, 6, 0.5),
    actionButton("run","▶ Run Simulation", class="btn-success btn-block")
  ),
  dashboardBody(
    tabItems(
      # ---- TAB 1: Patient Profile ----
      tabItem(tabName="tab_profile",
        fluidRow(
          box(width=12, title="CPPD Crystal Deposition Disease — Overview", status="primary", solidHeader=TRUE,
            fluidRow(
              column(6,
                h4("Epidemiology & Risk Factors"),
                tags$ul(
                  tags$li("Prevalence: ~4-7% of adults ≥60 yr; up to 35% of autopsies ≥80 yr"),
                  tags$li("Male = Female (contrast to gout); peak incidence 7th-8th decade"),
                  tags$li("Positively birefringent rhomboid crystals under compensated polarized light"),
                  tags$li("Associated: hyperparathyroidism, hypomagnesemia, hemochromatosis,",
                          "Wilson's disease, hypophosphatasia, familial CPPD (ANKH mutation)")
                ),
                h4("Key Pathophysiological Pathways"),
                tags$ul(
                  tags$li("ENPP1 → PPi generation; ANKH → PPi efflux; TNAP → PPi clearance"),
                  tags$li("When PPi > threshold + Ca²⁺↑ + Mg²⁺↓: CPPD nucleation in fibrocartilage"),
                  tags$li("Crystal shedding → SF → TLR2/4 signaling + NLRP3 inflammasome"),
                  tags$li("IL-1β → neutrophil recruitment → amplified inflammation"),
                  tags$li("MMP-3/13 + ADAMTS-4/5 → cartilage matrix degradation")
                )
              ),
              column(6,
                h4("Clinical Presentations"),
                tableOutput("profile_tbl"),
                h4("EULAR Recommendations (2018)"),
                tags$ul(
                  tags$li("Acute: Colchicine 0.5mg BID | NSAIDs | IA/systemic glucocorticoids"),
                  tags$li("Refractory: Anakinra 100mg/d SC × 3 days | Canakinumab 150mg SC"),
                  tags$li("Chronic prophylaxis: Colchicine 0.5mg/d (off-label) | Hydroxychloroquine"),
                  tags$li("No crystal-dissolving therapy currently available (unlike gout)")
                )
              )
            )
          )
        ),
        fluidRow(
          box(width=6, title="PPi Metabolism", status="info",
            tags$table(class="table table-bordered table-sm",
              tags$tr(tags$th("Enzyme/Transporter"),tags$th("Function"),tags$th("Direction")),
              tags$tr(tags$td("ENPP1/NPP1"),tags$td("ATP → AMP + PPi"),tags$td("PPi ↑")),
              tags$tr(tags$td("ANKH channel"),tags$td("PPi efflux (ICF→ECF)"),tags$td("PPi ↑")),
              tags$tr(tags$td("TNAP/ALPL"),tags$td("PPi → 2Pi (hydrolysis)"),tags$td("PPi ↓")),
              tags$tr(tags$td("NT5E/CD73"),tags$td("AMP → Adenosine"),tags$td("indirect")),
              tags$tr(tags$td("Mg²⁺"),tags$td("Inhibits CPPD nucleation"),tags$td("Crystal ↓"))
            )
          ),
          box(width=6, title="Drug Summary", status="warning",
            tableOutput("drug_tbl")
          )
        )
      ),

      # ---- TAB 2: Drug PK ----
      tabItem(tabName="tab_pk",
        fluidRow(
          box(width=6, title="Drug Plasma Concentrations", status="primary", plotlyOutput("pk_plot", height=350)),
          box(width=6, title="Receptor/Enzyme Occupancy", status="info",    plotlyOutput("occ_plot", height=350))
        ),
        fluidRow(
          box(width=12, title="PK Parameters Summary", status="primary",
            tableOutput("pk_params_tbl")
          )
        )
      ),

      # ---- TAB 3: Inflammation ----
      tabItem(tabName="tab_inflam",
        fluidRow(
          box(width=6, title="IL-1β Dynamics (pg/mL)",  status="danger",   plotlyOutput("il1b_plot", height=300)),
          box(width=6, title="NLRP3 Activation (0-1)",  status="warning",  plotlyOutput("nlrp3_plot", height=300))
        ),
        fluidRow(
          box(width=6, title="IL-6 (pg/mL)",            status="info",     plotlyOutput("il6_plot", height=300)),
          box(width=6, title="PGE2 in SF (pg/mL)",      status="warning",  plotlyOutput("pge2_plot", height=300))
        ),
        fluidRow(
          box(width=6, title="SF Neutrophils (cells/µL)",status="success",  plotlyOutput("neut_plot", height=300)),
          box(width=6, title="Lipoxin A4 (pro-resolving, pg/mL)", status="success", plotlyOutput("lipox_plot", height=300))
        )
      ),

      # ---- TAB 4: Crystal Dynamics ----
      tabItem(tabName="tab_crystal",
        fluidRow(
          box(width=6, title="Extracellular PPi (µM)", status="primary",   plotlyOutput("ppi_plot",   height=300)),
          box(width=6, title="CPPD in Cartilage (mg/g dry wt)", status="warning", plotlyOutput("cart_cryst_plot", height=300))
        ),
        fluidRow(
          box(width=6, title="CPPD Crystals in Synovial Fluid (µg/mL)", status="danger", plotlyOutput("sf_cryst_plot", height=300)),
          box(width=6, title="Crystal Dynamics Interpretation", status="info",
            tags$ul(
              tags$li("PPi threshold for nucleation: ~6 µM (adjust slider to modify)"),
              tags$li("Ca²⁺ elevation (hyperparathyroidism) accelerates crystal growth"),
              tags$li("Mg²⁺ deficiency removes inhibitory brake on nucleation"),
              tags$li("TNAP deficiency (hypophosphatasia, Fe/Cu toxicity) ↑ PPi"),
              tags$li("Crystal shedding into SF triggers acute attack via NLRP3"),
              tags$li("No current therapy reduces cartilage crystal burden directly")
            )
          )
        )
      ),

      # ---- TAB 5: Clinical Endpoints ----
      tabItem(tabName="tab_clinical",
        fluidRow(
          box(width=6, title="Pain VAS (0-10)", status="danger", plotlyOutput("pain_plot", height=300)),
          box(width=6, title="CRP (mg/L)",      status="warning",plotlyOutput("crp_plot",  height=300))
        ),
        fluidRow(
          box(width=6, title="Cartilage Integrity (%)", status="info",    plotlyOutput("ci_plot", height=300)),
          box(width=6, title="Key Time-Points Table",   status="primary", DTOutput("ktp_tbl"))
        )
      ),

      # ---- TAB 6: Scenario Comparison ----
      tabItem(tabName="tab_scenario",
        fluidRow(
          box(width=12, title="Treatment Scenario Comparison (all 6 regimens, fixed 10-day acute attack)",
              status="primary", solidHeader=TRUE,
            fluidRow(
              column(6, plotlyOutput("scen_il1b",  height=280)),
              column(6, plotlyOutput("scen_pain",  height=280))
            ),
            fluidRow(
              column(6, plotlyOutput("scen_crp",   height=280)),
              column(6, plotlyOutput("scen_nlrp3", height=280))
            )
          )
        ),
        fluidRow(
          box(width=12, title="Scenario Summary Table (Day 7)", status="info",
            DTOutput("scen_tbl")
          )
        )
      ),

      # ---- TAB 7: Biomarkers & Sensitivity ----
      tabItem(tabName="tab_bio",
        fluidRow(
          box(width=12, title="Biomarker Dashboard", status="primary", solidHeader=TRUE,
            fluidRow(
              column(4,
                h4("Clinical Biomarkers"),
                tableOutput("biomarker_tbl")
              ),
              column(8,
                plotlyOutput("biomarker_plot", height=350)
              )
            )
          )
        ),
        fluidRow(
          box(width=12, title="Sensitivity Analysis — Parameter Impact on Day-7 IL-1β & Pain VAS",
              status="warning",
            plotlyOutput("sens_plot", height=380)
          )
        )
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {

  sim_result <- eventReactive(input$run, {
    run_sim(
      drug         = input$drug,
      colch_dose   = input$colch_dose,
      indo_dose    = input$indo_dose,
      pred_dose    = input$pred_dose,
      ana_dose     = input$ana_dose,
      crystal_shed = input$crystal_shed,
      sim_days     = input$sim_days,
      ca_art       = input$ca_art,
      mg_art       = input$mg_art,
      ppith        = input$ppith
    )
  }, ignoreNULL = FALSE,
  ignoreInit  = FALSE)

  scen_result <- eventReactive(input$run, {
    scenario_sim(input$ca_art, input$mg_art, input$ppith)
  }, ignoreNULL = FALSE, ignoreInit = FALSE)

  # Helper: plotly line
  pl <- function(df, y, ylab, col="#2980b9") {
    plot_ly(df, x=~time_d, y=as.formula(paste0("~",y)),
            type="scatter", mode="lines",
            line=list(color=col, width=2)) %>%
      layout(xaxis=list(title="Time (days)"),
             yaxis=list(title=ylab))
  }

  # ---- TAB 1 ----
  output$profile_tbl <- renderTable({
    data.frame(
      Presentation = c("Acute pseudogout","Asymptomatic CPPD","Chronic CPP arthritis","Crowned dens syndrome"),
      Features     = c("Monoarthritis, SF WBC >25k, fever","Chondrocalcinosis on X-ray","OA-like symmetric arthritis","Neck pain, odontoid process crystals"),
      Treatment    = c("Colchicine/NSAID/GCS","Observation","MTX/HCQ (limited evidence)","Steroids/colchicine")
    )
  }, striped=TRUE, bordered=TRUE, hover=TRUE)

  output$drug_tbl <- renderTable({
    data.frame(
      Drug       = c("Colchicine","Indomethacin","Prednisolone","Anakinra","Canakinumab","HCQ","MTX"),
      Dose       = c("0.5mg BID","50mg TID","30mg/d→taper","100mg/d SC","150mg SC q3m","200-400mg/d","7.5-15mg/wk"),
      Target     = c("Tubulin/NLRP3","COX-1/2","GR/NF-κB","IL-1R","IL-1β mAb","Lysosomal pH","ATIC/DHFR"),
      `t½`       = c("26-31h","4.5h","2.5-4h","4-6h","26 days","40-50d","3-10h (MTX-PG ~days)")
    )
  }, striped=TRUE, bordered=TRUE, hover=TRUE)

  # ---- TAB 2: PK ----
  output$pk_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df) %>%
      add_lines(x=~time_d, y=~Cp_COLCH, name="Colchicine (ng/mL)", line=list(color="#2980b9")) %>%
      add_lines(x=~time_d, y=~Cp_PRED,  name="Prednisolone (ng/mL)", line=list(color="#e67e22")) %>%
      add_lines(x=~time_d, y=~Cp_ANA,   name="Anakinra (µg/mL)",    line=list(color="#8e44ad")) %>%
      add_lines(x=~time_d, y=~Cp_INDO,  name="Indomethacin (µg/mL)",line=list(color="#27ae60")) %>%
      layout(xaxis=list(title="Time (days)"), yaxis=list(title="Concentration"), legend=list(orientation="h"))
  })

  output$occ_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df) %>%
      add_lines(x=~time_d, y=~GR_occ,        name="GR Occupancy (Pred)", line=list(color="#e67e22")) %>%
      add_lines(x=~time_d, y=~E_COLCH_NLRP3, name="NLRP3 Inh (Colch)",  line=list(color="#2980b9")) %>%
      add_lines(x=~time_d, y=~E_INDO_COX,    name="COX Inh (Indo)",      line=list(color="#27ae60")) %>%
      add_lines(x=~time_d, y=~E_ANA,         name="IL-1R Block (Ana)",   line=list(color="#8e44ad")) %>%
      layout(xaxis=list(title="Time (days)"), yaxis=list(title="Fractional Effect (0-1)"),
             legend=list(orientation="h"))
  })

  output$pk_params_tbl <- renderTable({
    data.frame(
      Drug=c("Colchicine","Indomethacin","Prednisolone","Anakinra"),
      `F (%)`=c("45","98","82","95"),
      `Ka (/h)`=c(1.2,1.8,1.5,0.4),
      `Vd (L)`=c(1470,24,33,7),
      `CL (L/h)`=c(31,12,6.6,6.3),
      `t½ (h)`=c(31,4.5,3.2,5.5),
      Target=c("Tubulin/NLRP3","COX-1/2","GR (NF-κB)","IL-1R (Ka~0.5µg/mL)")
    )
  }, striped=TRUE, bordered=TRUE)

  # ---- TAB 3: Inflammation ----
  output$il1b_plot  <- renderPlotly({ df<-sim_result(); pl(df,"IL1b", "IL-1β (pg/mL)", "#c0392b") })
  output$nlrp3_plot <- renderPlotly({ df<-sim_result(); pl(df,"NLRP3_act","NLRP3 (0-1)","#e67e22") })
  output$il6_plot   <- renderPlotly({ df<-sim_result(); pl(df,"IL6","IL-6 (pg/mL)","#2980b9") })
  output$pge2_plot  <- renderPlotly({ df<-sim_result(); pl(df,"PGE2","PGE2 (pg/mL)","#f39c12") })
  output$neut_plot  <- renderPlotly({ df<-sim_result(); pl(df,"Neutrophil","SF Neutrophils (cells/µL)","#27ae60") })
  output$lipox_plot <- renderPlotly({ df<-sim_result(); pl(df,"LipoxA4","Lipoxin A4 (pg/mL)","#1abc9c") })

  # ---- TAB 4: Crystal ----
  output$ppi_plot      <- renderPlotly({ df<-sim_result(); pl(df,"PPi_ext","PPi (µM)","#8e44ad") })
  output$cart_cryst_plot<-renderPlotly({ df<-sim_result(); pl(df,"Cryst_cart","CPPD in Cartilage (mg/g)","#d35400") })
  output$sf_cryst_plot <- renderPlotly({ df<-sim_result(); pl(df,"Cryst_SF","CPPD in SF (µg/mL)","#c0392b") })

  # ---- TAB 5: Clinical ----
  output$pain_plot <- renderPlotly({ df<-sim_result(); pl(df,"PainVAS","Pain VAS (0-10)","#e74c3c") })
  output$crp_plot  <- renderPlotly({ df<-sim_result(); pl(df,"CRP","CRP (mg/L)","#e67e22") })
  output$ci_plot   <- renderPlotly({ df<-sim_result(); pl(df,"CartInteg","Cartilage Integrity (%)","#2980b9") })

  output$ktp_tbl <- renderDT({
    df <- sim_result()
    ktp <- df %>% filter(near(time_d, 1) | near(time_d, 3) | near(time_d, 7) |
                           near(time_d, round(input$sim_days,0))) %>%
      group_by(Day=round(time_d,0)) %>% slice_tail(n=1) %>% ungroup() %>%
      select(Day, IL1b, Neutrophil, IL6, CRP, PainVAS, CartInteg) %>%
      rename(`IL-1β (pg/mL)`=IL1b, `SF Neut (cells/µL)`=Neutrophil,
             `IL-6 (pg/mL)`=IL6, `CRP (mg/L)`=CRP,
             `Pain VAS`=PainVAS, `CartInteg (%)`=CartInteg) %>%
      mutate(across(where(is.numeric), ~round(.x,1)))
    datatable(ktp, options=list(pageLength=8, dom="t"), rownames=FALSE)
  })

  # ---- TAB 6: Scenarios ----
  output$scen_il1b <- renderPlotly({
    df <- scen_result()
    plot_ly(df, x=~time_d, y=~IL1b, color=~scenario, colors=scen_pal,
            type="scatter", mode="lines") %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="IL-1β (pg/mL)"),
             title="IL-1β", showlegend=TRUE)
  })
  output$scen_pain <- renderPlotly({
    df <- scen_result()
    plot_ly(df, x=~time_d, y=~PainVAS, color=~scenario, colors=scen_pal,
            type="scatter", mode="lines") %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="Pain VAS", range=c(0,10)),
             title="Pain VAS")
  })
  output$scen_crp <- renderPlotly({
    df <- scen_result()
    plot_ly(df, x=~time_d, y=~CRP, color=~scenario, colors=scen_pal,
            type="scatter", mode="lines") %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="CRP (mg/L)"), title="CRP")
  })
  output$scen_nlrp3 <- renderPlotly({
    df <- scen_result()
    plot_ly(df, x=~time_d, y=~NLRP3_act, color=~scenario, colors=scen_pal,
            type="scatter", mode="lines") %>%
      layout(xaxis=list(title="Days"), yaxis=list(title="NLRP3 (0-1)"), title="NLRP3 Activation")
  })

  output$scen_tbl <- renderDT({
    df <- scen_result() %>% filter(near(time_d,7)) %>%
      group_by(scenario) %>% slice_tail(n=1) %>% ungroup() %>%
      select(scenario, IL1b, Neutrophil, CRP, PainVAS, CartInteg) %>%
      rename(Scenario=scenario, `IL-1β`=IL1b, `SF Neut`=Neutrophil,
             `CRP`=CRP, Pain=PainVAS, `CartInteg%`=CartInteg) %>%
      mutate(across(where(is.numeric), ~round(.x,1)))
    datatable(df, options=list(dom="t", pageLength=8), rownames=FALSE)
  })

  # ---- TAB 7: Biomarkers & Sensitivity ----
  output$biomarker_tbl <- renderTable({
    df <- sim_result()
    row7 <- df %>% filter(near(time_d, min(7, max(time_d)))) %>% slice_tail(n=1)
    data.frame(
      Biomarker = c("IL-1β (pg/mL)","IL-6 (pg/mL)","CRP (mg/L)","SF Neutrophils (×10³/µL)","PGE2 (pg/mL)","Lipoxin A4 (pg/mL)","NLRP3 (0-1)","Pain VAS (0-10)","Cartilage Integrity (%)","CPPD in SF (µg/mL)","PPi (µM)"),
      `Day 7 Value` = round(c(row7$IL1b, row7$IL6, row7$CRP, row7$Neutrophil/1000, row7$PGE2, row7$LipoxA4, row7$NLRP3_act, row7$PainVAS, row7$CartInteg, row7$Cryst_SF, row7$PPi_ext),2),
      `Normal Range` = c("<10","<5","<5","<0.5","<100","10-50","0","<2","95-100","0","3-5")
    )
  }, striped=TRUE, bordered=TRUE)

  output$biomarker_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df) %>%
      add_lines(x=~time_d, y=~IL1b,   name="IL-1β",    line=list(color="#e74c3c")) %>%
      add_lines(x=~time_d, y=~IL6,    name="IL-6",     line=list(color="#3498db")) %>%
      add_lines(x=~time_d, y=~CRP,    name="CRP×10",   yaxis="y2", line=list(color="#f39c12",dash="dash")) %>%
      layout(
        xaxis=list(title="Time (days)"),
        yaxis =list(title="Cytokines (pg/mL)", side="left"),
        yaxis2=list(title="CRP (mg/L)×10", overlaying="y", side="right"),
        legend=list(orientation="h")
      )
  })

  output$sens_plot <- renderPlotly({
    base_params <- list(
      kIL1b_syn=2.0, kNLRP3_on=0.5, EC50_cryst_nlrp3=5.0,
      kNeut_rec=8.0, kCart_dam=2e-4, kIL1b_deg=0.4,
      kTNAP=0.20, kENPP1=0.15, kShed_base=3e-3, IL1Ra_endo=0.6
    )
    pnames <- names(base_params)
    tend <- 7 * 24
    tg   <- seq(0, tend, by=2)
    shed <- ev(time=0, cmt="Cryst_SF", amt=15)

    results <- map_dfr(pnames, function(p) {
      base_val <- base_params[[p]]
      run_p <- function(mult) {
        pv <- base_params
        pv[[p]] <- base_val * mult
        out <- mod %>% init(INIT_DEF) %>%
          param(pv) %>%
          mrgsim_e(shed, tgrid=tg, output="df") %>%
          filter(near(time, tend)) %>% slice_tail(n=1)
        data.frame(param=p, multiplier=mult, IL1b=out$IL1b, Pain=out$PainVAS)
      }
      bind_rows(run_p(0.5), run_p(1.0), run_p(2.0))
    })

    results %>% filter(multiplier != 1) %>%
      mutate(label=paste0(param," ×",multiplier),
             IL1b_pct = IL1b / filter(results, multiplier==1)$IL1b[match(param, filter(results,multiplier==1)$param)] - 1) -> res2

    # Bar chart of % deviation from baseline
    bl <- filter(results, multiplier==1) %>% select(param, IL1b_base=IL1b, Pain_base=Pain)
    res3 <- left_join(filter(results, multiplier!=1), bl, by="param") %>%
      mutate(IL1b_dev = (IL1b - IL1b_base)/IL1b_base * 100,
             mult_label = ifelse(multiplier==0.5,"×0.5 (−50%)","×2.0 (+100%)")) %>%
      arrange(abs(IL1b_dev))

    plot_ly(res3, y=~param, x=~IL1b_dev, color=~mult_label,
            colors=c("×0.5 (−50%)"="#2980b9","×2.0 (+100%)"="#e74c3c"),
            type="bar", orientation="h") %>%
      layout(barmode="group",
             xaxis=list(title="% Change in Day-7 IL-1β"),
             yaxis=list(title="Parameter"),
             title="Sensitivity: parameter ×0.5 vs ×2.0 → IL-1β change at Day 7")
  })
}

shinyApp(ui, server)
