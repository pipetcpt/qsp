##############################################################################
# Non-Hodgkin Lymphoma (DLBCL) — Shiny QSP Dashboard
# 6-Tab Interactive Application
#
# Tabs:
#   1. Patient Profile & Disease Classification
#   2. Drug PK — Rituximab & CHOP Components
#   3. Tumor Dynamics & Response
#   4. Clinical Endpoints (ORR, PFS, Waterfall)
#   5. Scenario Comparison (6 treatment arms)
#   6. Biomarkers & Toxicity (ANC, CRS, Resistance)
##############################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(scales)

# ── Inline mrgsolve model (core model from nhl_mrgsolve_model.R) ─────────────
MODEL_CODE <- '
$PARAM
CLR=0.008, VcR=3.1, QR=0.12, VpR=3.7,
konR=0.27, koffR=0.0003, kintR=0.05,
ksynCD20=0.012, kdegCD20=0.004, Bmax_CD20=0.35,
kelCPP=0.18, VdCPP=50,
CLDox=45, VcDox=25, QDox=15, VpDox=400,
kaVEN=0.2, FVEN=0.5, kelVEN=0.035, VdVEN=256,
kaIBR=0.5, FIBR=0.03, kelIBR=0.22, VdIBR=10000,
kg_tumor=0.012, kd_base=0.003, K_carry=1000, T0=100,
EC50_RTX=0.15, Emax_RTX=0.85,
EC50_CPP=0.8, Emax_CPP=0.80,
EC50_Dox=0.05, Emax_Dox=0.75,
EC50_VEN=0.3, Emax_VEN=0.80,
EC50_IBR=0.12, Emax_IBR=0.50, Hill_n=1.5,
kBCR_base=0.02, kBCR_deg=0.15, BCR_max=1.0, IBR_BCR_IC50=0.08,
BCL2_Kd=0.2, BCL2_0=1.0,
kNK_in=0.5, kNK_out=0.1, NK0=5, RTX_NK_stim=0.3,
kCD8_in=0.3, kCD8_out=0.05, CD8_0=10, kexhaust=0.01,
ANC_0=5, kANC_rec=0.08, kANC_CPP_kill=0.25, kANC_Dox_kill=0.15, ANC_min=0.1,
kresist=0.0003, resist_0=0.01,
is_ABC=0, is_DHL=0, DHL_growth_mult=1.5

$CMT RuxCent RuxPeriph CD20_free CD20_RTX CPP_active DoxCent DoxPeriph
     VEN_gut VEN_cent IBR_cent
     Tumor BCR_signal BCL2_occ NK_cells CD8_cells ANC Resistance
     Cum_RTX_dose Tumor_resp CRS_risk

$MAIN
double kg_eff = kg_tumor * (1 + is_DHL * (DHL_growth_mult - 1));
double CRux  = RuxCent / VcR;
double CCPP  = CPP_active / VdCPP;
double CDox  = DoxCent / VcDox;
double CVEN  = VEN_cent / VdVEN;
double CIBR  = IBR_cent / VdIBR;

$ODE
double koffR_eff = koffR + kintR;
dxdt_RuxCent   = -CLR/VcR*RuxCent - QR/VcR*RuxCent + QR/VpR*RuxPeriph
                  - konR*(RuxCent/VcR)*CD20_free + koffR*CD20_RTX;
dxdt_RuxPeriph = QR/VcR*RuxCent - QR/VpR*RuxPeriph;
double CD20_prod = ksynCD20*(Tumor/100.0);
dxdt_CD20_free = CD20_prod - kdegCD20*CD20_free
                  - konR*(RuxCent/VcR)*CD20_free + koffR*CD20_RTX;
dxdt_CD20_RTX  = konR*(RuxCent/VcR)*CD20_free - koffR_eff*CD20_RTX;
dxdt_CPP_active= -kelCPP*CPP_active;
dxdt_DoxCent   = -CLDox/VcDox*DoxCent - QDox/VcDox*DoxCent + QDox/VpDox*DoxPeriph;
dxdt_DoxPeriph = QDox/VcDox*DoxCent - QDox/VpDox*DoxPeriph;
dxdt_VEN_gut   = -kaVEN*VEN_gut;
dxdt_VEN_cent  = FVEN*kaVEN*VEN_gut - kelVEN*VEN_cent;
dxdt_IBR_cent  = FIBR*(-kelIBR*IBR_cent);
double IBR_BCR_inh = (is_ABC>0.5) ? (CIBR/(CIBR+IBR_BCR_IC50))*0.8 : 0.0;
dxdt_BCR_signal= kBCR_base*(BCR_max-BCR_signal)*(1-IBR_BCR_inh) - kBCR_deg*BCR_signal;
double VEN_occ = CVEN/(CVEN+BCL2_Kd);
dxdt_BCL2_occ  = kaVEN*(VEN_occ-BCL2_occ);
double NK_RTX_boost = 1.0 + RTX_NK_stim*(CD20_RTX/(CD20_RTX+0.05));
dxdt_NK_cells  = kNK_in*NK_RTX_boost - kNK_out*NK_cells;
double PD1_exhaust = kexhaust*(Tumor/100.0);
dxdt_CD8_cells = kCD8_in - (kCD8_out+PD1_exhaust)*CD8_cells;
double E_RTX = Emax_RTX*pow(CRux,Hill_n)/(pow(CRux,Hill_n)+pow(EC50_RTX,Hill_n));
double E_ADCC = E_RTX*(NK_cells/NK0)*1.2; if(E_ADCC>Emax_RTX) E_ADCC=Emax_RTX;
double E_CPP  = Emax_CPP*pow(CCPP,Hill_n)/(pow(CCPP,Hill_n)+pow(EC50_CPP,Hill_n));
double E_Dox  = Emax_Dox*pow(CDox,Hill_n)/(pow(CDox,Hill_n)+pow(EC50_Dox,Hill_n));
double E_VEN  = Emax_VEN*BCL2_occ;
double E_IBR  = (is_ABC>0.5) ? Emax_IBR*CIBR/(CIBR+EC50_IBR)*BCR_signal : 0.0;
double E_immune = 0.01*(NK_cells/NK0+CD8_cells/CD8_0)*(Tumor/100.0);
double resist_factor = 1-Resistance*0.8;
double E_total = (1-(1-E_ADCC)*(1-E_CPP)*(1-E_Dox)*(1-E_VEN)*(1-E_IBR))*resist_factor;
if(E_total>0.98) E_total=0.98;
double growth_term = kg_eff*Tumor*(1-Tumor/K_carry);
double kill_term = (kd_base+E_total)*Tumor + E_immune;
dxdt_Tumor = growth_term - kill_term;
if(Tumor<0.01) dxdt_Tumor=0;
dxdt_Resistance = kresist*Tumor*(1-Resistance);
double ANC_drive = kANC_rec*(ANC_0-ANC);
double ANC_CPP_kill = kANC_CPP_kill*CCPP*ANC/ANC_0;
double ANC_Dox_kill = kANC_Dox_kill*CDox*ANC/ANC_0;
dxdt_ANC = ANC_drive - ANC_CPP_kill - ANC_Dox_kill;
if(ANC<ANC_min) dxdt_ANC=0;
dxdt_Cum_RTX_dose = RuxCent/VcR;
dxdt_Tumor_resp = 0;
dxdt_CRS_risk = 0.1*kill_term*(1-CRS_risk) - 0.05*CRS_risk;

$TABLE
capture CRux_mgL = RuxCent/VcR;
capture CCPP_mgL = CPP_active/VdCPP;
capture CDox_mgL = DoxCent/VcDox;
capture CVEN_mgL = VEN_cent/VdVEN;
capture CIBR_mgL = IBR_cent/VdIBR;
capture SPD_pct  = Tumor/100.0*100.0;
capture CR_flag  = (Tumor/100.0 < 0.05) ? 1.0 : 0.0;
capture PR_flag  = (Tumor/100.0>=0.05 && Tumor/100.0<0.50) ? 1.0 : 0.0;
capture PD_flag  = (Tumor/100.0 > 1.50) ? 1.0 : 0.0;
capture ANC_out  = ANC;

$INIT
RuxCent=0, RuxPeriph=0, CD20_free=0.30, CD20_RTX=0, CPP_active=0,
DoxCent=0, DoxPeriph=0, VEN_gut=0, VEN_cent=0, IBR_cent=0,
Tumor=100, BCR_signal=0.133, BCL2_occ=0, NK_cells=5, CD8_cells=10,
ANC=5.0, Resistance=0.01, Cum_RTX_dose=0, Tumor_resp=100, CRS_risk=0
'

# Compile once
mod_global <- mcode("dlbcl_shiny", MODEL_CODE, quiet = TRUE)

# ── Helper functions ─────────────────────────────────────────────────────────
make_events_rchop <- function(n_cycles, rtx_mg, cpp_mg, dox_mg) {
  evs <- lapply(seq_len(n_cycles), function(i) {
    day_h <- (i - 1) * 21 * 24
    c(ev(cmt=1, time=day_h, amt=rtx_mg, rate=-2),
      ev(cmt=5, time=day_h, amt=cpp_mg, rate=-2),
      ev(cmt=6, time=day_h, amt=dox_mg, rate=-1))
  })
  do.call(c, evs)
}

make_events_ven <- function(dose_mg, days) {
  ev(cmt=8, time=0, amt=dose_mg, ii=24, addl=days-1)
}

make_events_ibr <- function(dose_mg, days) {
  ev(cmt=10, time=0, amt=dose_mg, ii=24, addl=days-1)
}

run_sim <- function(mod, params_list, events_obj, end_days=365, delta=6) {
  p <- do.call(param, c(list(mod), params_list))
  out <- mrgsim(p, events=events_obj, end=end_days*24, delta=delta)
  df  <- as.data.frame(out)
  df$time_day <- df$time / 24
  df
}

theme_dash <- theme_bw(base_size=13) +
  theme(legend.position="bottom",
        panel.grid.minor=element_blank(),
        strip.background=element_rect(fill="#E3F2FD"))

# ── UI ───────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "DLBCL QSP Dashboard", titleWidth = 280),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",     tabName="tab_patient",   icon=icon("user-md")),
      menuItem("Drug PK",             tabName="tab_pk",        icon=icon("pills")),
      menuItem("Tumor Dynamics",      tabName="tab_tumor",     icon=icon("chart-line")),
      menuItem("Clinical Endpoints",  tabName="tab_endpoints", icon=icon("stethoscope")),
      menuItem("Scenario Comparison", tabName="tab_scenarios", icon=icon("vials")),
      menuItem("Biomarkers & Toxicity",tabName="tab_bio",      icon=icon("flask"))
    ),
    hr(),
    h5("  Global Settings", style="color:#BBB;padding-left:12px"),
    sliderInput("sim_days",  "Simulation horizon (days)", 180, 730, 365, step=30),
    selectInput("subtype",   "DLBCL Subtype",
                choices = c("GCB-DLBCL"="gcb","ABC-DLBCL"="abc","Double-Hit (DHL)"="dhl"),
                selected = "gcb"),
    numericInput("bsa_m2", "BSA (m²)", value=1.73, min=1.2, max=2.5, step=0.05),
    hr(),
    h5("  R-CHOP Dosing", style="color:#BBB;padding-left:12px"),
    sliderInput("n_cycles",  "Number of R-CHOP cycles",  4, 8, 6, step=1),
    sliderInput("rtx_per_m2","Rituximab (mg/m²)",       250, 500, 375, step=25),
    sliderInput("cpp_per_m2","Cyclophosphamide (mg/m²)",500, 1000, 750, step=50),
    sliderInput("dox_per_m2","Doxorubicin (mg/m²)",     35, 75, 50, step=5),
    hr(),
    h5("  Add-On Agents", style="color:#BBB;padding-left:12px"),
    checkboxInput("add_venetoclax", "Add Venetoclax (BCL-2 inh.)", FALSE),
    conditionalPanel(
      condition = "input.add_venetoclax == true",
      sliderInput("ven_dose", "Venetoclax mg/day", 200, 1200, 800, step=100)
    ),
    checkboxInput("add_ibrutinib", "Add Ibrutinib (BTK inh.)", FALSE),
    conditionalPanel(
      condition = "input.add_ibrutinib == true",
      sliderInput("ibr_dose", "Ibrutinib mg/day", 280, 840, 560, step=140)
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML(
      ".content-wrapper,.right-side{background-color:#F5F7FA;}
       .box-header{background:#1565C0!important;color:white!important;}
       .value-box .icon{font-size:40px!important;}"
    ))),
    tabItems(

      # ── TAB 1: Patient Profile ──────────────────────────────────────────────
      tabItem(tabName="tab_patient",
        fluidRow(
          valueBoxOutput("vbox_subtype",  width=3),
          valueBoxOutput("vbox_ipi",      width=3),
          valueBoxOutput("vbox_bcl2_status", width=3),
          valueBoxOutput("vbox_prognosis",   width=3)
        ),
        fluidRow(
          box(title="Disease Characterization & IPI Score", width=6, solidHeader=TRUE,
            DTOutput("dt_patient_profile")
          ),
          box(title="DLBCL Subtype Biology", width=6, solidHeader=TRUE,
            p(strong("GCB-DLBCL (Germinal Center B-cell type):")),
            p("BCL6+, IRF4/MUM1−, CD10+. Associated with t(14;18)/BCL-2 rearrangement,
               EZH2 Y641 mutations. Better prognosis (~60% 5yr OS). R-CHOP is standard."),
            p(strong("ABC-DLBCL (Activated B-cell type):")),
            p("MUM1/IRF4+, FOXP1+. Characterized by constitutive NF-κB activation,
               MYD88 L265P mutation (~30%), CD79A/B ITAM mutations. Worse prognosis
               (~40% 5yr OS). Ibrutinib may add benefit in NF-κB driven cases."),
            p(strong("Double-Hit Lymphoma (DHL):")),
            p("MYC rearrangement + BCL-2 or BCL-6 rearrangement. High-grade B-cell
               lymphoma. Very aggressive (median OS ~12-18 months with R-CHOP alone).
               EPOCH-R or DA-EPOCH-R often used instead."),
            hr(),
            p(strong("Key Genomic Drivers by Subtype:")),
            DTOutput("dt_genomics")
          )
        ),
        fluidRow(
          box(title="Mechanistic Map Overview", width=12, solidHeader=TRUE,
            tags$img(src=NULL,
              style="max-width:100%;border:1px solid #ddd;border-radius:4px;",
              alt="DLBCL mechanistic map — see nhl_qsp_model.svg"),
            p(em("Mechanistic map: nhl_qsp_model.svg (14 clusters, 120+ nodes)"),
              style="color:#666;font-size:11px;margin-top:5px;")
          )
        )
      ),

      # ── TAB 2: Drug PK ──────────────────────────────────────────────────────
      tabItem(tabName="tab_pk",
        fluidRow(
          box(title="Rituximab — 2-Compartment PK with TMDD", width=6,
              solidHeader=TRUE,
            plotOutput("plot_rtx_pk", height=320)
          ),
          box(title="CHOP Component Concentrations", width=6, solidHeader=TRUE,
            plotOutput("plot_chop_pk", height=320)
          )
        ),
        fluidRow(
          box(title="CD20 Receptor Occupancy", width=6, solidHeader=TRUE,
            plotOutput("plot_cd20_occ", height=280)
          ),
          box(title="Add-On Agent PK (Venetoclax / Ibrutinib)", width=6,
              solidHeader=TRUE,
            plotOutput("plot_addon_pk", height=280)
          )
        ),
        fluidRow(
          box(title="PK Parameter Summary", width=12, solidHeader=TRUE,
            DTOutput("dt_pk_params")
          )
        )
      ),

      # ── TAB 3: Tumor Dynamics ───────────────────────────────────────────────
      tabItem(tabName="tab_tumor",
        fluidRow(
          valueBoxOutput("vbox_spd_d126",  width=3),
          valueBoxOutput("vbox_cr_d126",   width=3),
          valueBoxOutput("vbox_nadir_day", width=3),
          valueBoxOutput("vbox_resistance",width=3)
        ),
        fluidRow(
          box(title="Tumor Burden (SPD % of Baseline)", width=8, solidHeader=TRUE,
            plotOutput("plot_tumor_dyn", height=380)
          ),
          box(title="Response Classification", width=4, solidHeader=TRUE,
            plotOutput("plot_resp_pie", height=380)
          )
        ),
        fluidRow(
          box(title="NK Cell & CD8+ T-Cell Dynamics", width=6, solidHeader=TRUE,
            plotOutput("plot_immune_eff", height=280)
          ),
          box(title="BCR/NF-κB Signal Suppression", width=6, solidHeader=TRUE,
            plotOutput("plot_bcr_signal", height=280)
          )
        )
      ),

      # ── TAB 4: Clinical Endpoints ───────────────────────────────────────────
      tabItem(tabName="tab_endpoints",
        fluidRow(
          box(title="Response Rate Over Time", width=8, solidHeader=TRUE,
            plotOutput("plot_response_kinetics", height=320)
          ),
          box(title="Clinical Endpoint Summary", width=4, solidHeader=TRUE,
            DTOutput("dt_endpoints_sum")
          )
        ),
        fluidRow(
          box(title="Population Waterfall (N=50 Virtual Patients)", width=7,
              solidHeader=TRUE,
            plotOutput("plot_waterfall", height=320)
          ),
          box(title="Lugano Response Criteria", width=5, solidHeader=TRUE,
            DTOutput("dt_lugano")
          )
        ),
        fluidRow(
          box(title="PFS Simulation (Kaplan-Meier Approximation)", width=12,
              solidHeader=TRUE,
            plotOutput("plot_pfs_km", height=280)
          )
        )
      ),

      # ── TAB 5: Scenario Comparison ──────────────────────────────────────────
      tabItem(tabName="tab_scenarios",
        fluidRow(
          box(title="Treatment Scenarios: Tumor Burden Comparison", width=8,
              solidHeader=TRUE,
            plotOutput("plot_scenarios_tumor", height=380)
          ),
          box(title="Scenario Parameters", width=4, solidHeader=TRUE,
            DTOutput("dt_scenarios_table")
          )
        ),
        fluidRow(
          box(title="End-of-Treatment Response Comparison (Day 126)",
              width=6, solidHeader=TRUE,
            plotOutput("plot_eot_response", height=280)
          ),
          box(title="ANC Myelosuppression by Scenario", width=6, solidHeader=TRUE,
            plotOutput("plot_anc_scenarios", height=280)
          )
        ),
        fluidRow(
          box(title="Resistance Development Over Time", width=12, solidHeader=TRUE,
            plotOutput("plot_resistance", height=250)
          )
        )
      ),

      # ── TAB 6: Biomarkers & Toxicity ────────────────────────────────────────
      tabItem(tabName="tab_bio",
        fluidRow(
          valueBoxOutput("vbox_anc_nadir",  width=3),
          valueBoxOutput("vbox_anc_grade",  width=3),
          valueBoxOutput("vbox_crs_risk",   width=3),
          valueBoxOutput("vbox_bcl2_occ",   width=3)
        ),
        fluidRow(
          box(title="ANC Time Profile (Myelosuppression)", width=6,
              solidHeader=TRUE,
            plotOutput("plot_anc_biomarker", height=300)
          ),
          box(title="BCL-2 Occupancy (Venetoclax PD)", width=6, solidHeader=TRUE,
            plotOutput("plot_bcl2_occ", height=300)
          )
        ),
        fluidRow(
          box(title="CRS Risk Index Over Time", width=6, solidHeader=TRUE,
            plotOutput("plot_crs", height=280)
          ),
          box(title="Key Biomarker & Toxicity Table", width=6, solidHeader=TRUE,
            DTOutput("dt_biomarkers")
          )
        ),
        fluidRow(
          box(title="Cardiotoxicity Monitor (Cumulative Doxorubicin)", width=12,
              solidHeader=TRUE,
            plotOutput("plot_dox_cumulative", height=220)
          )
        )
      )

    ) # end tabItems
  ) # end dashboardBody
)

# ── SERVER ───────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Reactive: build parameter list
  params_reactive <- reactive({
    is_abc <- as.numeric(input$subtype == "abc")
    is_dhl <- as.numeric(input$subtype == "dhl")
    list(is_ABC = is_abc, is_DHL = is_dhl)
  })

  # Reactive: build dosing events
  events_reactive <- reactive({
    bsa   <- input$bsa_m2
    rtx_m  <- input$rtx_per_m2 * bsa
    cpp_m  <- input$cpp_per_m2 * bsa
    dox_m  <- input$dox_per_m2 * bsa
    evs   <- make_events_rchop(input$n_cycles, rtx_m, cpp_m, dox_m)
    if (input$add_venetoclax) {
      evs <- c(evs, make_events_ven(input$ven_dose, input$sim_days))
    }
    if (input$add_ibrutinib) {
      evs <- c(evs, make_events_ibr(input$ibr_dose, 21 * input$n_cycles))
    }
    evs
  })

  # Primary simulation
  sim_data <- reactive({
    run_sim(mod_global, params_reactive(), events_reactive(),
            end_days = input$sim_days, delta = 6)
  })

  # All 6 scenarios
  all_scenarios <- reactive({
    bsa  <- input$bsa_m2
    scn_list <- list(
      list(lbl="1. Untreated",         params=list(is_ABC=0,is_DHL=0),
           evs=ev(cmt=1,time=0,amt=0)),
      list(lbl="2. R-CHOP×6",          params=list(is_ABC=0,is_DHL=0),
           evs=make_events_rchop(6,375*bsa,750*bsa,50*bsa)),
      list(lbl="3. Pola-R-CHP×6",      params=list(is_ABC=0,is_DHL=0,Emax_CPP=0.85),
           evs=make_events_rchop(6,375*bsa,750*bsa,50*bsa)),
      list(lbl="4. R-CHOP+Ibrutinib",  params=list(is_ABC=1,is_DHL=0),
           evs=c(make_events_rchop(6,375*bsa,750*bsa,50*bsa),
                  make_events_ibr(560,126))),
      list(lbl="5. R-CHOP+Venetoclax", params=list(is_ABC=0,is_DHL=0),
           evs=c(make_events_rchop(6,375*bsa,750*bsa,50*bsa),
                  make_events_ven(800,126))),
      list(lbl="6. R-CHOP DHL",        params=list(is_ABC=0,is_DHL=1),
           evs=make_events_rchop(6,375*bsa,750*bsa,50*bsa))
    )
    bind_rows(lapply(scn_list, function(s) {
      df <- run_sim(mod_global, s$params, s$evs, input$sim_days, 12)
      df$Scenario <- s$lbl
      df
    }))
  })

  # Population sim (N=50)
  pop_sim <- reactive({
    bsa <- input$bsa_m2
    set.seed(123)
    idata <- data.frame(
      ID = 1:50,
      kg_tumor = rlnorm(50, log(0.012), 0.4),
      T0       = rlnorm(50, log(100), 0.5),
      EC50_RTX = rlnorm(50, log(0.15), 0.3)
    )
    out <- mrgsim(mod_global,
                  idata  = idata,
                  events = make_events_rchop(6,375*bsa,750*bsa,50*bsa),
                  end    = input$sim_days*24, delta=24)
    df  <- as.data.frame(out)
    df$time_day <- df$time/24
    df
  })

  # ── Tab 1: Patient Profile ────────────────────────────────────────────────
  output$vbox_subtype <- renderValueBox({
    subtype_lab <- switch(input$subtype,
      "gcb"="GCB-DLBCL","abc"="ABC-DLBCL","dhl"="Double-Hit (DHL)")
    valueBox(subtype_lab, "DLBCL Subtype", icon=icon("dna"), color="blue")
  })
  output$vbox_ipi <- renderValueBox({
    valueBox("IPI 3 (High-Intermediate)", "IPI Risk", icon=icon("chart-bar"), color="orange")
  })
  output$vbox_bcl2_status <- renderValueBox({
    bcl2_flag <- if (input$add_venetoclax) "BCL-2 High (VEN)" else "Standard"
    valueBox(bcl2_flag, "BCL-2 Status", icon=icon("circle"), color="purple")
  })
  output$vbox_prognosis <- renderValueBox({
    prog <- switch(input$subtype,
      "gcb"="Good (5yr OS ~60%)", "abc"="Moderate (5yr OS ~40%)",
      "dhl"="Poor (5yr OS ~25%)")
    valueBox(prog, "Prognosis (R-CHOP)", icon=icon("heartbeat"), color="red")
  })
  output$dt_patient_profile <- renderDT({
    data.frame(
      Parameter = c("Age","Performance Status","Stage","Extranodal sites",
                    "LDH","β2-Microglobulin","ECOG PS","Ki-67","IPI Score",
                    "COO (Hans algorithm)"),
      Value     = c("58 yrs","ECOG 1","Stage III-IV","1","Elevated (2× ULN)",
                    "3.5 mg/L (elevated)","1","85%","3 (High-Intermediate)",
                    switch(input$subtype,"gcb"="GCB","abc"="Non-GCB (ABC)",
                           "dhl"="Non-GCB (DHL)"))
    )
  }, options=list(pageLength=10,dom='t'), rownames=FALSE)

  output$dt_genomics <- renderDT({
    data.frame(
      Subtype   = c("GCB","GCB","ABC","ABC","DHL","DHL"),
      Alteration = c("t(14;18) BCL-2 (30-45%)","EZH2 Y641 (22%)",
                     "MYD88 L265P (30%)","CD79A/B ITAM (20%)",
                     "MYC rearrangement","BCL-2 rearrangement"),
      Frequency = c("30-45%","22%","30%","20%","5-10% of DLBCL","~5%"),
      Drug_Target = c("Venetoclax","EZH2 inhibitors","Ibrutinib",
                      "Ibrutinib","DA-EPOCH-R","Venetoclax")
    )
  }, options=list(pageLength=6,dom='t'), rownames=FALSE)

  # ── Tab 2: PK ─────────────────────────────────────────────────────────────
  output$plot_rtx_pk <- renderPlot({
    df <- sim_data()
    ggplot(df %>% filter(CRux_mgL > 0.001), aes(x=time_day, y=CRux_mgL)) +
      geom_line(color="#1565C0", linewidth=1.5) +
      scale_y_log10(labels=comma) +
      labs(title="Rituximab Central Compartment",
           x="Time (days)", y="Conc. (mg/L, log scale)") +
      theme_dash
  })
  output$plot_chop_pk <- renderPlot({
    df <- sim_data() %>%
      select(time_day, CCPP_mgL, CDox_mgL) %>%
      pivot_longer(cols=-time_day, names_to="drug", values_to="conc") %>%
      filter(conc > 1e-6) %>%
      mutate(drug = recode(drug, "CCPP_mgL"="4-OH-Cyclophos", "CDox_mgL"="Doxorubicin"))
    ggplot(df, aes(x=time_day, y=conc, color=drug)) +
      geom_line(linewidth=1.2) +
      scale_y_log10(labels=comma) +
      scale_color_brewer(palette="Set1") +
      labs(title="CHOP Active Concentrations", x="Time (days)",
           y="Conc. (mg/L, log)", color="Agent") +
      theme_dash
  })
  output$plot_cd20_occ <- renderPlot({
    df <- sim_data()
    df$CD20_total <- df$CD20_free + df$CD20_RTX
    df$Occ <- with(df, ifelse(CD20_total > 0, CD20_RTX/CD20_total, 0))
    ggplot(df, aes(x=time_day, y=Occ*100)) +
      geom_line(color="#388E3C", linewidth=1.5) +
      geom_hline(yintercept=80, linetype="dashed", color="orange") +
      labs(title="CD20 Receptor Occupancy by Rituximab",
           x="Time (days)", y="Occupancy (%)") +
      ylim(0,100) + theme_dash
  })
  output$plot_addon_pk <- renderPlot({
    df <- sim_data()
    plot_df <- data.frame(time_day=df$time_day,
                          Venetoclax=df$CVEN_mgL, Ibrutinib=df$CIBR_mgL) %>%
      pivot_longer(-time_day, names_to="drug", values_to="conc") %>%
      filter(conc > 1e-6)
    if (nrow(plot_df) == 0) {
      ggplot() + labs(title="No add-on agents selected", x="", y="") + theme_dash
    } else {
      ggplot(plot_df, aes(x=time_day, y=conc, color=drug)) +
        geom_line(linewidth=1.2) +
        scale_y_log10(labels=comma) +
        scale_color_manual(values=c("Venetoclax"="#7B1FA2","Ibrutinib"="#E65100")) +
        labs(title="Add-On Agent Concentrations",
             x="Time (days)", y="Conc. (mg/L, log)") +
        theme_dash
    }
  })
  output$dt_pk_params <- renderDT({
    data.frame(
      Drug     = c("Rituximab","Rituximab","Rituximab","4-OH-Cyclophosphamide",
                   "Doxorubicin","Doxorubicin","Venetoclax","Ibrutinib"),
      Parameter= c("CL (linear)","Vc","T½ (terminal ~21d)","T½ (active met.)",
                   "CL","Vc","T½","T½"),
      Value    = c("0.008 L/h","3.1 L","~21 days (FcRn)","~3.8h (4-OH-CPP)",
                   "45 L/h","25 L","~19h","~3h"),
      Source   = c("Berinstein 1998","Maloney 1999","Avivi 2014","Boddy 2003",
                   "Dobbs 2006","Dobbs 2006","Munakata 2022","de Claro 2013")
    )
  }, options=list(dom='t', pageLength=10), rownames=FALSE)

  # ── Tab 3: Tumor Dynamics ─────────────────────────────────────────────────
  output$vbox_spd_d126 <- renderValueBox({
    df <- sim_data()
    val <- df %>% filter(abs(time_day-126)<3) %>% slice_tail(n=1) %>% pull(SPD_pct)
    valueBox(paste0(round(val,1),"%"), "SPD at Day 126", icon=icon("ruler"), color="blue")
  })
  output$vbox_cr_d126 <- renderValueBox({
    df <- sim_data()
    cr_flag <- df %>% filter(abs(time_day-126)<3) %>% slice_tail(n=1) %>% pull(CR_flag)
    col <- if(cr_flag > 0.5) "green" else "orange"
    valueBox(if(cr_flag>0.5)"CR Achieved" else "No CR yet",
             "Day-126 Response", icon=icon("check-circle"), color=col)
  })
  output$vbox_nadir_day <- renderValueBox({
    df <- sim_data()
    nadir_day <- df$time_day[which.min(df$Tumor)]
    valueBox(paste0("Day ", round(nadir_day)), "Tumor Nadir", icon=icon("arrow-down"), color="purple")
  })
  output$vbox_resistance <- renderValueBox({
    df <- sim_data()
    resist_final <- tail(df$Resistance, 1)
    valueBox(paste0(round(resist_final*100,1),"%"), "Resistance at End",
             icon=icon("shield-alt"), color="red")
  })
  output$plot_tumor_dyn <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x=time_day, y=SPD_pct)) +
      geom_line(color="#1565C0", linewidth=1.8) +
      geom_hline(yintercept=c(5,50,150), linetype="dashed",
                 color=c("#4CAF50","#FF9800","#F44336"), linewidth=0.8) +
      annotate("text",x=input$sim_days*0.9,y=7,label="CR (<5%)",color="#4CAF50",size=3.5) +
      annotate("text",x=input$sim_days*0.9,y=52,label="PR (<50%)",color="#FF9800",size=3.5) +
      annotate("text",x=input$sim_days*0.9,y=152,label="PD (>150%)",color="#F44336",size=3.5) +
      scale_y_continuous(limits=c(0,NA)) +
      labs(title="Tumor Burden — SPD (% of Baseline)", x="Time (days)", y="SPD (%)") +
      theme_dash
  })
  output$plot_resp_pie <- renderPlot({
    df <- sim_data()
    latest <- tail(df, 1)
    resp_val <- latest$SPD_pct
    cat_lbl <- if(resp_val<5) "CR" else if(resp_val<50) "PR" else if(resp_val>150) "PD" else "SD"
    pie_df <- data.frame(resp=c("CR","PR","SD","PD"),
                         val=c(if(cat_lbl=="CR")1 else 0,
                                if(cat_lbl=="PR")1 else 0,
                                if(cat_lbl=="SD")1 else 0,
                                if(cat_lbl=="PD")1 else 0))
    ggplot(pie_df, aes(x="",y=val,fill=resp)) +
      geom_bar(stat="identity",width=1,color="white") +
      coord_polar("y") +
      scale_fill_manual(values=c("CR"="#4CAF50","PR"="#8BC34A","SD"="#FFC107","PD"="#F44336")) +
      labs(title=paste("Current response:", cat_lbl), fill="Response") +
      theme_void(base_size=13) + theme(legend.position="bottom")
  })
  output$plot_immune_eff <- renderPlot({
    df <- sim_data() %>%
      select(time_day, NK_cells, CD8_cells) %>%
      pivot_longer(-time_day, names_to="cell", values_to="count")
    ggplot(df, aes(x=time_day, y=count, color=cell)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=c("NK_cells"="#FF6F00","CD8_cells"="#1565C0")) +
      labs(title="Immune Effector Cells", x="Time (days)", y="Normalized count",
           color="Cell type") +
      theme_dash
  })
  output$plot_bcr_signal <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x=time_day, y=BCR_signal)) +
      geom_line(color="#6A1B9A", linewidth=1.5) +
      labs(title=paste("BCR/NF-κB Signal —", toupper(input$subtype)),
           x="Time (days)", y="Normalized signal (0-1)") +
      ylim(0,1) + theme_dash
  })

  # ── Tab 4: Clinical Endpoints ─────────────────────────────────────────────
  output$plot_response_kinetics <- renderPlot({
    df <- sim_data()
    resp_df <- df %>%
      mutate(resp = case_when(SPD_pct<5 ~ "CR", SPD_pct<50 ~ "PR",
                               SPD_pct>150 ~ "PD", TRUE ~ "SD")) %>%
      count(time_day, resp) %>%
      complete(time_day, resp, fill=list(n=0))
    ggplot(df, aes(x=time_day, y=100-SPD_pct)) +
      geom_line(color="#1565C0", linewidth=1.5) +
      geom_hline(yintercept=c(95,50), linetype="dashed",
                 color=c("#4CAF50","#FF9800")) +
      labs(title="Tumor Reduction (% from Baseline)", x="Time (days)",
           y="Reduction (%)") +
      theme_dash
  })
  output$dt_endpoints_sum <- renderDT({
    df <- sim_data()
    timepoints <- c(42,84,126,252,365)
    timepoints <- timepoints[timepoints <= input$sim_days]
    tab <- lapply(timepoints, function(tp) {
      row <- df %>% filter(abs(time_day-tp)<3) %>% slice_tail(n=1)
      data.frame(Day=tp,
                 SPD_pct=round(row$SPD_pct,1),
                 Response=if(row$SPD_pct<5)"CR" else if(row$SPD_pct<50)"PR" else if(row$SPD_pct>150)"PD" else "SD",
                 ANC=round(row$ANC_out,2))
    })
    do.call(rbind, tab)
  }, options=list(dom='t'), rownames=FALSE)
  output$plot_waterfall <- renderPlot({
    pop <- pop_sim()
    wf <- pop %>%
      filter(abs(time_day-126) < 1) %>%
      group_by(ID) %>% slice_tail(n=1) %>% ungroup() %>%
      arrange(SPD_pct) %>%
      mutate(rank=row_number(),
             resp=case_when(SPD_pct<5~"CR",SPD_pct<50~"PR",SPD_pct>150~"PD",TRUE~"SD"))
    ggplot(wf, aes(x=rank, y=SPD_pct-100, fill=resp)) +
      geom_bar(stat="identity") +
      scale_fill_manual(values=c("CR"="#4CAF50","PR"="#8BC34A","SD"="#FFC107","PD"="#F44336")) +
      labs(title="Waterfall Plot (N=50 Virtual Patients, Day 126)",
           x="Patient", y="% change SPD from baseline", fill="Response") +
      theme_dash
  })
  output$dt_lugano <- renderDT({
    data.frame(
      Response=c("Complete Response (CR)","Partial Response (PR)",
                 "Stable Disease (SD)","Progressive Disease (PD)"),
      Deauville=c("1-2","3-4","3 (no change)","4-5 (new)"),
      SPD_Change=c("<CR threshold (>95%)","≥50% reduction","<50% red, <50% inc",">50% increase"),
      Model_Threshold=c("<5% of T0","5-50% of T0","50-150%",">150% of T0")
    )
  }, options=list(dom='t'), rownames=FALSE)
  output$plot_pfs_km <- renderPlot({
    pop <- pop_sim()
    event_times <- pop %>%
      filter(PD_flag==1) %>%
      group_by(ID) %>%
      summarise(event_time=min(time_day), event=1, .groups="drop")
    all_ids <- data.frame(ID=1:50)
    km_df <- left_join(all_ids, event_times, by="ID") %>%
      mutate(event=replace_na(event,0),
             event_time=replace_na(event_time, max(pop$time_day)))
    times <- sort(unique(c(0, km_df$event_time)))
    surv <- sapply(times, function(t) mean(km_df$event_time >= t))
    ggplot(data.frame(time=times, surv=surv), aes(x=time,y=surv*100)) +
      geom_step(color="#1565C0", linewidth=1.5) +
      geom_hline(yintercept=50, linetype="dashed", color="red") +
      labs(title="Estimated PFS (KM-like, N=50 Virtual Patients, R-CHOP)",
           x="Time (days)", y="PFS (%)") +
      ylim(0,100) + theme_dash
  })

  # ── Tab 5: Scenario Comparison ───────────────────────────────────────────
  output$plot_scenarios_tumor <- renderPlot({
    df <- all_scenarios()
    ggplot(df, aes(x=time_day, y=SPD_pct, color=Scenario)) +
      geom_line(linewidth=1.2) +
      geom_hline(yintercept=5, linetype="dashed", color="#4CAF50", linewidth=0.8) +
      scale_color_brewer(palette="Dark2") +
      labs(title="Tumor Burden Comparison — 6 Treatment Scenarios",
           x="Time (days)", y="SPD (% of baseline)") +
      ylim(0,NA) + theme_dash
  })
  output$dt_scenarios_table <- renderDT({
    data.frame(
      Scenario=c("1.Untreated","2.R-CHOP×6","3.Pola-R-CHP","4.R-CHOP+Ibr","5.R-CHOP+VEN","6.R-CHOP DHL"),
      Subtype=c("GCB","GCB","GCB","ABC","GCB","DHL"),
      Key_Drug=c("None","Rituximab+CHOP","Pola+R-CHP","Ibrutinib+R-CHOP","Venetoclax+R-CHOP","R-CHOP"),
      Ref=c("—","Coiffier 2002","Tilly 2022","Younes 2019","Morschhauser 2021","Dunleavy 2013")
    )
  }, options=list(dom='t',pageLength=6), rownames=FALSE)
  output$plot_eot_response <- renderPlot({
    df <- all_scenarios() %>%
      filter(abs(time_day-126)<3) %>%
      group_by(Scenario) %>%
      summarise(CR=mean(CR_flag)*100, PR=mean(PR_flag)*100,
                PD=mean(PD_flag)*100, .groups="drop") %>%
      pivot_longer(-Scenario, names_to="Response", values_to="pct")
    ggplot(df, aes(x=Scenario, y=pct, fill=Response)) +
      geom_bar(stat="identity", position="stack") +
      scale_fill_manual(values=c("CR"="#4CAF50","PR"="#8BC34A","PD"="#F44336")) +
      labs(title="Day-126 Response Rates by Scenario",
           x="", y="Patients (%)") +
      coord_flip() + theme_dash
  })
  output$plot_anc_scenarios <- renderPlot({
    df <- all_scenarios() %>% filter(Scenario != "1. Untreated")
    ggplot(df, aes(x=time_day, y=ANC_out, color=Scenario)) +
      geom_line(linewidth=1.0) +
      geom_hline(yintercept=0.5, linetype="dashed", color="red") +
      scale_color_brewer(palette="Set2") +
      labs(title="ANC Dynamics by Scenario", x="Time (days)", y="ANC (×10⁹/L)") +
      theme_dash
  })
  output$plot_resistance <- renderPlot({
    df <- all_scenarios()
    ggplot(df, aes(x=time_day, y=Resistance*100, color=Scenario)) +
      geom_line(linewidth=1.0) +
      scale_color_brewer(palette="Dark2") +
      labs(title="Resistance Development Over Time",
           x="Time (days)", y="Resistance fraction (%)") +
      theme_dash
  })

  # ── Tab 6: Biomarkers & Toxicity ──────────────────────────────────────────
  output$vbox_anc_nadir <- renderValueBox({
    df <- sim_data()
    anc_min <- min(df$ANC_out)
    valueBox(round(anc_min,2), "ANC Nadir (×10⁹/L)", icon=icon("syringe"), color="red")
  })
  output$vbox_anc_grade <- renderValueBox({
    df <- sim_data()
    anc_min <- min(df$ANC_out)
    grade <- if(anc_min<0.5)"Grade 4" else if(anc_min<1.0)"Grade 3" else
             if(anc_min<1.5)"Grade 2" else "Grade 1"
    valueBox(grade, "Worst Neutropenia Grade", icon=icon("exclamation-triangle"), color="orange")
  })
  output$vbox_crs_risk <- renderValueBox({
    df <- sim_data()
    crs_max <- max(df$CRS_risk)
    valueBox(round(crs_max,3), "Max CRS Risk Index", icon=icon("fire"), color="purple")
  })
  output$vbox_bcl2_occ <- renderValueBox({
    df <- sim_data()
    if(input$add_venetoclax) {
      occ_max <- max(df$BCL2_occ, na.rm=TRUE)
      valueBox(paste0(round(occ_max*100,1),"%"), "Peak BCL-2 Occupancy",
               icon=icon("lock"), color="green")
    } else {
      valueBox("N/A", "BCL-2 Occupancy (Venetoclax not selected)",
               icon=icon("lock-open"), color="gray")
    }
  })
  output$plot_anc_biomarker <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x=time_day, y=ANC_out)) +
      geom_line(color="#1565C0", linewidth=1.5) +
      geom_ribbon(aes(ymin=pmin(ANC_out,0.5), ymax=0.5), alpha=0.15, fill="red") +
      geom_hline(yintercept=c(0.5,1.0,1.5), linetype="dashed",
                 color=c("red","orange","gold")) +
      labs(title="ANC — Myelosuppression (G-CSF support typically if <0.5)",
           x="Time (days)", y="ANC (×10⁹/L)") +
      theme_dash
  })
  output$plot_bcl2_occ <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x=time_day, y=BCL2_occ*100)) +
      geom_line(color="#7B1FA2", linewidth=1.5) +
      labs(title="BCL-2 Occupancy (Venetoclax PD)",
           x="Time (days)", y="BCL-2 Occupancy (%)") +
      ylim(0,100) + theme_dash
  })
  output$plot_crs <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(x=time_day, y=CRS_risk)) +
      geom_line(color="#B71C1C", linewidth=1.4) +
      labs(title="CRS Risk Index Over Time",
           x="Time (days)", y="CRS Risk Index (0-1)") +
      ylim(0,1) + theme_dash
  })
  output$dt_biomarkers <- renderDT({
    df <- sim_data()
    checkpoints <- c(14,28,42,84,126)
    checkpoints <- checkpoints[checkpoints <= input$sim_days]
    tab <- lapply(checkpoints, function(tp) {
      row <- df %>% filter(abs(time_day-tp)<3) %>% slice_tail(n=1)
      data.frame(Day=tp,
                 ANC=round(row$ANC_out,2),
                 NK=round(row$NK_cells,2),
                 CD8=round(row$CD8_cells,2),
                 BCL2_Occ_pct=round(row$BCL2_occ*100,1),
                 Resist_pct=round(row$Resistance*100,2),
                 CRS_idx=round(row$CRS_risk,4))
    })
    do.call(rbind, tab)
  }, options=list(dom='t'), rownames=FALSE)
  output$plot_dox_cumulative <- renderPlot({
    df <- sim_data()
    bsa <- input$bsa_m2
    n_cy <- input$n_cycles
    dox_per_cycle <- input$dox_per_m2 * bsa
    cum_dox <- data.frame(
      cycle = 1:n_cy,
      day   = (1:n_cy - 1) * 21,
      cum   = dox_per_cycle * (1:n_cy)
    )
    p <- ggplot(cum_dox, aes(x=day, y=cum)) +
      geom_step(color="#E65100", linewidth=1.5) +
      geom_hline(yintercept=400, linetype="dashed", color="red") +
      annotate("text", x=max(cum_dox$day)*0.8, y=410,
               label="Cardiotoxicity threshold (~400 mg/m²)", color="red", size=3.5) +
      labs(title="Cumulative Doxorubicin Dose (Cardiotoxicity Monitor)",
           x="Time (days)", y=paste0("Cumulative dose (mg, BSA=",bsa,"m²)")) +
      theme_dash
    print(p)
  })
}

# ── Launch ───────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
