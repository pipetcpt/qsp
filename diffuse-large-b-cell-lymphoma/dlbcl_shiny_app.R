################################################################################
# DLBCL QSP Shiny Dashboard
# Disease: Diffuse Large B-Cell Lymphoma (DLBCL)
# Tabs: (1) Patient Profile · (2) Drug PK · (3) Oncogenic Signaling ·
#        (4) Clinical Endpoints · (5) Scenario Comparison · (6) TME & Biomarkers
################################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

# ── Inline model code (identical kernel to mrgsolve model file) ────────────
dlbcl_code <- '
$PARAM
  kBCR_syn=0.10, kBCR_deg=0.05, kSYK_act=0.30, kSYK_deg=0.40,
  kBTK_act=0.25, kBTK_deg=0.35, BTK_base=1.0,
  kNFkB_act=0.20, kNFkB_deg=0.30, NFkB_base=0.5,
  kAKT_act=0.15, kAKT_deg=0.25, AKT_base=0.8,
  kmTOR_act=0.12, kmTOR_deg=0.18,
  kSTAT3_act=0.08, kSTAT3_deg=0.15, IL6_base=0.3,
  BCL2_base=1.5, kBCL2_syn=0.10, kBCL2_deg=0.04,
  kBIM_syn=0.05, kBIM_deg=0.08, kApop=0.30, kApop_base=0.002, BCL2_thresh=1.0,
  kGCB_grow=0.025, kGCB_death=0.005, GCB_cap=1000.0, GCB_0=50.0,
  kABC_grow=0.030, kABC_death=0.004, ABC_cap=1000.0, ABC_0=0.0,
  kCD8_base=0.02, kCD8_deg=0.015, kNK_base=0.015, kNK_deg=0.012,
  kKill_CD8=0.04, kKill_NK=0.03, kExhaust=0.005,
  kLDH_syn=0.05, kLDH_deg=0.03, LDH_norm=1.0,
  kctDNA_shed=0.02, kctDNA_cl=0.08,
  RTX_CL=0.010, RTX_Vc=2.80, RTX_Q=0.005, RTX_Vp=4.00,
  RTX_Emax=0.90, RTX_EC50=5.00,
  CYC_ka=0.80, CYC_CL=8.00, CYC_Vd=35.0, CYC_Emax=0.70, CYC_EC50=0.50,
  DOX_CL=4.00, DOX_Vd=700.0, DOX_Emax=0.80, DOX_EC50=0.05,
  VCR_CL=5.00, VCR_Vd=560.0, VCR_Emax=0.50, VCR_EC50=0.01,
  POLA_CL=0.015, POLA_Vc=3.50, POLA_Q=0.008, POLA_Vp=5.00,
  POLA_krel=0.020, POLA_Emax=0.85, POLA_EC50=2.00,
  VEN_ka=0.50, VEN_CL=7.00, VEN_Vd=256.0, VEN_Emax=0.75, VEN_EC50=0.30,
  CART_kexp=0.040, CART_kprol=0.030, CART_kdeath=0.008,
  CART_kexh=0.003, CART_Emax=0.95, CART_EC50=10.0,
  give_RCHOP=0, give_PolaRC=0, give_VEN=0, give_CART=0, give_IBRU=0,
  COO_ABC=0, DHL=0, TP53_mut=0, IPI=2

$CMT
  RTX1 RTX2 CYC_p CYC4OH DOX VCR POLA1 POLA2 POLA_MMAE VEN_gut VEN_plasma
  pSYK pBTK pNFkB pAKT pSTAT3
  BCL2_prot BIM_prot Apoptosis
  GCB_tumor ABC_tumor
  CD8_T NK_cell CART_cells
  LDH_level ctDNA_level

$MAIN
  double ABC_boost  = (COO_ABC==1) ? 1.4 : 1.0;
  double DHL_boost  = (DHL==1)     ? 1.6 : 1.0;
  double p53_factor = (TP53_mut==1)? 0.4 : 1.0;
  RTX1_0=0; RTX2_0=0; CYC_p_0=0; CYC4OH_0=0; DOX_0=0; VCR_0=0;
  POLA1_0=0; POLA2_0=0; POLA_MMAE_0=0; VEN_gut_0=0; VEN_plasma_0=0;
  pSYK_0=0.5; pBTK_0=BTK_base*ABC_boost; pNFkB_0=NFkB_base*ABC_boost;
  pAKT_0=AKT_base; pSTAT3_0=0.3*ABC_boost;
  BCL2_prot_0=BCL2_base*(1+0.5*(DHL==1)); BIM_prot_0=0.5; Apoptosis_0=0;
  GCB_tumor_0=GCB_0*(1-COO_ABC); ABC_tumor_0=ABC_0*COO_ABC;
  CD8_T_0=2.0; NK_cell_0=1.5; CART_cells_0=0;
  LDH_level_0=LDH_norm+0.05*GCB_0+0.05*ABC_0;
  ctDNA_level_0=0.01*(GCB_0+ABC_0);

$ODE
  double RTX_conc=RTX1/RTX_Vc; double CYC_act_conc=CYC4OH/CYC_Vd;
  double DOX_conc=DOX/DOX_Vd; double VCR_conc=VCR/VCR_Vd;
  double POLA_conc=POLA1/POLA_Vc; double MMAE_conc=POLA_MMAE;
  double VEN_conc=VEN_plasma/VEN_Vd;
  dxdt_RTX1=-(RTX_CL/RTX_Vc)*RTX1-(RTX_Q/RTX_Vc)*RTX1+(RTX_Q/RTX_Vp)*RTX2;
  dxdt_RTX2=(RTX_Q/RTX_Vc)*RTX1-(RTX_Q/RTX_Vp)*RTX2;
  dxdt_CYC_p=-CYC_ka*CYC_p;
  dxdt_CYC4OH=CYC_ka*CYC_p-(CYC_CL/CYC_Vd)*CYC4OH;
  dxdt_DOX=-(DOX_CL/DOX_Vd)*DOX;
  dxdt_VCR=-(VCR_CL/VCR_Vd)*VCR;
  dxdt_POLA1=-(POLA_CL/POLA_Vc)*POLA1-(POLA_Q/POLA_Vc)*POLA1+(POLA_Q/POLA_Vp)*POLA2;
  dxdt_POLA2=(POLA_Q/POLA_Vc)*POLA1-(POLA_Q/POLA_Vp)*POLA2;
  dxdt_POLA_MMAE=POLA_krel*POLA1-0.04*POLA_MMAE;
  dxdt_VEN_gut=-VEN_ka*VEN_gut;
  dxdt_VEN_plasma=VEN_ka*VEN_gut-(VEN_CL/VEN_Vd)*VEN_plasma;
  double IBRU_inh=(give_IBRU==1)?0.85:0.0;
  double ABC_b=(COO_ABC==1)?1.4:1.0;
  double DHL_b=(DHL==1)?1.6:1.0;
  double p53_f=(TP53_mut==1)?0.4:1.0;
  dxdt_pSYK=kSYK_act*(1.0+0.5*ABC_b)-kSYK_deg*pSYK;
  dxdt_pBTK=kBTK_act*pSYK*(1-IBRU_inh)-kBTK_deg*pBTK;
  dxdt_pNFkB=kNFkB_act*pBTK*ABC_b-kNFkB_deg*pNFkB+kNFkB_act*0.1*ABC_b;
  dxdt_pAKT=kAKT_act*(0.5+pSYK*0.5)-kAKT_deg*pAKT;
  double IL6_p=IL6_base*ABC_b*(1+0.2*pNFkB);
  dxdt_pSTAT3=kSTAT3_act*IL6_p*(1+0.5*ABC_b)-kSTAT3_deg*pSTAT3;
  double BCL2_d=kBCL2_syn*(1+0.5*pNFkB+0.3*pSTAT3)*DHL_b;
  double VEN_eff=(give_VEN==1)?(VEN_Emax*VEN_conc/(VEN_EC50+VEN_conc)):0.0;
  double EffBCL2=BCL2_prot*(1-VEN_eff);
  dxdt_BCL2_prot=BCL2_d-kBCL2_deg*BCL2_prot;
  dxdt_BIM_prot=kBIM_syn-kBIM_deg*BIM_prot;
  double freeApop=BIM_prot/(1+EffBCL2/BCL2_thresh);
  double apop_r=kApop*freeApop*p53_f+kApop_base;
  dxdt_Apoptosis=apop_r;
  double total_t=GCB_tumor+ABC_tumor;
  double RTX_k=(give_RCHOP==1||give_PolaRC==1)?RTX_Emax*RTX_conc/(RTX_EC50+RTX_conc):0.0;
  double CYC_k=(give_RCHOP==1||give_PolaRC==1)?CYC_Emax*CYC_act_conc/(CYC_EC50+CYC_act_conc):0.0;
  double DOX_k=(give_RCHOP==1)?DOX_Emax*DOX_conc/(DOX_EC50+DOX_conc):0.0;
  double VCR_k=(give_RCHOP==1||give_PolaRC==1)?VCR_Emax*VCR_conc/(VCR_EC50+VCR_conc):0.0;
  double POLA_k=(give_PolaRC==1)?POLA_Emax*MMAE_conc/(POLA_EC50+MMAE_conc):0.0;
  double VEN_s=VEN_eff*0.5;
  double dk=RTX_k+CYC_k+DOX_k+VCR_k+POLA_k+VEN_s;
  if(dk>0.99) dk=0.99;
  double CD8_kt=kKill_CD8*CD8_T*GCB_tumor/(GCB_tumor+1.0);
  double NK_kt=kKill_NK*NK_cell*GCB_tumor/(GCB_tumor+1.0);
  double CART_k=(give_CART==1)?CART_Emax*CART_cells/(CART_EC50+CART_cells):0.0;
  dxdt_GCB_tumor=kGCB_grow*DHL_b*GCB_tumor*(1-total_t/GCB_cap)
                 -(kGCB_death+dk+apop_r*0.5+CART_k)*GCB_tumor
                 -(CD8_kt+NK_kt)/GCB_cap*GCB_tumor;
  double CD8_ka=kKill_CD8*CD8_T*ABC_tumor/(ABC_tumor+1.0);
  double NK_ka=kKill_NK*NK_cell*ABC_tumor/(ABC_tumor+1.0);
  double IBRU_k=(give_IBRU==1)?0.40:0.0;
  dxdt_ABC_tumor=kABC_grow*DHL_b*ABC_b*ABC_tumor*(1-total_t/ABC_cap)
                 -(kABC_death+dk+apop_r*0.5+IBRU_k+CART_k)*ABC_tumor
                 -(CD8_ka+NK_ka)/ABC_cap*ABC_tumor;
  double PDL1=0.3+0.4*pNFkB+0.3*pSTAT3;
  dxdt_CD8_T=kCD8_base*(1+0.3*RTX_k)-kCD8_deg*CD8_T-kExhaust*CD8_T*PDL1;
  dxdt_NK_cell=kNK_base*(1+0.5*RTX_k)-kNK_deg*NK_cell;
  dxdt_CART_cells=(give_CART==1)?
    CART_kexp*CART_cells+CART_kprol*total_t*CART_cells/(total_t+10.0)
    -CART_kdeath*CART_cells-CART_kexh*CART_cells : 0.0;
  dxdt_LDH_level=kLDH_syn*total_t-kLDH_deg*(LDH_level-LDH_norm);
  dxdt_ctDNA_level=kctDNA_shed*total_t-kctDNA_cl*ctDNA_level;

$TABLE
  double CONC_RTX=RTX1/RTX_Vc; double CONC_CYC=CYC4OH/CYC_Vd;
  double CONC_DOX=DOX/DOX_Vd; double CONC_VCR=VCR/VCR_Vd;
  double CONC_POLA=POLA1/POLA_Vc; double CONC_VEN=VEN_plasma/VEN_Vd;
  double CONC_MMAE=POLA_MMAE;
  double TumorTotal=GCB_tumor+ABC_tumor;
  double PercentRedux=100.0*(1.0-TumorTotal/(GCB_0+ABC_0+1e-6));
  double TMTV_proxy=TumorTotal*2.5;

$CAPTURE CONC_RTX CONC_CYC CONC_DOX CONC_VCR CONC_POLA CONC_VEN CONC_MMAE
         pSYK pBTK pNFkB pAKT pSTAT3
         BCL2_prot BIM_prot Apoptosis
         GCB_tumor ABC_tumor TumorTotal PercentRedux TMTV_proxy
         CD8_T NK_cell CART_cells LDH_level ctDNA_level
'

mod_base <- mcode("dlbcl_shiny", dlbcl_code, quiet = TRUE)

# Helper: build events
make_events <- function(regimen, bsa, weight, n_cycles, iv_start, ven_days) {
  interval <- 21 * 24
  starts   <- seq(iv_start, iv_start + (n_cycles - 1) * interval, by = interval)
  evlist   <- list()

  rtx <- 375 * bsa; cyc <- 750 * bsa; dox <- 50 * bsa
  vcr <- 1.4 * bsa; pola <- 1.8 * weight

  if (regimen %in% c("R-CHOP", "Pola-R-CHP", "R-CHOP + VEN")) {
    evlist[[1]] <- tibble(time = starts, cmt = "RTX1",  amt = rtx, evid = 1, rate = -2)
    evlist[[2]] <- tibble(time = starts, cmt = "CYC_p", amt = cyc, evid = 1, rate = -2)
    evlist[[3]] <- tibble(time = starts, cmt = "VCR",   amt = vcr, evid = 1, rate = -2)
  }
  if (regimen %in% c("R-CHOP", "R-CHOP + VEN")) {
    evlist[[4]] <- tibble(time = starts, cmt = "DOX", amt = dox, evid = 1, rate = -2)
  }
  if (regimen == "Pola-R-CHP") {
    evlist[[5]] <- tibble(time = starts, cmt = "POLA1", amt = pola, evid = 1, rate = -2)
  }
  if (regimen %in% c("Venetoclax", "R-CHOP + VEN")) {
    ven_times <- seq(iv_start, iv_start + ven_days * 24 - 24, by = 24)
    evlist[[6]] <- tibble(time = ven_times, cmt = "VEN_gut", amt = 400, evid = 1)
  }
  if (regimen == "CAR-T") {
    evlist[[7]] <- tibble(time = iv_start, cmt = "CART_cells", amt = 1.4e4, evid = 1)
  }

  all_ev <- bind_rows(evlist) %>% arrange(time)
  if (nrow(all_ev) == 0) return(ev(time = 0, cmt = 1, amt = 0, evid = 2))
  as.ev(all_ev)
}

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "DLBCL QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "tab1", icon = icon("user-md")),
      menuItem("Drug PK",            tabName = "tab2", icon = icon("capsules")),
      menuItem("Oncogenic Signaling",tabName = "tab3", icon = icon("dna")),
      menuItem("Clinical Endpoints", tabName = "tab4", icon = icon("chart-line")),
      menuItem("Scenario Comparison",tabName = "tab5", icon = icon("exchange-alt")),
      menuItem("TME & Biomarkers",   tabName = "tab6", icon = icon("microscope"))
    ),
    hr(),
    h5("  Treatment", style = "color:#ccc; margin-left:10px"),
    selectInput("regimen", "Regimen",
                choices = c("R-CHOP", "Pola-R-CHP", "Venetoclax", "CAR-T", "R-CHOP + VEN", "Untreated"),
                selected = "R-CHOP"),
    sliderInput("n_cycles", "# of cycles (IV)", 1, 8, 6),
    sliderInput("ven_days", "Venetoclax duration (days)", 30, 365, 180),
    hr(),
    h5("  Patient", style = "color:#ccc; margin-left:10px"),
    selectInput("coo", "Cell-of-Origin",
                choices = c("GCB" = 0, "ABC" = 1), selected = 0),
    checkboxInput("dhl",     "Double-hit (MYC + BCL2/BCL6)", FALSE),
    checkboxInput("tp53",    "TP53 mutation", FALSE),
    checkboxInput("ibrutinib","Add Ibrutinib (ABC)", FALSE),
    sliderInput("ipi", "IPI Score", 0, 5, 2),
    sliderInput("tumor0", "Initial Tumor Burden (AU)", 10, 200, 50),
    sliderInput("bsa", "BSA (m²)", 1.2, 2.5, 1.8, step = 0.1),
    sliderInput("weight","Weight (kg)", 40, 120, 70),
    sliderInput("sim_days","Simulation (days)", 90, 365, 180),
    actionButton("run_sim", "Run Simulation", icon = icon("play"),
                 class = "btn-success", width = "90%")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #F4F6F9; }
      .box { border-top: 3px solid #6C3483; }
    "))),
    tabItems(

      # ── TAB 1: Patient Profile ──────────────────────────────────────────
      tabItem("tab1",
        fluidRow(
          valueBoxOutput("vb_regimen", width = 4),
          valueBoxOutput("vb_coo",     width = 4),
          valueBoxOutput("vb_ipi",     width = 4)
        ),
        fluidRow(
          box(title = "Disease Characteristics", width = 6, solidHeader = TRUE,
              status = "purple",
              tableOutput("pt_summary")),
          box(title = "Mechanistic Pathway Map — DLBCL", width = 6, solidHeader = TRUE,
              status = "purple",
              p(strong("Key pathways in DLBCL:")),
              tags$ul(
                tags$li("BCR → SYK → BTK → NF-κB (ABC: MYD88 L265P, constitutive)"),
                tags$li("PI3K/AKT/mTOR: cell survival & translation"),
                tags$li("JAK/STAT3: IL-6/IL-10 autocrine loop, PD-L1 induction"),
                tags$li("MYC/BCL6/EZH2: GCB epigenetic driver"),
                tags$li("BCL-2 family: BCL2 translocation t(14;18) in GCB;"),
                tags$li("TME: CD8 exhaustion by PD-L1; NK ADCC by rituximab"),
                tags$li("Double-hit: MYC + BCL2/BCL6 → R-CHOP resistance")
              ),
              p(em("Cell-of-Origin (COO): GCB better prognosis; ABC worse (NF-κB driven)."))
          )
        ),
        fluidRow(
          box(title = "IPI Risk Calculator", width = 6, solidHeader = TRUE,
              status = "warning",
              p("International Prognostic Index (IPI) components:"),
              tableOutput("ipi_table")),
          box(title = "DLBCL Subtypes", width = 6, solidHeader = TRUE, status = "info",
              plotOutput("coo_radar", height = "280px"))
        )
      ),

      # ── TAB 2: Drug PK ──────────────────────────────────────────────────
      tabItem("tab2",
        fluidRow(
          box(title = "Rituximab PK (2-Compartment)", width = 6, solidHeader = TRUE,
              status = "purple", plotlyOutput("pk_rtx", height = "300px")),
          box(title = "Cytotoxic Agents PK (1-Compartment)", width = 6, solidHeader = TRUE,
              status = "purple", plotlyOutput("pk_chemo", height = "300px"))
        ),
        fluidRow(
          box(title = "Polatuzumab Vedotin + MMAE Payload", width = 6, solidHeader = TRUE,
              status = "purple", plotlyOutput("pk_pola", height = "300px")),
          box(title = "Venetoclax Oral PK", width = 6, solidHeader = TRUE,
              status = "purple", plotlyOutput("pk_ven", height = "300px"))
        ),
        fluidRow(
          box(title = "PK Summary Table (Day 21 Trough/Peak)", width = 12,
              solidHeader = TRUE, status = "info", DTOutput("pk_table"))
        )
      ),

      # ── TAB 3: Oncogenic Signaling ───────────────────────────────────────
      tabItem("tab3",
        fluidRow(
          box(title = "BCR / BTK / NF-κB Dynamics", width = 6, solidHeader = TRUE,
              status = "purple", plotlyOutput("sig_nfkb", height = "300px")),
          box(title = "PI3K / AKT and JAK / STAT3", width = 6, solidHeader = TRUE,
              status = "purple", plotlyOutput("sig_akt", height = "300px"))
        ),
        fluidRow(
          box(title = "BCL-2 vs BIM Balance (Apoptosis)", width = 6, solidHeader = TRUE,
              status = "danger", plotlyOutput("sig_bcl2", height = "300px")),
          box(title = "Cumulative Apoptosis Index", width = 6, solidHeader = TRUE,
              status = "danger", plotlyOutput("sig_apop", height = "300px"))
        )
      ),

      # ── TAB 4: Clinical Endpoints ─────────────────────────────────────────
      tabItem("tab4",
        fluidRow(
          valueBoxOutput("vb_response",  width = 3),
          valueBoxOutput("vb_reduction", width = 3),
          valueBoxOutput("vb_tmtv",      width = 3),
          valueBoxOutput("vb_ldh",       width = 3)
        ),
        fluidRow(
          box(title = "Tumor Burden Over Time", width = 6, solidHeader = TRUE,
              status = "purple", plotlyOutput("ep_tumor", height = "300px")),
          box(title = "% Tumor Reduction (Waterfall at Day 90)", width = 6,
              solidHeader = TRUE, status = "purple",
              plotOutput("ep_waterfall", height = "300px"))
        ),
        fluidRow(
          box(title = "GCB vs ABC Subtype Tumor Dynamics", width = 6, solidHeader = TRUE,
              status = "info", plotlyOutput("ep_coo_split", height = "300px")),
          box(title = "Simulated Progression-Free Survival Proxy", width = 6,
              solidHeader = TRUE, status = "success",
              plotlyOutput("ep_pfs", height = "300px"))
        )
      ),

      # ── TAB 5: Scenario Comparison ──────────────────────────────────────
      tabItem("tab5",
        fluidRow(
          box(title = "Multi-Regimen Tumor Burden Comparison", width = 8,
              solidHeader = TRUE, status = "purple",
              plotlyOutput("sc_tumor", height = "350px")),
          box(title = "Day-180 Response Summary", width = 4, solidHeader = TRUE,
              status = "info", DTOutput("sc_table"))
        ),
        fluidRow(
          box(title = "BCR Signaling Comparison (R-CHOP vs Pola-R-CHP vs CAR-T)", width = 6,
              solidHeader = TRUE, status = "purple",
              plotlyOutput("sc_nfkb", height = "300px")),
          box(title = "% Reduction: R-CHOP vs Pola-R-CHP (POLARIX-like)", width = 6,
              solidHeader = TRUE, status = "success",
              plotlyOutput("sc_compare", height = "300px"))
        )
      ),

      # ── TAB 6: TME & Biomarkers ──────────────────────────────────────────
      tabItem("tab6",
        fluidRow(
          box(title = "CD8+ T cells and NK cells", width = 6, solidHeader = TRUE,
              status = "purple", plotlyOutput("tme_immune", height = "300px")),
          box(title = "CAR-T Cell Kinetics (axi-cel)", width = 6, solidHeader = TRUE,
              status = "warning", plotlyOutput("tme_cart", height = "300px"))
        ),
        fluidRow(
          box(title = "Serum LDH (surrogate for tumor burden)", width = 6,
              solidHeader = TRUE, status = "info",
              plotlyOutput("bm_ldh", height = "300px")),
          box(title = "ctDNA Dynamics", width = 6, solidHeader = TRUE,
              status = "info", plotlyOutput("bm_ctdna", height = "300px"))
        ),
        fluidRow(
          box(title = "PD-L1 Proxy (NF-κB + STAT3 driven)", width = 6, solidHeader = TRUE,
              status = "danger", plotlyOutput("bm_pdl1", height = "300px")),
          box(title = "BCL-2 / BIM Ratio over Time", width = 6, solidHeader = TRUE,
              status = "danger", plotlyOutput("bm_bcl2ratio", height = "300px"))
        )
      )
    )
  )
)

# ── SERVER ─────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Reactive: run simulation
  sim_result <- eventReactive(input$run_sim, {
    coo_val  <- as.numeric(input$coo)
    gcb_init <- if (coo_val == 0) input$tumor0 else 0
    abc_init <- if (coo_val == 1) input$tumor0 else 0

    ev_drug <- make_events(input$regimen, input$bsa, input$weight,
                           input$n_cycles, 0, input$ven_days)

    p_override <- list(
      COO_ABC = coo_val, DHL = as.numeric(input$dhl),
      TP53_mut = as.numeric(input$tp53), IPI = input$ipi,
      give_IBRU = as.numeric(input$ibrutinib),
      give_RCHOP  = as.numeric(input$regimen %in% c("R-CHOP", "R-CHOP + VEN")),
      give_PolaRC = as.numeric(input$regimen == "Pola-R-CHP"),
      give_VEN    = as.numeric(input$regimen %in% c("Venetoclax", "R-CHOP + VEN")),
      give_CART   = as.numeric(input$regimen == "CAR-T"),
      GCB_0 = gcb_init, ABC_0 = abc_init
    )
    p_mod <- do.call(param, c(list(mod_base), p_override))

    out <- mrgsim(p_mod, events = ev_drug, obsonly = TRUE,
                  end = input$sim_days * 24, delta = 6)
    as_tibble(out) %>% mutate(day = time / 24)
  }, ignoreNULL = FALSE)

  # Also run a multi-scenario comparison (fixed parameters)
  multi_sim <- eventReactive(input$run_sim, {
    coo_val  <- as.numeric(input$coo)
    gcb_init <- if (coo_val == 0) input$tumor0 else 0
    abc_init <- if (coo_val == 1) input$tumor0 else 0
    common <- list(COO_ABC = coo_val, DHL = as.numeric(input$dhl),
                   TP53_mut = as.numeric(input$tp53), IPI = input$ipi,
                   GCB_0 = gcb_init, ABC_0 = abc_init)

    run_sc <- function(label, reg) {
      ev <- make_events(reg, input$bsa, input$weight, input$n_cycles, 0, input$ven_days)
      ovr <- c(common, list(
        give_RCHOP  = as.numeric(reg %in% c("R-CHOP", "R-CHOP + VEN")),
        give_PolaRC = as.numeric(reg == "Pola-R-CHP"),
        give_VEN    = as.numeric(reg %in% c("Venetoclax", "R-CHOP + VEN")),
        give_CART   = as.numeric(reg == "CAR-T"),
        give_IBRU   = 0
      ))
      pm <- do.call(param, c(list(mod_base), ovr))
      mrgsim(pm, events = ev, obsonly = TRUE,
             end = input$sim_days * 24, delta = 24) %>%
        as_tibble() %>% mutate(day = time / 24, scenario = label)
    }

    bind_rows(
      run_sc("Untreated",   "Untreated"),
      run_sc("R-CHOP",      "R-CHOP"),
      run_sc("Pola-R-CHP",  "Pola-R-CHP"),
      run_sc("Venetoclax",  "Venetoclax"),
      run_sc("CAR-T",       "CAR-T")
    )
  }, ignoreNULL = FALSE)

  # Value boxes
  output$vb_regimen <- renderValueBox(
    valueBox(input$regimen, "Treatment Regimen", icon = icon("pills"), color = "purple"))
  output$vb_coo <- renderValueBox(
    valueBox(ifelse(input$coo == "0", "GCB", "ABC"), "Cell-of-Origin",
             icon = icon("dna"), color = if (input$coo == "0") "blue" else "red"))
  output$vb_ipi <- renderValueBox(
    valueBox(paste("IPI =", input$ipi),
             ifelse(input$ipi <= 1, "Low Risk", ifelse(input$ipi <= 2, "Low-Intermediate",
                    ifelse(input$ipi <= 3, "High-Intermediate", "High Risk"))),
             icon = icon("chart-bar"),
             color = ifelse(input$ipi <= 1, "green", ifelse(input$ipi <= 2, "yellow",
                            ifelse(input$ipi <= 3, "orange", "red")))))

  output$vb_response <- renderValueBox({
    d <- sim_result()
    last <- tail(d, 1)
    resp <- if (last$PercentRedux >= 95) "CR" else
            if (last$PercentRedux >= 50) "PR" else
            if (last$PercentRedux >= -25) "SD" else "PD"
    col  <- c(CR="green", PR="blue", SD="yellow", PD="red")[resp]
    valueBox(resp, "Best Response", icon = icon("check-circle"), color = col)
  })
  output$vb_reduction <- renderValueBox({
    d   <- sim_result()
    val <- round(max(d$PercentRedux, na.rm = TRUE), 1)
    valueBox(paste0(val, "%"), "Max Tumor Reduction", icon = icon("arrow-down"), color = "teal")
  })
  output$vb_tmtv <- renderValueBox({
    d   <- sim_result()
    val <- round(tail(d$TMTV_proxy, 1), 0)
    valueBox(paste0(val, " cm³"), "TMTV Proxy (Day-end)", icon = icon("cube"), color = "maroon")
  })
  output$vb_ldh <- renderValueBox({
    d   <- sim_result()
    val <- round(tail(d$LDH_level, 1), 2)
    valueBox(paste0(val, " AU"), "LDH (Day-end)", icon = icon("vial"), color = "orange")
  })

  # Patient summary table
  output$pt_summary <- renderTable({
    tibble(
      Parameter = c("COO", "Double-hit", "TP53 mutation", "Initial Tumor Burden",
                    "IPI Score", "Ibrutinib add-on", "BSA", "Weight"),
      Value = c(ifelse(input$coo == "0", "GCB", "ABC"),
                ifelse(input$dhl, "Yes", "No"),
                ifelse(input$tp53, "Yes", "No"),
                paste(input$tumor0, "AU"),
                paste(input$ipi, "/5"),
                ifelse(input$ibrutinib, "Yes (ABC benefit)", "No"),
                paste(input$bsa, "m²"),
                paste(input$weight, "kg"))
    )
  })

  # IPI table
  output$ipi_table <- renderTable({
    tibble(
      Factor = c("Age > 60", "Elevated LDH", "ECOG PS ≥ 2", "Stage III/IV", "Extranodal sites > 1"),
      Points = c(1, 1, 1, 1, 1),
      `Risk Group (IPI)` = c("0–1: Low (5yr OS ~73%)", "2: Low-int (~51%)",
                              "3: High-int (~43%)", "4–5: High (~26%)",
                              rep("", 1))
    )
  })

  # COO bar chart
  output$coo_radar <- renderPlot({
    tibble(
      Feature = c("NF-κB", "BCR dep.", "STAT3", "EZH2", "BCL2 transl.", "MYD88 L265P", "OS (5yr)"),
      GCB = c(0.3, 0.5, 0.4, 0.8, 0.7, 0.05, 73),
      ABC = c(0.9, 0.85, 0.75, 0.2, 0.15, 0.35, 40)
    ) %>%
      pivot_longer(-Feature) %>%
      mutate(value = ifelse(Feature == "OS (5yr)", value / 100, value)) %>%
      ggplot(aes(Feature, value, fill = name)) +
      geom_bar(stat = "identity", position = "dodge") +
      coord_flip() +
      scale_fill_manual(values = c(GCB = "#2980B9", ABC = "#C0392B")) +
      labs(title = "GCB vs ABC Feature Profile", y = "Relative score (0–1)",
           fill = "Subtype") +
      theme_minimal(base_size = 11)
  })

  # ── PK PLOTS ────────────────────────────────────────────────────────────
  pk_data <- reactive({
    sim_result() %>% filter(day <= min(42, input$sim_days))
  })

  output$pk_rtx <- renderPlotly({
    p <- pk_data() %>%
      ggplot(aes(day, CONC_RTX)) + geom_line(color = "#2E86C1", linewidth = 1.1) +
      labs(x = "Day", y = "Conc (μg/mL)", title = "Rituximab") + theme_bw()
    ggplotly(p)
  })
  output$pk_chemo <- renderPlotly({
    p <- pk_data() %>%
      select(day, CONC_CYC, CONC_DOX, CONC_VCR) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, color = name)) + geom_line(linewidth = 1) +
      labs(x = "Day", y = "Conc (AU)", title = "Cytotoxic Agents") + theme_bw()
    ggplotly(p)
  })
  output$pk_pola <- renderPlotly({
    p <- pk_data() %>%
      select(day, CONC_POLA, CONC_MMAE) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, color = name)) + geom_line(linewidth = 1) +
      labs(x = "Day", y = "Conc (AU/nM)", title = "Polatuzumab + MMAE") + theme_bw()
    ggplotly(p)
  })
  output$pk_ven <- renderPlotly({
    p <- pk_data() %>%
      ggplot(aes(day, CONC_VEN)) + geom_line(color = "#8E44AD", linewidth = 1.1) +
      labs(x = "Day", y = "Conc (μg/mL)", title = "Venetoclax") + theme_bw()
    ggplotly(p)
  })

  output$pk_table <- renderDT({
    sim_result() %>%
      filter(day %in% c(1, 7, 14, 21, 42)) %>%
      select(day, CONC_RTX, CONC_CYC, CONC_DOX, CONC_VCR, CONC_POLA, CONC_VEN) %>%
      mutate(across(-day, \(x) round(x, 4))) %>%
      datatable(options = list(pageLength = 5))
  })

  # ── SIGNALING PLOTS ──────────────────────────────────────────────────────
  output$sig_nfkb <- renderPlotly({
    p <- sim_result() %>%
      select(day, pSYK, pBTK, pNFkB) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, color = name)) + geom_line(linewidth = 1) +
      labs(x = "Day", y = "Phospho (AU)", title = "BCR → BTK → NF-κB") + theme_bw()
    ggplotly(p)
  })
  output$sig_akt <- renderPlotly({
    p <- sim_result() %>%
      select(day, pAKT, pSTAT3) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, color = name)) + geom_line(linewidth = 1) +
      labs(x = "Day", y = "Phospho (AU)") + theme_bw()
    ggplotly(p)
  })
  output$sig_bcl2 <- renderPlotly({
    p <- sim_result() %>%
      select(day, BCL2_prot, BIM_prot) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, color = name)) + geom_line(linewidth = 1) +
      labs(x = "Day", y = "Protein (AU)") + theme_bw()
    ggplotly(p)
  })
  output$sig_apop <- renderPlotly({
    p <- sim_result() %>%
      ggplot(aes(day, Apoptosis)) + geom_line(color = "#C0392B", linewidth = 1.2) +
      labs(x = "Day", y = "Cumulative Apoptosis (AU)") + theme_bw()
    ggplotly(p)
  })

  # ── CLINICAL ENDPOINTS ───────────────────────────────────────────────────
  output$ep_tumor <- renderPlotly({
    p <- sim_result() %>%
      ggplot(aes(day, TumorTotal)) + geom_line(color = "#6C3483", linewidth = 1.3) +
      geom_hline(yintercept = input$tumor0 * 0.05, linetype = "dashed", color = "green") +
      geom_hline(yintercept = input$tumor0 * 0.50, linetype = "dashed", color = "blue") +
      annotate("text", x = 5, y = input$tumor0*0.05, label = "CR", hjust=0, size=3, color="green") +
      annotate("text", x = 5, y = input$tumor0*0.50, label = "PR", hjust=0, size=3, color="blue") +
      labs(x = "Day", y = "Tumor Burden (AU)") + theme_bw()
    ggplotly(p)
  })
  output$ep_waterfall <- renderPlot({
    d <- sim_result() %>% filter(day == min(90, input$sim_days - 1))
    ggplot(data.frame(x = "Patient", y = d$PercentRedux[1]),
           aes(x, y, fill = y > 0)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = c(50, 100), linetype = "dashed") +
      scale_fill_manual(values = c("FALSE" = "#E74C3C", "TRUE" = "#2ECC71")) +
      labs(x = "", y = "% Reduction", title = "Waterfall (Day 90)") +
      theme_bw() + theme(legend.position = "none")
  })
  output$ep_coo_split <- renderPlotly({
    p <- sim_result() %>%
      select(day, GCB_tumor, ABC_tumor) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, color = name)) + geom_line(linewidth = 1) +
      labs(x = "Day", y = "Tumor Burden (AU)", color = "Subtype") + theme_bw()
    ggplotly(p)
  })
  output$ep_pfs <- renderPlotly({
    p <- sim_result() %>%
      mutate(surv_proxy = exp(-0.003 * TMTV_proxy)) %>%
      ggplot(aes(day, surv_proxy)) + geom_line(color = "#27AE60", linewidth = 1.3) +
      ylim(0, 1) +
      labs(x = "Day", y = "Survival Probability (proxy)") + theme_bw()
    ggplotly(p)
  })

  # ── SCENARIO COMPARISON ──────────────────────────────────────────────────
  output$sc_tumor <- renderPlotly({
    p <- multi_sim() %>%
      ggplot(aes(day, TumorTotal, color = scenario)) +
      geom_line(linewidth = 1.1) +
      labs(x = "Day", y = "Tumor Burden (AU)", color = "Regimen") + theme_bw()
    ggplotly(p)
  })
  output$sc_table <- renderDT({
    multi_sim() %>%
      group_by(scenario) %>%
      slice_tail(n = 1) %>%
      ungroup() %>%
      mutate(
        Response = case_when(
          PercentRedux >= 95 ~ "CR",
          PercentRedux >= 50 ~ "PR",
          PercentRedux >= -25 ~ "SD",
          TRUE ~ "PD"),
        `% Reduction` = round(PercentRedux, 1),
        `TMTV (cm³)` = round(TMTV_proxy, 0)
      ) %>%
      select(scenario, Response, `% Reduction`, `TMTV (cm³)`) %>%
      datatable(options = list(dom = "t"))
  })
  output$sc_nfkb <- renderPlotly({
    p <- multi_sim() %>%
      filter(scenario %in% c("R-CHOP", "Pola-R-CHP", "CAR-T", "Untreated")) %>%
      ggplot(aes(day, pNFkB, color = scenario)) + geom_line(linewidth = 1) +
      labs(x = "Day", y = "NF-κB (AU)", color = "Regimen") + theme_bw()
    ggplotly(p)
  })
  output$sc_compare <- renderPlotly({
    p <- multi_sim() %>%
      filter(scenario %in% c("R-CHOP", "Pola-R-CHP")) %>%
      ggplot(aes(day, PercentRedux, color = scenario)) + geom_line(linewidth = 1.2) +
      geom_hline(yintercept = 50, linetype = "dashed") +
      labs(x = "Day", y = "% Reduction", color = "Regimen",
           title = "R-CHOP vs Pola-R-CHP (POLARIX-like)") + theme_bw()
    ggplotly(p)
  })

  # ── TME & BIOMARKERS ─────────────────────────────────────────────────────
  output$tme_immune <- renderPlotly({
    p <- sim_result() %>%
      select(day, CD8_T, NK_cell) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, color = name)) + geom_line(linewidth = 1) +
      labs(x = "Day", y = "Cells (AU)", color = "Cell type") + theme_bw()
    ggplotly(p)
  })
  output$tme_cart <- renderPlotly({
    p <- sim_result() %>%
      ggplot(aes(day, CART_cells)) + geom_line(color = "#E67E22", linewidth = 1.3) +
      labs(x = "Day", y = "CAR-T Cells (AU)") + theme_bw()
    ggplotly(p)
  })
  output$bm_ldh <- renderPlotly({
    p <- sim_result() %>%
      ggplot(aes(day, LDH_level)) + geom_line(color = "#3498DB", linewidth = 1.2) +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "red") +
      annotate("text", x = 5, y = 1.05, label = "ULN", size = 3, color = "red") +
      labs(x = "Day", y = "LDH (AU)") + theme_bw()
    ggplotly(p)
  })
  output$bm_ctdna <- renderPlotly({
    p <- sim_result() %>%
      ggplot(aes(day, ctDNA_level)) + geom_line(color = "#8E44AD", linewidth = 1.2) +
      labs(x = "Day", y = "ctDNA (AU)") + theme_bw()
    ggplotly(p)
  })
  output$bm_pdl1 <- renderPlotly({
    p <- sim_result() %>%
      mutate(PDL1_proxy = 0.3 + 0.4 * pNFkB + 0.3 * pSTAT3) %>%
      ggplot(aes(day, PDL1_proxy)) + geom_line(color = "#E74C3C", linewidth = 1.1) +
      labs(x = "Day", y = "PD-L1 level (proxy AU)") + theme_bw()
    ggplotly(p)
  })
  output$bm_bcl2ratio <- renderPlotly({
    p <- sim_result() %>%
      mutate(ratio = BCL2_prot / (BIM_prot + 0.001)) %>%
      ggplot(aes(day, ratio)) + geom_line(color = "#C0392B", linewidth = 1.1) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(x = "Day", y = "BCL-2 / BIM ratio") + theme_bw()
    ggplotly(p)
  })
}

shinyApp(ui, server)
