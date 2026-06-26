## =============================================================================
## Colorectal Cancer (CRC) QSP — Interactive Shiny Dashboard
## 6 Tabs: Patient Profile · PK · PD Indicators · Clinical Endpoints
##         Scenario Comparison · Biomarkers
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)

## ─────────────────────────────────────────────────────────────────────────────
## INLINE MODEL (subset of full model for Shiny responsiveness)
## ─────────────────────────────────────────────────────────────────────────────
crc_code_shiny <- '
$PARAM
  CL_FU=112  V1_FU=11  V2_FU=28  Q_FU=30  Kic_FU=0.5  Kout_ic=0.8
  IC50_FU=0.05  Emax_FU=0.80  Hill_FU=1.5
  CL_OX=10.2  V1_OX=4.9  V2_OX=21  Q_OX=4.6
  k_adduct=0.018  k_repair=0.003
  IC50_OX=0.2  Emax_OX=0.70  Hill_OX=1.2
  CL_IRI=20  V1_IRI=15  V2_IRI=120  Q_IRI=14
  k_conv=0.15  k_gluc=0.22  UGT_eff=1.0  CL_SN38=4.0
  IC50_SN38=0.015  Emax_SN38=0.75  Hill_SN38=1.3
  CL_BEV=0.197  V1_BEV=2.92
  Kin_VEGF=0.05  Kout_VEGF=0.04
  Kon_BEV=0.12  Koff_BEV=0.0001  Kint=0.001
  VEGF0=1.25  EC50_BEV=0.4  Emax_BEV=0.6
  CL_CTX=0.384  V1_CTX=3.57  k_EGFR=0.08  Kd_CTX=0.3
  Emax_CTX=0.65  KRAS_mut=0
  CL_PEM=0.214  V1_PEM=5.0  Kd_PEM=0.0004
  MSI_H=0  CD8_base=0.3  CD8_max=0.9  k_CD8act=0.005
  kg_s=0.008  kg_r=0.006  k_kill_s=0.025  k_kill_r=0.004
  k_mutate=1e-6  k_immune=0.003  K_carry=1e10  Ts0=5e8
  kout_CEA=0.004  kin_CEA0=0.004
  kout_ctDNA=0.05

$CMT FU1 FU2 FU_ic OX1 OX_DNA IRI1 SN38 SN38G BEV1 BEV_VEGF CTX1 PEM1
     Ts Tr CEA ctDNA CD8eff VEGF_free EGFR_occ PD1_occ

$INIT FU1=0 FU2=0 FU_ic=0 OX1=0 OX_DNA=0 IRI1=0 SN38=0 SN38G=0
      BEV1=0 BEV_VEGF=0 CTX1=0 PEM1=0
      Ts=5e8 Tr=1e6 CEA=10 ctDNA=0.01 CD8eff=0.3
      VEGF_free=1.25 EGFR_occ=0 PD1_occ=0

$ODE
  double Cp_FU=FU1/V1_FU;
  double FU_eff=Emax_FU*pow(FU_ic/(IC50_FU+FU_ic+1e-9),Hill_FU);
  double OX_eff=Emax_OX*pow(OX_DNA/(IC50_OX+OX_DNA+1e-9),Hill_OX);
  double SN38_eff=Emax_SN38*pow(SN38/(IC50_SN38+SN38+1e-9),Hill_SN38);
  double VEGF_rel=VEGF_free/(EC50_BEV+VEGF_free+1e-9);
  double BEV_TGI=Emax_BEV*(1.0-VEGF_rel);
  double CTX_eff=Emax_CTX*EGFR_occ*(1.0-KRAS_mut);
  double CD8_kill=CD8eff*k_immune*MSI_H+CD8eff*k_immune*0.1*(1.0-MSI_H);
  double kill_s=k_kill_s*(FU_eff+OX_eff+SN38_eff+CTX_eff+BEV_TGI)+CD8_kill;
  if(kill_s>0.99) kill_s=0.99;
  double kill_r=k_kill_r*(FU_eff+OX_eff+SN38_eff)+CD8_kill;
  double Ttot=Ts+Tr;
  double gfrac=1.0-Ttot/K_carry;
  if(gfrac<0) gfrac=0;
  dxdt_FU1=-(CL_FU/V1_FU)*FU1-(Q_FU/V1_FU)*FU1+(Q_FU/V2_FU)*FU2-Kic_FU*FU1;
  dxdt_FU2=(Q_FU/V1_FU)*FU1-(Q_FU/V2_FU)*FU2;
  dxdt_FU_ic=Kic_FU*Cp_FU-Kout_ic*FU_ic;
  double Cp_OX=OX1/V1_OX;
  dxdt_OX1=-(CL_OX/V1_OX)*OX1-(Q_OX/V1_OX)*OX1;
  dxdt_OX_DNA=k_adduct*Cp_OX-k_repair*OX_DNA;
  double Cp_IRI=IRI1/V1_IRI;
  dxdt_IRI1=-(CL_IRI/V1_IRI)*IRI1-(Q_IRI/V1_IRI)*IRI1-k_conv*IRI1;
  dxdt_SN38=k_conv*Cp_IRI*V1_IRI-CL_SN38*SN38-k_gluc*UGT_eff*SN38;
  dxdt_SN38G=k_gluc*UGT_eff*SN38;
  double Cp_BEV=BEV1/V1_BEV;
  double Cp_BEV_nM=Cp_BEV/0.149;
  dxdt_BEV1=-(CL_BEV/V1_BEV)*BEV1-Kon_BEV*Cp_BEV_nM*VEGF_free*V1_BEV+Koff_BEV*BEV_VEGF;
  dxdt_BEV_VEGF=Kon_BEV*Cp_BEV_nM*VEGF_free-(Koff_BEV+Kint)*BEV_VEGF;
  double VEGF_prod=Kin_VEGF*(Ttot/Ts0);
  dxdt_VEGF_free=VEGF_prod-Kout_VEGF*VEGF_free-Kon_BEV*Cp_BEV_nM*VEGF_free+Koff_BEV*BEV_VEGF;
  double Cp_CTX=CTX1/V1_CTX;
  double Cp_CTX_nM=Cp_CTX/0.145;
  dxdt_CTX1=-(CL_CTX/V1_CTX)*CTX1;
  dxdt_EGFR_occ=k_EGFR*(Cp_CTX_nM/(Kd_CTX+Cp_CTX_nM)-EGFR_occ);
  double Cp_PEM=PEM1/V1_PEM;
  double Cp_PEM_nM=Cp_PEM/0.149;
  dxdt_PEM1=-(CL_PEM/V1_PEM)*PEM1;
  dxdt_PD1_occ=0.01*(Cp_PEM_nM/(Kd_PEM+Cp_PEM_nM)-PD1_occ);
  double CD8_target=CD8_base+(CD8_max-CD8_base)*PD1_occ*MSI_H;
  dxdt_CD8eff=k_CD8act*(CD8_target-CD8eff);
  dxdt_Ts=kg_s*Ts*gfrac-kill_s*Ts-k_mutate*Ts;
  dxdt_Tr=kg_r*Tr*gfrac-kill_r*Tr+k_mutate*Ts;
  double Ttot_norm=Ttot/(Ts0+1e6);
  dxdt_CEA=kin_CEA0*Ttot_norm-kout_CEA*CEA;
  dxdt_ctDNA=1e-9*(kill_s*Ts+kill_r*Tr)-kout_ctDNA*ctDNA;

$TABLE
  double TumDiam=35.0*pow((Ts+Tr)/(Ts0+1e6),0.333);
  double PctChange=(TumDiam-35.0)/35.0*100.0;
  double ResistFrac=Tr/(Ts+Tr+1e-6);
  double Cp_FU_out=FU1/V1_FU;
  double Cp_OX_out=OX1/V1_OX;
  double Cp_SN38_out=SN38;
  double Cp_BEV_out=BEV1/V1_BEV;
  double Cp_CTX_out=CTX1/V1_CTX;
  double Cp_PEM_out=PEM1/V1_PEM;

$CAPTURE Cp_FU_out Cp_OX_out Cp_SN38_out Cp_BEV_out Cp_CTX_out Cp_PEM_out
         TumDiam PctChange ResistFrac CEA ctDNA CD8eff VEGF_free EGFR_occ PD1_occ
'

mod_shiny <- mcode("crc_shiny", crc_code_shiny, quiet = TRUE)

## ─────────────────────────────────────────────────────────────────────────────
## HELPER: build dosing events from regimen
## ─────────────────────────────────────────────────────────────────────────────
build_events <- function(reg, bsa, wt, n_cycles, ugt_eff = 1.0) {
  evts <- data.frame(time = numeric(), cmt = character(), amt = numeric(),
                     rate = numeric(), evid = integer(), stringsAsFactors = FALSE)

  if (reg %in% c("FOLFOX6", "FOLFOX_BEV")) {
    for (i in seq_len(n_cycles)) {
      t0 <- (i - 1) * 14 * 24
      evts <- bind_rows(evts,
        data.frame(time = t0, cmt = "OX1", amt = 85*bsa, rate = 85*bsa/2, evid = 1),
        data.frame(time = t0, cmt = "FU1", amt = 400*bsa, rate = -2, evid = 1),
        data.frame(time = t0, cmt = "FU1", amt = 2400*bsa, rate = 2400*bsa/46, evid = 1))
    }
  }
  if (reg %in% c("FOLFIRI", "FOLFIRI_BEV", "FOLFIRI_CTX")) {
    for (i in seq_len(n_cycles)) {
      t0 <- (i - 1) * 14 * 24
      evts <- bind_rows(evts,
        data.frame(time = t0, cmt = "IRI1", amt = 180*bsa, rate = 180*bsa/1.5, evid = 1),
        data.frame(time = t0, cmt = "FU1",  amt = 400*bsa, rate = -2, evid = 1),
        data.frame(time = t0, cmt = "FU1",  amt = 2400*bsa, rate = 2400*bsa/46, evid = 1))
    }
  }
  if (reg %in% c("FOLFOX_BEV", "FOLFIRI_BEV")) {
    for (i in seq_len(n_cycles)) {
      t0 <- (i - 1) * 14 * 24
      evts <- bind_rows(evts,
        data.frame(time = t0, cmt = "BEV1", amt = 5*wt, rate = 5*wt/0.5, evid = 1))
    }
  }
  if (reg == "FOLFIRI_CTX") {
    for (i in seq_len(n_cycles * 2)) {  # weekly
      t0 <- (i - 1) * 7 * 24
      dose <- ifelse(i == 1, 400*bsa, 250*bsa)
      evts <- bind_rows(evts,
        data.frame(time = t0, cmt = "CTX1", amt = dose, rate = dose/2, evid = 1))
    }
  }
  if (reg == "Pembro") {
    for (i in seq_len(ceiling(n_cycles * 14 / 21))) {
      t0 <- (i - 1) * 21 * 24
      evts <- bind_rows(evts,
        data.frame(time = t0, cmt = "PEM1", amt = 200, rate = 200/0.5, evid = 1))
    }
  }
  evts
}

## ─────────────────────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "CRC QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile", tabName = "patient",    icon = icon("user")),
      menuItem("PK Profiles",     tabName = "pk",         icon = icon("pills")),
      menuItem("PD Indicators",   tabName = "pd",         icon = icon("chart-line")),
      menuItem("Clinical Endpoints", tabName = "clinical",icon = icon("hospital")),
      menuItem("Scenario Compare", tabName = "scenario",  icon = icon("layer-group")),
      menuItem("Biomarkers",      tabName = "biomarker",  icon = icon("dna"))
    )
  ),

  dashboardBody(
    tabItems(

      ## ── TAB 1: PATIENT PROFILE ──────────────────────────────────────────
      tabItem("patient",
        fluidRow(
          box(title = "Patient Characteristics", width = 4, status = "primary",
            numericInput("age",  "Age (years)", 62, 25, 85),
            numericInput("bsa",  "BSA (m²)", 1.73, 1.2, 2.5, step = 0.01),
            numericInput("wt",   "Body Weight (kg)", 70, 40, 130),
            selectInput("ecog",  "ECOG PS", c("0","1","2")),
            selectInput("stage", "Disease Stage",
              c("Stage III (adjuvant)" = "III",
                "Stage IV metastatic"  = "IV"))
          ),
          box(title = "Tumor Molecular Profile", width = 4, status = "warning",
            selectInput("kras",  "RAS Status",
              c("RAS Wild-Type" = "0", "RAS Mutant (KRAS G12x/G13x)" = "1")),
            selectInput("braf",  "BRAF Status",
              c("BRAF Wild-Type" = "0", "BRAF V600E Mutant" = "1")),
            selectInput("msi",   "MSI / MMR Status",
              c("MSS / pMMR" = "0", "MSI-H / dMMR" = "1")),
            selectInput("pik3ca","PIK3CA Status",
              c("Wild-Type" = "0", "PIK3CA Mutant" = "1")),
            selectInput("cms",   "CMS Subtype",
              c("CMS1 (MSI immune)" = "1", "CMS2 (canonical)" = "2",
                "CMS3 (metabolic)"  = "3", "CMS4 (mesenchymal)" = "4"))
          ),
          box(title = "Treatment Regimen", width = 4, status = "success",
            selectInput("regimen", "Select Regimen",
              c("No Treatment"         = "None",
                "FOLFOX6"              = "FOLFOX6",
                "FOLFIRI"              = "FOLFIRI",
                "FOLFOX + Bevacizumab" = "FOLFOX_BEV",
                "FOLFIRI + Bevacizumab"= "FOLFIRI_BEV",
                "FOLFIRI + Cetuximab"  = "FOLFIRI_CTX",
                "Pembrolizumab (MSI-H)"= "Pembro"),
              selected = "FOLFOX6"),
            numericInput("n_cycles", "Number of Cycles (Q2W)", 12, 1, 24),
            selectInput("ugt_poly", "UGT1A1 Polymorphism",
              c("*1/*1 Normal"      = "1.0",
                "*1/*28 Hetero"     = "0.75",
                "*28/*28 Poor glucuronidator" = "0.5")),
            actionButton("run_sim", "Run Simulation", class = "btn-success btn-lg",
                         icon = icon("play"))
          )
        ),
        fluidRow(
          box(title = "Biomarker Eligibility Guide", width = 12, status = "info",
            tableOutput("eligibility_tbl"))
        )
      ),

      ## ── TAB 2: PK PROFILES ──────────────────────────────────────────────
      tabItem("pk",
        fluidRow(
          box(title = "5-FU Plasma Concentration", width = 6, plotlyOutput("pk_fu")),
          box(title = "Oxaliplatin / SN-38 Concentration", width = 6, plotlyOutput("pk_ox_sn38"))
        ),
        fluidRow(
          box(title = "Bevacizumab Plasma & Free VEGF", width = 6, plotlyOutput("pk_bev")),
          box(title = "Cetuximab / EGFR Occupancy", width = 6, plotlyOutput("pk_ctx"))
        ),
        fluidRow(
          box(title = "Pembrolizumab Plasma & PD-1 Occupancy", width = 6, plotlyOutput("pk_pem")),
          box(title = "PK Summary Table (Cycle 1)", width = 6, DTOutput("pk_summary"))
        )
      ),

      ## ── TAB 3: PD INDICATORS ────────────────────────────────────────────
      tabItem("pd",
        fluidRow(
          box(title = "Tumor Diameter Change (RECIST)", width = 6, plotlyOutput("pd_recist")),
          box(title = "Tumor Cell Kinetics (Sensitive vs Resistant)", width = 6, plotlyOutput("pd_cells"))
        ),
        fluidRow(
          box(title = "CD8 T Cell Activity", width = 6, plotlyOutput("pd_cd8")),
          box(title = "Resistant Cell Fraction", width = 6, plotlyOutput("pd_resist"))
        )
      ),

      ## ── TAB 4: CLINICAL ENDPOINTS ───────────────────────────────────────
      tabItem("clinical",
        fluidRow(
          valueBoxOutput("vbox_bor"),
          valueBoxOutput("vbox_pfs"),
          valueBoxOutput("vbox_orr")
        ),
        fluidRow(
          box(title = "Tumor Response Over Time", width = 8, plotlyOutput("clin_response")),
          box(title = "RECIST Classification", width = 4,
            tableOutput("recist_tbl"),
            tags$hr(),
            tags$b("RECIST 1.1 Thresholds:"),
            tags$ul(
              tags$li("CR: Disappearance of all lesions"),
              tags$li("PR: ≥30% decrease from baseline SLD"),
              tags$li("PD: ≥20% increase from nadir"),
              tags$li("SD: Neither PR nor PD criteria")
            )
          )
        )
      ),

      ## ── TAB 5: SCENARIO COMPARISON ──────────────────────────────────────
      tabItem("scenario",
        fluidRow(
          box(title = "Select Scenarios to Compare", width = 12,
            checkboxGroupInput("compare_regs", "",
              choices = c("No Treatment" = "None",
                          "FOLFOX6" = "FOLFOX6",
                          "FOLFIRI"  = "FOLFIRI",
                          "FOLFOX + Bevacizumab" = "FOLFOX_BEV",
                          "FOLFIRI + Bevacizumab" = "FOLFIRI_BEV",
                          "FOLFIRI + Cetuximab (RAS-WT)" = "FOLFIRI_CTX",
                          "Pembrolizumab (MSI-H)" = "Pembro"),
              selected = c("FOLFOX6", "FOLFIRI_BEV", "Pembro"),
              inline = TRUE),
            actionButton("run_compare", "Compare Selected", class = "btn-primary",
                         icon = icon("balance-scale"))
          )
        ),
        fluidRow(
          box(title = "Tumor Diameter Comparison", width = 12, plotlyOutput("cmp_diam"))
        ),
        fluidRow(
          box(title = "CEA Kinetics Comparison", width = 6, plotlyOutput("cmp_cea")),
          box(title = "Summary Statistics", width = 6, DTOutput("cmp_summary"))
        )
      ),

      ## ── TAB 6: BIOMARKERS ───────────────────────────────────────────────
      tabItem("biomarker",
        fluidRow(
          box(title = "CEA Kinetics", width = 6, plotlyOutput("bm_cea")),
          box(title = "ctDNA (Liquid Biopsy)", width = 6, plotlyOutput("bm_ctdna"))
        ),
        fluidRow(
          box(title = "Free VEGF Suppression", width = 6, plotlyOutput("bm_vegf")),
          box(title = "PD-1 / EGFR Occupancy", width = 6, plotlyOutput("bm_occ"))
        ),
        fluidRow(
          box(title = "Biomarker Reference Ranges & Interpretation", width = 12,
            tableOutput("bm_ref_tbl"))
        )
      )
    )
  )
)

## ─────────────────────────────────────────────────────────────────────────────
## SERVER
## ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## ── Run main simulation ────────────────────────────────────────────────
  sim_result <- eventReactive(input$run_sim, {
    req(input$regimen, input$bsa, input$wt, input$n_cycles)
    bsa       <- input$bsa
    wt        <- input$wt
    n_cycles  <- input$n_cycles
    kras_val  <- as.numeric(input$kras)
    msi_val   <- as.numeric(input$msi)
    ugt_val   <- as.numeric(input$ugt_poly)

    if (input$regimen == "Pembro" && msi_val == 0) {
      showNotification("Pembrolizumab is primarily indicated for MSI-H tumors. MSI status set to H.",
                       type = "warning", duration = 5)
    }

    evts <- build_events(input$regimen, bsa, wt, n_cycles, ugt_val)
    sim_end <- max(n_cycles * 14 * 24 + 2000, 4000)

    if (nrow(evts) == 0) {
      # No treatment — just run blank
      evts <- data.frame(time = 0, cmt = "FU1", amt = 0, evid = 0, rate = 0)
    }

    mod_run <- mod_shiny %>%
      param(KRAS_mut = kras_val,
            MSI_H    = msi_val,
            UGT_eff  = ugt_val)

    out <- mod_run %>%
      ev(as_data_set(evts)) %>%
      mrgsim(end = sim_end, delta = 12, obsonly = TRUE) %>%
      as.data.frame() %>%
      mutate(time_days   = time / 24,
             time_months = time_days / 30.4)
    out
  })

  ## ── Comparison simulation ───────────────────────────────────────────────
  compare_result <- eventReactive(input$run_compare, {
    req(input$compare_regs, input$bsa, input$wt, input$n_cycles)
    bsa      <- input$bsa
    wt       <- input$wt
    n_cycles <- input$n_cycles
    kras_val <- as.numeric(input$kras)
    msi_val  <- as.numeric(input$msi)
    ugt_val  <- as.numeric(input$ugt_poly)
    all_out  <- data.frame()

    for (reg in input$compare_regs) {
      evts <- build_events(reg, bsa, wt, n_cycles, ugt_val)
      if (nrow(evts) == 0)
        evts <- data.frame(time = 0, cmt = "FU1", amt = 0, evid = 0, rate = 0)
      msi_use <- ifelse(reg == "Pembro", max(msi_val, 0), msi_val)
      out <- mod_shiny %>%
        param(KRAS_mut = ifelse(reg == "FOLFIRI_CTX", kras_val, kras_val),
              MSI_H    = msi_use,
              UGT_eff  = ugt_val) %>%
        ev(as_data_set(evts)) %>%
        mrgsim(end = max(n_cycles * 14 * 24 + 2000, 4000), delta = 24, obsonly = TRUE) %>%
        as.data.frame() %>%
        mutate(scenario = reg, time_months = time / 24 / 30.4)
      all_out <- bind_rows(all_out, out)
    }
    all_out
  })

  ## ── TAB 1: Eligibility Table ────────────────────────────────────────────
  output$eligibility_tbl <- renderTable({
    data.frame(
      Biomarker = c("RAS (KRAS/NRAS)", "BRAF V600E", "MSI / dMMR", "PIK3CA", "HER2"),
      Status    = c(c("WT","Mutant")[as.integer(input$kras) + 1],
                    c("WT","V600E")[as.integer(input$braf) + 1],
                    c("MSS","MSI-H")[as.integer(input$msi) + 1],
                    c("WT","Mutant")[as.integer(input$pik3ca) + 1],
                    "Not tested"),
      Implication = c(
        c("Eligible: anti-EGFR (cetuximab/panitumumab)",
          "RESIST anti-EGFR; eligible KRAS G12C → sotorasib")[as.integer(input$kras)+1],
        c("Standard chemo eligible",
          "BRAF V600E: FOLFOXIRI + bevacizumab or encorafenib + cetuximab")[as.integer(input$braf)+1],
        c("Standard chemo; immunotherapy not preferred",
          "FIRST-LINE pembrolizumab (KEYNOTE-177); nivolumab + ipilimumab")[as.integer(input$msi)+1],
        c("No specific therapy", "Possible PI3K inhibitor consideration")[as.integer(input$pik3ca)+1],
        "Trastuzumab combinations (HER2-amp CRC)"
      ),
      stringsAsFactors = FALSE
    )
  })

  ## ── TAB 2: PK plots ─────────────────────────────────────────────────────
  output$pk_fu <- renderPlotly({
    req(sim_result())
    d <- sim_result() %>% filter(time_months <= min(2, max(time_months)))
    p <- ggplot(d, aes(x = time * 1, y = Cp_FU_out)) +
      geom_line(color = "#2980B9", linewidth = 1) +
      labs(title = "5-FU Plasma Conc.", x = "Time (h)", y = "5-FU (μg/mL)") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_ox_sn38 <- renderPlotly({
    req(sim_result())
    d <- sim_result() %>%
      filter(time_months <= 2) %>%
      select(time, Cp_OX_out, Cp_SN38_out) %>%
      pivot_longer(-time, names_to = "drug", values_to = "Cp")
    p <- ggplot(d, aes(x = time, y = Cp, color = drug)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c("#A9CCE3","#7DCEA0"),
                         labels = c("Oxaliplatin (μg/mL)","SN-38 (μg/mL)")) +
      labs(title = "Oxaliplatin / SN-38 Concentration",
           x = "Time (h)", y = "Concentration (μg/mL)", color = "") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_bev <- renderPlotly({
    req(sim_result())
    d <- sim_result() %>% select(time_months, Cp_BEV_out, VEGF_free) %>%
      pivot_longer(-time_months, names_to = "var", values_to = "val")
    p <- ggplot(d, aes(x = time_months, y = val, color = var)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c("#F1948A","#FADBD8"),
                         labels = c("Bevacizumab (mg/mL)","Free VEGF (nM)")) +
      labs(title = "Bevacizumab & Free VEGF",
           x = "Time (months)", y = "Concentration", color = "") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_ctx <- renderPlotly({
    req(sim_result())
    d <- sim_result() %>% select(time_months, Cp_CTX_out, EGFR_occ)
    p <- ggplot(d) +
      geom_line(aes(x = time_months, y = Cp_CTX_out), color = "#D7BDE2", linewidth = 1) +
      geom_line(aes(x = time_months, y = EGFR_occ * 100), color = "#8E44AD",
                linewidth = 1, linetype = "dashed") +
      scale_y_continuous(name = "Cetuximab (mg/mL)",
                         sec.axis = sec_axis(~./100, name = "EGFR Occupancy (fraction)")) +
      labs(title = "Cetuximab PK & EGFR Occupancy", x = "Time (months)") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_pem <- renderPlotly({
    req(sim_result())
    d <- sim_result() %>% select(time_months, Cp_PEM_out, PD1_occ)
    p <- ggplot(d) +
      geom_line(aes(x = time_months, y = Cp_PEM_out), color = "#A9DFBF", linewidth = 1) +
      geom_line(aes(x = time_months, y = PD1_occ * 100), color = "#1ABC9C",
                linewidth = 1, linetype = "dashed") +
      scale_y_continuous(name = "Pembrolizumab (mg/mL)",
                         sec.axis = sec_axis(~./100, name = "PD-1 Occupancy (fraction)")) +
      labs(title = "Pembrolizumab PK & PD-1 Occupancy", x = "Time (months)") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_summary <- renderDT({
    req(sim_result())
    sim_result() %>%
      filter(time_months <= 24/30.4) %>%
      summarise(
        `5-FU Cmax (μg/mL)`  = max(Cp_FU_out),
        `Oxaliplatin Cmax`    = max(Cp_OX_out),
        `SN-38 Cmax (μg/mL)` = max(Cp_SN38_out),
        `Bevacizumab Cmax`    = max(Cp_BEV_out),
        `Cetuximab Cmax`      = max(Cp_CTX_out),
        `Pembrolizumab Cmax`  = max(Cp_PEM_out)
      ) %>%
      round(4) %>%
      datatable(options = list(dom = 't'), rownames = FALSE)
  })

  ## ── TAB 3: PD plots ─────────────────────────────────────────────────────
  output$pd_recist <- renderPlotly({
    req(sim_result())
    p <- ggplot(sim_result(), aes(x = time_months, y = PctChange)) +
      geom_line(color = "#2980B9", linewidth = 1.2) +
      geom_hline(yintercept = 20,  linetype = "dashed", color = "red") +
      geom_hline(yintercept = -30, linetype = "dashed", color = "darkgreen") +
      coord_cartesian(ylim = c(-100, max(200, max(sim_result()$PctChange) * 1.1))) +
      labs(title = "Tumor Diameter Change (RECIST 1.1)",
           x = "Time (months)", y = "Change from Baseline (%)") +
      theme_bw()
    ggplotly(p)
  })

  output$pd_cells <- renderPlotly({
    req(sim_result())
    d <- sim_result() %>%
      select(time_months, Ts, Tr) %>%
      pivot_longer(-time_months, names_to = "type", values_to = "cells")
    p <- ggplot(d, aes(x = time_months, y = cells / 1e8, color = type)) +
      geom_area(aes(fill = type), alpha = 0.4, position = "stack") +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c("#2980B9","#E74C3C"),
                         labels = c("Sensitive","Resistant")) +
      scale_fill_manual(values = c("#2980B9","#E74C3C"),
                        labels = c("Sensitive","Resistant")) +
      labs(title = "Tumor Cell Kinetics",
           x = "Time (months)", y = "Cells (×10⁸)", color = "", fill = "") +
      theme_bw()
    ggplotly(p)
  })

  output$pd_cd8 <- renderPlotly({
    req(sim_result())
    p <- ggplot(sim_result(), aes(x = time_months, y = CD8eff)) +
      geom_line(color = "#1ABC9C", linewidth = 1.3) +
      geom_hline(yintercept = 0.3, linetype = "dotted", color = "gray50") +
      coord_cartesian(ylim = c(0, 1)) +
      labs(title = "CD8 T Cell Effector Activity",
           x = "Time (months)", y = "Activity (0–1 scale)") +
      theme_bw()
    ggplotly(p)
  })

  output$pd_resist <- renderPlotly({
    req(sim_result())
    p <- ggplot(sim_result(), aes(x = time_months, y = ResistFrac * 100)) +
      geom_line(color = "#E74C3C", linewidth = 1.2) +
      labs(title = "Resistant Cell Fraction Over Time",
           x = "Time (months)", y = "Resistant Fraction (%)") +
      theme_bw()
    ggplotly(p)
  })

  ## ── TAB 4: Clinical Endpoints ────────────────────────────────────────────
  best_response <- reactive({
    req(sim_result())
    min(sim_result()$PctChange)
  })

  pfs_days <- reactive({
    req(sim_result())
    d <- sim_result()
    pd_t <- d$time_days[d$PctChange >= 20]
    if (length(pd_t) > 0) min(pd_t) else max(d$time_days)
  })

  output$vbox_bor <- renderValueBox({
    br <- round(best_response(), 1)
    resp_cat <- dplyr::case_when(
      br <= -100 ~ "CR", br <= -30 ~ "PR",
      br < 20    ~ "SD", TRUE       ~ "PD"
    )
    valueBox(paste0(br, "%"), paste("Best Response:", resp_cat),
             icon = icon("ruler"), color = "blue")
  })

  output$vbox_pfs <- renderValueBox({
    pfs_m <- round(pfs_days() / 30.4, 1)
    valueBox(paste0(pfs_m, " mo"), "Simulated PFS",
             icon = icon("calendar"), color = "green")
  })

  output$vbox_orr <- renderValueBox({
    br <- best_response()
    orr <- ifelse(br <= -30, "Responder (PR/CR)", "Non-responder (SD/PD)")
    valueBox(orr, "Response Classification",
             icon = icon("check-circle"), color = ifelse(br <= -30, "green","red"))
  })

  output$clin_response <- renderPlotly({
    req(sim_result())
    p <- ggplot(sim_result(), aes(x = time_months, y = PctChange)) +
      geom_ribbon(aes(ymin = -30, ymax = 20), fill = "lightblue", alpha = 0.2) +
      geom_line(color = "#2980B9", linewidth = 1.3) +
      geom_hline(yintercept = 20,  linetype = "dashed", color = "red",   alpha = 0.8) +
      geom_hline(yintercept = -30, linetype = "dashed", color = "green4",alpha = 0.8) +
      annotate("text", x = 0.5, y = 22,  label = "PD threshold (+20%)",
               hjust = 0, color = "red",    size = 3.5) +
      annotate("text", x = 0.5, y = -32, label = "PR threshold (-30%)",
               hjust = 0, color = "green4", size = 3.5) +
      coord_cartesian(ylim = c(-100, min(300, max(sim_result()$PctChange) * 1.2))) +
      labs(title = "RECIST Tumor Response",
           x = "Time (months)", y = "% Change from Baseline SLD") +
      theme_bw(base_size = 14)
    ggplotly(p)
  })

  output$recist_tbl <- renderTable({
    data.frame(
      Response = c("CR", "PR", "SD", "PD"),
      Threshold = c("Disappearance", "≥30% decrease", "Neither PR/PD", "≥20% increase"),
      `Simulated` = c(
        ifelse(best_response() <= -100, "✓",""),
        ifelse(best_response() <= -30 & best_response() > -100, "✓",""),
        ifelse(best_response() > -30  & best_response() < 20, "✓",""),
        ifelse(best_response() >= 20, "✓","")
      ), stringsAsFactors = FALSE
    )
  })

  ## ── TAB 5: Scenario Comparison ───────────────────────────────────────────
  output$cmp_diam <- renderPlotly({
    req(compare_result())
    p <- ggplot(compare_result(), aes(x = time_months, y = PctChange, color = scenario)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = c(20, -30), linetype = "dashed",
                 color = c("red","green4"), alpha = 0.7) +
      coord_cartesian(ylim = c(-100, 300)) +
      labs(title = "Tumor Diameter — Multi-Regimen Comparison",
           x = "Time (months)", y = "% Change from Baseline", color = "Regimen") +
      theme_bw(base_size = 13)
    ggplotly(p)
  })

  output$cmp_cea <- renderPlotly({
    req(compare_result())
    p <- ggplot(compare_result() %>% filter(time_months <= 18),
                aes(x = time_months, y = CEA, color = scenario)) +
      geom_line(linewidth = 1.1) +
      scale_y_log10() +
      labs(title = "CEA Kinetics Comparison", x = "Time (months)",
           y = "CEA (ng/mL, log)", color = "Regimen") +
      theme_bw()
    ggplotly(p)
  })

  output$cmp_summary <- renderDT({
    req(compare_result())
    compare_result() %>%
      group_by(scenario) %>%
      summarise(
        `Best Response (%)` = round(min(PctChange), 1),
        `PFS est. (mo)`     = round(
          suppressWarnings(min(time_months[PctChange >= 20], na.rm = TRUE)), 1),
        `CEA Nadir (ng/mL)` = round(min(CEA), 2),
        `Final Resist. Frac.` = round(last(ResistFrac) * 100, 2),
        .groups = "drop"
      ) %>%
      datatable(rownames = FALSE, options = list(pageLength = 10))
  })

  ## ── TAB 6: Biomarkers ───────────────────────────────────────────────────
  output$bm_cea <- renderPlotly({
    req(sim_result())
    p <- ggplot(sim_result(), aes(x = time_months, y = CEA)) +
      geom_line(color = "#E67E22", linewidth = 1.2) +
      geom_hline(yintercept = 5, linetype = "dashed", color = "gray40") +
      annotate("text", x = 0.5, y = 6, label = "ULN = 5 ng/mL",
               hjust = 0, color = "gray40", size = 3.5) +
      labs(title = "CEA Biomarker Kinetics",
           x = "Time (months)", y = "CEA (ng/mL)") +
      theme_bw()
    ggplotly(p)
  })

  output$bm_ctdna <- renderPlotly({
    req(sim_result())
    p <- ggplot(sim_result(), aes(x = time_months, y = ctDNA)) +
      geom_line(color = "#8E44AD", linewidth = 1.2) +
      scale_y_log10() +
      labs(title = "ctDNA Kinetics (Liquid Biopsy)",
           x = "Time (months)", y = "ctDNA (relative units, log)") +
      theme_bw()
    ggplotly(p)
  })

  output$bm_vegf <- renderPlotly({
    req(sim_result())
    p <- ggplot(sim_result(), aes(x = time_months, y = VEGF_free)) +
      geom_line(color = "#D68910", linewidth = 1.2) +
      geom_hline(yintercept = 1.25, linetype = "dotted", color = "gray50") +
      annotate("text", x = 0.5, y = 1.35, label = "Baseline VEGF",
               hjust = 0, size = 3.5, color = "gray50") +
      labs(title = "Free VEGF-A Level",
           x = "Time (months)", y = "VEGF (nM)") +
      theme_bw()
    ggplotly(p)
  })

  output$bm_occ <- renderPlotly({
    req(sim_result())
    d <- sim_result() %>%
      select(time_months, EGFR_occ, PD1_occ) %>%
      pivot_longer(-time_months, names_to = "target", values_to = "occupancy")
    p <- ggplot(d, aes(x = time_months, y = occupancy * 100, color = target)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("#D7BDE2","#A9DFBF"),
                         labels = c("EGFR Occupancy (%)","PD-1 Occupancy (%)")) +
      coord_cartesian(ylim = c(0, 100)) +
      labs(title = "Target Receptor Occupancy",
           x = "Time (months)", y = "Occupancy (%)", color = "") +
      theme_bw()
    ggplotly(p)
  })

  output$bm_ref_tbl <- renderTable({
    data.frame(
      Biomarker = c("CEA","CA 19-9","ctDNA","VEGF-A","CD8 T cells"),
      `Normal Range`  = c("<5 ng/mL","<37 U/mL","Undetectable","50-200 pg/mL","Variable"),
      `Clinical Use`  = c(
        "Response monitoring, recurrence detection",
        "Pancreatic/biliary involvement",
        "MRD detection, resistance monitoring",
        "Bevacizumab efficacy prediction",
        "Immunotherapy response (MSI-H)"
      ),
      `Key Note` = c(
        "Rise >25% suggests PD",
        "Less CRC-specific than CEA",
        "ctDNA clearance = favorable outcome",
        "High baseline = worse anti-VEGF benefit",
        "High CD8/Treg ratio = better IO response"
      ),
      stringsAsFactors = FALSE
    )
  })
}

## ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
