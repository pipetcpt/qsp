## =============================================================================
## Small Cell Lung Cancer (SCLC) — QSP Shiny Dashboard
## 6 Tabs: Patient Profile | Drug PK | Molecular Pathways |
##         Clinical Endpoints | Scenario Comparison | Biomarkers & TME
## =============================================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)

# ---- Helper: Simplified ODE solver (Euler, hourly) ----
simulate_sclc <- function(params, t_max_days = 365, dt_h = 6) {
  p <- params
  t_seq <- seq(0, t_max_days * 24, by = dt_h)
  n <- length(t_seq)

  # State vector
  sv <- list(
    Atezo_C=0,Atezo_P=0,Durva_C=0,Durva_P=0,
    Carbo_C=0,Etop_C=0,Etop_T=0,Lurbi_C=0,Topo_C=0,
    Tarlat_C=0,Tarlat_T=0,Tarlat_P=0,Trila_C=0,
    DNA=0,Ts=p$Tum_S0,Ta=0,Tr=0.002,
    CD8T=1.0,PDL1_occ=0,NK=1.0,IFNg=0.2
  )

  out <- data.frame(time_h=numeric(n), time_d=numeric(n),
                    Tumor_total=numeric(n), SLD=numeric(n),
                    DNA=numeric(n), CD8T=numeric(n),
                    PDL1_occ=numeric(n), Atezo_C=numeric(n),
                    Durva_C=numeric(n), Tarlat_C=numeric(n),
                    Trila_C=numeric(n), IFNg=numeric(n), NK=numeric(n))

  # Dosing schedule builder
  dose_events <- list()

  # CE dosing
  if (p$f_Carbo > 0 || p$f_Etop > 0) {
    for (cy in 0:(p$n_CE_cycles - 1)) {
      t0 <- cy * p$CE_tau_d * 24
      dose_events <- c(dose_events,
        list(list(t=t0,   cmt="Carbo_C", amt=p$Carbo_dose * p$f_Carbo)),
        list(list(t=t0,   cmt="Etop_C",  amt=p$Etop_dose  * p$f_Etop)),
        list(list(t=t0+24,cmt="Etop_C",  amt=p$Etop_dose  * p$f_Etop)),
        list(list(t=t0+48,cmt="Etop_C",  amt=p$Etop_dose  * p$f_Etop))
      )
    }
  }
  if (p$f_Trila > 0) {
    for (cy in 0:(p$n_CE_cycles - 1)) {
      t0 <- cy * p$CE_tau_d * 24
      dose_events <- c(dose_events,
        list(list(t=t0,   cmt="Trila_C", amt=p$Trila_dose * p$f_Trila)),
        list(list(t=t0+24,cmt="Trila_C", amt=p$Trila_dose * p$f_Trila)),
        list(list(t=t0+48,cmt="Trila_C", amt=p$Trila_dose * p$f_Trila))
      )
    }
  }
  if (p$f_Atezo > 0) {
    for (cy in 0:(p$n_IO_cycles - 1)) {
      t0 <- cy * 21 * 24
      dose_events <- c(dose_events, list(list(t=t0, cmt="Atezo_C", amt=1200*p$f_Atezo)))
    }
  }
  if (p$f_Durva > 0) {
    for (cy in 0:(p$n_IO_cycles - 1)) {
      t0 <- cy * 28 * 24
      dose_events <- c(dose_events, list(list(t=t0, cmt="Durva_C", amt=1500*p$f_Durva)))
    }
  }
  if (p$f_Tarlat > 0) {
    for (cy in 0:(p$n_Tarlat_cycles - 1)) {
      t0 <- cy * 14 * 24
      dose_events <- c(dose_events, list(list(t=t0, cmt="Tarlat_C", amt=p$Tarlat_dose*p$f_Tarlat)))
    }
  }
  if (p$f_Lurbi > 0) {
    for (cy in 0:(p$n_Lurbi_cycles - 1)) {
      t0 <- cy * 21 * 24
      dose_events <- c(dose_events, list(list(t=t0, cmt="Lurbi_C", amt=p$Lurbi_dose*p$f_Lurbi)))
    }
  }
  if (p$f_Topo > 0) {
    for (cy in 0:(p$n_Topo_cycles - 1)) {
      t0 <- cy * 21 * 24
      for (d in 0:4) dose_events <- c(dose_events,
        list(list(t=t0+d*24, cmt="Topo_C", amt=p$Topo_dose*p$f_Topo)))
    }
  }

  # PK rates
  CL_atezo  <- p$Atezo_CL  / (p$Atezo_V1 * 24)
  Q_atezo   <- p$Atezo_Q   / (p$Atezo_V1 * 24)
  R_atezo   <- p$Atezo_V1  / p$Atezo_V2
  CL_durva  <- p$Durva_CL  / (p$Durva_V1 * 24)
  Q_durva   <- p$Durva_Q   / (p$Durva_V1 * 24)
  R_durva   <- p$Durva_V1  / p$Durva_V2
  CL_carbo  <- p$Carbo_CL  / p$Carbo_V
  CL_etop   <- p$Etop_CL   / p$Etop_V1
  Q_etop    <- p$Etop_Q    / p$Etop_V1
  R_etop    <- p$Etop_V1   / p$Etop_V2
  CL_lurbi  <- p$Lurbi_CL  / p$Lurbi_V
  CL_topo   <- p$Topo_CL   / p$Topo_V
  CL_tarlat <- p$Tarlat_CL / (p$Tarlat_V1 * 24)
  Q_tarl1   <- p$Tarlat_Q1 / (p$Tarlat_V1 * 24)
  Q_tarl2   <- p$Tarlat_Q2 / (p$Tarlat_V1 * 24)
  R_tarl1   <- p$Tarlat_V1 / p$Tarlat_Vt
  R_tarl2   <- p$Tarlat_V1 / p$Tarlat_V2
  CL_trila  <- p$Trila_CL  / p$Trila_V

  for (i in seq_along(t_seq)) {
    t_h <- t_seq[i]

    # Apply doses at this time step
    for (ev in dose_events) {
      if (abs(ev$t - t_h) < dt_h / 2) {
        sv[[ev$cmt]] <- sv[[ev$cmt]] + ev$amt / switch(ev$cmt,
          Atezo_C=p$Atezo_V1, Durva_C=p$Durva_V1, Carbo_C=p$Carbo_V,
          Etop_C=p$Etop_V1, Lurbi_C=p$Lurbi_V, Topo_C=p$Topo_V,
          Tarlat_C=p$Tarlat_V1, Trila_C=p$Trila_V, 1.0)
      }
    }

    Ttot <- max(sv$Ts + sv$Ta + sv$Tr, 1e-6)
    Tnorm <- Ttot / p$Tum_S0

    Atezo_nM <- sv$Atezo_C * 1000 / 145000
    Durva_nM <- sv$Durva_C * 1000 / 148000
    PDL1_RO  <- (Atezo_nM + Durva_nM) / (p$Kd_PDL1 + Atezo_nM + Durva_nM)
    Tarlat_nM <- sv$Tarlat_T * 1000 / 146000

    E_DNA   <- p$DNA_kill  * sv$DNA / (p$EC50_DNA + sv$DNA)
    E_Tarl  <- p$f_Tarlat  * p$Tarlat_Emax * p$DLL3_expr * sv$CD8T *
               Tarlat_nM / (p$Tarlat_EC50 + Tarlat_nM + 1e-9)
    E_CD8k  <- p$kCD8_kill  * sv$CD8T  * PDL1_RO
    E_NK    <- p$kNK_act    * sv$NK

    gom     <- p$kTumor_g * log(pmax(p$Tum_cap / Ttot, 1.0001))
    Trila_E <- p$f_Trila * p$Trila_Emax * sv$Trila_C / (p$Trila_EC50 + sv$Trila_C + 1e-9)

    E_IO_CD8  <- (p$f_Atezo * sv$Atezo_C/(sv$Atezo_C+50) +
                  p$f_Durva * sv$Durva_C/(sv$Durva_C+50)) * PDL1_RO
    E_BITE_CD8 <- p$f_Tarlat * p$kBiTE_kill * Tarlat_nM /
                  (p$Tarlat_EC50 + Tarlat_nM + 1e-9) * p$DLL3_expr

    # ODE RHS
    dAtezoC <- -CL_atezo * sv$Atezo_C - Q_atezo * (sv$Atezo_C - sv$Atezo_P * R_atezo)
    dAtezoP <- Q_atezo * (sv$Atezo_C - sv$Atezo_P * R_atezo)
    dDurvaC <- -CL_durva * sv$Durva_C - Q_durva * (sv$Durva_C - sv$Durva_P * R_durva)
    dDurvaP <- Q_durva * (sv$Durva_C - sv$Durva_P * R_durva)
    dCarboC <- -CL_carbo * sv$Carbo_C
    dEtopC  <- -CL_etop * sv$Etop_C - Q_etop * (sv$Etop_C - sv$Etop_T * R_etop)
    dEtopT  <- Q_etop * (sv$Etop_C - sv$Etop_T * R_etop)
    dLurbiC <- -CL_lurbi * sv$Lurbi_C
    dTopoC  <- -CL_topo  * sv$Topo_C
    dTarlatC<- -(CL_tarlat+Q_tarl1+Q_tarl2)*sv$Tarlat_C +
               Q_tarl1*sv$Tarlat_T*R_tarl1 + Q_tarl2*sv$Tarlat_P*R_tarl2
    dTarlatT<- Q_tarl1*(sv$Tarlat_C - sv$Tarlat_T*R_tarl1) - E_Tarl*sv$Tarlat_T
    dTarlatP<- Q_tarl2*(sv$Tarlat_C - sv$Tarlat_P*R_tarl2)
    dTrilaC <- -CL_trila * sv$Trila_C

    DNA_in  <- (p$kDNA_form*(p$f_Etop*sv$Etop_T + p$f_Carbo*sv$Carbo_C*0.5 +
                              p$f_Lurbi*sv$Lurbi_C*0.8 + p$f_Topo*sv$Topo_C*0.6))
    dDNA    <- DNA_in - p$kDNA_rep * sv$DNA
    dTs     <- gom*sv$Ts - E_DNA*sv$Ts - E_CD8k*sv$Ts - E_NK*sv$Ts*0.5 -
               E_Tarl*sv$Ts - Trila_E*gom*sv$Ts - p$kResist*sv$Ts
    dTa     <- Trila_E*gom*sv$Ts - E_DNA*sv$Ta*0.5 - E_CD8k*sv$Ta - p$kTumor_g*sv$Ta
    dTr     <- p$kResist*sv$Ts + p$kTumor_g*0.8*log(pmax(p$Tum_cap/(sv$Tr+1e-6),1.0001))*sv$Tr -
               E_CD8k*sv$Tr*0.3 - E_Tarl*sv$Tr*0.4
    dCD8T   <- p$kCD8_stim*E_IO_CD8*(1-sv$PDL1_occ) + E_BITE_CD8 -
               p$kCD8_death*sv$CD8T - p$kCD8_kill*sv$CD8T*Tnorm*0.5
    dPDL1   <- p$k_PDL1_on*(PDL1_RO - sv$PDL1_occ)
    dNK     <- 0.005*sv$IFNg*(2-sv$NK) - 0.005*sv$NK*Tnorm
    dIFNg   <- p$kIFNg_prod*sv$CD8T*PDL1_RO - p$kIFNg_deg*sv$IFNg

    sv$Atezo_C  <- max(0, sv$Atezo_C  + dAtezoC  * dt_h)
    sv$Atezo_P  <- max(0, sv$Atezo_P  + dAtezoP  * dt_h)
    sv$Durva_C  <- max(0, sv$Durva_C  + dDurvaC  * dt_h)
    sv$Durva_P  <- max(0, sv$Durva_P  + dDurvaP  * dt_h)
    sv$Carbo_C  <- max(0, sv$Carbo_C  + dCarboC  * dt_h)
    sv$Etop_C   <- max(0, sv$Etop_C   + dEtopC   * dt_h)
    sv$Etop_T   <- max(0, sv$Etop_T   + dEtopT   * dt_h)
    sv$Lurbi_C  <- max(0, sv$Lurbi_C  + dLurbiC  * dt_h)
    sv$Topo_C   <- max(0, sv$Topo_C   + dTopoC   * dt_h)
    sv$Tarlat_C <- max(0, sv$Tarlat_C + dTarlatC * dt_h)
    sv$Tarlat_T <- max(0, sv$Tarlat_T + dTarlatT * dt_h)
    sv$Tarlat_P <- max(0, sv$Tarlat_P + dTarlatP * dt_h)
    sv$Trila_C  <- max(0, sv$Trila_C  + dTrilaC  * dt_h)
    sv$DNA      <- max(0, sv$DNA      + dDNA     * dt_h)
    sv$Ts       <- max(0, sv$Ts       + dTs      * dt_h)
    sv$Ta       <- max(0, sv$Ta       + dTa      * dt_h)
    sv$Tr       <- max(0, sv$Tr       + dTr      * dt_h)
    sv$CD8T     <- max(0.01, sv$CD8T  + dCD8T    * dt_h)
    sv$PDL1_occ <- pmin(1, max(0, sv$PDL1_occ + dPDL1 * dt_h))
    sv$NK       <- max(0.01, sv$NK    + dNK      * dt_h)
    sv$IFNg     <- max(0, sv$IFNg     + dIFNg    * dt_h)

    Ttot2 <- max(sv$Ts + sv$Ta + sv$Tr, 0)
    out[i,] <- list(
      time_h=t_h, time_d=t_h/24,
      Tumor_total=Ttot2,
      SLD=8 * (Ttot2 / p$Tum_S0)^0.333,
      DNA=sv$DNA, CD8T=sv$CD8T,
      PDL1_occ=sv$PDL1_occ,
      Atezo_C=sv$Atezo_C, Durva_C=sv$Durva_C,
      Tarlat_C=sv$Tarlat_C, Trila_C=sv$Trila_C,
      IFNg=sv$IFNg, NK=sv$NK
    )
  }
  out
}

# ---- Default parameters ----
default_params <- function() {
  list(
    Atezo_CL=0.200, Atezo_V1=3.28, Atezo_Q=0.496, Atezo_V2=2.75,
    Durva_CL=0.213, Durva_V1=3.61, Durva_Q=0.468, Durva_V2=3.16,
    Carbo_CL=7.0,   Carbo_V=18.0,  Carbo_dose=500,
    Etop_CL=1.21,   Etop_V1=7.20,  Etop_Q=0.85,   Etop_V2=14.0, Etop_dose=170,
    Lurbi_CL=16.2,  Lurbi_V=489.0, Lurbi_dose=5.44,
    Topo_CL=17.5,   Topo_V=87.0,   Topo_dose=2.55,
    Tarlat_CL=0.048,Tarlat_V1=2.20,Tarlat_Q1=0.190,Tarlat_Vt=0.80,
    Tarlat_Q2=0.310,Tarlat_V2=4.10,Tarlat_dose=10, Tarlat_EC50=0.50,
    Trila_CL=22.4,  Trila_V=110.0, Trila_dose=408, Trila_EC50=15.0, Trila_Emax=0.95,
    kTumor_g=0.0048, Tum_cap=3000, Tum_S0=60, kResist=0.00006,
    Kd_PDL1=0.30, PDL1_base=1.0, PDL1_IFN=0.40, k_PDL1_on=0.20,
    kDNA_form=0.050, kDNA_rep=0.025, DNA_kill=0.030, EC50_DNA=2.50,
    kCD8_stim=0.030, kCD8_kill=0.012, kCD8_death=0.008,
    kNK_act=0.025, kIFNg_prod=0.020, kIFNg_deg=0.040,
    Tarlat_Emax=0.90, DLL3_expr=1.0, kBiTE_kill=0.035,
    f_Atezo=0,f_Durva=0,f_Carbo=0,f_Etop=0,f_Lurbi=0,
    f_Topo=0,f_Tarlat=0,f_Trila=0,f_Nivo=0,
    n_CE_cycles=6, CE_tau_d=21, n_IO_cycles=16, n_Tarlat_cycles=24,
    n_Lurbi_cycles=6, n_Topo_cycles=6
  )
}

# ---- UI ----
ui <- dashboardPage(
  dashboardHeader(title = "SCLC QSP Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",   tabName="tab_patient",  icon=icon("user")),
      menuItem("Drug PK",           tabName="tab_pk",       icon=icon("pills")),
      menuItem("Molecular Pathways",tabName="tab_mol",      icon=icon("dna")),
      menuItem("Clinical Endpoints",tabName="tab_clinical", icon=icon("chart-line")),
      menuItem("Scenario Comparison",tabName="tab_compare", icon=icon("balance-scale")),
      menuItem("Biomarkers & TME",  tabName="tab_bio",      icon=icon("microscope"))
    )
  ),
  dashboardBody(
    tabItems(

      # ========== TAB 1: Patient Profile ==========
      tabItem(tabName="tab_patient",
        fluidRow(
          box(title="Patient Characteristics", width=4, status="primary",
            selectInput("stage", "Staging", choices=c("Extensive stage (ES-SCLC)","Limited stage (LS-SCLC)"),
                        selected="Extensive stage (ES-SCLC)"),
            sliderInput("ecog", "ECOG PS", 0, 2, 1, step=1),
            sliderInput("Tum_S0", "Initial Tumor Burden (x10^9 cells)", 10, 200, 60, step=5),
            selectInput("subtype","NE Subtype",
                        choices=c("ASCL1-high (NE-high)","NEUROD1 (NE-high)","YAP1 (NE-low)","POU2F3 (Tuft)"),
                        selected="ASCL1-high (NE-high)"),
            sliderInput("DLL3_expr","DLL3 Expression (0=absent, 1=high)", 0, 1, 1.0, step=0.1),
            numericInput("bsa", "Body Surface Area (m²)", value=1.7, min=1.2, max=2.5, step=0.1)
          ),
          box(title="Comorbidities & Biomarkers", width=4, status="warning",
            checkboxInput("smoker", "Smoker (pack-year >30)", TRUE),
            checkboxInput("brain_mets", "Brain Metastases", FALSE),
            checkboxInput("ldh_high", "LDH > ULN", FALSE),
            sliderInput("TMB","Tumor Mutational Burden (mut/Mb)", 0, 50, 8, step=1),
            sliderInput("PDL1_base","Baseline PD-L1 TPS (%)", 0, 100, 15, step=5),
            numericInput("GFR","eGFR (mL/min)", 80, 20, 130, step=5)
          ),
          box(title="Patient Summary", width=4, status="success",
            h4("Key Risk Factors"),
            verbatimTextOutput("patient_summary"),
            hr(),
            h4("Genomic Profile"),
            tableOutput("genomic_table")
          )
        )
      ),

      # ========== TAB 2: Drug PK ==========
      tabItem(tabName="tab_pk",
        fluidRow(
          box(title="Treatment Regimen", width=3, status="primary",
            h4("1L Options"),
            checkboxInput("use_Carbo", "Carboplatin (AUC5)", TRUE),
            checkboxInput("use_Etop",  "Etoposide (100 mg/m²)", TRUE),
            checkboxInput("use_Atezo", "Atezolizumab 1200 mg q3w", FALSE),
            checkboxInput("use_Durva", "Durvalumab 1500 mg q4w", FALSE),
            checkboxInput("use_Trila", "Trilaciclib (CDK4/6i)", FALSE),
            hr(),
            h4("2L Options"),
            checkboxInput("use_Tarlat","Tarlatamab 10 mg q2w", FALSE),
            checkboxInput("use_Lurbi", "Lurbinectedin 3.2 mg/m²", FALSE),
            checkboxInput("use_Topo",  "Topotecan 1.5 mg/m²", FALSE),
            hr(),
            sliderInput("sim_days_pk","Simulation Duration (days)", 30, 365, 84),
            actionButton("run_pk","Run Simulation", class="btn-primary btn-block")
          ),
          box(title="Plasma Concentration — Immunotherapy", width=9, status="info",
            plotOutput("pk_IO_plot", height="300px"),
            plotOutput("pk_Tarlat_plot", height="300px")
          )
        ),
        fluidRow(
          box(title="Plasma Concentration — Chemotherapy", width=6, status="warning",
            plotOutput("pk_chemo_plot", height="300px")),
          box(title="CDK4/6 Inhibitor (Trilaciclib)", width=6, status="success",
            plotOutput("pk_trila_plot", height="300px"))
        )
      ),

      # ========== TAB 3: Molecular Pathways ==========
      tabItem(tabName="tab_mol",
        fluidRow(
          box(title="DNA Damage Response", width=6, status="danger",
            plotOutput("mol_dna_plot", height="300px"),
            p("Etoposide & carboplatin generate DSBs; lurbinectedin alkylates DNA.
               Topotecan traps TOP1-cleavage complexes (SSBs). Repair kinetics
               modulate killing efficiency.")
          ),
          box(title="Tarlatamab BiTE Effect on DLL3 Pathway", width=6, status="primary",
            plotOutput("mol_bite_plot", height="300px"),
            p("DLL3 (aberrant Notch ligand) is expressed in ~80% of SCLC.
               Tarlatamab bridges DLL3+ tumor cells to CD3+ T cells,
               inducing granzyme B-mediated killing.")
          )
        ),
        fluidRow(
          box(title="Immune Checkpoint Occupancy (PD-L1 RO)", width=6, status="info",
            plotOutput("mol_pdl1_plot", height="250px")),
          box(title="IFN-gamma & NK Activation", width=6, status="success",
            plotOutput("mol_ifng_plot", height="250px"))
        )
      ),

      # ========== TAB 4: Clinical Endpoints ==========
      tabItem(tabName="tab_clinical",
        fluidRow(
          box(title="Tumor Burden (SLD, cm)", width=8, status="primary",
            plotOutput("clin_tumor_plot", height="350px")),
          box(title="Response Summary", width=4, status="success",
            tableOutput("response_table"),
            hr(),
            h5("Reference Clinical Data"),
            tableOutput("clin_ref_table")
          )
        ),
        fluidRow(
          box(title="Tumor Cell Dynamics (by compartment)", width=12, status="info",
            plotOutput("clin_tumor_comp_plot", height="280px"))
        )
      ),

      # ========== TAB 5: Scenario Comparison ==========
      tabItem(tabName="tab_compare",
        fluidRow(
          box(title="Select Scenarios", width=3, status="primary",
            checkboxGroupInput("scen_select", "Scenarios to Compare:",
              choices = c(
                "Untreated"            = "untreated",
                "CE alone"             = "CE",
                "CE + Atezolizumab"    = "CE_Atezo",
                "CE + Durvalumab"      = "CE_Durva",
                "Tarlatamab 2L"        = "Tarlat",
                "Lurbinectedin 2L"     = "Lurbi",
                "Topotecan 2L"         = "Topo",
                "CE + Trilaciclib"     = "CE_Trila"
              ),
              selected = c("untreated","CE","CE_Atezo","Tarlat")
            ),
            sliderInput("sim_days_comp","Duration (days)", 90, 730, 365),
            actionButton("run_compare","Run All Scenarios", class="btn-success btn-block")
          ),
          box(title="Tumor Burden Comparison", width=9, status="info",
            plotOutput("compare_tumor_plot", height="400px"))
        ),
        fluidRow(
          box(title="CD8+ T Cell Activity", width=6, status="success",
            plotOutput("compare_cd8_plot", height="280px")),
          box(title="Comparative Efficacy Table", width=6, status="warning",
            tableOutput("compare_efficacy_table"))
        )
      ),

      # ========== TAB 6: Biomarkers & TME ==========
      tabItem(tabName="tab_bio",
        fluidRow(
          box(title="CD8+ T Cells & NK Activity", width=6, status="primary",
            plotOutput("bio_immune_plot", height="300px")),
          box(title="IFN-gamma Dynamics", width=6, status="info",
            plotOutput("bio_ifng_plot", height="300px"))
        ),
        fluidRow(
          box(title="Neuroendocrine Biomarkers", width=6, status="warning",
            plotOutput("bio_ne_plot", height="280px"),
            p("ProGRP and NSE levels reflect NE subtype activity.
               Both serve as pharmacodynamic markers for SCLC treatment.")
          ),
          box(title="NE Subtype & DLL3 Impact on Efficacy", width=6, status="success",
            plotOutput("bio_dll3_plot", height="280px"),
            p("Higher DLL3 expression predicts greater tarlatamab response.
               DLL3 is expressed in ~80% ES-SCLC at primary and ~95% at relapse.")
          )
        )
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {

  p_base <- reactive({
    p <- default_params()
    p$Tum_S0   <- input$Tum_S0
    p$DLL3_expr<- input$DLL3_expr
    p
  })

  # Reactive simulation (single scenario)
  sim_result <- eventReactive(input$run_pk, {
    p <- p_base()
    p$f_Carbo  <- as.numeric(input$use_Carbo)
    p$f_Etop   <- as.numeric(input$use_Etop)
    p$f_Atezo  <- as.numeric(input$use_Atezo)
    p$f_Durva  <- as.numeric(input$use_Durva)
    p$f_Trila  <- as.numeric(input$use_Trila)
    p$f_Tarlat <- as.numeric(input$use_Tarlat)
    p$f_Lurbi  <- as.numeric(input$use_Lurbi)
    p$f_Topo   <- as.numeric(input$use_Topo)
    simulate_sclc(p, t_max_days = input$sim_days_pk)
  }, ignoreNULL = FALSE)

  # Scenario comparison
  compare_result <- eventReactive(input$run_compare, {
    scenarios <- list(
      untreated = list(f_Atezo=0,f_Durva=0,f_Carbo=0,f_Etop=0,f_Lurbi=0,f_Topo=0,f_Tarlat=0,f_Trila=0),
      CE        = list(f_Carbo=1,f_Etop=1),
      CE_Atezo  = list(f_Carbo=1,f_Etop=1,f_Atezo=1),
      CE_Durva  = list(f_Carbo=1,f_Etop=1,f_Durva=1),
      Tarlat    = list(f_Tarlat=1),
      Lurbi     = list(f_Lurbi=1),
      Topo      = list(f_Topo=1),
      CE_Trila  = list(f_Carbo=1,f_Etop=1,f_Trila=1)
    )
    labels <- c(
      untreated="Untreated", CE="CE alone", CE_Atezo="CE + Atezolizumab (IMpower133)",
      CE_Durva="CE + Durvalumab (CASPIAN)", Tarlat="Tarlatamab 2L (DeLLphi-301)",
      Lurbi="Lurbinectedin 2L", Topo="Topotecan 2L", CE_Trila="CE + Trilaciclib"
    )
    selected <- input$scen_select
    result_list <- lapply(selected, function(s) {
      p <- p_base()
      p$Tum_S0 <- if (s %in% c("Tarlat","Lurbi","Topo")) 30 else p$Tum_S0
      flags <- scenarios[[s]]
      for (nm in names(flags)) p[[nm]] <- flags[[nm]]
      df <- simulate_sclc(p, t_max_days = input$sim_days_comp)
      df$Scenario <- labels[[s]]
      df
    })
    bind_rows(result_list)
  }, ignoreNULL = FALSE)

  # ---- Tab 1: Patient summary ----
  output$patient_summary <- renderText({
    risk <- "Standard risk"
    if (input$ldh_high) risk <- "High risk (LDH↑)"
    if (input$brain_mets) risk <- "Very high risk (BM+)"
    paste0(
      "Stage: ", input$stage, "\n",
      "ECOG PS: ", input$ecog, "\n",
      "NE Subtype: ", input$subtype, "\n",
      "DLL3: ", round(input$DLL3_expr * 100), "% expression\n",
      "TMB: ", input$TMB, " mut/Mb\n",
      "Risk: ", risk, "\n",
      "BSA: ", input$bsa, " m²"
    )
  })

  output$genomic_table <- renderTable({
    data.frame(
      Alteration = c("TP53","RB1","MYC/L/N","CREBBP","PTEN"),
      Frequency  = c("~90%","~85%","~30%","~22%","~10%"),
      Therapeutic_Target = c("No","No","Indirect","No","PI3Ki")
    )
  })

  # ---- Tab 2: PK plots ----
  output$pk_IO_plot <- renderPlot({
    df <- sim_result()
    df_long <- df %>%
      select(time_d, Atezo_C, Durva_C) %>%
      pivot_longer(-time_d, names_to="Drug", values_to="Conc")
    ggplot(df_long, aes(x=time_d, y=Conc, color=Drug)) +
      geom_line(linewidth=1.1) +
      labs(title="Immunotherapy Plasma Concentration",
           x="Time (days)", y="Concentration (µg/mL)") +
      scale_color_manual(values=c(Atezo_C="#1565C0",Durva_C="#283593"),
                         labels=c("Atezolizumab","Durvalumab")) +
      theme_bw()
  })

  output$pk_Tarlat_plot <- renderPlot({
    df <- sim_result()
    df_long <- df %>%
      select(time_d, Tarlat_C) %>%
      pivot_longer(-time_d, names_to="Drug", values_to="Conc")
    ggplot(df_long, aes(x=time_d, y=Conc, color=Drug)) +
      geom_line(linewidth=1.1, color="#1A237E") +
      labs(title="Tarlatamab Central Plasma", x="Time (days)", y="Concentration (µg/mL)") +
      theme_bw()
  })

  output$pk_chemo_plot <- renderPlot({
    df <- sim_result()
    df_long <- df %>%
      select(time_d, Atezo_C, Durva_C, Tarlat_C) %>%
      { . } # placeholder — show DNA damage as chemo proxy
    ggplot(df, aes(x=time_d, y=DNA)) +
      geom_line(linewidth=1.1, color="#FF6F00") +
      labs(title="DNA Damage Index (Chemo PD Marker)",
           x="Time (days)", y="DNA Damage Index") +
      theme_bw()
  })

  output$pk_trila_plot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x=time_d, y=Trila_C)) +
      geom_line(linewidth=1.1, color="#E65100") +
      labs(title="Trilaciclib Plasma Concentration",
           x="Time (days)", y="Trilaciclib (ng/mL)") +
      theme_bw()
  })

  # ---- Tab 3: Molecular plots ----
  output$mol_dna_plot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x=time_d, y=DNA)) +
      geom_line(linewidth=1.2, color="#FF6F00") +
      geom_area(fill="#FF6F00", alpha=0.2) +
      labs(title="DNA Damage Index Over Time",
           x="Time (days)", y="DNA Damage Index") +
      theme_bw()
  })

  output$mol_bite_plot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x=time_d, y=Tarlat_C * input$DLL3_expr)) +
      geom_line(linewidth=1.2, color="#1A237E") +
      labs(title="Tarlatamab DLL3-binding Activity (tumor)",
           x="Time (days)", y="Effective BiTE Activity (DLL3-weighted conc.)") +
      theme_bw()
  })

  output$mol_pdl1_plot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x=time_d, y=PDL1_occ)) +
      geom_line(linewidth=1.2, color="#006064") +
      scale_y_continuous(limits=c(0,1), labels=scales::percent) +
      labs(title="PD-L1 Receptor Occupancy",
           x="Time (days)", y="PD-L1 RO (%)") +
      theme_bw()
  })

  output$mol_ifng_plot <- renderPlot({
    df <- sim_result()
    df_long <- df %>% select(time_d, IFNg, NK) %>%
      pivot_longer(-time_d, names_to="Signal", values_to="Value")
    ggplot(df_long, aes(x=time_d, y=Value, color=Signal)) +
      geom_line(linewidth=1.1) +
      scale_color_manual(values=c(IFNg="#00BCD4", NK="#006064")) +
      labs(title="IFN-γ & NK Cell Activity",
           x="Time (days)", y="Relative Activity") +
      theme_bw()
  })

  # ---- Tab 4: Clinical ----
  output$clin_tumor_plot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x=time_d, y=SLD)) +
      geom_line(linewidth=1.3, color="#1565C0") +
      geom_hline(yintercept=8, linetype="dashed", color="gray50") +
      annotate("text", x=max(df$time_d)*0.7, y=8.5, label="Baseline SLD", size=3) +
      labs(title="Tumor Sum of Longest Diameters (SLD)",
           x="Time (days)", y="SLD (cm)") +
      theme_bw()
  })

  output$response_table <- renderTable({
    df <- sim_result()
    last  <- tail(df$SLD, 1)
    nadir <- min(df$SLD)
    base  <- df$SLD[1]
    resp  <- (nadir - base) / base * 100
    resp_cat <- if (resp < -30) "PR/CR" else if (resp < 0) "SD" else "PD"
    data.frame(
      Metric = c("Baseline SLD","Nadir SLD","% Change","Best Response","Final SLD"),
      Value  = c(round(base,1), round(nadir,1),
                 paste0(round(resp,1),"%"), resp_cat, round(last,1))
    )
  })

  output$clin_ref_table <- renderTable({
    data.frame(
      Trial = c("IMpower133","CASPIAN","DeLLphi-301","ATLANTIS"),
      Regimen = c("CE+Atezo","CE+Durva","Tarlatamab","Lurbi"),
      mOS_mo = c("13.9","13.0","14.3*","5.3"),
      mPFS_mo = c("5.2","5.1","4.9","3.5")
    )
  })

  output$clin_tumor_comp_plot <- renderPlot({
    df <- sim_result()
    df_long <- df %>%
      select(time_d, Tumor_S, Tumor_A, Tumor_R) %>%
      rename(`Sensitive`=Tumor_S, `G1-Arrested`=Tumor_A, `Resistant`=Tumor_R) %>%
      pivot_longer(-time_d, names_to="Compartment", values_to="Cells_billion")
    ggplot(df_long, aes(x=time_d, y=Cells_billion, fill=Compartment)) +
      geom_area(alpha=0.75) +
      scale_fill_manual(values=c(Sensitive="#455A64",`G1-Arrested`="#607D8B",Resistant="#263238")) +
      labs(title="Tumor Cell Dynamics by Compartment",
           x="Time (days)", y="Tumor Cells (×10⁹)") +
      theme_bw()
  })

  # ---- Tab 5: Scenario Comparison ----
  output$compare_tumor_plot <- renderPlot({
    df <- compare_result()
    ggplot(df, aes(x=time_d, y=SLD, color=Scenario)) +
      geom_line(linewidth=1.1, alpha=0.85) +
      labs(title="Tumor Burden — Scenario Comparison",
           x="Time (days)", y="SLD (cm)") +
      scale_color_brewer(palette="Set1") +
      theme_bw() +
      theme(legend.position="bottom", legend.text=element_text(size=8))
  })

  output$compare_cd8_plot <- renderPlot({
    df <- compare_result()
    ggplot(df, aes(x=time_d, y=CD8T, color=Scenario)) +
      geom_line(linewidth=1.0, alpha=0.8) +
      labs(title="CD8+ T Cell Activity",
           x="Time (days)", y="CD8T (relative)") +
      scale_color_brewer(palette="Set1") +
      theme_bw() + theme(legend.position="bottom")
  })

  output$compare_efficacy_table <- renderTable({
    df <- compare_result()
    df %>%
      group_by(Scenario) %>%
      summarise(
        `Nadir SLD (cm)` = round(min(SLD, na.rm=TRUE), 1),
        `% Reduction`    = round((min(SLD)-head(SLD,1))/head(SLD,1)*100, 1),
        `Peak CD8T`      = round(max(CD8T, na.rm=TRUE), 2),
        .groups="drop"
      )
  })

  # ---- Tab 6: Biomarkers ----
  output$bio_immune_plot <- renderPlot({
    df <- sim_result()
    df_long <- df %>% select(time_d, CD8T, NK) %>%
      pivot_longer(-time_d, names_to="Cell", values_to="Activity")
    ggplot(df_long, aes(x=time_d, y=Activity, color=Cell)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=c(CD8T="#00ACC1", NK="#006064")) +
      labs(title="Immune Cell Activity", x="Time (days)", y="Relative Activity") +
      theme_bw()
  })

  output$bio_ifng_plot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(x=time_d, y=IFNg)) +
      geom_line(linewidth=1.2, color="#00BCD4") +
      geom_area(fill="#00BCD4", alpha=0.2) +
      labs(title="IFN-γ Level Over Time",
           x="Time (days)", y="IFN-γ (relative)") +
      theme_bw()
  })

  output$bio_ne_plot <- renderPlot({
    # Surrogate NE markers from tumor burden
    df <- sim_result()
    df$ProGRP <- df$Tumor_total * 2.5 + rnorm(nrow(df), 0, 0.5)
    df$NSE    <- df$Tumor_total * 1.2 + rnorm(nrow(df), 0, 0.3)
    df_long <- df %>% select(time_d, ProGRP, NSE) %>%
      pivot_longer(-time_d, names_to="Biomarker", values_to="Level")
    ggplot(df_long, aes(x=time_d, y=Level, color=Biomarker)) +
      geom_line(linewidth=1.1) +
      scale_color_manual(values=c(ProGRP="#9C27B0",NSE="#FF9800")) +
      labs(title="Neuroendocrine Biomarkers",
           x="Time (days)", y="Level (relative units)") +
      theme_bw()
  })

  output$bio_dll3_plot <- renderPlot({
    dll3_vals <- c(0.2, 0.5, 0.8, 1.0)
    labels    <- paste0("DLL3 ", dll3_vals*100, "%")
    res <- lapply(seq_along(dll3_vals), function(j) {
      p <- p_base()
      p$f_Tarlat  <- 1
      p$DLL3_expr <- dll3_vals[j]
      p$Tum_S0    <- 30
      df <- simulate_sclc(p, t_max_days=180)
      df$DLL3_level <- labels[j]
      df
    })
    df_all <- bind_rows(res)
    ggplot(df_all, aes(x=time_d, y=SLD, color=DLL3_level)) +
      geom_line(linewidth=1.1) +
      scale_color_manual(values=c("#FFCDD2","#EF9A9A","#EF5350","#B71C1C")) +
      labs(title="Tarlatamab Efficacy vs DLL3 Expression",
           x="Time (days)", y="SLD (cm)", color="DLL3 Expression") +
      theme_bw()
  })
}

shinyApp(ui, server)
