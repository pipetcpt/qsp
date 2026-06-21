## ============================================================
## Chronic Gastritis QSP – Interactive Shiny Dashboard
## 7 Tabs:
##   1. Patient Profile & Disease Stage
##   2. H. pylori Dynamics (Eradication)
##   3. Inflammatory Response (Cytokines)
##   4. Gastric Physiology (Acid / Mucus / Gastrin)
##   5. Correa Cascade (Long-term Progression)
##   6. Treatment Comparison
##   7. Biomarkers & Clinical Endpoints
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

## ---- Embed model code (same as mrgsolve file) ----
cgast_code <- '
$PARAM
Hp0=6.0 Hp_max=8.0 Hp_growth=0.08 Hp_death=0.02
kill_AMX=0.15 kill_CLR=0.20 kill_MTZ=0.12 kill_BSS=0.05
res_CLR=0.0 res_MTZ=0.0
kNFkB_on=0.30 kNFkB_deg=0.25 kIL8_prod=2.5 kIL8_deg=0.18
kIL1b_prod=1.2 kIL1b_deg=0.15 kTNFa_prod=0.80 kTNFa_deg=0.20
kIFNg_prod=0.50 kIFNg_deg=0.12 kIL10_prod=0.40 kIL10_deg=0.10
kNeut_recr=0.35 kNeut_deg=0.25 kTh1_diff=0.12 kTh1_deg=0.08
kTreg_ind=0.08 kTreg_deg=0.06 Treg0=0.10
Gastrin_base=30.0 kGastr_Hp=8.0 kGastr_deg=0.03
Acid_max=15.0 kAcid_stim=0.50 kAcid_deg=0.20 kAcid_PPI=0.80
Mucus_base=200.0 kMucus_deg=0.04 kMucus_prod=0.30
kAtrophy_prog=0.0003 kAtrophy_regr=0.0001
kIM_prog=0.0002 kIM_regr=0.00005
Atrophy0=0.0 IM0=0.0
Ka_PPI=0.50 Ke_PPI=0.40 F_PPI=0.65
Ka_AMX=1.00 Ke_AMX=0.70 F_AMX=0.90
Ka_CLR=0.60 Ke_CLR=0.09 F_CLR=0.55
Ka_MTZ=0.80 Ke_MTZ=0.12 F_MTZ=0.99
w_IL8=0.01 w_TNFa=0.03 w_Acid=0.20 w_Mucus_inv=0.015

$CMT Hp_mucosal NFkB IL8 IL1b TNFa IFNg IL10 Neutrophil Th1 Treg
     Gastrin Acid Mucus Atrophy IM_score PPI_gut PPI_plasma AMX CLR MTZ Symptom

$INIT
Hp_mucosal=6.0 NFkB=0.5 IL8=50.0 IL1b=5.0 TNFa=8.0 IFNg=3.0 IL10=15.0
Neutrophil=0.2 Th1=0.15 Treg=0.10 Gastrin=60.0 Acid=8.0 Mucus=160.0
Atrophy=0.0 IM_score=0.0 PPI_gut=0.0 PPI_plasma=0.0 AMX=0.0 CLR=0.0 MTZ=0.0 Symptom=4.0

$ODE
double Hp_norm=Hp_mucosal/Hp_max;
double Hp_lin=pow(10.0,Hp_mucosal-6.0);
double kill_a=kill_AMX*AMX/(AMX+0.5);
double kill_c=kill_CLR*CLR/(CLR+0.1)*(1.0-res_CLR);
double kill_m=kill_MTZ*MTZ/(MTZ+2.0)*(1.0-res_MTZ);
double kill_total=kill_a+kill_c+kill_m;
dxdt_Hp_mucosal=Hp_growth*(1.0-Hp_norm)-Hp_death-kill_total;
double Hp_stim=(Hp_mucosal>0.0)?Hp_lin:0.0;
dxdt_NFkB=kNFkB_on*Hp_stim*(1.0-NFkB)-kNFkB_deg*NFkB-0.15*IL10*NFkB/(IL10+10.0);
dxdt_IL8=kIL8_prod*NFkB*100.0-kIL8_deg*IL8;
dxdt_IL1b=kIL1b_prod*NFkB*50.0-kIL1b_deg*IL1b;
dxdt_TNFa=kTNFa_prod*NFkB*40.0+0.05*Th1*40.0-kTNFa_deg*TNFa;
dxdt_IFNg=kIFNg_prod*Th1*30.0-kIFNg_deg*IFNg;
dxdt_IL10=kIL10_prod*Treg*50.0+5.0-kIL10_deg*IL10;
dxdt_Neutrophil=kNeut_recr*IL8/(IL8+100.0)-kNeut_deg*Neutrophil;
dxdt_Th1=kTh1_diff*Hp_stim/(Hp_stim+1.0)*(1.0-Th1)*(1.0-0.5*Treg)-kTh1_deg*Th1;
dxdt_Treg=kTreg_ind*IL10/(IL10+20.0)*(1.0-Treg)-kTreg_deg*(Treg-Treg0);
dxdt_Gastrin=kGastr_Hp*Hp_lin+Gastrin_base*kGastr_deg-kGastr_deg*Gastrin;
double PPI_active=PPI_plasma;
double acid_inhibit=kAcid_PPI*PPI_active/(PPI_active+0.5);
double IL1b_inhib=0.10*IL1b/(IL1b+20.0);
double acid_stim=kAcid_stim*Gastrin/(Gastrin+30.0);
dxdt_Acid=Acid_max*acid_stim*(1.0-acid_inhibit-IL1b_inhib)-kAcid_deg*Acid;
double ROS_proxy=Neutrophil+0.5*IL1b/30.0;
dxdt_Mucus=kMucus_prod*(Mucus_base-Mucus)-kMucus_deg*ROS_proxy*Mucus;
double inflam_drive=(IL1b/10.0+TNFa/15.0+IFNg/5.0)/3.0;
double eradicated=(Hp_mucosal<2.0)?1.0:0.0;
dxdt_Atrophy=kAtrophy_prog*inflam_drive*(3.0-Atrophy)-kAtrophy_regr*eradicated*Atrophy;
dxdt_IM_score=kIM_prog*Atrophy*(3.0-IM_score)-kIM_regr*eradicated*IM_score;
dxdt_PPI_gut=-Ka_PPI*PPI_gut;
dxdt_PPI_plasma=Ka_PPI*PPI_gut-Ke_PPI*PPI_plasma;
dxdt_AMX=-Ke_AMX*AMX;
dxdt_CLR=-Ke_CLR*CLR;
dxdt_MTZ=-Ke_MTZ*MTZ;
double raw_score=w_Acid*Acid+w_IL8*IL8/10.0+w_TNFa*TNFa/10.0+w_Mucus_inv*(Mucus_base-Mucus);
dxdt_Symptom=0.5*(raw_score>10.0?10.0:raw_score)-0.5*Symptom;

$TABLE
double PGI_proxy=100.0-25.0*Atrophy;
double PGII_proxy=10.0+5.0*(IL8/50.0);
double PG_ratio=(PGII_proxy>0)?PGI_proxy/PGII_proxy:99.0;
double OLGA_stage=(Atrophy<0.5)?0.0:(Atrophy<1.5)?1.0:(Atrophy<2.5)?2.0:3.0;
double eradicated_flag=(Hp_mucosal<2.0)?1.0:0.0;

$CAPTURE Hp_mucosal NFkB IL8 IL1b TNFa IFNg IL10 Neutrophil Th1 Treg
         Gastrin Acid Mucus Atrophy IM_score PGI_proxy PGII_proxy PG_ratio
         OLGA_stage PPI_plasma AMX CLR MTZ Symptom eradicated_flag
'

cgast_mod <- mcode("cgast_shiny", cgast_code, quiet = TRUE)

## ---- Helper: build event list from UI ----
build_events <- function(use_PPI, use_AMX, use_CLR, use_MTZ,
                         treat_days, ppi_dose, amx_dose, clr_dose, mtz_dose) {
  evs <- list()
  if (use_PPI) evs[["PPI"]] <- ev(cmt = "PPI_gut",
    amt = ppi_dose * 0.65, ii = 12, addl = 2 * treat_days - 1)
  if (use_AMX) evs[["AMX"]] <- ev(cmt = "AMX",
    amt = amx_dose * 0.90 / 10.0, ii = 12, addl = 2 * treat_days - 1)
  if (use_CLR) evs[["CLR"]] <- ev(cmt = "CLR",
    amt = clr_dose * 0.55 / 15.0, ii = 12, addl = 2 * treat_days - 1)
  if (use_MTZ) evs[["MTZ"]] <- ev(cmt = "MTZ",
    amt = mtz_dose * 0.99 / 15.0, ii = 8, addl = 3 * treat_days - 1)
  if (length(evs) > 0) do.call(c, evs) else ev()
}

run_sim <- function(inputs, sim_days = 180) {
  ev_all <- build_events(inputs$use_PPI, inputs$use_AMX, inputs$use_CLR, inputs$use_MTZ,
                         inputs$treat_days, inputs$ppi_dose, inputs$amx_dose,
                         inputs$clr_dose, inputs$mtz_dose)
  ini <- list(
    Hp_mucosal = inputs$Hp_init,
    Atrophy    = inputs$atrophy_init,
    IM_score   = inputs$IM_init,
    Symptom    = 4.0
  )
  cgast_mod %>%
    param(res_CLR = inputs$res_CLR, res_MTZ = inputs$res_MTZ,
          kAcid_PPI = if (inputs$use_VPZ) 0.92 else 0.80) %>%
    do.call(init, c(list(cgast_mod), ini)) %>%
    ev(ev_all) %>%
    mrgsim(end = sim_days * 24, delta = 12) %>%
    as_tibble() %>%
    mutate(time_day = time / 24)
}

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Chronic Gastritis QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",       tabName = "tab_patient",   icon = icon("user")),
      menuItem("H. pylori Dynamics",    tabName = "tab_hp",        icon = icon("bacterium")),
      menuItem("Inflammatory Response", tabName = "tab_inflam",    icon = icon("fire")),
      menuItem("Gastric Physiology",    tabName = "tab_gastric",   icon = icon("flask")),
      menuItem("Correa Cascade",        tabName = "tab_correa",    icon = icon("arrow-trend-up")),
      menuItem("Treatment Comparison",  tabName = "tab_compare",   icon = icon("pills")),
      menuItem("Biomarkers & Endpoints",tabName = "tab_bio",       icon = icon("chart-bar"))
    ),
    hr(),
    h5("  Patient & Disease", style = "color:#aaa; margin-left:10px"),
    sliderInput("Hp_init",       "Initial Hp load (log₁₀ CFU)", 3, 8, 6, 0.5),
    sliderInput("atrophy_init",  "Initial Atrophy Score (0-3)",  0, 3, 0, 0.1),
    sliderInput("IM_init",       "Initial IM Score (0-3)",       0, 3, 0, 0.1),
    hr(),
    h5("  Treatment", style = "color:#aaa; margin-left:10px"),
    checkboxInput("use_PPI", "PPI (Omeprazole)", TRUE),
    checkboxInput("use_VPZ", "Use Vonoprazan instead", FALSE),
    checkboxInput("use_AMX", "Amoxicillin",   TRUE),
    checkboxInput("use_CLR", "Clarithromycin",TRUE),
    checkboxInput("use_MTZ", "Metronidazole", FALSE),
    sliderInput("treat_days", "Treatment Duration (days)", 7, 21, 14, 1),
    numericInput("ppi_dose",  "PPI Dose (mg/dose)", 20, 10, 40, 10),
    numericInput("amx_dose",  "AMX Dose (mg/dose)", 1000, 500, 1500, 250),
    numericInput("clr_dose",  "CLR Dose (mg/dose)", 500, 250, 1000, 250),
    numericInput("mtz_dose",  "MTZ Dose (mg/dose)", 500, 250, 1000, 250),
    hr(),
    h5("  Resistance", style = "color:#aaa; margin-left:10px"),
    sliderInput("res_CLR", "CLR Resistance (0=none, 1=full)", 0, 1, 0, 0.25),
    sliderInput("res_MTZ", "MTZ Resistance (0=none, 1=full)", 0, 1, 0, 0.25),
    hr(),
    sliderInput("sim_days", "Simulation Period (days)", 30, 365, 180, 30),
    actionButton("run_btn", "Run Simulation", class = "btn-success", width = "90%")
  ),

  dashboardBody(
    tabItems(
      ## ---- Tab 1: Patient Profile ----
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Disease Overview", width = 6, status = "primary",
            h4("Chronic Gastritis – H. pylori Mechanistic Overview"),
            p("H. pylori colonizes the gastric mucosa via urease-mediated pH buffering and
               flagellar-driven mucus penetration. CagA (via T4SS injection) and VacA
               toxin drive NF-κB activation, IL-8/IL-1β/TNF-α production, and
               neutrophil infiltration. Chronic inflammation → atrophic gastritis →
               intestinal metaplasia → dysplasia → gastric cancer (Correa cascade)."),
            h4("Correa Cascade Stages"),
            tableOutput("stage_table")
          ),
          box(title = "Initial Patient State", width = 6, status = "info",
            h4("Current Parameters"),
            verbatimTextOutput("patient_summary"),
            h4("Recommended Treatment"),
            uiOutput("recommended_treatment"),
            br(),
            h4("Risk Assessment"),
            uiOutput("risk_assessment")
          )
        ),
        fluidRow(
          box(title = "Virulence Factor Summary", width = 12, status = "warning",
            fluidRow(
              column(4,
                h5("CagA (Type IV secretion)"),
                p("Present in ~70% of Western strains. EPIYA-C motif → SHP2 binding → MAPK/ERK → proliferation, loss of polarity. Major risk factor for gastric cancer."),
                p(strong("Key evidence:"), " Covacci et al. 1999 (Science) – PMID 10428014")
              ),
              column(4,
                h5("VacA (Vacuolating cytotoxin)"),
                p("s1/m1 genotype: highest toxicity. Vacuolation, mitochondrial damage, NLRP3 inflammasome activation → IL-1β. Inhibits T-cell activation."),
                p(strong("Key evidence:"), " Cover et al. 1992 (J Exp Med) – PMID 1607999")
              ),
              column(4,
                h5("OipA / HopQ / BabA"),
                p("OipA → direct NF-κB; HopQ → CEA-mediated CagA injection efficiency; BabA → Lewis b binding → colonisation density determinant."),
                p(strong("Key evidence:"), " Yamaoka et al. 2000 (J Infect Dis) – PMID 10842219")
              )
            )
          )
        )
      ),

      ## ---- Tab 2: H. pylori Dynamics ----
      tabItem(tabName = "tab_hp",
        fluidRow(
          box(title = "H. pylori Mucosal Density", width = 9, status = "danger",
            plotlyOutput("plot_hp", height = "420px")
          ),
          box(title = "Eradication Status", width = 3, status = "info",
            valueBoxOutput("vbox_hp_day14", width = 12),
            valueBoxOutput("vbox_hp_day90", width = 12),
            valueBoxOutput("vbox_eradication", width = 12),
            hr(),
            h5("Eradication Threshold"),
            p("< 2 log₁₀ CFU/biopsy = eradicated"),
            h5("Expected Eradication Rates (ITT)"),
            tableOutput("erad_rate_table")
          )
        ),
        fluidRow(
          box(title = "Drug Concentrations (AMX / CLR / MTZ)", width = 12, status = "primary",
            plotlyOutput("plot_drug_conc", height = "320px")
          )
        )
      ),

      ## ---- Tab 3: Inflammatory Response ----
      tabItem(tabName = "tab_inflam",
        fluidRow(
          box(title = "Cytokine Dynamics", width = 8, status = "danger",
            plotlyOutput("plot_cytokines", height = "400px")
          ),
          box(title = "Immune Cell Populations", width = 4, status = "warning",
            plotlyOutput("plot_immune_cells", height = "400px")
          )
        ),
        fluidRow(
          box(title = "NF-κB Activation & Regulatory Balance", width = 12, status = "primary",
            fluidRow(
              column(6, plotlyOutput("plot_NFkB", height = "300px")),
              column(6, plotlyOutput("plot_Th1_Treg", height = "300px"))
            )
          )
        )
      ),

      ## ---- Tab 4: Gastric Physiology ----
      tabItem(tabName = "tab_gastric",
        fluidRow(
          box(title = "Gastric Acid Output", width = 4, status = "warning",
            plotlyOutput("plot_acid", height = "320px")
          ),
          box(title = "Serum Gastrin (G17)", width = 4, status = "info",
            plotlyOutput("plot_gastrin", height = "320px")
          ),
          box(title = "Mucus Layer Thickness", width = 4, status = "success",
            plotlyOutput("plot_mucus", height = "320px")
          )
        ),
        fluidRow(
          box(title = "H⁺/K⁺-ATPase Inhibition & Acid Suppression", width = 12, status = "primary",
            p("PPI (omeprazole): irreversible covalent binding to active proton pumps. Maximum suppression
               requires 3-5 days of dosing. Vonoprazan achieves faster, more complete suppression via
               potassium-competitive acid blockade (PCAB). Higher intragastric pH (>5) markedly improves
               amoxicillin stability and H. pylori susceptibility."),
            plotlyOutput("plot_PPI_pd", height = "280px")
          )
        )
      ),

      ## ---- Tab 5: Correa Cascade ----
      tabItem(tabName = "tab_correa",
        fluidRow(
          box(title = "Disease Progression Score (Correa Cascade)", width = 9, status = "danger",
            plotlyOutput("plot_correa", height = "420px")
          ),
          box(title = "OLGA Stage", width = 3, status = "warning",
            valueBoxOutput("vbox_olga", width = 12),
            br(),
            h5("OLGA Staging"),
            tableOutput("olga_table"),
            br(),
            h5("Cancer Risk by Stage"),
            tableOutput("cancer_risk_table")
          )
        ),
        fluidRow(
          box(title = "PGI / PGII Ratio – Atrophy Biomarker Trend", width = 12, status = "info",
            plotlyOutput("plot_PG_ratio", height = "300px")
          )
        )
      ),

      ## ---- Tab 6: Treatment Comparison ----
      tabItem(tabName = "tab_compare",
        fluidRow(
          box(title = "Treatment Comparison – H. pylori Eradication", width = 12, status = "primary",
            plotlyOutput("plot_compare_hp", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Treatment Comparison – Symptom Score", width = 6, status = "success",
            plotlyOutput("plot_compare_symptom", height = "350px")
          ),
          box(title = "End-of-treatment Summary Table", width = 6, status = "info",
            DTOutput("compare_table")
          )
        )
      ),

      ## ---- Tab 7: Biomarkers & Clinical Endpoints ----
      tabItem(tabName = "tab_bio",
        fluidRow(
          box(title = "GastroPanel® Biomarker Panel", width = 6, status = "primary",
            plotlyOutput("plot_gastropanel", height = "380px")
          ),
          box(title = "Clinical Endpoints Over Time", width = 6, status = "info",
            plotlyOutput("plot_endpoints", height = "380px")
          )
        ),
        fluidRow(
          box(title = "Biomarker Reference Ranges", width = 12, status = "success",
            DTOutput("biomarker_table")
          )
        )
      )
    )
  )
)

## ============================================================
## Server
## ============================================================
server <- function(input, output, session) {

  ## Reactive simulation
  sim_data <- eventReactive(input$run_btn, {
    inputs <- list(
      Hp_init       = input$Hp_init,
      atrophy_init  = input$atrophy_init,
      IM_init       = input$IM_init,
      use_PPI       = input$use_PPI,
      use_VPZ       = input$use_VPZ,
      use_AMX       = input$use_AMX,
      use_CLR       = input$use_CLR,
      use_MTZ       = input$use_MTZ,
      treat_days    = input$treat_days,
      ppi_dose      = input$ppi_dose,
      amx_dose      = input$amx_dose,
      clr_dose      = input$clr_dose,
      mtz_dose      = input$mtz_dose,
      res_CLR       = input$res_CLR,
      res_MTZ       = input$res_MTZ
    )
    withProgress(message = "Running ODE simulation...", {
      run_sim(inputs, sim_days = input$sim_days)
    })
  }, ignoreNULL = FALSE)

  ## Multi-scenario data (pre-computed for comparison tab)
  compare_data <- eventReactive(input$run_btn, {
    withProgress(message = "Computing all scenarios...", {
      base_inputs <- list(
        Hp_init = input$Hp_init, atrophy_init = input$atrophy_init,
        IM_init = input$IM_init, treat_days = input$treat_days,
        ppi_dose = 20, amx_dose = 1000, clr_dose = 500, mtz_dose = 500,
        res_CLR = input$res_CLR, res_MTZ = input$res_MTZ,
        use_VPZ = FALSE
      )
      scenarios <- list(
        list(label = "No Treatment",       use_PPI=F,use_AMX=F,use_CLR=F,use_MTZ=F),
        list(label = "PPI Only",           use_PPI=T,use_AMX=F,use_CLR=F,use_MTZ=F),
        list(label = "Triple (PPI+AMX+CLR)", use_PPI=T,use_AMX=T,use_CLR=T,use_MTZ=F),
        list(label = "Quad (PPI+AMX+MTZ)", use_PPI=T,use_AMX=T,use_CLR=F,use_MTZ=T),
        list(label = "Bismuth Quad",       use_PPI=T,use_AMX=T,use_CLR=T,use_MTZ=F)
      )
      lapply(scenarios, function(s) {
        inp <- c(base_inputs, s)
        run_sim(inp, sim_days = input$sim_days) %>%
          mutate(scenario = s$label)
      }) %>% bind_rows()
    })
  }, ignoreNULL = FALSE)

  ## ---- Tab 1 outputs ----
  output$stage_table <- renderTable({
    data.frame(
      Stage = c("Normal Mucosa", "Non-atrophic Gastritis", "Atrophic Gastritis",
                "Intestinal Metaplasia (IM)", "Dysplasia (LGD/HGD)", "Gastric Cancer (EGC)"),
      OLGA  = c("0", "0-I", "I-II", "II-III", "III-IV", "IV"),
      Annual_Progression = c("—", "~0.5%/yr", "~1%/yr", "~0.5%/yr→cancer", "~5%/yr→HGD", "—"),
      Key_Biomarker = c("Normal PGI/PGII >3", "Mild PGI/PGII ↓", "PGI <70 + PGI/PGII <3",
                         "PGI/PGII <3 + G17↓", "Biopsy ± p53+", "EGD + biopsy")
    )
  })

  output$patient_summary <- renderText({
    paste0(
      "Hp Load: ", input$Hp_init, " log₁₀ CFU/biopsy\n",
      "Atrophy: ", input$atrophy_init, " (0-3)\n",
      "IM Score: ", input$IM_init, " (0-3)\n",
      "CLR Resistance: ", input$res_CLR * 100, "%\n",
      "MTZ Resistance: ", input$res_MTZ * 100, "%"
    )
  })

  output$recommended_treatment <- renderUI({
    if (input$res_CLR >= 0.5) {
      tags$div(class = "alert alert-danger",
        "CLR resistance ≥50%: Use Bismuth quadruple therapy (PPI + AMX + MTZ + BSS × 14d) or culture-guided therapy")
    } else if (input$res_CLR >= 0.25) {
      tags$div(class = "alert alert-warning",
        "CLR resistance 25-50%: Consider Bismuth quadruple or concomitant quadruple")
    } else {
      tags$div(class = "alert alert-success",
        "CLR resistance <25%: Standard triple (PPI + AMX + CLR × 14d) or Vonoprazan-based triple")
    }
  })

  output$risk_assessment <- renderUI({
    risk_level <- if (input$IM_init >= 1.5) "HIGH" else if (input$atrophy_init >= 1.5) "MODERATE" else "LOW"
    color <- if (risk_level == "HIGH") "danger" else if (risk_level == "MODERATE") "warning" else "success"
    tags$div(class = paste0("alert alert-", color),
      paste0("Gastric cancer risk: ", risk_level,
             " (Atrophy: ", input$atrophy_init, ", IM: ", input$IM_init, ")"))
  })

  ## ---- Tab 2 outputs ----
  output$plot_hp <- renderPlotly({
    req(sim_data())
    p <- sim_data() %>%
      ggplot(aes(x = time_day, y = Hp_mucosal)) +
      geom_line(color = "#2E7D32", linewidth = 1.2) +
      geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
      annotate("text", x = max(sim_data()$time_day) * 0.1, y = 2.3,
               label = "Eradication threshold (10²)", color = "red", size = 3) +
      labs(x = "Time (days)", y = "H. pylori (log₁₀ CFU/biopsy)") +
      theme_bw()
    ggplotly(p)
  })

  output$vbox_hp_day14 <- renderValueBox({
    d <- sim_data() %>% filter(abs(time_day - 14) < 0.6)
    val <- if (nrow(d) > 0) round(mean(d$Hp_mucosal), 2) else "N/A"
    valueBox(val, "Hp at Day 14 (log₁₀)", icon = icon("bacterium"), color = "red")
  })

  output$vbox_hp_day90 <- renderValueBox({
    d <- sim_data() %>% filter(abs(time_day - 90) < 0.6)
    val <- if (nrow(d) > 0) round(mean(d$Hp_mucosal), 2) else "N/A"
    valueBox(val, "Hp at Day 90 (log₁₀)", icon = icon("bacterium"), color = "orange")
  })

  output$vbox_eradication <- renderValueBox({
    d <- sim_data() %>% filter(time_day >= 28)
    pct <- if (nrow(d) > 0) round(mean(d$eradicated_flag) * 100) else 0
    valueBox(paste0(pct, "%"), "Eradication (Day 28+)", icon = icon("check"), color = "green")
  })

  output$erad_rate_table <- renderTable({
    data.frame(
      Regimen = c("Triple (PPI+AMX+CLR)", "Bismuth Quad", "VPZ Triple", "MTZ Quad"),
      ITT_Rate = c("80-85%", "85-92%", "88-94%", "82-88%")
    )
  })

  output$plot_drug_conc <- renderPlotly({
    req(sim_data())
    p <- sim_data() %>%
      select(time_day, AMX, CLR, MTZ) %>%
      pivot_longer(c(AMX, CLR, MTZ), names_to = "drug", values_to = "conc") %>%
      ggplot(aes(x = time_day, y = conc, color = drug)) +
      geom_line(linewidth = 0.9) + facet_wrap(~drug, scales = "free_y") +
      labs(x = "Time (days)", y = "Plasma Conc (a.u.)") + theme_bw()
    ggplotly(p)
  })

  ## ---- Tab 3 outputs ----
  output$plot_cytokines <- renderPlotly({
    req(sim_data())
    p <- sim_data() %>%
      select(time_day, IL8, IL1b, TNFa, IFNg, IL10) %>%
      pivot_longer(c(IL8, IL1b, TNFa, IFNg, IL10), names_to = "cytokine", values_to = "conc") %>%
      ggplot(aes(x = time_day, y = conc, color = cytokine)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Time (days)", y = "Concentration (pg/mL)") + theme_bw()
    ggplotly(p)
  })

  output$plot_immune_cells <- renderPlotly({
    req(sim_data())
    p <- sim_data() %>%
      select(time_day, Neutrophil, Th1, Treg) %>%
      pivot_longer(c(Neutrophil, Th1, Treg), names_to = "cell", values_to = "value") %>%
      ggplot(aes(x = time_day, y = value, color = cell)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Time (days)", y = "Normalized Level (0-1)") + theme_bw()
    ggplotly(p)
  })

  output$plot_NFkB <- renderPlotly({
    req(sim_data())
    p <- sim_data() %>%
      ggplot(aes(x = time_day, y = NFkB)) +
      geom_line(color = "#E91E63", linewidth = 1.2) +
      labs(x = "Time (days)", y = "NF-κB Activation (0-1)") + theme_bw()
    ggplotly(p)
  })

  output$plot_Th1_Treg <- renderPlotly({
    req(sim_data())
    p <- sim_data() %>%
      select(time_day, Th1, Treg) %>%
      pivot_longer(c(Th1, Treg), names_to = "cell", values_to = "value") %>%
      ggplot(aes(x = time_day, y = value, color = cell)) +
      geom_line(linewidth = 1.1) +
      labs(x = "Time (days)", y = "Normalized (0-1)", title = "Th1 / Treg Balance") + theme_bw()
    ggplotly(p)
  })

  ## ---- Tab 4 outputs ----
  output$plot_acid <- renderPlotly({
    req(sim_data())
    p <- sim_data() %>%
      ggplot(aes(x = time_day, y = Acid)) +
      geom_line(color = "#FF6F00", linewidth = 1.2) +
      geom_hline(yintercept = 5, linetype = "dashed", color = "gray50") +
      labs(x = "Time (days)", y = "Acid Output (mmol/h)") + theme_bw()
    ggplotly(p)
  })

  output$plot_gastrin <- renderPlotly({
    req(sim_data())
    p <- sim_data() %>%
      ggplot(aes(x = time_day, y = Gastrin)) +
      geom_line(color = "#6A1B9A", linewidth = 1.2) +
      geom_hline(yintercept = 30, linetype = "dashed", color = "gray50") +
      annotate("text", x = 5, y = 28, label = "Normal baseline ~30 pg/mL", size = 3, color = "gray40") +
      labs(x = "Time (days)", y = "Serum Gastrin (pg/mL)") + theme_bw()
    ggplotly(p)
  })

  output$plot_mucus <- renderPlotly({
    req(sim_data())
    p <- sim_data() %>%
      ggplot(aes(x = time_day, y = Mucus)) +
      geom_line(color = "#00838F", linewidth = 1.2) +
      geom_hline(yintercept = 200, linetype = "dashed", color = "gray50") +
      annotate("text", x = 5, y = 202, label = "Normal 200 μm", size = 3, color = "gray40") +
      labs(x = "Time (days)", y = "Mucus Thickness (μm)") + theme_bw()
    ggplotly(p)
  })

  output$plot_PPI_pd <- renderPlotly({
    req(sim_data())
    p <- sim_data() %>%
      select(time_day, PPI_plasma, Acid) %>%
      pivot_longer(c(PPI_plasma, Acid), names_to = "variable", values_to = "value") %>%
      ggplot(aes(x = time_day, y = value, color = variable)) +
      geom_line(linewidth = 1) + facet_wrap(~variable, scales = "free_y") +
      labs(x = "Time (days)", y = "Value") + theme_bw()
    ggplotly(p)
  })

  ## ---- Tab 5 outputs ----
  output$plot_correa <- renderPlotly({
    req(sim_data())
    p <- sim_data() %>%
      select(time_day, Atrophy, IM_score) %>%
      pivot_longer(c(Atrophy, IM_score), names_to = "marker", values_to = "score") %>%
      ggplot(aes(x = time_day, y = score, color = marker)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c(Atrophy = "#E64A19", IM_score = "#7C4DFF"),
                         labels = c(Atrophy = "Atrophy Score (0-3)", IM_score = "IM Score (0-3)")) +
      labs(x = "Time (days)", y = "Disease Score (0-3)") + theme_bw()
    ggplotly(p)
  })

  output$vbox_olga <- renderValueBox({
    d <- sim_data() %>% tail(1)
    val <- if (nrow(d) > 0) round(d$OLGA_stage) else "N/A"
    valueBox(paste0("Stage ", val), "Current OLGA Stage", icon = icon("stethoscope"), color = "orange")
  })

  output$olga_table <- renderTable({
    data.frame(
      Stage = c("0", "I", "II", "III", "IV"),
      Corpus = c("None", "Mild", "Mod", "Severe", "Severe"),
      Antrum = c("None", "Mild", "Mod", "Mod", "Severe"),
      Cancer_Risk = c("<0.1%/yr", "0.3%/yr", "0.5%/yr", "1.2%/yr", "2.5%/yr")
    )
  })

  output$cancer_risk_table <- renderTable({
    data.frame(
      Feature = c("H. pylori +", "CagA +", "OLGA III-IV", "IM Type III", "Family Hx"),
      OR_Cancer = c("3-6×", "2-3×", "5-10×", "6-8×", "2-4×")
    )
  })

  output$plot_PG_ratio <- renderPlotly({
    req(sim_data())
    p <- sim_data() %>%
      ggplot(aes(x = time_day, y = PG_ratio)) +
      geom_line(color = "#1565C0", linewidth = 1.2) +
      geom_hline(yintercept = 3, linetype = "dashed", color = "red") +
      annotate("text", x = 5, y = 3.3, label = "PGI/PGII < 3 → atrophy screen", color = "red", size = 3) +
      labs(x = "Time (days)", y = "PGI/PGII Ratio") + theme_bw()
    ggplotly(p)
  })

  ## ---- Tab 6 outputs ----
  output$plot_compare_hp <- renderPlotly({
    req(compare_data())
    p <- compare_data() %>%
      ggplot(aes(x = time_day, y = Hp_mucosal, color = scenario)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
      labs(x = "Time (days)", y = "H. pylori (log₁₀ CFU/biopsy)", color = "Regimen") +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p) %>% layout(legend = list(orientation = "h", y = -0.2))
  })

  output$plot_compare_symptom <- renderPlotly({
    req(compare_data())
    p <- compare_data() %>%
      ggplot(aes(x = time_day, y = Symptom, color = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Time (days)", y = "Symptom Score (0-10)") +
      ylim(0, 10) + theme_bw()
    ggplotly(p)
  })

  output$compare_table <- renderDT({
    req(compare_data())
    compare_data() %>%
      group_by(scenario) %>%
      summarise(
        Hp_day14  = round(mean(Hp_mucosal[abs(time_day - 14) < 0.6], na.rm = TRUE), 2),
        Hp_end    = round(last(Hp_mucosal), 2),
        Eradicated = paste0(round(mean(eradicated_flag[time_day >= 28], na.rm = TRUE) * 100), "%"),
        Symptom_end = round(last(Symptom), 2),
        .groups = "drop"
      ) %>%
      datatable(options = list(pageLength = 5, scrollX = TRUE))
  })

  ## ---- Tab 7 outputs ----
  output$plot_gastropanel <- renderPlotly({
    req(sim_data())
    p <- sim_data() %>%
      select(time_day, PGI_proxy, PGII_proxy, Gastrin) %>%
      pivot_longer(c(PGI_proxy, PGII_proxy, Gastrin), names_to = "marker", values_to = "value") %>%
      ggplot(aes(x = time_day, y = value, color = marker)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c(PGI_proxy = "#1E88E5", PGII_proxy = "#E53935", Gastrin = "#7B1FA2"),
                         labels = c(PGI_proxy = "Pepsinogen I (μg/L)",
                                    PGII_proxy = "Pepsinogen II (μg/L)",
                                    Gastrin = "Gastrin-17 (pg/mL)")) +
      labs(x = "Time (days)", y = "Concentration", color = "Marker") + theme_bw()
    ggplotly(p)
  })

  output$plot_endpoints <- renderPlotly({
    req(sim_data())
    p <- sim_data() %>%
      select(time_day, Symptom, OLGA_stage, PG_ratio) %>%
      pivot_longer(c(Symptom, OLGA_stage, PG_ratio), names_to = "endpoint", values_to = "value") %>%
      ggplot(aes(x = time_day, y = value, color = endpoint)) +
      geom_line(linewidth = 1) + facet_wrap(~endpoint, scales = "free_y") +
      labs(x = "Time (days)", y = "Value") + theme_bw()
    ggplotly(p)
  })

  output$biomarker_table <- renderDT({
    data.frame(
      Biomarker   = c("Pepsinogen I", "Pepsinogen II", "PGI/PGII Ratio",
                      "Gastrin-17", "H. pylori IgG", "CagA IgG",
                      "¹³C-UBT", "Stool HpSA"),
      Normal      = c(">70 μg/L", "<10 μg/L", ">3.0", "1-7 pmol/L",
                      "Negative", "Negative", "<3.5 DPM", "Negative"),
      Atrophy     = c("<70 μg/L", ">10 μg/L", "<3.0 (corpus)", "<1 pmol/L (antrum)",
                      "Positive", "High risk if +", "N/A", "N/A"),
      Clinical_Use = c("Corpus function", "Inflammation", "Non-invasive atrophy screen",
                        "Antral function", "Hp status", "Virulence/cancer risk",
                        "Eradication test (4-6w post-Rx)", "Eradication test")
    ) %>%
      datatable(options = list(pageLength = 8, scrollX = TRUE))
  })
}

shinyApp(ui, server)
