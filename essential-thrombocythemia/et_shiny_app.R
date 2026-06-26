## ============================================================
## Essential Thrombocythemia (ET) — Interactive Shiny Dashboard
## 8 Tabs: Patient Profile | Drug PK | Platelet PD |
##         JAK2 Allele Burden | Thrombosis Risk | Scenario Compare |
##         Biomarker Panel | BM & Progression
## ============================================================

library(shiny)
library(shinydashboard)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)

## ── Simulation engine (self-contained Euler ODE) ─────────────
simulate_ET <- function(
    dose_hu    = 0,    # mg/day
    dose_ana   = 0,    # mg/day
    dose_rux   = 0,    # mg/day
    dose_pifn  = 0,    # µg/week
    dose_asa   = 0,    # binary 0/1
    sim_days   = 730,
    JAK2_init  = 0.55,
    PLT_init   = 900,
    ASXL1_pos  = 0,
    age_gt60   = 0,
    prior_thr  = 0,
    dt         = 1.0
) {
  ## Parameters
  k_HSC_prod = 0.010; k_HSC_diff = 0.008; k_MKP_diff = 0.15
  k_MK_mat   = 0.10;  k_PLT_prod = 8.0;   k_PLT_destr = 0.115
  k_TPO_prod = 0.50;  k_TPO_elim = 0.012; k_TPO_abs   = 0.00008
  k_JAK2_exp = 0.003; k_SPL_grow = 0.008; k_SPL_shrk  = 0.05
  phi_JAK2   = 3.5;   PLT_norm   = 250;   TPO_ss      = 85

  EC50_HU = 3.50; Emax_HU = 0.85; gam_HU = 1.20
  EC50_ANA = 0.025; Emax_ANA = 0.75; gam_ANA = 1.50
  EC50_RUX = 0.15;  Emax_RUX = 0.80; gam_RUX = 1.80
  EC50_pIFN = 0.008; Emax_pIFN = 0.70; gam_pIFN = 1.00

  lambda_T  = 0.0005; alpha_T  = 1.80; delta_JAK2 = 0.80
  lambda_MF = 0.00012; delta_ASXL = 1.50

  ## Drug steady-state plasma concentrations (simplified: CSS = F*D/CL/24)
  Css_HU   <- dose_hu  * 0.80 / (4.20 * 24) * 1e3 / 76.05  # µg/mL
  Css_ANA  <- dose_ana * 0.70 / (9.50 * 24) * 1e3           # ng/mL → µg/mL
  Css_ANA  <- Css_ANA / 1000
  Css_RUX  <- dose_rux * 0.95 / (22.0 * 24) * 1e3           # ng/mL → µg/mL
  Css_RUX  <- Css_RUX / 1000
  Css_pIFN <- dose_pifn / 7 * 0.84 / (0.038 * 24) / 1000    # µg/mL

  Hill <- function(C, EC50, Emax, gam)
    Emax * C^gam / (EC50^gam + C^gam)

  E_HU   <- Hill(Css_HU,   EC50_HU,   Emax_HU,   gam_HU)
  E_ANA  <- Hill(Css_ANA,  EC50_ANA,  Emax_ANA,  gam_ANA)
  E_RUX  <- Hill(Css_RUX,  EC50_RUX,  Emax_RUX,  gam_RUX)
  E_pIFN <- Hill(Css_pIFN, EC50_pIFN, Emax_pIFN, gam_pIFN)
  E_ASA  <- ifelse(dose_asa > 0, 0.80, 0)

  ## State vector
  n  <- round(sim_days / dt) + 1
  tm <- seq(0, sim_days, by = dt)

  HSC  <- numeric(n); MKP <- numeric(n); MK  <- numeric(n)
  PLT  <- numeric(n); TPO <- numeric(n); JAK <- numeric(n)
  SPL  <- numeric(n); RT  <- numeric(n); RM  <- numeric(n)

  HSC[1]=1.0; MKP[1]=1.2; MK[1]=1.5; PLT[1]=PLT_init
  TPO[1]=60;  JAK[1]=JAK2_init; SPL[1]=2.5; RT[1]=0; RM[1]=0

  for (i in 2:n) {
    tpo_eff  <- TPO[i-1] / (TPO[i-1] + TPO_ss) * 2.0
    spd      <- max(PLT[i-1]/PLT_norm - 1, 0)
    kd_spl   <- k_PLT_destr * (1 + SPL[i-1]/10)
    thr_mod  <- (1 + age_gt60*0.5) * (1 + prior_thr*0.8)
    mf_mod   <- 1 + ASXL1_pos * delta_ASXL

    dHSC <- k_HSC_prod*(1+JAK[i-1]*phi_JAK2)*1.0 - k_HSC_diff*HSC[i-1]
    dMKP <- k_HSC_diff*HSC[i-1]*tpo_eff - k_MKP_diff*(1+E_HU+E_RUX)*MKP[i-1]
    dMK  <- k_MKP_diff*MKP[i-1] - k_MK_mat*(1+E_ANA*0.6)*MK[i-1] - k_MK_mat*0.4*MK[i-1]
    dPLT <- k_PLT_prod*k_MK_mat*MK[i-1] - kd_spl*PLT[i-1]
    dTPO <- k_TPO_prod - k_TPO_elim*TPO[i-1] - k_TPO_abs*PLT[i-1]*TPO[i-1]
    dJAK <- k_JAK2_exp*JAK[i-1]*(1-JAK[i-1]) - E_pIFN*0.015*JAK[i-1]
    dSPL <- k_SPL_grow*spd - k_SPL_shrk*E_RUX*SPL[i-1]
    dRT  <- lambda_T * PLT[i-1]^alpha_T / PLT_norm^alpha_T *
            (1+delta_JAK2*JAK[i-1]) * thr_mod * (1-E_ASA*0.4)
    dRM  <- lambda_MF * JAK[i-1] * mf_mod

    HSC[i] <- max(HSC[i-1] + dHSC*dt, 0)
    MKP[i] <- max(MKP[i-1] + dMKP*dt, 0)
    MK[i]  <- max(MK[i-1]  + dMK*dt,  0)
    PLT[i] <- max(PLT[i-1] + dPLT*dt, 0)
    TPO[i] <- max(TPO[i-1] + dTPO*dt, 0)
    JAK[i] <- max(min(JAK[i-1] + dJAK*dt, 1), 0)
    SPL[i] <- max(SPL[i-1] + dSPL*dt, 0)
    RT[i]  <- RT[i-1] + dRT*dt
    RM[i]  <- RM[i-1] + dRM*dt
  }

  data.frame(
    time     = tm,
    PLT      = PLT,
    JAK2_AB  = JAK * 100,
    Spleen   = SPL,
    TPO_level= TPO,
    MK_pool  = MK,
    RISK_T   = RT,
    RISK_MF  = RM,
    Hazard_T = lambda_T * PLT^alpha_T/PLT_norm^alpha_T * (1+delta_JAK2*JAK),
    CHR      = as.integer(PLT <= 400),
    PHR      = as.integer(PLT <= 600)
  )
}

## ── Scenarios ────────────────────────────────────────────────
scenario_defs <- list(
  list(id="obs",  label="① No Treatment",                dose_hu=0,    dose_ana=0, dose_rux=0,  dose_pifn=0,  dose_asa=0, color="#616161"),
  list(id="asa",  label="② Aspirin only",                dose_hu=0,    dose_ana=0, dose_rux=0,  dose_pifn=0,  dose_asa=1, color="#ff8f00"),
  list(id="hu5",  label="③ Hydroxyurea 500mg + ASA",     dose_hu=500,  dose_ana=0, dose_rux=0,  dose_pifn=0,  dose_asa=1, color="#1565c0"),
  list(id="hu15", label="④ Hydroxyurea 1500mg + ASA",    dose_hu=1500, dose_ana=0, dose_rux=0,  dose_pifn=0,  dose_asa=1, color="#0288d1"),
  list(id="ana",  label="⑤ Anagrelide 2mg/d + ASA",      dose_hu=0,    dose_ana=2, dose_rux=0,  dose_pifn=0,  dose_asa=1, color="#6a1b9a"),
  list(id="rux",  label="⑥ Ruxolitinib 20mg/d + ASA",    dose_hu=0,    dose_ana=0, dose_rux=20, dose_pifn=0,  dose_asa=1, color="#00695c"),
  list(id="ifn",  label="⑦ Peg-IFN-α2a 90µg/wk + ASA",  dose_hu=0,    dose_ana=0, dose_rux=0,  dose_pifn=90, dose_asa=1, color="#c62828")
)

## ── UI ───────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = span("ET QSP Dashboard", style="font-size:16px")),

  dashboardSidebar(
    sidebarMenu(
      menuItem("① Patient Profile",   tabName="tab_patient",  icon=icon("user")),
      menuItem("② Drug PK",           tabName="tab_pk",       icon=icon("pills")),
      menuItem("③ Platelet Dynamics", tabName="tab_plt",      icon=icon("circle")),
      menuItem("④ JAK2 Allele Burden",tabName="tab_jak2",     icon=icon("dna")),
      menuItem("⑤ Thrombosis Risk",   tabName="tab_risk",     icon=icon("heart-pulse")),
      menuItem("⑥ Scenario Compare",  tabName="tab_compare",  icon=icon("chart-line")),
      menuItem("⑦ Biomarker Panel",   tabName="tab_bio",      icon=icon("microscope")),
      menuItem("⑧ BM & Progression",  tabName="tab_bm",       icon=icon("bone"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color:#f5f6fa; }
      .box-header { border-bottom:2px solid #e0e0e0; }
    "))),

    tabItems(

      ## ── Tab 1: Patient Profile ────────────────────────────
      tabItem("tab_patient",
        fluidRow(
          box(title="Patient Characteristics", width=4, solidHeader=TRUE, status="purple",
            numericInput("PLT_init", "Baseline Platelet Count (×10⁹/L)", 900, 400, 3000, 50),
            sliderInput("JAK2_init", "JAK2 Allele Burden (%)", 0, 100, 55, 1),
            checkboxInput("age_gt60",   "Age > 60 years",         FALSE),
            checkboxInput("prior_thr",  "Prior Thrombosis",        FALSE),
            checkboxInput("ASXL1_pos",  "ASXL1 Co-mutation",      FALSE),
            numericInput("sim_days", "Simulation Duration (days)", 730, 90, 1825, 30)
          ),
          box(title="Risk Stratification (IPSS-ET)", width=4, solidHeader=TRUE, status="warning",
            tableOutput("risk_table")
          ),
          box(title="WHO 2016 Diagnosis Criteria", width=4, solidHeader=TRUE, status="info",
            HTML("<b>Major criteria:</b><ul>
              <li>PLT ≥450×10⁹/L</li>
              <li>BM biopsy: MK proliferation, hypercellularity</li>
              <li>Not meeting WHO criteria for BCR-ABL1+ CML, PV, PMF, MDS, or other MPN</li>
              <li>Presence of JAK2, CALR, or MPL mutation</li></ul>
              <b>Minor criterion:</b>
              <ul><li>Presence of a clonal marker or exclusion of reactive thrombocytosis</li></ul>
              <hr/><b>Diagnosis:</b> All 4 major OR 3 major + 1 minor")
          )
        ),
        fluidRow(
          box(title="Baseline Summary", width=12, solidHeader=TRUE,
            verbatimTextOutput("baseline_summary")
          )
        )
      ),

      ## ── Tab 2: Drug PK ────────────────────────────────────
      tabItem("tab_pk",
        fluidRow(
          box(title="PK Parameters", width=3, solidHeader=TRUE, status="primary",
            h4("Hydroxyurea"),
            sliderInput("dose_hu", "Dose (mg/day)", 0, 2000, 500, 100),
            hr(),
            h4("Anagrelide"),
            sliderInput("dose_ana", "Dose (mg/day)", 0, 4, 0, 0.5),
            hr(),
            h4("Ruxolitinib"),
            sliderInput("dose_rux", "Dose (mg/day)", 0, 40, 0, 5),
            hr(),
            h4("Peg-IFN-α2a"),
            sliderInput("dose_pifn", "Dose (µg/week)", 0, 180, 0, 45),
            hr(),
            checkboxInput("dose_asa", "Low-dose Aspirin (81 mg/d)", TRUE)
          ),
          box(title="Plasma Concentration Profiles (Steady-State Approx)", width=9,
              solidHeader=TRUE, status="primary",
            plotlyOutput("pk_plot", height="500px")
          )
        )
      ),

      ## ── Tab 3: Platelet Dynamics ──────────────────────────
      tabItem("tab_plt",
        fluidRow(
          box(title="Platelet Count Over Time", width=9, solidHeader=TRUE, status="info",
            plotlyOutput("plt_plot", height="450px")
          ),
          box(title="Response Summary", width=3, solidHeader=TRUE, status="success",
            tableOutput("plt_response_table"),
            hr(),
            h5("ELN Response Definitions:"),
            HTML("<ul>
              <li><b>CHR:</b> PLT ≤400 + WBC ≤10 + normal spleen</li>
              <li><b>PHR:</b> PLT ≤600 or ≥50% reduction</li>
              <li><b>NR:</b> Less than PHR</li></ul>")
          )
        )
      ),

      ## ── Tab 4: JAK2 Allele Burden ─────────────────────────
      tabItem("tab_jak2",
        fluidRow(
          box(title="JAK2 V617F Allele Burden Over Time", width=9, solidHeader=TRUE, status="warning",
            plotlyOutput("jak2_plot", height="450px")
          ),
          box(title="Molecular Response", width=3, solidHeader=TRUE, status="warning",
            tableOutput("mol_response_table"),
            hr(),
            HTML("<b>Molecular Response Criteria:</b><ul>
              <li>CMR: JAK2 AB <1%</li>
              <li>PMR: ≥50% AB reduction</li>
              <li>NMR: <50% reduction</li></ul>"),
            hr(),
            HTML("<b>Note:</b> Only Peg-IFN-α reliably reduces JAK2 allele burden.
                  Hydroxyurea and anagrelide have minimal molecular effect.")
          )
        )
      ),

      ## ── Tab 5: Thrombosis Risk ────────────────────────────
      tabItem("tab_risk",
        fluidRow(
          box(title="Annual Thrombosis Hazard Over Time", width=8, solidHeader=TRUE, status="danger",
            plotlyOutput("risk_plot", height="400px")
          ),
          box(title="Cumulative Risk Estimator", width=4, solidHeader=TRUE, status="danger",
            plotlyOutput("cum_risk_plot", height="400px")
          )
        ),
        fluidRow(
          box(title="Thrombosis Risk Factors in ET", width=12, solidHeader=TRUE,
            DT::DTOutput("risk_factors_table")
          )
        )
      ),

      ## ── Tab 6: Scenario Comparison ─────────────────────────
      tabItem("tab_compare",
        fluidRow(
          box(title="Select Scenarios to Compare", width=3, solidHeader=TRUE, status="success",
            checkboxGroupInput("sel_scenarios",
              label = "Scenarios:",
              choices = setNames(sapply(scenario_defs, `[[`, "id"),
                                 sapply(scenario_defs, `[[`, "label")),
              selected = c("obs","hu15","ana","ifn")),
            hr(),
            radioButtons("compare_var", "Variable:",
                         choices = c("Platelet Count"="PLT",
                                     "JAK2 AB (%)"="JAK2_AB",
                                     "Thrombosis Hazard"="Hazard_T",
                                     "Spleen Size"="Spleen"),
                         selected = "PLT")
          ),
          box(title="Scenario Comparison Plot", width=9, solidHeader=TRUE, status="success",
            plotlyOutput("compare_plot", height="450px")
          )
        ),
        fluidRow(
          box(title="Response Rate Comparison at Key Timepoints", width=12,
              solidHeader=TRUE, status="success",
            DT::DTOutput("compare_table")
          )
        )
      ),

      ## ── Tab 7: Biomarker Panel ────────────────────────────
      tabItem("tab_bio",
        fluidRow(
          box(title="Multi-Biomarker Dashboard", width=12, solidHeader=TRUE, status="info",
            plotlyOutput("bio_panel", height="600px")
          )
        ),
        fluidRow(
          box(title="Biomarker Reference Ranges", width=12, solidHeader=TRUE,
            DT::DTOutput("bio_ref_table")
          )
        )
      ),

      ## ── Tab 8: BM & Progression ───────────────────────────
      tabItem("tab_bm",
        fluidRow(
          box(title="Disease Progression Model", width=8, solidHeader=TRUE, status="danger",
            plotlyOutput("progression_plot", height="450px")
          ),
          box(title="MF Transformation Risk Factors", width=4, solidHeader=TRUE, status="danger",
            HTML("<b>IPSS-ET High Risk Features:</b><ul>
              <li>Age >60 y (1 point each)</li>
              <li>Prior thrombosis (2 points)</li>
              <li>JAK2 V617F (1 point)</li>
              <li>Cardiovascular risk factors (1 point)</li>
              <li>Leukocytosis (WBC >11) (1 point)</li></ul>
              <hr/>
              <b>Risk of Post-ET MF:</b><ul>
              <li>Low risk: ~5% at 10 yr</li>
              <li>Int risk: ~12% at 10 yr</li>
              <li>High risk: ~20% at 10 yr</li></ul>
              <hr/>
              <b>Risk of AML transformation:</b> ~1–3% cumulative"),
            hr(),
            tableOutput("mf_risk_summary")
          )
        ),
        fluidRow(
          box(title="Molecular Pathways to Transformation", width=12, solidHeader=TRUE,
            HTML("<div style='display:flex;gap:30px'>
              <div><h5>Post-ET MF</h5><ul>
                <li>ASXL1, EZH2 co-mutations</li>
                <li>TGF-β driven fibroblast activation</li>
                <li>Reticulin → collagen fibrosis</li>
                <li>Extramedullary hematopoiesis</li></ul></div>
              <div><h5>Blast Phase / AML</h5><ul>
                <li>TP53 mutation (therapy-related)</li>
                <li>IDH1/IDH2 mutation</li>
                <li>Cytogenetic evolution</li>
                <li>Hydroxyurea (rare, long-term)</li></ul></div>
              <div><h5>Protective Effects</h5><ul>
                <li>Peg-IFN: reduces JAK2+ clone</li>
                <li>Ruxolitinib: ↓cytokines + fibrosis</li>
                <li>Careful monitoring: BM biopsy q3-5yr</li></ul></div>
            </div>")
          )
        )
      )
    )
  )
)

## ── Server ───────────────────────────────────────────────────
server <- function(input, output, session) {

  ## Reactive: run simulation with current inputs
  sim_data <- reactive({
    simulate_ET(
      dose_hu   = input$dose_hu,
      dose_ana  = input$dose_ana,
      dose_rux  = input$dose_rux,
      dose_pifn = input$dose_pifn,
      dose_asa  = as.numeric(input$dose_asa),
      sim_days  = input$sim_days,
      JAK2_init = input$JAK2_init / 100,
      PLT_init  = input$PLT_init,
      ASXL1_pos = as.numeric(input$ASXL1_pos),
      age_gt60  = as.numeric(input$age_gt60),
      prior_thr = as.numeric(input$prior_thr)
    )
  })

  ## All scenarios for comparison
  all_scenarios_data <- reactive({
    purrr::map_dfr(scenario_defs, function(sc) {
      simulate_ET(dose_hu=sc$dose_hu, dose_ana=sc$dose_ana,
                  dose_rux=sc$dose_rux, dose_pifn=sc$dose_pifn,
                  dose_asa=sc$dose_asa, sim_days=input$sim_days,
                  JAK2_init=input$JAK2_init/100, PLT_init=input$PLT_init,
                  ASXL1_pos=as.numeric(input$ASXL1_pos),
                  age_gt60=as.numeric(input$age_gt60),
                  prior_thr=as.numeric(input$prior_thr)) %>%
        mutate(scenario=sc$label, color=sc$color, id=sc$id)
    })
  })

  ## Risk table
  output$risk_table <- renderTable({
    pts <- 0
    if (input$age_gt60)  pts <- pts + 1
    if (input$prior_thr) pts <- pts + 2
    if (input$JAK2_init >= 50) pts <- pts + 1
    data.frame(
      Factor = c("Age >60","Prior Thrombosis","JAK2 V617F ≥50%","Total Points","Risk Category"),
      Value  = c(as.integer(input$age_gt60), as.integer(input$prior_thr)*2,
                 as.integer(input$JAK2_init >= 50), pts,
                 ifelse(pts==0,"Low",ifelse(pts<=2,"Intermediate","High")))
    )
  })

  output$baseline_summary <- renderText({
    plt   <- input$PLT_init
    jak2  <- input$JAK2_init
    pts   <- as.numeric(input$age_gt60) + as.numeric(input$prior_thr)*2 +
             as.numeric(jak2 >= 50)
    risk  <- ifelse(pts==0,"LOW",ifelse(pts<=2,"INTERMEDIATE","HIGH"))
    thrombo <- ifelse(risk=="LOW","~1-2%/yr",ifelse(risk=="INTERMEDIATE","~3-4%/yr","~4-5%/yr"))
    paste0(
      "=== Patient Summary ===\n",
      sprintf("Baseline PLT:       %d × 10⁹/L (%s)\n", plt,
              ifelse(plt<400,"Normal range",ifelse(plt<600,"Mild thrombocytosis","Severe thrombocytosis"))),
      sprintf("JAK2 Allele Burden: %d%%\n", jak2),
      sprintf("Risk Category:      %s (IPSS-ET score = %d)\n", risk, pts),
      sprintf("Est. Thrombosis:    %s\n", thrombo),
      sprintf("Age >60:            %s\n", ifelse(input$age_gt60,"Yes","No")),
      sprintf("Prior Thrombosis:   %s\n", ifelse(input$prior_thr,"Yes (High-risk)","No")),
      sprintf("ASXL1 mutation:     %s\n", ifelse(input$ASXL1_pos,"Yes (↑MF risk)","No"))
    )
  })

  ## PK plot
  output$pk_plot <- renderPlotly({
    t_hours <- seq(0, 24, by = 0.5)
    pk_profiles <- data.frame(time = t_hours)

    ## Simple 2-cmpt PK profiles
    biexp <- function(t, D, F, V1, CL, Q, V2, ka) {
      alpha <- ((V1+V2)*CL + Q*(V1+V2)) / (2*V1*V2)
      beta_  <- CL*Q / (V1*V2*alpha)
      D*F*ka/(V1*(ka-alpha)) * (exp(-alpha*t) - exp(-ka*t)) +
      D*F*ka/(V1*(ka-beta_)) * (exp(-beta_*t) - exp(-ka*t))
    }

    pk_profiles$HU  <- if(input$dose_hu>0)
      biexp(t_hours, input$dose_hu/2/1e3, 0.80, 28, 4.20, 1.80, 12, 1.40) else 0
    pk_profiles$ANA <- if(input$dose_ana>0)
      biexp(t_hours, input$dose_ana/2/1e3, 0.70, 18, 9.50, 2.10, 8, 4.60)*1000 else 0
    pk_profiles$RUX <- if(input$dose_rux>0)
      biexp(t_hours, input$dose_rux/2/1e3, 0.95, 72, 22.0, 8.50, 38, 2.30)*1000 else 0

    fig <- plot_ly() %>%
      add_lines(data=pk_profiles, x=~time, y=~HU,  name="Hydroxyurea (µg/mL)", line=list(color="#1565c0")) %>%
      add_lines(data=pk_profiles, x=~time, y=~ANA, name="Anagrelide (ng/mL)",   line=list(color="#6a1b9a")) %>%
      add_lines(data=pk_profiles, x=~time, y=~RUX, name="Ruxolitinib (ng/mL)",  line=list(color="#00695c")) %>%
      layout(title="Drug PK Profiles (Single-dose, Day 1)",
             xaxis=list(title="Time (hours)"), yaxis=list(title="Concentration"),
             legend=list(orientation="h"))
    fig
  })

  ## Platelet dynamics
  output$plt_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~time, y=~PLT, type='scatter', mode='lines',
            line=list(color='#1565c0', width=2.5), name="Platelets") %>%
      add_lines(x=c(0,max(df$time)), y=c(400,400), line=list(dash='dash',color='darkgreen'), name="CHR (400)") %>%
      add_lines(x=c(0,max(df$time)), y=c(600,600), line=list(dash='dash',color='orange'),    name="PHR (600)") %>%
      add_lines(x=c(0,max(df$time)), y=c(1500,1500),line=list(dash='dash',color='red'),      name="Hemorrhage risk (1500)") %>%
      layout(title="Platelet Count (×10⁹/L)",
             xaxis=list(title="Day"), yaxis=list(title="Platelets (×10⁹/L)"),
             legend=list(orientation="h"))
  })

  output$plt_response_table <- renderTable({
    df <- sim_data()
    timepoints <- c(90, 180, 365, 730)
    valid_tp   <- timepoints[timepoints <= input$sim_days]
    purrr::map_dfr(valid_tp, function(tp) {
      row <- df[which.min(abs(df$time - tp)), ]
      data.frame(Day=tp, PLT=round(row$PLT), CHR=row$CHR, PHR=row$PHR)
    })
  })

  ## JAK2 dynamics
  output$jak2_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~time, y=~JAK2_AB, type='scatter', mode='lines',
            line=list(color='#e65100', width=2.5), name="JAK2 AB") %>%
      add_lines(x=c(0,max(df$time)), y=c(1,1), line=list(dash='dash',color='purple'), name="CMR (<1%)") %>%
      add_lines(x=c(0,max(df$time)), y=c(input$JAK2_init*0.5, input$JAK2_init*0.5),
                line=list(dash='dot',color='steelblue'), name="PMR (50% reduction)") %>%
      layout(title="JAK2 V617F Allele Burden (%)",
             xaxis=list(title="Day"), yaxis=list(title="JAK2 AB (%)", range=c(0,100)),
             legend=list(orientation="h"))
  })

  output$mol_response_table <- renderTable({
    df <- sim_data()
    tp <- c(90,180,365,730)
    purrr::map_dfr(tp[tp<=input$sim_days], function(t) {
      row <- df[which.min(abs(df$time-t)),]
      resp <- if(row$JAK2_AB < 1) "CMR" else if(row$JAK2_AB < input$JAK2_init*0.5) "PMR" else "NMR"
      data.frame(Day=t, `JAK2 AB(%)`=round(row$JAK2_AB,1), Response=resp)
    })
  }, colnames=TRUE)

  ## Thrombosis risk
  output$risk_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~time, y=~(Hazard_T*365*100), type='scatter', mode='lines',
            line=list(color='#c62828', width=2.5)) %>%
      layout(title="Annual Thrombosis Hazard (%/year)",
             xaxis=list(title="Day"), yaxis=list(title="Thrombosis Hazard (%/yr)"))
  })

  output$cum_risk_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~time, y=~(RISK_T*100), type='scatter', mode='lines',
            fill='tozeroy', fillcolor='rgba(198,40,40,0.2)',
            line=list(color='#c62828')) %>%
      layout(title="Cumulative Thrombosis Risk (%)",
             xaxis=list(title="Day"), yaxis=list(title="Cumulative Risk (%)"))
  })

  output$risk_factors_table <- DT::renderDT({
    data.frame(
      `Risk Factor` = c("Age >60 years","Prior thrombosis","JAK2 V617F vs CALR",
                        "Platelet >1000","Leukocytosis (WBC>11)","Cardiovascular risk"),
      `OR/HR` = c("2.0","3.2","1.8","1.5","1.8","1.4"),
      `Impact` = c("2× thrombosis risk","3× arterial events","More arterial events",
                   "Acquired vWD → hemorrhage","↑thrombosis","Synergistic"),
      `Reference` = c("Barbui 2012","Carobbio 2008","Rumi 2014 Blood",
                      "Guglielmelli 2021","Barbui 2013","Falchi 2017")
    )
  }, options=list(pageLength=6), rownames=FALSE)

  ## Scenario comparison
  output$compare_plot <- renderPlotly({
    all_df <- all_scenarios_data()
    sel_ids <- input$sel_scenarios
    yvar    <- input$compare_var
    ylabel  <- switch(yvar, PLT="Platelets (×10⁹/L)", JAK2_AB="JAK2 AB (%)",
                            Hazard_T="Thrombosis Hazard/day", Spleen="Spleen (cm)")

    filtered <- all_df %>% filter(id %in% sel_ids)
    pal <- setNames(sapply(scenario_defs, `[[`, "color"),
                    sapply(scenario_defs, `[[`, "id"))

    fig <- plot_ly()
    for (sc_id in sel_ids) {
      sc_df <- filtered %>% filter(id == sc_id)
      sc_label <- scenario_defs[[which(sapply(scenario_defs,`[[`,"id")==sc_id)]]$label
      fig <- fig %>% add_lines(data=sc_df, x=~time, y=as.formula(paste0("~",yvar)),
                               name=sc_label, line=list(color=pal[sc_id], width=2))
    }
    fig %>% layout(xaxis=list(title="Day"), yaxis=list(title=ylabel),
                   legend=list(orientation="h"))
  })

  output$compare_table <- DT::renderDT({
    all_df <- all_scenarios_data()
    tp <- c(90,180,365,730)
    purrr::map_dfr(scenario_defs, function(sc) {
      df <- all_df %>% filter(id==sc$id)
      purrr::map_dfr(tp[tp<=input$sim_days], function(t) {
        r <- df[which.min(abs(df$time-t)),]
        data.frame(Scenario=sc$label, Day=t, PLT=round(r$PLT),
                   JAK2_AB=round(r$JAK2_AB,1),
                   CHR=ifelse(r$CHR==1,"Yes","No"),
                   Spleen=round(r$Spleen,2),
                   Thromb_Haz=round(r$Hazard_T*365*100,2))
      })
    })
  }, options=list(pageLength=10, scrollX=TRUE), rownames=FALSE)

  ## Biomarker panel
  output$bio_panel <- renderPlotly({
    df <- sim_data()
    p1 <- plot_ly(df, x=~time, y=~PLT,     name="PLT (×10⁹/L)",  type='scatter', mode='lines', line=list(color='#1565c0'))
    p2 <- plot_ly(df, x=~time, y=~JAK2_AB, name="JAK2 AB (%)",   type='scatter', mode='lines', line=list(color='#e65100'))
    p3 <- plot_ly(df, x=~time, y=~Spleen,  name="Spleen (cm)",    type='scatter', mode='lines', line=list(color='#388e3c'))
    p4 <- plot_ly(df, x=~time, y=~TPO_level,name="TPO (pg/mL)",  type='scatter', mode='lines', line=list(color='#6a1b9a'))
    p5 <- plot_ly(df, x=~time, y=~MK_pool, name="MK Pool (AU)",   type='scatter', mode='lines', line=list(color='#00695c'))
    p6 <- plot_ly(df, x=~time, y=~(RISK_MF*100), name="Cum MF Risk (%)", type='scatter', mode='lines', line=list(color='#c62828'))

    subplot(list(p1,p2,p3,p4,p5,p6), nrows=3, shareX=TRUE,
            titleX=TRUE, titleY=TRUE) %>%
      layout(title="Multi-Biomarker Panel",
             showlegend=FALSE,
             annotations=list(
               list(x=0.22,y=1.02,text="Platelets",xref="paper",yref="paper",showarrow=FALSE),
               list(x=0.78,y=1.02,text="JAK2 AB",xref="paper",yref="paper",showarrow=FALSE),
               list(x=0.22,y=0.65,text="Spleen",xref="paper",yref="paper",showarrow=FALSE),
               list(x=0.78,y=0.65,text="TPO",xref="paper",yref="paper",showarrow=FALSE),
               list(x=0.22,y=0.30,text="MK Pool",xref="paper",yref="paper",showarrow=FALSE),
               list(x=0.78,y=0.30,text="MF Risk",xref="paper",yref="paper",showarrow=FALSE)
             ))
  })

  output$bio_ref_table <- DT::renderDT({
    data.frame(
      Biomarker = c("Platelet Count","JAK2 V617F AB","WBC","Hemoglobin","LDH","Spleen Size","Ferritin"),
      `Normal Range` = c("150–400×10⁹/L","0–1%","4–11×10⁹/L","12–16 g/dL (F)","<250 U/L","<13 cm","12–300 ng/mL"),
      `ET Typical` = c("500–1500×10⁹/L","40–70%","Normal/↑","Normal","Normal/↑","13–20 cm","Variable"),
      `Clinical Significance` = c("Primary criterion","Prognosis, thrombosis risk","Leukemic phase indicator",
                                   "Concurrent PV feature","BM turnover","EMH indicator","Iron status")
    )
  }, options=list(pageLength=7), rownames=FALSE)

  ## BM progression
  output$progression_plot <- renderPlotly({
    df <- sim_data()
    plot_ly() %>%
      add_lines(data=df, x=~time, y=~(RISK_MF*100), name="Cum MF Risk (%)",
                line=list(color='#880e4f', width=2.5)) %>%
      add_lines(data=df, x=~time, y=~JAK2_AB, name="JAK2 AB (%)",
                line=list(color='#e65100', width=2, dash='dot')) %>%
      layout(title="BM Progression Risk & JAK2 Dynamics",
             xaxis=list(title="Day"),
             yaxis=list(title="Value (%)"),
             legend=list(orientation="h"))
  })

  output$mf_risk_summary <- renderTable({
    df <- sim_data()
    last <- df[nrow(df),]
    data.frame(
      Metric    = c("Simulated Duration","Final JAK2 AB","Cum MF Risk","Cum Thrombo Risk"),
      Value     = c(paste0(input$sim_days," days"),
                    sprintf("%.1f%%", last$JAK2_AB),
                    sprintf("%.2f%%", last$RISK_MF*100),
                    sprintf("%.2f%%", last$RISK_T*100))
    )
  })
}

shinyApp(ui, server)
