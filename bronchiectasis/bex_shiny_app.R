## ============================================================
## Bronchiectasis QSP Interactive Shiny Dashboard
## ============================================================
## Tabs:
##   1. Patient Profile & Disease Severity
##   2. Drug PK — Concentration-Time Profiles
##   3. PD Biomarkers (IL-8, NE, Bacterial Load)
##   4. Lung Function & Structural Outcomes
##   5. Treatment Scenario Comparison
##   6. Exacerbation & Long-term Prognosis
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)
library(DT)

## ---- mrgsolve model (inline) ----
bex_code <- '
$PARAM
Bmax=1e9, B0=1e7, kgrow=0.8, kbiofilm=0.15, BF0=0.3, kBF_clear=0.05,
Nbase=2.5e6, kN_recruit=0.6, kN_death=0.3, kN_bact=1e-9,
IL8_base=5.0, kIL8_prod=2.0, kIL8_deg=0.5,
NE_base=0.8, kNE_prod=1.5e-7, kNE_deg=0.4, A1AT_base=100,
MCC0=1.0, kMCC_NE=0.003, kMCC_recov=0.05, MCC_min=0.1,
AD0=0.3, kAD_NE=0.002, kAD_repair=0.01, AD_max=0.95,
FEV1_0=1800, kFEV1_AD=500, FEV1_min=400,
Ex_threshold=0.7, kEx_resolve=0.25,
Ka_AZM=0.5, F_AZM=0.37, Vd_AZM=31, CL_AZM=22, Vp_AZM=800, Q_AZM=5,
EC50_AZM_IL8=0.05, Emax_AZM=0.50, EC50_AZM_QS=0.02, Emax_AZM_QS=0.35,
F_Tobra_inh=0.10, Vd_Tobra=0.26, CL_Tobra=5, ke_Tobra_lung=0.8,
MIC_PA_Tobra=4, EC50_Tobra=16, Emax_Tobra=0.85, Hill_Tobra=1.5,
Ka_Cipro=1.2, F_Cipro=0.70, Vd_Cipro=2.5, CL_Cipro=25,
ELF_ratio=0.60, MIC_PA_Cipro=0.5, Emax_Cipro=0.80, EC50_Cipro=1.0,
Kd_DNase=0.3, Emax_DNase=0.60, EC50_DNase=0.5,
AZM_on=0, TIP_on=0, Cipro_on=0, DNase_on=0,
wt=70

$CMT GUT_AZM CENT_AZM PERI_AZM LUNG_Tobra CENT_Tobra GUT_Cipro CENT_Cipro LUNG_DNase
$CMT BACT BIOFILM NEUT IL8 NE MCC AD EXAC

$MAIN
double C_AZM_lung    = (CENT_AZM/(Vd_AZM*wt))*10;
double C_Tobra_lung  = LUNG_Tobra/60;
double C_Cipro_ELF   = (CENT_Cipro/(Vd_Cipro*wt))*ELF_ratio;
double C_DNase_lung  = LUNG_DNase/30;
BACT_0=0.7; BIOFILM_0=0.40; NEUT_0=5.0; IL8_0=20; NE_0=2.5;
MCC_0=0.45; AD_0=0.50; EXAC_0=0;

$ODE
double C_AZM_p = CENT_AZM/(Vd_AZM*wt);
double C_AZM_l = C_AZM_p*10;
double C_Tl = LUNG_Tobra/60;
double C_Ce = (CENT_Cipro/(Vd_Cipro*wt))*ELF_ratio;
double C_Dl = LUNG_DNase/30;

double E_AZM_IL8 = Emax_AZM*C_AZM_l/(EC50_AZM_IL8+C_AZM_l);
double E_AZM_QS  = Emax_AZM_QS*C_AZM_l/(EC50_AZM_QS+C_AZM_l);
double E_Tobra   = Emax_Tobra*pow(C_Tl,Hill_Tobra)/(pow(EC50_Tobra,Hill_Tobra)+pow(C_Tl,Hill_Tobra));
double E_Cipro   = Emax_Cipro*C_Ce/(EC50_Cipro+C_Ce);
double E_kill    = 1-(1-E_Tobra*TIP_on)*(1-E_Cipro*Cipro_on);
if(E_kill>0.95) E_kill=0.95;
double E_Dn      = Emax_DNase*C_Dl/(EC50_DNase+C_Dl);

dxdt_GUT_AZM  = -Ka_AZM*GUT_AZM;
dxdt_CENT_AZM = Ka_AZM*GUT_AZM*F_AZM -(CL_AZM/Vd_AZM)*CENT_AZM -(Q_AZM/Vd_AZM)*CENT_AZM +(Q_AZM/Vp_AZM)*PERI_AZM;
dxdt_PERI_AZM = (Q_AZM/Vd_AZM)*CENT_AZM -(Q_AZM/Vp_AZM)*PERI_AZM;
dxdt_LUNG_Tobra= -ke_Tobra_lung*LUNG_Tobra-(F_Tobra_inh*ke_Tobra_lung)*LUNG_Tobra;
dxdt_CENT_Tobra= F_Tobra_inh*ke_Tobra_lung*LUNG_Tobra-(CL_Tobra/(Vd_Tobra*wt))*CENT_Tobra;
dxdt_GUT_Cipro = -Ka_Cipro*GUT_Cipro;
dxdt_CENT_Cipro= Ka_Cipro*GUT_Cipro*F_Cipro-(CL_Cipro/(Vd_Cipro*wt))*CENT_Cipro;
dxdt_LUNG_DNase= -Kd_DNase*LUNG_DNase;

double Bn = BACT/(Bmax/1e7);
double kNk = kN_bact*NEUT*1e6;
dxdt_BACT = kgrow*BACT*(1-Bn) - kNk*BACT - MCC*0.4*BACT - E_kill*2.0*BACT;
dxdt_BIOFILM = kbiofilm*Bn*(1-BIOFILM)*(1-E_AZM_QS*AZM_on) - kBF_clear*BIOFILM*MCC;
double IS = IL8/(IL8+5);
dxdt_NEUT = kN_recruit*IS*(Nbase/1e6) - kN_death*NEUT + 0.1*(Nbase/1e6-NEUT);
double IL8p = kIL8_prod*Bn*(1-E_AZM_IL8*AZM_on);
dxdt_IL8 = IL8p+0.3*IL8_base-kIL8_deg*(IL8-IL8_base);
double Nf = A1AT_base/(A1AT_base+50);
dxdt_NE = kNE_prod*NEUT*1e6*(1-Nf)-kNE_deg*(NE-NE_base);
double MCCt = MCC0-kMCC_NE*NE+E_Dn*DNase_on;
if(MCCt<MCC_min) MCCt=MCC_min; if(MCCt>1.0) MCCt=1.0;
dxdt_MCC = kMCC_recov*(MCCt-MCC);
dxdt_AD = kAD_NE*NE*(1-AD/AD_max)-kAD_repair*(1-AD);
double Exd = (Bn>Ex_threshold)?(Bn-Ex_threshold)*5:0;
dxdt_EXAC = Exd*(1-EXAC)-kEx_resolve*EXAC;

$TABLE
double FEV1 = FEV1_0 - kFEV1_AD*AD;
if(FEV1<FEV1_min) FEV1=FEV1_min;
double FEV1_pct = FEV1/2400*100;
double log10_BACT = log10(BACT*1e7+1);
capture FEV1 FEV1_pct log10_BACT IL8 NE MCC AD BIOFILM EXAC
capture C_AZM_lung C_Tobra_lung C_Cipro_ELF C_DNase_lung
'

mod <- suppressMessages(mcode("bex_shiny", bex_code))

## ---- UI ----
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Bronchiectasis QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",       tabName = "tab_patient",  icon = icon("user")),
      menuItem("Drug PK Profiles",      tabName = "tab_pk",       icon = icon("pills")),
      menuItem("PD Biomarkers",         tabName = "tab_pd",       icon = icon("flask")),
      menuItem("Lung Function",         tabName = "tab_lung",     icon = icon("lungs")),
      menuItem("Scenario Comparison",   tabName = "tab_scenario", icon = icon("chart-bar")),
      menuItem("Exacerbation & Prognosis", tabName = "tab_prog",  icon = icon("heartbeat"))
    ),

    hr(),
    h5("  Patient Settings", style = "color:white; margin-left:10px;"),
    sliderInput("wt",    "Body weight (kg)", 40, 120, 70),
    sliderInput("fev1_0","Baseline FEV1 (mL)", 600, 3500, 1800),
    selectInput("severity", "BSI Severity",
                choices = c("Mild (BSI 0-4)" = "mild",
                            "Moderate (BSI 5-8)" = "mod",
                            "Severe (BSI ≥9)" = "sev"),
                selected = "mod"),

    hr(),
    h5("  Treatment", style = "color:white; margin-left:10px;"),
    checkboxInput("use_azm",   "Azithromycin (250mg 3×/wk)", FALSE),
    checkboxInput("use_tip",   "Inhaled Tobramycin (300mg BID)", FALSE),
    checkboxInput("use_cipro", "Ciprofloxacin 500mg BID (acute)", FALSE),
    checkboxInput("use_dnase", "Dornase Alfa (2.5mg daily)", FALSE),

    hr(),
    sliderInput("sim_days", "Simulation duration (days)", 30, 730, 365)
  ),

  dashboardBody(
    tabItems(

      ## ------ TAB 1: Patient Profile ------
      tabItem(tabName = "tab_patient",
        h2("Patient Profile & Baseline Characteristics"),
        fluidRow(
          valueBoxOutput("box_severity"),
          valueBoxOutput("box_fev1"),
          valueBoxOutput("box_exac_risk")
        ),
        fluidRow(
          box(title = "Baseline Disease State", status = "primary", solidHeader = TRUE,
              width = 6, DTOutput("tbl_baseline")),
          box(title = "BSI Score Components", status = "warning", solidHeader = TRUE,
              width = 6,
              plotlyOutput("plot_bsi", height = "320px"))
        ),
        fluidRow(
          box(title = "Vicious Cycle Summary", status = "danger", solidHeader = TRUE,
              width = 12,
              HTML('<div style="background:#fff3e0;padding:12px;border-radius:6px;">
                <b>The Bronchiectasis Vicious Cycle:</b><br>
                <span style="color:#c62828;">⟳ Infection (PA/HI colonization)</span> →
                <span style="color:#e53935;">Neutrophilic Inflammation (↑IL-8, ↑NE)</span> →
                <span style="color:#283593;">Protease Imbalance (NE/MMP ≫ α1-AT/SLPI)</span> →
                <span style="color="#0d47a1;">Structural Airway Damage</span> →
                <span style="color:#1b5e20;">MCC Failure</span> →
                <span style="color:#c62828;">Re-infection ⟲</span><br><br>
                Key: Interrupting any step (antibiotics → infection, azithromycin → inflammation,
                mucolytics/physiotherapy → MCC) can slow disease progression.
              </div>'))
        )
      ),

      ## ------ TAB 2: Drug PK ------
      tabItem(tabName = "tab_pk",
        h2("Drug PK — Concentration-Time Profiles"),
        fluidRow(
          box(title = "Select Drug for PK Display", status = "primary", solidHeader = TRUE,
              width = 3,
              selectInput("pk_drug", "Drug", choices = c(
                "Azithromycin (plasma)"  = "C_AZM_lung",
                "Inhaled Tobramycin (lung)" = "C_Tobra_lung",
                "Ciprofloxacin (ELF)"    = "C_Cipro_ELF",
                "Dornase Alfa (lung)"    = "C_DNase_lung"
              )),
              numericInput("pk_dose_azm",   "AZM dose (mg)",         250, min = 125, max = 1000),
              numericInput("pk_dose_tip",   "Tobramycin dose (mg)",  300, min = 100, max = 600),
              numericInput("pk_dose_cipro", "Ciprofloxacin dose (mg)",500, min = 250, max = 750),
              numericInput("pk_dose_dn",    "DNase dose (mg)",        2.5, min = 1, max = 5)
          ),
          box(title = "Concentration-Time Profile", status = "primary", solidHeader = TRUE,
              width = 9, plotlyOutput("plot_pk", height = "420px"))
        ),
        fluidRow(
          box(title = "PK Parameter Summary Table", status = "info", solidHeader = TRUE,
              width = 12, DTOutput("tbl_pk_params"))
        )
      ),

      ## ------ TAB 3: PD Biomarkers ------
      tabItem(tabName = "tab_pd",
        h2("Pharmacodynamic Biomarkers"),
        fluidRow(
          box(title = "Sputum IL-8 (ng/mL)", status = "danger", solidHeader = TRUE,
              width = 6, plotlyOutput("plot_il8", height = "300px")),
          box(title = "Neutrophil Elastase Activity (µg/mL)", status = "danger", solidHeader = TRUE,
              width = 6, plotlyOutput("plot_ne", height = "300px"))
        ),
        fluidRow(
          box(title = "Bacterial Load (log10 CFU/mL)", status = "warning", solidHeader = TRUE,
              width = 6, plotlyOutput("plot_bact", height = "300px")),
          box(title = "Biofilm Fraction (0-1)", status = "warning", solidHeader = TRUE,
              width = 6, plotlyOutput("plot_biofilm", height = "300px"))
        ),
        fluidRow(
          box(title = "Mucociliary Clearance Index (0-1)", status = "success", solidHeader = TRUE,
              width = 6, plotlyOutput("plot_mcc", height = "300px")),
          box(title = "Sputum Neutrophil Count (×10⁶/mL)", status = "info", solidHeader = TRUE,
              width = 6, plotlyOutput("plot_neut", height = "300px"))
        )
      ),

      ## ------ TAB 4: Lung Function ------
      tabItem(tabName = "tab_lung",
        h2("Lung Function & Structural Outcomes"),
        fluidRow(
          box(title = "FEV1 (% predicted) Over Time", status = "primary", solidHeader = TRUE,
              width = 8, plotlyOutput("plot_fev1", height = "380px")),
          box(title = "FEV1 Classification", status = "info", solidHeader = TRUE,
              width = 4,
              HTML('<table style="width:100%;font-size:12px;">
                <tr bgcolor="#FFCDD2"><td><b>Stage</b></td><td><b>FEV1 % pred</b></td><td><b>BSI Category</b></td></tr>
                <tr><td>Mild</td><td>≥70%</td><td>BSI 0-4</td></tr>
                <tr bgcolor="#FFF9C4"><td>Moderate</td><td>50-69%</td><td>BSI 5-8</td></tr>
                <tr bgcolor="#FFCDD2"><td>Severe</td><td>30-49%</td><td>BSI 9-14</td></tr>
                <tr bgcolor="#EF9A9A"><td>Very Severe</td><td>&lt;30%</td><td>BSI ≥15</td></tr>
              </table>'),
              br(),
              plotlyOutput("plot_fev1_gauge", height = "220px"))
        ),
        fluidRow(
          box(title = "Airway Damage Score & MCC Index", status = "warning", solidHeader = TRUE,
              width = 12, plotlyOutput("plot_structural", height = "320px"))
        )
      ),

      ## ------ TAB 5: Scenario Comparison ------
      tabItem(tabName = "tab_scenario",
        h2("Treatment Scenario Comparison"),
        fluidRow(
          box(title = "Scenarios to Compare", status = "primary", solidHeader = TRUE, width = 3,
              checkboxGroupInput("scen_select", "Select Scenarios:",
                choices = c(
                  "No Treatment"                  = "notreat",
                  "Azithromycin"                  = "azm",
                  "Inhaled Tobramycin (cycled)"   = "tip",
                  "AZM + Inhaled Tobramycin"      = "azm_tip",
                  "AZM + TIP + Dornase Alfa"      = "azm_tip_dn"
                ),
                selected = c("notreat","azm","azm_tip")
              ),
              selectInput("scen_endpoint", "Primary Endpoint",
                choices = c("FEV1 (% predicted)" = "FEV1_pct",
                            "Bacterial Load (log10)" = "log10_BACT",
                            "IL-8 (ng/mL)" = "IL8",
                            "NE Activity (µg/mL)" = "NE",
                            "MCC Index" = "MCC",
                            "Airway Damage" = "AD",
                            "Exacerbation State" = "EXAC")),
              actionButton("run_scen", "Run Comparison", class = "btn-primary", width = "100%")
          ),
          box(title = "Scenario Comparison Plot", status = "primary", solidHeader = TRUE,
              width = 9, plotlyOutput("plot_scen", height = "440px"))
        ),
        fluidRow(
          box(title = "Numerical Comparison Table", status = "info", solidHeader = TRUE,
              width = 12, DTOutput("tbl_scen"))
        )
      ),

      ## ------ TAB 6: Exacerbation & Prognosis ------
      tabItem(tabName = "tab_prog",
        h2("Exacerbation Dynamics & Long-term Prognosis"),
        fluidRow(
          box(title = "Exacerbation State Over Time", status = "danger", solidHeader = TRUE,
              width = 8, plotlyOutput("plot_exac", height = "360px")),
          box(title = "Exacerbation Risk Factors", status = "warning", solidHeader = TRUE,
              width = 4,
              HTML('<b>High-risk features (BSI ≥9):</b>
              <ul>
                <li>FEV1 &lt;50% predicted</li>
                <li>≥3 exacerbations/year</li>
                <li>PA colonization</li>
                <li>Radiological severity</li>
                <li>MRC Dyspnoea ≥3</li>
              </ul>
              <b>Annual mortality:</b><br>
              BSI mild (0-4): 0-2%/yr<br>
              BSI moderate (5-8): 5-10%/yr<br>
              BSI severe (≥9): 15-20%/yr<br>
              <hr>
              <b>EMBRACE trial (AZM 250mg 3×/wk):</b><br>
              Exacerbation rate: 0.59 vs 1.57 (RR 0.38)<br>
              Sputum volume ↓47%, FEV1 ↑5.4%<br>
              <i>Lancet 2012;380:660-667</i>'))
        ),
        fluidRow(
          box(title = "5-Year FEV1 Trajectory Projection", status = "primary", solidHeader = TRUE,
              width = 6, plotlyOutput("plot_longterm", height = "320px")),
          box(title = "Exacerbation Frequency Summary", status = "info", solidHeader = TRUE,
              width = 6, DTOutput("tbl_exac_summary"))
        )
      )
    )
  )
)

## ---- SERVER ----
server <- function(input, output, session) {

  ## Helper: severity-based initial conditions
  get_init <- reactive({
    switch(input$severity,
      mild = init(mod, BACT=0.2, BIOFILM=0.15, NEUT=2.5, IL8=8, NE=0.9, MCC=0.75, AD=0.20, EXAC=0),
      mod  = init(mod, BACT=0.7, BIOFILM=0.40, NEUT=5.0, IL8=20, NE=2.5, MCC=0.45, AD=0.50, EXAC=0),
      sev  = init(mod, BACT=1.5, BIOFILM=0.65, NEUT=10.0, IL8=55, NE=5.5, MCC=0.20, AD=0.75, EXAC=0)
    )
  })

  ## Helper: make dosing events
  make_ev <- function(azm = 0, tip = 0, cipro = 0, dnase = 0, days = 365) {
    evs <- list()
    if(azm > 0) {
      t_a <- seq(0, days, by = 7/3)
      evs[["azm"]] <- ev(amt = azm, cmt = "GUT_AZM", time = t_a)
    }
    if(tip > 0) {
      on_d <- unlist(lapply(seq(0, days-56, by=56), function(d) seq(d, d+27)))
      t_t <- sort(unique(c(on_d, on_d+0.5))); t_t <- t_t[t_t<=days]
      evs[["tip"]] <- ev(amt = tip, cmt = "LUNG_Tobra", time = t_t)
    }
    if(cipro > 0) {
      t_c <- seq(0, days, by = 0.5)
      evs[["cipro"]] <- ev(amt = cipro, cmt = "GUT_Cipro", time = t_c)
    }
    if(dnase > 0) {
      evs[["dnase"]] <- ev(amt = dnase, cmt = "LUNG_DNase", time = seq(0, days, by = 1))
    }
    if(length(evs) == 0) return(ev(amt=0, cmt=1, time=Inf))
    do.call(c, evs)
  }

  ## Reactive simulation
  sim_data <- reactive({
    p <- param(mod,
      AZM_on   = as.integer(input$use_azm),
      TIP_on   = as.integer(input$use_tip),
      Cipro_on = as.integer(input$use_cipro),
      DNase_on = as.integer(input$use_dnase),
      FEV1_0   = input$fev1_0,
      wt       = input$wt
    )
    ev_all <- make_ev(
      azm   = if(input$use_azm)   250 else 0,
      tip   = if(input$use_tip)   300 else 0,
      cipro = if(input$use_cipro) 500 else 0,
      dnase = if(input$use_dnase) 2.5 else 0,
      days  = input$sim_days
    )
    mrgsim(p, init = get_init(), events = ev_all,
           tgrid = seq(0, input$sim_days, by = 0.5)) %>%
      as_tibble()
  })

  ## TAB 1 outputs
  output$box_severity <- renderValueBox({
    sev_lab <- switch(input$severity,
                      mild = "Mild (BSI 0-4)", mod = "Moderate (BSI 5-8)", sev = "Severe (BSI ≥9)")
    valueBox(sev_lab, "Disease Severity", icon = icon("lungs"), color = "blue")
  })
  output$box_fev1 <- renderValueBox({
    pct <- round(input$fev1_0 / 2400 * 100)
    valueBox(paste0(pct, "% predicted"), paste0("FEV1 (", input$fev1_0, " mL)"),
             icon = icon("chart-line"), color = if(pct >= 70) "green" else if(pct >= 50) "yellow" else "red")
  })
  output$box_exac_risk <- renderValueBox({
    risk <- switch(input$severity, mild = "Low (~1/yr)", mod = "Moderate (2-3/yr)", sev = "High (≥3/yr)")
    valueBox(risk, "Annual Exacerbation Risk", icon = icon("exclamation-triangle"),
             color = switch(input$severity, mild="green", mod="yellow", sev="red"))
  })
  output$tbl_baseline <- renderDT({
    init_v <- switch(input$severity,
      mild = data.frame(Variable=c("Bacterial Load","Biofilm","IL-8","NE","MCC","Airway Damage"),
                        Value=c("2×10⁶ CFU/mL","15%","8 ng/mL","0.9 µg/mL","75%","20%")),
      mod  = data.frame(Variable=c("Bacterial Load","Biofilm","IL-8","NE","MCC","Airway Damage"),
                        Value=c("7×10⁶ CFU/mL","40%","20 ng/mL","2.5 µg/mL","45%","50%")),
      sev  = data.frame(Variable=c("Bacterial Load","Biofilm","IL-8","NE","MCC","Airway Damage"),
                        Value=c("1.5×10⁷ CFU/mL","65%","55 ng/mL","5.5 µg/mL","20%","75%"))
    )
    datatable(init_v, options=list(dom='t', pageLength=8), rownames=FALSE)
  })
  output$plot_bsi <- renderPlotly({
    bsi_df <- data.frame(
      Component = c("FEV1","Exacerbations","Hospital","Colonization","Radiology","MRC Dyspnoea"),
      Score = switch(input$severity,
        mild = c(0,1,0,0,1,1),
        mod  = c(1,2,0,1,2,1),
        sev  = c(3,4,1,2,3,2)
      )
    )
    plot_ly(bsi_df, x = ~Component, y = ~Score, type = "bar",
            marker = list(color = "#1E88E5")) %>%
      layout(title = "BSI Component Scores",
             xaxis = list(title = ""), yaxis = list(title = "Score"),
             showlegend = FALSE)
  })

  ## TAB 2 PK
  output$plot_pk <- renderPlotly({
    pk_params <- param(mod, AZM_on=1, TIP_on=1, Cipro_on=1, DNase_on=1, wt=input$wt)
    t_pk <- seq(0, 96, by=0.25)
    ev_pk <- make_ev(azm=input$pk_dose_azm, tip=input$pk_dose_tip,
                     cipro=input$pk_dose_cipro, dnase=input$pk_dose_dn, days=4)
    pk_sim <- mrgsim(pk_params, events=ev_pk, tgrid=t_pk) %>% as_tibble()
    drug_var <- input$pk_drug
    y_label <- switch(drug_var,
      C_AZM_lung    = "AZM Lung Conc (mg/L)",
      C_Tobra_lung  = "Tobramycin Lung Conc (mg/L)",
      C_Cipro_ELF   = "Cipro ELF Conc (mg/L)",
      C_DNase_lung  = "DNase Lung Conc (mg/L)"
    )
    y_vals <- pk_sim[[drug_var]]
    plot_ly(x=pk_sim$time, y=y_vals, type="scatter", mode="lines",
            line=list(color="#1565C0", width=2)) %>%
      layout(title=y_label, xaxis=list(title="Time (hours)"),
             yaxis=list(title=y_label))
  })
  output$tbl_pk_params <- renderDT({
    pk_df <- data.frame(
      Drug = c("Azithromycin","Inhaled Tobramycin","Ciprofloxacin","Dornase Alfa"),
      Route = c("Oral","Inhaled","Oral","Inhaled"),
      Dose  = c("250mg 3×/wk","300mg BID","500mg BID","2.5mg daily"),
      Vd    = c("31 L/kg","0.26 L/kg","2.5 L/kg","Lung local"),
      t_half = c("68h","2-3h (lung)","4.2h","~2h (lung)"),
      F_pct = c("37%","~10% systemic","70%","N/A"),
      Lung_ELF = c("~10× plasma",">>10× plasma","0.6× plasma","Local"),
      Emax = c("50% IL-8↓","85% kill","80% kill","60% viscosity↓"),
      EC50  = c("0.05 mg/L","16 mg/L","1.0 mg/L (ELF)","0.5 mg/L")
    )
    datatable(pk_df, options=list(dom='t', scrollX=TRUE), rownames=FALSE)
  })

  ## TAB 3 PD Biomarkers
  output$plot_il8 <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~IL8, type="scatter", mode="lines",
            line=list(color="#E53935",width=2)) %>%
      add_lines(y=~IL8_base, name="Normal (5 ng/mL)",
                line=list(color="gray",dash="dot"), inherit=FALSE,
                x=d$time, y=rep(5, nrow(d))) %>%
      layout(title="IL-8 (ng/mL)", xaxis=list(title="Days"),
             yaxis=list(title="IL-8 (ng/mL)"), showlegend=TRUE)
  })
  output$plot_ne <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~NE, type="scatter", mode="lines",
            line=list(color="#C62828",width=2)) %>%
      layout(title="Neutrophil Elastase (µg/mL)",
             xaxis=list(title="Days"), yaxis=list(title="NE (µg/mL)"))
  })
  output$plot_bact <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~log10_BACT, type="scatter", mode="lines",
            line=list(color="#F57F17",width=2)) %>%
      layout(title="Bacterial Load (log10 CFU/mL)",
             xaxis=list(title="Days"), yaxis=list(title="log10[CFU/mL]"))
  })
  output$plot_biofilm <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~BIOFILM, type="scatter", mode="lines",
            line=list(color="#7B1FA2",width=2)) %>%
      layout(title="Biofilm Fraction (0-1)",
             xaxis=list(title="Days"), yaxis=list(title="Biofilm Fraction"))
  })
  output$plot_mcc <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~MCC, type="scatter", mode="lines",
            line=list(color="#2E7D32",width=2)) %>%
      add_lines(x=d$time, y=rep(1,nrow(d)), name="Normal MCC",
                line=list(color="gray",dash="dot")) %>%
      layout(title="MCC Index (1=normal, 0=no clearance)",
             xaxis=list(title="Days"), yaxis=list(title="MCC Index",range=c(0,1.1)))
  })
  output$plot_neut <- renderPlotly({
    d <- sim_data()
    d <- d %>% mutate(Neut_plot = NEUT * 1e6 / 1e6)  # back to ×10⁶
    plot_ly(d, x=~time, y=~Neut_plot, type="scatter", mode="lines",
            line=list(color="#1565C0",width=2)) %>%
      layout(title="Sputum Neutrophil Count (×10⁶/mL)",
             xaxis=list(title="Days"), yaxis=list(title="Neutrophils ×10⁶/mL"))
  })

  ## TAB 4 Lung Function
  output$plot_fev1 <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~FEV1_pct, type="scatter", mode="lines",
            line=list(color="#1565C0",width=2.5), name="FEV1 % pred") %>%
      add_lines(x=c(0,max(d$time)), y=c(70,70), name="Mild cutoff (70%)",
                line=list(color="green",dash="dash")) %>%
      add_lines(x=c(0,max(d$time)), y=c(50,50), name="Moderate cutoff (50%)",
                line=list(color="orange",dash="dash")) %>%
      add_lines(x=c(0,max(d$time)), y=c(30,30), name="Severe cutoff (30%)",
                line=list(color="red",dash="dash")) %>%
      layout(title="FEV1 (% predicted) Trend",
             xaxis=list(title="Days"), yaxis=list(title="FEV1 % predicted"))
  })
  output$plot_fev1_gauge <- renderPlotly({
    d <- sim_data()
    final_fev1 <- tail(d$FEV1_pct, 1)
    plot_ly(type="indicator", mode="gauge+number",
            value=round(final_fev1, 1),
            title=list(text="Final FEV1 (% pred)"),
            gauge=list(
              axis=list(range=list(0,100)),
              bar=list(color="#1E88E5"),
              steps=list(
                list(range=c(0,30),color="#FFCDD2"),
                list(range=c(30,50),color="#FFECB3"),
                list(range=c(50,70),color="#DCEDC8"),
                list(range=c(70,100),color="#C8E6C9")
              )
            ))
  })
  output$plot_structural <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~AD, type="scatter", mode="lines",
            name="Airway Damage", line=list(color="#B71C1C",width=2)) %>%
      add_lines(x=d$time, y=d$MCC, name="MCC Index",
                line=list(color="#1B5E20",width=2)) %>%
      layout(title="Airway Damage (0-1) & MCC Index (0-1) Over Time",
             xaxis=list(title="Days"), yaxis=list(title="Index (0-1)",range=c(0,1.1)))
  })

  ## TAB 5 Scenario Comparison
  scen_results <- eventReactive(input$run_scen, {
    scens <- input$scen_select
    scen_params_list <- list(
      notreat  = list(AZM_on=0, TIP_on=0, Cipro_on=0, DNase_on=0, azm=0, tip=0, dn=0),
      azm      = list(AZM_on=1, TIP_on=0, Cipro_on=0, DNase_on=0, azm=250, tip=0, dn=0),
      tip      = list(AZM_on=0, TIP_on=1, Cipro_on=0, DNase_on=0, azm=0, tip=300, dn=0),
      azm_tip  = list(AZM_on=1, TIP_on=1, Cipro_on=0, DNase_on=0, azm=250, tip=300, dn=0),
      azm_tip_dn=list(AZM_on=1, TIP_on=1, Cipro_on=0, DNase_on=1, azm=250, tip=300, dn=2.5)
    )
    scen_labels <- c(notreat="No Treatment", azm="Azithromycin",
                     tip="Inhaled Tobramycin", azm_tip="AZM + TIP",
                     azm_tip_dn="AZM + TIP + DNase")
    bind_rows(lapply(scens, function(s) {
      sp <- scen_params_list[[s]]
      p_s <- param(mod, AZM_on=sp$AZM_on, TIP_on=sp$TIP_on,
                   Cipro_on=sp$Cipro_on, DNase_on=sp$DNase_on,
                   FEV1_0=input$fev1_0, wt=input$wt)
      ev_s <- make_ev(azm=sp$azm, tip=sp$tip, dnase=sp$dn, days=input$sim_days)
      mrgsim(p_s, init=get_init(), events=ev_s,
             tgrid=seq(0,input$sim_days,by=1)) %>%
        as_tibble() %>%
        mutate(scenario=scen_labels[[s]])
    }))
  })
  output$plot_scen <- renderPlotly({
    req(scen_results())
    d <- scen_results()
    y_col <- input$scen_endpoint
    y_lab <- input$scen_endpoint
    fig <- plot_ly()
    for(s in unique(d$scenario)) {
      ds <- d %>% filter(scenario==s)
      fig <- fig %>% add_lines(data=ds, x=~time, y=as.formula(paste0("~",y_col)),
                               name=s, line=list(width=2))
    }
    fig %>% layout(title=paste("Scenario Comparison —", y_lab),
                   xaxis=list(title="Days"), yaxis=list(title=y_lab),
                   legend=list(orientation="h", y=-0.25))
  })
  output$tbl_scen <- renderDT({
    req(scen_results())
    d <- scen_results()
    d %>% filter(time %in% c(0, 90, 180, 365)) %>%
      group_by(scenario, time) %>%
      summarise(FEV1_pct=round(mean(FEV1_pct),1),
                log10_BACT=round(mean(log10_BACT),2),
                IL8=round(mean(IL8),1),
                NE=round(mean(NE),2),
                MCC=round(mean(MCC),3),
                AD=round(mean(AD),3), .groups="drop") %>%
      mutate(Day=time) %>% select(-time) %>%
      datatable(options=list(dom='t', scrollX=TRUE), rownames=FALSE)
  })

  ## TAB 6 Exacerbation & Prognosis
  output$plot_exac <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~EXAC, type="scatter", mode="lines",
            fill="tozeroy", fillcolor="rgba(229,57,53,0.2)",
            line=list(color="#E53935",width=2)) %>%
      add_lines(x=c(0,max(d$time)), y=c(0.5,0.5), name="Exacerbation threshold",
                line=list(color="gray",dash="dash")) %>%
      layout(title="Exacerbation State Over Time (>0.5 = Active)",
             xaxis=list(title="Days"),
             yaxis=list(title="Exacerbation Score",range=c(0,1.1)))
  })
  output$plot_longterm <- renderPlotly({
    p_lt <- param(mod, AZM_on=as.integer(input$use_azm),
                  TIP_on=as.integer(input$use_tip),
                  Cipro_on=0, DNase_on=as.integer(input$use_dnase),
                  FEV1_0=input$fev1_0, wt=input$wt)
    ev_lt <- make_ev(azm=if(input$use_azm)250 else 0,
                     tip=if(input$use_tip)300 else 0,
                     dnase=if(input$use_dnase)2.5 else 0,
                     days=365*5)
    lt_sim <- mrgsim(p_lt, init=get_init(), events=ev_lt,
                     tgrid=seq(0,365*5,by=7)) %>% as_tibble()
    plot_ly(lt_sim, x=~(time/365), y=~FEV1_pct, type="scatter",mode="lines",
            line=list(color="#1565C0",width=2)) %>%
      add_lines(x=lt_sim$time/365, y=rep(50,nrow(lt_sim)),
                name="Severe cutoff", line=list(color="red",dash="dash")) %>%
      layout(title="5-Year FEV1 Trajectory",
             xaxis=list(title="Years"),
             yaxis=list(title="FEV1 % predicted"))
  })
  output$tbl_exac_summary <- renderDT({
    d <- sim_data()
    n_exac <- sum(diff(as.integer(d$EXAC > 0.5)) == 1)
    exac_dur <- mean(rle(d$EXAC > 0.5)$lengths[rle(d$EXAC > 0.5)$values] * 0.5)
    data.frame(
      Metric = c("Total exacerbation episodes (simulated)",
                 "Mean exacerbation duration (days)",
                 "Min FEV1 (% pred) during simulation",
                 "Final FEV1 (% pred)",
                 "Final Airway Damage Score",
                 "Final MCC Index"),
      Value = c(
        as.character(n_exac),
        ifelse(is.nan(exac_dur), "0", round(exac_dur, 1)),
        round(min(d$FEV1_pct), 1),
        round(tail(d$FEV1_pct, 1), 1),
        round(tail(d$AD, 1), 3),
        round(tail(d$MCC, 1), 3)
      )
    ) %>% datatable(options=list(dom='t'), rownames=FALSE)
  })
}

## ---- Run App ----
shinyApp(ui = ui, server = server)
