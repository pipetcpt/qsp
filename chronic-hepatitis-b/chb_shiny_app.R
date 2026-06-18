## ============================================================
## Chronic Hepatitis B (CHB) — Interactive Shiny QSP Dashboard
## Author: Claude Code Routine (CCR)  |  Date: 2026-06-18
## Tabs: Patient Profile · PK · Viral/HBsAg Dynamics ·
##       Clinical Endpoints · Scenario Comparison · Biomarkers
## ============================================================

library(shiny)
library(shinydashboard)
library(plotly)
library(dplyr)
library(tidyr)
library(DT)
library(mrgsolve)
library(ggplot2)

## ─────────────────────────────────────────────────────────────
## Inline mrgsolve CHB model
## ─────────────────────────────────────────────────────────────

chb_code <- '
$PROB CHB QSP Shiny Model

$PARAM @annotated
BWT:65:Body weight (kg)
T0:1e7:Baseline target hepatocytes
I0:1e5:Baseline infected hepatocytes
V0:1e7:Baseline HBV DNA (IU/mL)
ccc0:10:Baseline cccDNA
Ag0:1000:Baseline HBsAg (IU/mL)
ALT0:80:Baseline ALT (IU/L)
FIB0:1.5:Baseline fibrosis score
CTL0:0.5:Baseline CTL
ExhSS:0.4:Baseline T cell exhaustion
s_T:1e5:Target cell production
d_T:0.01:Target cell death rate
beta:2e-10:Infection rate
delta_I:0.15:Infected cell clearance
p_V:50:Virion production
c_V:0.67:Viral clearance
lambda_ccc:0.002:cccDNA replenishment
delta_ccc:0.003:cccDNA loss rate
k_Ag:0.05:HBsAg production
c_Ag:0.01:HBsAg clearance
alpha_CTL:0.8:CTL killing
k_CTL_exp:0.05:CTL expansion
d_CTL:0.03:CTL contraction
k_exhaust:0.02:T cell exhaustion rate
r_exhaust:0.005:Exhaustion recovery
k_IFNinn:0.1:Innate IFN induction
d_IFNinn:0.15:Innate IFN decay
IFN_antiv:0.4:IFN antiviral effect
k_ALT_I:50:ALT per infected cell
k_ALT_ret:0.05:ALT return rate
k_fibrosis:0.0005:Fibrosis progression
k_fib_reg:0.0002:Fibrosis regression
k_HCC:0.002:HCC risk rate
k_HSC:0.01:HSC activation
d_HSC:0.03:HSC deactivation
ETV_F:0.37:ETV bioavailability
ETV_ka:0.90:ETV absorption (/h)
ETV_Vc:73:ETV central volume (L)
ETV_CL:10.5:ETV clearance (L/h)
ETV_Kp:5.0:ETV liver partition
ETV_kTP:0.15:ETV phosphorylation
ETV_kTP_d:0.008:ETV-TP elimination
ETV_IC50:0.004:ETV-TP IC50 (µM)
ETV_Imax:0.99:ETV max inhibition
TDF_F:0.25:TDF bioavailability
TDF_ka:0.50:TDF absorption (/h)
TDF_Vc:18:TFV volume (L)
TDF_CL:11.2:TFV clearance (L/h)
TFV_kDP:0.10:TFV phosphorylation
TFV_kDP_d:0.012:TFV-DP elimination
TFV_IC50:0.500:TFV-DP IC50 (µM)
TFV_Imax:0.99:TFV max inhibition
IFN_F:0.80:PegIFN bioavailability
IFN_ka:0.04:PegIFN absorption (/h)
IFN_Vc:8.0:PegIFN volume (L)
IFN_CL:0.055:PegIFN clearance (L/h)
IFN_EC50:0.08:PegIFN EC50
IFN_Emax:0.80:PegIFN Emax antiviral
IFN_immE:0.50:PegIFN immunomodulation
IFN_cccE:0.30:PegIFN anti-cccDNA
IAGM_dose:0:siRNA active
siRNA_Emax:0.90:siRNA HBsAg reduction
USE_ETV:0:Use Entecavir
USE_TDF:0:Use TDF
USE_PIFN:0:Use Peg-IFN

$CMT ETV_gut ETV_C ETV_TP TDF_gut TDF_C TFV_DP IFN_SC IFN_C T_cell I_cell V_dna ccc_DNA HBsAg_C CTL_resp Exhaust IFN_inn ALT_val Fibrosis HCC_risk HSC_act

$MAIN
double ETV_Cp = ETV_C / ETV_Vc;
double TFV_Cp = TDF_C / TDF_Vc;
double IFN_Cp = IFN_C / IFN_Vc;
double ETV_eff = USE_ETV * ETV_Imax * ETV_TP / (ETV_TP + ETV_IC50);
double TFV_eff = USE_TDF * TFV_Imax * TFV_DP / (TFV_DP + TFV_IC50);
double NUC_eff = 1.0 - (1.0-ETV_eff)*(1.0-TFV_eff);
double IFN_av  = USE_PIFN * IFN_Emax * IFN_Cp / (IFN_Cp + IFN_EC50);
double IFN_im  = USE_PIFN * IFN_immE * IFN_Cp / (IFN_Cp + IFN_EC50);
double IFN_ccc_eff = USE_PIFN * IFN_cccE * IFN_Cp / (IFN_Cp + IFN_EC50);
double siRNA_eff = IAGM_dose * siRNA_Emax;
double INN_av  = IFN_antiv * IFN_inn / (IFN_inn + 0.5);
double CTL_eff = CTL_resp * (1.0 - Exhaust * 0.7);
if(NEWIND <= 1) {
  T_cell_0 = T0/1e6; I_cell_0 = I0/1e6; V_dna_0 = V0;
  ccc_DNA_0 = ccc0; HBsAg_C_0 = Ag0; CTL_resp_0 = CTL0;
  Exhaust_0 = ExhSS; IFN_inn_0 = 0.1; ALT_val_0 = ALT0;
  Fibrosis_0 = FIB0; HCC_risk_0 = 0.0; HSC_act_0 = 0.1;
}

$ODE
dxdt_ETV_gut = -ETV_ka * ETV_gut;
dxdt_ETV_C   =  ETV_F*ETV_ka*ETV_gut - (ETV_CL/ETV_Vc)*ETV_C - ETV_kTP*ETV_C*ETV_Kp/ETV_Vc;
dxdt_ETV_TP  =  ETV_kTP*ETV_C*ETV_Kp/ETV_Vc - ETV_kTP_d*ETV_TP;
dxdt_TDF_gut = -TDF_ka * TDF_gut;
dxdt_TDF_C   =  TDF_F*TDF_ka*TDF_gut - (TDF_CL/TDF_Vc)*TDF_C - TFV_kDP*TDF_C/TDF_Vc;
dxdt_TFV_DP  =  TFV_kDP*TDF_C/TDF_Vc - TFV_kDP_d*TFV_DP;
dxdt_IFN_SC  = -IFN_ka * IFN_SC;
dxdt_IFN_C   =  IFN_F*IFN_ka*IFN_SC - (IFN_CL/IFN_Vc)*IFN_C;
double T = T_cell; double I = I_cell; double V = V_dna; double ccc = ccc_DNA;
dxdt_T_cell  = s_T/1e6 - d_T*T - beta*T*V;
double CTL_kill = alpha_CTL * CTL_eff * I / (I + 0.1);
dxdt_I_cell  = beta*T*V - delta_I*I*(1.0+CTL_eff) - d_T*I;
dxdt_V_dna   = p_V*I*(1.0-NUC_eff)*(1.0-INN_av) - c_V*V;
double ccc_input = lambda_ccc*I*(1.0-NUC_eff)*(1.0-IFN_ccc_eff);
double ccc_loss  = delta_ccc*ccc + alpha_CTL*CTL_eff*0.2*ccc;
dxdt_ccc_DNA = ccc_input - ccc_loss;
dxdt_HBsAg_C = k_Ag*ccc - c_Ag*HBsAg_C*(1.0+siRNA_eff);
dxdt_CTL_resp = k_CTL_exp*V/(V+1e5)*(1.0+IFN_im)*(1.0-CTL_resp) - d_CTL*CTL_resp;
double Ag_load = HBsAg_C/(HBsAg_C+100.0);
dxdt_Exhaust  = k_exhaust*Ag_load*(1.0-Exhaust) - r_exhaust*Exhaust - IFN_im*0.1*Exhaust;
dxdt_IFN_inn  = k_IFNinn*0.05*V/(V+1e6) - d_IFNinn*IFN_inn + IFN_av*0.5;
dxdt_ALT_val  = delta_I*I*CTL_eff*k_ALT_I - k_ALT_ret*(ALT_val-ALT0);
double inflam_norm = (ALT_val-ALT0)/(ALT0+10.0);
dxdt_HSC_act  = k_HSC*inflam_norm*(1.0-HSC_act) - d_HSC*HSC_act*(V<100?1.5:1.0);
double prog_rate = k_fibrosis*HSC_act*(ALT_val/ALT0);
double reg_rate  = k_fib_reg*(V_dna<100?1.5:0.5);
dxdt_Fibrosis = prog_rate*(4.0-Fibrosis)/4.0 - reg_rate*Fibrosis/4.0;
dxdt_HCC_risk = k_HCC*Fibrosis*(1.0+0.5*(V_dna>1e4?1.0:0.0));

$TABLE
capture ETV_Cp  = ETV_C/ETV_Vc;
capture TFV_Cp  = TDF_C/TDF_Vc;
capture IFN_Cp  = IFN_C/IFN_Vc;
capture V_log10 = log10(V_dna+1.0);
capture Ag_log10= log10(HBsAg_C+0.01);
capture ccc_log = log10(ccc_DNA+0.001);
capture ALT_out = ALT_val;
capture Fib_out = Fibrosis;
capture HCC_out = HCC_risk;
capture CTL_out = CTL_resp;
capture Exh_out = Exhaust;
capture HSC_out = HSC_act;
capture NUC_eff = 1.0-(1.0-USE_ETV*ETV_Imax*ETV_TP/(ETV_TP+ETV_IC50))*(1.0-USE_TDF*TFV_Imax*TFV_DP/(TFV_DP+TFV_IC50));
capture IFN_AV  = USE_PIFN*IFN_Emax*IFN_Cp/(IFN_Cp+IFN_EC50);
'

chb_mod <- mcode("CHB_Shiny", chb_code, quiet = TRUE)

## ─────────────────────────────────────────────────────────────
## Simulation Helper
## ─────────────────────────────────────────────────────────────

run_sim <- function(V0, Ag0, ALT0, FIB0, CTL0,
                    use_etv, etv_dose,
                    use_tdf, tdf_dose,
                    use_pifn, pifn_dose, pifn_wks,
                    use_sirna, sim_yrs) {

  sim_end <- sim_yrs * 365
  params <- list(
    V0      = V0,    Ag0    = Ag0,   ALT0 = ALT0,
    FIB0    = FIB0,  CTL0   = CTL0,
    USE_ETV = as.numeric(use_etv),
    USE_TDF = as.numeric(use_tdf),
    USE_PIFN= as.numeric(use_pifn),
    IAGM_dose = as.numeric(use_sirna)
  )

  doses <- list()
  if(use_etv) doses[[length(doses)+1]] <-
    ev(cmt="ETV_gut", amt=etv_dose, ii=24, addl=sim_end-1)
  if(use_tdf) doses[[length(doses)+1]] <-
    ev(cmt="TDF_gut", amt=tdf_dose, ii=24, addl=sim_end-1)
  if(use_pifn) doses[[length(doses)+1]] <-
    ev(cmt="IFN_SC",  amt=pifn_dose, ii=168, addl=pifn_wks-1)
  if(length(doses)==0) doses[[1]] <- ev(amt=0, cmt=1, time=0)
  all_ev <- do.call(c, doses)

  chb_mod %>%
    param(params) %>%
    ev(all_ev) %>%
    mrgsim(end=sim_end, delta=1,
           recover=paste0("V_log10,Ag_log10,ccc_log,ALT_out,Fib_out,",
                          "HCC_out,CTL_out,Exh_out,HSC_out,",
                          "ETV_Cp,TFV_Cp,IFN_Cp,NUC_eff,IFN_AV")) %>%
    as_tibble() %>%
    mutate(time_yr = time/365)
}

## ─────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "CHB QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",   tabName="patient",   icon=icon("user")),
      menuItem("PK Profiles",       tabName="pk",        icon=icon("pills")),
      menuItem("Viral & HBsAg",     tabName="viral",     icon=icon("virus")),
      menuItem("Clinical Endpoints",tabName="clinical",  icon=icon("chart-line")),
      menuItem("Scenario Compare",  tabName="scenario",  icon=icon("exchange-alt")),
      menuItem("Biomarkers",        tabName="biomarker", icon=icon("flask"))
    )
  ),

  dashboardBody(
    tabItems(

      ## ── TAB 1: Patient Profile ──────────────────────────────
      tabItem(tabName="patient",
        fluidRow(
          box(title="Patient Demographics", width=4, status="danger",
            sliderInput("BWT",  "Body Weight (kg)",   40, 120, 65, 1),
            selectInput("phase","Disease Phase",
                        choices=c("Immune Tolerant (HBeAg+, low ALT)"="tolerant",
                                  "Immune Active (HBeAg+, high ALT)"="active",
                                  "Inactive Carrier (HBeAg−)"="inactive",
                                  "HBeAg-negative Hepatitis"="hbeag_neg")),
            numericInput("ALT0","Baseline ALT (IU/L)", 80, 10, 500, 10),
            numericInput("FIB0","Baseline Fibrosis (0–4)", 1.5, 0, 4, 0.5)
          ),
          box(title="Baseline Virological Parameters", width=4, status="warning",
            selectInput("V0_sel","Baseline HBV DNA",
                        choices=c("High: 10^8 IU/mL"  = "1e8",
                                  "Medium: 10^7 IU/mL" = "1e7",
                                  "Low: 10^5 IU/mL"   = "1e5",
                                  "Very Low: 10^3 IU/mL"= "1e3")),
            selectInput("Ag0_sel","Baseline HBsAg (IU/mL)",
                        choices=c(">10,000 (SVP excess)"    = "30000",
                                  "1,000–10,000 (moderate)" = "3000",
                                  "100–1,000 (lower)"       = "500",
                                  "<100 (functional cure)"  = "50")),
            sliderInput("CTL0","CTL Response Strength (0–1)", 0.1, 0.9, 0.5, 0.05)
          ),
          box(title="Treatment Selection", width=4, status="primary",
            checkboxInput("use_etv",   "Entecavir (ETV)", FALSE),
            conditionalPanel("input.use_etv",
              selectInput("etv_dose_sel","ETV Dose",
                          choices=c("0.5 mg QD (naïve)"=0.5, "1.0 mg QD (LamivR)"=1.0))),
            checkboxInput("use_tdf",   "TDF 300 mg QD", FALSE),
            checkboxInput("use_pifn",  "Peg-IFN-α2a 180 µg QW", FALSE),
            conditionalPanel("input.use_pifn",
              sliderInput("pifn_wks","Peg-IFN Duration (weeks)", 24, 72, 48, 4)),
            checkboxInput("use_sirna", "GalNAc-siRNA (HBsAg knockdown)", FALSE),
            sliderInput("sim_yrs","Simulation Duration (years)", 1, 10, 5, 1),
            actionButton("run_sim","▶  Run Simulation", class="btn-success btn-lg",
                         width="100%")
          )
        ),
        fluidRow(
          box(title="Baseline Parameter Summary", width=12, status="info",
            DTOutput("baseline_tbl")
          )
        )
      ),

      ## ── TAB 2: PK Profiles ─────────────────────────────────
      tabItem(tabName="pk",
        fluidRow(
          box(title="ETV Plasma Concentration", width=6, status="primary",
              plotlyOutput("pk_etv_plot", height="300px")),
          box(title="ETV-TP Intracellular (Active)", width=6, status="primary",
              plotlyOutput("pk_etvtp_plot", height="300px"))
        ),
        fluidRow(
          box(title="TFV Plasma Concentration", width=6, status="success",
              plotlyOutput("pk_tdf_plot", height="300px")),
          box(title="Peg-IFN-α2a Plasma", width=6, status="warning",
              plotlyOutput("pk_ifn_plot", height="300px"))
        ),
        fluidRow(
          box(title="Drug Effect (% RT Inhibition)", width=12, status="danger",
              plotlyOutput("pk_effect_plot", height="280px"))
        )
      ),

      ## ── TAB 3: Viral & HBsAg Dynamics ──────────────────────
      tabItem(tabName="viral",
        fluidRow(
          box(title="Serum HBV DNA Kinetics", width=6, status="danger",
              plotlyOutput("vir_vl_plot", height="340px")),
          box(title="HBsAg Kinetics — Path to Functional Cure", width=6, status="warning",
              plotlyOutput("vir_hbsag_plot", height="340px"))
        ),
        fluidRow(
          box(title="Intrahepatic cccDNA Kinetics", width=6, status="primary",
              plotlyOutput("vir_cccdna_plot", height="300px")),
          box(title="Virologic Response Summary", width=6, status="info",
              DTOutput("vir_response_tbl"))
        )
      ),

      ## ── TAB 4: Clinical Endpoints ───────────────────────────
      tabItem(tabName="clinical",
        fluidRow(
          box(title="ALT Kinetics — Biochemical Response", width=6, status="warning",
              plotlyOutput("clin_alt_plot", height="300px")),
          box(title="Fibrosis Score (Metavir)", width=6, status="danger",
              plotlyOutput("clin_fib_plot", height="300px"))
        ),
        fluidRow(
          box(title="Cumulative HCC Risk", width=6, status="danger",
              plotlyOutput("clin_hcc_plot", height="300px")),
          box(title="HSC Activation Index", width=6, status="primary",
              plotlyOutput("clin_hsc_plot", height="300px"))
        ),
        fluidRow(
          box(title="Clinical Milestones", width=12, status="success",
              DTOutput("clin_milestone_tbl"))
        )
      ),

      ## ── TAB 5: Scenario Comparison ─────────────────────────
      tabItem(tabName="scenario",
        fluidRow(
          box(title="Scenario Comparison Setup", width=4, status="primary",
            p("Predefined scenarios using current patient profile:"),
            checkboxGroupInput("scen_sel","Include Scenarios:",
              choices = c("Untreated"              = "untreated",
                          "Entecavir 0.5 mg QD"   = "etv",
                          "TDF 300 mg QD"          = "tdf",
                          "Peg-IFN × 48wks"        = "pifn",
                          "ETV + Peg-IFN combo"    = "etv_pifn",
                          "ETV + siRNA add-on"     = "etv_sirna"),
              selected = c("untreated","etv","tdf","pifn")),
            actionButton("run_compare","▶  Compare All", class="btn-primary",
                         width="100%")
          ),
          box(title="Year-1 Comparison: HBV DNA (log IU/mL)", width=8, status="danger",
              plotlyOutput("scen_vl_plot", height="320px"))
        ),
        fluidRow(
          box(title="Year-1 HBsAg by Scenario", width=6, status="warning",
              plotlyOutput("scen_ag_plot", height="300px")),
          box(title="Year-3 Fibrosis by Scenario", width=6, status="danger",
              plotlyOutput("scen_fib_plot", height="300px"))
        ),
        fluidRow(
          box(title="Scenario Summary Table", width=12, status="info",
              DTOutput("scen_summary_tbl"))
        )
      ),

      ## ── TAB 6: Biomarkers ──────────────────────────────────
      tabItem(tabName="biomarker",
        fluidRow(
          box(title="CTL Response & T-cell Exhaustion", width=6, status="primary",
              plotlyOutput("bio_ctl_plot", height="320px")),
          box(title="Innate IFN Signaling", width=6, status="success",
              plotlyOutput("bio_ifn_plot", height="320px"))
        ),
        fluidRow(
          box(title="Biomarker at Key Timepoints", width=12, status="info",
              DTOutput("bio_tbl"))
        ),
        fluidRow(
          box(title="Treatment Response Radar — Year 1", width=6, status="warning",
              plotlyOutput("bio_radar_plot", height="380px")),
          box(title="HBsAg Decline Rate (log IU/mL/yr)", width=6, status="danger",
              plotlyOutput("bio_hbsag_rate_plot", height="380px"))
        )
      )
    )
  )
)

## ─────────────────────────────────────────────────────────────
## SERVER
## ─────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # Reactive: run simulation when button pressed
  sim_data <- eventReactive(input$run_sim, {
    withProgress(message="Running CHB simulation...", value=0.5, {
      V0  <- as.numeric(input$V0_sel)
      Ag0 <- as.numeric(input$Ag0_sel)
      run_sim(
        V0      = V0,
        Ag0     = Ag0,
        ALT0    = input$ALT0,
        FIB0    = input$FIB0,
        CTL0    = input$CTL0,
        use_etv = input$use_etv,
        etv_dose= as.numeric(isolate(input$etv_dose_sel) %||% 0.5),
        use_tdf = input$use_tdf,
        tdf_dose= 300,
        use_pifn= input$use_pifn,
        pifn_dose=180,
        pifn_wks= isolate(input$pifn_wks) %||% 48,
        use_sirna=input$use_sirna,
        sim_yrs = input$sim_yrs
      )
    })
  })

  # Initialize with default simulation
  observe({
    if(is.null(isolate(sim_data()))) {
      isolate({ input$run_sim })
    }
  })

  default_data <- reactive({
    run_sim(V0=1e7, Ag0=1000, ALT0=80, FIB0=1.5, CTL0=0.5,
            use_etv=TRUE, etv_dose=0.5,
            use_tdf=FALSE, tdf_dose=300,
            use_pifn=FALSE, pifn_dose=180, pifn_wks=48,
            use_sirna=FALSE, sim_yrs=3)
  })

  get_data <- reactive({
    if(input$run_sim == 0) default_data() else sim_data()
  })

  ## ── TAB 1: Baseline table ─────────────────────────────────
  output$baseline_tbl <- renderDT({
    datatable(
      data.frame(
        Parameter = c("Body Weight","Baseline HBV DNA","Baseline HBsAg","Baseline ALT",
                      "Fibrosis Score","CTL Strength","Disease Phase",
                      "ETV","TDF","Peg-IFN","siRNA"),
        Value     = c(paste(input$BWT,"kg"),
                      input$V0_sel,
                      paste(input$Ag0_sel,"IU/mL"),
                      paste(input$ALT0,"IU/L"),
                      paste(input$FIB0,"(Metavir)"),
                      input$CTL0,
                      input$phase,
                      ifelse(input$use_etv, "Yes","No"),
                      ifelse(input$use_tdf, "Yes","No"),
                      ifelse(input$use_pifn,"Yes","No"),
                      ifelse(input$use_sirna,"Yes","No"))
      ), options=list(pageLength=15, dom='t'), rownames=FALSE
    )
  })

  ## ── TAB 2: PK Plots ───────────────────────────────────────
  output$pk_etv_plot <- renderPlotly({
    df <- get_data() %>% filter(time <= 30)
    plot_ly(df, x=~time, y=~ETV_Cp, type="scatter", mode="lines",
            line=list(color="#1E88E5",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="ETV (ng/mL)"),
             title="ETV Plasma (First 30 days)")
  })

  output$pk_etvtp_plot <- renderPlotly({
    df <- get_data() %>% filter(time <= 60)
    plot_ly(df, x=~time, y=~ETV_TP, type="scatter", mode="lines",
            line=list(color="#5C6BC0",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="ETV-TP (µM)"),
             title="ETV-TP Intracellular")
  })

  output$pk_tdf_plot <- renderPlotly({
    df <- get_data() %>% filter(time <= 30)
    plot_ly(df, x=~time, y=~TFV_Cp, type="scatter", mode="lines",
            line=list(color="#43A047",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="TFV (ng/mL)"),
             title="TFV Plasma (First 30 days)")
  })

  output$pk_ifn_plot <- renderPlotly({
    df <- get_data() %>% filter(time <= min(400, max(get_data()$time)))
    plot_ly(df, x=~time, y=~IFN_Cp, type="scatter", mode="lines",
            line=list(color="#FB8C00",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Peg-IFN (ng/mL equiv)"),
             title="Peg-IFN-α2a Plasma")
  })

  output$pk_effect_plot <- renderPlotly({
    df <- get_data()
    plot_ly(df, x=~time_yr) %>%
      add_lines(y=~NUC_eff*100, name="NUC (ETV/TDF) Effect", line=list(color="#1E88E5")) %>%
      add_lines(y=~IFN_AV*100,  name="Peg-IFN Antiviral",    line=list(color="#FB8C00")) %>%
      layout(xaxis=list(title="Time (years)"),
             yaxis=list(title="% RT Inhibition", range=c(0,105)),
             title="Antiviral Drug Effect on RT/Replication",
             legend=list(x=0.7, y=0.5))
  })

  ## ── TAB 3: Viral Dynamics ─────────────────────────────────
  output$vir_vl_plot <- renderPlotly({
    df <- get_data()
    plot_ly(df, x=~time_yr, y=~V_log10, type="scatter", mode="lines",
            line=list(color="#E53935",width=2)) %>%
      add_segments(x=0, xend=max(df$time_yr), y=log10(20), yend=log10(20),
                   line=list(dash="dash", color="gray")) %>%
      layout(xaxis=list(title="Time (years)"),
             yaxis=list(title="HBV DNA (log₁₀ IU/mL)"),
             title="Serum HBV DNA")
  })

  output$vir_hbsag_plot <- renderPlotly({
    df <- get_data()
    plot_ly(df, x=~time_yr, y=~Ag_log10, type="scatter", mode="lines",
            line=list(color="#FF8F00",width=2.5)) %>%
      add_segments(x=0, xend=max(df$time_yr), y=log10(0.05), yend=log10(0.05),
                   line=list(dash="dash", color="red")) %>%
      layout(xaxis=list(title="Time (years)"),
             yaxis=list(title="HBsAg (log₁₀ IU/mL)"),
             title="HBsAg — Functional Cure Threshold")
  })

  output$vir_cccdna_plot <- renderPlotly({
    df <- get_data()
    plot_ly(df, x=~time_yr, y=~ccc_log, type="scatter", mode="lines",
            line=list(color="#7B1FA2",width=2)) %>%
      layout(xaxis=list(title="Time (years)"),
             yaxis=list(title="cccDNA (log₁₀ copies/cell equiv)"),
             title="Intrahepatic cccDNA Pool")
  })

  output$vir_response_tbl <- renderDT({
    df <- get_data()
    tps <- c(0, 52, 104, 156)  # weeks
    tps_d <- tps * 7
    tab <- df %>%
      filter(time %in% tps_d | time == max(time)) %>%
      mutate(
        Week     = round(time/7),
        `HBV DNA (log₁₀)` = round(V_log10, 2),
        `HBsAg (log₁₀)`   = round(Ag_log10, 2),
        `ViroResp (<20 IU/mL)` = V_log10 < log10(21),
        `Func Cure (<0.05 IU/mL)` = Ag_log10 < log10(0.06)
      ) %>%
      select(Week, `HBV DNA (log₁₀)`, `HBsAg (log₁₀)`,
             `ViroResp (<20 IU/mL)`, `Func Cure (<0.05 IU/mL)`)
    datatable(tab, options=list(dom='t',pageLength=10), rownames=FALSE)
  })

  ## ── TAB 4: Clinical Endpoints ─────────────────────────────
  output$clin_alt_plot <- renderPlotly({
    df <- get_data()
    plot_ly(df, x=~time_yr, y=~ALT_out, type="scatter", mode="lines",
            line=list(color="#EF6C00",width=2)) %>%
      add_segments(x=0, xend=max(df$time_yr), y=40, yend=40,
                   line=list(dash="dash",color="green")) %>%
      layout(xaxis=list(title="Time (years)"),
             yaxis=list(title="ALT (IU/L)"),
             title="ALT Dynamics")
  })

  output$clin_fib_plot <- renderPlotly({
    df <- get_data()
    plot_ly(df, x=~time_yr, y=~Fib_out, type="scatter", mode="lines",
            line=list(color="#AD1457",width=2)) %>%
      layout(xaxis=list(title="Time (years)"),
             yaxis=list(title="Fibrosis (Metavir)", range=c(0,4)),
             title="Hepatic Fibrosis Progression")
  })

  output$clin_hcc_plot <- renderPlotly({
    df <- get_data()
    plot_ly(df, x=~time_yr, y=~HCC_out*100, type="scatter", mode="lines",
            line=list(color="#B71C1C",width=2)) %>%
      layout(xaxis=list(title="Time (years)"),
             yaxis=list(title="Cumulative HCC Risk (%)"),
             title="HCC Risk Accumulation")
  })

  output$clin_hsc_plot <- renderPlotly({
    df <- get_data()
    plot_ly(df, x=~time_yr, y=~HSC_out, type="scatter", mode="lines",
            line=list(color="#880E4F",width=2)) %>%
      layout(xaxis=list(title="Time (years)"),
             yaxis=list(title="HSC Activation (0–1)"),
             title="Hepatic Stellate Cell Activation")
  })

  output$clin_milestone_tbl <- renderDT({
    df <- get_data()
    last <- df %>% slice_tail(n=1)
    yr1  <- df %>% filter(abs(time_yr - 1) == min(abs(time_yr - 1))) %>% slice(1)
    milestones <- data.frame(
      Milestone = c("HBV DNA <2,000 IU/mL","HBV DNA <20 IU/mL (undetect.)",
                    "ALT Normalization (<40)","HBsAg <100 IU/mL",
                    "Functional Cure (<0.05 IU/mL)","Fibrosis Regression"),
      `Year 1`  = c(yr1$V_log10 < log10(2001),
                    yr1$V_log10 < log10(21),
                    yr1$ALT_out < 40,
                    yr1$Ag_log10 < 2,
                    yr1$Ag_log10 < log10(0.06),
                    yr1$Fib_out < last$FIB0),
      `Year End`= c(last$V_log10 < log10(2001),
                    last$V_log10 < log10(21),
                    last$ALT_out < 40,
                    last$Ag_log10 < 2,
                    last$Ag_log10 < log10(0.06),
                    last$Fib_out < df$Fib_out[1])
    )
    datatable(milestones, options=list(dom='t',pageLength=10), rownames=FALSE)
  })

  ## ── TAB 5: Scenario Comparison ────────────────────────────
  scen_data <- eventReactive(input$run_compare, {
    withProgress(message="Running scenario comparisons...", value=0.3, {
      V0  <- as.numeric(input$V0_sel)
      Ag0 <- as.numeric(input$Ag0_sel)
      sel <- input$scen_sel

      scen_list <- list(
        "untreated" = list(etv=FALSE, tdf=FALSE, pifn=FALSE, sirna=FALSE, label="Untreated"),
        "etv"       = list(etv=TRUE,  tdf=FALSE, pifn=FALSE, sirna=FALSE, label="Entecavir 0.5 mg QD"),
        "tdf"       = list(etv=FALSE, tdf=TRUE,  pifn=FALSE, sirna=FALSE, label="TDF 300 mg QD"),
        "pifn"      = list(etv=FALSE, tdf=FALSE, pifn=TRUE,  sirna=FALSE, label="Peg-IFN × 48wks"),
        "etv_pifn"  = list(etv=TRUE,  tdf=FALSE, pifn=TRUE,  sirna=FALSE, label="ETV + Peg-IFN"),
        "etv_sirna" = list(etv=TRUE,  tdf=FALSE, pifn=FALSE, sirna=TRUE,  label="ETV + siRNA")
      )

      bind_rows(lapply(sel, function(s) {
        sc <- scen_list[[s]]
        run_sim(V0, Ag0, input$ALT0, input$FIB0, input$CTL0,
                sc$etv, 0.5, sc$tdf, 300, sc$pifn, 180, 48, sc$sirna, 3) %>%
          mutate(Scenario = sc$label)
      }))
    })
  })

  output$scen_vl_plot <- renderPlotly({
    req(scen_data())
    df <- scen_data()
    p <- ggplot(df, aes(x=time_yr, y=V_log10, color=Scenario)) +
      geom_line(size=1) +
      geom_hline(yintercept=log10(20), linetype="dashed", color="gray40") +
      labs(x="Time (years)", y="HBV DNA (log₁₀ IU/mL)") +
      theme_bw()
    ggplotly(p)
  })

  output$scen_ag_plot <- renderPlotly({
    req(scen_data())
    df <- scen_data()
    p <- ggplot(df, aes(x=time_yr, y=Ag_log10, color=Scenario)) +
      geom_line(size=1) +
      geom_hline(yintercept=log10(0.05), linetype="dashed", color="red") +
      labs(x="Time (years)", y="HBsAg (log₁₀ IU/mL)") +
      theme_bw()
    ggplotly(p)
  })

  output$scen_fib_plot <- renderPlotly({
    req(scen_data())
    df <- scen_data()
    p <- ggplot(df, aes(x=time_yr, y=Fib_out, color=Scenario)) +
      geom_line(size=1) +
      scale_y_continuous(limits=c(0,4)) +
      labs(x="Time (years)", y="Fibrosis (Metavir)") +
      theme_bw()
    ggplotly(p)
  })

  output$scen_summary_tbl <- renderDT({
    req(scen_data())
    df <- scen_data()
    tab <- df %>%
      filter(abs(time_yr - 1) == min(abs(time_yr - 1))) %>%
      group_by(Scenario) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(
        `HBV DNA Y1 (log)` = round(V_log10, 2),
        `HBsAg Y1 (log)`   = round(Ag_log10, 2),
        `ALT Y1`           = round(ALT_out, 1),
        `DNA Undetect.`    = V_log10 < log10(21),
        `Func Cure`        = Ag_log10 < log10(0.06)
      ) %>%
      select(Scenario, `HBV DNA Y1 (log)`, `HBsAg Y1 (log)`, `ALT Y1`,
             `DNA Undetect.`, `Func Cure`)
    datatable(tab, options=list(dom='t', pageLength=8), rownames=FALSE)
  })

  ## ── TAB 6: Biomarkers ─────────────────────────────────────
  output$bio_ctl_plot <- renderPlotly({
    df <- get_data()
    plot_ly(df, x=~time_yr) %>%
      add_lines(y=~CTL_out, name="CTL Response", line=list(color="#1565C0",width=2)) %>%
      add_lines(y=~Exh_out, name="T cell Exhaustion", line=list(color="#78909C",width=2,dash="dash")) %>%
      layout(xaxis=list(title="Time (years)"),
             yaxis=list(title="Normalized Index (0–1)", range=c(0,1)),
             title="CD8+ CTL & T-cell Exhaustion",
             legend=list(x=0.65, y=0.9))
  })

  output$bio_ifn_plot <- renderPlotly({
    df <- get_data()
    plot_ly(df, x=~time_yr, y=~IFN_AV, type="scatter", mode="lines",
            line=list(color="#2E7D32",width=2)) %>%
      layout(xaxis=list(title="Time (years)"),
             yaxis=list(title="IFN Antiviral Effect (0–1)"),
             title="Innate IFN Signaling & Antiviral Effect")
  })

  output$bio_tbl <- renderDT({
    df <- get_data()
    tps_d <- c(0, 90, 180, 365, 730)
    tab <- df %>%
      filter(sapply(time, function(t) min(abs(t - tps_d)) < 2)) %>%
      group_by(time) %>% slice(1) %>% ungroup() %>%
      mutate(
        Day    = time,
        `CTL`  = round(CTL_out, 3),
        `Exhaustion` = round(Exh_out, 3),
        `IFN innate` = round(IFN_AV, 3),
        `HSC act.`   = round(HSC_out, 3),
        `cccDNA(log)`= round(ccc_log, 2),
        `HBsAg(log)` = round(Ag_log10, 2),
        `HBV DNA(log)`= round(V_log10, 2)
      ) %>%
      select(Day, `CTL`, `Exhaustion`, `IFN innate`, `HSC act.`,
             `cccDNA(log)`, `HBsAg(log)`, `HBV DNA(log)`)
    datatable(tab, options=list(dom='t',pageLength=8), rownames=FALSE)
  })

  output$bio_radar_plot <- renderPlotly({
    df <- get_data()
    yr1 <- df %>% filter(abs(time_yr-1) == min(abs(time_yr-1))) %>% slice(1)
    # Normalize metrics for radar (0=worst, 1=best)
    categories <- c("DNA Suppression","HBsAg Reduction",
                     "ALT Normalization","Fibrosis Control",
                     "CTL Recovery","Exh. Reduction")
    values <- c(
      max(0, min(1, 1 - yr1$V_log10/8)),
      max(0, min(1, 1 - yr1$Ag_log10/5)),
      max(0, min(1, 1 - (yr1$ALT_out-40)/100)),
      max(0, min(1, 1 - yr1$Fib_out/4)),
      yr1$CTL_out,
      1 - yr1$Exh_out
    )
    plot_ly(type="scatterpolar", mode="lines+markers",
            r=c(values, values[1]),
            theta=c(categories, categories[1]),
            fill="toself",
            line=list(color="#1565C0"),
            fillcolor="rgba(30,136,229,0.3)") %>%
      layout(polar=list(radialaxis=list(range=c(0,1))),
             title="Treatment Response Radar (Year 1)")
  })

  output$bio_hbsag_rate_plot <- renderPlotly({
    df <- get_data()
    df_rate <- df %>%
      mutate(dAg = c(NA, diff(Ag_log10)),
             dt  = c(NA, diff(time_yr))) %>%
      filter(!is.na(dAg), time_yr > 0) %>%
      mutate(rate = dAg / dt)
    plot_ly(df_rate, x=~time_yr, y=~rate, type="scatter", mode="lines",
            line=list(color="#FF8F00",width=2)) %>%
      add_segments(x=0, xend=max(df_rate$time_yr, na.rm=TRUE), y=0, yend=0,
                   line=list(dash="dash", color="gray")) %>%
      layout(xaxis=list(title="Time (years)"),
             yaxis=list(title="ΔHBsAg (log IU/mL/year)"),
             title="HBsAg Decline Rate Over Time")
  })
}

## ─────────────────────────────────────────────────────────────
## LAUNCH APP
## ─────────────────────────────────────────────────────────────

shinyApp(ui = ui, server = server)
