# ============================================================================
# CLL QSP Shiny App — Interactive Dashboard
# Tabs: Patient Profile | PK | PD Biomarkers | Clinical Endpoints |
#       Scenario Comparison | Genetic Risk
# ============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(tidyr)
library(DT)

# ── Embed mrgsolve model ──────────────────────────────────────────────────────
cll_code <- '
$PARAM
Ka_IB=0.50,Vd_IB=10000,CL_IB=980,F_IB=0.25,
kinact_BTK=0.10,Ki_BTK=1.5,kdeg_BTK=0.010,
Ka_VEN=0.30,V1_VEN=250,V2_VEN=500,CL_VEN=65,Q_VEN=10,F_VEN=0.50,
Ki_BCL2=0.01,BCL2_tot=100,kout_BCL2=0.10,
V1_OBI=3.4,V2_OBI=3.0,CL_OBI=0.020,Q_OBI=0.150,
Kd_CD20=0.001,kint_CD20=0.0003,ksyn_CD20=0.050,kdeg_CD20=0.005,CD20_0=100,
kprol_CLL=0.0030,Kmax_ALC=300,ALC_0=50,BM_0=70,LN_0=60,
Emax_BTK=0.70,EC50_BTK=50,Emax_BCL2=0.90,EC50_BCL2=35,Emax_CD20=0.75,EC50_CD20=50,
kegress=0.012,egress_thr=20,
kin_MCL1=0.008,kout_MCL1=0.050,MCL1_max=4.0,
kin_NK=0.005,kout_NK=0.020,NK_max=3.0

$CMT DEPOT_IB CENT_IB BTK_FREE BTK_OCC DEPOT_VEN CENT_VEN PERI_VEN
     BCL2_FREE BCL2_OCC CENT_OBI PERI_OBI CD20_FREE CD20_OCC
     ALC BM_CLL LN_CLL MCL1_ADAPT NK_ACT

$INIT DEPOT_IB=0,CENT_IB=0,BTK_FREE=100,BTK_OCC=0,
      DEPOT_VEN=0,CENT_VEN=0,PERI_VEN=0,BCL2_FREE=100,BCL2_OCC=0,
      CENT_OBI=0,PERI_OBI=0,CD20_FREE=100,CD20_OCC=0,
      ALC=50,BM_CLL=70,LN_CLL=60,MCL1_ADAPT=1.0,NK_ACT=1.0

$ODE
dxdt_DEPOT_IB=-Ka_IB*DEPOT_IB;
dxdt_CENT_IB=Ka_IB*F_IB*DEPOT_IB-(CL_IB/Vd_IB)*CENT_IB;
double C_IB_nM=(CENT_IB/Vd_IB)*(1000.0/440.5);
double k_inact=kinact_BTK*C_IB_nM/(Ki_BTK+C_IB_nM);
dxdt_BTK_FREE=kdeg_BTK*100.0-kdeg_BTK*BTK_FREE-k_inact*BTK_FREE;
dxdt_BTK_OCC=k_inact*BTK_FREE-kdeg_BTK*BTK_OCC;
double BTK_OCC_pct=BTK_OCC;
dxdt_DEPOT_VEN=-Ka_VEN*DEPOT_VEN;
dxdt_CENT_VEN=Ka_VEN*F_VEN*DEPOT_VEN-(CL_VEN+Q_VEN)/V1_VEN*CENT_VEN+Q_VEN/V2_VEN*PERI_VEN;
dxdt_PERI_VEN=Q_VEN/V1_VEN*CENT_VEN-Q_VEN/V2_VEN*PERI_VEN;
double C_VEN_nM=(CENT_VEN/V1_VEN)*(1000.0/868.4);
double BCL2_OCC_ss=BCL2_tot*C_VEN_nM/(Ki_BCL2+C_VEN_nM);
dxdt_BCL2_FREE=kout_BCL2*(BCL2_tot-BCL2_OCC_ss-BCL2_FREE);
dxdt_BCL2_OCC=kout_BCL2*(BCL2_OCC_ss-BCL2_OCC);
double BCL2_OCC_pct=(BCL2_FREE+BCL2_OCC>0.001)?BCL2_OCC/(BCL2_FREE+BCL2_OCC)*100.0:0;
double C_OBI_mgL=CENT_OBI/V1_OBI;
double CD20_OCC_ss=CD20_0*C_OBI_mgL/(Kd_CD20+C_OBI_mgL);
dxdt_CENT_OBI=-(CL_OBI+Q_OBI)/V1_OBI*CENT_OBI+Q_OBI/V2_OBI*PERI_OBI-kint_CD20*(CD20_OCC_ss-CD20_OCC)*V1_OBI*0.01;
dxdt_PERI_OBI=Q_OBI/V1_OBI*CENT_OBI-Q_OBI/V2_OBI*PERI_OBI;
dxdt_CD20_FREE=ksyn_CD20*CD20_0-kdeg_CD20*CD20_FREE-kint_CD20*(CD20_OCC_ss-CD20_OCC)*0.5;
dxdt_CD20_OCC=kint_CD20*(CD20_OCC_ss-CD20_OCC)*0.5-kdeg_CD20*CD20_OCC;
double CD20_OCC_pct=(CD20_FREE+CD20_OCC>0.001)?CD20_OCC/(CD20_FREE+CD20_OCC)*100.0:0;
double E_BTKi=Emax_BTK*BTK_OCC_pct/(EC50_BTK+BTK_OCC_pct);
double E_BCL2i=Emax_BCL2*BCL2_OCC_pct/(EC50_BCL2+BCL2_OCC_pct)/MCL1_ADAPT;
double E_CD20=Emax_CD20*CD20_OCC_pct*NK_ACT/(EC50_CD20+CD20_OCC_pct);
double kill_ALC=(E_BTKi*0.50+E_BCL2i*0.80+E_CD20*0.50)*ALC;
double kill_BM=(E_BTKi*0.30+E_BCL2i*0.90+E_CD20*0.40)*BM_CLL;
double kill_LN=(E_BTKi*0.60+E_BCL2i*0.70+E_CD20*0.70)*LN_CLL;
double do_eg=(BTK_OCC_pct>egress_thr)?1.0:0.0;
double egr_BM=do_eg*kegress*BM_CLL;
double egr_LN=do_eg*kegress*LN_CLL*0.6;
double ALC_pos=(ALC>0)?ALC:0;
dxdt_ALC=kprol_CLL*ALC_pos*(1.0-ALC_pos/Kmax_ALC)-kill_ALC+egr_BM+egr_LN;
if(ALC<0.001)dxdt_ALC=0;
double BM_pos=(BM_CLL>0)?BM_CLL:0;
dxdt_BM_CLL=kprol_CLL*0.8*BM_pos*(1.0-BM_pos/100.0)-kill_BM-egr_BM;
if(BM_CLL<0.001)dxdt_BM_CLL=0;
double LN_pos=(LN_CLL>0)?LN_CLL:0;
dxdt_LN_CLL=kprol_CLL*1.2*LN_pos*(1.0-LN_pos/100.0)-kill_LN-egr_LN;
if(LN_CLL<0.001)dxdt_LN_CLL=0;
double stim_MCL1=BCL2_OCC_pct/100.0;
dxdt_MCL1_ADAPT=kin_MCL1*stim_MCL1*(MCL1_max-MCL1_ADAPT)-kout_MCL1*(MCL1_ADAPT-1.0);
if(MCL1_ADAPT<1.0)dxdt_MCL1_ADAPT=0;
dxdt_NK_ACT=kin_NK*(CD20_OCC_pct/100.0)*(NK_max-NK_ACT)-kout_NK*(NK_ACT-1.0);
if(NK_ACT<1.0)dxdt_NK_ACT=0;

$TABLE
double C_IB_ngmL=(CENT_IB/Vd_IB)*1000.0;
double C_VEN_ngmL=(CENT_VEN/V1_VEN)*1000.0;
double C_OBI_ugmL=CENT_OBI/V1_OBI*1000.0;
double BTK_OCC_out=BTK_OCC;
double BCL2_OCC_out=(BCL2_FREE+BCL2_OCC>0.001)?BCL2_OCC/(BCL2_FREE+BCL2_OCC)*100.0:0;
double CD20_OCC_out=(CD20_FREE+CD20_OCC>0.001)?CD20_OCC/(CD20_FREE+CD20_OCC)*100.0:0;
double ALC_pch=(ALC_0>0)?(ALC-ALC_0)/ALC_0*100.0:0;
int CR_flag=(ALC<4.0&&BM_CLL<30.0&&LN_CLL<20.0)?1:0;
int PR_flag=(ALC_pch<-50.0&&!CR_flag)?1:0;
int MRD_neg=(ALC<0.1&&BM_CLL<5.0)?1:0;
double BURDEN=(ALC/Kmax_ALC*100.0+BM_CLL+LN_CLL)/3.0;

$CAPTURE ALC,BM_CLL,LN_CLL,BTK_OCC_out,BCL2_OCC_out,CD20_OCC_out,
         MCL1_ADAPT,NK_ACT,C_IB_ngmL,C_VEN_ngmL,C_OBI_ugmL,
         BURDEN,CR_flag,PR_flag,MRD_neg,ALC_pch
'

mod <- mcode("CLL_Shiny", cll_code, quiet = TRUE)

# ── Simulation helper ──────────────────────────────────────────────────────────
simulate_cll <- function(alc0, bm0, ln0, use_ib, use_ven, use_obi,
                          ib_dose, ven_peak_dose, end_days) {
  ev_list <- list()

  if (use_ib) {
    ev_list$IB <- ev(amt = ib_dose, cmt = "DEPOT_IB",
                     time = 0, ii = 24, addl = end_days - 1)
  }

  if (use_ven) {
    ramp_doses <- c(20, 50, 100, 200, min(ven_peak_dose, 400))
    ramp_times <- c(0, 168, 336, 504, 672)
    for (i in 1:5) {
      next_t <- if (i < 5) ramp_times[i + 1] else end_days * 24
      add_doses <- max(0, floor((next_t - ramp_times[i]) / 24) - 1)
      ev_list[[paste0("VEN", i)]] <- ev(
        amt = ramp_doses[i], cmt = "DEPOT_VEN",
        time = ramp_times[i], ii = 24, addl = add_doses
      )
    }
  }

  if (use_obi) {
    obi_t   <- c(0, 24, 336, 672, 1344, 2016, 2688, 3360)
    obi_amt <- c(100, 900, 1000, 1000, 1000, 1000, 1000, 1000)
    valid   <- obi_t <= end_days * 24
    for (i in which(valid)) {
      ev_list[[paste0("OBI", i)]] <- ev(amt = obi_amt[i], cmt = "CENT_OBI",
                                        time = obi_t[i])
    }
  }

  if (length(ev_list) == 0) {
    e <- ev(amt = 0, cmt = "DEPOT_IB", time = 1e9)
  } else {
    e <- do.call(c, ev_list)
  }

  p <- list(ALC_0 = alc0, BM_0 = bm0, LN_0 = ln0)
  mod2 <- param(mod, p)
  init_list <- list(ALC = alc0, BM_CLL = bm0, LN_CLL = ln0)
  mod2 <- init(mod2, init_list)

  mrgsim(mod2, ev = e, start = 0, end = end_days * 24, delta = 12) %>%
    as_tibble() %>%
    mutate(time_days = time / 24)
}

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "CLL QSP Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "profile",   icon = icon("user")),
      menuItem("PK Profiles",        tabName = "pk",        icon = icon("pills")),
      menuItem("PD Biomarkers",      tabName = "pd",        icon = icon("flask")),
      menuItem("Clinical Endpoints", tabName = "endpoints", icon = icon("chart-line")),
      menuItem("Scenario Comparison",tabName = "scenarios", icon = icon("balance-scale")),
      menuItem("Genetic Risk",       tabName = "genetics",  icon = icon("dna"))
    )
  ),
  dashboardBody(
    tabItems(

      # ── Tab 1: Patient Profile ──────────────────────────────────────────────
      tabItem("profile",
        fluidRow(
          box(title = "Patient Characteristics", status = "primary", solidHeader = TRUE,
              width = 4,
              numericInput("alc0",     "Baseline ALC (×10⁹/L)",    value = 50, min = 5,  max = 500),
              numericInput("bm0",      "BM Infiltration (%)",        value = 70, min = 10, max = 100),
              numericInput("ln0",      "LN Burden (relative %)",     value = 60, min = 10, max = 100),
              numericInput("end_days", "Simulation Duration (days)", value = 730, min = 30, max = 1095),
              hr(),
              selectInput("rai_stage", "Rai Stage",
                          choices = c("0 – Lymphocytosis only",
                                      "I – + Lymphadenopathy",
                                      "II – + Splenomegaly",
                                      "III – + Anemia",
                                      "IV – + Thrombocytopenia")),
              selectInput("ighv_status", "IGHV Status",
                          choices = c("Mutated (favourable)",
                                      "Unmutated (unfavourable)")),
              selectInput("cytogenetics", "Cytogenetics (FISH)",
                          choices = c("del(13q) only – low risk",
                                      "Normal – intermediate",
                                      "Trisomy 12 – intermediate",
                                      "del(11q) – high risk",
                                      "del(17p) – very high risk"))
          ),
          box(title = "CLL-IPI Risk Calculator", status = "warning", solidHeader = TRUE,
              width = 4,
              numericInput("age_ipi",   "Age (years)",       value = 68, min = 18, max = 100),
              numericInput("b2mg_ipi",  "β2-Microglobulin (mg/L)", value = 4.5, min = 0, max = 30),
              hr(),
              p(strong("CLL-IPI Score Components:")),
              p("• Age >65: +1"),
              p("• IGHV unmutated: +2"),
              p("• β2M >3.5 mg/L: +2"),
              p("• del(17p) or TP53 mut: +4"),
              p("• Clinical stage Binet B/C or Rai I-IV: +1"),
              hr(),
              verbatimTextOutput("ipi_score")
          ),
          box(title = "Disease Summary", status = "info", solidHeader = TRUE,
              width = 4,
              tableOutput("patient_summary")
          )
        )
      ),

      # ── Tab 2: PK Profiles ─────────────────────────────────────────────────
      tabItem("pk",
        fluidRow(
          box(title = "Treatment Selection", status = "primary", solidHeader = TRUE,
              width = 3,
              checkboxInput("use_ib",  "Ibrutinib",      value = TRUE),
              numericInput("ib_dose",  "Ibrutinib dose (mg QD)", value = 420, min = 140, max = 560),
              hr(),
              checkboxInput("use_ven", "Venetoclax",     value = FALSE),
              numericInput("ven_dose", "Peak VEN dose (mg QD)", value = 400, min = 20, max = 400),
              hr(),
              checkboxInput("use_obi", "Obinutuzumab",   value = FALSE),
              actionButton("run_sim", "Run Simulation", class = "btn-success")
          ),
          box(title = "Ibrutinib PK", status = "info", solidHeader = TRUE,
              width = 9, plotlyOutput("pk_ib_plot", height = "280px"))
        ),
        fluidRow(
          box(title = "Venetoclax PK", status = "info", solidHeader = TRUE,
              width = 6, plotlyOutput("pk_ven_plot", height = "280px")),
          box(title = "Obinutuzumab PK", status = "info", solidHeader = TRUE,
              width = 6, plotlyOutput("pk_obi_plot", height = "280px"))
        )
      ),

      # ── Tab 3: PD Biomarkers ───────────────────────────────────────────────
      tabItem("pd",
        fluidRow(
          box(title = "BTK Occupancy", status = "warning", solidHeader = TRUE,
              width = 6, plotlyOutput("btk_occ_plot", height = "300px")),
          box(title = "BCL-2 Occupancy", status = "warning", solidHeader = TRUE,
              width = 6, plotlyOutput("bcl2_occ_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "CD20 Receptor Occupancy", status = "warning", solidHeader = TRUE,
              width = 6, plotlyOutput("cd20_occ_plot", height = "300px")),
          box(title = "MCL-1 Resistance & NK Activation", status = "danger", solidHeader = TRUE,
              width = 6, plotlyOutput("adapt_plot", height = "300px"))
        )
      ),

      # ── Tab 4: Clinical Endpoints ──────────────────────────────────────────
      tabItem("endpoints",
        fluidRow(
          box(title = "Absolute Lymphocyte Count (ALC)", status = "success", solidHeader = TRUE,
              width = 8, plotlyOutput("alc_plot", height = "320px")),
          box(title = "IWCLL Response Summary", status = "success", solidHeader = TRUE,
              width = 4,
              tableOutput("response_table"),
              hr(),
              verbatimTextOutput("cr_time")
          )
        ),
        fluidRow(
          box(title = "Tumor Burden by Compartment", status = "primary", solidHeader = TRUE,
              width = 6, plotlyOutput("burden_plot", height = "300px")),
          box(title = "Composite Tumor Burden", status = "primary", solidHeader = TRUE,
              width = 6, plotlyOutput("composite_plot", height = "300px"))
        )
      ),

      # ── Tab 5: Scenario Comparison ─────────────────────────────────────────
      tabItem("scenarios",
        fluidRow(
          box(title = "Select Scenarios to Compare", status = "primary",
              solidHeader = TRUE, width = 12,
              checkboxGroupInput("compare_scen",
                label = NULL,
                choices = c("Ibrutinib monotherapy"      = "ib",
                            "Venetoclax monotherapy"     = "ven",
                            "Obinutuzumab monotherapy"   = "obi",
                            "Venetoclax + Obinutuzumab (CLL14)" = "ven_obi",
                            "Ibrutinib + Venetoclax"     = "ib_ven",
                            "Triple: IB + VEN + OBI"     = "triplet"),
                selected = c("ib", "ven_obi"),
                inline = TRUE
              ),
              actionButton("run_compare", "Compare Scenarios", class = "btn-primary")
          )
        ),
        fluidRow(
          box(title = "ALC Comparison", status = "info", solidHeader = TRUE,
              width = 6, plotlyOutput("cmp_alc", height = "320px")),
          box(title = "Composite Burden Comparison", status = "info", solidHeader = TRUE,
              width = 6, plotlyOutput("cmp_burden", height = "320px"))
        ),
        fluidRow(
          box(title = "BTK Occupancy Comparison", status = "warning", solidHeader = TRUE,
              width = 6, plotlyOutput("cmp_btk", height = "300px")),
          box(title = "Response Summary Table", status = "success", solidHeader = TRUE,
              width = 6, DTOutput("cmp_table"))
        )
      ),

      # ── Tab 6: Genetic Risk ────────────────────────────────────────────────
      tabItem("genetics",
        fluidRow(
          box(title = "Prognostic Marker Impact on Treatment Selection",
              status = "danger", solidHeader = TRUE, width = 12,
              p("Key genomic risk factors in CLL and their implications:"),
              DTOutput("genetics_table")
          )
        ),
        fluidRow(
          box(title = "del(17p)/TP53 — Ibrutinib vs Venetoclax Response",
              status = "warning", solidHeader = TRUE, width = 6,
              selectInput("del17p_flag", "Patient has del(17p)/TP53 mutation?",
                          choices = c("No", "Yes")),
              actionButton("run_del17p", "Simulate"),
              plotlyOutput("del17p_plot", height = "300px")
          ),
          box(title = "IGHV Mutation Status — Impact on BTKi Response",
              status = "warning", solidHeader = TRUE, width = 6,
              selectInput("ighv_sim", "IGHV Status for Simulation",
                          choices = c("Mutated (slower BCR)" = "mut",
                                      "Unmutated (faster BCR)" = "unmut")),
              actionButton("run_ighv", "Simulate"),
              plotlyOutput("ighv_plot", height = "300px")
          )
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Reactive: main simulation ───────────────────────────────────────────────
  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Running CLL simulation...", {
      simulate_cll(
        alc0 = input$alc0, bm0 = input$bm0, ln0 = input$ln0,
        use_ib = input$use_ib, use_ven = input$use_ven, use_obi = input$use_obi,
        ib_dose = input$ib_dose, ven_peak_dose = input$ven_dose,
        end_days = input$end_days
      )
    })
  }, ignoreNULL = FALSE)

  # ── CLL-IPI Score ───────────────────────────────────────────────────────────
  output$ipi_score <- renderText({
    score <- 0
    if (input$age_ipi > 65)             score <- score + 1
    if (grepl("Unmutated", input$ighv_status)) score <- score + 2
    if (input$b2mg_ipi > 3.5)           score <- score + 2
    if (grepl("del\\(17p\\)", input$cytogenetics)) score <- score + 4
    if (!grepl("0", input$rai_stage))   score <- score + 1
    risk <- dplyr::case_when(
      score <= 1 ~ "Low (10-yr OS ~92%)",
      score <= 3 ~ "Intermediate (10-yr OS ~79%)",
      score <= 6 ~ "High (10-yr OS ~61%)",
      TRUE       ~ "Very High (10-yr OS ~23%)"
    )
    paste0("CLL-IPI Score: ", score, "\nRisk group: ", risk)
  })

  # ── Patient summary ─────────────────────────────────────────────────────────
  output$patient_summary <- renderTable({
    data.frame(
      Parameter = c("ALC", "BM Infiltration", "LN Burden", "Rai Stage",
                    "IGHV", "Cytogenetics"),
      Value     = c(paste0(input$alc0, " ×10⁹/L"),
                    paste0(input$bm0, "%"),
                    paste0(input$ln0, "%"),
                    input$rai_stage,
                    input$ighv_status,
                    input$cytogenetics),
      stringsAsFactors = FALSE
    )
  })

  # ── PK plots ────────────────────────────────────────────────────────────────
  output$pk_ib_plot <- renderPlotly({
    d <- sim_data()
    if (max(d$C_IB_ngmL, na.rm = TRUE) < 0.01) {
      return(plotly_empty() %>% layout(title = "Ibrutinib not selected"))
    }
    plot_ly(d, x = ~time_days, y = ~C_IB_ngmL, type = "scatter", mode = "lines",
            line = list(color = "#1565C0", width = 2)) %>%
      layout(title = "Ibrutinib Concentration",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Concentration (ng/mL)"))
  })

  output$pk_ven_plot <- renderPlotly({
    d <- sim_data()
    if (max(d$C_VEN_ngmL, na.rm = TRUE) < 0.01) {
      return(plotly_empty() %>% layout(title = "Venetoclax not selected"))
    }
    plot_ly(d, x = ~time_days, y = ~C_VEN_ngmL, type = "scatter", mode = "lines",
            line = list(color = "#33691E", width = 2)) %>%
      layout(title = "Venetoclax Concentration",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Concentration (ng/mL)"))
  })

  output$pk_obi_plot <- renderPlotly({
    d <- sim_data()
    if (max(d$C_OBI_ugmL, na.rm = TRUE) < 0.01) {
      return(plotly_empty() %>% layout(title = "Obinutuzumab not selected"))
    }
    plot_ly(d, x = ~time_days, y = ~C_OBI_ugmL, type = "scatter", mode = "lines",
            line = list(color = "#BF360C", width = 2)) %>%
      layout(title = "Obinutuzumab Concentration",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Concentration (µg/mL)"))
  })

  # ── PD Biomarker plots ──────────────────────────────────────────────────────
  output$btk_occ_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~BTK_OCC_out, type = "scatter", mode = "lines",
            line = list(color = "#006064", width = 2)) %>%
      add_lines(y = 95, line = list(dash = "dash", color = "red", width = 1),
                name = "95% target") %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "BTK Occupancy (%)", range = c(0, 105)),
             showlegend = TRUE)
  })

  output$bcl2_occ_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~BCL2_OCC_out, type = "scatter", mode = "lines",
            line = list(color = "#33691E", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "BCL-2 Occupancy (%)", range = c(0, 105)))
  })

  output$cd20_occ_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~CD20_OCC_out, type = "scatter", mode = "lines",
            line = list(color = "#BF360C", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "CD20 Occupancy (%)", range = c(0, 105)))
  })

  output$adapt_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~MCL1_ADAPT, name = "MCL-1 (fold)",
            type = "scatter", mode = "lines",
            line = list(color = "#D32F2F", width = 2)) %>%
      add_lines(y = ~NK_ACT, name = "NK Activation (fold)",
                line = list(color = "#7B1FA2", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Fold change"), showlegend = TRUE)
  })

  # ── Clinical endpoint plots ─────────────────────────────────────────────────
  output$alc_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~ALC, type = "scatter", mode = "lines",
            line = list(color = "#1565C0", width = 2.5), name = "ALC") %>%
      add_lines(y = 4, line = list(dash = "dash", color = "grey", width = 1),
                name = "IWCLL CR (4×10⁹/L)") %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "ALC (×10⁹/L)"), showlegend = TRUE)
  })

  output$burden_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days) %>%
      add_lines(y = ~BM_CLL, name = "BM (%)",  line = list(color = "#E65100")) %>%
      add_lines(y = ~LN_CLL, name = "LN (rel%)", line = list(color = "#7B1FA2")) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Burden (%)"), showlegend = TRUE)
  })

  output$composite_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~BURDEN, type = "scatter", mode = "lines",
            line = list(color = "#263238", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Composite burden (%)"))
  })

  output$response_table <- renderTable({
    d <- sim_data()
    data.frame(
      Endpoint  = c("ALC Nadir (×10⁹/L)", "Nadir Day", "CR Achieved", "MRD-neg"),
      Value     = c(round(min(d$ALC, na.rm = TRUE), 2),
                    round(d$time_days[which.min(d$ALC)], 0),
                    ifelse(any(d$CR_flag == 1), "Yes", "No"),
                    ifelse(any(d$MRD_neg == 1), "Yes", "No")),
      stringsAsFactors = FALSE
    )
  })

  output$cr_time <- renderText({
    d <- sim_data()
    first_cr <- d$time_days[d$CR_flag == 1]
    if (length(first_cr) > 0) {
      paste0("CR first achieved at day ", round(min(first_cr), 0))
    } else {
      "CR not achieved in simulation period"
    }
  })

  # ── Scenario comparison ─────────────────────────────────────────────────────
  compare_data <- eventReactive(input$run_compare, {
    scen_configs <- list(
      ib      = list(use_ib=TRUE,  use_ven=FALSE, use_obi=FALSE, label="Ibrutinib"),
      ven     = list(use_ib=FALSE, use_ven=TRUE,  use_obi=FALSE, label="Venetoclax"),
      obi     = list(use_ib=FALSE, use_ven=FALSE, use_obi=TRUE,  label="Obinutuzumab"),
      ven_obi = list(use_ib=FALSE, use_ven=TRUE,  use_obi=TRUE,  label="VEN+OBI"),
      ib_ven  = list(use_ib=TRUE,  use_ven=TRUE,  use_obi=FALSE, label="IB+VEN"),
      triplet = list(use_ib=TRUE,  use_ven=TRUE,  use_obi=TRUE,  label="IB+VEN+OBI")
    )
    withProgress(message = "Running all scenarios...", {
      bind_rows(lapply(input$compare_scen, function(s) {
        cfg <- scen_configs[[s]]
        simulate_cll(input$alc0, input$bm0, input$ln0,
                     cfg$use_ib, cfg$use_ven, cfg$use_obi,
                     420, 400, input$end_days) %>%
          mutate(scenario = cfg$label)
      }))
    })
  })

  pal6 <- c("#1565C0","#33691E","#BF360C","#006064","#F57F17","#7B1FA2")

  output$cmp_alc <- renderPlotly({
    d <- compare_data()
    plot_ly(d, x = ~time_days, y = ~ALC, color = ~scenario,
            colors = pal6, type = "scatter", mode = "lines") %>%
      add_lines(y = 4, line = list(dash = "dash", color = "grey"), showlegend = FALSE) %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "ALC (×10⁹/L)"),
             legend = list(orientation = "h"))
  })

  output$cmp_burden <- renderPlotly({
    d <- compare_data()
    plot_ly(d, x = ~time_days, y = ~BURDEN, color = ~scenario,
            colors = pal6, type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Composite burden (%)"),
             legend = list(orientation = "h"))
  })

  output$cmp_btk <- renderPlotly({
    d <- compare_data()
    plot_ly(d, x = ~time_days, y = ~BTK_OCC_out, color = ~scenario,
            colors = pal6, type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "BTK Occupancy (%)"),
             legend = list(orientation = "h"))
  })

  output$cmp_table <- renderDT({
    d <- compare_data()
    tbl <- d %>%
      group_by(scenario) %>%
      summarise(
        ALC_nadir   = round(min(ALC, na.rm = TRUE), 2),
        Nadir_day   = round(time_days[which.min(ALC)], 0),
        CR          = ifelse(any(CR_flag == 1), "Yes", "No"),
        MRD_neg     = ifelse(any(MRD_neg == 1), "Yes", "No"),
        BM_final    = round(last(BM_CLL), 1),
        MCL1_peak   = round(max(MCL1_ADAPT, na.rm = TRUE), 2),
        .groups = "drop"
      )
    datatable(tbl, options = list(pageLength = 10, dom = "t"), rownames = FALSE)
  })

  # ── Genetics tab ────────────────────────────────────────────────────────────
  output$genetics_table <- renderDT({
    tbl <- data.frame(
      Marker        = c("del(13q14)", "Normal karyotype", "Trisomy 12",
                        "del(11q22)", "del(17p13)/TP53", "IGHV mutated",
                        "IGHV unmutated", "NOTCH1 mut", "SF3B1 mut",
                        "BTK C481S", "BCL-2 G101V"),
      Frequency     = c("55%","~18%","16%","12%","7-10%","40-50%","50-60%",
                        "10%","10%","acquired","acquired"),
      Risk          = c("Low","Intermediate","Intermediate","High","Very High",
                        "Favorable","Unfavorable","Poor","Poor","Resistance","Resistance"),
      Implication   = c("miR-15a/16-1 loss → BCL-2 overexpression",
                        "Standard prognosis",
                        "BCR signaling activated; rituximab response",
                        "ATM deletion; DNA damage checkpoint impaired",
                        "p53 pathway inactivated; chemo-resistant; BTKi/VEN preferred",
                        "Good response to all therapies; longer PFS",
                        "BCR-driven; BTKi highly effective",
                        "NOTCH pathway activated; chemotherapy less effective",
                        "Aberrant splicing; reduced apoptosis; VEN preferred",
                        "Ibrutinib resistance; switch to pirtobrutinib (non-covalent)",
                        "Venetoclax resistance; combination with BTKi"),
      stringsAsFactors = FALSE
    )
    datatable(tbl, options = list(pageLength = 15, dom = "t"), rownames = FALSE)
  })

  del17p_sim <- eventReactive(input$run_del17p, {
    # del(17p) reduces p53 → less apoptosis (modeled as Emax_BCL2 reduction)
    alc0 <- input$alc0; bm0 <- input$bm0; ln0 <- input$ln0
    end_d <- input$end_days
    # Normal
    d_norm <- simulate_cll(alc0, bm0, ln0, TRUE, TRUE, FALSE, 420, 400, end_d) %>%
      mutate(group = "No del(17p)")
    # del(17p): impaired p53 → reduced apoptosis under venetoclax
    mod_17p <- param(mod, list(Emax_BCL2 = if (input$del17p_flag == "Yes") 0.40 else 0.90,
                                Emax_BTK  = if (input$del17p_flag == "Yes") 0.65 else 0.70))
    e_17p <- c(ev(amt=420,cmt="DEPOT_IB",time=0,ii=24,addl=end_d-1),
               ev(amt=400,cmt="DEPOT_VEN",time=0,ii=24,addl=end_d-1))
    d_17p <- mrgsim(init(param(mod, list(ALC_0=alc0,BM_0=bm0,LN_0=ln0,
                                         Emax_BCL2=if(input$del17p_flag=="Yes")0.40 else 0.90)),
                         list(ALC=alc0,BM_CLL=bm0,LN_CLL=ln0)),
                    ev=e_17p, start=0, end=end_d*24, delta=12) %>%
      as_tibble() %>%
      mutate(time_days=time/24,
             group=if(input$del17p_flag=="Yes")"With del(17p)" else "No del(17p)")
    bind_rows(d_norm, d_17p)
  })

  output$del17p_plot <- renderPlotly({
    d <- del17p_sim()
    plot_ly(d, x = ~time_days, y = ~ALC, color = ~group,
            type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "ALC (×10⁹/L)"))
  })

  ighv_sim_data <- eventReactive(input$run_ighv, {
    alc0 <- input$alc0; bm0 <- input$bm0; ln0 <- input$ln0
    end_d <- input$end_days
    # IGHV unmutated: stronger BCR → higher kprol_CLL
    kprol_val <- if (input$ighv_sim == "unmut") 0.0045 else 0.0020
    e_ib <- ev(amt = 420, cmt = "DEPOT_IB", time = 0, ii = 24, addl = end_d - 1)
    d <- mrgsim(
      init(param(mod, list(ALC_0=alc0, BM_0=bm0, LN_0=ln0, kprol_CLL=kprol_val)),
           list(ALC=alc0, BM_CLL=bm0, LN_CLL=ln0)),
      ev = e_ib, start = 0, end = end_d * 24, delta = 12
    ) %>%
      as_tibble() %>%
      mutate(time_days = time / 24,
             IGHV = if (input$ighv_sim == "unmut") "IGHV Unmutated" else "IGHV Mutated")
    d
  })

  output$ighv_plot <- renderPlotly({
    d <- ighv_sim_data()
    plot_ly(d, x = ~time_days, y = ~ALC, color = ~IGHV,
            type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "ALC (×10⁹/L)"))
  })
}

# ── Launch ────────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
