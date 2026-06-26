################################################################################
# Kawasaki Disease QSP — Shiny Interactive Dashboard
# Tabs: Patient Profile · PK · Cytokines/Inflammation · Clinical Endpoints ·
#       Scenario Comparison · Biomarkers & Risk
################################################################################

library(shiny)
library(shinydashboard)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)
library(mrgsolve)

# ---- Inline mrgsolve model (minimal version for Shiny) ----
kd_code_shiny <- '
$PARAM
CL_IVIG=0.0033 Vc_IVIG=0.050 Vp_IVIG=0.048 Qp_IVIG=0.010
ka_FcRn=0.015 F_FcRn=0.60
ka_ASA=0.80 Vc_ASA=0.18 CL_ASA=0.60 f_SA=0.85 CL_SA=0.010 Vc_SA=0.16
CL_MP=0.48 Vc_MP=0.70 Vp_MP=0.55 Qp_MP=0.060
CL_IFX=0.0062 Vc_IFX=0.052 Vp_IFX=0.046 Qp_IFX=0.0038
CL_ANK=0.28 Vc_ANK=0.17 ka_ANK=0.30
kprod_IL1=0.12 kdeg_IL1=0.15 ksyn_IL1=0.25 IL1_base=0.8
kprod_IL6=0.28 kdeg_IL6=0.10 ksyn_IL6=0.40 IL6_base=2.5
kprod_TNF=0.15 kdeg_TNF=0.20 ksyn_TNF=0.30 TNF_base=1.2
Mac0=1.0 kact_Mac=0.08 kinact_Mac=0.05
kact_EC=0.10 kinact_EC=0.04 EC0=1.0
T_base=36.8 krise_T=0.35 kfall_T=0.10 T_max=40.5
kprod_CRP=3.5 kdeg_CRP=0.040 CRP_base=1.0
PLT_base=250 kprod_PLT=0.12 kdeg_PLT=0.006
Z0=0.0 krise_Z=0.0015 kfall_Z=0.0008 Z_max=12.0
EC50_IVIG=8.0 Emax_IVIG=0.80 n_IVIG=1.5
EC50_ASA=15.0 Emax_ASA=0.60 IC50_COX1=3.0 Emax_COX1=0.95
EC50_MP=0.25 Emax_MP=0.75
EC50_IFX=2.5 Emax_IFX=0.90 n_IFX=1.8
EC50_ANK=1.0 Emax_ANK=0.85
WT=15 IVIG_RES=0

$CMT A_IVIG_c A_IVIG_p A_ASA_gut A_ASA_c A_SA_c A_MP_c A_MP_p
     A_IFX_c A_IFX_p A_ANK_gut A_ANK_c
     IL1b IL6 TNFa Mac_act EC_act Fever CRP PLT_c CAL_Z

$INIT A_IVIG_c=0 A_IVIG_p=0 A_ASA_gut=0 A_ASA_c=0 A_SA_c=0
      A_MP_c=0 A_MP_p=0 A_IFX_c=0 A_IFX_p=0 A_ANK_gut=0 A_ANK_c=0
      IL1b=4.8 IL6=20.0 TNFa=6.0 Mac_act=3.5 EC_act=2.5
      Fever=39.5 CRP=45.0 PLT_c=250 CAL_Z=0.0

$ODE
double C_IVIG  = A_IVIG_c/(Vc_IVIG*WT);
double C_ASA   = A_ASA_c/(Vc_ASA*WT);
double C_SA    = A_SA_c/(Vc_SA*WT);
double C_MP    = A_MP_c/(Vc_MP*WT);
double C_IFX_u = A_IFX_c/(Vc_IFX*WT)*1000.0;
double C_ANK_u = A_ANK_c/(Vc_ANK*WT)*1000.0;
double E_IVIG  = Emax_IVIG*pow(C_IVIG,n_IVIG)/(pow(EC50_IVIG,n_IVIG)+pow(C_IVIG,n_IVIG));
if(IVIG_RES==1) E_IVIG=E_IVIG*0.30;
double E_ASA_p = Emax_ASA*C_SA/(EC50_ASA+C_SA);
double E_MP    = Emax_MP*C_MP/(EC50_MP+C_MP);
double E_IFX   = Emax_IFX*pow(C_IFX_u,n_IFX)/(pow(EC50_IFX,n_IFX)+pow(C_IFX_u,n_IFX));
double E_ANK   = Emax_ANK*C_ANK_u/(EC50_ANK+C_ANK_u);
double E_tot   = 1.0-(1.0-E_IVIG)*(1.0-E_MP)*(1.0-E_ANK);
double E_TNFd  = 1.0-(1.0-E_IFX)*(1.0-E_MP);
dxdt_A_IVIG_c  = -(CL_IVIG*WT)*C_IVIG-(Qp_IVIG*WT)*C_IVIG+(Qp_IVIG*WT)*(A_IVIG_p/(Vp_IVIG*WT))+F_FcRn*ka_FcRn*A_IVIG_c;
dxdt_A_IVIG_p  =  (Qp_IVIG*WT)*C_IVIG-(Qp_IVIG*WT)*(A_IVIG_p/(Vp_IVIG*WT));
dxdt_A_ASA_gut = -ka_ASA*A_ASA_gut;
dxdt_A_ASA_c   =  ka_ASA*A_ASA_gut-(CL_ASA*WT)*C_ASA;
dxdt_A_SA_c    =  f_SA*(CL_ASA*WT)*C_ASA-(CL_SA*WT)*C_SA;
dxdt_A_MP_c    = -(CL_MP*WT)*C_MP-(Qp_MP*WT)*C_MP+(Qp_MP*WT)*(A_MP_p/(Vp_MP*WT));
dxdt_A_MP_p    =  (Qp_MP*WT)*C_MP-(Qp_MP*WT)*(A_MP_p/(Vp_MP*WT));
dxdt_A_IFX_c   = -(CL_IFX*WT)*A_IFX_c/(Vc_IFX*WT)*(Vc_IFX*WT)-(Qp_IFX*WT)*C_IFX_u/1000.0+(Qp_IFX*WT)*(A_IFX_p/(Vp_IFX*WT));
dxdt_A_IFX_p   =  (Qp_IFX*WT)*A_IFX_c/(Vc_IFX*WT)-(Qp_IFX*WT)*(A_IFX_p/(Vp_IFX*WT));
dxdt_A_ANK_gut = -ka_ANK*A_ANK_gut;
dxdt_A_ANK_c   =  ka_ANK*A_ANK_gut-(CL_ANK*WT)*C_ANK_u/1000.0;
dxdt_Mac_act   =  kact_Mac*(1.0-E_tot)*Mac_act*(1.0-Mac_act/10.0)-kinact_Mac*Mac_act;
dxdt_IL1b      =  kprod_IL1+ksyn_IL1*Mac_act*(1.0-E_tot)*(1.0-E_ANK)*IL1_base-kdeg_IL1*IL1b;
dxdt_IL6       =  kprod_IL6+ksyn_IL6*(IL1b/IL1_base)*Mac_act*(1.0-E_tot)*IL6_base-kdeg_IL6*IL6;
dxdt_TNFa      =  kprod_TNF+ksyn_TNF*Mac_act*(1.0-E_TNFd)*TNF_base-kdeg_TNF*TNFa;
double EC_stim = (TNFa/TNF_base+IL1b/IL1_base)/2.0;
dxdt_EC_act    =  kact_EC*EC_stim*(1.0-E_tot)-kinact_EC*EC_act;
double ap      =  E_ASA_p+E_MP*0.5;
dxdt_Fever     =  krise_T*(IL1b/IL1_base+TNFa/TNF_base+IL6/IL6_base)/3.0*(1.0-ap)-kfall_T*(Fever-T_base)*(1.0+E_IVIG);
dxdt_CRP       =  kprod_CRP*(IL6/IL6_base)*(1.0-E_tot)-kdeg_CRP*CRP;
dxdt_PLT_c     =  kprod_PLT*(IL6/IL6_base)-kdeg_PLT*PLT_c;
double dp      =  E_IVIG+E_MP*0.5+E_IFX*0.3; dp=(dp>1.0)?1.0:dp;
dxdt_CAL_Z     =  krise_Z*(EC_act-EC0)*(TNFa/TNF_base)*(1.0-dp)*(Z_max-CAL_Z)-kfall_Z*CAL_Z*E_IVIG;

$TABLE
capture IVIG_gL  = A_IVIG_c/(Vc_IVIG*WT);
capture ASA_mgL  = A_ASA_c/(Vc_ASA*WT);
capture SA_mgL   = A_SA_c/(Vc_SA*WT);
capture MP_mgL   = A_MP_c/(Vc_MP*WT);
capture IFX_ugmL = A_IFX_c/(Vc_IFX*WT)*1000.0;
capture ANK_ugmL = A_ANK_c/(Vc_ANK*WT)*1000.0;
capture TempC    = Fever;
capture CRP_mgL  = CRP;
capture PLT_k    = PLT_c;
capture Zscore   = CAL_Z;
capture IL1_ng   = IL1b;
capture IL6_ng   = IL6;
capture TNF_ng   = TNFa;
capture Mac      = Mac_act;
capture EC       = EC_act;
'

kd_mod <- mcode("kd_shiny", kd_code_shiny, quiet = TRUE)

# ---- Simulation helper ----
simulate_kd <- function(wt, ivig_dose, asa_hi, add_steroids, add_ifx, add_anakinra,
                        ivig_resistant, sim_days = 56) {
  # Build events
  events <- NULL

  # IVIG
  ivig_amt  <- ivig_dose * wt * 1000
  ivig_rate <- ivig_amt / 12
  ev_ivig   <- ev(amt = ivig_amt, cmt = "A_IVIG_c", time = 0, rate = ivig_rate)
  events    <- if(is.null(events)) ev_ivig else events + ev_ivig

  # ASA high dose
  asa_each <- asa_hi * wt / 4
  ev_asa   <- ev(amt = asa_each, cmt = "A_ASA_gut",
                 time = seq(0, by = 6, length.out = 4 * 14))
  events   <- events + ev_asa

  # ASA low dose
  asa_lo_each <- 5 * wt
  ev_lo    <- ev(amt = asa_lo_each, cmt = "A_ASA_gut",
                 time = seq(14 * 24, by = 24, length.out = 42))
  events   <- events + ev_lo

  # Methylprednisolone
  if(add_steroids) {
    ev_mp <- ev(amt = 2 * wt, cmt = "A_MP_c",
                time = seq(0, by = 24, length.out = 5))
    events <- events + ev_mp
  }

  # Infliximab
  if(add_ifx) {
    ev_ifx <- ev(amt = 5 * wt, cmt = "A_IFX_c", time = 48,
                 rate = 5 * wt / 2)
    events <- events + ev_ifx
  }

  # Anakinra
  if(add_anakinra) {
    ev_ank <- ev(amt = 4 * wt, cmt = "A_ANK_gut",
                 time = seq(0, by = 24, length.out = 14))
    events <- events + ev_ank
  }

  res <- mrgsim(kd_mod, events = events,
                end = sim_days * 24, delta = 1,
                param = list(WT = wt, IVIG_RES = as.integer(ivig_resistant)),
                quiet = TRUE)
  as.data.frame(res) %>% mutate(Day = time / 24)
}

# ==================================================================
# UI
# ==================================================================
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "Kawasaki Disease QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "patient",   icon = icon("child")),
      menuItem("Pharmacokinetics",    tabName = "pk",        icon = icon("flask")),
      menuItem("Cytokines/Inflam.",   tabName = "cytokine",  icon = icon("virus")),
      menuItem("Clinical Endpoints",  tabName = "clinical",  icon = icon("heartbeat")),
      menuItem("Scenario Comparison", tabName = "scenario",  icon = icon("chart-bar")),
      menuItem("Biomarkers & Risk",   tabName = "biomarker", icon = icon("dna"))
    ),

    hr(),
    h5("Patient Parameters", style = "color:#fff; padding-left:15px"),
    sliderInput("wt",      "Body Weight (kg)",      min = 5,  max = 40, value = 15, step = 1),
    sliderInput("ivig_dose","IVIG Dose (g/kg)",      min = 1.0,max = 2.5,value = 2.0,step = 0.5),
    sliderInput("asa_hi",  "High-dose ASA (mg/kg/d)",min=30,  max = 100,value = 80, step = 10),
    checkboxInput("add_steroids",  "Add Methylprednisolone",    value = FALSE),
    checkboxInput("add_ifx",       "Add Infliximab (rescue)",    value = FALSE),
    checkboxInput("add_anakinra",  "Add Anakinra (rescue)",      value = FALSE),
    checkboxInput("ivig_resistant","IVIG-Resistant Patient",     value = FALSE),
    sliderInput("sim_days", "Simulation Duration (days)", min = 14, max = 120, value = 56),
    actionButton("run_sim", "Run Simulation", class = "btn-warning btn-block")
  ),

  dashboardBody(
    tabItems(

      # ---- Tab 1: Patient Profile ----
      tabItem("patient",
        fluidRow(
          box(width = 12, title = "Kawasaki Disease — Disease Overview",
              status = "danger", solidHeader = TRUE,
              fluidRow(
                column(6,
                  h4("Disease Biology"),
                  p("Kawasaki Disease (KD) is an acute self-limited febrile vasculitis of
                    childhood, predominantly affecting children under 5 years of age.
                    It is the leading cause of acquired heart disease in children in
                    developed countries, primarily through coronary artery aneurysm (CAA) formation."),
                  p(strong("Key pathophysiology:"), "Unknown trigger (possibly viral/bacterial)
                    → innate immune activation → NLRP3 inflammasome → cytokine storm
                    (IL-1β, IL-6, TNF-α) → endothelial activation → coronary arteritis
                    → aneurysm formation"),
                  h4("Diagnostic Criteria (AHA 2017)"),
                  tags$ul(
                    tags$li("Fever ≥ 38.5°C for ≥ 5 days"),
                    tags$li("Bilateral non-purulent conjunctivitis"),
                    tags$li("Polymorphous rash"),
                    tags$li("Oral changes (strawberry tongue, cracked lips)"),
                    tags$li("Extremity changes (erythema, desquamation)"),
                    tags$li("Cervical lymphadenopathy")
                  )
                ),
                column(6,
                  h4("Treatment Algorithm"),
                  tags$ol(
                    tags$li("Standard: IVIG 2 g/kg + High-dose aspirin 80–100 mg/kg/day → Low-dose 3–5 mg/kg/day"),
                    tags$li("Kobayashi high-risk: + Prednisolone/Methylprednisolone"),
                    tags$li("IVIG-resistant (~10–20%): 2nd IVIG OR Infliximab OR Steroids OR Cyclosporine"),
                    tags$li("Refractory: Anakinra (IL-1 blockade) / Infliximab (TNF-α blockade)")
                  ),
                  h4("Coronary Artery Classification (AHA)"),
                  tableOutput("cal_table")
                )
              )
          )
        ),
        fluidRow(
          valueBoxOutput("vbox_fever", width = 3),
          valueBoxOutput("vbox_crp",   width = 3),
          valueBoxOutput("vbox_zscore",width = 3),
          valueBoxOutput("vbox_plt",   width = 3)
        )
      ),

      # ---- Tab 2: PK ----
      tabItem("pk",
        fluidRow(
          box(width = 12, title = "Pharmacokinetics", status = "primary",
              solidHeader = TRUE,
              tabsetPanel(
                tabPanel("IVIG",   plotlyOutput("plot_ivig_pk",  height = 380)),
                tabPanel("Aspirin/Salicylate", plotlyOutput("plot_asa_pk", height = 380)),
                tabPanel("Methylprednisolone", plotlyOutput("plot_mp_pk",  height = 380)),
                tabPanel("Infliximab",  plotlyOutput("plot_ifx_pk",  height = 380)),
                tabPanel("Anakinra",    plotlyOutput("plot_ank_pk",  height = 380))
              )
          )
        ),
        fluidRow(
          box(width = 12, title = "PK Summary Table", status = "primary",
              DTOutput("pk_summary_tbl"))
        )
      ),

      # ---- Tab 3: Cytokines ----
      tabItem("cytokine",
        fluidRow(
          box(width = 6, title = "IL-1β Dynamics", status = "warning",
              plotlyOutput("plot_il1", height = 300)),
          box(width = 6, title = "IL-6 Dynamics", status = "warning",
              plotlyOutput("plot_il6", height = 300))
        ),
        fluidRow(
          box(width = 6, title = "TNF-α Dynamics", status = "danger",
              plotlyOutput("plot_tnf", height = 300)),
          box(width = 6, title = "Macrophage & Endothelial Activation", status = "info",
              plotlyOutput("plot_mac_ec", height = 300))
        )
      ),

      # ---- Tab 4: Clinical Endpoints ----
      tabItem("clinical",
        fluidRow(
          box(width = 6, title = "Body Temperature (Fever)", status = "danger",
              plotlyOutput("plot_fever", height = 320)),
          box(width = 6, title = "C-Reactive Protein", status = "warning",
              plotlyOutput("plot_crp", height = 320))
        ),
        fluidRow(
          box(width = 6, title = "Coronary Artery Z-score", status = "danger",
              plotlyOutput("plot_zscore", height = 320)),
          box(width = 6, title = "Platelet Count (Thrombocytosis)", status = "info",
              plotlyOutput("plot_plt", height = 320))
        )
      ),

      # ---- Tab 5: Scenario Comparison ----
      tabItem("scenario",
        fluidRow(
          box(width = 12, title = "5-Scenario Treatment Comparison",
              status = "primary", solidHeader = TRUE,
              p("Comparing: S1=IVIG+ASA, S2=+Steroids, S3=2nd IVIG(resistant),
                S4=Infliximab rescue, S5=Anakinra rescue"),
              tabsetPanel(
                tabPanel("Fever",   plotlyOutput("sc_fever",  height = 380)),
                tabPanel("CRP",     plotlyOutput("sc_crp",    height = 380)),
                tabPanel("Z-score", plotlyOutput("sc_zscore", height = 380)),
                tabPanel("Cytokines", plotlyOutput("sc_cyt",  height = 380))
              )
          )
        ),
        fluidRow(
          box(width = 12, title = "Scenario Summary (Day 7 & Peak Values)",
              status = "primary", DTOutput("sc_summary_tbl"))
        )
      ),

      # ---- Tab 6: Biomarkers & Risk ----
      tabItem("biomarker",
        fluidRow(
          box(width = 6, title = "Kobayashi Risk Score Components",
              status = "warning",
              sliderInput("kob_na",    "Sodium (mEq/L)", 125, 145, 132),
              sliderInput("kob_crp",   "CRP (mg/dL)",    0,   20,   8),
              sliderInput("kob_alb",   "Albumin (g/dL)", 2.0, 4.5,  3.0, step = 0.1),
              sliderInput("kob_alt",   "ALT (IU/L)",     0,   400, 80),
              sliderInput("kob_days",  "Days of Fever",  0,   10,  4),
              sliderInput("kob_age",   "Age (months)",   0,   60,  12),
              numericInput("kob_plt",  "Platelet (×10³/μL)", value = 280),
              verbatimTextOutput("kob_result")
          ),
          box(width = 6, title = "CAL Risk & IVIG Response Probability",
              status = "danger",
              plotlyOutput("biomarker_risk", height = 350),
              verbatimTextOutput("ivig_resist_prob")
          )
        ),
        fluidRow(
          box(width = 12, title = "Z-score Trajectory with CAL Thresholds",
              status = "danger",
              plotlyOutput("zscore_detail", height = 380))
        )
      )
    )
  )
)

# ==================================================================
# Server
# ==================================================================
server <- function(input, output, session) {

  # ---- Reactive simulation ----
  sim_result <- eventReactive(input$run_sim, {
    withProgress(message = "Running KD QSP simulation...", value = 0.5, {
      simulate_kd(
        wt            = input$wt,
        ivig_dose     = input$ivig_dose,
        asa_hi        = input$asa_hi,
        add_steroids  = input$add_steroids,
        add_ifx       = input$add_ifx,
        add_anakinra  = input$add_anakinra,
        ivig_resistant= input$ivig_resistant,
        sim_days      = input$sim_days
      )
    })
  }, ignoreNULL = FALSE)

  # ---- CAL classification table ----
  output$cal_table <- renderTable({
    data.frame(
      Category = c("No involvement", "Dilation only", "Small CAA",
                   "Medium CAA", "Large/Giant CAA"),
      `Z-score` = c("< 2.0", "2.0–2.49", "2.5–< 5",
                    "5–< 10", "≥ 10"),
      `Lumen diameter` = c("Normal", "< 3mm (<5yr)", "3–< 4mm",
                            "4–8mm", "> 8mm")
    )
  }, bordered = TRUE)

  # ---- Value boxes (Day 3 snapshot) ----
  get_day <- function(df, d) df[which.min(abs(df$Day - d)), ]

  output$vbox_fever <- renderValueBox({
    d <- get_day(sim_result(), 3)
    valueBox(paste0(round(d$TempC, 1), "°C"), "Temperature (Day 3)",
             icon = icon("thermometer-half"), color = "red")
  })
  output$vbox_crp <- renderValueBox({
    d <- get_day(sim_result(), 7)
    valueBox(paste0(round(d$CRP_mgL, 1), " mg/L"), "CRP (Day 7)",
             icon = icon("vial"), color = "orange")
  })
  output$vbox_zscore <- renderValueBox({
    z <- max(sim_result()$Zscore)
    col <- if(z >= 10) "red" else if(z >= 5) "orange" else if(z >= 2.5) "yellow" else "green"
    valueBox(round(z, 2), "Peak Coronary Z-score",
             icon = icon("heartbeat"), color = col)
  })
  output$vbox_plt <- renderValueBox({
    p <- max(sim_result()$PLT_k)
    valueBox(paste0(round(p, 0), "k"), "Peak Platelet Count",
             icon = icon("tint"), color = "purple")
  })

  # ---- PK plots ----
  mk_pk <- function(df, y, ylab, col) {
    plot_ly(df, x = ~Day, y = as.formula(paste0("~", y)),
            type = "scatter", mode = "lines",
            line = list(color = col, width = 2)) %>%
      layout(yaxis = list(title = ylab), xaxis = list(title = "Day"))
  }

  output$plot_ivig_pk <- renderPlotly({
    mk_pk(sim_result(), "IVIG_gL", "IVIG (g/L)", "#1565C0")
  })
  output$plot_asa_pk <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~Day) %>%
      add_trace(y = ~ASA_mgL, name = "Aspirin (mg/L)", mode = "lines",
                line = list(color = "#FFA000")) %>%
      add_trace(y = ~SA_mgL,  name = "Salicylate (mg/L)", mode = "lines",
                line = list(color = "#FF5722")) %>%
      layout(yaxis = list(title = "Concentration (mg/L)"),
             xaxis = list(title = "Day"))
  })
  output$plot_mp_pk <- renderPlotly({
    mk_pk(sim_result(), "MP_mgL", "Methylprednisolone (mg/L)", "#2E7D32")
  })
  output$plot_ifx_pk <- renderPlotly({
    mk_pk(sim_result(), "IFX_ugmL", "Infliximab (μg/mL)", "#6A1B9A")
  })
  output$plot_ank_pk <- renderPlotly({
    mk_pk(sim_result(), "ANK_ugmL", "Anakinra (μg/mL)", "#4527A0")
  })

  output$pk_summary_tbl <- renderDT({
    df <- sim_result()
    tbl <- data.frame(
      Drug       = c("IVIG","Aspirin","Salicylate","Methylprednisolone","Infliximab","Anakinra"),
      Cmax_unit  = c("g/L","mg/L","mg/L","mg/L","μg/mL","μg/mL"),
      Cmax       = c(max(df$IVIG_gL), max(df$ASA_mgL), max(df$SA_mgL),
                     max(df$MP_mgL),  max(df$IFX_ugmL), max(df$ANK_ugmL)),
      Tmax_d     = c(which.max(df$IVIG_gL), which.max(df$ASA_mgL),
                     which.max(df$SA_mgL),  which.max(df$MP_mgL),
                     which.max(df$IFX_ugmL),which.max(df$ANK_ugmL)) / 24
    )
    tbl$Cmax <- round(tbl$Cmax, 3)
    tbl$Tmax_d <- round(tbl$Tmax_d, 1)
    datatable(tbl, options = list(pageLength = 10, dom = "t"))
  })

  # ---- Cytokine plots ----
  output$plot_il1 <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~Day, y = ~IL1_ng, type = "scatter", mode = "lines",
            line = list(color = "#E91E63", width = 2)) %>%
      add_segments(x = 0, xend = max(df$Day), y = 0.8, yend = 0.8,
                   line = list(dash = "dash", color = "gray")) %>%
      layout(yaxis = list(title = "IL-1β (ng/mL)"))
  })
  output$plot_il6 <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~Day, y = ~IL6_ng, type = "scatter", mode = "lines",
            line = list(color = "#1B5E20", width = 2)) %>%
      layout(yaxis = list(title = "IL-6 (ng/mL)"))
  })
  output$plot_tnf <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~Day, y = ~TNF_ng, type = "scatter", mode = "lines",
            line = list(color = "#B71C1C", width = 2)) %>%
      layout(yaxis = list(title = "TNF-α (ng/mL)"))
  })
  output$plot_mac_ec <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~Day) %>%
      add_trace(y = ~Mac, name = "Macrophage Activation", mode = "lines",
                line = list(color = "#1565C0")) %>%
      add_trace(y = ~EC,  name = "Endothelial Activation", mode = "lines",
                line = list(color = "#FF8F00")) %>%
      layout(yaxis = list(title = "Activation (fold-baseline)"))
  })

  # ---- Clinical endpoint plots ----
  output$plot_fever <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~Day, y = ~TempC, type = "scatter", mode = "lines",
            line = list(color = "#C62828", width = 2.5)) %>%
      add_segments(x = 0, xend = max(df$Day), y = 38.5, yend = 38.5,
                   line = list(dash = "dash", color = "red")) %>%
      layout(yaxis = list(title = "Temperature (°C)", range = c(35, 42)))
  })
  output$plot_crp <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~Day, y = ~CRP_mgL, type = "scatter", mode = "lines",
            line = list(color = "#E65100", width = 2)) %>%
      add_segments(x = 0, xend = max(df$Day), y = 10, yend = 10,
                   line = list(dash = "dash", color = "gray")) %>%
      layout(yaxis = list(title = "CRP (mg/L)"))
  })
  output$plot_zscore <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~Day, y = ~Zscore, type = "scatter", mode = "lines",
            line = list(color = "#B71C1C", width = 2.5)) %>%
      add_segments(x = 0, xend = max(df$Day), y = 2.5, yend = 2.5,
                   line = list(dash = "dash", color = "orange")) %>%
      add_segments(x = 0, xend = max(df$Day), y = 10, yend = 10,
                   line = list(dash = "dash", color = "red")) %>%
      layout(yaxis = list(title = "Coronary Z-score"),
             annotations = list(
               list(x = 5, y = 2.5, text = "CAL threshold", showarrow = FALSE),
               list(x = 5, y = 10,  text = "Giant CAA",    showarrow = FALSE)
             ))
  })
  output$plot_plt <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~Day, y = ~PLT_k, type = "scatter", mode = "lines",
            line = list(color = "#6A1B9A", width = 2)) %>%
      add_segments(x = 0, xend = max(df$Day), y = 500, yend = 500,
                   line = list(dash = "dash", color = "purple")) %>%
      layout(yaxis = list(title = "Platelet (×10³/μL)"))
  })

  # ---- Scenario comparison ----
  sc_data <- reactive({
    wt <- input$wt
    bind_rows(
      simulate_kd(wt, 2.0, 80, FALSE, FALSE, FALSE, FALSE) %>%
        mutate(Scenario = "S1: IVIG+ASA"),
      simulate_kd(wt, 2.0, 80, TRUE,  FALSE, FALSE, FALSE) %>%
        mutate(Scenario = "S2: +Steroids"),
      simulate_kd(wt, 2.0, 80, FALSE, FALSE, FALSE, TRUE)  %>%
        mutate(Scenario = "S3: 2nd IVIG (resistant)"),
      simulate_kd(wt, 2.0, 80, FALSE, TRUE,  FALSE, TRUE)  %>%
        mutate(Scenario = "S4: Infliximab rescue"),
      simulate_kd(wt, 2.0, 80, FALSE, FALSE, TRUE,  TRUE)  %>%
        mutate(Scenario = "S5: Anakinra rescue")
    )
  })

  output$sc_fever <- renderPlotly({
    plot_ly(sc_data(), x = ~Day, y = ~TempC, color = ~Scenario,
            type = "scatter", mode = "lines") %>%
      add_segments(x = 0, xend = max(sc_data()$Day), y = 38.5, yend = 38.5,
                   line = list(dash = "dash", color = "red", width = 1)) %>%
      layout(yaxis = list(title = "Temperature (°C)"))
  })
  output$sc_crp <- renderPlotly({
    plot_ly(sc_data(), x = ~Day, y = ~CRP_mgL, color = ~Scenario,
            type = "scatter", mode = "lines") %>%
      layout(yaxis = list(title = "CRP (mg/L)"))
  })
  output$sc_zscore <- renderPlotly({
    plot_ly(sc_data(), x = ~Day, y = ~Zscore, color = ~Scenario,
            type = "scatter", mode = "lines") %>%
      layout(yaxis = list(title = "Coronary Z-score"))
  })
  output$sc_cyt <- renderPlotly({
    df_l <- sc_data() %>%
      select(Day, Scenario, IL1_ng, IL6_ng, TNF_ng) %>%
      pivot_longer(c(IL1_ng, IL6_ng, TNF_ng),
                   names_to = "Cytokine", values_to = "Conc")
    plot_ly(df_l, x = ~Day, y = ~Conc, color = ~Scenario,
            symbol = ~Cytokine, type = "scatter", mode = "lines") %>%
      layout(yaxis = list(title = "Cytokine (ng/mL)"))
  })

  output$sc_summary_tbl <- renderDT({
    sc <- sc_data() %>%
      group_by(Scenario) %>%
      summarise(
        Fever_Day3   = round(TempC[which.min(abs(Day - 3))], 1),
        Fever_Day7   = round(TempC[which.min(abs(Day - 7))], 1),
        CRP_Day7     = round(CRP_mgL[which.min(abs(Day - 7))], 1),
        CRP_Peak     = round(max(CRP_mgL), 1),
        Zscore_Peak  = round(max(Zscore), 2),
        PLT_Peak_k   = round(max(PLT_k), 0),
        .groups = "drop"
      )
    datatable(sc, options = list(dom = "t", pageLength = 10))
  })

  # ---- Biomarker / Kobayashi score ----
  kobayashi_score <- reactive({
    score <- 0
    if(input$kob_na <= 133)                score <- score + 2
    if(input$kob_crp >= 10)               score <- score + 1  # mg/dL
    if(input$kob_alb <= 3.0)              score <- score + 1
    if(input$kob_alt >= 100)              score <- score + 1
    if(input$kob_days <= 4)               score <- score + 2
    if(input$kob_age < 12)               score <- score + 1
    if(input$kob_plt <= 300)             score <- score + 1
    score
  })

  output$kob_result <- renderText({
    s <- kobayashi_score()
    risk <- if(s >= 4) "HIGH RISK (IVIG resistance ~50%)\n→ Consider primary steroid + IVIG" else
            if(s >= 2) "MODERATE RISK (IVIG resistance ~20%)" else
            "LOW RISK (IVIG resistance ~5%)"
    paste0("Kobayashi Score: ", s, " / 10\nRisk: ", risk)
  })

  output$biomarker_risk <- renderPlotly({
    df <- sim_result()
    df_risk <- df %>%
      mutate(
        p_giant = pmin(1, Zscore / 15),
        inflam  = (IL1_ng/0.8 + IL6_ng/2.5 + TNF_ng/1.2) / 3
      )
    plot_ly(df_risk, x = ~Day) %>%
      add_trace(y = ~p_giant, name = "Giant CAA probability",
                type = "scatter", mode = "lines",
                line = list(color = "red")) %>%
      add_trace(y = ~inflam / max(df_risk$inflam),
                name = "Inflammation index (norm.)",
                type = "scatter", mode = "lines",
                line = list(color = "orange", dash = "dash")) %>%
      layout(yaxis = list(title = "Probability / Index"))
  })

  output$ivig_resist_prob <- renderText({
    df <- sim_result()
    day3 <- df[which.min(abs(df$Day - 3)), ]
    p_res <- if(input$ivig_resistant) 0.85 else
             1 / (1 + exp(-(0.8 * (day3$Mac - 3.0))))
    paste0("IVIG resistance probability (model estimate):\n",
           round(p_res * 100, 1), "%\n",
           "Macrophage activation at Day 3: ", round(day3$Mac, 2),
           "\nKobayashi score (sidebar): ", kobayashi_score())
  })

  output$zscore_detail <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x = ~Day, y = ~Zscore, type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(183,28,28,0.15)",
            line = list(color = "#B71C1C", width = 2.5),
            name = "Coronary Z-score") %>%
      add_segments(x = 0, xend = max(df$Day), y = 2.5, yend = 2.5,
                   line = list(dash = "dash", color = "#FF9800"), name = "Small CAA (z=2.5)") %>%
      add_segments(x = 0, xend = max(df$Day), y = 5, yend = 5,
                   line = list(dash = "dash", color = "#F44336"), name = "Medium CAA (z=5)") %>%
      add_segments(x = 0, xend = max(df$Day), y = 10, yend = 10,
                   line = list(dash = "dash", color = "#B71C1C"), name = "Giant CAA (z=10)") %>%
      layout(
        yaxis = list(title = "Coronary Artery Z-score", range = c(0, 13)),
        xaxis = list(title = "Day"),
        title = "Coronary Artery Z-score with AHA Classification Thresholds"
      )
  })
}

# Run app
shinyApp(ui = ui, server = server)
