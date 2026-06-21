##############################################################################
# GERD QSP Shiny App
# Interactive simulator: PPI / H2RA / P-CAB / Prokinetic pharmacology
# 6 tabs: Patient Profile · Drug PK · Gastric PD · Clinical Endpoints ·
#         Scenario Comparison · Biomarker / Risk Assessment
# Author: Claude Code Routine (CCR) — 2026-06-18
##############################################################################

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mrgsolve)
library(DT)
library(plotly)
library(scales)

# ─── mrgsolve model ──────────────────────────────────────────────────────────
gerd_code <- '
$PARAM
DOSE_PPI=20 DOSE_H2RA=40 DOSE_PCAB=20 DOSE_PROK=10
KA_PPI=0.80 CL_PPI=18 V_PPI=20 F_PPI=0.65 CYP2C19=1.0
KA_H2RA=1.50 CL_H2RA=14 V_H2RA=90 F_H2RA=0.50
KA_PCAB=2.00 CL_PCAB=8 V_PCAB=300 F_PCAB=0.60
KA_PROK=1.20 CL_PROK=60 V_PROK=450 F_PROK=0.15
PUMP_ACT0=30 K_ACT=0.30 K_DEACT=0.30 K_SYN_PUMP=3.0 K_DEG_PUMP=0.03
IC50_PPI=0.15 HILL_PPI=1.5 IC50_PCAB=0.08 IC50_H2RA=0.04 EMAX_H2RA=0.60
ACID_BASE=3.5 ACID_MAX=12.0 HIST_STIM=0.35 PH_BUFF=1.5 GASTRIC_EMP=0.50
LES_P0=22 TLESR_BASE=6 K_PROK_LES=5 K_PROK_EMP=0.30 REFLUX_K=0.05
MUC_HEAL0=0.30 MUC_INJ_K=0.40 K_HEAL=0.20 MAX_DAMAGE=100
SYM_K=0.05 SYM_HEAL_K=0.30

$CMT PPI_GUT PPI_CENT H2RA_GUT H2RA_CENT PCAB_GUT PCAB_CENT PROK_GUT PROK_CENT
     PUMP_INACT PUMP_ACT PUMP_INH ACID_RATE GAS_pH AET MUC_DMG MUC_HEAL SYM_SCORE BE_RISK

$MAIN
double CP_PPI  = PPI_CENT  / V_PPI;
double CP_H2RA = H2RA_CENT / V_H2RA;
double CP_PCAB = PCAB_CENT / V_PCAB;
double CP_PROK = PROK_CENT / V_PROK;
double INH_PPI  = pow(CP_PPI, HILL_PPI) / (pow(IC50_PPI, HILL_PPI) + pow(CP_PPI, HILL_PPI));
double INH_PCAB = CP_PCAB / (IC50_PCAB + CP_PCAB);
double INH_H2RA = EMAX_H2RA * CP_H2RA / (IC50_H2RA + CP_H2RA);
double ACID_INH = 1.0 - (1.0 - INH_PPI) * (1.0 - INH_PCAB) * (1.0 - INH_H2RA * HIST_STIM);
if(ACID_INH > 0.99) ACID_INH = 0.99;
double PROK_EFF_LES = K_PROK_LES * CP_PROK / (0.005 + CP_PROK);
double LES_P = LES_P0 + PROK_EFF_LES;
double TLESR = TLESR_BASE * (LES_P0 / LES_P);
double ACID_EFF = ACID_BASE + (ACID_MAX - ACID_BASE) * (1.0 - ACID_INH);
double pH_calc = 1.0 + PH_BUFF * exp(-ACID_EFF / 2.0);
if(pH_calc < 1.0) pH_calc = 1.0;
if(pH_calc > 7.0) pH_calc = 7.0;
double REFLUX_RATE = TLESR * REFLUX_K * (pH_calc < 4.0 ? 1.0 : exp(-pH_calc + 4.0));
double AET_INST = REFLUX_RATE * 100.0;
double DMG_RATE = MUC_INJ_K * (AET / 100.0) * (1.0 - MUC_HEAL / MAX_DAMAGE);
double HEAL_RATE = K_HEAL * (MUC_HEAL0 / 24.0) * (1.0 + ACID_INH);

$ODE
dxdt_PPI_GUT   = -KA_PPI * PPI_GUT;
dxdt_PPI_CENT  = KA_PPI * F_PPI * PPI_GUT * CYP2C19 - (CL_PPI/V_PPI)/CYP2C19 * PPI_CENT;
dxdt_H2RA_GUT  = -KA_H2RA * H2RA_GUT;
dxdt_H2RA_CENT = KA_H2RA * F_H2RA * H2RA_GUT - (CL_H2RA/V_H2RA) * H2RA_CENT;
dxdt_PCAB_GUT  = -KA_PCAB * PCAB_GUT;
dxdt_PCAB_CENT = KA_PCAB * F_PCAB * PCAB_GUT - (CL_PCAB/V_PCAB) * PCAB_CENT;
dxdt_PROK_GUT  = -KA_PROK * PROK_GUT;
dxdt_PROK_CENT = KA_PROK * F_PROK * PROK_GUT - (CL_PROK/V_PROK) * PROK_CENT;
dxdt_PUMP_INACT = K_SYN_PUMP - K_DEG_PUMP*PUMP_INACT - K_ACT*PUMP_INACT;
dxdt_PUMP_ACT   = K_ACT*PUMP_INACT - K_DEACT*PUMP_ACT - K_DEG_PUMP*PUMP_ACT - INH_PPI*0.5*PUMP_ACT;
dxdt_PUMP_INH   = INH_PPI*0.5*PUMP_ACT - K_DEG_PUMP*PUMP_INH;
dxdt_ACID_RATE  = 0.1*(ACID_EFF*(PUMP_ACT/PUMP_ACT0) - ACID_RATE);
dxdt_GAS_pH     = 0.5*(pH_calc - GAS_pH);
dxdt_AET        = 0.2*(AET_INST - AET);
dxdt_MUC_DMG    = DMG_RATE - HEAL_RATE;
dxdt_MUC_HEAL   = -DMG_RATE + HEAL_RATE;
dxdt_SYM_SCORE  = 0.1*(15.0*(AET/30.0)+3.0*(MUC_DMG/MAX_DAMAGE) - SYM_SCORE);
dxdt_BE_RISK    = 0.0001*AET*(1.0 - BE_RISK);

$TABLE
capture CP_PPI=PPI_CENT/V_PPI; capture CP_H2RA=H2RA_CENT/V_H2RA;
capture CP_PCAB=PCAB_CENT/V_PCAB; capture CP_PROK=PROK_CENT/V_PROK;
capture ACID_INH_pct=ACID_INH*100; capture pH=GAS_pH; capture AET_pct=AET;
capture DMG=MUC_DMG; capture HEAL=MUC_HEAL; capture SYM=SYM_SCORE;
capture LES_pressure=LES_P0+K_PROK_LES*CP_PROK/(0.005+CP_PROK);
capture Barrett=BE_RISK; capture TLESR_rate=TLESR_BASE*(LES_P0/(LES_P0+K_PROK_LES*CP_PROK/(0.005+CP_PROK)));

$INIT
PPI_GUT=0 PPI_CENT=0 H2RA_GUT=0 H2RA_CENT=0 PCAB_GUT=0 PCAB_CENT=0
PROK_GUT=0 PROK_CENT=0 PUMP_INACT=70 PUMP_ACT=30 PUMP_INH=0
ACID_RATE=3.5 GAS_pH=1.8 AET=15.0 MUC_DMG=25.0 MUC_HEAL=75.0
SYM_SCORE=8.0 BE_RISK=0.0
'

base_mod <- mcode("GERD_QSP_shiny", gerd_code, quiet = TRUE)

run_sim <- function(drug, dose, dose2 = 0, cyp = 1.0,
                    duration_wk = 8, interval = 24) {
  ev_list <- list()
  p_list  <- list(CYP2C19 = cyp)

  times <- seq(0, duration_wk * 7 * 24 - 1, by = interval)

  if (drug == "PPI") {
    ev_list <- ev(amt = dose, cmt = "PPI_GUT", time = times)
  } else if (drug == "H2RA") {
    t_bid <- sort(c(times, times + 12))
    ev_list <- ev(amt = dose, cmt = "H2RA_GUT", time = t_bid)
  } else if (drug == "P-CAB") {
    ev_list <- ev(amt = dose, cmt = "PCAB_GUT", time = times)
  } else if (drug == "PPI+Prokinetic") {
    t_tid <- sort(c(times, times + 8, times + 16))
    ev_list <- c(ev(amt = dose, cmt = "PPI_GUT", time = times),
                 ev(amt = dose2, cmt = "PROK_GUT", time = t_tid))
  } else {
    ev_list <- ev(amt = 0, cmt = "PPI_GUT", time = 0)
  }

  m2 <- param(base_mod, p_list)
  out <- mrgsim(m2, ev_list,
                start = 0, end = duration_wk * 7 * 24,
                delta = 0.5, carry_out = "evid")
  as.data.frame(out)
}

# ─── UI ───────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "GERD QSP Simulator"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",      tabName = "patient",    icon = icon("user")),
      menuItem("Drug PK",              tabName = "pk",         icon = icon("pills")),
      menuItem("Gastric PD",           tabName = "pd",         icon = icon("chart-line")),
      menuItem("Clinical Endpoints",   tabName = "clinical",   icon = icon("stethoscope")),
      menuItem("Scenario Comparison",  tabName = "compare",    icon = icon("table")),
      menuItem("Biomarker & Risk",     tabName = "biomarker",  icon = icon("dna"))
    ),

    hr(),
    h5("Drug Settings", style = "color:white; margin-left:10px;"),
    selectInput("drug", "Drug Class",
                choices = c("No Treatment", "PPI", "H2RA", "P-CAB", "PPI+Prokinetic"),
                selected = "PPI"),
    conditionalPanel("input.drug == 'PPI' || input.drug == 'PPI+Prokinetic'",
      sliderInput("ppi_dose", "PPI Dose (mg)", 10, 80, 20, step = 10)
    ),
    conditionalPanel("input.drug == 'H2RA'",
      sliderInput("h2ra_dose", "H2RA Dose (mg BID)", 20, 80, 40, step = 20)
    ),
    conditionalPanel("input.drug == 'P-CAB'",
      sliderInput("pcab_dose", "P-CAB Dose (mg)", 10, 40, 20, step = 10)
    ),
    conditionalPanel("input.drug == 'PPI+Prokinetic'",
      sliderInput("prok_dose", "Prokinetic Dose (mg TID)", 5, 20, 10, step = 5)
    ),

    hr(),
    h5("Patient Factors", style = "color:white; margin-left:10px;"),
    selectInput("cyp", "CYP2C19 Phenotype",
                choices = c("Ultra-Rapid (UM)"  = "2.0",
                            "Extensive (EM)"     = "1.0",
                            "Intermediate (IM)"  = "0.6",
                            "Poor (PM)"          = "0.25"),
                selected = "1.0"),
    sliderInput("duration", "Duration (weeks)", 4, 24, 8, step = 4),

    hr(),
    actionButton("run", "Run Simulation", icon = icon("play"),
                 class = "btn-success", width = "90%")
  ),

  dashboardBody(
    tabItems(

      # ── Tab 1: Patient Profile ────────────────────────────────────────────
      tabItem("patient",
        fluidRow(
          box(title = "GERD Overview", width = 6, solidHeader = TRUE, status = "primary",
            h4("Disease Characteristics"),
            p("Gastroesophageal Reflux Disease (GERD) affects ~20% of the Western population.
               The hallmark is transient lower esophageal sphincter relaxation (TLESR) allowing
               acid/bile to contact the esophageal mucosa."),
            h5("Key Metrics"),
            tableOutput("patient_table")
          ),
          box(title = "Pathophysiology Summary", width = 6, solidHeader = TRUE, status = "info",
            plotlyOutput("patho_radar", height = "350px")
          )
        ),
        fluidRow(
          box(title = "Risk Factor Impact on TLESR Rate", width = 12, solidHeader = TRUE, status = "warning",
            plotlyOutput("risk_bar", height = "250px")
          )
        )
      ),

      # ── Tab 2: Drug PK ───────────────────────────────────────────────────
      tabItem("pk",
        fluidRow(
          box(title = "Plasma Concentration-Time", width = 8, solidHeader = TRUE, status = "primary",
            plotlyOutput("pk_plot", height = "400px")
          ),
          box(title = "PK Parameters Summary", width = 4, solidHeader = TRUE, status = "info",
            DTOutput("pk_table")
          )
        ),
        fluidRow(
          box(title = "CYP2C19 Phenotype Effect on PPI Exposure", width = 12, solidHeader = TRUE, status = "warning",
            plotlyOutput("cyp_plot", height = "300px")
          )
        )
      ),

      # ── Tab 3: Gastric PD ────────────────────────────────────────────────
      tabItem("pd",
        fluidRow(
          box(title = "Intragastric pH (24h profile)", width = 6, solidHeader = TRUE, status = "primary",
            plotlyOutput("ph_plot", height = "350px")
          ),
          box(title = "H+/K+-ATPase Pool Dynamics", width = 6, solidHeader = TRUE, status = "info",
            plotlyOutput("pump_plot", height = "350px")
          )
        ),
        fluidRow(
          box(title = "Acid Secretion Rate & Inhibition %", width = 6, solidHeader = TRUE, status = "warning",
            plotlyOutput("acid_plot", height = "300px")
          ),
          box(title = "LES Pressure & TLESR Rate", width = 6, solidHeader = TRUE, status = "success",
            plotlyOutput("les_plot", height = "300px")
          )
        )
      ),

      # ── Tab 4: Clinical Endpoints ─────────────────────────────────────────
      tabItem("clinical",
        fluidRow(
          box(title = "Acid Exposure Time (%) over Treatment", width = 6, solidHeader = TRUE, status = "primary",
            plotlyOutput("aet_plot", height = "350px")
          ),
          box(title = "Mucosal Damage & Healing", width = 6, solidHeader = TRUE, status = "danger",
            plotlyOutput("mucosal_plot", height = "350px")
          )
        ),
        fluidRow(
          box(title = "Symptom Score (GERD-Q) Trajectory", width = 6, solidHeader = TRUE, status = "warning",
            plotlyOutput("sym_plot", height = "300px")
          ),
          box(title = "Endoscopic Healing Rate Prediction", width = 6, solidHeader = TRUE, status = "success",
            plotlyOutput("healing_rate_plot", height = "300px")
          )
        )
      ),

      # ── Tab 5: Scenario Comparison ────────────────────────────────────────
      tabItem("compare",
        fluidRow(
          box(title = "All Treatment Scenarios — AET (%)", width = 12, solidHeader = TRUE, status = "primary",
            plotlyOutput("compare_aet", height = "350px")
          )
        ),
        fluidRow(
          box(title = "Comparative Efficacy Table (Week 8)", width = 12, solidHeader = TRUE, status = "info",
            DTOutput("compare_table")
          )
        )
      ),

      # ── Tab 6: Biomarker & Risk ───────────────────────────────────────────
      tabItem("biomarker",
        fluidRow(
          box(title = "Barrett's Esophagus Risk Accumulation", width = 6, solidHeader = TRUE, status = "danger",
            plotlyOutput("barrett_plot", height = "350px")
          ),
          box(title = "Dose-Response: AET at Week 8", width = 6, solidHeader = TRUE, status = "primary",
            plotlyOutput("dr_plot", height = "350px")
          )
        ),
        fluidRow(
          box(title = "Biomarker Risk Classification", width = 12, solidHeader = TRUE, status = "warning",
            DTOutput("risk_table")
          )
        )
      )
    )
  )
)

# ─── Server ───────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  sim_data <- eventReactive(input$run, {
    req(input$drug)
    dose  <- switch(input$drug,
                    "PPI"           = input$ppi_dose,
                    "H2RA"          = input$h2ra_dose,
                    "P-CAB"         = input$pcab_dose,
                    "PPI+Prokinetic"= input$ppi_dose,
                    0)
    dose2 <- ifelse(input$drug == "PPI+Prokinetic", input$prok_dose, 0)
    run_sim(input$drug, dose, dose2,
            cyp      = as.numeric(input$cyp),
            duration = input$duration)
  }, ignoreNULL = FALSE)

  # ── Tab 1: Patient ─────────────────────────────────────────────────────────
  output$patient_table <- renderTable({
    data.frame(
      Parameter = c("Normal AET", "GERD threshold AET", "Normal LES pressure",
                    "TLESR rate (controls)", "TLESR rate (GERD)",
                    "Baseline symptom score"),
      Value     = c("< 6%", "≥ 6%", "10–40 mmHg",
                    "~4/h", "~6–8/h", "8–12 (GERD-Q)")
    )
  })

  output$patho_radar <- renderPlotly({
    categories <- c("Acid Secretion", "LES Dysfunction", "Mucosal Injury",
                    "Inflammation", "Dysmotility", "Neurogenic Pain")
    vals <- c(0.7, 0.8, 0.6, 0.5, 0.6, 0.7)
    plot_ly(
      type = "scatterpolar", fill = "toself",
      r = c(vals, vals[1]), theta = c(categories, categories[1]),
      name = "Typical GERD", fillcolor = "rgba(33,150,243,0.3)",
      line = list(color = "#1976D2")
    ) %>%
      layout(polar = list(radialaxis = list(visible = TRUE, range = c(0, 1))),
             title = "GERD Pathophysiology Profile")
  })

  output$risk_bar <- renderPlotly({
    df <- data.frame(
      Factor = c("Obesity (BMI>30)", "Hiatus Hernia", "Pregnancy",
                 "Smoking", "High-Fat Diet", "Alcohol"),
      TLESR_increase = c(2.1, 1.8, 1.6, 1.3, 1.4, 1.2)
    )
    plot_ly(df, x = ~Factor, y = ~TLESR_increase, type = "bar",
            marker = list(color = "#EF5350")) %>%
      layout(yaxis = list(title = "TLESR rate multiplier"),
             title = "Risk Factor Impact on TLESR Rate")
  })

  # ── Tab 2: PK ──────────────────────────────────────────────────────────────
  output$pk_plot <- renderPlotly({
    d <- sim_data() %>% filter(time <= 24)
    p <- ggplot(d, aes(x = time)) +
      geom_line(aes(y = CP_PPI,  colour = "PPI (mg/L)"),  linewidth = 1) +
      geom_line(aes(y = CP_H2RA, colour = "H2RA (mg/L)"), linewidth = 1) +
      geom_line(aes(y = CP_PCAB, colour = "P-CAB (mg/L)"), linewidth = 1) +
      scale_colour_manual(values = c("PPI (mg/L)" = "#8E24AA",
                                     "H2RA (mg/L)" = "#1E88E5",
                                     "P-CAB (mg/L)" = "#00897B")) +
      labs(x = "Time (h)", y = "Plasma Concentration (mg/L)",
           colour = "Drug", title = "Day 1 PK Profile") +
      theme_classic()
    ggplotly(p)
  })

  output$pk_table <- renderDT({
    d_last <- sim_data()
    d_day1 <- d_last %>% filter(time <= 24)
    df <- data.frame(
      Parameter   = c("PPI Cmax (mg/L)", "H2RA Cmax (mg/L)", "P-CAB Cmax (mg/L)",
                       "PPI Tmax (h)", "P-CAB Tmax (h)", "Acid Inhib (%)"),
      Value       = c(
        round(max(d_day1$CP_PPI),  3),
        round(max(d_day1$CP_H2RA), 3),
        round(max(d_day1$CP_PCAB), 3),
        round(d_day1$time[which.max(d_day1$CP_PPI)], 1),
        round(d_day1$time[which.max(d_day1$CP_PCAB)], 1),
        round(mean(tail(d_last$ACID_INH_pct, 48)), 1)
      )
    )
    datatable(df, options = list(dom = "t", paging = FALSE), rownames = FALSE)
  })

  output$cyp_plot <- renderPlotly({
    cyp_vals <- c("UM" = 2.0, "EM" = 1.0, "IM" = 0.6, "PM" = 0.25)
    cols     <- c("UM" = "#C62828", "EM" = "#1976D2", "IM" = "#388E3C", "PM" = "#F57F17")
    sims <- lapply(names(cyp_vals), function(ph) {
      d <- run_sim("PPI", 20, cyp = cyp_vals[ph], duration = 7)
      d$Phenotype <- ph
      d
    })
    d_all <- bind_rows(sims) %>% filter(time <= 24)
    p <- ggplot(d_all, aes(x = time, y = CP_PPI, colour = Phenotype)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = cols) +
      labs(x = "Time (h, Day 1)", y = "PPI Plasma Conc (mg/L)",
           title = "CYP2C19 Phenotype Effect") +
      theme_classic()
    ggplotly(p)
  })

  # ── Tab 3: PD ──────────────────────────────────────────────────────────────
  output$ph_plot <- renderPlotly({
    d <- sim_data() %>% filter(time <= 24)
    p <- ggplot(d, aes(x = time, y = pH)) +
      geom_line(colour = "#E53935", linewidth = 1) +
      geom_hline(yintercept = 4, linetype = "dashed", colour = "grey50") +
      labs(x = "Time (h, Day 1)", y = "Intragastric pH",
           title = "24-h Gastric pH Profile") +
      ylim(0, 7) + theme_classic()
    ggplotly(p)
  })

  output$pump_plot <- renderPlotly({
    d <- sim_data() %>% filter(time %% 12 == 0)
    d_long <- d %>%
      select(time, PUMP_ACT, PUMP_INH, PUMP_INACT) %>%
      pivot_longer(-time, names_to = "Pool", values_to = "Count") %>%
      mutate(Pool = recode(Pool,
        PUMP_ACT   = "Active",
        PUMP_INH   = "Inhibited (PPI)",
        PUMP_INACT = "Inactive"))
    p <- ggplot(d_long, aes(x = time / 24, y = Count, fill = Pool)) +
      geom_area(alpha = 0.7) +
      scale_fill_manual(values = c("Active" = "#EF5350",
                                   "Inhibited (PPI)" = "#7B1FA2",
                                   "Inactive" = "#90A4AE")) +
      labs(x = "Time (days)", y = "Pump Units",
           title = "H+/K+-ATPase Pool Dynamics") +
      theme_classic()
    ggplotly(p)
  })

  output$acid_plot <- renderPlotly({
    d <- sim_data() %>% filter(time %% 6 == 0)
    p <- ggplot(d, aes(x = time / 24)) +
      geom_line(aes(y = ACID_INH_pct), colour = "#8E24AA", linewidth = 1) +
      geom_hline(yintercept = 80, linetype = "dashed", colour = "grey50") +
      labs(x = "Time (days)", y = "Acid Inhibition (%)",
           title = "Acid Secretion Inhibition") +
      ylim(0, 100) + theme_classic()
    ggplotly(p)
  })

  output$les_plot <- renderPlotly({
    d <- sim_data() %>% filter(time %% 6 == 0)
    p <- ggplot(d, aes(x = time / 24)) +
      geom_line(aes(y = LES_pressure, colour = "LES Pressure (mmHg)"), linewidth = 1) +
      geom_line(aes(y = TLESR_rate * 2, colour = "TLESR Rate ×2 (events/h)"), linewidth = 1, linetype = "dashed") +
      scale_colour_manual(values = c("LES Pressure (mmHg)" = "#2E7D32",
                                     "TLESR Rate ×2 (events/h)" = "#C62828")) +
      labs(x = "Time (days)", y = "Value", colour = "Variable",
           title = "LES Pressure & TLESR Rate") +
      theme_classic()
    ggplotly(p)
  })

  # ── Tab 4: Clinical ────────────────────────────────────────────────────────
  output$aet_plot <- renderPlotly({
    d <- sim_data() %>% filter(time %% 6 == 0)
    p <- ggplot(d, aes(x = time / 24, y = AET_pct)) +
      geom_line(colour = "#E53935", linewidth = 1.2) +
      geom_hline(yintercept = 6, linetype = "dashed", colour = "#1565C0") +
      annotate("text", x = max(d$time / 24) * 0.8, y = 7,
               label = "Lyon 2.0 threshold: 6%", colour = "#1565C0", size = 3.5) +
      labs(x = "Time (days)", y = "AET (%)", title = "Esophageal Acid Exposure Time") +
      theme_classic()
    ggplotly(p)
  })

  output$mucosal_plot <- renderPlotly({
    d <- sim_data() %>% filter(time %% 6 == 0)
    p <- ggplot(d, aes(x = time / 24)) +
      geom_line(aes(y = DMG,  colour = "Damage Score"),  linewidth = 1) +
      geom_line(aes(y = HEAL, colour = "Integrity Index"), linewidth = 1) +
      scale_colour_manual(values = c("Damage Score" = "#C62828",
                                     "Integrity Index" = "#2E7D32")) +
      labs(x = "Time (days)", y = "Score (0-100)", colour = "Measure",
           title = "Mucosal Damage & Integrity") +
      ylim(0, 100) + theme_classic()
    ggplotly(p)
  })

  output$sym_plot <- renderPlotly({
    d <- sim_data() %>% filter(time %% 6 == 0)
    p <- ggplot(d, aes(x = time / 24, y = SYM)) +
      geom_line(colour = "#F57F17", linewidth = 1.2) +
      geom_hline(yintercept = 8, linetype = "dashed", colour = "#C62828") +
      annotate("text", x = max(d$time / 24) * 0.8, y = 9,
               label = "GERD-Q ≥8 = GERD+", colour = "#C62828", size = 3.5) +
      labs(x = "Time (days)", y = "GERD-Q Score",
           title = "Symptom Score Trajectory") +
      ylim(0, 18) + theme_classic()
    ggplotly(p)
  })

  output$healing_rate_plot <- renderPlotly({
    d <- sim_data() %>%
      mutate(Healed = ifelse(DMG < 10, "Healed", "Not Healed"),
             week   = floor(time / (24 * 7)))
    heal_by_week <- d %>%
      group_by(week) %>%
      summarise(heal_pct = 100 * mean(DMG < 10), .groups = "drop")
    p <- ggplot(heal_by_week, aes(x = week, y = heal_pct)) +
      geom_col(fill = "#43A047", alpha = 0.8) +
      labs(x = "Week", y = "Predicted Healing Rate (%)",
           title = "Mucosal Healing Prediction by Week") +
      ylim(0, 100) + theme_classic()
    ggplotly(p)
  })

  # ── Tab 5: Scenario Comparison ─────────────────────────────────────────────
  all_scenarios <- reactive({
    scenarios <- list(
      list(drug = "No Treatment", dose = 0, name = "No Treatment"),
      list(drug = "PPI", dose = 20, name = "Omeprazole 20 mg QD"),
      list(drug = "PPI", dose = 40, name = "Esomeprazole 40 mg QD"),
      list(drug = "P-CAB", dose = 20, name = "Vonoprazan 20 mg QD"),
      list(drug = "H2RA", dose = 40, name = "Famotidine 40 mg BID"),
      list(drug = "PPI+Prokinetic", dose = 40, dose2 = 10,
           name = "Eso 40 mg + Domperidone")
    )
    lapply(scenarios, function(s) {
      d <- run_sim(s$drug, s$dose, s$dose2 %||% 0,
                   cyp = as.numeric(input$cyp), duration = input$duration)
      d$Scenario <- s$name
      d
    }) %>% bind_rows()
  })

  output$compare_aet <- renderPlotly({
    d <- all_scenarios() %>% filter(time %% 6 == 0)
    p <- ggplot(d, aes(x = time / 24, y = AET_pct, colour = Scenario)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 6, linetype = "dashed", colour = "black") +
      labs(x = "Time (days)", y = "AET (%)",
           title = "Acid Exposure Time — All Scenarios") +
      theme_classic()
    ggplotly(p)
  })

  output$compare_table <- renderDT({
    d <- all_scenarios() %>%
      filter(time == max(time)) %>%
      group_by(Scenario) %>%
      summarise(
        "pH (Week-End)"     = round(mean(pH, na.rm = TRUE), 2),
        "AET % (Week-End)"  = round(mean(AET_pct, na.rm = TRUE), 1),
        "Acid Inhib %"      = round(mean(ACID_INH_pct, na.rm = TRUE), 1),
        "Mucosal Damage"    = round(mean(DMG, na.rm = TRUE), 1),
        "Symptom Score"     = round(mean(SYM, na.rm = TRUE), 1),
        "Barrett Risk"      = round(mean(Barrett, na.rm = TRUE), 4),
        .groups = "drop"
      )
    datatable(d, options = list(dom = "t", paging = FALSE), rownames = FALSE) %>%
      formatStyle("AET % (Week-End)",
                  backgroundColor = styleInterval(6, c("lightgreen", "lightyellow")))
  })

  # ── Tab 6: Biomarker & Risk ────────────────────────────────────────────────
  output$barrett_plot <- renderPlotly({
    d <- sim_data() %>% filter(time %% 12 == 0)
    p <- ggplot(d, aes(x = time / 24, y = Barrett * 100)) +
      geom_line(colour = "#8E24AA", linewidth = 1.2) +
      geom_area(fill = "#CE93D8", alpha = 0.3) +
      labs(x = "Time (days)", y = "Cumulative Barrett Risk (%)",
           title = "Barrett's Esophagus Risk Accumulation") +
      theme_classic()
    ggplotly(p)
  })

  output$dr_plot <- renderPlotly({
    doses_ppi  <- c(5, 10, 20, 40, 80)
    doses_pcab <- c(5, 10, 20, 40)
    dr_ppi <- lapply(doses_ppi, function(d) {
      out <- run_sim("PPI", d, duration = 8)
      data.frame(Dose = d, AET = tail(out$AET_pct, 1), Drug = "PPI")
    }) %>% bind_rows()
    dr_pcab <- lapply(doses_pcab, function(d) {
      out <- run_sim("P-CAB", d, duration = 8)
      data.frame(Dose = d, AET = tail(out$AET_pct, 1), Drug = "P-CAB")
    }) %>% bind_rows()
    dr <- bind_rows(dr_ppi, dr_pcab)
    p <- ggplot(dr, aes(x = Dose, y = AET, colour = Drug, group = Drug)) +
      geom_line(linewidth = 1) + geom_point(size = 3) +
      scale_x_log10() +
      scale_colour_manual(values = c("PPI" = "#8E24AA", "P-CAB" = "#00897B")) +
      geom_hline(yintercept = 6, linetype = "dashed") +
      labs(x = "Dose (mg, log)", y = "AET at Week 8 (%)",
           title = "Dose-Response: AET at Week 8") +
      theme_classic()
    ggplotly(p)
  })

  output$risk_table <- renderDT({
    d <- sim_data() %>% filter(time == max(time))
    aet_val <- mean(d$AET_pct)
    dmg_val <- mean(d$DMG)
    sym_val <- mean(d$SYM)
    be_val  <- mean(d$Barrett)

    risk_df <- data.frame(
      Biomarker   = c("Acid Exposure Time", "Mucosal Damage", "Symptom Score", "Barrett Risk"),
      Value       = c(round(aet_val, 1), round(dmg_val, 1), round(sym_val, 1), round(be_val * 100, 3)),
      Unit        = c("%", "score", "GERD-Q", "%"),
      Normal      = c("< 6%", "< 10", "< 8", "< 0.1%"),
      Status      = c(
        ifelse(aet_val < 6, "Normal", ifelse(aet_val < 15, "Elevated", "High")),
        ifelse(dmg_val < 10, "Normal", ifelse(dmg_val < 50, "Moderate", "Severe")),
        ifelse(sym_val < 8, "Normal", ifelse(sym_val < 12, "Symptomatic", "Severe")),
        ifelse(be_val < 0.001, "Low", ifelse(be_val < 0.01, "Intermediate", "High"))
      )
    )
    datatable(risk_df, options = list(dom = "t", paging = FALSE), rownames = FALSE) %>%
      formatStyle("Status",
                  backgroundColor = styleEqual(
                    c("Normal", "Low", "Elevated", "Moderate", "Symptomatic", "Intermediate", "High", "Severe"),
                    c("lightgreen", "lightgreen", "lightyellow", "lightyellow",
                      "lightyellow", "orange", "tomato", "tomato")
                  ))
  })
}

# ─── Helper ───────────────────────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

shinyApp(ui, server)
