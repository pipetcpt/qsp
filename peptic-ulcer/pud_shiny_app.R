## ============================================================
## Peptic Ulcer Disease (PUD) – Interactive QSP Shiny App
## ============================================================
## Tabs:
##   1. Patient Profile & Disease State
##   2. Drug Pharmacokinetics (PK)
##   3. PD Key Metrics (pH, Pump Inhibition, PG)
##   4. Clinical Endpoints (Pain, Ulcer, Bleeding)
##   5. Scenario Comparison
##   6. Biomarkers & Inflammation
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(plotly)
library(DT)
library(shinycssloaders)

# ── Compile mrgsolve model ────────────────────────────────────
pud_code <- '
$PARAM
F_PPI=0.65 ka_PPI=1.2 Vc_PPI=14.0 Vp_PPI=8.0 CL_PPI=28.0 Q_PPI=6.0
EC50_PPI=0.3 Emax_PPI=0.98
F_H2RA=0.45 ka_H2RA=0.8 Vc_H2RA=25.0 CL_H2RA=10.0
EC50_H2RA=0.05 Emax_H2RA=0.90
F_AMX=0.85 ka_AMX=1.5 Vc_AMX=22.0 CL_AMX=18.0 MIC_AMX=0.125 Emax_AMX=0.85
F_CLR=0.55 ka_CLR=1.0 Vc_CLR=72.0 CL_CLR=30.0 MIC_CLR=0.25 Emax_CLR=0.90
ResistanceCLR=0.15
F_NSAID=0.80 ka_NSAID=2.0 Vc_NSAID=8.0 CL_NSAID=4.0
IC50_COX1=0.5 IC50_COX2=2.0 Imax_COX=0.95
HP0=1e7 k_grow_HP=0.02 K_HP=1e9 k_clear_HP=0.001
pH_kill_HP=3.0 k_pH_kill=0.05
MaxAcidOutput=40.0 BasalAcid=4.0 k_pumpTurn=0.12
pHbasal=1.5 pHmax=7.0 HP_acid_stim=0.3
MucusMax=1.0 k_mucusProd=0.10 k_mucusDeg=0.05 PG_protect=0.8
k_PG_base=0.10 k_PG_deg=0.08 MucusIC50=4.0
k_inflam=0.05 k_inflam_res=0.02 InflamMax=10.0 IL8_HP=0.8
k_damage=0.15 k_heal=0.08 k_recur=0.005 UlcerMax=100.0 EGF_heal=0.3

$CMT PPI_GUT PPI_CENT PPI_PERI H2RA_GUT H2RA_CENT
AMX_GUT AMX_CENT CLR_GUT CLR_CENT NSAID_GUT NSAID_CENT
HP_LOAD PUMP_ACTIVE pH_INTRAGASTRIC MUCUS PG_LEVEL INFLAM ULCER_AREA

$MAIN
double C_PPI   = PPI_CENT / Vc_PPI;
double C_H2RA  = H2RA_CENT / Vc_H2RA;
double C_AMX   = AMX_CENT / Vc_AMX;
double C_CLR   = CLR_CENT / Vc_CLR;
double C_NSAID = NSAID_CENT / Vc_NSAID;
double inh_PPI  = Emax_PPI * C_PPI / (EC50_PPI + C_PPI);
double inh_H2RA = Emax_H2RA * C_H2RA / (EC50_H2RA + C_H2RA);
double total_acid_inh = 1.0 - (1.0 - inh_PPI) * (1.0 - inh_H2RA);
double inh_COX1 = Imax_COX * C_NSAID / (IC50_COX1 + C_NSAID);
double PG_NSAID_factor = 1.0 - 0.7 * inh_COX1;
double kill_AMX = (C_AMX > MIC_AMX) ? Emax_AMX * (C_AMX - MIC_AMX) / (C_AMX - MIC_AMX + MIC_AMX) : 0.0;
double kill_CLR = (C_CLR > MIC_CLR && ResistanceCLR < 0.5) ?
  Emax_CLR * (C_CLR - MIC_CLR) / (C_CLR - MIC_CLR + MIC_CLR) : 0.0;
double kill_combo = 1.0 - (1.0 - kill_AMX) * (1.0 - kill_CLR);
double pH = pH_INTRAGASTRIC;
double acid_HP_kill = (pH < pH_kill_HP) ? k_pH_kill * (pH_kill_HP - pH) : 0.0;

$ODE
dxdt_PPI_GUT  = -ka_PPI * PPI_GUT;
dxdt_PPI_CENT =  ka_PPI * F_PPI * PPI_GUT - (CL_PPI + Q_PPI) * (PPI_CENT/Vc_PPI) + Q_PPI * (PPI_PERI/Vp_PPI);
dxdt_PPI_PERI =  Q_PPI * (PPI_CENT/Vc_PPI) - Q_PPI * (PPI_PERI/Vp_PPI);
dxdt_H2RA_GUT  = -ka_H2RA * H2RA_GUT;
dxdt_H2RA_CENT =  ka_H2RA * F_H2RA * H2RA_GUT - CL_H2RA * C_H2RA;
dxdt_AMX_GUT  = -ka_AMX * AMX_GUT;
dxdt_AMX_CENT =  ka_AMX * F_AMX * AMX_GUT - CL_AMX * C_AMX;
dxdt_CLR_GUT  = -ka_CLR * CLR_GUT;
dxdt_CLR_CENT =  ka_CLR * F_CLR * CLR_GUT - CL_CLR * C_CLR;
dxdt_NSAID_GUT  = -ka_NSAID * NSAID_GUT;
dxdt_NSAID_CENT =  ka_NSAID * F_NSAID * NSAID_GUT - CL_NSAID * C_NSAID;
double HP = HP_LOAD;
double HP_growth = k_grow_HP * HP * (1.0 - HP / K_HP);
double HP_kill_total = (kill_combo + acid_HP_kill + k_clear_HP) * HP;
dxdt_HP_LOAD = HP_growth - HP_kill_total;
double HP_stim = HP_acid_stim * (HP / (HP + 1e6));
double pump_activation = k_pumpTurn * (1.0 + HP_stim) * (1.0 - PUMP_ACTIVE);
double pump_inhibition_rate = inh_PPI * k_pumpTurn;
dxdt_PUMP_ACTIVE = pump_activation - pump_inhibition_rate * PUMP_ACTIVE;
double acid_rate = (MaxAcidOutput - BasalAcid) * PUMP_ACTIVE * (1.0 - total_acid_inh) + BasalAcid;
double target_pH = pHbasal + (pHmax - pHbasal) * (1.0 - acid_rate / MaxAcidOutput);
dxdt_pH_INTRAGASTRIC = 1.5 * (target_pH - pH_INTRAGASTRIC);
double mucus_prod = k_mucusProd * PG_LEVEL * PG_NSAID_factor;
double mucus_acid_deg = k_mucusDeg * (1.0 + (1.0 - PG_NSAID_factor)) * exp(-0.5 * (pH - 2.0));
double mucus_HP_deg = 0.03 * (HP / K_HP);
dxdt_MUCUS = mucus_prod * (MucusMax - MUCUS) - (mucus_acid_deg + mucus_HP_deg) * MUCUS;
double PG_prod = k_PG_base * PG_NSAID_factor + 0.02 * INFLAM;
double PG_deg  = k_PG_deg * PG_LEVEL;
dxdt_PG_LEVEL = PG_prod - PG_deg;
double inflam_drive = k_inflam * (HP / K_HP) * InflamMax;
double NSAID_inflam = 0.1 * inh_COX1 * InflamMax;
double inflam_res = k_inflam_res * INFLAM * (MUCUS + 0.1);
dxdt_INFLAM = inflam_drive + NSAID_inflam - inflam_res;
double mucosal_integrity = MUCUS * PG_LEVEL;
double acid_damage_factor = (pH < 3.0) ? (3.0 - pH) : 0.0;
double ulcer_damage = k_damage * acid_damage_factor * INFLAM / (mucosal_integrity + 0.1);
double EGF_boost = 1.0 + EGF_heal * (1.0 - ULCER_AREA / UlcerMax);
double ulcer_healing = k_heal * EGF_boost * PG_LEVEL * MUCUS * ULCER_AREA;
double ulcer_recur = k_recur * (HP / K_HP) * (UlcerMax - ULCER_AREA);
dxdt_ULCER_AREA = ulcer_damage * (UlcerMax - ULCER_AREA) + ulcer_recur - ulcer_healing;

$TABLE
double C_PPI_out   = PPI_CENT / Vc_PPI;
double C_H2RA_out  = H2RA_CENT / Vc_H2RA;
double C_AMX_out   = AMX_CENT / Vc_AMX;
double C_CLR_out   = CLR_CENT / Vc_CLR;
double C_NSAID_out = NSAID_CENT / Vc_NSAID;
double pump_inh_pct = (1.0 - PUMP_ACTIVE) * 100.0;
double log10_HP = (HP_LOAD > 1) ? log10(HP_LOAD) : 0.0;
double HP_erad_flag = (HP_LOAD < 100) ? 1.0 : 0.0;
double mucosal_prot = MUCUS * PG_LEVEL;
double pain_VAS = 10.0 * (ULCER_AREA / UlcerMax) * (0.5 + 0.5 * INFLAM / 10.0);
double bleeding_risk = (ULCER_AREA > 20.0) ? 0.01 * (ULCER_AREA - 20.0) / (UlcerMax - 20.0) : 0.0;
double endoscopy_score = (ULCER_AREA < 1) ? 0.0 : (ULCER_AREA < 10) ? 1.0 :
  (ULCER_AREA < 25) ? 2.0 : (ULCER_AREA < 50) ? 3.0 : 4.0;
double COX1_inh_pct = Imax_COX * C_NSAID_out / (IC50_COX1 + C_NSAID_out) * 100.0;
capture C_PPI_out C_H2RA_out C_AMX_out C_CLR_out C_NSAID_out
capture pH_INTRAGASTRIC pump_inh_pct log10_HP HP_erad_flag
capture mucosal_prot pain_VAS bleeding_risk endoscopy_score COX1_inh_pct
capture MUCUS PG_LEVEL INFLAM ULCER_AREA PUMP_ACTIVE

$INIT
PPI_GUT=0 PPI_CENT=0 PPI_PERI=0 H2RA_GUT=0 H2RA_CENT=0
AMX_GUT=0 AMX_CENT=0 CLR_GUT=0 CLR_CENT=0 NSAID_GUT=0 NSAID_CENT=0
HP_LOAD=1e7 PUMP_ACTIVE=0.95 pH_INTRAGASTRIC=1.5
MUCUS=0.6 PG_LEVEL=0.7 INFLAM=4.0 ULCER_AREA=25.0
'
mod <- mcode("PUD_shiny", pud_code, quiet = TRUE)

# ── Helper: run simulation ────────────────────────────────────
run_sim <- function(
    ppi_dose = 0, ppi_freq = 12,
    h2ra_dose = 0, h2ra_freq = 12,
    amx_dose  = 0, clr_dose  = 0,
    nsaid_dose = 0, nsaid_freq = 8,
    dur_days = 28,
    hp_init = 1e7,
    ulcer_init = 25,
    inflam_init = 4.0,
    resistance_clr = 0.15
) {
  e_list <- list()
  if (ppi_dose  > 0) e_list[["ppi"]]  <- ev(amt = ppi_dose,  cmt = "PPI_GUT",  ii = ppi_freq,  addl = floor(dur_days * 24 / ppi_freq) - 1, time = 0)
  if (h2ra_dose > 0) e_list[["h2ra"]] <- ev(amt = h2ra_dose, cmt = "H2RA_GUT", ii = h2ra_freq, addl = floor(dur_days * 24 / h2ra_freq) - 1, time = 0)
  if (amx_dose  > 0) e_list[["amx"]]  <- ev(amt = amx_dose,  cmt = "AMX_GUT",  ii = 12,         addl = floor(dur_days * 24 / 12) - 1, time = 0)
  if (clr_dose  > 0) e_list[["clr"]]  <- ev(amt = clr_dose,  cmt = "CLR_GUT",  ii = 12,         addl = floor(dur_days * 24 / 12) - 1, time = 0)
  if (nsaid_dose > 0) e_list[["nsaid"]] <- ev(amt = nsaid_dose, cmt = "NSAID_GUT", ii = nsaid_freq, addl = floor(dur_days * 24 / nsaid_freq) - 1, time = 0)

  events <- if (length(e_list) > 0) do.call(c, e_list) else ev(amt = 0, cmt = 1, time = 9999)

  mod %>%
    param(ResistanceCLR = resistance_clr) %>%
    init(HP_LOAD = hp_init, ULCER_AREA = ulcer_init, INFLAM = inflam_init) %>%
    ev(events) %>%
    mrgsim(end = dur_days * 24, delta = 1) %>%
    as.data.frame() %>%
    mutate(time_day = time / 24)
}

# ── UI Definition ─────────────────────────────────────────────
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "PUD QSP Model"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",   tabName = "profile",   icon = icon("user")),
      menuItem("Pharmacokinetics",  tabName = "pk",        icon = icon("chart-line")),
      menuItem("PD Key Metrics",    tabName = "pd",        icon = icon("flask")),
      menuItem("Clinical Endpoints",tabName = "clinical",  icon = icon("stethoscope")),
      menuItem("Scenario Comparison",tabName = "scenarios", icon = icon("balance-scale")),
      menuItem("Biomarkers",        tabName = "biomarkers", icon = icon("dna"))
    ),

    hr(),
    h4("  Treatment Settings", style = "color:white; padding-left:10px"),

    # PPI settings
    checkboxInput("use_ppi", "PPI (Omeprazole)", value = TRUE),
    conditionalPanel("input.use_ppi",
      sliderInput("ppi_dose", "PPI Dose (mg)", 10, 40, 20, step = 10),
      selectInput("ppi_freq", "PPI Frequency (hr)", choices = c("12 (BID)" = 12, "24 (QD)" = 24), selected = 12)
    ),

    # H2RA settings
    checkboxInput("use_h2ra", "H2RA (Famotidine)", value = FALSE),
    conditionalPanel("input.use_h2ra",
      sliderInput("h2ra_dose", "H2RA Dose (mg)", 10, 40, 20, step = 10)
    ),

    # H. pylori eradication
    checkboxInput("use_triple", "H. pylori Triple Therapy", value = FALSE),
    conditionalPanel("input.use_triple",
      sliderInput("amx_dose", "Amoxicillin (mg)", 500, 1500, 1000, step = 250),
      sliderInput("clr_dose", "Clarithromycin (mg)", 250, 500, 500, step = 250),
      sliderInput("clr_resist", "CLR Resistance (%)", 0, 50, 15, step = 5)
    ),

    # NSAID settings
    checkboxInput("use_nsaid", "NSAID (Ibuprofen)", value = FALSE),
    conditionalPanel("input.use_nsaid",
      sliderInput("nsaid_dose", "NSAID Dose (mg)", 200, 800, 400, step = 200),
      selectInput("nsaid_freq", "NSAID Frequency (hr)", choices = c("8 (TID)" = 8, "12 (BID)" = 12), selected = 8)
    ),

    hr(),
    h4("  Patient Characteristics", style = "color:white; padding-left:10px"),
    sliderInput("hp_severity", "H. pylori Load (log10 CFU)", 3, 9, 7, step = 1),
    sliderInput("ulcer_init", "Initial Ulcer Area (mm²)", 0, 80, 25, step = 5),
    sliderInput("inflam_init", "Baseline Inflammation (0-10)", 0, 10, 4, step = 1),
    sliderInput("dur_days", "Simulation Duration (days)", 7, 56, 28, step = 7),
    actionButton("run_btn", "Run Simulation", class = "btn-primary btn-block")
  ),

  dashboardBody(
    tabItems(

      # ── TAB 1: Patient Profile ──────────────────────────────
      tabItem(tabName = "profile",
        fluidRow(
          box(title = "Disease Overview: Peptic Ulcer Disease", width = 12, status = "danger",
            p("Peptic ulcer disease (PUD) results from an imbalance between aggressive factors
               (gastric acid, H. pylori infection, NSAIDs) and mucosal defense mechanisms
               (mucus layer, bicarbonate, prostaglandins)."),
            p("This QSP model integrates:"),
            tags$ul(
              tags$li("H. pylori colonization dynamics and eradication pharmacodynamics"),
              tags$li("Gastric acid secretion via proton pump and H2 receptor pathways"),
              tags$li("Mucosal defense (mucus, PGE2) and NSAID-mediated disruption"),
              tags$li("Ulcer formation, progression, and healing dynamics"),
              tags$li("Multiple treatment scenarios (PPI, H2RA, triple therapy, misoprostol)")
            )
          )
        ),
        fluidRow(
          valueBoxOutput("box_ulcer_area", width = 3),
          valueBoxOutput("box_hp_load", width = 3),
          valueBoxOutput("box_ph", width = 3),
          valueBoxOutput("box_pain", width = 3)
        ),
        fluidRow(
          box(title = "Patient State at End of Simulation", width = 6, status = "primary",
            withSpinner(DTOutput("table_summary"))
          ),
          box(title = "Disease Mechanism Schematic", width = 6, status = "info",
            tags$img(
              src = "https://upload.wikimedia.org/wikipedia/commons/thumb/c/cf/Peptic_ulcer_diagram.svg/640px-Peptic_ulcer_diagram.svg.png",
              width = "100%", style = "border-radius:8px"
            ),
            p("Peptic ulcer: breach of mucosal barrier by acid/pepsin", style = "font-size:11px; color:gray")
          )
        )
      ),

      # ── TAB 2: Pharmacokinetics ─────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "PPI Plasma Concentration-Time Profile", width = 6, status = "primary",
            withSpinner(plotlyOutput("plot_ppi_pk", height = 300))
          ),
          box(title = "H2RA / NSAID Plasma Concentration", width = 6, status = "primary",
            withSpinner(plotlyOutput("plot_other_pk", height = 300))
          )
        ),
        fluidRow(
          box(title = "Antibiotic Plasma Concentrations (AMX / CLR)", width = 6, status = "warning",
            withSpinner(plotlyOutput("plot_atb_pk", height = 300))
          ),
          box(title = "PK Summary Statistics", width = 6, status = "info",
            withSpinner(DTOutput("table_pk_stats"))
          )
        )
      ),

      # ── TAB 3: PD Key Metrics ───────────────────────────────
      tabItem(tabName = "pd",
        fluidRow(
          box(title = "Intragastric pH Over Time", width = 6, status = "primary",
            withSpinner(plotlyOutput("plot_ph", height = 300))
          ),
          box(title = "Proton Pump Inhibition (%)", width = 6, status = "primary",
            withSpinner(plotlyOutput("plot_pump_inh", height = 300))
          )
        ),
        fluidRow(
          box(title = "Mucus Layer Integrity (0-1)", width = 6, status = "warning",
            withSpinner(plotlyOutput("plot_mucus", height = 300))
          ),
          box(title = "Prostaglandin Level (0-1)", width = 6, status = "warning",
            withSpinner(plotlyOutput("plot_pg", height = 300))
          )
        )
      ),

      # ── TAB 4: Clinical Endpoints ───────────────────────────
      tabItem(tabName = "clinical",
        fluidRow(
          box(title = "Ulcer Area Dynamics (mm²)", width = 6, status = "danger",
            withSpinner(plotlyOutput("plot_ulcer", height = 300))
          ),
          box(title = "Epigastric Pain VAS (0-10)", width = 6, status = "danger",
            withSpinner(plotlyOutput("plot_pain", height = 300))
          )
        ),
        fluidRow(
          box(title = "H. pylori Load (log10 CFU/mL)", width = 6, status = "warning",
            withSpinner(plotlyOutput("plot_hp", height = 300))
          ),
          box(title = "GI Bleeding Risk (%)", width = 6, status = "danger",
            withSpinner(plotlyOutput("plot_bleed", height = 300))
          )
        )
      ),

      # ── TAB 5: Scenario Comparison ──────────────────────────
      tabItem(tabName = "scenarios",
        fluidRow(
          box(title = "Pre-defined Treatment Scenarios", width = 12, status = "primary",
            p("Compare standard-of-care regimens: No treatment, PPI, H2RA, H. pylori triple therapy, NSAID ± PPI"),
            withSpinner(plotlyOutput("plot_scenarios", height = 500))
          )
        ),
        fluidRow(
          box(title = "Scenario Endpoint Summary Table", width = 12, status = "info",
            withSpinner(DTOutput("table_scenarios"))
          )
        )
      ),

      # ── TAB 6: Biomarkers ───────────────────────────────────
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(title = "Inflammation Score Over Time", width = 6, status = "info",
            withSpinner(plotlyOutput("plot_inflam", height = 300))
          ),
          box(title = "COX-1 Inhibition % (NSAID Effect)", width = 6, status = "warning",
            withSpinner(plotlyOutput("plot_cox1", height = 300))
          )
        ),
        fluidRow(
          box(title = "Mucosal Protection Score (Mucus × PG)", width = 6, status = "success",
            withSpinner(plotlyOutput("plot_muc_prot", height = 300))
          ),
          box(title = "Endoscopy Score (Lanza 0-4)", width = 6, status = "info",
            withSpinner(plotlyOutput("plot_endo", height = 300))
          )
        )
      )
    )
  )
)

# ── Server Logic ─────────────────────────────────────────────
server <- function(input, output, session) {

  sim_data <- eventReactive(input$run_btn, {
    run_sim(
      ppi_dose   = if (input$use_ppi)   input$ppi_dose  else 0,
      ppi_freq   = as.numeric(input$ppi_freq),
      h2ra_dose  = if (input$use_h2ra)  input$h2ra_dose else 0,
      h2ra_freq  = 12,
      amx_dose   = if (input$use_triple) input$amx_dose  else 0,
      clr_dose   = if (input$use_triple) input$clr_dose  else 0,
      nsaid_dose = if (input$use_nsaid) input$nsaid_dose else 0,
      nsaid_freq = as.numeric(input$nsaid_freq),
      dur_days   = input$dur_days,
      hp_init    = 10^input$hp_severity,
      ulcer_init = input$ulcer_init,
      inflam_init = input$inflam_init,
      resistance_clr = if (input$use_triple) input$clr_resist / 100 else 0.15
    )
  }, ignoreNULL = FALSE)

  # ── Value boxes ───────────────────────────────────────────
  get_last <- function(col) {
    d <- sim_data()
    if (is.null(d) || nrow(d) == 0) return(NA)
    tail(d[[col]], 1)
  }

  output$box_ulcer_area <- renderValueBox({
    v <- round(get_last("ULCER_AREA"), 1)
    valueBox(paste0(v, " mm²"), "Ulcer Area", icon = icon("circle"), color = "red")
  })
  output$box_hp_load <- renderValueBox({
    v <- round(get_last("log10_HP"), 2)
    valueBox(paste0("10^", v), "H. pylori Load", icon = icon("bacteria"), color = "orange")
  })
  output$box_ph <- renderValueBox({
    v <- round(get_last("pH_INTRAGASTRIC"), 2)
    valueBox(v, "Intragastric pH", icon = icon("flask"), color = "blue")
  })
  output$box_pain <- renderValueBox({
    v <- round(get_last("pain_VAS"), 1)
    valueBox(v, "Pain VAS (0-10)", icon = icon("heartbeat"), color = "yellow")
  })

  # ── Summary table ─────────────────────────────────────────
  output$table_summary <- renderDT({
    d <- sim_data()
    if (is.null(d)) return(NULL)
    last <- tail(d, 1)
    df <- data.frame(
      Parameter = c("Ulcer Area (mm²)", "H. pylori (log10 CFU)", "Intragastric pH",
                    "Pain VAS", "Bleeding Risk (%)", "Endoscopy Score", "Mucosal Protection",
                    "Inflammation Score", "HP Eradicated"),
      Value = round(c(last$ULCER_AREA, last$log10_HP, last$pH_INTRAGASTRIC,
                      last$pain_VAS, last$bleeding_risk * 100, last$endoscopy_score,
                      last$mucosal_prot, last$INFLAM,
                      last$HP_erad_flag), 2)
    )
    datatable(df, options = list(pageLength = 10, dom = "t"), rownames = FALSE)
  })

  # ── PK plots ──────────────────────────────────────────────
  output$plot_ppi_pk <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_day, C_PPI_out)) +
      geom_line(color = "#E74C3C", linewidth = 1) +
      labs(x = "Time (days)", y = "Plasma Concentration (mg/L)", title = "PPI (Omeprazole)") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_other_pk <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d) +
      geom_line(aes(time_day, C_H2RA_out, color = "H2RA"), linewidth = 1) +
      geom_line(aes(time_day, C_NSAID_out, color = "NSAID"), linewidth = 1) +
      scale_color_manual(values = c("H2RA" = "#3498DB", "NSAID" = "#E67E22")) +
      labs(x = "Time (days)", y = "Concentration (mg/L)", color = "Drug") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_atb_pk <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d) +
      geom_line(aes(time_day, C_AMX_out, color = "Amoxicillin"), linewidth = 1) +
      geom_line(aes(time_day, C_CLR_out, color = "Clarithromycin"), linewidth = 1) +
      geom_hline(yintercept = 0.125, linetype = "dashed", color = "purple", alpha = 0.7) +
      geom_hline(yintercept = 0.25,  linetype = "dashed", color = "blue", alpha = 0.7) +
      scale_color_manual(values = c("Amoxicillin" = "#9B59B6", "Clarithromycin" = "#1ABC9C")) +
      labs(x = "Time (days)", y = "Concentration (mg/L)", color = "Drug",
           subtitle = "Dashed lines = MIC") +
      theme_bw()
    ggplotly(p)
  })

  output$table_pk_stats <- renderDT({
    d <- sim_data()
    if (is.null(d)) return(NULL)
    df <- data.frame(
      Drug = c("PPI", "H2RA", "Amoxicillin", "Clarithromycin", "NSAID"),
      Cmax = round(c(max(d$C_PPI_out), max(d$C_H2RA_out), max(d$C_AMX_out),
                     max(d$C_CLR_out), max(d$C_NSAID_out)), 3),
      AUC_24h = round(c(
        sum(d$C_PPI_out[d$time <= 24]),
        sum(d$C_H2RA_out[d$time <= 24]),
        sum(d$C_AMX_out[d$time <= 24]),
        sum(d$C_CLR_out[d$time <= 24]),
        sum(d$C_NSAID_out[d$time <= 24])
      ), 2)
    )
    datatable(df, options = list(dom = "t"), rownames = FALSE)
  })

  # ── PD plots ──────────────────────────────────────────────
  make_plotly <- function(d, y, ylab, color = "#2980B9", hline = NULL) {
    p <- ggplot(d, aes(time_day, .data[[y]])) +
      geom_line(color = color, linewidth = 1) +
      labs(x = "Time (days)", y = ylab) +
      theme_bw()
    if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = "red")
    ggplotly(p)
  }

  output$plot_ph      <- renderPlotly(make_plotly(sim_data(), "pH_INTRAGASTRIC", "pH", "#2ECC71", hline = 4.0))
  output$plot_pump_inh<- renderPlotly(make_plotly(sim_data(), "pump_inh_pct",    "Pump Inhibition (%)", "#E74C3C"))
  output$plot_mucus   <- renderPlotly(make_plotly(sim_data(), "MUCUS",            "Mucus Layer (0-1)", "#27AE60"))
  output$plot_pg      <- renderPlotly(make_plotly(sim_data(), "PG_LEVEL",         "Prostaglandin Level (0-1)", "#F39C12"))

  # ── Clinical endpoint plots ────────────────────────────────
  output$plot_ulcer   <- renderPlotly(make_plotly(sim_data(), "ULCER_AREA", "Ulcer Area (mm²)", "#C0392B"))
  output$plot_pain    <- renderPlotly(make_plotly(sim_data(), "pain_VAS",   "Pain VAS (0-10)",  "#E74C3C"))
  output$plot_hp      <- renderPlotly(make_plotly(sim_data(), "log10_HP",   "log10(HP CFU/mL)", "#8E44AD", hline = log10(100)))
  output$plot_bleed   <- renderPlotly({
    d <- sim_data() %>% mutate(bleed_pct = bleeding_risk * 100)
    make_plotly(d, "bleed_pct", "Bleeding Risk (%)", "#922B21")
  })

  # ── Scenario comparison ────────────────────────────────────
  scenarios_data <- reactive({
    base_args <- list(hp_init = 10^input$hp_severity, ulcer_init = input$ulcer_init,
                      inflam_init = input$inflam_init, dur_days = input$dur_days)

    sc <- list(
      "No Treatment"   = do.call(run_sim, base_args),
      "PPI 20mg BID"   = do.call(run_sim, c(base_args, list(ppi_dose = 20))),
      "H2RA 20mg BID"  = do.call(run_sim, c(base_args, list(h2ra_dose = 20))),
      "Triple Therapy" = do.call(run_sim, c(base_args, list(ppi_dose = 20, amx_dose = 1000, clr_dose = 500))),
      "NSAID Only"     = do.call(run_sim, c(base_args, list(hp_init = 100, ulcer_init = 0, inflam_init = 1, nsaid_dose = 400))),
      "NSAID + PPI"    = do.call(run_sim, c(base_args, list(hp_init = 100, ulcer_init = 0, inflam_init = 1, nsaid_dose = 400, ppi_dose = 20)))
    )
    bind_rows(lapply(names(sc), function(nm) sc[[nm]] %>% mutate(scenario = nm)))
  })

  output$plot_scenarios <- renderPlotly({
    d <- scenarios_data()
    p <- ggplot(d, aes(time_day, ULCER_AREA, color = scenario)) +
      geom_line(linewidth = 1) +
      labs(title = "Ulcer Area by Treatment Scenario",
           x = "Time (days)", y = "Ulcer Area (mm²)", color = "Scenario") +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$table_scenarios <- renderDT({
    d <- scenarios_data()
    d %>%
      group_by(scenario) %>%
      filter(time == max(time)) %>%
      summarise(
        pH         = round(mean(pH_INTRAGASTRIC), 2),
        log10_HP   = round(mean(log10_HP), 2),
        HP_Erad    = round(mean(HP_erad_flag) * 100, 1),
        Ulcer_mm2  = round(mean(ULCER_AREA), 1),
        Pain_VAS   = round(mean(pain_VAS), 2),
        Bleed_risk = round(mean(bleeding_risk) * 100, 2)
      ) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # ── Biomarker plots ────────────────────────────────────────
  output$plot_inflam   <- renderPlotly(make_plotly(sim_data(), "INFLAM", "Inflammation (0-10)", "#E67E22"))
  output$plot_cox1     <- renderPlotly(make_plotly(sim_data(), "COX1_inh_pct", "COX-1 Inhibition (%)", "#E74C3C"))
  output$plot_muc_prot <- renderPlotly(make_plotly(sim_data(), "mucosal_prot", "Mucosal Protection (0-1)", "#27AE60"))
  output$plot_endo     <- renderPlotly(make_plotly(sim_data(), "endoscopy_score", "Lanza Endoscopy Score", "#8E44AD"))
}

shinyApp(ui, server)
