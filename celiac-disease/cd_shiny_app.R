# ==============================================================================
# Celiac Disease QSP Shiny App
# Interactive Quantitative Systems Pharmacology Dashboard
# ==============================================================================
# 6 Tabs:
#   1. Patient Profile & Disease Summary
#   2. Drug PK (Pharmacokinetics)
#   3. Immune Biomarkers (IL-15, IFN-γ, IL-17, IL-21, IEL)
#   4. Intestinal Histopathology (Marsh Score, V:C Ratio, Villous/Crypt)
#   5. Clinical Endpoints (Anti-tTG IgA, Hgb, BMD, Absorption)
#   6. Treatment Scenario Comparison
# ==============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(shinycssloaders)

# ==============================================================================
# Embedded mrgsolve Model
# ==============================================================================
cd_model_code <- '
$PARAM
GlutenIntake_g_day : 12
GFD               : 0
GFD_leak          : 0
kabs_gliadin   : 0.25
kdeg_lumen     : 0.40
k_deamid       : 0.30
kdeg_DGP       : 0.35
k_tTG2_TGFb    : 0.05
k_tTG2_deg     : 0.04
IP_basal       : 1.0
k_IP_IL15      : 0.08
k_IP_restore   : 0.03
EC50_IP_IL15   : 2.5
k_IL15_prod    : 0.15
k_IL15_deg     : 0.12
EC50_IL15_glia : 0.8
IEL_normal     : 20
k_IEL_IL15     : 0.06
k_IEL_death    : 0.025
IEL_max        : 120
EC50_IEL_IL15  : 3.0
CD4T_basal     : 100
k_CD4_prol     : 0.07
k_CD4_death    : 0.03
EC50_CD4_DGP   : 0.5
k_IFNg_prod    : 0.18
k_IFNg_deg     : 0.35
k_IL17_prod    : 0.09
k_IL17_deg     : 0.30
k_IL21_prod    : 0.10
k_IL21_deg     : 0.25
k_Bprol        : 0.05
k_Bdeath       : 0.02
k_AntiTTG_prod : 0.06
k_AntiTTG_deg  : 0.005
AntiTTG_ULN   : 10
VH0            : 1.0
VH_min         : 0.05
k_VH_damage    : 0.004
k_VH_repair    : 0.002
EC50_VH_IFNg   : 4.0
CD0            : 1.0
k_CD_hyp       : 0.003
k_CD_norm      : 0.001
k_Iron_abs_hr  : 0.083
k_Iron_loss_hr : 0.042
k_BMD_loss     : 0.0001
k_BMD_repair   : 0.00005
F_oral         : 0.5
ka_drug        : 0.6
CL_drug_L_h   : 8.0
Vd_drug_L     : 6.0
Drug_type      : 0
EC50_lara_IP   : 50
EC50_ZED_tTG   : 80
EC50_AMG_IL15  : 15
Emax_drug      : 0.85

$INIT
GlutenLumen  : 0
GlutenLP     : 0
DGP          : 0
tTG2_act     : 1.0
IP           : 1.0
IL15         : 0.1
IEL          : 20
CD4T         : 100
IFNg         : 0.1
IL17         : 0.05
IL21         : 0.05
Bcell        : 50
AntiTTG      : 2.0
VH           : 1.0
CrD          : 1.0
AbsArea      : 1.0
IronStores   : 1000
BMD          : 1.0
DrugGut      : 0
DrugPlasma   : 0

$ODE
double gluten_rate_h = GlutenIntake_g_day * ((1.0-GFD) + GFD*GFD_leak) / 24.0;
double Cp = DrugPlasma;
double E_lara = (Drug_type==1) ? Emax_drug*Cp/(Cp+EC50_lara_IP)  : 0.0;
double E_ZED  = (Drug_type==2) ? Emax_drug*Cp/(Cp+EC50_ZED_tTG)  : 0.0;
double E_AMG  = (Drug_type==3) ? Emax_drug*Cp/(Cp+EC50_AMG_IL15) : 0.0;
double h_IL15_IP  = IL15/(IL15+EC50_IP_IL15);
double h_IL15_IEL = IL15/(IL15+EC50_IEL_IL15);
double h_glia_IL15 = GlutenLP/(GlutenLP+EC50_IL15_glia);
double h_DGP_CD4   = DGP/(DGP+EC50_CD4_DGP);
double h_IFNg_VH   = IFNg/(IFNg+EC50_VH_IFNg);
dxdt_GlutenLumen = gluten_rate_h - kdeg_lumen*GlutenLumen - kabs_gliadin*GlutenLumen*IP;
dxdt_GlutenLP = kabs_gliadin*GlutenLumen*IP - k_deamid*tTG2_act*GlutenLP*(1.0-E_ZED) - kdeg_DGP*GlutenLP;
dxdt_DGP = k_deamid*tTG2_act*GlutenLP*(1.0-E_ZED) - kdeg_DGP*DGP;
dxdt_tTG2_act = k_tTG2_TGFb*(1.0+0.3*IFNg/(IFNg+2.0)) - k_tTG2_deg*tTG2_act;
dxdt_IP = k_IP_IL15*h_IL15_IP*(1.0-E_lara) - k_IP_restore*(IP-IP_basal);
dxdt_IL15 = k_IL15_prod*h_glia_IL15*(1.0-E_AMG) - k_IL15_deg*IL15;
double IEL_prol = k_IEL_IL15*h_IL15_IEL*IEL*(1.0-IEL/IEL_max);
dxdt_IEL = IEL_prol - k_IEL_death*(IEL-IEL_normal);
dxdt_CD4T = k_CD4_prol*h_DGP_CD4*CD4T - k_CD4_death*(CD4T-CD4T_basal);
dxdt_IFNg = k_IFNg_prod*(CD4T/100.0)*0.4 - k_IFNg_deg*IFNg;
dxdt_IL17 = k_IL17_prod*(CD4T/100.0)*0.2*h_IL15_IP - k_IL17_deg*IL17;
dxdt_IL21 = k_IL21_prod*(CD4T/100.0)*0.3 - k_IL21_deg*IL21;
double Bcell_prol = k_Bprol*h_DGP_CD4*(IL21/(IL21+0.1))*Bcell;
dxdt_Bcell = Bcell_prol - k_Bdeath*Bcell;
dxdt_AntiTTG = k_AntiTTG_prod*(Bcell/50.0) - k_AntiTTG_deg*AntiTTG;
double VH_damage = k_VH_damage*h_IFNg_VH*(IL17/(IL17+1.0))*VH;
double VH_repair = k_VH_repair*(1.0-h_DGP_CD4)*(VH0-VH);
dxdt_VH = VH_repair - VH_damage;
dxdt_CrD = k_CD_hyp*h_IL15_IP*(2.5-CrD) - k_CD_norm*(CrD-CD0)*(GlutenLP<0.1?1.0:0.0);
double Abs_target = fmax(0.05, VH/(VH+0.2));
dxdt_AbsArea = 0.015*(Abs_target-AbsArea);
dxdt_IronStores = k_Iron_abs_hr*AbsArea - k_Iron_loss_hr;
dxdt_BMD = k_BMD_repair*AbsArea - k_BMD_loss*(1.0-AbsArea);
dxdt_DrugGut = -ka_drug*DrugGut;
dxdt_DrugPlasma = F_oral*ka_drug*DrugGut*1000.0/Vd_drug_L - (CL_drug_L_h/Vd_drug_L)*DrugPlasma;

$TABLE
capture VH_CD_ratio = VH/CrD;
capture Marsh_score = (VH_CD_ratio<0.3) ? 3 : (VH_CD_ratio<0.7) ? 2 : (VH_CD_ratio<1.0) ? 1 : 0;
capture Serology_pos = (AntiTTG > AntiTTG_ULN) ? 1 : 0;
capture Hgb_g_dL = 8.0 + 6.0*(IronStores/1000.0);
capture Ferritin_ug = 15.0*(IronStores/1000.0);
capture BMD_Tscore = (BMD-1.0)/0.1*(-1);
capture IEL_elevated = (IEL>25) ? 1 : 0;
$CAPTURE VH_CD_ratio Marsh_score Serology_pos Hgb_g_dL Ferritin_ug
         BMD_Tscore IEL_elevated AntiTTG VH CrD AbsArea IronStores BMD
         IL15 IFNg IL17 IL21 IEL CD4T Bcell DrugPlasma GlutenLumen DGP IP
'

# Build model once at startup
cd_mod <- mcode("CD_Shiny", cd_model_code)

# ==============================================================================
# Helper: Run Simulation
# ==============================================================================
run_sim <- function(params_list, drug_type, drug_dose_mg, dose_interval_h, sim_days) {
  sim_end <- sim_days * 24

  if (drug_dose_mg > 0 && drug_type > 0) {
    ev_times <- seq(0, sim_end - dose_interval_h, by = dose_interval_h)
    dose_ev <- ev(time = ev_times, amt = drug_dose_mg, cmt = "DrugGut")
  } else {
    dose_ev <- ev(time = 0, amt = 0, cmt = "DrugGut")
  }

  mod_updated <- do.call(param, c(list(cd_mod), params_list))

  mrgsim(
    x       = mod_updated,
    events  = dose_ev,
    end     = sim_end,
    delta   = 24,
    obsonly = TRUE
  ) %>%
    as.data.frame() %>%
    mutate(Time_days = time / 24)
}

# ==============================================================================
# UI
# ==============================================================================
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "Celiac Disease QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("① Patient Profile",     tabName = "profile",    icon = icon("user-md")),
      menuItem("② Drug PK",             tabName = "pk",         icon = icon("pills")),
      menuItem("③ Immune Biomarkers",   tabName = "immune",     icon = icon("microscope")),
      menuItem("④ Histopathology",      tabName = "histo",      icon = icon("chart-bar")),
      menuItem("⑤ Clinical Endpoints",  tabName = "clinical",   icon = icon("heartbeat")),
      menuItem("⑥ Scenario Comparison", tabName = "scenarios",  icon = icon("layer-group"))
    ),
    hr(),
    h4("  Global Controls", style = "color:white; padding-left:10px"),

    selectInput("disease_hla", "HLA Type",
      choices = c("HLA-DQ2.5 (homozygous, high risk)" = "dq2_homo",
                  "HLA-DQ2.5 (heterozygous)" = "dq2_het",
                  "HLA-DQ8" = "dq8",
                  "DQ2.2 (low risk)" = "dq22"),
      selected = "dq2_homo"),

    sliderInput("gluten_g", "Daily Gluten Intake (g/day)", 0, 20, 12, step = 1),

    radioButtons("gfd", "GFD Compliance",
      choices = c("None" = 0, "Partial (5% leak)" = "partial", "Strict" = 1),
      selected = 0),

    selectInput("drug_type", "Add Drug",
      choices = c("None" = 0,
                  "Larazotide (TID 2 mg)" = 1,
                  "ZED1227 TG2i (QD 300 mg)" = 2,
                  "AMG714 anti-IL-15 (150 mg SC wkly)" = 3),
      selected = 0),

    sliderInput("sim_days", "Simulation Duration (days)", 30, 730, 365, step = 30),

    actionButton("run_sim", "▶ Run Simulation",
                 class = "btn-success btn-block", style = "margin:10px")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box { border-radius: 8px; }
      .info-box { border-radius: 6px; }
    "))),

    tabItems(
      # ----------------------------------------------------------------
      # TAB 1: Patient Profile
      # ----------------------------------------------------------------
      tabItem(tabName = "profile",
        fluidRow(
          box(title = "Disease Overview", width = 4, status = "primary", solidHeader = TRUE,
            h4("Celiac Disease (CD)"),
            p("Celiac disease is an immune-mediated systemic illness triggered by gluten
              in genetically susceptible individuals (HLA-DQ2/DQ8). It causes small
              intestinal villous atrophy, leading to malabsorption."),
            p(strong("Prevalence: "), "~1% worldwide (Western populations)"),
            p(strong("Genetic Risk: "), "HLA-DQ2.5 (~95%), HLA-DQ8 (~5%)"),
            p(strong("Gold-standard Dx: "), "Duodenal biopsy (Marsh ≥2) + positive anti-tTG IgA"),
            p(strong("Standard Tx: "), "Strict Gluten-Free Diet (GFD)"),
            hr(),
            h5("Current Patient Settings"),
            verbatimTextOutput("patient_summary")
          ),
          box(title = "Mechanistic Cascade", width = 8, status = "info", solidHeader = TRUE,
            fluidRow(
              valueBoxOutput("vb_marsh",  width = 3),
              valueBoxOutput("vb_antittg", width = 3),
              valueBoxOutput("vb_hgb",    width = 3),
              valueBoxOutput("vb_bmd",    width = 3)
            ),
            plotlyOutput("profile_radar", height = 350)
          )
        ),
        fluidRow(
          box(title = "Marsh Classification Guide", width = 6, status = "warning",
            tableOutput("marsh_table")
          ),
          box(title = "Drug Mechanism Summary", width = 6, status = "success",
            tableOutput("drug_moa_table")
          )
        )
      ),

      # ----------------------------------------------------------------
      # TAB 2: PK
      # ----------------------------------------------------------------
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "PK Parameters", width = 4, status = "primary", solidHeader = TRUE,
            sliderInput("dose_mg", "Drug Dose (mg/dose)", 0.5, 500, 2, step = 0.5),
            numericInput("F_oral",  "Oral Bioavailability (F)", 0.5, min=0.01, max=1, step=0.05),
            numericInput("ka_drug", "Absorption Rate ka (h⁻¹)", 0.6, min=0.1, max=3, step=0.1),
            numericInput("CL_drug", "Clearance CL (L/h)", 8, min=0.5, max=100, step=0.5),
            numericInput("Vd_drug", "Volume of Distribution Vd (L)", 6, min=1, max=200, step=1),
            numericInput("dose_interval", "Dosing Interval (h)", 8, min=4, max=168, step=4)
          ),
          box(title = "Plasma Concentration–Time Profile", width = 8, status = "info", solidHeader = TRUE,
            withSpinner(plotlyOutput("pk_plot", height = 400))
          )
        ),
        fluidRow(
          box(title = "PK Parameter Table", width = 12, status = "success",
            DTOutput("pk_table")
          )
        )
      ),

      # ----------------------------------------------------------------
      # TAB 3: Immune Biomarkers
      # ----------------------------------------------------------------
      tabItem(tabName = "immune",
        fluidRow(
          box(title = "Innate Immune — IL-15 & IEL", width = 6, status = "danger", solidHeader = TRUE,
            withSpinner(plotlyOutput("il15_plot", height = 300))
          ),
          box(title = "IEL Count (per 100 enterocytes)", width = 6, status = "danger", solidHeader = TRUE,
            withSpinner(plotlyOutput("iel_plot", height = 300))
          )
        ),
        fluidRow(
          box(title = "Th1 Cytokines — IFN-γ & TNF-α", width = 4, status = "warning", solidHeader = TRUE,
            withSpinner(plotlyOutput("ifng_plot", height = 300))
          ),
          box(title = "Th17 Cytokines — IL-17A & IL-21", width = 4, status = "warning", solidHeader = TRUE,
            withSpinner(plotlyOutput("il17_plot", height = 300))
          ),
          box(title = "CD4+ T Cells & B Cells", width = 4, status = "info", solidHeader = TRUE,
            withSpinner(plotlyOutput("tcell_plot", height = 300))
          )
        )
      ),

      # ----------------------------------------------------------------
      # TAB 4: Histopathology
      # ----------------------------------------------------------------
      tabItem(tabName = "histo",
        fluidRow(
          box(title = "Villous Height (VH, normalized)", width = 4, status = "warning", solidHeader = TRUE,
            withSpinner(plotlyOutput("vh_plot", height = 300))
          ),
          box(title = "Crypt Depth (normalized)", width = 4, status = "warning", solidHeader = TRUE,
            withSpinner(plotlyOutput("cd_plot", height = 300))
          ),
          box(title = "Villi:Crypt Ratio", width = 4, status = "warning", solidHeader = TRUE,
            withSpinner(plotlyOutput("vc_plot", height = 300))
          )
        ),
        fluidRow(
          box(title = "Marsh Score Trajectory", width = 7, status = "danger", solidHeader = TRUE,
            withSpinner(plotlyOutput("marsh_plot", height = 350))
          ),
          box(title = "Intestinal Permeability Index", width = 5, status = "info", solidHeader = TRUE,
            withSpinner(plotlyOutput("ip_plot", height = 350)),
            p("IP Index: 1=Normal baseline; higher = leaky gut")
          )
        )
      ),

      # ----------------------------------------------------------------
      # TAB 5: Clinical Endpoints
      # ----------------------------------------------------------------
      tabItem(tabName = "clinical",
        fluidRow(
          box(title = "Anti-tTG IgA Serology (U/mL)", width = 6, status = "primary", solidHeader = TRUE,
            withSpinner(plotlyOutput("antittg_plot", height = 300)),
            p("Positive threshold: >10 U/mL (ULN). Gold-standard non-invasive biomarker.")
          ),
          box(title = "Hemoglobin (g/dL, Iron Deficiency Anemia)", width = 6, status = "danger", solidHeader = TRUE,
            withSpinner(plotlyOutput("hgb_plot", height = 300)),
            p("Normal: Female ≥12 g/dL, Male ≥13 g/dL. Main consequence of iron malabsorption.")
          )
        ),
        fluidRow(
          box(title = "Bone Mineral Density (normalized BMD)", width = 6, status = "warning", solidHeader = TRUE,
            withSpinner(plotlyOutput("bmd_plot", height = 300)),
            p("BMD=1 normal; T-score < -2.5 = osteoporosis. Due to Ca/Vit D malabsorption.")
          ),
          box(title = "Intestinal Absorption Surface Area", width = 6, status = "success", solidHeader = TRUE,
            withSpinner(plotlyOutput("abs_plot", height = 300)),
            p("AbsArea=1 normal; reduced by villous atrophy. Governs all nutrient absorption.")
          )
        ),
        fluidRow(
          box(title = "1-Year Endpoint Summary Table", width = 12, status = "info",
            DTOutput("endpoint_table")
          )
        )
      ),

      # ----------------------------------------------------------------
      # TAB 6: Scenario Comparison
      # ----------------------------------------------------------------
      tabItem(tabName = "scenarios",
        fluidRow(
          box(title = "Predefined Treatment Scenarios", width = 12, status = "primary", solidHeader = TRUE,
            p("Comparing 6 standard treatment strategies. Click 'Run Simulation' to refresh."),
            withSpinner(plotlyOutput("scenario_antittg", height = 300))
          )
        ),
        fluidRow(
          box(title = "V:C Ratio by Scenario", width = 6, status = "warning", solidHeader = TRUE,
            withSpinner(plotlyOutput("scenario_vc", height = 300))
          ),
          box(title = "Hemoglobin by Scenario", width = 6, status = "danger", solidHeader = TRUE,
            withSpinner(plotlyOutput("scenario_hgb", height = 300))
          )
        ),
        fluidRow(
          box(title = "Scenario Results at 1 Year", width = 12, status = "success",
            DTOutput("scenario_table")
          )
        )
      )
    )
  )
)

# ==============================================================================
# SERVER
# ==============================================================================
server <- function(input, output, session) {

  # Reactive: run simulation when button pressed
  sim_result <- eventReactive(input$run_sim, {
    gfd_val  <- as.numeric(input$gfd)
    leak_val <- if (input$gfd == "partial") 0.05 else 0.0
    if (input$gfd == 1) gfd_val <- 1

    drug_type <- as.integer(input$drug_type)
    dose_interval_h <- switch(as.character(drug_type),
      "1" = as.numeric(input$dose_interval),
      "2" = 24,
      "3" = 168,
      8
    )

    params_list <- list(
      GlutenIntake_g_day = input$gluten_g,
      GFD                = gfd_val,
      GFD_leak           = leak_val,
      Drug_type          = drug_type,
      F_oral             = input$F_oral,
      ka_drug            = input$ka_drug,
      CL_drug_L_h        = input$CL_drug,
      Vd_drug_L          = input$Vd_drug
    )

    run_sim(params_list, drug_type, input$dose_mg, dose_interval_h, input$sim_days)
  }, ignoreNULL = FALSE)

  # Scenario comparison (6 predefined)
  scenario_result <- eventReactive(input$run_sim, {
    scen_defs <- list(
      list(GFD=0, leak=0, drug=0, dose=0,   lbl="①No GFD, No Drug",    col="#E53935"),
      list(GFD=1, leak=0, drug=0, dose=0,   lbl="②Strict GFD",         col="#43A047"),
      list(GFD=1, leak=0.05, drug=0, dose=0,lbl="③Partial GFD (5%)",   col="#FB8C00"),
      list(GFD=1, leak=0.10, drug=1, dose=2,lbl="④GFD+Larazotide",     col="#1E88E5"),
      list(GFD=1, leak=0.10, drug=2, dose=300, lbl="⑤GFD+ZED1227",    col="#8E24AA"),
      list(GFD=1, leak=0.20, drug=3, dose=150, lbl="⑥GFD+AMG714 RCD",col="#00897B")
    )
    lapply(scen_defs, function(s) {
      di <- switch(as.character(s$drug), "1"=8, "2"=24, "3"=168, 8)
      run_sim(
        list(GlutenIntake_g_day=12, GFD=s$GFD, GFD_leak=s$leak, Drug_type=s$drug),
        s$drug, s$dose, di, input$sim_days
      ) %>% mutate(Scenario = s$lbl, Color = s$col)
    }) %>% bind_rows()
  }, ignoreNULL = FALSE)

  # Patient summary text
  output$patient_summary <- renderText({
    d <- sim_result()
    last <- d[nrow(d),]
    paste0(
      "HLA Type:    ", input$disease_hla, "\n",
      "Gluten:      ", input$gluten_g, " g/day\n",
      "GFD:         ", input$gfd, "\n",
      "Drug:        ", input$drug_type, "\n",
      "Sim Days:    ", input$sim_days, "\n",
      "--- 1-Year Outcomes ---\n",
      "Anti-tTG:    ", round(last$AntiTTG, 1), " U/mL",
        if (last$AntiTTG > 10) " [POSITIVE]" else " [negative]", "\n",
      "Marsh Score: ", last$Marsh_score, "\n",
      "V:C Ratio:   ", round(last$VH_CD_ratio, 2), "\n",
      "Hgb:         ", round(last$Hgb_g_dL, 1), " g/dL\n",
      "BMD:         ", round(last$BMD, 3)
    )
  })

  # Value boxes
  output$vb_marsh <- renderValueBox({
    d <- sim_result(); last <- d[nrow(d),]
    valueBox(paste("Marsh", last$Marsh_score), "Histopathology",
             color = c("green","yellow","orange","red")[last$Marsh_score+1],
             icon  = icon("layer-group"))
  })
  output$vb_antittg <- renderValueBox({
    d <- sim_result(); last <- d[nrow(d),]
    valueBox(paste0(round(last$AntiTTG,1)," U/mL"), "Anti-tTG IgA",
             color = if (last$AntiTTG > 10) "red" else "green",
             icon  = icon("vial"))
  })
  output$vb_hgb <- renderValueBox({
    d <- sim_result(); last <- d[nrow(d),]
    valueBox(paste0(round(last$Hgb_g_dL,1)," g/dL"), "Hemoglobin",
             color = if (last$Hgb_g_dL < 12) "red" else "green",
             icon  = icon("tint"))
  })
  output$vb_bmd <- renderValueBox({
    d <- sim_result(); last <- d[nrow(d),]
    valueBox(round(last$BMD,3), "BMD (normalized)",
             color = if (last$BMD < 0.9) "orange" else "green",
             icon  = icon("bone"))
  })

  # Radar chart (profile)
  output$profile_radar <- renderPlotly({
    d <- sim_result(); last <- d[nrow(d),]
    categories <- c("Serology\n(anti-tTG)", "Histopathology\n(Marsh)", "Immune\n(IFN-γ)",
                    "Malabsorption\n(1-Abs)", "BMD\nLoss", "Iron\nAnemia")
    vals <- c(
      min(last$AntiTTG / 30, 1),
      last$Marsh_score / 3,
      min(last$IFNg / 15, 1),
      1 - last$AbsArea,
      1 - last$BMD,
      1 - last$Hgb_g_dL / 14
    )
    plot_ly(
      type = "scatterpolar", mode = "lines+markers",
      r = c(vals, vals[1]), theta = c(categories, categories[1]),
      fill = "toself", fillcolor = "rgba(30,136,229,0.3)",
      line = list(color = "#1E88E5")
    ) %>%
      layout(polar = list(radialaxis = list(range = c(0, 1))),
             title = "Disease Severity Profile (0=normal, 1=max)")
  })

  # Marsh table
  output$marsh_table <- renderTable({
    tibble(
      `Marsh Score` = c("0","1","2","3a","3b","3c"),
      `Histology` = c("Normal","↑IEL (>25/100)","+ Crypt hyperplasia",
                      "Partial villous atrophy","Subtotal villous atrophy",
                      "Total villous atrophy"),
      `V:C Ratio` = c(">3:1","≥2:1","≥1:1","1–2:1","<1:1","~0"),
      `Absorption` = c("Normal","Mild↓","Moderate↓","↓↓","↓↓↓","Minimal")
    )
  })

  output$drug_moa_table <- renderTable({
    tibble(
      Drug = c("GFD","Larazotide","ZED1227","AMG714","Budesonide","Azathioprine"),
      Target = c("Gluten (dietary)","Zonulin/TJ","tTG2 enzyme","IL-15","GR (local)","IMPDH"),
      Mechanism = c("Removes antigen source","Blocks tight junction opening",
                    "Prevents gliadin deamidation","Neutralizes innate IL-15 signal",
                    "Anti-inflammatory (steroid)","Immunosuppression"),
      Phase = c("Standard","Phase III","Phase II","Phase II","Approved (RCD)","Off-label")
    )
  })

  # PK plot
  output$pk_plot <- renderPlotly({
    d <- sim_result()
    p <- ggplot(d, aes(Time_days, DrugPlasma)) +
      geom_line(color = "#1565C0", linewidth = 1.5) +
      labs(title = "Drug Plasma Concentration vs Time",
           x = "Time (days)", y = "Plasma Cp (ng/mL)") +
      theme_bw(base_size = 13)
    ggplotly(p)
  })

  output$pk_table <- renderDT({
    d <- sim_result()
    tbl <- tibble(
      Parameter = c("F_oral", "ka (h⁻¹)", "CL (L/h)", "Vd (L)",
                    "t½ (h)", "Cmax (ng/mL)", "AUC0-24 (ng·h/mL)"),
      Value = c(
        input$F_oral, input$ka_drug, input$CL_drug, input$Vd_drug,
        round(0.693 / (input$CL_drug / input$Vd_drug), 1),
        round(max(d$DrugPlasma, na.rm=TRUE), 2),
        round(sum(d$DrugPlasma[d$Time_days <= 1], na.rm=TRUE) * 24, 1)
      )
    )
    datatable(tbl, options = list(dom="t", pageLength=10))
  })

  # Immune plots
  mk_plotly <- function(d, y_col, y_lab, title, ref_line = NULL, ref_label = NULL, color="#E53935") {
    p <- ggplot(d, aes_string("Time_days", y_col)) +
      geom_line(color = color, linewidth = 1.4) +
      labs(title = title, x = "Time (days)", y = y_lab) +
      theme_bw(base_size = 12)
    if (!is.null(ref_line)) {
      p <- p + geom_hline(yintercept = ref_line, linetype = "dashed", color = "gray50") +
        annotate("text", x = max(d$Time_days)*0.6, y = ref_line*1.05,
                 label = ref_label, size = 3, color = "gray40")
    }
    ggplotly(p)
  }

  output$il15_plot   <- renderPlotly({ d <- sim_result(); mk_plotly(d,"IL15","ng/mL","IL-15","#1E88E5") })
  output$iel_plot    <- renderPlotly({ d <- sim_result(); mk_plotly(d,"IEL","IEL/100 enterocytes","IEL Count", ref_line=25, ref_label="Marsh 1 threshold", color="#C62828") })
  output$ifng_plot   <- renderPlotly({ d <- sim_result(); mk_plotly(d,"IFNg","ng/mL","IFN-γ","#F57F17") })
  output$il17_plot   <- renderPlotly({ d <- sim_result(); mk_plotly(d,"IL17","ng/mL","IL-17A","#7B1FA2") })
  output$tcell_plot  <- renderPlotly({
    d <- sim_result()
    plot_ly(d, x=~Time_days, y=~CD4T, name="CD4+ T cells", type="scatter", mode="lines",
            line=list(color="#9C27B0")) %>%
      add_lines(y=~Bcell, name="B cells", line=list(color="#2196F3")) %>%
      layout(title="T & B Cell Dynamics", xaxis=list(title="Time (days)"),
             yaxis=list(title="Cell count (AU)"), legend=list(x=0.7,y=0.9))
  })

  # Histopathology plots
  output$vh_plot  <- renderPlotly({ d <- sim_result(); mk_plotly(d,"VH","VH (normalized)","Villous Height","#FF6F00") })
  output$cd_plot  <- renderPlotly({ d <- sim_result(); mk_plotly(d,"CrD","CrD (normalized)","Crypt Depth","#795548") })
  output$vc_plot  <- renderPlotly({
    d <- sim_result()
    p <- ggplot(d, aes(Time_days, VH_CD_ratio)) +
      geom_line(color="#FF6F00", linewidth=1.4) +
      geom_hline(yintercept=3, linetype="dashed", color="steelblue") +
      geom_hline(yintercept=1, linetype="dashed", color="red") +
      labs(title="Villi:Crypt Ratio", x="Time (days)", y="V:C Ratio") +
      theme_bw(base_size=12)
    ggplotly(p)
  })
  output$marsh_plot <- renderPlotly({
    d <- sim_result()
    p <- ggplot(d, aes(Time_days, Marsh_score)) +
      geom_step(color="#D32F2F", linewidth=1.6) +
      scale_y_continuous(breaks=0:3, labels=c("0 Normal","1 ↑IEL","2 Crypt","3 Atrophy")) +
      labs(title="Marsh Score Over Time", x="Time (days)", y="Marsh Score") +
      theme_bw(base_size=12)
    ggplotly(p)
  })
  output$ip_plot <- renderPlotly({ d <- sim_result(); mk_plotly(d,"IP","IP Index","Intestinal Permeability","#0288D1") })

  # Clinical endpoint plots
  output$antittg_plot <- renderPlotly({
    d <- sim_result()
    p <- ggplot(d, aes(Time_days, AntiTTG)) +
      geom_line(color="#1565C0", linewidth=1.4) +
      geom_hline(yintercept=10, linetype="dashed", color="red") +
      annotate("text", x=max(d$Time_days)*0.6, y=11.5, label="Positive (>10 ULN)",
               color="red", size=3.5) +
      labs(title="Anti-tTG IgA", x="Time (days)", y="U/mL") +
      theme_bw(base_size=12)
    ggplotly(p)
  })
  output$hgb_plot <- renderPlotly({
    d <- sim_result()
    p <- ggplot(d, aes(Time_days, Hgb_g_dL)) +
      geom_line(color="#C62828", linewidth=1.4) +
      geom_hline(yintercept=12, linetype="dashed", color="tomato") +
      annotate("text", x=max(d$Time_days)*0.6, y=11.5, label="Anemia threshold",
               color="tomato", size=3.5) +
      labs(title="Hemoglobin (g/dL)", x="Time (days)", y="Hgb (g/dL)") +
      theme_bw(base_size=12)
    ggplotly(p)
  })
  output$bmd_plot <- renderPlotly({ d <- sim_result(); mk_plotly(d,"BMD","BMD (normalized)","Bone Mineral Density","#5D4037") })
  output$abs_plot <- renderPlotly({ d <- sim_result(); mk_plotly(d,"AbsArea","Absorption Area (normalized)","Intestinal Absorption Surface","#2E7D32") })

  output$endpoint_table <- renderDT({
    d <- sim_result()
    last <- d[nrow(d),]
    tbl <- tibble(
      Endpoint = c("Anti-tTG IgA (U/mL)","Marsh Score","V:C Ratio","IFN-γ (ng/mL)",
                   "IL-15 (ng/mL)","IEL Count","Hemoglobin (g/dL)","BMD","AbsArea","Iron Stores (mg)"),
      Value    = c(round(last$AntiTTG,1), last$Marsh_score, round(last$VH_CD_ratio,2),
                   round(last$IFNg,2), round(last$IL15,2), round(last$IEL,1),
                   round(last$Hgb_g_dL,1), round(last$BMD,3), round(last$AbsArea,3),
                   round(last$IronStores,0)),
      Normal   = c("<10","0",">3","<1","<0.5","<25",">12","~1.0","~1.0","800-1200"),
      Status   = c(
        if (last$AntiTTG > 10) "⚠ ABNORMAL" else "✓ Normal",
        if (last$Marsh_score > 0) "⚠ ABNORMAL" else "✓ Normal",
        if (last$VH_CD_ratio < 1) "⚠ Atrophy" else "✓ Normal",
        if (last$IFNg > 5) "⚠ Elevated" else "✓ Normal",
        if (last$IL15 > 3) "⚠ Elevated" else "✓ Normal",
        if (last$IEL > 25) "⚠ Elevated" else "✓ Normal",
        if (last$Hgb_g_dL < 12) "⚠ Anemia" else "✓ Normal",
        if (last$BMD < 0.9) "⚠ ↓BMD" else "✓ Normal",
        if (last$AbsArea < 0.7) "⚠ Malabsorption" else "✓ Normal",
        if (last$IronStores < 500) "⚠ Iron-deficient" else "✓ Normal"
      )
    )
    datatable(tbl, options=list(dom="t",pageLength=12),
              rownames=FALSE) %>%
      formatStyle("Status",
        color = styleEqual(c("✓ Normal","⚠ ABNORMAL","⚠ Atrophy","⚠ Elevated","⚠ Anemia","⚠ ↓BMD","⚠ Malabsorption","⚠ Iron-deficient"),
                           c("green","red","red","orange","red","orange","orange","orange")))
  })

  # Scenario comparison plots
  scen_color_map <- c(
    "①No GFD, No Drug"="#E53935", "②Strict GFD"="#43A047",
    "③Partial GFD (5%)"="#FB8C00", "④GFD+Larazotide"="#1E88E5",
    "⑤GFD+ZED1227"="#8E24AA", "⑥GFD+AMG714 RCD"="#00897B"
  )
  output$scenario_antittg <- renderPlotly({
    d <- scenario_result()
    p <- ggplot(d, aes(Time_days, AntiTTG, color=Scenario)) +
      geom_line(linewidth=1.2) +
      geom_hline(yintercept=10, linetype="dashed", color="gray30") +
      scale_color_manual(values=scen_color_map) +
      labs(title="Anti-tTG IgA by Treatment Scenario", x="Time (days)", y="U/mL", color=NULL) +
      theme_bw(base_size=12) + theme(legend.position="bottom")
    ggplotly(p)
  })
  output$scenario_vc <- renderPlotly({
    d <- scenario_result()
    p <- ggplot(d, aes(Time_days, VH_CD_ratio, color=Scenario)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=scen_color_map) +
      labs(title="V:C Ratio by Scenario", x="Time (days)", y="V:C Ratio", color=NULL) +
      theme_bw(base_size=12) + theme(legend.position="bottom")
    ggplotly(p)
  })
  output$scenario_hgb <- renderPlotly({
    d <- scenario_result()
    p <- ggplot(d, aes(Time_days, Hgb_g_dL, color=Scenario)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=scen_color_map) +
      labs(title="Hemoglobin by Scenario", x="Time (days)", y="g/dL", color=NULL) +
      theme_bw(base_size=12) + theme(legend.position="bottom")
    ggplotly(p)
  })
  output$scenario_table <- renderDT({
    d <- scenario_result()
    tbl <- d %>%
      group_by(Scenario) %>%
      slice_tail(n=1) %>%
      select(Scenario, AntiTTG, VH_CD_ratio, Marsh_score, IFNg, Hgb_g_dL, BMD, AbsArea) %>%
      mutate(across(where(is.numeric), ~round(.x, 2))) %>%
      rename(
        `Anti-tTG IgA (U/mL)` = AntiTTG,
        `V:C Ratio` = VH_CD_ratio,
        `Marsh Score` = Marsh_score,
        `IFN-γ (ng/mL)` = IFNg,
        `Hgb (g/dL)` = Hgb_g_dL,
        `BMD` = BMD,
        `AbsArea` = AbsArea
      )
    datatable(tbl, options=list(dom="t", pageLength=8), rownames=FALSE)
  })
}

# ==============================================================================
shinyApp(ui = ui, server = server)
