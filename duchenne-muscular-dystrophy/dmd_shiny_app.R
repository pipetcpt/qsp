## ============================================================
## Duchenne Muscular Dystrophy (DMD) QSP — Shiny Interactive Dashboard
## ============================================================
## Interactive simulation dashboard for DMD QSP model
## Tabs: Patient Profile | Drug PK | Dystrophin Dynamics |
##       Motor Function | Clinical Endpoints | Scenario Comparison | Biomarkers
##
## Usage: shiny::runApp("dmd_shiny_app.R")
## Requirements: shiny, mrgsolve, dplyr, ggplot2, plotly, bslib, DT, patchwork
## ============================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(bslib)
library(DT)
library(tidyr)

## ─────────────────────────────────────────────────────────────
## MODEL DEFINITION (inline)
## ─────────────────────────────────────────────────────────────

dmd_code <- '
$PARAM
CL_Ete=80.4, Vc_Ete=15.8, Vp_Ete=42.2, Q_Ete=30.1,
kUptake=0.035, kElim_m=0.005,
ka_DFZ=1.20, CL_DFZ=18.5, Vd_DFZ=62.0, F_DFZ=0.89,
kconv_DFZ=2.10, CL_aDFZ=8.20, Vd_aDFZ=45.0,
EC50_skip=0.45, Emax_skip=0.85, Hill_skip=1.5,
kDyst_syn=0.0045, kDyst_deg=0.0040,
GT_contribution=0.0,
EC50_GC=0.020, Emax_GC=0.75, Hill_GC=1.2,
kNecrosis=0.0025, kRegen=0.0060, kMaturation=0.0020,
kFibrosis_rate=0.00015,
Inflam_basal=25.0, kInflam_stim=0.15, kInflam_decay=0.010,
kFib_stim=0.0025, kFib_decay=0.00008,
Fibrosis_max=0.85,
kSC_regen=0.005, kSC_exhaust=0.0003,
kCK_release=0.08, kCK_elim=0.025, CK_basal=150,
kFVC_decline=0.000095, FVC_min=10.0,
kLVEF_decline=0.000050, LVEF_min=20.0,
k_NSAA_decline=0.00015, k_6MWD_decline=0.0004,
BW=20.0, Fiber_H_0=100.0,
GT_eff=0.0

$INIT
Ete_C1=0, Ete_C2=0, Ete_Muscle=0,
DFZ_Gut=0, DFZ_Plasma=0, Active_DFZ=0,
Dystrophin=0, Fiber_H=100.0, Fiber_N=0.001,
Fiber_R=0.0, Inflam=25.0, Fibrosis=0.05,
SC_Pool=1.0, CK_serum=15000,
FVC_pct=95.0, LVEF=62.0, NSAA=28.0, SixMWD=380.0

$ODE
double Cp_Ete = Ete_C1/Vc_Ete;
double Ct_Ete = Ete_C2/Vp_Ete;
dxdt_Ete_C1 = -CL_Ete*Cp_Ete - Q_Ete*(Cp_Ete-Ct_Ete) - kUptake*Ete_C1;
dxdt_Ete_C2 = Q_Ete*(Cp_Ete-Ct_Ete);
dxdt_Ete_Muscle = kUptake*Ete_C1 - kElim_m*Ete_Muscle;
dxdt_DFZ_Gut = -ka_DFZ*DFZ_Gut;
double Cp_DFZ  = DFZ_Plasma/Vd_DFZ;
double Cp_aDFZ = Active_DFZ/Vd_aDFZ;
dxdt_DFZ_Plasma = F_DFZ*ka_DFZ*DFZ_Gut - CL_DFZ*Cp_DFZ - kconv_DFZ*DFZ_Plasma;
dxdt_Active_DFZ = kconv_DFZ*DFZ_Plasma - CL_aDFZ*Cp_aDFZ;
double ASO_pmol = Ete_Muscle/1000.0;
double Eskip = Emax_skip*pow(ASO_pmol,Hill_skip)/(pow(EC50_skip,Hill_skip)+pow(ASO_pmol,Hill_skip));
double kDyst_syn_eff = kDyst_syn*(0.02 + Eskip + GT_contribution);
double GC_NF_inhib = Emax_GC*pow(Cp_aDFZ,Hill_GC)/(pow(EC50_GC,Hill_GC)+pow(Cp_aDFZ,Hill_GC));
double GC_protect = 1.0 - GC_NF_inhib;
double SC_regen_eff = kRegen*SC_Pool;
dxdt_Dystrophin = kDyst_syn_eff - kDyst_deg*Dystrophin;
double dyst_protect = 1.0 - 1.0/(1.0 + Dystrophin/10.0);
double Necrosis_rate = kNecrosis*Inflam*GC_protect*(1.0-dyst_protect)*Fiber_H/100.0;
double Maturation_rate = kMaturation*Fiber_R*SC_Pool;
dxdt_Fiber_H = -Necrosis_rate + Maturation_rate;
dxdt_Fiber_N = Necrosis_rate - SC_regen_eff*Fiber_N;
dxdt_Fiber_R = SC_regen_eff*Fiber_N - Maturation_rate;
dxdt_Inflam = kInflam_stim*Fiber_N - kInflam_decay*Inflam - GC_NF_inhib*kInflam_decay*Inflam;
double Fibrosis_stim = kFib_stim*Inflam/100.0*(1.0-Fibrosis/Fibrosis_max);
double GC_anti_fib = GC_NF_inhib*0.5;
dxdt_Fibrosis = (1.0-GC_anti_fib)*Fibrosis_stim - kFib_decay*Fibrosis;
dxdt_SC_Pool = kSC_regen*(1.0-SC_Pool) - kSC_exhaust*Fiber_N*SC_Pool;
dxdt_CK_serum = CK_basal + kCK_release*Fiber_N - kCK_elim*CK_serum;
dxdt_FVC_pct = -kFVC_decline*Fibrosis*100.0*(FVC_pct-FVC_min);
dxdt_LVEF = -kLVEF_decline*Fibrosis*100.0*(LVEF-LVEF_min);
dxdt_NSAA = -k_NSAA_decline*Fibrosis*100.0*(NSAA-0.0);
dxdt_SixMWD = -k_6MWD_decline*Fibrosis*100.0*(SixMWD-0.0);

$TABLE
double Cp_Ete_ugmL = Ete_C1/Vc_Ete;
double aDFZ_Cp = Active_DFZ/Vd_aDFZ;
double Exon_Skip_pct = Emax_skip*pow(Ete_Muscle/1000.0,Hill_skip)/
  (pow(EC50_skip,Hill_skip)+pow(Ete_Muscle/1000.0,Hill_skip))*100.0;
double Dyst_pct = Dystrophin;
double GC_Inhib_pct = Emax_GC*pow(aDFZ_Cp,Hill_GC)/
  (pow(EC50_GC,Hill_GC)+pow(aDFZ_Cp,Hill_GC))*100.0;
double Fibrosis_pct = Fibrosis*100.0;
double Inflam_idx = Inflam;
double SC_pct = SC_Pool*100.0;
double CK_kUL = CK_serum/1000.0;
double Ambulatory = (SixMWD>10.0)?1.0:0.0;

$CAPTURE
Cp_Ete_ugmL aDFZ_Cp Exon_Skip_pct
Dyst_pct GC_Inhib_pct Fibrosis_pct Inflam_idx SC_pct
Fiber_H Fiber_N Fiber_R
CK_serum CK_kUL FVC_pct LVEF NSAA SixMWD Ambulatory
'

mod_base <- mcode("dmd_shiny", dmd_code, quiet=TRUE)

## ─────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────

ui <- page_navbar(
  title = "DMD QSP Dashboard",
  theme = bs_theme(
    bootswatch = "darkly",
    primary    = "#9b59b6",
    secondary  = "#3498db",
    success    = "#2ecc71",
    danger     = "#e74c3c",
    warning    = "#f39c12"
  ),
  bg = "#1a1a2e",
  fillable = FALSE,

  ## ── Tab 1: Patient Profile ─────────────────────────────────
  nav_panel(
    title = "Patient Profile",
    icon  = icon("user"),
    fluidRow(
      column(4,
        card(
          card_header("Patient Demographics"),
          sliderInput("BW",     "Body Weight (kg):",  min=10, max=60, value=20, step=1),
          sliderInput("Age_st", "Starting Age (yr):", min=3,  max=10, value=5,  step=1),
          hr(),
          h5("Mutation Type"),
          radioButtons("mutation_type",
            "DMD Mutation:",
            choices = c(
              "Deletion/Duplication (exon-skippable, ~72%)" = "deletion",
              "Nonsense/PTC (~13%)"                          = "nonsense",
              "Point/Other (~15%)"                           = "other"
            ),
            selected = "deletion"
          ),
          conditionalPanel(
            "input.mutation_type == 'deletion'",
            radioButtons("exon_target",
              "Exon Skip Target:",
              choices = c(
                "Exon-51 (EXONDYS 51 / Eteplirsen, ~13% eligible)"  = "exon51",
                "Exon-45 (AMONDYS 45 / Casimersen, ~8% eligible)"   = "exon45",
                "Exon-53 (VYONDYS 53 / Golodirsen, ~8% eligible)"   = "exon53",
                "Other deletion (not skippable)"                     = "other_del"
              ),
              selected = "exon51"
            )
          )
        )
      ),
      column(4,
        card(
          card_header("Disease Severity at Baseline"),
          sliderInput("Dyst_base",    "Baseline Dystrophin (% normal):",   0, 5,    0,   0.5),
          sliderInput("Fiber_H_base", "Healthy Fibers (au, 0-100):",       50, 100, 100, 1),
          sliderInput("Fibrosis_base","Baseline Fibrosis (%):",            0,  30,  5,   1),
          sliderInput("NSAA_base",    "Baseline NSAA Score (0-34):",       5,  34,  28,  1),
          sliderInput("SixMWD_base",  "Baseline 6MWD (meters):",           50, 500, 380, 10),
          sliderInput("FVC_base",     "Baseline FVC% Predicted:",          40, 100, 95,  1),
          sliderInput("LVEF_base",    "Baseline LVEF (%):",                30, 70,  62,  1)
        )
      ),
      column(4,
        card(
          card_header("Clinical Background"),
          div(style="font-size:14px; color:#ccc;",
            p(strong("DMD Overview:"), "X-linked recessive neuromuscular disorder caused by mutations in the",
              em("DMD"), "gene (Xp21.2) encoding dystrophin (427 kDa)."),
            p(strong("Epidemiology:"), "~1 in 3,500 live male births; ~300,000 affected worldwide."),
            p(strong("Natural History:")),
            tags$ul(
              tags$li("Gross motor delay: ages 2-3yr"),
              tags$li("Gower sign, waddling gait: ages 3-6yr"),
              tags$li("Calf pseudohypertrophy"),
              tags$li("6MWD peaks age 8-9yr (~380-400m)"),
              tags$li("Loss of ambulation: median age 12-13yr (untreated)"),
              tags$li("Respiratory support: ~20yr"),
              tags$li("Cardiomyopathy: ~10-12yr onset"),
              tags$li("Median survival: ~26yr (untreated) → >40yr with modern Rx")
            ),
            p(strong("Serum CK:"), "Elevated 10-100× ULN even before symptoms"),
            p(strong("LVEF threshold:"), "<55% = Dilated CM (ACEi/ARB start)")
          )
        )
      )
    )
  ),

  ## ── Tab 2: Drug PK ─────────────────────────────────────────
  nav_panel(
    title = "Drug PK",
    icon  = icon("pills"),
    fluidRow(
      column(3,
        card(
          card_header("ASO (Exon-skipping) PK"),
          h6("Eteplirsen / Casimersen / Golodirsen"),
          sliderInput("dose_Ete",     "ASO Dose (mg/kg/wk):", 10, 60, 30, 5),
          sliderInput("CL_Ete",       "CL (L/h):",            40, 150, 80, 5),
          sliderInput("Vc_Ete",       "Vc (L):",              5,  40,  16, 1),
          sliderInput("PK_Vp_Ete",    "Vp (L):",              10, 80,  42, 2),
          hr(),
          h6("Deflazacort (EMFLAZA) PK"),
          sliderInput("dose_DFZ",     "DFZ Dose (mg/kg/d):",  0.3, 1.5, 0.9, 0.1),
          sliderInput("CL_DFZ_inp",   "CL_DFZ (L/h):",        5,  40,  18, 1),
          numericInput("sim_days_pk", "PK Sim Days:", 28, 1, 180, 7)
        )
      ),
      column(9,
        fluidRow(
          column(6, card(card_header("Eteplirsen Plasma PK"),
            plotlyOutput("plot_ete_pk", height="280px"))),
          column(6, card(card_header("Eteplirsen Muscle Uptake"),
            plotlyOutput("plot_ete_muscle", height="280px")))
        ),
        fluidRow(
          column(6, card(card_header("Deflazacort → 21-DFZ Plasma PK"),
            plotlyOutput("plot_dfz_pk", height="280px"))),
          column(6, card(card_header("GC Effect (NF-κB Inhibition)"),
            plotlyOutput("plot_gc_effect", height="280px")))
        )
      )
    )
  ),

  ## ── Tab 3: Dystrophin Dynamics ─────────────────────────────
  nav_panel(
    title = "Dystrophin",
    icon  = icon("dna"),
    fluidRow(
      column(3,
        card(
          card_header("Exon-Skipping Parameters"),
          sliderInput("Emax_skip",  "Emax Skip (max eff, 0-1):",   0.3, 1.0, 0.85, 0.05),
          sliderInput("EC50_skip",  "EC50 Skip (pmol/g):",         0.1, 2.0, 0.45, 0.05),
          sliderInput("Hill_skip",  "Hill Coefficient:",           0.5, 3.0, 1.5,  0.1),
          hr(),
          h6("Gene Therapy (Elevidys)"),
          checkboxInput("GT_active",  "Gene Therapy (Elevidys)", FALSE),
          sliderInput("GT_contribution_inp", "GT Dystrophin Boost (%):",
                      0, 80, 50, 5),
          hr(),
          h6("Dystrophin Kinetics"),
          sliderInput("kDyst_syn_inp", "kDyst_syn (1/h):", 0.001, 0.020, 0.0045, 0.0005),
          sliderInput("kDyst_deg_inp", "kDyst_deg (1/h):", 0.001, 0.020, 0.0040, 0.0005),
          numericInput("sim_days_dyst", "Simulation Days:", 365, 30, 730, 30)
        )
      ),
      column(9,
        fluidRow(
          column(6, card(card_header("Dystrophin Level Over Time"),
            plotlyOutput("plot_dystrophin", height="300px"))),
          column(6, card(card_header("Exon-Skipping Efficiency (dose-response)"),
            plotlyOutput("plot_eskip_dr", height="300px")))
        ),
        fluidRow(
          column(12, card(card_header("Muscle Fiber Compartments"),
            plotlyOutput("plot_fibers", height="300px")))
        )
      )
    )
  ),

  ## ── Tab 4: Motor Function ──────────────────────────────────
  nav_panel(
    title = "Motor Function",
    icon  = icon("person-walking"),
    fluidRow(
      column(3,
        card(
          card_header("Motor Function Parameters"),
          sliderInput("sim_years", "Simulation Duration (yr):", 5, 20, 8, 1),
          hr(),
          h6("Treatment Selection"),
          checkboxInput("use_dfz",  "Deflazacort (0.9 mg/kg/d)",      TRUE),
          checkboxInput("use_ete",  "Eteplirsen 30 mg/kg/wk (ex51)",  FALSE),
          checkboxInput("use_gt",   "Gene Therapy (Elevidys)",         FALSE),
          checkboxInput("use_vamo", "Vamorolone 6 mg/kg/d",           FALSE),
          hr(),
          h6("NSAA / 6MWD Parameters"),
          sliderInput("k_NSAA_inp",  "NSAA Decline Rate:", 0.00005, 0.0005, 0.00015, 0.00005),
          sliderInput("k_6MWD_inp",  "6MWD Decline Rate:", 0.00010, 0.0010, 0.00040, 0.00010)
        )
      ),
      column(9,
        fluidRow(
          column(6, card(card_header("6-Minute Walk Distance (6MWD)"),
            plotlyOutput("plot_6mwd", height="320px"))),
          column(6, card(card_header("NSAA Score"),
            plotlyOutput("plot_nsaa", height="320px")))
        ),
        fluidRow(
          column(6, card(card_header("Ambulatory Status Timeline"),
            plotlyOutput("plot_ambulatory", height="300px"))),
          column(6, card(card_header("Satellite Cell Pool"),
            plotlyOutput("plot_sc_pool", height="300px")))
        )
      )
    )
  ),

  ## ── Tab 5: Clinical Endpoints ──────────────────────────────
  nav_panel(
    title = "Clinical Endpoints",
    icon  = icon("heartbeat"),
    fluidRow(
      column(3,
        card(
          card_header("Endpoint Parameters"),
          sliderInput("sim_years_clin", "Simulation Duration (yr):", 5, 25, 15, 1),
          hr(),
          h6("Respiratory"),
          sliderInput("kFVC_inp", "FVC Decline Rate (/h):",
                      0.00003, 0.0003, 0.000095, 0.000005),
          hr(),
          h6("Cardiac"),
          sliderInput("kLVEF_inp", "LVEF Decline Rate (/h):",
                      0.00001, 0.0002, 0.000050, 0.000005),
          hr(),
          h6("Cardioprotection"),
          checkboxInput("use_ACEi", "ACEi/ARB (start LVEF<55%)", TRUE),
          checkboxInput("use_BB",   "Beta-blocker", FALSE)
        )
      ),
      column(9,
        fluidRow(
          column(6, card(card_header("FVC% Predicted (Respiratory)"),
            plotlyOutput("plot_fvc", height="300px"))),
          column(6, card(card_header("Left Ventricular EF (Cardiac)"),
            plotlyOutput("plot_lvef", height="300px")))
        ),
        fluidRow(
          column(6, card(card_header("Serum CK (Muscle Damage Marker)"),
            plotlyOutput("plot_ck", height="300px"))),
          column(6, card(card_header("Inflammation Index"),
            plotlyOutput("plot_inflam", height="300px")))
        )
      )
    )
  ),

  ## ── Tab 6: Scenario Comparison ────────────────────────────
  nav_panel(
    title = "Scenario Comparison",
    icon  = icon("chart-bar"),
    fluidRow(
      column(3,
        card(
          card_header("Comparison Settings"),
          numericInput("comp_years", "Simulation Duration (yr):", 8, 5, 20, 1),
          hr(),
          h6("Scenarios to compare:"),
          checkboxGroupInput("scenarios_sel",
            "Select Scenarios:",
            choices = list(
              "1. Natural History"                     = "S1",
              "2. Deflazacort 0.9 mg/kg/d"             = "S2",
              "3. Prednisone 0.75 mg/kg/d"             = "S3",
              "4. Eteplirsen+DFZ (ex51)"               = "S4",
              "5. Casimersen+DFZ (ex45)"               = "S5",
              "6. Gene Therapy (Elevidys)"             = "S6",
              "7. Vamorolone 6 mg/kg/d"                = "S7"
            ),
            selected = c("S1","S2","S4","S6")
          ),
          hr(),
          actionButton("run_comparison", "Run Comparison",
                       class="btn-primary", width="100%")
        )
      ),
      column(9,
        fluidRow(
          column(6, card(card_header("6MWD Comparison"),
            plotlyOutput("comp_6mwd", height="300px"))),
          column(6, card(card_header("NSAA Comparison"),
            plotlyOutput("comp_nsaa", height="300px")))
        ),
        fluidRow(
          column(6, card(card_header("Fibrosis Comparison"),
            plotlyOutput("comp_fibrosis", height="300px"))),
          column(6, card(card_header("FVC% Comparison"),
            plotlyOutput("comp_fvc", height="300px")))
        ),
        fluidRow(
          column(12,
            card(card_header("Summary Table (at Year 5 post-start)"),
              DTOutput("comp_table")
            )
          )
        )
      )
    )
  ),

  ## ── Tab 7: Biomarkers ─────────────────────────────────────
  nav_panel(
    title = "Biomarkers",
    icon  = icon("flask"),
    fluidRow(
      column(3,
        card(
          card_header("Biomarker Settings"),
          numericInput("bio_years", "Simulation Duration (yr):", 10, 5, 20, 1),
          hr(),
          selectInput("bio_treatment",
            "Select Treatment:",
            choices = list(
              "Untreated"                      = "untreated",
              "Deflazacort"                    = "dfz",
              "Eteplirsen + Deflazacort"       = "ete_dfz",
              "Gene Therapy (Elevidys)"        = "gt"
            ),
            selected = "dfz"
          ),
          hr(),
          h6("Key Biomarker Thresholds"),
          div(style="font-size:12px; color:#aaa;",
            p(icon("circle-xmark", style="color:red"), " CK >10,000 U/L: significant necrosis"),
            p(icon("circle-xmark", style="color:orange"), " FVC <50%: high exacerbation risk"),
            p(icon("circle-xmark", style="color:red"), " FVC <30%: nocturnal NIV indicated"),
            p(icon("circle-xmark", style="color:orange"), " LVEF <55%: DCM, start ACEi/ARB"),
            p(icon("circle-xmark", style="color:red"), " LVEF <40%: advanced heart failure"),
            p(icon("circle-check", style="color:green"), " Dyst >3%: potentially functional"),
            p(icon("circle-check", style="color:green"), " Dyst >15%: BMD-like phenotype")
          )
        )
      ),
      column(9,
        fluidRow(
          column(6, card(card_header("Serum CK Trajectory"),
            plotlyOutput("bio_ck", height="280px"))),
          column(6, card(card_header("Dystrophin Expression"),
            plotlyOutput("bio_dyst", height="280px")))
        ),
        fluidRow(
          column(6, card(card_header("Muscle Fiber Composition"),
            plotlyOutput("bio_fibers", height="280px"))),
          column(6, card(card_header("Satellite Cell Reserve"),
            plotlyOutput("bio_sc", height="280px")))
        ),
        fluidRow(
          column(12, card(card_header("Biomarker Summary Table"),
            DTOutput("bio_table")
          ))
        )
      )
    )
  )  # end biomarkers tab
)  # end page_navbar

## ─────────────────────────────────────────────────────────────
## SERVER
## ─────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  ## ── Helper: build model with user params ──────────────────
  make_mod <- function(extra_params = list()) {
    p_list <- c(
      list(BW        = input$BW,
           EC50_skip = input$EC50_skip,
           Emax_skip = input$Emax_skip,
           Hill_skip = input$Hill_skip,
           kDyst_syn = input$kDyst_syn_inp,
           kDyst_deg = input$kDyst_deg_inp),
      extra_params
    )
    tryCatch(param(mod_base, p_list), error = function(e) mod_base)
  }

  ## ── Helper: run simulation ────────────────────────────────
  run_sim <- function(model, events, sim_hours, delta = 168) {
    tryCatch({
      df <- as.data.frame(mrgsim(model, events=events, end=sim_hours, delta=delta))
      df$Age_yr <- input$Age_st + df$time/(365.25*24)
      df
    }, error = function(e) NULL)
  }

  ## ── Helper: build dosing events ──────────────────────────
  ev_build <- function(dose_ete_mg=0, dose_dfz_mg=0, sim_h=70080, gt=FALSE) {
    evs <- list()
    if (dose_ete_mg > 0)
      evs[[length(evs)+1]] <- ev(time=0, evid=1, cmt=1, amt=dose_ete_mg,
                                  ii=168, addl=floor(sim_h/168))
    if (dose_dfz_mg > 0)
      evs[[length(evs)+1]] <- ev(time=0, evid=1, cmt=4, amt=dose_dfz_mg/4,
                                  ii=6, addl=floor(sim_h/6))
    if (length(evs)==0) return(ev(time=0, cmt=1, amt=0))
    do.call(c, evs)
  }

  ## ─────────────────────────────────────
  ## TAB 2: DRUG PK
  ## ─────────────────────────────────────

  pk_sim <- reactive({
    sim_h    <- input$sim_days_pk * 24
    dose_ete <- input$dose_Ete * input$BW
    dose_dfz <- input$dose_DFZ * input$BW
    m <- param(mod_base, list(
      CL_Ete=input$CL_Ete, Vc_Ete=input$Vc_Ete, Vp_Ete=input$PK_Vp_Ete,
      CL_DFZ=input$CL_DFZ_inp, BW=input$BW
    ))
    events <- ev_build(dose_ete_mg=dose_ete, dose_dfz_mg=dose_dfz, sim_h=sim_h)
    run_sim(m, events, sim_h, delta=1)
  })

  output$plot_ete_pk <- renderPlotly({
    df <- pk_sim(); req(!is.null(df))
    plot_ly(df, x=~time/24, y=~Cp_Ete_ugmL, type='scatter', mode='lines',
            line=list(color='#9b59b6', width=2)) %>%
      layout(xaxis=list(title="Time (days)", color="white"),
             yaxis=list(title="Cp (µg/mL)", color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"))
  })

  output$plot_ete_muscle <- renderPlotly({
    df <- pk_sim(); req(!is.null(df))
    plot_ly(df, x=~time/24, y=~Ete_Muscle/1000, type='scatter', mode='lines',
            line=list(color='#8e44ad', width=2)) %>%
      layout(xaxis=list(title="Time (days)", color="white"),
             yaxis=list(title="ASO Muscle (pmol/g)", color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"))
  })

  output$plot_dfz_pk <- renderPlotly({
    df <- pk_sim(); req(!is.null(df))
    plot_ly(df, x=~time/24, y=~aDFZ_Cp, type='scatter', mode='lines',
            line=list(color='#d4ac0d', width=2),
            name="21-DFZ (active)") %>%
      add_trace(y=~DFZ_Plasma/62, name="DFZ (pro-drug)",
                line=list(color='#f39c12', width=2, dash='dash')) %>%
      layout(xaxis=list(title="Time (days)", color="white"),
             yaxis=list(title="Cp (mg/L)", color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$plot_gc_effect <- renderPlotly({
    df <- pk_sim(); req(!is.null(df))
    plot_ly(df, x=~time/24, y=~GC_Inhib_pct, type='scatter', mode='lines',
            line=list(color='#27ae60', width=2)) %>%
      layout(xaxis=list(title="Time (days)", color="white"),
             yaxis=list(title="NF-κB Inhibition (%)", range=c(0,100), color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"))
  })

  ## ─────────────────────────────────────
  ## TAB 3: DYSTROPHIN
  ## ─────────────────────────────────────

  dyst_sim <- reactive({
    sim_h <- input$sim_days_dyst * 24
    dose_ete <- 30 * input$BW
    dose_dfz <- 0.9 * input$BW
    gt_c <- if (input$GT_active) input$GT_contribution_inp/100 else 0

    m <- make_mod(list(GT_contribution=gt_c))
    events <- ev_build(dose_ete_mg=dose_ete, dose_dfz_mg=dose_dfz, sim_h=sim_h)
    df <- run_sim(m, events, sim_h, delta=24); req(!is.null(df))
    df
  })

  output$plot_dystrophin <- renderPlotly({
    df <- dyst_sim()
    plot_ly(df, x=~time/24, y=~Dyst_pct, type='scatter', mode='lines',
            line=list(color='#2ecc71', width=2)) %>%
      add_segments(x=0, xend=max(df$time/24), y=3, yend=3,
                   line=list(color='yellow', dash='dash', width=1),
                   name="3% threshold") %>%
      add_segments(x=0, xend=max(df$time/24), y=15, yend=15,
                   line=list(color='orange', dash='dot', width=1),
                   name="15% (BMD-like)") %>%
      layout(xaxis=list(title="Time (days)", color="white"),
             yaxis=list(title="Dystrophin (% normal)", color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$plot_eskip_dr <- renderPlotly({
    aso_conc <- seq(0, 3, by=0.05)
    skip_eff <- input$Emax_skip * aso_conc^input$Hill_skip /
                (input$EC50_skip^input$Hill_skip + aso_conc^input$Hill_skip) * 100
    plot_ly(x=aso_conc, y=skip_eff, type='scatter', mode='lines',
            line=list(color='#9b59b6', width=2)) %>%
      layout(xaxis=list(title="ASO Muscle Conc (pmol/g)", color="white"),
             yaxis=list(title="Exon-Skip Efficiency (%)", range=c(0,100), color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"))
  })

  output$plot_fibers <- renderPlotly({
    df <- dyst_sim()
    plot_ly(df, x=~time/24) %>%
      add_trace(y=~Fiber_H, name="Healthy", type='scatter', mode='lines',
                line=list(color='#2ecc71', width=2)) %>%
      add_trace(y=~Fiber_R, name="Regenerating", type='scatter', mode='lines',
                line=list(color='#f39c12', width=2)) %>%
      add_trace(y=~Fiber_N, name="Necrotic", type='scatter', mode='lines',
                line=list(color='#e74c3c', width=2)) %>%
      layout(xaxis=list(title="Time (days)", color="white"),
             yaxis=list(title="Muscle Fiber (au)", color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  ## ─────────────────────────────────────
  ## TAB 4: MOTOR FUNCTION
  ## ─────────────────────────────────────

  motor_sim <- reactive({
    sim_h <- input$sim_years * 365.25 * 24
    dose_ete <- if (input$use_ete) 30 * input$BW else 0
    dose_dfz <- if (input$use_dfz) 0.9 * input$BW else 0
    gt_c <- if (input$use_gt) 0.5 else 0

    p_extra <- list(
      k_NSAA_decline = input$k_NSAA_inp,
      k_6MWD_decline = input$k_6MWD_inp,
      GT_contribution = gt_c
    )
    if (input$use_vamo) {
      p_extra$EC50_GC <- 0.035
      p_extra$Emax_GC <- 0.72
      p_extra$kFib_stim <- 0.0022
    }
    m <- make_mod(p_extra)
    events <- ev_build(dose_ete_mg=dose_ete, dose_dfz_mg=dose_dfz, sim_h=sim_h)
    df <- run_sim(m, events, sim_h, delta=168); req(!is.null(df)); df
  })

  output$plot_6mwd <- renderPlotly({
    df <- motor_sim()
    plot_ly(df, x=~Age_yr, y=~SixMWD, type='scatter', mode='lines',
            line=list(color='#3498db', width=2)) %>%
      add_segments(x=min(df$Age_yr), xend=max(df$Age_yr), y=0, yend=0,
                   line=list(color='red', dash='dash', width=1), name="Loss of ambulation") %>%
      layout(xaxis=list(title="Age (years)", color="white"),
             yaxis=list(title="6MWD (meters)", color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$plot_nsaa <- renderPlotly({
    df <- motor_sim()
    plot_ly(df, x=~Age_yr, y=~NSAA, type='scatter', mode='lines',
            line=list(color='#9b59b6', width=2)) %>%
      layout(xaxis=list(title="Age (years)", color="white"),
             yaxis=list(title="NSAA Score (0-34)", range=c(0,34), color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"))
  })

  output$plot_ambulatory <- renderPlotly({
    df <- motor_sim()
    df$Ambulatory_Status <- ifelse(df$SixMWD > 10, "Ambulatory", "Non-ambulatory")
    plot_ly(df, x=~Age_yr, y=~Ambulatory, type='scatter', mode='lines',
            fill='tozeroy', fillcolor='rgba(46,204,113,0.3)',
            line=list(color='#2ecc71', width=2)) %>%
      layout(xaxis=list(title="Age (years)", color="white"),
             yaxis=list(title="Ambulatory (1=Yes, 0=No)", range=c(-0.1,1.2), color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"))
  })

  output$plot_sc_pool <- renderPlotly({
    df <- motor_sim()
    plot_ly(df, x=~Age_yr, y=~SC_pct, type='scatter', mode='lines',
            line=list(color='#27ae60', width=2)) %>%
      add_segments(x=min(df$Age_yr), xend=max(df$Age_yr), y=30, yend=30,
                   line=list(color='orange', dash='dash'), name="Critical threshold") %>%
      layout(xaxis=list(title="Age (years)", color="white"),
             yaxis=list(title="Satellite Cell Reserve (%)", range=c(0,105), color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  ## ─────────────────────────────────────
  ## TAB 5: CLINICAL ENDPOINTS
  ## ─────────────────────────────────────

  clin_sim <- reactive({
    sim_h <- input$sim_years_clin * 365.25 * 24
    dose_dfz <- 0.9 * input$BW
    m <- make_mod(list(kFVC_decline=input$kFVC_inp, kLVEF_decline=input$kLVEF_inp))
    events <- ev_build(dose_dfz_mg=dose_dfz, sim_h=sim_h)
    df <- run_sim(m, events, sim_h, delta=168); req(!is.null(df)); df
  })

  output$plot_fvc <- renderPlotly({
    df <- clin_sim()
    plot_ly(df, x=~Age_yr, y=~FVC_pct, type='scatter', mode='lines',
            line=list(color='#3498db', width=2)) %>%
      add_segments(x=min(df$Age_yr), xend=max(df$Age_yr), y=30, yend=30,
                   line=list(color='orange', dash='dash'), name="NIV threshold (30%)") %>%
      add_segments(x=min(df$Age_yr), xend=max(df$Age_yr), y=50, yend=50,
                   line=list(color='yellow', dash='dot'), name="50% threshold") %>%
      layout(xaxis=list(title="Age (years)", color="white"),
             yaxis=list(title="FVC% Predicted", range=c(0,105), color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$plot_lvef <- renderPlotly({
    df <- clin_sim()
    plot_ly(df, x=~Age_yr, y=~LVEF, type='scatter', mode='lines',
            line=list(color='#e74c3c', width=2)) %>%
      add_segments(x=min(df$Age_yr), xend=max(df$Age_yr), y=55, yend=55,
                   line=list(color='lightblue', dash='dash'), name="DCM threshold (55%)") %>%
      add_segments(x=min(df$Age_yr), xend=max(df$Age_yr), y=40, yend=40,
                   line=list(color='red', dash='dot'), name="Severe DCM (40%)") %>%
      layout(xaxis=list(title="Age (years)", color="white"),
             yaxis=list(title="LVEF (%)", range=c(15,70), color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$plot_ck <- renderPlotly({
    df <- clin_sim()
    plot_ly(df, x=~Age_yr, y=~CK_kUL, type='scatter', mode='lines',
            line=list(color='#e67e22', width=2)) %>%
      layout(xaxis=list(title="Age (years)", color="white"),
             yaxis=list(title="Serum CK (kU/L)", color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"))
  })

  output$plot_inflam <- renderPlotly({
    df <- clin_sim()
    plot_ly(df, x=~Age_yr, y=~Inflam_idx, type='scatter', mode='lines',
            line=list(color='#f39c12', width=2)) %>%
      layout(xaxis=list(title="Age (years)", color="white"),
             yaxis=list(title="Inflammation Index (0-100)", color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"))
  })

  ## ─────────────────────────────────────
  ## TAB 6: SCENARIO COMPARISON
  ## ─────────────────────────────────────

  comp_data <- eventReactive(input$run_comparison, {
    sim_h <- input$comp_years * 365.25 * 24
    bw    <- input$BW
    result_list <- list()

    sc_defs <- list(
      S1 = list(name="1. Natural History", ete=0,     dfz=0,          gt=FALSE, vamo=FALSE),
      S2 = list(name="2. Deflazacort",     ete=0,     dfz=0.9*bw,     gt=FALSE, vamo=FALSE),
      S3 = list(name="3. Prednisone",      ete=0,     dfz=0.75*bw,    gt=FALSE, vamo=FALSE),
      S4 = list(name="4. Ete+DFZ (ex51)",  ete=30*bw, dfz=0.9*bw,    gt=FALSE, vamo=FALSE),
      S5 = list(name="5. Casimersen+DFZ",  ete=30*bw, dfz=0.9*bw,    gt=FALSE, vamo=FALSE),
      S6 = list(name="6. Gene Therapy",    ete=0,     dfz=0,          gt=TRUE,  vamo=FALSE),
      S7 = list(name="7. Vamorolone",      ete=0,     dfz=6*bw,       gt=FALSE, vamo=TRUE)
    )

    for (sc_id in input$scenarios_sel) {
      sc <- sc_defs[[sc_id]]
      p_extra <- list()
      if (sc$vamo) { p_extra$EC50_GC <- 0.035; p_extra$kFib_stim <- 0.0022 }
      if (sc$gt)   { p_extra$GT_contribution <- 0.50 }
      m <- make_mod(p_extra)
      evs <- ev_build(dose_ete_mg=sc$ete, dose_dfz_mg=sc$dfz, sim_h=sim_h)
      df <- run_sim(m, evs, sim_h, delta=168)
      if (!is.null(df)) {
        df$Scenario <- sc$name
        result_list[[sc_id]] <- df
      }
    }
    bind_rows(result_list)
  }, ignoreNULL=FALSE)

  output$comp_6mwd <- renderPlotly({
    df <- comp_data(); req(nrow(df)>0)
    plot_ly(df, x=~Age_yr, y=~SixMWD, color=~Scenario, type='scatter', mode='lines') %>%
      layout(xaxis=list(title="Age (years)", color="white"),
             yaxis=list(title="6MWD (m)", color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"), legend=list(font=list(color="white",size=9)))
  })

  output$comp_nsaa <- renderPlotly({
    df <- comp_data(); req(nrow(df)>0)
    plot_ly(df, x=~Age_yr, y=~NSAA, color=~Scenario, type='scatter', mode='lines') %>%
      layout(xaxis=list(title="Age (years)", color="white"),
             yaxis=list(title="NSAA (0-34)", color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"), legend=list(font=list(color="white",size=9)))
  })

  output$comp_fibrosis <- renderPlotly({
    df <- comp_data(); req(nrow(df)>0)
    plot_ly(df, x=~Age_yr, y=~Fibrosis_pct, color=~Scenario, type='scatter', mode='lines') %>%
      layout(xaxis=list(title="Age (years)", color="white"),
             yaxis=list(title="Fibrosis (%)", color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"), legend=list(font=list(color="white",size=9)))
  })

  output$comp_fvc <- renderPlotly({
    df <- comp_data(); req(nrow(df)>0)
    plot_ly(df, x=~Age_yr, y=~FVC_pct, color=~Scenario, type='scatter', mode='lines') %>%
      layout(xaxis=list(title="Age (years)", color="white"),
             yaxis=list(title="FVC% Predicted", color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"), legend=list(font=list(color="white",size=9)))
  })

  output$comp_table <- renderDT({
    df <- comp_data(); req(nrow(df)>0)
    tgt_yr <- input$Age_st + 5
    tbl <- df %>%
      filter(abs(Age_yr - tgt_yr) < 0.2) %>%
      group_by(Scenario) %>%
      slice(1) %>%
      ungroup() %>%
      select(Scenario, SixMWD, NSAA, Dyst_pct, Fibrosis_pct, FVC_pct, LVEF, CK_kUL) %>%
      mutate(across(where(is.numeric), ~round(.x,1))) %>%
      rename(`6MWD (m)`=SixMWD, `NSAA`=NSAA, `Dyst (%)`=Dyst_pct,
             `Fibrosis (%)`=Fibrosis_pct, `FVC% pred`=FVC_pct,
             `LVEF (%)`=LVEF, `CK (kU/L)`=CK_kUL)
    datatable(tbl, options=list(dom='t', pageLength=10),
              rownames=FALSE,
              class="table-dark table-sm") %>%
      formatStyle(columns=1:ncol(tbl),
                  backgroundColor="#1a1a2e", color="white")
  })

  ## ─────────────────────────────────────
  ## TAB 7: BIOMARKERS
  ## ─────────────────────────────────────

  bio_sim <- reactive({
    sim_h <- input$bio_years * 365.25 * 24
    bw    <- input$BW
    treat <- input$bio_treatment
    dose_ete <- ifelse(treat=="ete_dfz", 30*bw, 0)
    dose_dfz <- ifelse(treat %in% c("dfz","ete_dfz"), 0.9*bw, 0)
    gt_c     <- ifelse(treat=="gt", 0.50, 0.0)
    m <- make_mod(list(GT_contribution=gt_c))
    events <- ev_build(dose_ete_mg=dose_ete, dose_dfz_mg=dose_dfz, sim_h=sim_h)
    df <- run_sim(m, events, sim_h, delta=168); req(!is.null(df)); df
  })

  output$bio_ck <- renderPlotly({
    df <- bio_sim()
    plot_ly(df, x=~Age_yr, y=~CK_serum, type='scatter', mode='lines',
            line=list(color='#e67e22', width=2)) %>%
      add_segments(x=min(df$Age_yr), xend=max(df$Age_yr), y=10000, yend=10000,
                   line=list(color='yellow', dash='dash'), name="10,000 U/L") %>%
      layout(xaxis=list(title="Age (years)", color="white"),
             yaxis=list(title="Serum CK (U/L)", type="log", color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$bio_dyst <- renderPlotly({
    df <- bio_sim()
    plot_ly(df, x=~Age_yr, y=~Dyst_pct, type='scatter', mode='lines',
            line=list(color='#2ecc71', width=2)) %>%
      add_segments(x=min(df$Age_yr), xend=max(df$Age_yr), y=3, yend=3,
                   line=list(color='yellow', dash='dash'), name=">3% functional") %>%
      add_segments(x=min(df$Age_yr), xend=max(df$Age_yr), y=15, yend=15,
                   line=list(color='lightgreen', dash='dot'), name="15% BMD-like") %>%
      layout(xaxis=list(title="Age (years)", color="white"),
             yaxis=list(title="Dystrophin (% normal)", color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$bio_fibers <- renderPlotly({
    df <- bio_sim()
    # stacked area chart
    plot_ly(df, x=~Age_yr) %>%
      add_trace(y=~Fiber_H/(Fiber_H+Fiber_N+Fiber_R+0.001)*100,
                name="Healthy", type='scatter', mode='lines',
                fill='tozeroy', line=list(color='#2ecc71')) %>%
      add_trace(y=~Fiber_R/(Fiber_H+Fiber_N+Fiber_R+0.001)*100,
                name="Regenerating", type='scatter', mode='lines',
                fill='tozeroy', line=list(color='#f39c12')) %>%
      add_trace(y=~Fiber_N/(Fiber_H+Fiber_N+Fiber_R+0.001)*100,
                name="Necrotic", type='scatter', mode='lines',
                fill='tozeroy', line=list(color='#e74c3c')) %>%
      layout(xaxis=list(title="Age (years)", color="white"),
             yaxis=list(title="Fiber composition (%)", range=c(0,100), color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$bio_sc <- renderPlotly({
    df <- bio_sim()
    plot_ly(df, x=~Age_yr, y=~SC_pct, type='scatter', mode='lines',
            fill='tozeroy', fillcolor='rgba(46,204,113,0.2)',
            line=list(color='#27ae60', width=2)) %>%
      add_segments(x=min(df$Age_yr), xend=max(df$Age_yr), y=30, yend=30,
                   line=list(color='orange', dash='dash'), name="Critical reserve (30%)") %>%
      layout(xaxis=list(title="Age (years)", color="white"),
             yaxis=list(title="Satellite Cell Reserve (%)", range=c(0,105), color="white"),
             paper_bgcolor="#2d2d4a", plot_bgcolor="#2d2d4a",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$bio_table <- renderDT({
    df <- bio_sim()
    checkpoints <- c(input$Age_st+2, input$Age_st+5, input$Age_st+8,
                     input$Age_st+12)
    tbl <- lapply(checkpoints, function(yr) {
      row <- df[which.min(abs(df$Age_yr - yr)),]
      data.frame(
        `Age (yr)` = round(yr,0),
        `6MWD (m)` = round(row$SixMWD,0),
        `NSAA`     = round(row$NSAA,1),
        `Dyst (%)`     = round(row$Dyst_pct,2),
        `Fibrosis (%)` = round(row$Fibrosis_pct,1),
        `CK (kU/L)`    = round(row$CK_kUL,1),
        `FVC% pred`    = round(row$FVC_pct,1),
        `LVEF (%)`     = round(row$LVEF,1),
        `SC Pool (%)`  = round(row$SC_pct,1),
        check.names=FALSE
      )
    })
    tbl_df <- bind_rows(tbl)
    datatable(tbl_df, options=list(dom='t', pageLength=10),
              rownames=FALSE, class="table-dark table-sm") %>%
      formatStyle(columns=1:ncol(tbl_df),
                  backgroundColor="#1a1a2e", color="white")
  })

}  # end server

## ─────────────────────────────────────────────────────────────
## RUN
## ─────────────────────────────────────────────────────────────
shinyApp(ui=ui, server=server)
