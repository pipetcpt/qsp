# =============================================================================
#  arvc_shiny_app.R
#  Interactive dashboard for the ARVC / arrhythmogenic cardiomyopathy QSP model
# =============================================================================
#
#  Run:
#      # in the same directory as arvc_mrgsolve_model.R
#      shiny::runApp("arvc_shiny_app.R")
#
#  The app is deliberately organised around the model's THREE COMMITMENTS
#  rather than around organ systems, because the point of the model is that
#  those three claims are what make it different from a linear cascade — and a
#  dashboard that hides them behind a tidy row of organ tabs would defeat the
#  purpose.
#
#    Commitment 1  the clock is cumulative mechanical work, not time
#                  -> tabs "The clock", "Falsifier"
#    Commitment 2  RV selectivity is Laplace, not biology
#                  -> tab "RV vs LV"
#    Commitment 3  two arrhythmia generators; a drug cannot outperform the
#                  generator it occupies
#                  -> tabs "Two generators", "Therapy bench"
#
#  TABS
#    1.  Patient          genotype, sex, body size, and the exercise history
#                         that IS the disease clock
#    2.  The clock        wall stress, fatigue dose rate, cumulative dose —
#                         the arithmetic before any ODE is solved
#    3.  Drug exposure    PK of every agent, with occupancy at its own target
#    4.  Structure        myocyte loss, fibrous and fatty replacement, wall
#                         thickness, chamber volumes, RVEF/LVEF
#    5.  Two generators   Generator I (trigger) and Generator II (re-entry)
#                         side by side on the same time axis
#    6.  ECG / imaging    TWI extent, epsilon wave, TAD, late potentials — all
#                         outputs, never inputs
#    7.  Diagnosis        the 2010 Task Force Criteria scored live, category by
#                         category, plus the risk-calculator inputs
#    8.  Outcomes         sustained VA hazard, VA-free survival, survival, ICD
#                         shock burden, amiodarone toxicity accrual
#    9.  Therapy bench    build up to four arms and compare them
#    10. RV vs LV         chamber selectivity, and what happens when KAPPA_LV
#                         is switched on
#    11. Penetrance       cohort view: penetrance by age and exercise stratum
#    12. Falsifier        PHI_EX = 0 side by side with PHI_EX = 1
#    13. Model notes      what is fitted, what is predicted, the stated miss
#
#  DISCLAIMER: educational / research tool.  Not for clinical decision-making.
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

source("arvc_mrgsolve_model.R", local = TRUE)

YR <- 365.25
AGE_MIN <- 12

# -----------------------------------------------------------------------------
#  helpers
# -----------------------------------------------------------------------------
a2d <- function(age) (age - AGE_MIN) * YR

theme_arvc <- function() {
  theme_bw(base_size = 12) +
    theme(strip.background = element_rect(fill = "#f0f0f0"),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")
}

# Build the event table for one arm from the UI inputs.
build_events <- function(inp, age_end) {
  d <- NULL
  add <- function(x) if (is.null(d)) x else c(d, x)
  if (isTRUE(inp$bb))   d <- add(chronic(CMT_BB, inp$bb_dose, 1,   inp$drug_start, age_end))
  if (isTRUE(inp$fl))   d <- add(chronic(CMT_FL, inp$fl_dose, 0.5, inp$drug_start, age_end))
  if (isTRUE(inp$so))   d <- add(chronic(CMT_SO, inp$so_dose, 0.5, inp$drug_start, age_end))
  if (isTRUE(inp$am))   d <- add(chronic(CMT_AM, inp$am_dose, 1,   inp$drug_start, age_end))
  if (isTRUE(inp$aav))  d <- add(ev(time = a2d(inp$aav_age), amt = inp$aav_td,
                                    cmt = CMT_VEC))
  if (!identical(inp$ablation, "none")) {
    frac <- if (identical(inp$ablation, "endo only")) 0.35 else 0.78
    d <- add(ev(time = a2d(inp$abl_age), amt = frac,
                cmt = which(cmt(mod) == "ABL_HOM")))
  }
  if (is.null(d)) d <- ev(time = 0, amt = 0, cmt = CMT_BB)
  d
}

build_params <- function(inp) {
  p <- GENOTYPE[[inp$genotype]]
  p$EX         <- inp$ex_dose
  p$MALE       <- ifelse(inp$sex == "male", 1, 0)
  p$SEX_K      <- ifelse(inp$sex == "male", 1, 0.72)
  p$BSA        <- inp$bsa
  p$FRAILTY    <- inp$frailty
  p$PHI_EX     <- ifelse(isTRUE(inp$falsify), 0, 1)
  p$ON_ICD     <- ifelse(isTRUE(inp$icd), 1, 0)
  p$ON_LOADRED <- ifelse(isTRUE(inp$loadred), 1, 0)
  p$ON_MRA     <- ifelse(isTRUE(inp$mra), 1, 0)
  p$ON_IL1     <- ifelse(isTRUE(inp$il1), 1, 0)
  p$ON_GC      <- ifelse(isTRUE(inp$gc), 1, 0)
  p$ON_GSK     <- ifelse(isTRUE(inp$gsk), 1, 0)
  if (!is.null(inp$kappa_lv_override) && inp$kappa_lv_override > 0)
    p$KAPPA_LV <- inp$kappa_lv_override
  # Exercise history: the SAME patient with two load histories is the
  # comparison the model exists to make, so the dose change is a parameter
  # (T_RESTR) inside the model rather than a spliced two-part simulation.
  if (isTRUE(inp$restrict)) {
    p$EX2 <- inp$restrict_to
    p$T_RESTR <- a2d(inp$restrict_age)
  } else {
    p$EX2 <- inp$ex_dose
    p$T_RESTR <- 1e9
  }
  p
}

run_one <- function(inp, age_end, delta = 30) {
  mod %>%
    param(build_params(inp)) %>%
    mrgsim(events = build_events(inp, age_end), end = a2d(age_end),
           delta = delta) %>%
    as_tibble()
}

TFC_LABEL <- c("none", "possible", "borderline", "definite")

# =============================================================================
#  UI
# =============================================================================
ui <- fluidPage(
  titlePanel(
    div(
      h3("Arrhythmogenic Right Ventricular Cardiomyopathy — QSP dashboard"),
      p(em(paste("The clock is cumulative mechanical work, not time ·",
                 "RV selectivity is Laplace, not biology ·",
                 "two arrhythmia generators on different clocks"))),
      p(strong("Educational / research model. Not for clinical use."),
        style = "color:#b71c1c")
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Patient"),
      selectInput("genotype", "Genotype",
                  choices = c("WT (gene-elusive)" = "WT",
                              "PKP2 truncating" = "PKP2tv",
                              "DSP truncating" = "DSPtv",
                              "FLNC truncating" = "FLNCtv"),
                  selected = "PKP2tv"),
      radioButtons("sex", "Sex", c("male", "female"), inline = TRUE),
      sliderInput("bsa", "Body surface area (m2)", 1.4, 2.4, 1.9, 0.05),
      sliderInput("frailty", "Between-subject frailty on the fatigue rate",
                  0.4, 2.5, 1.0, 0.05),
      sliderInput("age_end", "Simulate to age (y)", 30, 85, 60, 1),

      hr(),
      h4("The clock — exercise history"),
      helpText("This is not a lifestyle field. It is the model's clock."),
      sliderInput("ex_dose", "Exercise dose (MET-hours/week)", 3, 120, 60, 1),
      helpText("6 sedentary · 15 guideline · 60 competitive · 100 elite"),
      checkboxInput("restrict", "Restrict exercise at a given age", FALSE),
      conditionalPanel(
        "input.restrict",
        sliderInput("restrict_age", "Restrict from age (y)", 13, 60, 20, 1),
        sliderInput("restrict_to", "New dose (MET-h/wk)", 3, 60, 15, 1)
      ),

      hr(),
      h4("Antiarrhythmic / device"),
      sliderInput("drug_start", "Start drugs at age (y)", 13, 70, 38, 1),
      checkboxInput("bb", "Beta-blocker (nadolol) — Generator I", TRUE),
      conditionalPanel("input.bb",
                       sliderInput("bb_dose", "nadolol mg daily", 20, 240, 80, 20)),
      checkboxInput("fl", "Flecainide — Generator I via RyR2, Generator II via Nav1.5",
                    FALSE),
      conditionalPanel("input.fl",
                       sliderInput("fl_dose", "flecainide mg per dose (bd)",
                                   50, 150, 100, 25)),
      checkboxInput("so", "Sotalol — IKr, plus a proarrhythmic arm", FALSE),
      conditionalPanel("input.so",
                       sliderInput("so_dose", "sotalol mg per dose (bd)",
                                   80, 240, 160, 40)),
      checkboxInput("am", "Amiodarone — both generators, and a 58-day half-life",
                    FALSE),
      conditionalPanel("input.am",
                       sliderInput("am_dose", "amiodarone mg daily",
                                   100, 400, 200, 100)),
      selectInput("ablation", "Catheter ablation",
                  c("none", "endo only", "endo + epicardial")),
      conditionalPanel("input.ablation != 'none'",
                       sliderInput("abl_age", "Ablation at age (y)", 15, 70, 38, 1)),
      checkboxInput("icd", "ICD implanted", FALSE),

      hr(),
      h4("Upstream / disease-modifying"),
      checkboxInput("loadred",
                    "Load-reducing therapy (diuretic + nitrate / ARNI) — the CLOCK",
                    FALSE),
      checkboxInput("mra", "MRA — fibrogenesis", FALSE),
      checkboxInput("il1", "IL-1 blockade — the amplifier", FALSE),
      checkboxInput("gc", "Corticosteroid — the amplifier", FALSE),
      checkboxInput("gsk", "GSK-3beta inhibition — the adipogenic switch", FALSE),
      checkboxInput("aav", "AAV9-PKP2 gene therapy (single infusion)", FALSE),
      conditionalPanel(
        "input.aav",
        sliderInput("aav_age", "Infusion at age (y)", 13, 70, 26, 1),
        sliderInput("aav_td", "Transduced cardiomyocyte fraction", 0.1, 0.8, 0.62, 0.02)
      ),

      hr(),
      h4("Falsifier"),
      checkboxInput("falsify",
                    "PHI_EX = 0: make damage load-INDEPENDENT", FALSE),
      helpText(paste("With this on, the model becomes a calendar-clock model.",
                     "Exercise should stop mattering, penetrance should become",
                     "complete and age-fixed, and the RV should lose its",
                     "priority. If it does not, commitment 1 is decorative.")),
      sliderInput("kappa_lv_override",
                  "Override KAPPA_LV (0 = use genotype default)", 0, 5, 0, 0.1)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        # ------------------------------------------------------------- 1
        tabPanel(
          "1. Patient",
          br(),
          fluidRow(
            column(6, h4("Where this patient is now"),
                   tableOutput("summary_now")),
            column(6, h4("What the genotype does — and only this"),
                   verbatimTextOutput("geno_note"))
          ),
          hr(),
          h4("Lifetime trajectory of the four things a clinic measures"),
          plotOutput("p_overview", height = "460px")
        ),

        # ------------------------------------------------------------- 2
        tabPanel(
          "2. The clock",
          br(),
          h4("The fatigue arithmetic, before any ODE is solved"),
          p(paste("Fatigue dose rate = beats x (wall stress)^NMECH, averaged",
                  "over the week and normalised so a perfectly rested healthy",
                  "chamber scores 1.00. NMECH is fixed at 4 a priori; the two",
                  "exercise stress coefficients come from exercise CMR",
                  "(RV +125%, LV +14%).")),
          tableOutput("load_table"),
          hr(),
          h4("Wall stress and fatigue dose over a lifetime"),
          plotOutput("p_clock", height = "420px"),
          hr(),
          h4("Cumulative mechanical dose — the model's only clock"),
          plotOutput("p_dose", height = "300px")
        ),

        # ------------------------------------------------------------- 3
        tabPanel(
          "3. Drug exposure",
          br(),
          h4("Plasma concentrations"),
          plotOutput("p_pk", height = "380px"),
          hr(),
          h4("Occupancy at the drug's own target"),
          p(paste("Each drug is attached to a specific node with occupancy",
                  "computed from published Ki/IC50, unbound fraction and PK.",
                  "A drug cannot outperform the node it occupies.")),
          tableOutput("occ_table")
        ),

        # ------------------------------------------------------------- 4
        tabPanel(
          "4. Structure",
          br(),
          h4("Myocardial composition of the RV free wall"),
          plotOutput("p_composition", height = "340px"),
          hr(),
          h4("Chamber size and function"),
          plotOutput("p_chamber", height = "340px"),
          hr(),
          helpText(paste("Replacement fills the space the myocytes vacated —",
                         "it does not create it. Fibrous tissue bears load",
                         "(0.70) and fat essentially does not (0.35 bulk",
                         "only), which is why the wall thins even as it is",
                         "replaced, and why the process is autocatalytic."))
        ),

        # ------------------------------------------------------------- 5
        tabPanel(
          "5. Two generators",
          br(),
          h4("Generator I and Generator II on the same time axis"),
          plotOutput("p_generators", height = "380px"),
          hr(),
          fluidRow(
            column(6, h4("Generator I — trigger"),
                   plotOutput("p_gen1_detail", height = "300px")),
            column(6, h4("Generator II — re-entry"),
                   plotOutput("p_gen2_detail", height = "300px"))
          ),
          hr(),
          helpText(paste("Generator I exists in a structurally normal heart:",
                         "PKP2 loss strips Nav1.5 and Cx43 from the",
                         "intercalated disc and destabilises RyR2. Generator",
                         "II needs scar, and peaks at about 50/50",
                         "interdigitation — a homogeneous scar has no",
                         "conducting channels either."))
        ),

        # ------------------------------------------------------------- 6
        tabPanel(
          "6. ECG / imaging",
          br(),
          h4("Every item here is an output of the state variables"),
          plotOutput("p_ecg", height = "420px"),
          hr(),
          tableOutput("ecg_table")
        ),

        # ------------------------------------------------------------- 7
        tabPanel(
          "7. Diagnosis",
          br(),
          h4("2010 Task Force Criteria, scored live"),
          plotOutput("p_tfc", height = "300px"),
          hr(),
          fluidRow(
            column(6, h4("Category counts at the selected age"),
                   tableOutput("tfc_table")),
            column(6, h4("Risk-calculator inputs, emitted not scored"),
                   p(paste("The model deliberately does not reproduce a",
                           "published linear predictor; it emits the seven",
                           "inputs so an external, independently validated",
                           "score can be applied to its trajectories.")),
                   tableOutput("risk_inputs"))
          ),
          hr(),
          sliderInput("dx_age", "Read the criteria at age (y)", 13, 85, 40, 1)
        ),

        # ------------------------------------------------------------- 8
        tabPanel(
          "8. Outcomes",
          br(),
          h4("Sustained ventricular arrhythmia hazard and survival"),
          plotOutput("p_outcomes", height = "420px"),
          hr(),
          fluidRow(
            column(6, h4("Device burden"),
                   plotOutput("p_shocks", height = "280px")),
            column(6, h4("Amiodarone: the cost of the most effective drug"),
                   plotOutput("p_tox", height = "280px"))
          )
        ),

        # ------------------------------------------------------------- 9
        tabPanel(
          "9. Therapy bench",
          br(),
          h4("Compare the standard arms on the current patient"),
          p(paste("Each arm starts at the 'Start drugs at age' set in the",
                  "sidebar. Watch which arms move STRUCTURE and which move",
                  "only the ARRHYTHMIA — that separation is commitment 3.")),
          actionButton("run_bench", "Run bench", class = "btn-primary"),
          br(), br(),
          plotOutput("p_bench", height = "460px"),
          hr(),
          tableOutput("bench_table")
        ),

        # ------------------------------------------------------------- 10
        tabPanel(
          "10. RV vs LV",
          br(),
          h4("Chamber selectivity with no chamber-specific biology"),
          plotOutput("p_rvlv", height = "400px"),
          hr(),
          h4("RV:LV ratio against exercise dose"),
          p(paste("The gap should WIDEN with training. That is a prediction of",
                  "the +125% vs +14% wall-stress asymmetry raised to the",
                  "fourth power, not an input.")),
          plotOutput("p_rvlv_ratio", height = "320px"),
          hr(),
          helpText(paste("KAPPA_LV is the only chamber-specific term in the",
                         "model and it is switched on only for DSP and FLNC.",
                         "Use the sidebar override to see what it takes to",
                         "make a PKP2 patient left-dominant."))
        ),

        # ------------------------------------------------------------- 11
        tabPanel(
          "11. Penetrance",
          br(),
          h4("Penetrance is a property of a cohort, not of a trajectory"),
          p(paste("15 subjects per stratum with a log-normal spread on the",
                  "fatigue rate (SD 0.55, fixed quantiles). Nothing here is",
                  "fitted.")),
          actionButton("run_pen", "Run cohort", class = "btn-primary"),
          br(), br(),
          plotOutput("p_penetrance", height = "420px"),
          hr(),
          tableOutput("pen_table"),
          helpText(paste("Observed anchors: penetrance in genotype-positive",
                         "relatives is incomplete and roughly 35-50% by",
                         "mid-life; most carriers in an unselected biobank",
                         "have no phenotype; competitive athletes are the",
                         "group in which penetrance approaches complete."))
        ),

        # ------------------------------------------------------------- 12
        tabPanel(
          "12. Falsifier",
          br(),
          h4("PHI_EX = 1 against PHI_EX = 0, everything else identical"),
          actionButton("run_fals", "Run both", class = "btn-danger"),
          br(), br(),
          plotOutput("p_falsifier", height = "460px"),
          hr(),
          tableOutput("fals_table"),
          helpText(paste("With the load term removed the model predicts that",
                         "exercise is irrelevant, penetrance is complete and",
                         "age-fixed, gene-elusive ARVC is impossible, the",
                         "ventricles are affected equally, and load-reducing",
                         "therapy does nothing. All five are contradicted by",
                         "data the model was not fitted to."))
        ),

        # ------------------------------------------------------------- 13
        tabPanel(
          "13. Model notes",
          br(),
          h4("What is fitted"),
          tags$ul(
            tags$li(strong("K_INJ"), " — fatigue-failure scale, to the median age at definite Task Force diagnosis in PKP2 carriers on ordinary recreational activity"),
            tags$li(strong("H0_VA"), " — arrhythmia hazard scale, to ~10%/yr sustained VA in definite disease"),
            tags$li(strong("LAM2"), " — Generator II weight, to the amiodarone-versus-sotalol gap"),
            tags$li(strong("K_DIL"), " — dilatation gain, to the RVEDVi trajectory into overt disease")
          ),
          h4("What is predicted"),
          tags$ul(
            tags$li("the ~3x exercise hazard ratio and its dose-response"),
            tags$li("incomplete penetrance, and its dependence on exercise stratum"),
            tags$li("gene-elusive ARVC in extreme-dose athletes only"),
            tags$li("male predominance, from a single multiplier"),
            tags$li("RV before LV with no chamber-specific biology, and the RV/LV gap widening with training"),
            tags$li("load-reducing therapy preventing the phenotype"),
            tags$li("the sotalol null — its IKr and proarrhythmic arms cancel"),
            tags$li("flecainide helping early and hurting late: one drug, two signs"),
            tags$li("the endocardial-only versus endo-epicardial ablation gap"),
            tags$li("the ICD's outcome-only benefit: all of the mortality, none of the disease"),
            tags$li("AAV-PKP2 gene therapy in which timing dominates dose")
          ),
          h4("A stated miss"),
          p(paste("This model gives beta-blocker monotherapy a larger",
                  "arrhythmia reduction than the North American ARVC registry",
                  "observed. Almost all of Generator I is routed through",
                  "beta-adrenergic drive here; if the registry is right, part",
                  "of the RyR2 leak must be adrenergic-independent. That is",
                  "the cleanest place to try to falsify commitment 3.")),
          h4("A second miss, which is the same miss"),
          p(paste("Flecainide added to a beta-blocker comes out slightly",
                  "harmful here at every age, whereas the published series",
                  "found it reduced VT. The reason is that beta-blockade has",
                  "already driven Generator I to zero in this model, leaving",
                  "only flecainide's conduction cost. Both misses therefore",
                  "say the same thing: some of the RyR2 leak has to be",
                  "adrenergic-independent. One added term repairs both, which",
                  "makes it a testable prediction about the model rather than",
                  "an excuse for it.")),
          h4("Verification"),
          p("The dependency-free Python twin regenerates every number:"),
          tags$pre("python3 arvc_reference_impl.py --check"),
          h4("Files"),
          tags$ul(
            tags$li(code("arvc_qsp_model.dot"), " — mechanistic map (20 clusters)"),
            tags$li(code("arvc_mrgsolve_model.R"), " — 54-compartment ODE model"),
            tags$li(code("arvc_reference_impl.py"), " — numerical twin + self-checks"),
            tags$li(code("arvc_references.md"), " — every citation, PMID-resolved")
          )
        )
      )
    )
  )
)

# =============================================================================
#  SERVER
# =============================================================================
server <- function(input, output, session) {

  sim <- reactive({
    run_one(input, input$age_end, delta = 30)
  })

  at_age <- function(df, age) {
    df %>% filter(AGE <= age) %>% slice_tail(n = 1)
  }

  # ------------------------------------------------------------------ tab 1
  output$summary_now <- renderTable({
    o <- at_age(sim(), input$age_end)
    tibble(
      quantity = c("age (y)", "Task Force call", "RVEF", "RVEDVi (mL/m2)",
                   "LVEF", "residual myocytes (%)", "fibrofatty RV (fraction)",
                   "PVC / 24 h", "RV conduction velocity (cm/s)",
                   "NT-proBNP (pg/mL)", "sustained VA hazard (%/yr)",
                   "VA-free survival", "survival"),
      value = c(sprintf("%.1f", o$AGE),
                TFC_LABEL[o$TFC + 1],
                sprintf("%.2f", o$RVEF),
                sprintf("%.0f", o$RVEDVI),
                sprintf("%.2f", o$LVEF),
                sprintf("%.0f", o$RESIDMYO),
                sprintf("%.3f", o$FIBFATRV),
                sprintf("%.0f", o$PVC24),
                sprintf("%.0f", o$CV_RV),
                sprintf("%.0f", o$NTBNP),
                sprintf("%.1f", 100 * o$HVA_YR),
                sprintf("%.2f", o$VAFREE),
                sprintf("%.2f", o$SURV))
    )
  })

  output$geno_note <- renderText({
    g <- GENOTYPE[[input$genotype]]
    paste0(
      "genotype: ", input$genotype, "\n",
      "PKP2 at the intercalated disc : ", g$PKP2_SET, "\n",
      "desmoplakin at the disc       : ", g$DSP_SET, "\n",
      "KAPPA_LV (LV-specific term)   : ", g$KAPPA_LV,
      if (g$KAPPA_LV > 1) "   <- ON: this genotype is left-dominant"
      else "   <- OFF: RV priority must come from Laplace alone", "\n\n",
      "The genotype enters the model at exactly ONE place: it divides the\n",
      "fatigue rate by (desmosomal reserve)^GAMMA_G.  It does not set the\n",
      "amount of damage, the chamber, the arrhythmia or the age of onset.\n",
      "Those all come out of the load history."
    )
  })

  output$p_overview <- renderPlot({
    sim() %>%
      select(AGE, RVEF, RVEDVI, FIBFATRV, PVC24) %>%
      pivot_longer(-AGE) %>%
      ggplot(aes(AGE, value)) +
      geom_line(linewidth = 0.9, colour = "#6a1b9a") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "age (years)", y = NULL) +
      theme_arvc()
  })

  # ------------------------------------------------------------------ tab 2
  output$load_table <- renderTable({
    nmech <- 4; hr_r <- 62; hr_x <- 158; met <- 8
    lf <- function(ex, kex) {
      f <- min(0.35, ex / (met * 168))
      ((1 - f) * hr_r + f * hr_x * (1 + kex)^nmech) / hr_r
    }
    ex <- c(6, 15, 60, 100, input$ex_dose)
    tibble(
      exercise = c("sedentary", "guideline", "competitive", "elite",
                   "THIS PATIENT"),
      `MET-h/wk` = ex,
      LOAD_RV = sprintf("%.3f", sapply(ex, lf, kex = 1.25)),
      LOAD_LV = sprintf("%.3f", sapply(ex, lf, kex = 0.14)),
      `RV/LV` = sprintf("%.2f", sapply(ex, lf, kex = 1.25) /
                          sapply(ex, lf, kex = 0.14))
    )
  })

  output$p_clock <- renderPlot({
    sim() %>%
      select(AGE, SIGRV, SIGLV, LOADRVO, LOADLVO) %>%
      pivot_longer(-AGE) %>%
      mutate(panel = ifelse(grepl("^SIG", name), "wall stress (relative)",
                            "fatigue dose rate"),
             chamber = ifelse(grepl("RV", name), "RV", "LV")) %>%
      ggplot(aes(AGE, value, colour = chamber)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~panel, scales = "free_y") +
      scale_colour_manual(values = c(RV = "#8e24aa", LV = "#b39ddb")) +
      labs(x = "age (years)", y = NULL, colour = NULL) +
      theme_arvc()
  })

  output$p_dose <- renderPlot({
    sim() %>% select(AGE, D_RV, D_LV) %>% pivot_longer(-AGE) %>%
      ggplot(aes(AGE, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = c(D_RV = "#8e24aa", D_LV = "#b39ddb")) +
      labs(x = "age (years)", y = "cumulative dose (load-years)",
           colour = NULL) +
      theme_arvc()
  })

  # ------------------------------------------------------------------ tab 3
  output$p_pk <- renderPlot({
    s <- sim()
    s %>% transmute(AGE,
                    nadolol   = A_BB_C / 140,
                    flecainide = A_FL_C / 600,
                    sotalol   = A_SO_C / 100,
                    amiodarone = A_AM_C / 60) %>%
      pivot_longer(-AGE) %>% filter(value > 1e-12) %>%
      ggplot(aes(AGE, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "age (years)", y = "plasma concentration (mg/L)",
           colour = NULL) +
      theme_arvc() + theme(legend.position = "none")
  })

  output$occ_table <- renderTable({
    o <- at_age(sim(), input$age_end)
    c_bb <- 0.70 * (o$A_BB_C / 140) / 309.4 * 1e6
    c_fl <- 0.60 * (o$A_FL_C / 600) / 414.3 * 1e6
    c_so <- (o$A_SO_C / 100) / 272.4 * 1e6
    c_am <- (o$A_AM_C / 60) / 645.3 * 1e6
    tibble(
      drug = c("nadolol", "flecainide", "flecainide", "sotalol", "sotalol",
               "amiodarone", "amiodarone"),
      target = c("beta1 (Kd 4 nM)", "RyR2 open state (IC50 2.5 uM)",
                 "Nav1.5 (IC50 6 uM)", "beta1 (Kd 320 nM)",
                 "IKr (IC50 8 uM)", "Generator I (IC50 1.2 uM)",
                 "Generator II (IC50 1.6 uM)"),
      `free conc (nM)` = sprintf("%.0f", c(c_bb, c_fl, c_fl, c_so, c_so,
                                           c_am, c_am)),
      occupancy = sprintf("%.0f%%", 100 * c(
        c_bb / (c_bb + 4), c_fl / (c_fl + 2500), c_fl / (c_fl + 6000),
        c_so / (c_so + 320), c_so / (c_so + 8000),
        c_am / (c_am + 1200), c_am / (c_am + 1600))),
      `acts on` = c("Generator I", "Generator I (helps)",
                    "Generator II (HURTS: slows conduction)",
                    "Generator I", "Generator II (helps) + proarrhythmia",
                    "Generator I", "Generator II")
    )
  })

  # ------------------------------------------------------------------ tab 4
  output$p_composition <- renderPlot({
    sim() %>% select(AGE, myocyte = MYO_RV, fibrous = FIB_RV, fat = FAT_RV) %>%
      pivot_longer(-AGE) %>%
      ggplot(aes(AGE, value, fill = name)) +
      geom_area(alpha = 0.85) +
      scale_fill_manual(values = c(myocyte = "#ef5350", fibrous = "#a1887f",
                                   fat = "#ffd54f")) +
      labs(x = "age (years)", y = "fraction of the RV free wall", fill = NULL) +
      theme_arvc()
  })

  output$p_chamber <- renderPlot({
    sim() %>% select(AGE, RVEDVI, LVEDVI, RVEF, LVEF) %>%
      pivot_longer(-AGE) %>%
      ggplot(aes(AGE, value)) +
      geom_line(linewidth = 0.9, colour = "#00695c") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "age (years)", y = NULL) + theme_arvc()
  })

  # ------------------------------------------------------------------ tab 5
  output$p_generators <- renderPlot({
    sim() %>% select(AGE, `Generator I` = G1OUT, `Generator II` = G2OUT) %>%
      pivot_longer(-AGE) %>%
      ggplot(aes(AGE, value, colour = name)) +
      geom_line(linewidth = 1.0) +
      scale_colour_manual(values = c(`Generator I` = "#1565c0",
                                     `Generator II` = "#4e342e")) +
      labs(x = "age (years)", y = "generator drive (dimensionless)",
           colour = NULL) +
      theme_arvc()
  })

  output$p_gen1_detail <- renderPlot({
    sim() %>% select(AGE, RYR_LEAK, CA_DIA, PVC24, BETA1_D) %>%
      pivot_longer(-AGE) %>%
      ggplot(aes(AGE, value)) + geom_line(colour = "#1565c0", linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "age (years)", y = NULL) + theme_arvc()
  })

  output$p_gen2_detail <- renderPlot({
    sim() %>% select(AGE, CV_RV, SCARHETRV, FIBFATRV, ABL_HOM) %>%
      pivot_longer(-AGE) %>%
      ggplot(aes(AGE, value)) + geom_line(colour = "#4e342e", linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "age (years)", y = NULL) + theme_arvc()
  })

  # ------------------------------------------------------------------ tab 6
  output$p_ecg <- renderPlot({
    sim() %>% select(AGE, `TWI leads` = NTWI, `TAD (ms)` = TAD_MS,
                     `epsilon wave` = EPSILONW,
                     `late potentials` = LATEPOT,
                     `CV RV (cm/s)` = CV_RV, `PVC / 24 h` = PVC24) %>%
      pivot_longer(-AGE) %>%
      ggplot(aes(AGE, value)) +
      geom_line(colour = "#558b2f", linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "age (years)", y = NULL) + theme_arvc()
  })

  output$ecg_table <- renderTable({
    o <- at_age(sim(), input$age_end)
    tibble(
      finding = c("T-wave inversion, right precordial leads",
                  "terminal activation duration",
                  "epsilon wave", "SAECG late potentials",
                  "PVC burden"),
      value = c(sprintf("%.0f leads", o$NTWI),
                sprintf("%.0f ms", o$TAD_MS),
                ifelse(o$EPSILONW > 0.5, "present", "absent"),
                ifelse(o$LATEPOT > 0.5, "present", "absent"),
                sprintf("%.0f / 24 h", o$PVC24)),
      `Task Force weight` = c(
        ifelse(o$NTWI >= 3, "MAJOR", ifelse(o$NTWI >= 2, "minor", "-")),
        ifelse(o$TAD_MS >= 55, "minor", "-"),
        ifelse(o$EPSILONW > 0.5, "MAJOR", "-"),
        ifelse(o$LATEPOT > 0.5, "minor", "-"),
        ifelse(o$PVC24 > 500, "minor", "-"))
    )
  })

  # ------------------------------------------------------------------ tab 7
  output$p_tfc <- renderPlot({
    sim() %>% select(AGE, major = NMAJ, minor = NMIN, call = TFC) %>%
      pivot_longer(-AGE) %>%
      ggplot(aes(AGE, value, colour = name)) +
      geom_step(linewidth = 0.9) +
      labs(x = "age (years)",
           y = "criteria count / call (0 none .. 3 definite)",
           colour = NULL) +
      theme_arvc()
  })

  output$tfc_table <- renderTable({
    o <- at_age(sim(), input$dx_age)
    tibble(
      item = c("age", "major criteria", "minor criteria", "diagnostic call"),
      value = c(sprintf("%.1f y", o$AGE), sprintf("%.0f", o$NMAJ),
                sprintf("%.0f", o$NMIN), TFC_LABEL[o$TFC + 1])
    )
  })

  output$risk_inputs <- renderTable({
    o <- at_age(sim(), input$dx_age)
    tibble(
      predictor = c("age (y)", "sex", "recent syncope",
                    "NSVT", "PVC count / 24 h",
                    "leads with T-wave inversion", "RVEF"),
      value = c(sprintf("%.0f", o$AGE), input$sex,
                "not modelled explicitly",
                ifelse(o$HVA_YR > 0.045, "yes", "no"),
                sprintf("%.0f", o$PVC24), sprintf("%.0f", o$NTWI),
                sprintf("%.2f", o$RVEF))
    )
  })

  # ------------------------------------------------------------------ tab 8
  output$p_outcomes <- renderPlot({
    sim() %>% select(AGE, `VA hazard (%/yr)` = HVA_YR,
                     `VA-free survival` = VAFREE,
                     `survival` = SURV, `NT-proBNP` = NTBNP) %>%
      mutate(`VA hazard (%/yr)` = 100 * `VA hazard (%/yr)`) %>%
      pivot_longer(-AGE) %>%
      ggplot(aes(AGE, value)) +
      geom_line(linewidth = 0.9, colour = "#1b5e20") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "age (years)", y = NULL) + theme_arvc()
  })

  output$p_shocks <- renderPlot({
    sim() %>% ggplot(aes(AGE, ICD_SHK)) +
      geom_line(linewidth = 0.9, colour = "#ad1457") +
      labs(x = "age (years)", y = "cumulative ICD shocks") + theme_arvc()
  })

  output$p_tox <- renderPlot({
    sim() %>% ggplot(aes(AGE, AMIO_TOX)) +
      geom_line(linewidth = 0.9, colour = "#e65100") +
      labs(x = "age (years)",
           y = "cumulative amiodarone toxicity index") + theme_arvc()
  })

  # ------------------------------------------------------------------ tab 9
  bench <- eventReactive(input$run_bench, {
    arms <- list(
      "no therapy"                  = list(),
      "exercise restriction"        = list(restrict = TRUE,
                                           restrict_age = input$drug_start,
                                           restrict_to = 15),
      "load-reducing therapy"       = list(loadred = TRUE),
      "beta-blocker"                = list(bb = TRUE),
      "beta-blocker + flecainide"   = list(bb = TRUE, fl = TRUE),
      "sotalol"                     = list(so = TRUE),
      "amiodarone"                  = list(am = TRUE),
      "endo-only ablation"          = list(ablation = "endo only",
                                           abl_age = input$drug_start),
      "endo + epicardial ablation"  = list(ablation = "endo + epicardial",
                                           abl_age = input$drug_start),
      "ICD alone"                   = list(icd = TRUE)
    )
    base <- reactiveValuesToList(input)
    bind_rows(lapply(names(arms), function(nm) {
      inp <- modifyList(base, list(bb = FALSE, fl = FALSE, so = FALSE,
                                   am = FALSE, icd = FALSE, loadred = FALSE,
                                   restrict = FALSE, aav = FALSE,
                                   ablation = "none"))
      inp <- modifyList(inp, arms[[nm]])
      run_one(inp, input$age_end, delta = 60) %>% mutate(arm = nm)
    }))
  })

  output$p_bench <- renderPlot({
    bench() %>%
      select(AGE, arm, `VA hazard (%/yr)` = HVA_YR, RVEF,
             `fibrofatty RV` = FIBFATRV, `survival` = SURV) %>%
      mutate(`VA hazard (%/yr)` = 100 * `VA hazard (%/yr)`) %>%
      pivot_longer(-c(AGE, arm)) %>%
      ggplot(aes(AGE, value, colour = arm)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "age (years)", y = NULL, colour = NULL) + theme_arvc()
  })

  output$bench_table <- renderTable({
    bench() %>% group_by(arm) %>% slice_tail(n = 1) %>% ungroup() %>%
      transmute(arm,
                `VA hazard %/yr` = round(100 * HVA_YR, 2),
                RVEF = round(RVEF, 2),
                `fibrofatty RV` = round(FIBFATRV, 3),
                `VA-free surv` = round(VAFREE, 3),
                survival = round(SURV, 3),
                `TFC` = TFC_LABEL[TFC + 1]) %>%
      arrange(`VA hazard %/yr`)
  })

  # ------------------------------------------------------------------ tab 10
  output$p_rvlv <- renderPlot({
    sim() %>% select(AGE, RV = FIBFATRV, LV = FIBFATLV) %>%
      pivot_longer(-AGE) %>%
      ggplot(aes(AGE, value, colour = name)) +
      geom_line(linewidth = 1.0) +
      scale_colour_manual(values = c(RV = "#8e24aa", LV = "#00838f")) +
      labs(x = "age (years)", y = "fibrofatty replacement (fraction)",
           colour = "chamber") + theme_arvc()
  })

  output$p_rvlv_ratio <- renderPlot({
    nmech <- 4; hr_r <- 62; hr_x <- 158; met <- 8
    lf <- function(ex, kex) {
      f <- pmin(0.35, ex / (met * 168))
      ((1 - f) * hr_r + f * hr_x * (1 + kex)^nmech) / hr_r
    }
    ex <- seq(3, 120, by = 1)
    tibble(ex = ex, ratio = lf(ex, 1.25) / lf(ex, 0.14)) %>%
      ggplot(aes(ex, ratio)) +
      geom_line(linewidth = 1.0, colour = "#6a1b9a") +
      geom_vline(xintercept = input$ex_dose, linetype = 2) +
      labs(x = "exercise dose (MET-hours/week)",
           y = "RV : LV fatigue dose rate ratio") + theme_arvc()
  })

  # ------------------------------------------------------------------ tab 11
  pen <- eventReactive(input$run_pen, {
    run_penetrance(genotype = input$genotype, n = 15, age_end = 75)
  })

  output$p_penetrance <- renderPlot({
    pen() %>% group_by(exercise, sex, ID) %>%
      summarise(age_dx = suppressWarnings(min(AGE[TFC >= 3])), .groups = "drop") %>%
      tidyr::crossing(age = seq(15, 75, 1)) %>%
      group_by(exercise, sex, age) %>%
      summarise(pen = mean(age_dx <= age), .groups = "drop") %>%
      ggplot(aes(age, 100 * pen, colour = exercise, linetype = sex)) +
      geom_line(linewidth = 0.9) +
      labs(x = "age (years)", y = "penetrance (% definite TFC)",
           colour = NULL, linetype = NULL) + theme_arvc()
  })

  output$pen_table <- renderTable({
    pen() %>% group_by(exercise, sex, ID) %>%
      summarise(age_dx = suppressWarnings(min(AGE[TFC >= 3])), .groups = "drop") %>%
      group_by(exercise, sex) %>%
      summarise(`by 30` = sprintf("%.0f%%", 100 * mean(age_dx <= 30)),
                `by 40` = sprintf("%.0f%%", 100 * mean(age_dx <= 40)),
                `by 50` = sprintf("%.0f%%", 100 * mean(age_dx <= 50)),
                `by 60` = sprintf("%.0f%%", 100 * mean(age_dx <= 60)),
                .groups = "drop")
  })

  # ------------------------------------------------------------------ tab 12
  fals <- eventReactive(input$run_fals, {
    base <- reactiveValuesToList(input)
    bind_rows(lapply(c(1, 0), function(phi) {
      bind_rows(lapply(c(6, 60), function(ex) {
        inp <- modifyList(base, list(falsify = (phi == 0), ex_dose = ex,
                                     restrict = FALSE))
        run_one(inp, 80, delta = 90) %>%
          mutate(PHI_EX = phi,
                 exercise = ifelse(ex == 6, "sedentary", "competitive"))
      }))
    }))
  })

  output$p_falsifier <- renderPlot({
    fals() %>%
      mutate(model = ifelse(PHI_EX == 1, "PHI_EX = 1 (fatigue clock)",
                            "PHI_EX = 0 (calendar clock)")) %>%
      select(AGE, model, exercise, RV = FIBFATRV, LV = FIBFATLV, RVEF) %>%
      pivot_longer(c(RV, LV, RVEF)) %>%
      ggplot(aes(AGE, value, colour = exercise)) +
      geom_line(linewidth = 0.9) +
      facet_grid(name ~ model, scales = "free_y") +
      labs(x = "age (years)", y = NULL, colour = NULL) + theme_arvc()
  })

  output$fals_table <- renderTable({
    fals() %>% group_by(PHI_EX, exercise) %>%
      summarise(age_definite = suppressWarnings(min(AGE[TFC >= 3])),
                `RV:LV at 50` = {
                  i <- which.min(abs(AGE - 50))
                  sprintf("%.2f", FIBFATRV[i] / max(1e-9, FIBFATLV[i]))
                }, .groups = "drop") %>%
      mutate(age_definite = ifelse(is.finite(age_definite),
                                   sprintf("%.1f", age_definite), "never"),
             model = ifelse(PHI_EX == 1, "fatigue clock", "calendar clock")) %>%
      select(model, exercise, age_definite, `RV:LV at 50`)
  })
}

shinyApp(ui, server)
