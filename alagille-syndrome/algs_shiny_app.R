# =====================================================================
# Alagille Syndrome (ALGS) QSP — Shiny dashboard
#   Author : Claude Code Routine (2026-07-30)
#
#   Front end for algs_mrgsolve_model.R (44 ODEs).  Twelve tabs, each
#   tied to one of the seven axes in the model header:
#
#     1  Patient profile      — duct capacity is the patient, not sBA
#     2  Drug PK / ASBT       — lumen-acting, <1% absorbed
#     3  Bile acid system     — the compartments and the balance
#     4  Pruritus             — two components, only one is the drug
#     5  Trial replication    — ASSERT and ICONIC, both arms
#     6  Duct threshold       — Axis 2, and the GALA cut-point collision
#     7  Scenario comparison  — every agent on one patient
#     8  Nutrition & growth   — Axis 5, where in the gut a drug acts
#     9  Fibrosis & portal HT
#    10  Survival             — Axes 4, 6 and 7
#    11  Biomarkers
#    12  Model notes          — the six failures, stated in the app
#
#   Run:  shiny::runApp("algs_shiny_app.R")
# =====================================================================

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

# --- locate and load the model ---------------------------------------
.here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),
                  error = function(e) ".")
.mfile <- file.path(.here, "algs_mrgsolve_model.R")
if (!file.exists(.mfile)) .mfile <- "algs_mrgsolve_model.R"
if (!exists("mod_algs")) suppressMessages(source(.mfile))

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom")

# Reference lines drawn from published data, never from the model
REF <- list(
  assert_sba_drug = 149, assert_sba_pbo = 271, assert_sba_base = 237,
  assert_itch_drug = 1.1, assert_itch_pbo = 2.2, assert_itch_base = 2.8,
  gala_nls = data.frame(age = c(5, 10, 18), nls = c(66.8, 54.4, 40.3)),
  gala_cut_low = 5.0, gala_cut_high = 10.0,
  mrx_hr = 0.305
)

# =====================================================================
#  UI
# =====================================================================
ui <- fluidPage(
  titlePanel("Alagille Syndrome — Quantitative Systems Pharmacology"),
  tags$p(style = "color:#555;margin-top:-8px;",
         tags$em(paste("JAG1/NOTCH2 haploinsufficiency to bile duct paucity.",
                       "Duct capacity is the patient; serum bile acid is only",
                       "the symptom the trials happen to measure."))),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("dpr0", "Bile duct : portal tract ratio (normal 0.9-1.8)",
                  min = 0.10, max = 0.60, value = 0.228, step = 0.005),
      sliderInput("wt", "Body weight (kg)", 5, 60, 15, 1),
      sliderInput("age0", "Age at start (yr)", 0.5, 16, 3, 0.5),
      hr(),
      h4("Therapy"),
      selectInput("drug", "Agent",
                  choices = c("none", "maralixibat 380 ug/kg/d",
                              "odevixibat 120 ug/kg/d", "UDCA 20 mg/kg/d",
                              "cholestyramine", "rifampicin", "naltrexone",
                              "PEBD (surgical diversion)"),
                  selected = "odevixibat 120 ug/kg/d"),
      checkboxInput("mct", "MCT-enriched formula", FALSE),
      sliderInput("tstart_yr", "Treatment start (yr from baseline)",
                  0, 10, 0, 0.5),
      hr(),
      h4("Simulation"),
      sliderInput("years", "Horizon (yr)", 0.5, 17, 3, 0.5),
      sliderInput("pbo", "Placebo / observer itch effect (points)",
                  0, 1.5, 0, 0.1),
      hr(),
      h4("Hazard structure"),
      sliderInput("nexp", "Exposure-hazard exponent n",
                  1.0, 3.0, 1.6, 0.05),
      helpText(paste("GALA's own bilirubin strata imply n = 1.5-1.85.",
                     "Reaching the published maralixibat HR of 0.305 from",
                     "the drug's bile-acid effect needs n ~ 2.9, which then",
                     "over-predicts GALA's severe stratum 5.6-fold.")),
      hr(),
      helpText("Sweep tabs (6, 7) re-simulate many patients and take a few
                seconds.")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 Patient", br(),
                 fluidRow(column(12, DTOutput("tbl_profile"))),
                 br(), plotOutput("p_profile", height = "320px"),
                 helpText("PHI_EHC is the share of the hepatocyte's bile-acid
                           input that came back from the gut. It is the hard
                           ceiling on any IBAT inhibitor, it is set by duct
                           capacity, and it is 0.98 in an unobstructed liver.")),
        tabPanel("2 Drug PK", br(),
                 plotOutput("p_pk", height = "300px"),
                 plotOutput("p_inh", height = "260px"),
                 helpText("Maralixibat and odevixibat are <1% absorbed. The
                           pharmacologically active compartment is the drug
                           sitting in the terminal ileal lumen.")),
        tabPanel("3 Bile acids", br(),
                 plotOutput("p_ba", height = "340px"),
                 plotOutput("p_flux", height = "280px"),
                 helpText("At steady state synthesis must equal faecal plus
                           renal loss. Serum bile acid is whatever the ducts
                           could not carry.")),
        tabPanel("4 Pruritus", br(),
                 plotOutput("p_itch", height = "320px"),
                 DTOutput("tbl_itch"),
                 helpText("The reported 0-4 scale is a ceiling on the
                           instrument, not on the sensation. Forcing the
                           operating point onto a saturating 0-4 map makes the
                           trial demand a supralinear Stevens exponent.")),
        tabPanel("5 Trial replication", br(),
                 h4("ASSERT (odevixibat, n=52) — both arms"),
                 DTOutput("tbl_assert"),
                 h4("ASSERT as published"), DTOutput("tbl_assert_obs"),
                 h4("Itch per umol/L of bile acid"), DTOutput("tbl_slopes"),
                 helpText("The two single-arm slopes agree with each other and
                           are both about twice the placebo-controlled slope.")),
        tabPanel("6 Duct threshold", br(),
                 actionButton("go_sweep", "Run duct-capacity sweep",
                              class = "btn-primary"),
                 br(), br(), plotOutput("p_sweep", height = "360px"),
                 DTOutput("tbl_threshold"),
                 helpText("Dashed vertical lines are GALA's prognostic
                           bilirubin cut-points (5.0 and 10.0 mg/dL). This
                           model never saw them.")),
        tabPanel("7 Scenarios", br(),
                 actionButton("go_panel", "Run drug panel (24 weeks)",
                              class = "btn-primary"),
                 br(), br(), DTOutput("tbl_panel"),
                 plotOutput("p_panel", height = "340px")),
        tabPanel("8 Nutrition", br(),
                 plotOutput("p_nutr", height = "320px"),
                 plotOutput("p_growth", height = "260px"),
                 helpText("ASBT is in the terminal ileum, downstream of the
                           duodenal micellar window. A sequestrant binds bile
                           acid inside that window. Both remove bile acids;
                           only one causes steatorrhoea.")),
        tabPanel("9 Fibrosis", br(),
                 plotOutput("p_fib", height = "320px"),
                 plotOutput("p_pht", height = "260px")),
        tabPanel("10 Survival", br(),
                 plotOutput("p_surv", height = "340px"),
                 h4("Competing risks to age 18"), DTOutput("tbl_risk"),
                 helpText("Population curves use a gamma frailty of variance
                           2.94, without which no single trajectory fits GALA.
                           Individual-level curves are not comparable to a
                           cohort study.")),
        tabPanel("11 Biomarkers", br(),
                 plotOutput("p_bio", height = "420px")),
        tabPanel("12 Model notes", br(), uiOutput("notes"))
      )
    )
  )
)

# =====================================================================
#  SERVER
# =====================================================================
server <- function(input, output, session) {

  drug_par <- reactive({
    ts <- yr(input$tstart_yr); te <- yr(input$years)
    base <- list(TSTART = ts, TSTOP = te, MCTF = as.numeric(input$mct))
    switch(input$drug,
      "none"                     = base,
      "maralixibat 380 ug/kg/d"  = c(base, list(DOSEKG = 380, POT = 0.28)),
      "odevixibat 120 ug/kg/d"   = c(base, list(DOSEKG = 120, POT = 1.00)),
      "UDCA 20 mg/kg/d"          = c(base, list(UDCAF = 1)),
      "cholestyramine"           = c(base, list(CHOLB = 0.55)),
      "rifampicin"               = c(base, list(RIFF = 1)),
      "naltrexone"               = c(base, list(NALF = 1)),
      "PEBD (surgical diversion)"= c(base, list(PEBDF = 1)))
  })

  run <- reactive({
    end <- yr(input$years)
    do.call(sim_eq, c(list(end = end, delta = max(1, end / 400),
                           DPR0 = input$dpr0, WT = input$wt, AGE0 = input$age0,
                           NEXP = input$nexp, PBOMAX = input$pbo),
                      drug_par()))
  })

  run_ctrl <- reactive({
    end <- yr(input$years)
    sim_eq(end = end, delta = max(1, end / 400), DPR0 = input$dpr0,
           WT = input$wt, AGE0 = input$age0, NEXP = input$nexp,
           PBOMAX = input$pbo)
  })

  both <- reactive({
    bind_rows(mutate(run_ctrl(), arm = "untreated"),
              mutate(run(),      arm = input$drug))
  })

  # ---- 1 patient ----------------------------------------------------
  output$tbl_profile <- renderDT({
    d <- run_ctrl(); t0 <- 0
    data.frame(
      quantity = c("duct : portal tract ratio", "duct capacity JCAP (x synthesis)",
                   "PHI_EHC (drug-accessible fraction)",
                   "serum bile acid (umol/L)", "total bilirubin (mg/dL)",
                   "GGT (U/L)", "ALT (U/L)", "cholesterol (mg/dL)",
                   "itch score (0-4)", "fat absorption", "25-OH vitamin D (ng/mL)",
                   "height z-score", "fibrosis stage"),
      baseline = round(c(at_day(d, "DPR", t0), at_day(d, "JCAP_OUT", t0),
                         at_day(d, "PHI_OUT", t0), at_day(d, "SBA_OUT", t0),
                         at_day(d, "BILI", t0), at_day(d, "GGTX", t0),
                         at_day(d, "ALTX", t0), at_day(d, "TCHOL_OUT", t0),
                         at_day(d, "ITCH", t0), at_day(d, "FATABS_O", t0),
                         at_day(d, "VITD", t0), at_day(d, "HTZ", t0),
                         at_day(d, "FIB", t0)), 3),
      unobstructed_liver = c("0.9 - 1.8", "~55", "0.98", "< 10", "< 1.0",
                             "< 25", "< 40", "< 200", "0", "0.90", "> 30",
                             "0", "0")
    )
  }, options = list(dom = "t", pageLength = 15), rownames = FALSE)

  output$p_profile <- renderPlot({
    d <- both()
    d %>% select(time, arm, SBA_OUT, BILI, ITCH, JCAP_OUT) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = recode(name, SBA_OUT = "serum bile acid (umol/L)",
                           BILI = "bilirubin (mg/dL)", ITCH = "itch (0-4)",
                           JCAP_OUT = "duct capacity JCAP")) %>%
      ggplot(aes(time / 365.25 + input$age0, value, colour = arm)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "age (yr)", y = NULL, colour = NULL,
           title = "Core trajectories") + THEME
  })

  # ---- 2 PK ---------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- run()
    d %>% select(time, DEPOT, DGUT, DILE, DSYS) %>%
      pivot_longer(-time) %>%
      mutate(name = factor(name, c("DEPOT", "DGUT", "DILE", "DSYS"),
                           c("upper gut depot", "mid gut",
                             "terminal ileum (site of action)", "systemic"))) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "day", y = "drug amount (ug)", colour = NULL,
           title = "IBAT inhibitor disposition") + THEME
  })

  output$p_inh <- renderPlot({
    d <- both()
    ggplot(d, aes(time, INH_OUT, colour = arm)) + geom_line(linewidth = 0.9) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "day", y = "fractional ASBT blockade", colour = NULL,
           title = "Target engagement") + THEME
  })

  # ---- 3 bile acids -------------------------------------------------
  output$p_ba <- renderPlot({
    both() %>%
      select(time, arm, BAHEP, BABIL, BADUO, BAILE, BACOL, BASYS) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = factor(name,
        c("BAHEP", "BABIL", "BADUO", "BAILE", "BACOL", "BASYS"),
        c("hepatocyte", "bile / gallbladder", "duodenum", "terminal ileum",
          "colon", "systemic"))) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.85) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = "bile acid (umol)", colour = NULL,
           title = "Bile acid compartments") + THEME
  })

  output$p_flux <- renderPlot({
    both() %>% select(time, arm, PHI_OUT, SIG_OUT, EXPO_OUT) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = recode(name,
        PHI_OUT = "PHI_EHC (drug-accessible fraction)",
        SIG_OUT = "synthesis reserve engaged (fold)",
        EXPO_OUT = "hepatocyte BA, relative")) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL, colour = NULL,
           title = "Why the ceiling exists") + THEME
  })

  # ---- 4 pruritus ---------------------------------------------------
  output$p_itch <- renderPlot({
    both() %>% select(time, arm, ITCH, ITCH_OBS, ATX, SLEEPD) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = recode(name, ITCH = "itch, drug-attributable",
        ITCH_OBS = "itch as reported (incl. placebo)",
        ATX = "autotaxin / LPA arm", SLEEPD = "sleep disturbance")) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL, colour = NULL,
           title = "Pruritus has two components; the drug moves one") + THEME
  })

  output$tbl_itch <- renderDT({
    d <- run(); u <- run_ctrl(); e <- max(d$time)
    data.frame(
      quantity = c("serum bile acid (umol/L)", "itch, drug-attributable",
                   "itch as reported", "autotaxin / LPA arm",
                   "central sensitisation"),
      untreated = round(c(at_day(u, "SBA_OUT", e), at_day(u, "ITCH", e),
                          at_day(u, "ITCH_OBS", e), at_day(u, "ATX", e),
                          at_day(u, "SENS", e)), 3),
      treated = round(c(at_day(d, "SBA_OUT", e), at_day(d, "ITCH", e),
                        at_day(d, "ITCH_OBS", e), at_day(d, "ATX", e),
                        at_day(d, "SENS", e)), 3))
  }, options = list(dom = "t"), rownames = FALSE)

  # ---- 5 trials -----------------------------------------------------
  output$tbl_assert    <- renderDT(round_df(sc03_assert()),
                                   options = list(dom = "t"), rownames = FALSE)
  output$tbl_assert_obs<- renderDT(ASSERT_OBS, options = list(dom = "t"),
                                   rownames = FALSE)
  output$tbl_slopes    <- renderDT(round_df(sc04_slopes(), 5),
                                   options = list(dom = "t"), rownames = FALSE)

  # ---- 6 duct threshold ---------------------------------------------
  sweep <- eventReactive(input$go_sweep, {
    withProgress(message = "Sweeping duct capacity", value = 0.5,
                 sc05_duct_sweep(seq(0.10, 0.45, by = 0.025)))
  })

  output$p_sweep <- renderPlot({
    s <- sweep()
    ggplot(s, aes(BILI, pct_change)) +
      geom_hline(yintercept = c(-15, -30), linetype = 3, colour = "grey40") +
      geom_vline(xintercept = c(REF$gala_cut_low, REF$gala_cut_high),
                 linetype = 2, colour = "firebrick", linewidth = 0.8) +
      geom_line(linewidth = 1.1, colour = "#1565c0") +
      geom_point(size = 1.6, colour = "#1565c0") +
      annotate("text", x = REF$gala_cut_low, y = -88, hjust = -0.08,
               label = "GALA 5.0", colour = "firebrick", size = 3.6) +
      annotate("text", x = REF$gala_cut_high, y = -88, hjust = -0.08,
               label = "GALA 10.0", colour = "firebrick", size = 3.6) +
      labs(x = "baseline total bilirubin (mg/dL)",
           y = "24-week serum bile acid change (%)",
           title = paste("IBAT-inhibitor response collapses with duct",
                         "capacity; the boundaries land on GALA's cut-points"))+
      THEME
  })

  output$tbl_threshold <- renderDT(round_df(sc05_threshold(sweep())),
                                   options = list(dom = "t"), rownames = FALSE)

  # ---- 7 scenarios --------------------------------------------------
  panel <- eventReactive(input$go_panel, {
    withProgress(message = "Running drug panel", value = 0.5, sc08_drug_panel())
  })
  output$tbl_panel <- renderDT(round_df(panel()), options = list(dom = "t"),
                               rownames = FALSE)
  output$p_panel <- renderPlot({
    panel() %>% select(arm, sBA, itch, fat_abs, vitD) %>%
      pivot_longer(-arm) %>%
      ggplot(aes(reorder(arm, value), value, fill = name)) +
      geom_col(show.legend = FALSE) + coord_flip() +
      facet_wrap(~name, scales = "free_x") +
      labs(x = NULL, y = NULL,
           title = "24 weeks, one patient, every agent") + THEME
  })

  # ---- 8 nutrition --------------------------------------------------
  output$p_nutr <- renderPlot({
    both() %>% select(time, arm, CDUO_OUT, FATABS_O, VITD, VITE) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = recode(name, CDUO_OUT = "duodenal bile acid (mM)",
        FATABS_O = "fat absorption coefficient",
        VITD = "25-OH vitamin D (ng/mL)", VITE = "vitamin E (relative)")) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL, colour = NULL,
           title = "The micellar window") + THEME
  })

  output$p_growth <- renderPlot({
    both() %>% select(time, arm, HTZ, IGF1, INR) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = recode(name, HTZ = "height z-score",
                           IGF1 = "IGF-1 (relative)", INR = "INR")) %>%
      ggplot(aes(time / 365.25 + input$age0, value, colour = arm)) +
      geom_line(linewidth = 0.9) + facet_wrap(~name, scales = "free_y") +
      labs(x = "age (yr)", y = NULL, colour = NULL,
           title = "Growth and coagulation") + THEME
  })

  # ---- 9 fibrosis ---------------------------------------------------
  output$p_fib <- renderPlot({
    both() %>% select(time, arm, FIB, HSC, DUCTR, CUMBA) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = recode(name, FIB = "fibrosis stage",
        HSC = "activated stellate cells", DUCTR = "ductular reaction",
        CUMBA = "cumulative BA exposure")) %>%
      ggplot(aes(time / 365.25 + input$age0, value, colour = arm)) +
      geom_line(linewidth = 0.9) + facet_wrap(~name, scales = "free_y") +
      labs(x = "age (yr)", y = NULL, colour = NULL,
           title = "Fibrosis tracks the INTEGRAL of exposure") + THEME
  })

  output$p_pht <- renderPlot({
    both() %>% select(time, arm, PORTP_O, PLT) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = recode(name, PORTP_O = "portal pressure index (mmHg)",
                           PLT = "platelets (10^9/L)")) %>%
      ggplot(aes(time / 365.25 + input$age0, value, colour = arm)) +
      geom_line(linewidth = 0.9) + facet_wrap(~name, scales = "free_y") +
      labs(x = "age (yr)", y = NULL, colour = NULL,
           title = "Portal hypertension") + THEME
  })

  # ---- 10 survival --------------------------------------------------
  output$p_surv <- renderPlot({
    d <- both() %>%
      select(time, arm, NLS_POP, EFS_POP, NLS, EFS) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(level = ifelse(grepl("POP", name), "population (frailty)",
                            "individual"),
             what  = ifelse(grepl("^NLS", name), "native liver survival",
                            "event-free survival"),
             age   = time / 365.25 + input$age0)
    ggplot(d, aes(age, 100 * value, colour = arm, linetype = level)) +
      geom_line(linewidth = 0.9) + facet_wrap(~what) +
      geom_point(data = REF$gala_nls %>% mutate(what = "native liver survival"),
                 aes(age, nls), inherit.aes = FALSE, colour = "firebrick",
                 size = 2.8) +
      labs(x = "age (yr)", y = "%", colour = NULL, linetype = NULL,
           title = "Red points are GALA observed native liver survival") + THEME
  })

  output$tbl_risk <- renderDT(round_df(sc13_competing_risk()),
                              options = list(dom = "t"), rownames = FALSE)

  # ---- 11 biomarkers ------------------------------------------------
  output$p_bio <- renderPlot({
    both() %>%
      select(time, arm, SBA_OUT, BILI, GGTX, ALTX, TCHOL_OUT, XANTH, PLT, INR) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = recode(name, SBA_OUT = "serum bile acid (umol/L)",
        BILI = "total bilirubin (mg/dL)", GGTX = "GGT (U/L)",
        ALTX = "ALT (U/L)", TCHOL_OUT = "cholesterol (mg/dL)",
        XANTH = "xanthoma burden", PLT = "platelets", INR = "INR")) %>%
      ggplot(aes(time / 365.25 + input$age0, value, colour = arm)) +
      geom_line(linewidth = 0.85) + facet_wrap(~name, scales = "free_y",
                                               ncol = 3) +
      labs(x = "age (yr)", y = NULL, colour = NULL,
           title = "Biomarker panel") + THEME
  })

  # ---- 12 notes -----------------------------------------------------
  output$notes <- renderUI({
    HTML(paste0(
      "<h4>What this model does not do</h4><ol>",
      "<li><b>PHI_EHC is not separately identifiable from one trial.</b> ",
      "Duct capacity, achieved ASBT blockade and the synthesis reserve enter ",
      "the steady-state bile-acid ratio as one combination. Re-solving duct ",
      "capacity for each assumed blockade moves it 3.58 to 1.95 while ",
      "PHI_EHC moves 0.57 to 0.42.</li>",
      "<li><b>The bilirubin boundaries rest on one assumed anchor</b> (the ",
      "ASSERT population's mean baseline bilirubin, taken as 3.5 mg/dL). ",
      "Only the <i>ratio</i> of the two boundaries, 1.97 against GALA's 2.00, ",
      "is anchor-independent.</li>",
      "<li><b>The bile-acid-independent itch fraction is assumed, not ",
      "solved.</b> The ASSERT placebo arm quantifies the non-drug component ",
      "(47%); nothing published constrains the residual at zero bile acid.</li>",
      "<li><b>Spontaneous improvement is only half reproduced.</b> The duct ",
      "term alone lowers bilirubin 3.50 to 2.20 from age 1 to 18; in the full ",
      "model fibrosis outruns it and bilirubin rises to 4.48 instead.</li>",
      "<li><b>Hazards are calibrated to a cohort and applied to an ",
      "individual.</b> The frailty variance needed to fit GALA is 2.94, so the ",
      "cohort curve is largely a statement about between-patient variance.</li>",
      "<li><b>No cirrhotic pharmacology, and no placebo-arm progression.</b> ",
      "ASSERT's placebo bile acids rose 22 umol/L; this model holds its ",
      "untreated arm flat and therefore over-states the absolute on-drug ",
      "change.</li></ol>",
      "<p style='color:#777'>Research and teaching model. Not a clinical ",
      "decision tool, not validated for any patient-level use.</p>"))
  })
}

# Round every numeric column of a data frame for display
round_df <- function(d, digits = 3) {
  num <- vapply(d, is.numeric, logical(1))
  d[num] <- lapply(d[num], round, digits)
  d
}

shinyApp(ui, server)
