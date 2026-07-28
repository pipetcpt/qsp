# =============================================================================
#  ppd_shiny_app.R
#  Postpartum depression QSP model — interactive dashboard
#  산후우울증 QSP 모델 · Shiny 대시보드 (9 tabs)
#  ---------------------------------------------------------------------------
#  Run:
#      # from the repository root
#      shiny::runApp("postpartum-depression/ppd_shiny_app.R")
#
#  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
#
#  The app compiles the model defined in ppd_mrgsolve_model.R.  To avoid
#  re-running that file's scenario block on startup, the model code string is
#  extracted from it rather than sourced.
#
#  DESIGN INTENT
#  -------------
#  The dashboard is built so that the model's central claim is visible on the
#  FIRST tab a user opens, and so that the two things the model is most
#  uncertain about (KR, THR0) are exposed as sliders rather than buried:
#  a user should be able to break this model in under a minute.
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

# ---------------------------------------------------------------------------
# Build the model by extracting the code string from the model file
# ---------------------------------------------------------------------------
MODEL_FILE <- "ppd_mrgsolve_model.R"
if (!file.exists(MODEL_FILE))
  MODEL_FILE <- file.path("postpartum-depression", "ppd_mrgsolve_model.R")

build_model <- function(path = MODEL_FILE) {
  txt <- readLines(path, warn = FALSE)
  i0 <- grep("^ppd_code <- '", txt)[1]
  i1 <- grep("^'$", txt)
  i1 <- i1[i1 > i0][1]
  code <- paste(txt[(i0 + 1):(i1 - 1)], collapse = "\n")
  mcode_cache("ppd_shiny", code)
}
mod <- build_model()

TDEL <- 1000                                  # simulation time of delivery (h)

# ---------------------------------------------------------------------------
# Event builders (mirror of the helpers in ppd_mrgsolve_model.R)
# ---------------------------------------------------------------------------
ev_delivery <- function() ev(time = TDEL, amt = 22, cmt = "INFL")

ev_brexanolone <- function(top, wt, start) {
  rates <- c(30, min(60, top), top, 30); durs <- c(4, 20, 28, 8)
  t0 <- start; out <- NULL
  for (i in seq_along(rates)) {
    e <- ev(time = t0, amt = rates[i] * wt * durs[i],
            rate = rates[i] * wt, cmt = "BRX1")
    out <- if (is.null(out)) e else c(out, e); t0 <- t0 + durs[i]
  }
  out
}
ev_zuranolone <- function(dose, days, start)
  ev(time = start, amt = dose, cmt = "ZURA", ii = 24, addl = max(0, days - 1))
ev_sertraline <- function(start, days)
  c(ev(time = start, amt = 50, cmt = "SERA", ii = 24, addl = 6),
    ev(time = start + 7 * 24, amt = 100, cmt = "SERA", ii = 24,
       addl = max(0, days - 8)))
ev_esketamine <- function(mgkg, wt, start)
  ev(time = start, amt = mgkg * wt, rate = mgkg * wt / (2/3), cmt = "ESKC")

OFF <- 1e9
win <- function(care = NULL, protect = NULL, nsp = NULL, cbt = NULL) {
  p <- list(CARE_T0 = OFF, CARE_T1 = OFF, PROT_T0 = OFF, PROT_T1 = OFF,
            NSP_T0 = OFF, CBT_T0 = OFF)
  if (!is.null(care))    { p$CARE_T0 <- care[1];    p$CARE_T1 <- care[2] }
  if (!is.null(protect)) { p$PROT_T0 <- protect[1]; p$PROT_T1 <- protect[2] }
  if (!is.null(nsp))       p$NSP_T0  <- nsp
  if (!is.null(cbt))       p$CBT_T0  <- cbt
  p
}

# Published endpoints, for the comparison tab.  Every value verified on PubMed;
# see ppd_references.md.
TRIALS <- tibble::tribble(
  ~trial,                        ~arm,                  ~day, ~observed, ~pmid,
  "Kanes 2017 (phase 2)",        "brexanolone",          2.5, -21.0, "28619476",
  "Kanes 2017 (phase 2)",        "placebo (inpatient)",  2.5,  -8.8, "28619476",
  "Meltzer-Brody 2018 study 1",  "brexanolone 60",       2.5, -19.5, "30177236",
  "Meltzer-Brody 2018 study 1",  "brexanolone 90",       2.5, -17.7, "30177236",
  "Meltzer-Brody 2018 study 1",  "placebo (inpatient)",  2.5, -14.0, "30177236",
  "Meltzer-Brody 2018 study 2",  "brexanolone 90",       2.5, -14.6, "30177236",
  "Meltzer-Brody 2018 study 2",  "placebo (inpatient)",  2.5, -12.1, "30177236",
  "ROBIN 2021",                  "zuranolone 30 x 14 d",  15, -17.8, "34190962",
  "ROBIN 2021",                  "placebo (outpatient)",  15, -13.6, "34190962",
  "SKYLARK 2023",                "zuranolone 50 x 14 d",  15, -15.6, "37491938",
  "SKYLARK 2023",                "placebo (outpatient)",  15, -11.6, "37491938"
)

ARMS <- c("none (natural history)", "placebo (inpatient)",
          "placebo (outpatient)", "brexanolone 60", "brexanolone 90",
          "zuranolone 50 x 14 d", "zuranolone 30 x 14 d",
          "zuranolone 50 x 3 d", "sertraline 50-100 mg",
          "esketamine 0.25 mg/kg", "protected sleep", "CBT/IPT")

run_arm <- function(arm, V, enrol_day, wt, extra = list(), horizon_days = 46) {
  enrol <- TDEL + enrol_day * 24
  end   <- enrol + horizon_days * 24
  e <- ev_delivery(); p <- win()
  inpatient <- function() win(care = c(enrol, enrol + 60),
                              protect = c(enrol, enrol + 60), nsp = enrol)
  if (arm == "placebo (inpatient)")        p <- inpatient()
  if (arm == "placebo (outpatient)")       p <- win(nsp = enrol)
  if (arm == "brexanolone 60")  { p <- inpatient(); e <- c(e, ev_brexanolone(60, wt, enrol)) }
  if (arm == "brexanolone 90")  { p <- inpatient(); e <- c(e, ev_brexanolone(90, wt, enrol)) }
  if (arm == "zuranolone 50 x 14 d") { p <- win(nsp = enrol); e <- c(e, ev_zuranolone(50, 14, enrol)) }
  if (arm == "zuranolone 30 x 14 d") { p <- win(nsp = enrol); e <- c(e, ev_zuranolone(30, 14, enrol)) }
  if (arm == "zuranolone 50 x 3 d")  { p <- win(nsp = enrol); e <- c(e, ev_zuranolone(50, 3, enrol)) }
  if (arm == "sertraline 50-100 mg") { p <- win(nsp = enrol); e <- c(e, ev_sertraline(enrol, horizon_days)) }
  if (arm == "esketamine 0.25 mg/kg"){ p <- win(nsp = enrol); e <- c(e, ev_esketamine(0.25, wt, enrol)) }
  if (arm == "protected sleep")      p <- win(protect = c(enrol, end), nsp = enrol)
  if (arm == "CBT/IPT")              p <- win(nsp = enrol, cbt = enrol)

  mod %>%
    param(c(list(V = V, WT = wt, TDEL = TDEL), p, extra)) %>%
    ev(e) %>%
    mrgsim(end = end, delta = 1) %>%
    as_tibble() %>%
    mutate(arm = arm, enrol_time = enrol,
           day_from_enrol = (time - enrol) / 24)
}

d_hamd <- function(d, days) {
  base <- d$HAMD[which.min(abs(d$time - d$enrol_time[1]))]
  sapply(days, function(dd)
    d$HAMD[which.min(abs(d$time - (d$enrol_time[1] + dd * 24)))] - base)
}

# ===========================================================================
# UI
# ===========================================================================
ui <- fluidPage(
  titlePanel("산후우울증 QSP 모델 · Postpartum Depression QSP Dashboard"),
  tags$p(tags$b("The model's one idea:"),
         " tonic inhibition is a PRODUCT of a ligand that collapses in hours ",
         "and a receptor pool that re-expands over weeks. ",
         tags$i("A product of a fast fall and a slow rise troughs.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("1 · 환자 · Patient"),
      sliderInput("V", "취약성 V · vulnerability gain", 0.6, 2.2, 1.70, 0.05),
      helpText("prior PPD / prior MDD / low support raise V. ",
               "V scales symptom accrual, slows receptor recovery and ",
               "narrows the E/I reserve."),
      sliderInput("enrol_day", "등록 시점 (분만 후 일) · enrolment day",
                  3, 180, 21, 1),
      sliderInput("wt", "체중 kg · maternal weight", 45, 120, 70, 1),
      hr(),
      h4("2 · 치료군 · Arms"),
      checkboxGroupInput("arms", NULL, choices = ARMS,
                         selected = c("placebo (inpatient)",
                                      "brexanolone 60",
                                      "zuranolone 50 x 14 d",
                                      "sertraline 50-100 mg")),
      hr(),
      h4("3 · 모델을 부숴보기 · break the model"),
      helpText("These two parameters have no human measurement and dominate ",
               "the model. Move them."),
      sliderInput("KR_mult", "KR × (δ-pool recovery rate; 0 = frozen)",
                  0, 3, 1, 0.05),
      sliderInput("THR0", "THR0 (E/I threshold at V = 1)", 1.5, 4.5, 2.85, 0.05),
      sliderInput("EC50_PAM", "EC50 potentiation (nM brain ALLO-eq)",
                  40, 400, 120, 5),
      hr(),
      actionButton("reset", "reset to reference", class = "btn-sm")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",

        # ---------------------------------------------------------------
        tabPanel(
          "① 기전 · The trough",
          h4("Why a transient exists where no steady state does"),
          plotOutput("p_trough", height = "460px"),
          hr(),
          fluidRow(
            column(6, h5("Steady states all coincide"),
                   tableOutput("t_setpoint")),
            column(6, h5("Trough geometry"), verbatimTextOutput("txt_trough"))
          )
        ),

        # ---------------------------------------------------------------
        tabPanel(
          "② PK · 신경스테로이드",
          h4("Endogenous and exogenous PAM in one currency"),
          plotOutput("p_pk", height = "380px"),
          hr(),
          h5("Saturation of potentiation with brexanolone infusion rate"),
          tableOutput("t_saturation"),
          helpText("A 50 % increase from 60 to 90 µg/kg/h buys a single-digit ",
                   "percent increase in potentiation. Site-2 direct gating ",
                   "(sedation) has no such ceiling — which is the model's ",
                   "account of the flat 60-vs-90 dose-response in phase 3.")
        ),

        # ---------------------------------------------------------------
        tabPanel(
          "③ 수용체 · Receptor & reserve",
          h4("δ-pool, tonic conductance, deficit and the E/I reserve"),
          plotOutput("p_receptor", height = "460px"),
          helpText("RESERVE = threshold − EXC. Positive means sub-threshold: ",
                   "symptoms resolve at the accelerated rate. The model spends ",
                   "most of its interesting time near RESERVE = 0, which is ",
                   "deliberate (assumption A7) and is why trajectories near ",
                   "the threshold move slowly.")
        ),

        # ---------------------------------------------------------------
        tabPanel(
          "④ 임상 엔드포인트 · HAM-D / EPDS",
          h4("HAM-D17 and EPDS by arm"),
          plotOutput("p_endpoint", height = "400px"),
          hr(),
          h5("Change from enrolment"),
          DTOutput("t_endpoint")
        ),

        # ---------------------------------------------------------------
        tabPanel(
          "⑤ 시나리오 비교 · vs published trials",
          h4("Model versus the published endpoints it was calibrated to"),
          DTOutput("t_trials"),
          helpText("Observed values are from PMID 28619476, 30177236, ",
                   "34190962 and 37491938. The model was fitted to a subset ",
                   "of these; agreement on the rest is not independent ",
                   "validation, and the enrolment-day slider will move every ",
                   "model column (assumption A6)."),
          hr(),
          h5("Placebo response versus enrolment day (a trial-design prediction)"),
          plotOutput("p_enrol", height = "320px")
        ),

        # ---------------------------------------------------------------
        tabPanel(
          "⑥ 바이오마커 · Biomarkers",
          h4("HPA axis, MAO-A, inflammation, kynurenine, sleep, plasticity"),
          plotOutput("p_biomarker", height = "520px"),
          helpText("Cortisol falls after delivery in this model because ",
                   "hypothalamic CRH was suppressed by months of placental ",
                   "CRH-driven hypercortisolism and recovers with a ~29-day ",
                   "half-life (PMID 8626857).")
        ),

        # ---------------------------------------------------------------
        tabPanel(
          "⑦ 수유·영아 노출 · Lactation",
          h4("Maternal exposure and infant exposure are different questions"),
          plotOutput("p_milk", height = "380px"),
          hr(),
          tableOutput("t_milk"),
          helpText("The relative infant dose is computed from a milk ",
                   "compartment and an infant PK compartment, and the infant ",
                   "plasma concentration is compared with the model's own ",
                   "potentiation EC50 — not with the maternal dose. Measured ",
                   "brexanolone milk data exist (PMID 35869362) and should be ",
                   "used to falsify the RID column.")
        ),

        # ---------------------------------------------------------------
        tabPanel(
          "⑧ 브릿지 검정 · Bridge test",
          h4("Is durability the drug, or the receptor?"),
          helpText("Set the KR slider to 0 in the sidebar to freeze δ-pool ",
                   "plasticity. The prediction: the on-drug response survives ",
                   "and the post-washout benefit does not, with no change to ",
                   "any symptom-equation parameter."),
          plotOutput("p_bridge", height = "400px"),
          hr(),
          tableOutput("t_bridge")
        ),

        # ---------------------------------------------------------------
        tabPanel(
          "⑨ 가정과 한계 · Assumptions",
          h4("Read this before quoting a number"),
          tags$ol(
            tags$li(tags$b("A1"), " Neurosteroids are total plasma/brain ",
                    "concentrations and EC50 is in the same currency, so >99 % ",
                    "protein binding is lumped into EC50."),
            tags$li(tags$b("A2"), " Brexanolone IS allopregnanolone, so its ",
                    "Kp and potency are set by identity, not fitted."),
            tags$li(tags$b("A3"), " ZUR_EQ is the only parameter calibrated ",
                    "against a zuranolone endpoint."),
            tags$li(tags$b("A4"), " Term ALLO ≈ 80 nM and the postpartum floor ",
                    "≈ 2 nM; assays disagree severalfold (PMID 11238543) and a ",
                    "2025 IPD meta-analysis finds absolute concentrations do ",
                    "not separate cases from controls (PMID 39511449). The ",
                    "~40× RATIO is what was matched."),
            tags$li(tags$b("A5"), " K_CARE and K_NSP are non-mechanistic ",
                    "expectancy/contact terms fitted once to placebo arms."),
            tags$li(tags$b("A6"), " Reference enrolment is day 21 postpartum; ",
                    "real trials enrolled far more widely."),
            tags$li(tags$b("A7"), " The model operates near a bifurcation, so ",
                    "late endpoints are intrinsically sensitive to THR0 and KR."),
            tags$li(tags$b("A8"), " KR has no human measurement and is the ",
                    "single most influential parameter.")
          ),
          hr(),
          h5("Sensitivity of the day-15 zuranolone endpoint (±30 %)"),
          actionButton("run_sens", "run sensitivity (slow)", class = "btn-sm"),
          DTOutput("t_sens"),
          hr(),
          tags$p("Model files: ",
                 tags$code("ppd_qsp_model.dot"), " · ",
                 tags$code("ppd_mrgsolve_model.R"), " · ",
                 tags$code("ppd_reference_check.py"), " · ",
                 tags$code("ppd_references.md")),
          tags$p(tags$b("교육·연구 목적 전용. "),
                 "This is a research and teaching model. It has not been ",
                 "validated against patient-level data and must not inform ",
                 "clinical decisions.")
        )
      )
    )
  )
)

# ===========================================================================
# SERVER
# ===========================================================================
server <- function(input, output, session) {

  observeEvent(input$reset, {
    updateSliderInput(session, "V", value = 1.70)
    updateSliderInput(session, "enrol_day", value = 21)
    updateSliderInput(session, "KR_mult", value = 1)
    updateSliderInput(session, "THR0", value = 2.85)
    updateSliderInput(session, "EC50_PAM", value = 120)
  })

  extra <- reactive({
    kr0 <- as.numeric(param(mod)[["KR"]])
    list(KR = kr0 * input$KR_mult, THR0 = input$THR0,
         EC50_PAM = input$EC50_PAM)
  })

  sims <- reactive({
    req(length(input$arms) > 0)
    bind_rows(lapply(input$arms, run_arm, V = input$V,
                     enrol_day = input$enrol_day, wt = input$wt,
                     extra = extra()))
  })

  natural <- reactive(
    run_arm("none (natural history)", input$V, input$enrol_day, input$wt,
            extra())
  )

  # ---- ① trough ---------------------------------------------------------
  output$p_trough <- renderPlot({
    natural() %>%
      filter(DAYPP >= -3, DAYPP <= 60) %>%
      transmute(DAYPP,
                `plasma ALLO (nM)` = ALLOP,
                `brain PAM (nM)` = PAM,
                `δ pool R` = RD,
                `tonic G (set-point = 1)` = GTONIC,
                `tonic deficit (%)` = DEFICIT,
                `HAM-D17` = HAMD) %>%
      pivot_longer(-DAYPP) %>%
      ggplot(aes(DAYPP, value)) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
      geom_line(linewidth = 0.9, colour = "#1b7837") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "days postpartum", y = NULL) +
      theme_bw(base_size = 12)
  })

  output$t_setpoint <- renderTable({
    p <- as.list(param(mod)); p$EC50_PAM <- input$EC50_PAM
    pot <- function(pam) p$EMAX_PAM * pam^p$H_PAM /
      (p$EC50_PAM^p$H_PAM + pam^p$H_PAM)
    pam <- c(6, 60, 160, 4.8)
    tibble(state = c("non-pregnant", "2nd trimester", "term",
                     "day 3 postpartum"),
           `PAM (nM)` = pam,
           potentiation = round(pot(pam), 3),
           `R_δ set-point` = round(pmin(1, 1 / (1 + pot(pam))), 3),
           `G_tonic (ss)` = round(pmin(1, 1 / (1 + pot(pam))) *
                                    (1 + pot(pam)), 3))
  }, digits = 3)

  output$txt_trough <- renderPrint({
    d <- natural() %>% filter(DAYPP >= 0)
    worst <- d[which.max(d$DEFICIT), ]
    cat(sprintf("worst tonic deficit      : %.0f %% at day %.1f\n",
                worst$DEFICIT, worst$DAYPP))
    for (thr in c(25, 10, 5)) {
      sel <- d$DAYPP[d$DEFICIT > thr]
      cat(sprintf("deficit > %2d %% until      : day %.1f\n", thr,
                  if (length(sel)) max(sel) else NA))
    }
    pk <- d[which.max(d$HAMD), ]
    cat(sprintf("peak HAM-D               : %.1f at day %.1f\n",
                pk$HAMD, pk$DAYPP))
    cat(sprintf("peak EPDS                : %.1f  (screening cut-off 12/13)\n",
                max(d$EPDS)))
  })

  # ---- ② PK -------------------------------------------------------------
  output$p_pk <- renderPlot({
    sims() %>%
      filter(day_from_enrol >= -1, day_from_enrol <= 20) %>%
      transmute(day_from_enrol, arm,
                `brain PAM (nM ALLO-eq)` = PAM,
                `brexanolone (ng/mL)` = CBRX,
                `zuranolone (ng/mL)` = CZUR,
                `sertraline (ng/mL)` = CSER) %>%
      pivot_longer(c(-day_from_enrol, -arm)) %>%
      ggplot(aes(day_from_enrol, value, colour = arm)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "days from enrolment", y = NULL) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom")
  })

  output$t_saturation <- renderTable({
    p <- as.list(param(mod)); p$EC50_PAM <- input$EC50_PAM
    rate <- c(30, 60, 90, 120)
    css <- rate * input$wt / p$BRX_CL
    nM <- css * 1000 / p$MW_ALLO
    brain <- p$KP_ALLO * nM
    pot <- p$EMAX_PAM * brain^p$H_PAM / (p$EC50_PAM^p$H_PAM + brain^p$H_PAM)
    tibble(`rate µg/kg/h` = rate, `Css ng/mL` = round(css, 1),
           `Css nM` = round(nM, 0), `brain nM-eq` = round(brain, 0),
           potentiation = round(pot, 3),
           `% vs 60 µg/kg/h` = round(100 * (pot / pot[2] - 1), 1))
  }, digits = 3)

  # ---- ③ receptor -------------------------------------------------------
  output$p_receptor <- renderPlot({
    sims() %>%
      filter(day_from_enrol >= -1, day_from_enrol <= 46) %>%
      transmute(day_from_enrol, arm,
                `δ pool R` = RD, `R_δ set-point` = RDSET,
                `tonic G` = GTONIC, `effective inhibition` = GEFF,
                `E/I index EXC` = EXC, `reserve (THR − EXC)` = RESERVE) %>%
      pivot_longer(c(-day_from_enrol, -arm)) %>%
      ggplot(aes(day_from_enrol, value, colour = arm)) +
      geom_hline(yintercept = 0, linetype = 3, colour = "grey50") +
      geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "days from enrolment", y = NULL) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom")
  })

  # ---- ④ endpoints ------------------------------------------------------
  output$p_endpoint <- renderPlot({
    sims() %>%
      filter(day_from_enrol >= -1, day_from_enrol <= 46) %>%
      select(day_from_enrol, arm, HAMD, EPDS) %>%
      pivot_longer(c(HAMD, EPDS)) %>%
      ggplot(aes(day_from_enrol, value, colour = arm)) +
      geom_hline(data = tibble(name = c("HAMD", "EPDS"), y = c(7, 12.5)),
                 aes(yintercept = y), linetype = 3) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "days from enrolment", y = NULL,
           caption = "dotted: HAM-D remission (≤7) and EPDS cut-off (12/13)") +
      theme_bw(base_size = 12) + theme(legend.position = "bottom")
  })

  output$t_endpoint <- renderDT({
    days <- c(1, 2.5, 3, 7, 15, 28, 45)
    d <- sims()
    out <- lapply(split(d, d$arm), function(x) {
      base <- x$HAMD[which.min(abs(x$time - x$enrol_time[1]))]
      c(base = base, d_hamd(x, days))
    })
    m <- as.data.frame(do.call(rbind, out))
    names(m) <- c("baseline", paste0("d", days))
    datatable(round(m, 1), options = list(dom = "t", paging = FALSE))
  })

  # ---- ⑤ trials ---------------------------------------------------------
  output$t_trials <- renderDT({
    d <- sims()
    model_val <- function(arm, day) {
      x <- d[d$arm == arm, ]
      if (!nrow(x)) return(NA_real_)
      d_hamd(x, day)
    }
    tt <- TRIALS %>%
      rowwise() %>%
      mutate(model = model_val(arm, day)) %>%
      ungroup() %>%
      mutate(across(c(observed, model), ~ round(.x, 1)),
             `model − obs` = round(model - observed, 1))
    datatable(tt, options = list(dom = "t", paging = FALSE))
  })

  output$p_enrol <- renderPlot({
    days <- c(7, 14, 21, 42, 90, 180)
    res <- lapply(days, function(dd) {
      a <- run_arm("placebo (inpatient)", input$V, dd, input$wt, extra(), 20)
      b <- run_arm("placebo (outpatient)", input$V, dd, input$wt, extra(), 20)
      tibble(enrol_day = dd,
             `inpatient, 60 h` = d_hamd(a, 2.5),
             `outpatient, day 15` = d_hamd(b, 15))
    })
    bind_rows(res) %>%
      pivot_longer(-enrol_day) %>%
      ggplot(aes(enrol_day, value, colour = name)) +
      geom_line(linewidth = 0.9) + geom_point() +
      labs(x = "enrolment day postpartum", y = "placebo HAM-D change",
           colour = NULL) +
      theme_bw(base_size = 12)
  })

  # ---- ⑥ biomarkers -----------------------------------------------------
  output$p_biomarker <- renderPlot({
    sims() %>%
      filter(day_from_enrol >= -1, day_from_enrol <= 46) %>%
      transmute(day_from_enrol, arm,
                `cortisol (µg/dL)` = CORT,
                `hypothalamic CRH drive` = HCRH,
                `MAO-A binding (rel)` = MAOA,
                `5-HT tone (rel)` = FIVEHT,
                `IL-6 eq (pg/mL)` = INFL,
                `KYN/TRP (rel)` = KYNR,
                `sleep (h/night)` = SLP,
                `sleep debt (h)` = SDEBT,
                `BDNF (rel)` = BDNF) %>%
      pivot_longer(c(-day_from_enrol, -arm)) %>%
      ggplot(aes(day_from_enrol, value, colour = arm)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "days from enrolment", y = NULL) +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })

  # ---- ⑦ lactation ------------------------------------------------------
  output$p_milk <- renderPlot({
    sims() %>%
      filter(day_from_enrol >= -1, day_from_enrol <= 20) %>%
      transmute(day_from_enrol, arm,
                `milk concentration (ng/mL)` = CMILK,
                `cumulative infant dose (µg/kg)` = MILKD,
                `infant plasma (nM ALLO-eq)` = INFP) %>%
      pivot_longer(c(-day_from_enrol, -arm)) %>%
      ggplot(aes(day_from_enrol, value, colour = arm)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "days from enrolment", y = NULL) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom")
  })

  output$t_milk <- renderTable({
    d <- sims()
    ec50 <- input$EC50_PAM
    bind_rows(lapply(split(d, d$arm), function(x) {
      tibble(arm = x$arm[1],
             `max milk ng/mL` = max(x$CMILK),
             `infant dose µg/kg (14 d)` =
               x$MILKD[which.min(abs(x$day_from_enrol - 14))],
             `max infant plasma nM` = max(x$INFP),
             `× below EC50` = ifelse(max(x$INFP) > 0,
                                     ec50 / max(x$INFP), NA))
    }))
  }, digits = 3)

  # ---- ⑧ bridge ---------------------------------------------------------
  bridge <- reactive({
    kr0 <- as.numeric(param(mod)[["KR"]])
    base_extra <- extra()
    frozen <- base_extra; frozen$KR <- 0
    arms <- intersect(input$arms, c("brexanolone 60", "brexanolone 90",
                                    "zuranolone 50 x 14 d",
                                    "zuranolone 30 x 14 d",
                                    "zuranolone 50 x 3 d"))
    if (!length(arms)) arms <- "zuranolone 50 x 14 d"
    bind_rows(
      bind_rows(lapply(arms, run_arm, V = input$V, enrol_day = input$enrol_day,
                       wt = input$wt, extra = base_extra)) %>%
        mutate(receptor = "KR as set"),
      bind_rows(lapply(arms, run_arm, V = input$V, enrol_day = input$enrol_day,
                       wt = input$wt, extra = frozen)) %>%
        mutate(receptor = "KR = 0 (frozen)")
    )
  })

  output$p_bridge <- renderPlot({
    bridge() %>%
      filter(day_from_enrol >= -1, day_from_enrol <= 46) %>%
      ggplot(aes(day_from_enrol, HAMD, colour = arm, linetype = receptor)) +
      geom_hline(yintercept = 7, linetype = 3) +
      geom_line(linewidth = 0.9) +
      labs(x = "days from enrolment", y = "HAM-D17", colour = NULL,
           linetype = NULL) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom")
  })

  output$t_bridge <- renderTable({
    days <- c(2.5, 3, 15, 45)
    d <- bridge()
    bind_rows(lapply(split(d, list(d$arm, d$receptor), drop = TRUE),
                     function(x) {
                       v <- d_hamd(x, days)
                       tibble(arm = x$arm[1], receptor = x$receptor[1],
                              d2.5 = v[1], d3 = v[2], d15 = v[3], d45 = v[4])
                     }))
  }, digits = 1)

  # ---- ⑨ sensitivity ---------------------------------------------------
  sens <- eventReactive(input$run_sens, {
    keys <- c("KR", "EC50_PAM", "EMAX_PAM", "ZUR_EQ", "THR0", "KON", "KOFF",
              "KFAST", "W_SD", "A_SYMP", "KDEC_SD", "K_NSP", "V_KR")
    withProgress(message = "running sensitivity", {
      bind_rows(lapply(keys, function(k) {
        incProgress(1 / length(keys), detail = k)
        base <- as.numeric(param(mod)[[k]])
        v <- sapply(c(0.7, 1.3), function(mult) {
          e <- extra(); e[[k]] <- base * mult
          d_hamd(run_arm("zuranolone 50 x 14 d", input$V, input$enrol_day,
                         input$wt, e, 16), 15)
        })
        tibble(parameter = k, `−30 %` = round(v[1], 1),
               `+30 %` = round(v[2], 1), span = round(abs(v[2] - v[1]), 2))
      })) %>% arrange(desc(span))
    })
  })
  output$t_sens <- renderDT(
    datatable(sens(), options = list(dom = "t", paging = FALSE))
  )
}

shinyApp(ui, server)
