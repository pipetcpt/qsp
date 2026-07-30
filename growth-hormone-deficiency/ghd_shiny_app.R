##  Paediatric Growth Hormone Deficiency (GHD) — QSP Shiny Dashboard
##  ============================================================================
##  Interactive front end for ghd_mrgsolve_model.R (42 ODEs).
##
##  Run:
##    setwd("<this directory>")
##    shiny::runApp("ghd_shiny_app.R")
##
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##
##  Eight tabs:
##    1. 환자 프로파일      Patient & aetiology — where is the lesion?
##    2. GH PK             daily vs weekly exposure, free vs total, carrier pool
##    3. 신호전달·IGF 축     GHR / STAT5b / SOCS2 / IGF-1 / IGFBP-3 / ALS / free IGF-1
##    4. 성장판             resting-zone reserve, zones, senescence, fusion
##    5. 임상 엔드포인트     height, height SDS, height velocity, bone age, near-adult height
##    6. 사춘기·성호르몬     puberty, testosterone, oestradiol, GnRHa / aromatase inhibitor
##    7. 대사·체성분·안전성  fat/lean mass, HOMA-IR, IGF-1 SDS ceiling, free T4, ICP, ADA
##    8. 시나리오 비교       side-by-side comparison + summary table
##
##  Disclaimer: research / education / hypothesis generation only. Not for
##  clinical use. See the model file header for calibration and limitations.
##  ============================================================================

# NOTE: load mrgsolve BEFORE shiny would be wrong the other way round too —
# mrgsolve exports its own req(), which masks shiny::req(). Any req() used for
# input validation in this file is therefore written as shiny::req().
library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------
MOD <- mread_cache("ghd_mrgsolve_model.R", quiet = TRUE)

LAGH <- list(
  "lonapegsomatropin (weekly)" =
    list(KA_W = 0.45, KREL = 0.55, CLP_K = 0.08, POT_W = 1.00, FBIO_W = 0.72,
         mgkgwk = 0.24),
  "somatrogon (weekly)" =
    list(KA_W = 0.60, KREL = 0.58, CLP_K = 0.10, POT_W = 0.33, FBIO_W = 0.72,
         mgkgwk = 0.66),
  "somapacitan (weekly)" =
    list(KA_W = 0.50, KREL = 0.35, CLP_K = 0.05, POT_W = 1.45, FBIO_W = 0.72,
         mgkgwk = 0.16)
)

AETIOLOGY <- list(
  "Severe congenital isolated GHD (GH1/GHRHR)" =
    list(PITMAX = 0.06, GHRCAP = 1, STAT5CAP = 1, ALSCAP = 1, THYRCAP = 1, HYPOGON = 0, XRT_ON = 0),
  "Partial / idiopathic isolated GHD" =
    list(PITMAX = 0.35, GHRCAP = 1, STAT5CAP = 1, ALSCAP = 1, THYRCAP = 1, HYPOGON = 0, XRT_ON = 0),
  "Combined pituitary hormone deficiency (POU1F1/PROP1)" =
    list(PITMAX = 0.08, GHRCAP = 1, STAT5CAP = 1, ALSCAP = 1, THYRCAP = 0.45, HYPOGON = 1, XRT_ON = 0),
  "Post-irradiation, progressive somatotroph loss" =
    list(PITMAX = 0.60, GHRCAP = 1, STAT5CAP = 1, ALSCAP = 1, THYRCAP = 0.8, HYPOGON = 0, XRT_ON = 1),
  "GH insensitivity - Laron (GHR loss of function)" =
    list(PITMAX = 1.0, GHRCAP = 0.05, STAT5CAP = 1, ALSCAP = 1, THYRCAP = 1, HYPOGON = 0, XRT_ON = 0),
  "GH insensitivity - STAT5B defect" =
    list(PITMAX = 1.0, GHRCAP = 1, STAT5CAP = 0.12, ALSCAP = 1, THYRCAP = 1, HYPOGON = 0, XRT_ON = 0),
  "IGFALS (acid-labile subunit) deficiency" =
    list(PITMAX = 1.0, GHRCAP = 1, STAT5CAP = 1, ALSCAP = 0.05, THYRCAP = 1, HYPOGON = 0, XRT_ON = 0),
  "Healthy reference child (no deficiency)" =
    list(PITMAX = 1.0, GHRCAP = 1, STAT5CAP = 1, ALSCAP = 1, THYRCAP = 1, HYPOGON = 0, XRT_ON = 0)
)

REGIMENS <- c("No GH treatment (natural history)",
              "Daily somatropin",
              names(LAGH),
              "Mecasermin (rhIGF-1) BID")

THEME <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        panel.grid.minor = element_blank())

# ---------------------------------------------------------------------------
# Simulation helpers
# ---------------------------------------------------------------------------
ev_daily <- function(mgkgd, WT0, days, adherence = 1, grow = 0.10) {
  d <- seq(0, days - 1, by = 1)
  a <- mgkgd * WT0 * (1 + grow)^floor(d / 365)
  keep <- if (adherence >= 1) rep(TRUE, length(d)) else
            (seq_along(d) %% max(1, round(1 / adherence))) == 0
  ev(time = d[keep], amt = a[keep], cmt = "DEPD")
}
ev_weekly <- function(mgkgwk, WT0, days, grow = 0.10) {
  d <- seq(0, days - 1, by = 7)
  ev(time = d, amt = mgkgwk * WT0 * (1 + grow)^floor(d / 365), cmt = "DEPW")
}
ev_meca <- function(ugkg, WT0, days, grow = 0.10) {
  d <- seq(0, days - 0.5, by = 0.5)
  ev(time = d, amt = ugkg * WT0 * (1 + grow)^floor(d / 365) / 1000, cmt = "DEPM")
}

#' Build the parameter list and event table for one arm, then integrate.
simulate_arm <- function(aetiology, regimen, dose, age0, ba0, mphsds, sex,
                         years, adherence, wt0, leuprate, airate, gcdose,
                         lt4, delta = 1) {
  p <- AETIOLOGY[[aetiology]]
  p$AGE0 <- age0; p$BA0 <- ba0; p$MPHSDS <- mphsds; p$SEX <- sex
  p$LEUPRATE <- leuprate; p$AIRATE <- airate; p$GCDOSE <- gcdose; p$LT4 <- lt4
  days <- round(years * 365)

  if (regimen %in% names(LAGH)) {
    prod <- LAGH[[regimen]]
    p[c("KA_W", "KREL", "CLP_K", "POT_W", "FBIO_W")] <-
      prod[c("KA_W", "KREL", "CLP_K", "POT_W", "FBIO_W")]
    e <- ev_weekly(dose, wt0, days)
  } else if (regimen == "Daily somatropin") {
    e <- ev_daily(dose, wt0, days, adherence)
  } else if (regimen == "Mecasermin (rhIGF-1) BID") {
    e <- ev_meca(dose, wt0, days)
  } else {
    e <- ev(amt = 0, cmt = "DEPD")
  }

  mrgsim(param(MOD, p), e, end = days, delta = delta, hmax = 0.25) %>%
    as_tibble() %>%
    filter(!duplicated(time, fromLast = TRUE))
}

#' Annualised height velocity from the height trajectory (cm per calendar year).
annual_hv <- function(d) {
  d %>% mutate(yr = floor(time / 365) + 1) %>% group_by(yr) %>%
    summarise(HV = max(Height_cm) - min(Height_cm),
              IGF1_SDS = mean(IGF1_SDS),
              HtSDS_end = last(Height_SDS), .groups = "drop")
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Paediatric Growth Hormone Deficiency — QSP Dashboard"),
  p(tags$em(paste("42-ODE mrgsolve model. Height velocity is computed as the",
                  "PRODUCT of hypothalamic drive x somatotroph mass x GH PK x",
                  "GHR/STAT5b transduction x IGF-1 output and binding state x",
                  "growth-plate responsiveness x REMAINING proliferative",
                  "reserve. Research/education only — not for clinical use."))),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      selectInput("aetiology", "Aetiology / lesion site",
                  choices = names(AETIOLOGY), selected = names(AETIOLOGY)[1]),
      radioButtons("sex", "Sex", choices = c("Boy" = 1, "Girl" = 0),
                   selected = 1, inline = TRUE),
      sliderInput("age0", "Age at start (yr)", 2, 15, 5, step = 0.5),
      sliderInput("ba0", "Bone age at start (yr)", 1, 15, 4, step = 0.5),
      sliderInput("mphsds", "Mid-parental (target) height SDS", -3, 2, -0.5, step = 0.1),
      sliderInput("wt0", "Weight at start (kg)", 8, 70, 15.5, step = 0.5),
      hr(),
      h4("치료 (Treatment)"),
      selectInput("regimen", "Regimen", choices = REGIMENS,
                  selected = "Daily somatropin"),
      uiOutput("dose_ui"),
      sliderInput("adherence", "Adherence (daily product only)", 0.2, 1, 1, step = 0.05),
      sliderInput("years", "Horizon (yr)", 1, 8, 4, step = 1),
      hr(),
      h4("병용·배경 (Adjuvants & background)"),
      checkboxInput("gnrha", "GnRH analogue depot (delay puberty)", FALSE),
      checkboxInput("ai", "Aromatase inhibitor (anastrozole 1 mg/day)", FALSE),
      sliderInput("gcdose", "Glucocorticoid (pred-equiv mg/m2/day)", 0, 20, 0, step = 1),
      sliderInput("lt4", "Levothyroxine replacement adequacy", 0.2, 1, 1, step = 0.05),
      hr(),
      h4("비교 대상 (Comparator)"),
      selectInput("cmp_regimen", "Second arm", choices = REGIMENS,
                  selected = "No GH treatment (natural history)"),
      helpText("The comparator uses the same patient and the label dose of the",
               "selected product.")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. 환자 프로파일",
                 br(), fluidRow(column(6, DTOutput("tbl_profile")),
                                column(6, plotOutput("plt_gain", height = 320))),
                 br(), verbatimTextOutput("txt_profile")),
        tabPanel("2. GH PK",
                 br(), plotOutput("plt_pk_long", height = 260),
                 plotOutput("plt_pk_zoom", height = 300),
                 helpText("Top: whole horizon. Bottom: a single 3-week window at",
                          "steady state, showing the shape a daily injection",
                          "makes versus a weekly depot, and the carrier pool that",
                          "the weekly products release from.")),
        tabPanel("3. 신호전달·IGF 축",
                 br(), plotOutput("plt_signal", height = 300),
                 plotOutput("plt_igf", height = 320),
                 helpText("Total IGF-1 falls further than free IGF-1 in GHD",
                          "because the IGFBP-3/ALS carrier pool falls too: the",
                          "free fraction rises and clearance of what is made",
                          "accelerates. This is why growth is driven in the model",
                          "mainly by the LOCAL, GH-dependent IGF-1 arm.")),
        tabPanel("4. 성장판",
                 br(), plotOutput("plt_plate", height = 340),
                 plotOutput("plt_reserve", height = 300),
                 helpText("Reserve (RZ) multiplies the proliferation rate AND is",
                          "spent in proportion to it. An untreated child saves",
                          "reserve; starting rhGH spends it. That is the whole",
                          "mechanism behind catch-up growth and behind the fall",
                          "in height velocity from year 1 to year 3 on an",
                          "unchanged mg/kg dose.")),
        tabPanel("5. 임상 엔드포인트",
                 br(), plotOutput("plt_height", height = 340),
                 fluidRow(column(6, plotOutput("plt_hv", height = 300)),
                          column(6, plotOutput("plt_ba", height = 300))),
                 br(), DTOutput("tbl_endpoints")),
        tabPanel("6. 사춘기·성호르몬",
                 br(), plotOutput("plt_pub", height = 340),
                 plotOutput("plt_fusion", height = 300),
                 helpText("Oestrogen enters the model twice with opposite sign:",
                          "it amplifies GH pulse amplitude (the spurt, cm/yr up)",
                          "and it accelerates reserve consumption and epiphyseal",
                          "fusion (growing years down). Whether a GnRH analogue",
                          "or an aromatase inhibitor helps is therefore an",
                          "emergent trade-off, not a coded height gain.")),
        tabPanel("7. 대사·체성분·안전성",
                 br(), plotOutput("plt_body", height = 300),
                 plotOutput("plt_safety", height = 340),
                 helpText("The dashed line on IGF-1 SDS is the +2 SDS ceiling at",
                          "which guidelines require a dose reduction. Free T4",
                          "falls after starting rhGH (increased T4 to T3",
                          "conversion) — the fall unmasks central hypothyroidism",
                          "rather than causing it, so the growth-relevant",
                          "thyroid tone is plotted alongside it.")),
        tabPanel("8. 시나리오 비교",
                 br(), plotOutput("plt_cmp", height = 420),
                 br(), DTOutput("tbl_cmp"),
                 br(), downloadButton("dl", "Download both arms (CSV)"))
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  # dose control adapts to the selected product
  output$dose_ui <- renderUI({
    r <- input$regimen
    if (r == "Daily somatropin")
      sliderInput("dose", "Dose (mg/kg/day)", 0.010, 0.070, 0.034, step = 0.001)
    else if (r %in% names(LAGH))
      sliderInput("dose", "Dose (mg/kg/week)", 0.05, 1.20,
                  LAGH[[r]]$mgkgwk, step = 0.01)
    else if (r == "Mecasermin (rhIGF-1) BID")
      sliderInput("dose", "Dose (ug/kg per dose, BID)", 20, 160, 100, step = 10)
    else helpText("No dose — natural history.")
  })

  label_dose <- function(r) {
    if (r == "Daily somatropin") 0.034
    else if (r %in% names(LAGH)) LAGH[[r]]$mgkgwk
    else if (r == "Mecasermin (rhIGF-1) BID") 100
    else 0
  }

  common <- reactive({
    list(aetiology = input$aetiology,
         age0 = input$age0, ba0 = input$ba0, mphsds = input$mphsds,
         sex = as.numeric(input$sex), years = input$years,
         wt0 = input$wt0,
         leuprate = if (isTRUE(input$gnrha)) 0.134 else 0,
         airate   = if (isTRUE(input$ai)) 1.0 else 0,
         gcdose = input$gcdose, lt4 = input$lt4)
  })

  simA <- reactive({
    shiny::req(input$regimen)
    d <- input$dose
    if (is.null(d)) d <- label_dose(input$regimen)
    a <- common()
    do.call(simulate_arm, c(list(regimen = input$regimen, dose = d,
                                 adherence = input$adherence), a)) %>%
      mutate(arm = paste0("A: ", input$regimen))
  })

  simB <- reactive({
    a <- common()
    do.call(simulate_arm, c(list(regimen = input$cmp_regimen,
                                 dose = label_dose(input$cmp_regimen),
                                 adherence = 1), a)) %>%
      mutate(arm = paste0("B: ", input$cmp_regimen))
  })

  both <- reactive(bind_rows(simA(), simB()))

  # high-resolution window for the PK tab (last 3 weeks of the horizon)
  simA_fine <- reactive({
    d <- input$dose; if (is.null(d)) d <- label_dose(input$regimen)
    a <- common()
    do.call(simulate_arm, c(list(regimen = input$regimen, dose = d,
                                 adherence = input$adherence, delta = 0.02), a)) %>%
      mutate(arm = paste0("A: ", input$regimen))
  })
  simB_fine <- reactive({
    a <- common()
    do.call(simulate_arm, c(list(regimen = input$cmp_regimen,
                                 dose = label_dose(input$cmp_regimen),
                                 adherence = 1, delta = 0.02), a)) %>%
      mutate(arm = paste0("B: ", input$cmp_regimen))
  })

  gg_line <- function(d, y, ylab, hline = NULL) {
    p <- ggplot(d, aes(Age_yr, .data[[y]], colour = arm)) +
      geom_line(linewidth = 0.7) + labs(x = "Age (yr)", y = ylab) + THEME
    if (!is.null(hline))
      p <- p + geom_hline(yintercept = hline, linetype = 2, colour = "grey30")
    p
  }

  # ---- tab 1: patient profile ---------------------------------------------
  output$tbl_profile <- renderDT({
    p <- AETIOLOGY[[input$aetiology]]
    tibble(
      Parameter = c("Somatotroph capacity (PITMAX)", "GH-receptor capacity (GHRCAP)",
                    "STAT5b capacity (STAT5CAP)", "ALS capacity (ALSCAP)",
                    "Thyrotroph capacity (THYRCAP)", "Hypogonadotrophic (HYPOGON)",
                    "Progressive injury (XRT_ON)"),
      Value = c(p$PITMAX, p$GHRCAP, p$STAT5CAP, p$ALSCAP, p$THYRCAP,
                p$HYPOGON, p$XRT_ON)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t", paging = FALSE))
  })

  output$plt_gain <- renderPlot({
    d0 <- filter(simA(), time == 0)
    tibble(
      gain = c("GH exposure", "STAT5b transduction", "IGF-1 output",
               "Carrier (ternary) pool", "Growth-plate reserve"),
      value = c(d0$GH_active / 2.2, d0$pSTAT5b, d0$IGF1_total / 97,
                d0$Ternary_frac / 0.993, d0$Reserve_RZ / 0.56)
    ) %>%
      mutate(gain = factor(gain, levels = rev(gain))) %>%
      ggplot(aes(value, gain)) +
      geom_col(fill = "#4878a8", width = 0.62) +
      geom_vline(xintercept = 1, linetype = 2) +
      labs(x = "Fraction of the healthy reference child (dashed = 1.0)", y = NULL,
           title = "Where is the chain broken? (at time 0)") + THEME
  })

  output$txt_profile <- renderPrint({
    hv <- annual_hv(simA())
    cat("Arm A —", input$regimen, "\n")
    cat("Annualised height velocity (cm/yr):",
        paste(sprintf("yr%d %.2f", hv$yr, hv$HV), collapse = "  "), "\n")
    l <- tail(simA(), 1); f <- head(simA(), 1)
    cat(sprintf("Height SDS %.2f -> %.2f   |  bone age %.1f yr (BA/CA %.2f)\n",
                f$Height_SDS, l$Height_SDS, l$BoneAge_yr, l$BA_CA_ratio))
    cat(sprintf("IGF-1 SDS %.2f -> %.2f   |  reserve remaining %.2f  |  fusion %.2f\n",
                f$IGF1_SDS, l$IGF1_SDS, l$Reserve_RZ, l$Fusion_prog))
  })

  # ---- tab 2: PK ----------------------------------------------------------
  output$plt_pk_long <- renderPlot({
    both() %>% select(Age_yr, arm, `Total GH` = GH_total, `Active GH` = GH_active) %>%
      pivot_longer(-c(Age_yr, arm)) %>%
      ggplot(aes(Age_yr, value, colour = arm)) + geom_line(linewidth = 0.5) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Age (yr)", y = "GH (ug/L)",
           title = "GH exposure over the whole horizon (daily-resolution samples)") +
      THEME
  })

  output$plt_pk_zoom <- renderPlot({
    d <- bind_rows(simA_fine(), simB_fine())
    tmax <- max(d$time); w <- c(max(0, tmax - 21), tmax)
    d %>% filter(time >= w[1], time <= w[2]) %>%
      select(time, arm, `Total GH (ug/L)` = GH_total,
             `Carrier / prodrug (mg)` = GH_carrier,
             `Total IGF-1 (ng/mL)` = IGF1_total) %>%
      pivot_longer(-c(time, arm)) %>%
      ggplot(aes(time - w[1], value, colour = arm)) + geom_line(linewidth = 0.6) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Day within the final 3 weeks", y = NULL,
           title = "Steady-state shape: daily peak-and-trough vs weekly depot release") +
      THEME
  })

  # ---- tab 3: signalling and IGF axis ------------------------------------
  output$plt_signal <- renderPlot({
    both() %>% select(Age_yr, arm, `GHR density` = GHR_density,
                      `STAT5b signal` = pSTAT5b, `SOCS2 brake` = SOCS2) %>%
      pivot_longer(-c(Age_yr, arm)) %>%
      ggplot(aes(Age_yr, value, colour = arm)) + geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Age (yr)", y = "normalised (1 = healthy)",
           title = "Transduction: receptor, STAT5b signal, and the SOCS2 brake") + THEME
  })

  output$plt_igf <- renderPlot({
    both() %>% select(Age_yr, arm, `Total IGF-1 (ng/mL)` = IGF1_total,
                      `Free IGF-1 (ng/mL)` = IGF1_free,
                      `IGFBP-3 (mg/L)` = IGFBP3, `ALS (norm.)` = ALS_level,
                      `Ternary-complexed fraction` = Ternary_frac,
                      `IGF-1 SDS` = IGF1_SDS) %>%
      pivot_longer(-c(Age_yr, arm)) %>%
      ggplot(aes(Age_yr, value, colour = arm)) + geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Age (yr)", y = NULL, title = "IGF system and its binding state") + THEME
  })

  # ---- tab 4: growth plate ----------------------------------------------
  output$plt_plate <- renderPlot({
    both() %>% select(Age_yr, arm, `Proliferative zone` = Prolif_zone,
                      `Hypertrophic zone` = Hyper_zone,
                      `IGF-1 tone (0-1)` = IGF1_tone,
                      `Direct GH effect` = GH_drive_dir) %>%
      pivot_longer(-c(Age_yr, arm)) %>%
      ggplot(aes(Age_yr, value, colour = arm)) + geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Age (yr)", y = NULL, title = "Growth-plate zones and their drives") + THEME
  })

  output$plt_reserve <- renderPlot({
    both() %>% select(Age_yr, arm, `Resting-zone reserve (RZ)` = Reserve_RZ,
                      `Fusion program (0-1)` = Fusion_prog,
                      `Height velocity (cm/yr)` = HV_cm_yr) %>%
      pivot_longer(-c(Age_yr, arm)) %>%
      ggplot(aes(Age_yr, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Age (yr)", y = NULL,
           title = "Reserve is a consumable: spending it is what buys the cm") + THEME
  })

  # ---- tab 5: clinical endpoints ----------------------------------------
  output$plt_height <- renderPlot({
    d <- both()
    ggplot(d, aes(Age_yr, Height_SDS, colour = arm)) +
      geom_line(linewidth = 0.9) +
      geom_hline(aes(yintercept = Target_SDS), linetype = 3, colour = "grey20") +
      geom_hline(yintercept = -2, linetype = 2, colour = "firebrick") +
      labs(x = "Age (yr)", y = "Height SDS",
           title = "Height SDS (dotted = mid-parental target, dashed = -2 SDS)") +
      THEME
  })
  output$plt_hv <- renderPlot(gg_line(both(), "HV_cm_yr", "Height velocity (cm/yr)"))
  output$plt_ba <- renderPlot({
    both() %>% select(Age_yr, arm, `Bone age (yr)` = BoneAge_yr,
                      `Bone age / chronological age` = BA_CA_ratio) %>%
      pivot_longer(-c(Age_yr, arm)) %>%
      ggplot(aes(Age_yr, value, colour = arm)) + geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Age (yr)", y = NULL,
           title = "Skeletal maturation — advancing it too fast steals final height") +
      THEME
  })
  output$tbl_endpoints <- renderDT({
    bind_rows(mutate(annual_hv(simA()), arm = "A"),
              mutate(annual_hv(simB()), arm = "B")) %>%
      transmute(Arm = arm, Year = yr,
                `HV (cm/yr)` = round(HV, 2),
                `IGF-1 SDS (mean)` = round(IGF1_SDS, 2),
                `Height SDS (end of yr)` = round(HtSDS_end, 2)) %>%
      datatable(rownames = FALSE, options = list(dom = "t", paging = FALSE))
  })

  # ---- tab 6: puberty ---------------------------------------------------
  output$plt_pub <- renderPlot({
    both() %>% select(Age_yr, arm, `GnRH pulse drive` = GnRH_drive,
                      `Testosterone (ng/dL)` = Testosterone,
                      `Oestradiol (pg/mL)` = Estradiol) %>%
      pivot_longer(-c(Age_yr, arm)) %>%
      ggplot(aes(Age_yr, value, colour = arm)) + geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Age (yr)", y = NULL, title = "Pubertal axis and sex steroids") + THEME
  })
  output$plt_fusion <- renderPlot({
    both() %>% select(Age_yr, arm, `Fusion program` = Fusion_prog,
                      `Reserve (RZ)` = Reserve_RZ,
                      `Height velocity (cm/yr)` = HV_cm_yr,
                      `BMD Z-score` = BMD_Z) %>%
      pivot_longer(-c(Age_yr, arm)) %>%
      ggplot(aes(Age_yr, value, colour = arm)) + geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Age (yr)", y = NULL,
           title = "The oestrogen trade-off: cm/yr against years-of-growth (and bone)") +
      THEME
  })

  # ---- tab 7: metabolism, body composition, safety ----------------------
  output$plt_body <- renderPlot({
    both() %>% select(Age_yr, arm, `Fat mass (kg)` = Fat_mass_kg,
                      `Body fat (%)` = Fat_pct,
                      `Lean mass (kg)` = Lean_mass_kg,
                      `Weight (kg)` = Weight_kg) %>%
      pivot_longer(-c(Age_yr, arm)) %>%
      ggplot(aes(Age_yr, value, colour = arm)) + geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", nrow = 1) +
      labs(x = "Age (yr)", y = NULL, title = "Body composition") + THEME
  })
  output$plt_safety <- renderPlot({
    d <- both() %>%
      select(Age_yr, arm, `IGF-1 SDS (cap +2)` = IGF1_SDS, `HOMA-IR` = HOMA_IR,
             `Free T4 (ng/dL)` = FreeT4, `Thyroid tone (growth-relevant)` = Thyroid_tone,
             `LDL-C (mg/dL)` = LDL_C, `Intracranial-pressure index` = ICP_index,
             `ADA titre` = ADA_titre, `BMD Z-score` = BMD_Z) %>%
      pivot_longer(-c(Age_yr, arm))
    caps <- tibble(name = "IGF-1 SDS (cap +2)", y = 2)
    ggplot(d, aes(Age_yr, value, colour = arm)) + geom_line(linewidth = 0.7) +
      geom_hline(data = caps, aes(yintercept = y), linetype = 2,
                 colour = "firebrick", inherit.aes = FALSE) +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "Age (yr)", y = NULL, title = "Safety and monitoring read-outs") + THEME
  })

  # ---- tab 8: comparison ------------------------------------------------
  output$plt_cmp <- renderPlot({
    both() %>% select(Age_yr, arm, `Height SDS` = Height_SDS,
                      `Height velocity (cm/yr)` = HV_cm_yr,
                      `IGF-1 SDS` = IGF1_SDS, `Reserve (RZ)` = Reserve_RZ,
                      `Bone age (yr)` = BoneAge_yr, `Fat (%)` = Fat_pct) %>%
      pivot_longer(-c(Age_yr, arm)) %>%
      ggplot(aes(Age_yr, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Age (yr)", y = NULL, title = "Arm A versus arm B") + THEME
  })

  output$tbl_cmp <- renderDT({
    summ <- function(d, lab) {
      f <- head(d, 1); l <- tail(d, 1); hv <- annual_hv(d)
      tibble(Arm = lab,
             `Yr-1 HV (cm/yr)` = round(hv$HV[1], 2),
             `Final-yr HV` = round(tail(hv$HV, 1), 2),
             `Height SDS start` = round(f$Height_SDS, 2),
             `Height SDS end` = round(l$Height_SDS, 2),
             `delta Height SDS` = round(l$Height_SDS - f$Height_SDS, 2),
             `Final height (cm)` = round(l$Height_cm, 1),
             `IGF-1 SDS end` = round(l$IGF1_SDS, 2),
             `Peak IGF-1 SDS` = round(max(d$IGF1_SDS), 2),
             `Bone age end` = round(l$BoneAge_yr, 1),
             `Reserve left` = round(l$Reserve_RZ, 3),
             `Fusion` = round(l$Fusion_prog, 2),
             `BMD Z end` = round(l$BMD_Z, 2),
             `Fat % end` = round(l$Fat_pct, 1),
             `HOMA-IR end` = round(l$HOMA_IR, 2))
    }
    bind_rows(summ(simA(), unique(simA()$arm)),
              summ(simB(), unique(simB()$arm))) %>%
      datatable(rownames = FALSE,
                options = list(dom = "t", paging = FALSE, scrollX = TRUE))
  })

  output$dl <- downloadHandler(
    filename = function() "ghd_qsp_simulation.csv",
    content  = function(f) write.csv(both(), f, row.names = FALSE)
  )
}

shinyApp(ui, server)
