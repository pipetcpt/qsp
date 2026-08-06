## =============================================================================
## HDFN QSP dashboard -- Shiny front end for hdfn_mrgsolve_model.R
##
## The app is organised around the identity the model is built on:
##
##        destruction  =  A (delivered antibody)
##                      x f_ag (antigen-positive fraction)
##                      x M (clearance capacity)
##
## so that a user can see WHICH FACTOR a given intervention moves.  Tab 4 shows
## the three factors as separate traces; every therapeutic control in the
## sidebar is labelled with the factor it acts on.
##
## Twelve tabs:
##   1  Patient & protocol      the case, and what is being done about it
##   2  Maternal PK/PD          IgG, anti-D, nipocalimab, IVIG, plasmapheresis
##   3  The placental conveyor  capacity vs gestation, free FcRn, F/M ratio
##   4  The three factors       A x f_ag x M and the flux they multiply to
##   5  Fetal haematology       Hb, MoM, red cell pools, erythropoiesis
##   6  Surveillance            MCA-PSV MoM with the derived threshold
##   7  Transfusion course      volumes, intervals, the predicted decline
##   8  Hydrops                 Starling terms, ascites, cardiac reserve
##   9  Neonatal bilirubin      the clock that changes at birth
##  10  Neonatal haematology    late anaemia, top-ups, iron
##  11  Prophylaxis calculator  the stoichiometric race and the 72-hour window
##  12  Scenario comparison     the sixteen scenarios side by side
##
## NOTE: no R toolchain was available in the environment where this model was
## built, so this file has not been executed.  It mirrors the equations and the
## controller of hdfn_python_reference.py, which has.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

source("hdfn_mrgsolve_model.R")   # defines hdfn, simulate_hdfn, hdfn_scenarios

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

## ---------------------------------------------------------------- UI ---------
ui <- fluidPage(
  titlePanel("Hemolytic Disease of the Fetus and Newborn — QSP dashboard"),
  p(tags$b("destruction = A × f_ag × M."),
    "Every control below is labelled with the factor it moves. Nipocalimab, ",
    "IVIG and plasmapheresis act on A (what reaches the fetus); intrauterine ",
    "transfusion acts on f_ag (the antigen-positive fraction — a transfusion ",
    "dilutes the substrate, it does not only add haemoglobin); high-dose IVIG ",
    "also acts on M (FcγR blockade)."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("The case"),
      sliderInput("anti_d", "Maternal anti-D quantitation (IU/mL) — sets A",
                  min = 0, max = 300, value = 15, step = 1),
      helpText("UK risk bands: <4 low, 4–15 moderate, >15 high."),
      selectInput("antigen", "Antibody specificity",
                  c("anti-D (R1r, 10 000 sites)" = "d1",
                    "anti-D (R2R2, 30 000 sites)" = "d2",
                    "anti-c" = "c",
                    "anti-K (kills progenitors)" = "k",
                    "ABO (soluble antigen sink)" = "abo")),
      sliderInput("ga0", "Gestational age at entry (wk)", 11, 24, 13, 1),
      sliderInput("ga_del", "Planned delivery (wk)", 30, 39, 36, 0.5),

      h4("Surveillance"),
      selectInput("protocol", "Protocol",
                  c("weekly MCA-PSV, IUT if ≥ cut-off" = "mca",
                    "fixed 2-weekly IUT" = "fixed",
                    "observation only" = "none")),
      sliderInput("mca_cut", "MCA-PSV threshold (MoM)", 1.2, 1.8, 1.5, 0.05),
      helpText(tags$b("1.50 MoM is not empirical:"),
               "if cerebral oxygen delivery is defended then PSV MoM = 1/Hb MoM,",
               "so 1.50 is the arithmetic image of Hb 0.667 MoM — the",
               "definition of moderate anaemia."),
      sliderInput("iut_target", "IUT target haematocrit", 0.35, 0.55, 0.45, 0.01),

      h4("Therapy on factor A"),
      sliderInput("nip_dose", "Nipocalimab (mg/kg/week)", 0, 60, 0, 5),
      sliderInput("nip_start", "Nipocalimab start (wk)", 12, 30, 14, 1),
      sliderInput("nip_pen", "Placental FcRn penetration (fraction of plasma)",
                  0, 1, 0.15, 0.05),
      helpText("The single most consequential structural assumption in the",
               "model: at 1.0 no fetus would ever need a transfusion, which",
               "the UNITY trial contradicts."),
      sliderInput("ivig", "Maternal IVIG (g/kg/week)", 0, 2, 0, 0.5),
      checkboxInput("pheresis", "Plasmapheresis ×3 before IVIG", FALSE),

      h4("Neonatal care"),
      sliderInput("photo", "Phototherapy irradiance (µW/cm²/nm)", 0, 50, 30, 5),
      helpText("Neonatal IVIG is NOT separately implemented: the model has a",
               "fetal FcγR-blockade term driven by MATERNAL IVIG, and there is",
               "no neonatal dosing route for it. Tab 9 states this."),
      numericInput("seed", "Surveillance seed", 7, 1, 999, 1),
      actionButton("go", "Run", class = "btn-primary")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · Patient & protocol",
                 h4("Case summary"), tableOutput("tbl_case"),
                 h4("What happened"), verbatimTextOutput("txt_summary"),
                 plotOutput("p_overview", height = "420px")),
        tabPanel("2 · Maternal PK/PD",
                 plotOutput("p_mat", height = "560px"),
                 helpText("Nipocalimab blocks the FcRn that recycles IgG, so",
                          "maternal IgG half-life falls from 21 d to about 6 d",
                          "and the steady state to ~28% of baseline.",
                          "Plasmapheresis removes the STOCK; the long-lived",
                          "plasma cells that make it are untouched, which is",
                          "the rebound.")),
        tabPanel("3 · The conveyor",
                 plotOutput("p_conveyor", height = "520px"),
                 helpText("Transfer capacity grows e^0.268 per week — a factor",
                          "of about 210 from 16 to 36 weeks, which after",
                          "dilution into a growing fetus leaves the",
                          "fetal:maternal antibody ratio ~17-fold higher at",
                          "term than at 19.5 weeks. The same maternal titre",
                          "is therefore a different disease at each",
                          "gestational age, and 'early-onset severe HDFN' is",
                          "a titre high enough to matter while the conveyor",
                          "is still weak.")),
        tabPanel("4 · A × f_ag × M",
                 plotOutput("p_factors", height = "560px"),
                 helpText("Watch what a transfusion does: A is unchanged",
                          "(packed cells carry almost no plasma), f_ag",
                          "collapses to about a third, and the product falls",
                          "with it. That is why the effect of a transfusion",
                          "outlasts the haemoglobin it delivered.")),
        tabPanel("5 · Fetal haematology",
                 plotOutput("p_fetal", height = "560px"),
                 tableOutput("tbl_pools")),
        tabPanel("6 · Surveillance",
                 plotOutput("p_psv", height = "460px"),
                 h4("The threshold, inverted"), tableOutput("tbl_psv")),
        tabPanel("7 · Transfusion course",
                 tableOutput("tbl_iut"),
                 plotOutput("p_decline", height = "440px"),
                 helpText("Reported decline between the first two procedures:",
                          "0.40 g/dL/day (SD 0.25). The model predicts the",
                          "ORDER (each interval flatter than the last, because",
                          "f_ag halves again each time) and under-predicts the",
                          "LEVEL by about a third — closing the gap requires",
                          "donor cells to survive 30–45 d in the fetus rather",
                          "than 70.")),
        tabPanel("8 · Hydrops",
                 plotOutput("p_hydrops", height = "560px"),
                 helpText("Lymphatic capacity is DEFINED as baseline",
                          "filtration, so a healthy fetus cannot drift into",
                          "ascites and hydrops has to come from a named term:",
                          "hypoalbuminaemia from hepatic extramedullary",
                          "erythropoiesis, venous pressure, or permeability.")),
        tabPanel("9 · Neonatal bilirubin",
                 plotOutput("p_bili", height = "520px"),
                 helpText("Neonatal IVIG has no dosing route in this model.",
                          "Nothing about the antibody changes at delivery.",
                          "What changes is the clearance route: the placenta",
                          "(the mother's liver, k = 14/d) is replaced by a",
                          "liver with 6% of adult UGT1A1. Identical",
                          "haemolysis, two different diseases.")),
        tabPanel("10 · Neonatal haematology",
                 plotOutput("p_neo", height = "520px"),
                 helpText("The nadir arrives at five to six weeks, not in the",
                          "first week: maternal IgG persists with a 3-week",
                          "half-life while a marrow suppressed by months of",
                          "transfusion has nothing in reserve.")),
        tabPanel("11 · Prophylaxis",
                 fluidRow(
                   column(4,
                          numericInput("fmh", "Fetomaternal haemorrhage (mL fetal whole blood)",
                                       30, 0, 300, 1),
                          numericInput("rhig", "Anti-D dose (µg)", 300, 0, 3000, 50),
                          sliderInput("delay", "Delay to administration (days)",
                                      0, 21, 1, 1)),
                   column(8, tableOutput("tbl_prophy"))),
                 plotOutput("p_prophy", height = "420px"),
                 helpText("The dose rule is stoichiometric: 20 µg per mL of",
                          "fetal red cells (100 IU/mL), so 300 µg covers 15 mL",
                          "of fetal red cells = 30 mL of fetal whole blood.",
                          "The 72-hour window is derived, not assumed: an",
                          "uncoated fetal cell in the maternal circulation is",
                          "just a red cell and survives for months, so a",
                          "three-day delay leaks only a few per cent of the",
                          "antigen exposure integral.")),
        tabPanel("12 · Scenarios",
                 actionButton("run_sc", "Run all sixteen scenarios"),
                 tableOutput("tbl_sc"))
      )
    )
  )
)

## ------------------------------------------------------------- SERVER --------
server <- function(input, output, session) {

  site_of <- function(a) switch(a, d1 = 1.0, d2 = 3.0, c = 0.6, k = 0.2,
                                abo = 0.08, 1.0)
  kell_of <- function(a) if (identical(a, "k")) 0.55 else 0.0

  sim <- eventReactive(input$go, {
    m <- param(hdfn,
               site_D = site_of(input$antigen),
               kell_kill = kell_of(input$antigen),
               nip_pl_pen = input$nip_pen)
    r <- simulate_hdfn(m, ga0 = input$ga0, ga_deliver = input$ga_del,
                       postnatal_days = 70, anti_d_iu = input$anti_d,
                       protocol = input$protocol,
                       nip_dose = input$nip_dose, nip_start = input$nip_start,
                       ivig_dose = input$ivig,
                       ivig_start = if (input$ivig > 0) input$ga0 else 99,
                       mca_cut = input$mca_cut, iut_target = input$iut_target,
                       photo_irr = input$photo, seed = input$seed)
    r$out$ga <- r$out$time / 7
    r$out$phase <- ifelse(r$out$ga < input$ga_del, "fetal", "neonatal")
    r
  }, ignoreNULL = FALSE)

  long <- function(d, cols) {
    d %>% select(ga, all_of(cols)) %>%
      pivot_longer(-ga, names_to = "variable", values_to = "value")
  }

  ## --- 1 -------------------------------------------------------------------
  output$tbl_case <- renderTable({
    data.frame(
      field = c("anti-D quantitation", "specificity / antigen density",
                "protocol", "MCA-PSV threshold", "nipocalimab", "IVIG",
                "planned delivery"),
      value = c(paste(input$anti_d, "IU/mL"),
                paste(input$antigen, "· site factor",
                      site_of(input$antigen)),
                input$protocol, paste(input$mca_cut, "MoM"),
                paste(input$nip_dose, "mg/kg/wk from", input$nip_start, "wk"),
                paste(input$ivig, "g/kg/wk"),
                paste(input$ga_del, "wk")))
  })

  output$txt_summary <- renderPrint({
    r <- sim(); o <- r$out
    cat("intrauterine transfusions :", r$n_iut, "\n")
    cat("first transfusion         :",
        ifelse(is.null(r$iut), "none", sprintf("%.1f wk", min(r$iut$ga))), "\n")
    cat("hydrops at any point      :", any(o$HYDROPS > 0.5), "\n")
    cat("cumulative survival       :", sprintf("%.3f", min(o$SURV)), "\n")
    nb <- o[o$phase == "neonatal", ]
    if (nrow(nb) > 0) {
      cat("Hb at birth               :", sprintf("%.1f g/dL", nb$HB[1]), "\n")
      cat("peak neonatal bilirubin   :", sprintf("%.1f mg/dL", max(nb$TSB)), "\n")
      cat("Hb nadir (postnatal)      :", sprintf("%.1f g/dL", min(nb$HB)), "\n")
    }
  })

  output$p_overview <- renderPlot({
    o <- sim()$out
    ggplot(long(o, c("HB", "HBMOM", "PSVMOM", "TSB")), aes(ga, value)) +
      geom_line(linewidth = 0.8, colour = "#2b5797") +
      geom_vline(xintercept = input$ga_del, linetype = 2, colour = "grey40") +
      facet_wrap(~variable, scales = "free_y") +
      labs(x = "gestational age (weeks)", y = NULL,
           caption = "dashed line = delivery") + THEME
  })

  ## --- 2 -------------------------------------------------------------------
  output$p_mat <- renderPlot({
    o <- sim()$out
    ggplot(long(o, c("CIG", "CA1", "CNIP")), aes(ga, value)) +
      geom_line(linewidth = 0.9, colour = "#7030a0") +
      facet_wrap(~variable, scales = "free_y", ncol = 1,
                 labeller = as_labeller(c(
                   CIG = "maternal total IgG (g/L)",
                   CA1 = "maternal anti-D IgG1 (IU/mL)",
                   CNIP = "nipocalimab (mg/L)"))) +
      labs(x = "gestational age (weeks)", y = NULL) + THEME
  })

  ## --- 3 -------------------------------------------------------------------
  output$p_conveyor <- renderPlot({
    o <- sim()$out
    p <- as.list(param(hdfn))
    d <- data.frame(ga = seq(14, 40, 0.2))
    d$capacity <- exp(p$g_pl * (d$ga - 20))
    ggplot(d, aes(ga, capacity)) +
      geom_line(linewidth = 1, colour = "#c9a227") +
      scale_y_log10() +
      labs(x = "gestational age (weeks)",
           y = "placental transfer capacity (relative to 20 wk, log scale)",
           title = sprintf("cap(GA) = v0 · e^(%.3f · (GA − 20))", p$g_pl),
           subtitle = "calibrated to a fetal:maternal IgG ratio of 0.075 at 19.5 wk and 1.25 at 39 wk") +
      THEME
  })

  ## --- 4 -------------------------------------------------------------------
  output$p_factors <- renderPlot({
    o <- sim()$out
    o$flux <- as.list(param(hdfn))$kops * o$AEFF * o$FAG *
      (o$HB * o$Vfp_o / 100)
    ggplot(long(o, c("AEFF", "FAG", "flux")), aes(ga, value)) +
      geom_line(linewidth = 0.9, colour = "#c0504d") +
      facet_wrap(~variable, scales = "free_y", ncol = 1,
                 labeller = as_labeller(c(
                   AEFF = "A — delivered antibody (IU/mL equivalents)",
                   FAG = "f_ag — antigen-positive fraction",
                   flux = "destruction flux (g Hb/day)"))) +
      labs(x = "gestational age (weeks)", y = NULL) + THEME
  })

  ## --- 5 -------------------------------------------------------------------
  output$p_fetal <- renderPlot({
    o <- sim()$out
    ggplot(long(o, c("HB", "HBMOM", "HCT", "AEFF")), aes(ga, value)) +
      geom_line(linewidth = 0.9, colour = "#1f7a1f") +
      geom_hline(data = data.frame(variable = "HBMOM", y = c(0.84, 0.65, 0.55)),
                 aes(yintercept = y), linetype = 3) +
      facet_wrap(~variable, scales = "free_y") +
      labs(x = "gestational age (weeks)", y = NULL,
           caption = "dotted lines: mild 0.84, moderate 0.65, severe 0.55 MoM") +
      THEME
  })
  output$tbl_pools <- renderTable({
    o <- sim()$out
    idx <- round(seq(1, nrow(o), length.out = 10))
    o[idx, c("ga", "HB", "HBMOM", "HCT", "FAG", "AEFF", "SURV")]
  }, digits = 3)

  ## --- 6 -------------------------------------------------------------------
  output$p_psv <- renderPlot({
    o <- sim()$out
    ggplot(o, aes(ga, PSVMOM)) +
      geom_line(linewidth = 0.9, colour = "#2b5797") +
      geom_hline(yintercept = input$mca_cut, linetype = 2, colour = "#c00000") +
      labs(x = "gestational age (weeks)", y = "MCA-PSV (MoM)",
           title = "PSV MoM = 1 / Hb MoM when cerebral oxygen delivery is defended") +
      THEME
  })
  output$tbl_psv <- renderTable({
    mom <- c(1.00, 0.84, 0.75, 0.667, 0.65, 0.55, 0.45)
    data.frame(`Hb MoM` = mom, `implied PSV MoM` = round(1 / mom, 3),
              interpretation = c("normal", "mild anaemia threshold", "",
                                 "PSV 1.50 MoM lands here",
                                 "moderate anaemia (definition)",
                                 "severe anaemia", "hydrops range"),
              check.names = FALSE)
  })

  ## --- 7 -------------------------------------------------------------------
  output$tbl_iut <- renderTable({
    r <- sim()
    if (is.null(r$iut)) return(data.frame(message = "no transfusion given"))
    d <- r$iut
    d$interval_days <- c(NA, diff(d$ga) * 7)
    d
  }, digits = 2)
  output$p_decline <- renderPlot({
    r <- sim(); o <- r$out
    ggplot(o[o$phase == "fetal", ], aes(ga, HB)) +
      geom_line(linewidth = 0.9) +
      { if (!is.null(r$iut)) geom_vline(xintercept = r$iut$ga, linetype = 2,
                                       colour = "#1f7a1f") } +
      labs(x = "gestational age (weeks)", y = "fetal haemoglobin (g/dL)",
           caption = "green lines = transfusions") + THEME
  })

  ## --- 8 -------------------------------------------------------------------
  output$p_hydrops <- renderPlot({
    o <- sim()$out
    ggplot(long(o, c("Asc", "Alb", "Card", "HYDROPS")), aes(ga, value)) +
      geom_line(linewidth = 0.9, colour = "#c00000") +
      facet_wrap(~variable, scales = "free_y",
                 labeller = as_labeller(c(
                   Asc = "ascites (mL)", Alb = "albumin (g/dL)",
                   Card = "cardiac decompensation index",
                   HYDROPS = "overt hydrops (0/1)"))) +
      labs(x = "gestational age (weeks)", y = NULL) + THEME
  })

  ## --- 9 -------------------------------------------------------------------
  output$p_bili <- renderPlot({
    o <- sim()$out
    o$pna <- (o$ga - input$ga_del) * 7
    nb <- o[o$pna >= -14, ]
    ggplot(long(nb %>% mutate(ga = pna), c("TSB", "BA", "Ugt", "OD450")),
           aes(ga, value)) +
      geom_line(linewidth = 0.9, colour = "#c9a227") +
      geom_vline(xintercept = 0, linetype = 2) +
      facet_wrap(~variable, scales = "free_y",
                 labeller = as_labeller(c(
                   TSB = "total serum bilirubin (mg/dL)",
                   BA = "bilirubin/albumin molar ratio",
                   Ugt = "UGT1A1 activity (fraction of adult)",
                   OD450 = "amniotic ΔOD450"))) +
      labs(x = "days from delivery", y = NULL) + THEME
  })

  ## --- 10 ------------------------------------------------------------------
  output$p_neo <- renderPlot({
    o <- sim()$out
    o$pna <- (o$ga - input$ga_del) * 7
    nb <- o[o$pna >= 0, ]
    ggplot(long(nb %>% mutate(ga = pna), c("HB", "HBMOM", "CFA1", "Fer")),
           aes(ga, value)) +
      geom_line(linewidth = 0.9, colour = "#1f7a1f") +
      facet_wrap(~variable, scales = "free_y",
                 labeller = as_labeller(c(
                   HB = "haemoglobin (g/dL)", HBMOM = "Hb MoM",
                   CFA1 = "residual maternal anti-D (IU/mL)",
                   Fer = "transfusional iron (mg)"))) +
      labs(x = "postnatal age (days)", y = NULL) + THEME
  })

  ## --- 11 ------------------------------------------------------------------
  prophy <- reactive({
    p <- as.list(param(hdfn))
    rbc <- 0.42 * input$fmh                 # mL of fetal RED CELLS
    need <- p$iu_per_ml_rbc * rbc           # IU required to coat them
    given <- input$rhig * p$ad_iu_per_ug    # IU given
    covered <- min(1, given / max(need, 1e-9))
    leak <- 1 - exp(-p$k_sen_free * input$delay)
    list(rbc = rbc, need = need, given = given, covered = covered, leak = leak)
  })
  output$tbl_prophy <- renderTable({
    q <- prophy()
    data.frame(
      quantity = c("fetal red cells in the bleed (mL)",
                   "anti-D required at 100 IU/mL RBC (IU)",
                   "anti-D given (IU)",
                   "fraction of the antigen that can be coated",
                   "exposure leaked before the dose arrives",
                   "dose needed for full coverage (µg)"),
      value = c(sprintf("%.2f", q$rbc), sprintf("%.0f", q$need),
                sprintf("%.0f", q$given), sprintf("%.1f%%", 100 * q$covered),
                sprintf("%.2f%%", 100 * q$leak),
                sprintf("%.0f", q$need / as.list(param(hdfn))$ad_iu_per_ug)))
  })
  output$p_prophy <- renderPlot({
    p <- as.list(param(hdfn))
    d <- expand.grid(fmh = 10^seq(-1, 2.3, 0.05),
                     dose = c(100, 300, 600, 1500))
    d$covered <- pmin(1, d$dose * p$ad_iu_per_ug /
                        (p$iu_per_ml_rbc * 0.42 * d$fmh))
    ggplot(d, aes(fmh, covered, colour = factor(dose))) +
      geom_line(linewidth = 0.9) + scale_x_log10() +
      geom_vline(xintercept = 30, linetype = 2) +
      labs(x = "fetomaternal haemorrhage (mL fetal whole blood, log scale)",
           y = "fraction of antigen coated", colour = "dose (µg)",
           title = "300 µg turns over at exactly 30 mL of fetal whole blood",
           subtitle = "which is the stated coverage of that dose — the rule is stoichiometric") +
      THEME
  })

  ## --- 12 ------------------------------------------------------------------
  sc <- eventReactive(input$run_sc, hdfn_scenarios(hdfn))
  output$tbl_sc <- renderTable(sc(), digits = 2)
}

shinyApp(ui, server)
