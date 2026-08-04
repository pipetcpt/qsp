## =============================================================================
##  vl_shiny_app.R
##  Visceral Leishmaniasis (kala-azar) QSP model — interactive dashboard
##  내장 리슈만편모충증 QSP 모델 — 인터랙티브 대시보드
##
##  Ten tabs, each built around one thing the model claims:
##    1  환자 프로파일    Patient profile — set up the patient, see the disease
##    2  약동학           PK — plasma and INTRAMACROPHAGE concentrations
##    3  두 개의 적분     The two integrals — spleen kill vs free-plasma toxicity
##    4  기생충 부하      Parasite burden, by organ and by population
##    5  숙주 면역        Host immunity — the priming vs IL-10 race
##    6  임상 종말점      Clinical endpoints — spleen, Hb, platelets, albumin
##    7  독성             Toxicity ledger by organ system
##    8  시나리오 비교    Scenario comparison across all 20 regimens
##    9  분리선           The separatrix — cure threshold vs CD4
##   10  바이오마커       Biomarkers and the serology trap
##
##  Run with:  shiny::runApp("vl_shiny_app.R")
##  Requires vl_mrgsolve_model.R in the same directory.
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

source("vl_mrgsolve_model.R")

DRUG_COLS <- c(`L-AmB` = "#C2452D", `d-AmB` = "#E08A76",
               miltefosine = "#C98A2B", paromomycin = "#2E8B8B",
               antimony = "#8A8A2E")
ORGAN_COLS <- c(spleen = "#B5342A", liver = "#C98A2B",
                marrow = "#2E6FA8", skin = "#4C9A5E")

## -----------------------------------------------------------------------------
##  UI
## -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel(paste("내장 리슈만편모충증 QSP 모델 —",
                   "Visceral Leishmaniasis (kala-azar), 73 ODEs")),
  tags$p(style = "color:#555;",
         paste("The drug and the parasite never meet in plasma.",
               "Every kill term in this model is driven by an intramacrophage",
               "concentration; amphotericin nephrotoxicity is driven by the",
               "free plasma concentration of the same dose.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 / Patient"),
      sliderInput("WT", "Body weight (kg)", 5, 90, 50, step = 1),
      sliderInput("CD40", "Baseline CD4 (cells/uL)", 20, 1200, 700, step = 10),
      checkboxInput("HIV", "HIV co-infection", FALSE),
      checkboxInput("ART", "On suppressive ART", FALSE),
      sliderInput("MALNUT", "Malnutrition index (0-1)", 0, 0.8, 0, step = 0.05),
      sliderInput("P0SCALE", "Presenting burden multiplier", 0.1, 8, 1,
                  step = 0.1),

      h4("요법 / Regimen"),
      selectInput("scenario", "Preset scenario",
                  choices = c(
                    "S01 untreated natural history" = "S01_untreated",
                    "S02 SSG 20 mg/kg x 30 d (East Africa)" = "S02_ssg_africa",
                    "S03 SSG 20 mg/kg x 30 d (Bihar, resistant)" = "S03_ssg_bihar_res",
                    "S04 AmB deoxycholate 15 mg/kg total" = "S04_damb_alt30",
                    "S05 L-AmB 10 mg/kg single dose" = "S05_lamb_single10",
                    "S06 L-AmB 21 mg/kg multidose" = "S06_lamb_21_multi",
                    "S07 L-AmB 3 mg/kg x 5" = "S07_lamb_5x3",
                    "S08 miltefosine 28 d, adult" = "S08_mil28_adult",
                    "S09 miltefosine 28 d, child, linear mg/kg" = "S09_mil28_child_linear",
                    "S10 miltefosine 28 d, child, allometric" = "S10_mil28_child_allom",
                    "S11 paromomycin 21 d (India)" = "S11_pm21",
                    "S12 paromomycin 21 d (East Africa)" = "S12_pm21_africa",
                    "S13 L-AmB 5 + miltefosine 7 d" = "S13_lamb5_mil7",
                    "S14 L-AmB 5 + paromomycin 10 d" = "S14_lamb5_pm10",
                    "S15 miltefosine 10 d + paromomycin 10 d" = "S15_mil10_pm10",
                    "S16 SSG + paromomycin 17 d" = "S16_ssg_pm17",
                    "S17 HIV-VL, L-AmB 10 mg/kg, no ART yet" = "S17_hiv_lamb10",
                    "S18 HIV-VL, L-AmB 30 + miltefosine 28 d" = "S18_hiv_lamb30_mil28",
                    "S19 PKDL, miltefosine 12 weeks" = "S19_pkdl_mil84",
                    "S20 PKDL, L-AmB 5 mg/kg weekly x 4" = "S20_pkdl_lamb"),
                  selected = "S05_lamb_single10"),
      sliderInput("days", "Follow-up (days)", 60, 720, 360, step = 30),

      h4("약리 파라미터 / Pharmacology"),
      sliderInput("FREL", "Liposome leakage fraction FREL", 0.02, 1.0, 0.30,
                  step = 0.02),
      helpText(paste("FREL is the whole of claim 1: the fraction of liposome",
                     "clearance that leaks drug into plasma instead of being",
                     "phagocytosed. Set it to 1.0 and liposomal amphotericin",
                     "becomes deoxycholate.")),
      sliderInput("RES_SB", "Antimony efflux resistance RES_SB", 1, 15, 1,
                  step = 0.5),
      sliderInput("ADHER_PCT", "Oral adherence (% of doses taken)", 40, 100,
                  100, step = 5),
      sliderInput("KTP", "T-cell priming rate KTP (1/h)", 0.001, 0.015, 0.005,
                  step = 0.0005),
      actionButton("go", "Simulate", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 환자 프로파일 Profile",
                 br(), verbatimTextOutput("profile"),
                 plotOutput("plot_overview", height = "420px")),
        tabPanel("2 약동학 PK",
                 br(),
                 helpText(paste("Top row: plasma. Bottom row: the",
                                "intramacrophage compartments where the",
                                "kill terms actually live. Note the scale",
                                "difference.")),
                 plotOutput("plot_pk", height = "560px")),
        tabPanel("3 두 개의 적분 Two integrals",
                 br(),
                 helpText(paste("The same mg/kg of amphotericin splits into",
                                "two integrals that move in OPPOSITE",
                                "directions when you encapsulate it.")),
                 plotOutput("plot_integrals", height = "420px"),
                 tableOutput("tab_integrals")),
        tabPanel("4 기생충 부하 Burden",
                 br(), plotOutput("plot_burden", height = "420px"),
                 plotOutput("plot_burden_pop", height = "300px")),
        tabPanel("5 숙주 면역 Immunity",
                 br(),
                 helpText(paste("The race: TMEM builds while antigen is",
                                "present; IL-10 only switches on above a",
                                "burden of ~60 units. Whichever crosses",
                                "first decides the outcome.")),
                 plotOutput("plot_immune", height = "480px")),
        tabPanel("6 임상 종말점 Endpoints",
                 br(), plotOutput("plot_clin", height = "520px"),
                 tableOutput("tab_outcome")),
        tabPanel("7 독성 Toxicity",
                 br(), plotOutput("plot_tox", height = "520px"),
                 tableOutput("tab_tox")),
        tabPanel("8 시나리오 비교 Scenarios",
                 br(),
                 actionButton("go_all", "Run all 20 scenarios",
                              class = "btn-warning"),
                 br(), br(), tableOutput("tab_all"),
                 plotOutput("plot_all", height = "420px")),
        tabPanel("9 분리선 Separatrix",
                 br(),
                 helpText(paste("Cure is not zero parasites. It is crossing a",
                                "line whose position depends on CD4 count.")),
                 actionButton("go_sep", "Compute separatrix",
                              class = "btn-warning"),
                 br(), br(), tableOutput("tab_sep"),
                 plotOutput("plot_sep", height = "360px")),
        tabPanel("10 바이오마커 Biomarkers",
                 br(),
                 helpText(paste("Polyclonal IgG decays with a half-life of",
                                "~96 days while the burden falls by orders of",
                                "magnitude in days. That gap is why rK39",
                                "serology cannot be used to call relapse.")),
                 plotOutput("plot_bio", height = "460px"),
                 tableOutput("tab_bio"))
      )
    )
  )
)

## -----------------------------------------------------------------------------
##  Server
## -----------------------------------------------------------------------------
server <- function(input, output, session) {

  user_par <- reactive({
    p <- vl_scenario_par(input$scenario)
    ## explicit UI settings win over the scenario preset, except for the
    ## paediatric scenarios where the weight IS the scenario
    if (!input$scenario %in% c("S09_mil28_child_linear",
                              "S10_mil28_child_allom")) {
      p$WT <- input$WT
      p$MALNUT <- input$MALNUT
    }
    if (!input$scenario %in% c("S17_hiv_lamb10", "S18_hiv_lamb30_mil28")) {
      p$CD40 <- input$CD40
      p$HIV <- as.numeric(input$HIV)
      p$ART <- as.numeric(input$ART)
    }
    p$P0SCALE <- input$P0SCALE
    p$FREL <- input$FREL
    p$RES_SB <- if (is.null(p$RES_SB)) input$RES_SB else
      max(p$RES_SB, input$RES_SB)
    p$KTP <- input$KTP
    p
  })

  events <- reactive({
    ev <- vl_regimens(input$WT)[[input$scenario]]
    if (!is.null(ev) && input$ADHER_PCT < 100 && any(ev$cmt == CMT_MIL)) {
      set.seed(1)
      keep <- ev$cmt != CMT_MIL | runif(nrow(ev)) <= input$ADHER_PCT / 100
      ev <- ev[keep, , drop = FALSE]
    }
    ev
  })

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    m <- vl_mod
    if (grepl("^S(19|20)", input$scenario)) m <- init(vl_mod, vl_pkdl_init())
    withProgress(message = "Integrating 73 ODEs...", {
      vl_simulate(mod = m, dosing = events(), par = user_par(),
                  days = input$days, delta = 6)
    })
  })

  outc <- reactive(vl_outcome(sim(), vl_eot_day(events())))

  ## --- tab 1 ---------------------------------------------------------------
  output$profile <- renderPrint({
    o <- outc(); s <- sim()
    cat(sprintf("Scenario            : %s\n", input$scenario))
    cat(sprintf("Weight / CD4        : %.0f kg / %.0f cells/uL%s\n",
                input$WT, input$CD40, if (input$HIV) "  [HIV+]" else ""))
    cat(sprintf("Cumulative dose     : %.1f mg/kg\n",
                if (is.null(events())) 0 else sum(events()$amt) / input$WT))
    cat(sprintf("Presenting burden   : %.3g x1e6 amastigotes\n", o$P0))
    cat(sprintf("Nadir               : %.3g  on day %.0f\n",
                o$nadir, o$nadir_day))
    cat(sprintf("Burden at day %-5.0f : %.3g\n", input$days, o$P_end))
    cat(sprintf("Log10 drop to nadir : %.2f\n", o$logdrop))
    cat(sprintf("OUTCOME             : %s\n",
                if (o$cure) "definitive cure" else
                  if (o$relapse) "RELAPSE" else "PRIMARY FAILURE"))
    cat(sprintf("Modelled mortality  : %.1f %%\n", o$mort_pct))
    cat(sprintf("Peak PKDL lesion    : %.2f\n", o$pkdl_max))
  })

  output$plot_overview <- renderPlot({
    s <- sim()
    df <- data.frame(day = s$day,
                     `log10 burden` = log10(pmax(s$PTOT_OUT, 1e-7)),
                     `spleen (cm)` = s$SPL, `Hb (g/dL)` = s$HGB,
                     `temperature (C)` = s$TEMP, check.names = FALSE) |>
      pivot_longer(-day)
    ggplot(df, aes(day, value)) + geom_line(linewidth = 0.7, colour = "#B5342A") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Day", y = NULL) + theme_bw(base_size = 12)
  })

  ## --- tab 2 ---------------------------------------------------------------
  output$plot_pk <- renderPlot({
    s <- sim()
    plasma <- data.frame(day = s$day,
                         `AmB liposomal` = s$CLIP_OUT,
                         `AmB free` = s$CFRE_OUT,
                         `miltefosine` = s$CMIL_OUT,
                         check.names = FALSE) |> pivot_longer(-day)
    tissue <- data.frame(day = s$day,
                         `AmB spleen macrophage` = s$CASP_OUT,
                         `miltefosine spleen macrophage` = s$CMSP_OUT,
                         `paromomycin spleen macrophage` = s$CPSP_OUT,
                         `Sb(III) spleen macrophage` = s$CSSP_OUT,
                         check.names = FALSE) |> pivot_longer(-day)
    ## faceted rather than gridExtra::grid.arrange, to keep the app's
    ## dependencies to shiny + mrgsolve + the tidyverse trio
    plasma$panel <- "Plasma"
    tissue$panel <- "Intramacrophage - where the kill terms live"
    ggplot(rbind(plasma, tissue),
           aes(day, value + 1e-4, colour = name)) +
      geom_line(linewidth = 0.7) + scale_y_log10() +
      facet_wrap(~panel, ncol = 1, scales = "free_y") +
      labs(x = "Day", y = "mg/L", colour = NULL) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom")
  })

  ## --- tab 3 ---------------------------------------------------------------
  integrals <- reactive({
    wt <- input$WT
    a <- vl_simulate(dosing = vl_ev(CMT_LAMB, mgkg = 10, n = 1, wt = wt),
                     par = list(WT = wt, FREL = input$FREL), days = 180,
                     delta = 4)
    b <- vl_simulate(dosing = vl_ev(CMT_DAMB, mgkg = 1, n = 10, ii_h = 48,
                                    wt = wt),
                     par = list(WT = wt), days = 180, delta = 4)
    list(lamb = a, damb = b)
  })

  output$plot_integrals <- renderPlot({
    z <- integrals()
    df <- rbind(
      data.frame(day = z$lamb$day, spleen = z$lamb$AUCASP,
                 free = z$lamb$AUCFRE, form = "liposomal"),
      data.frame(day = z$damb$day, spleen = z$damb$AUCASP,
                 free = z$damb$AUCFRE, form = "deoxycholate")) |>
      pivot_longer(c(spleen, free), names_to = "integral")
    ggplot(df, aes(day, value + 1, colour = form, linetype = integral)) +
      geom_line(linewidth = 0.8) + scale_y_log10() +
      scale_colour_manual(values = c(liposomal = "#C2452D",
                                     deoxycholate = "#2E6FA8")) +
      labs(x = "Day", y = "Cumulative AUC (mg.h/L)", colour = NULL,
           linetype = NULL,
           title = "One dose, two integrals",
           subtitle = paste("spleen-macrophage AUC drives killing;",
                            "free-plasma AUC drives nephrotoxicity")) +
      theme_bw(base_size = 12)
  })

  output$tab_integrals <- renderTable({
    z <- integrals()
    f <- function(s, lbl) data.frame(
      form = lbl,
      `AUC spleen` = s$AUCASP[nrow(s)],
      `AUC free plasma` = s$AUCFRE[nrow(s)],
      `TI surrogate` = s$AUCASP[nrow(s)] / s$AUCFRE[nrow(s)],
      `peak creatinine` = max(s$SCR), `nadir K` = min(s$KSER),
      check.names = FALSE)
    out <- rbind(f(z$lamb, "liposomal 10 mg/kg x1"),
                 f(z$damb, "deoxycholate 10 mg/kg total"))
    out$`gain` <- c(round(out$`TI surrogate`[1] / out$`TI surrogate`[2], 1), NA)
    out
  }, digits = 3)

  ## --- tab 4 ---------------------------------------------------------------
  output$plot_burden <- renderPlot({
    s <- sim()
    df <- data.frame(day = s$day,
                     spleen = s$P_SP_S + s$P_SP_R,
                     liver  = s$P_LI_S + s$P_LI_R,
                     marrow = s$P_BM_S + s$P_BM_R,
                     skin   = s$P_SK_S + s$P_SK_R) |> pivot_longer(-day)
    ggplot(df, aes(day, pmax(value, 1e-7), colour = name)) +
      geom_line(linewidth = 0.8) + scale_y_log10() +
      scale_colour_manual(values = ORGAN_COLS) +
      geom_hline(yintercept = EXTINCT, linetype = 3) +
      labs(x = "Day", y = "Amastigotes (x1e6)", colour = "Organ",
           title = "Burden by organ",
           subtitle = "dotted line = one amastigote (extinction floor)") +
      theme_bw(base_size = 12)
  })

  output$plot_burden_pop <- renderPlot({
    s <- sim()
    df <- data.frame(day = s$day,
                     replicating = s$P_SP_S + s$P_LI_S + s$P_BM_S + s$P_SK_S,
                     quiescent = s$P_SP_R + s$P_LI_R + s$P_BM_R + s$P_SK_R) |>
      pivot_longer(-day)
    ggplot(df, aes(day, pmax(value, 1e-7), colour = name)) +
      geom_line(linewidth = 0.8) + scale_y_log10() +
      labs(x = "Day", y = "Amastigotes (x1e6)", colour = NULL,
           title = "Replicating vs quiescent — where relapse comes from") +
      theme_bw(base_size = 12)
  })

  ## --- tab 5 ---------------------------------------------------------------
  output$plot_immune <- renderPlot({
    s <- sim()
    df <- data.frame(day = s$day, TMEM = s$TMEM, `IFN-gamma` = s$IFNG,
                     `IL-10` = s$IL10, `TNF-alpha` = s$TNFA,
                     `TGF-beta` = s$TGFB,
                     `activated macrophage` = s$MPHA,
                     `CD4 (/100)` = s$CD4 / 100, check.names = FALSE) |>
      pivot_longer(-day)
    ggplot(df, aes(day, value, colour = name)) +
      geom_line(linewidth = 0.8) +
      labs(x = "Day", y = "Level (arbitrary units)", colour = NULL,
           title = "Host immunity: the priming vs IL-10 race") +
      theme_bw(base_size = 12)
  })

  ## --- tab 6 ---------------------------------------------------------------
  output$plot_clin <- renderPlot({
    s <- sim()
    df <- data.frame(day = s$day,
                     `spleen (cm)` = s$SPL, `liver (cm)` = s$LIVS,
                     `Hb (g/dL)` = s$HGB, `platelets (1e9/L)` = s$PLT,
                     `WBC (1e9/L)` = s$WBC, `albumin (g/dL)` = s$ALB,
                     `temperature (C)` = s$TEMP,
                     `weight loss (frac)` = s$BWTL,
                     `splenic aspirate grade` = s$SPLGRADE,
                     check.names = FALSE) |> pivot_longer(-day)
    ggplot(df, aes(day, value)) + geom_line(linewidth = 0.7, colour = "#2E6FA8") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Day", y = NULL) + theme_bw(base_size = 11)
  })

  output$tab_outcome <- renderTable({
    o <- outc()
    data.frame(metric = c("nadir burden", "nadir day", "burden at end",
                          "log10 drop", "outcome", "mortality %",
                          "spleen (cm)", "Hb (g/dL)", "platelets"),
               value = c(sprintf("%.3g", o$nadir), sprintf("%.0f", o$nadir_day),
                         sprintf("%.3g", o$P_end), sprintf("%.2f", o$logdrop),
                         if (o$cure) "cure" else
                           if (o$relapse) "relapse" else "failure",
                         sprintf("%.1f", o$mort_pct),
                         sprintf("%.2f", o$spleen_end),
                         sprintf("%.1f", o$hgb_end),
                         sprintf("%.0f", o$plt_end)))
  })

  ## --- tab 7 ---------------------------------------------------------------
  output$plot_tox <- renderPlot({
    s <- sim()
    df <- data.frame(day = s$day,
                     `creatinine (mg/dL)` = s$SCR,
                     `potassium (mEq/L)` = s$KSER,
                     `magnesium (mg/dL)` = s$MGSER,
                     `hearing shift (dB)` = s$HEAR,
                     `QTc (ms)` = s$QTC,
                     `lipase (U/L)` = s$LIPA,
                     `ALT (U/L)` = s$ALTX,
                     `GI intolerance index` = s$GITX,
                     `tubular injury index` = s$TUBI,
                     check.names = FALSE) |> pivot_longer(-day)
    ggplot(df, aes(day, value)) + geom_line(linewidth = 0.7, colour = "#B5342A") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Day", y = NULL,
           title = "Toxicity: each panel belongs to a different drug") +
      theme_bw(base_size = 11)
  })

  output$tab_tox <- renderTable({
    o <- outc()
    data.frame(
      endpoint = c("peak creatinine (mg/dL)", "nadir potassium (mEq/L)",
                   "nadir magnesium (mg/dL)", "hearing shift (dB)",
                   "peak QTc (ms)", "peak lipase (U/L)", "peak ALT (U/L)",
                   "GI intolerance index"),
      attributable_to = c("amphotericin B", "amphotericin B", "amphotericin B",
                          "paromomycin", "antimony", "antimony",
                          "antimony / miltefosine", "miltefosine"),
      value = c(o$scr_max, o$k_min, o$mg_min, o$hear_dB, o$qtc_max,
                o$lipase_max, o$alt_max, o$gi_max))
  }, digits = 2)

  ## --- tab 8 ---------------------------------------------------------------
  all_scen <- eventReactive(input$go_all, {
    withProgress(message = "Running 20 scenarios...", vl_run_all(input$WT,
                                                                input$days))
  })
  output$tab_all <- renderTable(all_scen(), digits = 3)
  output$plot_all <- renderPlot({
    d <- all_scen()
    d$scenario <- factor(d$scenario, levels = d$scenario[order(d$log_drop)])
    ggplot(d, aes(scenario, -log_drop, fill = outcome)) +
      geom_col() + coord_flip() +
      scale_fill_manual(values = c(cure = "#4C9A5E", relapse = "#C98A2B",
                                   failure = "#B5342A")) +
      labs(x = NULL, y = "log10 burden reduction to nadir") +
      theme_bw(base_size = 11)
  })

  ## --- tab 9 ---------------------------------------------------------------
  sep <- eventReactive(input$go_sep, {
    withProgress(message = "Bisecting the separatrix...",
                 vl_analysis_separatrix())
  })
  output$tab_sep <- renderTable(sep(), digits = 6)
  output$plot_sep <- renderPlot({
    d <- sep()
    ggplot(d, aes(CD4, critical_amastigotes)) +
      geom_line(linewidth = 0.9, colour = "#B5342A") +
      geom_point(size = 2.5) + scale_y_log10() +
      labs(x = "CD4 count (cells/uL)",
           y = "Critical residual burden (amastigotes)",
           title = "The finish line moves with CD4",
           subtitle = paste("residual burdens below the line are cleared by",
                            "the host; above it, the patient relapses")) +
      theme_bw(base_size = 12)
  })

  ## --- tab 10 --------------------------------------------------------------
  output$plot_bio <- renderPlot({
    s <- sim()
    b <- s$PTOT_OUT[1]
    df <- data.frame(
      day = s$day,
      `parasite burden (% of baseline)` = 100 * pmax(s$PTOT_OUT, 1e-9) / b,
      `polyclonal IgG (% of excess)` =
        100 * pmax(s$IGG - 1.1, 0) / max(s$IGG[1] - 1.1, 1e-9),
      `spleen size (% of baseline)` = 100 * s$SPL / s$SPL[1],
      check.names = FALSE) |> pivot_longer(-day)
    ggplot(df, aes(day, pmax(value, 1e-7), colour = name)) +
      geom_line(linewidth = 0.9) + scale_y_log10() +
      labs(x = "Day", y = "% of baseline (log scale)", colour = NULL,
           title = "The serology trap",
           subtitle = paste("burden falls by orders of magnitude in days;",
                            "IgG has a half-life of ~96 days")) +
      theme_bw(base_size = 12)
  })

  output$tab_bio <- renderTable({
    s <- sim()
    idx <- sapply(c(0, 30, 90, 180, 360), function(d) which.min(abs(s$day - d)))
    data.frame(day = s$day[idx],
               `log10 burden` = log10(pmax(s$PTOT_OUT[idx], 1e-9)),
               `splenic aspirate grade` = s$SPLGRADE[idx],
               `IgG (g/dL)` = s$IGG[idx],
               `spleen (cm)` = s$SPL[idx],
               `Hb (g/dL)` = s$HGB[idx],
               check.names = FALSE)
  }, digits = 2)
}

shinyApp(ui, server)
