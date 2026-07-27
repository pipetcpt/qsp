## ============================================================================
## Endometrial Carcinoma (EC) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 8 tabs: Patient & molecular profile · Drug PK · Endocrine & PI3K PD ·
##         Tumour dynamics · Immune compartment · Clinical endpoints ·
##         Scenario comparison · TMB continuum (the model's central claim) ·
##         References
## Dependencies: shiny, shinydashboard, mrgsolve, dplyr, tidyr, ggplot2, DT
## Run with:  shiny::runApp("emc_shiny_app.R")
##            (expects emc_mrgsolve_model.R in the same directory)
## ----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(DT)
})

## ---------- Lazy model load ----------
get_model <- function() {
  if (!exists(".EMC_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.EMC_MOD)) {
    assign(".EMC_MOD",
           mread_cache("emc", project = ".", file = "emc_mrgsolve_model.R"),
           envir = .GlobalEnv)
  }
  .GlobalEnv$.EMC_MOD
}

CLASS_LABELS <- c("1 · POLEmut (ultramutated)"      = 1,
                  "2 · MMRd / MSI-H (hypermutated)" = 2,
                  "3 · NSMP / copy-number-low"      = 3,
                  "4 · p53abn / copy-number-high"   = 4)

CLASS_SHORT <- c("POLEmut", "MMRd", "NSMP", "p53abn")

SCENARIOS <- c(
  "1. Untreated natural history",
  "2. Carboplatin-paclitaxel x6 (GOG-0209)",
  "3. Dostarlimab + chemo (RUBY)",
  "4. Pembrolizumab + chemo (NRG-GY018)",
  "5. Pembrolizumab + lenvatinib (KEYNOTE-775)",
  "6. Fertility-sparing LNG-IUS (feMMe)",
  "7. Letrozole + everolimus (Slomovitz 2015)",
  "8. Megestrol +/- tamoxifen alternation (GOG-153)"
)

## ---------- Dosing building blocks ----------
ev_chemo <- function(cycles = 6) {
  ev(amt = 625, cmt = "CARB_C", ii = 21, addl = cycles - 1) +
    ev(amt = 300, cmt = "PTX_C", ii = 21, addl = cycles - 1)
}
ev_dostarlimab <- function(end_d) {
  ev(amt = 500, cmt = "IO_C", ii = 21, addl = 3) +
    ev(amt = 1000, cmt = "IO_C", time = 84, ii = 42,
       addl = max(floor((end_d - 84) / 42), 0))
}
ev_pembro <- function(end_d) {
  ev(amt = 200, cmt = "IO_C", ii = 21, addl = max(floor(end_d / 21) - 1, 0))
}
ev_lenva <- function(end_d, dose = 20) {
  ev(amt = dose, cmt = "LEN_G", ii = 1, addl = max(end_d - 1, 0))
}
ev_megace <- function(end_d, dose = 160) {
  ev(amt = dose, cmt = "PRG_G", ii = 1, addl = max(end_d - 1, 0))
}
ev_ius <- function(end_d) {
  ev(amt = 0.020, cmt = "PRG_G", ii = 1, addl = max(end_d - 1, 0))
}
ev_letro <- function(end_d) ev(amt = 2.5, cmt = "LTZ_G", ii = 1, addl = max(end_d - 1, 0))
ev_evero <- function(end_d) ev(amt = 10, cmt = "EVE_G", ii = 1, addl = max(end_d - 1, 0))

## Returns list(ev = <event object>, p = <named list of parameter overrides>)
build_regimen <- function(scenario, end_d, cycles = 6, len_dose = 20) {
  switch(
    scenario,
    "1. Untreated natural history" =
      list(ev = ev(amt = 0, cmt = "CARB_C", time = 0), p = list()),
    "2. Carboplatin-paclitaxel x6 (GOG-0209)" =
      list(ev = ev_chemo(cycles), p = list()),
    "3. Dostarlimab + chemo (RUBY)" =
      list(ev = ev_chemo(cycles) + ev_dostarlimab(end_d), p = list()),
    "4. Pembrolizumab + chemo (NRG-GY018)" =
      list(ev = ev_chemo(cycles) + ev_pembro(end_d), p = list()),
    "5. Pembrolizumab + lenvatinib (KEYNOTE-775)" =
      list(ev = ev_pembro(end_d) + ev_lenva(end_d, len_dose), p = list()),
    "6. Fertility-sparing LNG-IUS (feMMe)" =
      list(ev = ev_ius(min(end_d, 360)),
           p = list(GRADE = 1, TUM0 = 12, POSTMENO = 0,
                    LOCAL_IUS = 1, ERPR_POS = 1)),
    "7. Letrozole + everolimus (Slomovitz 2015)" =
      list(ev = ev_letro(end_d) + ev_evero(end_d), p = list(ERPR_POS = 1)),
    "8. Megestrol +/- tamoxifen alternation (GOG-153)" =
      list(ev = ev_megace(min(end_d, 360)), p = list(ERPR_POS = 1)),
    list(ev = ev(amt = 0, cmt = "CARB_C"), p = list())
  )
}

run_sim <- function(mod, scenario, pars, end_d, cycles, len_dose) {
  reg <- build_regimen(scenario, end_d, cycles, len_dose)
  allp <- modifyList(pars, reg$p)
  m <- do.call(param, c(list(mod), allp))
  mrgsim(m, reg$ev, end = end_d, delta = 1) %>%
    as_tibble() %>%
    mutate(scenario = scenario)
}

theme_emc <- function() {
  theme_minimal(base_size = 13) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}

## ============================================================================
## UI
## ============================================================================
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "Endometrial Cancer QSP", titleWidth = 320),

  dashboardSidebar(
    width = 320,
    sidebarMenu(
      id = "tabs",
      menuItem("Patient & molecular profile", tabName = "profile", icon = icon("dna")),
      menuItem("Drug PK",                     tabName = "pk",      icon = icon("prescription-bottle")),
      menuItem("Endocrine & PI3K PD",         tabName = "pd",      icon = icon("venus")),
      menuItem("Tumour dynamics",             tabName = "tumour",  icon = icon("chart-line")),
      menuItem("Immune compartment",          tabName = "immune",  icon = icon("shield-virus")),
      menuItem("Clinical endpoints",          tabName = "endpoint", icon = icon("notes-medical")),
      menuItem("Scenario comparison",         tabName = "compare", icon = icon("layer-group")),
      menuItem("TMB continuum",               tabName = "tmb",     icon = icon("wave-square")),
      menuItem("References",                  tabName = "refs",    icon = icon("book"))
    ),
    hr(),
    selectInput("molclass", "TCGA molecular class",
                choices = CLASS_LABELS, selected = 2),
    selectInput("scenario", "Treatment scenario",
                choices = SCENARIOS, selected = SCENARIOS[3]),
    sliderInput("horizon", "Simulation horizon (days)",
                min = 90, max = 1095, value = 730, step = 30),
    sliderInput("bmi", "BMI (kg/m2)", min = 20, max = 55, value = 34, step = 1),
    sliderInput("tum0", "Baseline RECIST sum (mm)",
                min = 10, max = 200, value = 60, step = 5),
    sliderInput("grade", "FIGO grade", min = 1, max = 3, value = 2, step = 1),
    checkboxInput("erpr", "ER+/PR+ by IHC", value = TRUE),
    checkboxInput("postmeno", "Postmenopausal", value = TRUE),
    checkboxInput("her2amp", "HER2 amplified (serous)", value = FALSE),
    checkboxInput("her2tx", "Trastuzumab on board", value = FALSE),
    checkboxInput("b2m", "B2M loss (antigen-presentation escape)", value = FALSE),
    checkboxInput("tam", "Tamoxifen alternation (re-induces PR-B)", value = FALSE),
    sliderInput("cycles", "Chemotherapy cycles", min = 3, max = 8, value = 6, step = 1),
    sliderInput("lendose", "Lenvatinib dose (mg/day)",
                min = 4, max = 24, value = 20, step = 2),
    sliderInput("tmbovr", "TMB override (0 = use class default)",
                min = 0, max = 200, value = 0, step = 1)
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f7f5fa; }
      .small-box { border-radius: 6px; }
      .box { border-top-color: #7b4f9d; }
    "))),
    tabItems(

      ## ---------------- 1. Patient & molecular profile ----------------
      tabItem(
        "profile",
        fluidRow(
          valueBoxOutput("vb_class", width = 3),
          valueBoxOutput("vb_tmb", width = 3),
          valueBoxOutput("vb_immun", width = 3),
          valueBoxOutput("vb_excl", width = 3)
        ),
        fluidRow(
          box(width = 7, title = "The TMB → immunogenicity curve, and where this patient sits",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_hill", height = 380),
              helpText(HTML(
                "The four TCGA classes are <b>four points on one continuous axis</b>. ",
                "The dashed vertical line at TMB50 = 13 mut/Mb is where the clinical ",
                "dMMR/pMMR dichotomy effectively falls — on the steepest part of the ",
                "curve. That is why the binary biomarker works, and also why tumours ",
                "near the boundary are the ones it must misclassify."))),
          box(width = 5, title = "Molecular class reference card", status = "primary",
              solidHeader = TRUE, DTOutput("t_classes"))
        ),
        fluidRow(
          box(width = 12, title = "Baseline patient state", status = "info",
              solidHeader = TRUE, DTOutput("t_baseline"))
        )
      ),

      ## ---------------- 2. Drug PK ----------------
      tabItem(
        "pk",
        fluidRow(
          box(width = 6, title = "Cytotoxic chemotherapy plasma concentrations",
              status = "primary", solidHeader = TRUE, plotOutput("p_pk_chemo", height = 320)),
          box(width = 6, title = "Anti-PD-1 concentration and receptor blockade",
              status = "primary", solidHeader = TRUE, plotOutput("p_pk_io", height = 320))
        ),
        fluidRow(
          box(width = 6, title = "Lenvatinib concentration and VEGFR-2 inhibition",
              status = "primary", solidHeader = TRUE, plotOutput("p_pk_len", height = 320)),
          box(width = 6, title = "Endocrine agents (letrozole, everolimus)",
              status = "primary", solidHeader = TRUE, plotOutput("p_pk_endo", height = 320))
        )
      ),

      ## ---------------- 3. Endocrine & PI3K PD ----------------
      tabItem(
        "pd",
        fluidRow(
          box(width = 6, title = "Free estradiol", status = "primary",
              solidHeader = TRUE, plotOutput("p_e2", height = 300)),
          box(width = 6, title = "ERa transcriptional activity", status = "primary",
              solidHeader = TRUE, plotOutput("p_er", height = 300))
        ),
        fluidRow(
          box(width = 6, title = "PR-B expression — the acquired-resistance clock",
              status = "warning", solidHeader = TRUE, plotOutput("p_pr", height = 300),
              helpText("Sustained progestin occupancy erodes PR-B, which is the ",
                       "receptor the drug needs. Tamoxifen alternation re-induces it.")),
          box(width = 6, title = "PI3K/AKT/mTOR pathway activity", status = "primary",
              solidHeader = TRUE, plotOutput("p_pi3k", height = 300))
        )
      ),

      ## ---------------- 4. Tumour dynamics ----------------
      tabItem(
        "tumour",
        fluidRow(
          box(width = 8, title = "RECIST sum of target lesions", status = "primary",
              solidHeader = TRUE, plotOutput("p_tum", height = 360)),
          box(width = 4, title = "Response summary", status = "info",
              solidHeader = TRUE, DTOutput("t_response"))
        ),
        fluidRow(
          box(width = 6, title = "Sensitive vs resistant clone", status = "primary",
              solidHeader = TRUE, plotOutput("p_clones", height = 300)),
          box(width = 6, title = "Angiogenesis: VEGF-A and microvessel density",
              status = "primary", solidHeader = TRUE, plotOutput("p_angio", height = 300))
        )
      ),

      ## ---------------- 5. Immune compartment ----------------
      tabItem(
        "immune",
        fluidRow(
          box(width = 6, title = "Neoantigen load", status = "primary",
              solidHeader = TRUE, plotOutput("p_neo", height = 300)),
          box(width = 6, title = "CD8 effector vs exhausted T cells", status = "primary",
              solidHeader = TRUE, plotOutput("p_tcell", height = 300))
        ),
        fluidRow(
          box(width = 6, title = "Immunosuppressive compartment (Treg, MDSC)",
              status = "warning", solidHeader = TRUE, plotOutput("p_supp", height = 300),
              helpText("Lenvatinib's contribution in pMMR disease is modelled here ",
                       "and nowhere else — it has no direct tumour-kill term.")),
          box(width = 6, title = "PD-1 pathway blockade", status = "primary",
              solidHeader = TRUE, plotOutput("p_pd1", height = 300))
        )
      ),

      ## ---------------- 6. Clinical endpoints ----------------
      tabItem(
        "endpoint",
        fluidRow(
          valueBoxOutput("vb_nadir", width = 3),
          valueBoxOutput("vb_best", width = 3),
          valueBoxOutput("vb_surv", width = 3),
          valueBoxOutput("vb_anc", width = 3)
        ),
        fluidRow(
          box(width = 6, title = "Progression-free survival surrogate",
              status = "success", solidHeader = TRUE, plotOutput("p_surv", height = 300)),
          box(width = 6, title = "CA-125 and ctDNA", status = "primary",
              solidHeader = TRUE, plotOutput("p_bio", height = 300))
        ),
        fluidRow(
          box(width = 6, title = "Neutropenia (Friberg myelosuppression)",
              status = "danger", solidHeader = TRUE, plotOutput("p_anc", height = 300)),
          box(width = 6, title = "Systolic blood pressure (lenvatinib)",
              status = "danger", solidHeader = TRUE, plotOutput("p_sbp", height = 300))
        )
      ),

      ## ---------------- 7. Scenario comparison ----------------
      tabItem(
        "compare",
        fluidRow(
          box(width = 12, status = "primary", solidHeader = TRUE,
              title = "Select scenarios to overlay (all share the sidebar patient profile)",
              checkboxGroupInput("cmp_sc", NULL, choices = SCENARIOS,
                                 selected = SCENARIOS[c(1, 2, 3, 5)], inline = FALSE),
              checkboxInput("cmp_allclass",
                            "Instead compare the FOUR molecular classes under the single selected scenario",
                            value = FALSE))
        ),
        fluidRow(
          box(width = 6, title = "Tumour burden", status = "primary",
              solidHeader = TRUE, plotOutput("p_cmp_tum", height = 340)),
          box(width = 6, title = "Survival surrogate", status = "success",
              solidHeader = TRUE, plotOutput("p_cmp_surv", height = 340))
        ),
        fluidRow(
          box(width = 6, title = "CD8 effector density", status = "primary",
              solidHeader = TRUE, plotOutput("p_cmp_teff", height = 320)),
          box(width = 6, title = "Endpoint table", status = "info",
              solidHeader = TRUE, DTOutput("t_cmp"))
        )
      ),

      ## ---------------- 8. TMB continuum ----------------
      tabItem(
        "tmb",
        fluidRow(
          box(width = 12, status = "primary", solidHeader = TRUE,
              title = "The central claim: checkpoint efficacy is continuous in TMB",
              HTML("<p>Everything else is held fixed — same drug, same dose, same ",
                   "patient, same immune-exclusion penalty. Only the mutational ",
                   "burden moves. If the dMMR/pMMR dichotomy were a biological ",
                   "edge, this sweep would show a step. It shows a smooth ",
                   "sigmoid whose steepest region happens to be where the clinical ",
                   "cut-off was drawn.</p>"),
              sliderInput("tmb_lo", "Sweep from (mut/Mb)", min = 0.5, max = 10,
                          value = 1, step = 0.5),
              sliderInput("tmb_hi", "Sweep to (mut/Mb)", min = 20, max = 300,
                          value = 200, step = 10),
              actionButton("run_sweep", "Run TMB sweep", class = "btn-primary"))
        ),
        fluidRow(
          box(width = 6, title = "Best percentage change from baseline vs TMB",
              status = "primary", solidHeader = TRUE, plotOutput("p_sweep_resp", height = 340)),
          box(width = 6, title = "Survival surrogate at end of horizon vs TMB",
              status = "success", solidHeader = TRUE, plotOutput("p_sweep_surv", height = 340))
        ),
        fluidRow(
          box(width = 12, title = "Sweep results", status = "info",
              solidHeader = TRUE, DTOutput("t_sweep"))
        )
      ),

      ## ---------------- 9. References ----------------
      tabItem(
        "refs",
        fluidRow(
          box(width = 12, title = "Key calibration anchors", status = "primary",
              solidHeader = TRUE,
              HTML("
<ul>
<li><b>TCGA</b> — Kandoth 2013 <i>Nature</i> 497:67. The four molecular classes and their mutation rates.</li>
<li><b>ProMisE</b> — Talhouk 2015 <i>Br J Cancer</i> 113:299; 2017 <i>Cancer</i> 123:802. The clinically deployable surrogate classifier.</li>
<li><b>PORTEC-3 molecular</b> — Le&oacute;n-Castillo 2020 <i>J Clin Oncol</i> 38:3388. Adjuvant chemo benefit concentrated in p53abn; POLEmut excellent regardless.</li>
<li><b>RUBY</b> — Mirza 2023 <i>N Engl J Med</i> 388:2145. Dostarlimab + chemo; dMMR PFS HR 0.28.</li>
<li><b>NRG-GY018</b> — Eskander 2023 <i>N Engl J Med</i> 388:2159. Pembrolizumab + chemo; dMMR HR 0.30, pMMR HR 0.54.</li>
<li><b>GARNET</b> — Oaknin 2020 <i>JAMA Oncol</i> 6:1766. Dostarlimab monotherapy; ORR 45.5% dMMR vs 15.4% pMMR.</li>
<li><b>KEYNOTE-775</b> — Makker 2022 <i>N Engl J Med</i> 386:437. Lenvatinib + pembrolizumab in pMMR; PFS HR 0.60, OS HR 0.68.</li>
<li><b>GOG-0209</b> — Miller 2020 <i>J Clin Oncol</i> 38:3841. Carboplatin-paclitaxel reference doublet.</li>
<li><b>Everolimus + letrozole</b> — Slomovitz 2015 <i>J Clin Oncol</i> 33:930. ORR 32%.</li>
<li><b>feMMe / LNG-IUS</b> — Janda 2021 <i>Gynecol Oncol</i> 161:143; Westin 2021 <i>Am J Obstet Gynecol</i> 224:191.e1.</li>
<li><b>GOG-119 / GOG-153</b> — Fiorica 2004; Whitney 2004 <i>Gynecol Oncol</i> 92:10 / 92:4. Megestrol-tamoxifen alternation.</li>
<li><b>Friberg myelosuppression</b> — Friberg 2002 <i>J Clin Oncol</i> 20:4713.</li>
<li><b>RAINBO</b> — RAINBO Research Consortium 2023 <i>Int J Gynecol Cancer</i> 33:109. Molecular-class-directed de-escalation and escalation.</li>
</ul>
<p>The full annotated bibliography (75 entries with PubMed links) is in
<code>emc_references.md</code> in this directory.</p>")),
          box(width = 12, title = "Disclaimer", status = "warning", solidHeader = TRUE,
              HTML("<p>This is a qualitative-to-semi-quantitative QSP model built for
                    education, research and hypothesis generation. It has not been
                    independently validated or qualified and must not be used for
                    clinical decisions, prescribing, or regulatory submission.</p>"))
        )
      )
    )
  )
)

## ============================================================================
## SERVER
## ============================================================================
server <- function(input, output, session) {

  base_params <- reactive({
    list(
      MOLCLASS = as.numeric(input$molclass),
      BMI      = input$bmi,
      TUM0     = input$tum0,
      GRADE    = input$grade,
      ERPR_POS = as.numeric(input$erpr),
      POSTMENO = as.numeric(input$postmeno),
      HER2_AMP = as.numeric(input$her2amp),
      HER2_TX  = as.numeric(input$her2tx),
      B2M_LOSS = as.numeric(input$b2m),
      TAM_PULSE = as.numeric(input$tam),
      TMB_OVR  = input$tmbovr
    )
  })

  sim <- reactive({
    withProgress(message = "Simulating...", value = 0.5, {
      run_sim(get_model(), input$scenario, base_params(),
              input$horizon, input$cycles, input$lendose)
    })
  })

  ## ---------------- Value boxes ----------------
  output$vb_class <- renderValueBox({
    cl <- as.numeric(input$molclass)
    valueBox(CLASS_SHORT[cl], "TCGA molecular class",
             icon = icon("dna"),
             color = c("green", "aqua", "yellow", "red")[cl])
  })
  output$vb_tmb <- renderValueBox({
    d <- sim()
    valueBox(sprintf("%.1f", d$TMB_mutMb[1]), "TMB (mut/Mb)",
             icon = icon("wave-square"), color = "purple")
  })
  output$vb_immun <- renderValueBox({
    d <- sim()
    valueBox(sprintf("%.2f", d$Immunogenicity[1]),
             "Immunogenicity (Hill output, 0-1)",
             icon = icon("shield-virus"), color = "olive")
  })
  output$vb_excl <- renderValueBox({
    cl <- as.numeric(input$molclass)
    valueBox(ifelse(cl == 4, "0.45 (excluded)", "0.12 (baseline)"),
             "Immune-exclusion penalty",
             icon = icon("ban"), color = ifelse(cl == 4, "red", "light-blue"))
  })

  output$vb_nadir <- renderValueBox({
    d <- sim()
    valueBox(sprintf("%.0f mm", min(d$Tumor_total_mm)), "Tumour nadir",
             icon = icon("arrow-down"), color = "green")
  })
  output$vb_best <- renderValueBox({
    d <- sim()
    valueBox(sprintf("%+.0f%%", min(d$Pct_change_baseline)),
             "Best change from baseline", icon = icon("percent"), color = "aqua")
  })
  output$vb_surv <- renderValueBox({
    d <- sim()
    valueBox(sprintf("%.0f%%", 100 * tail(d$Survival_fraction, 1)),
             sprintf("Survival surrogate at day %d", input$horizon),
             icon = icon("heart-pulse"), color = "purple")
  })
  output$vb_anc <- renderValueBox({
    d <- sim()
    valueBox(sprintf("%.2f", min(d$ANC_10e9L)), "ANC nadir (10^9/L)",
             icon = icon("droplet"),
             color = ifelse(min(d$ANC_10e9L) < 1.0, "red", "yellow"))
  })

  ## ---------------- Tab 1: profile ----------------
  output$p_hill <- renderPlot({
    p <- get_model()@param@data
    tmb50 <- p$TMB50; h <- p$HILL_TMB
    grid <- tibble(TMB = 10^seq(log10(0.5), log10(300), length.out = 400)) %>%
      mutate(Immunogenicity = TMB^h / (tmb50^h + TMB^h))
    pts <- tibble(
      TMB   = c(p$TMB_P53, p$TMB_NSMP, p$TMB_MMRD, p$TMB_POLE),
      Class = factor(c("p53abn", "NSMP", "MMRd", "POLEmut"),
                     levels = c("p53abn", "NSMP", "MMRd", "POLEmut"))
    ) %>% mutate(Immunogenicity = TMB^h / (tmb50^h + TMB^h))
    here <- tibble(TMB = sim()$TMB_mutMb[1]) %>%
      mutate(Immunogenicity = TMB^h / (tmb50^h + TMB^h))

    ggplot(grid, aes(TMB, Immunogenicity)) +
      geom_vline(xintercept = tmb50, linetype = 2, colour = "grey40") +
      annotate("text", x = tmb50 * 1.15, y = 0.06, hjust = 0, size = 3.6,
               colour = "grey30",
               label = sprintf("TMB50 = %g mut/Mb\n(where the dMMR/pMMR\ncut effectively falls)", tmb50)) +
      geom_line(linewidth = 1.2, colour = "#7b4f9d") +
      geom_point(data = pts, aes(colour = Class), size = 5) +
      geom_point(data = here, shape = 21, size = 8, stroke = 1.4,
                 fill = NA, colour = "black") +
      scale_colour_manual(values = c("p53abn" = "#c0392b", "NSMP" = "#b7950b",
                                     "MMRd" = "#2471a3", "POLEmut" = "#1e8449")) +
      scale_x_log10(breaks = c(1, 2, 5, 13, 25, 50, 140, 300)) +
      labs(x = "Tumour mutational burden (mut/Mb, log scale)",
           y = "Immunogenicity (0-1)",
           title = "One curve, four classes, no per-class fitting") +
      theme_emc()
  })

  output$t_classes <- renderDT({
    datatable(
      data.frame(
        Class = CLASS_SHORT,
        Prevalence = c("7-12%", "25-30%", "30-40%", "15-25%"),
        `TMB (mut/Mb)` = c(140, 25, 2.9, 2.3),
        `5-y PFS` = c("~95-100%", "~75%", "~80%", "~50%"),
        `Therapeutic implication` = c(
          "De-escalation candidate (RAINBO POLE-BLUE)",
          "Checkpoint blockade; RUBY HR 0.28",
          "Hormonal therapy if ER+/PR+; lenvatinib+IO",
          "Chemotherapy benefit concentrated here; immune-excluded"),
        check.names = FALSE),
      rownames = FALSE, options = list(dom = "t", ordering = FALSE))
  })

  output$t_baseline <- renderDT({
    d <- sim() %>% slice(1)
    datatable(
      data.frame(
        Quantity = c("TMB (mut/Mb)", "Immunogenicity", "Free estradiol (pmol/L)",
                     "ERa activity", "PR-B expression", "PI3K/AKT/mTOR activity",
                     "Tumour burden (mm)", "CD8 effector density",
                     "CA-125 (U/mL)", "ctDNA (%VAF)"),
        Value = round(c(d$TMB_mutMb, d$Immunogenicity, d$Estradiol_free,
                        d$ERa_activity, d$PRB_expression, d$PI3K_activity,
                        d$Tumor_total_mm, d$CD8_effector,
                        d$CA125_UmL, d$ctDNA_VAF), 3)),
      rownames = FALSE, options = list(dom = "t", pageLength = 12))
  })

  ## ---------------- Tab 2: PK ----------------
  output$p_pk_chemo <- renderPlot({
    sim() %>%
      select(time, Carboplatin = Carboplatin_conc, Paclitaxel = Paclitaxel_conc) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Time (days)", y = "Concentration (mg/L)") + theme_emc()
  })
  output$p_pk_io <- renderPlot({
    sim() %>%
      select(time, `Concentration (mg/L)` = AntiPD1_conc,
             `PD-1 blockade (0-1)` = PD1_blockade) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Time (days)", y = NULL) + theme_emc()
  })
  output$p_pk_len <- renderPlot({
    sim() %>%
      select(time, `Lenvatinib (mg/L)` = Lenvatinib_conc,
             `VEGFR-2 inhibition (0-1)` = VEGFR2_inhibition) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Time (days)", y = NULL) + theme_emc()
  })
  output$p_pk_endo <- renderPlot({
    sim() %>%
      select(time, Letrozole = Letrozole_conc, Everolimus = Everolimus_conc) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Time (days)", y = "Concentration (mg/L)") + theme_emc()
  })

  ## ---------------- Tab 3: endocrine / PI3K ----------------
  simple_line <- function(col, ylab, colour = "#7b4f9d") {
    renderPlot({
      d <- sim()
      ggplot(d, aes(time, .data[[col]])) +
        geom_line(linewidth = 1.1, colour = colour) +
        labs(x = "Time (days)", y = ylab) + theme_emc()
    })
  }
  output$p_e2   <- simple_line("Estradiol_free", "Free estradiol (pmol/L)", "#8e44ad")
  output$p_er   <- simple_line("ERa_activity", "ERa transcriptional activity", "#8e44ad")
  output$p_pr   <- simple_line("PRB_expression", "PR-B expression (fraction)", "#b7950b")
  output$p_pi3k <- simple_line("PI3K_activity", "PI3K/AKT/mTOR activity", "#1e8449")
  output$p_pd1  <- simple_line("PD1_blockade", "PD-1 pathway blockade (0-1)", "#2471a3")
  output$p_neo  <- simple_line("Neoantigen_load", "Neoantigen load (normalised)", "#16a085")
  output$p_sbp  <- simple_line("SBP_mmHg", "Systolic BP (mmHg)", "#c0392b")

  ## ---------------- Tab 4: tumour ----------------
  output$p_tum <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, Tumor_total_mm)) +
      geom_hline(yintercept = d$Tumor_total_mm[1] * 1.2, linetype = 2,
                 colour = "#c0392b") +
      geom_hline(yintercept = d$Tumor_total_mm[1] * 0.7, linetype = 2,
                 colour = "#1e8449") +
      geom_line(linewidth = 1.2, colour = "#7b4f9d") +
      annotate("text", x = max(d$time) * 0.02, y = d$Tumor_total_mm[1] * 1.22,
               hjust = 0, size = 3.4, colour = "#c0392b",
               label = "RECIST progression threshold (+20%)") +
      annotate("text", x = max(d$time) * 0.02, y = d$Tumor_total_mm[1] * 0.66,
               hjust = 0, size = 3.4, colour = "#1e8449",
               label = "RECIST partial response (-30%)") +
      labs(x = "Time (days)", y = "Sum of target lesions (mm)") + theme_emc()
  })

  output$t_response <- renderDT({
    d <- sim()
    b <- d$Tumor_total_mm[1]
    nadir <- min(d$Tumor_total_mm)
    best <- min(d$Pct_change_baseline)
    ## RECIST progression is measured FROM the nadir, so only timepoints after
    ## the nadir can qualify. Scanning the whole series would flag the
    ## pre-treatment baseline of any responding patient as progression.
    nad_i <- which.min(d$Tumor_total_mm)
    post <- d[seq(nad_i, nrow(d)), ]
    pd_idx <- which(post$Tumor_total_mm >= nadir * 1.2 & post$Tumor_total_mm >= nadir + 5)
    ttp <- if (length(pd_idx)) post$time[pd_idx[1]] else NA_real_
    resp <- if (best <= -100 + 1e-6) "CR" else if (best <= -30) "PR" else
      if (best >= 20) "PD" else "SD"
    datatable(
      data.frame(
        Metric = c("Baseline (mm)", "Nadir (mm)", "Best change (%)",
                   "Best RECIST response", "Time to progression (d)",
                   "Survival surrogate at horizon"),
        Value = c(sprintf("%.0f", b), sprintf("%.0f", nadir),
                  sprintf("%+.1f", best), resp,
                  ifelse(is.na(ttp), "not reached", sprintf("%.0f", ttp)),
                  sprintf("%.1f%%", 100 * tail(d$Survival_fraction, 1)))),
      rownames = FALSE, options = list(dom = "t"))
  })

  output$p_clones <- renderPlot({
    sim() %>%
      select(time, Sensitive = Tumor_sensitive_mm, Resistant = Tumor_resistant_mm) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, fill = name)) +
      geom_area(alpha = 0.8) +
      scale_fill_manual(values = c(Sensitive = "#7b4f9d", Resistant = "#c0392b")) +
      labs(x = "Time (days)", y = "Tumour burden (mm)") + theme_emc()
  })

  output$p_angio <- renderPlot({
    sim() %>%
      select(time, `VEGF-A` = VEGFA_level, `Microvessel density` = Microvessel_density) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      labs(x = "Time (days)", y = "Normalised level") + theme_emc()
  })

  ## ---------------- Tab 5: immune ----------------
  output$p_tcell <- renderPlot({
    sim() %>%
      select(time, Effector = CD8_effector, Exhausted = CD8_exhausted) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c(Effector = "#1e8449", Exhausted = "#c0392b")) +
      labs(x = "Time (days)", y = "Normalised density") + theme_emc()
  })
  output$p_supp <- renderPlot({
    sim() %>%
      select(time, Treg = Treg_density, MDSC = MDSC_density) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      labs(x = "Time (days)", y = "Normalised density") + theme_emc()
  })

  ## ---------------- Tab 6: endpoints ----------------
  output$p_surv <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, 100 * Survival_fraction)) +
      geom_line(linewidth = 1.2, colour = "#1e8449") +
      ylim(0, 100) +
      labs(x = "Time (days)", y = "Progression-free surrogate (%)") + theme_emc()
  })
  output$p_bio <- renderPlot({
    sim() %>%
      select(time, `CA-125 (U/mL)` = CA125_UmL, `ctDNA (%VAF)` = ctDNA_VAF) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Time (days)", y = NULL) + theme_emc()
  })
  output$p_anc <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, ANC_10e9L)) +
      geom_hline(yintercept = 1.0, linetype = 2, colour = "#c0392b") +
      geom_hline(yintercept = 0.5, linetype = 3, colour = "#c0392b") +
      geom_line(linewidth = 1.1, colour = "#2471a3") +
      annotate("text", x = max(d$time) * 0.02, y = 1.06, hjust = 0, size = 3.4,
               colour = "#c0392b", label = "grade 3 threshold") +
      labs(x = "Time (days)", y = "ANC (10^9/L)") + theme_emc()
  })

  ## ---------------- Tab 7: comparison ----------------
  cmp <- reactive({
    mod <- get_model()
    p <- base_params()
    if (isTRUE(input$cmp_allclass)) {
      bind_rows(lapply(1:4, function(cl) {
        pp <- p; pp$MOLCLASS <- cl
        run_sim(mod, input$scenario, pp, input$horizon,
                input$cycles, input$lendose) %>%
          mutate(arm = CLASS_SHORT[cl])
      }))
    } else {
      req(length(input$cmp_sc) > 0)
      bind_rows(lapply(input$cmp_sc, function(sc) {
        run_sim(mod, sc, p, input$horizon, input$cycles, input$lendose) %>%
          mutate(arm = sc)
      }))
    }
  })

  output$p_cmp_tum <- renderPlot({
    ggplot(cmp(), aes(time, Tumor_total_mm, colour = arm)) +
      geom_line(linewidth = 1.05) +
      labs(x = "Time (days)", y = "Sum of target lesions (mm)") +
      theme_emc() + guides(colour = guide_legend(ncol = 1))
  })
  output$p_cmp_surv <- renderPlot({
    ggplot(cmp(), aes(time, 100 * Survival_fraction, colour = arm)) +
      geom_line(linewidth = 1.05) + ylim(0, 100) +
      labs(x = "Time (days)", y = "Progression-free surrogate (%)") +
      theme_emc() + guides(colour = guide_legend(ncol = 1))
  })
  output$p_cmp_teff <- renderPlot({
    ggplot(cmp(), aes(time, CD8_effector, colour = arm)) +
      geom_line(linewidth = 1.05) +
      labs(x = "Time (days)", y = "CD8 effector density") +
      theme_emc() + guides(colour = guide_legend(ncol = 1))
  })
  output$t_cmp <- renderDT({
    cmp() %>%
      group_by(arm) %>%
      summarise(
        `TMB` = round(first(TMB_mutMb), 1),
        `Nadir (mm)` = round(min(Tumor_total_mm), 1),
        `Best change (%)` = round(min(Pct_change_baseline), 1),
        `ANC nadir` = round(min(ANC_10e9L), 2),
        `Survival at horizon (%)` = round(100 * last(Survival_fraction), 1),
        .groups = "drop") %>%
      datatable(rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
  })

  ## ---------------- Tab 8: TMB sweep ----------------
  sweep <- eventReactive(input$run_sweep, {
    mod <- get_model()
    p <- base_params()
    grid <- unique(round(10^seq(log10(input$tmb_lo), log10(input$tmb_hi),
                                length.out = 16), 2))
    withProgress(message = "Sweeping TMB...", value = 0, {
      bind_rows(lapply(seq_along(grid), function(i) {
        incProgress(1 / length(grid))
        pp <- p; pp$TMB_OVR <- grid[i]
        run_sim(mod, input$scenario, pp, input$horizon,
                input$cycles, input$lendose) %>%
          mutate(TMB = grid[i])
      }))
    })
  })

  sweep_sum <- reactive({
    sweep() %>%
      group_by(TMB) %>%
      summarise(`Immunogenicity` = first(Immunogenicity),
                `Best change (%)` = min(Pct_change_baseline),
                `Nadir (mm)` = min(Tumor_total_mm),
                `Survival at horizon (%)` = 100 * last(Survival_fraction),
                .groups = "drop")
  })

  output$p_sweep_resp <- renderPlot({
    tmb50 <- get_model()@param@data$TMB50
    ggplot(sweep_sum(), aes(TMB, `Best change (%)`)) +
      geom_vline(xintercept = tmb50, linetype = 2, colour = "grey40") +
      geom_hline(yintercept = -30, linetype = 3, colour = "#1e8449") +
      geom_line(linewidth = 1.2, colour = "#7b4f9d") +
      geom_point(size = 2.4, colour = "#7b4f9d") +
      scale_x_log10() +
      labs(x = "TMB (mut/Mb, log scale)", y = "Best change from baseline (%)",
           title = "Same drug, same patient — only TMB moves") +
      theme_emc()
  })
  output$p_sweep_surv <- renderPlot({
    tmb50 <- get_model()@param@data$TMB50
    ggplot(sweep_sum(), aes(TMB, `Survival at horizon (%)`)) +
      geom_vline(xintercept = tmb50, linetype = 2, colour = "grey40") +
      geom_line(linewidth = 1.2, colour = "#1e8449") +
      geom_point(size = 2.4, colour = "#1e8449") +
      scale_x_log10() + ylim(0, 100) +
      labs(x = "TMB (mut/Mb, log scale)", y = "Progression-free surrogate (%)") +
      theme_emc()
  })
  output$t_sweep <- renderDT({
    sweep_sum() %>%
      mutate(across(where(is.numeric), ~round(.x, 2))) %>%
      datatable(rownames = FALSE, options = list(dom = "tp", pageLength = 16))
  })
}

shinyApp(ui, server)
