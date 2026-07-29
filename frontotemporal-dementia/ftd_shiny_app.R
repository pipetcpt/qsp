## ============================================================================
##  Frontotemporal Dementia (FTD) QSP Model — Shiny Dashboard
## ============================================================================
##  Companion to ftd_mrgsolve_model.R. Nine tabs, built around the one thing
##  this model exists to show: a drug can normalise its registrational
##  biomarker while doing nothing — or harm — to the pool that matters.
##
##  Tabs
##    1. Patient & Genotype   — genotype, modifiers, predicted onset age
##    2. Drug PK & Target     — antibody PK, sortilin occupancy, ASO depots
##    3. Progranulin (3 pools)— THE headline plot: plasma vs CSF vs LYSOSOMAL
##    4. The p Sweep          — the sign flip at p = f_CSF, and the ceiling
##    5. Molecular Pathology   — TDP-43, cryptic splicing, DPRs, tau
##    6. Neuroinflammation     — microglia, C1q, C3, synaptic pruning
##    7. Clinical Endpoints    — CDR+NACC-FTLD SB, NPI, FRS, survival
##    8. Biomarkers            — NfL, GFAP, BMP, poly-GP, volumetry
##    9. Scenario Comparison   — all 13 prebuilt scenarios side by side
##
##  Run:  shiny::runApp("ftd_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  Disclaimer: research / education only.
## ============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

## The model, helpers and scenarios all live in the mrgsolve file.
source("ftd_mrgsolve_model.R")

F_CSF_REF <- 0.78   # sortilin share of CSF PGRN clearance (the break-even p)

theme_ftd <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Frontotemporal Dementia — QSP Model Explorer"),
  tags$p(
    tags$em(
      "GRN / C9orf72 / MAPT -> TDP-43 or tau -> complement-mediated synaptic ",
      "pruning -> salience-network degeneration. Research and education only."
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Patient"),
      selectInput("genotype", "Genotype",
                  choices = c("Sporadic FTLD" = "sporadic",
                              "GRN (progranulin)" = "GRN",
                              "C9orf72 expansion" = "C9",
                              "MAPT (tau)" = "MAPT"),
                  selected = "GRN"),
      sliderInput("age0", "Age at simulation start (y)", 25, 55, 30, step = 1),
      checkboxInput("tmem", "Protective TMEM106B haplotype", FALSE),
      sliderInput("temp_shift", "Temporal shift (0 = bvFTD, 0.4 = svPPA)",
                  0, 0.6, 0, step = 0.05),
      checkboxInput("als", "FTD-ALS phenotype (adds mortality hazard)", FALSE),

      hr(),
      h4("Disease-modifying therapy"),
      checkboxInput("use_lato", "Latozinemab (anti-sortilin) IV q4w", TRUE),
      conditionalPanel(
        "input.use_lato",
        sliderInput("lato_mgkg", "Dose (mg/kg)", 5, 120, 60, step = 5)
      ),
      checkboxInput("use_aav", "AAV-GRN gene therapy (single dose)", FALSE),
      checkboxInput("use_aso9", "C9orf72 ASO, intrathecal q12w", FALSE),
      checkboxInput("use_asot", "MAPT ASO, intrathecal q12w", FALSE),

      hr(),
      h4("Symptomatic therapy"),
      checkboxInput("use_ssri", "SSRI (citalopram 30 mg/d)", FALSE),
      checkboxInput("use_trz", "Trazodone 150 mg/d", FALSE),
      checkboxInput("use_dnp", "Donepezil 10 mg/d (NULL / HARM arm)", FALSE),

      hr(),
      h4(tags$span(style = "color:#b00020",
                   "The decisive unmeasured parameter")),
      sliderInput("p_lys", "p = sortilin share of LYSOSOMAL delivery",
                  0, 0.95, 0.50, step = 0.05),
      helpText(
        HTML(paste0(
          "Break-even is <b>p = f<sub>CSF</sub> = 0.78</b>. Below it the drug ",
          "helps; above it the drug <b>raises the biomarker and lowers the ",
          "target pool</b>. No human measurement of p exists, and plasma PGRN ",
          "cannot measure it."
        ))
      ),
      sliderInput("w_act_lys", "Weight of LYSOSOMAL pool in PGRN action",
                  0, 1, 0.80, step = 0.05),
      helpText("Set near 0.2 to express the rival 'PGRN acts extracellularly'",
               "hypothesis: same drug, opposite prediction."),

      hr(),
      sliderInput("weeks", "Trial duration (weeks from onset)",
                  24, 208, 96, step = 4),
      actionButton("run", "Run simulation", class = "btn-primary",
                   width = "100%")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel(
          "1. Patient & Genotype",
          br(),
          fluidRow(
            column(5, h4("Predicted onset age"),
                   tableOutput("onset_tbl"),
                   helpText("Onset age is an OUTPUT of the age gate, not an",
                            "input. Observed values are cohort means from",
                            "Moore 2020 Lancet Neurol (PMID 31810826).")),
            column(7, h4("Prodromal trajectory to onset"),
                   plotOutput("prodrome_plot", height = "330px"))
          ),
          hr(),
          h4("Baseline state at the modelled trial entry"),
          tableOutput("baseline_tbl")
        ),

        tabPanel(
          "2. Drug PK & Target",
          br(),
          fluidRow(
            column(6, h5("Antibody PK (plasma vs CSF, log scale)"),
                   plotOutput("pk_plot", height = "300px")),
            column(6, h5("Sortilin occupancy, peripheral vs CNS"),
                   plotOutput("occ_plot", height = "300px"))
          ),
          hr(),
          fluidRow(
            column(6, h5("ASO CNS tissue depots"),
                   plotOutput("aso_plot", height = "280px")),
            column(6, h5("Symptomatic drug concentrations"),
                   plotOutput("sym_pk_plot", height = "280px"))
          )
        ),

        tabPanel(
          "3. Progranulin — three pools",
          br(),
          h4("The dissociation this model exists to show"),
          plotOutput("pgrn_plot", height = "400px"),
          helpText(
            "Plasma PGRN is the registrational biomarker. The LYSOSOMAL pool",
            "is where progranulin acts. Sortilin blockade raises the first and",
            "removes a delivery route into the third at the same time."
          ),
          hr(),
          fluidRow(
            column(6, h5("Fold change from baseline, by pool"),
                   tableOutput("pgrn_fold_tbl")),
            column(6, h5("Lysosomal function and CSF BMP"),
                   plotOutput("lyso_plot", height = "260px"))
          )
        ),

        tabPanel(
          "4. The p sweep",
          br(),
          h4("Sign flip at p = f_CSF = 0.78, at any dose"),
          plotOutput("sweep_plot", height = "380px"),
          helpText(
            HTML(paste0(
              "Note the plasma-PGRN trace: <b>flat across the whole sweep</b>. ",
              "The drug looks identical on its biomarker whether it is helping, ",
              "doing nothing, or making the target pool worse. Dose escalation ",
              "cannot fix this, because raising occupancy raises the ligand and ",
              "removes delivery capacity in lockstep. The ceiling as occupancy ",
              "-> 1 is (1-p)/(1-f<sub>CSF</sub>)."
            ))
          ),
          hr(),
          DTOutput("sweep_tbl")
        ),

        tabPanel(
          "5. Molecular Pathology",
          br(),
          fluidRow(
            column(6, h5("TDP-43: nuclear loss and aggregate burden"),
                   plotOutput("tdp_plot", height = "290px")),
            column(6, h5("Cryptic splicing: STMN2 and UNC13A"),
                   plotOutput("cryptic_plot", height = "290px"))
          ),
          hr(),
          fluidRow(
            column(6, h5("C9orf72: repeat RNA and dipeptide repeats"),
                   plotOutput("c9_plot", height = "290px")),
            column(6, h5("Tau: mRNA -> soluble -> phospho -> aggregate"),
                   plotOutput("tau_plot", height = "290px"))
          )
        ),

        tabPanel(
          "6. Neuroinflammation",
          br(),
          h4("Microglia -> C1q -> C3 -> synapse elimination"),
          plotOutput("inflam_plot", height = "340px"),
          helpText(
            "Anti-sortilin therapy lowers CSF C1q and C3 in this model, which",
            "matches what was observed clinically. Target engagement was never",
            "the problem."
          ),
          hr(),
          fluidRow(
            column(6, h5("Synaptic density and regional neuron pools"),
                   plotOutput("syn_plot", height = "280px")),
            column(6, h5("Regional grey-matter volume"),
                   plotOutput("vol_plot", height = "280px"))
          )
        ),

        tabPanel(
          "7. Clinical Endpoints",
          br(),
          fluidRow(
            column(6, h5("CDR plus NACC-FTLD sum of boxes (0-24)"),
                   plotOutput("cdr_plot", height = "300px")),
            column(6, h5("NPI total (0-144)"),
                   plotOutput("npi_plot", height = "300px"))
          ),
          hr(),
          fluidRow(
            column(6, h5("Model survival probability"),
                   plotOutput("surv_plot", height = "280px")),
            column(6, h5("Treatment effect vs untreated reference"),
                   tableOutput("effect_tbl"),
                   helpText("A fraction-of-one-percent CDR effect over 96",
                            "weeks is below trial resolution — this is the",
                            "model predicting a phase 3 miss."))
          )
        ),

        tabPanel(
          "8. Biomarkers",
          br(),
          fluidRow(
            column(6, h5("Neurofilament light (CSF and plasma)"),
                   plotOutput("nfl_plot", height = "290px")),
            column(6, h5("CSF BMP — the LYSOSOMAL readout"),
                   plotOutput("bmp_plot", height = "290px"))
          ),
          helpText(
            "CSF BMP tracks the lysosomal pool rather than the plasma pool, so",
            "it separates cases that plasma PGRN cannot. It is the measurement",
            "this model argues should sit alongside plasma PGRN in any",
            "progranulin-directed trial."
          ),
          hr(),
          fluidRow(
            column(6, h5("Plasma GFAP"),
                   plotOutput("gfap_plot", height = "260px")),
            column(6, h5("CSF poly-GP (C9orf72 target engagement)"),
                   plotOutput("gp_plot", height = "260px"))
          )
        ),

        tabPanel(
          "9. Scenario Comparison",
          br(),
          h4("All 13 prebuilt scenarios"),
          actionButton("run_all", "Run all scenarios (slow)",
                       class = "btn-warning"),
          br(), br(),
          plotOutput("scen_plot", height = "420px"),
          hr(),
          DTOutput("scen_tbl")
        )
      )
    )
  )
)

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  ## --- assemble the parameter set from the sidebar -------------------------
  pars <- reactive({
    p <- FTD_genotype(input$genotype)
    p$AGE0 <- input$age0
    p$TMEM106B_PROT <- as.numeric(input$tmem)
    p$TEMP_SHIFT <- input$temp_shift
    p$PHENO_TEMP <- as.numeric(input$temp_shift > 0.2)
    p$ALS_FLAG <- as.numeric(input$als || identical(input$genotype, "C9"))
    p$P_LYS_SORT <- input$p_lys
    p$W_ACT_LYS <- input$w_act_lys
    p$W_ACT_EXTRA <- 1 - input$w_act_lys
    p
  })

  ## --- burn-in to the model's own onset state ------------------------------
  burn <- eventReactive(input$run, {
    withProgress(message = "Solving the disease clock (burn-in)...", {
      FTD_burnin(pars())
    })
  }, ignoreNULL = FALSE)

  ## --- dosing events ------------------------------------------------------
  events <- reactive({
    wk <- input$weeks
    days <- wk * 7
    ev_list <- list()
    if (input$use_lato) {
      ev_list <- c(ev_list, list(
        ev_latozinemab(input$lato_mgkg, 70, n = ceiling(days / 28))))
    }
    if (input$use_aav)  ev_list <- c(ev_list, list(ev_aav_grn(1.0)))
    if (input$use_aso9) ev_list <- c(ev_list, list(
      ev_aso_c9(100, n = max(1, ceiling(days / 84)))))
    if (input$use_asot) ev_list <- c(ev_list, list(
      ev_aso_mapt(100, n = max(1, ceiling(days / 84)))))
    if (input$use_ssri) ev_list <- c(ev_list, list(ev_ssri(30, n = days)))
    if (input$use_trz)  ev_list <- c(ev_list, list(ev_trazodone(150, n = days)))
    if (input$use_dnp)  ev_list <- c(ev_list, list(ev_donepezil(10, n = days)))
    if (!length(ev_list)) return(NULL)
    Reduce(c, ev_list)
  })

  ## --- treated and untreated arms from the SAME baseline state ------------
  sim <- eventReactive(input$run, {
    bi <- burn()
    withProgress(message = "Simulating treated and untreated arms...", {
      treated <- FTD_from_onset(pars(), bi, events(), weeks = input$weeks,
                                delta = 1) |> mutate(arm = "Treated")
      control <- FTD_from_onset(pars(), bi, NULL, weeks = input$weeks,
                                delta = 1) |> mutate(arm = "Untreated")
    })
    list(both = bind_rows(control, treated), burn = bi)
  }, ignoreNULL = FALSE)

  both <- reactive(sim()$both)

  ## helper: tidy long-format plot of selected columns
  longp <- function(cols, labels = cols, ylab = "", title = NULL,
                    logy = FALSE) {
    d <- both() |>
      select(time, arm, all_of(cols)) |>
      pivot_longer(all_of(cols), names_to = "var", values_to = "value") |>
      mutate(var = factor(var, levels = cols, labels = labels),
             weeks = time / 7)
    g <- ggplot(d, aes(weeks, value, colour = var, linetype = arm)) +
      geom_line(linewidth = 0.8) +
      labs(x = "Weeks from trial entry", y = ylab, colour = NULL,
           linetype = NULL, title = title) +
      theme_ftd
    if (logy) g <- g + scale_y_log10()
    g
  }

  ## ===================== Tab 1: patient & genotype ======================
  output$onset_tbl <- renderTable({
    obs <- c(sporadic = 58.0, GRN = 61.3, C9 = 58.2, MAPT = 49.5)
    bi <- burn()
    data.frame(
      Quantity = c("Genotype", "Predicted onset age (y)",
                   "Observed cohort mean (y)", "Error (y)",
                   "Years simulated to onset"),
      Value = c(input$genotype,
                sprintf("%.1f", bi$onset_age),
                sprintf("%.1f", obs[[input$genotype]]),
                sprintf("%+.1f", bi$onset_age - obs[[input$genotype]]),
                sprintf("%.1f", bi$time / 365))
    )
  }, colnames = FALSE)

  output$prodrome_plot <- renderPlot({
    tr <- burn()$trajectory |>
      mutate(age = pars()$AGE0 + time / 365)
    tr |>
      select(age, CDR, NFLPL, VOLF, PGRNLYS_PCT) |>
      pivot_longer(-age) |>
      mutate(name = recode(name,
                           CDR = "CDR+NACC-FTLD SB",
                           NFLPL = "plasma NfL (pg/mL)",
                           VOLF = "frontal volume (%)",
                           PGRNLYS_PCT = "lysosomal PGRN (%WT)")) |>
      ggplot(aes(age, value)) +
      geom_line(linewidth = 0.9, colour = "#2c6fbb") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Age (years)", y = NULL,
           title = "Prodrome from a clinically silent carrier to onset") +
      theme_ftd
  })

  output$baseline_tbl <- renderTable({
    b <- both() |> filter(arm == "Untreated") |> slice(1)
    data.frame(
      Measure = c("CDR+NACC-FTLD SB", "NPI total", "plasma NfL (pg/mL)",
                  "CSF NfL (pg/mL)", "frontal volume (% baseline)",
                  "synaptic density", "lysosomal function",
                  "plasma PGRN (% control)", "CSF PGRN (% control)",
                  "lysosomal PGRN (% WT)"),
      Value = sprintf("%.2f", c(b$CDR, b$NPI, b$NFLPL, b$NFLCSF, b$VOLF,
                                b$SYN, b$LYSO, b$PGRNPL_PCT, b$PGRNCSF_PCT,
                                b$PGRNLYS_PCT))
    )
  })

  ## ===================== Tab 2: PK & target =============================
  output$pk_plot <- renderPlot({
    d <- both() |> filter(arm == "Treated") |>
      select(time, CLATO, CSFLATO) |>
      pivot_longer(-time) |>
      mutate(name = recode(name, CLATO = "plasma (nM)",
                           CSFLATO = "CSF (nM)"),
             value = pmax(value, 1e-4), weeks = time / 7)
    ggplot(d, aes(weeks, value, colour = name)) +
      geom_line(linewidth = 0.8) + scale_y_log10() +
      labs(x = "Weeks", y = "Concentration (nM, log)", colour = NULL) +
      theme_ftd
  })

  output$occ_plot <- renderPlot({
    longp(c("SORT_OCC_PL", "SORT_OCC_CNS"),
          c("peripheral occupancy", "CNS occupancy"),
          "Fractional sortilin occupancy") +
      geom_hline(yintercept = F_CSF_REF, linetype = "dotted") +
      coord_cartesian(ylim = c(0, 1))
  })

  output$aso_plot <- renderPlot({
    longp(c("ASO9TIS", "ASOTTIS"),
          c("C9orf72 ASO depot", "MAPT ASO depot"), "Tissue depot (au)")
  })

  output$sym_pk_plot <- renderPlot({
    longp(c("SSRIC", "TRZC", "DNPC"),
          c("SSRI", "trazodone", "donepezil"), "mg/L")
  })

  ## ===================== Tab 3: progranulin =============================
  output$pgrn_plot <- renderPlot({
    longp(c("PGRNPL_PCT", "PGRNCSF_PCT", "PGRNLYS_PCT"),
          c("PLASMA (the biomarker)", "CSF", "LYSOSOMAL (where it acts)"),
          "% of wild-type level for that pool",
          "Plasma PGRN normalises; the lysosomal pool barely moves") +
      geom_hline(yintercept = 100, linetype = "dashed", colour = "grey40")
  })

  output$pgrn_fold_tbl <- renderTable({
    d <- both()
    f <- function(arm_, col) {
      x <- d[[col]][d$arm == arm_]
      c(first = x[1], last = tail(x, 1))
    }
    cols <- c("PGRNPL", "PGRNCSF", "PGRNLYS")
    labs <- c("Plasma PGRN", "CSF PGRN", "LYSOSOMAL PGRN")
    out <- lapply(seq_along(cols), function(i) {
      u <- f("Untreated", cols[i]); t <- f("Treated", cols[i])
      data.frame(Pool = labs[i],
                 Baseline = sprintf("%.3g", u[["first"]]),
                 Untreated_end = sprintf("%.3g", u[["last"]]),
                 Treated_end = sprintf("%.3g", t[["last"]]),
                 Fold = sprintf("%.2fx", t[["last"]] / u[["last"]]))
    })
    bind_rows(out)
  })

  output$lyso_plot <- renderPlot({
    longp(c("LYSO", "BMPCSF"),
          c("lysosomal function (0-1)", "CSF BMP (au)"), "")
  })

  ## ===================== Tab 4: the p sweep =============================
  sweep <- eventReactive(input$run, {
    grid <- seq(0, 0.95, by = 0.1)
    withProgress(message = "Sweeping p (this runs a burn-in per value)...", {
      bind_rows(lapply(grid, function(pv) {
        incProgress(1 / length(grid))
        pp <- pars(); pp$P_LYS_SORT <- pv
        bi <- FTD_burnin(pp)
        n <- ceiling(input$weeks * 7 / 28)
        tx <- FTD_from_onset(pp, bi, ev_latozinemab(input$lato_mgkg, 70, n = n),
                             weeks = input$weeks)
        ct <- FTD_from_onset(pp, bi, NULL, weeks = input$weeks)
        tibble(
          p = pv,
          plasma_pct = tail(tx$PGRNPL_PCT, 1),
          lyso_pct = tail(tx$PGRNLYS_PCT, 1),
          lyso_ratio = tail(tx$PGRNLYS, 1) / tail(ct$PGRNLYS, 1),
          ceiling = (1 - pv) / (1 - F_CSF_REF),
          dCDR = tail(tx$CDR, 1) - tail(ct$CDR, 1),
          dBMP = tail(tx$BMPCSF, 1) - tail(ct$BMPCSF, 1)
        )
      }))
    })
  }, ignoreNULL = FALSE)

  output$sweep_plot <- renderPlot({
    s <- sweep()
    d <- s |>
      select(p, `plasma PGRN (% control)` = plasma_pct,
             `lysosomal PGRN (% WT)` = lyso_pct,
             `lysosomal delivery ratio` = lyso_ratio,
             `dCDR at end of trial` = dCDR) |>
      pivot_longer(-p)
    ggplot(d, aes(p, value)) +
      geom_line(linewidth = 0.9, colour = "#b00020") +
      geom_point(size = 1.6, colour = "#b00020") +
      geom_vline(xintercept = F_CSF_REF, linetype = "dashed") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "p  =  sortilin share of lysosomal PGRN delivery", y = NULL,
           title = "Dashed line: break-even p = f_CSF = 0.78") +
      theme_ftd
  })

  output$sweep_tbl <- renderDT({
    sweep() |>
      mutate(verdict = case_when(lyso_ratio > 1.02 ~ "BENEFIT",
                                 lyso_ratio < 0.98 ~ "HARM",
                                 TRUE ~ "NEUTRAL")) |>
      mutate(across(where(is.numeric), ~round(.x, 4))) |>
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 12))
  })

  ## ===================== Tab 5: molecular ===============================
  output$tdp_plot <- renderPlot(longp(c("TDPN", "TDPAGG"),
    c("nuclear TDP-43 (0-1)", "TDP-43 aggregate (au)"), ""))
  output$cryptic_plot <- renderPlot(longp(c("STMN2T", "UNC13AT"),
    c("truncated STMN2", "UNC13A cryptic"), "au"))
  output$c9_plot <- renderPlot(longp(c("C9RNA", "POLYGP", "POLYGR"),
    c("repeat RNA", "poly-GP", "poly-GR"), "au"))
  output$tau_plot <- renderPlot(longp(c("TAUM", "TAUS", "TAUP", "TAUAGG"),
    c("MAPT mRNA", "soluble tau", "phospho tau", "aggregated tau"), "au"))

  ## ===================== Tab 6: neuroinflammation =======================
  output$inflam_plot <- renderPlot(longp(c("MG", "C1Q", "C3", "AST"),
    c("microglial activation", "C1q", "C3", "reactive astrocytes"), ""))
  output$syn_plot <- renderPlot(longp(c("SYN", "NEURSN", "NEURTMP"),
    c("synaptic density", "salience-network neurons",
      "anterior-temporal neurons"), "fraction of baseline"))
  output$vol_plot <- renderPlot(longp(c("VOLF", "VOLT", "VOL_WB"),
    c("frontal", "temporal", "whole brain"), "% of baseline"))

  ## ===================== Tab 7: clinical ================================
  output$cdr_plot <- renderPlot({
    longp("CDR", "CDR+NACC-FTLD SB", "Sum of boxes (0-24)") +
      coord_cartesian(ylim = c(0, 24))
  })
  output$npi_plot <- renderPlot(longp("NPI", "NPI total", "0-144"))
  output$surv_plot <- renderPlot(longp("SURVIVAL", "survival probability", ""))

  output$effect_tbl <- renderTable({
    d <- both()
    last_of <- function(a, col) tail(d[[col]][d$arm == a], 1)
    first_of <- function(col) d[[col]][d$arm == "Untreated"][1]
    rows <- c("CDR", "NPI", "NFLPL", "VOLF", "SYN", "BMPCSF", "C3")
    labs <- c("CDR+NACC-FTLD SB", "NPI total", "plasma NfL", "frontal volume %",
              "synaptic density", "CSF BMP", "CSF C3")
    prog <- last_of("Untreated", "CDR") - first_of("CDR")
    tab <- data.frame(
      Endpoint = labs,
      Untreated = sprintf("%.3f", sapply(rows, last_of, a = "Untreated")),
      Treated = sprintf("%.3f", sapply(rows, last_of, a = "Treated")),
      Difference = sprintf("%+.3f", sapply(rows, last_of, a = "Treated") -
                             sapply(rows, last_of, a = "Untreated"))
    )
    dcdr <- last_of("Treated", "CDR") - last_of("Untreated", "CDR")
    rbind(tab, data.frame(
      Endpoint = "-- CDR slowing vs untreated progression --",
      Untreated = sprintf("%+.2f pts", prog),
      Treated = sprintf("%+.4f pts", dcdr),
      Difference = sprintf("%.2f%%", -100 * dcdr / prog)))
  })

  ## ===================== Tab 8: biomarkers ==============================
  output$nfl_plot <- renderPlot(longp(c("NFLCSF", "NFLPL"),
    c("CSF NfL", "plasma NfL"), "pg/mL", logy = TRUE))
  output$bmp_plot <- renderPlot(longp("BMPCSF", "CSF BMP", "au"))
  output$gfap_plot <- renderPlot(longp("GFAPPL", "plasma GFAP", "pg/mL"))
  output$gp_plot <- renderPlot(longp("POLYGP", "CSF poly-GP", "au"))

  ## ===================== Tab 9: scenarios ===============================
  scen <- eventReactive(input$run_all, {
    withProgress(message = "Running all 13 scenarios...", {
      FTD_simulate_scenarios(weeks = input$weeks)
    })
  })

  output$scen_plot <- renderPlot({
    s <- scen()
    s |>
      filter(grepl("^(0[6-9]|1[0-3])", scenario)) |>
      select(time, scenario, CDR, PGRNLYS_PCT, NFLPL, NPI) |>
      pivot_longer(c(CDR, PGRNLYS_PCT, NFLPL, NPI)) |>
      ggplot(aes(time / 7, value, colour = scenario)) +
      geom_line(linewidth = 0.75) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Weeks from onset state", y = NULL, colour = NULL) +
      theme_ftd +
      theme(legend.text = element_text(size = 7))
  })

  output$scen_tbl <- renderDT({
    scen() |>
      group_by(scenario) |>
      slice_tail(n = 1) |>
      ungroup() |>
      select(scenario, CDR, NPI, NFLPL, PGRNPL_PCT, PGRNCSF_PCT,
             PGRNLYS_PCT, VOLF, BMPCSF, SURVIVAL) |>
      mutate(across(where(is.numeric), ~round(.x, 3))) |>
      datatable(rownames = FALSE,
                options = list(pageLength = 13, scrollX = TRUE))
  })
}

shinyApp(ui, server)
