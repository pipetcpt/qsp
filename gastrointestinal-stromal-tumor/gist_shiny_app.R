##  ============================================================================
##  GIST QSP Model — Shiny Dashboard
##  ============================================================================
##  Companion interactive front-end for gist_mrgsolve_model.R.
##
##  Run with:
##     library(shiny); runApp("gist_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, tidyr, purrr, ggplot2, DT
##
##  Design intent
##  -------------
##  The clinical questions this dashboard exists to answer are all questions
##  about CLONAL COMPOSITION, not about tumor size:
##     * Which clone is actually driving this patient's disease right now?
##     * Is the resistant fraction already rising while the scan still looks
##       like a response?
##     * Given this patient's clonal make-up, which second-line drug wins?
##     * Is the trough adequate, and does fixing it buy anything?
##  Every tab is built around one of those.
##
##  8 tabs:
##     1. Patient & genotype      — set up the case, see the starting clone mix
##     2. Drug exposure (PK)      — concentrations, imatinib trough vs 1100 ng/mL
##     3. Target engagement       — the drug x clone inhibition matrix, live
##     4. Clonal dynamics         — the six sub-populations over time
##     5. Imaging endpoints       — RECIST vs Choi density vs FDG-PET SUV
##     6. ctDNA & molecular lead  — resistant-clone VAF vs radiologic PD
##     7. Scenario comparison     — line sequencing and the INTRIGUE question
##     8. Toxicity & tolerability — edema, ANC, TSH, BP, hand-foot
##  ============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(DT)

mod <- mread_cache("gist_mrgsolve_model.R")

GENOTYPES <- c(
  "KIT exon 11 (juxtamembrane, ~65-70%)"            = 1,
  "KIT exon 9 (extracellular, ~8-10%)"              = 2,
  "PDGFRA exon 18 D842V (~5-6%)"                    = 3,
  "SDH-deficient / quadruple wild-type"             = 4,
  "KIT exon 11 + BIM intron-2 deletion polymorphism" = 5
)

CLONE_LABELS <- c(
  Clone_sensitive  = "Primary-driver sensitive",
  Clone_quiescent  = "Quiescent persister",
  Clone_exon13_14  = "Exon 13/14 (ATP pocket)",
  Clone_exon17_18  = "Exon 17/18 (activation loop)",
  Clone_refractory = "Primary-refractory (D842V / SDH)",
  Clone_KIT_indep  = "KIT-independent bypass"
)

CLONE_COLS <- c(
  "Primary-driver sensitive"         = "#2c7fb8",
  "Quiescent persister"              = "#7fcdbb",
  "Exon 13/14 (ATP pocket)"          = "#fdae61",
  "Exon 17/18 (activation loop)"     = "#d7301f",
  "Primary-refractory (D842V / SDH)" = "#762a83",
  "KIT-independent bypass"           = "#525252"
)

## ---------------------------------------------------------------------------
## Dosing helpers — all oral once daily; sunitinib 4/2, regorafenib 3/1
## ---------------------------------------------------------------------------
ev_daily <- function(cmt, mg, start, days) {
  if (days <= 0 || mg <= 0) return(NULL)
  as.data.frame(ev(amt = mg, cmt = cmt, time = start, ii = 1, addl = days - 1))
}
ev_cyclic <- function(cmt, mg, start, on_days, cycle_days, n_cycles) {
  if (n_cycles <= 0 || mg <= 0) return(NULL)
  map_dfr(seq_len(n_cycles) - 1, function(k)
    as.data.frame(ev(amt = mg, cmt = cmt,
                     time = start + k * cycle_days, ii = 1, addl = on_days - 1)))
}
bind_ev <- function(...) {
  parts <- Filter(Negate(is.null), list(...))
  if (!length(parts)) return(ev(amt = 0, cmt = "IMA_DEPOT", time = 0))
  bind_rows(parts) %>% arrange(time)
}

run_sim <- function(pars, events, end_d, delta = 1) {
  m <- mod
  if (length(pars)) m <- do.call(param, c(list(m), pars))
  mrgsim(m, data = events, end = end_d, delta = delta) %>% as_tibble()
}

## Radiologic progression: >=20% RELATIVE increase in diameter over the nadir
pd_time <- function(df) {
  d <- 100 + df$RECIST_pct_change
  nad <- cummin(d)
  i <- which(d / nad >= 1.20)
  if (length(i)) df$time[i[1]] else NA_real_
}
first_time <- function(df, col, thresh) {
  i <- which(df[[col]] >= thresh)
  if (length(i)) df$time[i[1]] else NA_real_
}

theme_qsp <- function() {
  theme_bw(base_size = 12) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92"))
}

## ===========================================================================
## UI
## ===========================================================================
ui <- fluidPage(
  titlePanel("GIST QSP Model — polyclonal resistance, TKI sequencing and response assessment"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      selectInput("geno", "Driver genotype", choices = GENOTYPES, selected = 1),
      sliderInput("tb0", "Baseline tumor burden (mL)", 50, 2000, 500, step = 50),
      checkboxInput("adjuvant", "Adjuvant setting (R0 resected, micrometastatic)", FALSE),
      sliderInput("clfac", "Imatinib clearance multiplier (CYP3A4 / adherence)",
                  0.6, 3.0, 1.0, step = 0.1),
      hr(),
      h4("Clonal background"),
      helpText("Pre-existing resistant subclone fractions at presentation. ",
               "Setting one to zero models a patient whose tumor does not ",
               "generate that resistance family — the ctDNA-defined INTRIGUE subgroups."),
      numericInput("fa0", "Exon 13/14 (ATP pocket) fraction", 3.0e-4,
                   min = 0, max = 0.05, step = 1e-4),
      numericInput("fl0", "Exon 17/18 (activation loop) fraction", 2.4e-4,
                   min = 0, max = 0.05, step = 1e-4),
      hr(),
      h4("Therapy"),
      numericInput("ima_mg", "Imatinib mg/day", 400, min = 0, max = 800, step = 100),
      numericInput("ima_start", "Imatinib start (day)", 0, min = 0, step = 30),
      numericInput("ima_days", "Imatinib duration (days)", 1460, min = 0, step = 30),
      numericInput("sun_mg", "Sunitinib mg/day (4 wk on / 2 off)", 0, min = 0, max = 50, step = 12.5),
      numericInput("sun_start", "Sunitinib start (day)", 700, min = 0, step = 30),
      numericInput("reg_mg", "Regorafenib mg/day (3 wk on / 1 off)", 0, min = 0, max = 160, step = 40),
      numericInput("reg_start", "Regorafenib start (day)", 950, min = 0, step = 30),
      numericInput("rip_mg", "Ripretinib mg/day", 0, min = 0, max = 150, step = 50),
      numericInput("rip_start", "Ripretinib start (day)", 1120, min = 0, step = 30),
      numericInput("ava_mg", "Avapritinib mg/day", 0, min = 0, max = 300, step = 100),
      numericInput("ava_start", "Avapritinib start (day)", 0, min = 0, step = 30),
      hr(),
      sliderInput("end_d", "Simulation horizon (days)", 180, 2555, 1460, step = 30),
      actionButton("go", "Run simulation", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel(
          "1. Patient profile",
          br(),
          fluidRow(
            column(6, h4("Baseline clonal composition"), plotOutput("p_compo0", height = 300)),
            column(6, h4("Genotype implications"), htmlOutput("txt_geno"))
          ),
          hr(),
          h4("Effective IC50 matrix in force for this patient (nM)"),
          helpText("Rows = drug, columns = tumor sub-population. Lower = more potent. ",
                   "These are calibrated effective potencies, not raw biochemical Kd."),
          DTOutput("tbl_ic50")
        ),

        tabPanel(
          "2. Drug exposure (PK)",
          br(),
          plotOutput("p_pk", height = 380),
          hr(),
          h4("Imatinib trough vs the 1,100 ng/mL threshold"),
          helpText("Demetri 2009 (JCO 27:3141): patients in the lowest steady-state ",
                   "trough quartile had markedly shorter time to progression."),
          plotOutput("p_cmin", height = 300),
          hr(),
          h4("Alpha-1-acid glycoprotein and free-drug availability"),
          helpText("AAG is an acute-phase protein induced by tumor burden. A large ",
                   "tumor sequesters imatinib, so free drug is LOWER at the same ",
                   "total plasma concentration — and rises again as the tumor responds."),
          plotOutput("p_aag", height = 260)
        ),

        tabPanel(
          "3. Target engagement",
          br(),
          h4("Pathway activity by sub-population"),
          helpText("1.0 = fully active KIT/PDGFRA signaling; 0 = complete shutdown. ",
                   "Deep inhibition is cytocidal, partial inhibition is only cytostatic ",
                   "(Hill threshold at IC_AP = 0.60)."),
          plotOutput("p_engage", height = 340),
          hr(),
          h4("Downstream signaling and apoptotic priming"),
          plotOutput("p_signal", height = 320)
        ),

        tabPanel(
          "4. Clonal dynamics",
          br(),
          h4("Sub-population volumes (log scale)"),
          plotOutput("p_clones", height = 380),
          hr(),
          h4("Clonal composition over time (stacked fraction)"),
          helpText("The clinically important quantity: not how big the tumor is, ",
                   "but what it is made of."),
          plotOutput("p_compo", height = 320),
          hr(),
          h4("Milestones"),
          DTOutput("tbl_milestones")
        ),

        tabPanel(
          "5. Imaging endpoints",
          br(),
          h4("Three response readouts on one clock"),
          helpText("FDG-PET SUV falls within 24-48 h, CT density (Choi) over weeks, ",
                   "RECIST diameter over months. This ordering is why RECIST ",
                   "systematically under-calls GIST response."),
          plotOutput("p_imaging", height = 400),
          hr(),
          fluidRow(
            column(6, h4("Choi vs RECIST partial response"), tableOutput("tbl_choi")),
            column(6, h4("Tumor burden"), plotOutput("p_burden", height = 280))
          )
        ),

        tabPanel(
          "6. ctDNA & molecular lead time",
          br(),
          h4("Resistant-clone ctDNA vs radiologic progression"),
          helpText("Resistant-clone variant allele fraction becomes detectable while ",
                   "the scan still shows response. The gap between the two vertical ",
                   "lines is the window in which a treatment change could be made ",
                   "before radiologic progression."),
          plotOutput("p_ctdna", height = 400),
          hr(),
          tableOutput("tbl_lead")
        ),

        tabPanel(
          "7. Scenario comparison",
          br(),
          h4("Prebuilt comparisons"),
          selectInput("cmp", NULL, width = "70%", choices = c(
            "First line: imatinib 400 vs 800 mg (current genotype)" = "dose",
            "KIT exon 9: 400 vs 800 mg (MetaGIST)"                  = "ex9",
            "Trough: normal vs high clearance vs TDM escalation"    = "tdm",
            "Continue vs stop imatinib at 1 year (BFR14)"           = "stop",
            "Second line after imatinib: sunitinib vs ripretinib (INTRIGUE)" = "intrigue",
            "PDGFRA D842V: imatinib vs avapritinib"                 = "d842v",
            "SDH-deficient: imatinib vs sunitinib"                  = "sdh",
            "Adjuvant 3 years vs observation (SSGXVIII)"            = "adj"
          )),
          actionButton("go_cmp", "Run comparison", class = "btn-primary"),
          br(), br(),
          plotOutput("p_cmp", height = 400),
          hr(),
          h4("Endpoint summary"),
          DTOutput("tbl_cmp"),
          br(),
          htmlOutput("txt_cmp")
        ),

        tabPanel(
          "8. Toxicity",
          br(),
          h4("Organ-system effects"),
          plotOutput("p_tox", height = 460),
          hr(),
          h4("Worst on-treatment values"),
          tableOutput("tbl_tox"),
          helpText("Edema is imatinib-dominant (PDGFR-beta); hypothyroidism is ",
                   "essentially sunitinib-specific; hypertension tracks VEGFR ",
                   "blockade; hand-foot skin reaction is weighted to regorafenib.")
        )
      )
    )
  )
)

## ===========================================================================
## SERVER
## ===========================================================================
server <- function(input, output, session) {

  base_pars <- reactive({
    list(GENO = as.numeric(input$geno), TB0 = input$tb0,
         ADJUVANT = as.numeric(input$adjuvant), CLFAC_IMA = input$clfac,
         FRAC_A0 = input$fa0, FRAC_L0 = input$fl0)
  })

  base_events <- reactive({
    e <- input$end_d
    bind_ev(
      ev_daily("IMA_DEPOT", input$ima_mg, input$ima_start,
               min(input$ima_days, e - input$ima_start)),
      ev_cyclic("SUN_DEPOT", input$sun_mg, input$sun_start, 28, 42,
                max(0, floor((e - input$sun_start) / 42))),
      ev_cyclic("REG_DEPOT", input$reg_mg, input$reg_start, 21, 28,
                max(0, floor((e - input$reg_start) / 28))),
      ev_daily("RIP_DEPOT", input$rip_mg, input$rip_start, e - input$rip_start),
      ev_daily("AVA_DEPOT", input$ava_mg, input$ava_start, e - input$ava_start)
    )
  })

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    run_sim(base_pars(), base_events(), input$end_d)
  })

  clone_long <- reactive({
    sim() %>%
      select(time, all_of(names(CLONE_LABELS))) %>%
      pivot_longer(-time, names_to = "clone", values_to = "vol") %>%
      mutate(clone = factor(CLONE_LABELS[clone], levels = CLONE_LABELS))
  })

  ## ---- Tab 1 ---------------------------------------------------------------
  output$p_compo0 <- renderPlot({
    d <- clone_long() %>% filter(time == min(time))
    ggplot(d, aes("", vol, fill = clone)) +
      geom_col(width = 0.6) +
      scale_fill_manual(values = CLONE_COLS) +
      scale_y_continuous(trans = "sqrt") +
      labs(x = NULL, y = "Volume (mL, sqrt scale)") +
      theme_qsp() + guides(fill = guide_legend(ncol = 1))
  })

  output$txt_geno <- renderUI({
    g <- as.numeric(input$geno)
    msg <- switch(as.character(g),
      "1" = "<b>KIT exon 11.</b> The most imatinib-sensitive genotype. Expect a deep
             response, a Choi partial response within weeks, and secondary
             resistance driven by exon 13/14 and exon 17/18 mutants at roughly
             20-24 months. Standard dose is 400 mg/day.",
      "2" = "<b>KIT exon 9.</b> Imatinib IC50 is ~5-fold higher, so 400 mg gives only
             partial pathway suppression and a shallow response. 800 mg/day is the
             standard of care and the model reproduces the response-depth gain.
             Sunitinib is relatively <i>more</i> potent on exon 9 than on exon 11.",
      "3" = "<b>PDGFRA exon 18 D842V.</b> Imatinib is inactive — this genotype must be
             identified before therapy, not after failure. Avapritinib, a type I
             inhibitor that binds the active conformation, is the treatment.",
      "4" = "<b>SDH-deficient / quadruple wild-type.</b> KIT-independent. Succinate
             accumulation inhibits prolyl hydroxylases, stabilizes HIF-1&alpha; and drives
             VEGF — so this tumor's dependence is <i>vascular</i>, not kinase.
             KIT-directed therapy is inert; sunitinib slows growth through its
             anti-angiogenic arm rather than by target inhibition.",
      "5" = "<b>KIT exon 11 with the BIM intron-2 2,903-bp deletion polymorphism</b>
             (~12% of East Asians). Exposure and target inhibition are normal, but
             apoptotic competence is impaired, so the response is shallower at the
             same pathway shutdown — cytostasis without cytocidal effect."
    )
    HTML(paste0("<div style='padding:8px'>", msg, "</div>"))
  })

  output$tbl_ic50 <- renderDT({
    p <- as.list(param(mod))
    f9  <- if (as.numeric(input$geno) == 2) p$GF_EX9 else 1
    f9s <- if (as.numeric(input$geno) == 2) p$GF_EX9_SUN else 1
    avaD <- if (as.numeric(input$geno) == 4) p$IC_AVA_D_SDH else p$IC_AVA_D
    m <- data.frame(
      Drug = c("Imatinib", "Sunitinib", "Regorafenib", "Ripretinib", "Avapritinib"),
      `Primary driver`   = c(p$IC_IMA_S*f9, p$IC_SUN_S*f9s, p$IC_REG_S, p$IC_RIP_S, p$IC_AVA_S),
      `Exon 13/14`       = c(p$IC_IMA_A*f9, p$IC_SUN_A*f9s, p$IC_REG_A, p$IC_RIP_A, p$IC_AVA_A),
      `Exon 17/18`       = c(p$IC_IMA_L*f9, p$IC_SUN_L*f9s, p$IC_REG_L, p$IC_RIP_L, p$IC_AVA_L),
      `Refractory clone` = c(p$IC_IMA_D, p$IC_SUN_D, p$IC_REG_D, p$IC_RIP_D, avaD),
      check.names = FALSE
    )
    datatable(m, rownames = FALSE, options = list(dom = "t", ordering = FALSE)) %>%
      formatRound(2:5, 1)
  })

  ## ---- Tab 2 ---------------------------------------------------------------
  output$p_pk <- renderPlot({
    sim() %>%
      select(time, Imatinib = Imatinib_ng_mL, Sunitinib = Sunitinib_ng_mL,
             Regorafenib = Regorafenib_ng_mL, Ripretinib = Ripretinib_ng_mL,
             Avapritinib = Avapritinib_ng_mL) %>%
      pivot_longer(-time) %>% filter(max(value) > 0, .by = name) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.5) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Day", y = "Plasma concentration (ng/mL)") +
      theme_qsp() + theme(legend.position = "none")
  })

  output$p_cmin <- renderPlot({
    d <- sim() %>% filter(Imatinib_ng_mL > 0)
    if (!nrow(d)) return(NULL)
    trough <- d %>% mutate(day = floor(time)) %>%
      summarise(Cmin = min(Imatinib_ng_mL), .by = day)
    ggplot(trough, aes(day, Cmin)) +
      geom_line(colour = "#2c7fb8", linewidth = 0.7) +
      geom_hline(yintercept = 1100, linetype = "dashed", colour = "#d7301f") +
      annotate("text", x = max(trough$day) * 0.75, y = 1160,
               label = "1,100 ng/mL efficacy threshold", colour = "#d7301f") +
      labs(x = "Day", y = "Imatinib daily trough (ng/mL)") +
      theme_qsp()
  })

  output$p_aag <- renderPlot({
    sim() %>% select(time, `AAG (relative)` = AAG_level,
                     `Tumor burden (mL)` = Tumor_burden_mL) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#762a83", linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "Day", y = NULL) + theme_qsp()
  })

  ## ---- Tab 3 ---------------------------------------------------------------
  output$p_engage <- renderPlot({
    sim() %>% select(time, `Aggregate pKIT` = pKIT_activity,
                     `Vascular support` = Vascular_support) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 0.40, linetype = "dotted") +
      annotate("text", x = 0, y = 0.44, hjust = 0, size = 3.2,
               label = "residual activity 0.40 = inhibition 0.60 = cytocidal threshold") +
      scale_colour_manual(values = c("#2c7fb8", "#41ab5d")) +
      ylim(0, 1.05) + labs(x = "Day", y = "Relative activity") + theme_qsp()
  })

  output$p_signal <- renderPlot({
    sim() %>% select(time, AKT = AKT_activity, ERK = ERK_activity,
                     ETV1 = ETV1_level, BIM = BIM_level) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Day", y = "Relative level") +
      theme_qsp() + theme(legend.position = "none")
  })

  ## ---- Tab 4 ---------------------------------------------------------------
  output$p_clones <- renderPlot({
    ggplot(clone_long(), aes(time, pmax(vol, 1e-4), colour = clone)) +
      geom_line(linewidth = 0.8) +
      scale_y_log10() + scale_colour_manual(values = CLONE_COLS) +
      labs(x = "Day", y = "Sub-population volume (mL, log scale)") +
      theme_qsp()
  })

  output$p_compo <- renderPlot({
    ggplot(clone_long(), aes(time, vol, fill = clone)) +
      geom_area(position = "fill") +
      scale_fill_manual(values = CLONE_COLS) +
      labs(x = "Day", y = "Fraction of tumor burden") + theme_qsp()
  })

  output$tbl_milestones <- renderDT({
    d <- sim()
    ms <- tibble(
      Milestone = c("Best RECIST response (%)", "Day of nadir",
                    "Choi partial response (day)", "RECIST partial response (day)",
                    "Resistant ctDNA detectable (day)",
                    "Radiologic progression (day)",
                    "Resistant fraction at end of run (%)"),
      Value = c(
        round(min(d$RECIST_pct_change), 1),
        d$time[which.min(d$RECIST_pct_change)],
        first_time(d, "Choi_PR", 1), first_time(d, "RECIST_PR", 1),
        first_time(d, "ctDNA_res_detected", 1),
        pd_time(d),
        round(100 * tail(d$Resistant_frac, 1), 1)
      )
    )
    datatable(ms, rownames = FALSE, options = list(dom = "t", ordering = FALSE))
  })

  ## ---- Tab 5 ---------------------------------------------------------------
  output$p_imaging <- renderPlot({
    sim() %>%
      select(time,
             `RECIST (% diameter change)` = RECIST_pct_change,
             `Choi CT density (% change)` = Choi_density_pct,
             `FDG-PET SUVmax (% change)`  = SUV_pct_change) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = c(-30, -15), linetype = c("dashed", "dotted")) +
      scale_colour_manual(values = c("#d7301f", "#2c7fb8", "#41ab5d")) +
      labs(x = "Day", y = "% change from baseline",
           caption = "dashed = RECIST PR (-30%); dotted = Choi density PR (-15%)") +
      theme_qsp()
  })

  output$tbl_choi <- renderTable({
    d <- sim()
    data.frame(
      Criterion = c("Choi PR (density -15% or size -10%)", "RECIST PR (-30% diameter)"),
      `First met (day)` = c(first_time(d, "Choi_PR", 1), first_time(d, "RECIST_PR", 1)),
      check.names = FALSE
    )
  })

  output$p_burden <- renderPlot({
    ggplot(sim(), aes(time, Tumor_burden_mL)) +
      geom_line(colour = "#252525", linewidth = 0.8) +
      scale_y_log10() + labs(x = "Day", y = "Total burden (mL, log)") + theme_qsp()
  })

  ## ---- Tab 6 ---------------------------------------------------------------
  output$p_ctdna <- renderPlot({
    d <- sim()
    tpd  <- pd_time(d)
    tmol <- first_time(d, "ctDNA_res_detected", 1)
    d %>% select(time, `Total ctDNA VAF (%)` = ctDNA_VAF_total,
                 `Resistant-clone VAF (%)` = ctDNA_VAF_resistant) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, pmax(value, 1e-3), colour = name)) +
      geom_line(linewidth = 0.9) +
      { if (!is.na(tmol)) geom_vline(xintercept = tmol, linetype = "dashed",
                                     colour = "#2c7fb8") } +
      { if (!is.na(tpd))  geom_vline(xintercept = tpd, linetype = "dashed",
                                     colour = "#d7301f") } +
      scale_y_log10() +
      scale_colour_manual(values = c("#762a83", "#d7301f")) +
      labs(x = "Day", y = "Variant allele fraction (%, log scale)",
           caption = "blue = resistant clone first detectable; red = radiologic progression") +
      theme_qsp()
  })

  output$tbl_lead <- renderTable({
    d <- sim(); tpd <- pd_time(d); tmol <- first_time(d, "ctDNA_res_detected", 1)
    data.frame(
      `Resistant ctDNA detectable (day)` = tmol,
      `Radiologic progression (day)`     = tpd,
      `Molecular lead time (days)`       = if (is.na(tpd) || is.na(tmol)) NA else tpd - tmol,
      check.names = FALSE
    )
  })

  ## ---- Tab 7 ---------------------------------------------------------------
  cmp <- eventReactive(input$go_cmp, {
    e <- input$end_d
    bp <- base_pars()
    arms <- switch(input$cmp,
      "dose" = list(
        "Imatinib 400 mg" = list(bp, bind_ev(ev_daily("IMA_DEPOT", 400, 0, e))),
        "Imatinib 800 mg" = list(bp, bind_ev(ev_daily("IMA_DEPOT", 800, 0, e)))),
      "ex9" = list(
        "Exon 9, imatinib 400 mg" = list(modifyList(bp, list(GENO = 2)),
                                         bind_ev(ev_daily("IMA_DEPOT", 400, 0, e))),
        "Exon 9, imatinib 800 mg" = list(modifyList(bp, list(GENO = 2)),
                                         bind_ev(ev_daily("IMA_DEPOT", 800, 0, e)))),
      "tdm" = list(
        "Normal clearance, 400 mg"  = list(modifyList(bp, list(CLFAC_IMA = 1.0)),
                                           bind_ev(ev_daily("IMA_DEPOT", 400, 0, e))),
        "High clearance, 400 mg"    = list(modifyList(bp, list(CLFAC_IMA = 2.2)),
                                           bind_ev(ev_daily("IMA_DEPOT", 400, 0, e))),
        "High clearance, TDM 800 mg"= list(modifyList(bp, list(CLFAC_IMA = 2.2)),
                                           bind_ev(ev_daily("IMA_DEPOT", 800, 0, e)))),
      "stop" = list(
        "Continue imatinib"        = list(bp, bind_ev(ev_daily("IMA_DEPOT", 400, 0, e))),
        "Stop imatinib at day 365" = list(bp, bind_ev(ev_daily("IMA_DEPOT", 400, 0, 365)))),
      "intrigue" = list(
        "2L sunitinib, exon 13/14 only" =
          list(modifyList(bp, list(FRAC_A0 = 5.4e-4, FRAC_L0 = 0, MU_L = 0)),
               bind_ev(ev_daily("IMA_DEPOT", 400, 0, 640),
                       ev_cyclic("SUN_DEPOT", 50, 640, 28, 42, floor((e-640)/42)))),
        "2L ripretinib, exon 13/14 only" =
          list(modifyList(bp, list(FRAC_A0 = 5.4e-4, FRAC_L0 = 0, MU_L = 0)),
               bind_ev(ev_daily("IMA_DEPOT", 400, 0, 640),
                       ev_daily("RIP_DEPOT", 150, 640, e - 640))),
        "2L sunitinib, exon 17/18 only" =
          list(modifyList(bp, list(FRAC_A0 = 0, FRAC_L0 = 5.4e-4, MU_A = 0)),
               bind_ev(ev_daily("IMA_DEPOT", 400, 0, 640),
                       ev_cyclic("SUN_DEPOT", 50, 640, 28, 42, floor((e-640)/42)))),
        "2L ripretinib, exon 17/18 only" =
          list(modifyList(bp, list(FRAC_A0 = 0, FRAC_L0 = 5.4e-4, MU_A = 0)),
               bind_ev(ev_daily("IMA_DEPOT", 400, 0, 640),
                       ev_daily("RIP_DEPOT", 150, 640, e - 640)))),
      "d842v" = list(
        "D842V, imatinib 400 mg"    = list(modifyList(bp, list(GENO = 3)),
                                           bind_ev(ev_daily("IMA_DEPOT", 400, 0, e))),
        "D842V, avapritinib 300 mg" = list(modifyList(bp, list(GENO = 3)),
                                           bind_ev(ev_daily("AVA_DEPOT", 300, 0, e)))),
      "sdh" = list(
        "SDH-deficient, imatinib"        = list(modifyList(bp, list(GENO = 4)),
                                                bind_ev(ev_daily("IMA_DEPOT", 400, 0, e))),
        "SDH-deficient, sunitinib 37.5"  = list(modifyList(bp, list(GENO = 4)),
                                                bind_ev(ev_daily("SUN_DEPOT", 37.5, 0, e)))),
      "adj" = list(
        "Resected, observation"      = list(modifyList(bp, list(ADJUVANT = 1)),
                                            bind_ev(NULL)),
        "Resected, adjuvant 3 years" = list(modifyList(bp, list(ADJUVANT = 1)),
                                            bind_ev(ev_daily("IMA_DEPOT", 400, 0, min(1095, e)))))
    )
    imap_dfr(arms, function(a, lab)
      run_sim(a[[1]], a[[2]], e) %>% mutate(arm = lab))
  })

  output$p_cmp <- renderPlot({
    ggplot(cmp(), aes(time, Tumor_burden_mL, colour = arm)) +
      geom_line(linewidth = 0.9) + scale_y_log10() +
      labs(x = "Day", y = "Tumor burden (mL, log scale)") +
      theme_qsp() + guides(colour = guide_legend(ncol = 2))
  })

  output$tbl_cmp <- renderDT({
    cmp() %>% group_by(arm) %>% group_modify(function(d, k) {
      tibble(
        `Best RECIST (%)`       = round(min(d$RECIST_pct_change), 1),
        `Nadir day`             = d$time[which.min(d$RECIST_pct_change)],
        `Choi PR (day)`         = first_time(d, "Choi_PR", 1),
        `Radiologic PD (day)`   = pd_time(d),
        `Burden at end (mL)`    = round(tail(d$Tumor_burden_mL, 1), 1),
        `Resistant frac (%)`    = round(100 * tail(d$Resistant_frac, 1), 1)
      )
    }) %>% ungroup() %>%
      datatable(rownames = FALSE, options = list(dom = "t", ordering = FALSE))
  })

  output$txt_cmp <- renderUI({
    msg <- switch(input$cmp,
      "intrigue" = "<b>Reading this comparison.</b> The two clonal backgrounds are the
         ctDNA-defined subgroups of INTRIGUE. Sunitinib is potent on ATP-binding-pocket
         (exon 13/14) mutants and nearly inert on activation-loop (exon 17/18) mutants;
         ripretinib is good on both but not the best on either. The model therefore
         predicts sunitinib wins the exon 13/14 subgroup and ripretinib wins the
         exon 17/18 subgroup by a much larger margin — matching the direction of the
         published ctDNA analysis. Note the model does <i>not</i> reproduce the overall
         equipoise of the INTRIGUE primary endpoint; see the README limitations section.",
      "stop"     = "<b>BFR14.</b> Regrowth after interruption comes out of the quiescent
         persister pool, which is why it is fast and why the stopped arm does not
         return to the continuous-therapy curve on re-treatment.",
      "adj"      = "<b>SSGXVIII.</b> Adjuvant imatinib delays recurrence by roughly the
         duration of therapy rather than eradicating micrometastatic disease — the
         'shift, not cure' shape of the published curves.",
      "ex9"      = "<b>MetaGIST.</b> The dominant modeled effect of 800 mg in exon 9 is
         on response depth; the progression-delay component is present but smaller
         than the published PFS separation.",
      ""
    )
    if (nzchar(msg)) HTML(paste0("<div style='padding:8px;background:#f6f6f6'>", msg, "</div>"))
  })

  ## ---- Tab 8 ---------------------------------------------------------------
  output$p_tox <- renderPlot({
    sim() %>%
      select(time, `Edema score (0-4)` = Edema_score,
             `Neutrophils (10^9/L)` = Neutrophils,
             `TSH (mIU/L)` = TSH_level,
             `Systolic BP (mmHg)` = Systolic_BP,
             `Hand-foot score (0-3)` = HFSR_score) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#d7301f", linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "Day", y = NULL) + theme_qsp()
  })

  output$tbl_tox <- renderTable({
    d <- sim()
    data.frame(
      Endpoint = c("Peak edema score", "Nadir ANC (10^9/L)", "Days with ANC < 1.0",
                   "Peak TSH (mIU/L)", "Peak systolic BP (mmHg)",
                   "Days with SBP >= 160", "Peak hand-foot score"),
      Value = c(round(max(d$Edema_score), 2), round(min(d$Neutrophils), 2),
                sum(d$Neutropenia_G34 > 0), round(max(d$TSH_level), 2),
                round(max(d$Systolic_BP), 0), sum(d$Hypertension_G3 > 0),
                round(max(d$HFSR_score), 2))
    )
  })
}

shinyApp(ui, server)
