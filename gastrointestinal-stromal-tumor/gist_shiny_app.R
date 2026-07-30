## =============================================================================
##  GIST QSP model — Shiny dashboard
## =============================================================================
##  12 tabs over the model in gist_mrgsolve_model.R.  The organising idea is that
##  each tab isolates one of the three structural commitments so that the user can
##  try to break it:
##
##    [C1] the tumour is a POPULATION      -> tabs 3, 4, 8, 12
##    [C2] occupancy fast, killing slow    -> tabs 6, 9
##    [C3] exposure matters where the       -> tabs 2, 5, 10
##         genotype puts the EC50
##
##  Run:  shiny::runApp("gist_shiny_app.R")
##  Needs: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

MOD <- mread_cache("gist", "gist_mrgsolve_model.R")

GENOS <- c("KIT exon 11 (65-70%)" = 1, "KIT exon 9 (8-10%)" = 2,
           "PDGFRA D842V (5%)" = 3, "SDH-deficient / wild type (5-8%)" = 4)
AGENTS <- c("imatinib" = "IM", "sunitinib" = "SU", "regorafenib" = "RE",
            "ripretinib" = "RI", "avapritinib" = "AV")
DEFDOSE <- c(IM = 400, SU = 50, RE = 160, RI = 150, AV = 300)
DEFSCHED <- c(IM = "qd", SU = "4/2", RE = "3/1", RI = "qd", AV = "qd")
CLONE_LAB <- c(N_PRIM = "primary genotype only",
               N_ATP  = "exon 13/14 ATP-pocket",
               N_LOOP = "exon 17/18 activation-loop",
               N_BYP  = "KIT-independent bypass")

## ---- dosing event builders --------------------------------------------------
ev_qd <- function(cmt, amt, start = 0, days = 3000) {
  if (amt <= 0 || days <= 0) return(NULL)
  ev(time = start, amt = amt, ii = 1, addl = days - 1, cmt = cmt)
}
ev_cyc <- function(cmt, amt, on, off, start = 0, days = 3000) {
  if (amt <= 0 || days <= 0) return(NULL)
  ncyc <- ceiling(days / (on + off))
  Reduce(`+`, lapply(seq_len(ncyc) - 1, function(k)
    ev(time = start + k * (on + off), amt = amt, ii = 1, addl = on - 1, cmt = cmt)))
}
dose_events <- function(drug, dose, sched, start = 0, days = 3000) {
  cmt <- paste0("A_", drug)
  switch(sched,
         "qd"  = ev_qd(cmt, dose, start, days),
         "4/2" = ev_cyc(cmt, dose, 28, 14, start, days),
         "3/1" = ev_cyc(cmt, dose, 21,  7, start, days),
         stop("unknown schedule"))
}

## ---- endpoints -------------------------------------------------------------
## RECIST 1.1 progression relative to the nadir on the current line.
recist_pd <- function(d, t0 = 0, assess = 7) {
  d <- dplyr::filter(d, time >= t0)
  if (nrow(d) < 3) return(NA_real_)
  tt  <- seq(min(d$time) + assess, max(d$time), by = assess)
  s   <- stats::approx(d$time, d$SLD, tt)$y
  nad <- cummin(s)
  hit <- which(s >= 1.2 * nad & s >= nad + 5)
  if (!length(hit)) return(NA_real_)
  tt[hit[1]] - t0
}
best_resp <- function(d, t0 = 0) {
  d <- dplyr::filter(d, time >= t0)
  if (!nrow(d)) return(NA_real_)
  100 * (min(d$SLD) - d$SLD[1]) / d$SLD[1]
}
mo <- function(d) d / 30.4375

## Run a sequence of lines, switching at RECIST progression.
run_sequence <- function(pars, lines, end = 4000, setting = 1) {
  m <- do.call(param, c(list(MOD), pars, list(SETTING = setting)))
  t0 <- 0; evs <- NULL; res <- list()
  for (i in seq_along(lines)) {
    ln <- lines[[i]]
    e  <- NULL
    for (a in ln) {
      e2 <- dose_events(a$drug, a$dose, a$sched, start = t0, days = end - t0)
      e  <- if (is.null(e)) e2 else e + e2
    }
    evs_try <- if (is.null(evs)) e else evs + e
    out <- as.data.frame(mrgsim(m, if (is.null(evs_try)) ev(amt = 0, cmt = 1) else evs_try,
                                end = end, delta = 1))
    pd <- recist_pd(out, t0)
    res[[i]] <- list(t0 = t0, pfs = pd, best = best_resp(out, t0), out = out)
    if (is.na(pd)) { evs <- evs_try; break }
    ## truncate this line's dosing at progression and continue from there
    e_trunc <- NULL
    for (a in ln) {
      e2 <- dose_events(a$drug, a$dose, a$sched, start = t0, days = ceiling(pd))
      e_trunc <- if (is.null(e_trunc)) e2 else e_trunc + e2
    }
    evs <- if (is.null(evs)) e_trunc else evs + e_trunc
    t0 <- t0 + ceiling(pd)
  }
  m2  <- do.call(param, c(list(MOD), pars, list(SETTING = setting)))
  fin <- as.data.frame(mrgsim(m2, if (is.null(evs)) ev(amt = 0, cmt = 1) else evs,
                              end = end, delta = 1))
  list(sim = fin, lines = res)
}

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

## =============================================================================
## UI
## =============================================================================
ui <- fluidPage(
  titlePanel("Gastrointestinal stromal tumour — QSP explorer"),
  p(tags$b("The model's claim:"), " GIST behaves like a ", tags$i("population"),
    " of clones with different drug sensitivities, most of which are already ",
    "present at diagnosis; the drug is largely cytostatic on a quiescent ",
    "reservoir; and exposure only matters where the genotype puts the EC50. ",
    "Each tab is a place to try to break one of those claims."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient and tumour"),
      selectInput("geno", "Primary genotype", choices = GENOS, selected = 1),
      sliderInput("v0", "Baseline tumour volume (mL)", 20, 1500, 400, step = 20),
      sliderInput("turnover", "Pre-diagnostic divisions per surviving cell",
                  0, 250, 85, step = 5,
                  post = "  (sets how many resistant cells exist on day 1)"),
      hr(),
      h4("Pharmacokinetics"),
      sliderInput("clf", "Clearance multiplier", 0.4, 4, 1, step = 0.1),
      sliderInput("agpf", "AGP multiplier", 0.6, 2.5, 1, step = 0.1),
      sliderInput("indf", "CYP3A4 induction (rifampicin ~2.6)", 0.5, 3, 1, step = 0.1),
      checkboxInput("difix", "Hold dose intensity at 100% (toxicity feedback off)",
                    FALSE),
      hr(),
      h4("First line"),
      selectInput("d1", "Agent", choices = AGENTS, selected = "IM"),
      numericInput("dose1", "Dose (mg/day)", 400, min = 0, max = 1600, step = 50),
      selectInput("sched1", "Schedule", c("qd", "4/2", "3/1"), "qd"),
      checkboxInput("add1", "Add a second agent from day 0 (combination)", FALSE),
      conditionalPanel(
        "input.add1",
        selectInput("d1b", "Second agent", choices = AGENTS, selected = "RI"),
        numericInput("dose1b", "Dose (mg/day)", 100, min = 0, max = 400, step = 25)),
      hr(),
      h4("Later lines (switch at progression)"),
      selectInput("d2", "Second line", choices = c("none" = "", AGENTS), "SU"),
      selectInput("d3", "Third line", choices = c("none" = "", AGENTS), "RE"),
      selectInput("d4", "Fourth line",
                  choices = c("none / best supportive care" = "", AGENTS), "RI"),
      hr(),
      sliderInput("end", "Simulate to (days)", 400, 5000, 3000, step = 100),
      actionButton("go", "Run", class = "btn-primary btn-block")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",
        tabPanel("1 · Profile",       br(), DTOutput("profile"),
                 br(), verbatimTextOutput("profile_txt")),
        tabPanel("2 · PK",            br(), plotOutput("pk", height = 560),
                 helpText("Imatinib exposure falls over the first month as it induces its own CYP3A4 clearance, and the AGP that rises with tumour burden lowers the free fraction — both are in the model, and both move the apparent EC50 rather than the potency.")),
        tabPanel("3 · Target coverage", br(), plotOutput("cover", height = 340),
                 DTOutput("ec50"),
                 helpText("The relative KIT pathway output SIG per clone. A clone with SIG above about 0.37 grows. This table is where the INTRIGUE crossover comes from: sunitinib covers the ATP-pocket column, ripretinib and regorafenib cover the activation-loop column, and nothing covers both.")),
        tabPanel("4 · Clonal architecture", br(), plotOutput("clones", height = 560),
                 helpText("Resistant subclones are generated per cell DIVISION, so almost all of them exist before treatment starts. Set the pre-diagnostic division slider to 0 to see a tumour with no pre-existing resistance.")),
        tabPanel("5 · Response",      br(), plotOutput("resp", height = 560),
                 DTOutput("endpoints")),
        tabPanel("6 · Three time scales", br(), plotOutput("scales", height = 560),
                 helpText("PET reads signalling (hours), Ki67 reads the cycling fraction (days), CT reads a mass that is mostly non-viable tissue (months). A responding tumour that does not shrink and a shrinking tumour that is not cured are the same fact.")),
        tabPanel("7 · ctDNA",         br(), plotOutput("ctdna", height = 520),
                 verbatimTextOutput("lead")),
        tabPanel("8 · Second-line choice", br(), plotOutput("intrigue", height = 420),
                 DTOutput("intrigue_tab"),
                 helpText("Sunitinib versus ripretinib in a tumour whose secondary mutations are confined to one class. Observed in INTRIGUE ctDNA: exon 11+13/14 sunitinib 15.0 vs ripretinib 4.0 months; exon 11+17/18 ripretinib 14.2 vs sunitinib 1.5 months.")),
        tabPanel("9 · Stopping",      br(), plotOutput("stopplot", height = 480),
                 DTOutput("stoptab"),
                 helpText("BFR14: interruption versus continuation after 1, 3 and 5 years of imatinib. The reservoir is not killed by years of deep response; what changes is how deeply dormant it is and how much of the imaged mass is still viable.")),
        tabPanel("10 · Adjuvant",     br(), plotOutput("adj", height = 480),
                 helpText("Recurrence-free survival for 12, 36 and 60 months of adjuvant imatinib over a virtual population. The curves separate while the drug is given and converge afterwards, because the model has no cure term — only a delay.")),
        tabPanel("11 · Toxicity",     br(), plotOutput("tox", height = 620),
                 helpText("Tolerability is part of efficacy: a toxicity-driven dose reduction lowers coverage, which accelerates selection. Tick 'hold dose intensity' in the sidebar to separate the two.")),
        tabPanel("12 · Falsification", br(), plotOutput("fals", height = 420),
                 DTOutput("falstab"),
                 helpText("The fair competitor is one clone whose potency degrades under drug pressure, refitted to the same first-line PFS. Two tests are then decisive against it: the genotype crossover, which it cannot express at any parameter value, and fourth-line ripretinib. The imatinib-rechallenge benefit does NOT discriminate -- the drifting model produces it too, and larger."))
      )
    )
  )
)

## =============================================================================
## SERVER
## =============================================================================
server <- function(input, output, session) {

  pars <- reactive({
    list(GENO = as.numeric(input$geno), V0 = input$v0,
         TURNOVER = input$turnover, CLF = input$clf, AGPF = input$agpf,
         INDF = input$indf, DIFIX = as.numeric(input$difix))
  })

  lines <- reactive({
    l1 <- list(list(drug = input$d1, dose = input$dose1, sched = input$sched1))
    if (isTRUE(input$add1) && input$dose1b > 0)
      l1 <- c(l1, list(list(drug = input$d1b, dose = input$dose1b,
                            sched = DEFSCHED[[input$d1b]])))
    out <- list(l1)
    for (d in c(input$d2, input$d3, input$d4)) {
      if (identical(d, "")) {
        out <- c(out, list(list()))            # best supportive care
      } else {
        out <- c(out, list(list(list(drug = d, dose = DEFDOSE[[d]],
                                     sched = DEFSCHED[[d]]))))
      }
    }
    out
  })

  RUN <- eventReactive(input$go, {
    withProgress(message = "simulating", value = 0.4,
                 run_sequence(pars(), lines(), end = input$end))
  }, ignoreNULL = FALSE)

  ## ---- 1 · profile ---------------------------------------------------------
  output$profile <- renderDT({
    p <- pars()
    data.frame(Quantity = names(p), Value = unlist(p)) |>
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })
  output$profile_txt <- renderPrint({
    r <- RUN()
    cat("Lines of therapy and per-line endpoints\n")
    for (i in seq_along(r$lines)) {
      li <- r$lines[[i]]
      cat(sprintf("  line %d  start day %5.0f   PFS %8s   best SLD %+6.1f %%\n",
                  i, li$t0,
                  if (is.na(li$pfs)) "not reached" else sprintf("%.1f mo", mo(li$pfs)),
                  li$best))
    }
  })

  ## ---- 2 · PK --------------------------------------------------------------
  output$pk <- renderPlot({
    d <- RUN()$sim
    d |>
      select(time, C_IMAT, C_SUNI, C_REGO, C_RIPR, C_AVAP, ENZ, AGP, DI) |>
      pivot_longer(-time) |>
      ggplot(aes(time, value)) +
      geom_line(linewidth = 0.6, colour = "#2c6fb5") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "day", y = NULL,
           title = "Concentrations (ng/mL), enzyme pool, AGP and achieved dose intensity") +
      THEME
  })

  ## ---- 3 · coverage --------------------------------------------------------
  output$cover <- renderPlot({
    d <- RUN()$sim
    ## SIG per clone is reconstructed from the captured concentrations using the
    ## same competitive expression as the model.
    ec <- ec_table(as.numeric(input$geno))
    cc <- d |> select(time, IM = C_IMAT, SU = C_SUNI, RE = C_REGO,
                      RI = C_RIPR, AV = C_AVAP)
    sig <- lapply(seq_len(nrow(ec)), function(i) {
      a <- 1 / ((1 + cc$IM / ec$IM[i]) * (1 + cc$SU / ec$SU[i]) *
                (1 + cc$RE / ec$RE[i]) * (1 + cc$RI / ec$RI[i]) *
                (1 + cc$AV / ec$AV[i]))
      data.frame(time = cc$time, clone = ec$clone[i],
                 SIG = 1 - ec$kitdep[i] * (1 - a))
    }) |> bind_rows()
    ggplot(sig, aes(time, SIG, colour = clone)) +
      geom_hline(yintercept = 0.37, linetype = 2, colour = "grey40") +
      geom_line(linewidth = 0.7) +
      annotate("text", x = max(sig$time) * 0.6, y = 0.40, size = 3.4,
               label = "above this line the clone grows", colour = "grey30") +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "day", y = "relative KIT pathway output  SIG", colour = NULL,
           title = "Coverage per clone") +
      THEME
  })
  output$ec50 <- renderDT({
    ec_table(as.numeric(input$geno)) |>
      datatable(rownames = FALSE, options = list(dom = "t"),
                caption = "EC50 (ng/mL, total plasma). Rows are clones, columns are drugs.")
  })

  ## ---- 4 · clones ----------------------------------------------------------
  output$clones <- renderPlot({
    d <- RUN()$sim
    long <- d |>
      select(time, N_PRIM, N_ATP, N_LOOP, N_BYP) |>
      pivot_longer(-time, names_to = "clone", values_to = "cells") |>
      mutate(clone = factor(CLONE_LAB[clone], levels = CLONE_LAB))
    p1 <- ggplot(long, aes(time, pmax(cells, 1), colour = clone)) +
      geom_line(linewidth = 0.7) + scale_y_log10() +
      labs(x = NULL, y = "cells (log)", colour = NULL,
           title = "Clone burdens: the resistant clones are present on day 1") + THEME
    p2 <- ggplot(long, aes(time, cells, fill = clone)) +
      geom_area(position = "fill") +
      labs(x = "day", y = "fraction of viable cells", fill = NULL) + THEME
    gridExtra::grid.arrange(p1, p2, ncol = 1, heights = c(1.15, 1))
  })

  ## ---- 5 · response --------------------------------------------------------
  output$resp <- renderPlot({
    d <- RUN()$sim
    lin <- RUN()$lines
    sw <- vapply(lin, function(x) x$t0, numeric(1))
    p1 <- d |>
      select(time, SLD) |>
      ggplot(aes(time, SLD)) +
      geom_vline(xintercept = sw, linetype = 3, colour = "grey55") +
      geom_line(linewidth = 0.7, colour = "#b03030") +
      labs(x = NULL, y = "SLD (mm)", title = "Imaged size, and what it is made of") + THEME
    p2 <- d |>
      transmute(time, viable = V_VIAB, `non-viable` = VNEC) |>
      pivot_longer(-time) |>
      ggplot(aes(time, value, fill = name)) + geom_area() +
      labs(x = "day", y = "volume (mL)", fill = NULL) + THEME
    gridExtra::grid.arrange(p1, p2, ncol = 1)
  })
  output$endpoints <- renderDT({
    lin <- RUN()$lines
    data.frame(line = seq_along(lin),
               start_day = round(vapply(lin, function(x) x$t0, numeric(1))),
               PFS_months = round(mo(vapply(lin, function(x) x$pfs, numeric(1))), 1),
               best_SLD_pct = round(vapply(lin, function(x) x$best, numeric(1)), 1)) |>
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- 6 · three time scales ----------------------------------------------
  output$scales <- renderPlot({
    d <- RUN()$sim |> filter(time <= 180)
    d |>
      transmute(time,
                `PET SUV (hours)` = 100 * SUV / SUV[1],
                `Ki67 (days)` = 100 * KI67 / KI67[1],
                `SLD (months)` = 100 * SLD / SLD[1],
                `viable fraction of mass` = 100 * V_FRAC) |>
      pivot_longer(-time) |>
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.8) +
      labs(x = "day", y = "% of baseline", colour = NULL,
           title = "Occupancy, cell cycle and mass move on different clocks") + THEME
  })

  ## ---- 7 · ctDNA -----------------------------------------------------------
  output$ctdna <- renderPlot({
    d <- RUN()$sim
    pd <- RUN()$lines[[1]]$pfs
    ggplot(d, aes(time, 100 * VAF_R)) +
      geom_line(linewidth = 0.7, colour = "#7a5ba6") +
      geom_hline(yintercept = 1, linetype = 2) +
      { if (!is.na(pd)) geom_vline(xintercept = pd, linetype = 3, colour = "#b03030") } +
      scale_y_log10() +
      labs(x = "day", y = "resistance-mutation VAF (%)",
           title = "ctDNA crosses the assay threshold before RECIST calls progression") +
      THEME
  })
  output$lead <- renderPrint({
    d <- RUN()$sim; pd <- RUN()$lines[[1]]$pfs
    t1 <- d$time[which(d$VAF_R >= 0.01)[1]]
    cat(sprintf("VAF 1%% crossed on day %s; RECIST progression on day %s; lead time %s\n",
                ifelse(is.na(t1), "never", round(t1)),
                ifelse(is.na(pd), "not reached", round(pd)),
                ifelse(is.na(t1) || is.na(pd), "-",
                       sprintf("%.1f months", mo(pd - t1)))))
  })

  ## ---- 8 · the INTRIGUE crossover -----------------------------------------
  intrigue <- reactive({
    mu_t <- 5.0e-7
    grid <- expand.grid(
      subgroup = c("exon 11 + 13/14 only", "exon 11 + 17/18 only", "both classes"),
      second = c("sunitinib", "ripretinib"), stringsAsFactors = FALSE)
    res <- lapply(seq_len(nrow(grid)), function(i) {
      mu <- switch(grid$subgroup[i],
                   "exon 11 + 13/14 only" = list(MU_ATP = mu_t, MU_AL = 1e-11),
                   "exon 11 + 17/18 only" = list(MU_ATP = 1e-11, MU_AL = mu_t),
                   list(MU_ATP = 2e-7, MU_AL = 3e-7))
      d2 <- if (grid$second[i] == "sunitinib") "SU" else "RI"
      r <- run_sequence(c(pars()[c("V0", "TURNOVER", "CLF", "AGPF")],
                          list(GENO = 1), mu),
                        list(list(list(drug = "IM", dose = 400, sched = "qd")),
                             list(list(drug = d2, dose = DEFDOSE[[d2]],
                                       sched = DEFSCHED[[d2]]))),
                        end = 3400)
      data.frame(subgroup = grid$subgroup[i], drug = grid$second[i],
                 pfs1 = mo(r$lines[[1]]$pfs),
                 pfs2 = if (length(r$lines) > 1) mo(r$lines[[2]]$pfs) else NA)
    })
    bind_rows(res)
  })
  output$intrigue <- renderPlot({
    ggplot(intrigue(), aes(subgroup, pfs2, fill = drug)) +
      geom_col(position = position_dodge(0.7), width = 0.6) +
      labs(x = NULL, y = "second-line PFS (months)", fill = NULL,
           title = "The ranking of two drugs inverts between genotypes") + THEME
  })
  output$intrigue_tab <- renderDT({
    intrigue() |> mutate(across(where(is.numeric), ~round(.x, 1))) |>
      datatable(rownames = FALSE, options = list(dom = "t"),
                caption = "Observed (INTRIGUE ctDNA): 15.0 vs 4.0 and 1.5 vs 14.2 months.")
  })

  ## ---- 9 · stopping --------------------------------------------------------
  stopping <- reactive({
    pheno <- list(`1 year` = list(y = 1, p = list(TURNOVER = 85)),
                  `3 years` = list(y = 3, p = list(TURNOVER = 0)),
                  `5 years` = list(y = 5, p = list(TURNOVER = 0, KPMAX = 0.0225)))
    bind_rows(lapply(names(pheno), function(nm) {
      ph <- pheno[[nm]]; d <- ph$y * 365
      m <- do.call(param, c(list(MOD), ph$p, list(GENO = 1, V0 = input$v0)))
      s_stop <- as.data.frame(mrgsim(m, ev_qd("A_IM", 400, 0, d),
                                     end = d + 1500, delta = 1))
      s_go   <- as.data.frame(mrgsim(m, ev_qd("A_IM", 400, 0, d + 3000),
                                     end = d + 3000, delta = 1))
      data.frame(arm = nm,
                 interruption = mo(recist_pd(s_stop, d)),
                 continuation = mo(recist_pd(s_go, d)),
                 viable_fraction = s_stop$V_FRAC[which.min(abs(s_stop$time - d))])
    }))
  })
  output$stopplot <- renderPlot({
    stopping() |>
      select(arm, interruption, continuation) |>
      pivot_longer(-arm) |>
      ggplot(aes(arm, value, fill = name)) +
      geom_col(position = position_dodge(0.7), width = 0.6) +
      labs(x = "randomised after", y = "PFS from randomisation (months)", fill = NULL,
           title = "Stopping imatinib in a responder") + THEME
  })
  output$stoptab <- renderDT({
    stopping() |> mutate(across(where(is.numeric), ~round(.x, 2))) |>
      datatable(rownames = FALSE, options = list(dom = "t"),
                caption = "Observed (BFR14): 6.1 vs 27.8, 7.0 vs 67.0, 12.0 vs not reached months.")
  })

  ## ---- 10 · adjuvant -------------------------------------------------------
  output$adj <- renderPlot({
    set.seed(20260730)
    n <- 200
    pop <- data.frame(ID = seq_len(n),
                      N_MICRO = 10^rnorm(n, log10(2.2e6), 0.55),
                      TURNOVER = pmax(0, rnorm(n, 85, 45)),
                      KPMAX = rlnorm(n, log(0.029), 0.25))
    curves <- bind_rows(lapply(c(365, 1095, 1826), function(dur) {
      out <- as.data.frame(mrgsim(param(MOD, SETTING = 2, GENO = as.numeric(input$geno)),
                                  idata = pop, events = ev_qd("A_IM", 400, 0, dur),
                                  end = 3652, delta = 28))
      out |> group_by(time) |>
        summarise(RFS = mean(DETECT == 0), .groups = "drop") |>
        mutate(duration = paste0(round(dur / 30.4375), " months"))
    }))
    ggplot(curves, aes(time / 365, 100 * RFS, colour = duration)) +
      geom_step(linewidth = 0.8) +
      geom_vline(xintercept = c(1, 3, 5), linetype = 3, colour = "grey60") +
      labs(x = "years since resection", y = "recurrence-free (%)", colour = NULL,
           title = "Adjuvant imatinib delays recurrence; the curves converge") +
      THEME
  })

  ## ---- 11 · toxicity ------------------------------------------------------
  output$tox <- renderPlot({
    RUN()$sim |>
      transmute(time, ANC, `MAP (mmHg)` = MAP_TOT, HFSR, oedema = EDEM,
                TSH, `free T4` = FT4, `dose intensity` = DI,
                `cumulative toxicity` = TOXC) |>
      pivot_longer(-time) |>
      ggplot(aes(time, value)) +
      geom_line(linewidth = 0.6, colour = "#b06a6a") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "day", y = NULL, title = "Organ-system effects and the dose-intensity loop") +
      THEME
  })

  ## ---- 12 · falsification -------------------------------------------------
  output$fals <- renderPlot({
    lv <- c("INTRIGUE crossover", "4th-line ripretinib",
            "imatinib rechallenge", "escalation confined to exon 9")
    d <- data.frame(test = factor(lv, levels = lv),
                    ## 1 = reproduces the observation, 0 = does not
                    polyclonal   = c(1, 1, 1, 1),
                    single_clone = c(0, 0, 1, 1)) |>
      pivot_longer(-test)
    ggplot(d, aes(test, value, fill = name)) +
      geom_col(position = position_dodge(0.7), width = 0.6) +
      scale_y_continuous(breaks = c(0, 1), labels = c("no", "yes")) +
      labs(x = NULL, y = NULL, fill = NULL,
           title = "Which observations each structure can reproduce") +
      THEME + theme(axis.text.x = element_text(angle = 20, hjust = 1))
  })
  output$falstab <- renderDT({
    data.frame(
      quantity = c("first-line PFS (anchor, refitted)",
                   "second line, exon 11+13/14: sunitinib vs ripretinib",
                   "second line, exon 11+17/18: sunitinib vs ripretinib",
                   "fourth line ripretinib (INVICTUS, observed 6.3 mo)",
                   "fourth line: imatinib rechallenge vs supportive care",
                   "escalation gain, exon 9 / exon 11"),
      polyclonal = c("28.5 mo", "6.2 vs 2.8 mo", "2.5 vs 31.3 mo", "6.7 mo",
                     "2.5 vs 2.1 mo", "+279% / +3%"),
      single_clone = c("28.5 mo (refit)",
                       "2.5 vs 2.8 mo, identical in both subgroups",
                       "2.5 vs 2.8 mo, identical in both subgroups",
                       "2.1 mo -- DECISIVE against one clone",
                       "9.0 vs 1.8 mo -- does NOT discriminate",
                       "+290% / +27% -- discriminates only weakly")) |>
      datatable(rownames = FALSE, options = list(dom = "t"),
                caption = paste("Numbers from gist_python_twin.py, which runs the",
                                "identical system and asserts them."))
  })
}

## ---- EC50 table used by tab 3 ----------------------------------------------
ec_table <- function(geno) {
  prim <- switch(as.character(geno),
    "2" = c(IM = 1000, SU = 15, RE = 1200, RI = 230, AV = 180),
    "3" = c(IM = 20000, SU = 2800, RE = 30000, RI = 600, AV = 40),
    "4" = c(IM = 1e7, SU = 1e7, RE = 1e7, RI = 1e7, AV = 1e7),
          c(IM = 185, SU = 20, RE = 1400, RI = 210, AV = 170))
  kd0 <- if (as.character(geno) == "4") 0.15 else 1.0
  data.frame(
    clone = c("primary genotype", "exon 13/14 ATP-pocket",
              "exon 17/18 activation-loop", "KIT-independent bypass"),
    IM = c(prim[["IM"]], 4500, 6500, 1e7),
    SU = c(prim[["SU"]], 12, 380, 1e7),
    RE = c(prim[["RE"]], 9000, 950, 1e7),
    RI = c(prim[["RI"]], 2200, 165, 1e7),
    AV = c(prim[["AV"]], 2700, 130, 1e7),
    kitdep = c(kd0, 1, 1, 0.05))
}

shinyApp(ui, server)
