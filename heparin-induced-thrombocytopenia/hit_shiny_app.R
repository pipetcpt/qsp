# =====================================================================
# Heparin-Induced Thrombocytopenia (HIT) — Shiny dashboard
# Companion to hit_mrgsolve_model.R (59-ODE QSP model)
#
#   10 tabs:
#     1  Patient & exposure     — who, on what heparin, at what dose
#     2  Antigen (the bell curve) — PF4, chains, molar ratio, ULC
#     3  Heparin & drug PK      — every anticoagulant on one axis
#     4  Immunity              — B cells, plasmablasts, IgG, seroreversion
#     5  Platelet axis         — count, nadir, 4Ts, recovery
#     6  Thrombin axis         — TF, microparticles, thrombin, D-dimer
#     7  Clinical endpoints    — thrombosis and venous limb gangrene
#     8  Diagnostics           — ELISA OD, SRA at low / high / no heparin
#     9  The two traps         — warfarin mismatch, argatroban-INR
#    10  Scenario comparison   — any number of arms side by side
#
# Run:  shiny::runApp("hit_shiny_app.R")
# Requires: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
#
# EDUCATIONAL / RESEARCH USE ONLY. Not for clinical decision-making.
# =====================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

# ---------------------------------------------------------------------
# Load the model by sourcing the spec out of the model file, so the app
# and the model can never drift apart.
# ---------------------------------------------------------------------
MODEL_FILE <- "hit_mrgsolve_model.R"
stopifnot(file.exists(MODEL_FILE))
src  <- readLines(MODEL_FILE)
i0   <- grep("^hit_code <- '", src)[1]
i1   <- grep("^'$", src)
i1   <- i1[i1 > i0][1]
spec <- paste(src[(i0 + 1):(i1 - 1)], collapse = "\n")
mod  <- mcode_cache("hit_shiny", spec)

WT_DEF <- 70

theme_hit <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "grey35", size = 10),
          legend.position = "bottom", legend.title = element_blank(),
          strip.text = element_text(face = "bold"))
}
PAL <- c("#C62828", "#1565C0", "#2E7D32", "#E65100", "#6A1B9A",
         "#00838F", "#AD1457", "#37474F", "#827717", "#4527A0")

# ---------------------------------------------------------------------
# Build an event object from the UI state
# ---------------------------------------------------------------------
build_events <- function(inp) {
  wt <- inp$wt
  ev_list <- list()
  add <- function(e) ev_list[[length(ev_list) + 1]] <<- e

  # ---- index heparin exposure ----
  dur <- inp$hep_days * 24
  if (inp$hep_type == "UFH intravenous") {
    add(ev(amt = inp$ufh_bolus * wt, cmt = "HEP"))
    add(ev(amt = inp$ufh_rate * wt * dur, rate = inp$ufh_rate * wt, cmt = "HEP"))
  } else if (inp$hep_type == "UFH subcutaneous") {
    add(ev(amt = inp$ufh_sc_dose, cmt = "HEPD", ii = inp$ufh_sc_ii,
           addl = max(0, floor(dur/inp$ufh_sc_ii) - 1)))
  } else if (inp$hep_type == "LMWH (enoxaparin)") {
    add(ev(amt = inp$lmw_mg * 100, cmt = "LMWD", ii = inp$lmw_ii,
           addl = max(0, floor(dur/inp$lmw_ii) - 1)))
  } else if (inp$hep_type == "Fondaparinux") {
    add(ev(amt = inp$fon_mg * 1000, cmt = "FOND", ii = 24,
           addl = max(0, inp$hep_days - 1)))
  } else if (inp$hep_type == "Danaparoid") {
    add(ev(amt = inp$dan_u, cmt = "DANC", ii = 8,
           addl = max(0, floor(dur/8) - 1)))
  } else if (inp$hep_type == "Cardiopulmonary bypass") {
    add(ev(amt = inp$cpb_ukg * wt, cmt = "HEP"))
    add(ev(amt = 1, cmt = "SURG"))
    if (inp$protamine) add(ev(amt = 300, cmt = "PROT", time = inp$cpb_h))
    if (inp$hep_days > 1)
      add(ev(amt = 5000, cmt = "HEPD", time = 24, ii = 8,
             addl = max(0, floor((dur - 24)/8) - 1)))
  } else if (inp$hep_type == "None (VITT / spontaneous)") {
    add(ev(amt = 1, cmt = "SURG"))
    add(ev(amt = inp$vaccine_ag, cmt = "PF4DNA"))
  }

  # ---- treatment, started on the day of recognition ----
  t0  <- inp$dx_day * 24
  tdur<- max(24, (inp$tx_days) * 24)
  switch(inp$anticoag,
    "None"         = NULL,
    "Argatroban"   = add(ev(amt = inp$arg_rate * wt * 60 * tdur,
                            rate = inp$arg_rate * wt * 60, cmt = "ARGC", time = t0)),
    "Bivalirudin"  = add(ev(amt = inp$biv_rate * wt * 1000 * tdur,
                            rate = inp$biv_rate * wt * 1000, cmt = "BIVC", time = t0)),
    "Fondaparinux" = add(ev(amt = 7500, cmt = "FOND", time = t0, ii = 24,
                            addl = max(0, inp$tx_days - 1))),
    "Rivaroxaban"  = add(ev(amt = 15, cmt = "RIVD", time = t0, ii = 12,
                            addl = max(0, inp$tx_days * 2 - 1))),
    "Danaparoid"   = add(ev(amt = 1500, cmt = "DANC", time = t0, ii = 8,
                            addl = max(0, inp$tx_days * 3 - 1)))
  )
  if (inp$ivig)   add(ev(amt = inp$ivig_gkg * wt, cmt = "IVGC", time = t0,
                         ii = 24, addl = inp$ivig_days - 1))
  if (inp$plex)   add(ev(amt = 1, cmt = "PLEXA", time = t0, ii = 24, addl = inp$plex_n - 1))
  if (inp$rtx)    add(ev(amt = 700, cmt = "RTXC", time = t0, ii = 168, addl = 3))
  if (inp$plt_txn) add(ev(amt = 40, cmt = "PLT", time = t0, ii = 24, addl = 3))
  if (inp$warfarin) {
    tw <- inp$war_day * 24
    add(ev(amt = inp$war_load, cmt = "WARD", time = tw))
    add(ev(amt = inp$war_maint, cmt = "WARD", time = tw + 24, ii = 24, addl = 29))
    if (inp$vitk) add(ev(amt = 10, cmt = "VITK", time = tw + 24))
  }
  if (!length(ev_list)) return(ev(amt = 0, cmt = "HEP"))
  Reduce(c, ev_list)
}

build_params <- function(inp) {
  list(WT = inp$wt, PLT0 = inp$plt0, HIND = inp$hind, FCG = inp$fcg,
       FPATH = inp$fpath, BCELL0 = inp$bcell0, RENFN = inp$renfn,
       HEPFN = inp$hepfn, KTREG = inp$ktreg, STER = inp$ster,
       KRECALL = inp$krecall)
}

run_model <- function(inp, end_days = NULL) {
  end <- (if (is.null(end_days)) inp$horizon else end_days) * 24
  o <- as.data.frame(mrgsim(param(mod, build_params(inp)), build_events(inp),
                            end = end, delta = 2))
  o$day <- o$time/24
  o
}

# ---------------------------------------------------------------------
# 4Ts, computed from the trajectory
# ---------------------------------------------------------------------
score_4ts <- function(o) {
  fall <- max(o$PLTFALL); nad <- min(o$PLT)
  t1 <- if (fall > 50 && nad >= 20) 2 else if (fall >= 30 || (nad >= 10 && nad < 20)) 1 else 0
  idx <- which(o$PLTFALL > 30)
  d   <- if (length(idx)) o$day[idx[1]] else NA_real_
  t2 <- if (is.na(d)) 0 else if (d >= 5 && d <= 10) 2 else if (d > 10) 1 else 0
  te <- max(o$TEC)
  t3 <- if (te > 0.30) 2 else if (te > 0.10) 1 else 0
  data.frame(
    Component = c("T1 Thrombocytopenia", "T2 Timing of fall",
                  "T3 Thrombosis", "T4 oTher causes absent", "TOTAL"),
    Score = c(t1, t2, t3, 2, t1 + t2 + t3 + 2),
    Basis = c(sprintf("max fall %.1f%%, nadir %.0f", fall, nad),
              if (is.na(d)) "no fall >30%" else sprintf("30%% fall on day %.1f", d),
              sprintf("modelled event probability %.1f%%", 100*te),
              "assumed (single-mechanism simulation)",
              c("low (0-3)", "intermediate (4-5)", "high (6-8)")[
                findInterval(t1 + t2 + t3 + 2, c(0, 4, 6))])
  )
}

lp <- function(o, vars, labs, title, sub, ylab, logy = FALSE) {
  d <- o[, c("day", vars)] |> pivot_longer(-day)
  d$name <- factor(d$name, levels = vars, labels = labs)
  p <- ggplot(d, aes(day, value, colour = name)) +
    geom_line(linewidth = 0.85) +
    scale_colour_manual(values = PAL[seq_along(vars)]) +
    labs(title = title, subtitle = sub, x = "Day", y = ylab) + theme_hit()
  if (logy) p <- p + scale_y_log10()
  p
}

# =====================================================================
# UI
# =====================================================================
ui <- fluidPage(
  titlePanel("Heparin-Induced Thrombocytopenia — QSP dashboard"),
  tags$p(style = "color:#555;margin-top:-10px",
         tags$b("HIT is an immune complex disease in which the drug is one half of the antigen."),
         " The antigen is a stoichiometric complex, so it is gated by chain length, chain",
         " concentration and the PF4:heparin molar ratio — not by dose alone."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("wt", "Weight (kg)", 40, 140, WT_DEF, 1),
      sliderInput("plt0", "Baseline platelets (10^9/L)", 100, 450, 250, 10),
      sliderInput("horizon", "Simulation horizon (days)", 20, 200, 60, 5),
      selectInput("fcg", "FcgammaRIIa 131 genotype",
                  c("HH (high affinity) = 1.35" = 1.35, "HR (reference) = 1.0" = 1.0,
                    "RR (low affinity) = 0.75" = 0.75), selected = 1.0),
      sliderInput("renfn", "Renal function (CrCl/100)", 0.2, 1.4, 1.0, 0.1),
      sliderInput("hepfn", "Hepatic function multiplier", 0.2, 1.2, 1.0, 0.1),

      hr(), h4("Index heparin exposure"),
      selectInput("hep_type", "Exposure",
                  c("UFH intravenous", "UFH subcutaneous", "LMWH (enoxaparin)",
                    "Fondaparinux", "Danaparoid", "Cardiopulmonary bypass",
                    "None (VITT / spontaneous)")),
      sliderInput("hep_days", "Duration (days)", 1, 40, 14, 1),
      conditionalPanel("input.hep_type == 'UFH intravenous'",
        sliderInput("ufh_bolus", "Bolus (U/kg)", 0, 120, 80, 5),
        sliderInput("ufh_rate", "Infusion (U/kg/h)", 0, 30, 18, 1)),
      conditionalPanel("input.hep_type == 'UFH subcutaneous'",
        sliderInput("ufh_sc_dose", "Dose (U)", 2500, 15000, 5000, 500),
        selectInput("ufh_sc_ii", "Interval (h)", c(8, 12), selected = 8)),
      conditionalPanel("input.hep_type == 'LMWH (enoxaparin)'",
        sliderInput("lmw_mg", "Dose (mg)", 20, 120, 70, 5),
        selectInput("lmw_ii", "Interval (h)", c(12, 24), selected = 12)),
      conditionalPanel("input.hep_type == 'Fondaparinux'",
        sliderInput("fon_mg", "Dose (mg)", 1.5, 10, 2.5, 0.5)),
      conditionalPanel("input.hep_type == 'Danaparoid'",
        sliderInput("dan_u", "Dose (anti-Xa U q8h)", 250, 2000, 750, 50)),
      conditionalPanel("input.hep_type == 'Cardiopulmonary bypass'",
        sliderInput("cpb_ukg", "Bypass heparin (U/kg)", 200, 600, 350, 25),
        checkboxInput("protamine", "Protamine reversal", TRUE),
        sliderInput("cpb_h", "Reversal at (h)", 1, 12, 4, 1)),
      conditionalPanel("input.hep_type == 'None (VITT / spontaneous)'",
        sliderInput("vaccine_ag", "Exogenous polyanion (nM)", 0, 200, 40, 10)),

      hr(), h4("Variant biology"),
      sliderInput("hind", "HIND — heparin independence", 0, 1, 0, 0.05),
      tags$small(style = "color:#777",
        "0 = classic HIT · 0.3-0.7 = autoimmune HIT · 0.8-1.0 = VITT. ",
        "One coefficient turns classic HIT into every variant."),
      sliderInput("fpath", "Platelet-activating IgG fraction", 0.01, 0.6, 0.15, 0.01),
      sliderInput("bcell0", "Pre-existing anti-PF4 clone size", 0.2, 20, 1, 0.2),
      sliderInput("krecall", "Memory recall (0 = no boost, as observed)", 0, 0.02, 0, 0.001),
      sliderInput("ktreg", "Regulatory suppression", 0, 0.8, 0, 0.05),

      hr(), h4("Recognition & treatment"),
      sliderInput("dx_day", "HIT recognised on day", 3, 30, 7, 1),
      tags$small(style = "color:#777", "Index heparin is not extended past its duration; this sets when treatment starts."),
      selectInput("anticoag", "Non-heparin anticoagulant",
                  c("None", "Argatroban", "Bivalirudin", "Fondaparinux",
                    "Rivaroxaban", "Danaparoid"), selected = "Argatroban"),
      sliderInput("tx_days", "Treatment duration (days)", 1, 40, 20, 1),
      conditionalPanel("input.anticoag == 'Argatroban'",
        sliderInput("arg_rate", "Argatroban (ug/kg/min)", 0.1, 4, 2, 0.1)),
      conditionalPanel("input.anticoag == 'Bivalirudin'",
        sliderInput("biv_rate", "Bivalirudin (mg/kg/h)", 0.02, 0.5, 0.15, 0.01)),
      checkboxInput("ivig", "High-dose IVIG", FALSE),
      conditionalPanel("input.ivig",
        sliderInput("ivig_gkg", "IVIG (g/kg/day)", 0.25, 2, 1, 0.25),
        sliderInput("ivig_days", "Days", 1, 5, 2, 1)),
      checkboxInput("plex", "Therapeutic plasma exchange", FALSE),
      conditionalPanel("input.plex", sliderInput("plex_n", "Sessions", 1, 6, 3, 1)),
      checkboxInput("rtx", "Rituximab 375 mg/m2 weekly x4", FALSE),
      checkboxInput("plt_txn", "Platelet transfusion (anti-therapy)", FALSE),
      checkboxInput("ster", "Corticosteroid adjunct", FALSE),

      hr(), h4("Warfarin transition"),
      checkboxInput("warfarin", "Start a vitamin K antagonist", FALSE),
      conditionalPanel("input.warfarin",
        sliderInput("war_day", "Start on day", 4, 40, 18, 1),
        sliderInput("war_load", "Loading dose (mg)", 0, 15, 10, 1),
        sliderInput("war_maint", "Maintenance (mg/day)", 1, 10, 5, 0.5),
        checkboxInput("vitk", "Vitamin K rescue at 24 h", FALSE))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · Patient & exposure", br(),
                 fluidRow(column(7, verbatimTextOutput("summary")),
                          column(5, h5("4Ts score, computed from the trajectory"),
                                 tableOutput("fourts"))),
                 plotOutput("overview", height = 430)),
        tabPanel("2 · Antigen (bell curve)", br(),
                 tags$p(tags$b("The central non-linearity."), " The antigen is the ultralarge",
                        " complex, and its abundance is bell-shaped in heparin: too little and",
                        " there is nothing to bridge with, too much and every PF4 tetramer gets",
                        " its own chain."),
                 plotOutput("antigen", height = 330),
                 plotOutput("bellcurve", height = 330)),
        tabPanel("3 · Heparin & drug PK", br(),
                 plotOutput("pk_hep", height = 300), plotOutput("pk_drug", height = 300),
                 plotOutput("chain", height = 260)),
        tabPanel("4 · Immunity", br(),
                 tags$p("Plasmablasts dominate and long-lived plasma cells are near zero:",
                        " that is why the antibody is transient and why re-exposure gives no",
                        " anamnestic boost."),
                 plotOutput("immune", height = 320), plotOutput("igg", height = 320)),
        tabPanel("5 · Platelet axis", br(),
                 plotOutput("platelets", height = 340), plotOutput("megak", height = 280),
                 verbatimTextOutput("plt_txt")),
        tabPanel("6 · Thrombin axis", br(),
                 tags$p("The same signal that clears platelets also builds thrombin.",
                        " A direct thrombin inhibitor cuts this panel and leaves tab 5 alone."),
                 plotOutput("coag", height = 320), plotOutput("markers", height = 300)),
        tabPanel("7 · Clinical endpoints", br(),
                 plotOutput("endpoints", height = 340), plotOutput("hazard", height = 280),
                 verbatimTextOutput("ep_txt")),
        tabPanel("8 · Diagnostics", br(),
                 tags$p(tags$b("The high-heparin step is not a rule, it is the bell curve."),
                        " In classic HIT, release at 100 U/mL is abolished. Raise HIND and it",
                        " is not — which is why autoimmune HIT and VITT sera need",
                        " PF4-enhanced rather than heparin-dependent assays."),
                 plotOutput("assays", height = 340), plotOutput("iceberg", height = 300),
                 tableOutput("dx_tab")),
        tabPanel("9 · The two traps", br(),
                 tags$p(tags$b("Trap 1 — the half-life mismatch."),
                        " Protein C t1/2 8 h, prothrombin t1/2 60-72 h."),
                 plotOutput("trap_war", height = 320),
                 tags$p(tags$b("Trap 2 — the argatroban INR artefact."),
                        " Measured INR, true INR, and chromogenic factor X on one axis."),
                 plotOutput("trap_inr", height = 320)),
        tabPanel("10 · Scenario comparison", br(),
                 tags$p("Every arm uses the current patient and index exposure; only the",
                        " treatment differs."),
                 plotOutput("cmp_plt", height = 320), plotOutput("cmp_tec", height = 320),
                 DTOutput("cmp_tab"))
      )
    )
  )
)

# =====================================================================
# SERVER
# =====================================================================
server <- function(input, output, session) {

  o <- reactive(run_model(input))

  output$summary <- renderPrint({
    d <- o(); s4 <- score_4ts(d)
    cat(sprintf("Exposure                 : %s for %d days\n", input$hep_type, input$hep_days))
    cat(sprintf("Peak local PF4           : %8.0f nM   (plasma %.2f nM)\n",
                max(d$PF4Lo), max(d$PF4P)))
    cat(sprintf("Peak polyanion chains    : %8.0f nM\n", max(d$HEPTo)))
    cat(sprintf("PF4:chain ratio at peak  : %8.2f     (optimum >= 8)\n",
                d$RRATo[which.max(d$ULC)]))
    cat(sprintf("Peak ULC (the antigen)   : %8.1f nM\n", max(d$ULC)))
    cat(sprintf("Peak anti-PF4/H IgG      : %8.1f U/mL  (ELISA OD %.2f)\n",
                max(d$IGGP + d$IGGN), max(d$OD)))
    cat(sprintf("Platelet nadir           : %8.1f on day %.1f (fall %.1f%%)\n",
                min(d$PLT), d$day[which.min(d$PLT)], max(d$PLTFALL)))
    cat(sprintf("Peak thrombin activity   : %8.2f nM   (baseline 0.5)\n", max(d$THRACTo)))
    cat(sprintf("Thrombotic event risk    : %8.1f%% by day %d\n",
                100*tail(d$TEC, 1), input$horizon))
    cat(sprintf("Venous limb gangrene risk: %8.2f%%\n", 100*tail(d$NECR, 1)))
    cat(sprintf("SRA: %.0f%% at low heparin | %.0f%% at 100 U/mL | %.0f%% with none\n",
                max(d$SRALO), max(d$SRAHI), max(d$SRA0)))
    cat(sprintf("Assay pattern            : %s\n",
      if (max(d$SRAPOSA) > 0) "HEPARIN-INDEPENDENT (aHIT / VITT pattern)"
      else if (max(d$SRAPOS) > 0) "classic HIT (positive low, inhibited at 100 U/mL)"
      else if (max(d$SEROPOS) > 0) "seropositive, not platelet-activating"
      else "no antibody detected"))
    cat(sprintf("4Ts total                : %d\n", s4$Score[5]))
  })

  output$fourts  <- renderTable(score_4ts(o()), rownames = FALSE)

  output$overview <- renderPlot({
    d <- o()
    df <- data.frame(day = d$day,
                     `Platelets (10^9/L)` = d$PLT,
                     `ULC antigen (nM)` = d$ULC,
                     `anti-PF4 IgG (U/mL)` = d$IGGP + d$IGGN,
                     `Thrombin (nM)` = d$THRACTo,
                     `Thrombosis risk (%)` = 100*d$TEC,
                     check.names = FALSE) |> pivot_longer(-day)
    ggplot(df, aes(day, value)) + geom_line(colour = "#C62828", linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(title = "The disease in five quantities",
           subtitle = "antigen -> antibody -> platelet loss and thrombin, in series",
           x = "Day", y = NULL) + theme_hit()
  })

  # ---- 2. antigen -------------------------------------------------
  output$antigen <- renderPlot({
    d <- o()
    df <- data.frame(day = d$day,
                     `Local PF4 (nM)` = d$PF4Lo,
                     `Polyanion chains (nM)` = d$HEPTo,
                     `PF4:chain ratio` = d$RRATo,
                     `Antigenicity of ratio (0-1)` = d$FRATo,
                     `ULC (nM)` = d$ULC,
                     `Immune complex (nM)` = d$IC,
                     check.names = FALSE) |> pivot_longer(-day)
    ggplot(df, aes(day, value)) + geom_line(colour = "#1565C0", linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(title = "Assembling the antigen",
           subtitle = "the ratio panel is the gate: antigenicity collapses in heparin excess",
           x = "Day", y = NULL) + theme_hit()
  })

  output$bellcurve <- renderPlot({
    hs <- 10^seq(-3, 2.4, length.out = 40)
    pr <- build_params(input)
    d0 <- o()
    igg <- max(d0$IGGP)
    res <- sapply(hs, function(h) {
      s <- as.data.frame(mrgsim(param(mod, c(pr, list(HEPSRALO = h))),
                                build_events(input), end = input$horizon*24, delta = 24))
      max(s$SRALO)
    })
    df <- data.frame(heparin = hs, release = res)
    ggplot(df, aes(heparin, release)) +
      annotate("rect", xmin = 0.03, xmax = 2, ymin = -Inf, ymax = Inf,
               fill = "#2E7D32", alpha = 0.08) +
      annotate("rect", xmin = 30, xmax = Inf, ymin = -Inf, ymax = Inf,
               fill = "#C62828", alpha = 0.08) +
      geom_line(linewidth = 1.1, colour = "#C62828") +
      geom_vline(xintercept = 100, linetype = 2, colour = "grey40") +
      annotate("text", x = 100, y = 95, label = "SRA confirmatory\n100 U/mL",
               hjust = 1.1, size = 3, colour = "grey30") +
      annotate("text", x = 0.25, y = 15, label = "optimal ratio", size = 3.2, colour = "#2E7D32") +
      annotate("text", x = 90, y = 60, label = "heparin excess", size = 3.2, colour = "#C62828") +
      scale_x_log10() +
      labs(title = "The bell curve, swept over four decades of heparin",
           subtitle = paste0("in-vitro platelet activation by this patient's serum ",
                             "(platelet-activating IgG ", round(igg, 1), " U/mL)"),
           x = "Heparin in the assay (U/mL, log scale)",
           y = "Serotonin release (%)") + theme_hit()
  })

  # ---- 3. PK ------------------------------------------------------
  output$pk_hep <- renderPlot(
    lp(o(), c("CUFHo", "CLMWo", "CFONo", "CDANo"),
       c("UFH (U/mL)", "LMWH (anti-Xa IU/mL)", "Fondaparinux (ug/mL)", "Danaparoid (U/mL)"),
       "Heparin-family pharmacokinetics",
       "UFH clearance is saturable, so its half-life lengthens with dose", "Concentration"))
  output$pk_drug <- renderPlot(
    lp(o(), c("CARGo", "CBIVo", "CRIVo", "CWARo", "CRTXo"),
       c("Argatroban (ug/mL)", "Bivalirudin (ug/mL)", "Rivaroxaban (ng/mL)",
         "Warfarin (mg/L)", "Rituximab (ug/mL)"),
       "Non-heparin anticoagulants and immunomodulators",
       "argatroban is hepatically cleared; fondaparinux, danaparoid and DOACs are renal",
       "Concentration"))
  output$chain <- renderPlot({
    df <- data.frame(
      species = factor(c("UFH", "LMWH", "danaparoid", "fondaparinux", "DNA / polyP"),
                       levels = c("fondaparinux", "danaparoid", "LMWH", "UFH", "DNA / polyP")),
      saccharides = c(50, 15, 9, 5, 200))
    df$bridging <- df$saccharides^6/(df$saccharides^6 + 24^6)
    ggplot(df, aes(species, bridging, fill = species)) +
      geom_col(width = 0.6) + coord_flip() +
      scale_fill_manual(values = PAL, guide = "none") +
      geom_text(aes(label = sprintf("%.4f", bridging)), hjust = -0.1, size = 3.4) +
      ylim(0, 1.2) +
      labs(title = "Bridging competence is the class effect",
           subtitle = "a chain must reach two PF4 tetramers; fondaparinux is a structural null",
           x = NULL, y = "Fraction of chains able to bridge") + theme_hit()
  })

  # ---- 4. immunity ------------------------------------------------
  output$immune <- renderPlot(
    lp(o(), c("BCELL", "PBLAST", "LLPC", "MEMB"),
       c("B-cell clone", "Plasmablasts", "Long-lived plasma cells", "Memory B cells"),
       "The cellular response is a plasmablast burst, not a memory response",
       "long-lived plasma cells stay near zero: that is why the antibody disappears",
       "Cells (arbitrary units)"))
  output$igg <- renderPlot({
    d <- o()
    df <- data.frame(day = d$day,
                     `Platelet-activating IgG (U/mL)` = d$IGGP,
                     `Non-activating IgG (U/mL)` = d$IGGN,
                     `ELISA optical density` = d$OD,
                     check.names = FALSE) |> pivot_longer(-day)
    ggplot(df, aes(day, value)) + geom_line(colour = "#6A1B9A", linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      geom_hline(data = data.frame(name = "ELISA optical density", y = 0.4),
                 aes(yintercept = y), linetype = 2, colour = "grey40") +
      labs(title = "Antibody, and the assay that sees it",
           subtitle = "most anti-PF4 IgG is not platelet-activating: the mass of the iceberg",
           x = "Day", y = NULL) + theme_hit()
  })

  # ---- 5. platelets -----------------------------------------------
  output$platelets <- renderPlot({
    d <- o()
    ggplot(d, aes(day, PLT)) +
      geom_hline(yintercept = 150, linetype = 2, colour = "grey50") +
      geom_hline(yintercept = 0.5*input$plt0, linetype = 3, colour = "#C62828") +
      geom_line(linewidth = 1, colour = "#4527A0") +
      annotate("text", x = max(d$day), y = 150, hjust = 1, vjust = -0.5, size = 3,
               colour = "grey40", label = "150 x10^9/L") +
      annotate("text", x = max(d$day), y = 0.5*input$plt0, hjust = 1, vjust = 1.5,
               size = 3, colour = "#C62828", label = "50% fall — the diagnostic quantity") +
      labs(title = "The platelet axis",
           subtitle = "a >50% RELATIVE fall matters more than any absolute count; nadir is rarely <20",
           x = "Day", y = "Platelets (10^9/L)") + theme_hit()
  })
  output$megak <- renderPlot(
    lp(o(), c("PLTA", "TPO", "MK2", "MP"),
       c("Activated platelets (10^9/L)", "Thrombopoietin (x baseline)",
         "Megakaryocyte output", "Procoagulant microparticles"),
       "Consumption, compensation and the procoagulant debris",
       "thrombopoietin rises as platelets fall, but marrow output lags 3-4 days", "Level"))
  output$plt_txt <- renderPrint({
    d <- o()
    rec <- which(d$day > d$day[which.min(d$PLT)] & d$PLT > 150)
    cat(sprintf("Nadir %.1f x10^9/L on day %.1f (%.1f%% fall from %.0f)\n",
                min(d$PLT), d$day[which.min(d$PLT)], max(d$PLTFALL), input$plt0))
    cat(sprintf("Recovery above 150 x10^9/L: %s\n",
                if (length(rec)) sprintf("day %.1f", d$day[rec[1]]) else "not within horizon"))
    cat(sprintf("Peak FcgammaRIIa occupancy %.3f, crosslinking signal %.3f\n",
                max(d$OCCo), max(d$XLo)))
    if (input$plt_txn) cat("\nNOTE: platelet transfusion is switched on. Watch tab 6 and 7:\n",
                           "the count improves while the thrombin axis does not.\n")
  })

  # ---- 6. thrombin ------------------------------------------------
  output$coag <- renderPlot(
    lp(o(), c("THRACTo", "APCo", "FFVo", "APTT"),
       c("Effective thrombin (nM)", "Activated protein C", "Factor Va availability", "aPTT (s)"),
       "The thrombin axis", "HIT generates among the highest thrombin of any acquired thrombophilia",
       "Level"))
  output$markers <- renderPlot({
    d <- o()
    df <- data.frame(day = d$day, `D-dimer (ug/mL FEU)` = d$DDIM,
                     `Prothrombin fragment 1.2 (nmol/L)` = d$FRG,
                     `Fibrin burden` = d$FBN, `Antithrombin (%)` = d$AT,
                     `Tissue factor` = d$TF, `von Willebrand factor (%)` = d$VWF,
                     check.names = FALSE) |> pivot_longer(-day)
    ggplot(df, aes(day, value)) + geom_line(colour = "#880E4F", linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(title = "Coagulation biomarkers", x = "Day", y = NULL) + theme_hit()
  })

  # ---- 7. endpoints ----------------------------------------------
  output$endpoints <- renderPlot({
    d <- o()
    df <- data.frame(day = d$day,
                     `Thrombotic event probability (%)` = 100*d$TEC,
                     `Venous limb gangrene probability (%)` = 100*d$NECR,
                     check.names = FALSE) |> pivot_longer(-day)
    ggplot(df, aes(day, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = c("#BF360C", "#E65100")) +
      labs(title = "Clinical endpoints",
           subtitle = "cumulative, HIT-attributable (excess over the background rate)",
           x = "Day", y = "Probability (%)") + theme_hit()
  })
  output$hazard <- renderPlot(
    lp(o(), c("THRACTo", "XLo"),
       c("Thrombin route (nM)", "Platelet-activation route (0-1)"),
       "Two routes to thrombosis, cut by different drugs",
       "anticoagulants act on the first only; heparin cessation, IVIG and PLEX act on the second",
       "Level"))
  output$ep_txt <- renderPrint({
    d <- o()
    cat(sprintf("Thrombotic event probability : %.1f%% at day %d\n",
                100*tail(d$TEC, 1), input$horizon))
    cat(sprintf("Venous limb gangrene         : %.2f%%\n", 100*tail(d$NECR, 1)))
    cat(sprintf("Cumulative thrombin exposure : %.0f nM.h\n", tail(d$THRAUC, 1)))
    if (input$warfarin && min(d$PC) < 60 && max(d$THRACTo) > 5)
      cat("\nWARNING (model): a vitamin K antagonist was started while thrombin was still\n",
          "high. Protein C falls far faster than prothrombin, so the net effect is a\n",
          "transient INCREASE in thrombin. See tab 9.\n")
  })

  # ---- 8. diagnostics --------------------------------------------
  output$assays <- renderPlot(
    lp(o(), c("SRALO", "SRAHI", "SRA0"),
       c("SRA at 0.2 U/mL heparin", "SRA at 100 U/mL heparin", "SRA with NO heparin"),
       "The functional assay, simulated",
       "classic HIT: high at low heparin, abolished at 100 U/mL. Raise HIND and that inhibition is lost.",
       "Serotonin release (%)") +
      geom_hline(yintercept = 20, linetype = 2, colour = "grey40"))
  output$iceberg <- renderPlot({
    rates <- c(0.5, 1, 2, 4, 8, 12, 18, 25)
    pr <- build_params(input); wt <- input$wt
    res <- do.call(rbind, lapply(rates, function(r) {
      e <- c(ev(amt = 4*r*wt, cmt = "HEP"), ev(amt = r*wt*336, rate = r*wt, cmt = "HEP"))
      s <- as.data.frame(mrgsim(param(mod, pr), e, end = 1440, delta = 8))
      data.frame(rate = r,
                 `1 seropositive` = 100*max(s$SEROPOS),
                 `2 platelet-activating` = 100*max(s$SRAPOS),
                 `3 clinical HIT` = 100*as.integer(max(s$PLTFALL) > 50),
                 check.names = FALSE)
    })) |> pivot_longer(-rate)
    ggplot(res, aes(rate, value, colour = name)) +
      geom_step(linewidth = 1.1) + scale_x_log10() +
      scale_colour_manual(values = c("#1565C0", "#E65100", "#C62828")) +
      labs(title = "The iceberg is generated by exposure intensity",
           subtitle = "three successive waterlines as UFH infusion rate rises (this patient's genotype)",
           x = "UFH infusion rate (U/kg/h, log scale)", y = "Level reached (%)") + theme_hit()
  })
  output$dx_tab <- renderTable({
    d <- o()
    data.frame(Assay = c("PF4/heparin IgG ELISA (OD)", "Rapid immunoassay",
                         "SRA, 0.2 U/mL heparin", "SRA, 100 U/mL heparin",
                         "SRA, no heparin (PF4-enhanced)", "Interpretation"),
               Result = c(sprintf("%.2f", max(d$OD)),
                          if (max(d$SEROPOS) > 0) "positive" else "negative",
                          sprintf("%.0f%%", max(d$SRALO)),
                          sprintf("%.0f%%", max(d$SRAHI)),
                          sprintf("%.0f%%", max(d$SRA0)),
                          if (max(d$SRAPOSA) > 0) "heparin-INDEPENDENT activation: aHIT / VITT"
                          else if (max(d$SRAPOS) > 0) "classic heparin-dependent HIT"
                          else if (max(d$SEROPOS) > 0) "seropositive but not platelet-activating"
                          else "no anti-PF4 antibody"),
               Threshold = c(">0.40 positive", "negative rules out", ">20% positive",
                             "<20% confirms classic HIT", "raised only in aHIT / VITT", ""))
  }, rownames = FALSE)

  # ---- 9. traps ---------------------------------------------------
  output$trap_war <- renderPlot({
    pr <- build_params(input); wt <- input$wt
    base <- c(ev(amt = 80*wt, cmt = "HEP"),
              ev(amt = 18*wt*input$dx_day*24, rate = 18*wt, cmt = "HEP"))
    arms <- lapply(c(7, 11, 15, 21), function(sd) {
      e <- c(base, ev(amt = 10, cmt = "WARD", time = sd*24),
             ev(amt = 5, cmt = "WARD", time = sd*24 + 24, ii = 24, addl = 25))
      s <- as.data.frame(mrgsim(param(mod, pr), e, end = 1440, delta = 2))
      s$day <- s$time/24; s$arm <- paste0("warfarin from day ", sd); s
    })
    d <- bind_rows(arms)
    df <- data.frame(day = d$day, arm = d$arm, `Protein C (%)` = d$PC,
                     `Prothrombin (%)` = d$FII, `Thrombin (nM)` = d$THRACTo,
                     `Venous limb gangrene (%)` = 100*d$NECR, check.names = FALSE) |>
      pivot_longer(c(-day, -arm))
    ggplot(df, aes(day, value, colour = arm)) + geom_line(linewidth = 0.85) +
      facet_wrap(~name, scales = "free_y") + xlim(0, 45) +
      scale_colour_manual(values = PAL) +
      labs(title = "The half-life mismatch, and the rule it generates",
           subtitle = "protein C (t1/2 8 h) collapses long before prothrombin (t1/2 60 h)",
           x = "Day", y = NULL) + theme_hit()
  })
  output$trap_inr <- renderPlot({
    pr <- build_params(input); wt <- input$wt
    e <- c(ev(amt = 80*wt, cmt = "HEP"),
           ev(amt = 18*wt*input$dx_day*24, rate = 18*wt, cmt = "HEP"),
           ev(amt = 2*wt*60*360, rate = 2*wt*60, cmt = "ARGC", time = input$dx_day*24),
           ev(amt = 10, cmt = "WARD", time = 18*24),
           ev(amt = 5, cmt = "WARD", time = 19*24, ii = 24, addl = 25))
    s <- as.data.frame(mrgsim(param(mod, pr), e, end = 1200, delta = 2)); s$day <- s$time/24
    df <- data.frame(day = s$day, `Measured INR (with DTI artefact)` = s$INR,
                     `True INR (factor depletion only)` = s$INRTRUE,
                     `Chromogenic factor X (%)/25` = s$FXCHROM/25,
                     check.names = FALSE) |> pivot_longer(-day)
    ggplot(df, aes(day, value, colour = name)) + geom_line(linewidth = 0.95) +
      geom_hline(yintercept = c(2, 3), linetype = 2, colour = "grey55") +
      xlim(15, 40) + scale_colour_manual(values = c("#C62828", "#1565C0", "#2E7D32")) +
      labs(title = "An INR of 4 on argatroban does not mean the warfarin is therapeutic",
           subtitle = "factor X (green, scaled /25) is the correct surrogate during overlap",
           x = "Day", y = "INR  /  scaled factor X") + theme_hit()
  })

  # ---- 10. comparison --------------------------------------------
  cmp <- reactive({
    pr <- build_params(input); wt <- input$wt; t0 <- input$dx_day*24
    dur <- input$tx_days*24
    hep <- build_events(modifyList(reactiveValuesToList(input),
                                   list(anticoag = "None", ivig = FALSE, plex = FALSE,
                                        rtx = FALSE, plt_txn = FALSE, warfarin = FALSE)))
    arms <- list(
      "stop heparin only" = NULL,
      "argatroban"        = ev(amt = 2*wt*60*dur, rate = 2*wt*60, cmt = "ARGC", time = t0),
      "bivalirudin"       = ev(amt = 0.15*wt*1000*dur, rate = 0.15*wt*1000, cmt = "BIVC", time = t0),
      "fondaparinux"      = ev(amt = 7500, cmt = "FOND", time = t0, ii = 24, addl = input$tx_days - 1),
      "rivaroxaban"       = ev(amt = 15, cmt = "RIVD", time = t0, ii = 12, addl = input$tx_days*2 - 1),
      "argatroban + IVIG" = c(ev(amt = 2*wt*60*dur, rate = 2*wt*60, cmt = "ARGC", time = t0),
                              ev(amt = wt, cmt = "IVGC", time = t0, ii = 24, addl = 1)),
      "platelet transfusion" = ev(amt = 40, cmt = "PLT", time = t0, ii = 24, addl = 3))
    bind_rows(lapply(names(arms), function(n) {
      e <- if (is.null(arms[[n]])) hep else c(hep, arms[[n]])
      s <- as.data.frame(mrgsim(param(mod, pr), e, end = input$horizon*24, delta = 4))
      s$day <- s$time/24; s$arm <- n; s
    }))
  })
  output$cmp_plt <- renderPlot({
    ggplot(cmp(), aes(day, PLT, colour = arm)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 150, linetype = 2, colour = "grey50") +
      scale_colour_manual(values = PAL) +
      labs(title = "Platelet axis by treatment arm",
           subtitle = "anticoagulants barely move this panel: platelets recover on the immune system's schedule",
           x = "Day", y = "Platelets (10^9/L)") + theme_hit()
  })
  output$cmp_tec <- renderPlot({
    ggplot(cmp(), aes(day, 100*TEC, colour = arm)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(title = "Thrombotic risk by treatment arm",
           subtitle = "and here they separate sharply — the second of the two axes",
           x = "Day", y = "Cumulative thrombotic event probability (%)") + theme_hit()
  })
  output$cmp_tab <- renderDT({
    cmp() |> group_by(arm) |>
      summarise(`platelet nadir` = round(min(PLT), 1),
                `nadir day` = round(day[which.min(PLT)], 1),
                `peak thrombin (nM)` = round(max(THRACTo), 2),
                `peak crosslink signal` = round(max(XLo), 3),
                `thrombosis risk (%)` = round(100*max(TEC), 1),
                `VLG risk (%)` = round(100*max(NECR), 2), .groups = "drop") |>
      arrange(`thrombosis risk (%)`) |>
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 10))
  })
}

shinyApp(ui, server)
