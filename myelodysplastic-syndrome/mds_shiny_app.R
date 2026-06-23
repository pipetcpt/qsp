# =============================================================================
# MDS QSP Shiny Dashboard
# Myelodysplastic Syndrome — Quantitative Systems Pharmacology Model
# =============================================================================

library(shiny)
library(dplyr)
library(ggplot2)
library(plotly)
library(bslib)
library(DT)
library(tidyr)

# =============================================================================
# Helper functions
# =============================================================================

# BSA (Mosteller formula)
calc_bsa <- function(weight_kg, height_cm = 170) {
  sqrt((weight_kg * height_cm) / 3600)
}

# IPSS-R scoring
calc_ipssr <- function(blast_pct, cyto_cat, hgb, plt, anc) {
  # Blast score
  blast_score <- dplyr::case_when(
    blast_pct <= 2  ~ 0,
    blast_pct <= 4  ~ 1,
    blast_pct <= 10 ~ 2,
    TRUE            ~ 3
  )
  # Cytogenetics score
  cyto_score <- dplyr::case_when(
    cyto_cat == "Very Good"   ~ 0,
    cyto_cat == "Good"        ~ 1,
    cyto_cat == "Intermediate"~ 2,
    cyto_cat == "Poor"        ~ 3,
    cyto_cat == "Very Poor"   ~ 4,
    TRUE                      ~ 2
  )
  # Hemoglobin score
  hgb_score <- dplyr::case_when(
    hgb >= 10  ~ 0,
    hgb >= 8   ~ 1,
    TRUE        ~ 1.5
  )
  # Platelet score
  plt_score <- dplyr::case_when(
    plt >= 100 ~ 0,
    plt >= 50  ~ 0.5,
    TRUE        ~ 1
  )
  # ANC score
  anc_score <- ifelse(anc >= 0.8, 0, 0.5)

  total <- blast_score + cyto_score + hgb_score + plt_score + anc_score

  risk_group <- dplyr::case_when(
    total <= 1.5  ~ "Very Low",
    total <= 3.0  ~ "Low",
    total <= 4.5  ~ "Intermediate",
    total <= 6.0  ~ "High",
    TRUE          ~ "Very High"
  )
  list(score = round(total, 1), risk_group = risk_group)
}

# IPSS-M estimate (simplified based on mutation burden)
calc_ipssm <- function(mutations, ipssr_score) {
  # Each driver mutation adds risk
  mut_count <- sum(mutations)
  tp53 <- "TP53" %in% mutations
  runx1 <- "RUNX1" %in% mutations
  additional <- mut_count * 0.4 + (if (tp53) 1.5 else 0) + (if (runx1) 0.8 else 0)
  score <- ipssr_score * 0.3 + additional
  risk_group <- dplyr::case_when(
    score < 0.5  ~ "Very Low",
    score < 1.0  ~ "Low",
    score < 2.0  ~ "Moderate Low",
    score < 3.0  ~ "Moderate High",
    score < 4.0  ~ "High",
    TRUE         ~ "Very High"
  )
  list(score = round(score, 2), risk_group = risk_group)
}

# =============================================================================
# PK Parameter Tables
# =============================================================================
pk_params <- list(
  AZA = list(
    dose_unit   = "75 mg/m² SC",
    route       = "SC",
    t_half      = 0.6,      # hours
    Vd          = 76,       # L
    CL          = 89,       # L/h
    F           = 0.89,
    tmax        = 0.5,
    dose_mg_m2  = 75,
    duration    = 7         # days on
  ),
  DEC = list(
    dose_unit   = "20 mg/m² IV",
    route       = "IV",
    t_half      = 0.5,
    Vd          = 63,
    CL          = 125,
    F           = 1.0,
    tmax        = 0.08,
    dose_mg_m2  = 20,
    duration    = 5
  ),
  `Oral-DEC` = list(
    dose_unit   = "35 mg oral",
    route       = "oral",
    t_half      = 1.3,
    Vd          = 50,
    CL          = 40,
    F           = 0.93,
    tmax        = 1.0,
    dose_mg_m2  = 35 / 1.8,  # fixed dose
    duration    = 5
  ),
  Lenalidomide = list(
    dose_unit   = "10 mg QD",
    route       = "oral",
    t_half      = 3,
    Vd          = 67,
    CL          = 15,
    F           = 0.90,
    tmax        = 1.0,
    dose_mg_m2  = 10 / 1.8,
    duration    = 21
  ),
  Luspatercept = list(
    dose_unit   = "1.0 mg/kg SC",
    route       = "SC",
    t_half      = 168,     # ~7 days (protein)
    Vd          = 5.5,
    CL          = 0.022,
    F           = 0.67,
    tmax        = 72,
    dose_mg_m2  = NA,
    duration    = 1        # q21d single dose
  ),
  Darbepoetin = list(
    dose_unit   = "500 mcg SC q3w",
    route       = "SC",
    t_half      = 70,
    Vd          = 3.4,
    CL          = 0.034,
    F           = 0.37,
    tmax        = 48,
    dose_mg_m2  = NA,
    duration    = 1
  ),
  Venetoclax = list(
    dose_unit   = "400 mg QD",
    route       = "oral",
    t_half      = 17,
    Vd          = 256,
    CL          = 10.4,
    F           = 0.72,
    tmax        = 5,
    dose_mg_m2  = 400 / 1.8,
    duration    = 28
  )
)

# =============================================================================
# PK Simulation
# =============================================================================
simulate_pk <- function(drug, n_cycles, bsa, weight_kg) {
  p  <- pk_params[[drug]]
  ke <- log(2) / p$t_half

  # Dose in mg
  if (!is.na(p$dose_mg_m2)) {
    dose_mg <- p$dose_mg_m2 * bsa
  } else {
    dose_mg <- if (drug == "Luspatercept") 1.0 * weight_kg
               else if (drug == "Darbepoetin") 0.5  # 500 mcg
               else 400
  }

  cycle_len <- 28  # days per cycle
  dose_days <- p$duration
  total_days <- n_cycles * cycle_len

  dt   <- 0.1   # h
  tmax_h <- total_days * 24
  times <- seq(0, tmax_h, by = dt)
  conc <- numeric(length(times))

  # Dose events
  dose_times_h <- c()
  for (cyc in 0:(n_cycles - 1)) {
    for (d in 0:(dose_days - 1)) {
      dose_times_h <- c(dose_times_h, (cyc * cycle_len + d) * 24)
    }
  }

  A <- 0  # drug amount in central compartment (mg)

  for (j in 2:length(times)) {
    t_now <- times[j - 1]
    t_new <- times[j]
    dt_h  <- t_new - t_now

    # Check for dose
    dose_given <- sum(dose_times_h >= t_now & dose_times_h < t_new)
    if (dose_given > 0) {
      if (p$route == "IV") {
        A <- A + dose_mg * dose_given
      } else {
        # Approximate SC/oral absorption as bolus at peak
        A <- A + dose_mg * p$F * dose_given
      }
    }

    A <- A * exp(-ke * dt_h)
    conc[j] <- A / p$Vd
  }

  # Summarise hourly
  idx <- seq(1, length(times), by = 10)
  df <- data.frame(
    time_h    = times[idx],
    time_day  = times[idx] / 24,
    conc_ugml = conc[idx]
  )

  # PK metrics per cycle 1
  cyc1_idx <- df$time_day <= cycle_len
  df_c1    <- df[cyc1_idx, ]
  cmax  <- max(df_c1$conc_ugml)
  tmax_obs <- df_c1$time_h[which.max(df_c1$conc_ugml)]
  auc   <- sum(df_c1$conc_ugml) * (cycle_len * 24 / nrow(df_c1))
  css   <- if (n_cycles > 1) {
    last_cyc <- df$time_day > (n_cycles - 1) * cycle_len
    mean(df$conc_ugml[last_cyc])
  } else NA_real_

  list(
    data     = df,
    metrics  = data.frame(
      Parameter = c("Cmax (μg/mL)", "Tmax (h)", "AUC0-28d (μg·h/mL)", "Css (μg/mL)"),
      Value     = round(c(cmax, tmax_obs, auc, css), 3)
    )
  )
}

# =============================================================================
# Disease Dynamics Simulation
# =============================================================================
simulate_disease <- function(treatment, days = 365,
                              blast_init = 10, vaf_init = 30,
                              hgb_init = 8.5, plt_init = 80, anc_init = 1.0,
                              methylation_init = 0.65) {
  dt   <- 1
  t    <- seq(0, days, by = dt)
  n    <- length(t)

  # Treatment parameters
  params <- list(
    BSC          = list(eff_blast = 0,    eff_met = 0,    eff_vaf = 0,    eff_hgb = 0,    eff_plt = 0,    eff_anc = 0,    cycle = 28, dur = 0),
    AZA          = list(eff_blast = 0.25, eff_met = 0.30, eff_vaf = 0.15, eff_hgb = 0.12, eff_plt = 0.10, eff_anc = 0.08, cycle = 28, dur = 7),
    DEC          = list(eff_blast = 0.30, eff_met = 0.35, eff_vaf = 0.18, eff_hgb = 0.10, eff_plt = 0.08, eff_anc = 0.07, cycle = 28, dur = 5),
    `Oral-DEC`   = list(eff_blast = 0.28, eff_met = 0.33, eff_vaf = 0.17, eff_hgb = 0.10, eff_plt = 0.08, eff_anc = 0.07, cycle = 28, dur = 5),
    Lenalidomide = list(eff_blast = 0.20, eff_met = 0.10, eff_vaf = 0.25, eff_hgb = 0.25, eff_plt = 0.15, eff_anc = 0.10, cycle = 28, dur = 21),
    Luspatercept = list(eff_blast = 0.05, eff_met = 0.05, eff_vaf = 0.03, eff_hgb = 0.30, eff_plt = 0.05, eff_anc = 0.02, cycle = 21, dur = 1),
    `VEN+AZA`    = list(eff_blast = 0.45, eff_met = 0.35, eff_vaf = 0.30, eff_hgb = 0.08, eff_plt = 0.05, eff_anc = 0.05, cycle = 28, dur = 7)
  )

  p <- if (treatment %in% names(params)) params[[treatment]] else params[["BSC"]]

  # Natural history rates (per day)
  blast_growth  <- 0.002
  vaf_drift     <- 0.0005
  hgb_decay     <- 0.005
  plt_decay     <- 0.003
  anc_decay     <- 0.002
  met_drift     <- 0.0003

  # Initialise
  Blast <- numeric(n);  Blast[1]  <- blast_init
  VAF   <- numeric(n);  VAF[1]    <- vaf_init
  Hgb   <- numeric(n);  Hgb[1]   <- hgb_init
  Plt   <- numeric(n);  Plt[1]    <- plt_init
  ANC   <- numeric(n);  ANC[1]   <- anc_init
  Met   <- numeric(n);  Met[1]   <- methylation_init

  for (i in 2:n) {
    day <- t[i]

    # Is drug active today?
    cycle_day <- day %% p$cycle
    drug_on   <- (p$dur > 0) && (cycle_day < p$dur)
    eff       <- if (drug_on) 1 else 0

    # Blast dynamics
    net_blast <- blast_growth * Blast[i-1] * (1 - Blast[i-1] / 30) -
                 p$eff_blast * eff * Blast[i-1]
    Blast[i]  <- max(0.5, min(30, Blast[i-1] + net_blast * dt))

    # VAF dynamics
    net_vaf  <- vaf_drift * VAF[i-1] * (1 - VAF[i-1] / 60) -
                p$eff_vaf * eff * VAF[i-1] * 0.5
    VAF[i]   <- max(1, min(60, VAF[i-1] + net_vaf * dt))

    # Methylation
    net_met <- met_drift - p$eff_met * eff * Met[i-1]
    Met[i]  <- max(0.1, min(0.9, Met[i-1] + net_met * dt))

    # Hgb — erythropoiesis suppressed by blasts, improved by treatment
    hgb_target <- 12 - Blast[i] * 0.3
    net_hgb    <- -hgb_decay * (Hgb[i-1] - hgb_target * 0.5) +
                   p$eff_hgb * eff * (hgb_target - Hgb[i-1])
    Hgb[i]     <- max(4, min(14, Hgb[i-1] + net_hgb * dt))

    # Plt
    plt_target <- 200 - Blast[i] * 5
    net_plt    <- -plt_decay * (Plt[i-1] - plt_target * 0.3) +
                   p$eff_plt * eff * (plt_target - Plt[i-1]) * 0.02
    Plt[i]     <- max(5, min(500, Plt[i-1] + net_plt * dt))

    # ANC
    anc_target <- 3 - Blast[i] * 0.1
    net_anc    <- -anc_decay * (ANC[i-1] - anc_target * 0.3) +
                   p$eff_anc * eff * (anc_target - ANC[i-1]) * 0.1
    ANC[i]     <- max(0.05, min(10, ANC[i-1] + net_anc * dt))
  }

  data.frame(
    time        = t,
    Blast       = Blast,
    VAF         = VAF,
    Hgb         = Hgb,
    Plt         = Plt,
    ANC         = ANC,
    Methylation = Met,
    Treatment   = treatment
  )
}

# Determine response category
response_category <- function(sim_df, blast_init, hgb_init) {
  last  <- tail(sim_df, 1)
  early <- sim_df[sim_df$time <= 60, ]
  blast_nadir <- min(early$Blast, na.rm = TRUE)

  if (blast_nadir < 5 && last$Blast < 5)       return("CR (Complete Remission)")
  if (last$Blast < blast_init * 0.5)            return("PR (Partial Remission)")
  if (last$Blast <= blast_init * 1.25)          return("SD (Stable Disease)")
  return("PD (Progressive Disease)")
}

# Transfusion burden simulation
simulate_transfusion <- function(sim_df) {
  # One unit RBC given when Hgb < 8 g/dL
  monthly <- sim_df %>%
    mutate(month = ceiling(time / 30)) %>%
    group_by(month) %>%
    summarise(
      avg_hgb  = mean(Hgb, na.rm = TRUE),
      units_rbc = sum(Hgb < 8) * 2 / 30,  # ~2 units per episode, rough
      .groups = "drop"
    ) %>%
    mutate(units_rbc = round(units_rbc))
  monthly
}

# =============================================================================
# Biomarker Simulation
# =============================================================================
simulate_biomarkers <- function(treatment, days = 365, sf3b1 = TRUE,
                                 gdf11_init = 5.0, hepcidin_init = 80,
                                 iron_init = 500) {
  dt <- 1
  t  <- seq(0, days, by = dt)
  n  <- length(t)

  # Luspatercept suppresses GDF11; AZA/DEC reduce iron loading
  lusp_on <- treatment %in% c("Luspatercept")
  aza_dec  <- treatment %in% c("AZA", "DEC", "Oral-DEC")
  cycle_gdf11 <- 21
  dur_gdf11   <- 1

  GDF11    <- numeric(n);  GDF11[1]   <- gdf11_init
  Hepcidin <- numeric(n);  Hepcidin[1]<- hepcidin_init
  Iron     <- numeric(n);  Iron[1]    <- iron_init
  RS_pct   <- numeric(n);  RS_pct[1]  <- if (sf3b1) 30 else 5
  IPSSM_tr <- numeric(n);  IPSSM_tr[1]<- 2.5

  for (i in 2:n) {
    day       <- t[i]
    drug_on_g <- lusp_on && ((day %% cycle_gdf11) < dur_gdf11 + 14)

    # GDF11: elevated in RS-MDS, suppressed by luspatercept
    net_gdf11 <- if (drug_on_g) -0.03 * GDF11[i-1]
                 else 0.005 * (gdf11_init - GDF11[i-1])
    GDF11[i]  <- max(0.5, GDF11[i-1] + net_gdf11)

    # Hepcidin: rises with transfusion iron loading
    transfusion_boost <- if (Hepcidin[i-1] < 30) 0.05 else 0
    net_hep  <- -0.002 * (Hepcidin[i-1] - 30) + transfusion_boost -
                (if (aza_dec) 0.005 * Hepcidin[i-1] else 0)
    Hepcidin[i] <- max(5, min(400, Hepcidin[i-1] + net_hep))

    # Iron stores: accumulate ~250 mg per RBC unit; chelation reduces
    iron_input <- 1.5   # mg/day from chronic transfusion (approx)
    chelation  <- if (Iron[i-1] > 1000) 2.5 else 0
    Iron[i]    <- max(0, Iron[i-1] + iron_input - chelation)

    # Ring sideroblast % (SF3B1 driven, reduced by luspatercept)
    rs_effect  <- if (drug_on_g) -0.05 * RS_pct[i-1] else 0.001
    RS_pct[i]  <- max(0, min(80, RS_pct[i-1] + rs_effect))

    # IPSS-M trajectory
    net_ipssm   <- 0.001 - (if (aza_dec || lusp_on) 0.002 else 0)
    IPSSM_tr[i] <- max(0.5, IPSSM_tr[i-1] + net_ipssm)
  }

  data.frame(
    time     = t,
    GDF11    = GDF11,
    Hepcidin = Hepcidin,
    Iron     = Iron,
    RS_pct   = RS_pct,
    IPSSM    = IPSSM_tr,
    Treatment = treatment
  )
}

# =============================================================================
# Treatment list
# =============================================================================
treatment_choices <- c(
  "BSC (Best Supportive Care)"             = "BSC",
  "AZA 75 mg/m² SC d1-7 q28d"             = "AZA",
  "DEC 20 mg/m² IV d1-5 q28d"             = "DEC",
  "Oral-DEC/Cedazuridine 35/100 mg d1-5"  = "Oral-DEC",
  "Lenalidomide 10 mg QD d1-21 q28d"      = "Lenalidomide",
  "Luspatercept 1.0 mg/kg SC q21d"        = "Luspatercept",
  "VEN 400mg QD + AZA 75 mg/m² d1-7"     = "VEN+AZA"
)

mutation_choices <- c("SF3B1", "SRSF2", "U2AF1", "TET2", "DNMT3A",
                       "ASXL1", "TP53", "RUNX1")

who_subtypes <- c("MDS-LB (Low Blasts)", "MDS-IB1 (Increased Blasts-1)",
                  "MDS-IB2 (Increased Blasts-2)", "MDS-RS (Ring Sideroblasts)",
                  "MDS-5q (del 5q)", "MDS-EB (Excess Blasts)")

# =============================================================================
# UI
# =============================================================================
ui <- page_navbar(
  title = "MDS QSP Dashboard",
  theme = bs_theme(
    bootswatch = "darkly",
    primary    = "#4e9af1",
    font_scale = 0.9
  ),
  fillable = TRUE,

  # ── Tab 1: Patient Profile ──────────────────────────────────────────────────
  nav_panel(
    title = icon_text("Patient Profile"),
    layout_sidebar(
      sidebar = sidebar(
        width = 310,
        h5("IPSS-R Calculator", class = "text-primary"),
        sliderInput("bm_blast",  "BM Blast (%)",        0,  20, 5,   step = 0.5),
        selectInput("cyto_cat",  "Cytogenetic Category",
                    choices = c("Very Good", "Good", "Intermediate", "Poor", "Very Poor"),
                    selected = "Intermediate"),
        sliderInput("hgb_base",  "Hemoglobin (g/dL)",   4,  15, 8.5, step = 0.1),
        sliderInput("plt_base",  "Platelets (×10⁹/L)", 0, 500, 80,  step = 5),
        sliderInput("anc_base",  "ANC (×10⁹/L)",       0,   5, 1.0, step = 0.1),
        hr(),
        h5("WHO 2022 Subtype", class = "text-primary"),
        selectInput("who_subtype", NULL, choices = who_subtypes),
        hr(),
        h5("Mutation Profile", class = "text-primary"),
        checkboxGroupInput("mutations", NULL,
                           choices  = mutation_choices,
                           selected = c("TET2", "DNMT3A"),
                           inline   = FALSE)
      ),
      # main panel
      layout_columns(
        col_widths = c(4, 4, 4),
        # IPSS-R score card
        card(
          card_header("IPSS-R Score"),
          card_body(
            uiOutput("ipssr_display")
          )
        ),
        # IPSS-M score card
        card(
          card_header("IPSS-M Estimate"),
          card_body(
            uiOutput("ipssm_display")
          )
        ),
        # WHO subtype card
        card(
          card_header("WHO 2022 Classification"),
          card_body(
            uiOutput("who_display")
          )
        )
      ),
      card(
        card_header("Patient Summary"),
        card_body(
          DTOutput("patient_table")
        )
      )
    )
  ),

  # ── Tab 2: Drug PK ──────────────────────────────────────────────────────────
  nav_panel(
    title = icon_text("Drug PK"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        h5("Drug Selection", class = "text-primary"),
        selectInput("pk_drug",    "Drug",
                    choices  = names(pk_params),
                    selected = "AZA"),
        sliderInput("pk_cycles",  "Number of Cycles", 1, 12, 3),
        sliderInput("pk_weight",  "Body Weight (kg)", 40, 120, 70),
        sliderInput("pk_height",  "Height (cm)",     140, 200, 170),
        hr(),
        uiOutput("bsa_display"),
        hr(),
        h6("Reference PK Parameters", class = "text-info"),
        DTOutput("pk_ref_table")
      ),
      card(
        card_header("PK Concentration–Time Profile"),
        card_body(
          plotlyOutput("pk_plot", height = "340px")
        )
      ),
      card(
        card_header("PK Metrics (Cycle 1)"),
        card_body(
          DTOutput("pk_metrics_table")
        )
      )
    )
  ),

  # ── Tab 3: Disease Dynamics ─────────────────────────────────────────────────
  nav_panel(
    title = icon_text("Disease Dynamics"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        h5("Initial Disease State", class = "text-primary"),
        sliderInput("dd_blast",  "Initial BM Blast (%)",      0.5, 20,  10,   step = 0.5),
        sliderInput("dd_vaf",    "Initial Clonal VAF (%)",     1,   60,  30,   step = 1),
        sliderInput("dd_met",    "Initial Methylation Index",  0.1,  0.9, 0.65, step = 0.05),
        hr(),
        h5("Treatment", class = "text-primary"),
        selectInput("dd_trt", NULL,
                    choices = treatment_choices, selected = "AZA"),
        sliderInput("dd_horizon", "Time Horizon (days)", 30, 730, 365, step = 30),
        hr(),
        uiOutput("response_badge")
      ),
      layout_columns(
        col_widths = c(4, 4, 4),
        card(card_header("BM Blast %"),
             card_body(plotlyOutput("blast_plot", height = "250px"))),
        card(card_header("Clonal VAF (%)"),
             card_body(plotlyOutput("vaf_plot",   height = "250px"))),
        card(card_header("DNA Methylation Index"),
             card_body(plotlyOutput("met_plot",   height = "250px")))
      )
    )
  ),

  # ── Tab 4: Hematologic Endpoints ────────────────────────────────────────────
  nav_panel(
    title = icon_text("Hematologic Endpoints"),
    layout_sidebar(
      sidebar = sidebar(
        width = 260,
        h5("Simulation Settings", class = "text-primary"),
        selectInput("he_trt", "Treatment",
                    choices = treatment_choices, selected = "AZA"),
        sliderInput("he_horizon", "Time Horizon (days)", 60, 730, 365, step = 30),
        sliderInput("he_hgb",  "Initial Hgb (g/dL)",  4,   12, 8.5, step = 0.1),
        sliderInput("he_plt",  "Initial PLT (×10⁹/L)", 5, 400, 80,  step = 5),
        sliderInput("he_anc",  "Initial ANC (×10⁹/L)", 0.1, 5.0, 1.0, step = 0.1),
        sliderInput("he_blast","Initial Blast (%)",    0.5,  20, 10,  step = 0.5),
        hr(),
        h5("Hematologic Improvement Criteria", class = "text-info"),
        uiOutput("hi_summary")
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Hemoglobin (g/dL)"),
             card_body(plotlyOutput("hgb_plot",  height = "240px"))),
        card(card_header("Platelets (×10⁹/L)"),
             card_body(plotlyOutput("plt_plot",  height = "240px"))),
        card(card_header("ANC (×10⁹/L)"),
             card_body(plotlyOutput("anc_plot",  height = "240px"))),
        card(card_header("Monthly Transfusion Burden (RBC units)"),
             card_body(plotlyOutput("transfusion_plot", height = "240px")))
      )
    )
  ),

  # ── Tab 5: Scenario Comparison ──────────────────────────────────────────────
  nav_panel(
    title = icon_text("Scenario Comparison"),
    layout_sidebar(
      sidebar = sidebar(
        width = 260,
        h5("Base Patient Parameters", class = "text-primary"),
        sliderInput("sc_blast",   "Initial Blast (%)",     1, 20, 10,  step = 1),
        sliderInput("sc_hgb",     "Initial Hgb (g/dL)",   4, 12, 8.5, step = 0.5),
        sliderInput("sc_horizon", "Time Horizon (days)", 30, 730, 365, step = 30),
        hr(),
        h5("Scenarios", class = "text-info"),
        p("All 7 treatment arms simulated simultaneously", class = "text-muted small")
      ),
      card(
        card_header("Hgb over Time — All Scenarios"),
        card_body(plotlyOutput("sc_hgb_plot",   height = "280px"))
      ),
      card(
        card_header("BM Blasts over Time — All Scenarios"),
        card_body(plotlyOutput("sc_blast_plot", height = "280px"))
      ),
      card(
        card_header("AML Transformation Risk over Time"),
        card_body(plotlyOutput("sc_aml_plot",   height = "280px"))
      ),
      card(
        card_header("Summary Comparison Table"),
        card_body(DTOutput("sc_summary_table"))
      )
    )
  ),

  # ── Tab 6: Biomarker Tracker ─────────────────────────────────────────────────
  nav_panel(
    title = icon_text("Biomarker Tracker"),
    layout_sidebar(
      sidebar = sidebar(
        width = 270,
        h5("Biomarker Settings", class = "text-primary"),
        selectInput("bm_trt", "Primary Treatment",
                    choices  = treatment_choices,
                    selected = "Luspatercept"),
        checkboxInput("bm_overlay_bsc", "Overlay BSC for comparison", TRUE),
        checkboxInput("bm_sf3b1",       "SF3B1 mutation present",     TRUE),
        sliderInput("bm_gdf11",   "Initial GDF11 (ng/mL)",   0.5, 15, 5,   step = 0.5),
        sliderInput("bm_hep",     "Initial Hepcidin (ng/mL)", 10, 300, 80,  step = 10),
        sliderInput("bm_iron",    "Initial Iron Stores (mg)",  0, 3000, 500, step = 100),
        sliderInput("bm_horizon", "Time Horizon (days)",       30, 730, 365, step = 30),
        hr(),
        p("GDF11 elevated in MDS-RS; suppressed by luspatercept activin-receptor ligand trap.",
          class = "text-muted small")
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("GDF11 (ng/mL)"),
             card_body(plotlyOutput("bm_gdf11_plot",   height = "240px"))),
        card(card_header("Hepcidin (ng/mL)"),
             card_body(plotlyOutput("bm_hep_plot",     height = "240px"))),
        card(card_header("Iron Stores (mg)"),
             card_body(plotlyOutput("bm_iron_plot",    height = "240px"))),
        card(card_header("Ring Sideroblast % (SF3B1)"),
             card_body(plotlyOutput("bm_rs_plot",      height = "240px"))),
        card(card_header("IPSS-M Score Trajectory"),
             card_body(plotlyOutput("bm_ipssm_plot",   height = "240px")))
      )
    )
  )
)

# Helper to add icon labels to nav_panel titles
icon_text <- function(label) label  # Simple passthrough; icons need htmltools in full version

# =============================================================================
# Server
# =============================================================================
server <- function(input, output, session) {

  # ── Reactive: core patient data ────────────────────────────────────────────
  ipssr_r <- reactive({
    calc_ipssr(input$bm_blast, input$cyto_cat,
               input$hgb_base, input$plt_base, input$anc_base)
  })

  ipssm_r <- reactive({
    calc_ipssm(input$mutations, ipssr_r()$score)
  })

  # ── Tab 1: Patient Profile ─────────────────────────────────────────────────
  risk_color <- function(rg) {
    switch(rg,
      "Very Low"     = "success",
      "Low"          = "success",
      "Intermediate" = "warning",
      "Moderate Low" = "warning",
      "Moderate High"= "danger",
      "High"         = "danger",
      "Very High"    = "danger",
      "secondary"
    )
  }

  output$ipssr_display <- renderUI({
    r <- ipssr_r()
    tagList(
      tags$h2(r$score, class = "display-5 fw-bold"),
      span(r$risk_group,
           class = paste0("badge bg-", risk_color(r$risk_group), " fs-6"))
    )
  })

  output$ipssm_display <- renderUI({
    r <- ipssm_r()
    tagList(
      tags$h2(r$score, class = "display-5 fw-bold"),
      span(r$risk_group,
           class = paste0("badge bg-", risk_color(r$risk_group), " fs-6"))
    )
  })

  output$who_display <- renderUI({
    tagList(
      tags$h5(input$who_subtype, class = "mt-2"),
      tags$p(paste(length(input$mutations), "driver mutation(s) detected"),
             class = "text-muted small"),
      tags$p(paste("Mutations:", paste(input$mutations, collapse = ", ")),
             class = "small")
    )
  })

  output$patient_table <- renderDT({
    r    <- ipssr_r()
    rm   <- ipssm_r()
    df <- data.frame(
      Parameter = c("BM Blasts (%)", "Cytogenetics", "Hemoglobin (g/dL)",
                    "Platelets (×10⁹/L)", "ANC (×10⁹/L)",
                    "IPSS-R Score", "IPSS-R Risk Group",
                    "IPSS-M Score (est.)", "IPSS-M Risk Group",
                    "WHO 2022 Subtype", "Driver Mutations"),
      Value     = c(input$bm_blast, input$cyto_cat, input$hgb_base,
                    input$plt_base, input$anc_base,
                    r$score, r$risk_group,
                    rm$score, rm$risk_group,
                    input$who_subtype,
                    paste(input$mutations, collapse = ", "))
    )
    datatable(df, options = list(dom = "t", pageLength = 15, scrollX = TRUE),
              rownames = FALSE,
              class = "compact stripe")
  })

  # ── Tab 2: PK ──────────────────────────────────────────────────────────────
  bsa_r <- reactive({
    calc_bsa(input$pk_weight, input$pk_height)
  })

  output$bsa_display <- renderUI({
    tagList(
      tags$p(sprintf("BSA (Mosteller): %.2f m²", bsa_r()),
             class = "fw-bold text-info")
    )
  })

  pk_sim_r <- reactive({
    simulate_pk(input$pk_drug, input$pk_cycles, bsa_r(), input$pk_weight)
  })

  output$pk_plot <- renderPlotly({
    df   <- pk_sim_r()$data
    drug <- input$pk_drug
    p <- ggplot(df, aes(x = time_day, y = conc_ugml)) +
      geom_line(color = "#4e9af1", linewidth = 0.8) +
      labs(x = "Time (days)", y = "Concentration (μg/mL)",
           title = paste(drug, "— Multi-Cycle PK")) +
      theme_minimal(base_size = 11) +
      theme(panel.background = element_rect(fill = "#2b3035", color = NA),
            plot.background  = element_rect(fill = "#2b3035", color = NA),
            text             = element_text(color = "white"),
            axis.text        = element_text(color = "#cccccc"),
            panel.grid       = element_line(color = "#3d4852"))
    ggplotly(p) %>% layout(paper_bgcolor = "#2b3035", plot_bgcolor = "#2b3035",
                           font = list(color = "white"))
  })

  output$pk_metrics_table <- renderDT({
    df <- pk_sim_r()$metrics
    datatable(df, options = list(dom = "t", pageLength = 10),
              rownames = FALSE, class = "compact")
  })

  output$pk_ref_table <- renderDT({
    p  <- pk_params[[input$pk_drug]]
    df <- data.frame(
      Parameter = c("Route", "t½ (h)", "Vd (L)", "CL (L/h)", "F (%)"),
      Value     = c(p$route, p$t_half, p$Vd, p$CL, round(p$F * 100))
    )
    datatable(df, options = list(dom = "t", pageLength = 10),
              rownames = FALSE, class = "compact")
  })

  # ── Tab 3: Disease Dynamics ────────────────────────────────────────────────
  dd_sim_r <- reactive({
    simulate_disease(
      treatment    = input$dd_trt,
      days         = input$dd_horizon,
      blast_init   = input$dd_blast,
      vaf_init     = input$dd_vaf,
      hgb_init     = input$hgb_base,
      plt_init     = input$plt_base,
      anc_init     = input$anc_base,
      methylation_init = input$dd_met
    )
  })

  output$response_badge <- renderUI({
    sim <- dd_sim_r()
    rc  <- response_category(sim, input$dd_blast, input$hgb_base)
    color <- if (grepl("CR", rc)) "success"
             else if (grepl("PR", rc)) "primary"
             else if (grepl("SD", rc)) "warning"
             else "danger"
    tagList(
      h5("Predicted Response"),
      span(rc, class = paste0("badge bg-", color, " fs-6"))
    )
  })

  make_dd_plot <- function(var, ylab, color_hex) {
    sim <- dd_sim_r()
    p <- ggplot(sim, aes(x = time, y = .data[[var]])) +
      geom_line(color = color_hex, linewidth = 1) +
      labs(x = "Day", y = ylab) +
      theme_minimal(base_size = 10) +
      theme(panel.background = element_rect(fill = "#2b3035", color = NA),
            plot.background  = element_rect(fill = "#2b3035", color = NA),
            text             = element_text(color = "white"),
            axis.text        = element_text(color = "#cccccc"),
            panel.grid       = element_line(color = "#3d4852"))
    ggplotly(p) %>% layout(paper_bgcolor = "#2b3035", plot_bgcolor = "#2b3035",
                           font = list(color = "white"))
  }

  output$blast_plot <- renderPlotly(make_dd_plot("Blast",       "BM Blast (%)",  "#ff6b6b"))
  output$vaf_plot   <- renderPlotly(make_dd_plot("VAF",         "VAF (%)",        "#ffd166"))
  output$met_plot   <- renderPlotly(make_dd_plot("Methylation", "Methylation",    "#06d6a0"))

  # ── Tab 4: Hematologic Endpoints ───────────────────────────────────────────
  he_sim_r <- reactive({
    simulate_disease(
      treatment  = input$he_trt,
      days       = input$he_horizon,
      blast_init = input$he_blast,
      hgb_init   = input$he_hgb,
      plt_init   = input$he_plt,
      anc_init   = input$he_anc
    )
  })

  make_he_plot <- function(var, ylab, ref_line = NULL, color_hex = "#4e9af1") {
    sim <- he_sim_r()
    p <- ggplot(sim, aes(x = time, y = .data[[var]])) +
      geom_line(color = color_hex, linewidth = 1)
    if (!is.null(ref_line))
      p <- p + geom_hline(yintercept = ref_line, linetype = "dashed",
                          color = "#ff6b6b", alpha = 0.7)
    p <- p + labs(x = "Day", y = ylab) +
      theme_minimal(base_size = 10) +
      theme(panel.background = element_rect(fill = "#2b3035", color = NA),
            plot.background  = element_rect(fill = "#2b3035", color = NA),
            text             = element_text(color = "white"),
            axis.text        = element_text(color = "#cccccc"),
            panel.grid       = element_line(color = "#3d4852"))
    ggplotly(p) %>% layout(paper_bgcolor = "#2b3035", plot_bgcolor = "#2b3035",
                           font = list(color = "white"))
  }

  output$hgb_plot  <- renderPlotly(make_he_plot("Hgb", "Hgb (g/dL)", ref_line = 8, color_hex = "#ff6b6b"))
  output$plt_plot  <- renderPlotly(make_he_plot("Plt", "PLT (×10⁹/L)", ref_line = 50, color_hex = "#ffd166"))
  output$anc_plot  <- renderPlotly(make_he_plot("ANC", "ANC (×10⁹/L)", ref_line = 0.5, color_hex = "#06d6a0"))

  output$transfusion_plot <- renderPlotly({
    sim <- he_sim_r()
    monthly <- simulate_transfusion(sim)
    p <- ggplot(monthly, aes(x = month, y = units_rbc)) +
      geom_bar(stat = "identity", fill = "#e63946", alpha = 0.8) +
      labs(x = "Month", y = "RBC Units") +
      theme_minimal(base_size = 10) +
      theme(panel.background = element_rect(fill = "#2b3035", color = NA),
            plot.background  = element_rect(fill = "#2b3035", color = NA),
            text             = element_text(color = "white"),
            axis.text        = element_text(color = "#cccccc"),
            panel.grid       = element_line(color = "#3d4852"))
    ggplotly(p) %>% layout(paper_bgcolor = "#2b3035", plot_bgcolor = "#2b3035",
                           font = list(color = "white"))
  })

  output$hi_summary <- renderUI({
    sim <- he_sim_r()
    hgb_rise  <- max(sim$Hgb) - sim$Hgb[1]
    plt_rise  <- (max(sim$Plt) - sim$Plt[1]) / sim$Plt[1] * 100
    anc_rise  <- max(sim$ANC) - sim$ANC[1]
    transfusion_free <- mean(sim$Hgb > 9) * 100

    hi_e <- if (hgb_rise >= 1.5)  paste0("HI-E: YES (+", round(hgb_rise, 1), " g/dL)")
            else                   paste0("HI-E: NO (+",  round(hgb_rise, 1), " g/dL)")
    hi_p <- if (plt_rise  >= 50)  paste0("HI-P: YES (+", round(plt_rise,  0), "%)")
            else                   paste0("HI-P: NO (+",  round(plt_rise,  0), "%)")
    hi_n <- if (anc_rise  >= 0.1) paste0("HI-N: YES (+", round(anc_rise,  2), " ×10⁹/L)")
            else                   paste0("HI-N: NO (+",  round(anc_rise,  2), " ×10⁹/L)")

    tagList(
      tags$p(hi_e, class = if (hgb_rise >= 1.5) "text-success" else "text-danger"),
      tags$p(hi_p, class = if (plt_rise >= 50)  "text-success" else "text-danger"),
      tags$p(hi_n, class = if (anc_rise >= 0.1) "text-success" else "text-danger"),
      tags$p(sprintf("TI probability: %.0f%%", transfusion_free),
             class = "text-info fw-bold")
    )
  })

  # ── Tab 5: Scenario Comparison ─────────────────────────────────────────────
  sc_all_r <- reactive({
    trt_codes <- unname(treatment_choices)
    results <- lapply(trt_codes, function(trt) {
      simulate_disease(
        treatment  = trt,
        days       = input$sc_horizon,
        blast_init = input$sc_blast,
        hgb_init   = input$sc_hgb,
        plt_init   = 80,
        anc_init   = 1.0
      )
    })
    df <- bind_rows(results)
    # AML risk proxy: blast^1.5 / 100
    df$AML_risk <- pmin(100, df$Blast^1.5)
    df
  })

  sc_colors <- c(
    "BSC"          = "#8d99ae",
    "AZA"          = "#4e9af1",
    "DEC"          = "#06d6a0",
    "Oral-DEC"     = "#26c6da",
    "Lenalidomide" = "#ffd166",
    "Luspatercept" = "#ff6b6b",
    "VEN+AZA"      = "#e040fb"
  )

  make_sc_plot <- function(var, ylab) {
    df <- sc_all_r()
    p  <- ggplot(df, aes(x = time, y = .data[[var]],
                          color = Treatment, group = Treatment)) +
      geom_line(linewidth = 0.9) +
      scale_color_manual(values = sc_colors) +
      labs(x = "Day", y = ylab, color = NULL) +
      theme_minimal(base_size = 10) +
      theme(panel.background = element_rect(fill = "#2b3035", color = NA),
            plot.background  = element_rect(fill = "#2b3035", color = NA),
            text             = element_text(color = "white"),
            axis.text        = element_text(color = "#cccccc"),
            panel.grid       = element_line(color = "#3d4852"),
            legend.background= element_rect(fill = "#2b3035", color = NA))
    ggplotly(p) %>%
      layout(paper_bgcolor = "#2b3035", plot_bgcolor = "#2b3035",
             font = list(color = "white"),
             legend = list(bgcolor = "#2b3035"))
  }

  output$sc_hgb_plot   <- renderPlotly(make_sc_plot("Hgb",      "Hgb (g/dL)"))
  output$sc_blast_plot <- renderPlotly(make_sc_plot("Blast",    "BM Blast (%)"))
  output$sc_aml_plot   <- renderPlotly(make_sc_plot("AML_risk", "AML Risk Index"))

  output$sc_summary_table <- renderDT({
    df <- sc_all_r()
    summary_df <- df %>%
      group_by(Treatment) %>%
      summarise(
        `Median Hgb (g/dL)`    = round(median(Hgb,   na.rm = TRUE), 1),
        `Max Hgb (g/dL)`       = round(max(Hgb,      na.rm = TRUE), 1),
        `Min Blast (%)`        = round(min(Blast,    na.rm = TRUE), 1),
        `Final Blast (%)`      = round(last(Blast),  1),
        `AML Risk (end)`       = round(last(AML_risk), 1),
        `TI Days (%)`          = round(mean(Hgb > 9) * 100),
        .groups = "drop"
      ) %>%
      arrange(`Min Blast (%)`)
    datatable(summary_df,
              options = list(dom = "t", pageLength = 10, scrollX = TRUE),
              rownames = FALSE, class = "compact stripe")
  })

  # ── Tab 6: Biomarker Tracker ─────────────────────────────────────────────────
  bm_sim_r <- reactive({
    sim1 <- simulate_biomarkers(
      treatment   = input$bm_trt,
      days        = input$bm_horizon,
      sf3b1       = input$bm_sf3b1,
      gdf11_init  = input$bm_gdf11,
      hepcidin_init = input$bm_hep,
      iron_init   = input$bm_iron
    )
    if (input$bm_overlay_bsc) {
      sim2 <- simulate_biomarkers(
        treatment   = "BSC",
        days        = input$bm_horizon,
        sf3b1       = input$bm_sf3b1,
        gdf11_init  = input$bm_gdf11,
        hepcidin_init = input$bm_hep,
        iron_init   = input$bm_iron
      )
      return(bind_rows(sim1, sim2))
    }
    sim1
  })

  bm_palette <- c("BSC" = "#8d99ae", "AZA" = "#4e9af1", "DEC" = "#06d6a0",
                  "Oral-DEC" = "#26c6da", "Lenalidomide" = "#ffd166",
                  "Luspatercept" = "#ff6b6b", "VEN+AZA" = "#e040fb")

  make_bm_plot <- function(var, ylab, ref_line = NULL) {
    df <- bm_sim_r()
    p  <- ggplot(df, aes(x = time, y = .data[[var]],
                          color = Treatment, group = Treatment)) +
      geom_line(linewidth = 0.9) +
      scale_color_manual(values = bm_palette) +
      labs(x = "Day", y = ylab, color = NULL)
    if (!is.null(ref_line))
      p <- p + geom_hline(yintercept = ref_line, linetype = "dashed",
                          color = "#ffffff50")
    p <- p +
      theme_minimal(base_size = 10) +
      theme(panel.background = element_rect(fill = "#2b3035", color = NA),
            plot.background  = element_rect(fill = "#2b3035", color = NA),
            text             = element_text(color = "white"),
            axis.text        = element_text(color = "#cccccc"),
            panel.grid       = element_line(color = "#3d4852"),
            legend.background= element_rect(fill = "#2b3035", color = NA))
    ggplotly(p) %>%
      layout(paper_bgcolor = "#2b3035", plot_bgcolor = "#2b3035",
             font = list(color = "white"),
             legend = list(bgcolor = "#2b3035"))
  }

  output$bm_gdf11_plot  <- renderPlotly(make_bm_plot("GDF11",    "GDF11 (ng/mL)"))
  output$bm_hep_plot    <- renderPlotly(make_bm_plot("Hepcidin", "Hepcidin (ng/mL)"))
  output$bm_iron_plot   <- renderPlotly(make_bm_plot("Iron",     "Iron Stores (mg)", ref_line = 1000))
  output$bm_rs_plot     <- renderPlotly(make_bm_plot("RS_pct",   "Ring Sideroblast (%)"))
  output$bm_ipssm_plot  <- renderPlotly(make_bm_plot("IPSSM",    "IPSS-M Score"))
}

# =============================================================================
# Launch
# =============================================================================
shinyApp(ui = ui, server = server)
