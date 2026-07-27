## ============================================================================
## Anemia of Chronic Kidney Disease (ACKD) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 8 tabs:
##   1 Patient profile      — build the virtual CKD patient, see the derived
##                            baseline (REP pool, EPO, iron fluxes) before dosing
##   2 Drug PK              — ESA and HIF-PHI exposure profiles by agent/schedule
##   3 Erythropoiesis PD    — HIF -> EPO -> progenitor -> reticulocyte -> RBC chain
##   4 Iron & hepcidin      — TSAT, ferritin, hepcidin, ERFE, the iron gate and
##                            the full mg/day iron flux balance
##   5 Clinical endpoints   — Hb vs the KDIGO 10-11.5 band, FACIT-Fatigue, LVMI
##   6 Scenario comparison  — all 15 prebuilt arms side by side + summary table
##   7 Safety               — Hb overshoot, MAP, thrombotic hazard, NTBI,
##                            FGF23/phosphate after IV iron
##   8 Titration & docs     — closed-loop KDIGO titration simulator + references
##
## Dependencies: shiny, shinydashboard, mrgsolve, dplyr, tidyr, ggplot2, DT
## Run with:  shiny::runApp("ackd_shiny_app.R")
##            (keep ackd_mrgsolve_model.R in the same directory)
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

## ---------------------------------------------------------------- model load
## The model spec lives in ackd_mrgsolve_model.R (this library's naming
## convention). mrgsolve's mread() resolves <model>.cpp / .mod, so load it by
## FILE NAME via mread_file(), falling back to a temp .cpp copy on older
## mrgsolve builds that lack mread_file().
ACKD_SPEC <- "ackd_mrgsolve_model.R"

load_ackd <- function(spec = ACKD_SPEC) {
  stopifnot(file.exists(spec))
  if (exists("mread_file", where = asNamespace("mrgsolve"), inherits = FALSE)) {
    mrgsolve::mread_file(basename(spec), project = dirname(normalizePath(spec)),
                         soloc = tempdir(), quiet = TRUE)
  } else {
    tmp <- file.path(tempdir(), "ackd.cpp")
    file.copy(spec, tmp, overwrite = TRUE)
    mrgsolve::mread("ackd", project = tempdir(), soloc = tempdir(), quiet = TRUE)
  }
}

get_model <- function() {
  if (!exists(".ACKD_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.ACKD_MOD)) {
    assign(".ACKD_MOD", load_ackd(), envir = .GlobalEnv)
  }
  .GlobalEnv$.ACKD_MOD
}

## ------------------------------------------------------------- agent presets
ESA_AGENTS <- list(
  "Epoetin alfa (SC)"            = list(KA_ESA = 0.35, F_ESA = 0.30, KEL_ESA = 1.80),
  "Darbepoetin alfa (SC)"        = list(KA_ESA = 0.25, F_ESA = 0.45, KEL_ESA = 0.65),
  "Methoxy-PEG-epoetin beta (SC)"= list(KA_ESA = 0.08, F_ESA = 0.55, KEL_ESA = 0.13)
)
PHI_AGENTS <- list(
  "Roxadustat"  = list(KA_PHI = 3.0, F_PHI = 0.85, KEL_PHI = 0.55),
  "Daprodustat" = list(KA_PHI = 6.0, F_PHI = 0.65, KEL_PHI = 3.50),
  "Vadadustat"  = list(KA_PHI = 4.0, F_PHI = 0.80, KEL_PHI = 0.70)
)
IRON_IV <- list(
  "Ferric carboxymaltose (FGF23-inducing)" = list(FCM_FLAG = 1),
  "Iron sucrose / derisomaltose"           = list(FCM_FLAG = 0)
)
TIW <- 7 / 3

## --------------------------------------------------------- scenario library
## Each entry: list(param overrides, function(end_d) -> event object)
SCENARIOS <- list(
  "1. Untreated natural history (CKD G4)" = list(
    p = list(EGFR0 = 20, DIALYSIS = 0),
    e = function(d) ev(amt = 0, cmt = "ESA_SC")),
  "2. Epoetin alfa 75 U/kg SC TIW (maintenance)" = list(
    p = c(list(EGFR0 = 20, ESA_WK_UKG = 225), ESA_AGENTS[[1]]),
    e = function(d) ev(amt = 0.5, cmt = "ESA_SC", ii = TIW, addl = floor(d / TIW))),
  "3. Epoetin alfa 150 U/kg SC TIW (correction, overshoots)" = list(
    p = c(list(EGFR0 = 20, ESA_WK_UKG = 450), ESA_AGENTS[[1]]),
    e = function(d) ev(amt = 1.0, cmt = "ESA_SC", ii = TIW, addl = floor(d / TIW))),
  "4. Darbepoetin alfa Q2W" = list(
    p = c(list(EGFR0 = 20, ESA_WK_UKG = 450), ESA_AGENTS[[2]]),
    e = function(d) ev(amt = 4.0, cmt = "ESA_SC", ii = 14, addl = floor(d / 14))),
  "5. Methoxy-PEG-epoetin beta Q4W" = list(
    p = c(list(EGFR0 = 20, ESA_WK_UKG = 450), ESA_AGENTS[[3]]),
    e = function(d) ev(amt = 7.0, cmt = "ESA_SC", ii = 28, addl = floor(d / 28))),
  "6. Epoetin + proactive IV iron (HD, PIVOTAL)" = list(
    p = c(list(EGFR0 = 8, DIALYSIS = 1, ESA_WK_UKG = 300, FCM_FLAG = 0), ESA_AGENTS[[1]]),
    e = function(d) ev(amt = 1.0, cmt = "ESA_SC", ii = TIW, addl = floor(d / TIW)) +
                    ev(amt = 400, cmt = "FE_IV", ii = 28, addl = floor(d / 28))),
  "7. Roxadustat TIW (dialysis)" = list(
    p = c(list(EGFR0 = 8, DIALYSIS = 1), PHI_AGENTS[[1]]),
    e = function(d) ev(amt = 2.0, cmt = "PHI_GUT", ii = TIW, addl = floor(d / TIW))),
  "8. Daprodustat daily (non-dialysis)" = list(
    p = c(list(EGFR0 = 20), PHI_AGENTS[[2]]),
    e = function(d) ev(amt = 4.0, cmt = "PHI_GUT", ii = 1, addl = max(d - 1, 0))),
  "9. Oral iron alone (hepcidin-capped)" = list(
    p = list(EGFR0 = 20),
    e = function(d) ev(amt = 65, cmt = "FE_GUT", ii = 1, addl = max(d - 1, 0))),
  "10. Inflammatory hyporesponder, matched ESA dose" = list(
    p = c(list(EGFR0 = 20, IL6_DRIVE = 4, ESA_WK_UKG = 450), ESA_AGENTS[[1]]),
    e = function(d) ev(amt = 1.0, cmt = "ESA_SC", ii = TIW, addl = floor(d / TIW))),
  "11. Inflammatory hyporesponder, ESA escalated 3x" = list(
    p = c(list(EGFR0 = 20, IL6_DRIVE = 4, ESA_WK_UKG = 1350), ESA_AGENTS[[1]]),
    e = function(d) ev(amt = 3.0, cmt = "ESA_SC", ii = TIW, addl = floor(d / TIW))),
  "12. Inflammatory hyporesponder, roxadustat" = list(
    p = c(list(EGFR0 = 20, IL6_DRIVE = 4), PHI_AGENTS[[1]]),
    e = function(d) ev(amt = 1.5, cmt = "PHI_GUT", ii = TIW, addl = floor(d / TIW))),
  "13. Hb normalisation overshoot (safety)" = list(
    p = c(list(EGFR0 = 20, ESA_WK_UKG = 2250), ESA_AGENTS[[1]]),
    e = function(d) ev(amt = 5.0, cmt = "ESA_SC", ii = TIW, addl = floor(d / TIW)) +
                    ev(amt = 300, cmt = "FE_IV", ii = 28, addl = floor(d / 28))),
  "14. High-dose FCM -> FGF23 -> hypophosphatemia" = list(
    p = list(EGFR0 = 20, FCM_FLAG = 1),
    e = function(d) ev(amt = 750, cmt = "FE_IV", ii = 7, addl = 1)),
  "15. ESA interruption at day 120 (Hb cycling)" = list(
    p = c(list(EGFR0 = 20, ESA_WK_UKG = 450), ESA_AGENTS[[1]]),
    e = function(d) ev(amt = 1.0, cmt = "ESA_SC", ii = TIW,
                       addl = floor(min(120, d) / TIW)))
)

## param() rejects an empty list, so only apply overrides when there are some
with_param <- function(m, p) if (length(p)) param(m, p) else m

run_scenario <- function(name, end_d, extra = list()) {
  sc <- SCENARIOS[[name]]
  get_model() %>%
    with_param(c(sc$p, extra)) %>%
    mrgsim(events = sc$e(end_d), end = end_d, delta = 1) %>%
    as_tibble() %>%
    mutate(scenario = name)
}

## ------------------------------------------------- closed-loop KDIGO titration
titrate <- function(dose0 = 0.5, lo = 10, hi = 11.5, blocks = 13, block_d = 28,
                    p = list(), step = 0.25) {
  base <- get_model() %>% with_param(p)
  cmts <- names(init(get_model()))
  st <- NULL; dose <- dose0; out <- list()
  for (b in seq_len(blocks)) {
    ## [MAIN] re-derives the initial conditions on every call, so a resumed
    ## block MUST set INIT_FROM_PARAM = 0 or init() is silently overwritten.
    m <- if (is.null(st)) base else base %>% param(INIT_FROM_PARAM = 0) %>% init(st)
    e <- if (dose > 0) {
      ev(amt = dose, cmt = "ESA_SC", ii = TIW, addl = floor(block_d / TIW))
    } else ev(amt = 0, cmt = "ESA_SC")
    r  <- m %>% mrgsim(events = e, end = block_d, delta = 1) %>% as_tibble()
    hb <- tail(r$Hemoglobin, 1)
    out[[b]] <- r %>% mutate(block = b, time = time + (b - 1) * block_d,
                             dose_au = dose)
    ## KDIGO-style titration: hold above the ceiling, +-25% steps inside
    dose <- if (hb > hi) max(0, dose * (1 - step)) else
            if (hb < lo) dose * (1 + step) else dose
    st <- as.list(tail(r, 1))[cmts]
  }
  bind_rows(out)
}

## ================================================================== UI ======
ui <- dashboardPage(
  dashboardHeader(title = "Anemia of CKD — QSP", titleWidth = 300),
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("1. Patient profile",     tabName = "patient", icon = icon("user")),
      menuItem("2. Drug PK",             tabName = "pk",      icon = icon("pills")),
      menuItem("3. Erythropoiesis PD",   tabName = "pd",      icon = icon("dna")),
      menuItem("4. Iron & hepcidin",     tabName = "iron",    icon = icon("magnet")),
      menuItem("5. Clinical endpoints",  tabName = "endpoint",icon = icon("heart-pulse")),
      menuItem("6. Scenario comparison", tabName = "compare", icon = icon("layer-group")),
      menuItem("7. Safety",              tabName = "safety",  icon = icon("triangle-exclamation")),
      menuItem("8. Titration & docs",    tabName = "docs",    icon = icon("sliders"))
    ),
    hr(),
    h4("Virtual patient", style = "padding-left:15px"),
    sliderInput("EGFR0", "Baseline eGFR (mL/min/1.73m2)", 5, 60, 20, 1),
    sliderInput("BASE_HB", "Presenting Hb (g/dL)", 7, 12, 9.5, 0.1),
    sliderInput("FIB", "Interstitial fibrosis / REP loss (0-1)", 0, 1, 0.6, 0.05),
    sliderInput("BLUNT", "Blunting of renal hypoxic drive (0-1)", 0, 0.95, 0.7, 0.05),
    sliderInput("IL6_DRIVE", "Inflammatory burden (x baseline IL-6)", 0.5, 6, 1, 0.1),
    sliderInput("URE", "Uremic burden (RBC lifespan)", 0, 1, 0.5, 0.05),
    sliderInput("TSAT0", "Baseline TSAT (%)", 8, 45, 24, 1),
    sliderInput("FERRITIN0", "Baseline ferritin (ng/mL)", 30, 800, 147, 5),
    sliderInput("HEPC_BASE", "Baseline hepcidin-25 (ng/mL)", 10, 250, 60, 5),
    checkboxInput("DIALYSIS", "On maintenance haemodialysis", FALSE),
    checkboxInput("ANTI_IL6", "Anti-IL-6(R) adjunct on board", FALSE),
    hr(),
    sliderInput("end_d", "Simulation horizon (days)", 90, 730, 365, 5)
  ),
  dashboardBody(
    tags$style(HTML(".content-wrapper{background:#f7f9fb}
                     .small-note{font-size:12px;color:#666}")),
    tabItems(

      ## ---------------------------------------------------------- 1 patient
      tabItem("patient",
        fluidRow(
          box(width = 12, status = "primary", solidHeader = TRUE,
              title = "The virtual CKD patient and its self-calibrated baseline",
              p("The model back-calculates its production constants from the ",
                "descriptors on the left, so every patient you build STARTS IN ",
                "BALANCE. Only eGFR (and therefore the REP pool and EPO) drifts ",
                "in the untreated arm."),
              p(class = "small-note",
                "Caveat carried through from the model file: because baseline Hb ",
                "is pinned to your slider, raising the inflammatory burden does ",
                "NOT lower the starting Hb. Inflammation shows up as a smaller Hb ",
                "gain per unit of ESA and as a tighter iron gate. Compare ",
                "hyporesponders at MATCHED DOSE (scenario 10 vs 3)."))),
        fluidRow(
          valueBoxOutput("vb_hb", 3), valueBoxOutput("vb_epo", 3),
          valueBoxOutput("vb_rep", 3), valueBoxOutput("vb_hepc", 3)),
        fluidRow(
          valueBoxOutput("vb_tsat", 3), valueBoxOutput("vb_ferr", 3),
          valueBoxOutput("vb_ret", 3), valueBoxOutput("vb_facit", 3)),
        fluidRow(
          box(width = 6, title = "Untreated natural history", status = "info",
              solidHeader = TRUE, plotOutput("p_nat", height = 340)),
          box(width = 6, title = "Baseline state (all compartments)",
              status = "info", solidHeader = TRUE,
              DTOutput("t_baseline")))),

      ## --------------------------------------------------------------- 2 PK
      tabItem("pk",
        fluidRow(
          box(width = 12, status = "primary", solidHeader = TRUE,
              title = "Drug exposure",
              fluidRow(
                column(4, selectInput("esa_agent", "ESA molecule", names(ESA_AGENTS))),
                column(4, numericInput("esa_dose", "ESA dose (au; 1.0 = 150 U/kg epoetin)",
                                       0.5, 0, 10, 0.1)),
                column(4, numericInput("esa_ii", "ESA interval (days)", 2.333, 1, 28, 0.333))),
              fluidRow(
                column(4, selectInput("phi_agent", "HIF-PHI molecule", names(PHI_AGENTS))),
                column(4, numericInput("phi_dose", "HIF-PHI dose (au; 1.0 = 100 mg roxadustat)",
                                       0, 0, 10, 0.5)),
                column(4, numericInput("phi_ii", "HIF-PHI interval (days)", 2.333, 1, 7, 0.333))))),
        fluidRow(
          box(width = 6, title = "ESA plasma exposure (au)", status = "info",
              solidHeader = TRUE, plotOutput("p_esa_pk", height = 320)),
          box(width = 6, title = "HIF-PHI plasma exposure (au)", status = "info",
              solidHeader = TRUE, plotOutput("p_phi_pk", height = 320))),
        fluidRow(
          box(width = 12, title = "Total erythropoietic signal (endogenous EPO + ESA equivalent)",
              status = "info", solidHeader = TRUE, plotOutput("p_epo_sig", height = 300),
              p(class = "small-note",
                "Note the contrast the trials show: ESAs produce large supra-",
                "physiological EPO peaks, HIF-PHIs only a 3-5x rise."))) ),

      ## --------------------------------------------------------------- 3 PD
      tabItem("pd",
        fluidRow(
          box(width = 12, status = "primary", solidHeader = TRUE,
              title = "Custom regimen (uses the ESA / HIF-PHI settings from tab 2)",
              checkboxGroupInput("pd_add", NULL, inline = TRUE,
                choices = c("Oral iron 65 mg/day" = "po_fe",
                            "IV iron 400 mg q4w"  = "iv_fe"),
                selected = character(0)),
              selectInput("iv_agent", "IV iron product", names(IRON_IV)))),
        fluidRow(
          box(width = 6, title = "HIF-2alpha activity", status = "info",
              solidHeader = TRUE, plotOutput("p_hif", height = 280)),
          box(width = 6, title = "REP-cell pool & endogenous EPO", status = "info",
              solidHeader = TRUE, plotOutput("p_rep", height = 280))),
        fluidRow(
          box(width = 6, title = "Marrow chain: progenitors -> erythroblasts",
              status = "info", solidHeader = TRUE, plotOutput("p_marrow", height = 280)),
          box(width = 6, title = "Reticulocytes (the early response signal)",
              status = "info", solidHeader = TRUE, plotOutput("p_retic", height = 280))),
        fluidRow(
          box(width = 12, title = "Haemoglobin and the EPOR-signalling attenuation",
              status = "info", solidHeader = TRUE, plotOutput("p_hb_pd", height = 300)))),

      ## ------------------------------------------------------------- 4 iron
      tabItem("iron",
        fluidRow(
          box(width = 12, status = "primary", solidHeader = TRUE,
              title = "The hepcidin-ferroportin-erythroferrone axis",
              p("Baseline iron balance is closed by construction: obligate loss ",
                "equals baseline absorption, and recycling equals erythron ",
                "utilisation (~11 mg/day). Any deficit you see is mechanistic — ",
                "dialysis/GI loss, hepcidin trapping, or an ESA-driven demand surge."))),
        fluidRow(
          box(width = 6, title = "TSAT and the iron-availability gate",
              status = "info", solidHeader = TRUE, plotOutput("p_tsat", height = 300)),
          box(width = 6, title = "Ferritin & storage iron", status = "info",
              solidHeader = TRUE, plotOutput("p_ferr", height = 300))),
        fluidRow(
          box(width = 6, title = "Hepcidin-25 and erythroferrone", status = "info",
              solidHeader = TRUE, plotOutput("p_hepc", height = 300)),
          box(width = 6, title = "TIBC / transferrin", status = "info",
              solidHeader = TRUE, plotOutput("p_tibc", height = 300))),
        fluidRow(
          box(width = 12, title = "Functional vs absolute iron deficiency map",
              status = "info", solidHeader = TRUE, plotOutput("p_iron_map", height = 340),
              p(class = "small-note",
                "The KDIGO iron-sufficiency corner is TSAT >20% AND ferritin ",
                ">100 ng/mL (>200 on dialysis). Trajectories that drop out of ",
                "that corner while Hb is still rising are FUNCTIONAL iron ",
                "deficiency created by the drug.")))),

      ## -------------------------------------------------------- 5 endpoints
      tabItem("endpoint",
        fluidRow(
          box(width = 12, title = "Haemoglobin vs the KDIGO target band",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_hb_band", height = 340))),
        fluidRow(
          box(width = 4, title = "FACIT-Fatigue", status = "info",
              solidHeader = TRUE, plotOutput("p_facit", height = 280)),
          box(width = 4, title = "LV mass index", status = "info",
              solidHeader = TRUE, plotOutput("p_lvmi", height = 280)),
          box(width = 4, title = "eGFR trajectory", status = "info",
              solidHeader = TRUE, plotOutput("p_egfr", height = 280))),
        fluidRow(
          box(width = 12, title = "Endpoint summary for the selected regimen",
              status = "info", solidHeader = TRUE, DTOutput("t_endpoint")))),

      ## --------------------------------------------------------- 6 compare
      tabItem("compare",
        fluidRow(
          box(width = 12, status = "primary", solidHeader = TRUE,
              title = "Compare prebuilt scenarios",
              selectizeInput("sc_pick", NULL, choices = names(SCENARIOS),
                             multiple = TRUE, width = "100%",
                             selected = names(SCENARIOS)[c(1, 2, 3, 7, 13)]),
              p(class = "small-note",
                "Scenario patient descriptors override the sidebar sliders, so ",
                "each arm runs on its own intended patient (e.g. arms 6-7 are ",
                "haemodialysis patients)."))),
        fluidRow(
          box(width = 6, title = "Haemoglobin", status = "info",
              solidHeader = TRUE, plotOutput("c_hb", height = 320)),
          box(width = 6, title = "Hepcidin-25", status = "info",
              solidHeader = TRUE, plotOutput("c_hepc", height = 320))),
        fluidRow(
          box(width = 6, title = "TSAT", status = "info",
              solidHeader = TRUE, plotOutput("c_tsat", height = 320)),
          box(width = 6, title = "Cumulative thrombotic hazard", status = "info",
              solidHeader = TRUE, plotOutput("c_thr", height = 320))),
        fluidRow(
          box(width = 12, title = "Scenario summary at end of horizon",
              status = "info", solidHeader = TRUE, DTOutput("c_table")))),

      ## ---------------------------------------------------------- 7 safety
      tabItem("safety",
        fluidRow(
          box(width = 12, status = "danger", solidHeader = TRUE,
              title = "Safety axes",
              p("The Normal-Hematocrit, CHOIR, CREATE and TREAT trials all found ",
                "that pushing Hb toward normal produced no benefit and excess ",
                "thrombotic events. In this model that emerges from viscosity, ",
                "ESA exposure, MAP and labile-iron terms in the hazard."),
              selectizeInput("sf_pick", "Arms to overlay", choices = names(SCENARIOS),
                             multiple = TRUE, width = "100%",
                             selected = names(SCENARIOS)[c(2, 3, 13, 14)]))),
        fluidRow(
          box(width = 6, title = "Hb overshoot above 13 g/dL", status = "warning",
              solidHeader = TRUE, plotOutput("s_over", height = 300)),
          box(width = 6, title = "Mean arterial pressure", status = "warning",
              solidHeader = TRUE, plotOutput("s_map", height = 300))),
        fluidRow(
          box(width = 6, title = "Cumulative thrombotic hazard", status = "danger",
              solidHeader = TRUE, plotOutput("s_thr", height = 300)),
          box(width = 6, title = "Labile (non-transferrin-bound) iron",
              status = "warning", solidHeader = TRUE, plotOutput("s_ntbi", height = 300))),
        fluidRow(
          box(width = 6, title = "FGF23 after IV iron", status = "warning",
              solidHeader = TRUE, plotOutput("s_fgf", height = 300)),
          box(width = 6, title = "Serum phosphate (hypophosphatemia < 2.0 mg/dL)",
              status = "danger", solidHeader = TRUE, plotOutput("s_phos", height = 300)))),

      ## ------------------------------------------------------------ 8 docs
      tabItem("docs",
        fluidRow(
          box(width = 12, status = "primary", solidHeader = TRUE,
              title = "Closed-loop KDIGO titration simulator",
              p("The prebuilt scenarios are open-loop, which is why a correction ",
                "dose held for a year overshoots. This titrator re-evaluates Hb ",
                "every block and applies +-25% dose steps, holding above the ceiling."),
              fluidRow(
                column(3, numericInput("ti_dose", "Starting dose (au)", 0.5, 0, 10, 0.1)),
                column(3, numericInput("ti_lo", "Target floor (g/dL)", 10, 8, 12, 0.1)),
                column(3, numericInput("ti_hi", "Target ceiling (g/dL)", 11.5, 9, 14, 0.1)),
                column(3, numericInput("ti_block", "Block length (days)", 28, 7, 56, 7))),
              actionButton("ti_go", "Run titration", class = "btn-primary"))),
        fluidRow(
          box(width = 8, title = "Titrated Hb trajectory", status = "info",
              solidHeader = TRUE, plotOutput("ti_hb", height = 330)),
          box(width = 4, title = "Dose per block", status = "info",
              solidHeader = TRUE, plotOutput("ti_dose_plot", height = 330))),
        fluidRow(
          box(width = 12, title = "Model documentation", status = "info",
              solidHeader = TRUE, collapsible = TRUE,
              HTML("
<h4>Deliverables in this directory</h4>
<ul>
<li><code>ackd_qsp_model.dot / .svg / .png</code> — mechanistic map, 219 nodes in
    15 clusters</li>
<li><code>ackd_mrgsolve_model.R</code> — 30-compartment ODE model with a
    self-calibrating baseline</li>
<li><code>ackd_shiny_app.R</code> — this dashboard</li>
<li><code>ackd_references.md</code> — annotated bibliography</li>
<li><code>README.md</code> — mechanism narrative and the validation table</li>
</ul>
<h4>What the model is for</h4>
<p>Teaching and hypothesis generation about <b>why</b> renal anemia responds the
way it does: how much of a poor ESA response is EPO deficiency, how much is
hepcidin-driven iron restriction, and how much is cytokine attenuation of EPOR
signalling — and why HIF-PHIs, which act on all three at once, look different
in the trials.</p>
<h4>Key simplifications</h4>
<ul>
<li>Baseline Hb is <b>pinned</b> to the slider for every patient, so
    inflammation blunts the <i>response</i>, not the presenting Hb. Compare
    hyporesponders at matched dose.</li>
<li>Drug exposure is in dose-proportional <code>au</code>, first-order PK, no
    target-mediated disposition.</li>
<li>One lumped HIF activity for HIF-1&alpha; and HIF-2&alpha; across kidney,
    liver and gut.</li>
<li>The thrombotic hazard <b>ranks</b> arms; it is not a calibrated event rate.</li>
</ul>
<p><b>Not for clinical use.</b> Research and education only.</p>"))),
        fluidRow(
          box(width = 12, title = "Anchor trials & mechanism papers", status = "info",
              solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
              HTML("
<ul>
<li>Eschbach 1987 <i>NEJM</i> 316:73 — first correction of ESKD anemia with rHuEPO</li>
<li>Besarab 1998 <i>NEJM</i> 339:584 — Normal Hematocrit Trial</li>
<li>Singh 2006 <i>NEJM</i> 355:2085 (CHOIR) · Dr&uuml;eke 2006 <i>NEJM</i> 355:2071 (CREATE)</li>
<li>Pfeffer 2009 <i>NEJM</i> 361:2019 (TREAT)</li>
<li>Macdougall 2019 <i>NEJM</i> 380:447 (PIVOTAL)</li>
<li>Chen 2019 <i>NEJM</i> 381:1011 &amp; 381:1001 (roxadustat)</li>
<li>Singh 2021 <i>NEJM</i> 385:2325 (ASCEND-D) &amp; 385:2313 (ASCEND-ND)</li>
<li>Eckardt / Chertow 2021 <i>NEJM</i> 384:1601 / 1589 (INNO2VATE / PRO2TECT)</li>
<li>Nemeth 2004 <i>Science</i> 306:2090 &amp; <i>JCI</i> 113:1271 — hepcidin-ferroportin, IL-6</li>
<li>Kautz 2014 <i>Nat Genet</i> 46:678 — erythroferrone</li>
<li>Wolf 2020 <i>JAMA</i> 323:432 (PHOSPHARE-IDA) — FCM, FGF23, hypophosphatemia</li>
<li>KDIGO 2012 Anemia Guideline · Ku 2023 <i>Kidney Int</i> 104:655 (KDIGO Controversies)</li>
</ul>
<p>Full list with PubMed links: <code>ackd_references.md</code></p>")))
      )
    )
  )
)

## ============================================================== SERVER ======
server <- function(input, output, session) {

  patient <- reactive({
    list(EGFR0 = input$EGFR0, BASE_HB = input$BASE_HB, FIB = input$FIB,
         BLUNT = input$BLUNT, IL6_DRIVE = input$IL6_DRIVE, URE = input$URE,
         TSAT0 = input$TSAT0, FERRITIN0 = input$FERRITIN0,
         HEPC_BASE = input$HEPC_BASE,
         DIALYSIS = as.numeric(input$DIALYSIS),
         ANTI_IL6 = as.numeric(input$ANTI_IL6))
  })

  ## untreated natural history for the selected patient
  natural <- reactive({
    get_model() %>% param(patient()) %>%
      mrgsim(events = ev(amt = 0, cmt = "ESA_SC"), end = input$end_d, delta = 1) %>%
      as_tibble()
  })

  ## custom regimen assembled from tabs 2-3
  custom <- reactive({
    p <- c(patient(), ESA_AGENTS[[input$esa_agent]], PHI_AGENTS[[input$phi_agent]],
           IRON_IV[[input$iv_agent]])
    e <- ev(amt = 0, cmt = "ESA_SC")
    if (input$esa_dose > 0)
      e <- e + ev(amt = input$esa_dose, cmt = "ESA_SC", ii = input$esa_ii,
                  addl = floor(input$end_d / input$esa_ii))
    if (input$phi_dose > 0)
      e <- e + ev(amt = input$phi_dose, cmt = "PHI_GUT", ii = input$phi_ii,
                  addl = floor(input$end_d / input$phi_ii))
    if ("po_fe" %in% input$pd_add)
      e <- e + ev(amt = 65, cmt = "FE_GUT", ii = 1, addl = max(input$end_d - 1, 0))
    if ("iv_fe" %in% input$pd_add)
      e <- e + ev(amt = 400, cmt = "FE_IV", ii = 28, addl = floor(input$end_d / 28))
    get_model() %>% param(p) %>%
      mrgsim(events = e, end = input$end_d, delta = 1) %>% as_tibble()
  })

  compared <- reactive({
    req(length(input$sc_pick) > 0)
    bind_rows(lapply(input$sc_pick, run_scenario, end_d = input$end_d))
  })
  safety_runs <- reactive({
    req(length(input$sf_pick) > 0)
    bind_rows(lapply(input$sf_pick, run_scenario, end_d = input$end_d))
  })

  ## -------------------------------------------------------------- helpers
  base_theme <- theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          panel.grid.minor = element_blank())

  line1 <- function(d, y, ylab, hl = NULL) {
    g <- ggplot(d, aes(time, .data[[y]])) +
      geom_line(linewidth = 0.9, colour = "#2c6fb5") +
      labs(x = "Day", y = ylab) + base_theme
    if (!is.null(hl)) g <- g + geom_hline(yintercept = hl, linetype = 2,
                                          colour = "grey45")
    g
  }
  lineN <- function(d, y, ylab, hl = NULL) {
    g <- ggplot(d, aes(time, .data[[y]], colour = scenario)) +
      geom_line(linewidth = 0.85) + labs(x = "Day", y = ylab) + base_theme +
      guides(colour = guide_legend(ncol = 1))
    if (!is.null(hl)) g <- g + geom_hline(yintercept = hl, linetype = 2,
                                          colour = "grey45")
    g
  }

  ## ------------------------------------------------------------ 1 patient
  vb <- function(x, lab, icon_name, col) {
    valueBox(x, lab, icon = icon(icon_name), color = col)
  }
  output$vb_hb    <- renderValueBox(vb(sprintf("%.2f g/dL", natural()$Hemoglobin[1]),
                                       "Baseline Hb", "droplet", "red"))
  output$vb_epo   <- renderValueBox(vb(sprintf("%.1f mIU/mL", natural()$EPO_endogenous[1]),
                                       "Endogenous EPO", "bolt", "purple"))
  output$vb_rep   <- renderValueBox(vb(sprintf("%.0f%%", 100 * natural()$REP_pool[1]),
                                       "REP pool vs healthy", "kidneys", "orange"))
  output$vb_hepc  <- renderValueBox(vb(sprintf("%.0f ng/mL", natural()$Hepcidin[1]),
                                       "Hepcidin-25", "lock", "teal"))
  output$vb_tsat  <- renderValueBox(vb(sprintf("%.0f%%", natural()$TSAT_pct[1]),
                                       "TSAT", "magnet", "olive"))
  output$vb_ferr  <- renderValueBox(vb(sprintf("%.0f ng/mL", natural()$Ferritin[1]),
                                       "Ferritin", "box", "olive"))
  output$vb_ret   <- renderValueBox(vb(sprintf("%.2f%%", natural()$Reticulocyte_pct[1]),
                                       "Reticulocytes", "seedling", "maroon"))
  output$vb_facit <- renderValueBox(vb(sprintf("%.0f / 52", natural()$FACIT_fatigue[1]),
                                       "FACIT-Fatigue", "battery-half", "navy"))

  output$p_nat <- renderPlot({
    natural() %>%
      select(time, Hb = Hemoglobin, TSAT = TSAT_pct, Hepcidin, eGFR = eGFR_ml_min) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#2c6fb5") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Day", y = NULL) + base_theme
  })
  output$t_baseline <- renderDT({
    d <- natural()[1, ]
    keep <- c("Hemoglobin", "Hematocrit", "RBC_count", "Reticulocyte_pct",
              "EPO_endogenous", "EPOR_attenuation", "HIF_activity", "REP_pool",
              "Serum_iron", "TSAT_pct", "TIBC_level", "Ferritin",
              "Storage_iron_mg", "Iron_gate", "Hepcidin", "Erythroferrone",
              "IL6_level", "CRP_level", "Phosphate", "LDL_C", "MAP_mmHg",
              "FACIT_fatigue", "LVMI_g_m2", "eGFR_ml_min")
    tibble(Variable = keep,
           Value = round(as.numeric(d[keep]), 3)) %>%
      datatable(rownames = FALSE, options = list(pageLength = 12, dom = "tp"))
  })

  ## ----------------------------------------------------------------- 2 PK
  output$p_esa_pk  <- renderPlot(line1(custom(), "ESA_conc", "ESA (au)"))
  output$p_phi_pk  <- renderPlot(line1(custom(), "PHI_conc", "HIF-PHI (au)"))
  output$p_epo_sig <- renderPlot({
    custom() %>%
      select(time, `Endogenous EPO` = EPO_endogenous,
             `Total signal (EPO + ESA equiv.)` = EPO_total_equiv) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_y_log10() + labs(x = "Day", y = "mIU/mL (log scale)") + base_theme
  })

  ## ----------------------------------------------------------------- 3 PD
  output$p_hif   <- renderPlot(line1(custom(), "HIF_activity",
                                     "HIF-2a activity (rel.)", hl = 1))
  output$p_rep   <- renderPlot({
    custom() %>% select(time, `REP pool (x100)` = REP_pool,
                        `EPO (mIU/mL)` = EPO_endogenous) %>%
      mutate(`REP pool (x100)` = 100 * `REP pool (x100)`) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "Day", y = NULL) + base_theme
  })
  output$p_marrow <- renderPlot({
    custom() %>% select(time, PROG, EB1, EB2) %>% pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "Day", y = "10^9 /L") + base_theme
  })
  output$p_retic <- renderPlot({
    custom() %>% select(time, `Retic %` = Reticulocyte_pct,
                        `Retic absolute (10^9/L)` = Reticulocyte_abs) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#8e44ad") +
      facet_wrap(~name, scales = "free_y") + labs(x = "Day", y = NULL) + base_theme
  })
  output$p_hb_pd <- renderPlot({
    custom() %>% select(time, `Hb (g/dL)` = Hemoglobin,
                        `EPOR attenuation (fold)` = EPOR_attenuation) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#c0392b") +
      facet_wrap(~name, scales = "free_y") + labs(x = "Day", y = NULL) + base_theme
  })

  ## --------------------------------------------------------------- 4 iron
  output$p_tsat <- renderPlot({
    custom() %>% select(time, `TSAT (%)` = TSAT_pct,
                        `Iron gate (0-1)` = Iron_gate) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#7f8c33") +
      facet_wrap(~name, scales = "free_y") + labs(x = "Day", y = NULL) + base_theme
  })
  output$p_ferr <- renderPlot({
    custom() %>% select(time, `Ferritin (ng/mL)` = Ferritin,
                        `Storage iron (mg)` = Storage_iron_mg) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#b9770e") +
      facet_wrap(~name, scales = "free_y") + labs(x = "Day", y = NULL) + base_theme
  })
  output$p_hepc <- renderPlot({
    custom() %>% select(time, `Hepcidin-25 (ng/mL)` = Hepcidin,
                        `Erythroferrone (rel.)` = Erythroferrone) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#148f77") +
      facet_wrap(~name, scales = "free_y") + labs(x = "Day", y = NULL) + base_theme
  })
  output$p_tibc <- renderPlot(line1(custom(), "TIBC_level", "TIBC (ug/dL)"))
  output$p_iron_map <- renderPlot({
    compared() %>%
      ggplot(aes(TSAT_pct, Ferritin, colour = scenario)) +
      geom_path(linewidth = 0.8, alpha = 0.9) +
      geom_point(data = ~ group_by(.x, scenario) %>% slice_tail(n = 1), size = 2.5) +
      geom_vline(xintercept = 20, linetype = 2, colour = "grey40") +
      geom_hline(yintercept = 100, linetype = 2, colour = "grey40") +
      annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, size = 3.4,
               colour = "grey30", label = "KDIGO iron-sufficient corner") +
      labs(x = "TSAT (%)", y = "Ferritin (ng/mL)") + base_theme +
      guides(colour = guide_legend(ncol = 1))
  })

  ## ---------------------------------------------------------- 5 endpoints
  output$p_hb_band <- renderPlot({
    ggplot(custom(), aes(time, Hemoglobin)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 10, ymax = 11.5,
               fill = "#2ecc71", alpha = 0.16) +
      geom_hline(yintercept = 13, linetype = 2, colour = "#c0392b") +
      geom_line(linewidth = 1, colour = "#c0392b") +
      annotate("text", x = 0, y = 10.75, hjust = -0.05, size = 3.6,
               colour = "#1d7a45", label = "KDIGO target 10-11.5 g/dL") +
      labs(x = "Day", y = "Hb (g/dL)") + base_theme
  })
  output$p_facit <- renderPlot(line1(custom(), "FACIT_fatigue", "FACIT-Fatigue"))
  output$p_lvmi  <- renderPlot(line1(custom(), "LVMI_g_m2", "LVMI (g/m2)"))
  output$p_egfr  <- renderPlot(line1(custom(), "eGFR_ml_min", "eGFR (mL/min/1.73m2)"))
  output$t_endpoint <- renderDT({
    d <- custom(); n <- nrow(d)
    tibble(
      Endpoint = c("Hb start (g/dL)", "Hb end (g/dL)", "Hb change (g/dL)",
                   "Days in 10-11.5 band", "Days Hb > 13 g/dL",
                   "TSAT end (%)", "Ferritin end (ng/mL)", "Hepcidin end (ng/mL)",
                   "FACIT-Fatigue change", "LVMI change (g/m2)",
                   "Cumulative thrombotic hazard", "MAP end (mmHg)",
                   "LDL-C change (mg/dL)", "Lowest phosphate (mg/dL)"),
      Value = c(round(d$Hemoglobin[1], 2), round(d$Hemoglobin[n], 2),
                round(d$Hemoglobin[n] - d$Hemoglobin[1], 2),
                sum(d$In_target_10_115), sum(d$Hb_overshoot_13),
                round(d$TSAT_pct[n], 1), round(d$Ferritin[n], 0),
                round(d$Hepcidin[n], 1),
                round(d$FACIT_fatigue[n] - d$FACIT_fatigue[1], 1),
                round(d$LVMI_g_m2[n] - d$LVMI_g_m2[1], 1),
                round(d$Thrombotic_hazard[n], 2), round(d$MAP_mmHg[n], 1),
                round(d$LDL_C[n] - d$LDL_C[1], 1), round(min(d$Phosphate), 2))
    ) %>% datatable(rownames = FALSE, options = list(pageLength = 14, dom = "t"))
  })

  ## ------------------------------------------------------------ 6 compare
  output$c_hb   <- renderPlot(lineN(compared(), "Hemoglobin", "Hb (g/dL)",
                                    hl = c(10, 11.5, 13)))
  output$c_hepc <- renderPlot(lineN(compared(), "Hepcidin", "Hepcidin (ng/mL)"))
  output$c_tsat <- renderPlot(lineN(compared(), "TSAT_pct", "TSAT (%)", hl = 20))
  output$c_thr  <- renderPlot(lineN(compared(), "Thrombotic_hazard",
                                    "Cumulative hazard (event-equiv.)"))
  output$c_table <- renderDT({
    compared() %>% group_by(scenario) %>%
      summarise(`Hb start` = round(first(Hemoglobin), 2),
                `Hb end`   = round(last(Hemoglobin), 2),
                `dHb`      = round(last(Hemoglobin) - first(Hemoglobin), 2),
                `Days in band` = sum(In_target_10_115),
                `Days >13`     = sum(Hb_overshoot_13),
                `TSAT end`     = round(last(TSAT_pct), 1),
                `Ferritin end` = round(last(Ferritin), 0),
                `Hepcidin end` = round(last(Hepcidin), 1),
                `EPO peak`     = round(max(EPO_total_equiv), 0),
                `LDL end`      = round(last(LDL_C), 1),
                `MAP end`      = round(last(MAP_mmHg), 1),
                `Min phosphate`= round(min(Phosphate), 2),
                `Cum. hazard`  = round(last(Thrombotic_hazard), 2),
                `dFACIT`       = round(last(FACIT_fatigue) - first(FACIT_fatigue), 1),
                `dLVMI`        = round(last(LVMI_g_m2) - first(LVMI_g_m2), 1),
                .groups = "drop") %>%
      datatable(rownames = FALSE,
                options = list(pageLength = 15, dom = "tp", scrollX = TRUE))
  })

  ## ------------------------------------------------------------- 7 safety
  output$s_over <- renderPlot(lineN(safety_runs(), "Hemoglobin", "Hb (g/dL)", hl = 13))
  output$s_map  <- renderPlot(lineN(safety_runs(), "MAP_mmHg", "MAP (mmHg)"))
  output$s_thr  <- renderPlot(lineN(safety_runs(), "Thrombotic_hazard",
                                    "Cumulative hazard"))
  output$s_ntbi <- renderPlot(lineN(safety_runs(), "NTBI_level", "NTBI (rel.)"))
  output$s_fgf  <- renderPlot(lineN(safety_runs(), "FGF23_level", "FGF23 (rel.)"))
  output$s_phos <- renderPlot(lineN(safety_runs(), "Phosphate",
                                    "Phosphate (mg/dL)", hl = 2.0))

  ## --------------------------------------------------------------- 8 docs
  ti <- eventReactive(input$ti_go, {
    titrate(dose0 = input$ti_dose, lo = input$ti_lo, hi = input$ti_hi,
            blocks = max(1, floor(input$end_d / input$ti_block)),
            block_d = input$ti_block,
            p = c(patient(), ESA_AGENTS[[input$esa_agent]]))
  })
  output$ti_hb <- renderPlot({
    ggplot(ti(), aes(time, Hemoglobin)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = input$ti_lo,
               ymax = input$ti_hi, fill = "#2ecc71", alpha = 0.16) +
      geom_line(linewidth = 1, colour = "#c0392b") +
      labs(x = "Day", y = "Hb (g/dL)") + base_theme
  })
  output$ti_dose_plot <- renderPlot({
    ti() %>% group_by(block) %>% summarise(dose = first(dose_au), .groups = "drop") %>%
      ggplot(aes(block, dose)) + geom_step(linewidth = 0.9, colour = "#2c6fb5") +
      labs(x = "Titration block", y = "ESA dose (au)") + base_theme
  })
}

shinyApp(ui, server)
