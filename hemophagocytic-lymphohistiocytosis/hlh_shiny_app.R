## =============================================================================
## HLH QSP Shiny Dashboard
## Hemophagocytic Lymphohistiocytosis — Interactive Simulation
## =============================================================================
## Tabs: 1-Patient Profile | 2-Drug PK | 3-Cytokine Storm | 4-Clinical Endpoints
##       5-Scenario Comparison | 6-Biomarkers & HScore
## =============================================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)

## --------------------------------------------------------------------------
## Simulation engine (simplified ODE via Euler method for Shiny)
## --------------------------------------------------------------------------
run_hlh_sim <- function(
    tmax = 60, dt = 0.5,
    fNK_dysfunc = 0.3, kTrig = 0.8,
    kMacro_act = 0.5, EC50_IFNg_M = 200,
    # Drug regimens (TRUE/FALSE)
    use_dex = FALSE, use_etop = FALSE, use_csa = FALSE,
    use_emapa = FALSE, use_anakr = FALSE, use_ruxo = FALSE,
    # Drug doses
    dose_dex = 17,  # mg/day
    dose_etop = 270, # mg every 2 weeks
    dose_csa = 360,  # mg/day
    dose_emapa = 70, # mg every 2 weeks
    dose_anakr = 200, # mg/day SC
    dose_ruxo = 20,  # mg BID
    # Initial disease state
    init_IFNg = 50, init_IL6 = 20, init_FERR = 500,
    init_APC = 0.5, init_T_ACT = 0.05, init_MAC = 0.2
) {
    times <- seq(0, tmax, by = dt)
    n <- length(times)

    # State variables
    NK <- rep(1, n); CTL <- rep(1, n)
    T_ACT <- rep(0.01, n); APC <- rep(0.1, n)
    IFNg <- rep(5, n); IL6 <- rep(2, n); TNFa <- rep(1, n)
    IL18 <- rep(10, n); IL10 <- rep(3, n); IL12 <- rep(1, n)
    MAC_ACT <- rep(0.1, n); HEMOPHAG <- rep(0, n)
    FERR <- rep(150, n); BM_SUPP <- rep(0, n)
    LIVER <- rep(0, n); COAG <- rep(0, n)
    # Drug concentrations (simplified 1-cmt)
    C_DEX <- rep(0, n); C_ETOP <- rep(0, n); C_CSA <- rep(0, n)
    C_EMAPA <- rep(0, n); C_ANK <- rep(0, n); C_RUX <- rep(0, n)

    # Set initial disease state
    T_ACT[1] <- init_T_ACT; APC[1] <- init_APC; MAC_ACT[1] <- init_MAC
    IFNg[1] <- init_IFNg; IL6[1] <- init_IL6; FERR[1] <- init_FERR

    # Drug PK parameters (simplified)
    CL_DEX <- 18; V_DEX <- 50; ka_DEX <- 1.2     # DEX ng/mL
    CL_ETOP <- 20; V_ETOP <- 8                     # μg/mL
    CL_CSA <- 25; V_CSA <- 200; ka_CSA <- 0.5      # ng/mL
    CL_EMAPA <- 0.004 * 24; V_EMAPA <- 4           # μg/mL
    CL_ANK <- 2.5 * 24; V_ANK <- 7; ka_ANK <- 0.08 * 24
    CL_RUX <- 17.5 * 24; V_RUX <- 75; ka_RUX <- 1.8 * 24

    depot_DEX <- rep(0, n); depot_CSA <- rep(0, n)
    depot_ANK <- rep(0, n); depot_RUX <- rep(0, n)

    for (i in 2:n) {
        t <- times[i]
        dt_i <- dt

        # Drug dosing events (simplified)
        if (use_dex) {
            dose_interval <- t / 1  # daily
            if (abs(t - round(t)) < dt_i / 2) depot_DEX[i-1] <- depot_DEX[i-1] + dose_dex
        }
        if (use_csa) {
            if (abs(t - round(t)) < dt_i / 2) depot_CSA[i-1] <- depot_CSA[i-1] + dose_csa
        }
        if (use_etop && (abs(t %% 14) < dt_i / 2 || abs(t %% 14 - 14) < dt_i / 2)) {
            C_ETOP[i-1] <- C_ETOP[i-1] + dose_etop / V_ETOP
        }
        if (use_emapa && (abs(t %% 14) < dt_i / 2 || abs(t %% 14 - 14) < dt_i / 2)) {
            C_EMAPA[i-1] <- C_EMAPA[i-1] + dose_emapa / V_EMAPA
        }
        if (use_anakr) {
            if (abs(t - round(t)) < dt_i / 2) depot_ANK[i-1] <- depot_ANK[i-1] + dose_anakr
        }
        if (use_ruxo) {
            if (abs(t %% 0.5) < dt_i / 2) depot_RUX[i-1] <- depot_RUX[i-1] + dose_ruxo
        }

        # Drug PK (simplified Euler)
        depot_DEX[i] <- max(0, depot_DEX[i-1] - ka_DEX * dt_i * depot_DEX[i-1] * 24)
        C_DEX[i] <- max(0, C_DEX[i-1] + dt_i * (ka_DEX * 24 * depot_DEX[i-1] / V_DEX -
                            CL_DEX * 24 / V_DEX * C_DEX[i-1]))

        depot_CSA[i] <- max(0, depot_CSA[i-1] - ka_CSA * dt_i * 24 * depot_CSA[i-1])
        C_CSA[i] <- max(0, C_CSA[i-1] + dt_i * (ka_CSA * 24 * depot_CSA[i-1] / V_CSA -
                            CL_CSA * 24 / V_CSA * C_CSA[i-1]))

        C_ETOP[i] <- max(0, C_ETOP[i-1] - dt_i * CL_ETOP * 24 / V_ETOP * C_ETOP[i-1])

        C_EMAPA[i] <- max(0, C_EMAPA[i-1] - dt_i * CL_EMAPA * 24 / V_EMAPA * C_EMAPA[i-1])

        depot_ANK[i] <- max(0, depot_ANK[i-1] - ka_ANK * dt_i * depot_ANK[i-1])
        C_ANK[i] <- max(0, C_ANK[i-1] + dt_i * (ka_ANK * depot_ANK[i-1] * 0.95 / V_ANK -
                            CL_ANK / V_ANK * C_ANK[i-1]))

        depot_RUX[i] <- max(0, depot_RUX[i-1] - ka_RUX * dt_i * depot_RUX[i-1])
        C_RUX[i] <- max(0, C_RUX[i-1] + dt_i * (ka_RUX * depot_RUX[i-1] / V_RUX -
                            CL_RUX / V_RUX * C_RUX[i-1]))

        # Drug effects (inhibition)
        inh_DEX_IFNg  <- 0.80 * C_DEX[i] / (C_DEX[i] + 50)
        inh_DEX_Macro <- 0.75 * C_DEX[i] / (C_DEX[i] + 60)
        inh_ETOP      <- 0.90 * C_ETOP[i]^2 / (C_ETOP[i]^2 + 4)
        inh_CSA       <- 0.85 * C_CSA[i] / (C_CSA[i] + 300)
        inh_EMAPA     <- 0.95 * (C_EMAPA[i] * 1000) / (C_EMAPA[i] * 1000 + 50)
        inh_ANK       <- 0.80 * C_ANK[i] / (C_ANK[i] + 80)
        inh_RUX_n     <- 0.75 * C_RUX[i]^1.5 / (C_RUX[i]^1.5 + 200^1.5)

        inh_IFNg_total <- 1 - (1 - inh_DEX_IFNg) * (1 - inh_CSA) *
                              (1 - inh_EMAPA) * (1 - inh_RUX_n)

        # Disease ODEs (Euler)
        IL12_stim <- IL12[i-1] / (IL12[i-1] + 10)
        IL18_stim <- IL18[i-1] / (IL18[i-1] + 20)
        TNFa_stim <- TNFa[i-1] / (TNFa[i-1] + 50)
        IFNg_stim_M <- IFNg[i-1] / (IFNg[i-1] + EC50_IFNg_M)

        NK_cytotox <- (1 - fNK_dysfunc) * NK[i-1] * (1 - inh_ETOP)
        CTL_cytotox <- CTL[i-1] * (1 - inh_ETOP)
        APC_clear <- (NK_cytotox + CTL_cytotox) * 0.05

        dNK  <- 1.0 * (1 - fNK_dysfunc) - 0.1 * NK[i-1] - inh_ETOP * 0.05 * NK[i-1]
        dCTL <- 1.0 - 0.1 * CTL[i-1] - inh_ETOP * 0.05 * CTL[i-1]
        dAPC <- kTrig * 2 - 0.3 * APC[i-1] - APC_clear
        dT   <- 0.5 * APC[i-1] * (1 + IL12_stim + IL18_stim) -
                0.4 * T_ACT[i-1] * (1 - inh_ETOP) * (1 - inh_CSA) - 0.2 * T_ACT[i-1]

        dIFNg <- 100 * T_ACT[i-1] * (1 + IL18_stim) * (1 - inh_IFNg_total) + 5 -
                 2.0 * IFNg[i-1]
        dIL6  <- 50 * (APC[i-1] + MAC_ACT[i-1]) * (1 + TNFa_stim) * (1 - inh_RUX_n * 0.6) +
                 2.0 - 4.0 * IL6[i-1]
        dTNFa <- 80 * MAC_ACT[i-1] * (1 + IFNg_stim_M) * (1 - inh_DEX_Macro) + 1 -
                 6.0 * TNFa[i-1]
        dIL18 <- 30 * MAC_ACT[i-1] * (1 + IFNg_stim_M) + 10 - 1.5 * IL18[i-1]
        dIL10 <- 40 * T_ACT[i-1] * MAC_ACT[i-1] + 3 - 3.0 * IL10[i-1]
        dIL12 <- 20 * APC[i-1] * (1 - inh_DEX_Macro) + 1 - 2.0 * IL12[i-1]

        dMAC  <- kMacro_act * IFNg_stim_M * (1 - inh_DEX_Macro) - 0.2 * MAC_ACT[i-1]
        dHEMO <- 0.3 * MAC_ACT[i-1] * (1 - HEMOPHAG[i-1]) - 0.1 * HEMOPHAG[i-1] * (1 - inh_ETOP)
        Ferr_stim <- 1 + IL6[i-1] / (IL6[i-1] + 100) + TNFa[i-1] / (TNFa[i-1] + 100)
        dFERR <- 10 * MAC_ACT[i-1] * Ferr_stim - 0.05 * FERR[i-1]
        dBM   <- 0.4 * HEMOPHAG[i-1] * (1 - BM_SUPP[i-1]) - 0.1 * BM_SUPP[i-1]
        dLIVER <- 0.2 * MAC_ACT[i-1] * (TNFa[i-1] / 100) * (1 - LIVER[i-1]) -
                  0.15 * LIVER[i-1] * (1 - HEMOPHAG[i-1])
        dCOAG <- 0.15 * MAC_ACT[i-1] * HEMOPHAG[i-1] * (1 - COAG[i-1]) - 0.1 * COAG[i-1]

        NK[i]       <- max(0, NK[i-1] + dNK * dt_i)
        CTL[i]      <- max(0, CTL[i-1] + dCTL * dt_i)
        APC[i]      <- max(0, APC[i-1] + dAPC * dt_i)
        T_ACT[i]    <- max(0, T_ACT[i-1] + dT * dt_i)
        IFNg[i]     <- max(0, IFNg[i-1] + dIFNg * dt_i)
        IL6[i]      <- max(0, IL6[i-1] + dIL6 * dt_i)
        TNFa[i]     <- max(0, TNFa[i-1] + dTNFa * dt_i)
        IL18[i]     <- max(0, IL18[i-1] + dIL18 * dt_i)
        IL10[i]     <- max(0, IL10[i-1] + dIL10 * dt_i)
        IL12[i]     <- max(0, IL12[i-1] + dIL12 * dt_i)
        MAC_ACT[i]  <- max(0, MAC_ACT[i-1] + dMAC * dt_i)
        HEMOPHAG[i] <- pmin(1, pmax(0, HEMOPHAG[i-1] + dHEMO * dt_i))
        FERR[i]     <- max(0, FERR[i-1] + dFERR * dt_i)
        BM_SUPP[i]  <- pmin(1, pmax(0, BM_SUPP[i-1] + dBM * dt_i))
        LIVER[i]    <- pmin(1, pmax(0, LIVER[i-1] + dLIVER * dt_i))
        COAG[i]     <- pmin(1, pmax(0, COAG[i-1] + dCOAG * dt_i))
    }

    death_risk <- LIVER * 0.3 + COAG * 0.25 + BM_SUPP * 0.25 + HEMOPHAG * 0.2
    surv_prob  <- exp(-death_risk * 2.5)
    sCD25 <- 100 + T_ACT * 500
    Fibrinogen <- 3.5 * (1 - COAG * 0.8)
    Triglycerides_mM <- 1.5 + MAC_ACT * 4.0
    HScore_approx <- pmin(337, (IFNg / 2000 + IL6 / 1000 + FERR / 100000 +
                                HEMOPHAG + BM_SUPP) * 20)

    data.frame(
        time = times,
        NK_cells = NK, CTL_cells = CTL, T_ACT = T_ACT, APC = APC,
        IFNg = IFNg, IL6 = IL6, TNFa = TNFa, IL18 = IL18, IL10 = IL10, IL12 = IL12,
        MAC_ACT = MAC_ACT, HEMOPHAG = HEMOPHAG, FERR = FERR,
        BM_SUPP = BM_SUPP, LIVER = LIVER, COAG = COAG,
        C_DEX = C_DEX, C_ETOP = C_ETOP, C_CSA = C_CSA,
        C_EMAPA = C_EMAPA, C_ANK = C_ANK, C_RUX = C_RUX,
        Survival_prob = surv_prob, sCD25 = sCD25,
        Fibrinogen = Fibrinogen, Triglycerides_mM = Triglycerides_mM,
        HScore_approx = HScore_approx
    )
}

## --------------------------------------------------------------------------
## UI
## --------------------------------------------------------------------------
ui <- dashboardPage(
    skin = "red",
    dashboardHeader(
        title = span("HLH QSP Dashboard", style = "font-size:16px; font-weight:bold;")
    ),
    dashboardSidebar(
        width = 260,
        sidebarMenu(
            menuItem("Patient Profile",    tabName = "tab_patient",   icon = icon("user")),
            menuItem("Drug PK",            tabName = "tab_pk",        icon = icon("pills")),
            menuItem("Cytokine Storm",     tabName = "tab_cytokine",  icon = icon("fire")),
            menuItem("Clinical Endpoints", tabName = "tab_endpoints", icon = icon("chart-line")),
            menuItem("Scenario Comparison",tabName = "tab_compare",   icon = icon("code-compare")),
            menuItem("Biomarkers & HScore",tabName = "tab_biomarkers",icon = icon("microscope"))
        ),
        tags$hr(),
        tags$div(style = "padding: 10px; font-size: 12px; color: #aaa;",
                 "Simulation Settings"),

        selectInput("hlh_type", "HLH Subtype",
                    choices = c("Secondary HLH (viral)" = "secondary",
                                "Primary/Familial HLH" = "primary",
                                "MAS / sJIA-HLH" = "mas"),
                    selected = "secondary"),

        sliderInput("sim_days", "Simulation Duration (days)", 7, 90, 60),

        tags$div(style = "padding: 10px; font-size: 12px; color: #aaa;",
                 "Disease Parameters"),
        sliderInput("nk_dysfunc", "NK Dysfunction Fraction", 0, 1, 0.3, 0.05),
        sliderInput("trig_strength", "Trigger Strength", 0, 1, 0.8, 0.1),

        tags$div(style = "padding: 10px; font-size: 12px; color: #aaa;",
                 "Initial Conditions"),
        numericInput("init_IFNg", "Initial IFN-γ (pg/mL)", 50, 5, 500),
        numericInput("init_FERR", "Initial Ferritin (ng/mL)", 500, 100, 50000),

        tags$div(style = "padding: 10px; font-size: 12px; color: #aaa;",
                 "Treatment Selection"),
        checkboxInput("use_dex",   "Dexamethasone", TRUE),
        checkboxInput("use_etop",  "Etoposide", FALSE),
        checkboxInput("use_csa",   "Cyclosporine A", FALSE),
        checkboxInput("use_emapa", "Emapalumab (anti-IFN-γ)", FALSE),
        checkboxInput("use_anakr", "Anakinra (IL-1Ra)", FALSE),
        checkboxInput("use_ruxo",  "Ruxolitinib", FALSE),

        tags$div(style = "padding: 10px; font-size: 12px; color: #aaa;",
                 "Drug Doses"),
        numericInput("dose_dex",   "DEX dose (mg/day)", 17, 1, 40),
        numericInput("dose_etop",  "ETOP dose (mg Q2W)", 270, 50, 500),
        numericInput("dose_csa",   "CsA dose (mg/day)", 360, 100, 800),
        numericInput("dose_emapa", "Emapalumab (mg Q2W)", 70, 10, 200),
        numericInput("dose_anakr", "Anakinra (mg/day)", 200, 50, 600),
        numericInput("dose_ruxo",  "Ruxolitinib (mg BID)", 20, 5, 50)
    ),
    dashboardBody(
        tabItems(
            # ------------------------------------------------------------------
            # TAB 1: Patient Profile
            # ------------------------------------------------------------------
            tabItem("tab_patient",
                h2("HLH Patient Profile & Disease Overview"),
                fluidRow(
                    valueBoxOutput("vbox_ifng"),
                    valueBoxOutput("vbox_ferr"),
                    valueBoxOutput("vbox_surv")
                ),
                fluidRow(
                    box(title = "Disease Progression Overview", width = 8, status = "danger",
                        plotlyOutput("plot_overview", height = 350)),
                    box(title = "HLH Subtype Characteristics", width = 4, status = "info",
                        tableOutput("tbl_subtype"))
                ),
                fluidRow(
                    box(title = "What is HLH?", width = 12, status = "primary",
                        solidHeader = FALSE,
                        p("Hemophagocytic Lymphohistiocytosis (HLH) is a life-threatening
                          hyperinflammatory syndrome characterized by uncontrolled immune
                          activation leading to cytokine storm, macrophage activation,
                          and multi-organ damage. Without treatment, mortality exceeds 50%
                          within 5 weeks."),
                        p(strong("Primary HLH"), ": Genetic defects in NK/CTL cytotoxicity
                          (perforin, MUNC13-4, Rab27a, etc.) — requires HSCT for cure."),
                        p(strong("Secondary HLH / MAS"), ": Triggered by infections (EBV, CMV),
                          autoimmune diseases (sJIA, SLE), or malignancy — more often responsive
                          to immunosuppression alone."),
                        p(strong("Diagnostic Criteria"), ": HLH-2004 (≥5 of 8 criteria) or
                          HScore ≥169 (93% sensitivity). Key biomarkers: ferritin >10,000 ng/mL,
                          sCD25 >2,400 U/mL, NK activity ↓.")
                    )
                )
            ),

            # ------------------------------------------------------------------
            # TAB 2: Drug PK
            # ------------------------------------------------------------------
            tabItem("tab_pk",
                h2("Drug Pharmacokinetics"),
                fluidRow(
                    box(title = "Drug Concentration-Time Profiles", width = 8, status = "primary",
                        plotlyOutput("plot_pk", height = 400)),
                    box(title = "PK Parameters Summary", width = 4, status = "info",
                        tableOutput("tbl_pk_params"))
                ),
                fluidRow(
                    box(title = "Drug Mechanism of Action", width = 12, status = "warning",
                        tableOutput("tbl_moa"))
                )
            ),

            # ------------------------------------------------------------------
            # TAB 3: Cytokine Storm
            # ------------------------------------------------------------------
            tabItem("tab_cytokine",
                h2("Cytokine Storm Dynamics"),
                fluidRow(
                    box(title = "Pro-inflammatory Cytokines", width = 6, status = "danger",
                        plotlyOutput("plot_cytokine1", height = 350)),
                    box(title = "Macrophage Activation & Hemophagocytosis", width = 6,
                        status = "warning",
                        plotlyOutput("plot_cytokine2", height = 350))
                ),
                fluidRow(
                    box(title = "Cytokine Network — All Mediators", width = 12, status = "info",
                        plotlyOutput("plot_cytokine_all", height = 350))
                )
            ),

            # ------------------------------------------------------------------
            # TAB 4: Clinical Endpoints
            # ------------------------------------------------------------------
            tabItem("tab_endpoints",
                h2("Clinical Endpoints & Outcomes"),
                fluidRow(
                    box(title = "Survival Probability", width = 6, status = "danger",
                        plotlyOutput("plot_survival", height = 300)),
                    box(title = "Organ Damage Indices", width = 6, status = "warning",
                        plotlyOutput("plot_organ_damage", height = 300))
                ),
                fluidRow(
                    box(title = "Bone Marrow Suppression & Cytopenias", width = 6,
                        status = "info",
                        plotlyOutput("plot_bm", height = 300)),
                    box(title = "Key Endpoint Summary (Day 14 & 30)", width = 6,
                        status = "success",
                        tableOutput("tbl_endpoints_summary"))
                )
            ),

            # ------------------------------------------------------------------
            # TAB 5: Scenario Comparison
            # ------------------------------------------------------------------
            tabItem("tab_compare",
                h2("Multi-Scenario Treatment Comparison"),
                fluidRow(
                    box(title = "IFN-γ: All Treatment Scenarios", width = 6, status = "danger",
                        plotlyOutput("plot_comp_ifng", height = 320)),
                    box(title = "Ferritin: All Treatment Scenarios", width = 6,
                        status = "warning",
                        plotlyOutput("plot_comp_ferr", height = 320))
                ),
                fluidRow(
                    box(title = "Survival Probability Comparison", width = 6,
                        status = "primary",
                        plotlyOutput("plot_comp_surv", height = 320)),
                    box(title = "Hemophagocytosis Index", width = 6, status = "info",
                        plotlyOutput("plot_comp_hemo", height = 320))
                ),
                fluidRow(
                    box(title = "HScore Trajectory", width = 12, status = "success",
                        plotlyOutput("plot_comp_hscore", height = 280))
                )
            ),

            # ------------------------------------------------------------------
            # TAB 6: Biomarkers & HScore
            # ------------------------------------------------------------------
            tabItem("tab_biomarkers",
                h2("Biomarkers & HScore Monitoring"),
                fluidRow(
                    box(title = "HScore Trajectory", width = 6, status = "danger",
                        plotlyOutput("plot_hscore", height = 300)),
                    box(title = "Serum Ferritin (Log Scale)", width = 6, status = "warning",
                        plotlyOutput("plot_ferr_log", height = 300))
                ),
                fluidRow(
                    box(title = "sCD25 & NK Activity", width = 6, status = "info",
                        plotlyOutput("plot_scd25", height = 300)),
                    box(title = "Fibrinogen & Triglycerides", width = 6, status = "primary",
                        plotlyOutput("plot_coag_lipids", height = 300))
                ),
                fluidRow(
                    box(title = "HLH-2004 Diagnostic Criteria Tracker", width = 12,
                        status = "success",
                        tableOutput("tbl_hlh2004_criteria"))
                )
            )
        )
    )
)

## --------------------------------------------------------------------------
## SERVER
## --------------------------------------------------------------------------
server <- function(input, output, session) {

    ## Reactive: HLH type defaults
    hlh_defaults <- reactive({
        switch(input$hlh_type,
               "primary"   = list(nk = 0.85, trig = 0.9, init_IFNg = 150,
                                   init_FERR = 3000, init_MAC = 0.4),
               "secondary" = list(nk = 0.3,  trig = 0.8, init_IFNg = 50,
                                   init_FERR = 500,  init_MAC = 0.2),
               "mas"       = list(nk = 0.1,  trig = 0.5, init_IFNg = 30,
                                   init_FERR = 1500, init_MAC = 0.3))
    })

    ## Update sliders when HLH type changes
    observe({
        d <- hlh_defaults()
        updateSliderInput(session, "nk_dysfunc", value = d$nk)
        updateSliderInput(session, "trig_strength", value = d$trig)
        updateNumericInput(session, "init_IFNg", value = d$init_IFNg)
        updateNumericInput(session, "init_FERR", value = d$init_FERR)
    })

    ## Main simulation (current settings)
    sim_data <- reactive({
        run_hlh_sim(
            tmax = input$sim_days, dt = 0.5,
            fNK_dysfunc = input$nk_dysfunc,
            kTrig = input$trig_strength,
            use_dex = input$use_dex, use_etop = input$use_etop,
            use_csa = input$use_csa, use_emapa = input$use_emapa,
            use_anakr = input$use_anakr, use_ruxo = input$use_ruxo,
            dose_dex = input$dose_dex, dose_etop = input$dose_etop,
            dose_csa = input$dose_csa, dose_emapa = input$dose_emapa,
            dose_anakr = input$dose_anakr, dose_ruxo = input$dose_ruxo,
            init_IFNg = input$init_IFNg, init_IL6 = 20,
            init_FERR = input$init_FERR, init_APC = 0.5,
            init_T_ACT = 0.05, init_MAC = 0.2
        )
    })

    ## All 5 scenarios comparison data
    compare_data <- reactive({
        params_base <- list(tmax = input$sim_days, dt = 0.5,
                            fNK_dysfunc = input$nk_dysfunc,
                            kTrig = input$trig_strength,
                            init_IFNg = input$init_IFNg,
                            init_FERR = input$init_FERR)

        bind_rows(
            do.call(run_hlh_sim, c(params_base,
                list(use_dex=F,use_etop=F,use_csa=F,use_emapa=F,use_anakr=F,use_ruxo=F))) %>%
                mutate(Scenario = "Untreated"),
            do.call(run_hlh_sim, c(params_base,
                list(use_dex=T,use_etop=T,use_csa=T,use_emapa=F,use_anakr=F,use_ruxo=F))) %>%
                mutate(Scenario = "HLH-2004 (DEX+ETOP+CsA)"),
            do.call(run_hlh_sim, c(params_base,
                list(use_dex=T,use_etop=F,use_csa=F,use_emapa=T,use_anakr=F,use_ruxo=F))) %>%
                mutate(Scenario = "Emapalumab + DEX"),
            do.call(run_hlh_sim, c(params_base,
                list(use_dex=T,use_etop=F,use_csa=F,use_emapa=F,use_anakr=T,use_ruxo=F))) %>%
                mutate(Scenario = "Anakinra + DEX (MAS)"),
            do.call(run_hlh_sim, c(params_base,
                list(use_dex=T,use_etop=F,use_csa=F,use_emapa=F,use_anakr=F,use_ruxo=T))) %>%
                mutate(Scenario = "Ruxolitinib + DEX")
        )
    })

    scen_colors <- c(
        "Untreated" = "#E74C3C",
        "HLH-2004 (DEX+ETOP+CsA)" = "#2980B9",
        "Emapalumab + DEX" = "#8E44AD",
        "Anakinra + DEX (MAS)" = "#27AE60",
        "Ruxolitinib + DEX" = "#D35400"
    )

    ## ---- Value Boxes ----
    output$vbox_ifng <- renderValueBox({
        d <- sim_data()
        val <- round(tail(d$IFNg, 1), 1)
        color <- if (val > 500) "red" else if (val > 100) "orange" else "green"
        valueBox(paste0(val, " pg/mL"), "Final IFN-γ", icon = icon("fire"), color = color)
    })
    output$vbox_ferr <- renderValueBox({
        d <- sim_data()
        val <- round(tail(d$FERR, 1), 0)
        color <- if (val > 10000) "red" else if (val > 500) "orange" else "green"
        valueBox(paste0(format(val, big.mark=","), " ng/mL"), "Final Ferritin",
                 icon = icon("vials"), color = color)
    })
    output$vbox_surv <- renderValueBox({
        d <- sim_data()
        val <- round(tail(d$Survival_prob, 1) * 100, 1)
        color <- if (val < 40) "red" else if (val < 75) "orange" else "green"
        valueBox(paste0(val, "%"), "Survival Probability", icon = icon("heart-pulse"),
                 color = color)
    })

    ## ---- Overview Plot ----
    output$plot_overview <- renderPlotly({
        d <- sim_data()
        p <- ggplot(d, aes(time)) +
            geom_line(aes(y = IFNg / max(IFNg) * 100, color = "IFN-γ (scaled)"), linewidth=1) +
            geom_line(aes(y = FERR / max(FERR) * 100, color = "Ferritin (scaled)"), linewidth=1) +
            geom_line(aes(y = Survival_prob * 100, color = "Survival %"), linewidth=1.5,
                      linetype = "dashed") +
            geom_line(aes(y = HEMOPHAG * 100, color = "Hemophagocytosis Index (%)"),
                      linewidth=1) +
            scale_color_manual(values = c("IFN-γ (scaled)"="#E74C3C",
                                           "Ferritin (scaled)"="#E67E22",
                                           "Survival %"="#27AE60",
                                           "Hemophagocytosis Index (%)"="#8E44AD")) +
            labs(x="Time (days)", y="% of peak / %", color="Metric") +
            theme_bw()
        ggplotly(p)
    })

    ## ---- Subtype Table ----
    output$tbl_subtype <- renderTable({
        data.frame(
            Feature = c("Genetic basis", "NK dysfunction", "Trigger", "Primary treatment",
                         "HSCT needed", "Mortality (untreated)"),
            Primary_HLH = c("Yes (PRF1, UNC13D...)", "Severe (>70%)", "Often absent",
                             "DEX + ETOP + CsA", "Yes (curative)", ">90%"),
            Secondary_HLH = c("No", "Moderate (20-60%)", "Infection/malignancy",
                               "Treat trigger + DEX", "Rarely", "~50%"),
            MAS = c("No (predisposing)", "Mild (<30%)", "sJIA/SLE flare",
                    "High-dose corticosteroids ± IL-1/6 blockade", "Rarely", "10-20%")
        )
    })

    ## ---- Drug PK Plot ----
    output$plot_pk <- renderPlotly({
        d <- sim_data()
        pk_long <- d %>%
            select(time, C_DEX, C_ETOP, C_CSA, C_EMAPA, C_ANK, C_RUX) %>%
            pivot_longer(-time, names_to = "Drug", values_to = "Conc") %>%
            filter(Conc > 0.001)
        p <- ggplot(pk_long, aes(time, Conc, color = Drug)) +
            geom_line(linewidth = 1) +
            labs(x = "Time (days)", y = "Concentration (normalized units)", color = "Drug") +
            theme_bw()
        ggplotly(p)
    })

    output$tbl_pk_params <- renderTable({
        data.frame(
            Drug = c("Dexamethasone", "Etoposide", "Cyclosporine A",
                     "Emapalumab", "Anakinra", "Ruxolitinib"),
            Route = c("Oral/IV", "IV infusion", "Oral", "IV", "SC", "Oral"),
            t_half = c("3-5 h", "4-11 h", "8-12 h", "~22 days", "~6 h", "~3 h"),
            Target = c("GR (NF-κB↓)", "Topo-II (T cell↓)", "Calcineurin (IL-2↓)",
                       "IFN-γ (neutralize)", "IL-1R (block)", "JAK1/2 (STAT↓)")
        )
    })

    output$tbl_moa <- renderTable({
        data.frame(
            Drug = c("Dexamethasone", "Etoposide", "Cyclosporine A",
                     "Emapalumab", "Anakinra", "Ruxolitinib"),
            Primary_Target = c("GR/NF-κB", "Topoisomerase II", "Calcineurin/NFAT",
                                "IFN-γ", "IL-1 receptor", "JAK1/2"),
            Effect_on_HLH = c("↓ IFN-γ, ↓ TNF-α, ↓ macrophage activation",
                               "Depletes activated T/NK cells",
                               "↓ IL-2, IFN-γ; blocks T cell proliferation",
                               "Directly neutralizes free IFN-γ",
                               "Blocks IL-1β signaling (esp. MAS/NLRP3)",
                               "Inhibits JAK-STAT downstream of IFN-γ, IL-6"),
            FDA_Approval = c("Yes (off-label HLH)", "Yes (HLH-2004 protocol)",
                              "Yes (HLH-2004 protocol)", "Yes (primary/refractory HLH, 2018)",
                              "No (compassionate use, MAS)", "No (clinical trials)")
        )
    })

    ## ---- Cytokine Plots ----
    output$plot_cytokine1 <- renderPlotly({
        d <- sim_data()
        d_long <- d %>% select(time, IFNg, TNFa, IL18, IL12) %>%
            pivot_longer(-time)
        p <- ggplot(d_long, aes(time, value, color = name)) +
            geom_line(linewidth = 1) +
            scale_color_manual(values = c(IFNg="#E74C3C", TNFa="#CB4335",
                                           IL18="#E67E22", IL12="#F39C12")) +
            labs(x="Time (days)", y="Concentration (pg/mL)", color="Cytokine") + theme_bw()
        ggplotly(p)
    })
    output$plot_cytokine2 <- renderPlotly({
        d <- sim_data()
        p <- ggplot(d, aes(time)) +
            geom_line(aes(y = MAC_ACT * 100, color = "Macrophage Activation (%)"), linewidth=1) +
            geom_line(aes(y = HEMOPHAG * 100, color = "Hemophagocytosis Index (%)"), linewidth=1) +
            scale_color_manual(values = c("Macrophage Activation (%)"="#D35400",
                                           "Hemophagocytosis Index (%)"="#8E44AD")) +
            labs(x="Time (days)", y="Index (%)", color="Variable") + theme_bw()
        ggplotly(p)
    })
    output$plot_cytokine_all <- renderPlotly({
        d <- sim_data()
        d_long <- d %>% select(time, IFNg, IL6, TNFa, IL18, IL10, IL12) %>%
            pivot_longer(-time)
        p <- ggplot(d_long, aes(time, value, color = name)) +
            geom_line(linewidth = 1) +
            facet_wrap(~name, scales = "free_y") +
            labs(x="Time (days)", y="pg/mL") + theme_bw() + theme(legend.position = "none")
        ggplotly(p)
    })

    ## ---- Clinical Endpoints Plots ----
    output$plot_survival <- renderPlotly({
        d <- sim_data()
        p <- ggplot(d, aes(time, Survival_prob * 100)) +
            geom_line(color = "#27AE60", linewidth = 1.5) +
            geom_ribbon(aes(ymin = pmax(0, Survival_prob*100 - 10),
                             ymax = pmin(100, Survival_prob*100 + 10)),
                        alpha = 0.15, fill = "#27AE60") +
            geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
            labs(x="Time (days)", y="Survival Probability (%)") +
            coord_cartesian(ylim=c(0,100)) + theme_bw()
        ggplotly(p)
    })
    output$plot_organ_damage <- renderPlotly({
        d <- sim_data()
        d_long <- d %>% select(time, LIVER, COAG) %>%
            pivot_longer(-time, names_to = "Organ",
                         values_to = "Damage Index")
        p <- ggplot(d_long, aes(time, `Damage Index`, color = Organ)) +
            geom_line(linewidth = 1.2) +
            scale_color_manual(values = c(LIVER="#E74C3C", COAG="#8E44AD")) +
            labs(x="Time (days)", y="Damage Index (0-1)") + theme_bw()
        ggplotly(p)
    })
    output$plot_bm <- renderPlotly({
        d <- sim_data()
        p <- ggplot(d, aes(time, BM_SUPP * 100)) +
            geom_line(color = "#2980B9", linewidth = 1.2) +
            labs(x="Time (days)", y="BM Suppression Index (%)") + theme_bw()
        ggplotly(p)
    })
    output$tbl_endpoints_summary <- renderTable({
        d <- sim_data()
        d14 <- d %>% filter(abs(time - 14) < 0.3) %>% tail(1)
        d30 <- d %>% filter(abs(time - 30) < 0.3) %>% tail(1)
        data.frame(
            Timepoint   = c("Day 14", "Day 30"),
            IFNg_pgmL   = round(c(d14$IFNg, d30$IFNg), 1),
            Ferritin     = round(c(d14$FERR, d30$FERR), 0),
            HemophagIdx  = round(c(d14$HEMOPHAG, d30$HEMOPHAG), 3),
            Survival_pct = round(c(d14$Survival_prob, d30$Survival_prob) * 100, 1),
            HScore       = round(c(d14$HScore_approx, d30$HScore_approx), 1)
        )
    })

    ## ---- Scenario Comparison Plots ----
    output$plot_comp_ifng <- renderPlotly({
        d <- compare_data()
        p <- ggplot(d, aes(time, IFNg, color = Scenario)) +
            geom_line(linewidth = 1) +
            scale_color_manual(values = scen_colors) +
            labs(x="Time (days)", y="IFN-γ (pg/mL)") + theme_bw() +
            theme(legend.position = "bottom")
        ggplotly(p)
    })
    output$plot_comp_ferr <- renderPlotly({
        d <- compare_data()
        p <- ggplot(d, aes(time, FERR, color = Scenario)) +
            geom_line(linewidth = 1) + scale_y_log10() +
            scale_color_manual(values = scen_colors) +
            labs(x="Time (days)", y="Ferritin (ng/mL, log)") + theme_bw() +
            theme(legend.position = "bottom")
        ggplotly(p)
    })
    output$plot_comp_surv <- renderPlotly({
        d <- compare_data()
        p <- ggplot(d, aes(time, Survival_prob * 100, color = Scenario)) +
            geom_line(linewidth = 1.2) +
            scale_color_manual(values = scen_colors) +
            coord_cartesian(ylim = c(0, 100)) +
            labs(x="Time (days)", y="Survival Probability (%)") + theme_bw() +
            theme(legend.position = "bottom")
        ggplotly(p)
    })
    output$plot_comp_hemo <- renderPlotly({
        d <- compare_data()
        p <- ggplot(d, aes(time, HEMOPHAG, color = Scenario)) +
            geom_line(linewidth = 1) +
            scale_color_manual(values = scen_colors) +
            labs(x="Time (days)", y="Hemophagocytosis Index") + theme_bw() +
            theme(legend.position = "bottom")
        ggplotly(p)
    })
    output$plot_comp_hscore <- renderPlotly({
        d <- compare_data()
        p <- ggplot(d, aes(time, HScore_approx, color = Scenario)) +
            geom_line(linewidth = 1) +
            geom_hline(yintercept = 169, linetype = "dashed", color = "red") +
            scale_color_manual(values = scen_colors) +
            labs(x="Time (days)", y="HScore (0-337)") + theme_bw() +
            theme(legend.position = "bottom")
        ggplotly(p)
    })

    ## ---- Biomarker Plots ----
    output$plot_hscore <- renderPlotly({
        d <- sim_data()
        p <- ggplot(d, aes(time, HScore_approx)) +
            geom_line(color = "#8E44AD", linewidth = 1.5) +
            geom_hline(yintercept = 169, linetype = "dashed", color = "red") +
            annotate("text", x = 5, y = 175, label = "Diagnostic threshold (169)",
                     size = 3, color = "red") +
            labs(x="Time (days)", y="HScore (0-337)") + theme_bw()
        ggplotly(p)
    })
    output$plot_ferr_log <- renderPlotly({
        d <- sim_data()
        p <- ggplot(d, aes(time, FERR)) +
            geom_line(color = "#E67E22", linewidth = 1.5) +
            scale_y_log10() +
            geom_hline(yintercept = c(500, 10000), linetype = c("dashed","solid"),
                       color = c("orange","red"), alpha = 0.7) +
            labs(x="Time (days)", y="Ferritin (ng/mL, log scale)") + theme_bw()
        ggplotly(p)
    })
    output$plot_scd25 <- renderPlotly({
        d <- sim_data()
        p <- ggplot(d, aes(time)) +
            geom_line(aes(y = sCD25, color = "sCD25 (U/mL)"), linewidth = 1) +
            geom_line(aes(y = NK_cells * 100, color = "NK Function (%)"), linewidth = 1) +
            scale_color_manual(values = c("sCD25 (U/mL)"="#2980B9",
                                           "NK Function (%)"="#27AE60")) +
            labs(x="Time (days)", y="Value", color="Biomarker") + theme_bw()
        ggplotly(p)
    })
    output$plot_coag_lipids <- renderPlotly({
        d <- sim_data()
        p <- ggplot(d, aes(time)) +
            geom_line(aes(y = Fibrinogen, color = "Fibrinogen (g/L)"), linewidth = 1) +
            geom_line(aes(y = Triglycerides_mM, color = "Triglycerides (mmol/L)"),
                      linewidth = 1) +
            geom_hline(yintercept = 1.5, linetype = "dashed", color = "red", alpha = 0.5) +
            scale_color_manual(values = c("Fibrinogen (g/L)"="#E74C3C",
                                           "Triglycerides (mmol/L)"="#F39C12")) +
            labs(x="Time (days)", y="Value", color="Parameter") + theme_bw()
        ggplotly(p)
    })
    output$tbl_hlh2004_criteria <- renderTable({
        d <- sim_data()
        d_end <- tail(d, 1)
        data.frame(
            Criterion = c("Fever >38.5°C",
                          "Cytopenias (≥2 cell lines)",
                          "Hemophagocytosis in BM/spleen/LN",
                          "Hypertriglyceridemia ≥3 mmol/L",
                          "Hypofibrinogenemia <1.5 g/L",
                          "Low/absent NK cell activity",
                          "sCD25 ≥2400 U/mL",
                          "Splenomegaly"),
            Status = c(
                ifelse(d_end$MAC_ACT > 0.3, "PRESENT", "absent"),
                ifelse(d_end$BM_SUPP > 0.3, "PRESENT", "absent"),
                ifelse(d_end$HEMOPHAG > 0.3, "PRESENT", "absent"),
                ifelse(d_end$Triglycerides_mM > 3, "PRESENT", "absent"),
                ifelse(d_end$Fibrinogen < 1.5, "PRESENT", "absent"),
                ifelse(d_end$NK_cells * (1 - input$nk_dysfunc) < 0.5, "PRESENT", "absent"),
                ifelse(d_end$sCD25 > 2400, "PRESENT", "absent"),
                ifelse(d_end$HEMOPHAG > 0.4, "PRESENT (extrapolated)", "absent")
            ),
            `Reference Value` = c(">38.5°C for >7 days", "ANC <1000, Hgb <9, Plt <100",
                                    "≥3 phagocytic events per field", "≥3 mmol/L",
                                    "<1.5 g/L", "<Natural killer cell lytic units 10",
                                    "≥2400 U/mL", ">22 cm on ultrasound")
        )
    })
}

shinyApp(ui = ui, server = server)
