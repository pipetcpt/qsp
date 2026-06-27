## ============================================================
## Lambert-Eaton Myasthenic Syndrome (LEMS) QSP — Shiny App
## ============================================================
## Interactive dashboard for LEMS mechanistic PK/PD simulation
## Tabs: Patient Profile · PK · VGCC/Antibody · NMJ/CMAP ·
##       Clinical Endpoints · Scenario Comparison · Biomarkers
## ============================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(shinydashboard)
library(tidyr)

## ---- Embed model code (copy from lems_mrgsolve_model.R) ----
lems_code <- '
$PARAM
Ka_amif=0.693, CL_amif=18.0, V1_amif=35.0, V2_amif=60.0, Q_amif=5.0, F_amif=1.0,
Ka_pred=1.2, CL_pred=1.3, V1_pred=15.0, V2_pred=30.0, Q_pred=2.0, F_pred=0.82,
Emax_kch=0.85, EC50_kch=120.0, nH_kch=1.5,
Ab0=2000.0, kin_Ab=100.0, kout_Ab=0.05, kout_Ab_PE=0.0, inh_pred_Ab=0.7,
k_VGCC_block=0.001, k_VGCC_recov=0.015,
Ca_basal=0.1, kCa_in=5.0, kCa_out=2.0, Ca_thresh=1.5, Hill_Ca=3.0,
V_RRP0=100.0, V_reserve=1000.0, k_refill=0.1, k_deplete=0.02, Vmax_fuse=50.0,
EPP_base=40.0, EPP_thresh=10.0,
CMAP_norm=5.0, CMAP_min=0.3, kCMAP=0.5,
QMG_max=39.0, k_QMG=0.8,
Bcell0=1.0, kB_prolif=0.02, kB_death=0.02, IC50_pred_B=50.0,
Tumor0=0.0, kTumor_growth=0.05, kchemo=0.0, Tumor_max=10.0,
tau_facil=0.5

$CMT A_gut A_central A_periph P_gut P_central P_periph
     Ab_VGCC VGCC_free RRP EPP_amp CMAP QMG Bcell Tumor Facil

$MAIN
double Camif = A_central / V1_amif * 1000.0;
double Cpred = P_central / V1_pred  * 1000.0;
double Kblock = Emax_kch * pow(Camif, nH_kch) /
                (pow(EC50_kch, nH_kch) + pow(Camif, nH_kch));
double VGCC_effective = VGCC_free * (1.0 + 2.5 * Kblock);
if(VGCC_effective > 1.5) VGCC_effective = 1.5;
double Ca_pre = Ca_basal + kCa_in * VGCC_effective / kCa_out;
double Ca_facilitated = Ca_pre * (1.0 + 0.5 * Facil);
double Ca_above = (Ca_facilitated > Ca_basal) ? (Ca_facilitated - Ca_basal) : 0.0;
double F_ACh = Vmax_fuse * pow(Ca_above, Hill_Ca) /
               (pow(Ca_thresh - Ca_basal, Hill_Ca) + pow(Ca_above, Hill_Ca));
double EPP_ss = EPP_base * (F_ACh / Vmax_fuse) * (RRP / V_RRP0);
double EPP_ratio = EPP_ss / EPP_thresh;
double CMAP_frac = (EPP_ratio > 0.2) ?
    (1.0 / (1.0 + exp(-kCMAP * (EPP_ratio - 1.5)))) : 0.0;
double CMAP_ss = CMAP_min + (CMAP_norm - CMAP_min) * CMAP_frac;
double QMG_ss = QMG_max * (1.0 - CMAP_ss / CMAP_norm) * 0.8;
if(QMG_ss < 0) QMG_ss = 0;
if(QMG_ss > QMG_max) QMG_ss = QMG_max;
double Imax_pred = inh_pred_Ab * Cpred / (IC50_pred_B + Cpred);
double tumor_factor = (Tumor > 0) ? (1.0 + 0.5 * Tumor) : 1.0;
if(NEWIND <= 1) {
    A_gut_0=0; A_central_0=0; A_periph_0=0;
    P_gut_0=0; P_central_0=0; P_periph_0=0;
    Ab_VGCC_0=Ab0;
    VGCC_free_0=1.0-(Ab0*k_VGCC_block/(k_VGCC_block*Ab0+k_VGCC_recov));
    RRP_0=V_RRP0; EPP_amp_0=EPP_base;
    CMAP_0=CMAP_norm; QMG_0=0;
    Bcell_0=Bcell0; Tumor_0=Tumor0; Facil_0=0;
}

$ODE
dxdt_A_gut     = -Ka_amif * A_gut;
dxdt_A_central =  Ka_amif * A_gut * F_amif
                - (CL_amif/V1_amif)*A_central - (Q_amif/V1_amif)*A_central
                + (Q_amif/V2_amif)*A_periph;
dxdt_A_periph  =  (Q_amif/V1_amif)*A_central - (Q_amif/V2_amif)*A_periph;
dxdt_P_gut     = -Ka_pred * P_gut;
dxdt_P_central =  Ka_pred*P_gut*F_pred - (CL_pred/V1_pred)*P_central
                - (Q_pred/V1_pred)*P_central + (Q_pred/V2_pred)*P_periph;
dxdt_P_periph  =  (Q_pred/V1_pred)*P_central - (Q_pred/V2_pred)*P_periph;
double kin_eff  = kin_Ab * tumor_factor * Bcell * (1.0 - Imax_pred);
dxdt_Ab_VGCC   = kin_eff - (kout_Ab + kout_Ab_PE) * Ab_VGCC;
dxdt_VGCC_free = k_VGCC_recov*(1.0-VGCC_free) - k_VGCC_block*Ab_VGCC*VGCC_free;
dxdt_RRP       = k_refill*(V_RRP0-RRP) - k_deplete*F_ACh;
dxdt_EPP_amp   = 10.0*(EPP_ss - EPP_amp);
dxdt_CMAP      = 5.0*(CMAP_ss - CMAP);
dxdt_QMG       = 2.0*(QMG_ss - QMG);
dxdt_Bcell     = kB_prolif*Bcell*(1.0-Imax_pred) - kB_death*Bcell;
dxdt_Tumor     = kTumor_growth*Tumor*(1.0-Tumor/Tumor_max) - kchemo*Tumor;
dxdt_Facil     = -Facil/tau_facil;

$TABLE
double Camif_ngmL      = A_central/V1_amif*1000.0;
double Cpred_ngmL      = P_central/V1_pred*1000.0;
double VGCC_blocked_pct= (1.0-VGCC_free)*100.0;
double Ab_fold         = Ab_VGCC/Ab0;
double Safety_Factor   = EPP_amp/EPP_thresh;
double CMAP_pct_normal = CMAP/CMAP_norm*100.0;
double Kblock_frac     = Emax_kch*pow(Camif_ngmL,nH_kch)/
                         (pow(EC50_kch,nH_kch)+pow(Camif_ngmL,nH_kch));

$CAPTURE Camif_ngmL Cpred_ngmL VGCC_blocked_pct Ab_fold Safety_Factor
         CMAP_pct_normal Kblock_frac EPP_amp CMAP QMG Tumor Bcell VGCC_free Ab_VGCC
'

mod <- mcode_cache("LEMS_shiny", lems_code, quiet = TRUE)

## ---- Helper: run simulation ----
run_sim <- function(
    ab0, amif_dose, amif_interval, pred_dose,
    pe_sessions, chemo_effect, sim_days,
    paraneoplastic, tumor_growth
) {
    kout_pe <- ifelse(pe_sessions > 0, pe_sessions * 0.15, 0)

    # Build events
    ev_list <- list()
    if (amif_dose > 0) {
        times_amif <- seq(0, sim_days * 24, by = amif_interval)
        ev_list[[length(ev_list) + 1]] <- ev(
            cmt = 1, amt = amif_dose, time = times_amif)
    }
    if (pred_dose > 0) {
        times_pred <- seq(0, sim_days * 24, by = 24)
        ev_list[[length(ev_list) + 1]] <- ev(
            cmt = 4, amt = pred_dose, time = times_pred)
    }

    ev_use <- if (length(ev_list) == 0) {
        ev(time = 0, amt = 0, cmt = 1)
    } else {
        do.call(combine_ev, ev_list)
    }

    tumor0_val <- ifelse(paraneoplastic, 1.0, 0.0)

    mod %>%
        param(Ab0 = ab0, kin_Ab = ab0 * 0.05,
              kout_Ab_PE = kout_pe,
              kchemo = chemo_effect,
              Tumor0 = tumor0_val,
              kTumor_growth = tumor_growth) %>%
        mrgsim(ev = ev_use, end = sim_days * 24, delta = 1) %>%
        as.data.frame() %>%
        mutate(time_days = time / 24)
}

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
    skin = "purple",
    dashboardHeader(
        title = "LEMS QSP Dashboard",
        titleWidth = 300
    ),
    dashboardSidebar(
        width = 300,
        sidebarMenu(
            menuItem("Patient Profile",    tabName = "patient",   icon = icon("user")),
            menuItem("Drug PK",            tabName = "pk",        icon = icon("pills")),
            menuItem("VGCC & Antibody",    tabName = "vgcc",      icon = icon("dna")),
            menuItem("NMJ / CMAP",         tabName = "nmj",       icon = icon("bolt")),
            menuItem("Clinical Endpoints", tabName = "clinical",  icon = icon("stethoscope")),
            menuItem("Scenario Comparison",tabName = "scenario",  icon = icon("chart-line")),
            menuItem("Biomarkers",         tabName = "biomarker", icon = icon("vials"))
        )
    ),
    dashboardBody(
        tags$head(tags$style(HTML("
            .box { border-radius: 8px; }
            .skin-purple .main-header .logo { background-color: #6a1b9a; }
            .skin-purple .main-header .navbar { background-color: #7b1fa2; }
            .skin-purple .main-sidebar { background-color: #4a148c; }
        "))),
        tabItems(

            ## ---- TAB 1: Patient Profile ----
            tabItem(tabName = "patient",
                fluidRow(
                    box(title = "LEMS Disease Characteristics", width = 4, status = "purple",
                        h4("Anti-VGCC Antibody"),
                        sliderInput("ab0", "Baseline Anti-VGCC Ab (pmol/L):",
                                    min = 100, max = 10000, value = 2000, step = 100),
                        helpText("Normal < 30 pmol/L; LEMS typically 100–10,000 pmol/L"),
                        hr(),
                        h4("Paraneoplastic LEMS (SCLC)"),
                        checkboxInput("paraneoplastic", "SCLC-associated LEMS", FALSE),
                        conditionalPanel(
                            condition = "input.paraneoplastic == true",
                            sliderInput("tumor_growth", "Tumor Growth Rate:",
                                        min = 0.01, max = 0.1, value = 0.05, step = 0.005),
                            sliderInput("chemo_effect", "Chemotherapy Efficacy:",
                                        min = 0, max = 0.3, value = 0.08, step = 0.01)
                        )
                    ),
                    box(title = "Simulation Settings", width = 4, status = "purple",
                        sliderInput("sim_days", "Simulation Duration (days):",
                                    min = 30, max = 365, value = 180, step = 15),
                        hr(),
                        h4("Plasma Exchange"),
                        sliderInput("pe_sessions", "Number of PE Sessions (per week):",
                                    min = 0, max = 5, value = 0, step = 1),
                        helpText("PE rapidly removes ~50-70% of circulating IgG per session")
                    ),
                    box(title = "LEMS Diagnostic Criteria", width = 4, status = "info",
                        h4("Diagnostic Features"),
                        tags$ul(
                            tags$li("Proximal limb weakness (legs > arms)"),
                            tags$li("Areflexia / hyporeflexia"),
                            tags$li("Autonomic dysfunction"),
                            tags$li("Facilitation on repetitive stimulation"),
                            tags$li("Anti-VGCC antibodies (P/Q-type)")
                        ),
                        hr(),
                        h4("VGCC Subtypes Affected"),
                        tags$ul(
                            tags$li("P/Q-type (Cav2.1): Primary target"),
                            tags$li("N-type (Cav2.2): Secondary"),
                            tags$li("Autonomic VGCC: Cholinergic ganglia")
                        )
                    )
                ),
                fluidRow(
                    box(title = "Disease State Overview", width = 12, status = "warning",
                        plotlyOutput("overview_plot", height = "350px")
                    )
                )
            ),

            ## ---- TAB 2: Drug PK ----
            tabItem(tabName = "pk",
                fluidRow(
                    box(title = "Amifampridine (3,4-DAP)", width = 6, status = "purple",
                        sliderInput("amif_dose", "Dose per Administration (mg):",
                                    min = 0, max = 30, value = 15, step = 2.5),
                        radioButtons("amif_freq", "Dosing Frequency:",
                                     choices = c("BID (12h)" = 12,
                                                 "TID (8h)"  = 8,
                                                 "QID (6h)"  = 6),
                                     selected = 8, inline = TRUE),
                        helpText("Approved range: 5–25 mg per dose, up to 4× daily"),
                        hr(),
                        verbatimTextOutput("amif_pk_summary")
                    ),
                    box(title = "Prednisolone", width = 6, status = "danger",
                        sliderInput("pred_dose", "Daily Dose (mg):",
                                    min = 0, max = 80, value = 0, step = 5),
                        helpText("Typical LEMS dose: 20–60 mg/day; taper after response"),
                        hr(),
                        verbatimTextOutput("pred_pk_summary")
                    )
                ),
                fluidRow(
                    box(title = "Amifampridine PK Profile", width = 6, status = "purple",
                        plotlyOutput("pk_amif_plot", height = "300px")),
                    box(title = "Prednisolone PK Profile", width = 6, status = "danger",
                        plotlyOutput("pk_pred_plot", height = "300px"))
                ),
                fluidRow(
                    box(title = "K+ Channel Blockade vs. Amifampridine Concentration",
                        width = 12, status = "info",
                        plotlyOutput("kblock_plot", height = "300px"))
                )
            ),

            ## ---- TAB 3: VGCC & Antibody ----
            tabItem(tabName = "vgcc",
                fluidRow(
                    box(title = "Anti-VGCC Antibody Dynamics", width = 7, status = "primary",
                        plotlyOutput("ab_plot", height = "350px")),
                    box(title = "VGCC Functional Status", width = 5, status = "warning",
                        plotlyOutput("vgcc_plot", height = "350px"))
                ),
                fluidRow(
                    box(title = "Antibody Titer vs Clinical Effect", width = 12, status = "info",
                        plotlyOutput("ab_cmap_scatter", height = "350px"))
                )
            ),

            ## ---- TAB 4: NMJ / CMAP ----
            tabItem(tabName = "nmj",
                fluidRow(
                    box(title = "CMAP Amplitude Over Time", width = 8, status = "primary",
                        plotlyOutput("cmap_plot", height = "350px")),
                    box(title = "NMJ Safety Factor", width = 4, status = "warning",
                        plotlyOutput("safety_factor_plot", height = "350px"))
                ),
                fluidRow(
                    box(title = "EPP Amplitude", width = 6, status = "info",
                        plotlyOutput("epp_plot", height = "300px")),
                    box(title = "CMAP % of Normal", width = 6, status = "success",
                        plotlyOutput("cmap_pct_plot", height = "300px"))
                )
            ),

            ## ---- TAB 5: Clinical Endpoints ----
            tabItem(tabName = "clinical",
                fluidRow(
                    box(title = "QMG Score Over Time", width = 6, status = "danger",
                        plotlyOutput("qmg_plot", height = "300px"),
                        helpText("QMG: 0 = normal, 39 = severe; LEMS ≥ 10 = meaningful impairment")),
                    box(title = "CMAP Facilitation (Post-Exercise)", width = 6, status = "warning",
                        plotlyOutput("facilitation_plot", height = "300px"))
                ),
                fluidRow(
                    box(title = "Clinical Response Summary", width = 12, status = "info",
                        DTOutput("clinical_table"))
                )
            ),

            ## ---- TAB 6: Scenario Comparison ----
            tabItem(tabName = "scenario",
                fluidRow(
                    box(title = "Scenario Selection", width = 12, status = "purple",
                        checkboxGroupInput("scenarios",
                            "Select Scenarios to Compare:",
                            choices = c(
                                "No Treatment"                          = "none",
                                "Amifampridine Monotherapy (15mg TID)"  = "amif",
                                "Prednisolone 40mg/day"                 = "pred",
                                "Amifampridine + Prednisolone"          = "combo",
                                "PE (5 sessions/week) + Amifampridine"  = "pe",
                                "IVIG + Amifampridine"                  = "ivig"
                            ),
                            selected = c("none", "amif", "combo"),
                            inline = TRUE
                        )
                    )
                ),
                fluidRow(
                    box(title = "CMAP Comparison", width = 6, status = "primary",
                        plotlyOutput("scen_cmap", height = "350px")),
                    box(title = "QMG Score Comparison", width = 6, status = "danger",
                        plotlyOutput("scen_qmg", height = "350px"))
                ),
                fluidRow(
                    box(title = "Anti-VGCC Ab Comparison", width = 6, status = "info",
                        plotlyOutput("scen_ab", height = "300px")),
                    box(title = "VGCC Blockade Comparison", width = 6, status = "warning",
                        plotlyOutput("scen_vgcc", height = "300px"))
                )
            ),

            ## ---- TAB 7: Biomarkers ----
            tabItem(tabName = "biomarker",
                fluidRow(
                    box(title = "Antibody Titer → CMAP Relationship", width = 6, status = "primary",
                        plotlyOutput("bm_ab_cmap", height = "300px")),
                    box(title = "VGCC Blockade → Safety Factor", width = 6, status = "warning",
                        plotlyOutput("bm_vgcc_sf", height = "300px"))
                ),
                fluidRow(
                    box(title = "Dose–Response: Amifampridine", width = 6, status = "purple",
                        plotlyOutput("bm_dr_amif", height = "300px")),
                    box(title = "B Cell / Immune Dynamics", width = 6, status = "success",
                        plotlyOutput("bm_bcell", height = "300px"))
                ),
                fluidRow(
                    box(title = "Biomarker Summary Table", width = 12, status = "info",
                        DTOutput("bm_table"))
                )
            )
        )
    )
)

## ============================================================
## SERVER
## ============================================================
server <- function(input, output, session) {

    ## Reactive: main simulation
    sim_data <- reactive({
        run_sim(
            ab0           = input$ab0,
            amif_dose     = input$amif_dose,
            amif_interval = as.numeric(input$amif_freq),
            pred_dose     = input$pred_dose,
            pe_sessions   = input$pe_sessions,
            chemo_effect  = if (input$paraneoplastic) input$chemo_effect else 0,
            sim_days      = input$sim_days,
            paraneoplastic= input$paraneoplastic,
            tumor_growth  = input$tumor_growth
        )
    })

    ## ---- Overview Plot ----
    output$overview_plot <- renderPlotly({
        d <- sim_data()
        p <- ggplot(d, aes(x = time_days)) +
            geom_line(aes(y = CMAP, color = "CMAP (mV)"), size = 1.2) +
            geom_line(aes(y = Ab_fold * 2, color = "Ab Fold × 2"), size = 1, linetype = "dashed") +
            geom_line(aes(y = VGCC_blocked_pct / 20, color = "VGCC Blocked /20"), size = 1, linetype = "dotted") +
            scale_y_continuous(name = "CMAP (mV) / [scaled]") +
            scale_color_manual(values = c("CMAP (mV)" = "#7b1fa2",
                                          "Ab Fold × 2" = "#e53935",
                                          "VGCC Blocked /20" = "#1565c0")) +
            labs(title = "LEMS Disease State Overview", x = "Time (days)", color = "") +
            theme_bw()
        ggplotly(p)
    })

    ## ---- PK Plots ----
    output$pk_amif_plot <- renderPlotly({
        d <- sim_data() %>% filter(time_days <= 3)
        p <- ggplot(d, aes(x = time_days * 24, y = Camif_ngmL)) +
            geom_line(color = "#7b1fa2", size = 1.2) +
            geom_hline(yintercept = 120, linetype = "dashed", color = "#e53935") +
            annotate("text", x = 10, y = 130, label = "EC50 (120 ng/mL)",
                     size = 3, color = "#e53935") +
            labs(title = "Amifampridine PK (72h)", x = "Time (h)", y = "Concentration (ng/mL)") +
            theme_bw()
        ggplotly(p)
    })

    output$pk_pred_plot <- renderPlotly({
        d <- sim_data() %>% filter(time_days <= 3)
        p <- ggplot(d, aes(x = time_days * 24, y = Cpred_ngmL)) +
            geom_line(color = "#c62828", size = 1.2) +
            labs(title = "Prednisolone PK (72h)", x = "Time (h)", y = "Concentration (ng/mL)") +
            theme_bw()
        ggplotly(p)
    })

    output$kblock_plot <- renderPlotly({
        conc_seq <- seq(0, 500, length.out = 200)
        emax <- 0.85; ec50 <- 120; nh <- 1.5
        kb <- emax * conc_seq^nh / (ec50^nh + conc_seq^nh)
        df_kblock <- data.frame(Conc = conc_seq, Kblock = kb)
        p <- ggplot(df_kblock, aes(x = Conc, y = Kblock)) +
            geom_line(color = "#4a148c", size = 1.5) +
            geom_vline(xintercept = input$amif_dose * 1000 / 35,
                       linetype = "dashed", color = "#7b1fa2") +
            labs(title = "K+ Channel Blockade (Amifampridine Emax Model)",
                 x = "Amifampridine Concentration (ng/mL)",
                 y = "Fractional K+ Channel Blockade") +
            theme_bw()
        ggplotly(p)
    })

    output$amif_pk_summary <- renderText({
        d <- sim_data()
        peak <- max(d$Camif_ngmL[d$time_days <= 1])
        trough <- min(d$Camif_ngmL[d$time_days <= 2 & d$Camif_ngmL > 0])
        kb_peak <- 0.85 * peak^1.5 / (120^1.5 + peak^1.5)
        sprintf("Peak Camif: %.1f ng/mL\nTrough Camif: %.2f ng/mL\nPeak K-Block: %.1f%%\nT1/2: ~2.5 h",
                peak, trough, kb_peak * 100)
    })

    output$pred_pk_summary <- renderText({
        d <- sim_data()
        peak <- max(d$Cpred_ngmL[d$time_days <= 1])
        sprintf("Peak Cpred: %.1f ng/mL\nDose: %g mg/day\nT1/2: ~3-4 h",
                peak, input$pred_dose)
    })

    ## ---- VGCC & Ab Plots ----
    output$ab_plot <- renderPlotly({
        d <- sim_data()
        p <- ggplot(d, aes(x = time_days, y = Ab_fold)) +
            geom_line(color = "#0d47a1", size = 1.3) +
            geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
            labs(title = "Anti-VGCC Ab (Fold of Baseline)",
                 x = "Time (days)", y = "Ab / Baseline") +
            theme_bw()
        ggplotly(p)
    })

    output$vgcc_plot <- renderPlotly({
        d <- sim_data()
        p <- ggplot(d, aes(x = time_days, y = VGCC_blocked_pct)) +
            geom_line(color = "#e53935", size = 1.3) +
            geom_hline(yintercept = 50, linetype = "dashed", color = "grey50") +
            annotate("text", x = max(d$time_days) * 0.1, y = 53,
                     label = "50% blockade threshold", size = 3, color = "grey40") +
            labs(title = "VGCC Blocked (%)",
                 x = "Time (days)", y = "% VGCC Blocked") +
            ylim(0, 100) +
            theme_bw()
        ggplotly(p)
    })

    output$ab_cmap_scatter <- renderPlotly({
        # Cross-sectional Ab vs CMAP
        ab_seq <- seq(100, 8000, by = 200)
        df_ab <- map_dfr(ab_seq, function(ab) {
            s <- mod %>%
                param(Ab0 = ab, kin_Ab = ab * 0.05) %>%
                mrgsim(ev = ev(time = 0, amt = 0, cmt = 1), end = 30 * 24, delta = 24) %>%
                as.data.frame()
            data.frame(Ab0 = ab, CMAP_D30 = tail(s$CMAP, 1))
        })
        p <- ggplot(df_ab, aes(x = Ab0, y = CMAP_D30)) +
            geom_line(color = "#6a1b9a", size = 1.2) +
            geom_point(size = 2, color = "#9c27b0") +
            geom_hline(yintercept = 0.5 * 5, linetype = "dashed", color = "#e53935") +
            annotate("text", x = 1000, y = 2.8, label = "50% CMAP threshold", size = 3) +
            labs(title = "Anti-VGCC Ab Titer vs. CMAP Amplitude (Day 30)",
                 x = "Baseline Anti-VGCC Ab (pmol/L)",
                 y = "CMAP at Day 30 (mV)") +
            theme_bw()
        ggplotly(p)
    })

    ## ---- NMJ / CMAP ----
    output$cmap_plot <- renderPlotly({
        d <- sim_data()
        p <- ggplot(d, aes(x = time_days, y = CMAP)) +
            geom_ribbon(aes(ymin = 0, ymax = CMAP), fill = "#7e57c2", alpha = 0.3) +
            geom_line(color = "#4a148c", size = 1.3) +
            geom_hline(yintercept = 5.0, linetype = "dashed", color = "grey50") +
            geom_hline(yintercept = 0.5, linetype = "dotted", color = "#e53935") +
            annotate("text", x = max(d$time_days) * 0.05, y = 5.2,
                     label = "Normal (5 mV)", size = 3, color = "grey50") +
            labs(title = "CMAP Amplitude Over Time",
                 x = "Time (days)", y = "CMAP (mV)") +
            ylim(0, 6) +
            theme_bw()
        ggplotly(p)
    })

    output$safety_factor_plot <- renderPlotly({
        d <- sim_data()
        p <- ggplot(d, aes(x = time_days, y = Safety_Factor)) +
            geom_line(color = "#f57f17", size = 1.3) +
            geom_hline(yintercept = 1.0, linetype = "dashed", color = "#e53935") +
            annotate("text", x = max(d$time_days) * 0.05, y = 1.1,
                     label = "Threshold (SF = 1)", size = 3, color = "#e53935") +
            labs(title = "Neuromuscular Safety Factor",
                 x = "Time (days)", y = "EPP / Threshold") +
            theme_bw()
        ggplotly(p)
    })

    output$epp_plot <- renderPlotly({
        d <- sim_data()
        p <- ggplot(d, aes(x = time_days, y = EPP_amp)) +
            geom_line(color = "#0097a7", size = 1.2) +
            geom_hline(yintercept = 10, linetype = "dashed", color = "#e53935") +
            annotate("text", x = max(d$time_days) * 0.05, y = 11,
                     label = "Threshold EPP (10 mV)", size = 3) +
            labs(title = "Endplate Potential Amplitude",
                 x = "Time (days)", y = "EPP (mV)") +
            theme_bw()
        ggplotly(p)
    })

    output$cmap_pct_plot <- renderPlotly({
        d <- sim_data()
        p <- ggplot(d, aes(x = time_days, y = CMAP_pct_normal)) +
            geom_line(color = "#2e7d32", size = 1.2) +
            geom_hline(yintercept = 100, linetype = "dashed", color = "grey50") +
            labs(title = "CMAP as % of Normal",
                 x = "Time (days)", y = "CMAP (% of normal)") +
            ylim(0, 110) +
            theme_bw()
        ggplotly(p)
    })

    ## ---- Clinical Endpoints ----
    output$qmg_plot <- renderPlotly({
        d <- sim_data()
        p <- ggplot(d, aes(x = time_days, y = QMG)) +
            geom_line(color = "#b71c1c", size = 1.3) +
            geom_hline(yintercept = 10, linetype = "dashed", color = "#ff7043") +
            annotate("text", x = max(d$time_days) * 0.05, y = 11,
                     label = "Clinically meaningful (QMG ≥ 10)", size = 3) +
            labs(title = "QMG Score (Lower = Better)",
                 x = "Time (days)", y = "QMG Score") +
            ylim(0, 40) +
            theme_bw()
        ggplotly(p)
    })

    output$facilitation_plot <- renderPlotly({
        d <- sim_data()
        # Simulate facilitation ratio: post-exercise CMAP ~2-3× pre-exercise in LEMS
        d <- d %>%
            mutate(
                CMAP_pre  = CMAP,
                CMAP_post = CMAP * pmin(3.0, 1 + 1.5 * VGCC_blocked_pct / 100),
                Facil_ratio = CMAP_post / (CMAP_pre + 0.001)
            )
        p <- ggplot(d, aes(x = time_days)) +
            geom_line(aes(y = CMAP_pre, color = "Pre-exercise"), size = 1.2) +
            geom_line(aes(y = CMAP_post, color = "Post-exercise (10s)"), size = 1.2, linetype = "dashed") +
            scale_color_manual(values = c("Pre-exercise" = "#7b1fa2",
                                          "Post-exercise (10s)" = "#fb8c00")) +
            labs(title = "CMAP Facilitation (Characteristic of LEMS)",
                 x = "Time (days)", y = "CMAP (mV)", color = "") +
            theme_bw()
        ggplotly(p)
    })

    output$clinical_table <- renderDT({
        d <- sim_data()
        timepoints <- c(1, 7, 14, 30, 90, 180)
        tbl <- map_dfr(timepoints, function(tp) {
            idx <- which.min(abs(d$time_days - tp))
            row <- d[idx, ]
            data.frame(
                "Time (days)"    = tp,
                "CMAP (mV)"      = round(row$CMAP, 2),
                "CMAP % Normal"  = round(row$CMAP_pct_normal, 1),
                "QMG Score"      = round(row$QMG, 1),
                "Ab Fold"        = round(row$Ab_fold, 3),
                "VGCC Blocked %" = round(row$VGCC_blocked_pct, 1),
                "Safety Factor"  = round(row$Safety_Factor, 2),
                check.names = FALSE
            )
        })
        datatable(tbl, options = list(pageLength = 10, dom = "t"),
                  rownames = FALSE) %>%
            formatStyle("CMAP (mV)",
                        background = styleColorBar(range(tbl[["CMAP (mV)"]]), "#9c27b0"),
                        backgroundSize = "100% 90%", backgroundRepeat = "no-repeat",
                        backgroundPosition = "center")
    })

    ## ---- Scenario Comparison ----
    scen_data <- reactive({
        req(input$scenarios)
        scenarios_list <- list()

        if ("none" %in% input$scenarios) {
            s <- run_sim(input$ab0, 0, 8, 0, 0, 0, input$sim_days, FALSE, 0.05)
            scenarios_list[["No Treatment"]] <- s
        }
        if ("amif" %in% input$scenarios) {
            s <- run_sim(input$ab0, 15, 8, 0, 0, 0, input$sim_days, FALSE, 0.05)
            scenarios_list[["Amifampridine 15mg TID"]] <- s
        }
        if ("pred" %in% input$scenarios) {
            s <- run_sim(input$ab0, 0, 8, 40, 0, 0, input$sim_days, FALSE, 0.05)
            scenarios_list[["Prednisolone 40mg/day"]] <- s
        }
        if ("combo" %in% input$scenarios) {
            s <- run_sim(input$ab0, 15, 8, 40, 0, 0, input$sim_days, FALSE, 0.05)
            scenarios_list[["Amifampridine + Prednisolone"]] <- s
        }
        if ("pe" %in% input$scenarios) {
            s <- run_sim(input$ab0, 15, 8, 0, 5, 0, input$sim_days, FALSE, 0.05)
            scenarios_list[["PE + Amifampridine"]] <- s
        }
        if ("ivig" %in% input$scenarios) {
            s <- run_sim(input$ab0, 15, 8, 20, 2, 0, input$sim_days, FALSE, 0.05)
            scenarios_list[["IVIG + Amifampridine"]] <- s
        }

        bind_rows(lapply(names(scenarios_list), function(nm) {
            scenarios_list[[nm]] %>% mutate(scenario = nm)
        }))
    })

    plot_scenario <- function(yvar, ytitle, color_map) {
        d <- scen_data()
        p <- ggplot(d, aes_string(x = "time_days", y = yvar, color = "scenario")) +
            geom_line(size = 1.1) +
            scale_color_brewer(palette = "Set1") +
            labs(x = "Time (days)", y = ytitle, color = "Scenario") +
            theme_bw(base_size = 11) +
            theme(legend.position = "bottom",
                  legend.text = element_text(size = 8))
        ggplotly(p) %>% layout(legend = list(orientation = "h", y = -0.3))
    }

    output$scen_cmap  <- renderPlotly(plot_scenario("CMAP",             "CMAP (mV)",           NULL))
    output$scen_qmg   <- renderPlotly(plot_scenario("QMG",              "QMG Score",            NULL))
    output$scen_ab    <- renderPlotly(plot_scenario("Ab_fold",          "Ab / Baseline",        NULL))
    output$scen_vgcc  <- renderPlotly(plot_scenario("VGCC_blocked_pct", "VGCC Blocked (%)",     NULL))

    ## ---- Biomarkers ----
    output$bm_ab_cmap <- renderPlotly({
        ab_seq <- seq(100, 5000, by = 250)
        df <- map_dfr(ab_seq, function(ab) {
            s <- mod %>% param(Ab0 = ab, kin_Ab = ab * 0.05) %>%
                mrgsim(ev = ev(0, 0, 1), end = 30 * 24, delta = 24) %>%
                as.data.frame()
            data.frame(Ab = ab, CMAP = tail(s$CMAP, 1))
        })
        p <- ggplot(df, aes(x = Ab, y = CMAP)) +
            geom_line(color = "#6a1b9a", size = 1.2) +
            labs(title = "Ab Titer → CMAP (Day 30)",
                 x = "Anti-VGCC Ab (pmol/L)", y = "CMAP (mV)") +
            theme_bw()
        ggplotly(p)
    })

    output$bm_vgcc_sf <- renderPlotly({
        d <- sim_data()
        p <- ggplot(d, aes(x = VGCC_blocked_pct, y = Safety_Factor)) +
            geom_path(color = "#e65100", size = 1.2, arrow = arrow(length = unit(0.2, "cm"))) +
            geom_hline(yintercept = 1.0, linetype = "dashed") +
            labs(title = "VGCC Blockade → Safety Factor (trajectory)",
                 x = "VGCC Blocked (%)", y = "Safety Factor") +
            theme_bw()
        ggplotly(p)
    })

    output$bm_dr_amif <- renderPlotly({
        doses <- c(0, 5, 10, 15, 20, 25, 30)
        df <- map_dfr(doses, function(d) {
            ev_dr <- ev(cmt = 1, amt = d, time = seq(0, 14 * 24, by = 8))
            s <- mod %>%
                param(Ab0 = input$ab0, kin_Ab = input$ab0 * 0.05) %>%
                mrgsim(ev = ev_dr, end = 14 * 24, delta = 1) %>% as.data.frame()
            data.frame(Dose = d, CMAP_D14 = tail(s$CMAP, 1), QMG_D14 = tail(s$QMG, 1))
        })
        p <- ggplot(df, aes(x = Dose)) +
            geom_line(aes(y = CMAP_D14, color = "CMAP (mV)"), size = 1.2) +
            geom_line(aes(y = QMG_D14 / 8, color = "QMG /8"), size = 1.2, linetype = "dashed") +
            scale_color_manual(values = c("CMAP (mV)" = "#7b1fa2", "QMG /8" = "#e53935")) +
            labs(title = "Dose–Response: Amifampridine (Day 14)",
                 x = "Amifampridine Dose (mg/dose, TID)",
                 y = "CMAP (mV) / [QMG/8]", color = "") +
            theme_bw()
        ggplotly(p)
    })

    output$bm_bcell <- renderPlotly({
        d <- sim_data()
        p <- ggplot(d, aes(x = time_days, y = Bcell)) +
            geom_line(color = "#2e7d32", size = 1.2) +
            geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey50") +
            labs(title = "B Cell Dynamics (Prednisolone Effect)",
                 x = "Time (days)", y = "B Cell Count (normalized)") +
            theme_bw()
        ggplotly(p)
    })

    output$bm_table <- renderDT({
        d <- sim_data()
        tp_list <- c(1, 7, 14, 30, 90, 180)
        tbl <- map_dfr(tp_list, function(tp) {
            idx <- which.min(abs(d$time_days - tp))
            row <- d[idx, ]
            data.frame(
                Day            = tp,
                "CMAP (mV)"   = round(row$CMAP, 2),
                "QMG"         = round(row$QMG, 1),
                "Ab Fold"     = round(row$Ab_fold, 3),
                "VGCC Blk %" = round(row$VGCC_blocked_pct, 1),
                "Safety F"   = round(row$Safety_Factor, 2),
                "Camif ng/mL"= round(row$Camif_ngmL, 1),
                "Cpred ng/mL"= round(row$Cpred_ngmL, 1),
                check.names = FALSE
            )
        })
        datatable(tbl, options = list(dom = "t", pageLength = 10), rownames = FALSE)
    })
}

## ============================================================
## Launch
## ============================================================
shinyApp(ui = ui, server = server)
