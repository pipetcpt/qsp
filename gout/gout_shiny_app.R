## ============================================================
## Gout QSP Interactive Shiny Application
## Quantitative Systems Pharmacology Dashboard
## Version: 1.0  Date: 2026-06-17
## ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)
library(shinycssloaders)

## ============================================================
## Simulation Engine (lightweight version for Shiny)
## ============================================================

simulate_gout <- function(
    # Patient characteristics
    BW=80, AGE=50, SEX=1, food_score=0.5, ETOH=0, eGFR0=90,
    CKD=FALSE, ABCG2_Q141K=FALSE,
    # ULT Drugs
    Allo_dose=0, Febu_dose=0, Prob_dose=0, Lesi_dose=0,
    # Acute anti-inflammatory
    Colch_dose=0, NSAID_dose=0, GC_dose=0, Ana_dose=0, Cana_dose=0,
    # Flare induction
    induce_flare=FALSE, flare_crystal=5,
    # Simulation duration (days)
    sim_days=365
) {
    dt   <- 0.25  # time step (hours)
    nstep <- sim_days * 24 / dt

    # Parameters
    kprod_UA   <- 0.15 * (BW/70) * (1 + food_score*0.4 + ETOH*0.2)
    kprod_diet <- 0.03 * (1 + food_score)
    GFR        <- eGFR0 * ifelse(SEX==1, 1.0, 0.85)
    Vc_UA      <- 14; Vp_UA <- 28; Vsyn_UA <- 0.5
    ktp_UA     <- 0.08; kpt_UA <- 0.04
    ksyn_in    <- 0.06; ksyn_out <- 0.05

    # XO inhibitor PK parameters
    Ki_Oxy   <- 0.001;   Ki_Febu  <- 0.000001
    IC50_URAT1 <- 5.0;   IC50_Lesi <- 0.1
    IC50_Colch <- 0.0003; Emax_Colch <- 0.85
    IC50_Indo  <- 0.002
    IC50_Ana   <- 0.5

    kcryst    <- 0.002;  kdissolve <- 0.001; sUA_sat <- 6.8
    kIL1b     <- 0.5;    kIL1b_deg <- 0.3
    kTNFa     <- 0.3;    kTNFa_deg <- 0.4
    kPMN_rec  <- 0.2;    kPMN_deg  <- 0.15
    kpain_IL1b<- 0.04;   kpain_res <- 0.1
    CLr_base  <- GFR * 0.001 * 60 * (1 - 0.88 + 0.12 + 0.04)

    # Initial states
    UA_c  <- 6.0 * Vc_UA
    UA_p  <- 6.0 * Vp_UA * 0.8
    UA_s  <- 5.5 * Vsyn_UA
    Cryst <- 0.0
    Toph  <- 0.0
    IL1b  <- 1.0
    TNFa  <- 5.0
    PMN   <- 1.0
    Pain  <- 0.0
    eGFR  <- eGFR0
    JDmg  <- 0.0

    # Dosing schedules (daily dose converted to hourly amount delivered q24h)
    # We use proportional dosing effect for simplicity in this ODE solver
    allo_conc  <- function(t) Allo_dose * 0.80 / 30 * 0.3    # simplified steady-state
    oxy_conc   <- function(t) Allo_dose * 0.80 / 45 * 1.2
    febu_conc  <- function(t) Febu_dose * 0.49 / 12 * 0.5
    prob_conc  <- function(t) Prob_dose * 1.0  / 10 * 0.5
    lesi_conc  <- function(t) Lesi_dose * 0.95 / 14 * 0.5
    colch_conc <- function(t) Colch_dose * 0.45 / 100 * 0.3
    indo_conc  <- function(t) NSAID_dose * 0.90 / 20 * 0.4
    ana_conc   <- function(t) Ana_dose * 0.95 / 8 * 0.6
    cana_conc  <- function(t) Cana_dose * 0.70 / 6 * 0.5

    # Record output every 24h
    rec_int <- round(24 / dt)
    nrec    <- sim_days + 1
    out_day <- out_sUA <- out_cryst <- out_toph <-
        out_IL1b <- out_PMN <- out_pain <- out_eGFR <- numeric(nrec)
    idx <- 1

    for (i in seq_len(nstep)) {
        t <- (i - 1) * dt

        # Flare induction at t=0
        if (induce_flare && i == 1) Cryst <- Cryst + flare_crystal

        C_UA   <- UA_c / Vc_UA
        C_UAs  <- UA_s / Vsyn_UA

        Coxy   <- oxy_conc(t);   Cfebu  <- febu_conc(t)
        Cprob  <- prob_conc(t);  Clesi  <- lesi_conc(t)
        Ccolch <- colch_conc(t); Cindo  <- indo_conc(t)
        Cana   <- ana_conc(t);   Ccana  <- cana_conc(t)

        XO_inh  <- 1 - (1 - Coxy/(Coxy+Ki_Oxy)) * (1 - Cfebu/(Cfebu+Ki_Febu))
        XO_act  <- 1.0 - XO_inh
        URAT1_i <- 1 - (1 - Cprob/(Cprob+IC50_URAT1)) * (1 - Clesi/(Clesi+IC50_Lesi))
        CLr     <- GFR * 0.001 * 60 * (1 - 0.88*(1-URAT1_i) + 0.12 + 0.04)
        CLr     <- CLr * ifelse(CKD, 0.4, 1.0)
        ABCG2f  <- ifelse(ABCG2_Q141K, 0.3, 1.0)

        NLRP3_i  <- Emax_Colch * Ccolch / (Ccolch + IC50_Colch)
        COX_inh  <- Cindo / (Cindo + IC50_Indo)
        Ana_eff  <- Cana / (Cana + IC50_Ana)
        Cana_eff <- Ccana / (Ccana + 0.001)

        # Crystal
        Cf_cryst <- ifelse(C_UAs > sUA_sat, kcryst*(C_UAs - sUA_sat), 0)
        Cd_cryst <- kdissolve * Cryst * (1/(1 + C_UAs/sUA_sat))

        # NLRP3/Inflam
        NLRP3act <- 0.8 * Cryst/(1+Cryst) * (1 - NLRP3_i)
        IL1b_prod <- 0.5 * NLRP3act * (1 - Ana_eff) * (1 - Cana_eff)
        PGE2      <- (1 - COX_inh) * IL1b/(IL1b+10)

        # UA production
        UA_prod <- (kprod_UA + kprod_diet) * XO_act

        # ODEs (Euler)
        dUA_c <- UA_prod - CLr*C_UA - 0.015*ABCG2f*UA_c - ktp_UA*UA_c + kpt_UA*UA_p - ksyn_in*UA_c + ksyn_out*UA_s
        dUA_p <- ktp_UA*UA_c - kpt_UA*UA_p
        dUA_s <- ksyn_in*UA_c - ksyn_out*UA_s - Cf_cryst*Vsyn_UA + Cd_cryst*Vsyn_UA
        dCryst<- Cf_cryst - Cd_cryst - 0.005*Cryst*PMN
        dToph <- 0.0001*Cryst - 0.00005*Toph*max(0, 6-C_UA)
        dIL1b <- IL1b_prod - kIL1b_deg*IL1b
        dTNFa <- kTNFa*NLRP3act - kTNFa_deg*TNFa
        dPMN  <- kPMN_rec*(IL1b/(IL1b+5)) - kPMN_deg*PMN
        dPain <- kpain_IL1b*IL1b + 0.03*PGE2*10 + 0.02*PMN - kpain_res*Pain
        dJDmg <- 0.0001*PMN*Cryst
        deGFR <- -0.00002*max(0, C_UA-6)*eGFR

        UA_c  <- max(0, UA_c  + dUA_c  * dt)
        UA_p  <- max(0, UA_p  + dUA_p  * dt)
        UA_s  <- max(0, UA_s  + dUA_s  * dt)
        Cryst <- max(0, Cryst + dCryst * dt)
        Toph  <- max(0, Toph  + dToph  * dt)
        IL1b  <- max(0, IL1b  + dIL1b  * dt)
        TNFa  <- max(0, TNFa  + dTNFa  * dt)
        PMN   <- max(0, PMN   + dPMN   * dt)
        Pain  <- min(10, max(0, Pain + dPain * dt))
        JDmg  <- JDmg  + dJDmg  * dt
        eGFR  <- max(1,  eGFR  + deGFR  * dt)

        # Record
        if (i %% rec_int == 0 || i == 1) {
            out_day[idx]  <- t / 24
            out_sUA[idx]  <- UA_c / Vc_UA
            out_cryst[idx]<- Cryst
            out_toph[idx] <- Toph
            out_IL1b[idx] <- IL1b
            out_PMN[idx]  <- PMN
            out_pain[idx] <- Pain
            out_eGFR[idx] <- eGFR
            idx <- min(idx + 1, nrec)
        }
    }

    data.frame(
        Day   = out_day[1:idx-1],
        sUA   = out_sUA[1:idx-1],
        Crystal = out_cryst[1:idx-1],
        Tophus  = out_toph[1:idx-1],
        IL1b    = out_IL1b[1:idx-1],
        PMN     = out_PMN[1:idx-1],
        Pain    = out_pain[1:idx-1],
        eGFR    = out_eGFR[1:idx-1]
    )
}

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
    skin = "blue",
    dashboardHeader(title = "Gout QSP Model"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Patient Profile",        tabName="tab_patient",   icon=icon("user")),
            menuItem("Serum Urate (sUA)",      tabName="tab_sua",       icon=icon("tint")),
            menuItem("Inflammation & Pain",    tabName="tab_inflam",    icon=icon("fire")),
            menuItem("Clinical Endpoints",     tabName="tab_outcomes",  icon=icon("chart-line")),
            menuItem("Scenario Comparison",    tabName="tab_scen",      icon=icon("exchange-alt")),
            menuItem("Biomarkers & Genetics",  tabName="tab_bio",       icon=icon("dna")),
            menuItem("References",             tabName="tab_ref",       icon=icon("book"))
        )
    ),
    dashboardBody(
        tabItems(

## ---- TAB 1: Patient Profile ----
tabItem(tabName="tab_patient",
    fluidRow(
        box(title="Patient Demographics", status="primary", solidHeader=TRUE, width=4,
            sliderInput("BW",  "Body Weight (kg)", 40, 150, 80, step=5),
            sliderInput("AGE", "Age (years)",       18, 90,  50, step=1),
            radioButtons("SEX","Sex", choices=c("Male"=1,"Female"=0), inline=TRUE),
            sliderInput("eGFR0","Baseline eGFR (mL/min/1.73m²)", 10, 120, 90, step=5),
            checkboxInput("CKD","CKD present (eGFR correction)", FALSE),
            checkboxInput("ABCG2_Q141K","ABCG2 Q141K polymorphism", FALSE)
        ),
        box(title="Lifestyle & Diet", status="warning", solidHeader=TRUE, width=4,
            sliderInput("food_score","High-purine diet score (0=low, 1=high)", 0, 1, 0.5, step=0.05),
            sliderInput("ETOH","Alcohol (drinks/day)", 0, 10, 0, step=0.5),
            helpText("High-purine diet: red meat, organ meat, shellfish, fructose drinks"),
            helpText("Alcohol increases uric acid production AND decreases excretion.")
        ),
        box(title="Disease State", status="danger", solidHeader=TRUE, width=4,
            checkboxInput("induce_flare","Simulate Acute Flare at Day 0", FALSE),
            sliderInput("flare_crystal","Crystal load at flare", 1, 20, 5, step=1),
            sliderInput("sim_days","Simulation duration (days)", 30, 730, 365, step=30),
            helpText("Baseline serum urate ~6 mg/dL. Crystal flare triggers NLRP3 → IL-1β → Pain.")
        )
    ),
    fluidRow(
        box(title="Urate-Lowering Therapy (ULT)", status="success", solidHeader=TRUE, width=6,
            strong("XO Inhibitors:"),
            sliderInput("Allo_dose",  "Allopurinol dose (mg/day; 0=off)",   0, 800, 0,   step=100),
            sliderInput("Febu_dose",  "Febuxostat dose (mg/day; 0=off)",    0, 120, 0,   step=40),
            hr(),
            strong("Uricosurics:"),
            sliderInput("Prob_dose",  "Probenecid dose (mg/day; 0=off)",    0, 2000, 0,  step=500),
            sliderInput("Lesi_dose",  "Lesinurad dose (mg/day; 0=off)",     0, 400, 0,   step=200)
        ),
        box(title="Acute Anti-Inflammatory Drugs", status="danger", solidHeader=TRUE, width=6,
            strong("Prophylaxis / Acute:"),
            sliderInput("Colch_dose", "Colchicine (mg/day; 0=off)",   0, 1.8, 0, step=0.3),
            sliderInput("NSAID_dose", "Indomethacin (mg/day; 0=off)", 0, 150, 0, step=25),
            sliderInput("GC_dose",    "Prednisolone (mg/day; 0=off)", 0, 60,  0, step=10),
            hr(),
            strong("Biologics (IL-1 inhibitors):"),
            sliderInput("Ana_dose",   "Anakinra (mg/day; 0=off)",     0, 200, 0, step=50),
            sliderInput("Cana_dose",  "Canakinumab (mg/dose; 0=off)", 0, 300, 0, step=150)
        )
    )
),

## ---- TAB 2: Serum Urate ----
tabItem(tabName="tab_sua",
    fluidRow(
        valueBoxOutput("vbox_sUA",      width=3),
        valueBoxOutput("vbox_target",   width=3),
        valueBoxOutput("vbox_XO_inh",   width=3),
        valueBoxOutput("vbox_crystal",  width=3)
    ),
    fluidRow(
        box(title="Serum Urate Over Time", status="primary", solidHeader=TRUE, width=8,
            withSpinner(plotlyOutput("plot_sUA", height="380px"))
        ),
        box(title="sUA Distribution Summary", status="info", solidHeader=TRUE, width=4,
            withSpinner(plotlyOutput("plot_sUA_dist", height="380px"))
        )
    ),
    fluidRow(
        box(title="MSU Crystal Burden in Synovial Fluid", status="warning", solidHeader=TRUE, width=6,
            withSpinner(plotlyOutput("plot_crystal", height="300px"))
        ),
        box(title="Tophus Volume Dynamics", status="warning", solidHeader=TRUE, width=6,
            withSpinner(plotlyOutput("plot_tophus", height="300px"))
        )
    )
),

## ---- TAB 3: Inflammation & Pain ----
tabItem(tabName="tab_inflam",
    fluidRow(
        valueBoxOutput("vbox_IL1b",  width=3),
        valueBoxOutput("vbox_PMN",   width=3),
        valueBoxOutput("vbox_pain",  width=3),
        valueBoxOutput("vbox_dur",   width=3)
    ),
    fluidRow(
        box(title="IL-1β Dynamics (NLRP3 Inflammasome Output)", status="danger", solidHeader=TRUE, width=6,
            withSpinner(plotlyOutput("plot_IL1b", height="320px"))
        ),
        box(title="Neutrophil Influx to Joint", status="warning", solidHeader=TRUE, width=6,
            withSpinner(plotlyOutput("plot_PMN", height="320px"))
        )
    ),
    fluidRow(
        box(title="Acute Gout Pain Score (NRS 0-10)", status="danger", solidHeader=TRUE, width=12,
            withSpinner(plotlyOutput("plot_pain", height="280px"))
        )
    )
),

## ---- TAB 4: Clinical Endpoints ----
tabItem(tabName="tab_outcomes",
    fluidRow(
        box(title="eGFR Trajectory (Renal Protection)", status="primary", solidHeader=TRUE, width=6,
            withSpinner(plotlyOutput("plot_eGFR", height="320px"))
        ),
        box(title="Cumulative Joint Damage", status="warning", solidHeader=TRUE, width=6,
            withSpinner(plotlyOutput("plot_jdmg", height="320px"))
        )
    ),
    fluidRow(
        box(title="Annual Gout Flare Risk Estimate", status="danger", solidHeader=TRUE, width=4,
            withSpinner(plotlyOutput("plot_flare_risk", height="280px"))
        ),
        box(title="Clinical Outcomes Summary Table", status="success", solidHeader=TRUE, width=8,
            withSpinner(DTOutput("tbl_outcomes"))
        )
    )
),

## ---- TAB 5: Scenario Comparison ----
tabItem(tabName="tab_scen",
    fluidRow(
        box(title="Scenario Manager", status="primary", solidHeader=TRUE, width=12,
            helpText("Compare up to 5 pre-defined treatment scenarios."),
            checkboxGroupInput("scen_select", "Select Scenarios:",
                choices = c(
                    "1. Untreated (baseline)"              = "s1",
                    "2. Allopurinol 300mg/day"             = "s2",
                    "3. Febuxostat 80mg/day"               = "s3",
                    "4. Allo 300mg + Lesinurad 200mg"      = "s4",
                    "5. Febu 80mg + Colchicine prophylaxis"= "s5"
                ),
                selected = c("s1","s2","s3"),
                inline = TRUE
            )
        )
    ),
    fluidRow(
        box(title="sUA Comparison", status="primary", solidHeader=TRUE, width=6,
            withSpinner(plotlyOutput("plot_scen_sUA", height="350px"))
        ),
        box(title="Pain Score Comparison", status="danger", solidHeader=TRUE, width=6,
            withSpinner(plotlyOutput("plot_scen_pain", height="350px"))
        )
    ),
    fluidRow(
        box(title="Crystal Burden Comparison", status="warning", solidHeader=TRUE, width=6,
            withSpinner(plotlyOutput("plot_scen_cryst", height="320px"))
        ),
        box(title="eGFR Comparison", status="info", solidHeader=TRUE, width=6,
            withSpinner(plotlyOutput("plot_scen_eGFR", height="320px"))
        )
    )
),

## ---- TAB 6: Biomarkers & Genetics ----
tabItem(tabName="tab_bio",
    fluidRow(
        box(title="Genetic Risk Factors & Transporters", status="info", solidHeader=TRUE, width=6,
            h4("ABCG2 Q141K Polymorphism"),
            p("Reduces intestinal urate secretion by ~70%, increasing serum urate by ~1-1.5 mg/dL."),
            p("Present in ~10% Caucasians, ~25% East Asians."),
            p("Patients with Q141K respond poorly to uricosurics alone."),
            hr(),
            h4("SLC2A9 (GLUT9) Variants"),
            p("Major urate transporter in proximal tubule and gut."),
            p("Accounts for ~2-5 mg/dL sUA variance between individuals."),
            hr(),
            h4("SLC22A12 (URAT1) Variants"),
            p("Loss-of-function mutations cause hypouricemia and exercise-induced acute kidney injury."),
            hr(),
            h4("HLA-B*58:01"),
            p("Associated with allopurinol-induced severe cutaneous adverse reactions (SJS/DRESS)."),
            p("Prevalence: Han Chinese ~6-8%, Thais ~7%, Europeans <1%."),
            p("CPIC guidelines recommend genetic screening before allopurinol in high-risk populations.")
        ),
        box(title="Biomarker Trajectories", status="success", solidHeader=TRUE, width=6,
            withSpinner(plotlyOutput("plot_biomarkers", height="420px")),
            helpText("CRP correlates with IL-6, which is induced by IL-1β during acute flares.")
        )
    ),
    fluidRow(
        box(title="NLRP3 Inflammasome Pathway", status="danger", solidHeader=TRUE, width=12,
            h4("MSU Crystal → NLRP3 Activation Cascade"),
            tags$ol(
                tags$li("MSU crystals phagocytosed by macrophages/neutrophils"),
                tags$li("Lysosomal rupture releases Cathepsin B → NLRP3 signal 2"),
                tags$li("TLR2/TLR4 activation → NF-κB → pro-IL-1β mRNA (signal 1)"),
                tags$li("NLRP3 inflammasome assembly (ASC speck formation)"),
                tags$li("Caspase-1 autoactivation → cleaves pro-IL-1β and pro-IL-18"),
                tags$li("Gasdermin D pore formation → pyroptotic cell death"),
                tags$li("IL-1β release → systemic and local inflammation")
            ),
            h4("Drug Targets in this Cascade:"),
            tags$ul(
                tags$li("Colchicine: inhibits tubulin polymerization → ↓ NLRP3 assembly + ↓ neutrophil migration"),
                tags$li("Anakinra (IL-1Ra): blocks IL-1 receptor → prevents IL-1β signaling"),
                tags$li("Canakinumab (anti-IL-1β): neutralizes circulating IL-1β"),
                tags$li("Rilonacept (IL-1 Trap): decoy receptor for IL-1α and IL-1β")
            )
        )
    )
),

## ---- TAB 7: References ----
tabItem(tabName="tab_ref",
    fluidRow(
        box(title="Key References", status="primary", solidHeader=TRUE, width=12,
            h4("Clinical Trials:"),
            tags$ul(
                tags$li("Becker MA et al. NEJM 2005 — Febuxostat vs Allopurinol (CONFIRMS)"),
                tags$li("Sundy JS et al. JAMA 2011 — Pegloticase in Refractory Gout"),
                tags$li("Terkeltaub RA et al. Arthritis Rheum 2010 — Low-dose Colchicine (AGREE)"),
                tags$li("Saag KG et al. Arthritis Rheumatol 2017 — Lesinurad + XOI (CLEAR 1&2)"),
                tags$li("So A et al. Ann Rheum Dis 2010 — Canakinumab in Acute Gout"),
                tags$li("White WB et al. NEJM 2018 — CARES trial (Febuxostat CV safety)")
            ),
            h4("Mechanistic / QSP:"),
            tags$ul(
                tags$li("Martinon F et al. Nature 2006 — MSU crystals activate NLRP3"),
                tags$li("Dalbeth N et al. Lancet 2016 — Pathophysiology of Gout"),
                tags$li("Richette P & Bardin T. Lancet 2010 — Gout comprehensive review"),
                tags$li("Neogi T et al. Arthritis Rheum 2015 — 2015 Gout Classification Criteria"),
                tags$li("Kobylecki CJ et al. Eur Heart J 2017 — SUA and CV risk (Mendelian Randomization)"),
                tags$li("Stamp LK et al. Rheum 2019 — Oxypurinol PK in CKD patients")
            ),
            h4("Guidelines:"),
            tags$ul(
                tags$li("FitzGerald JD et al. ACR Guidelines 2020 — Management of Gout"),
                tags$li("Richette P et al. Ann Rheum Dis 2017 — EULAR Recommendations for Gout"),
                tags$li("CPIC Guideline — Allopurinol and HLA-B*58:01 (pharmgkb.org)")
            )
        )
    )
)

        ) # end tabItems
    )
)

## ============================================================
## SERVER
## ============================================================
server <- function(input, output, session) {

    sim_result <- reactive({
        simulate_gout(
            BW=input$BW, AGE=input$AGE, SEX=as.numeric(input$SEX),
            food_score=input$food_score, ETOH=input$ETOH, eGFR0=input$eGFR0,
            CKD=input$CKD, ABCG2_Q141K=input$ABCG2_Q141K,
            Allo_dose=input$Allo_dose, Febu_dose=input$Febu_dose,
            Prob_dose=input$Prob_dose, Lesi_dose=input$Lesi_dose,
            Colch_dose=input$Colch_dose, NSAID_dose=input$NSAID_dose,
            GC_dose=input$GC_dose, Ana_dose=input$Ana_dose, Cana_dose=input$Cana_dose,
            induce_flare=input$induce_flare, flare_crystal=input$flare_crystal,
            sim_days=input$sim_days
        )
    })

    # Pre-defined scenarios
    scen_data <- reactive({
        base <- list(BW=80, AGE=50, SEX=1, food_score=0.5, ETOH=1, eGFR0=90,
                     CKD=FALSE, ABCG2_Q141K=FALSE,
                     Colch_dose=0, NSAID_dose=0, GC_dose=0, Ana_dose=0, Cana_dose=0,
                     induce_flare=FALSE, flare_crystal=5, sim_days=input$sim_days)

        run_s <- function(params, label) {
            args <- c(base, params)
            do.call(simulate_gout, args) %>% mutate(Scenario=label)
        }

        all <- list()
        if ("s1" %in% input$scen_select)
            all[["s1"]] <- run_s(list(Allo_dose=0, Febu_dose=0, Prob_dose=0, Lesi_dose=0), "Untreated")
        if ("s2" %in% input$scen_select)
            all[["s2"]] <- run_s(list(Allo_dose=300, Febu_dose=0, Prob_dose=0, Lesi_dose=0), "Allopurinol 300mg")
        if ("s3" %in% input$scen_select)
            all[["s3"]] <- run_s(list(Allo_dose=0, Febu_dose=80, Prob_dose=0, Lesi_dose=0), "Febuxostat 80mg")
        if ("s4" %in% input$scen_select)
            all[["s4"]] <- run_s(list(Allo_dose=300, Febu_dose=0, Prob_dose=0, Lesi_dose=200), "Allo+Lesinurad")
        if ("s5" %in% input$scen_select)
            all[["s5"]] <- run_s(list(Allo_dose=0, Febu_dose=80, Prob_dose=0, Lesi_dose=0, Colch_dose=0.5), "Febu+Colch")
        bind_rows(all)
    })

    # Value boxes
    output$vbox_sUA <- renderValueBox({
        d <- sim_result()
        last <- tail(d$sUA, 1)
        col  <- ifelse(last < 6.0, "green", ifelse(last < 8.0, "yellow", "red"))
        valueBox(round(last, 2), "Final sUA (mg/dL)", icon=icon("tint"), color=col)
    })
    output$vbox_target <- renderValueBox({
        d   <- sim_result()
        pct <- round(mean(d$sUA < 6.0) * 100, 1)
        col <- ifelse(pct > 70, "green", ifelse(pct > 40, "yellow", "red"))
        valueBox(paste0(pct, "%"), "Time with sUA < 6 mg/dL", icon=icon("check"), color=col)
    })
    output$vbox_XO_inh <- renderValueBox({
        Coxy  <- ifelse(input$Allo_dose>0, input$Allo_dose*0.80/45*1.2, 0)
        Cfebu <- ifelse(input$Febu_dose>0, input$Febu_dose*0.49/12*0.5, 0)
        inh   <- round((1-(1-Coxy/(Coxy+0.001))*(1-Cfebu/(Cfebu+0.000001)))*100, 1)
        valueBox(paste0(inh, "%"), "XO Inhibition", icon=icon("ban"), color="blue")
    })
    output$vbox_crystal <- renderValueBox({
        d   <- sim_result()
        last <- round(tail(d$Crystal, 1), 3)
        col  <- ifelse(last < 0.5, "green", ifelse(last < 2, "yellow", "red"))
        valueBox(last, "Crystal Burden", icon=icon("snowflake"), color=col)
    })
    output$vbox_IL1b <- renderValueBox({
        last <- round(tail(sim_result()$IL1b, 1), 2)
        col  <- ifelse(last < 5, "green", ifelse(last < 20, "yellow", "red"))
        valueBox(last, "IL-1β (pg/mL)", icon=icon("fire"), color=col)
    })
    output$vbox_PMN <- renderValueBox({
        last <- round(tail(sim_result()$PMN, 1), 2)
        valueBox(last, "Neutrophil Influx", icon=icon("capsules"), color="orange")
    })
    output$vbox_pain <- renderValueBox({
        last <- round(tail(sim_result()$Pain, 1), 1)
        col  <- ifelse(last < 3, "green", ifelse(last < 6, "yellow", "red"))
        valueBox(last, "NRS Pain Score", icon=icon("bolt"), color=col)
    })
    output$vbox_dur <- renderValueBox({
        d   <- sim_result()
        dur <- sum(d$Pain > 3)  # days with moderate-severe pain
        valueBox(dur, "Days Pain > 3", icon=icon("calendar-times"), color="red")
    })

    # Plots
    make_plot <- function(df, x, y, ylab, title, refline=NULL, color="steelblue") {
        p <- ggplot(df, aes_string(x=x, y=y)) +
            geom_line(color=color, linewidth=1.0) +
            labs(title=title, x="Day", y=ylab) +
            theme_bw(base_size=12)
        if (!is.null(refline))
            p <- p + geom_hline(yintercept=refline, linetype="dashed", color="red")
        p
    }

    output$plot_sUA   <- renderPlotly({ ggplotly(make_plot(sim_result(), "Day","sUA","sUA (mg/dL)","Serum Urate", refline=6.0)) })
    output$plot_sUA_dist <- renderPlotly({
        d <- sim_result()
        p <- ggplot(d, aes(x=sUA)) + geom_histogram(bins=30, fill="steelblue", color="white") +
            geom_vline(xintercept=6.0, linetype="dashed", color="red") +
            labs(title="sUA Distribution", x="sUA (mg/dL)", y="Frequency") + theme_bw()
        ggplotly(p)
    })
    output$plot_crystal <- renderPlotly({ ggplotly(make_plot(sim_result(), "Day","Crystal","Crystal (normalized)","MSU Crystal Burden", color="#FF8F00")) })
    output$plot_tophus  <- renderPlotly({ ggplotly(make_plot(sim_result(), "Day","Tophus","Volume (cm³)","Tophus Volume", color="#9C27B0")) })
    output$plot_IL1b    <- renderPlotly({ ggplotly(make_plot(sim_result(), "Day","IL1b","IL-1β (pg/mL)","IL-1β Dynamics", color="#E53935")) })
    output$plot_PMN     <- renderPlotly({ ggplotly(make_plot(sim_result(), "Day","PMN","PMN (normalized)","Neutrophil Influx", color="#FF7043")) })
    output$plot_pain    <- renderPlotly({ ggplotly(make_plot(sim_result(), "Day","Pain","NRS Pain (0-10)","Gout Pain Score", refline=3, color="#D32F2F")) })
    output$plot_eGFR    <- renderPlotly({ ggplotly(make_plot(sim_result(), "Day","eGFR","eGFR (mL/min/1.73m²)","eGFR Trajectory", refline=60, color="#1565C0")) })
    output$plot_jdmg    <- renderPlotly({ ggplotly(make_plot(sim_result(), "Day","Pain","Damage (a.u.)","Joint Damage Accumulation", color="#795548")) })
    output$plot_flare_risk <- renderPlotly({
        d <- sim_result()
        monthly <- d %>% mutate(Month = ceiling(Day/30)) %>%
            group_by(Month) %>%
            summarise(avg_crystal=mean(Crystal), flare_prob=pmin(avg_crystal*0.15, 1.0))
        p <- ggplot(monthly, aes(x=Month, y=flare_prob*100)) +
            geom_bar(stat="identity", fill="#EF5350") +
            labs(title="Monthly Flare Risk (%)", x="Month", y="Probability (%)") + theme_bw()
        ggplotly(p)
    })

    output$tbl_outcomes <- renderDT({
        d <- sim_result()
        data.frame(
            Metric          = c("Final sUA (mg/dL)","Mean sUA","% Time sUA<6","Crystal Burden",
                                "Tophus Volume","Peak Pain (NRS)","Days Pain>3","Final eGFR"),
            Value           = c(round(tail(d$sUA,1),2), round(mean(d$sUA),2),
                                round(mean(d$sUA<6)*100,1), round(tail(d$Crystal,1),3),
                                round(tail(d$Tophus,1),4), round(max(d$Pain),1),
                                sum(d$Pain>3), round(tail(d$eGFR,1),1)),
            Target          = c("<6.0","<6.0",">80%","<0.1","Resolving","<3","<10","≥60"),
            Status          = c(
                ifelse(tail(d$sUA,1)<6,"✓ MET","✗ NOT MET"),
                ifelse(mean(d$sUA)<6,"✓ MET","✗ NOT MET"),
                ifelse(mean(d$sUA<6)*100>80,"✓ MET","✗ NOT MET"),
                ifelse(tail(d$Crystal,1)<0.1,"✓ MET","✗ NOT MET"),
                ifelse(tail(d$Tophus,1)<tail(head(d,2)$Tophus,1),"Reducing","Increasing"),
                ifelse(max(d$Pain)<3,"✓ MET","✗ NOT MET"),
                ifelse(sum(d$Pain>3)<10,"✓ MET","✗ NOT MET"),
                ifelse(tail(d$eGFR,1)>=60,"✓ PRESERVED","⚠ REDUCED")
            )
        )
    }, options=list(pageLength=10, dom='t'))

    output$plot_biomarkers <- renderPlotly({
        d <- sim_result()
        long <- d %>% select(Day, IL1b, PMN) %>%
            mutate(CRP=IL1b*0.8) %>%
            pivot_longer(-Day, names_to="Marker", values_to="Value")
        p <- ggplot(long, aes(x=Day, y=Value, color=Marker)) +
            geom_line(linewidth=0.9) +
            labs(title="Inflammatory Biomarkers", x="Day", y="Value (normalized/pg/mL)") +
            theme_bw()
        ggplotly(p)
    })

    # Scenario plots
    scen_plot <- function(var, ylab, title, refline=NULL) {
        d <- scen_data()
        if (nrow(d)==0) return(ggplotly(ggplot() + theme_bw()))
        p <- ggplot(d, aes_string(x="Day", y=var, color="Scenario")) +
            geom_line(linewidth=0.9) +
            labs(title=title, x="Day", y=ylab) +
            theme_bw(base_size=12) + theme(legend.position="bottom")
        if (!is.null(refline)) p <- p + geom_hline(yintercept=refline, linetype="dashed", color="red")
        ggplotly(p)
    }
    output$plot_scen_sUA   <- renderPlotly({ scen_plot("sUA",    "sUA (mg/dL)",       "sUA Comparison",    refline=6) })
    output$plot_scen_pain  <- renderPlotly({ scen_plot("Pain",   "NRS Pain",          "Pain Comparison",   refline=3) })
    output$plot_scen_cryst <- renderPlotly({ scen_plot("Crystal","Crystal (norm.)",   "Crystal Comparison") })
    output$plot_scen_eGFR  <- renderPlotly({ scen_plot("eGFR",   "eGFR (mL/min)",    "eGFR Comparison",   refline=60) })
}

## ============================================================
## Launch
## ============================================================
shinyApp(ui, server)
