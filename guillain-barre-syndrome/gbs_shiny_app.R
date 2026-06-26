# ==============================================================================
# Guillain-Barré Syndrome (GBS) QSP — Interactive Shiny Dashboard
# 7 tabs: Patient Profile · PK · Immunology · Nerve Pathology ·
#         Clinical Endpoints · Treatment Comparison · Biomarkers
# ==============================================================================
library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)

# ---- Simplified ODE solver (deSolve-based, no mrgsolve dependency) ----------
library(deSolve)

gbs_ode <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    # Pathogen
    LOS_stim <- ifelse(t >= t_inf_start & t <= t_inf_end, LOS_peak, 0)
    dPathogen <- LOS_stim - k_infect_clear * Pathogen

    # Innate
    dDC_act  <- k_DC_act * Pathogen - k_DC_decay * DC_act
    dMac_act <- k_Mac_recruit * (DC_act + C5a) - k_Mac_decay * Mac_act

    # B cells / Plasma
    BCR_signal <- Pathogen * (1 + DC_act)
    dBcell  <- k_Bcell_act * BCR_signal * (1 - Bcell/2) - 0.05 * Bcell
    dPlasma <- k_plasma_diff * Bcell - 0.03 * Plasma

    # Antibodies
    IVIG_antiidio <- ifelse(IVIG_C > 0.01, 0.5 * IVIG_C / (IVIG_C + 0.1), 0)
    PE_drain      <- PE_flag * PE_Ab_eff * 0.10

    dAb_GM1  <- k_Ab_prod * fGM1  * Plasma - k_Ab_decay * Ab_GM1
                - PE_drain * Ab_GM1 - IVIG_antiidio * Ab_GM1
    dAb_GD1a <- k_Ab_prod * fGD1a * Plasma - k_Ab_decay * Ab_GD1a
                - PE_drain * Ab_GD1a - IVIG_antiidio * Ab_GD1a
    dAb_GQ1b <- k_Ab_GQ1b * fGQ1b * Plasma - k_Ab_decay * Ab_GQ1b
                - PE_drain * Ab_GQ1b

    # T cells
    Treg_inh     <- Treg / (Treg + 0.3)
    IVIG_Treg    <- ifelse(IVIG_C > 0.01, 0.3 * IVIG_C / (IVIG_C + 0.2), 0)
    dTh1  <- k_Th1_diff  * DC_act * (1 - Treg_inh) - k_Tcell_decay * Th1
    dTh17 <- k_Th17_diff * DC_act * (1 - Treg_inh) - k_Tcell_decay * Th17
    dTreg <- k_Treg_form * IL10 - 0.08 * Treg * (Th1 + Th17) + IVIG_Treg

    # Cytokines
    dIL6  <- k_IL6_prod  * (Mac_act + Th17) - k_IL6_decay  * IL6
    dTNFa <- k_TNFa_prod * (Mac_act + Th1)  - k_TNFa_decay * TNFa
    dIL10 <- k_IL10_prod * (Treg + Mac_act) - k_IL10_decay * IL10

    # Complement
    Ab_total     <- Ab_GM1 + Ab_GD1a + Ab_GQ1b
    IVIG_c_inh   <- ifelse(IVIG_C > 0.01, 0.6 * IVIG_C / (IVIG_C + 0.15), 0)
    ECU_inh      <- ifelse(ECU_flag > 0.5 & ECU_C > 0.01, ECU_C / (ECU_C + ECU_IC50), 0)

    dC3b <- k_C3b_prod * Ab_total - k_C3b_decay * C3b - IVIG_c_inh * C3b
    dMAC <- k_MAC_form * C3b * (1 - ECU_inh) - k_MAC_decay * MAC
    dC5a <- k_C5a_prod * C3b * (1 - ECU_inh) - 0.8 * C5a

    # Nerve damage
    if (subtype < 1.5) {
      AIDP_drv <- MAC + Mac_act * 0.5 + IL6 * 0.2; AMAN_drv <- MAC * 0.1
    } else if (subtype < 2.5) {
      AIDP_drv <- MAC * 0.1; AMAN_drv <- MAC + C5a * 0.5
    } else {
      AIDP_drv <- Ab_GQ1b * 0.2; AMAN_drv <- Ab_GQ1b * 0.05
    }
    dMyelin_dmg <- k_myelin_dmg * AIDP_drv * (1 - Myelin_dmg) -
                   k_remyelin * (1 - Myelin_dmg) * Nerve_func
    dAxon_dmg   <- k_axon_dmg * AMAN_drv * (1 - Axon_dmg) -
                   k_axon_repair * (1 - Axon_dmg) * Nerve_func

    dmg_c <- 0.7 * Myelin_dmg + 0.3 * Axon_dmg
    dNerve_func <- -0.5 * dmg_c * Nerve_func +
                   0.03 * (1 - Nerve_func) * (IL10 * 0.05 + 0.05) * (1 - dmg_c)

    dGBS_score <- 0.3 * (6 * (1 - Nerve_func) - GBS_score)
    resp_dmg   <- max(0, (GBS_score - 3) * 0.25)
    dFVC_pct   <- -15 * resp_dmg * (FVC_pct/100) + 2 * Nerve_func * (1 - FVC_pct/100)

    # IVIG PK
    IVIG_rate <- ifelse(IVIG_dose > 0 & t >= IVIG_start & t < (IVIG_start + 5),
                        IVIG_dose * WT / 5, 0)
    CL_t <- IVIG_CL * 24 * WT; Q_t <- IVIG_Q * 24 * WT
    Vc_t <- IVIG_Vc * WT;      Vp_t <- IVIG_Vp * WT
    dIVIG_C <- IVIG_rate/Vc_t - (CL_t/Vc_t)*IVIG_C - (Q_t/Vc_t)*IVIG_C + (Q_t/Vp_t)*IVIG_P
    dIVIG_P <- (Q_t/Vc_t)*IVIG_C - (Q_t/Vp_t)*IVIG_P

    dPE_effect <- PE_flag * 0.1 * (5 - PE_effect)

    ECU_rate <- 0
    if (ECU_flag > 0.5 & t >= IVIG_start)
      ECU_rate <- ifelse(t < IVIG_start + 28, 900/7, 1200/14)
    dECU_C <- ECU_rate/(ECU_Vc*WT) - ECU_CL*24/(ECU_Vc*WT)*ECU_C

    list(c(dPathogen, dDC_act, dMac_act, dBcell, dPlasma,
           dAb_GM1, dAb_GD1a, dAb_GQ1b,
           dTh1, dTh17, dTreg, dIL6, dTNFa, dIL10,
           dC3b, dMAC, dC5a, dMyelin_dmg, dAxon_dmg, dNerve_func,
           dGBS_score, dFVC_pct, dIVIG_C, dIVIG_P, dPE_effect, dECU_C))
  })
}

state_names <- c("Pathogen","DC_act","Mac_act","Bcell","Plasma",
                 "Ab_GM1","Ab_GD1a","Ab_GQ1b",
                 "Th1","Th17","Treg","IL6","TNFa","IL10",
                 "C3b","MAC","C5a","Myelin_dmg","Axon_dmg","Nerve_func",
                 "GBS_score","FVC_pct","IVIG_C","IVIG_P","PE_effect","ECU_C")

default_state <- c(Pathogen=1,DC_act=0.01,Mac_act=0.01,Bcell=0.01,Plasma=0,
                   Ab_GM1=0,Ab_GD1a=0,Ab_GQ1b=0,
                   Th1=0.01,Th17=0.01,Treg=0.1,IL6=0.01,TNFa=0.01,IL10=0.05,
                   C3b=0,MAC=0,C5a=0,Myelin_dmg=0,Axon_dmg=0,Nerve_func=1,
                   GBS_score=0,FVC_pct=100,IVIG_C=0,IVIG_P=0,PE_effect=0,ECU_C=0)

run_gbs <- function(end_day=180, dt=0.5, subtype=1,
                    IVIG_dose=0, IVIG_start=7,
                    PE_flag=0, ECU_flag=0, WT=70) {
  fGM1  <- ifelse(subtype <= 1.5, 1.0, 0.1)
  fGD1a <- ifelse(subtype > 1.5 & subtype < 2.5, 1.0, 0.3)
  fGQ1b <- ifelse(subtype > 2.5, 1.0, 0.05)

  parms <- c(
    k_infect_clear=0.5, LOS_peak=1, t_inf_start=0, t_inf_end=14,
    k_DC_act=0.3, k_DC_decay=0.15, k_Mac_recruit=0.25, k_Mac_decay=0.12,
    k_Bcell_act=0.1, k_plasma_diff=0.08, k_Ab_prod=0.5, k_Ab_decay=0.035,
    k_Ab_GQ1b=0.3, fGM1=fGM1, fGD1a=fGD1a, fGQ1b=fGQ1b,
    k_Th1_diff=0.12, k_Th17_diff=0.08, k_Treg_form=0.05, k_Tcell_decay=0.1,
    k_IL6_prod=0.4, k_IL6_decay=0.5, k_TNFa_prod=0.35, k_TNFa_decay=0.6,
    k_IL10_prod=0.2, k_IL10_decay=0.4,
    k_C3b_prod=0.3, k_C3b_decay=0.2, k_MAC_form=0.25, k_MAC_decay=0.15, k_C5a_prod=0.2,
    k_myelin_dmg=0.08, k_remyelin=0.04, k_axon_dmg=0.06, k_axon_repair=0.015,
    IVIG_CL=0.0033, IVIG_Vc=0.05, IVIG_Vp=0.09, IVIG_Q=0.005, WT=WT,
    PE_flag=PE_flag, PE_Ab_eff=0.6,
    ECU_flag=ECU_flag, ECU_CL=0.0026, ECU_Vc=0.07, ECU_IC50=100,
    IVIG_dose=IVIG_dose, IVIG_start=IVIG_start, subtype=subtype
  )

  times <- seq(0, end_day, by=dt)
  out   <- ode(y=default_state, times=times, func=gbs_ode, parms=parms,
               method="lsoda")
  as.data.frame(out) %>%
    rename(time_day=time) %>%
    mutate(
      Hughes_grade  = case_when(GBS_score < 1 ~ 0, GBS_score < 2 ~ 1,
                                GBS_score < 3 ~ 2, GBS_score < 4 ~ 3,
                                GBS_score < 5 ~ 4, GBS_score < 5.5 ~ 5, TRUE ~ 6),
      Vent_risk     = case_when(FVC_pct < 25 ~ 1, FVC_pct < 40 ~ 0.7,
                                FVC_pct < 60 ~ 0.3, TRUE ~ 0),
      NCV_predicted = 50 * Nerve_func,
      MRC_sum       = 60 * Nerve_func,
      CSF_protein   = 50 + 200 * (1 - Nerve_func) * (Myelin_dmg + Axon_dmg),
      Ab_total      = Ab_GM1 + Ab_GD1a + Ab_GQ1b
    )
}

# ==============================================================================
# UI
# ==============================================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "GBS QSP Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName="tab_patient",  icon=icon("user")),
      menuItem("PK Profile",          tabName="tab_pk",       icon=icon("flask")),
      menuItem("Immunology",          tabName="tab_immune",   icon=icon("bacteria")),
      menuItem("Nerve Pathology",     tabName="tab_nerve",    icon=icon("brain")),
      menuItem("Clinical Endpoints",  tabName="tab_clinical", icon=icon("heartbeat")),
      menuItem("Treatment Comparison",tabName="tab_compare",  icon=icon("chart-bar")),
      menuItem("Biomarkers",          tabName="tab_biomarker",icon=icon("vial"))
    )
  ),
  dashboardBody(
    tabItems(

      # ---- TAB 1: Patient Profile ----
      tabItem("tab_patient",
        h2("Patient Profile & Disease Characteristics"),
        fluidRow(
          box(title="Patient Parameters", width=4, status="primary",
            sliderInput("age",   "Age (years):",       min=10, max=80, value=45, step=1),
            sliderInput("WT",    "Body Weight (kg):",  min=40, max=120, value=70, step=5),
            selectInput("sex",   "Sex:", choices=c("Male","Female")),
            selectInput("subtype","GBS Subtype:",
                        choices=list("AIDP (Demyelinating)"=1,
                                     "AMAN (Axonal Motor)"=2,
                                     "MFS (Miller Fisher)"=3), selected=1),
            selectInput("infection","Antecedent Infection:",
                        choices=c("Campylobacter jejuni","CMV","EBV",
                                  "Mycoplasma pneumoniae","SARS-CoV-2","Unknown"))
          ),
          box(title="Disease Severity", width=4, status="warning",
            sliderInput("days_onset","Days from Symptom Onset:",
                        min=1, max=28, value=7, step=1),
            sliderInput("GBS_baseline","Initial GBS Score (0-6):",
                        min=0, max=6, value=0, step=0.5),
            selectInput("severity","Disease Severity:",
                        choices=c("Mild (Grade 1-2)","Moderate (Grade 3)","Severe (Grade 4-5)")),
            checkboxInput("bulbar","Bulbar involvement", FALSE),
            checkboxInput("autonomic","Autonomic dysfunction", FALSE)
          ),
          box(title="Clinical Assessment", width=4, status="info",
            h4("GBS Prognosis Score (mEGOS)"),
            sliderInput("mEGOS_age",   "Age component (years):", min=10, max=80, value=45),
            sliderInput("mEGOS_GI",    "Preceding GI infection (1=yes):", min=0, max=1, step=1, value=0),
            sliderInput("mEGOS_MRC",   "MRC Sum Score (0-60):", min=0, max=60, value=48),
            h4(textOutput("mEGOS_score")),
            p("mEGOS ≥7: high risk of mechanical ventilation")
          )
        ),
        fluidRow(
          box(title="Disease Overview", width=12, status="success",
            h4("Guillain-Barré Syndrome (GBS) — Key Facts"),
            p("GBS is an acute immune-mediated polyradiculoneuropathy, typically triggered by a preceding
              infection. Campylobacter jejuni is the most common precipitant (~30% of cases). Molecular
              mimicry between microbial lipooligosaccharide (LOS) antigens and gangliosides (GM1, GD1a, GQ1b)
              drives production of pathogenic antibodies that activate complement and damage peripheral nerves."),
            tags$ul(
              tags$li("AIDP (Acute Inflammatory Demyelinating Polyneuropathy): most common in Europe/North America (~90%)"),
              tags$li("AMAN (Acute Motor Axonal Neuropathy): common in Asia/China; anti-GM1 and anti-GD1a antibodies"),
              tags$li("MFS (Miller Fisher Syndrome): ophthalmoplegia, ataxia, areflexia; anti-GQ1b antibodies (~90%)"),
              tags$li("Annual incidence: 1–2 per 100,000; bimodal age distribution (15–35 and 50–75 years)"),
              tags$li("Peak disability at 2–4 weeks; 80% walk independently at 6 months with treatment")
            )
          )
        )
      ),

      # ---- TAB 2: PK Profile ----
      tabItem("tab_pk",
        h2("Pharmacokinetics — IVIG & Eculizumab"),
        fluidRow(
          box(title="Treatment Selection", width=3, status="primary",
            selectInput("treat","Treatment:",
                        choices=c("Untreated","IVIG 2 g/kg","Plasma Exchange (5×)","IVIG + Eculizumab")),
            sliderInput("IVIG_start","IVIG Start Day:", min=1, max=21, value=7, step=1),
            numericInput("IVIG_dose","IVIG Dose (g/kg):", value=2, min=0.4, max=2.4, step=0.4),
            checkboxInput("ECU_on","Add Eculizumab", FALSE),
            p(em("IVIG standard: 2 g/kg over 5 days")),
            p(em("Eculizumab: 900 mg q7d × 4 → 1200 mg q14d"))
          ),
          box(title="IVIG Concentration (Central Compartment)", width=9, status="info",
            plotOutput("pk_plot", height="350px"))
        ),
        fluidRow(
          box(title="IVIG 2-Compartment PK Parameters", width=6, status="warning",
            tableOutput("pk_params_tbl")),
          box(title="PK Summary Statistics", width=6, status="success",
            tableOutput("pk_summary_tbl"))
        )
      ),

      # ---- TAB 3: Immunology ----
      tabItem("tab_immune",
        h2("Immune Cell & Cytokine Dynamics"),
        fluidRow(
          box(title="T Cell Populations (Th1 / Th17 / Treg)", width=6, status="primary",
            plotOutput("tcell_plot", height="300px")),
          box(title="Cytokines (IL-6, TNF-α, IL-10)", width=6, status="warning",
            plotOutput("cytokine_plot", height="300px"))
        ),
        fluidRow(
          box(title="Anti-Ganglioside Antibody Titers", width=6, status="danger",
            plotOutput("ab_plot", height="300px")),
          box(title="Complement Cascade (C3b / MAC / C5a)", width=6, status="info",
            plotOutput("complement_plot", height="300px"))
        )
      ),

      # ---- TAB 4: Nerve Pathology ----
      tabItem("tab_nerve",
        h2("Peripheral Nerve Pathology"),
        fluidRow(
          box(title="Nerve Function Index", width=6, status="primary",
            plotOutput("nerve_func_plot", height="300px")),
          box(title="Myelin vs Axon Damage (subtype-specific)", width=6, status="danger",
            plotOutput("damage_plot", height="300px"))
        ),
        fluidRow(
          box(title="Nerve Conduction Velocity (simulated)", width=6, status="info",
            plotOutput("ncv_plot", height="300px")),
          box(title="Macrophage Infiltration & Innate Immunity", width=6, status="warning",
            plotOutput("innate_plot", height="300px"))
        )
      ),

      # ---- TAB 5: Clinical Endpoints ----
      tabItem("tab_clinical",
        h2("Clinical Endpoints & Outcome Trajectories"),
        fluidRow(
          box(title="GBS Disability Score (Hughes Grade)", width=6, status="danger",
            plotOutput("gbs_score_plot", height="300px")),
          box(title="FVC% & Ventilation Risk", width=6, status="warning",
            plotOutput("fvc_plot", height="300px"))
        ),
        fluidRow(
          box(title="MRC Sum Score (Motor Function)", width=6, status="info",
            plotOutput("mrc_plot", height="300px")),
          box(title="Key Clinical Milestones", width=6, status="success",
            tableOutput("milestones_tbl"))
        )
      ),

      # ---- TAB 6: Treatment Comparison ----
      tabItem("tab_compare",
        h2("Treatment Scenario Comparison"),
        fluidRow(
          box(title="Simulation Duration", width=3, status="primary",
            sliderInput("sim_days","Simulation Days:", min=30, max=365, value=180, step=30),
            sliderInput("comp_subtype","Subtype for Comparison:",
                        min=1, max=3, step=1, value=1,
                        ticks=FALSE),
            p("1=AIDP  2=AMAN  3=MFS")
          ),
          box(title="GBS Score Comparison — All Scenarios", width=9, status="info",
            plotOutput("compare_plot", height="350px"))
        ),
        fluidRow(
          box(title="Clinical Endpoint Summary Table (Day 28 / 90 / 180)", width=12,
            tableOutput("compare_tbl"))
        )
      ),

      # ---- TAB 7: Biomarkers ----
      tabItem("tab_biomarker",
        h2("Prognostic Biomarkers"),
        fluidRow(
          box(title="CSF Protein (mg/dL)", width=6, status="primary",
            plotOutput("csf_plot", height="300px")),
          box(title="Serum Neurofilament Light Chain (NfL proxy)", width=6, status="warning",
            plotOutput("nfl_plot", height="300px"))
        ),
        fluidRow(
          box(title="Anti-Ganglioside Ab as Prognostic Marker", width=6, status="danger",
            plotOutput("ab_prog_plot", height="300px")),
          box(title="Biomarker Reference Ranges", width=6, status="info",
            tableOutput("biomarker_ref_tbl"))
        )
      )
    )
  )
)

# ==============================================================================
# SERVER
# ==============================================================================
server <- function(input, output, session) {

  # Reactive: run simulation based on user inputs
  sim_data <- reactive({
    subtype   <- as.numeric(input$subtype)
    IVIG_dose <- ifelse(input$treat %in% c("IVIG 2 g/kg","IVIG + Eculizumab"), input$IVIG_dose, 0)
    PE_flag   <- ifelse(input$treat == "Plasma Exchange (5×)", 1, 0)
    ECU_flag  <- ifelse(input$treat == "IVIG + Eculizumab" | input$ECU_on, 1, 0)

    run_gbs(end_day=180, dt=0.5,
            subtype=subtype, IVIG_dose=IVIG_dose,
            IVIG_start=input$IVIG_start, PE_flag=PE_flag,
            ECU_flag=ECU_flag, WT=input$WT)
  })

  # Reactive: comparison scenarios
  compare_data <- reactive({
    sub <- input$comp_subtype
    days <- input$sim_days

    sc <- list(
      list(label="Untreated",          IVIG_dose=0,   PE=0, ECU=0),
      list(label="IVIG early (day 7)", IVIG_dose=2,   PE=0, ECU=0, IVIG_start=7),
      list(label="IVIG late (day 14)", IVIG_dose=2,   PE=0, ECU=0, IVIG_start=14),
      list(label="Plasma Exchange",    IVIG_dose=0,   PE=1, ECU=0, IVIG_start=7),
      list(label="IVIG + Eculizumab",  IVIG_dose=2,   PE=0, ECU=1, IVIG_start=7)
    )

    purrr::map_dfr(sc, function(s) {
      d <- run_gbs(end_day=days, dt=1, subtype=sub,
                   IVIG_dose=s$IVIG_dose, IVIG_start=ifelse(is.null(s$IVIG_start), 7, s$IVIG_start),
                   PE_flag=s$PE, ECU_flag=s$ECU, WT=70)
      d$scenario <- s$label
      d
    })
  })

  # mEGOS score
  output$mEGOS_score <- renderText({
    score <- 0
    if (input$mEGOS_age >= 60) score <- score + 2
    else if (input$mEGOS_age >= 40) score <- score + 1
    if (input$mEGOS_GI == 1) score <- score + 2
    score <- score + max(0, round((60 - input$mEGOS_MRC) / 10))
    paste0("mEGOS Score: ", score, " / 12")
  })

  # ---- PK plots ----
  output$pk_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x=time_day, y=IVIG_C)) +
      geom_line(colour="#1F78B4", size=1.3) +
      geom_vline(xintercept=input$IVIG_start, linetype="dashed", colour="red") +
      annotate("text", x=input$IVIG_start+1, y=max(d$IVIG_C)*0.9,
               label=paste0("IVIG start\nday ", input$IVIG_start), colour="red", size=3.5) +
      labs(title="IVIG Central Compartment Concentration",
           x="Time (days)", y="IVIG_C (normalised g/L)") +
      theme_classic(base_size=13)
  })

  output$pk_params_tbl <- renderTable({
    data.frame(
      Parameter = c("CL (L/h/kg)","Vc (L/kg)","Vp (L/kg)","Q (L/h/kg)","T½ (days)","Dose"),
      Value     = c("0.0033","0.05","0.09","0.005","~21","2 g/kg over 5 days"),
      Source    = c("Gelfand 2012","Bhatt 2012","Bhatt 2012","Gelfand 2012","Gelfand 2012","Standard of care")
    )
  })

  output$pk_summary_tbl <- renderTable({
    d <- sim_data()
    data.frame(
      Metric   = c("Peak IVIG_C","Time to Peak","AUC (0-30d)","IVIG_C at Day 28"),
      Value    = c(round(max(d$IVIG_C),3),
                   round(d$time_day[which.max(d$IVIG_C)],1),
                   round(sum(d$IVIG_C[d$time_day <= 30]) * 0.5, 2),
                   round(d$IVIG_C[which.min(abs(d$time_day - 28))], 3))
    )
  })

  # ---- Immunology ----
  output$tcell_plot <- renderPlot({
    d <- sim_data() %>% select(time_day, Th1, Th17, Treg) %>%
      pivot_longer(c(Th1, Th17, Treg), names_to="Cell", values_to="Count")
    ggplot(d, aes(x=time_day, y=Count, colour=Cell)) +
      geom_line(size=1.1) +
      scale_colour_manual(values=c(Th1="#E41A1C", Th17="#FF7F00", Treg="#4DAF4A")) +
      labs(title="T Cell Dynamics", x="Time (days)", y="Count (relative)") +
      theme_bw(base_size=12)
  })

  output$cytokine_plot <- renderPlot({
    d <- sim_data() %>% select(time_day, IL6, TNFa, IL10) %>%
      pivot_longer(c(IL6, TNFa, IL10), names_to="Cytokine", values_to="Level")
    ggplot(d, aes(x=time_day, y=Level, colour=Cytokine)) +
      geom_line(size=1.1) +
      scale_colour_manual(values=c(IL6="#E41A1C", TNFa="#FF7F00", IL10="#4DAF4A")) +
      labs(title="Cytokine Network", x="Time (days)", y="Level (relative)") +
      theme_bw(base_size=12)
  })

  output$ab_plot <- renderPlot({
    d <- sim_data() %>% select(time_day, Ab_GM1, Ab_GD1a, Ab_GQ1b) %>%
      pivot_longer(c(Ab_GM1, Ab_GD1a, Ab_GQ1b), names_to="Ab", values_to="Titer")
    ggplot(d, aes(x=time_day, y=Titer, colour=Ab)) +
      geom_line(size=1.1) +
      scale_colour_manual(values=c(Ab_GM1="#E41A1C", Ab_GD1a="#377EB8", Ab_GQ1b="#4DAF4A")) +
      labs(title="Anti-Ganglioside Antibody Titers", x="Time (days)", y="Titer (relative)") +
      theme_bw(base_size=12)
  })

  output$complement_plot <- renderPlot({
    d <- sim_data() %>% select(time_day, C3b, MAC, C5a) %>%
      pivot_longer(c(C3b, MAC, C5a), names_to="Component", values_to="Level")
    ggplot(d, aes(x=time_day, y=Level, colour=Component)) +
      geom_line(size=1.1) +
      scale_colour_manual(values=c(C3b="#1B9E77", MAC="#D95F02", C5a="#7570B3")) +
      labs(title="Complement Cascade", x="Time (days)", y="Level (relative)") +
      theme_bw(base_size=12)
  })

  # ---- Nerve Pathology ----
  output$nerve_func_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x=time_day, y=Nerve_func)) +
      geom_line(colour="#1F78B4", size=1.3) +
      geom_hline(yintercept=0.5, linetype="dashed", colour="orange") +
      scale_y_continuous(limits=c(0,1)) +
      labs(title="Peripheral Nerve Function (0=loss, 1=normal)",
           x="Time (days)", y="Nerve Function Index") +
      theme_classic(base_size=12)
  })

  output$damage_plot <- renderPlot({
    d <- sim_data() %>% select(time_day, Myelin_dmg, Axon_dmg) %>%
      pivot_longer(c(Myelin_dmg, Axon_dmg), names_to="Type", values_to="Damage")
    ggplot(d, aes(x=time_day, y=Damage, colour=Type)) +
      geom_line(size=1.2) +
      scale_colour_manual(values=c(Myelin_dmg="#E41A1C", Axon_dmg="#377EB8")) +
      scale_y_continuous(limits=c(0,1)) +
      labs(title="Myelin vs Axon Damage", x="Time (days)", y="Damage Index (0-1)") +
      theme_bw(base_size=12)
  })

  output$ncv_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x=time_day, y=NCV_predicted)) +
      geom_line(colour="#984EA3", size=1.3) +
      geom_hline(yintercept=38, linetype="dashed", colour="red") +
      annotate("text", x=30, y=36, label="Lower limit normal (38 m/s)", colour="red", size=3) +
      labs(title="Nerve Conduction Velocity (simulated)",
           x="Time (days)", y="NCV (m/s)") +
      theme_classic(base_size=12)
  })

  output$innate_plot <- renderPlot({
    d <- sim_data() %>% select(time_day, DC_act, Mac_act, Pathogen) %>%
      pivot_longer(c(DC_act, Mac_act, Pathogen), names_to="Cell", values_to="Level")
    ggplot(d, aes(x=time_day, y=Level, colour=Cell)) +
      geom_line(size=1.1) +
      scale_colour_manual(values=c(DC_act="#FF7F00", Mac_act="#E41A1C", Pathogen="#1F78B4")) +
      labs(title="Innate Immune Response & Pathogen Clearance",
           x="Time (days)", y="Level (relative)") +
      theme_bw(base_size=12)
  })

  # ---- Clinical Endpoints ----
  output$gbs_score_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x=time_day, y=GBS_score)) +
      geom_line(colour="#E41A1C", size=1.3) +
      geom_hline(yintercept=3, linetype="dashed", colour="grey40") +
      annotate("text", x=10, y=3.2, label="Grade 3 (hospitalization)", size=3.5) +
      scale_y_continuous(limits=c(0,6), breaks=0:6) +
      labs(title="GBS Disability Score (Hughes Grade 0-6)",
           x="Time (days)", y="GBS Score") +
      theme_classic(base_size=12)
  })

  output$fvc_plot <- renderPlot({
    d <- sim_data() %>% select(time_day, FVC_pct, Vent_risk) %>%
      pivot_longer(c(FVC_pct, Vent_risk), names_to="Measure", values_to="Value") %>%
      mutate(Measure=ifelse(Measure=="FVC_pct","FVC (%)","Ventilation Risk"))
    ggplot(d, aes(x=time_day, y=Value, colour=Measure)) +
      geom_line(size=1.2) +
      facet_wrap(~Measure, scales="free_y") +
      scale_colour_manual(values=c("FVC (%)"="#377EB8","Ventilation Risk"="#E41A1C")) +
      labs(title="Respiratory Function", x="Time (days)", y="Value") +
      theme_bw(base_size=12) + theme(legend.position="none")
  })

  output$mrc_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x=time_day, y=MRC_sum)) +
      geom_line(colour="#4DAF4A", size=1.3) +
      geom_hline(yintercept=48, linetype="dashed", colour="grey40") +
      annotate("text", x=15, y=49, label="MRC 48 (mild weakness)", size=3.5) +
      scale_y_continuous(limits=c(0,60)) +
      labs(title="MRC Sum Score (Motor Strength, max 60)",
           x="Time (days)", y="MRC Sum Score") +
      theme_classic(base_size=12)
  })

  output$milestones_tbl <- renderTable({
    d <- sim_data()
    nadir_day <- d$time_day[which.min(d$Nerve_func)]
    walk_day  <- d %>% filter(time_day > nadir_day, GBS_score < 3) %>%
                  slice(1) %>% pull(time_day)
    indep_day <- d %>% filter(time_day > nadir_day, GBS_score < 2) %>%
                  slice(1) %>% pull(time_day)

    data.frame(
      Milestone             = c("Peak GBS Score","Nadir (worst function)",
                                "GBS Score < 3 (Grade 2)","GBS Score < 2 (Grade 1)",
                                "FVC Nadir (%)","ICU ventilation risk"),
      Value                 = c(round(max(d$GBS_score), 2),
                                round(min(d$Nerve_func), 3),
                                ifelse(length(walk_day)==0,"Not achieved",paste0("Day ",round(walk_day))),
                                ifelse(length(indep_day)==0,"Not achieved",paste0("Day ",round(indep_day))),
                                round(min(d$FVC_pct), 1),
                                round(max(d$Vent_risk), 2))
    )
  })

  # ---- Treatment Comparison ----
  output$compare_plot <- renderPlot({
    d <- compare_data()
    ggplot(d, aes(x=time_day, y=GBS_score, colour=scenario)) +
      geom_line(size=1.1) +
      geom_hline(yintercept=3, linetype="dashed", colour="grey50") +
      scale_colour_brewer(palette="Set1") +
      labs(title="GBS Disability Score — All Scenarios",
           x="Time (days)", y="GBS Score (0-6)", colour="Scenario") +
      theme_classic(base_size=12) +
      theme(legend.position="right")
  })

  output$compare_tbl <- renderTable({
    d <- compare_data()
    d %>% filter(time_day %in% c(14, 28, 90, 180)) %>%
      group_by(Scenario=scenario, Day=time_day) %>%
      summarise(
        GBS_Score   = round(mean(GBS_score), 2),
        Hughes      = round(mean(Hughes_grade), 1),
        FVC_pct     = round(mean(FVC_pct), 1),
        Nerve_Func  = round(mean(Nerve_func), 3),
        MRC_Sum     = round(mean(MRC_sum), 0),
        NCV_ms      = round(mean(NCV_predicted), 1),
        Vent_Risk   = round(mean(Vent_risk), 2),
        .groups="drop"
      )
  })

  # ---- Biomarkers ----
  output$csf_plot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x=time_day, y=CSF_protein)) +
      geom_line(colour="#E41A1C", size=1.3) +
      geom_hline(yintercept=45, linetype="dashed", colour="grey40") +
      annotate("text", x=10, y=47, label="ULN (45 mg/dL)", size=3.5) +
      labs(title="CSF Total Protein (mg/dL)",
           x="Time (days)", y="CSF Protein (mg/dL)") +
      theme_classic(base_size=12)
  })

  output$nfl_plot <- renderPlot({
    d <- sim_data() %>%
      mutate(NfL_proxy = 10 + 500 * Axon_dmg * (1 - Nerve_func))
    ggplot(d, aes(x=time_day, y=NfL_proxy)) +
      geom_line(colour="#FF7F00", size=1.3) +
      geom_hline(yintercept=50, linetype="dashed", colour="grey40") +
      annotate("text", x=15, y=52, label="Approx. threshold for axonal injury", size=3.5) +
      labs(title="Neurofilament Light Chain (NfL, pg/mL proxy)",
           x="Time (days)", y="NfL proxy (pg/mL)") +
      theme_classic(base_size=12)
  })

  output$ab_prog_plot <- renderPlot({
    d <- sim_data() %>% select(time_day, Ab_GM1, Ab_GD1a, Ab_GQ1b) %>%
      pivot_longer(c(Ab_GM1, Ab_GD1a, Ab_GQ1b), names_to="Ab", values_to="Titer")
    ggplot(d, aes(x=time_day, y=Titer, colour=Ab)) +
      geom_line(size=1.1) +
      scale_colour_manual(values=c(Ab_GM1="#E41A1C", Ab_GD1a="#377EB8", Ab_GQ1b="#4DAF4A")) +
      labs(title="Anti-Ganglioside Abs as Prognostic Markers",
           subtitle="Anti-GD1a→AMAN; Anti-GQ1b→MFS; Anti-GM1→poor prognosis AMAN",
           x="Time (days)", y="Titer (relative units)") +
      theme_bw(base_size=12)
  })

  output$biomarker_ref_tbl <- renderTable({
    data.frame(
      Biomarker    = c("CSF Protein","Anti-GM1 IgG","Anti-GD1a IgG","Anti-GQ1b IgG",
                       "NfL (serum)","GFAP (serum)","IL-6 (CSF)","Complement C3"),
      Normal_Range = c("<45 mg/dL","Negative","Negative","Negative",
                       "<10 pg/mL","<0.3 ng/mL","<5 pg/mL","70-150 mg/dL"),
      In_GBS       = c("↑↑ (50-200+)","↑ AIDP/AMAN","↑ AMAN","↑↑ MFS/AMSAN",
                       "↑ axonal GBS","↑ Schwann cell damage","↑ severe GBS","↓ consumed"),
      Prognostic   = c("Moderate","High (AMAN)","High (AMAN)","Dx MFS",
                       "Predicts disability","Research","Research","Moderate")
    )
  })
}

# Run the app
shinyApp(ui=ui, server=server)
