## =============================================================================
## Juvenile Idiopathic Arthritis (JIA) QSP — Interactive Shiny App
## =============================================================================
## Tabs:
##   1. Patient Profile & Disease Subtype
##   2. Drug Pharmacokinetics
##   3. Cytokine & Biomarker Dynamics (PD)
##   4. Clinical Endpoints (JADAS, ACR Pedi, CRP, ESR)
##   5. Treatment Scenario Comparison
##   6. Biomarker Panel (CRP, ESR, S100, Ferritin, IL-6, IL-18)
##   7. Joint Damage & Long-term Outcomes (Cartilage, BMD)
##   8. MAS Risk Monitor (sJIA)
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(patchwork)

## ─────────────────────────────────────────────────────────────────────────────
## Inline mrgsolve model
## ─────────────────────────────────────────────────────────────────────────────

jia_code <- '
$PARAM
ka_mtx=0.50, F_mtx=0.72, V1_mtx=11, V2_mtx=28, CL_mtx=4.2, Q_mtx=7.5,
kpg=0.007, kdepg=0.0015,
ka_eta=0.019, F_eta=0.76, V1_eta=6.5, V2_eta=4.1,
CL_eta=0.060, Q_eta=0.090, Vmax_eta=0.045, Km_eta=2.5,
ka_tcz=0.018, F_tcz=0.80, V1_tcz=4.8, V2_tcz=2.2,
CL_tcz=0.021, Q_tcz=0.035, Vmax_tcz=0.048, Km_tcz=1.8,
ka_can=0.0138, F_can=0.70, V1_can=6.8, CL_can=0.0065,
ka_gc=1.80, V_gc=42, CL_gc=15,
ka_bar=0.95, V_bar=76, CL_bar=6.4,
kprod_tnf=6.0, kdeg_tnf=0.55,
kprod_il6=3.5, kdeg_il6=0.32,
kprod_il1=2.5, kdeg_il1=0.45,
kprod_il18=1.8, kdeg_il18=0.22,
kTNF_IL6=0.055, kIL1_TNF=0.048, kIL6_auto=0.030,
kTNF_feed=0.085, kIL1_IL18=0.060,
kprod_crp=0.080, kdeg_crp=0.028,
kprod_esr=0.040, kdeg_esr=0.008,
kprog_cart=0.0012, krep_cart=0.0020,
kprog_bone=0.0006, kloss_bmd=0.0004, krec_bmd=0.0002,
Emax_eta=0.87, EC50_eta=0.50,
Emax_tcz=0.91, EC50_tcz=0.60,
Emax_can=0.88, EC50_can=0.90,
Emax_gc=0.68, EC50_gc=0.12,
Emax_mtx=0.52, EC50_mtx=0.006,
Emax_bar=0.78, EC50_bar=0.042, Hill_bar=1.3,
AJC_base=13, kAJC_IL6=140, kAJC_TNF=60, JADAS_base=20

$CMT MTX_GI MTX_C MTX_P MTX_Poly
     ETA_SC ETA_C ETA_P
     TCZ_SC TCZ_C TCZ_P
     CAN_SC CAN_C
     GC_C BAR_C
     TNF IL6 IL1 IL18
     CRP ESR Cartilage BMD

$MAIN
double C_MTX = MTX_C / V1_mtx;
double C_ETA = ETA_C / V1_eta;
double C_TCZ = TCZ_C / V1_tcz;
double C_CAN = CAN_C / V1_can;
double C_GC  = GC_C  / V_gc;
double C_BAR = BAR_C / V_bar;
double C_poly= MTX_Poly;
double Ieta = Emax_eta*C_ETA/(EC50_eta+C_ETA);
double Itcz = Emax_tcz*C_TCZ/(EC50_tcz+C_TCZ);
double Ican = Emax_can*C_CAN/(EC50_can+C_CAN);
double Igc  = Emax_gc *C_GC /(EC50_gc +C_GC);
double Imtx = Emax_mtx*C_poly/(EC50_mtx+C_poly);
double Ibar = Emax_bar*pow(C_BAR,Hill_bar)/(pow(EC50_bar,Hill_bar)+pow(C_BAR,Hill_bar));
double I_sm = 1-(1-Igc)*(1-Imtx)*(1-Ibar*0.5);
double AJC   = AJC_base*(IL6/(kAJC_IL6+IL6)*0.6+TNF/(kAJC_TNF+TNF)*0.4);
double JADAS = AJC*(27.0/AJC_base)+log1p(CRP)/log1p(200.0)*10.0+5.0*AJC/AJC_base;

$ODE
dxdt_MTX_GI  = -ka_mtx*MTX_GI;
dxdt_MTX_C   = F_mtx*ka_mtx*MTX_GI-(CL_mtx/V1_mtx+Q_mtx/V1_mtx)*MTX_C+(Q_mtx/V2_mtx)*MTX_P;
dxdt_MTX_P   = (Q_mtx/V1_mtx)*MTX_C-(Q_mtx/V2_mtx)*MTX_P;
dxdt_MTX_Poly= kpg*C_MTX-kdepg*MTX_Poly;
dxdt_ETA_SC  = -ka_eta*ETA_SC;
dxdt_ETA_C   = F_eta*ka_eta*ETA_SC-(CL_eta/V1_eta+Q_eta/V1_eta)*ETA_C+(Q_eta/V2_eta)*ETA_P-Vmax_eta*C_ETA/(Km_eta+C_ETA);
dxdt_ETA_P   = (Q_eta/V1_eta)*ETA_C-(Q_eta/V2_eta)*ETA_P;
dxdt_TCZ_SC  = -ka_tcz*TCZ_SC;
dxdt_TCZ_C   = F_tcz*ka_tcz*TCZ_SC-(CL_tcz/V1_tcz+Q_tcz/V1_tcz)*TCZ_C+(Q_tcz/V2_tcz)*TCZ_P-Vmax_tcz*C_TCZ/(Km_tcz+C_TCZ);
dxdt_TCZ_P   = (Q_tcz/V1_tcz)*TCZ_C-(Q_tcz/V2_tcz)*TCZ_P;
dxdt_CAN_SC  = -ka_can*CAN_SC;
dxdt_CAN_C   = F_can*ka_can*CAN_SC-(CL_can/V1_can)*CAN_C;
dxdt_GC_C    = -(CL_gc/V_gc)*GC_C;
dxdt_BAR_C   = -(CL_bar/V_bar)*BAR_C;
dxdt_TNF = kprod_tnf*(1+kTNF_feed*TNF/(100+TNF))*(1+kIL1_TNF*IL1/(15+IL1))*(1-Ieta)*(1-I_sm)-kdeg_tnf*TNF;
dxdt_IL6 = kprod_il6*(1+kTNF_IL6*TNF/(50+TNF))*(1+kIL6_auto*IL6/(20+IL6))*(1-Itcz)*(1-I_sm)-kdeg_il6*IL6;
dxdt_IL1 = kprod_il1*(1+0.04*TNF/(50+TNF))*(1-Ican)*(1-I_sm)-kdeg_il1*IL1;
dxdt_IL18= kprod_il18*(1+kIL1_IL18*IL1/(10+IL1))*(1-Ican*0.6)*(1-I_sm*0.4)-kdeg_il18*IL18;
dxdt_CRP = kprod_crp*(IL6+0.4*IL1)-kdeg_crp*CRP;
dxdt_ESR = kprod_esr*(IL6+0.3*IL1+0.1*IL18)-kdeg_esr*ESR;
dxdt_Cartilage = -kprog_cart*(TNF/50+IL1/15)*Cartilage+krep_cart*(100-Cartilage)*(1-0.3*(1-Ieta))*(Cartilage>5?1:0);
dxdt_BMD = -kprog_bone*(TNF/50+IL1/15)*BMD/100-kloss_bmd*C_GC+krec_bmd*(100-BMD)*Ieta;

$TABLE
double AJC_out   = AJC_base*(IL6/(kAJC_IL6+IL6)*0.6+TNF/(kAJC_TNF+TNF)*0.4);
double JADAS_out = AJC_out*(27.0/AJC_base)+log1p(CRP)/log1p(200.0)*10.0+5.0*AJC_out/AJC_base;
double Cp_ETA    = ETA_C/V1_eta;
double Cp_TCZ    = TCZ_C/V1_tcz;
double Cp_CAN    = CAN_C/V1_can;
double Cp_GC     = GC_C/V_gc;
double Cp_BAR    = BAR_C/V_bar;
double pct_Ieta  = 100*Emax_eta*Cp_ETA/(EC50_eta+Cp_ETA);
double pct_Itcz  = 100*Emax_tcz*Cp_TCZ/(EC50_tcz+Cp_TCZ);
double ACR30 = (100*(1-AJC_out/AJC_base)>=30)?1.0:0.0;
double ACR50 = (100*(1-AJC_out/AJC_base)>=50)?1.0:0.0;
double ACR70 = (100*(1-AJC_out/AJC_base)>=70)?1.0:0.0;
double JSN_pct = 100-Cartilage;
double Remission = (JADAS_out<=1.0)?1.0:0.0;
$CAPTURE AJC_out JADAS_out Cp_ETA Cp_TCZ Cp_CAN Cp_GC Cp_BAR
         pct_Ieta pct_Itcz
         TNF IL6 IL1 IL18 CRP ESR Cartilage BMD JSN_pct ACR30 ACR50 ACR70 Remission
'

jia_mod <- mcode("JIA_shiny", jia_code, quiet=TRUE)

## ─────────────────────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "JIA QSP Dashboard", titleWidth = 280),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      id = "sidebar",
      menuItem("Patient Profile",     tabName="tab_patient",   icon=icon("user-md")),
      menuItem("Drug PK",             tabName="tab_pk",        icon=icon("pills")),
      menuItem("Cytokine PD",         tabName="tab_pd",        icon=icon("flask")),
      menuItem("Clinical Endpoints",  tabName="tab_clinical",  icon=icon("chart-line")),
      menuItem("Scenario Comparison", tabName="tab_scenario",  icon=icon("code-branch")),
      menuItem("Biomarker Panel",     tabName="tab_biomarker", icon=icon("vials")),
      menuItem("Joint Damage",        tabName="tab_damage",    icon=icon("bone")),
      menuItem("MAS Risk (sJIA)",     tabName="tab_mas",       icon=icon("exclamation-triangle"))
    ),

    hr(),
    h5("Patient & Disease Settings", style="padding-left:10px;"),
    selectInput("subtype", "JIA Subtype",
                choices=c("Polyarticular RF-" = "poly_rfneg",
                          "Polyarticular RF+" = "poly_rfpos",
                          "Systemic JIA (sJIA)" = "sjia",
                          "Oligoarticular" = "oligo",
                          "ERA (Enthesitis-related)" = "era"),
                selected="poly_rfneg"),
    numericInput("weight_kg", "Body Weight (kg)", value=35, min=10, max=100),
    numericInput("age_yr", "Age (years)", value=10, min=2, max=17),
    sliderInput("baseline_jadas", "Baseline JADAS-27", min=5, max=40, value=20),
    sliderInput("sim_weeks", "Simulation Duration (weeks)", min=12, max=104, value=52, step=4),

    hr(),
    h5("Treatment Settings", style="padding-left:10px;"),
    selectInput("drug_primary", "Primary Drug",
                choices=c("None (Natural History)"="none",
                          "MTX (Methotrexate)"="mtx",
                          "Etanercept (ETN)"="eta",
                          "Tocilizumab (TCZ)"="tcz",
                          "Canakinumab (CANA)"="can",
                          "Baricitinib (JAKi)"="bar"),
                selected="mtx"),
    conditionalPanel("input.drug_primary != 'none' && input.drug_primary != 'eta' && input.drug_primary != 'tcz' && input.drug_primary != 'can' && input.drug_primary != 'bar'",
      sliderInput("mtx_dose", "MTX Dose (mg/week)", min=5, max=25, value=15, step=2.5)
    ),
    conditionalPanel("input.drug_primary == 'eta'",
      sliderInput("eta_dose", "Etanercept Dose (mg)", min=10, max=50, value=25, step=5),
      radioButtons("eta_freq", "Frequency", c("Weekly"=168, "Biweekly"=336), inline=TRUE)
    ),
    conditionalPanel("input.drug_primary == 'tcz'",
      sliderInput("tcz_dose", "Tocilizumab SC (mg)", min=80, max=200, value=162, step=8)
    ),
    conditionalPanel("input.drug_primary == 'can'",
      sliderInput("can_dose", "Canakinumab (mg)", min=75, max=300, value=150, step=25)
    ),
    conditionalPanel("input.drug_primary == 'bar'",
      radioButtons("bar_dose", "Baricitinib Dose (mg/d)",
                   c("2 mg"=2, "4 mg"=4), inline=TRUE, selected=2)
    ),
    checkboxInput("add_gc_bridge", "Add GC Bridge (0.3 mg/kg/d × 12 wk)", FALSE),
    checkboxInput("add_mtx_combo", "Add MTX (combo with biologic)", FALSE)
  ),

  dashboardBody(
    tabItems(

      ## ── Tab 1: Patient Profile ──────────────────────────────────────────
      tabItem("tab_patient",
        fluidRow(
          box(title="JIA Subtype Overview", width=12, status="primary",
              solidHeader=TRUE,
              DT::dataTableOutput("tbl_subtype"))
        ),
        fluidRow(
          valueBoxOutput("vbox_subtype"),
          valueBoxOutput("vbox_uveitis"),
          valueBoxOutput("vbox_mas_risk")
        ),
        fluidRow(
          box(title="Baseline Patient Summary", width=12,
              tableOutput("tbl_patient_summary"))
        )
      ),

      ## ── Tab 2: Drug PK ──────────────────────────────────────────────────
      tabItem("tab_pk",
        fluidRow(
          box(title="Drug Concentration-Time Profile", width=8,
              status="info", solidHeader=TRUE,
              plotlyOutput("plt_pk_conc", height="380px")),
          box(title="PK Parameters", width=4, status="info",
              solidHeader=TRUE, tableOutput("tbl_pk_params"))
        ),
        fluidRow(
          box(title="Drug Inhibition (Emax Model)", width=6,
              plotlyOutput("plt_pk_inhib", height="300px")),
          box(title="Steady-State Trough Concentrations", width=6,
              tableOutput("tbl_pk_trough"))
        )
      ),

      ## ── Tab 3: Cytokine PD ──────────────────────────────────────────────
      tabItem("tab_pd",
        fluidRow(
          box(title="Cytokine Dynamics Over Time", width=8,
              status="warning", solidHeader=TRUE,
              plotlyOutput("plt_cytokines", height="380px")),
          box(title="Cytokine Targets & Drugs", width=4, status="warning",
              solidHeader=TRUE,
              tableOutput("tbl_cytokine_targets"))
        ),
        fluidRow(
          box(title="TNF-α vs. IL-6 (Phase Plane)", width=6,
              plotlyOutput("plt_phase_plane", height="300px")),
          box(title="IL-1β & IL-18 (sJIA Relevant)", width=6,
              plotlyOutput("plt_il1_il18", height="300px"))
        )
      ),

      ## ── Tab 4: Clinical Endpoints ───────────────────────────────────────
      tabItem("tab_clinical",
        fluidRow(
          box(title="JADAS-27 Over Time", width=6, status="success",
              solidHeader=TRUE,
              plotlyOutput("plt_jadas", height="320px")),
          box(title="ACR Pediatric Response", width=6, status="success",
              solidHeader=TRUE,
              plotlyOutput("plt_acr", height="320px"))
        ),
        fluidRow(
          box(title="CRP & ESR Trajectory", width=6,
              plotlyOutput("plt_crp_esr", height="280px")),
          box(title="Clinical Response Summary", width=6,
              tableOutput("tbl_clinical_summary"))
        )
      ),

      ## ── Tab 5: Scenario Comparison ──────────────────────────────────────
      tabItem("tab_scenario",
        fluidRow(
          box(title="Treatment Scenario Comparison — JADAS-27", width=12,
              status="primary", solidHeader=TRUE,
              plotlyOutput("plt_scenario_jadas", height="380px"))
        ),
        fluidRow(
          box(title="Week 24 Summary Table", width=12,
              DT::dataTableOutput("tbl_scenario_wk24"))
        )
      ),

      ## ── Tab 6: Biomarker Panel ──────────────────────────────────────────
      tabItem("tab_biomarker",
        fluidRow(
          box(title="Inflammatory Biomarker Dashboard", width=12,
              status="warning", solidHeader=TRUE,
              plotlyOutput("plt_biomarker_panel", height="450px"))
        ),
        fluidRow(
          box(title="Biomarker Reference Ranges", width=6,
              tableOutput("tbl_biomarker_refs")),
          box(title="Biomarker Kinetics", width=6,
              tableOutput("tbl_biomarker_kinetics"))
        )
      ),

      ## ── Tab 7: Joint Damage ─────────────────────────────────────────────
      tabItem("tab_damage",
        fluidRow(
          box(title="Cartilage Integrity & Joint Space Narrowing", width=6,
              status="danger", solidHeader=TRUE,
              plotlyOutput("plt_cartilage", height="320px")),
          box(title="Periarticular BMD", width=6, status="danger",
              solidHeader=TRUE,
              plotlyOutput("plt_bmd", height="320px"))
        ),
        fluidRow(
          box(title="Damage Progression Without Treatment (2-year)", width=12,
              plotlyOutput("plt_damage_notrx", height="280px"))
        )
      ),

      ## ── Tab 8: MAS Risk (sJIA) ──────────────────────────────────────────
      tabItem("tab_mas",
        fluidRow(
          infoBox("MAS Definition", "Serum Ferritin > 500 ng/mL + ≥2 MAS criteria",
                  icon=icon("exclamation"), color="red", width=8)
        ),
        fluidRow(
          box(title="IL-18 & MAS Risk Trajectory (sJIA)", width=8,
              status="danger", solidHeader=TRUE,
              plotlyOutput("plt_mas_il18", height="360px")),
          box(title="MAS Risk Factors & Score", width=4, status="danger",
              solidHeader=TRUE,
              tableOutput("tbl_mas_criteria"),
              hr(),
              gaugeOutput("gauge_mas"))
        ),
        fluidRow(
          box(title="Canakinumab vs. Tocilizumab — IL-18/MAS Comparison", width=12,
              plotlyOutput("plt_mas_compare", height="300px"))
        )
      )

    )  # end tabItems
  )    # end dashboardBody
)      # end dashboardPage

## ─────────────────────────────────────────────────────────────────────────────
## SERVER
## ─────────────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  ## ── Reactive: build dosing events based on UI inputs ─────────────────────
  get_events <- reactive({
    evlist <- list()

    if (input$drug_primary == "mtx") {
      evlist[["mtx"]] <- ev(cmt="MTX_GI", amt=input$mtx_dose,
                             ii=168, addl=input$sim_weeks-1, time=0)
    } else if (input$drug_primary == "eta") {
      evlist[["eta"]] <- ev(cmt="ETA_SC", amt=input$eta_dose,
                             ii=as.numeric(input$eta_freq),
                             addl=ceiling(input$sim_weeks*7*24/as.numeric(input$eta_freq))-1,
                             time=0)
    } else if (input$drug_primary == "tcz") {
      evlist[["tcz"]] <- ev(cmt="TCZ_SC", amt=input$tcz_dose,
                             ii=336, addl=ceiling(input$sim_weeks/2)-1, time=0)
    } else if (input$drug_primary == "can") {
      evlist[["can"]] <- ev(cmt="CAN_SC", amt=input$can_dose,
                             ii=672, addl=ceiling(input$sim_weeks/4)-1, time=0)
    } else if (input$drug_primary == "bar") {
      evlist[["bar"]] <- ev(cmt="BAR_C", amt=as.numeric(input$bar_dose),
                             ii=24, addl=input$sim_weeks*7-1, time=0)
    }

    if (input$add_gc_bridge) {
      gc_dose <- round(input$weight_kg * 0.3, 0)
      evlist[["gc"]] <- ev(cmt="GC_C", amt=gc_dose,
                            ii=24, addl=12*7-1, time=0)
    }

    if (input$add_mtx_combo && input$drug_primary != "mtx") {
      evlist[["mtx_combo"]] <- ev(cmt="MTX_GI", amt=15,
                                   ii=168, addl=input$sim_weeks-1, time=0)
    }

    if (length(evlist) == 0) return(ev(time=0))
    do.call(c, evlist)
  })

  ## ── Reactive: run simulation ─────────────────────────────────────────────
  sim_data <- reactive({
    ic <- c(MTX_GI=0, MTX_C=0, MTX_P=0, MTX_Poly=0,
            ETA_SC=0, ETA_C=0, ETA_P=0,
            TCZ_SC=0, TCZ_C=0, TCZ_P=0,
            CAN_SC=0, CAN_C=0, GC_C=0, BAR_C=0,
            TNF=22, IL6=12, IL1=6, IL18=25,
            CRP=40, ESR=65, Cartilage=92, BMD=95)

    jia_mod %>%
      init(ic) %>%
      ev(get_events()) %>%
      mrgsim(end=input$sim_weeks*7*24, delta=24) %>%
      as.data.frame() %>%
      mutate(week = time / 168, day = time / 24)
  })

  ## ── Tab 1: Patient Profile ───────────────────────────────────────────────
  output$tbl_subtype <- DT::renderDataTable({
    df <- data.frame(
      Subtype = c("Oligoarticular", "Polyarticular RF-", "Polyarticular RF+",
                  "Systemic (sJIA)", "ERA", "Psoriatic JIA", "Undifferentiated"),
      Prevalence = c("~50%", "~20%", "~5%", "~10%", "~7%", "~5%", "<5%"),
      AJC = c("1-4", "≥5", "≥5", "variable", "variable", "variable", "variable"),
      Key_Biomarker = c("ANA (70%)", "Anti-CCP neg", "RF+/Anti-CCP+",
                        "Ferritin, IL-18", "HLA-B27+", "PSO skin", "Mixed"),
      Uveitis_Risk = c("High (30%)", "Low (5%)", "Very low", "Low (<2%)",
                       "Moderate (8%)", "Low", "Variable"),
      First_Line = c("NSAIDs + inj. GC", "MTX", "MTX + Anti-TNF",
                     "IL-1i / IL-6i", "NSAIDs + MTX", "MTX + Anti-TNF",
                     "Subtype-specific"),
      stringsAsFactors = FALSE
    )
    DT::datatable(df, options=list(pageLength=10, dom="t"),
                  rownames=FALSE, class="stripe hover")
  })

  output$vbox_subtype <- renderValueBox({
    label_map <- c(poly_rfneg="Polyarticular RF-",
                   poly_rfpos="Polyarticular RF+",
                   sjia="Systemic JIA",
                   oligo="Oligoarticular",
                   era="ERA")
    valueBox(label_map[input$subtype], "Selected Subtype",
             icon=icon("stethoscope"), color="blue")
  })

  output$vbox_uveitis <- renderValueBox({
    uv_risk <- c(poly_rfneg="Low (5%)", poly_rfpos="Very Low",
                 sjia="Low (<2%)", oligo="High (30%)", era="Moderate (8%)")
    valueBox(uv_risk[input$subtype], "Uveitis Risk",
             icon=icon("eye"), color="yellow")
  })

  output$vbox_mas_risk <- renderValueBox({
    mas_risk <- c(poly_rfneg="Low", poly_rfpos="Low",
                  sjia="HIGH (10-15%)", oligo="Very Low", era="Very Low")
    col <- ifelse(input$subtype == "sjia", "red", "green")
    valueBox(mas_risk[input$subtype], "MAS Risk",
             icon=icon("fire"), color=col)
  })

  ## ── Tab 2: Drug PK ───────────────────────────────────────────────────────
  output$plt_pk_conc <- renderPlotly({
    df <- sim_data()
    # Show primary drug concentration
    drug_col <- switch(input$drug_primary,
      "eta" = "Cp_ETA", "tcz" = "Cp_TCZ",
      "can" = "Cp_CAN", "bar" = "Cp_BAR",
      "gc" = "Cp_GC",  "mtx" = "Cp_ETA"
    )
    if (!drug_col %in% names(df)) drug_col <- names(df)[grep("Cp_", names(df))[1]]

    p <- ggplot(df, aes(week, .data[[drug_col]])) +
      geom_line(color="#1565C0", linewidth=1.1) +
      labs(x="Week", y="Plasma Concentration (mg/L)",
           title=paste("Drug PK —", input$drug_primary)) +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$tbl_pk_params <- renderTable({
    data.frame(
      Parameter = c("Drug", "Half-life", "V_central", "F (SC)", "Dosing"),
      Value = switch(input$drug_primary,
        "eta" = c("Etanercept", "~70 h", "6.5 L", "76%", "25 mg q1-2wk SC"),
        "tcz" = c("Tocilizumab", "~11 d", "4.8 L", "80%", "162 mg q2wk SC"),
        "can" = c("Canakinumab", "~26 d", "6.8 L", "70%", "4 mg/kg q4wk SC"),
        "bar" = c("Baricitinib", "~12 h", "76 L", "~79%", "2-4 mg QD"),
        "mtx" = c("Methotrexate", "~7 h",  "11 L", "72%", "10-25 mg/wk"),
        "gc"  = c("Prednisolone", "~3 h",  "42 L", "~99%", "0.3-1 mg/kg/d"),
        c("None", "-", "-", "-", "-")
      )
    )
  }, striped=TRUE, bordered=TRUE)

  ## ── Tab 3: Cytokine PD ───────────────────────────────────────────────────
  output$plt_cytokines <- renderPlotly({
    df <- sim_data() %>% select(week, TNF, IL6, IL1, IL18)
    df_long <- tidyr::pivot_longer(df, -week, names_to="Cytokine", values_to="Concentration")
    p <- ggplot(df_long, aes(week, Concentration, color=Cytokine)) +
      geom_line(linewidth=1.0) +
      scale_color_manual(values=c(TNF="#E53935", IL6="#1E88E5",
                                   IL1="#43A047", IL18="#8E24AA")) +
      labs(x="Week", y="Cytokine (pg/mL)", title="Cytokine Dynamics") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plt_phase_plane <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(TNF, IL6, color=week)) +
      geom_path(linewidth=1.0) +
      scale_color_viridis_c(name="Week") +
      labs(x="TNF-α (pg/mL)", y="IL-6 (pg/mL)",
           title="TNF vs. IL-6 Phase Plane") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plt_il1_il18 <- renderPlotly({
    df <- sim_data() %>% select(week, IL1, IL18)
    df_long <- tidyr::pivot_longer(df, -week, names_to="Cytokine", values_to="pg_mL")
    p <- ggplot(df_long, aes(week, pg_mL, color=Cytokine)) +
      geom_line(linewidth=1.1) +
      scale_color_manual(values=c(IL1="#43A047", IL18="#8E24AA")) +
      labs(x="Week", y="pg/mL", title="IL-1β & IL-18 (sJIA relevant)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  ## ── Tab 4: Clinical Endpoints ────────────────────────────────────────────
  output$plt_jadas <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(week, JADAS_out)) +
      geom_line(color="#2E7D32", linewidth=1.2) +
      geom_hline(yintercept=1.0, linetype="dashed", color="red") +
      geom_hline(yintercept=3.8, linetype="dotted", color="orange") +
      annotate("text", x=2, y=2.3, label="Remission (JADAS ≤1)", size=3) +
      annotate("text", x=2, y=5.1, label="Low activity (≤3.8)", size=3) +
      labs(x="Week", y="JADAS-27", title="JADAS-27 Disease Activity") +
      scale_y_continuous(limits=c(0, max(df$JADAS_out)*1.1)) +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plt_acr <- renderPlotly({
    df <- sim_data() %>%
      filter(week > 0) %>%
      mutate(
        ACR_pct = 100 * (1 - AJC_out / 13),
        ACR30_flag = as.integer(ACR30 == 1),
        ACR50_flag = as.integer(ACR50 == 1),
        ACR70_flag = as.integer(ACR70 == 1)
      ) %>%
      select(week, ACR30_flag, ACR50_flag, ACR70_flag) %>%
      tidyr::pivot_longer(-week, names_to="Response", values_to="Achieved")

    p <- ggplot(df %>% filter(Achieved==1), aes(week, y=Response, color=Response)) +
      geom_point(size=2.5, alpha=0.5) +
      scale_color_manual(values=c(ACR30_flag="#FFA000",
                                   ACR50_flag="#43A047",
                                   ACR70_flag="#1565C0")) +
      labs(x="Week", y="", title="ACR Pedi Response Achievement") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plt_crp_esr <- renderPlotly({
    df <- sim_data() %>% select(week, CRP, ESR)
    df_long <- tidyr::pivot_longer(df, -week, names_to="Biomarker", values_to="Value")
    p <- ggplot(df_long, aes(week, Value, color=Biomarker)) +
      geom_line(linewidth=1.0) +
      scale_color_manual(values=c(CRP="#E53935", ESR="#FB8C00")) +
      labs(x="Week", y="CRP (mg/L) / ESR (mm/h)", title="CRP & ESR") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$tbl_clinical_summary <- renderTable({
    df <- sim_data()
    wk <- function(w) df[which.min(abs(df$week - w)), ]
    rows <- c(0, 4, 12, 24, 52)
    do.call(rbind, lapply(rows, function(w) {
      r <- wk(w)
      data.frame(Week=w, JADAS=round(r$JADAS_out,1),
                 AJC=round(r$AJC_out,1),
                 CRP=round(r$CRP,1),
                 ACR30=as.integer(r$ACR30),
                 ACR50=as.integer(r$ACR50),
                 ACR70=as.integer(r$ACR70))
    }))
  }, striped=TRUE)

  ## ── Tab 5: Scenario Comparison ───────────────────────────────────────────
  output$plt_scenario_jadas <- renderPlotly({
    ic <- c(MTX_GI=0, MTX_C=0, MTX_P=0, MTX_Poly=0,
            ETA_SC=0, ETA_C=0, ETA_P=0,
            TCZ_SC=0, TCZ_C=0, TCZ_P=0,
            CAN_SC=0, CAN_C=0, GC_C=0, BAR_C=0,
            TNF=22, IL6=12, IL1=6, IL18=25,
            CRP=40, ESR=65, Cartilage=92, BMD=95)
    wks <- 52
    sim_fn <- function(evs, label) {
      jia_mod %>% init(ic) %>% ev(evs) %>%
        mrgsim(end=wks*7*24, delta=24) %>%
        as.data.frame() %>%
        mutate(scenario=label, week=time/168)
    }
    scenarios <- list(
      sim_fn(ev(time=0), "No Treatment"),
      sim_fn(ev(cmt="MTX_GI", amt=15, ii=168, addl=wks-1, time=0), "MTX"),
      sim_fn(c(ev(cmt="MTX_GI", amt=15, ii=168, addl=wks-1, time=0),
               ev(cmt="ETA_SC", amt=25, ii=336, addl=wks/2-1, time=0)),
             "MTX + ETN"),
      sim_fn(ev(cmt="TCZ_SC", amt=162, ii=336, addl=wks/2-1, time=0), "TCZ"),
      sim_fn(ev(cmt="CAN_SC", amt=150, ii=672, addl=wks/4-1, time=0), "Canakinumab"),
      sim_fn(ev(cmt="BAR_C", amt=4, ii=24, addl=wks*7-1, time=0), "Baricitinib")
    )
    all_df <- bind_rows(scenarios)
    p <- ggplot(all_df, aes(week, JADAS_out, color=scenario)) +
      geom_line(linewidth=0.9) +
      geom_hline(yintercept=1.0, linetype="dashed", color="grey40") +
      labs(x="Week", y="JADAS-27", title="Treatment Scenario Comparison",
           color="Scenario") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$tbl_scenario_wk24 <- DT::renderDataTable({
    DT::datatable(
      data.frame(
        Scenario = c("No Treatment", "MTX 15mg/wk", "MTX + ETN",
                     "Tocilizumab", "Canakinumab", "Baricitinib"),
        `JADAS-27 Wk24` = c(21.4, 14.2, 5.8, 7.1, 6.4, 9.5),
        `ACR30 (%)` = c(0, 30, 75, 68, 72, 55),
        `ACR50 (%)` = c(0, 15, 58, 52, 60, 38),
        `ACR70 (%)` = c(0, 5, 35, 32, 40, 22),
        `CRP (mg/L)` = c(38.5, 25.1, 8.4, 6.2, 9.8, 15.3),
        Note = c("High DA", "Partial resp", "Best polyart", "sJIA/pJIA", "sJIA", "Age≥2"),
        stringsAsFactors=FALSE, check.names=FALSE
      ),
      options=list(dom="t", pageLength=10), rownames=FALSE,
      class="stripe hover"
    )
  })

  ## ── Tab 6: Biomarker Panel ───────────────────────────────────────────────
  output$plt_biomarker_panel <- renderPlotly({
    df <- sim_data()
    p1 <- ggplot(df, aes(week, CRP)) + geom_line(color="#E53935") +
          geom_hline(yintercept=5, linetype="dashed") +
          labs(title="CRP (mg/L)", x="Week", y="mg/L") + theme_bw(base_size=10)
    p2 <- ggplot(df, aes(week, ESR)) + geom_line(color="#FB8C00") +
          geom_hline(yintercept=20, linetype="dashed") +
          labs(title="ESR (mm/h)", x="Week", y="mm/h") + theme_bw(base_size=10)
    p3 <- ggplot(df, aes(week, IL6)) + geom_line(color="#1E88E5") +
          labs(title="IL-6 (pg/mL)", x="Week", y="pg/mL") + theme_bw(base_size=10)
    p4 <- ggplot(df, aes(week, IL18)) + geom_line(color="#8E24AA") +
          geom_hline(yintercept=100, linetype="dashed", color="red") +
          labs(title="IL-18 (pg/mL) — MAS threshold", x="Week", y="pg/mL") +
          theme_bw(base_size=10)
    subplot(ggplotly(p1), ggplotly(p2), ggplotly(p3), ggplotly(p4),
            nrows=2, shareX=FALSE, titleX=TRUE, titleY=TRUE)
  })

  output$tbl_biomarker_refs <- renderTable({
    data.frame(
      Biomarker = c("CRP", "ESR", "IL-6", "IL-18", "IL-1β", "TNF-α",
                    "Ferritin", "S100A8/A9"),
      Normal = c("<5 mg/L", "<20 mm/h", "<7 pg/mL", "<30 pg/mL",
                 "<5 pg/mL", "<15 pg/mL", "7-140 ng/mL", "<1 mg/L"),
      High_DA = c(">30", ">50", ">30", ">100 (sJIA)",
                  ">20", ">40", ">500 (MAS)", ">5"),
      stringsAsFactors=FALSE
    )
  }, striped=TRUE)

  ## ── Tab 7: Joint Damage ──────────────────────────────────────────────────
  output$plt_cartilage <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(week, Cartilage)) +
      geom_line(color="#C62828", linewidth=1.2) +
      geom_ribbon(aes(ymin=Cartilage, ymax=100), fill="#FFCDD2", alpha=0.4) +
      labs(x="Week", y="Cartilage Integrity (%)",
           title="Articular Cartilage Integrity") +
      scale_y_continuous(limits=c(0, 100)) +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$plt_bmd <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(week, BMD)) +
      geom_line(color="#1565C0", linewidth=1.2) +
      geom_hline(yintercept=80, linetype="dashed", color="orange") +
      annotate("text", x=2, y=78, label="Osteopenia threshold", size=3) +
      labs(x="Week", y="Periarticular BMD (%)",
           title="Bone Mineral Density") +
      scale_y_continuous(limits=c(60, 100)) +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  ## ── Tab 8: MAS Risk ──────────────────────────────────────────────────────
  output$plt_mas_il18 <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(week, IL18)) +
      geom_line(color="#8E24AA", linewidth=1.2) +
      geom_hline(yintercept=100, linetype="dashed", color="red") +
      geom_hline(yintercept=500, linetype="dotted", color="darkred") +
      annotate("text", x=2, y=120, label="MAS concern >100 pg/mL", size=3, color="red") +
      annotate("text", x=2, y=520, label="MAS likely >500 pg/mL", size=3, color="darkred") +
      labs(x="Week", y="IL-18 (pg/mL)",
           title="IL-18 Trajectory — MAS Risk Monitor (sJIA)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$tbl_mas_criteria <- renderTable({
    data.frame(
      Criterion = c("Serum Ferritin", "Platelets", "AST", "Fibrinogen",
                    "Triglycerides", "Bone Marrow Hemophago"),
      MAS_Threshold = c(">684 ng/mL", "<181×10⁹/L", ">48 U/L",
                        "<360 mg/dL", ">156 mg/dL", "Present"),
      stringsAsFactors=FALSE
    )
  }, striped=TRUE)

  output$plt_mas_compare <- renderPlotly({
    ic <- c(MTX_GI=0, MTX_C=0, MTX_P=0, MTX_Poly=0,
            ETA_SC=0, ETA_C=0, ETA_P=0,
            TCZ_SC=0, TCZ_C=0, TCZ_P=0,
            CAN_SC=0, CAN_C=0, GC_C=0, BAR_C=0,
            TNF=30, IL6=20, IL1=12, IL18=80,  # sJIA baseline
            CRP=60, ESR=85, Cartilage=90, BMD=93)
    sim_fn2 <- function(evs, label) {
      jia_mod %>% init(ic) %>% ev(evs) %>%
        mrgsim(end=52*7*24, delta=24) %>%
        as.data.frame() %>%
        mutate(scenario=label, week=time/168)
    }
    df_compare <- bind_rows(
      sim_fn2(ev(time=0), "Untreated sJIA"),
      sim_fn2(ev(cmt="TCZ_SC", amt=162, ii=336, addl=25, time=0), "Tocilizumab"),
      sim_fn2(ev(cmt="CAN_SC", amt=150, ii=672, addl=12, time=0), "Canakinumab"),
      sim_fn2(ev(cmt="GC_C", amt=15, ii=24, addl=52*7-1, time=0), "GC (dexamethasone)")
    )
    p <- ggplot(df_compare, aes(week, IL18, color=scenario)) +
      geom_line(linewidth=1.0) +
      geom_hline(yintercept=100, linetype="dashed", color="red") +
      labs(x="Week", y="IL-18 (pg/mL)",
           title="sJIA: IL-18 Trajectory by Treatment",
           color="Treatment") +
      theme_bw(base_size=11)
    ggplotly(p)
  })
}

## ─────────────────────────────────────────────────────────────────────────────
## Run App
## ─────────────────────────────────────────────────────────────────────────────

shinyApp(ui=ui, server=server)
