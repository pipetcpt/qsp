## ============================================================
## Recurrent Pericarditis — Interactive QSP Shiny App
## ============================================================
## Tabs: 1-Patient/Scenario | 2-Drug PK | 3-Inflammasome/Cytokines
##       4-Pericardial Pathology | 5-Clinical Endpoints
##       6-Scenario Comparison | 7-Biomarkers
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)

# ---- Inline model code (mirrors rp_mrgsolve_model.R) ----
rp_model_code <- '
$PARAM
ka_colch=0.80, Vd_colch=450, CL_colch=18, k12_colch=0.20, k21_colch=0.10, F_colch=0.45
ka_nsaid=1.20, Vd_nsaid=15, CL_nsaid=2.50, F_nsaid=0.87
ka_cs=1.00, Vd_cs=50, CL_cs=8.00, F_cs=0.82
ka_ana=0.50, Vd_ana=6, CL_ana=1.20, F_ana=0.95
ka_rilo=0.04, Vd_rilo=4, CL_rilo=0.015, F_rilo=0.99
k_nlrp3_on=0.50, k_nlrp3_off=0.08, NLRP3_base=0.05, IC50_colch_nlrp3=0.50
k_il1b_prod=80, k_il1b_deg=0.25, IL1B_base=3, IC50_ana=0.10, IC50_rilo=0.05
k_il18_prod=30, k_il18_deg=0.18, IL18_base=50
k_tnf_prod=25, k_tnf_deg=0.60, TNF_base=5
k_il6_prod=20, k_il6_deg=0.35, IL6_base=2
k_neutro_in=1.50, k_neutro_out=0.12, NEUTRO_base=1, IC50_colch_neutro=0.30
k_m1_on=0.08, k_m1_off=0.04, M1_base=1
k_inflam_on=0.30, k_inflam_off=0.06, INFLAM_base=0.5, INFLAM_max=10
k_eff_prod=3.00, k_eff_resorp=0.025, EFF_base=30
k_fibrin_on=0.40, k_fibrin_off=0.10
k_fibro_prod=0.008, k_fibro_deg=0.004
k_crp_prod=12, k_crp_deg=0.045, CRP_base=1
k_pain_il1=0.80, k_pain_pge2=0.50, k_pain_off=0.20
IC50_nsaid_cox2=15, IC50_cs_nfkb=0.05
DISEASE_ON=1, TRIGGER_STRENGTH=1

$CMT
COLCH_GUT COLCH_CENT COLCH_PERI NSAID_GUT NSAID_CENT
CS_GUT CS_CENT ANA_SC ANA_CENT RILO_SC RILO_CENT
NLRP3_ACT IL1B IL18 TNF IL6 NEUTRO M1_MACRO
INFLAM EFFUSION FIBRIN FIBROSIS CRP PAIN

$INIT
NLRP3_ACT=0.05, IL1B=3, IL18=50, TNF=5, IL6=2, NEUTRO=1, M1_MACRO=1
INFLAM=0.5, EFFUSION=30, FIBRIN=0, FIBROSIS=0, CRP=1, PAIN=0.5

$MAIN
double C_colch = (COLCH_CENT/Vd_colch)*1000;
double C_nsaid = (NSAID_CENT/Vd_nsaid)*1000/206.28*1000;
double C_cs    = (CS_CENT/Vd_cs)*1000/360.44*1000;
double C_ana   = (ANA_CENT/Vd_ana)/17263.0*1e12;
double C_rilo  = (RILO_CENT/Vd_rilo)/251000.0*1e12;
double INH_colch_nlrp3 = C_colch/(IC50_colch_nlrp3+C_colch);
double INH_colch_neutro= C_colch/(IC50_colch_neutro+C_colch);
double INH_nsaid_cox2  = C_nsaid/(IC50_nsaid_cox2+C_nsaid);
double INH_cs_nfkb     = C_cs/(IC50_cs_nfkb+C_cs);
double INH_ana         = C_ana/(IC50_ana+C_ana);
double INH_rilo        = C_rilo/(IC50_rilo+C_rilo);
double IL1_BIO_INH = 1.0-(1.0-INH_ana)*(1.0-INH_rilo);
double NLRP3_stim = DISEASE_ON*TRIGGER_STRENGTH*(1.0-INH_colch_nlrp3);
double PGE2  = (1.0-INH_nsaid_cox2)*INFLAM/INFLAM_max;
double NFKB_act = (1.0-INH_cs_nfkb)*INFLAM/INFLAM_max;

$ODE
dxdt_COLCH_GUT  = -ka_colch*COLCH_GUT;
dxdt_COLCH_CENT = ka_colch*F_colch*COLCH_GUT-(CL_colch/Vd_colch+k12_colch)*COLCH_CENT+k21_colch*COLCH_PERI;
dxdt_COLCH_PERI = k12_colch*COLCH_CENT-k21_colch*COLCH_PERI;
dxdt_NSAID_GUT  = -ka_nsaid*NSAID_GUT;
dxdt_NSAID_CENT = ka_nsaid*F_nsaid*NSAID_GUT-(CL_nsaid/Vd_nsaid)*NSAID_CENT;
dxdt_CS_GUT  = -ka_cs*CS_GUT;
dxdt_CS_CENT = ka_cs*F_cs*CS_GUT-(CL_cs/Vd_cs)*CS_CENT;
dxdt_ANA_SC   = -ka_ana*ANA_SC;
dxdt_ANA_CENT = ka_ana*F_ana*ANA_SC-(CL_ana/Vd_ana)*ANA_CENT;
dxdt_RILO_SC   = -ka_rilo*RILO_SC;
dxdt_RILO_CENT = ka_rilo*F_rilo*RILO_SC-(CL_rilo/Vd_rilo)*RILO_CENT;
dxdt_NLRP3_ACT = k_nlrp3_on*NLRP3_stim*(1.0-NLRP3_ACT)-k_nlrp3_off*NLRP3_ACT;
double IL1B_prod = k_il1b_prod*NLRP3_ACT*(1.0-IL1_BIO_INH);
dxdt_IL1B = IL1B_prod - k_il1b_deg*IL1B;
dxdt_IL18 = k_il18_prod*NLRP3_ACT - k_il18_deg*IL18;
double TNF_prod = k_tnf_prod*NFKB_act*(1.0-INH_cs_nfkb);
dxdt_TNF = TNF_prod - k_tnf_deg*TNF;
double IL6_prod = k_il6_prod*IL1B/(IL1B_base+IL1B)*NFKB_act;
dxdt_IL6 = IL6_prod - k_il6_deg*IL6;
dxdt_NEUTRO = k_neutro_in*INFLAM/INFLAM_max*(1.0-INH_colch_neutro) - k_neutro_out*NEUTRO;
dxdt_M1_MACRO = k_m1_on*INFLAM/INFLAM_max - k_m1_off*M1_MACRO;
double INFLAM_drive = (IL1B/100.0)+(TNF/200.0)+(NEUTRO/5.0)+(M1_MACRO/10.0);
dxdt_INFLAM = k_inflam_on*INFLAM_drive*(INFLAM_max-INFLAM) - k_inflam_off*INFLAM;
dxdt_EFFUSION = k_eff_prod*INFLAM/INFLAM_max - k_eff_resorp*(EFFUSION-EFF_base);
dxdt_FIBRIN = k_fibrin_on*INFLAM/INFLAM_max*(1.0-FIBRIN) - k_fibrin_off*FIBRIN;
dxdt_FIBROSIS = k_fibro_prod*FIBRIN - k_fibro_deg*FIBROSIS;
double CRP_prod = k_crp_prod*IL6/(IL6_base+IL6);
dxdt_CRP = CRP_prod - k_crp_deg*(CRP-CRP_base);
double PAIN_drive = k_pain_il1*(IL1B/100.0)+k_pain_pge2*PGE2*INFLAM/INFLAM_max;
double PAIN_max10 = fmin(PAIN_drive*INFLAM_max, 10.0);
dxdt_PAIN = k_inflam_on*(PAIN_max10-PAIN) - k_pain_off*PAIN;

$CAPTURE C_colch C_nsaid C_cs C_ana C_rilo INH_colch_nlrp3 INH_nsaid_cox2 INH_cs_nfkb IL1_BIO_INH PGE2 NFKB_act
'

# Pre-compile model (done once at startup)
mod <- mcode("rp_shiny", rp_model_code, quiet = TRUE)

# ---- Scenario event builder ----
build_dose <- function(use_colch, dose_colch, use_nsaid, dose_nsaid,
                       use_cs, dose_cs, use_ana, use_rilo) {
  ev_list <- list()
  if (use_colch) ev_list[[length(ev_list)+1]] <- ev(amt=dose_colch, ii=12, addl=179, cmt="COLCH_GUT")
  if (use_nsaid) ev_list[[length(ev_list)+1]] <- ev(amt=dose_nsaid, ii=8,  addl=83,  cmt="NSAID_GUT")
  if (use_cs) {
    ev_list[[length(ev_list)+1]] <- ev(amt=dose_cs,         ii=24, addl=27, cmt="CS_GUT", time=0)
    ev_list[[length(ev_list)+1]] <- ev(amt=dose_cs*0.70,    ii=24, addl=27, cmt="CS_GUT", time=672)
    ev_list[[length(ev_list)+1]] <- ev(amt=dose_cs*0.40,    ii=24, addl=27, cmt="CS_GUT", time=1344)
    ev_list[[length(ev_list)+1]] <- ev(amt=dose_cs*0.15,    ii=24, addl=27, cmt="CS_GUT", time=2016)
  }
  if (use_ana)  ev_list[[length(ev_list)+1]] <- ev(amt=100, ii=24, addl=179, cmt="ANA_SC")
  if (use_rilo) {
    ev_list[[length(ev_list)+1]] <- ev(amt=320, time=0,   cmt="RILO_SC")
    ev_list[[length(ev_list)+1]] <- ev(amt=160, ii=168, addl=25, time=168, cmt="RILO_SC")
  }
  if (length(ev_list) == 0) return(ev(amt=0, time=0, cmt=1))
  do.call(c, ev_list)
}

# ---- Colours ----
SCEN_COLS <- c("Untreated"="#e74c3c","NSAID"="#f39c12","Colchicine"="#3498db",
               "Colch+NSAID"="#2ecc71","Prednisone"="#9b59b6","Anakinra"="#1abc9c","Rilonacept"="#e67e22")

# ---- UI ----
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "Recurrent Pericarditis QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("① Patient & Scenario",  tabName = "tab_patient",  icon = icon("user")),
      menuItem("② Drug PK",             tabName = "tab_pk",       icon = icon("flask")),
      menuItem("③ Inflammasome/Cytokines", tabName="tab_cytokine", icon=icon("fire")),
      menuItem("④ Pericardial Pathology", tabName="tab_pericardial", icon=icon("heart")),
      menuItem("⑤ Clinical Endpoints",  tabName = "tab_endpoints", icon = icon("chart-line")),
      menuItem("⑥ Scenario Comparison", tabName = "tab_compare",  icon = icon("layer-group")),
      menuItem("⑦ Biomarkers",          tabName = "tab_biomarkers",icon = icon("vials"))
    )
  ),

  dashboardBody(
    tabItems(

      ## ---- TAB 1: Patient Profile & Scenario ----
      tabItem("tab_patient",
        fluidRow(
          box(title="Patient Profile", status="danger", solidHeader=TRUE, width=4,
            numericInput("weight_kg", "Body Weight (kg)", 70, 40, 120),
            selectInput("sex", "Sex", c("Male","Female")),
            numericInput("age", "Age (years)", 45, 18, 80),
            selectInput("etiology",  "Presumed Etiology",
              c("Idiopathic","Post-viral","Post-cardiac injury (PCIS)","Autoinflammatory (CAPS/FMF)","Autoimmune (SLE/RA)","Post-COVID-19")),
            numericInput("n_prior", "Number of Prior Episodes", 0, 0, 10),
            sliderInput("trigger", "Trigger Strength (0-2)", 0.5, 0, 2, 0.1)
          ),
          box(title="Treatment Selection", status="warning", solidHeader=TRUE, width=4,
            checkboxInput("use_colch", "Colchicine", TRUE),
            conditionalPanel("input.use_colch",
              numericInput("dose_colch", "Colchicine Dose (mg, per dose)", 0.5, 0.25, 1)),
            checkboxInput("use_nsaid", "NSAID (Ibuprofen)", TRUE),
            conditionalPanel("input.use_nsaid",
              numericInput("dose_nsaid", "Ibuprofen Dose (mg, per dose)", 600, 200, 800)),
            checkboxInput("use_cs", "Prednisone", FALSE),
            conditionalPanel("input.use_cs",
              numericInput("dose_cs", "Prednisone start dose (mg/d)", 35, 5, 80)),
            checkboxInput("use_ana",  "Anakinra 100 mg/d SC", FALSE),
            checkboxInput("use_rilo", "Rilonacept (320→160 mg qw)", FALSE),
            actionButton("run_sim", "▶ Run Simulation", class="btn-danger btn-lg")
          ),
          box(title="ESC 2015 Diagnostic Criteria", status="info", solidHeader=TRUE, width=4,
            tags$p(tags$b("≥2 of the following:")),
            tags$ul(
              tags$li("Pericarditic chest pain"),
              tags$li("Pericardial friction rub"),
              tags$li("New ST-elevation or PR depression"),
              tags$li("Pericardial effusion (new / worsening)")
            ),
            tags$hr(),
            tags$p(tags$b("High-risk features (hospitalise):")),
            tags$ul(
              tags$li("Fever >38°C"),
              tags$li("Subacute onset"),
              tags$li("Large effusion / tamponade"),
              tags$li("Myopericarditis"),
              tags$li("Immunosuppressed host")
            ),
            tags$hr(),
            tags$p(tags$b("Recurrence definition:")),
            tags$p("New episode ≥4-6 wks after symptom-free interval"),
            tags$p(tags$b("Incessant:"), " >4-6 wks without remission")
          )
        ),
        fluidRow(
          box(title="Simulation Overview", status="success", solidHeader=TRUE, width=12,
            plotlyOutput("overview_plot", height="400px")
          )
        )
      ),

      ## ---- TAB 2: Drug PK ----
      tabItem("tab_pk",
        fluidRow(
          box(title="Drug Concentration-Time Profiles", status="primary", solidHeader=TRUE, width=12,
            plotlyOutput("pk_plot", height="500px")
          )
        ),
        fluidRow(
          box(title="Pharmacodynamic Inhibition (%)", status="info", solidHeader=TRUE, width=6,
            plotlyOutput("pd_inh_plot", height="350px")
          ),
          box(title="PK Parameter Summary", status="warning", solidHeader=TRUE, width=6,
            tableOutput("pk_table")
          )
        )
      ),

      ## ---- TAB 3: Inflammasome / Cytokines ----
      tabItem("tab_cytokine",
        fluidRow(
          box(title="NLRP3 Inflammasome Activity (%)", status="danger", solidHeader=TRUE, width=6,
            plotlyOutput("nlrp3_plot", height="300px")
          ),
          box(title="IL-1β Dynamics (pg/mL)", status="danger", solidHeader=TRUE, width=6,
            plotlyOutput("il1b_plot", height="300px")
          )
        ),
        fluidRow(
          box(title="IL-18 (pg/mL)", status="warning", solidHeader=TRUE, width=4,
            plotlyOutput("il18_plot", height="280px")
          ),
          box(title="TNF-α (pg/mL)", status="warning", solidHeader=TRUE, width=4,
            plotlyOutput("tnf_plot", height="280px")
          ),
          box(title="IL-6 (pg/mL)", status="warning", solidHeader=TRUE, width=4,
            plotlyOutput("il6_plot", height="280px")
          )
        ),
        fluidRow(
          box(title="Pericardial Neutrophils & M1 Macrophages (rel.)", status="info", solidHeader=TRUE, width=12,
            plotlyOutput("immune_cells_plot", height="280px")
          )
        )
      ),

      ## ---- TAB 4: Pericardial Pathology ----
      tabItem("tab_pericardial",
        fluidRow(
          box(title="Pericardial Inflammation Score (0-10)", status="danger", solidHeader=TRUE, width=6,
            plotlyOutput("inflam_plot", height="300px")
          ),
          box(title="Pericardial Effusion Volume (mL)", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("effusion_plot", height="300px")
          )
        ),
        fluidRow(
          box(title="Fibrin Deposition (0-1)", status="info", solidHeader=TRUE, width=6,
            plotlyOutput("fibrin_plot", height="280px")
          ),
          box(title="Fibrosis Index (0-1)", status="primary", solidHeader=TRUE, width=6,
            plotlyOutput("fibrosis_plot", height="280px")
          )
        ),
        fluidRow(
          box(title="Risk Stratification", status="danger", solidHeader=TRUE, width=12,
            tags$table(class="table table-striped",
              tags$thead(tags$tr(tags$th("Effusion Volume"),tags$th("Echo Finding"),tags$th("Risk"),tags$th("Action"))),
              tags$tbody(
                tags$tr(tags$td("<100 mL"),tags$td("Small / trace"),tags$td("Low"),tags$td("Outpatient management")),
                tags$tr(tags$td("100-250 mL"),tags$td("Moderate rim"),tags$td("Moderate"),tags$td("Close follow-up")),
                tags$tr(tags$td(">250 mL"),tags$td("Large / swinging"),tags$td("High"),tags$td("Hospitalise")),
                tags$tr(tags$td("Any + tamponade signs"),tags$td("IVC plethora, RV collapse"),tags$td("Emergency"),tags$td("Pericardiocentesis"))
              )
            )
          )
        )
      ),

      ## ---- TAB 5: Clinical Endpoints ----
      tabItem("tab_endpoints",
        fluidRow(
          box(title="Chest Pain VAS (0-10)", status="danger", solidHeader=TRUE, width=6,
            plotlyOutput("pain_plot", height="300px")
          ),
          box(title="CRP (mg/L) — Recurrence Risk Marker", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("crp_plot", height="300px")
          )
        ),
        fluidRow(
          box(title="Key Time-Point Summary Table", status="info", solidHeader=TRUE, width=12,
            DTOutput("endpoint_table")
          )
        ),
        fluidRow(
          box(title="Treatment Goals (ESC 2015)", status="success", solidHeader=TRUE, width=12,
            tags$ul(
              tags$li(tags$b("Pain:"), " VAS ≤1 within 1-2 weeks"),
              tags$li(tags$b("CRP:"), " Normalize (≤3 mg/L) before tapering → CORP trial criterion"),
              tags$li(tags$b("Effusion:"), " Resolve or reduce on serial echo"),
              tags$li(tags$b("Colchicine duration:"), " ≥3 months for first episode; ≥6 months for recurrence"),
              tags$li(tags$b("Steroid taper:"), " Only after CRP normalization; slow taper (2-4 mg/2 wks)"),
              tags$li(tags$b("Recurrence target:"), " <15% at 18 months with colchicine+NSAID")
            )
          )
        )
      ),

      ## ---- TAB 6: Scenario Comparison ----
      tabItem("tab_compare",
        fluidRow(
          box(title="All Scenarios — Inflammation Score", status="danger", solidHeader=TRUE, width=6,
            plotlyOutput("cmp_inflam", height="300px")
          ),
          box(title="All Scenarios — CRP", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("cmp_crp", height="300px")
          )
        ),
        fluidRow(
          box(title="All Scenarios — IL-1β", status="info", solidHeader=TRUE, width=6,
            plotlyOutput("cmp_il1b", height="300px")
          ),
          box(title="All Scenarios — Effusion", status="primary", solidHeader=TRUE, width=6,
            plotlyOutput("cmp_effusion", height="300px")
          )
        ),
        fluidRow(
          box(title="Clinical Trial Benchmarks", status="success", solidHeader=TRUE, width=12,
            DTOutput("trial_table")
          )
        )
      ),

      ## ---- TAB 7: Biomarkers ----
      tabItem("tab_biomarkers",
        fluidRow(
          box(title="Biomarker Dashboard at Selected Timepoint", status="info", solidHeader=TRUE, width=4,
            sliderInput("bm_day", "Select Day", 0, 180, 7, 1),
            tableOutput("bm_snapshot")
          ),
          box(title="Biomarker Trajectories", status="primary", solidHeader=TRUE, width=8,
            plotlyOutput("bm_traj", height="400px")
          )
        ),
        fluidRow(
          box(title="NLRP3 Inhibition vs Colchicine Concentration (Emax)", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("emax_colch", height="300px")
          ),
          box(title="Pericardial Risk Summary", status="danger", solidHeader=TRUE, width=6,
            plotlyOutput("risk_radar", height="300px")
          )
        )
      )

    ) # end tabItems
  ) # end dashboardBody
) # end dashboardPage

# ---- SERVER ----
server <- function(input, output, session) {

  # Pre-compile all 7 scenarios on startup
  all_scenarios <- reactive({
    list(
      list(label="Untreated",   colch=FALSE, nsaid=FALSE, cs=FALSE, ana=FALSE, rilo=FALSE),
      list(label="NSAID",       colch=FALSE, nsaid=TRUE,  cs=FALSE, ana=FALSE, rilo=FALSE),
      list(label="Colchicine",  colch=TRUE,  nsaid=FALSE, cs=FALSE, ana=FALSE, rilo=FALSE),
      list(label="Colch+NSAID", colch=TRUE,  nsaid=TRUE,  cs=FALSE, ana=FALSE, rilo=FALSE),
      list(label="Prednisone",  colch=FALSE, nsaid=FALSE, cs=TRUE,  ana=FALSE, rilo=FALSE),
      list(label="Anakinra",    colch=FALSE, nsaid=FALSE, cs=FALSE, ana=TRUE,  rilo=FALSE),
      list(label="Rilonacept",  colch=FALSE, nsaid=FALSE, cs=FALSE, ana=FALSE, rilo=TRUE)
    )
  })

  sim_all <- reactive({
    purrr::map_dfr(all_scenarios(), function(sc) {
      e <- build_dose(sc$colch, 0.5, sc$nsaid, 600, sc$cs, 35, sc$ana, sc$rilo)
      mod %>%
        param(DISEASE_ON=1, TRIGGER_STRENGTH=input$trigger) %>%
        mrgsim(events=e, end=4320, delta=2) %>%
        as_tibble() %>%
        mutate(scenario=sc$label, day=time/24)
    })
  })

  # User custom scenario
  sim_user <- eventReactive(input$run_sim, {
    e <- build_dose(input$use_colch, input$dose_colch,
                    input$use_nsaid, input$dose_nsaid,
                    input$use_cs,   input$dose_cs,
                    input$use_ana,  input$use_rilo)
    mod %>%
      param(DISEASE_ON=1, TRIGGER_STRENGTH=input$trigger) %>%
      mrgsim(events=e, end=4320, delta=2) %>%
      as_tibble() %>%
      mutate(scenario="Custom", day=time/24)
  })

  make_plotly <- function(data, var, ylab, col_map=NULL) {
    d <- data %>% select(day, scenario, value=all_of(var))
    p <- ggplot(d, aes(day, value, colour=scenario)) +
      geom_line(size=0.9) +
      labs(x="Day", y=ylab, colour="Treatment") +
      theme_classic(base_size=11) +
      theme(legend.position="bottom")
    if (!is.null(col_map)) p <- p + scale_colour_manual(values=col_map)
    ggplotly(p) %>% layout(legend=list(orientation="h"))
  }

  output$overview_plot <- renderPlotly({
    df <- sim_user() %>% select(day, INFLAM, CRP, PAIN, EFFUSION) %>%
      pivot_longer(-day)
    ggplotly(
      ggplot(df, aes(day, value, colour=name)) +
        geom_line(size=1) +
        facet_wrap(~name, scales="free_y") +
        labs(x="Day", y="Value", title="Custom Scenario — Key Outputs") +
        theme_classic(base_size=11) +
        theme(legend.position="none")
    )
  })

  output$pk_plot <- renderPlotly({
    d <- sim_user() %>%
      select(day, `Colchicine (ng/mL)`=C_colch, `Ibuprofen (μM)`=C_nsaid,
             `Prednisolone (μM)`=C_cs, `Anakinra (nM)`=C_ana, `Rilonacept (nM)`=C_rilo) %>%
      pivot_longer(-day)
    ggplotly(
      ggplot(d %>% filter(value>0), aes(day, value, colour=name)) +
        geom_line(size=0.9) +
        facet_wrap(~name, scales="free_y") +
        labs(x="Day", y="Concentration") +
        theme_classic(base_size=11) +
        theme(legend.position="none")
    )
  })

  output$pd_inh_plot <- renderPlotly({
    d <- sim_user() %>%
      mutate(
        `Colch→NLRP3 (%)` = INH_colch_nlrp3 * 100,
        `NSAID→COX-2 (%)` = INH_nsaid_cox2  * 100,
        `CS→NF-κB (%)`    = INH_cs_nfkb     * 100,
        `Bio→IL-1β (%)`   = IL1_BIO_INH     * 100
      ) %>%
      select(day, `Colch→NLRP3 (%)`, `NSAID→COX-2 (%)`, `CS→NF-κB (%)`, `Bio→IL-1β (%)`) %>%
      pivot_longer(-day)
    ggplotly(
      ggplot(d %>% filter(value>0.5), aes(day, value, colour=name)) +
        geom_line(size=0.9) +
        ylim(0,100) +
        labs(x="Day", y="Inhibition (%)") +
        theme_classic(base_size=11) +
        theme(legend.position="bottom")
    )
  })

  output$pk_table <- renderTable({
    tibble(
      Drug         = c("Colchicine","Ibuprofen","Prednisone","Anakinra","Rilonacept"),
      `Dose`       = c("0.5 mg BID","600 mg TID","0.5 mg/kg/d","100 mg/d SC","320→160 mg qw"),
      `Vd (L)`     = c(450,15,50,6,4),
      `CL (L/h)`   = c(18,2.5,8,1.2,0.015),
      `t½ (h)`     = c(27,6,3,5,206),
      `Bioavail.`  = c("45%","87%","82%","95%","99%")
    )
  })

  output$nlrp3_plot  <- renderPlotly(make_plotly(sim_all(), "NLRP3_ACT", "Activity (0-1)", SCEN_COLS))
  output$il1b_plot   <- renderPlotly(make_plotly(sim_all(), "IL1B",      "IL-1β (pg/mL)", SCEN_COLS))
  output$il18_plot   <- renderPlotly(make_plotly(sim_all(), "IL18",      "IL-18 (pg/mL)", SCEN_COLS))
  output$tnf_plot    <- renderPlotly(make_plotly(sim_all(), "TNF",       "TNF-α (pg/mL)", SCEN_COLS))
  output$il6_plot    <- renderPlotly(make_plotly(sim_all(), "IL6",       "IL-6 (pg/mL)",  SCEN_COLS))
  output$inflam_plot <- renderPlotly(make_plotly(sim_all(), "INFLAM",    "Score (0-10)",   SCEN_COLS))
  output$effusion_plot <- renderPlotly(make_plotly(sim_all(), "EFFUSION","Effusion (mL)",  SCEN_COLS))
  output$fibrin_plot   <- renderPlotly(make_plotly(sim_all(), "FIBRIN",  "Fibrin (0-1)",   SCEN_COLS))
  output$fibrosis_plot <- renderPlotly(make_plotly(sim_all(), "FIBROSIS","Fibrosis (0-1)", SCEN_COLS))
  output$pain_plot     <- renderPlotly(make_plotly(sim_all(), "PAIN",    "VAS (0-10)",     SCEN_COLS))
  output$crp_plot      <- renderPlotly(make_plotly(sim_all(), "CRP",     "CRP (mg/L)",     SCEN_COLS))

  output$immune_cells_plot <- renderPlotly({
    d <- sim_all() %>%
      select(day, scenario, Neutrophils=NEUTRO, M1_Macro=M1_MACRO) %>%
      pivot_longer(c(Neutrophils, M1_Macro))
    ggplotly(
      ggplot(d, aes(day, value, colour=scenario, linetype=name)) +
        geom_line(size=0.9) +
        labs(x="Day", y="Relative Units", colour="Scenario", linetype="Cell") +
        theme_classic(base_size=11) +
        theme(legend.position="bottom")
    )
  })

  output$endpoint_table <- renderDT({
    sim_all() %>%
      filter(day %in% c(0, 3, 7, 14, 30, 90, 180)) %>%
      group_by(scenario, day) %>%
      slice(1) %>%
      mutate(across(c(IL1B, CRP, INFLAM, EFFUSION, PAIN, FIBROSIS), round, 2)) %>%
      select(Treatment=scenario, Day=day, `IL-1β`=IL1B, CRP, Inflammation=INFLAM,
             `Effusion(mL)`=EFFUSION, `Pain VAS`=PAIN, Fibrosis=FIBROSIS) %>%
      arrange(Day, Treatment)
  }, options=list(pageLength=20, scrollX=TRUE))

  output$cmp_inflam   <- renderPlotly(make_plotly(sim_all(), "INFLAM",   "Score (0-10)", SCEN_COLS))
  output$cmp_crp      <- renderPlotly(make_plotly(sim_all(), "CRP",      "CRP (mg/L)",  SCEN_COLS))
  output$cmp_il1b     <- renderPlotly(make_plotly(sim_all(), "IL1B",     "IL-1β (pg/mL)", SCEN_COLS))
  output$cmp_effusion <- renderPlotly(make_plotly(sim_all(), "EFFUSION", "Effusion (mL)", SCEN_COLS))

  output$trial_table <- renderDT({
    tibble(
      Trial   = c("COPE","ICAP","CORP","CORP-2","AIRTRIP","RHAPSODY"),
      Year    = c(2005,2013,2011,2014,2016,2021),
      Drug    = c("Colchicine+ASA","Colchicine 0.5mg BID","Colchicine 2nd ep","Colchicine incessant","Anakinra 100mg/d","Rilonacept 320→160mg"),
      n       = c(120,240,96,100,22,86),
      `Ctrl %` = c(45,32.3,45.5,50,90.9,74.4),
      `Trt %` = c(24,16.7,19.2,26.0,18.2,8.8),
      RRR     = c("47%","48%","58%","48%","80%","88%"),
      Journal = c("Lancet","NEJM","Arch Int Med","Ann Int Med","NEJM","NEJM")
    )
  })

  output$bm_snapshot <- renderTable({
    d <- sim_user() %>%
      filter(abs(day - input$bm_day) < 0.1) %>%
      slice(1) %>%
      select(IL1B, IL18, TNF, IL6, CRP, PAIN, INFLAM, EFFUSION, FIBRIN, FIBROSIS) %>%
      pivot_longer(everything(), names_to="Biomarker", values_to="Value") %>%
      mutate(Value=round(Value,2))
    d
  })

  output$bm_traj <- renderPlotly({
    d <- sim_user() %>%
      select(day, IL1B, CRP, PAIN, EFFUSION, FIBROSIS) %>%
      pivot_longer(-day)
    ggplotly(
      ggplot(d, aes(day, value, colour=name)) +
        geom_line(size=0.9) +
        facet_wrap(~name, scales="free_y") +
        labs(x="Day", y="Value") +
        theme_classic(base_size=11) +
        theme(legend.position="none")
    )
  })

  output$emax_colch <- renderPlotly({
    conc <- seq(0, 5, 0.01)
    inh  <- conc / (0.5 + conc)
    ggplotly(
      ggplot(tibble(conc, inh), aes(conc, inh*100)) +
        geom_line(colour="steelblue", size=1.2) +
        geom_vline(xintercept=0.5, linetype="dashed", colour="red") +
        geom_hline(yintercept=50, linetype="dashed", colour="red") +
        annotate("text", x=0.55, y=5, label="IC50 = 0.5 ng/mL", hjust=0, size=3) +
        labs(x="Colchicine Concentration (ng/mL)", y="NLRP3 Inhibition (%)") +
        theme_classic(base_size=11)
    )
  })

  output$risk_radar <- renderPlotly({
    d <- sim_user() %>% filter(day==7) %>% slice(1)
    tibble(
      Category = c("IL-1β","CRP","Inflammation","Effusion","Pain","Fibrosis Risk"),
      Score    = c(
        pmin(d$IL1B/300, 1)*10,
        pmin(d$CRP/30, 1)*10,
        d$INFLAM,
        pmin(d$EFFUSION/500, 1)*10,
        d$PAIN,
        d$FIBROSIS*10
      )
    ) %>%
    plot_ly(type="scatterpolar", r=~Score, theta=~Category, fill="toself") %>%
    layout(title=paste("Risk Radar at Day 7 (custom scenario)"),
           polar=list(radialaxis=list(visible=TRUE, range=c(0,10))))
  })

}

shinyApp(ui, server)
