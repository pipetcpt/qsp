## ============================================================
## Marfan Syndrome QSP — Interactive Shiny Dashboard
## ============================================================
## Tabs:
##   1. Patient Profile & Disease Overview
##   2. Drug Pharmacokinetics
##   3. TGF-β / Molecular PD
##   4. Cardiovascular Endpoints
##   5. Treatment Scenario Comparison
##   6. Biomarkers & Monitoring
##   7. Surgical Decision Support
## ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)
library(mrgsolve)

# ─────────────────────────────────────────────────────────────
# Inline mrgsolve model (compact version for Shiny)
# ─────────────────────────────────────────────────────────────
mfs_code <- '
$PARAM
KA_ATN=1.06 F_ATN=0.50 CL_ATN=10.8 Vc_ATN=67.0 Vp_ATN=64.0 Q_ATN=15.0
KA_LOS=1.41 F_LOS=0.33 CL_LOS=55.0 Vc_LOS=34.0 FM_LOS=0.14 CL_EXP=18.0 Vc_EXP=12.0
kprod_TGFb=0.05 kdeg_TGFb=0.10 TGFb0=0.50 TGFb_MFS_mult=2.5
kact_SMAD=0.30 kdact_SMAD=0.08 kact_ERK=0.20 kdact_ERK=0.12
kprod_MMP=0.04 kdeg_MMP=0.15 MMP0=0.27
k_ao_growth=1.8e-4 Emax_MMP_ao=0.80 EC50_MMP_ao=0.30
Emax_SMAD_ao=0.60 EC50_SMAD_ao=0.50
SBP0=125 DBP0=75 HR0=78 dPdt0=1200
k_AR_prog=4e-5 AR0=0.20
IC50_ATN_HR=65 Emax_ATN_HR=0.28
IC50_ATN_dPdt=55 Emax_ATN_dPdt=0.30
IC50_LOS_TGFb=40 Emax_LOS_TGFb=0.55
IC50_LOS_SBP=45 Emax_LOS_SBP=0.18
WT=68 AGE=30 MFS_FBN1=1

$CMT DEPOT_ATN C1_ATN C2_ATN DEPOT_LOS C1_LOS C_EXP3174
     TGFb pSMAD pERK MMP Ao_Diam AR_Grade HR SBP dPdt
     NT_proBNP LVEDD TGFb_plasma_obs Systemic_score

$INIT DEPOT_ATN=0 C1_ATN=0 C2_ATN=0 DEPOT_LOS=0 C1_LOS=0 C_EXP3174=0
      TGFb=1.25 pSMAD=1.0 pERK=1.0 MMP=0.27 Ao_Diam=30.0 AR_Grade=0.20
      HR=78 SBP=125 dPdt=1200 NT_proBNP=50 LVEDD=46 TGFb_plasma_obs=1.25
      Systemic_score=3.0

$ODE
dxdt_DEPOT_ATN = -KA_ATN*DEPOT_ATN;
dxdt_C1_ATN    =  KA_ATN*DEPOT_ATN*F_ATN - (CL_ATN+Q_ATN)/Vc_ATN*C1_ATN + Q_ATN/Vp_ATN*C2_ATN;
dxdt_C2_ATN    =  Q_ATN/Vc_ATN*C1_ATN - Q_ATN/Vp_ATN*C2_ATN;
dxdt_DEPOT_LOS = -KA_LOS*DEPOT_LOS;
dxdt_C1_LOS    =  KA_LOS*DEPOT_LOS*F_LOS - CL_LOS/Vc_LOS*C1_LOS;
dxdt_C_EXP3174 =  FM_LOS*CL_LOS/Vc_LOS*C1_LOS - CL_EXP/Vc_EXP*C_EXP3174;
double TGFb_ss = TGFb0*(1+(TGFb_MFS_mult-1)*MFS_FBN1);
double Cp_EXP  = C_EXP3174/Vc_EXP*1000.0;
double E_LOS   = Emax_LOS_TGFb*Cp_EXP/(IC50_LOS_TGFb+Cp_EXP);
dxdt_TGFb      = kprod_TGFb*TGFb_ss*(1-E_LOS) - kdeg_TGFb*TGFb;
double TGFb_norm = TGFb/TGFb_ss;
dxdt_pSMAD     = kact_SMAD*TGFb_norm*(1-E_LOS) - kdact_SMAD*pSMAD;
dxdt_pERK      = kact_ERK*TGFb_norm*(1-E_LOS*0.6) - kdact_ERK*pERK;
double kp_MMP  = kprod_MMP*(0.5*pSMAD+0.5*pERK);
dxdt_MMP       = kp_MMP - kdeg_MMP*MMP;
double Cp_ATN  = C1_ATN/Vc_ATN*1000.0;
double E_dPdt  = Emax_ATN_dPdt*Cp_ATN/(IC50_ATN_dPdt+Cp_ATN);
double hstress = (dPdt/dPdt0)*(1-E_dPdt);
double E_MMPd  = Emax_MMP_ao*MMP/(EC50_MMP_ao+MMP);
double E_SMAD  = Emax_SMAD_ao*pSMAD/(EC50_SMAD_ao+pSMAD);
dxdt_Ao_Diam   = k_ao_growth*hstress*(1+E_MMPd+E_SMAD)*Ao_Diam;
double AR_d    = (Ao_Diam>37.0)?(Ao_Diam-37.0)/20.0:0.0;
dxdt_AR_Grade  = k_AR_prog*AR_d*(4.0-AR_Grade);
double E_HR    = Emax_ATN_HR*Cp_ATN/(IC50_ATN_HR+Cp_ATN);
double E_SBP   = Emax_LOS_SBP*Cp_EXP/(IC50_LOS_SBP+Cp_EXP);
dxdt_HR    = 0.5*(HR0*(1-E_HR) - HR);
dxdt_SBP   = 0.5*(SBP0*(1-E_SBP) - SBP);
dxdt_dPdt  = 0.5*(dPdt0*(1-E_dPdt) - dPdt);
dxdt_NT_proBNP = 0.02*(50*(1+AR_Grade*0.5) - NT_proBNP);
dxdt_LVEDD = 0.005*(AR_Grade/4.0)*(1+(Ao_Diam-30)/30.0);
dxdt_TGFb_plasma_obs = 0.5*(TGFb - TGFb_plasma_obs);
double tsc = 0.0 + 2.0 + 3.0 + 2.0 + ((Ao_Diam>=50)?2.0:(Ao_Diam>=42)?1.0:0.0);
dxdt_Systemic_score  = 0.001*(tsc - Systemic_score);

$TABLE
capture CpATN = C1_ATN/Vc_ATN*1000;
capture CpLOS = C1_LOS/Vc_LOS*1000;
capture CpEXP = C_EXP3174/Vc_EXP*1000;
capture TGFb_c = TGFb; capture pSMAD_f=pSMAD; capture pERK_f=pERK;
capture MMP_a=MMP; capture AoD=Ao_Diam; capture AR_g=AR_Grade;
capture HR_v=HR; capture SBP_v=SBP; capture dPdt_v=dPdt;
capture BNP=NT_proBNP; capture LVD=LVEDD; capture GS=Systemic_score;
capture Zscore=(Ao_Diam-23.0)/3.5;
capture growth_rate_yr = k_ao_growth*(dPdt/dPdt0)*
  (1-(Emax_ATN_dPdt*(C1_ATN/Vc_ATN*1000)/(IC50_ATN_dPdt+C1_ATN/Vc_ATN*1000)))*
  (1+(Emax_MMP_ao*MMP/(EC50_MMP_ao+MMP))+(Emax_SMAD_ao*pSMAD/(EC50_SMAD_ao+pSMAD)))*
  Ao_Diam*8760;
$CAPTURE CpATN CpLOS CpEXP TGFb_c pSMAD_f pERK_f MMP_a AoD AR_g HR_v SBP_v dPdt_v BNP LVD GS Zscore growth_rate_yr
'

mod <- mcode("mfs_shiny", mfs_code)

run_sim <- function(mod, dose_atn, dose_los, dur_yr,
                    aod0 = 30, tgfbm = 2.5, ageyr = 30) {
  end_h <- dur_yr * 8760
  ev_atn <- data.frame(time = seq(0, end_h - 24, by = 24),
                       amt = dose_atn, cmt = 1, evid = 1)
  ev_los <- data.frame(time = seq(0, end_h - 24, by = 24),
                       amt = dose_los, cmt = 4, evid = 1)
  ev_all <- bind_rows(ev_atn, ev_los) %>% arrange(time)
  mod %>%
    param(TGFb_MFS_mult = tgfbm, AGE = ageyr) %>%
    init(Ao_Diam = aod0, TGFb = 0.5 * tgfbm, pSMAD = 1, pERK = 1) %>%
    mrgsim_df(data = ev_all, end = end_h, delta = 24) %>%
    mutate(time_yr = time / 8760)
}

# ─────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "MFS QSP Dashboard", titleWidth = 280),
  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("① Patient Profile",    tabName = "patient",   icon = icon("user-md")),
      menuItem("② Drug PK",            tabName = "pk",        icon = icon("pills")),
      menuItem("③ TGF-β / Mol. PD",   tabName = "tgfb",      icon = icon("dna")),
      menuItem("④ Cardiovascular",     tabName = "cardio",    icon = icon("heartbeat")),
      menuItem("⑤ Scenario Comparison",tabName = "scenarios", icon = icon("chart-line")),
      menuItem("⑥ Biomarkers",         tabName = "biomarkers",icon = icon("vial")),
      menuItem("⑦ Surgical Decision",  tabName = "surgery",   icon = icon("scalpel"))
    ),
    hr(),
    h5("Patient Parameters", style = "padding-left:15px; color:#AAA"),
    sliderInput("aod0",   "Baseline Ao Root (mm)", 25, 45, 30, 1),
    sliderInput("ageyr",  "Age (years)",            5,  60, 30, 1),
    sliderInput("tgfbm",  "TGF-β multiplier",       1.0, 4.0, 2.5, 0.1),
    sliderInput("dur_yr", "Duration (years)",        1, 10, 5, 1),
    hr(),
    h5("Drug Doses", style = "padding-left:15px; color:#AAA"),
    sliderInput("dose_atn", "Atenolol (mg/day)",   0, 200, 50, 25),
    sliderInput("dose_los", "Losartan (mg/day)",   0, 200, 50, 25),
    actionButton("run", "▶ Run Simulation", class = "btn-primary btn-block",
                 style = "margin:10px")
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .box-header { background-color: #1565C0 !important; color: white !important; }
      .small-box { border-radius:8px; }
      body { font-family: 'Helvetica Neue', sans-serif; }
    "))),
    tabItems(
      # ── Tab 1: Patient Profile ─────────────────────────────
      tabItem(tabName = "patient",
        fluidRow(
          valueBoxOutput("vb_aod",  width = 3),
          valueBoxOutput("vb_zscore", width = 3),
          valueBoxOutput("vb_ghent", width = 3),
          valueBoxOutput("vb_ar",  width = 3)
        ),
        fluidRow(
          box(title = "Marfan Syndrome Pathophysiology Overview", width = 8, status = "primary",
            p(strong("Marfan syndrome (MFS)"), "is an autosomal dominant connective tissue disorder caused by mutations in the",
              strong("FBN1 gene"), "(Chr 15q21.1), which encodes fibrillin-1."),
            p("Defective fibrillin-1 disrupts extracellular microfibril networks, leading to:"),
            tags$ul(
              tags$li(strong("Impaired TGF-β sequestration"), "→ excess free TGF-β1/2 → SMAD2/3 + ERK/MAPK hyperactivation"),
              tags$li("ECM degradation via MMP-2/9 upregulation → aortic media degeneration"),
              tags$li("Progressive", strong("aortic root dilation"), "(Valsalva sinus) → AR, dissection risk"),
              tags$li("Skeletal features: tall stature, arachnodactyly, scoliosis, pectus"),
              tags$li("Ocular: ectopia lentis (60-70%), myopia, retinal detachment risk")
            ),
            p(strong("Treatment goals:"), "reduce aortic root growth rate to < 0.5 mm/year; prevent dissection; maintain diastolic BP 60-70 mmHg, HR < 70 bpm.")
          ),
          box(title = "Revised Ghent Criteria (2010)", width = 4, status = "info",
            tableOutput("ghent_table")
          )
        ),
        fluidRow(
          box(title = "Mechanistic Map", width = 12, status = "success",
            tags$img(src = "mfs_qsp_model.png", width = "100%",
                     onerror = "this.style.display='none'; this.nextSibling.style.display='block'"),
            tags$p("(PNG preview — open mfs_qsp_model.svg for full resolution)",
                   style = "display:none; color:#888")
          )
        )
      ),

      # ── Tab 2: Drug PK ─────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Atenolol PK — Central Concentration (ng/mL)", width = 6, status = "primary",
              plotOutput("pk_atn_plot", height = 300)),
          box(title = "Losartan & EXP-3174 PK (ng/mL)", width = 6, status = "success",
              plotOutput("pk_los_plot", height = 300))
        ),
        fluidRow(
          box(title = "PK Summary at Steady State", width = 12, status = "info",
              DTOutput("pk_summary"))
        ),
        fluidRow(
          box(title = "PK Parameter Reference", width = 12, status = "warning",
            tableOutput("pk_params_table")
          )
        )
      ),

      # ── Tab 3: TGF-β / Molecular PD ────────────────────────
      tabItem(tabName = "tgfb",
        fluidRow(
          box(title = "Plasma TGF-β1 Dynamics (ng/mL)", width = 6, status = "primary",
              plotOutput("tgfb_plot", height = 300)),
          box(title = "p-SMAD2/3 Activity (fold vs baseline)", width = 6, status = "success",
              plotOutput("psmad_plot", height = 300))
        ),
        fluidRow(
          box(title = "p-ERK1/2 Activity (fold vs baseline)", width = 6, status = "warning",
              plotOutput("perk_plot", height = 300)),
          box(title = "MMP Activity (U/mL)", width = 6, status = "danger",
              plotOutput("mmp_plot", height = 300))
        )
      ),

      # ── Tab 4: Cardiovascular Endpoints ────────────────────
      tabItem(tabName = "cardio",
        fluidRow(
          box(title = "Aortic Root Diameter (mm)", width = 6, status = "danger",
              plotOutput("aod_plot", height = 320)),
          box(title = "Aortic Root Z-Score", width = 6, status = "warning",
              plotOutput("zscore_plot", height = 320))
        ),
        fluidRow(
          box(title = "Aortic Regurgitation Grade (0-4)", width = 6, status = "primary",
              plotOutput("ar_plot", height = 300)),
          box(title = "Heart Rate & dP/dt (β-blocker effects)", width = 6, status = "success",
              plotOutput("hemo_plot", height = 300))
        ),
        fluidRow(
          box(title = "LV End-Diastolic Diameter (mm)", width = 6, status = "info",
              plotOutput("lvedd_plot", height = 280)),
          box(title = "Systolic Blood Pressure (mmHg)", width = 6, status = "info",
              plotOutput("sbp_plot", height = 280))
        )
      ),

      # ── Tab 5: Scenario Comparison ──────────────────────────
      tabItem(tabName = "scenarios",
        fluidRow(
          box(title = "All 6 Treatment Scenarios — Aortic Root Diameter", width = 12,
              status = "danger", plotOutput("scenario_aod", height = 380))
        ),
        fluidRow(
          box(title = "Summary at 5 Years", width = 12, status = "primary",
              DTOutput("scenario_table"))
        )
      ),

      # ── Tab 6: Biomarkers ───────────────────────────────────
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(title = "NT-proBNP (pg/mL)", width = 6, status = "warning",
              plotOutput("bnp_plot", height = 300)),
          box(title = "TGF-β1 Plasma (ng/mL)", width = 6, status = "info",
              plotOutput("biomarker_tgfb", height = 300))
        ),
        fluidRow(
          box(title = "Ghent Systemic Score (0-20)", width = 6, status = "success",
              plotOutput("ghent_plot", height = 280)),
          box(title = "Annual Aortic Growth Rate (mm/yr)", width = 6, status = "danger",
              plotOutput("growth_rate_plot", height = 280))
        ),
        fluidRow(
          box(title = "Monitoring Thresholds (clinical reference)", width = 12, status = "info",
              tableOutput("biomarker_ref"))
        )
      ),

      # ── Tab 7: Surgical Decision Support ───────────────────
      tabItem(tabName = "surgery",
        fluidRow(
          valueBoxOutput("vb_surg_now", width = 4),
          valueBoxOutput("vb_surg_yr", width = 4),
          valueBoxOutput("vb_dissect", width = 4)
        ),
        fluidRow(
          box(title = "Time to Surgical Threshold (50mm)", width = 6, status = "danger",
              plotOutput("surg_threshold_plot", height = 320)),
          box(title = "Surgical Decision Guidelines (AHA/ESC)", width = 6, status = "info",
              tableOutput("surg_guidelines"))
        ),
        fluidRow(
          box(title = "Treatment Effect on Time to Surgery", width = 12, status = "primary",
              plotOutput("treatment_effect_plot", height = 320))
        )
      )
    )
  )
)

# ─────────────────────────────────────────────────────────────
# SERVER
# ─────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Reactive simulation
  sim_data <- eventReactive(input$run, {
    withProgress(message = "Running QSP simulation...", {
      run_sim(mod,
              dose_atn = input$dose_atn,
              dose_los = input$dose_los,
              dur_yr   = input$dur_yr,
              aod0     = input$aod0,
              tgfbm    = input$tgfbm,
              ageyr    = input$ageyr)
    })
  }, ignoreNULL = FALSE)

  # Run all 6 scenarios (for scenario comparison tab)
  all6 <- eventReactive(input$run, {
    scen <- list(
      list(atn=0,   los=0,   label="1. Untreated"),
      list(atn=50,  los=0,   label="2. Atenolol 50mg"),
      list(atn=100, los=0,   label="3. Atenolol 100mg"),
      list(atn=0,   los=50,  label="4. Losartan 50mg"),
      list(atn=0,   los=100, label="5. Losartan 100mg"),
      list(atn=50,  los=50,  label="6. ATN+LOS combo")
    )
    lapply(scen, function(s) {
      run_sim(mod, s$atn, s$los, dur_yr = 5,
              aod0 = input$aod0, tgfbm = input$tgfbm, ageyr = input$ageyr) %>%
        mutate(scenario = s$label)
    }) %>% bind_rows() %>%
      mutate(scenario = factor(scenario, levels = sapply(scen, `[[`, "label")))
  }, ignoreNULL = FALSE)

  sc_colors <- c("1. Untreated"="#B71C1C","2. Atenolol 50mg"="#1565C0",
                 "3. Atenolol 100mg"="#0D47A1","4. Losartan 50mg"="#2E7D32",
                 "5. Losartan 100mg"="#1B5E20","6. ATN+LOS combo"="#4527A0")

  # ── Value boxes ──
  output$vb_aod <- renderValueBox({
    d <- tail(sim_data(), 1)
    valueBox(sprintf("%.1f mm", d$AoD), "Aortic Root Diameter", icon = icon("heart"), color = "red")
  })
  output$vb_zscore <- renderValueBox({
    d <- tail(sim_data(), 1)
    col <- if (d$Zscore > 3) "red" else if (d$Zscore > 2) "orange" else "green"
    valueBox(sprintf("%.2f", d$Zscore), "Z-Score", icon = icon("ruler"), color = col)
  })
  output$vb_ghent <- renderValueBox({
    d <- tail(sim_data(), 1)
    valueBox(sprintf("%.1f", d$GS), "Ghent Score", icon = icon("clipboard-check"), color = "blue")
  })
  output$vb_ar <- renderValueBox({
    d <- tail(sim_data(), 1)
    col <- if (d$AR_g > 2) "red" else if (d$AR_g > 1) "orange" else "green"
    valueBox(sprintf("%.2f", d$AR_g), "AR Grade", icon = icon("wave-square"), color = col)
  })

  output$vb_surg_now <- renderValueBox({
    d <- tail(sim_data(), 1)
    surg <- if (d$AoD >= 50) "NOW" else sprintf("%.1f mm to go", 50 - d$AoD)
    valueBox(surg, "Surgery Threshold (50mm)", icon = icon("cut"), color = "red")
  })
  output$vb_surg_yr <- renderValueBox({
    df <- sim_data()
    cross <- df %>% filter(AoD >= 50) %>% pull(time_yr)
    val <- if (length(cross) > 0) sprintf("Year %.1f", min(cross)) else "> sim horizon"
    valueBox(val, "Time to 50mm", icon = icon("clock"), color = "orange")
  })
  output$vb_dissect <- renderValueBox({
    d <- tail(sim_data(), 1)
    risk <- round((d$AoD / 50) ^ 4 * 2, 1)   # simplified risk estimate %/yr
    valueBox(sprintf("%.1f%%/yr", risk), "Estimated Dissection Risk", icon = icon("exclamation-triangle"), color = "red")
  })

  # Ghent table
  output$ghent_table <- renderTable({
    data.frame(
      Feature      = c("Aortic root Z ≥ 2", "Ectopia lentis", "FBN1 pathogenic variant",
                       "Systemic score ≥ 7", "Dural ectasia", "Wrist+thumb sign",
                       "Pectus carinatum", "Scoliosis ≥ 20°"),
      Points       = c(2, 2, 2, 1, 2, 3, 2, 1),
      Diagnostic   = c("Major", "Major", "Major", "Minor", "Minor", "Minor", "Minor", "Minor")
    )
  })

  # PK table
  output$pk_params_table <- renderTable({
    data.frame(
      Drug        = c("Atenolol","Atenolol","Losartan","EXP-3174"),
      Parameter   = c("Oral F","T½","Oral F","T½"),
      Value       = c("50%","6-7 h","33%","3-4 h"),
      Source      = c("Öhrvall 1994","Öhrvall 1994","McCrea 1996","McCrea 1996")
    )
  })

  # PK plots
  output$pk_atn_plot <- renderPlot({
    df <- sim_data() %>% filter(time_yr <= 0.5)
    ggplot(df, aes(time_yr * 365, CpATN)) +
      geom_line(colour = "#1565C0", linewidth = 1.2) +
      labs(x = "Day", y = "Atenolol [ng/mL]") + theme_bw()
  })
  output$pk_los_plot <- renderPlot({
    df <- sim_data() %>% filter(time_yr <= 0.5) %>%
      select(time_yr, CpLOS, CpEXP) %>%
      pivot_longer(c(CpLOS, CpEXP))
    ggplot(df, aes(time_yr * 365, value, colour = name)) +
      geom_line(linewidth = 1.2) +
      scale_colour_manual(values = c(CpLOS = "#2E7D32", CpEXP = "#4CAF50"),
                          labels = c("Losartan", "EXP-3174")) +
      labs(x = "Day", y = "Concentration [ng/mL]", colour = NULL) + theme_bw()
  })
  output$pk_summary <- renderDT({
    sim_data() %>% filter(time_yr > input$dur_yr - 0.1) %>%
      summarise(Atenolol_Ctrough_ng = round(mean(CpATN), 1),
                Losartan_Ctrough    = round(mean(CpLOS), 1),
                EXP3174_Ctrough     = round(mean(CpEXP), 1)) %>%
      datatable(options = list(dom = "t"), rownames = FALSE)
  })

  # TGF-β / molecular
  output$tgfb_plot <- renderPlot({
    ggplot(sim_data(), aes(time_yr, TGFb_c)) +
      geom_line(colour = "#2E7D32", linewidth = 1.2) +
      geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
      labs(x = "Years", y = "TGF-β1 [ng/mL]") + theme_bw()
  })
  output$psmad_plot <- renderPlot({
    ggplot(sim_data(), aes(time_yr, pSMAD_f)) +
      geom_line(colour = "#1B5E20", linewidth = 1.2) +
      geom_hline(yintercept = 1, linetype = "dotted") +
      labs(x = "Years", y = "p-SMAD2/3 (fold)") + theme_bw()
  })
  output$perk_plot <- renderPlot({
    ggplot(sim_data(), aes(time_yr, pERK_f)) +
      geom_line(colour = "#1565C0", linewidth = 1.2) +
      labs(x = "Years", y = "p-ERK1/2 (fold)") + theme_bw()
  })
  output$mmp_plot <- renderPlot({
    ggplot(sim_data(), aes(time_yr, MMP_a)) +
      geom_line(colour = "#880E4F", linewidth = 1.2) +
      labs(x = "Years", y = "MMP [U/mL]") + theme_bw()
  })

  # Cardiovascular
  output$aod_plot <- renderPlot({
    ggplot(sim_data(), aes(time_yr, AoD)) +
      geom_line(colour = "#B71C1C", linewidth = 1.3) +
      geom_hline(yintercept = 50, linetype = "dashed", colour = "red") +
      annotate("text", x = 0.1, y = 50.5, label = "Surgical threshold (50mm)",
               hjust = 0, colour = "red", size = 3.5) +
      labs(x = "Years", y = "Aortic Root Diameter (mm)") + theme_bw()
  })
  output$zscore_plot <- renderPlot({
    ggplot(sim_data(), aes(time_yr, Zscore)) +
      geom_line(colour = "#E65100", linewidth = 1.3) +
      geom_hline(yintercept = 2, linetype = "dashed", colour = "orange") +
      labs(x = "Years", y = "Aortic Root Z-Score") + theme_bw()
  })
  output$ar_plot <- renderPlot({
    ggplot(sim_data(), aes(time_yr, AR_g)) +
      geom_line(colour = "#1565C0", linewidth = 1.2) +
      scale_y_continuous(limits = c(0, 4)) +
      labs(x = "Years", y = "AR Grade (0-4)") + theme_bw()
  })
  output$hemo_plot <- renderPlot({
    sim_data() %>% select(time_yr, HR_v, dPdt_v) %>%
      pivot_longer(c(HR_v, dPdt_v)) %>%
      ggplot(aes(time_yr, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      facet_wrap(~name, scales = "free_y",
                 labeller = labeller(name = c(HR_v = "HR (bpm)", dPdt_v = "dP/dt (mmHg/s)"))) +
      scale_colour_manual(values = c(HR_v = "#1565C0", dPdt_v = "#880E4F"), guide = "none") +
      labs(x = "Years", y = "Value") + theme_bw()
  })
  output$lvedd_plot <- renderPlot({
    ggplot(sim_data(), aes(time_yr, LVD)) +
      geom_line(colour = "#4527A0", linewidth = 1.2) +
      geom_hline(yintercept = 56, linetype = "dashed", colour = "red") +
      labs(x = "Years", y = "LVEDD (mm)") + theme_bw()
  })
  output$sbp_plot <- renderPlot({
    ggplot(sim_data(), aes(time_yr, SBP_v)) +
      geom_line(colour = "#00695C", linewidth = 1.2) +
      geom_hline(yintercept = 120, linetype = "dashed", colour = "orange") +
      labs(x = "Years", y = "SBP (mmHg)") + theme_bw()
  })

  # Scenario comparison
  output$scenario_aod <- renderPlot({
    ggplot(all6(), aes(time_yr, AoD, colour = scenario)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = 50, linetype = "dashed", colour = "red") +
      scale_colour_manual(values = sc_colors) +
      labs(x = "Years", y = "Aortic Root (mm)", colour = "Scenario") +
      theme_bw() + theme(legend.position = "bottom")
  })
  output$scenario_table <- renderDT({
    all6() %>% filter(time_yr >= 4.9) %>%
      group_by(scenario) %>%
      summarise(AoD_5yr = round(mean(AoD), 1),
                Zscore_5yr = round(mean(Zscore), 2),
                AR_5yr = round(mean(AR_g), 2),
                TGFb_5yr = round(mean(TGFb_c), 3),
                HR_5yr = round(mean(HR_v), 1),
                SBP_5yr = round(mean(SBP_v), 1),
                .groups = "drop") %>%
      datatable(options = list(pageLength = 10), rownames = FALSE)
  })

  # Biomarkers
  output$bnp_plot <- renderPlot({
    ggplot(sim_data(), aes(time_yr, BNP)) +
      geom_line(colour = "#E65100", linewidth = 1.2) +
      labs(x = "Years", y = "NT-proBNP (pg/mL)") + theme_bw()
  })
  output$biomarker_tgfb <- renderPlot({
    ggplot(sim_data(), aes(time_yr, TGFb_c)) +
      geom_line(colour = "#2E7D32", linewidth = 1.2) +
      labs(x = "Years", y = "TGF-β1 [ng/mL]") + theme_bw()
  })
  output$ghent_plot <- renderPlot({
    ggplot(sim_data(), aes(time_yr, GS)) +
      geom_line(colour = "#1565C0", linewidth = 1.2) +
      geom_hline(yintercept = 7, linetype = "dashed", colour = "orange") +
      labs(x = "Years", y = "Ghent Systemic Score") + theme_bw()
  })
  output$growth_rate_plot <- renderPlot({
    ggplot(sim_data(), aes(time_yr, growth_rate_yr)) +
      geom_line(colour = "#880E4F", linewidth = 1.2) +
      geom_hline(yintercept = 0.5, linetype = "dashed", colour = "red") +
      annotate("text", x = 0.1, y = 0.55, label = "Target < 0.5 mm/yr",
               hjust = 0, colour = "red", size = 3.5) +
      labs(x = "Years", y = "Annual Growth Rate (mm/yr)") + theme_bw()
  })
  output$biomarker_ref <- renderTable({
    data.frame(
      Biomarker     = c("Plasma TGF-β1","MMP-9","NT-proBNP","Aortic root Z-score","Aortic root [mm]"),
      Normal        = c("< 1 ng/mL","< 0.25 U/mL","< 125 pg/mL","< 2","< 40 mm"),
      MFS_Elevated  = c("2-4× ↑","↑","↑ with AR/LV dilation","≥ 2","varies"),
      Clinical_Use  = c("TGF-β pathway activity","ECM remodelling","Cardiac stress","Paediatric diagnosis","Adult threshold (surgical consideration)")
    )
  })

  # Surgical decision
  output$surg_threshold_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(time_yr, AoD)) +
      geom_line(colour = "#B71C1C", linewidth = 1.4) +
      geom_hline(yintercept = 50, linetype = "solid",  colour = "red",   linewidth = 1) +
      geom_hline(yintercept = 45, linetype = "dashed", colour = "orange",linewidth = 0.8) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 50, ymax = Inf,
               alpha = 0.1, fill = "red") +
      annotate("text", x = 0.1, y = 50.8, label = "Class I: ≥ 50mm (elective repair)", hjust = 0, size = 3.5) +
      annotate("text", x = 0.1, y = 45.8, label = "Class IIa: ≥ 45mm + risk factors",  hjust = 0, size = 3.5) +
      labs(x = "Years", y = "Aortic Root (mm)",
           title = "Aortic Root Growth vs Surgical Thresholds") + theme_bw()
  })
  output$surg_guidelines <- renderTable({
    data.frame(
      Class = c("I","IIa","IIa","IIb"),
      Indication = c(
        "Ao root ≥ 50mm → elective repair",
        "Ao root 45-50mm + family hx dissection",
        "Ao root 45-50mm + rapid growth (≥ 3mm/yr)",
        "Ao root 40-45mm + severe AR"
      ),
      Procedure = c("Bentall / David","Bentall / David","Bentall / David","Valve-sparing / Bentall"),
      Source = c("ESC 2024","ESC 2024","AHA 2022","AHA 2022")
    )
  })
  output$treatment_effect_plot <- renderPlot({
    all6() %>%
      ggplot(aes(time_yr, AoD, colour = scenario)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = 50, linetype = "dashed", colour = "red") +
      geom_hline(yintercept = 45, linetype = "dotted", colour = "orange") +
      scale_colour_manual(values = sc_colors) +
      labs(x = "Years", y = "Aortic Root Diameter (mm)",
           colour = "Treatment",
           title  = "Treatment Effect on Time to Surgical Threshold") +
      theme_bw() + theme(legend.position = "bottom")
  })
}

# ─────────────────────────────────────────────────────────────
# LAUNCH
# ─────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
