# =============================================================================
# Vitiligo QSP — Interactive Shiny Dashboard
# =============================================================================
# Requires: shiny, bslib, plotly, DT, mrgsolve, dplyr, tidyr, ggplot2
# Run:      shiny::runApp("vitiligo/vit_shiny_app.R")
# =============================================================================

library(shiny)
library(bslib)
library(plotly)
library(DT)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---- Inline mrgsolve model --------------------------------------------------
vit_code <- '
$PARAM @annotated
ka_ruxo:1.50:Ruxo oral ka(/h)  Vc_ruxo:72.0:Vc oral(L)
CL_ruxo:17.7:CL oral(L/h)     F_ruxo:0.95:Oral F
ka_sk:0.120:Topical ka(/h)     ke_sk:0.028:Skin ke(/h)
ka_afam:0.018:Afam depot ka    CL_afam:0.350:Afam CL
Vc_afam:8.00:Afam Vc           JAK1_IC50:3.30:JAK1 IC50(nM)
n_hill:1.20:Hill coef          Emax_JAK:0.95:Max JAK inh
k_mel_birth:3.0e-4:Mel birth   k_mel_death:3.0e-4:Mel death
k_mel_kill:0.085:CD8 kill rate k_NKGL_up:0.055:NKGL up
k_NKGL_down:0.040:NKGL down    k_CD8_rec:0.0022:CD8 recruit
k_CD8_death:0.042:CD8 death    k_CD8_prol:0.055:CD8 prolif
k_Treg_in:0.008:Treg in        k_Treg_out:0.025:Treg out
k_IFNG_prod:2.20:IFNg prod     k_IFNG_deg:0.28:IFNg deg
Treg_sup_IFNG:0.35:Treg sup    IFNG0:12.0:IFNg baseline
k_CXCL10_prod:14.0:CXCL10 syn k_CXCL10_deg:0.080:CXCL10 deg
CXCL10_0:82.0:CXCL10 baseline  k_pST1_act:0.55:pSTAT1 act
k_pST1_deg:0.22:pSTAT1 deg     pSTAT1_0:0.35:pSTAT1 baseline
k_MITF_base:0.040:MITF base    k_MITF_cAMP:0.38:MITF cAMP
k_MITF_deg:0.055:MITF deg      MITF0:0.72:MITF baseline
k_MEL_prod:0.0038:Melanin syn  k_MEL_loss:0.0008:Melanin loss
MELANIN0:0.40:Melanin baseline  k_FHAIR_mob:0.0015:Follicle mob
k_FHAIR_ren:1.5e-4:Follicle ren HAIRFOL0:0.85:Follicle base
VASI0:48.0:Baseline VASI        k_VASI_prog:3.5e-4:VASI worsening
k_VASI_rep:5.0e-4:VASI improv  VASI_floor:0.5:VASI floor
k_repig_rate:5.0e-4:Repig rate  k_inflam_in:0.12:Inflam in
k_inflam_out:0.10:Inflam out    NBUVB_Treg_stim:0.012:NBUVB Treg
NBUVB_mel_stim:0.005:NBUVB mel  ruxo_cream_BID:0:Cream BID flag
ruxo_cream_QD:0:Cream QD flag   ruxo_oral_BID:0:Oral BID flag
afam_on:0:Afam flag             NBUVB_on:0:NBUVB flag
CD8E0:0.30:CD8 baseline         TREG0:0.32:Treg baseline
NKGL0:0.20:NKGL baseline        NKGD0:0.18:NKGD baseline

$CMT
RUXO_GUT RUXO_C RUXO_SK AFAM_D AFAM_C
MEL NKGL CD8E TREG IFNG CXCL10 PSTAT1 MITF_C MELANIN
HAIRFOL NKGD_ACT TREG_SKIN INFLAM VASI REPIG

$MAIN
RUXO_GUT_0=0; RUXO_C_0=0; RUXO_SK_0=0;
AFAM_D_0=0; AFAM_C_0=0;
MEL_0=1.0; NKGL_0=NKGL0; CD8E_0=CD8E0; TREG_0=TREG0;
IFNG_0=IFNG0; CXCL10_0=CXCL10_0; PSTAT1_0=pSTAT1_0;
MITF_C_0=MITF0; MELANIN_0=MELANIN0; HAIRFOL_0=HAIRFOL0;
NKGD_ACT_0=NKGD0; TREG_SKIN_0=TREG0*0.5;
INFLAM_0=(IFNG0/15.0+CXCL10_0/100.0)*0.5;
VASI_0=VASI0; REPIG_0=0;

$ODE
double ke_ruxo=CL_ruxo/Vc_ruxo;
dxdt_RUXO_GUT=-ka_ruxo*RUXO_GUT;
dxdt_RUXO_C=ka_ruxo*F_ruxo*RUXO_GUT/Vc_ruxo-ke_ruxo*RUXO_C;
dxdt_RUXO_SK=-ke_sk*RUXO_SK;
double ke_afam=CL_afam/Vc_afam;
dxdt_AFAM_D=-ka_afam*AFAM_D;
dxdt_AFAM_C=ka_afam*AFAM_D/Vc_afam-ke_afam*AFAM_C;
double ruxo_oral_nM=RUXO_C*ruxo_oral_BID*1000.0/306.4;
double ruxo_skin_nM=RUXO_SK*(ruxo_cream_BID+ruxo_cream_QD*0.6);
double ruxo_eff_nM=ruxo_skin_nM+ruxo_oral_nM;
double JAK_inh=Emax_JAK*pow(ruxo_eff_nM,n_hill)/(pow(JAK1_IC50,n_hill)+pow(ruxo_eff_nM,n_hill));
if(JAK_inh>Emax_JAK)JAK_inh=Emax_JAK;
if(JAK_inh<0)JAK_inh=0;
double afam_cAMP=afam_on*AFAM_C/(AFAM_C+0.30);
double nbuvb_eff=NBUVB_on*1.0;
double stress_signal=0.3+(1.0-MEL)*0.5;
dxdt_NKGL=k_NKGL_up*MEL*stress_signal-k_NKGL_down*NKGL;
dxdt_NKGD_ACT=k_NKGL_up*NKGL-k_NKGL_down*NKGD_ACT;
double CD8_recruit=k_CD8_rec*CXCL10*(CD8E0/(CD8E0+0.5));
double JAK_CD8_blk=JAK_inh*0.85;
dxdt_CD8E=CD8_recruit*(1.0-JAK_CD8_blk)+k_CD8_prol*CD8E*(1.0-CD8E/3.0)-k_CD8_death*CD8E;
if(CD8E<0.01)CD8E=0.01;
dxdt_TREG=k_Treg_in*(1.0+NBUVB_Treg_stim*nbuvb_eff*5.0)-k_Treg_out*TREG;
dxdt_TREG_SKIN=k_Treg_in*0.4*(1.0+nbuvb_eff*2.0)-k_Treg_out*TREG_SKIN;
double Treg_sup=Treg_sup_IFNG*TREG;
dxdt_IFNG=k_IFNG_prod*(CD8E+CD8E*0.6)-k_IFNG_deg*IFNG-Treg_sup*IFNG;
if(IFNG<0.1)IFNG=0.1;
dxdt_PSTAT1=k_pST1_act*IFNG*(1.0-JAK_inh)-k_pST1_deg*PSTAT1;
if(PSTAT1<0)PSTAT1=0; if(PSTAT1>1)PSTAT1=1;
dxdt_CXCL10=k_CXCL10_prod*PSTAT1-k_CXCL10_deg*CXCL10;
if(CXCL10<1)CXCL10=1;
double IFNg_MITF_sup=0.4*IFNG/(IFNG+20.0);
dxdt_MITF_C=k_MITF_base+k_MITF_cAMP*afam_cAMP-k_MITF_deg*MITF_C-IFNg_MITF_sup*MITF_C+0.06*JAK_inh*MITF_C;
if(MITF_C<0)MITF_C=0; if(MITF_C>1)MITF_C=1;
double UV_mob=NBUVB_on*NBUVB_mel_stim*2.0;
dxdt_HAIRFOL=k_FHAIR_ren*HAIRFOL-k_FHAIR_mob*(1.0+UV_mob)*HAIRFOL;
double mel_killing=k_mel_kill*CD8E*MEL+0.040*NKGD_ACT*MEL;
double mel_repop=k_mel_birth*MEL+NBUVB_mel_stim*nbuvb_eff*HAIRFOL*0.8+k_MITF_cAMP*afam_cAMP*0.12*HAIRFOL;
dxdt_MEL=mel_repop-k_mel_death*MEL-mel_killing;
if(MEL<0.001)MEL=0.001; if(MEL>1)MEL=1;
dxdt_MELANIN=k_MEL_prod*MITF_C*MEL-k_MEL_loss*MELANIN;
if(MELANIN<0)MELANIN=0; if(MELANIN>1)MELANIN=1;
dxdt_INFLAM=k_inflam_in*(IFNG/20.0+CXCL10/100.0+CD8E/2.0)/3.0-k_inflam_out*INFLAM;
double VASI_worse=k_VASI_prog*CD8E*(1.0-MEL)*VASI;
double VASI_better=k_VASI_rep*MELANIN*MITF_C*(VASI-VASI_floor);
dxdt_VASI=VASI_worse-VASI_better;
if(VASI<VASI_floor)VASI=VASI_floor; if(VASI>100)VASI=100;
double repig_rate=k_repig_rate*MEL*MITF_C*(VASI0-VASI+0.1);
if(repig_rate<0)repig_rate=0;
dxdt_REPIG=repig_rate;
if(REPIG>100)REPIG=100;

$TABLE
double JAK_inh_pct=JAK_inh*100;
double pSTAT1_inh_pct=(1.0-PSTAT1/pSTAT1_0)*100;
double VASI50_resp=(VASI<=VASI0*0.50)?1:0;
double VASI75_resp=(VASI<=VASI0*0.25)?1:0;

$CAPTURE
VASI VASI50_resp VASI75_resp REPIG
CXCL10 IFNG PSTAT1 CD8E TREG MEL MITF_C MELANIN
JAK_inh_pct pSTAT1_inh_pct INFLAM HAIRFOL
RUXO_SK RUXO_C AFAM_C
'

vit_mod <- mcode("vit_shiny", vit_code)

CREAM_DOSE_NM <- 1960  # nM per topical application
VASI0_DEFAULT <- 48.0

make_events <- function(tx, dur_wk) {
  dur_h <- dur_wk * 7 * 24
  if (tx == "placebo") return(ev(time = 0, amt = 0, cmt = 1))
  if (tx == "ruxo_cream_BID") {
    times <- seq(0, dur_h - 12, by = 12)
    return(ev(amt = CREAM_DOSE_NM, cmt = "RUXO_SK", time = times))
  }
  if (tx == "ruxo_cream_QD") {
    times <- seq(0, dur_h - 24, by = 24)
    return(ev(amt = CREAM_DOSE_NM, cmt = "RUXO_SK", time = times))
  }
  if (tx == "ruxo_oral") {
    times <- seq(0, dur_h - 12, by = 12)
    return(ev(amt = 10, cmt = "RUXO_GUT", time = times))
  }
  if (tx == "afam_nbuvb") {
    times_afam <- seq(0, 3 * 60 * 24, by = 60 * 24)
    return(ev(amt = 16000, cmt = "AFAM_D", time = times_afam))
  }
  ev(time = 0, amt = 0, cmt = 1)
}

tx_params <- list(
  placebo       = list(ruxo_cream_BID=0, ruxo_cream_QD=0, ruxo_oral_BID=0, afam_on=0, NBUVB_on=0),
  ruxo_cream_BID= list(ruxo_cream_BID=1, ruxo_cream_QD=0, ruxo_oral_BID=0, afam_on=0, NBUVB_on=0),
  ruxo_cream_QD = list(ruxo_cream_BID=0, ruxo_cream_QD=1, ruxo_oral_BID=0, afam_on=0, NBUVB_on=0),
  ruxo_oral     = list(ruxo_cream_BID=0, ruxo_cream_QD=0, ruxo_oral_BID=1, afam_on=0, NBUVB_on=0),
  afam_nbuvb    = list(ruxo_cream_BID=0, ruxo_cream_QD=0, ruxo_oral_BID=0, afam_on=1, NBUVB_on=1)
)
tx_labels <- c(
  placebo        = "① Placebo (Vehicle)",
  ruxo_cream_BID = "② Ruxo Cream 1.5% BID",
  ruxo_cream_QD  = "③ Ruxo Cream 1.5% QD",
  ruxo_oral      = "④ Ruxo Oral 10mg BID",
  afam_nbuvb     = "⑤ Afamelanotide + NB-UVB"
)
tx_colors <- c(
  placebo        = "#607D8B",
  ruxo_cream_BID = "#1565C0",
  ruxo_cream_QD  = "#64B5F6",
  ruxo_oral      = "#0D47A1",
  afam_nbuvb     = "#2E7D32"
)

run_sim <- function(mod_in, tx, dur_wk, extra_params = list()) {
  params_use <- modifyList(tx_params[[tx]], extra_params)
  mod_run    <- mod_in %>% param(params_use)
  ev_use     <- make_events(tx, dur_wk)
  out        <- mrgsim(mod_run, events = ev_use,
                       end   = dur_wk * 7 * 24,
                       delta = 6,
                       carry_out = "evid")
  as.data.frame(out) %>% mutate(time_wk = time / (7 * 24), scenario = tx)
}

# ============================================================
# UI
# ============================================================
ui <- page_navbar(
  title = "Vitiligo QSP Dashboard",
  theme = bs_theme(version = 5, bootswatch = "flatly",
                   primary = "#1565C0", secondary = "#2E7D32"),
  bg = "#1565C0",

  # ---- TAB 1: Patient Profile ----------------------------------------
  nav_panel(
    title = "① Patient Profile",
    icon  = icon("user"),
    layout_sidebar(
      sidebar = sidebar(
        title = "Patient Parameters",
        width = 300,
        h6("Demographics"),
        sliderInput("age",      "Age (years)",       min = 10,  max = 75,  value = 35),
        selectInput("sex",      "Sex",               choices = c("Female", "Male")),
        sliderInput("bsa_pct",  "BSA Involvement (%)", min = 1, max = 80,  value = 15),
        hr(),
        h6("Vitiligo Characteristics"),
        selectInput("vit_type", "Vitiligo Type",
                    choices = c("Non-segmental (NSV)" = "nsv",
                                "Segmental (SV)"      = "sv",
                                "Acrofacial"          = "acrofacial"),
                    selected = "nsv"),
        sliderInput("VASI0_user", "Baseline VASI Score", min = 5, max = 90, value = 48),
        selectInput("activity",  "Disease Activity",
                    choices = c("Active (Koebner +)" = "active",
                                "Stable (>1yr)"      = "stable")),
        sliderInput("CD8E0_user", "CD8+ T cell activity (rel.)",
                    min = 0.1, max = 1.0, value = 0.30, step = 0.05),
        sliderInput("CXCL10_0_user", "Baseline CXCL10 (pg/mL)",
                    min = 20, max = 200, value = 82),
        hr(),
        h6("Simulation Duration"),
        sliderInput("dur_wk", "Duration (weeks)", min = 12, max = 52, value = 24, step = 4)
      ),
      layout_columns(
        fill = FALSE,
        value_box(title = "VASI Score",   value = textOutput("vb_vasi"),
                  showcase = icon("percent"), theme = "primary"),
        value_box(title = "BSA Involved", value = textOutput("vb_bsa"),
                  showcase = icon("body"),   theme = "info"),
        value_box(title = "Disease Activity", value = textOutput("vb_act"),
                  showcase = icon("fire"),   theme = "warning"),
        value_box(title = "Baseline CXCL10", value = textOutput("vb_cxcl10"),
                  showcase = icon("flask"),  theme = "success")
      ),
      card(
        card_header("Patient Classification & Disease Severity"),
        layout_columns(
          col_widths = c(6, 6),
          plotlyOutput("plt_vasi_gauge",   height = "280px"),
          plotlyOutput("plt_cxcl10_gauge", height = "280px")
        )
      )
    )
  ),

  # ---- TAB 2: Drug PK ------------------------------------------------
  nav_panel(
    title = "② Drug PK",
    icon  = icon("pills"),
    layout_sidebar(
      sidebar = sidebar(
        title = "PK Parameters",
        width = 300,
        selectInput("pk_drug", "Select Drug",
                    choices = c("Ruxolitinib cream BID" = "ruxo_cream_BID",
                                "Ruxolitinib cream QD"  = "ruxo_cream_QD",
                                "Ruxolitinib oral 10mg BID" = "ruxo_oral",
                                "Afamelanotide SC 16mg q60d" = "afam_nbuvb")),
        hr(),
        h6("PK Parameter Override"),
        sliderInput("JAK1_IC50_user", "JAK1 IC50 (nM)",
                    min = 1, max = 10, value = 3.3, step = 0.1),
        hr(),
        card(
          card_header("Reference PK Parameters"),
          tableOutput("tbl_pk_params")
        )
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Drug Concentration vs Time"),
             plotlyOutput("plt_pk_conc", height = "320px")),
        card(card_header("PK Summary Statistics"),
             tableOutput("tbl_pk_summary"))
      ),
      card(
        card_header("JAK Inhibition vs Skin Concentration (Dose-Response)"),
        plotlyOutput("plt_pk_dr", height = "280px")
      )
    )
  ),

  # ---- TAB 3: PD Biomarkers ------------------------------------------
  nav_panel(
    title = "③ PD Biomarkers",
    icon  = icon("dna"),
    layout_sidebar(
      sidebar = sidebar(
        title = "Treatment Selection",
        width = 260,
        checkboxGroupInput("pd_tx_sel", "Select Treatments",
                           choices   = tx_labels,
                           selected  = c("placebo", "ruxo_cream_BID",
                                         "ruxo_oral", "afam_nbuvb")),
        hr(),
        h6("Reference Biomarker Ranges"),
        tableOutput("tbl_biomarker_ref")
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("IFN-γ (pg/mL)"),
             plotlyOutput("plt_pd_ifng",   height = "260px")),
        card(card_header("Serum CXCL10 (pg/mL) ★ Disease Activity Biomarker"),
             plotlyOutput("plt_pd_cxcl10", height = "260px"))
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("p-STAT1 Inhibition (%)"),
             plotlyOutput("plt_pd_pstat1", height = "260px")),
        card(card_header("CD8⁺ Effector T Cells (relative)"),
             plotlyOutput("plt_pd_cd8",    height = "260px"))
      )
    )
  ),

  # ---- TAB 4: Clinical Endpoints -------------------------------------
  nav_panel(
    title = "④ Clinical Endpoints",
    icon  = icon("chart-line"),
    layout_sidebar(
      sidebar = sidebar(
        title = "Treatment",
        width = 260,
        selectInput("ep_tx", "Select Treatment",
                    choices = tx_labels,
                    selected = "ruxo_cream_BID"),
        hr(),
        sliderInput("ep_wk", "Evaluation Time (weeks)",
                    min = 4, max = 52, value = 24, step = 4),
        hr(),
        card(
          card_header("Endpoint Summary"),
          tableOutput("tbl_endpoint_summary")
        )
      ),
      layout_columns(
        fill = FALSE,
        value_box("VASI at Selected Week", textOutput("vb_ep_vasi"),
                  showcase = icon("chart-bar"), theme = "primary"),
        value_box("VASI50 Response",       textOutput("vb_ep_v50"),
                  showcase = icon("check"),   theme = "success"),
        value_box("VASI75 Response",       textOutput("vb_ep_v75"),
                  showcase = icon("star"),    theme = "warning"),
        value_box("Repigmentation %",      textOutput("vb_ep_repig"),
                  showcase = icon("sparkles"),theme = "info")
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("VASI Score Trajectory"),
             plotlyOutput("plt_ep_vasi",    height = "300px")),
        card(card_header("Melanocyte Density & Melanin Content"),
             plotlyOutput("plt_ep_mel",     height = "300px"))
      )
    )
  ),

  # ---- TAB 5: Scenario Comparison ------------------------------------
  nav_panel(
    title = "⑤ Scenario Comparison",
    icon  = icon("table-columns"),
    layout_columns(
      col_widths = c(8, 4),
      card(
        card_header("VASI Score — All Scenarios"),
        plotlyOutput("plt_cmp_vasi",  height = "310px")
      ),
      card(
        card_header("CXCL10 — All Scenarios"),
        plotlyOutput("plt_cmp_cxcl10", height = "310px")
      )
    ),
    layout_columns(
      col_widths = c(8, 4),
      card(
        card_header("Melanocyte Density — All Scenarios"),
        plotlyOutput("plt_cmp_mel",  height = "290px")
      ),
      card(
        card_header("Repigmentation % — All Scenarios"),
        plotlyOutput("plt_cmp_repig", height = "290px")
      )
    ),
    card(
      card_header("Scenario Summary Table at Selected Week"),
      layout_columns(
        col_widths = c(3, 9),
        sliderInput("cmp_wk", "Evaluation Week", min = 4, max = 52, value = 24, step = 4),
        DTOutput("tbl_cmp_summary")
      )
    )
  ),

  # ---- TAB 6: Biomarker & Risk ---------------------------------------
  nav_panel(
    title = "⑥ Biomarker & Risk",
    icon  = icon("microscope"),
    layout_columns(
      col_widths = c(6, 6),
      card(
        card_header("CXCL10 vs VASI Correlation"),
        plotlyOutput("plt_bm_scatter", height = "320px")
      ),
      card(
        card_header("Disease Activity Risk Calculator"),
        sliderInput("risk_cxcl10", "Serum CXCL10 (pg/mL)", 10, 250, 82),
        sliderInput("risk_vasi",   "Current VASI",           0,  80, 30),
        sliderInput("risk_dur_yr", "Disease Duration (years)", 0.5, 30, 5),
        selectInput("risk_act",    "Koebner Phenomenon",
                    choices = c("Absent" = "no", "Present" = "yes")),
        hr(),
        card(card_header("Risk Assessment"),
             tableOutput("tbl_risk"))
      )
    ),
    card(
      card_header("Serum CXCL10 as Pharmacodynamic Biomarker — Literature Summary"),
      layout_columns(
        col_widths = c(6, 6),
        plotlyOutput("plt_bm_cxcl10_resp", height = "300px"),
        tableOutput("tbl_biomarker_lit")
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  # ---- Reactive: build extra user params ----
  user_extra_params <- reactive({
    list(
      VASI0    = input$VASI0_user,
      CD8E0    = input$CD8E0_user,
      CXCL10_0 = input$CXCL10_0_user,
      IFNG0    = 12 * (input$CD8E0_user / 0.30)
    )
  })

  # ---- Reactive: run single scenario (for tabs 2,4) ----
  sim_single <- reactive({
    req(input$ep_tx)
    tryCatch(
      run_sim(vit_mod, input$ep_tx, input$dur_wk, user_extra_params()),
      error = function(e) NULL
    )
  })

  # ---- Reactive: run PD subset ----
  sim_pd <- reactive({
    req(input$pd_tx_sel)
    txs <- input$pd_tx_sel
    bind_rows(lapply(txs, function(tx) {
      tryCatch(run_sim(vit_mod, tx, input$dur_wk, user_extra_params()),
               error = function(e) NULL)
    }))
  })

  # ---- Reactive: run all scenarios (for tab 5) ----
  sim_all <- reactive({
    bind_rows(lapply(names(tx_labels), function(tx) {
      tryCatch(run_sim(vit_mod, tx, input$dur_wk, user_extra_params()),
               error = function(e) NULL)
    }))
  })

  # ---- TAB 1: Value Boxes ----
  output$vb_vasi    <- renderText(paste0(input$VASI0_user, " pts"))
  output$vb_bsa     <- renderText(paste0(input$bsa_pct, "% BSA"))
  output$vb_act     <- renderText(if (input$activity == "active") "Active" else "Stable")
  output$vb_cxcl10  <- renderText(paste0(input$CXCL10_0_user, " pg/mL"))

  output$plt_vasi_gauge <- renderPlotly({
    v <- input$VASI0_user
    plot_ly(type = "indicator", mode = "gauge+number",
            value = v,
            title = list(text = "VASI Score"),
            gauge = list(
              axis  = list(range = list(0, 100)),
              steps = list(
                list(range = c(0,  25), color = "#A5D6A7"),
                list(range = c(25, 50), color = "#FFF176"),
                list(range = c(50, 75), color = "#FFAB40"),
                list(range = c(75,100), color = "#EF9A9A")),
              threshold = list(line = list(color="red", width = 4),
                               thickness = 0.75, value = v))) %>%
      layout(margin = list(t=60,b=20,l=20,r=20))
  })

  output$plt_cxcl10_gauge <- renderPlotly({
    v <- input$CXCL10_0_user
    plot_ly(type = "indicator", mode = "gauge+number",
            value = v,
            title = list(text = "Serum CXCL10 (pg/mL)"),
            gauge = list(
              axis  = list(range = list(0, 250)),
              steps = list(
                list(range = c(0,  40), color = "#A5D6A7"),
                list(range = c(40, 80), color = "#FFF176"),
                list(range = c(80,150), color = "#FFAB40"),
                list(range = c(150,250),color = "#EF9A9A")),
              threshold = list(line = list(color="red", width = 4),
                               thickness = 0.75, value = v))) %>%
      layout(margin = list(t=60,b=20,l=20,r=20))
  })

  # ---- TAB 2: PK ----
  output$tbl_pk_params <- renderTable({
    data.frame(
      Drug         = c("Ruxolitinib oral","Ruxolitinib topical","Afamelanotide SC"),
      `F (%)`      = c("95","~10 (skin)","~90 depot"),
      `t½ (h)`     = c("3","local ~25","22"),
      `IC50 nM`    = c("3.3 (JAK1)","3.3 (JAK1)","EC50~0.3ng/mL"),
      check.names  = FALSE
    )
  }, striped = TRUE, hover = TRUE, small = TRUE)

  pk_sim <- reactive({
    tx <- input$pk_drug
    run_sim(vit_mod, tx, input$dur_wk, user_extra_params())
  })

  output$plt_pk_conc <- renderPlotly({
    df <- pk_sim(); req(nrow(df) > 0)
    if (input$pk_drug == "afam_nbuvb") {
      p <- plot_ly(df, x = ~time_wk, y = ~AFAM_C, type = "scatter", mode = "lines",
                   line = list(color="#2E7D32", width=2),
                   name = "Afamelanotide [ng/mL]")
    } else {
      y_col <- if (input$pk_drug %in% c("ruxo_cream_BID","ruxo_cream_QD")) ~RUXO_SK
               else ~RUXO_C
      y_lab <- if (input$pk_drug %in% c("ruxo_cream_BID","ruxo_cream_QD")) "Skin [nM]"
               else "Plasma [μg/L]"
      p <- plot_ly(df, x = ~time_wk, y = y_col, type = "scatter", mode = "lines",
                   line = list(color="#1565C0", width=2),
                   name = y_lab)
    }
    p %>% layout(xaxis = list(title="Time (weeks)"), yaxis = list(title="Concentration"),
                 title = "PK Concentration vs Time")
  })

  output$tbl_pk_summary <- renderTable({
    df <- pk_sim(); req(nrow(df) > 0)
    y_col <- if (input$pk_drug == "afam_nbuvb") df$AFAM_C
             else if (input$pk_drug %in% c("ruxo_cream_BID","ruxo_cream_QD")) df$RUXO_SK
             else df$RUXO_C
    unit  <- if (input$pk_drug == "afam_nbuvb") "ng/mL"
             else if (input$pk_drug %in% c("ruxo_cream_BID","ruxo_cream_QD")) "nM"
             else "μg/L"
    data.frame(
      Metric = c(paste("Cmax (", unit,")"), paste("Cmin SS (", unit,")"),
                 "t½ (h)", "JAK Inh at Cmax (%)"),
      Value  = c(round(max(y_col, na.rm=T),2),
                 round(quantile(y_col[df$time_wk > input$dur_wk*0.8], 0.1, na.rm=T),2),
                 ifelse(input$pk_drug %in% c("ruxo_cream_BID","ruxo_cream_QD"),"~25h","3h"),
                 round(max(df$JAK_inh_pct, na.rm=T),1)),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE)

  output$plt_pk_dr <- renderPlotly({
    nM  <- 10^seq(-1, 3.5, length.out = 200)
    ic50 <- input$JAK1_IC50_user
    inh <- 0.95 * nM^1.2 / (ic50^1.2 + nM^1.2) * 100
    df  <- data.frame(nM = nM, inh = inh)
    plot_ly(df, x = ~nM, y = ~inh, type="scatter", mode="lines",
            line = list(color="#1565C0", width=2)) %>%
      add_segments(x = ic50, xend = ic50, y = 0, yend = 50,
                   line = list(color="red", dash="dash"), name = paste0("IC50=",ic50,"nM")) %>%
      layout(xaxis = list(title="Ruxolitinib [nM] (log)", type="log"),
             yaxis = list(title="p-STAT1 Inhibition (%)"),
             title = "Dose-Response: JAK Inhibition")
  })

  # ---- TAB 3: PD Biomarkers ----
  output$tbl_biomarker_ref <- renderTable({
    data.frame(
      Biomarker = c("CXCL10", "IFN-γ", "p-STAT1", "CD8+"),
      Normal    = c("<40 pg/mL","<5 pg/mL","low","low"),
      `Active VIT`= c("80-200","10-30","elevated","elevated"),
      check.names=FALSE
    )
  }, striped=T, small=T)

  make_pd_plot <- function(df, y_col, y_lab, title_str) {
    req(nrow(df) > 0, y_col %in% names(df))
    colors_used <- tx_colors[unique(df$scenario)]
    labels_used <- tx_labels[unique(df$scenario)]
    plot_ly() %>%
      {
        p <- .
        for (tx in unique(df$scenario)) {
          d <- df %>% filter(scenario == tx)
          p <- add_trace(p, data = d, x = ~time_wk, y = as.formula(paste0("~", y_col)),
                         type = "scatter", mode = "lines",
                         line = list(color = tx_colors[[tx]], width = 2),
                         name = tx_labels[[tx]])
        }
        p
      } %>%
      layout(xaxis = list(title="Weeks"), yaxis = list(title = y_lab),
             title = title_str, legend = list(orientation = "h", y = -0.3))
  }

  output$plt_pd_ifng   <- renderPlotly({ make_pd_plot(sim_pd(), "IFNG",  "IFN-γ (pg/mL)",      "IFN-γ") })
  output$plt_pd_cxcl10 <- renderPlotly({ make_pd_plot(sim_pd(), "CXCL10","CXCL10 (pg/mL)",     "CXCL10") })
  output$plt_pd_pstat1 <- renderPlotly({ make_pd_plot(sim_pd(), "pSTAT1_inh_pct","pSTAT1 Inh %","p-STAT1 Inhibition") })
  output$plt_pd_cd8    <- renderPlotly({ make_pd_plot(sim_pd(), "CD8E",  "CD8+ (relative)","CD8⁺ Effectors") })

  # ---- TAB 4: Clinical Endpoints ----
  ep_at_wk <- reactive({
    df <- sim_single(); req(!is.null(df))
    df %>% filter(abs(time_wk - input$ep_wk) == min(abs(time_wk - input$ep_wk))) %>%
      slice(1)
  })

  output$vb_ep_vasi  <- renderText({ r <- ep_at_wk(); paste0(round(r$VASI, 1)) })
  output$vb_ep_v50   <- renderText({ r <- ep_at_wk(); if (r$VASI50_resp == 1) "YES ✓" else "No" })
  output$vb_ep_v75   <- renderText({ r <- ep_at_wk(); if (r$VASI75_resp == 1) "YES ✓" else "No" })
  output$vb_ep_repig <- renderText({ r <- ep_at_wk(); paste0(round(r$REPIG, 1), "%") })

  output$tbl_endpoint_summary <- renderTable({
    r <- ep_at_wk()
    data.frame(
      Endpoint = c("VASI score","VASI50","VASI75","Repigmentation","CXCL10","Melanocyte %","MITF %"),
      Value    = c(round(r$VASI,1),
                   ifelse(r$VASI50_resp==1,"Responder","Non-responder"),
                   ifelse(r$VASI75_resp==1,"Responder","Non-responder"),
                   paste0(round(r$REPIG,1),"%"),
                   paste0(round(r$CXCL10,1)," pg/mL"),
                   paste0(round(r$MEL*100,1),"%"),
                   paste0(round(r$MITF_C*100,1),"%"))
    )
  }, striped=TRUE)

  output$plt_ep_vasi <- renderPlotly({
    df <- sim_single(); req(!is.null(df))
    v0 <- input$VASI0_user
    plot_ly(df, x=~time_wk, y=~VASI, type="scatter", mode="lines",
            line=list(color=tx_colors[[input$ep_tx]], width=2.5), name="VASI") %>%
      add_segments(x=0, xend=max(df$time_wk),
                   y=v0*0.50, yend=v0*0.50,
                   line=list(color="orange", dash="dash"), name="VASI50") %>%
      add_segments(x=0, xend=max(df$time_wk),
                   y=v0*0.25, yend=v0*0.25,
                   line=list(color="red",    dash="dash"), name="VASI75") %>%
      layout(xaxis=list(title="Weeks"), yaxis=list(title="VASI Score"),
             title=paste("VASI —", tx_labels[[input$ep_tx]]))
  })

  output$plt_ep_mel <- renderPlotly({
    df <- sim_single(); req(!is.null(df))
    df <- df %>% mutate(mel_pct = MEL*100, melanin_pct = MELANIN*100)
    plot_ly(df, x=~time_wk) %>%
      add_trace(y=~mel_pct, type="scatter", mode="lines",
                line=list(color="#7B1FA2",width=2), name="Melanocyte (%)") %>%
      add_trace(y=~melanin_pct, type="scatter", mode="lines",
                line=list(color="#2E7D32",width=2), name="Melanin (%)") %>%
      layout(xaxis=list(title="Weeks"), yaxis=list(title="% of normal"),
             title="Melanocyte Density & Melanin Content")
  })

  # ---- TAB 5: Scenario Comparison ----
  make_cmp_plotly <- function(y_col, y_lab, title_str) {
    df <- sim_all(); req(nrow(df)>0)
    p <- plot_ly()
    for (tx in names(tx_labels)) {
      d <- df %>% filter(scenario == tx)
      if (nrow(d) > 0) {
        p <- add_trace(p, data=d, x=~time_wk, y=as.formula(paste0("~",y_col)),
                       type="scatter", mode="lines",
                       line=list(color=tx_colors[[tx]],width=2),
                       name=tx_labels[[tx]])
      }
    }
    p %>% layout(xaxis=list(title="Weeks"), yaxis=list(title=y_lab),
                 title=title_str, legend=list(orientation="h", y=-0.35))
  }

  output$plt_cmp_vasi  <- renderPlotly({ make_cmp_plotly("VASI",    "VASI", "VASI Score") })
  output$plt_cmp_cxcl10<- renderPlotly({ make_cmp_plotly("CXCL10", "pg/mL","CXCL10") })
  output$plt_cmp_mel   <- renderPlotly({ make_cmp_plotly("MEL",     "rel.","Melanocyte Density") })
  output$plt_cmp_repig <- renderPlotly({ make_cmp_plotly("REPIG",   "%","Repigmentation %") })

  output$tbl_cmp_summary <- renderDT({
    df <- sim_all(); req(nrow(df)>0)
    tgt_wk <- input$cmp_wk
    summ <- df %>%
      group_by(scenario) %>%
      filter(abs(time_wk - tgt_wk) == min(abs(time_wk - tgt_wk))) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(
        VASI_chg = round((input$VASI0_user - VASI)/input$VASI0_user*100,1),
        VASI50_r = ifelse(VASI50_resp==1,"Yes","No"),
        VASI75_r = ifelse(VASI75_resp==1,"Yes","No")
      ) %>%
      select(scenario, VASI, VASI_chg, VASI50_r, VASI75_r,
             CXCL10, pSTAT1_inh_pct, REPIG, JAK_inh_pct) %>%
      mutate(
        scenario = tx_labels[scenario],
        across(c(VASI, VASI_chg, CXCL10, pSTAT1_inh_pct, REPIG, JAK_inh_pct), ~round(.x,1))
      ) %>%
      rename(`Scenario`=scenario, `VASI`=VASI, `%ΔVASI`=VASI_chg,
             `VASI50`=VASI50_r, `VASI75`=VASI75_r,
             `CXCL10\n(pg/mL)`=CXCL10, `pSTAT1\nInh%`=pSTAT1_inh_pct,
             `Repig%`=REPIG, `JAK\nInh%`=JAK_inh_pct)
    datatable(summ, options=list(dom="t", pageLength=10),
              rownames=FALSE, class="stripe hover compact")
  })

  # ---- TAB 6: Biomarker & Risk ----
  output$plt_bm_scatter <- renderPlotly({
    df <- sim_all(); req(nrow(df)>0)
    plot_ly() %>%
      {
        p <- .
        for (tx in names(tx_labels)) {
          d <- df %>% filter(scenario==tx)
          if (nrow(d)>0) {
            p <- add_trace(p, data=d, x=~CXCL10, y=~VASI,
                           type="scatter", mode="markers",
                           marker=list(color=tx_colors[[tx]], opacity=0.6, size=4),
                           name=tx_labels[[tx]])
          }
        }
        p
      } %>%
      layout(xaxis=list(title="Serum CXCL10 (pg/mL)"),
             yaxis=list(title="VASI Score"),
             title="CXCL10 vs VASI — Disease Activity Correlation",
             legend=list(orientation="h", y=-0.3))
  })

  output$tbl_risk <- renderTable({
    cxcl10 <- input$risk_cxcl10
    vasi   <- input$risk_vasi
    dur    <- input$risk_dur_yr
    act    <- input$risk_act == "yes"
    score  <- 0
    if (cxcl10 > 120) score <- score + 3
    else if (cxcl10 > 80) score <- score + 2
    else if (cxcl10 > 40) score <- score + 1
    if (vasi > 50) score <- score + 3
    else if (vasi > 20) score <- score + 2
    else score <- score + 1
    if (dur < 2) score <- score + 2
    else if (dur < 5) score <- score + 1
    if (act) score <- score + 2
    risk_level <- if (score >= 8) "HIGH — Aggressive JAK inhibition indicated"
                  else if (score >= 5) "MODERATE — Topical ruxolitinib appropriate"
                  else "LOW — Topical steroids / calcineurin inhibitors"
    data.frame(
      Parameter = c("Risk Score","Risk Level","Suggested 1st-line Tx","Serum CXCL10 Level"),
      Value     = c(score, risk_level,
                    if (score >= 8) "Ruxolitinib oral 10mg BID"
                    else if (score >= 5) "Ruxolitinib cream 1.5% BID"
                    else "Topical steroid / Tacrolimus",
                    if (cxcl10 > 80) paste0(cxcl10," (elevated, active disease)")
                    else paste0(cxcl10," (borderline/normal)"))
    )
  }, striped=TRUE)

  output$plt_bm_cxcl10_resp <- renderPlotly({
    df_resp <- data.frame(
      CXCL10_baseline = c(40,60,82,110,150,180),
      VASI50_ruxo_BID = c(38,44,50,55,58,60),
      VASI50_placebo  = c(12,11,10, 9, 8, 7)
    )
    plot_ly(df_resp, x=~CXCL10_baseline) %>%
      add_trace(y=~VASI50_ruxo_BID, type="scatter", mode="lines+markers",
                line=list(color="#1565C0",width=2), name="Ruxo BID VASI50%") %>%
      add_trace(y=~VASI50_placebo,  type="scatter", mode="lines+markers",
                line=list(color="#607D8B",width=2,dash="dash"), name="Placebo VASI50%") %>%
      layout(xaxis=list(title="Baseline CXCL10 (pg/mL)"),
             yaxis=list(title="VASI50 Response Rate (%)"),
             title="CXCL10 as Predictive Biomarker of Response")
  })

  output$tbl_biomarker_lit <- renderTable({
    data.frame(
      Biomarker    = c("Serum CXCL10","pSTAT1","CD8+ TRM","IFN-γ","Serum CXCL9"),
      Threshold    = c(">80 pg/mL",">0.3 rel.","dense infiltrate",">15 pg/mL",">120 pg/mL"),
      Clinical_sig = c("Active spreading","JAK target engagement","CD8 cytotoxic lesion","Active autoimmune","Active autoimmune"),
      Reference    = c("Liu et al JAAD 2019","Rashighi Sci TrMed 2014","Tulic JID 2019","Harris JID 2012","Liu et al JAAD 2019"),
      check.names  = FALSE
    )
  }, striped=TRUE, small=TRUE)
}

shinyApp(ui = ui, server = server)
