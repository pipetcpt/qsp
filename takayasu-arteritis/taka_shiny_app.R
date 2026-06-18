## =============================================================================
## Takayasu Arteritis QSP — Interactive Shiny Dashboard
## 7 tabs: Patient Profile · Drug PK · Cytokine Biomarkers ·
##         Vascular Remodeling · Disease Activity · Scenario Comparison ·
##         Risk Assessment
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)
library(patchwork)

# ─────────────────────────────────────────────────────────────────────────────
# Embedded mrgsolve model (condensed for Shiny)
# ─────────────────────────────────────────────────────────────────────────────
taka_code <- '
$PARAM @annotated
KA_PRED:1.20:Pred absorption (1/h) | F_PRED:0.82:Bioavailability
VC_PRED:17.0:Pred Vc (L) | VP_PRED:24.0:Pred Vp (L)
CL_PRED:6.60:Pred CL (L/h) | Q_PRED:5.00:Pred Q (L/h)
FUBP_PRED:0.23:Unbound fraction
VC_TCZ:3.50:TCZ Vc (L) | CL_TCZ:0.0078:TCZ CL (L/h)
Q_TCZ:0.107:TCZ Q (L/h) | VP_TCZ:2.90:TCZ Vp (L)
KON_TCZ:0.012:TCZ kon | KOFF_TCZ:0.0001:TCZ koff
KDEG_TCZ:0.004:TCZ-IL6R degrade | IL6R_BASE:8.0:IL6R baseline
EMAX_GC:0.90:GC Emax | EC50_GC:15.0:GC EC50 (ng/mL) | HILL_GC:1.5:GC Hill
EMAX_TCZ:0.95:TCZ Emax | EC50_TCZ:1.0:TCZ EC50 (ug/mL)
EMAX_IFX:0.88:IFX Emax | EC50_IFX:0.8:IFX EC50 (ug/mL)
EMAX_MTX:0.55:MTX Emax | EC50_MTX:0.5:MTX EC50
EMAX_AZA:0.50:AZA Emax | EC50_AZA:0.5:AZA EC50
KPR_TH1:0.012:Th1 prolif | KDTH_TH1:0.010:Th1 death
KPR_TH17:0.011:Th17 prolif | KDTH_TH17:0.009:Th17 death
KPR_TREG:0.008:Treg prolif | KDTH_TREG:0.007:Treg death
KPR_MACRO:0.015:Macro recruit | KDTH_MACRO:0.013:Macro death
KPROD_IL6:3.20:IL6 prod | KDEG_IL6:0.45:IL6 deg | BASE_IL6:7.0:IL6 base
KPROD_TNF:1.50:TNF prod | KDEG_TNF:0.65:TNF deg | BASE_TNF:5.0:TNF base
KPROD_IFNG:0.80:IFNG prod | KDEG_IFNG:0.55:IFNG deg | BASE_IFNG:3.0:IFNG base
KPROD_IL17:0.60:IL17 prod | KDEG_IL17:0.40:IL17 deg | BASE_IL17:5.0:IL17 base
KPROD_CRP:1.50:CRP prod | KDEG_CRP:0.035:CRP deg | BASE_CRP:3.0:CRP base
KPROD_ESR:0.40:ESR prod | KDEG_ESR:0.020:ESR deg | BASE_ESR:18.0:ESR base
KPROD_VEGF:0.50:VEGF prod | KDEG_VEGF:0.30:VEGF deg | BASE_VEGF:120.0:VEGF base
KGR_WALL:0.0003:Wall thicken | KDECR_WALL:0.0001:Wall regress | BASE_WALL:2.0:Wall base
KGR_STEN:0.0002:Stenosis prog | KDECR_STEN:0.00005:Sten regress | MAX_STEN:0.95:Max sten
KSYN_NIH:0.008:NIH synth | KDEG_NIH:0.006:NIH degrade
USE_PRED:0:Pred flag | USE_TCZ:0:TCZ flag | USE_MTX:0:MTX flag
USE_AZA:0:AZA flag | USE_IFX:0:IFX flag
MTX_CONC:0:MTX conc | AZA_CONC:0:AZA conc | IFX_CONC:0:IFX conc

$CMT PRED_C PRED_P TCZ_C TCZ_P TCZ_B TH1 TH17 TREG MACRO CD8T
     IL6 TNFA IFNG IL17 VEGFS CRP ESR_P WALL STEN NIH_S

$INIT PRED_C=0,PRED_P=0,TCZ_C=0,TCZ_P=0,TCZ_B=0
TH1=1,TH17=1,TREG=1,MACRO=1,CD8T=1
IL6=7,TNFA=5,IFNG=3,IL17=5,VEGFS=120,CRP=3,ESR_P=18
WALL=2,STEN=0,NIH_S=0

$MAIN
double PRED_free=FUBP_PRED*(PRED_C/VC_PRED);
double GC_occ=0;
if(USE_PRED>0.5) GC_occ=EMAX_GC*pow(PRED_free,HILL_GC)/(pow(EC50_GC,HILL_GC)+pow(PRED_free,HILL_GC));
double TCZ_Cc=TCZ_C/VC_TCZ;
double TCZ_eff=0;
if(USE_TCZ>0.5) TCZ_eff=EMAX_TCZ*TCZ_Cc/(EC50_TCZ+TCZ_Cc);
double MTX_eff=0;
if(USE_MTX>0.5) MTX_eff=EMAX_MTX*MTX_CONC/(EC50_MTX+MTX_CONC);
double AZA_eff=0;
if(USE_AZA>0.5) AZA_eff=EMAX_AZA*AZA_CONC/(EC50_AZA+AZA_CONC);
double IFX_eff=0;
if(USE_IFX>0.5) IFX_eff=EMAX_IFX*IFX_CONC/(EC50_IFX+IFX_CONC);
double SUP_LYMPH=1.0-0.7*GC_occ-0.5*MTX_eff-0.45*AZA_eff;
if(SUP_LYMPH<0.05) SUP_LYMPH=0.05;

$ODE
dxdt_PRED_C=-(CL_PRED+Q_PRED)*(PRED_C/VC_PRED)+Q_PRED*(PRED_P/VP_PRED);
dxdt_PRED_P=Q_PRED*(PRED_C/VC_PRED)-Q_PRED*(PRED_P/VP_PRED);
double TCZ_Cp=TCZ_P/VP_TCZ;
dxdt_TCZ_C=-(CL_TCZ+Q_TCZ)*TCZ_Cc+Q_TCZ*TCZ_Cp-KON_TCZ*TCZ_Cc*IL6R_BASE+KOFF_TCZ*(TCZ_B/VC_TCZ);
dxdt_TCZ_P=Q_TCZ*TCZ_Cc-Q_TCZ*TCZ_Cp;
dxdt_TCZ_B=KON_TCZ*TCZ_Cc*IL6R_BASE*VC_TCZ-(KOFF_TCZ+KDEG_TCZ)*TCZ_B;
double TH1s=1+0.4*(IL6/BASE_IL6-1)+0.3*(IFNG/BASE_IFNG-1);
if(TH1s<0.1)TH1s=0.1;
dxdt_TH1=KPR_TH1*TH1s*TH1*SUP_LYMPH-KDTH_TH1*TH1;
double TH17s=1+0.45*(IL6/BASE_IL6-1)+0.2*(IL17/BASE_IL17-1);
if(TH17s<0.1)TH17s=0.1;
dxdt_TH17=KPR_TH17*TH17s*TH17*SUP_LYMPH-KDTH_TH17*TH17-0.3*TREG*TH17*0.005;
double TRb=1+0.2*GC_occ;
dxdt_TREG=KPR_TREG*TRb*TREG-KDTH_TREG*TREG;
double MRec=1+0.5*(IL6/BASE_IL6-1)+0.3*(TNFA/BASE_TNF-1);
if(MRec<0.1)MRec=0.1;
dxdt_MACRO=KPR_MACRO*MRec*MACRO*(1-0.6*GC_occ)-KDTH_MACRO*MACRO;
double CD8s=1+0.35*(IFNG/BASE_IFNG-1)+0.2*(IL6/BASE_IL6-1);
if(CD8s<0.1)CD8s=0.1;
dxdt_CD8T=KPR_CD8*CD8s*CD8T*SUP_LYMPH-KDTH_CD8*CD8T;
double IL6p=MACRO*TH17*(1-0.85*GC_occ)*(1-0.10*TCZ_eff);
if(IL6p<0)IL6p=0;
dxdt_IL6=KPROD_IL6*IL6p-KDEG_IL6*IL6;
double TNFp=MACRO*TH1*(1-0.80*GC_occ)*(1-IFX_eff);
if(TNFp<0)TNFp=0;
dxdt_TNFA=KPROD_TNF*TNFp-KDEG_TNF*TNFA;
double IFNp=(TH1*1.2+CD8T*0.5)*(1-0.70*GC_occ);
if(IFNp<0)IFNp=0;
dxdt_IFNG=KPROD_IFNG*IFNp-KDEG_IFNG*IFNG;
double IL17p=TH17*(1-0.60*GC_occ-0.4*TCZ_eff);
if(IL17p<0)IL17p=0;
dxdt_IL17=KPROD_IL17*IL17p-KDEG_IL17*IL17;
double VEGFp=(1+0.8*(IL6/BASE_IL6-1))*(1-0.50*TCZ_eff)*(1-0.30*GC_occ);
if(VEGFp<0)VEGFp=0;
dxdt_VEGFS=KPROD_VEGF*BASE_VEGF*VEGFp-KDEG_VEGF*VEGFS;
double CRPs=(IL6/BASE_IL6)*(1-0.95*TCZ_eff)*(1-0.50*GC_occ);
if(CRPs<0)CRPs=0;
dxdt_CRP=KPROD_CRP*CRPs-KDEG_CRP*CRP;
double ESRs=(CRP/BASE_CRP+IL6/BASE_IL6)*0.5;
if(ESRs<0)ESRs=0;
dxdt_ESR_P=KPROD_ESR*ESRs-KDEG_ESR*ESR_P;
double LOAD=(IL6/BASE_IL6+TNFA/BASE_TNF+IFNG/BASE_IFNG+IL17/BASE_IL17)/4*MACRO;
dxdt_WALL=KGR_WALL*LOAD*(1-0.4*GC_occ-0.3*TCZ_eff)-KDECR_WALL*(WALL-BASE_WALL);
double STCAP=1-STEN/MAX_STEN;
dxdt_STEN=KGR_STEN*LOAD*STCAP*(1-0.3*GC_occ-0.2*TCZ_eff)-KDECR_STEN*STEN;
double NIHd=0.6*(CRP/10)+0.2*(IFNG/BASE_IFNG)+0.2*(IL6/BASE_IL6);
if(NIHd<0)NIHd=0;
dxdt_NIH_S=KSYN_NIH*NIHd-KDEG_NIH*NIH_S;

$TABLE
capture Pred_free=FUBP_PRED*PRED_C/VC_PRED;
capture TCZ_conc=TCZ_C/VC_TCZ;
capture GC_occ_out=EMAX_GC*pow(Pred_free,HILL_GC)/(pow(EC50_GC,HILL_GC)+pow(Pred_free,HILL_GC));
capture NIH_disc=(NIH_S<0.5)?0:(NIH_S<1.5)?1:(NIH_S<2.5)?2:(NIH_S<3.5)?3:4;
'

mod_shiny <- mcode("taka_shiny", taka_code, quiet = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# Active disease initial conditions (used as default)
# ─────────────────────────────────────────────────────────────────────────────
active_init <- list(
    TH1=2.8, TH17=3.2, TREG=0.7, MACRO=2.5, CD8T=2.0,
    IL6=85, TNFA=30, IFNG=20, IL17=35, VEGFS=450,
    CRP=55, ESR_P=68, WALL=5.0, STEN=0.18, NIH_S=2.5
)

# ─────────────────────────────────────────────────────────────────────────────
# Simulation helper
# ─────────────────────────────────────────────────────────────────────────────
run_sim <- function(scenario, bw_kg = 60, sim_days = 365,
                    il6_init = 85, crp_init = 55, sten_init = 0.18) {
    init_list <- active_init
    init_list$IL6  <- il6_init
    init_list$CRP  <- crp_init
    init_list$STEN <- sten_init

    params_map <- list(
        "Untreated" = c(USE_PRED=0, USE_TCZ=0, USE_MTX=0, USE_AZA=0, USE_IFX=0,
                        MTX_CONC=0, AZA_CONC=0, IFX_CONC=0),
        "Prednisolone (0.8 mg/kg)" = c(USE_PRED=1, USE_TCZ=0, USE_MTX=0, USE_AZA=0, USE_IFX=0,
                                        MTX_CONC=0, AZA_CONC=0, IFX_CONC=0),
        "Pred + MTX" = c(USE_PRED=1, USE_TCZ=0, USE_MTX=1, USE_AZA=0, USE_IFX=0,
                         MTX_CONC=1, AZA_CONC=0, IFX_CONC=0),
        "Tocilizumab IV q4w" = c(USE_PRED=0, USE_TCZ=1, USE_MTX=0, USE_AZA=0, USE_IFX=0,
                                  MTX_CONC=0, AZA_CONC=0, IFX_CONC=0),
        "Pred + Tocilizumab" = c(USE_PRED=1, USE_TCZ=1, USE_MTX=0, USE_AZA=0, USE_IFX=0,
                                  MTX_CONC=0, AZA_CONC=0, IFX_CONC=0),
        "Pred + AZA" = c(USE_PRED=1, USE_TCZ=0, USE_MTX=0, USE_AZA=1, USE_IFX=0,
                         MTX_CONC=0, AZA_CONC=1, IFX_CONC=0),
        "Infliximab" = c(USE_PRED=0, USE_TCZ=0, USE_MTX=0, USE_AZA=0, USE_IFX=1,
                         MTX_CONC=0, AZA_CONC=0, IFX_CONC=2.5)
    )

    ev_map <- list(
        "Untreated"              = ev(time=0, amt=0, cmt="PRED_C"),
        "Prednisolone (0.8 mg/kg)" = ev(cmt="PRED_C",
                                         amt=0.8*bw_kg*1e6*0.82,
                                         ii=24, addl=sim_days-1, time=0),
        "Pred + MTX"             = ev(cmt="PRED_C",
                                       amt=0.6*bw_kg*1e6*0.82,
                                       ii=24, addl=sim_days-1, time=0),
        "Tocilizumab IV q4w"     = ev(cmt="TCZ_C",
                                       amt=8*bw_kg*1000,
                                       ii=28*24, addl=floor(sim_days/28), time=0),
        "Pred + Tocilizumab"     = c(
            ev(cmt="PRED_C", amt=0.5*bw_kg*1e6*0.82, ii=24,
               addl=sim_days-1, time=0),
            ev(cmt="TCZ_C",  amt=8*bw_kg*1000, ii=28*24,
               addl=floor(sim_days/28), time=0)
        ),
        "Pred + AZA"             = ev(cmt="PRED_C",
                                       amt=0.5*bw_kg*1e6*0.82,
                                       ii=24, addl=sim_days-1, time=0),
        "Infliximab"             = ev(time=0, amt=0, cmt="PRED_C")
    )

    out <- tryCatch({
        mod_shiny %>%
            init(init_list) %>%
            param(params_map[[scenario]]) %>%
            mrgsim(events = ev_map[[scenario]],
                   end = sim_days * 24, delta = 12) %>%
            as.data.frame() %>%
            mutate(time_day = time / 24, scenario = scenario)
    }, error = function(e) {
        data.frame(time=0, time_day=0, scenario=scenario,
                   CRP=NA, IL6=NA, NIH_S=NA, WALL=NA, STEN=NA, VEGFS=NA,
                   Pred_free=NA, TCZ_conc=NA, GC_occ_out=NA, NIH_disc=NA)
    })
    out
}

# ─────────────────────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
    skin = "red",
    dashboardHeader(title = "Takayasu Arteritis QSP"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Patient Profile",     tabName = "tab_profile",  icon = icon("user")),
            menuItem("Drug PK",             tabName = "tab_pk",       icon = icon("capsules")),
            menuItem("Cytokine Biomarkers", tabName = "tab_cyto",     icon = icon("vials")),
            menuItem("Vascular Remodeling", tabName = "tab_vasc",     icon = icon("heartbeat")),
            menuItem("Disease Activity",    tabName = "tab_activity", icon = icon("chart-line")),
            menuItem("Scenario Comparison", tabName = "tab_compare",  icon = icon("table")),
            menuItem("Risk Assessment",     tabName = "tab_risk",     icon = icon("exclamation-triangle"))
        )
    ),
    dashboardBody(
        tabItems(

# ─── TAB 1: Patient Profile ────────────────────────────────────────────────
tabItem(tabName = "tab_profile",
    fluidRow(
        box(title = "Patient Demographics & Disease Parameters", width = 4,
            status = "danger", solidHeader = TRUE,
            numericInput("bw", "Body weight (kg)", value = 60, min = 30, max = 120),
            numericInput("age", "Age (years)", value = 32, min = 10, max = 80),
            selectInput("sex", "Sex", choices = c("Female", "Male"), selected = "Female"),
            selectInput("disease_phase", "Disease Phase (ACR)",
                choices = c("Phase I (pre-pulseless)" = "I",
                            "Phase IIa (vascular symptoms)" = "IIa",
                            "Phase IIb (mixed)" = "IIb",
                            "Phase III (pulseless)" = "III",
                            "Phase IV (+hypertension)" = "IV"),
                selected = "IIb"),
            sliderInput("il6_init", "Baseline IL-6 (pg/mL)", min = 5, max = 300, value = 85),
            sliderInput("crp_init", "Baseline CRP (mg/L)", min = 1, max = 200, value = 55),
            sliderInput("sten_init", "Baseline Stenosis Score", min = 0, max = 0.8, value = 0.18, step = 0.01),
            selectInput("treatment", "Treatment", choices = c(
                "Untreated",
                "Prednisolone (0.8 mg/kg)",
                "Pred + MTX",
                "Tocilizumab IV q4w",
                "Pred + Tocilizumab",
                "Pred + AZA",
                "Infliximab"
            )),
            sliderInput("sim_days", "Simulation duration (days)", 90, 730, 365, 30),
            actionButton("run_sim", "Run Simulation", class = "btn-danger", icon = icon("play"))
        ),
        box(title = "Disease Summary", width = 4, status = "danger",
            valueBoxOutput("vb_crp", width = 12),
            valueBoxOutput("vb_il6", width = 12),
            valueBoxOutput("vb_nih", width = 12),
            valueBoxOutput("vb_wall", width = 12)
        ),
        box(title = "ACR Diagnostic Criteria", width = 4, status = "warning",
            solidHeader = TRUE,
            HTML('<table class="table table-sm table-bordered">
                <thead><tr><th>Criterion</th><th>Feature</th></tr></thead>
                <tbody>
                <tr><td>1</td><td>Age at onset ≤ 40 years</td></tr>
                <tr><td>2</td><td>Claudication of extremity</td></tr>
                <tr><td>3</td><td>Decreased brachial artery pulse</td></tr>
                <tr><td>4</td><td>BP difference > 10 mmHg</td></tr>
                <tr><td>5</td><td>Subclavian or aortic bruit</td></tr>
                <tr><td>6</td><td>Arteriographic abnormality</td></tr>
                </tbody></table>
                <p><strong>Diagnosis: ≥3 of 6 criteria (Sens 90.5%, Spec 97.8%)</strong></p>')
        )
    ),
    fluidRow(
        box(title = "Vessel Distribution & Angiographic Classification", width = 6,
            status = "danger",
            HTML('<table class="table table-sm">
            <thead><tr><th>Type</th><th>Vessel Involvement</th><th>Frequency</th></tr></thead>
            <tbody>
            <tr><td>Type I</td><td>Aortic arch branches only</td><td>~8%</td></tr>
            <tr><td>Type IIa</td><td>Ascending aorta + arch + branches</td><td>~11%</td></tr>
            <tr><td>Type IIb</td><td>+Descending thoracic aorta</td><td>~11%</td></tr>
            <tr><td>Type III</td><td>Descending thoracic + abdominal aorta ± branches</td><td>~11%</td></tr>
            <tr><td>Type IV</td><td>Abdominal aorta + renal arteries</td><td>~18%</td></tr>
            <tr><td>Type V</td><td>Combined Type IIb + IV (most common)</td><td>~41%</td></tr>
            </tbody></table>'),
            p("Plus: ±pulmonary artery involvement (Type P)")
        ),
        box(title = "Epidemiology & Key Features", width = 6, status = "info",
            HTML('<ul>
            <li><strong>Prevalence:</strong> ~2.6/million (USA); ~40/million (Japan/Korea)</li>
            <li><strong>Age of onset:</strong> 10–40 years (peak 20s–30s)</li>
            <li><strong>Female predominance:</strong> 8:1 (F:M)</li>
            <li><strong>Geographic:</strong> Asia > Middle East > Latin America > Europe</li>
            <li><strong>Genetic:</strong> HLA-B*52:01 (Japan/Korea), RNF213 p.R4810K</li>
            <li><strong>Relapse rate:</strong> ~50–80% off treatment</li>
            <li><strong>5-year survival:</strong> >90% with modern treatment</li>
            </ul>')
        )
    )
),

# ─── TAB 2: Drug PK ────────────────────────────────────────────────────────
tabItem(tabName = "tab_pk",
    fluidRow(
        box(title = "Pharmacokinetic Parameters", width = 3, status = "primary",
            solidHeader = TRUE,
            h5("Prednisolone PK"),
            p("• ka = 1.20 h⁻¹, F = 82%"),
            p("• Vc = 17 L, Vp = 24 L"),
            p("• CL = 6.6 L/h → t½ ≈ 3.5 h"),
            p("• fu = 23% (CBG-bound 77%)"),
            hr(),
            h5("Tocilizumab PK (TMDD)"),
            p("• IV: Vc = 3.5 L, t½ ≈ 6–7 d"),
            p("• SC 162 mg: F = 80%, ka = 0.015 h⁻¹"),
            p("• TMDD: kon/koff/kdeg modeled"),
            p("• Concentration-dependent clearance"),
            hr(),
            h5("Infliximab PK"),
            p("• Vc ~ 3 L, t½ ≈ 8–10 d"),
            p("• Kd ~ 0.1 nM (anti-TNF)")
        ),
        box(title = "Prednisolone: Concentration–Time", width = 9, status = "primary",
            plotlyOutput("pk_pred_plot", height = "350px")
        )
    ),
    fluidRow(
        box(title = "Tocilizumab: Central Concentration", width = 6, status = "success",
            plotlyOutput("pk_tcz_plot", height = "300px")
        ),
        box(title = "GC Receptor Occupancy (%)", width = 6, status = "warning",
            plotlyOutput("pk_gcocc_plot", height = "300px")
        )
    )
),

# ─── TAB 3: Cytokine Biomarkers ────────────────────────────────────────────
tabItem(tabName = "tab_cyto",
    fluidRow(
        box(title = "IL-6 Dynamics (pg/mL)", width = 6, status = "danger",
            plotlyOutput("cyto_il6", height = "300px"),
            p("Reference range: < 7 pg/mL | Active disease: 50–300 pg/mL")
        ),
        box(title = "CRP (mg/L)", width = 6, status = "danger",
            plotlyOutput("cyto_crp", height = "300px"),
            p("Normal: < 5 mg/L | Active disease: 20–150 mg/L")
        )
    ),
    fluidRow(
        box(title = "TNF-α & IFN-γ (pg/mL)", width = 6, status = "warning",
            plotlyOutput("cyto_tnf_ifn", height = "300px")
        ),
        box(title = "ESR (mm/h) & Serum VEGF (pg/mL)", width = 6, status = "warning",
            plotlyOutput("cyto_esr_vegf", height = "300px")
        )
    ),
    fluidRow(
        box(title = "Cytokine Summary Table", width = 12, status = "info",
            DTOutput("cyto_table")
        )
    )
),

# ─── TAB 4: Vascular Remodeling ────────────────────────────────────────────
tabItem(tabName = "tab_vasc",
    fluidRow(
        box(title = "Aortic Wall Thickness (mm)", width = 6, status = "danger",
            solidHeader = TRUE,
            plotlyOutput("vasc_wall", height = "320px"),
            p("Normal: ~2 mm | Active disease: 4–8 mm by PET-CT/MRI")
        ),
        box(title = "Stenosis Score (0–100%)", width = 6, status = "danger",
            solidHeader = TRUE,
            plotlyOutput("vasc_sten", height = "320px"),
            p("Irreversible vascular damage accumulates over months–years")
        )
    ),
    fluidRow(
        box(title = "Vascular Remodeling vs Inflammatory Load", width = 6,
            status = "warning",
            plotlyOutput("vasc_inflam", height = "280px")
        ),
        box(title = "Vascular Remodeling Mechanism", width = 6, status = "info",
            HTML('<ul>
            <li><strong>Granuloma formation:</strong> Adventitia → Giant cells → MMP-9/12 release</li>
            <li><strong>Intimal hyperplasia:</strong> SMC proliferation + migration (PDGF, FGF-2)</li>
            <li><strong>Stenosis:</strong> Progressive lumen narrowing (irreversible)</li>
            <li><strong>Aneurysm:</strong> Elastin degradation → wall weakening</li>
            <li><strong>Neovascularization:</strong> VEGF → Vasa vasorum expansion</li>
            <li><strong>Calcification:</strong> Chronic granuloma → mural calcification</li>
            </ul>
            <p><strong>¹⁸F-FDG PET-CT</strong> detects active wall inflammation (FDG-avid granuloma)</p>
            <p>Wall thickness correlates with FDG SUVmax (r=0.64)</p>')
        )
    )
),

# ─── TAB 5: Disease Activity ───────────────────────────────────────────────
tabItem(tabName = "tab_activity",
    fluidRow(
        box(title = "NIH Activity Score (continuous)", width = 8,
            status = "danger", solidHeader = TRUE,
            plotlyOutput("activity_nih", height = "350px")
        ),
        box(title = "NIH Score Components", width = 4, status = "warning",
            HTML('<table class="table table-sm">
            <thead><tr><th>Feature</th><th>Score</th></tr></thead>
            <tbody>
            <tr><td>Systemic features (fever, MSK pain, rash)</td><td>+1</td></tr>
            <tr><td>Elevated ESR (>20 mm/h F; >15 mm/h M)</td><td>+1</td></tr>
            <tr><td>Vascular ischemia/inflammation</td><td>+1</td></tr>
            <tr><td>Angiographic features</td><td>+1</td></tr>
            </tbody></table>
            <hr>
            <p><strong>Active disease: NIH ≥ 2</strong></p>
            <p><strong>Remission: NIH = 0</strong></p>
            <p><strong>Relapse: ≥2 features after remission</strong></p>'),
            hr(),
            valueBoxOutput("vb_remission_rate", width = 12)
        )
    ),
    fluidRow(
        box(title = "Immune Cell Dynamics (fold over baseline)", width = 12,
            status = "info",
            plotlyOutput("activity_immune", height = "300px")
        )
    )
),

# ─── TAB 6: Scenario Comparison ────────────────────────────────────────────
tabItem(tabName = "tab_compare",
    fluidRow(
        box(title = "All Scenarios: CRP & IL-6", width = 12,
            status = "danger", solidHeader = TRUE,
            actionButton("run_all", "Run All 7 Scenarios", class = "btn-danger",
                         icon = icon("play")),
            plotlyOutput("compare_plot", height = "450px")
        )
    ),
    fluidRow(
        box(title = "12-Month Outcome Comparison", width = 12, status = "warning",
            DTOutput("compare_table")
        )
    ),
    fluidRow(
        box(title = "Clinical Trial Reference Data", width = 12, status = "info",
            HTML('<table class="table table-sm table-bordered">
            <thead class="thead-dark"><tr>
            <th>Treatment</th><th>Trial</th><th>N</th>
            <th>Remission 12mo</th><th>CRP norm</th><th>Relapse</th></tr></thead>
            <tbody>
            <tr><td>Prednisolone</td><td>Kerr 1994</td><td>60</td><td>33%</td><td>38%</td><td>53%</td></tr>
            <tr><td>Pred + MTX</td><td>Keser 2014 (RCT)</td><td>30</td><td>55%</td><td>52%</td><td>47%</td></tr>
            <tr><td>Tocilizumab IV</td><td>TAKT 2018 (RCT)</td><td>36</td><td>50.5%</td><td>61%</td><td>34%</td></tr>
            <tr><td>Pred + TCZ</td><td>Nakaoka 2018</td><td>36</td><td>52%</td><td>65%</td><td>28%</td></tr>
            <tr><td>Pred + AZA</td><td>Valsakumar 2003</td><td>82</td><td>42%</td><td>45%</td><td>50%</td></tr>
            <tr><td>Infliximab</td><td>Comarmond 2012</td><td>15</td><td>60%</td><td>70%</td><td>38%</td></tr>
            </tbody></table>')
        )
    )
),

# ─── TAB 7: Risk Assessment ────────────────────────────────────────────────
tabItem(tabName = "tab_risk",
    fluidRow(
        box(title = "5-Year Cumulative Risk Projection", width = 8,
            status = "danger", solidHeader = TRUE,
            plotlyOutput("risk_plot", height = "350px")
        ),
        box(title = "Risk Stratification", width = 4, status = "warning",
            h5("5-Year Complication Risk"),
            verbatimTextOutput("risk_text"),
            hr(),
            h5("Poor Prognosis Factors"),
            HTML('<ul>
            <li>Severe aortic regurgitation</li>
            <li>Pulmonary artery involvement</li>
            <li>Renovascular hypertension</li>
            <li>Cardiac complications</li>
            <li>Disease onset at age < 25 years</li>
            <li>Continuous relapses (monophasic → worse prognosis)</li>
            </ul>
            <hr>
            <p><strong>Mortality predictor (Ishikawa criteria):</strong></p>
            <p>Major: AR, hypertension, aortic regurgitation, aneurysm</p>
            <p>Minor: retinopathy, limb ischemia, carotidynia</p>')
        )
    ),
    fluidRow(
        box(title = "Organ-Specific Complication Rates by Phase", width = 12,
            status = "info",
            HTML('<table class="table table-sm table-bordered">
            <thead class="thead-dark"><tr>
            <th>Complication</th><th>Overall (%)</th><th>Phase I–II</th><th>Phase III–IV</th><th>Mechanism</th></tr></thead>
            <tbody>
            <tr><td>Renovascular HTN</td><td>33–83%</td><td>Low</td><td>High</td><td>Renal artery stenosis → RAAS</td></tr>
            <tr><td>Aortic regurgitation</td><td>14–24%</td><td>Low</td><td>Moderate</td><td>Aortic root dilation/AR</td></tr>
            <tr><td>Stroke/TIA</td><td>5–20%</td><td>Low</td><td>High</td><td>Carotid/vertebral stenosis</td></tr>
            <tr><td>Heart failure</td><td>2–15%</td><td>Rare</td><td>Moderate</td><td>HTN cardiomyopathy + AR</td></tr>
            <tr><td>Retinopathy</td><td>10–34%</td><td>Rare</td><td>Moderate</td><td>Reduced carotid flow</td></tr>
            <tr><td>Pulmonary HTN</td><td>5–17%</td><td>Low</td><td>Moderate</td><td>Pulmonary arteritis</td></tr>
            <tr><td>MI</td><td>2–10%</td><td>Rare</td><td>Moderate</td><td>Coronary ostia stenosis</td></tr>
            <tr><td>Death (5-yr)</td><td>3–10%</td><td>Low</td><td>Higher</td><td>Cardiac/cerebrovascular</td></tr>
            </tbody></table>')
        )
    )
)

# End tabItems
        )  # end tabItems
    )      # end dashboardBody
)          # end dashboardPage

# ─────────────────────────────────────────────────────────────────────────────
# Server
# ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

    sim_result <- eventReactive(input$run_sim, {
        run_sim(
            scenario  = input$treatment,
            bw_kg     = input$bw,
            sim_days  = input$sim_days,
            il6_init  = input$il6_init,
            crp_init  = input$crp_init,
            sten_init = input$sten_init
        )
    })

    all_results <- eventReactive(input$run_all, {
        sc_list <- c("Untreated","Prednisolone (0.8 mg/kg)","Pred + MTX",
                     "Tocilizumab IV q4w","Pred + Tocilizumab","Pred + AZA",
                     "Infliximab")
        withProgress(message = "Simulating all scenarios...", {
            bind_rows(lapply(sc_list, function(sc) {
                run_sim(sc, bw_kg = input$bw, sim_days = input$sim_days,
                        il6_init = input$il6_init,
                        crp_init = input$crp_init,
                        sten_init = input$sten_init)
            }))
        })
    })

    # ── Value boxes ────────────────────────────────────────────────────────
    output$vb_crp <- renderValueBox({
        d <- sim_result()
        v <- round(tail(d$CRP, 1), 1)
        clr <- if(v < 5) "green" else if(v < 20) "yellow" else "red"
        valueBox(paste0(v, " mg/L"), "CRP (final)", icon = icon("flask"), color = clr)
    })
    output$vb_il6 <- renderValueBox({
        d <- sim_result()
        v <- round(tail(d$IL6, 1), 1)
        clr <- if(v < 7) "green" else if(v < 30) "yellow" else "red"
        valueBox(paste0(v, " pg/mL"), "IL-6 (final)", icon = icon("dna"), color = clr)
    })
    output$vb_nih <- renderValueBox({
        d <- sim_result()
        v <- round(tail(d$NIH_S, 1), 2)
        clr <- if(v < 1) "green" else if(v < 2) "yellow" else "red"
        valueBox(round(v, 1), "NIH Score (final)", icon = icon("chart-bar"), color = clr)
    })
    output$vb_wall <- renderValueBox({
        d <- sim_result()
        v <- round(tail(d$WALL, 1), 2)
        clr <- if(v < 3) "green" else if(v < 5) "yellow" else "red"
        valueBox(paste0(v, " mm"), "Wall Thickness (final)", icon = icon("circle"), color = clr)
    })
    output$vb_remission_rate <- renderValueBox({
        d <- sim_result()
        rem_pct <- round(100 * mean(d$NIH_disc == 0, na.rm = TRUE), 0)
        clr <- if(rem_pct > 60) "green" else if(rem_pct > 30) "yellow" else "red"
        valueBox(paste0(rem_pct, "%"), "Time in remission", icon = icon("check-circle"), color = clr)
    })

    # ── PK plots ───────────────────────────────────────────────────────────
    output$pk_pred_plot <- renderPlotly({
        d <- sim_result()
        p <- ggplot(d[d$time_day <= 14,], aes(time_day, Pred_free)) +
            geom_line(color="#E63946", linewidth=1) +
            labs(x="Day", y="Free Prednisolone (ng/mL)",
                 title="Prednisolone Free Concentration (first 14 days)") +
            theme_classic()
        ggplotly(p)
    })
    output$pk_tcz_plot <- renderPlotly({
        d <- sim_result()
        p <- ggplot(d, aes(time_day, TCZ_conc)) +
            geom_line(color="#2A9D8F", linewidth=1) +
            geom_hline(yintercept = 1, linetype="dashed", color="gray60") +
            labs(x="Day", y="Tocilizumab (ug/mL)",
                 title="Tocilizumab Central Concentration") +
            theme_classic()
        ggplotly(p)
    })
    output$pk_gcocc_plot <- renderPlotly({
        d <- sim_result()
        d$GC_pct <- d$GC_occ_out * 100
        p <- ggplot(d[d$time_day <= 14,], aes(time_day, GC_pct)) +
            geom_line(color="#F4A261", linewidth=1) +
            geom_hline(yintercept = 50, linetype="dashed", color="gray60") +
            labs(x="Day", y="GR Occupancy (%)",
                 title="Glucocorticoid Receptor Occupancy") +
            theme_classic()
        ggplotly(p)
    })

    # ── Cytokine plots ─────────────────────────────────────────────────────
    output$cyto_il6 <- renderPlotly({
        d <- sim_result()
        p <- ggplot(d, aes(time_day, IL6)) +
            geom_line(color="#E63946", linewidth=1) +
            geom_hline(yintercept=7, linetype="dashed", color="gray60") +
            labs(x="Day", y="IL-6 (pg/mL)") + theme_classic()
        ggplotly(p)
    })
    output$cyto_crp <- renderPlotly({
        d <- sim_result()
        p <- ggplot(d, aes(time_day, CRP)) +
            geom_line(color="#E9C46A", linewidth=1) +
            geom_hline(yintercept=5, linetype="dashed", color="gray60") +
            labs(x="Day", y="CRP (mg/L)") + theme_classic()
        ggplotly(p)
    })
    output$cyto_tnf_ifn <- renderPlotly({
        d <- sim_result() %>%
            select(time_day, TNFA, IFNG) %>%
            pivot_longer(c(TNFA, IFNG), names_to="marker", values_to="conc")
        p <- ggplot(d, aes(time_day, conc, color=marker)) +
            geom_line(linewidth=1) +
            scale_color_manual(values=c(TNFA="#264653", IFNG="#457B9D")) +
            labs(x="Day", y="pg/mL", color=NULL) + theme_classic()
        ggplotly(p)
    })
    output$cyto_esr_vegf <- renderPlotly({
        d <- sim_result()
        p1 <- ggplot(d, aes(time_day, ESR_P)) +
            geom_line(color="#6D6875", linewidth=1) +
            geom_hline(yintercept=20, linetype="dashed", color="gray60") +
            labs(x="Day", y="ESR (mm/h)") + theme_classic()
        p2 <- ggplot(d, aes(time_day, VEGFS)) +
            geom_line(color="#B5838D", linewidth=1) +
            geom_hline(yintercept=240, linetype="dashed", color="gray60") +
            labs(x="Day", y="VEGF (pg/mL)") + theme_classic()
        subplot(ggplotly(p1), ggplotly(p2), nrows=1, shareX=TRUE)
    })
    output$cyto_table <- renderDT({
        d <- sim_result()
        tbl <- d %>%
            filter(time_day %in% c(0, 30, 90, 180, 365)) %>%
            group_by(time_day) %>% slice(1) %>% ungroup() %>%
            mutate(across(c(IL6, TNFA, IFNG, IL17, CRP, ESR_P, VEGFS), round, 1)) %>%
            select(`Day`=time_day, `IL-6`=IL6, `TNF-α`=TNFA,
                   `IFN-γ`=IFNG, `IL-17`=IL17,
                   `CRP`=CRP, `ESR`=ESR_P, `VEGF`=VEGFS)
        datatable(tbl, options=list(dom='t', pageLength=10),
                  rownames=FALSE) %>%
            formatStyle("CRP",
                backgroundColor = styleInterval(c(5, 20),
                    c("#d4edda","#fff3cd","#f8d7da")))
    })

    # ── Vascular plots ────────────────────────────────────────────────────
    output$vasc_wall <- renderPlotly({
        d <- sim_result()
        p <- ggplot(d, aes(time_day, WALL)) +
            geom_line(color="#6A1B9A", linewidth=1.2) +
            geom_hline(yintercept=2, linetype="dashed", color="gray60") +
            labs(x="Day", y="Wall thickness (mm)") + theme_classic()
        ggplotly(p)
    })
    output$vasc_sten <- renderPlotly({
        d <- sim_result()
        p <- ggplot(d, aes(time_day, STEN*100)) +
            geom_line(color="#AB47BC", linewidth=1.2) +
            ylim(0, 100) +
            labs(x="Day", y="Stenosis (%)") + theme_classic()
        ggplotly(p)
    })
    output$vasc_inflam <- renderPlotly({
        d <- sim_result()
        d$inflam_load <- (d$IL6/7 + d$TNFA/5 + d$IFNG/3 + d$IL17/5) / 4
        p <- ggplot(d, aes(inflam_load, WALL, color=time_day)) +
            geom_path(linewidth=1) +
            scale_color_viridis_c(name="Day") +
            labs(x="Cytokine Burden Index", y="Wall Thickness (mm)",
                 title="Vascular Remodeling Phase Portrait") +
            theme_classic()
        ggplotly(p)
    })

    # ── Disease Activity ──────────────────────────────────────────────────
    output$activity_nih <- renderPlotly({
        d <- sim_result()
        p <- ggplot(d, aes(time_day, NIH_S)) +
            geom_line(color="#E63946", linewidth=1.2) +
            geom_ribbon(aes(ymin=0, ymax=NIH_S), fill="#E63946", alpha=0.15) +
            geom_hline(yintercept=2, linetype="dashed", color="darkred",
                       linewidth=0.8) +
            annotate("text", x=max(d$time_day)*0.6, y=2.2,
                     label="Active disease (NIH≥2)", size=3.5, color="darkred") +
            labs(x="Day", y="NIH Activity Score (continuous)") +
            theme_classic()
        ggplotly(p)
    })
    output$activity_immune <- renderPlotly({
        d <- sim_result() %>%
            select(time_day, TH1, TH17, TREG, MACRO, CD8T) %>%
            pivot_longer(-time_day, names_to="cell", values_to="fold")
        p <- ggplot(d, aes(time_day, fold, color=cell)) +
            geom_line(linewidth=0.9) +
            geom_hline(yintercept=1, linetype="dashed", color="gray60") +
            scale_color_manual(values=c(TH1="#E63946",TH17="#2A9D8F",
                                        TREG="#264653",MACRO="#E9C46A",CD8T="#457B9D")) +
            labs(x="Day", y="Fold change vs baseline", color=NULL) +
            theme_classic()
        ggplotly(p)
    })

    # ── Scenario comparison ────────────────────────────────────────────────
    output$compare_plot <- renderPlotly({
        all <- all_results()
        p1 <- ggplot(all, aes(time_day, CRP, color=scenario)) +
            geom_line(linewidth=0.8) +
            geom_hline(yintercept=5, linetype="dashed", color="gray50") +
            labs(x="Day", y="CRP (mg/L)", color=NULL) + theme_classic()
        p2 <- ggplot(all, aes(time_day, IL6, color=scenario)) +
            geom_line(linewidth=0.8) +
            geom_hline(yintercept=7, linetype="dashed", color="gray50") +
            labs(x="Day", y="IL-6 (pg/mL)", color=NULL) + theme_classic()
        subplot(ggplotly(p1) %>% layout(showlegend=TRUE),
                ggplotly(p2) %>% layout(showlegend=FALSE),
                nrows=1, shareX=TRUE, titleY=TRUE)
    })
    output$compare_table <- renderDT({
        all <- all_results()
        tbl <- all %>%
            filter(time_day >= max(time_day) - 1) %>%
            group_by(scenario) %>%
            summarise(
                CRP_final   = round(mean(CRP, na.rm=TRUE), 1),
                IL6_final   = round(mean(IL6, na.rm=TRUE), 1),
                ESR_final   = round(mean(ESR_P, na.rm=TRUE), 1),
                NIH_final   = round(mean(NIH_S, na.rm=TRUE), 2),
                WALL_final  = round(mean(WALL, na.rm=TRUE), 2),
                STEN_pct    = round(mean(STEN*100, na.rm=TRUE), 1),
                Remission   = paste0(round(100*mean(NIH_disc==0, na.rm=TRUE), 0), "%"),
                .groups = "drop"
            )
        datatable(tbl, rownames=FALSE,
                  options=list(dom='t', pageLength=10)) %>%
            formatStyle("CRP_final",
                backgroundColor = styleInterval(c(5, 20),
                    c("#d4edda","#fff3cd","#f8d7da")))
    })

    # ── Risk assessment ────────────────────────────────────────────────────
    output$risk_plot <- renderPlotly({
        d <- sim_result()
        risk_df <- data.frame(
            day  = d$time_day,
            HTN_risk = pmin(100, (d$STEN * 100) * 0.6 + (d$WALL - 2) * 3),
            Stroke_risk = pmin(100, (d$STEN * 100) * 0.25),
            AR_risk = pmin(100, d$STEN * 80 * 0.15)
        ) %>% pivot_longer(-day, names_to="risk", values_to="pct")
        p <- ggplot(risk_df, aes(day, pct, color=risk)) +
            geom_line(linewidth=1) +
            scale_color_manual(
                values=c(HTN_risk="#E63946", Stroke_risk="#264653", AR_risk="#F4A261"),
                labels=c("Hypertension risk","Stroke risk","Aortic regurg risk")) +
            labs(x="Day", y="Cumulative risk (%)", color=NULL,
                 title="Organ Complication Risk Trajectory") +
            theme_classic()
        ggplotly(p)
    })
    output$risk_text <- renderText({
        d <- sim_result()
        crp_f  <- round(tail(d$CRP, 1), 1)
        nih_f  <- round(tail(d$NIH_S, 1), 2)
        sten_f <- round(tail(d$STEN, 1) * 100, 1)
        wall_f <- round(tail(d$WALL, 1), 2)

        risk_level <- if(nih_f < 1 && sten_f < 10) "LOW RISK" else
                      if(nih_f < 2 && sten_f < 30) "MODERATE RISK" else
                      "HIGH RISK"

        paste0(
            "Risk category: ", risk_level, "\n",
            "CRP: ", crp_f, " mg/L\n",
            "NIH Score: ", nih_f, "\n",
            "Stenosis: ", sten_f, "%\n",
            "Wall thickness: ", wall_f, " mm\n\n",
            "Est. 5-yr HTN risk: ",
            round(pmin(95, sten_f * 0.6 + wall_f * 5), 0), "%\n",
            "Est. stroke risk: ",
            round(pmin(50, sten_f * 0.25), 0), "%"
        )
    })
}

shinyApp(ui, server)
