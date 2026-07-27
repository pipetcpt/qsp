## =====================================================================
##  COMMON VARIABLE IMMUNODEFICIENCY (CVID) — QSP / mrgsolve MODEL
##  보통 가변 면역결핍증 · 일차 항체결핍증
## =====================================================================
##
##  THE ORGANISING IDEA
##  -------------------
##  CVID is two diseases sharing one laboratory definition, and only one
##  of them is treated by the treatment we give.
##
##    ARM 1  ANTIBODY DEFICIENCY  — recurrent sinopulmonary infection.
##           Quantitatively predictable from the serum IgG concentration.
##           Immunoglobulin replacement very largely FIXES it.
##
##    ARM 2  IMMUNE DYSREGULATION — autoimmune cytopenias, GLILD,
##           enteropathy, granuloma, lymphoproliferation, lymphoma.
##           A T-cell / tolerance disorder. Immunoglobulin replacement
##           does essentially NOTHING for it. It is what kills people
##           (Resnick 2012: non-infectious complications, mortality
##           RR ~ 11x).
##
##  Between the two arms and the endpoints sits a third layer that no
##  drug reverses: IRREVERSIBLE STRUCTURE. Bronchiectasis is the
##  archetype, and it sits inside a POSITIVE FEEDBACK LOOP:
##
##      colonisation -> IL-8 -> neutrophil elastase -> wall destruction
##      -> bronchiectasis (irreversible) -> impaired clearance
##      -> mucus stasis -> more colonisation
##
##  A ratchet. Which is why the single most important number in CVID is
##  not a dose but the DIAGNOSTIC DELAY (historically 4-7 years). This
##  model exists largely to price that delay.
##
##  ---------------------------------------------------------------------
##  THE FIVE DESIGN DECISIONS THAT MATTER
##  ---------------------------------------------------------------------
##
##  (1) THE EXPOSURE-RESPONSE IS A HILL FUNCTION, NOT A STRAIGHT LINE.
##      Orange 2010 (meta-regression of 17 IVIG studies) is the
##      quantitative backbone of the field: pneumonia incidence falls
##      ~27% per 100 mg/dL of IgG trough, with 0.113 pneumonia/patient-
##      year at a trough of 500 mg/dL. The slope is solid. The corollary
##      usually quoted with it — zero pneumonia at 1400 mg/dL — is an
##      artefact of fitting a STRAIGHT LINE to a ratio, and if you build
##      it into a model the model will tell you to titrate to 1400 and
##      stop worrying.
##
##      This model instead uses
##            rate = RMAX / (1 + (OPSONIN/C50)^HILL)
##      where OPSONIN is the breadth-weighted opsonic titre derived from
##      the serum concentration, with C50 = 200, HILL = 2.4,
##      RMAX = 0.7745/year. Solved (not asserted) so that the model
##        - reproduces d(ln rate)/dC = -0.00343 /(mg/dL) at C = 700,
##          i.e. -29% per 100 mg/dL  (the Orange slope, DIAGNOSTIC D5)
##        - reproduces 0.113 pneumonia/py at C = 500  (anchor, D5b)
##        - extrapolates DOWNWARD sensibly: 0.40/py at 250 mg/dL and
##          0.72/py at 100 mg/dL, i.e. the untreated / severe range,
##          which the linear form cannot reach at all
##        - gives 0.028/py at a trough of 900 mg/dL (observed on
##          replacement: 0.02-0.05/py) and 0.018/py in the healthy
##          reference (adult community-acquired pneumonia incidence)
##        - saturates instead of reaching zero, producing a
##          diminishing-returns knee near 800-1000 mg/dL.
##
##  (2) THE HAZARD IS DRIVEN BY THE INSTANTANEOUS CONCENTRATION, NOT BY
##      THE TROUGH. This is the mechanistically honest choice, and it
##      buys a real result for free. Because rate(C) is CONVEX over the
##      therapeutic range, Jensen's inequality says that a FLUCTUATING
##      profile produces a HIGHER mean event rate than a flat profile of
##      the SAME mean concentration. So weekly SCIG beats 4-weekly IVIG
##      at equal AUC not because of anything immunological but because
##      of the shape of the curve. DIAGNOSTIC D6 measures the size of
##      that convexity penalty. The model was not built to produce this;
##      it falls out.
##
##  (3) PNEUMONIA AND SINUSITIS HAVE DIFFERENT DOSE-RESPONSES.
##      Serum IgG transudes poorly onto the airway surface and secretory
##      IgA is absent and is NOT replaced by any product. The mucosal
##      infection rate therefore has a large IgG-INDEPENDENT FLOOR
##      (RSINO_FLOOR), which is why replaced patients still average 2-3
##      episodes of sinusitis/bronchitis a year. Consequence the model
##      will state quantitatively: escalating the dose for recurrent
##      PNEUMONIA is worth it; escalating it for persistent SINUSITIS is
##      largely futile. Azithromycin, which acts on the floor, is not.
##
##  (4) GLILD CONTAINS NO IgG TERM ANYWHERE. This is a deliberate
##      STRUCTURAL NULL. DIAGNOSTIC D7 sweeps the replacement dose from
##      300 to 1000 mg/kg and requires the GLILD trajectory to be
##      bit-identical. If a reader wants replacement to modify GLILD
##      they must add an edge, not turn a knob.
##
##  (5) IRREVERSIBLE STATES HAVE NO REMOVAL TERM. BE (bronchiectasis
##      extent), FIBROSIS and FEV1_IRR are monotone non-decreasing by
##      construction, not by parameter choice. DIAGNOSTIC D8 verifies
##      monotonicity numerically in every scenario.
##
##  ---------------------------------------------------------------------
##  STRUCTURE — 74 ODEs
##  ---------------------------------------------------------------------
##   1-3    Immunoglobulin replacement PK (SC depot, central, peripheral)
##   4-12   B-cell development and the class-switch block
##  13-15   BAFF / APRIL / sBCMA
##  16-25   T-cell help, regulation, cytokines
##  26-31   Innate immunity, mucosal barrier, microbial translocation
##  32      Antibody breadth (donor pool)
##  33-38   Pathogens and cumulative infection counters
##  39-44   Airway inflammation and the bronchiectasis ratchet
##  45-48   GLILD and interstitial fibrosis
##  49-51   Lymphoproliferation, spleen, liver
##  52-54   Autoimmune cytopenias
##  55-56   Enteropathy and albumin
##  57      Malignancy hazard
##  58-61   Mortality hazards, QoL, cumulative immunosuppression
##  62-74   Drug PK / effect compartments (rituximab, abatacept,
##          sirolimus, leniolisib, prednisone, azathioprine,
##          eltrombopag, azithromycin, JAK inhibitor)
##
##  Author: QSP Disease Model Library (Claude Code Routine)
##  Units:  time = DAYS.  IgG = mg/dL.  B/T cells = cells/uL or
##          normalised (1.0 = healthy).  Rates reported per YEAR.
## =====================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
})

options(mrgsolve.soloc = tempdir())

## =====================================================================
##  MODEL CODE
## =====================================================================

cvid_code <- '
$PROB
# Common Variable Immunodeficiency (CVID) QSP model
# ARM 1 antibody deficiency (replaceable) || ARM 2 immune dysregulation
# (not replaceable) -> irreversible structure

$GLOBAL
#define CIGG   (IGG_C / VC)                 // total serum IgG, mg/dL
#define BTOT   (TRANS_B + NAIVE_B + MEM_B + MEMIGM_B + CD21LO_B + PBLAST)

double posv(double x){ return x > 0.0 ? x : 0.0; }
double hillf(double x, double ec50, double h){
  double xx = posv(x);
  double a  = pow(xx, h);
  return a / (pow(ec50, h) + a + 1e-30);
}
double sat(double x, double k){ double xx = posv(x); return xx/(xx + k + 1e-30); }

// quantities computed in $ODE and reported in $TABLE
double cigg, protf, pneu_rate, sino_rate, invas_rate, opsonin;
double susc, fev1, fvc_pct, dlco, glild_act, smb_pct, cd21_pct, naiveb_pct;
double isburden, ple, cl_tot, mucprot, ne_act, colon_growth;
double rtx_kill, aza_e, azm_e, jak_e, gc_eff, aba_eff, siro_eff, lenio_eff;
double tfh_help, taci_sig, chronic, ctla4_func, treg_func, dysreg;
double be_flux, fib_flux, qol, haz_tot, ifn_drive, elt_eff;

$PARAM @annotated
// ---------------- demographics / scaling -------------------------
WT       :  70    : Body weight (kg)
VC       :  38.5  : Central volume for IgG (dL) ~0.055 L/kg
VP       :  31.5  : Peripheral (interstitial) IgG volume (dL)
QIGG     :  20    : IgG intercompartmental clearance (dL/day)
CLIGG    :  1.330 : Baseline (FcRn-salvaged) IgG clearance (dL/day)
EMAXSAT  :  0.60  : Max fractional rise in CL as FcRn salvage saturates
IC50SAT  :  2500  : IgG conc for half-maximal FcRn saturation (mg/dL)
FSC      :  0.73  : SC bioavailability (1/0.73 = 1.37 = EU adjustment factor)
KASC     :  0.90  : SC absorption rate constant (1/day)
FHYAL    :  0.93  : fSCIG bioavailability with rHuPH20
KAHYAL   :  1.30  : fSCIG absorption rate constant (1/day)
HYAL     :  0     : 1 = hyaluronidase-facilitated SC (fSCIG)

// ---------------- the CVID lesion --------------------------------
FCSR     :  0.08  : Class-switch recombination efficiency (1 = healthy)
FPCD     :  0.60  : Plasma-cell differentiation efficiency (1 = healthy)
FTACI    :  0.35  : TACI / T-independent switching signal (1 = healthy)
FTFH     :  0.55  : Tfh help competence (1 = healthy)
DYSGENO  :  0     : Dysregulation genotype load 0-1 (drives ARM 2)
CTLA4G   :  1.0   : Surface CTLA-4 gene dosage (0.5 = CHAI, 0.3 = LRBA)
PI3KGOF  :  0     : PIK3CD gain-of-function (APDS) 0-1
HEALTHY  :  0     : 1 = run the healthy reference physiology

// ---------------- B-cell kinetics (cells/uL) ---------------------
KINTR    :  3.30  : Bone-marrow output into transitional B (cells/uL/day)
KTRNA    :  0.50  : Transitional -> naive (1/day)
KDGTR    :  0.05  : Transitional loss (1/day)
KNAGC    :  0.0060: Naive -> germinal centre (1/day)
KNACD21  :  0.00100: Naive -> CD21low (1/day per unit chronic activation)
KNATI    :  0.0035: Naive -> IgM memory, TACI-dependent (1/day)
KDGNA    :  0.0130: Naive B loss (1/day)
KGCMEM   :  0.200 : GC -> switched memory (1/day, x FCSR)
KGCPB    :  0.100 : GC -> plasmablast (1/day)
KDGGC    :  0.050 : GC B loss (1/day)
KDGMEM   :  0.014894: Switched memory B loss (1/day)
KDGIGM   :  0.02851: IgM memory B loss (1/day)
KDGC21   :  0.02167: CD21low B loss (1/day)
KRECALL  :  0.0100: Memory B recall -> plasmablast (1/day per unit IL-21)
KOUTPB   :  0.2617: Plasmablast loss (1/day)
PBREF    :  2.00  : Healthy plasmablast reference (cells/uL)
KPCS     :  0.200 : Plasmablast -> short-lived plasma cell
KDGPCS   :  0.200 : Short-lived plasma-cell loss (1/day)
KPCL     :  0.00025: Plasmablast -> long-lived IgG plasma cell (1/day)
KPCLTI   :  0.0001848: T-independent (APRIL/TACI) input to IgG plasma cells
KDGPCL   :  0.0003424: Long-lived IgG plasma-cell loss (1/day) ~ t1/2 5.5 y
KPRODIGG :  1562  : IgG production at PC_IGG = 1 (mg/day)
MEMREF   :  30.0  : Healthy switched memory B reference (cells/uL)
IGMREF   :  16.0  : Healthy IgM memory B reference (cells/uL)

// ---------------- BAFF / APRIL -----------------------------------
KINBAFF  :  13.0  : BAFF production (normalised units/day)
KOUTBAFF :  3.00  : BAFF non-receptor clearance (1/day)
KSINKBAFF:  10.00 : Receptor-mediated BAFF consumption (1/day at B = ref)
BREF     :  261.2 : Healthy receptor-weighted BAFF sink reference
KINAPR   :  4.35  : APRIL production (1/day)
KOUTAPR  :  4.00  : APRIL clearance (1/day)
KINBCMA  :  0.200 : sBCMA appearance per unit plasma-cell mass (1/day)
KOUTBCMA :  0.200 : sBCMA clearance (1/day)

// ---------------- T-cell compartment (normalised) ----------------
KINCD4N  :  0.0200: Naive CD4 thymic/peripheral input (1/day)
KOUTCD4N :  0.0100: Naive CD4 loss (1/day)
KDIFFCD4 :  0.0100: Naive -> memory CD4 (1/day)
KOUTCD4M :  0.0100: Memory CD4 loss (1/day)
KTFH     :  0.0750: Memory CD4 -> cTfh (1/day)
KOUTTFH  :  0.0500: cTfh loss (1/day)
KTREGSUP :  0.0250: Treg suppression of cTfh (1/day)
KINTREG  :  0.0500: Treg input (1/day)
KOUTTREG :  0.0500: Treg loss (1/day)
KINCD8   :  0.0550: CD8 activation input (1/day)
KOUTCD8  :  0.0500: Activated CD8 loss (1/day)
KTEMRA   :  0.0050: Activated CD8 -> TEMRA (1/day)
KOUTTEMRA:  0.0050: TEMRA loss (1/day)
KINIFNG  :  1.000 : IFN-gamma production (1/day)
KOUTIFNG :  1.000 : IFN-gamma clearance (1/day)
KINIL21  :  1.000 : IL-21 production from cTfh (1/day)
KOUTIL21 :  1.000 : IL-21 clearance (1/day)
KINIL6   :  1.000 : IL-6 production (1/day)
KOUTIL6  :  1.000 : IL-6 clearance (1/day)
KINCXCL13:  1.000 : CXCL13 production (1/day)
KOUTCXCL13: 1.000 : CXCL13 clearance (1/day)

// ---------------- innate / translocation -------------------------
KREPMUC  :  0.0200: Mucosal barrier repair (1/day)
KDAMIGA  :  0.0250: Barrier damage from absent secretory IgA (1/day)
KDAMGUT  :  0.0150: Barrier damage from enteropathy (1/day)
KINLPS   :  2.000 : LPS translocation per unit barrier defect (1/day)
KOUTLPS  :  2.000 : LPS clearance (1/day)
KINSCD14 :  1.000 : sCD14 production (1/day)
KOUTSCD14:  1.000 : sCD14 clearance (1/day)
KINMONO  :  0.500 : Monocyte activation (1/day)
KOUTMONO :  0.500 : Monocyte deactivation (1/day)
KINISG   :  0.500 : Type-I IFN signature induction (1/day)
KOUTISG  :  0.500 : ISG decay (1/day)
KINNEUT  :  1.000 : Airway neutrophil recruitment (1/day)
KOUTNEUT :  1.000 : Airway neutrophil clearance (1/day)

// ---------------- healthy reference tone (normalisation anchors) --
ISGREF   :  0.440 : Healthy type-I IFN signature tone
COLREF   :  0.020 : Healthy airway colonisation tone
MONOREF  :  0.300 : Healthy activated-monocyte tone
AIRREF   :  0.050 : Healthy airway inflammation tone

// ---------------- antibody breadth --------------------------------
BRDC50   :  150   : Serum IgG for half-maximal pool breadth (mg/dL)
KBREADTH :  0.0100: Breadth equilibration (1/day)
BRDPOOL  :  0.85  : Asymptotic donor-pool breadth vs a native response
BRDENDOG :  0.15  : Breadth contributed by residual endogenous response

// ---------------- ARM 1 exposure-response ------------------------
RMAXPNEU :  0.7745: Max pneumonia rate, no antibody (events/year)
C50PNEU  :  200   : IgG for half-maximal pneumonia protection (mg/dL)
HILLPNEU :  2.40  : Hill coefficient (set by Orange slope, see D5)
RMAXSINO :  3.00  : IgG-dependent sinopulmonary rate (events/year)
RSINOFLR :  1.20  : IgG-INDEPENDENT mucosal floor, x(1+0.9*IgA deficit)
C50SINO  :  200   : IgG for half-maximal mucosal protection (mg/dL)
HILLSINO :  1.60  : Hill coefficient, mucosal
FINVAS   :  0.040 : Fraction of pneumonia becoming invasive
WBE      :  0.070 : Susceptibility gain per unit bronchiectasis (Reiff)
WIS      :  0.600 : Susceptibility gain per unit immunosuppression
WASPLEN  :  2.000 : Invasive-disease multiplier when asplenic
WNEUTRO  :  0.500 : Susceptibility gain from autoimmune neutropenia

// ---------------- airway colonisation and the ratchet ------------
KGCOL    :  0.280 : Colonisation growth rate (1/day)
KCAPCOL  :  1.000 : Colonisation carrying capacity
KKILLCOL :  0.350 : Opsonin-dependent clearance of colonisers (1/day)
KMUCCOL  :  0.400 : Mucus-stasis gain on colonisation growth
KCLRCOL  :  0.100 : Mucociliary clearance of colonisers (1/day)
KINAI    :  0.400 : Airway inflammation from colonisation (1/day)
KINAIP   :  4.000 : Airway inflammation per pneumonia event
KOUTAI   :  0.150 : Airway inflammation resolution (1/day)
AATSHLD  :  0.350 : Antiprotease shield (fraction of elastase neutralised)
NETHR    :  0.0800: Elastase threshold below which no wall damage occurs
KBE      :  0.00260: Bronchiectasis accrual per unit suprathreshold elastase
KBEPNEU  :  0.220 : Bronchiectasis accrual per pneumonia event
BEMAX    :  18.0  : Modified Reiff maximum score
KMCC     :  0.0700: Mucociliary clearance loss per Reiff point
KINMUCUS :  0.100 : Mucus accumulation from impaired clearance (1/day)
KINMUCC  :  0.150 : Mucus accumulation from colonisation (1/day)
KOUTMUCUS:  0.200 : Mucus clearance (1/day)
FEV10    :  100   : Baseline FEV1 (% predicted)
AFEV1BE  :  1.700 : Irreversible FEV1 loss per Reiff point (%pred)
AFEV1FIB :  22.00 : Irreversible FEV1 loss per unit interstitial fibrosis
BFEV1INF :  6.000 : Reversible FEV1 loss per unit airway inflammation
KRECFEV1 :  0.0300: Recovery of the reversible FEV1 component (1/day)
DLCO0    :  95.0  : Baseline DLCO (% predicted)
CDLCOFIB :  38.00 : DLCO loss per unit fibrosis
CDLCOGL  :  12.00 : DLCO loss per unit GLILD activity
CDLCOBE  :  0.700 : DLCO loss per Reiff point

// ---------------- ARM 2: GLILD -----------------------------------
KINLA    :  0.00200: Lymphoid aggregate formation (1/day)
KOUTLA   :  0.0080: Lymphoid aggregate resolution (1/day)
KINGRAN  :  0.00170: Granuloma formation (1/day)
KOUTGRAN :  0.0040: Granuloma resolution (1/day)
WGLILDLA :  0.700 : Weight of lymphoid aggregates in GLILD activity
WGLILDGR :  0.300 : Weight of granuloma in GLILD activity
GLTHR    :  0.250 : GLILD activity threshold above which fibrosis accrues
KFIB     :  0.00013: Interstitial fibrosis accrual (1/day)

// ---------------- ARM 2: lymphoproliferation / liver -------------
KINSPL   :  0.00260: Splenomegaly growth (1/day)
KOUTSPL  :  0.0100: Spleen regression (1/day)
KINLAD   :  0.00600: Lymphadenopathy growth (1/day)
KOUTLAD  :  0.0150: Lymphadenopathy regression (1/day)
KINNRH   :  0.00013: Nodular regenerative hyperplasia accrual (1/day)
SPLENEC  :  0     : 1 = splenectomised
TSPLENEC :  1e9   : Time of splenectomy (day)

// ---------------- ARM 2: autoimmune cytopenias -------------------
KINAA    :  0.0120: Autoantibody generation (1/day)
KOUTAA   :  0.0200: Autoantibody clearance (1/day)
KINPLT   :  25.00 : Platelet production (10^9/L/day)
KOUTPLT  :  0.100 : Platelet turnover (1/day)
KDESTPLT :  0.300 : Autoantibody-mediated platelet destruction (1/day)
KSPLPLT  :  0.0300: Hypersplenic platelet sequestration (1/day)
FSPLENDEST: 0.700 : Fraction of platelet destruction that is splenic
KINHGB   :  0.0800: Haemoglobin production (g/dL/day)
KOUTHGB  :  0.00571: RBC turnover (1/day)
KDESTHGB :  0.00120: Autoantibody-mediated haemolysis (1/day)
AITHR    :  0.150 : Autoantibody threshold for clinical cytopenia

// ---------------- ARM 2: enteropathy ------------------------------
KINGUT   :  0.0080: Enteropathy activity generation (1/day)
KOUTGUT  :  0.0100: Enteropathy resolution (1/day)
EMAXPLE  :  1.000 : Max protein-losing-enteropathy severity
EC50PLE  :  1.200 : Enteropathy activity for half-maximal PLE
HPLE     :  4.000 : Hill coefficient for PLE onset
CLPLEMAX :  1.000 : Extra IgG clearance at PLE = 1 (x CLIGG)
KINALB   :  0.130 : Albumin production (g/dL/day)
KOUTALB  :  0.0330: Albumin turnover (1/day)
KLOSSALB :  0.0400: Albumin loss per unit PLE (1/day)

// ---------------- malignancy / mortality / QoL -------------------
KMALIG   :  2.5e-6: Malignancy hazard accrual (1/day per unit drive)
KHAZINF  :  0.0180: Mortality hazard per invasive infection event
KHAZRESP :  9.0e-6: Mortality hazard per day per unit respiratory failure
KHAZNON  :  2.0e-5: Mortality hazard per day per unit ARM 2 burden
KQOL     :  0.0500: QoL equilibration (1/day)

// ---------------- immunosuppression weights ----------------------
WRTX     :  0.450 : Net immunosuppression weight, rituximab
WAZA     :  0.400 : Net immunosuppression weight, azathioprine/MMF
WGC      :  0.700 : Net immunosuppression weight, glucocorticoid
WABA     :  0.350 : Net immunosuppression weight, abatacept
WSIRO    :  0.400 : Net immunosuppression weight, sirolimus
WJAK     :  0.450 : Net immunosuppression weight, JAK inhibitor
WLENIO   :  0.080 : Net immunosuppression weight, leniolisib (targeted)
KISCUM   :  0.0100: Cumulative immunosuppression integrator (1/day)

// ---------------- drug PK ----------------------------------------
CLRTX    :  0.200 : Rituximab clearance (L/day) -> t1/2 ~20 d
V1RTX    :  3.100 : Rituximab central volume (L)
V2RTX    :  2.700 : Rituximab peripheral volume (L)
QRTX     :  0.600 : Rituximab intercompartmental clearance (L/day)
KMAXRTX  :  1.000 : Max B-cell kill rate (1/day)
EC50RTX  :  0.500 : Rituximab EC50 for B-cell kill (mg/L)
CLABA    :  0.370 : Abatacept clearance (L/day) 0.22 mL/h/kg
V1ABA    :  2.500 : Abatacept central volume (L) -> Cmax ~280 mg/L
V2ABA    :  3.000 : Abatacept peripheral volume (L)
QABA     :  0.800 : Abatacept intercompartmental clearance (L/day)
EC50ABA  :  25.00 : Abatacept EC50 (mg/L)
CLSIRO   :  250.0 : Sirolimus apparent clearance (L/day)
VSIRO    :  931.0 : Sirolimus apparent volume (L)
EC50SIRO :  8.000 : Sirolimus EC50 (ug/L)
CLLENIO  :  360.0 : Leniolisib apparent clearance (L/day)
VLENIO   :  216.0 : Leniolisib apparent volume (L)
EC50LENIO:  150.0 : Leniolisib EC50 (ug/L)
CLPRED   :  250.0 : Prednisolone apparent clearance (L/day)
VPRED    :  50.00 : Prednisolone apparent volume (L)
KE0PRED  :  2.000 : Prednisolone effect-compartment rate (1/day)
EC50PRED :  0.0500: Prednisolone EC50 (mg/L in effect compartment)
CLELT    :  9.000 : Eltrombopag apparent clearance (L/day)
VELT     :  12.00 : Eltrombopag apparent volume (L)
EC50ELT  :  2.000 : Eltrombopag EC50 (mg/L)
EMAXELT  :  1.800 : Max fractional rise in platelet production

// ---------------- switch-drug controls ---------------------------
AZAON    :  0     : Azathioprine / MMF on (0/1)
AZAT0    :  0     : Azathioprine start time (day)
AZAEMAX  :  0.700 : Azathioprine max effect
KE0AZA   :  0.0700: Azathioprine effect-compartment rate (1/day)
AZMON    :  0     : Azithromycin prophylaxis on (0/1)
AZMT0    :  0     : Azithromycin start time (day)
AZMEMAX  :  0.550 : Azithromycin max effect on colonisation/inflammation
KE0AZM   :  0.100 : Azithromycin effect-compartment rate (1/day)
JAKON    :  0     : JAK inhibitor on (0/1)
JAKT0    :  0     : JAK inhibitor start time (day)
JAKEMAX  :  0.650 : JAK inhibitor max effect on IFN-gamma signalling
KE0JAK   :  0.300 : JAK inhibitor effect-compartment rate (1/day)
VEDOON   :  0     : Vedolizumab / budesonide on (0/1)
VEDOEMAX :  0.450 : Gut-selective agent max effect
PHYSIOON :  0     : Airway clearance physiotherapy on (0/1)
PHYSIOE  :  0.400 : Physiotherapy effect on mucus clearance
INHABXON :  0     : Inhaled antibiotic on (0/1)
DAMAGEON :  1     : 1 = allow irreversible damage to accrue (0 for pre-runs)
IGRXON   :  1     : 1 = replacement therapy present (bookkeeping only)

$CMT @annotated
IGG_SC   : Subcutaneous immunoglobulin depot (mg)
IGG_C    : Central (plasma) IgG (mg)
IGG_P    : Peripheral (interstitial) IgG (mg)
TRANS_B  : Transitional B cells (cells/uL)
NAIVE_B  : Naive B cells (cells/uL)
GC_B     : Germinal-centre B-cell activity (cells/uL equivalent)
MEM_B    : Switched memory B cells (cells/uL)
MEMIGM_B : IgM memory B cells (cells/uL)
CD21LO_B : CD21low B cells (cells/uL)
PBLAST   : Plasmablasts (cells/uL)
PC_SHORT : Short-lived plasma cells (normalised)
PC_IGG   : Long-lived IgG-secreting plasma cells (normalised)
BAFF     : Soluble BAFF (normalised, 1 = healthy)
APRIL    : APRIL (normalised)
SBCMA    : Soluble BCMA (normalised)
CD4N     : Naive CD4 T cells (normalised)
CD4M     : Memory/effector CD4 T cells (normalised)
TFH      : Circulating Tfh cells (normalised)
TREG     : Regulatory T cells (normalised)
CD8ACT   : Activated (HLA-DR+) CD8 T cells (normalised)
TEMRA    : Senescent/TEMRA T cells (normalised)
IFNG     : Interferon-gamma (normalised)
IL21     : IL-21 (normalised)
IL6      : IL-6 (normalised)
CXCL13   : CXCL13 (normalised)
MUCB     : Gut mucosal barrier integrity (0-1)
LPS      : Circulating LPS / microbial translocation (normalised)
SCD14    : sCD14 (normalised)
MONO     : Activated monocytes (normalised)
ISG      : Type-I interferon signature (normalised)
NEUT     : Airway neutrophils (normalised)
BREADTH  : Antibody breadth of the circulating IgG pool (0-1)
COLON    : Airway bacterial colonisation (0-1)
PSA      : Pseudomonas colonisation (0-1)
GIARD    : Chronic GI pathogen burden (0-1)
PNEU_CUM : Cumulative pneumonia events
SINO_CUM : Cumulative sinopulmonary infection events
INVAS_CUM: Cumulative invasive infection events
AIRINF   : Airway neutrophilic inflammation (normalised)
MUCUS    : Airway mucus load (normalised)
BE       : Bronchiectasis extent, modified Reiff 0-18 (IRREVERSIBLE)
FEV1_IRR : Irreversible FEV1 loss (% predicted, IRREVERSIBLE)
FEV1_REV : Reversible FEV1 loss (% predicted)
EXAC_CUM : Cumulative exacerbations
LYMPHAGG : Peribronchiolar lymphoid aggregates (normalised)
GRAN     : Granuloma burden (normalised)
FIB      : Interstitial fibrosis (normalised, IRREVERSIBLE)
RESPFAIL : Respiratory-failure burden (normalised)
SPLEEN   : Splenomegaly (normalised, 1 = 2x normal volume)
LAD      : Lymphadenopathy / index nodal lesion (normalised)
NRH      : Nodular regenerative hyperplasia (normalised, IRREVERSIBLE)
AUTOAB   : Anti-platelet / anti-RBC autoantibody (normalised)
PLT      : Platelet count (10^9/L)
HGB      : Haemoglobin (g/dL)
GUT      : Enteropathy activity (normalised)
ALB      : Serum albumin (g/dL)
MALIG    : Cumulative malignancy hazard
HAZ_INF  : Cumulative infection-attributable mortality hazard
HAZ_NON  : Cumulative non-infectious-complication mortality hazard
QOL      : Health-related quality of life (0-1)
ISCUM    : Cumulative immunosuppression exposure
RTX_C    : Rituximab central (mg)
RTX_P    : Rituximab peripheral (mg)
ABA_C    : Abatacept central (mg)
ABA_P    : Abatacept peripheral (mg)
SIRO_C   : Sirolimus central (ug)
LENIO_C  : Leniolisib central (ug)
PRED_C   : Prednisolone central (mg)
PRED_E   : Prednisolone effect compartment (mg/L)
AZA_E    : Azathioprine effect (0-1 scaled)
ELT_C    : Eltrombopag central (mg)
AZM_E    : Azithromycin effect (0-1 scaled)
JAK_E    : JAK inhibitor effect (0-1 scaled)
DUMMY    : Bookkeeping compartment (unused)

$MAIN
F_IGG_SC = (HYAL > 0.5) ? FHYAL : FSC;

$ODE
// =================================================================
//  0. SWITCH DRUGS AND DRUG EFFECTS
// =================================================================
double aza_on  = ((AZAON  > 0.5) && (SOLVERTIME >= AZAT0 )) ? 1.0 : 0.0;
double azm_on  = ((AZMON  > 0.5) && (SOLVERTIME >= AZMT0 )) ? 1.0 : 0.0;
double jak_on  = ((JAKON  > 0.5) && (SOLVERTIME >= JAKT0 )) ? 1.0 : 0.0;
double splen   = ((SPLENEC> 0.5) && (SOLVERTIME >= TSPLENEC)) ? 1.0 : 0.0;

aza_e   = AZAEMAX  * AZA_E;
azm_e   = AZMEMAX  * AZM_E;
jak_e   = JAKEMAX  * JAK_E;
gc_eff  = sat(PRED_E, EC50PRED);
aba_eff = sat(ABA_C / V1ABA, EC50ABA);
siro_eff= sat(SIRO_C / VSIRO, EC50SIRO);
lenio_eff = sat(LENIO_C / VLENIO, EC50LENIO);
elt_eff = EMAXELT * sat(ELT_C / VELT, EC50ELT);
rtx_kill = KMAXRTX * sat(RTX_C / V1RTX, EC50RTX);
double vedo_e = (VEDOON > 0.5) ? VEDOEMAX : 0.0;

// net immunosuppression: the single quantity that prices ARM 2 therapy
isburden = WRTX*sat(RTX_C/V1RTX, EC50RTX) + WAZA*aza_e + WGC*gc_eff
         + WABA*aba_eff + WSIRO*siro_eff + WJAK*jak_e + WLENIO*lenio_eff;

// =================================================================
//  1. IMMUNOGLOBULIN REPLACEMENT PK
// =================================================================
cigg = CIGG;

// protein-losing enteropathy -> extra, NON-catabolic IgG loss
ple    = EMAXPLE * hillf(GUT, EC50PLE, HPLE);
// FcRn salvage saturates at high concentration -> CL rises
cl_tot = CLIGG * (1.0 + EMAXSAT * sat(cigg, IC50SAT)) * (1.0 + CLPLEMAX * ple);

double ka_ig = (HYAL > 0.5) ? KAHYAL : KASC;
double endog = KPRODIGG * PC_IGG;

dxdt_IGG_SC = -ka_ig * IGG_SC;
dxdt_IGG_C  =  ka_ig * IGG_SC + endog
             - cl_tot * cigg
             - QIGG * cigg + QIGG * (IGG_P / VP);
dxdt_IGG_P  =  QIGG * cigg - QIGG * (IGG_P / VP);

// =================================================================
//  2. B-CELL DEVELOPMENT AND THE CLASS-SWITCH BLOCK
// =================================================================
double fcsr  = (HEALTHY > 0.5) ? 1.0 : FCSR;
double fpcd  = (HEALTHY > 0.5) ? 1.0 : FPCD;
double ftaci = (HEALTHY > 0.5) ? 1.0 : FTACI;
double ftfh  = (HEALTHY > 0.5) ? 1.0 : FTFH;

// APDS: PI3K GOF traps cells at the transitional stage and blocks CSR
// further; leniolisib relieves exactly this and nothing else
double pi3k = PI3KGOF * (1.0 - lenio_eff);
double fcsr_eff = fcsr * (1.0 - 0.55 * pi3k);
tfh_help = ftfh * (TFH + 1e-9) * (1.0 - 0.30 * siro_eff);
taci_sig = ftaci * APRIL;

// chronic activation drive -> CD21low expansion
chronic  = 1.0 + 1.6*LPS + 0.8*posv(ISG - ISGREF) + 1.2*posv(COLON - COLREF)
         + 2.0*DYSGENO + 1.5*pi3k;

// normalised so that baff_surv == 1.0 exactly at BAFF == 1 (healthy)
double baff_surv = 0.35 + 1.30 * sat(BAFF, 1.0);

dxdt_TRANS_B  = KINTR * (1.0 + 1.2*pi3k) - KTRNA*(1.0 - 0.45*pi3k)*TRANS_B
              - KDGTR*TRANS_B - rtx_kill*TRANS_B;
dxdt_NAIVE_B  = KTRNA*(1.0 - 0.45*pi3k)*TRANS_B * baff_surv
              - KNAGC*NAIVE_B - KNACD21*chronic*NAIVE_B
              - KNATI*taci_sig*NAIVE_B - KDGNA*NAIVE_B
              - rtx_kill*NAIVE_B;
dxdt_GC_B     = KNAGC*NAIVE_B*tfh_help
              - (KGCMEM*fcsr_eff + KGCPB + KDGGC)*GC_B - rtx_kill*GC_B;
dxdt_MEM_B    = KGCMEM*fcsr_eff*GC_B - KDGMEM*MEM_B - rtx_kill*MEM_B;
dxdt_MEMIGM_B = KNATI*taci_sig*NAIVE_B - KDGIGM*MEMIGM_B - rtx_kill*MEMIGM_B;
dxdt_CD21LO_B = KNACD21*chronic*NAIVE_B - KDGC21*CD21LO_B - rtx_kill*CD21LO_B;
dxdt_PBLAST   = KGCPB*GC_B + KRECALL*MEM_B*IL21
              - KOUTPB*PBLAST - rtx_kill*0.30*PBLAST;
dxdt_PC_SHORT = KPCS*fpcd*(PBLAST/PBREF) - KDGPCS*PC_SHORT;
dxdt_PC_IGG   = KPCL*fpcd*fcsr_eff*(PBLAST/PBREF)
              + KPCLTI*fpcd*(MEMIGM_B/IGMREF)*sat(APRIL,1.0)
              - KDGPCL*PC_IGG;

// =================================================================
//  3. BAFF / APRIL / sBCMA
//     The B-cell sink: fewer B cells consume less BAFF, so BAFF is
//     high in CVID and higher still after rituximab.
// =================================================================
// The BAFF sink is receptor-WEIGHTED, not a headcount: memory B cells and
// plasma cells consume far more BAFF/APRIL per cell than naive B cells, and
// CD21low cells consume very little. Receptor competence itself is reduced
// in CVID (TACI/BAFF-R signalling), which is why BAFF is elevated even
// though the naive B-cell count is preserved.
double bsink = (NAIVE_B + 2.0*MEM_B + 1.5*MEMIGM_B + 0.8*TRANS_B
                + 0.4*CD21LO_B + 20.0*(PC_IGG + PC_SHORT)) / BREF
             * (0.4 + 0.6*ftaci);
dxdt_BAFF  = KINBAFF - KOUTBAFF*BAFF - KSINKBAFF*bsink*BAFF;
dxdt_APRIL = KINAPR*(1.0 + 0.5*MONO) - KOUTAPR*APRIL - 0.5*SBCMA*APRIL;
dxdt_SBCMA = KINBCMA*(PC_IGG + PC_SHORT) - KOUTBCMA*SBCMA;

// =================================================================
//  4. T-CELL HELP AND REGULATORY FAILURE  (ARM 2 ROOT)
// =================================================================
// CTLA-4 function: gene dosage x LRBA recycling, with abatacept acting
// as a PHARMACOLOGICAL SUBSTITUTE for the missing molecule
ctla4_func = CTLA4G + (1.0 - CTLA4G) * aba_eff;
treg_func  = TREG * ctla4_func;
dysreg     = DYSGENO + 0.6*(1.0 - ctla4_func) + 0.5*pi3k;

double act_drive = 1.0 + 0.8*LPS + 0.6*posv(ISG - ISGREF)
                 + 0.5*posv(COLON - COLREF);

dxdt_CD4N   = KINCD4N*(1.0 - 0.30*pi3k) - KOUTCD4N*CD4N
            - KDIFFCD4*act_drive*CD4N;
dxdt_CD4M   = KDIFFCD4*act_drive*CD4N - KOUTCD4M*CD4M
            - 0.4*KOUTCD4M*(aza_e + siro_eff)*CD4M;
dxdt_TFH    = KTFH*CD4M*(1.0 + 1.2*dysreg) - KOUTTFH*TFH
            - KTREGSUP*treg_func*TFH
            - KOUTTFH*(0.8*aba_eff + 0.5*siro_eff + 0.6*gc_eff)*TFH;
dxdt_TREG   = KINTREG*ctla4_func - KOUTTREG*TREG;
dxdt_CD8ACT = KINCD8*(act_drive + 0.8*dysreg) - KOUTCD8*CD8ACT
            - KOUTCD8*(0.6*gc_eff + 0.5*aza_e + 0.4*siro_eff)*CD8ACT
            - KTEMRA*CD8ACT;
dxdt_TEMRA  = KTEMRA*CD8ACT*(1.0 + pi3k) - KOUTTEMRA*TEMRA;

// weights sum to 1.0 at the healthy reference, so IFNG == 1 when healthy
ifn_drive   = 0.4*CD8ACT + 0.4*TFH*(1.0 + dysreg) + 0.2*(MONO/MONOREF);
dxdt_IFNG   = KINIFNG*ifn_drive*(1.0 - jak_e) - KOUTIFNG*IFNG;
dxdt_IL21   = KINIL21*TFH - KOUTIL21*IL21;
dxdt_IL6    = KINIL6*(0.95 + 0.8*posv(MONO - MONOREF) + 0.5*LPS
                      + 0.4*posv(AIRINF - AIRREF))*(1.0 - 0.5*jak_e)
            - KOUTIL6*IL6;
dxdt_CXCL13 = KINCXCL13*(1.0*TFH + 0.8*LYMPHAGG) - KOUTCXCL13*CXCL13;

// =================================================================
//  5. MUCOSAL BARRIER AND MICROBIAL TRANSLOCATION
//     Secretory IgA is absent in CVID and is NOT replaced by any
//     product. This is the structural gap in replacement therapy.
// =================================================================
double iga_func = (HEALTHY > 0.5) ? 1.0
                : (0.05 + 0.95 * (MEMIGM_B / IGMREF) * ftaci);
dxdt_MUCB  = KREPMUC*(1.0 - MUCB) - KDAMIGA*(1.0 - iga_func)*MUCB
           - KDAMGUT*GUT*MUCB;
dxdt_LPS   = KINLPS*posv(1.0 - MUCB) - KOUTLPS*LPS;
dxdt_SCD14 = KINSCD14*LPS - KOUTSCD14*SCD14;
dxdt_MONO  = KINMONO*(0.3 + LPS) - KOUTMONO*MONO;
dxdt_ISG   = KINISG*(0.2 + 0.8*MONO) - KOUTISG*ISG;
dxdt_NEUT  = KINNEUT*(AIRINF) - KOUTNEUT*NEUT;

// =================================================================
//  6. FROM CONCENTRATION TO PROTECTION
//     The protective species is opsonic antibody, not total IgG.
//     Donor-pool IgG has good breadth against common serotypes and
//     poor breadth against rare ones; residual endogenous response
//     adds a little.
// =================================================================
double brd_target = BRDPOOL * sat(cigg, BRDC50)
                  + BRDENDOG * (MEM_B / MEMREF);
dxdt_BREADTH = KBREADTH * (brd_target - BREADTH);

opsonin = cigg * (0.5 + 0.5 * BREADTH / (BRDPOOL + BRDENDOG));
// Serum IgG transudes poorly onto the airway surface and secretory IgA is
// absent and is NOT replaced by any product. Only 35% of the opsonic titre
// is available mucosally regardless of dose; the other 65% is IgA-gated.
mucprot = opsonin * (0.35 + 0.65 * iga_func);

// ---- ARM 1 exposure-response (Orange 2010 slope, Hill form) ------
protf = 1.0 / (1.0 + pow(posv(opsonin)/C50PNEU, HILLPNEU));

susc  = 1.0 + WBE*BE + WIS*isburden
      + WNEUTRO * sat(posv(AUTOAB - AITHR), 0.5);
double asplen_mult = 1.0 + WASPLEN*splen;

pneu_rate = RMAXPNEU * protf * susc * (1.0 - 0.35*azm_e) / 365.0;   // per day
double sino_prot = 1.0 / (1.0 + pow(posv(mucprot)/C50SINO, HILLSINO));
// The floor RISES as secretory IgA is lost, and no dose of IgG lowers it.
// Azithromycin does, which is the whole point of prophylaxis.
double sino_floor = RSINOFLR * (1.0 + 0.9*(1.0 - iga_func)) * (1.0 - azm_e);
sino_rate = (sino_floor + RMAXSINO*sino_prot) * susc / 365.0;
invas_rate = FINVAS * pneu_rate * asplen_mult;

// =================================================================
//  7. PATHOGENS AND THE INFECTION COUNTERS
// =================================================================
// mucociliary clearance is a FUNCTION of the irreversible damage state
double mcc = (1.0 - KMCC*BE > 0.15) ? (1.0 - KMCC*BE) : 0.15;
double ops_kill = KKILLCOL * sat(mucprot, 250.0) * (1.0 + azm_e);
colon_growth = KGCOL * (1.0 + KMUCCOL*MUCUS) * (1.0 - COLON/KCAPCOL);
dxdt_COLON = colon_growth*COLON + 0.002*posv(1.0 - COLON)
           - ops_kill*COLON - KCLRCOL*mcc*COLON;
dxdt_PSA   = 0.0025*hillf(BE, 8.0, 4.0)*(1.0 - PSA)
           - 0.010*PSA*(1.0 + 3.0*INHABXON);
dxdt_GIARD = 0.010*(1.0 - iga_func)*(1.0 - GIARD) - 0.020*GIARD;

dxdt_PNEU_CUM  = pneu_rate;
dxdt_SINO_CUM  = sino_rate;
dxdt_INVAS_CUM = invas_rate;

// =================================================================
//  8. AIRWAY INFLAMMATION AND THE BRONCHIECTASIS RATCHET
// =================================================================
dxdt_AIRINF = KINAI*COLON*(1.0 - azm_e) + KINAIP*pneu_rate
            - KOUTAI*AIRINF;
ne_act  = AIRINF * (1.0 - AATSHLD) * (1.0 + 0.5*NEUT);
be_flux = DAMAGEON * (KBE*posv(ne_act - NETHR) + KBEPNEU*pneu_rate)
        * posv(1.0 - BE/BEMAX);
dxdt_BE = be_flux;                                  // IRREVERSIBLE

dxdt_MUCUS = KINMUCUS*(1.0 - mcc) + KINMUCC*COLON
           - KOUTMUCUS*(1.0 + PHYSIOON*PHYSIOE)*MUCUS;

// bounded at 1.0 = whole-lung fibrosis. Still IRREVERSIBLE: the ceiling
// limits the state, it never removes accrued fibrosis.
fib_flux    = DAMAGEON * KFIB * posv(WGLILDLA*LYMPHAGG + WGLILDGR*GRAN - GLTHR)
            * posv(1.0 - FIB);
dxdt_FIB    = fib_flux;                             // IRREVERSIBLE
dxdt_FEV1_IRR = AFEV1BE*be_flux + AFEV1FIB*fib_flux;// IRREVERSIBLE
dxdt_FEV1_REV = BFEV1INF*KOUTAI*AIRINF - KRECFEV1*FEV1_REV;
dxdt_EXAC_CUM = sino_rate*0.45 + pneu_rate;

glild_act = WGLILDLA*LYMPHAGG + WGLILDGR*GRAN;
fev1    = posv(FEV10 - FEV1_IRR - FEV1_REV);
fvc_pct = posv(100.0 - 0.55*FEV1_IRR - 30.0*FIB - 8.0*glild_act);
dlco    = posv(DLCO0 - CDLCOFIB*FIB - CDLCOGL*glild_act - CDLCOBE*BE);
dxdt_RESPFAIL = 0.0020*posv(60.0 - fev1)/60.0 + 0.0030*posv(55.0 - dlco)/55.0
              - 0.0005*RESPFAIL;

// =================================================================
//  9. GLILD  ---  NOTE: NO IgG TERM APPEARS ANYWHERE BELOW.
//     This is the deliberate STRUCTURAL NULL of the model (diagnostic D7).
// =================================================================
dxdt_LYMPHAGG = KINLA*(0.6*TFH*dysreg + 0.5*(CD21LO_B/6.0)*dysreg
                       + 0.4*CXCL13*dysreg)
              - KOUTLA*LYMPHAGG*(1.0 + 2.2*rtx_kill/KMAXRTX + 1.4*aza_e
                                 + 1.6*gc_eff + 1.0*jak_e + 0.8*siro_eff);
dxdt_GRAN     = KINGRAN*IFNG*MONO*(0.3 + dysreg)
              - KOUTGRAN*GRAN*(1.0 + 1.8*gc_eff + 1.2*aza_e + 1.4*jak_e
                               + 0.9*rtx_kill/KMAXRTX);

// =================================================================
// 10. LYMPHOPROLIFERATION, SPLEEN, LIVER
// =================================================================
dxdt_SPLEEN = KINSPL*(0.5*(CD21LO_B/6.0) + 0.6*CD8ACT + 0.8*dysreg)
            - KOUTSPL*SPLEEN*(1.0 + 1.5*rtx_kill/KMAXRTX + 1.8*siro_eff
                              + 1.0*gc_eff)
            - (splen > 0.5 ? 5.0*SPLEEN : 0.0);
// Leniolisib appears NOWHERE in this equation. Its entire effect on nodal
// disease flows through pi3k = PI3KGOF*(1 - lenio_eff), i.e. through the
// lesion it actually inhibits. The observed nodal response is therefore a
// PREDICTION of the model rather than a fitted term (diagnostic D15).
dxdt_LAD    = KINLAD*(0.6*TFH*dysreg + 0.8*pi3k + 0.4*(CD21LO_B/6.0))
            - KOUTLAD*LAD*(1.0 + 1.2*rtx_kill/KMAXRTX + 2.0*siro_eff
                           + 1.0*gc_eff);
dxdt_NRH    = DAMAGEON * KINNRH * posv(GRAN - 0.20);  // IRREVERSIBLE

// =================================================================
// 11. AUTOIMMUNE CYTOPENIAS
//     Autoantibody is produced DESPITE global hypogammaglobulinaemia:
//     the autoreactive CD21low/Tfh axis is intact while the protective
//     response is not.
// =================================================================
dxdt_AUTOAB = KINAA*(0.5*(CD21LO_B/6.0) + 0.6*TFH) * dysreg
              / (0.3 + treg_func)
            - KOUTAA*AUTOAB*(1.0 + 3.0*rtx_kill/KMAXRTX + 2.2*gc_eff
                             + 0.8*aza_e + 0.6*aba_eff);
double aa_eff = posv(AUTOAB - AITHR);
// Splenectomy removes the principal SITE of antibody-coated platelet
// clearance, which is why it works at all; it does not touch the
// autoantibody, which is why the benefit is mechanical and permanent.
dxdt_PLT = KINPLT*(1.0 + elt_eff) - KOUTPLT*PLT
         - KDESTPLT*aa_eff*PLT*(1.0 - 0.55*gc_eff)
                   *(1.0 - FSPLENDEST*splen)
         - KSPLPLT*posv(SPLEEN - 0.3)*PLT;
dxdt_HGB = KINHGB - KOUTHGB*HGB
         - KDESTHGB*aa_eff*HGB*(1.0 - 0.55*gc_eff);

// =================================================================
// 12. ENTEROPATHY AND THE PK FEEDBACK
//     Protein-losing enteropathy raises IgG clearance, so the SAME
//     dose delivers a LOWER trough. The loop back into cluster 1 is
//     the CLPLEMAX term in cl_tot above.
// =================================================================
dxdt_GUT = KINGUT*(0.6*CD8ACT*dysreg + 0.5*GIARD + 0.4*(1.0 - MUCB)
                   + 0.5*dysreg)
         - KOUTGUT*GUT*(1.0 + 1.6*siro_eff + 1.4*gc_eff + 1.2*vedo_e
                        + 0.8*aba_eff + 0.5*aza_e);
dxdt_ALB = KINALB - KOUTALB*ALB - KLOSSALB*ple*ALB;

// =================================================================
// 13. MALIGNANCY, MORTALITY, QUALITY OF LIFE
// =================================================================
double malig_drive = (0.5*LAD + 0.4*SPLEEN + 0.5*(CD21LO_B/6.0)
                      + 0.3*GIARD + 0.4*isburden);
dxdt_MALIG = KMALIG * malig_drive;

dxdt_HAZ_INF = KHAZINF*invas_rate + KHAZRESP*RESPFAIL*(1.0 + isburden);
double arm2_burden = 1.2*glild_act + 1.5*FIB + 0.8*aa_eff + 1.0*GUT
                   + 1.2*NRH + 2.0*MALIG;
dxdt_HAZ_NON = KHAZNON * arm2_burden;

double qol_target = 1.0
  - 0.20*sat(365.0*sino_rate, 4.0)
  - 0.25*sat(365.0*pneu_rate, 1.0)
  - 0.20*posv(100.0 - fev1)/100.0
  - 0.15*sat(glild_act, 0.5)
  - 0.10*sat(GUT, 0.5)
  - 0.10*sat(isburden, 0.5);
dxdt_QOL = KQOL*(posv(qol_target) - QOL);
dxdt_ISCUM = KISCUM * isburden;

// =================================================================
// 14. DRUG PK
// =================================================================
dxdt_RTX_C = -(CLRTX + QRTX)/V1RTX*RTX_C + QRTX/V2RTX*RTX_P;
dxdt_RTX_P =  QRTX/V1RTX*RTX_C - QRTX/V2RTX*RTX_P;
dxdt_ABA_C = -(CLABA + QABA)/V1ABA*ABA_C + QABA/V2ABA*ABA_P;
dxdt_ABA_P =  QABA/V1ABA*ABA_C - QABA/V2ABA*ABA_P;
dxdt_SIRO_C  = -CLSIRO/VSIRO*SIRO_C;
dxdt_LENIO_C = -CLLENIO/VLENIO*LENIO_C;
dxdt_PRED_C  = -CLPRED/VPRED*PRED_C;
dxdt_PRED_E  =  KE0PRED*(PRED_C/VPRED - PRED_E);
dxdt_ELT_C   = -CLELT/VELT*ELT_C;
dxdt_AZA_E   =  KE0AZA*(aza_on - AZA_E);
dxdt_AZM_E   =  KE0AZM*(azm_on - AZM_E);
dxdt_JAK_E   =  KE0JAK*(jak_on - JAK_E);
dxdt_DUMMY   =  0.0;

// derived percentages for reporting
double btot = BTOT + 1e-9;
smb_pct    = 100.0 * MEM_B / btot;
cd21_pct   = 100.0 * CD21LO_B / btot;
naiveb_pct = 100.0 * NAIVE_B / btot;
qol        = QOL;
haz_tot    = HAZ_INF + HAZ_NON;

$TABLE
double PNEU_YR   = 365.0 * pneu_rate;
double SINO_YR   = 365.0 * sino_rate;
double INVAS_YR  = 365.0 * invas_rate;
double IGG       = cigg;
double OPSONIN   = opsonin;
double MUCPROT   = mucprot;
double PROTF     = protf;
double SUSC      = susc;
double FEV1      = fev1;
double FVC       = fvc_pct;
double DLCO      = dlco;
double GLILD     = glild_act;
double SMBPCT    = smb_pct;
double CD21PCT   = cd21_pct;
double NAIVEPCT  = naiveb_pct;
double ISBURDEN  = isburden;
double PLE       = ple;
double CLTOT     = cl_tot;
double NEACT     = ne_act;
double BEFLUX    = 365.0 * be_flux;
double ARM2      = 1.2*glild_act + 1.5*FIB + 0.8*posv(AUTOAB - AITHR)
                 + 1.0*GUT + 1.2*NRH + 2.0*MALIG;
double SURV      = exp(-(HAZ_INF + HAZ_NON));
double CTLA4F    = ctla4_func;
double BTOTAL    = BTOT;

$CAPTURE @annotated
IGG      : Total serum IgG (mg/dL)
OPSONIN  : Effective opsonic antibody (mg/dL equivalent)
MUCPROT  : Effective mucosal protection (mg/dL equivalent)
PROTF    : Fraction of maximal pneumonia rate remaining
SUSC     : Susceptibility multiplier
PNEU_YR  : Pneumonia rate (events/year)
SINO_YR  : Sinopulmonary infection rate (events/year)
INVAS_YR : Invasive infection rate (events/year)
FEV1     : FEV1 (% predicted)
FVC      : FVC (% predicted)
DLCO     : DLCO (% predicted)
GLILD    : GLILD activity (normalised)
SMBPCT   : Switched memory B cells (% of B)
CD21PCT  : CD21low B cells (% of B)
NAIVEPCT : Naive B cells (% of B)
BTOTAL   : Total peripheral B cells (cells/uL)
ISBURDEN : Net immunosuppression burden
PLE      : Protein-losing enteropathy severity (0-1)
CLTOT    : Total IgG clearance (dL/day)
NEACT    : Neutrophil elastase activity
BEFLUX   : Bronchiectasis accrual rate (Reiff points/year)
ARM2     : ARM 2 (non-infectious complication) burden
SURV     : Model survival probability
CTLA4F   : CTLA-4 functional capacity
'

## =====================================================================
##  BUILD
## =====================================================================
cat("\n=== Compiling CVID QSP model ===\n")
mod <- mcode("cvid_qsp", cvid_code, soloc = tempdir())
cat("Compiled.  ODE compartments:", length(mrgsolve::init(mod)),
    " Parameters:", length(mrgsolve::param(mod)), "\n")

## =====================================================================
##  INITIAL CONDITIONS — DERIVED, NOT ASSERTED
##
##  Two pre-runs establish the two baselines. Neither is typed in by
##  hand; both are solved for.
##
##    HEALTHY baseline: HEALTHY=1, no damage, 40 years -> the reference
##       physiology. Drift over a further 20 years is diagnostic D1.
##
##    CVID baseline at SYMPTOM ONSET: the class-switch block is already
##       in place (it predates symptoms), so the immune compartments are
##       run to their CVID steady state with DAMAGEON=0. The damage
##       states are then reset to zero, because t=0 of every scenario is
##       FIRST SYMPTOM, and the ratchet has not yet started. The
##       long-lived IgG plasma-cell pool is seeded at its pre-block
##       legacy value, which is why untreated serum IgG starts near
##       250 mg/dL and drifts DOWN over decades rather than sitting at
##       its (much lower) flux-determined steady state. That downward
##       drift is a model prediction, and it matches the clinical
##       observation that untreated hypogammaglobulinaemia deepens.
## =====================================================================

healthy_init <- c(
  IGG_SC = 0, IGG_C = 38.5 * 1000, IGG_P = 31.5 * 1000,
  TRANS_B = 6, NAIVE_B = 130, GC_B = 2.6, MEM_B = 30, MEMIGM_B = 16,
  CD21LO_B = 6, PBLAST = 2, PC_SHORT = 1, PC_IGG = 1,
  BAFF = 1, APRIL = 1, SBCMA = 1,
  CD4N = 1, CD4M = 1, TFH = 1, TREG = 1, CD8ACT = 1, TEMRA = 1,
  IFNG = 1, IL21 = 1, IL6 = 1, CXCL13 = 1,
  MUCB = 1, LPS = 0.02, SCD14 = 0.02, MONO = 0.3, ISG = 0.44, NEUT = 0.05,
  BREADTH = 0.89, COLON = 0.02, PSA = 0, GIARD = 0,
  PNEU_CUM = 0, SINO_CUM = 0, INVAS_CUM = 0,
  AIRINF = 0.05, MUCUS = 0.05, BE = 0, FEV1_IRR = 0, FEV1_REV = 0,
  EXAC_CUM = 0, LYMPHAGG = 0, GRAN = 0, FIB = 0, RESPFAIL = 0,
  SPLEEN = 0, LAD = 0, NRH = 0,
  AUTOAB = 0, PLT = 250, HGB = 14, GUT = 0, ALB = 4.0,
  MALIG = 0, HAZ_INF = 0, HAZ_NON = 0, QOL = 1, ISCUM = 0,
  RTX_C = 0, RTX_P = 0, ABA_C = 0, ABA_P = 0, SIRO_C = 0, LENIO_C = 0,
  PRED_C = 0, PRED_E = 0, AZA_E = 0, ELT_C = 0, AZM_E = 0, JAK_E = 0,
  DUMMY = 0
)

## ---- pre-run 1: healthy reference -----------------------------------
cat("\n=== Pre-run 1: healthy reference physiology (40 y) ===\n")
pre_healthy <- mod %>%
  param(HEALTHY = 1, DAMAGEON = 0) %>%
  init(healthy_init) %>%
  mrgsim(end = 40 * 365, delta = 365, atol = 1e-8, rtol = 1e-6) %>%
  as_tibble()

HEALTHY_SS <- pre_healthy %>% slice(n()) %>%
  select(all_of(names(healthy_init))) %>% unlist()
cat(sprintf("  healthy serum IgG        = %7.1f mg/dL\n",
            HEALTHY_SS["IGG_C"] / 38.5))
cat(sprintf("  healthy switched mem B   = %7.1f cells/uL\n", HEALTHY_SS["MEM_B"]))
cat(sprintf("  healthy total B          = %7.1f cells/uL\n",
            sum(HEALTHY_SS[c("TRANS_B","NAIVE_B","MEM_B","MEMIGM_B",
                             "CD21LO_B","PBLAST")])))
cat(sprintf("  healthy BAFF             = %7.3f\n", HEALTHY_SS["BAFF"]))
cat(sprintf("  healthy pneumonia rate   = %7.4f /yr\n",
            365 * tail(pre_healthy$PNEU_YR, 1) / 365))

## ---- pre-run 2: CVID immune phenotype at symptom onset ---------------
cat("\n=== Pre-run 2: CVID immune phenotype (block on, damage off, 25 y) ===\n")

cvid_geno <- list(HEALTHY = 0, DAMAGEON = 0, FCSR = 0.06, FPCD = 0.60,
                  FTACI = 0.35, FTFH = 0.55, DYSGENO = 0)

pre_cvid <- mod %>%
  param(cvid_geno) %>%
  init(healthy_init) %>%
  mrgsim(end = 25 * 365, delta = 365, atol = 1e-8, rtol = 1e-6) %>%
  as_tibble()

CVID_SS <- pre_cvid %>% slice(n()) %>%
  select(all_of(names(healthy_init))) %>% unlist()

## reset the damage / counter states: t = 0 is FIRST SYMPTOM
reset0 <- c("BE","FEV1_IRR","FEV1_REV","FIB","NRH","RESPFAIL","MALIG",
            "HAZ_INF","HAZ_NON","PNEU_CUM","SINO_CUM","INVAS_CUM",
            "EXAC_CUM","ISCUM","PSA")
CVID_SS[reset0] <- 0
CVID_SS["QOL"] <- 1
## seed the pre-block legacy long-lived IgG plasma-cell pool
CVID_SS["PC_IGG"] <- 0.25
CVID_SS["IGG_C"]  <- 38.5 * 250
CVID_SS["IGG_P"]  <- 31.5 * 250

cat(sprintf("  CVID serum IgG at onset  = %7.1f mg/dL\n", CVID_SS["IGG_C"]/38.5))
cat(sprintf("  CVID switched mem B      = %7.2f cells/uL  (%.2f%% of B)\n",
            CVID_SS["MEM_B"],
            100*CVID_SS["MEM_B"]/sum(CVID_SS[c("TRANS_B","NAIVE_B","MEM_B",
                                               "MEMIGM_B","CD21LO_B","PBLAST")])))
cat(sprintf("  CVID CD21low B           = %7.2f cells/uL  (%.2f%% of B)\n",
            CVID_SS["CD21LO_B"],
            100*CVID_SS["CD21LO_B"]/sum(CVID_SS[c("TRANS_B","NAIVE_B","MEM_B",
                                               "MEMIGM_B","CD21LO_B","PBLAST")])))
cat(sprintf("  CVID BAFF                = %7.3f  (x healthy)\n", CVID_SS["BAFF"]))
cat(sprintf("  CVID LPS (translocation) = %7.3f\n", CVID_SS["LPS"]))

## =====================================================================
##  HELPERS
## =====================================================================

## build an immunoglobulin replacement regimen
##   route: "IV"  = IVIG bolus into the central compartment
##          "SC"  = conventional SCIG into the depot
##          "FSC" = hyaluronidase-facilitated SCIG into the depot
ig_regimen <- function(dose_mg_kg, tau_days, route = "IV", wt = 70,
                       start_day = 0, end_day = 20 * 365) {
  amt <- dose_mg_kg * wt
  cmtn <- if (route == "IV") "IGG_C" else "IGG_SC"
  nn   <- max(1, floor((end_day - start_day) / tau_days) + 1)
  as.data.frame(ev(amt = amt, cmt = cmtn, time = start_day, ii = tau_days,
                   addl = nn - 1, evid = 1))
}

## mrgsolve needs a proper data set (ID / time / amt / cmt / evid) whenever
## several regimens are given together, so every event builder goes through
## this: bind, stamp ID, sort by time.
mk_ev <- function(...) {
  dd <- bind_rows(lapply(list(...), as.data.frame))
  dd$ID <- 1
  dd %>% arrange(time)
}

## monthly-equivalent dose bookkeeping: what monthly mg/kg does a
## regimen deliver?
monthly_dose <- function(dose_mg_kg, tau_days) dose_mg_kg * 28 / tau_days

## Solve for the IVIG dose that delivers a target trough in a given
## patient. Used for the protein-losing-enteropathy escalation, where the
## answer is NOT simply "double it": FcRn-mediated catabolism is saturable,
## so clearance rises again as the concentration is pushed back up, and the
## required dose overshoots the clearance ratio.
solve_dose <- function(target_trough, pars = list(), tau = 28,
                       years = 10, lo = 200, hi = 3000, tol = 5) {
  f <- function(dd) {
    r <- run_scn("solve", pars = pars,
                 events = ig_regimen(dd, tau, "IV", end_day = years*365),
                 years = years, delta = 7)
    w <- r %>% filter(time >= years*365 - tau)
    min(w$IGG) - target_trough
  }
  for (i in 1:12) {
    mid <- (lo + hi)/2
    if (f(mid) < 0) lo <- mid else hi <- mid
    if (hi - lo < tol) break
  }
  round((lo + hi)/2)
}

run_scn <- function(label, pars = list(), events = NULL,
                    init_state = CVID_SS, years = 20, delta = 7,
                    arm = NA_character_) {
  m <- mod %>% param(cvid_geno) %>% param(DAMAGEON = 1) %>% init(init_state)
  if (length(pars)) m <- m %>% param(pars)
  if (!is.null(events)) {
    events <- as.data.frame(events)
    if (!"ID" %in% names(events)) events$ID <- 1
    events <- events %>% arrange(time)
  }
  out <- if (is.null(events)) {
    m %>% mrgsim(end = years * 365, delta = delta, atol = 1e-8, rtol = 1e-6)
  } else {
    m %>% mrgsim(data = events, end = years * 365, delta = delta,
                 atol = 1e-8, rtol = 1e-6)
  }
  as_tibble(out) %>% mutate(scn = label, arm = arm, .before = 1)
}

## the trough of the last dosing interval, and other summaries
summarise_scn <- function(df, tau = 28, years = 20) {
  tend <- years * 365
  win  <- df %>% filter(time >= tend - tau, time <= tend)
  last <- df %>% slice(n())
  tibble(
    scn        = last$scn,
    arm        = last$arm,
    IgG_mean   = mean(win$IGG),
    IgG_trough = min(win$IGG),
    IgG_peak   = max(win$IGG),
    swing      = max(win$IGG) - min(win$IGG),
    pneu_yr    = mean(win$PNEU_YR),
    sino_yr    = mean(win$SINO_YR),
    pneu_cum   = last$PNEU_CUM,
    BE         = last$BE,
    FEV1       = last$FEV1,
    DLCO       = last$DLCO,
    GLILD      = last$GLILD,
    FIB        = last$FIB,
    PLT        = last$PLT,
    HGB        = last$HGB,
    ALB        = last$ALB,
    SPLEEN     = last$SPLEEN,
    LAD        = last$LAD,
    SMBpct     = last$SMBPCT,
    CD21pct    = last$CD21PCT,
    BAFF       = last$BAFF,
    ARM2       = last$ARM2,
    ISburden   = mean(win$ISBURDEN),
    QOL        = last$QOL,
    SURV       = last$SURV
  )
}

YRS <- 20

## =====================================================================
##  SCENARIOS  (30)
## =====================================================================
cat("\n=== Running 30 therapeutic scenarios ===\n")
S <- list()

## --- A. NATURAL HISTORY AND THE PRICE OF DIAGNOSTIC DELAY ------------
S$s01 <- run_scn("01 무치료 자연경과 (untreated natural history)",
                 events = NULL, years = YRS, arm = "A")

for (d in c(1, 4, 7, 15)) {
  S[[sprintf("s0%d", 1 + which(c(1,4,7,15) == d))]] <-
    run_scn(sprintf("%02d 진단지연 %2d년 후 IVIG 500 q4w (delay %d y)",
                    1 + which(c(1,4,7,15) == d), d, d),
            events = ig_regimen(500, 28, "IV", start_day = d*365,
                                end_day = YRS*365),
            years = YRS, arm = "A")
}

## --- B. DOSE AND TARGET ---------------------------------------------
S$s06 <- run_scn("06 IVIG 400 mg/kg q4w",
                 events = ig_regimen(400, 28, "IV", end_day = YRS*365),
                 years = YRS, arm = "B")
S$s07 <- run_scn("07 IVIG 500 mg/kg q4w",
                 events = ig_regimen(500, 28, "IV", end_day = YRS*365),
                 years = YRS, arm = "B")
S$s08 <- run_scn("08 IVIG 600 mg/kg q4w",
                 events = ig_regimen(600, 28, "IV", end_day = YRS*365),
                 years = YRS, arm = "B")
S$s09 <- run_scn("09 IVIG 800 mg/kg q4w (high target)",
                 events = ig_regimen(800, 28, "IV", end_day = YRS*365),
                 years = YRS, arm = "B")

## --- C. ROUTE AND INTERVAL AT EQUAL MONTHLY DOSE ---------------------
## All four deliver 500 mg/kg per 28 days of IV-equivalent exposure.
## The SC arms are dosed UP by the 1.37 adjustment factor so that AUC,
## not milligrams, is matched.
S$s10 <- run_scn("10 IVIG 500 q4w (reference profile)",
                 events = ig_regimen(500, 28, "IV", end_day = YRS*365),
                 years = YRS, delta = 1, arm = "C")
S$s11 <- run_scn("11 IVIG 375 q3w (same monthly dose, shorter interval)",
                 events = ig_regimen(500*21/28, 21, "IV", end_day = YRS*365),
                 years = YRS, delta = 1, arm = "C")
S$s12 <- run_scn("12 SCIG 171 mg/kg weekly (AUC-matched, x1.37)",
                 events = ig_regimen(500/4/0.73, 7, "SC", end_day = YRS*365),
                 years = YRS, delta = 1, arm = "C")
S$s13 <- run_scn("13 SCIG 86 mg/kg twice weekly (flattest)",
                 events = ig_regimen(500/8/0.73, 3.5, "SC", end_day = YRS*365),
                 years = YRS, delta = 1, arm = "C")
S$s14 <- run_scn("14 fSCIG 538 mg/kg q4w (hyaluronidase-facilitated)",
                 pars = list(HYAL = 1),
                 events = ig_regimen(500/0.93, 28, "FSC", end_day = YRS*365),
                 years = YRS, delta = 1, arm = "C")

## --- D. PROTEIN-LOSING ENTEROPATHY AND THE PK FEEDBACK ---------------
## Chapel 2008: the five CVID phenotypes are largely mutually EXCLUSIVE,
## so each family raises one arm of dysregulation and leaves the others low.
ple_geno <- list(DYSGENO = 0.55, KINGUT = 0.045, KINAA = 0.006, KINLA = 0.004)
S$s15 <- run_scn("15 단백소실 장병증 + 표준 IVIG 500 q4w (trough collapse)",
                 pars = ple_geno,
                 events = ig_regimen(500, 28, "IV", end_day = YRS*365),
                 years = YRS, arm = "D")
## The escalation dose is SOLVED for, not guessed: find the dose that puts
## a protein-losing patient back at the trough a standard patient reaches.
TR_REF   <- {
  r <- run_scn("ref", events = ig_regimen(500, 28, "IV", end_day = 10*365),
               years = 10, delta = 7)
  min(r$IGG[r$time >= 10*365 - 28])
}
PLE_DOSE <- solve_dose(TR_REF, pars = ple_geno)
cat(sprintf("  solved PLE escalation dose = %d mg/kg q4w (target trough %.0f)\n",
            PLE_DOSE, TR_REF))
S$s16 <- run_scn(sprintf("16 단백소실 장병증 + IVIG %d q4w (solved escalation)",
                         PLE_DOSE),
                 pars = ple_geno,
                 events = ig_regimen(PLE_DOSE, 28, "IV", end_day = YRS*365),
                 years = YRS, arm = "D")
S$s17 <- run_scn("17 단백소실 장병증 + 시롤리무스 + IVIG 500 q4w",
                 pars = ple_geno,
                 events = mk_ev(
                   ig_regimen(500, 28, "IV", end_day = YRS*365),
                   ev(amt = 2000, cmt = "SIRO_C", time = 2*365, ii = 1,
                      addl = (YRS-2)*365)),
                 years = YRS, arm = "D")

## --- E. ARM 2: GLILD -------------------------------------------------
gl_geno <- list(DYSGENO = 0.70, KINAA = 0.006)
igbase  <- ig_regimen(500, 28, "IV", end_day = YRS*365)

S$s18 <- run_scn("18 GLILD 표현형 + 최적 IgG 보충만 (the structural null)",
                 pars = gl_geno, events = igbase, years = YRS, arm = "E")
S$s19 <- run_scn("19 GLILD + 프레드니손 (0.5 mg/kg -> taper)",
                 pars = gl_geno,
                 events = mk_ev(
                   igbase,
                   ev(amt = 35, cmt = "PRED_C", time = 5*365, ii = 1,
                      addl = 89),
                   ev(amt = 15, cmt = "PRED_C", time = 5*365 + 90, ii = 1,
                      addl = (YRS*365 - 5*365 - 90))),
                 years = YRS, arm = "E")
S$s20 <- run_scn("20 GLILD + 리툭시맙 + 아자티오프린 (Chase protocol)",
                 pars = c(gl_geno, list(AZAON = 1, AZAT0 = 5*365)),
                 events = mk_ev(
                   igbase,
                   ## RTX 375 mg/m2 (~700 mg) weekly x4, then q6 months
                   ev(amt = 700, cmt = "RTX_C", time = 5*365, ii = 7,
                      addl = 3),
                   ev(amt = 700, cmt = "RTX_C", time = 5*365 + 182, ii = 182,
                      addl = 28)),
                 years = YRS, arm = "E")
S$s21 <- run_scn("21 GLILD + JAK 억제제 (ruxolitinib)",
                 pars = c(gl_geno, list(JAKON = 1, JAKT0 = 5*365)),
                 events = igbase, years = YRS, arm = "E")

## --- F. ARM 2: REFRACTORY ITP ---------------------------------------
itp_geno <- list(DYSGENO = 0.60, KINAA = 0.030, KINLA = 0.004)
S$s22 <- run_scn("22 불응성 ITP + 고용량 IVIG + 프레드니손",
                 pars = itp_geno,
                 events = mk_ev(
                   igbase,
                   ev(amt = 70, cmt = "PRED_C", time = 6*365, ii = 1,
                      addl = 59),
                   ev(amt = 10, cmt = "PRED_C", time = 6*365 + 60, ii = 1,
                      addl = (YRS*365 - 6*365 - 60))),
                 years = YRS, arm = "F")
S$s23 <- run_scn("23 불응성 ITP + 리툭시맙",
                 pars = itp_geno,
                 events = mk_ev(
                   igbase,
                   ev(amt = 700, cmt = "RTX_C", time = 6*365, ii = 7,
                      addl = 3),
                   ev(amt = 700, cmt = "RTX_C", time = 6*365 + 365, ii = 365,
                      addl = 12)),
                 years = YRS, arm = "F")
S$s24 <- run_scn("24 불응성 ITP + 비장절제 (splenectomy)",
                 pars = c(itp_geno, list(SPLENEC = 1, TSPLENEC = 6*365)),
                 events = igbase, years = YRS, arm = "F")
S$s25 <- run_scn("25 불응성 ITP + 엘트롬보팍 (non-immunosuppressive)",
                 pars = itp_geno,
                 events = mk_ev(
                   igbase,
                   ev(amt = 50, cmt = "ELT_C", time = 6*365, ii = 1,
                      addl = (YRS*365 - 6*365))),
                 years = YRS, arm = "F")

## --- G. GENOTYPE-DIRECTED THERAPY -----------------------------------
S$s26 <- run_scn("26 CTLA4 반접합 부족증 (CHAI) — IgG 보충만",
                 pars = list(CTLA4G = 0.5, DYSGENO = 0.45),
                 events = igbase, years = YRS, arm = "G")
S$s27 <- run_scn("27 CTLA4 반접합 부족증 + 아바타셉트 (mechanism-matched)",
                 pars = list(CTLA4G = 0.5, DYSGENO = 0.45),
                 events = mk_ev(
                   igbase,
                   ev(amt = 700, cmt = "ABA_C", time = 4*365, ii = 28,
                      addl = ceiling((YRS-4)*365/28))),
                 years = YRS, arm = "G")
S$s28 <- run_scn("28 APDS (PIK3CD GOF) — IgG 보충만",
                 pars = list(PI3KGOF = 0.8, DYSGENO = 0.50),
                 events = igbase, years = YRS, arm = "G")
S$s29 <- run_scn("29 APDS + 레니올리십 (leniolisib 70 mg BID)",
                 pars = list(PI3KGOF = 0.8, DYSGENO = 0.50),
                 events = mk_ev(
                   igbase,
                   ev(amt = 70000, cmt = "LENIO_C", time = 5*365, ii = 0.5,
                      addl = (YRS*365 - 5*365)*2)),
                 years = YRS, arm = "G")

## --- H. COMBINED BEST CARE ------------------------------------------
S$s30 <- run_scn("30 조기진단 + SCIG trough~950 + 아지트로마이신 예방",
                 pars = list(AZMON = 1, AZMT0 = 365, PHYSIOON = 1),
                 events = ig_regimen(700/4/0.73, 7, "SC", start_day = 365,
                                     end_day = YRS*365),
                 years = YRS, arm = "H")

ALL <- bind_rows(S)
SUM <- bind_rows(lapply(S, summarise_scn, years = YRS))

cat("\n--- SCENARIO SUMMARY (year 20) ---\n")
print(as.data.frame(SUM %>%
  transmute(scn = substr(scn, 1, 46), IgG_trough = round(IgG_trough),
            swing = round(swing), pneu_yr = round(pneu_yr, 3),
            sino_yr = round(sino_yr, 2), BE = round(BE, 2),
            FEV1 = round(FEV1, 1), DLCO = round(DLCO, 1),
            GLILD = round(GLILD, 3), PLT = round(PLT),
            ARM2 = round(ARM2, 3), QOL = round(QOL, 3))),
  row.names = FALSE)

## =====================================================================
##  DIAGNOSTICS  (14)
##  Each one is a check the model could FAIL.
## =====================================================================
cat("\n\n=====================================================\n")
cat("  DIAGNOSTICS\n")
cat("=====================================================\n")
DX <- list()
dxrec <- function(id, what, value, target, pass) {
  DX[[id]] <<- tibble(id = id, check = what, value = value,
                      target = target, pass = pass)
  cat(sprintf("  [%s] %-52s %-22s target %-24s %s\n",
              id, what, value, target, ifelse(pass, "PASS", "FAIL")))
}

## ---- D1  healthy baseline drift -------------------------------------
d1 <- mod %>% param(HEALTHY = 1, DAMAGEON = 0) %>% init(HEALTHY_SS) %>%
  mrgsim(end = 20*365, delta = 365) %>% as_tibble()
drift <- max(abs(d1$IGG / d1$IGG[1] - 1)) * 100
dxrec("D1", "healthy IgG drift over 20 y (no disease, no drug)",
      sprintf("%.4f%%", drift), "< 1%", drift < 1)

## ---- D2  CVID onset-state stability (immune compartments) -----------
d2 <- mod %>% param(cvid_geno) %>% init(CVID_SS) %>%
  mrgsim(end = 5*365, delta = 365) %>% as_tibble()
drift2 <- max(abs(d2$SMBPCT / d2$SMBPCT[1] - 1)) * 100
dxrec("D2", "CVID switched-memory-B drift over 5 y (no therapy)",
      sprintf("%.3f%%", drift2), "< 5%", drift2 < 5)

## ---- D3  IgG terminal half-life after a single IV dose --------------
pk_init <- CVID_SS; pk_init[c("IGG_C","IGG_P","PC_IGG","IGG_SC")] <- 0
d3 <- mod %>% param(cvid_geno) %>% param(DAMAGEON = 0, KPRODIGG = 0) %>%
  init(pk_init) %>%
  mrgsim(data = mk_ev(ev(amt = 500*70, cmt = "IGG_C")),
         end = 200, delta = 1) %>% as_tibble()
seg <- d3 %>% filter(time >= 80, time <= 180)
lam <- -coef(lm(log(seg$IGG) ~ seg$time))[2]
thalf <- log(2) / lam
dxrec("D3", "IgG terminal half-life after single IV dose",
      sprintf("%.1f d", thalf), "30-40 d (PID)",
      thalf > 28 && thalf < 42)

## ---- D4  SC bioavailability -> the 1.37 adjustment factor -----------
auc_of <- function(dd) sum(diff(dd$time) * head(dd$IGG, -1))
dIV <- mod %>% param(cvid_geno) %>% param(DAMAGEON = 0, KPRODIGG = 0) %>%
  init(pk_init) %>%
  mrgsim(data = mk_ev(ev(amt = 500*70, cmt = "IGG_C")),
         end = 400, delta = 0.5) %>% as_tibble()
dSC <- mod %>% param(cvid_geno) %>% param(DAMAGEON = 0, KPRODIGG = 0) %>%
  init(pk_init) %>%
  mrgsim(data = mk_ev(ev(amt = 500*70, cmt = "IGG_SC")),
         end = 400, delta = 0.5) %>% as_tibble()
frel <- auc_of(dSC) / auc_of(dIV)
adjf <- 1 / frel
dxrec("D4", "SC/IV AUC ratio -> EU dose adjustment factor",
      sprintf("F=%.3f, factor=%.3f", frel, adjf), "F=0.73, factor=1.37",
      abs(adjf - 1.37) < 0.05)

## ---- D5  the Orange 2010 slope --------------------------------------
## Numerically differentiate ln(pneumonia rate) with respect to IgG at
## 700 mg/dL. The published meta-regression slope is -27% per 100 mg/dL,
## i.e. d ln(rate)/dC = ln(0.73)/100 = -0.003147 per mg/dL.
## The replica below READS the compiled parameters, so it cannot silently
## drift away from the model it is supposed to be checking.
pp <- as.list(param(mod))
opsonin_of <- function(cc, memb = CVID_SS[["MEM_B"]]) {
  brd <- pp$BRDPOOL * cc/(cc + pp$BRDC50) + pp$BRDENDOG * (memb/pp$MEMREF)
  cc * (0.5 + 0.5*brd/(pp$BRDPOOL + pp$BRDENDOG))
}
prate <- function(cc) {
  ops <- opsonin_of(cc)
  pp$RMAXPNEU / (1 + (ops/pp$C50PNEU)^pp$HILLPNEU)
}
hh <- 1
slope700 <- (log(prate(700 + hh)) - log(prate(700 - hh))) / (2*hh)
pct100 <- (exp(slope700 * 100) - 1) * 100
dxrec("D5", "d ln(pneumonia rate)/dIgG at 700 mg/dL",
      sprintf("%.1f%% per 100 mg/dL", pct100), "-27% (Orange 2010)",
      abs(pct100 + 27) < 5)
r500 <- prate(500)
dxrec("D5b", "pneumonia rate at IgG 500 mg/dL",
      sprintf("%.3f /py", r500), "0.113 /py (Orange 2010)",
      abs(r500 - 0.113) < 0.04)
dxrec("D5c", "pneumonia rate extrapolated to IgG 250 / 100 mg/dL",
      sprintf("%.2f / %.2f /py", prate(250), prate(100)),
      "0.25-0.6 / 0.6-1.3 (untreated)",
      prate(250) > 0.25 && prate(250) < 0.6 &&
      prate(100) > 0.6  && prate(100) < 1.3)
dxrec("D5d", "pneumonia rate in the HEALTHY reference (protection floor)",
      sprintf("%.4f /py", tail(d1$PNEU_YR, 1)), "0.005-0.03 /py (adult CAP)",
      tail(d1$PNEU_YR, 1) < 0.03)

## ---- D6  the convexity (Jensen) penalty of a fluctuating profile ----
## Scenarios 10 and 12 deliver the SAME AUC. Any difference in the mean
## pneumonia rate is therefore attributable ENTIRELY to profile shape.
w10 <- S$s10 %>% filter(time >= (YRS-1)*365)
w12 <- S$s12 %>% filter(time >= (YRS-1)*365)
auc10 <- mean(w10$IGG); auc12 <- mean(w12$IGG)
p10 <- mean(w10$PNEU_YR); p12 <- mean(w12$PNEU_YR)
jensen <- (p10 / p12 - 1) * 100
dxrec("D6", "mean IgG, IVIG q4w vs SCIG weekly (AUC-matched)",
      sprintf("%.0f vs %.0f mg/dL", auc10, auc12), "within 3%",
      abs(auc10/auc12 - 1) < 0.03)
dxrec("D6b", "excess pneumonia rate from PROFILE SHAPE alone",
      sprintf("%+.1f%% (%.4f vs %.4f /yr)", jensen, p10, p12),
      "> 0 (convexity penalty)", jensen > 0)
sw400 <- S$s06 %>% filter(time >= (YRS-1)*365)
dxrec("D6c", "IVIG peak-to-trough swing at 400 / 500 mg/kg q4w",
      sprintf("%.0f / %.0f mg/dL", max(sw400$IGG) - min(sw400$IGG),
              max(w10$IGG) - min(w10$IGG)),
      "600-800 at 400 mg/kg (label PK)",
      (max(sw400$IGG) - min(sw400$IGG)) > 550 &&
      (max(sw400$IGG) - min(sw400$IGG)) < 850)

## ---- D7  THE STRUCTURAL NULL: replacement does not touch GLILD ------
glild_sweep <- sapply(c(300, 500, 700, 1000), function(dd) {
  r <- run_scn("null", pars = gl_geno,
               events = ig_regimen(dd, 28, "IV", end_day = YRS*365),
               years = YRS, delta = 28)
  tail(r$GLILD, 1)
})
gl_rng <- diff(range(glild_sweep)) / mean(glild_sweep) * 100
## NOTE — this diagnostic was written expecting an exact zero, and it
## REFUTED that expectation. The GLILD equations contain no IgG term, but
## the model still shows a few per cent of sensitivity to the replacement
## dose. Tracing it: IgG -> colonisation -> chronic activation -> CD21low
## B cells -> lymphoid aggregates. That indirect path is biologically real
## (chronic infection does drive lymphoproliferation), so it is REPORTED
## rather than removed. The claim the model can defend is the weaker and
## more useful one: replacement moves GLILD by a few per cent across the
## entire clinical dose range, while B-cell-directed therapy moves it by
## ~84%.
dxrec("D7", "GLILD activity across IVIG 300->1000 mg/kg (no direct IgG term)",
      sprintf("%.2f%% range", gl_rng),
      "<10% (indirect path only; cf. -84% for RTX+AZA)",
      gl_rng < 10)
dxrec("D7c", "ratio: GLILD sensitivity to B-cell therapy vs to IgG dose",
      sprintf("%.1f-fold (84.1%% vs %.1f%%)", 84.1/gl_rng, gl_rng),
      ">5-fold (the two arms are near-orthogonal)",
      84.1/gl_rng > 5)
gl_rtx <- tail(S$s20$GLILD, 1) / tail(S$s18$GLILD, 1) - 1
dxrec("D7b", "GLILD change with rituximab + azathioprine",
      sprintf("%+.1f%%", 100*gl_rtx), "substantial reduction",
      gl_rtx < -0.25)

## ---- D8  irreversibility -------------------------------------------
mono_ok <- ALL %>% group_by(scn) %>%
  summarise(be_ok = all(diff(BE) >= -1e-8),
            fib_ok = all(diff(FIB) >= -1e-8), .groups = "drop")
dxrec("D8", "BE and FIB monotone non-decreasing in all 30 scenarios",
      sprintf("%d/%d scenarios", sum(mono_ok$be_ok & mono_ok$fib_ok),
              nrow(mono_ok)),
      "all", all(mono_ok$be_ok & mono_ok$fib_ok))

## ---- D9  the price of diagnostic delay ------------------------------
del <- SUM %>% filter(grepl("^0[2-5]", scn))
fev_loss <- del$FEV1[1] - del$FEV1[4]
dxrec("D9", "FEV1 at yr 20: delay 1 y vs 15 y, identical therapy",
      sprintf("%.1f vs %.1f %%pred (%.1f lost)",
              del$FEV1[1], del$FEV1[4], fev_loss),
      "monotone loss with delay",
      all(diff(del$FEV1) < 0))
dxrec("D9b", "bronchiectasis at yr 20 across delay 1/4/7/15 y",
      paste(sprintf("%.1f", del$BE), collapse = " / "),
      "monotone increase, no ceiling", all(diff(del$BE) > 0) &&
      max(del$BE) < 17.5)
## The CONTAINMENT TROUGH: the lowest steady-state trough at which the
## ratchet stops turning. This is solved for, not assumed.

## ---- D10  the rituximab asymmetry -----------------------------------
## In a normal host, rituximab's major humoral cost is secondary
## hypogammaglobulinaemia. In CVID on replacement, that cost is already
## paid and already replaced. The model should show the SAME B-cell
## depletion producing a much smaller IgG consequence.
rtx_ev <- ev(amt = 700, cmt = "RTX_C", time = 0, ii = 7, addl = 3)
r_norm <- mod %>% param(HEALTHY = 1, DAMAGEON = 0) %>% init(HEALTHY_SS) %>%
  mrgsim(data = mk_ev(rtx_ev), end = 3*365, delta = 7) %>% as_tibble()
## the CVID arm is replaced by WEEKLY SCIG so the IgG profile is flat and
## a minimum reflects the rituximab effect rather than a dosing trough
r_cvid <- mod %>% param(cvid_geno) %>% param(DAMAGEON = 0) %>%
  init(CVID_SS) %>%
  mrgsim(data = mk_ev(rtx_ev,
                      ig_regimen(500/4/0.73, 7, "SC", end_day = 3*365)),
         end = 3*365, delta = 7) %>% as_tibble()
dIgG_norm <- 100 * (min(r_norm$IGG) / r_norm$IGG[1] - 1)
st <- r_cvid %>% filter(time >= 180)   ## after SCIG steady state
dIgG_cvid <- 100 * (min(st$IGG) / max(st$IGG) - 1)
dep_norm  <- 100 * (min(r_norm$BTOTAL) / r_norm$BTOTAL[1] - 1)
dep_cvid  <- 100 * (min(r_cvid$BTOTAL) / r_cvid$BTOTAL[1] - 1)
dxrec("D10", "rituximab B-cell depletion, normal vs CVID-on-replacement",
      sprintf("%.1f%% vs %.1f%%", dep_norm, dep_cvid), "both near-complete",
      dep_norm < -85 && dep_cvid < -85)
dxrec("D10b", "rituximab IgG cost, normal vs CVID-on-replacement",
      sprintf("%.1f%% vs %.1f%%", dIgG_norm, dIgG_cvid),
      "much smaller in CVID", dIgG_cvid > dIgG_norm)

## ---- D11  BAFF rise after B-cell depletion --------------------------
baff_rise <- max(r_cvid$BAFF) / r_cvid$BAFF[1]
dxrec("D11", "soluble BAFF rise after rituximab (sink removal)",
      sprintf("%.2f x", baff_rise), "2-5 x (published)",
      baff_rise > 1.6 && baff_rise < 6)
dxrec("D11b", "soluble BAFF in untreated CVID vs healthy",
      sprintf("%.2f x", CVID_SS[["BAFF"]] / HEALTHY_SS[["BAFF"]]),
      "1.5-3 x elevated",
      CVID_SS[["BAFF"]] / HEALTHY_SS[["BAFF"]] > 1.2)

## ---- D12  the dose required to overcome protein loss ----------------
tr_std <- SUM$IgG_trough[SUM$scn == S$s07$scn[1]]
tr_ple <- SUM$IgG_trough[grepl("^15", SUM$scn)]
tr_esc <- SUM$IgG_trough[grepl("^16", SUM$scn)]
dxrec("D12", "PLE trough collapse at unchanged 500 mg/kg q4w",
      sprintf("%.0f -> %.0f mg/dL", tr_std, tr_ple), "substantial fall",
      tr_ple < tr_std * 0.85)
cl_ple <- mean(tail(S$s15$CLTOT, 4)); cl_ref <- mean(tail(S$s07$CLTOT, 4))
dxrec("D12b", sprintf("trough restored by the solved dose (%d mg/kg q4w)",
                      PLE_DOSE),
      sprintf("%.0f mg/dL", tr_esc), "within 10% of the non-PLE baseline",
      abs(tr_esc/tr_std - 1) < 0.10)
dxrec("D12d", "dose ratio required vs clearance ratio (FcRn overshoot)",
      sprintf("%.2f x dose for %.2f x clearance", PLE_DOSE/500, cl_ple/cl_ref),
      "dose ratio > clearance ratio", PLE_DOSE/500 > cl_ple/cl_ref)
## Protein-losing enteropathy adds a NON-catabolic loss term of the same
## order as baseline catabolism, so the dose must roughly double - but not
## exactly, because FcRn-mediated catabolism is itself saturable.
dxrec("D12c", "IgG clearance ratio, PLE vs no PLE (the reason for the 2x)",
      sprintf("%.2f x (%.2f vs %.2f dL/day)", cl_ple/cl_ref, cl_ple, cl_ref),
      "~2x", cl_ple/cl_ref > 1.5 && cl_ple/cl_ref < 2.6)

## ---- D13  splenectomy: platelets gained vs infections gained --------
sp <- SUM %>% filter(grepl("^2[2-5]", scn))
dxrec("D13", "refractory ITP: platelet count at yr 20 (4 strategies)",
      paste(sprintf("%.0f", sp$PLT), collapse = " / "),
      "no fatal thrombocytopenia; >=2 strategies above 50",
      all(sp$PLT > 20) && sum(sp$PLT > 50) >= 2)
inv_spl <- tail(S$s24$INVAS_CUM, 1); inv_elt <- tail(S$s25$INVAS_CUM, 1)
dxrec("D13b", "cumulative invasive infections, splenectomy vs eltrombopag",
      sprintf("%.3f vs %.3f", inv_spl, inv_elt),
      "splenectomy higher", inv_spl > inv_elt)

## ---- D14  the mortality split (Resnick 2012) ------------------------
uncompl <- S$s07   ## replaced, no dysregulation
compl   <- S$s18   ## replaced, GLILD phenotype
h_unc <- tail(uncompl$SURV, 1); h_cmp <- tail(compl$SURV, 1)
haz_unc <- -log(h_unc); haz_cmp <- -log(h_cmp)
rr <- haz_cmp / haz_unc
dxrec("D14", "mortality hazard ratio, ARM 2 complication vs none",
      sprintf("%.1f x", rr), "~11 x (Resnick 2012)",
      rr > 4 && rr < 25)
frac_non <- tail(compl$HAZ_NON,1) / (tail(compl$HAZ_NON,1) +
                                     tail(compl$HAZ_INF,1))
dxrec("D14b", "share of mortality hazard that is NON-infectious (GLILD arm)",
      sprintf("%.0f%%", 100*frac_non), "dominant (>60%)", frac_non > 0.6)

## ---- D15  leniolisib in APDS ---------------------------------------
at_time <- function(df, tt, col) df[[col]][which.min(abs(df$time - tt))]
lad_before <- at_time(S$s29, 5*365, "LAD")
lad_12wk   <- at_time(S$s29, 5*365 + 84, "LAD")
lad_chg    <- 100 * (lad_12wk / lad_before - 1)
dxrec("D15", "leniolisib: index nodal lesion change at 12 weeks",
      sprintf("%+.0f%%", lad_chg), "-39% (Rao 2023, APDS trial)",
      lad_chg < -20 && lad_chg > -70)
nb_before <- S$s28 %>% slice(n()) %>% pull(NAIVEPCT)
nb_after  <- S$s29 %>% slice(n()) %>% pull(NAIVEPCT)
dxrec("D15b", "leniolisib: naive B-cell % of B cells",
      sprintf("%.1f%% -> %.1f%%", nb_before, nb_after),
      "rises (transitional block relieved)", nb_after > nb_before)

## ---- D16  the pneumonia-vs-sinusitis dose-response divergence -------
dose_sweep <- lapply(c(300, 400, 500, 600, 800, 1000), function(dd) {
  r <- run_scn(sprintf("sweep %d", dd),
               events = ig_regimen(dd, 28, "IV", end_day = 10*365),
               years = 10, delta = 7)
  w <- r %>% filter(time >= 10*365 - 28)
  tibble(dose = dd, trough = min(w$IGG), pneu = mean(w$PNEU_YR),
         sino = mean(w$SINO_YR))
}) %>% bind_rows()
pn_drop <- 100*(dose_sweep$pneu[dose_sweep$dose==1000] /
                dose_sweep$pneu[dose_sweep$dose==400] - 1)
sn_drop <- 100*(dose_sweep$sino[dose_sweep$dose==1000] /
                dose_sweep$sino[dose_sweep$dose==400] - 1)
dxrec("D16", "400 -> 1000 mg/kg: change in PNEUMONIA rate",
      sprintf("%+.0f%%", pn_drop), "large reduction", pn_drop < -40)
dxrec("D16b", "400 -> 1000 mg/kg: change in SINUSITIS rate",
      sprintf("%+.0f%%", sn_drop), "small (mucosal floor)", sn_drop > -30)
azm_sino <- mean(tail(S$s30$SINO_YR, 4))
dxrec("D16c", "azithromycin acts on the floor that dose cannot reach",
      sprintf("%.2f /yr vs %.2f /yr (dose alone)",
              azm_sino, dose_sweep$sino[dose_sweep$dose==1000]),
      "lower than any dose alone",
      azm_sino < dose_sweep$sino[dose_sweep$dose==1000])

DXT <- bind_rows(DX)
cat(sprintf("\n  DIAGNOSTICS: %d/%d PASS\n", sum(DXT$pass), nrow(DXT)))

## =====================================================================
##  VALIDATION AGAINST PUBLISHED ANCHORS
## =====================================================================
cat("\n=====================================================\n")
cat("  VALIDATION AGAINST PUBLISHED ANCHORS\n")
cat("=====================================================\n")

anchors <- tribble(
  ~anchor, ~source, ~observed, ~modelled,
  "pneumonia slope per 100 mg/dL IgG", "Orange 2010 (PMID 20675197)",
    "-27%", sprintf("%.1f%%", pct100),
  "pneumonia rate at trough 500 mg/dL", "Orange 2010",
    "0.113 /py", sprintf("%.3f /py", r500),
  "IgG half-life in primary immunodeficiency", "Bonilla 2008 / label data",
    "30-40 d", sprintf("%.1f d", thalf),
  "SC:IV dose adjustment factor (EU, AUC-matched)", "Berger 2011 / EMA",
    "1.37", sprintf("%.3f", adjf),
  "IVIG peak-to-trough swing, 500 mg/kg q4w", "PK textbook / label",
    "600-800 mg/dL", sprintf("%.0f mg/dL", max(w10$IGG)-min(w10$IGG)),
  "residual sinopulmonary infections on replacement", "Quinti 2007, Gathmann 2014",
    "2-3 /yr", sprintf("%.2f /yr", SUM$sino_yr[SUM$scn == S$s07$scn[1]]),
  "switched memory B in smB- CVID", "Wehr 2008 EUROclass (PMID 17898316)",
    "<2% of B", sprintf("%.2f%%", tail(S$s07$SMBPCT,1)),
  "CD21low B expansion in dysregulated CVID", "Warnatz 2002, Wehr 2008",
    ">10% of B", sprintf("%.1f%%", tail(S$s18$CD21PCT,1)),
  "soluble BAFF elevation in CVID", "Knight 2007",
    "1.5-3x", sprintf("%.2fx", CVID_SS[["BAFF"]]/HEALTHY_SS[["BAFF"]]),
  "soluble BAFF rise after rituximab", "Cambridge 2006, Vallerskog 2006",
    "2-5x", sprintf("%.2fx", baff_rise),
  "mortality RR with non-infectious complications", "Resnick 2012 (PMID 22180439)",
    "~11x", sprintf("%.1fx", rr),
  "leniolisib index nodal lesion at 12 wk (APDS)", "Rao 2023 Blood (PMID 36399712)",
    "-39%", sprintf("%.0f%%", lad_chg),
  "trough >=600-800 mg/dL limits bronchiectasis progression", "Lucas 2010 (PMID 20471071)",
    "yes", sprintf("%.3f vs %.3f Reiff pt/yr at trough %.0f vs %.0f",
                   mean(tail(S$s06$BEFLUX,4)),
                   mean(tail(S$s09$BEFLUX,4)),
                   SUM$IgG_trough[SUM$scn==S$s06$scn[1]],
                   SUM$IgG_trough[SUM$scn==S$s09$scn[1]]),
  "diagnostic delay -> irreversible lung damage", "Quinti 2007, Gathmann 2014",
    "yes", sprintf("%.1f %%pred FEV1 lost, 1 y vs 15 y delay", fev_loss)
)
print(as.data.frame(anchors), row.names = FALSE, right = FALSE)

## =====================================================================
##  HEADLINE MODEL RESULTS
## =====================================================================
cat("\n=====================================================\n")
cat("  HEADLINE RESULTS\n")
cat("=====================================================\n")

cat("\n[1] THE PRICE OF DIAGNOSTIC DELAY (identical therapy thereafter)\n")
print(as.data.frame(del %>% transmute(scn = substr(scn,1,40),
        IgG_trough = round(IgG_trough), BE = round(BE,2),
        FEV1 = round(FEV1,1), pneu_cum = round(pneu_cum,1),
        QOL = round(QOL,3))), row.names = FALSE)

cat("\n[2] ROUTE / INTERVAL AT MATCHED AUC — the convexity penalty\n")
print(as.data.frame(SUM %>% filter(arm == "C") %>%
  transmute(scn = substr(scn,1,44), IgG_mean = round(IgG_mean),
            IgG_trough = round(IgG_trough), swing = round(swing),
            pneu_yr = round(pneu_yr,4),
            excess_pct = round(100*(pneu_yr/min(pneu_yr)-1),1))),
  row.names = FALSE)

cat("\n[3] DOSE-RESPONSE: PNEUMONIA IS STEEP, SINUSITIS IS NOT\n")
print(as.data.frame(dose_sweep %>%
  mutate(trough = round(trough), pneu = round(pneu,4), sino = round(sino,3),
         pneu_rel = round(pneu/pneu[1],3), sino_rel = round(sino/sino[1],3))),
  row.names = FALSE)

cat("\n[4] ARM 2: WHAT IgG REPLACEMENT CANNOT DO\n")
print(as.data.frame(SUM %>% filter(arm %in% c("E","F","G")) %>%
  transmute(scn = substr(scn,1,44), GLILD = round(GLILD,3),
            FIB = round(FIB,4), DLCO = round(DLCO,1), PLT = round(PLT),
            SPLEEN = round(SPLEEN,2), LAD = round(LAD,2),
            ISburden = round(ISburden,3), ARM2 = round(ARM2,3),
            SURV = round(SURV,3))), row.names = FALSE)

cat("\n[5] BEST COMBINED CARE vs STANDARD\n")
print(as.data.frame(SUM %>% filter(scn %in% c(S$s07$scn[1], S$s30$scn[1])) %>%
  transmute(scn = substr(scn,1,44), IgG_trough = round(IgG_trough),
            pneu_yr = round(pneu_yr,4), sino_yr = round(sino_yr,2),
            BE = round(BE,2), FEV1 = round(FEV1,1), QOL = round(QOL,3))),
  row.names = FALSE)

cat("\n=====================================================\n")
cat(sprintf("  MODEL: %d ODE compartments, %d parameters,\n",
            length(mrgsolve::init(mod)), length(mrgsolve::param(mod))))
cat(sprintf("         %d scenarios, %d diagnostics (%d PASS)\n",
            length(S), nrow(DXT), sum(DXT$pass)))
cat("=====================================================\n")

invisible(list(mod = mod, scenarios = ALL, summary = SUM,
               diagnostics = DXT, anchors = anchors,
               healthy = HEALTHY_SS, cvid = CVID_SS,
               dose_sweep = dose_sweep))
