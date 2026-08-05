## =============================================================================
##  gsd1a_mrgsolve_model.R
##  Glycogen Storage Disease Type Ia (von Gierke) — QSP model for mrgsolve
##  (with a switchable GSD Ib / SLC37A4 branch)
## =============================================================================
##
##  WHAT THIS MODEL IS ABOUT
##  ------------------------
##  Glucose-6-phosphatase (G6PC1) is the single terminal step shared by BOTH
##  routes of endogenous glucose production.  Glycogenolysis and gluconeogenesis
##  both converge on glucose-6-phosphate, and only G6Pase releases free glucose
##  from it.  Losing it does not merely reduce hepatic glucose output — it turns
##  the liver from the body's glucose SOURCE into an obligate glucose SINK,
##  because glucokinase keeps running forward with nothing to reverse it.
##
##  Everything else in GSD Ia is the overflow.  The model is therefore written
##  as ONE mass balance on hepatic G6P with four exits, and the four classical
##  syndromes come out of the branch fractions rather than being four separately
##  parameterised submodels:
##
##      G6P -> glycogen   (G6P allosterically activates glycogen synthase)
##      G6P -> lactate    (glycolysis)
##      G6P -> lipid      (acetyl-CoA -> de novo lipogenesis)
##      G6P -> PRPP       (pentose phosphate pathway -> purines -> urate)
##
##  47 ODEs.  Independently re-implemented in Python/scipy as
##  `gsd1a_reference_model.py`; every number quoted in README.md is produced by
##  integrating one of the two and checking it against the other.
##
##  UNITS
##    time            h
##    glucose         mmol/L plasma, mmol/h flux
##    G6P, glycogen   umol/g liver
##    lactate         mmol/L
##    urate           mg/dL
##    triglyceride    mg/dL plasma, mg/g liver
##    conversion      1 mg/kg/min glucose = 0.33333 mmol/kg/h
##
##  DOSING COMPARTMENTS
##    AST     uncooked cornstarch / Glycosade, as mmol glucose equivalents
##            (1 g cornstarch ~ 0.97 g glucose equivalents ~ 5.39 mmol)
##    ALLOg   allopurinol, mg, oral
##    ACEIc   ACE inhibitor, ug/L, into the central compartment
##    EMPAc   empagliflozin, nmol/L, into the central compartment
##    AAVvg   AAV8-G6PC transduced load, dimensionless (0.22 ~ 1e13 GC/kg)
##    MRNA    LNP-mRNA derived G6Pase protein, dimensionless
##
##  CALIBRATION TARGETS (all reproduced by the shipped parameter set)
##    residual EGP in GSD Ia          0.9 mg/kg/min at age 1, 0.3 in adults
##                                    (~15 % of demand; see the limitation note
##                                    in README.md -- isotopic studies suggest
##                                    25-35 %, and the model sits below them)
##    hepatic glycogen                90-140 mg/g        (normal ~50)
##    liver volume                    1.8-2.2 x normal
##    blood lactate, treated          2-5 mmol/L; crisis 10-15
##    serum urate, treated            6-8 mg/dL; poor control 11-15
##    plasma triglyceride             400-1000 mg/dL treated
##    beta-hydroxybutyrate            INAPPROPRIATELY LOW for the FFA level
##    fasting tolerance               3.2 h at age 1, 5.9 h in adults
##    gain from DOUBLING the starch   1.60-1.64 h at every age (= ln2/kdis)
## =============================================================================

library(mrgsolve)

gsd1a_code <- '
$PROB
# GSD Ia (von Gierke) QSP model — hepatic G6P branch point with four overflow
# exits, dietary/adjunctive/gene/mRNA therapy, and a switchable GSD Ib arm.

$PARAM @annotated
// ---------------- patient ---------------------------------------------------
AGE      :  1.0  : age (years)
BW       : 10.0  : body weight (kg)
GENO     :  1    : genotype 0=control 1=GSD Ia 2=GSD Ib
RESID    :  0.0  : germline residual G6Pase activity (fraction of normal)

// ---------------- whole-body glucose demand ---------------------------------
GURA     :  2.10 : GUR asymptote in adults (mg/kg/min)
GURB     :  4.60 : GUR paediatric excess (mg/kg/min)
GURTAU   :  7.00 : GUR age decay constant (years)
FCNSA    :  0.42 : cerebral share of demand, adult asymptote
FCNSB    :  0.36 : cerebral share, paediatric excess
FCNSTAU  :  6.00 : cerebral share decay constant (years)
VGKG     :  0.20 : glucose distribution volume (L/kg)
VLACKG   :  0.50 : lactate distribution volume (L/kg)
VPLKG    :  0.045: plasma volume (L/kg)
VUAKG    :  0.30 : urate distribution volume (L/kg)
FLWA     :  0.024: liver mass fraction, adult
FLWB     :  0.014: liver mass fraction, paediatric excess
FLWTAU   :  4.00 : liver mass fraction decay constant (years)

// ---------------- cerebral fuel ---------------------------------------------
VLACCNS  :  0.90 : max lactate share of cerebral demand
KLACCNS  :  4.00 : lactate half-saturation for brain uptake (mmol/L)
KB       :  0.70 : apparent GLUT1 Km at the BBB (mmol/L)
FAITHR   :  0.830: fuel adequacy index below which neuroglycopenia begins

// ---------------- gut / starch ----------------------------------------------
KDIS     :  0.45 : starch amylolysis rate (1/h) UCCS 0.45, waxy-maize 0.28
FABS     :  0.75 : fraction of starch reaching portal vein as glucose
KA       :  3.00 : absorption of luminal free glucose (1/h)
KCOL     :  0.20 : colonic disposal of resistant starch (1/h)
RDRIP    :  0.0  : mean enteral glucose delivery (mmol/h)
RDRIPAMP :  0.0  : day-cycle amplitude, 0=drip, ~0.9=three meals
RDRIPPER : 24.0  : day-cycle period (h)

// ---------------- hepatic fluxes --------------------------------------------
FVGK     :  0.28 : glucokinase Vmax as multiple of GUR
KGK      :  8.00 : glucokinase Km (mmol/L)
NGK      :  1.70 : glucokinase Hill coefficient
GKINS    :  0.60 : GK activation by insulin (GKRP release)
FVG6P    :  6.20 : G6Pase Vmax as multiple of GUR (large: normal liver idles)
KMG6P    :  1.50 : G6Pase Km for G6P (umol/g)
FEGPEX   :  0.10 : extrahepatic (renal + intestinal) G6Pase share
GLYC0    :310.0  : normal hepatic glycogen (umol glucosyl/g)
KGP      :  0.048: phosphorylase rate constant (1/h)
FVGS     :  1.30 : glycogen synthase Vmax as multiple of GUR
KGS      :  1.10 : G6P allosteric activation constant for GS (umol/g)
FDEB     :  0.100: alpha-1,6 branch glucose released free by debrancher
KLYS     :  0.005: lysosomal (GAA) glycogen -> free glucose (1/h)
GLYCMAX  :950.0  : physical glycogen ceiling (umol/g)
FVPFK    :  0.95 : glycolysis Vmax as multiple of GUR
KPFK     :  3.00 : glycolysis Km for G6P (umol/g)
KACIDPFK :  9.00 : bicarbonate at which PFK-1 is half-inhibited (mmol/L)
FVPPP    :  0.10 : pentose phosphate Vmax as multiple of GUR
KPPP     :  3.00 : PPP Km for G6P (umol/g)
FVGNG    :  0.55 : gluconeogenesis Vmax as multiple of GUR
KGNG     :  2.00 : gluconeogenic substrate Km (mmol/L)
SUB0     :  1.20 : alanine + glycerol equivalent substrate supply (mmol/L)

// ---------------- lactate ---------------------------------------------------
FVLACOX  :  0.75 : peripheral lactate oxidation Vmax as multiple of GUR
KLACOX   :  3.50 : lactate oxidation Km (mmol/L)
KLACREN  :  0.55 : renal lactate clearance above threshold (L/h per 10 kg)
LACTHR   :  5.00 : renal lactate threshold (mmol/L)
FLACPROD :  0.50 : peripheral lactate production as multiple of GUR
FLACOUT  :  0.80 : glycolytic pyruvate leaving as lactate (vs PDH)
FGNGLAC  :  0.55 : gluconeogenic carbon coming from lactate
LPRATIO  : 14.0  : lactate:pyruvate ratio (NORMAL in GSD I — not a redox defect)

// ---------------- acid-base -------------------------------------------------
HCO3N    : 24.0  : normal bicarbonate (mmol/L)
KBUF     :  1.20 : bicarbonate equilibration rate (1/h)
CLACHCO3 :  1.25 : bicarbonate lost per mmol/L lactate

// ---------------- lipid -----------------------------------------------------
KDNL     :  0.55 : de novo lipogenesis coupling to PDH flux
KMC      :  1.60 : malonyl-CoA turnover (1/h)
KMCB     :  1.00 : CPT1 inhibition constant for beta-oxidation
KMCKET   :  0.60 : CPT1 inhibition constant for ketogenesis
VVLDL    :  0.625: ApoB-limited VLDL-TG export Vmax (mg/g/h)
KVLDLTG  : 12.0  : hepatic TG at half-maximal export (mg/g)
VLPL     :110.0  : saturable plasma TG clearance (mg/dL/h)
KLPL     :500.0  : LPL Km (mg/dL)
KLINTG   :  0.035: non-saturable remnant TG clearance (1/h)
FLPLDEF  :  0.52 : LPL activity in GSD Ia relative to normal
KBOX     :  0.030: hepatic beta-oxidation (1/h)
KFFATG   :  0.22 : hepatic FFA re-esterification to TG
KLIP     :  0.55 : lipolysis rate scaler
KFFAU    :  1.90 : FFA clearance (1/h)
KKET     :  3.00 : ketogenesis rate constant
KKBU     :  1.30 : ketone utilisation (1/h)

// ---------------- hormones --------------------------------------------------
IMAX     :420.0  : maximal insulin secretion (pmol/L)
KINSG    :  7.00 : glucose EC50 for insulin (mmol/L)
HINS     :  2.40 : insulin secretion Hill
KINS     :  4.20 : insulin elimination (1/h)
GCGMAX   : 42.0  : maximal glucagon (pmol/L)
KGCG     :  4.30 : glucose IC50 for glucagon (mmol/L)
HGCG     :  3.00 : glucagon Hill
KGCGEL   :  5.00 : glucagon elimination (1/h)
EPIMAX   :  3.20 : maximal epinephrine (nmol/L)
KEPI     :  3.20 : glucose IC50 for epinephrine (mmol/L)
HEPI     :  6.00 : epinephrine Hill
KEPIEL   :  8.00 : epinephrine elimination (1/h)
CORTMAX  :620.0  : maximal cortisol (nmol/L)
KCORT    :  3.10 : glucose IC50 for cortisol (mmol/L)
KCORTEL  :  0.45 : cortisol elimination (1/h)
GHMAX    :  9.00 : maximal GH tone (ug/L)
KGH      :  3.40 : glucose IC50 for GH (mmol/L)
KGHEL    :  0.90 : GH elimination (1/h)
IGF1N    : 25.0  : normal IGF-1 (nmol/L)
KIGF     :  0.03 : IGF-1 turnover (1/h)

// ---------------- purine / urate --------------------------------------------
FUAPROD  :  2.60 : baseline urate production (mg/h per 10 kg)
KPRPP    :  2.20 : PRPP formation coupling
KPRPPUSE :  2.20 : PRPP consumption (1/h)
NPRPP    :  0.45 : compressive exponent, PRPP pool vs PPP flux
KATPDEG  :  0.15 : ATP degradation from Pi trapping
G6PTHATP :  1.20 : G6P above which Pi trapping starts (umol/g)
FUACL    :  0.115: urate clearance (L/h per 10 kg)
KILACUA  :  7.00 : lactate Ki for urate reabsorption (mmol/L)
KIKBUA   :  3.00 : ketone Ki for urate reabsorption (mmol/L)
EMAXALLO :  0.62 : maximal xanthine oxidase inhibition
IC50OXY  :  5.50 : oxypurinol IC50 (mg/L)
KAALLO   :  2.20 : allopurinol absorption (1/h)
KALLO    :  0.55 : allopurinol elimination (1/h)
FOXY     :  0.76 : fraction converted to oxypurinol
KOXY     :  0.031: oxypurinol elimination (1/h, t1/2 ~22 h)
FVDALLO  :  1.60 : allopurinol Vd (L per 10 kg)

// ---------------- kidney ----------------------------------------------------
KHF      :  0.0016  : hyperfiltration adaptation (1/h)
KDAM     :  8.5e-7  : nephron damage accrual
KREP     :  4.0e-6  : nephron repair
WLAC     :  0.20    : lactate weight in the damage driver
WHF      :  0.85    : hyperfiltration weight
WUA      :  0.045   : urate weight
WTG      :  0.00035 : triglyceride weight
KUACR    :  0.004   : UACR equilibration (1/h)
EMAXACEI :  0.55    : maximal ACE inhibitor effect
IC50ACEI : 12.0     : ACE inhibitor IC50 (ug/L)
KACEI    :  0.058   : ACE inhibitor elimination (1/h)
KCIT     :  0.05    : urinary citrate equilibration (1/h)
KNCA     :  2.0e-5  : nephrocalcinosis accrual

// ---------------- long-term liver -------------------------------------------
KHCAINI  :  3.1e-7  : adenoma initiation
KHCAGRO  :  1.4e-5  : adenoma growth
HCAMAX   :  0.35    : adenoma burden ceiling (fraction of liver volume)
NMCI     :  2.60    : supralinearity of the metabolic control index
KHCC     :  2.4e-6  : malignant transformation hazard

// ---------------- growth / bone ---------------------------------------------
KGRO     :  6.0e-5  : height SDS relaxation rate (1/h)
HTSDSTAR :  -0.40   : attainable height SDS on perfect control
HTSDSMCI :  2.60    : height SDS lost at maximal MCI
AGEMAT   : 18.0     : skeletal maturity (years)
KBMD     :  4.0e-6  : BMD Z relaxation rate (1/h)
BMDZTAR  :  0.20    : attainable BMD Z on perfect control
BMDZACID :  2.60    : BMD Z lost per unit acidosis
BMDZMCI  :  1.60    : BMD Z lost at maximal MCI

// ---------------- gene / mRNA therapy ---------------------------------------
EPSAAV   :  1.00    : fractional G6Pase activity per unit transduced load
KEXPR    :  0.020   : transgene expression onset (1/h)
KSIL     :  6.0e-6  : slow episome loss beyond growth dilution (1/h)
KMRNA    :  0.0072  : LNP-mRNA protein decay (1/h, t1/2 ~4 d)
EPSMRNA  :  1.00    : fractional activity per unit mRNA protein
KABEL    :  2.0e-5  : anti-capsid antibody decay (1/h)
LWADULT  : 1680.0   : adult liver mass (g)
TAUGROW  : 48213.0  : liver growth time constant (h) = 5.5 y

// ---------------- GSD Ib ----------------------------------------------------
AG15IN   :  0.85    : dietary 1,5-anhydroglucitol input (ug/mL/h)
KAG15    :  0.055   : 1,5-AG renal handling (1/h)
KAGUP    :  0.030   : neutrophil 1,5-AG uptake
KAG6POUT :  0.200   : G6PT-dependent 1,5-AG6P efflux (absent in Ib)
KAG6PDEG :  0.020   : G6PT-independent 1,5-AG6P loss (1/h)
ANCPROD  :  3.90    : neutrophil production (10^9/L/h)
KAG6P    :  1.00    : 1,5-AG6P at half-maximal marrow suppression
HAG6P    :  3.20    : marrow suppression Hill (switch-like)
KANC     :  1.00    : neutrophil turnover (1/h)
EMPAEMAX :  3.40    : maximal fold rise in 1,5-AG clearance
EMPAEC50 : 60.0     : empagliflozin EC50 (nmol/L)
KEMPA    :  0.055   : empagliflozin elimination (1/h)
GCSF     :  0.0     : G-CSF on/off

// ---------------- switches --------------------------------------------------
FIBRATE  :  0.0  : fibrate on/off
STATIN   :  0.0  : statin on/off
CITRATE  :  0.0  : alkali citrate on/off
ILLNESS  :  0.0  : intercurrent catabolic stress 0-1
CHRONIC  :  1.0  : accelerate slow states for long-horizon runs

$CMT @annotated
AST      : slow starch in gut lumen (mmol glucose-eq)
AGL      : luminal free glucose (mmol)
ACOL     : starch escaping to colon (mmol glucose-eq)
Gp       : plasma glucose (mmol/L)
G6P      : hepatic glucose-6-phosphate (umol/g liver)
Glyc     : hepatic glycogen (umol glucosyl/g liver)
LV       : liver volume (mL)
Lac      : blood lactate (mmol/L)
Pyr      : blood pyruvate (mmol/L)
HCO3     : plasma bicarbonate (mmol/L)
MalCoA   : hepatic malonyl-CoA (relative to normal)
TGliv    : hepatic triglyceride (mg/g)
TGpl     : plasma triglyceride (mg/dL)
FFA      : plasma free fatty acid (mmol/L)
KBOH     : 3-hydroxybutyrate (mmol/L)
Ins      : plasma insulin (pmol/L)
Gcg      : plasma glucagon (pmol/L)
Epi      : plasma epinephrine (nmol/L)
Cort     : plasma cortisol (nmol/L)
GHc      : growth hormone tone (ug/L)
IGF1     : IGF-1 (nmol/L)
PRPP     : hepatic PRPP pool (relative to normal)
UA       : serum urate (mg/dL)
ALLOg    : allopurinol, gut (mg)
ALLOc    : allopurinol, plasma (mg/L)
OXY      : oxypurinol, plasma (mg/L)
GFRrel   : GFR relative to age-normal
Rdam     : accumulated nephron damage
UACR     : urine albumin:creatinine (mg/g)
CitU     : urinary citrate (relative to normal)
NCa      : nephrocalcinosis burden
HCA      : hepatocellular adenoma burden (relative volume)
HCChaz   : cumulative HCC hazard
HtSDS    : height standard-deviation score
BMDz     : lumbar BMD Z-score
TimeHypo : cumulative hours with glucose < 3.9 mmol/L
TimeNGC  : cumulative hours with fuel adequacy below threshold
AUCdef   : cumulative cerebral fuel deficit (mmol)
AAVvg    : transduced hepatocyte vector load
MRNAp    : LNP-mRNA derived G6Pase protein
G6Pact   : EXPRESSED fractional G6Pase activity
AntiAAV  : anti-capsid neutralising antibody
AG15     : plasma 1,5-anhydroglucitol (ug/mL)
AG6Pn    : neutrophil 1,5-AG6P (relative)
ANC      : absolute neutrophil count (10^9/L)
EMPAc    : empagliflozin (nmol/L)
ACEIc    : ACE inhibitor (ug/L)

$GLOBAL
#define HILL(x, K, n) (pow(fmax((x), 0.0), (n)) / (pow((K), (n)) + pow(fmax((x), 0.0), (n))))
#define POS(x) (fmax((x), 0.0))

$MAIN
// ---- allometry / age scaling (evaluated once per subject) -------------------
double GURMGKG = GURA + GURB * exp(-AGE / GURTAU);
double FCNS    = FCNSA + FCNSB * exp(-AGE / FCNSTAU);
double FLW     = FLWA + FLWB * exp(-AGE / FLWTAU);
double LW0     = FLW * BW * 1000.0;                 // g
double GUR     = GURMGKG * BW / 3.0;                // mmol/h
double VG      = VGKG * BW;
double VLAC    = VLACKG * BW;
double VPL     = VPLKG * BW;
double VUA     = VUAKG * BW;
double DCNS    = FCNS * GUR;
double UPER0   = 0.30 * (1.0 - FCNS) * GUR;
double SI      = 0.70 * (1.0 - FCNS) * GUR / 60.0;
double JPPPREF = FVPPP * GUR * (0.25 / (KPPP + 0.25));

// GENO: 0 = control, 1 = GSD Ia (catalytic subunit), 2 = GSD Ib (transporter)
double FG6PT   = (GENO == 2) ? 0.0 : 1.0;
double RESID_E = (GENO == 0) ? 1.0 : ((GENO == 2) ? fmax(RESID, 0.03) : RESID);
double FLPL    = (GENO == 0) ? 1.0 : FLPLDEF;

// ---- initial conditions -----------------------------------------------------
double CTRL = (GENO == 0) ? 1.0 : 0.0;
Gp_0     = CTRL ? 4.80 : 3.60;
G6P_0    = CTRL ? 0.22 : 3.60;
Glyc_0   = GLYC0 * (CTRL ? 1.0 : 2.1);
LV_0     = LW0 / 1.05 * (CTRL ? 1.0 : 2.0);
Lac_0    = CTRL ? 1.00 : 5.00;
Pyr_0    = Lac_0 / LPRATIO;
HCO3_0   = CTRL ? 24.0 : 19.0;
MalCoA_0 = CTRL ? 1.00 : 2.20;
TGliv_0  = CTRL ? 22.0 : 80.0;
TGpl_0   = CTRL ? 70.0 : 700.0;
FFA_0    = 0.45;
KBOH_0   = 0.15;
Ins_0    = 60.0;
Gcg_0    = 12.0;
Epi_0    = 0.40;
Cort_0   = 250.0;
GHc_0    = 2.00;
IGF1_0   = IGF1N * (CTRL ? 1.0 : 0.6);
PRPP_0   = CTRL ? 1.0 : 2.0;
UA_0     = CTRL ? 3.4 : 7.5;
GFRrel_0 = 1.0;
CitU_0   = CTRL ? 1.0 : 0.7;
HtSDS_0  = CTRL ? 0.0 : -1.6;
BMDz_0   = CTRL ? 0.0 : -1.1;
AG15_0   = AG15IN / KAG15;
AG6Pn_0  = (GENO == 2) ? 1.60 : 0.30;
ANC_0    = (GENO == 2) ? 0.50 : 3.90;

$ODE
double gp   = POS(Gp);
double g6p  = POS(G6P);
double gly  = POS(Glyc);
double lac  = POS(Lac);
double hco3 = POS(HCO3);
double LWg  = fmax(POS(LV) * 1.05, 1.0);
double LR   = LWg / LW0;                       // liver-size scaling of fluxes

// ================= hormones =================================================
double Iss  = IMAX * HILL(gp, KINSG, HINS);
double Gss  = GCGMAX * (1.0 - HILL(gp, KGCG, HGCG)) + 4.0;
double Ess  = EPIMAX * (1.0 - HILL(gp, KEPI, HEPI)) + 0.15 + 0.9 * ILLNESS;
double Css  = CORTMAX * (1.0 - HILL(gp, KCORT, 4.0)) + 180.0 + 220.0 * ILLNESS;
double GHss = GHMAX * (1.0 - HILL(gp, KGH, 4.0)) + 1.2;
dxdt_Ins  = KINS   * (Iss  - POS(Ins));
dxdt_Gcg  = KGCGEL * (Gss  - POS(Gcg));
dxdt_Epi  = KEPIEL * (Ess  - POS(Epi));
dxdt_Cort = KCORTEL* (Css  - POS(Cort));
dxdt_GHc  = KGHEL  * (GHss - POS(GHc));

double INSREL = POS(Ins) / 120.0;
double GCGREL = POS(Gcg) / 12.0;
double EPIREL = POS(Epi) / 0.5;

// ================= gut delivery =============================================
double Jdis = KDIS * POS(AST);
double Ragut = KA * POS(AGL);
double Rd = RDRIP * (1.0 + RDRIPAMP * sin(2.0 * M_PI * SOLVERTIME / RDRIPPER));
dxdt_AST  = -Jdis;
dxdt_AGL  = Jdis * FABS + fmax(Rd, 0.0) - Ragut;
dxdt_ACOL = Jdis * (1.0 - FABS) - KCOL * POS(ACOL);

// ================= hepatic branch point =====================================
// (a) glucokinase: glucose IN, one-way.  In GSD Ia nothing ever reverses it,
//     so the liver is a net glucose CONSUMER at every plasma concentration.
double Jgk = FVGK * GUR * (1.0 + GKINS * fmin(INSREL, 2.5))
             * HILL(gp, KGK, NGK) * LR;

// (b) glycogen phosphorylase.  Hill constants sit ABOVE resting hormone levels
//     so a glucagon challenge has somewhere to go.
double phos = (0.10 + 2.2 * HILL(GCGREL, 2.5, 2.0) + 0.90 * HILL(EPIREL, 2.5, 2.0))
              / (1.0 + 0.9 * HILL(gp, 6.0, 3.0)) / (1.0 + 0.6 * fmin(INSREL, 3.0));
double Jgp  = KGP * gly * phos * LWg / 1000.0;
double Jdeb = FDEB * Jgp;                       // alpha-1,6 -> FREE glucose
double Jlys = KLYS * gly * LWg / 1000.0;        // lysosomal -> FREE glucose

// (c) gluconeogenesis: still runs, still dumps into the SAME blocked node
double gngh = 0.45 + 0.75 * HILL(GCGREL, 1.0, 1.4) + 0.45 * (POS(Cort) / 400.0);
double Jgng = FVGNG * GUR * gngh * HILL(lac + 0.6 * POS(FFA) + SUB0, KGNG, 1.0);

// (d) G6Pase — the lesion.  a(t) = germline residual + expressed transgene.
//     In Ib the catalytic subunit is intact but substrate never reaches it.
double a_enz = fmin(RESID_E + POS(G6Pact), 1.0);
double a_eff = (GENO == 2) ? a_enz * FG6PT : a_enz;
double Jg6p  = FVG6P * GUR * a_eff * HILL(g6p, KMG6P, 1.0) * LR;
double Jextra = FEGPEX * FVG6P * GUR * a_eff * HILL(g6p, KMG6P, 1.0);

// (e) glycogen synthase: covalent (insulin) arm TIMES allosteric (G6P) arm.
//     The allosteric arm is why this is a STORAGE disease: G6P switches GS on
//     even when the covalent arm says stop.
double gscov = (0.15 + 1.6 * HILL(POS(Ins), 110.0, 1.5))
               / (1.0 + 0.8 * HILL(GCGREL, 1.2, 2.0));
double gsall = 0.22 + 0.78 * HILL(g6p, KGS, 2.0);
double Jgs = FVGS * GUR * gscov * gsall * (1.0 - HILL(gly, GLYCMAX, 6.0)) * LR;

// (f) the two overflow drains.  PFK-1 is strongly acidosis-inhibited: that is
//     the only negative feedback that terminates a lactic crisis.
double f26 = 1.0 + 0.55 * fmin(INSREL, 2.0) - 0.20 * HILL(GCGREL, 1.0, 1.5);
double acidbrake = HILL(hco3, KACIDPFK, 3.0) / HILL(24.0, KACIDPFK, 3.0);
double Jgly = FVPFK * GUR * HILL(g6p, KPFK, 1.0) * fmax(f26, 0.25) * acidbrake * LR;
double Jppp = FVPPP * GUR * HILL(g6p, KPPP, 1.0) * LR;

dxdt_G6P  = (Jgk + (Jgp - Jdeb) + Jgng - Jg6p - Jgs - Jgly - Jppp) * 1000.0 / LWg;
dxdt_Glyc = (Jgs - Jgp - Jlys) * 1000.0 / LWg;

// ================= cerebral fuel and systemic glucose =======================
double supLac  = VLACCNS * DCNS * HILL(lac, KLACCNS, 1.0);
double needGlc = fmax(DCNS - supLac, 0.0);
double gotGlc  = needGlc * HILL(gp, KB, 1.0);
double FAI     = (DCNS > 0.0) ? (supLac + gotGlc) / DCNS : 1.0;

double Uper = (UPER0 + SI * POS(Ins)) * HILL(gp, 2.5, 1.0) * (1.0 - 0.25 * ILLNESS);
double Gthr = 10.0 - 3.4 * HILL(POS(EMPAc), EMPAEC50, 1.0);
double Uren = 0.9 * BW / 10.0 * fmax(gp - Gthr, 0.0) * POS(GFRrel);
double EGP  = Jg6p + Jextra + Jdeb + Jlys;

dxdt_Gp = (Ragut + EGP - gotGlc - Uper - Jgk - Uren) / VG;

// ================= lactate / pyruvate / acid-base ===========================
double Jlacox = FVLACOX * GUR * HILL(lac, KLACOX, 1.0) * (1.0 + 0.25 * EPIREL / 3.0);
double Jlacren = KLACREN * BW / 10.0 * fmax(lac - LACTHR, 0.0) * POS(GFRrel);
double LacPer = FLACPROD * GUR * (1.0 + 0.5 * HILL(EPIREL, 2.0, 1.0))
                + 0.35 * FLACPROD * GUR * ILLNESS;
double JlacCNS = supLac * 2.0;
dxdt_Lac = (2.0 * Jgly * FLACOUT + LacPer - Jlacox - Jlacren
            - JlacCNS - 2.0 * Jgng * FGNGLAC) / VLAC;
dxdt_Pyr = (lac / LPRATIO - POS(Pyr)) * 6.0;
double HCO3ss = HCO3N - CLACHCO3 * fmax(lac - 1.0, 0.0) + 4.0 * CITRATE;
dxdt_HCO3 = KBUF * (fmax(HCO3ss, 6.0) - hco3);

// ================= lipid ====================================================
double chrebp = HILL(g6p, 2.0, 1.8);
// carbon accounting: pyruvate splits between LDH (lactate) and PDH (lipid).
// Counting the same pyruvate into both would inflate lactate AND triglyceride.
double Jdnl = KDNL * 2.0 * Jgly * (1.0 - FLACOUT) * (0.35 + 0.9 * chrebp)
              * (1.0 + 0.5 * fmin(INSREL, 2.0));
dxdt_MalCoA = KMC * ((0.25 + 2.1 * chrebp + 0.5 * fmin(INSREL, 2.0)) - POS(MalCoA));
double lipo = KLIP * (1.0 + 1.4 * HILL(EPIREL, 2.0, 1.0) + 0.5 * HILL(GCGREL, 1.5, 1.0))
              / (1.0 + 1.8 * fmin(INSREL, 3.0));
dxdt_FFA = lipo - KFFAU * POS(FFA);
// malonyl-CoA stays high because ChREBP is driven by G6P, not by insulin.
// Hence HYPOketotic hypoglycaemia — the bedside discriminator from GSD 0/III/VI.
dxdt_KBOH = KKET * POS(FFA) / (1.0 + POS(MalCoA) / KMCKET) - KKBU * POS(KBOH);

double Jvldl = VVLDL * HILL(POS(TGliv), KVLDLTG, 1.0);
double Jbox  = KBOX * POS(TGliv) / (1.0 + POS(MalCoA) / KMCB);
dxdt_TGliv = Jdnl * 180.0 * 0.42 / LWg + KFFATG * POS(FFA) * 12.0 - Jvldl - Jbox;
double lpl = FLPL * (1.0 + 0.55 * fmin(INSREL, 2.0)) * (1.0 + 1.35 * FIBRATE);
double Jlpl = (VLPL * HILL(POS(TGpl), KLPL, 1.0) + KLINTG * POS(TGpl)) * lpl;
dxdt_TGpl = Jvldl * LWg / (VPL * 10.0) * (1.0 - 0.22 * STATIN) - Jlpl;

// ================= liver volume =============================================
double LVtar = (LW0 / 1.05) * (1.0 + 0.62 * fmax(gly / GLYC0 - 1.0, 0.0)
               + 0.30 * fmax(POS(TGliv) / 25.0 - 1.0, 0.0)) * (1.0 + 0.8 * POS(HCA));
dxdt_LV = 0.0025 * (LVtar - POS(LV));

// ================= purine / urate ===========================================
double prpptar = pow(Jppp / fmax(JPPPREF, 1e-9), NPRPP);
dxdt_PRPP = KPRPP * prpptar - KPRPPUSE * POS(PRPP);
double ATPdeg = KATPDEG * fmax(g6p - G6PTHATP, 0.0);
double Eallo  = EMAXALLO * HILL(POS(OXY), IC50OXY, 1.0);
double UAprod = FUAPROD * BW / 10.0 * (0.45 + 0.75 * POS(PRPP) + ATPdeg) * (1.0 - Eallo);
// urate is a QUOTIENT: allopurinol scales the numerator, lactate and ketones
// scale the denominator.  The two levers therefore multiply, not add.
double UAcl = FUACL * BW / 10.0 * POS(GFRrel)
              / (1.0 + lac / KILACUA + POS(KBOH) / KIKBUA);
dxdt_UA = (UAprod - UAcl * POS(UA) * 10.0) / (VUA * 10.0);
dxdt_ALLOg = -KAALLO * POS(ALLOg);
dxdt_ALLOc = KAALLO * POS(ALLOg) / (FVDALLO * BW / 10.0) - KALLO * POS(ALLOc);
dxdt_OXY   = FOXY * KALLO * POS(ALLOc) - KOXY * POS(OXY);

// ================= kidney ===================================================
double Eacei = EMAXACEI * HILL(POS(ACEIc), IC50ACEI, 1.0);
double GFRtar = 1.0 + 0.30 * HILL(fmax(lac - 2.0, 0.0), 4.0, 1.0) - 0.9 * POS(Rdam);
dxdt_GFRrel = CHRONIC * KHF * (GFRtar - POS(GFRrel));
double drive = WLAC * fmax(lac - 2.0, 0.0) + WHF * fmax(POS(GFRrel) - 1.0, 0.0)
               + WUA * fmax(POS(UA) - 5.5, 0.0) + WTG * POS(TGpl);
dxdt_Rdam = CHRONIC * (KDAM * drive * (1.0 - Eacei) - KREP * POS(Rdam));
dxdt_UACR = CHRONIC * KUACR * (POS(Rdam) * 2600.0 * (1.0 - Eacei) - POS(UACR));
dxdt_CitU = CHRONIC * KCIT * (pow(hco3 / HCO3N, 2.2) * (1.0 + 0.6 * CITRATE) - POS(CitU));
dxdt_NCa  = CHRONIC * KNCA * fmax(1.0 - POS(CitU), 0.0) * (1.0 + 0.4 * fmax(POS(UA) - 6.0, 0.0));

// ================= long-term liver ==========================================
// One scalar index of chronic metabolic control, entered supralinearly: brief
// severe excursions therefore count for more than their duration alone.
double MCI = fmin(1.0, (0.34 * fmax(lac - 2.2, 0.0) / 5.0
                        + 0.26 * fmax(POS(TGpl) - 250.0, 0.0) / 900.0
                        + 0.22 * fmax(POS(UA) - 5.5, 0.0) / 5.0
                        + 0.18 * fmax(3.9 - gp, 0.0) / 2.0) * 1.9);
// logistic, not exponential: adenoma burden is a fraction of liver volume and
// cannot exceed it.  Left unbounded this term reaches absurd values by year 30.
dxdt_HCA    = CHRONIC * (KHCAINI + KHCAGRO * POS(HCA)
                         * fmax(1.0 - POS(HCA) / HCAMAX, 0.0)) * pow(MCI, NMCI) * 4.0;
dxdt_HCChaz = CHRONIC * KHCC * POS(HCA) * (1.0 + 2.0 * MCI);

// ================= growth / bone ============================================
dxdt_IGF1 = KIGF * (IGF1N * (1.0 - 0.55 * MCI)
                    * (0.55 + 0.45 * HILL(gp, 3.5, 3.0)) - POS(IGF1));
double acid = fmax(HCO3N - hco3, 0.0) / 8.0;
// Relaxation towards a control-dependent TARGET, not a free integral, and
// linear growth stops at skeletal maturity whatever the metabolic state.  A
// height SDS that integrates a constant negative rate for 30 years reaches -8,
// which is not a number a human being can have.
double agenow = AGE + SOLVERTIME / 8766.0;
double growing = 1.0 / (1.0 + exp((agenow - AGEMAT) / 1.5));
dxdt_HtSDS = CHRONIC * KGRO * growing * ((HTSDSTAR - HTSDSMCI * MCI) - HtSDS);
dxdt_BMDz  = CHRONIC * KBMD * ((BMDZTAR - BMDZACID * acid - BMDZMCI * MCI) - BMDz);

// ================= exposure integrators (smoothed) ==========================
dxdt_TimeHypo = 1.0 / (1.0 + exp(6.0 * (gp - 3.9)));
dxdt_TimeNGC  = 1.0 / (1.0 + exp(120.0 * (FAI - FAITHR)));
dxdt_AUCdef   = fmax(DCNS * FAITHR - (supLac + gotGlc), 0.0);

// ================= gene / mRNA therapy ======================================
// THE durability term.  AAV episomes do not replicate, so a growing liver
// dilutes them; in a small child this dominates every other loss process and
// it cannot be answered by redosing, because anti-capsid antibody forbids it.
double aged = AGE * 8766.0 + SOLVERTIME;
double LWt  = LWADULT - (LWADULT - LW0) * exp(-(aged - AGE * 8766.0) / TAUGROW);
double grow = fmax((LWADULT - LWt) / TAUGROW / fmax(LWt, 1.0), 0.0);
dxdt_AAVvg = -(grow + KSIL) * POS(AAVvg) * CHRONIC;
dxdt_MRNAp = -KMRNA * POS(MRNAp);
dxdt_G6Pact = KEXPR * (EPSAAV * POS(AAVvg) + EPSMRNA * POS(MRNAp) - POS(G6Pact));
dxdt_AntiAAV = -KABEL * POS(AntiAAV) * CHRONIC;

// ================= GSD Ib neutrophil arm ====================================
double empa = 1.0 + (EMPAEMAX - 1.0) * HILL(POS(EMPAc), EMPAEC50, 1.0);
dxdt_AG15 = AG15IN - KAG15 * empa * POS(AG15) * POS(GFRrel);
// 1,5-AG is phosphorylated by neutrophil hexokinase; only G6PT can export the
// phosphate ester.  In Ib it accumulates and competitively starves hexokinase.
dxdt_AG6Pn = KAGUP * POS(AG15) / 14.0 - KAG6POUT * FG6PT * POS(AG6Pn)
             - KAG6PDEG * POS(AG6Pn);
dxdt_ANC = ANCPROD * (1.0 - HILL(POS(AG6Pn), KAG6P, HAG6P)) * (1.0 + 1.8 * GCSF)
           - KANC * POS(ANC);
dxdt_EMPAc = -KEMPA * POS(EMPAc);
dxdt_ACEIc = -KACEI * POS(ACEIc);

$TABLE
double LWgT = fmax(POS(LV) * 1.05, 1.0);
double supLacT = VLACCNS * (FCNSA + FCNSB * exp(-AGE / FCNSTAU))
                 * ((GURA + GURB * exp(-AGE / GURTAU)) * BW / 3.0)
                 * HILL(POS(Lac), KLACCNS, 1.0);
double DCNST = (FCNSA + FCNSB * exp(-AGE / FCNSTAU))
               * ((GURA + GURB * exp(-AGE / GURTAU)) * BW / 3.0);
double gotGlcT = fmax(DCNST - supLacT, 0.0) * HILL(POS(Gp), KB, 1.0);

capture GLUCOSE     = POS(Gp);
capture GLUCOSE_MGDL = POS(Gp) * 18.0;
capture LACTATE     = POS(Lac);
capture URATE       = POS(UA);
capture TG          = POS(TGpl);
capture BOHB        = POS(KBOH);
capture BICARB      = POS(HCO3);
capture GLYCOGEN_MGG = POS(Glyc) * 0.162;
capture LIVER_RATIO = POS(LV) / ((FLWA + FLWB * exp(-AGE / FLWTAU)) * BW * 1000.0 / 1.05);
capture G6PASE_ACT  = fmin(((GENO == 0) ? 1.0 : ((GENO == 2) ? fmax(RESID, 0.03) : RESID))
                           + POS(G6Pact), 1.0);
capture FUEL_INDEX  = (DCNST > 0.0) ? (supLacT + gotGlcT) / DCNST : 1.0;
capture CNS_LAC_SHARE = (DCNST > 0.0) ? supLacT / DCNST : 0.0;
capture ANION_GAP   = 12.0 + fmax(24.0 - POS(HCO3), 0.0);
capture NEUTROPHILS = POS(ANC);
capture AG15_UGML   = POS(AG15);
capture UACR_MGG    = POS(UACR);
capture EGFR_REL    = POS(GFRrel);
capture ADENOMA     = POS(HCA);
capture HEIGHT_SDS  = HtSDS;
capture BMD_Z       = BMDz;
capture HOURS_HYPO  = TimeHypo;
capture HOURS_NGC   = TimeNGC;
'

mod_gsd1a <- mcode_cache("gsd1a", gsd1a_code)

## =============================================================================
##  PATIENT PROFILES
## =============================================================================
##  Glucose demand per kilogram falls ~3-fold from infancy to adulthood while
##  the cornstarch dose per kilogram stays roughly constant.  That single fact
##  is why an infant needs a pump and an adult can sleep through the night.

gsd1a_patient <- function(age = 1, bw = NULL, geno = 1, resid = 0) {
  if (is.null(bw)) {
    bw <- if (age < 1) 4 + 6 * age
          else if (age < 10) 9 + 2.1 * age
          else if (age < 18) 30 + 2.4 * (age - 10)
          else 68
  }
  list(AGE = age, BW = bw, GENO = geno, RESID = resid)
}

PATIENTS <- list(
  infant_6mo = gsd1a_patient(0.5,  7),
  child_1y   = gsd1a_patient(1.0, 10),
  child_5y   = gsd1a_patient(5.0, 18),
  child_8y   = gsd1a_patient(8.0, 26),
  adol_14y   = gsd1a_patient(14.0, 50),
  adult_30y  = gsd1a_patient(30.0, 70),
  ib_child   = gsd1a_patient(12.0, 38, geno = 2),
  control_1y = gsd1a_patient(1.0, 10, geno = 0, resid = 1)
)

## Cornstarch dose -> mmol glucose equivalents.
## 1 g cornstarch ~ 88 % starch x 1.11 g glucose per g starch ~ 0.97 g glucose.
cornstarch_mmol <- function(g_per_kg, bw) g_per_kg * bw * 0.97 * 1000 / 180

## Bring the slow pools (glycogen, liver volume, hepatic and plasma TG, urate)
## to a settled repeating-day state before any acute experiment is run on top.
gsd1a_equilibrate <- function(patient, delivery_frac = 0.80, amplitude = 0.0,
                              hours = 900, ...) {
  p <- c(patient, list(RDRIPAMP = amplitude), list(...))
  GUR <- (2.10 + 4.60 * exp(-patient$AGE / 7)) * patient$BW / 3
  p$RDRIP <- delivery_frac * GUR
  out <- mod_gsd1a |>
    param(p) |>
    mrgsim(end = hours, delta = 1, atol = 1e-8, rtol = 1e-6)
  list(state = as.list(tail(as.data.frame(out), 1)), param = p, GUR = GUR)
}

## =============================================================================
##  SCENARIOS
## =============================================================================
##  Each returns a data frame.  The point of the set is that the model is asked
##  the questions clinicians actually argue about, not just "does it run".

## ---- 1. Untreated natural history: what the disease does on its own ---------
scn_natural_history <- function(patient = PATIENTS$child_1y, hours = 12) {
  eq <- gsd1a_equilibrate(patient, 0.80)
  mod_gsd1a |> param(eq$param) |> param(RDRIP = 0) |>
    init(eq$state[names(eq$state) %in% names(init(mod_gsd1a))]) |>
    mrgsim(end = hours, delta = 0.05)
}

## ---- 2. Uncooked cornstarch, the standard regimen ---------------------------
scn_uccs <- function(patient = PATIENTS$child_5y, dose_gkg = 1.6,
                     interval = 4, hours = 24) {
  eq <- gsd1a_equilibrate(patient, 0.80)
  ev <- ev(time = seq(0, hours, by = interval), cmt = "AST",
           amt = cornstarch_mmol(dose_gkg, patient$BW))
  mod_gsd1a |> param(eq$param) |> param(RDRIP = 0) |>
    init(eq$state[names(eq$state) %in% names(init(mod_gsd1a))]) |>
    mrgsim(events = ev, end = hours, delta = 0.05)
}

## ---- 3. Extended-release waxy-maize starch (Glycosade) ----------------------
##  The mechanism that matters is kdis, not the dose.  Coverage ends when the
##  RELEASE RATE falls below the deficit, so a slower starch buys more time than
##  a bigger one: doubling the dose adds ln2/kdis hours and nothing more.
scn_glycosade <- function(patient = PATIENTS$adol_14y, dose_gkg = 2.0, hours = 12) {
  eq <- gsd1a_equilibrate(patient, 0.80)
  ev <- ev(time = 0, cmt = "AST", amt = cornstarch_mmol(dose_gkg, patient$BW))
  mod_gsd1a |> param(eq$param) |>
    param(RDRIP = 0, KDIS = 0.28, FABS = 0.78) |>
    init(eq$state[names(eq$state) %in% names(init(mod_gsd1a))]) |>
    mrgsim(events = ev, end = hours, delta = 0.05)
}

## ---- 4. Continuous nocturnal drip, and what happens when the pump stops -----
scn_drip_failure <- function(patient = PATIENTS$child_5y, rate_mgkgmin = 7,
                             fail_at = 5, hours = 12) {
  eq <- gsd1a_equilibrate(patient, 0.80)
  R <- rate_mgkgmin * patient$BW / 3 * 0.98
  st <- eq$state[names(eq$state) %in% names(init(mod_gsd1a))]
  ## RDRIP is a PARAMETER, not a dosing compartment, so the pump stopping
  ## cannot be expressed as an event -- the run has to be spliced in two.
  pre <- as.data.frame(
    mod_gsd1a |> param(eq$param) |> param(RDRIP = R) |> init(st) |>
      mrgsim(end = fail_at, delta = 0.05))
  st2 <- as.list(tail(pre, 1))
  post <- as.data.frame(
    mod_gsd1a |> param(eq$param) |> param(RDRIP = 0) |>
      init(st2[names(st2) %in% names(init(mod_gsd1a))]) |>
      mrgsim(end = hours - fail_at, delta = 0.05))
  post$time <- post$time + fail_at
  rbind(pre, post[-1, ])
}

## ---- 5. Intercurrent illness with vomiting, +/- IV dextrose rescue ----------
scn_illness <- function(patient = PATIENTS$child_5y, rescue = TRUE, hours = 24) {
  eq <- gsd1a_equilibrate(patient, 0.80)
  ev1 <- ev(time = seq(0, 8, by = 4), cmt = "AST",
            amt = cornstarch_mmol(1.6, patient$BW))
  m <- mod_gsd1a |> param(eq$param) |> param(RDRIP = 0) |>
    init(eq$state[names(eq$state) %in% names(init(mod_gsd1a))])
  ## illness raises catecholamines and cortisol AND blocks enteral absorption
  m <- m |> param(ILLNESS = 1, FABS = 0.10)
  if (rescue) m <- m |> param(RDRIP = 8 * patient$BW / 3)
  mrgsim(m, events = ev1, end = hours, delta = 0.05)
}

## ---- 6. Allopurinol vs intensified diet for hyperuricaemia ------------------
##  Urate is a quotient, so these two levers multiply rather than add and the
##  combination is markedly sub-additive.  A clinician who reads the shortfall
##  as non-adherence will be wrong.
scn_urate <- function(patient = PATIENTS$adol_14y, allopurinol_mg = 300,
                      intensify = FALSE, days = 28) {
  eq <- gsd1a_equilibrate(patient, if (intensify) 0.92 else 0.85,
                          amplitude = if (intensify) 0.20 else 0.95)
  ev <- if (allopurinol_mg > 0)
          ev(time = seq(0, days * 24, by = 24), cmt = "ALLOg", amt = allopurinol_mg)
        else ev(time = 0, cmt = "ALLOg", amt = 0)
  mod_gsd1a |> param(eq$param) |>
    init(eq$state[names(eq$state) %in% names(init(mod_gsd1a))]) |>
    mrgsim(events = ev, end = days * 24, delta = 2)
}

## ---- 7. Lipid therapy -------------------------------------------------------
scn_lipid <- function(patient = PATIENTS$adol_14y, fibrate = 1, statin = 1,
                      days = 90) {
  eq <- gsd1a_equilibrate(patient, 0.85, amplitude = 0.5)
  mod_gsd1a |> param(eq$param) |> param(FIBRATE = fibrate, STATIN = statin) |>
    init(eq$state[names(eq$state) %in% names(init(mod_gsd1a))]) |>
    mrgsim(end = days * 24, delta = 12)
}

## ---- 8. ACE inhibitor for microalbuminuria (10 years) -----------------------
scn_renal <- function(patient = PATIENTS$child_8y, acei = TRUE, years = 20) {
  eq <- gsd1a_equilibrate(patient, 0.88, amplitude = 0.6)
  m <- mod_gsd1a |> param(eq$param) |>
    init(eq$state[names(eq$state) %in% names(init(mod_gsd1a))])
  ## the renal endpoint responds to AVERAGE ACE inhibition over years, so the
  ## drug is entered at its steady-state average rather than as 7000 doses
  if (acei) m <- m |> param(KACEI = 0) |> init(ACEIc = 55 / (0.058 * 24))
  mrgsim(m, end = years * 8766, delta = 24 * 30)
}

## ---- 9. DTX401 gene therapy at different ages -------------------------------
##  The headline is NOT the peak activity — that is the same at every age.  It
##  is the RETENTION, which is set by how much liver growth remains.  And the
##  decision is irreversible: anti-capsid antibody rules out a second dose.
scn_gene_therapy <- function(patient = PATIENTS$child_5y, load = 0.22,
                             years = 20) {
  eq <- gsd1a_equilibrate(patient, 0.85)
  st <- eq$state[names(eq$state) %in% names(init(mod_gsd1a))]
  st$AAVvg <- load
  st$AntiAAV <- 1
  mod_gsd1a |> param(eq$param) |> init(st) |>
    mrgsim(end = years * 8766, delta = 24 * 30)
}

## ---- 10. mRNA-3745: repeat dosing, NOT diluted by growth --------------------
scn_mrna <- function(patient = PATIENTS$child_5y, dose = 0.22,
                     interval_d = 14, weeks = 24) {
  eq <- gsd1a_equilibrate(patient, 0.85)
  ev <- ev(time = seq(0, weeks * 168, by = interval_d * 24),
           cmt = "MRNAp", amt = dose)
  mod_gsd1a |> param(eq$param) |>
    init(eq$state[names(eq$state) %in% names(init(mod_gsd1a))]) |>
    mrgsim(events = ev, end = weeks * 168, delta = 6)
}

## ---- 11. GSD Ib: empagliflozin for neutropenia ------------------------------
scn_empagliflozin <- function(patient = PATIENTS$ib_child, dose_nM = 780,
                              days = 180, start_day = 30) {
  eq <- gsd1a_equilibrate(patient, 0.85, hours = 2000)
  ev <- ev(time = seq(start_day * 24, days * 24, by = 24),
           cmt = "EMPAc", amt = dose_nM)
  mod_gsd1a |> param(eq$param) |>
    init(eq$state[names(eq$state) %in% names(init(mod_gsd1a))]) |>
    mrgsim(events = ev, end = days * 24, delta = 6)
}

## ---- 12. Overtreatment: surplus glucose is not free -------------------------
##  A normal liver exports surplus glucose again.  A GSD Ia liver cannot, so
##  every extra gram is trapped as G6P and leaves as glycogen, lactate or fat.
scn_overtreatment <- function(patient = PATIENTS$child_8y,
                              fracs = c(0.60, 0.75, 0.90, 1.05, 1.25, 1.45)) {
  do.call(rbind, lapply(fracs, function(f) {
    eq <- gsd1a_equilibrate(patient, f, amplitude = 0.30, hours = 1200)
    data.frame(delivery_frac = f,
               glucose  = eq$state$GLUCOSE,
               lactate  = eq$state$LACTATE,
               urate    = eq$state$URATE,
               TG       = eq$state$TG,
               glycogen = eq$state$GLYCOGEN_MGG,
               liver_x  = eq$state$LIVER_RATIO)
  }))
}

## ---- 13. Thirty-year natural history, good vs poor control ------------------
scn_lifetime <- function(patient = PATIENTS$child_8y, quality = c("good", "poor"),
                         years = 30) {
  quality <- match.arg(quality)
  frac <- if (quality == "good") 0.92 else 0.85
  amp  <- if (quality == "good") 0.20 else 0.95
  eq <- gsd1a_equilibrate(patient, frac, amplitude = amp)
  mod_gsd1a |> param(eq$param) |>
    init(eq$state[names(eq$state) %in% names(init(mod_gsd1a))]) |>
    mrgsim(end = years * 8766, delta = 24 * 90)
}

## ---- 14. Residual-activity dose-response ------------------------------------
##  Daily starch requirement falls LINEARLY in restored activity while fasting
##  tolerance rises HYPERBOLICALLY, so there is a critical activity below which
##  the biochemistry improves and the patient's night does not.
scn_activity_response <- function(patient = PATIENTS$adol_14y,
                                  a = c(0, 0.01, 0.02, 0.03, 0.05, 0.08,
                                        0.12, 0.18, 0.25, 0.40)) {
  do.call(rbind, lapply(a, function(ai) {
    pat <- patient; pat$RESID <- ai
    eq <- gsd1a_equilibrate(pat, 0.85)
    st <- eq$state[names(eq$state) %in% names(init(mod_gsd1a))]
    sim <- mod_gsd1a |> param(eq$param) |> param(RDRIP = 0) |> init(st) |>
      mrgsim(events = ev(time = 0, cmt = "AST",
                         amt = cornstarch_mmol(1.45, pat$BW)),
             end = 30, delta = 0.05) |> as.data.frame()
    below <- sim$time[sim$GLUCOSE < 3.3]
    data.frame(activity = ai,
               fasting_tolerance_h = if (length(below)) min(below) else 30,
               lactate = eq$state$LACTATE,
               urate   = eq$state$URATE,
               TG      = eq$state$TG)
  }))
}

SCENARIOS <- list(
  natural_history   = scn_natural_history,
  uccs              = scn_uccs,
  glycosade         = scn_glycosade,
  drip_failure      = scn_drip_failure,
  illness           = scn_illness,
  urate             = scn_urate,
  lipid             = scn_lipid,
  renal             = scn_renal,
  gene_therapy      = scn_gene_therapy,
  mrna              = scn_mrna,
  empagliflozin     = scn_empagliflozin,
  overtreatment     = scn_overtreatment,
  lifetime          = scn_lifetime,
  activity_response = scn_activity_response
)

## =============================================================================
##  PARAMETER PROVENANCE
## =============================================================================
##  Glucose demand by age            Bier 1977 Diabetes; Haymond, Kalhan
##  Residual EGP in GSD I            Powell 1981; Tsalikian 1984; Kalderon 1989
##                                   (isotopic EGP ~1.0-1.5 vs 4.5-5.5 mg/kg/min)
##  alpha-1,6 branch fraction        7-10 % of hepatic glycogen glucosyl units
##  Lactate as cerebral fuel         van Hall 2009 J Cereb Blood Flow Metab;
##                                   Boumezbeur 2010 J Neurosci (~60 % at 7 mM)
##  Cornstarch dose and interval     Chen 1984 NEJM 310:171; Wolfsdorf 1990
##  Extended-release waxy maize      Correia 2008 Am J Clin Nutr; Bhattacharya 2007
##  LPL activity in GSD Ia           Forget 1974; Levy 1988
##  URAT1 lactate competition        Enomoto 2002 Nature; Lipkowitz 2012 review
##  Allopurinol / oxypurinol PK      Day 2007 Clin Pharmacokinet
##  Adenoma epidemiology             Franco 2005; Wang 2011; Kishnani 2014 ACMG
##  DTX401 phase 1/2 and phase 3     Weinstein 2021 Hum Gene Ther; NCT03517085;
##                                   GlucoGene NCT05139316 (cornstarch endpoint)
##  mRNA-3745                        NCT05095727
##  Empagliflozin in GSD Ib          Wortmann 2020 NEJM 383:1480;
##                                   Grunert 2022 Orphanet J Rare Dis (n = 112)
##  Management guideline             Kishnani 2014 Genet Med 16:e1 (ACMG)
##
##  Full citation list with PubMed links: gsd1a_references.md
## =============================================================================
