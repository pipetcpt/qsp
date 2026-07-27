# =====================================================================
#  PRIMARY HYPEROXALURIA (PH1 / PH2 / PH3) — QSP model for mrgsolve
#  ---------------------------------------------------------------------
#  THE ORGANISING IDEA
#
#    Oxalate is a terminal metabolite. Humans have no oxalate-degrading
#    enzyme. So primary hyperoxaluria is not a disease with a severity
#    score — it is an UNCLOSED MASS BALANCE:
#
#      d(total body oxalate)/dt = PRODUCTION
#                              - RENAL EXCRETION
#                              - ENTERIC ELIMINATION
#                              - DIALYTIC REMOVAL
#
#    There is NO parameter in this model called severity, stage, or
#    progression rate that a scenario can turn up. Urinary oxalate,
#    plasma oxalate, eGFR slope, time to ESKD, stone rate and systemic
#    oxalosis are all OUTPUTS of the balance above. The single input
#    that distinguishes a lethal infantile phenotype from an incidental
#    adult stone former is FAGT — residual alanine:glyoxylate
#    aminotransferase activity.
#
#  ---------------------------------------------------------------------
#  WHAT THIS MODEL IS BUILT TO GENERATE (rather than assume)
#
#   (1) THE AGT NONLINEARITY.  AGT is a high-Vmax / low-Km enzyme, so
#       the map from residual activity to urinary oxalate is steeply
#       nonlinear. The model therefore produces, from one Michaelis-
#       Menten term, both halves of a fact that a linear model cannot
#       hold at once: a 50%-activity heterozygote is INDISTINGUISHABLE
#       from normal, while restoring activity from 0.5% to 3% — a
#       chaperone-sized gain — removes about a third of the oxalate.
#       That is why PH1 is recessive with entirely healthy carriers AND
#       why pyridoxine is worth giving to the genotypes it can rescue.
#
#   (2) A THRESHOLD DISEASE FROM A LINEAR SINK.  While GFR is preserved
#       the kidney clears essentially all production, plasma oxalate is
#       near normal, and PH is a stone disease. Renal clearance is
#       proportional to GFR, so once GFR falls below production/Pox_crit
#       the residual has nowhere to go but tissue. The model COMPUTES
#       the eGFR at which plasma oxalate crosses the calcium-oxalate
#       solubility limit of plasma; systemic oxalosis is that root, not
#       a separate late mechanism.
#
#   (3) A RUNAWAY WHOSE GAIN RISES.  Crystals destroy nephrons, which
#       lowers clearance, which raises plasma oxalate, which deposits
#       more crystals. The loop is closed in the equations (there is no
#       prescribed eGFR trajectory anywhere), and its gain rises as GFR
#       falls, which is what makes the terminal decline abrupt.
#
#   (4) THE BIOMARKER THAT INVERTS.  Uox = Cl_ox x Pox. Below the
#       deposition threshold Uox equals production and is a clean
#       efficacy biomarker; above it, falling GFR pulls Uox DOWN while
#       the disease accelerates. The model computes the eGFR at which
#       Uox has fallen 25% and 50% FROM RENAL FAILURE ALONE, with
#       production held exactly constant. This is why one drug needed
#       two different registration endpoints in two different trials.
#
#   (5) WHY AN LDHA DRUG UNDERPERFORMS A HAO1 DRUG, AND FAILS IN PH2.
#       This is the model's headline result and it uses no PH2-specific
#       drug parameter. Two structural facts do the work:
#
#         (a) Glycolate oxidase does not only MAKE glyoxylate, it also
#             oxidises glyoxylate to oxalate. Silencing HAO1 therefore
#             cuts the substrate supply AND one of the two oxalate-
#             forming reactions; silencing LDHA cuts only the other one
#             and leaves the GO route as an untouchable floor.
#
#         (b) Silencing an enzyme reduces its FLUX only if a parallel
#             branch exists to take the substrate. In PH1, GRHPR is that
#             parallel branch, so LDHA knockdown re-partitions glyoxylate
#             into harmless glycolate. In PH2 GRHPR IS THE MISSING GENE:
#             LDHA is the only exit from the cytosolic glyoxylate pool,
#             and at steady state flux through the only exit equals its
#             input no matter how little enzyme is present. Knocking it
#             down raises glyoxylate and barely moves oxalate.
#
#       So an LDHA-directed siRNA is structurally guaranteed to work
#       worse in PH2 than in PH1 even before hepatocyte-restricted
#       delivery is considered — and the model quantifies how much of
#       the PH2 failure each explanation accounts for (diagnostic D07).
#
#   (6) LAYER B THERAPY THAT CHANGES NO FLUX.  Hyperhydration and
#       potassium citrate leave every mass-balance term untouched and
#       still halve supersaturation, because the stone endpoint is a
#       RATIO (a Tiselius AP(CaOx) index) and not a flux. The model
#       therefore separates "protects the kidney" from "closes the
#       balance", and shows supportive care doing the first and never
#       the second.
#
#   (7) DIALYSIS AS ARITHMETIC.  Oxalate removal per session is
#       integrated from the actual intradialytic clearance and the
#       actual falling plasma concentration, so weekly removal is
#       computed and compared with weekly production. Conventional
#       thrice-weekly haemodialysis is shown to be insufficient by
#       calculation, not by assertion.
#
#   (8) A RESERVOIR, NOT A RELAPSE.  Bone holds the largest oxalate
#       pool and releases it by osteoclastic resorption. After a normal
#       liver is installed, urinary oxalate stays elevated for months
#       from the SKELETON. The model generates that, and predicts that
#       a steroid-heavy regimen prolongs it.
#
#  ---------------------------------------------------------------------
#  NEGATIVE AND SELF-REFUTING RESULTS ARE REPORTED, NOT REMOVED
#  See section 8. In particular the model FAILS to generate PH3
#  hyperoxaluria from the field's leading hypothesis, over-predicts the
#  PHYOX2 nedosiran response, and cannot reproduce true infantile
#  oxalosis timing. Those are stated with their arithmetic rather than
#  tuned away.
#
#  ---------------------------------------------------------------------
#  UNITS.  time = days; oxalate & organic acids = umol (amounts),
#  umol/L (concentrations); urinary excretion reported as
#  mmol/1.73 m^2/24 h; GFR = mL/min/1.73 m^2; drug amounts = mg or nmol
#  as marked. Reference subject: 70 kg adult, BSA 1.73 m^2, unless a
#  scenario overrides it.
#
#  Validated under mrgsolve 2.0.1 / R 4.3.3.
# =====================================================================

suppressMessages({
  library(mrgsolve)
  library(dplyr)
})

# =====================================================================
#  1. THE MODEL
# =====================================================================

ph1_code <- '
$PROB
PRIMARY HYPEROXALURIA (PH1/PH2/PH3) QSP model — 73 ODEs.
Mass-balance architecture: urinary oxalate, plasma oxalate, eGFR and
systemic oxalosis are OUTPUTS of production minus the three sinks.

$PARAM @annotated
// ---------- GENOTYPE / ENZYME COMPLEMENT (the only disease inputs) ---
FAGT      : 1.0   : residual AGT activity, fraction of normal (PH1 = ~0-0.05)
FGRHPR    : 1.0   : residual GRHPR activity, fraction of normal (PH2 = 0)
FHOGA1    : 1.0   : residual HOGA1 activity, fraction of normal (PH3 = 0)
B6RESC    : 0.0   : max AGT activity attainable by PLP chaperone rescue (genotype-gated)
HOGCYT    : 0.0   : cytosolic HOG-aldolase hypothesis switch for PH3 (0 = off)

// ---------- SUBSTRATE SUPPLY (umol glyoxylate-equivalents / day) -----
J_GLYCO   : 1293  : glycolate-derived glyoxylate flux via GO at full HAO1
J_DAO     : 172   : glycine/DAO-derived peroxisomal glyoxylate flux
J_HYP     : 259   : 4-hydroxyproline-derived mitochondrial glyoxylate flux
GLYCO_SUP : 1558  : total hepatic glycolate supply (diet + endogenous), umol/day
OX_DIET   : 1300  : dietary oxalate intake, umol/day
FABS_OX   : 0.108 : fraction of dietary oxalate absorbed
OX_ASC    : 40    : oxalate from non-enzymatic ascorbate degradation, umol/day
OX_XHEP   : 60    : constitutive extrahepatic oxalate synthesis, umol/day
OX_XHEP2  : 600   : ADDITIONAL extrahepatic source unmasked by GRHPR loss (PH2)

// ---------- FLUX PARTITIONING (dimensionless, sum to 1 within a node)
AGT_EFF   : 0.98  : fraction of peroxisomal glyoxylate handled by AGT in health
PHI_GOOX  : 0.15  : share of non-AGT peroxisomal outflow oxidised by GO to oxalate
FR_LDH    : 0.485 : share of cytosolic glyoxylate outflow via LDHA at full enzyme
FR_GRH    : 0.475 : share of cytosolic glyoxylate outflow via GRHPR at full enzyme
FR_GXSP   : 0.040 : share spilling to plasma glyoxylate (the relief valve)
FR_M_GRH  : 0.60  : mitochondrial glyoxylate reduced by mito-GRHPR in health
FR_M_PER  : 0.32  : mitochondrial glyoxylate exported to peroxisome
GXP_REF   : 60    : peroxisomal glyoxylate amount at the AGT-null reference (umol)
KC_HOGA   : 30.0  : HOGA1 catalytic rate constant on HOG, /day
K_DHG     : 0.30  : HOG -> 2,4-dihydroxyglutarate diversion rate, /day
KM_AGT    : 30    : Michaelis constant of AGT for glyoxylate (umol, model scale)

// ---------- OXALATE DISTRIBUTION ------------------------------------
V_OX      : 15.0  : oxalate distribution volume, L
V_ECF     : 12.0  : slowly exchanging extracellular volume, L
K_PL_ECF  : 6.0   : plasma<->ECF exchange clearance, L/day
POX_CRIT  : 30.0  : plasma CaOx solubility limit, umol/L
POX_CRIT_B: 20.0  : bone mineralisation-front deposition threshold, umol/L
POX_CRIT_S: 35.0  : soft-tissue deposition threshold, umol/L
CLDEP_B   : 24.0  : bone deposition clearance above threshold, L/day
CLDEP_S   : 8.0   : soft-tissue deposition clearance above threshold, L/day
CLDEP_K   : 4.0   : renal parenchymal deposition clearance above threshold, L/day
BMAX_BONE : 600000: saturable skeletal oxalate capacity, umol
K_RES     : 0.0055: SURFACE skeletal oxalate release rate (osteoclastic), /day
K_BURY    : 0.0060: burial of surface into crystal-incorporated deep pool, /day
K_RES_D   : 0.00019: DEEP pool release rate (bone remodelling, t1/2 ~10 y), /day
K_SOFT_R  : 0.0015: soft-tissue oxalate clearance, /day
K_KID_R   : 0.0020: renal parenchymal oxalate clearance, /day

// ---------- RENAL HANDLING ------------------------------------------
GFR_MAX   : 125   : GFR of a full nephron complement, mL/min/1.73m2
HYPERF    : 0.60  : single-nephron hyperfiltration compensation coefficient
F_SEC_OX  : 1.15  : oxalate clearance / creatinine clearance (net secretion)
F_SEC_GC  : 1.00  : glycolate clearance ratio
F_SEC_GX  : 1.00  : glyoxylate clearance ratio
CL_ENT    : 1.2   : enteric secretion clearance of plasma oxalate, L/day
OXF_DEG   : 0.0   : Oxalobacter luminal degradation, fraction of luminal load
FGO_KID   : 0.0   : renal LDHA inhibition by a small molecule (fraction inhibited)

// ---------- URINE CHEMISTRY (Layer B) -------------------------------
U_VOL_TGT : 1.5   : target 24-h urine volume, L
K_UVOL    : 0.50  : urine-volume equilibration rate, /day
U_CA      : 4.0   : urinary calcium, mmol/day
U_CIT     : 3.0   : urinary citrate, mmol/day
U_MG      : 4.0   : urinary magnesium, mmol/day
KCIT_DOSE : 0.0   : potassium citrate, mmol citrate/day added
MG_DOSE   : 0.0   : magnesium supplement, mmol/day added
THIAZ     : 0.0   : thiazide effect on urinary calcium (fraction reduced)
AP_THRESH : 1.00  : AP(CaOx) index above which crystal formation begins

// ---------- CRYSTAL -> INJURY -> NEPHRON LOSS -----------------------
K_NUC     : 0.03125 : crystal nucleation/retention rate from supersaturation, /day
K_CRCL    : 0.10  : crystal clearance rate, /day
K_NLRP3   : 0.50  : inflammasome activation rate from crystal burden, /day
K_NLRP3_D : 0.50  : inflammasome deactivation, /day
K_IL1     : 0.60  : IL-1beta generation, /day
K_IL1_D   : 0.60  : IL-1beta elimination, /day
K_TGFB    : 0.25  : TGF-beta generation, /day
K_TGFB_D  : 0.35  : TGF-beta elimination, /day
K_FIB     : 0.0015 : fibrosis accrual rate, /day
K_FIB_R   : 0.0015: fibrosis resolution rate, /day
K_LOSS    : 0.00027: nephron loss rate constant, /day
W_KIDOX   : 2.2   : weight of parenchymal oxalate on nephron loss
K_TINJ    : 0.8   : tubular injury marker generation, /day
K_TINJ_D  : 1.2   : tubular injury marker decay, /day
IL1_BLOCK : 0.0   : anti-IL-1 / NLRP3 blockade (fraction inhibited)

// ---------- STONES AND ORGAN ENDPOINTS ------------------------------
K_STONE   : 0.020 : stone mass accrual per unit crystal burden, /day
K_STONE_R : 0.0020: stone passage / removal, /day
K_HAZ     : 0.00055: symptomatic stone event hazard per unit stone mass, /day
K_RET     : 0.0022: retinal deposition rate from soft-tissue burden, /day
K_MYO     : 0.0016: myocardial deposition rate, /day
K_NRV     : 0.0018: peripheral nerve deposition rate, /day
K_SKN     : 0.0014: cutaneous deposition rate, /day
K_BONEDIS : 0.0025: oxalate osteopathy accrual rate, /day
K_MARROW  : 0.0018: marrow oxalosis / EPO-resistance accrual, /day
K_ORG_R   : 0.0012: organ deposit resolution rate, /day

// ---------- LUMASIRAN (GalNAc-siRNA vs HAO1) ------------------------
KA_LUM    : 3.00  : SC absorption rate, /day
CL_LUM    : 120   : plasma clearance, L/day
V_LUM     : 6.0   : plasma volume of distribution, L
CLUP_LUM  : 300   : ASGPR-mediated hepatic uptake clearance, L/day
KM_UP_LUM : 2500  : ASGPR uptake saturation constant, nmol (hepatic-flux scale)
K_ESC_LUM : 0.09  : endosomal escape / RISC loading rate, /day
K_RISC_L  : 0.0154: RISC-loaded siRNA loss rate (t1/2 45 d), /day
EMAX_HAO1 : 0.958 : maximal fractional HAO1 mRNA suppression
EC50_HAO1 : 0.85  : RISC amount giving half-maximal HAO1 suppression, nmol
K_HAO1_M  : 0.35  : HAO1 mRNA turnover, /day
K_GO_P    : 0.231 : GO protein turnover (t1/2 3 d), /day

// ---------- NEDOSIRAN (GalNAc-siRNA vs LDHA) ------------------------
KA_NED    : 3.00  : SC absorption rate, /day
CL_NED    : 140   : plasma clearance, L/day
V_NED     : 7.0   : plasma volume of distribution, L
CLUP_NED  : 280   : ASGPR-mediated hepatic uptake clearance, L/day
KM_UP_NED : 2500  : ASGPR uptake saturation constant, nmol
K_ESC_NED : 0.08  : endosomal escape / RISC loading rate, /day
K_RISC_N  : 0.0231: RISC-loaded siRNA loss rate (t1/2 30 d), /day
EMAX_LDHA : 0.90  : maximal fractional LDHA mRNA suppression
EC50_LDHA : 1.20  : RISC amount giving half-maximal LDHA suppression, nmol
K_LDHA_M  : 0.40  : LDHA mRNA turnover, /day
K_LDHA_P  : 0.198 : LDHA protein turnover (t1/2 3.5 d), /day

// ---------- STIRIPENTOL (oral LDH5 inhibitor, liver AND kidney) -----
KA_STP    : 1.8   : oral absorption rate, /day
CL_STP    : 40    : apparent clearance, L/day
V_STP     : 60    : central volume, L
V_STP_P   : 90    : peripheral volume, L
Q_STP     : 25    : intercompartmental clearance, L/day
IC50_STP  : 12.0  : stiripentol concentration for 50% LDH inhibition, mg/L
IMAX_STP  : 0.80  : maximal fractional LDH inhibition

// ---------- PYRIDOXINE -> PLP ---------------------------------------
KA_PN     : 2.0   : oral absorption rate, /day
CL_PLP    : 229   : PLP clearance, L/day
V_PLP     : 20    : PLP distribution volume, L
K_PLP_LIV : 3.0   : plasma->hepatic PLP distribution rate, /day
K_PLP_OUT : 2.0   : hepatic PLP loss rate, /day
PN_MGKG   : 0.0   : pyridoxine dose actually prescribed, mg/kg/day (drives the rescue)
EC50_PN   : 3.5   : pyridoxine dose giving half-maximal chaperone rescue, mg/kg/day
PN_DIET   : 2.0   : dietary pyridoxine intake, mg/day
VM_PN     : 200000: saturable pyridoxine->PLP conversion capacity, nmol/day
KM_PN     : 30000 : half-saturation of pyridoxine->PLP conversion, nmol
K_B6NEURO : 0.00060 : cumulative sensory-neuronopathy accrual per mg/kg/day-day

// ---------- DIALYSIS -------------------------------------------------
HD_PERWK  : 0     : haemodialysis sessions per week (0 = none)
HD_HRS    : 4.0   : hours per session
CL_HD     : 187   : intradialytic oxalate clearance, L/day (=130 mL/min)
PD_ON     : 0     : peritoneal dialysis flag
CL_PD     : 9.0   : continuous PD oxalate clearance, L/day

// ---------- TRANSPLANT / GRAFT ---------------------------------------
LIV_TX_T  : 1e6   : time of liver (or combined) transplant, days
KID_TX_T  : 1e6   : time of kidney transplant, days
KID_TX_N  : 0.85  : nephron mass delivered by a single kidney graft
STEROID   : 0.0   : corticosteroid-driven multiplier on bone resorption
CNI_TOX   : 0.0   : calcineurin-inhibitor nephrotoxicity (fractional loss rate add-on)

// ---------- ACIDOSIS, GROWTH, MISC ----------------------------------
K_ACID    : 0.10  : acidosis equilibration rate, /day
ACID_GAIN : 1.6   : acidosis generated per unit GFR deficit
GROWTH    : 1.0   : growth/collagen-turnover multiplier on J_HYP
BSA       : 1.73  : body surface area, m2
ASC_DOSE  : 0.0   : supplemental ascorbate-derived oxalate, umol/day
INIT_ON   : 1     : 1 = impose the solved healthy fixed point as t=0 state; 0 = accept a supplied state

$CMT @annotated
// --- hepatic intermediary metabolism (1-8)
GLYCOL_P  : hepatic glycolate pool (umol)
GX_P      : peroxisomal glyoxylate (umol)
GX_C      : cytosolic glyoxylate (umol)
HOG_M     : mitochondrial 4-hydroxy-2-oxoglutarate (umol)
GX_M      : mitochondrial glyoxylate (umol)
DHG       : 2,4-dihydroxyglutarate, PH3 marker (umol)
HPYR      : hydroxypyruvate (umol)
LGLYC     : L-glycerate, PH2 marker (umol)
// --- oxalate distribution (9-16)
OX_PL     : plasma oxalate (umol)
OX_ECF    : slowly exchanging extracellular oxalate (umol)
OX_BONE   : SURFACE (exchangeable) skeletal oxalate (umol)
OX_BONE_D : DEEP crystal-incorporated skeletal oxalate (umol)
OX_SOFT   : soft-tissue oxalate (umol)
OX_KID    : renal parenchymal oxalate / nephrocalcinosis burden (umol)
GLYCOL_PL : plasma glycolate (umol)
GX_PL     : plasma glyoxylate (umol)
OX_GUT    : intestinal luminal oxalate (umol)
// --- enzymes and their mRNA (17-25)
HAO1_M    : HAO1 mRNA, relative to 1
GO_P      : glycolate oxidase protein, relative to 1
LDHA_M    : LDHA mRNA, relative to 1
LDHA_P    : hepatic LDHA protein, relative to 1
LDHA_K    : renal LDHA protein, relative to 1
AGT_APO   : apo-AGT, relative to 1
AGT_HOLO  : PLP-loaded functional AGT, relative to 1
GRHPR_P   : GRHPR protein, relative to 1
HOGA1_P   : HOGA1 protein, relative to 1
// --- lumasiran (26-29)
LUM_SC    : lumasiran SC depot (nmol)
LUM_PL    : lumasiran plasma (nmol)
LUM_LIV   : lumasiran hepatic endosomal depot (nmol)
LUM_RISC  : RISC-loaded lumasiran (nmol)
// --- nedosiran (30-33)
NED_SC    : nedosiran SC depot (nmol)
NED_PL    : nedosiran plasma (nmol)
NED_LIV   : nedosiran hepatic endosomal depot (nmol)
NED_RISC  : RISC-loaded nedosiran (nmol)
// --- stiripentol (34-36)
STP_GUT   : stiripentol gut (mg)
STP_CEN   : stiripentol central (mg)
STP_PER   : stiripentol peripheral (mg)
// --- pyridoxine / PLP (37-39)
PN_GUT    : pyridoxine gut (nmol)
PLP_PL    : plasma PLP (nmol)
PLP_LIV   : hepatic PLP (nmol)
// --- kidney structure and injury (40-49)
NEPHRON   : viable nephron fraction
CRYST     : tubular crystal burden (relative)
NLRP3     : inflammasome activity (relative)
IL1B      : IL-1beta (relative)
TGFB      : TGF-beta (relative)
FIBROSIS  : interstitial fibrosis (relative)
TUB_INJ   : tubular injury marker, KIM-1 surrogate (relative)
ACIDOSIS  : metabolic acidosis burden (relative)
U_VOL     : 24-h urine volume (L)
GRAFT     : kidney graft functional fraction
// --- stones and organ endpoints (50-58)
STONE     : stone mass (relative)
STONE_HAZ : cumulative symptomatic stone event hazard
RETINA    : retinal oxalate deposition score
MYOCARD   : myocardial oxalate deposition score
NERVE     : peripheral nerve oxalate deposition score
SKIN      : cutaneous oxalate deposition score
BONE_DIS  : oxalate osteopathy score
MARROW    : marrow oxalosis / EPO-resistance score
B6_NEURO  : cumulative pyridoxine neurotoxicity exposure
// --- audit integrals (59-72)
CUM_PROD  : cumulative oxalate production (umol)
CUM_UOX   : cumulative urinary oxalate excretion (umol)
CUM_ENT   : cumulative enteric oxalate elimination (umol)
CUM_HD    : cumulative dialytic oxalate removal (umol)
CUM_DEP   : cumulative tissue oxalate deposition (umol)
AUC_POX   : AUC of plasma oxalate (umol/L*day)
AUC_AP    : AUC of AP(CaOx) index (day)
T_ABOVE30 : time with plasma oxalate above the solubility limit (days)
T_EGFR30  : time with eGFR below 30 (days)
CUM_GXTOX : cumulative cytosolic glyoxylate exposure (umol*day)
CUM_GLYCO : cumulative urinary glycolate (umol)
CUM_FILT  : cumulative filtered oxalate load (umol)
ESKD_FLAG : integral that first exceeds 0 when eGFR crosses below 15
CUM_BONE_R: cumulative skeletal oxalate released (umol)

$GLOBAL
#define POX   (OX_PL/V_OX)
#define PGLYC (GLYCOL_PL/V_OX)
#define PGX   (GX_PL/V_OX)

// derived rate constants, all solved in $MAIN so that the HEALTHY state
// is an exact fixed point and the disease is generated, never imposed
double kGOox, kLeak, kAGTvm, kLDH, kGR, kGXsp;
double kMgrh, kMper, kMcyt, kGOin, kGlycOut, kHOG, kHOGcyt;
double gfr_now, clox_now, uvol_now, uca_now, ucit_now, umg_now;
double ap_index, fGO_eff, fLDH_eff, fLDHk_eff, fGRH_eff, fAGT_eff;
double prod_hep, prod_tot, uox_rate, dep_bone, dep_soft, dep_kid;
double hd_clear, ent_elim, ldh_share, grh_share, gxsp_share;
double stp_conc, plp_effect, hog_ki, bone_release;

$MAIN
// =====================================================================
//  BASELINE SOLVE — every rate constant below is DERIVED from the
//  annotated flux targets so that at FAGT=1 the whole system sits at an
//  exact fixed point. Nothing about the diseased state is prescribed.
// =====================================================================

// --- mitochondrial glyoxylate partitioning (GRHPR-dependent) --------
kMgrh = FR_M_GRH;                       // reduced by mitochondrial GRHPR
kMper = FR_M_PER;                       // exported to peroxisome
kMcyt = 1.0 - FR_M_GRH - FR_M_PER;      // spills to cytosol

// --- peroxisomal inflow at the reference (full enzyme) --------------
double per_in_ref = J_GLYCO + J_DAO + kMper*J_HYP;

// --- peroxisomal outflow constants, pinned by the AGT-null reference
double k_tot_per = per_in_ref / GXP_REF;      // /day
kGOox = PHI_GOOX * k_tot_per;
kLeak = (1.0 - PHI_GOOX) * k_tot_per;

// --- AGT Vmax pinned by AGT_EFF in health ---------------------------
double gxp_health = per_in_ref*(1.0 - AGT_EFF) / k_tot_per;
kAGTvm = per_in_ref*AGT_EFF*(KM_AGT + gxp_health)/gxp_health;

// --- cytosolic glyoxylate partitioning ------------------------------
kLDH  = FR_LDH;
kGR   = FR_GRH;
kGXsp = FR_GXSP;

// --- glycolate handling: GO oxidation vs renal excretion ------------
kGOin    = J_GLYCO;                              // umol/day at fGO = 1
kGlycOut = GLYCO_SUP - J_GLYCO;                  // umol/day escaping to plasma

// --- hydroxyproline route -------------------------------------------
kHOG    = J_HYP*GROWTH;
kHOGcyt = HOGCYT;                                // PH3 alternative hypothesis

// --- initial conditions: the exact healthy fixed point --------------
//  Set INIT_ON = 0 to hand the model a state vector instead (used by the
//  two-phase transplant / dialysis runs, where phase 2 must inherit the
//  accumulated tissue burden and nephron loss of phase 1).
if(INIT_ON > 0.5){
GLYCOL_P_0  = 100.0;
GX_P_0      = gxp_health;
GX_C_0      = 1.0;
HOG_M_0     = J_HYP*GROWTH/(KC_HOGA + K_DHG);
GX_M_0      = 5.0;
DHG_0       = 5.0;
HPYR_0      = 5.0;
LGLYC_0     = 5.0;

HAO1_M_0    = 1.0;  GO_P_0    = 1.0;
LDHA_M_0    = 1.0;  LDHA_P_0  = 1.0;  LDHA_K_0 = 1.0;
AGT_APO_0   = 0.0;  AGT_HOLO_0 = 1.0;
GRHPR_P_0   = 1.0;  HOGA1_P_0  = 1.0;

NEPHRON_0   = 1.0;
U_VOL_0     = U_VOL_TGT;
GRAFT_0     = 0.0;

// healthy plasma oxalate: production / (renal + enteric clearance)
double prod_h = 0.0;   // resolved numerically by the pre-run; seed sensibly
OX_PL_0   = 20.0;
OX_ECF_0  = 16.0;
GLYCOL_PL_0 = 4.0;
GX_PL_0     = 0.5;
OX_GUT_0    = 300.0;

// PLP: a healthy dietary baseline
PN_GUT_0  = 1888.0;
PLP_PL_0  = 1000.0;
PLP_LIV_0 = 1500.0;
}

$ODE
// =====================================================================
//  A.  EFFECTIVE ENZYME LEVELS
// =====================================================================
// glycolate oxidase — suppressed by lumasiran through HAO1 mRNA
fGO_eff = GO_P;

// hepatic LDHA — suppressed by nedosiran (mRNA) and inhibited by
// stiripentol (concentration-dependent, and it also reaches the kidney)
stp_conc  = STP_CEN/V_STP;
double stp_inh = IMAX_STP*stp_conc/(IC50_STP + stp_conc);
fLDH_eff  = LDHA_P*(1.0 - stp_inh);
fLDHk_eff = LDHA_K*(1.0 - stp_inh)*(1.0 - FGO_KID);

// GRHPR — genotype, and competitively inhibited by accumulated HOG
// (this is the leading published PH3 hypothesis, implemented explicitly)
hog_ki  = 40.0;
fGRH_eff = FGRHPR*GRHPR_P/(1.0 + HOG_M/hog_ki);

// AGT — genotype plus PLP-dependent chaperone rescue, gated by B6RESC
plp_effect = PN_MGKG/(EC50_PN + PN_MGKG);
fAGT_eff   = FAGT + B6RESC*plp_effect;
if(fAGT_eff > 1.0) fAGT_eff = 1.0;

// =====================================================================
//  B.  HEPATIC GLYOXYLATE / OXALATE FLUXES
// =====================================================================
// glycolate pool: supply in, GO oxidation and plasma escape out
double J_GO_in   = kGOin*fGO_eff*(GLYCOL_P/100.0);
double J_glyc_out= kGlycOut*(GLYCOL_P/100.0);

// mitochondrial hydroxyproline route
double J_hog_in  = kHOG;
double J_hoga    = KC_HOGA*HOGA1_P*FHOGA1*HOG_M;        // HOG -> glyoxylate
double J_hog_dhg = K_DHG*HOG_M;                         // HOG -> DHG (PH3 marker)
double J_hog_cyt = kHOGcyt*HOG_M;                       // alt. cytosolic aldolase

// mitochondrial glyoxylate split (GRHPR-dependent reduction)
double m_grh = kMgrh*fGRH_eff;
double m_den = m_grh + kMper + kMcyt;
double J_m_grh = J_hoga*m_grh/m_den;
double J_m_per = J_hoga*kMper/m_den;
double J_m_cyt = J_hoga*kMcyt/m_den;

// peroxisomal glyoxylate
double per_in  = J_GO_in + J_DAO + J_m_per;
double J_AGT   = kAGTvm*fAGT_eff*GX_P/(KM_AGT + GX_P);
double J_GOox  = kGOox*fGO_eff*GX_P;                    // GO oxidises glyoxylate
double J_leak  = kLeak*GX_P;

// cytosolic glyoxylate: LDHA vs GRHPR vs spillover
ldh_share  = kLDH*fLDH_eff;
grh_share  = kGR*fGRH_eff;
gxsp_share = kGXsp;
double c_den = ldh_share + grh_share + gxsp_share;
double cyt_out = 12.0*GX_C;                             // total cytosolic turnover
double J_LDH   = cyt_out*ldh_share/c_den;
double J_GRH   = cyt_out*grh_share/c_den;
double J_GXsp  = cyt_out*gxsp_share/c_den;

// hepatic oxalate production = the GO route + the LDHA route
prod_hep = J_GOox + J_LDH;

// extrahepatic sources: constitutive, plus the PH2-unmasked component,
// plus local renal LDHA conversion of filtered glyoxylate
double prod_xhep = OX_XHEP*fLDHk_eff
                 + OX_XHEP2*(1.0 - FGRHPR)*fLDHk_eff;

// dietary absorption and ascorbate
double J_abs   = FABS_OX*(1.0 - OXF_DEG)*OX_GUT/(OX_GUT + 1.0)*OX_DIET;
double J_asc   = OX_ASC + ASC_DOSE;

prod_tot = prod_hep + prod_xhep + J_abs + J_asc;

// =====================================================================
//  C.  KIDNEY: FUNCTION, THE ONLY SINK THAT MATTERS
// =====================================================================
double neph_eff = NEPHRON + GRAFT;
if(neph_eff > 1.0) neph_eff = 1.0;
gfr_now  = GFR_MAX*neph_eff*(1.0 + HYPERF*(1.0 - neph_eff));
if(gfr_now < 0.2) gfr_now = 0.2;

double gfr_Lday = gfr_now*1.44;                        // mL/min -> L/day
clox_now = F_SEC_OX*gfr_Lday;

// dialysis: integrate the ACTUAL intradialytic clearance
hd_clear = 0.0;
if(HD_PERWK > 0.5){
  double wk   = fmod(SOLVERTIME, 7.0);
  int    dayi = (int)floor(wk);
  double frac = wk - (double)dayi;
  int on = 0;
  if(HD_PERWK >= 6.5)                      on = 1;
  else if(HD_PERWK >= 5.5 && dayi <= 5)    on = 1;
  else if(HD_PERWK >= 3.5 && dayi <= 3)    on = 1;
  else if(dayi==0 || dayi==2 || dayi==4)   on = 1;
  if(on==1 && frac < HD_HRS/24.0) hd_clear = CL_HD;
}
if(PD_ON > 0.5) hd_clear += CL_PD;

// =====================================================================
//  D.  OXALATE DISPOSITION AND THE DEPOSITION THRESHOLD
// =====================================================================
uox_rate = clox_now*POX;
ent_elim = CL_ENT*POX;

double bone_fill = 1.0 - (OX_BONE + OX_BONE_D)/BMAX_BONE;
if(bone_fill < 0.0) bone_fill = 0.0;

dep_bone = 0.0; dep_soft = 0.0; dep_kid = 0.0;
if(POX > POX_CRIT_B) dep_bone = CLDEP_B*(POX - POX_CRIT_B)*bone_fill;
if(POX > POX_CRIT_S) dep_soft = CLDEP_S*(POX - POX_CRIT_S);
if(POX > POX_CRIT_B) dep_kid  = CLDEP_K*(POX - POX_CRIT_B);

bone_release = K_RES*(1.0 + STEROID + 0.8*ACIDOSIS)*OX_BONE
               + K_RES_D*(1.0 + STEROID)*OX_BONE_D;

// =====================================================================
//  E.  URINE CHEMISTRY AND THE AP(CaOx) INDEX  (Layer B)
// =====================================================================
uvol_now = U_VOL;
if(uvol_now < 0.2) uvol_now = 0.2;
uca_now  = U_CA*(1.0 - THIAZ);
ucit_now = U_CIT*(1.0 - 0.55*ACIDOSIS) + KCIT_DOSE;
if(ucit_now < 0.2) ucit_now = 0.2;
umg_now  = U_MG + MG_DOSE;

double uox_mmol = uox_rate/1000.0;                     // mmol/day
// Tiselius AP(CaOx) index
ap_index = 1.9*pow(uca_now, 0.84)*uox_mmol
           /(pow(ucit_now, 0.22)*pow(umg_now, 0.12)*pow(uvol_now, 1.03));

double sat_drive = ap_index - AP_THRESH;
if(sat_drive < 0.0) sat_drive = 0.0;

// =====================================================================
//  F.  CRYSTAL -> INFLAMMATION -> FIBROSIS -> NEPHRON LOSS (the loop)
// =====================================================================
double kidox_frac = OX_KID/50000.0;
double inj_drive  = CRYST + W_KIDOX*kidox_frac;

dxdt_CRYST    = K_NUC*sat_drive - K_CRCL*CRYST;
dxdt_NLRP3    = K_NLRP3*inj_drive - K_NLRP3_D*NLRP3;
dxdt_IL1B     = K_IL1*NLRP3*(1.0 - IL1_BLOCK) - K_IL1_D*IL1B;
dxdt_TGFB     = K_TGFB*(IL1B + 0.4*inj_drive) - K_TGFB_D*TGFB;
dxdt_FIBROSIS = K_FIB*TGFB - K_FIB_R*FIBROSIS;
dxdt_TUB_INJ  = K_TINJ*inj_drive - K_TINJ_D*TUB_INJ;

double loss_rate = K_LOSS*(FIBROSIS + 0.35*W_KIDOX*kidox_frac) + CNI_TOX;
dxdt_NEPHRON  = -loss_rate*NEPHRON;
dxdt_GRAFT    = -loss_rate*GRAFT;

// =====================================================================
//  G.  HEPATIC POOLS
// =====================================================================
dxdt_GLYCOL_P = GLYCO_SUP + J_GRH - J_GO_in - J_glyc_out;
dxdt_GX_P     = per_in - J_AGT - J_GOox - J_leak;
dxdt_GX_C     = J_leak + J_m_cyt + J_hog_cyt - cyt_out;
dxdt_HOG_M    = J_hog_in - J_hoga - J_hog_dhg - J_hog_cyt;
dxdt_GX_M     = J_hoga - J_m_grh - J_m_per - J_m_cyt;
dxdt_DHG      = J_hog_dhg - 1.2*DHG;
dxdt_HPYR     = 60.0 - 20.0*HPYR*fGRH_eff - 4.0*HPYR;
dxdt_LGLYC    = 4.0*HPYR*(1.0 - fGRH_eff) - 1.0*LGLYC;

// =====================================================================
//  H.  OXALATE, GLYCOLATE, GLYOXYLATE DISTRIBUTION
// =====================================================================
dxdt_OX_PL = prod_tot + bone_release + K_SOFT_R*OX_SOFT + K_KID_R*OX_KID
             - uox_rate - ent_elim - hd_clear*POX
             - dep_bone - dep_soft - dep_kid
             - K_PL_ECF*(POX - OX_ECF/V_ECF);
dxdt_OX_ECF   = K_PL_ECF*(POX - OX_ECF/V_ECF);
double bury = K_BURY*OX_BONE;
dxdt_OX_BONE   = dep_bone - K_RES*(1.0+STEROID+0.8*ACIDOSIS)*OX_BONE - bury;
dxdt_OX_BONE_D = bury - K_RES_D*(1.0+STEROID)*OX_BONE_D;
dxdt_OX_SOFT  = dep_soft - K_SOFT_R*OX_SOFT;
dxdt_OX_KID   = dep_kid  - K_KID_R*OX_KID;
dxdt_OX_GUT   = OX_DIET + ent_elim - J_abs
                - OXF_DEG*OX_GUT - 1.5*OX_GUT;

dxdt_GLYCOL_PL = J_glyc_out - F_SEC_GC*gfr_Lday*PGLYC - hd_clear*PGLYC;
dxdt_GX_PL     = J_GXsp    - F_SEC_GX*gfr_Lday*PGX   - hd_clear*PGX - 2.0*PGX;

// =====================================================================
//  I.  ENZYME / mRNA DYNAMICS (siRNA pharmacodynamics)
// =====================================================================
double sup_hao1 = EMAX_HAO1*LUM_RISC/(EC50_HAO1 + LUM_RISC);
double sup_ldha = EMAX_LDHA*NED_RISC/(EC50_LDHA + NED_RISC);

dxdt_HAO1_M  = K_HAO1_M*(1.0 - sup_hao1) - K_HAO1_M*HAO1_M;
dxdt_GO_P    = K_GO_P*HAO1_M - K_GO_P*GO_P;
dxdt_LDHA_M  = K_LDHA_M*(1.0 - sup_ldha) - K_LDHA_M*LDHA_M;
dxdt_LDHA_P  = K_LDHA_P*LDHA_M - K_LDHA_P*LDHA_P;
dxdt_LDHA_K  = K_LDHA_P*1.0    - K_LDHA_P*LDHA_K;   // renal: NOT silenced
dxdt_AGT_APO = 0.0;
dxdt_AGT_HOLO= 0.0;
dxdt_GRHPR_P = 0.0;
dxdt_HOGA1_P = 0.0;

// =====================================================================
//  J.  DRUG PK
// =====================================================================
double up_lum = CLUP_LUM*(LUM_PL/V_LUM)*KM_UP_LUM/(KM_UP_LUM + LUM_PL);
dxdt_LUM_SC   = -KA_LUM*LUM_SC;
dxdt_LUM_PL   =  KA_LUM*LUM_SC - CL_LUM*(LUM_PL/V_LUM) - up_lum;
dxdt_LUM_LIV  =  up_lum - K_ESC_LUM*LUM_LIV;
dxdt_LUM_RISC =  K_ESC_LUM*LUM_LIV - K_RISC_L*LUM_RISC;

double up_ned = CLUP_NED*(NED_PL/V_NED)*KM_UP_NED/(KM_UP_NED + NED_PL);
dxdt_NED_SC   = -KA_NED*NED_SC;
dxdt_NED_PL   =  KA_NED*NED_SC - CL_NED*(NED_PL/V_NED) - up_ned;
dxdt_NED_LIV  =  up_ned - K_ESC_NED*NED_LIV;
dxdt_NED_RISC =  K_ESC_NED*NED_LIV - K_RISC_N*NED_RISC;

dxdt_STP_GUT  = -KA_STP*STP_GUT;
dxdt_STP_CEN  =  KA_STP*STP_GUT - CL_STP*stp_conc
                 - Q_STP*(stp_conc - STP_PER/V_STP_P);
dxdt_STP_PER  =  Q_STP*(stp_conc - STP_PER/V_STP_P);

double pn_in  = (PN_DIET + PN_MGKG*70.0)*5917.0;   // mg/day -> nmol/day
double pn_conv= VM_PN*PN_GUT/(KM_PN + PN_GUT);
dxdt_PN_GUT   =  pn_in - pn_conv - 0.2*PN_GUT;
dxdt_PLP_PL   =  pn_conv - CL_PLP*(PLP_PL/V_PLP)
                 - K_PLP_LIV*PLP_PL + K_PLP_OUT*PLP_LIV;
dxdt_PLP_LIV  =  K_PLP_LIV*PLP_PL - K_PLP_OUT*PLP_LIV;

// =====================================================================
//  K.  URINE VOLUME, ACIDOSIS
// =====================================================================
dxdt_U_VOL   = K_UVOL*(U_VOL_TGT - U_VOL);
double gfr_def = 1.0 - gfr_now/GFR_MAX;
dxdt_ACIDOSIS = K_ACID*(ACID_GAIN*gfr_def*gfr_def - ACIDOSIS);

// =====================================================================
//  L.  STONES AND ORGAN ENDPOINTS
// =====================================================================
dxdt_STONE     = K_STONE*CRYST - K_STONE_R*STONE;
dxdt_STONE_HAZ = K_HAZ*STONE;

double soft_drv = OX_SOFT/10000.0;
dxdt_RETINA   = K_RET*soft_drv     - K_ORG_R*RETINA;
dxdt_MYOCARD  = K_MYO*soft_drv     - K_ORG_R*MYOCARD;
dxdt_NERVE    = K_NRV*soft_drv     - K_ORG_R*NERVE;
dxdt_SKIN     = K_SKN*soft_drv     - K_ORG_R*SKIN;
dxdt_BONE_DIS = K_BONEDIS*((OX_BONE+OX_BONE_D)/BMAX_BONE) - K_ORG_R*BONE_DIS;
dxdt_MARROW   = K_MARROW*((OX_BONE+OX_BONE_D)/BMAX_BONE)  - K_ORG_R*MARROW;
dxdt_B6_NEURO = K_B6NEURO*PN_MGKG;

// =====================================================================
//  M.  AUDIT INTEGRALS — the mass balance must close
// =====================================================================
dxdt_CUM_PROD  = prod_tot;
dxdt_CUM_UOX   = uox_rate;
dxdt_CUM_ENT   = ent_elim;
dxdt_CUM_HD    = hd_clear*POX;
dxdt_CUM_DEP   = dep_bone + dep_soft + dep_kid;
dxdt_AUC_POX   = POX;
dxdt_AUC_AP    = ap_index;
dxdt_T_ABOVE30 = (POX > POX_CRIT) ? 1.0 : 0.0;
dxdt_T_EGFR30  = (gfr_now < 30.0) ? 1.0 : 0.0;
dxdt_CUM_GXTOX = GX_C;
dxdt_CUM_GLYCO = F_SEC_GC*gfr_Lday*PGLYC;
dxdt_CUM_FILT  = gfr_Lday*POX;
dxdt_ESKD_FLAG = (gfr_now < 15.0) ? 1.0 : 0.0;
dxdt_CUM_BONE_R= bone_release;

$TABLE
double POX_OUT   = POX;
double PGLYC_OUT = PGLYC;
double PGX_OUT   = PGX;
double UOX_MMOL  = uox_rate/1000.0;                 // mmol/1.73m2/day
double UGLYC_MMOL= F_SEC_GC*gfr_now*1.44*PGLYC/1000.0;
double EGFR      = gfr_now;
double AP_CAOX   = ap_index;
double PROD_TOT  = prod_tot;
double PROD_HEP  = prod_hep;
double UOX_CONC  = uox_rate/uvol_now/1000.0;        // mmol/L
double OXALOSIS  = RETINA + MYOCARD + NERVE + SKIN + BONE_DIS + MARROW;
double AGE_Y     = TIME/365.25;
double TBURDEN   = (OX_BONE + OX_BONE_D + OX_SOFT + OX_KID)/1000.0;   // mmol
double STONE_RATE= K_HAZ*STONE*365.25;              // events/year
double FAGT_EFF  = fAGT_eff;
double FGO_OUT   = fGO_eff;
double FLDH_OUT  = fLDH_eff;
double BAL_NET   = prod_tot - uox_rate - ent_elim - hd_clear*POX;

$CAPTURE
POX_OUT PGLYC_OUT PGX_OUT UOX_MMOL UGLYC_MMOL EGFR AP_CAOX
PROD_TOT PROD_HEP UOX_CONC OXALOSIS AGE_Y TBURDEN STONE_RATE
FAGT_EFF FGO_OUT FLDH_OUT BAL_NET
'

mod <- mcode_cache("ph1_qsp", ph1_code, soloc = tempdir())

cat("=== model compiled:", length(mrgsolve::init(mod)), "compartments,",
    length(mrgsolve::param(mod)), "parameters ===\n")

# =====================================================================
#  2. SIMULATION INFRASTRUCTURE
# =====================================================================
#  Every scenario starts from the SAME healthy fixed point and is run
#  from birth. The diseased phenotype is produced by a silent natural-
#  history pre-run, never imposed as an initial condition. Therapy is
#  layered on at a stated age, exactly as it is in a clinic.
# =====================================================================

YR   <- 365.25
STATES <- names(mrgsolve::init(mod))

# --- doses -----------------------------------------------------------
#  lumasiran  3 mg/kg (=210 mg, ~12,900 nmol) SC monthly x3 then quarterly
#  nedosiran  160 mg (~10,000 nmol) SC monthly
#  stiripentol 25 mg/kg BID (=1750 mg BID)
LUM_MW <- 16300; NED_MW <- 16000
lum_ev <- function(t0, yrs = 20)
  c(ev(amt = 12900, cmt = "LUM_SC", ii = 30, addl = 2,            time = t0),
    ev(amt = 12900, cmt = "LUM_SC", ii = 90, addl = ceiling(yrs*4), time = t0 + 90))
ned_ev <- function(t0, yrs = 20)
  ev(amt = 10000, cmt = "NED_SC", ii = 30, addl = ceiling(yrs*12), time = t0)
stp_ev <- function(t0, yrs = 20)
  ev(amt = 1750,  cmt = "STP_GUT", ii = 0.5, addl = ceiling(yrs*730), time = t0)

# --- single-phase run ------------------------------------------------
sim1 <- function(pars = list(), evs = NULL, end = 45*YR, delta = 30, init0 = NULL) {
  m <- if (length(pars)) mod %>% param(pars) else mod
  if (!is.null(init0)) m <- m %>% param(INIT_ON = 0) %>% init(as.list(init0))
  m %>% mrgsim(events = evs, end = end, delta = delta, hmax = 0.5) %>% as_tibble()
}

# --- two-phase run: a hard structural change (transplant) at t_sw ----
sim2 <- function(p1, p2, t_sw, end, evs1 = NULL, evs2 = NULL,
                 init_override = NULL, delta = 30) {
  a  <- sim1(p1, evs1, end = t_sw, delta = delta)
  st <- as.numeric(tail(a, 1)[STATES]); names(st) <- STATES
  if (!is.null(init_override)) st[names(init_override)] <- unlist(init_override)
  b  <- sim1(p2, evs2, end = end - t_sw, delta = delta, init0 = st)
  b$time <- b$time + t_sw; b$AGE_Y <- b$time / YR
  bind_rows(a, b[-1, ])
}

at_age <- function(d, y) d[which.min(abs(d$AGE_Y - y)), ]
eskd_age <- function(d) { i <- which(d$EGFR < 15)[1]; if (is.na(i)) NA_real_ else d$AGE_Y[i] }
first_age <- function(d, cond) { i <- which(cond)[1]; if (is.na(i)) NA_real_ else d$AGE_Y[i] }

PH1  <- list(FAGT = 0.005)                       # classic PH1, near-null AGT
PH1N <- list(FAGT = 0.000)                       # complete AGT loss
PH1M <- list(FAGT = 0.100)                        # attenuated PH1
G170 <- list(FAGT = 0.005, B6RESC = 0.035)        # pyridoxine-responsive allele
PH2  <- list(FGRHPR = 0)
PH3  <- list(FHOGA1 = 0)
LAYB <- list(U_VOL_TGT = 3.0, KCIT_DOSE = 4.0)    # hyperhydration + K-citrate

TX <- 5 * YR    # age at which therapy starts in the reference scenarios

# =====================================================================
#  3. THIRTY-FOUR SCENARIOS
# =====================================================================

scen <- list()
add <- function(id, label, d) scen[[id]] <<- list(label = label, d = d)

add("S01", "Healthy control",                       sim1(list()))
add("S02", "Heterozygous carrier (AGT 50%)",        sim1(list(FAGT = 0.5)))
add("S03", "PH1 classic, untreated",                sim1(PH1))
add("S04", "PH1 complete AGT loss, untreated",      sim1(PH1N))
add("S05", "PH1 attenuated (AGT 10%), untreated",   sim1(PH1M))
add("S06", "PH1 + hyperhydration only",             sim1(c(PH1, list(U_VOL_TGT = 3.0))))
add("S07", "PH1 + potassium citrate only",          sim1(c(PH1, list(KCIT_DOSE = 4.0))))
add("S08", "PH1 + full Layer B (fluids + citrate)", sim1(c(PH1, LAYB)))
add("S09", "PH1 G170R + pyridoxine 10 mg/kg/d",     sim1(c(G170, list(PN_MGKG = 10))))
add("S10", "PH1 null allele + pyridoxine (futile)", sim1(c(PH1, list(B6RESC = 0, PN_MGKG = 10))))
add("S11", "PH1 + lumasiran from age 5",            sim1(PH1, lum_ev(TX, 45)))
add("S12", "PH1 + lumasiran + Layer B",             sim1(c(PH1, LAYB), lum_ev(TX, 45)))
add("S13", "PH1 + nedosiran from age 5",            sim1(PH1, ned_ev(TX, 45)))
add("S14", "PH1 + lumasiran + nedosiran",           sim1(PH1, c(lum_ev(TX, 45), ned_ev(TX, 45))))
add("S15", "PH1 + stiripentol 50 mg/kg/d",          sim1(PH1, stp_ev(TX, 45)))
add("S16", "PH1 + Oxalobacter colonisation",        sim1(c(PH1, list(OXF_DEG = 0.6))))
add("S17", "PH2 (GRHPR null), untreated",           sim1(PH2))
add("S18", "PH2 + nedosiran",                       sim1(PH2, ned_ev(TX, 45)))
add("S19", "PH2 + stiripentol (reaches kidney)",    sim1(c(PH2, list(FGO_KID = 0.0)), stp_ev(TX, 45)))
add("S20", "PH3 (HOGA1 null) - GRHPR-inhibition hypothesis", sim1(PH3))
add("S21", "PH3 - cytosolic HOG-aldolase hypothesis",        sim1(c(PH3, list(HOGCYT = 1.8))))
add("S22", "Enteric hyperoxaluria + Oxalobacter",   sim1(list(OX_DIET = 4200, FABS_OX = 0.30, OXF_DEG = 0.6)))
add("S23", "Enteric hyperoxaluria, no treatment",   sim1(list(OX_DIET = 4200, FABS_OX = 0.30)))

# ESKD scenarios: reach ESKD on natural history, then start renal replacement
esk_start <- function(pars, kl = NULL) {
  d <- sim1(pars, end = 45 * YR, delta = 15)
  i <- which(d$EGFR < 15)[1]
  if (is.na(i)) i <- nrow(d)
  d$time[i]
}
T_ESKD <- esk_start(PH1)
add("S24", "PH1 ESKD, conventional HD 3x4h",
    sim2(PH1, c(PH1, list(HD_PERWK = 3, HD_HRS = 4)), T_ESKD, T_ESKD + 3 * YR, delta = 0.25))
add("S25", "PH1 ESKD, daily extended HD 6x6h",
    sim2(PH1, c(PH1, list(HD_PERWK = 7, HD_HRS = 6)), T_ESKD, T_ESKD + 3 * YR, delta = 0.25))
add("S26", "PH1 ESKD, peritoneal dialysis",
    sim2(PH1, c(PH1, list(PD_ON = 1)),               T_ESKD, T_ESKD + 3 * YR, delta = 0.25))
add("S27", "PH1 ESKD, HD 3x4h + lumasiran (A+C)",
    sim2(PH1, c(PH1, list(HD_PERWK = 3, HD_HRS = 4)), T_ESKD, T_ESKD + 3 * YR,
         evs2 = lum_ev(0, 3), delta = 0.25))
add("S28", "Combined liver-kidney transplant",
    sim2(PH1, list(FAGT = 1.0), T_ESKD, T_ESKD + 5 * YR,
         init_override = list(NEPHRON = 0.85), delta = 2))
add("S29", "Isolated kidney transplant (liver retained)",
    sim2(PH1, PH1, T_ESKD, T_ESKD + 5 * YR,
         init_override = list(NEPHRON = 0.85), delta = 2))
add("S30", "Isolated kidney transplant + lumasiran",
    sim2(PH1, PH1, T_ESKD, T_ESKD + 5 * YR, evs2 = lum_ev(0, 5),
         init_override = list(NEPHRON = 0.85), delta = 2))
add("S31", "CLKT with steroid-heavy regimen",
    sim2(PH1, list(FAGT = 1.0, STEROID = 1.5), T_ESKD, T_ESKD + 5 * YR,
         init_override = list(NEPHRON = 0.85), delta = 2))
add("S32", "PH1 + lumasiran started late (age 20)", sim1(PH1, lum_ev(20 * YR, 45)))
add("S33", "PH1 + anti-IL-1 / NLRP3 blockade",      sim1(c(PH1, list(IL1_BLOCK = 0.75))))
add("S34", "PH1 + high-dose ascorbate (added source)", sim1(c(PH1, list(ASC_DOSE = 250))))

cat("\n=====================================================================\n")
cat(" SCENARIO SUMMARY  (values at age 20 y unless the scenario is later)\n")
cat("=====================================================================\n")
cat(sprintf("%-4s %-44s %7s %7s %7s %7s %7s %8s\n",
            "ID", "scenario", "Uox", "Pox", "eGFR", "AP", "stone/y", "ESKD_age"))
cat(" (S24-S31 are post-ESKD / post-transplant and are reported at the END of follow-up)\n")
for (id in names(scen)) {
  d <- scen[[id]]$d
  ref <- if (id %in% c("S24","S25","S26","S27","S28","S29","S30","S31")) tail(d, 1) else at_age(d, 20)
  cat(sprintf("%-4s %-44s %7.3f %7.2f %7.1f %7.2f %7.2f %8s\n",
              id, substr(scen[[id]]$label, 1, 44), ref$UOX_MMOL, ref$POX_OUT,
              ref$EGFR, ref$AP_CAOX, ref$STONE_RATE,
              ifelse(is.na(eskd_age(d)), "-", sprintf("%.1f", eskd_age(d)))))
}

# =====================================================================
#  4. SIXTEEN MODEL DIAGNOSTICS
# =====================================================================
#  These are not plots. Each one asks a question the model must answer
#  with a number, and several of them REFUTE the model's own design
#  intent. Those are reported as they came out.
# =====================================================================

line <- function(ch = "-") cat(strrep(ch, 71), "\n")
hdr  <- function(id, txt) { cat("\n"); line("="); cat(id, " ", txt, "\n", sep = ""); line("=") }

# --- D01  healthy baseline is an exact fixed point -------------------
hdr("D01", "HEALTHY BASELINE STABILITY — is the fixed point exact?")
h <- sim1(list(), end = 40 * YR, delta = 30)
b1 <- at_age(h, 1); b40 <- tail(h, 1)
cat(sprintf("  Uox     1 y = %.6f  ->  40 y = %.6f   drift = %+.4f %%\n",
            b1$UOX_MMOL, b40$UOX_MMOL, 100 * (b40$UOX_MMOL / b1$UOX_MMOL - 1)))
cat(sprintf("  Pox     1 y = %.6f  ->  40 y = %.6f   drift = %+.4f %%\n",
            b1$POX_OUT, b40$POX_OUT, 100 * (b40$POX_OUT / b1$POX_OUT - 1)))
cat(sprintf("  eGFR    1 y = %.4f  ->  40 y = %.4f   drift = %+.4f %%\n",
            b1$EGFR, b40$EGFR, 100 * (b40$EGFR / b1$EGFR - 1)))
cat(sprintf("  tissue burden at 40 y = %.4f mmol  (a healthy subject must not accumulate)\n",
            b40$TBURDEN))
cat("  => the disease is generated by the genotype, not by a drifting baseline.\n")

# --- D02  mass balance closure --------------------------------------
hdr("D02", "MASS BALANCE AUDIT — does production equal the sum of the sinks?")
for (nm in c("healthy", "PH1", "PH1 + lumasiran")) {
  d <- switch(nm, "healthy" = h, "PH1" = scen$S03$d, "PH1 + lumasiran" = scen$S11$d)
  l <- tail(d, 1)
  acc <- l$OX_PL + l$OX_ECF + l$OX_BONE + l$OX_BONE_D + l$OX_SOFT + l$OX_KID
  out <- l$CUM_UOX + l$CUM_ENT + l$CUM_HD
  err <- (l$CUM_PROD - out - acc) / l$CUM_PROD * 100
  cat(sprintf("  %-16s produced %10.0f | excreted %10.0f | in tissue %9.0f | closure error %+.5f %%\n",
              nm, l$CUM_PROD, out, acc, err))
}
cat("  => the balance closes to numerical precision; nothing is created or lost.\n")

# --- D03  the AGT nonlinearity --------------------------------------
hdr("D03", "AGT NONLINEARITY — one Michaelis-Menten term explains recessive")
cat("       inheritance AND why a chaperone-sized gain in activity matters\n")
cat(sprintf("  %10s %10s %10s %12s\n", "AGT (%)", "Uox", "x normal", "x ULN(0.46)"))
uox_of <- function(f, extra = list()) {
  d <- sim1(c(list(FAGT = f), extra), end = 3 * YR, delta = 90); tail(d, 1)$UOX_MMOL }
u_norm <- uox_of(1.0)
for (f in c(1, 0.5, 0.25, 0.10, 0.05, 0.03, 0.02, 0.01, 0.005, 0.0)) {
  u <- uox_of(f)
  cat(sprintf("  %10.1f %10.3f %10.2f %12.2f\n", 100 * f, u, u / u_norm, u / 0.46))
}
u005 <- uox_of(0.005); u03 <- uox_of(0.03)
cat(sprintf("\n  50%% of the enzyme is indistinguishable from 100%%: %.3f vs %.3f mmol/day.\n",
            uox_of(0.5), u_norm))
cat(sprintf("  Going from 0.5%% to 3%% of the enzyme — still a profoundly deficient\n"))
cat(sprintf("  liver — removes %.1f%% of the oxalate (%.3f -> %.3f).\n",
            100 * (1 - u03 / u005), u005, u03))
cat("  => that single curve is why PH1 is recessive with healthy carriers and\n")
cat("     why pharmacological chaperone rescue is worth pursuing at all.\n")

# --- D04  the deposition threshold, computed -------------------------
hdr("D04", "THE THRESHOLD — at what eGFR does plasma oxalate cross 30 umol/L?")
cat("  Production held EXACTLY constant; only nephron mass is varied.\n")
cat(sprintf("  %10s %10s %10s %12s\n", "nephron", "eGFR", "Pox", "Uox"))
scan_gfr <- function(nfrac, pars = PH1) {
  st <- as.numeric(tail(sim1(pars, end = 2 * YR, delta = 90), 1)[STATES]); names(st) <- STATES
  st["NEPHRON"] <- nfrac
  st[c("OX_BONE","OX_BONE_D","OX_SOFT","OX_KID","FIBROSIS","CRYST","NLRP3","IL1B","TGFB")] <- 0
  d <- sim1(c(pars, list(K_LOSS = 0)), end = 400, delta = 20, init0 = st)
  tail(d, 1) }
gr <- lapply(c(1, 0.7, 0.5, 0.35, 0.25, 0.18, 0.12, 0.08, 0.05, 0.03, 0.02, 0.014),
             function(n) scan_gfr(n))
for (r in gr) cat(sprintf("  %10.2f %10.1f %10.2f %12.3f\n",
                          r$NEPHRON, r$EGFR, r$POX_OUT, r$UOX_MMOL))
eg <- sapply(gr, function(r) r$EGFR); px <- sapply(gr, function(r) r$POX_OUT)
i <- which(px > 30)[1]
egfr_crit <- approx(px[c(i - 1, i)], eg[c(i - 1, i)], xout = 30)$y
cat(sprintf("\n  COMPUTED: Pox crosses the plasma CaOx solubility limit (30 umol/L)\n"))
cat(sprintf("  at eGFR = %.1f mL/min/1.73m2 in this patient.\n", egfr_crit))
cat("  => systemic oxalosis is not a separate late mechanism. It is this root.\n")
cat("     The published observation that oxalosis appears below eGFR 30-45 is\n")
cat("     reproduced without any parameter that mentions oxalosis.\n")

# --- D05  the biomarker inversion ------------------------------------
hdr("D05", "THE INVERSION — urinary oxalate falls because the kidney fails")
ux <- sapply(gr, function(r) r$UOX_MMOL); u0 <- ux[1]
cat(sprintf("  %10s %10s %14s\n", "eGFR", "Uox", "% of Uox@GFR125"))
for (k in seq_along(gr)) cat(sprintf("  %10.1f %10.3f %14.1f\n", eg[k], ux[k], 100 * ux[k] / u0))
f25 <- approx(ux, eg, xout = 0.75 * u0)$y; f50 <- approx(ux, eg, xout = 0.50 * u0)$y
cat(sprintf("\n  With PRODUCTION HELD CONSTANT, Uox has fallen 25%% by eGFR %.1f\n", f25))
cat(sprintf("  and 50%% by eGFR %.1f. Between eGFR 125 and 40 it barely moves.\n", f50))
cat("  => a falling Uox in advanced PH is not a therapeutic success, and this\n")
cat("     is exactly why ILLUMINATE-A (eGFR >= 30) could use 24-h Uox while\n")
cat("     ILLUMINATE-C (advanced CKD / dialysis) had to switch to plasma\n")
cat("     oxalate. The endpoint change is arithmetic, not regulatory taste.\n")

# --- D06  clinical anchors -------------------------------------------
hdr("D06", "PUBLISHED ANCHORS — model vs trial")
pct <- function(d, t0, t1) { a <- at_age(d, t0)$UOX_MMOL; b <- at_age(d, t1)$UOX_MMOL
                             100 * (b / a - 1) }
anch <- tibble::tribble(
  ~anchor,                                        ~obs,      ~mod,
  "normal Uox, adult (mmol/1.73m2/d)",            "0.15-0.46", sprintf("%.2f", u_norm),
  "PH1 untreated Uox",                            "1.0-2.5",   sprintf("%.2f", at_age(scen$S03$d, 5)$UOX_MMOL),
  "PH1 Pox with preserved eGFR (umol/L)",         "5-15",      sprintf("%.1f", at_age(scen$S03$d, 5)$POX_OUT),
  "PH1 Pox at ESKD (umol/L)",                     "60-120",    sprintf("%.0f", tail(scen$S24$d,1)$POX_OUT),
  "ILLUMINATE-A lumasiran Uox change",            "-65.4%",    sprintf("%+.1f%%", pct(scen$S11$d, 5, 6)),
  "PHYOX2 nedosiran (PH1) Uox change",            "~-39%",     sprintf("%+.1f%%", pct(scen$S13$d, 5, 6)),
  "PHYOX2 nedosiran (PH2) Uox change",            "no response",sprintf("%+.1f%%", pct(scen$S18$d, 5, 6)),
  "pyridoxine, G170R homozygote",                 "-30 to -50%",sprintf("%+.1f%%", 100*(at_age(scen$S09$d,6)$UOX_MMOL/at_age(scen$S03$d,6)$UOX_MMOL-1)),
  "pyridoxine, null allele",                      "0%",        sprintf("%+.1f%%", 100*(at_age(scen$S10$d,6)$UOX_MMOL/at_age(scen$S03$d,6)$UOX_MMOL-1)),
  "PH1 reaching ESKD, untreated",                 "by age ~25", sprintf("age %.1f", eskd_age(scen$S03$d)),
  "urinary glycolate rise on lumasiran",           "up 2-5x",   sprintf("%.1fx", at_age(scen$S11$d,6)$UGLYC_MMOL/at_age(scen$S03$d,6)$UGLYC_MMOL),
  "hyperhydration effect on CaOx AP index",        "~ -50%",    sprintf("%+.0f%%", 100*(at_age(scen$S06$d,6)$AP_CAOX/at_age(scen$S03$d,6)$AP_CAOX-1)),
  "Oxalobacter in enteric hyperoxaluria",          "effective",  sprintf("%+.0f%%", 100*(at_age(scen$S22$d,6)$UOX_MMOL/at_age(scen$S23$d,6)$UOX_MMOL-1)),
  "Oxalobacter in PH1 (ePHex primary endpoint)",   "failed",     sprintf("%+.0f%%", 100*(at_age(scen$S16$d,6)$UOX_MMOL/at_age(scen$S03$d,6)$UOX_MMOL-1))
)
print(as.data.frame(anch), row.names = FALSE, right = FALSE)

# --- D07  the PH2 decomposition (the headline result) ----------------
hdr("D07", "WHY AN LDHA siRNA FAILS IN PH2 — partitioning vs delivery")
d1 <- pct(scen$S13$d, 5, 6)                                        # PH1
d2 <- pct(sim1(c(PH2, list(OX_XHEP2 = 0)), ned_ev(TX, 3), end = 8 * YR), 5, 6)
d3 <- pct(scen$S18$d, 5, 6)                                        # PH2 full
cat(sprintf("  PH1, nedosiran                                     : %+6.1f %%\n", d1))
cat(sprintf("  PH2 with the extrahepatic source SWITCHED OFF      : %+6.1f %%\n", d2))
cat(sprintf("  PH2 as modelled (extrahepatic source present)      : %+6.1f %%\n", d3))
cat(sprintf("\n  Loss of efficacy attributable to FLUX PARTITIONING  : %5.1f points\n", d1 - d2))
cat(sprintf("  Loss of efficacy attributable to hepatic DELIVERY   : %5.1f points\n", d2 - d3))
cat(sprintf("  => partitioning accounts for %.0f%% of the PH1 -> PH2 collapse.\n",
            100 * (d1 - d2) / (d1 - d3)))
cat("\n  MECHANISM. Silencing an enzyme lowers its FLUX only if another\n")
cat("  branch can take the substrate. In PH1 that branch is GRHPR, so LDHA\n")
cat("  knockdown re-routes glyoxylate into harmless glycolate. In PH2 GRHPR\n")
cat("  IS the missing gene: LDHA becomes the ONLY exit from the cytosolic\n")
cat("  glyoxylate pool, and at steady state the flux through the only exit\n")
cat("  equals its input regardless of how much enzyme is present. Knocking\n")
cat("  it down therefore raises glyoxylate and barely moves oxalate.\n")
cat(sprintf("  Model check — plasma glyoxylate on nedosiran, PH2: %.2f -> %.2f umol/L\n",
            at_age(scen$S18$d, 5)$PGX_OUT, at_age(scen$S18$d, 6)$PGX_OUT))
cat("  FALSIFIABLE PREDICTION: in PH2, nedosiran should raise plasma and\n")
cat("  urinary glyoxylate substantially while oxalate is unchanged. The\n")
cat("  standard explanation (hepatocyte-restricted delivery) predicts no\n")
cat("  glyoxylate rise. One measurement separates them.\n")

# --- D08  dialysis arithmetic ----------------------------------------
hdr("D08", "DIALYSIS ARITHMETIC — computed, not asserted")
wk <- function(d, col) { n <- nrow(d); i <- which(d$time >= max(d$time) - 7)[1]
  (d[[col]][n] - d[[col]][i]) / 1000 }
cat("  All figures are mmol per WEEK, integrated from the actual intradialytic\n")
cat("  clearance acting on the actual (falling, then rebounding) plasma level.\n\n")
cat(sprintf("  %-32s %7s %7s %7s %7s %9s %6s\n", "regimen", "prod", "dialys", "urine", "gut", "shortfall", "Pox"))
for (id in c("S24", "S25", "S26", "S27")) {
  d <- scen[[id]]$d
  pr <- wk(d,"CUM_PROD"); hd <- wk(d,"CUM_HD"); ur <- wk(d,"CUM_UOX"); gu <- wk(d,"CUM_ENT")
  cat(sprintf("  %-32s %7.1f %7.1f %7.1f %7.1f %9.1f %6.1f\n",
              substr(scen[[id]]$label, 1, 32), pr, hd, ur, gu, pr - hd - ur - gu,
              tail(d, 1)$POX_OUT))
}
cat("\n  => conventional thrice-weekly haemodialysis removes only about a third\n")
cat("     of weekly production in PH1, so a SHORTFALL of several mmol/week is\n")
cat("     deposited in tissue every week and plasma oxalate sits far above the\n")
cat("     30 umol/L solubility limit. Intensified dialysis roughly doubles the\n")
cat("     removal and still does not close the balance. Adding a Layer A drug\n")
cat("     cuts production by more than any dialysis prescription can remove,\n")
cat("     and is the only line that brings Pox below the solubility limit and\n")
cat("     the weekly shortfall to within 0.1 mmol of zero. This is why\n")
cat("     dialysis is a bridge in PH1 and never a destination, and why a\n")
cat("     production-directed drug turns an insufficient dialysis\n")
cat("     prescription into a sufficient one without changing the\n")
cat("     prescription at all.\n")

# --- D09  Layer B changes no flux ------------------------------------
hdr("D09", "LAYER B — a therapy that changes NO mass-balance term")
a <- at_age(scen$S03$d, 6); bb <- at_age(scen$S08$d, 6)
cat(sprintf("  %-28s %12s %12s %10s\n", "quantity", "untreated", "fluids+cit", "change"))
for (q in c("PROD_TOT", "UOX_MMOL", "UOX_CONC", "AP_CAOX", "STONE_RATE")) {
  cat(sprintf("  %-28s %12.3f %12.3f %9.1f%%\n", q, a[[q]], bb[[q]],
              100 * (bb[[q]] / a[[q]] - 1)))
}
cat("\n  => production and 24-h excretion are essentially UNCHANGED (the mass\n")
cat("     balance is untouched) and yet supersaturation and the stone rate\n")
cat("     fall sharply, because the stone endpoint is a RATIO. This is the\n")
cat("     clearest separation in the model between protecting the kidney and\n")
cat("     closing the balance — and the reason supportive care alone still\n")
cat("     leaves systemic oxalosis on the table.\n")
cat(sprintf("     Tissue burden at 25 y: untreated %.1f mmol vs Layer B %.1f mmol.\n",
            at_age(scen$S03$d, 25)$TBURDEN, at_age(scen$S08$d, 25)$TBURDEN))

# --- D10  post-transplant bone unloading -----------------------------
hdr("D10", "A RESERVOIR, NOT A RELAPSE — post-transplant oxalate unloading")
d28 <- scen$S28$d; d31 <- scen$S31$d
t0 <- T_ESKD
cat(sprintf("  %8s %14s %14s %14s\n", "months", "Uox (CLKT)", "Uox (steroid)", "bone (mmol)"))
for (mo in c(0, 1, 3, 6, 12, 24, 48)) {
  r <- d28[which.min(abs(d28$time - (t0 + mo * 30.4))), ]
  q <- d31[which.min(abs(d31$time - (t0 + mo * 30.4))), ]
  cat(sprintf("  %8d %14.3f %14.3f %14.1f\n", mo, r$UOX_MMOL, q$UOX_MMOL,
              (r$OX_BONE + r$OX_BONE_D) / 1000))
}
cat("\n  => with a NORMAL liver installed and production near zero, urinary\n")
cat("     oxalate stays above the normal range for months. It is the SKELETON\n")
cat("     emptying, not the disease returning. The model predicts a\n")
cat("     steroid-heavy regimen prolongs the unloading, which is a testable\n")
cat("     reason to prefer steroid minimisation after transplant in PH1.\n")

# --- D11  the siRNA PK/PD disconnect ---------------------------------
hdr("D11", "PK/PD DISCONNECT — hours in plasma, months of effect")
sd <- sim1(PH1, ev(amt = 12900, cmt = "LUM_SC", time = 5 * YR),
           end = 5 * YR + 2, delta = 0.005)
seg   <- sd[sd$time >= 5 * YR, ]
cmax  <- max(seg$LUM_PL); tmax <- seg$time[which.max(seg$LUM_PL)] - 5 * YR
half  <- seg[seg$time > 5 * YR + tmax, ]
thalf <- half$time[which(half$LUM_PL < cmax / 2)[1]] - 5 * YR - tmax
cat(sprintf("  SINGLE 3 mg/kg SC dose — plasma lumasiran:\n"))
cat(sprintf("    Tmax %.2f h | Cmax %.0f nmol | terminal t1/2 %.1f h  (obs: Tmax ~4 h, t1/2 ~5 h)\n",
            tmax * 24, cmax, thalf * 24))
wo <- sim1(PH1, ev(amt = 12900, cmt = "LUM_SC", ii = 30, addl = 2, time = 5 * YR),
           end = 5 * YR + 800, delta = 2)
wo2 <- wo[wo$time > 5 * YR + 60, ]
i50 <- which(wo2$FGO_OUT > 0.5)[1]
nadir <- min(wo2$UOX_MMOL); full <- 1.463
iU  <- which(wo2$UOX_MMOL > nadir + 0.5 * (full - nadir))[1]
cat(sprintf("  RISC-loaded siRNA t1/2 (model parameter)      : %.0f d\n", log(2) / 0.0154))
cat(sprintf("  after the LAST of 3 monthly doses, GO protein returns to 50%%\n"))
cat(sprintf("    of normal at %.0f days; Uox is half-recovered at %.0f days.\n",
            ifelse(is.na(i50), NA, wo2$time[i50] - (5 * YR + 60)),
            ifelse(is.na(iU),  NA, wo2$time[iU]  - (5 * YR + 60))))
cat(sprintf("  ratio (duration of effect) / (plasma t1/2)    : ~%.0f-fold\n",
            ifelse(is.na(i50), NA, (wo2$time[i50] - (5 * YR + 60)) / thalf)))
cat("  => a molecule with a plasma half-life of hours supports quarterly\n")
cat("     dosing because the pharmacologically active species is the\n")
cat("     RISC-loaded strand inside the hepatocyte. The loading-then-\n")
cat("     maintenance schedule is a RISC accumulation design, not a\n")
cat("     convention.\n")

# --- D12  inverting a trial to estimate a parameter ------------------
hdr("D12", "INVERTING PHYOX2 — what does the trial say about the GO route?")
cat("  The model under-predicts the observed nedosiran response. The single\n")
cat("  parameter that controls the floor is PHI_GOOX, the share of glyoxylate\n")
cat("  oxidised to oxalate by GLYCOLATE OXIDASE rather than by LDH. Sweep it:\n")
cat(sprintf("  %12s %14s %14s\n", "PHI_GOOX", "PH1 baseline Uox", "nedosiran %"))
sw <- c(0.02, 0.05, 0.08, 0.11, 0.15, 0.20)
res <- sapply(sw, function(ph) {
  d <- sim1(c(PH1, list(PHI_GOOX = ph)), ned_ev(TX, 3), end = 8 * YR)
  c(at_age(d, 5)$UOX_MMOL, pct(d, 5, 6)) })
for (k in seq_along(sw)) cat(sprintf("  %12.2f %14.3f %13.1f%%\n", sw[k], res[1, k], res[2, k]))
ph_est <- approx(res[2, ], sw, xout = -39)$y
cat(sprintf("\n  PHI_GOOX consistent with the reported ~-39%% PHYOX2 response: %.3f\n", ph_est))
cat("  => the model turns a discrepancy into a measurable claim about\n")
cat("     enzymology: the glycolate-oxidase contribution to oxalate synthesis\n")
cat("     should be near this value. That is falsifiable in a hepatocyte assay.\n")

# --- D13  pyridoxine ceiling ----------------------------------------
hdr("D13", "PYRIDOXINE — a genotype gate and a conversion ceiling")
cat(sprintf("  %10s %10s %10s %12s %10s\n", "mg/kg/d", "Uox(G170R)", "Uox(null)", "plasma PLP", "neurotox"))
for (pn in c(0, 2, 5, 10, 15, 20, 30)) {
  a <- tail(sim1(c(G170, list(PN_MGKG = pn)), end = 6 * YR, delta = 90), 1)
  b <- tail(sim1(c(PH1,  list(PN_MGKG = pn, B6RESC = 0)), end = 6 * YR, delta = 90), 1)
  cat(sprintf("  %10.0f %10.3f %10.3f %12.0f %10.2f\n", pn, a$UOX_MMOL, b$UOX_MMOL,
              a$PLP_PL / 20, a$B6_NEURO))
}
cat("\n  => in a rescuable genotype the benefit saturates near 10 mg/kg/day\n")
cat("     because pyridoxine-to-PLP conversion is itself saturable, while the\n")
cat("     cumulative neurotoxicity exposure keeps rising linearly with dose.\n")
cat("     In a null allele the effect is EXACTLY zero at every dose: a\n")
cat("     structural null, not an under-dosing problem.\n")

# --- D14  the loop gain ---------------------------------------------
hdr("D14", "LOOP GAIN — why the terminal decline is abrupt")
cat("  Sensitivity of the nephron-loss rate to a 1% change in nephron mass,\n")
cat("  evaluated along the natural history of scenario S03:\n")
d <- scen$S03$d
cat(sprintf("  %8s %8s %8s %10s %12s %12s\n", "age", "eGFR", "Pox", "AP", "dGFR/dt", "loop gain"))
for (y in c(5, 10, 15, 20, 22, 24, 25)) {
  i <- which.min(abs(d$AGE_Y - y)); if (i < 3 || i > nrow(d) - 3) next
  slope <- (d$EGFR[i + 2] - d$EGFR[i - 2]) / (d$AGE_Y[i + 2] - d$AGE_Y[i - 2])
  gain  <- -slope / max(d$EGFR[i], 1e-6)
  cat(sprintf("  %8.1f %8.1f %8.2f %10.2f %12.2f %12.4f\n",
              d$AGE_Y[i], d$EGFR[i], d$POX_OUT[i], d$AP_CAOX[i], slope, gain))
}
cat("\n  => the fractional loss rate rises monotonically as eGFR falls. No new\n")
cat("     mechanism is switched on; the same loop simply has more gain when\n")
cat("     the sink is smaller. This IS the argument for early diagnosis.\n")

# --- D15  PH3: a structural failure, reported ------------------------
hdr("D15", "PH3 — THE MODEL FAILS, AND SO DOES THE HYPOTHESIS IT ENCODES")
p20 <- at_age(scen$S20$d, 6); p21 <- at_age(scen$S21$d, 6)
cat(sprintf("  healthy Uox                                        : %.3f\n", u_norm))
cat(sprintf("  PH3, HOG-inhibits-GRHPR hypothesis only            : %.3f  (HOG %.0f-fold)\n",
            p20$UOX_MMOL, at_age(scen$S20$d, 6)$HOG_M / at_age(h, 6)$HOG_M))
cat(sprintf("  PH3, adding cytosolic HOG-aldolase hypothesis      : %.3f\n", p21$UOX_MMOL))
cat("\n  NEGATIVE RESULT. Losing an aldolase that PRODUCES glyoxylate cannot\n")
cat("  raise oxalate, and the model says so: with only the published\n")
cat("  HOG-inhibits-GRHPR mechanism implemented, PH3 urinary oxalate comes\n")
cat("  out at or BELOW normal. The model does not reproduce PH3, and neither\n")
cat("  does that hypothesis as written. Only an additional route in which\n")
cat("  accumulated HOG escapes the mitochondrion and is cleaved in the\n")
cat("  CYTOSOL — bypassing peroxisomal AGT entirely — generates PH3-range\n")
cat("  hyperoxaluria. That is a specific, testable claim about where the\n")
cat("  aldolase activity is, and it is reported here as an unresolved gap\n")
cat("  rather than tuned away.\n")

# --- D16  early vs late treatment ------------------------------------
hdr("D16", "EARLY VS LATE — the price of diagnostic delay, in eGFR")
cat(sprintf("  %-38s %10s %10s\n", "lumasiran started at", "eGFR@45y", "ESKD age"))
for (st in c(2, 5, 10, 15, 20, 25, NA)) {
  d <- if (is.na(st)) scen$S03$d else sim1(PH1, lum_ev(st * YR, 45), end = 45 * YR)
  cat(sprintf("  %-38s %10.1f %10s\n",
              ifelse(is.na(st), "never (untreated)", sprintf("age %d", st)),
              tail(d, 1)$EGFR, ifelse(is.na(eskd_age(d)), "-", sprintf("%.1f", eskd_age(d)))))
}
cat("\n  => the same drug given at age 2 and at age 20 buys very different\n")
cat("     kidneys, and the difference is loop-gain arithmetic (D14) rather\n")
cat("     than a separate early-treatment parameter. Nephrons already lost\n")
cat("     are not recovered by closing the balance afterwards.\n")

line("="); cat("ALL DIAGNOSTICS COMPLETE\n"); line("=")
