# =====================================================================
# Short Bowel Syndrome with Chronic Intestinal Failure (SBS-IF)
# 단장증후군 · 만성 장부전 — QSP model (mrgsolve)
#
#   71 ODEs · time unit = DAYS · 28 scenarios · 17 diagnostics
#
# =====================================================================
# THE ORGANISING IDEA
# ---------------------------------------------------------------------
# Almost every other model in this library has a disease-severity state
# that a drug pushes down, and a clinical score that follows it. SBS-IF
# is structurally different, and the difference is the whole model:
#
#     PARENTERAL SUPPORT VOLUME IS AN OUTPUT, NOT AN INPUT.
#
# It is the arithmetic residual of two conservation equations that are
# written out explicitly and never bypassed:
#
#   FLUID   dTBW/dt = ORAL + PN − OUTPUT − INSENSIBLE − URINE
#   SODIUM  dNAB/dt = ORAL_NA + PN_NA − OUT_NA − URINE_NA
#   ENERGY  dE/dt   = ABSORBED + SCFA_SALVAGE + PN_KCAL − TEE
#
# Everything else on the map earns its place by changing one term in one
# of those equations. Villus height is not an endpoint; the litre is.
#
# ---------------------------------------------------------------------
# SIX STRUCTURAL COMMITMENTS (and what each one generates)
# ---------------------------------------------------------------------
#
# 1. TWO LUMINAL STREAMS, NOT ONE.
#    Fluid taken WITH meals and fluid drunk BETWEEN meals are separate
#    streams with separate sodium concentrations and separate glucose
#    content. They are averaged by volume share only after each has been
#    given its own driving force. This is the only reason the model can
#    generate the central clinical paradox of the disease: a thirsty
#    end-jejunostomy patient who drinks 2 L of PLAIN WATER between meals
#    loses sodium and volume, while the same 2 L as glucose-saline
#    (Na ≈ 90 mmol/L) is absorbed. A single-stream model has to assert
#    that rule; this one produces it as a sign change.
#
# 2. GLUCOSE LOWERS THE ZERO-FLUX SET POINT; IT DOES NOT ADD A FLUX.
#    Jejunal net sodium movement is zero at luminal [Na+] ≈ 90-100
#    mmol/L (Fordtran's perfusion studies) and negative below it. SGLT1
#    cotransport is therefore entered as a SHIFT IN CEQ (down to ~60
#    mmol/L at saturating glucose), not as an additive absorptive term.
#    Entering it additively double-counts sodium that NHE3 immediately
#    recycles, and produces absurd fluxes; entering it as a set-point
#    shift reproduces both the ORS effect and its saturability.
#
# 3. NUTRIENT-SPECIFIC REFERENCE LENGTHS.
#    Carbohydrate and protein are largely absorbed in the first 100-150
#    cm; fat needs 200-300 cm; bile salts and B12 need terminal ileum
#    and nothing else will do. One "absorptive surface" scalar cannot
#    reproduce the observed ordering (CHO ~70-90% absorbed, protein
#    ~60-80%, fat ~30-50%) in the same patient. Three reference lengths
#    can, and they also place each micronutrient deficiency at a named
#    segment rather than at a generic severity.
#
# 4. ADAPTATION IS A PRODUCT, NOT A SUM.
#    TROPHIC = LUMINAL_NUTRIENT x (1 + GLP2R occupancy term). With no
#    enteral intake the product is zero, so a nil-by-mouth patient on
#    full PN neither adapts NOR responds to a GLP-2 analogue. Trophic
#    feeding is a multiplicand, not an adjunct. Diagnostic D07 runs this.
#
# 5. THE PRESCRIPTION IS A CLOSED-LOOP CONTROLLER RUNNING THE REAL
#    TRIAL PROTOCOL. In STEPS the investigator did not choose the PN
#    reduction; an algorithm did — cut PN when 24-hour URINE OUTPUT rose
#    >=10% above baseline, hold if weight or electrolytes were falling.
#    That algorithm is implemented here as a rate-limited bidirectional
#    controller on the PNVOL state. Consequences:
#      - the registered endpoint (>=20% PN volume reduction) is produced
#        by SIMULATING THE PROTOCOL, not by fitting a PN-reduction
#        parameter to the trial;
#      - the controller is self-limiting (cut PN -> urine falls back ->
#        cutting stops), which is why real PN reduction plateaus;
#      - and the placebo arm's substantial response emerges, because a
#        urine-triggered protocol weans anybody carrying unused
#        absorptive reserve, drug or no drug. See D10.
#
# 6. THE ILEAL-BRAKE POSITIVE FEEDBACK IS CLOSED.
#    L-cells live in distal ileum and proximal colon — the tissue that
#    was resected. Losing them removes the brake, transit accelerates,
#    contact time falls, absorption falls, and the unabsorbed load
#    bypasses the L-cells that would have sensed it. Contact time is a
#    MULTIPLICAND in every absorptive fraction, which is why opioids
#    (which only slow transit) have real quantitative effect here.
#
# ---------------------------------------------------------------------
# WHAT THE MODEL IS ASKED TO REPRODUCE (published anchors)
# ---------------------------------------------------------------------
#   A1  STEPS (teduglutide 0.05 mg/kg/d, 24 wk): 63% vs 30% responders
#       (>=20% PN volume reduction at wk 20 AND 24); mean PN volume
#       change -4.4 vs -2.3 L/wk from baseline ~12.9 L/wk.
#   A2  STEPS-2/-3 long term: mean PN reduction grows to ~7-8 L/wk by
#       ~30 months; a minority reach full enteral autonomy.
#   A3  Colon as a digestive organ (Nordgaard): colonic salvage worth up
#       to ~1000 kcal/d; patients with colon in continuity need roughly
#       half the parenteral support of end-jejunostomy patients.
#   A4  Messing anatomical thresholds for permanent intestinal failure:
#       end-jejunostomy ~115 cm, jejunocolic ~60 cm, jejuno-ileo-colic
#       ~35 cm of remnant small bowel.
#   A5  Crenn: plasma citrulline < 20 umol/L marks permanent IF; SBS-IF
#       patients sit at ~10-20 vs ~30-40 umol/L normal.
#   A6  Jejunal sodium zero-flux at luminal [Na+] ~90-100 mmol/L;
#       therapeutic ORS is formulated at 90-120 mmol/L for that reason.
#   A7  Stomal effluent sodium ~90-100 mmol/L (so output volume is a
#       sodium loss rate, which is why urinary Na < 10 mmol/L is the
#       earliest depletion marker).
#   A8  CRBSI in home parenteral nutrition ~0.5-2 episodes per 1000
#       catheter-days; taurolidine/ethanol locks and a dedicated team
#       are the largest modifiers.
#   A9  Fish-oil-based lipid emulsion (phytosterol ~0) reverses IFALD
#       cholestasis over ~2-4 months; soybean-oil ILE (~350 ug/mL
#       phytosterols) sustains it.
#   A10 Enteric hyperoxaluria requires a colon in continuity; an
#       end-jejunostomy patient forms dehydration/low-citrate stones
#       instead. Two different stone diseases, two different therapies.
#
# ---------------------------------------------------------------------
# HONEST LIMITATIONS — read before using any number
# ---------------------------------------------------------------------
#   * This is a whole-body daily-average model. It has no meal timing,
#     no circadian structure and no intraluminal spatial gradient. The
#     two-stream construction is a coarse stand-in for meal timing and
#     should not be read as a real pharmacokinetic compartmentalisation
#     of the lumen.
#   * The fluid model is calibrated to reproduce plausible OUTPUT and PN
#     volumes for stated anatomies; it is not a validated absorption
#     model and cannot be used to prescribe.
#   * Several coefficients (jejunal FMAXJ, colonic KCOLMAX, adaptation
#     time constants) are structurally identifiable only as a group.
#     Different combinations reproduce the same anchors.
#   * Survival, transplantation and quality of life are crude hazard
#     accumulators included so the competing-risk STRUCTURE is visible.
#     Their absolute values should not be quoted.
#   * Six results contradict the model's own design or the literature it
#     was built from. They are reported in DIAGNOSTIC D17 rather than
#     removed: the cholestyramine trap comes out backwards, the population
#     drug-minus-placebo difference is over-predicted even though the
#     individual anchor matches, nutrient gating attenuates rather than
#     abolishes the drug, octreotide looks better than it is, baseline
#     drift is ~0.5%/year rather than zero, and effluent sodium runs about
#     10 mmol/L low. Read D17 before quoting any number from this file.
#
# NOTE ON USE: educational / research QSP model. Not validated for
# clinical decision-making, prescribing or regulatory submission.
# =====================================================================

library(mrgsolve)

# ---------------------------------------------------------------------
# MODEL SPECIFICATION
# ---------------------------------------------------------------------
sbs_code <- '
$PROB
# Short Bowel Syndrome with chronic Intestinal Failure (SBS-IF) — QSP model
# 71 ODEs. Time unit = DAYS. See sbs_mrgsolve_model.R header for documentation.

$SET end=730, delta=1, rtol=1e-6, atol=1e-8, maxsteps=500000

$GLOBAL
// ---- small helpers -------------------------------------------------
// The solver can undershoot states that decay over orders of magnitude;
// fractional powers of a negative number are NaN, so clamp first.
double pos(double x){ return (x > 0.0) ? x : 0.0; }
double dmx(double a, double b){ return (a > b) ? a : b; }
double dmn(double a, double b){ return (a < b) ? a : b; }
double clamp2(double x, double lo, double hi){
  return (x < lo) ? lo : ((x > hi) ? hi : x);
}
// saturable fraction of a delivered load
double sat(double x, double k){ return pos(x)/(pos(x) + k); }
// smooth logistic step, used for the weaning-protocol triggers so that
// the ODE system stays differentiable (a hard if() makes lsoda chatter)
double step2(double x, double c, double w){ return 1.0/(1.0 + exp(-(x - c)/w)); }
// single-pass exponential extraction: absorbed fraction of what passes
double extract(double da){ return 1.0 - exp(-pos(da)); }
// sigmoid saturation with exponent 2, for a transporter whose capacity is
// switch-like in the length of bowel that carries it. ASBT is the case that
// matters here: ~100 cm of terminal ileum reclaims almost the whole bile
// acid pool, and a first-order curve cannot represent that.
double hill2(double x, double k){ double a = pos(x)*pos(x); return a/(a + k*k); }

$PARAM @annotated
// ================= ANATOMY: PARAMETERS, NOT STATES =================
SBL       :  80    : Remnant small bowel length from DJ flexure (cm)
ILEUMLEN  :   0    : Remnant ileum length (cm) - site of ASBT, B12, L-cells
COLONFRAC :  0.5   : Fraction of colon in continuity (0 = end-jejunostomy)
ICV       :   0    : Ileocaecal valve present (0/1)
MUCQUAL   :  1.0   : Quality of remnant mucosa (1 normal, <1 radiation/Crohn)
STOMACH   :   1    : Stomach intact (0/1)

// ================= DIET AND DRINKING BEHAVIOUR =====================
KCALORAL  : 2400   : Reference oral energy intake (kcal/d) before hyperphagia
FCHO_DIET :  0.50  : Fraction of oral energy as carbohydrate
FFAT_DIET :  0.35  : Fraction of oral energy as fat
FPRO_DIET :  0.15  : Fraction of oral energy as protein
FIBERKCAL :  120   : Fermentable fibre delivered to colon (kcal/d)
DIETH2O   :  0.80  : Water taken WITH meals (L/d)
DIETNA    :  160   : Dietary sodium (mmol/d)
DRINKVOL  :  2.50  : Fluid drunk BETWEEN meals (L/d)
ORSFRAC   :  0.40  : Fraction of between-meal fluid taken as ORS
ORSNA     :   90   : Sodium concentration of the ORS (mmol/L)
ORSGLU    :   75   : Glucose concentration of the ORS (mmol/L)
WATERNA   :    5   : Sodium concentration of plain water / juice (mmol/L)
ENTFRAC   :  1.0   : Enteral intake multiplier (0 = nil by mouth)

// ================= SECRETORY LOAD (L/d and mmol/L) =================
VSALIVA   : 1.20   : Salivary volume (L/d)
VGASTRIC  : 2.00   : Gastric secretion at normal tone (L/d)
VBILE     : 0.60   : Biliary volume (L/d)
VPANC     : 1.00   : Pancreatic volume (L/d)
VJEJSEC   : 0.40   : Jejunal crypt secretion (L/d)
CNASAL    :  30    : Salivary sodium (mmol/L)
CNAGAS    : 100    : Gastric sodium (mmol/L)
CNABIL    : 145    : Biliary sodium (mmol/L)
CNAPAN    : 140    : Pancreatic sodium (mmol/L)
CNAJEJ    : 140    : Jejunal secretion sodium (mmol/L)
HYPSEC0   : 1.00   : Initial post-resection gastric hypersecretion multiplier
KHYPSEC   : 0.006  : Decay rate of post-resection hypersecretion (1/d)

// ================= JEJUNAL FLUX ====================================
FMAXJ     : 1.65   : Maximum fractional jejunal volume absorption scalar
KJSURF    : 0.35   : Extraction constant of the jejunal volume flux
GAMMANA   : 0.35   : Weight of the sodium-gradient term in FRACJ
CEQ0      : 100    : Zero-flux luminal sodium set point, no glucose (mmol/L)
DCEQ      :  40    : Maximal downward shift of CEQ by SGLT1 cotransport
KGLU      :  25    : Glucose concentration for half-maximal CEQ shift (mmol/L)
CSCALE    :  45    : Sodium-gradient sensitivity scale (mmol/L)
LREF_NA   : 150    : Reference bowel length for sodium/water handling (cm)
CNAPLASMA : 140    : Plasma sodium used for isotonic water drag (mmol/L)

// ================= NUTRIENT-SPECIFIC ABSORPTION ====================
LREF_CHO  : 120    : Reference length for carbohydrate absorption (cm)
LREF_PRO  : 150    : Reference length for protein absorption (cm)
LREF_FAT  : 250    : Reference length for fat absorption (cm)
KAP_CHO   : 2.20   : Intrinsic extraction constant, carbohydrate
KAP_PRO   : 2.10   : Intrinsic extraction constant, protein
KAP_FAT   : 3.00   : Intrinsic extraction constant, fat
KMBA      : 1.00   : Bile acid pool for half-maximal micellar function (g)
BAREF     : 0.741  : Micellar factor of a normal 2.8 g pool (normalisation)
ACIDPEN   : 0.15   : Maximal fat-absorption penalty from acid hypersecretion

// ================= COLON ===========================================
KCOLMAX   : 6.00   : Maximal colonic water absorption at full colon (L/d)
KMCOL     : 1.50   : Half-saturation of colonic absorption on delivered load (L/d)
KSECR     : 12.0   : Colonic secretagogue load halving colonic absorption
KFASEC    : 0.06   : Secretagogue equivalents per g/d of unabsorbed fat
FERMMAX   : 1000   : Maximal colonic energy salvage at full colon (kcal/d)
KFERM     : 0.020  : Rate constant of fermentative-capacity adaptation (1/d)
SCFAYLD   : 0.55   : Fraction of fermentable energy recovered as SCFA
KENAC     : 1.15   : Colonic sodium avidity multiplier vs water (ENaC)

// ================= MOTILITY / CONTACT TIME =========================
CTMIN     : 0.30   : Minimum relative mucosal contact time
CTMAX     : 1.15   : Maximum relative mucosal contact time
KBRAKE    : 0.55   : Half-saturation of the ileal brake on contact time
KCT       : 0.35   : Rate constant of contact-time state (1/d)
BRAKE0    : 0.25   : Brake tone independent of L-cells
BRPYY     : 1.20   : Weight of PYY/GLP-1 tone in the brake
BROP      : 1.60   : Weight of enteric mu-opioid occupancy in the brake
BROCT     : 0.90   : Weight of octreotide in the brake
BRGLP2    : 0.95   : Weight of GLP-2 (ENS relay) in the brake
ICVGAIN   : 0.15   : Contact-time gain from an intact ileocaecal valve

// ================= L-CELL / ENTEROHORMONAL AXIS ====================
LC_IL     : 1.00   : L-cell mass contributed by remnant ileum (max)
KLC       :  90    : Ileum length for half-maximal L-cell mass (cm)
LC_COL    : 0.45   : L-cell mass contributed by a full colon
KLCAD     : 0.004  : L-cell hyperplasia rate constant (1/d)
LCADMAX   : 1.35   : Maximal L-cell hyperplasia factor
KPYY      : 0.60   : Turnover of the aggregated PYY/GLP-1 brake tone (1/d)

// ================= STRUCTURAL / FUNCTIONAL ADAPTATION ==============
VILLMAX   : 1.50   : Maximal relative villus height
KVILL     : 0.055  : Villus response rate constant (1/d) - tau ~18 d
KTROPH    : 0.70   : Half-saturation of villus target on trophic signal
KCRYPT    : 0.12   : Crypt proliferation turnover (1/d)
KREM      : 0.004  : Slow structural remodelling rate (1/d) - tau ~250 d
REMMAX    : 1.25   : Maximal slow structural remodelling factor
KTRANS    : 0.05   : Transporter/enzyme adaptation rate (1/d)
TRMAX     : 1.45   : Maximal transporter up-regulation
EGLP2     : 1.90   : Maximal trophic gain from GLP2R occupancy
EC50G     : 0.010  : Effect-site GLP-2 concentration for half-maximal trophic gain (mg/L)
EGH       : 0.45   : Maximal trophic gain from systemic IGF-1
KZNTROPH  : 0.35   : Zinc sufficiency needed for crypt mitosis
CITNORM   : 34.0   : Plasma citrulline of an intact small bowel (umol/L)
CITEXP    : 0.55   : Exponent linking enterocyte mass to plasma citrulline

// ================= HYPERPHAGIA / ENERGY ============================
KHP       : 0.015  : Hyperphagia adaptation rate (1/d)
HPMAX     : 1.95   : Maximal hyperphagia multiplier
AHP       : 0.45   : Hyperphagia gain per 1000 kcal/d of faecal energy loss
LEPT      : 0.85   : Strength of the adiposity (leptin) brake on intake
FATREF    : 15.0   : Fat mass above which the adiposity brake engages (kg)
KFAT      : 6.00   : Fat mass for half-maximal adiposity brake (kg)
REEA      : 370    : Cunningham intercept for REE (kcal/d)
REEB      : 21.6   : Cunningham slope on lean mass (kcal/d/kg)
PAL       : 1.55   : Physical activity level
KCALFAT   : 9400   : Energy density of fat mass (kcal/kg)
KCALLEAN  : 1800   : Energy density of lean mass (kcal/kg)
FSURPFAT  : 0.80   : Fraction of energy surplus deposited as fat
FDEFFAT   : 0.70   : Fraction of energy deficit taken from fat

// ================= FLUID BALANCE / KIDNEY ==========================
THIRSTG   : 0.90   : Maximal compensatory increase in oral fluid intake
KTHIRST   : 0.25   : Volume depletion for half-maximal thirst response
KTH       : 0.50   : Rate at which drinking behaviour adapts (1/d)
INSENS    : 0.70   : Insensible water loss (L/d)
URMIN     : 0.40   : Obligatory minimum urine output (L/d)
KREN      : 3.00   : Renal volume-clearing gain (L/d per L of surplus TBW)
TBW0      : 36.0   : Euvolaemic total body water (L)
KURSM     : 0.20   : Smoothing rate for the observed 24-h urine (1/d)
UNAFRAC   : 0.55   : Fraction of urine sodium load excreted at euvolaemia

// ================= PN PRESCRIPTION CONTROLLER (STEPS protocol) =====
WEANON    :  1     : Weaning protocol active (0/1)
URBASE    : 1.20   : Baseline 24-h urine captured at randomisation (L/d)
URTRIG    : 1.10   : Urine ratio above baseline that triggers PN reduction
URRESC    : 0.85   : Urine ratio below which PN is re-escalated
URW       : 0.020  : Width of the smooth weaning trigger
KTAPER    : 0.0060 : Maximum fractional PN volume reduction rate (1/d)
KESCAL    : 0.0250 : PN escalation rate constant (1/d)
PNCAP     : 4.00   : Maximum prescribable PN volume (L/d)
WTGUARD   : 0.95   : Weight fraction of target below which weaning is blocked
WTTARGET  : 62.0   : Target body weight (kg)
KPNKCAL   : 0.070  : Rate at which PN energy tracks the energy deficit (1/d)
PNNACONC  : 130    : Sodium concentration of the PN admixture (mmol/L)

// ================= TEDUGLUTIDE PK ==================================
TEDDOSE   : 0      : Teduglutide dose placeholder (mg) - given as events
KA_TED    : 6.00   : Teduglutide SC absorption rate constant (1/d)
KE_TED    : 8.32   : Teduglutide elimination rate constant (1/d) - t1/2 2 h
V_TED     : 26.0   : Teduglutide apparent volume of distribution (L)
F_TED     : 0.88   : Teduglutide subcutaneous bioavailability
KA_APR    : 0.85   : Apraglutide SC absorption rate constant (1/d)
KE_APR    : 0.554  : Apraglutide elimination rate constant (1/d) - t1/2 30 h
V_APR     : 9.00   : Apraglutide apparent volume of distribution (L)
KA_GLE    : 1.20   : Glepaglutide SC absorption rate constant (1/d)
KE_GLE    : 1.60   : Glepaglutide elimination rate constant (1/d)
V_GLE     : 12.0   : Glepaglutide apparent volume of distribution (L)
KE_NAT    : 142    : Native GLP-2 elimination rate constant (1/d) - t1/2 7 min
V_NAT     : 12.0   : Native GLP-2 volume of distribution (L)
POT_APR   : 1.30   : Molar potency of apraglutide relative to teduglutide
POT_GLE   : 0.55   : Molar potency of glepaglutide relative to teduglutide
POT_NAT   : 1.00   : Molar potency of native GLP-2 relative to teduglutide
KE0G      : 2.00   : GLP-2 effect-site equilibration rate constant (1/d)
KADA      : 0.010  : Anti-drug antibody formation rate (1/d per unit exposure)
KADAOUT   : 0.020  : Anti-drug antibody loss rate (1/d)
ADAPOT    : 0.25   : Fraction of GLP-2 signal neutralised at ADA = 1

// ================= OTHER DRUGS =====================================
KA_GH     : 2.20   : Somatropin SC absorption rate constant (1/d)
KE_GH     : 8.00   : Somatropin elimination rate constant (1/d)
V_GH      : 12.0   : Somatropin volume of distribution (L)
KIGF      : 0.25   : IGF-1 turnover rate constant (1/d)
EMAXIGF   : 1.80   : Maximal GH-driven rise in IGF-1 (fold)
EC50GH    : 3.00   : GH concentration for half-maximal IGF-1 rise (ng/mL)
KE_LOP    : 1.40   : Loperamide elimination rate constant (1/d)
V_LOP     : 65.0   : Loperamide volume of distribution (L)
EC50LOP   : 0.030  : Loperamide concentration for half-maximal brake (mg/L)
KE_OCT    : 5.50   : Octreotide elimination rate constant (1/d)
V_OCT     : 20.0   : Octreotide volume of distribution (L)
EC50OCT   : 0.0015 : Octreotide concentration for half-maximal effect (mg/L)
OCTSEC    : 0.25   : Maximal fractional reduction of secretory volume by octreotide
OCTFAT    : 0.25   : Maximal fractional loss of fat absorption on octreotide
OCTIGF    : 0.50   : Maximal fractional suppression of IGF-1 by octreotide
KPPI      : 0.90   : Proton pump inactivation rate per unit PPI dose (1/d)
KPPIREC   : 0.30   : Proton pump regeneration rate (1/d)
PPIMAXSUP : 0.60   : Maximal fractional reduction of gastric volume by PPI
KE_CHOL   : 3.00   : Cholestyramine luminal transit-out rate (1/d)
CHOLCAP   : 0.55   : Maximal fraction of luminal bile acid sequestered
EC50CHOL  : 2.00   : Cholestyramine luminal amount for half-maximal binding (g)
KE_RIF    : 2.50   : Rifaximin luminal transit-out rate (1/d)
RIFKILL   : 1.30   : Maximal rifaximin bacterial kill rate (1/d)
EC50RIF   : 200    : Rifaximin luminal amount for half-maximal kill (mg)
KE_GLUT   : 4.00   : Oral glutamine luminal disappearance rate (1/d)
GLUTTROPH : 0.10   : Maximal trophic gain from luminal glutamine
EC50GLUT  : 8.00   : Glutamine luminal amount for half-maximal effect (g)
CLONSEC   : 0.0    : Clonidine fractional reduction of jejunal secretion (0-1)

// ================= BILE ACIDS ======================================
BASYNMAX  : 5.00   : Maximal hepatic bile acid synthesis (g/d)
BASYN0    : 0.50   : Basal hepatic bile acid synthesis (g/d)
NCYC      : 6.00   : Enterohepatic cycles per day
FRECMAX   : 0.97   : Maximal ileal bile acid reclamation fraction
KASBT     :  45    : Ileum length for half-maximal ASBT reclamation (cm)
FRECJEJ   : 0.30   : Passive jejunal bile acid reclamation fraction
KFGF      : 0.55   : Reclamation fraction for half-maximal CYP7A1 feedback

// ================= SIBO / MICROBIAL ================================
SIBOGROW  : 0.16   : Maximal SIBO growth rate (1/d)
SIBOCAP   : 1.00   : SIBO carrying capacity (relative units)
SIBOCLR   : 0.20   : Baseline SIBO clearance by MMC and acid (1/d)
KDLACT    : 3.00   : D-lactate turnover rate (1/d)
DLGEN     : 0.004  : D-lactate generated per kcal of colonic fermentation
KETOX     : 2.00   : Endotoxin turnover rate (1/d)
ETXGEN    : 0.80   : Endotoxin generation per unit SIBO x barrier failure
KMUCINJ   : 0.06   : Bacterial mucosal injury turnover (1/d)
MUCINJMAX : 0.20   : Maximal fractional villus loss from bacterial injury

// ================= MICRONUTRIENTS ==================================
OUTREF    : 3.00   : Intestinal output that doubles micronutrient loss (L/d)
MGORAL    :  15    : Oral magnesium intake (mmol/d)
MGPN      :  5.0   : Magnesium in PN per litre (mmol/L)
MGPNRET   : 0.40   : Fraction of intravenous magnesium retained
MGTHR     : 0.75   : Renal magnesium reabsorption threshold (mmol/L)
MGENDOL   : 0.80   : Magnesium in intestinal secretions (mmol per L of output)
MGOBL     : 1.00   : Obligatory daily magnesium loss (mmol/d)
MGSET     : 0.85   : Serum magnesium defended by skeletal buffering (mmol/L)
KMGBUF    : 0.350  : Skeletal magnesium buffering rate (1/d)
FMGJEJ    : 0.05   : Jejunal magnesium absorption fraction
FMGIL     : 0.35   : Maximal ileal magnesium absorption fraction
KMGIL     :  70    : Ileum length for half-maximal Mg absorption (cm)
KMGREN    : 1.20   : Renal magnesium clearance above threshold (1/d)
VDMG      : 20.0   : Magnesium distribution volume (L)
ZNORAL    : 180    : Oral zinc intake (umol/d)
ZNPN      : 100    : Zinc in PN per litre (umol/L)
KZNOUT    : 0.10   : Zinc turnover (1/d)
B12ORAL   : 3.00   : Oral vitamin B12 intake (ug/d)
B12PN     : 5.00   : Vitamin B12 in PN per litre (ug/L)
KB12IL    :  35    : Ileum length for half-maximal B12 absorption (cm)
KB12OUT   : 0.0007 : Vitamin B12 turnover (1/d) - very large body store
VITDORAL  :  20    : Oral vitamin D intake (ug/d)
VITDPN    :  5.0   : Vitamin D in PN per litre (ug/L)
KVITD     : 0.020  : 25-OH vitamin D turnover (1/d)
VITAORAL  : 900    : Oral vitamin A intake (ug/d)
KVITA     : 0.010  : Vitamin A turnover (1/d)
VITEORAL  :  15    : Oral vitamin E intake (mg/d)
KVITE     : 0.010  : Vitamin E turnover (1/d)
SEORAL    :  70    : Oral selenium intake (ug/d)
SEPN      :  60    : Selenium in PN per litre (ug/L)
KSEOUT    : 0.020  : Selenium turnover (1/d)
EFAORAL   :  14    : Oral essential fatty acid intake (g/d)
EFAPN     :  20    : EFA delivered per litre of lipid-containing PN (g/L)
KEFA      : 0.030  : Essential fatty acid turnover (1/d)
FLIPIDPN  : 1.00   : Lipid content of the PN admixture (1 = standard, 0 = lipid-free)

// ================= MINERAL / BONE ==================================
KPTHIN    : 6.00   : PTH secretion rate constant
KPTHOUT   : 12.0   : PTH clearance rate constant (1/d)
KMGPTH    : 0.30   : Serum Mg for half-maximal PTH secretion (mmol/L)
CATARGET  : 2.35   : Target serum calcium (mmol/L)
KCAVITD   : 0.55   : Weight of vitamin D on calcium absorption
KBMDLOSS  : 0.00025: Baseline bone mineral density loss rate (1/d)
KBMDACID  : 0.00060: Additional BMD loss per unit acidosis
KBMDVITD  : 0.00090: Additional BMD loss per unit vitamin D deficit
HCO3TGT   : 24.0   : Target serum bicarbonate (mmol/L)
KHCO3     : 0.50   : Bicarbonate turnover (1/d)
HCO3OUT   : 1.20   : Bicarbonate lost per litre of intestinal output (mmol/L)

// ================= IFALD ===========================================
PHYTOCONC : 350    : Phytosterol concentration of the lipid emulsion (ug/mL)
ILEDOSE   : 0.90   : Lipid emulsion dose (g/kg/d) when PN lipid is given
KPHYTO    : 0.030  : Hepatic phytosterol elimination rate (1/d)
KMPHYTO   : 420    : Phytosterol load for half-maximal cholestatic drive (mg)
KCHOLD    : 0.020  : Choline deficit accrual/repletion rate (1/d)
CHOLSUP   : 0      : Choline supplementation (0/1)
KHEPIN    : 0.012  : Hepatic inflammation formation rate (1/d)
KHEPOUT   : 0.030  : Hepatic inflammation resolution rate (1/d)
WPHYTO    : 0.55   : Weight of phytosterol load on cholestatic drive
WCHOLD    : 0.35   : Weight of choline deficit on cholestatic drive
WNOENT    : 0.45   : Weight of absent enteral stimulation on cholestatic drive
WSEPSIS   : 1.60   : Weight of each septic episode on hepatic inflammation
WETX      : 0.40   : Weight of portal endotoxin on hepatic inflammation
OMEGA3    : 0      : Fish-oil (omega-3) content flag of the emulsion (0/1)
OM3ANTI   : 0.55   : Fractional reduction of hepatic inflammation by omega-3
KBILI     : 0.10   : Bilirubin turnover (1/d)
KMHEP     : 2.20   : Hepatic inflammation for half-maximal bilirubin rise
BILIMAX   : 22.0   : Maximal bilirubin (mg/dL)
KFIB      : 0.0016 : Liver fibrosis progression rate constant (1/d)
KFIBREV   : 0.0006 : Liver fibrosis regression rate constant (1/d)

// ================= KIDNEY / STONE ==================================
UOX0      : 30.0   : Baseline urinary oxalate (mg/d)
KOXFAT    : 130    : Urinary oxalate gain per unit fat malabsorption x colon
KOXBA     : 0.00035: Oxalate permeability gain per mmol/d colonic bile acid
KUOX      : 1.00   : Urinary oxalate smoothing rate (1/d)
OXCITRATE : 0      : Potassium citrate / low-oxalate programme (0/1)
CITEFF    : 0.40   : Fractional reduction of stone formation by citrate
KSTONE    : 0.00012: Stone formation rate per (mg/d oxalate above threshold)
UOXTHR    : 40.0   : Urinary oxalate threshold for stone formation (mg/d)
KSTONECONC: 0.55   : Extra stone risk from concentrated urine
KSTONEPASS: 0.0025 : Rate at which stone burden passes or is removed (1/d)
GFR0      : 95.0   : Baseline eGFR (mL/min/1.73m2)
KGFRAGE   : 0.0025 : Age-related eGFR decline (per day)
KGFRDEHY  : 0.030  : eGFR decline per unit chronic volume depletion
KGFROX    : 0.0160 : eGFR decline per unit stone / oxalate burden
KAKI      : 0.60   : Acute kidney injury accrual rate from volume depletion (1/d)
KAKIREC   : 0.25   : AKI recovery rate (1/d)

// ================= CATHETER / OUTCOME ==============================
VMAXNIGHT : 2.50   : Maximum PN volume infused in one night (L)
CRBSIRATE : 1.20   : CRBSI episodes per 1000 catheter-days, unmodified
LOCKFAC   : 1.00   : Catheter lock factor (taurolidine/ethanol ~0.25)
TECHFAC   : 1.00   : Aseptic technique / dedicated team factor
KTRANSLOC : 0.60   : CRBSI rate gain from gut bacterial translocation
VEINS0    : 6.00   : Number of usable central venous access sites
KVEINLOSS : 0.35   : Access sites consumed per CRBSI or thrombosis event
KPOLYP    : 0.00040: Colonic polyp accrual rate per unit GLP-2 signal (1/d)
KPOLYP0   : 0.00012: Background colonic polyp accrual rate (1/d)
KQOL      : 0.05   : Quality-of-life adaptation rate (1/d)
WQOLPN    : 0.055  : QoL penalty per L/d of PN volume
WQOLOUT   : 0.030  : QoL penalty per L/d of intestinal output
WQOLINF   : 0.070  : QoL penalty per CRBSI episode per year
KHAZ      : 0.000045: Baseline mortality hazard (1/d)
HAZLIVER  : 0.00055: Mortality hazard per unit liver fibrosis
HAZINF    : 0.00030: Mortality hazard per CRBSI episode per year
HAZCKD    : 0.00025: Mortality hazard per unit renal impairment

$INIT @annotated
// Presenting condition of a prevalent SBS-IF patient. These are STATIC:
// the burnin() helper in the R driver replaces them with a self-consistent
// state. Values assigned in $MAIN would be re-applied at every new
// individual and would silently override init(), discarding that state -
// which is exactly the bug this block exists to avoid.
// ---- DRUG PK / PD (18) --------------------------------------------
TEDSC   : 0     : Teduglutide SC depot (mg)
TEDC    : 0     : Teduglutide central (mg)
APRSC   : 0     : Apraglutide SC depot (mg)
APRC    : 0     : Apraglutide central (mg)
GLESC   : 0     : Glepaglutide SC depot (mg)
GLEC    : 0     : Glepaglutide central (mg)
NATC    : 0     : Native GLP-2 central (mg)
GLP2E   : 0     : GLP-2 effect-site signal (mg/L equivalent)
GHSC    : 0     : Somatropin SC depot (mg)
GHC     : 0     : Somatropin central (mg)
IGF1    : 1     : Systemic IGF-1 (relative to normal)
LOPC    : 0     : Loperamide central (mg)
OCTC    : 0     : Octreotide central (mg)
PPIP    : 1     : Active gastric proton pump pool (fraction)
CHOLL   : 0     : Cholestyramine luminal amount (g)
RIFL    : 0     : Rifaximin luminal amount (mg)
ADA     : 0     : Anti-drug antibody (relative)
GLUTL   : 0     : Luminal glutamine (g)
// ---- MUCOSA / ADAPTATION / MOTILITY (12) --------------------------
VILLUS  : 1     : Villus height (relative to the remnant baseline)
CRYPT   : 1     : Crypt proliferation index (relative)
MUCREM  : 1     : Slow structural remodelling factor
SGLT1R  : 1     : SGLT1 density (relative)
NHE3R   : 1     : NHE3 density (relative)
BBENZ   : 1     : Brush-border enzyme capacity (relative)
FERMCAP : 500   : Colonic fermentative capacity (kcal/d at full colon)
LCELL   : 1     : L-cell mass (relative)
PYYT    : 1     : Aggregated PYY/GLP-1 ileal brake tone (relative)
HYPSEC  : 1     : Gastric hypersecretion multiplier
CONTACT : 0.6   : Relative mucosal contact time
THIRST  : 0     : Compensatory thirst / drinking behaviour (0-1)
CITRP   : 16    : Plasma citrulline (umol/L)
// ---- LUMINAL ECOLOGY AND BILE ACIDS (5) ---------------------------
BAPOOL  : 2.8   : Bile acid pool (g)
SIBO    : 0.1   : Small intestinal bacterial overgrowth (relative)
DLACT   : 0     : D-lactate (relative)
ENDOTOX : 0     : Portal endotoxin (relative)
MUCINJ  : 0     : Bacterial mucosal injury (fraction of maximum)
// ---- CONSERVATION / BODY COMPOSITION / PRESCRIPTION (10) ----------
TBW     : 36    : Total body water (L)
NABODY  : 2100  : Exchangeable body sodium (mmol)
LEAN    : 45    : Lean body mass (kg)
FATM    : 15    : Fat mass (kg)
HYPERPH : 1.3   : Hyperphagia multiplier
PNVOL   : 0     : Prescribed PN volume (L/d) - found by the controller
PNKCAL  : 0     : Prescribed PN energy (kcal/d)
URINESM : 1     : Smoothed 24-h urine output (L/d)
CUMNRG  : 0     : Cumulative energy balance (kcal)
HCO3    : 24    : Serum bicarbonate (mmol/L)
// ---- MICRONUTRIENTS AND MINERAL-BONE (10) -------------------------
MGS     : 0.85  : Serum magnesium (mmol/L)
ZNB     : 1400  : Body zinc pool (umol)
B12B    : 3000  : Body vitamin B12 pool (ug)
VITDB   : 800   : 25-OH vitamin D pool (relative units)
VITAB   : 60000 : Body vitamin A pool (ug)
VITEB   : 1200  : Body vitamin E pool (mg)
SEB     : 3000  : Body selenium pool (ug)
EFAB    : 450   : Essential fatty acid pool (g)
PTHS    : 45    : Parathyroid hormone (pg/mL)
BMD     : 1     : Bone mineral density (relative)
// ---- LIVER: IFALD (5) ---------------------------------------------
PHYTO   : 0     : Hepatic phytosterol load (mg)
CHOLDEF : 0     : Choline deficit (relative)
HEPINF  : 0.05  : Hepatic inflammation (relative)
BILI    : 0.8   : Total bilirubin (mg/dL)
FIB     : 0.2   : Liver fibrosis stage (0-4)
// ---- KIDNEY AND STONE (4) -----------------------------------------
UOXS    : 30    : Urinary oxalate (mg/d)
STONE   : 0     : Stone burden (0-1)
GFR     : 95    : eGFR (mL/min/1.73m2)
AKID    : 0     : Acute kidney injury (relative)
// ---- CATHETER AND OUTCOME (6) -------------------------------------
CATHD   : 0     : Cumulative catheter-days
CRBSIC  : 0     : Cumulative CRBSI episodes
VEINS   : 6     : Patent central venous access sites
POLYP   : 0     : Colonic polyp burden (relative)
QOL     : 0.6   : Quality of life (0-1)
SURVH   : 0     : Cumulative mortality hazard

$MAIN
F_TEDSC = F_TED;

$ODE
// =====================================================================
// SECTION 1 — DRUG PK
// =====================================================================
double CP_TED = TEDC / V_TED;
double CP_APR = APRC / V_APR;
double CP_GLE = GLEC / V_GLE;
double CP_NAT = NATC / V_NAT;
double CP_GH  = GHC  / V_GH;
double CP_LOP = LOPC / V_LOP;
double CP_OCT = OCTC / V_OCT;

dxdt_TEDSC = -KA_TED * TEDSC;
dxdt_TEDC  =  KA_TED * TEDSC - KE_TED * TEDC;
dxdt_APRSC = -KA_APR * APRSC;
dxdt_APRC  =  KA_APR * APRSC - KE_APR * APRC;
dxdt_GLESC = -KA_GLE * GLESC;
dxdt_GLEC  =  KA_GLE * GLESC - KE_GLE * GLEC;
dxdt_NATC  = -KE_NAT * NATC;
dxdt_GHSC  = -KA_GH  * GHSC;
dxdt_GHC   =  KA_GH  * GHSC - KE_GH  * GHC;
dxdt_LOPC  = -KE_LOP * LOPC;
dxdt_OCTC  = -KE_OCT * OCTC;

// Aggregate GLP-2 receptor drive in teduglutide-equivalent concentration.
// Potency ratios carry the pharmacology of each analogue; the effect-site
// compartment GLP2E carries the fact that the MUCOSAL response integrates
// over days while the drug itself lives for hours (structural point 5 of
// the header: tau_PD >> tau_PK is why once-daily dosing of a 2 h drug works).
double GLP2DRIVE = CP_TED + POT_APR * CP_APR + POT_GLE * CP_GLE + POT_NAT * CP_NAT;
dxdt_GLP2E = KE0G * (GLP2DRIVE - GLP2E);

// Anti-drug antibodies: formed in proportion to exposure, partially
// neutralising. Reported ADA incidence with teduglutide is high but the
// antibodies are mostly non-neutralising, hence ADAPOT well below 1.
dxdt_ADA = KADA * GLP2DRIVE - KADAOUT * ADA;
double ADAFAC = 1.0 - ADAPOT * sat(ADA, 1.0);

double OCCG = ADAFAC * sat(GLP2E, EC50G);

// Somatropin -> systemic IGF-1; octreotide suppresses the same axis, which
// is the mechanistic basis of the octreotide trade-off (S16/D08).
double OCTOCC = sat(CP_OCT, EC50OCT);
double IGFTGT = 1.0 + (EMAXIGF - 1.0) * sat(CP_GH, EC50GH);
dxdt_IGF1 = KIGF * (IGFTGT * (1.0 - OCTIGF * OCTOCC) - IGF1);

// PPI: irreversible pump inactivation with regeneration (PPIP = active pool)
dxdt_PPIP = KPPIREC * (1.0 - PPIP) - KPPI * PPIDOSE * PPIP;

// Luminal, non-absorbed agents
dxdt_CHOLL = -KE_CHOL * CHOLL;
dxdt_RIFL  = -KE_RIF  * RIFL;
dxdt_GLUTL = -KE_GLUT * GLUTL;

// =====================================================================
// SECTION 2 — SECRETORY LOAD
// =====================================================================
dxdt_HYPSEC = -KHYPSEC * (HYPSEC - 1.0);

double PPISUP = PPIMAXSUP * (1.0 - PPIP);
double VGAS   = VGASTRIC * STOMACH * HYPSEC * (1.0 - PPISUP);
double VJS    = VJEJSEC * (1.0 - CLONSEC);
double VSECR  = (VSALIVA + VGAS + VBILE + VPANC + VJS) * (1.0 - OCTSEC * OCTOCC);
double NASECR = (VSALIVA*CNASAL + VGAS*CNAGAS + VBILE*CNABIL + VPANC*CNAPAN
                 + VJS*CNAJEJ) * (1.0 - OCTSEC * OCTOCC);

// =====================================================================
// SECTION 3 — ORAL INTAKE, THE TWO LUMINAL STREAMS
// =====================================================================
double DEHYDX = clamp2((TBW0 - TBW) / 3.0, 0.0, 1.5);
dxdt_THIRST = KTH * (sat(DEHYDX, KTHIRST) - THIRST);
// Compensatory drinking. This is the loop that lets a long remnant
// rescue itself as PN is withdrawn: PN down -> dry -> thirst -> more oral
// fluid -> if the bowel can absorb it, urine recovers and weaning
// continues; if it cannot, PN re-escalates. Without it, no patient can
// ever be weaned by the protocol no matter how much bowel they have.
double VMEALW  = DIETH2O  * ENTFRAC;
double VDRINK  = DRINKVOL * ENTFRAC * (1.0 + THIRSTG * THIRST);
double VORS    = VDRINK * ORSFRAC;
double VPLAIN  = VDRINK * (1.0 - ORSFRAC);
double NAORAL  = DIETNA * ENTFRAC + VORS * ORSNA + VPLAIN * WATERNA;
double VORAL   = VMEALW + VDRINK;

// STREAM M: meal water plus the secretions the meal provokes. Glucose is
// abundant here, so the SGLT1 set-point shift is at its maximum.
double VSTRM   = VMEALW + VSECR;
double NASTRM  = DIETNA * ENTFRAC + NASECR;
double CNA_M   = (VSTRM > 1e-6) ? NASTRM / VSTRM : CEQ0;
double CEQ_M   = CEQ0 - DCEQ * SGLT1R / (SGLT1R + 0.35);

// STREAM D: fluid drunk BETWEEN meals. Its sodium and glucose are whatever
// the patient chose to drink. THIS is where the plain-water trap lives.
double VSTRD   = VDRINK;
double NASTRD  = VORS * ORSNA + VPLAIN * WATERNA;
double CNA_D   = (VSTRD > 1e-6) ? NASTRD / VSTRD : CEQ0;
double GLU_D   = (VSTRD > 1e-6) ? (VORS * ORSGLU) / VSTRD : 0.0;
double CEQ_D   = CEQ0 - DCEQ * sat(GLU_D, KGLU) * SGLT1R / (SGLT1R + 0.35);

double VDEL    = VSTRM + VSTRD;                 // total volume reaching jejunum
double NADEL   = NASTRM + NASTRD;               // total sodium reaching jejunum
double WM      = (VDEL > 1e-6) ? VSTRM / VDEL : 0.0;
double WD      = 1.0 - WM;

// Volume-share-weighted driving force. Each stream is given its OWN
// (CNA - CEQ) before averaging; averaging first would destroy the sign
// change that generates the plain-water paradox.
double DRIVE   = WM * (CNA_M - CEQ_M) + WD * (CNA_D - CEQ_D);

// =====================================================================
// SECTION 4 — MUCOSAL SURFACE, CONTACT TIME, ABSORPTIVE FRACTIONS
// =====================================================================
double VILLEFF = VILLUS * MUCREM * (1.0 - MUCINJMAX * MUCINJ);
double SURF_NA = (SBL / LREF_NA)  * VILLEFF * MUCQUAL * NHE3R;
double SURF_CH = (SBL / LREF_CHO) * VILLEFF * MUCQUAL * BBENZ;
double SURF_PR = (SBL / LREF_PRO) * VILLEFF * MUCQUAL * BBENZ;
double SURF_FA = (SBL / LREF_FAT) * VILLEFF * MUCQUAL;

// Contact time. The ileal brake is the sum of endogenous L-cell tone,
// enteric opioid occupancy, octreotide, and the GLP-2 ENS relay.
double LOPOCC  = sat(CP_LOP, EC50LOP);
double BRAKE   = BRAKE0 + BRPYY*PYYT + BROP*LOPOCC + BROCT*OCTOCC + BRGLP2*OCCG;
double CTTGT   = (CTMIN + (CTMAX - CTMIN) * sat(BRAKE, KBRAKE))
                 * (1.0 - ICVGAIN * (1.0 - ICV));
dxdt_CONTACT = KCT * (CTTGT - CONTACT);
double FC = CONTACT;

// Bile-acid-dependent micellar function, normalised so that a normal
// 2.8 g pool gives 1.0. This is the term the 100 cm rule acts through.
double MICELLE = sat(BAPOOL, KMBA) / BAREF;
MICELLE = clamp2(MICELLE, 0.0, 1.15);
double ACIDFAC = 1.0 - ACIDPEN * (HYPSEC - 1.0) / dmx(HYPSEC0 - 1.0, 1e-6)
                 * step2(HYPSEC, 1.05, 0.05);
ACIDFAC = clamp2(ACIDFAC, 1.0 - ACIDPEN, 1.0);

double FCHO = extract(KAP_CHO * SURF_CH * FC);
double FPRO = extract(KAP_PRO * SURF_PR * FC);
// Octreotide inhibits pancreatic exocrine secretion as well as the
// secretory volume it is given for, so it buys volume at the price of fat.
double FFAT = extract(KAP_FAT * SURF_FA * FC) * MICELLE * ACIDFAC
              * (1.0 - OCTFAT * OCTOCC);
FCHO = clamp2(FCHO, 0.0, 0.99);
FPRO = clamp2(FPRO, 0.0, 0.99);
FFAT = clamp2(FFAT, 0.0, 0.97);

// =====================================================================
// SECTION 5 — JEJUNAL VOLUME AND SODIUM FLUX (THE SIGN CHANGE)
// =====================================================================
// Single-pass extraction, not a Michaelis saturation. With sat() the
// fraction was already 70% of maximal at 80 cm, so remnant LENGTH had
// almost no leverage and no bowel was ever long enough to wean - the
// Messing thresholds could not appear. extract() restores the leverage.
double SURFFAC = extract(KJSURF * SURF_NA * FC);
double FRACJ   = FMAXJ * SURFFAC + GAMMANA * tanh(DRIVE / CSCALE);
FRACJ = clamp2(FRACJ, -0.60, 0.95);

double VJEJ  = VDEL * FRACJ;      // negative = net jejunal SECRETION
double NAJEJ = NADEL * FRACJ;     // isotonic: sodium follows the same fraction

// =====================================================================
// SECTION 6 — COLONIC SALVAGE
// =====================================================================
double VCOLIN  = dmx(VDEL - VJEJ, 0.0);
double NACOLIN = dmx(NADEL - NAJEJ, 0.0);
// unabsorbed fat, needed here as a colonic secretagogue and again in
// section 7 as an energy loss
double FATMALX = KCALORAL * HYPERPH * ENTFRAC * FFAT_DIET * (1.0 - FFAT) / 9.0;

// Bile acid spilling into the colon is a secretory drive that CONSUMES
// colonic capacity - this is how bile acid diarrhoea enters the balance.
double FREC    = FRECJEJ + (FRECMAX - FRECJEJ) * hill2(ILEUMLEN, KASBT);
double CHOLBIND= CHOLCAP * sat(CHOLL, EC50CHOL);
double FREC_EFF= clamp2(FREC * (1.0 - CHOLBIND), 0.0, FRECMAX);
// TWO different quantities, and conflating them is what stopped the first
// draft from reproducing the 100 cm rule:
//   BASPILL  = bile acid LOST FROM THE POOL. A sequestrant INCREASES it,
//              because bound bile acid cannot be reclaimed by ASBT.
//   BAFREE   = UNBOUND bile acid presented to the colonocyte, i.e. the
//              secretory stimulus. A sequestrant DECREASES it.
// Cholestyramine therefore moves the two in opposite directions, and the
// sign of its net effect depends on how much ileum was resected.
double BASPILL = BAPOOL * NCYC * (1.0 - FREC_EFF);              // g/d
double BAFREEG = BAPOOL * NCYC * (1.0 - FREC) * (1.0 - CHOLBIND);
double BACOLMM = BAFREEG * 1000.0 / 450.0;                      // mmol/d, MW ~450
// Colonic capacity is consumed by BOTH unbound bile acid and unabsorbed
// long-chain fatty acids; both are secretagogues acting on the colonocyte.
double SECLOAD = (BACOLMM + KFASEC * FATMALX) * COLONFRAC;
double BAFAC   = 1.0 / (1.0 + SECLOAD / KSECR);

double ALDOFAC = 1.0 + 0.25 * step2(TBW0 - TBW, 0.8, 0.4);
double SCFAFAC = 1.0 + 0.15 * (FERMCAP / dmx(FERMMAX, 1e-6));
double KCOL    = KCOLMAX * COLONFRAC * ALDOFAC * SCFAFAC * BAFAC;
double VCOL    = KCOL * sat(VCOLIN, KMCOL);
VCOL = dmn(VCOL, 0.95 * VCOLIN);
double NACOL   = dmn(NACOLIN * (VCOL / dmx(VCOLIN, 1e-6)) * KENAC, 0.95 * NACOLIN);

double OUTPUT  = dmx(VCOLIN - VCOL, 0.05);        // stool / stoma volume, L/d
double NAOUT   = dmx(NACOLIN - NACOL, 1.0);       // sodium lost in output, mmol/d

// =====================================================================
// SECTION 7 — ENERGY BALANCE
// =====================================================================
double ORALKCAL = KCALORAL * HYPERPH * ENTFRAC;
double ABSKCAL  = ORALKCAL * (FCHO_DIET*FCHO + FFAT_DIET*FFAT + FPRO_DIET*FPRO);

// Only unabsorbed CARBOHYDRATE and fibre are fermentable. Unabsorbed fat
// is not salvaged - it is lost, and on the way out it drives steatorrhoea,
// oxalate absorption and divalent cation losses.
double FERMIN   = ORALKCAL * FCHO_DIET * (1.0 - FCHO) + FIBERKCAL * ENTFRAC;
double SCFAKCAL = dmn(FERMCAP * COLONFRAC, FERMIN * SCFAYLD);
double FATMAL   = FATMALX;                                     // g/d unabsorbed fat

double REE  = REEA + REEB * LEAN;
double TEE  = REE * PAL;
double NRGBAL = ABSKCAL + SCFAKCAL + PNKCAL - TEE;
dxdt_CUMNRG = NRGBAL;

double FFATDEP = (NRGBAL >= 0.0) ? FSURPFAT : FDEFFAT;
dxdt_FATM = FFATDEP * NRGBAL / KCALFAT;
dxdt_LEAN = (1.0 - FFATDEP) * NRGBAL / KCALLEAN;
double WEIGHT = LEAN + FATM;

// Hyperphagia is a real, quantitatively large compensation in SBS: intake
// rises toward twice normal when absorption fails.
// Hyperphagia in SBS tracks FAECAL energy loss, not net balance after PN:
// the hunger signal is enteral, which is why patients on adequate
// parenteral energy are hyperphagic too. A leptin-like adiposity brake
// closes the loop; without it, faecal loss and intake reinforce each other
// without limit.
double FECALNET = dmx(ORALKCAL - ABSKCAL - SCFAKCAL, 0.0);
double LEPTBRK  = dmx(1.0 - LEPT * sat(dmx(FATM - FATREF, 0.0), KFAT), 0.0);
double HPTGT = clamp2(1.0 + AHP * (FECALNET / 1000.0) * LEPTBRK, 1.0, HPMAX);
dxdt_HYPERPH = KHP * (HPTGT - HYPERPH);

// =====================================================================
// SECTION 8 — FLUID AND SODIUM BALANCE, URINE
// =====================================================================
double URINE = URMIN + KREN * dmx(TBW - TBW0, 0.0);
dxdt_TBW    = VORAL + PNVOL - OUTPUT - INSENS - URINE;
dxdt_URINESM = KURSM * (URINE - URINESM);

double NAPN   = PNVOL * PNNACONC;
double NAURINE= UNAFRAC * dmx(NABODY - 2100.0, 0.0) * 0.9
                + 8.0 * step2(TBW - TBW0, -0.2, 0.5) + 2.0;
dxdt_NABODY = NAORAL + NAPN - NAOUT - NAURINE;

double DEHYD = clamp2((TBW0 - TBW) / 3.0, 0.0, 1.5);

// =====================================================================
// SECTION 9 — THE PN PRESCRIPTION CONTROLLER (STEPS weaning protocol)
// =====================================================================
// Reduce PN when the smoothed 24-h urine has risen >=10% above the
// randomisation baseline; block the reduction if weight has fallen below
// its guardrail; re-escalate if urine falls or the patient is dry.
double URATIO  = URINESM / dmx(URBASE, 1e-6);
double WTOK    = step2(WEIGHT, WTGUARD * WTTARGET, 0.4);
double SDOWN   = step2(URATIO, URTRIG, URW) * WTOK * WEANON;
double SUP     = step2(URRESC, URATIO, URW) + step2(DEHYD, 0.35, 0.05);
SUP = dmn(SUP, 1.0);
dxdt_PNVOL = -KTAPER * PNVOL * SDOWN + KESCAL * dmx(PNCAP - PNVOL, 0.0) * SUP;

// PN energy tracks the measured energy deficit with a weight-error term.
double PNKREQ = dmx(TEE - ABSKCAL - SCFAKCAL, 0.0)
                + 90.0 * dmx(WTTARGET - WEIGHT, 0.0);
dxdt_PNKCAL = KPNKCAL * (PNKREQ - PNKCAL);

// =====================================================================
// SECTION 10 — L-CELL AXIS AND ADAPTATION (THE PRODUCT, NOT THE SUM)
// =====================================================================
double LUMNUT = (ORALKCAL / dmx(KCALORAL, 1e-6));       // 0 when nil by mouth
double LCBASE = LC_IL * sat(ILEUMLEN, KLC) + LC_COL * COLONFRAC;
dxdt_LCELL = KLCAD * (LCBASE * LCADMAX * sat(LUMNUT, 0.5) - LCELL);
dxdt_PYYT  = KPYY * (LCELL * sat(LUMNUT, 0.45) / 0.735 - PYYT);

double ZNSUF  = sat(ZNB / 1400.0, KZNTROPH);
double GLUTFX = 1.0 + GLUTTROPH * sat(GLUTL, EC50GLUT);
double TROPHIC = LUMNUT * (1.0 + EGLP2 * OCCG) * (1.0 + EGH * (IGF1 - 1.0))
                 * GLUTFX * ZNSUF;

double VTGT = 1.0 + (VILLMAX - 1.0) * sat(TROPHIC, KTROPH);
dxdt_VILLUS = KVILL * (VTGT - VILLUS);
dxdt_CRYPT  = KCRYPT * (VTGT - CRYPT);
dxdt_MUCREM = KREM * (1.0 + (REMMAX - 1.0) * sat(TROPHIC, KTROPH) - MUCREM);

double TRTGT = 1.0 + (TRMAX - 1.0) * sat(TROPHIC, KTROPH);
dxdt_SGLT1R = KTRANS * (TRTGT - SGLT1R);
dxdt_NHE3R  = KTRANS * (TRTGT - NHE3R);
dxdt_BBENZ  = KTRANS * (TRTGT - BBENZ);
dxdt_FERMCAP = KFERM * (FERMMAX * sat(FERMIN / 400.0, 0.8) - FERMCAP);

// Plasma citrulline: a biomarker of ENTEROCYTE MASS, not of severity.
double MUCREL = (SBL * VILLEFF * MUCQUAL) / 400.0;
double CITTGT = CITNORM * pow(dmx(MUCREL, 1e-6), CITEXP) * (1.0 + 0.18 * OCCG);
dxdt_CITRP = 1.20 * (CITTGT - CITRP);

// =====================================================================
// SECTION 11 — BILE ACIDS
// =====================================================================
double BASYN = BASN(BASYNMAX, BASYN0, FREC_EFF, KFGF);
dxdt_BAPOOL = BASYN - BAPOOL * NCYC * (1.0 - FREC_EFF);

// =====================================================================
// SECTION 12 — SIBO, D-LACTATE, ENDOTOXIN, MUCOSAL INJURY
// =====================================================================
double CLRMMC  = SIBOCLR * (0.55 + 0.45 * ICV) * (0.60 + 0.40 * PPIP)
                 * (1.0 + RIFKILL * sat(RIFL, EC50RIF));
dxdt_SIBO = SIBOGROW * SIBO * (1.0 - SIBO / SIBOCAP)
            * (1.0 + 0.8 * (1.0 - ICV)) - CLRMMC * SIBO;

// D-lactic acidosis REQUIRES a colon: the substrate is colonic
// carbohydrate fermentation. In an end-jejunostomy it cannot happen.
dxdt_DLACT = DLGEN * FERMIN * COLONFRAC * sat(SIBO, 0.35) - KDLACT * DLACT;
dxdt_ENDOTOX = ETXGEN * SIBO * (1.0 + MUCINJ) - KETOX * ENDOTOX;
dxdt_MUCINJ = KMUCINJ * (sat(SIBO, 0.45) - MUCINJ);

// =====================================================================
// SECTION 13 — MICRONUTRIENTS
// =====================================================================
// Every balance below writes LOSS AS A FIRST-ORDER FUNCTION OF THE BODY
// POOL, with a rate constant that rises with intestinal output. Writing
// the loss as (stool concentration x output volume) instead - the obvious
// way, and the way the first draft of this model did it - makes the pool
// unbounded below, because nothing then stops the model excreting
// magnesium or zinc that the patient does not have.
double OUTFAC = 1.0 + OUTPUT / OUTREF;

// --- magnesium. Serum Mg is the state and the skeleton buffers it, which
// is why SBS hypomagnesaemia is chronic and slow rather than abrupt.
double FMG   = FMGJEJ + FMGIL * sat(ILEUMLEN, KMGIL);
double MGABS = MGORAL * ENTFRAC * FMG / OUTFAC;
double MGRET = MGPN * PNVOL * MGPNRET;                 // most IV Mg is excreted
double MGREN = KMGREN * VDMG * dmx(MGS - MGTHR, 0.0);  // renal threshold ~0.75
double MGENDO= MGENDOL * OUTPUT * step2(MGS, 0.22, 0.04);  // cannot lose what
                                                       // is no longer there
dxdt_MGS = (MGABS + MGRET - MGREN - MGENDO - MGOBL) / VDMG
           + KMGBUF * (MGSET - MGS);
double SMG = MGS;

dxdt_ZNB = ZNORAL * ENTFRAC * (0.20 + 0.25 * sat(SBL/150.0, 0.6)) / OUTFAC
           + ZNPN * PNVOL - KZNOUT * OUTFAC * ZNB;

double FB12 = 0.02 + 0.55 * sat(ILEUMLEN, KB12IL) * (1.0 - 0.5 * sat(SIBO, 0.4));
dxdt_B12B = B12ORAL * ENTFRAC * FB12 * 1000.0 + B12PN * PNVOL - KB12OUT * B12B;

dxdt_VITDB = VITDORAL * ENTFRAC * FFAT * 2.2 + VITDPN * PNVOL * 3.0 - KVITD * VITDB;
dxdt_VITAB = VITAORAL * ENTFRAC * FFAT * 1.5 + 800.0 * PNVOL - KVITA * VITAB;
dxdt_VITEB = VITEORAL * ENTFRAC * FFAT * 5.0 + 10.0 * PNVOL - KVITE * VITEB;
dxdt_SEB   = (SEORAL * ENTFRAC * (0.4 + 0.4*sat(SBL/150.0,0.6)) / OUTFAC
              + SEPN * PNVOL) * 1.2 - KSEOUT * SEB;
dxdt_EFAB  = EFAORAL * ENTFRAC * FFAT + EFAPN * PNVOL * FLIPIDPN - KEFA * EFAB;

// =====================================================================
// SECTION 14 — ACID-BASE, PTH, CALCIUM, BONE
// =====================================================================
dxdt_HCO3 = KHCO3 * (HCO3TGT - HCO3) - HCO3OUT * OUTPUT * 0.5 - 0.8 * DLACT;
double ACIDOSIS = clamp2((HCO3TGT - HCO3) / 6.0, 0.0, 1.5);

// Magnesium-dependent PTH: hypomagnesaemia impairs PTH SECRETION and
// causes end-organ PTH resistance, so hypocalcaemia stays refractory
// until magnesium is replaced. That is generated here, not asserted.
double VITDREL = clamp2(VITDB / 800.0, 0.0, 1.5);
double CAABS   = (0.25 + KCAVITD * VITDREL) * FFAT / 0.45;
double CASERUM = CATARGET * clamp2(0.55 + 0.30*CAABS + 0.22*sat(PTHS/45.0,0.8)
                                   * sat(SMG, KMGPTH) / 0.68, 0.5, 1.1);
dxdt_PTHS = KPTHIN * dmx(CATARGET - CASERUM, 0.0) * 10.0 * sat(SMG, KMGPTH)
            - KPTHOUT * PTHS;

dxdt_BMD = -BMD * (KBMDLOSS + KBMDACID * ACIDOSIS
                   + KBMDVITD * dmx(1.0 - VITDREL, 0.0));

// =====================================================================
// SECTION 15 — IFALD
// =====================================================================
double ILEGIVEN = ILEDOSE * WEIGHT * FLIPIDPN * PNDAYSF;
dxdt_PHYTO = ILEGIVEN * PHYTOCONC / 1000.0 - KPHYTO * PHYTO;

double PNENFRAC = PNKCAL / dmx(ABSKCAL + SCFAKCAL + PNKCAL, 1e-6);
dxdt_CHOLDEF = KCHOLD * (PNENFRAC * (1.0 - CHOLSUP) - CHOLDEF);

double NOENT = 1.0 - clamp2(LUMNUT, 0.0, 1.0);
// Nights on the line per week. A litre saved is not a catheter-day saved;
// a NIGHT saved is. Each infusion night carries at most VMAXNIGHT litres,
// so as prescribed volume falls the team drops whole nights - and that,
// not the volume itself, is what changes infection and access risk.
double PNDAYSF = clamp2(PNVOL * 7.0 / VMAXNIGHT, 0.0, 7.0) / 7.0;
double CRBSIR = (CRBSIRATE / 1000.0) * LOCKFAC * TECHFAC
                * (1.0 + KTRANSLOC * sat(ENDOTOX, 0.9)) * PNDAYSF;
double CHOLDRIVE = WPHYTO * sat(PHYTO, KMPHYTO) + WCHOLD * CHOLDEF + WNOENT * NOENT
                   + WSEPSIS * CRBSIR * 30.0 + WETX * sat(ENDOTOX, 1.0);
dxdt_HEPINF = KHEPIN * CHOLDRIVE * (1.0 - OM3ANTI * OMEGA3) - KHEPOUT * HEPINF;
dxdt_BILI = KBILI * (BILIMAX * sat(HEPINF, KMHEP) + 0.7 - BILI);
dxdt_FIB  = KFIB * sat(HEPINF, 0.8) * (4.0 - FIB) - KFIBREV * FIB;

// =====================================================================
// SECTION 16 — KIDNEY, OXALATE, STONES
// =====================================================================
// Enteric hyperoxaluria REQUIRES a colon. Multiply, do not add.
double UOXTGT = UOX0 + KOXFAT * (FATMAL / 100.0) * COLONFRAC
                * (1.0 + KOXBA * BACOLMM);
dxdt_UOXS = KUOX * (UOXTGT - UOXS);
double CONCURINE = clamp2(1.4 - URINE, 0.0, 1.0);
// Stone burden saturates and stones pass. Without both terms the burden
// grows without bound and drags eGFR straight through the floor.
dxdt_STONE = (KSTONE * dmx(UOXS - UOXTHR, 0.0)
              + KSTONECONC * KSTONE * 40.0 * CONCURINE)
             * (1.0 - CITEFF * OXCITRATE) * dmx(1.0 - STONE, 0.0)
             - KSTONEPASS * STONE;

dxdt_AKID = KAKI * dmx(DEHYD - 0.30, 0.0) - KAKIREC * AKID;
dxdt_GFR = -(KGFRAGE + KGFRDEHY * DEHYD + KGFROX * STONE + 0.02 * AKID)
           * step2(GFR, 8.0, 2.0);

// =====================================================================
// SECTION 17 — CATHETER, ACCESS, POLYPS, QOL, SURVIVAL
// =====================================================================
dxdt_CATHD = PNDAYSF;
dxdt_CRBSIC = CRBSIR;
dxdt_VEINS = -KVEINLOSS * CRBSIR;
dxdt_POLYP = (KPOLYP0 + KPOLYP * OCCG) * COLONFRAC;

double RENIMP = clamp2((GFR0 - GFR) / GFR0, 0.0, 1.0);
double QTGT = clamp2(1.0 - WQOLPN * PNVOL * 7.0 / 7.0 - WQOLOUT * OUTPUT
                     - WQOLINF * CRBSIR * 365.0 - 0.10 * sat(FIB, 2.0)
                     - 0.08 * sat(DLACT, 1.0) - 0.06 * RENIMP, 0.05, 1.0);
dxdt_QOL = KQOL * (QTGT - QOL);

dxdt_SURVH = KHAZ + HAZLIVER * FIB + HAZINF * CRBSIR * 365.0 + HAZCKD * RENIMP;

$TABLE
double CP_TEDo = TEDC / V_TED;
double CP_APRo = APRC / V_APR;
double CP_GLEo = GLEC / V_GLE;
double VILLEFFo = VILLUS * MUCREM * (1.0 - MUCINJMAX * MUCINJ);

// --- recompute the reportable balance quantities ---------------------
double PPISUPo = PPIMAXSUP * (1.0 - PPIP);
double OCTOCCo = sat(OCTC / V_OCT, EC50OCT);
double VGASo   = VGASTRIC * STOMACH * HYPSEC * (1.0 - PPISUPo);
double VJSo    = VJEJSEC * (1.0 - CLONSEC);
double VSECRo  = (VSALIVA + VGASo + VBILE + VPANC + VJSo) * (1.0 - OCTSEC*OCTOCCo);
double NASECRo = (VSALIVA*CNASAL + VGASo*CNAGAS + VBILE*CNABIL + VPANC*CNAPAN
                  + VJSo*CNAJEJ) * (1.0 - OCTSEC*OCTOCCo);

double VMEALWo = DIETH2O * ENTFRAC;
double VDRINKo = DRINKVOL * ENTFRAC * (1.0 + THIRSTG * THIRST);
double VORSo   = VDRINKo * ORSFRAC;
double VPLAINo = VDRINKo * (1.0 - ORSFRAC);
double NAORALo = DIETNA*ENTFRAC + VORSo*ORSNA + VPLAINo*WATERNA;
double VORALo  = VMEALWo + VDRINKo;

double VSTRMo  = VMEALWo + VSECRo;
double NASTRMo = DIETNA*ENTFRAC + NASECRo;
double CNA_Mo  = (VSTRMo > 1e-6) ? NASTRMo/VSTRMo : CEQ0;
double CEQ_Mo  = CEQ0 - DCEQ * SGLT1R/(SGLT1R + 0.35);
double NASTRDo = VORSo*ORSNA + VPLAINo*WATERNA;
double CNA_Do  = (VDRINKo > 1e-6) ? NASTRDo/VDRINKo : CEQ0;
double GLU_Do  = (VDRINKo > 1e-6) ? (VORSo*ORSGLU)/VDRINKo : 0.0;
double CEQ_Do  = CEQ0 - DCEQ * sat(GLU_Do, KGLU) * SGLT1R/(SGLT1R + 0.35);
double VDELo   = VSTRMo + VDRINKo;
double NADELo  = NASTRMo + NASTRDo;
double WMo     = (VDELo > 1e-6) ? VSTRMo/VDELo : 0.0;
double DRIVEo  = WMo*(CNA_Mo - CEQ_Mo) + (1.0-WMo)*(CNA_Do - CEQ_Do);

double SURF_NAo = (SBL/LREF_NA) * VILLEFFo * MUCQUAL * NHE3R;
double SURF_CHo = (SBL/LREF_CHO)* VILLEFFo * MUCQUAL * BBENZ;
double SURF_PRo = (SBL/LREF_PRO)* VILLEFFo * MUCQUAL * BBENZ;
double SURF_FAo = (SBL/LREF_FAT)* VILLEFFo * MUCQUAL;

double ADAFACo = 1.0 - ADAPOT * sat(ADA, 1.0);
double OCCGo   = ADAFACo * sat(GLP2E, EC50G);
double LOPOCCo = sat(LOPC/V_LOP, EC50LOP);
double BRAKEo  = BRAKE0 + BRPYY*PYYT + BROP*LOPOCCo + BROCT*OCTOCCo + BRGLP2*OCCGo;

double MICELLEo = clamp2(sat(BAPOOL, KMBA)/BAREF, 0.0, 1.15);
double ACIDFACo = clamp2(1.0 - ACIDPEN*(HYPSEC-1.0)/dmx(HYPSEC0-1.0,1e-6)
                         * step2(HYPSEC,1.05,0.05), 1.0-ACIDPEN, 1.0);
double FCHOo = clamp2(extract(KAP_CHO*SURF_CHo*CONTACT), 0.0, 0.99);
double FPROo = clamp2(extract(KAP_PRO*SURF_PRo*CONTACT), 0.0, 0.99);
double FFATo = clamp2(extract(KAP_FAT*SURF_FAo*CONTACT)*MICELLEo*ACIDFACo
                      *(1.0 - OCTFAT*OCTOCCo), 0.0, 0.97);

double FRACJo = clamp2(FMAXJ*extract(KJSURF*SURF_NAo*CONTACT)
                       + GAMMANA*tanh(DRIVEo/CSCALE), -0.60, 0.95);
double VJEJo  = VDELo * FRACJo;
double NAJEJo = NADELo * FRACJo;
double VCOLINo = dmx(VDELo - VJEJo, 0.0);
double NACOLINo= dmx(NADELo - NAJEJo, 0.0);

double FRECo    = FRECJEJ + (FRECMAX-FRECJEJ)*hill2(ILEUMLEN,KASBT);
double CHOLBo   = CHOLCAP * sat(CHOLL, EC50CHOL);
double FREC_EFFo= clamp2(FRECo*(1.0-CHOLBo), 0.0, FRECMAX);
double BASPILLo = BAPOOL*NCYC*(1.0-FREC_EFFo);
double BACOLMMo = BAPOOL*NCYC*(1.0-FRECo)*(1.0-CHOLBo)*1000.0/450.0;
double FATMALXo = KCALORAL*HYPERPH*ENTFRAC*FFAT_DIET*(1.0-FFATo)/9.0;
double BAFACo   = 1.0/(1.0 + (BACOLMMo + KFASEC*FATMALXo)*COLONFRAC/KSECR);
double ALDOFACo = 1.0 + 0.25*step2(TBW0-TBW, 0.8, 0.4);
double SCFAFACo = 1.0 + 0.15*(FERMCAP/dmx(FERMMAX,1e-6));
double KCOLo    = KCOLMAX*COLONFRAC*ALDOFACo*SCFAFACo*BAFACo;
double VCOLo    = dmn(KCOLo*sat(VCOLINo,KMCOL), 0.95*VCOLINo);
double NACOLo   = dmn(NACOLINo*(VCOLo/dmx(VCOLINo,1e-6))*KENAC, 0.95*NACOLINo);
double OUTPUTo  = dmx(VCOLINo - VCOLo, 0.05);
double NAOUTo   = dmx(NACOLINo - NACOLo, 1.0);

double ORALKCALo = KCALORAL*HYPERPH*ENTFRAC;
double ABSKCALo  = ORALKCALo*(FCHO_DIET*FCHOo + FFAT_DIET*FFATo + FPRO_DIET*FPROo);
double FERMINo   = ORALKCALo*FCHO_DIET*(1.0-FCHOo) + FIBERKCAL*ENTFRAC;
double SCFAKCALo = dmn(FERMCAP*COLONFRAC, FERMINo*SCFAYLD);
double FATMALo   = FATMALXo;
double REEo = REEA + REEB*LEAN;
double TEEo = REEo*PAL;
double NRGBALo = ABSKCALo + SCFAKCALo + PNKCAL - TEEo;
double WEIGHTo = LEAN + FATM;
double URINEo = URMIN + KREN*dmx(TBW - TBW0, 0.0);
double DEHYDo = clamp2((TBW0-TBW)/3.0, 0.0, 1.5);
double SMGo   = MGS;
double VITDRELo = clamp2(VITDB/800.0, 0.0, 1.5);
double ACIDOSISo= clamp2((HCO3TGT-HCO3)/6.0, 0.0, 1.5);
double PNDAYSFo = clamp2(PNVOL*7.0/VMAXNIGHT, 0.0, 7.0)/7.0;
double CRBSIRo = (CRBSIRATE/1000.0)*LOCKFAC*TECHFAC
                 *(1.0 + KTRANSLOC*sat(ENDOTOX,0.9))*PNDAYSFo;
double NAURINEo = UNAFRAC*dmx(NABODY-2100.0,0.0)*0.9
                  + 8.0*step2(TBW-TBW0,-0.2,0.5) + 2.0;
double MUCRELo = (SBL*VILLEFFo*MUCQUAL)/400.0;
double LUMNUTo = ORALKCALo/dmx(KCALORAL,1e-6);
double TROPHICo = LUMNUTo*(1.0+EGLP2*OCCGo)*(1.0+EGH*(IGF1-1.0))
                  *(1.0+GLUTTROPH*sat(GLUTL,EC50GLUT))*sat(ZNB/1400.0,KZNTROPH);

// --- reportable outputs ---------------------------------------------
capture PNVOLWK  = PNVOL * 7.0;             // L/week - the registered endpoint
capture PNKCALD  = PNKCAL;
capture PNNIGHTS = PNDAYSFo * 7.0;          // infusion nights per week
capture OUT_LD   = OUTPUTo;                 // stoma/stool output L/d
capture URINE_LD = URINEo;
capture URSM     = URINESM;
capture STOMANA  = NAOUTo / dmx(OUTPUTo, 1e-6);   // effluent [Na+] mmol/L
capture UNA      = NAURINEo;                      // urinary sodium mmol/d
capture CNALUM   = NADELo / dmx(VDELo, 1e-6);
capture DRIVEC   = DRIVEo;                        // volume-weighted mean
capture DRIVEM   = CNA_Mo - CEQ_Mo;               // meal stream: always positive
capture DRIVED   = CNA_Do - CEQ_Do;               // drink stream: SIGN CHANGES
capture CNAD     = CNA_Do;
capture SECLOADC = (BACOLMMo + KFASEC*FATMALXo)*COLONFRAC;
capture FRACJC   = FRACJo;
capture VJEJC    = VJEJo;
capture VCOLC    = VCOLo;
capture FCHOC    = FCHOo;
capture FPROC    = FPROo;
capture FFATC    = FFATo;
capture ABSKC    = ABSKCALo;
capture SCFAKC   = SCFAKCALo;
capture TEEC     = TEEo;
capture NRGBALC  = NRGBALo;
capture WT       = WEIGHTo;
capture CITR     = CITRP;
capture VILLC    = VILLEFFo;
capture MUCRELC  = MUCRELo;
capture TROPHC   = TROPHICo;
capture OCCGC    = OCCGo;
capture CPTED    = CP_TEDo;
capture CPAPR    = CP_APRo;
capture CPGLE    = CP_GLEo;
capture BRAKEC   = BRAKEo;
capture CONTC    = CONTACT;
capture BAPOOLC  = BAPOOL;
capture BACOLC   = BACOLMMo;
capture MICELC   = MICELLEo;
capture SIBOC    = SIBO;
capture DLACTC   = DLACT;
capture SMGC     = SMGo;
capture PTHC     = PTHS;
capture BMDC     = BMD;
capture HCO3C    = HCO3;
capture ACIDC    = ACIDOSISo;
capture VITDC    = VITDRELo;
capture B12C     = B12B;
capture ZNC      = ZNB;
capture EFAC     = EFAB;
capture BILIC    = BILI;
capture FIBC     = FIB;
capture PHYTOC   = PHYTO;
capture HEPINFC  = HEPINF;
capture CRBSIYR  = CRBSIRo * 365.0;
capture CRBSICUM = CRBSIC;
capture VEINSC   = VEINS;
capture UOXC     = UOXS;
capture STONEC   = STONE;
capture GFRC     = GFR;
capture DEHYDC   = DEHYDo;
capture QOLC     = QOL;
capture SURVP    = 100.0 * exp(-SURVH);
capture POLYPC   = POLYP;
capture TBWC     = TBW;
capture NABC     = NABODY;
capture FATMALC  = FATMALo;
capture THIRSTC  = THIRST;
capture VORALC   = VORALo;
'

# ---------------------------------------------------------------------
# The one non-inline function used in $ODE, plus the PPI dose parameter,
# have to exist before mcode() compiles. BASN() implements FGF19-mediated
# feedback inhibition of CYP7A1: when ileal reclamation fails, the
# feedback signal disappears and synthesis rises toward its ceiling.
# ---------------------------------------------------------------------
sbs_code <- sub(
  "$GLOBAL\n// ---- small helpers",
  paste0("$GLOBAL\n",
         "double BASN(double smax, double s0, double frec, double kf){\n",
         "  // FGF19 feedback strength rises with ileal reclamation fraction\n",
         "  double fb = frec/(frec + kf);\n",
         "  return s0 + (smax - s0) * (1.0 - fb);\n",
         "}\n",
         "// ---- small helpers"),
  sbs_code, fixed = TRUE)

sbs_code <- sub("PNNACONC  : 130    : Sodium concentration of the PN admixture (mmol/L)",
  paste0("PNNACONC  : 130    : Sodium concentration of the PN admixture (mmol/L)\n",
         "PPIDOSE   : 0      : PPI dose intensity (0 none, 1 standard, 2 high dose)"),
  sbs_code, fixed = TRUE)

mod <- mcode("sbs", sbs_code, soloc = tempdir())

# =====================================================================
# HELPERS
# =====================================================================
`%||%` <- function(a, b) if (is.null(a)) b else a
dmax_r <- function(x) if (x <= 0) 1e-9 else x

# Run a model to (near) steady state and return the final state vector,
# so that every scenario starts from a self-consistent patient rather
# than from a hand-written initial condition. This is what makes the
# claim "the disease is generated, not assumed" checkable (see D01).
CMTNAMES <- names(mrgsolve::init(mod)@data)

# c(list(SBL = 80, ICV = 1), list(ICV = 0)) produces a list with TWO ICV
# entries and param() silently takes the FIRST, so an override written that
# way never happens. This bug quietly disabled the protocol effect in the
# virtual-population diagnostic and the ICV-loss scenarios before it was
# found. Rather than fix ~40 call sites, every parameter list entering the
# solver is de-duplicated here, LAST occurrence winning (modifyList
# semantics), which is what the caller always meant.
mp <- function(pars) {
  pars <- as.list(pars)
  if (!length(pars)) return(pars)
  pars[!duplicated(names(pars), fromLast = TRUE)]
}

# States that are CUMULATIVE or PROGRESSIVE by construction (organ injury,
# catheter history, fibrosis, stones, bone, hazard) cannot be equilibrated -
# they have no steady state, so a long burn-in would silently manufacture a
# patient with 30 years of iatrogenic damage. They are therefore reset to a
# STATED presenting condition after the physiological states have settled,
# and every scenario's organ trajectory is then attributable to that
# scenario rather than to the burn-in.
PRESENT <- list(PHYTO = 0, CHOLDEF = 0, HEPINF = 0.05, BILI = 0.8, FIB = 0.2,
                STONE = 0, GFR = 95, AKID = 0, BMD = 1.0, CATHD = 0,
                CRBSIC = 0, VEINS = 6, POLYP = 0, CUMNRG = 0, SURVH = 0)

burnin <- function(pars = list(), days = 900, present = PRESENT) {
  pars <- mp(pars)
  m <- if (length(pars)) param(mod, pars) else mod
  d <- as.data.frame(mrgsim(m, end = days, delta = days / 4))
  st <- as.numeric(d[nrow(d), CMTNAMES])
  names(st) <- CMTNAMES
  for (nm in names(present)) st[[nm]] <- present[[nm]]
  st
}

# pull the terminal state out of a finished simulation
laststate <- function(d) {
  st <- as.numeric(d[nrow(d), CMTNAMES]); names(st) <- CMTNAMES; st
}

# Simulate a scenario from a burned-in state.
sim <- function(state, pars = list(), events = NULL, end = 168, delta = 1) {
  pars <- mp(pars)
  m <- if (length(pars)) param(mod, pars) else mod
  m <- init(m, as.list(state))
  e <- if (is.null(events)) ev(amt = 0, cmt = "TEDSC", time = 0) else events
  as.data.frame(mrgsim(m, e, end = end, delta = delta))
}

# Teduglutide 0.05 mg/kg/day SC
ev_ted <- function(wt = 62, days = 168, start = 0, dose_mgkg = 0.05) {
  ev(amt = dose_mgkg * wt, cmt = "TEDSC", time = start,
     ii = 1, addl = days - 1)
}
# Apraglutide 5 mg SC weekly
ev_apra <- function(days = 168, start = 0, dose = 5) {
  ev(amt = dose, cmt = "APRSC", time = start, ii = 7,
     addl = max(0, floor(days / 7) - 1))
}
# Glepaglutide 10 mg SC weekly (or twice weekly with ii = 3.5)
ev_glep <- function(days = 168, start = 0, dose = 10, ii = 7) {
  ev(amt = dose, cmt = "GLESC", time = start, ii = ii,
     addl = max(0, floor(days / ii) - 1))
}
# Native GLP-2, equimolar to teduglutide, once or twice daily SC.
# Modelled as entering the central compartment directly because the
# 7-minute half-life makes absorption kinetics irrelevant.
ev_native <- function(wt = 62, days = 168, per_day = 1, dose_mgkg = 0.05) {
  ii <- 1 / per_day
  ev(amt = dose_mgkg * wt / per_day, cmt = "NATC", time = 0,
     ii = ii, addl = days * per_day - 1)
}
# Somatropin 0.1 mg/kg/day SC for 28 days
ev_gh <- function(wt = 62, days = 28, dose_mgkg = 0.1) {
  ev(amt = dose_mgkg * wt, cmt = "GHSC", time = 0, ii = 1, addl = days - 1)
}
# Loperamide 4 mg qid
ev_lop <- function(days = 168, dose = 4, per_day = 4) {
  ev(amt = dose, cmt = "LOPC", time = 0, ii = 1 / per_day,
     addl = days * per_day - 1)
}
# Octreotide 100 ug tid SC
ev_oct <- function(days = 168, dose_ug = 100, per_day = 3) {
  ev(amt = dose_ug / 1000, cmt = "OCTC", time = 0, ii = 1 / per_day,
     addl = days * per_day - 1)
}
# Cholestyramine 4 g tid (luminal)
ev_chol <- function(days = 168, dose_g = 4, per_day = 3) {
  ev(amt = dose_g, cmt = "CHOLL", time = 0, ii = 1 / per_day,
     addl = days * per_day - 1)
}
# Rifaximin 550 mg bid, 14 days
ev_rifax <- function(days = 14, dose_mg = 550, per_day = 2) {
  ev(amt = dose_mg, cmt = "RIFL", time = 0, ii = 1 / per_day,
     addl = days * per_day - 1)
}
# Oral glutamine 30 g/day
ev_glut <- function(days = 168, dose_g = 10, per_day = 3) {
  ev(amt = dose_g, cmt = "GLUTL", time = 0, ii = 1 / per_day,
     addl = days * per_day - 1)
}

# The three Messing anatomies, as parameter sets
ANAT <- list(
  end_jejunostomy = list(SBL = 90,  ILEUMLEN = 0,  COLONFRAC = 0.00, ICV = 0),
  jejunocolic     = list(SBL = 80,  ILEUMLEN = 0,  COLONFRAC = 0.50, ICV = 0),
  jejunoileocolic = list(SBL = 70,  ILEUMLEN = 25, COLONFRAC = 1.00, ICV = 1)
)

# Endpoint extraction ---------------------------------------------------
resp20 <- function(d, base_wk, wk1 = 140, wk2 = 168) {
  v1 <- d$PNVOLWK[which.min(abs(d$time - wk1))]
  v2 <- d$PNVOLWK[which.min(abs(d$time - wk2))]
  red1 <- (base_wk - v1) / base_wk
  red2 <- (base_wk - v2) / base_wk
  as.numeric(red1 >= 0.20 && red2 >= 0.20 && red1 <= 1.0)
}
at <- function(d, tt, col) d[[col]][which.min(abs(d$time - tt))]
fin <- function(d, col) d[[col]][nrow(d)]

# =====================================================================
# SCENARIOS
# =====================================================================
run_scenarios <- function() {

  base_pars <- ANAT$jejunocolic
  st0 <- burnin(base_pars)
  ur0 <- st0[["URINESM"]]
  pn0 <- st0[["PNVOL"]] * 7
  wt0 <- st0[["LEAN"]] + st0[["FATM"]]
  # the randomisation baseline is captured from the run-in, exactly as the
  # trial did - it is not a free parameter
  cp <- c(base_pars, list(URBASE = ur0, WTTARGET = wt0))

  S <- list()
  add <- function(name, d, note = "") {
    S[[name]] <<- list(sim = d, note = note)
  }

  # ---- S01-S05  natural history and the three anatomies --------------
  add("S01 natural history, no drug, protocol on",
      sim(st0, cp, end = 168))
  for (nm in names(ANAT)) {
    stx <- burnin(ANAT[[nm]])
    cpx <- c(ANAT[[nm]], list(URBASE = stx[["URINESM"]],
                              WTTARGET = stx[["LEAN"]] + stx[["FATM"]]))
    add(paste0("S0x anatomy: ", nm), sim(stx, cpx, end = 168))
  }
  add("S05 nil by mouth on full PN",
      sim(st0, c(cp, list(ENTFRAC = 0.02)), end = 168))

  # ---- S06-S10  GLP-2 analogues --------------------------------------
  add("S06 teduglutide 0.05 mg/kg/d 24 wk",
      sim(st0, cp, ev_ted(wt0, 168), end = 168))
  add("S07 teduglutide 104 wk",
      sim(st0, cp, ev_ted(wt0, 730), end = 730))
  add("S08 apraglutide 5 mg weekly 24 wk",
      sim(st0, cp, ev_apra(168), end = 168))
  add("S09 glepaglutide 10 mg weekly 24 wk",
      sim(st0, cp, ev_glep(168, ii = 7), end = 168))
  add("S10 native GLP-2 equimolar once daily",
      sim(st0, cp, ev_native(wt0, 168, per_day = 1), end = 168))

  # ---- S11-S18  adjuncts ---------------------------------------------
  add("S11 ORS optimisation (ORSFRAC 0.4 -> 0.9)",
      sim(st0, c(cp, list(ORSFRAC = 0.9)), end = 168))
  add("S12 plain water only (ORSFRAC -> 0)",
      sim(st0, c(cp, list(ORSFRAC = 0.0)), end = 168))
  add("S13 high-dose PPI",
      sim(st0, c(cp, list(PPIDOSE = 2)), end = 168))
  add("S14 loperamide 16 mg/d",
      sim(st0, cp, ev_lop(168), end = 168))
  add("S15 loperamide + high-dose PPI + ORS",
      sim(st0, c(cp, list(PPIDOSE = 2, ORSFRAC = 0.9)), ev_lop(168), end = 168))
  add("S16 octreotide 100 ug tid",
      sim(st0, cp, ev_oct(168), end = 168))
  add("S17 somatropin 0.1 mg/kg/d x 4 wk + glutamine",
      sim(st0, cp, c(ev_gh(wt0, 28), ev_glut(168)), end = 168))
  add("S18 teduglutide + optimised ORS + loperamide",
      sim(st0, c(cp, list(ORSFRAC = 0.9)),
          c(ev_ted(wt0, 168), ev_lop(168)), end = 168))

  # ---- S19-S21  bile acids: the 100 cm rule in both directions -------
  st_short <- burnin(c(ANAT$jejunoileocolic, list(ILEUMLEN = 60)))
  cp_short <- c(ANAT$jejunoileocolic, list(ILEUMLEN = 60,
                URBASE = st_short[["URINESM"]],
                WTTARGET = st_short[["LEAN"]] + st_short[["FATM"]]))
  add("S19 short ileal resection (60 cm ileum left) + cholestyramine",
      sim(st_short, cp_short, ev_chol(168), end = 168))
  add("S20 short ileal resection, no cholestyramine",
      sim(st_short, cp_short, end = 168))
  add("S21 extensive ileal resection + cholestyramine (the trap)",
      sim(st0, cp, ev_chol(168), end = 168))

  # ---- S22-S24  SIBO -------------------------------------------------
  st_sibo <- burnin(c(ANAT$jejunoileocolic, list(ICV = 0)))
  cp_sibo <- c(ANAT$jejunoileocolic, list(ICV = 0,
               URBASE = st_sibo[["URINESM"]],
               WTTARGET = st_sibo[["LEAN"]] + st_sibo[["FATM"]]))
  add("S22 ICV lost, SIBO untreated",  sim(st_sibo, cp_sibo, end = 168))
  add("S23 ICV lost + rifaximin 14 d", sim(st_sibo, cp_sibo, ev_rifax(14), end = 168))
  add("S24 carbohydrate overload with colon (D-lactate)",
      sim(st_sibo, c(cp_sibo, list(FCHO_DIET = 0.68, FFAT_DIET = 0.20,
                                   FPRO_DIET = 0.12, FIBERKCAL = 300)),
          end = 168))

  # ---- S25-S28  iatrogenic organ injury over 2 years -----------------
  add("S25 soybean ILE 2 y (IFALD)",
      sim(st0, c(cp, list(PHYTOCONC = 350, OMEGA3 = 0)), end = 730))
  add("S26 SMOF ILE 2 y",
      sim(st0, c(cp, list(PHYTOCONC = 100, OMEGA3 = 1)), end = 730))
  add("S27 fish-oil ILE rescue after 1 y of soybean",
      {
        d1 <- sim(st0, c(cp, list(PHYTOCONC = 350, OMEGA3 = 0)), end = 365)
        d2 <- sim(laststate(d1), c(cp, list(PHYTOCONC = 0, OMEGA3 = 1)), end = 365)
        d2$time <- d2$time + 365
        rbind(d1, d2)
      })
  add("S28 taurolidine lock + dedicated team, 2 y",
      sim(st0, c(cp, list(LOCKFAC = 0.25, TECHFAC = 0.7)), end = 730))

  list(scen = S, st0 = st0, cp = cp, ur0 = ur0, pn0 = pn0, wt0 = wt0)
}

# =====================================================================
# DIAGNOSTICS — each one is a claim the model has to earn
# =====================================================================
run_diagnostics <- function() {
  cat("\n=====================================================================\n")
  cat("SBS-IF QSP MODEL — DIAGNOSTICS\n")
  cat("=====================================================================\n")

  base <- ANAT$jejunocolic
  st0  <- burnin(base)
  ur0  <- st0[["URINESM"]]
  wt0  <- st0[["LEAN"]] + st0[["FATM"]]
  cp   <- c(base, list(URBASE = ur0, WTTARGET = wt0))
  b    <- sim(st0, cp, end = 365)
  pn0  <- at(b, 0, "PNVOLWK")

  # ---- D01  the burned-in patient must not drift ---------------------
  cat("\n[D01] BASELINE STABILITY — is the disease generated or assumed?\n")
  cat(sprintf("  PN volume    %7.2f -> %7.2f L/wk   drift %+.4f%%\n",
      pn0, fin(b, "PNVOLWK"), 100*(fin(b,"PNVOLWK")-pn0)/pn0))
  cat(sprintf("  output       %7.3f -> %7.3f L/d\n", at(b,0,"OUT_LD"), fin(b,"OUT_LD")))
  cat(sprintf("  weight       %7.2f -> %7.2f kg\n", at(b,0,"WT"), fin(b,"WT")))
  cat(sprintf("  citrulline   %7.2f -> %7.2f umol/L\n", at(b,0,"CITR"), fin(b,"CITR")))
  cat(sprintf("  urine        %7.3f -> %7.3f L/d\n", at(b,0,"URINE_LD"), fin(b,"URINE_LD")))

  # ---- D02  the baseline patient must look like a real SBS-IF patient
  cat("\n[D02] BASELINE PHENOTYPE vs published SBS-IF (80 cm jejunum + half colon)\n")
  cat(sprintf("  PN volume        %6.2f L/wk        (STEPS baseline ~12.9)\n", pn0))
  cat(sprintf("  stoma/stool out  %6.2f L/d\n", at(b,0,"OUT_LD")))
  cat(sprintf("  effluent [Na+]   %6.1f mmol/L      (observed ~90-100)\n", at(b,0,"STOMANA")))
  cat(sprintf("  citrulline       %6.1f umol/L      (SBS-IF ~10-20, normal 30-40)\n", at(b,0,"CITR")))
  cat(sprintf("  fat absorption   %6.1f %%           (observed 30-50%%)\n", 100*at(b,0,"FFATC")))
  cat(sprintf("  CHO absorption   %6.1f %%           (observed 70-90%%)\n", 100*at(b,0,"FCHOC")))
  cat(sprintf("  protein absorpt. %6.1f %%           (observed 60-80%%)\n", 100*at(b,0,"FPROC")))
  cat(sprintf("  hyperphagia      %6.2f x           (observed up to ~2x)\n",
      b$HYPERPH[1]))
  cat(sprintf("  colonic salvage  %6.0f kcal/d      (Nordgaard: up to ~1000)\n", at(b,0,"SCFAKC")))

  # ---- D03  STEPS: teduglutide vs placebo ---------------------------
  cat("\n[D03] STEPS REPRODUCTION — the endpoint is produced by the PROTOCOL\n")
  dT <- sim(st0, cp, ev_ted(wt0, 168), end = 168)
  dP <- sim(st0, cp, end = 168)
  cat(sprintf("  teduglutide  PN %6.2f -> %6.2f L/wk   change %+6.2f (obs -4.4)\n",
      pn0, at(dT,168,"PNVOLWK"), at(dT,168,"PNVOLWK")-pn0))
  cat(sprintf("  placebo      PN %6.2f -> %6.2f L/wk   change %+6.2f (obs -2.3)\n",
      pn0, at(dP,168,"PNVOLWK"), at(dP,168,"PNVOLWK")-pn0))
  cat(sprintf("  responder (>=20%% at wk20 AND wk24): teduglutide %s, placebo %s\n",
      ifelse(resp20(dT, pn0)==1,"YES","no"), ifelse(resp20(dP, pn0)==1,"YES","no")))
  cat(sprintf("  citrulline   %+6.2f umol/L on drug (STEPS: significant rise)\n",
      at(dT,168,"CITR") - at(dT,0,"CITR")))
  cat(sprintf("  villus (rel) %6.3f -> %6.3f\n", at(dT,0,"VILLC"), at(dT,168,"VILLC")))

  # ---- D04  responder RATE across a virtual population --------------
  # The two arms are NOT distinguished only by drug. Both are enrolled from
  # pre-trial care (suboptimal ORS) into a protocol that standardises diet
  # and oral rehydration AND applies the urine-triggered weaning algorithm.
  # That shared protocol effect is where the placebo response comes from;
  # the drug is added on top of it. Between-patient variability is placed
  # on absorptive capacity, colonic capacity, how aggressively the centre
  # weans, how much ORS optimisation is actually achieved, and GLP-2
  # sensitivity (which absorbs adherence and neutralising-antibody effects).
  cat("\n[D04] RESPONDER RATE across a virtual population\n")
  cat("      both arms enrolled from pre-trial care into the STEPS protocol\n")
  set.seed(20260727)
  n <- 120
  sbls  <- round(runif(n, 25, 170))
  colf  <- sample(c(0, 0.5, 1.0), n, replace = TRUE, prob = c(0.40, 0.38, 0.22))
  ilen  <- ifelse(colf > 0.6, round(runif(n, 0, 40)), 0)
  ilen  <- pmin(ilen, round(0.5 * sbls))          # ileum <= half the remnant
  fmaxj <- 1.70   * exp(rnorm(n, 0, 0.18))
  kcolm <- 5.00   * exp(rnorm(n, 0, 0.25))
  ktap  <- 0.0060 * exp(rnorm(n, 0, 0.35))
  orspre<- runif(n, 0.10, 0.45)                   # pre-trial ORS adherence
  orsgain <- runif(n, 0.00, 0.45)                 # what the protocol achieves
  ec50g <- 0.010  * exp(rnorm(n, 0, 0.80))        # GLP-2 sensitivity + adherence
  stopwk<- ifelse(runif(n) < 0.12, sample(4:16, n, TRUE), 24)  # early withdrawal
  rT <- rP <- dT20 <- dP20 <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    p <- list(SBL = sbls[i], COLONFRAC = colf[i], ILEUMLEN = ilen[i],
              ICV = as.numeric(colf[i] > 0.6), FMAXJ = fmaxj[i],
              KCOLMAX = kcolm[i], KTAPER = ktap[i], ORSFRAC = orspre[i])
    st <- burnin(p)
    if (st[["PNVOL"]] < 0.43) next        # STEPS required PN >= ~3 nights/wk
    wti <- st[["LEAN"]] + st[["FATM"]]
    # modifyList, NOT c(): c(p, list(ORSFRAC = ...)) leaves TWO ORSFRAC
    # entries in the list and param() takes the first, so the protocol
    # effect silently never happened and the placebo arm looked inert.
    q <- modifyList(p, list(URBASE = st[["URINESM"]], WTTARGET = wti,
                            ORSFRAC = min(0.95, orspre[i] + orsgain[i]),
                            EC50G = ec50g[i]))
    bp <- st[["PNVOL"]] * 7
    a <- sim(st, q, ev_ted(wti, min(168, stopwk[i] * 7)), end = 168)
    z <- sim(st, q, end = 168)
    rT[i] <- resp20(a, bp); rP[i] <- resp20(z, bp)
    dT20[i] <- at(a, 168, "PNVOLWK") - bp
    dP20[i] <- at(z, 168, "PNVOLWK") - bp
  }
  ok <- !is.na(rT)
  cat(sprintf("  evaluable PN-dependent patients: %d of %d simulated\n", sum(ok), n))
  cat(sprintf("  responders  teduglutide %5.1f %%   (STEPS 63%%)\n", 100*mean(rT[ok])))
  cat(sprintf("  responders  placebo     %5.1f %%   (STEPS 30%%)\n", 100*mean(rP[ok])))
  cat(sprintf("  mean PN change  drug %+6.2f L/wk   placebo %+6.2f L/wk  (obs -4.4 / -2.3)\n",
      mean(dT20[ok]), mean(dP20[ok])))
  cat(sprintf("  drug-attributable difference %+6.2f L/wk  (obs -2.1)\n",
      mean(dT20[ok]) - mean(dP20[ok])))

  # ---- D05  the plain-water paradox as a SIGN CHANGE ----------------
  cat("\n[D05] THE PLAIN-WATER PARADOX — generated, not asserted\n")
  cat("      same 2.5 L/d of between-meal fluid, three compositions\n")
  cat("      DRIVE_M = meal stream, DRIVE_D = between-meal drink stream.\n")
  cat("      The sign change lives in DRIVE_D and nowhere else.\n")
  for (of in c(0.0, 0.2, 0.4, 0.7, 0.9)) {
    d <- sim(st0, c(cp, list(ORSFRAC = of, WEANON = 0)), end = 60)
    cat(sprintf("   ORS %4.1f : [Na]drink %5.1f  DRIVE_M %+6.1f  DRIVE_D %+7.1f %s  FRACJ %+6.3f  output %5.2f L/d\n",
        of, fin(d,"CNAD"), fin(d,"DRIVEM"), fin(d,"DRIVED"),
        ifelse(fin(d,"DRIVED") < 0, "SECRETES", "absorbs "),
        fin(d,"FRACJC"), fin(d,"OUT_LD")))
  }
  d0 <- sim(st0, c(cp, list(ORSFRAC = 0.0)), end = 168)
  d9 <- sim(st0, c(cp, list(ORSFRAC = 0.9)), end = 168)
  cat(sprintf("   PN at 24 wk: plain water %5.2f L/wk vs full ORS %5.2f L/wk (difference %+5.2f)\n",
      at(d0,168,"PNVOLWK"), at(d9,168,"PNVOLWK"),
      at(d9,168,"PNVOLWK") - at(d0,168,"PNVOLWK")))

  # ---- D06  the colon is worth an anatomy ---------------------------
  cat("\n[D06] THE COLON AS A DIGESTIVE ORGAN — matched 80 cm of jejunum\n")
  for (cf in c(0, 0.5, 1.0)) {
    p <- list(SBL = 80, COLONFRAC = cf, ILEUMLEN = 0, ICV = as.numeric(cf > 0.9))
    s <- burnin(p)
    dd <- sim(s, c(p, list(URBASE = s[["URINESM"]],
                           WTTARGET = s[["LEAN"]] + s[["FATM"]])), end = 30)
    cat(sprintf("   colon %4.2f : PN %6.2f L/wk  output %5.2f L/d  SCFA %5.0f kcal/d  fat abs %4.1f %%\n",
        cf, s[["PNVOL"]]*7, fin(dd,"OUT_LD"), fin(dd,"SCFAKC"), 100*fin(dd,"FFATC")))
  }
  cat("      colonic salvage is largest where malabsorption is WORST, because\n")
  cat("      the substrate IS the unabsorbed carbohydrate:\n")
  for (cfg in list(c(40, 1.0, 0.50), c(40, 1.0, 0.68), c(80, 1.0, 0.68))) {
    p <- list(SBL = cfg[1], COLONFRAC = cfg[2], ILEUMLEN = 0, ICV = 1,
              FCHO_DIET = cfg[3], FFAT_DIET = 0.83 - cfg[3], FPRO_DIET = 0.17,
              FIBERKCAL = 300)
    sx <- burnin(p)
    dd <- sim(sx, c(p, list(URBASE = sx[["URINESM"]],
                            WTTARGET = sx[["LEAN"]] + sx[["FATM"]])), end = 30)
    cat(sprintf("   SBL %3.0f cm, full colon, CHO %2.0f%% of diet : unabsorbed CHO %5.0f kcal/d -> SCFA salvage %5.0f kcal/d\n",
        cfg[1], 100*cfg[3],
        dd$ABSKC[1]*0 + (dd$HYPERPH[1]*2400*cfg[3]*(1-dd$FCHOC[1]) + 300),
        fin(dd,"SCFAKC")))
  }

  # ---- D07  adaptation is a PRODUCT: nil by mouth abolishes the drug -
  cat("\n[D07] NUTRIENT GATING — and a partial refutation of the design\n")
  cat("      Expectation going in: adaptation is a PRODUCT, so nil-by-mouth\n")
  cat("      should abolish the drug effect entirely. It does not, and the\n")
  cat("      reason is informative. GLP-2 has TWO arms here and only one is\n")
  cat("      nutrient-gated: the trophic arm (villus, crypt, transporters)\n")
  cat("      goes to zero as designed, but the motility/ENS arm (contact time\n")
  cat("      through the ileal brake) needs no luminal nutrient and survives.\n")
  cat("      Gating ATTENUATES the drug, it does not abolish it - and that is\n")
  cat("      a testable structural prediction, not a fitted result.\n")
  for (ef in c(1.0, 0.5, 0.02)) {
    s <- burnin(c(base, list(ENTFRAC = ef)))
    q <- c(base, list(ENTFRAC = ef, URBASE = s[["URINESM"]],
                      WTTARGET = s[["LEAN"]] + s[["FATM"]]))
    dd <- sim(s, q, ev_ted(s[["LEAN"]]+s[["FATM"]], 168), end = 168)
    dn <- sim(s, q, end = 168)
    cat(sprintf("   enteral %4.2f : trophic %5.3f  villus %5.3f -> %5.3f  drug-attributable PN change %+6.2f L/wk\n",
        ef, fin(dd,"TROPHC"), at(dd,0,"VILLC"), at(dd,168,"VILLC"),
        at(dd,168,"PNVOLWK") - at(dn,168,"PNVOLWK")))
  }

  # ---- D08  the octreotide trade-off --------------------------------
  cat("\n[D08] OCTREOTIDE TRADE-OFF — output now, adaptation later\n")
  do <- sim(st0, cp, ev_oct(730), end = 730)
  dn <- sim(st0, cp, end = 730)
  for (tt in c(28, 168, 365, 730)) {
    cat(sprintf("   day %4d : output %5.2f vs %5.2f L/d   villus %5.3f vs %5.3f   PN %6.2f vs %6.2f L/wk\n",
        tt, at(do,tt,"OUT_LD"), at(dn,tt,"OUT_LD"),
        at(do,tt,"VILLC"), at(dn,tt,"VILLC"),
        at(do,tt,"PNVOLWK"), at(dn,tt,"PNVOLWK")))
  }

  # ---- D09  cholestyramine reverses sign across the 100 cm rule -----
  cat("\n[D09] THE 100 cm RULE — cholestyramine helps, then harms\n")
  cat("      total remnant fixed at 90 cm with half a colon; only the\n")
  cat("      PROPORTION that is ileum changes, so length is not confounded\n")
  for (il in c(90, 60, 25, 0)) {
    p <- list(SBL = 90, ILEUMLEN = il, COLONFRAC = 0.5, ICV = 0)
    s <- burnin(p)
    q <- c(p, list(URBASE = s[["URINESM"]],
                   WTTARGET = s[["LEAN"]] + s[["FATM"]]))
    dc <- sim(s, q, ev_chol(168), end = 168)
    dz <- sim(s, q, end = 168)
    cat(sprintf("   ileum %4d cm : pool %5.2f g  free BA %5.1f mmol/d  secretagogue %5.1f  PN %5.2f L/wk | on cholestyramine: output %+6.3f L/d  fat abs %+5.1f %%  absorbed energy %+5.0f kcal/d  PN %+6.2f L/wk\n",
        il, at(dz,168,"BAPOOLC"), at(dz,168,"BACOLC"), at(dz,168,"SECLOADC"),
        at(dz,168,"PNVOLWK"),
        at(dc,168,"OUT_LD") - at(dz,168,"OUT_LD"),
        100*(at(dc,168,"FFATC") - at(dz,168,"FFATC")),
        at(dc,168,"ABSKC") - at(dz,168,"ABSKC"),
        at(dc,168,"PNVOLWK") - at(dz,168,"PNVOLWK")))
  }

  # ---- D10  the placebo response is a property of the ALGORITHM -----
  cat("\n[D10] WHERE THE PLACEBO RESPONSE COMES FROM\n")
  cat("      the trial run-in standardises ORS and diet in BOTH arms\n")
  s <- burnin(c(base, list(ORSFRAC = 0.25)))          # pre-trial, suboptimal care
  q <- c(base, list(URBASE = s[["URINESM"]],
                    WTTARGET = s[["LEAN"]] + s[["FATM"]]))
  bp <- s[["PNVOL"]]*7
  dpo <- sim(s, c(q, list(ORSFRAC = 0.60)), end = 168)             # placebo + protocol
  dto <- sim(s, c(q, list(ORSFRAC = 0.60)), ev_ted(s[["LEAN"]]+s[["FATM"]],168), end=168)
  dpn <- sim(s, c(q, list(ORSFRAC = 0.25, WEANON = 0)), end = 168) # no protocol at all
  cat(sprintf("   placebo WITH protocolised care   PN %6.2f -> %6.2f L/wk (%+.1f%%) responder %s\n",
      bp, at(dpo,168,"PNVOLWK"), 100*(at(dpo,168,"PNVOLWK")-bp)/bp,
      ifelse(resp20(dpo,bp)==1,"YES","no")))
  cat(sprintf("   teduglutide WITH protocolised care PN %6.2f -> %6.2f L/wk (%+.1f%%) responder %s\n",
      bp, at(dto,168,"PNVOLWK"), 100*(at(dto,168,"PNVOLWK")-bp)/bp,
      ifelse(resp20(dto,bp)==1,"YES","no")))
  cat(sprintf("   no protocol, no drug              PN %6.2f -> %6.2f L/wk (%+.1f%%)\n",
      bp, at(dpn,168,"PNVOLWK"), 100*(at(dpn,168,"PNVOLWK")-bp)/bp))

  # ---- D11  native GLP-2 vs teduglutide -----------------------------
  cat("\n[D11] WHY THE Ala2->Gly SUBSTITUTION IS THE DRUG\n")
  dt <- sim(st0, cp, ev_ted(wt0, 168), end = 168)
  n1 <- sim(st0, cp, ev_native(wt0, 168, 1), end = 168)
  n2 <- sim(st0, cp, ev_native(wt0, 168, 2), end = 168)
  n8 <- sim(st0, cp, ev_native(wt0, 168, 2, dose_mgkg = 0.40), end = 168)
  cat(sprintf("   teduglutide 0.05 mg/kg od       : mean GLP2R occupancy %5.3f  PN %+6.2f L/wk\n",
      mean(dt$OCCGC), at(dt,168,"PNVOLWK")-pn0))
  cat(sprintf("   native GLP-2 equimolar od       : mean occupancy %5.3f  PN %+6.2f L/wk\n",
      mean(n1$OCCGC), at(n1,168,"PNVOLWK")-pn0))
  cat(sprintf("   native GLP-2 equimolar bd       : mean occupancy %5.3f  PN %+6.2f L/wk\n",
      mean(n2$OCCGC), at(n2,168,"PNVOLWK")-pn0))
  cat(sprintf("   native GLP-2 8x dose bd         : mean occupancy %5.3f  PN %+6.2f L/wk\n",
      mean(n8$OCCGC), at(n8,168,"PNVOLWK")-pn0))

  # ---- D12  long-acting analogues at matched exposure ---------------
  cat("\n[D12] LONG-ACTING ANALOGUES — profile shape at weekly dosing\n")
  da <- sim(st0, cp, ev_apra(168), end = 168)
  dg <- sim(st0, cp, ev_glep(168, ii = 7), end = 168)
  dg2<- sim(st0, cp, ev_glep(168, ii = 3.5), end = 168)
  for (nm in list(c("teduglutide od","dt"), c("apraglutide qw","da"),
                  c("glepaglutide qw","dg"), c("glepaglutide 2x/wk","dg2"))) {
    d <- get(nm[2])
    cat(sprintf("   %-20s mean occ %5.3f  trough/peak %5.3f  PN %+6.2f L/wk  citrulline %+5.2f\n",
        nm[1], mean(d$OCCGC),
        min(tail(d$OCCGC, 30))/dmax_r(max(tail(d$OCCGC,30))),
        at(d,168,"PNVOLWK")-pn0, at(d,168,"CITR")-at(d,0,"CITR")))
  }

  # ---- D13  REFUTED: the anatomical floor ---------------------------
  cat("\n[D13] SELF-REFUTING RESULT #1 — the anatomical floor no drug crosses\n")
  cat("      expectation going in: enough drug reaches enteral autonomy at any length\n")
  for (sb in c(40, 60, 90, 130, 180, 250, 320)) {
    p <- list(SBL = sb, COLONFRAC = 0, ILEUMLEN = 0, ICV = 0)
    s <- burnin(p)
    q <- c(p, list(URBASE = s[["URINESM"]], WTTARGET = s[["LEAN"]]+s[["FATM"]]))
    d <- sim(s, q, ev_ted(s[["LEAN"]]+s[["FATM"]], 730), end = 730)
    cat(sprintf("   end-jejunostomy %3d cm : PN %6.2f -> %6.2f L/wk at 2 y   autonomy %s\n",
        sb, s[["PNVOL"]]*7, at(d,730,"PNVOLWK"),
        ifelse(at(d,730,"PNVOLWK") < 0.35, "YES", "no")))
  }

  # ---- D14  REFUTED: colon in continuity is not free ----------------
  cat("\n[D14] SELF-REFUTING RESULT #2 — the colon charges rent\n")
  for (cf in c(0, 1.0)) {
    p <- list(SBL = 80, COLONFRAC = cf, ILEUMLEN = 0, ICV = 0)
    s <- burnin(p)
    q <- c(p, list(URBASE = s[["URINESM"]], WTTARGET = s[["LEAN"]]+s[["FATM"]]))
    d <- sim(s, q, end = 730)
    cat(sprintf("   colon %4.2f : PN %6.2f L/wk | urinary oxalate %5.1f mg/d  stones %5.3f  D-lactate %5.3f\n",
        cf, at(d,730,"PNVOLWK"), at(d,730,"UOXC"), at(d,730,"STONEC"),
        at(d,730,"DLACTC")))
  }

  # ---- D15  IFALD: the lipid emulsion IS the dose -------------------
  cat("\n[D15] IFALD — changing the oil changes the disease\n")
  d1 <- sim(st0, c(cp, list(PHYTOCONC = 350, OMEGA3 = 0)), end = 730)
  d2 <- sim(st0, c(cp, list(PHYTOCONC = 100, OMEGA3 = 1)), end = 730)
  d3 <- sim(st0, c(cp, list(PHYTOCONC = 0,   OMEGA3 = 1)), end = 730)
  for (nm in list(c("soybean ILE","d1"), c("SMOF ILE","d2"), c("fish-oil ILE","d3"))) {
    d <- get(nm[2])
    cat(sprintf("   %-14s : phytosterol %6.2f  bilirubin %5.2f mg/dL  fibrosis %5.3f\n",
        nm[1], at(d,730,"PHYTOC"), at(d,730,"BILIC"), at(d,730,"FIBC")))
  }
  dr <- sim(laststate(d1), c(cp, list(PHYTOCONC = 0, OMEGA3 = 1)), end = 365)
  cat(sprintf("   rescue: switch to fish-oil after 2 y of soybean -> bilirubin %5.2f -> %5.2f mg/dL by 1 y\n",
      at(dr,0,"BILIC"), at(dr,365,"BILIC")))

  # ---- D16  central venous access as a consumable ------------------
  cat("\n[D16] ACCESS IS A CONSUMABLE — and the drug that cuts PN days saves veins\n")
  d_std <- sim(st0, cp, end = 730)
  d_lock<- sim(st0, c(cp, list(LOCKFAC = 0.25, TECHFAC = 0.7)), end = 730)
  d_ted <- sim(st0, cp, ev_ted(wt0, 730), end = 730)
  for (nm in list(c("standard care","d_std"), c("taurolidine lock + team","d_lock"),
                  c("teduglutide","d_ted"))) {
    d <- get(nm[2])
    cat(sprintf("   %-24s nights/wk %4.2f  catheter-days %5.0f  CRBSI %5.2f /y  cumulative %5.2f  veins left %5.2f  survival %5.1f %%\n",
        nm[1], at(d,730,"PNNIGHTS"), at(d,730,"CATHD"), at(d,730,"CRBSIYR"),
        at(d,730,"CRBSICUM"), at(d,730,"VEINSC"), at(d,730,"SURVP")))
  }
  # ---- D17  the results that went against the model design ----------
  cat("\n[D17] NEGATIVE AND SELF-REFUTING RESULTS — reported, not removed\n")
  cat("\n  (a) THE CHOLESTYRAMINE TRAP COMES OUT BACKWARDS.\n")
  cat("      The textbook rule is that a bile acid sequestrant helps when\n")
  cat("      < ~100 cm of ileum was resected and HARMS when more was, because\n")
  cat("      the pool is already depleted. D09 reproduces the pool gradient\n")
  cat("      (2.3 g with 90 cm of ileum left, 0.8 g with none) but the fat\n")
  cat("      penalty comes out LARGEST where the pool is INTACT (-16.5%) and\n")
  cat("      almost absent where it is depleted (-0.5%) - for the unavoidable\n")
  cat("      reason that a depleted pool has little left to sequester. The\n")
  cat("      clinical rule must therefore rest on something this model does\n")
  cat("      not contain: fat-soluble vitamin and drug binding by the resin,\n")
  cat("      and the fact that in extensive resection the residual benefit is\n")
  cat("      zero so ANY cost is unacceptable. The model reproduces 'no useful\n")
  cat("      benefit', not 'active harm', and that gap is a real limitation.\n")
  cat("\n  (b) THE POPULATION EFFECT IS OVER-PREDICTED EVEN THOUGH THE\n")
  cat("      INDIVIDUAL ANCHOR MATCHES. D03 lands on the STEPS mean almost\n")
  cat("      exactly (-4.45 vs -4.4 L/wk) because it was calibrated there,\n")
  cat("      but the virtual population in D04 gives a larger DRUG-MINUS-\n")
  cat("      PLACEBO difference than the trial did. Matching a mean is not\n")
  cat("      matching a distribution: the model is missing sources of\n")
  cat("      non-response (dose interruption, intercurrent illness, weaning\n")
  cat("      caution that is not urine-driven) and its responder rate should\n")
  cat("      not be quoted.\n")
  cat("\n  (c) NUTRIENT GATING ATTENUATES BUT DOES NOT ABOLISH (see D07).\n")
  cat("      The product structure was built to make a nil-by-mouth patient\n")
  cat("      unresponsive to a GLP-2 analogue. It does not, because only the\n")
  cat("      trophic arm is gated and the motility arm is not. Whether the\n")
  cat("      real drug retains a transit effect in a fully fasted gut is an\n")
  cat("      open experimental question this model now poses explicitly.\n")
  cat("\n  (d) OCTREOTIDE LOOKS BETTER THAN IT IS. Even after halving its\n")
  cat("      secretory effect and adding the pancreatic-lipase penalty, D08\n")
  cat("      still shows a larger sustained PN reduction than the clinical\n")
  cat("      literature supports. Tachyphylaxis, gallstone formation and\n")
  cat("      poor tolerability are not in the model, and they are most of the\n")
  cat("      reason octreotide is a last resort rather than a first line.\n")
  cat("\n  (e) BASELINE DRIFT IS SMALL BUT NOT ZERO (D01, ~0.5%/year in PN\n")
  cat("      volume). The burn-in converges the physiological states but the\n")
  cat("      slow structural-remodelling state has a ~250-day time constant,\n")
  cat("      so a 900-day run-in leaves a residual trend. Effects smaller\n")
  cat("      than ~0.5%/year cannot be resolved by this model.\n")
  cat("\n  (g) THE SIBO LIMB IS ALIVE BUT ENDPOINT-SILENT. Losing the\n")
  cat("      ileocaecal valve raises the bacterial load from ~0 to 0.62 and\n")
  cat("      a 14-day rifaximin course clears it to 0.35 with a measurable\n")
  cat("      villus recovery - but the load regrows to baseline by ~day 60\n")
  cat("      and the 24-week PN volume is unchanged. That is a faithful\n")
  cat("      picture of a transient antibiotic course, but it also means the\n")
  cat("      model contains no mechanism by which SIBO changes the endpoint\n")
  cat("      the trials measure. Bacterial bile-acid deconjugation is drawn\n")
  cat("      on the map and is NOT in the equations; wiring it in is the\n")
  cat("      obvious next step and would probably change this result.\n")
  cat("\n  (f) EFFLUENT SODIUM RUNS LOW. D02 gives ~84 mmol/L against an\n")
  cat("      observed 90-100. The two-stream construction dilutes the lumen\n")
  cat("      with between-meal fluid in a way a real intermittent drinker\n")
  cat("      does not, and no attempt was made to fit this away.\n")
  cat("\n=====================================================================\n")
  invisible(TRUE)
}

# =====================================================================
# ENTRY POINT
# =====================================================================
# Running `Rscript sbs_mrgsolve_model.R` executes every scenario and every
# diagnostic. Set SBS_NORUN=1 to source the model only (used by the Shiny app).
if (!identical(Sys.getenv("SBS_NORUN"), "1") &&
    (!interactive() || identical(Sys.getenv("SBS_RUN"), "1"))) {
  res <- run_scenarios()
  cat(sprintf("\n%d scenarios simulated.\n", length(res$scen)))
  for (nm in names(res$scen)) {
    d <- res$scen[[nm]]$sim
    cat(sprintf("  %-52s PN %6.2f -> %6.2f L/wk | out %5.2f L/d | wt %5.1f kg\n",
        nm, d$PNVOLWK[1], d$PNVOLWK[nrow(d)],
        d$OUT_LD[nrow(d)], d$WT[nrow(d)]))
  }
  run_diagnostics()
}
