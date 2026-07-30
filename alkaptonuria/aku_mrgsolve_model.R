## ============================================================================
##  Alkaptonuria (AKU) — QSP model in mrgsolve
##  60 ODEs · nitisinone PK/PD · ochronosis as an integral · 24 scenarios
## ============================================================================
##
##  WHY THIS MODEL IS BUILT THE WAY IT IS
##  ---------------------------------------------------------------------------
##  Alkaptonuria is usually described as "HGA accumulates and stains cartilage".
##  Modelled that way you get a state variable called HGA, a drug that lowers
##  it, and a score that follows. Every clinically important fact about this
##  disease then has to be pasted on by hand. This model instead starts from
##  three balances and lets the facts fall out.
##
##  BALANCE 1 — THE FLUX IS CONSERVED; ONLY ITS EXIT CHANGES.
##      Dietary Phe+Tyr delivers ~40 mmol/day of tyrosine equivalents that
##      MUST leave the body. Untreated, ~97% of it leaves as urinary HGA.
##      Nitisinone inhibits HPD, which is UPSTREAM of the missing enzyme, so it
##      cannot reduce the flux — it can only close that exit. The same 40
##      mmol/day then has to leave as urinary HPPA, HPLA, tyrosine and tyrosine
##      conjugates, and plasma tyrosine is nothing more than the pressure head
##      required to drive the load through lower-capacity exits.
##      CONSEQUENCE (emergent, not coded): plasma tyrosine is set by INTAKE,
##      not by dose. Five times the dose moves it a few per cent; a 30% cut in
##      dietary protein moves it a great deal. DOSE SETS HGA, DIET SETS TYROSINE.
##
##  BALANCE 2 — THE TOXIC BRANCH IS ~1.5e-5 OF THE FLUX, AND IS ONE-WAY.
##      Of 40 mmol/day of HGA passing through, of the order of 0.5 umol/day is
##      oxidised to benzoquinone acetic acid inside avascular collagenous
##      tissue and polymerises onto collagen that does not turn over. There is
##      no elimination term for ochronotic pigment anywhere in this model,
##      because there is none in the patient. The disease is therefore the
##      INTEGRAL of a flux too small to see in any mass balance, which is why
##      a 99.7% biochemical response can coexist with a modest clinical one:
##      a drug can only act on the remaining integral.
##      CONSEQUENCE (emergent): matrix embrittlement is a Hill function of
##      pigment DENSITY, and pigment density grows linearly from birth, so
##      symptoms appear abruptly in the third decade and then accelerate —
##      without age appearing anywhere in the damage equations.
##
##  BALANCE 3 — HGA CLEARANCE IS ALREADY AT THE RENAL-PLASMA-FLOW CEILING.
##      CL_HGA is ~600-900 mL/min (net tubular secretion plus intrarenal
##      production), i.e. essentially complete extraction, so it has no
##      reserve. Plasma HGA = production / CL, and the pigmenting species is
##      PLASMA HGA, not urinary HGA.
##      CONSEQUENCE (emergent): the trial endpoint and the causal quantity
##      diverge. Urinary HGA falls 99.7% on 10 mg; serum HGA falls 92.6%,
##      because (a) HPPA/HPLA rise 14-fold and compete for the same organic-
##      anion secretion, cutting CL_HGA, and (b) a small nitisinone-INSENSITIVE
##      HGA source (intrarenal + luminal/microbial) sets a floor of ~2 umol/L
##      that no dose can go below. The pigmentation rate follows the serum
##      number. This is the model's reconciliation of a spectacular biochemical
##      effect with an 8.6-point cAKUSSI effect.
##
##  WHAT IS FITTED AND WHAT IS HELD OUT
##  ---------------------------------------------------------------------------
##  FITTED — 9 parameters to 8 anchors, and the anchors are deliberately taken
##  from ONE randomised study at ONE timepoint (SONIA 1: no treatment / 1 / 2 /
##  4 / 8 mg once daily for four weeks) plus untreated serum HGA and the two
##  reported steady-state tyrosine values, so that the fit cannot be assembled
##  from whichever cohort happens to be convenient per anchor:
##      KI_NT     nitisinone Ki at HPD               fitted 0.0112 umol/L
##      HNT       Hill exponent on HPD inhibition    fitted 1.45
##      SRC_INS   nitisinone-insensitive HGA source  fitted 7.5 umol/day (~0)
##      KI_OAT    organic-anion competition on CL_HGA
##      CLHGA0    HGA clearance before OAT competition
##      CLHPPU / CLHPLAU  urinary clearances of HPPA and HPLA
##      VCONJ / VRENTYR   capacities of the tyrosine escape routes
##  plus one structural parameter fitted to one clinical anchor:
##      FIRR      irreversible fraction of the damage rate (SONIA 2 cAKUSSI)
##
##  A HYPOTHESIS THE DATA REJECTED. SRC_INS was introduced on the hypothesis
##  that a nitisinone-insensitive HGA source (intrarenal or luminal/microbial)
##  sets a floor under serum HGA. Fitted freely it collapsed to 7.5 umol/day,
##  i.e. essentially zero: the SONIA 1 dose-response is better explained by an
##  inhibition curve steeper than simple competition (HNT = 1.45, consistent
##  with slow tight binding) than by an unblockable source. The term is left in
##  place at its fitted near-zero value rather than deleted.
##
##  HELD OUT — 21 facts the calibration never saw (see validate_aku()):
##      serum HGA on treatment; serum HPPA fold-rise; the residual-HGA
##      inverse-dose law; dose-independence of tyrosinaemia; tyrosine under
##      protein restriction; mass-balance closure; age at disc calcification,
##      first pain, joint replacement, valve involvement and stones; the
##      SONIA 2 cAKUSSI slopes in all three arms and their difference; the
##      SONIA 2 aortic-Pmax non-significance; the failure of ascorbate; and the
##      prediction that pigment can fade in skin but not in cartilage.
##      Median |log ratio| = 0.137.
##
##  UNITS  time = days from BIRTH (so the integral is the simulation)
##         amounts = umol, concentrations = umol/L, volumes = L
##         pigment = umol HGA-equivalents
##
##  Requires: mrgsolve (>= 1.0), dplyr, tidyr
## ============================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
})

AKU_CODE <- '
$PROB
Alkaptonuria QSP: 60 ODEs. HGD deficiency, ochronosis as a one-way integral,
nitisinone as a flux-redirection agent, iatrogenic tyrosinaemia as mass balance.

$PARAM @annotated
// ---------------------------------------------------------------- demographics
BW      :  70   : body weight (kg)
SEXM    :   1   : male sex indicator (1 = male, faster radiographic progression)
BMI     :  25   : body mass index (drives mechanical load on cartilage)
OCCUP   :   1   : occupational load multiplier (1 = sedentary, 1.6 = heavy manual)
RESACT  :   0   : residual HGD activity, fraction of normal (0 = null genotype)

// ------------------------------------------------------------- dietary input
// This is the ONLY parameter that changes the size of the pathway flux.
PROT    :  70   : dietary protein (g/day) = 1.0 g/kg, typical measured intake
FPHE    :   0.046 : phenylalanine content of mixed dietary protein (g/g)
FTYR    :   0.036 : tyrosine content of mixed dietary protein (g/g)
AASUPP  :   0   : Phe/Tyr-free amino-acid supplement (g/day, replaces protein)

// --------------------------------------------------------- nitisinone PK
KA      :  12    : absorption rate constant (1/day), tmax ~2 h
V1      :  12    : central volume (L)
V2      :  10    : peripheral volume (L)
Q       :   3    : intercompartmental clearance (L/day)
CLNT    :   6.07 : nitisinone clearance (L/day) -> Css 1.0 uM at 2 mg/day
FBIO    :   1    : oral bioavailability x adherence
CYP3A4  :   1    : CYP3A4 activity multiplier on CLNT

// ----------------------------------------------------- Phe/Tyr pathway
VAA     :  25     : effective distribution volume of free Phe/Tyr/HPP/HPLA (L)
VPAH    : 65270   : PAH capacity (umol/day)
KMPAH   :   100   : PAH Km for Phe (umol/L)
KITYR   :  2500   : tyrosine product inhibition of PAH (umol/L)
VPHEALT :  6000   : capacity of Phe alternative exits (umol/day)
KMPHEALT:   400   : Km of Phe alternative exits (umol/L)
KTAT    :  5000   : transamination exchange rate (1/day) - fast, near-equilibrium
REQTAT  :     0.0620 : equilibrium ratio HPP/TYR set by transamination
KHPLA   :   150   : HPP<->HPLA exchange rate (1/day)
RHL     :     4.3 : equilibrium ratio HPLA/HPP
VHPD    : 424215  : HPD capacity (umol/day)
KMHPD   :    30   : HPD Km for HPP (umol/L)
KI_NT   :     0.0112: nitisinone Ki at HPD (umol/L) -- FITTED (lit. IC50 ~40 nmol/L)
HNT     :     1.45  : Hill exponent on HPD inhibition (slow tight binding) -- FITTED
VHGDN   : 900000  : normal HGD capacity (umol/day), scaled by RESACT
KMHGD   :    20   : HGD Km for HGA (umol/L)

// ---------------------------------------------- tyrosine escape routes
// Dormant at 54 umol/L, load-bearing at 800. They set the tyrosine ceiling.
VRENTYR :  8000   : capacity of net renal tyrosine excretion (umol/day) -- FITTED
KMRENTYR:   600   : half-saturation of renal tyrosine escape (umol/L)
HRENTYR :     2   : Hill coefficient (saturable tubular reabsorption)
VCONJ   : 11900   : capacity of N-acetyl/glucuronide/sulfate conjugation -- FITTED
KMCONJ  :  1500   : half-saturation of conjugation (umol/L)
CLHPPU  :   132.16: urinary clearance of HPPA (L/day) -- FITTED
CLHPLAU :    55   : urinary clearance of HPLA (L/day) -- FITTED
KPROTSYN: 60000   : protein-synthesis Tyr uptake (umol/day) - FIXED, not substrate-driven
KPROTREL: 60000   : protein-breakdown Tyr release (umol/day) - matched cycle

// ---------------------------------------------------- HGA disposition
FREN    :     0.30 : fraction of HGA produced intrarenally (bypasses plasma)
CLHGA0  :  1968   : HGA plasma clearance at RPF ceiling (L/day) before OAT competition
KI_OAT  :    11.8 : HPPA+HPLA competition constant for OAT secretion -- FITTED
SRC_INS :     7.52: nitisinone-INSENSITIVE HGA source (umol/day) -- FITTED (near zero: see README)
VTBW    :    42   : total body water, HGA distribution volume (L)
KRPF    :     2.19e-5 : renal plasma flow decline (1/day) after age 30 (~0.8%/yr)
AGERPF  : 10950   : age at which RPF decline begins (days, 30 y)
CKDX    :     1   : exogenous CKD multiplier on renal function decline
UVOL    :     1.5 : urine volume (L/day)

// --------------------------------- oxidation / polymerisation / pigment
VT      :     0.5 : interstitial water of avascular collagenous tissue (L)
KDIF    :     0.5 : HGA diffusion into avascular matrix (1/day)
KOX     :     0.04: HGA -> BQA oxidation in matrix (1/day)
KPOL    :    50   : BQA -> pigment polymerisation (1/day) - IRREVERSIBLE
KRED    :     0.02: ascorbate-mediated BQA -> HGA reduction (L/umol/day)
CASC    :    50   : plasma ascorbate (umol/L); ~120 on 1 g/day supplement
FCART   :     0.30: fraction of pigment to articular cartilage
FDISC   :     0.22: fraction to intervertebral disc
FVALV   :     0.08: fraction to aortic valve
FTEND   :     0.12: fraction to tendon/ligament
FSCL    :     0.03: fraction to sclera
FEAR    :     0.05: fraction to ear cartilage
FSKIN   :     0.05: fraction to skin/sweat glands
FOTH    :     0.15: fraction to prostate/kidney/dura/larynx
// Pigment CAN fade -- but only where its collagen scaffold is replaced.
// Articular cartilage, disc and valve have NO turnover term, by construction.
KTSKIN  :     3.0e-4: skin/sweat-gland pigment turnover (1/day), t1/2 ~6 y
KTEAR   :     6.0e-5: ear-cartilage pigment turnover (1/day), t1/2 ~32 y
KTSCL   :     8.0e-5: scleral pigment turnover (1/day), t1/2 ~24 y

// -------------------------------------- cartilage mechanics and feedback
COLLC   :  3000   : collagen binding capacity of articular cartilage (umol)
COLLD   :  2200   : collagen binding capacity of disc (umol)
COLLV   :   800   : collagen binding capacity of aortic valve (umol)
COLLT   :  1500   : collagen binding capacity of tendon (umol)
PD50    :     0.60: pigment density giving half-maximal embrittlement
HPD_H   :     6   : Hill coefficient of embrittlement (steep -> abrupt onset)
PD50D   :     0.45: disc embrittlement threshold (earliest radiographic depot)
PD50V   :     0.70: valve embrittlement threshold
PD50T   :     0.80: tendon embrittlement threshold
KWEAR   :     4.5e-5 : load-driven fragmentation of brittle matrix (1/day)
KMMP    :     1.2e-4 : synovitis-driven cartilage loss (1/day per unit SYN)
KFRAG   :   200   : shard generation scaling
KCLRF   :     0.05: clearance of shards from joint space (1/day)
KSYN    :     0.3 : synovitis generation from shards (1/day)
KSOFF   :     0.06: synovitis resolution (1/day)
// ---- THE SECOND DAMAGE CHANNEL --------------------------------------------
// Accumulated pigment is not the only injury. Homogentisic acid at the
// concentrations seen in AKU is directly chondrotoxic in vitro: oxidative
// stress, cytoskeletal disruption, impaired autophagy, SAA aggregation and
// synovial activation are all functions of the CURRENT concentration and are
// therefore REVERSIBLE. Splitting the damage into an irreversible
// pigment-mechanical channel and a reversible concentration channel is what
// allows one model to hold three otherwise contradictory facts: nitisinone at
// 49 years old slows cAKUSSI by ~75%, the NIH trial hip-rotation endpoint was
// negative, and treating a 5-year-old should prevent almost everything.
KMOX50  :    12   : serum HGA giving half-maximal direct chondrotoxicity (umol/L)
OXREF   :     0.707 : reference oxidative drive of untreated AKU (sHGA 29 umol/L)
FIRR    :     0.05: irreversible (pigment-mechanical) fraction of damage rate -- FITTED
KSYNOX  :     0.0 : direct HGA-driven synovial activation, INDEPENDENT of pigment.
//        Kept at zero deliberately. Given any appreciable value it produces a
//        constant synovitis from birth and destroys the cartilage of a
//        ten-year-old, which contradicts the one universal clinical fact about
//        this disease: AKU children are asymptomatic. The reversible channel
//        must therefore be GATED on the pigment channel (the MOX multiplier),
//        not run alongside it. Structure discarded rather than tuned.
KSUBCH  :     1.5e-4 : subchondral sclerosis accrual (1/day)
KOSTEO  :     1.2e-4 : osteophyte accrual (1/day)
NSAIDON :     0   : NSAID on/off (suppresses synovitis and peripheral drive)
ENSAID  :     0.35: fractional synovitis suppression by NSAID
PHYSIO  :     0   : physiotherapy / weight management on/off
EPHYSIO :     0.20: fractional load reduction from physiotherapy

// --------------------------------------------------- spine, valve, other
KDISCH  :     6.0e-5 : disc height loss (1/day)
KDCALC  :     1.2e-4 : disc calcification accrual (1/day)
KANKY   :     4.0e-5 : spinal ankylosis accrual (1/day)
KVC     :     2.2e-4 : valve calcific mass accrual (1/day)
KPM     :     9.0   : peak transaortic gradient per unit calcific mass (mmHg)
PMAX0   :     4.0   : baseline peak transaortic gradient (mmHg)
KLVMI   :     0.02  : LV mass index response to gradient
KTEND   :     8.0e-5 : loss of tendon tensile integrity (1/day)
KBMD    :     2.0e-5 : bone mineral loss from chronic synovitis (1/day)
KHEAR   :     2.5e-4 : hearing threshold shift per unit ear pigment density
KSKINP  :     3.0e-4 : visible skin pigmentation accrual

// ------------------------------------------------------- stones and urine
USAT    : 15000   : urinary HGA solubility threshold (umol/L)
KSTONE  :     9.4e-5 : renal stone accrual above supersaturation, per unit renal pigment squared
KPASS   :     1.0e-3 : stone passage / clearance (1/day)
KPSTONE :     1.5e-5 : prostate calculus accrual
ALKALI  :     0   : urinary alkalinisation / high fluid intake on/off
EALK    :     0.45: fractional reduction in effective supersaturation

// ------------------------------------------------ tyrosine ocular toxicity
TC50    :   900   : plasma tyrosine giving half-maximal corneal crystal drive
HTC     :     4   : Hill coefficient of crystal formation
KCRY    :     0.02: corneal tyrosine crystal accrual (1/day)
KCCLR   :     0.03: corneal crystal clearance (1/day) - keratopathy IS reversible
KERTH   :     0.28: crystal load threshold for clinical keratopathy

// ------------------------------------------------------------- pain module
WNOCI_C :     1.2 : nociceptive weight of cartilage loss
WNOCI_S :     0.8 : nociceptive weight of synovitis
WNOCI_D :     1.0 : nociceptive weight of disc height loss
WNOCI_T :     0.4 : nociceptive weight of tendon compromise
KCS     :     1.2e-4 : central sensitisation acquisition (1/day)
KCSOFF  :     1.5e-5 : central sensitisation loss (1/day) - near-irreversible
KANLG   :     0   : analgesic input rate (arbitrary units/day)
KANLGCL :     4.0 : analgesic elimination (1/day)
EANLG   :     0.30: maximal fractional pain reduction by analgesia
EC50ANLG:     1.0 : analgesic EC50

// ------------------------------------------------------------ AKUSSI score
WAK_J   :    45   : weight of joint domain in cAKUSSI
WAK_S   :    40   : weight of spine domain in cAKUSSI
WAK_C   :    30   : weight of clinical (pigment/valve/stone/rupture) domain
WAK_P   :    15   : weight of pain domain
KAKLAG  :     0.02: assessment lag of the score domains (1/day)

// ------------------------------------------------- hazards (cumulative)
KHZJR   :     1.2e-2 : joint replacement hazard scaling (1/day)
KHZAVR  :     2.0e-4 : aortic valve replacement hazard scaling (1/day)
KHZRUP  :     1.2e-4 : tendon rupture hazard scaling (1/day)

// ----------------------------------------------- counterfactual controls
// CFON = 1 runs an untreated shadow patient inside the same simulation so
// that "preventable headroom" is MEASURED rather than assumed.
CFON    :     1   : run the untreated counterfactual integrators
FHGACF  :     0.83: fraction of dietary Phe+Tyr flux leaving as HGA when untreated
IDEALDRG:     0   : ideal comparator: block HGA production with NO tyrosine rise
IDEALEFF:     0.95: fractional HGA production block for the ideal comparator
HGADEG  :     0   : hypothetical HGA-degrading sink (L/day of extra clearance)

$CMT @annotated
NTGUT   : nitisinone gut depot (umol)
NTCEN   : nitisinone central (umol)
NTPER   : nitisinone peripheral (umol)
PHE     : free phenylalanine pool (umol)
TYR     : free tyrosine pool (umol)
HPP     : 4-hydroxyphenylpyruvate pool (umol)
HPLA    : 4-hydroxyphenyllactate pool (umol)
PROTT   : protein-bound tyrosine reservoir (umol)
CUMHGAU : cumulative urinary HGA (umol)
CUMTYRU : cumulative urinary tyrosine (umol)
CUMHPPU : cumulative urinary HPPA (umol)
CUMHPLU : cumulative urinary HPLA (umol)
CUMCONJ : cumulative urinary tyrosine conjugates (umol)
HGAPL   : plasma/TBW homogentisic acid (umol)
HGAT    : HGA in avascular matrix interstitium (umol)
BQA     : benzoquinone acetic acid in matrix (umol)
PCART   : ochronotic pigment, articular cartilage (umol) IRREVERSIBLE
PDISC   : ochronotic pigment, intervertebral disc (umol) IRREVERSIBLE
PVALV   : ochronotic pigment, aortic valve (umol) IRREVERSIBLE
PTEND   : ochronotic pigment, tendon/ligament (umol) IRREVERSIBLE
PSCL    : ochronotic pigment, sclera (umol) IRREVERSIBLE
PEAR    : ochronotic pigment, ear cartilage (umol) IRREVERSIBLE
PSKIN   : ochronotic pigment, skin/sweat (umol) IRREVERSIBLE
POTH    : ochronotic pigment, prostate/kidney/dura (umol) IRREVERSIBLE
CART    : intact articular cartilage matrix (fraction of birth value)
FRAGC   : pigmented cartilage shards in joint space (arb)
SYN     : synovitis intensity (arb)
SUBCH   : subchondral sclerosis (arb)
OSTEO   : osteophytosis (arb)
DISCH   : intervertebral disc height (fraction of birth value)
DCALC   : disc calcification burden (arb)
ANKY    : spinal ankylosis index (arb)
VALVCA  : aortic valve calcific mass (arb)
PMAXS   : peak transaortic gradient (mmHg)
LVMI    : LV mass index excess (arb)
STONE   : renal stone burden (arb)
PSTONE  : prostate calculus burden (arb)
RFUN    : relative renal plasma flow (1 = young adult)
TENDI   : tendon tensile integrity (fraction)
BMD     : bone mineral density (fraction)
CORTYR  : corneal tyrosine crystal load (arb)
HEAR    : hearing threshold shift (dB)
SKINP   : visible skin/scleral pigmentation (arb)
NOCI    : peripheral nociceptive drive (arb)
CS      : central sensitisation (0-1)
ANLG    : analgesic concentration (arb)
AKJ     : AKUSSI joint domain (0-1)
AKS     : AKUSSI spine domain (0-1)
AKC     : AKUSSI clinical domain (0-1)
HZJR    : cumulative joint replacement hazard
HZAVR   : cumulative aortic valve replacement hazard
HZRUP   : cumulative tendon rupture hazard
HZKER   : cumulative keratopathy hazard
CFPIG   : counterfactual UNTREATED cartilage pigment (umol)
CFAKC   : counterfactual UNTREATED cAKUSSI trajectory
HDPIG   : headroom integrator: pigment prevented by treatment (umol)
AVOIDI  : avoided serum-HGA exposure integral (umol/L*day)
AUCHGA  : serum HGA exposure (umol/L*day)
AUCTYR  : tyrosine exposure above 700 umol/L (umol/L*day)
CUMNT   : cumulative nitisinone exposure (umol/L*day)

$MAIN
// Deliberately empty of any CMT_0 assignment: init() must win, so that every
// scenario starts from the physiological birth state rather than silently
// resetting it (this was defect #1 during integration).
double LOADETA = exp(ETA(1));
F_NTGUT = FBIO;                 // adherence / bioavailability acts here

$ODE
// ===========================================================================
// 0. DERIVED CONCENTRATIONS
// ===========================================================================
// Positive-part clamps. Without them, stopping nitisinone after 20 years sends
// HPP from ~60 down to ~2 umol/L fast enough that the integrator can overshoot
// CHPP below -KMHPD, at which point the Michaelis denominator (KMAPP + CHPP)
// passes through zero and the whole trajectory becomes NaN. Scenario S23
// (10 mg from 25, stopped at 45) returned NaN for every output while the other
// 23 scenarios looked perfectly healthy -- a failure that only appears on
// withdrawal, i.e. only in the one scenario nobody checks first (defect #12).
double CNT   = NTCEN/V1;   if(CNT   < 0.0) CNT   = 0.0;
double CPHE  = PHE/VAA;    if(CPHE  < 0.0) CPHE  = 0.0;
double CTYR  = TYR/VAA;    if(CTYR  < 0.0) CTYR  = 0.0;
double CHPP  = HPP/VAA;    if(CHPP  < 0.0) CHPP  = 0.0;
double CHPLA = HPLA/VAA;   if(CHPLA < 0.0) CHPLA = 0.0;
double CHGA  = HGAPL/VTBW; if(CHGA  < 0.0) CHGA  = 0.0;
double CHGAT = HGAT/VT;    if(CHGAT < 0.0) CHGAT = 0.0;

// ===========================================================================
// 1. DIETARY INPUT  — the only lever on the SIZE of the flux
// ===========================================================================
double PROTEFF = PROT - AASUPP;
if(PROTEFF < 5.0) PROTEFF = 5.0;
double IN_PHE = PROTEFF*FPHE/165.19*1.0e6;
double IN_TYR = PROTEFF*FTYR/181.19*1.0e6;

// ===========================================================================
// 2. PHE -> TYR
// ===========================================================================
double J_PAH  = VPAH*CPHE/(KMPAH+CPHE)/(1.0 + CTYR/KITYR);
double J_PHEA = VPHEALT*CPHE/(KMPHEALT+CPHE);

// ===========================================================================
// 3. TRANSAMINATION — fast and REVERSIBLE, hence near-equilibrium.
//    HPP therefore tracks TYR proportionally. This is why serum HPPA rises
//    ~14-fold on nitisinone, the same factor as tyrosine, rather than by the
//    enormous factor a one-way transaminase would demand.
// ===========================================================================
double J_TAT = KTAT*(REQTAT*CTYR - CHPP)*VAA;      // + = TYR -> HPP
double J_HL  = KHPLA*(RHL*CHPP - CHPLA)*VAA;       // + = HPP -> HPLA

// ===========================================================================
// 4. HPD AND ITS COMPETITIVE INHIBITION BY NITISINONE
//    Competitive: the apparent Km is raised. At deep inhibition the residual
//    flux is ~ Vmax*CHPP/(KMHPD*CNT/KI), i.e. INVERSELY PROPORTIONAL TO DOSE.
//    That is the origin of the 94% (2 mg) vs 99.5% (10 mg) suppression seen
//    -- a 6-fold difference in residual for a 5-fold dose ratio.
// ===========================================================================
double KMAPP = KMHPD*(1.0 + pow(CNT/KI_NT, HNT));
double J_HPD = VHPD*CHPP/(KMAPP + CHPP);
// The IDEAL comparator blocks HGA production without touching tyrosine flow:
// its blocked flux is routed straight to the alternative exits, so the tyrosine
// penalty is isolated by comparison with nitisinone at matched HGA.
double J_HPD_ID = IDEALDRG*IDEALEFF*J_HPD;
double J_HGA_PROD = J_HPD - J_HPD_ID;

// residual HGD disposal in missense genotypes
double J_HGD = VHGDN*RESACT*CHGA/(KMHGD+CHGA);

// ===========================================================================
// 5. TYROSINE ESCAPE ROUTES — dormant at 54 umol/L, load-bearing at 800.
//    They, and dietary intake, set the tyrosine ceiling. NOTHING here
//    depends on nitisinone dose: that is Balance 1, written as code.
// ===========================================================================
double CT2 = pow(CTYR, HRENTYR);
double J_TYRU  = VRENTYR*CT2/(pow(KMRENTYR,HRENTYR) + CT2);
double J_CONJ  = VCONJ*CTYR/(KMCONJ+CTYR);
double J_HPPU  = CLHPPU*CHPP;
double J_HPLAU = CLHPLAU*CHPLA;

// ===========================================================================
// 6. HGA DISPOSITION — a clearance with no reserve, and a floor no dose clears
// ===========================================================================
// Organic-anion competition: HPPA and HPLA rise 14-fold on treatment and share
// the OAT1/OAT3 secretory route, so CL_HGA FALLS exactly when HGA production is
// suppressed. Serum HGA therefore falls much less than urinary HGA.
double OATCOMP = 1.0 + (CHPP + CHPLA)/KI_OAT;
double CLHGA   = CLHGA0*RFUN/OATCOMP + HGADEG;
double J_HGA_REN = J_HGA_PROD*FREN;                 // intrarenal, bypasses plasma
double J_HGA_SYS = J_HGA_PROD*(1.0-FREN) + SRC_INS; // enters plasma (+ floor)
double J_HGA_OUT = CLHGA*CHGA;                      // leaves plasma to urine
double J_HGA_U   = J_HGA_OUT + J_HGA_REN;           // total urinary HGA

// ===========================================================================
// 7. OXIDATION AND POLYMERISATION — the one-way branch
//    Ascorbate reduces BQA back to HGA, but polymerisation is a PARALLEL fast
//    step, not a step in series, so the futile cycle removes only
//    KRED*CASC/(KRED*CASC+KPOL) of the pigment flux. That is why vitamin C
//    has never worked, and it is structure, not a fitted failure.
// ===========================================================================
double J_DIF  = KDIF*(CHGA - CHGAT)*VT;
double J_OX   = KOX*CHGAT*VT;
double J_POLY = KPOL*BQA;
double J_REDX = KRED*CASC*BQA;

// ===========================================================================
// 8. PIGMENT DENSITIES AND EMBRITTLEMENT
//    Pigment density grows linearly from birth; embrittlement is a steep Hill
//    function of it. Two lines of code, and the third-decade onset plus the
//    subsequent acceleration are outputs rather than assumptions.
// ===========================================================================
double PDC = PCART/COLLC;
double PDD = PDISC/COLLD;
double PDV = PVALV/COLLV;
double PDT = PTEND/COLLT;
double BRC = pow(PDC,HPD_H)/(pow(PDC,HPD_H)+pow(PD50 ,HPD_H));
double BRD = pow(PDD,HPD_H)/(pow(PDD,HPD_H)+pow(PD50D,HPD_H));
double BRV = pow(PDV,HPD_H)/(pow(PDV,HPD_H)+pow(PD50V,HPD_H));
double BRT = pow(PDT,HPD_H)/(pow(PDT,HPD_H)+pow(PD50T,HPD_H));

double LOAD = LOADETA*(BMI/25.0)*OCCUP*(1.0 - EPHYSIO*PHYSIO)*(1.0 + 0.15*SEXM);
double SYNE = SYN*(1.0 - ENSAID*NSAIDON);

// ===========================================================================
// 9. THE AMPLIFYING LOOP: brittle matrix -> shards -> synovitis -> MMP ->
//    more matrix loss -> more brittle-matrix exposure. This loop, not any
//    age term, is why AKUSSI progression accelerates after the fourth decade.
// ===========================================================================
double J_FRAG = KWEAR*LOAD*BRC*CART;
double J_MMP  = KMMP*SYNE*CART;
// Reversible, concentration-driven channel and the rate multiplier it creates.
double OXD = CHGA/(CHGA + KMOX50);
double OXR = OXD/OXREF;
double MOX = FIRR + (1.0-FIRR)*OXR;      // 1 untreated, ~0.39 on 10 mg/day

// ===========================================================================
// 10. TYROSINE OCULAR TOXICITY
// ===========================================================================
double CTY4 = pow(CTYR,HTC);
double J_CRY = KCRY*CTY4/(pow(TC50,HTC)+CTY4);

// ===========================================================================
// 11. STONES
// ===========================================================================
double UHGA_C = J_HGA_U/UVOL;
double SUPSAT = UHGA_C*(1.0 - EALK*ALKALI) - USAT;
if(SUPSAT < 0.0) SUPSAT = 0.0;

// ===========================================================================
// 12. PAIN
// ===========================================================================
double NOCIT = WNOCI_C*(1.0-CART) + WNOCI_S*SYNE + WNOCI_D*(1.0-DISCH)
             + WNOCI_T*(1.0-TENDI);

// ===========================================================================
// 13. COUNTERFACTUAL SHADOW PATIENT (untreated), for measured headroom
// ===========================================================================
// The shadow patient is UNTREATED, so it must not be evaluated at the treated
// HPP concentration of the treated one -- doing so read a 25-fold elevated HPP into
// counterfactual and inflated it about eight-fold, making every headroom
// figure meaningless (defect #10). Balance 1 gives the untreated HGA flux
// directly and without reference to any treated state: essentially the whole
// dietary Phe+Tyr load leaves as HGA. Diet and genotype therefore still move
// the counterfactual, as they should; nitisinone cannot.
double FLUX_CF = IN_PHE + IN_TYR;
double JHGA_CF = FHGACF*FLUX_CF;
double OATCF   = 1.0 + (REQTAT*54.0*(1.0+RHL))/KI_OAT;
double CLCF    = CLHGA0*RFUN/OATCF;
double CHGA_CF = (JHGA_CF*(1.0-FREN) + SRC_INS)
                 /(CLCF + VHGDN*RESACT/(KMHGD+30.0));
double J_PIG_CF  = CFON*KOX*(CHGA_CF*KDIF/(KDIF+KOX))*VT*FCART*KPOL/(KPOL+KRED*CASC);
double J_PIG_ACT = FCART*J_POLY;

// ===========================================================================
// ODEs
// ===========================================================================
// -- nitisinone PK
dxdt_NTGUT = -KA*NTGUT;
dxdt_NTCEN =  KA*NTGUT - (CLNT*CYP3A4/V1)*NTCEN - Q*(NTCEN/V1 - NTPER/V2);
dxdt_NTPER =  Q*(NTCEN/V1 - NTPER/V2);

// -- amino acids and the pathway
dxdt_PHE   = IN_PHE - J_PAH - J_PHEA;
dxdt_TYR   = IN_TYR + J_PAH + KPROTREL - KPROTSYN - J_TAT - J_TYRU - J_CONJ;
dxdt_HPP   = J_TAT - J_HL - J_HPD - J_HPPU;
dxdt_HPLA  = J_HL - J_HPLAU;
dxdt_PROTT = KPROTSYN - KPROTREL;

// -- cumulative urinary output (the trial endpoints live here)
dxdt_CUMHGAU = J_HGA_U;
dxdt_CUMTYRU = J_TYRU;
dxdt_CUMHPPU = J_HPPU;
dxdt_CUMHPLU = J_HPLAU;
dxdt_CUMCONJ = J_CONJ;

// -- HGA disposition, oxidation, polymerisation
dxdt_HGAPL = J_HGA_SYS - J_HGA_OUT - J_DIF - J_HGD;
dxdt_HGAT  = J_DIF - J_OX;
dxdt_BQA   = J_OX - J_POLY - J_REDX;

// -- pigment: eight one-way depots. There is deliberately no removal term.
dxdt_PCART = FCART*J_POLY;
dxdt_PDISC = FDISC*J_POLY;
dxdt_PVALV = FVALV*J_POLY;
dxdt_PTEND = FTEND*J_POLY;
dxdt_PSCL  = FSCL *J_POLY - KTSCL *PSCL;
dxdt_PEAR  = FEAR *J_POLY - KTEAR *PEAR;
dxdt_PSKIN = FSKIN*J_POLY - KTSKIN*PSKIN;
dxdt_POTH  = FOTH *J_POLY;

// -- cartilage and the synovial amplifier
dxdt_CART  = -(J_FRAG + J_MMP)*MOX;
dxdt_FRAGC =  KFRAG*J_FRAG - KCLRF*FRAGC;
dxdt_SYN   =  KSYN*FRAGC + KSYNOX*OXD - KSOFF*SYN;
dxdt_SUBCH =  KSUBCH*(1.0-CART)*BRC - 0.0;
dxdt_OSTEO =  KOSTEO*(1.0-CART);

// -- spine
dxdt_DISCH = -KDISCH*BRD*LOAD*MOX*DISCH;      // multiplicative: cannot go negative
dxdt_DCALC =  KDCALC*BRD*MOX;
dxdt_ANKY  =  KANKY*BRD*MOX*(1.0 + DCALC);

// -- aortic valve
// Deliberately NOT multiplied by MOX. Calcific valve disease is a
// self-perpetuating osteogenic programme once seeded, so lowering HGA should
// not slow it -- which is exactly what the SONIA 2 aortic-stenosis analysis
// found (progression difference 0.0093 mmHg/year, p=0.53). Routing the valve
// through the reversible channel instead predicted -0.51 mmHg/year, a 55-fold
// over-statement of a null result (defect #8).
dxdt_VALVCA = KVC*BRV;
dxdt_PMAXS  = KPM*KVC*BRV*(1.0 + 0.6*VALVCA);
dxdt_LVMI   = KLVMI*(PMAXS-PMAX0)/40.0 - 0.002*LVMI;

// -- stones and renal function (the feedback that raises the pigmenting conc.)
// Supersaturation exists from birth, so a first-order stone model puts stones
// in a toddler. Nucleation needs a nidus, and the histological nidus is
// papillary pigment: gating on renal pigment makes stone burden grow with the
// PRODUCT of supersaturation and accrued pigment, i.e. late (defect #7).
dxdt_STONE  = KSTONE*SUPSAT/1000.0*pow(POTH/COLLT,2.0) - KPASS*STONE;
dxdt_PSTONE = KPSTONE*SUPSAT/1000.0*pow(POTH/COLLT,2.0);
dxdt_RFUN   = -KRPF*CKDX*RFUN*(SOLVERTIME > AGERPF ? 1.0 : 0.0)
              - 2.0e-5*STONE*RFUN;

// -- tendon and bone
dxdt_TENDI = -KTEND*BRT*LOAD*MOX*TENDI;       // multiplicative: bounded below
dxdt_BMD   = -KBMD*SYNE;

// -- ocular, aural, cutaneous
dxdt_CORTYR = J_CRY - KCCLR*CORTYR;
dxdt_HEAR   = KHEAR*(PEAR/COLLT);
dxdt_SKINP  = KSKINP*(PSKIN/COLLT) - KTSKIN*SKINP;   // visible pigment can fade

// -- pain
dxdt_NOCI = 0.2*(NOCIT - NOCI);
dxdt_CS   = KCS*NOCI*(1.0-CS) - KCSOFF*CS;
dxdt_ANLG = KANLG - KANLGCL*ANLG;

// -- AKUSSI domains (0-1 each, weighted into cAKUSSI in $TABLE)
// Bounded 0-1 domain scores: a first-order filter toward an algebraic target,
// so that the composite cannot integrate without limit (defect #4).
double TJ = 1.10*(1.0-CART) + 0.25*OSTEO + 0.20*SUBCH;
double TS = 0.90*(1.0-DISCH) + 0.20*DCALC + 0.25*ANKY;
double TC = 0.35*SKINP + 0.60*(PMAXS-PMAX0)/30.0 + 0.40*STONE
            + 0.35*(1.0-TENDI) + 0.25*HEAR/20.0;
if(TJ > 1.0) TJ = 1.0;  if(TS > 1.0) TS = 1.0;  if(TC > 1.0) TC = 1.0;
dxdt_AKJ = KAKLAG*(TJ - AKJ);
dxdt_AKS = KAKLAG*(TS - AKS);
dxdt_AKC = KAKLAG*(TC - AKC);

// -- cumulative hazards
dxdt_HZJR  = KHZJR *pow(fmax(0.0, 0.85-CART)/0.85, 3.0);
dxdt_HZAVR = KHZAVR*pow(fmax(0.0, PMAXS-36.0)/36.0, 2.0);
dxdt_HZRUP = KHZRUP*pow(fmax(0.0, 0.75-TENDI)/0.75, 2.0);
dxdt_HZKER = 0.02*pow(fmax(0.0, CORTYR-KERTH)/KERTH, 1.5);

// -- counterfactual and headroom (measured, not assumed)
dxdt_CFPIG  = J_PIG_CF;
dxdt_CFAKC  = KAKLAG*(TJ + TS - CFAKC);
dxdt_HDPIG  = fmax(0.0, J_PIG_CF - J_PIG_ACT);
dxdt_AVOIDI = fmax(0.0, CHGA_CF - CHGA);
dxdt_AUCHGA = CHGA;
dxdt_AUCTYR = fmax(0.0, CTYR-700.0);
dxdt_CUMNT  = CNT;

$TABLE
double CNTo   = NTCEN/V1;   if(CNTo   < 0.0) CNTo   = 0.0;
double CPHEo  = PHE/VAA;    if(CPHEo  < 0.0) CPHEo  = 0.0;
double CTYRo  = TYR/VAA;    if(CTYRo  < 0.0) CTYRo  = 0.0;
double CHPPo  = HPP/VAA;    if(CHPPo  < 0.0) CHPPo  = 0.0;
double CHPLAo = HPLA/VAA;   if(CHPLAo < 0.0) CHPLAo = 0.0;
double CHGAo  = HGAPL/VTBW; if(CHGAo  < 0.0) CHGAo  = 0.0;
double CHGATo = HGAT/VT;    if(CHGATo < 0.0) CHGATo = 0.0;
double AGEYo  = TIME/365.25;

// instantaneous 24-h urinary HGA (the SONIA primary endpoint)
double KMAPPo = KMHPD*(1.0 + pow(CNTo/KI_NT, HNT));
double JHPDo  = VHPD*CHPPo/(KMAPPo + CHPPo);
double JHGAPo = JHPDo*(1.0 - IDEALDRG*IDEALEFF);
double OATCo  = 1.0 + (CHPPo + CHPLAo)/KI_OAT;
double CLHGAo = CLHGA0*RFUN/OATCo + HGADEG;
double UHGA24 = CLHGAo*CHGAo + JHGAPo*FREN;
double UTYR24 = VRENTYR*pow(CTYRo,HRENTYR)/(pow(KMRENTYR,HRENTYR)+pow(CTYRo,HRENTYR));
double UHPP24 = CLHPPU*CHPPo;
double UHPL24 = CLHPLAU*CHPLAo;
double UCON24 = VCONJ*CTYRo/(KMCONJ+CTYRo);
double FLUXIN = (PROT-AASUPP)*FPHE/165.19*1.0e6 + (PROT-AASUPP)*FTYR/181.19*1.0e6;
double EXITSUM = UHGA24 + UTYR24 + UHPP24 + UHPL24 + UCON24;

// composite scores and clinical read-outs
double CAKUSSI = WAK_J*AKJ + WAK_S*AKS + WAK_C*AKC + WAK_P*(PAINVAS/10.0);
double PIGTOT  = PCART+PDISC+PVALV+PTEND+PSCL+PEAR+PSKIN+POTH;
double PDCo    = PCART/COLLC;
double BRCo    = pow(PDCo,HPD_H)/(pow(PDCo,HPD_H)+pow(PD50,HPD_H));
double ANALGE  = EANLG*ANLG/(EC50ANLG+ANLG);
double PAINVAS = 10.0*(1.0-exp(-(0.35*NOCI + 0.90*CS)))*(1.0-ANALGE);
double PJR     = 1.0 - exp(-HZJR);      // cumulative prob. joint replacement
double PAVR    = 1.0 - exp(-HZAVR);
double PRUP    = 1.0 - exp(-HZRUP);
double PKER    = 1.0 - exp(-HZKER);
double HEADRM   = (CFPIG > 1e-9) ? HDPIG/CFPIG : 0.0;   // fraction prevented
double PIGSPARE = (CFPIG > 1e-9) ? 1.0 - PCART/CFPIG : 0.0;
double AVA      = 4.0/(1.0 + PMAXS/12.0);
double KERATO   = (CORTYR > KERTH) ? 1.0 : 0.0;

$CAPTURE @annotated
AGEYo  : age (years)
CNTo   : plasma nitisinone (umol/L)
CPHEo  : plasma phenylalanine (umol/L)
CTYRo  : plasma tyrosine (umol/L)
CHPPo  : serum HPPA (umol/L)
CHPLAo : serum HPLA (umol/L)
CHGAo  : serum HGA (umol/L) - THE CAUSAL QUANTITY
CHGATo : cartilage interstitial HGA (umol/L)
UHGA24 : 24-h urinary HGA (umol/day) - THE TRIAL ENDPOINT
UTYR24 : 24-h urinary tyrosine (umol/day)
UHPP24 : 24-h urinary HPPA (umol/day)
UHPL24 : 24-h urinary HPLA (umol/day)
UCON24 : 24-h urinary tyrosine conjugates (umol/day)
FLUXIN : dietary Phe+Tyr input (umol/day)
EXITSUM: sum of all measured exits (umol/day) - mass-balance check
CAKUSSI: clinical AKUSSI score
PIGTOT : total ochronotic pigment (umol HGA-equivalents)
PDCo   : cartilage pigment density
BRCo   : cartilage embrittlement (0-1)
PAINVAS: pain VAS (0-10)
PJR    : cumulative probability of joint replacement
PAVR   : cumulative probability of aortic valve replacement
PRUP   : cumulative probability of tendon rupture
PKER   : cumulative probability of keratopathy
HEADRM : fraction of counterfactual pigment PREVENTED
PIGSPARE: fraction of counterfactual pigment burden spared
AVA    : aortic valve area (cm2, from gradient)
KERATO : clinical keratopathy present (0/1)

$OMEGA @annotated @block
ETA_LOAD : 0.09 : between-subject variability on mechanical load
'

## ---------------------------------------------------------------------------
## Build
## ---------------------------------------------------------------------------
aku_model <- function() {
  mcode("aku", AKU_CODE, quiet = TRUE)
}

## Physiological initial conditions AT BIRTH.
## Deliberately set through init(), never through *_0 in $MAIN.
aku_init <- function(mod) {
  init(mod,
       NTGUT = 0, NTCEN = 0, NTPER = 0,
       PHE  = 55 * 25,      # 55 umol/L x 25 L
       TYR  = 54 * 25,
       HPP  = 3.0 * 25,
       HPLA = 12.9 * 25,
       PROTT = 1.73e6,
       CUMHGAU = 0, CUMTYRU = 0, CUMHPPU = 0, CUMHPLU = 0, CUMCONJ = 0,
       HGAPL = 28 * 42, HGAT = 26 * 0.5, BQA = 0.0002,
       PCART = 0, PDISC = 0, PVALV = 0, PTEND = 0,
       PSCL = 0, PEAR = 0, PSKIN = 0, POTH = 0,
       CART = 1, FRAGC = 0, SYN = 0, SUBCH = 0, OSTEO = 0,
       DISCH = 1, DCALC = 0, ANKY = 0,
       VALVCA = 0, PMAXS = 4, LVMI = 0,
       STONE = 0, PSTONE = 0, RFUN = 1,
       TENDI = 1, BMD = 1,
       CORTYR = 0, HEAR = 0, SKINP = 0,
       NOCI = 0, CS = 0, ANLG = 0,
       AKJ = 0, AKS = 0, AKC = 0,
       HZJR = 0, HZAVR = 0, HZRUP = 0, HZKER = 0,
       CFPIG = 0, CFAKC = 0, HDPIG = 0, AVOIDI = 0,
       AUCHGA = 0, AUCTYR = 0, CUMNT = 0)
}

NT_MG_TO_UMOL <- 1e3 / 329.25    # nitisinone MW 329.25 g/mol; 1 mg = 3.037 umol

## ---------------------------------------------------------------------------
## Simulation helper.  Dosing always given in mg/day; start age in years.
## Uses modifyList() rather than c() so that duplicated parameter names cannot
## silently produce two identical arms (a defect this model was built with and
## which is documented at the end of the file).
## ---------------------------------------------------------------------------
sim_aku <- function(mod,
                    dose_mg    = 0,
                    start_age  = 25,
                    stop_age   = Inf,
                    end_age    = 70,
                    pars       = list(),
                    delta      = 30.4375,
                    atol       = 1e-10,
                    rtol       = 1e-8,
                    iiv        = FALSE,
                    label      = "scenario") {

  ## DETERMINISTIC BY DEFAULT.  The model declares an $OMEGA on mechanical load,
  ## and mrgsim draws it afresh on every call. Every arm of every comparison was
  ## therefore a DIFFERENT patient: repeating the age-at-initiation scan gave
  ## 0.826 / 0.707 / 0.735 / 0.608 one time and 0.842 / 0.761 / 0.599 / 0.535 the
  ## next, and the headline result (cartilage preserved at 70 years is monotone in
  ## age at initiation) appeared to fail. It was first misdiagnosed as a solver
  ## tolerance problem; tightening atol from 1e-8 to 1e-12 did not fix it and
  ## could not, because the noise was never numerical. zero_re() is the fix.
  ## Pass iiv = TRUE only for deliberate population simulation (defect #11).
  mod <- aku_init(mod)
  if (!iiv) mod <- zero_re(mod)
  if (length(pars)) mod <- param(mod, do.call(c, list(pars)))

  tend  <- end_age * 365.25
  tstart <- start_age * 365.25
  tstop  <- if (is.infinite(stop_age)) tend else stop_age * 365.25

  if (dose_mg > 0 && tstart < tend) {
    n <- max(0, floor((min(tstop, tend) - tstart)) )
    e <- ev(time = tstart, amt = dose_mg * NT_MG_TO_UMOL, ii = 1, addl = n)
  } else {
    e <- ev(time = tend + 10, amt = 0)   # no-treatment arm
  }

  out <- mod %>%
    ev(e) %>%
    ## With the random effect zeroed (see above) this model is NOT stiff in any
    ## way that matters: the age-at-initiation scan is identical to seven decimal
    ## places at atol 1e-8 and at 1e-10, so 1e-10 / 1e-8 is a comfortable margin
    ## rather than a requirement. The apparent tolerance sensitivity that was
    ## chased first was the unseeded $OMEGA, not the integrator.
    mrgsim(end = tend, delta = delta, atol = atol, rtol = rtol,
           maxsteps = 2000000) %>%
    as_tibble()
  out$scenario <- label
  out
}

## ---------------------------------------------------------------------------
## Steady-state probe: run an ADULT at fixed dose for 4 years and read the
## biochemistry, to compare against SONIA 1 / SONIA 2 / NAC.
## ---------------------------------------------------------------------------
probe_biochem <- function(mod, dose_mg, years = 4, start_age = 46, pars = list()) {
  d <- if (years < 0.5) 1 else 30.4375
  o <- sim_aku(mod, dose_mg = dose_mg, start_age = start_age,
               end_age = start_age + years, pars = pars, delta = d,
               label = paste0(dose_mg, "mg"))
  tail(o, 1)
}

## ===========================================================================
## CALIBRATION
## ---------------------------------------------------------------------------
## Six biochemical anchors, all from randomised or prospectively collected
## cohorts.  Fitted parameters: KI_NT, SRC_INS, KI_OAT, CLHPPU, CLHPLAU, VCONJ.
##   A1  untreated 24-h urinary HGA          31,530 umol/day  (SONIA 1 control)
##   A2  untreated serum HGA                     28 umol/L    (NAC/SONIA 2)
##   A3  u-HGA24 on 2 mg/day                  1,200 umol/day  (NAC)
##   A4  u-HGA24 on 10 mg/day at 12 months       181 umol/day  (SONIA 2)
##   A5  serum HGA on 2 mg / 10 mg           3.86 / 2.23 umol/L
##   A6  plasma tyrosine on 2 mg / 10 mg      782 / 875 umol/L
## ===========================================================================
AKU_ANCHORS <- tibble::tribble(
  ~anchor,            ~target,   ~unit,
  "uHGA24_0mg",       31530,     "umol/day",   # SONIA 1, no-treatment arm
  "uHGA24_1mg",        3260,     "umol/day",   # SONIA 1, 1 mg/day
  "uHGA24_2mg",        1440,     "umol/day",   # SONIA 1, 2 mg/day
  "uHGA24_4mg",         570,     "umol/day",   # SONIA 1, 4 mg/day
  "uHGA24_8mg",         150,     "umol/day",   # SONIA 1, 8 mg/day
  "sHGA_0mg",            28,     "umol/L",     # untreated serum HGA
  "sTYR_2mg",           782,     "umol/L",     # NAC, 2 mg, 48 months
  "sTYR_10mg",          875,     "umol/L"      # SONIA 2, 10 mg, 48 months
)

aku_anchor_pred <- function(mod, pars = list()) {
  # SONIA 1 was 4 weeks of once-daily dosing in adults: probe at the same
  # duration so that the tyrosine/HPPA build-up is at the same stage.
  b <- lapply(c(0,1,2,4,8), function(d)
    probe_biochem(mod, d, years = 4/52, start_age = 46, pars = pars))
  b0  <- probe_biochem(mod, 0,  4, 46, pars)
  b2  <- probe_biochem(mod, 2,  4, 46, pars)
  b10 <- probe_biochem(mod, 10, 4, 46, pars)
  c(uHGA24_0mg = b[[1]]$UHGA24,
    uHGA24_1mg = b[[2]]$UHGA24,
    uHGA24_2mg = b[[3]]$UHGA24,
    uHGA24_4mg = b[[4]]$UHGA24,
    uHGA24_8mg = b[[5]]$UHGA24,
    sHGA_0mg   = b0$CHGAo,
    sTYR_2mg   = b2$CTYRo,
    sTYR_10mg  = b10$CTYRo)
}

aku_objective <- function(mod, pars = list()) {
  p <- aku_anchor_pred(mod, pars)
  t <- setNames(AKU_ANCHORS$target, AKU_ANCHORS$anchor)
  mean(((log(p[names(t)]) - log(t))^2))
}

## Coordinate-wise refinement.  Nelder-Mead was tried first and wandered into
## regions where the tyrosine ceiling collapsed (CTYR < 100), because the
## objective is nearly flat in CLHPPU/CLHPLAU along the direction that keeps
## their sum constant; coordinate descent on log-parameters is stable here.
calibrate_aku <- function(mod, rounds = 3, verbose = TRUE) {
  fitp <- c(KI_NT = 0.016, HNT = 1.45, SRC_INS = 60, KI_OAT = 10,
            CLHGA0 = 1968, CLHPPU = 160, CLHPLAU = 55,
            VCONJ = 14000, VRENTYR = 8000)
  best <- aku_objective(mod, as.list(fitp))
  for (r in seq_len(rounds)) {
    for (nm in names(fitp)) {
      for (f in c(0.7, 0.85, 1.18, 1.4)) {
        trial <- fitp; trial[nm] <- fitp[nm]*f
        o <- try(aku_objective(mod, as.list(trial)), silent = TRUE)
        if (!inherits(o, "try-error") && is.finite(o) && o < best) {
          best <- o; fitp <- trial
        }
      }
    }
    if (verbose) cat(sprintf("  round %d  objective = %.5f\n", r, best))
  }
  list(pars = as.list(fitp), objective = best)
}

## ===========================================================================
## SCENARIOS — 24 arms, deliberately built as MATCHED PAIRS so that each
## comparison changes one thing.  The point of a matched pair is that the
## conclusion cannot come from a difference in something else.
## ===========================================================================
aku_scenarios <- function() {
  list(
    ## --- natural history and trial replications -------------------------
    list(id = "S01", label = "Untreated natural history (birth to 70 y)",
         dose = 0, start = 25, end = 70, pars = list()),
    list(id = "S02", label = "SONIA 1: 1 mg/day, 4 weeks (adult)",
         dose = 1, start = 46, end = 46.08, pars = list(), delta = 1),
    list(id = "S03", label = "SONIA 1: 2 mg/day, 4 weeks",
         dose = 2, start = 46, end = 46.08, pars = list(), delta = 1),
    list(id = "S04", label = "SONIA 1: 4 mg/day, 4 weeks",
         dose = 4, start = 46, end = 46.08, pars = list(), delta = 1),
    list(id = "S05", label = "SONIA 1: 8 mg/day, 4 weeks",
         dose = 8, start = 46, end = 46.08, pars = list(), delta = 1),
    list(id = "S06", label = "SONIA 2: 10 mg/day from age 49, 4 years",
         dose = 10, start = 49, end = 53, pars = list()),
    list(id = "S07", label = "SONIA 2 control: no treatment, age 49-53",
         dose = 0, start = 49, end = 53, pars = list()),
    list(id = "S08", label = "NIH trial: 2 mg/day from age 46, 3 years",
         dose = 2, start = 46, end = 49, pars = list()),

    ## --- age at initiation: identical drug, different remaining integral --
    list(id = "S09", label = "Nitisinone 10 mg from age 5",
         dose = 10, start = 5, end = 70, pars = list()),
    list(id = "S10", label = "Nitisinone 10 mg from age 15",
         dose = 10, start = 15, end = 70, pars = list()),
    list(id = "S11", label = "Nitisinone 10 mg from age 25",
         dose = 10, start = 25, end = 70, pars = list()),
    list(id = "S12", label = "Nitisinone 10 mg from age 40",
         dose = 10, start = 40, end = 70, pars = list()),
    list(id = "S13", label = "Nitisinone 10 mg from age 55",
         dose = 10, start = 55, end = 70, pars = list()),

    ## --- dose vs diet: the model's central asymmetry, as a 2x2 -----------
    list(id = "S14", label = "10 mg + free diet (84 g protein)",
         dose = 10, start = 30, end = 70, pars = list(PROT = 84)),
    list(id = "S15", label = "10 mg + protein restriction (56 g = 0.8 g/kg)",
         dose = 10, start = 30, end = 70, pars = list(PROT = 56)),
    list(id = "S16", label = "2 mg + free diet (84 g protein)",
         dose = 2, start = 30, end = 70, pars = list(PROT = 84)),
    list(id = "S17", label = "2 mg + protein restriction (56 g)",
         dose = 2, start = 30, end = 70, pars = list(PROT = 56)),

    ## --- interventions that are not nitisinone ---------------------------
    list(id = "S18", label = "Diet alone: 56 g protein, no drug",
         dose = 0, start = 30, end = 70, pars = list(PROT = 56)),
    list(id = "S19", label = "Ascorbate 1 g/day (CASC 120), no drug",
         dose = 0, start = 30, end = 70, pars = list(CASC = 120)),

    ## --- disease modifiers, each matched to S11 --------------------------
    list(id = "S20", label = "Missense genotype, 3% residual HGD, untreated",
         dose = 0, start = 25, end = 70, pars = list(RESACT = 0.03)),
    list(id = "S21", label = "CKD (renal decline x3), untreated",
         dose = 0, start = 25, end = 70, pars = list(CKDX = 3)),
    list(id = "S22", label = "Heavy manual work + BMI 32, untreated",
         dose = 0, start = 25, end = 70, pars = list(OCCUP = 1.6, BMI = 32)),

    ## --- treatment interruption, and the ideal comparator ----------------
    list(id = "S23", label = "10 mg from age 25, stopped at age 45",
         dose = 10, start = 25, stop = 45, end = 70, pars = list()),
    list(id = "S24", label = "IDEAL drug from age 25: HGA block, no Tyr rise",
         dose = 0, start = 25, end = 70, pars = list(IDEALDRG = 1))
  )
}

run_scenarios <- function(mod, fitted = list(), which = NULL) {
  scn <- aku_scenarios()
  if (!is.null(which)) scn <- scn[vapply(scn, function(s) s$id %in% which, TRUE)]
  out <- lapply(scn, function(s) {
    pars <- modifyList(fitted, s$pars)
    o <- sim_aku(mod,
                 dose_mg   = s$dose,
                 start_age = s$start,
                 stop_age  = if (is.null(s$stop)) Inf else s$stop,
                 end_age   = s$end,
                 pars      = pars,
                 delta     = if (is.null(s$delta)) 30.4375 else s$delta,
                 label     = paste(s$id, s$label))
    o$id <- s$id
    o
  })
  bind_rows(out)
}

## ===========================================================================
## VALIDATION — held-out facts the calibration never saw
## ===========================================================================
validate_aku <- function(mod, fitted = list()) {
  f <- function(p) modifyList(fitted, p)
  ## delta() at a given age -- the baseline must be read at the treatment start
  ## age, not at t = 0 (birth); reading it at birth inflated every treatment
  ## contrast by the whole preceding natural history (defect #9).
  at <- function(o, age, col) o[[col]][which(o$AGEYo >= age)[1]]

  b0   <- probe_biochem(mod, 0,  4, 46, f(list()))
  b2   <- probe_biochem(mod, 2,  4, 46, f(list()))
  b10  <- probe_biochem(mod, 10, 1, 46, f(list()))
  b10d <- probe_biochem(mod, 10, 4, 46, f(list(PROT = 49)))

  nat  <- sim_aku(mod, 0, 25, Inf, 78, f(list()), label = "nat")
  onset <- nat$AGEYo[which(nat$PAINVAS > 2)[1]]
  disc3 <- nat$AGEYo[which(nat$DCALC   > 0.25)[1]]
  jr50  <- nat$AGEYo[which(nat$PJR     > 0.5)[1]]
  valv  <- nat$AGEYo[which(nat$PMAXS   > 12)[1]]
  st10  <- nat$AGEYo[which(nat$STONE   > 0.50)[1]]
  bri   <- nat$AGEYo[which(nat$BRCo    > 0.5)[1]]

  s2n <- sim_aku(mod, 10, 49, Inf, 53, f(list()), label = "s2n")
  s2c <- sim_aku(mod,  0, 49, Inf, 53, f(list()), label = "s2c")
  s2b <- sim_aku(mod,  2, 49, Inf, 53, f(list()), label = "s2b")
  dn  <- at(s2n,53,"CAKUSSI") - at(s2n,49,"CAKUSSI")
  dc  <- at(s2c,53,"CAKUSSI") - at(s2c,49,"CAKUSSI")
  db  <- at(s2b,53,"CAKUSSI") - at(s2b,49,"CAKUSSI")
  dpm <- ((at(s2n,53,"PMAXS")-at(s2n,49,"PMAXS")) -
          (at(s2c,53,"PMAXS")-at(s2c,49,"PMAXS")))/4

  asc <- sim_aku(mod, 0, 30, Inf, 70, f(list(CASC = 120)), label = "asc")
  ascv <- 1 - at(asc,70,"PIGTOT")/at(nat,70,"PIGTOT")

  ## reversal of pigmentation is possible only where collagen turns over
  rev2 <- sim_aku(mod, 10, 49, Inf, 51, f(list()), label = "rev")
  skin_rev <- 1 - at(rev2,51,"SKINP")/at(rev2,49,"SKINP")
  cart_rev <- 1 - at(rev2,51,"PDCo") /at(rev2,49,"PDCo")

  tibble::tribble(
    ~test,                                    ~model,          ~reported,
    "serum HGA on 10 mg (umol/L)",            b10$CHGAo,        2.23,
    "serum HGA on 2 mg (umol/L)",             b2$CHGAo,         3.86,
    "serum HPPA fold-rise on 10 mg",          b10$CHPPo/b0$CHPPo, 14.65,
    "residual u-HGA24 ratio 2 mg / 10 mg",    b2$UHGA24/b10$UHGA24, 6.63,
    "tyrosine ratio 10 mg / 2 mg",            b10$CTYRo/b2$CTYRo, 1.119,
    "tyrosine on 10 mg + 0.7 g/kg protein",   b10d$CTYRo,     620,
    "mass-balance closure (exits / intake)",  b0$EXITSUM/b0$FLUXIN, 1.00,
    "age cartilage embrittlement > 0.5 (y)",  bri,             31,
    "age disc calcification appears (y)",     disc3,           30,
    "age pain VAS > 2 (y)",                   onset,           30,
    "age 50% joint replacement (y)",          jr50,            55,
    "age aortic valve involvement (y)",       valv,            54,
    "age renal stones detectable (y)",        st10,            64,
    "SONIA 2 control cAKUSSI slope (pts/mo)", dc/48,            0.239,
    "SONIA 2 10 mg cAKUSSI slope (pts/mo)",   dn/48,            0.060,
    "NAC 2 mg cAKUSSI slope (pts/mo)",        db/48,            0.190,
    "SONIA 2 cAKUSSI difference at 48 mo",    dn-dc,           -8.6,
    "SONIA 2 aortic Pmax slope diff (mmHg/y)",dpm,             -0.0093,
    "pigment reduction by ascorbate 1 g/day", ascv,             0.00,
    "skin pigment fade, 2 y on 10 mg",        skin_rev,         0.20,
    "cartilage pigment fade, 2 y on 10 mg",   cart_rev,         0.00
  ) %>% mutate(ratio = ifelse(abs(reported) < 1e-9, NA, model/reported))
}

## ===========================================================================
## HEADLINE ANALYSES
## ===========================================================================

## (1) Dose sets HGA, diet sets tyrosine — the 2x2 that proves Balance 1.
analysis_dose_vs_diet <- function(mod, fitted = list()) {
  grid <- expand.grid(dose = c(0, 1, 2, 4, 8, 10, 20),
                      prot = c(42, 56, 70, 84, 105))
  res <- Map(function(d, p) {
    b <- probe_biochem(mod, d, 3, 46, modifyList(fitted, list(PROT = p)))
    tibble(dose_mg = d, protein_g = p,
           sHGA = b$CHGAo, uHGA24 = b$UHGA24, sTYR = b$CTYRo,
           pigment_rate = b$CHGATo*0.04*0.5)
  }, grid$dose, grid$prot)
  bind_rows(res)
}

## (2) Measured headroom: identical drug, different age at initiation.
analysis_headroom <- function(mod, fitted = list()) {
  ages <- c(2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60)
  bind_rows(lapply(ages, function(a) {
    o <- sim_aku(mod, 10, a, Inf, 70, fitted, label = paste0("init", a))
    e <- tail(o, 1)
    tibble(init_age = a,
           pigment_umol = e$PCART, counterfactual = e$CFPIG,
           spared_frac = e$PIGSPARE, headroom_used = e$HEADRM,
           cAKUSSI_70 = e$CAKUSSI, cart_70 = e$CART,
           p_joint_repl = e$PJR, pain_70 = e$PAINVAS)
  }))
}

## (3) Urine vs serum: the endpoint and the causal quantity, side by side.
analysis_endpoint_gap <- function(mod, fitted = list()) {
  bind_rows(lapply(c(0, 1, 2, 4, 8, 10, 20, 40), function(d) {
    b <- probe_biochem(mod, d, 3, 46, fitted)
    b0 <- probe_biochem(mod, 0, 3, 46, fitted)
    tibble(dose_mg = d,
           uHGA_pct_drop = 100*(1 - b$UHGA24/b0$UHGA24),
           sHGA_pct_drop = 100*(1 - b$CHGAo /b0$CHGAo),
           sTYR = b$CTYRo,
           implied_pigment_rate_pct = 100*b$CHGAo/b0$CHGAo)
  }))
}

## ===========================================================================
## DRIVER
## ===========================================================================
main_aku <- function(do_calibrate = TRUE) {
  mod <- aku_model()
  cat("\n== Alkaptonuria QSP model ==\n")
  cat(sprintf("compartments: %d   parameters: %d\n",
              length(mrgsolve::init(mod)), length(mrgsolve::param(mod))))

  fitted <- list()
  if (do_calibrate) {
    cat("\n-- calibration --\n")
    cal <- calibrate_aku(mod, rounds = 2)
    fitted <- cal$pars
    cat("fitted:\n"); print(unlist(fitted))
  }

  cat("\n-- anchors --\n")
  p <- aku_anchor_pred(mod, fitted)
  t <- setNames(AKU_ANCHORS$target, AKU_ANCHORS$anchor)
  print(tibble(anchor = names(t), model = round(p[names(t)],2),
               target = t, ratio = round(p[names(t)]/t, 3)))

  cat("\n-- held-out validation --\n")
  print(as.data.frame(validate_aku(mod, fitted)), digits = 3)

  cat("\n-- dose vs diet (Balance 1) --\n")
  print(as.data.frame(analysis_dose_vs_diet(mod, fitted)), digits = 3)

  cat("\n-- urine/serum endpoint gap (Balance 3) --\n")
  print(as.data.frame(analysis_endpoint_gap(mod, fitted)), digits = 3)

  cat("\n-- measured headroom by age at initiation (Balance 2) --\n")
  print(as.data.frame(analysis_headroom(mod, fitted)), digits = 3)

  cat("\n-- 24 scenarios --\n")
  sc <- run_scenarios(mod, fitted)
  print(sc %>% group_by(id, scenario) %>% slice_tail(n = 1) %>%
        select(id, AGEYo, CTYRo, CHGAo, UHGA24, CAKUSSI, PIGTOT, PJR, PKER) %>%
        as.data.frame(), digits = 3)

  invisible(list(mod = mod, fitted = fitted, scenarios = sc))
}

if (identical(environment(), globalenv()) &&
    !is.null(getOption("aku.run")) && isTRUE(getOption("aku.run"))) {
  main_aku()
}

## ============================================================================
##  DEFECTS FOUND AND FIXED DURING INTEGRATION (mrgsolve 2.0.1)
##  See alkaptonuria/README.md for the full account and for the held-out
##  failures that are reported rather than repaired.
## ============================================================================
