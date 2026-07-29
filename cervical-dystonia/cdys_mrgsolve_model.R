## =====================================================================
##  cdys_mrgsolve_model.R
##  Cervical dystonia (spasmodic torticollis) -- QSP model for mrgsolve
##
##  70 ODEs.  Botulinum neurotoxin chemodenervation PK/PD across eight neck
##  muscles, diffusive spread to the pharynx and to autonomic targets, humoral
##  immunity against two serotypes, the central dystonic-drive ratchet, three
##  oral adjuncts, and the TWSTRS subscales.
##
##  THIS FILE IS THE SECOND OF TWO PORTS OF THE SAME EQUATIONS.
##  The first is `cdys_reference_check.py` (numpy/scipy + LSODA).  This build
##  environment has no R toolchain, so the numbers quoted in README.md come from
##  the Python port.  The two files are meant to be diffed against each other:
##  same states, same parameters, same right-hand side.  If they disagree, one of
##  them is wrong, and that is the point of having two.
##
##  Quick start
##  -----------
##      library(mrgsolve); library(dplyr)
##      mod <- mread("cdys_mrgsolve_model.R")
##      out <- sim_scenario(mod, "S03_std_240U_q12wk")
##      plot(out, TWSTRS_TOTAL ~ time)
##
##  THE CENTRAL RESULT, in one line of algebra
##  ------------------------------------------
##      L(t) >= L_min = (1 - phi) * DCEN ,    phi = RHO * sum(w_m over injected)
##  No dose, potency, serotype or dilution term appears on the right.  phi is a
##  property of the injection PLAN and of the patient's anatomy.  RHO -- the
##  share of a targeted muscle the injectate actually reaches -- is identified by
##  the PLATEAU of the clinical dose-response curve, which is therefore a
##  measurement of geometry and not of pharmacology.  See A2 in the Python port.
## =====================================================================

$PROB
# Cervical dystonia QSP model
- 70 ODEs; 8 muscle compartments x 5 states
- Endpoints: TWSTRS severity (0-35), disability (0-30), pain (0-20), total (0-85)
- Interventions: 5 botulinum products, dose/dilution/target-list/placement,
  trihexyphenidyl, baclofen, clonazepam, GPi-DBS, selective peripheral denervation

$PLUGIN base

$GLOBAL
// E(S): transmission efficacy as a function of the intact SNARE fraction.
// Two nonlinearities, both real and both load-bearing:
//   (1) SNARE cooperativity -- several intact complexes are needed per fusion
//       event, so release capacity RCAP(S) is a Hill function of S;
//   (2) the neuromuscular SAFETY FACTOR -- a healthy terminal releases about
//       SF times more quanta than are needed to fire the fibre, so release
//       capacity must fall BELOW 1/SF before any weakness appears.
// Together they make E(S) a THRESHOLD function.  That is why half the substrate
// can be destroyed with no measurable weakness, why dose -> duration is
// logarithmic, and why the clinical dose-response plateaus.
#define RCAP(S) (pow(fmax((S),0.0), NSNARE) / (pow(fmax((S),0.0), NSNARE) + pow(S50, NSNARE)))
#define EEFF(S) ((SF * RCAP(S) / (1.0 + SF * RCAP(S))) / (SF / (1.0 + SF)))
#define CLAMP01(x) (fmin(fmax((x), 0.0), 1.0))

## ---------------------------------------------------------------------
## PARAMETERS
## ---------------------------------------------------------------------
$PARAM @annotated
// ---- toxin disposition in muscle ----
KBIND    :  4.20  : Productive binding to nerve terminals (1/d)
KCLEAR   :  1.00  : Local proteolysis and lymphatic clearance (1/d)
KDIFF0   :  0.90  : Diffusive efflux at the reference dilution (1/d)
GV       :  1.00  : Volume exponent for diffusive efflux (-)
KTRANS   :  0.60  : Endosomal translocation of the light chain (1/d)
KLC      :  0.02465: Loss of catalytically active light chain (1/d) CALIBRATED
BMAXM    :  260.0 : Saturable terminal binding capacity per unit mass (U)

// ---- SNARE cleavage and resynthesis ----
KSYN     :  0.140 : SNAP-25 resynthesis (1/d)
KCL      :  0.13613: Cleavage rate per (U-eq per unit mass) (1/d) CALIBRATED

// ---- neuromuscular safety factor ----
NSNARE   :  3.0   : SNARE cooperativity Hill coefficient (-)
S50      :  0.20  : Intact-SNARE fraction at half release capacity (-)
SF       :  3.0   : Neuromuscular safety factor (-)

// ---- terminal sprouting ----
KSP      :  0.050 : Sprout outgrowth (1/d)
KRG      :  0.050 : Sprout regression once the parent terminal recovers (1/d)
WQ       :  0.55  : Maximum share of the deficit sprouts can carry (-)

// ---- geometry: the only place phi is set ----
RHO      :  0.5463: Share of a TARGETED muscle the injectate reaches CALIBRATED
RHOMULT  :  1.0   : Multiplier on RHO for EMG/ultrasound guidance (-)
DENERV   :  0.0   : Fraction of NON-injected drive removed surgically (-)

// ---- spread compartments ----
THSCALE  :  1.00  : Global scale on the muscle-to-pharynx proximity index (-)
MASSSW   :  0.35  : Relative mass of the swallow compartment (-)
MASSAU   :  0.50  : Relative mass of the autonomic compartment (-)
THETAAU  :  0.06  : Efflux share reaching autonomic targets (-)
AUTOPREF :  1.0   : Autonomic tropism multiplier (serotype B ~4) (-)

// ---- clinical mapping ----
SEVMAX   :  35.0  : TWSTRS severity ceiling (points)
HS       :  1.60  : Hill coefficient, torque load to severity (-)
L50      :  0.8347: Torque load at half maximal severity (-)
KSEV     :  0.25  : Severity equilibration (1/d)
PAINMAX  :  20.0  : TWSTRS pain ceiling (points)
KPON     :  0.030 : Nociceptive sensitisation on-rate (1/d)
KPOFF    :  0.030 : Nociceptive sensitisation off-rate (1/d)
QP       :  1.50  : Torque exponent driving sensitisation (-)
PAINDIR  :  0.0   : Direct antinociceptive BoNT action (0 in base; see A13) (-)
DISMAX   :  30.0  : TWSTRS disability ceiling (points)
GDIS     :  0.80  : Disability scaling (-)
WDISSEV  :  0.55  : Weight of severity vs pain in disability (-)
KDIS     :  0.080 : Disability equilibration (1/d)

// ---- central and spinal ----
KD       :  0.0080: Maladaptive-plasticity time constant (1/d)
GA       :  0.80  : Afferent-drive exponent (-)
AFFFLOOR :  0.35  : Spindle output surviving complete intrafusal block (-)
RI0      :  0.45  : Reciprocal inhibition in CD, 1 = healthy (-)
SI0      :  0.40  : Surround inhibition in CD, 1 = healthy (-)
ARI      :  0.50  : Weight of the reciprocal-inhibition deficit on drive gain (-)
ASI      :  0.50  : Weight of the surround-inhibition deficit on drive gain (-)
KRI      :  0.050 : Reciprocal-inhibition equilibration (1/d)
KSI      :  0.050 : Surround-inhibition equilibration (1/d)
KCB      :  0.020 : Cerebellar-gain equilibration (1/d)
CB0      :  1.00  : Cerebellar gain set point (-)

// ---- adverse-effect tolerance distributions ----
DYSD50   :  1.7158: Swallow deficit at 50% dysphagia CALIBRATED
DYSK     :  0.5426: Logistic width, dysphagia CALIBRATED
NWD50    :  0.6536: Posterior deficit at 50% neck weakness CALIBRATED
NWK      :  0.130 : Logistic width, neck weakness (-)
DYSTHR   :  0.30  : Swallow deficit above which burden accrues (-)

// ---- immunogenicity ----
KAG      :  0.50  : Antigen depot clearance (1/d)
KB       :  0.000082: Naive memory-B priming per ng CALIBRATED
BETABOOST:  9.0   : Memory recall amplification (-)
KAAG     :  40.0  : Antigen-processing saturation (ng)
KBD      :  0.0019: Memory-B decay (1/d)
KNP      :  0.55  : Nab production per unit memory B (1/d)
KND      :  0.0347: IgG elimination (1/d)
NAB50    :  1.00  : Nab giving 50% loss of injected potency (-)
HN       :  1.60  : Hill coefficient of the potency gate (-)

// ---- oral adjunct PK/PD ----
THPKA    :  36.0  : Trihexyphenidyl absorption (1/d)
THPKE    :  4.50  : Trihexyphenidyl elimination (1/d)
THPV     :  560.0 : Trihexyphenidyl volume (L)
THPEMAX  :  0.30  : Maximal fractional reduction of central drive gain (-)
THPEC50  :  8.0   : Trihexyphenidyl EC50 (ng/mL)
BACKA    :  30.0  : Baclofen absorption (1/d)
BACKE    :  4.75  : Baclofen elimination (1/d)
BACV     :  60.0  : Baclofen volume (L)
BACEMAX  :  0.12  : Maximal fractional reduction of central drive gain (-)
BACEC50  :  200.0 : Baclofen EC50 (ng/mL)
CLZKA    :  24.0  : Clonazepam absorption (1/d)
CLZKE    :  0.555 : Clonazepam elimination (1/d)
CLZV     :  210.0 : Clonazepam volume (L)
CLZERI   :  0.36  : Fractional restoration of the RecInh deficit (-)
CLZESI   :  0.25  : Fractional restoration of the SurrInh deficit (-)
CLZEC50  :  12.0  : Clonazepam EC50 (ng/mL)
THPDOSE  :  0.0   : Trihexyphenidyl dose (mg/d)
BACDOSE  :  0.0   : Baclofen dose (mg/d)
CLZDOSE  :  0.0   : Clonazepam dose (mg/d)

// ---- GPi deep brain stimulation ----
DBSON    : -1.0   : Day of DBS activation, negative = off (d)
DBSEMAX  :  0.50  : Maximal fractional reduction of central drive gain (-)
DBSK     :  0.0154: DBS ramp rate (1/d)

// ---- bound mode: force transmission to zero in the injected muscles ----
BLOCK    :  0.0   : 1 = perfect permanent blockade of the injected muscles (-)

// ---- reference baseline for the accounting compartments ----
TW0REF   :  42.94 : Untreated baseline TWSTRS total (points)

$PARAM @annotated @covariates
// torque shares w_m; they must sum to 1
W1 : 0.20 : Sternocleidomastoid, contralateral (-)
W2 : 0.24 : Splenius capitis, ipsilateral (-)
W3 : 0.13 : Trapezius, ipsilateral (-)
W4 : 0.08 : Levator scapulae (-)
W5 : 0.12 : Semispinalis capitis, DEEP (-)
W6 : 0.06 : Scalene complex (-)
W7 : 0.10 : Obliquus capitis inferior, DEEP (-)
W8 : 0.07 : Longus colli, effectively unreachable (-)
// relative muscle mass; toxin acts on C/mass, a concentration
M1 : 1.00 : Relative mass, SCM (-)
M2 : 1.60 : Relative mass, splenius (-)
M3 : 2.20 : Relative mass, trapezius (-)
M4 : 0.80 : Relative mass, levator scapulae (-)
M5 : 1.80 : Relative mass, semispinalis (-)
M6 : 0.90 : Relative mass, scalenes (-)
M7 : 0.50 : Relative mass, obliquus capitis inferior (-)
M8 : 0.70 : Relative mass, longus colli (-)
// proximity of each muscle to the swallow compartment
T1 : 0.30 : Pharyngeal proximity index, SCM (-)
T2 : 0.08 : Pharyngeal proximity index, splenius (-)
T3 : 0.04 : Pharyngeal proximity index, trapezius (-)
T4 : 0.05 : Pharyngeal proximity index, levator scapulae (-)
T5 : 0.12 : Pharyngeal proximity index, semispinalis (-)
T6 : 0.28 : Pharyngeal proximity index, scalenes (-)
T7 : 0.18 : Pharyngeal proximity index, obliquus capitis inferior (-)
T8 : 0.60 : Pharyngeal proximity index, longus colli (-)
// posterior (head extensor) flag, for the neck-weakness read-out
P1 : 0 : Posterior flag, SCM (-)
P2 : 1 : Posterior flag, splenius (-)
P3 : 1 : Posterior flag, trapezius (-)
P4 : 1 : Posterior flag, levator scapulae (-)
P5 : 1 : Posterior flag, semispinalis (-)
P6 : 0 : Posterior flag, scalenes (-)
P7 : 1 : Posterior flag, obliquus capitis inferior (-)
P8 : 0 : Posterior flag, longus colli (-)
// injectate volume per muscle, mL.  Reference dilution = 50 U/mL, i.e.
// 0.02 mL per Unit.  V=0 means "use the reference dilution for that dose".
V1 : 0 : Injectate volume, SCM (mL)
V2 : 0 : Injectate volume, splenius (mL)
V3 : 0 : Injectate volume, trapezius (mL)
V4 : 0 : Injectate volume, levator scapulae (mL)
V5 : 0 : Injectate volume, semispinalis (mL)
V6 : 0 : Injectate volume, scalenes (mL)
V7 : 0 : Injectate volume, obliquus capitis inferior (mL)
V8 : 0 : Injectate volume, longus colli (mL)
// nominal Units per muscle, used only to derive the reference volume and to
// mark which muscles are on the target list (INJ_m)
D1 : 50 : Nominal dose, SCM (U)
D2 : 90 : Nominal dose, splenius (U)
D3 : 60 : Nominal dose, trapezius (U)
D4 : 40 : Nominal dose, levator scapulae (U)
D5 :  0 : Nominal dose, semispinalis (U)
D6 :  0 : Nominal dose, scalenes (U)
D7 :  0 : Nominal dose, obliquus capitis inferior (U)
D8 :  0 : Nominal dose, longus colli (U)

## ---------------------------------------------------------------------
## STATE VECTOR -- 70 ODEs
## ---------------------------------------------------------------------
$CMT @annotated
// --- muscle 1: sternocleidomastoid (contralateral) ---
A1_ : Free bioactive toxin, SCM (U)
B1_ : Membrane-bound / internalising toxin, SCM (U)
C1_ : Active light chain in cytosol, SCM (U-eq)
S1_ : Intact SNARE fraction, SCM (-)
Q1_ : Sprout-mediated release capacity, SCM (-)
// --- muscle 2: splenius capitis (ipsilateral) ---
A2_ : Free bioactive toxin, splenius (U)
B2_ : Internalising toxin, splenius (U)
C2_ : Active light chain, splenius (U-eq)
S2_ : Intact SNARE fraction, splenius (-)
Q2_ : Sprout capacity, splenius (-)
// --- muscle 3: trapezius ---
A3_ : Free bioactive toxin, trapezius (U)
B3_ : Internalising toxin, trapezius (U)
C3_ : Active light chain, trapezius (U-eq)
S3_ : Intact SNARE fraction, trapezius (-)
Q3_ : Sprout capacity, trapezius (-)
// --- muscle 4: levator scapulae ---
A4_ : Free bioactive toxin, levator scapulae (U)
B4_ : Internalising toxin, levator scapulae (U)
C4_ : Active light chain, levator scapulae (U-eq)
S4_ : Intact SNARE fraction, levator scapulae (-)
Q4_ : Sprout capacity, levator scapulae (-)
// --- muscle 5: semispinalis capitis (deep) ---
A5_ : Free bioactive toxin, semispinalis (U)
B5_ : Internalising toxin, semispinalis (U)
C5_ : Active light chain, semispinalis (U-eq)
S5_ : Intact SNARE fraction, semispinalis (-)
Q5_ : Sprout capacity, semispinalis (-)
// --- muscle 6: scalene complex ---
A6_ : Free bioactive toxin, scalenes (U)
B6_ : Internalising toxin, scalenes (U)
C6_ : Active light chain, scalenes (U-eq)
S6_ : Intact SNARE fraction, scalenes (-)
Q6_ : Sprout capacity, scalenes (-)
// --- muscle 7: obliquus capitis inferior (deep) ---
A7_ : Free bioactive toxin, OCI (U)
B7_ : Internalising toxin, OCI (U)
C7_ : Active light chain, OCI (U-eq)
S7_ : Intact SNARE fraction, OCI (-)
Q7_ : Sprout capacity, OCI (-)
// --- muscle 8: longus colli / deep prevertebral ---
A8_ : Free bioactive toxin, longus colli (U)
B8_ : Internalising toxin, longus colli (U)
C8_ : Active light chain, longus colli (U-eq)
S8_ : Intact SNARE fraction, longus colli (-)
Q8_ : Sprout capacity, longus colli (-)
// --- swallow compartment (pharyngeal constrictors) ---
ASW : Free toxin, swallow compartment (U)
BSW : Internalising toxin, swallow compartment (U)
CSW : Active light chain, swallow compartment (U-eq)
SSW : Intact SNARE fraction, swallow compartment (-)
// --- autonomic compartment (salivary / ganglionic) ---
AAU : Free toxin, autonomic compartment (U)
BAU : Internalising toxin, autonomic compartment (U)
CAU : Active light chain, autonomic compartment (U-eq)
SAU : Intact SNARE fraction, autonomic compartment (-)
// --- humoral immunity, serotype A ---
AGA  : Antigen depot, serotype A (ng)
BMEMA: Memory B cells, serotype A (-)
NABA : Neutralising IgG, serotype A (-)
// --- humoral immunity, serotype B ---
AGB  : Antigen depot, serotype B (ng)
BMEMB: Memory B cells, serotype B (-)
NABB : Neutralising IgG, serotype B (-)
// --- central and spinal ---
DCEN : Central dystonic drive, maladaptive plasticity state (-)
RECI : Reciprocal inhibition (-)
SURI : Surround inhibition (-)
CBLL : Cerebellar gain (-)
// --- clinical ---
SEV  : TWSTRS severity subscale (points)
PAIN : Nociceptive sensitisation fraction (-)
DISAB: TWSTRS disability subscale (points)
// --- oral adjunct PK ---
THPD : Trihexyphenidyl depot (ng)
THPC : Trihexyphenidyl central (ng)
BACD : Baclofen depot (ng)
BACC : Baclofen central (ng)
CLZD : Clonazepam depot (ng)
CLZC : Clonazepam central (ng)
// --- accounting ---
AUCBEN : Cumulative TWSTRS benefit-time (point.d)
AUCDYS : Cumulative swallow-deficit burden (deficit.d)
CUMU   : Cumulative Units administered (U)

## ---------------------------------------------------------------------
## INITIAL CONDITIONS
## ---------------------------------------------------------------------
$MAIN
// every SNARE pool starts intact
S1_0 = 1.0; S2_0 = 1.0; S3_0 = 1.0; S4_0 = 1.0;
S5_0 = 1.0; S6_0 = 1.0; S7_0 = 1.0; S8_0 = 1.0;
SSW_0 = 1.0; SAU_0 = 1.0;

DCEN_0 = 1.0;
RECI_0 = RI0;
SURI_0 = SI0;
CBLL_0 = CB0;

// baseline severity is the steady state of the clinical map at L = 1
double sev0  = SEVMAX * pow(1.0, HS) / (pow(1.0, HS) + pow(L50, HS));
double pain0 = KPON / (KPON + KPOFF);
SEV_0   = sev0;
PAIN_0  = pain0;
DISAB_0 = DISMAX * GDIS * (WDISSEV * sev0 / SEVMAX + (1.0 - WDISSEV) * pain0);

## ---------------------------------------------------------------------
## RIGHT-HAND SIDE
## ---------------------------------------------------------------------
$ODE

// ============ local helper values ==================================
double rho = fmin(RHO * RHOMULT, 1.0);
double vref = 0.02;                         // mL per Unit at 50 U/mL

// per-muscle geometry, dose flags and dilution-dependent efflux
double w[8]  = {W1,W2,W3,W4,W5,W6,W7,W8};
double ms[8] = {M1,M2,M3,M4,M5,M6,M7,M8};
double th[8] = {T1,T2,T3,T4,T5,T6,T7,T8};
double po[8] = {P1,P2,P3,P4,P5,P6,P7,P8};
double dn[8] = {D1,D2,D3,D4,D5,D6,D7,D8};
double vl[8] = {V1,V2,V3,V4,V5,V6,V7,V8};

double A[8] = {A1_,A2_,A3_,A4_,A5_,A6_,A7_,A8_};
double B[8] = {B1_,B2_,B3_,B4_,B5_,B6_,B7_,B8_};
double C[8] = {C1_,C2_,C3_,C4_,C5_,C6_,C7_,C8_};
double S[8] = {S1_,S2_,S3_,S4_,S5_,S6_,S7_,S8_};
double Q[8] = {Q1_,Q2_,Q3_,Q4_,Q5_,Q6_,Q7_,Q8_};

double dA[8], dB[8], dC[8], dS[8], dQ[8];
double fluxsw = 0.0, fluxau = 0.0;
double Lsum = 0.0, affsum = 0.0, postnum = 0.0, postden = 0.0;

for (int i = 0; i < 8; ++i) {
  double Si = fmin(fmax(S[i], 0.0), 1.0);
  double Qi = fmin(fmax(Q[i], 0.0), 1.0);
  double inj = (dn[i] > 0.0) ? 1.0 : 0.0;
  double blk = (BLOCK > 0.5) ? inj : 0.0;      // bound mode

  // dilution: efflux scales with injectate volume relative to the reference
  double vr = (dn[i] > 0.0 && vl[i] > 0.0) ? (vl[i] / (dn[i] * vref)) : 1.0;
  double kdiff = KDIFF0 * pow(fmax(vr, 1e-6), GV);

  // A -> B -> C cascade.  Note that the drug is essentially gone within a day;
  // the EFFECT is set by C, S and Q, which are all much slower.
  double occ = fmin(fmax(1.0 - B[i] / (BMAXM * ms[i]), 0.0), 1.0);
  dA[i] = -(KBIND * occ + KCLEAR + kdiff) * A[i];
  dB[i] =  KBIND * occ * A[i] - KTRANS * B[i];
  dC[i] =  KTRANS * B[i] - KLC * C[i];
  dS[i] =  KSYN * (1.0 - Si) - KCL * (C[i] / ms[i]) * Si;

  double Em = EEFF(Si) * (1.0 - blk);
  dQ[i] = KSP * (1.0 - Em) * (1.0 - Qi) - KRG * Em * Qi;

  // Sprouts grow from the SAME axon and load vesicles from the SAME cytosolic
  // SNARE pool the light chain has destroyed, so their release is gated by S
  // too.  Sprouting therefore ACCELERATES recovery; it does not cap blockade.
  // A5 shows this makes sprouting a minor contributor to wearing-off, which
  // contradicts the usual account.
  double Tm  = 1.0 - (1.0 - Em) * (1.0 - WQ * Qi * RCAP(Si) * (1.0 - blk));
  double Tef = 1.0 - rho * (1.0 - Tm);        // only rho of the muscle reached

  // selective peripheral denervation permanently removes NON-injected drive
  double weff = (DENERV > 0.0 && inj < 0.5) ? w[i] * (1.0 - DENERV) : w[i];
  Lsum += weff * Tef;

  // spindle afferent drive; complete blockade takes out the intrafusal
  // (gamma) terminals as well
  double spindle = (blk > 0.5) ? AFFFLOOR
                               : AFFFLOOR + (1.0 - AFFFLOOR) * Si;
  affsum += w[i] * (1.0 - rho * (1.0 - spindle));

  fluxsw += kdiff * th[i] * THSCALE * A[i];
  fluxau += kdiff * THETAAU * AUTOPREF * A[i];

  postnum += w[i] * po[i] * (1.0 - Tef);
  postden += w[i] * po[i];
}

dxdt_A1_ = dA[0]; dxdt_B1_ = dB[0]; dxdt_C1_ = dC[0]; dxdt_S1_ = dS[0]; dxdt_Q1_ = dQ[0];
dxdt_A2_ = dA[1]; dxdt_B2_ = dB[1]; dxdt_C2_ = dC[1]; dxdt_S2_ = dS[1]; dxdt_Q2_ = dQ[1];
dxdt_A3_ = dA[2]; dxdt_B3_ = dB[2]; dxdt_C3_ = dC[2]; dxdt_S3_ = dS[2]; dxdt_Q3_ = dQ[2];
dxdt_A4_ = dA[3]; dxdt_B4_ = dB[3]; dxdt_C4_ = dC[3]; dxdt_S4_ = dS[3]; dxdt_Q4_ = dQ[3];
dxdt_A5_ = dA[4]; dxdt_B5_ = dB[4]; dxdt_C5_ = dC[4]; dxdt_S5_ = dS[4]; dxdt_Q5_ = dQ[4];
dxdt_A6_ = dA[5]; dxdt_B6_ = dB[5]; dxdt_C6_ = dC[5]; dxdt_S6_ = dS[5]; dxdt_Q6_ = dQ[5];
dxdt_A7_ = dA[6]; dxdt_B7_ = dB[6]; dxdt_C7_ = dC[6]; dxdt_S7_ = dS[6]; dxdt_Q7_ = dQ[6];
dxdt_A8_ = dA[7]; dxdt_B8_ = dB[7]; dxdt_C8_ = dC[7]; dxdt_S8_ = dS[7]; dxdt_Q8_ = dQ[7];

// ============ spread compartments ==================================
// The pharyngeal constrictors are SMALL, so a small absolute amount of toxin
// produces a large concentration there.  That asymmetry, not any special
// affinity, is what makes dysphagia the dose-limiting toxicity.
double occsw = fmin(fmax(1.0 - BSW / (BMAXM * MASSSW), 0.0), 1.0);
dxdt_ASW = fluxsw - (KBIND * occsw + KCLEAR) * ASW;
dxdt_BSW = KBIND * occsw * ASW - KTRANS * BSW;
dxdt_CSW = KTRANS * BSW - KLC * CSW;
dxdt_SSW = KSYN * (1.0 - fmin(fmax(SSW,0.0),1.0)) - KCL * (CSW / MASSSW) * fmin(fmax(SSW,0.0),1.0);

double occau = fmin(fmax(1.0 - BAU / (BMAXM * MASSAU), 0.0), 1.0);
dxdt_AAU = fluxau - (KBIND * occau + KCLEAR) * AAU;
dxdt_BAU = KBIND * occau * AAU - KTRANS * BAU;
dxdt_CAU = KTRANS * BAU - KLC * CAU;
dxdt_SAU = KSYN * (1.0 - fmin(fmax(SAU,0.0),1.0)) - KCL * (CAU / MASSAU) * fmin(fmax(SAU,0.0),1.0);

// ============ humoral immunity, two independent serotype pools =====
// The memory compartment makes short intervals worse than a linear cumulative
// dose would, because recall is amplified while memory persists.  A7 shows
// that, at the OBSERVED antibody rate, this still does not produce an interior
// optimum -- so the q12wk convention is not explained by antibody risk.
double primeA = KB * AGA / (1.0 + AGA / KAAG) * (1.0 + BETABOOST * BMEMA);
dxdt_AGA   = -KAG * AGA;
dxdt_BMEMA = primeA - KBD * BMEMA;
dxdt_NABA  = KNP * BMEMA - KND * NABA;

double primeB = KB * AGB / (1.0 + AGB / KAAG) * (1.0 + BETABOOST * BMEMB);
dxdt_AGB   = -KAG * AGB;
dxdt_BMEMB = primeB - KBD * BMEMB;
dxdt_NABB  = KNP * BMEMB - KND * NABB;

// ============ oral adjunct PK =======================================
dxdt_THPD = THPDOSE * 1e6 - THPKA * THPD;
dxdt_THPC = THPKA * THPD - THPKE * THPC;
dxdt_BACD = BACDOSE * 1e6 - BACKA * BACD;
dxdt_BACC = BACKA * BACD - BACKE * BACC;
dxdt_CLZD = CLZDOSE * 1e6 - CLZKA * CLZD;
dxdt_CLZC = CLZKA * CLZD - CLZKE * CLZC;

double cthp = THPC / THPV / 1000.0;          // ng/mL
double cbac = BACC / BACV / 1000.0;
double cclz = CLZC / CLZV / 1000.0;

// ============ spinal and cortical gating ============================
double fclz = cclz / (cclz + CLZEC50);
dxdt_RECI = KRI * ((RI0 + (1.0 - RI0) * CLZERI * fclz) - RECI);
dxdt_SURI = KSI * ((SI0 + (1.0 - SI0) * CLZESI * fclz) - SURI);
dxdt_CBLL = KCB * (CB0 - CBLL);

// ============ central drive: the ratchet ============================
// This is the route from a PERIPHERAL injection to a SLOW CENTRAL benefit:
// BoNT blocks the intrafusal terminals too, lowering the abnormal spindle
// afferent drive that feeds the plasticity state.  A9 shows the consequence --
// troughs that fall across cycles while peak paralysis is unchanged.
double Graw = (1.0 + ARI * (1.0 - RECI)) * (1.0 + ASI * (1.0 - SURI));
double Gref = (1.0 + ARI * (1.0 - RI0))  * (1.0 + ASI * (1.0 - SI0));
double dbs  = (DBSON >= 0.0 && SOLVERTIME >= DBSON)
              ? DBSEMAX * (1.0 - exp(-DBSK * (SOLVERTIME - DBSON))) : 0.0;
double Gcen = (Graw / Gref) * (1.0 - dbs)
              * (1.0 - THPEMAX * cthp / (cthp + THPEC50))
              * (1.0 - BACEMAX * cbac / (cbac + BACEC50))
              * CBLL;
dxdt_DCEN = KD * (Gcen * pow(fmax(affsum, 1e-9), GA) - DCEN);

// ============ dystonic torque load and clinical scores =============
// THE BOUND:  L >= (1 - phi) * DCEN  with phi = rho * sum(w over injected).
// No drug parameter appears in it.
double L = DCEN * Lsum;

dxdt_SEV = KSEV * (SEVMAX * pow(fmax(L,0.0), HS)
                   / (pow(fmax(L,0.0), HS) + pow(L50, HS)) - SEV);

double Pn   = fmin(fmax(PAIN, 0.0), 1.0);
double anti = fmax(1.0 - PAINDIR * postnum, 0.0);
dxdt_PAIN = KPON * pow(fmax(L,0.0), QP) * anti * (1.0 - Pn) - KPOFF * Pn;

dxdt_DISAB = KDIS * (DISMAX * GDIS * (WDISSEV * SEV / SEVMAX
                                      + (1.0 - WDISSEV) * Pn) - DISAB);

// ============ accounting ============================================
double twnow  = SEV + DISAB + PAINMAX * Pn;
double defsw  = 1.0 - EEFF(fmin(fmax(SSW,0.0),1.0));
dxdt_AUCBEN = fmax(TW0REF - twnow, 0.0);
dxdt_AUCDYS = fmax(defsw - DYSTHR, 0.0);
dxdt_CUMU   = 0.0;

## ---------------------------------------------------------------------
## DERIVED OUTPUTS
## ---------------------------------------------------------------------
$TABLE
double rho_o = fmin(RHO * RHOMULT, 1.0);
double phi_o = rho_o * ( (D1>0?W1:0) + (D2>0?W2:0) + (D3>0?W3:0) + (D4>0?W4:0)
                       + (D5>0?W5:0) + (D6>0?W6:0) + (D7>0?W7:0) + (D8>0?W8:0) );

double TWSTRS_TOTAL = SEV + DISAB + PAINMAX * fmin(fmax(PAIN,0.0),1.0);
double TWSTRS_PAIN  = PAINMAX * fmin(fmax(PAIN,0.0),1.0);
double GAIN         = TW0REF - TWSTRS_TOTAL;

double DEF_SW = 1.0 - EEFF(fmin(fmax(SSW,0.0),1.0));
double DEF_AU = 1.0 - EEFF(fmin(fmax(SAU,0.0),1.0));

// posterior (head-extensor) deficit -> neck weakness / head drop
double pn_ = 0.0, pd_ = 0.0;
{
  double w_[8]  = {W1,W2,W3,W4,W5,W6,W7,W8};
  double po_[8] = {P1,P2,P3,P4,P5,P6,P7,P8};
  double dn_[8] = {D1,D2,D3,D4,D5,D6,D7,D8};
  double S_[8]  = {S1_,S2_,S3_,S4_,S5_,S6_,S7_,S8_};
  double Q_[8]  = {Q1_,Q2_,Q3_,Q4_,Q5_,Q6_,Q7_,Q8_};
  for (int i = 0; i < 8; ++i) {
    double Si = fmin(fmax(S_[i],0.0),1.0), Qi = fmin(fmax(Q_[i],0.0),1.0);
    double blk = (BLOCK > 0.5 && dn_[i] > 0.0) ? 1.0 : 0.0;
    double Em = EEFF(Si) * (1.0 - blk);
    double Tm = 1.0 - (1.0 - Em) * (1.0 - WQ * Qi * RCAP(Si) * (1.0 - blk));
    double Tef = 1.0 - rho_o * (1.0 - Tm);
    pn_ += w_[i] * po_[i] * (1.0 - Tef);
    pd_ += w_[i] * po_[i];
  }
}
double POST_DEF = (pd_ > 0.0) ? pn_ / pd_ : 0.0;

double P_DYSPHAGIA = 1.0 / (1.0 + exp(-(DEF_SW   - DYSD50) / DYSK));
double P_NECKWEAK  = 1.0 / (1.0 + exp(-(POST_DEF - NWD50)  / NWK));
double P_DRYMOUTH  = 1.0 / (1.0 + exp(-(DEF_AU   - DYSD50) / DYSK));

// potency remaining against each serotype
double GATE_A = 1.0 / (1.0 + pow(NABA / NAB50, HN));
double GATE_B = 1.0 / (1.0 + pow(NABB / NAB50, HN));

double S_SCM = S1_, S_SPLEN = S2_, Q_SCM = Q1_, C_SCM = C1_;
double E_SCM = EEFF(fmin(fmax(S1_,0.0),1.0));

$CAPTURE @annotated
TWSTRS_TOTAL : TWSTRS total score, 0-85 (points)
SEV          : TWSTRS severity subscale, 0-35 (points)
DISAB        : TWSTRS disability subscale, 0-30 (points)
TWSTRS_PAIN  : TWSTRS pain subscale, 0-20 (points)
GAIN         : TWSTRS improvement vs untreated baseline (points)
DCEN         : Central dystonic drive (-)
phi_o        : phi = rho * sum(w) over injected muscles (-)
S_SCM        : Intact SNARE fraction, sternocleidomastoid (-)
E_SCM        : Parent-terminal transmission efficacy, SCM (-)
Q_SCM        : Sprout capacity, SCM (-)
C_SCM        : Active light chain, SCM (U-eq)
S_SPLEN      : Intact SNARE fraction, splenius (-)
DEF_SW       : Swallow-compartment transmission deficit (-)
POST_DEF     : Posterior (head-extensor) transmission deficit (-)
P_DYSPHAGIA  : Modelled probability of dysphagia (-)
P_NECKWEAK   : Modelled probability of neck weakness (-)
P_DRYMOUTH   : Modelled probability of dry mouth (-)
NABA         : Neutralising IgG, serotype A (-)
NABB         : Neutralising IgG, serotype B (-)
GATE_A       : Fraction of injected serotype-A potency surviving antibody (-)
GATE_B       : Fraction of injected serotype-B potency surviving antibody (-)
AUCBEN       : Cumulative benefit-time (point.d)
AUCDYS       : Cumulative swallow-deficit burden (deficit.d)

## =====================================================================
$ENV
## ---------------------------------------------------------------------
## PRODUCTS
##   unit_scale : label Unit -> internal Allergan-equivalent Unit.  These are
##                clinical rules of thumb, NOT equipotency claims.
##   load       : neurotoxin-complex protein per label Unit (ng/U) = THE ANTIGEN
##   klc_mult   : multiplier on the active light-chain loss rate
## ---------------------------------------------------------------------
PRODUCTS <- list(
  incobotulinumtoxinA = list(serotype="A", unit_scale=1.00, load=0.44/100,
                             klc_mult=1.00, auto_pref=1.0),
  onabotulinumtoxinA  = list(serotype="A", unit_scale=1.00, load=5.00/100,
                             klc_mult=1.00, auto_pref=1.0),
  abobotulinumtoxinA  = list(serotype="A", unit_scale=0.34, load=4.35/500,
                             klc_mult=1.00, auto_pref=1.0),
  daxibotulinumtoxinA = list(serotype="A", unit_scale=1.00, load=0.50/100,
                             klc_mult=0.50, auto_pref=1.0),
  rimabotulinumtoxinB = list(serotype="B", unit_scale=0.03, load=5.00/5000,
                             klc_mult=1.45, auto_pref=4.0)
)

## Injection patterns, Allergan-equivalent Units per muscle.
## 240 U over four muscles is the pivotal-trial dose.
PATTERN_STD <- c(50, 90, 60, 40,  0, 0,  0, 0)   # 240 U, phi ~ 0.36
PATTERN_EXT <- c(45, 70, 40, 25, 40, 0, 20, 0)   # 240 U, phi ~ 0.48
CMT_A <- c("A1_","A2_","A3_","A4_","A5_","A6_","A7_","A8_")

## ---------------------------------------------------------------------
## build_dosing() -- turn an injection plan into an mrgsolve event table
##
## Every injection is FIVE things at once and the helper keeps them together:
##   (1) a bolus of bioactive toxin into each targeted muscle's A compartment,
##       scaled by the product's unit conversion AND by the surviving potency
##       after neutralising antibody;
##   (2) a bolus of ANTIGEN into the serotype-matched depot;
##   (3) the nominal per-muscle Units D1..D8, which mark the target list;
##   (4) the injectate volumes V1..V8, which set diffusive efflux;
##   (5) the product-specific light-chain decay rate.
##
## The antibody gate is applied at each injection time from the simulated NABA
## at that moment, so long horizons need the two-pass helper sim_scenario()
## below rather than a single call.
## ---------------------------------------------------------------------
build_dosing <- function(pattern, product = "incobotulinumtoxinA",
                         interval = 84, n_inj = 1, first = 0,
                         upermL = 50, gates = NULL) {
  pr <- PRODUCTS[[product]]
  stopifnot(!is.null(pr))
  times <- first + (seq_len(n_inj) - 1) * interval
  if (is.null(gates)) gates <- rep(1, n_inj)
  ev <- NULL
  for (k in seq_len(n_inj)) {
    for (i in which(pattern > 0)) {
      e <- data.frame(time = times[k], cmt = CMT_A[i],
                      amt = pattern[i] * pr$unit_scale * gates[k],
                      evid = 1)
      ev <- rbind(ev, e)
    }
    agcmt <- if (pr$serotype == "A") "AGA" else "AGB"
    ev <- rbind(ev, data.frame(time = times[k], cmt = agcmt,
                               amt = sum(pattern) * pr$load, evid = 1))
  }
  attr(ev, "params") <- c(
    stats::setNames(as.list(pattern), paste0("D", 1:8)),
    stats::setNames(as.list(pattern / upermL), paste0("V", 1:8)),
    list(KLC = 0.02465 * pr$klc_mult, AUTOPREF = pr$auto_pref)
  )
  attr(ev, "serotype") <- pr$serotype
  ev
}

## ---------------------------------------------------------------------
## SCENARIOS -- 14 of them
## ---------------------------------------------------------------------
SCENARIOS <- list(

  S01_placebo = list(
    label = "Placebo (no injection)",
    pattern = rep(0, 8), product = "incobotulinumtoxinA", n_inj = 0, end = 168,
    note = "reference trajectory; confirms the untreated baseline is a steady state"),

  S02_std_240U_single = list(
    label = "240 U, standard 4-muscle surface plan, single injection",
    pattern = PATTERN_STD, product = "incobotulinumtoxinA", n_inj = 1, end = 210,
    note = "THE CALIBRATION SCENARIO: week-4 and week-12 TWSTRS anchors"),

  S03_std_240U_q12wk = list(
    label = "240 U q12wk x 10 cycles",
    pattern = PATTERN_STD, product = "incobotulinumtoxinA",
    interval = 84, n_inj = 10, end = 840,
    note = "chronic therapy; shows the central ratchet in the TROUGHS (A9)"),

  S04_dose_120U = list(
    label = "120 U, half dose",
    pattern = PATTERN_STD * 0.5, product = "incobotulinumtoxinA", n_inj = 1, end = 210,
    note = "held-back check: response should be shallow and sub-proportional"),

  S05_dose_480U = list(
    label = "480 U, double dose",
    pattern = PATTERN_STD * 2, product = "incobotulinumtoxinA", n_inj = 1, end = 260,
    note = "the PLATEAU that identifies rho; note dysphagia keeps climbing"),

  S06_extended_target_list = list(
    label = "240 U redistributed over 6 muscles (semispinalis + OCI added)",
    pattern = PATTERN_EXT, product = "incobotulinumtoxinA", n_inj = 1, end = 210,
    note = "raises Sigma-w at constant Units; buys little ALONE (A2 point ii)"),

  S07_guided_placement = list(
    label = "240 U, standard list, PERFECT needle placement (rho -> 1)",
    pattern = PATTERN_STD, product = "incobotulinumtoxinA", n_inj = 1, end = 260,
    params = list(RHOMULT = 1 / 0.5463),
    note = "the single largest non-dose lever in the model (A2 lever 2)"),

  S08_extended_and_guided = list(
    label = "240 U, extended list AND perfect placement",
    pattern = PATTERN_EXT, product = "incobotulinumtoxinA", n_inj = 1, end = 260,
    params = list(RHOMULT = 1 / 0.5463),
    note = "the two geometry levers are COMPLEMENTS, not substitutes"),

  S09_bound_perfect_blockade = list(
    label = "THE BOUND: perfect permanent blockade of the injected muscles",
    pattern = PATTERN_STD, product = "incobotulinumtoxinA",
    interval = 84, n_inj = 8, end = 672,
    params = list(BLOCK = 1),
    note = "not a therapy; the upper limit of EVERY therapy at this phi"),

  S10_dilution_high_volume = list(
    label = "240 U at 12.5 U/mL (4x the reference volume)",
    pattern = PATTERN_STD, product = "incobotulinumtoxinA", n_inj = 1, end = 210,
    upermL = 12.5,
    note = "efficacy nearly flat, dysphagia risk is not: A6"),

  S11_product_comparison_daxi = list(
    label = "daxibotulinumtoxinA 240 U (slower light-chain loss)",
    pattern = PATTERN_STD, product = "daxibotulinumtoxinA", n_inj = 1, end = 400,
    note = "longer duration at the SAME nadir -- the ceiling does not move (A11)"),

  S12_serotype_rescue_B = list(
    label = "rimabotulinumtoxinB at A-equivalent dose",
    pattern = PATTERN_STD / 0.03, product = "rimabotulinumtoxinB", n_inj = 1, end = 300,
    note = "different substrate (VAMP), no cross-neutralisation, autonomic cost"),

  S13_adjunct_clonazepam = list(
    label = "240 U q12wk + clonazepam 2 mg/d",
    pattern = PATTERN_STD, product = "incobotulinumtoxinA",
    interval = 84, n_inj = 3, end = 252,
    params = list(CLZDOSE = 2),
    note = "GATING operator: restores part of the RecInh / SurrInh deficit"),

  S14_gpi_dbs = list(
    label = "240 U q12wk + GPi deep brain stimulation from day 0",
    pattern = PATTERN_STD, product = "incobotulinumtoxinA",
    interval = 84, n_inj = 3, end = 252,
    params = list(DBSON = 0),
    note = "CENTRAL DRIVE operator; the only strong one that is not a needle"),

  S15_denervation = list(
    label = "240 U q12wk + selective peripheral denervation (40% of unreached drive)",
    pattern = PATTERN_STD, product = "incobotulinumtoxinA",
    interval = 84, n_inj = 3, end = 252,
    params = list(DENERV = 0.40),
    note = "the only intervention here that permanently RELOCATES the ceiling")
)

## ---------------------------------------------------------------------
## sim_scenario() -- run one scenario
##
## Two passes.  The first pass simulates without the antibody gate to obtain
## NABA/NABB at each injection time; the second applies the gate.  For a single
## injection the two passes are identical, so the extra work only matters on the
## long chronic horizons where antibody actually accumulates.
## ---------------------------------------------------------------------
sim_scenario <- function(mod, name, delta = 0.5, passes = 2) {
  sc <- SCENARIOS[[name]]
  if (is.null(sc)) stop("unknown scenario: ", name)
  n_inj <- if (is.null(sc$n_inj)) 1 else sc$n_inj
  interval <- if (is.null(sc$interval)) 84 else sc$interval
  upermL <- if (is.null(sc$upermL)) 50 else sc$upermL
  end <- if (is.null(sc$end)) 210 else sc$end

  if (n_inj == 0 || all(sc$pattern <= 0)) {
    m <- mod
    if (!is.null(sc$params)) m <- mrgsolve::param(m, sc$params)
    return(mrgsolve::mrgsim(m, end = end, delta = delta))
  }

  gates <- rep(1, n_inj)
  out <- NULL
  for (pass in seq_len(passes)) {
    ev <- build_dosing(sc$pattern, sc$product, interval, n_inj,
                       upermL = upermL, gates = gates)
    m <- mrgsolve::param(mod, attr(ev, "params"))
    if (!is.null(sc$params)) m <- mrgsolve::param(m, sc$params)
    out <- mrgsolve::mrgsim(m, data = ev, end = end, delta = delta,
                            recover = "cmt")
    d <- as.data.frame(out)
    nabcol <- if (attr(ev, "serotype") == "A") "NABA" else "NABB"
    nab50 <- mrgsolve::param(m)$NAB50
    hn <- mrgsolve::param(m)$HN
    tt <- (seq_len(n_inj) - 1) * interval
    gates <- vapply(tt, function(t0) {
      nb <- stats::approx(d$time, d[[nabcol]], xout = t0, rule = 2)$y
      1 / (1 + (nb / nab50)^hn)
    }, numeric(1))
  }
  out
}

## ---------------------------------------------------------------------
## Read-outs matching the Python port
## ---------------------------------------------------------------------
tw_at <- function(out, day) {
  d <- as.data.frame(out)
  stats::approx(d$time, d$TWSTRS_TOTAL, xout = day, rule = 2)$y
}

nadir_tw <- function(out, lo = 0, hi = Inf) {
  d <- as.data.frame(out); d <- d[d$time >= lo & d$time <= hi, ]
  i <- which.min(d$TWSTRS_TOTAL)
  list(tw = d$TWSTRS_TOTAL[i], day = d$time[i])
}

## Duration of benefit: days until the gain falls back below the minimal
## clinically important change.  A PREDICTION, not a calibration anchor.
MCID_TWSTRS <- 4.5
duration_of_benefit <- function(out, mcid = MCID_TWSTRS) {
  d <- as.data.frame(out)
  g <- d$GAIN
  i <- which.max(g)
  if (g[i] <= mcid) return(0)
  j <- which(seq_along(g) > i & g <= mcid)
  if (!length(j)) return(max(d$time))
  d$time[min(j)] - d$time[1]
}

## Mean gain across a steady-state cycle -- the currency used in A2, because it
## is the only read-out on which a wearing-off therapy and a permanent one can
## be compared honestly.
cycle_mean_gain <- function(out, cycle = 8, interval = 84) {
  d <- as.data.frame(out)
  m <- d$time >= (cycle - 1) * interval & d$time <= cycle * interval
  c(mean = mean(d$GAIN[m]), best = max(d$GAIN[m]), worst = min(d$GAIN[m]),
    p_dysphagia = max(d$P_DYSPHAGIA), p_neckweak = max(d$P_NECKWEAK))
}

## phi, computed from a pattern -- the quantity the whole model turns on
phi_of <- function(pattern, rho = 0.5463,
                   w = c(0.20,0.24,0.13,0.08,0.12,0.06,0.10,0.07)) {
  rho * sum(w[pattern > 0])
}

## The bound itself, in closed form: no drug parameter appears.
tw_floor <- function(phi, p = list(SEVMAX=35, HS=1.6, L50=0.8347, PAINMAX=20,
                                   KPON=0.03, KPOFF=0.03, QP=1.5,
                                   DISMAX=30, GDIS=0.8, WDISSEV=0.55)) {
  L <- 1 - phi
  sev <- p$SEVMAX * L^p$HS / (L^p$HS + p$L50^p$HS)
  a <- p$KPON * L^p$QP
  pn <- a / (a + p$KPOFF)
  dis <- p$DISMAX * p$GDIS * (p$WDISSEV * sev / p$SEVMAX + (1 - p$WDISSEV) * pn)
  c(severity = sev, pain = p$PAINMAX * pn, disability = dis,
    total = sev + dis + p$PAINMAX * pn)
}

## Convenience accessor: mrgsolve exposes the $ENV block via env_get().
model_env <- function(mod) {
  if (exists("env_get", asNamespace("mrgsolve"))) return(mrgsolve::env_get(mod))
  mod@envir
}

run_all_scenarios <- function(mod, delta = 1) {
  res <- lapply(names(SCENARIOS), function(nm) {
    out <- sim_scenario(mod, nm, delta = delta)
    nd <- nadir_tw(out)
    data.frame(scenario = nm, label = SCENARIOS[[nm]]$label,
               phi = phi_of(SCENARIOS[[nm]]$pattern),
               wk4 = tw_at(out, 28) - 42.94,
               wk12 = tw_at(out, 84) - 42.94,
               nadir = nd$tw - 42.94, nadir_day = nd$day,
               duration_d = duration_of_benefit(out))
  })
  do.call(rbind, res)
}

## =====================================================================
## CALIBRATION NOTES
## =====================================================================
## Seven parameters were fitted to seven reported clinical quantities.  Every
## other number in this file is a literature-shaped order-of-magnitude choice,
## and every result quoted in README.md is a PREDICTION of these seven.
##
##   k_cl, k_LC, rho  <- TWSTRS-total change from baseline at week 4 (-10.5) and
##                       week 12 (-6.5) with 240 U, and the nadir at 480 U
##                       (-13.0), which is the dose-response PLATEAU
##   dys_d50, dys_k   <- dysphagia incidence, 240 U (11%) and 120 U (5.5%)
##   k_b              <- 5-year neutralising-antibody rate on onabotulinumtoxinA
##                       q12wk (about 2%)
##   nw_d50           <- neck-weakness incidence at 240 U (about 9%)
##
## rho deserves its own note.  A dose-response curve that flattens while
## adverse effects keep climbing is not pharmacological saturation -- the toxin
## has not run out of anything.  It is the NEEDLE running out of disease to
## reach.  Reading the plateau that way turns it from a nuisance into a
## measurement, and the measurement says that a standard surface injection
## affects only about 36% of the dystonic drive.
##
## Held back from the fit, as checks: the 60 U and 120 U points; the timing of
## peak effect; patient-perceived duration; and the continued rise of dysphagia
## above 480 U.  The first and last are reproduced.  The middle two are NOT --
## the model's effect-time curve is too square, and that failure is reported in
## A13(2) rather than fitted away.
##
## FIVE PLACES THIS MODEL DOES NOT FIT are enumerated in A13 of the Python
## port.  The most important is A13(1): the plateau identifies the PRODUCT
## phi = rho * Sigma-w, but NOT its factorisation.  The split between "better
## placement" and "longer target list" therefore rests on the assumed torque
## shares W1..W8 and is the weakest number in the whole construction.  The
## TOTAL non-dose headroom does not depend on that split, and it is roughly
## five times the dose headroom.
## =====================================================================
