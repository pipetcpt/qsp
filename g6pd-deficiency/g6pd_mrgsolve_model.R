## =====================================================================
##  G6PD DEFICIENCY — drug- and fava-induced acute haemolysis
##  포도당-6-인산 탈수소효소 결핍증 · QSP / PK-PD model for mrgsolve
##
##  48 ODE compartments:
##      12 drug PK      (primaquine + its oxidant equivalent, tafenoquine,
##                       dapsone, rasburicase + urate, methylene blue,
##                       fava divicine — gut and central)
##      12 RBC age bins (R1..R12, 10 days apiece, a 120-day Erlang chain)
##      12 damage bins  (Z1..Z12, hemichrome / Heinz-body burden per age bin)
##       4 erythropoiesis (3 marrow transit compartments + blood retics)
##       8 read-outs    (methaemoglobin, unconjugated bilirubin, free plasma
##                       haemoglobin, haptoglobin, tubular injury, creatinine,
##                       cumulative lysis, cumulative filtered haemoglobin)
##
##  TIME UNIT: DAYS.
##  UNITS: red cells 10^12/L whole blood · Hb g/dL · drug mg and mg/L ·
##         oxidant flux mmol H2O2 per L red cells per day · GSH-equivalent
##         damage mmol/L · bilirubin mg/dL · free Hb g/L plasma ·
##         haptoglobin g/L · creatinine mg/dL · urate mmol.
##
## =====================================================================
##
##  THE THESIS THIS MODEL IS BUILT AROUND
##  --------------------------------------
##  G6PD deficiency is normally written down as a number: "activity is
##  <10% of normal." This model asserts that the number is the wrong
##  object, that the disease is a FUNCTION, and that essentially every
##  clinically important behaviour of the condition — including three
##  that look like unrelated trivia — is a consequence of that one
##  substitution.
##
##  (1) THE ENZYME IS A FUNCTION OF CELL AGE, AND THE VARIANT IS ITS SLOPE.
##      A mature erythrocyte has no nucleus and no ribosome. It cannot
##      make G6PD. The enzyme is loaded once, in the erythroblast, and
##      from then on it only decays:
##
##            E(a) = E0 * exp(-ln2 * a / TAU)
##
##      TAU is the variant. That is the entire genotype in this model:
##
##            B (wild type)    E0 1.00   TAU 62 d
##            A- (Class III)   E0 0.55   TAU 13 d
##            Mediterranean    E0 0.05   TAU 12 d
##
##      Note what this does to a 115-day-old cell: it keeps 28% of its
##      enzyme if it is wild type, 0.1% if it is A-, and 0.006% if it is
##      Mediterranean. The "10% of normal" that the laboratory reports is
##      an AVERAGE OVER THIS CURVE, and averages hide slopes.
##
##  (2) HAEMOLYSIS IS A THRESHOLD SWEEPING AN AGE DISTRIBUTION.
##      A cell survives an oxidant load OX if its NADPH regeneration
##      capacity exceeds it: VMAXOX * E(a) > OX. Because E(a) is a clean
##      exponential, the critical age is available in closed form and the
##      model reports it at every time point as ASTAR:
##
##            a* = (TAU / ln2) * ln( E0 * VMAXOX / OXeff )
##
##      where OXeff = OX - KPIT*HZ50 is the flux net of the small burden the
##      cell can carry without triggering removal. Note what the formula is
##      NOT: there is no fitted "haemolysis rate" in it. It contains the
##      genotype (E0, TAU), a biochemical capacity (VMAXOX) and the exposure
##      (OX), and nothing else.
##
##      Cells older than a* cannot hold glutathione, oxidise their
##      haemoglobin, precipitate Heinz bodies and are removed. Haemolysis
##      is therefore not a rate constant multiplying a pool; it is a
##      KNIFE MOVING THROUGH AN AGE HISTOGRAM. Everything below follows.
##
##  (3) WHY THE SAME DRUG SELF-LIMITS IN ONE PATIENT AND NOT THE OTHER.
##      This is the model's headline result and it is scenario 1 vs 2.
##      Give primaquine 30 mg daily and hold it:
##
##        A- : the age-activity curve is STEEP, so a* lands INSIDE the
##             distribution. The oldest cells die; every survivor is young
##             and enzyme-replete; reticulocytes arrive on day 5-8 carrying
##             a FULL enzyme load and reset the age clock. Haemoglobin
##             falls, bottoms out around day 7-8, and then CLIMBS BACK
##             TOWARDS BASELINE WHILE THE DRUG IS STILL BEING TAKEN.
##             This is Dern's 1954 experiment, and the model reproduces it
##             without a single parameter that knows about it.
##
##        Mediterranean : the curve is FLAT AND LOW. a* lands far down in
##             the distribution, so the survivor pool is thin, the cells the
##             marrow supplies are themselves deficient and reach a* within
##             weeks, and a large fraction of the removal is so far past the
##             cell\'s remaining buffering capacity that it happens IN THE
##             CIRCULATION rather than in the spleen. The nadir is roughly
##             twice as deep, haemoglobinuria appears, and the plateau the
##             patient ends on is bought with near-maximal erythropoiesis
##             rather than with recovery.
##
##      Two patients, one drug, one dose, both labelled "severely
##      deficient", and the difference between "keep going, monitor" and
##      "stop and transfuse" is a single decay constant.
##
##  (4) THE HALF-LIFE IS THE TOXICOLOGY (scenarios 3 and 4).
##      Tafenoquine and primaquine are the same chemistry doing the same
##      damage. What differs is that primaquine's half-life is ~7 h and
##      tafenoquine's is ~15 DAYS. The reticulocyte rescue in (3) needs
##      4-7 days to arrive. So a fast-off drug delivers its oxidant insult
##      into a system that can answer it and can be withdrawn if it cannot;
##      a slow-off drug outruns the rescue and cannot be withdrawn at all.
##      Deliver the SAME total exposure both ways and the nadirs differ.
##      This is the mechanistic content of the FDA's decision to demand
##      QUANTITATIVE G6PD testing (>70%) before tafenoquine while allowing
##      weekly primaquine in mild deficiency with no test at all.
##
##  (5) WEEKLY DOSING IS NOT A SMALLER DOSE, IT IS A DIFFERENT MECHANISM
##      (scenario 5). WHO recommends primaquine 45 mg ONCE WEEKLY x 8
##      weeks for radical cure in G6PD-deficient P. vivax. In this model
##      the weekly pulse drives a* down hard for a few hours, kills the
##      thin oldest slice of the histogram, and then the 7-day gap lets
##      reticulocytes refill it before the next pulse. The haemoglobin
##      trace is a sawtooth around a stable mean instead of a slide. The
##      total dose is not much lower; the TIME STRUCTURE is what is doing
##      the work.
##
##  (6) ONE NADPH POOL, TWO CONSUMERS — SO METHYLENE BLUE IS A DERIVATION
##      (scenario 7). NADPH feeds glutathione reductase AND it feeds the
##      NADPH-methaemoglobin reductase that converts methylene blue to
##      leucomethylene blue, which is the species that actually reduces
##      Fe3+ back to Fe2+. The model therefore multiplies methylene blue's
##      efficacy by the same NADPH availability that governs haemolysis,
##      and adds unreduced methylene blue to the oxidant flux. Give
##      dapsone until methaemoglobin is 10%, then give the antidote: in a
##      normal subject methaemoglobin collapses within the hour; in a
##      deficient one it barely moves AND the haemolysis gets worse. The
##      contraindication is not memorised here, it is computed.
##
##  (7) THE OXIDANT DOSE OF RASBURICASE IS SET BY THE TUMOUR, NOT THE
##      PRESCRIPTION (scenario 8). Urate oxidase runs
##      urate + O2 + H2O -> allantoin + H2O2, one mole for one mole. So
##      the peroxide burden is stoichiometric in the urate pool being
##      destroyed. A patient with a big tumour lysis load receives a
##      bigger oxidant dose than a patient with a small one at exactly the
##      same 0.2 mg/kg. This is why the sensible manoeuvre when G6PD
##      status is unknown is not dose reduction but allopurinol, which
##      generates no peroxide at all.
##
##  (8) THE DIAGNOSTIC TEST LIES EXACTLY WHEN YOU USE IT (scenario 9).
##      The quantitative assay reports the age-weighted MEAN of E(a) over
##      the circulating population. Haemolysis has just deleted the old,
##      enzyme-poor cells and replaced them with reticulocytes carrying
##      1.5x a young cell's activity. So the assay READING RISES DURING
##      AND AFTER A HAEMOLYTIC EPISODE — the model outputs G6PDPCT so you
##      can watch a genuinely A- patient test in the normal range at the
##      nadir. Retest at three months, one full red cell lifespan.
##
##  (9) NEONATAL KERNICTERUS IS A PRODUCT, NOT A SUM (scenario 10).
##      Serum bilirubin is production over conjugation. G6PD deficiency
##      raises the numerator; the UGT1A1 (TA)7/7 Gilbert promoter lowers
##      the denominator. Each alone leaves the neonate in the safe range.
##      The model multiplies them, and only the double carrier crosses the
##      exchange-transfusion line. G6PD deficiency is the leading
##      identified cause of kernicterus worldwide and this is the shape of
##      the reason.
##
##  ---------------------------------------------------------------------
##  WHAT THIS MODEL IS NOT
##  ---------------------------------------------------------------------
##  * It is a mid-fidelity teaching / hypothesis-generating model, not a
##    validated clinical tool. Parameters are anchored to published
##    central values (see g6pd_references.md) but the model as a whole has
##    NOT been fitted to individual patient data.
##  * The age structure is an Erlang chain of 12 compartments, so a cell's
##    age within a bin is smeared (CV about 29% on the 120-day lifespan)
##    rather than sharp. This softens a* into a narrow band instead of a
##    knife edge, which is closer to reality anyway, but it does mean the
##    model cannot resolve age differences finer than about 10 days.
##  * Heterozygous females are TWO red cell populations, not one. Rather
##    than double the state vector the model is run twice and mixed —
##    see mix_mosaic() below. That is exact for everything except the
##    shared plasma read-outs (bilirubin, free Hb), which the helper
##    recombines by mass.
##  * Glutathione is treated as a quasi-steady-state algebraic function of
##    the oxidant load and the bin's NADPH capacity. GSH turnover is on
##    the order of minutes against a disease that runs over days, so this
##    is a safe reduction, but it means the model cannot say anything
##    about the first few minutes of an exposure.
## =====================================================================

library(mrgsolve)
suppressPackageStartupMessages({
  library(dplyr)
})

## =====================================================================
##  1. THE MODEL
## =====================================================================

g6pd_code <- '
$PROB
G6PD deficiency: age-structured erythrocyte redox model with drug PK/PD.

$GLOBAL
#define NBIN 12
// per-age-bin quantities, filled in $MAIN (they depend only on parameters)
double abin_[NBIN];   // nominal mean age of the bin (days)
double Ebin_[NBIN];   // G6PD activity, fraction of a normal young cell
double Vbin_[NBIN];   // peroxide disposal capacity, mmol/L RBC/day
// loop-body scratch. mrgsolve hoists every "double" to file scope, so each
// name must be declared exactly once in the whole model; declaring them here
// and only ASSIGNING inside the blocks keeps that rule easy to hold.
double zi, unbuf, zev, kev, omeg, omn, kiv, rk, inflow;
double R_[NBIN], Z_[NBIN], dR_[NBIN], dZ_[NBIN];
#define SAFEPOS(x) ((x) > 0.0 ? (x) : 0.0)

$PARAM @annotated
// ---- THE VARIANT: these two numbers are the genotype -----------------
E0     :  1.00 : G6PD activity of a brand new red cell (fraction of normal)
TAU    : 62.0  : in-vivo enzyme half-life inside the cell (d) -- THE SLOPE
FRET   :  1.50 : reticulocyte activity relative to a young mature cell

// ---- red cell age structure ------------------------------------------
LSPAN  : 120.0 : nominal red cell lifespan (d), = NBIN * bin width
TRET   :  1.20 : reticulocyte maturation time in blood (d)
MCH    : 29.70 : haemoglobin per cell (pg)
BV     :  5.00 : blood volume (L)
VPLAS  :  3.00 : plasma volume (L)
VRBCL  :  2.00 : packed red cell volume (L)

// ---- redox capacity and damage ---------------------------------------
VMAXOX : 60.0  : max peroxide disposal by a normal young cell (mmol/L RBC/d)
OXBASE :  0.060: baseline oxidant load from Hb autoxidation (mmol/L RBC/d)
KHZ    :  1.00 : hemichrome formed per unit unbuffered peroxide
KPIT   :  0.50 : splenic pitting / repair of Heinz burden (1/d)
KEVMAX :  0.80 : max extravascular removal rate (1/d)
HZ50   :  0.10 : Heinz burden for half-max extravascular removal
NEV    :  3.00 : Hill coefficient, extravascular
KIVMAX :  8.00 : max intravascular lysis rate (1/d)
OMEG50 : 20.00 : overwhelm ratio OX/(residual capacity) for half-max in-situ lysis
NIV    :  3.00 : Hill coefficient, intravascular
OXREF  :  0.02 : floor on residual buffering capacity (mmol/L RBC/d)
SPLEEN :  1.00 : spleen function multiplier (0 = splenectomised)

// ---- erythropoiesis ---------------------------------------------------
HB0    : 15.00 : haemoglobin set-point the kidney defends (g/dL)
EPO0   : 10.00 : baseline erythropoietin (mIU/mL)
GEPO   :  5.40 : log-linear EPO gain on the Hb deficit
EPOMAX : 6000. : ceiling on EPO (mIU/mL)
PROD0  :  0.041667 : baseline red cell production (10^12/L/d)
EMAXERY:  5.00 : maximum fold INCREASE in production above baseline
EC50ERY: 60.00 : EPO above baseline giving half-max expansion (mIU/mL)
TTM0   :  5.00 : marrow transit time at baseline (d)
FSHIFT :  0.60 : max fractional shortening of marrow transit
FERY   :  1.00 : marrow capacity multiplier (iron/folate/B12; 0 = aplasia)

// ---- primaquine -------------------------------------------------------
FPQ    :  0.96 : oral bioavailability
KAPQ   : 36.00 : absorption rate (1/d)
VPQ    : 240.0 : apparent volume of distribution (L)
CLPQ   : 576.0 : apparent clearance (L/d) -> t1/2 ~7 h
KMPQ   :  2.10 : formation/loss rate of the oxidant equivalent (1/d)
AS2D6  :  1.00 : CYP2D6 activity score (0 = poor metaboliser)
SPQ    :  6.00 : oxidant flux per unit primaquine oxidant equivalent

// ---- tafenoquine ------------------------------------------------------
FTQ    :  0.90 : oral bioavailability
KATQ   : 12.00 : absorption rate (1/d)
VTQ    : 1600. : volume of distribution (L)
CLTQ   : 73.90 : clearance (L/d) -> t1/2 ~15 d
STQ    :  6.00 : oxidant flux per mg/L tafenoquine

// ---- dapsone ----------------------------------------------------------
FDP    :  0.90 : oral bioavailability
KADP   : 12.00 : absorption rate (1/d)
VDP    : 70.00 : volume of distribution (L)
CLDP   : 46.60 : clearance (L/d) -> t1/2 ~25 h
FMNOH  :  0.05 : fraction present as the hydroxylamine at steady state
FNAT2  :  1.00 : NAT2 multiplier (1.0 fast, 1.5 slow acetylator)
SDAP   :  1.50 : oxidant flux per mg/L dapsone hydroxylamine
KMETDAP: 1000. : methaemoglobin formation per mg/L hydroxylamine (%/d)

// ---- rasburicase and urate --------------------------------------------
VRBX   :  6.00 : rasburicase volume of distribution (L)
CLRBX  :  4.75 : rasburicase clearance (L/d) -> t1/2 ~21 h
KCATRBX:  3.30 : urate turnover per mg/L rasburicase (1/d)
KMUR   :  0.10 : Michaelis constant for urate (mmol/L)
VUR    : 35.00 : urate volume of distribution (L)
KGENUR :  3.00 : urate generation (mmol/d); tumour lysis raises this
CLUR   :  9.00 : renal urate clearance (L/d)
SRBX   :  0.005: fraction of generated H2O2 reaching the red cell

// ---- methylene blue ----------------------------------------------------
VMB    : 1400. : volume of distribution (L)
CLMB   : 4400. : clearance (L/d) -> t1/2 ~5 h
KMBRED : 540.0 : MetHb reduction per mg/L leuco-MB (1/d)
SMB    :  8.00 : oxidant flux per mg/L of UNREDUCED methylene blue

// ---- fava beans (divicine / isouramil) ---------------------------------
FFV    :  0.30 : fraction of the aglycone absorbed
KAFV   : 12.00 : absorption rate (1/d)
VFV    : 40.00 : volume of distribution (L)
CLFV   : 333.0 : clearance (L/d) -> t1/2 ~2 h
SFV    :  0.80 : oxidant flux per mg/L divicine

// ---- infection / other non-drug oxidant stress -------------------------
INFON  :  0.00 : infection switch (0/1)
OXINF  :  0.15 : oxidant flux added by systemic infection

// ---- methaemoglobin -----------------------------------------------------
KCYB5  : 12.00 : cytochrome b5 reductase, NADH-dependent (1/d)
FMETBAS: 12.00 : baseline MetHb formation (%/d) -> 1% at steady state
KMETOX : 73.00 : MetHb formation per unit oxidant flux the cell cannot reduce

// ---- bilirubin ----------------------------------------------------------
BILHB  : 34.00 : mg bilirubin produced per g haemoglobin catabolised
VBIL   : 12.00 : bilirubin distribution volume (L)
CLBIL  :  2.92 : adult UGT1A1 conjugation rate constant (1/d)
FUGT   :  1.00 : UGT1A1 genotype multiplier (6/6 1.0, 6/7 0.7, 7/7 0.35)
NEO    :  0.00 : neonatal mode switch (0/1)
UGTMAT0:  0.03 : neonatal UGT1A1 activity at birth (fraction of adult)
KUGTMAT:  0.025: postnatal UGT1A1 maturation rate (1/d) -> t1/2 ~28 d
PNA0   :  0.00 : postnatal age at t = 0 (d)
KENT   :  0.00 : enterohepatic recirculation input (mg/dL/d)
PHOTO  :  0.00 : phototherapy clearance added to CLBIL (1/d)
ALBM   :  0.60 : albumin (mM); 0.60 adult, ~0.45 term neonate
KABIL  : 7.0e7 : bilirubin-albumin association constant (1/M)

// ---- intravascular haemolysis handling -----------------------------------
HP0    :  1.30 : baseline haptoglobin (g/L)
KSYNHP :  0.20 : haptoglobin turnover (1/d)
KBINDHP: 200.0 : Hb-haptoglobin association (L/g/d)
STOIHP :  1.50 : g haptoglobin consumed per g haemoglobin bound
CLRENHB: 60.00 : clearance of unbound plasma Hb (1/d) -> t1/2 ~17 min
KTUBI  :  0.005: tubular injury per g/L/d of filtered Hb
KTUBR  :  0.35 : tubular recovery (1/d)
KTUB50 :  1.00 : tubular injury giving 50% GFR loss
CRGEN  : 1000. : creatinine generation (mg/d)
CLCR   : 100.0 : creatinine clearance at normal GFR (L/d)
VCRD   : 42.00 : creatinine distribution volume (L)
LDH0   : 180.0 : baseline LDH (U/L)
KLDH   :  0.90 : LDH rise per unit relative lysis rate

$CMT @annotated
// -- drug PK
PQGUT  : primaquine in gut (mg)
PQC    : primaquine central (mg)
PQOXF  : primaquine oxidant equivalent (mg/L, kinetically filtered)
TQGUT  : tafenoquine in gut (mg)
TQC    : tafenoquine central (mg)
DPGUT  : dapsone in gut (mg)
DPC    : dapsone central (mg)
RBXC   : rasburicase central (mg)
URATE  : urate pool (mmol)
MBC    : methylene blue central (mg)
FVGUT  : fava aglycone in gut (mg)
FVC    : divicine central (mg)
// -- red cells, oldest last (each bin is LSPAN/12 = 10 days wide)
R1     : red cells age 0-10 d   (10^12/L)
R2     : red cells age 10-20 d  (10^12/L)
R3     : red cells age 20-30 d  (10^12/L)
R4     : red cells age 30-40 d  (10^12/L)
R5     : red cells age 40-50 d  (10^12/L)
R6     : red cells age 50-60 d  (10^12/L)
R7     : red cells age 60-70 d  (10^12/L)
R8     : red cells age 70-80 d  (10^12/L)
R9     : red cells age 80-90 d  (10^12/L)
R10    : red cells age 90-100 d (10^12/L)
R11    : red cells age 100-110 d(10^12/L)
R12    : red cells age 110-120 d(10^12/L)
// -- oxidative damage carried by each age bin
Z1     : Heinz/hemichrome burden, bin 1 (mmol/L)
Z2     : Heinz/hemichrome burden, bin 2
Z3     : Heinz/hemichrome burden, bin 3
Z4     : Heinz/hemichrome burden, bin 4
Z5     : Heinz/hemichrome burden, bin 5
Z6     : Heinz/hemichrome burden, bin 6
Z7     : Heinz/hemichrome burden, bin 7
Z8     : Heinz/hemichrome burden, bin 8
Z9     : Heinz/hemichrome burden, bin 9
Z10    : Heinz/hemichrome burden, bin 10
Z11    : Heinz/hemichrome burden, bin 11
Z12    : Heinz/hemichrome burden, bin 12
// -- erythropoiesis
M1     : marrow erythroid transit 1 (10^12/L)
M2     : marrow erythroid transit 2 (10^12/L)
M3     : marrow erythroid transit 3 (10^12/L)
RETB   : blood reticulocytes (10^12/L)
// -- plasma read-outs
METHB  : methaemoglobin (% of total Hb)
UCB    : unconjugated bilirubin (mg/dL)
FREEHB : free plasma haemoglobin (g/L)
HPTG   : haptoglobin (g/L)
TUB    : tubular injury burden (arbitrary)
CREA   : serum creatinine (mg/dL)
CUMLYS : cumulative cells lysed (10^12/L)
HBURIN : cumulative filtered haemoglobin (g/L equivalent)

$MAIN
// bioavailability on the oral depots
F_PQGUT = FPQ;
F_TQGUT = FTQ;
F_DPGUT = FDP;
F_FVGUT = FFV;

// ---------------------------------------------------------------------
// THE ONE LINE. Enzyme activity as a function of the bin nominal age,
// and the peroxide disposal capacity that follows from it.
// ---------------------------------------------------------------------
double wbin = LSPAN / (double) NBIN;
for (int i = 0; i < NBIN; ++i) {
  abin_[i] = wbin * ((double) i + 0.5);
  Ebin_[i] = E0 * exp(-log(2.0) * abin_[i] / TAU);
  Vbin_[i] = VMAXOX * Ebin_[i];
}

// steady-state initial condition: a flat age histogram, 1% reticulocytes
double rbin0 = PROD0 / (1.0 / wbin);
R1_0 = rbin0; R2_0 = rbin0; R3_0  = rbin0; R4_0  = rbin0;
R5_0 = rbin0; R6_0 = rbin0; R7_0  = rbin0; R8_0  = rbin0;
R9_0 = rbin0; R10_0 = rbin0; R11_0 = rbin0; R12_0 = rbin0;
M1_0 = PROD0 * TTM0 / 3.0;
M2_0 = PROD0 * TTM0 / 3.0;
M3_0 = PROD0 * TTM0 / 3.0;
RETB_0  = PROD0 * TRET;
METHB_0 = FMETBAS / KCYB5;
UCB_0   = 0.60;
HPTG_0  = HP0;
CREA_0  = 1.00;

$ODE
// ---------------------------------------------------------------------
// 0. gather the age-structured states
// ---------------------------------------------------------------------
R_[0]=R1; R_[1]=R2; R_[2]=R3;  R_[3]=R4;  R_[4]=R5;  R_[5]=R6;
R_[6]=R7; R_[7]=R8; R_[8]=R9;  R_[9]=R10; R_[10]=R11; R_[11]=R12;
Z_[0]=Z1; Z_[1]=Z2; Z_[2]=Z3;  Z_[3]=Z4;  Z_[4]=Z5;  Z_[5]=Z6;
Z_[6]=Z7; Z_[7]=Z8; Z_[8]=Z9;  Z_[9]=Z10; Z_[10]=Z11; Z_[11]=Z12;

double Rmat = 0.0;
double Esum = 0.0;
for (int j = 0; j < NBIN; ++j) {
  Rmat += SAFEPOS(R_[j]);
  Esum += SAFEPOS(R_[j]) * Ebin_[j];
}
double Rret = SAFEPOS(RETB);
double Rtot = Rmat + Rret;
Esum += Rret * E0 * FRET;

// The population-mean activity. This is what the laboratory assay
// measures, and it is ALSO the right scalar for whole-body NADPH-
// dependent processes such as methylene blue reduction.
double meanE = (Rtot > 1e-9) ? (Esum / Rtot) : 0.0;
// EREF is the same quantity in a normal subject at steady state: the mean of
// exp(-ln2*a/62) over a flat 0-120 d histogram plus 1% reticulocytes at 1.5x.
// Normalising by it makes NADPHAV read 1.00 in a normal subject, 0.16 in A-
// and 0.014 in Mediterranean -- i.e. the published activity classes.
#define EREF 0.5593
double NADPHAV = meanE / EREF;
if (NADPHAV > 1.0) NADPHAV = 1.0;

// ---------------------------------------------------------------------
// 1. drug concentrations
// ---------------------------------------------------------------------
double CPQ  = PQC  / VPQ;
double CTQ  = TQC  / VTQ;
double CDP  = DPC  / VDP;
double CRBX = RBXC / VRBX;
double CMB  = MBC  / VMB;
double CFV  = FVC  / VFV;
double CNHOH = FNAT2 * FMNOH * CDP;         // dapsone hydroxylamine, QSS

// urate oxidase runs one mole of peroxide per mole of urate destroyed
double urc   = SAFEPOS(URATE) / VUR;
double v_uox = KCATRBX * CRBX * urc / (KMUR + urc);   // mmol/L/d
double h2o2_mmol_d = v_uox * VUR;                     // mmol/d, whole body

// methylene blue only reduces MetHb through NADPH it cannot get here;
// whatever is not reduced behaves as a fresh oxidant.
double mb_reduced   = CMB * NADPHAV;
double mb_unreduced = CMB * (1.0 - NADPHAV);

// ---------------------------------------------------------------------
// 2. TOTAL OXIDANT FLUX -- every trigger enters the model right here
// ---------------------------------------------------------------------
double OX = OXBASE
          + SPQ  * SAFEPOS(PQOXF)
          + STQ  * CTQ
          + SDAP * CNHOH
          + SRBX * h2o2_mmol_d / VRBCL
          + SMB  * mb_unreduced
          + SFV  * CFV
          + OXINF * INFON;

// ---------------------------------------------------------------------
// 3. per-age-bin redox balance, damage and removal
//    A bin survives if VMAXOX * E(a) > OX. Everything above that line
//    is unbuffered peroxide and goes into haemoglobin.
// ---------------------------------------------------------------------
double kage = (double) NBIN / LSPAN;
double lysEV = 0.0;
double lysIV = 0.0;
double hz50n = pow(HZ50, NEV);
double omegn = pow(OMEG50, NIV);

for (int k = 0; k < NBIN; ++k) {
  zi    = SAFEPOS(Z_[k]);
  unbuf = SAFEPOS(OX - Vbin_[k]);
  dZ_[k] = KHZ * unbuf - KPIT * zi;

  zev = pow(zi, NEV);
  kev = KEVMAX * SPLEEN * zev / (hz50n + zev);
  // Intravascular lysis is not "more damage" -- it is being OVERWHELMED.
  // The governing quantity is the RATIO of the load to whatever buffering
  // the cell still has, so a cell whose enzyme is essentially zero ruptures
  // in the circulation while a merely old cell is pitted by the spleen.
  // This is what separates Mediterranean (haemoglobinuria) from A- (not).
  omeg = OX / (Vbin_[k] + OXREF);
  omn  = pow(omeg, NIV);
  kiv  = (zi > HZ50) ? (KIVMAX * omn / (omegn + omn)) : 0.0;

  rk = SAFEPOS(R_[k]);
  lysEV += kev * rk;
  lysIV += kiv * rk;

  inflow = (k == 0) ? (Rret / TRET) : (kage * SAFEPOS(R_[k-1]));
  dR_[k] = inflow - kage * rk - (kev + kiv) * rk;
}
double senesce = kage * SAFEPOS(R_[NBIN-1]);   // normal end-of-life removal

dxdt_R1  = dR_[0];  dxdt_R2  = dR_[1];  dxdt_R3  = dR_[2];  dxdt_R4  = dR_[3];
dxdt_R5  = dR_[4];  dxdt_R6  = dR_[5];  dxdt_R7  = dR_[6];  dxdt_R8  = dR_[7];
dxdt_R9  = dR_[8];  dxdt_R10 = dR_[9];  dxdt_R11 = dR_[10]; dxdt_R12 = dR_[11];
dxdt_Z1  = dZ_[0];  dxdt_Z2  = dZ_[1];  dxdt_Z3  = dZ_[2];  dxdt_Z4  = dZ_[3];
dxdt_Z5  = dZ_[4];  dxdt_Z6  = dZ_[5];  dxdt_Z7  = dZ_[6];  dxdt_Z8  = dZ_[7];
dxdt_Z9  = dZ_[8];  dxdt_Z10 = dZ_[9];  dxdt_Z11 = dZ_[10]; dxdt_Z12 = dZ_[11];

// ---------------------------------------------------------------------
// 4. erythropoiesis -- the arm that decides whether it all self-limits
// ---------------------------------------------------------------------
double Hb  = Rtot * MCH / 10.0;
double epo = EPO0 * pow(HB0 / (Hb > 1.0 ? Hb : 1.0), GEPO);
if (epo > EPOMAX) epo = EPOMAX;
double dEPO  = SAFEPOS(epo - EPO0);
double drive = dEPO / (dEPO + EC50ERY);
double prod  = PROD0 * FERY * (1.0 + EMAXERY * drive);
// above the set-point the kidney switches erythropoiesis DOWN, not merely
// back to baseline -- without this the rebound overshoots and never returns
if (epo < EPO0) prod = PROD0 * FERY * (epo / EPO0);
double ttm   = TTM0 * (1.0 - FSHIFT * drive);
double kM    = 3.0 / ttm;
double tretb = TRET * (1.0 + drive);        // shift reticulocytes mature slower

dxdt_M1 = prod  - kM * M1;
dxdt_M2 = kM*M1 - kM * M2;
dxdt_M3 = kM*M2 - kM * M3;
dxdt_RETB = kM * M3 - Rret / tretb;

// ---------------------------------------------------------------------
// 5. methaemoglobin -- the SECOND consumer of the same NADPH pool
// ---------------------------------------------------------------------
double kred  = KCYB5 + KMBRED * mb_reduced;
// Formation acts on the ferrous haemoglobin that is LEFT, so MetHb cannot
// exceed 100%. The third term is the reason rasburicase and favism cause
// methaemoglobinaemia and not just haemolysis: an oxidant flux the cell has
// no NADPH to answer oxidises the iron as well as the thiols.
double mform = FMETBAS + KMETDAP * CNHOH + KMETOX * OX * (1.0 - NADPHAV);
dxdt_METHB = mform * (1.0 - SAFEPOS(METHB) / 100.0) - kred * SAFEPOS(METHB);

// ---------------------------------------------------------------------
// 6. bilirubin -- production over conjugation, two independent hits
// ---------------------------------------------------------------------
double hb_lysed_g_d = (lysEV + lysIV + senesce) * BV * MCH;   // g Hb/d
double bilprod      = BILHB * hb_lysed_g_d;                   // mg/d
double ugtmat = 1.0;
if (NEO > 0.5) {
  ugtmat = UGTMAT0 + (1.0 - UGTMAT0) * (1.0 - exp(-KUGTMAT * (PNA0 + SOLVERTIME)));
}
double clb = CLBIL * FUGT * ugtmat + PHOTO;
dxdt_UCB = bilprod / (VBIL * 10.0) + KENT - clb * SAFEPOS(UCB);

// ---------------------------------------------------------------------
// 7. intravascular haemolysis: free Hb, haptoglobin, tubule, creatinine
// ---------------------------------------------------------------------
double iv_hb_g_d = lysIV * BV * MCH;               // g Hb/d released free
double fhb = SAFEPOS(FREEHB);
double hpg = SAFEPOS(HPTG);
double bind = KBINDHP * hpg * fhb;
double filt = CLRENHB * fhb;
dxdt_FREEHB = iv_hb_g_d / VPLAS - bind - filt;
dxdt_HPTG   = KSYNHP * (HP0 - hpg) - STOIHP * bind;
dxdt_HBURIN = filt;
dxdt_TUB    = KTUBI * filt - KTUBR * SAFEPOS(TUB);
double gfrf = 1.0 / (1.0 + SAFEPOS(TUB) / KTUB50);
dxdt_CREA   = (CRGEN / 10.0 - CLCR * gfrf * SAFEPOS(CREA)) / VCRD;

dxdt_CUMLYS = lysEV + lysIV;

// ---------------------------------------------------------------------
// 8. drug PK
// ---------------------------------------------------------------------
dxdt_PQGUT = -KAPQ * PQGUT;
dxdt_PQC   =  KAPQ * PQGUT - (CLPQ / VPQ) * PQC;
dxdt_PQOXF =  KMPQ * (AS2D6 * CPQ - SAFEPOS(PQOXF));
dxdt_TQGUT = -KATQ * TQGUT;
dxdt_TQC   =  KATQ * TQGUT - (CLTQ / VTQ) * TQC;
dxdt_DPGUT = -KADP * DPGUT;
dxdt_DPC   =  KADP * DPGUT - (CLDP / VDP) * DPC;
dxdt_RBXC  = -(CLRBX / VRBX) * RBXC;
dxdt_URATE =  KGENUR - (CLUR / VUR) * SAFEPOS(URATE) - v_uox * VUR;
dxdt_MBC   = -(CLMB / VMB) * MBC;
dxdt_FVGUT = -KAFV * FVGUT;
dxdt_FVC   =  KAFV * FVGUT - (CLFV / VFV) * FVC;

$TABLE
// ---- the age-structured state, summarised -----------------------------
double Rmat_o = R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12;
double Rtot_o = Rmat_o + RETB;
double HB     = Rtot_o * MCH / 10.0;
double RETPCT = (Rtot_o > 1e-9) ? 100.0 * RETB / Rtot_o : 0.0;

double Esum_o = R1*Ebin_[0] + R2*Ebin_[1] + R3*Ebin_[2]  + R4*Ebin_[3]
              + R5*Ebin_[4] + R6*Ebin_[5] + R7*Ebin_[6]  + R8*Ebin_[7]
              + R9*Ebin_[8] + R10*Ebin_[9]+ R11*Ebin_[10]+ R12*Ebin_[11]
              + RETB * E0 * FRET;
// What the quantitative laboratory assay actually reports: the age-weighted
// MEAN, as a percentage of a normal subject at steady state (0.5593 = mean of
// exp(-ln2*a/62) over a flat 0-120 d histogram plus 1% retics at 1.5x).
// This normalisation is not fitted, and it lands the untreated variants on
// their published activity classes: B 100%, A- 16%, Mediterranean 1.4%.
double MEANE   = (Rtot_o > 1e-9) ? Esum_o / Rtot_o : 0.0;
double G6PDPCT = 100.0 * MEANE / 0.5593;

// ---- oxidant flux and the critical age, both in closed form ------------
double CPQ_o  = PQC / VPQ;
double CTQ_o  = TQC / VTQ;
double CDP_o  = DPC / VDP;
double CNH_o  = FNAT2 * FMNOH * CDP_o;
double CMB_o  = MBC / VMB;
double CFV_o  = FVC / VFV;
double urc_o  = (URATE > 0 ? URATE : 0) / VUR;
double vuox_o = KCATRBX * (RBXC/VRBX) * urc_o / (KMUR + urc_o);
double OXFLUX = OXBASE + SPQ*(PQOXF > 0 ? PQOXF : 0) + STQ*CTQ_o + SDAP*CNH_o
              + SRBX * (vuox_o * VUR) / VRBCL + SMB * CMB_o * (1.0 - MEANE)
              + SFV * CFV_o + OXINF * INFON;

// a* = (TAU/ln2) * ln(E0 * VMAXOX / OX): the age above which a cell cannot
// hold glutathione. LSPAN means "no cell in the circulation is at risk".
// The removal Hill is half-maximal at a steady-state burden HZ50, which the
// bin reaches when its UNBUFFERED flux is KPIT*HZ50. So the flux that matters
// is OXFLUX minus that offset, not OXFLUX itself.
double OXEFF = OXFLUX - KPIT * HZ50;
double ASTAR = LSPAN;
if (OXEFF <= 0.0) {
  ASTAR = LSPAN;                       // no circulating cell is at risk
} else if (E0 * VMAXOX <= OXEFF) {
  ASTAR = 0.0;                         // not even a reticulocyte is safe
} else {
  ASTAR = (TAU / log(2.0)) * log(E0 * VMAXOX / OXEFF);
  if (ASTAR > LSPAN) ASTAR = LSPAN;
}
// fraction of the circulating red cell mass that sits above a*
double ATRISK = 100.0 * (1.0 - (ASTAR / LSPAN));

// ---- plasma chemistry ---------------------------------------------------
double TSB   = UCB;
// one-site bilirubin-albumin binding, solved exactly for the free species
double Btot  = TSB * 10.0 / 584.7 / 1000.0;   // mg/dL -> mol/L
double bq    = KABIL * (ALBM/1000.0 - Btot) + 1.0;
double BFREE = 1e9 * (-bq + sqrt(bq*bq + 4.0*KABIL*Btot)) / (2.0*KABIL);  // nM

double EPOOUT = EPO0 * pow(HB0 / (HB > 1.0 ? HB : 1.0), GEPO);
if (EPOOUT > EPOMAX) EPOOUT = EPOMAX;

double HPOUT  = HPTG;
double FHBOUT = FREEHB;
double HBURIA = (FREEHB > 0.10) ? 1.0 : 0.0;   // visibly dark urine
double GFRF   = 1.0 / (1.0 + (TUB > 0 ? TUB : 0) / KTUB50);
double METPCT = METHB;

$CAPTURE @annotated
HB      : haemoglobin (g/dL)
RETPCT  : reticulocytes (% of red cells)
G6PDPCT : G6PD activity AS THE ASSAY WOULD REPORT IT (% of normal)
OXFLUX  : total oxidant flux (mmol H2O2 / L RBC / d)
ASTAR   : critical cell age -- cells older than this lyse (d)
ATRISK  : % of the circulating red cell mass older than a*
METPCT  : methaemoglobin (%)
TSB     : total serum bilirubin (mg/dL)
BFREE   : free (albumin-unbound) bilirubin (nM)
HPOUT   : haptoglobin (g/L)
FHBOUT  : free plasma haemoglobin (g/L)
HBURIA  : haemoglobinuria present (0/1)
GFRF    : GFR as a fraction of normal
EPOOUT  : erythropoietin (mIU/mL)
CPQ_o   : primaquine (mg/L)
CTQ_o   : tafenoquine (mg/L)
CDP_o   : dapsone (mg/L)
'

g6pd_mod <- mcode("g6pd", g6pd_code)

## =====================================================================
##  2. VARIANTS -- the genotype is two numbers, E0 and TAU
## =====================================================================
##  E0  = activity of a brand-new red cell, as a fraction of normal
##  TAU = in-vivo half-life of the enzyme INSIDE the cell (days)
##
##  These are the model. Everything clinical is downstream of them.
##  Sources for the decay half-lives: Piomelli 1968 (A- 13 d vs B 62 d),
##  Morelli 1978, Luzzatto & Arese 2018 review. See g6pd_references.md.

VARIANTS <- list(
  ## WHO class IV -- not deficient
  normal        = list(E0 = 1.00, TAU = 62, label = "B (wild type)"),
  Aplus         = list(E0 = 0.90, TAU = 45, label = "A+ (A376G)"),
  ## WHO class III -- 10-60% activity, self-limiting drug-induced haemolysis
  Aminus        = list(E0 = 0.55, TAU = 13, label = "A- (A376G+G202A)"),
  Mahidol       = list(E0 = 0.35, TAU = 20, label = "Mahidol (G487A)"),
  ## WHO class II -- <10% activity, severe, includes favism
  Mediterranean = list(E0 = 0.05, TAU = 12, label = "Mediterranean (C563T)"),
  Canton        = list(E0 = 0.08, TAU = 13, label = "Canton (G1376T)"),
  ## WHO class I -- chronic non-spherocytic haemolytic anaemia at baseline
  ClassI        = list(E0 = 0.02, TAU =  7, label = "Class I (CNSHA)")
)

#' Build a parameter list for one patient.
#'
#' @param variant  name from VARIANTS
#' @param wt       body weight (kg); scales the volumes that should scale
#' @param ugt      UGT1A1 promoter genotype: "6/6", "6/7" or "7/7"
#' @param cyp2d6   CYP2D6 activity score (0 poor .. 2 ultrarapid)
#' @param nat2     "fast" or "slow" acetylator
#' @param spleen   1 = normal spleen, 0 = splenectomised
#' @param neonate  TRUE to switch on immature UGT1A1 and neonatal volumes
#' @param pna      postnatal age in days at t = 0 (neonates)
g6pd_patient <- function(variant = "normal", wt = 70, ugt = "6/6",
                         cyp2d6 = 1.0, nat2 = "fast", spleen = 1,
                         marrow = 1.0, neonate = FALSE, pna = 0) {
  v <- VARIANTS[[variant]]
  if (is.null(v)) stop("unknown variant: ", variant)
  fugt <- switch(ugt, "6/6" = 1.00, "6/7" = 0.70, "7/7" = 0.35,
                 stop("ugt must be 6/6, 6/7 or 7/7"))
  p <- list(E0 = v$E0, TAU = v$TAU,
            AS2D6 = cyp2d6, FNAT2 = if (nat2 == "slow") 1.5 else 1.0,
            FUGT = fugt, SPLEEN = spleen, FERY = marrow,
            BV = 0.0714 * wt, VPLAS = 0.0429 * wt, VRBCL = 0.0286 * wt)
  if (neonate) {
    ## a 3 kg term newborn: smaller everything, immature conjugation,
    ## a higher Hb set-point, and enterohepatic recirculation switched on
    p <- modifyList(p, list(
      ## the newborn red cell lives in a more oxidising environment than the
      ## adult one (higher pO2 after the first breath, HbF, low antioxidant
      ## enzymes), which is why OXBASE is raised here: it is what makes a
      ## G6PD-deficient NEWBORN haemolyse without any drug at all.
      NEO = 1, PNA0 = pna, ALBM = 0.45,
      BV = 0.085 * wt, VPLAS = 0.048 * wt, VRBCL = 0.037 * wt,
      VBIL = 0.55 * wt, HB0 = 17.0, LSPAN = 80, MCH = 34,
      PROD0 = 0.041667 * (120 / 80), KENT = 0.8, OXBASE = 0.15,
      CRGEN = 1000 * wt / 70, CLCR = 100 * wt / 70, VCRD = 42 * wt / 70))
  }
  p
}

## =====================================================================
##  3. REGIMENS
## =====================================================================
##  Doses are built as mrgsolve event objects. Amounts in mg; the FAVA
##  amount is mg of absorbable divicine-equivalent aglycone (a 250 g meal
##  of fresh broad beans is roughly 1000 mg).

rx_none        <- function() ev(amt = 0, cmt = "PQGUT", time = 0)
rx_primaquine  <- function(dose = 30, days = 14, start = 0, ii = 1)
  ev(amt = dose, cmt = "PQGUT", time = start, ii = ii,
     addl = max(0, ceiling(days / ii) - 1))
rx_pq_weekly   <- function(dose = 45, weeks = 8, start = 0)
  ev(amt = dose, cmt = "PQGUT", time = start, ii = 7, addl = weeks - 1)
rx_tafenoquine <- function(dose = 300, start = 0)
  ev(amt = dose, cmt = "TQGUT", time = start)
rx_dapsone     <- function(dose = 100, days = 60, start = 0)
  ev(amt = dose, cmt = "DPGUT", time = start, ii = 1, addl = days - 1)
rx_rasburicase <- function(dose_mg = 14, days = 3, start = 0)
  ev(amt = dose_mg, cmt = "RBXC", time = start, ii = 1, addl = days - 1)
rx_methyleneblue <- function(dose = 140, start = 0)
  ev(amt = dose, cmt = "MBC", time = start)
rx_fava        <- function(amt = 1000, start = 0)
  ev(amt = amt, cmt = "FVGUT", time = start)

#' Run one scenario.
#'
#' Always simulates a long drug-free lead-in first so that every
#' trajectory starts from ITS OWN steady state. This matters: a Class I or
#' Mediterranean patient is already mildly anaemic before anything is given
#' because a* sits inside their age distribution at BASELINE, and endpoints
#' must be measured against that, not against 15 g/dL. Two full red cell
#' lifespans are needed for a shortened-survival variant to settle.
run_g6pd <- function(mod = g6pd_mod, patient = g6pd_patient(),
                     rx = NULL, days = 60, lead_in = 250,
                     extra = list(), name = "scenario", delta = 0.05) {
  p <- modifyList(patient, extra)
  e <- if (is.null(rx)) rx_none() else rx
  e <- mutate(as.data.frame(e), time = time + lead_in)
  out <- mod %>%
    param(p) %>%
    mrgsim(events = e, end = days + lead_in, delta = delta,
           maxsteps = 200000, atol = 1e-10, rtol = 1e-8) %>%
    as.data.frame()
  out$time <- out$time - lead_in
  out$scenario <- name
  out
}

#' Summarise a run against its own pre-exposure baseline.
summarise_g6pd <- function(out) {
  base <- out[which.min(abs(out$time - 0)), ]
  post <- out[out$time >= 0, ]
  nadir_i <- which.min(post$HB)
  data.frame(
    scenario   = out$scenario[1],
    Hb_base    = round(base$HB, 2),
    Hb_nadir   = round(post$HB[nadir_i], 2),
    drop_pct   = round(100 * (base$HB - post$HB[nadir_i]) / base$HB, 1),
    nadir_day  = round(post$time[nadir_i], 1),
    Hb_end     = round(tail(post$HB, 1), 2),
    recovered  = (tail(post$HB, 1) - post$HB[nadir_i]) >
                 0.5 * (base$HB - post$HB[nadir_i]),
    astar_min  = round(min(post$ASTAR), 1),
    atrisk_max = round(max(post$ATRISK), 1),
    retic_max  = round(max(post$RETPCT), 1),
    assay_base = round(base$G6PDPCT, 1),
    assay_max  = round(max(post$G6PDPCT), 1),
    methb_max  = round(max(post$METPCT), 1),
    tsb_max    = round(max(post$TSB), 2),
    bfree_max  = round(max(post$BFREE), 1),
    hapto_min  = round(min(post$HPOUT), 2),
    freehb_max = round(max(post$FHBOUT), 3),
    ## mrgsolve emits duplicate time stamps at dose/observation boundaries,
    ## so take the modal positive step rather than the first difference
    hburia_d   = round(sum(post$HBURIA) *
                       median(diff(post$time)[diff(post$time) > 0]), 2),
    crea_max   = round(max(post$CREA), 2),
    row.names  = NULL
  )
}

#' Heterozygous females carry TWO red cell populations, not one.
#' Run the model twice and mix by the normal-cell fraction f.
mix_mosaic <- function(out_normal, out_deficient, f = 0.5) {
  stopifnot(nrow(out_normal) == nrow(out_deficient))
  m <- out_normal
  add <- c("HB", "TSB", "BFREE", "FHBOUT", "METPCT")
  wtd <- c("RETPCT", "G6PDPCT", "HPOUT", "GFRF", "CREA", "EPOOUT")
  for (v in add) m[[v]] <- f * out_normal[[v]] + (1 - f) * out_deficient[[v]]
  for (v in wtd) m[[v]] <- f * out_normal[[v]] + (1 - f) * out_deficient[[v]]
  ## a* is a property of each population, so report the deficient one
  m$ASTAR    <- out_deficient$ASTAR
  m$ATRISK   <- (1 - f) * out_deficient$ATRISK
  m$scenario <- paste0("mosaic f=", f)
  m
}

## =====================================================================
##  4. SCENARIOS
## =====================================================================
##  Every number quoted in the comments below is produced by running this
##  file; nothing is hand-entered. Regenerate with:  source(...); g6pd_all()

g6pd_all <- function() {

  ## -------------------------------------------------------------------
  ## 1  THE DERN EXPERIMENT (1954). A- male, primaquine 30 mg daily,
  ##    HELD FOR 60 DAYS. The haemoglobin must fall, bottom out around
  ##    day 7-9, and then climb back WHILE THE DRUG CONTINUES.
  ## -------------------------------------------------------------------
  ##    RESULT: Hb 15.00 -> 11.77 (-21.5%) at day 5.7, then back to 14.73 by
  ##    day 60 WITH THE DRUG STILL RUNNING. a* = 84 d, so 30% of the mass is
  ##    at risk and 70% is not. Reticulocytes 3.3%. This is Dern's curve.
  s01 <- run_g6pd(patient = g6pd_patient("Aminus"),
                  rx = rx_primaquine(30, days = 60), days = 60,
                  name = "1 A- - primaquine 30 mg daily x60 d")

  ## -------------------------------------------------------------------
  ## 2  THE SAME PRESCRIPTION IN A MEDITERRANEAN PATIENT.
  ##    Same drug, same dose, same "severely deficient" label. The
  ##    age-activity curve is flat, so there is no young survivor pool.
  ## -------------------------------------------------------------------
  ##    RESULT: baseline is ALREADY 14.06 (a* = 99 d at rest, so the oldest
  ##    bin is being trimmed before any drug). Nadir 8.74 (-37.9%) at day 4.1,
  ##    a* falls to 36 d and 70% of the mass is at risk, reticulocytes 8.9%,
  ##    haptoglobin 0.01 g/L and 0.25 d of frank haemoglobinuria. It ends at
  ##    12.55 -- a lower plateau held up by near-maximal erythropoiesis, not a
  ##    recovery. Same drug, same dose, twice the nadir and a different route.
  s02 <- run_g6pd(patient = g6pd_patient("Mediterranean"),
                  rx = rx_primaquine(30, days = 60), days = 60,
                  name = "2 Mediterranean - primaquine 30 mg daily x60 d")

  ## -------------------------------------------------------------------
  ## 3  A NORMAL SUBJECT ON THE SAME REGIMEN -- the negative control that
  ##    shows the drug is not intrinsically haemolytic.
  ## -------------------------------------------------------------------
  ##    RESULT: nothing happens. Hb 15.00 throughout, a* stays at 120 d
  ##    (i.e. no circulating cell of any age is at risk), 0% at risk.
  s03 <- run_g6pd(patient = g6pd_patient("normal"),
                  rx = rx_primaquine(30, days = 60), days = 60,
                  name = "3 normal - primaquine 30 mg daily x60 d")

  ## -------------------------------------------------------------------
  ## 4  TAFENOQUINE 300 mg SINGLE DOSE vs an equivalent primaquine course
  ##    in the same A- patient. The exposure is comparable; the half-life
  ##    is not. Nothing can be withdrawn from arm (a).
  ## -------------------------------------------------------------------
  ##    RESULT 4a: nadir 9.71 (-35.3%) at day 3.2, a* 66 d, 45% at risk,
  ##    reticulocytes 6.1%, free plasma Hb 0.48 g/L and 1.1 days of
  ##    haemoglobinuria. 4b (a THIRTY PERCENT smaller total dose, but
  ##    fast-off): nadir 12.96 (-13.6%), a* 97 d, 19% at risk, no
  ##    haemoglobinuria at all. The difference is not potency, it is that
  ##    tafenoquine is still there on day 15 and primaquine is not.
  s04a <- run_g6pd(patient = g6pd_patient("Aminus"),
                   rx = rx_tafenoquine(300), days = 60,
                   name = "4a A- - tafenoquine 300 mg single dose")
  s04b <- run_g6pd(patient = g6pd_patient("Aminus"),
                   rx = rx_primaquine(15, days = 14), days = 60,
                   name = "4b A- - primaquine 15 mg daily x14 d")

  ## -------------------------------------------------------------------
  ## 5  WHO WEEKLY REGIMEN. Primaquine 45 mg once weekly for 8 weeks in
  ##    the same A- patient -- a HIGHER single dose than scenario 1, and
  ##    a comparable total, but delivered so the marrow can answer it.
  ## -------------------------------------------------------------------
  ##    RESULT: nadir 12.45 (-17.0%) on the FIRST pulse and then a sawtooth
  ##    that never goes deeper -- each later pulse hits an age histogram the
  ##    reticulocytes have already refilled. Compare scenario 1: a bigger
  ##    single dose, a comparable total, and a shallower nadir with no slide.
  ##    An honest caveat the model does volunteer: the first weekly pulse is
  ##    nearly as costly as daily dosing. What weekly buys is the absence of
  ##    accumulation, not a gentler start.
  s05 <- run_g6pd(patient = g6pd_patient("Aminus"),
                  rx = rx_pq_weekly(45, weeks = 8), days = 70,
                  name = "5 A- - primaquine 45 mg WEEKLY x8")

  ## -------------------------------------------------------------------
  ## 6  FAVISM. A 20 kg child, Mediterranean variant, one meal of fresh
  ##    broad beans. Intravascular, hours not days, haptoglobin to zero.
  ## -------------------------------------------------------------------
  ##    RESULT: nadir 4.84 (-65.6%) at day 2.5 -- hours, not days. a* falls to
  ##    2.3 d, so 98% of the mass is at risk and even reticulocytes are not
  ##    safe. Haptoglobin 0.00, free plasma Hb 6.25 g/L, 0.4 d of
  ##    haemoglobinuria, methaemoglobin 12.1%, creatinine 1.19, reticulocytes
  ##    20.2%. This is a transfusion, and the model says so from the a* alone.
  s06 <- run_g6pd(patient = g6pd_patient("Mediterranean", wt = 20),
                  rx = rx_fava(1000), days = 40, delta = 0.02,
                  name = "6 Mediterranean child - fava bean meal")

  ## -------------------------------------------------------------------
  ## 7  DAPSONE, THEN THE ANTIDOTE. Dapsone 100 mg daily raises MetHb;
  ##    methylene blue 2 mg/kg is given on day 30. Run it in a normal
  ##    subject and in an A- subject and compare what the antidote does.
  ## -------------------------------------------------------------------
  rx_dap_mb <- function() c(rx_dapsone(100, days = 45),
                            rx_methyleneblue(140, start = 30))
  ##    RESULT -- and the summary table's MAXIMUM is the wrong statistic here,
  ##    so read the hours after the dose:
  ##      normal : MetHb 7.2% -> 2.2% within 2.4 h. The antidote works.
  ##      A-     : MetHb 7.9% -> 6.6% -> back up to 8.3% by 12 h, AND the
  ##               oxidant flux jumps 0.18 -> 0.68 mmol/L/d, because the dye
  ##               that cannot be reduced to leuco-MB is itself an oxidant.
  ##    Neither number was put in by hand. Both follow from multiplying
  ##    methylene blue's effect by the same NADPH availability that governs
  ##    the haemolysis.
  s07a <- run_g6pd(patient = g6pd_patient("normal"),
                   rx = rx_dap_mb(), days = 45, delta = 0.01,
                   name = "7a normal - dapsone 100 mg + MB on d30")
  s07b <- run_g6pd(patient = g6pd_patient("Aminus"),
                   rx = rx_dap_mb(), days = 45, delta = 0.01,
                   name = "7b A- - dapsone 100 mg + MB on d30")

  ## -------------------------------------------------------------------
  ## 8  RASBURICASE IN TUMOUR LYSIS SYNDROME. Urate oxidase runs
  ##    urate + O2 + H2O -> allantoin + H2O2, ONE MOLE FOR ONE MOLE, so the
  ##    peroxide dose is set by the urate pool being destroyed and not by
  ##    the prescription. Same 0.2 mg/kg across three genotypes and two
  ##    tumour burdens, plus the allopurinol arm that makes no peroxide.
  ## -------------------------------------------------------------------
  ##    RESULT: normal 0.0% | A- -31.2% | Mediterranean -49.6% | the SAME
  ##    Mediterranean patient with a small tumour (-29.9%) | allopurinol
  ##    -0.1%. Two readings. First, the genotype axis: identical prescription,
  ##    three outcomes. Second, the tumour axis: identical prescription and
  ##    identical genotype, and the nadir still moves 20 percentage points
  ##    because the urate pool IS the peroxide dose. Methaemoglobin rises to
  ##    5-6% in the deficient arms and not at all in the normal one, which is
  ##    the documented rasburicase signature.
  tls_hi <- list(KGENUR = 90, URATE_0 = 41.7)   # urate ~20 mg/dL, heavy load
  tls_lo <- list(KGENUR = 20, URATE_0 = 16.7)   # urate ~8 mg/dL
  s08a <- run_g6pd(patient = g6pd_patient("normal"),
                   rx = rx_rasburicase(14, days = 3), days = 30, extra = tls_hi,
                   name = "8a normal - rasburicase, HIGH urate")
  s08b <- run_g6pd(patient = g6pd_patient("Aminus"),
                   rx = rx_rasburicase(14, days = 3), days = 30, extra = tls_hi,
                   name = "8b A- - rasburicase, HIGH urate")
  s08c <- run_g6pd(patient = g6pd_patient("Mediterranean"),
                   rx = rx_rasburicase(14, days = 3), days = 30, extra = tls_hi,
                   name = "8c Mediterranean - rasburicase, HIGH urate")
  s08d <- run_g6pd(patient = g6pd_patient("Mediterranean"),
                   rx = rx_rasburicase(14, days = 3), days = 30, extra = tls_lo,
                   name = "8d Mediterranean - rasburicase, LOW urate")
  s08e <- run_g6pd(patient = g6pd_patient("Mediterranean"),
                   rx = NULL, days = 30, extra = tls_hi,
                   name = "8e Mediterranean - allopurinol (no uricase at all)")

  ## -------------------------------------------------------------------
  ## 9  THE DIAGNOSTIC TRAP. Same A- patient as scenario 1: watch
  ##    G6PDPCT, which is what the laboratory would report, while the
  ##    patient is actively haemolysing.
  ## -------------------------------------------------------------------
  ##    RESULT: the assay reads 16.5% before exposure -- correctly A-. During
  ##    the haemolysis it climbs to 29.5%, because the cells it was counting
  ##    have been deleted and replaced by reticulocytes carrying 1.5x a young
  ##    cell's activity. A patient tested at the nadir looks nearly twice as
  ##    normal as they are, and the number keeps drifting for months. Retest
  ##    at 3 months, one full red cell lifespan.
  s09 <- run_g6pd(patient = g6pd_patient("Aminus"),
                  rx = rx_primaquine(30, days = 10), days = 120,
                  name = "9 A- - 10-day course, assay followed to 120 d")

  ## -------------------------------------------------------------------
  ## 10 NEONATAL JAUNDICE: G6PD x UGT1A1 is a PRODUCT. Four newborns,
  ##    3 kg, followed from 12 h of life through the bilirubin peak.
  ## -------------------------------------------------------------------
  neo <- function(var, ugt, nm)
    run_g6pd(patient = g6pd_patient(var, wt = 3, ugt = ugt,
                                    neonate = TRUE, pna = 0.5),
             rx = NULL, days = 12, lead_in = 0, delta = 0.02, name = nm)
  ##    RESULT, peak total serum bilirubin (mg/dL) and peak FREE bilirubin (nM,
  ##    the species that actually crosses the blood-brain barrier):
  ##        neither hit    TSB  4.62   Bf    3.0
  ##        G6PD only      TSB 15.57   Bf   20.7
  ##        UGT1A1 only    TSB  8.58   Bf    6.9
  ##        BOTH           TSB 23.15   Bf  104.4
  ##    On the TSB scale the interaction looks nearly additive (19.5 predicted
  ##    vs 23.15 observed). On the FREE bilirubin scale it is 4.2x the
  ##    additive prediction (24.6 predicted vs 104.4 observed) -- because
  ##    albumin binding saturates, and the toxic species is the unbound one.
  ##    THE NON-LINEARITY LIVES EXACTLY WHERE THE DAMAGE DOES, and a model
  ##    that only tracked TSB would have missed it.
  s10a <- neo("normal",        "6/6", "10a neonate - neither hit")
  s10b <- neo("Mediterranean", "6/6", "10b neonate - G6PD only")
  s10c <- neo("normal",        "7/7", "10c neonate - UGT1A1 7/7 only")
  s10d <- neo("Mediterranean", "7/7", "10d neonate - BOTH hits")

  ## -------------------------------------------------------------------
  ## 11 HETEROZYGOUS FEMALE. Two populations in one circulation. The
  ##    whole-blood assay averages them into the normal range while a
  ##    real fraction of her red cells is at risk.
  ## -------------------------------------------------------------------
  hn <- run_g6pd(patient = g6pd_patient("normal"),
                 rx = rx_tafenoquine(300), days = 60, name = "het-normal-pop")
  hd <- run_g6pd(patient = g6pd_patient("Mediterranean"),
                 rx = rx_tafenoquine(300), days = 60, name = "het-def-pop")
  ##    RESULT: her whole-blood assay reads 60.9% -- under the FDA's 70%
  ##    tafenoquine threshold but comfortably inside what most laboratories
  ##    flag as normal, and far above the 30% line. Her haemoglobin still
  ##    falls 22.9%, because 40% of her red cells have a* = 19 d.
  s11 <- mix_mosaic(hn, hd, f = 0.6)
  s11$scenario <- "11 heterozygous female (60% normal) - tafenoquine 300 mg"

  ## -------------------------------------------------------------------
  ## 12 REMOVE THE RESCUE. Scenario 1 repeated with the marrow knocked
  ##    out (parvovirus B19 aplastic crisis, FERY = 0.15). Self-limitation
  ##    is not a property of the drug or the variant -- it is a property
  ##    of the marrow being able to answer.
  ## -------------------------------------------------------------------
  ##    RESULT: baseline is 10.65 (the marrow alone cannot hold 15). Nadir
  ##    8.24 (-22.7%) and it ends at 10.21 instead of 14.73. The DEPTH of the
  ##    nadir is unchanged -- a* does not know about the marrow -- but the
  ##    recovery is gone. Self-limitation was never a property of the drug or
  ##    of the variant; it was the reticulocyte arm answering.
  s12 <- run_g6pd(patient = g6pd_patient("Aminus", marrow = 0.15),
                  rx = rx_primaquine(30, days = 60), days = 60,
                  name = "12 A- + aplastic crisis - primaquine 30 mg")

  ## -------------------------------------------------------------------
  ## 13 CYP2D6 POOR METABOLISER. The enzyme that bioactivates primaquine to
  ##    the oxidant is the same one that bioactivates it to the cure. A PM
  ##    on 30 mg daily haemolyses not at all -- and is not cured either.
  ##    Note the assay still reads 16.5%: the genotype has not changed.
  ## -------------------------------------------------------------------
  ##    RESULT: Hb 15.00 flat, a* 120 d, 0% at risk -- and the assay still
  ##    reads 16.5%. The genotype is unchanged; only the bioactivation is
  ##    gone. The same CYP2D6 that makes the oxidant makes the cure, so this
  ##    patient is spared the haemolysis and also fails radical cure.
  s13 <- run_g6pd(patient = g6pd_patient("Aminus", cyp2d6 = 0),
                  rx = rx_primaquine(30, days = 60), days = 60,
                  name = "13 A- CYP2D6 poor metaboliser - primaquine 30 mg")

  runs <- list(s01, s02, s03, s04a, s04b, s05, s06, s07a, s07b,
               s08a, s08b, s08c, s08d, s08e, s09, s10a, s10b, s10c, s10d,
               s11, s12, s13)
  summ <- do.call(rbind, lapply(runs, summarise_g6pd))
  list(runs = runs, summary = summ)
}

## =====================================================================
##  5. THE TWO TABLES THAT MAKE THE ARGUMENT
## =====================================================================

#' Table 1. The genotype IS the age-activity curve.
#' Print E(a) for each variant at a range of cell ages, plus the critical
#' age a* under three oxidant loads. This table needs no simulation at
#' all -- it is arithmetic on E0 and TAU -- which is rather the point.
g6pd_age_table <- function(ages = c(0, 20, 40, 60, 80, 100, 120),
                           ox = c(baseline = 0.006, primaquine30 = 0.32,
                                  fava = 60)) {
  vmax <- 60
  do.call(rbind, lapply(names(VARIANTS), function(v) {
    p <- VARIANTS[[v]]
    e <- p$E0 * exp(-log(2) * ages / p$TAU)
    astar <- sapply(ox, function(o) {
      if (p$E0 * vmax <= o) return(0)
      min(120, (p$TAU / log(2)) * log(p$E0 * vmax / o))
    })
    df <- as.data.frame(c(
      list(variant = p$label, E0 = p$E0, TAU = p$TAU),
      setNames(as.list(signif(e, 3)), paste0("E_", ages, "d")),
      setNames(as.list(round(astar, 1)), paste0("astar_", names(ox)))
    ))
    df
  }))
}

#' Table 2. The headline comparison: everything in scenarios 1-3.
g6pd_headline <- function(res) {
  res$summary[1:3, c("scenario", "Hb_base", "Hb_nadir", "drop_pct",
                     "nadir_day", "Hb_end", "recovered", "astar_min",
                     "atrisk_max", "retic_max")]
}

if (identical(environment(), globalenv()) &&
    !is.null(getOption("g6pd.run")) && isTRUE(getOption("g6pd.run"))) {
  res <- g6pd_all()
  print(g6pd_age_table())
  print(res$summary)
}
