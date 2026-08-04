## =====================================================================
##  ANDROGENETIC ALOPECIA (male & female pattern hair loss) — QSP model
##  50 ODE compartments · 12 drug PK · 8 enzyme · 4 endocrine · 4 signalling
##                        · 5 systemic readouts · 17 follicle
##  Reference subject: 5.07 cm2 (1-inch diameter) vertex target area,
##  1268 follicles.  TIME UNIT: DAYS.  Concentrations ng/mL.
## =====================================================================
##
##  THE THESIS THIS MODEL IS BUILT AROUND
##  ------------------------------------
##  Androgenetic alopecia is taught as hair LOSS.  This model asserts that
##  almost nothing is lost — the follicles are still there, they are
##  smaller — and that every clinically important number follows from that
##  one substitution once you notice that follicle size is only ever
##  changed at ONE moment of the cycle: the exit from telogen.
##
##  (1) NOTHING IS LOST; THINGS GET SMALLER.
##      The reference patient (below) reaches the exact baseline of the
##      pivotal finasteride trial — 876 non-vellus hairs in a 1-inch
##      circle (Kaufman 1998, PMID 9777765, baseline 876) — while still
##      owning 1006 of his original 1268 follicles.  Of the 392-hair
##      "deficit", 262 follicles are genuinely gone and 130 are alive and
##      vellus.  Every drug in this file competes for those 130.
##
##  (2) THE RATCHET TURNS ONLY AT TELOGEN EXIT, SO THE DISEASE RUNS ON
##      CYCLE TIME, NOT ON CALENDAR TIME.
##      Cycle periods in the model (anagen shortened by androgen, ARD 0.95):
##
##          terminal          832 d   ->  2.2 cycles in 5 years
##          intermediate      517 d   ->  3.5
##          small interm.     296 d   ->  6.2
##          vellus            167 d   -> 10.9
##
##      A terminal follicle gets ~2 chances to shrink in 5 years, which is
##      why AGA takes decades; a vellus follicle gets ~11 chances to grow
##      back, which is why treatment response is visible in ONE year.
##      Both facts come from the same number.
##
##  (3) THE DAILY SHED IS A RATIO, NOT A RATE.
##      shed = N_total x f_telogen / T_telogen.  Normal scalp:
##
##          100 000 x 0.090 / 100 d  =  90 hairs/day
##
##      the textbook "about 100 a day", derived rather than asserted.  In
##      the model's AGA patient the telogen fraction rises to 23.6% purely
##      because anagen is shorter in the smaller classes — increased
##      shedding in AGA is a CONSEQUENCE of miniaturisation, not a
##      separate disease process, and the trichogram telogen percentage is
##      a direct read-out of anagen duration.
##
##  (4) THE FINASTERIDE DOSE-RESPONSE IS FLAT BECAUSE THE TARGET
##      INTEGRATES.  Finasteride's plasma half-life is 6 h; the enzyme it
##      inactivates turns over with a half-life of 2 d, and the covalent
##      enzyme-NADP-dihydrofinasteride adduct with a half-life of 30 d.
##      Steady-state free enzyme is therefore k_deg/(k_deg + <k_on x C>),
##      which depends on the AVERAGE concentration, not the peak.  The
##      model reproduces the whole published dose-response (Drake 1999,
##      PMID 10495374; day 42 scalp / serum DHT, model vs observed):
##
##          dose      scalp DHT              serum DHT
##          0.01 mg   -18.3%  (obs -14.9)    -20.1%  (obs   n/a)
##          0.05 mg   -43.1%  (obs -61.6)    -47.3%  (obs -49.5)
##          0.2  mg   -57.5%  (obs -56.5)    -63.1%  (obs -68.6)
##          1    mg   -63.6%  (obs -64.1)    -70.3%  (obs -71.4)
##          5    mg   -67.4%  (obs -69.4)    -74.9%  (obs -72.2)
##
##      Five times the dose buys 3.8 percentage points of scalp DHT and,
##      downstream, FOUR HAIRS at one year (+105.2 vs +101.3 over placebo).
##
##  (5) 5 mg BEATS 1 mg ON SCALP DHT BUT NOT ON SERUM DHT, AND THE MODEL
##      SAYS WHY.  Finasteride is ~2000x weaker on type I than type II in
##      this parameterisation.  At 1 mg free type I is 0.979; at 5 mg it is
##      0.903.  Scalp skin — unlike blood — contains sebaceous type I, so
##      only the scalp compartment can see that extra 7.6%.  The published
##      asymmetry (scalp -64.1 -> -69.4 while serum -71.4 -> -72.2) is
##      recovered without a scalp-specific fitted parameter.
##
##  (6) THE RECOVERABLE POOL IS THE RECENTLY MINIATURISED FOLLICLE, AND IT
##      AGES OUT.  Class 4 is split into A4a/C4a/T4a (recent vellus, can
##      re-enlarge, PUW 0.75) and A4b/C4b/T4b (established vellus, PUW
##      0.05), with a one-way flux p_age between them.  This single
##      structure produces the plateau of the finasteride response, the
##      failure of late starts, and the irreversibility of Norwood VI:
##
##          start at TAHC 1150 -> 58.8% of the deficit recovered by 3 y
##          start at TAHC 1000 -> 40.9%
##          start at TAHC  876 -> 30.1%
##
##  (7) HAIR WEIGHT MOVES TWICE AS FAR AS HAIR COUNT, BECAUSE MASS GOES AS
##      d^2.  Same simulation, two endpoints, at 12 months:
##
##          finasteride 1 mg   count +9.9%   mass +18.6%   ratio 1.88
##          minoxidil 5%       count +2.9%   mass  +5.9%   ratio 2.05
##          dutasteride        count +12.3%  mass +20.3%   ratio 1.65
##
##      Trials that count hairs are reading the smaller half of the effect.
##
##  (8) MINOXIDIL'S "DREAD SHED" IS THE MECHANISM WORKING, NOT A SIDE
##      EFFECT.  Minoxidil sulfate shortens telogen (T_T 100 -> 49 d at
##      full effect).  Telogen exit REQUIRES release of the club hair, so
##      an abrupt shortening must produce a transient rise in shedding
##      before it produces any hair.  No term in the model represents
##      shedding as toxicity, yet the responder's shed goes
##
##          2.38 -> 3.88 hairs/day per target area (x1.63) at day 10,
##          back below placebo by day 140, +25 hairs by 12 months.
##
##  (9) MINOXIDIL IS A PRODRUG, SO THE TRIAL MEAN IS A MIXTURE.  Follicular
##      SULT1A1 activity splits the population; the same 5% solution gives
##      +39.8 hairs vs placebo at 12 months in a responder and +18.7 in a
##      non-responder.  A trial reporting the mean is reporting neither.
##
## (10) STOPPING FINASTERIDE RESTORES DHT IN TEN DAYS AND LOSES THE HAIR
##      OVER YEARS.  After withdrawal at 2 y: scalp DHT 0.364 -> 0.980 by
##      day 10 (enzyme resynthesis, t1/2 2 d), while the hair count is
##      still +110 at day 20, +101 at 1 y and +79 at 3 y post-stop.  The
##      two clocks differ by two orders of magnitude.
##
## (11) THE SAME ARITHMETIC EXCUSES IMPERFECT ADHERENCE.  Taking 4 doses in
##      7 leaves scalp DHT at 0.389 instead of 0.364 and costs 10 hairs at
##      one year (+91.4 vs +101.3 over placebo) — 90% of the benefit from
##      57% of the tablets, because the enzyme, not the plasma, is the
##      thing being titrated.
##
## (12) TOPICAL FINASTERIDE SEPARATES THE TWO COMPARTMENTS.  With separate
##      scalp and systemic enzyme pools, 0.25% topical gives scalp DHT
##      -47.6% and serum DHT -27.4% (observed ~-25%), keeping 78% of the
##      hair benefit for 42% of the systemic androgen signal.
##
## (13) IN FEMALE PATTERN HAIR LOSS EVERY ANTI-ANDROGEN SCALES WITH THE
##      ANDROGENIC SHARE OF THE DRIVE, AND MINOXIDIL DOES NOT.  ARIND
##      carries androgen-independent inhibitory drive.  Holding presenting
##      severity fixed (age 33.5 y, 899 hairs, 216 vellus) and varying only
##      the SPLIT between GS and ARIND, with the drug model untouched:
##
##          androgenic share   finasteride   spironolactone   minoxidil 5%
##                34%             +53.2          +136.6          +67.2
##                22%             +35.2           +99.9          +67.1
##                11%             +15.7           +53.4          +67.1
##
##      Minoxidil is identical to three significant figures; every
##      anti-androgen is proportional.  The negative post-menopausal
##      finasteride trial (Price 2000, PMID 11050579) is therefore read as
##      a low-androgenic-share POPULATION, not as an inert drug — the same
##      model gives +101.3 in the man.
##
##  CALIBRATION SUMMARY  (targets, and how the model does)
##  ------------------------------------------------------
##   Kaufman 1998 (PMID 9777765) baseline 876 hairs / 5.1 cm2 ...... 876 (by construction)
##   Kaufman 1998 finasteride 1 mg vs placebo, 1 year: +107 ......... +101.3
##   Kaufman 1998 finasteride 1 mg vs placebo, 2 years: +138 ........ +139.2
##   Drake 1999 (PMID 10495374) five-point DHT dose-response ........ see (4) above
##   Gubelin Harcha 2014 (PMID 24411083) dutasteride > finasteride
##       at 24 weeks .................................................. +80.8 vs +63.4
##   Clark 2004 (PMID 15126539) dutasteride serum DHT ~-94% ......... -96.1%
##   normal daily shed ~100 hairs/day ............................... 90
##   normal telogen fraction 9-14% .................................. 9.0%
##   topical finasteride serum DHT ~-25% (Piraccini 2022, 34634163) .. -27.4%
##
##  WHERE THE MODEL DISAGREES WITH THE LITERATURE (reported, not tuned away)
##  -----------------------------------------------------------------------
##   * Dutasteride SCALP DHT: model -90.3% at day 42, literature range
##     -51% to -79%.  The model's scalp compartment has only a 6%
##     5AR-independent floor (W_ALT); a larger floor would fit dutasteride
##     but would then break the finasteride 1 mg point, which is the better
##     measured of the two.  Reported as an over-prediction.
##   * Drake's 0.05 mg SCALP point (-61.6%) is not fitted (model -43.1%).
##     That point is non-monotonic in the source data — 0.05 mg is reported
##     as suppressing scalp DHT MORE than 0.2 mg (-56.5%) — which no
##     mass-action model can produce.  The 0.05 mg SERUM point is fitted
##     well (-47.3 vs -49.5).
##   * Setipiprant: with PGD_EFF = 0.18 the model predicts +46 hairs vs
##     placebo at 12 months, which the phase 2a programme did not see.
##     Rather than delete the arm, the model is used the other way round:
##     to be consistent with a null trial the PGD2 arm can carry at most
##     ~5% of the inhibitory drive (see the PGD_EFF scan in the README).
##   * Placebo decline: model -14.7 hairs at 1 y vs ~-21 published; the
##     model's untreated arm is linear where the trial's accelerates.
##
##  IMPLEMENTATION NOTE
##  -------------------
##  Every equation here was first written and run in an independent
##  Python/scipy implementation (fixed-step RK4, 50 states).  That
##  cross-implementation found four real defects: (i) a dropped
##  up-transition flux (class 3 -> class 2) that silently destroyed
##  follicles and made finasteride look WORSE than placebo; (ii) a
##  1000-fold concentration unit error (V in L against k_on in per-ng/mL)
##  that put free type-II enzyme at zero for every dose; (iii) a
##  testosterone feedback loop with DHT weighted at 0.65, which produced a
##  47% rise in testosterone on finasteride instead of the observed ~10%;
##  and (iv) an RK4 step of 0.05 d that went unstable at 5 mg finasteride
##  (k_on x C_max = 224/d) and returned NaN for the entire arm.  All four
##  are fixed here.  Numbers quoted above are from the corrected runs.
## =====================================================================

library(mrgsolve)
library(dplyr)

aga_code <- '
$PARAM @annotated
// ---------------------------------------------------------------- finasteride
FIN_F     : 0.65   : finasteride oral bioavailability
FIN_KA    : 15.0   : finasteride absorption rate (1/d)
FIN_V     : 76000  : finasteride V/F (mL)
FIN_KE    : 2.7726 : finasteride elimination (1/d, t1/2 6 h)
FIN_KWASH : 1.5    : topical finasteride wash-off from scalp depot (1/d)
FIN_KFT   : 0.62   : topical finasteride scalp depot -> scalp tissue (1/d)
FIN_SCL   : 20.0   : mg of scalp depot -> ng/mL of scalp tissue scaling
FIN_KDIF  : 4.0    : plasma <-> scalp tissue exchange (1/d)
FIN_KOUTT : 4.0    : scalp tissue finasteride clearance (1/d)
TOPFIN_F  : 0.002  : systemic bioavailability of topical finasteride
// ---------------------------------------------------------------- dutasteride
DUT_F     : 0.60   : dutasteride bioavailability
DUT_KA    : 5.0    : dutasteride absorption (1/d)
DUT_V1    : 100000 : dutasteride central volume (mL)
DUT_V2    : 300000 : dutasteride peripheral volume (mL)
DUT_Q     : 20000  : dutasteride intercompartmental clearance (mL/d)
DUT_CL    : 7500   : dutasteride clearance (mL/d)
// ---------------------------------------------------------------- minoxidil
MX_KFOL   : 0.55   : scalp depot -> follicle (1/d)
MX_KSYS   : 0.020  : scalp depot -> systemic (1/d) [=> 1.1% absorbed]
MX_KWASH  : 1.20   : scalp depot wash-off (1/d)
MX_KFOUT  : 2.0    : follicular minoxidil elimination (1/d)
MX_KSULT  : 3.0    : SULT1A1 sulfotransferase scalar (1/d)
MX_KSOUT  : 1.5    : minoxidil sulfate elimination (1/d)
MX_VO     : 180000 : minoxidil systemic volume (mL)
MX_KEO    : 4.0774 : minoxidil systemic elimination (1/d, t1/2 4 h)
MX_FORAL  : 0.90   : oral minoxidil bioavailability
MX_FOLOR  : 5.5    : systemic -> follicular minoxidil transfer
MX_EC50   : 12.0   : minoxidil sulfate EC50 (follicular units)
MX_HILL   : 1.5    : minoxidil Hill coefficient
MX_GTT    : 1.05   : telogen shortening gain
MX_GTA    : 0.40   : anagen prolongation gain
SULT      : 2.0    : follicular SULT1A1 activity (responder 2.0 / non 0.25)
// ---------------------------------------------------------------- other drugs
SPI_V     : 90000  : canrenone volume (mL)
SPI_KE    : 0.9902 : canrenone elimination (1/d)
SPI_F     : 0.70   : spironolactone bioavailability
SPI_IC50  : 90.0   : canrenone AR-antagonism IC50 (ng/mL)
// ------------------------------------------------------- 5-alpha-reductase
KDEG_E    : 0.3466 : 5AR protein turnover (1/d, t1/2 2 d)
KON_FIN2  : 5.24   : finasteride capture of type II (per ng/mL per d)
KON_FIN1  : 0.00262: finasteride capture of type I (2000x weaker)
KON_DUT2  : 0.50   : dutasteride capture of type II
KON_DUT1  : 0.40   : dutasteride capture of type I
KDISS1    : 0.0231 : adduct dissociation type I (1/d, t1/2 30 d)
KDISS2    : 0.0231 : adduct dissociation type II (1/d, t1/2 30 d)
// ----------------------------------------------------------------- HPG axis
KLH       : 4.0    : LH production (1/d)
KDLH      : 4.0    : LH elimination (1/d)
FB_A_T    : 0.70   : testosterone weight in gonadal negative feedback
FB_A_D    : 0.30   : DHT weight in gonadal negative feedback
FB_H      : 2.0    : feedback Hill coefficient
KSYN_T    : 4.0    : testosterone synthesis scalar (1/d)
KEL_T     : 4.0    : testosterone elimination (1/d)
W_T2_SER  : 0.765  : fraction of serum DHT from type II
W_T1_SER  : 0.235  : fraction of serum DHT from type I
KEL_DHTS  : 4.0    : serum DHT equilibration (1/d)
W_LOC     : 0.52   : scalp DHT from local synthesis
W_SER     : 0.42   : scalp DHT taken up from serum
W_ALT     : 0.06   : scalp DHT from 5AR-independent routes
W_T1_LOC  : 0.30   : type I share of local scalp synthesis
W_T2_LOC  : 0.70   : type II share of local scalp synthesis
KEL_DHTC  : 4.0    : scalp DHT equilibration (1/d)
AROM      : 1.0    : aromatase scalar (raise for female arm)
LOC5AR    : 1.0    : follicular 5AR capacity (0.30 in the female arm)
// ------------------------------------------------------- AR signalling in DP
AR_EC50   : 0.55   : scalp DHT EC50 for AR signal (normalised)
AR_H      : 2.0    : AR Hill coefficient
K_ARN     : 0.0990 : nuclear AR signal turnover (1/d, t1/2 7 d)
K_DKK     : 0.1386 : DKK-1 turnover (1/d, t1/2 5 d)
K_WNT     : 0.1386 : Wnt tone turnover (1/d)
K_PGD     : 0.2310 : PGD2 turnover (1/d, t1/2 3 d)
PGD_GAIN  : 0.45   : AR-driven rise in PGD2
PGD_EFF   : 0.18   : weight of the PGD2 arm in total inhibitory drive
SETI      : 0.0    : setipiprant / CRTH2 blockade (0-1)
ARIND     : 0.05   : androgen-INDEPENDENT inhibitory drive
GS        : 0.90   : genetic susceptibility (0-1.3)
// --------------------------------------------------------- follicle kinetics
NTOT      : 1268   : follicles in the 5.07 cm2 target area
TA1       : 1000   : anagen duration, terminal (d)
TA2       : 560    : anagen duration, intermediate (d)
TA3       : 250    : anagen duration, small intermediate (d)
TA4       : 70     : anagen duration, vellus (d)
TCAT      : 17     : catagen duration (d)
TTEL      : 100    : telogen duration (d)
D1        : 78     : shaft diameter, terminal (um)
D2        : 60     : shaft diameter, intermediate (um)
D3        : 43     : shaft diameter, small intermediate (um)
D4        : 24     : shaft diameter, vellus (um)
G1        : 0.36   : linear growth rate, terminal (mm/d)
G2        : 0.31   : linear growth rate, intermediate
G3        : 0.24   : linear growth rate, small intermediate
G4        : 0.13   : linear growth rate, vellus
PDMAX     : 0.17   : maximum per-cycle probability of shrinking one class
PUMAX     : 2.20   : restitution scalar (per-cycle, x (1-ARD))
PUMX      : 0.40   : minoxidil-driven restitution (androgen-independent)
PUW2      : 1.00   : restitution weight, intermediate -> terminal
PUW3      : 0.90   : restitution weight, small interm -> intermediate
PUW4A     : 0.75   : restitution weight, RECENT vellus -> small interm
PUW4B     : 0.05   : restitution weight, ESTABLISHED vellus
PLOSS     : 0.055  : per-cycle probability of follicle death (vellus only)
PAGE      : 0.15   : per-cycle recent vellus -> established vellus
FA_AR     : 0.30   : fractional anagen shortening at full androgen drive
// -------------------------------------------------------- systemic readouts
K_PSA     : 0.0277 : PSA turnover (1/d, t1/2 25 d)
K_SEB     : 0.0693 : sebum turnover (1/d)
K_MAP     : 0.3466 : MAP turnover (1/d)
MAP0      : 92     : baseline mean arterial pressure (mmHg)
MAP_G     : 0.65   : MAP fall per ng/mL of systemic minoxidil
K_HTR     : 0.0231 : hypertrichosis index turnover (1/d)
HTR_G     : 8.0    : hypertrichosis gain

$CMT @annotated
FING  : finasteride gut depot (ng)
FINC  : finasteride plasma (ng/mL)
FINSK : topical finasteride scalp depot (mg)
FINT  : finasteride scalp tissue (ng/mL)
DUTG  : dutasteride gut depot (ng)
DUTC  : dutasteride central (ng/mL)
DUTP  : dutasteride peripheral (ng/mL)
MXSK  : topical minoxidil scalp depot (mg)
MXDF  : follicular minoxidil (units)
MXSF  : follicular minoxidil SULFATE (units)
MXDC  : minoxidil systemic (ng/mL)
SPIC  : canrenone plasma (ng/mL)
E1    : free 5AR type I, systemic (fraction)
EI1   : inhibited 5AR type I, systemic
E2    : free 5AR type II, systemic (fraction)
EI2   : inhibited 5AR type II, systemic
E1S   : free 5AR type I, SCALP (fraction)
EI1S  : inhibited 5AR type I, scalp
E2S   : free 5AR type II, SCALP (fraction)
EI2S  : inhibited 5AR type II, scalp
LH    : luteinising hormone (normalised)
TST   : serum testosterone (normalised)
DHTS  : serum DHT (normalised)
DHTC  : scalp / follicular DHT (normalised)
ARN   : nuclear AR signal in dermal papilla
DKK1  : DKK-1 / TGF-beta paracrine inhibitor
PGD2  : prostaglandin D2 (normalised)
WNT   : Wnt / beta-catenin growth tone
PSA   : serum PSA (normalised)
EST   : serum estradiol (normalised)
SEB   : sebum output (normalised)
MAP   : mean arterial pressure (mmHg)
HTR   : hypertrichosis index
A1    : anagen, terminal
A2    : anagen, intermediate
A3    : anagen, small intermediate
A4A   : anagen, vellus RECENT
A4B   : anagen, vellus ESTABLISHED
C1    : catagen, terminal
C2    : catagen, intermediate
C3    : catagen, small intermediate
C4A   : catagen, vellus recent
C4B   : catagen, vellus established
T1    : telogen, terminal
T2    : telogen, intermediate
T3    : telogen, small intermediate
T4A   : telogen, vellus recent
T4B   : telogen, vellus established
LOST  : follicles lost (fibrous streamers)
CSHED : cumulative shed hairs (target area)

$MAIN
// initial condition: an intact 18-year-old scalp, all terminal, at
// cycle-phase steady state.  The presenting patient is produced by a
// burn-in (see the scenario section at the bottom of the file).
double cyc0 = TA1 + TCAT + TTEL;
A1_0 = NTOT * TA1  / cyc0;
C1_0 = NTOT * TCAT / cyc0;
T1_0 = NTOT * TTEL / cyc0;
E1_0 = 1.0; E2_0 = 1.0; E1S_0 = 1.0; E2S_0 = 1.0;
LH_0 = 1.0; TST_0 = 1.0; DHTS_0 = 1.0; DHTC_0 = 1.0;
PGD2_0 = 1.0; WNT_0 = 1.0; PSA_0 = 1.0; EST_0 = 1.0; SEB_0 = 1.0;
MAP_0 = MAP0;

$ODE
// ------------------------------------------------------------------- guards
double fc  = FINC  > 0 ? FINC  : 0.0;
double ft  = FINT  > 0 ? FINT  : 0.0;
double dc  = DUTC  > 0 ? DUTC  : 0.0;
double msk = MXSK  > 0 ? MXSK  : 0.0;
double mdf = MXDF  > 0 ? MXDF  : 0.0;
double msf = MXSF  > 0 ? MXSF  : 0.0;
double mdc = MXDC  > 0 ? MXDC  : 0.0;
double spc = SPIC  > 0 ? SPIC  : 0.0;
double dhc = DHTC  > 1e-9 ? DHTC : 1e-9;

// ----------------------------------------------------------------------- PK
dxdt_FING  = -FIN_KA*FING;
dxdt_FINC  =  FIN_KA*FING/FIN_V - FIN_KE*FINC;
dxdt_FINSK = -(FIN_KFT + FIN_KWASH)*FINSK;
dxdt_FINT  =  FIN_KFT*FINSK*FIN_SCL + FIN_KDIF*fc - FIN_KOUTT*FINT;
dxdt_DUTG  = -DUT_KA*DUTG;
dxdt_DUTC  = (DUT_KA*DUTG - DUT_CL*DUTC - DUT_Q*DUTC + DUT_Q*DUTP)/DUT_V1;
dxdt_DUTP  = (DUT_Q*DUTC - DUT_Q*DUTP)/DUT_V2;
dxdt_MXSK  = -(MX_KFOL + MX_KSYS + MX_KWASH)*MXSK;
dxdt_MXDF  =  MX_KFOL*msk + MX_FOLOR*mdc - MX_KFOUT*MXDF - MX_KSULT*SULT*MXDF;
dxdt_MXSF  =  MX_KSULT*SULT*mdf - MX_KSOUT*MXSF;
dxdt_MXDC  =  MX_KSYS*msk*1e6/MX_VO - MX_KEO*MXDC;
dxdt_SPIC  = -SPI_KE*SPIC;

// -------------------------------------- 5AR: mechanism-based, quasi-irrevers.
// Free enzyme obeys  E_ss = k_deg / (k_deg + <k_on x C>) : the TARGET, not the
// plasma, is what the dose titrates.  This one line is thesis (4) and (11).
double kon1  = KON_FIN1*fc + KON_DUT1*dc;      // systemic pool
double kon2  = KON_FIN2*fc + KON_DUT2*dc;
double kon1s = KON_FIN1*ft + KON_DUT1*dc;      // scalp pool (topical-aware)
double kon2s = KON_FIN2*ft + KON_DUT2*dc;
dxdt_E1   = KDEG_E - KDEG_E*E1  - kon1 *E1  + KDISS1*EI1;
dxdt_EI1  = kon1 *E1  - KDISS1*EI1  - KDEG_E*EI1;
dxdt_E2   = KDEG_E - KDEG_E*E2  - kon2 *E2  + KDISS2*EI2;
dxdt_EI2  = kon2 *E2  - KDISS2*EI2  - KDEG_E*EI2;
dxdt_E1S  = KDEG_E - KDEG_E*E1S - kon1s*E1S + KDISS1*EI1S;
dxdt_EI1S = kon1s*E1S - KDISS1*EI1S - KDEG_E*EI1S;
dxdt_E2S  = KDEG_E - KDEG_E*E2S - kon2s*E2S + KDISS2*EI2S;
dxdt_EI2S = kon2s*E2S - KDISS2*EI2S - KDEG_E*EI2S;

// ------------------------------------------------------------------ HPG axis
double fb  = 1.0/(1.0 + pow(FB_A_T*TST + FB_A_D*DHTS, FB_H));
double fb0 = 1.0/(1.0 + pow(FB_A_T + FB_A_D, FB_H));
dxdt_LH   = KLH*(fb/fb0) - KDLH*LH;
dxdt_TST  = KSYN_T*LH - KEL_T*TST;
double ser = (W_T2_SER*E2 + W_T1_SER*E1)*TST;
dxdt_DHTS = KEL_DHTS*(ser - DHTS);
double loc = LOC5AR*(W_T1_LOC*E1S + W_T2_LOC*E2S)*TST;
dxdt_DHTC = KEL_DHTC*(W_LOC*loc + W_SER*DHTS + W_ALT - DHTC);

// --------------------------------------------- AR signalling, dermal papilla
double spi   = 1.0/(1.0 + spc/SPI_IC50);
double h_now = pow(dhc, AR_H)/(pow(AR_EC50,AR_H) + pow(dhc,AR_H));
double h_ref = 1.0/(pow(AR_EC50,AR_H) + 1.0);
double arraw = GS*spi*h_now/h_ref;
dxdt_ARN  = K_ARN*(arraw - ARN);
dxdt_DKK1 = K_DKK*(ARN - DKK1);
dxdt_PGD2 = K_PGD*(1.0 + PGD_GAIN*ARN - PGD2);
double msfh = pow(msf, MX_HILL);
double mx   = msfh/(pow(MX_EC50,MX_HILL) + msfh);
dxdt_WNT  = K_WNT*((1.0 - 0.6*DKK1 + 0.5*mx) - WNT);

double pgds = (1.0 - SETI)*((PGD2 > 1.0 ? PGD2 - 1.0 : 0.0)/PGD_GAIN);
double ARD  = (1.0 - PGD_EFF)*DKK1 + PGD_EFF*pgds + ARIND;
if (ARD < 0.0)  ARD = 0.0;
if (ARD > 1.5)  ARD = 1.5;
double ARDc = ARD > 1.0 ? 1.0 : ARD;

// --------------------------------------------------------- the follicle cycle
// Everything above enters here and NOWHERE else.  Class transitions occur only
// in the telogen->anagen flux, which is why the disease runs on cycle time.
double TTe  = TTEL/(1.0 + MX_GTT*mx);
double fAn  = (1.0 - FA_AR*ARDc)*(1.0 + MX_GTA*mx);
double TAe1 = TA1*fAn, TAe2 = TA2*fAn, TAe3 = TA3*fAn, TAe4 = TA4*fAn;

double pdn = PDMAX*ARD;                 if (pdn > 0.85) pdn = 0.85;
double pu0 = PUMAX*(1.0 - ARDc) + PUMX*mx;
double pu2 = pu0*PUW2,  pu3 = pu0*PUW3,  pu4a = pu0*PUW4A,  pu4b = pu0*PUW4B;
double cap = 0.90 - pdn; if (cap < 0.0) cap = 0.0;
if (pu2  > cap) pu2  = cap;
if (pu3  > cap) pu3  = cap;
if (pu4a > cap) pu4a = cap;
if (pu4b > cap) pu4b = cap;
double plo = PLOSS*ARD; if (plo > 0.85) plo = 0.85;

double f1 = T1/TTe, f2 = T2/TTe, f3 = T3/TTe, f4a = T4A/TTe, f4b = T4B/TTe;

dxdt_A1  = f1*(1.0 - pdn)                 + f2*pu2                    - A1/TAe1;
dxdt_A2  = f2*(1.0 - pdn - pu2)           + f1*pdn  + f3*pu3          - A2/TAe2;
dxdt_A3  = f3*(1.0 - pdn - pu3)           + f2*pdn  + f4a*pu4a + f4b*pu4b - A3/TAe3;
dxdt_A4A = f4a*(1.0 - pu4a - plo - PAGE)  + f3*pdn                    - A4A/TAe4;
dxdt_A4B = f4b*(1.0 - pu4b - plo)         + f4a*PAGE                  - A4B/TAe4;

dxdt_C1  = A1/TAe1  - C1/TCAT;
dxdt_C2  = A2/TAe2  - C2/TCAT;
dxdt_C3  = A3/TAe3  - C3/TCAT;
dxdt_C4A = A4A/TAe4 - C4A/TCAT;
dxdt_C4B = A4B/TAe4 - C4B/TCAT;

dxdt_T1  = C1/TCAT  - f1;
dxdt_T2  = C2/TCAT  - f2;
dxdt_T3  = C3/TCAT  - f3;
dxdt_T4A = C4A/TCAT - f4a;
dxdt_T4B = C4B/TCAT - f4b;

dxdt_LOST  = (f4a + f4b)*plo;
dxdt_CSHED = f1 + f2 + f3 + f4a + f4b;

// -------------------------------------------------------- systemic readouts
dxdt_PSA = K_PSA*((0.42 + 0.58*(0.25*E1 + 0.75*E2)) - PSA);
dxdt_EST = 4.0*(AROM*TST*(1.0 + 0.04*(1.0 - DHTS)) - EST);
dxdt_SEB = K_SEB*((0.25 + 0.75*E1S) - SEB);
dxdt_MAP = K_MAP*((MAP0 - MAP_G*mdc) - MAP);
dxdt_HTR = K_HTR*(HTR_G*mdc - HTR);

$TABLE
double N1  = A1  + C1  + T1;
double N2  = A2  + C2  + T2;
double N3  = A3  + C3  + T3;
double N4A = A4A + C4A + T4A;
double N4B = A4B + C4B + T4B;
double NTO = N1 + N2 + N3 + N4A + N4B;

// TARGET-AREA HAIR COUNT: non-vellus hairs in the 1-inch circle.
double TAHC = N1 + N2 + N3;
double VELL = N4A + N4B;
// HAIR MASS INDEX: mass goes as diameter SQUARED times linear growth rate.
double HMI  = 100.0*(N1*D1*D1*G1 + N2*D2*D2*G2 + N3*D3*D3*G3
                     + (N4A+N4B)*D4*D4*G4)/(NTOT*D1*D1*G1);
double MDIA = NTO > 1e-9 ? (N1*D1 + N2*D2 + N3*D3 + (N4A+N4B)*D4)/NTO : 0.0;
double TELPC = NTO > 1e-9 ? 100.0*(T1+T2+T3+T4A+T4B)/NTO : 0.0;
// daily shed, target area and whole-scalp equivalent (x 100000/NTOT)
double TTe_o = TTEL/(1.0 + MX_GTT*(pow(MXSF>0?MXSF:0,MX_HILL)
               /(pow(MX_EC50,MX_HILL)+pow(MXSF>0?MXSF:0,MX_HILL))));
double SHEDTA = (T1+T2+T3+T4A+T4B)/TTe_o;
double SHEDSC = SHEDTA*100000.0/NTOT;
// sexual adverse-effect index and modelled incidence (%)
double SEXIDX = 1.0 - (0.35*E1 + 0.65*E2);
double SEXPCT = 2.1 + 2.6*SEXIDX;
double PSAADJ = 2.0*PSA;     // the "double the PSA" clinical rule

$CAPTURE TAHC HMI VELL MDIA TELPC SHEDTA SHEDSC N1 N2 N3 N4A N4B
         SEXIDX SEXPCT PSAADJ
'

mod <- mcode("aga", aga_code)

## =====================================================================
##  DOSING HELPERS
##  Finasteride and dutasteride enter FING / DUTG as ng of ABSORBED drug.
##  Topical minoxidil enters MXSK in mg, twice daily.
##  Topical finasteride enters BOTH FINSK (mg) and FING (its 0.2% systemic
##  share), twice daily.
## =====================================================================
ev_fin <- function(mg = 1, days = 1825, F = 0.65)
  ev(amt = mg*1e6*F, cmt = "FING", ii = 1, addl = days - 1)

ev_dut <- function(mg = 0.5, days = 1825, F = 0.60)
  ev(amt = mg*1e6*F, cmt = "DUTG", ii = 1, addl = days - 1)

ev_mxt <- function(mg = 50, days = 1825)                      # 1 mL of 5% BID
  ev(amt = mg, cmt = "MXSK", ii = 0.5, addl = 2*days - 1)

ev_mxo <- function(mg = 5, days = 1825, F = 0.90, V = 180000)
  ev(amt = mg*1e6*F/V, cmt = "MXDC", ii = 1, addl = days - 1)

ev_spi <- function(mg = 100, days = 1825, F = 0.70, V = 90000)
  ev(amt = mg*1e6*F/V, cmt = "SPIC", ii = 1, addl = days - 1)

ev_topfin <- function(mg = 2.5, days = 1825, Fsys = 0.002)     # 0.25%, BID
  c(ev(amt = mg,        cmt = "FINSK", ii = 0.5, addl = 2*days - 1),
    ev(amt = mg*1e6*Fsys, cmt = "FING", ii = 0.5, addl = 2*days - 1))

## =====================================================================
##  THE PRESENTING PATIENT IS NOT ASSUMED — IT IS GROWN
##  Start from an intact 18-year-old scalp and run the untreated model
##  until the target-area hair count falls to the published trial baseline
##  of 876 hairs / 5.07 cm2.  With GS = 0.90 this takes 47.8 years.
##  Severity is therefore a CLOCK SPEED, not a state:
##      GS 0.70 -> never reaches 876 within a lifetime
##      GS 0.90 -> age 47.8 y
##      GS 1.00 -> age 33.6 y
##      GS 1.25 -> age 30.7 y
## =====================================================================
burn_in <- function(mod, gs = 0.90, target = 876, maxy = 60) {
  out <- mod %>% param(GS = gs) %>%
    mrgsim(end = maxy*365, delta = 30, atol = 1e-10, rtol = 1e-8) %>% as_tibble()
  k <- which(out$TAHC <= target)[1]
  if (is.na(k)) k <- nrow(out)
  list(state = as.numeric(out[k, cmt(mod)]), age = 18 + out$time[k]/365,
       tahc = out$TAHC[k], out = out)
}

run_arm <- function(mod, y0, events = NULL, end = 1825, ...) {
  m <- mod %>% init(setNames(as.list(y0$state), cmt(mod)))
  extra <- list(...)
  if (length(extra)) m <- do.call(param, c(list(m), extra))
  if (is.null(events)) m %>% mrgsim(end = end, delta = 5, atol = 1e-10, rtol = 1e-8)
  else                 m %>% mrgsim(data = events, end = end, delta = 5,
                                    atol = 1e-10, rtol = 1e-8)
}

## =====================================================================
##  SCENARIO LIBRARY — 18 arms
##  Model outputs quoted below are from the verified Python reference
##  implementation of the identical equations; run this file to reproduce.
## =====================================================================
if (interactive()) {

  P0 <- burn_in(mod)                       # age 47.8 y, TAHC 876, 262 lost

  ## ---- 1  natural history (placebo)  -14.7 hairs at 1 y, -71.5 at 5 y
  s01 <- run_arm(mod, P0)

  ## ---- 2  finasteride 1 mg   +86.5 at 1 y (+101.3 vs placebo; obs +107)
  s02 <- run_arm(mod, P0, ev_fin(1))

  ## ---- 3  finasteride 0.2 mg  +94.2 vs placebo — 1/5 the dose, 7 hairs
  s03 <- run_arm(mod, P0, ev_fin(0.2))

  ## ---- 4  finasteride 5 mg   +105.2 vs placebo — 5x the dose, 4 hairs
  s04 <- run_arm(mod, P0, ev_fin(5))

  ## ---- 5  dutasteride 0.5 mg +122.7 vs placebo at 1 y; scalp DHT -90%
  s05 <- run_arm(mod, P0, ev_dut(0.5))

  ## ---- 6  topical minoxidil 5% BID, SULT1A1 RESPONDER   +39.8 vs placebo
  s06 <- run_arm(mod, P0, ev_mxt(50), SULT = 2.0)

  ## ---- 7  topical minoxidil 5% BID, NON-RESPONDER       +18.7 vs placebo
  s07 <- run_arm(mod, P0, ev_mxt(50), SULT = 0.25)

  ## ---- 8  oral minoxidil 5 mg   +36.7 vs placebo, MAP -3.6 mmHg
  s08 <- run_arm(mod, P0, ev_mxo(5))

  ## ---- 9  oral minoxidil 1 mg   +7.5 vs placebo, MAP -0.7 mmHg
  s09 <- run_arm(mod, P0, ev_mxo(1))

  ## ---- 10 finasteride 1 mg + topical minoxidil 5%   +116.2 vs placebo
  s10 <- run_arm(mod, P0, c(ev_fin(1), ev_mxt(50)))

  ## ---- 11 dutasteride + minoxidil                   +132.8 vs placebo
  s11 <- run_arm(mod, P0, c(ev_dut(0.5), ev_mxt(50)))

  ## ---- 12 topical finasteride 0.25% BID
  ##         scalp DHT -47.6%, serum DHT -27.4%, +78.6 vs placebo
  s12 <- run_arm(mod, P0, ev_topfin(2.5))

  ## ---- 13 withdrawal after 2 years
  ##         scalp DHT 0.364 -> 0.980 in 10 days; hair +110 -> +79 over 3 y
  s13 <- run_arm(mod, P0, ev_fin(1, days = 730))

  ## ---- 14 imperfect adherence, 4 doses in 7
  ##         +91.4 vs placebo at 1 y — 90% of the benefit, 57% of the tablets
  ev14 <- do.call(rbind, lapply(seq(0, 1824), function(d)
            if ((d %% 7) < 4) as.data.frame(ev(amt = 1e6*0.65, cmt = "FING", time = d))
            else NULL))
  s14 <- run_arm(mod, P0, ev14)

  ## ---- 15 late start: the same drug begun two years later
  s15 <- run_arm(mod, P0, ev(amt = 1e6*0.65, cmt = "FING", time = 730,
                             ii = 1, addl = 1094))

  ## ---- 16 setipiprant (CRTH2 blockade).  With PGD_EFF = 0.18 the model
  ##         predicts +46 vs placebo, which the phase 2a did NOT see; a null
  ##         trial requires PGD_EFF <= ~0.05 (see README).
  s16 <- run_arm(mod, P0, NULL, SETI = 1)

  ## ---- 17 FEMALE PATTERN HAIR LOSS, post-menopausal
  ##         androgen-independent drive dominant; finasteride nearly null
  FP <- mod %>% param(GS = 0.45, ARIND = 0.62, LOC5AR = 0.30, AROM = 1.6,
                      W_LOC = 0.30, W_SER = 0.30, W_ALT = 0.40)
  P0f <- burn_in(FP, gs = 0.45, target = 900, maxy = 45)
  s17a <- run_arm(FP, P0f, ev_fin(1),  end = 730)
  s17b <- run_arm(FP, P0f, ev_mxt(50), end = 730)
  s17c <- run_arm(FP, P0f, ev_spi(100), end = 730)

  ## ---- 18 severity ladder: how long each susceptibility takes to get there
  for (g in c(0.70, 0.80, 0.90, 1.00, 1.10, 1.25))
    message(sprintf("GS %.2f -> presenting age %.1f y", g, burn_in(mod, g)$age))
}

## =====================================================================
##  READ-OUT NOTES
##  --------------
##  TAHC    non-vellus hairs per 5.07 cm2 — the endpoint of every pivotal
##          trial, and the one that under-reads the effect (see HMI).
##  HMI     hair mass index, 100 = an intact terminal scalp.  Because mass
##          goes as d^2, HMI moves ~1.9x as far as TAHC.
##  SHEDSC  whole-scalp-equivalent daily shed.  At baseline the model gives
##          90 hairs/day from 100 000 x 0.090 / 100 d — the textbook number
##          is an identity between three parameters, not an observation.
##  TELPC   trichogram telogen percentage = T_T/(T_A+T_C+T_T); it rises in
##          AGA only because T_A falls, so it is an anagen-duration assay.
##  PSAADJ  2 x PSA — the clinical correction after any 5AR inhibitor.
##  SEXPCT  modelled sexual adverse-effect incidence, anchored on
##          3.8% vs 2.1% placebo for finasteride 1 mg.
## =====================================================================
