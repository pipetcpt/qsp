## =============================================================================
##  myp_mrgsolve_model.R
##  Progressive (school) myopia and its pharmacological / optical control
##  Quantitative Systems Pharmacology model for mrgsolve
##
##  52 ODE compartments.  TIME UNIT = DAYS, because this disease forces a
##  2-minute tear-film half-life, a 5-day choroidal time constant, a 20-day
##  scleral creep time constant, a 90-day muscarinic receptor up-regulation
##  time constant and a 10-year growth trajectory to live in one system.
##
##  ---------------------------------------------------------------------------
##  THE ONE STRUCTURAL COMMITMENT
##  ---------------------------------------------------------------------------
##
##  Refraction is NOT a state variable.  It is an OUTPUT computed from three
##  independently regulated optical elements:
##
##      SER = optics( AXIAL LENGTH , CORNEAL POWER , CRYSTALLINE LENS POWER )
##                    ^^^^^^^^^^^^   ^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^
##                    element 1      element 2        element 3
##                    grows under    static after     falls 0.2-0.4 D/yr
##                    scleral creep  age 3            through childhood
##
##  The optics block is an exact paraxial two-thin-lens vergence trace, not a
##  fitted D-per-mm constant, and this is what makes the following OUTPUTS
##  rather than assumptions:
##
##    * dSER/dAL comes out at -2.87 D/mm at AL 23 and -2.58 D/mm at AL 30 at
##      the SPECTACLE plane.  At the CORNEAL plane the same eye falls from
##      -2.91 to -1.73, a 31% loss of sensitivity between AL 24 and AL 29 --
##      but the vertex-distance conversion very nearly cancels it (7%), which
##      is why the empirical "2.7 D per mm" rule of thumb works so well for
##      spectacle refraction and fails for anyone reasoning in contact-lens
##      or ocular terms;
##
##    * emmetropia in a normal growing child is not stasis, it is
##      CANCELLATION.  +0.10 mm/yr of axial growth is -0.28 D/yr of myopic
##      shift, and it is almost exactly offset by -0.30 D/yr of crystalline
##      lens power loss.  Take the lens term out and every emmetrope becomes
##      myopic;
##
##    * every myopia-control treatment ever licensed acts on element 1 and
##      only on element 1.  Element 3 keeps moving regardless, so it BIASES
##      any refractive endpoint.  Applying the model's own optics to the four
##      published LAMP arms recovers an implied lens compensation of
##      +0.353 +/- 0.053 D/yr that is INDEPENDENT of treatment arm (CV 15%) --
##      the internal-consistency test the two LAMP endpoints have to pass.
##      The same test applied ACROSS trials fails: reconciling ATOM1's placebo
##      arm with ATOM2's 0.01% arm requires the latter's lens to have lost
##      0.40 D/yr more power, which 0.8 D of accommodative loss cannot buy.
##
##  ---------------------------------------------------------------------------
##  TWO FURTHER COMMITMENTS
##  ---------------------------------------------------------------------------
##
##  (a) TWO ACTUATORS SHARE ONE MEASUREMENT.  A biometer measures cornea to
##      retinal pigment epithelium.  The choroid sits in that path, changes by
##      tens of microns in DAYS and is fully reversible; the sclera changes by
##      millimetres over MONTHS and is not.  So
##
##          AL_measured = AL_scleral - (CHT - CHT_baseline)/1000
##
##      and the model MEASURES what fraction of each treatment's apparent
##      first-year benefit is the reversible one: 6.3% for atropine 0.01%,
##      7.6% for 0.05%, 13.7% for 1%, 8.4% for DIMS, 8.8% for
##      orthokeratology, but 15.8% for repeated low-level red light.  The
##      model did NOT, however, confirm the author's expectation that this
##      makes red-light rebound visibly FASTER: in the 60 days after
##      withdrawal red light delivers 21% of its washout-year progression
##      and atropine 0.5% delivers 26%, so the choroidal step is too small
##      next to the sustained scleral term to separate the two by SHAPE.
##      Only the SIZE of the rebound differs, and that prediction stands.
##
##  (b) ATROPINE HITS THREE SITES WITH THREE DIFFERENT APPARENT AFFINITIES.
##      The iris sphincter (mydriasis) and the ciliary muscle (cycloplegia)
##      are anterior and heavily/lightly pigmented respectively; the
##      growth-control site is posterior and LOW affinity.  One drug, three
##      dose-response curves, which is why the therapeutic index cannot be
##      fixed by dose selection alone and why the model's efficacy-per-unit-
##      mydriasis ratio is WORST at 0.01%.
##
##  ---------------------------------------------------------------------------
##  WHAT IS FITTED  (13 constants; everything else is physiology or structure)
##  ---------------------------------------------------------------------------
##    KEXP     emmetropic axial growth at age 7 = 0.100 mm/yr
##    G0       LAMP placebo year 1 = 0.410 mm/yr
##    TAUAGE   the same child in year 4 (age 13.2) = 0.200 mm/yr
##    ED50ATR  \
##    HATR      >  atropine 0.01% / 0.05% / 1% axial arms (0.025% HELD OUT)
##    KDIR     /
##    EMAXRL   RLRL 1-year axial arm
##    KRUP     ATOM2 0.5% washout-year refraction (0.01%, 0.1% PREDICTED)
##    EMAXPD   LAMP pupil diameter
##    X1IRIS   LAMP pupil dose-response shape
##    X1CIL    LAMP accommodative-amplitude dose-response
##    AVI,BVI  Tideman 2016 lifetime visual-impairment strata (4 strata)
##  = 13 fitted constants.  KDA (the dopamine / outdoor-light brake) is NOT on
##  that list: it was set a priori and the He 2015 outdoor-time trial was then
##  used as a CHECK, which the model FAILS -- it predicts a 26% reduction in
##  3-year axial elongation from +40 min/day outdoors where that trial measured
##  roughly 11% on refraction.  The outdoor arm is the least trustworthy part
##  of this model and is labelled as such rather than tuned away.
##
##  ---------------------------------------------------------------------------
##  VERIFICATION
##  ---------------------------------------------------------------------------
##  No R runtime was available in the environment this model was written in, so
##  all 52 equations were FIRST implemented and integrated in dependency-free
##  Python RK4 (myp_reference_model.py) and calibrated there
##  (myp_calibration.py, output in myp_calibration_output.txt).  That exposed
##  four real defects, each fixed and recorded at the point of the fix:
##    D1  linear effector couplings drove TGF-beta and then TIMP-2 negative,
##        MMP2/TIMP2 reached 49, scleral creep reached 125 and the eye reached
##        -246 D in three years.  Every coupling is now a POWER LAW, which
##        cannot change sign;
##    D2  measured axial length was referenced to the analytic baseline
##        choroid rather than the settled one, so t = 0 did not reproduce the
##        prescribed baseline refraction (-1.78 D instead of -1.50 D);
##    D3  RK4 intermediate stages do not respect positivity, and a negative
##        base under a fractional exponent returns a COMPLEX number in Python;
##    D4  the atropine dose-response was written as a plain hyperbola in dose
##        (Hill = 1) and could NOT fit the 0.01 / 0.05 / 1% arms at once --
##        the fit pushed two constants onto their bounds and still left the
##        0.05% arm 9 percentage points short.  The published triplet implies
##        a Hill slope near 1.2 in BOTH dose intervals.
##
##  ---------------------------------------------------------------------------
##  A NOTE ON THE PK BLOCK
##  ---------------------------------------------------------------------------
##  Compartments 1-11 are the real topical/systemic atropine PK and are
##  integrated so that dosing events behave naturally and so that tissue and
##  plasma exposure are simulable (the model reproduces plasma C_max of
##  ~8 pg/mL at 0.01% and ~765 pg/mL at 1%, against a measured range of
##  300-900 pg/mL at 1%).  The PD, however, reads the DOSE-LEVEL periodic
##  steady-state summaries (ED50ATR / X1IRIS / X1CIL below) rather than the
##  instantaneous compartment amounts.  This is deliberate and it is what was
##  verified: the clinically measured side-effect readouts are taken ~14 h
##  after a nightly drop, i.e. they are a time-AVERAGE of a nonlinear
##  occupancy, and mean(Occ(C)) is not Occ(mean C).  The posterior/efficacy
##  site is the exception -- it sits far below its EC50, so occupancy there is
##  near-linear in dose and the two formulations coincide.  Do not rewire the
##  PD onto the PK compartments without re-deriving the summaries.
##
##  Usage:
##    library(mrgsolve)
##    mod <- mread("myp_mrgsolve_model.R")
##    # untreated progressing myope, 10 years
##    mod %>% mrgsim(end = 3650, delta = 7) %>% plot(SER + AL + CHT + CREEPS ~ time)
##    # atropine 0.05% nightly for 2 years then stop (see the scenarios at
##    # the bottom of this file)
## =============================================================================

$PROB
# Progressive myopia QSP model
- 52 ODE compartments; time in DAYS
- refraction is an OUTPUT of an exact paraxial optics block, not a state
- choroid (days, reversible) and sclera (months, permanent) share one
  measured axial length
- atropine acts at three sites with three different apparent affinities

$PARAM @annotated
// ---------------- patient descriptors ----------------------------------
AGE0    :  9.7  : baseline age (yr)
SER0    : -3.00 : baseline cycloplegic spherical equivalent (D) [informational]
NPAR    :  2    : number of myopic parents (0/1/2)
GRS     :  0.70 : polygenic risk score (0-1, population mean 0.5)
ETHN    :  1    : 1 = East Asian ancestry, 0 = European
OUTD    :  1.0  : time outdoors (h/day)
NEARD   :  3.0  : near work (h/day)
FUIRISR :  1.0  : iris-pigmentation multiplier on anterior free drug (1 = dark)
IOP0    : 15.5  : intra-ocular pressure set point (mmHg)

// ---------------- treatment inputs -------------------------------------
ATROPCT  : 0.0   : atropine strength (% w/v); 0.01, 0.025, 0.05, 0.1, 0.5, 1.0
TRTSTART : 0.0   : day treatment begins
ATROSTOP : 1e9   : day atropine stops
ATROTAP  : 0.0   : days over which atropine is tapered (0 = abrupt stop)
OPTSTOP  : 1e9   : day the optical/device treatment stops
TRTPD    : 0.15  : imposed peripheral defocus from spectacles/CL (D; -ve = myopic)
                 : SV +0.15, PAL -0.30, MiSight -0.90, DIMS -1.20, HAL -2.20
OKON     : 0     : orthokeratology on (0/1)
RLRLON   : 0     : repeated low-level red light on (0/1)
MXDOSE   : 0     : 7-methylxanthine (mg/day, oral; 400 in the Danish series)
ADHFIX   : -1    : fix adherence at this value; -1 = use the dynamic loop

// ---------------- optics -----------------------------------------------
NVIT    : 1.336  : vitreous refractive index
NAQ     : 1.336  : aqueous refractive index
VERTEX  : 0.012  : spectacle vertex distance (m)
CRAD0   : 7.80   : corneal radius at baseline (mm)
ACD0    : 3.60   : anterior chamber depth at baseline (mm)
LT0     : 3.45   : lens thickness at baseline (mm)
PLENS0  : 22.60  : crystalline lens equivalent power at age 7 (D)
AL0     : 24.291 : baseline axial length (mm); see solve_AL0() in the Python file
KCRAD   : 5.48e-6 : corneal radius drift (mm/day) = 0.002 mm/yr
KACD    : 3.29e-5 : ACD drift (mm/day) = 0.012 mm/yr
KLTAGE  : -5.48e-5 : lens thickness drift (mm/day) = -0.020 mm/yr
TAUAS   : 30     : anterior-segment adaptation time constant (day)
KACDA   : 0.16   : ACD deepening at full ciliary occupancy (mm)
KLTA    : -0.14  : lens thinning at full ciliary occupancy (mm)

// ---------------- crystalline lens (optical element 3) -----------------
KPL     : 1.370e-3 : lens power loss at age 7 (D/day) = 0.50 D/yr
TAUPL   : 7.0    : decay of the lens-power-loss rate (yr)
KMC     : 0.05   : myopia switches lens compensation off (per D)
KPLA    : 0.0    : extra atropine-driven lens flattening (D/day at full effect)

// ---------------- light / dopamine -------------------------------------
LUXREF  : 14500  : reference daily light dose (lux.h)
TAUDA   : 2.0    : retinal dopamine time constant (day)
KDA     : 0.90   : gain of the dopamine brake on the growth drive
PDA     : 0.35   : exponent on light dose

// ---------------- constitutive myopigenic drive ------------------------
G0      : 0.394985 : FITTED -- scale of the constitutive drive
KPAR    : 0.35   : per myopic parent
KGRS    : 1.20   : per unit polygenic risk score above the mean
KETH    : 0.60   : East Asian ancestry
PNEAR   : 0.50   : exponent on near-work hours
NEARREF : 3.0    : reference near work (h/day)

// ---------------- signed defocus response (ONE curve, all devices) -----
KH      : 0.55   : gain, hyperopic defocus -> GO
DH      : 1.20   : saturation of the GO limb (D)
KM      : 1.40   : gain, myopic defocus -> STOP (stronger than GO)
DM      : 0.80   : saturation of the STOP limb (D)
WPER    : 0.70   : weight on peripheral retina
WCEN    : 0.30   : weight on central retina
KSIG    : 1.00   : defocus-to-drive scale
RPR0    : -0.30  : relative peripheral refraction of an emmetrope (D)
KRPR    : 0.25   : relative peripheral hyperopia per D of central myopia
PRPR    : 0.80   : exponent, central myopia -> peripheral hyperopia
TAURPR  : 60     : eye-shape adaptation time constant (day)
ALAG0   : 0.65   : accommodative lag at near, untreated (D)
KLAGM   : 0.10   : extra lag per D of myopia
KLAGA   : 1.60   : extra lag at full ciliary occupancy (D)

// ---------------- effector cascade (all couplings are power laws) ------
TAURA   : 3.0    : retinoic acid time constant (day)
ARA     : 0.55   : drive -> retinoic acid gain
TAUNO   : 2.0    : nitric oxide time constant (day)
ANO     : 0.45   : drive -> NO gain (negative direction)
TAUTG   : 14     : TGF-beta time constant (day)
ATG     : 0.60   : retinoic acid -> TGF-beta exponent (suppressive)
TAUMM   : 10     : MMP-2 time constant (day)
AMM     : 0.70   : retinoic acid -> MMP-2 exponent
AMMH    : 0.40   : scleral hypoxia -> MMP-2 exponent
TAUTI   : 14     : TIMP-2 time constant (day)
ATI     : 0.90   : TGF-beta -> TIMP-2 exponent
KGS     : 0.060  : aggrecan/proteoglycan synthesis rate (1/day)
KGD     : 0.060  : aggrecan degradation rate (1/day)
AGD     : 0.80   : MMP2:TIMP2 exponent on aggrecan degradation
KCS     : 0.030  : collagen I synthesis rate (1/day)
KCD     : 0.030  : collagen I degradation rate (1/day)
ACOD    : 0.70   : MMP2:TIMP2 exponent on collagen degradation
TAULX   : 45     : lysyl-oxidase crosslink time constant (day)
ALX     : 0.70   : TGF-beta -> crosslink exponent
TAUSF   : 20     : myofibroblast index time constant (day)
ASF     : 0.80   : TGF-beta -> myofibroblast exponent
TAUHY   : 15     : scleral hypoxia time constant (day)
AHY     : 0.90   : choroidal perfusion -> scleral hypoxia exponent

// ---------------- scleral creep (the final common path) ----------------
TAUCR   : 20     : creep-state time constant (day)
PC1     : 0.50   : MMP2:TIMP2 exponent
PC2     : 0.35   : aggrecan exponent (protective)
PC3     : 0.50   : collagen I exponent (protective)
PC4     : 0.25   : crosslink exponent (protective)
PC5     : 0.30   : hypoxia exponent

// ---------------- growth -----------------------------------------------
KEXP    : 1.84749e-4 : FITTED -- emmetropic axial growth (mm/day) = 0.0674 mm/yr
TAUAGE  : 3.67856 : FITTED -- plasticity envelope decay (yr)
PEXP    : 1.50   : creep exponent on growth
PIOP    : 0.50   : IOP exponent on growth
IOPREF  : 15.0   : reference IOP (mmHg)
KTHIN   : 0.55   : scleral stretch-thinning coefficient
KDEGS   : 0.06   : scleral thinning from matrix turnover (um/day)
SCT0    : 1000   : baseline posterior scleral thickness (um)
KSTAPH  : 2.0e-4 : staphyloma formation rate
ALSTAPH : 26.5   : axial length above which staphyloma accrues (mm)
KLACQ   : 8.0e-4 : lacquer-crack / patchy-atrophy accrual rate

// ---------------- choroid ----------------------------------------------
CHT0    : 300    : choroidal thickness reference at AL 23.5 mm (um)
TAUCH   : 5.0    : choroidal time constant (day)
KCH     : 0.10   : drive -> choroidal thinning
KCHA    : 0.055  : atropine -> choroidal thickening
KCHR    : 0.075  : red light -> choroidal thickening
KCHAL   : 25     : choroidal thinning per mm of axial length (um/mm)
TAUCBF  : 7.0    : choroidal perfusion time constant (day)

// ---------------- atropine PD ------------------------------------------
EMAXATR : 1.00   : STRUCTURAL -- full blockade of the modifiable drive
ED50ATR : 0.0792 : FITTED -- dose halving the drive (% w/v)
HATR    : 1.2469 : FITTED -- Hill slope of the axial dose-response
KDIR    : 0.5448 : FITTED -- direct scleral suppression at full effect
X1IRIS  : 0.1880 : FITTED -- iris  C_avg/Kd at 0.01% (mydriasis)
X1CIL   : 0.0617 : FITTED -- ciliary C_avg/Kd at 0.01% (cycloplegia)
PUPD0   : 4.60   : photopic pupil diameter, untreated (mm)
EMAXPD  : 3.80   : FITTED -- maximal mydriasis (mm)
AAMP0   : 13.40  : amplitude of accommodation, untreated (D)
TAURUP  : 90     : muscarinic receptor up-regulation time constant (day)
KRUP    : 0.6757 : FITTED -- receptor up-regulation at full blockade

// ---------------- other treatments -------------------------------------
TAUOK   : 5.0    : orthokeratology corneal remodelling time constant (day)
KOKPD   : 1.80   : peripheral myopic defocus imposed by full ortho-K (D)
KOKCR   : 0.42   : central corneal flattening from full ortho-K (mm)
TAUPBA  : 1.5    : fast photobiomodulation time constant (day)
TAUPBC  : 30     : slow photobiomodulation time constant (day)
EMAXRL  : 0.4383 : FITTED -- red-light suppression of the growth drive
KAMX    : 6.0    : 7-MX absorption rate (1/day)
KELMX   : 2.4    : 7-MX elimination rate (1/day)
VMX     : 40     : 7-MX volume of distribution (L)
EMAXMX  : 0.38   : maximal 7-MX suppression of the growth drive
EC50MX  : 1.2    : 7-MX EC50 (umol/L)

// ---------------- adherence loop ---------------------------------------
ADH0    : 0.92   : baseline adherence
TAUADH  : 60     : adherence adaptation time constant (day)
KADHP   : 0.25   : photophobia penalty on adherence
KADHN   : 0.35   : near-blur penalty on adherence

// ---------------- atropine PK (exposure readouts; see the header note) --
MWATR   : 289.37 : atropine molecular weight (g/mol)
KDRAIN  : 500    : tear turnover (1/day), t1/2 = 2 min
KCOR    : 12     : tear -> cornea (1/day)
KTSC    : 1.5    : tear -> conjunctiva -> sclera, periocular route (1/day)
KCA     : 30     : cornea -> aqueous (1/day)
KAO     : 14.4   : aqueous outflow (1/day)
KAI     : 20     : aqueous -> iris/ciliary body (1/day)
KIA     : 8      : iris -> aqueous (1/day)
KAV     : 1.2    : aqueous -> vitreous (1/day)
KVA     : 0.6    : vitreous -> aqueous (1/day)
KVC     : 2.0    : vitreous -> choroid (1/day)
KCS2    : 6.0    : choroid -> sclera (1/day)
KSC     : 3.0    : sclera -> choroid (1/day)
KCP     : 24     : choroid -> plasma (1/day)
KSL     : 4.0    : sclera -> systemic (1/day)
KELP    : 5.5    : plasma elimination (1/day), t1/2 = 3 h
KCPER   : 2.0    : plasma -> peripheral (1/day)
KPERC   : 1.5    : peripheral -> plasma (1/day)
FSYS    : 0.08   : nasolacrimal systemic bioavailability
VPLA    : 30     : plasma volume of distribution (L)
VCHOR   : 0.10   : choroid volume (mL)
VSCLD   : 0.35   : sclera volume (mL)
VIRISD  : 0.30   : iris/ciliary volume (mL)

// ---------------- risk model (Tideman 2016 logistic) -------------------
AVI     : -2.060 : FITTED -- intercept, lifetime visual impairment
BVI     :  0.950 : FITTED -- slope per mm of axial length
KMMD    : 1.1e-5 : myopic maculopathy hazard scale (1/day)
BMMD    : 0.95   : maculopathy hazard exponent (per mm above 26.5)
KRD     : 2.2e-6 : retinal detachment hazard scale (1/day)
BRD     : 0.42   : detachment hazard exponent (per D of myopia)
KCNV    : 3.0e-6 : myopic CNV hazard scale (1/day)
BCNV    : 1.05   : CNV hazard exponent
KGLC    : 1.4e-5 : open-angle glaucoma hazard scale (1/day)
BGLC    : 0.30   : glaucoma hazard exponent
KCAT    : 2.0e-5 : nuclear cataract hazard scale (1/day)
BCAT    : 0.22   : cataract hazard exponent

$CMT @annotated
// ---- atropine PK (1-11) ----
ATEAR   : atropine, precorneal tear film (nmol)
ATCOR   : atropine, cornea (nmol)
ATAQ    : atropine, aqueous humour (nmol)
ATIRIS  : atropine, iris / ciliary body (nmol)
ATVIT   : atropine, vitreous (nmol)
ATCHOR  : atropine, choroid / RPE (nmol)
ATSCL   : atropine, sclera (nmol)
ATPLA   : atropine, plasma (nmol)
ATPER   : atropine, peripheral tissue (nmol)
RBIRIS  : iris muscarinic receptor occupancy (fraction)
RBCIL   : ciliary muscarinic receptor occupancy (fraction)
// ---- retinal / RPE signal cascade (12-18) ----
DA      : retinal dopamine tone (relative, 1 = normal)
RA      : RPE/choroid all-trans retinoic acid (relative)
NOSIG   : retinal + choroidal nitric oxide signal (relative)
TGFB    : scleral TGF-beta1/2 (relative)
MMP2    : active scleral MMP-2 (relative)
TIMP2   : scleral TIMP-2 (relative)
SHYP    : scleral hypoxia / HIF-1alpha index (relative)
// ---- sclera (19-25) ----
GAG     : scleral aggrecan / proteoglycan content (relative)
COL1    : scleral collagen I content (relative)
LOX     : scleral lysyl-oxidase crosslink density (relative)
SFIB    : scleral myofibroblast index (relative)
CREEPS  : scleral creep rate (relative, 1 = emmetropic reference)
SCT     : posterior scleral thickness (um)
PSTAPH  : posterior staphyloma index
// ---- choroid (26-29) ----
LACQ    : lacquer cracks / patchy chorioretinal atrophy index
CHT     : subfoveal choroidal thickness (um)
CBF     : choroidal perfusion (relative)
CHOX    : choroidal hypoxia index (relative)
// ---- biometry (30-35) ----
ALS     : SCLERAL axial length (mm) -- not what a biometer reports
ACD     : anterior chamber depth (mm)
LT      : lens thickness (mm)
PLENS   : crystalline lens equivalent power (D)
CRAD    : corneal radius of curvature (mm)
RPRS    : relative peripheral refraction (D)
// ---- receptor / behaviour (36-39) ----
MRUP    : muscarinic receptor up-regulation (relative, 1 = normal)
ADH     : adherence (fraction of prescribed drops instilled)
AAMPS   : amplitude of accommodation (D)
NEARWS  : near-work burden state (h/day)
// ---- device states (40-44) ----
OKEPI   : orthokeratology corneal epithelial remodelling (0-1)
PBMA    : fast photobiomodulation signal (0-1)
PBMC    : slow photobiomodulation signal (0-1)
MXGUT   : 7-methylxanthine, gut (umol)
MXPLA   : 7-methylxanthine, plasma (umol/L)
// ---- mechanics and cumulative outputs (45-52) ----
IOPS    : intra-ocular pressure (mmHg)
AXCUM   : cumulative SCLERAL axial elongation from baseline (mm)
SERAUC  : integrated myopia exposure (D.yr)
HMMD    : cumulative hazard, myopic maculopathy
HRD     : cumulative hazard, retinal detachment
HCNV    : cumulative hazard, myopic choroidal neovascularisation
HGLC    : cumulative hazard, open-angle glaucoma
HCAT    : cumulative hazard, nuclear cataract

$GLOBAL
#define YR 365.0

// Power with a clamped base.  DEFECT D3: the solver evaluates the right-hand
// side at intermediate states that do not respect positivity, and a negative
// base under a fractional exponent is undefined.  Every power-law coupling in
// the cascade goes through this.
double pw(double b, double e) {
  return pow(b > 1e-6 ? b : 1e-6, e);
}

// Exact paraxial two-thin-lens eye -> SPECTACLE-plane spherical equivalent.
// Vergence trace, backwards from the retina:
//   L2req = n_vit / d2                (vergence needed at the lens plane)
//   L2    = L2req - P_lens            (remove the lens)
//   L1    = L2 / (1 + d1*L2/n_aq)     (transfer back to the corneal plane)
//   F_c   = L1 - P_cornea             (ocular refraction, corneal plane)
//   SER   = F_c / (1 + v*F_c)         (refer to the spectacle plane)
double ser_of(double AL, double ACDx, double LTx, double PL, double CR,
              double naq, double nvit, double vtx) {
  double PCORN = (naq - 1.0) / (CR / 1000.0);
  double d1 = ACDx + LTx / 2.0;
  double d2 = AL - d1;
  if (d2 < 5.0) d2 = 5.0;
  double L2req = 1000.0 * nvit / d2;
  double L2 = L2req - PL;
  double L1 = L2 / (1.0 + (d1 / 1000.0 / naq) * L2);
  double FC = L1 - PCORN;
  return FC / (1.0 + vtx * FC);
}

// Signed, saturating, ASYMMETRIC defocus response.  D > 0 is hyperopic
// defocus (image behind the retina -> GO); D < 0 is myopic defocus (STOP).
// The STOP limb is deliberately the stronger one -- brief myopic defocus
// dominates prolonged hyperopic defocus in every animal model.  ONE curve
// serves the untreated eye and all five optical treatments.
double respD(double D, double kh, double dh, double km, double dm) {
  if (D > 0.0) return  kh * D / (1.0 + D / dh);
  double a = -D;
  return -km * a / (1.0 + a / dm);
}

$MAIN
// ---- initial conditions ----------------------------------------------
// Note on AL0: the Python reference model solves the baseline axial length
// from the prescribed baseline refraction by bisection on the optics block
// (solve_AL0), so that t = 0 reproduces SER0 exactly.  Pass the solved value
// in as AL0.  For SER0 = -3.00 D with the default anterior segment it is
// 24.291 mm; -1.50 D -> 23.759; -1.00 D -> 23.583; -0.75 D -> 23.495;
// +0.50 D -> 23.058; -4.70 D -> 24.898.
ALS_0    = AL0;
ACD_0    = ACD0;
LT_0     = LT0;
PLENS_0  = PLENS0;
CRAD_0   = CRAD0;
CHT_0    = CHT0 - KCHAL * (AL0 > 23.50 ? AL0 - 23.50 : 0.0);
SCT_0    = SCT0;
RPRS_0   = RPR0 + KRPR * pw(SER0 < 0 ? -SER0 : 0.0, PRPR);
DA_0     = 1.0;  RA_0    = 1.0;  NOSIG_0 = 1.0;  TGFB_0 = 1.0;
MMP2_0   = 1.0;  TIMP2_0 = 1.0;  SHYP_0  = 1.0;
GAG_0    = 1.0;  COL1_0  = 1.0;  LOX_0   = 1.0;  SFIB_0 = 1.0;
CREEPS_0 = 1.0;  CBF_0   = 1.0;  CHOX_0  = 1.0;  MRUP_0 = 1.0;
ADH_0    = ADH0;
AAMPS_0  = AAMP0;
NEARWS_0 = NEARD;
IOPS_0   = IOP0;

// The measured axial length is referenced to the choroidal thickness at
// baseline.  DEFECT D2: this must be the SETTLED choroid, so run the model
// with ATROPCT = 0 and no device for ~200 days first if you need t = 0 to
// reproduce SER0 to the third decimal place.
double CHTREF = CHT0 - KCHAL * (AL0 > 23.50 ? AL0 - 23.50 : 0.0);

$ODE
// =====================================================================
//  0.  time-dependent treatment inputs
// =====================================================================
double tnow = SOLVERTIME;
double pct = 0.0;
if (tnow >= TRTSTART) {
  if (tnow < ATROSTOP)            pct = ATROPCT;
  else if (ATROTAP > 0.0)         pct = ATROPCT * (1.0 - (tnow - ATROSTOP) / ATROTAP);
  if (pct < 0.0) pct = 0.0;
}
double devon = (tnow >= TRTSTART && tnow < OPTSTOP) ? 1.0 : 0.0;
double AGE = AGE0 + tnow / YR;
double adh = (ADHFIX >= 0.0) ? ADHFIX : ADH;

// =====================================================================
//  1.  atropine occupancies and effects
//      Three sites, three apparent affinities.  See the header note on why
//      the PD reads dose-level summaries rather than the PK compartments.
// =====================================================================
double dose_eff = pct * adh;
double xi  = X1IRIS * FUIRISR * (dose_eff / 0.01);
double xc  = X1CIL  * FUIRISR * (dose_eff / 0.01);
double OCCI = (dose_eff > 0.0) ? xi / (1.0 + xi) : 0.0;   // mydriasis
double OCCC = (dose_eff > 0.0) ? xc / (1.0 + xc) : 0.0;   // cycloplegia
// DEFECT D4: Hill slope, not a plain hyperbola.  1.2 is what the published
// 0.01 / 0.05 / 1% triplet implies, and it is the same in both intervals.
double dh_ = (dose_eff > 0.0) ? pw(dose_eff, HATR) : 0.0;
double EATR = (dose_eff > 0.0)
              ? EMAXATR * dh_ / (dh_ + pw(ED50ATR, HATR)) : 0.0;
double EMX  = EMAXMX * MXPLA / (EC50MX + MXPLA);
double ERL  = EMAXRL * (0.35 * PBMA + 0.65 * PBMC);

// =====================================================================
//  2.  OPTICS -- refraction is computed here, it is never integrated
// =====================================================================
double ALMEAS  = ALS - (CHT - CHTREF) / 1000.0;
double CRADEFF = CRAD + KOKCR * OKEPI;          // ortho-K flattens the cornea
double SERnow  = ser_of(ALMEAS, ACD, LT, PLENS, CRADEFF, NAQ, NVIT, VERTEX);
// ...and the refraction the eye WOULD have without the corneal change.  This
// is why orthokeratology trials cannot use a refractive endpoint at all.
double SERTRUE = ser_of(ALMEAS, ACD, LT, PLENS, CRAD,    NAQ, NVIT, VERTEX);
double MYO = (SERTRUE < 0.0) ? -SERTRUE : 0.0;

// =====================================================================
//  3.  defocus -> drive.  ONE response curve for every optical treatment.
// =====================================================================
double RPRT = RPR0 + KRPR * pw(MYO, PRPR);
double DPER = RPRS + devon * TRTPD - KOKPD * OKEPI;
double ALAG = ALAG0 + KLAGM * MYO + KLAGA * OCCC;
double SIGRAW = WPER * respD(DPER, KH, DH, KM, DM)
              + WCEN * (NEARD / 12.0) * respD(ALAG, KH, DH, KM, DM);

double GBASE = G0 * (1.0 + KPAR * NPAR) * (1.0 + KGRS * (GRS - 0.50))
             * (1.0 + KETH * ETHN) * pw(NEARD / NEARREF, PNEAR);
double FDA = 1.0 / (1.0 + KDA * (DA - 1.0));
double SIG = GBASE + KSIG * SIGRAW;
double SIGEFF = SIG * FDA * (1.0 - EATR) * (1.0 - EMX) * (1.0 - ERL) * MRUP;
if (SIGEFF < -0.90) SIGEFF = -0.90;

// =====================================================================
//  4.  effector cascade -- every coupling is a POWER LAW (defect D1)
// =====================================================================
double LUXH = OUTD * 10000.0 + (16.0 - OUTD) * 300.0;
dxdt_DA    = (pw(LUXH / LUXREF, PDA) - DA) / TAUDA;
dxdt_RA    = (exp(ARA * SIGEFF) - RA) / TAURA;
dxdt_NOSIG = (exp(-ANO * SIGEFF) - NOSIG) / TAUNO;
dxdt_TGFB  = (pw(RA, -ATG) - TGFB) / TAUTG;
dxdt_MMP2  = (pw(RA, AMM) * pw(SHYP, AMMH) - MMP2) / TAUMM;
dxdt_TIMP2 = (pw(TGFB, ATI) - TIMP2) / TAUTI;

double MTR = MMP2 / (TIMP2 > 1e-4 ? TIMP2 : 1e-4);   // MMP-2 : TIMP-2 balance
dxdt_GAG   = KGS * TGFB - KGD * GAG  * pw(MTR, AGD);
dxdt_COL1  = KCS * TGFB - KCD * COL1 * pw(MTR, ACOD);
dxdt_LOX   = (pw(TGFB, ALX) - LOX)  / TAULX;
dxdt_SFIB  = (pw(TGFB, ASF) - SFIB) / TAUSF;
dxdt_SHYP  = (pw(CBF, -AHY) - SHYP) / TAUHY;

// =====================================================================
//  5.  scleral creep -- the final common path for every treatment
// =====================================================================
double CRTGT = pw(MTR, PC1) * pw(GAG, -PC2) * pw(COL1, -PC3)
             * pw(LOX, -PC4) * pw(SHYP, PC5);
dxdt_CREEPS = (CRTGT - CREEPS) / TAUCR;

// =====================================================================
//  6.  axial growth.  PHI(age) is the plasticity envelope, and it is what
//      makes the cumulative benefit of ANY treatment plateau -- both arms
//      share it, so the protectable integral is front-loaded in age.
// =====================================================================
double PHI = exp(-(AGE - 7.0) / TAUAGE);
// Direct scleral limb: ATOM1 measured -0.02 mm over TWO YEARS on 1%
// atropine, which is below the EMMETROPIC growth rate.  Removing the
// myopigenic drive cannot get there; high-dose atropine must also act
// directly on the sclera.  KDIR carries that, fitted to the 1% arm alone.
double gsupp = 1.0 - KDIR * EATR;
if (gsupp < 0.02) gsupp = 0.02;
double dALS = KEXP * PHI * pw(CREEPS, PEXP) * pw(IOPS / IOPREF, PIOP) * gsupp;
dxdt_ALS   = dALS;
dxdt_AXCUM = dALS;

dxdt_SCT = -KTHIN * dALS * 1000.0 * (SCT / SCT0)
           - KDEGS * ((MTR > 1.0) ? MTR - 1.0 : 0.0);
dxdt_PSTAPH = KSTAPH * ((ALS > ALSTAPH) ? ALS - ALSTAPH : 0.0)
              * (SCT0 / (SCT > 100.0 ? SCT : 100.0));
dxdt_LACQ = KLACQ * PSTAPH * ((ALS > ALSTAPH) ? ALS - ALSTAPH : 0.0);

// =====================================================================
//  7.  choroid -- the fast, reversible actuator on the SAME measurement
// =====================================================================
double CHT0E = CHT0 - KCHAL * ((ALS > 23.50) ? ALS - 23.50 : 0.0);
double CHTTGT = CHT0E * (1.0 - KCH * SIGEFF + KCHA * EATR
                         + KCHR * (0.5 * PBMA + 0.5 * PBMC));
if (CHTTGT < 40.0) CHTTGT = 40.0;
dxdt_CHT  = (CHTTGT - CHT) / TAUCH;
dxdt_CBF  = (CHT / (CHT0E > 50.0 ? CHT0E : 50.0) - CBF) / TAUCBF;
dxdt_CHOX = (1.0 / (CBF > 0.2 ? CBF : 0.2) - CHOX) / TAUCH;

// =====================================================================
//  8.  crystalline lens (element 3) and anterior segment
// =====================================================================
double FMYO = 1.0 / (1.0 + KMC * MYO);   // compensation stops at onset
dxdt_PLENS = -KPL * exp(-(AGE - 7.0) / TAUPL) * FMYO - KPLA * EATR;
dxdt_ACD  = ((ACD0 + KACD  * YR * (AGE - 7.0) + KACDA * OCCC) - ACD) / TAUAS;
dxdt_LT   = ((LT0  + KLTAGE * YR * (AGE - 7.0) + KLTA  * OCCC) - LT)  / TAUAS;
dxdt_CRAD = KCRAD;
dxdt_RPRS = (RPRT - RPRS) / TAURPR;

// =====================================================================
//  9.  anterior-segment pharmacodynamics and the adherence loop
// =====================================================================
dxdt_AAMPS = (AAMP0 * (1.0 - OCCC) - AAMPS) / 3.0;
double NEARBLUR = 1.0 - AAMPS / AAMP0;
double ADHT = ADH0 * (1.0 - KADHP * OCCI - KADHN * NEARBLUR);
if (ADHT < 0.15) ADHT = 0.15;
dxdt_ADH    = (ADHT - ADH) / TAUADH;
dxdt_NEARWS = (NEARD - NEARWS) / 30.0;
// Chronic blockade up-regulates receptors, so withdrawal is supersensitive.
// This is the SLOW half of rebound; the choroidal collapse is the fast half.
dxdt_MRUP   = (1.0 + KRUP * EATR - MRUP) / TAURUP;

// =====================================================================
// 10.  device states
// =====================================================================
dxdt_OKEPI = (OKON   * devon - OKEPI) / TAUOK;
dxdt_PBMA  = (RLRLON * devon - PBMA)  / TAUPBA;
dxdt_PBMC  = (RLRLON * devon - PBMC)  / TAUPBC;
dxdt_MXGUT = -KAMX * MXGUT;
dxdt_MXPLA = KAMX * MXGUT / VMX - KELMX * MXPLA;
dxdt_IOPS  = (IOP0 - IOPS) / 30.0;

// =====================================================================
// 11.  atropine PK (exposure readouts -- see the header note)
// =====================================================================
dxdt_ATEAR  = -(KDRAIN + KCOR + KTSC) * ATEAR;
dxdt_ATCOR  = KCOR * ATEAR - KCA * ATCOR;
dxdt_ATAQ   = KCA * ATCOR - (KAO + KAI + KAV) * ATAQ + KIA * ATIRIS + KVA * ATVIT;
dxdt_ATIRIS = KAI * ATAQ - KIA * ATIRIS;
dxdt_ATVIT  = KAV * ATAQ - (KVA + KVC) * ATVIT;
dxdt_ATCHOR = KVC * ATVIT + KSC * ATSCL - (KCP + KCS2) * ATCHOR;
dxdt_ATSCL  = KTSC * ATEAR + KCS2 * ATCHOR - (KSC + KSL) * ATSCL;
dxdt_ATPLA  = FSYS * KDRAIN * ATEAR + KCP * ATCHOR + KSL * ATSCL
              - (KELP + KCPER) * ATPLA + KPERC * ATPER;
dxdt_ATPER  = KCPER * ATPLA - KPERC * ATPER;
// slowly reversible anterior receptor complexes (free drug drives them; the
// iris is pigmented, so FUIRISR gates how much of the iris drug is free)
double CIRF = 0.02 * FUIRISR * ATIRIS / VIRISD * 1000.0;   // nM
dxdt_RBIRIS = 0.0234 * CIRF * (1.0 - RBIRIS) - 0.35 * RBIRIS;
dxdt_RBCIL  = 0.0714 * CIRF * (1.0 - RBCIL)  - 0.35 * RBCIL;

// =====================================================================
// 12.  cumulative risk hazards
// =====================================================================
dxdt_SERAUC = MYO / YR;
dxdt_HMMD = KMMD * exp(BMMD * (ALMEAS - 26.5));
dxdt_HRD  = KRD  * exp(BRD * MYO);
dxdt_HCNV = KCNV * exp(BCNV * ((ALMEAS > 26.0) ? ALMEAS - 26.0 : 0.0));
dxdt_HGLC = KGLC * exp(BGLC * MYO);
dxdt_HCAT = KCAT * exp(BCAT * MYO);

$TABLE
double AGEy = AGE0 + TIME / YR;
double CHTREFo = CHT0 - KCHAL * (AL0 > 23.50 ? AL0 - 23.50 : 0.0);
double AL  = ALS - (CHT - CHTREFo) / 1000.0;         // what a biometer reads
double CRADo = CRAD + KOKCR * OKEPI;
double SER = ser_of(AL, ACD, LT, PLENS, CRADo, NAQ, NVIT, VERTEX);
double SERTR = ser_of(AL, ACD, LT, PLENS, CRAD, NAQ, NVIT, VERTEX);
double ALCR = AL / CRAD;                             // >3.0 signals myopia

// atropine exposure readouts
double CPLA_PG = ATPLA / VPLA * 1000.0 * MWATR;      // pg/mL
double CCHOR_NM = ATCHOR / VCHOR * 1000.0;           // nM
double CSCL_NM  = ATSCL  / VSCLD * 1000.0;           // nM

// anterior-segment effects
double adh2 = (ADHFIX >= 0.0) ? ADHFIX : ADH;
double pct2 = 0.0;
if (TIME >= TRTSTART) {
  if (TIME < ATROSTOP)     pct2 = ATROPCT;
  else if (ATROTAP > 0.0)  pct2 = ATROPCT * (1.0 - (TIME - ATROSTOP) / ATROTAP);
  if (pct2 < 0.0) pct2 = 0.0;
}
double de2 = pct2 * adh2;
double xi2 = X1IRIS * FUIRISR * (de2 / 0.01);
double OCCI2 = (de2 > 0.0) ? xi2 / (1.0 + xi2) : 0.0;
double PUPD = PUPD0 + EMAXPD * OCCI2;
double AAMP = AAMPS;
double dh2 = (de2 > 0.0) ? pw(de2, HATR) : 0.0;
double EATRo = (de2 > 0.0) ? EMAXATR * dh2 / (dh2 + pw(ED50ATR, HATR)) : 0.0;

// structural risk
double VIRISK = 1.0 / (1.0 + exp(-(AVI + BVI * (AL - 26.0))));  // lifetime
double PMMD = 1.0 - exp(-HMMD);
double PRD  = 1.0 - exp(-HRD);
double PCNV = 1.0 - exp(-HCNV);
double PGLC = 1.0 - exp(-HGLC);
double PCAT = 1.0 - exp(-HCAT);
double HIGHMYO = (SERTR <= -6.0) ? 1.0 : 0.0;
double AL26 = (AL >= 26.0) ? 1.0 : 0.0;
double MTRo = MMP2 / (TIMP2 > 1e-4 ? TIMP2 : 1e-4);

$CAPTURE @annotated
AGEy    : age (yr)
SER     : cycloplegic spherical equivalent, spectacle plane (D)
SERTR    : refraction with the ortho-K corneal change removed (D)
AL      : MEASURED axial length, cornea to RPE (mm)
ALCR    : axial length / corneal radius ratio
PUPD    : photopic pupil diameter (mm)
AAMP    : amplitude of accommodation (D)
EATRo   : fractional suppression of the growth drive by atropine
MTRo    : MMP-2 : TIMP-2 balance
CPLA_PG : atropine plasma concentration (pg/mL)
CCHOR_NM : atropine choroid concentration (nM)
CSCL_NM : atropine sclera concentration (nM)
VIRISK  : lifetime risk of irreversible visual impairment (Tideman logistic)
PMMD    : cumulative probability, myopic maculopathy
PRD     : cumulative probability, retinal detachment
PCNV    : cumulative probability, myopic CNV
PGLC    : cumulative probability, open-angle glaucoma
PCAT    : cumulative probability, nuclear cataract
HIGHMYO : 1 if the eye has reached -6.00 D
AL26    : 1 if the eye has reached 26 mm

## =============================================================================
##  TREATMENT SCENARIOS  (24)
## =============================================================================
##  Atropine is a nightly 30 uL drop into the tear compartment.  Because the
##  PD reads the dose-level summaries (see the header), the drop events are
##  needed only for the PK/exposure readouts -- ATROPCT drives the PD.
##
##  library(mrgsolve); library(dplyr)
##  mod <- mread("myp_mrgsolve_model.R")
##
##  drops <- function(pct, days = 3650)                     # 30 uL nightly
##    ev(time = 0, amt = pct * 10 * 30 / 289.37 * 1000,
##       cmt = "ATEAR", ii = 1, addl = days - 1)
##
##  run <- function(..., pct = 0, days = 3650)
##    mod %>% param(ATROPCT = pct, ...) %>%
##      mrgsim(events = drops(pct, days), end = days, delta = 7)
##
##  ---- natural history --------------------------------------------------
##   1  emmetropic child, age 7        param(AGE0=7, SER0=0.5, AL0=23.058,
##                                            NPAR=0, GRS=0.5, ETHN=0)
##   2  untreated myope, age 8 -1.5 D  param(AGE0=8, SER0=-1.5, AL0=23.759)
##   3  untreated myope, age 9.7 -3 D  param()                     [reference]
##   4  early onset, age 6 -0.75 D     param(AGE0=6, SER0=-0.75, AL0=23.495)
##   5  late onset, age 12 -1 D        param(AGE0=12, SER0=-1.0, AL0=23.583)
##   6  European, slow progressor      param(ETHN=0, GRS=0.5, NPAR=1)
##
##  ---- atropine ---------------------------------------------------------
##   7  atropine 0.01%                 run(pct = 0.01)
##   8  atropine 0.025%                run(pct = 0.025)
##   9  atropine 0.05%   (LAMP choice) run(pct = 0.05)
##  10  atropine 0.1%                  run(pct = 0.10)
##  11  atropine 1%      (ATOM1)       run(pct = 1.00)
##  12  0.05%, 2 yr then abrupt stop   run(pct=0.05, ATROSTOP=730)
##  13  0.5%,  2 yr then abrupt stop   run(pct=0.50, ATROSTOP=730)
##  14  0.5%,  2 yr then 1-yr taper    run(pct=0.50, ATROSTOP=730, ATROTAP=365)
##  15  0.05%, adherence fixed at 0.5  run(pct=0.05, ADHFIX=0.5)
##  16  0.05% in a light iris          run(pct=0.05, FUIRISR=4)
##  17  0.05% started 2 years late     run(pct=0.05, TRTSTART=730)
##
##  ---- optical and device ----------------------------------------------
##  18  DIMS spectacles                param(TRTPD=-1.20)
##  19  MiSight dual-focus SCL         param(TRTPD=-0.90)
##  20  HAL lenslet spectacles         param(TRTPD=-2.20)
##  21  orthokeratology                param(OKON=1, TRTPD=0)
##  22  ortho-K + atropine 0.01%       run(pct=0.01, OKON=1, TRTPD=0)
##  23  repeated low-level red light   param(RLRLON=1)
##  24  RLRL, 2 yr then stop           param(RLRLON=1, OPTSTOP=730)
##
##  ---- environment and pitfalls ---------------------------------------
##  25  +40 min/day outdoors           param(OUTD=1.67)
##  26  heavy near work                param(NEARD=5.0)
##  27  UNDER-correction by 0.75 D     param(TRTPD=0.90)
##  28  7-methylxanthine 400 mg/day    param(MXDOSE=400)  + a daily gut event
## =============================================================================
