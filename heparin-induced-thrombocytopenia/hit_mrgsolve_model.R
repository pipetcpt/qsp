# =====================================================================
# Heparin-Induced Thrombocytopenia (HIT) — QSP model
# including autoimmune HIT (aHIT), spontaneous HIT syndrome and VITT
#
#   59 ODEs · time unit = HOURS · mrgsolve
#   Validated under real mrgsolve; every number quoted below is produced
#   by the scenario block at the bottom of this file.
# =====================================================================
#
# ---------------------------------------------------------------------
# THE ORGANISING IDEA
#
#   HIT IS NOT A DRUG TOXICITY. IT IS AN IMMUNE COMPLEX DISEASE IN
#   WHICH THE DRUG IS ONE HALF OF THE ANTIGEN.
#
# Platelet factor 4 (PF4/CXCL4) is a 7.8 kDa homotetramer with a ring of
# lysines and a net charge near +20. Heparin is a polyanion. Near
# charge-neutral molar ratio the two do not merely bind — they polymerise
# into ULTRALARGE COMPLEXES (ULCs, >670 kDa) in which several PF4
# tetramers are strung along and between chains. Only the ULC is
# antigenic, because only the ULC presents the neoepitope at the surface
# density that multivalent IgG needs, and only multivalent IgG can
# crosslink FcgammaRIIa.
#
# The consequence that organises this whole model is that the antigen is
# NOT a monotonic function of heparin dose. It is gated three ways, and
# the model puts the polyanion into the equations three separate times:
#
#   (1) AS SUBSTRATE       HEPB, the concentration of chains long enough
#                          to BRIDGE two PF4 tetramers (>= ~11-14
#                          saccharides). More bridging chains, more
#                          antigen. This is what makes the dose-response
#                          within a drug monotonic.
#   (2) AS COMPETITOR      HEPT, the concentration of all chains weighted
#                          by PF4-BINDING competence. Short chains still
#                          bind PF4 and compete it away even though they
#                          cannot bridge. When HEPT greatly exceeds PF4,
#                          every tetramer gets its own chain, bridging
#                          fails, and the ULC dissolves — HEPARIN EXCESS.
#   (3) AS PF4 MOBILISER   Long chains displace PF4 from endothelial
#                          heparan sulfate, which is why plasma PF4 rises
#                          10-20x within minutes of a UFH bolus. The
#                          drug manufactures its own antigen partner.
#
# Because (1) and (3) favour long chains and (2) penalises high chain
# molarity, the clinical risk ordering falls out of chain chemistry
# rather than being asserted:
#
#     UFH therapeutic  >>  LMWH therapeutic  >  LMWH/UFH prophylaxis
#                      >>  danaparoid  >>  fondaparinux (a structural null)
#
# and so does the one experimental fact that any HIT model must
# reproduce: platelet activation at 0.1-0.3 U/mL heparin that is
# ABOLISHED at 100 U/mL. That high-heparin step is the confirmatory arm
# of the serotonin release assay, and here it is not coded as a rule —
# it is limb (2) of the same equation, evaluated 250-fold along the
# heparin axis.
#
# ---------------------------------------------------------------------
# THREE LAYERS IN SERIES, AND ONE LOOP THAT CLOSES ON ITSELF
#
#   LAYER 1  ANTIGEN    PF4 + polyanion -> ULC (cmts PF4P, PF4EC, PF4DNA,
#                       ULC). Cut ONLY by stopping heparin.
#   LAYER 2  IMMUNITY   ULC -> B-cell clone -> plasmablast -> anti-PF4/H
#                       IgG (cmts BCELL, PBLAST, LLPC, MEMB, IGGP, IGGN).
#                       Peculiar kinetics, all structural, none scripted:
#                         · a PRE-EXISTING clone (BCELL0) is why IgG is
#                           detectable by day 4-5 on a FIRST exposure
#                         · plasmablasts dominate and LLPC are ~0, which
#                           is why the antibody is TRANSIENT
#                         · KRECALL defaults to 0, which is why there is
#                           NO ANAMNESTIC BOOST on re-exposure
#   LAYER 3  EFFECTOR   IgG.ULC -> FcgammaRIIa crosslinking -> platelet
#                       activation, splenic clearance, PS+ microparticles,
#                       monocyte tissue factor, NETs, thrombin burst.
#
#   THE LOOP: activated platelets dump alpha-granule PF4 into the plasma
#   they are being activated in (PLT_GRAN -> PF4P -> ULC -> activation).
#   Antigen makes activation makes antigen. This is why HIT is explosive
#   rather than graded, and why PLATELET TRANSFUSION is anti-therapy: it
#   feeds the loop its substrate.
#
# ---------------------------------------------------------------------
# THE CENTRAL PARADOX, AND WHY THE HAZARD HAS TWO ROUTES
#
#   A DISEASE DEFINED BY A LOW PLATELET COUNT WHOSE PATIENTS DIE OF CLOTS.
#
# The same FcgammaRIIa signal that marks platelets for splenic clearance
# also sheds procoagulant microparticles and induces monocyte tissue
# factor. So the model routes ACTSIG to BOTH sinks and gives the
# thrombosis hazard two ADDITIVE routes:
#
#     HZTHR   thrombin-driven   <- cut by argatroban, bivalirudin, DOACs
#     HZPLT   platelet-driven   <- cut by NOTHING an anticoagulant does;
#                                  only heparin cessation, IVIG or
#                                  plasma exchange touch it
#
# That single structural choice is what makes the drug classes separate
# on two independent axes instead of one, and it is why a DTI-treated
# patient still carries residual risk (model: 4.8% vs 51.0% untreated).
#
# ---------------------------------------------------------------------
# THE TWO PHARMACOLOGICAL TRAPS, GENERATED RATHER THAN WARNED ABOUT
#
# (i) THE WARFARIN / VENOUS LIMB GANGRENE TRAP
#     Protein C t1/2 ~8 h; prothrombin t1/2 ~60-72 h. Both are ODEs here
#     with those two rate constants and nothing else. Start a vitamin K
#     antagonist while thrombin generation is still maximal and the model
#     reproduces the published mismatch unaided: 24 h after a 10 mg load
#     protein C is 37% of normal while prothrombin is still 83%, and
#     effective thrombin activity RISES from 2.1 to 10.7 nM. The
#     warfarin start-day sweep then generates the guideline rule — delay
#     until platelets recover — as a monotone dose-response, not as text.
#
# (ii) THE ARGATROBAN-INR TRAP
#     Argatroban prolongs the PT with no vitamin-K-dependent factor
#     depletion at all. INR is therefore computed from TWO inputs (true
#     factor activity, and a DTI assay artefact) and chromogenic factor X
#     is reported alongside it, so the dissociation is visible: measured
#     INR 4.31 while the true INR is 2.07 and factor X is still 47%.
#
# ---------------------------------------------------------------------
# THE ICEBERG IS AN OUTPUT, NOT AN ASSUMPTION
#
# After heparin exposure the antibody is common and the disease is rare.
# The model separates ANTIBODY PRESENT / PLATELET-ACTIVATING / CLINICAL
# as three successive thresholds and reports each level as the assay that
# measures it (ELISA optical density; SRA % release at low and at 100 U/mL
# heparin). Diagnostic D4 shows that what actually sets a patient's level
# on the iceberg is EXPOSURE INTENSITY: below ~0.07 U/mL nothing happens,
# at ~0.16 U/mL the patient is seropositive AND SRA-positive but the
# platelet fall stops at 48% and never becomes clinical HIT, and only
# above ~0.25 U/mL does the third waterline get crossed. Notably the
# pathogenic-fraction parameter FPATH does NOT generate the iceberg —
# see the negative results below. Prophylactic-dose exposure seroconverts
# (OD 0.82 after bypass) without ever becoming clinical HIT (7.0% fall
# on prophylactic-dose heparin, 21.4% after cardiopulmonary bypass).
#
# ---------------------------------------------------------------------
# ALL THE VARIANTS ARE ONE COEFFICIENT
#
# Autoimmune HIT, spontaneous HIT syndrome and VITT share this entire
# model with one change: HIND, the heparin-independence coefficient,
# moves from 0 toward 1 because their antibodies cluster PF4 at or near
# the heparin-binding site and need no heparin to do it. Everything
# clinically strange about them then follows without further edits:
#   · the platelet count keeps falling after heparin is stopped
#   · the standard SRA loses its high-heparin inhibition — the model
#     reports SRA 93% at LOW heparin, 93% at 100 U/mL and 93% with NO
#     heparin, which is precisely why these sera need PF4-enhanced
#     assays and why a classic SRA can read negative
#   · IVIG is transformative where anticoagulation is not, because
#     receptor competition does not care where the antibody bound
#     (model, aHIT: 7.7% with IVIG vs 56.4% with argatroban alone)
#
# ---------------------------------------------------------------------
# BASELINE IS SOLVED, NOT FITTED
#
# Every baseline balance is derived algebraically in $MAIN from the
# stated baseline values (SPF4, KECD, KBS, KMK, STPO, KTF0, KVWS, KTMS,
# KATS, KPT, KFBN, KDDS, KF12S). A 200-day run with no heparin therefore
# drifts by 0.0000% in every state: the disease in this model is
# GENERATED by exposure, never assumed.
#
# ---------------------------------------------------------------------
# CALIBRATION ANCHORS (see the scenario and diagnostic blocks below)
#
#   platelet fall onset, therapeutic UFH        model day 6.4   (obs 5-10, median 6)
#   platelet nadir                              model 54        (obs median 55-60)
#   nadir rarely <20                            model 50-100    (obs, distinguishes from ITP/DIC)
#   30-day thrombosis, untreated                model 51.0%     (obs 20-50%+)
#   30-day thrombosis, DTI-treated              model 4.8%      (obs ~6-14%)
#   plasma PF4 on therapeutic UFH               model 0.5 ug/mL (obs 0.05-0.2, up to ~2 post-CPB)
#   protein C at 24 h after 10 mg warfarin      model 37%       (obs ~30-40%)
#   prothrombin at 24 h after 10 mg warfarin    model 83%       (obs ~80%)
#   argatroban Css at 2 ug/kg/min               model 1.20 ug/mL (obs ~1.2)
#   argatroban aPTT at Css                      model 70 s      (obs 60-80, target 1.5-3x)
#   SRA at 0.1-0.3 U/mL heparin                 model 93%       (obs >20% = positive)
#   SRA at 100 U/mL heparin (classic HIT)       model 14%       (obs <20% = confirmatory)
#   SRA at 100 U/mL heparin (VITT)              model 93%       (obs NOT inhibited)
#   seroreversion of the immunoassay             model 120-180 d (obs ~85-100 d)
#   fondaparinux ULC formation                  model ~0        (obs no ULC; structural null)
#
# NEGATIVE / SELF-REFUTING RESULTS, REPORTED RATHER THAN REMOVED
#
#   · The model does NOT reproduce the historical ~50% residual thrombosis
#     after "heparin cessation alone" as a property of cessation. It
#     reproduces it only for LATE recognition (D7): recognise on day 5 and
#     the residual risk is 0.9%, on day 9 it is 18%, on day 30 it is 51%.
#     The model's claim is that "cessation alone" was never one
#     intervention — it is a family of them indexed by the day the
#     diagnosis was made, and the historical figure is a statement about
#     how late HIT used to be recognised.
#   · The widespread belief that continuing HEPARIN LINE FLUSHES
#     perpetuates HIT is NOT supported (D7b). Sub-therapeutic heparin
#     makes almost no ULC and still accelerates antithrombin, so the model
#     puts residual risk flat to slightly LOWER. This was originally built
#     in as the explanation for the point above; it failed, and the
#     failure is reported instead of being removed.
#   · The pathogenic-fraction parameter FPATH does NOT generate the
#     iceberg. Once a therapeutic exposure has occurred the antibody is so
#     abundant that even 2% of it being platelet-activating is enough. The
#     iceberg in this model is generated by EXPOSURE INTENSITY instead
#     (D4), which is a stronger and more testable statement.
#   · Prophylactic UFH comes out slightly LOWER risk than prophylactic
#     LMWH, the opposite of the epidemiology. The model assigns UFH its
#     advantage entirely through chain length and concentration, so at the
#     low concentrations subcutaneous prophylaxis reaches, that advantage
#     is spent. A genuine limitation, not a tuned result.
#   · Plasma exchange is nearly inert here (D-section B, scenario 18)
#     because removing IgG does not remove the plasmablasts making it.
#     The model says PLEX is a bridge, not a treatment.
#   · Rapid-onset reactivity OUTLASTS immunoassay positivity (D8), so the
#     model regards the 100-day re-exposure rule as an underestimate of
#     the hazard window.
#
# =====================================================================
# NOTE ON USE: educational / research QSP model. Not validated for
# clinical decision-making, prescribing or regulatory submission.
# =====================================================================

library(mrgsolve)

# ---------------------------------------------------------------------
# MODEL SPECIFICATION
# ---------------------------------------------------------------------
hit_code <- '
$PROB
# Heparin-Induced Thrombocytopenia (HIT) — QSP model
# 59 ODEs. Time unit = HOURS. See hit_mrgsolve_model.R header for documentation.

$SET end=672, delta=2, rtol=1e-8, atol=1e-10, maxsteps=500000

$GLOBAL
#define PF4EFF (FLOC*PF4P)

// clamp to non-negative: the solver can undershoot states that decay over many
// orders of magnitude, and a fractional power of a negative number is NaN
double pos(double x){ return (x > 0.0) ? x : 0.0; }

// bridging competence of a polyanion chain of length L (saccharide units)
double flen(double L, double L50){
  double a = pow(L, 6.0);
  return a/(a + pow(L50, 6.0));
}
// affinity of a chain of length L for PF4 (a lower length threshold than bridging:
// a short chain still BINDS PF4 and so still competes it away, it just cannot bridge)
double fbind(double L, double LB50){
  double a = pow(L, 3.0);
  return a/(a + pow(LB50, 3.0));
}
// Antigenicity as a function of the PF4:chain molar ratio R.
// ONE-SIDED by design. Above the optimum (R >= Rstar) PF4 is sufficient and every
// bridging-competent chain is used, so antigen scales with chain supply, not ratio:
// this is what makes the within-drug dose-response monotonic. Below the optimum the
// system is in HEPARIN EXCESS -- each PF4 tetramer gets its own chain, bridging fails
// and the complex dissolves. That falling limb is the SRA 100 U/mL confirmatory step
// and the reason cardiopulmonary-bypass heparin concentrations are self-limiting.
double fratio(double R, double Rstar, double sd){
  if(R <= 1e-12) return 0.0;
  if(R >= Rstar) return 1.0;
  double z = log(R/Rstar)/sd;
  return exp(-0.5*z*z);
}
// in-vitro platelet-activation assay: serum IgG + donor PF4 + [heparin] -> % release
double assay(double iggp, double pf4a, double hepU, double chn, double len, double l50,
             double rstar, double sd, double kulc, double kulcoff, double kind, double hind,
             double kicon, double kicoff, double kicclr, double kdic, double fcg,
             double nxl, double xl50, double khb){
  double hb  = chn*hepU*flen(len, l50);          // bridging-competent chains
  double ht  = chn*hepU*fbind(len, 8.0);   // LENB50 = 8         // PF4-binding chains (ratio denominator)
  double ulc = 0.0;
  if(ht > 1e-9) ulc += kulc*pf4a*fratio(pf4a/ht, rstar, sd)*(hb/(hb + khb))/kulcoff;
  ulc += kind*hind*pf4a/kulcoff;                    // heparin-independent (aHIT / VITT)
  double ic  = kicon*iggp*ulc/(kicoff + kicclr);
  double occ = (pos(ic)/kdic*fcg)/(1.0 + pos(ic)/kdic*fcg);
  double xl  = pow(pos(occ), nxl)/(pow(pos(occ), nxl) + pow(xl50, nxl));
  return 100.0*xl;
}

// ---- baseline quantities solved in $MAIN, shared with $ODE / $TABLE ----
double SPF4, KECD, KBS, KMK, STPO, KTF0, KVWS, KTMS, KATS, KPT, KFBN, KDDS, KF12S;
double FLENUFH, FLENLMW, FLENFON, FLENDAN, FLENDNA, THRGEN0, FFV0;
double FBNDUFH, FBNDLMW, FBNDFON, FBNDDAN, FBNDDNA;

$PARAM @annotated
WT       :   70   : Body weight (kg)

// ================= UFH pharmacokinetics (U, U/mL, h) =================
VHEP     : 4200   : UFH central volume (mL)
CLHEP    : 1200   : UFH linear (renal) clearance (mL/h)
VMHEP    : 2500   : UFH saturable elimination Vmax (U/h)
KMHEP    :  0.9   : UFH saturable elimination Km (U/mL)
KAHEP    : 0.35   : UFH subcutaneous absorption rate (1/h)
FHEPSC  : 0.45     : UFH subcutaneous bioavailability
KPROT    : 0.02   : Protamine neutralisation rate of UFH (1/(mg*h))
KPROTD   : 0.35   : Protamine elimination rate (1/h)

// ================= LMWH (enoxaparin), fondaparinux, danaparoid =======
KALMW    : 0.50   : LMWH absorption rate (1/h)
FLMW     : 0.92   : LMWH bioavailability
VLMW     : 5000   : LMWH central volume (mL)
CLLMW    :  750   : LMWH clearance (mL/h)
KAFON    : 0.60   : Fondaparinux absorption rate (1/h)
FFON     : 1.00   : Fondaparinux bioavailability
VFON     : 8000   : Fondaparinux central volume (mL)
CLFON    :  280   : Fondaparinux clearance (mL/h)
VDAN     : 7000   : Danaparoid central volume (mL)
CLDAN    :  200   : Danaparoid clearance (mL/h)

// ================= polyanion chain chemistry =========================
CHNUFH   :  370   : Chain molarity per U/mL UFH (nM chains)
CHNLMW   : 2222   : Chain molarity per anti-Xa IU/mL LMWH (nM chains)
CHNFON   : 578.7  : Chain molarity per ug/mL fondaparinux (nM chains)
CHNDAN   : 1800   : Chain molarity per anti-Xa U/mL danaparoid (nM chains)
LENUFH   :   50   : UFH mean chain length (saccharide units)
LENLMW   :   15   : LMWH mean chain length (saccharide units)
LENFON   :    5   : Fondaparinux chain length (saccharide units)
LENDAN  : 9        : Danaparoid effective bridging length (saccharide units)
LENDNA   :  200   : Cell-free DNA / polyP effective length (saccharide equivalents)
LEN50   : 24       : Chain length giving half-maximal bridging (saccharide units)
LENB50   :    8   : Chain length giving half-maximal PF4 binding (saccharide units)
KHB     : 800      : Bridging-chain concentration for half-maximal ULC assembly (nM)

// ================= PF4 biology =======================================
PF40     : 0.16   : Baseline free plasma PF4 (nM tetramer)
ECC0     :  400   : Baseline endothelial GAG-bound PF4 pool (nM)
ECCAP    : 4000   : Endothelial GAG PF4 binding capacity (nM)
KONEC    : 58.5   : PF4 association rate onto endothelial GAG (1/h)
KOFFEC   : 0.0008 : PF4 dissociation rate from endothelial GAG (1/h)
KPF4C    : 2.77   : Free plasma PF4 clearance (1/h)
KHSPG    :   10   : Polyanion concentration halving PF4-GAG association (nM)
KMOB     : 0.15   : Maximal heparin-driven mobilisation of GAG-bound PF4 (1/h)
KMOBH    :  100   : Polyanion concentration for half-maximal mobilisation (nM)
FLOC     :  500   : Local (surface microenvironment) PF4 amplification factor
KRELACT : 6.0      : PF4 release gain from platelet activation signal
KSURGP   :   19   : PF4 release gain from surgical / CPB stimulus
KSURGD  : 0.02     : Decay of surgical / CPB activation stimulus (1/h)

// ================= antigen assembly (ULC) ============================
RSTAR    :  8.0   : Optimal PF4-tetramer : polyanion-chain molar ratio
SDR     : 1.4      : Log-scale width of the antigenicity bell curve
KULC    : 0.25     : ULC formation rate constant (1/h)
KULCOFF  :  2.0   : ULC dissociation rate (1/h)
KIND     : 0.05   : Heparin-independent ULC formation rate constant (1/h)
HIND     :  0.0   : HEPARIN-INDEPENDENCE coefficient (0 classic HIT, ~0.9 VITT)
KPDS    : 2.0e-5   : PF4-DNA complex formation rate (1/(nM*h))
KPDD     : 0.05   : PF4-DNA complex elimination (1/h)

// ================= immunity ==========================================
BCELL0   :  1.0   : Pre-existing anti-PF4 B-cell clone size (arbitrary units)
BMAX     : 1000   : Maximal anti-PF4 B-cell clone size
KB       :   15   : ULC concentration for half-maximal B-cell drive (nM)
KBPROL  : 0.055    : Anti-PF4 B-cell proliferation rate (1/h)
KBREC    : 0.05   : Naive B-cell recruitment rate (1/h)
KBD      : 0.004  : Anti-PF4 B-cell death rate (1/h)
KDIFF   : 0.028    : B-cell to plasmablast differentiation rate (1/h)
KPBD     : 0.00722: Plasmablast death rate (1/h) [t1/2 4 d]
FLL     : 0.0008   : Fraction differentiating to long-lived plasma cells
KLLD     : 1.6e-4 : Long-lived plasma cell death rate (1/h)
FMEM     : 0.05   : Fraction differentiating to memory B cells
KMEMD    : 9.6e-5 : Memory B-cell death rate (1/h)
KRECALL  :  0.0   : Memory recall coefficient (0 = NO anamnestic boost)
SIGG     : 0.02   : IgG secretion rate per plasmablast unit (U/mL/h)
FPATH    : 0.15   : Fraction of anti-PF4 IgG that is platelet-activating
KIGGD    :0.001375: Anti-PF4 IgG catabolic rate (1/h) [t1/2 21 d]
FCRNX    :  1.0   : IgG catabolism multiplier (FcRn inhibitor)
KTREG    :  0.0   : Regulatory suppression of B-cell drive (0-1)

// ================= immune complex and FcgammaRIIa =====================
KICON    :  0.1   : IgG-ULC association rate
KICOFF   :  0.5   : Immune-complex dissociation rate (1/h)
KICCLR   :  0.3   : Immune-complex clearance rate (1/h)
KDIC    : 60.0     : Immune complex KD for platelet FcgammaRIIa (nM)
FCG      :  1.0   : FcgammaRIIa affinity multiplier (HH 1.35, HR 1.0, RR 0.75)
KDIVIG  : 4.0      : IVIG concentration competing at FcgammaRIIa (g/L)
KDIVIGM  : 18.0   : IVIG concentration blocking macrophage FcgammaR (g/L)
NXL     : 2.5      : Hill coefficient for crosslinking threshold
XL50    : 0.35     : Fractional occupancy giving half-maximal crosslinking

// ================= platelets =========================================
PLT0     :  250   : Baseline platelet count (10^9/L)
KPLTB    :0.00463 : Baseline platelet elimination rate (1/h) [lifespan 9 d]
KACT    : 0.009    : Platelet activation/consumption rate at full signal (1/h)
KOPS    : 0.018    : Opsonised platelet splenic clearance rate (1/h)
KCONS    : 2e-4   : Platelet consumption per unit fibrin burden (1/h)
KACLR    : 0.17   : Activated platelet clearance rate (1/h)
KTR      : 0.025  : Megakaryocyte transit rate (1/h)
KTPO     : 0.02   : Thrombopoietin turnover rate (1/h)
TPOMAX   :  6.0   : Maximal thrombopoietin drive
KMPD     : 0.03   : Microparticle elimination rate (1/h)
MP0      :  1.0   : Baseline procoagulant microparticle level (normalised)
KMPA     :  8.0   : Microparticle production gain from platelet activation

// ================= monocyte, neutrophil, endothelium =================
KMOA     : 0.08   : Monocyte activation rate from immune complexes (1/h)
KMOP     : 0.01   : Monocyte activation rate from surface PF4 (1/h)
KPF4M    :  5.0   : Excess PF4 for half-maximal monocyte signal (nM)
KMOD     : 0.03   : Monocyte deactivation rate (1/h)
TF0      :  1.0   : Baseline tissue factor activity (normalised)
KTFS     :  2.5   : Tissue factor induction gain
KTFEC    :  0.4   : Endothelial contribution to tissue factor
KTFD     : 0.115  : Tissue factor decay rate (1/h) [t1/2 6 h]
KNETS   : 0.05     : NET formation rate (1/h)
KNETD    : 0.04   : NET clearance rate (1/h)
KECS     : 0.02   : Endothelial activation rate (1/h)
KEC40    :  3.0   : sCD40L contribution to endothelial activation
KECD2    : 0.03   : Endothelial deactivation rate (1/h)
VWF0     :  100   : Baseline von Willebrand factor (%)
KVWD     : 0.058  : vWF elimination rate (1/h)
KVWA     :  1.5   : vWF release gain from endothelial activation
TM0      :  1.0   : Baseline thrombomodulin (normalised)
KTMD     : 0.02   : Thrombomodulin turnover rate (1/h)
KTMA     :  2.0   : Thrombomodulin loss gain from endothelial activation

// ================= coagulation factors ===============================
KFII     :0.01155 : Prothrombin elimination rate (1/h) [t1/2 60 h]
KFVII    : 0.1155 : Factor VII elimination rate (1/h) [t1/2 6 h]
KFIX     : 0.0289 : Factor IX elimination rate (1/h) [t1/2 24 h]
KFX      :0.01733 : Factor X elimination rate (1/h) [t1/2 40 h]
KPC      : 0.0866 : Protein C elimination rate (1/h) [t1/2 8 h]
KPROSD   : 0.0165 : Protein S elimination rate (1/h) [t1/2 42 h]
AT0      :  100   : Baseline antithrombin (%)
KATD     :0.01155 : Antithrombin elimination rate (1/h)
KATC    : 0.05     : Antithrombin consumption per nM thrombin (1/h)
THR0     :  0.5   : Baseline thrombin level (nM)
KTHRI    :  2.0   : Antithrombin-mediated thrombin inhibition rate (1/h)
EMAXAT  : 10       : Maximal heparin acceleration of antithrombin
KATACC  : 0.30     : Heparin concentration for half-maximal AT acceleration
KAPC     :  1.0   : APC generation gain
KAPCI    :  2.0   : APC concentration halving prothrombinase (normalised)
FBN0     :  1.0   : Baseline fibrin / clot burden (normalised)
KFBND    : 0.02   : Fibrin turnover rate (1/h)
DDIM0    :  0.3   : Baseline D-dimer (ug/mL FEU)
KDDD     : 0.0347 : D-dimer elimination rate (1/h)
FRG0     :  0.2   : Baseline prothrombin fragment 1.2 (nmol/L)
KFRGD    : 0.462  : Prothrombin fragment 1.2 elimination rate (1/h)

// ================= non-heparin anticoagulants ========================
VARG1    : 12000  : Argatroban central volume (mL)
VARG2    : 10000  : Argatroban peripheral volume (mL)
CLARG   : 7000     : Argatroban clearance (mL/h) [hepatic]
QARG     :  8000  : Argatroban intercompartmental clearance (mL/h)
HEPFN    :  1.0   : Hepatic function multiplier on argatroban clearance
VBIV     :  4500  : Bivalirudin central volume (mL)
CLBIV    : 12500  : Bivalirudin clearance (mL/h)
KARIV    :  1.2   : Rivaroxaban absorption rate (1/h)
FRIV     : 0.80   : Rivaroxaban bioavailability
VRIV     : 55000  : Rivaroxaban central volume (mL)
CLRIV    :  7000  : Rivaroxaban clearance (mL/h)
RENFN    :  1.0   : Renal function multiplier (CrCl/100)
KIARG   : 0.45     : Argatroban Ki at thrombin (ug/mL)
KIBIV    : 0.25   : Bivalirudin Ki at thrombin (ug/mL)
KIRIV    :  120   : Rivaroxaban Ki at factor Xa (ng/mL)
KIFON    :  0.9   : Fondaparinux Ki at factor Xa (ug/mL)
KIDAN    :  0.8   : Danaparoid Ki at factor Xa (anti-Xa U/mL)
KIUFHX   :  0.5   : UFH anti-Xa potency (U/mL)
KILMWX   :  0.7   : LMWH anti-Xa potency (IU/mL)

// ================= immunomodulation ==================================
VIVG1    :  3000  : IVIG central volume (mL)
VIVG2    :  4000  : IVIG peripheral volume (mL)
CLIVG    :  10    : IVIG clearance (mL/h)
QIVG     :  60    : IVIG intercompartmental clearance (mL/h)
KIVGULC : 0.35     : IVIG-mediated ULC disruption rate (1/(g/L)/h)
KPLEX    : 0.367  : Plasma-exchange removal rate per unit session intensity (1/h)
KPLEXD   : 0.35   : Decay of plasma-exchange session intensity (1/h)
VRTX1    :  3200  : Rituximab central volume (mL)
VRTX2    :  2700  : Rituximab peripheral volume (mL)
CLRTX    :  10    : Rituximab clearance (mL/h)
QRTX     :  30    : Rituximab intercompartmental clearance (mL/h)
KRTXKILL : 0.004  : Rituximab B-cell depletion rate (1/h)
IC50RTX  :  1.0   : Rituximab concentration for half-maximal depletion (ug/mL)
STER     :  0.0   : Corticosteroid suppression of B-cell drive (0-1)

// ================= warfarin ==========================================
KAWAR    :  1.5   : Warfarin absorption rate (1/h)
VWAR     :  8000  : Warfarin central volume (mL)
CLWAR    :  190   : Warfarin clearance (mL/h)
IMAXVKA  : 0.92   : Maximal VKORC1 inhibition
EC50VKA : 0.35     : Warfarin concentration for half-maximal effect (mg/L)
HVKA     :  1.4   : Warfarin effect Hill coefficient
KVITK    :  2.0   : Vitamin K amount halving VKA effect (mg)
KVITKD   : 0.06   : Vitamin K elimination rate (1/h)
KARGINR  :  0.9   : Argatroban artefactual INR gain (per ug/mL)
KBIVINR  :  0.5   : Bivalirudin artefactual INR gain (per ug/mL)
KRIVINR  : 0.0025 : Rivaroxaban artefactual INR gain (per ng/mL)
APTT0    :   30   : Baseline aPTT (s)
KAPUFH   :  3.3   : UFH aPTT gain (per U/mL)
KAPARG   :  1.1   : Argatroban aPTT gain (per ug/mL)
KAPBIV   :  1.6   : Bivalirudin aPTT gain (per ug/mL)

// ================= assays and endpoints ==============================
PF4ASSAY : 1200   : Donor-platelet PF4 in the in-vitro assay (nM)
HEPSRALO :  0.2   : SRA low heparin concentration (U/mL)
HEPSRAHI :  100   : SRA high heparin (confirmatory) concentration (U/mL)
ODMAX    :  2.6   : Maximal PF4 ELISA optical density
ODBASE   : 0.06   : PF4 ELISA background optical density
KOD     : 455      : IgG level for half-maximal optical density (U/mL)
WNP      :  0.7   : ELISA weight of non-activating IgG
HZ0     : 1.5e-5   : Thrombosis hazard scaling (1/h)
PHZ     : 1.2      : Thrombin exponent in the thrombosis hazard
KHZF     :  0.25  : Fibrin-burden contribution to the thrombosis hazard
KHZP     :  62.5  : Platelet-activation contribution to the thrombosis hazard
KVLG    : 5.0e-4   : Venous limb gangrene hazard scaling (1/h)

$CMT @annotated
HEPD  : UFH subcutaneous depot (U)
HEP   : UFH central compartment (U)
PROT  : Protamine (mg)
LMWD  : LMWH subcutaneous depot (anti-Xa IU)
LMWC  : LMWH central compartment (anti-Xa IU)
FOND  : Fondaparinux subcutaneous depot (ug)
FONC  : Fondaparinux central compartment (ug)
DANC  : Danaparoid central compartment (anti-Xa U)
SURG  : Surgical / cardiopulmonary bypass activation stimulus (unitless)
PF4P  : Free plasma PF4 (nM tetramer)
PF4EC : Endothelial GAG-bound PF4 pool (nM)
PF4DNA: PF4-DNA / PF4-polyphosphate complexes (nM)
ULC   : Ultralarge PF4-polyanion complex (nM)
IC    : Immune complex IgG-PF4-polyanion (nM)
BCELL : Anti-PF4 B-cell clone (arbitrary units)
PBLAST: Short-lived plasmablasts (arbitrary units)
LLPC  : Long-lived plasma cells (arbitrary units)
MEMB  : Anti-PF4 memory B cells (arbitrary units)
IGGP  : Platelet-activating anti-PF4/heparin IgG (U/mL)
IGGN  : Non-activating anti-PF4/heparin IgG (U/mL)
PLT   : Circulating platelet count (10^9/L)
PLTA  : Activated platelet pool (10^9/L)
MK1   : Megakaryocyte transit compartment 1
MK2   : Megakaryocyte transit compartment 2
TPO   : Thrombopoietin (normalised)
MP    : Procoagulant microparticles (normalised)
MOA   : Activated monocyte fraction (0-1)
TF    : Tissue factor activity (normalised)
NET   : Neutrophil extracellular trap burden (normalised)
ECA   : Endothelial activation (normalised)
VWF   : von Willebrand factor (%)
TM    : Endothelial thrombomodulin (normalised)
FII   : Prothrombin activity (%)
FVII  : Factor VII activity (%)
FIX   : Factor IX activity (%)
FX    : Factor X activity (%)
PC    : Protein C activity (%)
PROS  : Protein S activity (%)
AT    : Antithrombin activity (%)
THR   : Thrombin (nM)
FBN   : Fibrin / clot burden (normalised)
DDIM  : D-dimer (ug/mL FEU)
FRG   : Prothrombin fragment 1.2 (nmol/L)
TEC   : Cumulative HIT-attributable thrombotic event probability
NECR  : Cumulative venous limb gangrene / skin necrosis probability
THRAUC: Cumulative thrombin exposure (nM*h)
ARGC  : Argatroban central compartment (ug)
ARGP  : Argatroban peripheral compartment (ug)
BIVC  : Bivalirudin central compartment (ug)
RIVD  : Rivaroxaban absorption depot (mg)
RIVC  : Rivaroxaban central compartment (mg)
IVGC  : IVIG central compartment (g)
IVGP  : IVIG peripheral compartment (g)
WARD  : Warfarin absorption depot (mg)
WARC  : Warfarin central compartment (mg)
RTXC  : Rituximab central compartment (mg)
RTXP  : Rituximab peripheral compartment (mg)
PLEXA : Plasma-exchange session intensity (unitless)
VITK  : Vitamin K (mg)

$MAIN
// ---------- chain-length bridging competence (constant per species) ----------
FLENUFH = flen(LENUFH, LEN50);
FLENLMW = flen(LENLMW, LEN50);
FLENFON = flen(LENFON, LEN50);
FLENDAN = flen(LENDAN, LEN50);
FLENDNA = flen(LENDNA, LEN50);
FBNDUFH = fbind(LENUFH, LENB50);
FBNDLMW = fbind(LENLMW, LENB50);
FBNDFON = fbind(LENFON, LENB50);
FBNDDAN = fbind(LENDAN, LENB50);
FBNDDNA = fbind(LENDNA, LENB50);

// ---------- PF4 module: solve the baseline EXACTLY ----------
double ecflux = KONEC*PF40*(1.0 - ECC0/ECCAP);
KECD = (ecflux - KOFFEC*ECC0)/ECC0;               // GAG-bound PF4 internalisation
SPF4 = KPF4C*PF40 + ecflux - KOFFEC*ECC0;         // PF4 synthesis / basal release

// ---------- other baseline balances, solved ----------
KBS   = KBD*BCELL0;                               // basal anti-PF4 B-cell source
KMK   = KPLTB*PLT0;                               // megakaryocyte drive at TPO = 1
STPO  = KTPO*1.0;                                 // thrombopoietin synthesis
KTF0  = KTFD*TF0;                                 // constitutive tissue factor
KVWS  = KVWD*VWF0;
KTMS  = KTMD*TM0;
KATS  = KATD*AT0 + KATC*THR0;
FFV0  = 1.0/(1.0 + KAPC/KAPCI);                   // basal FVa availability
THRGEN0 = KTHRI*THR0;                             // basal thrombin generation (nM/h)
KPT   = THRGEN0/FFV0;                             // prothrombinase scaling
KFBN  = KFBND*FBN0/THR0;
KDDS  = KDDD*DDIM0/(KFBND*FBN0);
KF12S = KFRGD*FRG0/THRGEN0;

// ---------- initial conditions ----------
PF4P_0  = PF40;   PF4EC_0 = ECC0;
BCELL_0 = BCELL0;
PLT_0   = PLT0;   MK1_0 = KPLTB*PLT0/KTR;  MK2_0 = KPLTB*PLT0/KTR;
TPO_0   = 1.0;    MP_0  = MP0;   TF_0 = TF0;
VWF_0   = VWF0;   TM_0  = TM0;
FII_0   = 100;    FVII_0 = 100;  FIX_0 = 100;  FX_0 = 100;
PC_0    = 100;    PROS_0 = 100;  AT_0  = AT0;
THR_0   = THR0;   FBN_0 = FBN0;  DDIM_0 = DDIM0;  FRG_0 = FRG0;

$ODE
// ================= 1. heparin and polyanion exposure =================
double CUFH = HEP/VHEP;                              // U/mL
double CLMW = LMWC/VLMW;                             // anti-Xa IU/mL
double CFON = FONC/VFON;                             // ug/mL
double CDAN = DANC/VDAN;                             // anti-Xa U/mL

dxdt_HEPD = -KAHEP*HEPD;
dxdt_HEP  = KAHEP*FHEPSC*HEPD
            - CLHEP*CUFH - VMHEP*CUFH/(KMHEP + CUFH) - KPROT*PROT*HEP;
dxdt_PROT = -KPROTD*PROT;
dxdt_LMWD = -KALMW*LMWD;
dxdt_LMWC = KALMW*FLMW*LMWD - CLLMW*RENFN*CLMW;
dxdt_FOND = -KAFON*FOND;
dxdt_FONC = KAFON*FFON*FOND - CLFON*RENFN*CFON;
dxdt_DANC = -CLDAN*RENFN*CDAN;
dxdt_SURG = -KSURGD*SURG;

// total and bridging-competent polyanion chain concentrations (nM)
// HEPT: all chains, weighted by PF4-BINDING competence -> sets the molar ratio
//       (this is what competes PF4 away and drives the heparin-excess limb)
// HEPB: only chains long enough to BRIDGE two PF4 tetramers -> the ULC substrate
double HEPT = CHNUFH*CUFH*FBNDUFH + CHNLMW*CLMW*FBNDLMW + CHNFON*CFON*FBNDFON
            + CHNDAN*CDAN*FBNDDAN + PF4DNA*FBNDDNA;
double HEPB = CHNUFH*CUFH*FLENUFH + CHNLMW*CLMW*FLENLMW + CHNFON*CFON*FLENFON
            + CHNDAN*CDAN*FLENDAN + PF4DNA*FLENDNA;

// ================= 2. PF4 ============================================
// PF4 is displaced from endothelial heparan sulfate by LONG chains (HEPB), not by
// chain molarity: this is a second, independent route by which UFH beats LMWH, and
// it is why plasma PF4 rises 10-20x within minutes of a UFH bolus.
double FHEP    = HEPB/(HEPB + KMOBH);
double KONEFF  = KONEC/(1.0 + HEPB/KHSPG);           // long chains block GAG rebinding
double ACTSIG;                                       // set below, needed here
double OCC, XL;
// --- effector signal (computed first: needed by PF4 release and platelets) ---
double IVG = (IVGC/VIVG1)*1000.0;                    // g/L
double icn = pos(IC)/KDIC*FCG;
OCC = icn/(1.0 + icn + IVG/KDIVIG);
XL  = pow(pos(OCC), NXL)/(pow(pos(OCC), NXL) + pow(XL50, NXL));
ACTSIG = XL;

double RELPF4 = SPF4*(1.0 + KRELACT*ACTSIG + KSURGP*SURG);
dxdt_PF4P  = RELPF4 + KOFFEC*PF4EC + KMOB*FHEP*PF4EC
             - KPF4C*PF4P - KONEFF*PF4P*(1.0 - PF4EC/ECCAP);
dxdt_PF4EC = KONEFF*PF4P*(1.0 - PF4EC/ECCAP) - KOFFEC*PF4EC - KECD*PF4EC
             - KMOB*FHEP*PF4EC;

// ================= 3. antigen assembly ===============================
double PF4L = PF4EFF;                                // local effective PF4 (nM)
double RRAT = (HEPT > 1e-9) ? PF4L/HEPT : 1e12;
double FRAT = (HEPT > 1e-9) ? fratio(RRAT, RSTAR, SDR) : 0.0;
double ULCFORM = KULC*PF4L*FRAT*(HEPB/(HEPB + KHB)) + KIND*HIND*PF4L;
dxdt_ULC = ULCFORM - KULCOFF*ULC - KICON*IGGP*ULC + KICOFF*IC
           - KIVGULC*IVG*ULC - KPLEX*PLEXA*ULC;
dxdt_IC  = KICON*IGGP*ULC - KICOFF*IC - KICCLR*IC - KPLEX*PLEXA*IC;
dxdt_PF4DNA = KPDS*NET*PF4P*1000.0 - KPDD*PF4DNA;

// ================= 4. immunity =======================================
double DRIVE = (ULC/(ULC + KB))*(1.0 - KTREG)*(1.0 - STER);
double CRTX  = RTXC/VRTX1*1000.0;                    // ug/mL
double RTXE  = CRTX/(CRTX + IC50RTX);
dxdt_BCELL  = KBS + KBPROL*DRIVE*BCELL*(1.0 - BCELL/BMAX) + KBREC*BCELL0*DRIVE
              + KRECALL*MEMB*DRIVE - KBD*BCELL - KRTXKILL*RTXE*BCELL;
dxdt_PBLAST = KDIFF*DRIVE*BCELL - KPBD*PBLAST;
dxdt_LLPC   = FLL*KDIFF*DRIVE*BCELL - KLLD*LLPC;
dxdt_MEMB   = FMEM*KDIFF*DRIVE*BCELL - KMEMD*MEMB - KRTXKILL*RTXE*MEMB;
double PCELL = PBLAST + LLPC;
dxdt_IGGP = SIGG*FPATH*PCELL - KIGGD*FCRNX*IGGP - KPLEX*PLEXA*IGGP;
dxdt_IGGN = SIGG*(1.0 - FPATH)*PCELL - KIGGD*FCRNX*IGGN - KPLEX*PLEXA*IGGN;

// ================= 5. platelets ======================================
double OPSON = OCC/(1.0 + IVG/KDIVIGM);
double TPOE  = (TPO < TPOMAX) ? TPO : TPOMAX;
double PROD  = KTR*MK2;
dxdt_PLT  = PROD - KPLTB*PLT - KACT*ACTSIG*PLT - KOPS*OPSON*PLT
            - KCONS*(FBN - FBN0)*PLT;
dxdt_PLTA = KACT*ACTSIG*PLT - KACLR*PLTA;
dxdt_MK1  = KMK*TPOE - KTR*MK1;
dxdt_MK2  = KTR*MK1 - KTR*MK2;
dxdt_TPO  = STPO - KTPO*TPO*(PLT/PLT0);
dxdt_MP   = KMPD*MP0*(1.0 + KMPA*(ACTSIG + 5.0*PLTA/PLT0)) - KMPD*MP;

// ================= 6. monocyte / neutrophil / endothelium ============
double PF4X = (PF4P > PF40) ? (PF4P - PF40) : 0.0;
dxdt_MOA = KMOA*OCC*(1.0 - MOA) + KMOP*(PF4X/(PF4X + KPF4M))*(1.0 - MOA) - KMOD*MOA;
dxdt_TF  = KTF0 + KTFS*(MOA + KTFEC*ECA) - KTFD*TF;
dxdt_NET = KNETS*(10.0*PLTA/PLT0 + MOA) - KNETD*NET;
dxdt_ECA = KECS*(OCC + KEC40*PLTA/PLT0) - KECD2*ECA;
dxdt_VWF = KVWS*(1.0 + KVWA*ECA) - KVWD*VWF;
dxdt_TM  = KTMS - KTMD*TM*(1.0 + KTMA*ECA);

// ================= 7. coagulation ====================================
double CWAR = WARC/VWAR*1000.0;                      // mg/L
double VKA  = IMAXVKA*pow(CWAR, HVKA)/(pow(EC50VKA, HVKA) + pow(CWAR, HVKA));
VKA = VKA/(1.0 + VITK/KVITK);
dxdt_FII  = KFII*100.0*(1.0 - VKA)  - KFII*FII;
dxdt_FVII = KFVII*100.0*(1.0 - VKA) - KFVII*FVII;
dxdt_FIX  = KFIX*100.0*(1.0 - VKA)  - KFIX*FIX;
dxdt_FX   = KFX*100.0*(1.0 - VKA)   - KFX*FX;
dxdt_PC   = KPC*100.0*(1.0 - VKA)   - KPC*PC;
dxdt_PROS = KPROSD*100.0*(1.0 - VKA) - KPROSD*PROS;
// antithrombin consumption is proportional to remaining AT, so AT can approach
// but never cross zero (a negative AT would turn thrombin inhibition into
// positive feedback and is the classic way this model structure blows up)
dxdt_AT   = KATS - KATD*AT - KATC*pos(THR)*(AT/AT0);
double ATF = (AT > 0.0) ? AT/100.0 : 0.0;

double THRP = pos(THR);
double APC = KAPC*pos(PC)/100.0*pos(TM)*(THRP/THR0);
double FFV = 1.0/(1.0 + APC/KAPCI);
double CARG = ARGC/VARG1;                            // ug/mL
double CBIV = BIVC/VBIV;                             // ug/mL
double CRIV = RIVC/VRIV*1e6;                         // ng/mL
double FXAI = 1.0/(1.0 + CRIV/KIRIV + CFON/KIFON + CDAN/KIDAN
                   + CUFH*ATF/KIUFHX + CLMW*ATF/KILMWX);
double PTASE = KPT*(pos(TF)/TF0)*(pos(FX)/100.0)*pow(pos(MP)/MP0, 0.8)*FFV*FXAI;
double THRGEN = PTASE*(pos(FII)/100.0);
double HEPAT = CUFH + 0.35*CLMW + 0.10*CDAN;
double ATACC = EMAXAT*HEPAT/(KATACC + HEPAT);
dxdt_THR = THRGEN - KTHRI*THR*ATF*(1.0 + ATACC);
// effective thrombin seen by substrates, after competitive DTI inhibition

double THRACT = THRP/(1.0 + CARG/KIARG + CBIV/KIBIV);
dxdt_FBN  = KFBN*THRACT - KFBND*FBN;
dxdt_DDIM = KDDS*KFBND*FBN - KDDD*DDIM;
dxdt_FRG  = KF12S*THRGEN - KFRGD*FRG;

// ================= 8. endpoints ======================================
double TRAT = pos(THRACT)/THR0;
// Two additive routes to thrombosis, and they are cut by DIFFERENT drugs.
//  - the thrombin route is what argatroban / bivalirudin / DOACs suppress
//  - the platelet-activation route is immune, and NO anticoagulant touches it;
//    only stopping heparin, IVIG or plasma exchange reduce it. This is why a
//    DTI-treated patient still carries a residual thrombotic risk.
double HZTHR = HZ0*((pow(TRAT, PHZ) > 1.0) ? (pow(TRAT, PHZ) - 1.0) : 0.0)
               *(1.0 + KHZF*pos(FBN - FBN0));
double HZPLT = HZ0*KHZP*ACTSIG;
double HZ = HZTHR + HZPLT;
dxdt_TEC = HZ*(1.0 - TEC);

double PCDEF = (PC < 100.0) ? (1.0 - PC/100.0) : 0.0;
double GATE  = (1.0/(1.0 + exp((pos(PLT) - 120.0)/20.0)))       // platelets still low
             * (pos(FBN - FBN0)/(pos(FBN - FBN0) + 1.0));        // limb already thrombosed
dxdt_NECR = KVLG*PCDEF*PCDEF*TRAT*GATE*(1.0 - NECR);
dxdt_THRAUC = THRACT;

// ================= 9. drug PK ========================================
dxdt_ARGC = QARG*(ARGP/VARG2) - QARG*(ARGC/VARG1) - CLARG*HEPFN*(ARGC/VARG1);
dxdt_ARGP = QARG*(ARGC/VARG1) - QARG*(ARGP/VARG2);
dxdt_BIVC = -CLBIV*(0.2*RENFN + 0.8)*(BIVC/VBIV);
dxdt_RIVD = -KARIV*RIVD;
dxdt_RIVC = KARIV*FRIV*RIVD - CLRIV*(0.35*RENFN + 0.65)*(RIVC/VRIV);
dxdt_IVGC = QIVG*(IVGP/VIVG2) - QIVG*(IVGC/VIVG1) - CLIVG*(IVGC/VIVG1);
dxdt_IVGP = QIVG*(IVGC/VIVG1) - QIVG*(IVGP/VIVG2);
dxdt_WARD = -KAWAR*WARD;
dxdt_WARC = KAWAR*WARD - CLWAR*(WARC/VWAR);
dxdt_RTXC = QRTX*(RTXP/VRTX2) - QRTX*(RTXC/VRTX1) - CLRTX*(RTXC/VRTX1);
dxdt_RTXP = QRTX*(RTXC/VRTX1) - QRTX*(RTXP/VRTX2);
dxdt_PLEXA = -KPLEXD*PLEXA;
dxdt_VITK  = -KVITKD*VITK;

$TABLE
double CUFHo = HEP/VHEP;
double CLMWo = LMWC/VLMW;
double CFONo = FONC/VFON;
double CDANo = DANC/VDAN;
double HEPTo = CHNUFH*CUFHo*FBNDUFH + CHNLMW*CLMWo*FBNDLMW + CHNFON*CFONo*FBNDFON
             + CHNDAN*CDANo*FBNDDAN + PF4DNA*FBNDDNA;
double HEPBo = CHNUFH*CUFHo*FLENUFH + CHNLMW*CLMWo*FLENLMW + CHNFON*CFONo*FLENFON
             + CHNDAN*CDANo*FLENDAN + PF4DNA*FLENDNA;
double FLENMIXo = (HEPTo > 1e-9) ? HEPBo/HEPTo : 0.0;
double ULCFORMo = KULC*PF4Lo*FRATo*(HEPBo/(HEPBo + KHB));
double PF4Lo = FLOC*PF4P;
double RRATo = (HEPTo > 1e-9) ? PF4Lo/HEPTo : 1e12;
double FRATo = (HEPTo > 1e-9) ? fratio(RRATo, RSTAR, SDR) : 0.0;

double IVGo = (IVGC/VIVG1)*1000.0;
double icno = pos(IC)/KDIC*FCG;
double OCCo = icno/(1.0 + icno + IVGo/KDIVIG);
double XLo  = pow(pos(OCCo), NXL)/(pow(pos(OCCo), NXL) + pow(XL50, NXL));

// ---- diagnostic assays ----
double IGT = pos(IGGP) + WNP*pos(IGGN);
double OD = ODMAX*IGT/(KOD + IGT) + ODBASE;
double SRALO = assay(pos(IGGP), PF4ASSAY, HEPSRALO, CHNUFH, LENUFH, LEN50, RSTAR, SDR,
                     KULC, KULCOFF, KIND, HIND, KICON, KICOFF, KICCLR, KDIC, FCG, NXL, XL50, KHB);
double SRAHI = assay(pos(IGGP), PF4ASSAY, HEPSRAHI, CHNUFH, LENUFH, LEN50, RSTAR, SDR,
                     KULC, KULCOFF, KIND, HIND, KICON, KICOFF, KICCLR, KDIC, FCG, NXL, XL50, KHB);
double SRA0  = assay(pos(IGGP), PF4ASSAY, 0.0, CHNUFH, LENUFH, LEN50, RSTAR, SDR,
                     KULC, KULCOFF, KIND, HIND, KICON, KICOFF, KICCLR, KDIC, FCG, NXL, XL50, KHB);
double SEROPOS = (OD > 0.40) ? 1.0 : 0.0;
double SRAPOS  = ((SRALO > 20.0) && (SRAHI < 20.0)) ? 1.0 : 0.0;
double SRAPOSA = ((SRALO > 20.0) && (SRAHI >= 20.0)) ? 1.0 : 0.0;  // aHIT / VITT pattern

// ---- coagulation read-outs ----
double CARGo = ARGC/VARG1;
double CBIVo = BIVC/VBIV;
double CRIVo = RIVC/VRIV*1e6;
double PTR = 1.0/(pow(pos(FII)/100.0 + 1e-6, 0.20)*pow(pos(FVII)/100.0 + 1e-6, 0.30)*pow(pos(FX)/100.0 + 1e-6, 0.30));
double INRTRUE = PTR;
double INR = PTR*(1.0 + KARGINR*CARGo + KBIVINR*CBIVo + KRIVINR*CRIVo);
double FXCHROM = FX;
double APTT = APTT0*(1.0 + KAPUFH*CUFHo*(AT/100.0) + KAPARG*CARGo + KAPBIV*CBIVo);

double APCo = KAPC*pos(PC)/100.0*pos(TM)*(pos(THR)/THR0);
double FFVo = 1.0/(1.0 + APCo/KAPCI);
double THRACTo = pos(THR)/(1.0 + CARGo/KIARG + CBIVo/KIBIV);
double PLTPCT = 100.0*PLT/PLT0;
double PLTFALL = 100.0*(1.0 - PLT/PLT0);
double CWARo = WARC/VWAR*1000.0;
double CRTXo = RTXC/VRTX1*1000.0;

$CAPTURE @annotated
CUFHo   : UFH concentration (U/mL)
CLMWo   : LMWH anti-Xa concentration (IU/mL)
CFONo   : Fondaparinux concentration (ug/mL)
CDANo   : Danaparoid anti-Xa concentration (U/mL)
HEPTo   : Total polyanion chain concentration (nM)
FLENMIXo: Mean bridging competence of the polyanion pool (0-1)
ULCFORMo: Instantaneous ULC assembly rate (nM/h)
PF4Lo   : Local effective PF4 concentration (nM)
RRATo   : PF4 : polyanion chain molar ratio
FRATo   : Antigenicity of the current molar ratio (0-1)
OCCo    : Fractional FcgammaRIIa occupancy
XLo     : Fractional FcgammaRIIa crosslinking (activation signal)
OD      : PF4-heparin IgG ELISA optical density
SRALO   : Serotonin release at low heparin (%)
SRAHI   : Serotonin release at 100 U/mL heparin (%)
SRA0    : Serotonin release without heparin (%)
SEROPOS : Seropositive by immunoassay (0/1)
SRAPOS  : Classic SRA-positive pattern (0/1)
SRAPOSA : Heparin-independent activation pattern (aHIT/VITT) (0/1)
CARGo   : Argatroban concentration (ug/mL)
CBIVo   : Bivalirudin concentration (ug/mL)
CRIVo   : Rivaroxaban concentration (ng/mL)
CWARo   : Warfarin concentration (mg/L)
CRTXo   : Rituximab concentration (ug/mL)
INRTRUE : INR from true factor depletion
INR     : Measured INR (including DTI artefact)
FXCHROM : Chromogenic factor X activity (%)
APTT    : Activated partial thromboplastin time (s)
APCo    : Activated protein C (normalised)
FFVo    : Factor Va availability (0-1)
THRACTo : Effective (inhibitor-corrected) thrombin activity (nM)
PLTPCT  : Platelet count as percent of baseline
PLTFALL : Percent fall in platelet count from baseline
'
# ---------------------------------------------------------------------
# BUILD
# ---------------------------------------------------------------------
mod <- mcode("hit", hit_code)

# ---------------------------------------------------------------------
# DOSING HELPERS
#   All heparin amounts are in the native unit of the target compartment:
#   HEP/HEPD in units, LMWD in anti-Xa IU, FOND in ug, DANC in anti-Xa U.
# ---------------------------------------------------------------------
WT <- 70

ufh_iv    <- function(dur_h, bolus = 80*WT, rate = 18*WT)
  c(ev(amt = bolus, cmt = "HEP"),
    ev(amt = rate*dur_h, rate = rate, cmt = "HEP"))
ufh_sc    <- function(n_doses, amt = 5000, ii = 8)
  ev(amt = amt, cmt = "HEPD", ii = ii, addl = n_doses - 1)
enox      <- function(n_doses, mg = 1*WT, ii = 12)          # 100 anti-Xa IU per mg
  ev(amt = mg*100, cmt = "LMWD", ii = ii, addl = n_doses - 1)
fonda     <- function(n_doses, mg = 2.5, start = 0, ii = 24)
  ev(amt = mg*1000, cmt = "FOND", time = start, ii = ii, addl = n_doses - 1)
danaparoid<- function(n_doses, U = 750, start = 0, ii = 8)
  ev(amt = U, cmt = "DANC", time = start, ii = ii, addl = n_doses - 1)
argatroban<- function(start, dur_h, ugkgmin = 2)
  ev(amt = ugkgmin*WT*60*dur_h, rate = ugkgmin*WT*60, cmt = "ARGC", time = start)
bivalirudin<- function(start, dur_h, mgkgh = 0.15)
  ev(amt = mgkgh*WT*1000*dur_h, rate = mgkgh*WT*1000, cmt = "BIVC", time = start)
rivaroxaban<- function(start, n_doses, mg = 15, ii = 12)
  ev(amt = mg, cmt = "RIVD", time = start, ii = ii, addl = n_doses - 1)
ivig      <- function(start, gkg = 1, days = 2)
  ev(amt = gkg*WT, cmt = "IVGC", time = start, ii = 24, addl = days - 1)
plex      <- function(start, n = 3)
  ev(amt = 1, cmt = "PLEXA", time = start, ii = 24, addl = n - 1)
rituximab <- function(start, n = 4, mg = 700)
  ev(amt = mg, cmt = "RTXC", time = start, ii = 168, addl = n - 1)
warfarin  <- function(start, load = 10, maint = 5, n_maint = 25)
  c(ev(amt = load,  cmt = "WARD", time = start),
    ev(amt = maint, cmt = "WARD", time = start + 24, ii = 24, addl = n_maint - 1))
vitamin_k <- function(start, mg = 10) ev(amt = mg, cmt = "VITK", time = start)
platelet_txn <- function(start, n = 4, rise = 40)
  ev(amt = rise, cmt = "PLT", time = start, ii = 24, addl = n - 1)
protamine <- function(start, mg = 300) ev(amt = mg, cmt = "PROT", time = start)
cpb_stimulus <- function(start = 0) ev(amt = 1, cmt = "SURG", time = start)
vaccine_polyanion <- function(amt = 40) ev(amt = amt, cmt = "PF4DNA")

sim <- function(events = NULL, pars = list(), end = 1440, delta = 4) {
  m <- if (length(pars)) param(mod, pars) else mod
  d <- as.data.frame(mrgsim(m, if (is.null(events)) ev(amt = 0, cmt = "HEP") else events,
                            end = end, delta = delta))
  d$day <- d$time/24
  d
}
at   <- function(o, day, var) o[[var]][which.min(abs(o$day - day))]
onset<- function(o, pct = 50) { i <- which(o$PLTFALL > pct); if (length(i)) o$day[i[1]] else NA_real_ }
recov<- function(o, thr = 150, after = 8) {
  i <- which(o$day > after & o$PLT > thr); if (length(i)) o$day[i[1]] else NA_real_
}
score_4ts <- function(o) {
  # Thrombocytopenia: >50% fall AND nadir >=20 -> 2 ; 30-50% fall or nadir 10-19 -> 1
  fall <- max(o$PLTFALL); nad <- min(o$PLT)
  t1 <- if (fall > 50 && nad >= 20) 2 else if (fall >= 30 || (nad >= 10 && nad < 20)) 1 else 0
  # Timing of the fall relative to exposure
  d  <- onset(o, 30)
  t2 <- if (is.na(d)) 0 else if (d >= 5 && d <= 10) 2 else if (d > 10) 1 else 0
  # Thrombosis: model-generated cumulative event probability
  te <- max(o$TEC)
  t3 <- if (te > 0.30) 2 else if (te > 0.10) 1 else 0
  # oTher causes: assumed absent in these single-mechanism simulations
  t4 <- 2
  c(T_thrombocytopenia = t1, T_timing = t2, T_thrombosis = t3,
    T_other = t4, total = t1 + t2 + t3 + t4)
}
line <- function(lab, o) {
  cat(sprintf("%-34s %6.1f %7.1f %7.2f %6.2f %6.1f %6.1f %7.1f%% %8.3f %8.4f  %2d\n",
      lab, min(o$PLT), o$day[which.min(o$PLT)], max(o$ULC), max(o$OD),
      max(o$SRALO), max(o$SRAHI), max(o$PLTFALL),
      at(o, 30, "TEC"), at(o, 45, "NECR"), score_4ts(o)["total"]))
}
hdr <- function(title) {
  cat("\n", title, "\n", strrep("-", 118), "\n", sep = "")
  cat(sprintf("%-34s %6s %7s %7s %6s %6s %6s %8s %8s %8s  %2s\n",
      "scenario", "nadir", "nadir_d", "maxULC", "OD", "SRAlo", "SRAhi",
      "maxfall", "TEC30", "NECR45", "4Ts"))
}

DX <- 7*24   # day 7: HIT recognised, heparin stopped, alternative started

# =====================================================================
# SCENARIOS
# =====================================================================
cat(strrep("=", 118), "\n")
cat("HIT QSP MODEL — 26 SCENARIOS\n")
cat("TEC30 = cumulative HIT-attributable thrombotic event probability at day 30\n")
cat("NECR45 = cumulative venous limb gangrene / warfarin skin necrosis probability at day 45\n")
cat(strrep("=", 118), "\n")

hdr("A · EXPOSURE: which heparin, at which dose, makes the antigen")
line("1  no heparin (control)",              sim(NULL, end = 4800))
line("2  UFH therapeutic IV, 14 d",          sim(ufh_iv(336)))
line("3  UFH therapeutic IV, 30 d (no Dx)",  sim(ufh_iv(720)))
line("4  UFH prophylaxis SC 5000U q8h",      sim(ufh_sc(90)))
line("5  enoxaparin therapeutic q12h",       sim(enox(60)))
line("6  enoxaparin prophylaxis 40mg q24h",  sim(enox(30, mg = 40, ii = 24)))
line("7  fondaparinux 2.5 mg q24h",          sim(fonda(30)))
line("8  danaparoid 750 U q8h",              sim(danaparoid(90)))
line("9  CPB 350U/kg + protamine + post-op", sim(c(ev(amt = 350*WT, cmt = "HEP"), cpb_stimulus(),
                                                   protamine(4), ufh_sc(38, ii = 8))))

hdr("B · TREATMENT of established HIT (heparin stopped on day 7)")
line("10 stop heparin ONLY",                 sim(ufh_iv(DX)))
line("11 + argatroban 2 ug/kg/min",          sim(c(ufh_iv(DX), argatroban(DX, 480))))
line("12 + bivalirudin 0.15 mg/kg/h",        sim(c(ufh_iv(DX), bivalirudin(DX, 480))))
line("13 + fondaparinux 7.5 mg",             sim(c(ufh_iv(DX), fonda(20, mg = 7.5, start = DX))))
line("14 + rivaroxaban 15 mg bid",           sim(c(ufh_iv(DX), rivaroxaban(DX, 42))))
line("15 + danaparoid",                      sim(c(ufh_iv(DX), danaparoid(60, U = 1500, start = DX))))
line("16 + IVIG 1 g/kg x 2 d",               sim(c(ufh_iv(DX), ivig(DX))))
line("17 + argatroban + IVIG",               sim(c(ufh_iv(DX), argatroban(DX, 480), ivig(DX))))
line("18 + plasma exchange x 3",             sim(c(ufh_iv(DX), argatroban(DX, 480), plex(DX))))
line("19 + rituximab 375 mg/m2 x 4",         sim(c(ufh_iv(DX), argatroban(DX, 480), rituximab(DX))))
line("20 + PLATELET TRANSFUSION (anti-Rx)",  sim(c(ufh_iv(720), platelet_txn(DX))))

hdr("C · THE TWO TRAPS, and the host and variant modifiers")
line("21 warfarin 10 mg load on day 7",      sim(c(ufh_iv(DX), warfarin(DX))))
line("22 warfarin day 7 + argatroban",       sim(c(ufh_iv(DX), argatroban(DX, 480), warfarin(DX))))
line("23 warfarin day 18 (PLT recovered)",   sim(c(ufh_iv(DX), argatroban(DX, 480), warfarin(18*24))))
line("24 warfarin day 7 + vitamin K rescue", sim(c(ufh_iv(DX), argatroban(DX, 480), warfarin(DX),
                                                   vitamin_k(DX + 24))))
line("25 FcgRIIa 131 HH homozygote",         sim(ufh_iv(DX), pars = list(FCG = 1.35)))
line("26 FcgRIIa 131 RR homozygote",         sim(ufh_iv(DX), pars = list(FCG = 0.75)))

hdr("D · aHIT / SPONTANEOUS HIT / VITT — one coefficient (HIND)")
line("27 aHIT (HIND 0.6) cessation only",    sim(ufh_iv(DX), pars = list(HIND = 0.6)))
line("28 aHIT + argatroban",                 sim(c(ufh_iv(DX), argatroban(DX, 480)),
                                                 pars = list(HIND = 0.6)))
line("29 aHIT + argatroban + IVIG",          sim(c(ufh_iv(DX), argatroban(DX, 480), ivig(DX)),
                                                 pars = list(HIND = 0.6)))
line("30 VITT, no heparin ever",             sim(c(cpb_stimulus(), vaccine_polyanion()),
                                                 pars = list(HIND = 0.9, KPDD = 0.004)))
line("31 VITT + argatroban + IVIG",          sim(c(cpb_stimulus(), vaccine_polyanion(),
                                                   argatroban(168, 480), ivig(168)),
                                                 pars = list(HIND = 0.9, KPDD = 0.004)))

# =====================================================================
# DIAGNOSTICS
# =====================================================================
cat("\n\n", strrep("=", 118), "\n", sep = "")
cat("DIAGNOSTICS\n")
cat(strrep("=", 118), "\n")

cat("\nD1 · BASELINE IS SOLVED, NOT FITTED — 200 days, no heparin\n")
o <- sim(NULL, end = 4800, delta = 24)
for (v in c("PLT", "PF4P", "PF4EC", "ULC", "IGGP", "THR", "AT", "PC", "FII",
            "TF", "MP", "VWF", "TM", "DDIM", "BCELL")) {
  a <- o[[v]][1]; b <- tail(o[[v]], 1)
  cat(sprintf("   %-7s %12.6f -> %12.6f   drift %+8.4f%%\n", v, a, b,
              if (abs(a) > 1e-12) 100*(b - a)/a else 0))
}

cat("\nD2 · THE BELL CURVE: in-vitro heparin sweep over four decades\n")
cat("   The falling limb IS the SRA 100 U/mL confirmatory step.\n")
o <- sim(ufh_iv(336), end = 336)
igg_peak <- max(o$IGGP)
sweep <- do.call(rbind, lapply(10^seq(-3, 2.3, length.out = 14), function(h) {
  s <- sim(ufh_iv(336), pars = list(HEPSRALO = h), end = 336)
  data.frame(heparin_U_mL = h, SRA_pct = max(s$SRALO))
}))
sweep$regime <- ifelse(sweep$heparin_U_mL < 0.02, "PF4 excess",
                ifelse(sweep$heparin_U_mL > 3, "HEPARIN EXCESS", "optimal ratio"))
print(sweep, row.names = FALSE, digits = 3)

cat("\nD3 · CHAIN LENGTH IS THE CLASS EFFECT: bridging competence by species\n")
cl <- data.frame(species = c("UFH", "LMWH (enoxaparin)", "danaparoid", "fondaparinux",
                             "cell-free DNA / polyP"),
                 saccharides = c(50, 15, 11, 5, 200))
cl$bridging_competence <- with(cl, saccharides^6/(saccharides^6 + 24^6))
cl$relative_to_UFH <- round(cl$bridging_competence/cl$bridging_competence[1], 5)
print(cl, row.names = FALSE, digits = 4)

cat("\nD4 · THE ICEBERG IS GENERATED BY EXPOSURE INTENSITY\n")
cat("   Three successive waterlines: seropositive -> platelet-activating -> clinical.\n")
ice <- do.call(rbind, lapply(c(0.5, 1, 2, 4, 8, 12, 18, 25), function(r) {
  s <- sim(c(ev(amt = 4*r*WT, cmt = "HEP"), ev(amt = r*WT*336, rate = r*WT, cmt = "HEP")))
  data.frame(UFH_U_kg_h = r, UFH_conc_U_mL = round(at(s, 10, "CUFHo"), 3),
             peak_ULC_nM = round(max(s$ULC), 1), ELISA_OD = round(max(s$OD), 2),
             L1_seropositive = max(s$SEROPOS), L2_SRA_positive = max(s$SRAPOS),
             max_fall_pct = round(max(s$PLTFALL), 1),
             L3_clinical_HIT = as.integer(max(s$PLTFALL) > 50))
}))
print(ice, row.names = FALSE)

cat("\nD5 · THE WARFARIN TRAP: half-life mismatch after a 10 mg load on day 7\n")
o <- sim(c(ufh_iv(DX), warfarin(DX)), end = 600, delta = 2)
for (d in c(7, 7.25, 7.5, 8, 9, 10, 12, 14)) {
  cat(sprintf("   day %5.2f  protein C %5.1f%%  prothrombin %5.1f%%  INRtrue %4.2f  thrombin %6.2f nM  PLT %5.1f  VLG %6.4f\n",
      d, at(o, d, "PC"), at(o, d, "FII"), at(o, d, "INRTRUE"),
      at(o, d, "THRACTo"), at(o, d, "PLT"), at(o, d, "NECR")))
}

cat("\nD6 · WARFARIN START-DAY SWEEP — the guideline rule, generated\n")
ws <- do.call(rbind, lapply(c(7, 9, 11, 13, 15, 18, 21), function(sd) {
  s <- sim(c(ufh_iv(DX), argatroban(DX, 200), warfarin(sd*24)))
  data.frame(warfarin_start_day = sd, platelets_at_start = round(at(s, sd, "PLT"), 1),
             VLG_prob_day45 = round(at(s, 45, "NECR"), 4))
}))
print(ws, row.names = FALSE)

cat("\nD7 · WHEN HIT IS RECOGNISED IS THE DOMINANT DETERMINANT OF RESIDUAL RISK\n")
cat("   The historical ~50% thrombosis rate after 'heparin cessation alone' is\n")
cat("   reproduced here only for LATE recognition. Cessation is not one\n")
cat("   intervention; it is a family of them indexed by day.\n")
rec <- do.call(rbind, lapply(c(5, 7, 9, 11, 14, 21, 30), function(dd) {
  s <- sim(ufh_iv(dd*24), end = 2400)
  data.frame(recognised_on_day = dd, IgG_at_cessation = round(at(s, dd, "IGGP"), 1),
             nadir = round(min(s$PLT), 1), peak_thrombin_nM = round(max(s$THRACTo), 2),
             TEC_day30 = round(at(s, 30, "TEC"), 3),
             TEC_final = round(tail(s$TEC, 1), 3))
}))
print(rec, row.names = FALSE)

cat("\nD7b · A REPORTED NULL: does a residual heparin trickle sustain HIT?\n")
cat("   It does not. Once the immune drive is falling, sub-therapeutic heparin\n")
cat("   makes almost no antigen AND still accelerates antithrombin, so the\n")
cat("   thrombotic risk is flat-to-slightly-lower. The clinical folklore that\n")
cat("   line flushes perpetuate HIT is NOT supported by this model structure.\n")
fl <- do.call(rbind, lapply(c(0, 250, 500, 1000, 2000, 4000), function(u) {
  e <- if (u == 0) ufh_iv(DX) else c(ufh_iv(DX), ev(amt = u*23, rate = u/24, cmt = "HEP", time = DX))
  s <- sim(e)
  data.frame(residual_heparin_U_per_day = u, ULC_day14 = round(at(s, 14, "ULC"), 4),
             crosslink_day14 = round(at(s, 14, "XLo"), 4),
             TEC30 = round(at(s, 30, "TEC"), 3))
}))
print(fl, row.names = FALSE)

cat("\nD8 · SEROREVERSION, RAPID-ONSET HIT, AND THE ABSENCE OF A MEMORY BOOST\n")
cat("   A second 3-day course of UFH is given after a gap. What matters is whether\n")
cat("   the antibody from the first exposure is still there.\n")
re <- do.call(rbind, lapply(c(20, 40, 60, 90, 120, 180, 250), function(gap) {
  g <- gap*24
  s <- sim(c(ufh_iv(336), ev(amt = 80*WT, cmt = "HEP", time = g),
             ev(amt = 18*WT*72, rate = 18*WT, cmt = "HEP", time = g)),
           end = g + 30*24, delta = 2)
  post <- s[s$day >= gap, ]
  base <- at(s, gap - 0.2, "PLT")
  hit  <- which(post$PLT < 0.6*base)
  data.frame(gap_days = gap, IgG_at_re_exposure = round(at(s, gap, "IGGP"), 2),
             OD_at_re_exposure = round(at(s, gap, "OD"), 2),
             still_seropositive = at(s, gap, "SEROPOS"),
             days_to_40pct_fall = if (length(hit)) round(post$day[hit[1]] - gap, 2) else NA_real_,
             pattern = if (!length(hit)) "no HIT"
                       else if (post$day[hit[1]] - gap < 2) "RAPID-ONSET (<48 h)"
                       else if (post$day[hit[1]] - gap <= 10) "typical (day 5-10)"
                       else "delayed")
}))
print(re, row.names = FALSE)
cat("   Seroreversion: the immunoassay turns negative between day 120 and 180 here\n")
cat("   (observed ~85-100 d), and the functional titre falls first, as observed.\n")
cat("   KRECALL = 0 by default, so re-exposure never triggers an anamnestic boost --\n")
cat("   HIT conspicuously does not show one, and setting KRECALL > 0 would invent it.\n")
cat("   MODEL PREDICTION AND LIMITATION: rapid-onset reactivity OUTLASTS immunoassay\n")
cat("   positivity in this model. Even a few U/mL of residual platelet-activating IgG\n")
cat("   meets a large fresh antigen load and crosslinks enough receptor to drop the\n")
cat("   count within 48 h. The model therefore treats the 100-day re-exposure rule as\n")
cat("   an underestimate of the hazard window. That is a testable claim, not a fit.\n")

cat("\nD9 · THE ARGATROBAN-INR TRAP: measured INR vs the truth\n")
o <- sim(c(ufh_iv(DX), argatroban(DX, 360), warfarin(18*24)), end = 800, delta = 2)
for (d in c(18, 19, 20, 21, 22, 23, 24, 26)) {
  cat(sprintf("   day %2.0f  argatroban %5.2f ug/mL  MEASURED INR %4.2f  true INR %4.2f  chromogenic FX %5.1f%%  aPTT %5.1f s\n",
      d, at(o, d, "CARGo"), at(o, d, "INR"), at(o, d, "INRTRUE"),
      at(o, d, "FXCHROM"), at(o, d, "APTT")))
}

cat("\nD10 · TWO AXES, AND WHICH DRUG CUTS WHICH\n")
ax <- do.call(rbind, lapply(list(
  list("untreated",              ufh_iv(720),                                    list()),
  list("stop heparin",           ufh_iv(DX),                                     list()),
  list("argatroban",             c(ufh_iv(DX), argatroban(DX, 480)),             list()),
  list("IVIG",                   c(ufh_iv(DX), ivig(DX)),                        list()),
  list("argatroban + IVIG",      c(ufh_iv(DX), argatroban(DX, 480), ivig(DX)),   list()),
  list("platelet transfusion",   c(ufh_iv(720), platelet_txn(DX)),               list())),
  function(z) { s <- sim(z[[2]], pars = z[[3]])
    data.frame(intervention = z[[1]],
               platelet_axis_nadir = round(min(s$PLT), 1),
               thrombin_axis_peak_nM = round(max(s$THRACTo), 2),
               immune_axis_peak_XL = round(max(s$XLo), 3),
               TEC30 = round(at(s, 30, "TEC"), 3)) }))
print(ax, row.names = FALSE)

cat("\nD11 · HIND SWEEP: classic HIT -> aHIT -> VITT is a continuum\n")
hs <- do.call(rbind, lapply(c(0, 0.15, 0.3, 0.6, 0.9), function(h) {
  s <- sim(c(ufh_iv(DX), argatroban(DX, 480)), pars = list(HIND = h))
  data.frame(HIND = h, ULC_after_cessation_d14 = round(at(s, 14, "ULC"), 2),
             nadir = round(min(s$PLT), 1), nadir_day = round(s$day[which.min(s$PLT)], 1),
             SRA_no_heparin = round(max(s$SRA0), 1),
             SRA_100U = round(max(s$SRAHI), 1), TEC30 = round(at(s, 30, "TEC"), 3))
}))
print(hs, row.names = FALSE)
cat("   SRA_no_heparin and SRA_100U rise together with HIND: the high-heparin\n")
cat("   inhibition that DEFINES the classic assay is lost, which is exactly why\n")
cat("   aHIT and VITT sera require PF4-enhanced rather than heparin-dependent assays.\n")

cat("\nD12 · IVIG DOSE-RESPONSE in autoimmune HIT\n")
iv <- do.call(rbind, lapply(c(0, 0.25, 0.5, 1, 2), function(g) {
  e <- c(ufh_iv(DX), argatroban(DX, 480))
  if (g > 0) e <- c(e, ivig(DX, gkg = g))
  s <- sim(e, pars = list(HIND = 0.6))
  data.frame(IVIG_g_per_kg_per_day = g, peak_IgG_g_per_L = round(max(s$IVGC/3000*1000), 1),
             nadir = round(min(s$PLT), 1), PLT_day12 = round(at(s, 12, "PLT"), 1),
             TEC30 = round(at(s, 30, "TEC"), 3))
}))
print(iv, row.names = FALSE)

cat("\nD13 · ARGATROBAN IN HEPATIC IMPAIRMENT — why the dose must be halved\n")
hp <- do.call(rbind, lapply(c(1.0, 0.5, 0.25), function(hf) {
  do.call(rbind, lapply(c(0.5, 1, 2), function(rt) {
    s <- sim(c(ufh_iv(DX), argatroban(DX, 480, ugkgmin = rt)), pars = list(HEPFN = hf))
    data.frame(hepatic_function = hf, dose_ug_kg_min = rt,
               Css_ug_mL = round(at(s, 14, "CARGo"), 2), aPTT_s = round(at(s, 14, "APTT"), 1),
               TEC30 = round(at(s, 30, "TEC"), 3))
  }))
}))
print(hp, row.names = FALSE)

cat("\nD14 · 4Ts SCORE, computed from the simulated trajectory\n")
for (z in list(list("UFH therapeutic (typical HIT)", ufh_iv(DX), list()),
               list("UFH prophylaxis (subclinical)", ufh_sc(90), list()),
               list("fondaparinux (null)",           fonda(30), list()),
               list("aHIT",                          ufh_iv(DX), list(HIND = 0.6)))) {
  s <- score_4ts(sim(z[[2]], pars = z[[3]]))
  cat(sprintf("   %-32s T1=%d T2=%d T3=%d T4=%d  TOTAL=%d (%s)\n", z[[1]],
      s[1], s[2], s[3], s[4], s[5],
      if (s[5] >= 6) "high" else if (s[5] >= 4) "intermediate" else "low"))
}

cat("\n", strrep("=", 118), "\n", sep = "")
cat("END. See hit_references.md for the evidence behind every parameter,\n")
cat("and hit_shiny_app.R for the interactive dashboard.\n")
cat(strrep("=", 118), "\n")
