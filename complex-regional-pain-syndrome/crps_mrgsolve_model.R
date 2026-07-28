## =============================================================================
## Complex Regional Pain Syndrome (CRPS) — mrgsolve QSP model
## 복합부위통증증후군 정량적 시스템 약리학 모델
##
## 34 ODEs (14 PK + 20 disease/endpoint) · 9 treatment scenarios ·
## 8 analysis functions.  Time unit = hours.
##
## -----------------------------------------------------------------------------
## WHAT THIS MODEL ASSERTS (and what it turned out to say instead)
## -----------------------------------------------------------------------------
## CRPS is normally presented as a *list* of mechanisms: neurogenic inflammation,
## autoantibodies, sympatho-afferent coupling, oxidative stress, tissue hypoxia,
## spinal sensitisation, glial activation, cortical reorganisation, disuse
## osteopenia. This model imposes a TOPOLOGY on that list, and the topology is
## the scientific content:
##
##   (A) A FAST FEED-FORWARD PERIPHERAL NODE — SP/CGRP, cytokines, NGF, ROS,
##       oedema, alpha1 upregulation, osteoclast activation. Time constants of
##       days; driven by an injury input that decays with tau = 400 h. Left
##       alone it returns toward baseline. It cannot, by itself, be chronic.
##
##   (B) TWO NESTED POSITIVE-FEEDBACK LOOPS with month-scale time constants:
##
##       B1  behavioural-cortical ring
##           PAIN -(KFEAR, Hill n=3)-> DISUSE -(W_DIS_CTX, Hill n=3)->
##           CORTEX -(W_CTX_PAIN)-> PAIN
##           (with a second branch DISUSE -> hypoxia/bone resorption -> PSENS)
##
##       B2  neuroimmune latch
##           SSENS -> GLIA -(W_GLIA_SELF, Hill n=16 about GLIA50)-> SSENS
##           i.e. glial priming: above a threshold of glial activation, spinal
##           sensitisation supplies its own drive and no longer needs afferent
##           input. KOUT_GLIA = 8e-4 /h gives the latch a ~5-week memory.
##
## Because (A) decays while (B1)/(B2) latch, the model has THREE stable states
## for one parameter set, and every clinical statement below is a statement
## about which state you land in and when you can still change it. All numbers
## quoted are outputs of this model (verified against an independent Python
## transcription of the same equations, see README section 5), not literature
## values; the literature enters only through the calibration anchors listed
## further down.
##
##   attractor 0  resolution                     NRS 0.0, CSS 0,    BMD 1.00
##   attractor 1  behavioural-cortical chronic    NRS 5.5, CSS 4.6,  GLIA 0.23
##   attractor 2  neuroimmune-latched chronic     NRS 6.9, CSS 6.8,  GLIA 0.59
##
## FIVE COMPUTED CONSEQUENCES (functions, not prose):
##
##  1. CRPS_trait_bifurcation() — the fear-avoidance gain KFEAR, a *psychological*
##     parameter, is the switch. KFEAR 0.50 -> full resolution; 0.60 -> attractor
##     1 (NRS 4.4); 0.90 -> attractor 2 (NRS 6.9). Two successive bifurcations
##     in a trait, with the injury unchanged.
##
##  2. CRPS_insult_scan() — the insult magnitude does NOT decide chronicity. Over
##     INJ_AMP 0.10-0.90 the 3-year NRS is 5.50 in every case; at 0.95 it steps
##     to 6.92. Injury size selects the SEVERITY TIER; the trait selects whether
##     there is a tier at all.
##
##  3. CRPS_window_scan() — the identical multimodal package (prednisolone taper
##     + NAC + rehabilitation) gives complete resolution when started on day
##     3-90 and gives *nothing* from day 95 onward (NRS 6.92 = untreated).
##     The critical delay t* lies between 90 and 95 days and is set by when
##     GLIA crosses GLIA50, not by any drug property.
##
##  4. CRPS_arm_decomposition() — the window belongs to the rehabilitation arm.
##     Rehabilitation alone resolves at day 7/30/60. Prednisolone alone or NAC
##     alone NEVER resolve; their entire measurable benefit is preventing the
##     late escalation from attractor 1 to attractor 2 (5.50 vs 6.92) and only
##     if given within ~30 days. Consequence for trial design: an
##     anti-inflammatory arm is invisible both in patients who were never going
##     to escalate and in patients who already have.
##
##  5. CRPS_dose_vs_timing() — AGAINST THE HYPOTHESIS THIS MODEL WAS WRITTEN TO
##     EXPRESS. The model was built expecting ketamine + rehabilitation to be
##     supra-additive ("the analgesic window lets the rehabilitation work").
##     It is not. At the 700-day endpoint the interaction term is NEGATIVE
##     (-1.42 NRS): rehabilitation alone already reaches attractor 0, so the
##     drug can only be redundant at the endpoint, even though it clearly helps
##     early (day 34: 1.89 with the combination vs 4.05 for rehabilitation
##     alone). And on the dose axis, going from 22 to 88 mg/h — a 4-fold
##     increase, Cmax 258 -> 1034 ng/mL — changes the 2-year NRS by exactly
##     0.00, while moving the start day from 60 to 120 changes it by 1.42.
##     Ketamine's whole long-run contribution is a tier downgrade, and the
##     dose axis is flat above threshold. A stronger analgesic is not a
##     stronger disease-modifier.
##
## -----------------------------------------------------------------------------
## CALIBRATION ANCHORS (PMIDs in crps_references.md)
## -----------------------------------------------------------------------------
##  Epidemiology / course
##   - Incidence 26.2/100,000 py: de Mos 2007 Pain; 5.46/100,000: Sandroni 2003
##   - ~7% after distal radius fracture: Beerthuizen 2012 Pain
##   - Budapest criteria: Harden 2010 Pain; CRPS Severity Score: Harden 2017
##   - warm(acute) -> cold(chronic) phenotype shift: Bruehl 2016; Eberle 2009
##  Peripheral / immune
##   - blister-fluid IL-6, TNF-alpha up: Huygen 2002; SP/CGRP: Birklein 2001
##   - IENFD reduced ~29%: Oaklander 2006
##   - autoantibodies to beta2-AR/M2/AT1R: Kohr 2011; passive transfer: Tekus 2014
##   - IVIG positive small trial Goebel 2010; LIPS RCT NEGATIVE Goebel 2017
##     -> W_AAB_PS kept small and IVIG effect is nil late (model reproduces this)
##   - DMSO/NAC by warm/cold subtype: Perez 2003; vitamin C prophylaxis:
##     Zollinger 2007; ET-1 up / NO down: Groeneweg 2006; tadalafil: Groeneweg 2008
##  Bone
##   - neridronate 100 mg IV x4, VAS benefit persisting to 1 year: Varenna 2013,
##     2017 -> model gives -1.83 NRS at 1 year with CTX-I -61% and BMD preserved
##   - pamidronate: Robinson 2004; alendronate: Manicourt 2004
##  Central / NMDA
##   - S-ketamine 100-h infusion, NRS 7.2 -> 2.7 during infusion, benefit gone
##     by ~11 weeks: Sigtermans 2009; PK-PD Dahan 2011
##     -> model gives 6.04 -> 3.02 (day 34) and return to within 0.3 of
##        untreated by day 44
##   - glial priming / self-sustaining central sensitisation: Ji 2013, 2018
##  Cortical / rehabilitation
##   - graded motor imagery RCTs: Moseley 2004, 2006; S1 map normalises with
##     pain relief: Maihofner 2004; mirror therapy: McCabe 2003
##   - pain-exposure physiotherapy: Barnhoorn 2015
##   - pain-related fear predicts disability: de Jong 2005; Bean 2014
##  Neuromodulation
##   - SCS positive at 6 months: Kemler 2000 NEJM; between-group difference lost
##     at 5 years: Kemler 2008 -> SCS_HAB_K = 8e-5 /h gives +3.27 NRS at month 7
##     decaying to +0.15 at month 60
##   - DRG stimulation: Deer 2017 (ACCURATE)
##  Steroid
##   - prednisolone in early CRPS: Christensen 1982; Kalita 2006
## =============================================================================

$PROB
# CRPS QSP model — fast peripheral node + behavioural-cortical ring +
# neuroimmune (glial) latch + bone axis + PK/PD of 8 agents and 3 devices.

$PARAM
// ======================= KETAMINE PK (IV, 2-compartment) ====================
KET_V1     = 22      // L, central volume
KET_V2     = 200     // L, peripheral volume
KET_Q      = 60      // L/h, intercompartmental clearance
KET_CL     = 85      // L/h, clearance (high hepatic extraction)
KET_FM     = 0.80    // fraction of clearance forming norketamine
NORKET_V   = 70      // L
NORKET_CL  = 40      // L/h
NORKET_POT = 0.30    // norketamine potency relative to ketamine at NMDA

// =========================== PREDNISOLONE PK (oral) ========================
PRED_KA    = 2.0     // 1/h
PRED_F     = 0.80
PRED_V     = 40      // L
PRED_CL    = 9.0     // L/h (t1/2 ~3 h)

// ==================== NERIDRONATE PK (IV, bone-seeking) ====================
NER_V      = 20      // L, plasma volume
NER_CLR    = 6.0     // L/h, renal clearance
NER_KBONE  = 4.0     // L/h, plasma-to-bone uptake (~40% of dose to bone)
NER_KOFF   = 2e-5    // 1/h, release from bone (months-years residence)

// ============================ GABAPENTIN PK (oral) =========================
GBP_VMAX   = 45      // mg/h, saturable LAT1 absorption
GBP_KM     = 420     // mg in gut
GBP_F      = 0.90
GBP_V      = 58      // L
GBP_CL     = 7.5     // L/h

// ========================== AMITRIPTYLINE PK (oral) ========================
AMT_KA     = 1.0     // 1/h
AMT_F      = 0.50
AMT_V      = 900     // L (~12 L/kg)
AMT_CL     = 45      // L/h

// ======================= N-ACETYLCYSTEINE PK (oral) ========================
NAC_KA     = 1.2     // 1/h
NAC_F      = 0.10    // extensive first pass
NAC_V      = 35      // L
NAC_CL     = 15      // L/h

// =============================== IVIG PK ===================================
IVIG_V     = 5.0     // L
IVIG_KEL   = 0.0014  // 1/h (t1/2 ~3 weeks)

// ==================== INCITING EVENT / IMMOBILISATION ======================
INJ_AMP    = 1.0     // injury drive amplitude (1 = distal radius fracture)
INJ_TAU    = 400     // h, injury drive decay time constant (~17 d)
INJ_T0     = 0       // h, time of the inciting event
IMMOB_DUR  = 720     // h, cast immobilisation duration (30 d)
IMMOB_AMP  = 0.25    // disuse drive contributed by the cast

// ============ FAST PERIPHERAL NODE (kept in its linear regime) =============
KIN_NP     = 0.020   // 1/h, SP/CGRP release
KOUT_NP    = 0.050   // 1/h
NP_SENS_FB = 0.30    // afferent positive feedback onto neuropeptide release
KIN_CYT    = 0.020   // 1/h, IL-1b/IL-6/TNF composite
KOUT_CYT   = 0.045   // 1/h
W_NP_CYT   = 1.20    // neuropeptide -> keratinocyte/mast-cell cytokine drive
KIN_NGF    = 0.010   // 1/h
KOUT_NGF   = 0.015   // 1/h
KIN_EDEMA  = 0.020   // 1/h, plasma extravasation
KOUT_EDEMA = 0.030   // 1/h
KIN_ROS    = 0.050   // 1/h
KOUT_ROS   = 0.090   // 1/h
W_HYP_ROS  = 0.70    // hypoxia -> ROS
KIN_AAB    = 0.0015  // 1/h, autoantibody generation (slow)
KOUT_AAB   = 0.0018  // 1/h
W_AAB_PS   = 0.30    // autoantibody -> nociceptor gain (kept small: LIPS RCT)
KIN_A1     = 0.0025  // 1/h, alpha1 upregulation / sympathetic sprouting (slow)
KOUT_A1    = 0.006   // 1/h
W_NGF_A1   = 0.80    // NGF -> alpha1 / sprouting
SYMP_TONE  = 1.0     // patient covariate: sympathetic outflow scaling

// ================= VASOMOTOR / HYPOXIA (warm -> cold shift) ================
PERF_BASE  = 1.0     // symmetric perfusion index
KPERF      = 0.03    // 1/h
W_NP_PERF  = 3.0     // neurogenic vasodilation (WARM phenotype, acute)
W_A1_PERF  = 1.4     // alpha1 vasoconstriction (COLD phenotype)
W_ROS_PERF = 0.40    // endothelial dysfunction (ET-1 up / NO down)
W_SS_PERF  = 0.50    // central sympathetic vasoconstriction from spinal state
KIN_HYP    = 0.05    // 1/h
KOUT_HYP   = 0.04    // 1/h
W_DIS_HYP  = 0.20    // deconditioning -> nutritive-flow deficit
TEMP_GAIN  = 3.0     // degC temperature asymmetry per unit perfusion deviation

// ======================== PERIPHERAL SENSITISATION =========================
KIN_PS     = 0.040   // 1/h
KOUT_PS    = 0.070   // 1/h
W_NP_PS    = 0.35    // SP/CGRP + NGF (TRPV1, Nav1.7/1.8)
W_CYT_PS   = 0.40    // cytokines
W_ROS_PS   = 0.30    // ROS / lipid peroxidation
W_A1_PS    = 0.45    // sympatho-afferent coupling
W_HYP_PS   = 0.35    // acidosis / ASIC3
W_BONE_PS  = 0.30    // intraosseous acidosis from osteoclastic resorption

// ==================== SLOW CENTRAL RING AND GLIAL LATCH ====================
KIN_SS      = 0.012  // 1/h, spinal sensitisation build-up
KOUT_SS     = 0.010  // 1/h
W_PS_SS     = 1.60   // afferent barrage weight
KIN_GLIA    = 0.0016 // 1/h, microglial/astroglial activation
KOUT_GLIA   = 0.0008 // 1/h -> ~5-week glial memory
W_GLIA_SS   = 0.40   // glial amplification of afferent gain
W_GLIA_SELF = 1.30   // GLIAL SELF-DRIVE (the latch)
GLIA50      = 0.36   // latch threshold
HILL_GLIA   = 16     // latch steepness (glial priming is switch-like)
DINH_BASE   = 1.0    // descending inhibitory tone at baseline
KDINH       = 0.004  // 1/h
W_SS_DINH   = 0.55   // CPM/DNIC loss caused by spinal sensitisation
KIN_CTX     = 0.0025 // 1/h, cortical map degradation
KOUT_CTX    = 0.0015 // 1/h, spontaneous cortical recovery (slow)
W_SS_CTX    = 0.45   // nociceptive drive on cortical degradation
W_DIS_CTX   = 0.90   // non-use drive on cortical degradation
DIS50_CTX   = 0.45   // disuse threshold for use-dependent map loss
HILL_CTX    = 3      // steepness
KIN_DIS     = 0.006  // 1/h, disuse accumulation
KOUT_DIS    = 0.004  // 1/h, spontaneous return to use
KFEAR       = 1.00   // PAIN -> kinesiophobia -> DISUSE gain (TRAIT; the switch)
PAIN50_FEAR = 5.0    // NRS at which avoidance behaviour is half-maximal
HILL_FEAR   = 3      // steepness of the fear-avoidance response
W_ROM_DIS   = 0.15   // mechanical (contracture) contribution to disuse

// ====================== PAIN AND CLINICAL TRANSLATION ======================
PAIN_MAX     = 10.0  // NRS ceiling
KPAIN        = 0.25  // 1/h, pain state equilibration (~4 h)
PAIN50_DRIVE = 2.5   // drive giving half-maximal NRS
W_PS_PAIN    = 2.2   // peripheral sensitisation weight
W_SS_PAIN    = 3.6   // spinal sensitisation weight
W_CTX_PAIN   = 3.2   // cortical / body-perception weight (RING GAIN)
W_DIS_PAIN   = 1.6   // stiffness / deconditioning weight
KROM         = 0.010 // 1/h
ROM_MAX      = 1.0
W_EDEMA_ROM  = 0.35
W_PAIN_ROM   = 0.055 // per NRS point
W_DIS_ROM    = 0.45

// =============================== BONE AXIS =================================
KIN_OC    = 0.008    // 1/h, osteoclast activity adjustment
OC_BASE   = 1.0
W_NP_OC   = 0.70     // CGRP/SP -> RANKL/OPG
W_CYT_OC  = 0.50     // TNF/IL-6 -> osteoclastogenesis
W_DIS_OC  = 0.35     // loss of mechanical strain -> resorption
KBMD      = 0.00016  // 1/h, BMD turnover
BMD_REC   = 2.5      // loading-dependent BMD recovery multiplier
BMD_BASE  = 1.0      // regional BMD as fraction of contralateral limb
KCTX      = 0.05     // 1/h, serum CTX-I equilibration
CTX_BASE  = 0.35     // ng/mL
CTX_GAIN  = 0.30     // ng/mL per unit osteoclast activity above baseline

// ============================== DRUG PD ====================================
EC50_KET   = 130     // ng/mL, NMDA open-channel block
EMAX_KET   = 0.85    // maximal block of NMDA-dependent spinal gain
KET_GLIA_E = 0.30    // maximal ketamine effect on glial activation rate
EC50_PRED  = 25      // ng/mL
EMAX_PRED  = 0.75    // maximal NF-kB-driven cytokine suppression
PRED_EDEMA = 0.55    // maximal suppression of plasma extravasation
EMAX_NER   = 0.80    // maximal osteoclast inhibition (FPPS block)
BONE50_NER = 25      // mg bone-bound for half-maximal effect
EC50_GBP   = 4000    // ng/mL
EMAX_GBP   = 0.35    // alpha2delta-1 mediated reduction of glutamate release
EC50_AMT   = 60      // ng/mL
EMAX_AMT   = 0.60    // maximal restoration of descending inhibition
EC50_NAC   = 900     // ng/mL (also stands for topical DMSO 50%)
EMAX_NAC   = 0.55    // maximal ROS scavenging
EC50_IVIG  = 6.0     // g/L
EMAX_IVIG  = 0.65    // maximal neutralisation of autoantibody drive
VASODIL    = 0       // 0/1 PDE5-inhibitor-like vasodilator (cold phenotype)
VASODIL_E  = 0.45    // fractional reversal of vasoconstrictor drive

// ---- acute, state-independent analgesia (transmission block, not plasticity)
A_KET      = 0.75    // ketamine
A_GBP      = 0.30    // gabapentinoid
A_AMT      = 0.25    // TCA
A_SCS      = 0.80    // spinal cord stimulation (segmental gating)
A_MAX      = 0.85    // ceiling on combined acute analgesia

// =================== NON-PHARMACOLOGICAL INTERVENTIONS =====================
REHAB      = 0       // rehabilitation intensity 0-1 (PT + GMI + graded exposure)
REHAB_T0   = 0       // h, start
REHAB_DUR  = 4380    // h, duration (default 6 months)
REHAB_CTX  = 0.85    // use-dependent cortical remapping gain
REHAB_DIS  = 0.90    // graded-exposure reduction of the disuse state
REHAB_FEAR = 0.60    // CBT reduction of the PAIN->DISUSE gain itself
SCS_ON     = 0       // spinal cord stimulation switch
SCS_T0     = 1e9     // h, implantation time
SCS_EFF    = 0.55    // initial suppression of spinal sensitisation gain
SCS_HAB_K  = 8e-5    // 1/h, habituation (Kemler 5-year attenuation)
SYMPBLOCK  = 0       // sympathetic block switch
SB_T0      = 1e9     // h
SB_EFF     = 0.60    // fractional alpha1-drive suppression
SB_TAU     = 336     // h, block effect duration (2 weeks)

$CMT
// ---- PK (14) ----
KET_C1 KET_C2 NORKET
PRED_GUT PRED_C
NER_C NER_BONE
GBP_GUT GBP_C
AMT_GUT AMT_C
NAC_GUT NAC_C
IVIG_C
// ---- fast peripheral node (8) ----
NP CYT NGF EDEMA ROS AAB ALPHA1 PSENS
// ---- vasomotor / hypoxia (2) ----
PERF HYPOX
// ---- slow central ring + glial latch (5) ----
SSENS GLIA DINH CORTEX DISUSE
// ---- endpoints and bone (5) ----
PAIN ROM OC BMD CTXI

$MAIN
// healthy limb before the inciting event
PERF_0 = PERF_BASE;
DINH_0 = DINH_BASE;
ROM_0  = ROM_MAX;
OC_0   = OC_BASE;
BMD_0  = BMD_BASE;
CTXI_0 = CTX_BASE;

$ODE
// ==========================================================================
// 0. Concentrations and exposure-driven effects
// ==========================================================================
double ket_c     = 1000.0 * KET_C1 / KET_V1;          // ng/mL
double nork_c    = 1000.0 * NORKET / NORKET_V;        // ng/mL
double ket_eff_c = ket_c + NORKET_POT * nork_c;       // NMDA-equivalent
double pred_c    = 1000.0 * PRED_C / PRED_V;          // ng/mL
double gbp_c     = 1000.0 * GBP_C  / GBP_V;           // ng/mL
double amt_c     = 1000.0 * AMT_C  / AMT_V;           // ng/mL
double nac_c     = 1000.0 * NAC_C  / NAC_V;           // ng/mL
double ivig_c    = IVIG_C / IVIG_V;                   // g/L

double E_ket   = EMAX_KET   * ket_eff_c / (EC50_KET  + ket_eff_c);
double E_ketgl = KET_GLIA_E * ket_eff_c / (EC50_KET  + ket_eff_c);
double E_pred  = EMAX_PRED  * pred_c    / (EC50_PRED + pred_c);
double E_predE = PRED_EDEMA * pred_c    / (EC50_PRED + pred_c);
double E_gbp   = EMAX_GBP   * gbp_c     / (EC50_GBP  + gbp_c);
double E_amt   = EMAX_AMT   * amt_c     / (EC50_AMT  + amt_c);
double E_nac   = EMAX_NAC   * nac_c     / (EC50_NAC  + nac_c);
double E_ivig  = EMAX_IVIG  * ivig_c    / (EC50_IVIG + ivig_c);
double E_ner   = EMAX_NER   * NER_BONE  / (BONE50_NER + NER_BONE);

// ==========================================================================
// 1. Exogenous forcing: injury, cast, rehabilitation, devices
// ==========================================================================
double t_inj = SOLVERTIME - INJ_T0;
double INJ   = (t_inj >= 0.0) ? INJ_AMP * exp(-t_inj / INJ_TAU) : 0.0;
double IMMOB = (SOLVERTIME >= INJ_T0 && SOLVERTIME < INJ_T0 + IMMOB_DUR)
                 ? IMMOB_AMP : 0.0;
double rehab = (REHAB > 0.0 && SOLVERTIME >= REHAB_T0
                && SOLVERTIME < REHAB_T0 + REHAB_DUR) ? REHAB : 0.0;

double scs_t   = SOLVERTIME - SCS_T0;
double scs_eff = (SCS_ON > 0.0 && scs_t >= 0.0)
                   ? SCS_EFF * exp(-SCS_HAB_K * scs_t) : 0.0;
double sb_t    = SOLVERTIME - SB_T0;
double sb_eff  = (SYMPBLOCK > 0.0 && sb_t >= 0.0 && sb_t < SB_TAU)
                   ? SB_EFF * exp(-sb_t / (SB_TAU / 3.0)) : 0.0;

// ==========================================================================
// 2. PK
// ==========================================================================
dxdt_KET_C1 = KET_Q * (KET_C2 / KET_V2) - (KET_Q + KET_CL) * (KET_C1 / KET_V1);
dxdt_KET_C2 = KET_Q * (KET_C1 / KET_V1 - KET_C2 / KET_V2);
dxdt_NORKET = KET_FM * KET_CL * (KET_C1 / KET_V1) - NORKET_CL * (NORKET / NORKET_V);

dxdt_PRED_GUT = -PRED_KA * PRED_GUT;
dxdt_PRED_C   =  PRED_F * PRED_KA * PRED_GUT - PRED_CL * (PRED_C / PRED_V);

dxdt_NER_C    = -(NER_CLR + NER_KBONE) * (NER_C / NER_V) + NER_KOFF * NER_BONE;
dxdt_NER_BONE =  NER_KBONE * (NER_C / NER_V) - NER_KOFF * NER_BONE;

double gbp_abs = GBP_VMAX * GBP_GUT / (GBP_KM + GBP_GUT);   // saturable LAT1
dxdt_GBP_GUT = -gbp_abs;
dxdt_GBP_C   =  GBP_F * gbp_abs - GBP_CL * (GBP_C / GBP_V);

dxdt_AMT_GUT = -AMT_KA * AMT_GUT;
dxdt_AMT_C   =  AMT_F * AMT_KA * AMT_GUT - AMT_CL * (AMT_C / AMT_V);

dxdt_NAC_GUT = -NAC_KA * NAC_GUT;
dxdt_NAC_C   =  NAC_F * NAC_KA * NAC_GUT - NAC_CL * (NAC_C / NAC_V);

dxdt_IVIG_C  = -IVIG_KEL * IVIG_C;

// ==========================================================================
// 3. FAST PERIPHERAL NODE — feed-forward, days, all indices bounded 0-1
// ==========================================================================
dxdt_NP  = KIN_NP * (INJ + NP_SENS_FB * PSENS) * (1.0 - NP) - KOUT_NP * NP;

dxdt_CYT = KIN_CYT * (0.6 * INJ + W_NP_CYT * NP) * (1.0 - E_pred) * (1.0 - CYT)
           - KOUT_CYT * CYT;

dxdt_NGF = KIN_NGF * (0.5 * INJ + 0.8 * CYT) * (1.0 - NGF) - KOUT_NGF * NGF;

dxdt_EDEMA = KIN_EDEMA * (NP + 0.5 * CYT) * (1.0 - E_predE) * (1.0 - EDEMA)
             - KOUT_EDEMA * EDEMA * (1.0 + 0.5 * rehab);

dxdt_ROS = KIN_ROS * (CYT + W_HYP_ROS * HYPOX) * (1.0 - ROS)
           - KOUT_ROS * ROS * (1.0 + 3.0 * E_nac);

dxdt_AAB = KIN_AAB * (INJ + 0.4 * CYT) * (1.0 - AAB) - KOUT_AAB * AAB;

dxdt_ALPHA1 = KIN_A1 * (W_NGF_A1 * NGF + 0.3 * INJ) * SYMP_TONE * (1.0 - sb_eff)
              * (1.0 - ALPHA1) - KOUT_A1 * ALPHA1;

// Bone resorption above baseline acidifies bone and excites intraosseous
// nociceptors: this is the route by which a bisphosphonate becomes analgesic.
double bone_excess = (OC - OC_BASE > 0.0) ? (OC - OC_BASE) : 0.0;
double ps_in = W_NP_PS  * (NP + 0.6 * NGF)
             + W_CYT_PS * CYT
             + W_ROS_PS * ROS
             + W_A1_PS  * ALPHA1 * SYMP_TONE * (1.0 - sb_eff)
             + W_HYP_PS * HYPOX
             + W_AAB_PS * AAB * (1.0 - E_ivig)
             + W_BONE_PS * bone_excess;
dxdt_PSENS = KIN_PS * ps_in * (1.0 - PSENS) - KOUT_PS * PSENS;

// ==========================================================================
// 4. VASOMOTOR / HYPOXIA — acute warm limb, chronic cold limb
// ==========================================================================
double vaso_rev = (VASODIL > 0.0) ? VASODIL_E : 0.0;
double perf_target = PERF_BASE + W_NP_PERF * NP
                   - (W_A1_PERF * ALPHA1 * SYMP_TONE
                      + W_ROS_PERF * ROS
                      + W_SS_PERF * SSENS * SYMP_TONE) * (1.0 - vaso_rev);
dxdt_PERF = KPERF * (perf_target - PERF);

double perf_def = (PERF_BASE - PERF > 0.0) ? (PERF_BASE - PERF) : 0.0;
dxdt_HYPOX = KIN_HYP * (perf_def + W_DIS_HYP * DISUSE) * (1.0 - HYPOX)
             - KOUT_HYP * HYPOX;

// ==========================================================================
// 5. SLOW CENTRAL RING AND GLIAL LATCH
// ==========================================================================
// Glial self-drive: switch-like above GLIA50 (glial priming). This term is
// what makes spinal sensitisation able to outlive its afferent input.
double glia_self = W_GLIA_SELF * pow(GLIA, HILL_GLIA)
                   / (pow(GLIA50, HILL_GLIA) + pow(GLIA, HILL_GLIA));
double spinal_gain = (W_PS_SS * PSENS * (1.0 + W_GLIA_SS * GLIA) + glia_self)
                     * (1.0 - E_ket) * (1.0 - E_gbp) * (1.0 - scs_eff);
dxdt_SSENS = KIN_SS * spinal_gain * (1.0 - SSENS) - KOUT_SS * SSENS * DINH;

dxdt_GLIA = KIN_GLIA * SSENS * (1.0 - E_ketgl) * (1.0 - GLIA) - KOUT_GLIA * GLIA;

double dinh_target = DINH_BASE * (1.0 - W_SS_DINH * SSENS / (1.0 + SSENS))
                     * (1.0 + E_amt);
dxdt_DINH = KDINH * (dinh_target - DINH);

// Cortical map degradation: nociceptive drive + non-use, reversed only by
// use-dependent input (graded motor imagery / mirror therapy / loading).
double hill_dis = pow(DISUSE, HILL_CTX)
                  / (pow(DIS50_CTX, HILL_CTX) + pow(DISUSE, HILL_CTX));
dxdt_CORTEX = KIN_CTX * (W_SS_CTX * SSENS + W_DIS_CTX * hill_dis)
              * (1.0 - CORTEX)
              - (KOUT_CTX + REHAB_CTX * KIN_CTX * rehab * 8.0) * CORTEX;

// Disuse / fear-avoidance. CBT and graded exposure reduce the GAIN (KFEAR),
// not merely the state -- which is why they can move the attractor.
double fear_gain = KFEAR * (1.0 - REHAB_FEAR * rehab);
double hill_pain = pow(PAIN, HILL_FEAR)
                   / (pow(PAIN50_FEAR, HILL_FEAR) + pow(PAIN, HILL_FEAR));
double dis_in = fear_gain * hill_pain + W_ROM_DIS * (ROM_MAX - ROM) + IMMOB;
dxdt_DISUSE = KIN_DIS * dis_in * (1.0 - DISUSE)
              - (KOUT_DIS + REHAB_DIS * KIN_DIS * rehab * 6.0) * DISUSE;

// ==========================================================================
// 6. CLINICAL TRANSLATION
// ==========================================================================
double pain_drive = W_PS_PAIN * PSENS + W_SS_PAIN * SSENS
                  + W_CTX_PAIN * CORTEX + W_DIS_PAIN * DISUSE;
double acute = A_KET * E_ket + A_GBP * E_gbp + A_AMT * E_amt + A_SCS * scs_eff;
if(acute > A_MAX) acute = A_MAX;
double pain_target = PAIN_MAX * (pain_drive / (PAIN50_DRIVE + pain_drive))
                     * (1.0 - acute);
dxdt_PAIN = KPAIN * (pain_target - PAIN);

double rom_target = ROM_MAX * (1.0 - W_EDEMA_ROM * EDEMA
                                  - W_PAIN_ROM * PAIN
                                  - W_DIS_ROM * DISUSE);
if(rom_target < 0.05) rom_target = 0.05;
dxdt_ROM = KROM * (rom_target - ROM) * (1.0 + 2.0 * rehab);

// ==========================================================================
// 7. BONE AXIS
// ==========================================================================
double oc_target = OC_BASE * (1.0 + W_NP_OC * NP + W_CYT_OC * CYT
                              + W_DIS_OC * DISUSE) * (1.0 - E_ner);
dxdt_OC   = KIN_OC * (oc_target - OC);
dxdt_CTXI = KCTX * (CTX_BASE + CTX_GAIN * (OC - OC_BASE) - CTXI);
// BMD is expressed as a fraction of the contralateral limb, so only excess
// resorption can lose it and loading can only return it toward 1.
dxdt_BMD  = -KBMD * bone_excess * BMD
            + KBMD * BMD_REC * (1.0 - DISUSE) * (BMD_BASE - BMD);

$TABLE
double KET_ng  = 1000.0 * KET_C1 / KET_V1;
double NORK_ng = 1000.0 * NORKET / NORKET_V;
double PRED_ng = 1000.0 * PRED_C / PRED_V;
double GBP_ng  = 1000.0 * GBP_C / GBP_V;
double AMT_ng  = 1000.0 * AMT_C / AMT_V;
double NAC_ng  = 1000.0 * NAC_C / NAC_V;
double IVIG_gL = IVIG_C / IVIG_V;
double NER_mgL = NER_C / NER_V;

// Temperature asymmetry (degC): + warm limb (acute), - cold limb (chronic)
double TEMP_ASYM = TEMP_GAIN * (PERF - PERF_BASE);

// CRPS Severity Score surrogate (0-16; Harden 2010/2017 four-domain structure)
double css_pain  = 4.0 * (PAIN / PAIN_MAX);
double ta_abs    = fabs(TEMP_ASYM) / 3.0;
double css_vaso  = 4.0 * ((ta_abs > 1.0) ? 1.0 : ta_abs);
double css_edema = 4.0 * (EDEMA / (0.8 + EDEMA));
double css_motor = 4.0 * ((ROM_MAX - ROM) / ROM_MAX);
double CSS = css_pain + css_vaso + css_edema + css_motor;
if(CSS > 16.0) CSS = 16.0;

double ACTIVE_CRPS = (CSS >= 5.0 && PAIN >= 3.0) ? 1.0 : 0.0;
double REMISSION   = (CSS <  5.0 && PAIN <  3.0) ? 1.0 : 0.0;
double PHENOTYPE   = (TEMP_ASYM > 0.5) ? 1.0 : ((TEMP_ASYM < -0.5) ? -1.0 : 0.0);
// 1 = glial latch engaged (attractor 2), 0 = not engaged
double LATCHED     = (GLIA > GLIA50) ? 1.0 : 0.0;

// Local gain of the PAIN -> DISUSE -> CORTEX -> PAIN ring, evaluated as the
// product of the three steady-state sensitivities at the current point.
// > 1 means the ring can hold itself up without peripheral input.
double dDIS_dPAIN = (KIN_DIS / KOUT_DIS) * KFEAR * HILL_FEAR
                    * pow(PAIN50_FEAR, HILL_FEAR) * pow(PAIN, HILL_FEAR - 1.0)
                    / pow(pow(PAIN50_FEAR, HILL_FEAR)
                          + pow(PAIN, HILL_FEAR), 2.0);
double dCTX_dDIS  = (KIN_CTX / KOUT_CTX) * W_DIS_CTX * HILL_CTX
                    * pow(DIS50_CTX, HILL_CTX) * pow(DISUSE, HILL_CTX - 1.0)
                    / pow(pow(DIS50_CTX, HILL_CTX)
                          + pow(DISUSE, HILL_CTX), 2.0);
// pain_drive is local to $ODE, so it is recomputed here from the states
double drive_now  = W_PS_PAIN * PSENS + W_SS_PAIN * SSENS
                    + W_CTX_PAIN * CORTEX + W_DIS_PAIN * DISUSE;
double dPAIN_dCTX = PAIN_MAX * PAIN50_DRIVE * W_CTX_PAIN
                    / pow(PAIN50_DRIVE + drive_now, 2.0);
double RING_GAIN  = dDIS_dPAIN * dCTX_dDIS * dPAIN_dCTX;

$CAPTURE
KET_ng NORK_ng PRED_ng GBP_ng AMT_ng NAC_ng IVIG_gL NER_mgL
TEMP_ASYM CSS ACTIVE_CRPS REMISSION PHENOTYPE LATCHED RING_GAIN

$SET delta = 6, end = 17520     // 6-hourly output, 2-year default horizon

## =============================================================================
## R HARNESS — scenarios and analyses
##
##   library(mrgsolve)
##   mod <- mread("crps_mrgsolve_model.R")
##   e   <- mod@envir
##   e$run_scenarios(mod)            # 9 treatment scenarios, summary table
##   e$CRPS_trait_bifurcation(mod)   # KFEAR switch: 3 attractors
##   e$CRPS_insult_scan(mod)         # severity tier vs chronicity
##   e$CRPS_window_scan(mod)         # therapeutic window t*
##   e$CRPS_arm_decomposition(mod)   # which arm carries the window
##   e$CRPS_ketamine_washout(mod)    # Sigtermans-style washout
##   e$CRPS_dose_vs_timing(mod)      # timing beats dose (and the failed
##                                   #   supra-additivity hypothesis)
##   e$CRPS_phenotype_ordering(mod)  # warm-early vs cold-late arm ranking
##   e$CRPS_bone_axis(mod)           # neridronate / Varenna
##   e$CRPS_scs_habituation(mod)     # Kemler 2000 -> 2008
##   e$CRPS_report(mod)              # all of the above
##
## The harness deliberately uses only the most stable parts of the mrgsolve
## API (param(), ev(), as.ev() on a data frame, mrgsim(..., events=, end=,
## delta=)) so that it runs unchanged across mrgsolve versions.
## =============================================================================

$ENV

DAY  <- 24
YEAR <- 8760

## ---------------------------------------------------------------------------
## event helpers
## ---------------------------------------------------------------------------

## Combine several ev objects into one, ordered by time.
comb_ev <- function(...) {
  parts <- Filter(Negate(is.null), list(...))
  if (!length(parts)) return(NULL)
  d <- do.call(rbind, lapply(parts, as.data.frame))
  d <- d[order(d$time), , drop = FALSE]
  mrgsolve::as.ev(d)
}

## S-ketamine style continuous IV infusion (Sigtermans 2009): 100 h,
## titrated in four steps to `rate_max` mg/h (mean achieved ~22 mg/h).
ev_ketamine <- function(start = 30 * DAY, hours = 100, rate_max = 22) {
  steps <- c(0.25, 0.50, 0.75, 1.00)
  dur   <- hours / length(steps)
  d <- do.call(rbind, lapply(seq_along(steps), function(i) {
    r <- rate_max * steps[i]
    data.frame(time = start + (i - 1) * dur, cmt = 1, amt = r * dur,
               rate = r, evid = 1, ii = 0, addl = 0)
  }))
  mrgsolve::as.ev(d)
}

## Oral prednisolone 40 mg x 14 d, 20 mg x 7 d, 10 mg x 7 d
## (Kalita 2006 / Christensen 1982 structure)
ev_prednisolone <- function(start = 7 * DAY) {
  comb_ev(
    mrgsolve::ev(time = start,               amt = 40, cmt = "PRED_GUT", ii = 24, addl = 13),
    mrgsolve::ev(time = start + 14 * DAY,    amt = 20, cmt = "PRED_GUT", ii = 24, addl = 6),
    mrgsolve::ev(time = start + 21 * DAY,    amt = 10, cmt = "PRED_GUT", ii = 24, addl = 6))
}

## Neridronate 100 mg IV on days 0, 3, 6, 9 (Varenna 2013)
ev_neridronate <- function(start = 30 * DAY) {
  mrgsolve::ev(time = start, amt = 100, cmt = "NER_C", ii = 72, addl = 3)
}

## Gabapentin 600 mg TID (van de Vusse 2004)
ev_gabapentin <- function(start = 30 * DAY, days = 180, dose = 600) {
  mrgsolve::ev(time = start, amt = dose, cmt = "GBP_GUT", ii = 8,
               addl = 3 * days - 1)
}

## Amitriptyline 50 mg nightly
ev_amitriptyline <- function(start = 30 * DAY, days = 180, dose = 50) {
  mrgsolve::ev(time = start, amt = dose, cmt = "AMT_GUT", ii = 24,
               addl = days - 1)
}

## N-acetylcysteine 600 mg TID (Perez 2003; also stands in for DMSO 50% topical)
ev_nac <- function(start = 7 * DAY, days = 180, dose = 600) {
  mrgsolve::ev(time = start, amt = dose, cmt = "NAC_GUT", ii = 8,
               addl = 3 * days - 1)
}

## IVIG 0.5 g/kg (35 g) every 4 weeks x 6 (Goebel 2010 / LIPS 2017)
ev_ivig <- function(start = 30 * DAY, months = 6, g = 35) {
  mrgsolve::ev(time = start, amt = g, cmt = "IVIG_C", ii = 28 * DAY,
               addl = months - 1)
}

## ---------------------------------------------------------------------------
## simulation helpers
## ---------------------------------------------------------------------------
sim <- function(mod, pars = list(), events = NULL, end = 3 * 8760, delta = 6) {
  m <- mod
  if (length(pars)) m <- do.call(mrgsolve::param, c(list(m), pars))
  if (is.null(events)) {
    out <- mrgsolve::mrgsim(m, end = end, delta = delta)
  } else {
    out <- mrgsolve::mrgsim(m, events = events, end = end, delta = delta)
  }
  as.data.frame(out)
}

at_day <- function(d, day, col = "PAIN") {
  i <- which.min(abs(d$time - day * 24))
  d[[col]][i]
}

final_row <- function(d) d[nrow(d), ]

mean_nrs <- function(d, from_day, to_day) {
  s <- d[d$time >= from_day * 24 & d$time <= to_day * 24, ]
  mean(s$PAIN)
}

summarise_run <- function(d) {
  f <- final_row(d)
  data.frame(NRS = round(f$PAIN, 2), CSS = round(f$CSS, 2),
             ROM = round(f$ROM, 3), BMD = round(f$BMD, 3),
             GLIA = round(f$GLIA, 3), CORTEX = round(f$CORTEX, 3),
             DISUSE = round(f$DISUSE, 3), TEMP = round(f$TEMP_ASYM, 2),
             LATCHED = f$LATCHED, REMISSION = f$REMISSION)
}

## Multimodal package used throughout: prednisolone taper + NAC + rehabilitation
package <- function(start_day, arms = c("pred", "nac", "rehab")) {
  st <- start_day * DAY
  pars <- list(); evs <- list()
  if ("rehab" %in% arms)
    pars <- c(pars, list(REHAB = 1, REHAB_T0 = st, REHAB_DUR = 180 * DAY))
  if ("pred" %in% arms) evs <- c(evs, list(ev_prednisolone(st)))
  if ("nac"  %in% arms) evs <- c(evs, list(ev_nac(st, days = 180)))
  list(pars = pars, events = if (length(evs)) do.call(comb_ev, evs) else NULL)
}

## ---------------------------------------------------------------------------
## SCENARIOS (9)
## ---------------------------------------------------------------------------
run_scenarios <- function(mod, end = 3 * 8760) {

  out <- list()

  ## 1. untreated natural history: distal radius fracture, 30-day cast,
  ##    vulnerable trait (KFEAR = 1.0)
  out$untreated <- sim(mod, end = end)

  ## 2. EARLY multimodal package, day 7
  p <- package(7);   out$early_d7   <- sim(mod, p$pars, p$events, end)

  ## 3. LATE multimodal package, day 240 (identical package)
  p <- package(240); out$late_d240  <- sim(mod, p$pars, p$events, end)

  ## 4. Rehabilitation alone, day 30 (the arm that carries the window)
  out$rehab_alone <- sim(mod, list(REHAB = 1, REHAB_T0 = 30 * DAY,
                                   REHAB_DUR = 180 * DAY), NULL, end)

  ## 5. Prednisolone alone, day 7 (anti-inflammatory arm only)
  out$steroid_alone <- sim(mod, list(), ev_prednisolone(7 * DAY), end)

  ## 6. Neridronate 100 mg IV x4 from day 30 (Varenna 2013)
  out$neridronate <- sim(mod, list(), ev_neridronate(30 * DAY), end)

  ## 7. Ketamine 100-h infusion alone, day 30 (Sigtermans 2009)
  out$ketamine <- sim(mod, list(), ev_ketamine(30 * DAY), end)

  ## 8. Ketamine + rehabilitation started in the analgesic window
  out$ketamine_rehab <- sim(mod, list(REHAB = 1, REHAB_T0 = 30 * DAY,
                                      REHAB_DUR = 180 * DAY),
                            ev_ketamine(30 * DAY), end)

  ## 9. Cold, longstanding phenotype: high sympathetic tone, SCS at month 6
  ##    plus antioxidant and a PDE5-inhibitor-like vasodilator
  out$cold_scs <- sim(mod, list(SYMP_TONE = 1.7, VASODIL = 1,
                                SCS_ON = 1, SCS_T0 = 180 * DAY),
                      ev_nac(180 * DAY, days = 180), end)

  tab <- do.call(rbind, lapply(names(out),
                               function(n) cbind(scenario = n, summarise_run(out[[n]]))))
  cat("\n=== 9 SCENARIOS, state at 3 years ===\n"); print(tab)
  attr(out, "summary") <- tab
  invisible(out)
}

## ---------------------------------------------------------------------------
## 1. TRAIT BIFURCATION — the fear-avoidance gain is the switch
## ---------------------------------------------------------------------------
CRPS_trait_bifurcation <- function(mod, kfears = c(0.2, 0.35, 0.5, 0.6, 0.65,
                                                   0.7, 0.75, 0.8, 0.9, 1.0,
                                                   1.2, 1.4),
                                   end = 3 * 8760) {
  rows <- lapply(kfears, function(k) {
    d <- sim(mod, list(KFEAR = k), NULL, end, delta = 24)
    cbind(KFEAR = k, summarise_run(d))
  })
  res <- do.call(rbind, rows)
  cat("\n=== 1. BIFURCATION IN THE FEAR-AVOIDANCE TRAIT (3-year state) ===\n")
  print(res)
  jumps <- which(diff(res$NRS) > 1)
  cat("Step transitions between KFEAR = ",
      paste(sprintf("%.2f->%.2f", res$KFEAR[jumps], res$KFEAR[jumps + 1]),
            collapse = ", "), "\n")
  cat("The injury is identical in every row. A psychological gain, not a\n")
  cat("nociceptive one, decides both whether the disease persists and which\n")
  cat("of the two chronic states it settles into.\n")
  invisible(res)
}

## ---------------------------------------------------------------------------
## 2. INSULT MAGNITUDE — sets the tier, not the chronicity
## ---------------------------------------------------------------------------
CRPS_insult_scan <- function(mod, amps = c(0.1, 0.2, 0.35, 0.5, 0.7, 0.85,
                                           0.9, 0.95, 1.0, 1.2, 1.5),
                             end = 3 * 8760) {
  rows <- lapply(amps, function(a) {
    d <- sim(mod, list(INJ_AMP = a), NULL, end, delta = 24)
    cbind(INJ_AMP = a, peak_NRS = round(max(d$PAIN), 2), summarise_run(d))
  })
  res <- do.call(rbind, rows)
  cat("\n=== 2. INSULT-MAGNITUDE SCAN (3-year state) ===\n")
  print(res)
  cat("Chronicity is flat in the insult amplitude; only the severity tier\n")
  cat("steps (attractor 1 -> attractor 2) once the afferent drive is large\n")
  cat("enough to push GLIA through GLIA50.\n")
  invisible(res)
}

## ---------------------------------------------------------------------------
## 3. THERAPEUTIC WINDOW — identical package, varying start day
## ---------------------------------------------------------------------------
CRPS_window_scan <- function(mod, delays = c(3, 7, 14, 21, 30, 45, 60, 75, 90,
                                             95, 100, 120, 180, 240, 365),
                             end = 3 * 8760) {
  rows <- lapply(delays, function(dd) {
    p <- package(dd)
    d <- sim(mod, p$pars, p$events, end, delta = 24)
    cbind(start_day = dd, summarise_run(d))
  })
  res <- do.call(rbind, rows)
  cat("\n=== 3. THERAPEUTIC WINDOW (identical multimodal package) ===\n")
  print(res)
  bad <- which(res$NRS > (min(res$NRS) + max(res$NRS)) / 2)
  tstar <- if (length(bad)) res$start_day[bad[1]] else NA
  cat(sprintf("Critical delay t* between %s and %s days.\n",
              ifelse(is.na(tstar), "-", res$start_day[max(bad[1] - 1, 1)]),
              ifelse(is.na(tstar), "-", tstar)))
  cat("Nothing about the package changes across rows. What changes is whether\n")
  cat("GLIA has already crossed its threshold when the package arrives.\n")
  invisible(res)
}

## ---------------------------------------------------------------------------
## 4. WHICH ARM CARRIES THE WINDOW?
## ---------------------------------------------------------------------------
CRPS_arm_decomposition <- function(mod, days = c(7, 30, 60, 120, 240),
                                   end = 3 * 8760) {
  arms <- list(prednisolone = "pred", NAC = "nac", rehabilitation = "rehab",
               pred_plus_nac = c("pred", "nac"),
               full_package = c("pred", "nac", "rehab"))
  m <- matrix(NA_real_, nrow = length(days), ncol = length(arms),
              dimnames = list(paste0("d", days), names(arms)))
  for (i in seq_along(days)) for (j in seq_along(arms)) {
    p <- package(days[i], arms[[j]])
    d <- sim(mod, p$pars, p$events, end, delta = 24)
    m[i, j] <- round(final_row(d)$PAIN, 2)
  }
  cat("\n=== 4. ARM DECOMPOSITION (3-year NRS) ===\n"); print(m)
  cat("Rehabilitation is the only arm that reaches resolution. The\n")
  cat("anti-inflammatory arms buy a tier downgrade and only if given early:\n")
  cat("their measurable benefit is the prevention of a LATE transition, so a\n")
  cat("trial that enrols patients who will never make that transition, or who\n")
  cat("already have, must read out as null however potent the drug is.\n")
  invisible(m)
}

## ---------------------------------------------------------------------------
## 5. KETAMINE WASHOUT (Sigtermans 2009 reproduction)
## ---------------------------------------------------------------------------
CRPS_ketamine_washout <- function(mod, start_day = 30, end = 2 * 8760) {
  st <- start_day * DAY
  runs <- list(
    untreated  = sim(mod, list(), NULL, end),
    ketamine   = sim(mod, list(), ev_ketamine(st), end),
    rehab      = sim(mod, list(REHAB = 1, REHAB_T0 = st, REHAB_DUR = 180 * DAY),
                     NULL, end),
    ket_rehab  = sim(mod, list(REHAB = 1, REHAB_T0 = st, REHAB_DUR = 180 * DAY),
                     ev_ketamine(st), end))
  days <- c(29, 31, 32, 34, 37, 44, 58, 79, 107, 210, 395, 700)
  tab <- data.frame(day = days)
  for (n in names(runs))
    tab[[n]] <- vapply(days, function(d) round(at_day(runs[[n]], d), 2), numeric(1))
  cat("\n=== 5. KETAMINE 100-h INFUSION: NRS trajectory ===\n"); print(tab)
  cat(sprintf("ketamine Cmax %.0f ng/mL, norketamine Cmax %.0f ng/mL\n",
              max(runs$ketamine$KET_ng), max(runs$ketamine$NORK_ng)))
  e_k <- tail(tab$untreated, 1) - tail(tab$ketamine, 1)
  e_r <- tail(tab$untreated, 1) - tail(tab$rehab, 1)
  e_kr <- tail(tab$untreated, 1) - tail(tab$ket_rehab, 1)
  cat(sprintf("700-day effects: ketamine %+.2f, rehab %+.2f, both %+.2f, interaction %+.2f\n",
              e_k, e_r, e_kr, e_kr - (e_k + e_r)))
  cat("The combination is clearly better EARLY and exactly equal LATE: the\n")
  cat("interaction at the endpoint is negative because rehabilitation alone\n")
  cat("already reaches resolution. The drug buys time, not a better ceiling.\n")
  invisible(tab)
}

## ---------------------------------------------------------------------------
## 6. DOSE vs TIMING
## ---------------------------------------------------------------------------
CRPS_dose_vs_timing <- function(mod, rates = c(5.5, 11, 22, 44, 88),
                                delays = c(14, 30, 60, 120, 240),
                                end = 2 * 8760) {
  m <- matrix(NA_real_, nrow = length(rates), ncol = length(delays),
              dimnames = list(paste0(rates, "mg/h"), paste0("d", delays)))
  cmax <- numeric(length(rates))
  for (i in seq_along(rates)) for (j in seq_along(delays)) {
    d <- sim(mod, list(), ev_ketamine(delays[j] * DAY, rate_max = rates[i]),
             end, delta = 24)
    m[i, j] <- round(final_row(d)$PAIN, 2)
    if (j == 1) cmax[i] <- round(max(d$KET_ng), 0)
  }
  cat("\n=== 6. DOSE vs TIMING (2-year NRS) ===\n")
  print(cbind(m, Cmax_ng_mL = cmax))
  cat(sprintf("dose axis at the best start day: range %.2f NRS over a %gx dose span\n",
              diff(range(m[, 1])), max(rates) / min(rates)))
  cat(sprintf("timing axis at the highest dose: range %.2f NRS over d%g -> d%g\n",
              diff(range(m[nrow(m), ])), min(delays), max(delays)))
  cat("Above threshold the dose axis is flat and the timing axis is not.\n")
  invisible(m)
}

## ---------------------------------------------------------------------------
## 7. PHENOTYPE / TIMING ORDERING OF EVERY ARM
## ---------------------------------------------------------------------------
CRPS_phenotype_ordering <- function(mod, end = 2 * 8760) {
  arms <- c("none", "steroid", "antioxidant", "vasodilator", "bisphosphonate",
            "gabapentin", "amitriptyline", "rehab_GMI", "sympathetic_block",
            "IVIG", "ketamine_100h", "SCS")
  build <- function(arm, st) {
    pars <- list(); ev <- NULL
    if (arm == "steroid")            ev <- ev_prednisolone(st)
    else if (arm == "antioxidant")   ev <- ev_nac(st, days = 180)
    else if (arm == "vasodilator")   pars <- list(VASODIL = 1)
    else if (arm == "bisphosphonate")ev <- ev_neridronate(st)
    else if (arm == "gabapentin")    ev <- ev_gabapentin(st, days = 180)
    else if (arm == "amitriptyline") ev <- ev_amitriptyline(st, days = 180)
    else if (arm == "rehab_GMI")     pars <- list(REHAB = 1, REHAB_T0 = st,
                                                  REHAB_DUR = 180 * DAY)
    else if (arm == "sympathetic_block") pars <- list(SYMPBLOCK = 1, SB_T0 = st)
    else if (arm == "IVIG")          ev <- ev_ivig(st)
    else if (arm == "ketamine_100h") ev <- ev_ketamine(st)
    else if (arm == "SCS")           pars <- list(SCS_ON = 1, SCS_T0 = st)
    list(pars = pars, ev = ev)
  }
  conds <- list(warm_early = list(pars = list(SYMP_TONE = 0.8), day = 14),
                cold_late  = list(pars = list(SYMP_TONE = 1.7), day = 240))
  res <- list()
  for (cn in names(conds)) {
    cc <- conds[[cn]]; st <- cc$day * DAY
    vals <- lapply(arms, function(a) {
      b <- build(a, st)
      d <- sim(mod, utils::modifyList(cc$pars, b$pars), b$ev, end, delta = 24)
      c(mean180 = mean_nrs(d, cc$day, cc$day + 180),
        NRS_2y = final_row(d)$PAIN)
    })
    names(vals) <- arms
    base <- vals[["none"]]
    tab <- do.call(rbind, lapply(arms[-1], function(a)
      data.frame(arm = a,
                 d_mean_NRS_180d = round(base["mean180"] - vals[[a]]["mean180"], 2),
                 d_NRS_2y = round(base["NRS_2y"] - vals[[a]]["NRS_2y"], 2))))
    tab <- tab[order(-tab$d_mean_NRS_180d), ]
    rownames(tab) <- NULL
    cat(sprintf("\n=== 7. ARM ORDERING, %s (untreated mean NRS %.2f, 2-y %.2f) ===\n",
                cn, base["mean180"], base["NRS_2y"]))
    print(tab)
    res[[cn]] <- tab
  }
  cat("\nThe ordering is a property of the system state, not of the molecules:\n")
  cat("the loop-directed and device arms lead early, and in the longstanding\n")
  cat("cold limb every pharmacological arm collapses to symptomatic-only.\n")
  invisible(res)
}

## ---------------------------------------------------------------------------
## 8. BONE AXIS (neridronate, Varenna 2013/2017)
## ---------------------------------------------------------------------------
CRPS_bone_axis <- function(mod, end = 8760) {
  runs <- list(untreated = sim(mod, list(), NULL, end, delta = 12),
               nerid_d30 = sim(mod, list(), ev_neridronate(30 * DAY), end, delta = 12),
               nerid_d120 = sim(mod, list(), ev_neridronate(120 * DAY), end, delta = 12))
  days <- c(30, 40, 70, 120, 160, 250, 365)
  tab <- data.frame(day = days)
  for (n in names(runs))
    tab[[paste0("NRS_", n)]] <-
      vapply(days, function(d) round(at_day(runs[[n]], d), 2), numeric(1))
  bone <- data.frame(
    day = days,
    CTXI_untreated = vapply(days, function(d) round(at_day(runs$untreated, d, "CTXI"), 3), numeric(1)),
    CTXI_nerid     = vapply(days, function(d) round(at_day(runs$nerid_d30, d, "CTXI"), 3), numeric(1)),
    BMD_untreated  = vapply(days, function(d) round(at_day(runs$untreated, d, "BMD"), 3), numeric(1)),
    BMD_nerid      = vapply(days, function(d) round(at_day(runs$nerid_d30, d, "BMD"), 3), numeric(1)),
    bone_bound_mg  = vapply(days, function(d) round(at_day(runs$nerid_d30, d, "NER_BONE"), 1), numeric(1)))
  cat("\n=== 8. BONE AXIS (neridronate 100 mg IV x4) ===\n")
  print(tab); print(bone)
  cat("The analgesia here is not a bone endpoint: suppressing osteoclastic\n")
  cat("acidification removes W_BONE_PS from the peripheral drive. It behaves\n")
  cat("like the other peripheral-node agents, which is why the trials were\n")
  cat("positive in short-duration disease.\n")
  invisible(list(pain = tab, bone = bone))
}

## ---------------------------------------------------------------------------
## 9. SCS HABITUATION (Kemler 2000 positive -> 2008 lost)
## ---------------------------------------------------------------------------
CRPS_scs_habituation <- function(mod, implant_day = 180, end = 5 * 8760) {
  a <- sim(mod, list(), NULL, end, delta = 24)
  b <- sim(mod, list(SCS_ON = 1, SCS_T0 = implant_day * DAY), NULL, end, delta = 24)
  months <- c(6, 7, 9, 12, 18, 24, 36, 48, 60)
  tab <- data.frame(
    month = months,
    NRS_control = vapply(months, function(m) round(at_day(a, m * 30), 2), numeric(1)),
    NRS_SCS     = vapply(months, function(m) round(at_day(b, m * 30), 2), numeric(1)))
  tab$difference <- round(tab$NRS_control - tab$NRS_SCS, 2)
  cat("\n=== 9. SCS WITH HABITUATION ===\n"); print(tab)
  cat("Same device, same setting: a large early difference that decays as\n")
  cat("SCS_EFF * exp(-SCS_HAB_K * t) while the underlying attractor is\n")
  cat("untouched -- the shape of the Kemler 6-month vs 5-year discordance.\n")
  invisible(tab)
}

## ---------------------------------------------------------------------------
## everything
## ---------------------------------------------------------------------------
CRPS_report <- function(mod) {
  invisible(list(
    scenarios   = run_scenarios(mod),
    trait       = CRPS_trait_bifurcation(mod),
    insult      = CRPS_insult_scan(mod),
    window      = CRPS_window_scan(mod),
    arms        = CRPS_arm_decomposition(mod),
    ketamine    = CRPS_ketamine_washout(mod),
    dose_timing = CRPS_dose_vs_timing(mod),
    ordering    = CRPS_phenotype_ordering(mod),
    bone        = CRPS_bone_axis(mod),
    scs         = CRPS_scs_habituation(mod)))
}
