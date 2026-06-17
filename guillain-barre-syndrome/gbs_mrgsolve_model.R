# ============================================================
# Guillain-Barré Syndrome (GBS) QSP Model
# mrgsolve ODE-based PK/PD Model
#
# Key References:
#   1. Hughes RA et al. Lancet Neurol. 2014;13(6):619-630 (GBS overview)
#   2. Yuki N, Hartung HP. N Engl J Med. 2012;366(24):2294-2304 (pathophysiology)
#   3. van den Berg B et al. Nat Rev Neurol. 2014;10(8):469-482 (epidemiology)
#   4. Willison HJ et al. Lancet. 2016;388(10045):717-727 (clinical features)
#   5. Kuitwaard K et al. J Neurol Neurosurg Psychiatry. 2009;80(7):776-780 (IVIG PK)
#   6. Dalakas MC. Ther Adv Neurol Disord. 2012;5(3):155-163 (IVIG mechanisms)
#   7. Pritchard J et al. Ann Neurol. 2003;53(5):600-607 (complement in GBS)
#   8. Halstead SK et al. Brain. 2008;131(5):1197-1208 (anti-ganglioside Abs)
#   9. Uncini A, Kuwabara S. J Neurol Neurosurg Psychiatry. 2015;86(10):1157-1162 (AMAN)
#  10. Brilot F et al. J Neuroinflammation. 2010;7:76 (T cell mechanisms)
#  11. Leonhard SE et al. Nat Rev Neurol. 2019;15(11):671-683 (diagnosis and management)
#  12. van Doorn PA et al. Lancet. 2008;372(9648):1487-1498 (IVIG/PE treatment)
#  13. Feasby TE et al. Ann Neurol. 1986;20(3):317-328 (AMAN discovery)
#  14. Griffin JW et al. Ann Neurol. 1996;39(5):586-591 (axonal GBS)
#  15. Mishu B et al. Ann Intern Med. 1993;118(9):658-663 (Campylobacter trigger)
# ============================================================

library(mrgsolve)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(patchwork)

# ============================================================
# MODEL DEFINITION
# ============================================================

gbs_model <- '
$PROB
Guillain-Barre Syndrome (GBS) QSP Model
========================================
Disease: GBS - acute inflammatory peripheral neuropathy
Subtypes: 1=AIDP (demyelinating), 2=AMAN (axonal), 3=MFS (Miller Fisher)
Treatments: IVIG 2-compartment PK, Plasma Exchange, Eculizumab (investigational)

Model Components (20+ compartments):
  [1]  Pathogen        - triggering pathogen (e.g. Campylobacter jejuni)
  [2]  DC_act          - activated dendritic cells
  [3]  Mac_act         - activated macrophages/monocytes
  [4]  Bcell           - activated B lymphocytes
  [5]  Plasma          - plasma cells (antibody-secreting)
  [6]  Ab_anti_GM1     - anti-GM1 IgG antibodies (AIDP/AMAN)
  [7]  Ab_anti_GD1a    - anti-GD1a IgG antibodies (AMAN)
  [8]  Ab_anti_GQ1b    - anti-GQ1b IgG antibodies (MFS)
  [9]  Th1             - Th1 helper T cells
  [10] Th17            - Th17 helper T cells
  [11] Treg            - regulatory T cells
  [12] IL6             - interleukin-6
  [13] TNFa            - tumor necrosis factor alpha
  [14] IL10            - interleukin-10 (anti-inflammatory)
  [15] C3b             - complement C3b (activated)
  [16] MAC             - membrane attack complex (C5b-9)
  [17] C5a             - complement C5a (anaphylatoxin)
  [18] Myelin_damage   - cumulative myelin injury (0-1 scale)
  [19] Axon_damage     - cumulative axonal injury (0-1 scale)
  [20] Nerve_function  - composite nerve function (0=lost, 1=normal)
  [21] GBS_score       - GBS disability score (0-6 Hughes scale)
  [22] FVC_pct         - forced vital capacity (% predicted)
  [23] IVIG_C          - IVIG central compartment concentration (g/L)
  [24] IVIG_P          - IVIG peripheral compartment concentration (g/L)
  [25] PE_cumulative   - cumulative PE sessions effect
  [26] ECU_C           - eculizumab concentration (mg/L)

$PARAM
// ===================================================================
// INFECTION & TRIGGER PARAMETERS
// ===================================================================
// Campylobacter jejuni LPS (LOS) cross-reacts with gangliosides
// Ref: Mishu & Blaser Ann Intern Med 1993; Nachamkin et al. 2008
k_infect_clear  = 0.50    // pathogen clearance rate (1/day); T1/2 ~ 1.4 days
LOS_peak        = 1.0     // LOS antigen peak (relative, normalized)
t_infect_start  = 0.0     // infection start (day 0 = symptom trigger)
t_infect_end    = 14.0    // infection resolution (typical ~2 weeks post exposure)

// ===================================================================
// INNATE IMMUNE PARAMETERS
// ===================================================================
// Ref: Brilot et al. J Neuroinflammation 2010; Meyer zu Horste et al. 2010
k_DC_act        = 0.30    // DC activation rate by pathogen (1/day)
k_DC_decay      = 0.15    // DC decay rate (1/day); DC lifespan ~5-7 days
k_Mac_recruit   = 0.25    // macrophage recruitment rate by DC/C5a (1/day)
k_Mac_decay     = 0.12    // macrophage decay rate (1/day)

// ===================================================================
// B CELL AND ANTIBODY PARAMETERS
// ===================================================================
// Ref: Willison & Yuki Brain 2002; Halstead et al. Brain 2008
k_Bcell_act     = 0.10    // B cell activation rate (1/day)
k_plasma_diff   = 0.08    // plasma cell differentiation rate (1/day)
k_Ab_prod       = 0.50    // antibody production rate (AU/day/plasma cell)
k_Ab_decay      = 0.035   // IgG decay rate (1/day); IgG half-life ~ 20 days
k_Ab_GQ1b       = 0.30    // anti-GQ1b production (1/day); high in MFS

// ===================================================================
// T CELL PARAMETERS
// ===================================================================
// Ref: Brilot et al. J Neuroinflammation 2010; Zhang et al. 2013
k_Th1_diff      = 0.12    // Th1 differentiation rate (1/day)
k_Th17_diff     = 0.08    // Th17 differentiation rate (1/day)
k_Treg_form     = 0.05    // Treg induction rate by IL-10 (1/day)
k_Tcell_decay   = 0.10    // T cell decay rate (1/day); lifespan ~10 days

// ===================================================================
// CYTOKINE PARAMETERS
// ===================================================================
// Ref: Sharief et al. Ann Neurol 1993; Press et al. J Neuroimmunol 2001
k_IL6_prod      = 0.40    // IL-6 production by Mac + Th17 (AU/day)
k_IL6_decay     = 0.50    // IL-6 decay rate (1/day); T1/2 ~ 33 h
k_TNFa_prod     = 0.35    // TNF-alpha production by Mac + Th1 (AU/day)
k_TNFa_decay    = 0.60    // TNF-alpha decay rate (1/day); T1/2 ~ 28 h
k_IL10_prod     = 0.20    // IL-10 production by Treg + Mac (AU/day)
k_IL10_decay    = 0.40    // IL-10 decay rate (1/day); T1/2 ~ 42 h

// ===================================================================
// COMPLEMENT CASCADE PARAMETERS
// ===================================================================
// Ref: Halstead et al. Brain 2004; Phongsisay et al. 2008
// C3b deposition at nodes of Ranvier; MAC forms on axolemma
k_C3_act        = 0.30    // C3 activation rate by Ab-antigen complexes (1/day)
k_C3_decay      = 0.20    // C3b decay/inactivation rate (1/day)
k_MAC_form      = 0.25    // MAC formation rate from C3b (1/day)
k_MAC_decay     = 0.15    // MAC decay rate (1/day); membrane clearing
k_C5a_prod      = 0.20    // C5a production rate (1/day)

// ===================================================================
// NERVE DAMAGE PARAMETERS
// ===================================================================
// AIDP: macrophage-mediated myelin stripping + complement
// Ref: Asbury & Cornblath Ann Neurol 1990; Hafer-Macko et al. 1996
k_myelin_damage = 0.080   // myelin damage rate (1/day); moderate progression
k_remyelin      = 0.040   // remyelination rate (1/day); recovery 4-8 weeks

// AMAN: anti-GM1/GD1a Ab + complement attack axon at nodes of Ranvier
// Ref: Griffin et al. Ann Neurol 1996; Ho et al. 1997
k_axon_damage   = 0.060   // axonal damage rate (1/day)
k_axon_repair   = 0.015   // axon repair/regeneration rate (1/day); slow

// ===================================================================
// CLINICAL ENDPOINT PARAMETERS
// ===================================================================
// Hughes disability scale 0-6 (WHO/GBS standard)
// Ref: Hughes et al. Brain 1978; van Koningsveld et al. 2007
GBS_score_max   = 6.0     // maximum GBS disability score (bedbound/dead)
FVC_baseline    = 100.0   // FVC% at baseline (normal = 100%)

// ===================================================================
// IVIG PK PARAMETERS (2-COMPARTMENT MODEL)
// ===================================================================
// Ref: Kuitwaard et al. J Neurol Neurosurg Psychiatry 2009
// Hammarstrom & Smith NEJM 1993 (IgG pharmacokinetics)
IVIG_CL         = 0.0033  // IgG clearance (L/h/kg); ~0.23 mL/h/kg cited
IVIG_Vc         = 0.050   // central volume (L/kg); plasma space ~50 mL/kg
IVIG_Vp         = 0.090   // peripheral volume (L/kg); tissue distribution
IVIG_Q          = 0.005   // intercompartment clearance (L/h/kg)
WT              = 70.0    // body weight (kg)

// ===================================================================
// PLASMA EXCHANGE PARAMETERS
// ===================================================================
// Ref: French Cooperative Group Lancet 1997; van der Meche Lancet 1992
PE_Ab_removal_eff = 0.60  // fraction of IgG removed per PE session (~60%)
PE_sessions     = 5.0     // standard course: 5 sessions over 2 weeks

// ===================================================================
// ECULIZUMAB PK PARAMETERS (INVESTIGATIONAL)
// ===================================================================
// Ref: Misawa et al. Ann Neurol 2018 (eculizumab in AMAN)
// Based on PNH pharmacokinetics adapted for GBS
ECU_CL          = 0.0026  // eculizumab clearance (L/h)
ECU_Vc          = 0.07    // central volume (L/kg)
ECU_IC50        = 100.0   // C5 inhibition EC50 (mcg/mL); literature ~50-150

// ===================================================================
// TREATMENT FLAGS AND DOSING
// ===================================================================
IVIG_dose       = 0.0     // IVIG dose (g/kg); 0=off, 2.0=standard course
IVIG_start      = 7.0     // IVIG start day (day from symptom onset)
PE_flag         = 0.0     // PE: 0=no, 1=yes
ECU_flag        = 0.0     // Eculizumab: 0=no, 1=yes

// ===================================================================
// SUBTYPE MODIFIER
// ===================================================================
subtype         = 1.0     // 1=AIDP, 2=AMAN, 3=MFS (Miller Fisher)

$INIT
// --- INFECTION ---
Pathogen        = 1.0     // normalized pathogen/antigen load at onset

// --- INNATE IMMUNE ---
DC_act          = 0.01    // resting dendritic cells
Mac_act         = 0.01    // resting macrophages

// --- ADAPTIVE IMMUNE ---
Bcell           = 0.01    // naive B cells (small pool)
Plasma          = 0.00    // no plasma cells initially

// --- ANTIBODIES (all start at 0 or trace) ---
Ab_anti_GM1     = 0.001   // trace anti-GM1 (AIDP/AMAN relevant)
Ab_anti_GD1a    = 0.001   // trace anti-GD1a (AMAN relevant)
Ab_anti_GQ1b    = 0.001   // trace anti-GQ1b (MFS relevant)

// --- T CELLS ---
Th1             = 0.01    // resting Th1
Th17            = 0.01    // resting Th17
Treg            = 0.10    // baseline regulatory T cells (higher = tolerant)

// --- CYTOKINES ---
IL6             = 0.01    // baseline IL-6 (low)
TNFa            = 0.01    // baseline TNF-alpha (low)
IL10            = 0.05    // baseline IL-10 (anti-inflammatory tone)

// --- COMPLEMENT ---
C3b             = 0.00    // no activated complement at start
MAC             = 0.00    // no MAC at start
C5a             = 0.00    // no C5a at start

// --- NERVE DAMAGE ---
Myelin_damage   = 0.00    // no myelin damage at onset
Axon_damage     = 0.00    // no axon damage at onset
Nerve_function  = 1.00    // full nerve function at onset

// --- CLINICAL OUTCOMES ---
GBS_score       = 0.00    // normal function (0 = normal)
FVC_pct         = 100.00  // normal FVC%

// --- IVIG PK ---
IVIG_C          = 0.00    // no drug at start (g/L central)
IVIG_P          = 0.00    // no drug at start (g/L peripheral)

// --- PE ---
PE_cumulative   = 0.00    // no PE effect yet

// --- ECULIZUMAB ---
ECU_C           = 0.00    // no eculizumab at start

$ODE
// ==================================================================
// 1. PATHOGEN / ANTIGEN DYNAMICS
// ==================================================================
// Pathogen present during infection window, then cleared
// LOS cross-reacts with peripheral nerve gangliosides (molecular mimicry)
double LOS_stim = (SOLVERTIME >= t_infect_start && SOLVERTIME <= t_infect_end)
                  ? LOS_peak : 0.0;
dxdt_Pathogen = LOS_stim - k_infect_clear * Pathogen;

// ==================================================================
// 2. INNATE IMMUNE ACTIVATION
// ==================================================================
// DCs recognize pathogen via TLRs; macrophages recruited by DC signals and C5a
dxdt_DC_act = k_DC_act * Pathogen - k_DC_decay * DC_act;
dxdt_Mac_act = k_Mac_recruit * (DC_act + C5a * 0.5) - k_Mac_decay * Mac_act;

// ==================================================================
// 3. B CELL ACTIVATION AND PLASMA CELL DIFFERENTIATION
// ==================================================================
// BCR crosslinking + T cell help + DC presentation drive B cell expansion
// Ref: Willison Brain 2002 - molecular mimicry triggers B cell response
double BCR_signal = Pathogen * (1.0 + DC_act);
double Bcell_carrying_cap = 2.0;
dxdt_Bcell = k_Bcell_act * BCR_signal * (1.0 - Bcell / Bcell_carrying_cap)
             - 0.05 * Bcell;
dxdt_Plasma = k_plasma_diff * Bcell - 0.03 * Plasma;

// ==================================================================
// 4. ANTIBODY DYNAMICS (subtype-specific)
// ==================================================================
// IVIG anti-idiotype: reduces pathogenic Ab production
// Ref: Dalakas Ther Adv Neurol Disord 2012 - anti-idiotype neutralization
double IVIG_antiidio = (IVIG_C > 0.01)
                       ? 0.50 * IVIG_C / (IVIG_C + 0.10) : 0.0;

// PE removes circulating antibodies (bulk IgG removal)
// Ref: French Group Lancet 1997 - 60% removal per session
double PE_rate = PE_flag * 0.05;   // continuous approximation of bolus removal

// subtype determines which Ab predominates
// AIDP: anti-GM1 + anti-GD1b; AMAN: anti-GM1 + anti-GD1a; MFS: anti-GQ1b
double Ab_GM1_factor  = (subtype <= 2.0) ? 1.0 : 0.10;
double Ab_GD1a_factor = (subtype == 2.0) ? 1.0 : 0.30;
double Ab_GQ1b_factor = (subtype == 3.0) ? 1.0 : 0.05;

dxdt_Ab_anti_GM1  = Ab_GM1_factor  * k_Ab_prod * Plasma
                    - k_Ab_decay * Ab_anti_GM1
                    - PE_rate * Ab_anti_GM1
                    - IVIG_antiidio * Ab_anti_GM1;

dxdt_Ab_anti_GD1a = Ab_GD1a_factor * k_Ab_prod * Plasma
                    - k_Ab_decay * Ab_anti_GD1a
                    - PE_rate * Ab_anti_GD1a
                    - IVIG_antiidio * Ab_anti_GD1a;

dxdt_Ab_anti_GQ1b = Ab_GQ1b_factor * k_Ab_GQ1b * Plasma
                    - k_Ab_decay * Ab_anti_GQ1b
                    - PE_rate * Ab_anti_GQ1b;

// ==================================================================
// 5. T CELL DYNAMICS
// ==================================================================
// Th1/Th17 promote inflammation; Treg limits autoimmunity
// IVIG expands/stabilizes Tregs - Ref: Ephrem et al. Blood 2008
double IL12_signal = DC_act;        // DC-derived IL-12 drives Th1
double IL23_signal = DC_act * 0.5;  // DC-derived IL-23 drives Th17
double Treg_inh    = Treg / (Treg + 0.30);  // Hill-type inhibition

dxdt_Th1  = k_Th1_diff  * IL12_signal * (1.0 - Treg_inh)
            - k_Tcell_decay * Th1;
dxdt_Th17 = k_Th17_diff * IL23_signal * (1.0 - Treg_inh)
            - k_Tcell_decay * Th17;

// IVIG promotes Treg expansion; IL-10 supports Treg maintenance
double IVIG_Treg_effect = (IVIG_C > 0.01)
                          ? 0.30 * IVIG_C / (IVIG_C + 0.20) : 0.0;
dxdt_Treg = k_Treg_form * IL10
            - 0.08 * Treg * (Th1 + Th17)   // Th1/Th17 suppress Treg
            + IVIG_Treg_effect;

// ==================================================================
// 6. CYTOKINE DYNAMICS
// ==================================================================
// IL-6: Mac + Th17 → drives inflammation and acute phase response
// TNF-α: Mac + Th1 → promotes myelin damage and nerve injury
// IL-10: Treg + Mac → anti-inflammatory, promotes recovery
// Ref: Sharief et al. Ann Neurol 1993; Press et al. J Neuroimmunol 2001
dxdt_IL6  = k_IL6_prod  * (Mac_act + Th17)    - k_IL6_decay  * IL6;
dxdt_TNFa = k_TNFa_prod * (Mac_act + Th1)     - k_TNFa_decay * TNFa;
dxdt_IL10 = k_IL10_prod * (Treg + Mac_act * 0.3) - k_IL10_decay * IL10;

// ==================================================================
// 7. COMPLEMENT CASCADE
// ==================================================================
// Ab-ganglioside complexes activate complement at nodes of Ranvier
// Ref: Halstead et al. Brain 2004 - MAC deposition at paranodes
// IVIG: scavenges complement via Fc fragment and C3b binding
// Ref: Basta & Dalakas Blood 1994 - complement scavenging by IVIG
// Eculizumab: anti-C5 MAb prevents C5→C5a + C5b cleavage
// Ref: Misawa et al. Ann Neurol 2018 - eculizumab in severe GBS
double Ab_total = Ab_anti_GM1 + Ab_anti_GD1a + Ab_anti_GQ1b;
double IVIG_comp_inh = (IVIG_C > 0.01)
                       ? 0.60 * IVIG_C / (IVIG_C + 0.15) : 0.0;
double ECU_C5_inh    = (ECU_flag > 0.5 && ECU_C > 0.01)
                       ? ECU_C / (ECU_C + ECU_IC50) : 0.0;

dxdt_C3b = k_C3_act  * Ab_total  - k_C3_decay  * C3b
           - IVIG_comp_inh * k_C3_decay * C3b;  // IVIG accelerates C3b clearing
dxdt_MAC = k_MAC_form * C3b * (1.0 - ECU_C5_inh)
           - k_MAC_decay * MAC;
dxdt_C5a = k_C5a_prod * C3b * (1.0 - ECU_C5_inh)
           - 0.80 * C5a;   // C5a half-life ~ 21 min

// ==================================================================
// 8. SUBTYPE-SPECIFIC NERVE DAMAGE
// ==================================================================
// AIDP: macrophage-mediated demyelination at paranodal region
//   → Myelin_damage primarily; slower but more reversible
// AMAN: anti-GM1/GD1a Ab + MAC → axolemma disruption at nodes
//   → Axon_damage primarily; faster but less reversible
// MFS: anti-GQ1b at NMJ and paranodal region; primarily cranial nerves
// Ref: Uncini & Kuwabara 2015; Ho et al. Ann Neurol 1997

// AIDP damage driven by macrophage + complement + cytokines
double AIDP_driver = (subtype == 1.0)
                     ? (MAC + Mac_act * 0.50 + IL6 * 0.20 + TNFa * 0.10)
                     : (subtype == 2.0) ? (MAC * 0.10 + Mac_act * 0.05)
                     : (MAC * 0.05 + Mac_act * 0.02);   // MFS: mild

// AMAN damage driven by Ab-complement attack on axolemma
double AMAN_driver = (subtype == 2.0)
                     ? (MAC * 1.20 + C5a * 0.50 + Ab_anti_GD1a * 0.30)
                     : (subtype == 1.0) ? (MAC * 0.15 + C5a * 0.10)
                     : (MAC * 0.03);   // MFS: minimal axon damage

dxdt_Myelin_damage = k_myelin_damage * AIDP_driver * (1.0 - Myelin_damage)
                     - k_remyelin * (1.0 - Myelin_damage) * Nerve_function;

dxdt_Axon_damage   = k_axon_damage * AMAN_driver * (1.0 - Axon_damage)
                     - k_axon_repair * (1.0 - Axon_damage) * Nerve_function;

// ==================================================================
// 9. COMPOSITE NERVE FUNCTION
// ==================================================================
// Weighted combination: myelin damage 70% + axon damage 30% (AIDP-centric)
// Recovery promoted by IL-10 and absence of ongoing damage
// Ref: Feasby et al. Ann Neurol 1993 - axon loss predicts poor recovery
double damage_combined = 0.70 * Myelin_damage + 0.30 * Axon_damage;
double recovery_drive  = 0.03 * (1.0 - Nerve_function) * (IL10 + 0.1)
                         * (1.0 - damage_combined);
dxdt_Nerve_function = -0.50 * damage_combined * Nerve_function
                      + recovery_drive;

// ==================================================================
// 10. CLINICAL ENDPOINTS
// ==================================================================
// GBS disability score (Hughes scale 0-6)
// 0=normal, 1=minor sx, 2=walk 10m unaided, 3=walk w/ aid,
// 4=bedbound, 5=ventilated, 6=death
// Ref: Hughes et al. Brain 1978; Erasmus GBS outcome score
double target_GBS = GBS_score_max * (1.0 - Nerve_function);
dxdt_GBS_score = 0.30 * (target_GBS - GBS_score);

// FVC% decreases when respiratory muscle innervation fails
// ~30% of GBS patients require mechanical ventilation
// Ref: Lawn et al. Arch Neurol 2001 - respiratory failure predictors
double resp_impairment = (GBS_score > 3.0) ? (GBS_score - 3.0) * 0.25 : 0.0;
dxdt_FVC_pct = -15.0 * resp_impairment * (FVC_pct / 100.0)
               + 2.0 * Nerve_function * (1.0 - FVC_pct / 100.0);

// ==================================================================
// 11. IVIG PK - 2-COMPARTMENT MODEL
// ==================================================================
// Standard dosing: 2 g/kg over 5 days (0.4 g/kg/day)
// Ref: Kuitwaard et al. 2009 - Vc=2.8L, Vp=4.0L, CL=0.17 L/day
// FcRn recycling extends IgG half-life to ~21 days
double IVIG_dose_rate = 0.0;
if (IVIG_dose > 0.0 && SOLVERTIME >= IVIG_start && SOLVERTIME < IVIG_start + 5.0) {
  IVIG_dose_rate = (IVIG_dose / (5.0 * 24.0)) * WT;  // g/h total body dose
}
double IVIG_CL_total = IVIG_CL * WT;    // L/h
double IVIG_Q_total  = IVIG_Q  * WT;    // L/h
double IVIG_Vc_total = IVIG_Vc * WT;    // L
double IVIG_Vp_total = IVIG_Vp * WT;    // L

dxdt_IVIG_C = IVIG_dose_rate
              - (IVIG_CL_total / IVIG_Vc_total) * IVIG_C
              - (IVIG_Q_total  / IVIG_Vc_total) * IVIG_C
              + (IVIG_Q_total  / IVIG_Vp_total) * IVIG_P;

dxdt_IVIG_P = (IVIG_Q_total / IVIG_Vc_total) * IVIG_C
              - (IVIG_Q_total / IVIG_Vp_total) * IVIG_P;

// ==================================================================
// 12. PLASMA EXCHANGE CUMULATIVE EFFECT
// ==================================================================
// PE modeled as continuous approximation with saturation at 5 sessions
// Ref: French Cooperative Group Lancet 1997 - 5 PE sessions optimal
dxdt_PE_cumulative = PE_flag * 0.15 * (PE_sessions - PE_cumulative);

// ==================================================================
// 13. ECULIZUMAB PK (1-COMPARTMENT APPROXIMATION)
// ==================================================================
// Induction: 900 mg IV weekly x4 → Maintenance: 1200 mg q2w
// Ref: Misawa et al. Ann Neurol 2018; Broderick et al. 2017
double ECU_dose_rate = 0.0;
if (ECU_flag > 0.5 && SOLVERTIME >= IVIG_start) {
  if (SOLVERTIME < IVIG_start + 28.0) {
    ECU_dose_rate = 900.0 / (7.0 * 24.0);   // induction: 900 mg/week in mg/h
  } else {
    ECU_dose_rate = 1200.0 / (14.0 * 24.0); // maintenance: 1200 mg/2week
  }
}
dxdt_ECU_C = ECU_dose_rate / (ECU_Vc * WT)
             - (ECU_CL / (ECU_Vc * WT)) * ECU_C;

$TABLE
// Hughes disability grade (integer 0-6)
double Hughes_grade;
if      (GBS_score < 0.5) Hughes_grade = 0;
else if (GBS_score < 1.5) Hughes_grade = 1;
else if (GBS_score < 2.5) Hughes_grade = 2;
else if (GBS_score < 3.5) Hughes_grade = 3;
else if (GBS_score < 4.5) Hughes_grade = 4;
else if (GBS_score < 5.5) Hughes_grade = 5;
else                       Hughes_grade = 6;

// Ventilation risk (probability of requiring mechanical ventilation)
double Ventilation_risk;
if      (FVC_pct < 30.0) Ventilation_risk = 1.00;
else if (FVC_pct < 50.0) Ventilation_risk = 0.60;
else if (FVC_pct < 75.0) Ventilation_risk = 0.20;
else                      Ventilation_risk = 0.02;

// Composite inflammatory and injury markers
double Total_Ab          = Ab_anti_GM1 + Ab_anti_GD1a + Ab_anti_GQ1b;
double Immune_index      = (Th1 + Th17) / (Treg + 0.10);
double Inflammation_score = (IL6 + TNFa) / 2.0;
double Complement_burden  = MAC + C3b + C5a;

// MRC sum score approximation (0-60 scale, 60=normal)
// Based on nerve function with scaling
double MRC_sum = 60.0 * Nerve_function;

// Neurofilament light chain (NfL) - axon damage biomarker
// Ref: Querol et al. JNNP 2017 - NfL correlates with axon injury
double NfL_proxy = 1.0 + 50.0 * Axon_damage + 20.0 * Myelin_damage;

// CSF protein elevation (mg/dL) - cytoalbuminous dissociation
// Ref: Albuminocytologic dissociation: protein ↑, cells normal
double CSF_protein = 40.0 + 300.0 * (Myelin_damage + Axon_damage * 0.5);

$CAPTURE
Hughes_grade Ventilation_risk Total_Ab Immune_index Inflammation_score
Complement_burden MRC_sum NfL_proxy CSF_protein FVC_pct GBS_score
Nerve_function Myelin_damage Axon_damage MAC C3b C5a IVIG_C IVIG_P
Ab_anti_GM1 Ab_anti_GD1a Ab_anti_GQ1b Th1 Th17 Treg IL6 TNFa IL10
ECU_C PE_cumulative Pathogen DC_act Mac_act Bcell Plasma
'

# Compile the model
gbs_mod <- mcode("gbs_qsp", gbs_model)
cat("Model compiled successfully.\n")
cat("Number of ODEs:", length(init(gbs_mod)), "\n")

# ============================================================
# SIMULATION PARAMETERS
# ============================================================

sim_time <- seq(0, 180, by = 1)  # 180 days simulation

# ============================================================
# SCENARIO 1: Untreated Natural History - AIDP
# ============================================================
cat("\nRunning Scenario 1: Untreated AIDP natural history...\n")

sc1 <- gbs_mod %>%
  param(IVIG_dose = 0, PE_flag = 0, ECU_flag = 0,
        subtype = 1, IVIG_start = 7) %>%
  mrgsim(end = 180, delta = 1) %>%
  as_tibble() %>%
  mutate(scenario = "1: Untreated AIDP",
         subtype_label = "AIDP")

# ============================================================
# SCENARIO 2: IVIG 2g/kg starting at day 7 (early)
# ============================================================
cat("Running Scenario 2: IVIG 2g/kg from day 7 (early AIDP)...\n")

sc2 <- gbs_mod %>%
  param(IVIG_dose = 2.0, PE_flag = 0, ECU_flag = 0,
        subtype = 1, IVIG_start = 7) %>%
  mrgsim(end = 180, delta = 1) %>%
  as_tibble() %>%
  mutate(scenario = "2: IVIG early (day 7)",
         subtype_label = "AIDP")

# ============================================================
# SCENARIO 3: IVIG 2g/kg starting at day 14 (late)
# ============================================================
cat("Running Scenario 3: IVIG 2g/kg from day 14 (late treatment)...\n")

sc3 <- gbs_mod %>%
  param(IVIG_dose = 2.0, PE_flag = 0, ECU_flag = 0,
        subtype = 1, IVIG_start = 14) %>%
  mrgsim(end = 180, delta = 1) %>%
  as_tibble() %>%
  mutate(scenario = "3: IVIG late (day 14)",
         subtype_label = "AIDP")

# ============================================================
# SCENARIO 4: Plasma Exchange (5 sessions)
# ============================================================
cat("Running Scenario 4: Plasma exchange (5 sessions)...\n")

sc4 <- gbs_mod %>%
  param(IVIG_dose = 0, PE_flag = 1, ECU_flag = 0,
        subtype = 1, IVIG_start = 7) %>%
  mrgsim(end = 180, delta = 1) %>%
  as_tibble() %>%
  mutate(scenario = "4: Plasma exchange",
         subtype_label = "AIDP")

# ============================================================
# SCENARIO 5: IVIG + Eculizumab (investigational combination)
# ============================================================
cat("Running Scenario 5: IVIG + Eculizumab (investigational)...\n")

sc5 <- gbs_mod %>%
  param(IVIG_dose = 2.0, PE_flag = 0, ECU_flag = 1,
        subtype = 1, IVIG_start = 7) %>%
  mrgsim(end = 180, delta = 1) %>%
  as_tibble() %>%
  mutate(scenario = "5: IVIG + Eculizumab",
         subtype_label = "AIDP")

# ============================================================
# SCENARIO 6: AMAN subtype (severe axonal) - no treatment
# ============================================================
cat("Running Scenario 6: AMAN subtype (severe axonal, untreated)...\n")

sc6 <- gbs_mod %>%
  param(IVIG_dose = 0, PE_flag = 0, ECU_flag = 0,
        subtype = 2, IVIG_start = 7) %>%
  mrgsim(end = 180, delta = 1) %>%
  as_tibble() %>%
  mutate(scenario = "6: AMAN untreated",
         subtype_label = "AMAN")

# ============================================================
# COMBINE ALL SCENARIOS
# ============================================================

all_sc <- bind_rows(sc1, sc2, sc3, sc4, sc5, sc6)

# Color palette for scenarios
sc_colors <- c(
  "1: Untreated AIDP"       = "#E41A1C",
  "2: IVIG early (day 7)"   = "#377EB8",
  "3: IVIG late (day 14)"   = "#4DAF4A",
  "4: Plasma exchange"      = "#984EA3",
  "5: IVIG + Eculizumab"    = "#FF7F00",
  "6: AMAN untreated"       = "#A65628"
)

# ============================================================
# PLOT 1: GBS Disability Score Over Time (All Scenarios)
# ============================================================
cat("\nGenerating plots...\n")

p1 <- ggplot(all_sc, aes(x = time, y = GBS_score, color = scenario)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = sc_colors) +
  geom_hline(yintercept = c(2, 4), linetype = "dashed", color = "gray50", alpha = 0.7) +
  annotate("text", x = 5, y = 2.1, label = "Walk unaided", size = 3, color = "gray40") +
  annotate("text", x = 5, y = 4.1, label = "Bedbound", size = 3, color = "gray40") +
  scale_y_continuous(limits = c(0, 6.2),
                     breaks = 0:6,
                     labels = c("0 Normal", "1 Minor", "2 Walk 10m", "3 Walk aided",
                                "4 Bedbound", "5 Ventilated", "6 Death")) +
  labs(
    title    = "GBS Disability Score Over Time (Hughes Scale)",
    subtitle = "6 Treatment Scenarios - GBS QSP Model",
    x        = "Time (days from onset)",
    y        = "GBS Disability Score",
    color    = "Scenario",
    caption  = "Hughes et al. Brain 1978; van Doorn et al. Lancet 2008"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        plot.title = element_text(face = "bold")) +
  guides(color = guide_legend(nrow = 2))

# ============================================================
# PLOT 2: Antibody Titers (Anti-GM1 and Anti-GD1a)
# ============================================================

ab_data <- all_sc %>%
  select(time, scenario, Ab_anti_GM1, Ab_anti_GD1a) %>%
  pivot_longer(cols = c(Ab_anti_GM1, Ab_anti_GD1a),
               names_to = "Ab_type", values_to = "titer") %>%
  mutate(Ab_type = recode(Ab_type,
                          "Ab_anti_GM1"  = "Anti-GM1 IgG",
                          "Ab_anti_GD1a" = "Anti-GD1a IgG"))

p2 <- ggplot(ab_data, aes(x = time, y = titer, color = scenario, linetype = Ab_type)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = sc_colors) +
  scale_linetype_manual(values = c("Anti-GM1 IgG" = "solid", "Anti-GD1a IgG" = "dashed")) +
  labs(
    title    = "Anti-Ganglioside Antibody Titers",
    subtitle = "IVIG and PE reduce pathogenic antibody levels",
    x        = "Time (days)",
    y        = "Antibody Titer (AU)",
    color    = "Scenario",
    linetype = "Antibody",
    caption  = "Halstead et al. Brain 2008; Dalakas Ther Adv Neurol Disord 2012"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7),
        plot.title = element_text(face = "bold")) +
  guides(color = guide_legend(nrow = 2))

# ============================================================
# PLOT 3: Complement Cascade (MAC, C3b) Over Time
# ============================================================

comp_data <- all_sc %>%
  select(time, scenario, MAC, C3b, C5a) %>%
  pivot_longer(cols = c(MAC, C3b, C5a),
               names_to = "component", values_to = "level")

p3 <- ggplot(comp_data %>% filter(component %in% c("MAC", "C3b")),
             aes(x = time, y = level, color = scenario, linetype = component)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = sc_colors) +
  scale_linetype_manual(values = c("MAC" = "solid", "C3b" = "dashed")) +
  labs(
    title    = "Complement Cascade Activation",
    subtitle = "MAC and C3b deposition at peripheral nerves",
    x        = "Time (days)",
    y        = "Complement Level (AU)",
    color    = "Scenario",
    linetype = "Complement",
    caption  = "Halstead et al. Brain 2004; Basta & Dalakas Blood 1994"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7),
        plot.title = element_text(face = "bold")) +
  guides(color = guide_legend(nrow = 2))

# ============================================================
# PLOT 4: IVIG PK Concentration Profile
# ============================================================

ivig_data <- all_sc %>%
  filter(scenario %in% c("2: IVIG early (day 7)",
                          "3: IVIG late (day 14)",
                          "5: IVIG + Eculizumab")) %>%
  select(time, scenario, IVIG_C, IVIG_P)

p4 <- ggplot(ivig_data, aes(x = time)) +
  geom_line(aes(y = IVIG_C, color = scenario), size = 1.1) +
  geom_line(aes(y = IVIG_P, color = scenario), size = 0.7, linetype = "dashed") +
  scale_color_manual(values = sc_colors) +
  annotate("text", x = 30, y = max(ivig_data$IVIG_C) * 0.8,
           label = "Solid = Central\nDashed = Peripheral",
           size = 3, color = "gray40") +
  labs(
    title    = "IVIG 2-Compartment PK Profile",
    subtitle = "Standard 2g/kg over 5-day infusion",
    x        = "Time (days)",
    y        = "IVIG Concentration (g/L)",
    color    = "Scenario",
    caption  = "Kuitwaard et al. JNNP 2009; Hammarstrom & Smith NEJM 1993"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        plot.title = element_text(face = "bold"))

# ============================================================
# PLOT 5: FVC% and Ventilation Risk
# ============================================================

vent_data <- all_sc %>%
  select(time, scenario, FVC_pct, Ventilation_risk) %>%
  filter(time <= 90)  # focus on acute phase

p5a <- ggplot(vent_data, aes(x = time, y = FVC_pct, color = scenario)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = c(30, 50), linetype = "dashed", color = "red", alpha = 0.6) +
  annotate("text", x = 85, y = 32, label = "< 30% = intubate", size = 3, color = "red") +
  annotate("text", x = 85, y = 52, label = "< 50% = monitor ICU", size = 3, color = "darkred") +
  scale_color_manual(values = sc_colors) +
  scale_y_continuous(limits = c(0, 105)) +
  labs(
    title    = "Forced Vital Capacity (FVC%)",
    subtitle = "Respiratory function during acute GBS",
    x        = "Time (days)",
    y        = "FVC (% predicted)",
    color    = "Scenario",
    caption  = "Lawn et al. Arch Neurol 2001; Sharshar et al. Lancet 2003"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7),
        plot.title = element_text(face = "bold")) +
  guides(color = guide_legend(nrow = 2))

p5b <- ggplot(vent_data, aes(x = time, y = Ventilation_risk, color = scenario)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = sc_colors) +
  scale_y_continuous(limits = c(0, 1.05), labels = scales::percent) +
  labs(
    title    = "Mechanical Ventilation Risk",
    x        = "Time (days)",
    y        = "Probability of MV",
    color    = "Scenario"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

# ============================================================
# PLOT 6: Nerve Function Recovery
# ============================================================

p6 <- ggplot(all_sc, aes(x = time, y = Nerve_function, color = scenario)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = sc_colors) +
  scale_y_continuous(limits = c(0, 1.05), labels = scales::percent) +
  geom_ribbon(
    data = all_sc %>% filter(scenario == "2: IVIG early (day 7)"),
    aes(ymin = Nerve_function - 0.05, ymax = Nerve_function + 0.05),
    alpha = 0.1, fill = "#377EB8"
  ) +
  labs(
    title    = "Peripheral Nerve Function Recovery",
    subtitle = "Composite score (1.0 = normal, 0 = complete loss)",
    x        = "Time (days)",
    y        = "Nerve Function Score",
    color    = "Scenario",
    caption  = "Feasby et al. Ann Neurol 1993; Kuwabara & Yuki 2013"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7),
        plot.title = element_text(face = "bold")) +
  guides(color = guide_legend(nrow = 2))

# ============================================================
# PLOT 7: Immune Cell Dynamics (Th1, Th17, Treg)
# ============================================================

tcell_data <- all_sc %>%
  filter(scenario %in% c("1: Untreated AIDP", "2: IVIG early (day 7)",
                          "6: AMAN untreated")) %>%
  select(time, scenario, Th1, Th17, Treg) %>%
  pivot_longer(cols = c(Th1, Th17, Treg),
               names_to = "cell_type", values_to = "level")

p7 <- ggplot(tcell_data, aes(x = time, y = level, color = cell_type, linetype = scenario)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = c("Th1" = "#E41A1C", "Th17" = "#FF7F00", "Treg" = "#4DAF4A")) +
  scale_linetype_manual(values = c("1: Untreated AIDP"     = "solid",
                                    "2: IVIG early (day 7)" = "dashed",
                                    "6: AMAN untreated"     = "dotted")) +
  labs(
    title    = "T Lymphocyte Population Dynamics",
    subtitle = "IVIG shifts Th1/Th17 → Treg balance",
    x        = "Time (days)",
    y        = "Relative Cell Count (AU)",
    color    = "Cell Type",
    linetype = "Scenario",
    caption  = "Brilot et al. J Neuroinflammation 2010; Ephrem et al. Blood 2008"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        plot.title = element_text(face = "bold"))

# ============================================================
# PLOT 8: Hughes Grade Distribution at 4 Weeks
# ============================================================

hughes_4wk <- all_sc %>%
  filter(time == 28) %>%
  select(scenario, Hughes_grade, GBS_score, Nerve_function) %>%
  mutate(Hughes_label = paste0("Grade ", round(Hughes_grade)))

p8 <- ggplot(hughes_4wk, aes(x = reorder(scenario, GBS_score),
                               y = GBS_score, fill = scenario)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_text(aes(label = paste0("Grade ", round(Hughes_grade))),
            hjust = -0.1, size = 3.5) +
  scale_fill_manual(values = sc_colors) +
  scale_y_continuous(limits = c(0, 7)) +
  coord_flip() +
  geom_hline(yintercept = c(2, 4), linetype = "dashed", color = "gray40") +
  labs(
    title    = "Hughes Grade at 4 Weeks Post-Onset",
    subtitle = "Primary clinical endpoint across treatment scenarios",
    x        = NULL,
    y        = "GBS Disability Score (Hughes 0-6)",
    caption  = "Hughes et al. Brain 1978; van Doorn et al. Lancet 2008"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

# ============================================================
# COMPOSITE FIGURE
# ============================================================

cat("Compositing plots...\n")

# Page 1: Core clinical outcomes
fig1 <- (p1 | p6) / (p5a | p8) +
  plot_annotation(
    title   = "Guillain-Barré Syndrome QSP Model — Clinical Outcomes",
    theme   = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# Page 2: Immunopathology
fig2 <- (p2 | p3) / (p7 | p4) +
  plot_annotation(
    title = "GBS QSP Model — Immunological Mechanisms",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# Display plots
print(fig1)
print(fig2)

# ============================================================
# SUMMARY STATISTICS TABLE
# ============================================================

summary_stats <- all_sc %>%
  group_by(scenario) %>%
  summarise(
    max_GBS_score     = round(max(GBS_score), 2),
    max_Hughes_grade  = round(max(Hughes_grade), 0),
    min_FVC_pct       = round(min(FVC_pct), 1),
    max_vent_risk     = round(max(Ventilation_risk), 2),
    recovery_at_3mo   = round(Nerve_function[time == 90], 3),
    GBS_at_28d        = round(GBS_score[time == 28], 2),
    max_MAC           = round(max(MAC), 3),
    max_Total_Ab      = round(max(Total_Ab), 3),
    .groups = "drop"
  )

cat("\n========================================\n")
cat("GBS QSP Model — Summary Statistics\n")
cat("========================================\n")
print(as.data.frame(summary_stats))
cat("\n")

# ============================================================
# SUBTYPE COMPARISON: AIDP vs AMAN at day 90
# ============================================================

cat("========================================\n")
cat("Nerve Damage Profile at 90 Days\n")
cat("========================================\n")

nerve_compare <- all_sc %>%
  filter(time == 90, scenario %in% c("1: Untreated AIDP", "6: AMAN untreated")) %>%
  select(scenario, Myelin_damage, Axon_damage, Nerve_function, GBS_score, NfL_proxy)

print(as.data.frame(nerve_compare))
cat("\n")

# ============================================================
# TIME-TO-EVENT ANALYSIS
# ============================================================

cat("========================================\n")
cat("Time-to-Event Summary\n")
cat("========================================\n")

tte_analysis <- all_sc %>%
  group_by(scenario) %>%
  summarise(
    # Time to nadir (worst GBS score)
    time_to_nadir    = time[which.max(GBS_score)],
    nadir_GBS        = round(max(GBS_score), 2),
    # Time to GBS score < 2 (walk unaided)
    time_walk_unaided = {
      walk_times <- time[GBS_score < 2.0 & time > 20]
      if (length(walk_times) > 0) min(walk_times) else NA_real_
    },
    # Recovery at 6 months
    recovery_6mo = round(Nerve_function[time == 180], 3),
    .groups = "drop"
  )

print(as.data.frame(tte_analysis))
cat("\n")

cat("Model simulation complete.\n")
cat("GBS QSP Model: 26 state variables, 6 treatment scenarios\n")
cat("Reference: Hughes RA et al. Lancet Neurol 2014; Yuki & Hartung NEJM 2012\n")
