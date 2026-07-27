# =============================================================================
#  Clostridioides difficile Infection (CDI) — QSP model for mrgsolve
# =============================================================================
#
#  WHAT THIS MODEL IS FOR
#  ----------------------
#  CDI is unusual among infections in that killing the organism is the easy
#  part.  Oral vancomycin sterilises the stool of vegetative C. difficile
#  within days and cures ~80% of episodes — and then one in four patients
#  relapses, because the same drug that cleared the organism also kept the
#  microbial guild that provides colonization resistance flat on the floor
#  while a drug-proof spore reservoir sat waiting.  Every therapeutic advance
#  of the last 15 years (fidaxomicin, bezlotoxumab, FMT, SER-109, RBX2660)
#  is an attempt to act on that gap rather than on the kill rate.
#
#  So this model is built around the recurrence loop rather than around the
#  bug.  It carries three things simultaneously:
#
#    (1) an ECOLOGICAL layer   — six microbial guilds under drug-specific
#        kill, and the two currencies through which they suppress
#        C. difficile: secondary bile acids and the luminal nutrient niche;
#    (2) a PATHOGEN layer      — spore -> germination -> vegetative outgrowth
#        -> sporulation, with a mucosa-associated spore reservoir that is
#        the actual seed of relapse;
#    (3) a HOST layer          — PaLoc-regulated toxin output, Rho-GTPase
#        driven barrier failure, Wnt/FZD blockade of epithelial renewal,
#        neutrophilic inflammation, and the clinical read-outs that trials
#        actually measure (unformed stools/day, WBC, creatinine, albumin).
#
#  The therapeutic point of the structure: KILLCD (how fast the organism
#  dies) and KILLSBA (how hard the restoring guild is hit) are separate
#  quantities.  Agents that score well on the first and badly on the second
#  cure the episode and buy the relapse.
#
#  STRUCTURE
#  ---------
#    61 ODE compartments:
#      microbiota guilds ........  6   (MB_SBA MB_BUT MB_BAC MB_BIF MB_ENT MB_ENC)
#      bile acids ...............  5   (BA_TCA BA_CA BA_CDCA BA_DCA BA_LCA)
#      luminal nutrients / SCFA .  3   (NUT_SIA NUT_AA SCFA_BUT)
#      C. difficile .............  4   (CD_SPORE_L CD_VEG CD_MUC CD_SPORE_B)
#      toxins ...................  6   (TCDA TCDB TOX_CPLX TCDA_MUC TCDB_MUC CDT)
#      epithelium / barrier .....  5   (EPI EPI_SC EPI_TJ EPI_MUCUS EPI_PERM)
#      immunity .................  7   (IM_IL8 IM_IL1B IM_TNF IM_NEUT IM_IL22 IM_PSM AB_IGG)
#      clinical .................  5   (EPI_H2O STOOL WBC ALB CRE)
#      drug PK .................. 20   (vancomycin 3, fidaxomicin 4, metronidazole 3,
#                                       bezlotoxumab 3, rifaximin 2, ridinilazole 2,
#                                       index antibiotic 2, live biotherapeutic 1)
#
#  SELF-CALIBRATING HEALTHY BASELINE
#  ---------------------------------
#  Rather than hard-coding initial conditions and hoping they are consistent
#  with the rate constants, $MAIN solves the healthy colonic steady state
#  from the parameters themselves.  The bile-acid initial conditions come out
#  of an exact algebraic solution of the three-step microbial cascade
#  (conjugated pool -> BSH deconjugation -> bai 7alpha-dehydroxylation),
#  including the closed-form root of the Michaelis-Menten balance for
#  cholate and chenodeoxycholate.  Change k7A or the BSH parameters and the
#  baseline moves with them instead of silently drifting.
#
#  UNITS
#  -----
#    time            days
#    guilds          relative abundance of total colonic community (0-1)
#    bile acids      micromol/L of colonic water
#    NUT_SIA/NUT_AA  mmol/L        SCFA_BUT  mmol/L
#    C. difficile    10^6 CFU/g faeces   (so CD_VEG = 100  <=>  10^8 CFU/g)
#    toxins          ng/mL faecal supernatant (mucosal compartments: ng/mL-equiv)
#    EPI/TJ/MUCUS/SC fraction of normal (0-1)
#    cytokines       arbitrary normalised units (0 at health)
#    drug amounts    mg  (colonic concentrations reported as ug/g faeces)
#    bezlotoxumab    mg central/peripheral;  BEZ_GUT in mg/L (= ug/mL) of lumen
#    STOOL           unformed stools / 24 h        WBC  10^9/L
#    ALB             g/dL                          CRE  mg/dL
#
#  CALIBRATION ANCHORS (see cdi_references.md for full citations)
#  --------------------------------------------------------------
#    * Louie 2011 NEJM / Cornely 2012 Lancet ID — fidaxomicin vs vancomycin:
#      clinical cure ~88% both arms, recurrence 15.4% vs 25.3%.
#    * Johnson 2014 Clin Infect Dis (POLYMER pooled) — vancomycin cure 81.1%
#      vs metronidazole 72.7%.
#    * Guery 2018 Lancet ID (EXTEND) — extended-pulsed fidaxomicin,
#      sustained cure day 90 70% vs 59%.
#    * Wilcox 2017 NEJM (MODIFY I/II) — bezlotoxumab recurrence 16.5% vs 26.6%.
#    * van Nood 2013 NEJM — donor faeces infusion cure 81% (94% incl. retreat).
#    * Feuerstadt 2022 NEJM (ECOSPOR III) — SER-109 recurrence 12% vs 40%.
#    * Khanna 2022 Drugs (PUNCH CD3) — RBX2660 sustained response 70.6% vs 57.5%.
#    * Kyne 2001 Lancet — low day-3 serum anti-toxin A IgG predicts recurrence.
#    * Theriot 2014 Nat Commun / Buffie 2015 Nature — antibiotic loss of
#      secondary bile acids and of bai+ Clostridium scindens.
#    * Sorg & Sonenshein 2008/2009 J Bacteriol — taurocholate germination,
#      chenodeoxycholate competitive inhibition.
#    * Ng 2013 Nature — microbiota-liberated sialic acid fuels expansion.
#    * Louie 2009 AAC / Sears 2012 CID — faecal fidaxomicin and OP-1118.
#    * Bolton 1986 Gut — faecal metronidazole falls as diarrhoea resolves.
#
#  DISCLAIMER
#  ----------
#  Educational / research QSP model.  Semi-quantitative: parameters are
#  literature-anchored order-of-magnitude estimates, not a fitted population
#  model.  Not for clinical or regulatory use.
#
#  Usage:
#    Rscript cdi_mrgsolve_model.R          # builds, runs 16 scenarios, prints report
#    mod <- mread("cdi", "cdi_mrgsolve_model.R")   # if extracted to .cpp
# =============================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
})

cdi_code <- '
$PROB
# Clostridioides difficile infection — ecology / pathogen / host QSP model

$PARAM @annotated
// ---------------- patient & strain covariates ----------------------------
WT      :  70   : Body weight (kg)
RT027   :   0   : Hypervirulent ribotype 027/NAP1 flag (0/1)
VRE     :   0   : Vancomycin-resistant Enterococcus flag (0/1)
IGG0    :   1.0 : Baseline anti-TcdB IgG (1 = population median)
IMMCOMP :   1.0 : Humoral response capacity (0.2 immunosuppressed .. 1.5 robust)
FIBER   :   1.0 : Fermentable fibre intake (1 = normal, 0.3 = tube feed)
PPI     :   0   : Proton-pump inhibitor flag (raises effective spore inoculum)
TIGE    :   0   : IV tigecycline salvage flag (0/1)

// ---------------- microbiota guilds --------------------------------------
SBA0    : 0.10  : Healthy abundance, bai+ 7a-dehydroxylating Clostridia
BUT0    : 0.30  : Healthy abundance, butyrogenic Lachnospiraceae/Ruminococcaceae
BAC0    : 0.42  : Healthy abundance, Bacteroidetes
BIF0    : 0.08  : Healthy abundance, Bifidobacterium/Actinobacteria
ENT0    : 0.004 : Healthy abundance, Enterobacteriaceae
ENC0    : 0.002 : Healthy abundance, Enterococcus
MUSBA   : 0.20  : Community-limited recovery rate, MB_SBA (1/d)
MUBUT   : 0.25  : Community-limited recovery rate, MB_BUT (1/d)
MUBAC   : 0.55  : Community-limited recovery rate, MB_BAC (1/d)
MUBIF   : 0.30  : Community-limited recovery rate, MB_BIF (1/d)
MUENT   : 2.20  : Intrinsic growth rate, MB_ENT (1/d)
MUENC   : 1.60  : Intrinsic growth rate, MB_ENC (1/d)
EPSMB   : 1e-7  : Relative immigration floor per guild (1/d)
TRECOL  : 20.0  : Re-establishment time constant of a depleted guild (d)
KRECK   : 0.50  : Guild-specific kill rate halving that guilds re-seeding
KENTMAX : 0.28  : Max Enterobacteriaceae abundance in an emptied niche
KENCMAX : 0.20  : Max Enterococcus abundance in an emptied niche
HDYS    : 1.30  : Hill exponent on niche emptiness for opportunist bloom
KI22ENT : 0.55  : IL-22 level halving the opportunist carrying capacity

// ---------------- bile acid metabolism ----------------------------------
TCAIN   : 727   : Colonic delivery of conjugated primary bile acids (uM/d)
KDEC    : 60    : Max BSH deconjugation rate constant (1/d)
BSHMIN  : 0.08  : Residual BSH activity when the BSH guilds are gone
WENCBSH : 3.0   : Per-abundance BSH weight of blooming Enterococcus/Lactobacillus
FCA     : 0.60  : Fraction of deconjugated pool that is cholate
K7A     : 10590 : bai 7a-dehydroxylation Vmax for cholate (uM/d per abundance)
K7C     : 9640  : bai 7a-dehydroxylation Vmax for CDCA (uM/d per abundance)
KM7     : 100   : Michaelis constant of 7a-dehydroxylation (uM)
KABST   : 0.60  : Loss rate of conjugated pool (1/d)
KABSCA  : 1.40  : Loss rate of cholate (1/d)
KABSCD  : 1.086 : Loss rate of chenodeoxycholate (1/d)
KABSDCA : 0.80  : Loss rate of deoxycholate (1/d)
KABSLCA : 0.833 : Loss rate of lithocholate (1/d)
UDCA    : 0     : Ursodiol-derived luminal UDCA (uM, germination competitor)

// ---------------- luminal nutrient niche --------------------------------
JSIA0   : 0.20  : Host-derived free sialic acid input (mM/d)
JSIAM   : 1.00  : Microbial sialidase-liberated sialic acid input (mM/d)
FSIAH   : 0.35  : Fraction of sialidase input independent of Bacteroidetes
KCSIA   : 23.5  : Fermenter consumption of sialic acid (1/d at healthy niche)
KCSIACD : 0.05  : C. difficile sialic acid consumption (1/d per 10^6 CFU/g)
KOUTSIA : 0.50  : Washout of free sialic acid (1/d)
JAA     : 6.00  : Stickland amino-acid input (mM/d)
KCAA    : 5.50  : Fermenter consumption of Stickland amino acids (1/d, healthy)
KCAACD  : 0.06  : C. difficile amino-acid consumption (1/d per 10^6 CFU/g)
KOUTAA  : 0.50  : Washout of Stickland amino acids (1/d)
KBUTP   : 45    : Butyrate production (mM/d at healthy MB_BUT and fibre)
KABSBUT : 3.00  : Butyrate absorption/consumption (1/d)

// ---------------- C. difficile life cycle --------------------------------
INOC    : 0.05  : Ingested spore inoculum (10^6 CFU/g equivalent)
KGERM   : 6.0   : Maximum germination rate constant (1/d)
KTCA    : 150   : Taurocholate germination K50 (uM)
KICDCA  : 400   : CDCA competitive inhibition constant for germination (uM)
KIUDCA  : 500   : UDCA competitive inhibition constant for germination (uM)
KI2GERM : 250   : Secondary bile-acid inhibition of germination (uM)
MUCDMAX : 14.0  : Maximum C. difficile growth rate (1/d)
KSIA    : 0.35  : Sialic acid K50 for growth (mM)
KAA     : 3.00  : Stickland amino-acid K50 for growth (mM)
KI2GROW : 120   : Secondary bile-acid inhibition of outgrowth (uM)
KIBUTCD : 25    : Butyrate inhibition of outgrowth (mM)
KCDMAX  : 200   : Maximum luminal carrying capacity (10^6 CFU/g)
LAMNICHE: 5.50  : Steepness of niche-based colonization resistance
WNICHEF : 0.65  : Weight of the fermenter guilds in colonization resistance
KNSUP   : 1.20  : Neutrophil level halving commensal regrowth
KSSUP   : 6.00  : Excess stools/day halving commensal regrowth
KDCD    : 0.35  : Baseline vegetative death rate (1/d)
KOUTCD  : 0.50  : Vegetative washout with colonic transit (1/d)
KOUTSP  : 0.90  : Spore washout with colonic transit (1/d)
KADH    : 0.55  : Adhesion to mucus/mucosa (1/d)
KDET    : 0.45  : Detachment from mucosa (1/d)
FMUCADV : 1.60  : Growth advantage of the mucus-associated niche
KDMUC   : 0.60  : Death rate of the mucosa-associated population (1/d)
KMUCMAX : 12    : Carrying capacity of the mucosal niche (10^6 CFU/g)
KSPOR   : 0.25  : Baseline sporulation rate (1/d)
RTSPOR  : 0.80  : Extra sporulation in ribotype 027
KI2SPOR : 300   : Secondary bile-acid suppression of sporulation (uM)
FSPORB  : 0.45  : Fraction of sporulation entering the mucosal reservoir
KGERMB  : 0.55  : Re-germination rate of the mucosal reservoir (1/d)
KCLRB   : 0.20  : Immune/mechanical clearance of the reservoir (1/d)

// ---------------- PaLoc regulation and toxins ----------------------------
KTOXB   : 20.0  : TcdB production (ng/mL/d per 10^6 CFU/g at full induction)
KTOXA   : 50.0  : TcdA production (ng/mL/d per 10^6 CFU/g at full induction)
WMUC    : 5.0   : Toxin-delivery weight of mucosa-adherent cells
FTOXB   : 1.0   : Baseline PaLoc induction scalar
KNTOX   : 2.50  : CodY/CcpA nutrient repression constant (mM)
KIBUTOX : 30    : Butyrate damping of tcdA/tcdB (mM)
RTTOXF  : 2.20  : Toxin multiplier for ribotype 027 (truncated tcdC)
IC50TF  : 40    : Fidaxomicin concentration halving toxin output (ug/g)
IC50SF  : 30    : Fidaxomicin concentration halving sporulation (ug/g)
KDEGB   : 2.00  : Luminal TcdB degradation (1/d)
KDEGA   : 2.40  : Luminal TcdA degradation (1/d)
KOUTTOX : 0.80  : Toxin washout (1/d)
KTRB    : 0.40  : TcdB transfer to the mucosal compartment (1/d)
KTRA    : 0.30  : TcdA transfer to the mucosal compartment (1/d)
KINTB   : 3.00  : Internalisation/turnover of mucosal TcdB (1/d)
KINTA   : 3.40  : Internalisation/turnover of mucosal TcdA (1/d)
KIAB    : 3.00  : Anti-toxin IgG level halving mucosal toxin delivery
WBEZ    : 0.90  : Weight of luminal bezlotoxumab in the antibody term
ENTSAT  : 0.05  : Enterobacteriaceae abundance saturating the LPS signal
KONB    : 4.50  : Bezlotoxumab-TcdB capture rate (per (ug/mL) antibody per day)
KOFFB   : 0.05  : Bezlotoxumab-TcdB dissociation (1/d)
KDEGCX  : 1.50  : Clearance of the antibody-toxin complex (1/d)
KCDTP   : 3.0   : CDT binary toxin production (ng/mL/d per 10^6 CFU/g)
KDEGCDT : 2.00  : CDT degradation (1/d)

// ---------------- epithelium and barrier ---------------------------------
WA      : 0.35  : Potency weight of TcdA relative to TcdB
WCDT    : 0.45  : Potency weight of CDT on epithelial injury
KCDT50  : 8.0   : CDT K50 (ng/mL)
KDMG    : 1.20  : Maximum toxin-driven colonocyte loss (1/d)
KDMG50  : 160   : Mucosal toxin load producing half-maximal injury (ng/mL)
KAPOP   : 0.05  : Cytokine-driven colonocyte apoptosis (1/d per unit)
KREG    : 1.20  : Stem-cell-driven epithelial regeneration (1/d)
KSCR    : 0.35  : Stem/progenitor pool recovery (1/d)
KSCL    : 0.60  : Toxin-driven stem/progenitor loss (1/d)
KSC50   : 60    : Mucosal TcdB K50 for stem-cell loss (ng/mL)
KFZD    : 45    : Mucosal TcdB halving Wnt/FZD renewal signalling (ng/mL)
KTJR    : 1.60  : Tight-junction repair (1/d)
KTJT    : 2.60  : Rho-glucosylation-driven junction loss (1/d)
KTJ50   : 76    : Mucosal toxin load K50 for junction loss (ng/mL)
KTJI    : 0.12  : Cytokine-driven junction loss (1/d per unit)
KPON    : 2.20  : Permeability build-up (1/d)
KPOFF   : 2.00  : Permeability resolution (1/d)
KMUCR   : 0.70  : Mucus restoration (1/d)
KMUCL   : 0.55  : Mucus loss from neutrophils and toxin (1/d)

// ---------------- immunity -----------------------------------------------
KIL8    : 1.60  : IL-8/CXCL8 production (units/d)
KI850   : 75    : Toxin load K50 for IL-8 induction (ng/mL)
KDIL8   : 1.60  : IL-8 elimination (1/d)
KIL1    : 3.60  : IL-1beta production via pyrin/NLRP3 (units/d)
KI150   : 74    : Toxin load K50 for IL-1beta induction (ng/mL)
KDIL1   : 1.50  : IL-1beta elimination (1/d)
KTNF    : 1.40  : TNF-alpha/IL-6 production (units/d)
KTN50   : 3.0   : Mucosal C. difficile K50 for TLR signalling (10^6 CFU/g)
KDTNF   : 1.40  : TNF elimination (1/d)
KNREC   : 1.10  : Neutrophil recruitment (1/d per cytokine unit)
KNDEC   : 1.00  : Neutrophil turnover (1/d)
KIL22   : 0.60  : IL-22 production by ILC3/Th17 (units/d)
KDIL22  : 1.20  : IL-22 elimination (1/d)
KPSM    : 0.50  : Pseudomembrane formation (1/d)
KPSMR   : 0.55  : Pseudomembrane resolution (1/d)
KABP    : 0.26  : Anti-toxin IgG production (units/d)
KABA50  : 120   : Luminal toxin K50 for IgG induction (ng/mL)
ABMAX   : 4.5   : Maximum anti-toxin IgG level
KABD    : 0.008 : Anti-toxin IgG decay (1/d)

// ---------------- clinical read-outs ------------------------------------
KH2O    : 1.80  : Luminal fluid accumulation (1/d)
WSEC    : 1.60  : Secretory (CFTR/NKCC1) weight of toxin load
KSEC50  : 75    : Toxin load K50 for secretion (ng/mL)
KH2OA   : 2.20  : Colonic water reabsorption (1/d)
STOOL0  : 0.80  : Baseline unformed stools per day
SMAX    : 14.0  : Maximum unformed stools per day
KSTL    : 1.10  : Luminal fluid producing half-maximal stool frequency
TAUSTL  : 0.60  : Stool-frequency transduction time (d)
KILEUS  : 8.00  : Pseudomembrane/ileus level halving measured stool output
WBC0    : 7.0   : Baseline WBC (10^9/L)
WBCMAX  : 16.0  : Maximum WBC increment (10^9/L)
KWBC    : 2.20  : Cytokine level giving half-maximal leukocytosis
TAUWBC  : 0.50  : WBC transduction time (d)
ALB0    : 4.20  : Baseline serum albumin (g/dL)
KALBD   : 0.055 : Albumin turnover (1/d)
KALBL   : 0.018 : Permeability-driven albumin loss (1/d per permeability unit)
KALBI   : 8.00  : TNF level halving hepatic albumin synthesis
CRE0    : 0.90  : Baseline serum creatinine (mg/dL)
KDEH    : 1.30  : Maximal creatinine rise from volume depletion (fold)
KSTC    : 6.00  : Excess stools/day giving half-maximal creatinine rise
TAUCRE  : 1.00  : Creatinine transduction time (d)

// ---------------- antibacterial PK --------------------------------------
KTRV    : 6.00  : Oral vancomycin gastrointestinal transit (1/d)
KEXV    : 2.20  : Vancomycin faecal elimination (1/d)
VCOLG   : 200   : Mass of colonic content used for faecal concentration (g)
KWASH   : 0.15  : Extra drug washout per excess stool/day
KTRF    : 6.00  : Oral fidaxomicin transit (1/d)
KEXF    : 2.00  : Fidaxomicin faecal elimination (1/d)
KMETF   : 1.00  : Fidaxomicin -> OP-1118 hydrolysis (1/d)
FMET    : 0.90  : Molar yield of OP-1118
KEXOP   : 1.80  : OP-1118 faecal elimination (1/d)
WOP     : 0.50  : Relative potency of OP-1118
KAM     : 12.0  : Metronidazole absorption (1/d)
FMTZ    : 0.9995: Metronidazole oral bioavailability (essentially complete)
KEM     : 2.08  : Metronidazole systemic elimination (1/d)
VMTZ    : 50    : Metronidazole volume of distribution (L)
KSECM   : 0.060 : Metronidazole secretion into inflamed colon (1/d)
KEXM    : 4.00  : Metronidazole faecal elimination (1/d)
KTRR    : 6.00  : Rifaximin transit (1/d)
KEXR    : 1.60  : Rifaximin faecal elimination (1/d)
KTRD    : 6.00  : Ridinilazole transit (1/d)
KEXD    : 2.00  : Ridinilazole faecal elimination (1/d)
KTRX    : 6.00  : Index antibiotic delivery to colon (1/d)
KEXX    : 1.00  : Index antibiotic faecal elimination (1/d)

// ---------------- bezlotoxumab PK ---------------------------------------
CLB     : 0.317 : Bezlotoxumab clearance (L/d)
VCB     : 3.08  : Bezlotoxumab central volume (L)
VPB     : 4.10  : Bezlotoxumab peripheral volume (L)
QB      : 0.43  : Bezlotoxumab intercompartmental clearance (L/d)
KTRAB   : 0.450 : Transudation of IgG into the lumen (1/d per permeability unit)
KDGUT   : 1.20  : Luminal antibody degradation/washout (1/d)

// ---------------- live biotherapeutic / FMT ------------------------------
KENG    : 0.70  : Engraftment rate of the administered consortium (1/d)
KABXENG : 0.35  : Antibiotic pressure halving engraftment
FE_SBA  : 0.00  : Fraction of the dose that is bai+ Clostridia
FE_BUT  : 0.00  : Fraction of the dose that is butyrogenic Firmicutes
FE_BAC  : 0.00  : Fraction of the dose that is Bacteroidetes
FE_BIF  : 0.00  : Fraction of the dose that is Bifidobacterium/Actinobacteria

// ---------------- drug PD (kill and collateral damage) -------------------
EMAXV   : 14.0  : Vancomycin maximal kill of C. difficile (1/d)
EC50V   : 45    : Vancomycin faecal EC50 vs C. difficile (ug/g)
EMAXF   : 13.0  : Fidaxomicin maximal kill of C. difficile (1/d)
EC50F   : 18    : Fidaxomicin faecal EC50 vs C. difficile (ug/g)
EMAXM   : 12.0  : Metronidazole maximal kill of C. difficile (1/d)
EC50M   : 4.00  : Metronidazole faecal EC50 vs C. difficile (ug/g)
EMAXD   : 13.0  : Ridinilazole maximal kill of C. difficile (1/d)
EC50D   : 22    : Ridinilazole faecal EC50 vs C. difficile (ug/g)
EMAXR   : 6.00  : Rifaximin maximal kill of C. difficile (1/d)
EC50R   : 350   : Rifaximin faecal EC50 vs C. difficile (ug/g)
EMAXT   : 1.10  : Tigecycline systemic contribution to kill (1/d)
FMUCPEN : 0.65  : Penetration of luminal drug to the mucosal population
KV_SBA  : 3.40  : Vancomycin kill of MB_SBA (1/d)
KV_BUT  : 2.70  : Vancomycin kill of MB_BUT (1/d)
KV_BAC  : 0.10  : Vancomycin kill of MB_BAC (1/d)
KV_BIF  : 1.95  : Vancomycin kill of MB_BIF (1/d)
KV_ENC  : 2.00  : Vancomycin kill of MB_ENC (1/d, abolished if VRE)
KF_SBA  : 0.52  : Fidaxomicin kill of MB_SBA (1/d)
KF_BUT  : 0.55  : Fidaxomicin kill of MB_BUT (1/d)
KF_BAC  : 0.05  : Fidaxomicin kill of MB_BAC (1/d)
KF_BIF  : 0.37  : Fidaxomicin kill of MB_BIF (1/d)
KF_ENC  : 0.50  : Fidaxomicin kill of MB_ENC (1/d)
KM_SBA  : 1.30  : Metronidazole kill of MB_SBA (1/d)
KM_BUT  : 1.20  : Metronidazole kill of MB_BUT (1/d)
KM_BAC  : 1.50  : Metronidazole kill of MB_BAC (1/d)
KM_BIF  : 0.60  : Metronidazole kill of MB_BIF (1/d)
KD_SBA  : 0.28  : Ridinilazole kill of MB_SBA (1/d)
KD_BUT  : 0.35  : Ridinilazole kill of MB_BUT (1/d)
KD_BAC  : 0.04  : Ridinilazole kill of MB_BAC (1/d)
KD_BIF  : 0.35  : Ridinilazole kill of MB_BIF (1/d)
KD_ENC  : 0.25  : Ridinilazole kill of MB_ENC (1/d)
KR_SBA  : 0.52  : Rifaximin kill of MB_SBA (1/d)
KR_BUT  : 0.55  : Rifaximin kill of MB_BUT (1/d)
KR_BAC  : 0.30  : Rifaximin kill of MB_BAC (1/d)
KR_BIF  : 0.52  : Rifaximin kill of MB_BIF (1/d)
KX_SBA  : 1.00  : Index antibiotic kill of MB_SBA (1/d)
KX_BUT  : 1.05  : Index antibiotic kill of MB_BUT (1/d)
KX_BAC  : 0.40  : Index antibiotic kill of MB_BAC (1/d)
KX_BIF  : 0.90  : Index antibiotic kill of MB_BIF (1/d)
EC50X   : 0.30  : Index antibiotic exposure EC50 (exposure units)

$CMT @annotated
// --- 1-6 microbiota guilds
MB_SBA     : bai+ 7a-dehydroxylating Clostridia (relative abundance)
MB_BUT     : Butyrogenic Lachnospiraceae/Ruminococcaceae (relative abundance)
MB_BAC     : Bacteroidetes (relative abundance)
MB_BIF     : Bifidobacterium/Actinobacteria (relative abundance)
MB_ENT     : Enterobacteriaceae (relative abundance)
MB_ENC     : Enterococcus (relative abundance)
// --- 7-11 bile acids
BA_TCA     : Conjugated primary bile acids (uM)
BA_CA      : Cholate (uM)
BA_CDCA    : Chenodeoxycholate (uM)
BA_DCA     : Deoxycholate (uM)
BA_LCA     : Lithocholate (uM)
// --- 12-14 nutrients
NUT_SIA    : Free sialic acid + succinate (mM)
NUT_AA     : Stickland amino acids (mM)
SCFA_BUT   : Luminal butyrate (mM)
// --- 15-18 C. difficile
CD_SPORE_L : Luminal spores (10^6 CFU/g)
CD_VEG     : Luminal vegetative cells (10^6 CFU/g)
CD_MUC     : Mucosa-adherent vegetative cells (10^6 CFU/g)
CD_SPORE_B : Mucosa/biofilm spore reservoir (10^6 CFU/g)
// --- 19-24 toxins
TCDA       : Luminal TcdA (ng/mL)
TCDB       : Luminal TcdB (ng/mL)
TOX_CPLX   : Bezlotoxumab-TcdB complex (ng/mL)
TCDA_MUC   : Mucosa-bound TcdA (ng/mL-equiv)
TCDB_MUC   : Mucosa-bound TcdB (ng/mL-equiv)
CDT        : Binary toxin CDT (ng/mL)
// --- 25-29 epithelium
EPI        : Viable colonocyte fraction
EPI_SC     : Lgr5+ stem/progenitor pool (fraction of normal)
EPI_TJ     : Tight junction integrity (fraction of normal)
EPI_MUCUS  : MUC2 mucus layer (fraction of normal)
EPI_PERM   : Paracellular permeability index
// --- 30-36 immunity
IM_IL8     : IL-8/CXCL8 (normalised units)
IM_IL1B    : IL-1beta (normalised units)
IM_TNF     : TNF-alpha/IL-6 composite (normalised units)
IM_NEUT    : Mucosal neutrophil density (normalised units)
IM_IL22    : IL-22 (normalised units)
IM_PSM     : Pseudomembrane burden (normalised units)
AB_IGG     : Endogenous neutralising anti-TcdB IgG (normalised units)
// --- 37-41 clinical
EPI_H2O    : Luminal fluid load (normalised units)
STOOL      : Unformed stools per 24 h
WBC        : Peripheral WBC (10^9/L)
ALB        : Serum albumin (g/dL)
CRE        : Serum creatinine (mg/dL)
// --- 42-44 oral vancomycin
VAN_D      : Vancomycin, gastric/dosing compartment (mg)
VAN_T      : Vancomycin, small-bowel transit compartment (mg)
VAN_COL    : Vancomycin, colonic content (mg)
// --- 45-48 fidaxomicin + OP-1118
FDX_D      : Fidaxomicin, gastric/dosing compartment (mg)
FDX_T      : Fidaxomicin, small-bowel transit compartment (mg)
FDX_COL    : Fidaxomicin, colonic content (mg)
OP_COL     : OP-1118 metabolite, colonic content (mg)
// --- 49-51 metronidazole
MTZ_A      : Metronidazole, absorption site (mg)
MTZ_C      : Metronidazole, systemic (mg)
MTZ_COL    : Metronidazole, colonic content (mg)
// --- 52-54 bezlotoxumab
BEZ_C      : Bezlotoxumab, central plasma (mg)
BEZ_P      : Bezlotoxumab, peripheral (mg)
BEZ_GUT    : Bezlotoxumab, luminal (mg/L)
// --- 55-56 rifaximin
RFX_D      : Rifaximin, dosing compartment (mg)
RFX_COL    : Rifaximin, colonic content (mg)
// --- 57-58 ridinilazole
RID_D      : Ridinilazole, dosing compartment (mg)
RID_COL    : Ridinilazole, colonic content (mg)
// --- 59-60 index (precipitating) antibiotic
ABX_D      : Index antibiotic, dosing compartment (exposure units)
ABX_COL    : Index antibiotic, colonic exposure (exposure units)
// --- 61 live biotherapeutic / FMT
LBP        : Administered microbial consortium, unengrafted dose

$GLOBAL
#define pos(x) ((x) > 0.0 ? (x) : 0.0)

// faecal / systemic concentrations
#define CVAN   (1000.0 * VAN_COL / VCOLG)
#define CFDX   (1000.0 * FDX_COL / VCOLG)
#define COP    (1000.0 * OP_COL  / VCOLG)
#define CMTZ   (1000.0 * MTZ_COL / VCOLG)
#define CRFX   (1000.0 * RFX_COL / VCOLG)
#define CRID   (1000.0 * RID_COL / VCOLG)
#define CFDXE  (CFDX + WOP * COP)
#define CBEZP  (BEZ_C / VCB)

// fractional target engagement of each agent (0-1)
#define FV     (CVAN  / (EC50V + CVAN))
#define FF     (CFDXE / (EC50F + CFDXE))
#define FM     (CMTZ  / (EC50M + CMTZ))
#define FD     (CRID  / (EC50D + CRID))
#define FR     (CRFX  / (EC50R + CRFX))
#define FX     (ABX_COL / (EC50X + ABX_COL))

double ANA0;
double BSHREF;
double TCASS;
double CASS;
double CDCASS;
double DCASS;
double LCASS;
double SIASS;
double AASS;
double BUTSS;
double bq;
double disc;
double PCA;
double PCD;

$MAIN
// ---------------------------------------------------------------------------
// Self-calibrating healthy colonic steady state.
// Solved from the parameters so that the "no infection, no drug" run is a
// true fixed point rather than an assumed one.
// ---------------------------------------------------------------------------
ANA0   = SBA0 + BUT0 + BAC0 + BIF0;                 // obligate-anaerobe niche
BSHREF = BAC0 + BIF0 + BUT0 + WENCBSH * ENC0;       // reference BSH capacity

// conjugated pool: input balanced by deconjugation + loss
TCASS  = TCAIN / (KDEC + KABST);

// cholate: production = 7a-dehydroxylation (Michaelis-Menten) + loss
//   KABSCA*CA^2 + (K7A*SBA0 + KABSCA*KM7 - P)*CA - P*KM7 = 0
PCA  = FCA * KDEC * TCASS;
bq   = K7A * SBA0 + KABSCA * KM7 - PCA;
disc = bq * bq + 4.0 * KABSCA * PCA * KM7;
CASS = (-bq + sqrt(disc)) / (2.0 * KABSCA);
DCASS = (K7A * SBA0 * CASS / (KM7 + CASS)) / KABSDCA;
// chenodeoxycholate and lithocholate: same algebra
PCD  = (1.0 - FCA) * KDEC * TCASS;
bq   = K7C * SBA0 + KABSCD * KM7 - PCD;
disc = bq * bq + 4.0 * KABSCD * PCD * KM7;
CDCASS = (-bq + sqrt(disc)) / (2.0 * KABSCD);
LCASS  = (K7C * SBA0 * CDCASS / (KM7 + CDCASS)) / KABSLCA;
SIASS = (JSIA0 + JSIAM) / (KCSIA + KOUTSIA);
AASS  = JAA / (KCAA + KOUTAA);
BUTSS = KBUTP * FIBER / KABSBUT;

MB_SBA_0  = SBA0;   MB_BUT_0 = BUT0;   MB_BAC_0 = BAC0;
MB_BIF_0  = BIF0;   MB_ENT_0 = ENT0;   MB_ENC_0 = ENC0;
BA_TCA_0  = TCASS;  BA_CA_0  = CASS;   BA_CDCA_0 = CDCASS;
BA_DCA_0  = DCASS;  BA_LCA_0 = LCASS;
NUT_SIA_0 = SIASS;  NUT_AA_0 = AASS;   SCFA_BUT_0 = BUTSS;

EPI_0 = 1.0;  EPI_SC_0 = 1.0;  EPI_TJ_0 = 1.0;  EPI_MUCUS_0 = 1.0;
AB_IGG_0 = IGG0;
STOOL_0 = STOOL0;  WBC_0 = WBC0;  ALB_0 = ALB0;  CRE_0 = CRE0;

$ODE
// ===========================================================================
// 1. ECOLOGY — guild dynamics under drug-specific kill
// ===========================================================================
double sba = pos(MB_SBA);
double but = pos(MB_BUT);
double bac = pos(MB_BAC);
double bif = pos(MB_BIF);
double ent = pos(MB_ENT);
double enc = pos(MB_ENC);
double ANA  = sba + but + bac + bif;
double fANA = ANA / ANA0;                        // 1 at health, ->0 in dysbiosis
double fdys = pos(1.0 - fANA);
// Bacteroidetes LIBERATE sialic acid rather than consume it (Ng 2013), and
// Stickland fermentation is a Clostridia/Actinobacteria trait.  So the guild
// set that competes with C. difficile for nutrients excludes MB_BAC — which
// is exactly why vancomycin (Bacteroides-sparing) still leaves the niche open.
double FERM  = sba + but + bif;
double fFERM = FERM / (SBA0 + BUT0 + BIF0);

double fVRE = (1.0 - VRE);
double KILL_SBA = KV_SBA*FV + KF_SBA*FF + KM_SBA*FM + KD_SBA*FD + KR_SBA*FR + KX_SBA*FX;
double KILL_BUT = KV_BUT*FV + KF_BUT*FF + KM_BUT*FM + KD_BUT*FD + KR_BUT*FR + KX_BUT*FX;
double KILL_BAC = KV_BAC*FV + KF_BAC*FF + KM_BAC*FM + KD_BAC*FD + KR_BAC*FR + KX_BAC*FX;
double KILL_BIF = KV_BIF*FV + KF_BIF*FF + KM_BIF*FM + KD_BIF*FD + KR_BIF*FR + KX_BIF*FX;
double KILL_ENC = KV_ENC*FV*fVRE + KF_ENC*FF + KD_ENC*FD;

// engraftment of an administered consortium, blocked while antibiotics persist
double ABXPRESS = FV + FF + FM + FD + FR + FX;
double ENGFLUX  = KENG * pos(LBP) / (1.0 + ABXPRESS / KABXENG);

// accelerated transit and mucosal inflammation hold the commensals down, so
// the infection perpetuates the very dysbiosis that permitted it
double ECOSUP = 1.0/(1.0 + pos(IM_NEUT)/KNSUP + pos(STOOL - STOOL0)/KSSUP);

// Re-establishment flux.  Guilds do not climb back purely by multiplying what
// survived: they are re-seeded from the ileum, the biofilm and the environment.
// This zeroth-order limb is why recovery takes weeks almost independently of
// how deep the nadir was, and why the DEPTH of the nadir under drug is set by
// the kill/re-seeding balance rather than by an exponential decay.
// Incoming propagules are killed by whatever is in the lumen, so the brake on
// re-seeding is the GUILD-SPECIFIC kill, not total antibacterial activity.
// This is the whole narrow-vs-broad-spectrum argument in one line: fidaxomicin
// engages its target as completely as vancomycin does, yet leaves the door
// open for the Clostridia to come back.
double RECOL_SBA = (SBA0/TRECOL)*ECOSUP/(1.0 + KILL_SBA/KRECK);
double RECOL_BUT = (BUT0/TRECOL)*ECOSUP/(1.0 + KILL_BUT/KRECK);
double RECOL_BAC = (BAC0/TRECOL)*ECOSUP/(1.0 + KILL_BAC/KRECK);
double RECOL_BIF = (BIF0/TRECOL)*ECOSUP/(1.0 + KILL_BIF/KRECK);
dxdt_MB_SBA = MUSBA*ECOSUP*sba*(1.0 - sba/SBA0) + RECOL_SBA*pos(1.0 - sba/SBA0)
              - KILL_SBA*sba + EPSMB*SBA0 + FE_SBA*ENGFLUX;
dxdt_MB_BUT = MUBUT*ECOSUP*but*(1.0 - but/BUT0) + RECOL_BUT*pos(1.0 - but/BUT0)
              - KILL_BUT*but + EPSMB*BUT0 + FE_BUT*ENGFLUX;
dxdt_MB_BAC = MUBAC*ECOSUP*bac*(1.0 - bac/BAC0) + RECOL_BAC*pos(1.0 - bac/BAC0)
              - KILL_BAC*bac + EPSMB*BAC0 + FE_BAC*ENGFLUX;
dxdt_MB_BIF = MUBIF*ECOSUP*bif*(1.0 - bif/BIF0) + RECOL_BIF*pos(1.0 - bif/BIF0)
              - KILL_BIF*bif + EPSMB*BIF0 + FE_BIF*ENGFLUX;

// opportunists: their carrying capacity IS the emptiness of the anaerobe niche
double KENT_EFF = (ENT0 + (KENTMAX - ENT0)*pow(fdys, HDYS)) / (1.0 + IM_IL22/KI22ENT);
double KENC_EFF = (ENC0 + (KENCMAX - ENC0)*pow(fdys, HDYS)) / (1.0 + IM_IL22/KI22ENT);
dxdt_MB_ENT = MUENT*ent*(1.0 - ent/KENT_EFF) + EPSMB*ENT0;
dxdt_MB_ENC = MUENC*enc*(1.0 - enc/KENC_EFF) - KILL_ENC*enc + EPSMB*ENC0;

// ===========================================================================
// 2. BILE ACIDS — the germination switch
// ===========================================================================
double BSHACT = BSHMIN + (1.0 - BSHMIN) *
                fmin(1.0, (bac + bif + but + WENCBSH*enc) / BSHREF);
double DECON  = KDEC * BSHACT * pos(BA_TCA);
double V7CA   = K7A * sba * pos(BA_CA)   / (KM7 + pos(BA_CA));
double V7CD   = K7C * sba * pos(BA_CDCA) / (KM7 + pos(BA_CDCA));

dxdt_BA_TCA  = TCAIN - DECON - KABST*BA_TCA;
dxdt_BA_CA   = FCA*DECON - V7CA - KABSCA*BA_CA;
dxdt_BA_CDCA = (1.0-FCA)*DECON - V7CD - KABSCD*BA_CDCA;
dxdt_BA_DCA  = V7CA - KABSDCA*BA_DCA;
dxdt_BA_LCA  = V7CD - KABSLCA*BA_LCA;

double SBA2 = pos(BA_DCA) + 0.6*pos(BA_LCA);     // combined secondary-BA signal

// ===========================================================================
// 3. NUTRIENT NICHE AND SCFA
// ===========================================================================
double veg = pos(CD_VEG);
double muc = pos(CD_MUC);
double SIAIN = JSIA0 + JSIAM*(FSIAH + (1.0-FSIAH)*fmin(1.0, bac/BAC0));

dxdt_NUT_SIA = SIAIN - KCSIA*fFERM*NUT_SIA - KCSIACD*veg*NUT_SIA - KOUTSIA*NUT_SIA;
dxdt_NUT_AA  = JAA   - KCAA *fFERM*NUT_AA  - KCAACD *veg*NUT_AA  - KOUTAA *NUT_AA;
dxdt_SCFA_BUT= KBUTP*FIBER*(but/BUT0) - KABSBUT*SCFA_BUT;

// ===========================================================================
// 4. C. DIFFICILE LIFE CYCLE
// ===========================================================================
// germination: taurocholate agonist, CDCA/UDCA competitive antagonists,
// secondary bile acids non-competitively damping the whole step
double GERMR = KGERM * pos(BA_TCA) /
               (KTCA*(1.0 + pos(BA_CDCA)/KICDCA + UDCA/KIUDCA) + pos(BA_TCA)) /
               (1.0 + SBA2/KI2GERM);

// colonization resistance is mostly, but not only, a fermenter-guild property
double KCD_EFF = KCDMAX * exp(-LAMNICHE * ((1.0-WNICHEF)*fANA + WNICHEF*fFERM));
double MUCD = MUCDMAX * (pos(NUT_SIA)/(KSIA + pos(NUT_SIA)))
                      * (pos(NUT_AA) /(KAA  + pos(NUT_AA)))
                      / (1.0 + SBA2/KI2GROW)
                      / (1.0 + pos(SCFA_BUT)/KIBUTCD);

double KILLCD = EMAXV*FV + EMAXF*FF + EMAXM*FM + EMAXD*FD + EMAXR*FR + EMAXT*TIGE;
double SPORF  = KSPOR * (1.0 + RT027*RTSPOR) / (1.0 + SBA2/KI2SPOR)
                      / (1.0 + CFDXE/IC50SF);

double ADH = KADH * (0.25 + 0.75*(1.0 - pos(EPI_MUCUS)));   // mucus shields mucosa

dxdt_CD_SPORE_L = SPORF*(1.0-FSPORB)*veg - GERMR*pos(CD_SPORE_L)
                  - KOUTSP*(1.0 + KWASH*pos(STOOL-STOOL0))*pos(CD_SPORE_L);

dxdt_CD_VEG = GERMR*pos(CD_SPORE_L) + MUCD*veg*(1.0 - veg/KCD_EFF)
              - KDCD*veg - SPORF*veg - ADH*veg + KDET*muc
              - KILLCD*veg - KOUTCD*(1.0 + KWASH*pos(STOOL-STOOL0))*veg;

// the mucus-associated population is nutrient-privileged but is NOT exempt
// from colonization resistance: it grows on the same MUCD, times FMUCADV
double MUMUCE = MUCD * FMUCADV;
dxdt_CD_MUC = ADH*veg - KDET*muc + MUMUCE*muc*(1.0 - muc/KMUCMAX)
              - KDMUC*muc - KILLCD*FMUCPEN*muc - SPORF*muc
              - 0.20*pos(AB_IGG)*muc;

dxdt_CD_SPORE_B = SPORF*FSPORB*veg + SPORF*muc
                  - KGERMB*GERMR/KGERM*pos(CD_SPORE_B) - KCLRB*pos(CD_SPORE_B);

// the reservoir re-seeds the lumen when germination conditions return
double RESEED = KGERMB*GERMR/KGERM*pos(CD_SPORE_B);
dxdt_CD_VEG += RESEED;

// ===========================================================================
// 5. PaLoc REGULATION AND TOXINS
// ===========================================================================
// CodY/CcpA nutrient repression: toxin rises as the population exhausts its
// own niche, so toxin lags the bacterial peak (as observed in vivo).
double FTOX = FTOXB * (KNTOX/(KNTOX + pos(NUT_AA)))
                    / (1.0 + pos(SCFA_BUT)/KIBUTOX)
                    * (1.0 + RT027*(RTTOXF - 1.0))
                    / (1.0 + CFDXE/IC50TF);

double TOXSRC = veg + WMUC*muc;
// Neutralising capacity at the mucosal surface: endogenous anti-toxin IgG plus
// transudated bezlotoxumab.  BEZ_GUT is in ug/mL, converted to IgG-equivalent
// units by WBEZ, so a single 10 mg/kg dose contributes on the same scale as a
// good natural humoral response.
double ABEQ   = pos(AB_IGG) + WBEZ*pos(BEZ_GUT);
double fnAB   = 1.0 / (1.0 + ABEQ/KIAB);
// TcdA/TcdB attack an intact epithelium; rising permeability then amplifies
// delivery, which is the feed-forward limb of the injury loop
double PERMF  = 1.0 + pos(EPI_PERM);

dxdt_TCDA = KTOXA*TOXSRC*FTOX - KDEGA*TCDA - KOUTTOX*TCDA - KTRA*TCDA*PERMF;
double BIND = KONB*pos(TCDB)*pos(BEZ_GUT);       // ng/mL/d of TcdB captured
dxdt_TCDB = KTOXB*TOXSRC*FTOX - KDEGB*TCDB - KOUTTOX*TCDB - KTRB*TCDB*PERMF
            - BIND + KOFFB*pos(TOX_CPLX);
dxdt_TOX_CPLX = BIND - KOFFB*pos(TOX_CPLX) - KDEGCX*pos(TOX_CPLX);
dxdt_TCDA_MUC = KTRA*TCDA*PERMF*fnAB - KINTA*TCDA_MUC;
dxdt_TCDB_MUC = KTRB*TCDB*PERMF*fnAB - KINTB*TCDB_MUC;
dxdt_CDT      = KCDTP*RT027*TOXSRC*FTOX - KDEGCDT*CDT;

double TOXLOAD = pos(TCDB_MUC) + WA*pos(TCDA_MUC);
double CDTF    = 1.0 + WCDT*pos(CDT)/(KCDT50 + pos(CDT));

// ===========================================================================
// 6. EPITHELIUM AND BARRIER
// ===========================================================================
double DMG  = KDMG * TOXLOAD/(KDMG50 + TOXLOAD) * CDTF;
double FWNT = 1.0/(1.0 + pos(TCDB_MUC)/KFZD);
double CYT  = pos(IM_IL1B) + pos(IM_TNF);

dxdt_EPI      = KREG*pos(EPI_SC)*(1.0 - EPI) - DMG*EPI - KAPOP*CYT*EPI
                - 0.030*pos(IM_NEUT)*EPI + 0.040*pos(IM_IL22)*(1.0 - EPI);
dxdt_EPI_SC   = KSCR*(1.0 - EPI_SC)*FWNT
                - KSCL*(pos(TCDB_MUC)/(KSC50 + pos(TCDB_MUC)))*EPI_SC;
dxdt_EPI_TJ   = KTJR*(1.0 - EPI_TJ)*pos(EPI)
                - KTJT*(TOXLOAD/(KTJ50 + TOXLOAD))*EPI_TJ
                - KTJI*CYT*EPI_TJ
                + 0.15*(pos(SCFA_BUT)/BUTSS)*(1.0 - EPI_TJ);
dxdt_EPI_MUCUS= KMUCR*pos(EPI)*(0.4 + 0.6*pos(SCFA_BUT)/BUTSS)*(1.0 - EPI_MUCUS)
                - KMUCL*(pos(IM_NEUT) + 0.5*TOXLOAD/(KTJ50 + TOXLOAD))*EPI_MUCUS;
dxdt_EPI_PERM = KPON*(pos(1.0 - EPI_TJ) + 0.6*pos(1.0 - EPI)) - KPOFF*EPI_PERM;

// ===========================================================================
// 7. IMMUNITY
// ===========================================================================
// translocated microbial product: bounded on both factors so the composite
// cytokine scale stays interpretable (0 at health, ~1-3 in severe disease)
double ENTF = ent / (ent + ENTSAT);
double LPSX = fmin(1.0, pos(EPI_PERM)/2.0) * ENTF;

dxdt_IM_IL8  = KIL8*(TOXLOAD/(KI850 + TOXLOAD) + 0.5*pos(IM_IL1B)
                     + 0.3*pos(IM_TNF) + 0.4*LPSX) - KDIL8*IM_IL8;
dxdt_IM_IL1B = KIL1*(TOXLOAD/(KI150 + TOXLOAD))*(1.0 + 0.5*ENTF) - KDIL1*IM_IL1B;
dxdt_IM_TNF  = KTNF*(muc/(KTN50 + muc) + 0.5*LPSX) - KDTNF*IM_TNF;
dxdt_IM_NEUT = KNREC*(pos(IM_IL8) + 0.6*pos(IM_IL1B)) - KNDEC*IM_NEUT;
dxdt_IM_IL22 = KIL22*(pos(IM_IL1B) + 0.3*pos(IM_TNF))
                     *(1.0 + 0.4*pos(SCFA_BUT)/BUTSS) - KDIL22*IM_IL22;
dxdt_IM_PSM  = KPSM*pos(IM_NEUT)*pos(1.0 - EPI) - KPSMR*IM_PSM*(0.2 + pos(EPI));

double TOXAG = pos(TCDA) + pos(TCDB);
dxdt_AB_IGG  = KABP*IMMCOMP*(TOXAG/(KABA50 + TOXAG))*pos(1.0 - AB_IGG/ABMAX)
               - KABD*AB_IGG;

// ===========================================================================
// 8. CLINICAL READ-OUTS
// ===========================================================================
dxdt_EPI_H2O = KH2O*(pos(EPI_PERM) + WSEC*TOXLOAD/(KSEC50 + TOXLOAD)
                     + 0.4*pos(IM_IL1B)) - KH2OA*EPI_H2O*(0.25 + 0.75*pos(EPI));

double ILEUSF = 1.0/(1.0 + pos(IM_PSM)/KILEUS);
double STARG  = STOOL0 + (SMAX - STOOL0)*(pos(EPI_H2O)/(KSTL + pos(EPI_H2O)))*ILEUSF;
dxdt_STOOL    = (STARG - STOOL)/TAUSTL;

double CYT2   = pos(IM_IL8) + pos(IM_IL1B);
double WTARG  = WBC0 + WBCMAX*CYT2/(KWBC + CYT2);
dxdt_WBC      = (WTARG - WBC)/TAUWBC;

dxdt_ALB      = KALBD*ALB0/(1.0 + pos(IM_TNF)/KALBI) - KALBD*ALB
                - KALBL*pos(EPI_PERM)*ALB;

double EXS    = pos(STOOL - STOOL0);
double CTARG  = CRE0*(1.0 + KDEH*EXS/(KSTC + EXS));
dxdt_CRE      = (CTARG - CRE)/TAUCRE;

// ===========================================================================
// 9. DRUG PHARMACOKINETICS
// ===========================================================================
double WASHF = 1.0 + KWASH*EXS;                  // diarrhoeal washout of drug

dxdt_VAN_D = -KTRV*VAN_D;
dxdt_VAN_T =  KTRV*VAN_D - KTRV*VAN_T;
dxdt_VAN_COL = KTRV*VAN_T - KEXV*WASHF*VAN_COL;

dxdt_FDX_D = -KTRF*FDX_D;
dxdt_FDX_T =  KTRF*FDX_D - KTRF*FDX_T;
dxdt_FDX_COL = KTRF*FDX_T - (KEXF*WASHF + KMETF)*FDX_COL;
dxdt_OP_COL  = FMET*KMETF*FDX_COL - KEXOP*WASHF*OP_COL;

dxdt_MTZ_A = -KAM*MTZ_A;
dxdt_MTZ_C =  KAM*MTZ_A*FMTZ - KEM*MTZ_C;
// Metronidazole reaches the colonic lumen almost entirely by secretion across
// INFLAMED mucosa, so its faecal concentration collapses once the mucosa heals
// (Bolton 1986: 9.3 ug/g in watery stool, 1.2 ug/g in formed stool).  That is
// the mechanistic reason metronidazole cannot hold a cure it has achieved.
dxdt_MTZ_COL = KAM*MTZ_A*(1.0-FMTZ) + KSECM*MTZ_C*(0.02 + pos(EPI_PERM))
               - KEXM*WASHF*MTZ_COL;

dxdt_BEZ_C = -(CLB/VCB)*BEZ_C - (QB/VCB)*BEZ_C + (QB/VPB)*BEZ_P;
dxdt_BEZ_P =  (QB/VCB)*BEZ_C - (QB/VPB)*BEZ_P;
// antibody consumed stoichiometrically by complex formation; BIND is in ng/mL
// of TcdB (270 kDa) so the IgG (150 kDa) mass equivalent is BIND*150/270,
// converted from ng/mL to mg/L by /1000.  Negligible at the molar excess
// achieved by a 10 mg/kg dose, but kept so the mass balance is closed.
dxdt_BEZ_GUT = KTRAB*CBEZP*(0.05 + pos(EPI_PERM)) - KDGUT*BEZ_GUT
               - BIND*(150.0/270.0)/1000.0;

dxdt_RFX_D = -KTRR*RFX_D;
dxdt_RFX_COL = KTRR*RFX_D - KEXR*WASHF*RFX_COL;

dxdt_RID_D = -KTRD*RID_D;
dxdt_RID_COL = KTRD*RID_D - KEXD*WASHF*RID_COL;

dxdt_ABX_D = -KTRX*ABX_D;
dxdt_ABX_COL = KTRX*ABX_D - KEXX*ABX_COL;

dxdt_LBP = -KENG*pos(LBP);

$TABLE
// derived, human-readable outputs
double fANAo  = (MB_SBA + MB_BUT + MB_BAC + MB_BIF) / ANA0;
double SBA2o  = pos(BA_DCA) + 0.6*pos(BA_LCA);
double TOXLo  = pos(TCDB_MUC) + WA*pos(TCDA_MUC);
double CVANo  = 1000.0*VAN_COL/VCOLG;
double CFDXo  = 1000.0*FDX_COL/VCOLG;
double COPo   = 1000.0*OP_COL /VCOLG;
double CMTZFo = 1000.0*MTZ_COL/VCOLG;
double CMTZSo = MTZ_C/VMTZ;
double CRIDo  = 1000.0*RID_COL/VCOLG;
double CRFXo  = 1000.0*RFX_COL/VCOLG;
double CBEZo  = BEZ_C/VCB;

// log10 burdens (floored so plots stay finite)
double LGVEG   = log10(fmax(pos(CD_VEG)   , 1e-8)) + 6.0;
double LGMUC   = log10(fmax(pos(CD_MUC)   , 1e-8)) + 6.0;
double LGSPL   = log10(fmax(pos(CD_SPORE_L), 1e-8)) + 6.0;
double LGSPB   = log10(fmax(pos(CD_SPORE_B), 1e-8)) + 6.0;
double LGSBA   = log10(fmax(pos(MB_SBA)   , 1e-12));

// Shannon diversity over the six modelled guilds
double tot = pos(MB_SBA)+pos(MB_BUT)+pos(MB_BAC)+pos(MB_BIF)+pos(MB_ENT)+pos(MB_ENC);
double SHAN = 0.0;
double pk1 = pos(MB_SBA)/fmax(tot,1e-12);
double pk2 = pos(MB_BUT)/fmax(tot,1e-12);
double pk3 = pos(MB_BAC)/fmax(tot,1e-12);
double pk4 = pos(MB_BIF)/fmax(tot,1e-12);
double pk5 = pos(MB_ENT)/fmax(tot,1e-12);
double pk6 = pos(MB_ENC)/fmax(tot,1e-12);
if(pk1>1e-12) SHAN -= pk1*log(pk1);
if(pk2>1e-12) SHAN -= pk2*log(pk2);
if(pk3>1e-12) SHAN -= pk3*log(pk3);
if(pk4>1e-12) SHAN -= pk4*log(pk4);
if(pk5>1e-12) SHAN -= pk5*log(pk5);
if(pk6>1e-12) SHAN -= pk6*log(pk6);

// IDSA/SHEA severity classification
double SEVERE   = ((WBC >= 15.0) || (CRE >= 1.5)) ? 1.0 : 0.0;
double FULMIN   = ((IM_PSM >= 2.2) && (WBC >= 20.0)) ? 1.0 : 0.0;
double SYMPT    = (STOOL >= 3.0) ? 1.0 : 0.0;

// Recurrence risk index: the three things that actually predict relapse —
// a surviving spore reservoir, an unrecovered restoring guild, and no
// neutralising antibody.  Mapped to a probability in the R driver.
double ABTOT = pos(AB_IGG) + WBEZ*pos(BEZ_GUT);
double RRI = (pos(CD_SPORE_B)/(0.15 + pos(CD_SPORE_B)))
             * pos(1.0 - pos(MB_SBA)/SBA0)
             / (1.0 + ABTOT/KIAB)
             / (1.0 + SBA2o/200.0);

$CAPTURE @annotated
fANAo  : Anaerobe niche occupancy (fraction of healthy)
SHAN   : Shannon diversity over modelled guilds
SBA2o  : Combined secondary bile-acid signal (uM)
TOXLo  : Mucosal toxin load (ng/mL-equiv)
LGVEG  : log10 luminal C. difficile (CFU/g)
LGMUC  : log10 mucosa-adherent C. difficile (CFU/g)
LGSPL  : log10 luminal spores (CFU/g)
LGSPB  : log10 mucosal spore reservoir (CFU/g)
LGSBA  : log10 relative abundance of bai+ Clostridia
CVANo  : Faecal vancomycin (ug/g)
CFDXo  : Faecal fidaxomicin (ug/g)
COPo   : Faecal OP-1118 (ug/g)
CMTZFo : Faecal metronidazole (ug/g)
CMTZSo : Plasma metronidazole (mg/L)
CRIDo  : Faecal ridinilazole (ug/g)
CRFXo  : Faecal rifaximin (ug/g)
CBEZo  : Plasma bezlotoxumab (mg/L)
SEVERE : IDSA/SHEA severe CDI flag
FULMIN : Fulminant CDI flag
SYMPT  : Symptomatic (>=3 unformed stools/day) flag
RRI    : Recurrence risk index
'

# =============================================================================
#  BUILD
# =============================================================================
mod <- mcode_cache("cdi", cdi_code, soloc = tempdir())

# -----------------------------------------------------------------------------
#  Dosing helpers.  Compartment names are used directly so the event tables
#  stay readable and survive any renumbering of $CMT.
# -----------------------------------------------------------------------------
ev_index_abx <- function(start = 0, dur = 7, amt = 1) {
  # precipitating antibiotic: daily exposure units for `dur` days
  ev(time = start, amt = amt, cmt = "ABX_D", ii = 1, addl = dur - 1)
}
ev_spore <- function(time = 2, amt = 0.05, ppi = 0) {
  ev(time = time, amt = amt * (1 + 1.5 * ppi), cmt = "CD_SPORE_L")
}
ev_vanco <- function(start, days = 10, dose = 125, ii = 0.25) {
  n <- days / ii
  ev(time = start, amt = dose, cmt = "VAN_D", ii = ii, addl = n - 1)
}
ev_vanco_taper <- function(start) {
  # IDSA taper/pulse for recurrent CDI: 125 qid x14d, tid x7, bid x7, qd x7,
  # then every 2-3 days x 3 weeks
  c(ev(time = start,      amt = 125, cmt = "VAN_D", ii = 0.25,  addl = 55),
    ev(time = start + 14, amt = 125, cmt = "VAN_D", ii = 1/3,   addl = 20),
    ev(time = start + 21, amt = 125, cmt = "VAN_D", ii = 0.5,   addl = 13),
    ev(time = start + 28, amt = 125, cmt = "VAN_D", ii = 1,     addl = 6),
    ev(time = start + 35, amt = 125, cmt = "VAN_D", ii = 2.5,   addl = 8))
}
ev_fidaxo <- function(start, days = 10, dose = 200) {
  ev(time = start, amt = dose, cmt = "FDX_D", ii = 0.5, addl = days * 2 - 1)
}
ev_fidaxo_extend <- function(start) {
  # EXTEND: 200 mg bid days 1-5, then 200 mg once every other day days 7-25
  c(ev(time = start,     amt = 200, cmt = "FDX_D", ii = 0.5, addl = 9),
    ev(time = start + 6, amt = 200, cmt = "FDX_D", ii = 2,   addl = 9))
}
ev_metro <- function(start, days = 10, dose = 500) {
  ev(time = start, amt = dose, cmt = "MTZ_A", ii = 1/3, addl = days * 3 - 1)
}
ev_bezlo <- function(time, wt = 70) {
  ev(time = time, amt = 10 * wt, cmt = "BEZ_C")
}
ev_rifax <- function(start, days = 14, dose = 400) {
  ev(time = start, amt = dose, cmt = "RFX_D", ii = 1/3, addl = days * 3 - 1)
}
ev_ridini <- function(start, days = 10, dose = 200) {
  ev(time = start, amt = dose, cmt = "RID_D", ii = 0.5, addl = days * 2 - 1)
}
ev_lbp <- function(time, amt = 1, n = 1, ii = 1) {
  ev(time = time, amt = amt, cmt = "LBP", ii = ii, addl = n - 1)
}

# consortium composition presets (fractions of the dose that engraft per guild)
LBP_FMT <- list(FE_SBA = 0.130, FE_BUT = 0.320, FE_BAC = 0.450, FE_BIF = 0.090)
LBP_SER <- list(FE_SBA = 0.100, FE_BUT = 0.250, FE_BAC = 0.000, FE_BIF = 0.020)
LBP_RBX <- list(FE_SBA = 0.075, FE_BUT = 0.210, FE_BAC = 0.260, FE_BIF = 0.045)
LBP_NONE<- list(FE_SBA = 0,     FE_BUT = 0,     FE_BAC = 0,     FE_BIF = 0)

# =============================================================================
#  SCENARIOS
#  Timeline convention: index antibiotic runs days 0-6, spores are ingested on
#  day 3, symptoms appear ~day 5-6, therapy starts on day 7 (the clinical
#  presentation), and follow-up continues to day 90 so recurrence is visible.
# =============================================================================
TEND <- 90
TSTART <- 9          # day therapy begins (2-3 days of diarrhoea by then)

base_insult <- function(abx_days = 7, inoc = 0.05, ppi = 0) {
  c(ev_index_abx(0, abx_days), ev_spore(2, inoc, ppi))
}

scen <- list()

scen[["S01 healthy control"]] <- list(
  ev = ev(time = 0, amt = 0, cmt = "LBP"), par = LBP_NONE,
  note = "no antibiotic, no spore exposure — verifies the self-calibrated baseline")

scen[["S02 dysbiosis, no exposure"]] <- list(
  ev = ev_index_abx(0, 7), par = LBP_NONE,
  note = "antibiotic alone: the permissive state without the organism")

scen[["S03 untreated CDI"]] <- list(
  ev = base_insult(), par = LBP_NONE,
  note = "natural history, no anti-C. difficile therapy")

scen[["S04 asymptomatic carriage"]] <- list(
  ev = base_insult(), par = c(LBP_NONE, list(IGG0 = 3.2)),
  note = "same exposure, high pre-existing anti-toxin IgG (Kyne 2000)")

scen[["S05 metronidazole 500 tid x10d"]] <- list(
  ev = c(base_insult(), ev_metro(TSTART, 10)), par = LBP_NONE,
  note = "inferior comparator; faecal drug falls as the mucosa heals")

scen[["S06 vancomycin 125 qid x10d"]] <- list(
  ev = c(base_insult(), ev_vanco(TSTART, 10)), par = LBP_NONE,
  note = "standard of care first episode")

scen[["S07 fidaxomicin 200 bid x10d"]] <- list(
  ev = c(base_insult(), ev_fidaxo(TSTART, 10)), par = LBP_NONE,
  note = "narrow spectrum + sporulation and toxin suppression")

scen[["S08 fidaxomicin extended-pulsed"]] <- list(
  ev = c(base_insult(), ev_fidaxo_extend(TSTART)), par = LBP_NONE,
  note = "EXTEND regimen (Guery 2018)")

scen[["S09 vancomycin + bezlotoxumab"]] <- list(
  ev = c(base_insult(), ev_vanco(TSTART, 10), ev_bezlo(TSTART + 1)), par = LBP_NONE,
  note = "MODIFY I/II strategy: kill the organism, neutralise the toxin")

scen[["S10 fidaxomicin + bezlotoxumab"]] <- list(
  ev = c(base_insult(), ev_fidaxo(TSTART, 10), ev_bezlo(TSTART + 1)), par = LBP_NONE,
  note = "both non-ecological and ecological protection")

scen[["S11 vancomycin then FMT"]] <- list(
  ev = c(base_insult(), ev_vanco(TSTART, 10), ev_lbp(TSTART + 10.5)), par = LBP_FMT,
  note = "van Nood 2013: re-seed the guild 2 days after stopping the drug")

scen[["S12 vancomycin then SER-109"]] <- list(
  ev = c(base_insult(), ev_vanco(TSTART, 10), ev_lbp(TSTART + 11.5, amt = 0.40, n = 3)),
  par = LBP_SER,
  note = "ECOSPOR III: purified Firmicutes spores x3 days")

scen[["S13 vancomycin then RBX2660"]] <- list(
  ev = c(base_insult(), ev_vanco(TSTART, 10), ev_lbp(TSTART + 10.5, amt = 1.0)),
  par = LBP_RBX,
  note = "PUNCH CD3: single rectal microbiota suspension")

scen[["S14 ridinilazole 200 bid x10d"]] <- list(
  ev = c(base_insult(), ev_ridini(TSTART, 10)), par = LBP_NONE,
  note = "Ri-CoDIFy: narrow spectrum, microbiome largely spared")

scen[["S15 vancomycin + rifaximin chaser"]] <- list(
  ev = c(base_insult(), ev_vanco(TSTART, 10), ev_rifax(TSTART + 10, 14)), par = LBP_NONE,
  note = "rifaximin 'chaser' after vancomycin")

scen[["S16 vancomycin taper/pulse"]] <- list(
  ev = c(base_insult(), ev_vanco_taper(TSTART)), par = LBP_NONE,
  note = "IDSA vancomycin taper/pulse for multiply recurrent disease")

scen[["S17 RT027 fulminant, delayed rx"]] <- list(
  ev = c(base_insult(abx_days = 9, inoc = 0.30),
         ev_vanco(TSTART + 5, 14, dose = 500), ev_metro(TSTART + 5, 14)),
  par = c(LBP_NONE, list(RT027 = 1, IGG0 = 0.30, IMMCOMP = 0.35, TIGE = 1)),
  note = "hypervirulent strain, immunosuppressed host, 5 days late to therapy")

scen[["S18 index antibiotic never stopped"]] <- list(
  ev = c(ev_index_abx(0, 45), ev_spore(2), ev_vanco(TSTART, 10)), par = LBP_NONE,
  note = "worst case: continued dysbiotic pressure through and past therapy")

# -----------------------------------------------------------------------------
#  Runner
# -----------------------------------------------------------------------------
run_scen <- function(s, mod. = mod, tend = TEND) {
  p <- s$par
  m <- mod.
  if (length(p)) m <- param(m, p)
  mrgsim_e(m, s$ev, end = tend, delta = 0.25, atol = 1e-8, rtol = 1e-6,
           maxsteps = 200000) %>% as_tibble()
}

# -----------------------------------------------------------------------------
#  Endpoint extraction: clinical cure, time to resolution, recurrence
# -----------------------------------------------------------------------------
# Last day on anti-C. difficile therapy, read straight off the event table so
# taper/pulse and extended-pulsed regimens are handled without bookkeeping.
# Intervention compartments: the assessment point is after BOTH the last
# antibacterial dose and the last microbiome-restoration dose, because the
# ecological state at that moment is what the recurrence index is about.
RX_CMT <- c("VAN_D", "FDX_D", "MTZ_A", "RID_D", "RFX_D", "LBP", "BEZ_C")
rx_end <- function(e) {
  d <- as.data.frame(e)
  d <- d[d$cmt %in% RX_CMT & d$amt > 0, , drop = FALSE]
  if (!nrow(d)) return(NA_real_)
  ii   <- if ("ii"   %in% names(d)) ifelse(is.na(d$ii),   0, d$ii)   else 0
  addl <- if ("addl" %in% names(d)) ifelse(is.na(d$addl), 0, d$addl) else 0
  max(d$time + ii * addl)
}

# Clinical cure, per the trial definition: <=3 unformed stools per day for two
# consecutive days.  Only meaningful in a patient who was symptomatic to begin
# with, so an asymptomatic run reports cure_day = NA rather than "cured on day 1".
endpoints <- function(d, tstart = TSTART, t_rx_end = NA_real_,
                      stool_thresh = 3.0, tox_thresh = 3.0) {
  day <- d %>% filter(time >= 0, abs(time - round(time)) < 1e-6) %>%
    group_by(time) %>% slice(1) %>% ungroup()
  tt <- day$time; stool <- day$STOOL

  symptomatic <- any(stool >= stool_thresh)
  onset_day <- if (symptomatic) tt[which(stool >= stool_thresh)[1]] else NA_real_

  ok <- stool <= stool_thresh
  cure_day <- NA_real_
  if (symptomatic) {
    idx <- which(tt >= max(tstart, onset_day) + 1)
    for (i in idx) if (i < length(ok) && ok[i] && ok[i + 1]) { cure_day <- tt[i]; break }
  }
  # Recurrence, as the trials define it: recurrent diarrhoea AND a positive
  # toxin test -- so two consecutive symptomatic days with toxin present, not a
  # single-day blip while a live biotherapeutic is still engrafting.
  recur_day <- NA_real_
  if (!is.na(cure_day)) {
    bad <- stool > stool_thresh & day$TOXLo > tox_thresh
    after <- which(tt > cure_day + 2)
    for (i in after) if (i < length(bad) && bad[i] && bad[i + 1]) { recur_day <- tt[i]; break }
  }
  at <- function(x, d0) x[which.min(abs(tt - d0))]
  assess <- if (is.na(t_rx_end)) 21 else min(max(t_rx_end, tstart) + 3, TEND)

  tibble(
    symptomatic  = symptomatic,
    onset_day    = onset_day,
    # peak during the index episode vs peak over the whole 90 days: a relapse
    # is often worse than the episode that preceded it, and conflating the two
    # makes an effective drug look like a harmful one
    peak_stool   = if (is.na(t_rx_end)) max(stool)
                   else max(stool[tt <= min(t_rx_end + 2, TEND)]),
    peak_stool_90= max(stool),
    peak_WBC     = max(day$WBC),
    min_ALB      = min(day$ALB),
    max_CRE      = max(day$CRE),
    min_EPI      = min(day$EPI),
    peak_toxin   = max(day$TCDB),
    severe_days  = sum(day$SEVERE),
    cure_day     = cure_day,
    ttrod        = if (is.na(cure_day)) NA_real_ else cure_day - tstart,
    recur_day    = recur_day,
    recurred     = !is.na(recur_day),
    sick_days    = sum(stool >= stool_thresh),
    rx_end       = t_rx_end,
    assess_day   = assess,
    # ecology at the end of therapy: the state that decides what happens next
    sba_frac_rx  = at(day$MB_SBA, assess) / 0.10,
    dca_rx       = at(day$BA_DCA, assess),
    tca_rx       = at(day$BA_TCA, assess),
    shan_rx      = at(day$SHAN,   assess),
    spb_rx       = at(day$LGSPB,  assess),
    igg_rx       = at(day$AB_IGG, assess),
    RRI_rx       = at(day$RRI,    assess),
    sba_frac_d40 = at(day$MB_SBA, 40) / 0.10,
    dca_d40      = at(day$BA_DCA, 40)
  )
}

# ---------------------------------------------------------------------------
#  From the mechanistic recurrence risk index to an 8-week recurrence rate.
#
#  The deterministic trajectory either relapses or it does not; a real cohort
#  splits.  RRI_rx is the model quantity that should order the arms -- surviving
#  spore reservoir x unrecovered restoring guild / available neutralising
#  antibody -- and this logistic is fitted (below, at run time) to the pooled
#  phase 3 rates so the mapping is explicit and auditable rather than hidden.
# ---------------------------------------------------------------------------
# Caveat kept explicit rather than buried: the FMT and SER-109 trials enrolled
# recurrent-CDI populations whose untreated recurrence risk is higher than the
# first-episode arms above, so these two anchors are conservative (they would
# look even better against their own controls).  RBX2660 (PUNCH CD3, 29.4%
# recurrence against a 42.5% control) is deliberately NOT used as an anchor for
# the same reason -- its control arm is not the same population as Louie 2011.
RRI_ANCHORS <- tibble::tribble(
  ~scenario,                        ~observed, ~source,
  "S05 metronidazole 500 tid x10d",     0.230, "Johnson 2014 CID pooled",
  "S06 vancomycin 125 qid x10d",        0.253, "Louie 2011 / Cornely 2012",
  "S07 fidaxomicin 200 bid x10d",       0.154, "Louie 2011 / Cornely 2012",
  "S09 vancomycin + bezlotoxumab",      0.165, "Wilcox 2017 MODIFY I/II",
  "S11 vancomycin then FMT",            0.090, "van Nood 2013 / Kelly 2016",
  "S12 vancomycin then SER-109",        0.120, "Feuerstadt 2022 ECOSPOR III"
)

fit_rri_logistic <- function(res) {
  d <- dplyr::inner_join(res, RRI_ANCHORS, by = "scenario")
  y <- log(d$observed / (1 - d$observed))
  fit <- stats::lm(y ~ d$RRI_rx)
  list(b0 = unname(coef(fit)[1]), b1 = unname(coef(fit)[2]),
       n = nrow(d), r2 = summary(fit)$r.squared, data = d)
}
rri_to_prob <- function(rri, b0, b1) 1 / (1 + exp(-(b0 + b1 * rri)))

# =============================================================================
#  EXECUTE
# =============================================================================
if (!exists("CDI_SOURCE_ONLY")) {

  cat("\n=====================================================================\n")
  cat(" Clostridioides difficile infection QSP model\n")
  cat(" ", length(names(mrgsolve::init(mod))), "compartments |",
      length(names(param(mod))), "parameters\n")
  cat("=====================================================================\n\n")

  # ---- 1. baseline verification -------------------------------------------
  b <- run_scen(scen[["S01 healthy control"]])
  bl  <- b %>% filter(time == TEND) %>% slice(1)
  cat("--- Healthy steady state (day 90 vs day 0, self-calibrated) ---\n")
  b0r <- b %>% filter(time == 0) %>% slice(1)
  chk <- tibble(
    state = c("MB_SBA","BA_TCA","BA_CA","BA_CDCA","BA_DCA","BA_LCA",
              "NUT_SIA","NUT_AA","SCFA_BUT","EPI","STOOL","WBC","ALB","CRE"),
    d0 = c(b0r$MB_SBA,b0r$BA_TCA,b0r$BA_CA,b0r$BA_CDCA,b0r$BA_DCA,b0r$BA_LCA,
           b0r$NUT_SIA,b0r$NUT_AA,b0r$SCFA_BUT,b0r$EPI,b0r$STOOL,b0r$WBC,
           b0r$ALB,b0r$CRE),
    d90 = c(bl$MB_SBA,bl$BA_TCA,bl$BA_CA,bl$BA_CDCA,bl$BA_DCA,bl$BA_LCA,
            bl$NUT_SIA,bl$NUT_AA,bl$SCFA_BUT,bl$EPI,bl$STOOL,bl$WBC,
            bl$ALB,bl$CRE))
  chk$drift_pct <- 100 * (chk$d90 - chk$d0) / pmax(abs(chk$d0), 1e-9)
  print(as.data.frame(chk), digits = 4, row.names = FALSE)
  cat("\nMax |drift| over 90 drug-free days:",
      sprintf("%.3f%%\n\n", max(abs(chk$drift_pct))))

  # ---- 2. all scenarios ---------------------------------------------------
  res <- lapply(names(scen), function(nm) {
    d <- run_scen(scen[[nm]])
    e <- endpoints(d, t_rx_end = rx_end(scen[[nm]]$ev))
    e$scenario <- nm
    e
  }) %>% bind_rows() %>% relocate(scenario)

  lg <- fit_rri_logistic(res)
  res$p_recur <- rri_to_prob(res$RRI_rx, lg$b0, lg$b1)

  cat("--- Acute-phase severity ---\n")
  print(as.data.frame(res %>% select(scenario, symptomatic, onset_day, peak_stool,
                                     peak_stool_90, peak_WBC, max_CRE, min_ALB,
                                     min_EPI, peak_toxin, severe_days)),
        digits = 3, row.names = FALSE)

  cat("\n--- Response, relapse trajectory, and mapped 8-week recurrence rate ---\n")
  print(as.data.frame(res %>% select(scenario, cure_day, ttrod, sick_days,
                                     recurred, recur_day, RRI_rx, p_recur)),
        digits = 3, row.names = FALSE)

  cat("\n--- Ecology at the end of therapy (the state that decides the sequel) ---\n")
  print(as.data.frame(res %>% select(scenario, rx_end, sba_frac_rx, dca_rx,
                                     tca_rx, shan_rx, spb_rx, igg_rx,
                                     sba_frac_d40, dca_d40)),
        digits = 3, row.names = FALSE)

  cat("\n--- Calibration of the RRI -> recurrence-probability mapping ---\n")
  cat(sprintf("  logit(p) = %.3f + %.3f * RRI    (n = %d anchors, R2 = %.3f)\n",
              lg$b0, lg$b1, lg$n, lg$r2))
  print(as.data.frame(lg$data %>% transmute(scenario, RRI_rx, observed,
                        predicted = rri_to_prob(RRI_rx, lg$b0, lg$b1), source)),
        digits = 3, row.names = FALSE)

  cat("\n--- Faecal / plasma drug exposure check ---\n")
  pk <- lapply(c("S06 vancomycin 125 qid x10d","S07 fidaxomicin 200 bid x10d",
                 "S05 metronidazole 500 tid x10d","S14 ridinilazole 200 bid x10d",
                 "S09 vancomycin + bezlotoxumab","S15 vancomycin + rifaximin chaser"),
    function(nm) {
      d <- run_scen(scen[[nm]]) %>% filter(time > TSTART + 3, time < TSTART + 9)
      tibble(scenario = nm,
             van_ugg  = mean(d$CVANo), fdx_ugg = mean(d$CFDXo),
             op_ugg   = mean(d$COPo),  mtz_fec = mean(d$CMTZFo),
             mtz_plas = mean(d$CMTZSo), rid_ugg = mean(d$CRIDo),
             rfx_ugg  = mean(d$CRFXo), bez_mgL = mean(d$CBEZo))
    }) %>% bind_rows()
  print(as.data.frame(pk), digits = 3, row.names = FALSE)

  cat("\n--- Metronidazole: faecal drug follows mucosal inflammation, not dose ---\n")
  dm <- run_scen(scen[["S05 metronidazole 500 tid x10d"]])
  print(as.data.frame(dm %>% filter(time %in% c(7,8,10,12,14,16)) %>%
          group_by(time) %>% slice(1) %>% ungroup() %>%
          select(time, EPI_PERM, CMTZFo, CMTZSo, LGVEG)),
        digits = 3, row.names = FALSE)

  cat("\n=====================================================================\n")
  cat(" Reference targets: vancomycin recurrence ~25%, fidaxomicin ~15%,\n")
  cat(" bezlotoxumab add-on ~16%, FMT/SER-109 ~10-12%, metronidazole ~27%.\n")
  cat("=====================================================================\n")

  invisible(res)
}
