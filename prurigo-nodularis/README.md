# Prurigo Nodularis (PN) — QSP Model

> Epidermal barrier dysfunction → keratinocyte alarmins (TSLP/IL-33/IL-25)
> → Th2/Th17/Th22 skewing → **IL-31 → IL-31RA/OSMRβ** peripheral nerve
> sensitization → central (spinal GRPR) sensitization + MOR:KOR opioid
> imbalance → conscious itch → scratch → mechanical trauma → TGF-β/
> fibroblast-driven **fibrotic nodule formation**, closing a self-
> reinforcing itch-scratch-fibrosis cycle. Modeled therapies: dupilumab
> (anti-IL-4Rα), nemolizumab (anti-IL-31RA), gabapentin, nalbuphine ER
> (MOR antagonist/KOR agonist), and off-label oral JAK1 inhibition.

---

## 1. Disease overview

Prurigo nodularis (PN) is a chronic neuroimmune skin disease defined by
multiple intensely pruritic, hyperkeratotic nodules arising from a
self-perpetuating **itch-scratch cycle**:

1. **Barrier dysfunction & alarmins** — impaired filaggrin/lipid barrier
   raises transepidermal water loss (TEWL) and triggers keratinocyte
   release of TSLP, IL-33, and IL-25, priming ILC2 and Th2 differentiation.
2. **Type-2/Th17 cytokine skewing** — IL-4, IL-13, IL-5, IL-17A, and IL-22
   drive epidermal hyperplasia, while Th2/Th22 cells and keratinocytes
   secrete **IL-31**, the dominant pruritogenic cytokine in PN.
3. **Peripheral neuro-sensitization** — IL-31 signals through the
   IL-31RA/OSMRβ heterodimer on cutaneous C-fibers; TRPV1/TRPA1 channels,
   NGF-driven nerve sprouting, and loss of the axon-repellent Semaphorin
   3A produce a biphasic nerve fiber phenotype (early sprouting →
   late intraepidermal nerve fiber density (IENFD) depletion within
   chronic nodules, i.e. "denervated but hypersensitive" skin).
4. **Central sensitization & opioid imbalance** — chronic peripheral
   input sensitizes spinal GRPR+/Nppb+ itch-relay neurons; an increased
   µ-opioid : κ-opioid receptor tone (amplified in CKD-associated
   pruritus) further raises central itch gain.
5. **Itch-scratch-fibrosis cycle** — conscious itch perception drives
   scratching, causing mechanical trauma that (a) re-injures the barrier
   (feedback loop), and (b) activates dermal fibroblasts via TGF-β,
   depositing collagen that forms the pathognomonic **fibrotic nodules**
   and lichenification.

PN disproportionately affects older, female, and Black patients, and
frequently co-occurs with atopic disease, CKD-associated pruritus,
cholestatic liver disease, HIV, and psychiatric comorbidity (anxiety,
depression) driven by the chronic itch-sleep-quality-of-life burden.

## 2. Files in this directory

| File | Description |
|------|-------------|
| `pn_qsp_model.dot` / `.svg` / `.png` | Graphviz mechanistic map (113 nodes / 11 clusters) |
| `pn_mrgsolve_model.R` | mrgsolve ODE QSP model (16 disease PD ODEs + 11 PK compartments, 8 scenarios) |
| `pn_shiny_app.R` | Shiny dashboard (8 tabs) |
| `pn_references.md` | 35+ PubMed-linked references |
| `README.md` | This file |

## 3. Mechanistic map clusters

| # | Cluster | Key nodes |
|---|---------|-----------|
| ① | Etiology & Predisposing Factors | Atopic diathesis · CKD-aP · cholestasis · psychiatric comorbidity · filaggrin variants · female/older age · Black race disparity · HIV · DM · neuropathic overlap |
| ② | Epidermal Barrier & Alarmins | TEWL↑ · filaggrin/loricrin↓ · lipid lamellae disorder · S. aureus dysbiosis · TSLP/IL-33/IL-25 · hyperkeratosis · acanthosis · NGF↑ · Semaphorin 3A↓ · PGD2/E2 |
| ③ | Peripheral Sensory Nerve & Itch Transduction | C-fiber pruriceptor · TRPV1/TRPA1 · IL-31RA/OSMRβ · H1R · MRGPRX2/4 · PAR-2 · IENFD (early↑/late↓) · substance P · CGRP · mast cell · DRG · alloknesis |
| ④ | Type-2 (Th2) Immune Axis | Naive CD4 → Th2(GATA3) · IL-4/IL-13/IL-5/IL-31 · ILC2 · eosinophils · IgE-mast cell · IL-4Rα · STAT6 |
| ⑤ | Th17/Th22 Axis | Th17/Th22 differentiation · IL-17A · IL-22 · IFN-γ · neutrophils · keratinocyte proliferation · antimicrobial peptides |
| ⑥ | Central Sensitization & CNS Processing | Spinal GRPR+/Nppb+ neurons · gastrin-releasing peptide · central sensitization · thalamus · cortex (itch percept) · reward circuit · HPA-stress axis · stigma feedback |
| ⑦ | Opioid Receptor Balance | µ-opioid receptor (pruritogenic) · κ-opioid receptor (antipruritic) · endogenous opioid tone · MOR:KOR ratio · difelikefalin target |
| ⑧ | Itch-Scratch Cycle & Nodule Fibrosis | Itch perception · scratch behavior · mechanical trauma · alloknesis · pain · fibroblast activation · TGF-β · collagen · nodule · lichenification · secondary infection · barrier feedback |
| ⑨ | Drug PK/PD — Biologics | Dupilumab (anti-IL-4Rα) · nemolizumab (anti-IL-31RA) · vixarelimab (anti-OSMRβ) · immunogenicity/ADA |
| ⑩ | Drug PK/PD — Small Molecules | Gabapentin/pregabalin (α2δ) · abrocitinib (JAK1i) · topical corticosteroid/calcineurin inhibitor · antihistamine · thalidomide · naltrexone · nalbuphine ER (MOR antag./KOR agonist) · SSRI |
| ⑪ | Clinical Endpoints & QoL | Worst-Itch NRS · PN-IGA · nodule count · sleep disturbance · DLQI · HADS · composite responder · IL-31/TARC biomarker · flare rate |

Total: **113 nodes**, **11 clusters**, **~140 edges** (causal / feedback / drug-inhibition).

## 4. mrgsolve model — ODE structure

* **PK compartments (11):** dupilumab (SC depot + central), nemolizumab
  (SC depot + central), gabapentin (GI depot + central), nalbuphine ER
  (GI depot + central + peripheral, 2-cpt), oral JAK1 inhibitor
  (GI depot + central).
* **Disease/PD compartments (16):** barrier dysfunction, TSLP, Th2
  composite, Th17 composite, IL-31, peripheral sensitization, IENFD,
  central sensitization, opioid tone, scratch intensity, nodule burden,
  cumulative scratch exposure, Worst-Itch NRS, PN-IGA, sleep disturbance,
  DLQI.
* **Time horizon:** 24 weeks (4032 h) by default, adjustable to 52 weeks.
* **8 scenarios:** natural history (topical standard-of-care only),
  dupilumab (600 mg LD / 300 mg Q2W SC), nemolizumab (60 mg LD / 30 mg
  Q4W SC), gabapentin 300 mg TID, nalbuphine ER 162 mg BID, abrocitinib
  200 mg QD (off-label), dupilumab + gabapentin combination, and
  nalbuphine ER + adjunct peripheral KOR-agonist for CKD-associated
  pruritus overlap.

## 5. Drug PK/PD calibration (key)

| Drug | Dose | Disposition | EC₅₀ | Calibration source |
|------|------|-------------|------|--------------------|
| Dupilumab | 600 mg LD / 300 mg Q2W SC | CL 0.0021 L/h·kg · V 4.6 L · F 0.64 | 50 mg/L (IL-4Rα blockade) | Yosipovitch 2023 *Nat Med* (PRIME/PRIME2) |
| Nemolizumab | 60 mg LD / 30 mg Q4W SC | CL 0.014 L/h · V 4.6 L · F 0.68 | 2.5 mg/L (IL-31RA blockade) | Ständer 2020 *NEJM*; Kwatra 2023 *NEJM* (OLYMPIA 2) |
| Gabapentin | 300 mg TID oral | CL 8.5 L/h · V 58 L · F 0.55 | 4 mg/L (central damping) | Gooding 2010 systematic review |
| Nalbuphine ER | 162 mg BID oral | CL 115 L/h · V2 180 L/V3 260 L · F 0.12 | 15 ng/mL (opioid correction) | Weisshaar 2022 *JEADV*; Eudy-Byrne 2023 *Br J Clin Pharmacol* |
| Abrocitinib (off-label) | 200 mg QD oral | CL 42 L/h · V 100 L · F 0.65 | 150 ng/mL (JAK1 blockade) | Case series / extrapolated AD PK |

## 6. Clinical-trial anchors

* **PRIME / PRIME2** (Yosipovitch et al. 2023 *Nat Med*) — dupilumab
  Phase 3 in PN: ~60% achieved ≥4-point WI-NRS improvement at week 24
  vs ~18% placebo; IGA success ~48% vs 18%.
* **Ständer 2020** *NEJM* — nemolizumab Phase 2 in PN: rapid itch
  reduction from week 4, larger effect size than most PN agents to date.
* **OLYMPIA 2** (Kwatra/Ständer 2023 *NEJM*) and **OLYMPIA 1** (Ständer
  2025 *JAMA Dermatol*) — nemolizumab Phase 3: ~56-60% WI-NRS responders
  at week 16 vs ~21-27% placebo.
* **Weisshaar 2022** *JEADV* — oral nalbuphine ER Phase 2 + open-label
  extension: ~44% achieved ≥30% pruritus reduction at 162 mg BID by
  week 10 vs ~36% placebo, with further gains during the extension.
* **KALM-1** (Fishbane 2020 *NEJM*) — difelikefalin (IV KOR agonist)
  in hemodialysis pruritus; mechanistic anchor for the opioid-tone
  correction term used in the CKD-aP overlap scenario.

## 7. Shiny dashboard (8 tabs)

1. **Patient Profile** — baseline itch/nodule trajectory & summary table.
2. **Drug PK** — concentration curves for all 5 active agents.
3. **Barrier / Alarmin / Th2-Th17** — TSLP, Th2, Th17, IL-31 cascade.
4. **Nerve Sensitization** — peripheral sensitization, IENFD, central
   sensitization, and opioid (MOR:KOR) tone.
5. **Itch-Scratch-Fibrosis Cycle** — scratch intensity and nodule burden.
6. **Clinical Endpoints** — WI-NRS, PN-IGA, sleep, DLQI trajectories.
7. **Scenario Comparison** — overlay of WI-NRS, IGA, and nodule count
   across all 8 protocols.
8. **Biomarker Summary** — normalized heat-map of all state variables
   across the simulation horizon.

## 8. Running

```r
# Render map
dot -Tsvg pn_qsp_model.dot -o pn_qsp_model.svg
dot -Tpng -Gdpi=150 pn_qsp_model.dot -o pn_qsp_model.png

# Run model from R
source("pn_mrgsolve_model.R")
head(pn_sim_all)

# Launch dashboard
shiny::runApp("pn_shiny_app.R")
```

## 9. Limitations

This QSP model is **pedagogical and illustrative**: cytokines are
aggregated into composite Th2/Th17 signals rather than individually
resolved, nodule counts are treated as a continuous fibrotic-burden
variable rather than discrete lesion-level state, and opioid receptor
tone is a simplified net-effect term rather than a full receptor
occupancy model. Parameters are order-of-magnitude estimates
cross-referenced to published trial effect sizes, not a fitted
population-PK/PD model — quantitative predictions should not be used
for clinical decision-making.
