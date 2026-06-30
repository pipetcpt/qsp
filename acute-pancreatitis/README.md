# Acute Pancreatitis (AP) — QSP Model

> Trypsinogen → trypsin auto-activation · DAMP/TLR4 → NF-κB cytokine storm
> · acinar necrosis · capillary leak · gut translocation · SIRS → MODS
> with Atlanta 2012 severity classification, WATERFALL fluid resuscitation,
> ERCP-PEP indomethacin prophylaxis, protease inhibitors (gabexate · nafamostat ·
> ulinastatin), octreotide, meropenem for infected necrosis, anakinra,
> fentanyl PCA, and early enteral nutrition.

---

## 1. Disease overview

Acute pancreatitis (AP) is an abrupt inflammatory disorder of the exocrine
pancreas, characterised by:

1. **Premature intra-acinar trypsinogen activation** driven by sustained
   cytosolic Ca²⁺ elevation, mitochondrial permeability transition pore
   opening, ER stress, impaired autophagy and zymogen/lysosome
   co-localisation (cathepsin B catalyses trypsinogen → trypsin in
   pathological compartments).
2. **Pancreatic auto-digestion** by trypsin, chymotrypsin, elastase,
   phospholipase A2 (PLA2) and carboxypeptidases.
3. **Local & systemic inflammatory cascade** through NF-κB / AP-1 /
   NLRP3 inflammasome activation, generating TNF-α, IL-1β, IL-6, IL-8,
   PAF, HMGB1 and complement (C3a/C5a).
4. **Microcirculatory failure** — endothelial dysfunction, glycocalyx
   shedding, capillary leak, micro-thrombi, ischaemia–reperfusion.
5. **Gut barrier loss & bacterial translocation** generating endotoxemia
   and infected pancreatic necrosis in the late phase.
6. **SIRS → MODS** — ARDS, AKI, hepatic dysfunction, encephalopathy and
   distributive shock; SOFA / APACHE-II / BISAP track severity.

Atlanta 2012 classifies AP as mild (no organ failure), moderately severe
(transient organ failure ≤48 h) or severe (persistent organ failure
>48 h). Severe AP carries 15–30 % mortality, largely from infected
necrosis and persistent MODS.

## 2. Files in this directory

| File | Description |
|------|-------------|
| `ap_qsp_model.dot` / `.svg` / `.png` | Graphviz mechanistic map (151 nodes / 14 clusters) |
| `ap_mrgsolve_model.R` | mrgsolve ODE QSP model (25 disease ODEs + 11 PK compartments, 10 scenarios) |
| `ap_shiny_app.R` | Shiny dashboard (8 tabs) |
| `ap_references.md` | 80+ PubMed-linked references |
| `README.md` | This file |

## 3. Mechanistic map clusters

| # | Cluster | Key nodes |
|---|---------|-----------|
| ① | Etiology / TIGAR-O | Gallstones · Alcohol · HTG · Post-ERCP · Drugs · Trauma · Hypercalcemia · IPMN · AIP · Idiopathic · PRSS1/SPINK1/CFTR/CTRC/CASR/CLDN2 · Viral |
| ② | Acinar Injury | Ca²⁺ overload · SERCA · IP3R · STIM1/ORAI1 · MPTP · ROS · ER stress (PERK/IRE1/ATF6) · impaired autophagy · zymogen co-localisation |
| ③ | Trypsin Cascade | Trypsinogen → Trypsin (cathepsin B) · SPINK1 brake · CTRC degradation · chymo · elastase · PLA2 · PAR-2 |
| ④ | Inflammation | NF-κB · AP-1 · JNK · TLR4 · NLRP3 · caspase-1 · GSDMD · IL-1β · IL-18 · TNF-α · IL-6 · IL-8 · MCP-1 · PAF · HMGB1 |
| ⑤ | Innate Immunity | Neutrophils · NETs · MPO · M1/M2 macrophages · DAMPs · Complement · Kallikrein–Kinin · Resolvins |
| ⑥ | Death Modes | Apoptosis · Necroptosis (RIPK1/3/MLKL) · Pyroptosis · Ferroptosis · Coagulative necrosis · Fat necrosis (Ca-soap) |
| ⑦ | Vasculature | Endothelium · capillary leak · ICAM/VCAM · thrombin / DIC · ischemia · NO/iNOS · ET-1 · glycocalyx |
| ⑧ | Local Complications | APFC · pseudocyst · ANC · WON · infected necrosis · splenic/portal thrombosis · pseudoaneurysm · colonic necrosis |
| ⑨ | Gut Translocation | barrier · ZO-1/claudin · dysbiosis · BT · LPS · ileus · IAH · ACS |
| ⑩ | SIRS / MODS | SIRS · sepsis · ARDS · AKI · liver · CNS · shock · lactate · SOFA · MODS |
| ⑪ | Pain | visceral nociceptors · TRPV1 · PGE2 · substance-P · CGRP · spinal sensitisation · μ-OR · VAS |
| ⑫ | Drugs | LR · NS · indomethacin · diclofenac · octreotide · gabexate · nafamostat · ulinastatin · meropenem · ciprofloxacin · fentanyl PCA · morphine · PPI · plasmapheresis · ERCP · drainage · EN · pirfenidone · anakinra |
| ⑬ | Clinical Endpoints | lipase · amylase · CRP · APACHE-II · BISAP · Ranson · CTSI · Atlanta · mortality · LOS · recurrence · post-AP DM |
| ⑭ | Recovery / CARS | IL-10 · TGF-β · Treg · CARS · late immunosuppression · exocrine / endocrine insufficiency · regeneration |

Total: **151 nodes**, **14 clusters**, **~150 edges** (causal / inhibition / drug effect).

## 4. mrgsolve model — ODE structure

* **PK compartments (11):** indomethacin (PR depot + central), octreotide
  (SC depot + central), gabexate, nafamostat, ulinastatin, meropenem,
  anakinra (SC depot + central), fentanyl.
* **Disease compartments (25):** trypsinogen pool, active trypsin, PLA2,
  Ca²⁺, ROS, NF-κB, TNF-α, IL-1β, IL-6, IL-8, CRP, neutrophils, necrosis,
  DAMPs, permeability, gut barrier, bacterial translocation, PaO₂/FiO₂,
  creatinine, bilirubin, MAP, GCS, SOFA, VAS, cumulative mortality
  hazard.
* **Time horizon:** 14 days (336 h) with hourly output.
* **10 scenarios:** supportive only, LR aggressive (10 mL/kg/h),
  LR moderate (5 mL/kg/h), indomethacin PR (PEP), octreotide SC,
  gabexate IV, nafamostat infusion, ulinastatin q8h, meropenem for
  infected necrosis, anakinra for severe AP (HTG-induced).

## 5. Drug PK / PD calibration (key)

| Drug | Dose | Disposition | EC₅₀ / Ki | Calibration source |
|------|------|-------------|-----------|--------------------|
| Indomethacin | 100 mg PR ×1 | CL 5 L/h · V 60 L · F 0.9 | EC₅₀ 1.2 mg/L | Elmunzer 2012 NEJM |
| Octreotide   | 100 µg SC q8h | CL 9 L/h · V 20 L | EC₅₀ 1.5 ng/mL | Chanson 1993 / Uhl 1999 |
| Gabexate     | 600 mg IV q6h | CL 30 L/h · V 20 L | EC₅₀ 3.0 mg/L | Chen 2000 |
| Nafamostat   | 1.67 mg/h IV continuous | CL 60 L/h · V 16 L | EC₅₀ 0.5 mg/L | Yoshikawa 1996 |
| Ulinastatin  | 200 000 U IV q8h | CL 5 L/h · V 8 L | EC₅₀ 50 U/L | Tsujino 2005 / Park 2014 |
| Meropenem    | 1 g IV q8h | CL 12 L/h · V 18 L | EC₅₀ 8 mg/L | Buchler 2000 / PROCAP |
| Anakinra     | 100 mg SC q24h | CL 0.4 L/h · V 17 L · F 0.95 | Ki 1.0 mg/L | Akinosoglou 2024 |
| Fentanyl PCA | 50 µg q1h | CL 50 L/h · V 250 L | EC₅₀ 1.5 ng/mL | (analgesia) |
| Lactated Ringer's | 5 mL/kg/h | volume effect (no PK) | — | WATERFALL 2022 |

## 6. Clinical-trial anchors

* **Atlanta 2012** revision (Banks 2013 *Gut*) — severity classification.
* **WATERFALL** (de-Madaria 2022 *NEJM*) — aggressive LR (10 mL/kg/h)
  *not* superior to moderate (5 mL/kg/h), with more fluid overload.
* **Elmunzer 2012** *NEJM* — rectal indomethacin halves post-ERCP
  pancreatitis incidence in high-risk patients.
* **PANTER** (van Santvoort 2010 *NEJM*) and **TENSION** (van Brunschot
  2018 *Lancet*) — endoscopic / surgical step-up drainage for infected
  necrosis is preferred over open necrosectomy.
* **PROCAP** (Dellinger 2007 *Ann Surg*) — meropenem in predicted
  severe necrotising AP did not reduce mortality vs placebo (used here
  for documented infection only).
* **PYTHON** (Bakker 2014 *NEJM*) — early nasoenteric feeding
  equivalent to on-demand oral feeding; supports early EN.
* **Anakinra in SAP** (Akinosoglou 2024 *J Clin Med*) — IL-1Ra
  reduced organ failure duration in severe AP.

## 7. Shiny dashboard (8 tabs)

1. **Patient Profile** — etiology, genetics, baseline TG; hemodynamic
   summary.
2. **Drug PK** — concentration curves for all 8 active agents.
3. **Trypsin / Cytokines** — active trypsin and TNF / IL-1 / IL-6 / IL-8
   / CRP trajectories.
4. **Necrosis & DAMPs** — pancreatic necrosis %, permeability, gut
   barrier and BT.
5. **Organ Failure** — SOFA composite plus PF ratio, creatinine,
   bilirubin, MAP, GCS.
6. **Severity & Survival** — BISAP / SOFA / mortality hazard /
   summary score table.
7. **Scenario Comparison** — overlay of SOFA, necrosis, and survival
   for all 10 scenarios.
8. **Biomarker Heat-map** — peak trypsin / TNF / IL-6 / CRP /
   necrosis / SOFA / BT / VAS across scenarios.

## 8. Running

```r
# Render map
dot -Tsvg ap_qsp_model.dot -o ap_qsp_model.svg
dot -Tpng -Gdpi=150 ap_qsp_model.dot -o ap_qsp_model.png

# Run model from R
source("ap_mrgsolve_model.R")
df <- run_all_scenarios()
head(df)

# Launch dashboard
shiny::runApp("ap_shiny_app.R")
```

## 9. Limitations

The QSP model is **pedagogical and illustrative**: simplified ODE
representation aggregates many parallel pathways (e.g. Ca²⁺ subspecies,
specific complement components, individual neutrophil sub-populations).
Quantitative predictions should be cross-checked against patient-level
data and validated trial cohorts before any clinical interpretation.
