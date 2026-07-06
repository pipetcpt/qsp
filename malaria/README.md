# Malaria QSP Model

**Disease**: Malaria (*Plasmodium falciparum* / *P. vivax* infection)
**Category**: Infectious Disease — Vector-Borne Parasitic
**Date Added**: 2026-07-06
**Model Abbreviation**: `mal`

---

## Mechanistic Map

[![Malaria QSP Mechanistic Map](mal_qsp_model.png)](mal_qsp_model.svg)

> 11 subgraph clusters · 152 nodes · cross-cluster regulation edges

---

## Overview

Malaria is caused by *Plasmodium* spp. protozoa transmitted by female *Anopheles*
mosquitoes. *P. falciparum* causes the great majority of severe disease and death
(cerebral malaria, severe anemia, multi-organ failure) via cytoadherent
sequestration of mature-stage parasites in the microvasculature. *P. vivax* and
*P. ovale* form dormant hepatic **hypnozoites** that cause relapse weeks to
months after apparent cure unless an 8-aminoquinoline **radical cure** is given.
Artemisinin-based combination therapy (ACT) is first-line treatment worldwide;
emerging partial artemisinin resistance (K13 propeller mutations) threatens
control in the Greater Mekong Subregion and, increasingly, East Africa.

### Key Pathophysiology

| Process | Mechanism | Clinical Consequence |
|---------|-----------|----------------------|
| Liver Stage | Sporozoite → hepatocyte → exo-erythrocytic schizogony (6–16d) | Pre-patent period; hypnozoites (P. vivax/ovale) → relapse |
| Erythrocytic Cycle | 48h ring→trophozoite→schizont→rupture (P. falciparum) | Cyclic fever (paroxysm), exponential parasitemia rise |
| Sequestration | PfEMP1 (var genes) binds ICAM-1/CD36/EPCR/CSA | Cerebral malaria, placental malaria; parasites invisible on smear |
| Hemolysis | Schizont rupture destroys host RBC + splenic clearance | Severe anemia, blackwater fever (hemoglobinuria) |
| Cytokine Storm | Hemozoin/GPI-PAMPs → TLR2/9 → TNF-α/IL-6/IL-1β | Fever paroxysm, hypoglycemia, acidosis |
| Naturally-Acquired Immunity | Repeated exposure → anti-PfEMP1/MSP1/AMA1 antibodies | Premunition (asymptomatic parasitemia) in endemic adults |

---

## Antimalarial Drug Classes

| Class | Mechanism | Key Drugs | Notable PK |
|-------|-----------|-----------|------------|
| Endoperoxide (artemisinin) | Heme/Fe²⁺-activated free-radical alkylation; broad, dominant ring-stage kill | Artesunate → DHA | t½ ~0.5–1h; PRR ~10⁴/cycle |
| Aryl-amino alcohol | Heme detoxification inhibition | Lumefantrine, Mefloquine | t½ 3–21 days (post-treatment prophylaxis) |
| Bisquinoline | Hemozoin formation inhibition | Piperaquine | t½ 20–33 days |
| 4-Aminoquinoline | Heme polymerization inhibition | Amodiaquine → Desethyl-AQ | t½ (DEAQ) ~9–18 days |
| 8-Aminoquinoline | Hypnozoiticidal + gametocytocidal (CYP2D6-dependent) | Primaquine, Tafenoquine | Oxidative hemolysis risk in G6PD deficiency |

---

## Recommended WHO-Guideline Regimens

| Regimen | Use Case | Advantages | Limitations |
|---------|----------|-----------|--------------|
| Artemether-Lumefantrine (AL) | First-line uncomplicated P.f. | Widely available, well-tolerated | Twice-daily × 3d, food required |
| Artesunate-Amodiaquine (ASAQ) | First-line (Africa) | Once daily | QT risk, AQ-related neutropenia (rare) |
| Dihydroartemisinin-Piperaquine (DP) | First-line, longest prophylaxis | Once daily × 3d | QTc prolongation |
| IV Artesunate | Severe malaria (WHO preferred over quinine) | Superior mortality vs quinine (AQUAMAT/SEAQUAMAT) | Requires ≥24h parenteral therapy then oral ACT |
| ACT + Primaquine (14d) | P. vivax/ovale radical cure | Prevents relapse | Requires G6PD testing |

---

## mrgsolve Model Structure

### Compartments (22)

| # | Compartment | Description | Units |
|---|-------------|-------------|-------|
| 1 | AS_GUT | Artesunate gut depot | mg |
| 2 | AS_PLASMA | Artesunate plasma | mg/L |
| 3 | DHA_PLASMA | Dihydroartemisinin (active metabolite) | mg/L |
| 4 | LUM_GUT | Lumefantrine gut depot | mg |
| 5 | LUM_CENTRAL | Lumefantrine central | mg/L |
| 6 | LUM_PERIPH | Lumefantrine peripheral | mg/L |
| 7 | PPQ_GUT | Piperaquine gut depot | mg |
| 8 | PPQ_CENTRAL | Piperaquine central | mg/L |
| 9 | PPQ_PERIPH | Piperaquine peripheral | mg/L |
| 10 | AQ_GUT | Amodiaquine gut depot | mg |
| 11 | DEAQ_PLASMA | Desethylamodiaquine (active metabolite) | mg/L |
| 12 | PQ_GUT | Primaquine gut depot | mg |
| 13 | PQ_PLASMA | Primaquine plasma | mg/L |
| 14 | RBC_U | Uninfected RBC pool | cells/µL |
| 15 | PRBC_RING | Ring-stage infected RBC | parasites/µL |
| 16 | PRBC_TROPH | Trophozoite-stage infected RBC | parasites/µL |
| 17 | PRBC_SCHIZONT | Schizont-stage infected RBC (sequestered) | parasites/µL |
| 18 | LIVER_PARASITE | Liver-stage/hypnozoite burden | arbitrary units |
| 19 | GAMETOCYTE | Mature (Stage V) gametocyte density | gametocytes/µL |
| 20 | HB | Hemoglobin | g/dL |
| 21 | TEMP | Body temperature | °C |
| 22 | IMMUNITY | Acquired immunity index | 0–1 |

### Core ODE System (age-structured erythrocytic cycle)

```
dRing/dt      = new_infections − k_RT·Ring − E_AS·Ring
dTroph/dt     = k_RT·Ring − k_TS·Troph − (0.6·E_AS + E_partner)·Troph − gam_commit·k_RT·Ring
dSchizont/dt  = k_TS·Troph − k_SR·Schizont − 0.8·E_partner·Schizont
new_infections = (Liver_egress + burst·k_SR·Schizont) · inv_eff0·(1−immune_block·IMMUNITY)
E_drug = Emax · C^h / (EC50_eff^h + C^h)          [Hill/Emax]
EC50_AS,eff = EC50_AS · (1 + (Kres_shift−1)·K13_RES)   [artemisinin resistance]
dHB/dt   = −hb_per_rupture·(k_SR·Schizont) + erythropoietic compensation
dTEMP/dt = (37 + fever_drive(rupture flux) − TEMP) · k_temp_decay
```

Calibration notes: artesunate/DHA PK from Morris (2011, *Malar J*) population PK;
lumefantrine and piperaquine PK from Tarning (2012 *CPT*; 2008 *AAC*); ring-stage
killing kinetics from Saralamba (2011, *PNAS*) and the WWARN Parasite Clearance
Estimator methodology; K13/ring-stage-survival resistance parameterization from
Ashley (2014, *NEJM*).

---

## Treatment Scenarios (9)

| # | Scenario | Regimen | Notes |
|---|----------|---------|-------|
| 1 | Untreated, non-immune | — | Natural history, traveler/naive host |
| 2 | Untreated, semi-immune | — | High baseline immunity (premunition) |
| 3 | Artemether-Lumefantrine | AL, 6-dose 3-day | WHO first-line |
| 4 | Artesunate-Amodiaquine | ASAQ, 3-day | First-line (Africa) |
| 5 | Dihydroartemisinin-Piperaquine | DP, 3-day | Longest post-treatment prophylaxis |
| 6 | Severe malaria | IV Artesunate → oral AL | WHO-preferred severe malaria pathway |
| 7 | Artemisinin resistance | K13 C580Y + standard AL | Delayed clearance/recrudescence risk |
| 8 | P. vivax + radical cure | ACT + 14-day Primaquine | Prevents relapse |
| 9 | P. vivax, no radical cure | ACT only | Hypnozoite-driven relapse |

---

## Shiny App Tabs

| Tab | Content |
|-----|---------|
| ① Patient/Infection Profile | Species, baseline parasitemia/immunity, WHO severe-malaria criteria |
| ② Drug PK | DHA, partner-drug (LUM/PPQ/DEAQ), and primaquine plasma concentrations |
| ③ Parasite Dynamics | Peripheral parasitemia, stage composition, gametocyte carriage |
| ④ Clinical Endpoints | Parasite/fever clearance time, Day-28 ACPR outcome |
| ⑤ Scenario Comparison | Multi-regimen parasitemia clearance comparison |
| ⑥ Biomarkers | Hemoglobin, temperature, immunity index, liver/hypnozoite burden |
| ⑦ Resistance (K13) | Wild-type vs K13-mutant clearance kinetics |

---

## Key Clinical Trial Evidence

| Trial | Regimen | n | Key Finding |
|-------|---------|---|-------------|
| SEAQUAMAT (2005) | IV Artesunate vs Quinine (Asia, severe malaria) | 1461 | Mortality 15% vs 22% (RRR 34.7%) |
| AQUAMAT (2010) | IV Artesunate vs Quinine (Africa, pediatric severe malaria) | 5425 | Mortality 8.5% vs 10.9% (RRR 22.5%) |
| Ashley et al. (2014, TRAC) | K13 genotype vs parasite clearance half-life | 1241 | K13 mutations → delayed clearance across GMS |
| WWARN DP pooled analysis | DHA-piperaquine efficacy | >7000 | Day-42 PCR-adjusted efficacy >95% (most sites) |
| CDC/WHO Radical Cure guidance | Primaquine/Tafenoquine + G6PD testing | — | Point-of-care G6PD testing recommended before radical cure |

---

## Files

| File | Description |
|------|-------------|
| `mal_qsp_model.dot` | Graphviz DOT mechanistic map source (11 clusters, 152 nodes) |
| `mal_qsp_model.svg` | SVG rendered mechanistic map |
| `mal_qsp_model.png` | PNG rendered mechanistic map (150 dpi) |
| `mal_mrgsolve_model.R` | mrgsolve ODE model (22 CMT, 9 scenarios) |
| `mal_shiny_app.R` | Shiny dashboard (7 tabs, interactive PK/PD simulation) |
| `mal_references.md` | PubMed references (10 sections) |
| `README.md` | This file |
