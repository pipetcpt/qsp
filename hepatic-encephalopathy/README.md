# Hepatic Encephalopathy (HE) — QSP Model

> **Gut-Liver-Muscle-Brain Ammonia Axis · NH3-Glutamine-Glutamate · Astrocyte Swelling · Neuroinflammation Synergy · West Haven 0–IV**

A comprehensive Quantitative Systems Pharmacology model for **Hepatic Encephalopathy (HE)** — a neuropsychiatric syndrome arising from acute or chronic liver failure (cirrhosis · TIPS shunt · ALF) characterized by hyperammonemia, astrocyte osmotic swelling (Norenberg's glutamine "Trojan horse"), GABA-A hyperactivation by neurosteroids, manganese basal-ganglia accumulation, and systemic-inflammation synergy (Shawcross 2004) that converge on the West Haven 0–IV continuum and 30/90-day mortality.

---

## Mechanism (one-paragraph)

Liver dysfunction (cirrhosis MELD↑, ALF, portosystemic shunting after TIPS) reduces periportal **CPS1-driven ureagenesis** so colonic urease+ bacteria (Enterobacteriaceae, *Streptococcus salivarius*, *Klebsiella*) generate NH3 that bypasses the failing liver via the **shunt fraction (fps)** and reaches plasma at >50 µmol/L. Non-ionized NH3 diffuses across the BBB and is trapped by **astrocyte glutamine synthetase (GS)**, producing intracellular glutamine that acts as a "**Trojan horse**" osmolyte (Albrecht & Norenberg 2006) → mitochondrial permeability transition → ROS → low-grade cerebral edema. NH3 also stimulates **TSPO/PBR**, releasing the neurosteroids **allopregnanolone and THDOC** that potentiate **GABA-A α1γ2** receptors, while LPS translocation from a leaky gut amplifies neurotoxicity through TNF-α/IL-1β/IL-6 (NH3-inflammation **synergy**, Shawcross 2004). Decreased Fischer ratio (**BCAA/AAA**) lets aromatic amino acids cross via LAT-1 and generate false neurotransmitters; manganese deposits in the globus pallidus (T1 hyperintensity) impair dopaminergic tone (parkinsonism). The integrated insult presents as the **West Haven grade 0–IV** continuum (PHES, CFF, ICT, asterixis, GCS) with hospitalization, driving-test failure, and a 30-day mortality HR ≈2–4.

---

## Deliverables

| File | Description |
|------|-------------|
| [`he_qsp_model.dot`](he_qsp_model.dot) | Mechanistic map (12 clusters · 150+ nodes) |
| [`he_qsp_model.svg`](he_qsp_model.svg) | Vector rendering |
| [`he_qsp_model.png`](he_qsp_model.png) | Raster rendering (150 dpi) |
| [`he_mrgsolve_model.R`](he_mrgsolve_model.R) | 25-compartment ODE QSP model (mrgsolve) |
| [`he_shiny_app.R`](he_shiny_app.R) | Shiny dashboard (8 tabs) |
| [`he_references.md`](he_references.md) | 62 PubMed-linked references |

---

## Drug PK/PD (10 agents)

| # | Drug | Mechanism | Dose | Trial anchor |
|---|------|-----------|------|--------------|
| 1 | **Lactulose** | Non-absorbable disaccharide → ↓gut pH (5.5) → NH3→NH4+ trap; laxative | 30–90 g/d PO | Sharma 2009, Gluud 2016 meta |
| 2 | **Rifaximin** | Rifamycin (F<0.4%) → ↓urease+ Enterobacteriaceae | 550 mg BID | **Bass 2010 NEJM RFHE3001** (22% vs 46% breakthrough HE) |
| 3 | **LOLA** | L-ornithine L-aspartate → ↑urea cycle (CPS1) + ↑muscle GS | 20 g IV/d | **Kircheis 1997 Hepatology** (NH3 −38%) |
| 4 | **BCAA** | Leu/Ile/Val → Fischer ratio, muscle ammonia detox | 12 g TID | **Les 2011 AJG** (HE −44%) |
| 5 | **Na-Benzoate** | Alternative N excretion via hippurate | 5–10 g TID | Sushma 1992 |
| 6 | **Glycerol Phenylbutyrate** | PAGN excretion (HPN-100) | 6 mL TID | **Rockey 2014 Hepatology** |
| 7 | **PEG 3350** | Rapid catharsis | 4 L | **Rahimi 2014 JAMA IM (HELP)** |
| 8 | **Probiotic / FMT** | Microbiome shift (Bacteroidetes ↑) | VSL#3 450B / FMT | **Bajaj 2019 PROFIT** |
| 9 | **Albumin 20%** | Anti-inflammatory, LPS binding | 50 g IV | Caraceni 2018 ANSWER |
| 10 | **Flumazenil** | GABA-A competitive antagonist (rescue) | 0.5–1 mg IV | Goulenok 2002 meta |

---

## 25 ODE Compartments

Gut/portal/systemic NH3 · plasma glutamine · brain NH3 · brain glutamine (osmolyte) · astrocyte swelling · plasma LPS · TNFα · GABA-A PAM tone · plasma BCAA/AAA · brain Mn deposit · LOLA PK · BCAA PK · benzoate PK · GPB PK · albumin PK · flumazenil PK · probiotic effect decay · West Haven continuous surrogate · cumulative mortality hazard.

## 9 Treatment Scenarios

S0 untreated · S1 lactulose · S2 rifaximin · S3 lactulose+rifaximin (Sharma 2013) · S4 +LOLA · S5 +BCAA · S6 +Na-Benzoate · S7 +Probiotic/FMT · S8 ALF protocol (Albumin + Flumazenil rescue, fps=0.6).

## Shiny app (8 tabs)

1. Patient Profile (MELD, shunt, sarcopenia, hepatocyte mass)
2. Drug PK (10 agents)
3. Gut → Portal → Systemic NH3 kinetics
4. Brain NH3 / Glutamine / Astrocyte swelling
5. Neuroinflammation (LPS · TNFα · GABA-A tone · Mn)
6. Clinical Endpoints (West Haven, PHES surrogate, mortality hazard)
7. Scenario comparison (9 ladders)
8. Biomarkers (Fischer ratio, plasma/brain glutamine, trial overlays)

---

## Reproduce

```bash
# Mechanistic map
dot -Tsvg he_qsp_model.dot -o he_qsp_model.svg
dot -Tpng -Gdpi=150 he_qsp_model.dot -o he_qsp_model.png

# R / mrgsolve
Rscript -e "source('he_mrgsolve_model.R')"

# Shiny
Rscript -e "shiny::runApp('he_shiny_app.R', launch.browser=TRUE)"
```

## Key references (selected)

- Vilstrup 2014 AASLD/EASL guideline (PMID 25042402)
- Montagnese 2022 EASL update (PMID 35724930)
- Bass 2010 NEJM Rifaximin RFHE3001 (PMID 20335583)
- Sharma 2013 AJG Lactulose+Rifaximin (PMID 23877348)
- Kircheis 1997 Hepatology LOLA pivotal RCT (PMID 9185752)
- Les 2011 AJG BCAA (PMID 21326220)
- Rahimi 2014 JAMA IM HELP — PEG (PMID 25243839)
- Bajaj 2017/2019 Hepatology FMT (PMID 28586116, 31038755)
- Shawcross 2004 J Hepatol NH3-inflammation synergy (PMID 14739095)
- Albrecht & Norenberg 2006 Hepatology Trojan-horse glutamine (PMID 17006913)
- Häussinger 2000 J Hepatol astrocyte swelling (PMID 10898326)
- Mullen 1990 Lancet endogenous BZD (PMID 1975326)
- Olde Damink 2003 Hepatology renal NH3 (PMID 12774005)
- Tandon 2012 Liver Transpl sarcopenia (PMID 22740290)
