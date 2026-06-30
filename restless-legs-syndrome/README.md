# Restless Legs Syndrome (RLS / Willis-Ekbom Disease) — QSP Model

Comprehensive Quantitative Systems Pharmacology (QSP) model integrating brain-iron deficiency, dopaminergic, adenosinergic (A1R), α2δ, opioid, and glutamatergic signaling that drive sensorimotor symptoms and periodic limb movements in sleep (PLMS) of restless legs syndrome — and the response to first- and second-line pharmacotherapies.

> 분류 · **만성 신경계** | 약자 · **RLS / WED** | 최근 갱신 · 2026-06-30

## Quick links

| File | Purpose |
|------|---------|
| [`rls_qsp_model.dot`](rls_qsp_model.dot) · [`.svg`](rls_qsp_model.svg) · [`.png`](rls_qsp_model.png) | 100+ node mechanistic map (12 subgraph clusters) |
| [`rls_mrgsolve_model.R`](rls_mrgsolve_model.R) | 19-compartment ODE / PK-PD QSP model |
| [`rls_shiny_app.R`](rls_shiny_app.R) | 8-tab interactive dashboard |
| [`rls_references.md`](rls_references.md) | 50+ PubMed-indexed references |

## Pathophysiology (mechanistic summary)

```
GWAS risk loci (BTBD9, MEIS1, MAP2K5/SKOR1, PTPRD, TOX3)
    ↓ ↘
Low serum / CSF ferritin → ↓ brain iron in SN, putamen, thalamus  (Connor 2003 Brain)
    ↓
↓ tyrosine hydroxylase activity (Fe²⁺ cofactor) + paradoxical ↑ DA synthesis with ↓ post-syn D2/D3
    ↓
↓ A1 adenosine tone (Ferré 2018) → disinhibited glutamate + dopamine release at striatum/spine
    ↓
A11 → spinal dorsal horn dopamine tone drops at night → sensorimotor hyperexcitability
    ↓
Periodic limb movements in sleep (PLMS) + urge-to-move + IRLS ↑
    ↓
Augmentation hazard rises with chronic high-dose DA agonist exposure
```

## Model components

### 1. Mechanistic map (`.dot`)

- **12 clusters**: genetics · iron · dopamine · adenosine · glutamate/GABA · opioid · spinal cord/PNS · circadian/sleep · augmentation · drug PK · pharmacotherapies · clinical endpoints
- **>100 nodes** with edge annotations citing primary literature (Earley, Connor, Allen, Ferré, Trenkwalder)
- All current and emerging drug classes mapped to their molecular targets

### 2. mrgsolve ODE model (`.R`)

- **19 ODE compartments**: 6 drug-class PKs + IV iron PK + ferritin + brain iron + DA augmentation + IRLS + PLMS
- **6 drug classes** included: pramipexole, ropinirole, rotigotine (patch), gabapentin enacarbil, pregabalin, oxycodone/naloxone PR, IV ferric carboxymaltose
- **8 pre-built treatment scenarios**: untreated · pramipexole 0.25/0.5 mg · rotigotine 2 mg · gabapentin enacarbil 600 mg · pregabalin 300 mg · oxycodone/naloxone PR · IV FCM · augmentation scenario (1 yr)
- Calibration anchors: Allen 2003 (IRLS validation), Trenkwalder 2008 (rotigotine), Garcia-Borreguero 2010 (gabapentin enacarbil PIVOT-RLS), Allen 2014 (pregabalin head-to-head vs pramipexole), Trenkwalder 2013 (oxycodone-naloxone RELOXYN), IRLSSG augmentation criteria (Garcia-Borreguero 2016)
- Patient covariates: weight, eGFR (gabapentinoid CL), sex, pregnancy, ESRD, SSRI trigger, daily caffeine

### 3. Shiny app (`.R`)

8 tabs:
1. **Patient profile** — covariate inputs & risk snapshot
2. **Drug PK** — concentration-time curves
3. **Iron / ferritin** — serum ferritin & brain iron trajectories
4. **Network tones** — DA · A1 · α2δ · MOR effects
5. **IRLS / PLMS / Sleep** — symptom dynamics
6. **Augmentation** — cumulative hazard
7. **Scenario comparison** — side-by-side IRLS / endpoints
8. **Biomarkers & QoL** — CGI-I, ICD hazard, constipation surrogate

### 4. References (`.md`)

50+ PubMed-indexed references grouped by:
- Epidemiology & IRLSSG criteria
- GWAS / risk loci
- Brain iron / ferritin
- Dopaminergic & adenosinergic pathophysiology
- Drug-specific RCTs (DA agonists, α2δ, opioids, IV iron)
- Augmentation / safety / impulse-control disorders
- Special populations (ESRD, pregnancy)
- Pharmacokinetic primary sources

## Clinical / regulatory anchors

| Trial | Drug | Endpoint | Effect |
|-------|------|----------|--------|
| Trenkwalder 2008 *Lancet Neurol* | Rotigotine 3 mg/24h patch | ΔIRLS | −13.7 vs −7.4 placebo @ 6 mo |
| Garcia-Borreguero 2010 (PIVOT-RLS) | Gabapentin enacarbil 1200 mg | ΔIRLS | −13.2 vs −8.8 @ 12 wk |
| Allen 2014 *NEJM* | Pregabalin 300 mg vs pramipexole 0.5 mg | ΔIRLS @ 12 wk · augmentation @ 52 wk | −14.6 vs −12.7 · 1.7% vs 7.7% |
| Trenkwalder 2013 (RELOXYN) | Oxycodone/Naloxone PR | ΔIRLS | −16.5 vs −9.4 @ 12 wk |
| Allen 2011 | IV FCM 1000 mg | IRLS responders | >50% vs ~30% placebo |

## How to run

```r
# (1) Compile the model
source("rls_mrgsolve_model.R")

# (2) Launch the Shiny dashboard (R ≥ 4.3, mrgsolve ≥ 1.5, shiny, DT)
shiny::runApp("rls_shiny_app.R")
```

## Notes / limitations

- Brain iron pool is modeled as a phenomenological index (au, 100 = healthy reference); volume-of-interest MRI R2* could be substituted.
- The augmentation hazard accumulates monotonically with effective DA-agonist tone above a threshold; clinical augmentation is multifactorial and individual susceptibility varies.
- A1 adenosine tone is included for emerging dipyridamole/A1 PAM strategies but is not yet a routine therapeutic target.
- Sleep efficiency is a linear surrogate; full polysomnographic architecture (REM/N3 staging) is not simulated.

## License

Same as the parent repository (MIT) — see [`../LICENSE`](../LICENSE).
