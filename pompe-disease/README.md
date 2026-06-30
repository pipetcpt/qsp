# Pompe Disease (GSDII) — QSP Model

> Glycogen storage disease type II · acid α-glucosidase (GAA) deficiency
> Infantile-onset (IOPD) cardiomyopathy + late-onset (LOPD) limb-girdle/diaphragmatic myopathy
> ERT (alglucosidase / avalglucosidase / cipaglucosidase + miglustat) · AAV9-GAA gene therapy

![Mechanistic map](pompe_qsp_model.png)

[🗺️ Full SVG](pompe_qsp_model.svg) · [📐 DOT source](pompe_qsp_model.dot) · [⚙️ mrgsolve](pompe_mrgsolve_model.R) · [💻 Shiny](pompe_shiny_app.R) · [📚 References](pompe_references.md)

## 1. Disease snapshot

* **Genetics.** Autosomal recessive *GAA* (17q25.3) mutations. >600 reported variants.
  Common LOPD allele: `c.-32-13T>G`. CRIM-negative IOPD when both alleles abolish protein.
* **Pathogenesis.** GAA deficiency → lysosomal accumulation of α-1,4 and α-1,6 glycogen →
  enlarged lysosomes (LAMP1+) → autophagic block (LC3-II, p62) → cytoplasmic glycogen
  spillover, mitochondrial dysfunction, sarcomere disruption, satellite-cell exhaustion.
* **Tissue tropism.** IOPD: cardiomyocyte storage drives hypertrophic cardiomyopathy,
  short PR/tall QRS, LVOT obstruction. LOPD: diaphragm, paraspinals, proximal limb girdle,
  with progressive respiratory failure.

## 2. Model deliverables

| File | Description |
|------|-------------|
| `pompe_qsp_model.dot` / `.svg` / `.png` | 12-cluster mechanistic map (≥100 nodes) covering genetics, lysosomal biology, multi-organ pathology, ERT, gene therapy, ADA, biomarkers, endpoints |
| `pompe_mrgsolve_model.R` | 22-compartment ODE QSP model with PK for alglucosidase/avalglucosidase/cipaglucosidase + miglustat + AAV9-GAA + rituximab ITI |
| `pompe_shiny_app.R` | 8-tab Shiny dashboard: patient · PK · PD · clinical endpoints · scenario comparison · biomarkers · ADA · references |
| `pompe_references.md` | 80 curated references with PubMed links |

## 3. mrgsolve compartments

Drug PK (8): `ALGLU_C`, `ALGLU_P`, `AVAL_C`, `AVAL_P`, `CIPA_C`, `CIPA_P`, `MIG_A`, `MIG_C`
Tissue enzyme (3): `GAA_M`, `GAA_C`, `GAA_D`
Lysosomal glycogen (3): `GLYC_M`, `GLYC_C`, `GLYC_D`
Biomarker / immune (3): `HEX4`, `ADA_T`, `RTX_C`
Physiology (4): `LVMI`, `MM_IDX`, `DIAPH_F`, `FVC_UP`
Gene therapy (1): `AAV_X`

Captured derived endpoints include 6-minute walk distance, ventilator-failure hazard,
SF-36 PCS, NT-proBNP, LV ejection fraction, serum CK, and ADA-mediated uptake block.

## 4. Treatment scenarios shipped

| ID | Description |
|----|-------------|
| `no_tx` | Untreated natural history (LOPD or IOPD switchable) |
| `alglu` | Alglucosidase alfa 20 mg/kg IV q2w (Myozyme/Lumizyme) |
| `aval` | Avalglucosidase alfa 20 mg/kg IV q2w (Nexviazyme, COMET) |
| `cipa_mig` | Cipaglucosidase alfa 20 mg/kg IV q2w + Miglustat 195 mg PO (PROPEL) |
| `aav_gt` | One-time AAV9-GAA gene therapy bolus |
| `alglu_iti` | Alglucosidase + Rituximab-based ITI (CRIM-negative IOPD prophylaxis) |

## 5. Calibration anchors

* **Alglucosidase alfa PK** — CL ≈ 0.27 mL/min/kg; V_ss ≈ 100 mL/kg (Hahn 2008).
* **Avalglucosidase alfa** — ~15-fold higher M6P content, ~3-fold higher CI-MPR-mediated
  uptake; COMET trial: +2.4% FVC vs alglucosidase at 49 weeks (Diaz-Manera 2021).
* **Cipaglucosidase + Miglustat** — PROPEL trial: +14 m 6MWT vs alglucosidase (Schoser 2021).
* **IOPD ERT** — 1-year ventilator-free survival ≈ 88% vs <12-month mortality untreated
  (Kishnani 2007, 2009).
* **CRIM-negative ADA** — HSAT ≥51,200 → loss of clinical response; ITI restores tolerance
  (Banugaria 2011; Messinger 2012; Kazi 2017).
* **Hex4 biomarker** — declines on effective ERT (Young 2009; An 2005).
* **AAV9-GAA** — preclinical and Phase I/II results from SPK-3006, AT845, and AAVrh74-MyoAAV
  programs.

## 6. Quick start

```r
# Render the mechanistic map
system("dot -Tsvg pompe_qsp_model.dot -o pompe_qsp_model.svg")

# Simulate a 3-year LOPD adult on avalglucosidase
source("pompe_mrgsolve_model.R")
out <- pompe_run("aval", iopd = FALSE, years = 3)
plot(out$time, out$FVC_UP, type = "l", xlab = "Days", ylab = "FVC upright (%)")

# Launch the Shiny dashboard
shiny::runApp("pompe_shiny_app.R")
```

## 7. Limitations

* Parameter values are **illustrative**, anchored to published ranges but **not** validated
  against individual-patient data. Use as a teaching / hypothesis-generation tool.
* CNS involvement (white-matter changes, cognitive decline in long-term ERT-treated IOPD)
  is captured qualitatively in the map but not explicitly modelled in the ODE system.
* The AAV gene-therapy module is a deliberately simple two-state (induction + decay)
  abstraction; real PD requires capsid- and route-specific PBPK plus immune dynamics.
* Anti-drug-antibody dynamics treat ADAs as a single titre; isotype, epitope spreading,
  and FcRn recycling are not resolved.

See `pompe_references.md` for the full source list.
