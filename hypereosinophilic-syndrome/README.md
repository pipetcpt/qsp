# Hypereosinophilic Syndrome (HES, 과호산구증후군) — QSP Model

> Integrated Quantitative Systems Pharmacology model of hypereosinophilic
> syndrome, linking molecular-genotype-specific eosinophilopoiesis —
> constitutively-active fusion tyrosine kinase signaling in myeloid/
> neoplastic HES (FIP1L1-PDGFRA, PDGFRB-rearranged, JAK2-rearranged) or
> clonal-T-cell IL-5 overproduction in lymphocytic-variant HES — to
> circulating/tissue eosinophil burden, degranulation-mediated organ
> damage (cardiac Loeffler endocarditis as the anchor severity endpoint,
> plus pulmonary/dermatologic/neurologic/GI involvement and a
> hypercoagulable state), and genotype-matched pharmacologic response
> (imatinib, corticosteroids, anti-IL-5-axis biologics [mepolizumab,
> benralizumab], hydroxyurea, ruxolitinib).

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`hes_qsp_model.dot`](hes_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`hes_qsp_model.svg`](hes_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`hes_qsp_model.png`](hes_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`hes_mrgsolve_model.R`](hes_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`hes_shiny_app.R`](hes_shiny_app.R) |
| 📚 References             | [`hes_references.md`](hes_references.md) |

---

## 1. Disease in one paragraph

Hypereosinophilic syndrome is defined by persistent peripheral blood
eosinophilia (≥1500/µL on ≥2 occasions ≥1 month apart) with eosinophil-
mediated organ damage after secondary causes are excluded (Valent 2012
JACI). Roughly 10-20% of cases carry a cryptic del(4)(q12q12) that fuses
**FIP1L1** to **PDGFRA**, producing a constitutively active tyrosine
kinase that drives clonal eosinophil-progenitor expansion via STAT5/
PI3K-AKT/RAS-MAPK signaling and is exquisitely sensitive to low-dose
**imatinib** (Cools 2003 NEJM) — a paradigm subsequently extended to rarer
**PDGFRB** and **FGFR1/JAK2** rearrangements (myeloid/neoplastic HES,
M-HES). A distinct **lymphocytic variant (L-HES)** arises from an
aberrant clonal CD3-CD4+ T-cell population that chronically secretes
IL-5 (and IL-4/IL-13), driving eosinophilopoiesis through the IL-5
receptor/JAK2/STAT5 axis without a kinase-fusion driver (Simon 1999 NEJM);
many cases remain **idiopathic**. Whatever the driver, marrow-derived
eosinophils enter the circulation, home to tissue via the eotaxin/CCR3
axis, and degranulate — releasing major basic protein, eosinophil
cationic protein, eosinophil peroxidase, and eosinophil-derived
neurotoxin that directly damage endothelium and parenchyma. The most
feared consequence is **Loeffler endocarditis**: acute eosinophilic
endocardial necrosis progresses through mural thrombus formation to
endomyocardial fibrosis, entrapping the papillary muscles/chordae to
produce valvular regurgitation and restrictive cardiomyopathy (Parrillo
1979 Am J Med); pulmonary infiltrates, urticaria/angioedema, peripheral
neuropathy, eosinophilic gastroenteritis, and a hypercoagulable
thromboembolic state (including intracardiac-thrombus embolization to
the CNS) complete the multi-organ picture (Ogbogu 2009 JACI). Treatment
is genotype-matched: **imatinib** achieves rapid, durable hematologic and
molecular remission in FIP1L1-PDGFRA+/PDGFRB-rearranged M-HES (resistance
emerging via the T674I kinase-domain mutation); **corticosteroids** remain
first-line for idiopathic/lymphocytic-variant disease via eosinophil
apoptosis induction and IL-5 transcriptional suppression, but relapse on
taper is common; **mepolizumab** (anti-IL-5) and **benralizumab**
(anti-IL-5Rα, afucosylated Fc driving NK-cell ADCC eosinophil depletion)
provide steroid-sparing options for steroid-refractory/dependent disease;
and **hydroxyurea**, **ruxolitinib** (JAK2-rearranged HES), interferon-α,
or alemtuzumab are reserved for refractory or kinase-inhibitor-resistant
cases.

## 2. Mechanistic clusters (13 in the DOT map, 126 nodes)

1. Etiology & molecular classification (diagnostic criteria, FIP1L1-PDGFRA/PDGFRB/FGFR1/JAK2/KIT, L-HES, idiopathic, familial, overlap)
2. Bone marrow eosinophilopoiesis & cytokine drivers (CD34+ progenitor, GATA-1/PU.1/C/EBPα, IL-5/IL-3/GM-CSF, JAK2-STAT5)
3. Myeloid/neoplastic kinase signaling (constitutive FIP1L1-PDGFRA/PDGFRB/FGFR1 TK, STAT5/PI3K-AKT/RAS-MAPK, T674I resistance)
4. Lymphocytic-variant (L-HES) T-cell pathway (CD3-CD4+ clone, TCR rearrangement, Th2 skew, IL-5/IL-4/IL-13, lymphoma-transformation risk)
5. Eosinophil activation & tissue effector mediators (eotaxin/CCR3, adhesion/transmigration, MBP/ECP/EPO/EDN degranulation, ROS, tissue factor)
6. Cardiac involvement — Loeffler endocarditis (necrosis → mural thrombus → fibrosis → valvular regurgitation/restrictive cardiomyopathy)
7. Pulmonary, dermatologic & neurologic involvement (infiltrates/fibrosis, urticaria/angioedema, mononeuritis multiplex, CNS embolism/encephalopathy)
8. GI/hepatosplenic involvement & hypercoagulable state (eosinophilic gastroenteritis, hepatosplenomegaly, platelet activation, VTE risk)
9. Drug PK/PD — imatinib (ATP-competitive PDGFRA/B kinase inhibition, hematologic/molecular response, T674I resistance)
10. Drug PK/PD — corticosteroids (glucocorticoid receptor, NF-κB inhibition, eosinophil apoptosis, IL-5 suppression, taper relapse)
11. Drug PK/PD — anti-IL-5-axis biologics (mepolizumab IL-5 neutralization; benralizumab IL-5Rα blockade + ADCC depletion)
12. Drug PK/PD — hydroxyurea / ruxolitinib / IFN-α / alemtuzumab (nonselective cytoreduction, JAK1/2 inhibition, T-cell immunomodulation/depletion)
13. Clinical endpoints, biomarkers & monitoring (AEC, ECP, IL-5, tryptase, composite symptom/flare score, organ-damage composite, survival)

## 3. mrgsolve model (21 ODE compartments)

* **Drug PK (6 drug classes, 12 compartments)** — imatinib (gut + central),
  prednisone/prednisolone (gut + central), mepolizumab (SC depot +
  central), benralizumab (SC depot + central), hydroxyurea (gut +
  central), ruxolitinib (gut + central).
* **Disease/PD network (9 compartments)** — constitutive kinase-signaling
  index (genotype-dependent baseline, imatinib/ruxolitinib-inhibitable),
  endogenous IL-5 index (corticosteroid-suppressible, mepolizumab-
  neutralizable), marrow eosinophil-production drive, circulating
  absolute eosinophil count (AEC; corticosteroid apoptosis + benralizumab
  ADCC depletion), tissue eosinophil burden, cardiac (Loeffler) damage
  index (slow, partly irreversible fibrotic accrual), cardiac troponin,
  serum ECP, and a composite HES symptom/flare score.
* A `GENOTYPE` flag (0=idiopathic/L-HES, 1=FIP1L1-PDGFRA+, 2=PDGFRB-
  rearranged, 3=JAK2-rearranged) sets the baseline kinase-signaling vs.
  IL-5 drive and gates drug sensitivity exactly as in clinical practice
  (imatinib/ruxolitinib act only on kinase-driven genotypes; anti-IL-5-
  axis biologics and corticosteroids act on the IL-5/apoptosis axis
  regardless of genotype).

### 10 prebuilt scenarios

| # | Scenario | Calibration anchor |
|---|---|---|
| 1 | Untreated FIP1L1-PDGFRA+ M-HES (natural history) | Gotlib 2004 Blood |
| 2 | FIP1L1-PDGFRA+ M-HES + imatinib 100 mg PO QD | Cools 2003 NEJM; Klion 2004 Blood |
| 3 | Imatinib-resistant (T674I) → switch to hydroxyurea | Cools 2004 Cancer Cell |
| 4 | Idiopathic HES + prednisone (steroid-responsive) | Butterfield 2007 Immunol Allergy Clin North Am |
| 5 | Idiopathic HES steroid-refractory + mepolizumab 300 mg SC Q4W | Roufosse 2020 JACI |
| 6 | Lymphocytic-variant HES + benralizumab 30 mg SC Q4W | Kuang 2019 NEJM |
| 7 | PDGFRB-rearranged HES + low-dose imatinib 100 mg QD | David 2007 Blood |
| 8 | Cardiac Loeffler endocarditis + high-dose steroid pulse | Parrillo 1979 Am J Med; Ogbogu 2007 |
| 9 | JAK2-rearranged HES + ruxolitinib 20 mg BID | Reiter & Gotlib 2017 Blood |
| 10 | Steroid taper + mepolizumab steroid-sparing combination | Kuang 2018 JACI Pract |

## 4. Shiny dashboard (8 tabs)

1. **Patient profile** — genotype/mutation flags, age/weight/baseline AEC.
2. **Drug PK** — log-scale plasma concentrations for the six tracked drugs.
3. **Kinase/IL-5/marrow PD** — constitutive kinase-signaling index, endogenous IL-5, marrow production drive.
4. **Eosinophil & tissue burden** — circulating AEC (log scale, diagnostic threshold line) and tissue eosinophil burden.
5. **Cardiac & clinical endpoints** — Loeffler damage index, composite symptom/flare score, flare-event table.
6. **Biomarkers** — serum ECP and cardiac troponin.
7. **Scenario comparison** — runs all 10 scenarios with the chosen profile.
8. **References** — key citations and link to the full bibliography.

## 5. How to run

```bash
# 1) Render the mechanistic map
dot -Tsvg hes_qsp_model.dot -o hes_qsp_model.svg
dot -Tpng -Gdpi=150 hes_qsp_model.dot -o hes_qsp_model.png
```

```r
# 2) Simulate scenarios in R
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny",
                   "shinydashboard","DT"))
library(mrgsolve)
mod <- mread_cache("hes_mrgsolve_model.R")
out <- mod %>% param(GENOTYPE=1) %>%
  mrgsim(ev(amt=100, cmt="GUT_IMA", ii=24, addl=180), end=4320, delta=1)
plot(out, c("AEC","TISSUE","CARDIAC","SYMPTOM"))

# 3) Launch the dashboard
shiny::runApp("hes_shiny_app.R")
```

## 6. Key clinical anchors used during calibration

| Endpoint | Comparator | Expected effect |
|---|---|---|
| Imatinib hematologic response | FIP1L1-PDGFRA+ M-HES, 100-400 mg/d | Complete hematologic response within 1-4 weeks; molecular remission over months (Cools 2003 NEJM; Klion 2004 Blood) |
| T674I resistance | Imatinib EC50 shift | ~10-15x rightward EC50 shift, analogous to BCR-ABL T315I (Cools 2004 Cancer Cell) |
| Corticosteroid response | Idiopathic/L-HES, 0.5-1 mg/kg/d | Rapid AEC suppression via apoptosis + IL-5 transcriptional block; relapse common on taper (Butterfield 2007) |
| Mepolizumab flare reduction | Steroid-dependent/refractory HES, 300 mg Q4W | Reduced flare frequency, enables steroid-sparing (Roufosse 2020 JACI; Kuang 2022 NEJM) |
| Benralizumab depletion | PDGFRA-negative HES, 30 mg Q4-8W | Near-complete eosinophil depletion via ADCC, independent of IL-5 status (Kuang 2019 NEJM) |
| Loeffler endocarditis progression | Uncontrolled eosinophilia | Necrotic → thrombotic → fibrotic stage over weeks-months; fibrotic stage largely irreversible (Parrillo 1979 Am J Med) |

## 7. Caveats

* Designed for **research, education, and hypothesis generation** — not
  clinical decision support, hematology/immunology practice guidance, or
  a substitute for validated organ-damage risk-prediction tools.
* The kinase-signaling and IL-5 indices are lumped arbitrary-unit
  surrogates for genotype-specific molecular drive, not mechanistic
  receptor-binding/transcriptional kinetic models; real FIP1L1-PDGFRA/
  PDGFRB/JAK2 signaling involves additional feedback loops (e.g.,
  mTOR, negative regulators of cytokine signaling) not explicitly
  represented here.
* Cardiac (Loeffler) damage is modeled as a single lumped fibrotic-
  accrual index; it does not resolve the distinct necrotic/thrombotic/
  fibrotic histopathologic stages or chamber-specific (mitral vs.
  tricuspid) valvular mechanics.
* Inter-individual variability (genotype penetrance, organ-damage
  susceptibility, CYP3A4 polymorphism affecting imatinib clearance)
  needs `omega()` blocks for population simulations; all parameters here
  are typical-value point estimates.

## 8. License

Inherits the repository [LICENSE](../LICENSE).
