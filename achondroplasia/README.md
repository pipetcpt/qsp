# Achondroplasia (연골무형성증) — QSP Disease Model

> **Achondroplasia (ACH)** is the most common form of disproportionate short
> stature, caused by a **gain-of-function mutation in FGFR3** (G380R in the
> transmembrane domain, ~97% of cases). Constitutive FGFR3 dimerization
> hyperactivates the downstream **RAS-RAF-MEK-ERK (MAPK)** cascade in the
> growth-plate, which is a *negative* regulator of chondrocyte proliferation
> and hypertrophic differentiation — so the mutation suppresses endochondral
> bone growth, producing rhizomelic limb shortening, foramen magnum stenosis,
> spinal stenosis, and characteristic craniofacial features.

The model captures: **FGFR3 G380R** → ligand-independent constitutive
dimerization → RAS-RAF-MEK-**ERK**/STAT1 hyperactivation → growth-plate
chondrocyte cell-cycle arrest (p21^CIP1) + disorganized columns + blunted
hypertrophic differentiation → **growth-plate height ↓** → rhizomelic short
stature, premature cranial-base synchondrosis fusion (**foramen magnum
stenosis**, cervicomedullary compression), short pedicles (**spinal
stenosis**), midface hypoplasia (Eustachian dysfunction → **otitis media**,
**OSA**). The **CNP-NPR2-cGMP-PKGII** axis physiologically counter-regulates
FGFR3-MAPK signaling in the growth plate — this is the molecular target of
**vosoritide** (CNP analog) and **TransCon CNP/navepegritide**; **infigratinib**
instead blocks the FGFR3 kinase domain directly.

---

## Deliverables

| File | Purpose |
|------|---------|
| `acdp_qsp_model.dot` / `.svg` / `.png` | Graphviz mechanistic map — **15 clusters, 160 nodes** |
| `acdp_mrgsolve_model.R` | **20-ODE** mrgsolve model (6 PK + 14 disease/clinical) with 10 scenarios |
| `acdp_shiny_app.R` | Interactive **8-tab** Shiny dashboard |
| `acdp_references.md` | **50** PubMed citations grouped by section |

---

## Mechanistic Map — Cluster Index

1. **Genetics & etiology** — FGFR3 4p16.3, G380R (~97%), de novo (~80%, paternal age), autosomal dominant, homozygous ACH (lethal), HCH↔ACH↔TD spectrum
2. **FGFR3 receptor signaling** — ligand-independent constitutive dimerization, autophosphorylation, FRS2-RAS-RAF1-MEK-ERK, STAT1, PI3K-AKT, blunted SOCS3 feedback
3. **Growth-plate chondrocyte biology** — resting/proliferative/prehypertrophic/hypertrophic zones, p21^CIP1 cell-cycle arrest, SOX9/RUNX2, PTHrP-Ihh feedback, column disorganization
4. **CNP-NPR2 counter-regulatory axis (drug target)** — CNP-NPR-B-cGMP-PKGII-RAF1 inhibitory phosphorylation; NPR-C/neprilysin clearance; vosoritide binds NPR-B
5. **Craniofacial & skull base** — cranial-base synchondroses, premature fusion, foramen magnum stenosis, midface hypoplasia, jugular foramen/venous outflow, ICP, ventriculomegaly
6. **Spine** — short pedicles, spinal canal narrowing, thoracolumbar kyphosis (infancy), lumbar lordosis, progressive adult spinal stenosis, neurogenic claudication
7. **Neurologic complications** — cervicomedullary compression, hydrocephalus, syringomyelia, central sleep apnea, hyperreflexia, sudden-infant-death risk, motor delay
8. **ENT / respiratory** — Eustachian tube dysfunction, recurrent otitis media, conductive hearing loss, adenotonsillar hypertrophy, OSA/AHI, restrictive lung disease
9. **Orthopedic / limb** — rhizomelic shortening, trident hand, genu varum, fibular overgrowth, joint laxity, osteotomy, limb lengthening, adult height (~124-131 cm)
10. **Metabolic / systemic comorbidities** — reduced mobility, obesity risk, cardiometabolic risk, chronic pain, QoL, psychosocial impact, life expectancy
11. **Drug PK** — vosoritide (SC QD), TransCon CNP/navepegritide (SC QW sustained release), infigratinib (PO QD, CYP3A4), growth hormone, adherence/immunogenicity
12. **Drug PD** — NPR-B occupancy → cGMP → pERK inhibition → chondrocyte rescue → growth-plate rescue; direct FGFR3 kinase inhibition; off-target FGFR1 (hyperphosphatemia); CNP vasodilation/reflex tachycardia
13. **Clinical endpoints** — annualized growth velocity (AGV), height Z-score, upper/lower ratio, foramen magnum area, spinal canal Z, BMI-Z, motor milestones, sleep AHI, otitis frequency, PROMs
14. **Surgical / interventional** — foramen magnum decompression, spinal decompression/fusion, adenotonsillectomy, ear tubes, limb osteotomy/lengthening, AAP-guideline age-specific monitoring
15. **Safety / adverse events** — injection-site reactions, transient hypotension, vomiting, ear pain, corneal deposits/nail toxicity (FGFR-TKI class), growth-plate & long-term safety surveillance

---

## mrgsolve Model

### ODE Compartments (20)
**PK (6):** VOS_DEPOT, VOS_CP (vosoritide), TCNP_DEPOT, TCNP_CP (TransCon CNP), INFIG_GUT, INFIG_CP (infigratinib)

**Disease / PD / clinical (14):** PERK, CGMP_SIG, CHONDRO, HEIGHT_CM, HEIGHTZ,
FMAREA, SPCANALZ, AHI, OTITIS, BMIZ, MAP_BP, HR, PHOS (+ derived AGV_CALC)

### Treatment Scenarios (10)
1. **Untreated** — natural history reference (AGV ≈ 3.9 cm/yr)
2. **Vosoritide 15 µg/kg SC QD** — approved dose (Savarirayan 2020 Lancet, ΔAGV +1.57 cm/yr)
3. **Vosoritide 2.5 µg/kg SC QD** — phase 2 low dose
4. **Vosoritide 7.5 µg/kg SC QD** — phase 2 intermediate dose
5. **Vosoritide 30 µg/kg SC QD** — supratherapeutic, efficacy plateau vs 15 µg/kg
6. **TransCon CNP (navepegritide) SC QW** — sustained-release prodrug; pivotal ApproaCH phase 3, FDA-approved Feb-2026 as YUVIWEL (age ≥2 yr)
7. **Infigratinib PO QD** — FGFR1-3 TKI, direct kinase blockade (PROPEL 2/3 phase 2/3, positive; NDA planned Q3-2026, not yet approved)
8. **Growth hormone (off-label)** — historical, small/limited effect
9. **Vosoritide + FMD surgery** — combined pharmacologic + surgical decompression subgroup
10. **Vosoritide, 60% adherence** — real-world injection-adherence sensitivity

### Calibration Anchors
| Trial / source | Endpoint | Modeled target |
|---|---|---|
| Savarirayan 2020 Lancet (PMID 32891212, phase 3) | ΔAGV vosoritide 15 µg/kg vs placebo | +1.57 cm/yr |
| Savarirayan 2019 NEJM (PMID 31269546, phase 2) | Dose-response 2.5/7.5/15/30 µg/kg | plateau at 15 µg/kg |
| Horton 1978 J Pediatr (PMID 690757) | Untreated ACH growth velocity | ~3.9 cm/yr (age 5-14) |
| Ascendis ApproaCH phase 3 (topline Sep-2024) | ΔAGV TransCon CNP vs placebo | +1.49 cm/yr |
| BridgeBio PROPEL 3 (NEJM 2026) | Infigratinib off-target FGFR1 signal | mild/transient hyperphosphatemia ~4%, no discontinuations |
| Hunter 1998 J Med Genet (PMID 9733026) | Foramen magnum / cervicomedullary natural history | slow structural remodeling |

---

## Shiny App — 8 Tabs

1. **Patient & Overview** — covariate sidebar + mechanistic-map schematic
2. **Drug PK** — vosoritide/TransCon CNP and infigratinib concentration-time
3. **Pathway PD** — pERK activity vs. cGMP counter-signal; chondrocyte rescue index
4. **Growth endpoints** — AGV, height Z-score, cumulative height trajectories
5. **Structural & biomarker endpoints** — foramen magnum area, spinal canal Z, OSA-AHI, otitis, BMI-Z
6. **Scenario comparison** — all 10 regimens overlaid + endpoint table
7. **Safety** — hemodynamics (MAP/HR) and serum phosphate (FGFR1 off-target)
8. **References** — key trial citations

---

## Usage

```bash
# Render the mechanistic map
dot -Tsvg acdp_qsp_model.dot -o acdp_qsp_model.svg
dot -Tpng -Gdpi=150 acdp_qsp_model.dot -o acdp_qsp_model.png
```

```r
# Run a scenario in R
library(mrgsolve); library(dplyr); library(ggplot2)
source("acdp_mrgsolve_model.R")           # builds `acdp_mod` + `scenarios`
res <- run_scenario("2_Vosoritide_15ugkg_QD", scenarios[["2_Vosoritide_15ugkg_QD"]])
plot(res$time/24, res$HEIGHTZ, type = "l")

# Launch the dashboard
shiny::runApp("acdp_shiny_app.R")
```

---

## Disease-specific quick-reference (Korean / English)

| 한국어 | English |
|---|---|
| 연골무형성증 | Achondroplasia (ACH) |
| 섬유아세포성장인자수용체 3 | Fibroblast growth factor receptor 3 (FGFR3) |
| 성장판 | Growth plate |
| 근위지 단축 | Rhizomelic shortening |
| 대후두공 협착 | Foramen magnum stenosis |
| 척추관 협착 | Spinal stenosis |
| C형 나트륨이용성 펩타이드 | C-type natriuretic peptide (CNP) |
| 연간 성장 속도 | Annualized growth velocity (AGV) |
| 폐쇄성 수면무호흡 | Obstructive sleep apnea (OSA) |
| 재발성 중이염 | Recurrent otitis media |

---

*Built by Claude Code Routine on 2026-07-01 as part of the QSP Disease Model
Library. See root [README.md](../README.md) for the full model gallery.*
