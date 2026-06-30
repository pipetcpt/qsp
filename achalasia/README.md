# Achalasia (식도이완불능증) — QSP Disease Model

> **Primary idiopathic achalasia** is an esophageal motility disorder in which
> loss of inhibitory NO/VIP-ergic neurons of the myenteric (Auerbach's) plexus
> causes (i) **failure of lower esophageal sphincter (LES) relaxation** and
> (ii) **aperistalsis of the esophageal body**, leading to dysphagia for solids
> and liquids, regurgitation of undigested food, chest pain, and progressive
> weight loss.

The model captures: HLA-DQA1*0103/DQB1*0603 + putative HSV-1/VZV trigger →
T-cell mediated myenteric ganglionitis → **nNOS/VIP deficit** → cGMP/cAMP
shortage in LES smooth muscle → high resting LES tone (>35 mmHg) + failed
deglutitive relaxation (IRP > 15 mmHg on HRM, Chicago Classification v4.0)
+ esophageal aperistalsis/dilatation/sigmoid evolution → **Eckardt 0–12**,
**TBE**, HRQoL, esophageal cancer risk; treated with **isosorbide dinitrate /
nifedipine / sildenafil / tadalafil / intra-LES botulinum toxin A / pneumatic
dilation / Heller myotomy / POEM**.

---

## Deliverables

| File | Purpose |
|------|---------|
| `ach_qsp_model.dot` / `.svg` / `.png` | Graphviz mechanistic map — **14 clusters, 110+ nodes** |
| `ach_mrgsolve_model.R` | **23-ODE** mrgsolve model (10 PK + 12 disease/clinical + 2 AE) with 10 scenarios |
| `ach_shiny_app.R` | Interactive **8-tab** Shiny dashboard |
| `ach_references.md` | **82** PubMed citations grouped by section |

---

## Mechanistic Map — Cluster Index

1. **Etiology & triggers** — HLA-DQA1*0103/DQB1*0603 risk, HSV-1/VZV reactivation, Chagas T. cruzi (exclude), Allgrove syndrome, Down syndrome, familial form (<5%)
2. **Autoimmune ganglionitis** — anti-myenteric IgG, CD8⁺ T cytotoxicity, Th17/IL-23, Treg loss, mast cells, lymphocytic fibrosis, ICC loss
3. **Neurochemistry** — nNOS → NO → sGC → cGMP; VIP/PACAP → VPAC1 → cAMP; counter-balance of ACh/SP excitatory; GABA inhibitory
4. **LES smooth muscle** — Cav1.2 → Ca²⁺ → MLCK; MLCP regulation by cGMP/cAMP; RhoA/ROCK Ca²⁺-sensitization; K-ATP hyperpolarization
5. **Esophageal body** — swallow → propulsive peristalsis → deglutitive inhibition; aperistalsis types I/II/III (DCI, pan-esophageal pressurization, spasm)
6. **Chicago Classification v4.0** — HRM IRP > 15 mmHg, EndoFLIP distensibility, TBE column at 1/2/5 min, EGD birdbeak, provocative maneuvers
7. **Drug PK** — ISDN SL, nifedipine SL, sildenafil PO, tadalafil PO, intra-LES botulinum toxin A; CYP3A4 dependence
8. **Drug PD** — NO donation (ISDN), Cav1.2 blockade (nifedipine), PDE5 inhibition (sildenafil/tadalafil), SNAP-25 cleavage (botox); Bliss-independent combination
9. **Procedures** — pneumatic dilation (graded 30→35→40 mm), Heller myotomy + Dor/Toupet, POEM (long for Type III)
10. **Clinical endpoints** — dysphagia, regurgitation, chest pain, weight loss → **Eckardt 0–12**, TBE 5-min column, residual IRP, GERD-HRQL
11. **Long-term complications** — sigmoid megaesophagus, esophageal SCC (16-fold), aspiration pneumonia, malnutrition, esophageal candidiasis
12. **Patient covariates** — age, sex, BMI, CYP3A4 phenotype, cardiac history, disease duration, Chicago subtype, surgical candidacy
13. **Pivotal trials & guidelines** — Boeckxstaens 2011 NEJM, Werner 2019 NEJM, Ponds 2019 JAMA, Pasricha 1995, Triadafilopoulos 1991, Bortolotti 2000, ACG 2020, ISDE 2018
14. **Safety / AE** — nitrate headache/hypotension/tolerance, CCB edema, PDE5 NAION, botox heartburn, PD perforation 1–3%, POEM reflux 30–50%

---

## mrgsolve Model

### ODE Compartments (23)
**PK (10):** ISDN_GUT, ISDN_CEN, ISMN_CEN (active metabolite), NIF_GUT, NIF_CEN, SIL_GUT, SIL_CEN, TAD_GUT, TAD_CEN, BOTOX_LES (local depot)

**Disease physiology (11):** LES_PRESS, IRP_C, ESO_DIL, ESO_STAS, PERIST,
DYS_S, REG_S, CP_S, WT_S, TBE5, QOL

**Adverse events (2):** AE_HA (headache), AE_HYPO (postural hypotension)

### Treatment Scenarios (10)
1. **S01_NoTx** — natural history reference
2. **S02_ISDN** — Isosorbide dinitrate 10 mg SL TID
3. **S03_Nifedipine** — 20 mg SL before meals
4. **S04_Sildenafil** — 50 mg TID
5. **S05_Tadalafil** — 20 mg QD (long t½ off-label)
6. **S06_Botox** — Botulinum toxin A 100 U intra-LES q6 mo
7. **S07_PneumaticDilation** — 35 mm balloon (graded)
8. **S08_HellerMyotomy** — Laparoscopic + Dor fundoplication
9. **S09_POEM** — Peroral endoscopic myotomy (Type III: long myotomy)
10. **S10_Combo** — ISDN 10 mg + Sildenafil 25 mg SL TID

### Calibration Anchors
| Trial / source | Endpoint | Modeled target |
|---|---|---|
| Boeckxstaens 2011 NEJM | PD vs LHM 2 y success | 86% vs 90% |
| Werner 2019 NEJM | POEM vs LHM 2 y success | 83% vs 81% |
| Ponds 2019 JAMA | POEM vs PD 2 y success | 92% vs 54% |
| Pasricha 1995 NEJM | Botox response 1 mo | ~70%, falls to ~32% at 6 mo |
| Triadafilopoulos 1991 DDS | Nifedipine SL | ΔLES −30–40% |
| Bortolotti 2000 GE | Sildenafil 50 mg | ΔLES −35% × 2 h |
| Eckardt 1992 GE | Success ≤ 3 | Calibration threshold |
| Yadlapati 2021 NGM | Chicago v4.0 | IRP > 15 = abnormal |

---

## Shiny App — 8 Tabs

1. **Patient profile** — covariate-driven baseline trajectory plot + summary
2. **Drug PK** — log-scale Cp tracks for ISDN, NIF, SIL, TAD + Botox local depot
3. **LES / IRP** — manometric endpoint trajectories + subtype comparison
4. **Clinical endpoints** — Eckardt components + total, TBE column, QoL
5. **Scenario comparison** — multi-arm Eckardt + LES pressure
6. **Procedure outcomes** — PD vs LHM vs POEM 2-y / 5-y / GERD bar plot + table
7. **Safety** — headache + postural hypotension severity over time
8. **Calibration** — searchable table of trial anchors and how they're used

---

## Usage

```bash
# Render the mechanistic map
dot -Tsvg ach_qsp_model.dot -o ach_qsp_model.svg
dot -Tpng -Gdpi=110 ach_qsp_model.dot -o ach_qsp_model.png
```

```r
# Run a scenario in R
library(mrgsolve); library(dplyr); library(ggplot2)
source("ach_mrgsolve_model.R")           # builds `mod_ach` + `scenarios`
res <- run_scenario(scenarios[[9]],      # S09_POEM
                    subtype = 3,         # Type III spastic
                    end_days = 720)
plot(res$time, res$ECKARDT, type = "l")

# Launch the dashboard
shiny::runApp("ach_shiny_app.R")
```

---

## Disease-specific quick-reference (Korean / English)

| 한국어 | English |
|---|---|
| 식도이완불능증 | Achalasia |
| 하부식도괄약근 | Lower esophageal sphincter (LES) |
| 통합 이완 압력 | Integrated relaxation pressure (IRP) |
| 시간 바륨 식도조영 | Timed barium esophagram (TBE) |
| 풍선확장술 | Pneumatic dilation (PD) |
| 헬러근절개술 | Heller myotomy (LHM) |
| 경구 내시경 근절개술 | Per-oral endoscopic myotomy (POEM) |
| 미엔테릭 신경총 | Myenteric (Auerbach's) plexus |
| 카할 간질세포 | Interstitial cells of Cajal (ICC) |
| 거대식도증 | Megaesophagus |

---

*Built by Claude Code Routine on 2026-06-30 as part of the QSP Disease Model
Library. See repository root for full gallery.*
