# Cluster Headache (CH) — QSP Model

> Hypothalamic circadian generator → trigeminovascular CGRP/PACAP release → parasympathetic outflow (SPG) and partial Horner's → 15–180 min attacks of unilateral periorbital VAS 9–10/10 pain occurring in stereotyped circadian/circannual bouts. Episodic in 80–90 %, chronic in 10–20 %. Male:female ≈ 3:1.

## Files
| File | Description |
|------|-------------|
| [`ch_qsp_model.dot`](ch_qsp_model.dot) | Graphviz source — 12 clusters · 104 nodes |
| [`ch_qsp_model.svg`](ch_qsp_model.svg) · [`ch_qsp_model.png`](ch_qsp_model.png) | Rendered mechanistic map |
| [`ch_mrgsolve_model.R`](ch_mrgsolve_model.R) | 25-compartment ODE QSP, 7 drugs + O2 + GON block |
| [`ch_shiny_app.R`](ch_shiny_app.R) | 8-tab Shiny dashboard |
| [`ch_references.md`](ch_references.md) | 72 PubMed references grouped by topic |

## Mechanistic map — 12 clusters
1. **Genetic & chronobiologic susceptibility** — *HCRTR2, ADCYAP1, MTNR1A, CLOCK, BMAL1, ADH4*, family Hx ×5–18.
2. **Hypothalamic generator** — posterior hypothalamic gray (May 1998 *Lancet*), SCN / PVN / orexin / pineal melatonin, blunted cortisol & testosterone, in-bout vs remission gate.
3. **Trigeminovascular activation** — V1 trigeminal ganglion → dura · ICA · pial vessels; CGRP (jugular ↑ in attack, Goadsby 1994), PACAP-38, VIP, substance P, NO.
4. **Parasympathetic limb (trigeminal-autonomic reflex)** — SSN → GSPN → **SPG** → lacrimation, rhinorrhea, conjunctival injection, forehead sweating.
5. **Sympathetic disruption** — pericarotid plexus → SCG → ipsilateral ptosis + miosis (partial Horner's).
6. **Central pain processing** — TCC (NMDA/AMPA) → ipsilateral thalamic VPM/Po → S1·insula·ACC; restlessness / agitation distinguishes CH from migraine.
7. **Receptor pharmacology targets** — 5-HT 1B/1D/1F, L-type Ca²⁺, GSK-3β/IMPase (Li), OX2R, SSTR2/5, TRPV1, CGRP/CGRP-R mAbs.
8. **Drugs** — acute (O2, sumatriptan SC, zolmitriptan IN, lidocaine IN, octreotide), transitional (prednisone, GON block), preventive (verapamil, lithium, topiramate, melatonin, **galcanezumab CGRP-mAb**, erenumab off-label, civamide IN), devices (SPG-stim, nVNS, ONS, posterior hypothalamic DBS), research (psilocybin / LSA).
9. **Pharmacokinetics** — sumatriptan SC (CL 18 L/h, t½ ≈ 2 h), zolmitriptan IN, verapamil PR (CYP3A4, active norverapamil), lithium (renal CL ≈25 mL/min), topiramate (CL 1.2 L/h, t½ 21 h), galcanezumab (CL 0.008 L/h, t½ 27 d), prednisolone.
10. **Clinical endpoints** — attacks/week (primary), VAS pain, pain-free 15 min, ≥50 % responder rate, serum CGRP / VIP, AE bands (PR/QT, lithium tox, topiramate cognition, mAb injection-site).
11. **Triggers / lifestyle** — alcohol (bout-on only), histamine SC, nitroglycerin, nitrites/MSG, altitude/flights, odors, heat, daytime naps.
12. **Comorbidity & burden** — depression 55–65 %, suicide ideation 20–25 % (highest among pain disorders), CV risk on verapamil, lost work, HIT-6/CH-NDI.

## ODE model — 25 compartments
- **Drug PK (16)**: sumatriptan SC (1-cmt), zolmitriptan IN (2-cmt + peripheral), verapamil PR (2-cmt), lithium (1-cmt, renal-CL adjusted by CrCL), topiramate (1-cmt), galcanezumab SC (1-cmt with first-order SC absorption, linear FcRn-recycled), prednisolone.
- **Disease PD (9)**: hypothalamic drive (circadian + bout gate), CGRP tone, PACAP tone, pial effect site, attack hazard (smoothed), cumulative attacks, bout timer, O2 effect compartment, GON-block effect compartment (decays t½ ≈ 14 d).
- **Composite preventive effect** is a multiplicative escape: `1 − ∏(1 − Eᵢ)` over verapamil, lithium, topiramate, galcanezumab, prednisolone, GON.
- **Acute abortive** is a similar product over sumatriptan, zolmitriptan and O2 effect.
- **Calibrations**: Cohen 2009 *JAMA* (O2 78 % pain-free@15 min), Ekbom 1991 *NEJM* (SC sumatriptan 74 %), Leone 2000 *Neurology* (verapamil 240 mg/d), Steiner 1997 *Cephalalgia* (lithium), **Goadsby 2019 *NEJM* / Dodick 2020 *Cephalalgia* (galcanezumab ECH/CCH)**, Obermann 2021 *Lancet Neurol* (prednisone bridge), Leroux 2011 *Lancet Neurol* (GON block).

## Six treatment scenarios (12-week simulation horizon)
| Scenario | Description |
|----------|-------------|
| **S0** | No treatment (natural bout) |
| **S1** | O2 + sumatriptan 6 mg SC for an indexed attack |
| **S2** | Verapamil 240 mg PR BID |
| **S3** | Verapamil + lithium 300 mg BID (chronic CH) |
| **S4** | Galcanezumab 300 mg SC q4w |
| **S5** | Prednisone 60 → 40 → 20 mg taper + verapamil |
| **S6** | GON block (one-shot 0.65 effect) + verapamil |

## Shiny dashboard — 8 tabs
1. Patient & disease profile (chronic vs episodic, hypothalamic drive)
2. Drug PK (acute / preventive / mAb)
3. Pathway PD (hypothalamic drive · CGRP/PACAP · pial)
4. Clinical endpoints (attacks/week, hazard, trial-anchored table)
5. Scenario comparison (S0–S6)
6. Biomarkers (CGRP, PACAP, pial)
7. Safety (verapamil PR/QT, lithium therapeutic band, galcanezumab exposure)
8. References (renders `ch_references.md`)

## Run
```bash
# Render map
dot -Tsvg ch_qsp_model.dot -o ch_qsp_model.svg
dot -Tpng -Gdpi=150 ch_qsp_model.dot -o ch_qsp_model.png

# Simulate (R)
Rscript -e 'source("ch_mrgsolve_model.R"); print(simulate_all())'

# Dashboard
Rscript -e 'shiny::runApp("ch_shiny_app.R")'
```

*Built by Claude Code Routine — 2026-06-30.*
