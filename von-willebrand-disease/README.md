# Von Willebrand Disease (폰 빌레브란트병) — QSP Disease Model

> **Von Willebrand Disease (VWD)** is the most common inherited bleeding
> disorder, caused by **quantitative** (Type 1 partial deficiency, Type 3
> total deficiency) or **qualitative** (Type 2A/2B/2M/2N) defects of the
> **VWF gene (12p13.31)**. VWF is a multimeric adhesive glycoprotein that
> mediates **platelet-collagen/GPIbα adhesion** (primary hemostasis) and
> **stabilizes circulating FVIII** (secondary hemostasis); its size-dependent
> function is regulated by the shear-activated protease **ADAMTS13**. Loss of
> function produces mucocutaneous bleeding, menorrhagia, and — in severe
> quantitative or multimer-depleted qualitative subtypes — GI angiodysplasia
> bleeding and postpartum hemorrhage.

The model captures: **VWF gene defect** → altered biosynthesis/multimerization
(Weibel-Palade bodies) or altered **GPIbα (A1)/collagen (A3)/FVIII (D'D3)**
binding → impaired **platelet adhesion** and/or **FVIII stabilization** →
mucocutaneous bleeding, menorrhagia, GI bleeding (Type 2A/3/acquired),
postpartum hemorrhage. **ADAMTS13** cleaves shear-unfolded high-molecular-weight
multimers (HMWM); Type 2B/platelet-type GPIbα gain-of-function causes
spontaneous platelet aggregation/clearance (thrombocytopenia). **Desmopressin
(DDAVP)** acutely releases endogenous VWF/FVIII/VWFpp from Weibel-Palade
bodies via a V2-like endothelial receptor (tachyphylaxis with repeated
dosing; contraindicated in Type 2B/platelet-type). **Recombinant VWF**
(vonicog alfa) and **plasma-derived VWF/FVIII concentrate** (Humate-P/Wilate-
like) restore multimers directly; **tranexamic acid** and **hormonal
therapy** address menorrhagia/surgical bleeding.

---

## Deliverables

| File | Purpose |
|------|---------|
| `vwd_qsp_model.dot` / `.svg` / `.png` | Graphviz mechanistic map — **12 clusters, 124 nodes** |
| `vwd_mrgsolve_model.R` | **21-ODE** mrgsolve model (6 PK + 15 disease/clinical) with 10 scenarios |
| `vwd_shiny_app.R` | Interactive **8-tab** Shiny dashboard |
| `vwd_references.md` | **61** PubMed citations grouped by section (all verified) |

---

## Mechanistic Map — Cluster Index

1. **Genetics & classification** — VWF 12p13.31, autosomal dominant (Type 1/2A/2B/2M) vs recessive (Type 3/2N), platelet-type pseudo-VWD, acquired VWS, ABO modifier
2. **VWF biosynthesis & multimerization** — pre-pro-VWF, ER dimerization, Golgi multimerization, furin cleavage (VWFpp), ultra-large VWF (ULVWF), Weibel-Palade body storage, constitutive/regulated secretion
3. **ADAMTS13 & multimer-size regulation** — shear-induced A2-domain unfolding, Tyr1605-Met1606 cleavage, HMWM loss, acquired anti-ADAMTS13 antibody, Heyde syndrome (aortic stenosis shear)
4. **Platelet adhesion / primary hemostasis** — subendothelial collagen, VWF A3-collagen and A1-GPIbα binding, GPIbα gain-of-function (Type 2B/platelet-type), GPIIb/IIIa activation, aggregation, Type 2B RES clearance
5. **FVIII stabilization / secondary hemostasis** — VWF D'D3-FVIII binding (Type 2N defect), FVIII protection from APC/protein C, tenase complex, thrombin generation, clot strength
6. **Clinical bleeding phenotype** — epistaxis, gingival/oral bleeding, easy bruising, menorrhagia, GI angiodysplasia bleeding, post-surgical/dental bleeding, postpartum hemorrhage, ISTH-BAT score
7. **Diagnostic biomarkers** — VWF:Ag, VWF:RCo/GPIbM, VWF:CB, FVIII:C, RCo/Ag ratio, multimer gel electrophoresis, RIPA, VWFpp:Ag ratio, PFA closure time, blood-type-O modulation
8. **Desmopressin (DDAVP) PK/PD** — V2 renal (AQP2/hyponatremia) vs V2-like endothelial WPB release, tachyphylaxis, Type 2B contraindication, trial-of-response testing
9. **VWF/FVIII replacement PK/PD** — plasma-derived concentrate (co-formulated FVIII, immediate rise) vs recombinant VWF (VWF-only, delayed endogenous FVIII stabilization), PK-guided perioperative dosing
10. **Adjunctive therapies** — tranexamic acid/aminocaproic acid (antifibrinolytic), hormonal therapy (COC/LNG-IUS), iron repletion, fibrin sealant, antiplatelet avoidance
11. **Special populations / modulators** — pregnancy 3rd-trimester VWF/FVIII rise, postpartum fall (delayed PPH risk), acute-phase/exercise/age modulation, Type 2B pregnancy risk
12. **Safety / adverse events** — DDAVP hyponatremia/seizure, thrombotic risk from factor overcorrection, allergic reactions, rare alloantibody (Type 3), volume overload, pathogen risk (historical)

---

## mrgsolve Model

### ODE Compartments (21)
**PK (6):** DDAVP_DEPOT, DDAVP_CP (desmopressin), RVWF_CP (recombinant VWF),
PDVWF_CP (plasma-derived VWF/FVIII concentrate), TXA_DEPOT, TXA_CP (tranexamic acid)

**Disease / PD / clinical (15):** WPB_STORE, VWFPP, VWF_AG, VWF_RCO, HMWM,
ADAMTS13_ACT, FVIII_C, PLT_COUNT, BLEED_SCORE, MENS_LOSS, GI_LOSS, HB,
NA_SERUM, THROMB_RISK (+ derived VWF_AG_TOTAL/VWF_RCO_TOTAL/FVIII_C_TOTAL/HMWM_EFFECTIVE)

### Treatment Scenarios (10)
1. **Type 1, mild, untreated** — natural-history reference
2. **Type 1, DDAVP IV single dose (0.3 mcg/kg)** — trial-of-response
3. **Type 1, DDAVP intranasal, repeated (q12h x3d)** — tachyphylaxis demonstration
4. **Type 2B, DDAVP** — illustrates worsened thrombocytopenia (contraindicated)
5. **Type 3, severe, recombinant VWF (vonicog alfa-like)** — delayed endogenous FVIII rise
6. **Type 3, severe, plasma-derived VWF/FVIII concentrate** — immediate co-formulated FVIII correction
7. **Menorrhagia — tranexamic acid + hormonal therapy** — heavy menstrual bleeding management
8. **Acquired VWS — transient DDAVP response** — shear/lymphoproliferative-associated
9. **Pregnancy, Type 1, peripartum** — 3rd-trimester physiologic boost, postpartum fall
10. **Major surgery — PK-guided plasma-derived concentrate (q12h)** — perioperative trough targeting

### Calibration Anchors
| Trial / source | Endpoint | Modeled target |
|---|---|---|
| James/Connell 2021 Blood Adv, PMID 33570651 (ASH/ISTH/NHF/WFH diagnosis guideline) | VWF:RCo/Ag diagnostic thresholds | Type 1 <30-50 IU/dL; Type 3 <3 IU/dL |
| Connell/Flood 2021 Blood Adv, PMID 33570647 (ASH/ISTH/NHF/WFH management guideline) | DDAVP trial-of-response, perioperative PK-guided dosing | scenario design (2, 3, 10) |
| Mannucci 1997 Blood, PMID 9326215 (DDAVP review) | Tachyphylaxis with repeated dosing | K_WPBDEPLETE vs slow K_WPBREGEN |
| Federici 2004 Blood, PMID 14630825 (multicenter biologic-response study) | DDAVP response magnitude/duration by VWD type | EC50_DDAVP / EMAX_WPBREL |
| Mannucci 2013 Blood (PMID 23777763); Gill 2015 Blood (PMID 26239086) (rVWF/vonicog alfa) | VWF:RCo t1/2 ~21h; delayed endogenous FVIII climb | KE_RVWF; FVIII_STAB_DRIVE routing |
| Dobrkovska 1998 Haemophilia, PMID 10028316 (Humate-P PK) | VWF:RCo t1/2 ~12-20h; FVIII:C t1/2 ~8-12h | KE_PDVWF; VWF_FVIII_RATIO=2.4 |
| Leebeek & Eikenboom 2016 NEJM, PMID 27959741 (comprehensive review) | Overall disease/PK-PD framework | model scope & structure |

---

## Shiny App — 8 Tabs

1. **Patient & Overview** — genotype/regimen sidebar + mechanistic-map schematic
2. **Drug PK** — DDAVP/TXA and rVWF/PD-concentrate concentration-time
3. **Hemostatic biomarkers** — VWF:Ag/RCo/HMWM and FVIII:C/platelets/ADAMTS13
4. **Clinical endpoints** — bleeding score, menstrual & GI blood loss trajectories
5. **Scenario comparison** — all 7 regimens overlaid + endpoint table
6. **Biomarkers & diagnostics** — VWFpp acute marker, VWF:RCo/Ag ratio discriminator
7. **Safety** — serum sodium (DDAVP) and thrombotic-risk index (overcorrection)
8. **References** — key trial/guideline citations

---

## Usage

```bash
# Render the mechanistic map
dot -Tsvg vwd_qsp_model.dot -o vwd_qsp_model.svg
dot -Tpng -Gdpi=150 vwd_qsp_model.dot -o vwd_qsp_model.png
```

```r
# Run a scenario in R
library(mrgsolve); library(dplyr); library(ggplot2)
source("vwd_mrgsolve_model.R")           # builds `vwd_mod` + `scenarios`
res <- run_scenario("5_Type3_Severe_rVWF_Vonvendi", scenarios[["5_Type3_Severe_rVWF_Vonvendi"]], p_type3, end = 168)
plot(res$time, res$VWF_RCO_TOTAL, type = "l")

# Launch the dashboard
shiny::runApp("vwd_shiny_app.R")
```

---

## Disease-specific quick-reference (Korean / English)

| 한국어 | English |
|---|---|
| 폰 빌레브란트병 | Von Willebrand Disease (VWD) |
| 폰 빌레브란트 인자 | Von Willebrand Factor (VWF) |
| 혈소판 부착 | Platelet adhesion |
| 당단백질 Ib알파 | Glycoprotein Ibα (GPIbα) |
| 8번 인자 안정화 | Factor VIII stabilization |
| 고분자량 멀티머 | High-molecular-weight multimer (HMWM) |
| 데스모프레신 | Desmopressin (DDAVP) |
| 월경과다 | Menorrhagia / heavy menstrual bleeding |
| 위장관 출혈 | Gastrointestinal bleeding |
| 산후 출혈 | Postpartum hemorrhage |

---

*Built by Claude Code Routine on 2026-07-01 as part of the QSP Disease Model
Library. See root [README.md](../README.md) for the full model gallery.*
