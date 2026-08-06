# 대퇴골두 무혈성 괴사 (Osteonecrosis of the Femoral Head) — QSP 모델

**The necrotic bone never fails. The repair does.**

Dead bone in the femoral head is fully mineralised and as stiff as it ever was.
It has no living cells, so it never resorbs, never remodels, and never gets
weaker. A femoral head full of it can carry body weight indefinitely — and
sometimes does, silently, for years.

What breaks the hip is creeping substitution: osteoclasts excavate the dead
trabeculae *before* osteoblasts refill them, so the interface between dead and
living bone is mechanically weaker while it heals than it was while it was
merely dead. And because the front keeps moving, the active zone is
continuously reset to mid-cycle, so the weakness lasts for as long as the front
is still crossing.

In this model the interface loses **two thirds of its stiffness** at month 4–6
and has it all back by year 3, while the stress across it does not move at all
(1.101 → 1.105 MPa). Collapse is a division whose numerator is constant.

---

## Three clocks

| | clock | timescale | what it does |
|---|---|---|---|
| **1** | **perfusion** | hours–days | Marrow adipogenesis and microthrombosis raise intraosseous pressure until the Starling resistor closes and the head infarcts. **This clock has stopped before the patient has a symptom.** Everything aimed at it — statin, anticoagulant, prostacyclin, and the venting half of a core decompression — can only work as prophylaxis. |
| **2** | **repair** | months | Creeping substitution. Resorption precedes formation, so healing *costs* strength before it returns any. Sets the *rate* of clock 3. |
| **3** | **fatigue** | months–years | Microdamage. Living bone erases it by targeted remodelling, which osteocytes signal. Necrotic bone has none, so its damage is permanent. **This clock alone sets the endpoint.** |

---

## The mechanical core: one cone, two areas

The lesion is a cone with its apex at the centre of the head, half-angle α, axis
tilted φ from the hip joint resultant force. Two areas follow from that and
nothing else:

```
loaded cap        A_cap = 2 π R² (1 − cos α)      grows like α²
conical interface A_int =   π R²   sin α          grows like α
```

Every newton that lands on the necrotic cap has to **leave** through the conical
interface, because the apex is a point. So

```
σ_int = L_eff / A_int          an equilibrium traction
S_int = S₀ · β²                β = 1 − CAV − 0.65·NB·(1 − MINZ)
```

and `A_cap/A_int = 2 tan(α/2)` rises without bound:

| α | Kerboull CNA | A_cap (mm²) | A_int (mm²) | A_cap/A_int |
|---:|---:|---:|---:|---:|
| 15° | 60° | 108 | 412 | 0.263 |
| 35° | 140° | 575 | 912 | 0.631 |
| 55° | 220° | 1356 | 1303 | 1.041 |
| 75° | 300° | 2358 | 1536 | 1.535 |

**That quotient is the whole of ONFH staging.** Kerboull's combined necrotic
angle measures α (CNA = 4α for a circular cone). The JIC types measure *where*
the cap sits relative to the acetabular contact patch, which is what sets
`L_eff`. Steinberg's percent-volume measures a mixture of the two — which is
why it predicts collapse worse than either.

The JIC class boundaries are not model parameters. The weight-bearing surface
spans ±θ_c = ±50°, so its medial-third and medial-two-thirds lines sit at
−16.7° and +16.7° by arithmetic. A/B/C1/C2 fall out of where the lateral lesion
boundary θ_L = φ + α lands.

---

## Deliverables

| file | what it is |
|---|---|
| [`onfh_qsp_model.dot`](onfh_qsp_model.dot) · [SVG](onfh_qsp_model.svg) · [PNG](onfh_qsp_model.png) | mechanistic map: 159 nodes, 18 clusters, 219 edges |
| [`onfh_mrgsolve_model.R`](onfh_mrgsolve_model.R) | 49-ODE mrgsolve model, 20 annotated scenarios |
| [`onfh_shiny_app.R`](onfh_shiny_app.R) | 12-tab dashboard, one tab per subsystem |
| [`onfh_references.md`](onfh_references.md) | 152 references, every PMID retrieved from PubMed by [`mkrefs.py`](mkrefs.py) rather than recalled |
| [`onfh_python_reference.py`](onfh_python_reference.py) | the independent numpy/scipy implementation — the source of every number below |
| [`onfh_reference_output.txt`](onfh_reference_output.txt) | its verbatim output |

---

## What is fitted

**One number is fitted to an outcome.** `k_dmg = 0.90`, the microdamage rate
constant, anchored on a single quantity: the five-year collapse rate of
untreated JIC C1 hips. Two more (`h0_tha`, `k_pain`) scale the arthroplasty
hazard and the pain axis and touch nothing mechanical.

Two structural constants — the fatigue exponent `m_fat = 4` and the
contact-to-trabecular spreading factor `f_spread = 0.30` — were chosen so that
(a) normal bone does not fatigue-fail in a heavy, active person and (b) latency
to collapse spans months rather than weeks. Those are qualitative requirements,
not fits to any collapse rate. `m_fat` is a *structural* exponent: trabecular
bone **specimens** give 12–16, but a structure sheds load from failing elements
to their neighbours, and σ_int here is a lumped equilibrium traction rather than
a local tissue stress. It is the least certain constant in the model and it
tops the sensitivity table.

Everything else is geometry, published material properties, or physiology.
Peak contact pressure, for example, is *derived*: 2.80 × 70 kg over a 50°
contact patch gives **p₀ = 2.468 MPa**, against in-vivo instrumented-endo-
prosthesis gait values of 2–5 MPa.

---

## Results

### 1. The vulnerable window

Untreated JIC C1 (α = 45°, φ = −12°). This hip collapses at month 8.4; rows
after that are omitted because the post-collapse incongruity factor inflates
every stress by up to 5.4×.

| month | XF (mm) | CAV | NB | MINZ | β | E_int (MPa) | S_int (MPa) | σ_int (MPa) | σ/S |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.50 | 0.000 | 0.000 | 0.000 | 1.000 | 620.0 | 9.500 | 1.101 | 0.116 |
| 3 | 2.64 | 0.258 | 0.202 | 0.541 | 0.710 | 263.4 | 4.790 | 1.101 | 0.230 |
| 6 | 6.18 | 0.252 | 0.211 | 0.521 | 0.668 | 226.3 | 4.242 | 1.105 | 0.260 |

Minimum at **month 4.4**: β = 0.661, modulus 220 MPa — **35.5% of normal**.
By month 60 it is back to β = 0.994.

The stress column is the point. σ_int moves by 0.4% across the whole window.
Everything that happens to σ/S happens in the denominator.

### 2. Which structure fails first

Pre-collapse (month 6), stress/strength at five probes:

| probe | σ/S | is its damage cleared? |
|---|---:|---|
| **reparative interface** | **0.260** | yes — but the tissue is weak |
| subchondral plate | 0.198 | partly |
| adjacent living bone (control) | 0.078 | yes |
| acetabular rim | 0.068 | only if the rim sits on living bone |
| untouched necrotic core | 0.053 | **no** — no osteocytes to signal it |

The necrotic core is the one place where damage is permanent, and it is the
*least* stressed place in the head. That is the model's answer to why dead bone
does not simply crumble: it is strong, and it is protected by its own size —
the cap area it presents to the joint grows faster than the load it collects.

### 3. Collapse rates by JIC class — only C1 was fitted

360 virtual hips, geometry and habitus sampled; **k_dmg was anchored on the C1
row alone.**

| JIC | n | 2 y | 3 y | 5 y | median month among collapsers | published (JIC / Mont 2010) |
|---|---:|---:|---:|---:|---:|---|
| A | 19 | 0.000 | 0.000 | 0.000 | — | rarely collapses |
| B | 120 | 0.025 | 0.025 | 0.025 | 7.7 | ~15–25% |
| **C1** | 165 | 0.448 | 0.461 | **0.479** | 5.8 | ~40–60% ← the anchor |
| C2 | 56 | 0.750 | 0.768 | **0.768** | 3.2 | ~70–85% |

**A and C2 are predictions and both land inside the reported ranges. B does
not** — see the misses below.

### 4. The hazard is not monotone, and nobody told it to be

Because the interface loses strength only while the front is crossing and
regains it afterwards, the collapse hazard must rise, peak and then vanish.
Conditional probability of collapsing in the *next* 12 months:

| still intact at | 0 mo | 12 mo | 24 mo | 36 mo | 48 mo |
|---|---:|---:|---:|---:|---:|
| P(collapse within 12 months) | 0.294 | 0.051 | 0.012 | 0.008 | 0.004 |

85% of all collapses occur inside 12 months, 95% inside 24. Clinically, a head
still spherical at four years is very unlikely to fail. Nothing was fitted to
that; it follows from the fact that creeping substitution *finishes*.

### 5. Same volume, opposite outcome

Eight lesions with the **identical** cone half-angle — identical Steinberg
percent-volume — differing only in where the cone points:

| φ | θ_L | JIC | fraction of joint load on the lesion | collapse |
|---:|---:|---|---:|---|
| −55° | −10° | B | 0.289 | never |
| −35° | +10° | B | 0.541 | never |
| −25° | +20° | C1 | 0.668 | never |
| −15° | +30° | C1 | 0.785 | 11.9 mo |
| −5° | +40° | C1 | 0.873 | 5.6 mo |
| +15° | +60° | C2 | 0.785 | 5.3 mo |

This is why a volume-based stage predicts collapse worse than a location-based
one, stated as arithmetic rather than as clinical lore.

### 6. The Kerboull threshold, derived rather than assumed

Sweeping the combined necrotic angle with the lesion axis on the load axis (the
worst case, θ_L = α), the deterministic collapse threshold falls **between CNA
125° and 140°**. In the virtual population:

| CNA band | n | 5-year collapse |
|---|---:|---:|
| < 150° | 190 | 0.126 |
| 150–190° | 84 | 0.512 |
| 190–240° | 62 | 0.613 |
| > 240° | 24 | 0.833 |

Kerboull's published risk bands are <190 low, 190–240 intermediate, >240 high.
The model contains no number resembling 190; it contains a cone, a cosine
pressure field and a power-law S-N curve, and it produces a threshold and a
monotone gradient across exactly those bands — but it puts the transition
**lower than Kerboull does**, around 150° rather than 190°. That gap is real
and is reported, not tuned away.

### 7. Bisphosphonate: the ceiling is structural

Alendronate acts on **one limb** of the interface. It suppresses the excavation
that makes the interface porous; it cannot make bone. And because osteoblast
recruitment follows resorbed surface, suppressing resorption also removes part
of the stimulus for the formation it is protecting — visible below as new bone
that starts *ahead* of control at month 6 and finishes *behind* it at month 36.

| month | CAV ctrl | CAV ALN | NB ctrl | NB ALN | β ctrl | β ALN |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 0.258 | 0.190 | 0.202 | 0.149 | 0.710 | 0.903 |
| 6 | 0.252 | 0.164 | 0.211 | 0.250 | 0.668 | 0.835 |
| 24 | 0.101 | 0.116 | 0.361 | 0.315 | 0.864 | 0.843 |
| 36 | 0.043 | 0.141 | 0.419 | 0.315 | 0.944 | 0.810 |

Minimum β rises 0.661 → 0.793 (+20.0%); peak σ/S falls 0.508 → 0.184 (−63.8%).
On the matched C1 hip that **delays** collapse from month 8.4 to month 27.2 —
it does not prevent it. At trial scale, in paired virtual patients enrolled as
both trials enrolled (large lesions, CNA ≥ 160°):

| | collapse at 24 months |
|---|---:|
| model, control | 0.649 |
| model, alendronate | 0.541 → **RR 0.83** |
| Chen 2012 (RCT, double-blind, ARCO IIC/IIIC) | 0.500 vs 0.515 → **RR 0.97** |
| Lai 2005 (Steinberg IIC/IIIC) | 0.069 vs 0.760 → **RR 0.09** |

The two randomised trials disagree by a factor of ten. The model **cannot
produce Lai's effect size at any degree of osteoclast blockade** — pushing
inhibition from 0.85 to 0.999 moves minimum β from 0.793 to 0.777, because
blockade cannot lift the interface above the strength it would have had with no
resorption at all, and the front still has to cross. The model lands close to
Chen. That is a falsifiable structural claim about the mechanism's ceiling, and
it is the model's actual position on the controversy rather than an average of
the two arms.

### 8. Core decompression changes sign

It does three things at once: vents the compartment (good — but clock 1 has
stopped), opens a channel for repair (good), and removes load-bearing bone
(bad). Only the third scales with the lesion.

| CNA | no CD | with CD | change | verdict |
|---:|---|---|---:|---|
| 125° | never | never | — | both fine |
| 135° | 21.7 mo | never | **rescued** | HELP |
| 145° | 10.7 mo | 14.8 mo | +4.0 mo | help |
| 170° | 5.7 mo | 6.1 mo | +0.4 mo | help |
| 200° | 4.0 mo | 4.0 mo | −0.0 mo | **harm** |
| 250° | 4.2 mo | 4.2 mo | +0.1 mo | irrelevant |

The benefit decays monotonically and the sign flips near CNA 200°, which is why
core decompression reads as effective in small pre-collapse lesions and as
futile in large ones without either result being wrong.

### 9. Peak steroid dose beats cumulative steroid dose

Cumulative prednisolone-equivalent held **constant at 4650 mg**; only the
schedule changes.

| daily mg | days | peak marrow exposure | P_io max | Q min | necrotic fraction | CNA | outcome |
|---:|---:|---:|---:|---:|---:|---:|---|
| 100 | 46 | 0.390 | 80.8 | 0.000 | 0.420 | 323° | collapse, month 4.6 |
| 60 | 78 | 0.235 | 77.5 | 0.000 | 0.420 | 323° | collapse, month 4.6 |
| 40 | 116 | 0.157 | 67.9 | 0.000 | 0.420 | 323° | collapse, month 4.6 |
| 25 | 186 | 0.098 | 50.8 | 0.069 | 0.391 | 309° | collapse, month 4.6 |
| 15 | 310 | 0.059 | 34.6 | 0.336 | 0.043 | 96° | never |
| 10 | 465 | 0.039 | 26.8 | 0.586 | 0.009 | 43° | never |
| 6 | 775 | 0.024 | 18.4 | 0.868 | 0.003 | — | **no lesion at all** |

Identical total dose; the schedule alone spans "hip destroyed in five months"
and "nothing happened". The model was given one Hill coefficient of 2.5 on the
adipogenic drive with its EC50 at the exposure a 40 mg/d schedule produces; the
dominance of peak over cumulative is then arithmetic. The transition sits
between 15 and 25 mg/d, i.e. somewhat *lower* than the 40 mg/d the epidemiology
emphasises.

### 10. Prophylaxis, and the only thing that works

| arm | P_io max | Q min | necrotic fraction | CNA | outcome |
|---|---:|---:|---:|---:|---|
| 60 mg/d × 30 d then taper | 73.5 | 0.000 | 0.420 | 323° | collapse, 4.6 mo |
| + rosuvastatin from day 0 | 61.7 | 0.000 | 0.420 | 323° | collapse, 4.6 mo |
| + enoxaparin × 12 weeks | 71.0 | 0.000 | 0.420 | 323° | collapse, 4.9 mo |
| + iloprost × 5 d at day 14 | 73.5 | 0.000 | 0.420 | 323° | collapse, 4.6 mo |
| + statin started 6 months later | 73.5 | 0.000 | 0.420 | 323° | collapse, 4.6 mo |
| **same total dose at 12.4 mg/d × 375 d** | 30.5 | 0.455 | 0.019 | 63° | **never** |

Statin measurably lowers intraosseous pressure (73.5 → 61.7 mmHg) and still
does not save the head, because the compartment closes anyway. In this model
the only intervention that prevents the infarct is **not giving the peak dose**.
That is a strong claim and it is the direct consequence of clock 1 being over
before anyone knows it started.

### 11. Ten arms on one matched C1 hip

| arm | collapse | P(THA) at 5 y | minimum β |
|---|---|---:|---:|
| untreated | 8.4 mo | 0.997 | 0.661 |
| iloprost | 8.4 mo | 0.997 | 0.661 |
| core decompression (single track) | 10.3 mo | 0.996 | 0.655 |
| multiple drilling + BMAC | 14.8 mo | 0.994 | 0.660 |
| denosumab | 16.7 mo | 0.992 | 0.682 |
| alendronate | 27.2 mo | 0.974 | 0.793 |
| alendronate → teriparatide | 34.0 mo | 0.943 | 0.793 |
| teriparatide | 37.2 mo | 0.922 | 0.736 |
| protected weight bearing alone | never | 0.156 | 0.661 |
| structural graft + BMAC | never | 0.000 | 0.660 |

Note the two columns that disagree. Alendronate raises minimum β the most, but
teriparatide delays collapse the longest — because the anabolic agent shortens
the window rather than propping up its floor. And note the bottom two rows,
which the model gets wrong; see below.

### 12. Sensitivity: the model is a mechanics model

±20%, output is months to collapse (censored at 60). Five-year depression is
useless as a sensitivity output because it saturates for every hip that
collapses at all — that was itself a defect the first run exposed.

| parameter | −20% | base | +20% | normalised sensitivity |
|---|---:|---:|---:|---:|
| m_fat | 4.67 | 8.41 | 60.00 | **+16.4** |
| theta_c | 4.69 | 8.41 | 60.00 | +16.4 |
| alpha_deg | 60.00 | 8.41 | 5.13 | −16.3 |
| f_hip / BW | 60.00 | 8.41 | 5.50 | −16.2 |
| S_tr0 | 5.14 | 8.41 | 49.62 | +13.2 |
| RESP_max | 22.09 | 8.41 | 5.41 | −5.0 |
| **k_dmg** (the fitted one) | 10.24 | 8.41 | 7.38 | **−0.85** |

Geometry, loading and material strength dominate; the repair-biology constants
are an order of magnitude weaker; and the **one number fitted to an outcome is
nearly the least influential parameter in the model.** `E_tr0` and `f_spread`
show exactly zero sensitivity because they act only on probes that never fail.

---

## Known misses, stated rather than fitted away

1. **Type B collapses far too rarely** — 2.5% against a reported 15–25%. Real
   JIC typing has measurement error and real lesions are not cones; both would
   move mass from B into C1 in the published series in a way this model does not
   reproduce. This is the largest quantitative failure.
2. **The Kerboull threshold lands near 150°, not 190°.** The gradient across
   his bands is right; the intercept is not.
3. **Protected weight bearing alone prevents collapse.** Clinically it does not.
   Cycles enter the damage law linearly and the model has no compliance, no
   muscle loading, and no notion that a patient cannot stay non-weight-bearing
   for the 6–12 months the window lasts. The model is describing a treatment
   nobody can actually take.
4. **A structural graft abolishes collapse entirely.** The parallel load path is
   modelled as a fixed 35% offload with no failure mode of its own.
5. **The two alendronate trials cannot both be reproduced.** The model sits near
   Chen and says Lai's effect size is above the mechanistic ceiling.
6. **Damage is tracked at five fixed probes, not as a field.** A real crescent
   propagates; here it appears.
7. **Post-collapse numbers are not comparable across probes.** Once the head
   depresses, the incongruity factor multiplies every stress by up to 5.4, so the
   "damage at 5 years" column in the reference output compares a collapsed hip
   with an intact one. The pre-collapse row is the fair comparison and is what
   section 2 above uses. Rows 2–4 of the scenario table in `[1]` similarly show
   depression and hazard for the incident-lesion arms computed on placeholder
   geometry; those arms are resolved properly in section `[10]`.
8. **No R toolchain was available**, so `onfh_mrgsolve_model.R` mirrors the
   executed Python reference equation for equation but has not itself been run.

---

## Verification: what the second implementation caught

Every equation was written twice — once in mrgsolve form, once independently in
Python/scipy — and the executed version was checked against physiological
requirements a correct model must satisfy (a healthy hip must be a fixed point;
normal bone must not fatigue-fail; a dose schedule must be able to change an
outcome). That found and fixed **ten defects**:

1. **Every hip collapsed, including type A.** The crescent term used a logistic
   whose value at zero damage was 3.4 × 10⁻⁴ rather than zero, so a floor of
   crescent propagation drove every head to the 8 mm depression cap.
2. **The disease disappeared.** The reparative-interface trough was cancelled by
   its own new-bone term: modelled as one synchronised remodelling cohort, the
   interface healed in weeks. Creeping substitution is a *moving* front that
   continuously resets the active zone to un-excavated; adding that renewal term
   is what makes the vulnerable window last as long as the front is crossing.
3. **Normal bone fatigue-failed in heavy people.** The contact probes were fed
   raw joint contact pressure, which put healthy trabecular bone *at* the
   endurance limit; a 100 kg patient's intact bone accumulated damage of 1273.
   Fixed by separating equilibrium traction (σ_int, unscaled) from local contact
   stress (scaled by the load-spreading factor).
4. **A normal femoral head slowly infarcted itself.** Perfusion was normalised
   to venous pressure, so a healthy hip sat at Q = 0.925, which fed the oedema
   term, which raised intraosseous pressure, which lowered Q further. The loop
   gain exceeded one. Consequence: *every* steroid schedule produced the
   identical maximal lesion, silently destroying the peak-versus-cumulative
   result.
5. **Healthy hips lost their cartilage.** A baseline cartilage-loss term
   independent of collapse cost an intact hip 45% of its cartilage over five
   years and scored it Harris Hip 53.9.
6. **Asymptomatic hips were replaced.** The arthroplasty hazard had an additive
   constant, giving a spherical, painless, pain-free hip a 34% five-year
   probability of arthroplasty. The hazard is now proportional to disease, and a
   healthy hip is now an exact fixed point: Q = 1.000, P_io = 16.0 mmHg, oedema
   0, pain 0, Harris Hip 100.0, P(THA) 0.000.
7. **The steroid schedule could not affect the outcome.** Incident scenarios
   never fed the computed necrotic volume back into the lesion geometry, so they
   inherited whatever `alpha_deg` was in the parameter list. Fixed with an
   explicit two-stage simulation: the insult decides how much infarcts, then
   *that* cone carries the load.
8. **Alendronate looked 85% protective.** Drug comparisons used unrestricted
   maxima that included post-collapse stress amplification, comparing a collapsed
   control with an intact treated hip. Restricting peaks to before collapse
   changed the effect from −84.9% to −63.8% and the rate ratio from 0.001 to
   0.017.
9. **The sensitivity analysis returned zero for twelve of eighteen parameters**,
   because five-year depression saturates at its cap for every hip that collapses.
10. **Time-to-collapse carried no information.** The crack-coalescence term was
    strong enough to make the damage law explosive, so every hip above the
    endurance limit failed within weeks. Weakening it separates the two
    questions the model is asked: *whether* a hip collapses is the integral of
    the damage rate over the window; *when* is its magnitude.

---

## Running it

```bash
# mechanistic map
dot -Tsvg onfh_qsp_model.dot -o onfh_qsp_model.svg
dot -Tpng -Gdpi=150 onfh_qsp_model.dot -o onfh_qsp_model.png

# the executed reference implementation (writes onfh_reference_output.txt)
pip install numpy scipy
python3 onfh_python_reference.py

# regenerate the reference list from PubMed
python3 mkrefs.py
```

```r
# the mrgsolve model
install.packages(c("mrgsolve", "dplyr", "ggplot2"))
library(mrgsolve)
mod <- mread("onfh_mrgsolve_model.R")
mod |> param(alpha_deg = 45, phi_deg = -12) |> mrgsim(end = 1826, delta = 1) |>
  plot(DEPR + SR_INT + BETA_EFF ~ time)

# the dashboard
install.packages(c("shiny", "tidyr", "DT"))
shiny::runApp("onfh_shiny_app.R")
```

---

## ⚠️ 면책 조항

본 모델은 **교육 및 연구 목적의 QSP 모델**입니다. 공개 문헌을 근거로
구성되었으나 환자 데이터에 대해 독립적으로 검증되지 않았으며, **임상 의사결정,
처방, 규제 제출에 사용해서는 안 됩니다.** 위의 "Known misses" 절에 모델이
문헌과 어긋나는 지점을 그대로 적어 두었습니다.
