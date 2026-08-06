#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build gad_references.md from PubMed.

Every title, journal, year, author and PMID below is resolved live from NCBI
esummary -- nothing in the bibliography is written from memory.  Each entry is
a (section, intent, query) triple; the top PubMed hit for the query is taken,
duplicates are dropped, and the resolved record is written out together with
the INTENT so a reader can see immediately whether the retrieved paper is
actually the one the model leans on.
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import OrderedDict

B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
CACHE = "refs_raw.json"


def _get(url):
    for k in range(6):
        try:
            return urllib.request.urlopen(url, timeout=45).read().decode("utf-8", "replace")
        except Exception:
            time.sleep(2 + 2 * k)
    return ""


def esearch(term, n=1):
    # NOTE: E-utilities defaults to date order, NOT relevance.  Without
    # sort=relevance the top hit is simply the newest paper matching any
    # of the terms, which is how an earlier build of this file ended up
    # citing a 2026 alopecia meta-analysis for "anxiety and cardiovascular
    # risk".  Do not remove it.
    u = B + ("esearch.fcgi?db=pubmed&retmode=json&sort=relevance&retmax=%d&term=%s"
             % (n, urllib.parse.quote(term)))
    try:
        return json.loads(_get(u))["esearchresult"]["idlist"]
    except Exception:
        return []


def esummary(pmids):
    out = {}
    for i in range(0, len(pmids), 150):
        chunk = pmids[i:i + 150]
        u = B + "esummary.fcgi?db=pubmed&retmode=json&id=" + ",".join(chunk)
        try:
            r = json.loads(_get(u))["result"]
        except Exception:
            continue
        for pid in chunk:
            if pid in r:
                out[pid] = r[pid]
        time.sleep(0.34)
    return out


SECTIONS = OrderedDict([
 ("1", ("The disorder itself, and why it is a gain and not a state",
        "Sustained, unpredictable-threat anxiety with a self-reinforcing worry "
        "engine. The BNST rather than the central amygdala is the node that "
        "matches the phenomenology, which is why the model gives E_amy a slow "
        "time constant and lets CBT rather than a drug move it.")),
 ("2", ("The threat circuit: E_amy in the numerator, C_pfc in the denominator",
        "The two circuit factors of Phi. Everything in this section constrains "
        "either the amygdala/BNST drive or the prefrontal control capacity.")),
 ("3", ("Serotonin: transporter occupancy, the 5-HT1A autoreceptor gate, and why the clock is slow",
        "The single most important calibration anchor in the model. SERT "
        "occupancy is >80% within days; the clinical effect is not. The gap is "
        "the autoreceptor.")),
 ("4", ("Noradrenaline and NET occupancy: the arm that helps and harms at once",
        "Venlafaxine's measured NET occupancy sets the SNRI dose-response, and "
        "the same NE rise produces the early activation that drives dropout.")),
 ("5", ("GABA-A: subunit pharmacology, benzodiazepine-site occupancy, tolerance and rebound",
        "Where I_gaba, R_a1, R_a2 and DEPEND come from. The measured lorazepam "
        "occupancy EC50 is the anchor; the alpha1/alpha5-versus-alpha2/alpha3 "
        "split is why sedation tolerates in days and anxiolysis does not.")),
 ("6", ("Glutamate, alpha2delta-1 and the pregabalin mechanism",
        "Where S_glu and its ~2 day trafficking delay come from.")),
 ("7", ("HPA axis and cortisol in anxiety",
        "The loop the model had to bound. Cortisol in GAD is modestly raised, "
        "not grossly so, and the model is held to that.")),
 ("8", ("Neuroplasticity, BDNF and extinction learning",
        "The mechanism by which a transporter block becomes a change in "
        "prefrontal control weeks later, and the reason CBT's gain outlasts "
        "the drug's.")),
 ("9", ("Autonomic function, heart-rate variability and sleep",
        "The somatic half of HAM-A, and the sleep-to-amygdala loop that is "
        "inside the disease rather than downstream of it.")),
 ("10", ("SSRI and SNRI randomised trials in GAD",
         "Out-of-sample prediction targets. None of these were fitted.")),
 ("11", ("Pregabalin randomised trials in GAD",
         "Rickels 2005 supplies two of the five fitted numbers (pregabalin "
         "300 mg and the placebo arm at week 4). The 450 and 600 mg arms, and "
         "the psychic/somatic split, are predictions.")),
 ("12", ("Benzodiazepines: efficacy, tolerance, dependence and discontinuation",
         "Rickels 2005's alprazolam arm supplies the fifth fitted number. "
         "Everything about tolerance and rebound is a prediction constrained "
         "by this literature.")),
 ("13", ("Quetiapine, buspirone, hydroxyzine and the rest of the armamentarium",
         "The quetiapine XR dose-response is the model's sharpest out-of-sample "
         "test, because the reported curve is non-monotonic.")),
 ("14", ("Relapse prevention and long-term treatment",
         "The Allgulander randomised-withdrawal design is the money test: 19% "
         "versus 56% relapse, predicted from an acute-phase calibration.")),
 ("15", ("Placebo response, assay sensitivity and trial methodology",
         "The fifth clock. The model treats expectancy as real top-down "
         "control, which is why assay sensitivity falls out of the algebra.")),
 ("16", ("Network meta-analyses, comparative efficacy and guidelines",
         "The benchmark the whole set of predictions is scored against.")),
 ("17", ("Psychological treatment",
         "CBT enters the model at three places at once and is the only "
         "intervention that moves E_amy.")),
 ("18", ("The Hamilton Anxiety Rating Scale as a measuring instrument",
         "The readout is not transparent. Its psychic/somatic structure and "
         "its ceiling are load-bearing parts of the model.")),
 ("19", ("Clinical pharmacokinetics of the seven modelled drugs",
         "Every PK parameter in the model traces to this section.")),
 ("20", ("Quantitative systems pharmacology and PK/PD methodology",
         "The modelling machinery itself.")),
])

# (section, intent, query)
ENTRIES = [
 # --- 1 ------------------------------------------------------------------
 ("1", "Definition, epidemiology and course of GAD",
  "generalized anxiety disorder epidemiology course review Lancet"),
 ("1", "Lifetime prevalence and comorbidity in national surveys",
  "Lifetime prevalence age-of-onset distributions DSM-IV disorders National Comorbidity Survey Replication"),
 ("1", "Worry as the diagnostic core; intolerance of uncertainty",
  "Dugas intolerance of uncertainty problem orientation worry generalized anxiety disorder"),
 ("1", "Sustained versus phasic threat: BNST and the anxiety/fear distinction",
  "Alvarez phasic and sustained fear humans bed nucleus of the stria terminalis amygdala neuroimaging"),
 ("1", "Disability and economic burden of GAD",
  "Wittchen generalized anxiety disorder prevalence burden cost of illness"),
 ("1", "Diagnostic boundary with major depression",
  "Kendler major depression generalised anxiety disorder same genes partly different environments"),
 # --- 2 ------------------------------------------------------------------
 ("2", "Amygdala hyperreactivity in GAD",
  "amygdala hyperactivation generalized anxiety disorder emotional faces fMRI"),
 ("2", "Amygdala-prefrontal connectivity deficits in GAD",
  "amygdala prefrontal functional connectivity generalized anxiety disorder"),
 ("2", "vmPFC top-down inhibition of amygdala output",
  "ventromedial prefrontal cortex inhibition amygdala intercalated cells extinction"),
 ("2", "BNST and unpredictable threat anticipation in humans",
  "bed nucleus of the stria terminalis anticipation unpredictable threat human neuroimaging"),
 ("2", "Insula and interoception in anxiety",
  "anterior insula interoception anxiety disorders neuroimaging"),
 ("2", "Default mode network and worry",
  "default mode network worry rumination generalized anxiety disorder resting state"),
 ("2", "Startle potentiation to unpredictable threat as a GAD laboratory model",
  "startle potentiation unpredictable threat anxiety disorders NPU"),
 ("2", "Structural changes in amygdala and prefrontal cortex in GAD",
  "grey matter volume amygdala prefrontal generalized anxiety disorder voxel-based"),
 # --- 3 ------------------------------------------------------------------
 ("3", "SERT occupancy versus plasma concentration for SSRIs (the 80% threshold)",
  "Hart lessons from PET imaging therapeutic plasma concentrations antidepressants serotonin transporter occupancy"),
 ("3", "Escitalopram SERT occupancy by [11C]DASB and its regional differences",
  "Kim regional differences in serotonin transporter occupancy by escitalopram PK-PD DASB"),
 ("3", "Dose-occupancy relationship for SSRIs at clinical doses",
  "Meyer serotonin transporter occupancies in patients with major depression administered paroxetine citalopram DASB"),
 ("3", "5-HT1A somatodendritic autoreceptor desensitisation with chronic SSRI",
  "chronic SSRI 5-HT1A somatodendritic autoreceptor desensitization raphe firing"),
 ("3", "Time course of extracellular 5-HT rise with chronic SSRI (microdialysis)",
  "chronic fluoxetine increases extracellular serotonin frontal cortex microdialysis rat"),
 ("3", "Pindolol augmentation as a test of the autoreceptor hypothesis",
  "pindolol augmentation selective serotonin reuptake inhibitors meta-analysis depression"),
 ("3", "Delayed onset of antidepressant/anxiolytic action: the pattern to be explained",
  "onset of action antidepressant delayed weeks pattern early improvement meta-analysis"),
 ("3", "5-HT1A receptor binding in anxiety disorders",
  "5-HT1A receptor binding PET anxiety disorder panic reduced"),
 ("3", "5-HT3 mediated nausea with SSRIs",
  "adverse effects of selective serotonin reuptake inhibitors nausea sexual dysfunction review"),
 # --- 4 ------------------------------------------------------------------
 ("4", "Venlafaxine ER norepinephrine transporter occupancy in patients (PET)",
  "venlafaxine extended release blocks norepinephrine transporter PET FMeNER major depressive"),
 ("4", "Duloxetine SERT and NET occupancy",
  "duloxetine serotonin norepinephrine transporter occupancy PET humans"),
 ("4", "Locus coeruleus and noradrenergic contribution to anxiety",
  "locus coeruleus noradrenergic system anxiety stress arousal review"),
 ("4", "Prefrontal noradrenaline inverted-U and alpha1/alpha2 receptor actions",
  "noradrenaline prefrontal cortex inverted U alpha1 alpha2 receptor working memory Arnsten"),
 ("4", "SNRI early activation / jitteriness and blood pressure",
  "venlafaxine sustained hypertension dose related blood pressure increase"),
 # --- 5 ------------------------------------------------------------------
 ("5", "Lorazepam plasma EC50 for GABA-A benzodiazepine-site occupancy",
  "lorazepam occupancy GABA-A receptors flumazenil micro-positron emission tomography rat"),
 ("5", "Benzodiazepine receptor occupancy at clinical doses in humans",
  "benzodiazepine receptor occupancy humans PET clinical doses flumazenil displacement"),
 ("5", "alpha2/alpha3 subunits mediate anxiolysis; alpha1 mediates sedation",
  "alpha2 GABAA receptors mediate anxiolytic action of benzodiazepines knock-in point mutation"),
 ("5", "alpha1 subunit mediates the sedative action of diazepam",
  "alpha1 GABAA receptor sedative action diazepam knock-in histidine arginine"),
 ("5", "alpha5-GABA-A receptors are required for tolerance to sedation",
  "requirement alpha5 GABAA receptors development tolerance sedative diazepam"),
 ("5", "Differential tolerance to sedative versus anxiolytic benzodiazepine effects",
  "differential tolerance sedative anxiolytic effects benzodiazepines chronic treatment"),
 ("5", "Cortical GABA concentration in anxiety disorders (1H-MRS)",
  "magnetic resonance spectroscopy cortical GABA anxiety disorder reduced"),
 ("5", "Neurosteroids, allopregnanolone and GABA-A tonic inhibition",
  "allopregnanolone neurosteroid GABAA receptor tonic inhibition anxiety"),
 ("5", "Benzodiazepine physical dependence and the withdrawal syndrome",
  "benzodiazepine physical dependence withdrawal syndrome mechanisms review"),
 # --- 6 ------------------------------------------------------------------
 ("6", "alpha2delta-1 as the pregabalin binding site",
  "identification of the alpha2-delta-1 subunit of voltage-dependent calcium channels as a molecular target for pain"),
 ("6", "alpha2delta trafficking as the mechanism of delayed onset",
  "gabapentin alpha2delta trafficking calcium channel inhibition trafficking mechanism"),
 ("6", "Pregabalin reduces excitatory neurotransmitter release",
  "pregabalin inhibits release excitatory neurotransmitter glutamate cortical synaptosomes"),
 ("6", "Glutamate abnormalities in anxiety disorders",
  "glutamate anxiety disorders proton magnetic resonance spectroscopy meta-analysis"),
 ("6", "Stress upregulates alpha2delta-1 and presynaptic release probability",
  "acute stress enhances glutamate release prefrontal cortex corticosterone Popoli"),
 # --- 7 ------------------------------------------------------------------
 ("7", "Cortisol in generalized anxiety disorder",
  "cortisol generalized anxiety disorder salivary diurnal elevated"),
 ("7", "CRH, CRHR1 and anxiety behaviour",
  "corticotropin releasing hormone CRHR1 amygdala anxiety behaviour"),
 ("7", "Glucocorticoid receptor resistance under chronic stress",
  "chronic stress glucocorticoid receptor resistance and disease risk Cohen"),
 ("7", "Dexamethasone suppression and HPA feedback in anxiety",
  "dexamethasone suppression test panic disorder generalized anxiety nonsuppression"),
 ("7", "Cortisol and prefrontal function",
  "cortisol administration impairs working memory prefrontal cortex humans"),
 # --- 8 ------------------------------------------------------------------
 ("8", "BDNF and antidepressant-induced plasticity",
  "Castren neuronal plasticity and network plasticity antidepressant action"),
 ("8", "Fear extinction learning: circuitry and consolidation",
  "neuroscience of fear extinction implications for assessment and treatment of fear-based anxiety disorders"),
 ("8", "BDNF Val66Met and impaired extinction",
  "BDNF Val66Met polymorphism impaired fear extinction humans"),
 ("8", "Chronic stress causes opposite dendritic remodelling in PFC and amygdala",
  "chronic stress induces dendritic hypertrophy amygdala atrophy medial prefrontal cortex Vyas McEwen"),
 ("8", "Inflammation, kynurenine and reduced BDNF",
  "inflammation interleukin-6 kynurenine pathway BDNF depression anxiety"),
 # --- 9 ------------------------------------------------------------------
 ("9", "Reduced heart rate variability in anxiety disorders (meta-analysis)",
  "heart rate variability anxiety disorders meta-analysis reduced vagal"),
 ("9", "Autonomic inflexibility in GAD",
  "Thayer autonomic characteristics generalized anxiety disorder worry"),
 ("9", "Insomnia and anxiety: bidirectional relationship",
  "insomnia as a predictor of subsequent anxiety disorder longitudinal epidemiological study"),
 ("9", "Sleep deprivation amplifies amygdala reactivity",
  "the human emotional brain without sleep a prefrontal amygdala disconnect"),
 ("9", "Anxiety and cardiovascular risk",
  "Roest anxiety and risk of incident coronary heart disease meta-analysis"),
 # --- 10 -----------------------------------------------------------------
 ("10", "Escitalopram fixed-dose efficacy in GAD",
  "Baldwin escitalopram paroxetine generalised anxiety disorder randomised placebo-controlled double-blind"),
 ("10", "Paroxetine in GAD",
  "Rickels paroxetine treatment of generalized anxiety disorder double-blind placebo-controlled flexible-dose"),
 ("10", "Sertraline in GAD",
  "Allgulander sertraline treatment of generalized anxiety disorder randomized placebo-controlled"),
 ("10", "Venlafaxine ER in GAD: short and long term",
  "efficacy of venlafaxine extended-release capsules in nondepressed outpatients with generalized anxiety disorder six-month"),
 ("10", "Venlafaxine ER in older adults with GAD",
  "venlafaxine ER treatment generalized anxiety disorder older adults pooled analysis"),
 ("10", "Duloxetine in GAD",
  "Carter duloxetine a review of its use in the treatment of generalized anxiety disorder"),
 ("10", "Duloxetine versus venlafaxine in GAD",
  "duloxetine venlafaxine generalized anxiety disorder active comparator randomized"),
 ("10", "Antidepressant dose-response is flat for SSRIs",
  "Furukawa optimal dose of selective serotonin reuptake inhibitors dose-response meta-analysis"),
 # --- 11 -----------------------------------------------------------------
 ("11", "Rickels 2005: pregabalin 300/450/600 versus alprazolam and placebo",
  "pregabalin treatment generalized anxiety disorder 4-week multicenter double-blind alprazolam Rickels"),
 ("11", "Pregabalin review: onset within one week, psychic and somatic",
  "Frampton pregabalin a review of its use in adults with generalized anxiety disorder"),
 ("11", "Pregabalin efficacy and safety profile in generalized anxiety",
  "pregabalin efficacy safety tolerability profile generalized anxiety"),
 ("11", "Pregabalin in elderly patients with GAD",
  "Montgomery pregabalin elderly patients with generalized anxiety disorder double-blind placebo-controlled"),
 ("11", "Pregabalin relapse prevention in GAD",
  "pregabalin long-term treatment relapse prevention generalized anxiety disorder"),
 ("11", "Pregabalin augmentation of SSRI/SNRI partial responders",
  "pregabalin adjunctive treatment generalized anxiety disorder partial response antidepressant"),
 # --- 12 -----------------------------------------------------------------
 ("12", "Alprazolam clinical pharmacology, efficacy and behavioural toxicity",
  "clinical pharmacology clinical efficacy behavioral toxicity alprazolam review"),
 ("12", "Benzodiazepines versus antidepressants in GAD",
  "benzodiazepines antidepressants generalized anxiety disorder comparison efficacy meta-analysis"),
 ("12", "Long-term benzodiazepine use: benefits and harms",
  "long-term use of benzodiazepines in anxiety disorders risks benefits review"),
 ("12", "Benzodiazepine discontinuation strategies and rebound anxiety",
  "discontinuation of long-term benzodiazepine use taper rebound anxiety randomized controlled"),
 ("12", "Benzodiazepines, falls and fractures in older adults",
  "benzodiazepine use and risk of hip fracture in older adults"),
 ("12", "Cognitive and psychomotor impairment with benzodiazepines",
  "benzodiazepines cognitive psychomotor impairment driving meta-analysis"),
 # --- 13 -----------------------------------------------------------------
 ("13", "Khan 2011: quetiapine XR 50/150/300 mg monotherapy in GAD",
  "randomized double-blind once-daily extended release quetiapine fumarate monotherapy generalized anxiety disorder Khan"),
 ("13", "Quetiapine XR in GAD: meta-analysis",
  "quetiapine monotherapy acute treatment generalized anxiety disorder systematic review meta-analysis"),
 ("13", "Buspirone in generalized anxiety",
  "buspirone in the treatment of generalized anxiety disorder efficacy review"),
 ("13", "Azapirones for generalised anxiety disorder (Cochrane)",
  "azapirones versus placebo for generalised anxiety disorder Cochrane"),
 ("13", "Hydroxyzine for generalised anxiety disorder",
  "hydroxyzine for generalised anxiety disorder Cochrane review"),
 ("13", "Agomelatine in GAD",
  "agomelatine generalized anxiety disorder placebo controlled active comparator"),
 ("13", "Vortioxetine in GAD",
  "vortioxetine Lu AA21004 generalized anxiety disorder randomized double-blind placebo-controlled"),
 ("13", "Antipsychotic augmentation in refractory GAD",
  "quetiapine augmentation refractory generalized anxiety disorder randomized"),
 # --- 14 -----------------------------------------------------------------
 ("14", "Allgulander 2006: escitalopram relapse prevention in GAD",
  "Allgulander prevention of relapse in generalized anxiety disorder by escitalopram treatment"),
 ("14", "Duloxetine relapse prevention in GAD",
  "duloxetine prevention relapse generalized anxiety disorder randomized withdrawal"),
 ("14", "Randomised-withdrawal designs in anxiety disorders: meta-analysis",
  "Batelaan risk of relapse after antidepressant discontinuation anxiety disorders meta-analysis"),
 ("14", "How long should treatment continue after remission",
  "optimal duration of continuation treatment anxiety disorders after remission"),
 ("14", "Discontinuation symptoms after stopping SSRIs and SNRIs",
  "incidence of antidepressant discontinuation symptoms systematic review meta-analysis"),
 # --- 15 -----------------------------------------------------------------
 ("15", "Placebo response in anxiety disorder trials",
  "placebo response in generalized anxiety disorder trials predictors magnitude"),
 ("15", "Rising placebo response and failed trials in psychiatry",
  "increasing placebo response in antidepressant trials assay sensitivity failed trials"),
 ("15", "Number of study visits predicts placebo response",
  "frequency of study visits and placebo response in antidepressant clinical trials"),
 ("15", "Regression to the mean and baseline severity inflation in trials",
  "baseline score inflation clinical trials antidepressant eligibility criteria"),
 ("15", "Expectancy mechanisms of placebo analgesia and anxiolysis",
  "Benedetti neurobiological mechanisms of the placebo effect expectation"),
 ("15", "Baseline severity and drug-placebo difference",
  "initial severity and antidepressant benefit meta-analysis Kirsch Fournier"),
 # --- 16 -----------------------------------------------------------------
 ("16", "Network meta-analysis of first-line drugs for GAD",
  "comparative efficacy acceptability first-line drugs generalized anxiety disorder network meta-analysis"),
 ("16", "Remission rates and tolerability network meta-analysis in GAD",
  "comparative remission rates tolerability drugs generalised anxiety disorder network meta-analysis"),
 ("16", "Anxiolytic drugs across anxiety disorders: network meta-analysis",
  "comparative efficacy acceptability anxiolytic drugs anxiety disorders systematic review network meta-analysis"),
 ("16", "WFSBP guidelines for anxiety disorders",
  "Bandelow WFSBP guidelines pharmacological treatment anxiety obsessive-compulsive post-traumatic stress disorders"),
 ("16", "NICE / national guidance on GAD management",
  "Baldwin evidence-based pharmacological treatment anxiety disorders British Association Psychopharmacology guidelines"),
 # --- 17 -----------------------------------------------------------------
 ("17", "CBT for GAD: meta-analysis of effect size",
  "cognitive behavioural therapy generalised anxiety disorder randomised controlled trials meta-analysis Hanrahan Cuijpers"),
 ("17", "CBT for pathological worry",
  "meta-analysis of CBT for pathological worry among clients with GAD"),
 ("17", "Combined psychotherapy and pharmacotherapy in anxiety disorders",
  "combined psychotherapy and antidepressants for anxiety disorders meta-analysis Cuijpers"),
 ("17", "Internet-delivered CBT for anxiety",
  "guided internet-delivered cognitive behaviour therapy for anxiety disorders meta-analysis Andersson"),
 ("17", "Long-term durability of CBT gains in anxiety",
  "long-term outcome of cognitive behaviour therapy for anxiety disorders follow-up"),
 ("17", "Exercise for anxiety disorders",
  "meta-analysis aerobic exercise treatment anxiety disorders"),
 # --- 18 -----------------------------------------------------------------
 ("18", "The Hamilton Anxiety Rating Scale: original description",
  "Hamilton M the assessment of anxiety states by rating 1959"),
 ("18", "HAM-A psychometric properties and factor structure",
  "Maier Hamilton Anxiety Scale reliability validity sensitivity to change anxiety disorders"),
 ("18", "Clinically meaningful change on the HAM-A",
  "Hamilton Anxiety Rating Scale clinically meaningful change threshold generalized anxiety disorder"),
 ("18", "GAD-7 development and validation",
  "Spitzer a brief measure for assessing generalized anxiety disorder the GAD-7"),
 ("18", "Penn State Worry Questionnaire",
  "Meyer development and validation of the Penn State Worry Questionnaire"),
 ("18", "Limits of change-score endpoints in psychiatry trials",
  "psychometric limitations of the Hamilton scales unidimensionality item response"),
 # --- 19 -----------------------------------------------------------------
 ("19", "Escitalopram clinical pharmacokinetics",
  "escitalopram clinical pharmacokinetics review steady state half-life"),
 ("19", "CYP2C19 genotype and escitalopram exposure",
  "CYP2C19 genotype escitalopram serum concentration poor metabolizer"),
 ("19", "Venlafaxine and O-desmethylvenlafaxine pharmacokinetics; CYP2D6",
  "venlafaxine O-desmethylvenlafaxine pharmacokinetics CYP2D6 metabolizer status"),
 ("19", "Duloxetine pharmacokinetics",
  "duloxetine clinical pharmacokinetics review CYP1A2 CYP2D6 half-life"),
 ("19", "Pregabalin pharmacokinetics and renal dose adjustment",
  "pregabalin pharmacokinetics renal impairment dose adjustment creatinine clearance"),
 ("19", "Lorazepam pharmacokinetics",
  "Greenblatt clinical pharmacokinetics of lorazepam"),
 ("19", "Alprazolam pharmacokinetics",
  "Greenblatt alprazolam pharmacokinetics metabolism single dose humans"),
 ("19", "Buspirone pharmacokinetics and 1-PP metabolite",
  "buspirone pharmacokinetics 1-pyrimidinylpiperazine metabolite plasma concentrations humans"),
 ("19", "Quetiapine and norquetiapine pharmacokinetics",
  "quetiapine and norquetiapine plasma concentrations therapeutic drug monitoring"),
 ("19", "Therapeutic reference ranges for antidepressants (TDM consensus)",
  "Hiemke consensus guidelines for therapeutic drug monitoring in neuropsychopharmacology"),
 # --- 20 -----------------------------------------------------------------
 ("20", "Quantitative systems pharmacology: scope and role in drug development",
  "quantitative systems pharmacology models current perspective drug development white paper"),
 ("20", "mrgsolve for ODE-based PK/PD simulation in R",
  "Elmokadem quantitative systems pharmacology and physiologically-based pharmacokinetic modeling with mrgsolve"),
 ("20", "Model-informed drug development in CNS",
  "model-informed drug development regulatory review pharmacometrics FDA"),
 ("20", "Receptor occupancy modelling and the link to clinical effect",
  "receptor occupancy imaging dose selection CNS drug development PET translational"),
 ("20", "Clinical trial simulation with placebo and dropout models",
  "clinical trial simulation placebo response dropout model longitudinal exposure response"),
 ("20", "Disease progression modelling in neuropsychiatry",
  "disease progression model pharmacometric analysis Alzheimer ADAS-cog longitudinal"),
 ("20", "Item-response / latent-variable models for psychiatric rating scales",
  "item response theory pharmacometric modelling rating scale total score information"),
]


def main():
    if os.path.exists(CACHE) and "--refresh" not in sys.argv:
        raw = json.load(open(CACHE))
    else:
        raw = {}
    todo = [e for e in ENTRIES if e[2] not in raw]
    print("resolving %d queries (%d cached)" % (len(todo), len(ENTRIES) - len(todo)))
    for sec, intent, q in todo:
        ids = esearch(q, 1)
        raw[q] = ids
        print("  %-70s -> %s" % (q[:70], ids))
        time.sleep(0.34)
    json.dump(raw, open(CACHE, "w"), indent=0)

    pmids = []
    for _, _, q in ENTRIES:
        for p in raw.get(q, []):
            if p not in pmids:
                pmids.append(p)
    meta = esummary(pmids)
    json.dump(meta, open("refs_meta.json", "w"))

    # ---- write markdown -------------------------------------------------
    seen = set()
    lines = []
    lines.append("# Generalized Anxiety Disorder QSP model — references\n")
    lines.append(
        "Every entry below was resolved live from NCBI PubMed (`esearch` + "
        "`esummary`) by [`mkrefs.py`](mkrefs.py); titles, journals, years, "
        "author lists and PMIDs are therefore machine-transcribed rather than "
        "recalled. Each entry carries the **intent** — what the model actually "
        "takes from it — so that a retrieved paper that does not match its "
        "intent is visible at a glance rather than hidden.\n")
    lines.append(
        "> **Five numbers in this repository were fitted; everything else was "
        "predicted.** The fitted five are the placebo-arm trajectory of Khan "
        "2011 (week 1 −5.94, week 8 −11.10), and the week-4 endpoints of "
        "Rickels 2005 (placebo −8.4, pregabalin 300 mg −12.2, alprazolam "
        "1.5 mg −10.9). Sections 10–16 are, with those exceptions, "
        "out-of-sample prediction targets.\n")
    n = 0
    for sec, (title, blurb) in SECTIONS.items():
        rows = [e for e in ENTRIES if e[0] == sec]
        lines.append("\n## %s. %s\n" % (sec, title))
        lines.append("*%s*\n" % blurb)
        for _, intent, q in rows:
            ids = raw.get(q, [])
            if not ids:
                lines.append("- **%s** — *no PubMed record retrieved for this "
                             "query; stated as an unresolved gap rather than "
                             "cited from memory.*" % intent)
                continue
            pid = ids[0]
            m = meta.get(pid)
            if not m:
                lines.append("- **%s** — PMID %s (metadata unavailable)" % (intent, pid))
                continue
            auth = m.get("authors", [])
            if len(auth) == 0:
                a = ""
            elif len(auth) == 1:
                a = auth[0]["name"]
            elif len(auth) <= 3:
                a = ", ".join(x["name"] for x in auth)
            else:
                a = auth[0]["name"] + " et al."
            yr = (m.get("pubdate", "") or "")[:4]
            src = m.get("source", "")
            ttl = re.sub(r"\s+", " ", m.get("title", "")).strip().rstrip(".")
            dup = " *(also cited above)*" if pid in seen else ""
            seen.add(pid)
            n += 1
            lines.append("%d. **%s**  \n   %s%s %s. *%s* %s. "
                         "[PMID %s](https://pubmed.ncbi.nlm.nih.gov/%s/)%s"
                         % (n, intent, a, "" if a.endswith(".") else ".", ttl, src, yr, pid, pid, dup))
        lines.append("")
    lines.append("\n---\n")
    lines.append("*%d bibliography entries, %d distinct PubMed records. "
                 "Regenerate with `python3 mkrefs.py` (add `--refresh` to "
                 "re-run every query).*" % (n, len(seen)))
    open("gad_references.md", "w").write("\n".join(lines) + "\n")
    print("wrote gad_references.md: %d entries, %d distinct PMIDs" % (n, len(seen)))


if __name__ == "__main__":
    main()
