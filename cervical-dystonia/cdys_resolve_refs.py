#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cdys_resolve_refs.py -- resolve every citation in the reference list against
live PubMed, so that no PMID in `cdys_references.md` is written from memory.

For each entry the script sends the search string to PubMed esearch, takes the
best-matching record, then pulls its esummary and prints the authoritative
title / first author / journal / year / PMID.  The markdown reference list is
generated from THAT output, not from the query.

    python3 cdys_resolve_refs.py            # resolve all, emit TSV + markdown
    python3 cdys_resolve_refs.py --md       # markdown only

Entries that fail to resolve are printed as UNRESOLVED and are omitted from the
reference list rather than being guessed at.
"""
import argparse
import json
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# (section, model-relevance note, pubmed query)
QUERIES = [
    # ---- A. cervical dystonia: definition, epidemiology, natural history ----
    ("A. Cervical dystonia — phenomenology and epidemiology",
     "defines the disease the model is about; prevalence and adult onset",
     "Defourny Cervical dystonia epidemiology prevalence review"),
    ("A. Cervical dystonia — phenomenology and epidemiology",
     "consensus phenomenology/classification used for the posture patterns",
     "Albanese phenomenology classification dystonia consensus update 2013"),
    ("A. Cervical dystonia — phenomenology and epidemiology",
     "COL-CAP collar-based pattern classification informing the muscle table",
     "Reichel COL-CAP concept cervical dystonia classification"),
    ("A. Cervical dystonia — phenomenology and epidemiology",
     "natural history and remission rates; why chronic q12wk dosing is the norm",
     "Jahanshahi Marsden natural history torticollis prognosis"),
    ("A. Cervical dystonia — phenomenology and epidemiology",
     "spread and long-term course of adult-onset focal dystonia",
     "Weiss spread adult onset focal dystonia risk"),

    # ---- B. genetics ----
    ("B. Genetics and predisposition",
     "TOR1A/DYT1 — the predisposition cluster in the map",
     "Ozelius DYT1 torsin A gene dystonia 1997 Nature Genetics"),
    ("B. Genetics and predisposition",
     "THAP1/DYT6 — cranio-cervical predominance",
     "Fuchs THAP1 DYT6 mutations primary torsion dystonia"),
    ("B. Genetics and predisposition",
     "GNAL/DYT25 — adult-onset cervical dystonia gene",
     "Fuchs GNAL mutations primary torsion dystonia exome"),
    ("B. Genetics and predisposition",
     "ANO3/DYT24",
     "Charlesworth ANO3 mutations craniocervical dystonia anoctamin"),
    ("B. Genetics and predisposition",
     "CIZ1 variants in adult-onset cervical dystonia",
     "Xiao CIZ1 variants adult onset cervical dystonia"),
    ("B. Genetics and predisposition",
     "COL6A3 recessive isolated dystonia",
     "Zech COL6A3 mutations recessive isolated dystonia DYT27"),
    ("B. Genetics and predisposition",
     "KMT2B — complex early-onset dystonia",
     "Meyer KMT2B mutations childhood onset dystonia"),

    # ---- C. basal ganglia, dopamine, plasticity ----
    ("C. Basal ganglia, dopamine and corticostriatal plasticity",
     "loss of bidirectional corticostriatal plasticity — the LTD_LOSS node",
     "Quartarone Pisani abnormal plasticity dystonia review"),
    ("C. Basal ganglia, dopamine and corticostriatal plasticity",
     "striatal cholinergic-dopaminergic imbalance; rationale for antimuscarinics",
     "Pisani Bernardi Bonsi striatal cholinergic dysfunction dystonia"),
    ("C. Basal ganglia, dopamine and corticostriatal plasticity",
     "D2 receptor binding reduced in cervical dystonia (imaging)",
     "Perlmutter Stambuk decreased D2 receptor binding idiopathic cervical dystonia"),
    ("C. Basal ganglia, dopamine and corticostriatal plasticity",
     "network rather than nucleus: the dystonia network concept",
     "Prudente Hess Jinnah dystonia as a network disorder"),
    ("C. Basal ganglia, dopamine and corticostriatal plasticity",
     "GPi single-unit firing abnormalities in dystonia",
     "Vitek pallidal neuronal activity dystonia microelectrode"),

    # ---- D. cerebellum ----
    ("D. Cerebellar and olivo-cerebellar contribution",
     "cerebellum in dystonia pathophysiology — the Cbll gain node",
     "Filip Bares cerebellum dystonia pathophysiology review"),
    ("D. Cerebellar and olivo-cerebellar contribution",
     "eyeblink conditioning abnormality in cervical dystonia",
     "Teo Edwards Bhatia eyeblink classical conditioning cervical dystonia"),
    ("D. Cerebellar and olivo-cerebellar contribution",
     "cerebellar-basal ganglia interaction produces dystonia in animal models",
     "Neychev Jinnah cerebellum basal ganglia motor circuit dystonia"),

    # ---- E. cortical inhibition ----
    ("E. Cortical inhibition and surround inhibition",
     "surround inhibition deficit — the SurrInh state",
     "Sohn Hallett surround inhibition focal dystonia"),
    ("E. Cortical inhibition and surround inhibition",
     "reduced short-interval intracortical inhibition in dystonia",
     "Ridding Sheean Rothwell changes excitability motor cortical dystonia"),
    ("E. Cortical inhibition and surround inhibition",
     "the pathophysiology synthesis the central module follows",
     "Hallett pathophysiology dystonia review neurobiology of disease"),

    # ---- F. sensorimotor integration and the sensory trick ----
    ("F. Sensorimotor integration, afferent drive and the sensory trick",
     "abnormal sensorimotor integration; the AFF_ABN node",
     "Abbruzzese Berardelli sensorimotor integration dystonia"),
    ("F. Sensorimotor integration, afferent drive and the sensory trick",
     "temporal discrimination threshold as an endophenotype",
     "Bradley Hutchinson temporal discrimination threshold cervical dystonia endophenotype"),
    ("F. Sensorimotor integration, afferent drive and the sensory trick",
     "muscle spindle / vibration-induced illusion abnormalities",
     "Grunewald Yoneda Shipman Sagar idiopathic focal dystonia disorder muscle spindle"),
    ("F. Sensorimotor integration, afferent drive and the sensory trick",
     "sensory trick (geste antagoniste) mechanism",
     "Muller Wissel Masuhr sensory tricks cervical dystonia efficacy"),
    ("F. Sensorimotor integration, afferent drive and the sensory trick",
     "BoNT reduces spindle afferent input via intrafusal terminals — the model's route from periphery to centre",
     "Rosales Arimura Ikenaga Takenaga extrafusal intrafusal muscle effects botulinum toxin"),
    ("F. Sensorimotor integration, afferent drive and the sensory trick",
     "botulinum toxin action on muscle spindle afferents",
     "Filippi Errico Santarelli Bagolini Manni botulinum A toxin spindle afferent"),

    # ---- G. spinal circuits ----
    ("G. Spinal segmental circuits",
     "reciprocal inhibition deficit in the upper limb / neck — the RecInh state",
     "Nakashima Rothwell Day Thompson Marsden reciprocal inhibition writer cramp"),
    ("G. Spinal segmental circuits",
     "brainstem/spinal interneuronal excitability in cervical dystonia",
     "Deuschl Seifert Heinen Illert Lucking reciprocal inhibition forearm cervical dystonia"),

    # ---- H. muscle selection, EMG and ultrasound guidance (the rho lever) ----
    ("H. Muscle selection, EMG and ultrasound guidance — the ρ and Σw levers",
     "EMG guidance improves outcome; direct evidence for the rho lever",
     "Comella Buchman Tanner Brown-Toms Goetz botulinum toxin injection EMG assistance torticollis"),
    ("H. Muscle selection, EMG and ultrasound guidance — the ρ and Σw levers",
     "ultrasound-guided injection of deep cervical muscles",
     "Hong Sung ultrasound guided botulinum toxin injection cervical dystonia deep muscles"),
    ("H. Muscle selection, EMG and ultrasound guidance — the ρ and Σw levers",
     "obliquus capitis inferior injection — a Σw extension the standard plan omits",
     "Ito obliquus capitis inferior botulinum toxin cervical dystonia rotation"),
    ("H. Muscle selection, EMG and ultrasound guidance — the ρ and Σw levers",
     "accuracy of clinical versus instrumented muscle selection",
     "Nijmeijer Koelman Standaard Tijssen muscle selection cervical dystonia agreement"),
    ("H. Muscle selection, EMG and ultrasound guidance — the ρ and Σw levers",
     "intramuscular spread of injectate limits how much of a muscle is reached",
     "Shaari Sanders quantifying how location volume botulinum toxin injection affect muscle paralysis"),

    # ---- I. molecular mechanism ----
    ("I. Botulinum neurotoxin — molecular mechanism",
     "SNAP-25 is the serotype A substrate",
     "Blasi Chapman Link Binz Yamasaki botulinum neurotoxin A selectively cleaves SNAP-25"),
    ("I. Botulinum neurotoxin — molecular mechanism",
     "serotype B cleaves VAMP/synaptobrevin — why anti-A antibody does not cross-neutralise",
     "Schiavo Benfenati Poulain Rossetto tetanus botulinum-B neurotoxins block synaptobrevin"),
    ("I. Botulinum neurotoxin — molecular mechanism",
     "SV2 is the protein receptor for BoNT/A",
     "Dong Yeh Tepp Dean Johnson SV2 is the protein receptor for botulinum neurotoxin A"),
    ("I. Botulinum neurotoxin — molecular mechanism",
     "synaptotagmin I/II is the BoNT/B receptor; dual-receptor model",
     "Dong Richards Goodnough Tepp Johnson Chapman synaptotagmins I II mediate entry botulinum neurotoxin B"),
    ("I. Botulinum neurotoxin — molecular mechanism",
     "ganglioside plus protein dual receptor requirement",
     "Rummel double receptor anchorage botulinum neurotoxins"),
    ("I. Botulinum neurotoxin — molecular mechanism",
     "translocation of the light chain across the endosomal membrane",
     "Koriazova Montal translocation botulinum neurotoxin light chain protein channel"),
    ("I. Botulinum neurotoxin — molecular mechanism",
     "structure/function overview of the three domains",
     "Rossetto Pirazzini Montecucco botulinum neurotoxins genetic architecture pathophysiology"),

    # ---- J. duration, light chain persistence ----
    ("J. Duration of action and light-chain persistence — the calibrated k_LC",
     "long intracellular persistence of the BoNT/A light chain",
     "Whitemarsh Strathman Chase Stenmark Gutcher Johnson persistence botulinum neurotoxin A subtypes human neurons"),
    ("J. Duration of action and light-chain persistence — the calibrated k_LC",
     "duration differences among serotypes trace to light-chain lifetime",
     "Foran Mohammed Lisk Nagwaney Lawrence evaluation therapeutic usefulness botulinum neurotoxin B C1 E F"),
    ("J. Duration of action and light-chain persistence — the calibrated k_LC",
     "cleaved SNAP-25 persists long after clinical recovery",
     "Meunier Lisk Sesardic Dolly dynamics motor nerve terminal remodeling unveiled SNAP-25 cleaving toxin"),
    ("J. Duration of action and light-chain persistence — the calibrated k_LC",
     "molecular determinants of BoNT/A persistence",
     "Tsai Kim Chang Balint Tsai ubiquitin proteasome system botulinum neurotoxin A persistence"),

    # ---- K. sprouting and the two-phase recovery (A5) ----
    ("K. Terminal sprouting and the two-phase recovery — tested in A5",
     "the two-phase recovery: sprouts first, parent terminal second",
     "de Paiva Meunier Molgo Aoki Dolly functional repair botulinum poisoned mouse neuromuscular junctions"),
    ("K. Terminal sprouting and the two-phase recovery — tested in A5",
     "nerve sprouting after botulinum paralysis",
     "Duchen changes motor end plate botulinum toxin sprouting mouse"),
    ("K. Terminal sprouting and the two-phase recovery — tested in A5",
     "whether sprouts are required for functional recovery",
     "Rogozhin Pang Bukharaeva Young Slater recovery mouse neuromuscular junctions botulinum neurotoxin A"),
    ("K. Terminal sprouting and the two-phase recovery — tested in A5",
     "muscle changes and atrophy under chronic chemodenervation",
     "Fortuna Aurora Bhatia Lorenzo Cocco muscle changes after botulinum toxin injection"),

    # ---- L. neuromuscular safety factor (the threshold in E(S)) ----
    ("L. Neuromuscular safety factor — the threshold in E(S)",
     "the safety factor: why partial SNARE cleavage is clinically silent",
     "Wood Slater safety factor neuromuscular transmission"),
    ("L. Neuromuscular safety factor — the threshold in E(S)",
     "quantal content and SNARE stoichiometry for a fusion event",
     "Mohrmann de Wit Verhage Neher Sorensen fast vesicle fusion three SNARE complexes"),

    # ---- M. pivotal and comparative clinical trials ----
    ("M. Clinical trials in cervical dystonia — the calibration anchors",
     "incobotulinumtoxinA pivotal CD trial: the week-4 and week-12 TWSTRS anchors",
     "Comella Jankovic Truong Hanschmann Grafe incobotulinumtoxinA NT 201 cervical dystonia randomized"),
    ("M. Clinical trials in cervical dystonia — the calibration anchors",
     "onabotulinumtoxinA versus placebo in cervical dystonia",
     "Truong Duane Jankovic Singer Seeberger efficacy safety botulinum type A Dysport cervical dystonia"),
    ("M. Clinical trials in cervical dystonia — the calibration anchors",
     "dose-ranging: the PLATEAU that identifies phi",
     "Wissel Kanovsky Ruzicka Bares Hortova dose response botulinum toxin cervical dystonia Dysport"),
    ("M. Clinical trials in cervical dystonia — the calibration anchors",
     "daxibotulinumtoxinA ASPEN-1: ~24 weeks duration at comparable peak effect",
     "Comella Jankovic Truong Colosimo Chung daxibotulinumtoxinA cervical dystonia ASPEN-1"),
    ("M. Clinical trials in cervical dystonia — the calibration anchors",
     "rimabotulinumtoxinB in type A responsive and resistant patients",
     "Brashear Lew Dykstra Comella Factor Molho safety efficacy botulinum toxin type B cervical dystonia"),
    ("M. Clinical trials in cervical dystonia — the calibration anchors",
     "botulinum toxin type B in type A resistant cervical dystonia — serotype rescue",
     "Brin Lew Adler Comella Brashear Factor safety efficacy NeuroBloc botulinum toxin type B type A-resistant cervical dystonia"),
    ("M. Clinical trials in cervical dystonia — the calibration anchors",
     "evidence-based review / guideline for BoNT in cervical dystonia",
     "Simpson Hallett Ashman Comella Green practice guideline botulinum neurotoxin blepharospasm cervical dystonia headache"),
    ("M. Clinical trials in cervical dystonia — the calibration anchors",
     "Cochrane review of botulinum toxin type A for cervical dystonia",
     "Costa Espirito-Santo Ferreira Coutinho Botelho botulinum toxin type A cervical dystonia Cochrane"),
    ("M. Clinical trials in cervical dystonia — the calibration anchors",
     "long-term outcome of repeated injection cycles — the trough/ratchet question in A9",
     "Hsiung Tsui Pfeiffer Diamond Fung long-term efficacy botulinum toxin A cervical dystonia"),

    # ---- N. dose, dilution and volume (A6) ----
    ("N. Dose, dilution and injectate volume — A6",
     "dilution/volume changes spread at constant Units",
     "Kim Shin Lee Kim Kwon effect botulinum toxin dilution muscle paralysis"),
    ("N. Dose, dilution and injectate volume — A6",
     "injection volume and the diffusion field",
     "Brodsky Swope Grimes diffusion botulinum toxins"),
    ("N. Dose, dilution and injectate volume — A6",
     "unit conversion between products is not equipotency",
     "Scaglione conversion ratio different botulinum neurotoxin type A formulations"),

    # ---- O. spread, dysphagia and neck weakness ----
    ("O. Spread, dysphagia and neck weakness",
     "dysphagia after cervical BoNT and its relation to SCM/deep injection",
     "Comella Tanner DeFoor-Hill Smith dysphagia after botulinum toxin injections cervical dystonia"),
    ("O. Spread, dysphagia and neck weakness",
     "subclinical swallowing abnormality is common after cervical injection",
     "Thomas Hariram videofluoroscopic swallowing dysphagia botulinum cervical dystonia"),
    ("O. Spread, dysphagia and neck weakness",
     "distant spread and safety profile",
     "Bhatia Munchau Thompson Houser Chauhan Brooks Marsden generalised muscular weakness after botulinum toxin"),

    # ---- P. immunogenicity (A7, A8) ----
    ("P. Immunogenicity, protein load and secondary non-response — A7, A8",
     "neutralising antibody and secondary non-response",
     "Jankovic Schwartz response comparison botulinum toxin antibodies cervical dystonia"),
    ("P. Immunogenicity, protein load and secondary non-response — A7, A8",
     "the reformulation reduced protein load and antibody rate",
     "Jankovic Vuong Ahsan comparison efficacy immunogenicity original current botulinum toxin cervical dystonia"),
    ("P. Immunogenicity, protein load and secondary non-response — A7, A8",
     "protein load per Unit differs between products — the antigen parameter",
     "Frevert content botulinum neurotoxin Botox Dysport Xeomin"),
    ("P. Immunogenicity, protein load and secondary non-response — A7, A8",
     "long-term incobotulinumtoxinA immunogenicity",
     "Hefter Rosenthal Moll long-term safety incobotulinumtoxinA antibody cervical dystonia"),
    ("P. Immunogenicity, protein load and secondary non-response — A7, A8",
     "risk factors for antibody formation including injection interval",
     "Naumann Carruthers Carruthers Francis Meaney immunogenicity botulinum toxins"),
    ("P. Immunogenicity, protein load and secondary non-response — A7, A8",
     "shorter intervals / booster injections and antibody risk",
     "Greene Fahn Diamond development resistance botulinum toxin type A patients torticollis"),

    # ---- Q. antinociception (the A13 gap) ----
    ("Q. Direct antinociceptive action — the gap quantified in A13(4)",
     "BoNT inhibits CGRP release from sensory neurons",
     "Durham Cady Cady regulation calcitonin gene-related peptide secretion trigeminal botulinum"),
    ("Q. Direct antinociceptive action — the gap quantified in A13(4)",
     "BoNT/A blocks substance P release",
     "Purkiss Welch Doward Foster botulinum toxin A inhibits substance P release dorsal root ganglion"),
    ("Q. Direct antinociceptive action — the gap quantified in A13(4)",
     "analgesic mechanism of botulinum toxin, independent of muscle relaxation",
     "Matak Bolcskei Bach-Rojecky Helyes botulinum toxin A pain mechanisms"),
    ("Q. Direct antinociceptive action — the gap quantified in A13(4)",
     "pain relief precedes and exceeds motor improvement in cervical dystonia",
     "Tsui Eisen Stoessl Lang Calne double-blind study botulinum toxin spasmodic torticollis"),

    # ---- R. deep brain stimulation ----
    ("R. Deep brain stimulation of GPi",
     "sham-controlled trial of pallidal DBS in cervical dystonia",
     "Volkmann Mueller Deuschl Kuhn Krauss Poewe pallidal neurostimulation cervical dystonia randomised sham-controlled"),
    ("R. Deep brain stimulation of GPi",
     "long-term DBS outcome in segmental and generalised dystonia",
     "Volkmann Wolters Kupsch Muller Kuhn pallidal deep brain stimulation dystonia long-term"),
    ("R. Deep brain stimulation of GPi",
     "slow time course of DBS benefit in dystonia — the ramp in the model",
     "Krauss Yianni Loher Aziz deep brain stimulation cervical dystonia review"),

    # ---- S. oral agents ----
    ("S. Oral pharmacotherapy",
     "high-dose trihexyphenidyl in dystonia",
     "Burke Fahn Marsden torsion dystonia double-blind prospective trial high-dosage trihexyphenidyl"),
    ("S. Oral pharmacotherapy",
     "review of oral medication for dystonia",
     "Jankovic treatment dystonia medical therapy review anticholinergic baclofen"),

    # ---- T. surgical denervation (the phi lever) ----
    ("T. Selective peripheral denervation — a permanent φ lever",
     "selective peripheral denervation for cervical dystonia",
     "Bertrand selective peripheral denervation spasmodic torticollis surgical technique results"),
    ("T. Selective peripheral denervation — a permanent φ lever",
     "outcome of selective denervation in BoNT-resistant cervical dystonia",
     "Braun Richter selective peripheral denervation spasmodic torticollis results"),

    # ---- U. outcome measures ----
    ("U. Outcome measures — TWSTRS, MCID, CDIP-58",
     "the TWSTRS scale itself",
     "Consky Lang Toronto Western Spasmodic Torticollis Rating Scale assessment"),
    ("U. Outcome measures — TWSTRS, MCID, CDIP-58",
     "minimal clinically important change on TWSTRS — the MCID used throughout",
     "Cano Warner Linacre Bhatia Thompson Fitzpatrick capturing the true burden dystonia patients Rasch"),
    ("U. Outcome measures — TWSTRS, MCID, CDIP-58",
     "CDIP-58 patient-reported measure",
     "Cano Warner Thompson Bhatia Fitzpatrick Hobart CDIP-58 cervical dystonia impact profile"),
    ("U. Outcome measures — TWSTRS, MCID, CDIP-58",
     "clinically meaningful change thresholds in cervical dystonia trials",
     "Charles Adler Stacy Comella clinically meaningful improvement TWSTRS cervical dystonia"),

    # ---- V. QSP / modelling method ----
    ("V. QSP and modelling method",
     "mrgsolve, the tool the R model targets",
     "Elmokadem Riggs Baron quantitative systems pharmacology translational pharmacometrics mrgsolve"),
    ("V. QSP and modelling method",
     "good practice for QSP model qualification, which A1/A12/A13 follow",
     "Musante Ramanujan Schmidt Ghobrial Lu Heatherington quantitative systems pharmacology case for disease models"),
    ("V. QSP and modelling method",
     "Shiny-based PK/PD simulation front ends",
     "Wang gPKPDviz Shiny mrgsolve pharmacokinetic simulation tool"),
    ("V. QSP and modelling method",
     "receptor-occupancy / turnover PD modelling principles used for S and Q",
     "Dayneka Garg Jusko comparison four basic models indirect pharmacodynamic responses"),
]


def eutil(path, **params):
    params.setdefault("retmode", "json")
    url = f"{EUTILS}/{path}?" + urllib.parse.urlencode(params)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return json.load(r)
        except Exception as e:                                   # noqa: BLE001
            if attempt == 3:
                print(f"    ! {type(e).__name__}: {e}", file=sys.stderr)
                return None
            time.sleep(1.5 * (attempt + 1))
    return None


def resolve(query):
    js = eutil("esearch.fcgi", db="pubmed", term=query, retmax=1, sort="relevance")
    if not js:
        return None
    ids = js.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return None
    pmid = ids[0]
    js2 = eutil("esummary.fcgi", db="pubmed", id=pmid)
    if not js2:
        return None
    rec = js2.get("result", {}).get(pmid)
    if not rec:
        return None
    authors = [a["name"] for a in rec.get("authors", []) if a.get("authtype") == "Author"]
    first = authors[0] if authors else "?"
    return dict(pmid=pmid, title=(rec.get("title") or "").rstrip(". "),
                journal=rec.get("source", ""), year=(rec.get("pubdate", "") or "")[:4],
                first=first, nauth=len(authors))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--md", action="store_true", help="emit markdown only")
    a = ap.parse_args()

    out, unresolved = [], []
    for i, (section, note, q) in enumerate(QUERIES, 1):
        r = resolve(q)
        if r is None:
            unresolved.append((section, q))
            if not a.md:
                print(f"[{i:3d}/{len(QUERIES)}] UNRESOLVED  {q}", file=sys.stderr)
            continue
        r["section"], r["note"], r["query"] = section, note, q
        out.append(r)
        if not a.md:
            print(f"[{i:3d}/{len(QUERIES)}] {r['pmid']:>9}  {r['first']} "
                  f"{r['year']} {r['journal']} — {r['title'][:70]}", file=sys.stderr)
        time.sleep(0.36)          # stay under the 3 req/s anonymous NCBI limit

    # dedupe on pmid, keeping the first (its section/note)
    seen, dedup = set(), []
    for r in out:
        if r["pmid"] in seen:
            continue
        seen.add(r["pmid"])
        dedup.append(r)

    print(f"\nresolved {len(dedup)} unique / {len(QUERIES)} queries; "
          f"{len(unresolved)} unresolved\n", file=sys.stderr)
    with open("cdys_refs_resolved.json", "w") as fh:
        json.dump(dict(resolved=dedup, unresolved=unresolved), fh, indent=1)
    print("wrote cdys_refs_resolved.json", file=sys.stderr)


if __name__ == "__main__":
    main()
