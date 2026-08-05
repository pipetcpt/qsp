#!/usr/bin/env python3
"""
fetch_refs.py -- build the reference list for the myopia QSP model from PubMed
records that were ACTUALLY RETURNED by NCBI E-utilities.

Nothing is written from memory: every entry below is a SEARCH, and only the
record NCBI hands back (its real PMID, title, journal, year, authors) is
written into the reference file.  Searches that return nothing are reported as
misses rather than filled in.
"""
import json, time, urllib.request, urllib.parse, sys, re

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def _get(url):
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            if attempt == 3:
                print("   ! %s" % e, file=sys.stderr)
                return None
            time.sleep(1.5 * (attempt + 1))


def esearch(term, n=1):
    u = EUTILS + "esearch.fcgi?db=pubmed&retmode=json&sort=relevance&retmax=%d&term=%s" \
        % (n, urllib.parse.quote(term))
    d = _get(u)
    if not d:
        return []
    return d.get("esearchresult", {}).get("idlist", [])


def esummary(pmids):
    if not pmids:
        return {}
    u = EUTILS + "esummary.fcgi?db=pubmed&retmode=json&id=" + ",".join(pmids)
    d = _get(u)
    if not d:
        return {}
    return d.get("result", {})


# (section, search term, must-contain word in the title or "" for no check)
QUERIES = [
 # ---------------- 1. epidemiology and burden ----------------
 ("epi", "Holden global prevalence myopia high myopia 2000 through 2050 Ophthalmology", "prevalence"),
 ("epi", "Morgan Ohno-Matsui Saw myopia Lancet seminar", "myopia"),
 ("epi", "Baird myopia Nature Reviews Disease Primers", "myopia"),
 ("epi", "Dong prevalence myopia children meta-analysis worldwide", "myopia"),
 ("epi", "Rudnicka global variations myopia prevalence children ethnic differences meta-analysis", "myopia"),
 ("epi", "Vitale increased prevalence myopia United States 1971 2004 Arch Ophthalmol", "myopia"),
 ("epi", "Sankaridurg IMI impact myopia", "myopia"),
 ("epi", "Wang progression myopia home confinement COVID-19 school children JAMA Ophthalmology", "myopia"),
 ("epi", "Tricard progression myopia children COVID-19 lockdown", "myopia"),
 ("epi", "Jiang progression myopia Chinese children 2019 2021 photoscreening", "myopia"),
 # ---------------- 2. natural history, biometry, refraction ----------------
 ("hist", "Mutti Zadnik axial growth changes crystalline lens power emmetropization CLEERE", "lens"),
 ("hist", "Mutti refractive error axial length components onset juvenile myopia", "myopia"),
 ("hist", "Iribarren crystalline lens power refractive error review", "lens"),
 ("hist", "Zadnik ocular component development children Optom Vis Sci", "ocular"),
 ("hist", "Rozema Ocular biometry emmetropization review", "emmetropization"),
 ("hist", "Truckenbrod Longitudinal analysis axial length growth German children", "axial"),
 ("hist", "Sanz Diez growth curves axial length children percentile", "axial"),
 ("hist", "Hou Axial length percentile curves myopia children reference", "axial"),
 ("hist", "Chamberlain Non-adjusted axial length growth emmetropic children", "axial"),
 ("hist", "Gwiazda Thorn Bauer Held emmetropization accommodation myopic children", "myopi"),
 ("hist", "Atchison eye shape myopia MRI", "shape"),
 ("hist", "Verkicharla Ocular shape peripheral refraction myopia", "refraction"),
 ("hist", "Chakraborty diurnal variations axial length choroidal thickness human eye", "diurnal"),
 ("hist", "Stone Nickla diurnal rhythms eye growth ocular", "rhythm"),
 # ---------------- 3. emmetropization, defocus and visual signals ----------------
 ("sig", "Wallman Winawer Homeostasis eye growth question myopia Neuron", "growth"),
 ("sig", "Troilo IMI report emmetropization animal models experimental myopia", "myopia"),
 ("sig", "Smith Hung Huang peripheral vision refractive development rhesus monkeys", "refractive"),
 ("sig", "Smith Ramamirtham effects foveal ablation emmetropization monkeys", "emmetropization"),
 ("sig", "Zhu Winawer Wallman potency myopic defocus brief periods", "defocus"),
 ("sig", "Schaeffel Feldkaemper animal models myopia research review", "myopia"),
 ("sig", "Norton Siegwart Light levels refractive development tree shrew", "refractive"),
 ("sig", "Fischer McGuire Stell Light Ohlsson ZENK glucagon amacrine ocular growth", "ZENK"),
 ("sig", "Feldkaemper Schaeffel glucagon signalling eye growth chick", "glucagon"),
 ("sig", "Vessey Rushforth Stell Glucagon receptor agonists antagonists form deprivation chick", "glucagon"),
 ("sig", "Rohrer Stell Basic fibroblast growth factor transforming growth factor beta chick eye growth", "growth"),
 ("sig", "Rucker Wallman Chick eyes compensate chromatic simulations hyperopic myopic defocus", "defocus"),
 ("sig", "Ashby Schaeffel effect bright light form deprivation myopia chickens", "light"),
 ("sig", "Karouta Ashby Correlation light levels myopia protection", "light"),
 ("sig", "Read Collins Vincent Light exposure eye growth childhood", "light"),
 ("sig", "Landis Myopia protective effect objectively measured light exposure children", "light"),
 ("sig", "Norton Siegwart Light levels myopia tree shrew ambient", "myopia"),
 ("sig", "Flitcroft contribution optical defocus dioptric environment myopia", "defocus"),
 ("sig", "Charman Radhakrishnan peripheral refraction myopia development review", "peripheral"),
 ("sig", "Mutti Sinnott Peripheral refraction ocular shape children CLEERE", "refraction"),
 ("sig", "Atchison Peripheral refraction along horizontal vertical visual fields myopia", "refraction"),
 # ---------------- 4. dopamine, adenosine and retinal neurochemistry ----------------
 ("dopa", "Feldkaemper Schaeffel dopamine mechanism myopia review", "dopamine"),
 ("dopa", "Zhou Pardue Iuvone Qu Dopamine signaling myopia development what are the key challenges", "Dopamine"),
 ("dopa", "Stone Lin Laties Iuvone Retinal dopamine form-deprivation myopia chicks", "dopamine"),
 ("dopa", "Iuvone Tigges Stone Lambert Laties Effects apomorphine dopamine agonist ocular refraction", "dopamine"),
 ("dopa", "Ashby Schaeffel dopamine antagonist spiperone bright light protective myopia", "dopamine"),
 ("dopa", "Nickla Wildsoet Wallman nitric oxide choroidal thickness compensation defocus", "nitric"),
 ("dopa", "Nickla Wilken Lytle Yom Wildsoet Nitric oxide synthase inhibitors choroidal recovery", "nitric"),
 ("dopa", "Cui Trier Munk Nitric oxide myopia review", "nitric"),
 ("dopa", "Trier Munk Ribel-Madsen Cui Systemic 7-methylxanthine myopia progression children", "methylxanthine"),
 ("dopa", "Trier Cui Ribel-Madsen Guggenheim Oral administration caffeine 7-methylxanthine sclera", "methylxanthine"),
 ("dopa", "Cui Trier Zeng effects 7-methylxanthine posterior sclera guinea pig myopia", "methylxanthine"),
 ("dopa", "Stone Pendrak Sugimoto Adenosine receptor antagonists myopia", "myopia"),
 ("dopa", "Chakraborty Ostrin Nickla Iuvone Pardue circadian rhythms refractive development myopia", "circadian"),
 ("dopa", "Stone McGlinn Baldwin Melatonin ocular growth", "melatonin"),
 # ---------------- 5. choroid ----------------
 ("chor", "Nickla Wallman multifunctional choroid review", "choroid"),
 ("chor", "Wildsoet Wallman Choroidal scleral mechanisms compensation ocular defocus chick", "choroidal"),
 ("chor", "Read Collins Vincent Alonso-Caneiro Choroidal thickness childhood", "Choroidal"),
 ("chor", "Read Alonso-Caneiro Vincent Collins Longitudinal changes choroidal thickness children myopia", "choroidal"),
 ("chor", "Jiang Zhu Choroidal thickness myopia control children meta-analysis", "choroid"),
 ("chor", "Wang Zhou Choroidal thickness axial length refractive error meta-analysis", "choroidal"),
 ("chor", "Ye Zhao Choroidal thickness change low-concentration atropine children myopia", "atropine"),
 ("chor", "Zhang Wildsoet Choroid retinal pigment epithelium eye growth signal cascade", "choroid"),
 ("chor", "Summers choroid as sclera-controlling tissue review", "choroid"),
 ("chor", "Nickla choroidal thickness compensation defocus role dopamine", "choroid"),
 ("chor", "Hoseini-Yazdi Read Collins Diurnal variations choroidal thickness young adults", "choroidal"),
 ("chor", "Ostrin Choroidal thickness ocular biometry short-term light", "choroidal"),
 # ---------------- 6. sclera: matrix and biomechanics ----------------
 ("scl", "Rada Shelton Norton sclera myopia review extracellular matrix", "sclera"),
 ("scl", "McBrien Gentle Role sclera development pathological myopia review", "sclera"),
 ("scl", "Norton Siegwart Light Increased scleral creep tree shrew induced myopia", "creep"),
 ("scl", "Siegwart Norton Regulation scleral extracellular matrix remodeling tree shrew", "scleral"),
 ("scl", "Guggenheim McBrien Form-deprivation myopia gelatinase MMP-2 tree shrew sclera", "MMP"),
 ("scl", "Rada Perry Slover Achromatic scleral matrix metalloproteinase TIMP myopia", "matrix"),
 ("scl", "Jobling Nguyen McBrien Gentle Auditory transforming growth factor beta scleral remodelling myopia", "growth"),
 ("scl", "McBrien Lawlor Gentle Scleral remodeling during development recovery myopia", "Scleral"),
 ("scl", "Rada Nickla Troilo Decreased proteoglycan synthesis chick sclera myopic defocus", "proteoglycan"),
 ("scl", "Wu Chen Zhou Scleral hypoxia HIF-1alpha myopia development PNAS", "myopia"),
 ("scl", "Zhao Wu Zhou Scleral HIF-1alpha myopia progression", "scleral"),
 ("scl", "Backhouse Gentle Scleral remodelling myopia recovery review", "scleral"),
 ("scl", "Metlapally Wildsoet Scleral mechanisms underlying ocular growth myopia", "Scleral"),
 ("scl", "Wang Li Zhang scleral crosslinking genipin form-deprivation myopia guinea pig", "crosslink"),
 ("scl", "Dotan Kremer Livnat Scleral collagen cross-linking riboflavin ultraviolet", "cross-linking"),
 ("scl", "Levy Roberts Boote Structure sclera myopia collagen fibril", "sclera"),
 ("scl", "Grytz Siegwart Changing material properties tree shrew sclera during minus lens", "sclera"),
 # ---------------- 7. atropine and muscarinic pharmacology ----------------
 ("atro", "Yam Jiang Tang Low-Concentration Atropine Myopia Progression LAMP Study randomized", "Atropine"),
 ("atro", "Yam Li Tang Low-Concentration Atropine Myopia Progression LAMP Phase 2 report", "Atropine"),
 ("atro", "Yam Zhang Tang Three-year clinical trial Low-Concentration Atropine Myopia Progression LAMP washout", "Atropine"),
 ("atro", "Chua Balakrishnan Chan Atropine treatment myopia ATOM Ophthalmology 2006", "Atropine"),
 ("atro", "Chia Chua Cheung Atropine treatment myopia ATOM2 0.5% 0.1% 0.01% doses", "Atropine"),
 ("atro", "Chia Lu Tan Five-year clinical trial atropine myopia 2 myopia control atropine 0.01 eye drops", "Atropine"),
 ("atro", "Zadnik Schulman Flitcroft Efficacy safety 0.01% 0.02% atropine CHAMP randomized clinical trial", "Atropine"),
 ("atro", "Loughman Flitcroft Myopia Outcome Study Atropine Children MOSAIC randomised", "Atropine"),
 ("atro", "Lee Lingham Yazar Atropine 0.01% myopia control Western Australia", "Atropine"),
 ("atro", "Repka Weise Chandler Low-Dose 0.01% Atropine Eye Drops vs Placebo Myopia Progression PEDIG", "Atropine"),
 ("atro", "Hansen Atropine 0.1% 0.01% myopia control Danish children", "atropine"),
 ("atro", "Ha Kim Shim Comparison atropine concentrations myopia progression network meta-analysis", "atropine"),
 ("atro", "Lanca Repka Grzybowski Topical review atropine myopia control efficacy", "atropine"),
 ("atro", "Tan Atropine myopia mechanism muscarinic receptor review", "atropine"),
 ("atro", "McBrien Stell Carr How does atropine exert antimyopia effects", "atropine"),
 ("atro", "Barathi Beuerman Molecular mechanisms muscarinic receptors mouse myopia", "muscarinic"),
 ("atro", "Arumugam McBrien Muscarinic antagonist control experimental myopia selective M1", "myopia"),
 ("atro", "Siatkowski Cotter Miller Pirenzepine ophthalmic gel myopia progression children", "Pirenzepine"),
 ("atro", "Tan Lam Chua Pirenzepine ophthalmic gel myopia Asian children", "Pirenzepine"),
 ("atro", "Lind Chew Marzani Stell Muscarinic acetylcholine receptor antagonists chick scleral chondrocytes", "muscarinic"),
 ("atro", "Kaymak Fricke Mirshahi rebound effect atropine cessation myopia", "atropine"),
 ("atro", "Cui Yao Zhou rebound myopia progression after cessation atropine meta-analysis", "atropine"),
 ("atro", "Wu Chen Atropine pupil size accommodation amplitude adverse effects low concentration", "atropine"),
 ("atro", "Chiang Wang Systemic absorption plasma concentration topical atropine children", "atropine"),
 ("atro", "Cooper Eisenberg Schulman Maximum atropine dose without clinical signs symptoms", "atropine"),
 # ---------------- 8. optical interventions ----------------
 ("opt", "Lam Tang Tse Defocus Incorporated Multiple Segments DIMS spectacle lenses myopia randomised", "Defocus"),
 ("opt", "Lam Tang Qi Myopia control DIMS spectacle lenses 6-year results", "DIMS"),
 ("opt", "Chamberlain Peixoto-de-Matos Logan 3-year randomized clinical trial MiSight dual-focus contact lens", "dual-focus"),
 ("opt", "Chamberlain Bradley Arumugam Long-term effect dual-focus contact lenses myopia 6 year", "dual-focus"),
 ("opt", "Bao Huang Li Spectacle lenses highly aspherical lenslets myopia control 2-year", "lenslets"),
 ("opt", "Bao Yang Huang One-year myopia control efficacy spectacle lenses highly aspherical lenslets", "lenslets"),
 ("opt", "Cho Cheung Retardation Myopia Orthokeratology ROMIO randomized clinical trial", "Orthokeratology"),
 ("opt", "Cho Cheung Edwards longitudinal orthokeratology research children LORIC", "orthokeratology"),
 ("opt", "Sun Xu Zhang Orthokeratology myopia control meta-analysis axial length", "Orthokeratology"),
 ("opt", "Walline Greiner McVey Jones-Jordan Multifocal contact lens myopia control", "Multifocal"),
 ("opt", "Walline Walker Mutti BLINK Effect High Add Power Multifocal Contact Lenses Myopia Progression", "Multifocal"),
 ("opt", "Berntsen Sinnott Mutti Zadnik randomized trial progressive addition lenses accommodative lag myopia", "lenses"),
 ("opt", "Gwiazda Hyman Hussein COMET randomized clinical trial progressive addition lenses myopia", "progressive"),
 ("opt", "Chung Mohidin O'Leary Undercorrection myopia enhances rather than inhibits myopia progression", "Undercorrection"),
 ("opt", "Adler Millodot possible effect undercorrection myopic progression children", "undercorrection"),
 ("opt", "Kinoshita Konno Hikoya Additive effects orthokeratology low-dose atropine axial elongation", "orthokeratology"),
 ("opt", "Tan Zhou Combined atropine orthokeratology myopia control randomized", "orthokeratology"),
 ("opt", "Sankaridurg IMI Clinical Management Guidelines myopia control", "myopia"),
 ("opt", "Lawrenson Shah Huntjens Interventions myopia progression children Cochrane network meta-analysis", "myopia"),
 ("opt", "Prousali Efficacy interventions myopia control children systematic review network meta-analysis", "myopia"),
 # ---------------- 9. red light, environment interventions ----------------
 ("rlrl", "Jiang Zhu Tan Effect Repeated Low-level Red-light Therapy Myopia Control randomized", "Red-light"),
 ("rlrl", "Dong Zhi Xu Repeated low-level red-light therapy myopia control high myopia children", "red-light"),
 ("rlrl", "Xiong Zhu Jiang rebound effect after cessation repeated low-level red-light therapy", "red-light"),
 ("rlrl", "He Wang Xie Effect Repeated Low-Level Red Light myopia incidence premyopia", "Red Light"),
 ("rlrl", "Ostrin Schill Red light instruments myopia safety retinal", "light"),
 ("rlrl", "He Xiang Zeng Outdoor Activity Myopia Among Rural Chinese Schoolchildren randomized", "Myopia"),
 ("rlrl", "Wu Chen Tsai Outdoor activity myopia onset progression school children", "Outdoor"),
 ("rlrl", "Wu Chen Yang Myopia Prevention Outdoor Light Intervention Taiwan", "Outdoor"),
 ("rlrl", "Jin Hua Effect outdoor activity myopia children meta-analysis", "outdoor"),
 ("rlrl", "Xiong Sankaridurg Naduvilath Time spent outdoors myopia incidence progression meta-analysis", "outdoors"),
 # ---------------- 10. genetics ----------------
 ("gen", "Hysi Choquet Khawaja meta-analysis refractive error myopia genome-wide", "refractive"),
 ("gen", "Tedja Wojciechowski Hysi Genome-wide association meta-analysis refractive error CREAM", "refractive"),
 ("gen", "Tedja Haarman Meester-Smoor IMI Myopia Genetics Report", "Genetics"),
 ("gen", "Verhoeven Hysi Saw Genome-wide meta-analyses multiancestry refractive error myopia", "myopia"),
 ("gen", "Mountjoy Davies Plotnikov Education myopia Mendelian randomization", "myopia"),
 ("gen", "Guggenheim Northstone Williams Time outdoors physical activity myopia Mendelian", "myopia"),
 ("gen", "Ghorbani Mojarrad Plotnikov Williams Association polygenic risk score myopia refractive error", "polygenic"),
 ("gen", "Jiang Ohno-Matsui Wang ARR3 mutation early-onset high myopia female", "myopia"),
 # ---------------- 11. pathologic myopia and complications ----------------
 ("path", "Tideman Snabel Tedja Association Axial Length Risk Irreversible Visual Impairment", "Axial"),
 ("path", "Haarman Enthoven Tideman Complications Myopia Systematic Review Meta-Analysis", "Myopia"),
 ("path", "Ohno-Matsui Kawasaki Jonas International photographic classification grading system myopic maculopathy", "maculopathy"),
 ("path", "Ohno-Matsui Wu Yamashiro IMI Pathologic Myopia", "Myopia"),
 ("path", "Fang Yokoi Nagaoka Progression Myopic Maculopathy During 18-Year Follow-up", "Maculopathy"),
 ("path", "Hayashi Ohno-Matsui Shimada Long-term pattern progression myopic maculopathy", "myopic"),
 ("path", "Flitcroft Complex interactions myopia risk factors ocular disease", "myopia"),
 ("path", "Bullimore Brennan Myopia Control Why Each Diopter Matters", "Diopter"),
 ("path", "Bullimore Ritchey Shah Lopes Overestimation risk myopia control", "myopia"),
 ("path", "Marcus de Vries Junoy Montolio Jansonius Myopia glaucoma meta-analysis", "glaucoma"),
 ("path", "Ohno-Matsui Lai Cheung Updates pathologic myopia review", "myopia"),
 ("path", "Wong Ferreira Hughes Carter Mitchell Epidemiology disease burden pathologic myopia", "myopia"),
 ("path", "Verkicharla Ohno-Matsui Saw Current epidemiology myopia prevalence pathologic", "myopia"),
 # ---------------- 12. endpoints, methods, modelling ----------------
 ("meth", "Brennan Toubouti Cheng Bullimore Efficacy myopia control cumulative absolute reduction axial elongation", "myopia"),
 ("meth", "Chamberlain Lanca Bullimore axial length targets myopia control", "axial"),
 ("meth", "Wolffsohn Kollbaum Berntsen IMI Clinical Measurements myopia report", "myopia"),
 ("meth", "Wolffsohn Whayeb Logan Lingham IMI Industry Guidelines Ethical Considerations myopia control", "Myopia"),
 ("meth", "Hughes Jones Read Cumulative absolute reduction axial elongation limitations", "axial"),
 ("meth", "Kaphle Varnas Sankaridurg Predicting refraction axial length myopia model", "myopia"),
 ("meth", "Guo Tan Wang mathematical model emmetropization myopia development simulation", "myopia"),
 ("meth", "Rozema Dankert Iribarren Axial growth crystalline lens model refractive development", "model"),
 ("meth", "Nti Berntsen Optical changes myopia control lenses peripheral defocus measurement", "defocus"),
 ("meth", "Queiros Villa-Collar Gonzalez-Meijome peripheral refraction orthokeratology myopic defocus", "orthokeratology"),
]


def main():
    seen, recs, misses = set(), [], []
    for i, (sec, term, must) in enumerate(QUERIES):
        ids = esearch(term, n=3)
        got = None
        if ids:
            res = esummary(ids)
            for pid in ids:
                r = res.get(pid)
                if not r or "title" not in r:
                    continue
                if must and must.lower() not in r["title"].lower():
                    continue
                if pid in seen:
                    continue
                got = (pid, r)
                break
            # relax the title filter if nothing matched
            if got is None:
                for pid in ids:
                    r = res.get(pid)
                    if r and "title" in r and pid not in seen:
                        got = (pid, r)
                        break
        if got is None:
            misses.append((sec, term))
            print("MISS  %-6s %s" % (sec, term[:64]), file=sys.stderr)
        else:
            pid, r = got
            seen.add(pid)
            au = r.get("authors", [])
            first = au[0]["name"] if au else "?"
            recs.append(dict(sec=sec, pmid=pid, title=r["title"].rstrip("."),
                             journal=r.get("source", ""),
                             year=(r.get("pubdate", "") or "")[:4],
                             first=first, nauth=len(au)))
        time.sleep(0.36)
        if (i + 1) % 25 == 0:
            print("  ... %d/%d" % (i + 1, len(QUERIES)), file=sys.stderr)
    json.dump(dict(recs=recs, misses=misses), open("refs.json", "w"), indent=1)
    print("\nfetched %d unique records, %d misses" % (len(recs), len(misses)))


if __name__ == "__main__":
    main()
