"""
Part 2: Does protein structure explain which missense variants EVEE misses?

For ClinVar missense variants (1+ star, not conflicting) in genes with enough disease variants:
  1. UniProt  -> is the residue in a binding site, active site or other functional site?
  2. AlphaFold -> is the residue buried in the protein core or on the surface?
                 (relative solvent accessibility), and is the region well-folded (pLDDT)?
Then compare EVEE's miss rate and false-alarm rate across those groups, with AlphaMissense alongside.

Hypothesis: EVEE catches buried changes that break folding, but misses surface and
binding-site changes that alter what the protein does.

RUN THIS IN PYCHARM ON YOUR MAC (it needs the internet).
First time:  pip install duckdb biopython requests
For a quick test, set MAX_GENES = 30; leave it as None for the full run.
Downloads are cached in structure_cache/, so reruns are fast.

Outputs (structure_results/):
  variants_with_structure.csv   one row per variant with its structural context
  by_structure.csv, by_site.csv, by_plddt.csv, by_structure_and_size.csv
  report.txt                    everything printed below
"""

import duckdb
import requests
import csv
import json
import os
import re
import time
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

# ---------------- SETTINGS ----------------
FOLDER = os.environ.get("EVEE_DIR", os.path.expanduser("~/Downloads/EVEE"))
SHARDS = FOLDER + "/clean_shard_*.parquet"
CACHE = "structure_cache"
OUT = "structure_results"
CUTOFF = 0.5               # EVEE / AlphaMissense benign vs pathogenic line
MIN_STARS = 1
MIN_PATHOGENIC = 10        # only genes with at least this many pathogenic missense
MAX_GENES = None           # e.g. 30 for a quick test run; None = all genes (30-60 minutes the first time)
BURIED_RSA = 0.20          # relative solvent accessibility below this = buried
CONFIDENT_PLDDT = 70       # AlphaFold confidence; below this = likely disordered / unreliable
HIGH_MISS_FILE = "followup_results/genes_high_miss_rate.txt"   # the 80 genes, if present

os.makedirs(CACHE, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

THREE_TO_ONE = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
                "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
                "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
PDB_THREE = {"ALA": "Ala", "ARG": "Arg", "ASN": "Asn", "ASP": "Asp", "CYS": "Cys", "GLN": "Gln",
             "GLU": "Glu", "GLY": "Gly", "HIS": "His", "ILE": "Ile", "LEU": "Leu", "LYS": "Lys",
             "MET": "Met", "PHE": "Phe", "PRO": "Pro", "SER": "Ser", "THR": "Thr", "TRP": "Trp",
             "TYR": "Tyr", "VAL": "Val"}
# Maximum solvent-accessible area per residue, theoretical (Tien et al. 2013), square angstroms
MAX_ASA = {"Ala": 129, "Arg": 274, "Asn": 195, "Asp": 193, "Cys": 167, "Gln": 225, "Glu": 223,
           "Gly": 104, "His": 224, "Ile": 197, "Leu": 201, "Lys": 236, "Met": 224, "Phe": 240,
           "Pro": 159, "Ser": 155, "Thr": 172, "Trp": 285, "Tyr": 263, "Val": 174}
# Side-chain volume (Zamyatnin 1972), for the size-change breakdown
VOLUME = {"Gly": 60.1, "Ala": 88.6, "Ser": 89.0, "Cys": 108.5, "Asp": 111.1, "Pro": 112.7,
          "Asn": 114.1, "Thr": 116.1, "Glu": 138.4, "Val": 140.0, "Gln": 143.8, "His": 153.2,
          "Met": 162.9, "Ile": 166.7, "Leu": 166.7, "Lys": 168.6, "Arg": 173.4, "Phe": 189.9,
          "Tyr": 193.6, "Trp": 227.8}

report = []


def say(line=""):
    print(line)
    report.append(line)


# ======================================================================
# Step 1: missense variants from the EVEE dataset
# ======================================================================
def load_variants():
    con = duckdb.connect()
    con.execute("SET threads = 4")
    rows = con.execute(f"""
        SELECT gene_name, hgvsp_short, rs_id, chrom, pos, pathogenicity, gt_alphamissense_c,
               CASE
                   WHEN significance ILIKE '%pathogenic%' AND significance NOT ILIKE '%benign%'
                    AND significance NOT ILIKE '%conflicting%' THEN 'pathogenic'
                   WHEN significance ILIKE '%benign%' AND significance NOT ILIKE '%pathogenic%'
                    AND significance NOT ILIKE '%conflicting%' THEN 'benign'
                   ELSE 'other'
               END AS clin
        FROM read_parquet('{SHARDS}')
        WHERE consequence_display = 'Missense' AND stars >= {MIN_STARS}
          AND pathogenicity IS NOT NULL AND gene_name IS NOT NULL
    """).fetchall()
    pattern = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")
    variants = []
    for gene, prot, rs, chrom, gpos, evee, am, clin in rows:
        if clin == "other" or prot is None:
            continue
        m = pattern.match(prot)
        if not m:
            continue
        ref, pos, alt = m.group(1), int(m.group(2)), m.group(3)
        if ref not in THREE_TO_ONE or alt not in THREE_TO_ONE:
            continue
        variants.append({"gene": gene, "prot": prot, "ref": ref, "pos": pos, "alt": alt,
                         "evee": evee, "am": am, "clin": clin,
                         "rs_id": rs, "chrom": chrom, "genomic_pos": gpos})
    return variants


# ======================================================================
# Step 2: UniProt (accession, sequence, functional sites)
# ======================================================================
SITE_TYPES = {"Binding site": "binding site", "Active site": "active site",
              "Site": "other functional site", "Disulfide bond": "disulfide bond"}


def get_json(url, cache_file):
    path = os.path.join(CACHE, cache_file)
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 404:
                data = None
                break
            r.raise_for_status()
            data = r.json()
            break
        except requests.RequestException:
            time.sleep(2 * (attempt + 1))
            data = "error"
    if data == "error":
        return None   # don't cache network failures
    with open(path, "w") as f:
        json.dump(data, f)
    time.sleep(0.2)
    return data


def uniprot_for_gene(gene):
    url = ("https://rest.uniprot.org/uniprotkb/search?query=gene_exact:" + gene +
           "+AND+organism_id:9606+AND+reviewed:true"
           "&fields=accession,gene_primary,sequence,ft_binding,ft_act_site,ft_site,ft_disulfid"
           "&format=json&size=5")
    data = get_json(url, f"uniprot_{gene}.json")
    if not data or not data.get("results"):
        return None
    hit = None
    for entry in data["results"]:
        names = []
        for g in entry.get("genes", []):
            if "geneName" in g:
                names.append(g["geneName"]["value"])
        if gene in names:
            hit = entry
            break
    if hit is None:
        return None
    sites = {}   # residue number -> set of site labels
    for feat in hit.get("features", []):
        label = SITE_TYPES.get(feat.get("type"))
        if label is None:
            continue
        start = feat["location"]["start"].get("value")
        end = feat["location"]["end"].get("value")
        if start is None or end is None:
            continue
        if label == "disulfide bond":
            residues = [start, end]          # the two bonded cysteines
        else:
            residues = list(range(start, end + 1))
        for res in residues:
            if res not in sites:
                sites[res] = set()
            sites[res].add(label)
    return {"accession": hit["primaryAccession"], "sequence": hit["sequence"]["value"], "sites": sites}


# ======================================================================
# Step 3: AlphaFold (burial and confidence per residue)
# ======================================================================
def alphafold_residues(accession):
    """Returns {residue number: (three-letter aa, relative solvent accessibility, pLDDT)}."""
    out_path = os.path.join(CACHE, f"af_{accession}_residues.json")
    if os.path.exists(out_path):
        with open(out_path) as f:
            raw = json.load(f)
        result = {}
        for k in raw:
            result[int(k)] = tuple(raw[k])
        return result

    meta = get_json(f"https://alphafold.ebi.ac.uk/api/prediction/{accession}", f"af_meta_{accession}.json")
    if not meta:
        return None
    pdb_url = None
    for entry in meta:
        if entry.get("uniprotAccession") == accession and entry.get("pdbUrl"):
            pdb_url = entry["pdbUrl"]
            break
    if pdb_url is None:
        return None

    pdb_path = os.path.join(CACHE, f"af_{accession}.pdb")
    if not os.path.exists(pdb_path):
        try:
            r = requests.get(pdb_url, timeout=60)
            r.raise_for_status()
        except requests.RequestException:
            return None
        with open(pdb_path, "w") as f:
            f.write(r.text)

    result = residues_from_pdb(pdb_path)
    with open(out_path, "w") as f:
        json.dump({str(k): list(v) for k, v in result.items()}, f)
    os.remove(pdb_path)   # keep the cache small; the per-residue summary is saved
    return result


def residues_from_pdb(pdb_path):
    structure = PDBParser(QUIET=True).get_structure("af", pdb_path)
    model = structure[0]
    ShrakeRupley().compute(model, level="R")
    result = {}
    for chain in model:
        for res in chain:
            name = PDB_THREE.get(res.get_resname())
            if name is None:
                continue
            num = res.id[1]
            rsa = res.sasa / MAX_ASA[name]
            plddt = None
            if "CA" in res:
                plddt = res["CA"].get_bfactor()   # AlphaFold stores pLDDT in the B-factor column
            result[num] = (name, round(min(rsa, 1.5), 3), plddt)
    return result


# ======================================================================
# Step 4: classify each variant
# ======================================================================
def structural_context(rsa, plddt):
    if plddt is None or rsa is None:
        return None
    if plddt < CONFIDENT_PLDDT:
        return "3. Disordered / low confidence"
    if rsa < BURIED_RSA:
        return "1. Buried (protein core)"
    return "2. Surface"


def plddt_class(plddt):
    if plddt is None:
        return None
    if plddt >= 90:
        return "1. Very confident (90+)"
    if plddt >= 70:
        return "2. Confident (70-90)"
    if plddt >= 50:
        return "3. Low (50-70)"
    return "4. Very low (under 50)"


def size_class(ref, alt):
    dv = abs(VOLUME[alt] - VOLUME[ref])
    if dv < 30:
        return "small size change (<30)"
    return "big size change (30+)"


# ======================================================================
# Step 5: tables
# ======================================================================
def table(variants, title, key_func, filename):
    groups = {}
    for v in variants:
        k = key_func(v)
        if k is None:
            continue
        if k not in groups:
            groups[k] = [0, 0, 0, 0, 0, 0, 0, 0]
        g = groups[k]
        if v["clin"] == "pathogenic":
            g[0] += 1
            if v["evee"] < CUTOFF:
                g[1] += 1
            if v["am"] is not None:
                g[2] += 1
                if v["am"] < CUTOFF:
                    g[3] += 1
        else:
            g[4] += 1
            if v["evee"] >= CUTOFF:
                g[5] += 1
            if v["am"] is not None:
                g[6] += 1
                if v["am"] >= CUTOFF:
                    g[7] += 1

    def pct(a, b):
        return "" if b == 0 else round(100 * a / b, 1)

    say(f"=== {title} ===")
    say(f"{'Group':<52}{'Path n':>8}{'EVEE miss%':>11}{'AM miss%':>10}{'Benign n':>10}{'EVEE FP%':>10}{'AM FP%':>8}")
    rows = []
    for k in sorted(groups):
        g = groups[k]
        row = [k, g[0], pct(g[1], g[0]), pct(g[3], g[2]), g[4], pct(g[5], g[4]), pct(g[7], g[6])]
        rows.append(row)
        say(f"{row[0]:<52}{row[1]:>8,}{str(row[2]):>11}{str(row[3]):>10}{row[4]:>10,}{str(row[5]):>10}{str(row[6]):>8}")
    say()
    with open(os.path.join(OUT, filename), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Group", "Pathogenic n", "EVEE miss %", "AlphaMissense miss %",
                    "Benign n", "EVEE false-alarm %", "AlphaMissense false-alarm %"])
        for row in rows:
            w.writerow(row)


def main():
    print("Loading missense variants from the EVEE dataset...")
    variants = load_variants()
    print(f"  {len(variants):,} missense variants")

    path_count = {}
    for v in variants:
        if v["clin"] == "pathogenic":
            path_count[v["gene"]] = path_count.get(v["gene"], 0) + 1
    genes = []
    for g, n in sorted(path_count.items(), key=lambda p: -p[1]):
        if n >= MIN_PATHOGENIC:
            genes.append(g)

    high_miss = set()
    if os.path.exists(HIGH_MISS_FILE):
        with open(HIGH_MISS_FILE) as f:
            for line in f:
                if line.strip():
                    high_miss.add(line.strip())
    if MAX_GENES is not None:
        # test run: half high-miss genes, half other genes, so there is something to compare
        hm = []
        other = []
        for g in genes:
            if g in high_miss:
                hm.append(g)
            else:
                other.append(g)
        n_hm = min(len(hm), MAX_GENES // 2)
        genes = hm[:n_hm] + other[:MAX_GENES - n_hm]
    print(f"  Analysing {len(genes)} genes")

    by_gene = {}
    for v in variants:
        by_gene.setdefault(v["gene"], []).append(v)

    kept = []
    stats = {"no_uniprot": 0, "no_alphafold": 0, "aa_mismatch": 0, "no_residue": 0}
    for i, gene in enumerate(genes):
        print(f"  [{i + 1}/{len(genes)}] {gene}", end="", flush=True)
        up = uniprot_for_gene(gene)
        if up is None:
            stats["no_uniprot"] += 1
            print(" - no UniProt entry")
            continue
        af = alphafold_residues(up["accession"])
        if af is None:
            stats["no_alphafold"] += 1
            print(f" - no AlphaFold model ({up['accession']})")
            continue
        n_kept = 0
        for v in by_gene[gene]:
            pos = v["pos"]
            seq = up["sequence"]
            # Check that EVEE's protein numbering matches UniProt's sequence
            if pos > len(seq) or seq[pos - 1] != THREE_TO_ONE[v["ref"]]:
                stats["aa_mismatch"] += 1
                continue
            if pos not in af or af[pos][0] != v["ref"]:
                stats["no_residue"] += 1
                continue
            name, rsa, plddt = af[pos]
            v["accession"] = up["accession"]
            v["rsa"] = rsa
            v["plddt"] = plddt
            v["context"] = structural_context(rsa, plddt)
            site = up["sites"].get(pos)
            v["site"] = ", ".join(sorted(site)) if site else ""
            v["in_site"] = "1. In a functional site (UniProt)" if site else "2. Not in an annotated site"
            v["high_miss_gene"] = "yes" if gene in high_miss else ""
            kept.append(v)
            n_kept += 1
        print(f" - {up['accession']}, {n_kept} variants")

    say(f"Genes analysed: {len(genes)}")
    say(f"Skipped genes: {stats['no_uniprot']} without a UniProt entry, {stats['no_alphafold']} without a single AlphaFold model "
        "(very long proteins are split into fragments and skipped)")
    say(f"Skipped variants: {stats['aa_mismatch']} where EVEE's amino acid didn't match UniProt's sequence, "
        f"{stats['no_residue']} not in the AlphaFold model")
    say(f"Variants with structure: {len(kept):,}")
    say()

    table(kept, "Buried vs surface (AlphaFold)", lambda v: v["context"], "by_structure.csv")
    table(kept, "Functional sites (UniProt: binding, active, other sites, disulfide bonds)",
          lambda v: v["in_site"], "by_site.csv")
    table(kept, "AlphaFold confidence (pLDDT)", lambda v: plddt_class(v["plddt"]), "by_plddt.csv")

    def combo(v):
        if v["context"] is None or v["context"].startswith("3."):
            return None
        place = "Buried" if v["context"].startswith("1.") else "Surface"
        return f"{place}, {size_class(v['ref'], v['alt'])}"
    table(kept, "Location x size of the swap (well-folded regions only)", combo, "by_structure_and_size.csv")

    if high_miss:
        table([v for v in kept if v["high_miss_gene"] == "yes"], "Buried vs surface: the 80 high-miss genes only",
              lambda v: v["context"], "by_structure_high_miss_genes.csv")

    with open(os.path.join(OUT, "variants_with_structure.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Gene", "UniProt", "Protein change", "rs ID", "Chrom", "Position",
                    "ClinVar", "EVEE score", "AlphaMissense score",
                    "Relative solvent accessibility", "pLDDT", "Structural context", "UniProt site",
                    "Size change (cubic angstroms)", "Size change class", "High-miss gene"])
        for v in kept:
            size_change = round(abs(VOLUME[v["alt"]] - VOLUME[v["ref"]]), 1)
            w.writerow([v["gene"], v["accession"], v["prot"], v["rs_id"], v["chrom"], v["genomic_pos"],
                        v["clin"], v["evee"], v["am"], v["rsa"],
                        v["plddt"], v["context"], v["site"],
                        size_change, size_class(v["ref"], v["alt"]), v["high_miss_gene"]])

    say("Notes: AlphaFold models are single chains, so residues at interfaces with other proteins count as 'surface'.")
    say(f"Buried = relative solvent accessibility under {BURIED_RSA}; disordered = pLDDT under {CONFIDENT_PLDDT}.")
    with open(os.path.join(OUT, "report.txt"), "w") as f:
        f.write("\n".join(report) + "\n")
    print(f"\nSaved to {OUT}/")


if __name__ == "__main__":
    main()
