"""
Part 1: Does EVEE only catch dramatic amino acid changes?

For every ClinVar missense variant (1+ star, not conflicting), measure how different
the new amino acid is from the old one:
  - Grantham distance (overall chemical difference; Grantham 1974)
  - Size change (side-chain volume, cubic angstroms; Zamyatnin 1972)
  - Hydrophobicity change (Kyte-Doolittle 1982)
Then compare EVEE's miss rate (pathogenic called benign) and false-alarm rate
(benign called pathogenic) across those, with AlphaMissense alongside.

Outputs (aa_properties_results/):
  by_grantham.csv, by_size_change.csv, by_hydrophobicity_change.csv,
  swap_grid_miss_rate.csv, swap_grid_counts.csv, report.txt
"""

import duckdb
import csv
import os
import re

FOLDER = os.environ.get("EVEE_DIR", os.path.expanduser("~/Downloads/EVEE"))
SHARDS = FOLDER + "/clean_shard_*.parquet"
OUT = "aa_properties_results"
CUTOFF = 0.5
MIN_STARS = 1
os.makedirs(OUT, exist_ok=True)

# Grantham (1974): composition c, polarity p, volume v
GRANTHAM_CPV = {
    "Ser": (1.42, 9.2, 32), "Arg": (0.65, 10.5, 124), "Leu": (0.00, 4.9, 111),
    "Pro": (0.39, 8.0, 32.5), "Thr": (0.71, 8.6, 61), "Ala": (0.00, 8.1, 31),
    "Val": (0.00, 5.9, 84), "Gly": (0.74, 9.0, 3), "Ile": (0.00, 5.2, 111),
    "Phe": (0.00, 5.2, 132), "Tyr": (0.20, 6.2, 136), "Cys": (2.75, 5.5, 55),
    "His": (0.58, 10.4, 96), "Gln": (0.89, 10.5, 85), "Asn": (1.33, 11.6, 56),
    "Lys": (0.33, 11.3, 119), "Asp": (1.38, 13.0, 54), "Glu": (0.92, 12.3, 83),
    "Met": (0.00, 5.7, 105), "Trp": (0.13, 5.4, 170),
}
# Side-chain volume in cubic angstroms (Zamyatnin 1972)
VOLUME = {
    "Gly": 60.1, "Ala": 88.6, "Ser": 89.0, "Cys": 108.5, "Asp": 111.1, "Pro": 112.7,
    "Asn": 114.1, "Thr": 116.1, "Glu": 138.4, "Val": 140.0, "Gln": 143.8, "His": 153.2,
    "Met": 162.9, "Ile": 166.7, "Leu": 166.7, "Lys": 168.6, "Arg": 173.4, "Phe": 189.9,
    "Tyr": 193.6, "Trp": 227.8,
}
# Kyte-Doolittle hydropathy (higher = more hydrophobic / oily)
KD = {
    "Ile": 4.5, "Val": 4.2, "Leu": 3.8, "Phe": 2.8, "Cys": 2.5, "Met": 1.9, "Ala": 1.8,
    "Gly": -0.4, "Thr": -0.7, "Ser": -0.8, "Trp": -0.9, "Tyr": -1.3, "Pro": -1.6,
    "His": -3.2, "Glu": -3.5, "Gln": -3.5, "Asp": -3.5, "Asn": -3.5, "Lys": -3.9, "Arg": -4.5,
}
AAS = ["Gly", "Ala", "Ser", "Cys", "Thr", "Pro", "Asp", "Asn", "Glu", "Gln",
       "Val", "His", "Met", "Ile", "Leu", "Lys", "Arg", "Phe", "Tyr", "Trp"]


def grantham(a, b):
    c1, p1, v1 = GRANTHAM_CPV[a]
    c2, p2, v2 = GRANTHAM_CPV[b]
    d = (1.833 * (c1 - c2) ** 2 + 0.1018 * (p1 - p2) ** 2 + 0.000399 * (v1 - v2) ** 2) ** 0.5
    return d * 50.723   # scaling Grantham used so the average distance is 100


def grantham_class(d):
    # Li et al. (1984) classes
    if d <= 50:
        return "1. Conservative (0-50)"
    if d <= 100:
        return "2. Moderately conservative (51-100)"
    if d <= 150:
        return "3. Moderately radical (101-150)"
    return "4. Radical (over 150)"


def size_class(dv):
    dv = abs(dv)
    if dv < 20:
        return "1. Small (under 20)"
    if dv < 50:
        return "2. Medium (20-50)"
    if dv < 80:
        return "3. Large (50-80)"
    return "4. Very large (80+)"


def hydro_class(dh):
    dh = abs(dh)
    if dh < 1:
        return "1. Little change (under 1)"
    if dh < 3:
        return "2. Some change (1-3)"
    if dh < 5:
        return "3. Big change (3-5)"
    return "4. Very big change (5+)"


con = duckdb.connect()
con.execute("SET threads = 2")
print("Loading missense variants...")
rows = con.execute(f"""
    SELECT hgvsp_short, pathogenicity, gt_alphamissense_c,
           CASE
               WHEN significance ILIKE '%pathogenic%' AND significance NOT ILIKE '%benign%'
                AND significance NOT ILIKE '%conflicting%' THEN 'pathogenic'
               WHEN significance ILIKE '%benign%' AND significance NOT ILIKE '%pathogenic%'
                AND significance NOT ILIKE '%conflicting%' THEN 'benign'
               ELSE 'other'
           END AS clin
    FROM read_parquet('{SHARDS}')
    WHERE consequence_display = 'Missense' AND stars >= {MIN_STARS}
      AND pathogenicity IS NOT NULL
""").fetchall()

pattern = re.compile(r"^p\.([A-Z][a-z]{2})\d+([A-Z][a-z]{2})$")
variants = []
for prot, evee, am, clin in rows:
    if clin == "other" or prot is None:
        continue
    m = pattern.match(prot)
    if not m:
        continue
    ref = m.group(1)
    alt = m.group(2)
    if ref not in VOLUME or alt not in VOLUME or ref == alt:
        continue
    variants.append({
        "ref": ref, "alt": alt, "clin": clin, "evee": evee, "am": am,
        "grantham": grantham(ref, alt),
        "dvol": VOLUME[alt] - VOLUME[ref],
        "dkd": KD[alt] - KD[ref],
    })
print(f"Usable missense variants: {len(variants):,}")

report = []


def say(line=""):
    print(line)
    report.append(line)


def table(title, key_func, filename):
    groups = {}
    for v in variants:
        k = key_func(v)
        if k not in groups:
            groups[k] = {"pn": 0, "pmiss": 0, "pam": 0, "pam_miss": 0,
                         "bn": 0, "bfp": 0, "bam": 0, "bam_fp": 0}
        g = groups[k]
        if v["clin"] == "pathogenic":
            g["pn"] += 1
            if v["evee"] < CUTOFF:
                g["pmiss"] += 1
            if v["am"] is not None:
                g["pam"] += 1
                if v["am"] < CUTOFF:
                    g["pam_miss"] += 1
        else:
            g["bn"] += 1
            if v["evee"] >= CUTOFF:
                g["bfp"] += 1
            if v["am"] is not None:
                g["bam"] += 1
                if v["am"] >= CUTOFF:
                    g["bam_fp"] += 1

    def pct(a, b):
        return "" if b == 0 else round(100 * a / b, 1)

    say(f"=== {title} ===")
    say(f"{'Group':<38}{'Path n':>8}{'EVEE miss%':>11}{'AM miss%':>10}{'Benign n':>10}{'EVEE FP%':>10}{'AM FP%':>8}")
    out = []
    for k in sorted(groups):
        g = groups[k]
        row = [k, g["pn"], pct(g["pmiss"], g["pn"]), pct(g["pam_miss"], g["pam"]),
               g["bn"], pct(g["bfp"], g["bn"]), pct(g["bam_fp"], g["bam"])]
        out.append(row)
        say(f"{row[0]:<38}{row[1]:>8,}{str(row[2]):>11}{str(row[3]):>10}{row[4]:>10,}{str(row[5]):>10}{str(row[6]):>8}")
    say()
    with open(os.path.join(OUT, filename), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Group", "Pathogenic n", "EVEE miss %", "AlphaMissense miss %",
                    "Benign n", "EVEE false-alarm %", "AlphaMissense false-alarm %"])
        for row in out:
            w.writerow(row)


table("By Grantham distance (how chemically different the swap is)",
      lambda v: grantham_class(v["grantham"]), "by_grantham.csv")
table("By size change (side-chain volume, cubic angstroms)",
      lambda v: size_class(v["dvol"]), "by_size_change.csv")
table("By hydrophobicity change (Kyte-Doolittle)",
      lambda v: hydro_class(v["dkd"]), "by_hydrophobicity_change.csv")

# Average Grantham distance: missed vs caught pathogenic
missed = []
caught = []
for v in variants:
    if v["clin"] == "pathogenic":
        if v["evee"] < CUTOFF:
            missed.append(v["grantham"])
        else:
            caught.append(v["grantham"])
say("=== Average Grantham distance of disease-causing missense ===")
say(f"Missed by EVEE: {sum(missed) / len(missed):.1f}  (n = {len(missed):,})")
say(f"Caught by EVEE: {sum(caught) / len(caught):.1f}  (n = {len(caught):,})")
say()

# 20 x 20 swap grid: EVEE miss rate for pathogenic missense
grid_n = {}
grid_m = {}
for v in variants:
    if v["clin"] != "pathogenic":
        continue
    key = (v["ref"], v["alt"])
    grid_n[key] = grid_n.get(key, 0) + 1
    if v["evee"] < CUTOFF:
        grid_m[key] = grid_m.get(key, 0) + 1

with open(os.path.join(OUT, "swap_grid_miss_rate.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["From \\ To"] + AAS)
    for a in AAS:
        line = [a]
        for b in AAS:
            n = grid_n.get((a, b), 0)
            line.append("" if n < 20 else round(100 * grid_m.get((a, b), 0) / n, 1))
        w.writerow(line)
with open(os.path.join(OUT, "swap_grid_counts.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["From \\ To"] + AAS)
    for a in AAS:
        line = [a]
        for b in AAS:
            line.append(grid_n.get((a, b), 0))
        w.writerow(line)

# Most and least missed swaps (20+ pathogenic variants)
swaps = []
for key in grid_n:
    n = grid_n[key]
    if n >= 20:
        swaps.append((grid_m.get(key, 0) / n, n, key))
swaps.sort(reverse=True)
say("=== Swaps EVEE misses most (20+ disease variants) ===")
for rate, n, (a, b) in swaps[:10]:
    say(f"  {a}->{b}: {100 * rate:.0f}% missed (n = {n}, Grantham {grantham(a, b):.0f})")
say("=== Swaps EVEE misses least ===")
for rate, n, (a, b) in swaps[-10:]:
    say(f"  {a}->{b}: {100 * rate:.0f}% missed (n = {n}, Grantham {grantham(a, b):.0f})")
say("Swap grids: blank = fewer than 20 disease variants for that swap.")

with open(os.path.join(OUT, "report.txt"), "w") as f:
    f.write("\n".join(report) + "\n")
print(f"\nSaved to {OUT}/")
