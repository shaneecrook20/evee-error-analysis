"""
EVEE vs ClinVar: list of false negatives and false positives.

  False negative = ClinVar pathogenic, EVEE says benign   (EVEE score < CUTOFF)
  False positive = ClinVar benign,     EVEE says pathogenic (EVEE score >= CUTOFF)

"Confident" errors use stricter cutoffs (EVEE very sure and still wrong).

Outputs:
  false_negatives.csv
  false_positives.csv
"""

import duckdb
import csv
import os

# ---------------- SETTINGS ----------------
FOLDER = os.environ.get("EVEE_DIR", os.path.expanduser("~/Downloads/EVEE"))
ALL_SHARDS = FOLDER + "/clean_shard_*.parquet"
CUTOFF = 0.5            # EVEE's benign / pathogenic line
CONFIDENT_FN = 0.1      # pathogenic variant scored below this = confident miss
CONFIDENT_FP = 0.9      # benign variant scored above this = confident false alarm
MIN_STARS = 1           # ClinVar confidence (2 = stricter)
GENES = []              # e.g. ["APOE", "APOL1", "APOB"]; empty = all genes
SHOW_VARIANTS = 30      # how many to print per list (CSVs have all)
TOP_GENES = 15

con = duckdb.connect()
con.execute("SET threads = 4")
con.execute("SET preserve_insertion_order = false")

gene_filter = ""
if len(GENES) > 0:
    gene_list = ", ".join("'" + g + "'" for g in GENES)
    gene_filter = f"AND gene_name IN ({gene_list})"

# ---------- Step 1: pull confident ClinVar variants ----------
print("Searching all shards (may take a minute)...")
rows = con.execute(f"""
    SELECT gene_name, hgvsp_short, rs_id, chrom, pos, consequence_display,
           significance, stars, pathogenicity, gt_alphamissense_c,
           CASE
               WHEN significance ILIKE '%pathogenic%'
                AND significance NOT ILIKE '%benign%'
                AND significance NOT ILIKE '%conflicting%' THEN 'pathogenic'
               WHEN significance ILIKE '%benign%'
                AND significance NOT ILIKE '%pathogenic%'
                AND significance NOT ILIKE '%conflicting%' THEN 'benign'
               ELSE 'other'
           END AS clin_class
    FROM read_parquet('{ALL_SHARDS}')
    WHERE stars >= {MIN_STARS}
      AND pathogenicity IS NOT NULL
      {gene_filter}
""").fetchall()

# ---------- Step 2: sort into right / wrong ----------
n_path = 0
n_benign = 0
false_negatives = []
false_positives = []

# Per mutation type: [pathogenic total, EVEE misses, benign total, EVEE false alarms,
#                     AM pathogenic scored, AM misses, AM benign scored, AM false alarms]
by_type = {}

for row in rows:
    gene, protein, rsid, chrom, pos, cons, sig, stars, evee, am, clin_class = row
    variant = {
        "gene": gene, "protein": protein, "rsid": rsid, "chrom": chrom,
        "pos": pos, "consequence": cons, "clinvar": sig, "stars": stars,
        "evee": evee, "alphamissense": am,
    }
    if clin_class == "other":
        continue
    if cons not in by_type:
        by_type[cons] = [0, 0, 0, 0, 0, 0, 0, 0]
    counts = by_type[cons]

    if clin_class == "pathogenic":
        n_path += 1
        counts[0] += 1
        if evee < CUTOFF:
            false_negatives.append(variant)
            counts[1] += 1
        if am is not None:
            counts[4] += 1
            if am < CUTOFF:
                counts[5] += 1
    elif clin_class == "benign":
        n_benign += 1
        counts[2] += 1
        if evee >= CUTOFF:
            false_positives.append(variant)
            counts[3] += 1
        if am is not None:
            counts[6] += 1
            if am >= CUTOFF:
                counts[7] += 1

# Most confident errors first
false_negatives.sort(key=lambda v: v["evee"])
false_positives.sort(key=lambda v: v["evee"], reverse=True)

confident_fn = 0
for v in false_negatives:
    if v["evee"] < CONFIDENT_FN:
        confident_fn += 1
confident_fp = 0
for v in false_positives:
    if v["evee"] > CONFIDENT_FP:
        confident_fp += 1

# ---------- Step 3: summary ----------
print("\n=== Summary ===")
print(f"ClinVar pathogenic: {n_path:,}   ClinVar benign: {n_benign:,}")
if n_path > 0:
    print(f"False negatives: {len(false_negatives):,} "
          f"({100 * len(false_negatives) / n_path:.1f}% of pathogenic), "
          f"confident (< {CONFIDENT_FN}): {confident_fn:,}")
if n_benign > 0:
    print(f"False positives: {len(false_positives):,} "
          f"({100 * len(false_positives) / n_benign:.1f}% of benign), "
          f"confident (> {CONFIDENT_FP}): {confident_fp:,}")


# ---------- Step 3b: error RATE per mutation type ----------
def pct(part, whole):
    if whole == 0:
        return "   -  "
    return f"{100 * part / whole:5.1f}%"


print("\n=== Error rate by mutation type (EVEE, with AlphaMissense where it has scores) ===")
print(f"{'Mutation type':<22}{'Path n':>8}{'EVEE FN%':>10}{'AM FN%':>9}"
      f"{'Benign n':>10}{'EVEE FP%':>10}{'AM FP%':>9}")
type_rows = list(by_type.items())
type_rows.sort(key=lambda p: p[1][0] + p[1][2], reverse=True)
for cons, c in type_rows:
    print(f"{str(cons):<22}{c[0]:>8,}{pct(c[1], c[0]):>10}{pct(c[5], c[4]):>9}"
          f"{c[2]:>10,}{pct(c[3], c[2]):>10}{pct(c[7], c[6]):>9}")
print("FN% = % of pathogenic EVEE called benign.  FP% = % of benign EVEE called pathogenic.")
print("AM only scores missense, so '-' elsewhere is expected.")

# Save the rate table too
with open("error_rates_by_type.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Mutation type", "Pathogenic n", "EVEE misses", "EVEE FN rate",
                     "AM pathogenic n", "AM misses", "Benign n", "EVEE false alarms",
                     "EVEE FP rate", "AM benign n", "AM false alarms"])
    for cons, c in type_rows:
        fn_rate = round(c[1] / c[0], 4) if c[0] > 0 else ""
        fp_rate = round(c[3] / c[2], 4) if c[2] > 0 else ""
        writer.writerow([cons, c[0], c[1], fn_rate, c[4], c[5],
                         c[2], c[3], fp_rate, c[6], c[7]])


# ---------- Helpers ----------
def count_by(variants, key):
    counts = {}
    for v in variants:
        value = v[key]
        if value in counts:
            counts[value] += 1
        else:
            counts[value] = 1
    pairs = list(counts.items())
    pairs.sort(key=lambda p: p[1], reverse=True)
    return pairs


def print_variants(variants, title):
    shown = min(SHOW_VARIANTS, len(variants))
    print(f"\n=== {title} (showing {shown} of {len(variants):,}) ===")
    print(f"{'Gene':<10}{'Protein change':<20}{'rsID':<14}{'Consequence':<24}"
          f"{'EVEE':<8}{'AlphaMiss':<11}ClinVar")
    count = 0
    for v in variants:
        if count >= SHOW_VARIANTS:
            break
        am_text = "NA" if v["alphamissense"] is None else format(v["alphamissense"], ".3f")
        print(f"{str(v['gene']):<10}{str(v['protein']):<20}{str(v['rsid']):<14}"
              f"{str(v['consequence']):<24}{v['evee']:<8.3f}{am_text:<11}"
              f"{v['clinvar']} ({v['stars']} star)")
        count += 1


def print_breakdowns(variants, title):
    print(f"\n--- {title}: mutation types ---")
    for cons, n in count_by(variants, "consequence"):
        print(f"  {cons}: {n:,}")
    print(f"\n--- {title}: top {TOP_GENES} genes ---")
    count = 0
    for gene, n in count_by(variants, "gene"):
        if count >= TOP_GENES:
            break
        names = []
        for v in variants:
            if v["gene"] == gene:
                names.append(str(v["protein"]) + " (" + format(v["evee"], ".2f") + ")")
        shown = names[:8]
        extra = f" ... +{len(names) - 8} more" if len(names) > 8 else ""
        print(f"  {gene} [{n}]: " + ", ".join(shown) + extra)
        count += 1


# ---------- Step 4: the lists ----------
print_variants(false_negatives, "FALSE NEGATIVES: ClinVar pathogenic, EVEE benign (lowest score first)")
print_breakdowns(false_negatives, "False negatives")

print_variants(false_positives, "FALSE POSITIVES: ClinVar benign, EVEE pathogenic (highest score first)")
print_breakdowns(false_positives, "False positives")


# ---------- Step 5: save ----------
def save(variants, filename):
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Gene", "Protein change", "rsID", "Chrom", "Pos",
                         "Consequence", "ClinVar", "Stars", "EVEE", "AlphaMissense"])
        for v in variants:
            writer.writerow([v["gene"], v["protein"], v["rsid"], v["chrom"], v["pos"],
                             v["consequence"], v["clinvar"], v["stars"],
                             v["evee"], v["alphamissense"]])
    print(f"Saved {len(variants):,} rows to {filename}")


print()
save(false_negatives, "false_negatives.csv")
save(false_positives, "false_positives.csv")