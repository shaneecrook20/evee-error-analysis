"""
Chi-square tests for the protein structure results.

Reads structure_results/variants_with_structure.csv (made by structure_analysis.py)
and tests whether EVEE's miss rate (disease-causing variants scored below 0.5)
differs between groups:
  - buried vs surface
  - in vs outside a UniProt functional site
  - small vs big swaps (in the core, and on the surface)
  - buried vs surface (for small swaps, and for big swaps)
It also tests EVEE's false-alarm rate (harmless variants scored 0.5 or above)
for buried vs surface.

Uses scipy.stats.chi2_contingency, which applies Yates' correction to 2 x 2 tables.

Output: structure_results/chi_square_tests.txt
"""

import csv
import os
from scipy.stats import chi2_contingency

IN_FILE = os.path.join("structure_results", "variants_with_structure.csv")
OUT_FILE = os.path.join("structure_results", "chi_square_tests.txt")
CUTOFF = 0.5


def location(row):
    context = row["Structural context"]
    if context.startswith("1."):
        return "buried"
    if context.startswith("2."):
        return "surface"
    if context.startswith("3."):
        return "disordered"
    return None


def site(row):
    if row["UniProt site"].strip() != "":
        return "in site"
    return "not in site"


def size(row):
    if row["Size change class"].startswith("big"):
        return "big"
    return "small"


def location_and_size(row):
    place = location(row)
    if place == "buried" or place == "surface":
        return place + ", " + size(row)
    return None


def count(rows, group_func, clinvar_label):
    """Return {group: [wrong, right]} for one ClinVar label."""
    counts = {}
    for row in rows:
        if row["ClinVar"] != clinvar_label:
            continue
        group = group_func(row)
        if group is None:
            continue
        if group not in counts:
            counts[group] = [0, 0]
        score = float(row["EVEE score"])
        if clinvar_label == "pathogenic":
            wrong = score < CUTOFF
        else:
            wrong = score >= CUTOFF
        if wrong:
            counts[group][0] += 1
        else:
            counts[group][1] += 1
    return counts


def test(counts, group_a, group_b):
    table = [counts[group_a], counts[group_b]]
    chi2, p, dof, expected = chi2_contingency(table)
    rate_a = 100 * counts[group_a][0] / sum(counts[group_a])
    rate_b = 100 * counts[group_b][0] / sum(counts[group_b])
    return rate_a, rate_b, chi2, p


rows = []
with open(IN_FILE, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        rows.append(row)
print(f"Loaded {len(rows):,} variants")

tests = []
missed_by_location = count(rows, location, "pathogenic")
tests.append(("Missed: buried vs surface", missed_by_location, "buried", "surface"))
false_alarm_by_location = count(rows, location, "benign")
tests.append(("False alarms: buried vs surface", false_alarm_by_location, "buried", "surface"))
missed_by_site = count(rows, site, "pathogenic")
tests.append(("Missed: in vs outside a functional site", missed_by_site, "in site", "not in site"))
missed_by_combo = count(rows, location_and_size, "pathogenic")
tests.append(("Missed: buried small vs buried big", missed_by_combo, "buried, small", "buried, big"))
tests.append(("Missed: surface small vs surface big", missed_by_combo, "surface, small", "surface, big"))
tests.append(("Missed: buried small vs surface small", missed_by_combo, "buried, small", "surface, small"))
tests.append(("Missed: buried big vs surface big", missed_by_combo, "buried, big", "surface, big"))

lines = []
lines.append(f"{'Comparison':<42}{'Group A %':>11}{'Group B %':>11}{'Chi-square':>12}{'p-value':>12}")
for name, counts, a, b in tests:
    rate_a, rate_b, chi2, p = test(counts, a, b)
    lines.append(f"{name:<42}{rate_a:>10.1f}%{rate_b:>10.1f}%{chi2:>12.1f}{p:>12.2e}")
lines.append("")
lines.append("Group A is the first group named, group B the second. % = share EVEE got wrong.")

for line in lines:
    print(line)
with open(OUT_FILE, "w") as f:
    f.write("\n".join(lines) + "\n")
print(f"Saved {OUT_FILE}")
