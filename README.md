# EVEE error analysis

Code for comparing EVEE (Evo Variant Effect Explorer) pathogenicity predictions with ClinVar labels, to find which kinds of variants EVEE gets right and wrong. These scripts produce the results in Table 1 of the paper.

## Data

The EVEE dataset is not included because of its size (five Parquet files, about 35 GB). Download it from Zenodo:

Pearce et al. (2026). *EVEE: Interpretable variant effect prediction from genomic foundation model embeddings* (v1). https://doi.org/10.5281/zenodo.19701997

Put the five files (`clean_shard_0.parquet` to `clean_shard_4.parquet`) in one folder. The scripts look in `~/Downloads/EVEE` by default. To use a different folder, set the `EVEE_DIR` environment variable:

```
export EVEE_DIR=/path/to/EVEE
```

## Setup

Python 3 is required (3.14 was used).

```
pip install -r requirements.txt
```

## Scripts

Run them in this order from the repository folder.

| Order | Script | What it does | Output |
|---|---|---|---|
| 1 | `EVEE vs Clinvar Errors.py` | Labels ClinVar variants as disease-causing or harmless (1+ star, conflicting excluded), applies the EVEE 0.5 cutoff, and counts right and wrong by mutation type | `error_rates_by_type.csv`, `false_negatives.csv`, `false_positives.csv` |
| 2 | `amino_acid_properties.py` | For missense variants, calculates the Grantham distance of each amino acid swap and counts right and wrong for EVEE and AlphaMissense in each Grantham class | `aa_properties_results/` |
| 3 | `structure_analysis.py` | Downloads UniProt entries and AlphaFold structures, checks each variant's amino acid matches UniProt, calculates relative solvent accessibility, and labels each variant buried, surface or disordered, small or big swap, and in or outside a UniProt functional site. Needs an internet connection; downloads are cached in `structure_cache/` | `structure_results/` |
| 4 | `chi_square_tests.py` | Chi-square tests comparing groups from step 3 | `structure_results/chi_square_tests.txt` |

## Settings

Each script has its settings near the top:

- `CUTOFF = 0.5`: EVEE score at or above which a variant is called disease-causing
- `MIN_STARS = 1`: minimum ClinVar review stars

## Notes

- The ClinVar labels are those included in EVEE dataset v1.
- Buried = relative solvent accessibility under 0.2; disordered = AlphaFold pLDDT under 70; big swap = side-chain volume change of 30 cubic angstroms or more.
