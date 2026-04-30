# Neuropixels Population Correlation Architecture: Layer 2H Pairing-Family Decomposition

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Allen Brain Observatory](https://img.shields.io/badge/Data-DANDI%20000021-orange)](https://dandiarchive.org/dandiset/000021)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)

**Manuscript:** *Pairing-family decomposition of population correlation architecture during visual stimulation reveals signed reorganization dominance and continuous-discrete dissociation in mouse cortex Neuropixels recordings*

**Author:** Drake H. Harbert · [ORCID 0009-0007-7740-3616](https://orcid.org/0009-0007-7740-3616)
**Affiliation:** Inner Architecture LLC, Canton, OH 44721, USA
**Contact:** Drake@InnerArchitectureLLC.com

---

## Summary

This repository contains the complete analytical pipeline for a Layer 2H pairing-family architectural decomposition of population correlation reorganization in mouse visual cortex Neuropixels recordings. Five of six Layer 2H pairing types are empirically demonstrated on this neural substrate: Type 1 (continuous-discrete dissociation), Type 2 (sign-treatment dominance), Type 4 (cross-session consensus), Type 5 (local-global by brain area), and Type 6 (substrate-projection linearity validation). Type 3 (asymmetric-symmetric) is documented as inapplicable under percentile-threshold set construction.

**Methodology lineage:** Layer 2H framework introduced in Harbert (2026), Frontiers in Pharmacology, [doi:10.3389/fphar.2026.1830847](https://doi.org/10.3389/fphar.2026.1830847).

---

## Key Findings

| Layer 2H Type | Measurement | Result |
|---|---|---|
| Type 1 (Continuous-Discrete) | WJ_unsigned − BJ_top5 | gap = +0.179 ± 0.071 (n = 2,052; 30 sessions) |
| Type 2 (Sign-Treatment) | WJ_signed − WJ_unsigned | gap = 0.527; sign inversions = 91.92% |
| Type 3 (Asymmetric-Symmetric) | Tversky asymmetry on top-5% sets | Inapplicable (equal set cardinality by construction) |
| Type 4 (Set-vs-Multiset) | Cross-session sign-inversion consensus | 29/29 sessions ≥85%; 29/29 ≥90%; SD = 0.62% |
| Type 5 (Local-Global) | Per-area − full-population WJ | 25/29 areas show local exceeds global |
| Type 6 (Substrate-Projection) | WJ_spearman − WJ_pearson | gap ≈ −0.002; |gap|<0.05 in 100% of comparisons |

---

## Repository Structure

```
neuropixels-correlation-architecture-wj/
│
├── README.md                              (this file)
├── LICENSE                                (MIT)
├── requirements.txt
├── .gitignore
│
├── neuropixels_wj_pipeline.py             (Type 2 main pipeline; per-session WJ)
├── area_decomposition.py                  (per-brain-area decomposition)
├── area_wj_analysis.py                    (per-area aggregate analysis)
├── layer2h_analysis.py                    (Types 2/4/5/6 from existing summaries)
├── layer2h_type1_type3.py                 (Type 1 + Type 3 from raw NWB; 30-session)
│
├── sensitivity_bin_size.py
├── sensitivity_block_permutation.py
├── sensitivity_quality_filter.py
├── sensitivity_running_speed.py
├── sensitivity_sign_inversion.py
│
├── download_all_sessions.py               (download NWB files from DANDI)
├── generate_manuscript.py                 (manuscript builder)
│
└── results/
    ├── area_wj_summary.csv                (Type 5 per-area)
    ├── clean_summary.json
    ├── final_clean_summary.json           (Type 2 + Type 4 summary)
    ├── layer2h_type1_type3_per_comparison.csv  (n = 2,052 Type 1 comparisons)
    ├── layer2h_type1_type3_summary.csv
    ├── layer2h_summary_provenance.json
    ├── layer2h_type1_type3_provenance.json
    ├── sign_inversion_sensitivity.csv
    ├── quality_filter_effect.csv
    ├── bin_size_sensitivity.csv
    ├── running_speed_*.csv
    └── block_permutation_comparison.csv
```

---

## Reproduction

### Step 1: Download data

```bash
python download_all_sessions.py
```

Downloads all 30 NWB sessions from DANDI archive (accession 000021). Total ≈70 GB.

### Step 2: Run main WJ pipeline (Type 2)

```bash
python neuropixels_wj_pipeline.py
```

Computes per-session pairwise WJ for all condition pairs. Produces `results/ses-*_summary.json`. Runtime: ~30 min/session × 30 = ~15 hours.

### Step 3: Run per-area decomposition (Type 5)

```bash
python area_decomposition.py
python area_wj_analysis.py
```

Produces `results/area_wj_summary.csv`. Runtime: ~10 min.

### Step 4: Run Layer 2H Type 1 + Type 3 (NWB-direct compute)

```bash
python layer2h_type1_type3.py
```

Computes Type 1 (Continuous-Discrete) and Type 3 (Asymmetric-Symmetric) on all 30 sessions from raw NWB. Runtime: ~10.4 hours sequential.

### Step 5: Run Layer 2H summary analysis (Types 2/4/5/6)

```bash
python layer2h_analysis.py
```

Reads existing summary JSONs and produces consolidated Layer 2H output. Runtime: ~30 seconds.

### Step 6: Run sensitivity analyses

```bash
python sensitivity_quality_filter.py
python sensitivity_bin_size.py
python sensitivity_running_speed.py
python sensitivity_sign_inversion.py
python sensitivity_block_permutation.py
```

Reproducibility: all pipelines use `RANDOM_SEED = 42` and `FORCE_RECOMPUTE = True`.

---

## Citation

When the manuscript is published:

```
Harbert DH. Pairing-family decomposition of population correlation architecture
during visual stimulation reveals signed reorganization dominance and
continuous-discrete dissociation in mouse cortex Neuropixels recordings.
PLOS Computational Biology. (Submitted). [DOI when available]
```

Methodology reference (Layer 2H framework):

```
Harbert DH. Sigma-1 and Sigma-2 receptors exhibit divergent genome-wide
co-expression architectures in human brain despite shared subcellular
localization. Frontiers in Pharmacology. 2026.
https://doi.org/10.3389/fphar.2026.1830847
```

---

## Data Availability

| Source | Use |
|--------|-----|
| Allen Brain Observatory Visual Coding Neuropixels | DANDI archive accession 000021 |
| 30 sessions; 1,305-2,395 spike-sorted units per session | https://dandiarchive.org/dandiset/000021 |
| 6-8 brain regions (visual cortex, thalamus, hippocampus) | Public access |

---

## License

MIT License — see [LICENSE](LICENSE).

## Acknowledgments

Thanks to the Allen Institute for Brain Science for generating and publicly
releasing the Visual Coding Neuropixels dataset, and to the DANDI archive team
for maintaining open access to the data.
