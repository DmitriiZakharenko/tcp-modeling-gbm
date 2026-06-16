# TCP Modeling in Glioblastoma — Project Plan

## Repository

**Name:** `tcp-modeling-gbm`

**Description:**
> Tumor Control Probability modeling in glioblastoma using the CFB-GBM open radiotherapy dataset (TCIA). Reproducible Python workflow for NIfTI-based dosimetric feature extraction and dose-response modeling.

---

## Repository Structure

```
tcp-modeling-gbm/
├── data/
│   ├── raw/          # NIfTI files (RTDOSE, GTV) — excluded from git via .gitignore
│   └── processed/    # TSV/CSV feature tables and cohort files — versioned
├── notebooks/
│   ├── 01_cohort_overview.ipynb
│   ├── 02_feature_extraction.ipynb
│   ├── 03_tcp_models.ipynb
│   ├── 04_parameter_estimation.ipynb
│   └── 05_survival_analysis.ipynb
├── src/
│   ├── config.py         # Central paths (pathlib.Path)
│   ├── data/             # Cohort builder, NIfTI loader, DVH calculator, feature builder
│   ├── models/           # Base class, Poisson, Logistic, Probit, EUD-based TCP, bootstrap CI
│   └── utils/            # DVH and TCP plotting functions
├── reports/              # Written report sections and literature table
├── figures/              # Exported publication-quality figures (300 dpi PNG + PDF)
├── requirements.txt
├── README.md
├── .gitignore
└── TASKS.md
```

---

## Dataset — CFB-GBM (TCIA)

| Property | Value |
|---|---|
| Source | The Cancer Imaging Archive (TCIA) |
| DOI | https://doi.org/10.7937/V9PN-2F72 |
| Total patients | 264 |
| Patients with RTDOSE | 194 |
| Patients with GTV segmentation | 191 |
| Patients with unknown dose | ~70 (26.5%) — excluded from TCP modeling |
| Image format | NIfTI (.nii.gz), pre-processed (skull-stripped, co-registered) |
| RTDOSE resolution | 2.5 × 2.5 × 5 mm |
| Clinical data | TSV files (~50 KB total) — versioned in `data/processed/` |
| Full archive size | 208 GB (all MRI sequences + CT + RTDOSE) |

**Download strategy:** Only RTDOSE and GTV NIfTI files are required (~1–1.5 GB). Download via IBM Aspera Connect with folder-level filtering. Full 208 GB archive is not needed and should not be downloaded. HPC is not required.

**Fractionation:** Two schemes are present — 60 Gy / 30 fractions (standard Stupp protocol) and 40.05 Gy / 15 fractions (elderly / poor performance status). EQD2 correction (α/β = 10 Gy) is applied to pool both cohorts.

---

## Four-Week Timeline

### Week 1 — Infrastructure and Data Acquisition

**Goals:** set up repository, download clinical TSVs, build cohort table, start literature review.

- Initialize repository structure and `requirements.txt`
- Write `src/config.py` with central paths
- Download clinical and treatment TSV files from TCIA (available immediately, no Aspera needed)
- Write `src/data/cohort_builder.py`: merge TSVs, apply inclusion/exclusion criteria, export `cohort.csv`
- Run `notebooks/01_cohort_overview.ipynb`: cohort size, demographics, dose scheme breakdown
- Begin literature search on PubMed (TCP modeling in GBM, ≥8 papers)
- Write report section 1: Introduction and Clinical Background

**Output:** `data/processed/cohort.csv`, cohort summary notebook, literature table started

### Week 2 — Dosimetric Feature Extraction

**Goals:** parse NIfTI files, compute DVH metrics, build feature table.

- Write `src/data/nifti_loader.py`: load RTDOSE and GTV NIfTI files via `nibabel`
- Write `src/data/dvh_calculator.py`: compute cumulative DVH; extract D95, D98, D50, D2, Dmean, Dmax, GTV volume
- Write `src/data/feature_builder.py`: iterate over cohort, extract features, export `features.csv`
- Add EQD2 column: `EQD2 = D_total * (d_fraction + 10) / (2 + 10)`
- Run `notebooks/02_feature_extraction.ipynb`: descriptive statistics, DVH overlay plots
- Write report sections 2–4: Dataset Description, Curation Methodology, Feature Extraction

**Output:** `data/processed/features.csv`, DVH plots, report sections 2–4

### Week 3 — TCP Model Implementation and Statistical Evaluation

**Goals:** implement all TCP models, fit parameters, compute confidence intervals, run survival analysis.

- Write `src/models/base_model.py`: abstract `TCPModel` class with `fit`, `predict`, `log_likelihood`, `summary`
- Write `src/models/poisson_tcp.py`: Poisson TCP, MLE via `scipy.optimize.minimize`
- Write `src/models/logistic_tcp.py`: logistic dose-response, MLE
- Write `src/models/probit_tcp.py`: probit model, MLE
- Write `src/models/eud_tcp.py`: gEUD-based TCP with volume effect parameter `a`
- Write `src/models/bootstrap_ci.py`: bootstrap resampling (n=1000) for any `TCPModel` subclass
- Write `src/models/model_comparison.py`: AIC, BIC, log-likelihood, ROC AUC, calibration plots
- Write `src/models/survival_analysis.py`: Kaplan-Meier (OS), Cox PH regression with DVH covariates
- Run notebooks 03–05
- Write report sections 5–6: Mathematical Description of TCP Models, Parameter Estimation Methodology

**Output:** fitted models, parameter tables with 95% CI, model comparison table, survival curves

### Week 4 — Literature Comparison, Report, and Presentation

**Goals:** complete literature review, finalize report and repository, prepare oral presentation.

- Finalize literature table; compare estimated D50 / γ50 with published values
- Write report sections 9–11: Literature Review, Critical Discussion, Conclusion
- Assemble full report (10–15 pages); format references (Vancouver style)
- Export all publication-quality figures (300 dpi PNG + PDF)
- Final reproducibility check: all notebooks run top-to-bottom on clean environment
- Pin package versions in `requirements.txt`
- Build oral presentation (15 min, ~15 slides)

**Output:** final report PDF, figures, presentation, clean repository

---

## TCP Models

| Model | Parameters | Notes |
|---|---|---|
| Poisson TCP | D50, γ50 | Mechanistic; clonogen-based |
| Logistic | D50, k | Empirical; common baseline |
| Probit | μ, σ | Equivalent to logistic in practice |
| EUD-based TCP | D50, γ50, a | Uses gEUD from full DVH; volume effect |

All models fitted by MLE. Parameters D50 and γ50 are reported for all models (reparameterized where needed). 95% CI via bootstrap resampling (n=1000, seed=42).

Model comparison metrics: AIC, BIC, log-likelihood, ROC AUC, calibration plot.

---

## Team Roles

| Role | Responsibilities |
|---|---|
| Developer(s) | All coding tasks: data pipeline, models, notebooks, figures |
| Literature Reviewer | PubMed search, literature table, sections 9–10 (Literature, Discussion) |
| Report Writer | Sections 1–4, 7–8, 11 (Introduction, Methods, Results, Conclusion), assembly, presentation |

Detailed task breakdown with status tracking: see `TASKS.md`.

---

## Coding Standards

- All code in **English**: variable names, function names, docstrings, comments, commit messages
- Every function must have a **NumPy-style docstring**: Parameters, Returns, Raises
- No hardcoded paths — use `src/config.py` with `pathlib.Path`
- Notebooks must run **top-to-bottom without errors** (`Kernel → Restart & Run All`)
- Git commits: imperative mood, lowercase, specific — e.g. `add eud tcp model with gEUD from DVH`
- One feature per branch; merge via pull request

---

## Target Journals

| Journal | IF (2024) | Quartile | Format |
|---|---|---|---|
| Physics in Medicine & Biology | 3.4 | Q1 | Technical Note |
| Medical Physics | 3.5 | Q1 | Technical Note |
| Radiotherapy and Oncology | 5.3 | Q1 | Short Communication |

Publication requires going beyond course requirements: EUD-based TCP model, survival analysis with Cox regression, and honest uncertainty quantification are the key differentiators.

---

## References

Okunieff, P., Morgan, D., Niemierko, A., & Suit, H. D. (1995). Radiation dose-response of human tumors. *IJROBP*, 32(4), 1227–1237.

Moreau, N. N., et al. (2025). Pre and post treatment MRI and radiotherapy plans of patients with glioblastoma: the CFB-GBM cohort. *The Cancer Imaging Archive*. https://doi.org/10.7937/V9PN-2F72

Jarrett, D., Stride, E., Vallis, K., & Gooding, M. J. (2019). Applications and limitations of machine learning in radiation oncology. *British Journal of Radiology*, 92(1100), 20190001.
