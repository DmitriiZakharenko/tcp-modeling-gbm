# Task Board — TCP Modeling in GBM

## How to use this board

Checkboxes below are **clickable directly on GitHub** (no terminal needed):
- Open this file on [github.com/DmitriiZakharenko/tcp-modeling-gbm](https://github.com/DmitriiZakharenko/tcp-modeling-gbm/blob/main/TASKS.md)
- Click any `[ ]` to mark it done — GitHub saves automatically
- **Mirror issue:** [GitHub Issue #1](https://github.com/DmitriiZakharenko/tcp-modeling-gbm/issues/1) — sync via `bash scripts/sync_github_issue_1.sh` (requires `gh auth login`)

Roles: **[CODE]** · **[LIT]** literature & clinical review · **[WRITE]** report writing

**Last synced:** 2026-06-28 · **Git:** `cd967bf` · **Report:** `make report` → `reports/RESULTS.md`  
**Group guide:** [`reports/group_glossary_guide.md`](reports/group_glossary_guide.md) — plain-language terms and presentation FAQ

---

## Dataset Facts (CFB-GBM v3, TCIA) — verified 2026-06-28

| Metric | Value | Source |
|---|---|---|
| Total patients in TSV | 264 | `cohort.csv` |
| Has RTDOSE (imaging flag) | 194 | `cohort.csv` |
| Has GTV segmentation | 191 | `cohort.csv` |
| Both RTDOSE + GTV on disk | 191 | `verify_raw_data.py` |
| Unknown RT dose (excluded) | 70 | `cohort.csv` exclusion_reason |
| **Included after cohort rules** | **190** | `cohort.csv` (`included=True`) |
| **Modeling table (post DVH QC)** | **190 × 58** | `modeling_table.csv` (v3 RANO + clinical) |
| DVH QC exclusion | 1 (patient 32, Dmean=0 Gy) | `dvh_qc.py` |
| Primary fractionation | 120 × 60 Gy/30 fr; 61 × 40.05 Gy/15 fr | `modeling_table.csv` |
| RANO t0→t1 labels | **137 / 190** | `modeling_table.csv` |
| RANO in 40 Gy arm | **34** | `modeling_table.csv` |
| Median age (modeling) | 70 yr | `modeling_table.csv` |
| Sex M/F (modeling) | 117 / 73 | `modeling_table.csv` |
| Median OS (modeling) | 51 wk | `modeling_table.csv` |

### Key analysis results (auto-generated)

| Analysis | Result |
|---|---|
| Pooled TCP OS proxy (EQD2, n=190) | AUC ≈ **0.68** |
| Pooled EQD2 → RANO (n=137) | AUC ≈ **0.43** (confounded) |
| Pooled volume + clinical + scheme → RANO | AUC ≈ **0.72**; LOOCV ≈ **0.64** |
| 40 Gy volume → RANO (n=34) | in-sample AUC ≈ **0.90**; LOOCV ≈ **0.74** |
| PyRadiomics top-5 → RANO (t1gd GTV t0) | in-sample AUC ≈ **0.78**; nested CV ≈ **0.74** vs DVH volume **0.71** / **0.70** |
| Cox OS ~ RANO (n=137) | HR ≈ **0.48**, p ≈ 0.0009 |

- Format: **NIfTI (.nii.gz)**, pre-processed (skull-stripped, co-registered) — not DICOM
- Two fractionation schemes: 60 Gy / 30 fr (Stupp) and 40.05 Gy / 15 fr (elderly / poor PS)
- v3 adds: RANO, MRI/CT availability, PyRadiomics TSV, RT delay, BMI
- Clinical outcome: OS (weeks); RANO non-PD at t1 as secondary imaging endpoint

## Data Download Strategy

**Step 1 — Clinical TSVs (v3, ~200 KB, no Aspera):** `python -m src.data.download_clinical_data`

| File | Notes |
|---|---|
| `CFB-GBM_clinical_data_v03_*.tsv` | Demographics, OS |
| `CFB-GBM_rano_criteria_v03_*.tsv` | RANO t0/t1/t2 |
| `CFB_GBM_features_extraction_pyradiomics_v03_*.tsv` | Author radiomics |
| Treatment / imaging availability / dictionary | v2 + v3 merged in `cohort_builder.py` |

**Step 2 — RTDOSE + GTV t0 NIfTI (~52 GB, IBM Aspera):** `make download-rt`

**Step 3 — Follow-up GTV t1 NIfTI (optional, ~160 patients):** `make download-t1-gtv`  
_Status: blocked in CI/sandbox (TCP 33001); run locally with Aspera Connect._

---

## WEEK 1 — Infrastructure and Data Acquisition

### Coding

- [x] `[INFRA]` Initialize repository structure: `src/`, `notebooks/`, `data/`, `figures/`, `reports/`
- [x] `[INFRA]` Write `requirements.txt`
- [x] `[INFRA]` Write `.gitignore`
- [x] `[INFRA]` Write `README.md`
- [x] `[INFRA]` Write `src/config.py` with central paths and constants
- [x] `[INFRA]` Write `src/models/base_model.py` — abstract `TCPModel` base class
- [x] `[CODE]` Write `src/data/cohort_builder.py` — merge TSVs, apply inclusion criteria, EQD2, export `cohort.csv`
- [x] `[CODE]` Write stubs with docstrings: `src/data/nifti_loader.py`, `src/data/dvh_calculator.py`
- [x] `[CODE]` Download clinical TSVs → `data/processed/` (`download_clinical_data.py`)
- [x] `[CODE]` Run `python -m src.data.cohort_builder` → 264 total, 190 included, 74 excluded
- [x] `[CODE]` Write `notebooks/01_cohort_overview.ipynb`: demographics, fractionation, exclusion summary
- [x] `[CODE]` Download RT NIfTI (191 patients) — `download_rt_connect.py`, `organize_raw_data.py`
- [x] `[CODE]` Verify raw data completeness — `verify_raw_data.py`
- [x] `[CODE]` Results reporting — `src/reporting/update_results.py` → `reports/RESULTS.md` + `reports/metrics/*.csv` (`make report`)

### Literature and Writing

- [x] `[LIT]` Search PubMed for TCP modeling in GBM (≥8 papers; see `reports/literature_table.csv`, 12 entries)
- [x] `[LIT]` Build `reports/literature_table.csv` (18 refs, verified DOIs, ref_id 1–18)
- [x] `[WRITE]` In-text citations throughout `manuscript_draft.md` (Abstract–Discussion–slides)
- [x] `[LIT]` Read CFB-GBM data descriptor (Moreau et al. 2025, Front Oncol; PMID 39949753)
- [x] `[WRITE]` Draft Introduction — in `reports/manuscript_draft.md` §1 (full polish pending)

---

## WEEK 2 — Dosimetric Feature Extraction

### Coding

- [x] `[CODE]` Implement `src/data/nifti_loader.py`: `load_rtdose()` and `load_gtv_mask()` using `nibabel`
- [x] `[CODE]` Implement `src/data/dvh_calculator.py`: `compute_dvh()` and `extract_dvh_metrics()` (D95, D98, D50, D2, Dmean, Dmax, volume, Vx, gEUD, HI)
- [x] `[CODE]` Write `src/data/feature_builder.py`: iterate cohort, extract DVH features, merge, export `features.csv` (191 rows; gitignored)
- [x] `[CODE]` DVH quality control — `dvh_qc.py` (exclude Dmean < 1 Gy → patient 32)
- [x] `[CODE]` Export modeling dataset — `export_modeling_dataset.py` → `modeling_table.csv` (190 × 58 after v3)
- [x] `[CODE]` Write `src/utils/plot_dvh.py`: `plot_dvh_overlay()`, 300 dpi, publication style
- [x] `[CODE]` Write `notebooks/02_feature_extraction.ipynb`: descriptive stats, DVH overlay plots, distribution figures

### Literature and Writing

- [x] `[LIT]` Identify outcome definitions across papers vs CFB-GBM (OS / LC / PFS / RANO); documented in manuscript §4.1 and confounding audit
- [x] `[WRITE]` Draft Dataset Description — `manuscript_draft.md` §2.1 + `RESULTS.md` §1
- [x] `[WRITE]` Draft Data Acquisition — `manuscript_draft.md` §2.1 + README Data section
- [x] `[WRITE]` Draft Feature Extraction — `manuscript_draft.md` §2.2 + `RESULTS.md` DVH metrics

---

## WEEK 3 — TCP Models, Parameter Estimation, Survival Analysis

### Coding

- [x] `[CODE]` Write `src/models/poisson_tcp.py` — Poisson TCP model, MLE via `scipy.optimize.minimize`
- [x] `[CODE]` Write `src/models/logistic_tcp.py` — logistic dose-response, MLE
- [x] `[CODE]` Write `src/models/probit_tcp.py` — probit model, MLE
- [x] `[CODE]` Write `src/models/eud_tcp.py` — gEUD-based TCP with volume effect parameter `a`
- [x] `[CODE]` Write `src/models/bootstrap_ci.py` — bootstrap resampling (n=1000, seed=42) for any `TCPModel`
- [x] `[CODE]` Write `src/models/model_comparison.py` — AIC, BIC, log-likelihood, ROC AUC, calibration plot
- [x] `[CODE]` Write `src/models/survival_analysis.py` — Kaplan-Meier (OS), Cox PH regression with DVH covariates
- [x] `[CODE]` Write `notebooks/03_tcp_models.ipynb`: OS proxy + RANO v3, within-arm analysis
- [x] `[CODE]` Write `notebooks/04_parameter_estimation.ipynb`: bootstrap CI, model comparison table
- [x] `[CODE]` Write `notebooks/05_survival_analysis.ipynb`: KM curves, Cox, clinical stratification
- [x] `[CODE]` `src/analysis/confounding_audit.py` — dose heterogeneity, confounding correlations
- [x] `[CODE]` `src/analysis/stratified_analysis.py` — clinical prognosis, within-arm Spearman
- [x] `[CODE]` `src/analysis/within_arm_rano_tcp.py` — per-scheme DVH → RANO + Cox OS~RANO
- [x] `[CODE]` `src/analysis/rano_tcp_comparison.py` — OS vs RANO endpoint comparison
- [x] `[CODE]` `src/analysis/rano_multivariable.py` — 40 Gy multivariable logistic + volume validation
- [x] `[CODE]` `src/analysis/rano_prediction_suite.py` — pooled RANO, LOOCV, PyRadiomics comparison
- [x] `[CODE]` `src/data/rano_loader.py` + v3 cohort merge in `cohort_builder.py`
- [x] `[CODE]` Write `notebooks/06_rano_multivariable_40gy.ipynb`
- [x] `[CODE]` `scripts/build_analysis_notebooks.py` — regenerate notebooks 03–06

### Literature and Writing

- [x] `[WRITE]` Draft TCP equations — `manuscript_draft.md` §2.3 + `reports/manuscript_equations.tex`
- [x] `[WRITE]` Draft Parameter Estimation — `manuscript_draft.md` §2.3–2.5
- [x] `[WRITE]` Draft Results — `manuscript_draft.md` §3 + `reports/RESULTS.md`
- [x] `[WRITE]` Draft Discussion — `manuscript_draft.md` §4 (expanded 2026-06-28)
- [x] `[WRITE]` Draft Conclusion — `manuscript_draft.md` §5
- [x] `[WRITE]` Manuscript skeleton + presentation outline — `reports/manuscript_draft.md`

---

## WEEK 4 — Literature Comparison, Report, Presentation

### Coding

- [x] `[CODE]` Export analysis figures to `figures/` (22 PNG, 300 dpi; PDF export pending)
- [ ] `[CODE]` Download t1 GTV NIfTI + validate `size_t1_cm3` (`make download-t1-gtv`; Aspera required)
- [x] `[CODE]` Final check: all notebooks 01–06 run top-to-bottom without errors (`make check-notebooks`; 01, 03–06 verified 2026-06-28)
- [x] `[CODE]` Pin package versions in `requirements.txt`
- [x] `[CODE]` Nested CV for PyRadiomics feature selection (`pyradiomics_nested_cv_rano.csv`)
- [x] `[LIT]` Literature TCP D50 comparison table (`literature_tcp_d50_comparison.csv`)
- [x] `[WRITE]` LaTeX/Word/PDF export (`scripts/export_manuscript.sh` → manuscript.tex/docx/pdf)

### Literature and Writing

- [x] `[LIT]` Literature table with TCP/GBM/radiomics comparators — `reports/literature_table.csv`
- [x] `[LIT]` Critical Discussion draft — `manuscript_draft.md` §4 (limitations, confounding, Moreau comparison)
- [x] `[LIT]` Literature comparison in Discussion — D50 bootstrap vs literature; feasibility framing vs Ohri/Maitre
- [x] `[WRITE]` Merge manuscript sections into exportable document (`manuscript_draft.md` → `.docx` / `.pdf`)
- [x] `[WRITE]` Figure captions — `reports/figure_captions.md`
- [x] `[WRITE]` Assignment-style formal report — `assignment_report.md` → `.docx`/`.pdf` (`make export-assignment`)
- [x] `[LIT]` DOI/PubMed link verification — `scripts/verify_literature_dois.py` → `literature_doi_check.md`
- [x] `[WRITE]` Final polish of report document (10–15 pages; layout, figures embedded)
- [ ] `[LIT]` Format references Vancouver style in `manuscript.docx` (done in `assignment_report.pdf`)
- [ ] `[WRITE]` Build oral presentation slides from outline (15 min; outline in `manuscript_draft.md`)

---

## Remaining team deliverables

Shared resources: [`reports/RESULTS.md`](reports/RESULTS.md) (verified numbers) · [`reports/group_glossary_guide.md`](reports/group_glossary_guide.md) · [`figures/`](figures/)

### `[WRITE]` — report & presentation

- [ ] Build 15-min slide deck from presentation outline (`manuscript_draft.md`, bottom)
- [x] `[WRITE]` Write figure captions for all figures used in slides and final report — `figure_captions.md`
- [ ] Proofread and polish Introduction, Discussion, and Abstract in `manuscript_draft.md`
- [x] `[WRITE]` Assemble formal assignment report PDF — `assignment_report.pdf` (`make export-assignment`)

### `[LIT]` — literature & references

- [x] Format reference list Vancouver style in assignment report (`assignment_report.pdf` / `.docx`)
- [x] `[LIT]` Verify PubMed / DOI links — `make verify-dois` → `literature_doi_check.md`
- [ ] Draft one short paragraph per key reference for slides or appendix (optional)

### All members — presentation prep

- [ ] Read `group_glossary_guide.md` Parts 1–3 (project story) and Part 13 (FAQ)
- [ ] Assign presentation sections (cohort / TCP / RANO / radiomics / conclusions)
- [ ] Rehearse using talking points in glossary Part 15

---

## Remaining / Optional

- [ ] `[CODE]` External validation on second TCIA GBM cohort
- [ ] `[CODE]` Δvolume (t0→t1 NIfTI) → RANO after t1 GTV download
- [ ] `[WRITE]` Submit short paper / ESTRO poster (target: *Phys Med* or *Med Phys* technical note)

---

## Dependency Map

```
INFRA (repo, config) ✓
  └── DATA (clinical TSVs v3 → cohort_builder + rano_loader) ✓
        └── DOWNLOAD (RT NIfTI t0, verify) ✓
              └── EXTRACT (nifti_loader → dvh_calculator → feature_builder) ✓
                    └── QC + EXPORT (dvh_qc → modeling_table.csv 190×58) ✓
                          ├── PLOT (plot_dvh.py) ✓
                          ├── NOTEBOOK 01–02 ✓
                          ├── MODEL-01..05 (TCP models) ✓
                          │     └── STATS (bootstrap CI, model comparison) ✓
                          ├── STATS (survival analysis) ✓
                          ├── v3 RANO pipeline ✓
                          │     ├── rano_tcp_comparison, within_arm_rano_tcp ✓
                          │     ├── rano_multivariable (40 Gy) ✓
                          │     └── rano_prediction_suite (pooled, LOOCV, PyRadiomics) ✓
                          ├── NOTEBOOK 03–06 ✓
                          └── REPORT (update_results → RESULTS.md) ✓

Optional: download t1 GTV → validate_rano_volumes (t1) — pending Aspera

LIT tasks ── literature_table ✓; Vancouver refs in assignment report ✓; manuscript.docx refs pending
WRITE ── manuscript export ✓; final report polish + slides pending
Group ── glossary guide ✓; presentation prep pending
```

---

## Coding Standards

- All code in **English**: variable names, functions, docstrings, comments, commit messages
- Every function: **NumPy-style docstring** with Parameters, Returns, Raises
- No hardcoded paths — use `src/config.py` with `pathlib.Path`
- Notebooks must run **top-to-bottom** (`Kernel → Restart & Run All`)
- Commits: imperative, lowercase — e.g. `implement dvh_calculator with D95 D50 D2 extraction`
- **After each analysis commit:** run `make report` and include updated `reports/RESULTS.md` + `reports/metrics/*.csv`
- One feature per branch; merge via pull request
