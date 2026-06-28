# Task Board — TCP Modeling in GBM

## How to use this board

Checkboxes below are **clickable directly on GitHub** (no terminal needed):
- Open this file on [github.com/DmitriiZakharenko/tcp-modeling-gbm](https://github.com/DmitriiZakharenko/tcp-modeling-gbm/blob/main/TASKS.md)
- Click any `[ ]` to mark it done — GitHub saves automatically
- **Coders:** check off tasks as you push the code
- **Literature / Writing team:** check off your tasks as you complete each section

Roles: **[CODE]** · **[LIT]** literature & clinical review · **[WRITE]** report writing

---

## Dataset Facts (CFB-GBM, TCIA) — verified 2026-06-26

| Metric | Value | Source |
|---|---|---|
| Total patients in TSV | 264 | `cohort.csv` |
| Has RTDOSE (imaging flag) | 194 | `cohort.csv` |
| Has GTV segmentation | 191 | `cohort.csv` |
| Both RTDOSE + GTV on disk | 191 | `verify_raw_data.py` |
| Unknown RT dose (excluded) | 70 | `cohort.csv` exclusion_reason |
| **Included after cohort rules** | **190** | `cohort.csv` (`included=True`) |
| **Modeling table (post DVH QC)** | **190 × 33** | `modeling_table.csv` |
| DVH QC exclusion | 1 (patient 32, Dmean=0 Gy) | `dvh_qc.py` |
| Primary fractionation | 120 × 60 Gy/30 fr; 61 × 40.05 Gy/15 fr | `modeling_table.csv` |
| Median age (modeling) | 70 yr | `modeling_table.csv` |
| Sex M/F (modeling) | 117 / 73 | `modeling_table.csv` |
| Median OS (modeling) | 51 wk | `modeling_table.csv` |

- Format: **NIfTI (.nii.gz)**, pre-processed (skull-stripped, co-registered) — not DICOM
- Two fractionation schemes dominate: 60 Gy / 30 fr (standard Stupp) and 40.05 Gy / 15 fr (elderly / poor PS)
- ~26.5% of patients have unknown dose → excluded from TCP modeling
- Clinical outcome: overall survival in weeks, available for all 264 patients in TSV

## Data Download Strategy

**Step 1 — Clinical TSVs (~50 KB total, no login, no Aspera): download now**

Save all files to `data/processed/`:

| File | Size | Direct URL |
|---|---|---|
| Clinical data | 12.4 KB | https://www.cancerimagingarchive.net/wp-content/uploads/CFB-GBM_clinical_data_v02_20260129.tsv |
| Treatment data | 6.86 KB | https://www.cancerimagingarchive.net/wp-content/uploads/CFB-GBM_treatment_data_v02_20260129.tsv |
| Imaging availability | 5.47 KB | https://www.cancerimagingarchive.net/wp-content/uploads/CFB-GBM_treatment_imaging_availability_v02_20260129.tsv |
| Data dictionary | 4.09 KB | https://www.cancerimagingarchive.net/wp-content/uploads/CFB-GBM_columns_description_new_v02_20260129.tsv |

**Step 2 — RTDOSE and GTV NIfTI files (~52 GB on disk for 191 patients, requires IBM Aspera Connect)**

1. Install [IBM Aspera Connect](https://www.ibm.com/products/aspera/downloads)
2. Go to [CFB-GBM on TCIA](https://www.cancerimagingarchive.net/collection/cfb-gbm/)
3. Use `python -m src.data.download_rt_connect` (Faspex OAuth + Connect API) or manual Aspera dialog
4. Select only RTDOSE and GTV subfolders
5. Download to `data/raw/`; verify with `python -m src.data.verify_raw_data`

Do **not** download MRI sequences — they are not needed for this project.

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
- [x] `[CODE]` Download clinical TSVs → `data/processed/` (4 files; `download_clinical_data.py`)
- [x] `[CODE]` Run `python -m src.data.cohort_builder` → 264 total, 190 included, 74 excluded
- [x] `[CODE]` Write `notebooks/01_cohort_overview.ipynb`: demographics, fractionation, exclusion summary
- [x] `[CODE]` Download RT NIfTI (191 patients) — `download_rt_connect.py`, `organize_raw_data.py`
- [x] `[CODE]` Verify raw data completeness — `verify_raw_data.py`
- [x] `[CODE]` Results reporting — `src/reporting/update_results.py` → `reports/RESULTS.md` + `reports/metrics/*.csv` (run via `make report` before each analysis commit)

### Literature and Writing (start now, independent of coding)

- [ ] `[LIT]` Search PubMed for TCP modeling in GBM.
  Query: `("tumor control probability" OR "TCP model") AND ("glioblastoma" OR "brain tumor") AND "radiotherapy"`
  Filter: after 2000, peer-reviewed, English. Target ≥8 papers.
- [ ] `[LIT]` Build `reports/literature_table.csv`:
  columns `[authors, year, journal, n_patients, dataset, model_type, D50_Gy, gamma50, CI_method, outcome_definition, fractionation, notes]`
- [ ] `[LIT]` Read CFB-GBM data descriptor paper (Moreau et al. 2025, Frontiers in Oncology);
  summarize dataset limitations relevant to TCP modeling
- [ ] `[WRITE]` Write report section 1: Introduction and Clinical Background (~1 page)

---

## WEEK 2 — Dosimetric Feature Extraction

### Coding

- [x] `[CODE]` Implement `src/data/nifti_loader.py`: `load_rtdose()` and `load_gtv_mask()` using `nibabel`
- [x] `[CODE]` Implement `src/data/dvh_calculator.py`: `compute_dvh()` and `extract_dvh_metrics()` (D95, D98, D50, D2, Dmean, Dmax, volume, Vx, gEUD, HI)
- [x] `[CODE]` Write `src/data/feature_builder.py`: iterate cohort, extract DVH features, merge, export `data/processed/features.csv` (191 rows; local/gitignored)
- [x] `[CODE]` DVH quality control — `dvh_qc.py` (exclude Dmean < 1 Gy → patient 32)
- [x] `[CODE]` Export modeling dataset — `export_modeling_dataset.py` → `modeling_table.csv` (190 × 33, in git)
- [x] `[CODE]` Write `src/utils/plot_dvh.py`: `plot_dvh_overlay(patient_ids, dvh_data, save_path)`, 300 dpi, publication style
- [x] `[CODE]` Write `notebooks/02_feature_extraction.ipynb`: descriptive stats, DVH overlay plots, distribution figures

### Literature and Writing

- [ ] `[LIT]` Identify how outcome is defined across papers (local control vs. OS vs. PFS);
  flag mismatches with CFB-GBM (OS in weeks only — no local control endpoint)
- [ ] `[WRITE]` Write report section 2: Dataset Description and Patient Cohort (~1 page)
- [ ] `[WRITE]` Write report section 3: Data Acquisition and Curation Methodology (~0.5 page)
- [ ] `[WRITE]` Write report section 4: Dosimetric and Clinical Feature Extraction (~1 page)

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
- [x] `[CODE]` `src/analysis/within_arm_rano_tcp.py` — per-scheme DVH → RANO TCP + Cox OS~RANO
- [x] `[CODE]` Write `notebooks/06_rano_multivariable_40gy.ipynb`: multivariable logistic 40 Gy arm
- [x] `[WRITE]` Draft `reports/manuscript_draft.md` + `reports/literature_table.csv`

### Literature and Writing

- [ ] `[WRITE]` Write report section 5: Mathematical Description of TCP Models (~1.5 pages, include equations)
- [ ] `[WRITE]` Write report section 6: Parameter Estimation Methodology (~1 page)

---

## WEEK 4 — Literature Comparison, Report, Presentation

### Coding

- [ ] `[CODE]` Export all figures: 300 dpi PNG + PDF, consistent style, save to `figures/`
- [ ] `[CODE]` Final check: all notebooks run top-to-bottom without errors (`Kernel → Restart & Run All`)
- [ ] `[CODE]` Pin package versions: `pip freeze > requirements.txt`

### Literature and Writing

- [ ] `[LIT]` Write report section 9: Literature Review and Comparison
  (compare D50/γ50 estimates with literature table; ~2 pages)
- [ ] `[LIT]` Write report section 10: Critical Discussion
  (dataset limitations, model assumptions, fractionation heterogeneity, generalizability; ~2 pages)
- [ ] `[WRITE]` Write report sections 7–8: Results and Visualizations, Model Comparison (~2 pages)
- [ ] `[WRITE]` Write report section 11: Conclusion (~0.5 page)
- [ ] `[WRITE]` Assemble full report (10–15 pages), format references (Vancouver style), export PDF
- [ ] `[WRITE]` Build oral presentation (15 min, ~15 slides)

---

## Dependency Map

```
INFRA (repo, config) ✓
  └── DATA (clinical TSVs → cohort_builder) ✓
        └── DOWNLOAD (RT NIfTI, verify) ✓
              └── EXTRACT (nifti_loader → dvh_calculator → feature_builder) ✓
                    └── QC + EXPORT (dvh_qc → modeling_table.csv) ✓
                          ├── PLOT (plot_dvh.py) — in progress
                          ├── NOTEBOOK 02 — next
                          ├── MODEL-01..05 (TCP models)
                          │     └── STATS (bootstrap CI, model comparison)
                          └── STATS (survival analysis)

LIT tasks ── independent, start Week 1
WRITE 1–4 ── independent, start Week 1
WRITE 5–6 ── after models are designed (Week 3 start)
WRITE 7–11 ── after results available (Week 4)
```

---

## Coding Standards

- All code in **English**: variable names, functions, docstrings, comments, commit messages
- Every function: **NumPy-style docstring** with Parameters, Returns, Raises
- No hardcoded paths — use `src/config.py` with `pathlib.Path`
- Notebooks must run **top-to-bottom** (`Kernel → Restart & Run All`)
- Commits: imperative, lowercase — e.g. `implement dvh_calculator with D95 D50 D2 extraction`
- **After each analysis commit:** run `python -m src.reporting.update_results` (or `make report`) and include updated `reports/RESULTS.md` + `reports/metrics/*.csv` in the commit
- One feature per branch; merge via pull request
