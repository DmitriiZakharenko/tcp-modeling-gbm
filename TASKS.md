# Task Board — TCP Modeling in GBM

Tasks are labeled by type and role. Roles: **[CODE]** (any developer), **[LIT]** (literature/clinical review), **[WRITE]** (report writing).

Status: `[ ]` open · `[~]` in progress · `[x]` done

---

## Dataset Facts (CFB-GBM, TCIA)

- 264 patients; 194 have RTDOSE; 191 have GTV segmentation
- **Format: NIfTI (.nii.gz), NOT DICOM** — images pre-processed, skull-stripped, co-registered
- RTDOSE resolution: 2.5 × 2.5 × 5 mm, ~167 × 95 × 168 voxels
- Two fractionation schemes: 60 Gy / 30 fr (standard Stupp) and 40.05 Gy / 15 fr (elderly / poor PS)
- ~26.5% of patients have unknown dose — exclude from TCP modeling
- Clinical outcome (overall survival in weeks): available in TSV, all 264 patients
- **Clinical TSV files total ~50 KB** — download immediately, independent of images

## Data Strategy — No HPC Required

The 208 GB figure is for the **entire** collection: all MRI sequences (T1Gd, T2-FLAIR, T2*, ADC, etc.) at three time points, plus CT and RTDOSE.

For TCP modeling we need only:
- `RTDOSE` NIfTI files — estimated ~3–5 MB per patient × 194 ≈ **~0.6–1 GB**
- `GTV` segmentation NIfTI masks — estimated ~0.5–1 MB per patient × 191 ≈ **~100–200 MB**
- All clinical/treatment TSV files — **~50 KB total**, download now

**Total required: ~1–1.5 GB. A laptop is sufficient. HPC not needed.**

Download via IBM Aspera with selective file filtering (RTDOSE and GTV folders only). See task `[DATA-01]`.

---

## Labels

- `[DATA]` — data acquisition and organization
- `[EXTRACT]` — dosimetric and clinical feature extraction
- `[MODEL]` — TCP model implementation
- `[STATS]` — parameter estimation, CI, model comparison
- `[VIZ]` — figures and plots
- `[INFRA]` — repository infrastructure
- `[LIT]` — literature and clinical review (no coding required)
- `[WRITE]` — report sections (no coding required)

---

## WEEK 1 — Infrastructure, Data, and Literature Start

### Coding

- [ ] `[INFRA]` **[CODE]** Initialize repository structure:
  ```
  src/data/, src/models/, src/utils/
  notebooks/, data/raw/, data/processed/, figures/, reports/
  ```

- [ ] `[INFRA]` **[CODE]** Write `requirements.txt` with pinned versions:
  `nibabel`, `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`,
  `scikit-learn`, `lifelines`, `notebook`, `tqdm`

- [ ] `[INFRA]` **[CODE]** Write `.gitignore`:
  exclude `data/raw/`, `data/processed/*.nii.gz`, `__pycache__`, `.env`, `.DS_Store`

- [ ] `[INFRA]` **[CODE]** Write `src/config.py`: central paths using `pathlib.Path`
  (`DATA_RAW`, `DATA_PROCESSED`, `FIGURES_DIR`); read from environment or defaults

- [ ] `[INFRA]` **[CODE]** Write `README.md`: project summary, dataset citation, install
  instructions, data download note, usage

- [ ] `[DATA-01]` **[CODE]** Write `src/data/download_instructions.md` and
  `src/data/aspera_filter.txt`: document how to download only RTDOSE and GTV
  folders via IBM Aspera Connect from TCIA (selective download, ~1 GB total)

- [ ] `[DATA-02]` **[CODE]** Download all clinical TSV files from TCIA
  (clinical_data, treatment_data, treatment_imaging_availability, ct_availability, data_dictionary)
  and commit to `data/processed/` — these are tiny (~50 KB), safe to version-control

- [ ] `[DATA-03]` **[CODE]** Write `src/data/cohort_builder.py`:
  load clinical and treatment TSVs; merge on patient ID; apply inclusion criteria
  (has RTDOSE=True, has GTV=True, dose not unknown); export `data/processed/cohort.csv`
  with columns `[patient_id, rt_dose_gy, n_fractions, survival_weeks, age, sex, who_status, included, exclusion_reason]`

- [ ] `[DATA-04]` **[CODE]** Write `notebooks/01_cohort_overview.ipynb`:
  load `cohort.csv`; print inclusion/exclusion summary; plot cohort demographics
  (age distribution, sex, WHO status, dose scheme breakdown); save figures

### Literature and Writing (start immediately, independent of coding)

- [ ] `[LIT-01]` **[LIT]** Search PubMed for TCP modeling in GBM / brain tumors.
  Query: `("tumor control probability" OR "TCP model") AND ("glioblastoma" OR "brain tumor") AND "radiotherapy"`
  Filter: after 2000, peer-reviewed, English. Target: ≥8 papers (5 required, 3 buffer).

- [ ] `[LIT-02]` **[LIT]** Build literature table `reports/literature_table.csv`:
  columns `[authors, year, journal, n_patients, dataset, model_type, D50_Gy, gamma50, CI_method, outcome_definition, fractionation, notes]`

- [ ] `[LIT-03]` **[LIT]** Read TCIA CFB-GBM data descriptor paper (Moreau et al. 2025,
  Frontiers in Oncology) and summarize dataset limitations relevant to TCP modeling

- [ ] `[WRITE-01]` **[WRITE]** Write report section 1: "Introduction and Clinical Background"
  (GBM epidemiology, Stupp protocol, rationale for TCP modeling, project aim; ~1 page)

---

## WEEK 2 — Dosimetric Feature Extraction

### Coding

- [ ] `[EXTRACT-01]` **[CODE]** Write `src/data/nifti_loader.py`:
  load RTDOSE and GTV NIfTI files for a single patient using `nibabel`;
  return dose array (Gy), GTV mask array, and affine/voxel spacing metadata

- [ ] `[EXTRACT-02]` **[CODE]** Write `src/data/dvh_calculator.py`:
  compute cumulative DVH from dose array and binary mask; extract
  D95, D98, D50, D2, Dmean, Dmax, and GTV volume (cm³); return `dict`

- [ ] `[EXTRACT-03]` **[CODE]** Write `src/data/feature_builder.py`:
  iterate over all included patients; call `nifti_loader` and `dvh_calculator`;
  merge DVH features with `cohort.csv`; export `data/processed/features.csv`

- [ ] `[EXTRACT-04]` **[CODE]** Write `notebooks/02_feature_extraction.ipynb`:
  load `features.csv`; print descriptive statistics; check missing values;
  plot DVH curves for a sample of patients (overlay, one color per patient);
  plot distributions of D95, Dmean, GTV volume; save all figures

- [ ] `[VIZ-01]` **[CODE]** Write `src/utils/plot_dvh.py`:
  reusable `plot_dvh_overlay(patient_ids, dvh_data, save_path)` function;
  publication-quality style (300 dpi, consistent font size 11pt, axis labels in Gy and %)

### Literature and Writing

- [ ] `[LIT-04]` **[LIT]** Identify how outcome is defined across literature papers
  (local control vs. OS vs. PFS) and flag mismatches with CFB-GBM (OS in weeks only)

- [ ] `[WRITE-02]` **[WRITE]** Write report section 2: "Dataset Description and Patient Cohort"
  (CFB-GBM provenance, cohort demographics table, inclusion/exclusion criteria; ~1 page)

- [ ] `[WRITE-03]` **[WRITE]** Write report section 3: "Data Acquisition and Curation Methodology"
  (TCIA access, Aspera download, NIfTI format note, selective download rationale; ~0.5 page)

---

## WEEK 3 — TCP Models, Parameter Estimation, Survival Analysis

### Coding

- [ ] `[MODEL-01]` **[CODE]** Write `src/models/base_model.py`:
  abstract base class `TCPModel` with methods:
  `fit(doses, outcomes)`, `predict(doses)`, `log_likelihood(params, doses, outcomes)`, `summary()`

- [ ] `[MODEL-02]` **[CODE]** Write `src/models/poisson_tcp.py`:
  Poisson TCP: `TCP = exp(-N0 * exp(-alpha * d))`;
  MLE fitting via `scipy.optimize.minimize` (L-BFGS-B);
  parameters: D50, gamma50 (reparameterized for stability)

- [ ] `[MODEL-03]` **[CODE]** Write `src/models/logistic_tcp.py`:
  logistic dose-response: `TCP = 1 / (1 + exp(-k*(d - D50)))`;
  MLE fitting; parameters: D50, k (→ gamma50)

- [ ] `[MODEL-04]` **[CODE]** Write `src/models/probit_tcp.py`:
  probit model using `scipy.stats.norm.cdf`;
  MLE fitting; parameters: mu, sigma (→ D50, gamma50)

- [ ] `[MODEL-05]` **[CODE]** Write `src/models/eud_tcp.py`:
  gEUD from DVH: `EUD = (sum(v_i * d_i^a))^(1/a)`;
  combine with Poisson TCP on EUD instead of mean dose;
  parameters: D50, gamma50, a (volume effect)

- [ ] `[STATS-01]` **[CODE]** Write `src/models/bootstrap_ci.py`:
  bootstrap resampling (n=1000, seed=42) for any `TCPModel` subclass;
  return `pd.DataFrame` with columns `[param, estimate, ci_lower_95, ci_upper_95]`

- [ ] `[STATS-02]` **[CODE]** Write `src/models/model_comparison.py`:
  compute AIC, BIC, log-likelihood, ROC AUC for a list of fitted models;
  return comparison `pd.DataFrame`; generate calibration plot (predicted vs. observed)

- [ ] `[STATS-03]` **[CODE]** Write `src/models/survival_analysis.py`:
  Kaplan-Meier curve for OS using `lifelines.KaplanMeierFitter`;
  Cox PH regression with DVH features (D95, Dmean, GTV volume, D50) as covariates;
  export hazard ratios with 95% CI as `pd.DataFrame`; forest plot

- [ ] `[EXTRACT-05]` **[CODE]** Handle fractionation heterogeneity:
  add EQD2 correction column to `features.csv`
  (EQD2 = D_total * (d_fraction + alpha_beta) / (2 + alpha_beta), alpha/beta = 10 Gy for GBM)
  to allow pooling of 60Gy/30fr and 40.05Gy/15fr cohorts

- [ ] `[CODE]` Write `notebooks/03_tcp_models.ipynb`:
  fit all four models; print parameter tables; plot TCP curves per model on one figure

- [ ] `[CODE]` Write `notebooks/04_parameter_estimation.ipynb`:
  bootstrap CI for all models; comparison table AIC/BIC/AUC; calibration plots

- [ ] `[CODE]` Write `notebooks/05_survival_analysis.ipynb`:
  KM curves; Cox regression table; forest plot of hazard ratios

### Literature and Writing

- [ ] `[WRITE-04]` **[WRITE]** Write report section 5: "Mathematical Description of TCP Models"
  (equations for all four models with notation; reference Okunieff 1995; ~1.5 pages)

- [ ] `[WRITE-05]` **[WRITE]** Write report section 6: "Parameter Estimation Methodology"
  (MLE rationale, bootstrap procedure, EQD2 correction, fractionation pooling; ~1 page)

---

## WEEK 4 — Literature Comparison, Report, Presentation

### Coding

- [ ] `[INFRA]` **[CODE]** Final reproducibility check: all notebooks run top-to-bottom
  (`Kernel → Restart & Run All`) without errors on a clean environment

- [ ] `[VIZ-02]` **[CODE]** Export all publication-quality figures:
  300 dpi PNG + PDF, consistent style (font 11pt, color-blind-safe palette),
  save to `figures/`; update figure captions in `figures/README.md`

- [ ] `[INFRA]` **[CODE]** Pin all versions in `requirements.txt` via `pip freeze > requirements.txt`;
  add `## Reproducibility` section to README with Python version and OS

### Literature and Writing

- [ ] `[LIT-05]` **[LIT]** Write report section 9: "Literature Review and Comparison"
  (compare our D50/gamma50 with values from literature table; table format; ~2 pages)

- [ ] `[LIT-06]` **[LIT]** Write report section 10: "Critical Discussion":
  — CFB-GBM limitations: single-center, OS-only outcome (no local control), ~26.5% missing dose data
  — Model assumptions: homogeneous cell radiosensitivity, Poisson statistics, independence of clonogens
  — Fractionation heterogeneity and EQD2 correction validity
  — Small effective cohort size → wide CI → impact on generalizability
  — How our results compare to published TCP parameters and why they might differ

- [ ] `[WRITE-06]` **[WRITE]** Write report sections 4, 7, 8:
  "DVH Feature Extraction" (methods), "Results and Visualizations" (describe figures),
  "Model Comparison and Statistical Evaluation" (interpret AIC/BIC/AUC table)

- [ ] `[WRITE-07]` **[WRITE]** Write report section 11: "Conclusion" (~0.5 page)

- [ ] `[WRITE-08]` **[WRITE]** Assemble full report from all sections;
  format references (Vancouver style, numbered); export as PDF (10–15 pages)

- [ ] `[WRITE-09]` **[WRITE]** Build oral presentation (15 min, ~15 slides):
  data pipeline → cohort → DVH → TCP models → results → literature comparison → clinical interpretation;
  include key figures from notebooks

---

## Dependency Map

```
INFRA (repo, config)
  └── DATA-02 (clinical TSVs — download now)
        └── DATA-03 (cohort_builder.py)
              └── EXTRACT-01/02/03 (NIfTI loader → DVH calculator → feature_builder)
                    ├── EXTRACT-05 (EQD2 correction)
                    ├── MODEL-01 → MODEL-02/03/04/05 (TCP models)
                    │     └── STATS-01/02 (bootstrap CI, model comparison)
                    └── STATS-03 (survival analysis)

LIT-01 → LIT-02/03/04 ── independent, start Week 1
WRITE-01/02/03 ── independent, start Week 1
WRITE-04/05 ── after models are designed (Week 3 start)
WRITE-06/07/08/09 ── after results are available (Week 4)
```

---

## Coding Standards

- All code in **English**: variable names, function names, docstrings, comments, commit messages
- Every function must have a **NumPy-style docstring**: Parameters, Returns, Raises sections
- No hardcoded paths — use `src/config.py` with `pathlib.Path`
- Notebooks must run **top-to-bottom without errors** (`Kernel → Restart & Run All`)
- Git commits: imperative, lowercase, specific — e.g. `add eud tcp model with gEUD from DVH`
- One feature per branch; merge via pull request
