# CFB-GBM TCP Project — Plain-Language Guide & Glossary

**Audience:** Group members who are not programmers or statisticians.  
**Purpose:** Explain what we did, what every important term means, and how to talk about our results with clinicians, supervisors, or examiners.  
**Last updated:** 2026-06-28 · Numbers from [`reports/RESULTS.md`](RESULTS.md) (n=190 modeling cohort).

---

## Part 1 — The story in one page

### What is this project about?

We studied **glioblastoma (GBM)** — an aggressive brain tumour — using a **public patient dataset** from a French cancer centre (CFB-GBM, shared on TCIA). The clinical question behind the assignment was:

> *Can we build a mathematical model that links **radiotherapy dose** to **tumour control**?*

That type of model is called **TCP — Tumour Control Probability**.

### What we actually built

1. We downloaded and cleaned patient data (clinical tables + radiation dose maps + tumour outlines).
2. For each patient we computed **how dose was distributed inside the tumour** (DVH metrics).
3. We fitted **four standard TCP models** (Poisson, Logistic, Probit, gEUD) and compared them.
4. We tested whether dose predicts outcome — first using **overall survival (OS)**, then using **RANO imaging response** (a newer, better endpoint from dataset Version 3).
5. When classical dose–response failed, we looked at **tumour size (GTV volume)** and **radiomics** (texture features from MRI).
6. We used honest validation methods (**LOOCV**, **nested cross-validation**) so results are not inflated by overfitting.

### The main conclusion (in plain words)

- **Classical TCP on pooled dose does not work well on this dataset** — not because our code is wrong, but because the data were collected in routine care, not in a dose-escalation trial. Almost everyone got the planned dose; there is almost no dose variation *within* each treatment schedule.
- **Tumour burden (volume) and radiomics can predict early imaging response (RANO)** if we account for the fact that patients on 60 Gy and 40 Gy schedules are different (age, prognosis).
- We deliver a **reproducible open pipeline** and a **feasibility checklist** for future TCP studies.

---

## Part 2 — Step-by-step: what we did (workflow)

| Step | What happened | Why it matters |
|------|---------------|----------------|
| **1. Cohort selection** | Started with 264 patients in CFB-GBM; kept 190 with valid dose maps and tumour masks | Ensures every analysis uses trustworthy input |
| **2. DVH extraction** | Combined RTDOSE (3D dose grid) + GTV mask → dose statistics per patient | Converts imaging into numbers models can use |
| **3. TCP fitting** | Fit Poisson / Logistic / Probit / gEUD models; estimate D50, γ50; bootstrap CIs | Core assignment deliverable (Part IV–VI) |
| **4. Survival analysis** | Kaplan–Meier curves, Cox regression for OS | Shows known clinical prognostic factors (scheme, WHO PS) |
| **5. RANO integration** | Merged v3 TSV with RANO labels (137/190 patients) | Better endpoint than OS for early tumour response |
| **6. Confounding audit** | Checked age vs dose, scheme vs outcome | Explains why pooled dose models mislead |
| **7. Within-arm analysis** | Split by 60 Gy vs 40 Gy; test dose and volume separately | Dose spread exists *between* arms, not *within* |
| **8. Multivariable models** | Logistic regression: volume + age + PS + scheme → RANO | Adjusts for confounders |
| **9. Cross-validation** | LOOCV (leave-one-out) and nested 5-fold CV for radiomics | Honest performance estimates |
| **10. PyRadiomics comparison** | Compared author-provided MRI texture features vs DVH volume | Links our work to original CFB-GBM AI papers |
| **11. Reporting** | `make report` → `RESULTS.md`, figures, manuscript draft | Everything traceable and regeneratable |

---

## Part 3 — Clinical & disease terms

### Glioblastoma (GBM)
The most common malignant primary brain tumour in adults. Standard care: surgery + radiotherapy + chemotherapy (temozolomide).

### Overall Survival (OS)
Time from diagnosis (or start of treatment) until death. Measured in **weeks** in our cohort. Median OS = **51 weeks** (half of patients lived longer, half shorter).

### OS median-split proxy
We sometimes converted OS into a yes/no label: *Did the patient survive at least as long as the cohort median (51 weeks)?*  
This was an **exploratory stand-in** for “good outcome” when testing the TCP pipeline — **not** the same as true tumour control.

### WHO Performance Status (PS)
A clinician’s rating of how well a patient can perform daily activities (0 = fully active … 4 = bedridden). Higher PS → worse prognosis. In our Cox models, PS was a strong predictor (HR ≈ 1.42 per step).

### Fractionation / treatment scheme
How total radiation dose is split into daily sessions (fractions):

| Scheme | Meaning | n in cohort |
|--------|---------|-------------|
| **60 Gy / 30 fractions** | Standard curative schedule (~2 Gy per day) | 120 patients |
| **40.05 Gy / 15 fractions** | Hypofractionated schedule (~2.67 Gy per day), often for older or frail patients | 61 patients |

These two groups differ in age and prognosis — not a random split.

### Temozolomide (TMZ)
Chemotherapy drug given with radiotherapy. Mentioned in literature; our TCP models focus on radiotherapy dose metrics.

### RANO (Response Assessment in Neuro-Oncology)
Standardised rules for judging tumour response on MRI after treatment. Categories include:

| Code | Name | Plain meaning |
|------|------|---------------|
| **CR** | Complete Response | Tumour appears gone |
| **PR** | Partial Response | Clear shrinkage |
| **MR** | Minor Response | Small shrinkage |
| **SD** | Stable Disease | No meaningful change |
| **PD** | Progressive Disease | Tumour grew or new lesions |

### RANO non-PD (our binary endpoint)
**Yes** = SD, MR, PR, or CR at follow-up time **t1** (early post-treatment scan).  
**No** = PD.  
Event rate ≈ **77%** non-PD in the RANO subset (n=137). This is closer to “tumour controlled on imaging” than OS, but **still not the same as formal local control** (pathology or long-term recurrence).

### t0 and t1
- **t0** = baseline (before or at start of radiotherapy)
- **t1** = first follow-up imaging time point after treatment (early response assessment)

### Local Control (LC)
True endpoint for classical TCP: the tumour did not recur in the treated region. **Not available** in CFB-GBM as a labelled outcome — a major limitation for TCP validation.

### Confounding
When a third factor (e.g. age or treatment scheme) distorts the link between dose and outcome. Example: older patients got 40 Gy and had shorter OS — dose looks “protective” in pooled analysis, but age and intent of treatment drive the association.

---

## Part 4 — Data & imaging terms

### TCIA (The Cancer Imaging Archive)
Public repository where researchers share de-identified cancer imaging datasets. CFB-GBM DOI: [10.7937/v9pn-2f72](https://doi.org/10.7937/v9pn-2f72).

### CFB-GBM cohort
Dataset from Centre François Baclesse (France): MRI, radiotherapy dose, tumour contours, clinical tables. Version 3 (2025/2026) added RANO labels and PyRadiomics features.

### DICOM vs NIfTI
- **DICOM** — hospital standard format for CT/MRI/RT files (many files per patient).
- **NIfTI (.nii.gz)** — simplified 3D image format used in research.

**Important for our project:** TCIA provides **pre-processed NIfTI**, not full DICOM-RT plans. We do **not** have RTPLAN or RTSTRUCT DICOM in the public package.

### RTDOSE
A 3D map showing **how many Gray (Gy)** were delivered to each voxel (3D pixel) in the patient’s head/body. Think of it as a coloured heat map of radiation dose.

### GTV (Gross Tumour Volume)
The visible tumour outline drawn by radiation oncologists — the region that definitely contains tumour and must receive full dose.

### CTV / PTV
- **CTV** — Clinical Target Volume (GTV + margin for microscopic spread)
- **PTV** — Planning Target Volume (CTV + setup margin)

**Not included** in CFB-GBM public contours. All our DVH metrics refer to **GTV only**.

### Segmentation mask
A binary 3D image: inside GTV = 1, outside = 0. Used to tell the computer “only measure dose inside the tumour.”

### Skull-stripped / co-registered
Pre-processing steps already done by dataset authors: brain images cleaned of skull; MRI and dose maps aligned to the same coordinate system.

### Patient ID / modeling_table.csv
Our master spreadsheet: **190 rows × 58 columns** — one row per patient, columns for clinical data, DVH metrics, RANO labels, etc.

---

## Part 5 — Dosimetry terms (dose mathematics)

### Gray (Gy)
Unit of absorbed radiation dose. 1 Gy = 1 joule per kilogram of tissue.

### Dmean (mean dose to GTV)
Average dose inside the tumour volume. In the 60 Gy arm, Dmean varies only ~0.28 Gy between patients — essentially everyone received the planned dose.

### D95, D2, D98, Vx
DVH-derived metrics (see DVH below):
- **D95** — dose received by at least 95% of tumour volume (coverage metric)
- **D2** — dose to the hottest 2% (near-maximum dose)
- **Vx** — percentage of volume receiving at least x Gy

### DVH (Dose–Volume Histogram)
A curve answering: *“What fraction of the tumour received at least X Gy?”*  
Built by:
1. Take RTDOSE grid
2. Keep only voxels inside GTV mask
3. Count how many voxels exceed each dose level

**Input:** RTDOSE + GTV mask  
**Output:** curve + scalar summaries (Dmean, D95, volume, …)  
**Shows:** whether dose coverage is uniform or heterogeneous

### EQD2 (Equivalent Dose in 2 Gy fractions)
Converts different fractionation schedules to a **common scale**, as if dose were given in 2 Gy daily fractions. Uses the **linear-quadratic (LQ) model** with α/β = 10 Gy for GBM.

**Why we use it:** Lets us compare 60 Gy/30 fr and 40 Gy/15 fr on one axis.  
**Formula idea:** Bigger fractions biologically “weigh more” per Gy than 2 Gy fractions.

### α/β ratio
Radiobiology parameter describing how sensitive tissue is to fraction size. Tumour (GTV) often use α/β ≈ 10 Gy in GBM literature.

### gEUD (generalised Equivalent Uniform Dose)
Collapses the full DVH into **one number** that accounts for hot and cold spots inside the tumour. Parameter **a** controls volume effect:
- **a = 1** → similar to mean dose
- **a = −10** → emphasises cold (under-dosed) regions
- **a = +10** → emphasises hot spots

**Input:** full DVH  
**Output:** single dose metric per patient  
**Shows:** effective dose under different biological assumptions

### HI (Homogeneity Index)
Measures how uniform dose is inside GTV. Low variance in our cohort → not useful for predicting outcome.

### volume_cc
GTV volume in cubic centimetres (cc). Median ≈ **36 cc**. Strong predictor of RANO in several analyses.

---

## Part 6 — TCP models (the four we implemented)

**Common idea:** All TCP models predict **probability of “control”** (here: good binary outcome) as a function of dose (or gEUD). They differ in the **mathematical curve shape**.

**What they take as input:** A dose metric per patient (usually EQD2 or Dmean) + binary outcome (OS proxy or RANO non-PD).

**What they output:**
- Predicted TCP between 0 and 1 for each patient
- Fitted parameters (especially **D50** and **γ50** or equivalents)
- Goodness-of-fit statistics (AIC, AUC, …)

### D50
Dose at which predicted TCP = **50%**. Often reported in Gy.  
Our pooled Poisson estimate: **≈53 Gy** (bootstrap 95% CI ≈ 50–57 Gy) — interpret cautiously because endpoint is OS proxy, not true LC.

### γ50 (gamma50)
**Steepness** of the dose–response curve at D50. Higher γ → small dose changes cause large changes in TCP. Our estimate ≈ **3.3**.

### Poisson TCP
Based on **Poisson statistics of surviving clonogenic cells**. Assumes each cell has independent probability of being killed.

- **Parameters:** D50, γ50  
- **Typical use:** Classic radiobiology textbooks  
- **Our result:** AUC ≈ 0.68 for OS proxy; AUC ≈ 0.43 for RANO pooled

### Logistic TCP
Uses a **sigmoid (S-shaped) curve** — same family as logistic regression.

- **Parameters:** D50, k (slope)  
- **Our result:** Nearly identical to Poisson on this dataset (AIC within 0.1)

### Probit TCP
Uses the **normal cumulative distribution** (probit link) instead of logistic.

- **Parameters:** D50, σ (spread)  
- **Our result:** Best AIC among four models — but difference is tiny; all tell the same story

### gEUD TCP
First computes gEUD from DVH, then applies Poisson-style TCP on that single number.

- **Extra parameter:** a (volume effect)  
- **Our result:** Slightly higher AIC (3 parameters) — not clearly better

### Maximum Likelihood Estimation (MLE)
Standard way to find D50, γ50, etc.: choose parameters that make the observed outcomes **most probable** under the model.

### Bootstrap confidence interval
Resample patients with replacement 1000 times, refit model each time, take 2.5th–97.5th percentile of D50 and γ50. Shows uncertainty without assuming normal distribution.

---

## Part 7 — Statistical & machine-learning methods

### Kaplan–Meier (KM) curve
Non-parametric estimate of **survival over time**. Shows what fraction of patients are still alive at each week.

**Input:** survival time + event (died yes/no)  
**Output:** step curve + median survival  
**Our use:** Compare 60 Gy vs 40 Gy (median OS 60 vs 28 weeks, p ≈ 3×10⁻⁶)

### Log-rank test
Tests whether two KM curves differ significantly. Used for fractionation scheme comparison.

### Cox proportional hazards model
Regression for **time-to-event** data. Estimates **Hazard Ratio (HR)** — how much a factor increases/decreases instantaneous risk of death.

**Input:** survival time + covariates (age, PS, EQD2, RANO, …)  
**Output:** HR and p-value per covariate  
**Example:** RANO non-PD HR ≈ 0.48 → non-PD patients have ~half the hazard of PD patients, holding other variables fixed.

### Concordance index (C-index)
Cox model’s discrimination score (like AUC for survival): probability that patient who died earlier had higher predicted risk.

### Logistic regression
Predicts **probability of a yes/no outcome** (e.g. RANO non-PD) from one or more predictors.

**Our models:**
- Volume only
- Volume + age + WHO PS
- + treatment scheme indicator
- + volume × scheme interaction (does volume effect differ by arm?)

**Output:** coefficients, predicted probabilities, AUC

### Spearman correlation (ρ)
Measures **monotonic association** between two continuous variables (rank-based). Robust to outliers.

**Example:** GTV volume vs RANO non-PD in 40 Gy arm: ρ ≈ 0.41, p ≈ 0.016

### Mann–Whitney U test
Compares two groups on a continuous outcome (non-parametric). Used for OS difference between schemes.

### Likelihood ratio (LR) test
Compares full TCP model vs “dose has no effect” null model. Small p → dose adds significant information.

### Stratified / within-arm analysis
Analyse 60 Gy and 40 Gy **separately** so dose variation from different protocols does not confuse results.

### Pooled analysis
All patients together — useful for volume + covariate models, **misleading** for raw dose–TCP when schemes differ.

---

## Part 8 — Model evaluation metrics (how good is a model?)

### ROC curve & AUC (Area Under the Curve)
**ROC** plots true-positive rate vs false-positive rate at all classification thresholds.

**AUC** summary:
| AUC | Plain interpretation |
|-----|---------------------|
| 0.50 | Random guessing — useless |
| 0.60–0.70 | Weak discrimination |
| 0.70–0.80 | Acceptable / moderate |
| 0.80–0.90 | Good (watch for overfitting) |
| > 0.90 | Excellent — suspect small sample or in-sample optimism |

**Key results:**

| Analysis | In-sample AUC | Honest CV AUC |
|----------|---------------|---------------|
| Pooled EQD2 → OS proxy | 0.68 | 0.68 (5-fold) |
| Pooled EQD2 → RANO | 0.43 | 0.42 |
| Pooled volume + clinical + scheme → RANO | 0.72 | **0.64 (LOOCV)** |
| 40 Gy volume → RANO | 0.83–0.90 | **0.74 (LOOCV)** |
| PyRadiomics top-5 → RANO | 0.78 | **0.74 (nested CV)** |
| DVH volume → RANO | 0.71 | **0.70 (nested CV)** |

### Brier score
Mean squared error between predicted probability and actual 0/1 outcome. **Lower is better.** 0.25 = null model for 50% event rate.

### Sensitivity / Specificity
At threshold TCP = 0.5:
- **Sensitivity** — fraction of true “good outcomes” correctly flagged
- **Specificity** — fraction of true “bad outcomes” correctly flagged

### AIC / BIC
Balance model fit vs complexity. **Lower is better.** Used to compare Poisson vs Logistic vs Probit vs gEUD.

### McFadden pseudo-R²
Analogue of R² for logistic-type models. ~0.10 here → dose explains modest fraction of outcome variance.

### Hosmer–Lemeshow (HL) test
Calibration test: are predicted probabilities matched by observed frequencies in bins? Non-significant p → acceptable calibration.

### Optimism / optimism delta (Δ)
Difference between **in-sample AUC** (model tested on same data used to train) and **cross-validated AUC**. Large Δ → overfitting. PyRadiomics Δ ≈ 0.04–0.07.

---

## Part 9 — Cross-validation (honest performance)

### Why cross-validation?
If you train and test on the same patients, the model can **memorise noise**. CV holds out some patients during training and evaluates on them.

### 5-fold cross-validation
Split data into 5 parts; train on 4, test on 1; rotate. Report mean ± SD of AUC.

### LOOCV (Leave-One-Out Cross-Validation)
Special case: n−1 patients train, **1 patient** tests; repeat for every patient.  
**Pros:** Maximum data for training; deterministic  
**Cons:** Expensive; high variance for small n  
**Our use:** 40 Gy arm (n=34) → LOOCV AUC 0.74 for volume

### Nested cross-validation
When you **select features** (e.g. top-5 radiomics), selection must happen **inside** each training fold — otherwise information leaks from test set.

**Our nested 5-fold setup:**
1. Outer loop: 5 test folds
2. Inside each training fold: rank all radiomics features by univariate AUC → pick top 5 → train logistic model
3. Predict held-out fold
4. Aggregate AUC across all patients

**Result:** Radiomics nested AUC **0.74** vs volume **0.70** — real but modest advantage.

---

## Part 10 — PyRadiomics & imaging AI terms

### Radiomics
High-throughput extraction of **quantitative features** from medical images (shape, intensity, texture) — converts MRI into hundreds of numbers.

### PyRadiomics
Open-source Python library for radiomics (IBSI-compliant). **We used pre-computed features** from CFB-GBM v3 TSV (authors already extracted them).

### t1gd sequence
T1-weighted MRI **after gadolinium contrast** — bright enhancing tumour. Features computed on **GTV at t0**.

### Feature selection (top-5)
Rank all radiomics features by how well each alone predicts RANO (univariate AUC); keep best 5; combine in logistic regression.

### IBSI (Image Biomarker Standardisation Initiative)
Community standard so radiomics features are comparable across software and centres.

### dvh_volume vs pyro_mesh_volume
Both measure tumour size — nearly identical AUC (~0.71) because radiomics mesh volume correlates with DVH volume.

---

## Part 11 — Software & reproducibility terms

### Python pipeline
Scripts in `src/` automate download → DVH → modeling → figures. Not a black box: each step is a module you can run separately.

### Jupyter notebook
Interactive document mixing code, plots, and text. Notebooks `01`–`06` walk through the analysis story.

### `make report`
One command recomputes all numbers in `reports/RESULTS.md` from raw tables — prevents copy-paste errors.

### Random seed (42)
Fixed seed so random splits (CV, bootstrap) give **the same results every run**.

### Git / commit
Version control — every analysis change is traceable. `RESULTS.md` records git commit hash when generated.

### Manuscript export
`bash scripts/export_manuscript.sh` → Word (`.docx`), LaTeX (`.tex`), PDF from `manuscript_draft.md`.

---

## Part 12 — Key figures (what each picture shows)

| File | What to say when presenting |
|------|----------------------------|
| `01_demographics.png` | Who is in the cohort (age, sex, scheme) |
| `01_survival.png` | OS differs by fractionation |
| `03_tcp_curves_os_proxy.png` | Fitted TCP curves vs EQD2; S-shaped dose–response for OS proxy |
| `04_clinical_prognosis.png` | KM by scheme; WHO PS effect |
| `05_rano_vs_os_tcp_auc.png` | Same TCP model works for OS proxy, fails for RANO — side-by-side |
| `06_within_arm_rano_tcp.png` | Within each arm: dose flat, volume matters in 40 Gy |
| `07_rano_logistic_roc_40gy.png` | ROC for 40 Gy multivariable model |
| `07_rano_volume_validation_40gy.png` | Our DVH volume matches author RANO volume (ρ=1.0) |
| `08_pooled_rano_roc.png` | Pooled volume+clinical model ROC |
| `08_pyradiomics_vs_volume_auc.png` | Bar chart: radiomics vs DVH volume AUC |
| `08_pyradiomics_nested_cv_auc.png` | In-sample vs nested CV — optimism check |

---

## Part 13 — FAQ: how to explain to someone who asks

**Q: Did you prove that higher dose controls GBM better?**  
A: No. In pooled data, higher EQD2 correlates with **longer OS** but also with **more early PD on RANO** because 60 Gy patients are younger and on a curative schedule. Within each arm, dose barely varies — we cannot test dose–response properly.

**Q: What is the best predictor you found?**  
A: For early imaging response (RANO), **GTV volume** and **selected radiomics features** — after adjusting for age, PS, and treatment scheme. LOOCV/nested AUC ≈ 0.64–0.74 depending on model.

**Q: Is AUC 0.90 in the 40 Gy arm reliable?**  
A: That is **in-sample**. With only 34 patients, LOOCV drops to **0.74** — still interesting, but not clinical-grade prediction.

**Q: Why four TCP models if they all agree?**  
A: Assignment requirement and best practice — show that conclusions do not depend on one arbitrary curve shape.

**Q: Can we quote D50 = 53 Gy as GBM literature?**  
A: Only with heavy caveats: endpoint is OS median-split, not local control; pooled schemes confound dose. We compare qualitatively to literature (~40–80 Gy range) in §4k / Part VI table.

**Q: What would a proper TCP study need?**  
A: (1) Single protocol or randomised dose levels, (2) ≥1 Gy spread in GTV Dmean, (3) true LC endpoint, (4) enough events, (5) pre-specified validation cohort.

**Q: Did we need to download 52 GB of NIfTI?**  
A: Processed tables in the repo allow most analyses without local NIfTI. Full NIfTI needed to **regenerate** DVH from scratch.

---

## Part 14 — Alphabetical glossary (quick lookup)

| Term | One-line definition |
|------|---------------------|
| **AIC** | Model fit penalised for number of parameters |
| **AUC** | Discrimination ability; 0.5 = random |
| **Bootstrap** | Resampling uncertainty estimation |
| **Brier score** | Prediction accuracy for probabilities |
| **CFB-GBM** | Open French GBM radiotherapy dataset |
| **Confounding** | Hidden factor distorts exposure–outcome link |
| **Cox model** | Survival regression; outputs hazard ratios |
| **CR/PR/SD/PD/MR** | RANO response categories |
| **Cross-validation** | Train/test on different patients |
| **D50** | Dose for 50% tumour control probability |
| **Dmean** | Mean dose to GTV |
| **DVH** | Dose–volume histogram |
| **EQD2** | Dose metric normalised to 2 Gy fractions |
| **gEUD** | Single dose summarising full DVH |
| **GTV** | Gross tumour volume contour |
| **Gy (Gray)** | Unit of radiation dose |
| **HI** | Dose homogeneity index |
| **HR** | Hazard ratio (Cox model) |
| **IBSI** | Radiomics standardisation initiative |
| **KM curve** | Kaplan–Meier survival plot |
| **LOOCV** | Leave-one-out cross-validation |
| **LR test** | Likelihood ratio significance test |
| **MLE** | Maximum likelihood estimation |
| **Nested CV** | CV with feature selection inside training folds |
| **NIfTI** | Research 3D medical image format |
| **OS** | Overall survival |
| **Optimism Δ** | In-sample AUC minus CV AUC |
| **PD** | Progressive disease (RANO) |
| **Poisson TCP** | TCP model from clonogenic cell survival stats |
| **Probit TCP** | TCP model using normal CDF |
| **PyRadiomics** | Radiomics feature library / our feature table |
| **RANO** | Neuro-oncology response criteria on MRI |
| **ROC** | Receiver operating characteristic curve |
| **RTDOSE** | 3D radiotherapy dose map |
| **Spearman ρ** | Rank correlation coefficient |
| **TCP** | Tumour Control Probability |
| **TCIA** | The Cancer Imaging Archive |
| **t0 / t1** | Baseline / follow-up imaging time points |
| **TMZ** | Temozolomide chemotherapy |
| **volume_cc** | GTV volume in cubic centimetres |
| **WHO PS** | Performance status (0–4) |
| **γ50 (gamma50)** | TCP curve steepness at D50 |

---

## Part 15 — Suggested talking points for group presentation

1. **Problem:** TCP models need dose variation + proper endpoint — routine GBM data rarely provide both.
2. **Data:** 190 patients, two fractionation schemes, Version 3 RANO for 137.
3. **Methods:** DVH → four TCP models → survival → RANO logistic models → LOOCV/nested CV → radiomics benchmark.
4. **Negative result (important):** Pooled dose–TCP fails for RANO (AUC 0.43) — publishable honesty.
5. **Positive signal:** Tumour volume and radiomics predict RANO when scheme is adjusted (CV AUC ~0.64–0.74).
6. **Deliverables:** Open code, auto-generated report, manuscript draft, glossary (this document).

---

*For numeric details always cite `reports/RESULTS.md`. For clinical references see `reports/literature_table.csv` and `reports/manuscript_draft.md`.*
