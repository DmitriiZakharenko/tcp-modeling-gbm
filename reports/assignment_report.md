# Final Project Report: Tumor Control Probability Modeling in Glioblastoma

**Course:** TCP Modeling in GBM — Final Project  
**Dataset:** CFB-GBM (TCIA Version 3, DOI 10.7937/v9pn-2f72)  
**Institution:** Centre François Baclesse (open cohort via The Cancer Imaging Archive)  
**Modeling cohort:** n = 190 (after DVH quality control)  
**RANO imaging subset:** n = 137 (t0 to t1 labels)  
**Report date:** 2026-06-28

---

## Abstract

We implemented a reproducible Python pipeline for Tumor Control Probability (TCP) modeling on the open-access CFB-GBM glioblastoma radiotherapy cohort. The workflow covers clinical data curation, NIfTI-based DVH feature extraction, four TCP model families (Poisson, logistic, probit, gEUD-based), bootstrap confidence intervals, survival analysis, and Version 3 RANO imaging endpoints with PyRadiomics comparison.

Pooled EQD2 discriminated an exploratory overall survival (OS) median-split proxy (ROC AUC approximately 0.68) but not RANO non-progressive disease when dose was the sole predictor (AUC approximately 0.43). After adjusting for fractionation scheme and clinical covariates, pre-treatment GTV volume and author-provided PyRadiomics features predicted early RANO with moderate cross-validated discrimination (pooled LOOCV AUC approximately 0.64; hypofractionated arm LOOCV approximately 0.74; PyRadiomics nested CV approximately 0.74).

**Conclusion:** classical pooled TCP dose-response validation is not supported on routine CFB-GBM data because of limited within-protocol dose spread and endpoint confounding. The project delivers a complete open pipeline, formal negative result, and tumour-burden benchmarks for future TCP feasibility studies.

**Keywords:** glioblastoma; tumor control probability; DVH; EQD2; RANO; PyRadiomics; CFB-GBM; open data

---

## 1. Introduction and Clinical Background

### 1.1 Glioblastoma and standard radiotherapy

Glioblastoma (GBM) is the most common malignant primary brain tumour in adults. Maximal safe resection followed by radiotherapy with concurrent temozolomide remains the standard of care for fit patients (Stupp protocol: 60 Gy delivered in 30 fractions of 2 Gy) [1]. Despite a uniform nominal protocol, observed outcomes vary widely because of age, performance status, molecular factors, and tumour burden at baseline [2].

In elderly or frail patients, shorter hypofractionated schedules are widely used. The CFB-GBM cohort includes 40.05 Gy in 15 fractions, consistent with pragmatic GBM management described in Nordic and short-course trials [3,15]. These two fractionation schemes differ not only in total dose and fraction size but also in patient selection: in our cohort, 60 Gy patients were younger (median age 65 vs 75 years) and had longer OS (median 60 vs 28 weeks).

### 1.2 Tumor Control Probability (TCP) modeling

TCP models estimate the probability that all clonogenic tumour cells are sterilised after radiotherapy, as a function of delivered dose and radiobiological parameters such as D50 (dose for 50% control probability) and gamma50 (steepness of the dose-response curve) [4,5]. TCP models are used in treatment-planning research, protocol comparison, and hypothetical dose escalation studies.

Valid TCP validation requires three conditions emphasized in the radiotherapy modeling literature [7]:

1. **A tumour-control endpoint** (local control or a validated imaging surrogate), not a distant-survival proxy alone.
2. **Meaningful inter-patient dose heterogeneity** within the structure of interest (typically GTV Dmean spread greater than approximately 1 Gy in TCP audit workflows).
3. **Adequate sample size** and pre-specified cross-validation when predictive models are reported.

### 1.3 CFB-GBM and project aims

The CFB-GBM dataset (Centre François Baclesse, n = 264 in the public TSV release) provides pre- and post-treatment MRI, RTDOSE maps, GTV segmentations, clinical tables, and—since Version 3—RANO response labels and PyRadiomics feature tables [8–10]. The dataset is distributed as pre-processed NIfTI rather than full DICOM-RT plans.

**Primary aim:** build a reproducible TCP modeling pipeline on this open cohort, estimate D50 and gamma50 with uncertainty, compare models, and relate our parameters to published TCP literature.

**Secondary aim (Version 3 extension):** determine whether early RANO imaging response and PyRadiomics features rescue dose-response signal when classical pooled TCP fails, and report honest out-of-sample performance (LOOCV and nested cross-validation).

---

## 2. Dataset Description and Cohort Construction

### 2.1 Source data and imaging format

| Item | Value |
|:-----|:------|
| Data repository | The Cancer Imaging Archive (TCIA) |
| Collection | CFB-GBM Version 3 |
| DOI | 10.7937/v9pn-2f72 |
| Total patients in clinical TSV | 264 |
| Patients with RTDOSE flag | 194 |
| Patients with GTV segmentation | 191 |
| **Included for TCP modeling** | **190** |
| Primary exclusion | Unknown RT dose (70 patients) |
| DVH QC exclusion | 1 (patient 32, GTV Dmean = 0 Gy) |
| Fractionation scheme A | 120 patients: 60 Gy / 30 fractions |
| Fractionation scheme B | 61 patients: 40.05 Gy / 15 fractions |
| Median age (modeling cohort) | 70 years |
| Median overall survival | 51 weeks |
| Sex (M / F) | 117 / 73 |
| Median GTV volume | 36.2 cc |
| Median GTV Dmean | 59.84 Gy |
| RANO t0 to t1 labels available | 137 / 190 (72%) |
| RANO labels in 40 Gy arm | 34 patients |

**Structure set note:** all dosimetric metrics refer to **GTV at t0** only. CTV and PTV contours are not included in the public CFB-GBM release (GTV-only segmentation). RTPLAN and RTSTRUCT DICOM objects are not distributed on TCIA Version 3; the cohort provides co-registered **RTDOSE** and **GTV NIfTI** volumes. D95, D98, D50, D2, gEUD, EQD2, and GTV volume were computed from voxel dose grids and binary masks—the same quantities typically exported from a treatment-planning system after structure delineation and dose calculation.

### 2.2 Cohort flow (264 to 190)

Starting from 264 patients in the clinical table, 70 were excluded because prescribed RT dose was unknown or missing. Among imaging-eligible patients, one patient failed DVH quality control (zero dose to GTV on RTDOSE). The final modeling table contains 190 rows and 58 columns (`modeling_table.csv`), including clinical covariates, DVH scalars, EQD2, RANO categories, and PyRadiomics merge keys. Notebook 01 (`01_cohort_overview.ipynb`) documents the full exclusion breakdown and demographic tables.

### 2.3 Outcome definitions

**Primary TCP endpoint (exploratory):** OS greater than or equal to cohort median (51 weeks), used as a binary outcome for dose–response modeling. This is not equivalent to formal tumour control.

**Secondary endpoint (Version 3):** RANO non-progressive disease at t1, defined as stable disease, minor response, partial response, or complete response versus progressive disease [12,13]. RANO non-PD rate in the modeling subset was approximately 77%.

---

## 3. Data Acquisition, Curation, and Reproducibility

### 3.1 Pipeline overview

1. Download clinical, treatment, RANO, and PyRadiomics TSV files (`python -m src.data.download_clinical_data`).
2. Merge tables and apply inclusion rules (`src/data/cohort_builder.py`); export `cohort.csv`.
3. Download RTDOSE and GTV t0 NIfTI for 191 patients (IBM Aspera Connect; approximately 52 GB local storage).
4. Verify file completeness (`verify_raw_data.py`).
5. Compute cumulative DVH and scalar metrics (`feature_builder.py`, `dvh_calculator.py`).
6. Apply DVH QC (`dvh_qc.py`); export `modeling_table.csv`.
7. Fit TCP models, survival models, and RANO prediction suite; export summary tables and figures.

All analysis code lives in `src/`. Random seed 42 is fixed for bootstrap and cross-validation splits. Notebooks 01–06 reproduce the analysis end-to-end after `pip install -r requirements.txt`.

### 3.2 Imaging format and DVH inputs

| Typical DICOM-RT workflow | CFB-GBM Version 3 | This study |
|:--------------------------|:------------------|:-----------|
| RTPLAN / RTSTRUCT | Not in public download | RTDOSE + GTV NIfTI |
| CTV / PTV structures | Not released | GTV t0 only |
| RTDOSE + structure | RTDOSE + GTV NIfTI | `dvh_calculator.py`, `feature_builder.py` |
| Scalar DVH metrics | D95, D98, D50, D2, volume, gEUD, EQD2 | Exported to `modeling_table.csv` |

Given RTDOSE and GTV masks, cumulative DVH and standard metrics match those derived from a DICOM plan export.

### 3.3 EQD2 correction

To compare the two fractionation schemes on a common biologically effective scale, we computed equivalent dose in 2 Gy fractions (EQD2) using the linear-quadratic model with alpha/beta = 10 Gy for GTV [6,11]:

$$\mathrm{EQD2} = D_{\mathrm{total}} \cdot \frac{d_{\mathrm{fx}} + \alpha/\beta}{2 + \alpha/\beta}$$

---

## 4. Dosimetric Feature Extraction

### 4.1 DVH computation

For each patient, cumulative dose-volume histograms (DVH) were computed inside the t0 GTV mask from registered RTDOSE grids. Scalar metrics exported per patient include:

| Metric | Description |
|:-------|:------------|
| Dmean | Mean dose to GTV (Gy) |
| D95, D98, D50, D2 | Dose received by at least x% of GTV volume |
| Dmax | Maximum dose to GTV |
| volume_cc | GTV volume (cm³) |
| Vx | Percentage of GTV receiving at least x Gy |
| gEUD_a | Generalised equivalent uniform dose for a = -10, 1, 10 |
| HI | Homogeneity index inside GTV |
| eqd2_gy | EQD2 as defined above |

![Sample cumulative DVH curves for three patients (GTV t0, RTDOSE).](figures/02_dvh_overlay_sample.png){ width=90% }

**Figure (DVH).** Representative cumulative dose–volume histograms inside the t0 GTV mask, illustrating inter-patient dose distribution at fixed protocol dose (60 Gy arm) versus hypofractionation (40 Gy arm). Full cohort mean DVH: `figures/02_dvh_cohort_mean.png`.

### 4.2 Dose heterogeneity audit

A central feasibility finding is that **within-arm GTV Dmean variation is negligible**. In the 60 Gy arm, GTV Dmean standard deviation was only **0.28 Gy**. Therefore, classical dose-response fitting *within* a single protocol arm lacks statistical power regardless of TCP model choice [7]. Between-arm dose differences reflect different fractionation schedules and patient selection, not planned dose escalation.

---

## 5. Mathematical Description of TCP Models

We implemented four binary-outcome TCP formulations following standard radiobiology references [4–7]:

| Model | Parameters | Formula (conceptual) |
|:------|:-------------|:--------------------|
| Poisson TCP | D50, gamma50 | Clonogenic survival with exponential shoulder |
| Logistic TCP | D50, k | Sigmoid dose-response |
| Probit TCP | D50, sigma | Normal cumulative distribution link |
| gEUD TCP | D50, gamma50, a | DVH-collapsed dose metric then Poisson-style TCP |

Full equations are documented in `reports/manuscript_equations_fragment.tex`. Parameters were estimated by **maximum likelihood** using `scipy.optimize.minimize` with method **`L-BFGS-B`** and box constraints on D50 and gamma50 (`src/models/poisson_tcp.py`). This is the standard choice for smooth, low-dimensional TCP likelihoods: fast convergence, explicit bounds, and stable gradients.

### 5.1 Maximum-likelihood optimization

Poisson TCP parameters were estimated by constrained maximum likelihood (`scipy.optimize.minimize`, method L-BFGS-B). For this smooth two-parameter likelihood, bounded quasi-Newton methods are appropriate: they converge quickly, respect clinically plausible bounds (D50 30–80 Gy, gamma50 0.5–8), and remain stable across bootstrap resamples.

We cross-checked the MLE with TNC, SLSQP, Powell, Nelder-Mead, and differential evolution on the full cohort (n = 190). All methods recovered D50 = 53.20 Gy within 0.002 Gy; Nelder-Mead and differential evolution required 2–20× more computation time without improving the fit.

Model comparison used AIC, BIC, ROC AUC, Brier score, and Hosmer-Lemeshow calibration.

### 5.2 Uncertainty quantification (bootstrap and profile likelihood)

Bootstrap 95% confidence intervals (1000 resamples, seed 42) and profile-likelihood 95% intervals were computed for Poisson D50 and gamma50 [18]. Bootstrap intervals reflect resampling uncertainty; profile likelihood uses the chi-square(1) threshold on one-dimensional profile curves, consistent with standard MLE theory.

---

## 6. Parameter Estimation and Four-Model Comparison

### 6.1 Pooled TCP on OS proxy (EQD2, n = 190)

| Model | k | AIC | BIC | ROC AUC | Brier | HL p-value |
|:------|--:|----:|----:|--------:|------:|-----------:|
| probit_tcp | 2 | 241.73 | 248.23 | 0.684 | 0.218 | 0.101 |
| logistic_tcp | 2 | 241.75 | 248.24 | 0.684 | 0.218 | 0.101 |
| poisson_tcp | 2 | 241.82 | 248.31 | 0.684 | 0.218 | 0.100 |
| eud_tcp (gEUD) | 3 | 243.69 | 253.43 | 0.685 | 0.218 | 0.784 |

All three two-parameter models are nearly identical on this cohort. The Poisson fit is reported for literature comparison.

### 6.2 Bootstrap confidence intervals (Poisson, EQD2)

| Parameter | Estimate (Gy) | 95% Bootstrap CI | Bootstrap SD |
|:----------|-------------:|:-----------------|-------------:|
| D50 | 53.20 | [49.54, 56.75] | 1.87 |
| gamma50 | 3.32 | [2.06, 4.69] | 0.67 |

### 6.3 Profile-likelihood confidence intervals (Poisson, EQD2)

| Parameter | Estimate | 95% Profile-LI CI |
|:----------|--------:|:------------------|
| D50 (Gy) | 53.20 | [48.91, 58.15] |
| gamma50 | 3.32 | [1.76, 4.85] |

Profile and bootstrap intervals are similar in width; both exclude uninformative values (gamma50 > 1).

Likelihood-ratio test vs null model: p approximately 3 x 10^-6. Five-fold cross-validation AUC: 0.68 +/- 0.10.

### 6.4 Calibration — observed vs predicted

![Calibration plot: predicted TCP (Poisson, EQD2) vs observed OS proxy rate by decile.](figures/03_model_calibration.png){ width=88% }

**Figure (calibration).** Hosmer-Lemeshow groups (10 bins): observed fraction with OS >= median vs mean predicted TCP per bin. HL p approximately 0.10; no systematic over- or under-prediction across dose deciles.

---

## 7. Results

### 7.1 Clinical prognosis and confounding

![Figure 1. Kaplan-Meier overall survival by fractionation scheme and WHO performance status.](figures/04_clinical_prognosis.png){ width=92% }

**Figure 1.** Kaplan-Meier OS by fractionation scheme and WHO performance status (n = 190). Sixty Gy patients had median OS 60 weeks vs 28 weeks for 40 Gy (log-rank p approximately 3 x 10^-6). Cox model (age + sex + WHO PS + scheme): scheme HR approximately 0.54 (p approximately 0.0007), WHO PS HR approximately 1.42 (p approximately 0.001).

| OS by WHO PS | n | Median OS (weeks) | IQR (weeks) |
|:-------------|--:|------------------:|------------:|
| PS 0 | 34 | 59 | 46–79 |
| PS 1 | 101 | 54 | 32–77 |
| PS 2 | 47 | 29 | 17–55 |

Kruskal-Wallis across PS groups: p approximately 1.6 x 10^-4. Confounding audit: Pearson r(age, EQD2) = -0.57; older patients disproportionately received hypofractionation.

### 7.2 Pooled TCP — OS proxy

![Figure 2. TCP dose-response curves for EQD2 vs OS median-split proxy.](figures/03_tcp_curves_os_proxy.png){ width=88% }

**Figure 2.** Fitted Poisson, logistic, probit, and gEUD TCP curves (EQD2 vs OS >= median, n = 190). In-sample AUC approximately 0.68. Interpretation: association partly reflects fractionation and patient selection, not within-protocol dose escalation.

### 7.3 Pooled TCP — RANO endpoint failure

![Figure 3. OS vs RANO TCP endpoint comparison on the same patients.](figures/05_rano_vs_os_tcp_auc.png){ width=88% }

**Figure 3.** Poisson TCP (EQD2) on 137 RANO-labelled patients: OS proxy AUC approximately 0.62; **RANO non-PD AUC approximately 0.43** (LR p = 1.0 vs null). Higher EQD2 correlated with more early PD at t1 (27% PD in 60 Gy arm vs 15% in 40 Gy arm) despite better long-term OS.

| Endpoint | n | Poisson AUC | 5-fold CV AUC |
|:---------|--:|------------:|--------------:|
| OS median-split | 137 | 0.622 | 0.621 |
| RANO non-PD | 137 | 0.426 | 0.422 |

### 7.4 Within-arm stratified analysis

![Figure 4. Within-arm DVH metrics predicting RANO non-PD.](figures/06_within_arm_rano_tcp.png){ width=92% }

**Figure 4.** Within-arm Spearman and Poisson TCP AUC by fractionation scheme. Sixty Gy arm (n = 96): GTV Dmean SD = 0.28 Gy; no dose-response. Forty Gy arm (n = 34): GTV volume predicted RANO (Spearman rho approximately 0.41, p approximately 0.016; Poisson AUC = 0.83, LR p approximately 0.037).

### 7.5 Hypofractionated arm — multivariable logistic regression

![Figure 5. ROC curves for 40 Gy arm RANO models.](figures/07_rano_logistic_roc_40gy.png){ width=85% }

**Figure 5.** Multivariable logistic regression in the 40.05 Gy arm (n = 34): volume-only AUC = 0.83; volume + age + WHO PS AUC = 0.90 in-sample.

| Model (40 Gy arm) | n | In-sample AUC | LOOCV AUC |
|:------------------|--:|--------------:|----------:|
| volume only | 34 | 0.834 | 0.745 |
| volume + age + PS | 34 | 0.903 | 0.697 |

LOOCV confirms signal beyond chance but reduces optimism relative to in-sample AUC = 0.90.

### 7.6 Volume validation

![Figure 6. DVH GTV volume vs author RANO t0 volume.](figures/07_rano_volume_validation_40gy.png){ width=85% }

**Figure 6.** Agreement between DVH-derived GTV volume and author-reported RANO `size_t0_cm3` (Spearman rho = 1.00, n = 141). Validates the feature extraction pipeline.

### 7.7 Pooled RANO models (volume + clinical + scheme)

![Figure 7. Pooled logistic ROC for RANO non-PD.](figures/08_pooled_rano_roc.png){ width=92% }

**Figure 7.** Pooled models on n = 137: EQD2 alone AUC approximately 0.57; GTV volume alone approximately 0.71; volume + age + WHO PS + scheme approximately 0.72.

| Pooled model | n | AUC | LOOCV AUC |
|:-------------|--:|----:|----------:|
| volume + age + PS + scheme | 137 | 0.716 | 0.642 |
| volume + age + PS + scheme + interaction | 137 | 0.722 | 0.636 |

### 7.8 PyRadiomics vs DVH volume

![Figure 8. In-sample AUC: PyRadiomics top-5 vs DVH volume.](figures/08_pyradiomics_vs_volume_auc.png){ width=88% }

**Figure 8.** Author-provided PyRadiomics features (t1gd, GTV t0, Version 3 TSV) vs DVH GTV volume for RANO non-PD: top-5 radiomics AUC approximately 0.78 vs volume approximately 0.71 (n = 137).

![Figure 9. Nested cross-validation optimism check.](figures/08_pyradiomics_nested_cv_auc.png){ width=88% }

**Figure 9.** Nested 5-fold CV (top-5 feature selection on training folds only): radiomics nested AUC approximately 0.74 vs volume approximately 0.70; optimism delta approximately 0.04–0.07.

| Model (nested CV) | In-sample AUC | Nested CV AUC | Optimism delta |
|:------------------|-------------:|--------------:|---------------:|
| DVH volume only | 0.707 | 0.697 | 0.010 |
| PyRadiomics top-5 | 0.783 | 0.738 | 0.045 |
| PyRadiomics top-5 + clinical | 0.788 | 0.715 | 0.073 |

---

## 8. Survival Analysis

Kaplan-Meier and Cox proportional hazards models (lifelines) were fitted for overall survival.

### 8.1 Cox model with RANO (n = 137)

| Covariate | Hazard ratio | p-value |
|:----------|-------------:|--------:|
| RANO non-PD at t1 | 0.48 | 0.0009 |
| WHO performance status | 1.57 | 0.0005 |
| EQD2 (per Gy) | 0.97 | 0.009 |
| Age (per year) | 1.01 | 0.22 |

Concordance index: 0.667. RANO non-PD predicts longer OS independently of EQD2 and performance status, supporting the clinical relevance of the imaging endpoint even when dose-only TCP fails.

### 8.2 Clinical prognosis model (n = 190)

| Covariate | Hazard ratio | p-value |
|:----------|-------------:|--------:|
| 60 Gy scheme (vs 40 Gy) | 0.54 | 0.0007 |
| WHO performance status | 1.43 | 0.001 |
| Age | 1.01 | 0.23 |
| Sex (male) | 1.03 | 0.83 |

---

## 9. Literature Review and TCP Parameter Comparison

Literature sources are catalogued in `reports/literature_table.csv` (18 references). Link verification: `reports/literature_doi_check.md`.

### 9.1 Comparison of D50 and gamma50

| Source | Model | Endpoint | n | D50 (Gy) | gamma50 | Comparable to this study? |
|:-------|:------|:---------|--:|---------:|--------:|:--------------------------|
| **This study (CFB-GBM)** | Poisson TCP (EQD2) | OS >= median (proxy) | 190 | **53.2** | **3.32** | reference |
| Maitre et al. 2020 | TCP review | LC / NTCP (multiple) | review | 40–80 | 1–5 | partial |
| Okunieff et al. 1995 | Poisson meta-analysis | Tumour control (mixed sites) | review | 50–70 | 1–4 | partial |
| Ohri et al. 2017 | TCP/NTCP review | LC preferred | review | — | — | design guidance |
| Embring et al. 2020 | DVH prognostic | OS (GBM) | 120 | — | — | endpoint differs |
| Gardner et al. 2024 | Radiobiology review | Multi-scale | review | context | context | partial |

Our bootstrap D50 (53 Gy) falls within published broad ranges but **must not be interpreted as a GBM tumour-control D50** because (i) the endpoint is OS proxy, not local control; (ii) dose spread within protocol is less than 1 Gy; (iii) two fractionation schemes are pooled.

### 9.2 Relation to CFB-GBM author radiomics work

Moreau et al. demonstrated imaging AI and PyRadiomics for treatment efficacy prediction on the same TCIA resource [9,10]. Our head-to-head comparison shows that author-provided t1gd GTV t0 radiomics modestly outperform scalar DVH volume for early RANO (nested CV AUC 0.74 vs 0.70), while DVH volume remains fully reproducible from RTDOSE and GTV without MRI preprocessing [17].

---

## 10. Discussion

### 10.1 Why pooled dose-TCP fails on CFB-GBM

CFB-GBM documents routine clinical practice, not a dose-escalation trial. Within the 60 Gy arm, every patient received essentially the planned protocol dose (GTV Dmean SD = 0.28 Gy). TCP models mathematically require variation in the predictor; without it, fitted D50 and gamma50 reflect cohort mixing rather than a true dose-response law [7].

Using OS >= median as a TCP endpoint further confounds interpretation: higher EQD2 patients live longer primarily because they were younger and treated with curative-intent 60 Gy, not because incremental GTV dose within protocol controlled tumour. When RANO at t1 is used instead, pooled EQD2 alone performs worse than random (AUC approximately 0.43) because early imaging PD was more frequent in the 60 Gy arm despite better OS.

This negative result is scientifically valuable: it demonstrates that **endpoint and cohort design dominate model sophistication** in open-data TCP studies.

### 10.2 Where prognostic signal remains

When fractionation scheme and clinical covariates are included, **pre-treatment GTV volume** predicts RANO non-PD with moderate discrimination (pooled AUC approximately 0.72; LOOCV approximately 0.64). In the hypofractionated subgroup, the effect is stronger but the sample is small (n = 34; LOOCV AUC approximately 0.74 for volume only).

**PyRadiomics features** from the Version 3 TSV provide a modest incremental gain over volume (nested CV AUC approximately 0.74 vs 0.70). Optimism bias between in-sample and nested CV (delta approximately 0.04–0.07) underscores the need for pre-specified validation in future work.

### 10.3 Limitations

- No formal local control labels; RANO at t1 is an imaging surrogate, not pathological control [12,13].
- Single institution, retrospective, deceased-patients-only public release [8].
- Small hypofractionated RANO subset (n = 34); wide bootstrap CIs in 40 Gy arm.
- t1 GTV NIfTI masks not validated locally; t1 volumes taken from author TSV.
- PyRadiomics features were author-provided; we did not re-extract locally (intentional parity with published baseline [9]).

### 10.4 Feasibility checklist for future open TCP studies

1. Report within-arm GTV Dmean IQR before any TCP fit [7].
2. Prefer local control or validated imaging response over OS alone [12].
3. Adjust for fractionation, age, and performance status [15].
4. Report LOOCV or nested CV alongside in-sample AUC [18].
5. Compare against author-provided radiomics baselines when available [8,9,17].

---

## 11. Conclusion

We implemented a reproducible TCP modeling workflow on the CFB-GBM cohort: NIfTI-based DVH extraction (GTV only), four TCP model families, bootstrap and profile-likelihood uncertainty intervals, model comparison (AIC, BIC, ROC, Brier), calibration assessment, survival analysis, literature comparison of D50 and gamma50, and RANO/PyRadiomics prediction with nested cross-validation.

**Classical pooled TCP dose-response validation is not feasible** on this routine-care dataset. With confounding addressed, **tumour burden metrics** (GTV volume and selected PyRadiomics features) predict early RANO with moderate out-of-sample discrimination.

---

## References

1. Stupp R, Mason WP, van den Bent MJ, Weller M, Fisher B, Taphoorn MJB, et al. Effects of radiotherapy with concomitant and adjuvant temozolomide versus radiotherapy alone on survival in glioblastoma: a randomised phase III trial. Lancet Oncol. 2009;10(5):459-466. doi:10.1016/S1470-2045(09)70025-7

2. Minniti G, Niyazi M, Alongi F, Navarria P, Belka C. Current status and recent advances in reirradiation of glioblastoma. Radiat Oncol. 2021;16:36. doi:10.1186/s13014-021-01767-9

3. Malmstrom A, Glimelius B, Marosi C, Stupp R, Frappaz D, Schultz H, et al. Temozolomide versus standard 6-week radiotherapy versus hypofractionated radiotherapy in patients older than 60 years with glioblastoma: the Nordic randomised, open-label, phase 3 trial. Lancet Oncol. 2012;13(9):916-926. doi:10.1016/S1470-2045(12)70265-6

4. Maitre A, Maitre M, Boda-Heggemann J, Lartigau E. Construction of radiobiological models TCP and NTCP. Cancer Radiother. 2020;24(6-7):564-569. doi:10.1016/j.canrad.2019.08.002

5. Gardner LL, McMahon SJ, Butterworth KT, Prise KM. Modelling radiobiology. Phys Med Biol. 2024;69(20):20TR01. doi:10.1088/1361-6560/ad70f0

6. Niemierko A. Reporting and analyzing dose distributions: a concept of equivalent uniform dose. Med Phys. 1997;24(1):103-110. doi:10.1118/1.598061

7. Ohri N, Dicker AP, Showalter TN. Increasing power of tumor control and normal tissue complication probability modeling. Transl Cancer Res. 2017;6(S1):S92-S103. doi:10.21037/tcr.2017.02.01

8. Moreau NN, Fournier L, Darrigues L, et al. Pre and post treatment MRI and radiotherapy plans of patients with glioblastoma: the CFB-GBM cohort. The Cancer Imaging Archive; 2025. doi:10.7937/v9pn-2f72

9. Moreau NN, Fournier L, Darrigues L, et al. Early characterization and prediction of glioblastoma treatment efficacy using radiomics and AI. Front Oncol. 2025;15:1497195. doi:10.3389/fonc.2025.1497195

10. Moreau NN, Fournier L, Darrigues L, et al. AI-Driven Prediction of Treatment Efficacy in Glioblastoma Using Medical Imaging. In: PRIME 2025. LNCS. Springer; 2025. doi:10.1007/978-3-032-07904-6_3

11. Fowler JF. 21 years of biologically effective dose. Br J Radiol. 2010;83(994):554-568. doi:10.1259/bjr/31372149

12. Wen PY, Macdonald DR, Reardon DA, Cloughesy TF, Sorensen AG, Galanis E, et al. Updated response assessment criteria for high-grade gliomas: response assessment in neuro-oncology working group. J Clin Oncol. 2010;28(11):1963-1972. doi:10.1200/JCO.2009.26.3541

13. van Dijk WT, van den Bent MJ, Vogelbaum MA, et al. RANO 2.0 update for response assessment in neuro-oncology. Lancet Oncol. 2021;22(12):e503-e508. doi:10.1016/S1470-2045(21)00489-7

14. Cox DR. Regression models and life-tables. J R Stat Soc Series B. 1972;34(2):187-220. doi:10.1111/j.2517-6161.1972.tb00899.x

15. Perry JR, Laperriere N, O'Callaghan CJ, Brandes AA, Menten J, Phillips C, et al. Short-course radiation plus temozolomide in elderly patients with glioblastoma. N Engl J Med. 2017;376(11):1027-1037. doi:10.1056/NEJMoa1611977

16. Embring A, Glimelius B, Söderberg J, et al. DVH parameters and survival in glioblastoma patients treated with chemoradiation. Radiother Oncol. 2020;144:177-183. doi:10.1016/j.radonc.2019.11.014

17. Zwanenburg A, Vallières M, Abdalah MA, et al. The Image Biomarker Standardisation Initiative: standardized quantitative radiomics for high-throughput image-based phenotyping. Radiother Oncol. 2020;145:75-82. doi:10.1016/j.radonc.2019.11.009

18. Efron B, Tibshirani RJ. An Introduction to the Bootstrap. New York: Chapman and Hall; 1994.

---

## Appendix A — Reproducibility commands

```text
pip install -r requirements.txt
make report
make verify-dois
make check-notebooks
make export-assignment
```

Plain-language glossary for group presentations: `reports/group_glossary_guide.md`

## Appendix B — Repository outputs

| Output | Path |
|:-------|:-----|
| Modeling table | data/processed/modeling_table.csv |
| Summary results | reports/RESULTS.md |
| This report (PDF) | reports/assignment_report.pdf |
| Scientific manuscript | reports/manuscript.pdf |
| Figure captions | reports/figure_captions.md |
| Literature table | reports/literature_table.csv |
