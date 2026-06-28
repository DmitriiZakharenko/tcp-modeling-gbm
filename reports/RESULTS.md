# Current Results (auto-generated)

**Last updated:** 2026-06-28  
**Git commit:** `5c14e42`  
**Regenerate:** `python -m src.reporting.update_results`

> **Outcome caveat:** Primary TCP models still use an exploratory OS median-split proxy.
> **CFB-GBM v3** adds RANO imaging response (non-PD vs PD at t1) — see §4b below.
> RANO ≠ formal local control but is closer to tumor response than OS alone.

---

## 1. Cohort (verified from `modeling_table.csv`)

| Metric | Value |
|---|---|
| Total patients (TSV) | **264** |
| Included (cohort rules) | **190** |
| Modeling table rows | **190** |
| Excluded | **74** |
| 60 Gy / 30 fr | **120** |
| 40.05 Gy / 15 fr | **61** |
| Median age | **70** yr |
| Median OS | **51** wk |
| Sex — male | **117** |
| Sex — female | **73** |
| Median GTV Dmean | **59.84** Gy |
| Median GTV volume | **36.2** cc |

## 2. Survival by fractionation (descriptive)

| Scheme | n | Median OS (wk) | IQR (wk) |
|---|---:|---:|---:|
| 60Gy_30fr | 120 | 60 | 43–79 |
| 40Gy_15fr | 61 | 28 | 16–51 |

Mann–Whitney U (60 Gy vs 40.05 Gy): **p = 1.13e-07**

## 3. Dose–outcome association

| Pair | r | p |
|---|---:|---:|
| EQD2_vs_outcome_proxy | 0.3607 | 3.19e-07 |
| Dmean_vs_outcome_proxy | 0.3604 | 3.26e-07 |
| EQD2_vs_OS | 0.3389 | 1.72e-06 |
| Dmean_vs_OS | 0.3387 | 1.75e-06 |

## 4. TCP model quality metrics

### poisson_tcp — dose = `eqd2_gy`

Outcome proxy: OS >= median (51 wk)

| Metric | Value |
|---|---|
| D50 (Gy) | 53.200 |
| gamma50 | 3.3237 |
| NLL (model / null) | 118.91 / 131.69 |
| AIC (model / null) | 241.82 / 263.37 |
| BIC (model / null) | 248.31 / 263.37 |
| LR test p-value | 2.82e-06 |
| McFadden pseudo-R² | 0.0970 |
| ROC AUC (in-sample) | 0.6835 |
| ROC AUC (5-fold CV) | 0.6848 ± 0.0978 |
| Brier (model / null) | 0.2176 / 0.2500 |
| Accuracy @ TCP=0.5 | 0.6737 |
| Sensitivity / Specificity | 0.823 / 0.521 |

### logistic_tcp — dose = `eqd2_gy`

Outcome proxy: OS >= median (51 wk)

| Metric | Value |
|---|---|
| D50 (Gy) | 53.968 |
| k | 0.0929 |
| NLL (model / null) | 118.87 / 131.69 |
| AIC (model / null) | 241.75 / 263.37 |
| BIC (model / null) | 248.24 / 263.37 |
| LR test p-value | 2.72e-06 |
| McFadden pseudo-R² | 0.0973 |
| ROC AUC (in-sample) | 0.6835 |
| ROC AUC (5-fold CV) | 0.6848 ± 0.0978 |
| Brier (model / null) | 0.2175 / 0.2500 |
| Accuracy @ TCP=0.5 | 0.6737 |
| Sensitivity / Specificity | 0.823 / 0.521 |

### probit_tcp — dose = `eqd2_gy`

Outcome proxy: OS >= median (51 wk)

| Metric | Value |
|---|---|
| D50 (Gy) | 53.903 |
| sigma_gy | 17.4479 |
| NLL (model / null) | 118.87 / 131.69 |
| AIC (model / null) | 241.73 / 263.37 |
| BIC (model / null) | 248.23 / 263.37 |
| LR test p-value | 2.70e-06 |
| McFadden pseudo-R² | 0.0974 |
| ROC AUC (in-sample) | 0.6835 |
| ROC AUC (5-fold CV) | 0.6848 ± 0.0978 |
| Brier (model / null) | 0.2175 / 0.2500 |
| Accuracy @ TCP=0.5 | 0.6737 |
| Sensitivity / Specificity | 0.823 / 0.521 |

### eud_tcp — dose = `gEUD_a-10`

Outcome proxy: OS >= median (51 wk)

| Metric | Value |
|---|---|
| D50 (Gy) | 52.348 |
| gamma50 | 2.9076 |
| geud_a | -10.0 |
| NLL (model / null) | 118.85 / 131.69 |
| AIC (model / null) | 243.69 / 263.37 |
| BIC (model / null) | 253.43 / 263.37 |
| LR test p-value | 1.11e-05 |
| McFadden pseudo-R² | 0.0975 |
| ROC AUC (in-sample) | 0.6850 |
| ROC AUC (5-fold CV) | 0.6818 ± 0.0860 |
| Brier (model / null) | 0.2176 / 0.2500 |
| Accuracy @ TCP=0.5 | 0.6684 |
| Sensitivity / Specificity | 0.823 / 0.511 |

### poisson_tcp — dose = `Dmean_gy`

Outcome proxy: OS >= median (51 wk)

| Metric | Value |
|---|---|
| D50 (Gy) | 52.394 |
| gamma50 | 2.9042 |
| NLL (model / null) | 118.87 / 131.69 |
| AIC (model / null) | 241.73 / 263.37 |
| BIC (model / null) | 248.23 / 263.37 |
| LR test p-value | 2.71e-06 |
| McFadden pseudo-R² | 0.0974 |
| ROC AUC (in-sample) | 0.6864 |
| ROC AUC (5-fold CV) | 0.6857 ± 0.0856 |
| Brier (model / null) | 0.2176 / 0.2500 |
| Accuracy @ TCP=0.5 | 0.6684 |
| Sensitivity / Specificity | 0.823 / 0.511 |


## 4b. RANO endpoint vs OS proxy (same patients, EQD2 models)

Subset with RANO t0→t1 labels: **n = 137**

| Model | Endpoint | Event rate | AUC (in-sample) | AUC (5-fold CV) | LR p |
|---|---|---:|---:|---:|---:|
| poisson_tcp | OS median | 0.52 | 0.6222 | 0.6213 ± 0.0840 | 1.02e-02 |
| poisson_tcp | RANO non-PD | 0.77 | 0.4257 | 0.4221 ± 0.0531 | 1.00e+00 |
| logistic_tcp | OS median | 0.52 | 0.6222 | 0.6213 ± 0.0840 | 1.04e-02 |
| logistic_tcp | RANO non-PD | 0.77 | 0.4257 | 0.4221 ± 0.0531 | 1.00e+00 |
| probit_tcp | OS median | 0.52 | 0.6222 | 0.6213 ± 0.0840 | 1.04e-02 |
| probit_tcp | RANO non-PD | 0.77 | 0.4257 | 0.4221 ± 0.0531 | 1.00e+00 |

Dose–outcome on RANO subset:

| Pair | r | p |
|---|---:|---:|
| EQD2_vs_os_median_proxy | 0.2555 | 0.0026 |
| EQD2_vs_rano_non_pd | -0.1273 | 0.1384 |

RANO categories (modeling subset):

- Stable Disease (SD): 69
- Progressive Disease (PD): 32
- Minor Response (MR): 22
- Partial Response (PR): 13
- Complete Response (CR): 1


## 4c. Within-arm DVH → RANO (Poisson TCP + Spearman)

| Scheme | n | Metric | std | Spearman ρ | p | Poisson AUC | LR p |
|---|---:|---|---:|---:|---:|---:|---:|
| 60Gy_30fr | 96 | Dmean_gy | 0.287 | -0.041 | 0.6945 | 0.474 | 1.000 |
| 60Gy_30fr | 96 | D95_gy | 0.615 | 0.068 | 0.5124 | 0.544 | 0.973 |
| 60Gy_30fr | 96 | volume_cc | 27.999 | 0.239 | 0.0188 | 0.655 | 0.100 |
| 60Gy_30fr | 96 | gEUD_a10_gy | 0.284 | -0.037 | 0.7188 | 0.476 | 1.000 |
| 60Gy_30fr | 96 | HI_gy (low variance) | 0.020 | -0.025 | 0.8061 | — | — |
| 40Gy_15fr | 34 | Dmean_gy | 0.091 | 0.013 | 0.9432 | 0.510 | 0.973 |
| 40Gy_15fr | 34 | D95_gy | 0.244 | -0.190 | 0.2806 | 0.345 | 1.000 |
| 40Gy_15fr | 34 | volume_cc | 40.205 | 0.411 | 0.0159 | 0.834 | 0.037 |
| 40Gy_15fr | 34 | gEUD_a10_gy | 0.094 | 0.072 | 0.6859 | 0.559 | 0.949 |
| 40Gy_15fr | 34 | HI_gy (low variance) | 0.015 | 0.300 | 0.0842 | — | — |


## 4d. Cox OS ~ age + sex + WHO PS + EQD2 + RANO non-PD (n=137)

**Concordance index:** 0.6669

| Covariate | HR | p |
|---|---:|---:|
| age | 1.013 | 0.2244 |
| sex_M | 1.014 | 0.9365 |
| who_status | 1.574 | 0.0005 |
| eqd2_gy | 0.965 | 0.0093 |
| rano_controlled_t1 | 0.485 | 0.0009 |


## 4e. Multivariable logistic — RANO non-PD (40 Gy arm)

| Model | n | AUC | Brier |
|---|---:|---:|---:|
| intercept_only | 34 | — | 0.125 |
| volume_only | 34 | 0.834 | 0.106 |
| volume_age_ps | 34 | 0.903 | 0.093 |

Bootstrap AUC (volume+age+PS): **0.927** [0.814, 1.000] (n=995 resamples)


## 4f. Pooled logistic — RANO non-PD (n≈137, volume + clinical + scheme)

| Model | n | AUC | Brier |
|---|---:|---:|---:|
| eqd2_only | 137 | 0.574 | 0.176 |
| volume_only | 137 | 0.707 | 0.162 |
| volume_age_ps | 137 | 0.710 | 0.162 |
| volume_age_ps_scheme | 137 | 0.716 | 0.160 |
| volume_age_ps_scheme_interaction | 137 | 0.722 | 0.160 |


## 4g. LOOCV — 40 Gy arm (honest AUC)

| Model | n | LOOCV AUC | LOOCV Brier |
|---|---:|---:|---:|
| volume_only | 34 | 0.745 | 0.120 |
| volume_age_ps | 34 | 0.697 | 0.127 |


## 4h. LOOCV — pooled RANO models

| Model | n | LOOCV AUC | LOOCV Brier |
|---|---:|---:|---:|
| volume_age_ps_scheme | 137 | 0.642 | 0.173 |
| volume_age_ps_scheme_interaction | 137 | 0.636 | 0.174 |


## 4i. PyRadiomics vs DVH volume — RANO endpoint (t1gd GTV t0)

| Model | n | AUC | Brier |
|---|---:|---:|---:|
| dvh_volume_only | 137 | 0.707 | 0.162 |
| dvh_volume_clinical | 137 | 0.716 | 0.160 |
| pyro_mesh_volume | 137 | 0.707 | 0.162 |
| pyro_top5 | 137 | 0.783 | 0.143 |
| pyro_top5_clinical | 137 | 0.788 | 0.143 |
| dvh_volume_plus_pyro_top1 | 137 | 0.774 | 0.147 |


## 4j. PyRadiomics nested 5-fold CV — RANO (feature selection on train only)

| Model | n | In-sample AUC | Nested CV AUC | Optimism Δ |
|---|---:|---:|---:|---:|
| dvh_volume_only | 137 | 0.707 | 0.697 | 0.010 |
| dvh_volume_clinical | 137 | 0.716 | 0.674 | 0.042 |
| pyro_top5_nested | 137 | 0.783 | 0.738 | 0.045 |
| pyro_top5_clinical_nested | 137 | 0.788 | 0.715 | 0.073 |


## 4k. Literature comparison — TCP parameters (assignment Part VI)

| Source | Model | Endpoint | n | D50 (Gy) | γ50 | Comparable? |
|---|---|---|---:|---|---|---|
| This study (CFB-GBM) | Poisson TCP (EQD2) | OS ≥ median (exploratory proxy) | 190 | 53.2 | 3.32 | reference |
| Maitre et al. 2020 | Poisson / LQ TCP (review) | Various (LC, NTCP) | review | 40–80 | 1–5 | partial |
| Ohri et al. 2017 | TCP/NTCP (review) | LC preferred over OS | review | — | — | no |
| Embring et al. 2020 | DVH metrics (not TCP fit) | OS | 120 | — | — | partial |
| Gardner et al. 2024 | Radiobiology review | Multi-scale modelling | review | context-dependent | context-dependent | partial |
| Okunieff et al. 1995 | Poisson TCP (meta-analysis) | Tumour control (mixed sites) | review | ~50–70 | ~1–4 | partial |


## 5. Bootstrap 95% CI (Poisson TCP, EQD2)

| Parameter | Estimate | 95% CI | Bootstrap SD |
|---|---:|---|---:|
| D50_gy | 53.200 | [49.540, 56.749] | 1.870 |
| gamma50 | 3.324 | [2.058, 4.689] | 0.668 |



## 6. Four-model comparison (EQD2 / gEUD, sorted by AIC)

| Model | k | AIC | BIC | ROC AUC | Brier | HL p-value |
|---|---:|---:|---:|---:|---:|---:|
| probit_tcp | 2 | 241.73 | 248.23 | 0.6835 | 0.2175 | 0.1009 |
| logistic_tcp | 2 | 241.75 | 248.24 | 0.6835 | 0.2175 | 0.1008 |
| poisson_tcp | 2 | 241.82 | 248.31 | 0.6835 | 0.2176 | 0.1003 |
| eud_tcp | 3 | 243.69 | 253.43 | 0.6850 | 0.2176 | 0.7842 |



## 7. Cox PH regression (OS ~ EQD2 + Dmean + age + sex)

**Concordance index:** 0.6319

| Covariate | HR | 95% CI | p |
|---|---:|---|---:|
| eqd2_gy | 1.121 | [0.724, 1.735] | 0.6083 |
| Dmean_gy | 0.876 | [0.594, 1.292] | 0.5054 |
| age | 1.012 | [0.994, 1.031] | 0.1776 |
| sex_M | 0.949 | [0.705, 1.277] | 0.7304 |

KM log-rank (EQD2>=50.0_vs_<50.0): χ² = 21.75, p = 3.10e-06


## 8. Clinical prognosis (Cox: OS ~ age + sex + WHO PS + scheme)

**Concordance index:** 0.6560

| Covariate | HR | p |
|---|---:|---:|
| age | 1.011 | 0.2349 |
| sex_M | 1.033 | 0.8294 |
| who_status | 1.425 | 0.0009 |
| scheme_60gy | 0.541 | 0.0007 |


## 9. OS by WHO performance status

| WHO PS | n | Median OS (wk) | IQR (wk) |
|---:|---:|---:|---:|
| 0 | 34 | 59 | 46–79 |
| 1 | 101 | 54 | 32–77 |
| 2 | 47 | 29 | 17–55 |
| 3 | 8 | 28 | 19–40 |

Kruskal–Wallis: **p = 1.58e-04**


## 10. Within-arm DVH vs OS (Spearman)

| Scheme | n | Metric | ρ | p |
|---|---:|---|---:|---:|
| 60Gy_30fr | 120 | Dmean_gy | -0.091 | 0.3204 |
| 60Gy_30fr | 120 | D95_gy | 0.028 | 0.7634 |
| 60Gy_30fr | 120 | volume_cc | -0.123 | 0.1817 |
| 60Gy_30fr | 120 | HI_gy | -0.086 | 0.3499 |
| 60Gy_30fr | 120 | gEUD_a10_gy | -0.100 | 0.2756 |
| 40Gy_15fr | 61 | Dmean_gy | 0.172 | 0.1839 |
| 40Gy_15fr | 61 | D95_gy | -0.135 | 0.2994 |
| 40Gy_15fr | 61 | volume_cc | -0.265 | 0.0392 |
| 40Gy_15fr | 61 | HI_gy | 0.164 | 0.2071 |
| 40Gy_15fr | 61 | gEUD_a10_gy | 0.209 | 0.1058 |


## 11. Hypofractionated arm: Cox OS ~ volume + covariates

n = 61, C-index = 0.641

GTV volume HR = **1.0073**/cc, p = **0.0446** (exploratory)


## 12. TCP feasibility: dose heterogeneity within arm

| Scheme | n | Dmean SD (Gy) | min–max (Gy) | DVH-TCP feasible? |
|---|---:|---:|---|---|
| 60.00 Gy / 30 fr | 120 | 0.28 | 59.3–61.8 | no |
| 40.05 Gy / 15 fr | 61 | 0.08 | 39.9–40.4 | no |


## 13. Confounding correlations

| Pair | r | p |
|---|---:|---:|
| eqd2_vs_os | 0.3389 | 1.72e-06 |
| scheme60_vs_os | 0.3378 | 1.87e-06 |
| age_vs_os | -0.2799 | 9.16e-05 |
| age_vs_eqd2 | -0.5670 | 1.48e-17 |
| age_vs_dmean | -0.5645 | 2.20e-17 |


## 14. Unused TSV fields (association with OS)

| Field | In modeling table? | ρ vs OS | p |
|---|---|---:|---:|
| rt_delay_wk | yes | 0.259 | 3.11e-04 |
| bmi | yes | 0.117 | 1.09e-01 |
| mri_t1_weeks | yes | 0.394 | 9.60e-07 |


## 15. TCP feasibility verdict

- **endpoint:** OS (weeks to death) always available. CFB-GBM v3 adds RANO response: 137/190 modeling patients with t0→t1 label. RANO is imaging response (non-PD vs PD), not formal local control.
- **dose heterogeneity:** Within 60 Gy/30 fr, GTV Dmean SD = 0.28 Gy (threshold for DVH-TCP = 1.0 Gy). DVH-based TCP within standard arm is underpowered.
- **pooled tcp:** Pooled EQD2–TCP on OS proxy is confounded: r(age, EQD2) ≈ −0.57; scheme and age drive OS. Compare §4b: RANO endpoint on same patients.
- **recommendation:** RANO v3 enables tumor-response endpoint (137/190 with t0→t1). Pooled EQD2–RANO AUC ≈ 0.43 (worse than OS proxy 0.62 on same patients) because 60 Gy has higher PD rate than 40 Gy at t1 despite better OS. Within-arm dose-TCP still limited by Dmean homogeneity; exploratory signal: GTV volume vs RANO in 40 Gy arm only.

## 16. Figures

- [`figures/01_demographics.png`](../figures/01_demographics.png)
- [`figures/01_exclusion_reasons.png`](../figures/01_exclusion_reasons.png)
- [`figures/01_fractionation.png`](../figures/01_fractionation.png)
- [`figures/01_survival.png`](../figures/01_survival.png)
- [`figures/02_dmean_by_fractionation.png`](../figures/02_dmean_by_fractionation.png)
- [`figures/02_dvh_cohort_mean.png`](../figures/02_dvh_cohort_mean.png)
- [`figures/02_dvh_correlation.png`](../figures/02_dvh_correlation.png)
- [`figures/02_dvh_overlay_median.png`](../figures/02_dvh_overlay_median.png)
- [`figures/02_dvh_overlay_sample.png`](../figures/02_dvh_overlay_sample.png)
- [`figures/02_volume_hi.png`](../figures/02_volume_hi.png)
- [`figures/03_cox_forest.png`](../figures/03_cox_forest.png)
- [`figures/03_kaplan_meier_eqd2.png`](../figures/03_kaplan_meier_eqd2.png)
- [`figures/03_model_calibration.png`](../figures/03_model_calibration.png)
- [`figures/03_tcp_curves_os_proxy.png`](../figures/03_tcp_curves_os_proxy.png)
- [`figures/04_clinical_prognosis.png`](../figures/04_clinical_prognosis.png)
- [`figures/05_rano_vs_os_tcp_auc.png`](../figures/05_rano_vs_os_tcp_auc.png)
- [`figures/06_within_arm_rano_tcp.png`](../figures/06_within_arm_rano_tcp.png)
- [`figures/07_rano_logistic_roc_40gy.png`](../figures/07_rano_logistic_roc_40gy.png)
- [`figures/07_rano_volume_validation_40gy.png`](../figures/07_rano_volume_validation_40gy.png)
- [`figures/08_pooled_rano_roc.png`](../figures/08_pooled_rano_roc.png)
- [`figures/08_pyradiomics_nested_cv_auc.png`](../figures/08_pyradiomics_nested_cv_auc.png)
- [`figures/08_pyradiomics_vs_volume_auc.png`](../figures/08_pyradiomics_vs_volume_auc.png)

---

## Interpretation snapshot

| Question | Current answer |
|---|---|
| Data pipeline complete? | Yes — 190-patient modeling table with 21 DVH metrics |
| Strong clinical signal? | Yes — 60 vs 40 Gy OS (p ≈ 3×10⁻⁶); WHO PS (p ≈ 1.6×10⁻⁴); Cox scheme HR≈0.54 |
| DVH-TCP within standard arm? | No — 60 Gy Dmean SD = 0.28 Gy (below 1 Gy threshold) |
| TCP model beats null? | Yes — LR p ≈ 3×10⁻⁶ (Poisson, EQD2) but confounded by scheme/age |
| Good discrimination (AUC ≥ 0.7)? | No — in-sample AUC ≈ 0.68, CV ≈ 0.68 ± 0.10 |
| True TCP validation? | Partial — RANO non-PD available (v3); still not formal LC |
| RANO improves AUC vs OS on same n? | No — pooled RANO AUC ≈ 0.43 vs OS ≈ 0.62 (n=137) |
| Within-arm volume → RANO (40 Gy)? | Yes — Poisson AUC ≈ 0.83, LR p ≈ 0.037 (n=34); LOOCV AUC ≈ 0.74 |
| Pooled volume + scheme → RANO? | Yes — in-sample AUC ≈ 0.72 (n=137); LOOCV ≈ 0.64 |
| PyRadiomics beats DVH volume for RANO? | Exploratory — top-5 in-sample AUC ≈ 0.78 vs volume 0.71; nested CV see §4j |
| Within-arm volume → RANO (60 Gy)? | Exploratory — AUC ≈ 0.66, Spearman p ≈ 0.019 (n=96) |
| Calibration fixes ranking? | No — Platt scaling does not change AUC on same data |

