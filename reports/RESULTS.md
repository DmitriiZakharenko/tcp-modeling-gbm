# Current Results (auto-generated)

**Last updated:** 2026-06-28  
**Git commit:** `2113673`  
**Regenerate:** `python -m src.reporting.update_results`

> **Outcome caveat:** CFB-GBM provides overall survival (weeks) only.
> TCP models below use an exploratory binary proxy (OS ≥ cohort median),
> **not** true local tumour control. Interpret accordingly.

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

## 5. Figures

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

---

## Interpretation snapshot

| Question | Current answer |
|---|---|
| Data pipeline complete? | Yes — 190-patient modeling table with 21 DVH metrics |
| Strong clinical signal? | Yes — OS differs by fractionation (p ≪ 0.001) |
| TCP model beats null? | Yes — LR p ≈ 3×10⁻⁶ (Poisson, EQD2) |
| Good discrimination (AUC ≥ 0.7)? | No — in-sample AUC ≈ 0.68, CV ≈ 0.68 ± 0.10 |
| True TCP validation? | No — OS proxy only; local control unavailable |

