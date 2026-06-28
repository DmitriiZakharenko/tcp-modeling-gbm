# Figure Captions — CFB-GBM TCP Project

Use these captions in the **final report**, **slides**, and **manuscript export**.  
All numbers match [`RESULTS.md`](RESULTS.md) (regenerate: `make report`).

| Fig | File | Short title (slides) |
|-----|------|----------------------|
| 1 | `figures/04_clinical_prognosis.png` | Clinical prognosis by scheme and WHO PS |
| 2 | `figures/03_tcp_curves_os_proxy.png` | Pooled TCP dose–response (OS proxy) |
| 3 | `figures/05_rano_vs_os_tcp_auc.png` | OS vs RANO endpoint comparison |
| 4 | `figures/06_within_arm_rano_tcp.png` | Within-arm DVH metrics → RANO |
| 5 | `figures/07_rano_logistic_roc_40gy.png` | Multivariable ROC — 40 Gy arm |
| 6 | `figures/07_rano_volume_validation_40gy.png` | GTV volume validation (t0) |
| 7 | `figures/08_pooled_rano_roc.png` | Pooled RANO logistic models |
| 8 | `figures/08_pyradiomics_vs_volume_auc.png` | PyRadiomics vs DVH volume |
| 9 | `figures/08_pyradiomics_nested_cv_auc.png` | PyRadiomics nested CV optimism check |

---

## Figure 1 — Clinical prognosis

**File:** `figures/04_clinical_prognosis.png`

**Caption:** Kaplan–Meier overall survival by fractionation scheme (60 Gy/30 fractions vs 40.05 Gy/15 fractions) and by WHO performance status in the CFB-GBM modeling cohort (n=190). Median OS was 60 weeks (60 Gy arm) vs 28 weeks (40 Gy arm; log-rank p≈3×10⁻⁶). Higher WHO PS was associated with shorter survival (Kruskal–Wallis p≈1.6×10⁻⁴).

**Slide bullet:** Two fractionation schemes differ strongly in OS — confounding must be addressed in dose–response models.

---

## Figure 2 — Pooled TCP curves (OS proxy)

**File:** `figures/03_tcp_curves_os_proxy.png`

**Caption:** Fitted Poisson, logistic, probit, and gEUD-based TCP curves relating EQD2 to the exploratory OS median-split outcome proxy (OS ≥ 51 weeks, n=190). In-sample ROC AUC≈0.68 for all parametric forms; bootstrap Poisson D50≈53 Gy [50–57], γ50≈3.3. Association is partly driven by fractionation-scheme confounding (r(age, EQD2)=−0.57).

**Slide bullet:** TCP pipeline works technically, but pooled dose–response reflects mixed protocols, not within-trial dose escalation.

---

## Figure 3 — OS vs RANO TCP comparison

**File:** `figures/05_rano_vs_os_tcp_auc.png`

**Caption:** Head-to-head comparison of Poisson TCP (EQD2) discrimination for OS median-split vs RANO non-progressive disease at t1 on the same 137 patients with imaging labels. OS proxy AUC≈0.62; RANO non-PD AUC≈0.43 (LR p=1.0 vs null). Pooled EQD2 alone does not predict early imaging response across mixed fractionation schemes.

**Slide bullet:** Same dose metric succeeds on OS proxy but fails on RANO when schemes are pooled.

---

## Figure 4 — Within-arm analysis

**File:** `figures/06_within_arm_rano_tcp.png`

**Caption:** Within-arm Spearman correlations and Poisson TCP AUC for DVH metrics predicting RANO non-PD, stratified by fractionation scheme. In the 60 Gy arm (n=96), GTV Dmean SD=0.28 Gy — insufficient dose spread for dose–response. In the 40 Gy arm (n=34), pre-treatment GTV volume predicted RANO (Spearman ρ≈0.41, p≈0.016; Poisson AUC=0.83, LR p≈0.037).

**Slide bullet:** Signal appears for tumour volume in the hypofractionated arm; dose is flat within each protocol.

---

## Figure 5 — Multivariable ROC (40 Gy arm)

**File:** `figures/07_rano_logistic_roc_40gy.png`

**Caption:** Receiver operating characteristic curves for logistic regression predicting RANO non-PD in the 40.05 Gy/15-fraction arm (n=34). GTV volume alone: AUC=0.83; volume + age + WHO PS: AUC=0.90. Leave-one-out cross-validation AUC=0.74 (volume only), indicating moderate generalisation despite small sample size.

**Slide bullet:** Strong in-sample fit in n=34; LOOCV AUC≈0.74 is the honest estimate.

---

## Figure 6 — Volume validation

**File:** `figures/07_rano_volume_validation_40gy.png`

**Caption:** Agreement between GTV volume computed from DVH/RTDOSE+NIfTI masks at t0 and author-reported RANO `size_t0_cm3` from CFB-GBM Version 3 (Spearman ρ=1.00, n=141). Validates the dosimetric feature extraction pipeline against the published supplementary table.

**Slide bullet:** Our DVH volumes match the dataset authors' RANO segmentation volumes exactly.

---

## Figure 7 — Pooled RANO models

**File:** `figures/08_pooled_rano_roc.png`

**Caption:** ROC curves for pooled logistic models of RANO non-PD (n=137): EQD2 alone (AUC≈0.57), GTV volume alone (AUC≈0.71), and volume + age + WHO PS + fractionation scheme (AUC≈0.72). LOOCV AUC≈0.64 for the full clinical model.

**Slide bullet:** Adjusting for scheme and clinical covariates, tumour volume predicts RANO with moderate discrimination.

---

## Figure 8 — PyRadiomics vs DVH volume

**File:** `figures/08_pyradiomics_vs_volume_auc.png`

**Caption:** In-sample ROC AUC comparison on n=137 (t1gd GTV t0, author-provided PyRadiomics TSV): top-5 radiomics features AUC≈0.78 vs DVH GTV volume AUC≈0.71 for RANO non-PD. Combined volume + top radiomics feature: AUC≈0.77.

**Slide bullet:** Author-provided radiomics modestly outperforms scalar DVH volume in-sample.

---

## Figure 9 — PyRadiomics nested cross-validation

**File:** `figures/08_pyradiomics_nested_cv_auc.png`

**Caption:** In-sample vs nested 5-fold stratified cross-validation AUC for RANO non-PD models (n=137). Top-5 PyRadiomics features (selected on training folds only): nested CV AUC≈0.74 vs DVH volume nested CV≈0.70; optimism Δ≈0.04–0.07 relative to in-sample AUC≈0.78.

**Slide bullet:** Radiomics advantage persists out-of-sample but optimism bias is non-negligible.

---

## Supplementary figures (optional in report appendix)

| File | Suggested caption |
|------|-------------------|
| `figures/03_kaplan_meier_eqd2.png` | Kaplan–Meier OS by EQD2 dichotomy (≥50 vs <50 Gy EQD2). |
| `figures/03_model_calibration.png` | Hosmer–Lemeshow calibration of Poisson TCP (EQD2, OS proxy). |
| `figures/03_cox_forest.png` | Forest plot: Cox model OS ~ EQD2 + Dmean + age + sex. |
| `figures/01_demographics.png` | Cohort demographics and fractionation breakdown (264→190). |
| `figures/01_survival.png` | Overall survival distribution in modeling cohort. |
