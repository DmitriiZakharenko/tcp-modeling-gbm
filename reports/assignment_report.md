# Final Project Report: Tumor Control Probability Modeling in Glioblastoma

**Course:** TCP Modeling in GBM — Final Project  
**Dataset:** CFB-GBM (TCIA v3, DOI [10.7937/v9pn-2f72](https://doi.org/10.7937/v9pn-2f72))  
**Modeling cohort:** n=190 (after DVH QC) · **RANO subset:** n=137  
**Report date:** 2026-06-28 · Numbers from [`RESULTS.md`](RESULTS.md) (`make report`)

---

## Abstract

We implemented a reproducible Python pipeline for Tumor Control Probability (TCP) modeling on the open CFB-GBM glioblastoma cohort. Four TCP models (Poisson, logistic, probit, gEUD) were fitted with bootstrap confidence intervals. Version 3 supplementary data provided RANO imaging response labels and author PyRadiomics features. Pooled EQD2 discriminated an OS median-split proxy (AUC≈0.68) but not RANO non-progressive disease (AUC≈0.43) when used as a dose-only predictor. After adjusting for fractionation scheme, GTV volume and radiomics features predicted early RANO with moderate cross-validated discrimination (LOOCV AUC≈0.64–0.74). Classical pooled TCP dose–response validation is not supported on this routine-care dataset; we provide a feasibility checklist for future open-cohort TCP studies.

**Keywords:** glioblastoma; TCP; DVH; RANO; CFB-GBM; open data

---

## 1. Introduction and Clinical Background

Glioblastoma (GBM) is the most common malignant primary brain tumour in adults. Standard treatment includes maximal safe resection, radiotherapy, and concurrent temozolomide (Stupp protocol: 60 Gy in 30 fractions) [1]. In elderly or poor-performance patients, hypofractionated schedules (e.g. 40 Gy in 15 fractions) are commonly used [2,3,15].

**Tumor Control Probability (TCP)** models estimate the probability of sterilising clonogenic tumour cells as a function of dose and radiobiological parameters [4,5]. Valid TCP studies require (i) a well-defined tumour control endpoint, (ii) inter-patient dose variation, and (iii) adequate sample size [7].

The CFB-GBM dataset (Centre François Baclesse, n=264) provides MRI, RTDOSE, GTV contours, clinical covariates, and RANO response labels [8–10]. **Project aim:** build a standard TCP pipeline on this open cohort and assess whether classical dose–response validation is feasible; if not, identify where prognostic signal remains.

---

## 2. Dataset Description

| Item | Value |
|------|-------|
| Source | TCIA CFB-GBM Version 3 |
| Total patients (TSV) | 264 |
| Included for modeling | **190** |
| Excluded | 74 (unknown dose, missing imaging, etc.) |
| DVH QC exclusion | 1 (patient 32, Dmean=0 Gy) |
| Fractionation | 120 × 60 Gy/30 fr; 61 × 40.05 Gy/15 fr |
| Median age | 70 years |
| Median OS | 51 weeks |
| RANO t0→t1 labels | 137 / 190 |

**Imaging format:** Pre-processed **NIfTI** (skull-stripped, co-registered), not raw DICOM-RT. RTPLAN/RTSTRUCT DICOM are not in the public package. **GTV-only** contours at t0; CTV/PTV not provided.

---

## 3. Data Acquisition and Curation

1. Downloaded clinical/treatment/RANO/PyRadiomics TSVs (`python -m src.data.download_clinical_data`).
2. Built `cohort.csv` with inclusion rules (`cohort_builder.py`).
3. Downloaded RTDOSE + GTV t0 NIfTI for 191 patients (IBM Aspera).
4. Verified completeness (`verify_raw_data.py`); extracted DVH metrics (`feature_builder.py`).
5. Merged into `modeling_table.csv` (190 × 58 columns).

EQD2 computed with α/β = 10 Gy for GTV [6,11].

---

## 4. Dosimetric Feature Extraction

From registered RTDOSE and t0 GTV masks we computed:

- Cumulative **DVH** and scalars: D95, D98, D50, D2, Dmean, Dmax, GTV **volume_cc**, Vx, **gEUD** (a = −10, 1, 10), homogeneity index (HI).
- **EQD2** for cross-scheme comparison.

Within the 60 Gy arm, GTV Dmean SD = **0.28 Gy** — essentially no dose escalation within protocol.

---

## 5. Mathematical Description of TCP Models

Four binary-outcome TCP formulations [4–7]:

- **Poisson TCP:** clonogenic cell survival; parameters D50, γ50.
- **Logistic TCP:** sigmoid dose–response; D50, k.
- **Probit TCP:** normal CDF link; D50, σ.
- **gEUD TCP:** dose metric from full DVH, then Poisson-style TCP; extra parameter a.

Full equations: `reports/manuscript_equations_fragment.tex`.

**Estimation:** maximum likelihood; **95% CI:** bootstrap (1000 resamples, seed 42).

---

## 6. Parameter Estimation and Model Comparison

| Model | AIC | ROC AUC | Brier | Notes |
|-------|----:|--------:|------:|-------|
| Probit TCP | 241.73 | 0.684 | 0.218 | Best AIC (marginal) |
| Logistic TCP | 241.75 | 0.684 | 0.218 | ≈ Poisson |
| Poisson TCP | 241.82 | 0.684 | 0.218 | D50≈53 Gy, γ50≈3.3 |
| gEUD TCP | 243.69 | 0.685 | 0.218 | 3 parameters |

Bootstrap Poisson (EQD2, OS proxy): **D50 = 53.2 Gy** [49.5–56.7]; **γ50 = 3.32** [2.06–4.69].

---

## 7. Results

### 7.1 Clinical prognosis

![Figure 1](figures/04_clinical_prognosis.png)

**Figure 1.** Kaplan–Meier OS by fractionation scheme and WHO performance status. Sixty Gy patients had median OS 60 weeks vs 28 weeks for 40 Gy (log-rank p≈3×10⁻⁶). Cox model: scheme HR≈0.54 (p≈0.0007), WHO PS HR≈1.42 (p≈0.001).

### 7.2 Pooled TCP — OS proxy

![Figure 2](figures/03_tcp_curves_os_proxy.png)

**Figure 2.** TCP dose–response curves for EQD2 vs OS ≥ median (51 weeks, n=190). AUC≈0.68; LR p≈3×10⁻⁶. Interpretation confounded by fractionation (r(age, EQD2)=−0.57).

### 7.3 Pooled TCP — RANO endpoint

![Figure 3](figures/05_rano_vs_os_tcp_auc.png)

**Figure 3.** Same Poisson TCP on 137 RANO-labelled patients: OS proxy AUC≈0.62; **RANO non-PD AUC≈0.43** (no signal). Higher EQD2 associated with more early PD at t1 (27% vs 15% by scheme).

### 7.4 Within-arm analysis

![Figure 4](figures/06_within_arm_rano_tcp.png)

**Figure 4.** Stratified by fractionation: dose metrics flat within arm; **GTV volume** predicts RANO in 40 Gy arm (n=34; Poisson AUC=0.83, p≈0.037).

### 7.5 Multivariable and pooled RANO models

![Figure 5](figures/07_rano_logistic_roc_40gy.png)

**Figure 5.** 40 Gy arm multivariable logistic ROC (volume + age + PS; in-sample AUC=0.90).

![Figure 6](figures/07_rano_volume_validation_40gy.png)

**Figure 6.** DVH volume validation vs author RANO volumes (ρ=1.00, n=141).

![Figure 7](figures/08_pooled_rano_roc.png)

**Figure 7.** Pooled models (n=137): volume + clinical + scheme AUC≈0.72; LOOCV≈0.64.

### 7.6 PyRadiomics comparison

![Figure 8](figures/08_pyradiomics_vs_volume_auc.png)

**Figure 8.** Author PyRadiomics (t1gd GTV t0) top-5 features vs DVH volume (in-sample AUC 0.78 vs 0.71).

![Figure 9](figures/08_pyradiomics_nested_cv_auc.png)

**Figure 9.** Nested 5-fold CV reduces optimism: radiomics nested AUC≈0.74 vs volume 0.70.

---

## 8. Survival Analysis

Kaplan–Meier and Cox proportional hazards models (lifelines) included EQD2, Dmean, age, sex, WHO PS, fractionation scheme, and RANO non-PD.

| Covariate (Cox, n=137) | HR | p |
|------------------------|---:|---:|
| RANO non-PD | 0.48 | 0.0009 |
| WHO PS | 1.57 | 0.0005 |
| EQD2 | 0.97 | 0.009 |

RANO non-PD predicts OS independently of dose and performance status.

---

## 9. Literature Review and TCP Parameter Comparison (Part VI)

Literature table: `reports/literature_table.csv` (18 references). DOI verification: `reports/literature_doi_check.md`.

| Source | D50 (Gy) | γ50 | Endpoint | Comparable? |
|--------|----------|-----|----------|-------------|
| **This study (CFB-GBM)** | **53.2** | **3.32** | OS ≥ median (proxy) | reference |
| Maitre et al. 2020 (review) | 40–80 | 1–5 | LC / NTCP | partial |
| Okunieff et al. 1995 (meta) | ~50–70 | ~1–4 | Tumour control | partial |
| Ohri et al. 2017 | — | — | LC preferred | design note |

Our D50 falls in published ranges but is **not interpretable as GBM tumour-control D50** without local control endpoint and single-protocol dose spread.

---

## 10. Discussion

**Why pooled dose–TCP fails:** CFB-GBM reflects routine practice, not dose escalation. Within-arm GTV Dmean variation <1 Gy. OS and early RANO capture different biology; fractionation is confounded with age and treatment intent [7,15].

**Where signal remains:** GTV volume and PyRadiomics predict RANO when scheme is adjusted (pooled LOOCV AUC≈0.64; 40 Gy LOOCV≈0.74; PyRadiomics nested CV≈0.74). This reframes the project from classical TCP to **tumour burden → imaging response** — clinically interpretable and aligned with Moreau et al. radiomics work on the same dataset [9,10].

**Limitations:** No formal local control; single institution; small 40 Gy RANO subset (n=34); t1 GTV NIfTI not validated locally.

**Feasibility checklist for future open TCP studies [7]:** report within-arm Dmean IQR; prefer LC or RANO over OS; adjust for fractionation; report LOOCV/nested CV; compare against author radiomics baselines.

---

## 11. Conclusion

We delivered an open, reproducible TCP modeling pipeline on CFB-GBM and demonstrated that **classical pooled TCP validation is not feasible** on this cohort. With confounding addressed, **GTV volume and PyRadiomics features predict early RANO** with moderate cross-validated discrimination. All code, figures, and auto-generated metrics regenerate via `make report`.

---

## References

1. Stupp R et al. Lancet Oncol. 2009. doi:10.1016/S1470-2045(09)70025-7  
2. Minniti G et al. Expert Rev Anticancer Ther. 2021. doi:10.1080/14737140.2021.1919470  
3. Malmström A et al. Lancet Oncol. 2012. doi:10.1016/S1470-2045(12)70265-6  
4. Maitre A et al. Cancer Radiother. 2020. doi:10.1016/j.canrad.2019.08.002  
5. Gardner LL et al. Phys Med Biol. 2024. doi:10.1088/1361-6560/ad70f0  
6. Niemierko A. Med Phys. 1997. doi:10.1118/1.598061  
7. Ohri N et al. Transl Cancer Res. 2017. doi:10.21037/tcr.2017.02.01  
8. Moreau NN et al. TCIA. 2025. doi:10.7937/v9pn-2f72  
9. Moreau NN et al. Front Oncol. 2025. doi:10.3389/fonc.2025.1497195  
10. Moreau NN et al. PRIME 2025 LNCS. doi:10.1007/978-3-032-07904-6_3  
11. Fowler JF. Br J Radiol. 2010. doi:10.1259/bjr/86359321  
12. Wen PY et al. J Clin Oncol. 2010. doi:10.1200/JCO.2009.25.5874  
13. van Dijk WT et al. Lancet Oncol. 2021. doi:10.1016/S1470-2045(21)00489-7  
14. Cox DR. J R Stat Soc B. 1972. doi:10.1111/j.2517-6161.1972.tb00999.x  
15. Perry JR et al. N Engl J Med. 2017. doi:10.1056/NEJMoa1611977  
16. Embring A et al. Radiother Oncol. 2020. doi:10.1016/j.radonc.2019.11.014  
17. Zwanenburg A et al. Radiother Oncol. 2020. doi:10.1016/j.radonc.2019.11.009  
18. Efron B, Tibshirani RJ. An Introduction to the Bootstrap. Chapman & Hall; 1994.

---

## Appendix — Reproducibility

```bash
pip install -r requirements.txt
make report                              # regenerate RESULTS.md
python scripts/verify_literature_dois.py # check reference links
bash scripts/run_notebooks_check.sh      # execute notebooks 01, 03–06
bash scripts/export_assignment_report.sh # this document → DOCX/PDF
bash scripts/export_manuscript.sh        # scientific manuscript export
```

Group glossary (plain language): [`group_glossary_guide.md`](group_glossary_guide.md)
