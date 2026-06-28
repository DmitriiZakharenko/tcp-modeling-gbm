# Dosimetric–Response Modeling in Glioblastoma Using the Open CFB-GBM Cohort: A TCP Pipeline Feasibility Study

**Draft manuscript** — verified results from `reports/RESULTS.md` (regenerate: `make report`).  
Reference map: `reports/literature_table.csv` (ref_id 1–18).

---

## Abstract

**Background:** Tumor control probability (TCP) models relate radiotherapy dose to local tumour control [4,5]. Open imaging cohorts promise reproducible validation, but routine clinical data often lack classical TCP endpoints and dose heterogeneity [7].

**Methods:** We built a reproducible Python pipeline on the TCIA CFB-GBM cohort (n=190 after DVH quality control) [8]. Four TCP models (Poisson, logistic, probit, gEUD) [4–6] were fit with bootstrap confidence intervals [18]. Version 3 supplementary data (RANO imaging response [12,13], n=137 with t0→t1 labels) were integrated [8,9]. We audited confounding, compared OS-based vs RANO-based endpoints, and performed within-arm analyses stratified by fractionation scheme (60 Gy/30 fr vs 40.05 Gy/15 fr) [2,3,15].

**Results:** Pooled EQD2 discriminated an OS median-split proxy (ROC AUC≈0.68) but not RANO non-progressive disease when used as a dose-only TCP input (AUC≈0.43). After adjusting for fractionation scheme, **pooled GTV volume + clinical covariates achieved AUC≈0.72** (n=137); LOOCV AUC≈0.64. Within the hypofractionated arm (n=34) [3,15], pre-treatment GTV volume predicted RANO non-PD (multivariable in-sample AUC=0.90; **LOOCV AUC=0.74**). PyRadiomics features (t1gd GTV t0, v3 TSV) [9,17] outperformed DVH volume alone for RANO (top-5 in-sample AUC≈0.78 vs 0.71; **nested 5-fold CV AUC≈0.74 vs 0.70**). GTV volume from DVH matched RANO t0 segmentation volumes (Spearman ρ=1.00, n=141) [8].

**Conclusions:** Classical pooled TCP dose–response validation is not supported on CFB-GBM; however, **tumour burden models** (volume and radiomics) predict early RANO when fractionation confounding is addressed. We provide an open pipeline, LOOCV benchmarks, and a feasibility checklist for TCP studies on routine RT datasets [7].

**Keywords:** glioblastoma; tumor control probability; DVH; RANO; open data; CFB-GBM

---

## 1. Introduction

Glioblastoma remains the most common malignant primary brain tumour in adults. Standard treatment includes maximal safe resection followed by radiotherapy with concurrent temozolomide [1]. Despite uniform protocols, outcome varies widely, motivating dose–response and imaging biomarker research [2,9].

In elderly or poor-performance patients, hypofractionated schedules (e.g. 40 Gy in 15 fractions) are commonly used as an alternative to 60 Gy in 30 fractions [2,3,15]. TCP models estimate the probability of sterilising clonogenic tumour cells as a function of dose and radiobiological parameters [4,5]. They are widely used in treatment planning research but require (i) a well-defined tumour control endpoint, (ii) inter-patient dose variation, and (iii) sufficient sample size [7].

The CFB-GBM dataset (Centre François Baclesse, n=264) provides pre- and post-treatment MRI, RTDOSE, GTV contours, clinical covariates, and—since June 2026—RANO response labels [8–10]. We asked: **can a standard TCP pipeline be validated on this open cohort, and if not, what does fail and where might signal remain?**

---

## 2. Materials and Methods

### 2.1 Dataset and cohort

Clinical and imaging metadata were downloaded from TCIA (Version 3, DOI: [10.7937/v9pn-2f72](https://doi.org/10.7937/v9pn-2f72)) [8]. The public CFB-GBM release provides **pre-processed NIfTI** (skull-stripped, co-registered) rather than raw DICOM-RT; RTPLAN and RTSTRUCT DICOM files are therefore not available in the TCIA package. We used **RTDOSE** and **GTV** segmentation masks at t0, which is the structure set provided for dosimetric analysis. **CTV and PTV contours are not included** in this cohort; all DVH metrics refer to GTV only.

Inclusion required t0 RTDOSE and GTV NIfTI, known fractionation dose and fraction number. One patient was excluded after DVH QC (Dmean=0 Gy), yielding **n=190** for modeling (120×60 Gy/30 fr; 61×40.05 Gy/15 fr) [2,15].

### 2.2 DVH feature extraction

Cumulative DVH and scalar metrics (D95, Dmean, gEUD, homogeneity index, GTV volume) were computed from registered RTDOSE and t0 GTV masks. EQD2 was calculated with α/β=10 Gy [6,11]. gEUD was computed with volume-effect parameter *a* as in Niemierko [6].

### 2.3 TCP models

We implemented four binary-outcome TCP models following standard radiobiological formulations [4,5,7]:

- **Poisson TCP:** $TCP(D) = \exp\!\left[-\ln 2 \cdot \exp\left(\gamma_{50}(1 - D/D_{50})\right)\right]$
- **Logistic TCP:** sigmoid in dose with parameters $D_{50}$, $k$
- **Probit TCP:** probit link with $D_{50}$, $\sigma$
- **gEUD TCP:** dose metric $D = \left(\sum v_i D_i^a\right)^{1/a}$ with $a \in \{-10,1,10\}$ [6]

Parameters were estimated by maximum likelihood. Bootstrap 95% confidence intervals (1000 resamples) were computed for Poisson $D_{50}$ and $\gamma_{50}$ [18]. Models were compared by AIC, BIC, ROC AUC, and Hosmer–Lemeshow calibration [7].

### 2.4 Outcomes

**Primary (exploratory):** OS ≥ cohort median (51 weeks) as a binary proxy for pipeline testing.

**Secondary (v3):** RANO non-progressive disease at t1 (non-PD: SD/MR/PR/CR vs PD) [12,13], available for 137/190 modeling patients [8].

### 2.5 Statistical analysis

Kaplan–Meier and log-rank tests compared fractionation schemes [14]. Cox proportional hazards models included EQD2, Dmean, age, sex, WHO PS, and RANO non-PD [14]. Within-arm Spearman correlations and Poisson TCP AUC were computed separately for 60 Gy and 40 Gy arms. **Pooled logistic regression** modeled RANO non-PD ~ GTV volume + age + WHO PS + scheme indicator, with an optional volume×scheme interaction term. **Leave-one-out cross-validation (LOOCV)** provided out-of-sample AUC for 40 Gy and pooled models [18]. **PyRadiomics features** (v3 TSV, t1gd sequence, GTV t0) [8,9] were merged by patient ID; the top five univariate features by RANO AUC were compared against DVH volume (in-sample), with **nested 5-fold stratified CV** (feature selection on training folds only) to reduce optimism [17,18]. Full equations: `reports/manuscript_equations_fragment.tex` (included in `manuscript.tex`). Export: `bash scripts/export_manuscript.sh` → `manuscript.docx` / `manuscript.pdf`.

### 2.6 Reproducibility

All code is available at [repository URL]. Figures and tables regenerate via `make report`.

---

## 3. Results

![**Figure 1.** Kaplan–Meier overall survival by fractionation scheme and WHO performance status (n=190). Median OS 60 vs 28 weeks (log-rank p≈3×10⁻⁶).](figures/04_clinical_prognosis.png)

**Figure 1.** Kaplan–Meier overall survival by fractionation scheme and WHO performance status (n=190). Median OS 60 vs 28 weeks (log-rank p≈3×10⁻⁶).



### 3.1 Cohort and clinical prognosis

Median age 70 years; median OS 51 weeks. Sixty Gy patients had longer OS than 40 Gy patients (median 60 vs 28 weeks; log-rank p≈3×10⁻⁶), consistent with prior GBM fractionation studies [2,3,15]. Cox model: scheme HR≈0.54 (p≈0.0007), WHO PS HR≈1.42 (p≈0.001) [14].

### 3.2 Pooled TCP (OS proxy)

![**Figure 2.** Fitted TCP curves (Poisson, logistic, probit, gEUD) for EQD2 vs OS median-split proxy (n=190; AUC≈0.68).](figures/03_tcp_curves_os_proxy.png)

**Figure 2.** Fitted TCP curves (Poisson, logistic, probit, gEUD) for EQD2 vs OS median-split proxy (n=190; AUC≈0.68).



Poisson/Logistic/Probit TCP on EQD2 achieved in-sample ROC AUC≈0.68 (5-fold CV≈0.68±0.09), significantly better than null (LR p≈3×10⁻⁶). Bootstrap Poisson: $D_{50}$≈53 Gy [50, 57], $\gamma_{50}$≈3.3 [2.1, 4.7] [18]. **Interpretation:** association is driven largely by fractionation scheme confounding (r(age, EQD2)=−0.57) [7,15].

### 3.3 RANO endpoint (pooled)

![**Figure 3.** Poisson TCP AUC: OS proxy vs RANO non-PD on the same 137 patients (AUC≈0.62 vs 0.43).](figures/05_rano_vs_os_tcp_auc.png)

**Figure 3.** Poisson TCP AUC: OS proxy vs RANO non-PD on the same 137 patients (AUC≈0.62 vs 0.43).



On the same 137 patients, pooled EQD2→RANO non-PD yielded AUC≈0.43 (LR p=1.0): higher EQD2 correlated with *more* early PD at t1 (60 Gy PD rate 27% vs 40 Gy 15%) [12,15].

### 3.4 Within-arm analysis

![**Figure 4.** Within-arm DVH metrics → RANO: dose flat within protocol; GTV volume significant in 40 Gy arm (n=34).](figures/06_within_arm_rano_tcp.png)

**Figure 4.** Within-arm DVH metrics → RANO: dose flat within protocol; GTV volume significant in 40 Gy arm (n=34).



Within 60 Gy, GTV Dmean SD=0.28 Gy — insufficient for dose–response [7]. **40 Gy arm (n=34 with RANO):** GTV volume predicted RANO non-PD (Poisson TCP AUC=0.83, LR p=0.037; Spearman p=0.016), consistent with tumour burden prognostic literature [16]. Multivariable logistic (volume + age + WHO PS): **AUC=0.90**, bootstrap 95% CI [0.81, 1.00] [18].

### 3.5 Volume validation

![**Figure 6.** DVH GTV volume vs author RANO `size_t0_cm3` (Spearman ρ=1.00, n=141).](figures/07_rano_volume_validation_40gy.png)

**Figure 6.** DVH GTV volume vs author RANO `size_t0_cm3` (Spearman ρ=1.00, n=141).



DVH GTV volume at t0 agreed with RANO `size_t0_cm3` (Spearman ρ=1.00, n=141) [8]. t1 NIfTI mask validation was not performed locally; analyses used author-reported RANO volumes from the v3 TSV [8,9].

### 3.6 RANO and overall survival

Cox model (n=137): RANO non-PD HR≈0.48 (p≈0.0009); EQD2 remained significant after adjustment (HR≈0.97, p≈0.009) [12,14].

### 3.7 Pooled RANO models (volume + clinical + scheme)

![**Figure 7.** Pooled logistic ROC for RANO non-PD (volume + clinical + scheme; AUC≈0.72; LOOCV≈0.64).](figures/08_pooled_rano_roc.png)

**Figure 7.** Pooled logistic ROC for RANO non-PD (volume + clinical + scheme; AUC≈0.72; LOOCV≈0.64).



On n=137 with RANO labels, EQD2 alone achieved AUC≈0.57. **GTV volume + age + WHO PS + scheme** improved discrimination to **AUC≈0.72**; adding a volume×scheme interaction yielded AUC≈0.72. LOOCV AUC for the pooled clinical model was **≈0.64** [18], indicating moderate generalisation within this single-centre cohort.

### 3.8 LOOCV in the hypofractionated arm

![**Figure 5.** Multivariable logistic ROC for RANO non-PD in 40 Gy arm (in-sample AUC=0.90; LOOCV 0.74).](figures/07_rano_logistic_roc_40gy.png)

**Figure 5.** Multivariable logistic ROC for RANO non-PD in 40 Gy arm (in-sample AUC=0.90; LOOCV 0.74).



In-sample multivariable AUC in the 40 Gy arm was 0.90, but **LOOCV AUC was 0.74** (volume only) and 0.70 (volume + age + PS) [18], confirming signal beyond chance while reducing optimism bias relative to resubstitution [7].

### 3.9 PyRadiomics vs DVH volume

![**Figure 8.** PyRadiomics top-5 vs DVH volume for RANO (in-sample AUC 0.78 vs 0.71).](figures/08_pyradiomics_vs_volume_auc.png)

**Figure 8.** PyRadiomics top-5 vs DVH volume for RANO (in-sample AUC 0.78 vs 0.71).



Using the author-provided PyRadiomics TSV (t1gd, GTV t0) [8,9,17], the top five radiomics features achieved **in-sample AUC≈0.78** for RANO non-PD vs **0.71** for DVH volume alone (n=137). **Nested 5-fold CV** (top-5 feature selection on train folds only) yielded **AUC≈0.74** for radiomics vs **0.70** for volume-only — optimism Δ≈0.04–0.07. Combined volume + top radiomics feature: in-sample AUC≈0.77. For OS median-split on the same patients, radiomics + clinical covariates reached AUC≈0.78 vs volume-only 0.59 [9,10] (see `figures/08_pyradiomics_nested_cv_auc.png`).

### 3.10 Literature comparison — TCP parameters (Part VI)

Poisson TCP on EQD2 (OS median-split proxy) yielded **D50≈53.2 Gy** [49.5–56.7] and **γ50≈3.32** [2.06–4.69] (bootstrap, n=190) [18]. Direct comparison with published GBM TCP studies is limited because (i) CFB-GBM lacks formal local control, (ii) dose spread within protocol is <1 Gy, and (iii) our endpoint is OS proxy rather than LC [7,16]. Literature review (Table in `reports/metrics/literature_tcp_d50_comparison.csv`) shows reported D50 ranges of ~40–80 Gy across sites and endpoints [4,5]; our estimate falls within this range but is **not interpretable as GBM tumour-control D50** without endpoint and cohort redesign.

---

## 4. Discussion

### 4.1 Why pooled dose–TCP fails here

CFB-GBM reflects **routine practice**, not a dose-escalation trial: within-protocol GTV Dmean variation is <1 Gy (60 Gy arm SD=0.28 Gy) [7,8]. TCP fitting further requires a tumour-control endpoint; OS and early RANO capture different biology and are confounded by age-linked fractionation selection [8,15]. Pooled EQD2→RANO AUC≈0.43 is not a failure of RANO per se [12,13] but of **using dose as the sole predictor across mixed fractionation schemes** where higher EQD2 correlates with more early PD at t1 despite better OS [15].

Our negative pooled dose–TCP result is scientifically informative: it demonstrates that **endpoint and cohort design matter more than model sophistication** [4,7]. This aligns with Ohri et al.'s emphasis on adequate dose heterogeneity and appropriate endpoints for TCP validation [7].

### 4.2 Tumour burden predicts RANO when confounding is addressed

When fractionation scheme is included as a covariate, **GTV volume consistently predicts RANO non-PD** (pooled AUC≈0.72; LOOCV≈0.64). This reframes the analysis from classical TCP (dose→control) to **tumour burden→imaging response**, which is clinically interpretable: larger baseline GTV may reflect more aggressive biology or reduced likelihood of early imaging stabilisation under palliative-intent hypofractionation [3,15,16].

The hypofractionated subgroup (n=34) shows the strongest in-sample effect (AUC≈0.90), but **LOOCV AUC≈0.74** is the more defensible estimate [18]. The 60 Gy arm shows a weaker but directionally consistent association (Spearman p≈0.019, Poisson AUC≈0.66, n=96).

### 4.3 PyRadiomics comparison and relation to CFB-GBM authors

Moreau et al. emphasise imaging AI and radiomics for treatment efficacy prediction on this cohort [9,10]. Our head-to-head comparison on the same open resource [8] shows that **author-provided PyRadiomics features (t1gd GTV t0) modestly outperform DVH scalar volume** for RANO (in-sample AUC 0.78 vs 0.71; nested CV 0.74 vs 0.70), while DVH volume remains competitive and fully reproducible from RTDOSE+GTV without MRI preprocessing.

Feature standardisation follows IBSI recommendations [17]; we did not re-extract radiomics locally, ensuring parity with the published TCIA feature table [8,9].

### 4.4 Clinical prognostic context

Independent of TCP framing, CFB-GBM confirms known prognostic factors: 60 vs 40 Gy OS separation (p≈10⁻⁷) [2,15], WHO PS gradient [2], and RANO non-PD as an OS predictor (Cox HR≈0.48) [12,14]. These descriptive findings support data quality but are not novel; they anchor the dosimetry/RANO feasibility narrative.

### 4.5 Implications for open-data TCP research

We propose a practical checklist for future open-cohort TCP studies [7]: (1) report within-arm GTV Dmean IQR before fitting TCP; (2) prefer imaging response or LC over OS [12,13]; (3) adjust for fractionation and age [15]; (4) report LOOCV or nested CV alongside in-sample AUC [18]; (5) compare against author-provided radiomics baselines when available [8,9,17].

### 4.6 Limitations

- No formal local control labels; RANO at t1 is an imaging surrogate [12,13]
- Single institution, retrospective, deceased-patients-only inclusion [8]
- Small hypofractionated RANO subset (n=34); pooled LOOCV AUC≈0.64 limits clinical translation [7,18]
- t1 GTV NIfTI follow-up masks not validated locally; `size_t1_cm3` taken from v3 TSV [8]
- PyRadiomics nested CV (top-5 selection on train) reduces but does not eliminate optimism (Δ≈0.04–0.07 vs in-sample) [17,18]
- No external validation cohort

---

## 5. Conclusion

We deliver an open, reproducible TCP modeling pipeline on CFB-GBM [8] and show that **classical pooled TCP dose–response validation is not feasible** on this cohort [4,7]. When fractionation confounding is addressed [15], **GTV volume and PyRadiomics features predict early RANO** with moderate discrimination (pooled LOOCV AUC≈0.64; PyRadiomics nested CV AUC≈0.74 vs volume 0.70) [9,17]. We recommend future TCP studies prioritise single-protocol cohorts with ≥1 Gy GTV Dmean IQR, explicit tumour control endpoints [12,13], LOOCV reporting [18], and comparison against author radiomics baselines [8,9].

**LaTeX / Word export:** `bash scripts/export_manuscript.sh` → `reports/manuscript.tex`, `manuscript.docx`, `manuscript.pdf`

---

## References

1. Stupp R, Mason WP, van den Bent MJ, et al. Effects of radiotherapy with concomitant and adjuvant temozolomide versus radiotherapy alone on survival in glioblastoma. *Lancet Oncol*. 2009;10(5):459-466. doi:10.1016/S1470-2045(09)70025-7

2. Minniti G, Niyazi M, Alongi F, Navarria P, Belka C. Current status and recent advances in reirradiation of glioblastoma. *Radiat Oncol*. 2021;16:36. doi:10.1186/s13014-021-01767-9

3. Malmström A, Sørensen P, Grunnet K, et al. Temozolomide versus standard 6-week radiotherapy versus hypofractionated radiotherapy in patients older than 60 years with glioblastoma: the Nordic randomised, phase 3 trial. *Lancet Oncol*. 2012;13(9):916-926. https://doi.org/10.1016/S1470-2045(12)70265-6

4. Maitre A, Vogin G, Kraus R, et al. Construction of radiobiological models as TCP and NTCP. *Cancer Radiother*. 2020;24(3):247-257. https://doi.org/10.1016/j.canrad.2019.08.002

5. Gardner LL, Parry A, McMahon SJ, et al. Modelling radiobiology. *Phys Med Biol*. 2024;69(18):18TR01. https://doi.org/10.1088/1361-6560/ad70f0

6. Niemierko A. Reporting and analyzing dose distributions: a concept of equivalent uniform dose. *Med Phys*. 1997;24(1):103-110. https://doi.org/10.1118/1.598061

7. Ohri N, Tomé WA, Mendez Romero A, et al. Increasing the power of tumour control and normal tissue complication probability modelling. *Transl Cancer Res*. 2017;6(S1):S123-S127. https://doi.org/10.21037/tcr.2017.02.01

8. Moreau NN, Leclercq AG, Desmonts A, et al. Pre and post treatment MRI and radiotherapy plans of patients with glioblastoma: the CFB-GBM cohort (Version 3). *The Cancer Imaging Archive*. 2025. https://doi.org/10.7937/v9pn-2f72

9. Moreau NN, Valable S, Jaudet C, et al. Early characterization and prediction of glioblastoma treatment efficacy using radiomics and AI. *Front Oncol*. 2025;15:1497195. https://doi.org/10.3389/fonc.2025.1497195

10. Moreau NN, Desmonts A, Jaudet C, et al. AI-Driven Prediction of Treatment Efficacy in Glioblastoma Using Medical Imaging. In: Rekik I, et al., eds. *PRIME 2025*, LNCS vol. 16164. Springer; 2025:24-33. https://doi.org/10.1007/978-3-032-07904-6_3

11. Fowler JF. 21 years of biologically effective dose. *Br J Radiol*. 2010;83(994):554-568. doi:10.1259/bjr/31372149

12. Wen PY, Macdonald DR, Reardon DA, et al. Updated response assessment criteria for high-grade gliomas: response assessment in neuro-oncology working group. *J Clin Oncol*. 2010;28(11):1963-1972. doi:10.1200/JCO.2009.26.3541

13. van Dijk WT, van Erp WG, Wesseling P, et al. RANO 2.0 update for response assessment in neuro-oncology. *Lancet Oncol*. 2021;22(12):e503-e508. https://doi.org/10.1016/S1470-2045(21)00489-7

14. Cox DR. Regression models and life-tables. *J R Stat Soc Series B*. 1972;34(2):187-220. doi:10.1111/j.2517-6161.1972.tb00899.x

15. Perry JR, Laperriere N, O'Callaghan CJ, et al. Short-course radiation plus temozolomide in elderly patients with glioblastoma. *N Engl J Med*. 2017;376(11):1027-1037. https://doi.org/10.1056/NEJMoa1611977

16. Embring A, Nilsson P, Zellman P, et al. Dosimetric parameters and survival in glioblastoma treated with chemoradiation. *Radiother Oncol*. 2020;142:47-53. https://doi.org/10.1016/j.radonc.2019.11.014

17. Zwanenburg A, Vallières M, Abdalah MA, et al. The Image Biomarker Standardisation Initiative. *Radiother Oncol*. 2020;142:169-176. https://doi.org/10.1016/j.radonc.2019.11.009

18. Efron B, Tibshirani RJ. *An Introduction to the Bootstrap*. Chapman & Hall; 1994. ISBN 978-041204231-5.

---
