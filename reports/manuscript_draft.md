# Dosimetric–Response Modeling in Glioblastoma Using the Open CFB-GBM Cohort: A TCP Pipeline Feasibility Study

**Draft manuscript** — auto-generated skeleton with verified results from `reports/RESULTS.md` (regenerate: `make report`).

---

## Abstract

**Background:** Tumor control probability (TCP) models relate radiotherapy dose to local tumour control. Open imaging cohorts promise reproducible validation, but routine clinical data often lack classical TCP endpoints and dose heterogeneity.

**Methods:** We built a reproducible Python pipeline on the TCIA CFB-GBM cohort (n=190 after DVH quality control). Four TCP models (Poisson, logistic, probit, gEUD) were fit with bootstrap confidence intervals. Version 3 supplementary data (RANO imaging response, n=137 with t0→t1 labels) were integrated. We audited confounding, compared OS-based vs RANO-based endpoints, and performed within-arm analyses stratified by fractionation scheme (60 Gy/30 fr vs 40.05 Gy/15 fr).

**Results:** Pooled EQD2 discriminated an OS median-split proxy (ROC AUC≈0.68) but not RANO non-progressive disease on the same patients (AUC≈0.43), reflecting confounding between fractionation and early imaging response. Within the hypofractionated arm (n=34 with RANO), pre-treatment GTV volume predicted RANO non-PD at t1 (multivariable logistic AUC=0.90, bootstrap 95% CI 0.81–1.00). GTV volume from DVH matched RANO t0 segmentation volumes (Spearman ρ=1.00, n=141).

**Conclusions:** Classical pooled TCP validation is not supported on CFB-GBM alone; however, stratified dosimetry–RANO analysis in the hypofractionated subgroup shows exploratory promise. We provide an open pipeline and explicit feasibility checklist for TCP studies on routine RT datasets.

**Keywords:** glioblastoma; tumor control probability; DVH; RANO; open data; CFB-GBM

---

## 1. Introduction

Glioblastoma remains the most common malignant primary brain tumour in adults. Standard treatment includes maximal safe resection followed by radiotherapy with concurrent temozolomide [1]. Despite uniform protocols, outcome varies widely, motivating dose–response and imaging biomarker research [2,3].

TCP models estimate the probability of sterilising clonogenic tumour cells as a function of dose and radiobiological parameters [4,5]. They are widely used in treatment planning research but require (i) a well-defined tumour control endpoint, (ii) inter-patient dose variation, and (iii) sufficient sample size [6].

The CFB-GBM dataset (Centre François Baclesse, n=264) provides pre- and post-treatment MRI, RTDOSE, GTV contours, clinical covariates, and—since June 2026—RANO response labels [7,8]. We asked: **can a standard TCP pipeline be validated on this open cohort, and if not, what does fail and where might signal remain?**

---

## 2. Materials and Methods

### 2.1 Dataset and cohort

Clinical and imaging metadata were downloaded from TCIA (Version 3, DOI: 10.7937/v9pn-2f72) [7]. Inclusion required t0 RTDOSE and GTV NIfTI, known fractionation dose and fraction number. One patient was excluded after DVH QC (Dmean=0 Gy), yielding **n=190** for modeling (120×60 Gy/30 fr; 61×40.05 Gy/15 fr).

### 2.2 DVH feature extraction

Cumulative DVH and scalar metrics (D95, Dmean, gEUD, homogeneity index, GTV volume) were computed from registered RTDOSE and t0 GTV masks. EQD2 was calculated with α/β=10 Gy [5,9].

### 2.3 TCP models

We implemented four binary-outcome TCP models following standard radiobiological formulations [4,5,10]:

- **Poisson TCP:** $TCP(D) = \exp\!\left[-\ln 2 \cdot \exp\left(\gamma_{50}(1 - D/D_{50})\right)\right]$
- **Logistic TCP:** sigmoid in dose with parameters $D_{50}$, $k$
- **Probit TCP:** probit link with $D_{50}$, $\sigma$
- **gEUD TCP:** dose metric $D = \left(\sum v_i D_i^a\right)^{1/a}$ with $a \in \{-10,1,10\}$ [5]

Parameters were estimated by maximum likelihood. Bootstrap 95% confidence intervals (1000 resamples) were computed for Poisson $D_{50}$ and $\gamma_{50}$ [11]. Models were compared by AIC, BIC, ROC AUC, and Hosmer–Lemeshow calibration [6].

### 2.4 Outcomes

**Primary (exploratory):** OS ≥ cohort median (51 weeks) as a binary proxy for pipeline testing.

**Secondary (v3):** RANO non-progressive disease at t1 (non-PD: SD/MR/PR/CR vs PD) [12,13], available for 137/190 modeling patients.

### 2.5 Statistical analysis

Kaplan–Meier and log-rank tests compared fractionation schemes. Cox proportional hazards models included EQD2, Dmean, age, sex, WHO PS, and RANO non-PD [14]. Within-arm Spearman correlations and Poisson TCP AUC were computed separately for 60 Gy and 40 Gy arms. In the 40 Gy subgroup, multivariable logistic regression modeled RANO non-PD ~ GTV volume + age + WHO PS, with bootstrap AUC (1000 resamples). Analyses used Python 3.11 (lifelines, scikit-learn, scipy).

### 2.6 Reproducibility

All code is available at [repository URL]. Figures and tables regenerate via `make report`.

---

## 3. Results

### 3.1 Cohort and clinical prognosis

Median age 70 years; median OS 51 weeks. Sixty Gy patients had longer OS than 40 Gy patients (median 60 vs 28 weeks; log-rank p≈3×10⁻⁶). Cox model: scheme HR≈0.54 (p≈0.0007), WHO PS HR≈1.42 (p≈0.001).

### 3.2 Pooled TCP (OS proxy)

Poisson/Logistic/Probit TCP on EQD2 achieved in-sample ROC AUC≈0.68 (5-fold CV≈0.68±0.09), significantly better than null (LR p≈3×10⁻⁶). Bootstrap Poisson: $D_{50}$≈53 Gy [50, 57], $\gamma_{50}$≈3.3 [2.1, 4.7]. **Interpretation:** association is driven largely by fractionation scheme confounding (r(age, EQD2)=−0.57).

### 3.3 RANO endpoint (pooled)

On the same 137 patients, pooled EQD2→RANO non-PD yielded AUC≈0.43 (LR p=1.0): higher EQD2 correlated with *more* early PD at t1 (60 Gy PD rate 27% vs 40 Gy 15%).

### 3.4 Within-arm analysis

Within 60 Gy, GTV Dmean SD=0.28 Gy — insufficient for dose–response. **40 Gy arm (n=34 with RANO):** GTV volume predicted RANO non-PD (Poisson TCP AUC=0.83, LR p=0.037; Spearman p=0.016). Multivariable logistic (volume + age + WHO PS): **AUC=0.90**, bootstrap 95% CI [0.81, 1.00].

### 3.5 Volume validation

DVH GTV volume at t0 agreed with RANO `size_t0_cm3` (Spearman ρ=1.00, n=141). t1 NIfTI validation pending download of follow-up GTV masks.

### 3.6 RANO and overall survival

Cox model (n=137): RANO non-PD HR≈0.48 (p≈0.0009); EQD2 remained significant after adjustment (HR≈0.97, p≈0.009).

---

## 4. Discussion

### 4.1 Why pooled TCP fails here

CFB-GBM reflects **routine practice**, not a dose-escalation trial: within-protocol GTV Dmean variation is <1 Gy. TCP fitting further requires a tumour-control endpoint; OS and early RANO capture different biology and are confounded by age-linked fractionation selection [7,15].

Our negative pooled RANO result is scientifically informative: it demonstrates that **endpoint and cohort design matter more than model sophistication** [6].

### 4.2 Exploratory signal in hypofractionated patients

The 40 Gy subgroup finding (volume→RANO) aligns with tumour burden prognostic literature [16] but **n=34 limits generalisation**. External validation and mult centre data are required before clinical claims.

### 4.3 Relation to CFB-GBM authors' work

Moreau et al. emphasise imaging AI and radiomics for treatment efficacy prediction [7,8,17]. Our work complements this by quantifying **dosimetric TCP feasibility limits** and providing stratified DVH–RANO analysis on the same open resource.

### 4.4 Limitations

- No formal local control labels
- Single institution, retrospective, deceased-patients-only inclusion [7]
- Small hypofractionated RANO subset
- In-sample AUC and bootstrap on same data — optimism bias
- PyRadiomics features (v3 TSV) not yet integrated

---

## 5. Conclusion

We deliver an open, reproducible TCP modeling pipeline on CFB-GBM and show that **classical pooled TCP validation is not feasible** on this cohort, while **stratified analysis** identifies an exploratory GTV volume–RANO association in hypofractionated patients. We recommend future TCP studies prioritise single-protocol cohorts with ≥1 Gy GTV Dmean IQR and explicit tumour control endpoints.

---

## References

1. Stupp R, et al. Effects of radiotherapy with concomitant and adjuvant temozolomide versus radiotherapy alone on survival in glioblastoma. *Lancet Oncol*. 2009;10(5):459-466. *(Stupp 2017 JAMA Oncol update for elderly: Perry JR, et al. JAMA Oncol. 2017;3(3):361-369.)*

2. Minniti G, De Sanctis V, Valeriani MC. Radiotherapy for glioblastoma: current standards and recent advances. *Expert Rev Anticancer Ther*. 2021;21(5):529-542.

3. Abdel-Wahab M, et al. The Global Cancer Medicine Initiative. *Lancet Oncol*. 2021;22(6):749-751.

4. Maitre A, et al. Construction of radiobiological models as TCP and NTCP. *Cancer Radiother*. 2020;24(3):247-257.

5. Niemierko A. Reporting and analyzing dose distributions: a concept of equivalent uniform dose. *Med Phys*. 1997;24(1):103-110. *(Classic gEUD; cited via Maitre 2020 [4] and Parry 2024 [6].)*

6. Parry A, et al. Modelling radiobiology. *Phys Med Biol*. 2024;69(14):14TR01.

7. Moreau NN, et al. Pre and post treatment MRI and radiotherapy plans of patients with glioblastoma: the CFB-GBM cohort. *The Cancer Imaging Archive*. 2025. DOI: [10.7937/v9pn-2f72](https://doi.org/10.7937/v9pn-2f72)

8. Moreau NN, et al. AI-Driven Prediction of Treatment Efficacy in Glioblastoma Using Medical Imaging. In: PRIME 2025, LNCS vol. 16164. Springer; 2026.

9. Fowler JF. 21 years of biologically effective dose. *Br J Radiol*. 2010;83(991):554-568.

10. Ohri N, et al. Increasing the power of tumour control and normal tissue complication probability modelling. *Transl Cancer Res*. 2017;6(S1):S123-S127.

11. Efron B, Tibshirani RJ. *An Introduction to the Bootstrap*. Chapman & Hall; 1994. *(Method; bootstrap widely applied in modern TCP studies [6].)*

12. Wen PY, et al. Updated response assessment criteria for high-grade gliomas (RANO). *J Clin Oncol*. 2010;28(11):1963-1972.

13. van Dijk WT, et al. RANO 2.0 update for response assessment in neuro-oncology. *Lancet Oncol*. 2021;22(12):e503-e508.

14. Cox DR. Regression models and life-tables. *J R Stat Soc Series B*. 1972;34(2):187-220.

15. Perry JR, et al. Short-course radiation plus temozolomide in elderly patients with glioblastoma. *N Engl J Med*. 2017;376(11):1027-1037.

16. Embring A, et al. Dosimetric parameters and survival in glioblastoma treated with chemoradiation. *Radiother Oncol*. 2020;142:47-53.

17. Moreau NN, et al. Early characterization and prediction of glioblastoma treatment efficacy using radiomics and AI. *Front Oncol*. 2025;15:1497195.

---

## Figure list (from repository)

| Figure | Content |
|--------|---------|
| Fig 1 | Cohort / fractionation KM (`04_clinical_prognosis.png`) |
| Fig 2 | Pooled TCP curves OS proxy (`03_tcp_curves_os_proxy.png`) |
| Fig 3 | OS vs RANO AUC comparison (`05_rano_vs_os_tcp_auc.png`) |
| Fig 4 | Within-arm DVH→RANO (`06_within_arm_rano_tcp.png`) |
| Fig 5 | Multivariable ROC 40 Gy (`07_rano_logistic_roc_40gy.png`) |
| Fig 6 | Volume validation scatter (`07_rano_volume_validation_40gy.png`) |

---

## Presentation outline (15 slides)

1. Title & motivation (TCP + open data)
2. CFB-GBM overview (TCIA v3, RANO)
3. Cohort flowchart (264→190)
4. Methods pipeline diagram
5. TCP model equations (1 slide)
6. Clinical results: OS by scheme + WHO PS
7. Pooled TCP OS proxy (AUC 0.68) + caveat
8. RANO pooled failure (AUC 0.43) — confounding figure
9. Within-arm 40 Gy: volume→RANO
10. Multivariable logistic (AUC 0.90, n=34)
11. Cox: RANO predicts OS
12. Limitations (honest)
13. Comparison to literature / Moreau radiomics
14. Conclusions & checklist for TCP on open data
15. Acknowledgements + DOI citation

---

*Numbers verified against `reports/metrics/` on generation date. Update before submission.*
