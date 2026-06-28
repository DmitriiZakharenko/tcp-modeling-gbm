**Task board mirror** — synced from `TASKS.md` on 2026-06-28 (`14e985e`).  
Full board: [TASKS.md](https://github.com/DmitriiZakharenko/tcp-modeling-gbm/blob/main/TASKS.md) · Results: [`reports/RESULTS.md`](https://github.com/DmitriiZakharenko/tcp-modeling-gbm/blob/main/reports/RESULTS.md)

---

## Status summary

| Area | Progress |
|---|---|
| Data pipeline (v3 RANO + PyRadiomics) | ✅ Done |
| TCP models + survival | ✅ Done |
| RANO / pooled / LOOCV / PyRadiomics analysis | ✅ Done |
| Manuscript draft + equations + literature table | ✅ Draft done |
| t1 GTV download + PDF/slides | ⏳ Pending |

**Cohort:** 190 modeling patients · **RANO t0→t1:** 137 · **40 Gy + RANO:** 34

---

## Coding

### Infrastructure & data ✅
- [x] Repository structure, config, requirements, README
- [x] Clinical TSVs v3 download + cohort builder + RANO loader
- [x] RT NIfTI t0 download + verify (191 patients)
- [x] DVH extraction, QC, modeling_table (190×58)
- [x] `make report` → RESULTS.md + metrics CSVs

### TCP models ✅
- [x] Poisson / Logistic / Probit / EUD TCP + bootstrap CI
- [x] Model comparison (AIC, BIC, ROC, calibration)
- [x] Survival analysis (KM, Cox)
- [x] Notebooks 01–06 (03–06 via `scripts/build_analysis_notebooks.py`)

### RANO & v3 analysis ✅
- [x] `confounding_audit.py`, `stratified_analysis.py`
- [x] `rano_tcp_comparison.py` — OS vs RANO AUC
- [x] `within_arm_rano_tcp.py` — per-scheme DVH→RANO
- [x] `rano_multivariable.py` — 40 Gy logistic + volume validation
- [x] `rano_prediction_suite.py` — pooled models, LOOCV, PyRadiomics vs DVH
- [x] `validate_rano_volumes.py` — t0 NIfTI vs RANO TSV (t1 pending download)

### Remaining coding ⏳
- [ ] Download t1 GTV NIfTI (`make download-t1-gtv` — requires Aspera TCP/UDP 33001)
- [ ] All notebooks 01–06: Restart & Run All without errors
- [ ] `pip freeze > requirements.txt`
- [ ] Nested CV for PyRadiomics feature selection
- [ ] Export figures as PDF (PNG done, 22 files)

---

## Literature Review

- [x] PubMed search + `reports/literature_table.csv` (12 papers, PubMed IDs, 2010–2025)
- [x] CFB-GBM descriptor (Moreau 2025, PMID 39949753)
- [x] Outcome definition audit (OS vs LC vs RANO) — in manuscript §4.1
- [x] Literature comparison + critical discussion — manuscript §4 draft
- [ ] Final Vancouver-formatted reference list in Word/LaTeX

---

## Report Writing

- [x] `reports/manuscript_draft.md` — Abstract, Methods, Results, Discussion, Conclusion
- [x] `reports/manuscript_equations.tex` — TCP + logistic + Cox equations
- [x] Presentation outline (15 slides) in manuscript draft
- [ ] Merge to Word/LaTeX (10–15 pages)
- [ ] Export PDF
- [ ] Build slide deck

---

## Key results (2026-06-28)

| Analysis | AUC / metric |
|---|---|
| Pooled EQD2 → OS proxy (n=190) | 0.68 |
| Pooled EQD2 → RANO (n=137) | 0.43 |
| Pooled volume + clinical + scheme → RANO | 0.72 (LOOCV 0.64) |
| 40 Gy volume → RANO (n=34) | 0.90 in-sample; LOOCV 0.74 |
| PyRadiomics top-5 → RANO | 0.78 vs DVH volume 0.71 |
| Cox RANO → OS | HR 0.48, p≈0.0009 |

**Verdict:** Classical pooled TCP dose–response not validated; tumor burden (volume/radiomics) predicts RANO when scheme confounding is addressed.
