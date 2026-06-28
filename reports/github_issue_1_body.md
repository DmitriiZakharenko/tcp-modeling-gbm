**Task board mirror** — synced from `TASKS.md` on 2026-06-28 (`cd967bf`).  
Full board: [TASKS.md](https://github.com/DmitriiZakharenko/tcp-modeling-gbm/blob/main/TASKS.md) · Results: [`reports/RESULTS.md`](https://github.com/DmitriiZakharenko/tcp-modeling-gbm/blob/main/reports/RESULTS.md) · Group guide: [`group_glossary_guide.md`](https://github.com/DmitriiZakharenko/tcp-modeling-gbm/blob/main/reports/group_glossary_guide.md)

---

## Status summary

| Area | Progress |
|---|---|
| Data pipeline (v3 RANO + PyRadiomics) | ✅ Done |
| TCP models + survival | ✅ Done |
| RANO / pooled / LOOCV / nested CV / PyRadiomics | ✅ Done |
| Manuscript export (Word / LaTeX / PDF draft) | ✅ Done |
| Final report polish + slides + presentation prep | ⏳ In progress |
| t1 GTV download (optional) | ⏳ Pending Aspera |

**Cohort:** 190 modeling patients · **RANO t0→t1:** 137 · **40 Gy + RANO:** 34

---

## Coding — remaining

- [ ] Download t1 GTV NIfTI (`make download-t1-gtv` — optional; requires Aspera)
- [ ] All notebooks 01–06: Restart & Run All without errors

---

## Team deliverables — `[WRITE]`

- [ ] Build 15-min slide deck (outline in `manuscript_draft.md`)
- [ ] Write figure captions for slides and final report
- [ ] Proofread Introduction, Discussion, Abstract
- [ ] Assemble final report PDF (10–15 pp): `manuscript.docx` + figures

## Team deliverables — `[LIT]`

- [ ] Format references Vancouver style in `manuscript.docx`
- [ ] Verify PubMed / DOI links in `literature_table.csv` (18 refs)
- [ ] One short paragraph per key reference for slides/appendix (optional)

## All members

- [ ] Read `group_glossary_guide.md` (Parts 1–3, 13)
- [ ] Assign presentation sections and rehearse (Part 15 talking points)

---

## Key results (2026-06-28)

| Analysis | AUC / metric |
|---|---|
| Pooled EQD2 → OS proxy (n=190) | 0.68 |
| Pooled EQD2 → RANO (n=137) | 0.43 |
| Pooled volume + clinical + scheme → RANO | 0.72 (LOOCV 0.64) |
| 40 Gy volume → RANO (n=34) | 0.90 in-sample; LOOCV 0.74 |
| PyRadiomics top-5 → RANO | 0.78 in-sample; nested CV 0.74 vs volume 0.70 |
| Cox RANO → OS | HR 0.48, p≈0.0009 |

**Verdict:** Classical pooled TCP dose–response not validated; tumor burden (volume/radiomics) predicts RANO when scheme confounding is addressed.
