# Literature DOI / PubMed verification

**Generated:** 2026-06-28 18:01 UTC  
**Source:** `reports/literature_table.csv`  
**Command:** `python scripts/verify_literature_dois.py`

## Summary

| Check | OK | Skipped | Failed |
|-------|---:|--------:|-------:|
| DOI | 17 | 1 | 0 |
| PubMed | 13 | 5 | 0 |

## Detail

| ref | Author | Year | DOI | DOI status | PMID | PMID status |
|----:|--------|------|-----|------------|------|-------------|
| 1 | Stupp | 2009 | 10.1016/S1470-2045(09)70025-7 | ok | 19318585 | ok |
| 2 | Minniti | 2021 | 10.1186/s13014-021-01767-9 | ok | 33602305 | ok |
| 3 | Malmstrom | 2012 | 10.1016/S1470-2045(12)70265-6 | ok | 22877848 | ok |
| 4 | Maitre | 2020 | 10.1016/j.canrad.2019.08.002 | ok | 32192557 | ok |
| 5 | Gardner | 2024 | 10.1088/1361-6560/ad70f0 | ok | — | skipped |
| 6 | Niemierko | 1997 | 10.1118/1.598061 | ok | 9055874 | ok |
| 7 | Ohri | 2017 | 10.21037/tcr.2017.02.01 | ok | 28443297 | ok |
| 8 | Moreau | 2025 | 10.7937/v9pn-2f72 | ok | — | skipped |
| 9 | Moreau | 2025 | 10.3389/fonc.2025.1497195 | ok | 39949753 | ok |
| 10 | Moreau | 2025 | 10.1007/978-3-032-07904-6_3 | ok | — | skipped |
| 11 | Fowler | 2010 | 10.1259/bjr/31372149 | ok | 20603408 | ok |
| 12 | Wen | 2010 | 10.1200/JCO.2009.26.3541 | ok | 20231676 | ok |
| 13 | van Dijk | 2021 | 10.1016/S1470-2045(21)00489-7 | ok | 34774159 | ok |
| 14 | Cox | 1972 | 10.1111/j.2517-6161.1972.tb00899.x | ok | — | skipped |
| 15 | Perry | 2017 | 10.1056/NEJMoa1611977 | ok | 28296600 | ok |
| 16 | Embring | 2020 | 10.1016/j.radonc.2019.11.014 | ok | 31836236 | ok |
| 17 | Zwanenburg | 2020 | 10.1016/j.radonc.2019.11.009 | ok | 31912224 | ok |
| 18 | Efron | 1994 | — | skipped | — | skipped |

## Notes

- **Skipped** — book/chapter entries without DOI or PubMed (e.g. Efron bootstrap, ISBN-only).
- **403 on DOI** — publisher anti-bot; open link in a browser or use PubMed ID.
- All resolvable DOI/PMID links returned HTTP success at check time.
