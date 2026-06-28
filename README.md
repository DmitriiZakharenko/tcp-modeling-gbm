# TCP Modeling in Glioblastoma

Reproducible Python workflow for Tumor Control Probability (TCP) modeling in glioblastoma using the open-access CFB-GBM radiotherapy dataset from TCIA.

The pipeline covers DICOM-free NIfTI data curation, DVH feature extraction, Poisson / Logistic / Probit / EUD-based TCP model fitting with bootstrap confidence intervals, model comparison (AIC, BIC, ROC), and survival analysis (Kaplan-Meier, Cox regression).

**Dataset:** Moreau et al. 2025, *The Cancer Imaging Archive* — [doi:10.7937/V9PN-2F72](https://doi.org/10.7937/V9PN-2F72)  
**Cohort:** 190 / 264 patients for modeling after DVH QC (191 with RTDOSE + GTV; 1 excluded for invalid RTDOSE)  
**Task board:** [GitHub Issues #1](https://github.com/DmitriiZakharenko/tcp-modeling-gbm/issues/1)

---

## Project Structure

```
tcp-modeling-gbm/
├── data/
│   ├── raw/          # NIfTI files (RTDOSE + GTV) — not versioned
│   └── processed/    # Cohort CSV, feature tables, clinical TSVs
├── notebooks/        # Analysis notebooks (run in order: 01 → 05)
├── src/
│   ├── config.py     # Central paths and constants
│   ├── data/         # Cohort builder, NIfTI loader, DVH calculator, downloader
│   ├── models/       # TCP model classes, bootstrap CI, model comparison
│   └── utils/        # Plotting helpers
├── figures/          # Exported figures (300 dpi PNG + PDF)
├── reports/          # Report sections and literature table
└── requirements.txt
```

## Setup

```bash
git clone https://github.com/DmitriiZakharenko/tcp-modeling-gbm.git
cd tcp-modeling-gbm
pip install -r requirements.txt
```

## Data

### Included in the repository (no NIfTI required)

| File | Description |
|------|-------------|
| `data/processed/cohort.csv` | 264 patients; 190 modeling-eligible after DVH QC |
| `data/processed/modeling_table.csv` | Included cohort + clinical fields + DVH metrics (190×33) |
| `data/processed/CFB-GBM_*.tsv` | Source clinical / treatment tables from TCIA |

`modeling_table.csv` columns: `patient_id`, `rt_dose_gy`, `n_fractions`, `eqd2_gy`, `survival_weeks`, `age`, `sex`, `who_status`, DVH metrics (`D2_gy` … `D98_gy`, `Vx_pct`, `gEUD_*`, `HI_gy`, `volume_cc`).

```python
import pandas as pd
df = pd.read_csv("data/processed/modeling_table.csv")
```

### Local-only (large; not in git)

| Path | Size | Role |
|------|------|------|
| `data/raw/` | ~52 GB | RTDOSE + GTV NIfTI per patient |
| `data/processed/dvh_curves/`, `dose_slices/`, `dvh_curves_all.npz` | ~130 MB | Full DVH curves and axial slices (optional; regenerate from raw) |

### Regenerate processed tables from NIfTI

```bash
make process
# verify-rt → feature_builder → export modeling_table.csv
```

Requires `data/raw/{patient_id}/t0/*_t0_rtdose.nii.gz` and `*_t0_gtv.nii.gz` for all modeling-eligible patients.

### Download NIfTI from TCIA

When direct `ascp` (port 33001) is blocked, use IBM Aspera Connect:

```bash
python -m src.data.download_rt_connect
python -m src.data.organize_raw_data --import-from /path/to/downloads --move  # if needed
```

Fallback: `python -m src.data.download_rt_files` or manual download from [CFB-GBM on TCIA](https://www.cancerimagingarchive.net/collection/cfb-gbm/) (RTDOSE + GTV NIfTI only, ~70–80 GB).

Full pipeline:

```bash
make data              # TSVs + cohort + download + verify
python -m src.data.setup_data --skip-rt --features --workers 4
```

## Usage

```bash
# Build cohort table
python -m src.data.cohort_builder

# Run notebooks in order
jupyter notebook notebooks/
```

## Reproducibility

All notebooks are designed to run top-to-bottom (`Kernel → Restart & Run All`) without errors.  
Random seed: `42`. Python 3.9+.

## Citation

```
Moreau, N. N., et al. (2025). Pre and post treatment MRI and radiotherapy plans
of patients with glioblastoma: the CFB-GBM cohort. The Cancer Imaging Archive.
https://doi.org/10.7937/V9PN-2F72
```

## License

Code: MIT. Dataset: CC BY 4.0 (see TCIA data usage policy).
