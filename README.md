# TCP Modeling in Glioblastoma

Reproducible Python workflow for Tumor Control Probability (TCP) modeling in glioblastoma using the open-access CFB-GBM radiotherapy dataset from TCIA.

The pipeline covers DICOM-free NIfTI data curation, DVH feature extraction, Poisson / Logistic / Probit / EUD-based TCP model fitting with bootstrap confidence intervals, model comparison (AIC, BIC, ROC), and survival analysis (Kaplan-Meier, Cox regression).

**Dataset:** Moreau et al. 2025, *The Cancer Imaging Archive* — [doi:10.7937/V9PN-2F72](https://doi.org/10.7937/V9PN-2F72)  
**Cohort:** 191 / 264 patients with RTDOSE + GTV segmentation and known treatment dose  
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

Clinical and treatment TSVs (~50 KB, no login required):

```bash
cd data/processed
curl -O https://www.cancerimagingarchive.net/wp-content/uploads/CFB-GBM_clinical_data_v02_20260129.tsv
curl -O https://www.cancerimagingarchive.net/wp-content/uploads/CFB-GBM_treatment_data_v02_20260129.tsv
curl -O https://www.cancerimagingarchive.net/wp-content/uploads/CFB-GBM_treatment_imaging_availability_v02_20260129.tsv
curl -O https://www.cancerimagingarchive.net/wp-content/uploads/CFB-GBM_columns_description_new_v02_20260129.tsv
```

RTDOSE and GTV NIfTI files (~1.5 GB, requires IBM Aspera Connect):

```bash
python -m src.data.download_rt_files
```

Full 208 GB archive (all MRI sequences) is **not required** for this project.

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
