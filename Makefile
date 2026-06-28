.PHONY: data data-clinical cohort verify-rt features export-modeling download-rt help

help:
	@echo "Targets:"
	@echo "  make data           Full pipeline: TSVs + cohort + NIfTI + verify"
	@echo "  make data-clinical  Download clinical TSVs and build cohort only"
	@echo "  make cohort         Rebuild cohort.csv from existing TSVs"
	@echo "  make download-rt    Download RTDOSE + GTV NIfTI files (Aspera)"
	@echo "  make verify-rt      Check NIfTI completeness against cohort.csv"
	@echo "  make features       Extract DVH features (requires raw NIfTI data)"
	@echo "  make export-modeling  Merge cohort + features → modeling_table.csv"
	@echo "  make report           Regenerate reports/RESULTS.md and metrics CSVs"

data:
	python -m src.data.setup_data

data-clinical:
	python -m src.data.setup_data --skip-rt

cohort:
	python -m src.data.cohort_builder

download-rt:
	python -m src.data.download_rt_files

download-rt-connect:
	python -m src.data.download_rt_connect

verify-rt:
	python -m src.data.verify_raw_data --write-manifest

features:
	python -m src.data.feature_builder --workers 4

export-modeling:
	python -m src.data.export_modeling_dataset

process: verify-rt features export-modeling

report:
	python -m src.reporting.update_results
