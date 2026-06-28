"""Generate notebooks 03–05 with cell ids for nbformat."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NB = ROOT / "notebooks"


def cell_md(text: str, cid: str) -> dict:
    return {
        "cell_type": "markdown",
        "id": cid,
        "metadata": {},
        "source": [line + "\n" for line in text.strip().split("\n")],
    }


def cell_code(text: str, cid: str) -> dict:
    return {
        "cell_type": "code",
        "id": cid,
        "metadata": {},
        "outputs": [],
        "execution_count": None,
        "source": [line + "\n" for line in text.strip().split("\n")],
    }


def write_nb(name: str, cells: list) -> None:
    nb = {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python"},
        },
        "cells": cells,
    }
    path = NB / name
    path.write_text(json.dumps(nb, indent=1))
    print(f"Wrote {path}")


def build_03() -> None:
    write_nb(
        "03_tcp_models.ipynb",
        [
            cell_md(
                """# 03 — TCP Models (OS proxy + RANO v3)

- Pooled TCP curves (OS ≥ median)
- RANO non-PD vs OS on same patients
- Within-arm DVH → RANO (60 / 40 Gy)
- Cox OS ~ EQD2 + RANO""",
                "03-title",
            ),
            cell_code(
                """import sys
from pathlib import Path
sys.path.insert(0, str(Path.cwd().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.config import DATA_PROCESSED, FIGURES_DIR, REPORTS_DIR
from src.models.poisson_tcp import PoissonTCPModel
from src.models.logistic_tcp import LogisticTCPModel
from src.models.probit_tcp import ProbitTCPModel
from src.models.model_comparison import run_model_comparison
from src.analysis.rano_tcp_comparison import run_rano_tcp_comparison
from src.analysis.within_arm_rano_tcp import run_within_arm_rano_analysis""",
                "03-imports",
            ),
            cell_code(
                """modeling = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
median_os = modeling["survival_weeks"].median()
outcomes = (modeling["survival_weeks"] >= median_os).astype(float).to_numpy()
doses = modeling["eqd2_gy"].to_numpy()
print(f"n={len(modeling)}, OS median={median_os:.0f} wk, RANO labels={modeling['rano_controlled_t1'].notna().sum()}")""",
                "03-load",
            ),
            cell_code(
                """models = [
    ("Poisson", PoissonTCPModel(d50_init=55, gamma50_init=1.5)),
    ("Logistic", LogisticTCPModel(d50_init=53, k_init=0.1)),
    ("Probit", ProbitTCPModel(d50_init=53, sigma_init=10)),
]
d_grid = np.linspace(35, 65, 200)
fig, ax = plt.subplots(figsize=(7, 4.5))
for name, m in models:
    m.fit(doses, outcomes)
    ax.plot(d_grid, m.predict(d_grid), label=name, lw=2)
ax.set_xlabel("EQD2 (Gy)")
ax.set_ylabel("TCP (OS ≥ median proxy)")
ax.set_title(f"Pooled TCP (n={len(modeling)})")
ax.legend(frameon=False)
plt.tight_layout()
fig.savefig(FIGURES_DIR / "03_tcp_curves_os_proxy.png")
plt.close(fig)

comparison_df, _ = run_model_comparison(doses, outcomes, frame=modeling)
comparison_df[["model", "aic", "roc_auc", "hl_p_value"]]""",
                "03-pooled",
            ),
            cell_code(
                """metrics_dir = REPORTS_DIR / "metrics"
rano_cmp, _, _ = run_rano_tcp_comparison(modeling, metrics_dir)
rano_cmp[["model", "endpoint", "roc_auc_insample", "roc_auc_cv_mean", "lr_p_value"]]""",
                "03-rano",
            ),
            cell_code(
                """within_df, cox_rano = run_within_arm_rano_analysis(modeling, metrics_dir)
within_df[["scheme", "metric", "n_rano", "spearman_rho", "spearman_p", "poisson_auc", "poisson_lr_p"]]""",
                "03-within",
            ),
            cell_code(
                """cox_rano[["term", "hazard_ratio", "p_value"]]""",
                "03-cox",
            ),
        ],
    )


def build_04() -> None:
    write_nb(
        "04_parameter_estimation.ipynb",
        [
            cell_md("# 04 — Bootstrap CI (Poisson TCP, EQD2 / OS proxy)", "04-title"),
            cell_code(
                """import sys
from pathlib import Path
sys.path.insert(0, str(Path.cwd().parent))

import pandas as pd
from src.config import DATA_PROCESSED, REPORTS_DIR
from src.models.bootstrap_ci import bootstrap_tcp_params
from src.models.poisson_tcp import PoissonTCPModel
from src.models.model_comparison import run_model_comparison""",
                "04-imports",
            ),
            cell_code(
                """frame = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
median_os = frame["survival_weeks"].median()
outcomes = (frame["survival_weeks"] >= median_os).astype(float).to_numpy()
doses = frame["eqd2_gy"].to_numpy()

ci = bootstrap_tcp_params(
    PoissonTCPModel, doses, outcomes, n_bootstrap=1000,
    model_factory=lambda: PoissonTCPModel(d50_init=55.0, gamma50_init=1.5),
)
table, _ = run_model_comparison(doses, outcomes, frame=frame)
ci""",
                "04-run",
            ),
            cell_code("""table.sort_values("aic")[["model", "k_params", "aic", "bic", "roc_auc"]]""", "04-table"),
        ],
    )


def build_05() -> None:
    write_nb(
        "05_survival_analysis.ipynb",
        [
            cell_md("# 05 — Survival Analysis (KM, Cox, clinical stratification)", "05-title"),
            cell_code(
                """import sys
from pathlib import Path
sys.path.insert(0, str(Path.cwd().parent))

import pandas as pd
from src.config import DATA_PROCESSED, REPORTS_DIR
from src.models.survival_analysis import run_survival_analysis
from src.analysis.stratified_analysis import run_stratified_analysis""",
                "05-imports",
            ),
            cell_code(
                """frame = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
km, cox, _, _ = run_survival_analysis(frame)
clinical, who, within, hypo = run_stratified_analysis(frame, REPORTS_DIR / "metrics")
km""",
                "05-km",
            ),
            cell_code("""cox[["term", "hazard_ratio", "p"]]""", "05-cox-dose"),
            cell_code("""clinical""", "05-cox-clinical"),
            cell_code("""who""", "05-who"),
        ],
    )


def build_06() -> None:
    write_nb(
        "06_rano_multivariable_40gy.ipynb",
        [
            cell_md(
                """# 06 — Multivariable logistic: RANO ~ volume + age + PS (40 Gy arm)

Hypofractionated cohort (40.05 Gy / 15 fr) with RANO t0→t1 labels.
Compares volume-only vs adjusted model; validates DVH volume vs RANO TSV.""",
                "06-title",
            ),
            cell_code(
                """import sys
from pathlib import Path
sys.path.insert(0, str(Path.cwd().parent))

import matplotlib
matplotlib.use("Agg")
import pandas as pd

from src.config import DATA_PROCESSED, REPORTS_DIR
from src.analysis.rano_multivariable import run_rano_multivariable_40gy
from src.analysis.validate_rano_volumes import run_volume_validation""",
                "06-imports",
            ),
            cell_code(
                """modeling = pd.read_csv(DATA_PROCESSED / "modeling_table.csv")
metrics_dir = REPORTS_DIR / "metrics"
comp, coefs, boot, val = run_rano_multivariable_40gy(modeling, metrics_dir)
nifti_val = run_volume_validation(metrics_dir)
comp""",
                "06-run",
            ),
            cell_code("""coefs""", "06-coefs"),
            cell_code("""boot""", "06-boot"),
            cell_code("""val""", "06-val"),
            cell_code("""nifti_val""", "06-nifti"),
        ],
    )


if __name__ == "__main__":
    build_03()
    build_04()
    build_05()
    build_06()
