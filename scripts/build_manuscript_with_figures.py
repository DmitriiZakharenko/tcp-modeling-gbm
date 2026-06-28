#!/usr/bin/env python3
"""Insert numbered figures into manuscript draft for pandoc export."""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DRAFT = ROOT / "reports" / "manuscript_draft.md"
OUT = ROOT / "reports" / "manuscript_with_figures.md"

FIGURES = [
    (
        "## 3. Results",
        "figures/04_clinical_prognosis.png",
        "**Figure 1.** Kaplan–Meier overall survival by fractionation scheme and WHO "
        "performance status (n=190). Median OS 60 vs 28 weeks (log-rank p≈3×10⁻⁶).",
    ),
    (
        "### 3.2 Pooled TCP (OS proxy)",
        "figures/03_tcp_curves_os_proxy.png",
        "**Figure 2.** Fitted TCP curves (Poisson, logistic, probit, gEUD) for EQD2 vs "
        "OS median-split proxy (n=190; AUC≈0.68).",
    ),
    (
        "### 3.3 RANO endpoint (pooled)",
        "figures/05_rano_vs_os_tcp_auc.png",
        "**Figure 3.** Poisson TCP AUC: OS proxy vs RANO non-PD on the same 137 patients "
        "(AUC≈0.62 vs 0.43).",
    ),
    (
        "### 3.4 Within-arm analysis",
        "figures/06_within_arm_rano_tcp.png",
        "**Figure 4.** Within-arm DVH metrics → RANO: dose flat within protocol; GTV "
        "volume significant in 40 Gy arm (n=34).",
    ),
    (
        "### 3.8 LOOCV in the hypofractionated arm",
        "figures/07_rano_logistic_roc_40gy.png",
        "**Figure 5.** Multivariable logistic ROC for RANO non-PD in 40 Gy arm "
        "(in-sample AUC=0.90; LOOCV 0.74).",
    ),
    (
        "### 3.5 Volume validation",
        "figures/07_rano_volume_validation_40gy.png",
        "**Figure 6.** DVH GTV volume vs author RANO `size_t0_cm3` (Spearman ρ=1.00, n=141).",
    ),
    (
        "### 3.7 Pooled RANO models (volume + clinical + scheme)",
        "figures/08_pooled_rano_roc.png",
        "**Figure 7.** Pooled logistic ROC for RANO non-PD (volume + clinical + scheme; "
        "AUC≈0.72; LOOCV≈0.64).",
    ),
    (
        "### 3.9 PyRadiomics vs DVH volume",
        "figures/08_pyradiomics_vs_volume_auc.png",
        "**Figure 8.** PyRadiomics top-5 vs DVH volume for RANO (in-sample AUC 0.78 vs 0.71).",
    ),
    (
        "### 3.9 PyRadiomics vs DVH volume",
        "figures/08_pyradiomics_nested_cv_auc.png",
        "**Figure 9.** Nested 5-fold CV: in-sample vs out-of-sample AUC (optimism Δ≈0.04–0.07).",
        True,  # second insert at same section — append after first fig 8 block
    ),
]


def _block(path: str, caption: str) -> str:
    return (
        f"\n\n![{caption}]({path})\n\n"
        f"{caption}\n\n"
    )


def main() -> None:
    text = DRAFT.read_text(encoding="utf-8")
    # Remove old figure list table section from export (keep in draft)
    if "## Figure list (from repository)" in text:
        text = text.split("## Figure list (from repository)")[0].rstrip() + "\n"

    inserted_second_39 = False
    for item in FIGURES:
        marker, path, caption = item[0], item[1], item[2]
        append = len(item) > 3 and item[3]
        block = _block(path, caption)
        if append:
            # insert after Figure 8 block if already present
            needle = "figures/08_pyradiomics_vs_volume_auc.png"
            if needle in text and "figures/08_pyradiomics_nested_cv_auc.png" not in text:
                idx = text.index(needle)
                end = text.find("\n\n", text.find("**Figure 8.**", idx))
                text = text[: end + 2] + block + text[end + 2 :]
            inserted_second_39 = True
            continue
        if marker in text and path not in text:
            text = text.replace(marker, marker + block, 1)

    if "### 3.10" in text and "figures/08_pyradiomics_nested_cv_auc.png" not in text:
        block9 = _block(
            "figures/08_pyradiomics_nested_cv_auc.png",
            "**Figure 9.** Nested 5-fold CV: in-sample vs out-of-sample AUC (optimism Δ≈0.04–0.07).",
        )
        text = text.replace("### 3.10", block9 + "### 3.10", 1)

    OUT.write_text(text, encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
