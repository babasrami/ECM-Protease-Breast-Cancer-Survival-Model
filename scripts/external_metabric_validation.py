"""
METABRIC External Validation Script
===================================

External validation of a locked TCGA-BRCA trained Random Survival Forest model
on the independent METABRIC breast cancer cohort.

No data leakage:
- The RSF model is loaded as a locked model.
- No retraining or hyperparameter tuning is performed on METABRIC.
- METABRIC normalization is fitted only on METABRIC samples.
- Risk-group cutoffs are derived from the original training cohort.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import warnings
from pathlib import Path

import joblib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from sklearn.preprocessing import StandardScaler
from sksurv.util import Surv

warnings.filterwarnings("ignore")

GENE_ALIAS = {"MMP23": "MMP23B"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate locked RSF model on METABRIC.")
    parser.add_argument("--base-dir", default="results_survival_pipeline")
    parser.add_argument("--output-dir", default="Results")
    parser.add_argument("--model-path", default="results_survival_pipeline/models/final_integrative_rsf.joblib")
    parser.add_argument("--train-data", default="Data.xlsx")
    parser.add_argument("--metabric-dir", default="brca_metabric")
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    parser.add_argument("--random-seed", type=int, default=1000)
    parser.add_argument("--fig-dpi", type=int, default=300)
    return parser.parse_args()


def find_file(root: Path, filename: str) -> Path:
    """Find file at root/filename or recursively below root."""
    direct = root / filename
    if direct.is_file():
        return direct
    matches = list(root.rglob(filename))
    if not matches:
        return direct
    return matches[0]


def make_dirs(base_dir: Path, output_dir: Path) -> None:
    for folder in [
        base_dir / "tables",
        base_dir / "figures",
        output_dir / "tables",
        output_dir / "figures",
    ]:
        folder.mkdir(parents=True, exist_ok=True)


def check_files(paths: dict[str, Path]) -> None:
    print("\n" + "=" * 60)
    print("FILE EXISTENCE CHECK")
    print("=" * 60)

    missing = []
    for label, path in paths.items():
        exists = path.is_file()
        print(f"  {'[OK]' if exists else '[MISSING]'} {label}: {path}")
        if not exists:
            missing.append(f"{label}: {path}")

    if missing:
        raise FileNotFoundError("Missing required files:\n" + "\n".join(missing))

    print("\nAll files present. Proceeding.\n")


def load_model(model_path: Path):
    artifact = joblib.load(model_path)

    required_keys = {"model", "features", "ptpm_cols"}
    missing = required_keys.difference(artifact.keys())
    if missing:
        raise KeyError(f"Model artifact missing keys: {sorted(missing)}")

    rsf = artifact["model"]
    features = artifact["features"]
    ptpm_cols = artifact["ptpm_cols"]

    print(f"[Model] Loaded: {model_path}")
    print(f"[Model] Features used ({len(features)}): {features}")
    return rsf, features, ptpm_cols


def load_metabric_clinical(patient_file: Path, sample_file: Path) -> pd.DataFrame:
    clin_p = pd.read_csv(patient_file, sep="\t", comment="#")
    clin_s = pd.read_csv(sample_file, sep="\t", comment="#")

    clin_p["event"] = clin_p["OS_STATUS"].str.startswith("1").fillna(False).astype(int)
    clin_p["time"] = pd.to_numeric(clin_p["OS_MONTHS"], errors="coerce") * 30.4375
    clin_p["age"] = pd.to_numeric(clin_p["AGE_AT_DIAGNOSIS"], errors="coerce")

    stage_map = {1: 1, 2: 2, 3: 3, 4: 4}
    clin_s = clin_s.copy()
    clin_s["stage_ordinal"] = (
        pd.to_numeric(clin_s["TUMOR_STAGE"], errors="coerce")
        .map(stage_map)
        .fillna(0)
        .astype(float)
    )

    merged = clin_p.merge(
        clin_s[["PATIENT_ID", "stage_ordinal"]],
        on="PATIENT_ID",
        how="left",
    )
    merged["stage_ordinal"] = merged["stage_ordinal"].fillna(0)
    merged["Subtype"] = merged["CLAUDIN_SUBTYPE"].fillna("Unknown").astype(str)

    merged = merged.dropna(subset=["time", "event"])
    merged = merged[merged["time"] > 0]
    merged = merged.set_index("PATIENT_ID")

    print(f"[Clinical] METABRIC patients with valid survival: {len(merged)}")
    print(f"[Clinical] Events: {int(merged['event'].sum())}")
    print(f"[Clinical] Censored: {int((merged['event'] == 0).sum())}")
    return merged


def load_metabric_expression(expr_file: Path, needed_base_names: list[str]) -> tuple[pd.DataFrame, list[str]]:
    """Read only required protease genes from large METABRIC expression file."""
    metabric_lookup = {}
    for base in needed_base_names:
        metabric_lookup[GENE_ALIAS.get(base, base)] = base

    print(f"[Expression] Scanning {expr_file} for {len(metabric_lookup)} protease genes")
    print("[Expression] Alias applied: MMP23 -> MMP23B")

    collected = {}
    with expr_file.open("r", encoding="utf-8") as fh:
        header = fh.readline().strip().split("\t")
        sample_ids = header[2:]

        for line in fh:
            parts = line.rstrip("\n").split("\t")
            gene = parts[0]
            if gene in metabric_lookup:
                training_name = metabric_lookup[gene]
                if training_name not in collected:
                    collected[training_name] = dict(zip(sample_ids, parts[2:]))

    data = {}
    for base_name, sample_values in collected.items():
        data[f"{base_name}_pTPM"] = pd.to_numeric(pd.Series(sample_values), errors="coerce")

    expr_df = pd.DataFrame(data)
    expr_df.index.name = "PATIENT_ID"

    missing = [gene for gene in needed_base_names if gene not in collected]
    print(f"[Expression] Genes found: {len(collected)}/{len(needed_base_names)}")
    print(f"[Expression] Missing genes: {missing if missing else 'None'}")
    print(f"[Expression] Expression samples: {len(expr_df)}")
    return expr_df, missing


def preprocess_metabric(df: pd.DataFrame, final_features: list[str]) -> pd.DataFrame:
    out = df.copy()
    protease_features = [f for f in final_features if f.endswith("_pTPM")]

    for col in protease_features:
        if col not in out.columns:
            out[col] = 0.0
            print(f"[WARNING] Missing feature filled with 0: {col}")
        out[col] = np.log1p(pd.to_numeric(out[col], errors="coerce").fillna(0))

    if "age" in final_features:
        out["age"] = pd.to_numeric(out["age"], errors="coerce").fillna(out["age"].median())

    if "stage_ordinal" in final_features:
        out["stage_ordinal"] = out["stage_ordinal"].fillna(0).astype(float)

    scale_cols = [f for f in final_features if f in out.columns]
    scaler = StandardScaler()
    out[scale_cols] = scaler.fit_transform(out[scale_cols])
    return out


def get_training_cutoffs(train_data: Path, rsf_model, final_features: list[str]) -> tuple[float, float]:
    """Re-derive training risk-score tertile cutoffs without using METABRIC outcomes."""
    raw = pd.read_excel(train_data)

    status_map = {
        "dead": 1, "died": 1, "deceased": 1, "1": 1,
        "alive": 0, "living": 0, "0": 0,
    }
    raw["event"] = raw["status"].astype(str).str.lower().str.strip().map(status_map)
    raw["time"] = pd.to_numeric(raw["days"], errors="coerce")
    raw = raw.dropna(subset=["time", "event"])
    raw = raw[raw["time"] > 0].copy()

    raw["age"] = pd.to_numeric(raw["age"], errors="coerce").fillna(raw["age"].median())

    roman_map = {"I": 1, "II": 2, "III": 3, "IV": 4}
    extracted = raw["stage"].astype(str).str.extract(r"Stage\s*([IVX]+)", flags=re.IGNORECASE)[0]
    raw["stage_ordinal"] = extracted.str.upper().map(roman_map).fillna(0).astype(float)

    ptpm_cols = [
        c for c in raw.columns
        if ("MMP" in str(c) or "ADAM" in str(c)) and "_pTPM" in str(c)
    ]
    for col in ptpm_cols:
        raw[col] = np.log1p(pd.to_numeric(raw[col], errors="coerce").fillna(0))

    scale_cols = ptpm_cols + ["age", "stage_ordinal"]
    scaler = StandardScaler()
    raw[scale_cols] = scaler.fit_transform(raw[scale_cols])

    missing = [f for f in final_features if f not in raw.columns]
    if missing:
        raise ValueError(f"Training data is missing model features: {missing}")

    train_scores = rsf_model.predict(raw[final_features].values)
    q33, q67 = np.percentile(train_scores, [33.33, 66.67])
    print(f"[Cutoffs] Training tertiles: 33rd={q33:.4f}, 67th={q67:.4f}")
    return float(q33), float(q67)


def bootstrap_external_cindex(rsf_model, X: np.ndarray, y, n_boot: int, seed: int) -> tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    scores = []

    for _ in range(n_boot):
        idx = rng.integers(0, len(X), size=len(X))
        y_b = y[idx]
        if y_b["event"].sum() < 1:
            continue
        try:
            scores.append(float(rsf_model.score(X[idx], y_b)))
        except Exception:
            continue

    if not scores:
        return float("nan"), float("nan"), float("nan")

    arr = np.asarray(scores)
    return float(arr.mean()), float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))


def plot_external_km(df: pd.DataFrame, logrank_p: float, out_path: Path, dpi: int) -> None:
    colors = {"Low": "#2A9D8F", "Medium": "#E9C46A", "High": "#E76F51"}
    fig, ax = plt.subplots(figsize=(7, 5.2))

    for group in ["Low", "Medium", "High"]:
        mask = df["risk_group"].astype(str) == group
        if mask.sum() == 0:
            continue
        kmf = KaplanMeierFitter()
        kmf.fit(df.loc[mask, "time"], df.loc[mask, "event"], label=f"{group} (n={mask.sum()})")
        kmf.plot_survival_function(ax=ax, color=colors[group], linewidth=1.8, ci_show=True, ci_alpha=0.12)

    p_text = f"Log-rank p = {logrank_p:.4f}" if logrank_p >= 0.0001 else "Log-rank p < 0.0001"
    ax.text(
        0.97, 0.97, p_text,
        transform=ax.transAxes,
        fontsize=9,
        ha="right",
        va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#aaaaaa", alpha=0.85),
    )

    ax.set_title("METABRIC External Validation — KM Curves by RSF Risk Group")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Overall Survival Probability")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[Figure] Saved: {out_path}")


def save_table(df: pd.DataFrame, filename: str, base_dir: Path, output_dir: Path) -> None:
    p1 = base_dir / "tables" / filename
    p2 = output_dir / "tables" / filename
    df.to_csv(p1, index=False)
    shutil.copy2(p1, p2)
    print(f"[Table] Saved: {p2}")


def run(args: argparse.Namespace) -> dict:
    base_dir = Path(args.base_dir)
    output_dir = Path(args.output_dir)
    metabric_dir = Path(args.metabric_dir)
    model_path = Path(args.model_path)
    train_data = Path(args.train_data)

    expr_file = find_file(metabric_dir, "data_mrna_illumina_microarray.txt")
    patient_file = find_file(metabric_dir, "data_clinical_patient.txt")
    sample_file = find_file(metabric_dir, "data_clinical_sample.txt")

    make_dirs(base_dir, output_dir)
    check_files({
        "Locked RSF model": model_path,
        "Training data": train_data,
        "METABRIC expression": expr_file,
        "METABRIC clinical patient": patient_file,
        "METABRIC clinical sample": sample_file,
    })

    rsf, final_features, ptpm_cols = load_model(model_path)
    protease_base_names = [f.replace("_pTPM", "") for f in ptpm_cols]

    clin_df = load_metabric_clinical(patient_file, sample_file)
    expr_df, missing_genes = load_metabric_expression(expr_file, protease_base_names)
    metabric_df = clin_df.join(expr_df, how="inner")

    print(f"[Merge] Patients with clinical + expression: {len(metabric_df)}")
    print(f"[METABRIC] Events: {int(metabric_df['event'].sum())}")
    print(f"[METABRIC] Missing protease genes: {missing_genes if missing_genes else 'None'}")

    processed = preprocess_metabric(metabric_df, final_features)
    for feature in final_features:
        if feature not in processed.columns:
            processed[feature] = 0.0
            print(f"[WARNING] Added missing model feature as zero: {feature}")

    X_ext = processed[final_features].values
    y_ext = Surv.from_arrays(
        event=processed["event"].astype(bool).values,
        time=processed["time"].values,
    )

    ext_cindex = float(rsf.score(X_ext, y_ext))
    ci_mean, ci_lo, ci_hi = bootstrap_external_cindex(
        rsf, X_ext, y_ext, args.n_bootstrap, args.random_seed
    )

    print(f"\n[External C-index] {ext_cindex:.4f}")
    print(f"[Bootstrap 95% CI] {ci_lo:.4f} - {ci_hi:.4f} (mean={ci_mean:.4f})")

    q33, q67 = get_training_cutoffs(train_data, rsf, final_features)
    processed["risk_score"] = rsf.predict(X_ext)
    processed["risk_group"] = pd.cut(
        processed["risk_score"],
        bins=[-np.inf, q33, q67, np.inf],
        labels=["Low", "Medium", "High"],
    )

    valid = processed.dropna(subset=["risk_group", "time", "event"]).copy()
    valid = valid[valid["risk_group"].astype(str).isin(["Low", "Medium", "High"])]

    lr = multivariate_logrank_test(
        event_durations=valid["time"],
        groups=valid["risk_group"].astype(str),
        event_observed=valid["event"],
    )
    logrank_p = float(lr.p_value)
    logrank_chi2 = float(lr.test_statistic)
    p_label = f"{logrank_p:.4f}" if logrank_p >= 0.0001 else "< 0.0001"
    print(f"[Log-rank] chi2={logrank_chi2:.4f}, p={p_label}")

    km_path = base_dir / "figures" / "external_metabric_km_by_risk_group.png"
    plot_external_km(valid, logrank_p, km_path, args.fig_dpi)
    shutil.copy2(km_path, output_dir / "figures" / "external_metabric_km_by_risk_group.png")

    val_summary = pd.DataFrame([{
        "Cohort": "METABRIC (External)",
        "N_Patients": len(processed),
        "N_Events": int(processed["event"].sum()),
        "N_Censored": int((processed["event"] == 0).sum()),
        "N_Features_Used": len(final_features),
        "Missing_Protease_Genes": ";".join(missing_genes) if missing_genes else "None",
        "External_C_Index": round(ext_cindex, 4),
        "Bootstrap_Mean_C_Index": round(ci_mean, 4),
        "Bootstrap_95CI_Low": round(ci_lo, 4),
        "Bootstrap_95CI_High": round(ci_hi, 4),
        "Logrank_Chi2": round(logrank_chi2, 4),
        "Logrank_P_Value": round(logrank_p, 6),
        "Logrank_Significant": logrank_p < 0.05,
        "Risk_Cutoff_33rd_Training": round(q33, 4),
        "Risk_Cutoff_67th_Training": round(q67, 4),
    }])
    save_table(val_summary, "external_metabric_validation.csv", base_dir, output_dir)

    pred_out = processed[["time", "event", "risk_score", "risk_group"]].copy()
    pred_out.index.name = "PATIENT_ID"
    save_table(pred_out.reset_index(), "external_metabric_predictions.csv", base_dir, output_dir)

    missing_df = pd.DataFrame({"Missing_Feature": missing_genes if missing_genes else []})
    save_table(missing_df, "external_metabric_missing_features.csv", base_dir, output_dir)

    summary_path = output_dir / "pipeline_summary.json"
    summary = {}
    if summary_path.is_file():
        with summary_path.open("r", encoding="utf-8") as f:
            summary = json.load(f)

    summary["external_validation_metabric"] = {
        "n_patients": len(processed),
        "n_events": int(processed["event"].sum()),
        "missing_genes": missing_genes,
        "external_c_index": round(ext_cindex, 4),
        "bootstrap_mean": round(ci_mean, 4),
        "bootstrap_95ci_low": round(ci_lo, 4),
        "bootstrap_95ci_high": round(ci_hi, 4),
        "logrank_chi2": round(logrank_chi2, 4),
        "logrank_p_value": round(logrank_p, 6),
        "logrank_significant": logrank_p < 0.05,
    }
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 60)
    print("EXTERNAL VALIDATION COMPLETE")
    print("=" * 60)
    print(f"  METABRIC patients matched : {len(processed)}")
    print(f"  METABRIC events           : {int(processed['event'].sum())}")
    print(f"  Missing protease genes    : {missing_genes if missing_genes else 'None'}")
    print(f"  External C-index          : {ext_cindex:.4f}")
    print(f"  Bootstrap 95% CI          : {ci_lo:.4f} - {ci_hi:.4f}")
    print(f"  Log-rank p-value          : {p_label}")
    print("=" * 60)

    return summary["external_validation_metabric"]


if __name__ == "__main__":
    run(parse_args())
