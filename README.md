# METABRIC External Validation of a TCGA-Trained Breast Cancer Survival Model

This repository contains the external validation script used to test a locked Random Survival Forest (RSF) survival model on the independent METABRIC breast cancer cohort.

## Gap

The original survival model was developed using TCGA-BRCA data. Internal validation alone cannot prove that the model generalizes to an independent cohort.

## Need

An external validation step was required to check whether the trained model still carries prognostic signal in a different breast cancer cohort generated using a different expression platform.

## Goal

Validate the locked TCGA-trained integrative RSF model on METABRIC without data leakage, retraining, feature re-selection, or hyperparameter tuning.

## What we did

The script:

1. Loads the locked RSF model trained only on TCGA-BRCA.
2. Reads METABRIC clinical and Illumina microarray expression data.
3. Extracts the same ECM protease features used during training.
4. Applies METABRIC-only preprocessing and independent z-score normalization.
5. Predicts patient-level risk scores using the locked model.
6. Applies training-derived risk-group cutoffs to METABRIC.
7. Computes external C-index with 1,000-bootstrap 95% confidence interval.
8. Performs Kaplan-Meier survival analysis and multivariate log-rank testing.
9. Saves validation tables, predictions, missing-feature checks, and the KM figure.

## Main achieved result

External validation in METABRIC showed a statistically detectable but modest prognostic signal:

- External cohort: METABRIC
- Evaluable patients: 1,979
- Events: 1,143
- External C-index: 0.581
- Bootstrap 95% CI: 0.562-0.598
- Log-rank test across predicted risk groups: p < 0.0001

These results support that the model retained some external prognostic signal, but the modest C-index means the model is not yet clinically confirmatory. Further validation in uniformly annotated prospective cohorts is required.

## Repository structure

```text
.
├── README.md
├── requirements.txt
├── .gitignore
├── scripts/
│   └── external_metabric_validation.py
├── data/
│   └── README.md
├── Results/
│   ├── tables/
│   │   └── README.md
│   └── figures/
│       └── README.md
└── results_survival_pipeline/
    └── models/
        └── README.md
```

## Required input files

This repository does **not** include large public datasets or local model/training files.

You need:

```text
Data.xlsx
results_survival_pipeline/models/final_integrative_rsf.joblib
brca_metabric/data_mrna_illumina_microarray.txt
brca_metabric/data_clinical_patient.txt
brca_metabric/data_clinical_sample.txt
```

## Download METABRIC data

```bash
wget https://datahub.assets.cbioportal.org/brca_metabric.tar.gz
mkdir -p brca_metabric
tar -xzf brca_metabric.tar.gz -C brca_metabric --strip-components=1
```

## Installation

```bash
pip install -r requirements.txt
```

If `scikit-survival` fails with pip, install it with conda:

```bash
conda install -c conda-forge scikit-survival
```

## Run

```bash
python scripts/external_metabric_validation.py
```

Optional custom paths:

```bash
python scripts/external_metabric_validation.py \
  --model-path results_survival_pipeline/models/final_integrative_rsf.joblib \
  --train-data Data.xlsx \
  --metabric-dir brca_metabric \
  --output-dir Results
```

## Output files

```text
Results/tables/external_metabric_validation.csv
Results/tables/external_metabric_predictions.csv
Results/tables/external_metabric_missing_features.csv
Results/figures/external_metabric_km_by_risk_group.png
Results/pipeline_summary.json
```

## Interpretation

The METABRIC validation should be interpreted cautiously. The lower external C-index may reflect:

- TCGA vs METABRIC cohort differences
- RNA-seq pTPM vs Illumina microarray platform differences
- Batch effects
- Clinical-variable mismatch
- Endpoint and follow-up differences
- Feature-harmonization constraints
- Biological heterogeneity
- Possible cohort-specific overfitting in the internally strong RSF model

Therefore, this external validation supports a modest independent prognostic signal, not definitive clinical readiness.

## Data source

METABRIC data were obtained from cBioPortal/DataHub:

```text
https://datahub.assets.cbioportal.org/brca_metabric.tar.gz
```

Original METABRIC publications should also be cited in the manuscript/reference list.
