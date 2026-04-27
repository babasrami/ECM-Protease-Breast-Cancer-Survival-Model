# ECM Protease Breast Cancer Survival Model: TCGA Training and METABRIC External Validation

This repository contains the full analysis workflow used to develop and externally validate a breast cancer survival prediction model based on extracellular matrix (ECM) protease transcript signatures and clinical variables.

The project starts from TCGA-BRCA data preprocessing, protease feature screening, model training, internal validation, patient risk stratification, and ends with independent external validation in the METABRIC cohort.

---

## Project gap

Traditional breast cancer prognostic tools rely mainly on clinical variables such as age, tumour stage, tumour size, lymph node status, or histological grade. However, patients with similar clinical profiles can still have very different survival outcomes because breast cancer is molecularly heterogeneous.

ECM remodelling is strongly involved in tumour invasion, metastasis, and disease progression. Protease families such as MMP, ADAM, and ADAMTS may carry useful prognostic information, but their combined transcriptomic value with clinical variables needs to be tested in a survival-machine-learning framework.

---

## Project need

A complete and reproducible pipeline was needed to answer two main questions:

1. Do ECM protease transcript signatures improve survival prediction beyond clinical variables alone?
2. Does the final trained model retain prognostic signal when tested in an independent external cohort?

Internal validation is useful, but it is not enough. A model can perform well inside the discovery dataset and still lose performance when applied to a different cohort, platform, or patient population. Therefore, external validation was required.

---

## Project goal

The goal of this project was to build an integrative breast cancer survival prediction model using:

- ECM protease transcript features
- Age
- Clinical stage
- Machine learning survival models

The final model was trained on TCGA-BRCA and then tested on METABRIC as an independent external validation cohort without data leakage, retraining, feature re-selection, or hyperparameter tuning.

---

## What we did

The full workflow included:

1. Loaded and processed TCGA-BRCA clinical and transcriptomic data.
2. Harmonized survival time and event status.
3. Processed clinical variables, including age and clinical stage.
4. Extracted ECM protease transcript features from MMP, ADAM, and ADAMTS families.
5. Applied expression preprocessing, including log1p transformation and z-score scaling.
6. Performed univariate Cox screening to identify top prognostic protease candidates.
7. Performed subtype-specific survival association analysis.
8. Compared different survival model families:
   - Cox Proportional Hazards
   - Penalized Cox / CoxNet
   - Random Survival Forest
   - Gradient Boosting Survival
9. Compared protease-only models against integrated clinical + protease models.
10. Selected the best-performing integrative model.
11. Trained the final Random Survival Forest model on TCGA-BRCA.
12. Generated patient-level risk scores.
13. Stratified patients into Low, Medium, and High-risk groups.
14. Evaluated internal performance using held-out testing and cross-validation.
15. Downloaded and processed METABRIC clinical and expression data.
16. Applied the locked TCGA-trained model to METABRIC.
17. Used training-derived cutoffs to assign METABRIC risk groups.
18. Computed the external C-index with 1,000-bootstrap 95% confidence interval.
19. Performed Kaplan-Meier survival analysis and multivariate log-rank testing.
20. Saved validation tables, patient-level predictions, missing-feature checks, figures, and summary outputs.

---

## Main achieved results

The integrative TCGA-BRCA model showed strong internal performance, with the Random Survival Forest model achieving the best discrimination among the tested model families.

In external validation, the locked TCGA-trained model was tested in METABRIC.

External validation results:

- External cohort: METABRIC
- Evaluable patients: 1,979
- Events: 1,143
- External C-index: 0.581
- Bootstrap 95% CI: 0.562-0.598
- Log-rank test across predicted risk groups: p < 0.0001

These results show that the model retained a statistically detectable but modest independent prognostic signal in METABRIC. However, the modest external C-index means that the model is not yet clinically confirmatory. Further validation in uniformly annotated prospective cohorts is required.

---

## Repository structure

```text
.
├── README.md
├── requirements.txt
├── .gitignore
├── Code.py
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
    ├── models/
    │   └── README.md
    ├── tables/
    │   └── README.md
    └── figures/
        └── README.md
```

`Code.py` contains the main TCGA-BRCA pipeline for preprocessing, feature screening, model training, internal validation, and final model generation.

`scripts/external_metabric_validation.py` contains the METABRIC external validation workflow.

---

## Required input files

This repository does not include large public datasets or generated model files.

For the main TCGA-BRCA training pipeline, you need:

```text
Data.xlsx
```

For the METABRIC external validation step, you need:

```text
results_survival_pipeline/models/final_integrative_rsf.joblib
brca_metabric/data_mrna_illumina_microarray.txt
brca_metabric/data_clinical_patient.txt
brca_metabric/data_clinical_sample.txt
```

The file `final_integrative_rsf.joblib` is generated by running the main TCGA-BRCA training pipeline first.

---

## Download METABRIC data

The METABRIC dataset can be downloaded from cBioPortal/DataHub:

```bash
wget https://datahub.assets.cbioportal.org/brca_metabric.tar.gz
mkdir -p brca_metabric
tar -xzf brca_metabric.tar.gz -C brca_metabric --strip-components=1
```

If the extraction creates a nested folder, move the following files into `brca_metabric/`:

```text
data_mrna_illumina_microarray.txt
data_clinical_patient.txt
data_clinical_sample.txt
```

---

## Installation

Install the required Python packages:

```bash
pip install -r requirements.txt
```

If `scikit-survival` fails with pip, install it using conda:

```bash
conda install -c conda-forge scikit-survival
```

---

## Run the full workflow

### 1. Run the main TCGA-BRCA training pipeline

```bash
python Code.py
```

This step performs the main analysis and generates the trained model:

```text
results_survival_pipeline/models/final_integrative_rsf.joblib
```

It also generates internal validation tables and figures.

### 2. Run METABRIC external validation

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

---

## Main output files

Training and internal validation outputs include:

```text
results_survival_pipeline/models/final_integrative_rsf.joblib
results_survival_pipeline/tables/
results_survival_pipeline/figures/
```

External validation outputs include:

```text
Results/tables/external_metabric_validation.csv
Results/tables/external_metabric_predictions.csv
Results/tables/external_metabric_missing_features.csv
Results/figures/external_metabric_km_by_risk_group.png
Results/pipeline_summary.json
```

---

## External validation outputs

### `external_metabric_validation.csv`

Summary table containing:

- Number of METABRIC patients
- Number of events
- Number of censored patients
- Number of model features used
- Missing protease genes
- External C-index
- Bootstrap 95% confidence interval
- Log-rank statistic
- Log-rank p-value
- Training-derived risk cutoffs

### `external_metabric_predictions.csv`

Patient-level METABRIC predictions containing:

- Patient ID
- Survival time
- Event status
- Predicted risk score
- Assigned risk group

### `external_metabric_missing_features.csv`

List of protease features missing from METABRIC, if any.

### `external_metabric_km_by_risk_group.png`

Kaplan-Meier plot showing overall survival separation between predicted Low, Medium, and High-risk groups in METABRIC.

---

## Interpretation

The METABRIC validation should be interpreted cautiously.

The lower external C-index may reflect:

- TCGA vs METABRIC cohort differences
- RNA-seq pTPM vs Illumina microarray platform differences
- Batch effects
- Clinical-variable mismatch
- Endpoint and follow-up differences
- Feature harmonization constraints
- Biological heterogeneity
- Possible cohort-specific overfitting in the internally strong model

Therefore, this project supports ECM protease transcripts as promising candidate prognostic markers, but it does not establish a clinically ready model. Prospective validation using standardized clinical annotation, treatment information, breast cancer-specific endpoints, and ideally proteomic or activity-based protease measurements is still required.

---

## Data source

TCGA-BRCA data were obtained from public TCGA/UCSC Xena resources.

METABRIC data were obtained from cBioPortal/DataHub:

```text
https://datahub.assets.cbioportal.org/brca_metabric.tar.gz
```

The original METABRIC publications should also be cited in the manuscript/reference list.

---

## Notes

- The METABRIC validation uses the locked TCGA-trained model.
- No METABRIC data are used for model training.
- No METABRIC-based hyperparameter tuning is performed.
- Risk-group cutoffs are derived from the TCGA training cohort and applied blindly to METABRIC.
- MMP23 in the TCGA training data corresponds to MMP23B in METABRIC and is mapped automatically.

---

## License

Add your preferred license here.

Example:

```text
MIT License
```

---

## Citation

If you use this repository, please cite the related manuscript:

```text
Babas, R.; Vynios, D.H.; Kompothrekas, A.; Boutsinas, B.; Karamanos, N.K.
Integrative Survival Prediction in Breast Cancer Using Extracellular Matrix Protease Transcript Signatures and Clinical Variables: A Machine Learning Approach.
Cancers, 2026.
```
