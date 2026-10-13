
# Incremental value and conditional payer consequences of a residual disease prediction rule in HR+/HER2− breast cancer

Analysis code accompanying the manuscript *"Incremental value and conditional payer consequences of a residual disease prediction rule in HR+/HER2− breast cancer"*.

The code was written during a doctoral thesis, so its comments, variable names and printed output are in Turkish. This README provides an English map of the pipeline, a Turkish–English glossary, and a table linking each script to the tables and figures in the manuscript. The wider set of scripts from the doctoral thesis is at https://github.com/kdagdeviren/Doktora-Tez; this repository contains only the code specific to the present manuscript.

---

## What is and is not here

The scripts reproduce the numbers, tables and figures of the manuscript and its supplementary materials. Three checks are worth stating, because they are what a reader would want to verify first: running `04_economics/markov_base.py` unchanged reproduces the published economic figures exactly; running `03_validation/bootstrap_dca.py` reproduces the main-pipeline classification results exactly (108 flagged, 98 concordant, 10 discordant, positive predictive value 0.907); and running `06_erpor/fn_subtype_split.py` reproduces the fully nested results exactly (113 flagged, 97 concordant, 16 discordant, positive predictive value 0.858, area under the curve 0.694 in the eligible population). **The patient-level dataset is not included.** Individual records contain clinical information that could permit identification, and the ethics approval does not cover public release; the anonymised dataset is available from the corresponding author on reasonable request. Scripts expect the cohort file at the path given at the top of each file.

## Directory structure
01_cohort/ data correction, audit and baseline table
02_model/ model architecture, layer comparison, alternative specifications
03_validation/ cross-validation, nesting, bootstrap, decision curves, thresholds
04_economics/ drug and treatment costs, decision tree, Markov model, scenarios
05_figures/ figures for the manuscript and the supplementary materials
06_erpor/ the fully nested re-runs that produce the primary results, and the figures redrawn from them
outputs/ CSV results written by the scripts
figures/ figures written by the scripts


## Pipeline order

Scripts are meant to be run in this order; later steps consume files written by earlier ones.

| Step | Script | What it does |
|---|---|---|
| 1 | `01_cohort/kohort_duzeltme.py` | Builds the analysis cohort from the raw file and logs every correction with a version number |
| 2 | `01_cohort/kohort_denetimi.py`, `kohort_denetimi_tur2.py`, `klinik_denetim.py` | Structural, clinical and cross-consistency checks, including the outcome-leakage test on the two radiology coding versions |
| 3 | `01_cohort/baseline_table.py` | Table 1 |
| 4 | `02_model/katman_analizi.py` | Fits the two layers and their combinations; writes the out-of-fold probabilities of the main pipeline |
| 5 | `02_model/alt_spec.py`, `duyarlilik_modelleri.py`, `radyoloji_kademeli.py` | Radiology-free, pathology-only and low-missingness specifications; the radiology ablation sequence |
| 6 | `03_validation/bootstrap_dca.py` | Base-case results, patient-level cluster bootstrap, decision curve analysis |
| 7 | `03_validation/icice_cv.py` | Nested cross-validation for backbone selection |
| 8 | `03_validation/fold_imputasyon.py` | Fold-internal modal and model-based imputation; dummy-coded boosting |
| 9 | `03_validation/tam_icice.py` | **Fully nested protocol** — imputation, coding, block selection and variable selection all inside the outer fold. Produces the primary estimate (AUC 0.694) |
| 10 | `03_validation/simple_model_rival.py` | Regularised logistic regression alone through the identical protocol |
| 11 | `03_validation/esik_egrisi.py` | Threshold sweep from 0.60 to 0.90 with bootstrap intervals |
| 12 | `04_economics/ilac_maliyet.py`, `maliyet_modeli.py` | Payer drug prices and itemised treatment costs |
| 13 | `04_economics/markov_base.py` | Decision tree plus lifetime Markov model; base case, PSA, tornado, acceptability curve |
| 14 | `04_economics/yapisal_senaryolar.py` | Seventeen structural scenarios, including adjuvant CDK4/6 uptake |
| 15 | `04_economics/markov_evpi.py` | Expected value of perfect information; saves the probabilistic draws |
| 16 | `04_economics/false_flag_sensitivity.py` | Harm per wrong flag (*H*) sweep and break-even values |
| 17 | `05_figures/*.py` | Figure panels |
| 18 | `06_erpor/*.py` | The fully nested re-runs that produce the primary results, and the figures redrawn from them |

## Which script produces which result

The manuscript reports the **fully nested** analysis as primary and the conventional main pipeline alongside it as a measure of optimism. Where both appear in one table, both scripts are listed.

| Manuscript item | Script |
|---|---|
| Table 1, baseline characteristics | `01_cohort/baseline_table.py` |
| Table 2, fully nested column | `06_erpor/fn_subtype_split.py` (runs `03_validation/tam_icice.py`) |
| Table 2, main pipeline column | `03_validation/bootstrap_dca.py` |
| Table 3, incremental classification value | `06_erpor/fn_subtype_split.py` and `03_validation/bootstrap_dca.py` |
| Table 4, modeled costs and outcomes | `06_erpor/fn_economics.py` (runs `04_economics/markov_base.py`) |
| Table 5, structural scenarios | `04_economics/yapisal_senaryolar.py` |
| Table 5, CDK4/6 and adherence rows | `06_erpor/fn_economics.py` |
| Table 5, false-flag rows | `04_economics/false_flag_sensitivity.py` |
| Figure 1, patient flow | `05_figures/fn_fig1.py` |
| Figure 2, model structure | `05_figures/fn_fig2.py` |
| Figure 3, calibration and decision curve | `05_figures/calibration_and_dca_figures.py` |
| Figure 4a, tornado | `04_economics/mb_fn.py` |
| Figure 4b, acceptability curve, three strategies | `06_erpor/ceac_three_strategy.py` |
| Figure 5, net saving against false-flag harm | `04_economics/false_flag_sensitivity.py` |
| Section 3.2, range across outer partitions | `06_erpor/fn_subtype_split.py` |
| Supplementary Figures S1 and S2 | `06_erpor/supp_fig_s1_s2.py` |
| Supplementary Figures S3 and S4 | `06_erpor/supp_fig_s3_s4.py` |
| Supplementary material 3, data quality | `01_cohort/kohort_denetimi.py`, `kohort_denetimi_tur2.py`, `klinik_denetim.py` |
| Supplementary material 4, validation layers | `03_validation/tam_icice.py`, `fold_imputasyon.py` |
| Supplementary material 4, single-layer comparison | `03_validation/simple_model_rival.py` |
| Supplementary material 4, threshold tables | `03_validation/esik_egrisi.py` |
| Supplementary material 5, probabilistic parameters | `04_economics/psa_tablo.py` |
| Supplementary material 8, post-hoc combined strategy | `03_validation/s3_test.py` |

Two things this repository does **not** contain, stated here rather than left to be discovered:

- Supplementary Figures S5 and S6 have no script of their own. They were drawn from the output of `03_validation/esik_egrisi.py` and `03_validation/s3_test.py` respectively, but the plotting code was not kept.
- The adherence analysis had no script of its own. It is now derived explicitly in `06_erpor/fn_economics.py` from the cost structure of `markov_base.py`, and at 108 flagged that derivation reproduces the previously published break-even of 58.5%.

## Reproducibility

Every operation involving randomness uses a fixed seed.

| Operation | Seed |
|---|---|
| Cross-validation (5 repeats) | 300–304 |
| Nested cross-validation (3 repeats) | 500–502 |
| Fully nested pipeline (3 repeats) | 700–702 |
| Patient-level cluster bootstrap | 2026 |
| Probabilistic sensitivity analysis | 2026 |

All cross-validation is **grouped at the patient level**: sixteen patients have bilateral disease and contribute two records each, and their records are never split between training and test folds.

Written for Python 3.12 with `scikit-learn`, `pandas`, `numpy`, `scipy`, `matplotlib` and `openpyxl`.

## Turkish–English glossary

| Turkish | English |
|---|---|
| kohort, kayıt, hasta | cohort, record, patient |
| uygun popülasyon | eligible population (HR+/HER2−) |
| işaretlenen / uyumlu / uyumsuz | flagged / concordant / discordant |
| eşik | threshold |
| katlama, iç içe, tam iç içe | fold, nested, fully nested |
| imputasyon, eksiklik | imputation, missingness |
| ayırt gücü, kalibrasyon eğimi | discrimination, calibration slope |
| net fayda, karar eğrisi | net benefit, decision curve |
| duyarlılık / özgüllük | sensitivity / specificity |
| maliyet, tasarruf, ödeyici | cost, saving, payer |
| taban senaryo | base case |
| yapısal senaryo | structural scenario |
| uyum | adherence |
| ömür boyu, indirim oranı | lifetime, discount rate |
| kabul edilebilirlik eğrisi | acceptability curve |
| alt tip kuralı | subtype rule (S1) |
| öngörü kuralı | prediction rule (S2) |

## A note on the economic parameters

Two utility inputs were corrected late in the preparation of the manuscript and the values in this repository are the corrected ones: the locoregional recurrence utility is 0.73 (mapped to state R of the source meta-analysis) and the febrile neutropenia disutility is 0.150 (the value reported in the source). Capecitabine is excluded from the base case and examined as a scenario. Running `04_economics/markov_base.py` unchanged reproduces the main-pipeline figures: 499,580 / 495,316 / 492,185 TRY and 14.88344 / 14.88509 / 14.88635 QALYs for S0 / S1 / S2, a saving of 7,395 TRY against the cohort-observed pathway and 3,131 TRY against the subtype rule. The manuscript reports the fully nested flag count of 113 as primary, which `06_erpor/fn_economics.py` computes from the same model: 491,835 TRY and 14.88649 QALYs for S2, a saving of 7,745 TRY and 3,480 TRY respectively. Differences are computed before rounding, so they do not always equal the difference of the rounded figures.

## Licence

MIT.

## Citation

Please cite the article. Until it appears, cite this repository by its Zenodo concept identifier, https://doi.org/10.5281/zenodo.23051985, which always resolves to the most recent archived version.

## Running the scripts

Place the cohort file in `data/` and run each script from its own directory, for example:

```bash
cd code/03_validation && python tam_icice.py
cd code/06_erpor && python fn_subtype_split.py
```

Scripts read the cohort from `../../data/` and write results to `../../outputs/` and `../../figures/`.

Three input files are expected in `data/` and are not distributed here: the raw cohort
(`Kagan_TEZ.xlsx`), the corrected analysis cohort written by `01_cohort/kohort_duzeltme.py`
(`Kohort_v17.xlsx`), and the two national price lists used by `04_economics/ilac_maliyet.py`.
The price lists are public documents published by the Turkish Medicines and Medical Devices
Agency and the Social Security Institution; the cohort files are not shareable.
