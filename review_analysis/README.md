# Peer-review analyses

Run the commands from `review_analysis/`. Directory and output names match the review report. The revised manuscript uses adult visits in 2005–2025 and the corrected hepatic criteria. The model comparison also uses the historical cohort.

## Private inputs

The study team's `Dati/visite.csv.gz` and `Dati/dexa.csv.gz` are the analysis inputs. To rebuild them, use the historical visit panel from the original analysis and the updated clinical workbook with both INSTI generation fields.

```sh
python Paper/Scripts/prepare_inputs.py ../original_analysis/0_OriginalData/REFRAME_data.csv "/path/to/updated/MHMC_database_20251126 final 4.xlsx"
```

Model fitting uses `Dati/valori_iniziali_msm.rds` when supplied. Otherwise the scripts derive the same starting quantities from the four historical fits in `../original_analysis/6_Models/`. A team member with the saved review fits can go directly to the commands that export results.

Python dependencies are NumPy, pandas ≥2, SciPy, matplotlib, threadpoolctl and openpyxl. The model scripts use R with msm and VIM. Set `PYTHON` if the interpreter on your path is named something other than `python`.

## Revised manuscript

```sh
python Paper/Scripts/descrittive.py
python Paper/Scripts/concordanza.py
for model in 0 1 2 3; do Rscript Paper/Scripts/stima_msm.R "$model"; done
Rscript Paper/Scripts/risultati_msm.R
python Paper/Scripts/figure.py
```

## Reviewer requests

| Question | Script | Result |
|---|---|---|
| E1 | `E1_prevalenza/Scripts/prevalenze.py` | Crude, conditional and age/sex-standardised prevalence, using a random or the last visit per period. |
| R2.2 | `R2_2_eta/Scripts/confronto.py` | All three pairwise differences in clinical obesity prevalence. Run E1 first. |
| R2.1 and R3.4 | `R2_1_R3_4_tre_periodi/Scripts/presenze.py` and `sottogruppi.py` | Period attendance, baseline characteristics and period prevalence in 515 participants with continuous follow-up and 990 allowing re-entry. |
| R3.1 | `R3_1_BMI/Scripts/bmi_normale.py` | REFRAME states among participants with normal BMI. |
| R3.2 | `R3_2_DEXA/Scripts/dexa.py` | Paired VAT and total body fat, agreement, baseline reclassification and a figure of repeated examinations. |
| R3.9 | `R3_9_HIV/Scripts/hiv.py` | Full manuscript Tables 1 and 2, including CD4/CD8 and viral suppression. |
| R3.10 | `R3_10_menopausa/Scripts/menopausa.py` | Menopausal status and missingness at the visit. |
| R3.14 | `R3_14_INSTI/Scripts/insti.py` | INSTI generation at baseline and its denominator. |
| Sex-specific analysis | `Analisi_per_sesso/Scripts/per_sesso.py` | Baseline descriptions and observed transitions by sex. |
| E3e | `E3e_multistato/Scripts/follow_up.py` | Follow-up intervals and imputation counts. |

Run each Python script with `python`. E1 and the two three-period subgroup analyses fit their bootstrap models on the first run and reuse them on subsequent runs. In each analysis the standardisation population gives one third of the weight to each period.

### Four multistate models (E3e)

```sh
for model in M1 M2 M3 M4; do Rscript E3e_multistato/Scripts/stima.R "$model"; done
Rscript E3e_multistato/Scripts/risultati.R
```

| Model | Cohort | Effect of sex on mortality |
|---|---|---|
| M1 | Historical, from 2002 | Estimated separately for each origin state. |
| M2 | Historical, from 2002 | Common to all three death transitions. |
| M3 | Methods, 2005–2025 | Separate by origin, with the preclinical obesity to death HR fixed at 1. |
| M4 | Methods, 2005–2025 | Common to all three death transitions. |

Each analysis writes to its own `Tables/`, `Figures/` or `Modelli/` directory, as applicable.
