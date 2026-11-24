# Original analysis

Run the numbered R scripts from this directory. They retain the original filenames, variable names and statistical specifications. The historical cohort includes visits from 2002. The historical hepatic thresholds are recorded in the repository README.

The original workbooks belong in `0_OriginalData/`. Use the historical copy of `MHMC_database_20251126 final 4.xlsx` and `QoL_obesity_Roma.xlsx`. The updated workbook used for peer review has the same MHMC filename and additional drug fields. The data preparation script writes both the R panel and a CSV that can be used for the review.

The scripts use R with tidyverse, readxl, lubridate, knitr, kableExtra, openxlsx, patchwork, msm and VIM.

```sh
Rscript 1_Script/R/1_DataPreparation.R
Rscript 1_Script/R/2_Descriptive_Statistics.R
Rscript 1_Script/R/3_Multistate_Model.R
Rscript 1_Script/R/4_Agreement_Analysis.R
Rscript 1_Script/R/5_Make_Tables_From_Descriptives.R
```

Script 3 fits the four nested manuscript models (`reframe_fit0` through `reframe_fit3`). Its archived exploratory quadratic-age model can be run with `--quadratic-age`. Script 4 compares BMI and REFRAME state, with bootstrap resampling by participant. Scripts 4 and 5 can run as soon as script 2 finishes.
