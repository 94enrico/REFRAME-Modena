# Compare BMI and REFRAME at all visits, a random visit, the first visit and the last visit per period.

DATA_P <- "0_OriginalData"

# Main tables folder
TABS_P <- "3_Tables"
R_P <- "1_Script/R"

# Output folder for this analysis
AGREEMENT_P <- file.path(TABS_P, "Agreement_results")
LATEX_P <- file.path(AGREEMENT_P, "latex")
DATA_OUT_P <- file.path(AGREEMENT_P, "csv_xlsx")

# Source helper
source(file.path(R_P, "0_Helper_agreement.R"))

# Packages
library(dplyr)
library(tidyr)
library(tibble)
library(lubridate)
library(openxlsx)
library(kableExtra)

# User parameters
B_BOOT <- 2000
SEED_MASTER <- 20260419

# Output folders
dir.create(AGREEMENT_P, showWarnings = FALSE, recursive = TRUE)
dir.create(LATEX_P, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_OUT_P, showWarnings = FALSE, recursive = TRUE)

# Loading data
load(file = file.path(DATA_P, "REFRAME_agreement_data.RData"))
df_msm <- df_agreement_ready

# Data preparation

# Robust recoding of sex:
# female = 0, male = 1
df_msm <- df_msm %>%
  mutate(
    sex = case_when(
      sex %in% c("F", "Female", "female", 0, "0") ~ 0,
      sex %in% c("M", "Male", "male", 1, "1") ~ 1,
      TRUE ~ NA_real_
    )
  )

# Drop patients with inconsistent sex across visits
ids_inconsistent_sex <- df_msm %>%
  filter(!is.na(sex)) %>%
  group_by(patient_id) %>%
  summarise(n_sex = n_distinct(sex), .groups = "drop") %>%
  filter(n_sex > 1) %>%
  pull(patient_id)

if (length(ids_inconsistent_sex) > 0) {
  df_msm <- df_msm %>%
    filter(!patient_id %in% ids_inconsistent_sex)
}

# Mild access_date coercion
if (!inherits(df_msm$access_date, "Date")) {
  if (inherits(df_msm$access_date, c("POSIXct", "POSIXlt"))) {
    df_msm <- df_msm %>%
      mutate(access_date = as.Date(access_date))
  } else if (is.numeric(df_msm$access_date)) {
    df_msm <- df_msm %>%
      mutate(access_date = as.Date(access_date, origin = "1899-12-30"))
  } else {
    df_msm <- df_msm %>%
      mutate(access_date = as.Date(access_date))
  }
}

df_agreement <- df_msm %>%
  filter(state != 4) %>%
  select(patient_id, access_date, access_age, state, bmi_cat, sex) %>%
  arrange(patient_id, access_date, access_age) %>%
  group_by(patient_id) %>%
  mutate(visit_n = row_number()) %>%
  ungroup() %>%
  mutate(
    period = make_period(access_date)
  )

# Missingness diagnostics
missing_counts <- colSums(is.na(df_agreement))
print(missing_counts)

# General complete-case dataset for analyses 1-3
df_agreement_nona <- df_agreement %>%
  drop_na(patient_id, access_date, state, bmi_cat, sex)

# Period-complete dataset for analysis 4 only
df_agreement_period_nona <- df_agreement_nona %>%
  drop_na(period)

cat(
  "Rows removed due to missingness for analyses 1-3:",
  nrow(df_agreement) - nrow(df_agreement_nona), "\n"
)
cat(
  "Patients lost due to missingness for analyses 1-3:",
  length(unique(df_agreement$patient_id)) - length(unique(df_agreement_nona$patient_id)),
  "\n"
)
cat("Final number of rows for analyses 1-3:", nrow(df_agreement_nona), "\n")
cat("Final number of patients for analyses 1-3:", n_distinct(df_agreement_nona$patient_id), "\n")

cat(
  "Rows excluded because historical period is missing:",
  nrow(df_agreement_nona) - nrow(df_agreement_period_nona), "\n"
)
cat(
  "Patients excluded from period-specific analysis because historical period is missing:",
  length(unique(df_agreement_nona$patient_id)) - length(unique(df_agreement_period_nona$patient_id)),
  "\n"
)
cat("Final number of rows for period-specific analysis:", nrow(df_agreement_period_nona), "\n")
cat("Final number of patients for period-specific analysis:", n_distinct(df_agreement_period_nona$patient_id), "\n")

print(
  df_agreement_nona %>%
    distinct(patient_id, sex) %>%
    count(sex)
)

print(
  df_agreement_period_nona %>%
    count(period)
)

# Basic confusion matrices on the full all-visits dataset
write_confusion_matrix_tex(
  data = df_agreement_nona,
  file_path = file.path(LATEX_P, "confusion_matrix_overall.tex"),
  caption = "Confusion matrix between obesity state and BMI category (overall population, all visits)."
)

write_confusion_matrix_tex(
  data = df_agreement_nona %>% filter(sex == 1),
  file_path = file.path(LATEX_P, "confusion_matrix_male.tex"),
  caption = "Confusion matrix between obesity state and BMI category (men, all visits)."
)

write_confusion_matrix_tex(
  data = df_agreement_nona %>% filter(sex == 0),
  file_path = file.path(LATEX_P, "confusion_matrix_female.tex"),
  caption = "Confusion matrix between obesity state and BMI category (women, all visits)."
)

# Requested analysis datasets
set.seed(SEED_MASTER)

# 1) All visits
analysis_all_visits <- df_agreement_nona

# 2) Random visit per patient:
#    point estimate will use one random visit per patient,
#    and bootstrap will re-do the random visit selection inside each replicate
analysis_random_visit_source <- df_agreement_nona

# 3) First visit per patient
analysis_first_visit <- df_agreement_nona %>%
  arrange(patient_id, access_date, access_age, visit_n) %>%
  group_by(patient_id) %>%
  slice_head(n = 1) %>%
  ungroup()

# 4) Historical periods: last visit within each period for each patient
analysis_period_last_within_period <- df_agreement_period_nona %>%
  arrange(patient_id, period, access_date, access_age, visit_n) %>%
  group_by(patient_id, period) %>%
  slice_tail(n = 1) %>%
  ungroup()

# Helper to save outputs
save_setting_outputs <- function(setting_results,
                                 base_name,
                                 latex_caption_estimates,
                                 latex_caption_pairwise = NULL,
                                 latex_caption_kappa_pairwise = NULL) {
  xlsx_path <- file.path(DATA_OUT_P, paste0(base_name, ".xlsx"))
  csv_long_path <- file.path(DATA_OUT_P, paste0(base_name, "_estimates_long.csv"))
  csv_wide_path <- file.path(DATA_OUT_P, paste0(base_name, "_estimates_wide.csv"))

  csv_diff_long_path <- if (!is.null(setting_results$pairwise_long)) {
    file.path(DATA_OUT_P, paste0(base_name, "_pairwise_diff_long.csv"))
  } else {
    NULL
  }

  csv_diff_wide_path <- if (!is.null(setting_results$pairwise_wide)) {
    file.path(DATA_OUT_P, paste0(base_name, "_pairwise_diff_wide.csv"))
  } else {
    NULL
  }

  write_summary_workbook(
    long_df = setting_results$estimates_long,
    wide_df = setting_results$estimates_wide,
    xlsx_path = xlsx_path,
    csv_path_long = csv_long_path,
    csv_path_wide = csv_wide_path,
    diff_long_df = setting_results$pairwise_long,
    diff_wide_df = setting_results$pairwise_wide,
    csv_path_diff_long = csv_diff_long_path,
    csv_path_diff_wide = csv_diff_wide_path
  )

  write_latex_table(
    data = setting_results$estimates_wide,
    file_path = file.path(LATEX_P, paste0(base_name, "_summary.tex")),
    caption = latex_caption_estimates
  )

  if (!is.null(setting_results$pairwise_wide)) {
    write_latex_table(
      data = setting_results$pairwise_wide,
      file_path = file.path(LATEX_P, paste0(base_name, "_pairwise_diff_summary.tex")),
      caption = latex_caption_pairwise
    )

    kappa_pairwise_wide <- make_kappa_pairwise_wide(setting_results$pairwise_long)

    utils::write.csv(
      kappa_pairwise_wide,
      file = file.path(DATA_OUT_P, paste0(base_name, "_kappa_pairwise_diff.csv")),
      row.names = FALSE
    )

    write_latex_table(
      data = kappa_pairwise_wide,
      file_path = file.path(LATEX_P, paste0(base_name, "_kappa_pairwise_diff.tex")),
      caption = latex_caption_kappa_pairwise
    )
  }
}

# Run the four requested analyses
results_all_visits <- run_standard_setting(
  data = analysis_all_visits,
  setting_name = "all_visits",
  analysis_type = "all_visits",
  B = B_BOOT,
  seed_base = 1000
)

results_random_visit <- run_standard_setting(
  data = analysis_random_visit_source,
  setting_name = "random_visit",
  analysis_type = "random_visit",
  B = B_BOOT,
  seed_base = 2000
)

results_first_visit <- run_standard_setting(
  data = analysis_first_visit,
  setting_name = "first_visit",
  analysis_type = "first_visit",
  B = B_BOOT,
  seed_base = 3000
)

results_period_last_within_period <- run_period_last_setting(
  data = analysis_period_last_within_period,
  setting_name = "period_last_within_period",
  B = B_BOOT,
  seed_base = 4000
)

# Save outputs
save_setting_outputs(
  setting_results = results_all_visits,
  base_name = "agreement_all_visits",
  latex_caption_estimates = paste(
    "Agreement measures between BMI category and obesity classification.",
    "All visits included.",
    "Patients are resampled as clusters.",
    "Bootstrap 95\\% confidence intervals."
  )
)

save_setting_outputs(
  setting_results = results_random_visit,
  base_name = "agreement_random_visit",
  latex_caption_estimates = paste(
    "Agreement measures between BMI category and obesity classification.",
    "One random visit per patient.",
    "In each bootstrap replicate, patients are resampled and one visit is sampled at random within each resampled patient.",
    "Bootstrap 95\\% confidence intervals."
  )
)

save_setting_outputs(
  setting_results = results_first_visit,
  base_name = "agreement_first_visit",
  latex_caption_estimates = paste(
    "Agreement measures between BMI category and obesity classification.",
    "First visit per patient only.",
    "Bootstrap 95\\% confidence intervals."
  )
)

save_setting_outputs(
  setting_results = results_period_last_within_period,
  base_name = "agreement_period_last_within_period",
  latex_caption_estimates = paste(
    "Agreement measures between BMI category and obesity classification by historical period.",
    "Within each historical period, the last visit per patient is retained.",
    "Bootstrap 95\\% confidence intervals."
  ),
  latex_caption_pairwise = paste(
    "Pairwise differences across historical periods in agreement measures between BMI category and obesity classification.",
    "Within each historical period, the last visit per patient is retained.",
    "Difference defined as period 2 minus period 1.",
    "Bootstrap 95\\% confidence intervals."
  ),
  latex_caption_kappa_pairwise = paste(
    "Pairwise differences in weighted kappa across historical periods.",
    "Within each historical period, the last visit per patient is retained.",
    "Difference defined as period 2 minus period 1.",
    "Bootstrap 95\\% confidence intervals."
  )
)
