# Produce the descriptive tables and inputs for the agreement and multistate analyses.

data_p <- "0_OriginalData"
tabs_p <- "3_Tables"
dir.create(tabs_p, recursive = TRUE, showWarnings = FALSE)

figs_p <- "2_Figures"
R_p <- "1_Script/R"

# Packages
library(tidyverse)
library(lubridate)
library(openxlsx)
library(tidyverse)
library(readxl)
library(knitr)
library(kableExtra)

source(file.path(R_p, "0_Helper_descriptive_statistics.R"))
source(file.path(R_p, "0_Helper_msm_diagnostics.R"))

# Hyperparameters
study_end_date <- as.Date("2025-11-26")

# Thresholds
lower_dt <- 30 / 365.25 # 30 days in years
upper_dt <- 3 # 3 years

# Loading data
load(file.path(data_p, "REFRAME_data.RData"))

# Drop raw date columns + related intermediates
cols_drop <- c(
  "birth_date", "death_date", "death_age_exact",
  "HIV_positivity_date", "state_123"
)

msm_data <- msm_data %>%
  select(-any_of(cols_drop))
colnames(msm_data)

# Restrict gender to male/female categories (only one patient falls outside them)
msm_data <- msm_data %>%
  filter(sex %in% c("F", "M"))

# Study entry age
msm_data <- msm_data %>%
  group_by(patient_id) %>%
  mutate(
    study_entry_age = first(access_age)
  ) %>%
  ungroup()

# Time intervals and pruning

# Split dataset
visits_per_patient <- msm_data %>%
  dplyr::count(patient_id, name = "n_visits")

single_ids <- visits_per_patient %>%
  filter(n_visits <= 1) %>%
  pull(patient_id)

longitudinal_ids <- visits_per_patient %>%
  filter(n_visits >= 2) %>%
  pull(patient_id)

msm_data_single <- msm_data %>%
  filter(patient_id %in% single_ids)

msm_data_long <- msm_data %>%
  filter(patient_id %in% longitudinal_ids)

# Patients with single visit
(n_single_patients <- n_distinct(msm_data_single$patient_id))

# Build intervals
df_intervals <- build_intervals(
  msm_data_long,
  id = patient_id,
  time = access_age,
  state = state
)

# Extract intervals to trim
intervals_short <- df_intervals %>%
  filter(dt < lower_dt)

intervals_long <- df_intervals %>%
  filter(dt > upper_dt)

# Save trimmed observations
write.csv(
  intervals_short,
  file = file.path(tabs_p, "intervals_under_30days.csv"),
  row.names = FALSE
)

write.csv(
  intervals_long,
  file = file.path(tabs_p, "intervals_over_3years.csv"),
  row.names = FALSE
)

# Sequential trimming
msm_data_trimmed <- trim_panel_sequential(
  msm_data_long,
  id = patient_id,
  time = access_age,
  lower = lower_dt,
  upper = upper_dt,
  truncate_after_upper = TRUE
)

# Visits statistics
visits_long <- msm_data_long %>%
  count(patient_id, name = "n_visits")

visits_trim <- msm_data_trimmed %>%
  count(patient_id, name = "n_visits")

visits_stats_long <- visits_long %>%
  summarise(
    mean = mean(n_visits),
    sd   = sd(n_visits),
    min  = min(n_visits),
    p10  = quantile(n_visits, .10),
    q1   = quantile(n_visits, .25),
    med  = median(n_visits),
    q3   = quantile(n_visits, .75),
    p90  = quantile(n_visits, .90),
    max  = max(n_visits)
  )

visits_stats_trim <- visits_trim %>%
  summarise(
    mean = mean(n_visits),
    sd   = sd(n_visits),
    min  = min(n_visits),
    p10  = quantile(n_visits, .10),
    q1   = quantile(n_visits, .25),
    med  = median(n_visits),
    q3   = quantile(n_visits, .75),
    p90  = quantile(n_visits, .90),
    max  = max(n_visits)
  )

tab_visits <- bind_rows(
  visits_stats_long %>% mutate(dataset = "Longitudinal"),
  visits_stats_trim %>% mutate(dataset = "Trimmed")
) %>%
  select(dataset, everything())

tab_visits_latex <- kbl(
  tab_visits,
  format = "latex",
  booktabs = TRUE,
  digits = 3,
  linesep = "",
  caption = "Number of visits per patient",
  label = "visits_summary"
)

cat(
  tab_visits_latex,
  file = file.path(tabs_p, "tab_visits_summary.tex")
)

# Follow-up statistics
follow_long <- msm_data_long %>%
  group_by(patient_id) %>%
  summarise(
    followup = max(access_age) - min(access_age),
    .groups = "drop"
  )

follow_trim <- msm_data_trimmed %>%
  group_by(patient_id) %>%
  summarise(
    followup = max(access_age) - min(access_age),
    .groups = "drop"
  )

follow_stats_long <- follow_long %>%
  summarise(
    mean = mean(followup),
    sd   = sd(followup),
    min  = min(followup),
    p10  = quantile(followup, .10),
    q1   = quantile(followup, .25),
    med  = median(followup),
    q3   = quantile(followup, .75),
    p90  = quantile(followup, .90),
    max  = max(followup)
  )

follow_stats_trim <- follow_trim %>%
  summarise(
    mean = mean(followup),
    sd   = sd(followup),
    min  = min(followup),
    p10  = quantile(followup, .10),
    q1   = quantile(followup, .25),
    med  = median(followup),
    q3   = quantile(followup, .75),
    p90  = quantile(followup, .90),
    max  = max(followup)
  )

tab_follow <- bind_rows(
  follow_stats_long %>% mutate(dataset = "Longitudinal"),
  follow_stats_trim %>% mutate(dataset = "Trimmed")
) %>%
  select(dataset, everything())

tab_follow_latex <- kbl(
  tab_follow,
  format = "latex",
  booktabs = TRUE,
  digits = 3,
  linesep = "",
  caption = "Patient follow-up (years)",
  label = "followup_summary"
)

cat(
  tab_follow_latex,
  file = file.path(tabs_p, "tab_followup_summary.tex")
)

# Interval statistics
interval_stats_long <- df_intervals %>%
  summarise(
    mean = mean(dt),
    sd   = sd(dt),
    min  = min(dt),
    p10  = quantile(dt, .10),
    q1   = quantile(dt, .25),
    med  = median(dt),
    q3   = quantile(dt, .75),
    p90  = quantile(dt, .90),
    max  = max(dt)
  )

df_intervals_trim <- build_intervals(
  msm_data_trimmed,
  id = patient_id,
  time = access_age,
  state = state
)

interval_stats_trim <- df_intervals_trim %>%
  summarise(
    mean = mean(dt),
    sd   = sd(dt),
    min  = min(dt),
    p10  = quantile(dt, .10),
    q1   = quantile(dt, .25),
    med  = median(dt),
    q3   = quantile(dt, .75),
    p90  = quantile(dt, .90),
    max  = max(dt)
  )

tab_intervals <- bind_rows(
  interval_stats_long %>% mutate(dataset = "Longitudinal"),
  interval_stats_trim %>% mutate(dataset = "Trimmed")
) %>%
  select(dataset, everything())

tab_intervals_latex <- kbl(
  tab_intervals,
  format = "latex",
  booktabs = TRUE,
  digits = 3,
  linesep = "",
  caption = "Inter-visit intervals (years)",
  label = "interval_summary"
)

cat(
  tab_intervals_latex,
  file = file.path(tabs_p, "tab_interval_summary.tex")
)

# Rebuild dataset for exploratory analysis
msm_data_exploratory <- bind_rows(
  msm_data_single,
  msm_data_trimmed
) %>%
  arrange(patient_id, access_age)

# Missing data diagnostics
vars_state <- c("state", "bmi_cat")

vars_covariates <- c(
  "sex",
  "INSTI",
  "TAF",
  "IP",
  "hiv_duration_yrs",
  "nadir_cd4",
  "phys_act_bin",
  "diff_cal",
  "eq5d5l_bin",
  "comorbidity_Dislipidemia",
  "liver_disease",
  "insulin_resistance",
  "sppb_low"
)

vars_state_def <- setdiff(
  colnames(msm_data_exploratory),
  c(vars_state, vars_covariates)
)

vars_miss <- c(vars_state, vars_covariates, vars_state_def)

# Dataset ≥2 visits
msm_data_2vis <- msm_data_exploratory %>%
  group_by(patient_id) %>%
  filter(n() >= 2) %>%
  ungroup()

### First visit datasets
first_visit_full <- msm_data_exploratory %>%
  arrange(patient_id, access_age) %>%
  group_by(patient_id) %>%
  slice_min(access_age, n = 1, with_ties = FALSE) %>%
  ungroup()

first_visit_2vis <- msm_data_2vis %>%
  arrange(patient_id, access_age) %>%
  group_by(patient_id) %>%
  slice_min(access_age, n = 1, with_ties = FALSE) %>%
  ungroup()

# Full dataset
missing_first_full <- build_missing(first_visit_full)
missing_rows_full <- build_missing(msm_data_exploratory)

# ≥2 visits dataset
missing_first_2vis <- build_missing(first_visit_2vis)
missing_rows_2vis <- build_missing(msm_data_2vis)

# Merge tables
missing_first <- merge_missing(missing_first_full, missing_first_2vis)
missing_rows <- merge_missing(missing_rows_full, missing_rows_2vis)

# Table — cohort entry
tab_missing_first <- missing_first %>%
  select(
    Group = group,
    Variable = variable,
    `Missing full (n, \\%)` = missing_full,
    `Missing ≥2 visits (n, \\%)` = missing_2vis
  )

tab_missing_first_latex <- kbl(
  tab_missing_first,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Missing values at cohort entry (first visit per patient)",
  label = "missing_first_visit"
)

cat(
  tab_missing_first_latex,
  file = file.path(tabs_p, "tab_missing_first_visit.tex")
)

# Table — longitudinal dataset
tab_missing_rows <- missing_rows %>%
  select(
    Group = group,
    Variable = variable,
    `Missing full (n, \\%)` = missing_full,
    `Missing ≥2 visits (n, \\%)` = missing_2vis
  )

tab_missing_rows_latex <- kbl(
  tab_missing_rows,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Missing values in the longitudinal dataset",
  label = "missing_longitudinal"
)

cat(
  tab_missing_rows_latex,
  file = file.path(tabs_p, "tab_missing_longitudinal.tex")
)

# Variable selection for descriptive analyses

# All candidate variables, excluding ID and raw date variables
vars_desc <- names(msm_data_exploratory) %>%
  setdiff("patient_id") %>%
  discard(~ str_ends(.x, "date"))

# Variables to force as discrete
vars_discrete_forced <- c(
  # Outcomes / states
  "state",
  "bmi_cat",

  # Demographics / lifestyle
  "sex",
  "smoking",

  # HIV / therapy
  "HIV_viral_load_undetectable",
  "NRTI",
  "NNRTI",
  "IP",
  "INSTI",
  "TAF",

  # Comorbidities / derived indicators
  "comorbidity_COPD",
  "comorbidity_DM",
  "comorbidity_HTN",
  "comorbidity_Kidn",
  "comorbidity_CVD",
  "comorbidity_Dislipidemia",
  "phys_act_bin",
  "liver_disease",
  "insulin_resistance",
  "sppb_low"
)

# Variables to force as continuous
vars_cont_forced <- c(
  "Nadir_CD4__cells/uL",
  "CD4_absolute__cells/uL",
  "CD4/CD8__ratio"
)

# Automatic type-based classification
vars_discrete_auto <- vars_desc %>%
  keep(~ is.character(msm_data_exploratory[[.x]]) ||
    is.logical(msm_data_exploratory[[.x]]))

vars_cont_auto <- vars_desc %>%
  discard(~ is.character(msm_data_exploratory[[.x]]) ||
    is.logical(msm_data_exploratory[[.x]]))

# Final discrete / continuous sets
vars_discrete <- union(vars_discrete_forced, vars_discrete_auto)
vars_cont <- union(vars_cont_forced, vars_cont_auto)

# Ensure no overlap and remove time variable from descriptives
vars_cont <- setdiff(vars_cont, vars_discrete)
vars_cont <- setdiff(vars_cont, "access_age")

# Output folders
base_out <- file.path(tabs_p, "Descriptives_excel")
ensure_dir(base_out)

base_all <- ensure_dir(file.path(base_out, "all_patients"))
base_multi <- ensure_dir(file.path(base_out, "multi_visit_only"))

create_tree <- function(root) {
  list(
    visits = ensure_dir(file.path(root, "visits")),
    visits_strat = ensure_dir(file.path(root, "visits_stratified")),
    visits_by_sex = ensure_dir(file.path(root, "visits_by_sex")),
    years = ensure_dir(file.path(root, "calendar_year")),
    years_strat = ensure_dir(file.path(root, "calendar_year_stratified")),
    periods = ensure_dir(file.path(root, "historical_periods")),
    periods_strat = ensure_dir(file.path(root, "historical_periods_stratified")),
    periods_first = ensure_dir(file.path(root, "historical_periods_first")),
    periods_first_strat = ensure_dir(file.path(root, "historical_periods_first_stratified")),
    periods_last = ensure_dir(file.path(root, "historical_periods_last")),
    periods_last_strat = ensure_dir(file.path(root, "historical_periods_last_stratified")),
    sex = ensure_dir(file.path(root, "sex")),
    sex_strat = ensure_dir(file.path(root, "sex_stratified"))
  )
}

out_all <- create_tree(base_all)
out_multi <- create_tree(base_multi)

# Base filters
df_base_all <- msm_data_exploratory %>%
  mutate(access_date = as.Date(access_date)) %>%
  filter(!is.na(access_date))

ids_multi <- df_base_all %>%
  count(patient_id) %>%
  filter(n > 1) %>%
  pull(patient_id)

df_base_multi <- df_base_all %>%
  filter(patient_id %in% ids_multi)

# Sanity checks
cat("Patients (all):   ", n_distinct(df_base_all$patient_id), "\n")
cat("Patients (multi): ", n_distinct(df_base_multi$patient_id), "\n\n")

cat("Rows (all):   ", nrow(df_base_all), "\n")
cat("Rows (multi): ", nrow(df_base_multi), "\n\n")

# MASTER RUNNER
run_all_descriptives <- function(df_base, out) {
  cat("Running descriptives...\n")

  # Visits
  cat(" -> Visits\n")

  df_visit <- df_base %>%
    add_visit_index(id = "patient_id", date_var = "access_date")

  run_descriptive_pipeline(
    df = df_visit,
    group_var = "visit_index",
    vars_cont = vars_cont,
    vars_disc = vars_discrete,
    output_dir = out$visits,
    file_prefix = "desc_visit",
    stratify_by = FALSE,
    id = "patient_id",
    strat_var = "state"
  )

  run_descriptive_pipeline(
    df = df_visit,
    group_var = "visit_index",
    vars_cont = vars_cont,
    vars_disc = vars_discrete,
    output_dir = out$visits_strat,
    file_prefix = "desc_visit",
    stratify_by = TRUE,
    id = "patient_id",
    strat_var = "state"
  )

  # Visits stratified by sex
  cat(" -> Visits (stratified by sex)\n")

  df_visit_sex <- df_visit %>%
    mutate(
      sex = case_when(
        sex %in% c("F", 0) ~ "Female",
        sex %in% c("M", 1) ~ "Male",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(sex))

  run_descriptive_pipeline(
    df = df_visit_sex,
    group_var = "visit_index",
    vars_cont = vars_cont,
    vars_disc = vars_discrete,
    output_dir = out$visits_by_sex,
    file_prefix = "desc_visit_sex",
    stratify_by = TRUE,
    id = "patient_id",
    strat_var = "sex"
  )

  # Calendar year
  cat(" -> Calendar year\n")

  df_year <- df_base %>%
    add_calendar_year(date_var = "access_date")

  run_descriptive_pipeline(
    df = df_year,
    group_var = "year",
    vars_cont = vars_cont,
    vars_disc = vars_discrete,
    output_dir = out$years,
    file_prefix = "desc_year",
    stratify_by = FALSE,
    id = "patient_id",
    strat_var = "state"
  )

  run_descriptive_pipeline(
    df = df_year,
    group_var = "year",
    vars_cont = vars_cont,
    vars_disc = vars_discrete,
    output_dir = out$years_strat,
    file_prefix = "desc_year",
    stratify_by = TRUE,
    id = "patient_id",
    strat_var = "state"
  )

  # Sex
  cat(" -> Sex\n")

  df_sex <- df_base %>%
    mutate(
      sex = case_when(
        sex %in% c("F", 0) ~ "Female",
        sex %in% c("M", 1) ~ "Male",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(sex))

  run_descriptive_pipeline(
    df = df_sex,
    group_var = "sex",
    vars_cont = vars_cont,
    vars_disc = vars_discrete,
    output_dir = out$sex,
    file_prefix = "desc_sex",
    stratify_by = FALSE,
    id = "patient_id",
    strat_var = "state"
  )

  # Historical periods
  cat(" -> Historical periods\n")

  df_period <- df_base %>%
    add_period(date_var = "access_date") %>%
    filter(!is.na(period))

  run_descriptive_pipeline(
    df = df_period,
    group_var = "period",
    vars_cont = vars_cont,
    vars_disc = vars_discrete,
    output_dir = out$periods,
    file_prefix = "desc_period",
    stratify_by = FALSE,
    id = "patient_id",
    strat_var = "state"
  )

  run_descriptive_pipeline(
    df = df_period,
    group_var = "period",
    vars_cont = vars_cont,
    vars_disc = vars_discrete,
    output_dir = out$periods_strat,
    file_prefix = "desc_period",
    stratify_by = TRUE,
    id = "patient_id",
    strat_var = "state"
  )

  # First visit per patient within period
  cat(" -> Historical periods (first visit per patient)\n")

  df_period_first <- df_period %>%
    arrange(patient_id, access_date) %>%
    group_by(patient_id, period) %>%
    slice_min(access_date, n = 1, with_ties = FALSE) %>%
    ungroup()

  run_descriptive_pipeline(
    df = df_period_first,
    group_var = "period",
    vars_cont = vars_cont,
    vars_disc = vars_discrete,
    output_dir = out$periods_first,
    file_prefix = "desc_period_first",
    stratify_by = FALSE,
    id = "patient_id",
    strat_var = "state"
  )

  run_descriptive_pipeline(
    df = df_period_first,
    group_var = "period",
    vars_cont = vars_cont,
    vars_disc = vars_discrete,
    output_dir = out$periods_first_strat,
    file_prefix = "desc_period_first",
    stratify_by = TRUE,
    id = "patient_id",
    strat_var = "state"
  )

  # Last visit per patient within period
  cat(" -> Historical periods (last visit per patient)\n")

  df_period_last <- df_period %>%
    arrange(patient_id, access_date) %>%
    group_by(patient_id, period) %>%
    slice_max(access_date, n = 1, with_ties = FALSE) %>%
    ungroup()

  run_descriptive_pipeline(
    df = df_period_last,
    group_var = "period",
    vars_cont = vars_cont,
    vars_disc = vars_discrete,
    output_dir = out$periods_last,
    file_prefix = "desc_period_last",
    stratify_by = FALSE,
    id = "patient_id",
    strat_var = "state"
  )

  run_descriptive_pipeline(
    df = df_period_last,
    group_var = "period",
    vars_cont = vars_cont,
    vars_disc = vars_discrete,
    output_dir = out$periods_last_strat,
    file_prefix = "desc_period_last",
    stratify_by = TRUE,
    id = "patient_id",
    strat_var = "state"
  )

  cat("Done.\n\n")
}

# Run both populations
cat("DESCRIPTIVES — ALL PATIENTS\n")
run_all_descriptives(df_base_all, out_all)

cat("DESCRIPTIVES — MULTI VISIT ONLY\n")
run_all_descriptives(df_base_multi, out_multi)

# MSM dataset: variable selection and ordering
vars_msm <- c(
  # ID
  "patient_id",

  # Time
  "access_age",

  # State
  "state",

  # State definition
  "bmi_cat",

  # Demographics
  "sex",

  # HIV / therapy
  "INSTI",
  "TAF",
  "IP",
  "hiv_duration_yrs",
  "nadir_cd4",

  # Lifestyle / metabolic
  "phys_act_bin",
  "diff_cal",

  # QoL
  "eq5d5l_bin",

  # Comorbidities
  "comorbidity_Dislipidemia",
  "liver_disease",
  "insulin_resistance",

  # Functional status
  "sppb_low"
)

df_msm <- msm_data_exploratory %>%
  select(any_of(vars_msm)) %>%
  arrange(patient_id, access_age)

# Save MSM-ready dataset
save(df_msm, file = file.path(data_p, "REFRAME_data_clean.RData"))

write.csv(
  df_msm,
  file = file.path(data_p, "REFRAME_data_clean.csv"),
  row.names = FALSE
)

# Save agreement-analysis dataset
df_agreement_ready <- msm_data_exploratory %>%
  select(any_of(c("patient_id", "access_age", "access_date", "state", "bmi_cat", "sex"))) %>%
  arrange(patient_id, access_age)

save(df_agreement_ready, file = file.path(data_p, "REFRAME_agreement_data.RData"))
write.csv(df_agreement_ready, file = file.path(data_p, "REFRAME_agreement_data.csv"), row.names = FALSE)
