# Build the historical visit and death panel from the original workbooks.

data_p <- "0_OriginalData"
tabs_p <- "3_Tables"
dir.create(tabs_p, recursive = TRUE, showWarnings = FALSE)

figs_p <- "2_Figures"
R_p <- "1_Script/R"

# Packages

# Data preparation
library(tidyverse)
library(readxl)
library(lubridate)

# Exploratory analysis
library(msm)
library(knitr)
library(kableExtra)

source(file.path(R_p, "0_Helper_merge.R"))

# Hyperparameters
study_end_date <- as.Date("2025-11-26")
tolerance_days <- 30

# Main data import
df <- read_excel(
  file.path(data_p, "MHMC_database_20251126 final 4.xlsx"),
  na = c("NULL", "")
)

# Fix duplicated columns
df <- df %>%
  select(-c(comorbidity_CVD...391, WHtR...388)) %>%
  rename(
    WHtR = WHtR...20,
    comorbidity_CVD = comorbidity_CVD...387
  )

# QoL data import
df_QoL <- read_excel(
  file.path(data_p, "QoL_obesity_Roma.xlsx"),
  na = c("NULL", "")
)

# Date handling + exact ages
df <- df %>%
  mutate(
    across(c(
      birth_date, access_date, death_date,
      antropo_date, nutritional_date, IPAQ_date,
      DEXA_date, fibroscan_echosens_date, exams_date,
      HIV_positivity_date
    ), as.Date),
    access_age_exact =
      as.numeric(difftime(access_date, birth_date, units = "days")) / 365.25,
    death_age_exact =
      as.numeric(difftime(death_date, birth_date, units = "days")) / 365.25,
    study_end_age =
      as.numeric(difftime(study_end_date, birth_date, units = "days")) / 365.25
  )

# Helper: binary harmonization
to_logical01 <- function(x) {
  case_when(
    x %in% c(TRUE, "TRUE", 1, "1", "SI") ~ 1,
    x %in% c(FALSE, "FALSE", 0, "0", "NO") ~ 0,
    TRUE ~ NA
  )
}

df <- df %>%
  mutate(
    sex = as.character(sex),

    # Smoking -> binary yes/no
    smoking = case_when(
      smoking %in% c("Intense (>10day)", "Moderate (<=10day)") ~ "SI",
      smoking %in% c("No") ~ "NO",
      TRUE ~ NA_character_
    ),
    across(any_of(c(
      # HIV / therapy / smoking
      "HIV_viral_load_undetectable",
      "NRTI",
      "IP",
      "NNRTI",
      "INSTI",
      "TAF",
      "smoking",

      # Comorbidities
      "comorbidity_COPD",
      "comorbidity_DM",
      "comorbidity_HTN",
      "comorbidity_Dislipidemia",
      "comorbidity_Kidn",
      "comorbidity_CVD"
    )), to_logical01)
  )

# Base panel: one row per patient per access_date
vis_base <- df %>%
  filter(
    !is.na(patient_id),
    !is.na(access_date),
    !is.na(access_age_exact)
  ) %>%
  transmute(
    patient_id,
    sex,
    birth_date,
    access_date,
    access_age = access_age_exact,
    death_date,
    death_age_exact,
    study_end_age,

    # HIV history
    `Nadir_CD4__cells/uL`,
    `CD4_absolute__cells/uL`,
    `CD4/CD8__ratio`,
    HIV_viral_load_undetectable,
    HIV_positivity_date,

    # Therapy (visit-based)
    NRTI,
    IP,
    NNRTI,
    INSTI,
    TAF,

    # Smoking
    smoking,

    # Metabolic / physical performance
    HOMA,
    SPPB = sppb_score,

    # Chronic diseases / organ-damage components
    comorbidity_COPD,
    comorbidity_DM,
    comorbidity_HTN,
    comorbidity_Kidn,
    comorbidity_CVD,
    comorbidity_Dislipidemia
  ) %>%
  arrange(patient_id, access_date) %>%
  distinct(patient_id, access_date, .keep_all = TRUE)

# Anthropometry (fallback to access_date)
anthro_tbl <- df %>%
  transmute(
    patient_id,
    access_date,
    meas_date = antropo_date,
    WHtR,
    waist_circumference__cm,
    BMI
  )

# Diet
diet_tbl <- df %>%
  transmute(
    patient_id,
    access_date,
    meas_date = nutritional_date,
    diet__kcal
  )

# Physical activity (IPAQ)
ipaq_tbl <- df %>%
  transmute(
    patient_id,
    access_date,
    meas_date = IPAQ_date,
    phys_act_bin = case_when(
      # 1) precedence to physical_activity
      !is.na(physical_activity) & physical_activity == "No" ~ 0L,
      !is.na(physical_activity) ~ 1L,

      # 2) fallback to ipaq_result
      is.na(physical_activity) & !is.na(ipaq_result) & ipaq_result == "LOW" ~ 0L,
      is.na(physical_activity) & !is.na(ipaq_result) ~ 1L,

      # 3) both missing
      TRUE ~ NA_integer_
    )
  )

# DEXA
dexa_tbl <- df %>%
  transmute(
    patient_id,
    access_date,
    meas_date = DEXA_date,
    VAT_dexa,

    # BMR estimation
    BMR_kcal_day = case_when(
      # 1) Priority: DEXA lean mass (Katch–McArdle)
      !is.na(dexa_total_lean_mass__g) ~
        370 + 21.6 * (dexa_total_lean_mass__g / 1000),

      # 2) Fallback: Mifflin–St Jeor (only if ALL needed)
      is.na(dexa_total_lean_mass__g) &
        !is.na(weight__kg) &
        !is.na(height__cm) &
        !is.na(access_age) &
        !is.na(sex) ~
        case_when(
          sex %in% c("M", "Male", "male") ~
            10 * weight__kg + 6.25 * height__cm - 5 * access_age + 5,
          sex %in% c("F", "Female", "female") ~
            10 * weight__kg + 6.25 * height__cm - 5 * access_age - 161,
          TRUE ~ NA_real_
        ),

      # 3) Otherwise: missing
      TRUE ~ NA_real_
    )
  )

# Liver / Fibroscan (fallback to exams_date)
liver_tbl <- df %>%
  transmute(
    patient_id,
    access_date,
    meas_date = coalesce(fibroscan_echosens_date, exams_date),
    cap_med,
    e_med,
    FEGMIL,
    FIB4
  )

# QoL (external dataset)
qol_tbl <- df_QoL %>%
  transmute(
    patient_id = IdPaz,
    meas_date  = as.Date(DataEvento),
    eq5d5l_bin = if_else(eq5d5l_score > 0.8, 1L, 0L, missing = NA_integer_)
  )

# Merge covariate blocks to panel (LOCF with 30-day forward tolerance)

# Anthropometry (visit-based, fallback on access_date)
vis_anthro <- merge_locf_tolerance(
  panel_df = vis_base,
  meas_df = anthro_tbl,
  tolerance_days = tolerance_days,
  type = "visit"
)

# Diet (visit-based, fallback on access_date)
vis_diet <- merge_locf_tolerance(
  panel_df = vis_base,
  meas_df = diet_tbl,
  tolerance_days = tolerance_days,
  type = "visit"
)

# Physical activity (IPAQ; visit-based)
vis_ipaq <- merge_locf_tolerance(
  panel_df = vis_base,
  meas_df = ipaq_tbl,
  tolerance_days = tolerance_days,
  type = "visit"
)

# DEXA (visit-based)
vis_dexa <- merge_locf_tolerance(
  panel_df = vis_base,
  meas_df = dexa_tbl,
  tolerance_days = tolerance_days,
  type = "visit"
)

# Liver / Fibroscan (visit-based, fallback already handled upstream)
vis_liver <- merge_locf_tolerance(
  panel_df = vis_base,
  meas_df = liver_tbl,
  tolerance_days = tolerance_days,
  type = "visit"
)

# QoL (external event-based dataset, NO access_date)
vis_qol <- merge_locf_tolerance(
  panel_df = vis_base,
  meas_df = qol_tbl,
  tolerance_days = tolerance_days,
  type = "event"
)

# Merge all covariate blocks into final panel
vis_final <- vis_base %>%
  left_join(
    vis_anthro %>%
      select(patient_id, access_date, WHtR, waist_circumference__cm, BMI),
    by = c("patient_id", "access_date")
  ) %>%
  left_join(
    vis_diet %>%
      select(patient_id, access_date, diet__kcal),
    by = c("patient_id", "access_date")
  ) %>%
  left_join(
    vis_ipaq %>%
      select(patient_id, access_date, phys_act_bin),
    by = c("patient_id", "access_date")
  ) %>%
  left_join(
    vis_dexa %>%
      select(patient_id, access_date, VAT_dexa, BMR_kcal_day),
    by = c("patient_id", "access_date")
  ) %>%
  left_join(
    vis_liver %>%
      select(patient_id, access_date, cap_med, e_med, FEGMIL, FIB4),
    by = c("patient_id", "access_date")
  ) %>%
  left_join(
    vis_qol %>%
      select(patient_id, access_date, eq5d5l_bin),
    by = c("patient_id", "access_date")
  )

# Indicators for state construction + BMI categorization

# Helper: waist circumference threshold by sex
wc_threshold <- function(sex_vec) {
  if_else(sex_vec %in% c("M", "Male", "male"), 102, 88)
}

vis_ind <- vis_final %>%
  mutate(
    # HIV duration (years)
    hiv_duration_yrs =
      as.numeric(difftime(access_date, HIV_positivity_date, units = "days")) / 365.25,

    # BMI categorical classification (1..3)
    bmi_cat = case_when(
      is.na(BMI) ~ NA_integer_,
      BMI < 25 ~ 1L,
      BMI < 30 ~ 2L,
      BMI >= 30 ~ 3L
    ),

    # Adiposity criteria
    wc_thr = wc_threshold(sex),
    crit_vat = !is.na(VAT_dexa) & VAT_dexa >= 150,
    crit_whtr = !is.na(WHtR) & WHtR >= 0.50,
    crit_wc = !is.na(waist_circumference__cm) &
      waist_circumference__cm >= wc_thr,

    # Adiposity undefined if:
    # VAT is NA AND (WHtR is NA OR WC is NA)
    adiposity_undef =
      is.na(VAT_dexa) & (is.na(WHtR) | is.na(waist_circumference__cm)),
    excess_adiposity = case_when(
      adiposity_undef ~ NA,
      crit_vat ~ TRUE,
      (!is.na(WHtR) & !is.na(waist_circumference__cm)) ~ (crit_whtr & crit_wc),
      TRUE ~ NA
    ),

    # Liver disease (for descriptive purposes)
    # Definition:
    # Historical thresholds. The revision uses stiffness >= 8 and FIB-4 > 1.45.
    liver_disease = case_when(
      !is.na(cap_med) & !is.na(e_med) &
        cap_med >= 248 & e_med < 8 ~ 1L,
      !is.na(FEGMIL) & !is.na(FIB4) &
        FEGMIL < 1.1 & FIB4 <= 1.45 ~ 1L,
      (!is.na(cap_med) & !is.na(e_med)) |
        (!is.na(FEGMIL) & !is.na(FIB4)) ~ 0L,
      TRUE ~ NA_integer_
    ),

    # Organ damage (for state definition)
    # If ALL components are NA -> organ_damage = NA
    # Otherwise, TRUE if at least one component is present
    any_component_observed = !(
      is.na(comorbidity_COPD) &
        is.na(comorbidity_DM) &
        is.na(comorbidity_HTN) &
        is.na(comorbidity_Kidn) &
        is.na(comorbidity_CVD) &
        is.na(liver_disease)),
    organ_damage = case_when(
      !any_component_observed ~ NA,
      (coalesce(comorbidity_COPD, 0) == 1) |
        (coalesce(comorbidity_DM, 0) == 1) |
        (coalesce(comorbidity_HTN, 0) == 1) |
        (coalesce(comorbidity_Kidn, 0) == 1) |
        (coalesce(comorbidity_CVD, 0) == 1) |
        (coalesce(liver_disease, 0) == 1) ~ TRUE,
      TRUE ~ FALSE
    ),

    # State definition (1–3)
    state_123 = case_when(
      is.na(excess_adiposity) ~ NA_integer_,
      is.na(organ_damage) ~ NA_integer_,
      excess_adiposity == FALSE ~ 1L,
      excess_adiposity == TRUE & organ_damage == FALSE ~ 2L,
      excess_adiposity == TRUE & organ_damage == TRUE ~ 3L,
      TRUE ~ NA_integer_
    ),

    # Caloric balance
    # diff_cal = intake - expenditure (BMR only)
    diff_cal = if_else(
      is.na(diet__kcal) | is.na(BMR_kcal_day),
      NA_real_,
      diet__kcal - BMR_kcal_day
    ),

    # Descriptive-only indicators
    insulin_resistance = case_when(
      is.na(HOMA) ~ NA_integer_,
      HOMA >= 2.5 ~ 1L,
      TRUE ~ 0L
    ),
    sppb_low = case_when(
      is.na(SPPB) ~ NA_integer_,
      SPPB < 10 ~ 1L,
      TRUE ~ 0L
    )
  ) %>%
  select(
    -wc_thr, -crit_vat, -crit_whtr, -crit_wc,
    -adiposity_undef, -any_component_observed
  )

# NA Diagnostics for the state

# Total rows and patients
(n_rows_total <- nrow(vis_ind))
(n_patients_total <- vis_ind %>%
  distinct(patient_id) %>%
  nrow())

# Rows with NA state
(n_rows_na <- vis_ind %>%
  summarise(n = sum(is.na(state_123))) %>%
  pull(n))

(perc_rows_na <- round(100 * n_rows_na / n_rows_total, 2))

# Patients with at least one NA state
(n_patients_na <- vis_ind %>%
  filter(is.na(state_123)) %>%
  distinct(patient_id) %>%
  nrow())

(perc_patients_na <- round(100 * n_patients_na / n_patients_total, 2))

# Save NA observations
na_state_df <- vis_ind %>%
  filter(is.na(state_123)) %>%
  select(patient_id, access_date) %>%
  arrange(patient_id, access_date)

write_csv(
  na_state_df,
  file.path(tabs_p, "obesity_NA_visits.csv")
)

# Remove NA rows
panel_visits <- vis_ind %>%
  filter(!is.na(state_123)) %>%
  arrange(patient_id, access_age)

nrow(panel_visits)
panel_visits %>%
  distinct(patient_id) %>%
  nrow()

# Multistate dataset construction (MSM-style, keep ALL variables)

# Last visit per patient
last_visit <- panel_visits %>%
  group_by(patient_id) %>%
  slice_max(access_age, n = 1, with_ties = FALSE) %>%
  ungroup()

# Death rows
death_rows <- last_visit %>%
  filter(!is.na(death_date) & !is.na(death_age_exact)) %>%
  mutate(
    access_age = death_age_exact,
    state = 4L
  )

# Panel rows (states 1–3)
panel_rows <- panel_visits %>%
  mutate(
    state = state_123
  )

# Final MSM-ready dataset
msm_data <- bind_rows(panel_rows, death_rows) %>%
  arrange(patient_id, access_age) %>%
  group_by(patient_id) %>%
  distinct(access_age, .keep_all = TRUE) %>%
  ungroup()

# Save dataset
save(
  msm_data,
  file = file.path(data_p, "REFRAME_data.RData")
)
# The review links the historical panel to the updated workbook by person and visit date.
write.csv(msm_data, file.path(data_p, "REFRAME_data.csv"), row.names = FALSE)
