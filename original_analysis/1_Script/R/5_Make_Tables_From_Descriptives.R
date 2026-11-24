# Assemble the manuscript tables from the descriptive workbooks.

# Packages
library(tidyverse)
library(readxl)
library(openxlsx)
library(lubridate)

data_p <- "0_OriginalData"
out_dir <- "3_Tables/Descriptives_excel"

# Period definitions
period_levels <- c("2005-2014", "2015-2019", "2020+")

# Inputs

# Table 1 inputs
file_t1_overall <- file.path(
  out_dir,
  "multi_visit_only",
  "visits",
  "desc_visit_1.xlsx"
)

file_t1_state <- file.path(
  out_dir,
  "multi_visit_only",
  "visits_stratified",
  "desc_visit_1_by_state.xlsx"
)

# Table 3 inputs: last observation within historical periods
table3_overall_files <- c(
  "2005-2014" = file.path(
    out_dir,
    "multi_visit_only",
    "historical_periods_last",
    "desc_period_last_2005-2014.xlsx"
  ),
  "2015-2019" = file.path(
    out_dir,
    "multi_visit_only",
    "historical_periods_last",
    "desc_period_last_2015-2019.xlsx"
  ),
  "2020+" = file.path(
    out_dir,
    "multi_visit_only",
    "historical_periods_last",
    "desc_period_last_2020_.xlsx"
  )
)

table3_state_files <- c(
  "2005-2014" = file.path(
    out_dir,
    "multi_visit_only",
    "historical_periods_last_stratified",
    "desc_period_last_2005-2014_by_state.xlsx"
  ),
  "2015-2019" = file.path(
    out_dir,
    "multi_visit_only",
    "historical_periods_last_stratified",
    "desc_period_last_2015-2019_by_state.xlsx"
  ),
  "2020+" = file.path(
    out_dir,
    "multi_visit_only",
    "historical_periods_last_stratified",
    "desc_period_last_2020__by_state.xlsx"
  )
)

# Outputs (Excel only)
out_t1_xlsx <- file.path(out_dir, "table 1.xlsx")
out_t2_xlsx <- file.path(out_dir, "table 2.xlsx")
out_t3_xlsx <- file.path(out_dir, "table 3.xlsx")

# Variable specification
table_spec <- tribble(
  ~row_order, ~variable, ~label, ~type, ~yes_level,

  # Demographics
  1.0, "study_entry_age", "Age, years", "cont", NA_character_,
  2.0, "sex", "Female sex (yes)", "disc", "F",

  # HIV history
  3.0, "hiv_duration_yrs", "HIV duration, years", "cont", NA_character_,
  4.0, "Nadir_CD4__cells/uL", "Nadir CD4, cells/µL", "cont", NA_character_,
  5.0, "CD4_absolute__cells/uL", "CD4 absolute, cells/µL", "cont", NA_character_,
  6.0, "CD4/CD8__ratio", "CD4/CD8 ratio", "cont", NA_character_,
  7.0, "HIV_viral_load_undetectable", "Undetectable HIV viral load (yes)", "disc", "1",

  # Antiretroviral therapy
  8.0, "NRTI", "NRTI (yes)", "disc", "1",
  9.0, "NNRTI", "NNRTI (yes)", "disc", "1",
  10.0, "IP", "Protease inhibitors (yes)", "disc", "1",
  11.0, "INSTI", "INSTI (yes)", "disc", "1",
  12.0, "TAF", "TAF (yes)", "disc", "1",

  # Lifestyle / adiposity
  13.0, "smoking", "Smoking (yes)", "disc", "1",
  14.0, "BMI", "BMI, kg/m²", "cont", NA_character_,
  15.0, "waist_circumference__cm", "Waist circumference, cm", "cont", NA_character_,
  16.0, "WHtR", "Waist-to-height ratio", "cont", NA_character_,
  17.0, "VAT_dexa", "Visceral adipose tissue, cm²", "cont", NA_character_,

  # Comorbidities / clinical indicators
  18.0, "comorbidity_DM", "Diabetes (yes)", "disc", "1",
  19.0, "comorbidity_HTN", "Hypertension (yes)", "disc", "1",
  20.0, "comorbidity_Dislipidemia", "Dyslipidemia (yes)", "disc", "1",
  21.0, "liver_disease", "Liver disease (yes)", "disc", "1",
  22.0, "comorbidity_Kidn", "Chronic kidney disease (yes)", "disc", "1",
  23.0, "comorbidity_CVD", "Cardiovascular disease (yes)", "disc", "1",
  24.0, "phys_act_bin", "Physical activity (yes)", "disc", "1"
)

# Helpers
extract_p_from_string <- function(x) {
  x <- as.character(x)
  out <- stringr::str_extract(x, "p<1e-4|p=NA|p=\\d+\\.\\d+")
  out <- stringr::str_remove(out, "^p=")
  out
}

format_p_value <- function(p) {
  if (length(p) == 0 || is.null(p) || is.na(p)) {
    return(NA_character_)
  }
  if (p < 1e-4) {
    return("<1e-4")
  }
  sprintf("%.4f", p)
}

format_mean_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  sprintf("%.2f ± %.2f", mean(x), sd(x))
}

format_numden_yes <- function(x, yes_level) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  n_yes <- sum(x == yes_level)
  sprintf("%d/%d (%.2f%%)", n_yes, length(x), 100 * n_yes / length(x))
}

format_numden_level <- function(x, level) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  n_level <- sum(x == level)
  sprintf("%d/%d (%.2f%%)", n_level, length(x), 100 * n_level / length(x))
}

check_columns_present <- function(df, cols, df_name) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop(
      "Missing columns in ", df_name, ": ",
      paste(miss, collapse = ", ")
    )
  }
}

check_variable_values_present <- function(df, vars, df_name) {
  check_columns_present(df, "variable", df_name)

  miss <- setdiff(vars, unique(df$variable))
  if (length(miss) > 0) {
    stop(
      "Missing variables in ", df_name, ": ",
      paste(miss, collapse = ", ")
    )
  }
}

read_desc_workbook <- function(path) {
  list(
    continuous = read_excel(path, sheet = "continuous"),
    discrete   = read_excel(path, sheet = "discrete")
  )
}

build_prevalence_table <- function(df, var, levels, labels, periods) {
  out <- tibble(
    Variable = labels,
    level = as.character(levels)
  )

  out$Overall <- purrr::map_chr(
    out$level,
    ~ format_numden_level(df[[var]], .x)
  )

  for (p in periods) {
    df_p <- df %>% filter(period == p)
    out[[p]] <- purrr::map_chr(
      out$level,
      ~ format_numden_level(df_p[[var]], .x)
    )
  }

  out %>%
    select(-level)
}

save_single_sheet_xlsx <- function(df, path, sheet_name) {
  wb <- createWorkbook()
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, df)

  header_style <- createStyle(textDecoration = "bold")
  addStyle(
    wb,
    sheet = sheet_name,
    style = header_style,
    rows = 1,
    cols = 1:ncol(df),
    gridExpand = TRUE
  )

  setColWidths(wb, sheet_name, cols = 1:ncol(df), widths = "auto")
  saveWorkbook(wb, path, overwrite = TRUE)
}

build_table_from_workbooks <- function(cont_overall,
                                       disc_overall,
                                       cont_state,
                                       disc_state,
                                       spec,
                                       use_access_age = FALSE,
                                       period_label = NULL) {
  available_cont_vars <- union(
    unique(cont_overall$variable),
    unique(cont_state$variable)
  )

  age_var_for_workbook <- "study_entry_age"

  if (use_access_age) {
    if ("access_age" %in% available_cont_vars) {
      age_var_for_workbook <- "access_age"
    } else if ("study_end_age" %in% available_cont_vars) {
      age_var_for_workbook <- "study_end_age"
    } else if ("study_entry_age" %in% available_cont_vars) {
      age_var_for_workbook <- "study_entry_age"
    } else {
      stop("No age variable found in continuous workbook.")
    }
  }

  spec_use <- spec %>%
    mutate(
      variable = ifelse(
        variable == "study_entry_age",
        age_var_for_workbook,
        variable
      )
    )

  spec_cont <- spec_use %>% filter(type == "cont")
  spec_disc <- spec_use %>% filter(type == "disc")

  check_variable_values_present(
    cont_overall,
    spec_cont$variable,
    "continuous overall workbook"
  )

  check_variable_values_present(
    cont_state,
    spec_cont$variable,
    "continuous state workbook"
  )

  check_columns_present(
    disc_overall,
    c("variable", "level", "n", "perc_NA", "count", "prop"),
    "discrete overall workbook"
  )

  check_columns_present(
    disc_state,
    c(
      "variable", "level",
      "n_lv1", "perc_NA_lv1", "count_lv1", "prop_lv1",
      "n_lv2", "perc_NA_lv2", "count_lv2", "prop_lv2",
      "n_lv3", "perc_NA_lv3", "count_lv3", "prop_lv3",
      "prop_test"
    ),
    "discrete state workbook"
  )

  cont_overall_tab <- cont_overall %>%
    filter(variable %in% spec_use$variable) %>%
    select(variable, mean, sd) %>%
    mutate(Overall = sprintf("%.2f ± %.2f", mean, sd)) %>%
    select(variable, Overall)

  cont_state_tab <- cont_state %>%
    filter(variable %in% spec_use$variable) %>%
    select(variable, mean_lv1, sd_lv1, mean_lv2, sd_lv2, mean_lv3, sd_lv3, kw) %>%
    mutate(
      `Metabolically Healthy` = sprintf("%.2f ± %.2f", mean_lv1, sd_lv1),
      `Preclinical Obesity`   = sprintf("%.2f ± %.2f", mean_lv2, sd_lv2),
      `Clinical Obesity`      = sprintf("%.2f ± %.2f", mean_lv3, sd_lv3),
      `P value`               = extract_p_from_string(kw)
    ) %>%
    select(
      variable,
      `Metabolically Healthy`,
      `Preclinical Obesity`,
      `Clinical Obesity`,
      `P value`
    )

  disc_overall_tab <- disc_overall %>%
    mutate(level = as.character(level)) %>%
    inner_join(
      spec_disc %>% select(variable, yes_level),
      by = "variable"
    ) %>%
    filter(level == yes_level) %>%
    mutate(
      denom = round(n * (1 - perc_NA / 100)),
      Overall = sprintf("%d/%d (%.2f%%)", count, denom, 100 * prop)
    ) %>%
    select(variable, Overall)

  disc_state_tab <- disc_state %>%
    mutate(level = as.character(level)) %>%
    inner_join(
      spec_disc %>% select(variable, yes_level),
      by = "variable"
    ) %>%
    filter(level == yes_level) %>%
    select(
      variable,
      n_lv1, perc_NA_lv1, count_lv1, prop_lv1,
      n_lv2, perc_NA_lv2, count_lv2, prop_lv2,
      n_lv3, perc_NA_lv3, count_lv3, prop_lv3,
      prop_test
    ) %>%
    mutate(
      denom_lv1 = round(n_lv1 * (1 - perc_NA_lv1 / 100)),
      denom_lv2 = round(n_lv2 * (1 - perc_NA_lv2 / 100)),
      denom_lv3 = round(n_lv3 * (1 - perc_NA_lv3 / 100)),
      `Metabolically Healthy` = sprintf("%d/%d (%.2f%%)", count_lv1, denom_lv1, 100 * prop_lv1),
      `Preclinical Obesity` = sprintf("%d/%d (%.2f%%)", count_lv2, denom_lv2, 100 * prop_lv2),
      `Clinical Obesity` = sprintf("%d/%d (%.2f%%)", count_lv3, denom_lv3, 100 * prop_lv3),
      `P value` = extract_p_from_string(prop_test)
    ) %>%
    select(
      variable,
      `Metabolically Healthy`,
      `Preclinical Obesity`,
      `Clinical Obesity`,
      `P value`
    )

  out_cont <- spec_use %>%
    filter(type == "cont") %>%
    select(row_order, variable, label) %>%
    left_join(cont_overall_tab, by = "variable") %>%
    left_join(cont_state_tab, by = "variable")

  out_disc <- spec_use %>%
    filter(type == "disc") %>%
    select(row_order, variable, label) %>%
    left_join(disc_overall_tab, by = "variable") %>%
    left_join(disc_state_tab, by = "variable")

  out <- bind_rows(out_cont, out_disc) %>%
    arrange(row_order)

  if (is.null(period_label)) {
    out <- out %>%
      select(
        Variable = label,
        Overall,
        `Metabolically Healthy`,
        `Preclinical Obesity`,
        `Clinical Obesity`,
        `P value`
      )
  } else {
    out <- out %>%
      mutate(Period = period_label) %>%
      select(
        Period,
        Variable = label,
        Overall,
        `Metabolically Healthy`,
        `Preclinical Obesity`,
        `Clinical Obesity`,
        `P value`
      )
  }

  out
}

# TABLE 1
t1_overall <- read_desc_workbook(file_t1_overall)
t1_state <- read_desc_workbook(file_t1_state)

table1 <- build_table_from_workbooks(
  cont_overall   = t1_overall$continuous,
  disc_overall   = t1_overall$discrete,
  cont_state     = t1_state$continuous,
  disc_state     = t1_state$discrete,
  spec           = table_spec,
  use_access_age = FALSE,
  period_label   = NULL
)

print(table1, n = Inf)
save_single_sheet_xlsx(table1, out_t1_xlsx, "Table 1")

# TABLE 2
load(file.path(data_p, "REFRAME_data.RData"))

cols_drop <- c(
  "birth_date", "death_date", "death_age_exact",
  "HIV_positivity_date", "state_123"
)

msm_data <- msm_data %>%
  select(-any_of(cols_drop)) %>%
  filter(sex %in% c("F", "M")) %>%
  mutate(access_date = as.Date(access_date))

ids_multi <- msm_data %>%
  count(patient_id) %>%
  filter(n > 1) %>%
  pull(patient_id)

df_multi <- msm_data %>%
  filter(patient_id %in% ids_multi)

df_multi <- df_multi %>%
  filter(state %in% c(1, 2, 3))

df_multi <- df_multi %>%
  mutate(
    period = case_when(
      access_date >= as.Date("2005-01-01") & access_date < as.Date("2015-01-01") ~ "2005-2014",
      access_date >= as.Date("2015-01-01") & access_date < as.Date("2020-01-01") ~ "2015-2019",
      access_date >= as.Date("2020-01-01") ~ "2020+",
      TRUE ~ NA_character_
    ),
    period = factor(period, levels = period_levels)
  ) %>%
  filter(!is.na(period), !is.na(access_date))

df_last_period <- df_multi %>%
  arrange(patient_id, access_date) %>%
  group_by(patient_id, period) %>%
  slice_max(access_date, n = 1, with_ties = FALSE) %>%
  ungroup()

table2_spec <- table_spec %>%
  mutate(variable = ifelse(variable == "study_entry_age", "access_age", variable))

check_columns_present(df_last_period, table2_spec$variable, "df_last_period")

cont_tab <- table2_spec %>%
  filter(type == "cont") %>%
  mutate(
    Overall = NA_character_,
    `2005-2014` = NA_character_,
    `2015-2019` = NA_character_,
    `2020+` = NA_character_,
    `P value` = NA_character_
  )

for (i in seq_len(nrow(cont_tab))) {
  v <- cont_tab$variable[i]

  x_all <- df_last_period[[v]]
  x_05014 <- df_last_period[df_last_period$period == "2005-2014", v, drop = TRUE]
  x_1519 <- df_last_period[df_last_period$period == "2015-2019", v, drop = TRUE]
  x_2020 <- df_last_period[df_last_period$period == "2020+", v, drop = TRUE]

  cont_tab$Overall[i] <- format_mean_sd(x_all)
  cont_tab$`2005-2014`[i] <- format_mean_sd(x_05014)
  cont_tab$`2015-2019`[i] <- format_mean_sd(x_1519)
  cont_tab$`2020+`[i] <- format_mean_sd(x_2020)

  d_test <- df_last_period %>%
    transmute(period = period, x = .data[[v]]) %>%
    filter(!is.na(period), !is.na(x))

  if (n_distinct(d_test$period) >= 2 && nrow(d_test) >= 3) {
    kw <- tryCatch(
      kruskal.test(x ~ period, data = d_test),
      error = function(e) NULL
    )

    if (!is.null(kw)) {
      cont_tab$`P value`[i] <- format_p_value(kw$p.value)
    }
  }
}

disc_tab <- table2_spec %>%
  filter(type == "disc") %>%
  mutate(
    Overall = NA_character_,
    `2005-2014` = NA_character_,
    `2015-2019` = NA_character_,
    `2020+` = NA_character_,
    `P value` = NA_character_
  )

for (i in seq_len(nrow(disc_tab))) {
  v <- disc_tab$variable[i]
  yes <- disc_tab$yes_level[i]

  x_all <- df_last_period[[v]]
  x_05014 <- df_last_period[df_last_period$period == "2005-2014", v, drop = TRUE]
  x_1519 <- df_last_period[df_last_period$period == "2015-2019", v, drop = TRUE]
  x_2020 <- df_last_period[df_last_period$period == "2020+", v, drop = TRUE]

  disc_tab$Overall[i] <- format_numden_yes(x_all, yes)
  disc_tab$`2005-2014`[i] <- format_numden_yes(x_05014, yes)
  disc_tab$`2015-2019`[i] <- format_numden_yes(x_1519, yes)
  disc_tab$`2020+`[i] <- format_numden_yes(x_2020, yes)

  d_test <- df_last_period %>%
    transmute(period = period, x = as.character(.data[[v]])) %>%
    filter(!is.na(period), !is.na(x))

  if (nrow(d_test) > 0 && n_distinct(d_test$period) >= 2 && n_distinct(d_test$x) >= 2) {
    tab_test <- table(d_test$period, d_test$x)

    test_out <- tryCatch(
      suppressWarnings(chisq.test(tab_test)),
      error = function(e) NULL
    )

    if (!is.null(test_out)) {
      disc_tab$`P value`[i] <- format_p_value(test_out$p.value)
    } else {
      test_out2 <- tryCatch(
        fisher.test(tab_test),
        error = function(e) NULL
      )

      if (!is.null(test_out2)) {
        disc_tab$`P value`[i] <- format_p_value(test_out2$p.value)
      }
    }
  }
}

table2 <- bind_rows(
  cont_tab %>% select(row_order, label, Overall, `2005-2014`, `2015-2019`, `2020+`, `P value`),
  disc_tab %>% select(row_order, label, Overall, `2005-2014`, `2015-2019`, `2020+`, `P value`)
) %>%
  arrange(row_order) %>%
  select(
    Variable = label,
    Overall,
    `2005-2014`,
    `2015-2019`,
    `2020+`,
    `P value`
  )

print(table2, n = Inf)

bmi_prev <- build_prevalence_table(
  df      = df_last_period,
  var     = "bmi_cat",
  levels  = c("1", "2", "3"),
  labels  = c("BMI category 1", "BMI category 2", "BMI category 3"),
  periods = period_levels
)

print(bmi_prev, n = Inf)

state_prev <- build_prevalence_table(
  df      = df_last_period %>% filter(state %in% c(1, 2, 3)),
  var     = "state",
  levels  = c("1", "2", "3"),
  labels  = c("Metabolically Healthy", "Preclinical Obesity", "Clinical Obesity"),
  periods = period_levels
)

print(state_prev, n = Inf)

wb_t2 <- createWorkbook()

addWorksheet(wb_t2, "Table 2")
writeData(wb_t2, "Table 2", table2)

addWorksheet(wb_t2, "BMI prevalence")
writeData(wb_t2, "BMI prevalence", bmi_prev)

addWorksheet(wb_t2, "State prevalence")
writeData(wb_t2, "State prevalence", state_prev)

header_style <- createStyle(textDecoration = "bold")

addStyle(wb_t2, "Table 2", header_style, rows = 1, cols = 1:ncol(table2), gridExpand = TRUE)
addStyle(wb_t2, "BMI prevalence", header_style, rows = 1, cols = 1:ncol(bmi_prev), gridExpand = TRUE)
addStyle(wb_t2, "State prevalence", header_style, rows = 1, cols = 1:ncol(state_prev), gridExpand = TRUE)

setColWidths(wb_t2, "Table 2", cols = 1:ncol(table2), widths = "auto")
setColWidths(wb_t2, "BMI prevalence", cols = 1:ncol(bmi_prev), widths = "auto")
setColWidths(wb_t2, "State prevalence", cols = 1:ncol(state_prev), widths = "auto")

saveWorkbook(wb_t2, out_t2_xlsx, overwrite = TRUE)

# TABLE 3
# Same variables as Table 1, but last observation within each
# historical period, with columns Overall + obesity state.
# Rows are triplicated by period, and a Period column is added.
table3_list <- vector("list", length(period_levels))
names(table3_list) <- period_levels

for (p in period_levels) {
  wb_overall <- read_desc_workbook(table3_overall_files[[p]])
  wb_state <- read_desc_workbook(table3_state_files[[p]])

  table3_list[[p]] <- build_table_from_workbooks(
    cont_overall   = wb_overall$continuous,
    disc_overall   = wb_overall$discrete,
    cont_state     = wb_state$continuous,
    disc_state     = wb_state$discrete,
    spec           = table_spec,
    use_access_age = TRUE,
    period_label   = p
  )
}

table3 <- bind_rows(table3_list)

print(table3, n = Inf)
save_single_sheet_xlsx(table3, out_t3_xlsx, "Table 3")
