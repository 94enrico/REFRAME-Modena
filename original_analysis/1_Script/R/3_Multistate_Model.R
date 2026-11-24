# Fit the original multistate models with attained age as the time scale.

data_p <- "0_OriginalData"
tabs_p <- "3_Tables"
figs_p <- "2_Figures"
R_p <- "1_Script/R"
expl_p <- "5_Exploratory_Analysis"
modl_p <- "6_Models"
for (path in c(tabs_p, figs_p, expl_p, modl_p)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# Packages
library(tidyverse)
library(readxl)
library(lubridate)
library(knitr)
library(kableExtra)
library(msm)
library(VIM)

source(file.path(R_p, "0_Helper_msm_diagnostics.R"))

# Loading data
load(file = file.path(data_p, "REFRAME_data_clean.RData"))
str(df_msm)

# Keeping only patients with >=2 visits

# Save original dataset
df_msm_full <- df_msm

df_msm <- df_msm %>%
  group_by(patient_id) %>%
  mutate(n_visits = n()) %>%
  ungroup() %>%
  filter(n_visits >= 2) %>%
  select(-n_visits)

# Intervals

# Intervals
df_intervals <- build_intervals(df_msm,
  id = patient_id,
  time = access_age, state = state
)

# Stats on intervals
df_intervals %>%
  summarise(
    n_gaps = n(),
    mean = mean(dt),
    sd = sd(dt),
    min = min(dt),
    p5 = quantile(dt, 0.05),
    p10 = quantile(dt, 0.10),
    q1 = quantile(dt, 0.25),
    med = median(dt),
    q3 = quantile(dt, 0.75),
    p90 = quantile(dt, 0.90),
    p95 = quantile(dt, 0.95),
    p97 = quantile(dt, 0.97),
    p99 = quantile(dt, 0.99),
    max = max(dt)
  )

# Missing covariates
colSums(is.na(df_msm))
colSums(is.na(df_msm)) / nrow(df_msm)

# Drop variables with excessive missingness
df_msm <- df_msm %>%
  select(-diff_cal, -eq5d5l_bin, -bmi_cat)

# Binary clinical indicators: missing -> absence
df_msm <- df_msm %>%
  mutate(
    liver_disease = replace_na(liver_disease, 0),
    insulin_resistance = replace_na(insulin_resistance, 0),
    sppb_low = replace_na(sppb_low, 0)
  )

# HIV duration: fix value at study entry
df_msm <- df_msm %>%
  arrange(patient_id, access_age) %>%
  group_by(patient_id) %>%
  mutate(
    hiv_duration_yrs = {
      if (all(is.na(hiv_duration_yrs))) {
        hiv_duration_yrs
      } else {
        hiv_duration_yrs[which(!is.na(hiv_duration_yrs))[1]]
      }
    }
  ) %>%
  ungroup()

# Hot-deck / KNN imputations
# Physical activity
df_msm <- kNN(
  df_msm,
  variable = "phys_act_bin",
  k = 5,
  imp_var = FALSE
)

# HIV duration (remaining missing)
df_msm <- kNN(
  df_msm,
  variable = "hiv_duration_yrs",
  k = 5,
  imp_var = FALSE
)

# Check
colSums(is.na(df_msm))
colSums(is.na(df_msm)) / nrow(df_msm)

# Setting up the dataset for the analysis

# Dataset
# The descriptive cohort includes participants recorded as female or male.
df_msm <- df_msm %>%
  mutate(
    sex = case_when(
      sex == "F" ~ 0,
      sex == "M" ~ 1,
      TRUE ~ NA_real_
    )
  )

str(df_msm)

# Baseline HIV duration and tertiles
df_msm <- df_msm %>%
  arrange(patient_id, access_age) %>%
  group_by(patient_id) %>%
  mutate(hiv_dur_baseline_ter = first(hiv_duration_yrs)) %>%
  ungroup()

# Compute tertile cutpoints
cuts <- quantile(
  df_msm$hiv_dur_baseline_ter,
  probs = c(0, 1 / 3, 2 / 3, 1),
  na.rm = TRUE
)

df_msm <- df_msm %>%
  mutate(
    hiv_tertile = cut(
      hiv_dur_baseline_ter,
      breaks = cuts,
      include.lowest = TRUE,
      labels = c("Short", "Medium", "Long")
    )
  ) %>%
  select(-hiv_dur_baseline_ter)

# Transition tables

# State table
statetable.msm(state, patient_id, data = df_msm)

# State tables by covariate
covariates <- setdiff(
  names(df_msm),
  c("patient_id", "access_age", "state", "hiv_duration_yrs")
)

state_tables <- list()

for (v in covariates) {
  cat("Covariate:", v, "\n")

  tab <- statetable_by_interval(
    df_msm,
    id = patient_id,
    state = state,
    time = access_age,
    covariate = !!sym(v)
  )$statetable

  print(tab)

  state_tables[[v]] <- tab
}

# HIV duration
plot(density(df_msm$hiv_duration_yrs, na.rm = TRUE))
statetable_by_interval(df_msm,
  id = patient_id, state = state,
  time = access_age, covariate = hiv_tertile
)

# Controlling risk sets

# Transition by age class
age_cuts <- c(0, 40, 65, 120)
df_intervals <- build_intervals(df_msm,
  id = patient_id,
  time = access_age, state = state
)

diag <- transition_diagnostics(df_intervals, age_cuts,
  min_exposure = 20,
  min_events = 5,
  min_eff_exp = 8
)

print(diag, n = 12)

# Transition diagnostics by covariate
dir.create(expl_p, showWarnings = FALSE, recursive = TRUE)

# Covariates (ORDER = PDF NUMBER)
covariates <- list(
  list(name = "sex", type = "discrete"),
  list(name = "phys_act_bin", type = "discrete"),
  list(name = "INSTI", type = "discrete"),
  list(name = "TAF", type = "discrete"),
  list(name = "IP", type = "discrete"),
  list(name = "hiv_duration_yrs", type = "continuous", bins = 5),
  list(name = "hiv_tertile", type = "discrete"),
  list(name = "comorbidity_Dislipidemia", type = "discrete"),
  list(name = "liver_disease", type = "discrete"),
  list(name = "insulin_resistance", type = "discrete"),
  list(name = "sppb_low", type = "discrete")
)

# Extract transitions
states <- get_states(df_intervals)

transitions <- expand_grid(
  from = states,
  to   = states
) %>%
  filter(from != to) %>%
  filter(from != 6) %>%
  filter(!(from == 5 & to %in% 1:4))

print(transitions, n = 30)

# MAIN LOOP
for (i in seq_len(nrow(transitions))) {
  from_i <- transitions$from[i]
  to_i <- transitions$to[i]

  cat("\n====================================\n")
  cat("Processing transition:", from_i, "→", to_i, "\n")
  cat("====================================\n")

  # Create folder
  trans_folder <- file.path(
    expl_p,
    paste0("transition_", from_i, "_", to_i)
  )

  dir.create(trans_folder,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # TXT summary
  write_transition_txt(
    trans_folder,
    from_i,
    to_i,
    diag,
    df_msm
  )

  # 0 — DIAGNOSTIC BASELINE PLOT
  p_diag <- plot_transition_baseline(
    df_intervals,
    from = from_i,
    to = to_i,
    age_cuts = age_cuts,
    bandwidth = 5,
    zoom_empirical = TRUE
  )

  if (!is.null(p_diag)) {
    ggsave(
      filename = paste0("0_trans", from_i, to_i, "_diagnostic.pdf"),
      plot = p_diag,
      path = trans_folder,
      width = 7,
      height = 5
    )
  }

  # LOOP COVARIATES
  for (j in seq_along(covariates)) {
    cv <- covariates[[j]]
    cv_name <- cv$name

    cat("Covariate:", cv_name, "\n")

    if (cv$type == "baseline") next

    bins <- cv$bins %||% 4

    p <- plot_transition_baseline_cov(
      df_intervals,
      from_i,
      to_i,
      age_cuts,
      covariate = cv_name,
      cov_type = cv$type,
      cov_bins = bins
    )

    if (!is.null(p)) {
      ggsave(
        filename = paste0(j, "_trans", from_i, to_i, "_", cv_name, ".pdf"),
        plot = p,
        path = trans_folder,
        width = 9,
        height = 6
      )
    }
  }
}

# Warning expeted.

# State table overall
trans_tab <- statetable.msm(state, patient_id, data = df_msm) %>%
  as.data.frame.matrix()

trans_tab <- tibble::rownames_to_column(trans_tab, "From")

tab_state <- kbl(
  trans_tab,
  format = "latex",
  booktabs = TRUE,
  digits = 0,
  linesep = "",
  caption = "Observed transition counts in the CKM multistate dataset"
) %>%
  kable_styling(
    latex_options = "hold_position",
    font_size = 9
  )

tab_state <- sub(
  "\\\\caption\\{",
  "\\\\caption{Observed transition counts in the CKM multistate dataset}\\\\label{tab:state_table}\n",
  tab_state
)

cat(
  tab_state,
  file = file.path(tabs_p, "tab_state_table.tex")
)

# State tables by selected covariates
cov_list <- c("liver_disease", "sppb_low", "TAF")

for (cov in cov_list) {
  cov_tab <- statetable_by_interval(
    df_msm,
    id = patient_id,
    state = state,
    time = access_age,
    covariate = !!sym(cov)
  )$statetable

  for (lvl in dimnames(cov_tab)$cov0) {
    tab_lvl <- cov_tab[, , lvl] %>%
      as.data.frame.matrix() %>%
      tibble::rownames_to_column("From")

    core <- kbl(
      tab_lvl,
      format = "latex",
      booktabs = TRUE,
      digits = 0,
      linesep = ""
    )

    caption_txt <- paste0(
      "Observed transition counts by ", cov,
      " (level = ", lvl, ")"
    )

    latex_tab <- paste0(
      "\\begin{table}[!htbp]
\\centering
\\caption{", caption_txt, "}
\\label{tab:", cov, "_", lvl, "}

\\fontsize{9}{11}\\selectfont
",
      core,
      "
\\end{table}"
    )

    cat(
      latex_tab,
      file = file.path(
        tabs_p,
        paste0("tab_state_", cov, "_", lvl, ".tex")
      )
    )
  }
}

# Model fitting

# These two covariates had too few observations in some transition groups.
df_msm <- df_msm %>%
  select(-c(liver_disease, sppb_low))

# State 4 is absorbing.
Q <- matrix(1, nrow = 4, ncol = 4)
Q[4, ] <- 0

# Starting values
crude_int <- crudeinits.msm(state ~ access_age,
  subject = patient_id,
  data = df_msm,
  qmatrix = Q
)

# Null model
reframe_fit0 <- msm(state ~ access_age,
  subject = patient_id,
  data = df_msm,
  qmatrix = crude_int,
  gen.inits = FALSE,
  deathexact = TRUE
)
saveRDS(
  reframe_fit0,
  file = file.path(modl_p, "reframe_fit0.rds")
)

# Reframe 1
cov_list_1 <- list(
  "1-2" = ~access_age,
  "1-3" = ~access_age,
  "1-4" = ~ sex + access_age,
  "2-1" = ~access_age,
  "2-3" = ~access_age,
  "2-4" = ~ sex + access_age,
  "3-1" = ~access_age,
  "3-2" = ~access_age,
  "3-4" = ~ sex + access_age
)

constraints_1 = list(access_age = c(1, 2, 3, 4, 5, 3, 6, 7, 3))

## Starting values
qmatrix_init_1 <- qmatrix.msm(reframe_fit0)$estimates

## Model fit
reframe_fit1 <- msm(
  state ~ access_age,
  subject = patient_id,
  data = df_msm,
  covariates = cov_list_1,
  constraint = constraints_1,
  qmatrix = qmatrix_init_1,
  gen.inits = FALSE,
  deathexact = TRUE,
  opt.method = "optim",
  control = list(
    fnscale = reframe_fit0$minus2loglik,
    maxit   = 5000,
    trace   = 1,
    REPORT  = 10
  ),
  hessian = TRUE
)

saveRDS(
  reframe_fit1,
  file = file.path(modl_p, "reframe_fit1.rds")
)

# Reframe 2
cov_list_2 <- list(
  "1-2" = ~ access_age + hiv_tertile + phys_act_bin,
  "1-3" = ~ access_age + hiv_tertile + phys_act_bin,
  "1-4" = ~ sex + access_age,
  "2-1" = ~ access_age + hiv_tertile + phys_act_bin,
  "2-3" = ~ access_age + hiv_tertile + phys_act_bin,
  "2-4" = ~ sex + access_age,
  "3-1" = ~ access_age + hiv_tertile + phys_act_bin,
  "3-2" = ~ access_age + hiv_tertile + phys_act_bin,
  "3-4" = ~ sex + access_age
)

constraints_2 = list(access_age = c(1, 2, 3, 4, 5, 3, 6, 7, 3))

## Starting values
qmatrix_init_2 <- qmatrix.msm(reframe_fit1)$estimates

## Model fit
reframe_fit2 <- msm(
  state ~ access_age,
  subject = patient_id,
  data = df_msm,
  covariates = cov_list_2,
  constraint = constraints_2,
  qmatrix = qmatrix_init_2,
  gen.inits = FALSE,
  deathexact = TRUE,
  opt.method = "optim",
  control = list(
    fnscale = reframe_fit1$minus2loglik,
    maxit   = 5000,
    trace   = 1,
    REPORT  = 10
  ),
  hessian = TRUE
)

saveRDS(
  reframe_fit2,
  file = file.path(modl_p, "reframe_fit2.rds")
)

# Reframe 3
cov_list_3 <- list(
  "1-2" = ~ access_age + hiv_tertile + phys_act_bin + INSTI + TAF + IP + comorbidity_Dislipidemia,
  "1-3" = ~ access_age + hiv_tertile + phys_act_bin + INSTI + TAF + IP + comorbidity_Dislipidemia,
  "1-4" = ~ sex + access_age,
  "2-1" = ~ access_age + hiv_tertile + phys_act_bin + INSTI + TAF + IP + comorbidity_Dislipidemia,
  "2-3" = ~ access_age + hiv_tertile + phys_act_bin + INSTI + TAF + IP + comorbidity_Dislipidemia,
  "2-4" = ~ sex + access_age,
  "3-1" = ~ access_age + hiv_tertile + phys_act_bin + INSTI + TAF + IP + comorbidity_Dislipidemia,
  "3-2" = ~ access_age + hiv_tertile + phys_act_bin + INSTI + TAF + IP + comorbidity_Dislipidemia,
  "3-4" = ~ sex + access_age
)

constraints_3 = list(access_age = c(1, 2, 3, 4, 5, 3, 6, 7, 3))

## Starting values
qmatrix_init_3 <- qmatrix.msm(reframe_fit2)$estimates

## Model fit
reframe_fit3 <- msm(
  state ~ access_age,
  subject = patient_id,
  data = df_msm,
  covariates = cov_list_3,
  constraint = constraints_3,
  qmatrix = qmatrix_init_3,
  gen.inits = FALSE,
  deathexact = TRUE,
  opt.method = "optim",
  control = list(
    fnscale = reframe_fit2$minus2loglik,
    maxit   = 5000,
    trace   = 1,
    REPORT  = 10
  ),
  hessian = TRUE
)

saveRDS(
  reframe_fit3,
  file = file.path(modl_p, "reframe_fit3.rds")
)

# The published models end at fit 3. The archived script also explored quadratic age.
run_quadratic_age <- "--quadratic-age" %in% commandArgs(trailingOnly = TRUE)
if (run_quadratic_age) {
  df_msm$access_age2 <- I(df_msm$access_age)^2

  cov_list_4 <- list(
    "1-2" = ~ access_age + hiv_tertile + phys_act_bin + INSTI + TAF + IP + comorbidity_Dislipidemia,
    "1-3" = ~ access_age + hiv_tertile + phys_act_bin + INSTI + TAF + IP + comorbidity_Dislipidemia,
    "1-4" = ~ sex + access_age + access_age2,
    "2-1" = ~ access_age + hiv_tertile + phys_act_bin + INSTI + TAF + IP + comorbidity_Dislipidemia,
    "2-3" = ~ access_age + hiv_tertile + phys_act_bin + INSTI + TAF + IP + comorbidity_Dislipidemia,
    "2-4" = ~ sex + access_age + access_age2,
    "3-1" = ~ access_age + hiv_tertile + phys_act_bin + INSTI + TAF + IP + comorbidity_Dislipidemia,
    "3-2" = ~ access_age + hiv_tertile + phys_act_bin + INSTI + TAF + IP + comorbidity_Dislipidemia,
    "3-4" = ~ sex + access_age + access_age2
  )

  constraints_4 = list(
    access_age = c(1, 2, 3, 4, 5, 3, 6, 7, 3),
    access_age2 = c(1, 1, 1)
  )

  ## Starting values
  qmatrix_init_4 <- qmatrix.msm(reframe_fit3)$estimates

  ## Model fit
  reframe_fit4 <- msm(
    state ~ access_age,
    subject = patient_id,
    data = df_msm,
    covariates = cov_list_4,
    constraint = constraints_4,
    qmatrix = qmatrix_init_4,
    gen.inits = FALSE,
    deathexact = TRUE,
    opt.method = "optim",
    control = list(
      fnscale = reframe_fit3$minus2loglik,
      maxit   = 5000,
      trace   = 1,
      REPORT  = 10
    ),
    hessian = TRUE
  )

  saveRDS(
    reframe_fit4,
    file = file.path(modl_p, "reframe_fit4.rds")
  )
}

# Results for the report

# Model comparison (AIC)
AIC(reframe_fit0, reframe_fit1, reframe_fit2, reframe_fit3)
# AIC(reframe_fit0, reframe_fit1, reframe_fit2, reframe_fit3, reframe_fit4)

# Check predicted state prevalence.
plot.prevalence.msm(reframe_fit0)
plot.prevalence.msm(reframe_fit1)
plot.prevalence.msm(reframe_fit2)
plot.prevalence.msm(reframe_fit3)
# plot.prevalence.msm(reframe_fit4)

# Hazard ratios — Model 2
hzr_reframe_fit2 <- hazard.msm(reframe_fit2)

res_reframe_fit2 <- hazard_table_msm(reframe_fit2)

file_out <- file.path(tabs_p, "hazard_reframe_fit2.tex")

writeLines(
  as.character(res_reframe_fit2),
  file_out
)

# Hazard ratios — Model 3
hzr_reframe_fit3 <- hazard.msm(reframe_fit3)

res_reframe_fit3 <- hazard_table_msm(reframe_fit3)

file_out <- file.path(tabs_p, "hazard_reframe_fit3.tex")

writeLines(
  as.character(res_reframe_fit3),
  file_out
)

if (run_quadratic_age) {
  hzr_reframe_fit4 <- hazard.msm(reframe_fit4)

  res_reframe_fit4 <- hazard_table_msm(reframe_fit4)

  file_out <- file.path(tabs_p, "hazard_reframe_fit4.tex")

  writeLines(
    as.character(res_reframe_fit4),
    file_out
  )
}
