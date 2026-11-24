# Descriptive tables for continuous and categorical variables.

# Small helpers

# Ensure that a directory exists and return its path invisibly.
ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

# Convert a vector to character, trim spaces, and turn empty strings into NA.
as_clean_char <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    x <- stringr::str_trim(x)
    x <- dplyr::na_if(x, "")
  }
  as.character(x)
}

# Return a compact file-safe label.
sanitize_label <- function(x) {
  x <- as.character(x)
  stringr::str_replace_all(x, "[^A-Za-z0-9_\\-]+", "_")
}

# Summarize the amount of missingness in a vector.
na_stats <- function(x) {
  n_tot <- length(x)
  n_na <- sum(is.na(x))
  tibble::tibble(
    n_NA    = n_na,
    perc_NA = ifelse(n_tot == 0, NA_real_, 100 * n_na / n_tot)
  )
}

# Safely compute a quantile, returning NA if the input has no observed values.
q_safe <- function(x, p) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(stats::quantile(x, probs = p, names = FALSE, type = 7))
}

# Summarize a continuous vector using common descriptive statistics.
summ_cont <- function(x) {
  x0 <- x[!is.na(x)]
  if (length(x0) == 0) {
    return(tibble::tibble(
      mean = NA_real_, sd = NA_real_, median = NA_real_, mad = NA_real_,
      min = NA_real_, p5 = NA_real_, p10 = NA_real_, q1 = NA_real_,
      q3 = NA_real_, p90 = NA_real_, p95 = NA_real_, max = NA_real_
    ))
  }
  tibble::tibble(
    mean   = mean(x0),
    sd     = stats::sd(x0),
    median = stats::median(x0),
    mad    = stats::mad(x0),
    min    = min(x0),
    p5     = q_safe(x0, 0.05),
    p10    = q_safe(x0, 0.10),
    q1     = q_safe(x0, 0.25),
    q3     = q_safe(x0, 0.75),
    p90    = q_safe(x0, 0.90),
    p95    = q_safe(x0, 0.95),
    max    = max(x0)
  )
}

# Format a p-value into a compact human-readable string.
fmt_p <- function(p) {
  if (is.na(p)) return("p=NA")
  if (p < 1e-4) return("p<1e-4")
  sprintf("p=%.4f", p)
}

# Build a missingness table for a given set of variables.
build_missing <- function(df, variables = vars_miss) {
  n_tot <- nrow(df)

  tibble::tibble(variable = variables) %>%
    dplyr::mutate(
      n_missing = purrr::map_int(variable, ~ {
        if (!.x %in% names(df)) return(NA_integer_)
        sum(is.na(df[[.x]]))
      }),
      perc_missing = dplyr::if_else(
        is.na(n_missing) | n_tot == 0,
        NA_real_,
        100 * n_missing / n_tot
      ),
      missing = dplyr::if_else(
        is.na(n_missing),
        NA_character_,
        sprintf("%d (%.1f\\%%)", n_missing, perc_missing)
      )
    ) %>%
    dplyr::select(variable, missing, perc_missing)
}

# Merge two missingness tables and attach a variable group label.
merge_missing <- function(tab_full,
                          tab_2vis,
                          state_vars = vars_state,
                          covariate_vars = vars_covariates,
                          state_definition_vars = vars_state_def) {
  dplyr::left_join(
    tab_full,
    tab_2vis,
    by = "variable",
    suffix = c("_full", "_2vis")
  ) %>%
    dplyr::mutate(
      group = dplyr::case_when(
        variable %in% state_vars ~ "State variables",
        variable %in% covariate_vars ~ "Covariates",
        variable %in% state_definition_vars ~ "State definition variables",
        TRUE ~ "Other"
      )
    ) %>%
    dplyr::arrange(group, dplyr::desc(perc_missing_full))
}

# Grouping builders

# Add a visit index within each patient after ordering by the chosen date variable.
add_visit_index <- function(df, id = "patient_id", date_var = "access_date") {
  df %>%
    dplyr::arrange(.data[[id]], .data[[date_var]]) %>%
    dplyr::group_by(.data[[id]]) %>%
    dplyr::mutate(visit_index = dplyr::row_number()) %>%
    dplyr::ungroup()
}

# Add the calendar year extracted from a date variable.
add_calendar_year <- function(df, date_var = "access_date") {
  df %>%
    dplyr::mutate(year = lubridate::year(.data[[date_var]]))
}

# Add a historical period factor derived from a date variable.
add_period <- function(df, date_var = "access_date") {
  df %>%
    dplyr::mutate(
      period = dplyr::case_when(
        .data[[date_var]] >= as.Date("2005-01-01") & .data[[date_var]] < as.Date("2015-01-01") ~ "2005-2014",
        .data[[date_var]] >= as.Date("2015-01-01") & .data[[date_var]] < as.Date("2020-01-01") ~ "2015-2019",
        .data[[date_var]] >= as.Date("2020-01-01") ~ "2020+",
        TRUE ~ NA_character_
      ),
      period = factor(period, levels = c("2005-2014", "2015-2019", "2020+"))
    )
}

# Stratified tests

# Safely run a one-way ANOVA for a continuous variable across strata.
safe_anova <- function(df, xvar, strat_var = "state") {
  d <- df %>%
    dplyr::transmute(
      x = .data[[xvar]],
      g = as_clean_char(.data[[strat_var]])
    ) %>%
    dplyr::filter(!is.na(x), !is.na(g))

  k <- dplyr::n_distinct(d$g)
  n <- nrow(d)
  if (k < 2 || n < 3) return(NA_character_)

  tryCatch(
    {
      fit <- stats::aov(x ~ as.factor(g), data = d)
      sm <- summary(fit)[[1]]
      Fv <- as.numeric(sm$`F value`[1])
      df1 <- as.integer(sm$Df[1])
      df2 <- as.integer(sm$Df[2])
      pv <- as.numeric(sm$`Pr(>F)`[1])
      sprintf("F(%d,%d)=%.3f, %s, n=%d, k=%d", df1, df2, Fv, fmt_p(pv), n, k)
    },
    error = function(e) NA_character_
  )
}

# Safely run a Kruskal-Wallis test for a continuous variable across strata.
safe_kw <- function(df, xvar, strat_var = "state") {
  d <- df %>%
    dplyr::transmute(
      x = .data[[xvar]],
      g = as_clean_char(.data[[strat_var]])
    ) %>%
    dplyr::filter(!is.na(x), !is.na(g))

  k <- dplyr::n_distinct(d$g)
  n <- nrow(d)
  if (k < 2 || n < 3) return(NA_character_)

  tryCatch(
    {
      kt <- stats::kruskal.test(x ~ as.factor(g), data = d)
      H <- unname(kt$statistic)
      dfw <- unname(kt$parameter)
      pv <- unname(kt$p.value)
      sprintf("H(%d)=%.3f, %s, n=%d, k=%d", as.integer(dfw), as.numeric(H), fmt_p(pv), n, k)
    },
    error = function(e) NA_character_
  )
}

# Safely test association between a discrete variable and a stratification variable.
safe_prop_test <- function(df, xvar, strat_var = "state", simulate_B = 5000) {
  d <- df %>%
    dplyr::transmute(
      x = as_clean_char(.data[[xvar]]),
      g = as_clean_char(.data[[strat_var]])
    ) %>%
    dplyr::filter(!is.na(x), !is.na(g))

  if (nrow(d) == 0) return(NA_character_)
  if (dplyr::n_distinct(d$g) < 2) return(NA_character_)
  if (dplyr::n_distinct(d$x) < 2) return(NA_character_)

  tab <- table(d$g, d$x)

  tryCatch(
    {
      ct <- suppressWarnings(stats::chisq.test(tab))

      if (any(ct$expected < 5, na.rm = TRUE)) {
        ct2 <- suppressWarnings(stats::chisq.test(tab, simulate.p.value = TRUE, B = simulate_B))
        sprintf(
          "chisq_sim=%.3f, %s, n=%d",
          as.numeric(ct2$statistic), fmt_p(ct2$p.value), sum(tab)
        )
      } else {
        sprintf(
          "chisq(df=%d)=%.3f, %s, n=%d",
          as.integer(ct$parameter), as.numeric(ct$statistic), fmt_p(ct$p.value), sum(tab)
        )
      }
    },
    error = function(e) {
      tryCatch(
        {
          ft <- stats::fisher.test(tab)
          sprintf("fisher, %s, n=%d", fmt_p(ft$p.value), sum(tab))
        },
        error = function(e2) NA_character_
      )
    }
  )
}

# Summaries: non-stratified

# Build a continuous summary table for one grouped subset of data.
summ_cont_table <- function(df_g, vars_cont, id = "patient_id") {
  n_group <- nrow(df_g)
  sbj_group <- df_g %>%
    dplyr::distinct(.data[[id]]) %>%
    nrow()

  purrr::map_dfr(vars_cont, function(v) {
    x <- df_g[[v]]
    na <- na_stats(x)
    s <- summ_cont(x)

    tibble::tibble(
      variable = v,
      n        = n_group,
      sbj      = sbj_group,
      perc_NA  = na$perc_NA
    ) %>%
      dplyr::bind_cols(s)
  })
}

# Build a discrete summary table for one grouped subset of data.
summ_disc_table <- function(df_g, vars_disc, id = "patient_id") {
  n_group <- nrow(df_g)
  sbj_group <- df_g %>%
    dplyr::distinct(.data[[id]]) %>%
    nrow()

  purrr::map_dfr(vars_disc, function(v) {
    x <- as_clean_char(df_g[[v]])
    na <- na_stats(x)

    denom <- sum(!is.na(x))
    if (denom == 0) {
      return(tibble::tibble(
        variable = v,
        level    = NA_character_,
        n        = n_group,
        sbj      = sbj_group,
        perc_NA  = na$perc_NA,
        count    = NA_integer_,
        prop     = NA_real_
      ))
    }

    tab <- table(x, useNA = "no")

    tibble::tibble(
      variable = v,
      level    = names(tab),
      n        = n_group,
      sbj      = sbj_group,
      perc_NA  = na$perc_NA,
      count    = as.integer(tab),
      prop     = as.numeric(tab) / denom
    )
  })
}

# Summaries: stratified

# Build a wide continuous summary table stratified by a given variable.
summ_cont_table_by_strata <- function(df_g,
                                      vars_cont,
                                      id = "patient_id",
                                      strat_var = "state") {
  df2 <- df_g %>% dplyr::filter(!is.na(.data[[strat_var]]))

  purrr::map_dfr(vars_cont, function(v) {
    per_stratum <- df2 %>%
      dplyr::group_by(.data[[strat_var]]) %>%
      dplyr::group_modify(~ {
        x <- .x[[v]]
        na <- na_stats(x)
        s <- summ_cont(x)

        tibble::tibble(
          n       = nrow(.x),
          sbj     = .x %>% dplyr::distinct(.data[[id]]) %>% nrow(),
          perc_NA = na$perc_NA
        ) %>%
          dplyr::bind_cols(s)
      }) %>%
      dplyr::ungroup() %>%
      dplyr::rename(stratum = .data[[strat_var]]) %>%
      dplyr::mutate(stratum = sanitize_label(stratum))

    wide <- per_stratum %>%
      tidyr::pivot_wider(
        names_from  = stratum,
        values_from = c(n, sbj, perc_NA, mean, sd, median, mad, min, p5, p10, q1, q3, p90, p95, max),
        names_glue  = "{.value}_lv{stratum}"
      )

    tibble::tibble(variable = v) %>%
      dplyr::bind_cols(wide) %>%
      dplyr::mutate(
        anova = safe_anova(df2, v, strat_var = strat_var),
        kw    = safe_kw(df2, v, strat_var = strat_var)
      )
  })
}

# Build a wide discrete summary table stratified by a given variable.
summ_disc_table_by_strata <- function(df_g,
                                      vars_disc,
                                      id = "patient_id",
                                      strat_var = "state") {
  df2 <- df_g %>% dplyr::filter(!is.na(.data[[strat_var]]))
  strata_levels <- df2 %>%
    dplyr::distinct(.data[[strat_var]]) %>%
    dplyr::arrange(.data[[strat_var]]) %>%
    dplyr::pull() %>%
    as.character()

  purrr::map_dfr(vars_disc, function(v) {
    x_all <- as_clean_char(df2[[v]])
    levels_all <- sort(unique(x_all[!is.na(x_all)]))

    if (length(levels_all) == 0) {
      return(tibble::tibble(
        variable  = v,
        level     = NA_character_,
        prop_test = safe_prop_test(df2, v, strat_var = strat_var)
      ))
    }

    base <- tibble::tibble(variable = v, level = levels_all)

    out <- purrr::reduce(strata_levels, function(acc, lv) {
      d_lv <- df2 %>% dplyr::filter(as.character(.data[[strat_var]]) == lv)
      x <- as_clean_char(d_lv[[v]])

      n_lv <- nrow(d_lv)
      sbj_lv <- d_lv %>%
        dplyr::distinct(.data[[id]]) %>%
        nrow()
      na_lv <- na_stats(x)
      denom <- sum(!is.na(x))
      lv_lab <- sanitize_label(lv)

      tab <- if (denom == 0) {
        stats::setNames(rep(0L, length(levels_all)), levels_all)
      } else {
        tmp <- table(x, useNA = "no")
        full <- stats::setNames(rep(0L, length(levels_all)), levels_all)
        full[names(tmp)] <- as.integer(tmp)
        full
      }

      props <- if (denom == 0) {
        rep(NA_real_, length(levels_all))
      } else {
        as.integer(tab) / denom
      }

      df_join <- tibble::tibble(
        variable = v,
        level = levels_all,
        !!paste0("n_lv", lv_lab) := n_lv,
        !!paste0("sbj_lv", lv_lab) := sbj_lv,
        !!paste0("perc_NA_lv", lv_lab) := na_lv$perc_NA,
        !!paste0("count_lv", lv_lab) := as.integer(tab),
        !!paste0("prop_lv", lv_lab) := props
      )

      dplyr::left_join(acc, df_join, by = c("variable", "level"))
    }, .init = base)

    out %>%
      dplyr::mutate(prop_test = safe_prop_test(df2, v, strat_var = strat_var)) %>%
      dplyr::relocate(prop_test, .after = level)
  })
}

# Excel writer

# Write continuous and discrete summary tables into a two-sheet Excel workbook.
write_group_workbook <- function(cont_tbl, disc_tbl, file_path) {
  flatten_df <- function(df) {
    df %>%
      dplyr::mutate(dplyr::across(where(is.list), ~ purrr::map_chr(.x, paste, collapse = ";"))) %>%
      as.data.frame()
  }

  cont_tbl <- flatten_df(cont_tbl)
  disc_tbl <- flatten_df(disc_tbl)

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "continuous")
  openxlsx::addWorksheet(wb, "discrete")
  openxlsx::writeData(wb, "continuous", cont_tbl)
  openxlsx::writeData(wb, "discrete", disc_tbl)

  ensure_dir(dirname(file_path))
  openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
}

# Main runner

# Build and export descriptive workbooks for each level of a grouping variable.
run_descriptive_pipeline <- function(df,
                                     group_var,
                                     vars_cont,
                                     vars_disc,
                                     output_dir,
                                     file_prefix,
                                     stratify_by = FALSE,
                                     id = "patient_id",
                                     strat_var = "state") {
  ensure_dir(output_dir)

  df0 <- df %>% dplyr::filter(!is.na(.data[[group_var]]))

  groups <- df0 %>%
    dplyr::distinct(.data[[group_var]]) %>%
    dplyr::arrange(.data[[group_var]]) %>%
    dplyr::pull()

  purrr::walk(groups, function(gv) {
    df_g <- df0 %>% dplyr::filter(.data[[group_var]] == gv)
    g_lab <- sanitize_label(gv)

    if (!stratify_by) {
      cont_tbl <- summ_cont_table(df_g, vars_cont, id = id)
      disc_tbl <- summ_disc_table(df_g, vars_disc, id = id)
      fpath <- file.path(output_dir, sprintf("%s_%s.xlsx", file_prefix, g_lab))
      write_group_workbook(cont_tbl, disc_tbl, fpath)
    } else {
      cont_tbl <- summ_cont_table_by_strata(df_g, vars_cont, id = id, strat_var = strat_var)
      disc_tbl <- summ_disc_table_by_strata(df_g, vars_disc, id = id, strat_var = strat_var)
      fpath <- file.path(
        output_dir,
        sprintf("%s_%s_by_%s.xlsx", file_prefix, g_lab, sanitize_label(strat_var))
      )
      write_group_workbook(cont_tbl, disc_tbl, fpath)
    }
  })

  invisible(TRUE)
}
