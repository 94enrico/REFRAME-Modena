# Agreement between BMI and REFRAME state, with participant-level bootstrap intervals.

# Constants

AGREEMENT_PERIOD_LEVELS <- c("2005-2014", "2015-2019", "2020+")

AGREEMENT_STAT_NAMES <- local({
  out <- c()

  for (j in 1:3) {
    for (k in 1:3) {
      out <- c(out, paste0("P_O", j, "_B", k))
    }
  }

  for (j in 1:3) {
    for (k in 1:3) {
      out <- c(out, paste0("P_O", j, "_given_B", k))
    }
  }

  for (k in 1:3) {
    for (j in 1:3) {
      out <- c(out, paste0("P_B", k, "_given_O", j))
    }
  }

  out <- c(out, "p_agree", "kappa_w", "p_O_gt_B", "p_O_lt_B")
  out
})

AGREEMENT_STATE_LABELS <- c(
  "1" = "metabolically healthy",
  "2" = "preclinical obesity",
  "3" = "clinical obesity"
)

AGREEMENT_BMI_LABELS <- c(
  "1" = "normal BMI",
  "2" = "overweight BMI",
  "3" = "obesity BMI"
)

AGREEMENT_STAT_LABELS <- local({
  labels <- setNames(rep(NA_character_, length(AGREEMENT_STAT_NAMES)), AGREEMENT_STAT_NAMES)

  for (j in 1:3) {
    for (k in 1:3) {
      stat_name <- paste0("P_O", j, "_B", k)
      labels[stat_name] <- paste0(
        "P(obesity = ", AGREEMENT_STATE_LABELS[as.character(j)],
        ", BMI = ", AGREEMENT_BMI_LABELS[as.character(k)], ")"
      )
    }
  }

  for (j in 1:3) {
    for (k in 1:3) {
      stat_name <- paste0("P_O", j, "_given_B", k)
      labels[stat_name] <- paste0(
        "P(obesity = ", AGREEMENT_STATE_LABELS[as.character(j)],
        " | BMI = ", AGREEMENT_BMI_LABELS[as.character(k)], ")"
      )
    }
  }

  for (k in 1:3) {
    for (j in 1:3) {
      stat_name <- paste0("P_B", k, "_given_O", j)
      labels[stat_name] <- paste0(
        "P(BMI = ", AGREEMENT_BMI_LABELS[as.character(k)],
        " | obesity = ", AGREEMENT_STATE_LABELS[as.character(j)], ")"
      )
    }
  }

  labels["p_agree"] <- "Overall agreement"
  labels["kappa_w"] <- "Weighted kappa (quadratic weights)"
  labels["p_O_gt_B"] <- "P(obesity state > BMI category)"
  labels["p_O_lt_B"] <- "P(obesity state < BMI category)"

  labels
})

# Put P(clinical obesity | BMI) first
AGREEMENT_STAT_ORDER <- c(
  "P_O1_given_B1",
  "P_O2_given_B1",
  "P_O3_given_B1",
  "P_O1_given_B2",
  "P_O2_given_B2",
  "P_O3_given_B2",
  "P_O1_given_B3",
  "P_O2_given_B3",
  "P_O3_given_B3",
  "p_agree",
  "kappa_w",
  "p_O_gt_B",
  "p_O_lt_B",
  setdiff(
    AGREEMENT_STAT_NAMES,
    c(
      "P_O1_given_B1",
      "P_O2_given_B1",
      "P_O3_given_B1",
      "P_O1_given_B2",
      "P_O2_given_B2",
      "P_O3_given_B2",
      "P_O1_given_B3",
      "P_O2_given_B3",
      "P_O3_given_B3",
      "p_agree",
      "kappa_w",
      "p_O_gt_B",
      "p_O_lt_B"
    )
  )
)

AGREEMENT_SUBGROUP_LEVELS <- c("overall", "male", "female")

# Utility functions

make_period <- function(access_date) {
  year_vec <- lubridate::year(access_date)

  period <- dplyr::case_when(
    !is.na(year_vec) & year_vec >= 2005 & year_vec <= 2014 ~ "2005-2014",
    !is.na(year_vec) & year_vec >= 2015 & year_vec <= 2019 ~ "2015-2019",
    !is.na(year_vec) & year_vec >= 2020 ~ "2020+",
    TRUE ~ NA_character_
  )

  factor(period, levels = AGREEMENT_PERIOD_LEVELS)
}

safe_weighted_kappa_from_tab <- function(tab) {
  n <- sum(tab)

  if (n <= 1) {
    return(NA_real_)
  }

  O <- tab / n
  r <- rowSums(O)
  c <- colSums(O)
  E <- outer(r, c)

  K <- nrow(tab)
  w <- outer(
    1:K, 1:K,
    FUN = function(i, j) ((i - j) / (K - 1))^2
  )

  den <- sum(w * E)

  if (!is.finite(den) || den <= 0) {
    return(NA_real_)
  }

  1 - sum(w * O) / den
}

agreement_stats <- function(data) {
  if (is.null(data) || nrow(data) == 0) {
    out <- rep(NA_real_, length(AGREEMENT_STAT_NAMES))
    names(out) <- AGREEMENT_STAT_NAMES
    return(out)
  }

  tab <- table(
    factor(data$state, levels = 1:3),
    factor(data$bmi_cat, levels = 1:3)
  )

  tab <- as.matrix(tab)
  n <- sum(tab)

  if (n == 0) {
    out <- rep(NA_real_, length(AGREEMENT_STAT_NAMES))
    names(out) <- AGREEMENT_STAT_NAMES
    return(out)
  }

  p_joint <- tab / n

  col_tot <- colSums(tab)
  p_OgB <- sweep(tab, 2, col_tot, "/")
  p_OgB[, col_tot == 0] <- NA_real_

  row_tot <- rowSums(tab)
  p_BgO <- sweep(tab, 1, row_tot, "/")
  p_BgO[row_tot == 0, ] <- NA_real_

  p_agree <- sum(diag(tab)) / n

  idx <- expand.grid(O = 1:3, B = 1:3)
  probs_vec <- as.vector(p_joint)

  p_O_gt_B <- sum(probs_vec[idx$O > idx$B])
  p_O_lt_B <- sum(probs_vec[idx$O < idx$B])

  kappa_w <- safe_weighted_kappa_from_tab(tab)

  out <- c()

  for (j in 1:3) {
    for (k in 1:3) {
      out[paste0("P_O", j, "_B", k)] <- p_joint[j, k]
    }
  }

  for (j in 1:3) {
    for (k in 1:3) {
      out[paste0("P_O", j, "_given_B", k)] <- p_OgB[j, k]
    }
  }

  for (k in 1:3) {
    for (j in 1:3) {
      out[paste0("P_B", k, "_given_O", j)] <- p_BgO[j, k]
    }
  }

  out["p_agree"] <- p_agree
  out["kappa_w"] <- kappa_w
  out["p_O_gt_B"] <- p_O_gt_B
  out["p_O_lt_B"] <- p_O_lt_B

  out[AGREEMENT_STAT_NAMES]
}

safe_boot_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE))
}

format_estimate_ci <- function(estimate, ci_low, ci_high, digits = 3) {
  out <- rep(NA_character_, length(estimate))
  ok <- is.finite(estimate) & is.finite(ci_low) & is.finite(ci_high)

  out[ok] <- paste0(
    sprintf(paste0("%.", digits, "f"), estimate[ok]),
    " (",
    sprintf(paste0("%.", digits, "f"), ci_low[ok]),
    ", ",
    sprintf(paste0("%.", digits, "f"), ci_high[ok]),
    ")"
  )

  out
}

order_agreement_results <- function(data, diff_table = FALSE) {
  data <- data %>%
    dplyr::mutate(
      statistic = factor(statistic, levels = AGREEMENT_STAT_ORDER),
      label = factor(label, levels = unname(AGREEMENT_STAT_LABELS[AGREEMENT_STAT_ORDER]))
    ) %>%
    dplyr::arrange(statistic)

  if (diff_table) {
    data <- data %>%
      dplyr::mutate(
        comparison = paste0(period_2, " - ", period_1)
      )
  }

  data
}

make_summary_wide <- function(data, by_period = FALSE) {
  data <- data %>%
    dplyr::mutate(
      subgroup = factor(subgroup, levels = AGREEMENT_SUBGROUP_LEVELS)
    )

  if (by_period) {
    out <- data %>%
      dplyr::mutate(
        period = factor(period, levels = AGREEMENT_PERIOD_LEVELS),
        value = format_estimate_ci(estimate, ci_low, ci_high)
      ) %>%
      dplyr::select(label, statistic, subgroup, period, value) %>%
      tidyr::pivot_wider(
        names_from = c(subgroup, period),
        values_from = value,
        names_glue = "{subgroup} | {period}"
      )
  } else {
    out <- data %>%
      dplyr::mutate(
        value = format_estimate_ci(estimate, ci_low, ci_high)
      ) %>%
      dplyr::select(label, statistic, subgroup, value) %>%
      tidyr::pivot_wider(
        names_from = subgroup,
        values_from = value
      )
  }

  out %>%
    dplyr::mutate(statistic = factor(statistic, levels = AGREEMENT_STAT_ORDER)) %>%
    dplyr::arrange(statistic) %>%
    dplyr::mutate(label = as.character(label))
}

make_pairwise_wide <- function(data) {
  data %>%
    dplyr::mutate(
      subgroup = factor(subgroup, levels = AGREEMENT_SUBGROUP_LEVELS),
      comparison = factor(
        paste0(period_2, " - ", period_1),
        levels = c(
          "2015-2019 - 2005-2014",
          "2020+ - 2005-2014",
          "2020+ - 2015-2019"
        )
      ),
      value = format_estimate_ci(diff_estimate, ci_low, ci_high)
    ) %>%
    dplyr::select(label, statistic, subgroup, comparison, value) %>%
    tidyr::pivot_wider(
      names_from = c(subgroup, comparison),
      values_from = value,
      names_glue = "{subgroup} | {comparison}"
    ) %>%
    dplyr::mutate(statistic = factor(statistic, levels = AGREEMENT_STAT_ORDER)) %>%
    dplyr::arrange(statistic) %>%
    dplyr::mutate(label = as.character(label))
}

make_kappa_pairwise_wide <- function(data) {
  data %>%
    dplyr::filter(statistic == "kappa_w") %>%
    dplyr::mutate(
      subgroup = factor(subgroup, levels = AGREEMENT_SUBGROUP_LEVELS),
      comparison = factor(
        paste0(period_2, " - ", period_1),
        levels = c(
          "2015-2019 - 2005-2014",
          "2020+ - 2005-2014",
          "2020+ - 2015-2019"
        )
      ),
      value = format_estimate_ci(diff_estimate, ci_low, ci_high)
    ) %>%
    dplyr::select(subgroup, comparison, value) %>%
    tidyr::pivot_wider(
      names_from = comparison,
      values_from = value
    ) %>%
    dplyr::arrange(subgroup)
}

write_latex_table <- function(data, file_path, caption, align = NULL) {
  if (is.null(align)) {
    align <- paste0("l", strrep("c", ncol(data) - 1))
  }

  tab_latex <- kableExtra::kbl(
    data,
    format = "latex",
    booktabs = TRUE,
    align = align,
    caption = caption
  ) %>%
    kableExtra::kable_styling(latex_options = c("hold_position", "scale_down"))

  writeLines(as.character(tab_latex), con = file_path)
}

write_summary_workbook <- function(long_df,
                                   wide_df,
                                   xlsx_path,
                                   csv_path_long,
                                   csv_path_wide,
                                   diff_long_df = NULL,
                                   diff_wide_df = NULL,
                                   csv_path_diff_long = NULL,
                                   csv_path_diff_wide = NULL) {
  wb <- openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, "estimates_long")
  openxlsx::writeData(wb, "estimates_long", long_df)

  openxlsx::addWorksheet(wb, "estimates_wide")
  openxlsx::writeData(wb, "estimates_wide", wide_df)

  if (!is.null(diff_long_df)) {
    openxlsx::addWorksheet(wb, "pairwise_diff_long")
    openxlsx::writeData(wb, "pairwise_diff_long", diff_long_df)
  }

  if (!is.null(diff_wide_df)) {
    openxlsx::addWorksheet(wb, "pairwise_diff_wide")
    openxlsx::writeData(wb, "pairwise_diff_wide", diff_wide_df)
  }

  for (sheet in wb$sheet_names) {
    openxlsx::setColWidths(wb, sheet = sheet, cols = 1:100, widths = "auto")
  }

  openxlsx::saveWorkbook(wb, file = xlsx_path, overwrite = TRUE)

  utils::write.csv(long_df, file = csv_path_long, row.names = FALSE)
  utils::write.csv(wide_df, file = csv_path_wide, row.names = FALSE)

  if (!is.null(diff_long_df) && !is.null(csv_path_diff_long)) {
    utils::write.csv(diff_long_df, file = csv_path_diff_long, row.names = FALSE)
  }

  if (!is.null(diff_wide_df) && !is.null(csv_path_diff_wide)) {
    utils::write.csv(diff_wide_df, file = csv_path_diff_wide, row.names = FALSE)
  }
}

write_confusion_matrix_tex <- function(data, file_path, caption) {
  tab <- with(data, table(Obesity = state, BMI = bmi_cat))

  tab_df <- as.data.frame.matrix(tab)
  tab_df <- tibble::rownames_to_column(tab_df, var = "Obesity")

  write_latex_table(
    data = tab_df,
    file_path = file_path,
    caption = caption,
    align = paste0("l", strrep("c", ncol(tab_df) - 1))
  )
}

# Fast bootstrap machinery

make_cluster_rows <- function(data) {
  split(seq_len(nrow(data)), data$patient_id, drop = TRUE)
}

sample_cluster_rows_all_visits <- function(cluster_rows) {
  if (length(cluster_rows) == 0) {
    return(integer(0))
  }
  sampled_clusters <- sample.int(length(cluster_rows), size = length(cluster_rows), replace = TRUE)
  unlist(cluster_rows[sampled_clusters], use.names = FALSE)
}

sample_one_row_per_cluster <- function(cluster_rows) {
  if (length(cluster_rows) == 0) {
    return(integer(0))
  }
  vapply(cluster_rows, function(idx) sample(idx, size = 1), integer(1))
}

sample_cluster_rows_random_visit <- function(cluster_rows) {
  if (length(cluster_rows) == 0) {
    return(integer(0))
  }
  sampled_clusters <- sample.int(length(cluster_rows), size = length(cluster_rows), replace = TRUE)
  vapply(sampled_clusters, function(i) sample(cluster_rows[[i]], size = 1), integer(1))
}

summarise_bootstrap_estimates <- function(theta_hat, boot_mat = NULL) {
  if (is.null(boot_mat)) {
    ci_low <- rep(NA_real_, length(theta_hat))
    ci_high <- rep(NA_real_, length(theta_hat))
  } else {
    ci_low <- apply(boot_mat, 2, safe_boot_quantile, prob = 0.025)
    ci_high <- apply(boot_mat, 2, safe_boot_quantile, prob = 0.975)
  }

  tibble::tibble(
    statistic = AGREEMENT_STAT_NAMES,
    label     = unname(AGREEMENT_STAT_LABELS[AGREEMENT_STAT_NAMES]),
    estimate  = as.numeric(theta_hat[AGREEMENT_STAT_NAMES]),
    ci_low    = as.numeric(ci_low[AGREEMENT_STAT_NAMES]),
    ci_high   = as.numeric(ci_high[AGREEMENT_STAT_NAMES])
  )
}

bootstrap_from_sampler <- function(data,
                                   point_data,
                                   sampler_fun,
                                   B = 2000,
                                   seed = 12345) {
  theta_hat <- agreement_stats(point_data)

  if (nrow(point_data) == 0 || B <= 0) {
    return(list(
      estimates  = summarise_bootstrap_estimates(theta_hat, boot_mat = NULL),
      boot_mat   = NULL,
      point_data = point_data
    ))
  }

  set.seed(seed)

  boot_mat <- matrix(
    NA_real_,
    nrow = B,
    ncol = length(AGREEMENT_STAT_NAMES),
    dimnames = list(NULL, AGREEMENT_STAT_NAMES)
  )

  for (b in seq_len(B)) {
    row_idx <- sampler_fun()
    boot_mat[b, ] <- agreement_stats(data[row_idx, , drop = FALSE])[AGREEMENT_STAT_NAMES]
  }

  list(
    estimates  = summarise_bootstrap_estimates(theta_hat, boot_mat = boot_mat),
    boot_mat   = boot_mat,
    point_data = point_data
  )
}

bootstrap_all_visits <- function(data,
                                 B = 2000,
                                 seed = 12345) {
  point_data <- data

  if (nrow(data) == 0) {
    return(bootstrap_from_sampler(
      data = data,
      point_data = point_data,
      sampler_fun = function() integer(0),
      B = B,
      seed = seed
    ))
  }

  cluster_rows <- make_cluster_rows(data)

  sampler_fun <- function() {
    sample_cluster_rows_all_visits(cluster_rows)
  }

  bootstrap_from_sampler(
    data = data,
    point_data = point_data,
    sampler_fun = sampler_fun,
    B = B,
    seed = seed
  )
}

bootstrap_random_visit <- function(data,
                                   B = 2000,
                                   seed = 12345,
                                   point_seed = 54321) {
  if (nrow(data) == 0) {
    return(bootstrap_from_sampler(
      data = data,
      point_data = data,
      sampler_fun = function() integer(0),
      B = B,
      seed = seed
    ))
  }

  cluster_rows <- make_cluster_rows(data)

  set.seed(point_seed)
  point_idx <- sample_one_row_per_cluster(cluster_rows)
  point_data <- data[point_idx, , drop = FALSE]

  sampler_fun <- function() {
    sample_cluster_rows_random_visit(cluster_rows)
  }

  bootstrap_from_sampler(
    data = data,
    point_data = point_data,
    sampler_fun = sampler_fun,
    B = B,
    seed = seed
  )
}

bootstrap_first_visit <- function(data,
                                  B = 2000,
                                  seed = 12345) {
  point_data <- data

  if (nrow(data) == 0) {
    return(bootstrap_from_sampler(
      data = data,
      point_data = point_data,
      sampler_fun = function() integer(0),
      B = B,
      seed = seed
    ))
  }

  n <- nrow(data)

  sampler_fun <- function() {
    sample.int(n, size = n, replace = TRUE)
  }

  bootstrap_from_sampler(
    data = data,
    point_data = point_data,
    sampler_fun = sampler_fun,
    B = B,
    seed = seed
  )
}

bootstrap_period_last <- function(data,
                                  B = 2000,
                                  seed = 12345) {
  period_pairs <- utils::combn(AGREEMENT_PERIOD_LEVELS, 2, simplify = FALSE)

  point_stats <- lapply(
    AGREEMENT_PERIOD_LEVELS,
    function(p) agreement_stats(data[data$period == p, , drop = FALSE])
  )
  names(point_stats) <- AGREEMENT_PERIOD_LEVELS

  point_diffs <- lapply(
    period_pairs,
    function(pair) point_stats[[pair[2]]] - point_stats[[pair[1]]]
  )

  if (nrow(data) == 0 || B <= 0) {
    estimates_long <- dplyr::bind_rows(lapply(AGREEMENT_PERIOD_LEVELS, function(p) {
      tibble::tibble(
        statistic = AGREEMENT_STAT_NAMES,
        label     = unname(AGREEMENT_STAT_LABELS[AGREEMENT_STAT_NAMES]),
        period    = p,
        estimate  = as.numeric(point_stats[[p]][AGREEMENT_STAT_NAMES]),
        ci_low    = NA_real_,
        ci_high   = NA_real_
      )
    }))

    pairwise_long <- dplyr::bind_rows(lapply(seq_along(period_pairs), function(i) {
      pair <- period_pairs[[i]]
      diff_hat <- point_diffs[[i]]

      tibble::tibble(
        statistic     = AGREEMENT_STAT_NAMES,
        label         = unname(AGREEMENT_STAT_LABELS[AGREEMENT_STAT_NAMES]),
        period_1      = pair[1],
        period_2      = pair[2],
        diff_estimate = as.numeric(diff_hat[AGREEMENT_STAT_NAMES]),
        ci_low        = NA_real_,
        ci_high       = NA_real_
      )
    }))

    return(list(
      estimates_long = estimates_long,
      pairwise_long  = pairwise_long
    ))
  }

  cluster_rows <- make_cluster_rows(data)

  boot_estimates <- lapply(AGREEMENT_PERIOD_LEVELS, function(p) {
    matrix(
      NA_real_,
      nrow = B,
      ncol = length(AGREEMENT_STAT_NAMES),
      dimnames = list(NULL, AGREEMENT_STAT_NAMES)
    )
  })
  names(boot_estimates) <- AGREEMENT_PERIOD_LEVELS

  boot_diffs <- lapply(period_pairs, function(pair) {
    matrix(
      NA_real_,
      nrow = B,
      ncol = length(AGREEMENT_STAT_NAMES),
      dimnames = list(NULL, AGREEMENT_STAT_NAMES)
    )
  })

  set.seed(seed)

  for (b in seq_len(B)) {
    row_idx <- sample_cluster_rows_all_visits(cluster_rows)
    boot_data <- data[row_idx, , drop = FALSE]

    boot_stats <- lapply(
      AGREEMENT_PERIOD_LEVELS,
      function(p) agreement_stats(boot_data[boot_data$period == p, , drop = FALSE])
    )
    names(boot_stats) <- AGREEMENT_PERIOD_LEVELS

    for (p in AGREEMENT_PERIOD_LEVELS) {
      boot_estimates[[p]][b, ] <- boot_stats[[p]][AGREEMENT_STAT_NAMES]
    }

    for (i in seq_along(period_pairs)) {
      pair <- period_pairs[[i]]
      boot_diffs[[i]][b, ] <- (boot_stats[[pair[2]]] - boot_stats[[pair[1]]])[AGREEMENT_STAT_NAMES]
    }
  }

  estimates_long <- dplyr::bind_rows(lapply(AGREEMENT_PERIOD_LEVELS, function(p) {
    ci_low <- apply(boot_estimates[[p]], 2, safe_boot_quantile, prob = 0.025)
    ci_high <- apply(boot_estimates[[p]], 2, safe_boot_quantile, prob = 0.975)

    tibble::tibble(
      statistic = AGREEMENT_STAT_NAMES,
      label     = unname(AGREEMENT_STAT_LABELS[AGREEMENT_STAT_NAMES]),
      period    = p,
      estimate  = as.numeric(point_stats[[p]][AGREEMENT_STAT_NAMES]),
      ci_low    = as.numeric(ci_low[AGREEMENT_STAT_NAMES]),
      ci_high   = as.numeric(ci_high[AGREEMENT_STAT_NAMES])
    )
  }))

  pairwise_long <- dplyr::bind_rows(lapply(seq_along(period_pairs), function(i) {
    pair <- period_pairs[[i]]
    diff_hat <- point_diffs[[i]]

    ci_low <- apply(boot_diffs[[i]], 2, safe_boot_quantile, prob = 0.025)
    ci_high <- apply(boot_diffs[[i]], 2, safe_boot_quantile, prob = 0.975)

    tibble::tibble(
      statistic     = AGREEMENT_STAT_NAMES,
      label         = unname(AGREEMENT_STAT_LABELS[AGREEMENT_STAT_NAMES]),
      period_1      = pair[1],
      period_2      = pair[2],
      diff_estimate = as.numeric(diff_hat[AGREEMENT_STAT_NAMES]),
      ci_low        = as.numeric(ci_low[AGREEMENT_STAT_NAMES]),
      ci_high       = as.numeric(ci_high[AGREEMENT_STAT_NAMES])
    )
  }))

  list(
    estimates_long = estimates_long,
    pairwise_long  = pairwise_long
  )
}

# Wrappers for the requested analyses

run_standard_setting <- function(data,
                                 setting_name,
                                 analysis_type = c("all_visits", "random_visit", "first_visit"),
                                 B = 2000,
                                 seed_base = 1000) {
  analysis_type <- match.arg(analysis_type)

  subgroup_defs <- list(
    overall = function(x) x,
    male    = function(x) dplyr::filter(x, sex == 1),
    female  = function(x) dplyr::filter(x, sex == 0)
  )

  estimate_list <- list()
  counter <- 0

  for (subgroup_name in names(subgroup_defs)) {
    counter <- counter + 1
    subgroup_data <- subgroup_defs[[subgroup_name]](data)

    if (analysis_type == "all_visits") {
      boot_out <- bootstrap_all_visits(
        data = subgroup_data,
        B = B,
        seed = seed_base + counter
      )
    } else if (analysis_type == "random_visit") {
      boot_out <- bootstrap_random_visit(
        data = subgroup_data,
        B = B,
        seed = seed_base + counter,
        point_seed = seed_base + 50000 + counter
      )
    } else if (analysis_type == "first_visit") {
      boot_out <- bootstrap_first_visit(
        data = subgroup_data,
        B = B,
        seed = seed_base + counter
      )
    }

    estimate_list[[length(estimate_list) + 1]] <- boot_out$estimates %>%
      dplyr::mutate(
        setting = setting_name,
        subgroup = subgroup_name,
        n_obs = nrow(boot_out$point_data),
        n_patients = dplyr::n_distinct(boot_out$point_data$patient_id)
      )
  }

  estimates_long <- dplyr::bind_rows(estimate_list) %>%
    order_agreement_results(diff_table = FALSE)

  estimates_wide <- make_summary_wide(estimates_long, by_period = FALSE)

  list(
    estimates_long = estimates_long,
    estimates_wide = estimates_wide,
    pairwise_long = NULL,
    pairwise_wide = NULL
  )
}

run_period_last_setting <- function(data,
                                    setting_name,
                                    B = 2000,
                                    seed_base = 1000) {
  subgroup_defs <- list(
    overall = function(x) x,
    male    = function(x) dplyr::filter(x, sex == 1),
    female  = function(x) dplyr::filter(x, sex == 0)
  )

  estimate_list <- list()
  pairwise_list <- list()
  counter <- 0

  for (subgroup_name in names(subgroup_defs)) {
    counter <- counter + 1
    subgroup_data <- subgroup_defs[[subgroup_name]](data)

    boot_out <- bootstrap_period_last(
      data = subgroup_data,
      B = B,
      seed = seed_base + counter
    )

    period_counts <- subgroup_data %>%
      dplyr::group_by(period) %>%
      dplyr::summarise(
        n_obs = dplyr::n(),
        n_patients = dplyr::n_distinct(patient_id),
        .groups = "drop"
      )

    estimate_list[[length(estimate_list) + 1]] <- boot_out$estimates_long %>%
      dplyr::left_join(period_counts, by = "period") %>%
      dplyr::mutate(
        setting = setting_name,
        subgroup = subgroup_name
      )

    pairwise_list[[length(pairwise_list) + 1]] <- boot_out$pairwise_long %>%
      dplyr::mutate(
        setting = setting_name,
        subgroup = subgroup_name
      )
  }

  estimates_long <- dplyr::bind_rows(estimate_list) %>%
    order_agreement_results(diff_table = FALSE)

  estimates_wide <- make_summary_wide(estimates_long, by_period = TRUE)

  pairwise_long <- dplyr::bind_rows(pairwise_list) %>%
    order_agreement_results(diff_table = TRUE) %>%
    dplyr::mutate(
      comparison = paste0(period_2, " - ", period_1)
    )

  pairwise_wide <- make_pairwise_wide(pairwise_long)

  list(
    estimates_long = estimates_long,
    estimates_wide = estimates_wide,
    pairwise_long = pairwise_long,
    pairwise_wide = pairwise_wide
  )
}
