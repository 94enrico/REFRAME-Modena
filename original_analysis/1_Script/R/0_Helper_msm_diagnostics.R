# Explore baseline transition hazards in panel data.

# Utilities

# Extract state space

get_states <- function(df_intervals) {
  sort(unique(c(df_intervals$s0, df_intervals$s1)))
}

# Parse band limits from cut labels

band_lr <- function(band) {
  x <- gsub("\\[|\\(|\\]|\\)|\\s", "", band)
  ab <- strsplit(x, ",")[[1]]
  c(left = as.numeric(ab[1]), right = as.numeric(ab[2]))
}

# Compute midpoint of a band

band_mid <- function(band) {
  x <- gsub("\\[|\\(|\\]|\\)|\\s", "", band)
  ab <- strsplit(x, ",")[[1]]
  mean(as.numeric(ab))
}

# Predict Gompertz hazard on a grid

predict_gompertz <- function(fit, grid) {
  if (is.null(fit))
    return(NULL)

  pr <- tryCatch(
    predict(
      fit,
      newdata = tibble(t_mid = grid, dt = 1),
      se.fit = TRUE,
      type = "link"
    ),
    error = function(e) NULL
  )

  if (is.null(pr))
    return(NULL)

  tibble(
    age = grid,
    hazard = exp(pr$fit),
    lo = exp(pr$fit - 1.96 * pr$se.fit),
    hi = exp(pr$fit + 1.96 * pr$se.fit)
  )
}

# Data preparation

# Build interval dataset from panel observations

build_intervals <- function(df,
                            id = patient_id,
                            time = age,
                            state = CKM) {
  # Converts panel data into counting-process style intervals.
  # Each row represents exposure in state s0 over (t0, t1].

  df %>%
    arrange({{ id }}, {{ time }}) %>%
    group_by({{ id }}) %>%
    mutate(
      t0 = {{ time }},
      t1 = lead({{ time }}),
      s0 = {{ state }},
      s1 = lead({{ state }}),
      dt = t1 - t0
    ) %>%
    ungroup() %>%
    filter(!is.na(t1), dt > 0)
}

# Identifiability diagnostics

transition_diagnostics <- function(df_intervals,
                                   age_cuts,
                                   min_exposure = 20,
                                   min_events = 5,
                                   min_eff_exp = 8) {
  # Stable baseline estimation requires:
  #   • sufficient exposure
  #   • enough events
  #   • adequate effective sample size

  effective_n <- function(x) {
    s <- sum(x)
    if (s == 0) return(0)
    s^2 / sum(x^2)
  }

  df <- df_intervals %>%
    mutate(
      age_band = cut(
        t0,
        breaks = age_cuts,
        right = FALSE,
        include.lowest = TRUE
      ),
      transition = paste0(s0, "→", s1)
    )

  info_map <- df %>%
    group_by(transition, age_band) %>%
    summarise(
      events = n(),
      exposure = sum(dt),
      eff_exp = effective_n(dt),
      .groups = "drop"
    ) %>%
    mutate(
      rate = ifelse(exposure > 0, events / exposure, NA_real_),
      flag = case_when(
        exposure == 0 ~ "NO RISK SET",
        exposure < min_exposure |
          events < min_events |
          eff_exp < min_eff_exp ~ "LOW INFO",
        TRUE ~ ""
      ),
      exposure = round(exposure, 1),
      eff_exp = round(eff_exp, 1),
      rate = round(rate, 3),
      info = case_when(
        flag == "" ~
          paste0(
            events, " / ", exposure,
            " (", rate, ", eff=", eff_exp, ")"
          ),
        TRUE ~
          paste0(
            events, " / ", exposure,
            " (", rate, ", eff=", eff_exp, ") ⚠ ", flag
          )
      )
    ) %>%
    select(-events, -exposure, -rate, -eff_exp, -flag) %>%
    pivot_wider(
      names_from = age_band,
      values_from = info,
      values_fill = "0 / 0 (eff=0) ⚠ NO RISK SET"
    ) %>%
    arrange(
      as.numeric(sub("→.*", "", transition)),
      as.numeric(sub(".*→", "", transition))
    )

  info_map
}

# Baseline estimators

# Rolling empirical hazard (quasi-nonparametric)

rolling_empirical_hazard <- function(df_intervals,
                                     from,
                                     to,
                                     bandwidth = 5,
                                     step = 0.5,
                                     min_exposure = 50,
                                     min_events = 5) {
  df <- df_intervals %>%
    filter(s0 == from) %>%
    mutate(
      event = as.integer(s1 == to),
      t_mid = (t0 + t1) / 2
    )

  # Hard risk-set check
  if (nrow(df) == 0 || sum(df$dt) == 0)
    return(NULL)

  grid <- seq(floor(min(df$t_mid)),
    ceiling(max(df$t_mid)),
    by = step
  )

  out <- purrr::map_dfr(grid, function(a) {
    sub <- df %>%
      filter(
        t_mid >= a - bandwidth / 2,
        t_mid < a + bandwidth / 2
      )

    exposure <- sum(sub$dt)
    events <- sum(sub$event)

    if (exposure < min_exposure || events < min_events)
      return(NULL)

    ci <- poisson.test(events, T = exposure)$conf.int

    tibble(
      age = a,
      hazard = events / exposure,
      lo = ci[1],
      hi = ci[2]
    )
  })

  if (nrow(out) == 0) return(NULL)

  out
}

# Piecewise baseline via cloglog surrogate

pw_baseline <- function(df_intervals, from, to, age_cuts) {
  df <- df_intervals %>%
    filter(s0 == from) %>%
    mutate(
      event = as.integer(s1 == to),
      age_mid = (t0 + t1) / 2,
      age_band = cut(age_mid,
        breaks = age_cuts,
        right = FALSE,
        include.lowest = TRUE
      )
    ) %>%
    drop_na(age_band)

  if (nrow(df) == 0 || sum(df$dt) == 0)
    return(NULL)

  df$age_band <- droplevels(df$age_band)

  fit <- tryCatch(
    glm(
      event ~ age_band + offset(log(dt)),
      family = binomial(link = "cloglog"),
      data = df
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  bands <- levels(df$age_band)

  pred_df <- tibble(age_band = bands, dt = 1)

  pr <- predict(fit,
    newdata = pred_df,
    se.fit = TRUE,
    type = "link"
  )

  pred <- tibble(age_band = bands) %>%
    mutate(
      left = purrr::map_dbl(age_band, ~ band_lr(.x)["left"]),
      right = purrr::map_dbl(age_band, ~ band_lr(.x)["right"]),
      hazard = exp(pr$fit),
      lo = exp(pr$fit - 1.96 * pr$se.fit),
      hi = exp(pr$fit + 1.96 * pr$se.fit)
    ) %>%
    arrange(left)

  list(fit = fit, pred = pred)
}

# Gompertz baseline

gompertz_baseline <- function(df_intervals, from, to) {
  df <- df_intervals %>%
    filter(s0 == from) %>%
    mutate(
      event = as.integer(s1 == to),
      t_mid = (t0 + t1) / 2
    )

  if (nrow(df) == 0 || sum(df$dt) == 0)
    return(NULL)

  tryCatch(
    glm(
      event ~ t_mid + offset(log(dt)),
      family = poisson(link = "log"),
      data = df
    ),
    error = function(e) NULL
  )
}

# Master diagnostic plot

plot_transition_baseline <- function(df_intervals,
                                     from,
                                     to,
                                     age_cuts,
                                     bandwidth = 5,
                                     zoom_empirical = FALSE) {
  df_sub <- df_intervals %>% filter(s0 == from)

  if (nrow(df_sub) == 0 || sum(df_sub$dt) == 0) {
    message("No risk set for ", from, "->", to)
    return(NULL)
  }

  grid <- seq(floor(min(df_sub$t0)),
    ceiling(max(df_sub$t1)),
    by = 0.25
  )

  empirical <- rolling_empirical_hazard(df_intervals, from, to, bandwidth)
  pw <- pw_baseline(df_intervals, from, to, age_cuts)

  gfit <- gompertz_baseline(df_intervals, from, to)
  gompertz <- if (!is.null(gfit)) predict_gompertz(gfit, grid) else NULL

  if (is.null(empirical) & is.null(pw) & is.null(gompertz))
    return(NULL)

  g <- ggplot() +
    labs(
      title = paste0("Transition ", from, " -> ", to),
      y = "Hazard (per year)",
      x = "Age"
    ) +
    scale_y_log10() +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "left",
      panel.grid.minor = element_blank()
    )

  if (!is.null(gompertz)) {
    g <- g +
      geom_ribbon(
        data = gompertz,
        aes(age, ymin = lo, ymax = hi, fill = "Gompertz CI"),
        alpha = .18
      ) +
      geom_line(
        data = gompertz,
        aes(age, hazard, color = "Gompertz"),
        linewidth = 1.2
      )
  }

  if (!is.null(pw)) {
    pw_step <- pw$pred %>%
      select(left, right, hazard) %>%
      tidyr::pivot_longer(c(left, right),
        names_to = "which",
        values_to = "age"
      ) %>%
      arrange(age)

    g <- g +
      geom_step(
        data = pw_step,
        aes(x = age, y = hazard, color = "Piecewise"),
        linewidth = 1.1
      )
  }

  if (!is.null(empirical)) {
    g <- g +
      geom_ribbon(
        data = empirical,
        aes(age,
          ymin = lo, ymax = hi,
          fill = "Empirical CI"
        ),
        alpha = .15
      ) +
      geom_line(
        data = empirical,
        aes(age, hazard, color = "Empirical"),
        linewidth = 1.3
      )
  }

  g_final <- g +
    scale_color_manual(
      values = c(
        "Empirical" = "black",
        "Piecewise" = "#1B9E77",
        "Gompertz"  = "#D95F02"
      ),
      name = NULL
    ) +
    scale_fill_manual(
      values = c(
        "Empirical CI" = "black",
        "Gompertz CI" = "#D95F02"
      ),
      name = NULL
    )

  if (zoom_empirical && !is.null(empirical)) {
    xmin <- min(empirical$age)
    xmax <- max(empirical$age)

    g_final <- g_final +
      coord_cartesian(xlim = c(xmin, xmax))
  }

  g_final
}

# Covariate stratification

# Prepare grouping variable (auto-detect discrete vs continuous)

prepare_covariate <- function(df, covariate, cov_type, cov_bins) {
  v <- df[[covariate]]

  if (cov_type == "auto") {
    if (is.numeric(v) && dplyr::n_distinct(v, na.rm = TRUE) > 8) {
      cov_type <- "continuous"
    } else {
      cov_type <- "discrete"
    }
  }

  if (cov_type == "continuous") {
    qs <- quantile(
      v,
      probs = seq(0, 1, length.out = cov_bins + 1),
      na.rm = TRUE
    )

    qs <- unique(qs)

    # se tutti i quantili coincidono, cut() fallisce
    if (length(qs) < 2) {
      return(
        df %>%
          mutate(.cov_group = factor("all"))
      )
    }

    df %>%
      mutate(
        .cov_group = cut(
          .data[[covariate]],
          breaks = qs,
          include.lowest = TRUE
        )
      )
  } else {
    df %>%
      mutate(
        .cov_group = factor(.data[[covariate]])
      )
  }
}

# Stratified diagnostic plot wrapper

plot_transition_baseline_cov <- function(df_intervals,
                                         from,
                                         to,
                                         age_cuts,
                                         covariate = NULL,
                                         cov_type = "auto",
                                         cov_bins = 4,
                                         bandwidth = 5) {
  if (is.null(covariate)) {
    return(
      plot_transition_baseline(
        df_intervals,
        from,
        to,
        age_cuts,
        bandwidth
      )
    )
  }

  df <- prepare_covariate(
    df_intervals,
    covariate,
    cov_type,
    cov_bins
  ) %>%
    dplyr::filter(!is.na(.cov_group))

  groups <- unique(df$.cov_group)

  plots <- purrr::map(groups, function(g) {
    p <- plot_transition_baseline(
      df %>% dplyr::filter(.cov_group == g),
      from,
      to,
      age_cuts,
      bandwidth
    )

    if (!is.null(p)) {
      p + ggtitle(
        paste0(
          "Transition ", from, " -> ", to,
          " | covariate = ", covariate,
          " | group = ", g
        )
      )
    } else {
      NULL
    }
  })

  plots <- purrr::compact(plots)

  if (length(plots) == 0) return(NULL)

  patchwork::wrap_plots(plots)
}

# Trim panel

trim_panel_sequential <- function(df,
                                  id,
                                  time,
                                  lower = NULL,
                                  upper = NULL,
                                  truncate_after_upper = TRUE) {
  id <- rlang::enquo(id)
  time <- rlang::enquo(time)

  df %>%
    arrange(!!id, !!time) %>%
    group_by(!!id) %>%
    group_modify(~ {
      patient <- .x
      times <- dplyr::pull(patient, !!time)
      n <- length(times)

      # nothing to do if <2 visits
      if (n < 2) return(patient[0, ])

      keep <- logical(n)
      keep[1] <- TRUE

      last_kept <- 1

      for (i in 2:n) {
        dt <- times[i] - times[last_kept]

        # LOWER bound
        if (!is.null(lower) && dt < lower) {
          next
        }

        # UPPER bound
        if (!is.null(upper) && dt > upper) {
          # administratively censor follow-up
          if (truncate_after_upper) {
            break
          } else {
            next
          }
        }

        keep[i] <- TRUE
        last_kept <- i
      }

      patient[keep, ]
    }) %>%
    ungroup() %>%
    group_by(!!id) %>%
    filter(n() >= 2) %>%
    ungroup()
}

# Statetable by interval

statetable_by_interval <- function(data,
                                   id,
                                   state,
                                   time,
                                   covariate) {
  id <- rlang::ensym(id)
  st <- rlang::ensym(state)
  tm <- rlang::ensym(time)
  cv <- rlang::ensym(covariate)

  df_intervals <- data %>%
    arrange(!!id, !!tm) %>%
    group_by(!!id) %>%
    mutate(
      s0 = !!st,
      s1 = lead(!!st),
      cov0 = !!cv, # valore a inizio intervallo
      t1 = lead(!!tm)
    ) %>%
    filter(!is.na(s1)) %>%
    ungroup()

  tab <- with(
    df_intervals,
    table(s0, s1, cov0)
  )

  return(list(
    intervals = df_intervals,
    statetable = tab
  ))
}

################################################################################
# Helper: write TXT summary
write_transition_txt <- function(path,
                                 from,
                                 to,
                                 diag,
                                 df_panel) {
  file <- file.path(path, "transition_summary.txt")

  sink(file)

  cat("====================================\n")
  cat("Transition:", from, "->", to, "\n")
  cat("====================================\n\n")

  cat("IDENTIFIABILITY\n")
  print(
    diag %>%
      filter(transition == paste0(from, "→", to)),
    width = Inf
  )

  cat("\n\nSTATE TABLE\n")
  print(
    statetable.msm(
      state = state,
      subject = patient_id,
      data = df_panel
    )
  )

  sink()
}

################################################################################
# plot

hazard_table_msm <- function(fit, digits = 2) {
  haz <- hazard.msm(fit)

  df <- imap_dfr(haz, function(tab, cov) {
    tab %>%
      as.data.frame() %>%
      rownames_to_column("transition") %>% # <<< FIX VERO
      mutate(covariate = cov)
  })

  # remove baseline-only transitions
  df <- df %>%
    filter(!(HR == 1 & L == 1 & U == 1)) %>%
    arrange(transition) %>%
    relocate(covariate, .before = HR)

  # nicer transition format
  df <- df %>%
    mutate(
      transition = str_replace(
        transition,
        "State ([0-9]+) - State ([0-9]+)",
        "\\1$\\\\rightarrow$\\2"
      ),
      sig = ifelse(L > 1 | U < 1, "*", ""),
      HR_CI = sprintf(
        paste0("%.", digits, "f (%.", digits, "f--%.", digits, "f)"),
        HR, L, U
      ),
      covariate = paste0("\\texttt{", covariate, "}")
    ) %>%
    relocate(HR_CI, .before = sig)

  tab <- df %>%
    select(transition, covariate, HR_CI, sig) %>%
    kable(
      format = "latex",
      booktabs = TRUE,
      longtable = TRUE,
      escape = FALSE,
      col.names = c(
        "Transition",
        "Covariate",
        "HR (95\\% CI)",
        "Sig."
      )
    ) %>%
    collapse_rows(columns = 1, latex_hline = "major") %>%
    kable_styling(latex_options = "repeat_header")

  return(tab)
}
