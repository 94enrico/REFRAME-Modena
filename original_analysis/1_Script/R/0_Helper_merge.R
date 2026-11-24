merge_locf_tolerance <- function(
  panel_df,
  meas_df,
  tolerance_days = 30,
  type = c("visit", "event")
) {
  type <- match.arg(type)

  # Visit-based measurements use the encounter date when the measurement date is missing.

  if (type == "visit") {
    if (!all(c("patient_id", "access_date", "meas_date") %in% names(meas_df))) {
      stop("visit-type meas_df must contain: patient_id, access_date, meas_date")
    }

    meas_df <- meas_df %>%
      mutate(
        eff_meas_date = if_else(is.na(meas_date), access_date, meas_date)
      )

    panel_dates <- panel_df %>%
      select(patient_id, access_date) %>%
      rename(panel_access_date = access_date)

    joined <- panel_dates %>%
      left_join(
        meas_df %>% select(-access_date),
        by = "patient_id",
        relationship = "many-to-many"
      ) %>%
      mutate(
        diff_days = as.numeric(eff_meas_date - panel_access_date)
      )

    joined_valid <- joined %>%
      filter(!is.na(eff_meas_date), diff_days <= tolerance_days)

    joined_locf <- joined_valid %>%
      group_by(patient_id, panel_access_date) %>%
      slice_max(eff_meas_date, n = 1, with_ties = FALSE) %>%
      ungroup()

    out <- panel_df %>%
      rename(panel_access_date = access_date) %>%
      left_join(
        joined_locf %>%
          select(-diff_days, -meas_date),
        by = c("patient_id", "panel_access_date")
      ) %>%
      rename(access_date = panel_access_date)

    return(out)
  }

  # Event-based measurements (quality of life) have their own date.

  if (type == "event") {
    if (!all(c("patient_id", "meas_date") %in% names(meas_df))) {
      stop("event-type meas_df must contain: patient_id, meas_date")
    }

    panel_dates <- panel_df %>%
      select(patient_id, access_date) %>%
      rename(panel_access_date = access_date)

    joined <- panel_dates %>%
      left_join(
        meas_df,
        by = "patient_id",
        relationship = "many-to-many"
      ) %>%
      mutate(
        diff_days = as.numeric(meas_date - panel_access_date)
      )

    joined_valid <- joined %>%
      filter(!is.na(meas_date), diff_days <= tolerance_days)

    joined_locf <- joined_valid %>%
      group_by(patient_id, panel_access_date) %>%
      slice_max(meas_date, n = 1, with_ties = FALSE) %>%
      ungroup()

    out <- panel_df %>%
      rename(panel_access_date = access_date) %>%
      left_join(
        joined_locf %>% select(-diff_days),
        by = c("patient_id", "panel_access_date")
      ) %>%
      rename(access_date = panel_access_date)

    return(out)
  }
}
