# Export the four nested paper fits and Supplementary Tables 3–5.
suppressPackageStartupMessages(library(msm))
out <- "Paper/Tables"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
models <- "Paper/Modelli"
edges <- c("1-2", "1-3", "1-4", "2-1", "2-3", "2-4", "3-1", "3-2", "3-4")
fits <- lapply(0:3, function(i) readRDS(file.path(models, sprintf("reframe_fit%d.rds", i))))
fit <- fits[[4]]
mf <- fit$data$mf
d <- data.frame(patient_id = mf[["(subject)"]], state = mf[["(state)"]], sex = mf$sex)
next.row <- which(diff(d$patient_id) == 0)
counts <- table(factor(d$state[next.row], levels = 1:3), factor(d$state[next.row + 1], levels = 1:4))
stopifnot(sum(counts) == nrow(d) - length(unique(d$patient_id)), sum(counts[, 4]) == sum(d$state == 4))
counts <- data.frame(From = 1:3, unclass(counts), check.names = FALSE)
names(counts) <- c("From", "To 1", "To 2", "To 3", "To death")
write.csv(counts, file.path(out, "Supplementary_table_3.csv"), row.names = FALSE)

# Display each active covariate at every transition on which it acts.
vlabels <- c(
  access_age = "Age (per year)", hiv_tertileMedium = "HIV duration, middle vs lowest tertile",
  hiv_tertileLong = "HIV duration, highest vs lowest tertile", phys_act_bin = "Physical activity",
  INSTI = "Current INSTI use", TAF = "Current TAF use", IP = "Current PI use",
  comorbidity_Dislipidemia = "Dyslipidemia", sex = "Male vs female sex"
)

get_hr <- function(i) {
  fit <- fits[[i + 1]]
  haz <- msm::hazard.msm(fit)
  idx <- seq_along(edges)
  rows <- do.call(rbind, lapply(names(haz), function(v) {
    j <- idx[fit$qcmodel$cri[, v] == 1]
    values <- haz[[v]][j, , drop = FALSE]
    estimable <- is.finite(values[, "U"]) & values[, "L"] > 0
    data.frame(
      Model = i, Transition = edges[j], Variable = v, Covariate = unname(vlabels[v]),
      HR = ifelse(estimable, values[, "HR"], NA),
      Lower95 = ifelse(estimable, values[, "L"], NA), Upper95 = ifelse(estimable, values[, "U"], NA),
      Status = ifelse(estimable, "Estimated", "Not estimable with finite Wald bounds"), row.names = NULL
    )
  }))
  rows
}

all_hr <- do.call(rbind, lapply(1:3, get_hr))
write.csv(all_hr, file.path(out, "HR.csv"), row.names = FALSE, na = "")
rows <- subset(all_hr, Model == 3)
age_at_death <- subset(rows, Variable == "access_age" & endsWith(Transition, "4"))
stopifnot(nrow(rows) == 54, max(abs(age_at_death$HR - age_at_death$HR[1])) < 1e-10)

render <- function(transitions, variables, filename) {
  result <- data.frame(Covariate = unname(vlabels[variables]))
  for (edge in transitions) {
    estimates <- rows[rows$Transition == edge, ]
    index <- match(variables, estimates$Variable)
    active <- which(!is.na(index))
    matched <- estimates[index[active], ]
    value <- rep("—", length(variables))
    value[active] <- ifelse(
      is.na(matched$HR), "NE",
      sprintf("%.2f (%.2f–%.2f)", matched$HR, matched$Lower95, matched$Upper95)
    )
    result[[gsub("-", "→", edge)]] <- value
  }
  write.csv(result, file.path(out, filename), row.names = FALSE, fileEncoding = "UTF-8")
}

render(c("1-2", "1-3", "1-4", "2-3", "2-4", "3-4"), names(vlabels), "Supplementary_table_4.csv")
render(c("2-1", "3-1", "3-2"), setdiff(names(vlabels), "sex"), "Supplementary_table_5.csv")
