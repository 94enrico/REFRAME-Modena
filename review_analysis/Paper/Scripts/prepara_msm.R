# Prepare the same data used by the fitted models. Run from review_analysis/.
prepara <- function(coorte) {
  file <- tempfile(fileext = ".csv")
  on.exit(unlink(file))
  python <- Sys.getenv("PYTHON", unset = "python")
  status <- system2(python, c("Paper/Scripts/coorti.py", coorte, shQuote(file)))
  stopifnot(status == 0L)
  d <- read.csv(file, check.names = FALSE)
  d <- d[order(d$patient_id, d$access_age), ]
  for (v in c("liver_disease", "insulin_resistance", "sppb_low")) d[[v]][is.na(d[[v]])] <- 0
  # VIM's original Gower-distance calculation used a numeric patient identifier,
  # numeric clinical fields and a character sex field. Keep those data types here.
  d$patient_id <- as.numeric(d$patient_id)
  d$hiv_duration_yrs <- ave(d$hiv_duration_yrs, d$patient_id, FUN = function(x) {
    observed <- x[!is.na(x)]
    if (length(observed)) rep(observed[1], length(x)) else x
  })
  set.seed(52001)
  d <- VIM::kNN(d, variable = "phys_act_bin", k = 5, imp_var = FALSE)
  d <- VIM::kNN(d, variable = "hiv_duration_yrs", k = 5, imp_var = FALSE)
  d$sex <- as.integer(d$sex == "M")
  baseline <- ave(d$hiv_duration_yrs, d$patient_id, FUN = function(x) rep(x[1], length(x)))
  cuts <- quantile(baseline, c(0, 1 / 3, 2 / 3, 1))
  d$hiv_tertile <- cut(baseline, cuts, include.lowest = TRUE, labels = c("Short", "Medium", "Long"))
  stopifnot(!anyNA(d), all(table(d$patient_id) >= 2), all(d$state %in% 1:4))
  d
}

# Either use the compact private starts or derive them from the matching historical fit.
starting_values <- function(model) {
  stored <- "Dati/valori_iniziali_msm.rds"
  if (file.exists(stored)) return(readRDS(stored)[[model + 1L]])

  path <- file.path("../original_analysis/6_Models", sprintf("reframe_fit%d.rds", model))
  original <- readRDS(path)
  list(
    Qmatrices = original$Qmatrices,
    qcmodel = original$qcmodel,
    scale = sqrt(diag(solve(original$opt$hessian)))
  )
}
