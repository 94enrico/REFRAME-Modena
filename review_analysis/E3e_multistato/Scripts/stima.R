# Optional refit from review_analysis/, e.g. Rscript E3e_multistato/Scripts/stima.R M3
suppressPackageStartupMessages(library(msm))
source("Paper/Scripts/prepara_msm.R")
model <- commandArgs(trailingOnly = TRUE)
stopifnot(length(model) == 1L, model %in% paste0("M", 1:4))
d <- prepara(if (model %in% c("M1", "M2")) "storica" else "metodi")
archived <- starting_values(3L)
common <- model %in% c("M2", "M4")
fixed <- model == "M3"
stopifnot(!is.unsorted(d$patient_id), all(table(d$patient_id) >= 2))

edges <- c("1-2", "1-3", "1-4", "2-1", "2-3", "2-4", "3-1", "3-2", "3-4")
live <- ~ access_age + hiv_tertile + phys_act_bin + INSTI + TAF + IP + comorbidity_Dislipidemia
forms <- setNames(lapply(edges, function(e) if (endsWith(e, "4")) ~ sex + access_age else live), edges)
constraints <- list(access_age = c(1, 2, 3, 4, 5, 3, 6, 7, 3))
if (common) constraints$sex <- c(1, 1, 1)
if (fixed) forms[["2-4"]] <- ~access_age

# Use the original fit only to initialise the baseline intensities and slopes.
ij <- do.call(rbind, strsplit(edges, "-"))
index <- cbind(as.integer(ij[, 1]), as.integer(ij[, 2]))
included <- archived$qcmodel$cri == 1
if (fixed) included["2-4", "sex"] <- FALSE
covfrom <- function(v) unname(archived$Qmatrices[[v]][index][included[, v]])
initial <- setNames(lapply(colnames(included), covfrom), colnames(included))
if (common) initial$sex[] <- 0
rate <- unname(archived$Qmatrices$baseline)

# Scale each coefficient by its numerical curvature to avoid stopping early
# on the almost flat likelihood for mortality from preclinical obesity.
scales <- archived$scale
stopifnot(length(scales) == 61)
if (common) scales <- c(head(scales, -3), 0.25)
if (fixed) scales <- scales[-60]
fit_once <- function(shift, hessian) {
  covinit <- initial
  covinit$sex <- covinit$sex + shift
  msm(state ~ access_age,
    subject = patient_id, data = d, qmatrix = rate,
    covariates = forms, constraint = constraints, covinits = covinit,
    gen.inits = FALSE, deathexact = TRUE, opt.method = "optim", method = "BFGS",
    hessian = hessian, control = list(
      parscale = scales, fnscale = 1,
      reltol = 1e-13, maxit = 700, ndeps = rep(0.001, length(scales))
    )
  )
}
check <- fit_once(0.5, FALSE)
fit <- fit_once(0, TRUE)
score <- max(abs(fit$deriv) * scales)
discrepancy <- abs(fit$minus2loglik - check$minus2loglik)
eig <- eigen(fit$opt$hessian, symmetric = TRUE, only.values = TRUE)$values
stopifnot(
  check$opt$convergence == 0, fit$opt$convergence == 0, isTRUE(fit$foundse),
  min(eig) > 0, score < 0.005, discrepancy < 1e-5,
  max(abs(fit$opt$par - check$opt$par)) < 0.005
)
fit$reframe_checks <- c(
  max_scaled_score = score, min_hessian_eigenvalue = min(eig),
  difference_two_starts = discrepancy
)
dir.create("E3e_multistato/Modelli", recursive = TRUE, showWarnings = FALSE)
saveRDS(fit, file.path("E3e_multistato/Modelli", paste0(model, ".rds")))
cat(
  model, "people", length(unique(d$patient_id)), "observations", nrow(d),
  "deaths", sum(d$state == 4), "parameters", fit$paramdata$nopt,
  "-2logL", format(fit$minus2loglik, digits = 12), "AIC", format(AIC(fit), digits = 12), "\n"
)
print(hazard.msm(fit)$sex[c(3, 6, 9), , drop = FALSE])
