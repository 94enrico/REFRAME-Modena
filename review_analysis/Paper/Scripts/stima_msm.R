# Optional refit from review_analysis/, e.g. Rscript Paper/Scripts/stima_msm.R 3
suppressPackageStartupMessages(library(msm))
source("Paper/Scripts/prepara_msm.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, args %in% as.character(0:3))
n <- as.integer(args[1])
d <- prepara("metodi")
old <- starting_values(n)
edges <- c("1-2", "1-3", "1-4", "2-1", "2-3", "2-4", "3-1", "3-2", "3-4")
alive <- switch(as.character(n),
  "0" = NULL,
  "1" = ~access_age,
  "2" = ~ access_age + hiv_tertile + phys_act_bin,
  "3" = ~ access_age + hiv_tertile + phys_act_bin + INSTI + TAF + IP + comorbidity_Dislipidemia
)
forms <- if (n == 0) NULL else setNames(lapply(edges, function(e) if (endsWith(e, "4")) ~ sex + access_age else alive), edges)
constraints <- if (n == 0) NULL else list(access_age = c(1, 2, 3, 4, 5, 3, 6, 7, 3))
q <- unname(old$Qmatrices$baseline)
covinit <- NULL
if (n > 0) {
  covinit <- lapply(old$qcmodel$covlabels, function(v) {
    mat <- old$Qmatrices[[v]]
    idx <- do.call(rbind, strsplit(edges, "-"))
    values <- mat[cbind(as.integer(idx[, 1]), as.integer(idx[, 2]))]
    # Supply initial slopes only for intensities on which this covariate is estimated.
    as.numeric(values[old$qcmodel$cri[, v] == 1])
  })
  names(covinit) <- old$qcmodel$covlabels
}
scale <- old$scale
stopifnot(length(scale) == c(9, 19, 37, 61)[n + 1], all(is.finite(scale)), all(scale > 0))
ctl <- list(
  parscale = scale, fnscale = 1, maxit = 700, reltol = 1e-12,
  trace = 1, REPORT = 5, ndeps = rep(1e-3, length(scale))
)
fit <- msm(state ~ access_age,
  subject = patient_id, data = d, qmatrix = q,
  covariates = forms, covinits = covinit, constraint = constraints,
  gen.inits = FALSE, deathexact = TRUE, opt.method = "optim", method = "BFGS",
  control = ctl, hessian = TRUE
)
stopifnot(
  fit$opt$convergence == 0, isTRUE(fit$foundse),
  is.finite(fit$minus2loglik), all(is.finite(fit$opt$hessian))
)
eig <- eigen(fit$opt$hessian, symmetric = TRUE, only.values = TRUE)$values
scaled.score <- max(abs(fit$deriv) * scale)
stopifnot(min(eig) > 0, scaled.score < 0.05)
print(c(
  model = n, minus2loglik = fit$minus2loglik, AIC = AIC(fit), parameters = fit$paramdata$nopt,
  min_hessian_eigenvalue = min(eig), max_scaled_score = scaled.score
))
dir.create("Paper/Modelli", recursive = TRUE, showWarnings = FALSE)
saveRDS(fit, file.path("Paper/Modelli", sprintf("reframe_fit%d.rds", n)))
