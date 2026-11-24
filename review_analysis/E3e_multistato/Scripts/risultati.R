# Extract all four models from the saved fits. Run from review_analysis/.
suppressPackageStartupMessages(library(msm))
dest <- "E3e_multistato/Tables"
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
ids <- paste0("M", 1:4)
fits <- setNames(lapply(ids, function(id) readRDS(file.path("E3e_multistato/Modelli", paste0(id, ".rds")))), ids)
edges <- c("1-2", "1-3", "1-4", "2-1", "2-3", "2-4", "3-1", "3-2", "3-4")
cohort <- setNames(rep(c("Full follow-up", "Methods cohort"), each = 2), ids)
adjustment <- c(
  M1 = "Separate at each origin", M2 = "Common across origins",
  M3 = "Separate (preclinical fixed at 1)", M4 = "Common across origins"
)
labels <- c(
  access_age = "Age, per year", hiv_tertileMedium = "HIV duration, middle vs lowest tertile",
  hiv_tertileLong = "HIV duration, highest vs lowest tertile", phys_act_bin = "Physical activity",
  INSTI = "Current INSTI use", TAF = "Current TAF use", IP = "Current PI use",
  comorbidity_Dislipidemia = "Dyslipidemia", sex = "Men vs women"
)

check_pair <- function(i, j) {
  columns <- c("(subject)", "(state)", "(time)", "(obstype)", names(fits[[i]]$data$mf)[1:8])
  stopifnot(identical(fits[[i]]$data$mf[columns], fits[[j]]$data$mf[columns]))
}
check_pair("M1", "M2")
check_pair("M3", "M4")
model_summary <- do.call(rbind, lapply(ids, function(id) {
  f <- fits[[id]]
  mf <- f$data$mf
  data.frame(
    Model = id, Cohort = cohort[id], Sex_mortality = adjustment[id],
    People = length(unique(mf[["(subject)"]])), Observations = nrow(mf),
    Deaths = sum(mf[["(state)"]] == 4), Parameters = f$paramdata$nopt,
    Minus2logL = f$minus2loglik, AIC = AIC(f), row.names = NULL
  )
}))
model_summary$Delta_AIC_within_cohort <- with(model_summary, AIC - ave(AIC, Cohort, FUN = min))
stopifnot(
  all(model_summary$Parameters == c(61, 59, 60, 59)),
  all(abs(model_summary$Minus2logL + 2 * model_summary$Parameters - model_summary$AIC) < 1e-6)
)

all_hr <- do.call(rbind, lapply(ids, function(id) {
  f <- fits[[id]]
  hazards <- hazard.msm(f)
  do.call(rbind, lapply(names(labels), function(v) {
    idx <- if (v == "sex") which(endsWith(edges, "4")) else which(f$qcmodel$cri[, v] == 1)
    x <- hazards[[v]][idx, , drop = FALSE]
    is_fixed <- v == "sex" & id == "M3" & edges[idx] == "2-4"
    shared <- endsWith(edges[idx], "4") & (v == "access_age" | (v == "sex" & id %in% c("M2", "M4")))
    x[is_fixed, c("L", "U")] <- NA_real_
    stopifnot(
      all(x[is_fixed, "HR"] == 1), all(is.finite(x[!is_fixed, ])),
      all(x[!is_fixed, "L"] > 0), all(x[!is_fixed, "L"] <= x[!is_fixed, "HR"]),
      all(x[!is_fixed, "HR"] <= x[!is_fixed, "U"])
    )
    data.frame(
      Model = id, Transition = edges[idx], Covariate = labels[[v]],
      HR = x[, "HR"], Lower95 = x[, "L"], Upper95 = x[, "U"],
      Status = ifelse(is_fixed, "Fixed at 1", ifelse(shared, "Common to three death transitions", "Estimated separately")),
      row.names = NULL
    )
  }))
}))
stopifnot(nrow(all_hr) == 216, sum(all_hr$Status == "Fixed at 1") == 1)
for (id in ids) {
  common <- subset(all_hr, Model == id & Status == "Common to three death transitions")
  for (v in unique(common$Covariate)) {
    arr <- common[common$Covariate == v, c("HR", "Lower95", "Upper95")]
    stopifnot(nrow(unique(arr)) == 1)
  }
}

for (id in ids) {
  tab <- subset(all_hr, Model == id, select = -Model)
  write.csv(tab, file.path(dest, paste0("HR_", id, ".csv")), row.names = FALSE, na = "")
}
guide <- c(
  "Hazard ratios dei quattro modelli multistato",
  "",
  "Gli stati 1, 2, 3 e 4 indicano nessun eccesso di adiposità, obesità preclinica, obesità clinica e morte.",
  "Ogni riga riporta un HR con intervallo di Wald al 95% per una covariata e una transizione.",
  "Un HR superiore a 1 indica un’intensità di transizione maggiore per il confronto riportato.",
  "Status distingue gli effetti liberi, comuni alle tre transizioni alla morte o fissati.",
  "Per M3, l’HR uomini/donne per 2→4 è fissato a 1. Le celle del suo intervallo sono vuote.",
  "",
  "M1 e M2 usano la coorte storica. M3 e M4 usano la coorte dei Metodi (2005–2025, età ≥18 anni).",
  "M1 ha tre effetti del sesso sulla morte. In M3 sono liberi gli effetti 1→4 e 3→4.",
  "M2 e M4 hanno ciascuno un effetto del sesso comune alle tre transizioni alla morte.",
  "Tutti e quattro usano il criterio epatico corretto. Gli oggetti R sono nella cartella Modelli.",
  "",
  "Numerosità e AIC (confrontabile entro la stessa coorte)"
)
for (i in seq_len(nrow(model_summary))) {
  x <- model_summary[i, ]
  guide <- c(guide, sprintf(
    "%s   %d persone, %d osservazioni, %d decessi, %d parametri. AIC %.2f.",
    x$Model, x$People, x$Observations, x$Deaths, x$Parameters, x$AIC
  ))
}
guide <- c(
  guide, "",
  "Nella coorte dei Metodi sono state imputate 31 osservazioni per attività fisica e 339 per durata dell’HIV.",
  "L’intervallo mediano tra osservazioni è 364 giorni (Q1–Q3, 275–407). Il controllo è in Scripts/follow_up.py.",
  "Il calcolo storico degli intervalli in anni lascia fuori dieci visite a 30 giorni esatti per arrotondamento.",
  "Contare direttamente i giorni aggiungerebbe dieci osservazioni a parità di pazienti e decessi.",
  "Descrizione della verosimiglianza in msm. https://chjackson.github.io/msm/reference/msm.html"
)
writeLines(enc2utf8(guide), file.path(dest, "table_legend.txt"), useBytes = TRUE)
