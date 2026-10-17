#!/usr/bin/env Rscript
# ============================================================================
# v28.4_sa3_sa7.R
# SA3: DAG13（剔除 sepsis3），m=20 Rubin，SEQTaRget hazard，与主分析同口径
# SA7: 完整病例（原始未插补 DAG14 齐全者），无 MI，SEQTaRget hazard + bootstrap
#
# 用法：
#   Rscript v28.4_sa3_sa7.R --which dag13 --imps 1,...,20 --nboot 200
#   Rscript v28.4_sa3_sa7.R --which cca   --nboot 200
#   Rscript v28.4_sa3_sa7.R --which both  ...
# 版本 V28.4  日期 2026-10-10
# ============================================================================
CODE_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r"
source(file.path(CODE_DIR, "v28.4_common.R"))
suppressPackageStartupMessages(library(SEQTaRget))

parse_args <- function(argv) {
  getv <- function(key, def = NULL) {
    i <- match(key, argv); if (is.na(i)) def else argv[i + 1]
  }
  list(which = getv("--which", "both"),
       imps  = as.integer(strsplit(getv("--imps", "1"), ",")[[1]]),
       nboot = as.integer(getv("--nboot", 200)))
}
ARGS <- parse_args(commandArgs(trailingOnly = TRUE))

# ---- 原始（未插补）stay 唯一试验行 ----
load_raw_trials <- function() {
  aug <- as.data.table(read_parquet(AUG_FILE))
  setorder(aug, stay_id, trial_num)
  dt <- aug[, .SD[1], by = stay_id]
  out <- dt[, .(stay_id, trial_num, time_zero, itt_group,
                death_28day_num, surv_time)]
  out[, origin_day := as.integer(as.Date(time_zero) -
                                   min(as.Date(time_zero))), by = stay_id]
  out[, arm := fifelse(itt_group == "treatment", 1L, 0L)]
  for (v in DAG_FROM_IMP) out[, (v) := dt[[v]]]
  for (v in DAG_FROM_AUG) out[, (v) := dt[[v]]]
  for (v in BIN_COLS)    out[, (v) := to01(out[[v]])]
  out[, gender := factor(gender, levels = GENDER_LEV)]
  out[, admission_type := factor(admission_type, levels = ADMIT_LEV)]
  out[, nutrition_support_type := factor(nutrition_support_type,
                                         levels = NUTR_LEV)]
  out
}

run_seq <- function(trials, covset, label) {
  daily <- build_daily_input(trials)
  check_daily_input(daily, label = label)
  opts <- SEQopts(hazard = TRUE, km.curves = FALSE,
                  bootstrap = TRUE, bootstrap.nboot = ARGS$nboot,
                  bootstrap.sample = 1, followup.max = TAU,
                  seed = SEED, bootstrap.CI_method = "se")
  set.seed(SEED)
  model <- SEQuential(daily, id.col = "stay_id", time.col = "day",
                      eligible.col = "eligible", treatment.col = "tx_init",
                      outcome.col = "outcome", time_varying.cols = covset,
                      fixed.cols = list(), method = "ITT",
                      options = opts, verbose = TRUE)
  h0 <- hazard_ratio(model)
  cat(label, "HR:\n"); print(h0)
  extract_log(h0)
}

OUT_DIR <- file.path(V28_BASE,
                     "TTE/v28.4/output/tables/nonprimary_redo")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
SA3_PER_F <- file.path(OUT_DIR, "v28.4_sa3_dag13_per_imp.csv")

# ---- SA3 ----
if (ARGS$which %in% c("dag13", "both")) {
  DAG13 <- setdiff(DAG14, "sepsis3")
  IMPS <- ARGS$imps
  done_imps <- if (file.exists(SA3_PER_F)) {
    x <- fread(SA3_PER_F, colClasses = list(imp = "integer"))$imp
    sort(unique(x))
  } else integer(0)
  todo <- setdiff(IMPS, done_imps)
  if (length(done_imps) > 0 && !setequal(done_imps, intersect(done_imps, IMPS))) {
    cat(sprintf("[checkpoint] SA3 PER_F imps=%s ⊄ IMPS → reset\n",
                paste(done_imps, collapse = ",")))
    file.remove(SA3_PER_F); done_imps <- integer(0); todo <- IMPS
  }
  cat(sprintf("[SA3 checkpoint] todo=%s done=%s\n",
              if (length(todo)) paste(todo, collapse = ",") else "(none)",
              if (length(done_imps)) paste(done_imps, collapse = ",") else "(none)"))
  for (im in todo) {
    cat(sprintf("\n==== SA3 imp %d ====\n", im))
    r <- run_seq(load_trials(im), DAG13, paste0("dag13_imp", im))
    per_row <- data.table(imp = im, loghr = r$loghr, se = r$se,
                          hr = r$hr, lci = r$lci, uci = r$uci)
    fwrite(per_row, SA3_PER_F, append = file.exists(SA3_PER_F))
    rm(r, per_row); invisible(gc())
  }
  per_dt <- fread(SA3_PER_F, colClasses = list(imp = "integer"))
  stopifnot(setequal(sort(per_dt$imp), IMPS))
  loghr <- per_dt$loghr; sev <- per_dt$se
  pooled <- rubin_pool(loghr, sev, m = length(IMPS))
  cat(sprintf("\nSA3 DAG13 HR = %.3f (%.3f-%.3f), p=%.4f\n",
              pooled$hr, pooled$lci, pooled$uci, pooled$p))
  data.table(analysis = "SA3 DAG13", m = length(IMPS),
             n_boot = ARGS$nboot, hr = pooled$hr, lci = pooled$lci,
             uci = pooled$uci, p = pooled$p) |>
    fwrite(file.path(OUT_DIR, "v28.4_sa3_dag13.csv"))
}

# ---- SA7 ----
if (ARGS$which %in% c("cca", "both")) {
  raw <- load_raw_trials()
  comp <- raw[complete.cases(raw[, ..DAG14])]
  cat("\nSA7 complete-case n =", nrow(comp),
      "events =", sum(comp$death_28day_num), "\n")
  r <- run_seq(comp, DAG14, "cca")
  data.table(analysis = "SA7 complete-case", n = nrow(comp),
             n_event = sum(comp$death_28day_num), n_boot = ARGS$nboot,
             hr = r$hr, lci = r$lci, uci = r$uci,
             p = 2 * pnorm(-abs(r$loghr / r$se))) |>
    fwrite(file.path(OUT_DIR, "v28.4_sa7_cca.csv"))
  cat(sprintf("SA7 CCA HR = %.3f (%.3f-%.3f)\n", r$hr, r$lci, r$uci))
}
