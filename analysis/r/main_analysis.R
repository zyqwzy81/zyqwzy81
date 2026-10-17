#!/usr/bin/env Rscript
# ============================================================================
# v28.4_main_analysis.R
# 主分析：SEQTaRget ITT，hazard=TRUE，DAG14 校正，估计 IV 磷剂对 28d 死亡的 HR
#
# 用法：
#   pilot:  Rscript v28.4_main_analysis.R --pilot --nstays 3000 --imps 1,2 --nboot 20
#   全量:   Rscript v28.4_main_analysis.R --imps 1,2,...,20 --nboot 200
#   分批并行 (默认): 同全量参数，另加 --batch-size 50 --ncores 4
#
#   分批原理：点估计 1×SEQuential(bootstrap=FALSE) + n_batch×batch_size 次 bootstrap
#             各批 seed 递增偏移，SE_combined = sqrt(mean(SE_batch²))
#             实验已验证 imp01 B=20 vs 4×B=5: SE 差异仅 1.3%
#
# 须以 LD_LIBRARY_PATH=/opt/miniforge3/lib 运行（SEQTaRget 依赖 MKL）
# 版本 V28.4  日期 2026-10-09
# ============================================================================

CODE_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r"
source(file.path(CODE_DIR, "v28.4_common.R"))
suppressPackageStartupMessages({
  library(SEQTaRget)
  library(parallel)
})

# ---- 命令行参数 ----
parse_args <- function(argv) {
  getv <- function(key, def = NULL) {
    i <- match(key, argv)
    if (is.na(i)) def else argv[i + 1]
  }
  list(pilot      = "--pilot" %in% argv,
       nstays     = as.integer(getv("--nstays", NA)),
       imps       = as.integer(strsplit(getv("--imps", "1"), ",")[[1]]),
       nboot      = as.integer(getv("--nboot", 200)),
       batch_size = as.integer(getv("--batch-size", 50)),
       ncores     = as.integer(getv("--ncores", 4)),
       no_batch   = "--no-batch" %in% argv,
       tag        = getv("--tag", "main"))
}
ARGS <- parse_args(commandArgs(trailingOnly = TRUE))

run_one <- function(imp_id) {
  cat("\n================ imp", imp_id, "================\n")
  trials <- load_trials(imp_id)
  if (!is.na(ARGS$nstays)) {
    keep <- sort(sample(unique(trials$stay_id), ARGS$nstays, replace = FALSE))
    trials <- trials[stay_id %in% keep]
  }
  daily <- build_daily_input(trials)
  check_daily_input(daily, label = paste0("imp", imp_id))

  if (ARGS$no_batch) {
    # === 原方案: 单次 SEQuential(bootstrap=TRUE, B=nboot) ===
    cat("[mode] single-call bootstrap B=", ARGS$nboot, "\n")
    opts <- SEQopts(hazard = TRUE, km.curves = FALSE, bootstrap = TRUE,
                    bootstrap.nboot = ARGS$nboot, bootstrap.sample = 1,
                    followup.max = TAU, seed = SEED, bootstrap.CI_method = "se")
    set.seed(SEED)
    model <- SEQuential(daily, id.col = "stay_id", time.col = "day",
                        eligible.col = "eligible", treatment.col = "tx_init",
                        outcome.col = "outcome",
                        time_varying.cols = DAG14, fixed.cols = list(),
                        method = "ITT", options = opts, verbose = TRUE)
    h0 <- hazard_ratio(model)
    cat("imp", imp_id, "HR result:\n"); print(h0)
    extract_log(h0)
  } else {
    # === 分批并行方案 ===
    run_one_imp_batch_bootstrap(
      daily, DAG14 = DAG14,
      n_boot = ARGS$nboot, batch_size = ARGS$batch_size,
      base_seed = SEED, ncores = ARGS$ncores, verbose = TRUE)
  }
}

set.seed(SEED)
per_imp <- lapply(ARGS$imps, run_one)
loghr <- vapply(per_imp, function(x) x$loghr, numeric(1))
sev   <- vapply(per_imp, function(x) x$se, numeric(1))

pooled <- rubin_pool(loghr, sev, m = length(ARGS$imps))
cat("\n========== RUBIN POOLED (m =", length(ARGS$imps), ") ==========\n")
cat(sprintf("HR = %.3f (%.3f - %.3f), p = %.4f\n",
            pooled$hr, pooled$lci, pooled$uci, pooled$p))

# ---- 落盘 ----
OUT_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v28.4/output/tables/main"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
suffix <- if (ARGS$pilot) paste0("_pilot_", ARGS$tag) else ""
if (!ARGS$no_batch) suffix <- paste0(suffix, "_batch", ARGS$batch_size)
out_csv <- file.path(OUT_DIR, paste0("v28.4_main", suffix, ".csv"))

per_dt <- data.frame(imp = ARGS$imps,
                     loghr = loghr, se = sev,
                     hr = vapply(per_imp, function(x) x$hr, numeric(1)),
                     lci = vapply(per_imp, function(x) x$lci, numeric(1)),
                     uci = vapply(per_imp, function(x) x$uci, numeric(1)))
write.csv(per_dt, out_csv, row.names = FALSE)

mode_tag <- if (ARGS$no_batch) "single" else paste0("batch-", ARGS$batch_size, "x", ARGS$nboot/ARGS$batch_size)
pool_dt <- data.frame(version = "V28.4-batch", date = format(Sys.Date()),
                      mode = mode_tag, pilot = ARGS$pilot, m = length(ARGS$imps),
                      n_boot = ARGS$nboot, batch_size = ARGS$batch_size,
                      ncores = if (ARGS$no_batch) 1 else ARGS$ncores,
                      hr = pooled$hr, lci = pooled$lci, uci = pooled$uci,
                      p = pooled$p, B = pooled$B, df = pooled$df)
write.csv(pool_dt, file.path(OUT_DIR, paste0("v28.4_pooled", suffix, ".csv")),
          row.names = FALSE)
cat("\nresults saved to:", out_csv, "\n")
