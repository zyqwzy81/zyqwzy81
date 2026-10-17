#!/usr/bin/env Rscript
# ============================================================================
# v28.4_lasso_batch.R  — MITs-Lasso 协变量敏感性分析 (SA2)
# 复用 v28.4_resume_batch.R 架构, ncores=2 安全模式
# 仅差异: DAG14 → LASSO 31v (v26_final_covariates.json 32v - icu_type 缺失)
# 以及 trials merge v26 parquet 18 个数值型变量 (stay_id+trial_num)
# 用法: LD_LIBRARY_PATH=/opt/miniforge3/lib Rscript v28.4_lasso_batch.R 1,2,...,20
# ============================================================================
CODE_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r"
source(file.path(CODE_DIR, "v28.4_common.R"))
suppressPackageStartupMessages({library(SEQTaRget); library(parallel); library(arrow)})

BASE     <- "/data/小米云盘/文档/论文/血磷"
LASSO_JSON <- file.path(BASE, "sql/result/v26_final_covariates.json")
V26_PARQ_FMT <- file.path(BASE, "sql/result/mice_imputed",
                          "v26_mice_imputed_imp%02d.parquet")
OUT_DIR  <- file.path(BASE, "TTE/v28.4/output/tables/lasso")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_CSV  <- file.path(OUT_DIR, "v28.4_lasso31_batch50.csv")

NCORES_SAFE <- 2L

# ---- LASSO 32 全集 (从 JSON) ----
LASSO_FULL <- c("age","baseline_sofa","sofa_t0","sapsii","charlson","sepsis3",
  "baseline_phosphate","egfr","sodium_t0","albumin_min","bun_t0","inr_max",
  "ptt_max","ecg_abnormal","invasive_vent_at_t0","sf_ratio","vasoactive_at_t0",
  "steroid_at_t0","insulin_at_t0","nahco3_infusion_at_t0","transfusion_at_t0",
  "refeeding_risk","net_fluid_balance_24h","ever_crrt","lactate_t0",
  "admission_type","nutrition_support_type","icu_type","gender","race_group",
  "liver_disease","aki_kdigo2_plus")

# 从 v26 parquet merge 的数值变量 (非 factor, admission_type 等 trials 已有)
V26_MERGE_VARS <- c("sofa_t0","sapsii","sodium_t0","albumin_min","bun_t0","inr_max",
                    "ptt_max","ecg_abnormal","sf_ratio","steroid_at_t0","insulin_at_t0",
                    "nahco3_infusion_at_t0","transfusion_at_t0","net_fluid_balance_24h",
                    "ever_crrt","race_group","liver_disease","aki_kdigo2_plus")

# ============================================================================
# load_trials_lasso: 在 load_trials() 上 merge v26 parquet 18 变量
# ============================================================================
load_trials_lasso <- function(imp_id) {
  trials <- load_trials(imp_id)
  parquet_path <- sprintf(V26_PARQ_FMT, as.integer(imp_id))
  v26 <- as.data.table(read_parquet(parquet_path))
  # 只 merge LASSO 真正缺的 + v26 有的
  v26_has <- intersect(V26_MERGE_VARS, names(v26))
  setkey(trials, stay_id, trial_num)
  setkey(v26, stay_id, trial_num)
  trials_lasso <- merge(trials,
                        v26[, c("stay_id","trial_num", v26_has), with=FALSE],
                        by = c("stay_id","trial_num"), all.x = TRUE)
  # 最终 LASSO 可用变量 (排除 icu_type)
  LASSO_AVAIL <- intersect(LASSO_FULL, names(trials_lasso))
  # 排除 outcome 等非协变量
  LASSO_COV   <- setdiff(LASSO_AVAIL, c("outcome","stay_id","arm","origin_day",
                                         "death_28day_num","surv_time",
                                         "trial_num","itt_group","time_zero"))
  attr(trials_lasso, "lasso_covs") <- LASSO_COV
  trials_lasso
}

# ============================================================================
# build_daily_input_lasso: 同 v28.4_common.R 的 build_daily_input,
# 但协变量集换成 LASSO_COV
# ============================================================================
build_daily_input_lasso <- function(trials_lasso) {
  LASSO_COV <- attr(trials_lasso, "lasso_covs")
  TAU <- 28L
  dt <- trials_lasso
  # v28.4.1（2026-10-14）：统一使用 v28.4_common.R 的绝对日历日共识死亡日
  dd_agg <- stay_death_day(dt)
  dt_m <- merge(dt, dd_agg, by = "stay_id", all.x = TRUE)
  dt_m <- dt_m[is.na(death_day) | origin_day <= death_day]
  grid <- dt_m[, .(max_origin = max(origin_day),
                   dday = fifelse(any(death_28day_num == 1L),
                                  death_day[1], NA_integer_)),
               by = stay_id]
  grid[, last_day := max_origin + TAU]
  daily <- grid[, .(day = seq(0L, last_day[1])), by = stay_id]
  daily <- merge(daily, grid[, .(stay_id, dday)], by = "stay_id", all.x = TRUE)
  daily[, outcome := 0L]
  daily[!is.na(dday) & day == dday, outcome := 1L]
  daily[, dday := NULL]
  orig <- unique(dt_m[, .(stay_id, origin_day, arm)])
  daily[, eligible := 0L]; daily[, tx_init := 0L]
  daily[orig, on = c(stay_id = "stay_id", day = "origin_day"),
        `:=`(eligible = 1L, tx_init = i.arm)]

  for (v in LASSO_COV) {
    daily[, (v) := na_like(dt_m[[v]])]
    daily[dt_m, on = c(stay_id = "stay_id", day = "origin_day"),
          (v) := get(paste0("i.", v))]
    setorder(daily, stay_id, day)
    daily[, (v) := zoo_locf(get(v)), by = stay_id]
  }
  setorder(daily, stay_id, day)
  # 把 LASSO 协变量挂在 daily 属性, 后续 batch bootstrap 用
  attr(daily, "lasso_covs") <- LASSO_COV
  check_daily_input(daily, label = "lasso")
  daily
}

# ============================================================================
# argv
# ============================================================================
argv <- commandArgs(trailingOnly = TRUE)
IMPS <- as.integer(strsplit(argv[1], ",")[[1]])

# 检查已有 imp, 跳过
existing <- data.frame()
if (file.exists(OUT_CSV)) {
  existing <- read.csv(OUT_CSV)
  done_ids <- existing$imp
} else {
  done_ids <- integer(0)
}
to_run <- setdiff(IMPS, done_ids)
cat(sprintf("=== v28.4 LASSO 31v batch bootstrap ===\n"))
cat(sprintf("待跑 imp: %s (已完成 %s)\n",
            paste(to_run, collapse=","),
            if (length(done_ids)) paste(done_ids, collapse=",") else "none"))

for (imp_id in to_run) {
  cat(sprintf("\n===== imp %d · start =====\n", imp_id))
  t0 <- Sys.time()

  trials_lasso <- load_trials_lasso(imp_id)
  LASSO_COV <- attr(trials_lasso, "lasso_covs")
  cat(sprintf("LASSO covariates = %d vars\n", length(LASSO_COV)))

  daily <- build_daily_input_lasso(trials_lasso)
  rss_mb <- as.numeric(system("ps -o rss= -p $PPID", intern=TRUE))/1024
  cat(sprintf("daily: %d rows, %d stays, RSS=%.0f MB\n",
              nrow(daily), uniqueN(daily$stay_id), rss_mb))

  res <- run_one_imp_batch_bootstrap(
    daily, DAG14 = LASSO_COV,   # ← 关键: DAG14 参数传入 LASSO 协变量集
    n_boot = 200L, batch_size = 50L,
    base_seed = SEED, ncores = NCORES_SAFE, verbose = TRUE)

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  rss_mb  <- as.numeric(system("ps -o rss= -p $PPID", intern=TRUE))/1024

  row <- data.frame(
    imp   = imp_id,
    loghr = res$loghr,
    se    = res$se,
    hr    = res$hr,
    lci   = res$lci,
    uci   = res$uci
  )

  if (!file.exists(OUT_CSV)) {
    write.csv(row, OUT_CSV, row.names = FALSE)
  } else {
    write.table(row, OUT_CSV, append = TRUE, sep = ",",
                row.names = FALSE, col.names = FALSE)
  }

  cat(sprintf("imp %d · done · %.1fs · RSS=%.0f MB · HR=%.4f\n",
              imp_id, elapsed, rss_mb, res$hr))
  rm(trials_lasso, daily); gc(verbose = FALSE)
}

# ---- Rubin pool ----
all <- read.csv(OUT_CSV)
loghr <- all$loghr; se <- all$se; m <- nrow(all)

rubin_pool <- function(lh, se, m) {
  W <- mean(se^2); B <- var(lh); V <- W + B + B/m
  hr_pool <- exp(mean(lh)); se_pool <- sqrt(V)
  z <- qnorm(0.975)
  data.table::data.table(hr=hr_pool, lci=exp(mean(lh)-z*se_pool),
             uci=exp(mean(lh)+z*se_pool),
             p=2*(1-pnorm(abs(mean(lh)/se_pool))), B=B, df=m-1)
}
pooled <- rubin_pool(loghr, se, m)

cat(sprintf("\n===== RUBIN POOLED m=%d · LASSO 31v =====\n", m))
cat(sprintf("HR = %.4f (%.4f - %.4f), p = %.4f\n",
            pooled$hr, pooled$lci, pooled$uci, pooled$p))
cat(sprintf("vs DAG14 主分析 HR=0.972 (0.892-1.058): %.2f%% 差异\n",
            100*abs(pooled$hr - 0.972)/0.972))

# 写 pooled CSV
pool_dt <- data.frame(version="V28.4-LASSO31", date=format(Sys.Date()),
                      cov="lasso31", m=m, n_boot=200, ncores=NCORES_SAFE,
                      hr=pooled$hr, lci=pooled$lci, uci=pooled$uci,
                      p=pooled$p, B=pooled$B, df=pooled$df)
OUT_POOL <- file.path(OUT_DIR, "v28.4_lasso31_pooled_batch50.csv")
write.csv(pool_dt, OUT_POOL, row.names = FALSE)
cat(sprintf("Pooled → %s\n", OUT_POOL))
cat(sprintf("Per-imp CSV → %s\n", OUT_CSV))
