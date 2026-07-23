# ============================================================================
# V26.9 Bootstrap 比较分析 — nboot=100 vs nboot=200
#
# 目的：
#   在第 1 个 MICE 插补集上，分别以 nboot=100 和 nboot=200
#   运行 SEQTaRget ITT（DAG14），系统比较：
#     1. HR、95% CI、p 值等结果指标的差异
#     2. 运行时间、内存占用等计算效率指标
#   为后续分析选择最优 bootstrap 数量提供依据。
#
# 注意：nboot=200 结果已存在于 v26_9_seqtarget_per_imp_results.csv（imp01），
# 此处仅运行 nboot=100，然后与已有数据进行对比。
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(survival)
  library(SEQTaRget)
})

# ============================================================================
# 1. 参数配置
# ============================================================================

# 输出目录（与 V26.9 主分析共用）
OUT_DIR   <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9"
TBL_DIR   <- file.path(OUT_DIR, "tables")
LOG_DIR   <- file.path(OUT_DIR, "logs")
for (d in c(OUT_DIR, TBL_DIR, LOG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# 数据路径
DATA_FILE            <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_imputed_v26.parquet"
MICE_PARQUET_IMP01   <- "/data/小米云盘/文档/论文/血磷/sql/result/mice_imputed/v26_mice_imputed_imp01.parquet"

# DAG14 协变量集
DAG14 <- c(
  "age", "gender", "baseline_sofa", "sepsis3", "charlson",
  "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
  "admission_type", "baseline_phosphate", "lactate_t0", "refeeding_risk",
  "corrected_calcium_mmol_min", "nutrition_support_type"
)

# nboot=100 配置
BOOT_TEST <- 100L
SEED      <- 20260629L
set.seed(SEED)
options(mc.cores = 1)

# ============================================================================
# 2. 日志系统
# ============================================================================
log_file <- file.path(LOG_DIR, sprintf("v26_9_bootstrap_comp_%s.log",
                                        format(Sys.time(), "%Y%m%d_%H%M%S")))
log_con  <- file(log_file, open = "wt")

log_info <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg),
      file = log_con, append = TRUE)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
  flush(stdout())
}
log_section <- function(title) {
  cat(sprintf("\n========== %s ==========\n", title),
      file = log_con, append = TRUE)
  cat(sprintf("\n========== %s ==========\n", title))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

log_info("V26.9 Bootstrap 比较分析启动")
log_info(sprintf("种子: %d, nboot=%d, ncores=1", SEED, BOOT_TEST))
log_info(sprintf("输出目录: %s", OUT_DIR))

# ============================================================================
# 3. 数据加载（第 1 插补集）
# ============================================================================
log_section("数据加载 (imp01)")

# ---- 3a. 加载主数据 key 列 ----
dt_full <- as.data.table(read_parquet(DATA_FILE))
dt_full[, `:=`(
  tx_int   = as.integer(itt_group == "treatment"),
  event    = death_28day_num,
  id_int   = as.integer(factor(stay_id))
)]
setorder(dt_full, id_int, trial_num)

key_cols <- c("stay_id", "subject_id", "trial_num", "time_zero",
              "itt_group", "gender", "admission_type",
              "surv_time", "death_time_from_t0", "dischtime_from_t0",
              "race_group", "icu_type",
              "death_28day_num", "nutrition_support_type")
key_cols <- intersect(key_cols, names(dt_full))
dt_key <- dt_full[, ..key_cols]
log_info(sprintf("主数据: %d 行, %d 事件, key cols=%d",
                 nrow(dt_full), sum(dt_full$event), length(key_cols)))
rm(dt_full); gc(TRUE)

# ---- 3b. 加载 MICE imp01 Parquet ----
dt_imp <- as.data.table(read_parquet(MICE_PARQUET_IMP01))
keep_cols <- intersect(c("stay_id", "trial_num", DAG14), names(dt_imp))
dt_imp <- dt_imp[, ..keep_cols]

# DAG14 在 MICE 中的列
dag14_in_imp <- intersect(DAG14, names(dt_imp))
dag14_not_imp <- setdiff(DAG14, names(dt_imp))
log_info(sprintf("MICE imp01: %d 行, %d 列", nrow(dt_imp), ncol(dt_imp)))
log_info(sprintf("DAG14 在 MICE 中: %s", paste(dag14_in_imp, collapse = ", ")))
log_info(sprintf("DAG14 从主数据获取: %s", paste(dag14_not_imp, collapse = ", ")))

# ---- 3c. Merge 对齐 ----
dag14_remove_from_key <- intersect(dag14_in_imp, names(dt_key))
by_cols <- c("stay_id", "trial_num")

dt_comb <- merge(
  dt_key[, setdiff(names(dt_key), dag14_remove_from_key), with = FALSE],
  dt_imp,
  by = by_cols, all.x = TRUE, sort = FALSE
)
log_info(sprintf("Merge 后: %d 行, %d 列", nrow(dt_comb), ncol(dt_comb)))
rm(dt_key, dt_imp); gc(TRUE)

# ---- 3d. 字符型变量转 factor ----
char_covs <- intersect(DAG14, names(dt_comb))
char_covs <- char_covs[sapply(char_covs, function(v) is.character(dt_comb[[v]]))]
for (v in char_covs) {
  dt_comb[[v]] <- as.factor(dt_comb[[v]])
}
if (length(char_covs) > 0) {
  log_info(sprintf("字符型→factor: %s", paste(char_covs, collapse = ", ")))
}

# ============================================================================
# 4. SEQTaRget ITT 运行函数（含性能监控）
# ============================================================================
log_section("SEQTaRget ITT 运行")

run_seqtarget_imp01 <- function(dt, dag14_covs, boot) {
  dt_sub <- copy(dt)
  dt_sub[, `:=`(
    id_int       = as.integer(factor(stay_id)),
    period_int   = as.integer(trial_num) - 1L,
    a_int        = as.integer(itt_group == "treatment"),
    y_int        = as.integer(death_28day_num),
    event        = as.integer(death_28day_num),
    eligible_int = 1L
  )]
  setorder(dt_sub, id_int, period_int)

  # 按 period 过滤
  dt_sub[, `:=`(n_a = sum(a_int), n_y = sum(y_int)), by = period_int]
  dt_sub <- dt_sub[n_a > 0 & n_a < .N & n_y >= 5]
  n_before <- nrow(dt_sub)

  if (nrow(dt_sub) < 100) {
    return(list(error = "样本量不足"))
  }

  # 协变量检查
  covs <- intersect(dag14_covs, names(dt_sub))
  covs <- covs[sapply(covs, function(v) {
    x <- dt_sub[[v]]
    if (is.factor(x))    return(nlevels(x) > 1)
    if (is.character(x)) return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]
  if (length(covs) < 2) {
    return(list(error = "协变量不足"))
  }

  # 字符型转 factor
  for (v in covs) {
    if (is.character(dt_sub[[v]])) dt_sub[, (v) := as.factor(get(v))]
  }

  t0 <- Sys.time()
  result <- tryCatch({
    SEQuential(
      data            = as.data.frame(dt_sub),
      id.col          = "id_int",
      time.col        = "period_int",
      eligible.col    = "eligible_int",
      treatment.col   = "a_int",
      outcome.col     = "y_int",
      fixed.cols      = covs,
      time_varying.cols = character(0),
      method          = "ITT",
      options = SEQopts(
        bootstrap       = TRUE,
        bootstrap.nboot = boot,
        bootstrap.CI    = 0.95,
        ncores          = 1L,
        fastglm.method  = 3L,
        hazard          = TRUE,
        trial.include    = TRUE,
        followup.include = TRUE,
        data.return      = FALSE
      ),
      verbose = FALSE
    )
  }, error = function(e) {
    list(error = e$message)
  })
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (is.list(result) && !is.null(result$error)) {
    return(list(error = result$error, runtime = runtime))
  }

  hr_vec <- tryCatch(hazard_ratio(result), error = function(e) c(NA, NA, NA))
  rm(result); gc(TRUE)

  return(list(
    hr        = hr_vec[1],
    ci_lower  = hr_vec[2],
    ci_upper  = hr_vec[3],
    runtime   = runtime,
    n         = nrow(dt_sub)
  ))
}

# ============================================================================
# 5. 运行 nboot=100
# ============================================================================
log_info(sprintf("\n>>>> 启动 nboot=%d <<<<", BOOT_TEST))

r_100 <- run_seqtarget_imp01(dt_comb, DAG14, BOOT_TEST)

if (!is.null(r_100$error)) {
  log_info(sprintf("  [nboot=%d] 错误: %s", BOOT_TEST, r_100$error))
} else {
  log_info(sprintf("  [nboot=%d] HR=%.4f (%.4f–%.4f), n=%d, runtime=%.1f min (%.0f sec)",
                   BOOT_TEST, r_100$hr, r_100$ci_lower, r_100$ci_upper,
                   r_100$n, r_100$runtime / 60, r_100$runtime))
}

# ============================================================================
# 6. 载入 nboot=200（从已有结果）
# ============================================================================
log_section("载入 nboot=200 参考结果")

imp01_boot200 <- list(
  hr        = 1.04531132883746,
  ci_lower  = 0.944486975558084,
  ci_upper  = 1.15689872118171,
  runtime   = 4838.05685424805,  # 秒
  n         = 32325
)
imp01_boot200$log_hr     <- log(imp01_boot200$hr)
imp01_boot200$se_log     <- (log(imp01_boot200$ci_upper) - log(imp01_boot200$ci_lower)) / (2 * 1.96)
imp01_boot200$ci_width   <- imp01_boot200$ci_upper - imp01_boot200$ci_lower

log_info(sprintf("  [nboot=200] HR=%.4f (%.4f–%.4f), n=%d, runtime=%.1f min (%.0f sec)",
                 imp01_boot200$hr, imp01_boot200$ci_lower, imp01_boot200$ci_upper,
                 imp01_boot200$n, imp01_boot200$runtime / 60, imp01_boot200$runtime))
log_info(sprintf("  [nboot=200] log(HR)=%.6f, SE(log)=%.6f, CI宽度=%.4f",
                 imp01_boot200$log_hr, imp01_boot200$se_log, imp01_boot200$ci_width))

# ============================================================================
# 7. 比较分析
# ============================================================================
log_section("Bootstrap 比较分析")

if (is.null(r_100$error)) {
  r_100$log_hr   <- log(r_100$hr)
  r_100$se_log   <- (log(r_100$ci_upper) - log(r_100$ci_lower)) / (2 * 1.96)
  r_100$ci_width <- r_100$ci_upper - r_100$ci_lower

  # --- 7a. 点估计差异 ---
  hr_diff     <- abs(r_100$hr - imp01_boot200$hr)
  hr_diff_pct <- hr_diff / imp01_boot200$hr * 100

  # --- 7b. 置信区间宽度差异 ---
  width_diff     <- abs(r_100$ci_width - imp01_boot200$ci_width)
  width_diff_pct <- width_diff / imp01_boot200$ci_width * 100

  # --- 7c. 对数尺度 SE 差异 ---
  se_diff     <- abs(r_100$se_log - imp01_boot200$se_log)
  se_diff_pct <- se_diff / imp01_boot200$se_log * 100

  # --- 7d. 运行时间 ---
  runtime_ratio <- r_100$runtime / imp01_boot200$runtime

  # --- 7e. 效率比 ---
  # 效率 = 1 / (SE^2 * runtime)
  eff_100 <- 1 / (r_100$se_log^2 * r_100$runtime)
  eff_200 <- 1 / (imp01_boot200$se_log^2 * imp01_boot200$runtime)
  eff_ratio <- eff_100 / eff_200

  # --- 输出比较表 ---
  cat(c(
    "\n",
    "============================================================================\n",
    "  Bootstrap 参数比较报告 (imp01, DAG14 ITT)\n",
    "============================================================================\n",
    "\n",
    sprintf("  %-30s %12s %12s %12s\n", "指标", "nboot=100", "nboot=200", "差异"),
    sprintf("  %-30s %12s %12s %12s\n", "---", "---", "---", "---"),
    sprintf("  %-30s %12.4f %12.4f %10.2f%%\n", "HR", r_100$hr, imp01_boot200$hr, hr_diff_pct),
    sprintf("  %-30s %12.4f %12.4f %10.2f%%\n", "CI 下限", r_100$ci_lower, imp01_boot200$ci_lower, 
            abs(r_100$ci_lower - imp01_boot200$ci_lower) / imp01_boot200$ci_lower * 100),
    sprintf("  %-30s %12.4f %12.4f %10.2f%%\n", "CI 上限", r_100$ci_upper, imp01_boot200$ci_upper,
            abs(r_100$ci_upper - imp01_boot200$ci_upper) / imp01_boot200$ci_upper * 100),
    sprintf("  %-30s %12.4f %12.4f %10.2f%%\n", "CI 宽度", r_100$ci_width, imp01_boot200$ci_width, width_diff_pct),
    sprintf("  %-30s %12.6f %12.6f %10.2f%%\n", "log(HR)", r_100$log_hr, imp01_boot200$log_hr,
            abs(r_100$log_hr - imp01_boot200$log_hr) / abs(imp01_boot200$log_hr) * 100),
    sprintf("  %-30s %12.6f %12.6f %10.2f%%\n", "SE(log(HR))", r_100$se_log, imp01_boot200$se_log, se_diff_pct),
    sprintf("  %-30s %12.0f %12.0f %10.2fx\n", "运行时间 (秒)", r_100$runtime, imp01_boot200$runtime, runtime_ratio),
    sprintf("  %-30s %12.1f %12.1f %10s\n", "运行时间 (分钟)", r_100$runtime / 60, imp01_boot200$runtime / 60, ""),
    "\n",
    "  --- 效率分析 ---\n",
    "\n",
    sprintf("  %-30s %12s %12s\n", "效率指标", "nboot=100", "nboot=200"),
    sprintf("  %-30s %12.4f %12.4f\n", "效率 (1/(SE^2×runtime))", eff_100, eff_200),
    sprintf("  %-30s %12.4f (nboot=100/nboot=200)\n", "相对效率", eff_ratio),
    sprintf("  %-30s %12.4f\n", "SE 比率 (100/200)", se_diff / imp01_boot200$se_log + 1),
    sprintf("  %-30s %12.2f\n", "时间比率 (100/200)", runtime_ratio),
    "\n"
  ), file = log_con, append = TRUE)

  # 直接输出到控制台
  cat(sprintf("\n%s\n", paste(rep("=", 70), collapse = "")))
  cat("  Bootstrap 参数比较报告 (imp01, DAG14 ITT)\n")
  cat(sprintf("%s\n\n", paste(rep("=", 70), collapse = "")))
  
  cat(sprintf("  %-25s %12s %12s %10s\n", "指标", "nboot=100", "nboot=200", "差异"))
  cat(sprintf("  %s\n", paste(rep("-", 60), collapse = "")))
  cat(sprintf("  %-25s %12.4f %12.4f %8.2f%%\n", "HR", r_100$hr, imp01_boot200$hr, hr_diff_pct))
  cat(sprintf("  %-25s %12.4f %12.4f %8.2f%%\n", "CI 下限", r_100$ci_lower, imp01_boot200$ci_lower,
              abs(r_100$ci_lower - imp01_boot200$ci_lower) / imp01_boot200$ci_lower * 100))
  cat(sprintf("  %-25s %12.4f %12.4f %8.2f%%\n", "CI 上限", r_100$ci_upper, imp01_boot200$ci_upper,
              abs(r_100$ci_upper - imp01_boot200$ci_upper) / imp01_boot200$ci_upper * 100))
  cat(sprintf("  %-25s %12.4f %12.4f %8.2f%%\n", "CI 宽度", r_100$ci_width, imp01_boot200$ci_width, width_diff_pct))
  cat(sprintf("  %-25s %12.6f %12.6f %8.2f%%\n", "SE(log(HR))", r_100$se_log, imp01_boot200$se_log, se_diff_pct))
  cat(sprintf("  %-25s %12.0f 秒 %12.0f 秒  %6.2fx\n", "运行时间", r_100$runtime, imp01_boot200$runtime, runtime_ratio))
  cat(sprintf("  %-25s %12.1f 分 %12.1f 分\n", "", r_100$runtime / 60, imp01_boot200$runtime / 60))
  
  cat(sprintf("\n  --- 效率分析 ---\n"))
  cat(sprintf("  %-25s %12.4f %12.4f\n", "效率(1/(SE²×T))", eff_100, eff_200))
  cat(sprintf("  %-25s %12.4f\n", "相对效率(100/200)", eff_ratio))
  
  # --- 7f. 结论 ---
  cat(sprintf("\n  --- 综合评估 ---\n"))
  if (hr_diff_pct < 1 && width_diff_pct < 5 && se_diff_pct < 5) {
    cat(sprintf("  ✓ nboot=100 与 nboot=200 结果高度一致\n"))
    cat(sprintf("  ✓ 点估计差异 < 1%%, CI 宽度差异 < 5%%\n"))
    cat(sprintf("  ✓ 时间节省: %.1f× 加速\n", 1/runtime_ratio))
    cat(sprintf("  → 建议: 对探索性分析可使用 nboot=100，主分析建议 nboot=200\n"))
  } else if (hr_diff_pct < 2 && width_diff_pct < 10) {
    cat(sprintf("  △ 结果基本一致，但部分指标差异中等\n"))
    cat(sprintf("  → 建议: 主分析应使用 nboot=200, 敏感性分析可用 nboot=100\n"))
  } else {
    cat(sprintf("  ✗ 两种设置存在实质性差异\n"))
    cat(sprintf("  → 建议: 主分析必须使用 nboot=200\n"))
  }

  # --- 保存比较结果 CSV ---
  dt_comp <- data.table(
    metric   = c("HR", "CI_lower", "CI_upper", "CI_width", "log_HR", "SE_log", 
                 "runtime_sec", "runtime_min", "n"),
    boot_100 = c(r_100$hr, r_100$ci_lower, r_100$ci_upper, r_100$ci_width, 
                 r_100$log_hr, r_100$se_log, r_100$runtime, r_100$runtime/60, r_100$n),
    boot_200 = c(imp01_boot200$hr, imp01_boot200$ci_lower, imp01_boot200$ci_upper, 
                 imp01_boot200$ci_width, imp01_boot200$log_hr, imp01_boot200$se_log,
                 imp01_boot200$runtime, imp01_boot200$runtime/60, imp01_boot200$n),
    diff_pct = c(hr_diff_pct, 
                 abs(r_100$ci_lower - imp01_boot200$ci_lower) / abs(imp01_boot200$ci_lower) * 100,
                 abs(r_100$ci_upper - imp01_boot200$ci_upper) / abs(imp01_boot200$ci_upper) * 100,
                 width_diff_pct,
                 abs(r_100$log_hr - imp01_boot200$log_hr) / abs(imp01_boot200$log_hr) * 100,
                 se_diff_pct,
                 abs(r_100$runtime - imp01_boot200$runtime) / imp01_boot200$runtime * 100,
                 abs(r_100$runtime - imp01_boot200$runtime) / imp01_boot200$runtime * 100,
                 0)
  )
  fwrite(dt_comp, file.path(TBL_DIR, "v26_9_bootstrap_comparison.csv"))
  log_info(sprintf("比较结果已保存: %s",
                   file.path(TBL_DIR, "v26_9_bootstrap_comparison.csv")))
}

# ============================================================================
# 8. 清理与收尾
# ============================================================================
close(log_con)
cat(sprintf("\n日志文件: %s\n", log_file))
cat("Bootstrap 比较分析完成 ✓\n")
