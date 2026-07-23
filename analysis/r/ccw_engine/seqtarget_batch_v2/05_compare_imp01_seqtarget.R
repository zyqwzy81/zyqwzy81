#!/usr/bin/env Rscript
# ============================================================================
# 05_compare_imp01_seqtarget.R
#
# 目的：在首个 MICE 插补集（imp01）上，直接对比 V2 batch engine
#      （SEQexpand + pooled logistic）与原生 SEQTaRget 的 ITT 结果。
#
# 设计原理：
#   SEQTaRget ITT (method="ITT") = 克隆 + SEQexpand + pooled logistic
#                               + followup² + trial² + 协变量
#   - 无 IPTW（ITT 不做人工删失，无需删失权重）
#   - CI 通过 bootstrap 计算（nboot=200）
#
#   V2 batch engine ITT = 相同的 SEQexpand + pooled logistic + cluster SE
#   - 无权重（与 ITT 一致）
#   - CI 通过 cluster-robust sandwich estimator 计算（远快于 bootstrap）
#
# 预期：两方法 HR 应高度一致，V2 的 cluster SE 略小于 bootstrap SE
# ============================================================================

start_time <- Sys.time()

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(sandwich)
  library(lmtest)
  library(jsonlite)
  library(SEQTaRget)     # ← 必须显式加载
})

source("config_v2.R")

# ---- 额外路径（绝对路径，与 v26_9_main_analysis.R 一致） ----
IMP01_PARQUET <- "/data/小米云盘/文档/论文/血磷/sql/result/mice_imputed/v26_mice_imputed_imp01.parquet"
MAIN_PARQUET  <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_imputed_v26.parquet"
OUTPUT_CSV    <- file.path(RESULT_DIR, "v2_vs_seqtarget_imp01_comparison.csv")
OUTPUT_JSON   <- file.path(RESULT_DIR, "v2_vs_seqtarget_imp01_report.json")

# ---- DAG14 协变量集（与 v26_9_main_analysis.R 完全一致） ----
DAG14 <- c(
  "age", "gender", "baseline_sofa", "sepsis3", "charlson",
  "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
  "admission_type", "baseline_phosphate", "lactate_t0", "refeeding_risk",
  "corrected_calcium_mmol_min", "nutrition_support_type"
)

# 预检 imp01 列
IMP01_COLS_AVAIL <- intersect(DAG14, names(as.data.table(arrow::read_parquet(IMP01_PARQUET, as_data_frame = FALSE))))

module_log_section("V2 vs SEQTaRget — 首个 MICE 插补集 (imp01) 对比")
module_log_info(sprintf("SEQTaRget 参考 (imp01): HR=1.0453 (0.9445 - 1.1569), beta=0.0443, SE=0.0517"))

# ---- 1. 数据加载与合并 ----
t_load <- Sys.time()

dt_main <- as.data.table(read_parquet(MAIN_PARQUET))
module_log_info(sprintf("主 Parquet: %d 行, %d 列", nrow(dt_main), ncol(dt_main)))

dt_imp01 <- as.data.table(read_parquet(IMP01_PARQUET))
module_log_info(sprintf("imp01 Parquet: %d 行, %d 列", nrow(dt_imp01), ncol(dt_imp01)))

# DAG14 变量来源：
# - MICE Parquet (imp01) 包含部分 DAG14 协变量（插补后的）
# - 主 Parquet 包含 remaining 部分（如 nutrition_support_type 不在 MICE 中）
dag14_in_imp <- intersect(DAG14, names(dt_imp01))
dag14_from_main <- setdiff(DAG14, dag14_in_imp)

key_from_main <- c("stay_id", "trial_num", "itt_group", "death_28day_num",
                   "surv_time", dag14_from_main)
key_from_main <- intersect(key_from_main, names(dt_main))
dt_key <- dt_main[, ..key_from_main]

from_imp <- c("stay_id", "trial_num", dag14_in_imp)
dt_imp_sub <- dt_imp01[, ..from_imp]

dt_merge <- merge(dt_key, dt_imp_sub,
                  by = c("stay_id", "trial_num"),
                  all.x = TRUE, sort = FALSE)

# 过滤 + 准备 SEQTaRget 必需列
dt_merge[, `:=`(
  id_int       = as.integer(factor(stay_id)),
  period_int   = as.integer(trial_num) - 1L,
  a_int        = as.integer(itt_group == "treatment"),
  y_int        = as.integer(death_28day_num),
  eligible_int = 1L
)]
setorder(dt_merge, id_int, period_int)

# 按 period 过滤
dt_merge[, `:=`(n_a = sum(a_int), n_y = sum(y_int)), by = period_int]
dt_use <- dt_merge[n_a > 0 & n_a < .N & n_y >= 5]

module_log_info(sprintf("合并后: %d 行, %d 期 (过滤后 %d 行)",
                        nrow(dt_merge), uniqueN(dt_merge$period_int), nrow(dt_use)))

# 变量检查与 factor 转换
covs_use <- intersect(DAG14, names(dt_use))
covs_use <- covs_use[sapply(covs_use, function(v) {
  x <- dt_use[[v]]
  if (is.factor(x))   return(nlevels(x) > 1)
  if (is.character(x)) return(length(unique(x)) > 1)
  return(sd(x, na.rm = TRUE) > 1e-8)
})]
for (v in covs_use) {
  if (is.character(dt_use[[v]])) dt_use[, (v) := as.factor(get(v))]
}
module_log_info(sprintf("DAG14 可用协变量: %d/%d", length(covs_use), length(DAG14)))
module_log_info(sprintf("最终数据: %d 行, %d 例患者, %d 事件",
                        nrow(dt_use), uniqueN(dt_use$stay_id), sum(dt_use$y_int)))

rm(dt_main, dt_imp01, dt_key, dt_imp_sub, dt_merge); gc(verbose = FALSE)

load_time <- as.numeric(difftime(Sys.time(), t_load, units = "secs"))

# ============================================================
# 2. V2 Batch: SEQexpand + Pooled Logistic (cluster SE)
# ============================================================
t_v2 <- Sys.time()
module_log_section("V2 Batch Engine (SEQexpand + Pooled Logistic, cluster SE)")

# SEQexpand: 日粒度展开
ds <- copy(dt_use)
ds[, `:=`(
  tx_init    = a_int,
  day_max    = as.integer(floor(pmax(surv_time, 0)))
)]
ds[day_max > SEQEXPAND_MAX_DAY, day_max := SEQEXPAND_MAX_DAY]
ds[day_max < 0, day_max := 0L]

# 展开
expanded <- ds[rep(seq_len(.N), day_max + 1L)]
expanded[, day := seq_len(.N) - 1L, by = .(stay_id, trial_num)]
expanded[, y_daily := 0L]
expanded[day == day_max & y_int == 1L, y_daily := 1L]
expanded[, `:=`(
  followup    = as.numeric(day),
  followup_sq = as.numeric(day^2),
  trial       = as.numeric(trial_num),
  trial_sq    = as.numeric(trial_num^2)
)]

n_expanded <- nrow(expanded)
n_events_daily <- sum(expanded$y_daily)
module_log_info(sprintf("SEQexpand: %d 原始 → %d 日记录 (%.1f 倍)",
                        nrow(ds), n_expanded, n_expanded / nrow(ds)))
module_log_info(sprintf("展开后事件数: %d", n_events_daily))

# Pooled logistic（无权重，= ITT）
base_f <- "y_daily ~ tx_init + followup + followup_sq + trial + trial_sq"
cov_f  <- paste(covs_use, collapse = " + ")

cc_vars <- c("y_daily", "tx_init", "followup", "followup_sq",
             "trial", "trial_sq", "stay_id", covs_use)
cc_idx <- complete.cases(expanded[, ..cc_vars])
exp_cc <- expanded[cc_idx]
module_log_info(sprintf("完整案例: %d/%d (%.1f%%)",
                        nrow(exp_cc), n_expanded, 100 * nrow(exp_cc) / n_expanded))

v2_fit <- glm(as.formula(paste(base_f, cov_f, sep = " + ")),
              data = exp_cc, family = binomial(),
              control = list(maxit = 100))

v2_coef_tx <- coef(v2_fit)["tx_init"]
v2_se_cl   <- sqrt(diag(vcovCL(v2_fit, cluster = exp_cc$stay_id)))["tx_init"]

v2_results <- list(
  method = "V2_batch",
  hr_raw       = exp(v2_coef_tx),
  hr           = exp(v2_coef_tx),
  hr_lower     = exp(v2_coef_tx - 1.96 * v2_se_cl),
  hr_upper     = exp(v2_coef_tx + 1.96 * v2_se_cl),
  beta         = v2_coef_tx,
  se           = v2_se_cl,
  ci_method    = "cluster_sandwich",
  n_expanded   = n_expanded,
  n_events     = n_events_daily,
  run_time_sec = as.numeric(difftime(Sys.time(), t_v2, units = "secs"))
)

module_log_info(sprintf("V2: HR=%.4f (%.4f–%.4f), beta=%.6f, SE=%.6f",
                        v2_results$hr, v2_results$hr_lower, v2_results$hr_upper,
                        v2_results$beta, v2_results$se))

rm(ds, expanded, exp_cc, v2_fit); gc(verbose = FALSE)

# ============================================================
# 3. 原生 SEQTaRget ITT (nboot=200)
# ============================================================
t_seq <- Sys.time()
module_log_section("原生 SEQTaRget ITT (DAG14, nboot=200)")

# 手动转换为干净 data.frame（避免 data.table 属性问题）
dt_seq <- as.data.frame(dt_use)
module_log_info(sprintf("SEQTaRget 输入: data.frame %d × %d", nrow(dt_seq), ncol(dt_seq)))

seq_result <- tryCatch({
  SEQuential(
    data  = dt_seq,
    id.col       = "id_int",
    time.col     = "period_int",
    eligible.col = "eligible_int",
    treatment.col = "a_int",
    outcome.col   = "y_int",
    fixed.cols          = covs_use,
    time_varying.cols   = character(0),
    method = "ITT",
    options = SEQopts(
      bootstrap = TRUE,
      bootstrap.nboot = 200L,
      bootstrap.CI    = 0.95,
      ncores = 1L,
      fastglm.method = 3L,
      hazard = TRUE,
      trial.include    = TRUE,
      followup.include = TRUE,
      data.return      = FALSE
    ),
    verbose = TRUE
  )
}, error = function(e) {
  module_log_warn(sprintf("SEQTaRget 异常: %s", e$message))
  # 尝试无 bootstrap（仅点估计）
  tryCatch({
    module_log_info("  回退: 无 bootstrap 重试...")
    SEQuential(
      data  = dt_seq,
      id.col       = "id_int",
      time.col     = "period_int",
      eligible.col = "eligible_int",
      treatment.col = "a_int",
      outcome.col   = "y_int",
      fixed.cols          = covs_use,
      time_varying.cols   = character(0),
      method = "ITT",
      options = SEQopts(
        hazard = TRUE,
        bootstrap = FALSE,
        trial.include    = TRUE,
        followup.include = TRUE,
        data.return      = FALSE,
        fastglm.method = 3L
      ),
      verbose = TRUE
    )
  }, error = function(e2) {
    module_log_warn(sprintf("SEQTaRget 回退也失败: %s", e2$message))
    NULL
  })
})

if (is.null(seq_result)) {
  seq_results <- list(
    method = "SEQTaRget_ITT",
    hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
    beta = NA_real_, se = NA_real_, converged = FALSE,
    run_time_sec = as.numeric(difftime(Sys.time(), t_seq, units = "secs"))
  )
} else {
  hr_vec <- hazard_ratio(seq_result)
  seq_hr    <- hr_vec[1]
  seq_lower <- hr_vec[2]
  seq_upper <- hr_vec[3]
  seq_beta  <- log(seq_hr)
  seq_se    <- (log(seq_upper) - log(seq_lower)) / (2 * 1.96)

  seq_results <- list(
    method = "SEQTaRget_ITT",
    hr_raw       = seq_hr,
    hr           = seq_hr,
    hr_lower     = seq_lower,
    hr_upper     = seq_upper,
    beta         = seq_beta,
    se           = seq_se,
    ci_method    = "bootstrap_200",
    n_expanded   = NA_integer_,
    n_events     = sum(dt_use$y_int),
    converged    = TRUE,
    run_time_sec = as.numeric(difftime(Sys.time(), t_seq, units = "secs"))
  )

  module_log_info(sprintf("SEQTaRget: HR=%.4f (%.4f–%.4f), beta=%.6f, SE=%.6f",
                          seq_hr, seq_lower, seq_upper, seq_beta, seq_se))
}

rm(seq_result); gc(verbose = FALSE)

# ============================================================
# 4. 对比分析
# ============================================================
module_log_section("V2 vs SEQTaRget 对比")

# 4a: 差异计算
hr_diff <- abs(v2_results$hr - seq_results$hr)
hr_ratio <- v2_results$hr / seq_results$hr
beta_diff <- abs(v2_results$beta - seq_results$beta)
se_ratio <- v2_results$se / seq_results$se

comparison <- data.table(
  metric = c("HR", "HR_lower", "HR_upper", "beta(log_HR)", "SE", "run_time_min"),
  V2_batch = c(v2_results$hr, v2_results$hr_lower, v2_results$hr_upper,
               v2_results$beta, v2_results$se,
               round(v2_results$run_time_sec / 60, 2)),
  SEQTaRget = c(seq_results$hr, seq_results$hr_lower, seq_results$hr_upper,
                seq_results$beta, seq_results$se,
                round(seq_results$run_time_sec / 60, 2))
)

# 添加差异列
comparison[, diff := abs(V2_batch - SEQTaRget)]
comparison[, pct_diff := round(diff / abs(SEQTaRget) * 100, 2)]

# 打印
for (i in seq_len(nrow(comparison))) {
  cat(sprintf("  %-15s: V2=%.6f, SEQ=%.6f, |差异|=%.6f (%.2f%%)\n",
              comparison$metric[i],
              comparison$V2_batch[i],
              comparison$SEQTaRget[i],
              comparison$diff[i],
              comparison$pct_diff[i]))
}

# 4b: 收敛判定
convergence <- list(
  hr_relative_diff_pct = abs(v2_results$hr - seq_results$hr) / seq_results$hr * 100,
  beta_absolute_diff   = abs(v2_results$beta - seq_results$beta),
  se_ratio             = se_ratio,
  ci_overlap           = TRUE,  # 若两CI重叠则为TRUE
  hr_direction_consistent = sign(log(v2_results$hr)) == sign(log(seq_results$hr)),
  runtime_speedup      = round(seq_results$run_time_sec / v2_results$run_time_sec, 1)
)

convergence$ci_overlap <- !(
  v2_results$hr_upper < seq_results$hr_lower ||
  seq_results$hr_upper < v2_results$hr_lower
)

module_log_section("收敛判定")
module_log_info(sprintf("HR 相对差异: %.4f%%", convergence$hr_relative_diff_pct))
module_log_info(sprintf("Beta 绝对差异: %.6f", convergence$beta_absolute_diff))
module_log_info(sprintf("SE 比率 (V2/SEQ): %.4f", convergence$se_ratio))
module_log_info(sprintf("CI 重叠: %s", convergence$ci_overlap))
module_log_info(sprintf("方向一致: %s", convergence$hr_direction_consistent))
module_log_info(sprintf("V2 速度提升: %.1f 倍", convergence$runtime_speedup))

# ============================================================
# 5. 保存输出
# ============================================================
module_log_section("保存输出")

# CSV
fwrite(comparison, OUTPUT_CSV)
module_log_info(sprintf("对比表: %s", OUTPUT_CSV))

# JSON 报告
report <- list(
  title = "V2 Batch Engine vs Native SEQTaRget — imp01 DAG14 ITT 对比",
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  data = list(
    source_parquet = MAIN_PARQUET,
    imp01_parquet = IMP01_PARQUET,
    n_patients = uniqueN(dt_use$stay_id),
    n_records_pre_expand = nrow(dt_use),
    n_records_expanded = n_expanded,
    n_events = n_events_daily,
    dag14_used = covs_use,
    n_covariates = length(covs_use)
  ),
  v2 = v2_results[c("method", "hr", "hr_lower", "hr_upper", "beta", "se",
                    "ci_method", "run_time_sec")],
  seqtarget = seq_results[c("method", "hr", "hr_lower", "hr_upper", "beta", "se",
                           "ci_method", "run_time_sec")],
  comparison = list(
    hr_relative_diff_pct = round(convergence$hr_relative_diff_pct, 4),
    beta_absolute_diff   = round(convergence$beta_absolute_diff, 6),
    se_ratio             = round(convergence$se_ratio, 4),
    ci_overlap           = convergence$ci_overlap,
    direction_consistent = convergence$hr_direction_consistent,
    runtime_speedup      = convergence$runtime_speedup,
    v2_less_than_5pct    = convergence$hr_relative_diff_pct < 5
  ),
  runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2),
  peak_mem_mb = get_peak_mem_mb()
)

writeLines(toJSON(report, pretty = TRUE, auto_unbox = TRUE), OUTPUT_JSON)
module_log_info(sprintf("JSON 报告: %s", OUTPUT_JSON))

# ============================================================
# 6. 摘要
# ============================================================
module_log_section("对比完成")
cat("\n")
cat("==================================================================\n")
cat("  V2 Batch vs Native SEQTaRget 对比 — 首个 MICE 插补集 (imp01)\n")
cat("==================================================================\n\n")

cat(sprintf("  V2 (SEQexpand + pooled logistic, cluster SE):\n"))
cat(sprintf("    HR = %.4f (%.4f - %.4f)\n",
            v2_results$hr, v2_results$hr_lower, v2_results$hr_upper))
cat(sprintf("    beta = %.6f, SE = %.6f\n", v2_results$beta, v2_results$se))
cat(sprintf("    用时: %.1f 分钟\n\n", v2_results$run_time_sec / 60))

cat(sprintf("  SEQTaRget (native ITT, boot=200):\n"))
cat(sprintf("    HR = %.4f (%.4f - %.4f)\n",
            seq_results$hr, seq_results$hr_lower, seq_results$hr_upper))
cat(sprintf("    beta = %.6f, SE = %.6f\n", seq_results$beta, seq_results$se))
cat(sprintf("    用时: %.1f 分钟\n\n", seq_results$run_time_sec / 60))

cat(sprintf("  HR 相对差异: %s %s\n",
            if (is.na(convergence$hr_relative_diff_pct)) "NA" else sprintf("%.4f%%", convergence$hr_relative_diff_pct),
            if (isTRUE(convergence$hr_relative_diff_pct < 5)) "✓" else if (is.na(convergence$hr_relative_diff_pct)) "?" else "✗"))
cat(sprintf("  方向一致: %s\n",
            if (isTRUE(convergence$hr_direction_consistent)) "✓" else if (isTRUE(!convergence$hr_direction_consistent)) "✗" else "?"))
cat(sprintf("  CI 重叠: %s\n",
            if (isTRUE(convergence$ci_overlap)) "✓" else if (isTRUE(!convergence$ci_overlap)) "✗" else "?"))
cat(sprintf("  V2 速度提升: %s 倍\n\n",
            if (is.finite(convergence$runtime_speedup)) sprintf("%.1f", convergence$runtime_speedup) else "NA"))

cat(sprintf("  总运行时间: %.1f 分钟\n",
            as.numeric(difftime(Sys.time(), start_time, units = "secs")) / 60))
cat("==================================================================\n")
