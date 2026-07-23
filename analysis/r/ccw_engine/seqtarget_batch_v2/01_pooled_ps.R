# ============================================================================
# 01_pooled_ps.R — Phase 1: Pooled PS + IPTW
# 
# 核心改进: 全量数据拟合一个 pooled PS 模型（而非 V1 的 per-period 模型）
# 与 SEQTaRget weight.preexpansion=TRUE 对齐:
#   - 所有 period 合并，在原始 patient-period 数据上拟合 PS
#   - 计算稳定化 IPTW
#   - 将权重映射回各 period（为后续日粒度展开做准备）
#
# 输出:
#   - pooled_ps_model.rds: PS 模型对象
#   - pre_ps_weights.csv: 全量权重摘要
#   - expanded/weights_by_period/: 按 period 分区的权重 Parquet
# ============================================================================

source("config_v2.R")

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(jsonlite)
})

module_log_section("Phase 1: Pooled PS + IPTW（全量数据一次性拟合）")

t_start_total <- Sys.time()

# ---- 1. 读取全量数据 ----
dt <- as.data.table(arrow::read_parquet(INPUT_PARQUET))
set.seed(RANDOM_SEED)
setorder(dt, stay_id, trial_num)

module_log_info(sprintf("原始数据: %d 行, %d 列, %d 例患者",
                        nrow(dt), ncol(dt), length(unique(dt[[ID_COL]]))))

# ---- 2. 数据预处理 ----
# 2a: 治疗变量编码
dt[, tx_init := as.integer(get(TREATMENT_COL) == TREAT_LABEL)]
dt[, outcome := as.integer(get(OUTCOME_COL))]

# 2b: PS 协变量预处理
prepare_covariates <- function(dt, covs) {
  covs_avail <- intersect(covs, names(dt))
  keep <- c()
  for (col in covs_avail) {
    x <- dt[[col]]
    if (is.character(x)) {
      x <- as.factor(x)
      data.table::set(dt, j = col, value = x)
    }
    if (is.numeric(x)) {
      if (sd(x, na.rm = TRUE) > 1e-12) keep <- c(keep, col)
    } else if (is.factor(x)) {
      if (length(levels(x)) >= 2) keep <- c(keep, col)
    } else {
      keep <- c(keep, col)
    }
  }
  keep
}

covs_use <- prepare_covariates(dt, PS_COVARIATES)
dropped <- setdiff(PS_COVARIATES, covs_use)
if (length(dropped) > 0) {
  module_log_warn(sprintf("剔除零方差协变量: %s", paste(dropped, collapse = ", ")))
}
module_log_info(sprintf("PS 模型使用 %d 个协变量", length(covs_use)))

# 2c: 过滤无效期（事件数/样本量不足的 period）
dt[, period_events := sum(outcome, na.rm = TRUE), by = trial_num]
dt[, period_n := .N, by = trial_num]
dt[, period_tx_var := length(unique(tx_init)) > 1, by = trial_num]
dt_valid <- dt[period_events >= PERIOD_FILTER_MIN_EVENTS &
               period_n >= PERIOD_FILTER_MIN_N &
               period_tx_var == TRUE]

valid_periods <- sort(unique(dt_valid$trial_num))
module_log_info(sprintf("有效 period: %d 个 (%s)", length(valid_periods),
                         paste(valid_periods, collapse = ", ")))
dropped_periods <- setdiff(unique(dt$trial_num), valid_periods)
if (length(dropped_periods) > 0) {
  module_log_warn(sprintf("排除的 period（事件/样本不足或无治疗变异）: %s",
                          paste(dropped_periods, collapse = ", ")))
}

# ---- 3. Pooled PS 拟合 ----
# 使用完整案例（complete cases），与 SEQTaRget 默认一致
cc_idx <- complete.cases(dt_valid[, ..covs_use])
dt_cc <- dt_valid[cc_idx]
module_log_info(sprintf("完整案例: %d/%d (%.1f%%)",
                        nrow(dt_cc), nrow(dt_valid),
                        100 * nrow(dt_cc) / nrow(dt_valid)))

formula_str <- paste("tx_init ~", paste(covs_use, collapse = " + "))
t_start_ps <- Sys.time()

ps_model <- tryCatch(
  glm(as.formula(formula_str), data = dt_cc, family = binomial()),
  error = function(e) stop(sprintf("Pooled PS 拟合失败: %s", e$message))
)

ps_runtime <- as.numeric(difftime(Sys.time(), t_start_ps, units = "secs"))
module_log_info(sprintf("Pooled PS 拟合完成: 时间=%.2fs, AIC=%.1f",
                        ps_runtime, AIC(ps_model)))

# ---- 4. 计算 propensity scores ----
ps <- predict(ps_model, type = "response")
ps <- pmax(pmin(ps, PS_CLIP[2]), PS_CLIP[1])

# ---- 5. 稳定化 IPTW ----
p_treat <- mean(dt_cc$tx_init, na.rm = TRUE)
iptw <- ifelse(dt_cc$tx_init == TX_LEVEL,
               p_treat / ps,
               (1 - p_treat) / (1 - ps))

# ---- 6. 双侧截断 ----
q_upper <- quantile(iptw, IPTW_TRUNC_UPPER, na.rm = TRUE)
q_lower <- quantile(iptw, IPTW_TRUNC_LOWER, na.rm = TRUE)
iptw_trunc <- pmax(pmin(iptw, q_upper), q_lower)

dt_cc[, `:=`(ps = ps, iptw = iptw, iptw_trunc = iptw_trunc)]
module_log_info(sprintf("IPTW: 均值=%.4f, SD=%.4f, 截断后[%.4f, %.4f]",
                        mean(iptw_trunc), sd(iptw_trunc), min(iptw_trunc), max(iptw_trunc)))

# ESS
ess <- sum(iptw_trunc)^2 / sum(iptw_trunc^2)
module_log_info(sprintf("有效样本量(ESS): %.1f / %d (%.1f%%)",
                        ess, length(iptw_trunc), 100 * ess / length(iptw_trunc)))

# ---- 7. 按 period 分区保存权重 ----
weight_summaries <- list()
for (p in valid_periods) {
  sub <- dt_cc[trial_num == p]
  out_path <- file.path(WEIGHT_DIR, sprintf("trial_period=%d", p))
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  
  # 保存核心字段：stay_id + 权重 + PS
  weight_cols <- c("stay_id", "trial_num", "trial_period" = "trial_num",
                   "tx_init", "outcome", "surv_time",
                   "ps", "iptw", "iptw_trunc",
                   covs_use)
  # 重命名 trial_num → trial_period（与展开步骤对齐）
  sub_out <- data.table::copy(sub)
  data.table::setnames(sub_out, "trial_num", "trial_period")
  
  arrow::write_parquet(
    sub_out[, .(stay_id, trial_period, tx_init, outcome, surv_time,
                ps, iptw, iptw_trunc, ..covs_use)],
    file.path(out_path, "part-0.parquet")
  )
  
  summ <- data.table(
    trial_period = p,
    n_patients = nrow(sub),
    n_events = sum(sub$outcome, na.rm = TRUE),
    treatment_rate = round(mean(sub$tx_init), 4),
    ps_mean = round(mean(sub$ps), 4),
    ps_sd = round(sd(sub$ps), 4),
    iptw_trunc_mean = round(mean(sub$iptw_trunc), 4),
    iptw_trunc_sd = round(sd(sub$iptw_trunc), 4),
    ess = round(sum(sub$iptw_trunc)^2 / sum(sub$iptw_trunc^2), 1)
  )
  weight_summaries[[as.character(p)]] <- summ
  rm(sub, sub_out); gc(verbose = FALSE)
}

weight_summary_dt <- rbindlist(weight_summaries, use.names = TRUE, fill = TRUE)
weight_summary_csv <- file.path(RESULT_DIR, "pooled_ps_weights_summary.csv")
fwrite(weight_summary_dt, weight_summary_csv)
module_log_info(sprintf("权重摘要已保存: %s", weight_summary_csv))

# ---- 8. 保存 Pooled PS 模型 ----
ps_model_path <- file.path(RESULT_DIR, "pooled_ps_model.rds")
saveRDS(ps_model, ps_model_path)
module_log_info(sprintf("Pooled PS 模型已保存: %s", ps_model_path))

# ---- 9. 保存 PS 模型诊断 ----
ps_diagnostics <- data.table(
  n_training = nrow(dt_cc),
  n_covariates = length(covs_use),
  aic = AIC(ps_model),
  converged = ps_model$converged,
  p_treat = p_treat,
  iptw_mean = mean(iptw_trunc),
  iptw_sd = sd(iptw_trunc),
  ess = ess,
  runtime_sec = round(ps_runtime, 2),
  peak_mem_mb = get_peak_mem_mb()
)
ps_diag_csv <- file.path(RESULT_DIR, "pooled_ps_diagnostics.csv")
fwrite(ps_diagnostics, ps_diag_csv)
module_log_info(sprintf("PS 诊断已保存: %s", ps_diag_csv))

# ---- 10. 保存全量元数据 JSON ----
metadata <- list(
  phase = "01_pooled_ps",
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  input_file = INPUT_PARQUET,
  n_patients_total = length(unique(dt[[ID_COL]])),
  n_records_total = nrow(dt),
  n_records_valid = nrow(dt_valid),
  n_records_cc = nrow(dt_cc),
  n_periods_valid = length(valid_periods),
  valid_periods = as.list(valid_periods),
  n_covariates_input = length(PS_COVARIATES),
  n_covariates_used = length(covs_use),
  covariates_dropped = as.list(dropped),
  iptw_upper_quantile = q_upper,
  iptw_lower_quantile = q_lower,
  ess = ess,
  runtime_sec = round(as.numeric(difftime(Sys.time(), t_start_total, units = "secs")), 2),
  peak_mem_mb = get_peak_mem_mb()
)
metadata_json <- file.path(RESULT_DIR, "01_pooled_ps_metadata.json")
writeLines(toJSON(metadata, pretty = TRUE, auto_unbox = TRUE), metadata_json)

module_log_section("Phase 1 完成")
module_log_info(sprintf("总运行时间: %.2f 秒", metadata$runtime_sec))
module_log_info(sprintf("峰值内存: %s MB", ifelse(is.na(metadata$peak_mem_mb), "NA", as.character(metadata$peak_mem_mb))))

# 清理
rm(dt, dt_valid, dt_cc, ps_model); gc(verbose = FALSE)
