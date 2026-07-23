# ============================================================================
# 03_batched_logistic.R — Phase 3: Batched Pooled Logistic
#
# 核心验证机制: 累加式分批拟合 pooled logistic
#   按 trial_period 逐批读取展开数据 → 累加到 cumulative 数据集 →
#   每批拟合 complete-case pooled logistic（含 followup² + trial²）
#   → 记录每期的: β_tx, SE, HR, 95% CI
#   → 关键: 验证 HR 随 period 累加是否收敛到全量 HR
#
# 与 SEQTaRget ITT 对齐:
#   logit(P) = β₀ + β₁*tx_init + β₂*followup + β₃*followup²
#            + β₄*trial + β₅*trial² + Σβᵢ*covᵢ
#   Cluster-robust SE (stay_id)
#
# 输出:
#   - cumulative_hr.csv: 逐期累加 HR（收敛曲线）
#   - period_specific_hr.csv: 每期独立 HR（若事件数足够）
#   - final_model.rds: 全量数据最终模型
# ============================================================================

source("config_v2.R")

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(sandwich)
  library(jsonlite)
  if (requireNamespace("ggplot2", quietly = TRUE)) library(ggplot2)
})

module_log_section("Phase 3: Batched Pooled Logistic（累加式验证）")

t_start_total <- Sys.time()

# ---- 1. 发现展开分区 ----
period_dirs <- list.dirs(EXPAND_DIR, recursive = FALSE)
periods <- sort(as.integer(gsub("trial_period=", "", basename(period_dirs))))
if (length(periods) == 0) {
  stop(sprintf("未在 %s 发现展开分区，请先运行 02_batched_seqexpand.R", EXPAND_DIR))
}
module_log_info(sprintf("发现 %d 个展开分区: %s", length(periods),
                         paste(periods, collapse = ", ")))

# ---- 2. 辅助函数 ----

# 准备协变量（剔除零方差/单水平因子）
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

# 拟合加权 pooled logistic 并提取治疗效应
fit_weighted_logistic <- function(data, covs_use, label = "") {
  # 构建公式
  base_formula <- "y_daily ~ tx_init + followup + followup_sq + trial + trial_sq"
  if (length(covs_use) > 0) {
    cov_formula <- paste(covs_use, collapse = " + ")
    formula_str <- paste(base_formula, cov_formula, sep = " + ")
  } else {
    formula_str <- base_formula
  }
  
  # 完整案例
  cc_vars <- c("y_daily", "tx_init", "followup", "followup_sq",
               "trial", "trial_sq", "iptw_trunc", covs_use)
  cc_idx <- complete.cases(data[, ..cc_vars])
  data_cc <- data[cc_idx]
  
  if (nrow(data_cc) < 100 || sum(data_cc$y_daily) < 10 ||
      length(unique(data_cc$tx_init)) < 2) {
    module_log_warn(sprintf("  %s: 数据不足（n=%d, events=%d）, 跳过",
                            label, nrow(data_cc), sum(data_cc$y_daily)))
    return(NULL)
  }
  
  # 拟合加权 GLM
  fit <- tryCatch(
    glm(as.formula(formula_str), data = data_cc,
        family = binomial(), weights = iptw_trunc,
        control = list(maxit = 100)),
    error = function(e) {
      module_log_warn(sprintf("  %s: GLM 失败: %s", label, e$message))
      NULL
    }
  )
  if (is.null(fit) || !fit$converged) return(NULL)
  
  # Cluster-robust SE
  coef_tx <- coef(fit)["tx_init"]
  vcov_cl <- tryCatch(
    sandwich::vcovCL(fit, cluster = data_cc[[POOL_LOGISTIC_CLUSTER_BY]]),
    error = function(e) {
      module_log_warn(sprintf("  %s: vcovCL 失败, 回退到 naive SE", label))
      vcov(fit)
    }
  )
  se_tx <- sqrt(diag(vcov_cl))["tx_init"]
  
  # 如果 se 不合法，回退到 naive SE
  if (is.na(se_tx) || se_tx <= 0 || !is.finite(se_tx)) {
    se_tx <- sqrt(diag(vcov(fit)))["tx_init"]
  }
  
  hr <- exp(coef_tx)
  ci_lower <- exp(coef_tx - 1.96 * se_tx)
  ci_upper <- exp(coef_tx + 1.96 * se_tx)
  z_val <- coef_tx / se_tx
  p_val <- 2 * pnorm(-abs(z_val))
  
  list(
    label = label,
    n_rows = nrow(data_cc),
    n_events = sum(data_cc$y_daily),
    n_patients = length(unique(data_cc[[POOL_LOGISTIC_CLUSTER_BY]])),
    n_periods = length(unique(data_cc$trial_period)),
    log_hr = coef_tx,
    se = se_tx,
    hr = hr,
    hr_lower = ci_lower,
    hr_upper = ci_upper,
    p_value = p_val,
    converged = fit$converged
  )
}

# ---- 3. 逐期累加拟合 ----
cumulative_results <- list()
cumulative_data <- NULL

for (i in seq_along(periods)) {
  p <- periods[i]
  t_start <- Sys.time()
  
  # 3a: 读取该期展开数据
  in_path <- file.path(EXPAND_DIR, sprintf("trial_period=%d", p), "part-0.parquet")
  if (!file.exists(in_path)) {
    module_log_warn(sprintf("period=%d: 展开文件不存在，跳过", p))
    next
  }
  
  batch_data <- as.data.table(arrow::read_parquet(in_path))
  module_log_info(sprintf("period=%d: 读取 %d 行展开数据", p, nrow(batch_data)))
  
  # 3b: 累加到 cumulative 数据集
  if (is.null(cumulative_data)) {
    cumulative_data <- batch_data
  } else {
    cumulative_data <- rbind(cumulative_data, batch_data, use.names = TRUE, fill = TRUE)
  }
  
  # 3c: 检查 cumulative 数据是否有足够变异
  covs_use <- prepare_covariates(cumulative_data, PS_COVARIATES)
  
  n_events_total <- sum(cumulative_data$y_daily, na.rm = TRUE)
  n_tx_var <- length(unique(cumulative_data$tx_init))
  
  module_log_info(sprintf(
    "  cumulative period=1-%d: %d 行, %d 事件, 治疗变异=%d, 协变量=%d",
    p, nrow(cumulative_data), n_events_total, n_tx_var, length(covs_use)
  ))
  
  # 3d: 拟合累积 pooled logistic
  cum_result <- fit_weighted_logistic(
    cumulative_data, covs_use,
    label = sprintf("cumulative_1-%d", p)
  )
  
  if (!is.null(cum_result)) {
    cum_result$trial_period_cumulative <- p
    cum_result$runtime_sec <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 3)
    cumulative_results[[as.character(p)]] <- as.data.table(cum_result)
    
    module_log_info(sprintf(
      "  >>> period=1-%d: HR=%.4f (%.4f-%.4f), p=%.4f, n_events=%d, SE=%.4f",
      p, cum_result$hr, cum_result$hr_lower, cum_result$hr_upper,
      cum_result$p_value, cum_result$n_events, cum_result$se
    ))
  } else {
    module_log_warn(sprintf("  >>> period=1-%d: 模型未收敛", p))
  }
  
  rm(batch_data); gc(verbose = FALSE)
}

# ---- 4. 所有 period 汇总（全量模型）----
module_log_section("全量数据最终模型")
covs_use_full <- prepare_covariates(cumulative_data, PS_COVARIATES)
final_result <- fit_weighted_logistic(
  cumulative_data, covs_use_full,
  label = "full"
)

module_log_info(sprintf(
  "最终模型: HR=%.4f (%.4f-%.4f), p=%.4f, n_events=%d, n_periods=%d",
  final_result$hr, final_result$hr_lower, final_result$hr_upper,
  final_result$p_value, final_result$n_events, final_result$n_periods
))

# ---- 5. 保存结果 ----

# 5a: 累积 HR CSV
cumulative_dt <- rbindlist(cumulative_results, use.names = TRUE, fill = TRUE)
cumulative_csv <- file.path(RESULT_DIR, "cumulative_hr_by_period.csv")
fwrite(cumulative_dt, cumulative_csv)
module_log_info(sprintf("累积 HR 已保存: %s", cumulative_csv))

# 5b: 最终 HR 摘要
final_summary <- data.table(
  model = "cumulative_full",
  n_periods = final_result$n_periods,
  n_rows = final_result$n_rows,
  n_events = final_result$n_events,
  n_patients = final_result$n_patients,
  hr = final_result$hr,
  hr_lower = final_result$hr_lower,
  hr_upper = final_result$hr_upper,
  log_hr = final_result$log_hr,
  se = final_result$se,
  p_value = final_result$p_value
)
final_csv <- file.path(RESULT_DIR, "final_hr_summary.csv")
fwrite(final_summary, final_csv)
module_log_info(sprintf("最终 HR 摘要: %s", final_csv))

# 5c: 保存最终模型 RDS（可选）
final_model_path <- file.path(RESULT_DIR, "final_pooled_logistic_model.rds")
# 不保存全量 GLM 对象（太大），仅保存关键结果
saveRDS(final_result, final_model_path)

# ---- 6. 收敛性可视化（如可用）----
if (requireNamespace("ggplot2", quietly = TRUE) && nrow(cumulative_dt) > 0) {
  p <- ggplot(cumulative_dt, aes(x = trial_period_cumulative, y = hr)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    geom_ribbon(aes(ymin = hr_lower, ymax = hr_upper), alpha = 0.2, fill = "#0072B5") +
    geom_line(color = "#BC3C29", linewidth = 1) +
    geom_point(color = "#BC3C29", size = 2.5) +
    scale_x_continuous(breaks = cumulative_dt$trial_period_cumulative) +
    labs(
      title = "Cumulative HR Convergence by Trial Period",
      subtitle = sprintf("Pooled logistic with IPTW (DAG12 covariates, n_events=%d)",
                         final_result$n_events),
      x = "Cumulative trial period",
      y = "Hazard Ratio (95% CI)",
      caption = "Each point = model fitted on all periods up to that point"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey40", hjust = 0)
    )
  
  ggsave(file.path(RESULT_DIR, "Fig_cumulative_hr_convergence.pdf"), p, width = 8, height = 5)
  ggsave(file.path(RESULT_DIR, "Fig_cumulative_hr_convergence.png"), p, width = 8, height = 5, dpi = 150)
  module_log_info("收敛性图已保存")
}

# ---- 7. 元数据 ----
metadata <- list(
  phase = "03_batched_logistic",
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  expand_dir = EXPAND_DIR,
  n_periods_used = length(periods),
  n_cumulative_fits = nrow(cumulative_dt),
  final_hr = final_result$hr,
  final_hr_lower = final_result$hr_lower,
  final_hr_upper = final_result$hr_upper,
  final_p_value = final_result$p_value,
  final_n_events = final_result$n_events,
  final_n_patients = final_result$n_patients,
  runtime_sec = round(as.numeric(difftime(Sys.time(), t_start_total, units = "secs")), 2),
  peak_mem_mb = get_peak_mem_mb(),
  output_files = list(
    cumulative_hr = cumulative_csv,
    final_summary = final_csv,
    final_model = final_model_path
  )
)
metadata_json <- file.path(RESULT_DIR, "03_logistic_metadata.json")
writeLines(toJSON(metadata, pretty = TRUE, auto_unbox = TRUE), metadata_json)

module_log_section("Phase 3 完成")
module_log_info(sprintf("总运行时间: %.2f 秒", metadata$runtime_sec))
module_log_info(sprintf("峰值内存: %s MB", ifelse(is.na(metadata$peak_mem_mb), "NA", as.character(metadata$peak_mem_mb))))

# 清理
rm(cumulative_data, cumulative_results); gc(verbose = FALSE)
