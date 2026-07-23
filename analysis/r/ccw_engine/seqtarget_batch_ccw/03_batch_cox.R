# ============================================================================
# 03_batch_cox.R
# SEQTaRget 分批引擎 —— 步骤三：按 trial_period 分批拟合 IPTW-Cox
# 输出：每期 HR + 汇总 HR 的 CSV 文件
# ============================================================================

source("config.R")

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(survival)
})

get_peak_mem_mb <- function() {
  if (file.exists("/proc/self/status")) {
    lines <- readLines("/proc/self/status", warn = FALSE)
    peak_line <- grep("^VmPeak:", lines, value = TRUE)
    if (length(peak_line) > 0) {
      kb <- as.numeric(gsub("[^0-9]", "", peak_line))
      return(round(kb / 1024, 2))
    }
  }
  NA_real_
}

# ---- 固定效应逆方差合并 ----
pool_fixed <- function(log_hr, se, weights = NULL) {
  if (is.null(weights)) weights <- 1 / se^2
  weights <- weights[!is.na(log_hr) & !is.na(se) & se > 0]
  log_hr <- log_hr[!is.na(log_hr) & !is.na(se) & se > 0]
  se <- se[!is.na(log_hr) & !is.na(se) & se > 0]
  pooled <- sum(log_hr * weights) / sum(weights)
  pooled_se <- sqrt(1 / sum(weights))
  ci_l <- pooled - 1.96 * pooled_se
  ci_u <- pooled + 1.96 * pooled_se
  list(log_hr = pooled, se = pooled_se,
       hr = exp(pooled), hr_lower = exp(ci_l), hr_upper = exp(ci_u))
}

# ---- 随机效应 DerSimonian-Laird ----
pool_random <- function(log_hr, se) {
  if (length(log_hr) < 2) return(pool_fixed(log_hr, se))
  # 若 meta 包可用则使用 metagen，否则手动 DL
  if (requireNamespace("meta", quietly = TRUE)) {
    res <- meta::metagen(log_hr, se, sm = "HR", method.tau = "DL")
    return(list(
      log_hr = res$TE.random,
      se = res$seTE.random,
      hr = exp(res$TE.random),
      hr_lower = exp(res$lower.random),
      hr_upper = exp(res$upper.random)
    ))
  }
  # 手动 DL 实现
  wi <- 1 / se^2
  theta_fe <- sum(log_hr * wi) / sum(wi)
  Q <- sum(wi * (log_hr - theta_fe)^2)
  df <- length(log_hr) - 1
  C <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- max(0, (Q - df) / C)
  wi_star <- 1 / (se^2 + tau2)
  theta_re <- sum(log_hr * wi_star) / sum(wi_star)
  se_re <- sqrt(1 / sum(wi_star))
  list(log_hr = theta_re, se = se_re,
       hr = exp(theta_re), hr_lower = exp(theta_re - 1.96 * se_re),
       hr_upper = exp(theta_re + 1.96 * se_re))
}

module_log_section("步骤三：分批 Cox 拟合")

period_dirs <- list.dirs(WEIGHTED_DIR, recursive = FALSE)
periods <- sort(as.integer(gsub("trial_period=", "", basename(period_dirs))))
if (length(periods) == 0) {
  stop(sprintf("未在 %s 发现加权批次，请先运行 02_batch_ps_iptw.R", WEIGHTED_DIR))
}
module_log_info(sprintf("发现 %d 个加权批次", length(periods)))

period_results <- list()
t_start_total <- Sys.time()

for (p in periods) {
  t_start <- Sys.time()
  in_path <- file.path(WEIGHTED_DIR, sprintf("trial_period=%d", p), "part-0.parquet")
  if (!file.exists(in_path)) next

  batch <- as.data.table(arrow::read_parquet(in_path))
  batch <- batch[eligible == 1L & followup_time > 0]
  if (nrow(batch) == 0 || sum(batch$outcome, na.rm = TRUE) < 5 ||
      length(unique(batch$tx_init)) < 2) {
    module_log_warn(sprintf("period=%d 数据不足以拟合 Cox，跳过", p))
    rm(batch); gc(verbose = FALSE)
    next
  }

  fit <- tryCatch({
    coxph(Surv(followup_time, outcome) ~ tx_init,
          data = batch,
          weights = iptw_trunc,
          cluster = stay_id,
          model = FALSE,
          x = FALSE,
          y = FALSE)
  }, error = function(e) {
    module_log_warn(sprintf("period=%d Cox 拟合失败: %s", p, e$message))
    NULL
  })

  if (is.null(fit)) {
    rm(batch); gc(verbose = FALSE)
    next
  }

  s <- summary(fit)
  coef_idx <- "tx_init"
  log_hr <- s$coefficients[coef_idx, "coef"]
  se <- s$coefficients[coef_idx, "robust se"]
  hr <- exp(log_hr)
  ci <- confint(fit)
  ci_l <- exp(ci[coef_idx, 1])
  ci_u <- exp(ci[coef_idx, 2])
  pval <- s$coefficients[coef_idx, "Pr(>|z|)"]

  res <- data.table(
    trial_period = p,
    n_clones = nrow(batch),
    n_patients = length(unique(batch$stay_id)),
    n_events = sum(batch$outcome, na.rm = TRUE),
    log_hr = log_hr,
    se = se,
    hr = hr,
    hr_lower = ci_l,
    hr_upper = ci_u,
    p_value = pval,
    runtime_sec = round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 2),
    peak_mem_mb = get_peak_mem_mb()
  )
  period_results[[as.character(p)]] <- res

  module_log_info(sprintf(
    "  period=%d: HR=%.3f (%.3f-%.3f), events=%d, n=%d, time=%.2fs",
    p, hr, ci_l, ci_u, res$n_events, res$n_clones, res$runtime_sec
  ))

  rm(batch, fit, res)
  gc(verbose = FALSE)
}

results_dt <- rbindlist(period_results, use.names = TRUE, fill = TRUE)
period_csv <- file.path(RESULTS_DIR, "seqtarget_batch_hr_period_specific.csv")
fwrite(results_dt, period_csv)

# ---- 按配置方法汇总 ----
pool_list <- list()
if (nrow(results_dt) > 0) {
  if (POOL_METHOD == "period_specific") {
    pool_list[["period_specific"]] <- data.table(
      method = "period_specific",
      n_periods = nrow(results_dt),
      log_hr = NA_real_, se = NA_real_,
      hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
      note = "详见 seqtarget_batch_hr_period_specific.csv"
    )
  } else if (POOL_METHOD == "fixed") {
    pooled <- pool_fixed(results_dt$log_hr, results_dt$se,
                         weights = results_dt$n_events)
    pool_list[["fixed"]] <- data.table(
      method = "fixed",
      n_periods = nrow(results_dt),
      log_hr = pooled$log_hr,
      se = pooled$se,
      hr = pooled$hr,
      hr_lower = pooled$hr_lower,
      hr_upper = pooled$hr_upper,
      note = sprintf("以 n_events 为权重的逆方差固定效应合并")
    )
  } else if (POOL_METHOD == "random") {
    pooled <- pool_random(results_dt$log_hr, results_dt$se)
    pool_list[["random"]] <- data.table(
      method = "random",
      n_periods = nrow(results_dt),
      log_hr = pooled$log_hr,
      se = pooled$se,
      hr = pooled$hr,
      hr_lower = pooled$hr_lower,
      hr_upper = pooled$hr_upper,
      note = "DerSimonian-Laird 随机效应合并"
    )
  } else {
    module_log_warn(sprintf("未知的 POOL_METHOD=%s，默认使用 period_specific", POOL_METHOD))
    pool_list[["period_specific"]] <- data.table(
      method = "period_specific",
      n_periods = nrow(results_dt),
      log_hr = NA_real_, se = NA_real_,
      hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
      note = "POOL_METHOD 未识别，默认返回按期表格"
    )
  }
}

summary_dt <- rbindlist(pool_list, use.names = TRUE, fill = TRUE)
summary_csv <- file.path(RESULTS_DIR, "seqtarget_batch_hr_summary.csv")
fwrite(summary_dt, summary_csv)

module_log_section("Cox 拟合步骤完成")
module_log_info(sprintf("成功拟合 %d 个 period", nrow(results_dt)))
if (nrow(summary_dt) > 0 && !is.na(summary_dt$hr[1])) {
  module_log_info(sprintf("汇总 HR (%s): %.3f (%.3f-%.3f)",
                          summary_dt$method[1], summary_dt$hr[1],
                          summary_dt$hr_lower[1], summary_dt$hr_upper[1]))
}
module_log_info(sprintf("峰值内存: %s MB, 运行时间: %.2f 秒",
                        ifelse(is.na(get_peak_mem_mb()), "NA", as.character(get_peak_mem_mb())),
                        as.numeric(difftime(Sys.time(), t_start_total, units = "secs"))))
module_log_info(sprintf("按期结果 CSV: %s", period_csv))
module_log_info(sprintf("汇总结果 CSV: %s", summary_csv))

rm(results_dt, summary_dt, period_results); gc(verbose = FALSE)
