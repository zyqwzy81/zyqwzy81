# ============================================================================
# 02_batch_ps_iptw.R
# SEQTaRget 分批引擎 —— 步骤二：按 trial_period 分批拟合 PS 并计算 IPTW
# 输出：加权后的分批 Parquet + 每批次 PS/IPTW 摘要 CSV
# ============================================================================

source("config.R")

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
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

# ---- 确保字符型协变量可进入 glm 模型 ----
# glm 会自动将 factor/character 转为 dummy；这里统一转为 factor
# 同时剔除该批次内零方差或单水平的协变量
prepare_covariates <- function(dt, covs) {
  covs_avail <- intersect(covs, names(dt))
  keep <- c()
  for (col in covs_avail) {
    x <- dt[[col]]
    if (is.character(x)) {
      x <- as.factor(x)
      dt[[col]] <- x
    }
    # 剔除数值型零方差或因子型单水平
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

module_log_section("步骤二：分批 PS/IPTW")

# ---- 发现所有 period 分区 ----
period_dirs <- list.dirs(EXPANDED_DIR, recursive = FALSE)
periods <- sort(as.integer(gsub("trial_period=", "", basename(period_dirs))))
if (length(periods) == 0) {
  stop(sprintf("未在 %s 发现扩展数据分区，请先运行 01_expand_only.R", EXPANDED_DIR))
}
module_log_info(sprintf("发现 %d 个 period 分区", length(periods)))

batch_summaries <- list()
t_start_total <- Sys.time()

for (p in periods) {
  t_start <- Sys.time()
  in_path <- file.path(EXPANDED_DIR, sprintf("trial_period=%d", p), "part-0.parquet")
  if (!file.exists(in_path)) {
    module_log_warn(sprintf("period=%d 数据文件不存在，跳过", p))
    next
  }

  batch <- as.data.table(arrow::read_parquet(in_path))
  batch <- batch[eligible == 1L]
  if (nrow(batch) < 30 || sum(batch$outcome, na.rm = TRUE) < 5 ||
      length(unique(batch$tx_init)) < 2) {
    module_log_warn(sprintf("period=%d 样本量/事件数/治疗变异不足，跳过", p))
    rm(batch)
    gc(verbose = FALSE)
    next
  }

  covs_use <- prepare_covariates(batch, PS_COVARIATES)
  dropped <- setdiff(PS_COVARIATES, covs_use)
  if (length(dropped) > 0) {
    module_log_warn(sprintf("period=%d 剔除零方差协变量: %s", p, paste(dropped, collapse = ", ")))
  }
  if (length(covs_use) == 0) {
    module_log_warn(sprintf("period=%d 无可用协变量，跳过", p))
    rm(batch)
    gc(verbose = FALSE)
    next
  }
  # 仅保留完整案例用于 PS 模型（与 SEQTaRget 默认一致）
  cc_idx <- complete.cases(batch[, ..covs_use])
  batch_cc <- batch[cc_idx]
  if (nrow(batch_cc) == 0) {
    module_log_warn(sprintf("period=%d 完整案例数为 0，跳过", p))
    rm(batch, batch_cc)
    gc(verbose = FALSE)
    next
  }

  # ---- 拟合 logistic PS 模型 ----
  formula_str <- paste("tx_init ~", paste(covs_use, collapse = " + "))
  ps_model <- glm(as.formula(formula_str), data = batch_cc, family = binomial())
  ps <- predict(ps_model, type = "response")
  ps <- pmax(pmin(ps, PS_CLIP[2]), PS_CLIP[1])

  # ---- 稳定化 IPTW ----
  p_treat <- mean(batch_cc$tx_init, na.rm = TRUE)
  iptw <- ifelse(batch_cc$tx_init == 1,
                 p_treat / ps,
                 (1 - p_treat) / (1 - ps))

  # ---- 双侧截断 ----
  q_upper <- quantile(iptw, IPTW_TRUNC_UPPER, na.rm = TRUE)
  q_lower <- quantile(iptw, IPTW_TRUNC_LOWER, na.rm = TRUE)
  iptw_trunc <- pmax(pmin(iptw, q_upper), q_lower)

  batch_cc[, ps := ps]
  batch_cc[, iptw := iptw]
  batch_cc[, iptw_trunc := iptw_trunc]

  # ---- 保存加权批次 ----
  out_path <- file.path(WEIGHTED_DIR, sprintf("trial_period=%d", p))
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  arrow::write_parquet(batch_cc, file.path(out_path, "part-0.parquet"))

  # ---- 批次摘要 ----
  summ <- data.table(
    trial_period = p,
    n_clones = nrow(batch_cc),
    n_patients = length(unique(batch_cc$stay_id)),
    n_events = sum(batch_cc$outcome, na.rm = TRUE),
    treatment_rate = round(mean(batch_cc$tx_init, na.rm = TRUE), 4),
    ps_mean = round(mean(batch_cc$ps, na.rm = TRUE), 4),
    ps_sd = round(sd(batch_cc$ps, na.rm = TRUE), 4),
    iptw_mean = round(mean(batch_cc$iptw, na.rm = TRUE), 4),
    iptw_sd = round(sd(batch_cc$iptw, na.rm = TRUE), 4),
    iptw_trunc_mean = round(mean(batch_cc$iptw_trunc, na.rm = TRUE), 4),
    iptw_trunc_sd = round(sd(batch_cc$iptw_trunc, na.rm = TRUE), 4),
    iptw_min = round(min(batch_cc$iptw_trunc, na.rm = TRUE), 4),
    iptw_max = round(max(batch_cc$iptw_trunc, na.rm = TRUE), 4),
    ess = round(sum(batch_cc$iptw_trunc)^2 / sum(batch_cc$iptw_trunc^2), 1),
    runtime_sec = round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 2),
    peak_mem_mb = get_peak_mem_mb(),
    ps_model_converged = ps_model$converged
  )
  batch_summaries[[as.character(p)]] <- summ

  module_log_info(sprintf(
    "  period=%d: n=%d, events=%d, ps_mean=%.3f, iptw_mean=%.3f, ess=%.1f, time=%.2fs",
    p, summ$n_clones, summ$n_events, summ$ps_mean, summ$iptw_trunc_mean,
    summ$ess, summ$runtime_sec
  ))

  rm(batch, batch_cc, ps_model, ps, iptw, iptw_trunc, summ)
  gc(verbose = FALSE)
}

summary_dt <- rbindlist(batch_summaries, use.names = TRUE, fill = TRUE)
summary_csv <- file.path(RESULTS_DIR, "batch_ps_iptw_summary.csv")
fwrite(summary_dt, summary_csv)

module_log_section("PS/IPTW 步骤完成")
module_log_info(sprintf("共处理 %d 个批次", nrow(summary_dt)))
module_log_info(sprintf("总体 ps_mean: %.3f, 总体 iptw_trunc_mean: %.3f",
                        weighted.mean(summary_dt$ps_mean, summary_dt$n_clones, na.rm = TRUE),
                        weighted.mean(summary_dt$iptw_trunc_mean, summary_dt$n_clones, na.rm = TRUE)))
module_log_info(sprintf("峰值内存: %s MB, 运行时间: %.2f 秒",
                        ifelse(is.na(get_peak_mem_mb()), "NA", as.character(get_peak_mem_mb())),
                        as.numeric(difftime(Sys.time(), t_start_total, units = "secs"))))
module_log_info(sprintf("批次摘要 CSV: %s", summary_csv))

rm(summary_dt, batch_summaries); gc(verbose = FALSE)
