# ============================================================================
# 04_compare_with_seqtarget.R
# SEQTaRget 分批引擎 —— 步骤四：小样本与官方 SEQTaRget 对比验证
# 输出：对比 CSV + HR 比较图（PDF/PNG）
# ============================================================================

source("config.R")

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(survival)
  if (requireNamespace("sandwich", quietly = TRUE)) library(sandwich)
  if (requireNamespace("ggplot2", quietly = TRUE)) library(ggplot2)
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

# ---- 协变量预处理：剔除批次内零方差或单水平因子 ----
prepare_covariates <- function(dt, covs) {
  covs_avail <- intersect(covs, names(dt))
  keep <- c()
  for (col in covs_avail) {
    x <- dt[[col]]
    if (is.character(x)) {
      x <- as.factor(x)
      dt[[col]] <- x
    }
    if (is.numeric(x)) {
      x_sd <- sd(x, na.rm = TRUE)
      if (!is.na(x_sd) && x_sd > 1e-12) keep <- c(keep, col)
    } else if (is.factor(x)) {
      if (length(levels(x)) >= 2) keep <- c(keep, col)
    } else {
      keep <- c(keep, col)
    }
  }
  keep
}

# ---- 运行本引擎（01-03 步骤）在指定子集上 ----
run_batch_engine <- function(dt_subset, out_tag) {
  t_start <- Sys.time()

  cov_cols <- intersect(c(PS_COVARIATES, "gender", "admission_type"), names(dt_subset))
  expanded <- data.table(
    stay_id = dt_subset[[ID_COL]],
    subject_id = if ("subject_id" %in% names(dt_subset)) dt_subset$subject_id else NA_integer_,
    trial_period = dt_subset[[PERIOD_COL]],
    time_zero = if ("time_zero" %in% names(dt_subset)) dt_subset$time_zero else NA_character_,
    baseline_phosphate = if ("baseline_phosphate" %in% names(dt_subset)) dt_subset$baseline_phosphate else NA_real_,
    tx_init = as.integer(dt_subset[[TREATMENT_SRC_COL]] == TREATMENT_LABEL),
    outcome = as.integer(dt_subset[[OUTCOME_COL]]),
    followup_time = as.numeric(dt_subset[[SURV_TIME_COL]]),
    eligible = 1L
  )
  for (col in cov_cols) {
    data.table::set(expanded, j = col, value = dt_subset[[col]])
  }
  expanded[, clone_id := .I]
  expanded[followup_time <= 0 & outcome == 0, eligible := 0L]

  periods <- sort(unique(expanded$trial_period))
  period_results <- list()
  overall_iptw <- c()
  overall_ps <- c()

  for (p in periods) {
    batch <- expanded[trial_period == p & eligible == 1L]
    if (nrow(batch) == 0) next

    covs_use <- prepare_covariates(batch, PS_COVARIATES)
    dropped <- setdiff(PS_COVARIATES, covs_use)
    if (length(dropped) > 0) {
      module_log_warn(sprintf("batch engine period=%d 剔除零方差协变量: %s",
                              p, paste(dropped, collapse = ", ")))
    }
    if (length(covs_use) == 0) {
      module_log_warn(sprintf("batch engine period=%d 无可用协变量，跳过", p))
      rm(batch)
      gc(verbose = FALSE)
      next
    }
    cc_idx <- complete.cases(batch[, ..covs_use])
    batch_cc <- batch[cc_idx]
    if (nrow(batch_cc) == 0) {
      module_log_warn(sprintf("batch engine period=%d 完整案例数为 0，跳过", p))
      rm(batch, batch_cc)
      gc(verbose = FALSE)
      next
    }

    formula_str <- paste("tx_init ~", paste(covs_use, collapse = " + "))
    ps_model <- tryCatch(
      glm(as.formula(formula_str), data = batch_cc, family = binomial()),
      error = function(e) {
        module_log_warn(sprintf("batch engine period=%d PS 模型失败: %s", p, e$message))
        NULL
      }
    )
    if (is.null(ps_model)) {
      rm(batch, batch_cc)
      gc(verbose = FALSE)
      next
    }
    ps <- predict(ps_model, type = "response")
    ps <- pmax(pmin(ps, PS_CLIP[2]), PS_CLIP[1])
    p_treat <- mean(batch_cc$tx_init, na.rm = TRUE)
    iptw <- ifelse(batch_cc$tx_init == 1, p_treat / ps, (1 - p_treat) / (1 - ps))
    q_upper <- quantile(iptw, IPTW_TRUNC_UPPER, na.rm = TRUE)
    q_lower <- quantile(iptw, IPTW_TRUNC_LOWER, na.rm = TRUE)
    iptw_trunc <- pmax(pmin(iptw, q_upper), q_lower)

    batch_cc[, iptw_trunc := iptw_trunc]
    overall_ps <- c(overall_ps, ps)
    overall_iptw <- c(overall_iptw, iptw_trunc)

    if (sum(batch_cc$outcome, na.rm = TRUE) < 5 ||
        length(unique(batch_cc$tx_init)) < 2) next

    fit <- tryCatch(
      coxph(Surv(followup_time, outcome) ~ tx_init,
            data = batch_cc, weights = iptw_trunc, cluster = stay_id),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    s <- summary(fit)
    log_hr <- s$coefficients["tx_init", "coef"]
    se <- s$coefficients["tx_init", "robust se"]
    ci <- confint(fit)
    period_results[[as.character(p)]] <- data.table(
      trial_period = p,
      n_clones = nrow(batch_cc),
      n_events = sum(batch_cc$outcome, na.rm = TRUE),
      log_hr = log_hr,
      se = se,
      hr = exp(log_hr),
      hr_lower = exp(ci["tx_init", 1]),
      hr_upper = exp(ci["tx_init", 2])
    )
  }

  results_dt <- rbindlist(period_results, use.names = TRUE, fill = TRUE)
  pooled <- NULL
  if (nrow(results_dt) > 0) {
    weights <- results_dt$n_events
    if (all(is.na(weights)) || all(weights == 0, na.rm = TRUE)) {
      weights <- results_dt$n_clones
    }
    if (all(is.na(weights)) || all(weights == 0, na.rm = TRUE)) {
      weights <- rep(1, nrow(results_dt))
    }
    pooled <- .pool_fixed_internal(results_dt$log_hr, results_dt$se, weights)
  }

  runtime <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  list(
    method = "batch_engine",
    tag = out_tag,
    hr = ifelse(is.null(pooled), NA_real_, pooled$hr),
    hr_lower = ifelse(is.null(pooled), NA_real_, pooled$hr_lower),
    hr_upper = ifelse(is.null(pooled), NA_real_, pooled$hr_upper),
    se = ifelse(is.null(pooled), NA_real_, pooled$se),
    mean_iptw = ifelse(length(overall_iptw) > 0, mean(overall_iptw, na.rm = TRUE), NA_real_),
    ps_mean = ifelse(length(overall_ps) > 0, mean(overall_ps, na.rm = TRUE), NA_real_),
    runtime_sec = runtime,
    peak_mem_mb = get_peak_mem_mb(),
    n_periods = nrow(results_dt)
  )
}

.pool_fixed_internal <- function(log_hr, se, weights = NULL) {
  valid <- !is.na(log_hr) & !is.na(se) & se > 0
  log_hr <- log_hr[valid]; se <- se[valid]
  if (is.null(weights)) weights <- 1 / se^2 else weights <- weights[valid]
  pooled <- sum(log_hr * weights) / sum(weights)
  pooled_se <- sqrt(1 / sum(weights))
  list(log_hr = pooled, se = pooled_se,
       hr = exp(pooled),
       hr_lower = exp(pooled - 1.96 * pooled_se),
       hr_upper = exp(pooled + 1.96 * pooled_se))
}

# ---- 运行 SEQTaRget 展开式 Batch Engine（日粒度 pooled logistic） ----
# 差异说明：
# 原始 batch engine（run_batch_engine）：每 patient-period = 1 clone → 每期独立 Cox → Meta 分析
# SEQTaRget（SEQexpand）：每 patient-period 展开为日粒度 → pooled logistic（含 followup² + trial²）
#
# 本函数模拟 SEQTaRget 的 SEQexpand：
# - 对每 patient-period，从 day=0 到 floor(surv_time) 生成日记录
# - 死亡日 outcome=1，其余日 outcome=0
# - 治疗分配固定为 ITT 赋值
# - PS/IPTW 按 period 分别拟合（与原始 batch engine 一致）
# - 最终模型：pooled logistic + followup² + trial² + cluster-robust SE
run_batch_engine_seqexpand <- function(dt_subset, out_tag) {
  t_start <- Sys.time()

  cov_cols <- intersect(c(PS_COVARIATES, "gender", "admission_type"), names(dt_subset))
  expanded <- data.table(
    stay_id = dt_subset[[ID_COL]],
    trial_period = dt_subset[[PERIOD_COL]],
    tx_init = as.integer(dt_subset[[TREATMENT_SRC_COL]] == TREATMENT_LABEL),
    outcome = as.integer(dt_subset[[OUTCOME_COL]]),
    surv_time = as.numeric(dt_subset[[SURV_TIME_COL]]),
    eligible = 1L
  )
  for (col in cov_cols) {
    data.table::set(expanded, j = col, value = dt_subset[[col]])
  }
  expanded[surv_time <= 0 & outcome == 0, eligible := 0L]
  expanded <- expanded[eligible == 1L]

  if (nrow(expanded) == 0) {
    return(list(method = "batch_seqexpand", tag = out_tag,
                hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
                se = NA_real_, mean_iptw = NA_real_, ps_mean = NA_real_,
                runtime_sec = NA_real_, peak_mem_mb = NA_real_, n_periods = 0L))
  }

  # ---- 日粒度展开（模拟 SEQexpand 逻辑） ----
  # 对每 patient-period，从 day=0 到 floor(surv_time) 生成日记录
  expand_seqexpand <- function(dt) {
    dt[, day_max := as.integer(floor(pmax(surv_time, 1)))]
    dt[day_max > 28L, day_max := 28L]           # 最多展开到 28 天
    result <- dt[rep(seq_len(.N), day_max + 1L)]  # 按日展开
    result[, day := seq_len(.N) - 1L, by = .(stay_id, trial_period)]
    result[, y_daily := 0L]
    result[day == day_max & outcome == 1, y_daily := 1L]
    result[, c("surv_time", "day_max", "outcome") := NULL]
    result[]
  }

  person_dt <- expand_seqexpand(expanded)
  rm(expanded); gc(verbose = FALSE)

  # ---- 按 period 拟合 PS + IPTW（与原文 run_batch_engine 一致） ----
  periods <- sort(unique(person_dt$trial_period))
  overall_iptw <- c()
  overall_ps <- c()

  for (p in periods) {
    batch_rows <- which(person_dt$trial_period == p)
    batch <- person_dt[batch_rows]
    covs_use <- prepare_covariates(batch, PS_COVARIATES)
    if (length(covs_use) == 0) next

    cc_idx <- complete.cases(batch[, ..covs_use])
    cc_rows <- batch_rows[cc_idx]
    batch_cc <- batch[cc_idx]
    if (nrow(batch_cc) == 0 || length(unique(batch_cc$tx_init)) < 2) next

    formula_str <- paste("tx_init ~", paste(covs_use, collapse = " + "))
    ps_model <- tryCatch(
      glm(as.formula(formula_str), data = batch_cc, family = binomial()),
      error = function(e) NULL
    )
    if (is.null(ps_model)) next

    ps <- predict(ps_model, type = "response")
    ps <- pmax(pmin(ps, PS_CLIP[2]), PS_CLIP[1])
    p_treat <- mean(batch_cc$tx_init, na.rm = TRUE)
    iptw <- ifelse(batch_cc$tx_init == 1, p_treat / ps, (1 - p_treat) / (1 - ps))
    q_upper <- quantile(iptw, IPTW_TRUNC_UPPER, na.rm = TRUE)
    q_lower <- quantile(iptw, IPTW_TRUNC_LOWER, na.rm = TRUE)
    iptw_trunc <- pmax(pmin(iptw, q_upper), q_lower)
    batch_cc[, iptw_trunc := iptw_trunc]
    overall_ps <- c(overall_ps, ps)
    overall_iptw <- c(overall_iptw, iptw_trunc)

    # 使用保留的原始行索引回写（避免 subset 索引错位）
    person_dt[cc_rows, iptw_trunc := iptw_trunc]
    rm(batch, batch_cc, ps_model); gc(verbose = FALSE)
  }

  person_dt <- person_dt[!is.na(iptw_trunc)]

  if (nrow(person_dt) < 100 || sum(person_dt$y_daily) < 10 ||
      length(unique(person_dt$tx_init)) < 2) {
    return(list(method = "batch_seqexpand", tag = out_tag,
                hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
                se = NA_real_, mean_iptw = if (length(overall_iptw) > 0) mean(overall_iptw, na.rm = TRUE) else NA_real_,
                ps_mean = if (length(overall_ps) > 0) mean(overall_ps, na.rm = TRUE) else NA_real_,
                runtime_sec = as.numeric(difftime(Sys.time(), t_start, units = "secs")),
                peak_mem_mb = get_peak_mem_mb(), n_periods = 0L))
  }

  # ---- 拟合 pooled logistic（与 SEQTaRget ITT 内部一致） ----
  # logit(P) = β₀ + β₁*tx + β₂*followup + β₃*followup² + β₄*trial + β₅*trial² + 协变量
  person_dt[, `:=`(followup = as.numeric(day),
                   followup_sq = as.numeric(day^2),
                   trial = as.numeric(trial_period),
                   trial_sq = as.numeric(trial_period^2))]

  covs_use_all <- prepare_covariates(person_dt, PS_COVARIATES)
  cov_formula <- paste(covs_use_all, collapse = " + ")
  plogit_formula <- as.formula(
    sprintf("y_daily ~ tx_init + followup + followup_sq + trial + trial_sq + %s",
            cov_formula)
  )
  plogit_fit <- tryCatch(
    glm(plogit_formula, data = person_dt,
        family = binomial(), weights = iptw_trunc),
    error = function(e) NULL
  )

  if (is.null(plogit_fit)) {
    plogit_simple <- as.formula(sprintf("y_daily ~ tx_init + %s", cov_formula))
    plogit_fit <- tryCatch(
      glm(plogit_simple, data = person_dt,
          family = binomial(), weights = iptw_trunc),
      error = function(e) NULL
    )
  }

  if (is.null(plogit_fit)) {
    return(list(method = "batch_seqexpand", tag = out_tag,
                hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
                se = NA_real_, mean_iptw = mean(overall_iptw, na.rm = TRUE),
                ps_mean = mean(overall_ps, na.rm = TRUE),
                runtime_sec = as.numeric(difftime(Sys.time(), t_start, units = "secs")),
                peak_mem_mb = get_peak_mem_mb(), n_periods = 0L))
  }

  # ---- cluster-robust SE ----
  coef_tx <- coef(plogit_fit)["tx_init"]
  se_tx <- sqrt(diag(vcovCL(plogit_fit, cluster = person_dt$stay_id)))["tx_init"]

  list(
    method = "batch_seqexpand",
    tag = out_tag,
    hr = exp(coef_tx),
    hr_lower = exp(coef_tx - 1.96 * se_tx),
    hr_upper = exp(coef_tx + 1.96 * se_tx),
    se = se_tx,
    mean_iptw = mean(overall_iptw, na.rm = TRUE),
    ps_mean = mean(overall_ps, na.rm = TRUE),
    runtime_sec = as.numeric(difftime(Sys.time(), t_start, units = "secs")),
    peak_mem_mb = get_peak_mem_mb(),
    n_periods = length(periods)
  )
}

# ---- Pooled PS 版 Batch Seqexpand（在展开前拟合 pooled PS，与 SEQTaRget 对齐） ----
# 差异说明 vs 原始 batch_seqexpand：
# 原始 batch_seqexpand：展开后 → per-period PS → per-period IPTW → pooled logistic
# Pooled 版：        展开前 pooled PS → pooled IPTW → 展开 → 展开后应用权重 → pooled logistic
#
# 这更接近 SEQTaRget 的逻辑（weight.preexpansion = TRUE，pooled PS on pre-expansion data）
run_batch_engine_seqexpand_pooled <- function(dt_subset, out_tag) {
  t_start <- Sys.time()

  cov_cols <- intersect(c(PS_COVARIATES, "gender", "admission_type"), names(dt_subset))

  # ---- 第一步：在展开前（patient-period 级别）拟合 Pooled PS ----
  # 与 SEQTaRget (weight.preexpansion=TRUE) 对齐：
  # 所有 period 合并，在原始数据上拟合一个 pooled logistic PS 模型
  pre_dt <- data.table(
    stay_id = dt_subset[[ID_COL]],
    trial_period = dt_subset[[PERIOD_COL]],
    tx_init = as.integer(dt_subset[[TREATMENT_SRC_COL]] == TREATMENT_LABEL),
    outcome = as.integer(dt_subset[[OUTCOME_COL]]),
    surv_time = as.numeric(dt_subset[[SURV_TIME_COL]])
  )
  for (col in cov_cols) {
    data.table::set(pre_dt, j = col, value = dt_subset[[col]])
  }
  pre_dt <- pre_dt[surv_time > 0]
  if (nrow(pre_dt) == 0) {
    return(list(method = "batch_seqexpand_pooled_ps", tag = out_tag,
                hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
                se = NA_real_, mean_iptw = NA_real_, ps_mean = NA_real_,
                runtime_sec = NA_real_, peak_mem_mb = NA_real_, n_periods = 0L))
  }

  # ---- 剔除无治疗变异的 period（与 SEQTaRget 对齐） ----
  pre_dt[, period_tx_var := length(unique(tx_init)) > 1, by = trial_period]
  pre_dt <- pre_dt[period_tx_var == TRUE]
  pre_dt[, period_tx_var := NULL]

  if (nrow(pre_dt) < 50 || length(unique(pre_dt$tx_init)) < 2) {
    return(list(method = "batch_seqexpand_pooled_ps", tag = out_tag,
                hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
                se = NA_real_, mean_iptw = NA_real_, ps_mean = NA_real_,
                runtime_sec = NA_real_, peak_mem_mb = NA_real_, n_periods = 0L))
  }

  # ---- Pooled PS 拟合（所有 period 合并） ----
  covs_use <- prepare_covariates(pre_dt, PS_COVARIATES)
  if (length(covs_use) == 0) {
    return(list(method = "batch_seqexpand_pooled_ps", tag = out_tag,
                hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
                se = NA_real_, mean_iptw = NA_real_, ps_mean = NA_real_,
                runtime_sec = NA_real_, peak_mem_mb = NA_real_, n_periods = 0L))
  }

  cc_idx <- complete.cases(pre_dt[, ..covs_use])
  pre_cc <- pre_dt[cc_idx]
  if (nrow(pre_cc) < 50 || length(unique(pre_cc$tx_init)) < 2) {
    return(list(method = "batch_seqexpand_pooled_ps", tag = out_tag,
                hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
                se = NA_real_, mean_iptw = NA_real_, ps_mean = NA_real_,
                runtime_sec = NA_real_, peak_mem_mb = NA_real_, n_periods = 0L))
  }

  formula_str <- paste("tx_init ~", paste(covs_use, collapse = " + "))
  ps_model <- tryCatch(
    glm(as.formula(formula_str), data = pre_cc, family = binomial()),
    error = function(e) NULL
  )
  if (is.null(ps_model)) {
    return(list(method = "batch_seqexpand_pooled_ps", tag = out_tag,
                hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
                se = NA_real_, mean_iptw = NA_real_, ps_mean = NA_real_,
                runtime_sec = NA_real_, peak_mem_mb = NA_real_, n_periods = 0L))
  }

  ps <- predict(ps_model, type = "response")
  ps <- pmax(pmin(ps, PS_CLIP[2]), PS_CLIP[1])
  p_treat <- mean(pre_cc$tx_init, na.rm = TRUE)
  iptw <- ifelse(pre_cc$tx_init == 1, p_treat / ps, (1 - p_treat) / (1 - ps))
  q_upper <- quantile(iptw, IPTW_TRUNC_UPPER, na.rm = TRUE)
  q_lower <- quantile(iptw, IPTW_TRUNC_LOWER, na.rm = TRUE)
  iptw_trunc <- pmax(pmin(iptw, q_upper), q_lower)
  pre_cc[, iptw_trunc_pooled := iptw_trunc]

  mean_iptw <- mean(iptw_trunc, na.rm = TRUE)
  ps_mean <- mean(ps, na.rm = TRUE)

  # ---- 第二步：将 pooled PS 权重映射回原始数据 ----
  # 使用 stay_id + trial_period 作为 key 进行左连接
  # 这样每个 patient-period 的权重会被展开后的所有 day-level 行继承
  dt_sub <- data.table::as.data.table(dt_subset)
  # 将原始 period 列名重命名为 trial_period 以便 merge
  data.table::setnames(dt_sub, old = PERIOD_COL, new = "trial_period", skip_absent = TRUE)
  dt_sub <- merge(dt_sub, pre_cc[, .(stay_id, trial_period, iptw_trunc_pooled)],
                  by = c("stay_id", "trial_period"), all.x = TRUE)

  # ---- 第三步：日粒度展开（模拟 SEQexpand 逻辑） ----
  expanded <- data.table(
    stay_id = dt_sub[[ID_COL]],
    trial_period = dt_sub$trial_period,
    tx_init = as.integer(dt_sub[[TREATMENT_SRC_COL]] == TREATMENT_LABEL),
    outcome = as.integer(dt_sub[[OUTCOME_COL]]),
    surv_time = as.numeric(dt_sub[[SURV_TIME_COL]]),
    iptw_trunc_pooled = dt_sub$iptw_trunc_pooled,
    eligible = 1L
  )
  for (col in cov_cols) {
    data.table::set(expanded, j = col, value = dt_sub[[col]])
  }

  expanded[surv_time <= 0 & outcome == 0, eligible := 0L]
  expanded <- expanded[eligible == 1L & !is.na(iptw_trunc_pooled)]

  if (nrow(expanded) == 0) {
    return(list(method = "batch_seqexpand_pooled_ps", tag = out_tag,
                hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
                se = NA_real_, mean_iptw = mean_iptw, ps_mean = ps_mean,
                runtime_sec = NA_real_, peak_mem_mb = NA_real_, n_periods = 0L))
  }

  # ---- 日粒度展开 ----
  expand_seqexpand <- function(dt) {
    dt[, day_max := as.integer(floor(pmax(surv_time, 1)))]
    dt[day_max > 28L, day_max := 28L]
    result <- dt[rep(seq_len(.N), day_max + 1L)]
    result[, day := seq_len(.N) - 1L, by = .(stay_id, trial_period)]
    result[, y_daily := 0L]
    result[day == day_max & outcome == 1, y_daily := 1L]
    result[, c("surv_time", "day_max", "outcome") := NULL]
    result[]
  }

  person_dt <- expand_seqexpand(expanded)
  rm(expanded, dt_sub, pre_dt, pre_cc, ps_model); gc(verbose = FALSE)

  if (nrow(person_dt) < 100 || sum(person_dt$y_daily) < 10 ||
      length(unique(person_dt$tx_init)) < 2) {
    return(list(method = "batch_seqexpand_pooled_ps", tag = out_tag,
                hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
                se = NA_real_, mean_iptw = mean_iptw, ps_mean = ps_mean,
                runtime_sec = as.numeric(difftime(Sys.time(), t_start, units = "secs")),
                peak_mem_mb = get_peak_mem_mb(), n_periods = 0L))
  }

  # ---- 拟合 pooled logistic（与 SEQTaRget ITT 一致） ----
  person_dt[, `:=`(followup = as.numeric(day),
                   followup_sq = as.numeric(day^2),
                   trial = as.numeric(trial_period),
                   trial_sq = as.numeric(trial_period^2))]

  covs_use_all <- prepare_covariates(person_dt, PS_COVARIATES)
  cov_formula <- paste(covs_use_all, collapse = " + ")
  plogit_formula <- as.formula(
    sprintf("y_daily ~ tx_init + followup + followup_sq + trial + trial_sq + %s",
            cov_formula)
  )
  plogit_fit <- tryCatch(
    glm(plogit_formula, data = person_dt,
        family = binomial(), weights = iptw_trunc_pooled),
    error = function(e) NULL
  )

  if (is.null(plogit_fit)) {
    plogit_simple <- as.formula(sprintf("y_daily ~ tx_init + %s", cov_formula))
    plogit_fit <- tryCatch(
      glm(plogit_simple, data = person_dt,
          family = binomial(), weights = iptw_trunc_pooled),
      error = function(e) NULL
    )
  }

  if (is.null(plogit_fit)) {
    return(list(method = "batch_seqexpand_pooled_ps", tag = out_tag,
                hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
                se = NA_real_, mean_iptw = mean_iptw, ps_mean = ps_mean,
                runtime_sec = as.numeric(difftime(Sys.time(), t_start, units = "secs")),
                peak_mem_mb = get_peak_mem_mb(), n_periods = length(unique(person_dt$trial_period))))
  }

  # ---- cluster-robust SE ----
  coef_tx <- coef(plogit_fit)["tx_init"]
  se_tx <- sqrt(diag(vcovCL(plogit_fit, cluster = person_dt$stay_id)))["tx_init"]

  list(
    method = "batch_seqexpand_pooled_ps",
    tag = out_tag,
    hr = exp(coef_tx),
    hr_lower = exp(coef_tx - 1.96 * se_tx),
    hr_upper = exp(coef_tx + 1.96 * se_tx),
    se = se_tx,
    mean_iptw = mean_iptw,
    ps_mean = ps_mean,
    runtime_sec = as.numeric(difftime(Sys.time(), t_start, units = "secs")),
    peak_mem_mb = get_peak_mem_mb(),
    n_periods = length(unique(person_dt$trial_period))
  )
}

# ---- 运行官方 SEQTaRget ----
run_seqtarget_subset <- function(dt_subset) {
  t_start <- Sys.time()
  if (!requireNamespace("SEQTaRget", quietly = TRUE)) {
    stop("未安装 SEQTaRget 包")
  }

  dt <- data.table::as.data.table(dt_subset)
  dt[, id_int := as.integer(as.factor(get(ID_COL)))]
  dt[, period_int := as.integer(get(PERIOD_COL)) - 1L]
  dt[, y_int := as.integer(get(OUTCOME_COL))]
  dt[, a_int := as.integer(get(TREATMENT_SRC_COL) == TREATMENT_LABEL)]

  dt[, eligible_int := 1L]
  if ("baseline_phosphate" %in% names(dt)) {
    dt[baseline_phosphate > 0.80 | is.na(baseline_phosphate), eligible_int := 0L]
  }
  setorder(dt, id_int, period_int)
  dt[, cum_event := cumsum(y_int), by = id_int]
  dt[cum_event > 0 & y_int == 0 & eligible_int == 1, eligible_int := 0L]
  dt[, cum_event := NULL]

  # ---- 按 period 进一步剔除无治疗变异或事件不足的期 ----
  # SEQTaRget 内部按 period 拟合 PS，若某期全为治疗/对照或事件数为 0
  # 会导致 glm 无法找到有效起始值；提前标记为不合格可避免崩溃
  dt[, period_n := .N, by = period_int]
  dt[, period_tx_var := length(unique(a_int)) > 1, by = period_int]
  dt[, period_events := sum(y_int, na.rm = TRUE), by = period_int]
  dt[period_tx_var == FALSE | period_events < 5, eligible_int := 0L]
  dt[, c("period_n", "period_tx_var", "period_events") := NULL]

  covs_use <- intersect(COMPARE_COVARIATES, names(dt))
  for (col in covs_use) {
    if (is.character(dt[[col]]) || is.factor(dt[[col]])) {
      dt[[col]] <- as.numeric(factor(as.character(dt[[col]])))
    }
  }
  keep_cols <- c("id_int", "period_int", "a_int", "y_int", "eligible_int",
                 SURV_TIME_COL, "death_28day_num", "dischtime_from_t0",
                 covs_use)
  dt_sub <- dt[, ..keep_cols]
  dt_sub <- dt_sub[complete.cases(dt_sub)]

  if (nrow(dt_sub) == 0 || sum(dt_sub$y_int, na.rm = TRUE) < 5 ||
      length(unique(dt_sub$a_int)) < 2) {
    stop("SEQTaRget 子集：剔除不合格期后无足够变异或事件")
  }

  # ---- 创建 LTFU（失访）列与条件启用 cense ----
  # 原理解析：
  # 在 CCW 完全随访数据中，所有患者都有完整的结局记录（死亡或28天存活出院），
  # 因此实际上不存在传统意义上的失访（LTFU ≈ 0%）。
  #
  # SEQTaRget 内部通过 pooled logistic 模型拟合删失概率（IPCW），
  # 当 LTFU=0 或接近 0 时，存在完全分离问题，glm/fastglm
  # 无法找到有效起始值（"cannot find valid starting values"）。
  #
  # 修复策略：
  # 1. LTFU 仅标记幸存者的出院（不标记死亡），在方法学上更合理
  # 2. 判断 LTFU 率是否足够（≥0.5%）才启用 cense
  # 3. 若 LTFU 率不足或 cense 模型失败，自动回退到无 cense 模型

  # ---- 定义 LTFU ----
  # 标记规则：
  # - 死亡患者 (death_28day_num == 1)：不标记为 LTFU（死亡是结局事件）
  # - 幸存者 (death_28day_num == 0)：
  #   - 若出院时间 < 随访时间：标记为 LTFU=1（出院早于随访终点，信息性删失）
  #   - 否则：LTFU=0（随访完成）
  dt_sub[, ltfu := 0L]
  dt_sub[death_28day_num == 0 & !is.na(dischtime_from_t0) &
           dischtime_from_t0 < get(SURV_TIME_COL) - 0.5, ltfu := 1L]
  dt_sub[, eligible_cense := 1L]  # 所有患者均可被删失
  # 死亡患者不应该在删失模型中
  dt_sub[death_28day_num == 1, eligible_cense := 0L]

  ltfu_rate <- mean(dt_sub$ltfu, na.rm = TRUE)
  module_log_info(sprintf("  LTFU 率: %.2f%% (%d/%d)",
                          100 * ltfu_rate, sum(dt_sub$ltfu, na.rm = TRUE), nrow(dt_sub)))

  # 仅在 LTFU 足够时启用 cense（≥0.5%），否则保持关闭
  CENSE_ENABLE <- ltfu_rate >= 0.005

  time_varying <- intersect(COVARIATE_TIME_VARYING, names(dt_sub))
  fixed_covs <- setdiff(covs_use, time_varying)
  fixed_covs <- setdiff(fixed_covs, c("ltfu", "eligible_cense",
                                       "death_28day_num", "dischtime_from_t0"))

  # ---- 构建 SEQopts（条件启用 cense） ----
  # 若 LTFU 不足则设置 cense=NA 避免模型崩溃
  opts <- SEQTaRget::SEQopts(
    bootstrap = COMPARE_BOOTSTRAP > 0L,
    bootstrap.nboot = COMPARE_BOOTSTRAP,
    bootstrap.sample = 1.0,
    weighted = TRUE,
    weight.p99 = TRUE,
    hazard = TRUE,
    parallel = COMPARE_NTHREADS > 1L,
    ncores = COMPARE_NTHREADS,
    seed = RANDOM_SEED,
    fastglm.method = COMPARE_FASTGLM_METHOD,
    followup.include = COMPARE_FOLLOWUP_INCLUDE,
    trial.include = COMPARE_TRIAL_INCLUDE,
    cense = if (CENSE_ENABLE) "ltfu" else NA_character_,
    cense.eligible = if (CENSE_ENABLE) "eligible_cense" else NA_character_
  )

  if (CENSE_ENABLE) {
    module_log_info("  条件启用 cense（LTFU 率 ≥ 0.5%），IPCW 已激活")
  } else {
    module_log_info("  LTFU 率不足，cense 已关闭（IPCW 未激活）")
  }

  # ---- 多層自动回退拟合 ----
  # 原因分析：
  # a) fastglm.method=1L (unpivoted QR) 在 bootstrap 样本中遇到秩亏设计矩阵
  #    时崩溃，3L (LDLT Cholesky) 数值更稳定。
  # b) cense（IPCW）激活后，若 LTFU 率不足以支持删失模型拟合，
  #    fastglm/glm 因完全分离无法找到有效起始值
  #    （"cannot find valid starting values"）。
  #
  # 回退优先级：1L+boot → 3L+boot(串行) → 关闭cense+1L+boot
  #           → 关闭cense+3L+boot → 关闭cense+1L+boot
  #           → 关闭bootstrap(仅点估计)
  run_seq_with_fallback <- function(opts_input, label = "") {
    result <- tryCatch({
      SEQTaRget::SEQuential(
        data = as.data.frame(dt_sub),
        id.col = "id_int",
        time.col = "period_int",
        eligible.col = "eligible_int",
        treatment.col = "a_int",
        outcome.col = "y_int",
        fixed.cols = fixed_covs,
        time_varying.cols = time_varying,
        method = "ITT",
        options = opts_input,
        verbose = FALSE
      )
    }, error = function(e) {
      module_log_warn(sprintf("SEQTaRget%s失败: %s", label, e$message))
      NULL
    })
    result
  }

  has_ci <- function(res) {
    if (is.null(res)) return(FALSE)
    hr_vec <- tryCatch(SEQTaRget::hazard_ratio(res), error = function(e) NULL)
    if (is.null(hr_vec)) return(FALSE)
    has_lci <- "LCI" %in% names(hr_vec) && is.finite(hr_vec["LCI"])
    has_uci <- "UCI" %in% names(hr_vec) && is.finite(hr_vec["UCI"])
    has_lci && has_uci
  }

  # ---- 第1层：用户指定 method + bootstrap ----
  result <- run_seq_with_fallback(opts, "")

  # ---- 第2层：若失败，尝试 3L+bootstrap+关闭并行 ----
  if (!has_ci(result) && COMPARE_BOOTSTRAP > 0L) {
    if (is.null(result)) {
      module_log_warn("SEQTaRget → 回退[2]: 3L + bootstrap（关闭并行避免 globals 限制）")
    } else {
      module_log_info("  SEQTaRget 1L 无 CI → 尝试 3L+bootstrap（串行）")
    }
    opts_3l_ser <- opts
    opts_3l_ser@fastglm.method <- 3L
    opts_3l_ser@bootstrap <- TRUE
    opts_3l_ser@bootstrap.nboot <- COMPARE_BOOTSTRAP
    opts_3l_ser@parallel <- FALSE   # 关闭并行，避免 future globals 限制
    opts_3l_ser@ncores <- 1L
    result_3l <- run_seq_with_fallback(opts_3l_ser, " (3L-boot-serial)")
    if (!is.null(result_3l)) result <- result_3l
  }

  # ---- 第3层：若 cense 已启用但失败，尝试关闭 cense + 1L+boot ----
  # 适用于 LTFU=0 导致的 "cannot find valid starting values"
  if (!has_ci(result) && CENSE_ENABLE && COMPARE_BOOTSTRAP > 0L) {
    module_log_warn("SEQTaRget → 回退[3]: 关闭 cense（LTFU 不足导致完全分离）+ 1L+boot")
    opts_no_cense <- opts
    opts_no_cense@cense <- NA_character_
    opts_no_cense@cense.eligible <- NA_character_
    result_nc <- run_seq_with_fallback(opts_no_cense, " (no-cense)")
    if (!is.null(result_nc)) result <- result_nc
  }

  # ---- 第4层：若 cense 已启用但失败，尝试关闭 cense + 3L+boot ----
  if (!has_ci(result) && CENSE_ENABLE && COMPARE_BOOTSTRAP > 0L) {
    module_log_warn("SEQTaRget → 回退[4]: 关闭 cense + 3L+boot（串行）")
    opts_nc_3l <- opts
    opts_nc_3l@cense <- NA_character_
    opts_nc_3l@cense.eligible <- NA_character_
    opts_nc_3l@fastglm.method <- 3L
    opts_nc_3l@bootstrap <- TRUE
    opts_nc_3l@bootstrap.nboot <- COMPARE_BOOTSTRAP
    opts_nc_3l@parallel <- FALSE
    opts_nc_3l@ncores <- 1L
    result_nc_3l <- run_seq_with_fallback(opts_nc_3l, " (no-cense-3L)")
    if (!is.null(result_nc_3l)) result <- result_nc_3l
  }

  # ---- 最终回退：关闭 bootstrap 获取点估计 ----
  if (!has_ci(result)) {
    module_log_warn("SEQTaRget → 最终回退: 关闭 bootstrap 仅获取点估计")
    opts_no_boot <- opts
    opts_no_boot@bootstrap <- FALSE
    opts_no_boot@bootstrap.nboot <- 0L
    result_nb <- run_seq_with_fallback(opts_no_boot, " (no-boot)")
    if (!is.null(result_nb)) result <- result_nb
  }

  if (is.null(result)) {
    stop("SEQTaRget 最终模型结果为空（所有回退层均失败）")
  }

  hr_vec <- SEQTaRget::hazard_ratio(result)
  hr <- hr_vec[which(is.finite(hr_vec))[1]]
  if (is.na(hr)) hr <- NA_real_

  ci_lower <- NA_real_
  ci_upper <- NA_real_
  if ("LCI" %in% names(hr_vec) && is.finite(hr_vec["LCI"])) {
    ci_lower <- hr_vec["LCI"]
  }
  if ("UCI" %in% names(hr_vec) && is.finite(hr_vec["UCI"])) {
    ci_upper <- hr_vec["UCI"]
  }
  # 若 LCI/UCI 为空，尝试 stats::confint
  if (is.na(ci_lower) || is.na(ci_upper)) {
    ci <- tryCatch(stats::confint(result), error = function(e) NULL)
    if (!is.null(ci) && length(ci) >= 2) {
      if (is.matrix(ci)) {
        ci_row <- ci[which(is.finite(ci[, 1])), , drop = FALSE]
        if (nrow(ci_row) > 0) {
          ci_lower <- ci_row[1, 1]
          ci_upper <- ci_row[1, 2]
        }
      } else {
        ci_lower <- ci[1]
        ci_upper <- ci[2]
      }
    }
  }

  weights <- tryCatch(result@weights, error = function(e) NULL)
  mean_iptw <- NA_real_
  ps_mean <- NA_real_
  if (!is.null(weights) && length(weights) > 0) {
    mean_iptw <- mean(as.numeric(unlist(weights)), na.rm = TRUE)
  }

  # 导出使用的最终 fastglm.method（用于诊断）
  final_method <- tryCatch(result@params@fastglm.method, error = function(e) NA_integer_)
  used_cense <- tryCatch(!is.na(result@params@cense), error = function(e) FALSE)

  list(
    method = "SEQTaRget",
    tag = sprintf("SEQTaRget(fgm=%d,cense=%s)", final_method, used_cense),
    hr = hr,
    hr_lower = ci_lower,
    hr_upper = ci_upper,
    se = NA_real_,
    mean_iptw = mean_iptw,
    ps_mean = ps_mean,
    runtime_sec = as.numeric(difftime(Sys.time(), t_start, units = "secs")),
    peak_mem_mb = get_peak_mem_mb(),
    n_periods = length(unique(dt_sub$period_int))
  )
}

# ---- 主流程 ----
module_log_section(sprintf("步骤四：与 SEQTaRget 对比（fastglm.method=%d, followup.include=%s, trial.include=%s）",
                           COMPARE_FASTGLM_METHOD, COMPARE_FOLLOWUP_INCLUDE, COMPARE_TRIAL_INCLUDE))

set.seed(RANDOM_SEED)
dt_full <- as.data.table(arrow::read_parquet(INPUT_PARQUET_AUGMENTED))
all_ids <- unique(dt_full[[ID_COL]])
module_log_info(sprintf("全量数据: %d 例患者，%d 行记录", length(all_ids), nrow(dt_full)))

comparison_rows <- list()

for (n_sub in COMPARE_SAMPLE_SIZES) {
  n_sub <- min(n_sub, length(all_ids))
  sampled_ids <- sample(all_ids, n_sub)
  dt_sub <- dt_full[get(ID_COL) %in% sampled_ids]

  tmp_path <- file.path(TEMP_DIR, sprintf("tmp_seqtarget_batch_compare_n%d.parquet", n_sub))
  arrow::write_parquet(dt_sub, tmp_path)
  module_log_info(sprintf("子集 n=%d 已保存至: %s", n_sub, tmp_path))

  # 诊断：确认 COMPARE_SKIP_SEQTARGET 值
  module_log_info(sprintf("COMPARE_SKIP_SEQTARGET=%s (class=%s)",
                          COMPARE_SKIP_SEQTARGET, class(COMPARE_SKIP_SEQTARGET)))

  if (!COMPARE_SKIP_SEQTARGET) {
    # SEQTaRget
    seq_out <- tryCatch(run_seqtarget_subset(dt_sub), error = function(e) {
      module_log_warn(sprintf("SEQTaRget n=%d 失败: %s", n_sub, e$message))
      list(method = "SEQTaRget", tag = sprintf("n%d", n_sub),
           hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
           se = NA_real_, mean_iptw = NA_real_, ps_mean = NA_real_,
           runtime_sec = NA_real_, peak_mem_mb = NA_real_, n_periods = NA_integer_)
    })
    seq_out$tag <- sprintf("n%d", n_sub)
    comparison_rows[[paste0("seq_", n_sub)]] <- as.data.table(seq_out)
  }

  # Batch engine
  batch_out <- tryCatch(run_batch_engine(dt_sub, sprintf("n%d", n_sub)), error = function(e) {
    module_log_warn(sprintf("Batch engine n=%d 失败: %s", n_sub, e$message))
    list(method = "batch_engine", tag = sprintf("n%d", n_sub),
         hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
         se = NA_real_, mean_iptw = NA_real_, ps_mean = NA_real_,
         runtime_sec = NA_real_, peak_mem_mb = NA_real_, n_periods = NA_integer_)
  })
  comparison_rows[[paste0("batch_", n_sub)]] <- as.data.table(batch_out)

  # Batch engine SEQexpand（日粒度展开 + pooled logistic，用于方法学对比）
  if (n_sub == 5000L) {
    module_log_info("--- 运行 SEQTaRget 展开式 Batch Engine（n=5000） ---")
    batch_seq_out <- tryCatch(run_batch_engine_seqexpand(dt_sub, sprintf("n%d", n_sub)),
      error = function(e) {
        module_log_warn(sprintf("Batch seqexpand n=%d 失败: %s", n_sub, e$message))
        list(method = "batch_seqexpand", tag = sprintf("n%d", n_sub),
             hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
             se = NA_real_, mean_iptw = NA_real_, ps_mean = NA_real_,
             runtime_sec = NA_real_, peak_mem_mb = NA_real_, n_periods = NA_integer_)
      })
    comparison_rows[[paste0("batch_seq_", n_sub)]] <- as.data.table(batch_seq_out)
  }

  # Batch engine SEQexpand Pooled PS（展开前 pooled PS + 日粒度展开，与 SEQTaRget 对齐）
  if (n_sub == 5000L) {
    module_log_info("--- 运行 Pooled PS 版 Batch Seqexpand（n=5000） ---")
    batch_seq_pooled_out <- tryCatch(run_batch_engine_seqexpand_pooled(dt_sub, sprintf("n%d", n_sub)),
      error = function(e) {
        module_log_warn(sprintf("Batch seqexpand pooled PS n=%d 失败: %s", n_sub, e$message))
        list(method = "batch_seqexpand_pooled_ps", tag = sprintf("n%d", n_sub),
             hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
             se = NA_real_, mean_iptw = NA_real_, ps_mean = NA_real_,
             runtime_sec = NA_real_, peak_mem_mb = NA_real_, n_periods = NA_integer_)
      })
    comparison_rows[[paste0("batch_seq_pooled_", n_sub)]] <- as.data.table(batch_seq_pooled_out)
  }

  rm(dt_sub); gc(verbose = FALSE)
}

# ---- 全量数据对比 ----
if (COMPARE_FULL_DATA) {
  module_log_section("全量数据对比（所有患者）")
  n_full <- length(all_ids)
  module_log_info(sprintf("全量数据: n=%d 例患者", n_full))

  if (!COMPARE_SKIP_SEQTARGET) {
    # SEQTaRget 全量
    seq_out_full <- tryCatch(run_seqtarget_subset(dt_full), error = function(e) {
      module_log_warn(sprintf("SEQTaRget 全量失败: %s", e$message))
      list(method = "SEQTaRget", tag = "full",
           hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
           se = NA_real_, mean_iptw = NA_real_, ps_mean = NA_real_,
           runtime_sec = NA_real_, peak_mem_mb = NA_real_, n_periods = NA_integer_)
    })
    seq_out_full$tag <- "full"
    comparison_rows[["seq_full"]] <- as.data.table(seq_out_full)
  }

  # Batch engine 全量
  batch_out_full <- tryCatch(run_batch_engine(dt_full, "full"), error = function(e) {
    module_log_warn(sprintf("Batch engine 全量失败: %s", e$message))
    list(method = "batch_engine", tag = "full",
         hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
         se = NA_real_, mean_iptw = NA_real_, ps_mean = NA_real_,
         runtime_sec = NA_real_, peak_mem_mb = NA_real_, n_periods = NA_integer_)
  })
  comparison_rows[["batch_full"]] <- as.data.table(batch_out_full)

  # Batch engine SEQexpand Pooled PS 全量（展开前 pooled PS + 日粒度展开）
  module_log_info("--- 运行 Pooled PS 版 Batch Seqexpand（全量） ---")
  batch_seq_pooled_full <- tryCatch(run_batch_engine_seqexpand_pooled(dt_full, "full"),
    error = function(e) {
      module_log_warn(sprintf("Batch seqexpand pooled PS 全量失败: %s", e$message))
      list(method = "batch_seqexpand_pooled_ps", tag = "full",
           hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
           se = NA_real_, mean_iptw = NA_real_, ps_mean = NA_real_,
           runtime_sec = NA_real_, peak_mem_mb = NA_real_, n_periods = NA_integer_)
    })
  comparison_rows[["batch_seq_pooled_full"]] <- as.data.table(batch_seq_pooled_full)
}

comparison_dt <- rbindlist(comparison_rows, use.names = TRUE, fill = TRUE)
comparison_csv <- file.path(RESULTS_DIR, "seqtarget_vs_batch_engine_comparison.csv")
fwrite(comparison_dt, comparison_csv)

module_log_info(sprintf("对比结果已保存: %s", comparison_csv))
print(comparison_dt)

# ---- 可视化 ----
if (requireNamespace("ggplot2", quietly = TRUE) && nrow(comparison_dt) > 0) {
  plot_dt <- comparison_dt[!is.na(hr)]
  if (nrow(plot_dt) > 0) {
    p <- ggplot(plot_dt, aes(x = tag, y = hr, color = method)) +
      geom_point(position = position_dodge(width = 0.4), size = 3) +
      geom_errorbar(aes(ymin = hr_lower, ymax = hr_upper),
                    width = 0.2, position = position_dodge(width = 0.4)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
      scale_color_manual(values = c("SEQTaRget" = "#BC3C29",
                                     "batch_engine" = "#0072B5",
                                     "batch_seqexpand" = "#E69F00",
                                     "batch_seqexpand_pooled_ps" = "#2E8B57")) +
      labs(title = sprintf("SEQTaRget vs Batch Engine HR (fastglm.method=%d, followup.include=%s, trial.include=%s)",
                           COMPARE_FASTGLM_METHOD, COMPARE_FOLLOWUP_INCLUDE, COMPARE_TRIAL_INCLUDE),
           subtitle = sprintf("Covariates: %s", paste(PS_COVARIATES, collapse = ", ")),
           x = "Sample size", y = "Hazard Ratio", color = "Method") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")

    ggsave(file.path(RESULTS_DIR, "FigS_batch_vs_seqtarget_hr.pdf"), p, width = 7, height = 5)
    ggsave(file.path(RESULTS_DIR, "FigS_batch_vs_seqtarget_hr.png"), p, width = 7, height = 5, dpi = 300)
    module_log_info("对比图已保存: FigS_batch_vs_seqtarget_hr.pdf/png")
  }
}

module_log_section("对比步骤完成")
rm(dt_full, comparison_dt, comparison_rows); gc(verbose = FALSE)
