# ============================================================================
# V26 统一辅助函数集
# 所有分析脚本通过 source("_helpers_v26.R") 加载
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
})

# ---- 日志工具 ----
log_section <- function(msg) cat(sprintf("\n===== %s =====\n", msg))
log_info <- function(msg) cat(sprintf("  [INFO] %s\n", msg))
log_warn <- function(msg) cat(sprintf("  [WARN] %s\n", msg))
log_key <- function(lbl, val) cat(sprintf("  [KEY] %s: %s\n", lbl, val))

# ---- 数据加载 ----
load_data <- function(path) {
  if (grepl("\\.parquet$", path)) {
    if (!requireNamespace("arrow", quietly = TRUE)) stop("arrow package required")
    as.data.table(arrow::read_parquet(path))
  } else if (grepl("\\.csv$", path)) {
    fread(path)
  } else if (grepl("\\.rds$", path)) {
    readRDS(path)
  } else stop("Unsupported file format: ", path)
}

# ---- SMD 计算 ----
compute_smd <- function(x_treat, x_ctrl, w_treat = NULL, w_ctrl = NULL) {
  if (is.null(w_treat)) w_treat <- rep(1, length(x_treat))
  if (is.null(w_ctrl)) w_ctrl <- rep(1, length(x_ctrl))
  m1 <- weighted.mean(x_treat, w_treat, na.rm = TRUE)
  m0 <- weighted.mean(x_ctrl, w_ctrl, na.rm = TRUE)
  v1 <- sum(w_treat * (x_treat - m1)^2, na.rm = TRUE) / sum(w_treat, na.rm = TRUE)
  v0 <- sum(w_ctrl * (x_ctrl - m0)^2, na.rm = TRUE) / sum(w_ctrl, na.rm = TRUE)
  sd_pool <- sqrt((v1 + v0) / 2)
  if (sd_pool < 1e-8) 0 else abs(m1 - m0) / sd_pool
}

# ---- calc_smd 兼容包装（管道使用的旧版接口） ----
# x: 变量值向量, treat: 治疗指示向量, w: IPTW权重(NUll=无权重), var_name: 变量名(日志用)
calc_smd <- function(x, treat, w = NULL, var_name = "") {
  x_treat <- x[treat == 1]; x_ctrl <- x[treat == 0]
  if (!is.null(w)) {
    w_treat <- w[treat == 1]; w_ctrl <- w[treat == 0]
  } else {
    w_treat <- NULL; w_ctrl <- NULL
  }
  compute_smd(x_treat, x_ctrl, w_treat, w_ctrl)
}

# ---- E-value 计算 ----
compute_evalue <- function(hr, ci_lower, ci_upper) {
  if (is.na(hr) || is.na(ci_lower) || is.na(ci_upper)) {
    return(list(e_value_point = NA, e_value_ci = NA, grade = "无法计算"))
  }
  # E-value for point estimate
  if (hr <= 1) {
    rr <- 1 / hr
    e_value_point <- rr + sqrt(rr * (rr - 1))
  } else {
    e_value_point <- hr + sqrt(hr * (hr - 1))
  }

  # E-value for CI lower bound
  if (ci_lower <= 1 && ci_upper <= 1) {
    rr_ci <- 1 / ci_upper
    e_value_ci <- rr_ci + sqrt(rr_ci * (rr_ci - 1))
  } else if (ci_lower >= 1 && ci_upper >= 1) {
    e_value_ci <- ci_lower + sqrt(ci_lower * (ci_lower - 1))
  } else {
    e_value_ci <- 1.0
  }

  # Grade
  if (is.na(e_value_point)) {
    grade <- "无法计算"
  } else if (e_value_point >= 2.0) {
    grade <- "强稳健性"
  } else if (e_value_point >= 1.5) {
    grade <- "中等稳健性"
  } else {
    grade <- "弱稳健性"
  }

  list(e_value_point = round(e_value_point, 2),
       e_value_ci = round(e_value_ci, 2),
       grade = grade)
}

# ---- IPTW 权重 ----
compute_iptw <- function(data, covariates, treat_col = "treatment",
                          trunc_p = 0.99, lower_p = 0.01,
                          stabilize = TRUE, bilateral = TRUE) {
  formula_str <- paste(treat_col, "~", paste(covariates, collapse = " + "))
  ps_model <- glm(as.formula(formula_str), data = data, family = binomial())
  ps <- predict(ps_model, type = "response")
  ps <- pmax(pmin(ps, 0.999), 0.001)
  tx <- data[[treat_col]]
  p_tx <- mean(tx)
  w <- if (stabilize) ifelse(tx == 1, p_tx / ps, (1 - p_tx) / (1 - ps))
       else ifelse(tx == 1, 1 / ps, 1 / (1 - ps))
  if (bilateral) {
    w99 <- quantile(w, trunc_p, na.rm = TRUE)
    w01 <- quantile(w, lower_p, na.rm = TRUE)
    w <- pmax(pmin(w, w99), w01)
  }
  w[is.na(w) | is.infinite(w)] <- 1.0
  w
}

# ---- 协变量获取 ----
get_covariates <- function(data, pool) intersect(pool, names(data))

# ---- 提取 Cox HR ----
extract_cox_hr <- function(fit) {
  s <- summary(fit)
  ci <- exp(confint(fit))
  data.frame(hr = exp(s$coefficients[1,1]), ci_lower = ci[1,1],
             ci_upper = ci[1,2], p_value = s$coefficients[1, ncol(s$coefficients)],
             hr_formatted = sprintf("%.3f [%.3f-%.3f]", exp(s$coefficients[1,1]), ci[1,1], ci[1,2]))
}

# ---- Rubin 规则 ----
rubin_pool <- function(estimates, ses) {
  m <- length(estimates)
  pooled_est <- mean(estimates)
  within_var <- mean(ses^2)
  between_var <- var(estimates)
  total_var <- within_var + (1 + 1/m) * between_var
  pooled_se <- sqrt(total_var)
  list(estimate = pooled_est, se = pooled_se,
       ci_lower = pooled_est - 1.96 * pooled_se,
       ci_upper = pooled_est + 1.96 * pooled_se)
}

rubin_pool_hr <- function(log_hrs, log_ses) {
  pooled <- rubin_pool(log_hrs, log_ses)
  list(hr = exp(pooled$estimate), ci_lower = exp(pooled$ci_lower),
       ci_upper = exp(pooled$ci_upper))
}

# ---- 权重诊断 ----
weight_diagnostics <- function(weights, label = "") {
  w <- weights[!is.na(weights) & is.finite(weights)]
  ess <- sum(w)^2 / sum(w^2)
  list(label = label, n = length(w), mean = mean(w), sd = sd(w),
       min = min(w), p01 = quantile(w, 0.01), p99 = quantile(w, 0.99),
       max = max(w), ess = ess, ess_ratio = ess / length(w),
       extreme_pct = mean(w > quantile(w, 0.99) | w < quantile(w, 0.01)) * 100)
}

# ---- 格式化 ----
format_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
format_hr <- function(hr, ci_l, ci_u) sprintf("%.3f (%.3f-%.3f)", hr, ci_l, ci_u)

# ---- JAMA 配色 ----
JAMA_RED  <- "#BC3C29"
JAMA_BLUE <- "#0072B5"
JAMA_GREY <- "#BDC3C7"

cat("_helpers_v26.R loaded\n")
