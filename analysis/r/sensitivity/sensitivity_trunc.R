# ============================================================================
# V26.9 权重截断敏感性分析 — Manual IPTW-Cox (SEQTaRget 框架内)
#
# 目的：在 SEQTaRget 框架下评估不同权重截断水平对主分析结果的影响
#   截断水平：99%/97.5%/95%/90% (双侧对称截断)
#   
# 方法：Manual IPTW-Cox × 5 imp, DAG14 协变量, 4 个截断水平
#   注：SEQTaRget::SEQuential() 内部处理 IPTW 权重截断，不暴露截断参数
#   此处采用与 SEQTaRget 相同的数据预处理和协变量集，使用手动 IPTW-Cox
#   以评估截断水平的敏感性
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(survival)
})

# ============================================================================
# 1. 参数配置
# ============================================================================
OUT_DIR  <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9"
RDS_DIR  <- file.path(OUT_DIR, "rds")
TBL_DIR  <- file.path(OUT_DIR, "tables")
LOG_DIR  <- file.path(OUT_DIR, "logs")
for (d in c(OUT_DIR, RDS_DIR, TBL_DIR, LOG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

DATA_FILE           <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_imputed_v26.parquet"
MICE_PARQUET_TEMPL  <- "/data/小米云盘/文档/论文/血磷/sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet"

N_IMPUTE     <- 5L        # 前 5 个插补集
SEED         <- 20260705L
TRUNC_LEVELS <- c(0.99, 0.975, 0.95, 0.90)  # 双侧对称截断
set.seed(SEED)
options(mc.cores = 1)

DAG14 <- c(
  "age", "gender", "baseline_sofa", "sepsis3", "charlson",
  "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
  "admission_type", "baseline_phosphate", "lactate_t0", "refeeding_risk",
  "corrected_calcium_mmol_min", "nutrition_support_type"
)

# ============================================================================
# 2. 日志系统
# ============================================================================
log_file <- file.path(LOG_DIR, sprintf("v26_9_sensitivity_trunc_%s.log",
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

log_info("V26.9 权重截断敏感性分析 — Manual IPTW-Cox")
log_info(sprintf("种子=%d, 截断水平: %s", SEED, paste(TRUNC_LEVELS, collapse = ",")))
log_info(sprintf("插补集: 1-%d", N_IMPUTE))

# ============================================================================
# 3. 数据加载
# ============================================================================
log_section("数据加载")
dt_full <- as.data.table(read_parquet(DATA_FILE))
dt_full[, `:=`(
  tx_int = as.integer(itt_group == "treatment"),
  id_int = as.integer(factor(stay_id))
)]
setorder(dt_full, id_int, trial_num)

key_cols <- c("stay_id", "subject_id", "trial_num", "time_zero",
              "itt_group", "surv_time", "death_28day_num")
dt_key <- dt_full[, intersect(key_cols, names(dt_full)), with = FALSE]
log_info(sprintf("主数据 key 表: %d 行, %d 列", nrow(dt_key), ncol(dt_key)))
rm(dt_full); gc(TRUE)

# ============================================================================
# 4. 辅助函数：Manual IPTW-Cox
# ============================================================================
log_section("辅助函数定义")

run_manual_iptw_cox <- function(dt_comb, dag14_covs, trunc_upper, trunc_lower = NULL) {
  dt_sub <- copy(dt_comb)
  if (is.null(trunc_lower)) trunc_lower <- 1 - trunc_upper

  dt_sub[, `:=`(
    id_int       = as.integer(factor(stay_id)),
    period_int   = as.integer(trial_num) - 1L,
    a_int        = as.integer(itt_group == "treatment"),
    y_int        = as.integer(death_28day_num)
  )]
  setorder(dt_sub, id_int, period_int)

  # 按 period 过滤
  dt_sub[, `:=`(n_a = sum(a_int), n_y = sum(y_int)), by = period_int]
  dt_sub <- dt_sub[n_a > 0 & n_a < .N & n_y >= 5]
  if (nrow(dt_sub) < 100) return(NULL)

  # 协变量准备
  covs <- intersect(dag14_covs, names(dt_sub))
  covs <- covs[sapply(covs, function(v) {
    x <- dt_sub[[v]]
    if (is.factor(x))    return(nlevels(x) > 1)
    if (is.character(x)) return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]
  if (length(covs) < 2) return(NULL)
  for (v in covs) {
    if (is.character(dt_sub[[v]])) dt_sub[[v]] <- as.factor(dt_sub[[v]])
  }

  # 逐期 PS 模型 + IPTW
  all_periods <- sort(unique(dt_sub$period_int))
  weighted_list <- list()

  for (p in all_periods) {
    period_data <- dt_sub[period_int == p]
    if (nrow(period_data) < 20) next

    tx <- period_data$a_int
    n_tx <- sum(tx); n_ct <- sum(1 - tx)
    if (n_tx < 2 || n_ct < 2) next

    # PS 模型
    ps_formula <- as.formula(paste("a_int ~", paste(covs, collapse = " + ")))
    ps_fit <- tryCatch(glm(ps_formula, data = period_data, family = binomial), error = function(e) NULL)
    if (is.null(ps_fit)) next

    ps <- predict(ps_fit, type = "response")
    # 稳定化 IPTW
    p_tx <- mean(tx)
    iptw <- ifelse(tx == 1, p_tx / ps, (1 - p_tx) / (1 - ps))

    # 双侧截断
    q_upper <- quantile(iptw, trunc_upper, na.rm = TRUE)
    q_lower <- quantile(iptw, trunc_lower, na.rm = TRUE)
    iptw_trunc <- pmax(pmin(iptw, q_upper), q_lower)

    period_data[, `:=`(iptw_trunc = iptw_trunc, period_ps = ps)]
    weighted_list[[length(weighted_list) + 1]] <- period_data
  }

  if (length(weighted_list) == 0) return(NULL)
  dt_w <- rbindlist(weighted_list)

  # 加权 Cox 回归
  cox_formula <- as.formula(paste("Surv(surv_time, y_int) ~ a_int +", paste(covs, collapse = " + ")))
  cox_fit <- tryCatch(
    coxph(cox_formula, data = dt_w, weights = iptw_trunc, cluster = stay_id),
    error = function(e) NULL
  )
  if (is.null(cox_fit)) return(NULL)

  cox_sum <- summary(cox_fit)
  coef_idx <- which(names(coef(cox_fit)) == "a_int")
  if (length(coef_idx) == 0) return(NULL)

  hr <- exp(coef(cox_fit)[coef_idx])
  ci <- tryCatch(exp(confint(cox_fit)[coef_idx, ]), error = function(e) c(NA, NA))
  p_val <- tryCatch(cox_sum$coefficients[coef_idx, "Pr(>|z|)"], error = function(e) NA)

  return(list(
    hr = unname(hr),
    ci_lower = unname(ci[1]),
    ci_upper = unname(ci[2]),
    p = unname(p_val),
    n = nrow(dt_w),
    n_events = sum(dt_w$y_int == 1, na.rm = TRUE),
    n_periods = length(weighted_list),
    mean_iptw = mean(dt_w$iptw_trunc, na.rm = TRUE),
    ess = sum(dt_w$iptw_trunc)^2 / sum(dt_w$iptw_trunc^2)
  ))
}

rubin_combine <- function(betas, ses) {
  m  <- length(betas)
  bm <- mean(betas)
  W  <- mean(ses^2)
  B  <- if (m > 1) sum((betas - bm)^2) / (m - 1) else 0
  se_pool <- sqrt(W + (1 + 1 / m) * B)
  t_stat  <- bm / se_pool
  df <- if (m > 1) (m - 1) * (1 + W / ((1 + 1 / m) * B))^2 else 1
  p  <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)
  list(beta = bm, se = se_pool, t = t_stat, df = df, p = p,
       hr = exp(bm),
       ci_lower = exp(bm - 1.96 * se_pool),
       ci_upper = exp(bm + 1.96 * se_pool),
       W = W, B = B)
}

log_info("辅助函数定义完成")

# ============================================================================
# 5. 核心分析：Manual IPTW-Cox × 5 imp × 4 truncation levels
# ============================================================================
log_section("Manual IPTW-Cox × 5 插补 × 4 截断水平")

all_results <- list()

for (trunc in TRUNC_LEVELS) {
  trunc_label <- sprintf("%.0f%%", trunc * 100)
  log_info(sprintf("\n======== 截断水平: %s ========", trunc_label))

  imp_betas  <- numeric()
  imp_ses    <- numeric()
  imp_hrs    <- numeric()
  imp_ci_l   <- numeric()
  imp_ci_u   <- numeric()
  imp_p      <- numeric()
  imp_success <- integer()
  imp_fail   <- integer()

  for (i_imp in 1L:N_IMPUTE) {
    t_start <- Sys.time()
    log_info(sprintf("--- [%s] 插补集 %02d/%02d ---", trunc_label, i_imp, N_IMPUTE))

    parquet_file <- sprintf(MICE_PARQUET_TEMPL, i_imp)
    if (!file.exists(parquet_file)) {
      imp_fail <- c(imp_fail, i_imp)
      next
    }
    dt_imp <- as.data.table(read_parquet(parquet_file))
    keep_cols <- intersect(c("stay_id", "trial_num", DAG14), names(dt_imp))
    dt_imp <- dt_imp[, ..keep_cols]
    gc(TRUE)

    covs_in_imp <- intersect(DAG14, names(dt_imp))
    covs_remove <- intersect(covs_in_imp, names(dt_key))
    by_cols <- c("stay_id", "trial_num")
    dt_comb <- merge(
      dt_key[, setdiff(names(dt_key), covs_remove), with = FALSE],
      dt_imp,
      by = by_cols, all.x = TRUE, sort = FALSE
    )
    for (v in intersect(names(dt_comb), names(dt_imp))) {
      if (is.character(dt_comb[[v]])) dt_comb[[v]] <- as.factor(dt_comb[[v]])
    }

    r <- tryCatch(
      run_manual_iptw_cox(dt_comb, DAG14, trunc_upper = trunc),
      error = function(e) NULL
    )

    t_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    if (is.null(r)) {
      log_info(sprintf("  [%s imp%02d] 失败, %.1f min", trunc_label, i_imp, t_elapsed / 60))
      imp_fail <- c(imp_fail, i_imp)
      next
    }

    imp_success <- c(imp_success, i_imp)
    imp_betas   <- c(imp_betas, log(r$hr))
    imp_ses     <- c(imp_ses, (log(r$ci_upper) - log(r$ci_lower)) / (2 * 1.96))
    imp_hrs     <- c(imp_hrs, r$hr)
    imp_ci_l    <- c(imp_ci_l, r$ci_lower)
    imp_ci_u    <- c(imp_ci_u, r$ci_upper)
    imp_p       <- c(imp_p, r$p)

    log_info(sprintf("  [%s imp%02d] HR=%.4f (%.4f-%.4f), p=%.4f, ESS=%.0f, runtime=%.0fs",
                     trunc_label, i_imp, r$hr, r$ci_lower, r$ci_upper, r$p, r$ess, r$runtime))
  }

  # Rubin 合并
  N_success <- length(imp_success)
  log_info(sprintf("  [%s] 成功: %d/%d", trunc_label, N_success, N_IMPUTE))

  if (N_success >= 2) {
    rubin <- rubin_combine(imp_betas, imp_ses)
    all_results[[trunc_label]] <- list(
      hr = rubin$hr, ci_lower = rubin$ci_lower, ci_upper = rubin$ci_upper,
      p = rubin$p, W = rubin$W, B = rubin$B,
      n_success = N_success
    )
    log_info(sprintf("  [%s] Rubin HR=%.4f (%.4f-%.4f), p=%.4f",
                     trunc_label, rubin$hr, rubin$ci_lower, rubin$ci_upper, rubin$p))
    log_info(sprintf("  [%s] beta=%.6f, SE=%.6f, W=%.6f, B=%.6f",
                     trunc_label, rubin$beta, rubin$se, rubin$W, rubin$B))
  } else if (N_success == 1) {
    all_results[[trunc_label]] <- list(
      hr = imp_hrs[1], ci_lower = imp_ci_l[1], ci_upper = imp_ci_u[1],
      p = imp_p[1], W = NA, B = NA,
      n_success = 1
    )
    log_info(sprintf("  [%s] 单集 HR=%.4f (%.4f-%.4f)", trunc_label, imp_hrs[1], imp_ci_l[1], imp_ci_u[1]))
  } else {
    all_results[[trunc_label]] <- list(
      hr = NA, ci_lower = NA, ci_upper = NA, p = NA, W = NA, B = NA,
      n_success = 0
    )
    log_info(sprintf("  [%s] 全部失败", trunc_label))
  }
}

# ---- 导出汇总 ----
log_section("权重截断敏感性汇总")
summary_rows <- list()
for (trunc in TRUNC_LEVELS) {
  label <- sprintf("%.0f%%", trunc * 100)
  r <- all_results[[label]]
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    truncation = label,
    hr = r$hr,
    ci_lower = r$ci_lower,
    ci_upper = r$ci_upper,
    p = r$p,
    n_success = r$n_success
  )
  log_info(sprintf("截断 %s: HR=%.4f (%.4f-%.4f), p=%.4f, success=%d/%d",
                   label, r$hr %||% NA, r$ci_lower %||% NA, r$ci_upper %||% NA,
                   r$p %||% NA, r$n_success, N_IMPUTE))
}
dt_summary <- rbindlist(summary_rows)
fwrite(dt_summary, file.path(TBL_DIR, "v26_9_sensitivity_trunc_summary.csv"))

log_info("权重截断敏感性分析完成")
close(log_con)