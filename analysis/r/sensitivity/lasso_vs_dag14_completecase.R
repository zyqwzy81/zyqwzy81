# ============================================================================
# V26.9 MITS-LASSO vs DAG14 一致性对比（完整数据集，快速版）
#
# 目的：使用完整病例数据快速比较两种协变量集下 28 天死亡率的 HR 一致性
#   - Model A: DAG14（14 个理论驱动变量）
#   - Model B: LASSO ≥50%（44 个数据驱动变量）
#
# 方法：两阶段 IPTW + 加权 Cox（避免 SEQTaRget 依赖问题）
#   阶段1: Logistic PS 模型 → 稳定化 IPTW 权重（双侧截断）
#   阶段2: 加权 Cox 回归 + cluster-robust SE
#   阶段3: 对比 HR (95% CI)
#
# 注意：完整个案 | 快速 | 无 bootstrap | 敏感性分析用途
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(survival)
  library(sandwich)
  library(lmtest)
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

AUG_FILE <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_augmented.parquet"
LASSO_FREQ_FILE <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/lasso/v26_lasso_selection_frequency.csv"

DAG14 <- c(
  "age", "gender", "baseline_sofa", "sepsis3", "charlson",
  "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
  "admission_type", "baseline_phosphate", "lactate_t0", "refeeding_risk",
  "corrected_calcium_mmol_min", "nutrition_support_type"
)

IPTW_TRUNC_UPPER <- 0.99
IPTW_TRUNC_LOWER <- 0.01
SEED <- 20260629L
set.seed(SEED)

# ============================================================================
# 2. 日志系统
# ============================================================================
log_file <- file.path(LOG_DIR, sprintf("v26_9_lasso_vs_dag14_cc_%s.log",
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

# ============================================================================
# 3. 数据准备
# ============================================================================
log_section("数据准备")

# ---- 3a. LASSO 变量集 ----
dt_lasso <- fread(LASSO_FREQ_FILE)
LASSO_50 <- dt_lasso[Selection_Pct >= 50, sort(Variable)]
log_info(sprintf("LASSO ≥50%% 变量: %d 个", length(LASSO_50)))
log_info(sprintf("DAG14 变量: %d 个", length(DAG14)))

cat("\nDAG14 变量的 LASSO 选择频率:\n")
for (v in DAG14) {
  freq <- if (v %in% dt_lasso$Variable) dt_lasso[Variable == v, Selection_Pct] else 0L
  cat(sprintf("  %-35s %3d%% %s\n", v, freq, if (freq >= 50) "[入选]" else "[未入选]"))
}

# ---- 3b. 加载增强数据集 ----
log_info("加载增强数据集...")
dt_raw <- as.data.table(read_parquet(AUG_FILE))
log_info(sprintf("增强数据集: %d 行, %d 列", nrow(dt_raw), ncol(dt_raw)))

all_covs <- unique(c(DAG14, LASSO_50))
all_covs <- intersect(all_covs, names(dt_raw))
log_info(sprintf("可用候选变量: %d / %d", length(all_covs), length(unique(c(DAG14, LASSO_50)))))

# ---- 3c. 筛选完整个案 ----
required_cols <- c("stay_id", "trial_num", "itt_group",
                   "death_28day_num", "surv_time", all_covs)
required_cols <- intersect(required_cols, names(dt_raw))

dt_cc <- dt_raw[, ..required_cols]
dt_cc <- dt_cc[complete.cases(dt_cc)]
log_info(sprintf("完整个案: %d 行 (原始 %d 行, 保留 %.1f%%)",
                 nrow(dt_cc), nrow(dt_raw), 100 * nrow(dt_cc) / nrow(dt_raw)))

if (nrow(dt_cc) < 1000) {
  close(log_con)
  stop("完整个案不足 1000 行")
}

n_events_cc <- sum(dt_cc$death_28day_num == 1)
log_info(sprintf("28 天死亡事件: %d (%.1f%%)",
                 n_events_cc, 100 * n_events_cc / nrow(dt_cc)))

# 创建治疗指示变量
dt_cc[, tx := as.integer(itt_group == "treatment")]
log_info(sprintf("治疗组: %d (%.1f%%), 对照组: %d (%.1f%%)",
                 sum(dt_cc$tx == 1), 100 * mean(dt_cc$tx == 1),
                 sum(dt_cc$tx == 0), 100 * mean(dt_cc$tx == 0)))

# ============================================================================
# 4. 辅助函数
# ============================================================================
log_section("辅助函数")

run_iptw_cox <- function(dt_data, cov_set, label = "") {
  # ---- 4a. 协变量准备 ----
  covs <- intersect(cov_set, names(dt_data))
  covs <- covs[sapply(covs, function(v) {
    x <- dt_data[[v]]
    if (is.factor(x))    return(nlevels(x) > 1)
    if (is.character(x)) return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]
  if (length(covs) < 2) {
    log_info(sprintf("  [%s] 可用协变量不足: %d", label, length(covs)))
    return(NULL)
  }

  # 字符型转 factor
  dt_sub <- copy(dt_data)
  dt_sub[, `:=`(trunc_wt = 1.0)]  # 默认权重
  for (v in covs) {
    if (is.character(dt_sub[[v]])) dt_sub[[v]] <- as.factor(dt_sub[[v]])
  }

  # ---- 4b. Logistic PS 模型 ----
  ps_formula <- as.formula(paste("tx ~", paste(covs, collapse = " + ")))
  t0 <- Sys.time()

  ps_fit <- tryCatch(
    glm(ps_formula, data = dt_sub, family = binomial(link = "logit"),
        control = list(maxit = 100)),
    error = function(e) NULL
  )

  if (is.null(ps_fit)) {
    log_info(sprintf("  [%s] PS 模型失败", label))
    return(NULL)
  }

  # ---- 4c. 稳定化 IPTW ----
  dt_sub[, ps := predict(ps_fit, type = "response")]
  dt_sub[, iptw := fifelse(tx == 1, 1 / ps, 1 / (1 - ps))]

  # 稳定化
  p_tx <- mean(dt_sub$tx)
  dt_sub[, iptw_stable := fifelse(tx == 1, p_tx / ps, (1 - p_tx) / (1 - ps))]

  # 双侧截断
  trunc_u <- quantile(dt_sub$iptw_stable, IPTW_TRUNC_UPPER, na.rm = TRUE)
  trunc_l <- quantile(dt_sub$iptw_stable, IPTW_TRUNC_LOWER, na.rm = TRUE)
  dt_sub[, iptw_trunc := pmin(pmax(iptw_stable, trunc_l), trunc_u)]

  ps_runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  log_info(sprintf("  [%s] PS 模型: %.1f sec, IPTW range: %.2f–%.2f (trunc: %.2f–%.2f)",
                   label, ps_runtime,
                   min(dt_sub$iptw_stable), max(dt_sub$iptw_stable),
                   trunc_l, trunc_u))

  # ---- 4d. 加权 Cox 回归 ----
  t0 <- Sys.time()

  cox_fit <- tryCatch(
    coxph(Surv(surv_time, death_28day_num) ~ tx, data = dt_sub,
          weights = iptw_trunc),
    error = function(e) NULL
  )

  if (is.null(cox_fit)) {
    log_info(sprintf("  [%s] Cox 模型失败", label))
    return(NULL)
  }

  cox_runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # 优先使用 cluster-robust SE
  cox_se <- tryCatch(sqrt(diag(vcovCL(cox_fit, cluster = dt_sub$stay_id)))[1],
                     error = function(e) NULL)

  if (is.null(cox_se)) {
    # 回退到 model-based SE（唯一系数 "tx"）
    cox_se <- sqrt(diag(vcov(cox_fit)))[1]
    se_type <- "model-based"
  } else {
    se_type <- "cluster-robust"
  }

  beta_tx <- coef(cox_fit)[1]
  hr <- exp(beta_tx)
  ci_l <- exp(beta_tx - 1.96 * cox_se)
  ci_u <- exp(beta_tx + 1.96 * cox_se)
  p_val <- 2 * pnorm(abs(beta_tx / cox_se), lower.tail = FALSE)

  log_info(sprintf("  [%s] HR=%.4f (%.4f–%.4f), p=%.6f, se=%s, events=%d, %.1f sec",
                   label, hr, ci_l, ci_u, p_val, se_type,
                   sum(dt_sub$death_28day_num == 1), cox_runtime))

  rm(ps_fit, dt_sub); gc(TRUE)

  return(list(
    hr = hr, ci_lower = ci_l, ci_upper = ci_u,
    p_value = p_val, se_type = se_type,
    n = nrow(dt_data), n_events = sum(dt_data$death_28day_num == 1),
    n_covs = length(covs),
    runtime = ps_runtime + cox_runtime
  ))
}

log_info("辅助函数定义完成")

# ============================================================================
# 5. 核心分析
# ============================================================================
log_section("核心分析：DAG14 vs LASSO (完整个案)")

cov_sets <- list(
  DAG14       = DAG14,
  LASSO_50pct = LASSO_50
)

all_results <- list()

for (model_name in names(cov_sets)) {
  cov_set <- cov_sets[[model_name]]
  log_info(sprintf("\n--- 模型: %s (%d 变量) ---", model_name, length(cov_set)))

  r <- run_iptw_cox(dt_cc, cov_set, label = model_name)

  if (!is.null(r)) {
    all_results[[model_name]] <- r
  }
}

# ============================================================================
# 6. 结果对比
# ============================================================================
log_section("结果对比")

cat("\n")
cat("============================================================================\n")
cat("  协变量策略一致性对比：28天死亡率（完整个案，IPTW + Cox）\n")
cat("============================================================================\n")
cat(sprintf("  样本量: %d 行, %d 事件\n", nrow(dt_cc), n_events_cc))
cat("----------------------------------------------------------------------\n")
cat(sprintf("%-15s %5s %8s %8s %8s %10s %5s\n",
            "Model", "Covs", "HR", "CI_low", "CI_up", "p_value", "Events"))
cat("----------------------------------------------------------------------\n")

dt_compare_list <- list()

for (model_name in names(all_results)) {
  r <- all_results[[model_name]]
  line <- sprintf("  %-15s %5d %8.4f %8.4f %8.4f %10.6f %5d",
                  model_name, r$n_covs, r$hr, r$ci_lower, r$ci_upper,
                  r$p_value, r$n_events)
  cat(line, "\n")
  cat(line, "\n", file = log_con, append = TRUE)

  dt_compare_list[[model_name]] <- data.table(
    model        = model_name,
    n_covariates = r$n_covs,
    hr           = r$hr,
    ci_lower     = r$ci_lower,
    ci_upper     = r$ci_upper,
    p_value      = r$p_value,
    se_type      = r$se_type,
    n_events     = r$n_events,
    n_complete   = r$n,
    runtime_sec  = r$runtime
  )
}

# 一致性评估
if (length(all_results) == 2) {
  hr_a <- all_results[["DAG14"]]$hr
  hr_b <- all_results[["LASSO_50pct"]]$hr
  hr_ratio <- hr_a / hr_b
  hr_pct_diff <- abs(hr_a - hr_b) / ((hr_a + hr_b) / 2) * 100

  ci_a_l <- all_results[["DAG14"]]$ci_lower
  ci_a_u <- all_results[["DAG14"]]$ci_upper
  ci_b_l <- all_results[["LASSO_50pct"]]$ci_lower
  ci_b_u <- all_results[["LASSO_50pct"]]$ci_upper
  overlap <- min(ci_a_u, ci_b_u) - max(ci_a_l, ci_b_l)
  ci_union_len <- (ci_a_u - ci_a_l) + (ci_b_u - ci_b_l)
  overlap_ratio <- if (ci_union_len > 0) overlap / (ci_union_len - overlap) else 0

  same_direction <- (hr_a > 1 && hr_b > 1) || (hr_a < 1 && hr_b < 1)

  cat("----------------------------------------------------------------------\n")
  cat(sprintf("  HR 比值 (DAG14/LASSO): %.4f\n", hr_ratio))
  cat(sprintf("  HR 相对差异: %.2f%%\n", hr_pct_diff))
  cat(sprintf("  CI 重叠率: %.1f%%\n", overlap_ratio * 100))
  cat(sprintf("  方向一致性: %s\n", if (same_direction) "一致" else "不一致"))
  cat(sprintf("  总体判断: %s\n",
              if (hr_pct_diff < 10 && same_direction)
                "协变量选择策略支持结果稳健性"
              else "协变量选择策略需关注"))

  log_info(sprintf("HR 比值: %.4f, 差异: %.2f%%, CI重叠: %.1f%%, 方向一致: %s",
                   hr_ratio, hr_pct_diff, overlap_ratio * 100,
                   if (same_direction) "Yes" else "No"))
}

# 保存对比表
if (length(dt_compare_list) > 0) {
  dt_compare <- rbindlist(dt_compare_list)
  fwrite(dt_compare, file.path(TBL_DIR, "v26_9_lasso_vs_dag14_cc_comparison.csv"))
  log_info(sprintf("对比表: %s",
                   file.path(TBL_DIR, "v26_9_lasso_vs_dag14_cc_comparison.csv")))
}

# ============================================================================
# 7. 变量交叉表
# ============================================================================
log_section("变量交叉表")

cross_list <- list()

for (v in DAG14) {
  freq <- if (v %in% dt_lasso$Variable) dt_lasso[Variable == v, Selection_Pct] else 0L
  cross_list[[v]] <- data.table(
    variable = v, in_dag14 = TRUE,
    lasso_freq_pct = freq,
    in_lasso_50 = freq >= 50,
    source = if (freq >= 50) "Both" else "DAG14 only"
  )
}

lasso_only <- setdiff(LASSO_50, DAG14)
for (v in lasso_only) {
  freq <- dt_lasso[Variable == v, Selection_Pct]
  cross_list[[v]] <- data.table(
    variable = v, in_dag14 = FALSE,
    lasso_freq_pct = freq,
    in_lasso_50 = TRUE,
    source = "LASSO only"
  )
}

dt_cross <- rbindlist(cross_list)
setorder(dt_cross, -lasso_freq_pct, variable)
fwrite(dt_cross, file.path(TBL_DIR, "v26_9_lasso_dag14_cross_table.csv"))

n_both <- sum(dt_cross$source == "Both")
n_dag_only <- sum(dt_cross$source == "DAG14 only")
n_lasso_only <- sum(dt_cross$source == "LASSO only")

# ============================================================================
# 8. 保存 RDS
# ============================================================================
log_section("保存 RDS")

saveRDS(list(
  version    = "V26.9",
  script     = "v26_9_lasso_vs_dag14_completecase.R",
  timestamp  = as.character(Sys.time()),
  method     = "IPTW + Cox (complete case, cluster-robust SE)",
  dag14      = DAG14,
  lasso_50   = LASSO_50,
  n_complete = nrow(dt_cc),
  n_events   = n_events_cc,
  comparison = dt_compare,
  cross_table = dt_cross,
  per_model  = all_results
), file.path(RDS_DIR, "v26_9_lasso_vs_dag14_cc.rds"))

# ============================================================================
# 9. 自审核
# ============================================================================
log_section("流程自审核")

for (ac in list(
  c("数据加载", nrow(dt_cc) > 1000, sprintf("n=%d", nrow(dt_cc))),
  c("DAG14模型", "DAG14" %in% names(all_results), "OK"),
  c("LASSO模型", "LASSO_50pct" %in% names(all_results), "OK"),
  c("对比表", exists("dt_compare") && nrow(dt_compare) > 0,
    sprintf("n=%d", if (exists("dt_compare")) nrow(dt_compare) else 0))
)) {
  pass <- isTRUE(ac[2]) || ac[2] == "TRUE"
  log_info(sprintf("  [%s] %s: %s", if (pass) "PASS" else "FAIL", ac[1], ac[3]))
}

close(log_con)

cat("\n")
cat("============================================================================\n")
cat("  V26.9 LASSO vs DAG14 一致性对比完成\n")
cat("============================================================================\n")
cat(sprintf("  完整个案: %d 行, %d 事件\n", nrow(dt_cc), n_events_cc))
cat(sprintf("  DAG14: %d 变量 | LASSO: %d 变量\n", length(DAG14), length(LASSO_50)))
cat(sprintf("  共同: %d | DAG14独有: %d | LASSO独有: %d\n",
            n_both, n_dag_only, n_lasso_only))
cat(sprintf("  日志: %s\n", log_file))
cat("============================================================================\n")
