# ============================================================================
# V26.9 MITS-LASSO vs DAG14 协变量集一致性对比（28天死亡率）
#
# 目的：比较两种协变量选择策略下 SEQTaRget ITT 分析结果的一致性
#   - Model A: DAG14（理论驱动，14变量）— 当前主分析
#   - Model B: LASSO ≥50% 频率变量（数据驱动，44变量）— 敏感性分析
#
# 方法：SEQTaRget CCW-ITT ×20 MICE + Rubin 合并
#   - 28天死亡率为主要结局
#   - 使用已有 MICE 插补数据
#   - nboot=100（敏感性分析，略低于主分析的 200）
#
# 输出：
#   - 两种模型的 HR (95% CI) 及 p 值对比表
#   - 包含 DAG14 变量在 LASSO 中的选择频率交叉表
#   - 纳入论文补充材料 "Sensitivity Analysis: covariate selection strategy"
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

OUT_DIR  <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9"
RDS_DIR  <- file.path(OUT_DIR, "rds")
TBL_DIR  <- file.path(OUT_DIR, "tables")
LOG_DIR  <- file.path(OUT_DIR, "logs")
for (d in c(OUT_DIR, RDS_DIR, TBL_DIR, LOG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

DATA_FILE           <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_imputed_v26.parquet"
MICE_PARQUET_TEMPL  <- "/data/小米云盘/文档/论文/血磷/sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet"
LASSO_FREQ_FILE     <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/lasso/v26_lasso_selection_frequency.csv"

N_IMPUTE  <- 20L       # 全部 20 个插补集
BOOT      <- 100L      # 敏感性分析，nboot=100
SEED      <- 20260629L
set.seed(SEED)
options(mc.cores = 1)

# ---- DAG14 协变量集 ----
DAG14 <- c(
  "age", "gender", "baseline_sofa", "sepsis3", "charlson",
  "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
  "admission_type", "baseline_phosphate", "lactate_t0", "refeeding_risk",
  "corrected_calcium_mmol_min", "nutrition_support_type"
)

# ---- 读取 LASSO 选择频率，筛选 ≥50% 变量 ----
dt_lasso <- fread(LASSO_FREQ_FILE)
LASSO_50 <- dt_lasso[Selection_Pct >= 50, sort(Variable)]
cat(sprintf("LASSO ≥50%% 变量数: %d\n", length(LASSO_50)))
cat(sprintf("DAG14 变量数: %d\n", length(DAG14)))

# 交叉表：DAG14 变量在 LASSO 中的选择频率
dag_in_lasso <- intersect(DAG14, dt_lasso$Variable)
cat("\nDAG14 变量的 LASSO 选择频率:\n")
for (v in DAG14) {
  freq <- if (v %in% dag_in_lasso) {
    dt_lasso[Variable == v, Selection_Pct]
  } else {
    0L
  }
  cat(sprintf("  %-35s %3d%% %s\n", v, freq,
              if (freq >= 50) "[入选]" else "[未入选]"))
}

# ============================================================================
# 2. 日志系统
# ============================================================================
log_file <- file.path(LOG_DIR, sprintf("v26_9_lasso_vs_dag14_%s.log",
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

log_info("V26.9 LASSO vs DAG14 一致性对比分析启动")
log_info(sprintf("种子=%d, nboot=%d, ncores=1, 插补集: 1-%d", SEED, BOOT, N_IMPUTE))
log_info(sprintf("模型A (DAG14): %d 变量", length(DAG14)))
log_info(sprintf("模型B (LASSO≥50%%): %d 变量", length(LASSO_50)))

# ============================================================================
# 3. 数据加载
# ============================================================================
log_section("数据加载")

# 加载主数据（仅 key 列 + 28天死亡信息）
dt_full <- as.data.table(read_parquet(DATA_FILE))
dt_full[, `:=`(
  tx_int = as.integer(itt_group == "treatment"),
  id_int = as.integer(factor(stay_id))
)]
setorder(dt_full, id_int, trial_num)

key_cols <- c("stay_id", "subject_id", "trial_num", "time_zero",
              "itt_group",
              "surv_time", "death_28day_num")
key_cols <- intersect(key_cols, names(dt_full))
dt_key   <- dt_full[, ..key_cols]
log_info(sprintf("主数据 key 表: %d 行, %d 列", nrow(dt_key), ncol(dt_key)))
rm(dt_full); gc(TRUE)

# ============================================================================
# 4. 辅助函数
# ============================================================================
log_section("辅助函数定义")

# 通用 SEQTaRget ITT 单次运行
run_seqtarget_covset <- function(dt_key, cov_set, label = "") {
  dt_sub <- copy(dt_key)

  dt_sub[, `:=`(
    id_int       = as.integer(factor(stay_id)),
    period_int   = as.integer(trial_num) - 1L,
    a_int        = as.integer(itt_group == "treatment"),
    y_int        = as.integer(death_28day_num),
    eligible_int = 1L
  )]
  setorder(dt_sub, id_int, period_int)

  # 按 period 过滤
  dt_sub[, `:=`(n_a = sum(a_int), n_y = sum(y_int)), by = period_int]
  dt_sub <- dt_sub[n_a > 0 & n_a < .N & n_y >= 5]
  if (nrow(dt_sub) < 100) {
    log_info(sprintf("  [%s] 过滤后样本不足: %d", label, nrow(dt_sub)))
    return(NULL)
  }

  # 协变量准备：可用交集 + 排除零方差
  covs <- intersect(cov_set, names(dt_sub))
  covs <- covs[sapply(covs, function(v) {
    x <- dt_sub[[v]]
    if (is.factor(x))    return(nlevels(x) > 1)
    if (is.character(x)) return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]
  if (length(covs) < 2) {
    log_info(sprintf("  [%s] 可用协变量不足: %d", label, length(covs)))
    return(NULL)
  }
  for (v in covs) {
    if (is.character(dt_sub[[v]])) dt_sub[[v]] <- as.factor(dt_sub[[v]])
  }

  t0 <- Sys.time()
  result <- tryCatch({
    SEQuential(
      data              = as.data.frame(dt_sub),
      id.col            = "id_int",
      time.col          = "period_int",
      eligible.col      = "eligible_int",
      treatment.col     = "a_int",
      outcome.col       = "y_int",
      fixed.cols        = covs,
      time_varying.cols = character(0),
      method            = "ITT",
      options = SEQopts(
        bootstrap       = TRUE,
        bootstrap.nboot = BOOT,
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
  }, error = function(e) NULL)
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (is.null(result)) {
    log_info(sprintf("  [%s] SEQTaRget 返回 NULL", label))
    return(NULL)
  }

  hr_vec <- tryCatch(hazard_ratio(result), error = function(e) c(NA, NA, NA))
  rm(result); gc(TRUE)

  return(list(
    hr       = hr_vec[1],
    ci_lower = hr_vec[2],
    ci_upper = hr_vec[3],
    runtime  = runtime,
    n        = nrow(dt_sub),
    n_events = sum(dt_sub$y_int == 1, na.rm = TRUE)
  ))
}

# Rubin 合并
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
# 5. 核心分析：遍历两种协变量集
# ============================================================================
log_section("核心分析：两种协变量集对比")

cov_sets <- list(
  DAG14       = DAG14,
  LASSO_50pct = LASSO_50
)

all_models <- list()

for (model_name in names(cov_sets)) {
  cov_set <- cov_sets[[model_name]]
  log_info(sprintf("\n======== 模型: %s (%d 变量) ========", model_name, length(cov_set)))

  imp_betas   <- numeric()
  imp_ses     <- numeric()
  imp_hrs     <- numeric()
  imp_ci_l    <- numeric()
  imp_ci_u    <- numeric()
  imp_run     <- numeric()
  imp_nevt    <- numeric()
  imp_success <- integer()
  imp_fail    <- integer()

  t_model_start <- Sys.time()

  for (i_imp in seq_len(N_IMPUTE)) {
    t_imp_start <- Sys.time()

    # ---- 5a. 加载 MICE Parquet ----
    parquet_file <- sprintf(MICE_PARQUET_TEMPL, i_imp)
    if (!file.exists(parquet_file)) {
      log_info(sprintf("  [imp%02d/%s] 文件不存在", i_imp, model_name))
      imp_fail <- c(imp_fail, i_imp)
      next
    }
    dt_imp <- as.data.table(read_parquet(parquet_file))

    # 保留该模型需要的协变量
    keep_cols <- intersect(c("stay_id", "trial_num", cov_set), names(dt_imp))
    dt_imp <- dt_imp[, ..keep_cols]
    gc(TRUE)

    # ---- 5b. Merge 对齐 ----
    covs_in_imp <- intersect(cov_set, names(dt_imp))
    covs_remove <- intersect(covs_in_imp, names(dt_key))
    by_cols <- c("stay_id", "trial_num")

    dt_comb <- merge(
      dt_key[, setdiff(names(dt_key), covs_remove), with = FALSE],
      dt_imp,
      by = by_cols, all.x = TRUE, sort = FALSE
    )

    # 字符型→factor
    for (v in intersect(names(dt_comb), names(dt_imp))) {
      if (is.character(dt_comb[[v]])) dt_comb[[v]] <- as.factor(dt_comb[[v]])
    }

    # ---- 5c. 运行 SEQTaRget ----
    r <- tryCatch(
      run_seqtarget_covset(dt_comb, cov_set, label = sprintf("%s_imp%02d", model_name, i_imp)),
      error = function(e) NULL
    )

    t_imp_end <- Sys.time()
    imp_elapsed <- as.numeric(difftime(t_imp_end, t_imp_start, units = "secs"))

    if (is.null(r)) {
      log_info(sprintf("  [imp%02d/%s] 失败, %.1f min", i_imp, model_name, imp_elapsed / 60))
      imp_fail <- c(imp_fail, i_imp)
      rm(dt_comb, dt_imp); gc(TRUE)
      next
    }

    log_info(sprintf("  [imp%02d/%s] HR=%.4f (%.4f–%.4f), events=%d, %.1f min",
                     i_imp, model_name, r$hr, r$ci_lower, r$ci_upper,
                     r$n_events, r$runtime / 60))

    imp_betas   <- c(imp_betas,   log(r$hr))
    imp_ses     <- c(imp_ses,     (log(r$ci_upper) - log(r$ci_lower)) / (2 * 1.96))
    imp_hrs     <- c(imp_hrs,     r$hr)
    imp_ci_l    <- c(imp_ci_l,    r$ci_lower)
    imp_ci_u    <- c(imp_ci_u,    r$ci_upper)
    imp_run     <- c(imp_run,     r$runtime)
    imp_nevt    <- c(imp_nevt,    r$n_events)
    imp_success <- c(imp_success, i_imp)

    rm(dt_comb, dt_imp, r); gc(TRUE)
  }

  t_model_end <- Sys.time()
  model_elapsed <- as.numeric(difftime(t_model_end, t_model_start, units = "secs"))

  m_actual <- length(imp_success)

  # ---- 5d. Rubin 合并 ----
  if (m_actual >= 2) {
    rubin <- rubin_combine(imp_betas, imp_ses)

    log_info(sprintf("\n  [合并/%s] HR=%.4f (%.4f–%.4f), p=%.6f, m=%d",
                     model_name, rubin$hr, rubin$ci_lower, rubin$ci_upper,
                     rubin$p, m_actual))

    all_models[[model_name]] <- list(
      cov_set_name = model_name,
      n_covs       = length(cov_set),
      rubin        = rubin,
      m_actual     = m_actual,
      imp_hrs      = imp_hrs,
      imp_betas    = imp_betas,
      imp_ses      = imp_ses,
      imp_nevt     = imp_nevt,
      imp_success  = imp_success,
      imp_fail     = imp_fail,
      runtime_sec  = model_elapsed
    )
  } else {
    log_info(sprintf("  [失败/%s] 不足2个成功插补集", model_name))
  }
}

# ============================================================================
# 6. 结果比较
# ============================================================================
log_section("结果比较：DAG14 vs LASSO")

cat("\n")
cat("============================================================================\n")
cat("  协变量策略一致性对比：28天死亡率\n")
cat("============================================================================\n")

dt_compare_list <- list()

for (model_name in names(all_models)) {
  m <- all_models[[model_name]]
  r <- m$rubin

  # 一致性标记
  hr_diff_text <- "-"

  line <- sprintf("  %-15s %5d %8.4f %8.4f %8.4f %10.6f %4d %8.0f %s",
                  model_name, m$n_covs,
                  r$hr, r$ci_lower, r$ci_upper,
                  r$p, m$m_actual,
                  median(m$imp_nevt, na.rm = TRUE),
                  hr_diff_text)
  cat(line, "\n", file = log_con, append = TRUE)
  cat(line, "\n")

  dt_compare_list[[model_name]] <- data.table(
    model        = model_name,
    n_covariates = m$n_covs,
    hr           = r$hr,
    ci_lower     = r$ci_lower,
    ci_upper     = r$ci_upper,
    p_value      = r$p,
    m_success    = m$m_actual,
    median_events = median(m$imp_nevt, na.rm = TRUE),
    runtime_min  = m$runtime_sec / 60
  )
}

# 若两个模型均成功，计算 HR 相对差异
if (length(all_models) == 2) {
  hr_a <- all_models[[1]]$rubin$hr
  hr_b <- all_models[[2]]$rubin$hr
  hr_ratio <- hr_a / hr_b
  hr_pct_diff <- abs(hr_a - hr_b) / ((hr_a + hr_b) / 2) * 100

  log_info(sprintf("\nHR 相对差异: %.2f%%", hr_pct_diff))
  log_info(sprintf("HR 比值 (DAG14/LASSO): %.4f", hr_ratio))
  log_info(sprintf("一致性判断: %s",
                   if (hr_pct_diff < 5) "一致 (<5%)"
                   else if (hr_pct_diff < 10) "基本一致 (5-10%)"
                   else "不一致 (≥10%)"))

  cat(sprintf("\n  HR 相对差异: %.2f%%\n", hr_pct_diff))
  cat(sprintf("  HR 比值: %.4f\n", hr_ratio))
  cat(sprintf("  一致性判断: %s\n",
              if (hr_pct_diff < 5) "一致 (<5%)"
              else if (hr_pct_diff < 10) "基本一致 (5-10%)"
              else "不一致 (≥10%)"))
}

# 保存对比表
if (length(dt_compare_list) > 0) {
  dt_compare <- rbindlist(dt_compare_list)
  fwrite(dt_compare, file.path(TBL_DIR, "v26_9_lasso_vs_dag14_comparison.csv"))
  log_info(sprintf("对比表已保存: %s",
                   file.path(TBL_DIR, "v26_9_lasso_vs_dag14_comparison.csv")))
}

# ============================================================================
# 7. 生成 DAG14 vs LASSO 变量交叉表（用于论文补充材料）
# ============================================================================
log_section("变量交叉表")

cross_list <- list()

# DAG14 变量
for (v in DAG14) {
  freq <- if (v %in% dt_lasso$Variable) {
    dt_lasso[Variable == v, Selection_Pct]
  } else {
    0L
  }
  in_lasso <- freq >= 50
  cross_list[[v]] <- data.table(
    variable       = v,
    in_dag14       = TRUE,
    lasso_freq_pct = freq,
    in_lasso_50    = in_lasso,
    source         = if (in_lasso) "Both" else "DAG14 only"
  )
}

# LASSO-only 变量（不在 DAG14 中的）
lasso_only <- setdiff(LASSO_50, DAG14)
for (v in lasso_only) {
  freq <- dt_lasso[Variable == v, Selection_Pct]
  cross_list[[v]] <- data.table(
    variable       = v,
    in_dag14       = FALSE,
    lasso_freq_pct = freq,
    in_lasso_50    = TRUE,
    source         = "LASSO only"
  )
}

dt_cross <- rbindlist(cross_list)
setorder(dt_cross, -lasso_freq_pct, variable)
fwrite(dt_cross, file.path(TBL_DIR, "v26_9_lasso_dag14_cross_table.csv"))
log_info(sprintf("变量交叉表已保存: %s",
                 file.path(TBL_DIR, "v26_9_lasso_dag14_cross_table.csv")))

# 汇总统计
n_both <- sum(dt_cross$source == "Both")
n_dag_only <- sum(dt_cross$source == "DAG14 only")
n_lasso_only <- sum(dt_cross$source == "LASSO only")

log_info(sprintf("变量来源汇总:"))
log_info(sprintf("  DAG14 + LASSO 共同: %d", n_both))
log_info(sprintf("  DAG14 独有 (LASSO未选中): %d", n_dag_only))
log_info(sprintf("  LASSO 独有 (DAG14未含): %d", n_lasso_only))

# ============================================================================
# 8. 保存 RDS
# ============================================================================
log_section("保存 RDS")

rds_out <- list(
  version       = "V26.9",
  script        = "v26_9_lasso_vs_dag14_sensitivity.R",
  timestamp     = as.character(Sys.time()),
  seed          = SEED,
  nboot         = BOOT,
  n_impute      = N_IMPUTE,
  dag14         = DAG14,
  lasso_50_vars = LASSO_50,
  comparison    = dt_compare,
  cross_table   = dt_cross,
  per_model     = all_models
)
saveRDS(rds_out, file.path(RDS_DIR, "v26_9_lasso_vs_dag14.rds"))
log_info(sprintf("RDS 已保存: %s",
                 file.path(RDS_DIR, "v26_9_lasso_vs_dag14.rds")))

# ============================================================================
# 9. 流程自审核
# ============================================================================
log_section("流程自审核")

audit_checks <- list(
  c("数据加载", nrow(dt_key) > 10000, sprintf("n=%d", nrow(dt_key))),
  c("DAG14模型成功", "DAG14" %in% names(all_models),
    sprintf("DAG14存在=%s", "DAG14" %in% names(all_models))),
  c("LASSO模型成功", "LASSO_50pct" %in% names(all_models),
    sprintf("LASSO_50pct存在=%s", "LASSO_50pct" %in% names(all_models))),
  c("对比表已生成", exists("dt_compare") && nrow(dt_compare) > 0,
    sprintf("n=%d", if (exists("dt_compare")) nrow(dt_compare) else 0)),
  c("变量交叉表", nrow(dt_cross) > 0,
    sprintf("n=%d", nrow(dt_cross)))
)

for (ac in audit_checks) {
  pass <- ac[2] == TRUE || ac[2] == "TRUE"
  log_info(sprintf("  [%s] %s: %s",
                   if (pass) "PASS" else "FAIL", ac[1], ac[3]))
}

log_info(sprintf("\n日志文件: %s", log_file))
log_info(sprintf("对比表: %s",
                 file.path(TBL_DIR, "v26_9_lasso_vs_dag14_comparison.csv")))
log_info(sprintf("变量交叉表: %s",
                 file.path(TBL_DIR, "v26_9_lasso_dag14_cross_table.csv")))

close(log_con)

cat("\n")
cat("============================================================================\n")
cat("  V26.9 LASSO vs DAG14 一致性对比完成\n")
cat("============================================================================\n")
cat(sprintf("  模型A (DAG14):       %d 变量\n", length(DAG14)))
cat(sprintf("  模型B (LASSO ≥50%%): %d 变量\n", length(LASSO_50)))
cat(sprintf("  共同变量:           %d\n", n_both))
cat(sprintf("  DAG14独有:          %d (gender, admission_type, nutrition_support_type)\n", n_dag_only))
cat(sprintf("  LASSO独有:          %d\n", n_lasso_only))
cat(sprintf("  日志: %s\n", log_file))
cat(sprintf("  对比表: %s\n",
            file.path(TBL_DIR, "v26_9_lasso_vs_dag14_comparison.csv")))
cat("============================================================================\n")
