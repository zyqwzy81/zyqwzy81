# ============================================================================
# V26.9 敏感性分析 — DAG13（排除 sepsis3）
# 审稿 M-5 响应: 评估 sepsis3 缺失对结果的影响
# 在 MIMIC 全数据上使用 DAG13（DAG14 - sepsis3）运行 SEQTaRget ITT
# 对比 DAG14（完整） vs DAG13（排除 sepsis3）的 HR
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(survival)
  library(SEQTaRget)
})

# ---- 路径 ----
DATA_FILE       <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_imputed_v26.parquet"
MICE_IMP01      <- "/data/小米云盘/文档/论文/血磷/sql/result/mice_imputed/v26_mice_imputed_imp01.parquet"
OUT_DIR         <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9"
TBL_DIR         <- file.path(OUT_DIR, "tables")
FIG_DIR         <- file.path(OUT_DIR, "figures")
LOG_DIR         <- file.path(OUT_DIR, "logs")
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

LOGFILE <- file.path(LOG_DIR, sprintf("v26_9_sensitivity_dag13_%s.log",
                                       format(Sys.time(), "%Y%m%d_%H%M%S")))

log_con <- file(LOGFILE, open = "wt")
log_info <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg), file = log_con, append = TRUE)
  flush(stdout())
}
log_section <- function(title) {
  cat(sprintf("\n========== %s ==========\n", title))
  cat(sprintf("\n========== %s ==========\n", title), file = log_con, append = TRUE)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
BOOT <- 200L
SEED <- 20260629L
set.seed(SEED)
options(mc.cores = 1)

# ---- DAG14 vs DAG13 ----
DAG14 <- c("age", "gender", "baseline_sofa", "sepsis3", "charlson",
           "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
           "admission_type", "baseline_phosphate", "lactate_t0", "refeeding_risk",
           "corrected_calcium_mmol_min", "nutrition_support_type")
DAG13 <- setdiff(DAG14, "sepsis3")

log_info(sprintf("DAG14: %d 变量", length(DAG14)))
log_info(sprintf("DAG13: %d 变量 (排除 sepsis3)", length(DAG13)))

# ============================================================================
# 1. 数据加载
# ============================================================================
log_section("数据加载")

dt_full <- as.data.table(read_parquet(DATA_FILE))
dt_full[, `:=`(
  tx_int = as.integer(itt_group == "treatment"),
  event  = death_28day_num,
  id_int = as.integer(factor(stay_id))
)]
setorder(dt_full, id_int, trial_num)
log_info(sprintf("主 Parquet: %d 行, %d 事件", nrow(dt_full), sum(dt_full$event)))

key_cols <- c("stay_id", "subject_id", "trial_num", "time_zero",
              "itt_group", "gender", "admission_type",
              "surv_time", "death_time_from_t0", "dischtime_from_t0",
              "race_group", "icu_type",
              "death_28day_num", "nutrition_support_type")
key_cols <- intersect(key_cols, names(dt_full))
dt_key <- dt_full[, ..key_cols]

# MICE Parquet imp01
dt_imp <- as.data.table(read_parquet(MICE_IMP01))
keep <- intersect(c("stay_id", "trial_num", DAG14), names(dt_imp))
dt_imp <- dt_imp[, ..keep]

# merge 对齐
dag14_in_imp <- intersect(DAG14, names(dt_imp))
remove_from_key <- intersect(dag14_in_imp, names(dt_key))
dt_comb <- merge(
  dt_key[, setdiff(names(dt_key), remove_from_key), with = FALSE],
  dt_imp, by = c("stay_id", "trial_num"), all.x = TRUE, sort = FALSE
)
log_info(sprintf("merge 后: %d 行, DAG14 可用=%d/14",
                 nrow(dt_comb), length(intersect(DAG14, names(dt_comb)))))

# ============================================================================
# 2. SEQTaRget 运行函数
# ============================================================================
run_seqtarget_sens <- function(dt_comb, cov_set, label, boot = BOOT) {
  dt_sub <- copy(dt_comb)
  dt_sub[, `:=`(
    id_int       = as.integer(factor(stay_id)),
    period_int   = as.integer(trial_num) - 1L,
    a_int        = as.integer(itt_group == "treatment"),
    y_int        = as.integer(death_28day_num),
    event        = as.integer(death_28day_num),
    eligible_int = 1L
  )]
  setorder(dt_sub, id_int, period_int)
  dt_sub[, `:=`(n_a = sum(a_int), n_y = sum(y_int)), by = period_int]
  dt_sub <- dt_sub[n_a > 0 & n_a < .N & n_y >= 5]

  covs <- intersect(cov_set, names(dt_sub))
  covs <- covs[sapply(covs, function(v) {
    x <- dt_sub[[v]]
    if (is.factor(x))   return(nlevels(x) > 1)
    if (is.character(x)) return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]
  for (v in covs) {
    if (is.character(dt_sub[[v]])) dt_sub[, (v) := as.factor(get(v))]
  }

  if (nrow(dt_sub) < 100 || length(covs) < 2) return(NULL)

  t0 <- Sys.time()
  result <- tryCatch({
    SEQuential(
      data = as.data.frame(dt_sub),
      id.col = "id_int", time.col = "period_int",
      eligible.col = "eligible_int",
      treatment.col = "a_int", outcome.col = "y_int",
      fixed.cols = covs, time_varying.cols = character(0),
      method = "ITT",
      options = SEQopts(
        bootstrap = TRUE, bootstrap.nboot = boot, bootstrap.CI = 0.95,
        ncores = 1L, fastglm.method = 3L, hazard = TRUE,
        trial.include = TRUE, followup.include = TRUE, data.return = FALSE
      ),
      verbose = FALSE
    )
  }, error = function(e) {
    log_info(sprintf("  [%s] 错误: %s", label, e$message))
    NULL
  })
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (is.null(result)) return(NULL)
  hr_vec <- hazard_ratio(result)
  list(
    hr = hr_vec[1], ci_lower = hr_vec[2], ci_upper = hr_vec[3],
    runtime = runtime, n = nrow(dt_sub),
    n_events = sum(dt_sub$event),
    label = label, covs_used = length(covs)
  )
}

# ============================================================================
# 3. 运行 DAG14（参考）和 DAG13（敏感性）
# ============================================================================
log_section("DAG14 参考分析 (boot=200)")
r14 <- run_seqtarget_sens(dt_comb, DAG14, "DAG14_full", boot = BOOT)
if (!is.null(r14)) {
  log_info(sprintf("  DAG14: HR=%.4f (%.4f–%.4f), n=%d, events=%d, %.1f min",
                   r14$hr, r14$ci_lower, r14$ci_upper,
                   r14$n, r14$n_events, r14$runtime / 60))
} else {
  log_info("  DAG14: 返回 NULL")
}

log_section("DAG13 敏感性分析 (排除 sepsis3, boot=200)")
r13 <- run_seqtarget_sens(dt_comb, DAG13, "DAG13_no_sepsis3", boot = BOOT)
if (!is.null(r13)) {
  log_info(sprintf("  DAG13: HR=%.4f (%.4f–%.4f), n=%d, events=%d, %.1f min",
                   r13$hr, r13$ci_lower, r13$ci_upper,
                   r13$n, r13$n_events, r13$runtime / 60))
} else {
  log_info("  DAG13: 返回 NULL")
}

# ============================================================================
# 4. 对比与输出
# ============================================================================
log_section("结果对比")

if (!is.null(r14) && !is.null(r13)) {
  # 差异计算
  delta_hr <- abs(r13$hr - r14$hr) / r14$hr * 100
  ci_overlap <- (r13$ci_lower <= r14$ci_upper) & (r13$ci_upper >= r14$ci_lower)
  same_direction <- (r14$hr - 1) * (r13$hr - 1) > 0

  log_info("===== DAG14 vs DAG13 对比 =====")
  log_info(sprintf("  DAG14  HR=%.4f (%.4f–%.4f)", r14$hr, r14$ci_lower, r14$ci_upper))
  log_info(sprintf("  DAG13  HR=%.4f (%.4f–%.4f) [排除 sepsis3]", r13$hr, r13$ci_lower, r13$ci_upper))
  log_info(sprintf("  HR 相对差异: %.1f%%", delta_hr))
  log_info(sprintf("  方向一致性: %s", if (same_direction) "一致" else "不一致"))
  ci_overlap_label <- if (ci_overlap) "overlap" else "no overlap"
    log_info(sprintf("  95%% CI overlap: %s", ci_overlap_label))

  if (delta_hr < 10 && same_direction && ci_overlap) {
    log_info("  结论: DAG13 结果与 DAG14 一致，sepsis3 缺失影响有限")
  } else {
    log_info("  结论: sepsis3 缺失对结果有较大影响，需在论文中报告")
  }

  # 保存结果
  res <- data.table(
    analysis = c("DAG14_full", "DAG13_no_sepsis3"),
    covariates = c(14L, 13L),
    hr = c(r14$hr, r13$hr),
    ci_lower = c(r14$ci_lower, r13$ci_lower),
    ci_upper = c(r14$ci_upper, r13$ci_upper),
    n = c(r14$n, r13$n),
    n_events = c(r14$n_events, r13$n_events),
    runtime_min = round(c(r14$runtime, r13$runtime) / 60, 1)
  )

  # 差异行
  res <- rbind(res, data.table(
    analysis = "差异", covariates = NA_integer_,
    hr = delta_hr, ci_lower = NA_real_, ci_upper = NA_real_,
    n = NA_integer_, n_events = NA_integer_, runtime_min = NA_real_
  ), fill = TRUE)

  out_file <- file.path(TBL_DIR, "v26_9_sensitivity_dag13_comparison.csv")
  fwrite(res, out_file)
  log_info(sprintf("对比结果已保存: %s", out_file))
} else {
  log_info("对比失败: DAG14 或 DAG13 分析结果为 NULL")
}

log_section("敏感性分析完成")
close(log_con)
cat(sprintf("\n日志: %s\n", LOGFILE))
