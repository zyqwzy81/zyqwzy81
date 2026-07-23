# ============================================================================
# V26.9 DAG11 eICU 外部验证
#
# 使用 eICU 数据运行 SEQTaRget ITT (DAG11协变量集)
# DAG11: 与 eICU 共享的 11 个协变量
#    age, gender, baseline_sofa, charlson, egfr,
#    invasive_vent_at_t0, vasoactive_at_t0, admission_type,
#    baseline_phosphate, lactate_t0, refeeding_risk
#
# 策略:
#   1. eICU 无 MICE 插补 — 使用 CCA (完整案例分析)
#   2. SEQTaRget ITT, boot=100, DAG11
#   3. 与 MIMIC DAG11/DAG14 结果对比
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(survival)
  library(SEQTaRget)
  library(jsonlite)
})

# ============================================================================
# 1. 参数配置
# ============================================================================
OUT_DIR   <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9"
RDS_DIR   <- file.path(OUT_DIR, "rds")
TBL_DIR   <- file.path(OUT_DIR, "tables")
LOG_DIR   <- file.path(OUT_DIR, "logs")
for (d in c(OUT_DIR, RDS_DIR, TBL_DIR, LOG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

EICU_FILE <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/eicu/eicu_sequential_trials_enhanced.parquet"

# ---- DAG11 (eICU 共享协变量集) ----
DAG11 <- c("age", "gender", "baseline_sofa", "charlson",
           "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
           "admission_type", "baseline_phosphate", "lactate_t0",
           "refeeding_risk")

# ---- 分析参数 ----
BOOT      <- 100L
SEED      <- 20260629L
set.seed(SEED)

# ---- 日志系统 ----
log_file <- file.path(LOG_DIR, sprintf("v26_9_eicu_dag11_%s.log",
                                        format(Sys.time(), "%Y%m%d_%H%M%S")))
log_con  <- file(log_file, open = "wt")
log_info <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg),
      file = log_con, append = TRUE)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}
log_section <- function(title) {
  cat(sprintf("\n========== %s ==========\n", title),
      file = log_con, append = TRUE)
  cat(sprintf("\n========== %s ==========\n", title))
  flush(stdout())
}

log_info("V26.9 DAG11 eICU 外部验证启动")
log_info(sprintf("种子: %d, Bootstrap: %d", SEED, BOOT))
log_info(sprintf("DAG11 (%d 变量): %s", length(DAG11), paste(DAG11, collapse = ", ")))
log_info(sprintf("eICU 数据: %s", EICU_FILE))
log_info(sprintf("输出: %s", OUT_DIR))

# ============================================================================
# 2. 数据加载与准备
# ============================================================================
log_section("eICU 数据加载")

dt_raw <- as.data.table(read_parquet(EICU_FILE))
log_info(sprintf("原始 eICU 数据: %d 行, %d 列", nrow(dt_raw), ncol(dt_raw)))
log_info(sprintf("  治疗组: %d (%.1f%%), 对照组: %d (%.1f%%)",
                 dt_raw[itt_group == "treatment", .N],
                 dt_raw[itt_group == "treatment", .N] / nrow(dt_raw) * 100,
                 dt_raw[itt_group != "treatment", .N],
                 dt_raw[itt_group != "treatment", .N] / nrow(dt_raw) * 100))
log_info(sprintf("  28天事件: %d (%.1f%%)",
                 dt_raw$death_28day_num |> sum(na.rm = TRUE),
                 dt_raw$death_28day_num |> sum(na.rm = TRUE) / nrow(dt_raw) * 100))

# ---- 变量准备 ----
dt_raw[, `:=`(
  tx_int = as.integer(itt_group == "treatment"),
  event  = death_28day_num,
  id_int = as.integer(factor(stay_id))
)]
setorder(dt_raw, id_int, trial_num)

# ---- 完整案例筛选 (DAG11) ----
dt_cc <- dt_raw[complete.cases(dt_raw[, ..DAG11])]
log_info(sprintf("完整案例 (DAG11): %d 行, %d 事件 (%.1f%%)",
                 nrow(dt_cc), sum(dt_cc$event),
                 sum(dt_cc$event) / nrow(dt_cc) * 100))
log_info(sprintf("  治疗组: %d, 对照组: %d",
                 dt_cc[tx_int == 1, .N], dt_cc[tx_int == 0, .N]))

# 补充: 排除治疗无变异或事件不足的期
dt_cc[, `:=`(n_a = sum(tx_int), n_y = sum(event)), by = trial_num]
dt_cc <- dt_cc[n_a > 0 & n_a < .N & n_y >= 5]
log_info(sprintf("过滤不可识别期后: %d 行, %d 期",
                 nrow(dt_cc), uniqueN(dt_cc$trial_num)))

if (nrow(dt_cc) < 500) {
  log_info("[ERROR] eICU 完整案例样本量不足，退出")
  close(log_con)
  stop("样本量不足")
}

# ============================================================================
# 3. SEQTaRget ITT 分析 (DAG11)
# ============================================================================
log_section("SEQTaRget ITT × DAG11 (boot=100)")

dt_sub <- copy(dt_cc)
dt_sub[, `:=`(
  period_int   = as.integer(trial_num) - 1L,
  a_int        = tx_int,
  y_int        = as.integer(event),
  eligible_int = 1L
)]
setorder(dt_sub, id_int, period_int)

# 去重 (SEQTaRget 要求 id_int+period_int 唯一)
n_before <- nrow(dt_sub)
dt_sub <- unique(dt_sub, by = c("id_int", "period_int"))
n_dup <- n_before - nrow(dt_sub)
if (n_dup > 0) log_info(sprintf("  去重: 移除 %d 重复行 (%.1f%%)", n_dup, n_dup/n_before*100))

# 可用协变量检查
covs <- DAG11
covs <- covs[sapply(covs, function(v) {
  x <- dt_sub[[v]]
  if (is.factor(x))   return(nlevels(x) > 1)
  if (is.character(x)) return(length(unique(x)) > 1)
  return(sd(x, na.rm = TRUE) > 1e-8)
})]
log_info(sprintf("实际可用协变量 (%d): %s", length(covs), paste(covs, collapse = ", ")))

# 字符型→因子
for (v in covs) {
  if (is.character(dt_sub[[v]])) {
    dt_sub[, (v) := as.factor(get(v))]
  }
}

# 精简列 (避免 SEQTaRget 警告)
keep_cols <- unique(c("id_int", "period_int", "a_int", "y_int", "eligible_int", covs,
                      "surv_time", "event", "tx_int", "trial_num"))
extra_cols <- setdiff(names(dt_sub), keep_cols)
if (length(extra_cols) > 0) dt_sub[, (extra_cols) := NULL]

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
      trial.include   = TRUE,
      followup.include = TRUE,
      data.return     = FALSE
    ),
    verbose = FALSE
  )
}, error = function(e) {
  log_info(sprintf("[SEQ ERROR] %s", e$message))
  NULL
})
runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# ---- 回退方案: CoxPH ----
fallback <- NULL
if (is.null(result)) {
  log_info("SEQTaRget 失败，尝试 CoxPH 回退...")
  cov_formula <- paste(covs, collapse = " + ")
  f <- as.formula(sprintf("Surv(surv_time, event) ~ tx_int + %s", cov_formula))
  fb <- tryCatch({
    coxph(f, data = dt_sub, cluster = id_int)
  }, error = function(e) {
    log_info(sprintf("CoxPH 回退也失败: %s", e$message))
    NULL
  })
  if (!is.null(fb)) {
    s <- summary(fb)
    hr <- exp(coef(fb)["tx_int"])
    ci <- exp(confint(fb)["tx_int", ])
    p  <- s$coefficients["tx_int", "Pr(>|z|)"]
    fallback <- list(hr = hr, ci_lower = ci[1], ci_upper = ci[2], p = p)
    log_info(sprintf("CoxPH 回退: HR=%.4f (%.4f-%.4f), p=%.4f", hr, ci[1], ci[2], p))
  }
}

# ---- 提取 HR ----
if (!is.null(result)) {
  hr_vec <- hazard_ratio(result)
  hr_main     <- hr_vec[1]
  ci_lower    <- hr_vec[2]
  ci_upper    <- hr_vec[3]
  method_used <- "SEQTaRget_ITT"
  log_info(sprintf("SEQTaRget ITT: HR=%.4f (%.4f-%.4f), runtime=%.0fs",
                   hr_main, ci_lower, ci_upper, runtime))
} else if (!is.null(fallback)) {
  hr_main     <- fallback$hr
  ci_lower    <- fallback$ci_lower
  ci_upper    <- fallback$ci_upper
  method_used <- "CoxPH_cluster"
} else {
  hr_main <- NA_real_; ci_lower <- NA_real_; ci_upper <- NA_real_
  method_used <- "FAILED"
  log_info("[ERROR] 所有方法均失败")
}

# ============================================================================
# 4. 结果汇总与输出
# ============================================================================
log_section("结果保存")

# ---- 4a. 单个结果 CSV ----
result_single <- data.table(
  analysis      = "eICU External Validation (DAG11 CCA)",
  data          = "eICU-CRD",
  method        = method_used,
  n             = nrow(dt_cc),
  n_events      = sum(dt_cc$event),
  n_cov         = length(covs),
  cov_set       = "DAG11",
  hr            = hr_main,
  ci_lower      = ci_lower,
  ci_upper      = ci_upper,
  boot          = BOOT,
  runtime_sec   = runtime,
  missing_lactate_pct = round(sum(is.na(dt_raw$lactate_t0)) / nrow(dt_raw) * 100, 1),
  complete_case_pct   = round(nrow(dt_cc) / nrow(dt_raw) * 100, 1)
)
fwrite(result_single, file.path(TBL_DIR, "v26_9_eicu_dag11_result.csv"))
log_info(sprintf("结果已保存: v26_9_eicu_dag11_result.csv"))

# ---- 4b. DAG11 vs DAG14 对比表 ----
log_section("MIMIC vs eICU 对比")

# 读取 MIMIC DAG14 结果
dag14_file <- file.path(TBL_DIR, "v26_9_seqtarget_rubin_results.csv")
dag11_file <- file.path(TBL_DIR, "v26_9_s11_covariate_rubin.csv")

comparison_rows <- list()

# eICU DAG11 (当前结果)
comparison_rows[["eICU_DAG11"]] <- data.table(
  source    = "eICU-CRD",
  cov_set   = "DAG11",
  n         = nrow(dt_cc),
  n_events  = sum(dt_cc$event),
  hr        = sprintf("%.4f", hr_main),
  ci_lower  = sprintf("%.4f", ci_lower),
  ci_upper  = sprintf("%.4f", ci_upper),
  method    = method_used,
  notes     = "CCA (no MICE)"
)

# MIMIC DAG11 (S11)
if (file.exists(dag11_file)) {
  d11 <- fread(dag11_file)
  comparison_rows[["MIMIC_DAG11"]] <- data.table(
    source    = "MIMIC-IV",
    cov_set   = "DAG11",
    n         = d11$n[1],
    n_events  = d11$n_events[1],
    hr        = sprintf("%.4f", d11$hr[1]),
    ci_lower  = sprintf("%.4f", d11$ci_lower[1]),
    ci_upper  = sprintf("%.4f", d11$ci_upper[1]),
    method    = "SEQTaRget_ITT_Rubin×5",
    notes     = "MICE m=5, Rubin merged"
  )
} else {
  comparison_rows[["MIMIC_DAG11"]] <- data.table(
    source = "MIMIC-IV", cov_set = "DAG11",
    hr = "N/A", ci_lower = "N/A", ci_upper = "N/A",
    method = "NOT_RUN", notes = "Run S11 first"
  )
}

# MIMIC DAG14 (主分析)
if (file.exists(dag14_file)) {
  d14 <- fread(dag14_file)
  hr14_col <- grep("^hr$|^HR$", names(d14), value = TRUE, ignore.case = TRUE)[1]
  lo14_col <- grep("lower|lo", names(d14), value = TRUE, ignore.case = TRUE)[1]
  hi14_col <- grep("upper|hi", names(d14), value = TRUE, ignore.case = TRUE)[1]
  comparison_rows[["MIMIC_DAG14"]] <- data.table(
    source    = "MIMIC-IV",
    cov_set   = "DAG14",
    n         = NA_integer_,
    n_events  = NA_integer_,
    hr        = sprintf("%.4f", as.numeric(d14[[hr14_col]][1])),
    ci_lower  = sprintf("%.4f", as.numeric(d14[[lo14_col]][1])),
    ci_upper  = sprintf("%.4f", as.numeric(d14[[hi14_col]][1])),
    method    = "SEQTaRget_ITT_Rubin×20",
    notes     = "Primary analysis (MICE m=20)"
  )
}

dt_comp <- rbindlist(comparison_rows, use.names = TRUE, fill = TRUE)
fwrite(dt_comp, file.path(TBL_DIR, "v26_9_eicu_dag11_comparison.csv"))

# ---- 打印结果 ----
cat("\n")
cat("=", rep("=", 78), "\n", sep = "")
cat("  DAG11 eICU 外部验证结果\n")
cat("=", rep("=", 78), "\n", sep = "")
cat(sprintf("  eICU (DAG11 CCA):  HR=%.4f (%.4f–%.4f)\n", hr_main, ci_lower, ci_upper))
cat(sprintf("  方法: %s, n=%d, events=%d, boot=%d\n",
            method_used, nrow(dt_cc), sum(dt_cc$event), BOOT))
cat(sprintf("  运行时间: %.1f 分钟\n\n", runtime / 60))

log_info(sprintf("eICU DAG11: HR=%.4f (%.4f-%.4f), n=%d, events=%d, method=%s",
                 hr_main, ci_lower, ci_upper, nrow(dt_cc), sum(dt_cc$event), method_used))

# ============================================================================
# 5. 收尾
# ============================================================================
log_section("完成")

# 日志中输出对比
log_info("MIMIC vs eICU 对比:")
for (i in seq_len(nrow(dt_comp))) {
  r <- dt_comp[i]
  log_info(sprintf("  %s [%s]: HR=%s (%s-%s), method=%s",
                   r$source, r$cov_set, r$hr, r$ci_lower, r$ci_upper, r$method))
}

close(log_con)
cat(sprintf("\n✅ eICU DAG11 外部验证完成\n"))
cat(sprintf("   结果: %s\n", file.path(TBL_DIR, "v26_9_eicu_dag11_result.csv")))
cat(sprintf("   对比: %s\n", file.path(TBL_DIR, "v26_9_eicu_dag11_comparison.csv")))
cat(sprintf("   日志: %s\n", log_file))
