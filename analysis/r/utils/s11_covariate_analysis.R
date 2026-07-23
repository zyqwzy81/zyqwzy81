# ============================================================================
# Supplementary Analysis S16: 11-Covariate MIMIC-IV Reanalysis
# 使用 eICU 共享的 DAG11 协变量集重跑 SEQTaRget ITT
# 目的: 量化协变量减少（排除 sepsis3、corrected_calcium、nutrition_support_type）
#       对 HR 估计的影响
# 分析方案: SEQTaRget ITT × 5 MICE 插补集 × Rubin 合并 (DAG11, boot=50)
# ============================================================================
# DAG11 (eICU shared): age, gender, baseline_sofa, charlson, egfr,
#                       invasive_vent_at_t0, vasoactive_at_t0,
#                       admission_type, baseline_phosphate, lactate_t0,
#                       refeeding_risk
# DAG14 排除: sepsis3, corrected_calcium_mmol_min, nutrition_support_type
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(survival)
  library(SEQTaRget)
})

source("/data/小米云盘/文档/论文/血磷/TTE/v26/code/r/config.R")

# ---- 输出目录 ----
OUT_DIR     <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9"
RDS_DIR     <- file.path(OUT_DIR, "rds")
TBL_DIR     <- file.path(OUT_DIR, "tables")
FIG_DIR     <- file.path(OUT_DIR, "figures")
LOG_DIR     <- file.path(OUT_DIR, "logs")
for (d in c(OUT_DIR, RDS_DIR, TBL_DIR, FIG_DIR, LOG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- 数据路径 ----
DATA_FILE        <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_imputed_v26.parquet"
MICE_PARQUET_TEMPLATE <- file.path(
  "/data/小米云盘/文档/论文/血磷/sql/result/mice_imputed",
  "v26_mice_imputed_imp%02d.parquet"
)

# ---- 分析参数 ----
BOOT_S16  <- 50L       # 敏感性分析 boot=50 以节省时间
N_IMPUTE  <- 5L        # 与其他敏感性分析一致，使用 m=5
SEED      <- 20260629L
set.seed(SEED)
options(mc.cores = 1)

# ---- DAG11: 匹配 eICU 可用协变量 ----
DAG14_FULL <- c("age", "gender", "baseline_sofa", "sepsis3", "charlson",
                "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
                "admission_type", "baseline_phosphate", "lactate_t0",
                "refeeding_risk", "corrected_calcium_mmol_min",
                "nutrition_support_type")
DAG11 <- c("age", "gender", "baseline_sofa", "charlson",
           "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
           "admission_type", "baseline_phosphate", "lactate_t0",
           "refeeding_risk")
EXCLUDED_VARS <- setdiff(DAG14_FULL, DAG11)

# ---- 日志系统 ----
log_file <- file.path(LOG_DIR, sprintf("v26_9_s11_covariate_%s.log",
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

log_info("Supplementary Analysis S16: DAG11 MIMIC-IV Reanalysis (matching eICU)")
log_info(sprintf("种子: %d, Bootstrap: %d, 插补集: %d", SEED, BOOT_S16, N_IMPUTE))
log_info(sprintf("DAG11 变量 (%d 个): %s", length(DAG11),
                 paste(DAG11, collapse = ", ")))
log_info(sprintf("DAG14 排除变量 (%d 个): %s", length(EXCLUDED_VARS),
                 paste(EXCLUDED_VARS, collapse = ", ")))
log_info(sprintf("输出目录: %s", OUT_DIR))

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
log_info(sprintf("主 Parquet: %d 行, %d 事件, %d 列",
                 nrow(dt_full), sum(dt_full$event), ncol(dt_full)))

key_cols <- c("stay_id", "subject_id", "trial_num", "time_zero",
              "itt_group", "gender", "admission_type",
              "surv_time", "death_time_from_t0", "dischtime_from_t0",
              "race_group", "icu_type",
              "death_28day_num", "nutrition_support_type")
key_cols <- intersect(key_cols, names(dt_full))
dt_key   <- dt_full[, ..key_cols]

# ---- DAG11 可用性检查 ----
mice_p1 <- as.data.table(read_parquet(sprintf(MICE_PARQUET_TEMPLATE, 1L)))
dag11_in_mice  <- intersect(DAG11, names(mice_p1))
dag11_not_mice <- setdiff(DAG11, names(mice_p1))
log_info(sprintf("DAG11 在 MICE Parquet 中 (%d): %s",
                 length(dag11_in_mice), paste(dag11_in_mice, collapse = ", ")))
if (length(dag11_not_mice) > 0) {
  log_info(sprintf("DAG11 从主 Parquet 获取 (%d): %s",
                   length(dag11_not_mice), paste(dag11_not_mice, collapse = ", ")))
}
rm(mice_p1); gc(TRUE)

# ============================================================================
# 2. SEQTaRget ITT 封装函数 (DAG11)
# ============================================================================
run_seqtarget_dag11 <- function(dt_comb, dag11_covs, boot = BOOT_S16) {
  dt_sub <- dt_comb

  # 基础列
  # death_28day_num 在 MICE Parquet 中为字符串/categorical, as.numeric(factor) 会返回 1/2 而非 0/1
  y_vals <- dt_sub[["death_28day_num"]]
  if (is.factor(y_vals) || is.character(y_vals)) {
    y_int_vals <- as.integer(as.character(y_vals))
  } else {
    y_int_vals <- as.integer(y_vals)
  }
  y_int_vals <- pmax(0L, pmin(1L, y_int_vals, na.rm = TRUE))
  dt_sub[, `:=`(
    id_int       = as.integer(factor(stay_id)),
    period_int   = as.integer(trial_num) - 1L,
    a_int        = as.integer(itt_group == "treatment"),
    y_int        = y_int_vals,
    eligible_int = 1L
  )]
  rm(y_vals, y_int_vals)
  setorder(dt_sub, id_int, period_int)

  # 过滤每期：需同时有治疗/对照，且至少5个事件
  dt_sub[, `:=`(n_a = sum(a_int), n_y = sum(y_int)), by = period_int]
  dt_sub <- dt_sub[n_a > 0 & n_a < .N & n_y >= 5]
  if (nrow(dt_sub) < 100) return(NULL)

  # 可用协变量检查（剔除零方差/单水平因子）
  covs <- intersect(dag11_covs, names(dt_sub))
  covs <- covs[sapply(covs, function(v) {
    x <- dt_sub[[v]]
    if (is.factor(x))   return(nlevels(x) > 1)
    if (is.character(x)) return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]
  if (length(covs) < 2) return(NULL)

  # 字符型→因子
  for (v in covs) {
    if (is.character(dt_sub[[v]])) {
      dt_sub[, (v) := as.factor(get(v))]
    }
  }

  keep_cols <- unique(c("id_int","period_int","a_int","y_int","eligible_int", covs))
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
        bootstrap.nboot = boot,
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
    log_info(sprintf("  [SEQ ERROR] imp=%s: %s",
                     tryCatch(imp_i, error = function(e2) "?"), e$message))
    NULL
  })
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (is.null(result)) return(NULL)

  hr_vec <- hazard_ratio(result)
  out <- list(
    hr        = hr_vec[1],
    ci_lower  = hr_vec[2],
    ci_upper  = hr_vec[3],
    runtime   = runtime,
    n         = nrow(dt_sub),
    n_events  = sum(dt_sub$event, na.rm = TRUE)
  )
  rm(result, dt_sub); gc(TRUE)
  return(out)
}

# ============================================================================
# 3. Rubin 规则合并
# ============================================================================
rubin_combine <- function(betas, ses) {
  m <- length(betas)
  beta_bar <- mean(betas)
  W <- mean(ses^2)
  B <- if (m > 1) sum((betas - beta_bar)^2) / (m - 1) else 0
  se_pool <- sqrt(W + (1 + 1 / m) * B)
  t_stat  <- beta_bar / se_pool
  df <- if (m > 1) {
    (m - 1) * (1 + W / ((1 + 1 / m) * B))^2
  } else {
    1
  }
  p <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)
  list(
    beta      = beta_bar,
    se        = se_pool,
    t         = t_stat,
    df        = df,
    p         = p,
    hr        = exp(beta_bar),
    ci_lower  = exp(beta_bar - 1.96 * se_pool),
    ci_upper  = exp(beta_bar + 1.96 * se_pool),
    W = W, B = B
  )
}

# ============================================================================
# 4. 主循环：SEQTaRget ITT ×5 插补 (DAG11)
# ============================================================================
log_section("SEQTaRget DAG11 ×5 插补遍历 (boot=50, ncores=4)")

s16_results <- list()
success_count <- 0

for (imp_i in 1L:N_IMPUTE) {
  log_info(sprintf("--- DAG11 插补集 %02d/%02d ---", imp_i, N_IMPUTE))

  mice_file <- sprintf(MICE_PARQUET_TEMPLATE, imp_i)
  if (!file.exists(mice_file)) {
    log_info(sprintf("[WARN] MICE Parquet 不存在: %s", mice_file))
    next
  }
  mice_dt <- as.data.table(read_parquet(mice_file))

  # Merge 对齐
  by_cols <- c("stay_id", "trial_num")
  overlapped_remove <- intersect(names(mice_dt), names(dt_key))
  overlapped_remove <- setdiff(overlapped_remove, by_cols)
  dt_key_sub <- dt_key[, setdiff(names(dt_key), overlapped_remove), with = FALSE]
  dt_comb <- merge(dt_key_sub, mice_dt, by = by_cols, all.x = TRUE, sort = FALSE)
  rm(mice_dt, dt_key_sub)

  stopifnot("itt_group" %in% names(dt_comb),
            "death_28day_num" %in% names(dt_comb))

  res <- run_seqtarget_dag11(dt_comb, DAG11, boot = BOOT_S16)
  rm(dt_comb); gc(TRUE)

  if (is.null(res)) {
    log_info(sprintf("[WARN] 插补 %02d DAG11 分析返回 NULL", imp_i))
    next
  }

  success_count <- success_count + 1
  res$imp <- imp_i
  s16_results[[success_count]] <- res
  log_info(sprintf("  插补 %02d: HR=%.4f, CI=%.4f-%.4f, runtime=%.0fs, n=%d, events=%d",
                   imp_i, res$hr, res$ci_lower, res$ci_upper,
                   res$runtime, res$n, res$n_events))
}

log_section("DAG11 分析汇总")
log_info(sprintf("成功插补数: %d/%d", success_count, N_IMPUTE))

if (success_count == 0) {
  log_info("[ERROR] 所有 DAG11 分析均失败！")
  close(log_con)
  stop("DAG11 分析全部失败")
}

# ---- Rubin 合并 ----
betas <- sapply(s16_results, function(x) log(x$hr))
ses   <- sapply(s16_results, function(x) (log(x$ci_upper) - log(x$ci_lower)) / (2 * 1.96))

rubin <- rubin_combine(betas, ses)
log_info(sprintf("Rubin 合并结果 (DAG11):"))
log_info(sprintf("  HR = %.4f (%.4f-%.4f), p = %.4f",
                 rubin$hr, rubin$ci_lower, rubin$ci_upper, rubin$p))
log_info(sprintf("  beta = %.6f, SE = %.6f, W = %.6f, B = %.6f",
                 rubin$beta, rubin$se, rubin$W, rubin$B))
log_info(sprintf("  HR 范围: %.4f-%.4f, SD = %.6f",
                 min(sapply(s16_results, function(x) x$hr)),
                 max(sapply(s16_results, function(x) x$hr)),
                 sd(sapply(s16_results, function(x) x$hr))))

# ============================================================================
# 5. 导出结果
# ============================================================================
log_section("结果导出")

# Per-imp 结果
s16_per_imp <- data.table(
  imp      = sapply(s16_results, function(x) x$imp),
  hr       = sapply(s16_results, function(x) x$hr),
  ci_lower = sapply(s16_results, function(x) x$ci_lower),
  ci_upper = sapply(s16_results, function(x) x$ci_upper),
  runtime  = sapply(s16_results, function(x) x$runtime),
  n        = sapply(s16_results, function(x) x$n),
  n_events = sapply(s16_results, function(x) x$n_events)
)
fwrite(s16_per_imp, file.path(TBL_DIR, "v26_9_s11_covariate_per_imp.csv"))

# Rubin 合并结果
s16_rubin <- data.table(
  analysis = "S16: DAG11 MIMIC-IV (matching eICU)",
  cov_set  = "DAG11",
  n_cov    = length(DAG11),
  hr       = rubin$hr,
  ci_lower = rubin$ci_lower,
  ci_upper = rubin$ci_upper,
  p        = rubin$p,
  n        = s16_results[[1]]$n,
  n_events = s16_results[[1]]$n_events,
  W        = rubin$W,
  B        = rubin$B,
  fmi      = (1 + 1/success_count) * rubin$B / (rubin$W + (1 + 1/success_count) * rubin$B),
  boot     = BOOT_S16,
  n_imp    = success_count,
  excluded_vars = paste(EXCLUDED_VARS, collapse = ";")
)
fwrite(s16_rubin, file.path(TBL_DIR, "v26_9_s11_covariate_rubin.csv"))

# ============================================================================
# 6. 生成对比表：DAG11 vs DAG14 (参考)
# ============================================================================
log_section("DAG11 vs DAG14 对比")

# 读取已有的 DAG14 Rubin 结果
dag14_rubin_file <- file.path(TBL_DIR, "v26_9_seqtarget_rubin_results.csv")
if (file.exists(dag14_rubin_file)) {
  dag14_rubin <- fread(dag14_rubin_file)
  # 兼容不同列名
  hr14_col <- if ("hr" %in% names(dag14_rubin)) "hr" else
    grep("^hr$|^HR$", names(dag14_rubin), value = TRUE, ignore.case = TRUE)[1]
  lo14_col <- if ("ci_lower" %in% names(dag14_rubin)) "ci_lower" else
    grep("lower|lo|ci_l", names(dag14_rubin), value = TRUE, ignore.case = TRUE)[1]
  hi14_col <- if ("ci_upper" %in% names(dag14_rubin)) "ci_upper" else
    grep("upper|hi|ci_u", names(dag14_rubin), value = TRUE, ignore.case = TRUE)[1]

  if (length(hr14_col) == 1 && length(lo14_col) == 1 && length(hi14_col) == 1) {
    hr14  <- dag14_rubin[[hr14_col]][1]
    lo14  <- dag14_rubin[[lo14_col]][1]
    hi14  <- dag14_rubin[[hi14_col]][1]

    # 计算相对差异
    rel_diff <- (rubin$hr - hr14) / hr14 * 100

    # DAG13 结果（已有）
    dag13_file <- file.path(TBL_DIR, "v26_9_sensitivity_dag13_comparison.csv")
    dag13_hr <- NA; dag13_lo <- NA; dag13_hi <- NA
    if (file.exists(dag13_file)) {
      dag13_dt <- fread(dag13_file)
      # 找 DAG13 行
      dag13_row <- dag13_dt[grepl("DAG13", .SD[[1]], ignore.case = TRUE)]
      if (nrow(dag13_row) > 0) {
        dag13_hr <- dag13_row[[grep("hr|HR", names(dag13_row), value = TRUE, ignore.case = TRUE)[1]]]
        dag13_lo <- dag13_row[[grep("lower|lo", names(dag13_row), value = TRUE, ignore.case = TRUE)[1]]]
        dag13_hi <- dag13_row[[grep("upper|hi", names(dag13_row), value = TRUE, ignore.case = TRUE)[1]]]
      }
    }

    # 三组对比表
    comparison <- data.table(
      covariate_set = c("DAG14 (Primary)", "DAG13 (excl. sepsis3)", "DAG11 (eICU shared)"),
      n_covariates  = c(length(DAG14_FULL), length(DAG14_FULL) - 1, length(DAG11)),
      hr            = c(sprintf("%.3f", hr14),
                        if (is.na(dag13_hr)) "N/A" else sprintf("%.3f", dag13_hr),
                        sprintf("%.3f", rubin$hr)),
      ci_lower      = c(sprintf("%.3f", lo14),
                        if (is.na(dag13_lo)) "N/A" else sprintf("%.3f", dag13_lo),
                        sprintf("%.3f", rubin$ci_lower)),
      ci_upper      = c(sprintf("%.3f", hi14),
                        if (is.na(dag13_hi)) "N/A" else sprintf("%.3f", dag13_hi),
                        sprintf("%.3f", rubin$ci_upper)),
      boot          = c("200", "200", sprintf("%d", BOOT_S16)),
      n_impute      = c("20", "1 (imp01)", sprintf("%d", success_count))
    )

    log_info(sprintf("DAG11 vs DAG14 对比:"))
    log_info(sprintf("  DAG14 (参考): HR=%.3f (%.3f-%.3f)", hr14, lo14, hi14))
    if (!is.na(dag13_hr)) {
      log_info(sprintf("  DAG13:         HR=%.3f (%.3f-%.3f)", dag13_hr, dag13_lo, dag13_hi))
    }
    log_info(sprintf("  DAG11 (S16):   HR=%.3f (%.3f-%.3f)",
                     rubin$hr, rubin$ci_lower, rubin$ci_upper))
    log_info(sprintf("  相对差异 (DAG11 vs DAG14): %.1f%%", rel_diff))

    fwrite(comparison, file.path(TBL_DIR, "v26_9_s11_covariate_comparison.csv"))
    log_info("  对比表已导出: v26_9_s11_covariate_comparison.csv")
  } else {
    log_info("[WARN] DAG14 Rubin 结果列名不匹配，跳过对比")
  }
} else {
  log_info("[WARN] DAG14 Rubin 结果文件不存在，跳过对比")
}

# ============================================================================
# 7. 完成
# ============================================================================
log_section("S16 分析完成")
log_info(sprintf("DAG11 (eICU shared, %d covariates):", length(DAG11)))
log_info(sprintf("  Rubin HR=%.4f, 95%% CI (%.4f-%.4f), p=%.4f",
                 rubin$hr, rubin$ci_lower, rubin$ci_upper, rubin$p))
log_info(sprintf("  各插补 HR 范围: %.4f-%.4f",
                 min(s16_per_imp$hr), max(s16_per_imp$hr)))
log_info("完整输出路径:")
log_info(sprintf("  Per-imp:   %s", file.path(TBL_DIR, "v26_9_s11_covariate_per_imp.csv")))
log_info(sprintf("  Rubin:     %s", file.path(TBL_DIR, "v26_9_s11_covariate_rubin.csv")))
log_info(sprintf("  Comparison: %s", file.path(TBL_DIR, "v26_9_s11_covariate_comparison.csv")))

close(log_con)
cat("\n✅ S16 DAG11 分析完成!\n")
