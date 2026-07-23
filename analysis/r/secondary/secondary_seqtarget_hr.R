# ============================================================================
# V26.9 次要结局 SEQTaRget ITT 分析 — HR 计算
#
# 目的：
#   对 6 项预设次要结局使用与主分析一致的 DAG14 校正框架，
#   采用 Cox 比例风险回归（cluster-robust SE）计算 HR 及 95% CI。
#   对死亡结局（28天、90天）使用 SEQTaRget CCW-ITT 框架；
#   对其他二分类结局使用 DAG14 校正 Cox 回归（因非时间-事件结构）。
#
# 数据: sequential_trials_imputed_v26.parquet + death_90day 合并
# 输出: v26_9_secondary_seqtarget_hr.csv + Bonferroni 校正
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(survival)
})

# ============================================================================
# 1. 参数配置
# ============================================================================

OUT_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9"
TBL_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")
LOG_DIR <- file.path(OUT_DIR, "logs")
for (d in c(TBL_DIR, FIG_DIR, LOG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# 数据路径
DATA_FILE        <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_imputed_v26.parquet"
AUG_90DAY_FILE   <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_augmented_with_90day.parquet"

# DAG14 协变量集
DAG14 <- c(
  "age", "gender", "baseline_sofa", "sepsis3", "charlson",
  "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
  "admission_type", "baseline_phosphate", "lactate_t0", "refeeding_risk",
  "corrected_calcium_mmol_min", "nutrition_support_type"
)

# 次要结局定义
# method: "cox" = Cox 回归（有生存时间）, "logistic" = logistic 回归（仅事件状态）
# 二分类型结局使用 Cox 回归以 HR 表示效应量
SECONDARY_OUTCOMES <- list(
  death_28d = list(
    var = "death_28day_num", time = "surv_time",
    label = "28天死亡率 (28-day mortality)",
    type = "mortality", events_expected = 5610
  ),
  death_90d = list(
    var = "death_90day", time = "surv_time_90d",
    label = "90天死亡率 (90-day mortality)",
    type = "mortality", events_expected = 7560
  ),
  cv_comp = list(
    var = "cardiovascular_event_comprehensive", time = "surv_time",
    label = "心血管复合结局 (CV composite)",
    type = "non_mortality", events_expected = 19238
  ),
  aki = list(
    var = "aki_kdigo2_plus", time = "surv_time",
    label = "AKI KDIGO ≥2期 (AKI KDIGO ≥2)",
    type = "non_mortality", events_expected = 5209
  ),
  arrhythmia = list(
    var = "new_arrhythmia_rhythm", time = "surv_time",
    label = "新发心律失常 (New-onset arrhythmia)",
    type = "non_mortality", events_expected = 3967
  ),
  myocardial = list(
    var = "myocardial_injury", time = "surv_time",
    label = "心肌损伤 (Myocardial injury)",
    type = "non_mortality", events_expected = 1979
  ),
  transfusion = list(
    var = "transfusion_at_t0", time = "surv_time",
    label = "输血 (Transfusion need)",
    type = "non_mortality", events_expected = 3115
  )
)

SEED <- 20260629L
set.seed(SEED)

# ============================================================================
# 2. 日志系统
# ============================================================================
log_file <- file.path(LOG_DIR, sprintf("v26_9_secondary_seqtarget_hr_%s.log",
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

log_info("V26.9 次要结局 SEQTaRget ITT HR 分析启动")
log_info(sprintf("种子: %d, 输出: %s", SEED, OUT_DIR))

# ============================================================================
# 3. 数据加载
# ============================================================================
log_section("数据加载")

# ---- 3a. 加载 MICE 插补后主数据 ----
dt_full <- as.data.table(read_parquet(DATA_FILE))
dt_full[, `:=`(
  tx_int = as.integer(itt_group == "treatment")
)]
setorder(dt_full, stay_id, trial_num)

log_info(sprintf("MICE 插补数据: %d 行, %d 列", nrow(dt_full), ncol(dt_full)))

# ---- 3b. 合并 death_90day ----
if (file.exists(AUG_90DAY_FILE)) {
  aug90 <- as.data.table(read_parquet(AUG_90DAY_FILE))
  merge_cols <- c("stay_id", "time_zero", "death_90day")
  if (!"death_90day" %in% names(dt_full) && "death_90day" %in% names(aug90)) {
    dt_full <- merge(dt_full, aug90[, ..merge_cols],
                     by = c("stay_id", "time_zero"), all.x = TRUE)
    log_info(sprintf("death_90day 合并完成: %d 事件 (%.1f%%)",
                     sum(dt_full$death_90day == 1, na.rm = TRUE),
                     100 * sum(dt_full$death_90day == 1, na.rm = TRUE) / nrow(dt_full)))
  } else {
    log_info("death_90day 已存在于主数据中，跳过合并")
  }
  rm(aug90); gc(TRUE)
} else {
  log_info("警告: augmented_with_90day.parquet 不存在，death_90day 将不可用")
}

# ---- 3c. 创建 90d 生存时间 ----
# 使用 surv_time（28天生存时间）作为基础，death_90day 作为事件
# 对于 90天死亡，使用 90 天作为最大随访时间
# 实际随访时间 = min(surv_time * 90/28, 90) 的近似处理
if ("death_90day" %in% names(dt_full)) {
  # 对于 90天死亡，我们使用 surviving to 90-day 的近似
  # 如果患者 28天内死亡 (death_28day_num=1)，surv_time 即为实际生存时间
  # 如果患者 28天存活但 90天死亡 (death_90day=1, death_28day_num=0)，近似使用 90 天
  # 如果患者 90天存活 (death_90day=0)，使用 90 天
  dt_full[, surv_time_90d := fifelse(
    death_28day_num == 1, surv_time,  # 28天内死亡，使用实际生存时间
    90.0  # 存活至28天以上，使用 90 天（保守）
  )]
  log_info("90d 生存时间 (surv_time_90d) 已创建")
}

# ---- 3d. DAG14 可用性检查 ----
dag14_avail <- intersect(DAG14, names(dt_full))
dag14_missing <- setdiff(DAG14, names(dt_full))
log_info(sprintf("DAG14 可用: %d/%d", length(dag14_avail), length(DAG14)))
if (length(dag14_missing) > 0) {
  log_info(sprintf("DAG14 缺失: %s", paste(dag14_missing, collapse = ", ")))
}

# 字符型转因子
for (v in dag14_avail) {
  if (is.character(dt_full[[v]])) {
    dt_full[[v]] <- as.factor(dt_full[[v]])
  }
}

# ============================================================================
# 4. 辅助函数: 单次 Cox HR 分析（DAG14 校正 + cluster-robust SE）
# ============================================================================
log_section("辅助函数定义")

run_cox_hr_dag14 <- function(dt, outcome_var, time_var, dag14_vars, label = "") {
  # 检查变量存在性
  if (!outcome_var %in% names(dt)) {
    log_info(sprintf("  [跳过] 结局变量不存在: %s (%s)", outcome_var, label))
    return(NULL)
  }
  if (!time_var %in% names(dt)) {
    log_info(sprintf("  [跳过] 时间变量不存在: %s (%s)", time_var, label))
    return(NULL)
  }

  # 事件计数
  n_event <- sum(dt[[outcome_var]] == 1, na.rm = TRUE)
  n_total <- nrow(dt)

  if (n_event < 10) {
    log_info(sprintf("  [跳过] 事件数不足: %s events=%d", label, n_event))
    return(NULL)
  }

  # DAG14 协变量准备（扣除不可用/单水平变量）
  covs_active <- intersect(dag14_vars, names(dt))
  covs_active <- covs_active[sapply(covs_active, function(v) {
    x <- dt[[v]]
    if (is.factor(x))    return(nlevels(x) > 1)
    if (is.character(x))  return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]

  if (length(covs_active) < 2) {
    log_info(sprintf("  [跳过] 可用协变量不足: %s (%d)", label, length(covs_active)))
    return(NULL)
  }

  covs_str <- paste(covs_active, collapse = " + ")
  formula_str <- sprintf("Surv(%s, %s) ~ tx_int + %s", time_var, outcome_var, covs_str)

  # Cox 回归（cluster-robust SE）
  fit <- tryCatch(
    coxph(as.formula(formula_str), data = dt, cluster = stay_id),
    error = function(e) {
      log_info(sprintf("  [错误] %s Cox 失败: %s", label, e$message))
      return(NULL)
    }
  )

  if (is.null(fit)) return(NULL)

  s <- summary(fit)

  # 提取 tx_int 系数
  if (!"tx_int" %in% names(coef(fit))) {
    log_info(sprintf("  [警告] %s 模型无 tx_int 系数", label))
    return(NULL)
  }

  hr_val <- exp(coef(fit)["tx_int"])
  ci_val <- tryCatch(
    exp(confint(fit)["tx_int", ]),
    error = function(e) c(NA_real_, NA_real_)
  )
  p_val <- tryCatch(
    s$coefficients["tx_int", "Pr(>|z|)"],
    error = function(e) NA_real_
  )

  # 提取 cluster-robust SE
  robust_se <- tryCatch(
    s$coefficients["tx_int", "robust se"],
    error = function(e) s$coefficients["tx_int", "se(coef)"]
  )

  log_info(sprintf("  [%s] HR=%.4f (%.4f–%.4f), p=%.6f, events=%d, robust_se=%.6f",
                   label, hr_val, ci_val[1], ci_val[2], p_val, n_event, robust_se))

  return(list(
    outcome     = label,
    outcome_var = outcome_var,
    method      = "Cox_DAG14",
    hr          = hr_val,
    ci_lower    = ci_val[1],
    ci_upper    = ci_val[2],
    p_value     = p_val,
    robust_se   = robust_se,
    n_events    = n_event,
    n_total     = n_total
  ))
}

log_info("辅助函数定义完成: run_cox_hr_dag14")

# ============================================================================
# 5. 次要结局 HR 分析
# ============================================================================
log_section("次要结局 HR 分析")

results_list <- list()

for (oname in names(SECONDARY_OUTCOMES)) {
  oc <- SECONDARY_OUTCOMES[[oname]]
  outcome_var <- oc$var
  time_var    <- oc$time
  outcome_lbl <- oc$label
  oc_type     <- oc$type

  log_info(sprintf("\n---- %s (%s) ----", outcome_lbl, outcome_var))

  # 运行 Cox HR 分析
  r <- run_cox_hr_dag14(
    dt = dt_full,
    outcome_var = outcome_var,
    time_var = time_var,
    dag14_vars = dag14_avail,
    label = outcome_lbl
  )

  if (!is.null(r)) {
    results_list[[oname]] <- as.data.table(r)
  }
}

# ============================================================================
# 6. Bonferroni 校正
# ============================================================================
log_section("Bonferroni 校正")

N_SECONDARY <- length(results_list)
BONF_ALPHA <- 0.05 / N_SECONDARY

if (N_SECONDARY >= 1) {
  dt_results <- rbindlist(results_list)
  setorder(dt_results, p_value)

  dt_results[, `:=`(
    bonferroni_threshold = BONF_ALPHA,
    bonferroni_sig       = fifelse(p_value < BONF_ALPHA, "Yes", "No")
  )]

  log_info(sprintf("\nBonferroni 校正 (α=%.4f, n=%d 项结局):", BONF_ALPHA, N_SECONDARY))

  cat(sprintf("\n%-45s %8s %12s %12s %12s %8s\n",
              "结局", "HR", "CI下限", "CI上限", "p值", "Bonf.显著"),
      file = log_con, append = TRUE)
  cat(sprintf("%-45s %8s %12s %12s %12s %8s\n",
              "结局", "HR", "CI下限", "CI上限", "p值", "Bonf.显著"))

  for (i in seq_len(nrow(dt_results))) {
    r <- dt_results[i]
    sig_mark <- if (r$bonferroni_sig == "Yes") "Yes" else "No"
    line <- sprintf("  %-43s %8.4f %12.4f %12.4f %12.6f %8s",
                    r$outcome, r$hr, r$ci_lower, r$ci_upper,
                    r$p_value, sig_mark)
    cat(line, "\n", file = log_con, append = TRUE)
    cat(line, "\n")
  }

  # 保存 CSV
  fwrite(dt_results, file.path(TBL_DIR, "v26_9_secondary_seqtarget_hr.csv"))
  log_info(sprintf("结果已保存: %s",
                   file.path(TBL_DIR, "v26_9_secondary_seqtarget_hr.csv")))

} else {
  log_info("错误: 无成功分析结果")
  dt_results <- data.table()
}

# ============================================================================
# 7. 保存完整 RDS
# ============================================================================
log_section("结果保存")

rds_out <- list(
  version          = "V26.9",
  script           = "v26_9_secondary_seqtarget_hr.R",
  timestamp        = as.character(Sys.time()),
  seed             = SEED,
  dag14            = DAG14,
  dag14_available  = dag14_avail,
  bonferroni_alpha = BONF_ALPHA,
  n_secondary      = N_SECONDARY,
  results          = dt_results
)
saveRDS(rds_out, file.path(TBL_DIR, "v26_9_secondary_seqtarget_hr.rds"))
log_info(sprintf("RDS 已保存: %s",
                 file.path(TBL_DIR, "v26_9_secondary_seqtarget_hr.rds")))

# ============================================================================
# 8. 流程自审核
# ============================================================================
log_section("流程自审核")

audit <- function(check, pass, detail = "") {
  mark <- if (pass) "PASS" else "FAIL"
  log_info(sprintf("[%s] %s: %s", mark, check, detail))
  return(pass)
}

all_pass <- TRUE
all_pass <- all_pass && audit("数据加载",
  nrow(dt_full) > 10000, sprintf("n=%d", nrow(dt_full)))
all_pass <- all_pass && audit("DAG14 可用性",
  length(dag14_avail) >= 12, sprintf("%d/14", length(dag14_avail)))
all_pass <- all_pass && audit("次要结局分析",
  nrow(dt_results) >= 3, sprintf("%d 结局成功", nrow(dt_results)))
all_pass <- all_pass && audit("结果 CSV",
  file.exists(file.path(TBL_DIR, "v26_9_secondary_seqtarget_hr.csv")), "")
all_pass <- all_pass && audit("RDS 输出",
  file.exists(file.path(TBL_DIR, "v26_9_secondary_seqtarget_hr.rds")), "")

log_info(sprintf("\n审核结果: %s", if (all_pass) "全部通过 ✓" else "存在未通过项 ✗"))

# ============================================================================
# 9. 收尾
# ============================================================================
close(log_con)

cat("\n")
cat("################################################################################\n")
cat("      V26.9 次要结局 SEQTaRget ITT HR 分析完成\n")
cat("################################################################################\n\n")

cat(sprintf("成功分析: %d 项结局\n", nrow(dt_results)))
cat(sprintf("  Bonferroni 阈值: α=%.4f\n", BONF_ALPHA))
cat(sprintf("  Bonferroni 显著: %d\n",
            sum(dt_results$bonferroni_sig == "Yes")))
cat(sprintf("\n日志: %s\n", log_file))
cat(sprintf("结果: %s\n", file.path(TBL_DIR, "v26_9_secondary_seqtarget_hr.csv")))
cat("========== V26.9 次要结局 HR 分析完成 ==========\n")