# ============================================================================
# V26.9 eICU Sepsis3 缺失敏感性分析
#
# 审稿要求 M5: 在 MIMIC-IV 全数据上，排除 sepsis3 后重新运行 SEQTaRget ITT
# 比较 DAG14（完整）vs DAG13（排除 sepsis3）的 HR
#
# 分析设计:
#   - 非脓毒症亚组: sepsis3 == 0, 使用 DAG13 协变量集
#   - 脓毒症亚组:   sepsis3 == 1, 使用 DAG13 协变量集
#
# 数据: sequential_trials_imputed_v26.parquet (MIMIC-IV)
# 输出: TTE/v26/output/v26_9/
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(survival)
  library(SEQTaRget)
  library(jsonlite)
})

# ============================================================================
# 1. 加载配置 + 参数定义
# ============================================================================
source("/data/小米云盘/文档/论文/血磷/TTE/v26/code/r/config.R")

# ---- 输出目录 ----
OUT_DIR     <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9"
RDS_DIR     <- file.path(OUT_DIR, "rds")
TBL_DIR     <- file.path(OUT_DIR, "tables")
LOG_DIR     <- file.path(OUT_DIR, "logs")
for (d in c(OUT_DIR, RDS_DIR, TBL_DIR, LOG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- 数据路径 ----
DATA_FILE   <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_imputed_v26.parquet"

# ---- 分析参数 ----
BOOT_SUB    <- 50L       # 亚组分析 bootstrap 次数
SEED        <- 20260629L # 固定种子
set.seed(SEED)

# ---- DAG13 协变量集 = DAG14 minus "sepsis3" ----
DAG13 <- setdiff(DAG14, "sepsis3")

# ============================================================================
# 2. 日志系统
# ============================================================================
log_file <- file.path(LOG_DIR, sprintf("v26_9_eicu_sepsis3_sensitivity_%s.log",
                                        format(Sys.time(), "%Y%m%d_%H%M%S")))
log_con  <- file(log_file, open = "wt")

log_info <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg),
      file = log_con, append = TRUE)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}

log_section <- function(title) {
  sep <- paste(rep("=", 70), collapse = "")
  cat(sprintf("\n%s\n  %s\n%s\n", sep, title, sep),
      file = log_con, append = TRUE)
  cat(sprintf("\n%s\n  %s\n%s\n", sep, title, sep))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

log_info("V26.9 eICU Sepsis3 缺失敏感性分析启动")
log_info(sprintf("种子: %d, Bootstrap: %d, ncores: 1", SEED, BOOT_SUB))
log_info(sprintf("DAG14 (%d 变量): %s", length(DAG14), paste(DAG14, collapse = ", ")))
log_info(sprintf("DAG13 (%d 变量): %s", length(DAG13), paste(DAG13, collapse = ", ")))

# ============================================================================
# 3. 辅助函数
# ============================================================================
log_section("辅助函数定义")

# ---- 3a. SEQTaRget ITT 封装（通用版） ----
run_seqtarget_subgroup <- function(dt, label, boot = BOOT_SUB, cov_set = NULL) {
  dt_sub <- copy(dt)
  dt_sub[, `:=`(
    id_int       = as.integer(factor(stay_id)),
    period_int   = as.integer(trial_num) - 1L,
    a_int        = as.integer(tx_int),
    y_int        = as.integer(event),
    eligible_int = 1L
  )]
  setorder(dt_sub, id_int, period_int)

  # 剔除不可识别期：需同时有治疗/对照组、且至少 5 个事件
  dt_sub[, `:=`(n_a = sum(a_int), n_y = sum(y_int)), by = period_int]
  dt_sub <- dt_sub[n_a > 0 & n_a < .N & n_y >= 5]

  # 可用协变量检查（在当前子集中）
  use_cov <- intersect(cov_set, names(dt_sub))
  use_cov <- use_cov[sapply(use_cov, function(v) {
    x <- dt_sub[[v]]
    if (is.factor(x))   return(nlevels(x) > 1)
    if (is.character(x)) return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]

  if (nrow(dt_sub) < 100 || length(use_cov) < 2) {
    log_info(sprintf("[%s] 样本量或协变量不足 (n=%d, cov=%d)", label, nrow(dt_sub), length(use_cov)))
    return(NULL)
  }

  # 字符型变量转为 factor（SEQTaRget 要求）
  for (v in use_cov) {
    if (is.character(dt_sub[[v]])) {
      dt_sub[, (v) := as.factor(get(v))]
    }
  }

  t0 <- Sys.time()
  result <- tryCatch({
    SEQuential(
      data = as.data.frame(dt_sub),
      id.col = "id_int", time.col = "period_int",
      eligible.col = "eligible_int", treatment.col = "a_int",
      outcome.col = "y_int",
      fixed.cols = use_cov, time_varying.cols = character(0),
      method = "ITT",
      options = SEQopts(
        bootstrap = TRUE, bootstrap.nboot = boot, bootstrap.CI = 0.95,
        ncores = 1L, fastglm.method = 3L, hazard = TRUE,
        trial.include = TRUE, followup.include = TRUE,
        data.return = FALSE
      ),
      verbose = FALSE
    )
  }, error = function(e) {
    log_info(sprintf("[%s] SEQTaRget 错误: %s", label, e$message))
    return(NULL)
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
    n_events  = sum(dt_sub$event, na.rm = TRUE),
    n_cov     = length(use_cov),
    covariates = use_cov
  )
  rm(result, dt_sub); gc()
  return(out)
}

# ---- 3b. 回退 Cox PH 估计（SEQTaRget 失败时） ----
fallback_coxph <- function(dt, label, cov_set = NULL) {
  dt_sub <- copy(dt)
  dt_sub[, `:=`(
    id_int     = as.integer(factor(stay_id)),
    period_int = as.integer(trial_num) - 1L
  )]
  setorder(dt_sub, id_int, period_int)

  # 剔除不可识别期
  dt_sub[, `:=`(n_a = sum(tx_int), n_y = sum(event)), by = period_int]
  dt_sub <- dt_sub[n_a > 0 & n_a < .N & n_y >= 5]

  # 可用协变量检查
  use_cov <- intersect(cov_set, names(dt_sub))
  use_cov <- use_cov[sapply(use_cov, function(v) {
    x <- dt_sub[[v]]
    if (is.factor(x))   return(nlevels(x) > 1)
    if (is.character(x)) return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]

  if (nrow(dt_sub) < 100 || length(use_cov) < 2) return(NULL)

  # 构建公式
  cov_formula <- paste(use_cov, collapse = " + ")
  formula_str <- sprintf("Surv(surv_time, event) ~ tx_int + %s", cov_formula)
  formula_obj <- as.formula(formula_str)

  t0 <- Sys.time()
  fit <- tryCatch({
    coxph(formula_obj, data = dt_sub, cluster = stay_id)
  }, error = function(e) {
    log_info(sprintf("[%s] CoxPH 回退失败: %s", label, e$message))
    NULL
  })
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (is.null(fit)) return(NULL)

  s <- summary(fit)
  hr <- exp(coef(fit)["tx_int"])
  ci <- tryCatch({
    exp(confint(fit)["tx_int", ])
  }, error = function(e) c(NA_real_, NA_real_))
  p <- tryCatch({
    s$coefficients["tx_int", "Pr(>|z|)"]
  }, error = function(e) NA_real_)

  out <- list(
    hr        = hr,
    ci_lower  = ci[1],
    ci_upper  = ci[2],
    p         = p,
    runtime   = runtime,
    n         = nrow(dt_sub),
    n_events  = sum(dt_sub$event, na.rm = TRUE),
    n_cov     = length(use_cov),
    method    = "coxph_fallback",
    covariates = use_cov
  )
  rm(fit, dt_sub); gc()
  return(out)
}

# ---- 3c. 运行分析并记录 ----
run_and_log <- function(dt, label, cov_set, boot = BOOT_SUB) {
  log_info(sprintf("开始分析: %s", label))
  log_info(sprintf("  n=%d, events=%d, cov_set=%d", nrow(dt), sum(dt$event), length(cov_set)))

  # 完整案例分析（删除含 NA 协变量的行）
  dt_cc <- dt[complete.cases(dt[, ..cov_set])]
  log_info(sprintf("  完整案例: %d 行, %d 事件", nrow(dt_cc), sum(dt_cc$event)))

  if (nrow(dt_cc) < 200 || sum(dt_cc$event) < 10) {
    log_info(sprintf("  [%s] 样本量或事件数不足 (n=%d, events=%d), 跳过", label, nrow(dt_cc), sum(dt_cc$event)))
    return(NULL)
  }

  # 主分析：SEQTaRget ITT
  r <- run_seqtarget_subgroup(dt_cc, label = label, boot = boot, cov_set = cov_set)

  # 回退：如果 SEQTaRget 失败，使用 CoxPH
  if (is.null(r)) {
    log_info(sprintf("  [%s] SEQTaRget 失败，回退到 CoxPH", label))
    r <- fallback_coxph(dt_cc, label = label, cov_set = cov_set)
    if (!is.null(r)) {
      log_info(sprintf("  [%s] CoxPH 回退: HR=%.4f (%.4f–%.4f), n=%d, events=%d, %.1f min",
                       label, r$hr, r$ci_lower, r$ci_upper, r$n, r$n_events, r$runtime / 60))
    } else {
      log_info(sprintf("  [%s] CoxPH 回退也失败", label))
    }
  } else {
    log_info(sprintf("  [%s] SEQTaRget ITT: HR=%.4f (%.4f–%.4f), n=%d, events=%d, %.1f min, cov=%d",
                     label, r$hr, r$ci_lower, r$ci_upper, r$n, r$n_events, r$runtime / 60, r$n_cov))
  }

  return(r)
}

log_info("辅助函数定义完成: run_seqtarget_subgroup, fallback_coxph, run_and_log")

# ============================================================================
# 4. 数据加载与准备
# ============================================================================
log_section("数据加载")

dt_full <- as.data.table(read_parquet(DATA_FILE))
dt_full[, `:=`(
  tx_int = as.integer(itt_group == "treatment"),
  event  = death_28day_num
)]
setorder(dt_full, stay_id, trial_num)

log_info(sprintf("Parquet 数据: %d 行, %d 事件, %d 列",
                 nrow(dt_full), sum(dt_full$event), ncol(dt_full)))

# 检查 sepsis3 变量是否存在
if (!"sepsis3" %in% names(dt_full)) {
  log_info("错误: 数据中不存在 sepsis3 变量，无法进行亚组分析")
  log_info("请检查数据源是否包含 sepsis3 变量")
  close(log_con)
  stop("缺少 sepsis3 变量，退出")
}

# 统计 sepsis3 分布
n_sepsis   <- dt_full[sepsis3 == 1, .N]
n_sepsis_ev <- dt_full[sepsis3 == 1, sum(event)]
n_nonsepsis   <- dt_full[sepsis3 == 0, .N]
n_nonsepsis_ev <- dt_full[sepsis3 == 0, sum(event)]

log_info(sprintf("Sepsis3 分布: %d 脓毒症 (%d 事件), %d 非脓毒症 (%d 事件)",
                 n_sepsis, n_sepsis_ev, n_nonsepsis, n_nonsepsis_ev))
log_info(sprintf("  sepsis3 缺失: %d", dt_full[is.na(sepsis3), .N]))

# ============================================================================
# 5. 分层亚组分析
# ============================================================================
log_section("分层亚组分析 (Sepsis3 分层)")

results_list <- list()

# ---- 5a. 非脓毒症亚组 (sepsis3 == 0, DAG13) ----
log_section("  5a. 非脓毒症亚组 (sepsis3 == 0, DAG13)")
dt_nonsepsis <- dt_full[sepsis3 == 0]
r_nonsepsis <- run_and_log(dt_nonsepsis,
  label = "NonSepsis_DAG13",
  cov_set = DAG13,
  boot = BOOT_SUB
)
results_list[["NonSepsis_DAG13"]] <- r_nonsepsis

# ---- 5b. 脓毒症亚组 (sepsis3 == 1, DAG13) ----
log_section("  5b. 脓毒症亚组 (sepsis3 == 1, DAG13)")
dt_sepsis <- dt_full[sepsis3 == 1]
r_sepsis <- run_and_log(dt_sepsis,
  label = "Sepsis_DAG13",
  cov_set = DAG13,
  boot = BOOT_SUB
)
results_list[["Sepsis_DAG13"]] <- r_sepsis

# ============================================================================
# 6. 结果汇总与对比表格
# ============================================================================
log_section("结果汇总与对比")

# ---- 6a. 构建对比表格 ----
comparison_rows <- list()

# 全队列 DAG14 结果（占位 NA，由主分析结果填充）
comparison_rows[["Full_DAG14"]] <- data.table(
  analysis   = "全队列-DAG14",
  hr         = NA_real_,
  ci_lower   = NA_real_,
  ci_upper   = NA_real_,
  p          = NA_real_,
  n          = NA_integer_,
  n_events   = NA_integer_,
  n_cov      = NA_integer_,
  method     = "SEQTaRget_ITT",
  cov_set    = "DAG14",
  runtime_min = NA_real_
)

# 非脓毒症 DAG13
if (!is.null(r_nonsepsis)) {
  comparison_rows[["NonSepsis_DAG13"]] <- data.table(
    analysis   = "非脓毒症-DAG13",
    hr         = r_nonsepsis$hr,
    ci_lower   = r_nonsepsis$ci_lower,
    ci_upper   = r_nonsepsis$ci_upper,
    p          = r_nonsepsis$p %||% NA_real_,
    n          = r_nonsepsis$n,
    n_events   = r_nonsepsis$n_events,
    n_cov      = r_nonsepsis$n_cov,
    method     = r_nonsepsis$method %||% "SEQTaRget_ITT",
    cov_set    = "DAG13",
    runtime_min = round(r_nonsepsis$runtime / 60, 2)
  )
} else {
  comparison_rows[["NonSepsis_DAG13"]] <- data.table(
    analysis = "非脓毒症-DAG13",
    hr = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
    p = NA_real_, n = NA_integer_, n_events = NA_integer_,
    n_cov = NA_integer_, method = "FAILED", cov_set = "DAG13",
    runtime_min = NA_real_
  )
}

# 脓毒症 DAG13
if (!is.null(r_sepsis)) {
  comparison_rows[["Sepsis_DAG13"]] <- data.table(
    analysis   = "脓毒症-DAG13",
    hr         = r_sepsis$hr,
    ci_lower   = r_sepsis$ci_lower,
    ci_upper   = r_sepsis$ci_upper,
    p          = r_sepsis$p %||% NA_real_,
    n          = r_sepsis$n,
    n_events   = r_sepsis$n_events,
    n_cov      = r_sepsis$n_cov,
    method     = r_sepsis$method %||% "SEQTaRget_ITT",
    cov_set    = "DAG13",
    runtime_min = round(r_sepsis$runtime / 60, 2)
  )
} else {
  comparison_rows[["Sepsis_DAG13"]] <- data.table(
    analysis = "脓毒症-DAG13",
    hr = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
    p = NA_real_, n = NA_integer_, n_events = NA_integer_,
    n_cov = NA_integer_, method = "FAILED", cov_set = "DAG13",
    runtime_min = NA_real_
  )
}

# 合并为 data.table
dt_comparison <- rbindlist(comparison_rows, use.names = TRUE, fill = TRUE)

# ---- 6b. 打印对比表格 ----
cat("\n")
cat("================================================================================\n")
cat("  eICU Sepsis3 缺失敏感性分析 — 结果对比\n")
cat("================================================================================\n\n")
cat(sprintf("%-22s %8s %14s %10s %10s %10s\n",
            "分析", "HR", "95% CI", "n", "events", "cov_set"))
cat(rep("-", 80), "\n", sep = "")
for (i in seq_len(nrow(dt_comparison))) {
  row <- dt_comparison[i]
  if (is.na(row$hr)) {
    cat(sprintf("%-22s %8s %14s %10s %10s %10s\n",
                row$analysis, "NA", "NA", "NA", "NA", row$cov_set))
  } else {
    ci_str <- sprintf("%.4f–%.4f", row$ci_lower, row$ci_upper)
    cat(sprintf("%-22s %8.4f %14s %10d %10d %10s\n",
                row$analysis, row$hr, ci_str, row$n, row$n_events, row$cov_set))
  }
}
cat(rep("-", 80), "\n", sep = "")
cat("\n")

# 写入日志
log_info("结果对比表格:")
for (i in seq_len(nrow(dt_comparison))) {
  row <- dt_comparison[i]
  if (is.na(row$hr)) {
    log_info(sprintf("  %s: NA", row$analysis))
  } else {
    log_info(sprintf("  %s: HR=%.4f (%.4f–%.4f), n=%d, events=%d, cov_set=%s",
                     row$analysis, row$hr, row$ci_lower, row$ci_upper,
                     row$n, row$n_events, row$cov_set))
  }
}

# ============================================================================
# 7. 输出保存
# ============================================================================
log_section("结果保存")

# ---- 7a. 对比 CSV ----
fwrite(dt_comparison, file.path(TBL_DIR, "v26_9_eicu_sepsis3_comparison.csv"))
log_info(sprintf("对比 CSV 已保存: %s",
                 file.path(TBL_DIR, "v26_9_eicu_sepsis3_comparison.csv")))

# ---- 7b. 完整结果 RDS ----
rds_out <- list(
  version           = "V26.9",
  timestamp         = as.character(Sys.time()),
  seed              = SEED,
  boot              = BOOT_SUB,
  dag14             = DAG14,
  dag13             = DAG13,
  n_total           = nrow(dt_full),
  n_total_events    = sum(dt_full$event),
  n_sepsis          = n_sepsis,
  n_sepsis_events   = n_sepsis_ev,
  n_nonsepsis       = n_nonsepsis,
  n_nonsepsis_events = n_nonsepsis_ev,
  results = list(
    nonsepsis_DAG13 = r_nonsepsis,
    sepsis_DAG13    = r_sepsis
  ),
  comparison        = dt_comparison
)
saveRDS(rds_out, file.path(RDS_DIR, "v26_9_eicu_sepsis3_sensitivity.rds"))
log_info(sprintf("完整 RDS 已保存: %s",
                 file.path(RDS_DIR, "v26_9_eicu_sepsis3_sensitivity.rds")))

# ---- 7c. 对比 JSON（供主分析脚本集成） ----
json_out <- list(
  version   = "V26.9",
  analysis  = "eICU Sepsis3 缺失敏感性分析",
  timestamp = as.character(Sys.time()),
  comparison = lapply(seq_len(nrow(dt_comparison)), function(i) {
    row <- dt_comparison[i]
    list(
      analysis   = row$analysis,
      hr         = if (is.na(row$hr)) NA else row$hr,
      ci_lower   = if (is.na(row$ci_lower)) NA else row$ci_lower,
      ci_upper   = if (is.na(row$ci_upper)) NA else row$ci_upper,
      p          = if (is.na(row$p)) NA else row$p,
      n          = if (is.na(row$n)) NA else row$n,
      n_events   = if (is.na(row$n_events)) NA else row$n_events,
      method     = row$method,
      cov_set    = row$cov_set
    )
  })
)
write_json(json_out, file.path(TBL_DIR, "v26_9_eicu_sepsis3_summary.json"),
           pretty = TRUE, auto_unbox = TRUE)
log_info(sprintf("对比 JSON 已保存: %s",
                 file.path(TBL_DIR, "v26_9_eicu_sepsis3_summary.json")))

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
all_pass <- all_pass && audit("DAG13 定义正确",
  length(DAG13) == 13 && !"sepsis3" %in% DAG13,
  sprintf("%d 变量, sepsis3 已排除", length(DAG13)))
all_pass <- all_pass && audit("非脓毒症亚组分析",
  !is.null(r_nonsepsis),
  if (!is.null(r_nonsepsis)) sprintf("HR=%.4f", r_nonsepsis$hr) else "失败")
all_pass <- all_pass && audit("脓毒症亚组分析",
  !is.null(r_sepsis),
  if (!is.null(r_sepsis)) sprintf("HR=%.4f", r_sepsis$hr) else "失败")
all_pass <- all_pass && audit("对比 CSV 输出",
  file.exists(file.path(TBL_DIR, "v26_9_eicu_sepsis3_comparison.csv")), "")
all_pass <- all_pass && audit("RDS 输出",
  file.exists(file.path(RDS_DIR, "v26_9_eicu_sepsis3_sensitivity.rds")), "")
all_pass <- all_pass && audit("JSON 输出",
  file.exists(file.path(TBL_DIR, "v26_9_eicu_sepsis3_summary.json")), "")
all_pass <- all_pass && audit("日志文件",
  file.exists(log_file), "")

log_info(sprintf("\n审核结果: %s", if (all_pass) "全部通过 ✓" else "存在未通过项 ✗"))

# ============================================================================
# 9. 收尾
# ============================================================================
close(log_con)

cat("\n")
cat("################################################################################\n")
cat("              V26.9 eICU Sepsis3 缺失敏感性分析完成\n")
cat("################################################################################\n\n")

cat(sprintf("DAG14 (%d): %s\n", length(DAG14), paste(DAG14, collapse = ", ")))
cat(sprintf("DAG13 (%d): %s\n", length(DAG13), paste(DAG13, collapse = ", ")))
cat(sprintf("非脓毒症 (n=%d, events=%d): ", n_nonsepsis, n_nonsepsis_ev))
if (!is.null(r_nonsepsis)) {
  cat(sprintf("HR=%.4f (%.4f–%.4f)\n", r_nonsepsis$hr, r_nonsepsis$ci_lower, r_nonsepsis$ci_upper))
} else {
  cat("失败\n")
}
cat(sprintf("脓毒症 (n=%d, events=%d): ", n_sepsis, n_sepsis_ev))
if (!is.null(r_sepsis)) {
  cat(sprintf("HR=%.4f (%.4f–%.4f)\n", r_sepsis$hr, r_sepsis$ci_lower, r_sepsis$ci_upper))
} else {
  cat("失败\n")
}
cat(sprintf("\n日志: %s\n", log_file))
cat(sprintf("对比 CSV: %s\n", file.path(TBL_DIR, "v26_9_eicu_sepsis3_comparison.csv")))
cat(sprintf("RDS: %s\n", file.path(RDS_DIR, "v26_9_eicu_sepsis3_sensitivity.rds")))
cat(sprintf("JSON: %s\n", file.path(TBL_DIR, "v26_9_eicu_sepsis3_summary.json")))
cat("========== V26.9 eICU Sepsis3 缺失敏感性分析完成 ==========\n")
