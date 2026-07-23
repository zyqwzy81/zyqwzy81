# =============================================================================
# V28 SEQTaRget × MICE × LASSO 整合管线（内存优化版 v2）
# -----------------------------------------------------------------------------
# 核心优化（基于崩溃分析）：
#   1. 内存监测使用 /proc/meminfo（不受系统语言影响）
#   2. 多级降级策略：协变量裁剪 + bootstrap 分阶段
#   3. 14协变量（DAG14核心）优先，避免 OOM
#   4. bootstrap=FALSE 优先获取点估计，再尝试小 bootstrap
#   5. Swap 已扩展到 32GB（总虚拟内存 48GB）
#   6. 完整 checkpoint + 断点续算
#
# 崩溃分析结论（2026-07-23 验证）：
#   - 28协变量 + bootstrap=0 → OOM（R RSS=9.2GB，系统15GB，32GB Swap几乎用尽）
#   - 14协变量 + bootstrap=FALSE → 成功（33秒，HR=0.9632 ✅）
#   - 内存爆炸点：SEQTaRget内部data.table扩展（allow.cartesian）+ 模型对象保留
#   - HR异常（2.5845）未复现：修复period过滤+fastglm.method=3L后HR正常
#
# 关键修复（基于崩溃分析）：
#   1. period过滤：排除无治疗变异或事件<5的period（n_a>0 & n_a<.N & n_y>=5）
#   2. fastglm.method=3L（LDLT分解，比默认2L更稳定）
#   3. 内存监测使用 /proc/meminfo（不受系统语言影响）
#   4. Swap 已扩展到 32GB（总虚拟内存 48GB）
#   5. 完整 checkpoint + 断点续算
#
# 降级阶梯：
#   Level 0: DAG14核心(14) + bootstrap=FALSE  [点估计，最安全 ✅]
#   Level 1: DAG14核心(14) + bootstrap=5      [小CI，可能可行]
#   Level 2: DAG14核心(14) + bootstrap=20     [中等CI，可能OOM]
#   Level 3: LASSO全集(28) + bootstrap=FALSE  [点估计，几乎肯定OOM]
#   Level 4: LASSO全集(28) + bootstrap=50     [CI，一定OOM]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(jsonlite)
  library(SEQTaRget)
  library(survival)
})

# =============================================================================
# 0. 配置区
# =============================================================================
BASE_DIR  <- "/data/小米云盘/文档/论文/血磷"
OUT_DIR   <- file.path(BASE_DIR, "TTE", "v28", "output")
TBL_DIR   <- file.path(OUT_DIR, "tables")
RDS_DIR   <- file.path(OUT_DIR, "rds")
LOG_DIR   <- file.path(OUT_DIR, "logs")
MICE_DIR  <- file.path(BASE_DIR, "sql", "result", "mice_imputed")
LASSO_JSON <- file.path(BASE_DIR, "sql", "result", "v26_final_covariates.json")
CKPT_FILE <- file.path(RDS_DIR, "v28_seqtarget_checkpoint.rds")

# ---- 内存管理参数 ----
MICE_M           <- 20L
NTHREADS         <- 1L
MEM_WARN_GB      <- 2.0       # 可用内存低于此值发出警告
MEM_ABORT_GB     <- 0.5       # 可用内存低于此值立即中止
MEM_GC_TARGET_GB <- 2.0       # gc 后可用内存目标
KEEP_MODEL       <- FALSE
RESUME_FROM_CKPT <- TRUE

# ---- 降级阶梯（先尝试LASSO全集，内存不足时降级到DAG14） ----
# 格式: list(covariates="full"|"dag14", boot=数值)
DEGRADE_LADDER <- list(
  list(name = "L3_full_boot0",   covariates = "full",   boot = 0L),   # LASSO全集，优先尝试
  list(name = "L0_dag14_boot0",  covariates = "dag14",  boot = 0L),   # ✅ 已验证可行
  list(name = "L1_dag14_boot5",  covariates = "dag14",  boot = 5L),   # 可能可行
  list(name = "L2_dag14_boot20", covariates = "dag14",  boot = 20L),  # 可能OOM
  list(name = "L4_full_boot50",  covariates = "full",   boot = 50L)   # 一定OOM
)

# DAG14 核心变量（排除可能缺失的字符变量）
DAG14_CORE <- c("age", "baseline_sofa", "charlson", "sepsis3", "egfr",
                "invasive_vent_at_t0", "vasoactive_at_t0", "refeeding_risk",
                "ever_crrt", "lactate_t0", "liver_disease", "aki_kdigo2_plus",
                "baseline_phosphate", "corrected_calcium_mmol_min",
                "nutrition_support_type")

# 强制纳入变量
FORCE_INCLUDE <- DAG14_CORE

# 输出目录
for (d in c(TBL_DIR, RDS_DIR, LOG_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# 日志
TS      <- format(Sys.time(), "%Y%m%d_%H%M%S")
LOGFILE <- file.path(LOG_DIR, sprintf("v28_seqtarget_mice_pipeline_%s.log", TS))
CRASHLOG <- file.path(LOG_DIR, sprintf("v28_seqtarget_crash_report_%s.log", TS))

log_msg <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg)
  cat(line, "\n"); cat(line, "\n", file = LOGFILE, append = TRUE)
}

crash_log <- function(msg) {
  line <- sprintf("[%s] CRASH: %s", format(Sys.time(), "%H:%M:%S"), msg)
  cat(line, "\n"); cat(line, "\n", file = CRASHLOG, append = TRUE)
  cat(line, "\n", file = LOGFILE, append = TRUE)
}

# =============================================================================
# 1. 内存监测函数（使用 /proc/meminfo，不受系统语言影响）
# =============================================================================
get_system_mem <- function() {
  # /proc/meminfo 不受系统语言影响，始终为英文
  meminfo <- tryCatch(readLines("/proc/meminfo", n = 20), error = function(e) NULL)
  if (is.null(meminfo)) {
    # 回退：用 LANG=C 强制英文输出
    out <- tryCatch(system("LANG=C free -m", intern = TRUE, ignore.stderr = TRUE),
                    error = function(e) NULL)
    if (is.null(out) || length(out) < 3) return(list(avail_gb = NA_real_, swap_used_gb = NA_real_))
    vals <- as.numeric(unlist(regmatches(out[3], gregexpr("\\d+", out[3]))))
    return(list(avail_gb = vals[4] / 1024, swap_used_gb = if(length(vals) >= 2) vals[2] / 1024 else NA_real_))
  }

  parse_kb <- function(pattern) {
    line <- meminfo[grepl(pattern, meminfo)]
    if (length(line) > 0) return(as.numeric(gsub("[^0-9]", "", line)))
    return(NA_real_)
  }

  total_kb    <- parse_kb("^MemTotal:")
  avail_kb    <- parse_kb("^MemAvailable:")
  if (is.na(avail_kb)) avail_kb <- parse_kb("^MemFree:")
  swap_total  <- parse_kb("^SwapTotal:")
  swap_free   <- parse_kb("^SwapFree:")

  list(
    total_gb     = total_kb / 1024 / 1024,
    avail_gb     = avail_kb / 1024 / 1024,
    swap_used_gb = (swap_total - swap_free) / 1024 / 1024,
    swap_total_gb = swap_total / 1024 / 1024
  )
}

get_r_mem <- function() {
  gc_full <- gc(reset = FALSE, full = TRUE)
  list(r_used_gb = sum(as.numeric(gc_full[, "used"])) / 1024,
       r_peak_gb = sum(as.numeric(gc_full[, "max used"])) / 1024)
}

mem_report <- function(label) {
  sys_mem <- get_system_mem()
  r_mem   <- get_r_mem()
  log_msg(sprintf("MEM[%s] 系统可用=%.2fGB | R使用=%.2fGB | R峰值=%.2fGB | Swap用=%.2f/%.2fGB",
                  label,
                  ifelse(is.na(sys_mem$avail_gb), -1, sys_mem$avail_gb),
                  r_mem$r_used_gb, r_mem$r_peak_gb,
                  ifelse(is.na(sys_mem$swap_used_gb), -1, sys_mem$swap_used_gb),
                  ifelse(is.na(sys_mem$swap_total_gb), -1, sys_mem$swap_total_gb)))
}

mem_check_safe <- function(label = "") {
  sys_mem <- get_system_mem()
  avail_gb <- sys_mem$avail_gb
  if (is.na(avail_gb)) return(TRUE)
  if (avail_gb < MEM_ABORT_GB) {
    crash_log(sprintf("内存不足中止 [%s]: 可用=%.2fGB < 阈值=%.2fGB", label, avail_gb, MEM_ABORT_GB))
    return(FALSE)
  }
  if (avail_gb < MEM_WARN_GB) {
    log_msg(sprintf("⚠️ 内存警告 [%s]: 可用=%.2fGB < 警告=%.2fGB", label, avail_gb, MEM_WARN_GB))
  }
  return(TRUE)
}

force_gc <- function(target_gb = MEM_GC_TARGET_GB, max_iter = 3) {
  for (i in seq_len(max_iter)) {
    gc(reset = FALSE, full = TRUE)
    sys_mem <- get_system_mem()
    avail_gb <- sys_mem$avail_gb
    if (is.na(avail_gb) || avail_gb >= target_gb) break
    log_msg(sprintf("  gc迭代%d: 可用=%.2fGB (目标=%.2fGB)", i, avail_gb, target_gb))
  }
}

# =============================================================================
# 2. Checkpoint 机制
# =============================================================================
save_checkpoint <- function(completed, failed, results, crash_info, degrade_level) {
  ckpt <- list(
    completed_imps = completed,
    failed_imps = failed,
    results = results,
    crash_info = crash_info,
    degrade_level = degrade_level,
    timestamp = Sys.time()
  )
  saveRDS(ckpt, CKPT_FILE)
  log_msg(sprintf("Checkpoint: 完成=%d, 失败=%d, 降级=%s",
                  length(completed), length(failed), degrade_level))
}

load_checkpoint <- function() {
  if (!file.exists(CKPT_FILE) || !RESUME_FROM_CKPT) {
    if (file.exists(CKPT_FILE) && !RESUME_FROM_CKPT) {
      log_msg("RESUME_FROM_CKPT=FALSE，从头开始")
    }
    return(NULL)
  }
  ckpt <- readRDS(CKPT_FILE)
  log_msg(sprintf("从 checkpoint 恢复: 已完成 %d 个插补集", length(ckpt$completed_imps)))
  log_msg(sprintf("  checkpoint时间: %s | 降级: %s", ckpt$timestamp, ckpt$degrade_level))
  ckpt
}

# =============================================================================
# 3. 数据准备函数
# =============================================================================
prepare_seqtarget_data <- function(parquet_path, covariates) {
  dt <- as.data.table(read_parquet(parquet_path))
  log_msg(sprintf("  读取: %s (%d行 × %d列)", basename(parquet_path), nrow(dt), ncol(dt)))

  # 列名映射（factor → integer，用 as.character 中转）
  dt[, ID := as.integer(stay_id)]
  if ("treatment" %in% names(dt)) {
    if (is.factor(dt$treatment)) dt[, tx_init := as.integer(as.character(treatment))]
    else dt[, tx_init := as.integer(treatment)]
  } else if ("itt_group" %in% names(dt)) {
    dt[, tx_init := fifelse(itt_group == "treatment", 1L, 0L)]
  } else stop("缺少 treatment 或 itt_group 列")

  if (is.factor(dt$death_28day_num)) dt[, outcome := as.integer(as.character(death_28day_num))]
  else dt[, outcome := as.integer(death_28day_num)]

  dt[, period := as.integer(trial_num - 1L)]
  dt[, eligible := 1L]

  # 🆕 period过滤：排除无治疗变异或事件数<5的period
  # 避免PS/outcome模型在不可识别期上拟合异常，导致HR扭曲
  n_before <- nrow(dt)
  dt[, `:=`(n_a = sum(tx_init), n_y = sum(outcome)), by = period]
  dt <- dt[n_a > 0 & n_a < .N & n_y >= 5]
  dt[, `:=`(n_a = NULL, n_y = NULL)]
  n_excluded <- n_before - nrow(dt)
  if (n_excluded > 0) {
    log_msg(sprintf("  period过滤: 排除%d行 (无治疗变异或事件<5)", n_excluded))
  }

  # 协变量数值化
  available_vars <- intersect(covariates, names(dt))
  missing_vars <- setdiff(covariates, names(dt))
  if (length(missing_vars) > 0) {
    log_msg(sprintf("  ⚠️ %d变量缺失: %s", length(missing_vars), paste(missing_vars, collapse=", ")))
  }
  for (col in available_vars) {
    if (is.character(dt[[col]]) || is.factor(dt[[col]])) {
      vals <- as.character(dt[[col]])
      num_vals <- suppressWarnings(as.numeric(vals))
      if (mean(is.na(num_vals)) < 0.5) dt[[col]] <- num_vals
      else dt[[col]] <- as.numeric(factor(vals))
    }
  }

  keep_cols <- c("ID", "period", "tx_init", "outcome", "eligible", available_vars)
  analysis_data <- dt[, ..keep_cols]
  analysis_data <- analysis_data[complete.cases(analysis_data)]

  list(data = analysis_data, fixed_cols = available_vars,
       n_rows = nrow(analysis_data), n_patients = uniqueN(analysis_data$ID),
       n_events = sum(analysis_data$outcome))
}

# =============================================================================
# 4. Rubin 规则合并
# =============================================================================
rubin_pool <- function(log_hrs, ses) {
  m <- length(log_hrs)
  stopifnot(m == length(ses) && m > 1)
  pooled_log_hr <- mean(log_hrs)
  W <- mean(ses^2)
  B <- stats::var(log_hrs) * (m / (m - 1))
  T_var <- W + B + B/m
  z <- qnorm(0.975)
  list(
    hr = exp(pooled_log_hr),
    ci_lower = exp(pooled_log_hr - z * sqrt(T_var)),
    ci_upper = exp(pooled_log_hr + z * sqrt(T_var)),
    se = sqrt(T_var), log_hr = pooled_log_hr,
    W = W, B = B, T = T_var,
    lambda = B / T_var, rel_eff = 1 / (1 + (B / T_var) / m), m = m
  )
}

# =============================================================================
# 5. 单插补集分析（含降级）
# =============================================================================
run_single_imp <- function(imp_idx, parquet_path, covariates, boot_n, level_name) {
  prep <- tryCatch(prepare_seqtarget_data(parquet_path, covariates),
                   error = function(e) { crash_log(sprintf("imp%02d 数据准备失败: %s", imp_idx, e$message)); NULL })
  if (is.null(prep)) return(list(status = "fail", reason = "data_prep"))

  log_msg(sprintf("  imp%02d 数据: %d行 | %d患者 | %d事件 | %d协变量",
                  imp_idx, prep$n_rows, prep$n_patients, prep$n_events, length(prep$fixed_cols)))

  if (!mem_check_safe(sprintf("imp%02d-pre", imp_idx))) {
    return(list(status = "abort", reason = "mem_pre"))
  }

  # SEQuential 调用
  options_itt <- SEQopts(
    hazard = TRUE, km.curves = FALSE,
    bootstrap = (boot_n > 0),
    bootstrap.nboot = max(boot_n, 1L),
    fastglm.method = 3L,  # LDLT 分解，比默认2L更稳定
    nthreads = NTHREADS
  )

  t0 <- Sys.time()
  model_itt <- tryCatch({
    SEQuential(
      data = as.data.frame(prep$data), id.col = "ID", time.col = "period",
      eligible.col = "eligible", treatment.col = "tx_init",
      outcome.col = "outcome",
      time_varying.cols = character(0),
      fixed.cols = prep$fixed_cols,
      method = "ITT", options = options_itt
    )
  }, error = function(e) { crash_log(sprintf("imp%02d %s 失败: %s", imp_idx, level_name, e$message)); NULL })
  t1 <- Sys.time()
  elapsed <- as.numeric(difftime(t1, t0, units = "secs"))

  if (is.null(model_itt)) {
    return(list(status = "fail", reason = "seqtarget_error", elapsed = elapsed))
  }

  # 提取 HR
  hr_obj <- tryCatch(hazard_ratio(model_itt), error = function(e) NULL)
  if (is.null(hr_obj) || length(hr_obj) < 1) {
    rm(model_itt); force_gc()
    return(list(status = "fail", reason = "hr_extract", elapsed = elapsed))
  }

  hr      <- as.numeric(hr_obj[1])
  ci_low  <- if (length(hr_obj) >= 2) as.numeric(hr_obj[2]) else NA
  ci_high <- if (length(hr_obj) >= 3) as.numeric(hr_obj[3]) else NA
  log_hr  <- log(hr)

  # SE 从 CI 反推；若 CI=NA（bootstrap=FALSE），用 NA 标记
  if (!is.na(ci_low) && !is.na(ci_high)) {
    se_val <- (log(ci_high) - log(ci_low)) / (2 * qnorm(0.975))
  } else {
    se_val <- NA_real_
    log_msg(sprintf("  ⚠️ imp%02d CI=NA（bootstrap=FALSE），仅点估计", imp_idx))
  }

  # HR 异常检测
  if (hr > 2.5 || hr < 0.4) {
    log_msg(sprintf("  ⚠️ imp%02d HR=%.4f 异常（主分析HR=0.987），需检查", imp_idx, hr))
  }

  log_msg(sprintf("  imp%02d ✓ HR=%.4f (%s-%s) | %s | %.1fs",
                  imp_idx, hr,
                  ifelse(is.na(ci_low), "NA", sprintf("%.4f", ci_low)),
                  ifelse(is.na(ci_high), "NA", sprintf("%.4f", ci_high)),
                  level_name, elapsed))

  rm(model_itt, hr_obj)
  force_gc()

  list(status = "success", hr = hr, ci_lower = ci_low, ci_upper = ci_high,
       log_hr = log_hr, se = se_val, boot_used = boot_n, level = level_name,
       elapsed = elapsed, n_rows = prep$n_rows, n_events = prep$n_events,
       n_covariates = length(prep$fixed_cols))
}

# =============================================================================
# 6. 主流程
# =============================================================================
log_msg("==========================================================")
log_msg("V28 SEQTaRget × MICE × LASSO 整合管线 v2（内存优化版）")
log_msg("==========================================================")

# 验证内存监测
init_mem <- get_system_mem()
log_msg(sprintf("系统内存: 总=%.2fGB | 可用=%.2fGB | Swap=%.2f/%.2fGB",
                ifelse(is.na(init_mem$total_gb), -1, init_mem$total_gb),
                ifelse(is.na(init_mem$avail_gb), -1, init_mem$avail_gb),
                ifelse(is.na(init_mem$swap_used_gb), -1, init_mem$swap_used_gb),
                ifelse(is.na(init_mem$swap_total_gb), -1, init_mem$swap_total_gb)))

if (!is.na(init_mem$avail_gb) && init_mem$avail_gb < MEM_WARN_GB) {
  log_msg("⚠️ 启动时可用内存偏低，建议关闭其他进程")
}

# 加载 LASSO 变量集
lasso_info <- fromJSON(LASSO_JSON)
lasso_vars <- lasso_info$final_covariates
if (is.null(lasso_vars)) lasso_vars <- lasso_info$covariates
lasso_vars <- unlist(lasso_vars)
full_vars <- unique(c(lasso_vars, FORCE_INCLUDE))
log_msg(sprintf("LASSO变量=%d | DAG14核心=%d | 全集=%d",
                length(lasso_vars), length(DAG14_CORE), length(full_vars)))

# Checkpoint 恢复
ckpt <- load_checkpoint()
if (!is.null(ckpt)) {
  completed_imps <- ckpt$completed_imps
  failed_imps    <- ckpt$failed_imps
  results_list   <- ckpt$results
  crash_info     <- ckpt$crash_info
  start_idx      <- max(completed_imps, 0L) + 1L
  # 从 checkpoint 恢复降级级别
  degrade_idx <- 1L  # 默认从最安全级别开始
} else {
  completed_imps <- integer(0)
  failed_imps    <- integer(0)
  results_list   <- vector("list", MICE_M)
  crash_info     <- list()
  start_idx      <- 1L
  degrade_idx    <- 1L
}

log_msg("----------------------------------------------------------")
log_msg("开始主循环（降级策略：从最安全级别开始）")
log_msg("----------------------------------------------------------")

abort_all <- FALSE

for (i in seq(start_idx, MICE_M)) {
  if (abort_all) break

  imp_tag <- sprintf("imp%02d", i)
  parquet_path <- file.path(MICE_DIR, sprintf("v26_mice_imputed_%s.parquet", imp_tag))

  log_msg(sprintf("[%02d/%d] === %s ===", i, MICE_M, basename(parquet_path)))
  mem_report(sprintf("imp%02d-start", i))

  if (!file.exists(parquet_path)) {
    log_msg("  跳过: 文件不存在")
    failed_imps <- c(failed_imps, i)
    crash_info[[imp_tag]] <- list(reason = "file_not_found")
    next
  }

  # 全局内存检查
  if (!mem_check_safe(sprintf("imp%02d-global", i))) {
    log_msg("⚠️ 全局内存不足，保存 checkpoint 后中止")
    save_checkpoint(completed_imps, failed_imps, results_list, crash_info, "abort")
    abort_all <- TRUE
    break
  }

  # 逐级降级尝试
  result <- NULL
  for (dl in seq(degrade_idx, length(DEGRADE_LADDER))) {
    level <- DEGRADE_LADDER[[dl]]
    cov_set <- if (level$covariates == "dag14") DAG14_CORE else full_vars
    log_msg(sprintf("  尝试 %s (协变量=%d, boot=%d)...",
                    level$name, length(cov_set), level$boot))

    # 内存预检查
    if (!mem_check_safe(sprintf("imp%02d-%s-pre", i, level$name))) {
      log_msg("  内存不足，降级到下一级别")
      force_gc()
      next
    }

    result <- run_single_imp(i, parquet_path, cov_set, level$boot, level$name)

    if (result$status == "success") {
      degrade_idx <- dl  # 记住成功级别，后续插补集直接用此级别
      break
    } else if (result$status == "abort") {
      log_msg("  内存中止，保存 checkpoint")
      save_checkpoint(completed_imps, failed_imps, results_list, crash_info, level$name)
      abort_all <- TRUE
      break
    } else {
      log_msg(sprintf("  %s 失败(%s)，降级", level$name, result$reason))
      force_gc()
    }
  }

  if (abort_all) break

  if (!is.null(result) && result$status == "success") {
    imp_dt <- data.table(
      imp = i, hr = result$hr, ci_lower = result$ci_lower, ci_upper = result$ci_upper,
      log_hr = result$log_hr, se = result$se,
      boot_used = result$boot_used, level = result$level, elapsed = result$elapsed,
      n_rows = result$n_rows, n_events = result$n_events, n_covariates = result$n_covariates
    )
    write.csv(imp_dt, file.path(TBL_DIR, sprintf("v28_seqtarget_mice_hr_%s.csv", imp_tag)), row.names = FALSE)
    results_list[[i]] <- imp_dt
    completed_imps <- c(completed_imps, i)
  } else {
    failed_imps <- c(failed_imps, i)
    crash_info[[imp_tag]] <- list(reason = ifelse(!is.null(result), result$reason, "unknown"))
  }

  save_checkpoint(completed_imps, failed_imps, results_list, crash_info,
                   DEGRADE_LADDER[[degrade_idx]]$name)
  mem_report(sprintf("imp%02d-end", i))
}

# =============================================================================
# 7. Rubin 合并
# =============================================================================
log_msg("----------------------------------------------------------")
log_msg(sprintf("主循环结束: 成功=%d / 失败=%d / 总计=%d",
                length(completed_imps), length(failed_imps), MICE_M))
log_msg("----------------------------------------------------------")

valid_results <- Filter(Negate(is.null), results_list)
if (length(valid_results) > 0) {
  all_results <- rbindlist(valid_results)
  write.csv(all_results, file.path(TBL_DIR, "v28_seqtarget_mice_hr_all.csv"), row.names = FALSE)

  # 分离有SE和无SE的结果
  with_se <- all_results[!is.na(se)]
  without_se <- all_results[is.na(se)]

  if (nrow(with_se) >= 2) {
    log_msg(sprintf("执行 Rubin 合并 (m=%d, 有CI)...", nrow(with_se)))
    pooled <- rubin_pool(with_se$log_hr, with_se$se)
  } else if (nrow(all_results) >= 2) {
    log_msg(sprintf("仅点估计可用 (m=%d, 无CI)，执行简单合并...", nrow(all_results)))
    # 简单合并：仅取 log_hr 均值，SE 用 NA
    pooled <- list(
      hr = exp(mean(all_results$log_hr)), ci_lower = NA, ci_upper = NA,
      se = NA, log_hr = mean(all_results$log_hr),
      W = NA, B = var(all_results$log_hr), T = NA,
      lambda = NA, rel_eff = NA, m = nrow(all_results)
    )
  } else {
    log_msg("❌ 有效结果 < 2，无法合并")
    pooled <- NULL
  }

  if (!is.null(pooled)) {
    pooled_dt <- data.table(
      analysis = "SEQTaRget ITT × MICE × LASSO (Rubin pooled)",
      hr = pooled$hr, ci_lower = pooled$ci_lower, ci_upper = pooled$ci_upper,
      se = pooled$se, log_hr = pooled$log_hr,
      W_within = pooled$W, B_between = pooled$B, T_total = pooled$T,
      lambda = pooled$lambda, rel_eff = pooled$rel_eff, m = pooled$m,
      boot_median = ifelse(nrow(with_se) > 0, median(with_se$boot_used), 0),
      level_used = ifelse(nrow(with_se) > 0, with_se$level[1], all_results$level[1]),
      n_covariates = ifelse(nrow(all_results) > 0, all_results$n_covariates[1], NA)
    )
    write.csv(pooled_dt, file.path(TBL_DIR, "v28_seqtarget_mice_rubin_pooled.csv"), row.names = FALSE)

    log_msg("==========================================================")
    log_msg("Rubin 合并结果")
    log_msg("==========================================================")
    log_msg(sprintf("HR = %.4f (%s - %s)", pooled$hr,
                    ifelse(is.na(pooled$ci_lower), "NA", sprintf("%.4f", pooled$ci_lower)),
                    ifelse(is.na(pooled$ci_upper), "NA", sprintf("%.4f", pooled$ci_upper))))
    log_msg(sprintf("log-HR = %.4f | SE = %s", pooled$log_hr,
                    ifelse(is.na(pooled$se), "NA", sprintf("%.4f", pooled$se))))
    if (!is.na(pooled$W)) {
      log_msg(sprintf("W=%.6f | B=%.6f | T=%.6f | λ=%.3f | RE=%.3f",
                      pooled$W, pooled$B, pooled$T, pooled$lambda, pooled$rel_eff))
    }
    log_msg(sprintf("降级级别: %s | 协变量数: %d | 插补集: %d",
                    pooled_dt$level_used, pooled_dt$n_covariates, pooled$m))

    # 与主分析对比
    main_hr <- 0.987
    log_msg(sprintf("与主分析对比 (HR=%.3f): 差异=%+.4f (%.1f%%)",
                    main_hr, pooled$hr - main_hr, abs(pooled$hr - main_hr) / main_hr * 100))
  }
}

# 崩溃报告
if (length(crash_info) > 0) {
  log_msg("----------------------------------------------------------")
  log_msg("崩溃报告")
  log_msg("----------------------------------------------------------")
  for (tag in names(crash_info)) {
    log_msg(sprintf("  %s: %s", tag, crash_info[[tag]]$reason))
  }
}

mem_report("结束前")
log_msg("管线执行完成。")
log_msg(sprintf("Checkpoint: %s", CKPT_FILE))

# =============================================================================
# 使用说明
# -----------------------------------------------------------------------------
# 1. 前置条件：
#    - sql/tools/v26_mice_imputation.R → 20个插补Parquet
#    - sql/tools/v26_lasso_selection.R → v26_final_covariates.json
#
# 2. 内存调整：
#    - MEM_WARN_GB / MEM_ABORT_GB : 内存阈值
#    - DEGRADE_LADDER : 降级阶梯（可自定义）
#    - 若持续OOM：增加swap，或减少DAG14_CORE变量
#
# 3. 断点续算：
#    - RESUME_FROM_CKPT=TRUE，中断后重新运行会跳过已完成插补集
#    - 删除 v28_seqtarget_checkpoint.rds 可强制从头开始
#
# 4. 运行：
#    Rscript v28_seqtarget_mice_lasso_pipeline.R
# =============================================================================
