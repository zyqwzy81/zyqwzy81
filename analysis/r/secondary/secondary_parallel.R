# ============================================================================
# V26.9 统一次要结局 — 并行版本（mclapply 3核 × 9结局）
#
# 与 v26_9_secondary_unified_seqtarget.R 结构相同，
# 仅将结局遍历改为 mclapply 并行（3核，每核处理3结局~6h）。
#
# 检查点：每结局保存 checkpoint_rds 文件，避免重跑
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

OUT_DIR   <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9"
RDS_DIR   <- file.path(OUT_DIR, "rds")
TBL_DIR   <- file.path(OUT_DIR, "tables")
LOG_DIR   <- file.path(OUT_DIR, "logs")
for (d in c(OUT_DIR, RDS_DIR, TBL_DIR, LOG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

DATA_FILE           <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_imputed_v26.parquet"
OUTCOME_FILE        <- "/data/小米云盘/文档/论文/血磷/sql/result/M0_07_outcome_extraction.parquet"
AUG_90DAY_FILE      <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_augmented_with_90day.parquet"
MICE_PARQUET_TEMPL  <- "/data/小米云盘/文档/论文/血磷/sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet"

# DAG14 协变量集
DAG14 <- c(
  "age", "gender", "baseline_sofa", "sepsis3", "charlson",
  "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
  "admission_type", "baseline_phosphate", "lactate_t0", "refeeding_risk",
  "corrected_calcium_mmol_min", "nutrition_support_type"
)

N_IMPUTE  <- 5L
BOOT      <- 25L
SEED      <- 20260629L
set.seed(SEED)
options(mc.cores = 1)

OUTCOME_DEFS <- list(
  death_icu = list(
    var = "death_icu", time = "surv_time",
    label = "ICU死亡率 (ICU mortality)",
    group = "primary", order = 1,
    source = "outcome_data"
  ),
  death_hosp = list(
    var = "death_hospital", time = "surv_time",
    label = "住院死亡率 (In-hospital mortality)",
    group = "primary", order = 2,
    source = "outcome_data"
  ),
  cv_comp = list(
    var = "cardiovascular_event_comprehensive", time = "surv_time",
    label = "心血管复合结局 (CV composite)",
    group = "primary", order = 3,
    source = "imputed"
  ),
  aki = list(
    var = "aki_kdigo2_plus", time = "surv_time",
    label = "AKI KDIGO ≥2期 (AKI KDIGO ≥2)",
    group = "primary", order = 4,
    source = "imputed"
  ),
  transfusion = list(
    var = "transfusion_at_t0", time = "surv_time",
    label = "输血需求 (Transfusion need)",
    group = "primary", order = 5,
    source = "imputed"
  ),
  phos_norm = list(
    var = "phosphate_normalized", time = "surv_time",
    label = "血磷达标 (Phosphate normalization)",
    group = "primary", order = 6,
    source = "outcome_data"
  ),
  death_90d = list(
    var = "death_90day", time = "surv_time_90d",
    label = "90天死亡率 (90-day mortality)",
    group = "primary", order = 7,
    source = "outcome_data"
  ),
  arrhythmia = list(
    var = "new_arrhythmia_rhythm", time = "surv_time",
    label = "新发心律失常（CV子组件）",
    group = "sub", order = 1,
    source = "imputed"
  ),
  myocardial = list(
    var = "myocardial_injury", time = "surv_time",
    label = "心肌损伤（CV子组件）",
    group = "sub", order = 2,
    source = "imputed"
  ),
  vasopressor = list(
    var = "vasopressor_escalation", time = "surv_time",
    label = "血管活性药升级（CV子组件）",
    group = "sub", order = 3,
    source = "imputed"
  )
)

# ============================================================================
# 2. 日志系统
# ============================================================================
log_file <- file.path(LOG_DIR, sprintf("v26_9_secondary_parallel_%s.log",
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

log_info("V26.9 次要结局并行分析启动 (mclapply 3核)")
log_info(sprintf("种子: %d, nboot=%d, ncores=1/worker, 插补集: 1-%d, 结局数: %d",
                 SEED, BOOT, N_IMPUTE, length(OUTCOME_DEFS)))

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

IMP_OUTCOME_VARS <- c("cardiovascular_event_comprehensive", "aki_kdigo2_plus",
                       "transfusion_at_t0", "new_arrhythmia_rhythm",
                       "myocardial_injury", "vasopressor_escalation", "ecg_abnormal")
key_cols <- c("stay_id", "subject_id", "trial_num", "time_zero",
              "itt_group", "gender", "admission_type",
              "surv_time", "death_time_from_t0", "dischtime_from_t0",
              "death_28day_num",
              "nutrition_support_type",
              IMP_OUTCOME_VARS)
key_cols <- intersect(key_cols, names(dt_full))
dt_key   <- dt_full[, ..key_cols]
log_info(sprintf("主数据: %d 行, key cols=%d", nrow(dt_key), length(key_cols)))
rm(dt_full); gc(TRUE)

# 合并 outcome_data
dt_outcome <- as.data.table(read_parquet(OUTCOME_FILE))
outcome_cols <- c("stay_id", "time_zero",
                  "death_icu", "death_hospital", "death_28day", "death_90day",
                  "vfd_28", "total_mv_days",
                  "aki_stage2_plus",
                  "phosphate_normalized", "phosphate_normalization_days",
                  "icu_los_from_t0_days", "hospital_los_from_t0_days")
outcome_cols <- intersect(outcome_cols, names(dt_outcome))
dt_outcome <- dt_outcome[, ..outcome_cols]
setorder(dt_outcome, stay_id, time_zero)
dt_outcome[, trial_order := rowid(stay_id)]

n_before <- nrow(dt_key)
dt_key[, trial_order := trial_num]
dt_key <- merge(dt_key, dt_outcome[, setdiff(names(dt_outcome), c("time_zero", "death_28day")), with = FALSE],
                by = c("stay_id", "trial_order"), all.x = TRUE, sort = FALSE)
dt_outcome <- NULL; gc(TRUE)
log_info(sprintf("outcome_data 合并: %d→%d 行", n_before, nrow(dt_key)))

# 合并 death_90day 补全
if (file.exists(AUG_90DAY_FILE)) {
  dt_90 <- as.data.table(read_parquet(AUG_90DAY_FILE))
  setorder(dt_90, stay_id, time_zero)
  dt_90[, trial_order := rowid(stay_id)]
  dt_90 <- unique(dt_90[, .(stay_id, trial_order, death_90day)])
  if (sum(is.na(dt_key$death_90day)) > 0) {
    dt_key[dt_90, on = c("stay_id", "trial_order"),
           death_90day := fifelse(is.na(death_90day), i.death_90day, death_90day)]
  }
  rm(dt_90); gc(TRUE)
}
dt_key[, trial_order := NULL]

# 90天生存时间
dt_key[, surv_time_90d := fifelse(death_28day_num == 1, surv_time, 90.0)]

# 变量存在性检查
log_section("结局变量存在性检查")
for (oname in names(OUTCOME_DEFS)) {
  def <- OUTCOME_DEFS[[oname]]
  exists_var <- def$var %in% names(dt_key)
  exists_time <- def$time %in% names(dt_key)
  if (exists_var) {
    n_events <- sum(dt_key[[def$var]] == 1, na.rm = TRUE)
    log_info(sprintf("  [%s] %s var=%s time=%s events=%d",
                     oname, def$label,
                     if (exists_var) "✓" else "✗",
                     if (exists_time) "✓" else "✗",
                     n_events))
  } else {
    log_info(sprintf("  [%s] %s: 变量不存在!", oname, def$label))
  }
}

# ============================================================================
# 4. 辅助函数
# ============================================================================
log_section("辅助函数定义")

run_seqtarget_generic <- function(dt_comb, dag14_covs, outcome_var, time_var,
                                  boot, label = "") {
  dt_sub <- copy(dt_comb)
  if (!outcome_var %in% names(dt_sub)) return(NULL)
  n_events <- sum(dt_sub[[outcome_var]] == 1, na.rm = TRUE)
  if (n_events < 10) return(NULL)
  
  dt_sub[, `:=`(
    id_int       = as.integer(factor(stay_id)),
    period_int   = as.integer(trial_num) - 1L,
    a_int        = as.integer(itt_group == "treatment"),
    y_int        = as.integer(get(outcome_var)),
    eligible_int = 1L
  )]
  setorder(dt_sub, id_int, period_int)
  
  dt_sub[, `:=`(n_a = sum(a_int), n_y = sum(y_int)), by = period_int]
  dt_sub <- dt_sub[n_a > 0 & n_a < .N & n_y >= 5]
  if (nrow(dt_sub) < 100) return(NULL)
  
  covs <- intersect(dag14_covs, names(dt_sub))
  covs <- covs[sapply(covs, function(v) {
    x <- dt_sub[[v]]
    if (is.factor(x))    return(nlevels(x) > 1)
    if (is.character(x)) return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]
  if (length(covs) < 2) return(NULL)
  for (v in covs) {
    if (is.character(dt_sub[[v]])) dt_sub[, (v) := as.factor(get(v))]
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
        bootstrap       = boot > 0,
        bootstrap.nboot = boot,
        bootstrap.CI    = 0.95,
        ncores          = 1L,
        fastglm.method  = 3L,
        hazard          = TRUE,
        trial.include   = TRUE,
        followup.include = TRUE,
        data.return      = FALSE
      ),
      verbose = FALSE
    )
  }, error = function(e) NULL)
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  if (is.null(result)) return(NULL)
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

rubin_combine <- function(betas, ses) {
  m   <- length(betas)
  W   <- mean(ses^2)
  B   <- var(betas)
  VT  <- W + (1 + 1/m) * B
  se  <- sqrt(VT)
  beta_pooled <- mean(betas)
  t_stat <- beta_pooled / se
  df    <- (m - 1) * (1 + W / ((1 + 1/m) * B))^2
  p_val <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)
  list(
    beta = beta_pooled, se = se, t = t_stat, df = df, p = p_val,
    hr = exp(beta_pooled),
    ci_lower = exp(beta_pooled - 1.96 * se),
    ci_upper = exp(beta_pooled + 1.96 * se),
    W = W, B = B, VT = VT
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

log_info("辅助函数定义完成")

# ============================================================================
# 5. mclapply 并行处理所有结局
# ============================================================================
log_section(sprintf("并行处理 %d 个结局 (mclapply, 3 cores)", length(OUTCOME_DEFS)))

process_one_outcome <- function(oname) {
  def   <- OUTCOME_DEFS[[oname]]
  ovar  <- def$var
  tvar  <- def$time
  olbl  <- def$label
  ogrp  <- def$group

  # 检查检查点
  cp_file <- file.path(RDS_DIR, sprintf("secondary_checkpoint_%s.rds", oname))
  if (file.exists(cp_file)) {
    return(readRDS(cp_file))
  }

  if (!ovar %in% names(dt_key) || !tvar %in% names(dt_key)) {
    return(list(oname = oname, def = def, status = "skipped",
                reason = "variable_missing"))
  }

  imp_betas   <- numeric()
  imp_ses     <- numeric()
  imp_hrs     <- numeric()
  imp_ci_l    <- numeric()
  imp_ci_u    <- numeric()
  imp_run     <- numeric()
  imp_nevt    <- numeric()
  imp_success <- integer()
  imp_fail    <- integer()

  for (i_imp in seq_len(N_IMPUTE)) {
    parquet_file <- sprintf(MICE_PARQUET_TEMPL, i_imp)
    if (!file.exists(parquet_file)) {
      imp_fail <- c(imp_fail, i_imp)
      next
    }
    dt_imp <- as.data.table(read_parquet(parquet_file))
    keep_cols <- intersect(c("stay_id", "trial_num", DAG14), names(dt_imp))
    dt_imp <- dt_imp[, ..keep_cols]

    dag14_in_imp <- intersect(DAG14, names(dt_imp))
    dag14_remove <- intersect(dag14_in_imp, names(dt_key))
    dt_comb <- merge(
      dt_key[, setdiff(names(dt_key), dag14_remove), with = FALSE],
      dt_imp, by = c("stay_id", "trial_num"), all.x = TRUE, sort = FALSE
    )

    r <- tryCatch(
      run_seqtarget_generic(dt_comb, DAG14, ovar, tvar, BOOT, label = olbl),
      error = function(e) NULL
    )

    if (is.null(r)) {
      imp_fail <- c(imp_fail, i_imp)
      rm(dt_comb, dt_imp); gc(TRUE)
      next
    }

    imp_betas   <- c(imp_betas,   log(r$hr))
    imp_ses     <- c(imp_ses,     (log(r$ci_upper) - log(r$ci_lower)) / (2 * 1.96))
    imp_hrs     <- c(imp_hrs,     r$hr)
    imp_ci_l    <- c(imp_ci_l,    r$ci_lower)
    imp_ci_u    <- c(imp_ci_u,    r$ci_upper)
    imp_run     <- c(imp_run,     r$runtime)
    imp_nevt    <- c(imp_nevt,    r$n_events)
    imp_success <- c(imp_success, i_imp)
    rm(dt_comb, dt_imp); gc(TRUE)
  }

  m_actual <- length(imp_success)
  result <- list(
    oname = oname, def = def, status = "complete",
    m_actual = m_actual, imp_success = imp_success, imp_fail = imp_fail,
    imp_hrs = imp_hrs, imp_betas = imp_betas, imp_ses = imp_ses,
    imp_nevt = imp_nevt, runtime = sum(imp_run) / 60
  )

  if (m_actual >= 2) {
    rubin <- rubin_combine(imp_betas, imp_ses)
    result$rubin <- rubin
    result$hr_sd <- sd(imp_hrs)
  } else {
    result$rubin <- NULL
    result$hr_sd <- NA_real_
  }

  # 保存检查点
  saveRDS(result, cp_file)
  result
}

# 检查已有检查点
existing_cps <- list.files(RDS_DIR, pattern = "secondary_checkpoint_.*\\.rds")
completed <- gsub("secondary_checkpoint_|\\.rds", "", existing_cps)
all_names <- names(OUTCOME_DEFS)
remaining <- setdiff(all_names, completed)

if (length(existing_cps) >= 2) {
  log_info(sprintf("已有 %d/%d 结局检查点，跳过已完成项", length(existing_cps), length(OUTCOME_DEFS)))
  log_info(sprintf("已完成: %s", paste(completed, collapse = ", ")))
  log_info(sprintf("剩余: %s", paste(remaining, collapse = ", ")))
}

if (length(remaining) > 0) {
  # 只处理剩余结局
  names_subset <- remaining
} else {
  # 全部已完成
  names_subset <- character(0)
  log_info("所有结局已有检查点，跳过计算")
}

log_info(sprintf("启动并行计算: %d 个结局, 3 核", length(names_subset)))

if (length(names_subset) >= 2) {
  new_results <- parallel::mclapply(
    names_subset,
    process_one_outcome,
    mc.cores = min(2L, length(names_subset)),
    mc.preschedule = FALSE
  )
  log_info(sprintf("并行计算完成，处理了 %d 个结局", length(new_results)))
} else if (length(names_subset) == 1) {
  new_results <- list(process_one_outcome(names_subset[1]))
  log_info("单结局直接运行完成")
} else {
  new_results <- list()
}

# ============================================================================
# 6. 收集所有结果（从检查点恢复）
# ============================================================================
log_section("收集结果（从检查点恢复）")

all_results <- list()
for (oname in all_names) {
  cp_file <- file.path(RDS_DIR, sprintf("secondary_checkpoint_%s.rds", oname))
  if (file.exists(cp_file)) {
    r <- readRDS(cp_file)
    if (r$status == "complete" && !is.null(r$rubin)) {
      all_results[[oname]] <- r
      log_info(sprintf("  [%s] HR=%.4f (%.4f–%.4f), p=%.6f, m=%d",
                       oname, r$rubin$hr, r$rubin$ci_lower, r$rubin$ci_upper,
                       r$rubin$p, r$m_actual))
    }
  }
}

log_info(sprintf("收集完成: %d/%d 结局有有效结果", length(all_results), length(OUTCOME_DEFS)))

# ============================================================================
# 7. Bonferroni 校正 + 保存
# ============================================================================
log_section("Bonferroni 校正")

primary_results <- all_results[sapply(all_results, function(r) r$def$group == "primary")]
N_PRIMARY <- length(primary_results)
BONF_ALPHA <- 0.05 / max(N_PRIMARY, 1)
log_info(sprintf("主分析结局数: %d, Bonferroni α=%.6f", N_PRIMARY, BONF_ALPHA))

dt_results_list <- list()

for (oname in names(primary_results)) {
  r <- primary_results[[oname]]
  rubin <- r$rubin
  def   <- r$def
  dt_row <- data.table(
    outcome       = def$label,
    outcome_var   = def$var,
    group         = "primary",
    hr            = rubin$hr,
    ci_lower      = rubin$ci_lower,
    ci_upper      = rubin$ci_upper,
    p_value       = rubin$p,
    n_impute      = r$m_actual,
    n_events      = median(r$imp_nevt, na.rm = TRUE),
    hr_sd         = r$hr_sd,
    bonferroni    = fifelse(rubin$p < BONF_ALPHA, "Yes", "No"),
    bonf_alpha    = BONF_ALPHA,
    runtime_min   = r$runtime
  )
  dt_results_list[[oname]] <- dt_row
}

# CV 子组件
sub_results <- all_results[sapply(all_results, function(r) r$def$group == "sub")]
for (oname in names(sub_results)) {
  r <- sub_results[[oname]]
  rubin <- r$rubin
  def   <- r$def
  dt_row <- data.table(
    outcome       = def$label,
    outcome_var   = def$var,
    group         = "sub",
    hr            = rubin$hr,
    ci_lower      = rubin$ci_lower,
    ci_upper      = rubin$ci_upper,
    p_value       = rubin$p,
    n_impute      = r$m_actual,
    n_events      = median(r$imp_nevt, na.rm = TRUE),
    hr_sd         = r$hr_sd,
    bonferroni    = "N/A (sub-component)",
    bonf_alpha    = NA_real_,
    runtime_min   = r$runtime
  )
  dt_results_list[[oname]] <- dt_row
}

if (length(dt_results_list) > 0) {
  dt_all <- rbindlist(dt_results_list)
  setorder(dt_all, group, p_value)

  cat(sprintf("%-45s %8s %12s %12s %12s %10s %8s\n",
              "结局", "HR", "CI下限", "CI上限", "p值", "Bonf.显著", "事件数",
              file = log_con, append = TRUE))
  cat(sprintf("%-45s %8s %12s %12s %12s %10s %8s\n",
              "结局", "HR", "CI下限", "CI上限", "p值", "Bonf.显著", "事件数"))

  for (i in seq_len(nrow(dt_all))) {
    r <- dt_all[i]
    sig_txt <- if (r$group == "sub") "  -" else sprintf("%4s", r$bonferroni)
    line <- sprintf("  %-43s %8.4f %12.4f %12.4f %12.6f %10s %8.0f",
                    r$outcome, r$hr, r$ci_lower, r$ci_upper,
                    r$p_value, sig_txt, r$n_events)
    cat(line, "\n", file = log_con, append = TRUE)
    cat(line, "\n")
  }

  fwrite(dt_all, file.path(TBL_DIR, "v26_9_secondary_parallel_results.csv"))
  log_info(sprintf("结果已保存: %s",
                   file.path(TBL_DIR, "v26_9_secondary_parallel_results.csv")))
}

# ============================================================================
# 8. 连续结局
# ============================================================================
log_section("连续结局 — 手动IPTW加权线性回归")

continuous_results <- list()
continuous_outcomes <- list(
  vfd    = list(var = "vfd_28", label = "28天无呼吸机天数 (VFD-28)"),
  icu_los = list(var = "icu_los_from_t0_days", label = "ICU住院时长 (ICU LOS)")
)

# 简化版：仅用 imp01 做连续结局
for (coname in names(continuous_outcomes)) {
  convar <- continuous_outcomes[[coname]]$var
  conlbl <- continuous_outcomes[[coname]]$label

  if (!convar %in% names(dt_key)) {
    log_info(sprintf("  [%s] 变量不存在，跳过", conlbl))
    next
  }

  parquet_file <- sprintf(MICE_PARQUET_TEMPL, 1L)
  if (!file.exists(parquet_file)) {
    log_info(sprintf("  [连续/%s] MICE 文件不存在", conlbl))
    next
  }
  dt_imp <- as.data.table(read_parquet(parquet_file))
  keep_cols <- intersect(c("stay_id", "trial_num", DAG14), names(dt_imp))
  dt_imp <- dt_imp[, ..keep_cols]
  dag14_remove <- intersect(intersect(DAG14, names(dt_imp)), names(dt_key))
  dt_comb <- merge(
    dt_key[, setdiff(names(dt_key), dag14_remove), with = FALSE],
    dt_imp, by = c("stay_id", "trial_num"), all.x = TRUE, sort = FALSE
  )
  rm(dt_imp); gc(TRUE)

  # PS 模型
  dt_comb[, tx := as.integer(itt_group == "treatment")]
  ps_formula <- as.formula(paste("tx ~", paste(DAG14, collapse = " + ")))
  ps_model <- tryCatch(glm(ps_formula, data = dt_comb, family = binomial), error = function(e) NULL)
  if (is.null(ps_model)) {
    log_info(sprintf("  [%s] PS 模型失败", conlbl))
    next
  }
  dt_comb[, ps := predict(ps_model, type = "response")]
  dt_comb[, iptw := fifelse(tx == 1, 1/ps, 1/(1 - ps))]
  dt_comb[, iptw_trunc := pmin(iptw, quantile(iptw, 0.99, na.rm = TRUE))]

  # 加权线性回归
  lm_fit <- tryCatch(
    lm(as.formula(paste(convar, "~ tx")), data = dt_comb, weights = iptw_trunc),
    error = function(e) NULL
  )

  if (!is.null(lm_fit)) {
    coef_tx <- coef(lm_fit)["tx"]
    ci <- tryCatch(confint(lm_fit)["tx", ], error = function(e) c(NA, NA))
    pval <- tryCatch(summary(lm_fit)$coefficients["tx", "Pr(>|t|)"], error = function(e) NA)

    continuous_results[[coname]] <- list(
      var = convar, label = conlbl,
      coef = coef_tx, ci_lower = ci[1], ci_upper = ci[2],
      p_value = pval
    )

    log_info(sprintf("  [%s] Coef=%.3f (%.3f–%.3f), p=%.6f",
                     conlbl, coef_tx, ci[1], ci[2], pval))
  }
}

if (length(continuous_results) > 0) {
  dt_cont <- rbindlist(lapply(continuous_results, function(r) {
    data.table(
      outcome     = r$label,
      outcome_var = r$var,
      group       = "continuous",
      hr          = r$coef,
      ci_lower    = r$ci_lower,
      ci_upper    = r$ci_upper,
      p_value     = r$p_value,
      n_impute    = 1L,
      n_events    = NA_integer_,
      hr_sd       = NA_real_,
      bonferroni  = "N/A (continuous)",
      bonf_alpha  = NA_real_,
      runtime_min = NA_real_
    )
  }))
  fwrite(dt_cont, file.path(TBL_DIR, "v26_9_secondary_parallel_continuous.csv"))
  log_info(sprintf("连续结局已保存: %s",
                   file.path(TBL_DIR, "v26_9_secondary_parallel_continuous.csv")))
}

# ============================================================================
# 9. 清理检查点
# ============================================================================
log_section("清理检查点")
old_cps <- list.files(RDS_DIR, pattern = "secondary_checkpoint_.*\\.rds", full.names = TRUE)
if (length(old_cps) > 0) {
  file.remove(old_cps)
  log_info(sprintf("已清理 %d 个检查点文件", length(old_cps)))
}

log_info("V26.9 次要结局并行分析完成")
