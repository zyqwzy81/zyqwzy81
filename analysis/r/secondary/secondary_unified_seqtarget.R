# ============================================================================
# V26.9 统一次要结局 SEQTaRget ITT 分析 — 合并 outcome_data + 批量运行
#
# 策略：
#   1. 将 outcome_data 所有可用变量合并到 imputed 数据
#   2. 对每项二分类次要结局运行 SEQTaRget ×5 + nboot=100 + DAG14
#   3. Rubin 合并 5 个插补集的 log(HR)
#   4. Bonferroni 校正（按结局数）
#
# 结局结构（根据团队决策 2026-07-03）：
#   主分析结局（7项，Bonferroni 校正）：
#     1. ICU死亡率 (death_icu)
#     2. 住院死亡率 (death_hospital)
#     3. 心血管复合结局 (cardiovascular_event_comprehensive)
#     4. AKI KDIGO≥2 (aki_kdigo2_plus)
#     5. 输血需求 (transfusion_at_t0)
#     6. 血磷达标 (phosphate_normalized)
#     7. 90天死亡率 (death_90day)
#   CV复合子组件（描述性，不纳入主校正）：
#     - 新发心律失常 (new_arrhythmia_rhythm)
#     - 心肌损伤 (myocardial_injury)
#     - 血管活性药升级 (vasopressor_escalation)
#   连续结局（替代方法，补充材料）：
#     - VFD-28
#     - ICU住院时长
#     这些连续结局在本脚本末尾单独处理。
#
# 参数：
#   N_IMPUTE = 5（前5个MICE插补集）
#   BOOT     = 100
#   协变量   = DAG14
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
BOOT      <- 25L      # nboot=25（nboot=50 卡死在 AKI）
SEED      <- 20260629L
set.seed(SEED)
options(mc.cores = 1)

# 注意：结局变量来源分段：
#   - 来源A（imputed 数据已有）：cardiovascular_event_comprehensive, aki_kdigo2_plus,
#     transfusion_at_t0, new_arrhythmia_rhythm, myocardial_injury, vasopressor_escalation
#   - 来源B（outcome_data）：death_icu, death_hospital, death_90day, vfd_28,
#     phosphate_normalized, icu_los_from_t0_days
#   - 来源C（augmented_with_90day 补全）：death_90day
# 次要结局定义
# - var: outcome 变量名
# - time: 生存时间变量名（SEQTaRget 需要 binary outcome，这里 time 用于过滤/日志）
# - label: 显示标签
# - group: "primary" = 主校正组, "sub" = CV子组件, "continuous" = 连续结局
# - primary_order: 主分析中的展示顺序
OUTCOME_DEFS <- list(
  # ---- 主分析结局（7项） ----
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
  # ---- CV 复合子组件（描述性，不纳入主校正） ----
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
log_file <- file.path(LOG_DIR, sprintf("v26_9_secondary_unified_seqtarget_%s.log",
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

log_info("V26.9 统一次要结局 SEQTaRget ITT 分析启动")
log_info(sprintf("种子: %d, nboot=%d, ncores=1, 插补集: 1-%d", SEED, BOOT, N_IMPUTE))

# ============================================================================
# 3. 数据加载 — 合并 outcome_data
# ============================================================================
log_section("数据加载")

# ---- 3a. 主数据 key 列 ----
dt_full <- as.data.table(read_parquet(DATA_FILE))
dt_full[, `:=`(
  tx_int = as.integer(itt_group == "treatment"),
  id_int = as.integer(factor(stay_id))
)]
setorder(dt_full, id_int, trial_num)

# 保留 imputed 数据中已有的结局变量 + 核心列
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

# ---- 3b. 合并 outcome_data（ICU/住院死亡、VFD、AKI、磷达标、LOS） ----
# 注意：两个数据集的 time_zero 编码不同，按 stay_id + 试验顺序号合并
dt_outcome <- as.data.table(read_parquet(OUTCOME_FILE))
outcome_cols <- c("stay_id", "time_zero",
                  "death_icu", "death_hospital", "death_28day", "death_90day",
                  "vfd_28", "total_mv_days",
                  "aki_stage2_plus",
                  "phosphate_normalized", "phosphate_normalization_days",
                  "icu_los_from_t0_days", "hospital_los_from_t0_days")
outcome_cols <- intersect(outcome_cols, names(dt_outcome))
dt_outcome <- dt_outcome[, ..outcome_cols]

# 为 outcome_data 添加试验顺序号（按 stay_id 内部 time_zero 排序）
setorder(dt_outcome, stay_id, time_zero)
dt_outcome[, trial_order := rowid(stay_id)]

# 为 dt_key 添加顺序号（trial_num 即为顺序号）
n_before <- nrow(dt_key)
dt_key[, trial_order := trial_num]

# 按 stay_id + trial_order 合并
dt_key <- merge(dt_key, dt_outcome[, setdiff(names(dt_outcome), c("time_zero", "death_28day")), with = FALSE],
                by = c("stay_id", "trial_order"), all.x = TRUE, sort = FALSE)
dt_outcome <- NULL; gc(TRUE)

log_info(sprintf("outcome_data 合并: %d→%d 行", n_before, nrow(dt_key)))
for (v in intersect(c("death_icu","death_hospital","death_90day",
                       "vfd_28","aki_stage2_plus",
                       "phosphate_normalized","icu_los_from_t0_days"),
                     names(dt_key))) {
  if (is.numeric(dt_key[[v]])) {
    log_info(sprintf("  %s: mean=%.3f, NA=%d", v,
                     mean(dt_key[[v]], na.rm = TRUE), sum(is.na(dt_key[[v]]))))
  }
}

# ---- 3c. 合并 death_90day（从 augmented_with_90day 补全） ----
if (file.exists(AUG_90DAY_FILE)) {
  dt_90 <- as.data.table(read_parquet(AUG_90DAY_FILE))
  setorder(dt_90, stay_id, time_zero)
  dt_90[, trial_order := rowid(stay_id)]
  dt_90 <- unique(dt_90[, .(stay_id, trial_order, death_90day)])

  # 如果 dt_key 的 death_90day 缺失，从 dt_90 补
  if (sum(is.na(dt_key$death_90day)) > 0) {
    dt_key[dt_90, on = c("stay_id", "trial_order"),
           death_90day := fifelse(is.na(death_90day), i.death_90day, death_90day)]
    log_info(sprintf("death_90day 补全: NA 剩余 %d",
                     sum(is.na(dt_key$death_90day))))
  }
  rm(dt_90); gc(TRUE)
}

# 清理辅助列
dt_key[, trial_order := NULL]

# ---- 3d. 创建 90天生存时间 ----
dt_key[, surv_time_90d := fifelse(
  death_28day_num == 1, surv_time, 90.0
)]
log_info("surv_time_90d 已创建")

# ---- 3e. 变量存在性检查 ----
log_section("结局变量存在性检查")
for (oname in names(OUTCOME_DEFS)) {
  def <- OUTCOME_DEFS[[oname]]
  exists_var <- def$var %in% names(dt_key)
  exists_time <- def$time %in% names(dt_key)

  if (exists_var) {
    n_events <- sum(dt_key[[def$var]] == 1, na.rm = TRUE)
    log_info(sprintf("  [%s] %s: %s %s events=%d",
                     oname, def$label,
                     if (exists_var) "var✓" else "var✗",
                     if (exists_time) "time✓" else "time✗",
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

  # 检查结局变量存在性
  if (!outcome_var %in% names(dt_sub)) {
    log_info(sprintf("  [%s] 结局变量不存在: %s", label, outcome_var))
    return(NULL)
  }

  # 检查事件数
  n_events <- sum(dt_sub[[outcome_var]] == 1, na.rm = TRUE)
  if (n_events < 10) {
    log_info(sprintf("  [%s] 事件数不足: %d < 10", label, n_events))
    return(NULL)
  }

  dt_sub[, `:=`(
    id_int       = as.integer(factor(stay_id)),
    period_int   = as.integer(trial_num) - 1L,
    a_int        = as.integer(itt_group == "treatment"),
    y_int        = as.integer(get(outcome_var)),
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

  # 协变量准备
  covs <- intersect(dag14_covs, names(dt_sub))
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
        trial.include    = TRUE,
        followup.include = TRUE,
        data.return      = FALSE
      ),
      verbose = FALSE
    )
  }, error = function(e) {
    log_info(sprintf("  [%s] SEQTaRget 错误: %s", label, e$message))
    NULL
  })
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
       ci_upper = exp(bm + 1.96 * se_pool), W = W, B = B)
}

log_info("辅助函数定义完成: run_seqtarget_generic, rubin_combine")

# ============================================================================
# 5. 遍历所有结局
# ============================================================================
log_section("多重结局 SEQTaRget 遍历")

all_results <- list()

for (oname in names(OUTCOME_DEFS)) {
  def   <- OUTCOME_DEFS[[oname]]
  ovar  <- def$var
  tvar  <- def$time
  olbl  <- def$label
  ogrp  <- def$group

  log_info(sprintf("\n======== 结局: %s (%s) ========", olbl, oname))

  # 跳过缺失结局变量
  if (!ovar %in% names(dt_key)) {
    log_info(sprintf("  跳过: %s 变量不存在", ovar))
    next
  }
  if (!tvar %in% names(dt_key)) {
    log_info(sprintf("  跳过: %s 时间变量不存在", tvar))
    next
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

  t_out_start <- Sys.time()

  for (i_imp in seq_len(N_IMPUTE)) {
    t_imp_start <- Sys.time()

    # ---- 5a. 加载 MICE Parquet ----
    parquet_file <- sprintf(MICE_PARQUET_TEMPL, i_imp)
    if (!file.exists(parquet_file)) {
      log_info(sprintf("  [imp%02d] 文件不存在", i_imp))
      imp_fail <- c(imp_fail, i_imp)
      next
    }
    dt_imp <- as.data.table(read_parquet(parquet_file))
    keep_cols <- intersect(c("stay_id", "trial_num", DAG14), names(dt_imp))
    dt_imp <- dt_imp[, ..keep_cols]
    gc(TRUE)

    # ---- 5b. Merge 对齐 ----
    dag14_in_imp <- intersect(DAG14, names(dt_imp))
    dag14_remove <- intersect(dag14_in_imp, names(dt_key))
    by_cols <- c("stay_id", "trial_num")

    dt_comb <- merge(
      dt_key[, setdiff(names(dt_key), dag14_remove), with = FALSE],
      dt_imp,
      by = by_cols, all.x = TRUE, sort = FALSE
    )

    # 字符型→factor
    for (v in intersect(DAG14, names(dt_comb))) {
      if (is.character(dt_comb[[v]])) dt_comb[[v]] <- as.factor(dt_comb[[v]])
    }

    # ---- 5c. 运行 SEQTaRget ----
    r <- tryCatch(
      run_seqtarget_generic(dt_comb, DAG14, ovar, tvar, BOOT, label = olbl),
      error = function(e) NULL
    )

    t_imp_end <- Sys.time()
    imp_elapsed <- as.numeric(difftime(t_imp_end, t_imp_start, units = "secs"))

    if (is.null(r)) {
      log_info(sprintf("  [imp%02d/%s] 失败, %.1f min", i_imp, oname, imp_elapsed / 60))
      imp_fail <- c(imp_fail, i_imp)
      rm(dt_comb, dt_imp); gc(TRUE)
      next
    }

    log_info(sprintf("  [imp%02d/%s] HR=%.4f (%.4f–%.4f), events=%d, %.1f min",
                     i_imp, oname, r$hr, r$ci_lower, r$ci_upper,
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

  t_out_end <- Sys.time()
  out_elapsed <- as.numeric(difftime(t_out_end, t_out_start, units = "secs"))

  m_actual <- length(imp_success)

  # ---- 5d. Rubin合并 ----
  if (m_actual >= 2) {
    rubin <- rubin_combine(imp_betas, imp_ses)
    hr_sd <- if (m_actual >= 2) sd(imp_hrs) else NA_real_

    log_info(sprintf("  [合并/%s] HR=%.4f (%.4f–%.4f), p=%.6f, m=%d, events=%.0f",
                     oname, rubin$hr, rubin$ci_lower, rubin$ci_upper,
                     rubin$p, m_actual, median(imp_nevt)))

    all_results[[oname]] <- list(
      def       = def,
      rubin     = rubin,
      hr_sd     = hr_sd,
      m_actual  = m_actual,
      imp_hrs   = imp_hrs,
      imp_betas = imp_betas,
      imp_ses   = imp_ses,
      imp_nevt  = imp_nevt,
      imp_success = imp_success,
      imp_fail  = imp_fail,
      runtime   = out_elapsed
    )
  } else {
    log_info(sprintf("  [跳过/%s] 不足2个成功插补集", oname))
  }
}

# ============================================================================
# 6. Bonferroni 校正（仅 primary 组）
# ============================================================================
log_section("Bonferroni 校正")

primary_results <- all_results[sapply(all_results, function(r) r$def$group == "primary")]
N_PRIMARY <- length(primary_results)
BONF_ALPHA <- 0.05 / N_PRIMARY

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
    runtime_min   = r$runtime / 60
  )
  dt_results_list[[oname]] <- dt_row
}

# CV 子组件（sub 组）
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
    runtime_min   = r$runtime / 60
  )
  dt_results_list[[oname]] <- dt_row
}

if (length(dt_results_list) > 0) {
  dt_all <- rbindlist(dt_results_list)
  setorder(dt_all, group, p_value)

  # 打印结果表
  cat("\n")
  cat(sprintf("%-45s %8s %12s %12s %12s %10s %8s\n",
              "结局", "HR", "CI下限", "CI上限", "p值", "Bonf.显著", "事件数"),
      file = log_con, append = TRUE)
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

  # 保存 CSV
  fwrite(dt_all, file.path(TBL_DIR, "v26_9_secondary_unified_results.csv"))
  log_info(sprintf("结果已保存: %s",
                   file.path(TBL_DIR, "v26_9_secondary_unified_results.csv")))
}

# ============================================================================
# 7. 连续结局处理（VFD-28, ICU LOS）— 补充材料
# ============================================================================
log_section("连续结局 — 手动IPTW加权线性回归")

# 手动 PS（Logistic）+ IPTW + 加权线性回归 + cluster-robust SE
# 不使用 SEQTaRget（SEQTaRget 不支持连续结局输出权重）

continuous_results <- list()
continuous_outcomes <- list(
  vfd    = list(var = "vfd_28", label = "28天无呼吸机天数 (VFD-28)"),
  icu_los = list(var = "icu_los_from_t0_days", label = "ICU住院天数 (ICU LOS)")
)

for (coname in names(continuous_outcomes)) {
  cdef <- continuous_outcomes[[coname]]
  cvar <- cdef$var
  clbl <- cdef$label

  if (!cvar %in% names(dt_key)) {
    log_info(sprintf("  [连续/%s] 变量不存在: %s", coname, cvar))
    next
  }

  # 用第 1 个插补集合并 DAG14 协变量
  parquet_file <- sprintf(MICE_PARQUET_TEMPL, 1L)
  if (!file.exists(parquet_file)) {
    log_info(sprintf("  [连续/%s] MICE 文件不存在", coname))
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

  # 准备
  dt_comb[, a_int := as.integer(itt_group == "treatment")]

  # 非 NA 结局
  dt_comb <- dt_comb[!is.na(get(cvar))]

  # 准备协变量
  covs <- intersect(DAG14, names(dt_comb))
  covs <- covs[sapply(covs, function(v) {
    x <- dt_comb[[v]]
    if (is.factor(x)) return(nlevels(x) > 1)
    if (is.character(x)) return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]
  if (length(covs) < 2) next
  for (v in covs) {
    if (is.character(dt_comb[[v]])) dt_comb[[v]] <- as.factor(dt_comb[[v]])
  }

  # 手动 Logistic PS
  ps_formula <- as.formula(paste("a_int ~", paste(covs, collapse = " + ")))
  ps_fit <- tryCatch(glm(ps_formula, data = dt_comb, family = binomial(), maxit = 100),
                     error = function(e) NULL)
  if (is.null(ps_fit)) next

  # 稳定化 IPTW + 双侧截断
  dt_comb[, ps := predict(ps_fit, type = "response")]
  stabilizer <- mean(dt_comb$a_int, na.rm = TRUE)
  dt_comb[, iptw_stable := fifelse(a_int == 1, stabilizer / ps, (1 - stabilizer) / (1 - ps))]
  tr_l <- quantile(dt_comb$iptw_stable, 0.01, na.rm = TRUE)
  tr_u <- quantile(dt_comb$iptw_stable, 0.99, na.rm = TRUE)
  dt_comb[, iptw_trunc := pmin(pmax(iptw_stable, tr_l), tr_u)]

  dt_lm <- dt_comb[!is.na(iptw_trunc) & is.finite(iptw_trunc)]
  if (nrow(dt_lm) < 100) next

  log_info(sprintf("  [连续/%s] n=%d, IPTW range: %.2f–%.2f",
                   coname, nrow(dt_lm), min(dt_lm$iptw_trunc), max(dt_lm$iptw_trunc)))

  # 加权线性回归（含基础协变量校正）
  formula_lm <- as.formula(sprintf(
    "%s ~ a_int + age + gender + baseline_sofa + sepsis3 + charlson + egfr + invasive_vent_at_t0 + vasoactive_at_t0 + baseline_phosphate",
    cvar
  ))
  fit_lm <- tryCatch(lm(formula_lm, data = dt_lm, weights = iptw_trunc),
                     error = function(e) NULL)
  if (is.null(fit_lm)) next

  # Cluster-robust SE
  cluster_vcov <- tryCatch(vcovCL(fit_lm, cluster = dt_lm$stay_id),
                           error = function(e) NULL)
  if (is.null(cluster_vcov)) next

  coef_mat <- coeftest(fit_lm, vcov. = cluster_vcov)
  tx_idx <- grep("a_int", rownames(coef_mat))
  if (length(tx_idx) != 1) next

  cn <- colnames(coef_mat)
  beta  <- coef_mat[tx_idx, "Estimate"]
  se    <- coef_mat[tx_idx, "Std. Error"]
  p_col <- if ("Pr(>|z|)" %in% cn) "Pr(>|z|)" else "Pr(>|t|)"
  p_val <- coef_mat[tx_idx, p_col]
  ci_l  <- beta - 1.96 * se
  ci_u  <- beta + 1.96 * se

  log_info(sprintf("  [连续/%s] Beta=%.4f (%.4f–%.4f), p=%.6f, n=%d",
                   coname, beta, ci_l, ci_u, p_val, nrow(dt_lm)))

  continuous_results[[coname]] <- data.table(
    outcome     = clbl,
    outcome_var = cvar,
    group       = "continuous",
    beta        = beta,
    ci_lower    = ci_l,
    ci_upper    = ci_u,
    p_value     = p_val,
    n           = nrow(dt_lm)
  )
  rm(dt_comb, dt_lm, ps_fit, fit_lm); gc(TRUE)
}

if (length(continuous_results) > 0) {
  dt_cont <- rbindlist(continuous_results)
  fwrite(dt_cont, file.path(TBL_DIR, "v26_9_continuous_outcomes.csv"))
  log_info(sprintf("连续结局结果已保存: %s",
                   file.path(TBL_DIR, "v26_9_continuous_outcomes.csv")))
}

# ============================================================================
# 8. 保存完整结果 RDS
# ============================================================================
log_section("保存 RDS")

rds_out <- list(
  version            = "V26.9",
  script             = "v26_9_secondary_unified_seqtarget.R",
  timestamp          = as.character(Sys.time()),
  seed               = SEED,
  dag14              = DAG14,
  n_impute           = N_IMPUTE,
  bootstrap          = BOOT,
  n_primary          = N_PRIMARY,
  bonferroni_alpha   = BONF_ALPHA,
  primary_results    = dt_all[group == "primary"],
  sub_results        = dt_all[group == "sub"],
  continuous_results = if (exists("dt_cont")) dt_cont else data.table(),
  per_outcome_detail = all_results
)
saveRDS(rds_out, file.path(RDS_DIR, "v26_9_secondary_unified_seqtarget.rds"))
log_info(sprintf("RDS 已保存: %s",
                 file.path(RDS_DIR, "v26_9_secondary_unified_seqtarget.rds")))

# ============================================================================
# 9. 流程自审核
# ============================================================================
log_section("流程自审核")

audit_checks <- list(
  c("数据加载", nrow(dt_key) > 10000, sprintf("n=%d", nrow(dt_key))),
  c("outcome_data合并", "death_icu" %in% names(dt_key),
    sprintf("death_icu存在=%s", "death_icu" %in% names(dt_key))),
  c("VFD可用", "vfd_28" %in% names(dt_key),
    sprintf("vfd_28存在=%s", "vfd_28" %in% names(dt_key))),
  c("主分析结局数", length(primary_results) >= 5,
    sprintf("n=%d", length(primary_results))),
  c("CV子组件", length(sub_results) >= 2,
    sprintf("n=%d", length(sub_results)))
)

for (ac in audit_checks) {
  pass <- ac[2] == TRUE || ac[2] == "TRUE"
  log_info(sprintf("  [%s] %s: %s",
                   if (pass) "PASS" else "FAIL", ac[1], ac[3]))
}

log_info(sprintf("\n日志文件: %s", log_file))
log_info(sprintf("结果CSV: %s",
                 file.path(TBL_DIR, "v26_9_secondary_unified_results.csv")))

close(log_con)

cat("\n")
cat("============================================================================\n")
cat("  V26.9 统一次要结局 SEQTaRget ITT 分析完成\n")
cat("============================================================================\n")
cat(sprintf("  主分析结局: %d 项 (Bonferroni α=%.4f)\n", N_PRIMARY, BONF_ALPHA))
cat(sprintf("  CV子组件:   %d 项 (描述性)\n", length(sub_results)))
cat(sprintf("  连续结局:   %d 项 (IPTW加权LM)\n", length(continuous_results)))
cat(sprintf("  日志: %s\n", log_file))
cat(sprintf("  结果: %s\n", file.path(TBL_DIR, "v26_9_secondary_unified_results.csv")))
cat("============================================================================\n")
