# ============================================================================
# V26.9 宽限期敏感性分析 — SEQTaRget CCW-ITT (24h Grace Period)
#
# 目的：在 SEQTaRget 框架下评估 24h 宽限期对主分析结果的影响
#   主分析使用 12h 宽限期，此处扩展至 24h
#   
# 方法：使用 iv_within_24h 列重新定义治疗分配，SEQTaRget ITT × 5 imp × nboot=50
#   协变量集：DAG14
#   结局：28天死亡率
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

N_IMPUTE  <- 5L        # 前 5 个插补集
BOOT      <- 50L       # nboot=50
SEED      <- 20260705L
set.seed(SEED)
options(mc.cores = 1)

DAG14 <- c(
  "age", "gender", "baseline_sofa", "sepsis3", "charlson",
  "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
  "admission_type", "baseline_phosphate", "lactate_t0", "refeeding_risk",
  "corrected_calcium_mmol_min", "nutrition_support_type"
)

# ============================================================================
# 2. 日志系统
# ============================================================================
log_file <- file.path(LOG_DIR, sprintf("v26_9_sensitivity_grace_%s.log",
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

log_info("V26.9 宽限期敏感性分析 (24h Grace Period) — SEQTaRget CCW-ITT")
log_info(sprintf("种子=%d, nboot=%d, 插补集: 1-%d", SEED, BOOT, N_IMPUTE))
log_info(sprintf("DAG14: %d 变量", length(DAG14)))

# ============================================================================
# 3. 数据加载
# ============================================================================
log_section("数据加载")
dt_full <- as.data.table(read_parquet(DATA_FILE))
dt_full[, `:=`(
  tx_int = as.integer(itt_group == "treatment"),
  id_int = as.integer(factor(stay_id))
)]

# ---- 创建 24h 宽限期治疗分配 ----
# 当前 itt_group 基于 12h 宽限期。对于 24h 宽限期：
# 将 iv_within_24h==1 的患者重新定义为治疗组，
# 将 iv_within_24h==0 且曾在任何时间接受 IV 磷的患者重新定义为对照组
# 从未接受 IV 磷的患者保持对照组
dt_full[, tx_24h := as.integer(iv_within_24h == 1)]
log_info(sprintf("24h 治疗分配: 治疗=%d, 对照=%d (vs 12h: 治疗=%d, 对照=%d)",
                 sum(dt_full$tx_24h), sum(dt_full$tx_24h == 0),
                 sum(dt_full$tx_int), sum(dt_full$tx_int == 0)))

setorder(dt_full, id_int, trial_num)
key_cols <- c("stay_id", "subject_id", "trial_num", "time_zero",
              "itt_group", "tx_24h",
              "surv_time", "death_28day_num")
key_cols <- intersect(key_cols, names(dt_full))
dt_key   <- dt_full[, ..key_cols]
log_info(sprintf("主数据 key 表: %d 行, %d 列", nrow(dt_key), ncol(dt_key)))
rm(dt_full); gc(TRUE)

# ============================================================================
# 4. 辅助函数
# ============================================================================
log_section("辅助函数定义")

run_seqtarget_grace <- function(dt_key, dag14_covs, boot = BOOT) {
  dt_sub <- copy(dt_key)
  dt_sub[, `:=`(
    id_int       = as.integer(factor(stay_id)),
    period_int   = as.integer(trial_num) - 1L,
    a_int        = tx_24h,              # ★ 24h 宽限期治疗分配
    y_int        = as.integer(death_28day_num),
    eligible_int = 1L
  )]
  setorder(dt_sub, id_int, period_int)

  # 按 period 过滤
  dt_sub[, `:=`(n_a = sum(a_int), n_y = sum(y_int)), by = period_int]
  dt_sub <- dt_sub[n_a > 0 & n_a < .N & n_y >= 5]
  if (nrow(dt_sub) < 100) return(NULL)

  # 协变量准备
  covs <- intersect(dag14_covs, names(dt_sub))
  covs <- covs[sapply(covs, function(v) {
    x <- dt_sub[[v]]
    if (is.factor(x))    return(nlevels(x) > 1)
    if (is.character(x)) return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]
  if (length(covs) < 2) return(NULL)
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
# 5. 核心分析：SEQTaRget × 5 imp × 24h Grace Period
# ============================================================================
log_section("SEQTaRget 24h Grace Period × 5 插补")

imp_betas  <- numeric()
imp_ses    <- numeric()
imp_hrs    <- numeric()
imp_ci_l   <- numeric()
imp_ci_u   <- numeric()
imp_run    <- numeric()
imp_nevt   <- numeric()
imp_success <- integer()
imp_fail   <- integer()

for (i_imp in 1L:N_IMPUTE) {
  t_start <- Sys.time()
  log_info(sprintf("--- 插补集 %02d/%02d ---", i_imp, N_IMPUTE))

  parquet_file <- sprintf(MICE_PARQUET_TEMPL, i_imp)
  if (!file.exists(parquet_file)) {
    log_info(sprintf("  [imp%02d] 文件不存在: %s", i_imp, parquet_file))
    imp_fail <- c(imp_fail, i_imp)
    next
  }
  dt_imp <- as.data.table(read_parquet(parquet_file))

  # 保留 DAG14 协变量
  keep_cols <- intersect(c("stay_id", "trial_num", DAG14), names(dt_imp))
  dt_imp <- dt_imp[, ..keep_cols]
  gc(TRUE)

  # Merge 对齐
  covs_in_imp <- intersect(DAG14, names(dt_imp))
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

  # 运行 SEQTaRget (24h grace)
  r <- tryCatch(
    run_seqtarget_grace(dt_comb, DAG14, boot = BOOT),
    error = function(e) NULL
  )

  t_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  if (is.null(r)) {
    log_info(sprintf("  [imp%02d] 失败, %.1f min", i_imp, t_elapsed / 60))
    imp_fail <- c(imp_fail, i_imp)
    next
  }

  imp_success <- c(imp_success, i_imp)
  imp_betas   <- c(imp_betas, log(r$hr))
  imp_ses     <- c(imp_ses, (log(r$ci_upper) - log(r$ci_lower)) / (2 * 1.96))
  imp_hrs     <- c(imp_hrs, r$hr)
  imp_ci_l    <- c(imp_ci_l, r$ci_lower)
  imp_ci_u    <- c(imp_ci_u, r$ci_upper)
  imp_run     <- c(imp_run, r$runtime)
  imp_nevt    <- c(imp_nevt, r$n_events)

  log_info(sprintf("  [imp%02d] HR=%.4f (%.4f-%.4f), runtime=%.0fs, n=%d, events=%d",
                   i_imp, r$hr, r$ci_lower, r$ci_upper, r$runtime, r$n, r$n_events))
}

# 汇总
log_section("24h Grace Period 分析汇总")
N_success <- length(imp_success)
log_info(sprintf("成功插补数: %d/%d", N_success, N_IMPUTE))

if (N_success == 0) {
  log_error("所有分析均失败！")
  stop("全部失败")
}

# Rubin 合并
rubin <- rubin_combine(imp_betas, imp_ses)
log_info(sprintf("Rubin 合并结果 (24h Grace):"))
log_info(sprintf("  HR = %.4f (%.4f-%.4f), p = %.4f",
                 rubin$hr, rubin$ci_lower, rubin$ci_upper, rubin$p))
log_info(sprintf("  beta = %.6f, SE = %.6f, W = %.6f, B = %.6f",
                 rubin$beta, rubin$se, rubin$W, rubin$B))
log_info(sprintf("  HR 范围: %.4f-%.4f",
                 min(imp_hrs), max(imp_hrs)))

# ---- 导出 ----
per_imp <- data.table(
  imp      = imp_success,
  hr       = imp_hrs,
  ci_lower = imp_ci_l,
  ci_upper = imp_ci_u,
  runtime  = imp_run,
  n_events = imp_nevt
)
fwrite(per_imp, file.path(TBL_DIR, "v26_9_sensitivity_grace24h_per_imp.csv"))

rubin_tbl <- data.table(
  analysis = "Grace Period 24h (SEQTaRget ITT, DAG14, boot=50)",
  hr       = rubin$hr,
  ci_lower = rubin$ci_lower,
  ci_upper = rubin$ci_upper,
  p        = rubin$p,
  W        = rubin$W,
  B        = rubin$B
)
fwrite(rubin_tbl, file.path(TBL_DIR, "v26_9_sensitivity_grace24h_rubin.csv"))

log_info("24h Grace Period 敏感性分析完成")
close(log_con)