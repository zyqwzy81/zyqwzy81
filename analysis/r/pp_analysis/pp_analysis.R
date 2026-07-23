# ============================================================================
# V26.9 PP（Per-Protocol）分析
# SEQTaRget ×20 插补 × method="censoring" × Rubin 合并 (DAG14)
#
# 使用方法: SEQTaRget::SEQuential(method = "censoring")
# - ITT (method="ITT"): 保持初始分组，忽略方案偏离
# - PP (method="censoring"): CCCW 三步骤(Cloning→Censoring→IPTW)，
#   剔除方案偏离（治疗开始前死亡/出院/延迟补磷等）
# Bootstrap: 50 次
# 协变量集: DAG14
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(survival)
  library(SEQTaRget)
  library(ggplot2)
})

source("/data/小米云盘/文档/论文/血磷/TTE/v26/code/r/config.R")

# ---- 输出目录 ----
OUT_DIR     <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9"
RDS_DIR     <- file.path(OUT_DIR, "rds")
TBL_DIR     <- file.path(OUT_DIR, "tables")
LOG_DIR     <- file.path(OUT_DIR, "logs")
for (d in c(OUT_DIR, RDS_DIR, TBL_DIR, LOG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

DATA_FILE   <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_imputed_v26.parquet"
BOOT_PP     <- 50L
N_IMPUTE    <- 20L
SEED        <- 20260629L
set.seed(SEED)

DAG14 <- c(
  "age", "gender", "baseline_sofa", "sepsis3", "charlson",
  "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
  "admission_type", "baseline_phosphate", "lactate_t0", "refeeding_risk",
  "corrected_calcium_mmol_min", "nutrition_support_type"
)

MICE_PARQUET_TEMPLATE <- file.path(
  "/data/小米云盘/文档/论文/血磷/sql/result/mice_imputed",
  "v26_mice_imputed_imp%02d.parquet"
)

# ---- 日志系统 ----
log_file <- file.path(LOG_DIR, sprintf("v26_9_pp_analysis_%s.log",
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
log_warn  <- function(msg) log_info(sprintf("[WARN] %s", msg))
log_error <- function(msg) log_info(sprintf("[ERROR] %s", msg))

log_info(sprintf("V26.9 PP 分析启动 — SEQTaRget ×%d method=censoring (DAG14)", N_IMPUTE))
log_info(sprintf("种子: %d, Bootstrap: %d", SEED, BOOT_PP))
log_info(sprintf("输出目录: %s", OUT_DIR))

# ---- 数据加载 ----
log_section("数据加载")
dt_full <- as.data.table(read_parquet(DATA_FILE))
dt_full[, `:=`(
  tx_int   = as.integer(itt_group == "treatment"),
  event    = death_28day_num,
  id_int   = as.integer(factor(stay_id))
)]
setorder(dt_full, id_int, trial_num)
log_info(sprintf("Parquet 数据: %d 行, %d 事件, %d 列",
                 nrow(dt_full), sum(dt_full$event), ncol(dt_full)))

key_cols <- c("stay_id", "subject_id", "trial_num", "time_zero",
              "itt_group", "gender", "admission_type",
              "surv_time", "death_time_from_t0", "dischtime_from_t0",
              "race_group", "icu_type",
              "death_28day_num", "nutrition_support_type")
key_cols <- intersect(key_cols, names(dt_full))
dt_key <- dt_full[, ..key_cols]

# ---- DAG14 可用性检查 ----
mice_p1 <- as.data.table(read_parquet(sprintf(MICE_PARQUET_TEMPLATE, 1L)))
dag14_in_mice   <- intersect(DAG14, names(mice_p1))
dag14_not_mice  <- setdiff(DAG14, names(mice_p1))
log_info(sprintf("DAG14 在 MICE Parquet 中: %s", paste(dag14_in_mice, collapse = ", ")))
log_info(sprintf("DAG14 从主 Parquet 获取: %s", paste(dag14_not_mice, collapse = ", ")))
rm(mice_p1); gc(TRUE)

# ---- SEQTaRget PP 封装函数（method = "censoring", DAG14） ----
run_seqtarget_pp_dag14 <- function(dt_comb, dag14_covs, boot = BOOT_PP) {
  dt_sub <- dt_comb
  # death_28day_num 在 MICE Parquet 中为 category 类型，
  # as.numeric(factor) 会返回 1/2 而非 0/1
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

  # 按 period 过滤：需同时有治疗/对照组、且至少 5 个事件
  dt_sub[, `:=`(n_a = sum(a_int), n_y = sum(y_int)), by = period_int]
  dt_sub <- dt_sub[n_a > 0 & n_a < .N & n_y >= 5]
  if (nrow(dt_sub) < 100) return(NULL)

  # 可用协变量检查
  covs <- intersect(dag14_covs, names(dt_sub))
  covs <- covs[sapply(covs, function(v) {
    x <- dt_sub[[v]]
    if (is.factor(x))   return(nlevels(x) > 1)
    if (is.character(x)) return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]
  if (length(covs) < 2) return(NULL)

  for (v in covs) {
    if (is.character(dt_sub[[v]])) {
      dt_sub[, (v) := as.factor(get(v))]
    }
  }

  # 添加 event 列用于日志统计（不在 dt_comb 中，从 y_int 推导）
  dt_sub[, event := y_int]

  # SEQTaRget 会自动处理无关列（"Non-required columns provided and pruned"）
  # 不手动 prune 列，避免 data.table setdiff 作用域问题
  
  t0 <- Sys.time()
  result <- tryCatch({
    SEQuential(
      data  = as.data.frame(dt_sub),
      id.col          = "id_int",
      time.col        = "period_int",
      eligible.col    = "eligible_int",
      treatment.col   = "a_int",
      outcome.col     = "y_int",
      fixed.cols            = covs,
      time_varying.cols     = character(0),
      method = "censoring",  # ★ PP 分析：clone-censor-weight
      options = SEQopts(
        bootstrap = TRUE,
        bootstrap.nboot = boot,
        bootstrap.CI    = 0.95,
        ncores = 1L,
        fastglm.method = 3L,
        hazard = TRUE,
        trial.include    = TRUE,
        followup.include = TRUE,
        data.return      = FALSE
      ),
      verbose = FALSE
    )
  }, error = function(e) {
    cat(sprintf("[%s]  [PP ERROR] imp=%s: %s\n", format(Sys.time(), "%H:%M:%S"), 
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

# ---- Rubin 规则合并 ----
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
    beta      = beta_bar, se = se_pool, t = t_stat, df = df, p = p,
    hr        = exp(beta_bar),
    ci_lower  = exp(beta_bar - 1.96 * se_pool),
    ci_upper  = exp(beta_bar + 1.96 * se_pool),
    W = W, B = B
  )
}

# ============================================================================
# 主循环：SEQTaRget PP ×20 插补
# ============================================================================
log_section("SEQTaRget PP ×20 插补遍历")
pp_results <- list()
success_count <- 0

for (imp_i in 1L:N_IMPUTE) {
  log_info(sprintf("--- PP 插补集 %02d/%02d ---", imp_i, N_IMPUTE))
  
  # 加载 MICE 插补 Parquet
  mice_file <- sprintf(MICE_PARQUET_TEMPLATE, imp_i)
  if (!file.exists(mice_file)) {
    log_info(sprintf("[WARN] MICE Parquet 不存在: %s", mice_file))
    next
  }
  mice_dt <- as.data.table(read_parquet(mice_file))
  
  # merge 对齐（替代 cbind）按 stay_id + trial_num 主键对齐
  by_cols <- c("stay_id", "trial_num")
  # 移除 dt_key 中也在 mice_dt 出现的所有列（含 itt_group、gender 等）
  overlapped_remove <- intersect(names(mice_dt), names(dt_key))
  overlapped_remove <- setdiff(overlapped_remove, by_cols)  # 保留主键
  dt_key_sub <- dt_key[, setdiff(names(dt_key), overlapped_remove), with = FALSE]
  dt_comb <- merge(
    dt_key_sub,
    mice_dt,
    by = by_cols,
    all.x = TRUE,
    sort = FALSE
  )
  rm(mice_dt, dt_key_sub)
  # 确保关键列存在
  stopifnot("itt_group" %in% names(dt_comb),
            "death_28day_num" %in% names(dt_comb))
  
  # 运行 PP SEQTaRget
  res <- run_seqtarget_pp_dag14(dt_comb, DAG14, boot = BOOT_PP)
  rm(dt_comb); gc(TRUE)
  
  if (is.null(res)) {
    log_info(sprintf("[WARN] 插补 %02d PP 分析返回 NULL", imp_i))
    next
  }
  
  success_count <- success_count + 1
  res$imp <- imp_i
  pp_results[[success_count]] <- res
  log_info(sprintf("  插补 %02d: HR=%.4f, CI=%.4f-%.4f, runtime=%.0fs, n=%d, events=%d",
                   imp_i, res$hr, res$ci_lower, res$ci_upper, res$runtime, res$n, res$n_events))
}

log_section("PP 分析汇总")
log_info(sprintf("成功插补数: %d/%d", success_count, N_IMPUTE))

if (success_count == 0) {
  log_info("[ERROR] 所有 PP 分析均失败！")
  stop("PP分析全部失败")
}

# ---- Rubin 规则合并 ----
betas <- sapply(pp_results, function(x) log(x$hr))
ses   <- sapply(pp_results, function(x) (log(x$ci_upper) - log(x$ci_lower)) / (2 * 1.96))

rubin <- rubin_combine(betas, ses)
log_info(sprintf("Rubin 合并结果:"))
log_info(sprintf("  HR = %.4f (%.4f-%.4f), p = %.4f",
                 rubin$hr, rubin$ci_lower, rubin$ci_upper, rubin$p))
log_info(sprintf("  beta = %.6f, SE = %.6f, W = %.6f, B = %.6f",
                 rubin$beta, rubin$se, rubin$W, rubin$B))
log_info(sprintf("  HR 范围: %.4f-%.4f, SD = %.6f",
                 min(sapply(pp_results, function(x) x$hr)),
                 max(sapply(pp_results, function(x) x$hr)),
                 sd(sapply(pp_results, function(x) x$hr))))

# ---- 导出结果 ----
pp_summary <- data.table(
  imp        = sapply(pp_results, function(x) x$imp),
  hr         = sapply(pp_results, function(x) x$hr),
  ci_lower   = sapply(pp_results, function(x) x$ci_lower),
  ci_upper   = sapply(pp_results, function(x) x$ci_upper),
  runtime    = sapply(pp_results, function(x) x$runtime),
  n          = sapply(pp_results, function(x) x$n),
  n_events   = sapply(pp_results, function(x) x$n_events)
)
fwrite(pp_summary, file.path(TBL_DIR, "v26_9_pp_per_imp_results.csv"))

rubin_tbl <- data.table(
  analysis = "PP CCW+IPTW (DAG14, boot=50)",
  hr       = rubin$hr,
  ci_lower = rubin$ci_lower,
  ci_upper = rubin$ci_upper,
  p        = rubin$p,
  n        = pp_results[[1]]$n,
  n_events = pp_results[[1]]$n_events,
  W        = rubin$W,
  B        = rubin$B,
  hr_min   = min(sapply(pp_results, function(x) x$hr)),
  hr_max   = max(sapply(pp_results, function(x) x$hr))
)
fwrite(rubin_tbl, file.path(TBL_DIR, "v26_9_pp_rubin_results.csv"))

log_section("PP 分析完成")
log_info(sprintf("Per-Protocol (DAG14, method=censoring, boot=%d):",
                 BOOT_PP))
log_info(sprintf("  Rubin HR=%.4f, 95%% CI (%.4f-%.4f), p=%.4f",
                 rubin$hr, rubin$ci_lower, rubin$ci_upper, rubin$p))
log_info(sprintf("  各插补 HR 范围: %.4f-%.4f",
                 min(pp_summary$hr), max(pp_summary$hr)))
close(log_con)
