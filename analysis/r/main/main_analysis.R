# ============================================================================
# V26.9 主分析脚本 — SEQTaRget ×20 Rubin 合并 (DAG14)
#
# 审稿人要求 SEQTaRget 必须在全部 20 个 MICE 插补集上运行，
# 然后通过 Rubin 规则合并 log(HR) 及其标准误。
# 协变量集: DAG14（DAG12 + corrected_calcium_mmol_min + nutrition_support_type）
# Bootstrap: 200 次
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(survival)
  library(SEQTaRget)
  library(ggplot2)
})

# ============================================================================
# 1. 加载配置 + 参数定义
# ============================================================================
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
DATA_FILE   <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_imputed_v26.parquet"

# ---- 分析参数 ----
BOOT_MAIN   <- 200L      # SEQTaRget bootstrap 次数
N_IMPUTE    <- 20L       # MICE 插补集数量
SEED        <- 20260629L # SEQTaRget bootstrap 种子（MICE 用 config.R 的 MICE_SEED=12345）
set.seed(SEED)
options(mc.cores = 1)  # 防止 data.table 和其他包自动启用并行

# ---- DAG14 协变量集（V26.9 — 审稿 M-1.1 响应） ----
# DAG12 + 校正钙 + 营养支持类型
DAG14 <- c(
  "age", "gender", "baseline_sofa", "sepsis3", "charlson",
  "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
  "admission_type", "baseline_phosphate", "lactate_t0", "refeeding_risk",
  "corrected_calcium_mmol_min", "nutrition_support_type"
)
# nutrition_support_type 为字符型变量，在 SEQTaRget 中需要转为 factor

# ============================================================================
# 2. 日志系统
# ============================================================================
log_file <- file.path(LOG_DIR, sprintf("v26_9_main_analysis_%s.log",
                                        format(Sys.time(), "%Y%m%d_%H%M%S")))
log_con  <- file(log_file, open = "wt")

log_info <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg),
      file = log_con, append = TRUE)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
  flush(stdout())  # 强制刷新 stdout 缓冲，避免 nohup 下长时间无输出
}

log_section <- function(title) {
  cat(sprintf("\n========== %s ==========\n", title),
      file = log_con, append = TRUE)
  cat(sprintf("\n========== %s ==========\n", title))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

log_info(sprintf("V26.9 主分析启动 — SEQTaRget ×%d Rubin 合并 (DAG14)", N_IMPUTE))
log_info(sprintf("种子: %d, Bootstrap: %d, ncores: 1", SEED, BOOT_MAIN))
log_info(sprintf("输出目录: %s", OUT_DIR))

# ============================================================================
# 3. 数据加载
# ============================================================================
log_section("数据加载")

# ---- 3a. 加载 Parquet 主数据 ----
dt_full <- as.data.table(read_parquet(DATA_FILE))
dt_full[, `:=`(
  tx_int   = as.integer(itt_group == "treatment"),
  event    = death_28day_num,
  id_int   = as.integer(factor(stay_id))
)]
setorder(dt_full, id_int, trial_num)

log_info(sprintf("Parquet 数据: %d 行, %d 事件, %d 列",
                 nrow(dt_full), sum(dt_full$event), ncol(dt_full)))

# ---- 3b. 提取关键列（用于 cbind 到 MICE 插补集） ----
# 这些列在 MICE 插补过程中保持完整，从 Parquet 直接提取
# 必须包含 nutrition_support_type（DAG14 变量，不在 MICE 中）
key_cols <- c("stay_id", "subject_id", "trial_num", "time_zero",
              "itt_group", "gender", "admission_type",
              "surv_time", "death_time_from_t0", "dischtime_from_t0",
              "race_group", "icu_type",
              "death_28day_num", "nutrition_support_type")
key_cols <- intersect(key_cols, names(dt_full))
dt_key <- dt_full[, ..key_cols]

log_info(sprintf("关键列 (%d 个): %s", length(key_cols),
                 paste(key_cols, collapse = ", ")))

# ---- 3c. MICE Parquet 文件配置 ----
# 使用 sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet
# 这些文件包含 stay_id + trial_num + 所有插补变量
# 比 MICE RDS + complete() 更可靠（无行顺序问题）
MICE_PARQUET_TEMPLATE <- file.path(
  "/data/小米云盘/文档/论文/血磷/sql/result/mice_imputed",
  "v26_mice_imputed_imp%02d.parquet"
)

# ---- 3d. DAG14 可用性检查（从第1插补集 Parquet 检测） ----
dag14_from_config <- DAG14  # 来自 config.R
mice_p1 <- as.data.table(read_parquet(sprintf(MICE_PARQUET_TEMPLATE, 1L)))
dag14_in_mice   <- intersect(dag14_from_config, names(mice_p1))
dag14_not_mice  <- setdiff(dag14_from_config, names(mice_p1))
log_info(sprintf("DAG14 在 MICE Parquet 中: %s",
                 paste(dag14_in_mice, collapse = ", ")))
log_info(sprintf("DAG14 从主 Parquet 获取: %s",
                 paste(dag14_not_mice, collapse = ", ")))
rm(mice_p1); gc(TRUE)

# ---- 3e. 标识 DAG14 中的字符型变量（需转为 factor） ----
char_covariates <- intersect(dag14_from_config, names(dt_full))
char_covariates <- char_covariates[sapply(char_covariates, function(v) {
  is.character(dt_full[[v]])
})]
if (length(char_covariates) > 0) {
  log_info(sprintf("DAG14 中字符型变量（SEQTaRget 需转为 factor）: %s",
                   paste(char_covariates, collapse = ", ")))
}

# ============================================================================
# 4. 辅助函数
# ============================================================================
log_section("辅助函数定义")

# ---- 4a. SEQTaRget ITT 封装函数（DAG14 版） ----
# 输入: dt_comb — cbind 后的数据（key_cols + MICE 插补变量）
#        dag14_covs — DAG14 协变量集
#        boot — bootstrap 次数
# 输出: list(hr, ci_lower, ci_upper, runtime, n, n_events) 或 NULL
run_seqtarget_dag14 <- function(dt_comb, dag14_covs, boot = BOOT_MAIN) {
  dt_sub <- dt_comb  # 不创建副本：dt_comb 在调用后立即被 rm()
  dt_sub[, `:=`(
    id_int       = as.integer(factor(stay_id)),
    period_int   = as.integer(trial_num) - 1L,
    a_int        = as.integer(itt_group == "treatment"),
    y_int        = as.integer(death_28day_num),
    event        = as.integer(death_28day_num),
    eligible_int = 1L
  )]
  setorder(dt_sub, id_int, period_int)

  # 按 period 过滤：需同时有治疗/对照组、且至少 5 个事件
  dt_sub[, `:=`(n_a = sum(a_int), n_y = sum(y_int)), by = period_int]
  dt_sub <- dt_sub[n_a > 0 & n_a < .N & n_y >= 5]
  if (nrow(dt_sub) < 100) {
    return(NULL)
  }

  # 可用协变量检查
  covs <- intersect(dag14_covs, names(dt_sub))
  covs <- covs[sapply(covs, function(v) {
    x <- dt_sub[[v]]
    if (is.factor(x))   return(nlevels(x) > 1)
    if (is.character(x)) return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]

  # 字符型变量转为 factor（SEQTaRget 要求）
  for (v in covs) {
    if (is.character(dt_sub[[v]])) {
      dt_sub[, (v) := as.factor(get(v))]
    }
  }

  if (length(covs) < 2) {
    return(NULL)
  }

  t0 <- Sys.time()
  result <- tryCatch({
    SEQuential(
      data  = as.data.frame(dt_sub),
      id.col       = "id_int",
      time.col     = "period_int",
      eligible.col = "eligible_int",
      treatment.col = "a_int",
      outcome.col   = "y_int",
      fixed.cols          = covs,
      time_varying.cols   = character(0),
      method = "ITT",
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

# ---- 4b. Rubin 规则合并函数 ----
# 输入: betas — log(HR) 向量
#       ses   — log(HR) 标准误向量
# 输出: list(beta, se, t, df, p, hr, ci_lower, ci_upper, W, B)
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
    W         = W,
    B         = B
  )
}

log_info("辅助函数定义完成: run_seqtarget_dag14, rubin_combine")

# ============================================================================
# 5. SEQTaRget ×20 并行插补遍历 (mclapply, 4 cores)
# ============================================================================
log_section(sprintf("SEQTaRget ×%d 并行插补遍历 (DAG14, boot=%d, ncores=4)", N_IMPUTE, BOOT_MAIN))

# 释放 dt_full 节省内存（已提取 dt_key）
rm(dt_full); gc(TRUE)

# ---- 单插补集处理函数 (被 mclapply 调用) ----
process_one_imp <- function(i_imp) {
  setDTthreads(1)  # 防止 data.table 内部并行
  options(mc.cores = 1)

  # 加载 MICE Parquet
  parquet_file <- sprintf(MICE_PARQUET_TEMPLATE, i_imp)
  if (!file.exists(parquet_file)) return(NULL)

  dt_imp <- as.data.table(read_parquet(parquet_file))
  gc(TRUE)

  keep_cols <- intersect(c("stay_id", "trial_num", dag14_from_config), names(dt_imp))
  dt_imp <- dt_imp[, ..keep_cols]

  # merge 对齐
  dag14_in_dt_imp <- intersect(dag14_from_config, names(dt_imp))
  dag14_remove_from_key <- intersect(dag14_in_dt_imp, names(dt_key))
  by_cols <- c("stay_id", "trial_num")

  dt_comb <- merge(
    dt_key[, setdiff(names(dt_key), dag14_remove_from_key), with = FALSE],
    dt_imp,
    by = by_cols,
    all.x = TRUE,
    sort = FALSE
  )

  # 运行 SEQTaRget
  r <- tryCatch({
    run_seqtarget_dag14(dt_comb, dag14_covs = dag14_from_config, boot = BOOT_MAIN)
  }, error = function(e) NULL)

  if (is.null(r)) return(NULL)

  # SE 提取（优先 direct，回退 CI 反算）
  imp_se_val <- if (!is.null(r$se) && is.finite(r$se) && r$se > 0) {
    r$se
  } else if (!is.null(r$se_beta) && is.finite(r$se_beta) && r$se_beta > 0) {
    r$se_beta
  } else {
    (log(r$ci_upper) - log(r$ci_lower)) / (2 * 1.96)
  }
  imp_se_source <- if (!is.null(r$se) && is.finite(r$se) && r$se > 0) "direct" else
                   if (!is.null(r$se_beta) && is.finite(r$se_beta) && r$se_beta > 0) "direct_beta" else "ci_backcalc"

  # 保存检查点
  checkpoint_file <- file.path(RDS_DIR, sprintf("seqtarget_checkpoint_imp%02d.rds", i_imp))
  checkpoint_data <- list(
    imp = i_imp, hr = r$hr, ci_lower = r$ci_lower, ci_upper = r$ci_upper,
    beta = log(r$hr), se = imp_se_val, se_source = imp_se_source,
    runtime_s = r$runtime, n = r$n, n_events = r$n_events,
    timestamp = Sys.time()
  )
  saveRDS(checkpoint_data, checkpoint_file)

  rm(dt_comb, dt_imp, r); gc(TRUE)

  list(
    imp = i_imp, hr = r$hr, ci_lower = r$ci_lower, ci_upper = r$ci_upper,
    beta = log(r$hr), se = imp_se_val, se_source = imp_se_source,
    runtime_s = r$runtime, n = r$n, n_events = r$n_events
  )
}

time_start_all <- Sys.time()

# 检查已有的检查点文件（避免重跑 20 个 SEQTaRget）
existing_cps <- list.files(RDS_DIR, pattern = "seqtarget_checkpoint_imp.*\\.rds", full.names = TRUE)
if (length(existing_cps) >= 2) {
  log_info(sprintf("检测到 %d 个已有检查点文件，跳过 SEQTaRget 重跑，直接从检查点恢复", length(existing_cps)))
  results_list <- NULL
} else {
  # 并行执行（4 核，每个 worker 处理独立插补集）
  log_info("无检查点文件，启动 SEQTaRget ×20 并行计算 (4 cores)...")
  results_list <- parallel::mclapply(
    seq_len(N_IMPUTE),
    process_one_imp,
    mc.cores = 4L,
    mc.preschedule = FALSE
  )
}

time_end_all <- Sys.time()
total_runtime <- as.numeric(difftime(time_end_all, time_start_all, units = "secs"))

# 收集结果（从检查点恢复，绕过 mclapply 返回值序列化问题）
imp_betas      <- numeric()
imp_ses        <- numeric()
imp_se_sources <- character()
imp_hrs        <- numeric()
imp_ci_l       <- numeric()
imp_ci_u       <- numeric()
imp_run        <- numeric()
imp_n          <- integer()
imp_nevt       <- integer()
imp_success    <- integer()
imp_fail       <- integer()
per_imp_results <- list()

for (i_imp in seq_len(N_IMPUTE)) {
  cp_file <- file.path(RDS_DIR, sprintf("seqtarget_checkpoint_imp%02d.rds", i_imp))
  if (!file.exists(cp_file)) {
    log_info(sprintf("  [插补%d] 检查点不存在，视为失败", i_imp))
    imp_fail <- c(imp_fail, i_imp)
    next
  }
  r <- readRDS(cp_file)
  # 检查点文件是列表，直接包含结果字段
  log_info(sprintf("  [插补%d] HR=%.4f (%.4f–%.4f), SE=%.6f [来源:%s], n=%d, events=%d, %.1f min",
                   i_imp, r$hr, r$ci_lower, r$ci_upper,
                   r$se, r$se_source, r$n, r$n_events, r$runtime_s / 60))

  imp_betas <- c(imp_betas, r$beta)
  imp_ses   <- c(imp_ses,   r$se)
  imp_se_sources <- c(imp_se_sources, r$se_source)
  imp_hrs   <- c(imp_hrs,   r$hr)
  imp_ci_l  <- c(imp_ci_l,  r$ci_lower)
  imp_ci_u  <- c(imp_ci_u,  r$ci_upper)
  imp_run   <- c(imp_run,   r$runtime_s)
  imp_n     <- c(imp_n,     r$n)
  imp_nevt  <- c(imp_nevt,  r$n_events)
  imp_success <- c(imp_success, i_imp)

  per_imp_results[[as.character(i_imp)]] <- r
}

# 清理检查点文件
log_info("清理检查点文件...")
old_cps <- list.files(RDS_DIR, pattern = "seqtarget_checkpoint_imp.*\\.rds", full.names = TRUE)
if (length(old_cps) > 0) file.remove(old_cps)

# 遍历汇总
m_actual <- length(imp_success)
log_info(sprintf("\n遍历完成: %d/%d 成功, %d 失败, 总用时 %.1f min",
                 m_actual, N_IMPUTE, length(imp_fail), total_runtime / 60))
if (length(imp_fail) > 0) {
  log_info(sprintf("失败插补集编号: %s", paste(imp_fail, collapse = ", ")))
}

# ============================================================================
# 6. Rubin 规则合并
# ============================================================================
log_section("Rubin 规则合并")

if (m_actual >= 2) {
  rubin <- rubin_combine(imp_betas, imp_ses)

  log_info(sprintf("Rubin 合并结果 (m=%d):", m_actual))
  log_info(sprintf("  HR=%.4f (%.4f–%.4f)", rubin$hr, rubin$ci_lower, rubin$ci_upper))
  log_info(sprintf("  beta=%.6f, se=%.6f, t=%.4f, df=%.1f", rubin$beta, rubin$se, rubin$t, rubin$df))
  log_info(sprintf("  p=%.4f", rubin$p))
  log_info(sprintf("  项内方差 W=%.6f, 项间方差 B=%.6f", rubin$W, rubin$B))

  # 保存合并结果
  results_all <- data.table(
    category   = "SEQTaRget_x20_Rubin_DAG14",
    method     = sprintf("SEQTaRget_ITT_boot%d", BOOT_MAIN),
    hr         = rubin$hr,
    ci_lower   = rubin$ci_lower,
    ci_upper   = rubin$ci_upper,
    beta       = rubin$beta,
    se         = rubin$se,
    p          = rubin$p,
    n          = median(imp_n, na.rm = TRUE),
    n_events   = median(imp_nevt, na.rm = TRUE),
    m_actual   = m_actual,
    W          = rubin$W,
    B          = rubin$B,
    runtime_min = total_runtime / 60
  )
} else {
  log_info("Rubin 合并失败: 不足 2 个成功插补集!")
  results_all <- data.table(
    category = "SEQTaRget_x20_Rubin_DAG14",
    method   = sprintf("SEQTaRget_ITT_boot%d", BOOT_MAIN),
    hr       = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
    beta     = NA_real_, se = NA_real_, p = NA_real_,
    n        = NA_integer_, n_events = NA_integer_,
    m_actual = m_actual,
    W        = NA_real_, B = NA_real_,
    runtime_min = total_runtime / 60
  )
}

# ============================================================================
# 7. HR 分布报告
# ============================================================================
log_section("HR 分布报告")

if (m_actual >= 2) {
  hr_mean   <- mean(imp_hrs)
  hr_median <- median(imp_hrs)
  hr_min    <- min(imp_hrs)
  hr_max    <- max(imp_hrs)
  hr_iqr    <- IQR(imp_hrs)
  hr_sd     <- sd(imp_hrs)

  log_info(sprintf("20 个插补集的 HR 分布 (m=%d):", m_actual))
  log_info(sprintf("  均值   = %.4f", hr_mean))
  log_info(sprintf("  中位数 = %.4f", hr_median))
  log_info(sprintf("  范围   = %.4f – %.4f", hr_min, hr_max))
  log_info(sprintf("  IQR    = %.4f", hr_iqr))
  log_info(sprintf("  SD     = %.4f", hr_sd))

  # ---- 7a. HR 分布箱线图 ----
  dt_hr_dist <- data.table(
    imp     = imp_success,
    hr      = imp_hrs,
    ci_low  = imp_ci_l,
    ci_up   = imp_ci_u,
    runtime = imp_run
  )

  p_box <- ggplot(dt_hr_dist, aes(x = "DAG14", y = hr)) +
    geom_boxplot(fill = "#0072B5", alpha = 0.6, outlier.color = "#BC3C29",
                 outlier.size = 2, width = 0.5) +
    geom_jitter(width = 0.15, size = 2.5, alpha = 0.7, color = "#BC3C29") +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
    geom_hline(yintercept = rubin$hr, linetype = "solid",
               color = "#2E8B57", linewidth = 1.0) +
    annotate("text", x = 0.7, y = rubin$hr,
             label = sprintf("Rubin 合并 HR=%.3f", rubin$hr),
             color = "#2E8B57", hjust = 0, size = 3.8) +
    annotate("text", x = 0.7, y = 1.01,
             label = "HR=1.0 (无效线)",
             color = "grey40", hjust = 0, size = 3.3) +
    labs(
      title    = "SEQTaRget HR 在 20 个 MICE 插补集上的分布",
      subtitle = sprintf("DAG14 · Bootstrap=%d · %d/%d 成功 · 虚线: HR=1.0 · 绿线: Rubin 合并",
                         BOOT_MAIN, m_actual, N_IMPUTE),
      x        = NULL,
      y        = "HR (95% CI)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0, size = 13),
      plot.subtitle = element_text(color = "grey40", hjust = 0, size = 9),
      axis.text.x   = element_blank(),
      axis.ticks.x  = element_blank()
    )

  ggsave(file.path(FIG_DIR, "v26_9_hr_distribution_boxplot.pdf"),
         plot = p_box, width = 7, height = 5, device = "pdf")
  ggsave(file.path(FIG_DIR, "v26_9_hr_distribution_boxplot.png"),
         plot = p_box, width = 7, height = 5, dpi = 300, device = "png")
  log_info(sprintf("HR 分布箱线图已保存: %s", file.path(FIG_DIR, "v26_9_hr_distribution_boxplot.pdf")))
} else {
  log_info("HR 分布报告跳过: 不足 2 个成功插补集")
}

# ============================================================================
# 8. 输出保存
# ============================================================================
log_section("结果保存")

# ---- 8a. 保存合并结果 CSV ----
fwrite(results_all, file.path(TBL_DIR, "v26_9_seqtarget_rubin_results.csv"))
log_info(sprintf("Rubin 合并结果已保存: %s",
                 file.path(TBL_DIR, "v26_9_seqtarget_rubin_results.csv")))

# ---- 8b. 保存逐次结果 CSV ----
if (m_actual >= 2) {
  dt_per_imp <- data.table(
    imp        = imp_success,
    hr         = imp_hrs,
    ci_lower   = imp_ci_l,
    ci_upper   = imp_ci_u,
    beta       = imp_betas,
    se         = imp_ses,
    se_source  = imp_se_sources,
    runtime_s  = imp_run,
    n          = imp_n,
    n_events   = imp_nevt
  )
  fwrite(dt_per_imp, file.path(TBL_DIR, "v26_9_seqtarget_per_imp_results.csv"))
  log_info(sprintf("逐次结果已保存: %s",
                   file.path(TBL_DIR, "v26_9_seqtarget_per_imp_results.csv")))

  # ---- 8c. 保存 RDS ----
  rds_out <- list(
    version        = "V26.9",
    timestamp      = as.character(Sys.time()),
    seed           = SEED,
    dag14          = dag14_from_config,
    dag14_available = dag14_in_mice,
    bootstrap      = BOOT_MAIN,
    n_impute       = N_IMPUTE,
    m_actual       = m_actual,
    imp_success    = imp_success,
    imp_fail       = imp_fail,
    per_imp_betas  = imp_betas,
    per_imp_ses    = imp_ses,
    per_imp_hrs    = imp_hrs,
    rubin_combine  = rubin,
    results_summary = results_all,
    hr_distribution = list(
      mean   = if (m_actual >= 2) mean(imp_hrs) else NA_real_,
      median = if (m_actual >= 2) median(imp_hrs) else NA_real_,
      min    = if (m_actual >= 2) min(imp_hrs) else NA_real_,
      max    = if (m_actual >= 2) max(imp_hrs) else NA_real_,
      iqr    = if (m_actual >= 2) IQR(imp_hrs) else NA_real_,
      sd     = if (m_actual >= 2) sd(imp_hrs) else NA_real_
    ),
    runtime_secs = total_runtime
  )
  saveRDS(rds_out, file.path(RDS_DIR, "v26_9_seqtarget_rubin_full.rds"))
  log_info(sprintf("完整 RDS 已保存: %s",
                   file.path(RDS_DIR, "v26_9_seqtarget_rubin_full.rds")))
}

# ============================================================================
# 9. 流程自审核
# ============================================================================
log_section("流程自审核")

audit <- function(check, pass, detail = "") {
  mark <- if (pass) "PASS" else "FAIL"
  log_info(sprintf("[%s] %s: %s", mark, check, detail))
  return(pass)
}

all_pass <- TRUE
all_pass <- all_pass && audit("Parquet 数据加载",
  exists("dt_full") && nrow(dt_full) > 10000,
  if (exists("dt_full")) sprintf("n=%d", nrow(dt_full)) else "跳过（内存已释放）")
all_pass <- all_pass && audit("MICE Parquet 文件",
  file.exists(sprintf(MICE_PARQUET_TEMPLATE, 1L)) &&
  file.exists(sprintf(MICE_PARQUET_TEMPLATE, N_IMPUTE)),
  sprintf("imp01–imp%02d", N_IMPUTE))
all_pass <- all_pass && audit("DAG14 定义来自 config.R",
  length(dag14_from_config) == 14,
  sprintf("%d 变量", length(dag14_from_config)))
all_pass <- all_pass && audit("SEQTaRget 遍历",
  m_actual >= 2, sprintf("%d/%d 成功", m_actual, N_IMPUTE))
all_pass <- all_pass && audit("Rubin 合并",
  exists("rubin") && !is.na(rubin$hr),
  sprintf("HR=%.4f", if (exists("rubin")) rubin$hr else NA))
all_pass <- all_pass && audit("HR 分布箱线图",
  file.exists(file.path(FIG_DIR, "v26_9_hr_distribution_boxplot.pdf")), "")
all_pass <- all_pass && audit("合并结果 CSV",
  file.exists(file.path(TBL_DIR, "v26_9_seqtarget_rubin_results.csv")), "")
all_pass <- all_pass && audit("逐次结果 CSV",
  file.exists(file.path(TBL_DIR, "v26_9_seqtarget_per_imp_results.csv")), "")
all_pass <- all_pass && audit("完整 RDS 输出",
  file.exists(file.path(RDS_DIR, "v26_9_seqtarget_rubin_full.rds")), "")
all_pass <- all_pass && audit("日志文件",
  file.exists(log_file), "")

log_info(sprintf("\n审核结果: %s", if (all_pass) "全部通过 ✓" else "存在未通过项 ✗"))

# ============================================================================
# 10. 收尾
# ============================================================================
close(log_con)

cat("\n")
cat("################################################################################\n")
cat("              V26.9 主分析完成 — SEQTaRget ×20 Rubin 合并\n")
cat("################################################################################\n\n")

cat(sprintf("DAG14: %s\n", paste(dag14_from_config, collapse = ", ")))
cat(sprintf("Bootstrap: %d\n", BOOT_MAIN))

if (m_actual >= 2) {
  cat(sprintf("成功/总插补集: %d/%d\n", m_actual, N_IMPUTE))
  cat(sprintf("Rubin 合并 HR: %.4f (%.4f–%.4f)\n",
              rubin$hr, rubin$ci_lower, rubin$ci_upper))
  cat(sprintf("Rubin 合并 p 值: %.4f\n", rubin$p))
  cat(sprintf("HR 范围: %.4f–%.4f (IQR=%.4f)\n",
              min(imp_hrs), max(imp_hrs), IQR(imp_hrs)))
}
cat(sprintf("总运行时间: %.1f min\n", total_runtime / 60))
cat(sprintf("\n日志: %s\n", log_file))
cat(sprintf("结果: %s\n", file.path(TBL_DIR, "v26_9_seqtarget_rubin_results.csv")))
cat(sprintf("逐次: %s\n", file.path(TBL_DIR, "v26_9_seqtarget_per_imp_results.csv")))
cat(sprintf("箱线图: %s\n", file.path(FIG_DIR, "v26_9_hr_distribution_boxplot.pdf")))
cat(sprintf("RDS: %s\n", file.path(RDS_DIR, "v26_9_seqtarget_rubin_full.rds")))
cat("========== V26.9 主分析完成 ==========\n")
