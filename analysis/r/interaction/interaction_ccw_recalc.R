# ============================================================================
# V28 CCW 交互检验重新计算脚本（SEQTaRget CCW 框架）
#
# 目的：使用 SEQTaRget CCW 扩展后的 MICE 插补数据重新计算 10 项交互检验
#       替代原 R 脚本中的硬编码值（v26_18_fix_figures.R:166-214）
#
# 方法：Cox 比例风险模型 + treatment × variable 交互项 + cluster-robust SE (cluster=stay_id)
#       DAG14 协变量集 + 连续变量均值中心化
#       FDR 校正：Benjamini-Hochberg (q<0.10)
#
# 数据流：
#   MICE 插补 Parquet (imp01-20) + sequential_trials_imputed_v26.parquet (key_cols)
#     → cbind 合并 → Cox + 交互项 + cluster-robust SE
#     → Rubin 规则合并 20 个插补集的 log(HR_interaction) 及 SE
#     → FDR 校正
#
# 输出：
#   - v26_9_interaction_ccw_results.csv (主结果)
#   - v26_9_interaction_ccw_rubin_details.csv (Rubin 合并细节)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(survival)
})

# ============================================================================
# 1. 参数配置
# ============================================================================

OUT_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v28/output/tables"
LOG_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v28/output/logs"
for (d in c(OUT_DIR, LOG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# 数据路径
DATA_FILE         <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_imputed_v26.parquet"
MICE_PARQUET_TEMPL <- "/data/小米云盘/文档/论文/血磷/sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet"

# DAG14 协变量集（与主分析一致）
DAG14 <- c(
  "age", "gender", "baseline_sofa", "sepsis3", "charlson",
  "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
  "admission_type", "baseline_phosphate", "lactate_t0", "refeeding_risk",
  "corrected_calcium_mmol_min", "nutrition_support_type"
)

# 10 项预设交互检验变量
INTERACTION_VARS <- c(
  "vasoactive_at_t0", "baseline_sofa", "lactate_t0",
  "invasive_vent_at_t0", "charlson", "age",
  "sepsis3", "baseline_phosphate", "refeeding_risk", "egfr"
)

# 交互检验变量标签（与稿件 Table 8 一致）
INTERACTION_LABELS <- c(
  "vasoactive_at_t0"      = "Vasoactive agents",
  "baseline_sofa"         = "Baseline SOFA (per unit)",
  "lactate_t0"             = "Lactate (per unit)",
  "invasive_vent_at_t0"   = "Invasive ventilation",
  "charlson"              = "Charlson (per unit)",
  "age"                   = "Age (per year)",
  "sepsis3"               = "Sepsis-3",
  "baseline_phosphate"    = "Baseline phosphate (per unit)",
  "refeeding_risk"        = "Refeeding risk",
  "egfr"                  = "eGFR (per unit)"
)

N_IMPUTE <- 20L
SEED     <- 20260629L
set.seed(SEED)

# ============================================================================
# 2. 日志系统
# ============================================================================
log_file <- file.path(LOG_DIR, sprintf("v26_9_interaction_ccw_recalc_%s.log",
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

log_info("V28 CCW 交互检验重新计算启动 (SEQTaRget CCW 框架)")
log_info(sprintf("种子: %d, MICE 插补集数: %d", SEED, N_IMPUTE))
log_info(sprintf("交互检验变量数: %d", length(INTERACTION_VARS)))

# ============================================================================
# 3. 数据加载
# ============================================================================
log_section("数据加载")

# 加载主数据（提取 key_cols + 字符型 DAG14 变量）
dt_full <- as.data.table(read_parquet(DATA_FILE))
dt_full[, `:=`(
  tx_int   = as.integer(itt_group == "treatment"),
  event    = death_28day_num,
  id_int   = as.integer(factor(stay_id))
)]
setorder(dt_full, stay_id, trial_num)

log_info(sprintf("主数据: %d 行, %d 事件, %d 列",
                 nrow(dt_full), sum(dt_full$event), ncol(dt_full)))

# 提取关键列（key_cols 用于与 MICE 插补集合并）
# 必须包含 tx_int（在合并后保留，用于 Cox 公式）
key_cols <- c("stay_id", "subject_id", "trial_num", "time_zero",
              "itt_group", "tx_int", "gender", "admission_type",
              "surv_time", "death_time_from_t0", "dischtime_from_t0",
              "race_group", "icu_type",
              "death_28day_num", "nutrition_support_type")
key_cols <- intersect(key_cols, names(dt_full))
dt_key <- dt_full[, ..key_cols]

# DAG14 在主数据中的字符型变量
char_covariates <- intersect(DAG14, names(dt_full))
char_covariates <- char_covariates[sapply(char_covariates, function(v) {
  is.character(dt_full[[v]])
})]
if (length(char_covariates) > 0) {
  log_info(sprintf("DAG14 字符型变量: %s", paste(char_covariates, collapse = ", ")))
}

# 释放 dt_full
rm(dt_full); gc(TRUE)

# ============================================================================
# 4. 辅助函数：单插补集交互检验
# ============================================================================
log_section("辅助函数定义")

# 单插补集交互检验
# 输入：dt_comb (合并后的数据), var (交互变量)
# 输出：list(beta=log(HR_int), se, hr, ci_lower, ci_upper, p_value) 或 NULL
run_interaction_test <- function(dt_comb, var, dag14_covs) {
  # 数据准备
  dt_sub <- copy(dt_comb)
  dt_sub[, `:=`(
    id_int     = as.integer(factor(stay_id)),
    period_int = as.integer(trial_num) - 1L
  )]

  # 交互变量存在性检查
  if (!var %in% names(dt_sub)) return(NULL)

  # 连续变量均值中心化
  if (is.numeric(dt_sub[[var]])) {
    dt_sub[, (var) := scale(get(var), scale = FALSE)]
  }

  # 字符型 DAG14 变量转因子
  for (v in dag14_covs) {
    if (v %in% names(dt_sub) && is.character(dt_sub[[v]])) {
      dt_sub[, (v) := as.factor(get(v))]
    }
  }

  # 可用协变量检查
  covs <- intersect(dag14_covs, names(dt_sub))
  covs <- covs[sapply(covs, function(v) {
    x <- dt_sub[[v]]
    if (is.factor(x))    return(nlevels(x) > 1)
    if (is.character(x))  return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]

  if (length(covs) < 2) return(NULL)

  # 构建公式：treatment * var + DAG14 (排除交互变量本身避免共线)
  other_covs <- setdiff(covs, var)
  if (length(other_covs) == 0) return(NULL)

  formula_str <- sprintf(
    "Surv(surv_time, death_28day_num) ~ tx_int * %s + %s",
    var, paste(other_covs, collapse = " + ")
  )

  # Cox 回归 + cluster-robust SE
  fit <- tryCatch(
    coxph(as.formula(formula_str), data = dt_sub, cluster = stay_id),
    error = function(e) {
      log_info(sprintf("  [错误] %s Cox 失败: %s", var, e$message))
      return(NULL)
    }
  )

  if (is.null(fit)) return(NULL)

  # 提取交互项系数
  # 交互项命名为 "tx_int:var" 或 "tx_int:varValue"
  coef_names <- names(coef(fit))
  int_pattern <- paste0("^tx_int:", var)
  int_idx <- grep(int_pattern, coef_names)

  if (length(int_idx) == 0) {
    # 尝试反向命名
    int_pattern <- paste0("^", var, ":tx_int")
    int_idx <- grep(int_pattern, coef_names)
  }

  if (length(int_idx) == 0) {
    log_info(sprintf("  [警告] %s 未找到交互项系数", var))
    return(NULL)
  }

  s <- summary(fit)
  int_name <- coef_names[int_idx[1]]

  beta <- coef(fit)[int_name]
  # cluster-robust SE
  robust_se <- tryCatch(
    s$coefficients[int_name, "robust se"],
    error = function(e) s$coefficients[int_name, "se(coef)"]
  )

  return(list(
    variable  = var,
    beta      = beta,
    se        = robust_se,
    hr        = exp(beta),
    ci_lower  = exp(beta - 1.96 * robust_se),
    ci_upper  = exp(beta + 1.96 * robust_se),
    p_value   = 2 * pnorm(abs(beta / robust_se), lower.tail = FALSE)
  ))
}

# Rubin 规则合并
rubin_combine <- function(betas, ses) {
  m <- length(betas)
  beta_bar <- mean(betas)
  W <- mean(ses^2)
  B <- if (m > 1) sum((betas - beta_bar)^2) / (m - 1) else 0
  se_pool <- sqrt(W + (1 + 1 / m) * B)
  t_stat  <- beta_bar / se_pool
  df <- if (m > 1 && B > 0) {
    (m - 1) * (1 + W / ((1 + 1 / m) * B))^2
  } else {
    Inf
  }
  p <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)
  list(
    beta     = beta_bar,
    se       = se_pool,
    t        = t_stat,
    df       = df,
    p        = p,
    hr       = exp(beta_bar),
    ci_lower = exp(beta_bar - 1.96 * se_pool),
    ci_upper = exp(beta_bar + 1.96 * se_pool),
    W        = W,
    B        = B
  )
}

log_info("辅助函数定义完成: run_interaction_test, rubin_combine")

# ============================================================================
# 5. SEQTaRget ×20 交互检验遍历
# ============================================================================
log_section(sprintf("SEQTaRget CCW ×%d 交互检验遍历 (DAG14, 10 变量)", N_IMPUTE))

# 存储所有插补集的所有变量结果
# 结构：results_storage[[var]][[imp_idx]] = list(beta, se, hr, ...)
results_storage <- list()

for (var in INTERACTION_VARS) {
  results_storage[[var]] <- list()
  log_info(sprintf("\n---- 交互检验: %s (%s) ----", var, INTERACTION_LABELS[var]))
}

# 遍历 20 个 MICE 插补集
for (i_imp in seq_len(N_IMPUTE)) {
  log_info(sprintf("\n==== MICE 插补集 %d/%d ====", i_imp, N_IMPUTE))

  # 加载 MICE Parquet
  parquet_file <- sprintf(MICE_PARQUET_TEMPL, i_imp)
  if (!file.exists(parquet_file)) {
    log_info(sprintf("  [跳过] 文件不存在: %s", parquet_file))
    next
  }

  dt_imp <- as.data.table(read_parquet(parquet_file))
  gc(TRUE)

  # 合并 key_cols + DAG14
  keep_cols <- intersect(c("stay_id", "trial_num", DAG14), names(dt_imp))
  dt_imp <- dt_imp[, ..keep_cols]

  dag14_in_dt_imp <- intersect(DAG14, names(dt_imp))
  dag14_remove_from_key <- intersect(dag14_in_dt_imp, names(dt_key))
  by_cols <- c("stay_id", "trial_num")

  dt_comb <- merge(
    dt_key[, setdiff(names(dt_key), dag14_remove_from_key), with = FALSE],
    dt_imp,
    by = by_cols,
    all.x = TRUE,
    sort = FALSE
  )

  rm(dt_imp); gc(TRUE)

  log_info(sprintf("  合并后数据: %d 行, %d 列", nrow(dt_comb), ncol(dt_comb)))

  # 对每个交互变量运行检验
  for (var in INTERACTION_VARS) {
    r <- tryCatch({
      run_interaction_test(dt_comb, var, DAG14)
    }, error = function(e) {
      log_info(sprintf("  [错误] imp%d %s: %s", i_imp, var, e$message))
      NULL
    })

    if (!is.null(r)) {
      results_storage[[var]][[i_imp]] <- r
      if (i_imp == 1) {
        log_info(sprintf("  imp%d %s: HR=%.4f (%.4f-%.4f), p=%.4f",
                        i_imp, var, r$hr, r$ci_lower, r$ci_upper, r$p_value))
      }
    }
  }

  rm(dt_comb); gc(TRUE)
}

# ============================================================================
# 6. Rubin 规则合并 + FDR 校正
# ============================================================================
log_section("Rubin 规则合并 + FDR 校正")

# 存储合并结果
rubin_results <- list()
imp_success_count <- list()

for (var in INTERACTION_VARS) {
  imp_results <- results_storage[[var]]

  if (length(imp_results) == 0) {
    log_info(sprintf("  [警告] %s 无有效插补集结果", var))
    next
  }

  imp_success_count[[var]] <- length(imp_results)

  # 提取 betas 和 ses
  betas <- sapply(imp_results, function(x) x$beta)
  ses   <- sapply(imp_results, function(x) x$se)

  # Rubin 合并
  pooled <- rubin_combine(betas, ses)

  rubin_results[[var]] <- list(
    variable     = var,
    label        = INTERACTION_LABELS[var],
    n_imp        = length(imp_results),
    hr           = pooled$hr,
    ci_lower     = pooled$ci_lower,
    ci_upper     = pooled$ci_upper,
    p_value      = pooled$p,
    W            = pooled$W,
    B            = pooled$B,
    df           = pooled$df,
    beta_pooled  = pooled$beta,
    se_pooled    = pooled$se
  )

  log_info(sprintf("  %s: HR=%.4f (%.4f-%.4f), p=%.4f, n_imp=%d",
                   INTERACTION_LABELS[var], pooled$hr, pooled$ci_lower,
                   pooled$ci_upper, pooled$p, length(imp_results)))
}

# FDR 校正（Benjamini-Hochberg）
p_values <- sapply(rubin_results, function(x) x$p_value)
q_values <- p.adjust(p_values, method = "BH")

for (var in names(rubin_results)) {
  rubin_results[[var]]$q_value <- q_values[var]
}

# ============================================================================
# 7. 输出结果
# ============================================================================
log_section("输出结果")

# 主结果 CSV
dt_main <- rbindlist(lapply(rubin_results, function(x) {
  data.table(
    variable     = x$variable,
    label        = x$label,
    n_imp        = x$n_imp,
    hr           = round(x$hr, 4),
    ci_lower     = round(x$ci_lower, 4),
    ci_upper     = round(x$ci_upper, 4),
    p_value      = signif(x$p_value, 4),
    q_value     = signif(x$q_value, 4),
    significant  = ifelse(x$p_value < 0.05, "Raw p<0.05",
                   ifelse(x$q_value < 0.10, "FDR q<0.10", "NS"))
  )
}))

main_csv <- file.path(OUT_DIR, "v26_9_interaction_ccw_results.csv")
fwrite(dt_main, main_csv)
log_info(sprintf("主结果已保存: %s", main_csv))

# Rubin 合并细节 CSV
dt_details <- rbindlist(lapply(rubin_results, function(x) {
  data.table(
    variable     = x$variable,
    label        = x$label,
    n_imp        = x$n_imp,
    beta_pooled  = x$beta_pooled,
    se_pooled    = x$se_pooled,
    W            = x$W,
    B            = x$B,
    df           = x$df,
    hr           = x$hr,
    ci_lower     = x$ci_lower,
    ci_upper     = x$ci_upper,
    p_value      = x$p_value,
    q_value      = x$q_value
  )
}), fill = TRUE)

details_csv <- file.path(OUT_DIR, "v26_9_interaction_ccw_rubin_details.csv")
fwrite(dt_details, details_csv)
log_info(sprintf("Rubin 细节已保存: %s", details_csv))

# 打印结果摘要
log_section("交互检验结果摘要")
cat("\n")
cat(sprintf("%-35s %8s %15s %10s %10s %10s\n",
            "Variable", "HR", "95% CI", "p", "q", "Sig"))
cat(paste(rep("-", 95), collapse = ""), "\n")
for (var in INTERACTION_VARS) {
  if (var %in% names(rubin_results)) {
    r <- rubin_results[[var]]
    sig <- ifelse(r$p_value < 0.05 & r$q_value < 0.10, "FDR q<0.10",
           ifelse(r$p_value < 0.05, "Raw p<0.05", "NS"))
    cat(sprintf("%-35s %8.4f %7.4f-%.4f %10.4f %10.4f %10s\n",
                r$label, r$hr, r$ci_lower, r$ci_upper, r$p_value, r$q_value, sig))
  }
}

# 关闭日志
close(log_con)

cat("\n============================================================\n")
cat("  V28 CCW 交互检验重新计算完成\n")
cat("============================================================\n")
cat(sprintf("  输出目录: %s\n", OUT_DIR))
cat(sprintf("  主结果:   %s\n", main_csv))
cat(sprintf("  细节:     %s\n", details_csv))
cat(sprintf("  日志:     %s\n", log_file))
