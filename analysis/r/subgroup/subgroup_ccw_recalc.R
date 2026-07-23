# ============================================================================
# V28 CCW 亚组分析重新计算脚本（SEQTaRget CCW 框架）
#
# 目的：使用 SEQTaRget CCW 扩展后的 MICE 插补数据重新计算 18 项亚组分析
#       替代原 R 脚本中的硬编码值（v26_9_beautify_figures.R:274-295）
#
# 方法：Cox 比例风险模型 + DAG14 + cluster-robust SE (cluster=stay_id)
#       每个亚组单独拟合 Cox 模型
#       18 个亚组 × 20 个 MICE 插补集 → Rubin 规则合并
#
# 数据流：
#   MICE 插补 Parquet (imp01-20) + sequential_trials_imputed_v26.parquet (key_cols)
#     → cbind 合并 → 按亚组分层 Cox + DAG14 + cluster-robust SE
#     → Rubin 规则合并 20 个插补集的 log(HR) 及 SE
#
# 输出：
#   - v26_9_subgroup_ccw_results.csv (主结果)
#   - v26_9_subgroup_ccw_rubin_details.csv (Rubin 合并细节)
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

# DAG14 协变量集
DAG14 <- c(
  "age", "gender", "baseline_sofa", "sepsis3", "charlson",
  "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
  "admission_type", "baseline_phosphate", "lactate_t0", "refeeding_risk",
  "corrected_calcium_mmol_min", "nutrition_support_type"
)

# 18 项预设亚组定义
# 每个 subgroup: name, filter_expr, n_expected (用于校验)
SUBGROUPS <- list(
  list(name = "Phosphate < 0.50 mmol/L", var = "baseline_phosphate", op = "<", value = 0.50),
  list(name = "Phosphate >= 0.50 mmol/L", var = "baseline_phosphate", op = ">=", value = 0.50),
  list(name = "Sepsis-3", var = "sepsis3", op = "==", value = 1),
  list(name = "Non Sepsis-3", var = "sepsis3", op = "==", value = 0),
  list(name = "SOFA >= 8", var = "baseline_sofa", op = ">=", value = 8),
  list(name = "SOFA < 8", var = "baseline_sofa", op = "<", value = 8),
  list(name = "Age >= 65 yr", var = "age", op = ">=", value = 65),
  list(name = "Age < 65 yr", var = "age", op = "<", value = 65),
  list(name = "Male", var = "gender", op = "==", value = "M"),
  list(name = "Female", var = "gender", op = "==", value = "F"),
  list(name = "Invasive ventilation", var = "invasive_vent_at_t0", op = "==", value = 1),
  list(name = "No invasive ventilation", var = "invasive_vent_at_t0", op = "==", value = 0),
  list(name = "Vasoactive agents", var = "vasoactive_at_t0", op = "==", value = 1),
  list(name = "No vasoactive agents", var = "vasoactive_at_t0", op = "==", value = 0),
  list(name = "eGFR < 60", var = "egfr", op = "<", value = 60),
  list(name = "eGFR >= 60", var = "egfr", op = ">=", value = 60),
  list(name = "Refeeding risk", var = "refeeding_risk", op = "==", value = 1),
  list(name = "No refeeding risk", var = "refeeding_risk", op = "==", value = 0)
)

N_IMPUTE <- 20L
SEED     <- 20260629L
set.seed(SEED)

# ============================================================================
# 2. 日志系统
# ============================================================================
log_file <- file.path(LOG_DIR, sprintf("v26_9_subgroup_ccw_recalc_%s.log",
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

log_info("V28 CCW 亚组分析重新计算启动 (SEQTaRget CCW 框架)")
log_info(sprintf("种子: %d, MICE 插补集数: %d, 亚组数: %d", SEED, N_IMPUTE, length(SUBGROUPS)))

# ============================================================================
# 3. 数据加载
# ============================================================================
log_section("数据加载")

dt_full <- as.data.table(read_parquet(DATA_FILE))
dt_full[, `:=`(
  tx_int   = as.integer(itt_group == "treatment"),
  event    = death_28day_num,
  id_int   = as.integer(factor(stay_id))
)]
setorder(dt_full, stay_id, trial_num)

log_info(sprintf("主数据: %d 行, %d 事件, %d 列",
                 nrow(dt_full), sum(dt_full$event), ncol(dt_full)))

# key_cols 必须包含 tx_int
key_cols <- c("stay_id", "subject_id", "trial_num", "time_zero",
              "itt_group", "tx_int", "gender", "admission_type",
              "surv_time", "death_time_from_t0", "dischtime_from_t0",
              "race_group", "icu_type",
              "death_28day_num", "nutrition_support_type")
key_cols <- intersect(key_cols, names(dt_full))
dt_key <- dt_full[, ..key_cols]

char_covariates <- intersect(DAG14, names(dt_full))
char_covariates <- char_covariates[sapply(char_covariates, function(v) {
  is.character(dt_full[[v]])
})]
if (length(char_covariates) > 0) {
  log_info(sprintf("DAG14 字符型变量: %s", paste(char_covariates, collapse = ", ")))
}

rm(dt_full); gc(TRUE)

# ============================================================================
# 4. 辅助函数
# ============================================================================
log_section("辅助函数定义")

# 应用亚组过滤器
apply_subgroup_filter <- function(dt, sg) {
  var <- sg$var
  if (!var %in% names(dt)) return(NULL)

  op <- sg$op
  val <- sg$value

  if (op == "<") {
    dt_sub <- dt[get(var) < val]
  } else if (op == ">=") {
    dt_sub <- dt[get(var) >= val]
  } else if (op == "==") {
    dt_sub <- dt[get(var) == val]
  } else if (op == ">") {
    dt_sub <- dt[get(var) > val]
  } else if (op == "<=") {
    dt_sub <- dt[get(var) <= val]
  } else {
    return(NULL)
  }

  return(dt_sub)
}

# 单亚组 Cox 分析
# 输入：dt_sub (亚组子集), dag14_covs
# 输出：list(beta, se, hr, ci_lower, ci_upper, p_value, n, n_events) 或 NULL
run_subgroup_cox <- function(dt_sub, dag14_covs) {
  # 治疗组和对照组都需要有样本
  if (sum(dt_sub$tx_int == 1) < 10 || sum(dt_sub$tx_int == 0) < 10) return(NULL)

  # 事件数检查
  n_events <- sum(dt_sub$death_28day_num == 1, na.rm = TRUE)
  if (n_events < 20) return(NULL)

  # 字符型 DAG14 变量转因子
  for (v in dag14_covs) {
    if (v %in% names(dt_sub) && is.character(dt_sub[[v]])) {
      dt_sub[, (v) := as.factor(get(v))]
    }
  }

  # 可用协变量检查（亚组子集中需有变异）
  covs <- intersect(dag14_covs, names(dt_sub))
  covs <- covs[sapply(covs, function(v) {
    x <- dt_sub[[v]]
    if (is.factor(x))    return(nlevels(x) > 1)
    if (is.character(x))  return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]

  if (length(covs) < 2) return(NULL)

  formula_str <- sprintf("Surv(surv_time, death_28day_num) ~ tx_int + %s",
                         paste(covs, collapse = " + "))

  fit <- tryCatch(
    coxph(as.formula(formula_str), data = dt_sub, cluster = stay_id),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)
  if (!"tx_int" %in% names(coef(fit))) return(NULL)

  s <- summary(fit)
  beta <- coef(fit)["tx_int"]
  robust_se <- tryCatch(
    s$coefficients["tx_int", "robust se"],
    error = function(e) s$coefficients["tx_int", "se(coef)"]
  )

  return(list(
    beta     = beta,
    se       = robust_se,
    hr       = exp(beta),
    ci_lower = exp(beta - 1.96 * robust_se),
    ci_upper = exp(beta + 1.96 * robust_se),
    p_value  = 2 * pnorm(abs(beta / robust_se), lower.tail = FALSE),
    n        = nrow(dt_sub),
    n_events = n_events
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
    beta = beta_bar, se = se_pool, t = t_stat, df = df, p = p,
    hr = exp(beta_bar),
    ci_lower = exp(beta_bar - 1.96 * se_pool),
    ci_upper = exp(beta_bar + 1.96 * se_pool),
    W = W, B = B
  )
}

log_info("辅助函数定义完成: apply_subgroup_filter, run_subgroup_cox, rubin_combine")

# ============================================================================
# 5. SEQTaRget ×20 亚组分析遍历
# ============================================================================
log_section(sprintf("SEQTaRget CCW ×%d 亚组分析遍历 (DAG14, %d 亚组)", N_IMPUTE, length(SUBGROUPS)))

# 存储所有插补集的所有亚组结果
results_storage <- list()
for (sg in SUBGROUPS) {
  results_storage[[sg$name]] <- list()
}

# 遍历 20 个 MICE 插补集
for (i_imp in seq_len(N_IMPUTE)) {
  log_info(sprintf("\n==== MICE 插补集 %d/%d ====", i_imp, N_IMPUTE))

  parquet_file <- sprintf(MICE_PARQUET_TEMPL, i_imp)
  if (!file.exists(parquet_file)) {
    log_info(sprintf("  [跳过] 文件不存在: %s", parquet_file))
    next
  }

  dt_imp <- as.data.table(read_parquet(parquet_file))
  gc(TRUE)

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

  # 对每个亚组运行 Cox
  for (sg in SUBGROUPS) {
    dt_sub <- apply_subgroup_filter(dt_comb, sg)
    if (is.null(dt_sub) || nrow(dt_sub) < 100) {
      if (i_imp == 1) {
        log_info(sprintf("  [跳过] %s (样本不足)", sg$name))
      }
      next
    }

    r <- tryCatch({
      run_subgroup_cox(dt_sub, DAG14)
    }, error = function(e) NULL)

    if (!is.null(r)) {
      results_storage[[sg$name]][[i_imp]] <- r
      if (i_imp == 1) {
        log_info(sprintf("  imp1 %s: HR=%.4f (%.4f-%.4f), p=%.4f, n=%d, events=%d",
                        sg$name, r$hr, r$ci_lower, r$ci_upper, r$p_value, r$n, r$n_events))
      }
    }

    rm(dt_sub); gc(TRUE)
  }

  rm(dt_comb); gc(TRUE)
}

# ============================================================================
# 6. Rubin 规则合并
# ============================================================================
log_section("Rubin 规则合并")

rubin_results <- list()

for (sg in SUBGROUPS) {
  imp_results <- results_storage[[sg$name]]

  if (length(imp_results) == 0) {
    log_info(sprintf("  [警告] %s 无有效插补集结果", sg$name))
    next
  }

  betas <- sapply(imp_results, function(x) x$beta)
  ses   <- sapply(imp_results, function(x) x$se)

  pooled <- rubin_combine(betas, ses)

  # 取中位数 n 和 n_events
  n_vec <- sapply(imp_results, function(x) x$n)
  n_events_vec <- sapply(imp_results, function(x) x$n_events)

  rubin_results[[sg$name]] <- list(
    name     = sg$name,
    n_imp    = length(imp_results),
    n         = median(n_vec),
    n_events = median(n_events_vec),
    hr        = pooled$hr,
    ci_lower  = pooled$ci_lower,
    ci_upper  = pooled$ci_upper,
    p_value   = pooled$p,
    W         = pooled$W,
    B         = pooled$B,
    df        = pooled$df,
    beta      = pooled$beta,
    se        = pooled$se
  )

  log_info(sprintf("  %s: HR=%.4f (%.4f-%.4f), p=%.4f, n=%s, events=%s, n_imp=%s",
                   sg$name, pooled$hr, pooled$ci_lower, pooled$ci_upper,
                   pooled$p, as.integer(median(n_vec)),
                   as.integer(median(n_events_vec)),
                   as.integer(length(imp_results))))
}

# ============================================================================
# 7. 输出结果
# ============================================================================
log_section("输出结果")

# 主结果 CSV
dt_main <- rbindlist(lapply(rubin_results, function(x) {
  data.table(
    subgroup    = x$name,
    n_imp       = x$n_imp,
    n           = x$n,
    n_events    = x$n_events,
    hr          = round(x$hr, 4),
    ci_lower    = round(x$ci_lower, 4),
    ci_upper    = round(x$ci_upper, 4),
    p_value     = signif(x$p_value, 4),
    significant = ifelse(x$p_value < 0.05, "p<0.05", "NS")
  )
}), fill = TRUE)

main_csv <- file.path(OUT_DIR, "v26_9_subgroup_ccw_results.csv")
fwrite(dt_main, main_csv)
log_info(sprintf("主结果已保存: %s", main_csv))

# Rubin 合并细节 CSV
dt_details <- rbindlist(lapply(rubin_results, function(x) {
  data.table(
    subgroup = x$name,
    n_imp    = x$n_imp,
    n        = x$n,
    n_events = x$n_events,
    beta     = x$beta,
    se       = x$se,
    W        = x$W,
    B        = x$B,
    df       = x$df,
    hr       = x$hr,
    ci_lower = x$ci_lower,
    ci_upper = x$ci_upper,
    p_value  = x$p_value
  )
}), fill = TRUE)

details_csv <- file.path(OUT_DIR, "v26_9_subgroup_ccw_rubin_details.csv")
fwrite(dt_details, details_csv)
log_info(sprintf("Rubin 细节已保存: %s", details_csv))

# 打印结果摘要
log_section("亚组分析结果摘要")
cat("\n")
cat(sprintf("%-35s %8s %15s %10s %10s\n",
            "Subgroup", "HR", "95% CI", "p", "Sig"))
cat(paste(rep("-", 85), collapse = ""), "\n")
for (sg_name in names(rubin_results)) {
  r <- rubin_results[[sg_name]]
  sig <- ifelse(r$p_value < 0.05, "p<0.05", "NS")
  cat(sprintf("%-35s %8.4f %7.4f-%.4f %10.4f %10s\n",
              r$name, r$hr, r$ci_lower, r$ci_upper, r$p_value, sig))
}

# 异质性检验（Cochran's Q）
log_section("异质性检验 (Cochran's Q)")
# 对成对亚组计算 Q 统计量
pair_results <- list()
pairs <- list(
  c("Phosphate < 0.50 mmol/L", "Phosphate >= 0.50 mmol/L"),
  c("Sepsis-3", "Non Sepsis-3"),
  c("SOFA >= 8", "SOFA < 8"),
  c("Age >= 65 yr", "Age < 65 yr"),
  c("Male", "Female"),
  c("Invasive ventilation", "No invasive ventilation"),
  c("Vasoactive agents", "No vasoactive agents"),
  c("eGFR < 60", "eGFR >= 60"),
  c("Refeeding risk", "No refeeding risk")
)

for (pair in pairs) {
  if (all(pair %in% names(rubin_results))) {
    r1 <- rubin_results[[pair[1]]]
    r2 <- rubin_results[[pair[2]]]

    # Q = (beta1 - beta2)^2 / (se1^2 + se2^2)
    q_stat <- (r1$beta - r2$beta)^2 / (r1$se^2 + r2$se^2)
    q_p <- pchisq(q_stat, df = 1, lower.tail = FALSE)

    pair_results[[paste(pair, collapse = " vs ")]] <- list(
      pair = pair,
      hr1 = r1$hr, hr2 = r2$hr,
      q_stat = q_stat,
      q_p = q_p
    )

    log_info(sprintf("  %s vs %s: HR1=%.4f, HR2=%.4f, Q=%.4f, p=%.4f",
                     pair[1], pair[2], r1$hr, r2$hr, q_stat, q_p))
  }
}

# 异质性检验结果 CSV
if (length(pair_results) > 0) {
  dt_het <- rbindlist(lapply(pair_results, function(x) {
    data.table(
      pair = paste(x$pair, collapse = " vs "),
      hr1 = round(x$hr1, 4),
      hr2 = round(x$hr2, 4),
      Q_statistic = round(x$q_stat, 4),
      p_value = signif(x$q_p, 4),
      significant = ifelse(x$q_p < 0.05, "Heterogeneity p<0.05", "NS")
    )
  }))
  het_csv <- file.path(OUT_DIR, "v26_9_subgroup_heterogeneity.csv")
  fwrite(dt_het, het_csv)
  log_info(sprintf("异质性检验已保存: %s", het_csv))
}

# 关闭日志
close(log_con)

cat("\n============================================================\n")
cat("  V28 CCW 亚组分析重新计算完成\n")
cat("============================================================\n")
cat(sprintf("  输出目录: %s\n", OUT_DIR))
cat(sprintf("  主结果:   %s\n", main_csv))
cat(sprintf("  细节:     %s\n", details_csv))
cat(sprintf("  日志:     %s\n", log_file))
