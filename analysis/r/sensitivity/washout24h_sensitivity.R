#!/usr/bin/env Rscript
# ============================================================
# 24h 洗脱期敏感性分析 — SEQTaRget ITT (boot=50)
# 使用 MICE m=5, 仅首期 (trial_num == 1)
# 输出: 合并 HR + 95% CI
# ============================================================

library(SEQTaRget)
library(mice)
library(data.table)
library(arrow)

# ---- 配置 ----
DATA_FILE <- "/data/小米云盘/文档/论文/血磷/TTE/v28/output/washout24h/sequential_trials_washout24h.parquet"
OUTPUT_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v28/output/washout24h"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

# ---- DAG14 协变量 ----
DAG14 <- c(
  "age", "gender", "baseline_sofa", "sepsis3", "charlson",
  "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
  "admission_type", "baseline_phosphate", "lactate_t0", "refeeding_risk",
  "corrected_calcium_mmol_min", "nutrition_support_type"
)

# ---- 加载数据 ----
message("=== 加载数据 ===")
dt <- as.data.table(read_parquet(DATA_FILE))
message(sprintf("  总行数: %d", nrow(dt)))
message(sprintf("  death_28day_num: %d (%.1f%%)", 
    sum(dt$death_28day_num), 100*mean(dt$death_28day_num)))

# 仅首期 (避免 CCW 展开 OOM)
dt <- dt[trial_num == 1]
message(sprintf("  首期行数: %d", nrow(dt)))
message(sprintf("  首期死亡: %d (%.1f%%)", 
    sum(dt$death_28day_num), 100*mean(dt$death_28day_num)))

# SEQTaRget 所需列
dt[, `:=`(
  id_int = as.integer(as.factor(stay_id)),
  period_int = 0L,
  eligible_int = 1L,
  a_int = fifelse(itt_group == "treatment", 1L, 0L),
  y_int = as.integer(death_28day_num)
)]

# 确认可用协变量
avail_covs <- intersect(DAG14, names(dt))
message(sprintf("  可用协变量 (%d): %s", length(avail_covs), 
    paste(avail_covs, collapse = ", ")))

# 因子化
char_cols <- names(dt)[sapply(dt, is.character)]
for (v in intersect(char_cols, avail_covs)) {
  if (v %in% names(dt)) {
    dt[[v]] <- as.factor(dt[[v]])
  }
}

# ---- MICE 多重插补 (m=5) ----
message("\n=== MICE 插补 (m=5) ===")
imp_vars <- avail_covs
imp_df <- dt[, c("id_int", "y_int", "surv_time", imp_vars), with = FALSE]

# 排除零方差变量
keep <- sapply(imp_df, function(x) {
  if (all(is.na(x))) return(FALSE)
  if (is.factor(x) && nlevels(x) <= 1) return(FALSE)
  if (is.numeric(x) && var(x, na.rm = TRUE) == 0) return(FALSE)
  TRUE
})
keep["id_int"] <- FALSE
keep["y_int"] <- FALSE
keep["surv_time"] <- FALSE
imp_vars <- names(keep[keep])
message(sprintf("  MICE 插补变量 (%d): %s", length(imp_vars), 
    paste(imp_vars, collapse = ", ")))

imp <- mice(imp_df[, c("y_int", "surv_time", imp_vars), with = FALSE], 
    m = 5, maxit = 10, seed = 12345, printFlag = FALSE)

# ---- SEQTaRget 分析 ----
message("\n=== SEQTaRget ITT (boot=50) ===")
boot_n <- 50

run_seqtarget <- function(i) {
  dt_imp <- complete(imp, i)
  dt_sub <- dt[, c("id_int", "period_int", "eligible_int", "a_int", "y_int", "surv_time", imp_vars), with = FALSE]
  
  # 合并插补值
  for (v in imp_vars) {
    set(dt_sub, j = v, value = dt_imp[[v]])
  }
  
  # 剔除 NA
  dt_sub <- na.omit(dt_sub)
  if (nrow(dt_sub) < 100) return(NULL)
  
  # 剔除零方差
  covs_use <- imp_vars
  for (v in imp_vars) {
    if (is.factor(dt_sub[[v]]) && nlevels(dt_sub[[v]]) <= 1) {
      covs_use <- setdiff(covs_use, v)
    }
  }
  
  message(sprintf("   插补 %d: N=%d, Events=%d, Covs=%d", 
      i, nrow(dt_sub), sum(dt_sub$y_int), length(covs_use)))
  
  tryCatch({
    res <- SEQuential(
      data = as.data.frame(dt_sub),
      id.col = "id_int",
      time.col = "period_int",
      eligible.col = "eligible_int",
      treatment.col = "a_int",
      outcome.col = "y_int",
      fixed.cols = covs_use,
      time_varying.cols = character(0),
      method = "ITT",
      options = SEQopts(
        bootstrap = TRUE,
        bootstrap.nboot = boot_n,
        bootstrap.CI = 0.95,
        ncores = 1L,
        fastglm.method = 3L,
        hazard = TRUE,
        trial.include = TRUE,
        followup.include = TRUE,
        data.return = FALSE
      ),
      verbose = FALSE
    )
    return(res)
  }, error = function(e) {
    message(sprintf("   插补 %d 失败: %s", i, e$message))
    return(NULL)
  })
}

# 运行 5 个插补
results <- lapply(1:5, run_seqtarget)

# ---- Rubin 规则合并 ----
message("\n=== Rubin 规则合并 ===")
valid <- which(!sapply(results, is.null))
if (length(valid) < 2) {
  message("错误: 有效结果不足 2 个")
  quit(status = 1)
}

# 从 @hazard[[1]] 提取 HR, LCI, UCI
hrs <- sapply(results[valid], function(r) r@hazard[[1]]["Hazard ratio"])
lcis <- sapply(results[valid], function(r) r@hazard[[1]]["LCI"])
ucis <- sapply(results[valid], function(r) r@hazard[[1]]["UCI"])

# 计算 log(HR) 和 se(log(HR))
betas <- log(hrs)
ses <- (log(ucis) - log(lcis)) / (2 * 1.96)
m_actual <- length(valid)

beta_bar <- mean(betas)
W <- mean(ses^2)
B <- var(betas)
total_var <- W + (1 + 1/m_actual) * B
se_total <- sqrt(total_var)
df_rubin <- (m_actual - 1) * (1 + W / ((1 + 1/m_actual) * B))^2
t_val <- qt(0.975, df = df_rubin)
ci_lower <- beta_bar - t_val * se_total
ci_upper <- beta_bar + t_val * se_total
p_val <- 2 * pt(-abs(beta_bar / se_total), df = df_rubin)

hr <- exp(beta_bar)
hr_lower <- exp(ci_lower)
hr_upper <- exp(ci_upper)

message(sprintf("  HR = %.3f (95%% CI %.3f-%.3f)", hr, hr_lower, hr_upper))
message(sprintf("  beta = %.4f, se = %.4f", beta_bar, se_total))
message(sprintf("  p = %.4f", p_val))
message(sprintf("  m_actual = %d", m_actual))

# 各插补结果
message("\n=== 各插补结果 ===")
for (i in seq_along(valid)) {
  message(sprintf("  Imp %d: HR=%.3f (%.3f-%.3f)", 
      valid[i], hrs[i], lcis[i], ucis[i]))
}

# 与主分析比较
message(sprintf("\n=== 与主分析 (48h washout) 比较 ==="))
message(sprintf("  48h washout (首期 DAG14 boot=200): HR = 0.987 (95%% CI 0.893-1.092)"))
message(sprintf("  24h washout (首期 DAG14 boot=50):  HR = %.3f (95%% CI %.3f-%.3f)", hr, hr_lower, hr_upper))

# 保存结果
result_dt <- data.table(
  analysis = "washout24h_sensitivity",
  n = nrow(dt),
  n_events = sum(dt$death_28day_num),
  hr, ci_lower = hr_lower, ci_upper = hr_upper,
  beta = beta_bar, se = se_total, p = p_val,
  m_actual, W, B,
  boot_n = boot_n,
  dag = "DAG14", method = "SEQTaRget_ITT_boot50",
  ref = "48h washout: HR 0.987 (0.893-1.092)"
)
fwrite(result_dt, file.path(OUTPUT_DIR, "washout24h_seqtarget_result.csv"))
message(sprintf("\n结果已保存: %s", file.path(OUTPUT_DIR, "washout24h_seqtarget_result.csv")))