#!/usr/bin/env Rscript
# ============================================================================
# V2-pp 全流程分析管线
# 使用 person-period pooled logistic + sandwich SE
# 替代 SEQTaRget 进行 ITT 主分析、PP 分析、敏感性分析
# 数据源: 最新 MICE 插补集 (v26_mice_imputed_imp01-20.parquet)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(arrow); library(sandwich); library(lmtest)
})

# ============================================================================
# 1. 配置
# ============================================================================
OUT_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9"
TBL_DIR <- file.path(OUT_DIR, "tables")
LOG_DIR <- file.path(OUT_DIR, "logs")
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

MAIN_PARQUET <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_imputed_v26.parquet"
MICE_TEMPLATE <- "/data/小米云盘/文档/论文/血磷/sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet"
N_IMPUTE <- 20L
OUTCOME <- "death_28day_num"  # 28天死亡率

DAG14 <- c("age","gender","baseline_sofa","sepsis3","charlson",
           "egfr","invasive_vent_at_t0","vasoactive_at_t0",
           "admission_type","baseline_phosphate","lactate_t0",
           "refeeding_risk","corrected_calcium_mmol_min",
           "nutrition_support_type")

# 日志
log_file <- file.path(LOG_DIR, sprintf("v2pp_full_pipeline_%s.log",
                                        format(Sys.time(), "%Y%m%d_%H%M%S")))
log_con <- file(log_file, open = "wt")
log_info <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg),
      file = log_con, append = TRUE)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}
log_section <- function(title) {
  cat(sprintf("\n========== %s ==========\n", title),
      file = log_con, append = TRUE)
  cat(sprintf("\n========== %s ==========\n", title))
}

# ============================================================================
# 2. 辅助函数
# ============================================================================

# 加载并合并第 i 个 MICE 插补集
load_imp_data <- function(i_imp) {
  dt_main <- as.data.table(read_parquet(MAIN_PARQUET))
  dt_imp  <- as.data.table(read_parquet(sprintf(MICE_TEMPLATE, i_imp)))

  # 从 main 提取 key 列（结局、治疗、不在 MICE 中的 DAG14）
  dag14_in_imp  <- intersect(DAG14, names(dt_imp))
  dag14_from_main <- setdiff(DAG14, dag14_in_imp)
  key_from_main <- c("stay_id","trial_num","itt_group",
                     "death_28day_num","surv_time", dag14_from_main)
  key_from_main <- intersect(key_from_main, names(dt_main))
  dt_key <- dt_main[, ..key_from_main]
  dt_imp_sub <- dt_imp[, c("stay_id","trial_num", dag14_in_imp), with = FALSE]

  dt <- merge(dt_key, dt_imp_sub, by = c("stay_id", "trial_num"),
              all.x = TRUE, sort = FALSE)
  dt[, a_int := as.integer(itt_group == "treatment")]
  dt[, y_int := as.integer(get(OUTCOME))]
  return(dt)
}

# 准备协变量（排除零方差、因子化）
prepare_covariates <- function(dt, covs) {
  avail <- intersect(covs, names(dt))
  avail <- avail[sapply(avail, function(v) {
    x <- dt[[v]]
    if (is.factor(x))   return(nlevels(x) > 1)
    if (is.character(x)) return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]
  for (v in avail) if (is.character(dt[[v]])) dt[, (v) := as.factor(get(v))]
  return(avail)
}

# V2-pp ITT 拟合
v2pp_itt <- function(dt, covs_use) {
  # period 过滤
  dt[, `:=`(n_a = sum(a_int), n_y = sum(y_int)), by = trial_num]
  dt <- dt[n_a > 0 & n_a < .N & n_y >= 5]

  # person-period 格式
  dt[, `:=`(followup = surv_time, followup_sq = surv_time^2,
            trial = as.numeric(trial_num), trial_sq = as.numeric(trial_num^2))]
  base_f <- "y_int ~ a_int + followup + followup_sq + trial + trial_sq"
  cc_cols <- c("y_int","a_int","followup","followup_sq",
               "trial","trial_sq","stay_id", covs_use)
  cc <- dt[complete.cases(dt[, ..cc_cols])]

  fit <- glm(as.formula(paste(base_f, paste(covs_use, collapse = " + "),
                              sep = " + ")),
             data = cc, family = binomial(), control = list(maxit = 100))
  beta <- as.numeric(coef(fit)["a_int"])
  se   <- as.numeric(sqrt(diag(vcovCL(fit, cluster = cc$stay_id)))["a_int"])
  return(c(beta = beta, se = se, n = nrow(cc), nevt = sum(cc$y_int)))
}

# Rubin 规则合并
rubin_combine <- function(betas, ses) {
  m <- length(betas)
  beta_bar <- mean(betas)
  W <- mean(ses^2)
  B <- if (m > 1) sum((betas - beta_bar)^2) / (m - 1) else 0
  se_pool <- sqrt(W + (1 + 1 / m) * B)
  t_stat  <- beta_bar / se_pool
  df <- if (m > 1) (m - 1) * (1 + W / ((1 + 1 / m) * B))^2 else 1
  p <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)
  return(c(beta = beta_bar, se = se_pool, t = t_stat, df = df, p = p,
           hr = exp(beta_bar),
           ci_l = exp(beta_bar - 1.96 * se_pool),
           ci_u = exp(beta_bar + 1.96 * se_pool),
           W = W, B = B))
}

# ============================================================================
# 3. ITT 主分析 (V2-pp ×20)
# ============================================================================
log_section("ITT 主分析 — V2-pp ×20 Rubin 合并")

time_start <- Sys.time()
imp_results <- list()
betas <- ses <- numeric()

for (i in seq_len(N_IMPUTE)) {
  t0 <- Sys.time()
  log_info(sprintf("  [imp %02d/%02d] 加载数据...", i, N_IMPUTE))
  dt <- load_imp_data(i)
  covs_use <- prepare_covariates(dt, DAG14)
  log_info(sprintf("  [imp %02d] 协变量: %d/%d, 患者: %d, 事件: %d",
                   i, length(covs_use), length(DAG14),
                   uniqueN(dt$stay_id), sum(dt$y_int)))
  
  res <- v2pp_itt(dt, covs_use)
  betas <- c(betas, res["beta"])
  ses   <- c(ses,   res["se"])
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  imp_results[[i]] <- data.frame(
    imp = i, beta = res["beta"], se = res["se"],
    hr = exp(res["beta"]),
    ci_l = exp(res["beta"] - 1.96 * res["se"]),
    ci_u = exp(res["beta"] + 1.96 * res["se"]),
    n = res["n"], nevt = res["nevt"],
    runtime_sec = runtime
  )
  log_info(sprintf("  [imp %02d] HR=%.4f (%.4f-%.4f), β=%+.6f, SE=%.6f, %.1fs",
                   i, exp(res["beta"]), exp(res["beta"]-1.96*res["se"]),
                   exp(res["beta"]+1.96*res["se"]),
                   res["beta"], res["se"], runtime))
}

rubin <- rubin_combine(betas, ses)
total_time <- as.numeric(difftime(Sys.time(), time_start, units = "mins"))

log_info(sprintf("\n===== ITT Rubin 合并结果 ====="))
log_info(sprintf("HR=%.4f (%.4f-%.4f), β=%+.6f, SE=%.6f, p=%.4f",
                 rubin["hr"], rubin["ci_l"], rubin["ci_u"],
                 rubin["beta"], rubin["se"], rubin["p"]))
log_info(sprintf("W=%.8f (组内), B=%.8f (组间)", rubin["W"], rubin["B"]))
log_info(sprintf("总时间: %.1f min (%.1f sec/imp)", total_time,
                 total_time * 60 / N_IMPUTE))

# 保存逐次结果
per_imp_df <- do.call(rbind, imp_results)
fwrite(per_imp_df, file.path(TBL_DIR, "v2pp_itt_per_imp.csv"))
log_info(sprintf("逐次结果已保存: v2pp_itt_per_imp.csv"))

# 保存 Rubin 合并结果
rubin_df <- data.frame(t(rubin))
fwrite(rubin_df, file.path(TBL_DIR, "v2pp_itt_rubin_merged.csv"))
log_info(sprintf("Rubin合并已保存: v2pp_itt_rubin_merged.csv"))

# ============================================================================
# 4. 分段 HR 分析 (0-3, 3-7, 7-14, 14-28 天)
# ============================================================================
log_section("分段 HR 分析 — V2-pp ×20")

PERIODS <- list(c(0, 3), c(3, 7), c(7, 14), c(14, 28))
PERIOD_LABELS <- c("0-3d", "3-7d", "7-14d", "14-28d")

segment_results <- list()

for (seg_i in seq_along(PERIODS)) {
  seg <- PERIODS[[seg_i]]
  seg_label <- PERIOD_LABELS[seg_i]
  log_info(sprintf("  分段 %s...", seg_label))
  
  betas_s <- ses_s <- numeric()
  for (i in seq_len(N_IMPUTE)) {
    dt <- load_imp_data(i)
    dt <- dt[surv_time > seg[1] & surv_time <= seg[2] | 
             (death_28day_num == 1 & surv_time > seg[1])]
    covs_use <- prepare_covariates(dt, DAG14)
    
    # 用 surv_time_in_segment 作为 followup
    dt[, surv_in_seg := pmin(surv_time, seg[2]) - seg[1]]
    dt[surv_in_seg < 0, surv_in_seg := 0]
    dt[, `:=`(followup = surv_in_seg, followup_sq = surv_in_seg^2,
              trial = as.numeric(trial_num), trial_sq = as.numeric(trial_num^2))]
    dt[, `:=`(n_a = sum(a_int), n_y = sum(y_int)), by = trial_num]
    dt <- dt[n_a > 0 & n_a < .N & n_y >= 3]
    
    cc_cols <- c("y_int","a_int","followup","followup_sq",
                 "trial","trial_sq","stay_id", covs_use)
    cc <- dt[complete.cases(dt[, ..cc_cols])]
    base_f <- "y_int ~ a_int + followup + followup_sq + trial + trial_sq"
    
    if (nrow(cc) < 50) { betas_s <- c(betas_s, NA); ses_s <- c(ses_s, NA); next }
    fit <- glm(as.formula(paste(base_f, paste(covs_use, collapse = " + "),
                                sep = " + ")),
               data = cc, family = binomial(), control = list(maxit = 100))
    betas_s <- c(betas_s, as.numeric(coef(fit)["a_int"]))
    ses_s <- c(ses_s, as.numeric(sqrt(diag(vcovCL(fit, cluster=cc$stay_id)))["a_int"]))
  }
  
  rubin_s <- rubin_combine(betas_s[!is.na(betas_s)], ses_s[!is.na(ses_s)])
  segment_results[[seg_i]] <- data.frame(
    segment = seg_label, t(sapply(rubin_s, identity))
  )
  log_info(sprintf("  %s: HR=%.4f (%.4f-%.4f)", seg_label,
                   rubin_s["hr"], rubin_s["ci_l"], rubin_s["ci_u"]))
}

seg_df <- do.call(rbind, segment_results)
fwrite(seg_df, file.path(TBL_DIR, "v2pp_segmented_hr.csv"))
log_info(sprintf("分段HR已保存: v2pp_segmented_hr.csv"))

# ============================================================================
# 5. PP IPTW 分析
# ============================================================================
log_section("PP IPTW 分析 — V2-pp ×20")

# PS 模型拟合 IPTW
v2pp_pp <- function(dt, covs_use) {
  dt <- dt[complete.cases(dt[, c("a_int", covs_use), with = FALSE])]
  ps_formula <- as.formula(paste("a_int ~", paste(covs_use, collapse = " + ")))
  ps_fit <- glm(ps_formula, data = dt, family = binomial(), control = list(maxit = 100))
  dt[, ps := predict(ps_fit, type = "response")]
  dt[, iptw := (a_int / ps) + ((1 - a_int) / (1 - ps))]
  # 稳定化 IPTW
  p_treat <- mean(dt$a_int)
  dt[, iptw_stab := iptw * (a_int * p_treat + (1 - a_int) * (1 - p_treat))]
  # 截断
  q99 <- quantile(dt$iptw_stab, 0.99, na.rm = TRUE)
  q01 <- quantile(dt$iptw_stab, 0.01, na.rm = TRUE)
  dt[, iptw_trunc := pmin(pmax(iptw_stab, q01), q99)]
  
  # Period 过滤
  dt[, `:=`(n_a = sum(a_int), n_y = sum(y_int)), by = trial_num]
  dt <- dt[n_a > 0 & n_a < .N & n_y >= 5]
  
  dt[, `:=`(followup = surv_time, followup_sq = surv_time^2,
            trial = as.numeric(trial_num), trial_sq = as.numeric(trial_num^2))]
  base_f <- "y_int ~ a_int + followup + followup_sq + trial + trial_sq"
  cc_cols <- c("y_int","a_int","followup","followup_sq",
               "trial","trial_sq","stay_id","iptw_trunc", covs_use)
  cc <- dt[complete.cases(dt[, ..cc_cols])]
  
  fit <- glm(as.formula(paste(base_f, paste(covs_use, collapse = " + "),
                              sep = " + ")),
             data = cc, weights = iptw_trunc,
             family = binomial(), control = list(maxit = 100))
  beta <- as.numeric(coef(fit)["a_int"])
  se   <- as.numeric(sqrt(diag(vcovCL(fit, cluster = cc$stay_id)))["a_int"])
  return(c(beta = beta, se = se))
}

betas_pp <- ses_pp <- numeric()
for (i in seq_len(N_IMPUTE)) {
  dt <- load_imp_data(i)
  covs_use <- prepare_covariates(dt, DAG14)
  covs_use <- setdiff(covs_use, "a_int")
  res <- v2pp_pp(dt, covs_use)
  betas_pp <- c(betas_pp, res["beta"])
  ses_pp   <- c(ses_pp,   res["se"])
}

rubin_pp <- rubin_combine(betas_pp, ses_pp)
log_info(sprintf("PP IPTW: HR=%.4f (%.4f-%.4f), β=%+.6f, SE=%.6f",
                 rubin_pp["hr"], rubin_pp["ci_l"], rubin_pp["ci_u"],
                 rubin_pp["beta"], rubin_pp["se"]))
pp_df <- data.frame(t(rubin_pp))
fwrite(pp_df, file.path(TBL_DIR, "v2pp_pp_iptw.csv"))

# ============================================================================
# 6. 汇总
# ============================================================================
log_section("汇总")
cat("\n", file = log_con, append = TRUE)
cat("========================================================\n", file = log_con, append = TRUE)
cat("V2-pp 全流程分析完成\n", file = log_con, append = TRUE)
cat(sprintf("完成时间: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = log_con, append = TRUE)
cat(sprintf("总运行时间: %.1f min\n", total_time), file = log_con, append = TRUE)
cat("========================================================\n", file = log_con, append = TRUE)
cat(sprintf("\nITT:  HR=%.4f (%.4f-%.4f), p=%.4f\n", rubin["hr"], rubin["ci_l"], rubin["ci_u"], rubin["p"]))
cat(sprintf("PP:   HR=%.4f (%.4f-%.4f)\n", rubin_pp["hr"], rubin_pp["ci_l"], rubin_pp["ci_u"]))
cat("\n分段HR:\n")
for (seg_i in seq_along(PERIOD_LABELS)) {
  s <- segment_results[[seg_i]]
  cat(sprintf("  %s: HR=%.4f (%.4f-%.4f)\n", s$segment, s$hr, s$ci_l, s$ci_u))
}

close(log_con)
cat("全流程完成! 结果保存至:", TBL_DIR, "\n")
