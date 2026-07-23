# ============================================================================
# V26.9 统一次要结局 — 分批 CCW 引擎（手动 CCW + 逐期 IPTW/Cox）
#
# 方法学：
#   CCW 克隆-删失-加权框架，手动实现（非 SEQTaRget 包）：
#   1. 对每项结局：克隆扩展（所有 trial_period）
#   2. 逐期拟合 Logistic PS 模型（DAG14）→ 稳定化 IPTW（双侧 1%/99% 截断）
#   3. 逐期拟合加权 Cox 模型（cluster-robust SE via sandwich::vcovCL）
#   4. 固定效应逆方差加权合并各期 HR
#   5. 跨 5 个 MICE 插补集 → Rubin 规则合并
#   6. 主分析 7 项结局 → Bonferroni 校正
#
# 内存友好原因：
#   - 无 bootstrap（SE 来自 cluster-robust Cox，非重抽样）
#   - CCW 扩展后逐期处理，单期数据大小约 1/30 总量
#   - 协变量仅 DAG14（14 个），PS 模型简洁
#   - 不使用 SEQTaRget 包（其 OOM 因 bootstrap 多期叠加）
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(survival)
  library(sandwich)
  library(lmtest)
})

# ============================================================================
# 1. 参数
# ============================================================================

OUT_DIR   <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9"
RDS_DIR   <- file.path(OUT_DIR, "rds")
TBL_DIR   <- file.path(OUT_DIR, "tables")
LOG_DIR   <- file.path(OUT_DIR, "logs")
for (d in c(OUT_DIR, RDS_DIR, TBL_DIR, LOG_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

DATA_FILE          <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_imputed_v26.parquet"
OUTCOME_FILE       <- "/data/小米云盘/文档/论文/血磷/sql/result/M0_07_outcome_extraction.parquet"
AUG_90DAY_FILE     <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_augmented_with_90day.parquet"
MICE_PARQUET_TEMPL <- "/data/小米云盘/文档/论文/血磷/sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet"

DAG14 <- c("age","gender","baseline_sofa","sepsis3","charlson","egfr",
           "invasive_vent_at_t0","vasoactive_at_t0","admission_type",
           "baseline_phosphate","lactate_t0","refeeding_risk",
           "corrected_calcium_mmol_min","nutrition_support_type")

N_IMPUTE <- 5L
SEED     <- 20260629L
set.seed(SEED)
options(mc.cores = 1)

IPTW_TRUNC_LOWER <- 0.01
IPTW_TRUNC_UPPER <- 0.99

OUTCOME_DEFS <- list(
  death_icu   = list(var="death_icu",           time="surv_time",   label="ICU死亡率 (ICU mortality)",              group="primary", src="outcome_data"),
  death_hosp  = list(var="death_hospital",       time="surv_time",   label="住院死亡率 (In-hospital mortality)",      group="primary", src="outcome_data"),
  cv_comp     = list(var="cardiovascular_event_comprehensive", time="surv_time", label="心血管复合结局 (CV composite)", group="primary", src="imputed"),
  aki         = list(var="aki_kdigo2_plus",      time="surv_time",   label="AKI KDIGO ≥2期 (AKI KDIGO ≥2)",           group="primary", src="imputed"),
  transfusion = list(var="transfusion_at_t0",    time="surv_time",   label="输血需求 (Transfusion need)",             group="primary", src="imputed"),
  phos_norm   = list(var="phosphate_normalized", time="surv_time",   label="血磷达标 (Phosphate normalization)",      group="primary", src="outcome_data"),
  death_90d   = list(var="death_90day",          time="surv_time_90d",label="90天死亡率 (90-day mortality)",          group="primary", src="outcome_data"),
  arrhythmia  = list(var="new_arrhythmia_rhythm",time="surv_time",   label="新发心律失常（CV子组件）",                group="sub",   src="imputed"),
  myocardial  = list(var="myocardial_injury",    time="surv_time",   label="心肌损伤（CV子组件）",                    group="sub",   src="imputed"),
  vasopressor = list(var="vasopressor_escalation",time="surv_time",   label="血管活性药升级（CV子组件）",              group="sub",   src="imputed")
)

# ============================================================================
# 2. 日志
# ============================================================================
log_file <- file.path(LOG_DIR, sprintf("v26_9_secondary_unified_batch_ccw_%s.log",
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
log_info("V26.9 统一次要结局 分批 CCW 引擎启动")
log_info(sprintf("种子=%d, 插补集=1-%d, DAG14=%d个, IPTW截断=%.0f%%/%.0f%%",
                 SEED, N_IMPUTE, length(DAG14),
                 (1-IPTW_TRUNC_UPPER)*100, (1-IPTW_TRUNC_LOWER)*100))

# ============================================================================
# 3. 数据加载 — key cols
# ============================================================================
log_section("数据加载")
dt_full <- as.data.table(read_parquet(DATA_FILE))
IMP_OUTCOME_VARS <- c("cardiovascular_event_comprehensive","aki_kdigo2_plus",
                       "transfusion_at_t0","new_arrhythmia_rhythm",
                       "myocardial_injury","vasopressor_escalation","ecg_abnormal")
key_cols <- intersect(c("stay_id","subject_id","trial_num","time_zero","itt_group",
                         "surv_time","death_time_from_t0","dischtime_from_t0",
                         "death_28day_num","nutrition_support_type",
                         IMP_OUTCOME_VARS), names(dt_full))
dt_key <- dt_full[, ..key_cols]
rm(dt_full); gc(TRUE)
log_info(sprintf("主数据: %d 行, key cols=%d", nrow(dt_key), length(key_cols)))

# 合并 outcome_data
dt_outcome <- as.data.table(read_parquet(OUTCOME_FILE))
outcome_cols <- intersect(c("stay_id","time_zero","death_icu","death_hospital",
                             "death_28day","death_90day","vfd_28","aki_stage2_plus",
                             "phosphate_normalized","icu_los_from_t0_days"), names(dt_outcome))
setorder(dt_outcome, stay_id, time_zero)
dt_outcome[, trial_order := rowid(stay_id)]
dt_key[, trial_order := trial_num]
dt_key <- merge(dt_key, dt_outcome[, setdiff(names(dt_outcome),"time_zero"), with=FALSE],
                by=c("stay_id","trial_order"), all.x=TRUE, sort=FALSE)
rm(dt_outcome); gc(TRUE)
log_info(sprintf("outcome_data 合并: %d→%d 行", nrow(dt_key), nrow(dt_key)))

# death_90day 补全
if (file.exists(AUG_90DAY_FILE)) {
  dt_90 <- as.data.table(read_parquet(AUG_90DAY_FILE)); setorder(dt_90,stay_id,time_zero)
  dt_90[, trial_order := rowid(stay_id)]
  dt_90 <- unique(dt_90[, .(stay_id, trial_order, death_90day)])
  if (sum(is.na(dt_key$death_90day)) > 0) {
    dt_key[dt_90, on=c("stay_id","trial_order"),
           death_90day := fifelse(is.na(death_90day), i.death_90day, death_90day)]
  }
  rm(dt_90); gc(TRUE)
}
dt_key[, `:=`(trial_order=NULL, a_int = as.integer(itt_group=="treatment"))]
dt_key[, surv_time_90d := fifelse(death_28day_num==1, surv_time, 90.0)]

# 验证合并结果
for (v in intersect(c("death_icu","death_hospital","death_90day",
                        "vfd_28","aki_stage2_plus","phosphate_normalized",
                        "icu_los_from_t0_days"), names(dt_key))) {
  if (is.numeric(dt_key[[v]]))
    log_info(sprintf("  %s: mean=%.3f, NA=%d", v, mean(dt_key[[v]],na.rm=TRUE), sum(is.na(dt_key[[v]]))))
}

# ============================================================================
# 4. 核心函数
# ============================================================================
log_section("核心函数")

# ---- 4a. CCW 单期扩展 ----
DAG14_COLS <- c("age","gender","baseline_sofa","sepsis3","charlson","egfr",
                "invasive_vent_at_t0","vasoactive_at_t0","admission_type",
                "baseline_phosphate","lactate_t0","refeeding_risk",
                "corrected_calcium_mmol_min","nutrition_support_type")
ccw_expand_period <- function(dt_core, period_int, outcome_var, time_var) {
  clone <- dt_core[trial_num == period_int]
  if (nrow(clone) == 0) return(data.table())
  clone[, `:=`(clone_id = .I, trial_period = period_int)]
  covs_avail <- intersect(DAG14_COLS, names(dt_core))
  sel_cols <- c("stay_id", "trial_period", "clone_id", "a_int",
                outcome_var, time_var, covs_avail)
  sel_cols <- intersect(sel_cols, names(dt_core))
  clone <- clone[, ..sel_cols]
  setnames(clone, "a_int", "tx")
  setnames(clone, outcome_var, "outcome")
  setnames(clone, time_var, "surv")
  clone <- clone[!is.na(outcome) & !is.na(surv) & surv > 0]
  clone
}

# ---- 4b. 逐期拟合 Cox ----
run_period_cox <- function(dt_clone_eval, cov_set, outcome_name) {
  # dt_clone_eval: 单个 trial_period 的克隆数据
  # 剔除零方差协变量
  covs <- intersect(cov_set, names(dt_clone_eval))
  keep_covs <- c()
  for (v in covs) {
    x <- dt_clone_eval[[v]]
    if (is.character(x)) { dt_clone_eval[[v]] <- as.factor(x); x <- dt_clone_eval[[v]] }
    if (is.numeric(x) && sd(x, na.rm=TRUE) > 1e-12) keep_covs <- c(keep_covs, v)
    else if (is.factor(x) && nlevels(x) >= 2) keep_covs <- c(keep_covs, v)
  }
  if (length(keep_covs) < 2) return(NULL)

  # PS 模型
  ps_formula <- as.formula(paste("tx ~", paste(keep_covs, collapse="+")))
  ps_fit <- tryCatch(glm(ps_formula, data=dt_clone_eval, family=binomial(), maxit=100),
                     error=function(e) NULL)
  if (is.null(ps_fit)) return(NULL)

  # IPTW
  dt_clone_eval[, ps := predict(ps_fit, type="response")]
  stabilizer <- mean(dt_clone_eval$tx, na.rm=TRUE)
  dt_clone_eval[, iptw := fifelse(tx==1, stabilizer/ps, (1-stabilizer)/(1-ps))]
  tr_l <- quantile(dt_clone_eval$iptw, IPTW_TRUNC_LOWER, na.rm=TRUE)
  tr_u <- quantile(dt_clone_eval$iptw, IPTW_TRUNC_UPPER, na.rm=TRUE)
  dt_clone_eval[, iptw_trunc := pmin(pmax(iptw, tr_l), tr_u)]

  # 加权 Cox
  cox_formula <- as.formula(paste("Surv(surv, outcome) ~ tx"))
  cox_fit <- tryCatch(coxph(cox_formula, data=dt_clone_eval, weights=iptw_trunc),
                      error=function(e) NULL)
  if (is.null(cox_fit)) return(NULL)

  # cluster-robust SE
  se <- tryCatch(sqrt(diag(vcovCL(cox_fit, cluster=dt_clone_eval$stay_id)))[1],
                 error=function(e) sqrt(diag(vcov(cox_fit)))[1])

  beta <- coef(cox_fit)[1]
  n_ev <- sum(dt_clone_eval$outcome, na.rm=TRUE)
  list(beta=beta, se=se, hr=exp(beta), n_events=n_ev, n=nrow(dt_clone_eval))
}

# ---- 4c. 固定效应逆方差合并 ----
pool_fixed <- function(log_hr, se) {
  weights <- 1/se^2
  ok <- !is.na(weights) & !is.na(log_hr) & is.finite(weights) & weights > 0
  if (sum(ok) < 1) return(NULL)
  pw <- sum(log_hr[ok] * weights[ok]) / sum(weights[ok])
  pw_se <- sqrt(1/sum(weights[ok]))
  list(beta=pw, se=pw_se, hr=exp(pw), ci_lower=exp(pw-1.96*pw_se), ci_upper=exp(pw+1.96*pw_se))
}

# ---- 4d. Rubin 合并 ----
rubin_combine <- function(betas, ses) {
  m <- length(betas); if (m < 2) return(NULL)
  bm <- mean(betas); W <- mean(ses^2); B <- sum((betas-bm)^2)/(m-1)
  se_pool <- sqrt(W + (1+1/m)*B)
  df <- (m-1) * (1 + W/((1+1/m)*B))^2
  p <- 2*pt(abs(bm/se_pool), df=df, lower.tail=FALSE)
  list(beta=bm, se=se_pool, hr=exp(bm),
       ci_lower=exp(bm-1.96*se_pool), ci_upper=exp(bm+1.96*se_pool), p=p, W=W, B=B)
}

log_info("核心函数定义完成: ccw_expand_period, run_period_cox, pool_fixed, rubin_combine")

# ============================================================================
# 5. 遍历结局
# ============================================================================
log_section("分批 CCW 引擎 — 遍历结局")

all_results <- list()
all_imp_results <- list()

for (oname in names(OUTCOME_DEFS)) {
  def  <- OUTCOME_DEFS[[oname]]
  ovar <- def$var; tvar <- def$time; olbl <- def$label; ogrp <- def$group
  log_info(sprintf("\n======== 结局: %s (%s) ========", olbl, oname))

  if (!ovar %in% names(dt_key)) { log_info("  变量不存在，跳过"); next }
  if (!tvar %in% names(dt_key)) { log_info("  时间变量不存在，跳过"); next }

  imp_betas <- imp_ses <- imp_hrs <- imp_n <- imp_ne <- numeric()
  imp_results <- list()

  for (i_imp in seq_len(N_IMPUTE)) {
    t_imp <- Sys.time()
    parquet_file <- sprintf(MICE_PARQUET_TEMPL, i_imp)
    if (!file.exists(parquet_file)) { log_info(sprintf("  [imp%02d] 文件不存在", i_imp)); next }

    dt_imp <- as.data.table(read_parquet(parquet_file))
    dag14_in <- intersect(DAG14, names(dt_imp))
    dt_imp <- dt_imp[, c("stay_id","trial_num", dag14_in), with=FALSE]

    # 合并 key 列
    dag_rem <- intersect(dag14_in, names(dt_key))
    dt_core <- merge(dt_key[, setdiff(names(dt_key), dag_rem), with=FALSE],
                     dt_imp, by=c("stay_id","trial_num"), all.x=TRUE, sort=FALSE)
    rm(dt_imp); gc(TRUE)

    # 获取最大 period
    max_period <- dt_core[, max(trial_num, na.rm=TRUE)]

    # 对每个 trial_period 进行 CCW 扩展 + 建模
    period_log_hr <- period_se <- period_n <- period_ne <- numeric()

    for (tp in seq_len(max_period)) {
      # CCW 扩展：以 tp 为入组期，提取所有当期患者的克隆（1行/患者）
      clone_dt <- ccw_expand_period(dt_core, tp, ovar, tvar)
      if (nrow(clone_dt) < 50) { rm(clone_dt); next }

      cox_res <- run_period_cox(clone_dt, DAG14, oname)
      rm(clone_dt); gc(TRUE)

      if (!is.null(cox_res) && cox_res$n_events >= 5 && is.finite(cox_res$se) && cox_res$se > 0) {
        period_log_hr <- c(period_log_hr, cox_res$beta)
        period_se <- c(period_se, cox_res$se)
        period_n <- c(period_n, cox_res$n)
        period_ne <- c(period_ne, cox_res$n_events)
      }
    }
    rm(dt_core); gc(TRUE)

    # 固定效应合并各期
    period_pool <- pool_fixed(period_log_hr, period_se)
    t_elapsed <- as.numeric(difftime(Sys.time(), t_imp, units="secs"))

    if (!is.null(period_pool) && !is.na(period_pool$hr)) {
      log_info(sprintf("  [imp%02d/%s] HR=%.4f (%.4f–%.4f), periods=%d, events=%.0f, %.1f sec",
                       i_imp, oname, period_pool$hr, period_pool$ci_lower,
                       period_pool$ci_upper, length(period_log_hr), median(period_ne),
                       t_elapsed))
      imp_betas <- c(imp_betas, period_pool$beta)
      imp_ses   <- c(imp_ses,   period_pool$se)
      imp_hrs   <- c(imp_hrs,   period_pool$hr)
      imp_n     <- c(imp_n,     median(period_n))
      imp_ne    <- c(imp_ne,    median(period_ne))
      imp_results[[i_imp]] <- period_pool
    } else {
      log_info(sprintf("  [imp%02d/%s] 合并失败 (periods=%d, events=%.0f), %.1f sec",
                       i_imp, oname, length(period_log_hr), if (length(period_ne)>0) median(period_ne) else 0, t_elapsed))
    }
  }

  # Rubin 合并
  if (length(imp_betas) >= 2) {
    rubin <- rubin_combine(imp_betas, imp_ses)
    log_info(sprintf("  [合并/%s] HR=%.4f (%.4f–%.4f), p=%.6f, m=%d",
                     oname, rubin$hr, rubin$ci_lower, rubin$ci_upper,
                     rubin$p, length(imp_betas)))
    all_results[[oname]] <- list(def=def, rubin=rubin, m=length(imp_betas),
                                  imp_betas=imp_betas, imp_ses=imp_ses, imp_hrs=imp_hrs,
                                  per_imp = imp_results)
  } else {
    log_info(sprintf("  [跳过/%s] 不足2个成功插补集 (m=%d)", oname, length(imp_betas)))
  }
}

# ============================================================================
# 6. Bonferroni 校正 + 输出
# ============================================================================
log_section("Bonferroni 校正")
primary_results <- all_results[sapply(all_results, function(r) r$def$group=="primary")]
N_PRIMARY <- length(primary_results)
BONF_ALPHA <- 0.05 / N_PRIMARY
log_info(sprintf("主分析结局: %d 项, Bonferroni α=%.6f", N_PRIMARY, BONF_ALPHA))

dt_all <- rbindlist(lapply(names(all_results), function(oname) {
  r <- all_results[[oname]]
  data.table(
    outcome     = r$def$label,
    outcome_var = oname,
    group       = r$def$group,
    hr          = r$rubin$hr,
    ci_lower    = r$rubin$ci_lower,
    ci_upper    = r$rubin$ci_upper,
    p_value     = r$rubin$p,
    n_impute    = r$m,
    bonferroni  = if (r$def$group=="primary") fifelse(r$rubin$p < BONF_ALPHA, "Yes", "No") else "N/A (sub)",
    bonf_alpha  = if (r$def$group=="primary") BONF_ALPHA else NA_real_
  )
}))

cat(sprintf("\n%-45s %8s %12s %12s %12s %10s\n",
            "结局","HR","CI下限","CI上限","p值","Bonf."),
    file=log_con, append=TRUE)
for (i in seq_len(nrow(dt_all))) {
  r <- dt_all[i]
  cat(sprintf("  %-43s %8.4f %12.4f %12.4f %12.6f %10s\n",
              r$outcome, r$hr, r$ci_lower, r$ci_upper, r$p_value,
              if (r$group=="sub") "-" else r$bonferroni),
      file=log_con, append=TRUE)
}
fwrite(dt_all, file.path(TBL_DIR, "v26_9_secondary_unified_batch_ccw_results.csv"))

# ============================================================================
# 7. 保存 RDS
# ============================================================================
saveRDS(list(version="V26.9-batch-ccw", timestamp=as.character(Sys.time()),
             seed=SEED, dag14=DAG14, n_impute=N_IMPUTE,
             bonferroni_alpha=BONF_ALPHA, n_primary=N_PRIMARY,
             results=dt_all, per_outcome=all_results),
        file.path(RDS_DIR, "v26_9_secondary_unified_batch_ccw.rds"))

log_info(sprintf("结果CSV: %s", file.path(TBL_DIR, "v26_9_secondary_unified_batch_ccw_results.csv")))
log_info(sprintf("RDS: %s", file.path(RDS_DIR, "v26_9_secondary_unified_batch_ccw.rds")))

close(log_con)
cat("\n===== V26.9 分批 CCW 引擎 完成 =====\n")
cat(sprintf("  结局数: %d 项\n", nrow(dt_all)))
cat(sprintf("  日志: %s\n", log_file))
cat(sprintf("  结果: %s\n", file.path(TBL_DIR, "v26_9_secondary_unified_batch_ccw_results.csv")))
