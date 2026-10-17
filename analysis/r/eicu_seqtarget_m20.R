# ============================================================================
# v28.4_eicu_seqtarget_m20.R
# eICU-CRD v2.0 DAG14 外部验证 —— MICE m=20 + SEQTaRget B=50 + Rubin 合并
#
# 与 MIMIC v28.4 主分析对齐的口径：
#   * MICE m=20 (MIMIC 亦为 m=20)，仅插补 4 个原始变量
#   * SEQTaRget trial-period ITT，time.col=自然日历 day
#   * Rubin 合并 logHR 尺度 (同 v28.4_common.R rubin_pool)
#   * B=50 bootstrap（与之前 imp1 结果的 B=50 对齐）
#
# 缺失率审计 → 仅 4 变量需 MICE：
#   lactate_max (48.9%), albumin_min (25.6%), baseline_sofa_t0 (8.6%), egfr (0.2%)
# 其余变量 0–1.5% 缺失，无需插补
# arm/itt_group 缺失 7% (2,226 rows) → 直接 drop（不能插补治疗分配）
#
# 派生逻辑（每个插补集重建）：
#   sepsis3        = suspicion_of_infection (stay层) ∧ baseline_sofa≥2
#   refeeding_risk = hypokalemia | hypomagnesemia
#   corrected_calcium = albumin_min (代理)
#   nutrition_support_type = iv/po 给药路径
#   invasive_vent_at_t0    = sofa_t0.parquet intubated (无缺失)
#
# 版本: V28.4   日期: 2026-10-14
# ============================================================================

suppressPackageStartupMessages({
  library(arrow); library(data.table); library(SEQTaRget)
  library(mice); library(EValue)
})

# ---- 路径 ----
BASE    <- "/data/小米云盘/文档/论文/血磷"
SQL_DIR <- file.path(BASE, "sql/result")
OUT_DIR <- file.path(BASE, "TTE/v28.4/output/tables/eicu_seqtarget")
IMPDIR  <- file.path(SQL_DIR, "eicu_mice_imputed")
dir.create(OUT_DIR, showWarnings=F, recursive=T)
dir.create(IMPDIR, showWarnings=F, recursive=T)

# ---- 参数 ----
SEED <- 20260629L; set.seed(SEED)
TAU  <- 28L; BOOT <- 50L; M <- 20L

# ---- DAG14 ----
DAG14 <- c("age","gender","admission_type","baseline_phosphate",
           "baseline_sofa","charlson","egfr",
           "vasoactive_at_t0","invasive_vent_at_t0","sepsis3",
           "lactate_t0","refeeding_risk","corrected_calcium",
           "nutrition_support_type")

log_info <- function(msg) cat(sprintf("[eICU m=20] %s\n", msg))

# ---- v28.4_common.R 的 rubin_pool 复制 ----
rubin_pool <- function(loghr, se, m=length(loghr)) {
  qbar <- mean(loghr); ubar <- mean(se^2)
  bv <- if (m>1) var(loghr) else 0
  tv <- ubar + bv*(1+1/m)
  dfm <- if (m>1 && bv>0)
    (m-1)*(1 + ubar/(bv*(1+1/m)))^2 else Inf
  zt <- if (is.finite(dfm)) qt(0.975, dfm) else qnorm(0.975)
  list(loghr=qbar, hr=exp(qbar),
       lci=exp(qbar-zt*sqrt(tv)), uci=exp(qbar+zt*sqrt(tv)),
       p=2*(if (is.finite(dfm)) pt(-abs(qbar/sqrt(tv)), dfm)
            else pnorm(-abs(qbar/sqrt(tv)))),
       T=tv, B=bv, df=dfm)
}

zoo_locf <- function(x) {
  if (is.factor(x)) { lev<-levels(x); ii<-zoo_locf(as.integer(x));
    return(factor(ii, levels=seq_along(lev), labels=lev)) }
  nas<-which(is.na(x)); if (length(nas)==0) return(x)
  last<-NA; out<-x
  for (i in seq_along(x)) { if (!is.na(x[i])) last<-x[i]
    else if (!is.na(last)) out[i]<-last }
  out
}

# ============================================================================
# 0. 加载全部源数据
# ============================================================================
bc    <- as.data.table(read_parquet(file.path(SQL_DIR,"eicu_baseline_covariates.parquet")))
sofa0 <- as.data.table(read_parquet(file.path(SQL_DIR,"eicu_sofa_t0.parquet")))
exp   <- as.data.table(read_parquet(file.path(SQL_DIR,"eicu_exposure.parquet")))
tr    <- as.data.table(read_parquet(file.path(SQL_DIR,"eicu_sequential_trials.parquet")))
outc  <- as.data.table(read_parquet(file.path(SQL_DIR,"eicu_outcomes.parquet")))
soi   <- as.data.table(read_parquet(file.path(SQL_DIR,"eicu_suspicion_of_infection.parquet")))
log_info(sprintf("源数据加载完成"))

# ============================================================================
# 1. 合并 trial-level 表 (sofa0 主合并)
# ============================================================================
log_info("Step 1: 合并 trial-level 表")

trials <- sofa0[, .(stay_id=patientunitstayid, tzo=time_zero_offset,
                    sofa_t0, invasive_vent_at_t0=intubated)]
trials <- merge(trials, bc[, .(stay_id=patientunitstayid, tzo=time_zero_offset,
    age, sex, admission_type, baseline_phosphate_mmol_l, baseline_sofa_t0,
    charlson_comorbidity_index, egfr, vasoactive_at_t0, lactate_max,
    hypokalemia, hypomagnesemia, albumin_min)],
    by=c("stay_id","tzo"), all.x=TRUE)
trials <- merge(trials, exp[, .(stay_id=patientunitstayid, tzo=time_zero_offset,
    has_oral=has_oral_ng_phosphate_48h,
    cum_iv=cumulative_iv_phosphate_mmol_48h)],
    by=c("stay_id","tzo"), all.x=TRUE)
trials <- merge(trials, tr[, .(stay_id=patientunitstayid, tzo=time_zero_offset,
    itt_group, time_zero, subject_id, trial_number)],
    by=c("stay_id","tzo"), all.x=TRUE)
trials <- unique(trials, by=c("stay_id","tzo"))
log_info(sprintf("  合并后: %d trials", nrow(trials)))

# drop itt_group==NA (2,226 rows = 7%，合并缺口，不能插补治疗分配)
n_drop <- nrow(trials[is.na(itt_group)])
trials <- trials[!is.na(itt_group)]
log_info(sprintf("  drop itt_group==NA: -%d → %d trials", n_drop, nrow(trials)))

# arm + outcome
trials[, arm := as.integer(itt_group=="treatment")]
outc_s <- unique(outc, by="patientunitstayid")
setnames(outc_s, "patientunitstayid", "stay_id")
trials <- merge(trials, outc_s[, .(stay_id, mortality_28d, time_to_death_days)],
                by="stay_id", all.x=TRUE)
trials[, death_28day_num := as.integer(mortality_28d==1 | mortality_28d==TRUE)]

log_info(sprintf("  arm=1: %d, arm=0: %d, 死亡: %d (%.1f%%)",
                 trials[arm==1,.N], trials[arm==0,.N],
                 sum(trials$death_28day_num,na.rm=T),
                 sum(trials$death_28day_num,na.rm=T)/nrow(trials)*100))

# stay 层共识死亡日（同 v28.4_common.R stay_death_day）
trials[, origin_day := as.integer(tzo / 24)]
trials[, min_date := min(as.Date(time_zero)), by=stay_id]
trials[, origin_day := as.integer(as.Date(time_zero) - min_date)]
trials[, surv_time := time_to_death_days]
dth <- trials[death_28day_num==1L,
              .(cand=pmin(TAU, pmax(0L, origin_day + as.integer(round(surv_time))))),
              by=stay_id]
dth[, death_day := as.integer(names(sort(table(cand), decreasing=T))[1]), by=stay_id]
dth <- dth[!duplicated(stay_id)]
trials <- merge(trials, dth[, .(stay_id, death_day)], by="stay_id", all.x=T)
trials <- trials[is.na(death_day) | origin_day <= death_day]

# suspicion_of_infection stay 层
soi_stay <- soi[, .(sus = max(as.integer(suspected_infection=="true"|suspected_infection==1), na.rm=T)),
                by=patientunitstayid]
setnames(soi_stay, "patientunitstayid", "stay_id")
trials <- merge(trials, soi_stay, by="stay_id", all.x=T)
trials[, soi_flag := as.integer(sus==1)]

log_info(sprintf("  final trials: %d (excl. death<origin)", nrow(trials)))

# ============================================================================
# 2. MICE 插补 (m=20, 仅 4 个原始变量)
# ============================================================================
log_info("Step 2: MICE m=20（仅 4 个变量: lactate_max, albumin_min, baseline_sofa_t0, egfr）")

# baseline_sofa 合并 (sofa_t0 为主，bc 填补，-1 sentinel→NA)
trials[, baseline_sofa_raw := fifelse(!is.na(sofa_t0) & sofa_t0>=0, sofa_t0,
                                      fifelse(!is.na(baseline_sofa_t0) & baseline_sofa_t0>=0,
                                              baseline_sofa_t0, NA_real_))]

# 构建 MICE 输入表：仅 4 个待插补变量 + 辅助变量（arm, outcome）
mice_in <- trials[, .(
  lactate_max, albumin_min, baseline_sofa_raw, egfr,
  arm, death_28day_num, time_to_death_days,
  age, sex, admission_type, baseline_phosphate_mmol_l,
  charlson_comorbidity_index, vasoactive_at_t0, invasive_vent_at_t0,
  hypokalemia, hypomagnesemia, soi_flag, has_oral, cum_iv
)]
setnames(mice_in, c("lactate_max","albumin_min","baseline_sofa_raw","egfr"),
                  c("lactate_max","albumin_min","baseline_sofa","egfr"))

# 分类变量转 factor
mice_in[, sex := factor(sex, levels=c("F","M"))]
mice_in[, admission_type := factor(admission_type)]
mice_in[, vasoactive_at_t0 := as.integer(vasoactive_at_t0)]
mice_in[, invasive_vent_at_t0 := as.integer(invasive_vent_at_t0)]
mice_in[, arm := as.integer(arm)]
mice_in[, death_28day_num := as.integer(death_28day_num)]
mice_in[, hypokalemia := as.integer(hypokalemia)]
mice_in[, hypomagnesemia := as.integer(hypomagnesemia)]
mice_in[, soi_flag := as.integer(soi_flag)]

log_info(sprintf("  MICE input: %d rows, cols=%d", nrow(mice_in), ncol(mice_in)))
log_info(sprintf("  待插补缺失率: lactate=%.1f%%, albumin=%.1f%%, sofa=%.1f%%, egfr=%.1f%%",
                 sum(is.na(mice_in$lactate_max))/nrow(mice_in)*100,
                 sum(is.na(mice_in$albumin_min))/nrow(mice_in)*100,
                 sum(is.na(mice_in$baseline_sofa))/nrow(mice_in)*100,
                 sum(is.na(mice_in$egfr))/nrow(mice_in)*100))

# MICE: method 自动选 (pmm for num, logreg for binomial, polyreg for multi)
t0 <- Sys.time()
mice_model <- mice(mice_in, m=M, maxit=20, seed=SEED, printFlag=FALSE)
t1 <- Sys.time()
log_info(sprintf("  MICE m=20 完成, 耗时 %.1f min", as.numeric(difftime(t1,t0,units="min"))))

# 保存每个插补集
for (i in seq_len(M)) {
  imp <- as.data.table(complete(mice_model, i))
  # complete() 返回 data.frame，arm/outcome 已在里面（未被插补）
  write_parquet(imp, file.path(IMPDIR, sprintf("eicu_mice_imp%02d.parquet", i)))
}
log_info(sprintf("  20 个插补集已保存到 %s", IMPDIR))

# ============================================================================
# 3. 20 个 imp × SEQTaRget B=50 循环
# ============================================================================
log_info("Step 3: 20 imp × SEQTaRget B=50")

extract_log <- function(h0) {
  z <- qnorm(0.975)
  list(hr=as.numeric(h0["Hazard ratio"]),
       lci=as.numeric(h0["LCI"]), uci=as.numeric(h0["UCI"]),
       loghr=log(as.numeric(h0["Hazard ratio"])),
       se=(log(as.numeric(h0["UCI"])) - log(as.numeric(h0["LCI"])))/(2*z))
}

rubin_loghr <- numeric(M); rubin_se <- numeric(M)
all_results <- list()

for (i in seq_len(M)) {
  t_imp0 <- Sys.time()
  imp <- as.data.table(read_parquet(file.path(IMPDIR, sprintf("eicu_mice_imp%02d.parquet", i))))

  # ---- 3a. 构建该 imp 的 trial-level DAG14 ----
  dt <- copy(trials)
  dt[, age := imp$age]
  dt[, gender := factor(imp$sex, levels=c("F","M"))]
  dt[, admission_type := factor(imp$admission_type)]
  dt[, baseline_phosphate := imp$baseline_phosphate_mmol_l]
  dt[, baseline_sofa := imp$baseline_sofa]
  dt[, charlson := imp$charlson_comorbidity_index]
  dt[, egfr := imp$egfr]
  dt[, vasoactive_at_t0 := as.integer(imp$vasoactive_at_t0)]
  dt[, invasive_vent_at_t0 := as.integer(imp$invasive_vent_at_t0)]
  dt[, lactate_t0 := imp$lactate_max]
  dt[, corrected_calcium := imp$albumin_min]
  dt[, refeeding_risk := as.integer(imp$hypokalemia==1 | imp$hypomagnesemia==1)]
  # sepsis3 从 stay 层 soi 重建
  dt[, sepsis3 := as.integer(imp$soi_flag==1 & baseline_sofa>=2)]
  # nutrition_support_type 重建
  dt[, nutrition_support_type := fifelse(
    imp$cum_iv>0 & imp$has_oral>0, "mixed",
    fifelse(imp$cum_iv>0, "parenteral",
      fifelse(imp$has_oral>0, "enteral", "none")))]
  dt[, nutrition_support_type := factor(nutrition_support_type,
      levels=c("none","enteral","parenteral","mixed"))]

  # 确保无 NA
  for (v in DAG14) {
    if (any(is.na(dt[[v]]))) {
      if (is.factor(dt[[v]])) dt[is.na(get(v)), (v) := names(table(dt[[v]]))[1]]
      else dt[is.na(get(v)), (v) := median(dt[[v]], na.rm=T)]
    }
  }

  # ---- 3b. 构建 daily-long 表 (同 v28.4_common.R 框架) ----
  # stay 网格
  grid <- dt[, .(max_o=max(origin_day),
                 dday=fifelse(any(death_28day_num==1L), death_day[1], NA_integer_)),
              by=stay_id]
  grid[, last_day := max_o + TAU]
  # 展开
  daily_list <- mapply(function(s,m,d,l)
    data.table(stay_id=s, day=0L:l, dday=d),
    grid$stay_id, grid$max_o, grid$dday, grid$last_day,
    SIMPLIFY=FALSE)
  daily <- rbindlist(daily_list)
  daily[, outcome := 0L]
  daily[!is.na(dday) & day==dday, outcome := 1L]
  daily[, dday := NULL]

  # eligible + tx_init + DAG14 起点日值
  orig <- unique(dt[, c("stay_id","origin_day","arm",DAG14), with=FALSE],
                 by=c("stay_id","origin_day"))
  setnames(orig, "origin_day", "day")
  daily <- merge(daily, orig, by=c("stay_id","day"), all.x=TRUE)
  daily[!is.na(arm), eligible := 1L]
  daily[eligible==1L, tx_init := arm]
  daily[, arm := NULL]

  # DAG14 LOCF + NOCB (起点日写入已通过 merge 完成)
  setorder(daily, stay_id, day)
  for (v in DAG14) {
    tcol <- dt[[v]]
    if (!is.factor(tcol)) {
      daily[, (v) := zoo_locf(get(v)), by=stay_id]
    } else {
      daily[, (v) := factor(rev(zoo_locf(rev(as.integer(get(v))))),
                            levels=seq_along(levels(tcol)), labels=levels(tcol)), by=stay_id]
    }
    # 全局兜底
    gval <- if (is.factor(tcol)) names(table(tcol))[1] else median(tcol, na.rm=T)
    daily[is.na(get(v)), (v) := gval]
  }
  # 最终 NA 检查
  na_cols <- names(daily)[vapply(daily, function(x) any(is.na(x)), logical(1))]
  if (length(na_cols) > 0) {
    # 紧急兜底：对仍有 NA 的 factor 用 first level, numeric 用 median
    for (v in na_cols) {
      col <- daily[[v]]
      if (is.factor(col)) daily[is.na(get(v)), (v) := names(table(col))[1]]
      else daily[is.na(get(v)), (v) := median(col, na.rm=T)]
    }
    if (i==1) log_info(sprintf("  imp1 仍有 NA 列: %s（已紧急兜底）",
                               paste(na_cols, collapse=",")))
  } else if (i==1) {
    log_info("  imp1 DAG14 无 NA ✓")
  }

  # ---- 3c. SEQTaRget ----
  opts <- SEQopts(hazard=TRUE, km.curves=FALSE, bootstrap=TRUE,
                  bootstrap.nboot=BOOT, bootstrap.sample=1,
                  followup.max=TAU, seed=SEED+i*7L,
                  bootstrap.CI_method="se")
  mod <- SEQuential(
    as.data.frame(daily),
    id.col="stay_id", time.col="day",
    eligible.col="eligible", treatment.col="tx_init",
    outcome.col="outcome",
    time_varying.cols=DAG14, fixed.cols=list(),
    method="ITT", options=opts, verbose=FALSE)

  h <- extract_log(hazard_ratio(mod))
  rubin_loghr[i] <- h$loghr
  rubin_se[i] <- h$se
  all_results[[i]] <- h
  t_imp1 <- Sys.time()
  log_info(sprintf("  imp%02d: HR=%.4f (%.4f–%.4f), %.1fs",
                   i, h$hr, h$lci, h$uci,
                   as.numeric(difftime(t_imp1,t_imp0,units="secs"))))
}

# ============================================================================
# 4. Rubin 合并
# ============================================================================
log_info("Step 4: Rubin 合并 m=20 logHR+SE")

# 逐 imp 点估计范围
hrs <- vapply(all_results, function(x) x$hr, numeric(1))
loghrs <- rubin_loghr; ses <- rubin_se
log_info(sprintf("  逐 imp HR: min=%.4f, max=%.4f, range=%.4f",
                 min(hrs), max(hrs), diff(range(hrs))))
log_info(sprintf("  imputation var B=%.6f, within var Ubar=%.6f",
                 var(loghrs), mean(ses^2)))

pooled <- rubin_pool(loghrs, ses, m=M)
cat("\n=== Rubin 合并结果 ===\n")
cat(sprintf("HR = %.4f (95%% CI %.4f–%.4f), p = %.4f\n",
            pooled$hr, pooled$lci, pooled$uci, pooled$p))
cat(sprintf("Barnard-Rubin df = %.1f\n", pooled$df))
cat(sprintf("  B (between var) = %.6f\n", pooled$B))
cat(sprintf("  T (total var) = %.6f\n", pooled$T))

# E-value
ev_pt <- EValue::evalues.HR(pooled$hr, rare=FALSE)
log_info(sprintf("  E(point)=%.2f, E(bound)=%.2f", ev_pt$est, ev_pt$lo))

# ============================================================================
# 5. 保存 CSV
# ============================================================================
log_info("Step 5: 保存")

# 合并结果
events <- trials[death_28day_num==1, uniqueN(stay_id)]
fwrite(data.table(
  analysis = sprintf("eICU DAG14 (完整14v, m=%d, B=%d, Rubin 合并)", M, BOOT),
  framework = "SEQTaRget trial-period ITT + mice() + Rubin",
  dag = "DAG14 (14 vars, m=20 MICE imputation)",
  m_imputations = M, bootstrap_n = BOOT, seed = SEED,
  n_trial_periods = nrow(trials),
  n_stays = uniqueN(trials$stay_id),
  n_events = events,
  mortality_pct = round(events / uniqueN(trials$stay_id) * 100, 1),
  hr_rubin = pooled$hr,
  lci_rubin = pooled$lci,
  uci_rubin = pooled$uci,
  p_value = pooled$p,
  rubin_df = round(pooled$df, 1),
  rubin_B = pooled$B, rubin_T = pooled$T,
  evalue_point = round(ev_pt$est, 2), evalue_bound = round(ev_pt$lo, 2)
), file.path(OUT_DIR, sprintf("eicu_seqtarget_dag14_m%d_b%d_rubin.csv", M, BOOT)))

# 逐 imp 明细
fwrite(rbindlist(lapply(seq_len(M), function(i)
  data.table(imp=i, hr=all_results[[i]]$hr, lci=all_results[[i]]$lci,
             uci=all_results[[i]]$uci, loghr=all_results[[i]]$loghr,
             se=all_results[[i]]$se))),
  file.path(OUT_DIR, sprintf("eicu_seqtarget_dag14_m%d_b%d_by_imp.csv", M, BOOT)))

log_info(sprintf("DONE → %s", file.path(OUT_DIR, sprintf("eicu_seqtarget_dag14_m%d_b%d_rubin.csv", M, BOOT))))
cat("\n========== 最终对比 ==========\n")
cat(sprintf("eICU DAG14 m=20+B=50+Rubin: HR=%.4f (%.4f–%.4f)\n", pooled$hr, pooled$lci, pooled$uci))
cat(sprintf("eICU DAG14 imp1+B=50 (旧): HR=1.109 (0.828–1.485)\n"))
cat(sprintf("MIMIC DAG14 m=20+B=200+Rubin: HR=0.972 (0.892–1.058)\n"))
