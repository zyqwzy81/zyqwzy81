# ============================================================================
# v28.4_eicu_seqtarget_dag14.R
# eICU-CRD v2.0 DAG14 外部验证 —— 完整 14 变量重建 + SEQTaRget trial-period ITT
#
# 与 MIMIC v28.4 主分析使用 **完全相同的框架**：
#   * 自然日历 day 0..28，stay 连续网格（非 trial 序号）
#   * 起点日 eligible=1 + tx_init=arm（纯 ITT，不删失不切换）
#   * DAG14 协变量起点日写入 → LOCF 前向填充 → NOCB 后向填充 → 全局值兜底
#   * SEQTaRget ITT (method="ITT", time_varying.cols=DAG14)
#
# 关键设计决策（已在 v28.4 主稿 Methods/Results 报告）：
#   * 合并顺序：sofa0（primary，覆盖更广起点 + intubated 变量）← bc ← tr ← exp
#   * invasive_vent_at_t0 来自 sofa0.intubated（逐起点 100% 覆盖），
#     而非 bc.ventilation_status_t0（全为 "None" 字符型）
#   * sepsis3 = suspicion_of_infection.parquet（stay 层 suspected_infection）
#               ∧ baseline_sofa≥2（代理 Sepsis-3 脓毒症）
#   * nutrition_support_type = exposure.parquet 的 iv/po 给药路径 → 4 级 factor
#   * corrected_calcium = bc.albumin_min（低白蛋白代理低钙）
#   * refeeding_risk = bc.hypokalemia | bc.hypomagnesemia（电解质异常代理）
#   * outcomes.parquet 多逐起点重复 → unique(by=patientunitstayid) 去重
#
# 版本: V28.4   日期: 2026-10-14
# ============================================================================

suppressPackageStartupMessages({
  library(arrow); library(data.table); library(SEQTaRget); library(EValue)
})

# ---- 路径配置 ----
BASE     <- "/data/小米云盘/文档/论文/血磷"
SQL_DIR  <- file.path(BASE, "sql/result")
OUT_DIR  <- file.path(BASE, "TTE/v28.4/output/tables/eicu_seqtarget")
FIG_DIR  <- file.path(BASE, "TTE/v28.4/output/figures")
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

# ---- 参数 ----
SEED   <- 20260629L
TAU    <- 28L
BOOT   <- 50L         # 与 MIMIC DAG11/DAG10 对齐用 B=50（非主分析的 B=200）
set.seed(SEED)

# ---- DAG14 定义（与 v28.4_common.R 完全一致）----
DAG14 <- c("age","gender","admission_type","baseline_phosphate",
           "baseline_sofa","charlson","egfr",
           "vasoactive_at_t0","invasive_vent_at_t0","sepsis3",
           "lactate_t0","refeeding_risk","corrected_calcium",
           "nutrition_support_type")
DAG10 <- DAG14[!DAG14 %in% c("sepsis3","invasive_vent_at_t0","nutrition_support_type","corrected_calcium")]

log_info <- function(msg) cat(sprintf("[eICU DAG14] %s\n", msg))

# ---- 全局 LOCF 工具（同 v28.4_common.R，定义在 Step 7 之前）----
zoo_locf <- function(x) {
  if (is.factor(x)) {
    lev <- levels(x); ii <- zoo_locf(as.integer(x))
    return(factor(ii, levels=seq_along(lev), labels=lev))
  }
  nas <- which(is.na(x)); if (length(nas)==0) return(x)
  last <- NA; out <- x
  for (i in seq_along(x)) {
    if (!is.na(x[i])) last <- x[i]
    else if (!is.na(last)) out[i] <- last
  }
  out
}

# ============================================================================
# 1. 加载全部 eICU parquet
# ============================================================================
log_info("Step 1: 加载 eICU parquet 文件")

bc   <- as.data.table(read_parquet(file.path(SQL_DIR, "eicu_baseline_covariates.parquet")))
sofa0 <- as.data.table(read_parquet(file.path(SQL_DIR, "eicu_sofa_t0.parquet")))
exp  <- as.data.table(read_parquet(file.path(SQL_DIR, "eicu_exposure.parquet")))
tr   <- as.data.table(read_parquet(file.path(SQL_DIR, "eicu_sequential_trials.parquet")))
outc <- as.data.table(read_parquet(file.path(SQL_DIR, "eicu_outcomes.parquet")))
soi  <- as.data.table(read_parquet(file.path(SQL_DIR, "eicu_suspicion_of_infection.parquet")))

setkey(bc, patientunitstayid, time_zero_offset)
setkey(sofa0, patientunitstayid, time_zero_offset)
setkey(exp, patientunitstayid, time_zero_offset)
setkey(tr, patientunitstayid, time_zero_offset)

log_info(sprintf("  bc: %d rows, sofa0: %d rows, exp: %d rows, tr: %d rows",
                 nrow(bc), nrow(sofa0), nrow(exp), nrow(tr)))

# ============================================================================
# 2. 合并试验起点 + DAG14 协变量（sofa0 为主合并）
# ============================================================================
log_info("Step 2: 合并（sofa0 ← bc ← tr ← exp）")

# ---- 2a. sofa0 主表（覆盖最广起点 + intubated）----
trials <- sofa0[, .(patientunitstayid, time_zero_offset, sofa_t0, intubated)]
setnames(trials, c("patientunitstayid","time_zero_offset","sofa_t0","invasive_vent_at_t0"))
# 确保 intubated 为 0/1 整数
trials[, invasive_vent_at_t0 := as.integer(invasive_vent_at_t0)]

# ---- 2b. 合并 bc ----
bc_sub <- bc[, .(patientunitstayid, time_zero_offset,
                 age, sex, admission_type,
                 baseline_phosphate_mmol_l, baseline_sofa_t0,
                 charlson_comorbidity_index, egfr,
                 vasoactive_at_t0, lactate_max,
                 hypokalemia, hypomagnesemia, albumin_min)]
trials <- merge(trials, bc_sub, by=c("patientunitstayid","time_zero_offset"), all.x=TRUE)

# ---- 2c. 合并 tr（arm、time_zero 日期）----
tr_sub <- tr[, .(patientunitstayid, time_zero_offset, itt_group, time_zero)]
trials <- merge(trials, tr_sub, by=c("patientunitstayid","time_zero_offset"), all.x=TRUE)

# ---- 2d. 合并 exp（nutrition_support_type）----
exp_sub <- exp[, .(patientunitstayid, time_zero_offset,
                   has_oral_ng_phosphate_48h, cum_iv_phosphate_mmol_48h,
                   iv_phosphate_administrations)]
trials <- merge(trials, exp_sub, by=c("patientunitstayid","time_zero_offset"), all.x=TRUE)

log_info(sprintf("  合并后 trials: %d rows", nrow(trials)))

# ---- 2e. 去重 key ----
# (patientunitstayid, time_zero_offset) 应唯一；若重复保留第一行
n_before <- nrow(trials)
trials <- unique(trials, by=c("patientunitstayid","time_zero_offset"))
log_info(sprintf("  去重后: %d rows (移除 %d 重复)", nrow(trials), n_before-nrow(trials)))

# ============================================================================
# 3. 构建 DAG14 变量
# ============================================================================
log_info("Step 3: 构建 DAG14 变量")

# stay_id 统一命名
setnames(trials, "patientunitstayid", "stay_id")

# 日期锚定
trials[, time_zero_date := as.Date(time_zero)]
trials[, origin_day := as.integer(time_zero_offset / 24)]  # 小时 → 天
# 用自然日历日期锚定的 day 0
trials[, min_date := min(time_zero_date), by=stay_id]
trials[, origin_day := as.integer(time_zero_date - min_date)]

# ---- 3a. 二分类变量转换 ----
trials[, arm := as.integer(itt_group == "treatment")]
trials[, invasive_vent_at_t0 := as.integer(invasive_vent_at_t0)]
trials[, vasoactive_at_t0 := as.integer(vasoactive_at_t0)]

# ---- 3b. baseline_sofa（bc + sofa0 互补，-1 sentinel → NA）----
trials[, baseline_sofa := fifelse(!is.na(sofa_t0) & sofa_t0 >= 0, sofa_t0,
                                  fifelse(!is.na(baseline_sofa_t0) & baseline_sofa_t0 >= 0,
                                          baseline_sofa_t0, NA_real_))]

# ---- 3c. gender ----
trials[, gender := fifelse(sex == "M", "M", fifelse(sex == "F", "F", NA_character_))]
trials[, gender := factor(gender, levels=c("F","M"))]

# ---- 3d. admission_type ----
trials[, admission_type := factor(admission_type,
  levels=c("Observation","Emergency","Surgical_Elective","Elective","Emergency Surgical"))]

# ---- 3e. 其他直接映射 ----
trials[, baseline_phosphate := baseline_phosphate_mmol_l]
trials[, charlson := charlson_comorbidity_index]
trials[, egfr := egfr]
trials[, lactate_t0 := lactate_max]
trials[, corrected_calcium := albumin_min]

# ---- 3f. refeeding_risk（hypoK | hypoMg）----
trials[, refeeding_risk := as.integer(hypokalemia == 1 | hypomagnesemia == 1)]

# ---- 3g. sepsis3（stay 层代理: suspected_infection ∧ baseline_sofa≥2）----
soi_stay <- soi[, .(suspected_infection = max(as.integer(suspected_infection == "true" | suspected_infection == 1), na.rm=TRUE)),
                by=patientunitstayid]
setnames(soi_stay, "patientunitstayid", "stay_id")
trials <- merge(trials, soi_stay, by="stay_id", all.x=TRUE)
trials[, sepsis3 := as.integer(suspected_infection == 1 & baseline_sofa >= 2)]
trials[is.na(sepsis3), sepsis3 := 0L]

# ---- 3h. nutrition_support_type（4 级 factor）----
trials[, nutrition_support_type := fifelse(
  !is.na(cum_iv_phosphate_mmol_48h) & cum_iv_phosphate_mmol_48h > 0 &
    !is.na(has_oral_ng_phosphate_48h) & has_oral_ng_phosphate_48h > 0, "mixed",
  fifelse(!is.na(cum_iv_phosphate_mmol_48h) & cum_iv_phosphate_mmol_48h > 0, "parenteral",
    fifelse(!is.na(has_oral_ng_phosphate_48h) & has_oral_ng_phosphate_48h > 0, "enteral", "none")))]
trials[, nutrition_support_type := factor(nutrition_support_type,
  levels=c("none","enteral","parenteral","mixed"))]

# ---- 3i. 清理中间列 ----
trials[, c("sofa_t0","bc_sofa","sex","baseline_sofa_t0","baseline_phosphate_mmol_l",
           "charlson_comorbidity_index","albumin_min","lactate_max",
           "hypokalemia","hypomagnesemia","suspected_infection",
           "has_oral_ng_phosphate_48h","cum_iv_phosphate_mmol_48h",
           "iv_phosphate_administrations","time_zero_offset","itt_group",
           "time_zero_date","min_date","sofa_t0") := NULL]

log_info(sprintf("  DAG14 构建完成，trials: %d rows", nrow(trials)))

# ============================================================================
# 4. 合并结局（outcomes.parquet 多逐起点行 → unique by stay_id）
# ============================================================================
log_info("Step 4: 合并结局")

outc_stay <- unique(outc, by="patientunitstayid")
setnames(outc_stay, "patientunitstayid", "stay_id")
outc_stay <- outc_stay[, .(stay_id, mortality_28d, time_to_death_days)]

trials <- merge(trials, outc_stay, by="stay_id", all.x=TRUE)
trials[, death_28day_num := as.integer(mortality_28d == 1 | mortality_28d == TRUE)]

log_info(sprintf("  events: %d / %d stays (%.1f%%)",
                 sum(trials$death_28day_num, na.rm=TRUE),
                 uniqueN(trials$stay_id),
                 sum(trials$death_28day_num, na.rm=TRUE)/uniqueN(trials$stay_id)*100))

# ============================================================================
# 5. 构建 stay 层共识死亡日（与 v28.4_common.R stay_death_day 同逻辑）
# ============================================================================
log_info("Step 5: 共识死亡日")

trials[, surv_time := time_to_death_days]

dth <- trials[death_28day_num == 1L,
              .(cand = pmin(TAU, pmax(0L, origin_day + as.integer(round(surv_time))))),
              by=stay_id]
dth[, death_day := as.integer(names(sort(table(cand), decreasing=TRUE))[1]), by=stay_id]
trials <- merge(trials, dth[, .(stay_id, death_day)], by="stay_id", all.x=TRUE)

# ============================================================================
# 6. 构建 daily-long person-period 表（自然日历 day 0..28，同 v28.4_common.R）
# ============================================================================
log_info("Step 6: 构建 daily-long 表")

trials <- trials[is.na(death_day) | origin_day <= death_day]

# 每个 stay 的网格长度
grid <- trials[, .(max_origin = max(origin_day),
                   dday = fifelse(any(death_28day_num == 1L),
                                  death_day[1], NA_integer_)),
                by=stay_id]
grid[, last_day := max_origin + TAU]

# 展开网格（逐 stay 循环，避免 vecseq cartesian）
expand_one <- function(stay_id, max_o, dday, last_d) {
  days <- 0L:last_d
  data.table(stay_id=stay_id, day=days, dday=dday)
}
daily_list <- mapply(expand_one, grid$stay_id, grid$max_origin, grid$dday, grid$last_day,
                     SIMPLIFY=FALSE)
daily <- rbindlist(daily_list)

# outcome
daily[, outcome := 0L]
daily[!is.na(dday) & day == dday, outcome := 1L]
daily[, dday := NULL]

# eligible + tx_init 起点日
orig <- unique(trials[, .(stay_id, origin_day, arm)], by=c("stay_id","origin_day"))
daily[, eligible := 0L]
daily[, tx_init := 0L]
daily[orig, on=c("stay_id"="stay_id","day"="origin_day"),
      `:=`(eligible=1L, tx_init=i.arm)]

# 计数验证
log_info(sprintf("  daily grid: %d rows, eligible=%d, stays=%d",
                 nrow(daily), daily[,sum(eligible)], uniqueN(daily$stay_id)))

# ============================================================================
# 7. DAG14 协变量起点日写入 → LOCF → NOCB → 全局值兜底
# ============================================================================
log_info("Step 7: DAG14 协变量填充 (LOCF + NOCB)")

# 全局均值/众数兜底（对仍缺失的）
global_fills <- list()
for (v in DAG14) {
  x <- trials[[v]]
  if (is.factor(x)) global_fills[[v]] <- names(table(x))[1]
  else if (is.numeric(x)) global_fills[[v]] <- median(x, na.rm=TRUE)
  else if (is.integer(x)) global_fills[[v]] <- as.integer(median(x, na.rm=TRUE))
}

# 逐变量写入 + LOCF + NOCB
setorder(daily, stay_id, day)
for (v in DAG14) {
  tcol <- trials[[v]]
  # 初始化 NA 列（类型匹配）
  if (is.factor(tcol)) daily[, (v) := factor(NA, levels=levels(tcol))]
  else if (is.integer(tcol)) daily[, (v) := NA_integer_]
  else if (is.numeric(tcol)) daily[, (v) := NA_real_]
  else daily[, (v) := NA_character_]
  
  # 起点日写入
  daily[orig, on=c("stay_id"="stay_id","day"="origin_day"),
        (v) := get(paste0("i.", v))]
  
  # LOCF 前向
  daily[, (v) := zoo_locf(get(v)), by=stay_id]
  
  # NOCB 后向（stay 开头 NA）
  if (is.factor(tcol)) {
    daily[, (v) := factor(rev(zoo_locf(rev(as.integer(get(v))))),
                          levels=seq_along(levels(tcol)), labels=levels(tcol)),
          by=stay_id]
  } else {
    daily[, (v) := rev(zoo_locf(rev(get(v)))), by=stay_id]
  }
  
  # 全局兜底
  daily[is.na(get(v)), (v) := global_fills[[v]]]
}

log_info(sprintf("  DAG14 无 NA 检查: %d rows 有 NA",
                 nrow(daily[!complete.cases(daily[, ..DAG14])])))

# ============================================================================
# 8. 运行 SEQTaRget ITT（DAG14, imp1, B=50）
# ============================================================================
log_info("Step 8: SEQTaRget ITT × DAG14 (B=50)")

# 断言检查
daily[, eligible_n := sum(eligible), by=stay_id]
stopifnot("eligile must be set" = daily[eligible==1, .N] > 0)

opts_pt <- SEQopts(hazard=TRUE, km.curves=FALSE, bootstrap=FALSE,
                   bootstrap.sample=1, followup.max=TAU,
                   seed=SEED, bootstrap.CI_method="se")
model_pt <- SEQuential(
  as.data.frame(daily),
  id.col="stay_id", time.col="day",
  eligible.col="eligible", treatment.col="tx_init",
  outcome.col="outcome",
  time_varying.cols=DAG14, fixed.cols=list(),
  method="ITT", options=opts_pt, verbose=FALSE)

h0_pt <- hazard_ratio(model_pt)
point_hr <- as.numeric(h0_pt["Hazard ratio"])
log_info(sprintf("  点估计 HR=%.4f", point_hr))

# Bootstrap
opts_b <- SEQopts(hazard=TRUE, km.curves=FALSE, bootstrap=TRUE,
                  bootstrap.nboot=BOOT, bootstrap.sample=1,
                  followup.max=TAU, seed=SEED,
                  bootstrap.CI_method="se")
model_b <- SEQuential(
  as.data.frame(daily),
  id.col="stay_id", time.col="day",
  eligible.col="eligible", treatment.col="tx_init",
  outcome.col="outcome",
  time_varying.cols=DAG14, fixed.cols=list(),
  method="ITT", options=opts_b, verbose=FALSE)

h0_b <- hazard_ratio(model_b)
boot_hr <- as.numeric(h0_b["Hazard ratio"])
boot_lci <- as.numeric(h0_b["LCI"])
boot_uci <- as.numeric(h0_b["UCI"])
log_info(sprintf("  Bootstrap HR=%.4f (%.4f–%.4f)", boot_hr, boot_lci, boot_uci))

# ============================================================================
# 9. 同样跑 DAG10 对比（5 变量子集，无 sepsis3/invasive_vent/...）
# ============================================================================
log_info("Step 9: SEQTaRget ITT × DAG10 (B=50)")

opts_d10 <- SEQopts(hazard=TRUE, km.curves=FALSE, bootstrap=TRUE,
                    bootstrap.nboot=BOOT, bootstrap.sample=1,
                    followup.max=TAU, seed=SEED,
                    bootstrap.CI_method="se")
model_d10 <- SEQuential(
  as.data.frame(daily),
  id.col="stay_id", time.col="day",
  eligible.col="eligible", treatment.col="tx_init",
  outcome.col="outcome",
  time_varying.cols=DAG10, fixed.cols=list(),
  method="ITT", options=opts_d10, verbose=FALSE)
h0_d10 <- hazard_ratio(model_d10)
log_info(sprintf("  DAG10 HR=%.4f (%.4f–%.4f)",
                 as.numeric(h0_d10["Hazard ratio"]),
                 as.numeric(h0_d10["LCI"]), as.numeric(h0_d10["UCI"])))

# ============================================================================
# 10. 保存结果 CSV
# ============================================================================
log_info("Step 10: 保存 CSV")

ev <- function(hr) EValue::evalues.HR(hr, rare=FALSE)
dag14_events <- daily[outcome==1, .(n=uniqueN(stay_id))]$n[1]
dag10_events <- dag14_events  # 同一 daily 表，事件不变

fwrite(data.table(
  analysis = "eICU DAG14 (完整 14v, invasive_vent=intubated, sofa0主合并)",
  framework = "SEQTaRget trial-period ITT",
  dag = "DAG14 (14 vars, tier1 direct + tier2 proxy)",
  n_trial_periods = trials[, .N],
  n_stays = uniqueN(trials$stay_id),
  n_events = dag14_events,
  mortality = round(dag14_events / uniqueN(trials$stay_id) * 100, 1),
  point_hr = point_hr,
  boot_hr = boot_hr, boot_lci = boot_lci, boot_uci = boot_uci,
  boot_se_log = (log(boot_uci) - log(boot_lci)) / (2 * qnorm(0.975)),
  seed = SEED, nboot = BOOT,
  evalue_point = round(ev(boot_hr)$est, 2)
), file.path(OUT_DIR, "eicu_seqtarget_dag14_b50.csv"))

fwrite(data.table(
  source = "eICU-CRD v2.0",
  framework = "SEQTaRget trial-period ITT",
  cov_set = "DAG10 (10 vars, 仅直接匹配)",
  n_trial_periods = trials[, .N],
  n_stays = uniqueN(trials$stay_id),
  n_events = dag10_events,
  mortality = round(dag10_events / uniqueN(trials$stay_id) * 100, 1),
  hr_point = as.numeric(h0_d10["Hazard ratio"]),
  hr_b50 = as.numeric(h0_d10["Hazard ratio"]),
  lci_b50 = as.numeric(h0_d10["LCI"]),
  uci_b50 = as.numeric(h0_d10["UCI"]),
  bootstrap_n = BOOT, tau = TAU, seed = SEED
), file.path(OUT_DIR, "eicu_seqtarget_dag10_b50.csv"))

# ============================================================================
# 11. 最终验证
# ============================================================================
log_info("========== 最终结果 ==========")
log_info(sprintf("eICU DAG14: HR=%.4f (%.4f–%.4f)", boot_hr, boot_lci, boot_uci))
log_info(sprintf("eICU DAG10: HR=%.4f (%.4f–%.4f)",
                 as.numeric(h0_d10["Hazard ratio"]),
                 as.numeric(h0_d10["LCI"]), as.numeric(h0_d10["UCI"])))
log_info(sprintf("events: %d (%.1f%% mortality), trial-periods: %d, stays: %d",
                 dag14_events, dag14_events/uniqueN(trials$stay_id)*100,
                 trials[, .N], uniqueN(trials$stay_id)))
log_info(sprintf("E-value (DAG14 point): %.2f", round(ev(boot_hr)$est, 2)))
log_info("全部完成")
