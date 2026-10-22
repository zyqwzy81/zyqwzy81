# ============================================================================
# v28.4_eicu_dag14_m20.R — eICU 外部验证：MICE m=20 + SEQTaRget B=200 + Rubin
# 输入: sql/result/eicu_dag14_full_stay.parquet (自包含 SQL 提取, ICU入院24h窗口)
# 输出: eicu_seqtarget_dag14_m20_b200_rubin.csv + 20 个 imp 中间 CSV
# 时间轴: 自然日历 day 0..28, stay 连续网格, 起点日 eligible=1 + tx_init=arm
# ============================================================================
suppressPackageStartupMessages({library(data.table);library(arrow);library(SEQTaRget);library(mice);library(EValue)})
BASE   <- "/data/小米云盘/文档/论文/血磷"
OUTDIR <- file.path(BASE, "TTE/v28.4/output/tables/eicu_seqtarget")
IMPDIR <- file.path(BASE, "TTE/v28.4/output/tables/eicu_seqtarget/imps")
for (d in c(OUTDIR, IMPDIR)) dir.create(d, showWarnings=FALSE, recursive=TRUE)
TAU   <- 28L; SEED <- 20260629L; BOOT <- 200L; M <- 20L
DAG14 <- c("age","gender","admission_type","baseline_phosphate","baseline_sofa",
           "charlson","vasoactive_at_t0","invasive_vent_at_t0","sepsis3",
           "lactate_t0","refeeding_risk","corrected_calcium","nutrition_support_type")

zoo_locf <- function(x) {
  if (is.factor(x)) {lev<-levels(x);ii<-zoo_locf(as.integer(x));
    return(factor(ii, levels=seq_along(lev), labels=lev)) }
  nas<-which(is.na(x)); if (length(nas)==0) return(x)
  last<-NA; out<-x
  for (i in seq_along(x)) { if (!is.na(x[i])) last<-x[i]
    else if (!is.na(last)) out[i]<-last }
  out
}
rubin_pool <- function(loghr, se, m=length(loghr)) {
  qbar<-mean(loghr); ubar<-mean(se^2); bv<-if(m>1) var(loghr) else 0
  tv<-ubar+bv*(1+1/m); dfm<-if(m>1&&bv>0)(m-1)*(1+ubar/(bv*(1+1/m)))^2 else Inf
  zt<-if(is.finite(dfm)) qt(0.975,dfm) else qnorm(0.975)
  list(hr=exp(qbar), lci=exp(qbar-zt*sqrt(tv)), uci=exp(qbar+zt*sqrt(tv)),
       p=2*(if(is.finite(dfm)) pt(-abs(qbar/sqrt(tv)),dfm) else pnorm(-abs(qbar/sqrt(tv)))),
       df=dfm, loghr=qbar, se_pool=sqrt(tv))
}

# ---- 1. 加载 DAG14 parquet + 构建 cohort ----
cat("=== Step 1: 加载 DAG14 parquet ===\n")
t <- as.data.table(read_parquet(file.path(BASE,"sql/result/eicu_dag14_full_stay.parquet")))
t <- t[!is.na(itt_group)]  # drop 无治疗分配
cat(sprintf("  stays=%d treatment=%d conservative=%d\n", nrow(t),
    sum(t$itt_group=="treatment"), sum(t$itt_group=="conservative")))

# time_zero 从 full parquet 里已经有 (datetime64 → R POSIXct)
# 但有 2055 NA — 需要从 sequential_trials 补
if (sum(is.na(t$time_zero)) > 0) {
  seq_tb <- as.data.table(read_parquet(file.path(BASE,"sql/result/eicu_sequential_trials.parquet")))
  setnames(seq_tb, "patientunitstayid", "stay_id")
  seq_u <- unique(seq_tb[!is.na(time_zero), .(stay_id, time_zero_seq=time_zero)], by="stay_id")
  t <- merge(t, seq_u, by="stay_id", all.x=TRUE)
  t[, time_zero := fifelse(is.na(time_zero), time_zero_seq, time_zero)]
  t[, time_zero_seq := NULL]
}
t[, arm := as.integer(itt_group=="treatment")]
t[, tz_date := as.Date(time_zero)]
min_date <- min(t$tz_date, na.rm=TRUE)
t[, origin_day := as.integer(tz_date - min_date)]

# death_day
dth <- t[!is.na(mortality_28d) & mortality_28d==1L]
dth[, death_day := pmin(TAU, pmax(0L, origin_day + as.integer(round(time_to_death_days))))]
dth <- dth[, .(death_day = as.integer(names(sort(table(death_day), decreasing=T))[1])), by=stay_id]
t <- merge(t, dth, by="stay_id", all.x=TRUE)
# 保留 origin_day ≤ death_day 的 stay（安全过滤）
t <- t[is.na(death_day) | origin_day <= death_day]
t[, death_28day_num := as.integer(mortality_28d==1)]
cat(sprintf("  after filter: stays=%d deaths=%d origin_day range=%d..%d\n",
    nrow(t), sum(t$death_28day_num,na.rm=T), min(t$origin_day), max(t$origin_day)))

# ---- 2. 构建 R 数据 (所有 DAG14 变量) ----
cat("\n=== Step 2: 构建 DAG14 变量 ===\n")
t[, age := as.numeric(age)]
t[, gender := factor(fifelse(sex=="M","M","F"), levels=c("F","M"))]
t[, admission_type := factor(admission_type)]
t[, baseline_phosphate := as.numeric(baseline_phosphate_mmol_l)]
t[, charlson := as.numeric(charlson)]
t[, baseline_sofa := as.numeric(baseline_sofa)]
t[, vasoactive_at_t0 := as.integer(vasoactive_at_t0)]
t[, invasive_vent_at_t0 := as.integer(invasive_vent_at_t0)]
t[, sepsis3 := as.integer(sepsis3)]
t[, lactate_t0 := as.numeric(lactate_max)]
t[, corrected_calcium := as.numeric(albumin_min)]  # 代理 (eICU 无 real Ca)
t[, refeeding_risk := as.integer(is.na(sodium))]   # 占位 (需 K/Mg)
t[, nutrition_support_type := fifelse(
    cumulative_iv_phosphate_mmol_48h>0 & has_oral_ng_phosphate_48h>0, "mixed",
    fifelse(cumulative_iv_phosphate_mmol_48h>0, "parenteral",
      fifelse(has_oral_ng_phosphate_48h>0, "enteral", "none")))]
t[, nutrition_support_type := factor(nutrition_support_type,
    levels=c("none","enteral","parenteral","mixed"))]

# ---- 3. MICE m=20 ----
cat("\n=== Step 3: MICE m=20 ===\n")
mice_in <- data.table(
  lactate_t0 = t$lactate_t0,
  corrected_calcium = t$corrected_calcium,
  baseline_phosphate = t$baseline_phosphate,
  baseline_sofa = t$baseline_sofa,
  age = t$age,
  arm = t$arm,
  death_28day_num = t$death_28day_num,
  sepsis3 = t$sepsis3,
  vasoactive_at_t0 = t$vasoactive_at_t0,
  invasive_vent_at_t0 = t$invasive_vent_at_t0,
  charlson = t$charlson
)
method_vec <- c(lactate_t0="pmm", corrected_calcium="pmm", baseline_phosphate="pmm",
                baseline_sofa="", age="", arm="", death_28day_num="",
                sepsis3="", vasoactive_at_t0="", invasive_vent_at_t0="", charlson="")
t0 <- Sys.time()
mm <- mice(mice_in, m=M, maxit=15, seed=SEED, method=method_vec, printFlag=FALSE)
cat(sprintf("  MICE done %.1f min\n", as.numeric(difftime(Sys.time(),t0,units="min"))))
for (i in seq_len(M)) write_parquet(complete(mm, i), file.path(IMPDIR, sprintf("imp%02d.parquet",i)))
cat(sprintf("  saved %d imp parquets → %s\n", M, IMPDIR))

# ---- 4. 20 imp × SEQTaRget B=200 (分批 + 每 imp 单独 CSV + GC) ----
cat("\n=== Step 4: 20 imp × SEQTaRget B=200 ===\n")
extract_log <- function(h0) {
  z<-qnorm(0.975); hr<-as.numeric(h0["Hazard ratio"])
  list(hr=hr, loghr=log(hr), se=(log(as.numeric(h0["UCI"]))-log(as.numeric(h0["LCI"])))/(2*z))
}
hrs_all <- data.table(imp=integer(), hr=numeric(), loghr=numeric(), se=numeric())

for (i in seq_len(M)) {
  t0_imp <- Sys.time()
  imp <- as.data.table(read_parquet(file.path(IMPDIR, sprintf("imp%02d.parquet",i))))

  # 覆盖 MICE 插补值
  dt <- copy(t)
  dt[, lactate_t0 := imp$lactate_t0]
  dt[, corrected_calcium := imp$corrected_calcium]
  dt[, baseline_phosphate := imp$baseline_phosphate]

  # daily-long 展开
  grid <- dt[!is.na(origin_day), .(origin_day = origin_day[1],
    dday = fifelse(!is.na(death_day), death_day[1], NA_integer_)), by=stay_id]
  grid[, last_day := origin_day + TAU, by=stay_id]
  daily_list <- mapply(function(s,o,d,l)
    data.table(stay_id=s, day=seq.int(o, min(l, o+TAU-1L), by=1L), dday=d),
    grid$stay_id, grid$origin_day, grid$dday, grid$last_day, SIMPLIFY=FALSE)
  daily <- rbindlist(daily_list, fill=TRUE)
  daily[, outcome := 0L]
  daily[!is.na(dday) & day==dday, outcome := 1L]

  # 写协变量 (起点日) → LOCF → NOCB → 全局兜底
  orig <- unique(dt[, c("stay_id","origin_day","arm",DAG14), with=FALSE],
                 by=c("stay_id","origin_day"))
  setnames(orig, "origin_day", "day")
  daily <- merge(daily, orig, by=c("stay_id","day"), all.x=TRUE)
  daily[, eligible := as.integer(!is.na(arm))]
  daily[, tx_init := fifelse(eligible==1L, arm, 0L)]
  daily[, arm := NULL]

  setorder(daily, stay_id, day)
  for (v in DAG14) {
    tcol <- orig[[v]]; if (is.null(tcol)) next
    if (is.factor(tcol)) {
      daily[, (v) := factor(zoo_locf(as.integer(get(v))),
            levels=seq_along(levels(tcol)), labels=levels(tcol)), by=stay_id]
      daily[, (v) := factor(rev(zoo_locf(rev(as.integer(get(v))))),
            levels=seq_along(levels(tcol)), labels=levels(tcol)), by=stay_id]
    } else {
      daily[, (v) := zoo_locf(get(v)), by=stay_id]
      daily[, (v) := rev(zoo_locf(rev(get(v)))), by=stay_id]
    }
    gval <- if (is.factor(tcol)) names(table(tcol))[1] else median(tcol, na.rm=TRUE)
    daily[is.na(get(v)), (v) := gval]
  }

  # 检查 NA
  na_remaining <- sum(is.na(daily[, eligible]) | is.na(daily[, tx_init]))
  if (na_remaining > 0) daily[is.na(eligible), eligible := 0L]
  daily <- as.data.frame(daily)

  # SEQTaRget
  opts <- SEQopts(hazard=TRUE, km.curves=FALSE, bootstrap=TRUE,
    bootstrap.nboot=BOOT, bootstrap.sample=1,
    followup.max=TAU, seed=SEED+i*7L, bootstrap.CI_method="se")
  mod <- SEQuential(daily,
    id.col="stay_id", time.col="day",
    eligible.col="eligible", treatment.col="tx_init",
    outcome.col="outcome",
    time_varying.cols=DAG14, fixed.cols=list(),
    method="ITT", options=opts, verbose=FALSE)
  h <- extract_log(hazard_ratio(mod))
  hrs_all <- rbindlist(list(hrs_all, data.table(imp=i, hr=h$hr, loghr=h$loghr, se=h$se)))
  cat(sprintf("  imp%02d: HR=%.4f (%.1f s)\n", i, h$hr,
      as.numeric(difftime(Sys.time(),t0_imp,units="sec"))))

  # 写中间 CSV + GC
  fwrite(data.table(imp=i, hr=h$hr, loghr=h$loghr, se=h$se),
    file.path(IMPDIR, sprintf("imp%02d_result.csv",i)))
  rm(dt, grid, daily_list, daily, orig, mod, imp); gc(verbose=FALSE)
}

# ---- 5. Rubin 合并 + 输出 ----
cat("\n=== Step 5: Rubin 合并 ===\n")
pooled <- rubin_pool(hrs_all$loghr, hrs_all$se, m=M)
ev_pt <- EValue::evalues.HR(pooled$hr, rare=FALSE)
events <- sum(t$death_28day_num, na.rm=TRUE)

cat(sprintf("HR = %.4f (95%% CI %.4f–%.4f)\np  = %.4f\ndf = %.1f\n",
    pooled$hr, pooled$lci, pooled$uci, pooled$p, pooled$df))
cat(sprintf("events=%d stays=%d (%.1f%% mortality)\n",
    events, uniqueN(t$stay_id), 100*events/uniqueN(t$stay_id)))
cat(sprintf("E(point)=%.2f, E(bound)=%.2f\n", ev_pt$est, ev_pt$lo))
cat(sprintf("逐 imp HR: range=%.4f..%.4f, sd=%.4f\n",
    min(hrs_all$hr), max(hrs_all$hr), sd(hrs_all$hr)))

out_csv <- file.path(OUTDIR, sprintf("eicu_seqtarget_dag14_m%d_b%d_rubin.csv", M, BOOT))
fwrite(data.table(
  analysis = sprintf("eICU DAG14 m=%d B=%d Rubin (自包含SQL提取)", M, BOOT),
  n_stays = uniqueN(t$stay_id), n_events = events,
  mortality = round(events/uniqueN(t$stay_id)*100,1),
  n_treat = sum(t$itt_group=="treatment"),
  hr_rubin = pooled$hr, lci_rubin = pooled$lci, uci_rubin = pooled$uci,
  p_value = pooled$p, rubin_df = round(pooled$df, 1),
  evalue_point = round(ev_pt$est, 2), evalue_bound = round(ev_pt$lo, 2),
  se_imps = paste(round(hrs_all$hr,4), collapse=";")),
  out_csv)
cat(sprintf("\n✅ 保存: %s\n", out_csv))
