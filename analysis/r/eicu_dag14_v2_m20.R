# ============================================================================
# v28.4_eicu_dag14_v2_m20.R — eICU DAG14 v2 (refeeding_risk 真正修复版)
# 输入: eicu_dag14_full_stay_v2.parquet (SQL 已派生 refeeding_risk=低钾<3.5 OR 低镁<0.8)
# 修复: refeeding_risk 从占位符 is.na(sodium) 改为真正的临床变量 K/Mg
# ============================================================================
suppressPackageStartupMessages({library(data.table);library(arrow);library(SEQTaRget);library(mice);library(EValue)})
BASE   <- "/data/小米云盘/文档/论文/血磷"
OUTDIR <- file.path(BASE, "TTE/v28.4/output/tables/eicu_seqtarget")
IMPDIR <- file.path(BASE, "TTE/v28.4/output/tables/eicu_seqtarget/imps_v2")
for (d in c(OUTDIR, IMPDIR)) dir.create(d, showWarnings=FALSE, recursive=TRUE)
TAU   <- 28L; SEED <- 20260629L; BOOT <- 200L; M <- 20L

# DAG14 变量集
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

# ---- 1. 加载 v2 DAG14 parquet ----
cat("=== Step 1: 加载 v2 DAG14 parquet (含 refeeding_risk SQL 派生) ===\n")
t <- as.data.table(read_parquet(file.path(BASE,"sql/result/eicu_dag14_full_stay_v2.parquet")))
t <- t[!is.na(itt_group)]
cat(sprintf("  stays=%d treatment=%d conservative=%d\n", nrow(t),
    sum(t$itt_group=="treatment"), sum(t$itt_group=="conservative")))

# 检查 v2 DAG14 关键变量覆盖
cat(sprintf("\n  refeeding_risk: non-NA=%d, =1: %d (%.1f%%)\n",
    sum(!is.na(t$refeeding_risk)), sum(t$refeeding_risk==1,na.rm=T),
    100*mean(t$refeeding_risk==1,na.rm=T)))
cat(sprintf("  potassium_t0: non-NA=%d (%.1f%%)\n",
    sum(!is.na(t$potassium_t0)), 100*mean(!is.na(t$potassium_t0))))
cat(sprintf("  magnesium_t0: non-NA=%d (%.1f%%)\n",
    sum(!is.na(t$magnesium_t0)), 100*mean(!is.na(t$magnesium_t0))))

# time_zero
t[, arm := as.integer(itt_group=="treatment")]
t[, tz_date := as.Date(time_zero)]
min_date <- min(t$tz_date, na.rm=TRUE)
t[, origin_day := as.integer(tz_date - min_date)]

# death_day —— v28.4.1（2026-10-14）绝对日历日共识，与 MIMIC 主管线同源：
# 死亡日历日 as.Date(time_zero + time_to_death_days*86400) 减 time_zero 日历日。
# eICU 每 stay 单起点，共识众数仅为防御性保留。
# 旧规则 origin_day + round(time_to_death_days) 存在非整点舍入误差。
dth <- t[!is.na(mortality_28d) & mortality_28d==1L]
dth[, death_abs  := as.integer(as.Date(time_zero + time_to_death_days*86400))]
dth[, origin_abs := as.integer(as.Date(time_zero))]
dth[, death_day_new := origin_day + (death_abs - origin_abs)]
dth[, death_day_old := origin_day + as.integer(round(time_to_death_days))]
n_shift <- sum(dth$death_day_new != dth$death_day_old, na.rm=TRUE)
cat(sprintf("  death-day rule: new-vs-old shifted stays=%d/%d\n",
            n_shift, nrow(dth)))
dth[, death_day := pmin(TAU, pmax(0L, death_day_new))]
dth <- dth[, .(death_day = as.integer(names(sort(table(death_day), decreasing=T))[1])), by=stay_id]
t <- merge(t, dth, by="stay_id", all.x=TRUE)
t <- t[is.na(death_day) | origin_day <= death_day]
t[, death_28day_num := as.integer(mortality_28d==1)]
cat(sprintf("\n  after filter: stays=%d deaths=%d origin_day=%d..%d\n",
    nrow(t), sum(t$death_28day_num,na.rm=T), min(t$origin_day,na.rm=TRUE), max(t$origin_day,na.rm=TRUE)))

# ---- 2. 构建 DAG14 变量 (v2 版) ----
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
t[, corrected_calcium := as.numeric(calcium_min)]  # SQL 已派生 (94.3% 覆盖, vs albumin_min 59%)

# ✅ ref eeding_risk: SQL 已派生 (K<3.5 OR Mg<0.8) — 直接使用
t[, refeeding_risk := as.integer(refeeding_risk)]
cat(sprintf("  refeeding_risk=1: %d (%.1f%%)\n", sum(t$refeeding_risk==1,na.rm=T),
    100*mean(t$refeeding_risk==1,na.rm=T)))

# nutrition_support_type (cumulative_iv 全 0 → 全是 none)
t[, nutrition_support_type := "none"]
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

# ---- 4. 20 imp × SEQTaRget B=200 ----
cat("\n=== Step 4: 20 imp × SEQTaRget B=200 ===\n")
extract_log <- function(h0) {
  z<-qnorm(0.975); hr<-as.numeric(h0["Hazard ratio"])
  list(hr=hr, loghr=log(hr), se=(log(as.numeric(h0["UCI"]))-log(as.numeric(h0["LCI"])))/(2*z))
}
hrs_all <- data.table(imp=integer(), hr=numeric(), loghr=numeric(), se=numeric())

for (i in seq_len(M)) {
  t0_imp <- Sys.time()
  imp <- as.data.table(read_parquet(file.path(IMPDIR, sprintf("imp%02d.parquet",i))))

  dt <- copy(t)
  dt[, lactate_t0 := imp$lactate_t0]
  dt[, corrected_calcium := imp$corrected_calcium]
  dt[, baseline_phosphate := imp$baseline_phosphate]

  grid <- dt[!is.na(origin_day), .(origin_day = origin_day[1],
    dday = fifelse(!is.na(death_day), death_day[1], NA_integer_)), by=stay_id]
  grid[, last_day := origin_day + TAU, by=stay_id]
  daily_list <- mapply(function(s,o,d,l)
    data.table(stay_id=s, day=seq.int(o, min(l, o+TAU-1L), by=1L), dday=d),
    grid$stay_id, grid$origin_day, grid$dday, grid$last_day, SIMPLIFY=FALSE)
  daily <- rbindlist(daily_list, fill=TRUE)
  daily[, outcome := 0L]
  daily[!is.na(dday) & day==dday, outcome := 1L]

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

  na_remaining <- sum(is.na(daily[, eligible]) | is.na(daily[, tx_init]))
  if (na_remaining > 0) daily[is.na(eligible), eligible := 0L]
  daily <- as.data.frame(daily)

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

  fwrite(data.table(imp=i, hr=h$hr, loghr=h$loghr, se=h$se),
    file.path(IMPDIR, sprintf("imp%02d_result.csv",i)))
  rm(dt, grid, daily_list, daily, orig, mod, imp); gc(verbose=FALSE)
}

# ---- 5. Rubin 合并 ----
cat("\n=== Step 5: Rubin 合并 ===\n")
pooled <- rubin_pool(hrs_all$loghr, hrs_all$se, m=M)
ev <- evalues.HR(pooled$hr, pooled$lci, pooled$uci, rare=FALSE)
ev_pt <- if (is.data.frame(ev)) ev$est[1] else ev[[1]]
ev_lo <- if (is.data.frame(ev)) ev$lo[1] else ev[[2]]
events <- sum(t$death_28day_num, na.rm=TRUE)

cat(sprintf("HR = %.4f (95%% CI %.4f–%.4f)\np  = %.4f\ndf = %.1f\n",
    pooled$hr, pooled$lci, pooled$uci, pooled$p, pooled$df))
cat(sprintf("events=%d stays=%d (%.1f%% mortality)\n",
    events, uniqueN(t$stay_id), 100*events/uniqueN(t$stay_id)))
cat(sprintf("E(point)=%.2f, E(bound)=%.2f\n", ev_pt, ev_lo))
cat(sprintf("逐 imp HR: range=%.4f..%.4f, sd=%.4f\n",
    min(hrs_all$hr), max(hrs_all$hr), sd(hrs_all$hr)))

out_csv <- file.path(OUTDIR, sprintf("eicu_seqtarget_dag14_v2_m%d_b%d_rubin.csv", M, BOOT))
fwrite(data.table(
  analysis = sprintf("eICU DAG14-v2 m=%d B=%d (refeeding_risk=K<3.5|Mg<0.8, SQL 派生)", M, BOOT),
  n_stays = uniqueN(t$stay_id), n_events = events,
  mortality = round(events/uniqueN(t$stay_id)*100,1),
  n_treat = sum(t$itt_group=="treatment"),
  hr_rubin = pooled$hr, lci_rubin = pooled$lci, uci_rubin = pooled$uci,
  p_value = pooled$p, rubin_df = round(pooled$df, 1),
  evalue_point = round(ev_pt, 2), evalue_bound = round(ev_lo, 2),
  hr_imps = paste(round(hrs_all$hr,4), collapse="; ")),
  out_csv)
cat(sprintf("\n✅ 保存: %s\n", out_csv))
