# ============================================================================
# v28.4_eicu_seqtarget_m20_v2.R
# 完整修复版：MICE 仅 3 变量, SEQTaRget 每 imp 后 gc(), Rubin 合并
# ============================================================================
suppressPackageStartupMessages({library(data.table);library(arrow);library(SEQTaRget);library(mice);library(EValue)})
BASE <- "/data/小米云盘/文档/论文/血磷"
IMPDIR <- file.path(BASE, "sql/result/eicu_mice_imputed")
OUT_DIR <- file.path(BASE, "TTE/v28.4/output/tables/eicu_seqtarget")
dir.create(OUT_DIR, showWarnings=F, recursive=T)
TAU <- 28L; SEED <- 20260629L; BOOT <- 50L; M <- 20L
DAG14 <- c("age","gender","admission_type","baseline_phosphate","baseline_sofa",
           "charlson","egfr","vasoactive_at_t0","invasive_vent_at_t0","sepsis3",
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
       df=dfm)
}

# ---- 1. 加载 + 合并 trial-level ----
bc <- as.data.table(read_parquet(file.path(BASE,"sql/result/eicu_baseline_covariates.parquet")))
sofa0 <- as.data.table(read_parquet(file.path(BASE,"sql/result/eicu_sofa_t0.parquet")))
exp <- as.data.table(read_parquet(file.path(BASE,"sql/result/eicu_exposure.parquet")))
tr <- as.data.table(read_parquet(file.path(BASE,"sql/result/eicu_sequential_trials.parquet")))
outc <- as.data.table(read_parquet(file.path(BASE,"sql/result/eicu_outcomes.parquet")))
soi <- as.data.table(read_parquet(file.path(BASE,"sql/result/eicu_suspicion_of_infection.parquet")))

t <- sofa0[, .(stay_id=patientunitstayid, tzo=time_zero_offset, sofa_t0, inv=intubated)]
t <- merge(t, bc[, .(stay_id=patientunitstayid, tzo=time_zero_offset,
    age, sex, admission_type, baseline_phosphate_mmol_l, baseline_sofa_t0,
    charlson_comorbidity_index, egfr, vasoactive_at_t0, lactate_max,
    hypokalemia, hypomagnesemia, albumin_min)], by=c("stay_id","tzo"), all.x=T)
t <- merge(t, exp[, .(stay_id=patientunitstayid, tzo=time_zero_offset,
    has_oral=has_oral_ng_phosphate_48h, cum_iv=cumulative_iv_phosphate_mmol_48h)],
    by=c("stay_id","tzo"), all.x=T)
t <- merge(t, tr[, .(stay_id=patientunitstayid, tzo=time_zero_offset,
    itt_group, time_zero)], by=c("stay_id","tzo"), all.x=T)
t <- unique(t[!is.na(itt_group)], by=c("stay_id","tzo"))
t[, arm := as.integer(itt_group=="treatment")]

outc_s <- unique(outc, by="patientunitstayid")
setnames(outc_s, "patientunitstayid", "stay_id")
t <- merge(t, outc_s[, .(stay_id, mortality_28d, time_to_death_days)],
           by="stay_id", all.x=T)
t[, death_28day_num := as.integer(mortality_28d==1)]
t[, origin_day := as.integer(tzo/24)]
t[, min_date := min(as.Date(time_zero)), by=stay_id]
t[, origin_day := as.integer(as.Date(time_zero) - min_date)]
dth <- t[death_28day_num==1L,
         .(cand=pmin(TAU,pmax(0L,origin_day+as.integer(round(time_to_death_days))))),
         by=stay_id]
dth[, death_day := as.integer(names(sort(table(cand), decreasing=T))[1]), by=stay_id]
dth <- dth[!duplicated(stay_id)]
t <- merge(t, dth[, .(stay_id, death_day)], by="stay_id", all.x=T)
t <- t[is.na(death_day) | origin_day <= death_day]

soi_stay <- soi[, .(sus = max(as.integer(suspected_infection=="true"|suspected_infection==1),na.rm=T)),
                by=patientunitstayid]
setnames(soi_stay, "patientunitstayid", "stay_id")
t <- merge(t, soi_stay, by="stay_id", all.x=T)
t[, soi_flag := as.integer(sus==1)]
t[, baseline_sofa := fifelse(!is.na(sofa_t0) & sofa_t0>=0, sofa_t0, baseline_sofa_t0)]
cat(sprintf("trials=%d stays=%d deaths=%d\n", nrow(t), uniqueN(t$stay_id),
            sum(t$death_28day_num,na.rm=T)))

# ---- 2. MICE: 仅 3 变量需插补 ----
mice_in <- data.table(
  lactate_max = t$lactate_max,
  albumin_min = t$albumin_min,
  egfr = t$egfr,
  # 辅助变量（不插补，仅预测）
  baseline_sofa = t$baseline_sofa,
  age = t$age,
  arm = t$arm,
  death_28day_num = t$death_28day_num,
  time_to_death_days = t$time_to_death_days,
  inv = as.integer(t$inv),
  vasoactive_at_t0 = as.integer(t$vasoactive_at_t0),
  hypokalemia = as.integer(t$hypokalemia),
  hypomagnesemia = as.integer(t$hypomagnesemia),
  charlson_comorbidity_index = t$charlson_comorbidity_index
)
# 显式 method: 仅 3 变量 pmm，其余 "" 完全不碰
method_vec <- c(lactate_max="pmm", albumin_min="pmm", egfr="pmm",
                baseline_sofa="", age="", arm="", death_28day_num="",
                time_to_death_days="", inv="", vasoactive_at_t0="",
                hypokalemia="", hypomagnesemia="", charlson_comorbidity_index="")
cat("MICE method:", paste(names(method_vec), method_vec, sep="=", collapse=","), "\n")

t0 <- Sys.time()
mm <- mice(mice_in, m=M, maxit=15, seed=SEED, method=method_vec, printFlag=FALSE)
cat(sprintf("MICE done %.1f min\n", as.numeric(difftime(Sys.time(),t0,units="min"))))

# 保存 20 个插补集（纯 MICE 结果，后面追加固定列）
for (i in seq_len(M)) write_parquet(complete(mm, i),
    file.path(IMPDIR, sprintf("eicu_mice_imp%02d.parquet", i)))

# ---- 3. 20 imp × SEQTaRget B=50 ----
extract_log <- function(h0) {
  z <- qnorm(0.975)
  list(hr=as.numeric(h0["Hazard ratio"]),
       loghr=log(as.numeric(h0["Hazard ratio"])),
       se=(log(as.numeric(h0["UCI"])) - log(as.numeric(h0["LCI"])))/(2*z))
}
loghrs <- numeric(M); ses <- numeric(M); hrs_out <- numeric(M)
for (i in seq_len(M)) {
  imp <- as.data.table(read_parquet(file.path(IMPDIR, sprintf("eicu_mice_imp%02d.parquet", i))))

  dt <- copy(t)
  dt[, lactate_t0 := imp$lactate_max]
  dt[, corrected_calcium := imp$albumin_min]
  dt[, egfr := imp$egfr]
  dt[, baseline_sofa := imp$baseline_sofa]
  # DAG14 派生
  dt[, age := imp$age]
  dt[, gender := factor(fifelse(dt$sex=="M", "M", "F"), levels=c("F","M"))]
  dt[, admission_type := factor(admission_type)]
  dt[, baseline_phosphate := baseline_phosphate_mmol_l]
  dt[, charlson := charlson_comorbidity_index]
  dt[, vasoactive_at_t0 := as.integer(vasoactive_at_t0)]
  dt[, invasive_vent_at_t0 := as.integer(inv)]
  dt[, refeeding_risk := as.integer(hypokalemia==1 | hypomagnesemia==1)]
  dt[, sepsis3 := as.integer(soi_flag==1 & baseline_sofa>=2)]
  dt[, nutrition_support_type := fifelse(
      cum_iv>0 & has_oral>0, "mixed",
      fifelse(cum_iv>0, "parenteral",
        fifelse(has_oral>0, "enteral", "none")))]
  dt[, nutrition_support_type := factor(nutrition_support_type,
      levels=c("none","enteral","parenteral","mixed"))]

  # daily-long
  grid <- dt[, .(max_o=max(origin_day),
                 dday=fifelse(any(death_28day_num==1L), death_day[1], NA_integer_)),
              by=stay_id]
  grid[, last_day := max_o + TAU]
  daily_list <- mapply(function(s,m,d,l)
    data.table(stay_id=s, day=0L:l, dday=d),
    grid$stay_id, grid$max_o, grid$dday, grid$last_day, SIMPLIFY=FALSE)
  daily <- rbindlist(daily_list)
  daily[, outcome := 0L]
  daily[!is.na(dday) & day==dday, outcome := 1L]
  daily[, dday := NULL]

  orig <- unique(dt[, c("stay_id","origin_day","arm",DAG14), with=FALSE],
                 by=c("stay_id","origin_day"))
  setnames(orig, "origin_day", "day")
  daily <- merge(daily, orig, by=c("stay_id","day"), all.x=TRUE)
  daily[!is.na(arm), eligible := 1L]
  daily[eligible==1L, tx_init := arm]
  daily[, arm := NULL]
  daily[is.na(eligible), eligible := 0L]
  daily[is.na(tx_init), tx_init := 0L]

  setorder(daily, stay_id, day)
  for (v in DAG14) {
    tcol <- dt[[v]]
    if (!is.factor(tcol)) {
      daily[, (v) := zoo_locf(get(v)), by=stay_id]
      daily[, (v) := rev(zoo_locf(rev(get(v)))), by=stay_id]
    } else {
      daily[, (v) := factor(rev(zoo_locf(rev(as.integer(get(v))))),
                            levels=seq_along(levels(tcol)), labels=levels(tcol)), by=stay_id]
    }
    gval <- if (is.factor(tcol)) names(table(tcol))[1] else median(tcol, na.rm=T)
    daily[is.na(get(v)), (v) := gval]
  }
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
  hrs_out[i] <- h$hr; loghrs[i] <- h$loghr; ses[i] <- h$se
  cat(sprintf("imp%02d: HR=%.4f\n", i, h$hr))

  rm(dt, grid, daily_list, daily, orig, mod); gc(verbose=FALSE)
}

# ---- 4. Rubin 合并 ----
pooled <- rubin_pool(loghrs, ses, m=M)
cat(sprintf("\n=== Rubin 合并 ===\n"))
cat(sprintf("HR = %.4f (95%% CI %.4f–%.4f)\np  = %.4f\ndf = %.1f\n",
            pooled$hr, pooled$lci, pooled$uci, pooled$p, pooled$df))
cat(sprintf("逐 imp HR: min=%.4f, max=%.4f, range=%.4f, sd=%.4f\n",
            min(hrs_out), max(hrs_out), diff(range(hrs_out)), sd(hrs_out)))
ev_pt <- EValue::evalues.HR(pooled$hr, rare=FALSE)
cat(sprintf("E(point)=%.2f, E(bound)=%.2f\n", ev_pt$est, ev_pt$lo))

events <- uniqueN(t[death_28day_num==1, stay_id])
fwrite(data.table(
  analysis = sprintf("eICU DAG14 m=%d B=%d Rubin", M, BOOT),
  n_trial_periods = nrow(t), n_stays = uniqueN(t$stay_id),
  n_events = events, mortality = round(events/uniqueN(t$stay_id)*100,1),
  hr_rubin = pooled$hr, lci_rubin = pooled$lci, uci_rubin = pooled$uci,
  p_value = pooled$p, rubin_df = round(pooled$df, 1),
  evalue_point = round(ev_pt$est, 2), evalue_bound = round(ev_pt$lo, 2)
), file.path(OUT_DIR, sprintf("eicu_seqtarget_dag14_m%d_b%d_rubin.csv", M, BOOT)))
cat(sprintf("saved → %s\n", file.path(OUT_DIR, sprintf("eicu_seqtarget_dag14_m%d_b%d_rubin.csv", M, BOOT))))
