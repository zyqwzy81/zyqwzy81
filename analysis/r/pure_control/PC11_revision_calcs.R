# =====================================================================
# PC11: Post-review revision calculations (independent review 2026-10-11)
# Covers: M1 measurement-model c-statistic
#         M2 weighted time-averaged phosphate concentration + bootstrap CI
#         M5 safety event RR/RD with B=200 bootstrap CI (C1 + C3)
#         M6 pre-T0 slope difference CI + GEE time*group interaction p
#         M7a RCS nonlinearity test
#         TR-8.3 first-24h vs cumulative-48h dose metric comparison
#         M8 measured-vs-unmeasured baseline SMD (pre/post IPCW)
#         M3 A_IV vs primary-analysis treatment concordance table
# POST-HOC / hypothesis-generating addenda. Does NOT touch main analysis.
# DB password via env var: MIMIC_PG_PASSWORD (never hard-code).
# =====================================================================
suppressMessages({library(DBI); library(RPostgres); library(dplyr); library(tidyr)
  library(geepack); library(rms); library(readr)}
)
set.seed(20260629L)
B <- 200L
OUT <- "/data/小米云盘/文档/论文/血磷/TTE/v28.4/output/pure_control"
MMPASS <- Sys.getenv("MIMIC_PG_PASSWORD")
con <- dbConnect(Postgres(), host="localhost", port=5432, dbname="mimiciv",
                 user="postgres", password=MMPASS)
win <- dbGetQuery(con, "SELECT * FROM public.pure_control_phos_windows")
co  <- dbGetQuery(con, "SELECT * FROM public.pure_control_cohort")
lng  <- dbGetQuery(con, "SELECT stay_id, hph, phos FROM public.pure_control_phos_long")
saf  <- dbGetQuery(con, "SELECT * FROM public.pure_control_safety")
pre  <- dbGetQuery(con, "SELECT stay_id,win,actual_hph,phos FROM public.pure_control_prephos")
dose <- dbGetQuery(con, "SELECT stay_id,dose24,dose48,n_doses FROM public.pure_control_dose")
conc <- dbGetQuery(con, "
  WITH scs AS (SELECT stay_id, MIN(itt_group_48h) AS itt_group_48h
               FROM public.study_cohort GROUP BY stay_id)
  SELECT c.exposure_w48, scs.itt_group_48h, COUNT(*) AS n
  FROM public.pure_control_cohort c JOIN scs USING (stay_id)
  GROUP BY 1,2 ORDER BY 1,2")
dbDisconnect(con)

co <- co %>% mutate(sexM=as.integer(gender=="M"),
  sepsis3=as.integer(sepsis3), vent_dim3=as.integer(vent_dim3),
  vasoactive_at_t0=as.integer(vasoactive_at_t0), refeeding_risk=as.integer(refeeding_risk),
  adm_elective=as.integer(grepl("ELECTIVE|SURGICAL|AMBULATORY",admission_type)),
  nutr_any=as.integer(nutrition_support_type!="none"))
DAG14 <- c("age","sexM","baseline_sofa","sepsis3","charlson","egfr","vent_dim3",
           "vasoactive_at_t0","admission_type","baseline_phosphate","lactate_t0",
           "refeeding_risk","corrected_calcium_mmol_min","nutrition_support_type")
SMDV <- c("age","sexM","baseline_sofa","sepsis3","charlson","egfr","vent_dim3",
          "vasoactive_at_t0","baseline_phosphate","lactate_t0","refeeding_risk",
          "corrected_calcium_mmol_min","adm_elective","nutr_any")
wcent <- c(w6=6,w12=12,w24=24,w48=48,w72=72)
MMCOV <- c("age","baseline_sofa","sepsis3","egfr","vent_dim3",
           "vasoactive_at_t0","baseline_phosphate","charlson")

# stabilized IPTW weights for a complete-case two-arm subset
psw <- function(d, Acol="A"){
  fit <- suppressWarnings(glm(as.formula(paste(Acol,"~",paste(DAG14,collapse="+"))),
                              data=d, family=binomial))
  e <- predict(fit, type="response"); pA <- mean(d[[Acol]])
  list(w=ifelse(d[[Acol]]==1, pA/e, (1-pA)/(1-e)), e=e, pA=pA)
}
q025 <- function(x,...) quantile(x, .025, na.rm=TRUE, ...)
q975 <- function(x,...) quantile(x, .975, na.rm=TRUE, ...)

# =====================================================================
# M1 + M8: C1 balanced measurement panel (exact PC02 construction)
# =====================================================================
wlong <- win %>% mutate(time=wcent[win]) %>%
  select(stay_id, exposure_w48, time, actual_hph, phos) %>%
  filter(exposure_w48 %in% c("A_IV","D_none")) %>%
  inner_join(co %>% select(stay_id, all_of(DAG14), all_of(SMDV)), by="stay_id") %>%
  filter(complete.cases(across(all_of(DAG14))))
d0 <- co %>% filter(exposure_w48 %in% c("A_IV","D_none")) %>%
  select(stay_id, exposure_w48, baseline_phosphate, all_of(DAG14), all_of(setdiff(SMDV,DAG14))) %>%
  filter(complete.cases(across(all_of(DAG14))))
grid <- expand.grid(stay_id=unique(d0$stay_id), time=c(6,12,24,48,72),
                    stringsAsFactors=FALSE) %>%
  left_join(d0, by="stay_id") %>%
  left_join(wlong %>% select(stay_id,time,phos), by=c("stay_id","time")) %>%
  mutate(measured=as.integer(!is.na(phos)),
         A=as.integer(exposure_w48=="A_IV"), timef=factor(time))
ed <- grid %>% distinct(stay_id,.keep_all=TRUE)
ps1 <- glm(as.formula(paste("A~",paste(DAG14,collapse="+"))), data=ed, family=binomial)
ed$e <- predict(ps1, type="response"); pA <- mean(ed$A)
ed$iptw <- ifelse(ed$A==1, pA/ed$e, (1-pA)/(1-ed$e))
mm  <- suppressWarnings(glm(as.formula(paste(
        "measured ~ A*timef +",paste(MMCOV,collapse="+"))),
        data=grid, family=binomial))
mm0 <- glm(measured ~ A*timef, data=grid, family=binomial)
grid$pm  <- predict(mm, type="response")
grid$pm0 <- predict(mm0, type="response")
grid <- grid %>% left_join(ed %>% select(stay_id,iptw), by="stay_id")

auc_mw <- function(y, p){
  ok <- !is.na(p); y<-y[ok]; p<-p[ok]
  n1<-sum(y==1); n0<-sum(y==0)
  if(n1==0||n0==0) return(NA_real_)
  r<-rank(p, ties.method="average")
  as.numeric((sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0))
}
cstat_pooled <- auc_mw(grid$measured, grid$pm)
ever <- grid %>% group_by(stay_id) %>% filter(max(measured)==1) %>% ungroup()
cstat_ever  <- auc_mw(ever$measured, ever$pm)
st <- grid %>% group_by(stay_id) %>% summarise(y=max(measured),
  p_any=1-prod(1-pm), p_mean=mean(pm), .groups="drop")
cstat_any_prod <- auc_mw(st$y, st$p_any)
cstat_any_mean <- auc_mw(st$y, st$p_mean)
cstat_tab <- data.frame(
  level=c("pooled_stay_x_window_rows_all","pooled_rows_evermeasured_stays",
          "stay_any_measured_1minusProd","stay_any_measured_meanP"),
  n_units=c(nrow(grid), nrow(ever), nrow(st), nrow(st)),
  n_events=c(sum(grid$measured), sum(ever$measured), sum(st$y), sum(st$y)),
  c_statistic=c(cstat_pooled, cstat_ever, cstat_any_prod, cstat_any_mean),
  model=paste0("measured ~ A*timef + ",paste(MMCOV,collapse=" + ")),
  note=c("Mann-Whitney AUC of fitted measurement probability; full balanced panel",
         "same, panel rows restricted to stays with >=1 measured window",
         "stay score = 1-prod_windows(1-pm) (illustrative aggregation)",
         "stay score = mean_windows(pm) (illustrative aggregation)"))
write.csv(cstat_tab, file.path(OUT,"pc_measurement_cstat.csv"), row.names=FALSE)
cat("=== M1 MEASUREMENT MODEL C-STATISTIC ===\n");print(cstat_tab)

# ---- M8: measured-vs-unmeasured SMD per window, before/after weighting ----
smd2 <- function(x1, x0, w1=NULL, w0=NULL){
  if(is.null(w1)){ m1<-mean(x1); v1<-var(x1); m0<-mean(x0); v0<-var(x0)
  } else { m1<-weighted.mean(x1,w1); v1<-sum(w1*(x1-m1)^2)/sum(w1)
           m0<-weighted.mean(x0,w0); v0<-sum(w0*(x0-m0)^2)/sum(w0) }
  (m1-m0)/sqrt((v1+v0)/2)
}
smd_rows <- list(); k<-1
for(t in c(6,12,24,48,72)){
  g <- grid %>% filter(time==t) %>% left_join(
         ed %>% select(stay_id,iptw), by="stay_id", suffix=c("",".ed"))
  # panel IPCW at this window, 1/99 trimmed (mirrors PC02)
  g$ipcw <- g$pm0/g$pm
  g$w <- g$iptw*g$ipcw
  wq <- quantile(g$w, c(.01,.99), na.rm=TRUE); g$wt <- pmin(pmax(g$w,wq[1]),wq[2])
  for(v in SMDV){
    x <- g[[v]]; mflag <- g$measured==1
    su  <- smd2(x[mflag], x[!mflag])
    # weighted: measured rows carry joint weight; target = all stays with IPTW
    sw  <- smd2(x[mflag], x, g$wt[mflag], g$iptw)
    smd_rows[[k]] <- data.frame(window=t, variable=v,
                                smd_unweighted=su, smd_weighted=sw); k<-k+1
  }
}
smd_tab <- bind_rows(smd_rows)
write.csv(smd_tab, file.path(OUT,"pc_measurement_smd.csv"), row.names=FALSE)
smd_sum <- smd_tab %>% group_by(window) %>%
  summarise(max_abs_smd_unweighted=max(abs(smd_unweighted)),
            max_abs_smd_weighted=max(abs(smd_weighted)), .groups="drop")
cat("\n=== M8 MEASUREMENT SMD SUMMARY (max |SMD|) ===\n");print(smd_sum)

# =====================================================================
# M2: weighted time-averaged phosphate concentration 0-48h
# (trapezoid AUC / last sampling time, mmol/L; NOT mmol*L^-1*h)
# =====================================================================
auc <- lng %>% filter(hph<=48) %>% arrange(stay_id,hph) %>% group_by(stay_id) %>%
  summarise(npts=n(), t0h=min(hph), tnh=max(hph),
            tac=sum(diff(hph)*(head(phos,-1)+tail(phos,-1))/2)/max(hph), .groups="drop") %>%
  filter(npts>=3, t0h<=12, tnh>=30) %>%
  inner_join(co %>% select(stay_id,exposure_w48,all_of(DAG14)), by="stay_id")
auc_crude <- auc %>% filter(exposure_w48 %in% c("A_IV","B_PO_only","D_none")) %>%
  group_by(exposure_w48) %>% summarise(n=n(), crude_mean=mean(tac), .groups="drop")

wmean_boot <- function(lab){
  d <- co %>% filter(exposure_w48 %in% c(lab,"D_none")) %>%
    filter(complete.cases(across(all_of(DAG14)))) %>% mutate(A=as.integer(exposure_w48==lab))
  z <- psw(d); d$w <- z$w
  el <- auc %>% inner_join(d %>% select(stay_id,A,w), by="stay_id")
  pt <- el %>% group_by(A) %>% summarise(n=n(), wm=weighted.mean(tac,w), .groups="drop")
  R <- matrix(NA_real_, B, 3)
  uid <- unique(d$stay_id)
  for(b in seq_len(B)){
    db <- d[sample(seq_len(nrow(d)), nrow(d), replace=TRUE), ]
    fb <- try(psw(db), silent=TRUE); if(inherits(fb,"try-error")) next
    db$wb <- fb$w
    eb <- el %>% select(stay_id,tac) %>% inner_join(db %>% select(stay_id,wb), by="stay_id") %>%
      inner_join(db %>% select(stay_id,Ab=A), by="stay_id")
    R[b,1] <- weighted.mean(eb$tac[eb$Ab==1], eb$wb[eb$Ab==1])
    R[b,2] <- weighted.mean(eb$tac[eb$Ab==0], eb$wb[eb$Ab==0])
    R[b,3] <- R[b,1]-R[b,2]
  }
  data.frame(comparison=ifelse(lab=="A_IV","C1_IV_vs_none","C3_PO_vs_none"),
             arm=c(paste0(lab,"_weighted"),"D_none_weighted","diff_minus_none"),
             n=c(pt$n[pt$A==1], pt$n[pt$A==0], NA_integer_),
             est=c(pt$wm[pt$A==1], pt$wm[pt$A==0],
                   pt$wm[pt$A==1]-pt$wm[pt$A==0]),
             lcl=c(q025(R[,1]),q025(R[,2]),q025(R[,3])),
             ucl=c(q975(R[,1]),q975(R[,2]),q975(R[,3])),
             ci_method=paste0("bootstrap percentile, B=",B))
}
aucw <- bind_rows(wmean_boot("A_IV"), wmean_boot("B_PO_only"))
write.csv(aucw, file.path(OUT,"pc_auc_weighted.csv"), row.names=FALSE)
cat("\n=== M2 WEIGHTED TIME-AVERAGED CONCENTRATION (mmol/L) ===\n")
print(as.data.frame(aucw)); print(as.data.frame(auc_crude))

# =====================================================================
# M5: safety / normalization weighted RR + RD, B=200 bootstrap CI
# =====================================================================
events <- c("hyperphos","hypok","hypomg","hypoca","hypoca_severe")
denmap <- c(hyperphos="phos_n48",hypok="k_n",hypomg="mg_n",hypoca="ca_n",hypoca_severe="ca_n")
safety_ci <- function(trt, comp){
  # ONE stabilized PS model on the full complete-case two-arm cohort (as PC03);
  # denominator subsets only select rows for rate computation
  d0 <- co %>% filter(exposure_w48 %in% c(trt,"D_none")) %>%
    inner_join(saf, by="stay_id") %>%
    filter(complete.cases(across(all_of(DAG14)))) %>%
    mutate(A=as.integer(exposure_w48==trt))
  w48 <- win %>% filter(win=="w48") %>% select(stay_id, phos48=phos)
  dn <- d0 %>% left_join(w48, by="stay_id") %>% mutate(norm48=as.integer(phos48>=0.80))
  z0 <- psw(dn); dn$w <- z0$w
  specs <- c(as.list(denmap), list(norm48=NA_character_))
  sub <- lapply(names(specs), function(v){
    if(v=="norm48") dn %>% filter(!is.na(phos48))
    else dn %>% filter(.data[[specs[[v]]]]>0) })
  names(sub) <- names(specs)
  pt <- lapply(names(sub), function(v){ dd<-sub[[v]]; y<-dd[[v]]
    r1<-weighted.mean(y[dd$A==1],dd$w[dd$A==1]); r0<-weighted.mean(y[dd$A==0],dd$w[dd$A==0])
    data.frame(outcome=v,n_trt=sum(dd$A==1),n_ctrl=sum(dd$A==0),
               rate_trt=r1,rate_ctrl=r0,RD=r1-r0,RR=r1/r0) }) %>% bind_rows()
  R <- array(NA_real_, dim=c(B, length(sub), 4),
             dimnames=list(NULL, names(sub), c("rate1","rate0","RD","RR")))
  for(b in seq_len(B)){
    db <- dn[sample(seq_len(nrow(dn)), nrow(dn), replace=TRUE), ]
    fb <- try(psw(db), silent=TRUE); if(inherits(fb,"try-error")) next
    db$wb <- fb$w
    for(v in names(sub)){
      dd <- if(v=="norm48") db %>% filter(!is.na(phos48)) else db %>% filter(.data[[specs[[v]]]]>0)
      m1<-weighted.mean(dd[[v]][dd$A==1],dd$wb[dd$A==1])
      m0<-weighted.mean(dd[[v]][dd$A==0],dd$wb[dd$A==0])
      R[b,v,] <- c(m1,m0,m1-m0,m1/m0)
    }
  }
  ci <- lapply(names(sub), function(v) data.frame(
    outcome=v, RD_lcl=q025(R[,v,"RD"]), RD_ucl=q975(R[,v,"RD"]),
    RR_lcl=q025(R[,v,"RR"]), RR_ucl=q975(R[,v,"RR"]))) %>% bind_rows()
  pt %>% left_join(ci, by="outcome") %>% mutate(
    comparison=comp, ci_method=paste0("bootstrap percentile, B=",B,
      "; PS refit per replicate on full complete-case cohort")) %>%
    select(comparison,outcome,n_trt,n_ctrl,rate_trt,rate_ctrl,
           RD,RD_lcl,RD_ucl,RR,RR_lcl,RR_ucl,ci_method)
}
saf_ci <- bind_rows(safety_ci("A_IV","C1_IV_vs_none"),
                    safety_ci("B_PO_only","C3_PO_vs_none"))
write.csv(saf_ci, file.path(OUT,"pc_safety_weighted_ci.csv"), row.names=FALSE)
cat("\n=== M5 SAFETY WEIGHTED RR/RD WITH BOOTSTRAP CI ===\n")
print(as.data.frame(saf_ci %>% select(comparison,outcome,rate_trt,rate_ctrl,
          RD,RD_lcl,RD_ucl,RR,RR_lcl,RR_ucl)), digits=3)

# =====================================================================
# M6: NC2 pre-T0 slope difference CI + GEE A*time interaction p
# =====================================================================
pw <- pre %>% pivot_wider(names_from=win, values_from=c(phos,actual_hph)) %>%
  filter(!is.na(phos_pre24), !is.na(phos_pre6)) %>%
  mutate(slope=(phos_pre6-phos_pre24)/(actual_hph_pre6-actual_hph_pre24)) %>%
  inner_join(co, by="stay_id") %>%
  filter(exposure_w48 %in% c("A_IV","D_none"),
         complete.cases(across(all_of(DAG14)))) %>% mutate(A=as.integer(exposure_w48=="A_IV"))
z <- psw(pw); pw$w <- z$w
sl <- pw %>% group_by(A) %>% summarise(m=weighted.mean(slope,w), .groups="drop")
R <- matrix(NA_real_, B, 3)
for(b in seq_len(B)){
  db <- pw[sample(seq_len(nrow(pw)), nrow(pw), replace=TRUE), ]
  fb <- try(psw(db), silent=TRUE); if(inherits(fb,"try-error")) next
  wb <- fb$w
  R[b,1]<-weighted.mean(db$slope[db$A==1],wb[db$A==1])
  R[b,2]<-weighted.mean(db$slope[db$A==0],wb[db$A==0])
  R[b,3]<-R[b,1]-R[b,2]
}
long <- bind_rows(
  pw %>% transmute(stay_id,A,w,time=actual_hph_pre24,phos=phos_pre24,baseline_phosphate),
  pw %>% transmute(stay_id,A,w,time=actual_hph_pre6, phos=phos_pre6, baseline_phosphate))
gfit <- geeglm(phos ~ A*time + baseline_phosphate, id=stay_id, weights=w,
               data=long, corstr="exchangeable")
av <- anova(gfit); int_p <- as.numeric(av[grepl("A:time",rownames(av)), "P(>|Chi|)"])
nc2ci <- data.frame(
  metric=c("slope_IV","slope_none","slope_diff_IV_minus_none"),
  est=c(sl$m[sl$A==1], sl$m[sl$A==0], sl$m[sl$A==1]-sl$m[sl$A==0]),
  lcl=c(q025(R[,1]),q025(R[,2]),q025(R[,3])),
  ucl=c(q975(R[,1]),q975(R[,2]),q975(R[,3])),
  n_IV=sum(pw$A==1), n_none=sum(pw$A==0),
  gee_A_time_interaction_p=int_p,
  ci_method=paste0("slope: bootstrap percentile B=",B,
                   "; interaction: weighted GEE Wald"))
write.csv(nc2ci, file.path(OUT,"pc_nc2_pretrend_ci.csv"), row.names=FALSE)
cat("\n=== M6 PRE-T0 SLOPE CI + GEE INTERACTION p ===\n");print(nc2ci)

# =====================================================================
# M7a: RCS nonlinearity tests (exact PC06 specification)
# =====================================================================
w12 <- win %>% filter(win=="w12") %>% select(stay_id,p12=phos)
w48 <- win %>% filter(win=="w48") %>% select(stay_id,p48=phos)
d <- co %>% filter(exposure_w48=="A_IV") %>% inner_join(dose, by="stay_id") %>%
  left_join(w12, by="stay_id") %>% left_join(w48, by="stay_id") %>%
  mutate(d12=p12-baseline_phosphate, d48=p48-baseline_phosphate)
cap <- quantile(d$dose24, 0.99, na.rm=TRUE)
d <- d %>% mutate(dose=pmin(dose24,cap))
dd <<- datadist(d); options(datadist="dd")
nl_p <- sapply(c("d12","d48"), function(out){
  da <- d %>% filter(!is.na(.data[[out]]), complete.cases(across(all_of(DAG14))))
  f  <- ols(as.formula(paste(out,"~ rcs(dose, c(7.5,15,30,60)) +",
                             paste(DAG14,collapse="+"))), data=da)
  a  <- anova(f)
  rn <- trimws(rownames(a))
  c(nonlinear_p=as.numeric(a[grep("Nonlinear",rn)[1], "P"]),
    total_p=as.numeric(a[grep("TOTAL",rn)[1], "P"]))
})
rcs_aov <- data.frame(outcome=colnames(nl_p), t(nl_p),
  test="Wald test of nonlinear spline component vs adjusted model")
write.csv(rcs_aov, file.path(OUT,"pc_dose_rcs_anova.csv"), row.names=FALSE)
cat("\n=== M7a RCS NONLINEARITY ===\n");print(rcs_aov)

# ---- TR-8.3: first-24h vs cumulative-48h dose metric ----------------
d48tab <- d %>% filter(!is.na(d48))
q48 <- d48tab %>% mutate(q=ntile(dose48,4)) %>% group_by(q) %>%
  summarise(n=n(), dose48_med=median(dose48), d48=mean(d48),
            norm48=100*mean(p48>=0.80), .groups="drop")
late0 <- d48tab %>% filter(dose24==0)
mc <- bind_rows(
  data.frame(section="counts", metric="n_AIV",                value=nrow(d)),
  data.frame(section="counts", metric="n_dose24_eq_0",         value=sum(d$dose24==0,na.rm=TRUE)),
  data.frame(section="counts", metric="n_dose24_0_and_dose48_gt0", value=sum(d$dose24==0 & d$dose48>0,na.rm=TRUE)),
  data.frame(section="counts", metric="n_dose48_eq_0",         value=sum(d$dose48==0,na.rm=TRUE)),
  data.frame(section="quantiles", metric="dose24_median",      value=median(d$dose24,na.rm=TRUE)),
  data.frame(section="quantiles", metric="dose24_p99",         value=quantile(d$dose24,.99,na.rm=TRUE)),
  data.frame(section="quantiles", metric="dose48_median",      value=median(d$dose48,na.rm=TRUE)),
  data.frame(section="quantiles", metric="dose48_p99",         value=quantile(d$dose48,.99,na.rm=TRUE)),
  data.frame(section="quantiles", metric="spearman_dose24_dose48",
             value=suppressWarnings(cor(d$dose24,d$dose48,use="complete.obs",method="spearman"))),
  data.frame(section="late0_first24h_none", metric="n",            value=nrow(late0)),
  data.frame(section="late0_first24h_none", metric="mean_d48_delta", value=mean(late0$d48)),
  data.frame(section="late0_first24h_none", metric="mean_p48",     value=mean(late0$p48)),
  data.frame(section="late0_first24h_none", metric="norm48_pct",   value=100*mean(late0$p48>=0.80)))
write.csv(mc, file.path(OUT,"pc_dose_metric_compare.csv"), row.names=FALSE)
write.csv(q48, file.path(OUT,"pc_dose48_quartiles.csv"), row.names=FALSE)
cat("\n=== TR-8.3 DOSE METRIC COMPARISON ===\n");print(mc);print(q48)

# =====================================================================
# M3: concordance with primary-analysis ITT treatment assignment
# =====================================================================
write.csv(conc, file.path(OUT,"pc_concordance_primary.csv"), row.names=FALSE)
cat("\n=== M3 CONCORDANCE WITH PRIMARY ITT ASSIGNMENT ===\n");print(conc)
cat("\nDONE.\n")
