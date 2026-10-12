#!/usr/bin/env Rscript
# ============================================================================
# V28.4 Figure 2 — 两套候选估计量（stay 唯一队列，每 stay 取首个试验=t0）
#   (1) 观察口径：每臂真实观察到的 28d 死亡风险 + crude RMST
#   (2) 回归标准化 g-computation：pooled logistic 离散风险模型预测"若全程
#       随访到28天"的臂特异累计发病曲线（不依赖出院无信息删失）
# 用法:
#   验证:   Rscript fig2_estimands.R --imps 1 --nboot 0 --tag imp1
#   全量:   Rscript fig2_estimands.R --imps 1..20 --nboot 200
# 须 LD_LIBRARY_PATH=/opt/miniforge3/lib   版本 V28.4  2026-10-10
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(survival)
})

BASE <- "/data/小米云盘/文档/论文/血磷"
AUG_F <- file.path(BASE,"TTE/v26/output/sequential_trials_augmented.parquet")
MICE_F <- file.path(BASE,"sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet")
OUT_DIR <- file.path(BASE,"TTE/v28.4/output/tables/fig2")
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)
TAU <- 28L; SEED <- 20260629L

argv <- commandArgs(trailingOnly=TRUE)
gv <- function(k,d){i<-match(k,argv);if(is.na(i))d else argv[i+1]}
IMPS <- as.integer(strsplit(gv("--imps","1"),",")[[1]])
NBOOT <- as.integer(gv("--nboot","0"))
TAG <- gv("--tag","run")

BIN <- c("sepsis3","invasive_vent_at_t0","vasoactive_at_t0","refeeding_risk")
CONT <- c("age","baseline_sofa","charlson","egfr","baseline_phosphate",
          "lactate_t0","corrected_calcium_mmol_min")
FACT <- c("gender","admission_type","nutrition_support_type")
LEV <- list(gender=c("F","M"),
            admission_type=c("Observation","Emergency","Surgical_Elective"),
            nutrition_support_type=c("none","enteral","parenteral","mixed"))

# ---- 构建某 imp 的 stay 唯一 person 表 ----
load_persons <- function(imp_id){
  aug <- as.data.table(read_parquet(AUG_F))
  imp <- as.data.table(read_parquet(sprintf(MICE_F,imp_id)))
  setkey(aug,stay_id,trial_num); setkey(imp,stay_id,trial_num)
  a <- aug[order(stay_id,trial_num), .SD[1], by=stay_id]
  m <- imp[order(stay_id,trial_num), .SD[1], by=stay_id]
  stopifnot(identical(a$stay_id, m$stay_id))
  p <- a[, .(stay_id, subject_id, itt_group, death=as.integer(death_28day_num),
             t=as.numeric(surv_time))]
  p[, arm := as.integer(as.character(m$treatment))]
  for(v in c(CONT,BIN)) p[, (v) := m[[v]]]
  for(v in BIN) p[, (v) := as.integer(as.character(get(v)))]
  for(v in FACT) p[, (v) := factor(a[[v]], levels=LEV[[v]])]
  p[, tday := pmin(TAU, pmax(0L, as.integer(ceiling(t))))]
  p
}

# ---- 观察口径：风险 + crude RMST（死亡/出院时间） ----
observed <- function(p){
  risk <- p[, .(n=.N, deaths=sum(death), risk=mean(death)), by=arm]
  rmst <- sapply(c(1,0), function(aa){
    sf <- survfit(Surv(pmin(t,TAU), death)~1, data=p[arm==aa])
    summary(sf, rmean=TAU)$table["rmean"]
  })
  list(risk=risk, rmst=data.table(arm=c(1,0), rmst=as.numeric(rmst)))
}

# ---- 构建离散风险 person-day ----
# 每人含 day=0..last：存活者 last=tday-1（出院，全为0）；
# 死亡者 last=tday（死亡日当天 y=1）
make_long <- function(p){
  last <- ifelse(p$death==1L, p$tday, p$tday - 1L)
  long <- data.table(stay_id=rep(p$stay_id, last+1L),
                     day=as.integer(unlist(lapply(last, seq, from=0L))),
                     y=0L)
  ev <- p[death==1, .(stay_id, d=tday)]
  long[ev, on="stay_id", y := as.integer(day==i.d)]
  cols <- c("stay_id","arm",CONT,BIN,FACT)
  long <- merge(long, p[, cols, with=FALSE], by="stay_id")
  long
}

# ---- g-computation：拟合离散风险模型，标准化 ----
# RMST(tau) = sum_{k=0..tau-1} S(k)（日为整数，S(0)=1）
gcomp <- function(p){
  long <- make_long(p)
  long[, tf := factor(day)]
  covs <- paste(c(CONT,BIN,FACT), collapse=" + ")
  fml <- as.formula(paste("y ~ arm + tf +", covs, "+ arm:tf"))
  fit <- suppressWarnings(glm(fml, data=long, family=binomial(),
                              control=list(maxit=50)))
  base <- p[, c("stay_id",CONT,BIN,FACT), with=FALSE]
  day_seq <- seq(0L, TAU - 1L)
  out <- list()
  for(aa in c(1L,0L)){
    n <- nrow(base)
    pr <- base[rep(seq_len(n), each=TAU)]
    pr[, arm := aa]
    pr[, day := rep(day_seq, n)]
    pr[, tf := factor(day, levels=levels(long$tf))]
    pr[, h := pmin(pmax(predict(fit,newdata=pr,type="response"),1e-6),1-1e-6)]
    setorder(pr,stay_id,day)
    pr[, Sday := cumprod(1-h), by=stay_id]
    ci <- pr[, .(S=mean(Sday)), by=day][order(day)]
    ci[, `:=`(cif=1-S, arm=aa)]
    out[[length(out)+1]] <- ci
  }
  curve <- rbindlist(out)
  rmst <- curve[, .(rmst=sum(S)), by=arm]
  list(curve=curve[,.(arm,day,cif,S)], rmst=rmst)
}

# ---- 一次计算，返回结构化向量（供点估计与 bootstrap 复用） ----
one_run <- function(p){
  o <- observed(p); g <- gcomp(p)
  list(
    obs_risk = o$risk,
    obs_rmst = o$rmst,
    g_curve  = g$curve,
    g_rmst   = g$rmst)
}

# 从 one_run 结果提取关键对比（IV - control）
extract_contrasts <- function(r){
  or <- r$obs_risk
  c(obs_risk_t = or[arm==1]$risk, obs_risk_c = or[arm==0]$risk,
     obs_rd   = or[arm==1]$risk - or[arm==0]$risk,
     obs_rmst_t = r$obs_rmst[arm==1]$rmst, obs_rmst_c = r$obs_rmst[arm==0]$rmst,
     obs_rmst_d = r$obs_rmst[arm==1]$rmst - r$obs_rmst[arm==0]$rmst,
     g_risk_t = r$g_curve[arm==1 & day==TAU-1L]$cif,
     g_risk_c = r$g_curve[arm==0 & day==TAU-1L]$cif,
     g_rd     = r$g_curve[arm==1 & day==TAU-1L]$cif - r$g_curve[arm==0 & day==TAU-1L]$cif,
     g_rmst_t = r$g_rmst[arm==1]$rmst, g_rmst_c = r$g_rmst[arm==0]$rmst,
     g_rmst_d = r$g_rmst[arm==1]$rmst - r$g_rmst[arm==0]$rmst)
}

# ============ 执行（点估计：跨 imp 平均） ============
set.seed(SEED)
runs <- list()
for(i in seq_along(IMPS)) runs[[i]] <- one_run(load_persons(IMPS[i]))

# 汇总
mk <- function(elem, col){
  rbindlist(lapply(runs, `[[`, elem))[, setNames(list(mean(get(col))), col), by=arm]
}
obs_risk <- rbindlist(lapply(runs, `[[`, "obs_risk"))[,
              .(n=mean(n), deaths=mean(deaths), risk=mean(risk)), by=arm]
obs_rmst <- mk("obs_rmst","rmst")
gc_curve <- rbindlist(lapply(runs, `[[`,"g_curve"))[, .(cif=mean(cif),S=mean(S)), by=.(arm,day)]
gc_rmst  <- mk("g_rmst","rmst")

cat("\n===== (1) 观察口径 =====\n"); print(obs_risk); print(obs_rmst)
cat("\n===== (2) g-computation 日28 =====\n"); print(gc_curve[day==TAU-1L]); print(gc_rmst)

# ---- Bootstrap：MI-bootstrap（随机 imp × subject 重抽） ----
person_pool <- lapply(IMPS, function(id) load_persons(id))
if(NBOOT > 0){
  boot_mat <- matrix(NA_real_, nrow=NBOOT, ncol=12)
  k <- 0L
  while(k < NBOOT){
    ii <- sample(length(person_pool),1)
    p0 <- person_pool[[ii]]
    pick <- sample(unique(p0$subject_id), uniqueN(p0$subject_id), replace=TRUE)
    pb <- p0[data.table(subject_id=pick), on="subject_id", allow.cartesian=TRUE]
    r <- tryCatch(extract_contrasts(one_run(pb)), error=function(e) NULL)
    if(is.null(r) || any(!is.finite(r))) next
    k <- k+1L; boot_mat[k,] <- r
    if(k%%50==0) cat(sprintf("  boot %d/%d\n",k,NBOOT))
  }
  cnm <- names(extract_contrasts(one_run(person_pool[[1]])))
  colnames(boot_mat) <- cnm
  ci <- apply(boot_mat,2,quantile,c(.025,.975),na.rm=TRUE)
  cat("\n===== Bootstrap 95% CI（B=",NBOOT,"） =====\n",sep="")
  print(round(t(rbind(est=colMeans(boot_mat), lci=ci[1,], uci=ci[2,])),4))
  fwrite(as.data.table(boot_mat, keep.rownames=FALSE),
         file.path(OUT_DIR,paste0("fig2_boot_",TAG,".csv")))
}

fwrite(gc_curve, file.path(OUT_DIR,paste0("fig2_gcomp_curve_",TAG,".csv")))
fwrite(obs_risk, file.path(OUT_DIR,paste0("fig2_obs_risk_",TAG,".csv")))
cat("\nsaved tag=",TAG,"\n")
