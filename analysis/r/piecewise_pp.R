#!/usr/bin/env Rscript
# ============================================================================
# V28.4 分段HR + Per-Protocol(PP) 分析（stay 唯一队列，DAG14）
#
# 口径（2026-10-10 实证）：
#   * 治疗臂 3249 人全部 12h 内启动 IV -> PP 全依从，无删失
#   * 对照臂 19436 人中 2291 人随后交叉到 IV（中位 23.4h，<=48h 窗）-> PP 偏离，
#     在交叉日做人工删失；8936 完全保守、8209 仅口服（仅口服不违反 IV 策略）
#   * 死亡日：精确 2552 落真实日；仅知死亡 1076 在[出院日,28]区间落点
#
# 方法：
#   1) DAG14 估 PS -> 稳定化治疗权重（1/99 截断）
#   2) PP：交叉(对照臂启动IV)在交叉日删失，DAG14×interval  pooled logistic
#      估条件删失概率，累积得 IPC 权重，乘 PS 权重
#   3) 分段 HR：0-3/3-7/7-14/14-28，加权 Cox（段内 entry/exit，风险集=段初存活）
#   4) PP 总体 HR：0..28 加权 Cox（含交叉删失）
#   m=20 点估计 Rubin 合并；B=200 MI-bootstrap（subject 重抽 × 区间落点 ×
#      重估 PS 与删失权重），bootstrap 百分位 CI
#
# 用法:
#   验证: Rscript piecewise_pp.R --imps 1 --nboot 0 --tag imp1
#   全量: Rscript piecewise_pp.R --imps 1..20 --nboot 200
# 须 LD_LIBRARY_PATH=/opt/miniforge3/lib   V28.4  2026-10-10
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(survival)
})

BASE   <- "/data/小米云盘/文档/论文/血磷"
AUG_F  <- file.path(BASE,"TTE/v26/output/sequential_trials_augmented.parquet")
MICE_F <- file.path(BASE,"sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet")
OUT_T  <- file.path(BASE,"TTE/v28.4/output/tables/piecewise_pp")
dir.create(OUT_T, recursive=TRUE, showWarnings=FALSE)
TAU<-28L; SEED<-20260629L

argv <- commandArgs(trailingOnly=TRUE)
gv <- function(k,d){i<-match(k,argv);if(is.na(i))d else argv[i+1]}
IMPS  <- if("--imps" %in% argv) as.integer(strsplit(gv("--imps","1"),"\\.\\.")[[1]]) |>
           (\(x) if(length(x)==2) seq(x[1],x[2]) else as.integer(strsplit(gv("--imps","1"),",")[[1]]))() else 1L
NBOOT <- as.integer(gv("--nboot","0")); TAG <- gv("--tag","run")
N_PLACE<- 20L

BIN <- c("sepsis3","invasive_vent_at_t0","vasoactive_at_t0","refeeding_risk")
CONT<- c("age","baseline_sofa","charlson","egfr","baseline_phosphate",
         "lactate_t0","corrected_calcium_mmol_min")
FACT<- c("gender","admission_type","nutrition_support_type")
LEV <- list(gender=c("F","M"),
            admission_type=c("Observation","Emergency","Surgical_Elective"),
            nutrition_support_type=c("none","enteral","parenteral","mixed"))
DAG <- c(CONT,BIN,FACT)
IEDGES <- c(0L,3L,7L,14L,28L)
IIDS   <- 1:4
ILAB   <- c("0-3d","3-7d","7-14d","14-28d")

# ---- person 表（stay 唯一） ----
load_persons <- function(imp_id){
  aug <- as.data.table(read_parquet(AUG_F))
  imp <- as.data.table(read_parquet(sprintf(MICE_F,imp_id)))
  a <- aug[order(stay_id,trial_num),.SD[1],by=stay_id]
  m <- imp[order(stay_id,trial_num),.SD[1],by=stay_id]
  stopifnot(identical(a$stay_id,m$stay_id))
  p <- a[,.(stay_id,subject_id,death=as.integer(death_28day_num))]
  p[, `:=`(arm=as.integer(a$itt_group=="treatment"),
           exact_t=as.numeric(a$death_time_from_t0),
           disch=as.numeric(a$dischtime_from_t0),
           cross_h=as.numeric(a$first_iv_hours_from_tz))]
  for(v in c(CONT,BIN)) p[,(v):=m[[v]]]
  for(v in BIN) p[,(v):=as.integer(as.character(get(v)))]
  for(v in FACT) p[,(v):=factor(a[[v]],levels=LEV[[v]])]
  # 死亡日：先建列，精确死亡落真实日
  p[, death_day := NA_integer_]
  p[, has_exact := death==1L & !is.na(exact_t)]
  p[has_exact==TRUE, death_day := pmin(TAU,pmax(0L,as.integer(ceiling(exact_t))))]
  p[, dlower := as.integer(pmin(TAU,pmax(0L,ceiling(disch))))]
  p[is.na(dlower), dlower:=0L]
  # 交叉日（仅对照臂会交叉）：first_iv_hours_from_tz 换算到日(ceil)
  p[, cross_day := NA_integer_]
  ix <- which(p$arm==0L & !is.na(p$cross_h))
  p$cross_day[ix] <- pmin(TAU,pmax(1L,as.integer(ceiling(p$cross_h[ix]/24))))
  p
}

# ---- 死亡日区间落点 ----
assign_placement <- function(p){
  pp<-copy(p)
  idx<-which(pp$death==1L & !pp$has_exact)
  if(length(idx)){
    lo<-pp$dlower[idx]
    pp$death_day[idx]<-lo+as.integer(floor(runif(length(idx))*(TAU-lo+1L)))
  }
  pp$death_day[is.na(pp$death_day)]<-TAU+1L
  pp
}

# ---- DAG14 稳定化治疗权重 ----
add_ps <- function(p){
  fml<-as.formula(paste("arm ~",paste(DAG,collapse=" + ")))
  ps<-predict(suppressWarnings(glm(fml,data=p,family=binomial(),
                                   control=list(maxit=100))),type="response")
  ps[is.na(ps)]<-mean(p$arm); pa<-mean(p$arm)
  w<-ifelse(p$arm==1,pa/ps,(1-pa)/(1-ps))
  p[,wps:=w]
  p
}

# ---- 构造 person-interval 长表 + PP 交叉删失权重 ----
build_long <- function(p){
  pp<-copy(p)
  pp[, pid:=.I]
  pp[, cen_day:=pmin(death_day,TAU)]
  cc<-which(pp$arm==0L & !is.na(pp$cross_day) & pp$cross_day < pp$death_day)
  pp$cen_day[cc]<-pp$cross_day[cc]
  pp[, censored:=0L]; pp$censored[cc]<-1L
  rows<-lapply(seq_len(nrow(pp)),function(k){
    cut<-pp$cen_day[k]
    ii<-which(IEDGES[-length(IEDGES)] < cut)
    if(!length(ii)) ii<-which.min(IEDGES[-length(IEDGES)])
    data.table(pid=k,iid=ii,L=IEDGES[ii],R=IEDGES[ii+1L],cut=cut,
               cens=pp$censored[k],arm=pp$arm[k])
  })
  L<-rbindlist(rows)
  L[, c_event:=as.integer(cens==1L & cut>L & cut<=R)]
  for(v in DAG) L[, (v):=pp[[v]][pid]]
  L[, ifac:=factor(iid,levels=IIDS)]
  L[, Sc:=1]
  if(any(L$arm==0L & L$c_event==1L)){
    d0<-L[arm==0L & L<cut]
    cf<-as.formula(paste("c_event ~ ifac +",paste(DAG,collapse=" + ")))
    fit<-suppressWarnings(glm(cf,data=d0,family=binomial(),control=list(maxit=100)))
    d0[, hc:=predict(fit,newdata=d0,type="response")]
    setorder(d0,pid,iid)
    d0[, Sc:=c(1,cumprod(1-head(hc,-1))),by=pid]
    L<-merge(L,d0[,.(pid,iid,Sc_map=Sc)],by=c("pid","iid"),all.x=TRUE,sort=FALSE)
    L[!is.na(Sc_map), Sc:=Sc_map]; L[, Sc_map:=NULL]
  }
  base<-pp[,.(pid,stay_id,subject_id,arm,death_day,wps,cen_day,censored)]
  ML<-merge(base,L,by=c("pid","arm"),allow.cartesian=TRUE)
  ML[, wfinal:=wpos<-wps/Sc]
  qs<-quantile(ML$wfinal,c(.01,.99),na.rm=TRUE)
  ML[, wfinal:=pmin(pmax(wfinal,qs[1]),qs[2])]
  ML
}

# ---- PP 总体 Cox（0..28，含交叉删失） ----
pp_overall<-function(ML){
  setorder(ML,stay_id,iid)
  # 每 stay 取随访覆盖的最后一个 interval 的累积删失存活 Sc
  d<-ML[, .SD[.N], by=stay_id][,.(stay_id,subject_id,arm,death_day,
                                  cen_day,censored,wps,Sc)]
  d[, evt:=as.integer(death_day<=cen_day)]
  wv<-d$wps/d$Sc
  fit<-coxph(Surv(cen_day,evt)~arm,data=d,weights=wv,
             cluster=subject_id,robust=TRUE,control=coxph.control(iter.max=50))
  s<-summary(fit)
  c(beta=s$coefficients[1,"coef"],se=s$coefficients[1,"se(coef)"])
}

# ---- 分段 Cox（4 段，段内 entry/exit，风险集=段初存活） ----
piecewise<-function(ML){
  out<-list()
  for(ii in IIDS){
    l<-IEDGES[ii]; r<-IEDGES[ii+1L]
    d<-copy(ML[iid==ii])
    # 风险集：interval 起始 L 仍存活（death_day>L）
    d<-d[death_day>l]
    d[, t_in:=pmin(pmax(l,0),r-1L)]
    d[, t_out:=pmin(death_day,cen_day,r)]
    d[, evt:=as.integer(death_day<=pmin(cen_day,r))]
    wv<-d$wps/d$Sc
    qs<-quantile(wv,c(.01,.99),na.rm=TRUE); wv<-pmin(pmax(wv,qs[1]),qs[2])
    fit<-tryCatch(coxph(Surv(t_in,t_out,evt)~arm,data=d,weights=wv,
                  cluster=subject_id,robust=TRUE,
                  control=coxph.control(iter.max=50)),error=function(e)NULL)
    if(is.null(fit)){out[[ii]]<-c(beta=NA_real_,se=NA_real_,nrisk=NA_integer_,nevt=NA_integer_)}
    else{s<-summary(fit)
      out[[ii]]<-c(beta=s$coefficients[1,"coef"],se=s$coefficients[1,"se(coef)"],
                   nrisk=nrow(d),nevt=sum(d$evt))}
  }
  do.call(rbind,out)
}

# ---- Rubin 合并 logHR ----
rubin<-function(betas,ses){
  k<-!is.na(betas)&!is.na(ses);betas<-betas[k];ses<-ses[k];m<-length(betas)
  if(m==0)return(c(hr=NA,ci_l=NA,ci_u=NA,p=NA))
  bb<-mean(betas);W<-mean(ses^2);B<-if(m>1)sum((betas-bb)^2/(m-1))else 0
  sep<-sqrt(W+(1+1/m)*B);t<-bb/sep
  df<-if(m>1)(m-1)*(1+W/((1+1/m)*B))^2 else 1
  p<-2*pt(abs(t),df,lower.tail=FALSE)
  c(hr=exp(bb),ci_l=exp(bb-1.96*sep),ci_u=exp(bb+1.96*sep),p=p)
}

# ============ 点估计（跨 imp × 落点平均） ============
set.seed(SEED)
persons<-lapply(IMPS,load_persons)
cat("loaded m =",length(persons),"\n")

est_one<-function(p0){
  pw<-add_ps(p0)
  pp<-assign_placement(pw)
  ML<-build_long(pp)
  list(pp=pp_overall(ML),pw=piecewise(ML))
}
runs<-list()
for(p0 in persons){
  for(s in seq_len(N_PLACE)) runs[[length(runs)+1]]<-est_one(p0)
}
pp_b<-sapply(runs,function(z)z$pp["beta"]); pp_s<-sapply(runs,function(z)z$pp["se"])
pw_b<-sapply(runs,function(z)z$pw[,"beta"]); pw_s<-sapply(runs,function(z)z$pw[,"se"])
# 风险集/事件数取一次代表性落点（runs[[1]]），不跨落点平均
nrisk<-runs[[1]]$pw[,"nrisk"]; nevt<-runs[[1]]$pw[,"nevt"]

cat("\n===== 点估计 PP 总体（Rubin over imp×placement） =====\n")
pp_pt<-rubin(pp_b,pp_s); print(round(pp_pt,4))
cat("\n===== 点估计 分段 HR =====\n")
pw_pt<-t(sapply(1:4,function(ii)rubin(pw_b[ii,],pw_s[ii,])))
rownames(pw_pt)<-ILAB
pw_dt<-data.table(segment=ILAB,t_start=IEDGES[-5],t_end=IEDGES[-1],
                  n_at_risk=nrisk,n_events=nevt,as.data.table(pw_pt))
print(pw_dt)

# ============ MI-bootstrap ============
if(NBOOT>0){
  Bm<-matrix(NA_real_,NBOOT,5); colnames(Bm)<-c("pp_beta","b1","b2","b3","b4")
  nb<-0L
  while(nb<NBOOT){
    p0<-persons[[sample(length(persons),1)]]
    pick<-sample(unique(p0$subject_id),uniqueN(p0$subject_id),replace=TRUE)
    pb<-p0[data.table(subject_id=pick),on="subject_id",allow.cartesian=TRUE]
    r<-tryCatch(est_one(add_ps(pb)),error=function(e)NULL)
    if(is.null(r)||any(!is.finite(r$pp[1]),!is.finite(r$pw[,"beta"])))next
    nb<-nb+1L
    Bm[nb,]<-c(r$pp["beta"],r$pw[,"beta"])
    if(nb%%50==0)cat(sprintf("  boot %d/%d\n",nb,NBOOT))
  }
  ci<-apply(Bm,2,quantile,c(.025,.975),na.rm=TRUE)
  # PP bootstrap 结果
  pp_final<-data.table(analysis="Per-protocol (IPCW-crossover, DAG14)",
                       hr=exp(mean(Bm[,"pp_beta"])),
                       ci_lower=exp(ci[1,"pp_beta"]),ci_upper=exp(ci[2,"pp_beta"]))
  # 分段 bootstrap
  seg_final<-rbindlist(lapply(1:4,function(ii){
    data.table(segment=ILAB[ii],t_start=IEDGES[ii],t_end=IEDGES[ii+1],
               n_at_risk=nrisk[ii],n_events=nevt[ii],
               hr=exp(mean(Bm[,1+ii])),
               ci_lower=exp(ci[1,1+ii]),ci_upper=exp(ci[2,1+ii]))
  }))
  cat("\n===== PP 总体 Bootstrap 95% CI =====\n"); print(pp_final)
  cat("\n===== 分段 HR Bootstrap 95% CI =====\n"); print(seg_final)
  fwrite(pp_final,file.path(OUT_T,paste0("pp_",TAG,".csv")))
  fwrite(seg_final,file.path(OUT_T,paste0("piecewise_",TAG,".csv")))
  fwrite(as.data.table(Bm),file.path(OUT_T,paste0("piecewise_pp_bootraw_",TAG,".csv")))
}else{
  fwrite(data.table(analysis="Per-protocol",t(pp_pt)),
         file.path(OUT_T,paste0("pp_",TAG,".csv")))
  fwrite(pw_dt,file.path(OUT_T,paste0("piecewise_",TAG,".csv")))
}
cat("\ndone tag =",TAG,"\n")
