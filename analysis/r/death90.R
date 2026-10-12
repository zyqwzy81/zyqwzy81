#!/usr/bin/env Rscript
# ============================================================================
# V28.4 90 天死亡分析（stay 唯一队列，DAG14 稳定化 IPTW）
#
# 口径（2026-10-10 实证）：
#   * 死亡时间源 = mimiciv_hosp.patients.dod（可独立核验、院内外统一死亡日期）
#     死亡日 = dod - DATE(time_zero)；仅取落在 [0,90] 的死亡
#   * parquet 的 death_90day（4890）为更宽口径：多出 91 例其 dod 落在
#     91..96 天/更晚/缺失（52），本分析不采用，改以 patients.dod 为准 = 4799
#   * 未在 90 天内记录死亡者：第 90 日行政删失
# 方法：DAG14 估 PS -> 稳定化 IPTW（1/99 截断）->
#        加权组间风险 / RD / RR / OR + 加权 Cox 总体 HR（cluster=subject）
#        m=20 点估计 Rubin 合并；B=200 MI-bootstrap（subject 重抽 × 重估 PS）
# 用法:
#   验证: Rscript death90.R --imps 1 --nboot 0 --tag imp1
#   全量: Rscript death90.R --imps 1..20 --nboot 200
# 须 LD_LIBRARY_PATH=/opt/miniforge3/lib   V28.4  2026-10-10
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(survival)
})

BASE   <- "/data/小米云盘/文档/论文/血磷"
AUG_F  <- file.path(BASE,"TTE/v26/output/sequential_trials_augmented.parquet")
MICE_F <- file.path(BASE,"sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet")
SRC_F  <- file.path(BASE,"TTE/v28.4/output/tables/death90/death90_time_source.csv")
OUT_T  <- file.path(BASE,"TTE/v28.4/output/tables/death90")
dir.create(OUT_T, recursive=TRUE, showWarnings=FALSE)
TAU<-90L; SEED<-20260629L

argv <- commandArgs(trailingOnly=TRUE)
gv <- function(k,d){i<-match(k,argv);if(is.na(i))d else argv[i+1]}
IMPS  <- if("--imps" %in% argv) as.integer(strsplit(gv("--imps","1"),"\\.\\.")[[1]]) |>
           (\(x) if(length(x)==2) seq(x[1],x[2]) else as.integer(strsplit(gv("--imps","1"),",")[[1]]))() else 1L
NBOOT <- as.integer(gv("--nboot","0")); TAG <- gv("--tag","run")

BIN <- c("sepsis3","invasive_vent_at_t0","vasoactive_at_t0","refeeding_risk")
CONT<- c("age","baseline_sofa","charlson","egfr","baseline_phosphate",
         "lactate_t0","corrected_calcium_mmol_min")
FACT<- c("gender","admission_type","nutrition_support_type")
LEV <- list(gender=c("F","M"),
            admission_type=c("Observation","Emergency","Surgical_Elective"),
            nutrition_support_type=c("none","enteral","parenteral","mixed"))
DAG <- c(CONT,BIN,FACT)

src <- fread(SRC_F)[,.(stay_id,d90day=day_from_pdod)]

load_persons <- function(imp_id){
  aug <- as.data.table(read_parquet(AUG_F))
  imp <- as.data.table(read_parquet(sprintf(MICE_F,imp_id)))
  a <- aug[order(stay_id,trial_num),.SD[1],by=stay_id]
  m <- imp[order(stay_id,trial_num),.SD[1],by=stay_id]
  stopifnot(identical(a$stay_id,m$stay_id))
  p <- a[,.(stay_id,subject_id,arm=as.integer(a$itt_group=="treatment"))]
  p <- merge(p,src,by="stay_id",sort=FALSE)
  # 90 天内死亡事件 + 死亡日；否则行政删失到 90
  p[, evt:=as.integer(!is.na(d90day) & d90day>=0 & d90day<=TAU)]
  p[, time:=ifelse(evt==1L,d90day,TAU)]
  for(v in c(CONT,BIN)) p[,(v):=m[[v]]]
  for(v in BIN) p[,(v):=as.integer(as.character(get(v)))]
  for(v in FACT) p[,(v):=factor(a[[v]],levels=LEV[[v]])]
  p
}

add_weight <- function(p){
  fml<-as.formula(paste("arm ~",paste(DAG,collapse=" + ")))
  ps<-predict(suppressWarnings(glm(fml,data=p,family=binomial(),
                                   control=list(maxit=100))),type="response")
  ps[is.na(ps)]<-mean(p$arm); pa<-mean(p$arm)
  w<-ifelse(p$arm==1,pa/ps,(1-pa)/(1-ps))
  p[,wgt:=pmin(pmax(w,quantile(w,.01)),quantile(w,.99))]
  p
}

# 加权风险 / RD / RR / OR
risk_contrast <- function(p){
  d<-p[,.(risk=sum(wgt*evt)/sum(wgt)),by=arm][order(-arm)]
  rt<-d[arm==1]$risk; rc<-d[arm==0]$risk
  c(risk_t=rt,risk_c=rc,rd=rt-rc,rr=rt/rc,or=(rt/(1-rt))/(rc/(1-rc)))
}

# 加权 Cox 总体 HR
cox_hr <- function(p){
  fit<-coxph(Surv(time,evt)~arm,data=p,weights=wgt,cluster=subject_id,
             robust=TRUE,control=coxph.control(iter.max=50))
  s<-summary(fit)
  c(beta=s$coefficients[1,"coef"],se=s$coefficients[1,"se(coef)"])
}

set.seed(SEED)
persons<-lapply(IMPS,load_persons)
cat("loaded m =",length(persons),"\n")

# 点估计
pt_risk<-rowMeans(sapply(persons,function(p)risk_contrast(add_weight(p))))
cx_l<-lapply(persons,function(p)cox_hr(add_weight(p)))
betas<-sapply(cx_l,function(z)z["beta"]); ses<-sapply(cx_l,function(z)z["se"])
rubin<-function(b,s){m<-length(b);bb<-mean(b);W<-mean(s^2)
  B<-if(m>1)sum((b-bb)^2/(m-1))else 0;sep<-sqrt(W+(1+1/m)*B);t<-bb/sep
  df<-if(m>1)(m-1)*(1+W/((1+1/m)*B))^2 else 1
  p<-2*pt(abs(t),df,lower.tail=FALSE)
  c(hr=exp(bb),ci_l=exp(bb-1.96*sep),ci_u=exp(bb+1.96*sep),p=p)}
hr_pt<-rubin(betas,ses)

cat("\n===== 点估计 加权风险（m=",length(IMPS),"） =====\n",sep="");print(round(pt_risk,4))
cat("\n===== 点估计 90d HR（Rubin） =====\n");print(round(hr_pt,4))

if(NBOOT>0){
  Bm<-matrix(NA_real_,NBOOT,6); colnames(Bm)<-c("risk_t","risk_c","rd","rr","or","beta")
  nb<-0L
  while(nb<NBOOT){
    p0<-persons[[sample(length(persons),1)]]
    pick<-sample(unique(p0$subject_id),uniqueN(p0$subject_id),replace=TRUE)
    pb<-p0[data.table(subject_id=pick),on="subject_id",allow.cartesian=TRUE]
    r<-tryCatch({pw<-add_weight(pb);list(rk=risk_contrast(pw),cx=cox_hr(pw))},
                error=function(e)NULL)
    if(is.null(r)||any(!is.finite(r$rk),!is.finite(r$cx[1])))next
    nb<-nb+1L
    Bm[nb,]<-c(r$rk,r$cx["beta"])
    if(nb%%50==0)cat(sprintf("  boot %d/%d\n",nb,NBOOT))
  }
  ci<-apply(Bm,2,quantile,c(.025,.975),na.rm=TRUE)
  # 风险/RD/RR/OR
  risk_final<-data.table(
    metric=colnames(Bm)[1:5],est=colMeans(Bm)[1:5],
    lci=ci[1,1:5],uci=ci[2,1:5])
  # HR
  hr_final<-data.table(analysis="90-day death (DAG14 stabilized IPTW)",
                       hr=exp(mean(Bm[,"beta"])),
                       ci_lower=exp(ci[1,"beta"]),ci_upper=exp(ci[2,"beta"]))
  cat("\n===== Bootstrap 95% CI 风险指标 =====\n"); print(risk_final)
  cat("\n===== Bootstrap 95% CI 90d HR =====\n"); print(hr_final)
  fwrite(risk_final,file.path(OUT_T,paste0("death90_risk_",TAG,".csv")))
  fwrite(hr_final,file.path(OUT_T,paste0("death90_hr_",TAG,".csv")))
  fwrite(as.data.table(Bm),file.path(OUT_T,paste0("death90_bootraw_",TAG,".csv")))
}else{
  fwrite(data.table(metric=names(pt_risk),est=pt_risk),
         file.path(OUT_T,paste0("death90_risk_",TAG,".csv")))
  fwrite(data.table(analysis="90-day death",t(hr_pt)),
         file.path(OUT_T,paste0("death90_hr_",TAG,".csv")))
}
cat("\ndone tag =",TAG,"\n")
