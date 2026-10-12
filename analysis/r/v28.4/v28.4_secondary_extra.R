#!/usr/bin/env Rscript
# ============================================================================
# V28.4 额外次要结局重算（stay 唯一队列，DAG14 稳定化 IPTW）
#   二分类（6）：
#     * ICU 死亡            death_icu              (M0_07)
#     * 院内死亡            death_hospital         (M0_07)
#     * 3 组分 CV 复合      cv_3comp               (M0_14d: 心律失常 OR
#                            血管加压升级 OR 心肌损伤)
#     * 新发心律失常        new_arrhythmia         (M0_14a)
#     * 心肌损伤            myocardial_injury      (M0_14c)
#     * 血管加压升级        vasopressor_escalation (M0_14b)
#   连续（2）：
#     * 无呼吸机日 VFD-28   vfd_28                 (M0_07)
#     * ICU 住院时长        icu_los                (M0_07, 自 t0 起)
#   口径与 v28.4_secondary_outcomes.R 完全一致：每 stay 取首行；DAG14 估 PS
#   -> 稳定化 IPTW（1/99 截断）；m=20，B=200 MI-bootstrap
#     （随机 imp × subject 重抽 × 重估 PS）。
# 用法:
#   验证: Rscript v28.4_secondary_extra.R --imps 1 --nboot 0 --tag imp1
#   全量: Rscript v28.4_secondary_extra.R --imps 1..20 --nboot 200
# 须 LD_LIBRARY_PATH=/opt/miniforge3/lib   V28.4  2026-10-10
# ============================================================================
suppressPackageStartupMessages({library(data.table);library(arrow)})

BASE<-"/data/小米云盘/文档/论文/血磷"
AUG_F<-file.path(BASE,"TTE/v26/output/sequential_trials_augmented.parquet")
MICE_F<-file.path(BASE,"sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet")
RRES<-file.path(BASE,"sql/result")
OUT_T<-file.path(BASE,"TTE/v28.4/output/tables/secondary")
dir.create(OUT_T,recursive=TRUE,showWarnings=FALSE)
SEED<-20260629L

argv<-commandArgs(trailingOnly=TRUE)
gv<-function(k,d){i<-match(k,argv);if(is.na(i))d else argv[i+1]}
IMPS<-as.integer(strsplit(gv("--imps","1"),",")[[1]])
NBOOT<-as.integer(gv("--nboot","0")); TAG<-gv("--tag","run")

BIN<-c("sepsis3","invasive_vent_at_t0","vasoactive_at_t0","refeeding_risk")
CONT<-c("age","baseline_sofa","charlson","egfr","baseline_phosphate",
        "lactate_t0","corrected_calcium_mmol_min")
FACT<-c("gender","admission_type","nutrition_support_type")
LEV<-list(gender=c("F","M"),
          admission_type=c("Observation","Emergency","Surgical_Elective"),
          nutrition_support_type=c("none","enteral","parenteral","mixed"))

# 结局列映射（数据列名 -> 分析名）
BIN_Y<-c(death_icu="death_icu",
         death_hospital="death_hospital",
         cv_3comp="cardiovascular_event_3comp",
         new_arrhythmia="new_arrhythmia_rhythm",
         myocardial_injury="myocardial_injury",
         vasopressor_escalation="vasopressor_escalation")
CONT_Y<-c(vfd_28="vfd_28", icu_los="icu_los_from_t0_days")

# ---- 只读一次的结局表（不随 imp 变），按 stay 去重 ----
read_y<-function(f,col){
  d<-as.data.table(read_parquet(file.path(RRES,f)))
  d<-d[order(stay_id,time_zero),.SD[1],by=stay_id]
  setnames(d,col,"yv"); d[,.(stay_id,yv=as.integer(yv))]
}
read_yc<-function(f,col){
  d<-as.data.table(read_parquet(file.path(RRES,f)))
  d<-d[order(stay_id,time_zero),.SD[1],by=stay_id]
  d[,.(stay_id,yv=get(col))]
}

# 二分类结局（M0_07 含 death_icu/death_hospital）
o07<-as.data.table(read_parquet(file.path(RRES,"M0_07_outcome_extraction.parquet")))
o07<-o07[order(stay_id,time_zero),.SD[1],by=stay_id]
yA<-o07[,.(stay_id, death_icu=as.integer(death_icu),
           death_hospital=as.integer(death_hospital))]
yB<-read_y("M0_14d_composite.parquet","cardiovascular_event_3comp")
  setnames(yB,"yv",BIN_Y[["cv_3comp"]])
yC<-read_y("M0_14a_arrhythmia.parquet","new_arrhythmia_rhythm")
  setnames(yC,"yv",BIN_Y[["new_arrhythmia"]])
yD<-read_y("M0_14c_myocardial.parquet","myocardial_injury")
  setnames(yD,"yv",BIN_Y[["myocardial_injury"]])
yE<-read_y("M0_14b_hemodynamic.parquet","vasopressor_escalation")
  setnames(yE,"yv",BIN_Y[["vasopressor_escalation"]])
# 连续结局
ycV<-read_yc("M0_07_outcome_extraction.parquet","vfd_28"); setnames(ycV,"yv",CONT_Y[["vfd_28"]])
ycL<-read_yc("M0_07_outcome_extraction.parquet","icu_los_from_t0_days"); setnames(ycL,"yv",CONT_Y[["icu_los"]])

merge_all<-function(p){
  p<-merge(p,yA,by="stay_id")
  p<-merge(p,yB,by="stay_id"); p<-merge(p,yC,by="stay_id")
  p<-merge(p,yD,by="stay_id"); p<-merge(p,yE,by="stay_id")
  p<-merge(p,ycV,by="stay_id"); p<-merge(p,ycL,by="stay_id")
  p
}

load_persons<-function(imp_id){
  aug<-as.data.table(read_parquet(AUG_F))
  imp<-as.data.table(read_parquet(sprintf(MICE_F,imp_id)))
  a<-aug[order(stay_id,trial_num),.SD[1],by=stay_id]
  m<-imp[order(stay_id,trial_num),.SD[1],by=stay_id]
  stopifnot(identical(a$stay_id,m$stay_id))
  p<-a[,.(stay_id,subject_id)]
  p[,arm:=as.integer(as.character(m$treatment))]
  for(v in c(CONT,BIN))p[,(v):=m[[v]]]
  for(v in BIN)p[,(v):=as.integer(as.character(get(v)))]
  for(v in FACT)p[,(v):=factor(a[[v]],levels=LEV[[v]])]
  merge_all(p)
}

# 稳定化权重
add_weight<-function(p){
  fml<-as.formula(paste("arm ~",paste(c(CONT,BIN,FACT),collapse=" + ")))
  ps<-predict(suppressWarnings(glm(fml,data=p,family=binomial(),
                                  control=list(maxit=100))),type="response")
  ps[is.na(ps)]<-mean(p$arm); pa<-mean(p$arm)
  w<-ifelse(p$arm==1,pa/ps,(1-pa)/(1-ps))
  p[,wgt:=pmin(pmax(w,quantile(w,.01)),quantile(w,.99))]
  p
}

# 二分类对比
contrast_bin<-function(p,yvar){
  d<-p[,.(risk=sum(wgt*get(yvar))/sum(wgt)),by=arm][order(-arm)]
  rt<-d[arm==1]$risk; rc<-d[arm==0]$risk
  c(risk_t=rt,risk_c=rc,rd=rt-rc,rr=rt/rc,or=(rt/(1-rt))/(rc/(1-rc)))
}
# 连续对比（加权均值差）
contrast_cont<-function(p,yvar){
  d<-p[,.(m=sum(wgt*get(yvar))/sum(wgt)),by=arm][order(-arm)]
  mt<-d[arm==1]$m; mc<-d[arm==0]$m
  c(mean_t=mt,mean_c=mc,md=mt-mc)
}

set.seed(SEED)
persons<-lapply(IMPS,load_persons)
cat("loaded m =",length(persons)," stays =",nrow(persons[[1]]),"\n")

# ---- 点估计（跨imp平均） ----
ptb<-sapply(BIN_Y,function(yv) rowMeans(sapply(persons,function(p)contrast_bin(add_weight(p),yv))))
ptc<-sapply(CONT_Y,function(yv) rowMeans(sapply(persons,function(p)contrast_cont(add_weight(p),yv))))
cat("\n===== 二分类点估计（m=",length(IMPS),"） =====\n",sep="");print(round(ptb,4))
cat("\n===== 连续点估计（m=",length(IMPS),"） =====\n",sep="");print(round(ptc,4))

# ---- MI-bootstrap ----
if(NBOOT>0){
  ab<-array(NA_real_,c(NBOOT,5,length(BIN_Y)),
            dimnames=list(NULL,c("risk_t","risk_c","rd","rr","or"),names(BIN_Y)))
  ac<-array(NA_real_,c(NBOOT,3,length(CONT_Y)),
            dimnames=list(NULL,c("mean_t","mean_c","md"),names(CONT_Y)))
  nb<-0L
  while(nb<NBOOT){
    p0<-persons[[sample(length(persons),1)]]
    pick<-sample(unique(p0$subject_id),uniqueN(p0$subject_id),replace=TRUE)
    pb<-p0[data.table(subject_id=pick),on="subject_id",allow.cartesian=TRUE]
    pb<-add_weight(pb)
    rb<-sapply(BIN_Y,function(yv)tryCatch(contrast_bin(pb,yv),error=function(e)rep(NA_real_,5)))
    rc_<-sapply(CONT_Y,function(yv)tryCatch(contrast_cont(pb,yv),error=function(e)rep(NA_real_,3)))
    if(all(is.finite(rb))&&all(is.finite(rc_))){
      ab[nb+1,,]<-rb; ac[nb+1,,]<-rc_; nb<-nb+1L
      if(nb%%50==0)cat(sprintf("  boot %d/%d\n",nb,NBOOT))
    }
  }
  build<-function(arr){
    rbindlist(lapply(dimnames(arr)[[3]],function(j){
      m<-arr[,,j]; ci<-apply(m,2,quantile,c(.025,.975),na.rm=TRUE)
      data.table(outcome=j,metric=colnames(m),est=colMeans(m),lci=ci[1,],uci=ci[2,])
    }))
  }
  fb<-build(ab); fc<-build(ac)
  cat("\n===== 二分类 Bootstrap 95% CI (B=",NBOOT,") =====\n",sep="");print(fb)
  cat("\n===== 连续 Bootstrap 95% CI (B=",NBOOT,") =====\n",sep="");print(fc)
  fwrite(fb,file.path(OUT_T,paste0("v28.4_extra_bin_",TAG,".csv")))
  fwrite(fc,file.path(OUT_T,paste0("v28.4_extra_cont_",TAG,".csv")))
}
cat("\ndone tag =",TAG,"\n")
