#!/usr/bin/env Rscript
# ============================================================================
# V28.4 次要结局重算（stay 唯一队列，DAG14 稳定化 IPTW）
#   结局（二分类）：
#     * AKI stage>=2 : aki_stage2_plus  (来源 outcome parquet；time_zero 起
#                       28天内 KDIGO stage 最大>=2，含 t0 基线已有 AKI)
#     * 输血需求     : transfusion_at_t0 (来源 MICE)
#   方法：DAG14 估 PS -> 稳定化 IPTW（1/99 截断）-> 加权组间风险、RD、RR；
#         m=20，B=200 MI-bootstrap（随机 imp × subject 重抽 × 重估 PS）
# 用法:
#   验证: Rscript secondary_outcomes.R --imps 1 --nboot 0 --tag imp1
#   全量: Rscript secondary_outcomes.R --imps 1..20 --nboot 200
# 须 LD_LIBRARY_PATH=/opt/miniforge3/lib   V28.4  2026-10-10
# ============================================================================
suppressPackageStartupMessages({library(data.table);library(arrow)})

BASE<-"/data/小米云盘/文档/论文/血磷"
AUG_F<-file.path(BASE,"TTE/v26/output/sequential_trials_augmented.parquet")
MICE_F<-file.path(BASE,"sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet")
OC_F<-file.path(BASE,"sql/result/M0_07_outcome_extraction.parquet")
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
OUTCOMES<-c(aki_stage2_plus="aki_stage2_plus", transfusion="transfusion_at_t0",
            aki_incident="aki_incident")
# 各结局是否限制到"基线无 AKI"人群（按列名=OUTCOMES 的值索引）
RESTRICT_NOBASE<-c(aki_stage2_plus=FALSE, transfusion_at_t0=FALSE, aki_incident=TRUE)
INC_FLAGS<-file.path(OUT_T,"aki_incident_flags.csv")

load_persons<-function(imp_id){
  aug<-as.data.table(read_parquet(AUG_F))
  imp<-as.data.table(read_parquet(sprintf(MICE_F,imp_id)))
  oc <-as.data.table(read_parquet(OC_F))
  fl <-fread(INC_FLAGS)
  a<-aug[order(stay_id,trial_num),.SD[1],by=stay_id]
  m<-imp[order(stay_id,trial_num),.SD[1],by=stay_id]
  o<-oc[order(stay_id,time_zero),.SD[1],by=stay_id]
  f<-fl[order(stay_id,time_zero),.SD[1],by=stay_id]
  stopifnot(identical(a$stay_id,m$stay_id),identical(a$stay_id,o$stay_id),
            identical(a$stay_id,f$stay_id))
  p<-a[,.(stay_id,subject_id)]
  p[,`:=`(arm=as.integer(as.character(m$treatment)),
          aki_stage2_plus=as.integer(o$aki_stage2_plus),
          transfusion_at_t0=as.integer(as.character(m$transfusion_at_t0)),
          aki_incident=as.integer(f$aki_incident), aki_pre=as.integer(f$aki_pre))]
  for(v in c(CONT,BIN))p[,(v):=m[[v]]]
  for(v in BIN)p[,(v):=as.integer(as.character(get(v)))]
  for(v in FACT)p[,(v):=factor(a[[v]],levels=LEV[[v]])]
  p
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

# 加权组间风险 + RD + RR（subset=逻辑向量或 NULL）
contrast<-function(p,yvar,subset=NULL){
  if(!is.null(subset)) p<-p[subset]
  d<-p[,.(risk=sum(wgt*get(yvar))/sum(wgt)),by=arm][order(-arm)]
  rt<-d[arm==1]$risk; rc<-d[arm==0]$risk
  c(risk_t=rt,risk_c=rc,rd=rt-rc,rr=rt/rc,or=(rt/(1-rt))/(rc/(1-rc)))
}

set.seed(SEED)
persons<-lapply(IMPS,load_persons)
cat("loaded m =",length(persons),"\n")

# 限制子集
sub_for<-function(p,yv) if(RESTRICT_NOBASE[[yv]]) (p$aki_pre==0L) else NULL

# 点估计（跨imp平均）
pt<-sapply(OUTCOMES,function(yv){
  rowMeans(sapply(persons,function(p)contrast(add_weight(p),yv,sub_for(p,yv))))
})
cat("\n===== 点估计（m=",length(IMPS),"） =====\n",sep="");print(round(pt,4))

# MI-bootstrap
if(NBOOT>0){
  boot_arr<-array(NA_real_,c(NBOOT,5,length(OUTCOMES)),
                  dimnames=list(NULL,c("risk_t","risk_c","rd","rr","or"),names(OUTCOMES)))
  nb<-0L
  while(nb<NBOOT){
    p0<-persons[[sample(length(persons),1)]]
    pick<-sample(unique(p0$subject_id),uniqueN(p0$subject_id),replace=TRUE)
    pb<-p0[data.table(subject_id=pick),on="subject_id",allow.cartesian=TRUE]
    pb<-add_weight(pb)
    ok<-TRUE
    for(j in seq_along(OUTCOMES)){
            r<-tryCatch(contrast(pb,OUTCOMES[j],sub_for(pb,OUTCOMES[j])),error=function(e)NULL)
      if(is.null(r)||any(!is.finite(r))){ok<-FALSE;break}
      boot_arr[nb+1,,j]<-r
    }
    if(ok){nb<-nb+1L;if(nb%%50==0)cat(sprintf("  boot %d/%d\n",nb,NBOOT))}
  }
  res_l<-list()
  for(j in seq_along(OUTCOMES)){
    m<-boot_arr[,,j]; ci<-apply(m,2,quantile,c(.025,.975),na.rm=TRUE)
    tab<-data.table(outcome=names(OUTCOMES)[j],metric=colnames(m),
                    est=colMeans(m),lci=ci[1,],uci=ci[2,])
    res_l[[j]]<-tab
  }
  final<-rbindlist(res_l)
  cat("\n===== Bootstrap 95% CI (B=",NBOOT,") =====\n",sep="")
  print(final)
  fwrite(final,file.path(OUT_T,paste0("secondary_",TAG,".csv")))
  fwrite(as.data.table(as.data.frame(boot_arr[,,,drop=TRUE])),
         file.path(OUT_T,paste0("secondary_bootraw_",TAG,".csv")))
}
cat("\ndone tag =",TAG,"\n")
