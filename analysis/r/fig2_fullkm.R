#!/usr/bin/env Rscript
# ============================================================================
# V28.4 Figure 2 — 完整 28 天时间线上的稳定化加权 KM（推荐方案）
#   * 每 stay 一条完整 day0..28 时间线，不在出院处截断（结局对全员已知）
#   * 精确死亡 2552 人落真实死亡日；仅知死亡/不知日期 1076 人在已知区间
#     [出院日, 28] 内做区间落点（多重落点 + MI-bootstrap 平均）
#   * DAG14 稳定化治疗权重（1/99 截断）；因全程随访，无删失，加权 KM 即
#     加权累计风险，日28 自然落在 stay 层真实死亡水平（~16%）
#   * B=200 MI-bootstrap（subject 重抽 × 区间重落点 × 重估 PS）
# 用法:
#   验证: Rscript fig2_fullkm.R --imps 1 --nboot 0 --tag imp1
#   全量: Rscript fig2_fullkm.R --imps 1..20 --nboot 200
# 须 LD_LIBRARY_PATH=/opt/miniforge3/lib   V28.4  2026-10-10
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(arrow)
  library(ggplot2); library(scales); library(cowplot)
})

BASE <- "/data/小米云盘/文档/论文/血磷"
AUG_F <- file.path(BASE,"TTE/v26/output/sequential_trials_augmented.parquet")
MICE_F <- file.path(BASE,"sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet")
OUT_T <- file.path(BASE,"TTE/v28.4/output/tables/fig2")
OUT_F <- file.path(BASE,"TTE/v28.4/output/figures")
dir.create(OUT_T, recursive=TRUE, showWarnings=FALSE)
dir.create(OUT_F, recursive=TRUE, showWarnings=FALSE)
TAU<-28L; SEED<-20260629L

argv <- commandArgs(trailingOnly=TRUE)
gv <- function(k,d){i<-match(k,argv);if(is.na(i))d else argv[i+1]}
IMPS <- as.integer(strsplit(gv("--imps","1"),",")[[1]])
NBOOT <- as.integer(gv("--nboot","0")); TAG <- gv("--tag","run")
N_PLACE <- 20L
RISK_TIMES <- c(0L,7L,14L,21L,28L)

BIN <- c("sepsis3","invasive_vent_at_t0","vasoactive_at_t0","refeeding_risk")
CONT<- c("age","baseline_sofa","charlson","egfr","baseline_phosphate",
         "lactate_t0","corrected_calcium_mmol_min")
FACT<- c("gender","admission_type","nutrition_support_type")
LEV <- list(gender=c("F","M"),
            admission_type=c("Observation","Emergency","Surgical_Elective"),
            nutrition_support_type=c("none","enteral","parenteral","mixed"))

# ---- person 表（stay 唯一，含死亡日落点所需信息） ----
load_persons <- function(imp_id){
  aug <- as.data.table(read_parquet(AUG_F))
  imp <- as.data.table(read_parquet(sprintf(MICE_F,imp_id)))
  a <- aug[order(stay_id,trial_num),.SD[1],by=stay_id]
  m <- imp[order(stay_id,trial_num),.SD[1],by=stay_id]
  stopifnot(identical(a$stay_id,m$stay_id))
  p <- a[,.(stay_id,subject_id,death=as.integer(death_28day_num))]
  p[, `:=`(arm=as.integer(as.character(m$treatment)),
           exact_t=as.numeric(a$death_time_from_t0),
           disch=as.numeric(a$dischtime_from_t0))]
  for(v in c(CONT,BIN)) p[,(v):=m[[v]]]
  for(v in BIN) p[,(v):=as.integer(as.character(get(v)))]
  for(v in FACT) p[,(v):=factor(a[[v]],levels=LEV[[v]])]
  # 精确死亡日
  p[, has_exact := death==1L & !is.na(exact_t)]
  p[, death_day := NA_integer_]
  p[has_exact==TRUE, death_day := pmin(TAU,pmax(0L,as.integer(ceiling(exact_t))))]
  # 区间下界（出院日；缺失用0）
  p[, dlower := as.integer(pmin(TAU,pmax(0L,ceiling(disch))))]
  p[is.na(dlower), dlower := 0L]
  p
}

# ---- 区间死亡随机落点 ----
assign_placement <- function(p){
  pp <- copy(p)
  idx <- which(pp$death==1L & !pp$has_exact)
  if(length(idx)){
    lo <- pp$dlower[idx]
    pp$death_day[idx] <- lo + as.integer(floor(runif(length(idx))*(TAU-lo+1L)))
  }
  # 存活者（及任何未赋值）在整个0..28窗口内存活：落到 TAU+1
  pp$death_day[is.na(pp$death_day)] <- TAU+1L
  pp
}

# ---- DAG14 稳定化治疗权重（1/99 截断） ----
add_ps_weight <- function(p){
  fml <- as.formula(paste("arm ~",paste(c(CONT,BIN,FACT),collapse=" + ")))
  ps <- predict(suppressWarnings(glm(fml,data=p,family=binomial(),
                                     control=list(maxit=100))),type="response")
  ps[is.na(ps)] <- mean(p$arm)
  pa <- mean(p$arm)
  w <- ifelse(p$arm==1, pa/ps, (1-pa)/(1-ps))
  p[, wgt := pmin(pmax(w,quantile(w,.01)),quantile(w,.99))]
  p
}

# ---- 给定落点+权重，算日网格加权累计风险（每臂） ----
# 完全随访：t 日累计风险 = 加权 I(death_day<=t)/总权重
fullkm <- function(p){
  rbindlist(lapply(c(1L,0L), function(aa){
    d <- p[arm==aa]
    tot <- sum(d$wgt)
    data.table(arm=aa, day=0:TAU,
               cif=sapply(0:TAU, function(t) sum(d$wgt[d$death_day<=t])/tot))
  }))
}
risk_tab <- function(p){
  rbindlist(lapply(c(1L,0L), function(aa){
    d <- p[arm==aa]
    data.table(arm=aa, day=RISK_TIMES,
               nrisk=sapply(RISK_TIMES, function(t) sum(d$wgt[d$death_day>=t])))
  }))
}

# ============ 点估计（跨 imp × 多重落点 平均） ============
set.seed(SEED)
persons <- lapply(IMPS, load_persons)
cat("person 表加载完成 m =",length(persons),"\n")

curve_l <- list(); risk_l <- list()
for(p0 in persons){
  pw <- add_ps_weight(p0)
  for(s in seq_len(N_PLACE)){
    pp <- assign_placement(pw)
    curve_l[[length(curve_l)+1]] <- fullkm(pp)
  }
  risk_l[[length(risk_l)+1]] <- risk_tab(assign_placement(pw))
}
pt_curve <- rbindlist(curve_l)[,.(cif=mean(cif)),by=.(arm,day)]
pt_risk  <- rbindlist(risk_l)[,.(nrisk=mean(nrisk)),by=.(arm,day)]

cat("\n===== 点估计：日28 加权累计风险 =====\n")
print(pt_curve[day==TAU])
pt_curve[, S:=1-cif]
cat("RMST:"); print(pt_curve[,.(rmst=sum(S[day<TAU])),by=arm])

# ============ MI-bootstrap CI ============
if(NBOOT>0){
  nb <- 0L; Bmat <- matrix(NA_real_,NBOOT,8); Carray<-array(NA_real_,c(NBOOT,2,TAU+1L))
  colnames(Bmat)<-c("risk_t","risk_c","rd","rmst_t","rmst_c","rmst_d","hr_anchor","x")
  while(nb<NBOOT){
    p0 <- persons[[sample(length(persons),1)]]
    pick <- sample(unique(p0$subject_id),uniqueN(p0$subject_id),replace=TRUE)
    pb <- p0[data.table(subject_id=pick),on="subject_id",allow.cartesian=TRUE]
    pb <- assign_placement(add_ps_weight(pb))
    cv <- fullkm(pb)
    g_t <- cv[arm==1 & day==TAU]$cif; g_c<-cv[arm==0 & day==TAU]$cif
    cv[, S:=1-cif]
    r_t <- cv[arm==1,sum(S[day<TAU])]; r_c<-cv[arm==0,sum(S[day<TAU])]
    nb<-nb+1L
    Bmat[nb,]<-c(g_t,g_c,g_t-g_c,r_t,r_c,r_t-r_c, g_t/g_c, NA_real_)
    Carray[nb,1,]<-cv[arm==1][order(day)]$cif
    Carray[nb,2,]<-cv[arm==0][order(day)]$cif
    if(nb%%50==0) cat(sprintf("  boot %d/%d\n",nb,NBOOT))
  }
  ci<-apply(Bmat,2,quantile,c(.025,.975),na.rm=TRUE)
  cat("\n===== Bootstrap 95% CI =====\n")
  print(round(t(rbind(est=colMeans(Bmat),lci=ci[1,],uci=ci[2,])),4))
  # 曲线 CI
  cc <- rbindlist(lapply(c(1L,0L),function(aa){
    ai<-if(aa==1)1 else 2
    lo<-apply(Carray[,ai,],2,quantile,.025); hi<-apply(Carray[,ai,],2,quantile,.975)
    data.table(arm=aa,day=0:TAU,cif_lo=lo,cif_hi=hi)
  }))
  pt_curve <- merge(pt_curve,cc,by=c("arm","day"))
  fwrite(as.data.table(Bmat),file.path(OUT_T,paste0("fig2_fullkm_boot_",TAG,".csv")))
} else {
  pt_curve[, `:=`(cif_lo=cif,cif_hi=cif)]
}

fwrite(pt_curve,file.path(OUT_T,paste0("fig2_fullkm_curve_",TAG,".csv")))
fwrite(pt_risk,file.path(OUT_T,paste0("fig2_fullkm_risk_",TAG,".csv")))

# ================= 绘图 =================
FONT<-"Liberation Sans"; COL_T<-"#BC3C29"; COL_C<-"#0072B5"
pt_curve[, grp:=factor(arm,c(1,0),labels=c("IV phosphate","Conservative management"))]
pt_risk[, grp:=factor(arm,c(1,0),labels=c("IV phosphate","Conservative management"))]
pt_risk[, `:=`(label=comma(round(nrisk)), yn=ifelse(arm==1,2,1))]
pr<-pt_risk[day%in%RISK_TIMES]

theme_fig<-theme_classic(base_size=9,base_family=FONT)+
  theme(axis.line=element_line(linewidth=.4,color="black"),axis.ticks=element_line(linewidth=.4,color="black"),
        axis.title=element_text(size=9.5,face="bold"),axis.text=element_text(size=9,color="black"),
        plot.title=element_text(face="bold",size=10.5,hjust=0),legend.position=c(.99,.25),
        legend.justification=c("right","top"),legend.background=element_rect("white",color=NA),
        legend.title=element_text(size=8.5,face="bold"),legend.text=element_text(size=8.5),
        legend.key=element_blank(),legend.key.size=unit(.35,"cm"),
        plot.margin=margin(4,8,0,4),panel.grid=element_blank())

p_main<-ggplot(pt_curve,aes(day))+
  geom_ribbon(aes(ymin=cif_lo,ymax=cif_hi,fill=grp),alpha=.22,linewidth=0)+
  geom_step(aes(y=cif,colour=grp),linewidth=.8)+
  scale_colour_manual(values=c("IV phosphate"=COL_T,"Conservative management"=COL_C),name="Treatment strategy")+
  scale_fill_manual(values=c("IV phosphate"=COL_T,"Conservative management"=COL_C),guide="none")+
  coord_cartesian(c(0,TAU),c(0,.2))+
  scale_x_continuous(breaks=seq(0,TAU,7),expand=expansion(mult=.04))+
  scale_y_continuous(breaks=seq(0,.2,.05),labels=function(x)sprintf("%.2f",x))+
  labs(title="Adjusted Cumulative Incidence of 28-Day Death (full-horizon weighted KM)",
       x="Follow-up (days)",y="Cumulative incidence of death")+theme_fig

p_riskp<-ggplot(pr,aes(day,yn))+
  geom_text(aes(label=label),size=2.6,family=FONT)+
  scale_x_continuous(limits=c(0,TAU),breaks=seq(0,TAU,7),expand=expansion(mult=.04))+
  scale_y_continuous(limits=c(.6,2.4),breaks=c(1,2),labels=c("Conservative","IV phosphate"))+
  labs(title="Number at risk (weighted)")+
  theme_minimal(base_family=FONT)+
  theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.title=element_blank(),
        axis.text.y=element_text(size=8.5,face="bold",color="black"),
        plot.title=element_text(size=8.5,hjust=0,face="bold"),plot.margin=margin(10,8,2,4))

p_final<-suppressWarnings(plot_grid(p_main,p_riskp,ncol=1,rel_heights=c(.82,.18),align="v",axis="lr"))
ggsave(file.path(OUT_F,"Fig2_FullHorizon_Weighted_KM.png"),p_final,width=7.1,height=6.3,dpi=600,bg="white")
ggsave(file.path(OUT_F,"Fig2_FullHorizon_Weighted_KM.pdf"),p_final,width=7.1,height=6.3,device=cairo_pdf,bg="white")
ggsave(file.path(OUT_F,"Fig2_FullHorizon_Weighted_KM.svg"),p_final,width=7.1,height=6.3,bg="white")
cat("\nFigure 2 已保存至",OUT_F,"\n")
