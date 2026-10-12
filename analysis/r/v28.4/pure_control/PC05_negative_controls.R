# =====================================================================
# PC05: Negative controls
#  NC1 sodium weighted GEE (A vs D) — expect ~null, clinically trivial
#  NC2 pre-T0 phosphate slope (pre24 -> pre6) weighted comparison
# =====================================================================
suppressMessages({library(DBI);library(RPostgres);library(dplyr);library(tidyr);library(geepack);library(readr)})
set.seed(20260629L)
OUT <- "/data/小米云盘/文档/论文/血磷/TTE/v28.4/output/pure_control"
con <- dbConnect(Postgres(),host="localhost",port=5432,dbname="mimiciv",user="postgres",password=Sys.getenv("MIMIC_PG_PASSWORD"))
co <- dbGetQuery(con,"SELECT * FROM public.pure_control_cohort")
na_ <- dbGetQuery(con,"SELECT stay_id,win,actual_hph,sodium FROM public.pure_control_sodium_windows")
pre<- dbGetQuery(con,"SELECT stay_id,win,actual_hph,phos FROM public.pure_control_prephos")
dbDisconnect(con)
co <- co %>% mutate(sexM=as.integer(gender=="M"),sepsis3=as.integer(sepsis3),
  vent_dim3=as.integer(vent_dim3),vasoactive_at_t0=as.integer(vasoactive_at_t0),
  refeeding_risk=as.integer(refeeding_risk))
DAG14 <- c("age","sexM","baseline_sofa","sepsis3","charlson","egfr","vent_dim3",
           "vasoactive_at_t0","admission_type","baseline_phosphate","lactate_t0",
           "refeeding_risk","corrected_calcium_mmol_min","nutrition_support_type")
wcent <- c(w6=6,w12=12,w24=24,w48=48,w72=72)

# ---------- NC1: sodium weighted GEE ----------
sl <- na_ %>% mutate(time=wcent[win]) %>% select(stay_id,time,sodium)
d0 <- co %>% filter(exposure_w48 %in% c("A_IV","D_none")) %>%
  select(stay_id,exposure_w48,sodium_t0,all_of(DAG14)) %>%
  filter(complete.cases(across(all_of(DAG14))))
grid <- expand.grid(stay_id=unique(d0$stay_id),time=c(6,12,24,48,72),stringsAsFactors=FALSE) %>%
  left_join(d0,by="stay_id") %>%
  left_join(sl,by=c("stay_id","time")) %>%
  mutate(measured=as.integer(!is.na(sodium)),A=as.integer(exposure_w48=="A_IV"),timef=factor(time))
ed <- grid %>% distinct(stay_id,.keep_all=TRUE)
ps <- glm(as.formula(paste("A~",paste(DAG14,collapse="+"))),data=ed,family=binomial)
ed$e<-predict(ps,type="response");p<-mean(ed$A)
ed$iptw<-ifelse(ed$A==1,p/ed$e,(1-p)/(1-ed$e))
mm <- glm(measured~A*timef+age+baseline_sofa+sepsis3+egfr+vent_dim3+vasoactive_at_t0+baseline_phosphate+charlson,data=grid,family=binomial)
mm0<-glm(measured~A*timef,data=grid,family=binomial)
grid$pm<-predict(mm,type="response");grid$pm0<-predict(mm0,type="response")
grid <- grid%>%left_join(ed%>%select(stay_id,iptw),by="stay_id")%>%mutate(w=iptw*pm0/pm)
wq<-quantile(grid$w,c(.01,.99),na.rm=TRUE);grid$wt<-pmin(pmax(grid$w,wq[1]),wq[2])
obs<-grid%>%filter(measured==1,!is.na(sodium_t0))
fit<-geeglm(sodium~factor(A)*timef+sodium_t0,id=stay_id,weights=wt,data=obs,corstr="exchangeable")
le<-expand.grid(A=c(1,0),timef=factor(c(6,12,24,48,72)),sodium_t0=mean(obs$sodium_t0,na.rm=TRUE))
X<-model.matrix(~factor(A)*timef+sodium_t0,data=le);eta<-as.numeric(X%*%coef(fit));se<-sqrt(diag(X%*%vcov(fit)%*%t(X)))
nc1<-data.frame(control="NC1_sodium",time=as.numeric(as.character(le$timef)),A=le$A,est=eta,se=se)%>%
  pivot_wider(names_from=A,values_from=c(est,se),names_sep="")%>%
  mutate(diff=est1-est0,diff_se=sqrt(se1^2+se0^2),lcl=diff-1.96*diff_se,ucl=diff+1.96*diff_se)%>%
  select(control,time,diff,lcl,ucl)
write.csv(nc1,file.path(OUT,"pc_nc1_sodium.csv"),row.names=FALSE)

# ---------- NC2: pre-T0 phosphate slope ----------
pw <- pre %>% select(stay_id,win,actual_hph,phos) %>% pivot_wider(names_from=win,values_from=c(phos,actual_hph))
pw <- pw %>% filter(!is.na(phos_pre24),!is.na(phos_pre6)) %>%
  mutate(slope=(phos_pre6-phos_pre24)/(actual_hph_pre6-actual_hph_pre24),
         drop=phos_pre6-phos_pre24)
pd <- pw %>% inner_join(co,by="stay_id") %>% filter(exposure_w48 %in% c("A_IV","D_none")) %>%
  filter(complete.cases(across(all_of(DAG14)))) %>% mutate(A=as.integer(exposure_w48=="A_IV"))
psf<-glm(as.formula(paste("A~",paste(DAG14,collapse="+"))),data=pd,family=binomial)
pd$e<-predict(psf,type="response");pp<-mean(pd$A)
pd$w<-ifelse(pd$A==1,pp/pd$e,(1-pp)/(1-pd$e))
wm <- function(dd,col){dd%>%group_by(A)%>%summarise(m=weighted.mean(.data[[col]],w),.groups="drop")}
g_sl<-wm(pd,"slope");g_dr<-wm(pd,"drop");g_p24<-wm(pd,"phos_pre24");g_p6<-wm(pd,"phos_pre6")
nc2<-data.frame(
  metric=c("slope_mmol_per_h","drop_pre24_to_pre6","level_pre24","level_pre6"),
  A_IV=c(g_sl$m[g_sl$A==1],g_dr$m[g_dr$A==1],g_p24$m[g_p24$A==1],g_p6$m[g_p6$A==1]),
  D_none=c(g_sl$m[g_sl$A==0],g_dr$m[g_dr$A==0],g_p24$m[g_p24$A==0],g_p6$m[g_p6$A==0]),
  n_IV=nrow(pd%>%filter(A==1)),n_none=nrow(pd%>%filter(A==0))) %>%
  mutate(diff=round(A_IV-D_none,4))
write.csv(nc2,file.path(OUT,"pc_nc2_pretrend.csv"),row.names=FALSE)

cat("=== NC1 SODIUM weighted GEE contrast (IV - none), mmol/L ===\n")
write.table(as.data.frame(nc1),quote=FALSE,row.names=FALSE,sep="\t")
cat("\n=== NC2 PRE-T0 PHOSPHATE TREND ===\n")
write.table(as.data.frame(nc2),quote=FALSE,row.names=FALSE,sep="\t")
