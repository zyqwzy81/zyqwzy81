# =====================================================================
# PC12: Clinical-outcome uncertainty addendum (review M9)
# One B=200 subject-level MI-free bootstrap that re-estimates the PS
# every replicate and returns:
#   CCW HR + day-28 RD (clones rebuilt per replicate)
#   IPTW Cox HR, Adjusted Cox HR, 48-h landmark HR
#   AIPW RD + RR; RMST28 difference
#   48-h transfusion RD/RR; 7-d AKI>=2 RD/RR
#   nutrition=none subset HR (subset PS refit per replicate)
# Point estimates reproduce pc_clinical_main.csv / pc_clinical_secondary.csv.
# POST-HOC / hypothesis-generating. DB password: MIMIC_PG_PASSWORD env.
# =====================================================================
suppressMessages({library(DBI);library(RPostgres);library(dplyr);library(tidyr)
  library(survival);library(arrow);library(readr);library(data.table)})
set.seed(20260629L)
B <- 200L
OUT <- "/data/小米云盘/文档/论文/血磷/TTE/v28.4/output/pure_control"
BASE <- "/data/小米云盘/文档/论文/血磷"
con <- dbConnect(Postgres(),host="localhost",port=5432,dbname="mimiciv",
                 user="postgres",password=Sys.getenv("MIMIC_PG_PASSWORD"))
co <- dbGetQuery(con,"SELECT * FROM public.pure_control_cohort")
mor<- dbGetQuery(con,"SELECT stay_id,death28,surv_days FROM public.pure_control_mortality")
dbDisconnect(con)
co <- co %>% mutate(sexM=as.integer(gender=="M"),sepsis3=as.integer(sepsis3),
  vent_dim3=as.integer(vent_dim3),vasoactive_at_t0=as.integer(vasoactive_at_t0),
  refeeding_risk=as.integer(refeeding_risk))
DAG14 <- c("age","sexM","baseline_sofa","sepsis3","charlson","egfr","vent_dim3",
           "vasoactive_at_t0","admission_type","baseline_phosphate","lactate_t0",
           "refeeding_risk","corrected_calcium_mmol_min","nutrition_support_type")
prbc <- read_parquet(file.path(BASE,"TTE/v28.4/output/extraction/prbc_events.parquet"))
kd   <- read_parquet(file.path(BASE,"TTE/v28.4/output/extraction/kdigo_stages_ts.parquet"))
t0v <- co %>% select(stay_id,t0)
prbc_ev <- prbc %>% inner_join(t0v,by="stay_id") %>%
  mutate(hh=as.numeric(difftime(starttime,t0,units="hours"))) %>%
  filter(hh>0,hh<=48) %>% distinct(stay_id) %>% mutate(prbc48=1L)
aki_ev <- kd %>% filter(!is.na(aki_stage),aki_stage>=2) %>% inner_join(t0v,by="stay_id") %>%
  mutate(dd=as.numeric(difftime(charttime,t0,units="days"))) %>%
  filter(dd>0,dd<=7) %>% distinct(stay_id) %>% mutate(aki7=1L)
an <- co %>% filter(exposure_w48 %in% c("A_IV","D_none")) %>%
  inner_join(mor,by="stay_id") %>% left_join(prbc_ev,by="stay_id") %>%
  left_join(aki_ev,by="stay_id") %>%
  mutate(prbc48=tidyr::replace_na(prbc48,0L),aki7=tidyr::replace_na(aki7,0L),
         A=as.integer(exposure_w48=="A_IV"),
         iv_day=as.numeric(difftime(iv_first_w48,t0,units="days")),
         admission_type=ifelse(grepl("ELECTIVE|SURGICAL|AMBULATORY",admission_type),
                               "elective","emergency"),
         nutrition_support_type=ifelse(nutrition_support_type=="none","none","any"))
cc <- as.data.table(an %>% filter(complete.cases(across(all_of(DAG14)))))
psf <- glm(as.formula(paste("A~",paste(DAG14,collapse="+"))),data=cc,family=binomial)
cc$e <- predict(psf,type="response"); pA<-mean(cc$A)
cc$w <- ifelse(cc$A==1,pA/cc$e,(1-pA)/(1-cc$e))
wq <- quantile(cc$w,c(.01,.99)); cc$wt <- pmin(pmax(cc$w,wq[1]),wq[2])
DAGnn <- setdiff(DAG14,"nutrition_support_type")

build_ccw <- function(d){
  grid <- d[, .(day=seq(0,floor(pmin(surv_days[1],28)))), by=stay_id]
  pd <- grid[, .(z=c(1L,0L)), by=.(stay_id,day)]
  setkey(pd,stay_id)
  meta <- d[, .(stay_id,A,death28,surv_days,iv_day,
                age,sexM,baseline_sofa,sepsis3,charlson,egfr,vent_dim3,
                vasoactive_at_t0,admission_type,baseline_phosphate,lactate_t0,
                refeeding_risk,corrected_calcium_mmol_min,nutrition_support_type)]
  pd <- merge(pd,meta,by="stay_id",allow.cartesian=TRUE)
  pd[, `:=`(iv_day=ifelse(is.na(iv_day),99,iv_day),
    cens=fifelse(z==1 & A==0 & day>=2,1L,
          fifelse(z==0 & A==1 & iv_day<=day,1L,0L)))]
  pd[, firstc := min(c(day[cens==1],29)), by=.(stay_id,z)]
  pd[, evday := as.integer(death28==1 & day==floor(surv_days) & day<=firstc)]
  fit0 <- glm(cens~z+poly(day,2),data=pd,family=binomial)
  fit1 <- glm(as.formula(paste("cens~z*poly(day,2)+",paste(DAG14,collapse="+"))),data=pd,family=binomial)
  pd[, `:=`(p0=pmin(pmax(predict(fit0,pd,type="response"),1e-4),1-1e-4),
            p1=pmin(pmax(predict(fit1,pd,type="response"),1e-4),0.99))]
  pd <- pd[day<firstc | (day==firstc & cens==0)]
  pd[, sw := cumprod((1-p0)/(1-p1)), by=.(stay_id,z)]
  pd[, sw:=sw[.N], by=.(stay_id,z)]
  pd[, .(sw=sw[1],
         time=fifelse(firstc[1]<29,firstc[1],pmin(surv_days[1],28)),
         event=fifelse(firstc[1]<29,0L,as.integer(death28[1]==1 & surv_days[1]<=28)),
         z=z[1],censored=as.integer(firstc[1]<29)), by=.(stay_id,z)]
}
rmst28 <- function(d,g,wcol="wt"){
  dd<-d[A==g]; sf<-survfit(Surv(surv_days,death28)~1,data=dd,weights=dd[[wcol]])
  t<-c(0,sf$time);S<-c(1,sf$surv);k<-t<=28;sum(diff(c(t[k],28))*S[k])
}
risk28 <- function(cl0,wcol="swb"){
  s1<-survfit(Surv(time,event)~1,data=cl0[z==1],weights=cl0[z==1][[wcol]])
  s0<-survfit(Surv(time,event)~1,data=cl0[z==0],weights=cl0[z==0][[wcol]])
  c(r1=1-summary(s1,times=28,extend=TRUE)$surv, r0=1-summary(s0,times=28,extend=TRUE)$surv)
}
gethr <- function(fit,term="A"){
  if(inherits(fit,"try-error")||is.null(fit)) return(c(NA,NA,NA))
  c(unname(exp(coef(fit))[term]), as.numeric(exp(confint(fit))[term,]))
}

# ---- point estimates (original data) ----
cl <- build_ccw(cc)
cwq <- quantile(cl$sw,c(.01,.99)); cl[,swb:=pmin(pmax(sw,cwq[1]),cwq[2])]
ccw_p <- gethr(coxph(Surv(time,event)~z,data=cl,weights=swb,robust=TRUE),"z")
ccw_risk_p <- risk28(cl)
iptw_p <- gethr(coxph(Surv(surv_days,death28)~A,data=cc,weights=wt,robust=TRUE))
adj_p  <- gethr(coxph(as.formula(paste("Surv(surv_days,death28)~A+",
                    paste(DAG14,collapse="+"))),data=cc,weights=wt,robust=TRUE))
lmk0 <- cc[surv_days>2]; lmk0[,t2:=surv_days-2]
lmk_p <- gethr(coxph(Surv(t2,death28)~A,data=lmk0,weights=wt,robust=TRUE))
m1<-glm(as.formula(paste("death28~",paste(DAG14,collapse="+"))),data=cc[A==1],family=binomial)
m0<-glm(as.formula(paste("death28~",paste(DAG14,collapse="+"))),data=cc[A==0],family=binomial)
cc$q1<-predict(m1,cc,type="response");cc$q0<-predict(m0,cc,type="response")
cc$ip1<-cc$q1+cc$A/cc$e*(cc$death28-cc$q1);cc$ip0<-cc$q0+(1-cc$A)/(1-cc$e)*(cc$death28-cc$q0)
r1p<-weighted.mean(cc$ip1,cc$wt);r0p<-weighted.mean(cc$ip0,cc$wt)
rm_p <- c(rmst28(cc,1),rmst28(cc,0))
binp <- function(col){y<-cc[[col]];x1<-weighted.mean(y[cc$A==1],cc$wt[cc$A==1]);
  x0<-weighted.mean(y[cc$A==0],cc$wt[cc$A==0]);c(x1,x0,x1-x0,x1/x0)}
pr_p<-binp("prbc48");ak_p<-binp("aki7")
nn0<-cc[nutrition_support_type=="none"]
fnn<-glm(as.formula(paste("A~",paste(DAGnn,collapse="+"))),data=nn0,family=binomial)
nn0$en<-predict(fnn,type="response");pn<-mean(nn0$A)
nn0$wn<-pmin(pmax(ifelse(nn0$A==1,pn/nn0$en,(1-pn)/(1-nn0$en)),wq[1]),wq[2])
nn_p <- c(nrow(nn0), gethr(coxph(Surv(surv_days,death28)~A,data=nn0,weights=wn,robust=TRUE)))

# ---- bootstrap ----
n<-nrow(cc)
M <- c("CCW_HR","CCW_RD28","IPTW_HR","Adjusted_HR","Landmark48_HR",
       "AIPW_RD","AIPW_RR","RMST28_diff",
       "PRBC48_RD","PRBC48_RR","AKI7_RD","AKI7_RR","NutrNone_HR")
R<-matrix(NA_real_,B,length(M));colnames(R)<-M
t0<-Sys.time()
for(b in seq_len(B)){
  db<-cc[sample(seq_len(n),n,replace=TRUE)]
  fb<-try(glm(as.formula(paste("A~",paste(DAG14,collapse="+"))),data=db,family=binomial),silent=TRUE)
  if(inherits(fb,"try-error"))next
  eb<-predict(fb,type="response");pb<-mean(db$A)
  db$wb<-pmin(pmax(ifelse(db$A==1,pb/eb,(1-pb)/(1-eb)),wq[1]),wq[2])
  # CCW (clones rebuilt; original full-data sw trim points)
  clb<-try(build_ccw(db),silent=TRUE)
  if(!inherits(clb,"try-error")){
    clb[,swb:=pmin(pmax(sw,cwq[1]),cwq[2])]
    hb<-try(coxph(Surv(time,event)~z,data=clb,weights=swb,robust=TRUE),silent=TRUE)
    if(!inherits(hb,"try-error")) R[b,"CCW_HR"]<-unname(exp(coef(hb))["z"])
    rr<-try(risk28(clb),silent=TRUE)
    if(!inherits(rr,"try-error")) R[b,"CCW_RD28"]<-rr["r1"]-rr["r0"]
  }
  h1<-try(coxph(Surv(surv_days,death28)~A,data=db,weights=wb,robust=TRUE),silent=TRUE)
  if(!inherits(h1,"try-error"))R[b,"IPTW_HR"]<-unname(exp(coef(h1))["A"])
  h2<-try(coxph(as.formula(paste("Surv(surv_days,death28)~A+",paste(DAG14,collapse="+"))),
                data=db,weights=wb,robust=TRUE),silent=TRUE)
  if(!inherits(h2,"try-error"))R[b,"Adjusted_HR"]<-unname(exp(coef(h2))["A"])
  lb<-db[surv_days>2];lb[,t2:=surv_days-2]
  h3<-try(coxph(Surv(t2,death28)~A,data=lb,weights=wb,robust=TRUE),silent=TRUE)
  if(!inherits(h3,"try-error"))R[b,"Landmark48_HR"]<-unname(exp(coef(h3))["A"])
  g1<-try(glm(as.formula(paste("death28~",paste(DAG14,collapse="+"))),data=db[A==1],family=binomial),silent=TRUE)
  g0<-try(glm(as.formula(paste("death28~",paste(DAG14,collapse="+"))),data=db[A==0],family=binomial),silent=TRUE)
  if(!inherits(g1,"try-error")&&!inherits(g0,"try-error")){
    q1<-predict(g1,db,type="response");q0<-predict(g0,db,type="response")
    ip1<-q1+db$A/eb*(db$death28-q1);ip0<-q0+(1-db$A)/(1-eb)*(db$death28-q0)
    x1<-weighted.mean(ip1,db$wb);x0<-weighted.mean(ip0,db$wb)
    R[b,c("AIPW_RD","AIPW_RR")]<-c(x1-x0,x1/x0)}
  R[b,"RMST28_diff"]<-try(rmst28(db,1,"wb")-rmst28(db,0,"wb"),silent=TRUE)
  for(v in c("prbc48","aki7")){
    y<-db[[v]];x1<-weighted.mean(y[db$A==1],db$wb[db$A==1]);x0<-weighted.mean(y[db$A==0],db$wb[db$A==0])
    R[b,paste0(ifelse(v=="prbc48","PRBC48","AKI7"),"_RD")]<-x1-x0
    R[b,paste0(ifelse(v=="prbc48","PRBC48","AKI7"),"_RR")]<-x1/x0}
  nb<-db[nutrition_support_type=="none"]
  if(nrow(nb)>=100 && length(unique(nb$A))==2){
    fn<-try(glm(as.formula(paste("A~",paste(DAGnn,collapse="+"))),data=nb,family=binomial),silent=TRUE)
    if(!inherits(fn,"try-error")){
      en<-predict(fn,type="response");pnb<-mean(nb$A)
      nb$wn<-pmin(pmax(ifelse(nb$A==1,pnb/en,(1-pnb)/(1-en)),wq[1]),wq[2])
      hn<-try(coxph(Surv(surv_days,death28)~A,data=nb,weights=wn,robust=TRUE),silent=TRUE)
      if(!inherits(hn,"try-error"))R[b,"NutrNone_HR"]<-unname(exp(coef(hn))["A"])}}
  if(b%%20==0)cat("b=",b," elapsed=",round(as.numeric(Sys.time()-t0,units="mins"),1),"min\n")
}
qci <- function(x) quantile(x,c(.025,.975),na.rm=TRUE)
boot_ci <- t(apply(R,2,qci))
out <- data.frame(
  metric=c("CCW_HR","CCW_RD28_risk1","CCW_RD28_risk0","CCW_RD28",
           "IPTW_HR","Adjusted_HR","Landmark48_HR",
           "AIPW_risk1","AIPW_risk0","AIPW_RD","AIPW_RR",
           "RMST28_IV","RMST28_none","RMST28_diff",
           "PRBC48_rate1","PRBC48_rate0","PRBC48_RD","PRBC48_RR",
           "AKI7_rate1","AKI7_rate0","AKI7_RD","AKI7_RR",
           "NutrNone_HR"),
  est=c(ccw_p[1],ccw_risk_p["r1"],ccw_risk_p["r0"],ccw_risk_p["r1"]-ccw_risk_p["r0"],
        iptw_p[1],adj_p[1],lmk_p[1],r1p,r0p,r1p-r0p,r1p/r0p,
        rm_p[1],rm_p[2],rm_p[1]-rm_p[2],
        pr_p[1],pr_p[2],pr_p[3],pr_p[4],ak_p[1],ak_p[2],ak_p[3],ak_p[4],nn_p[2]),
  boot_lcl=c(boot_ci["CCW_HR",1],NA,NA,boot_ci["CCW_RD28",1],
        boot_ci["IPTW_HR",1],boot_ci["Adjusted_HR",1],boot_ci["Landmark48_HR",1],
        NA,NA,boot_ci["AIPW_RD",1],boot_ci["AIPW_RR",1],
        NA,NA,boot_ci["RMST28_diff",1],
        NA,NA,boot_ci["PRBC48_RD",1],boot_ci["PRBC48_RR",1],
        NA,NA,boot_ci["AKI7_RD",1],boot_ci["AKI7_RR",1],boot_ci["NutrNone_HR",1]),
  boot_ucl=c(boot_ci["CCW_HR",2],NA,NA,boot_ci["CCW_RD28",2],
        boot_ci["IPTW_HR",2],boot_ci["Adjusted_HR",2],boot_ci["Landmark48_HR",2],
        NA,NA,boot_ci["AIPW_RD",2],boot_ci["AIPW_RR",2],
        NA,NA,boot_ci["RMST28_diff",2],
        NA,NA,boot_ci["PRBC48_RD",2],boot_ci["PRBC48_RR",2],
        NA,NA,boot_ci["AKI7_RD",2],boot_ci["AKI7_RR",2],boot_ci["NutrNone_HR",2]),
  model_wald_lcl=c(ccw_p[2],NA,NA,NA,iptw_p[2],adj_p[2],lmk_p[2],
        rep(NA,15),nn_p[3]),
  model_wald_ucl=c(ccw_p[3],NA,NA,NA,iptw_p[3],adj_p[3],lmk_p[3],
        rep(NA,15),nn_p[4]),
  ci_type=c("bootstrap B=200 + robust Wald shown","","","bootstrap B=200",
        "bootstrap B=200 + robust Wald","bootstrap B=200 + robust Wald","bootstrap B=200 + robust Wald",
        "","","bootstrap B=200","bootstrap B=200","","","bootstrap B=200",
        "","","bootstrap B=200","bootstrap B=200","","","bootstrap B=200","bootstrap B=200",
        "bootstrap B=200 + robust Wald"),
  n_complete=nrow(cc), n_landmark=nrow(lmk0), n_nutrnnone=as.integer(nn_p[1]))
write.csv(out,file.path(OUT,"pc_clinical_ci.csv"),row.names=FALSE)
cat("\n=== POINT + BOOTSTRAP CI ===\n");print(out,row.names=FALSE,digits=4)
cat("\nclones=",nrow(cl),"censored=",cl[,sum(censored)]," sw trim=",round(cwq,3),
    " n_landmark=",nrow(lmk0)," (",lmk0[,sum(A==1)],"/",lmk0[,sum(A==0)],")\n")
