# PC07b: diagnose clinical-outcome discrepancy
suppressMessages({library(DBI);library(RPostgres);library(dplyr);library(survival);
  library(arrow);library(data.table)})
set.seed(20260629L)
con<-dbConnect(Postgres(),host="localhost",port=5432,dbname="mimiciv",user="postgres",password=Sys.getenv("MIMIC_PG_PASSWORD"))
co<-dbGetQuery(con,"SELECT * FROM public.pure_control_cohort")
mor<-dbGetQuery(con,"SELECT stay_id,death28,surv_days FROM public.pure_control_mortality");dbDisconnect(con)
co<-co%>%mutate(sexM=as.integer(gender=="M"),sepsis3=as.integer(sepsis3),
 vent_dim3=as.integer(vent_dim3),vasoactive_at_t0=as.integer(vasoactive_at_t0),
 refeeding_risk=as.integer(refeeding_risk),
 admission_type=ifelse(grepl("ELECTIVE|SURGICAL|AMBULATORY",admission_type),"elective","emergency"),
 nutrition_support_type=ifelse(nutrition_support_type=="none","none","any"))
DAG14<-c("age","sexM","baseline_sofa","sepsis3","charlson","egfr","vent_dim3","vasoactive_at_t0",
 "admission_type","baseline_phosphate","lactate_t0","refeeding_risk","corrected_calcium_mmol_min","nutrition_support_type")
an<-co%>%filter(exposure_w48%in%c("A_IV","D_none"))%>%inner_join(mor,by="stay_id")%>%
  mutate(A=as.integer(exposure_w48=="A_IV"))
cc<-as.data.table(an%>%filter(complete.cases(across(all_of(DAG14)))))
ps<-glm(as.formula(paste("A~",paste(DAG14,collapse="+"))),data=cc,family=binomial)
cc$e<-predict(ps,type="response");p<-mean(cc$A)
cc$w<-ifelse(cc$A==1,p/cc$e,(1-p)/(1-cc$e))
wq<-quantile(cc$w,c(.01,.99));cc$wt<-pmin(pmax(cc$w,wq[1]),wq[2])

cat("=== crude / weighted 28d risk ===\n")
cc[,.(n=.N,deaths=sum(death28),risk=mean(death28)),by=A]
cc[,.(wrisk=weighted.mean(death28,wt)),by=A]

cat("\n=== weighted KM at 28d + RMST (survfit) ===\n")
sf<-survfit(Surv(surv_days,death28)~A,data=cc,weights=wt)
print(summary(sf,times=28)$table)
print(setNames(summary(sf,times=28)$surv,c("A0","A1")))
for(g in 1:0){dd<-cc[A==g];s<-survfit(Surv(surv_days,death28)~1,data=dd,weights=wt)
 t<-c(0,s$time);S<-c(1,s$surv);k<-t<=28;cat("A=",g,"RMST28=",sum(diff(c(t[k],28))*S[k]),"\n")}

cat("\n=== PH assumption ===\n")
cx<-coxph(Surv(surv_days,death28)~A,data=cc,weights=wt,robust=TRUE)
print(summary(cx)$conf.int[c("A"),])
print(cox.zph(coxph(Surv(surv_days,death28)~A,data=cc,weights=wt)))

cat("\n=== weighted risk by time band (hazard timing) ===\n")
for(bd in list(c(0,3),c(3,7),c(7,28))){
  lo<-bd[1];hi<-bd[2]
  cc[,ev:=as.integer(death28==1&surv_days>lo&surv_days<=hi)]
  cc[,exposed:=as.integer(surv_days>lo)]
  cat(lo,"-",hi,"d: A rate%=",100*cc[A==1&exposed==1,weighted.mean(ev,wt)],
      " D rate%=",100*cc[A==0&exposed==1,weighted.mean(ev,wt)],"\n")
}

cat("\n=== weight by A ===\n");cc[,.(wmean=mean(wt),wmax=max(wt),p99=quantile(w,.99)),by=A]
cat("\n=== e distribution tails by death ===\n")
cc[,.(emean=mean(e),elo=quantile(e,.05),ehi=quantile(e,.95)),by=.(A,death28)]

cat("\n=== D_none: signals of withdrawal / very short stay ===\n")
cc[,.(n=.N,risk=mean(death28),med_los=median(surv_days),
      pct_die_d1=mean(death28==1&surv_days<=1),
      age=mean(age),sofa=mean(baseline_sofa),vent=mean(vent_dim3),
      pall=mean(baseline_sofa>=15)),by=A]
