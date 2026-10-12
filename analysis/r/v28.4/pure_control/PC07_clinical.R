# =====================================================================
# PC07: Post-hoc target-trial clinical outcomes, C1 = A_IV vs D_none
# POST-HOC / hypothesis-generating; does NOT replace main analysis.
# =====================================================================
suppressMessages({library(DBI);library(RPostgres);library(dplyr);library(tidyr)
  library(survival);library(arrow);library(readr);library(data.table)})
set.seed(20260629L)
B <- 200L
OUT <- "/data/小米云盘/文档/论文/血磷/TTE/v28.4/output/pure_control"
BASE <- "/data/小米云盘/文档/论文/血磷"

con <- dbConnect(Postgres(),host="localhost",port=5432,dbname="mimiciv",user="postgres",password=Sys.getenv("MIMIC_PG_PASSWORD"))
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
aki_ev <- kd %>% filter(!is.na(aki_stage),aki_stage>=2) %>%
  inner_join(t0v,by="stay_id") %>%
  mutate(dd=as.numeric(difftime(charttime,t0,units="days"))) %>%
  filter(dd>0,dd<=7) %>% distinct(stay_id) %>% mutate(aki7=1L)

an <- co %>% filter(exposure_w48 %in% c("A_IV","D_none")) %>%
  inner_join(mor,by="stay_id") %>%
  left_join(prbc_ev,by="stay_id") %>% left_join(aki_ev,by="stay_id") %>%
  mutate(prbc48=tidyr::replace_na(prbc48,0L),aki7=tidyr::replace_na(aki7,0L),
         A=as.integer(exposure_w48=="A_IV"),
         iv_day=as.numeric(difftime(iv_first_w48,t0,units="days")))
# recode factors to binary categories populated in BOTH arms (avoids unseen
# levels / separation in per-arm outcome models and censor models)
an$admission_type <- ifelse(grepl("ELECTIVE|SURGICAL|AMBULATORY",an$admission_type),
                            "elective","emergency")
an$nutrition_support_type <- ifelse(an$nutrition_support_type=="none","none","any")
cc <- as.data.table(an %>% filter(complete.cases(across(all_of(DAG14)))))
cat("C1 n=",nrow(an)," complete-case=",nrow(cc),
    " deaths A/D=",cc[A==1,sum(death28)],"/",cc[A==0,sum(death28)],"\n")
cat("prbc48 A/D=",cc[A==1,sum(prbc48)],"/",cc[A==0,sum(prbc48)],
    " aki7 A/D=",cc[A==1,sum(aki7)],"/",cc[A==0,sum(aki7)],"\n")

psf <- glm(as.formula(paste("A~",paste(DAG14,collapse="+"))),data=cc,family=binomial)
cc$e <- predict(psf,type="response"); pA <- mean(cc$A)
cc$w  <- ifelse(cc$A==1,pA/cc$e,(1-pA)/(1-cc$e))
wq <- quantile(cc$w,c(.01,.99)); cc$wt <- pmin(pmax(cc$w,wq[1]),wq[2])

# ---- CCW: clone at T0; 48h(=2d) grace; clone-level cumulative IPCW ----
build_ccw <- function(d){
  grid <- d[, .(day=seq(0,floor(pmin(surv_days[1],28)))), by=stay_id]
  pd <- grid[, .(z=c(1L,0L)), by=.(stay_id,day)]
  setkey(pd,stay_id)
  meta <- d[, .(stay_id,A,death28,surv_days,iv_day,
                age,sexM,baseline_sofa,sepsis3,charlson,egfr,vent_dim3,
                vasoactive_at_t0,admission_type,baseline_phosphate,lactate_t0,
                refeeding_risk,corrected_calcium_mmol_min,nutrition_support_type)]
  pd <- merge(pd,meta,by="stay_id",allow.cartesian=TRUE)
  pd[, `:=`(
    iv_day=ifelse(is.na(iv_day),99,iv_day),
    cens=fifelse(z==1 & A==0 & day>=2, 1L,
          fifelse(z==0 & A==1 & iv_day<=day, 1L, 0L)))]
  pd[, firstc := min(c(day[cens==1],29)), by=.(stay_id,z)]
  # uncensored-death event indicator only counts before censor
  pd[, evday := as.integer(death28==1 & day==floor(surv_days) & day<=firstc)]
  fit0 <- glm(cens~z+poly(day,2),data=pd,family=binomial)
  fit1 <- glm(as.formula(paste("cens~z*poly(day,2)+",paste(DAG14,collapse="+"))),data=pd,family=binomial)
  pd[, `:=`(p0=pmin(pmax(predict(fit0,pd,type="response"),1e-4),1-1e-4),
            p1=pmin(pmax(predict(fit1,pd,type="response"),1e-4),0.99))]
  pd <- pd[day<firstc | (day==firstc & cens==0)]   # drop the censor day itself
  pd[, sw := cumprod((1-p0)/(1-p1)), by=.(stay_id,z)]
  pd[, sw:=sw[.N], by=.(stay_id,z)]                # clone-level cumulative weight
  cl <- pd[, .(sw=sw[1],
               time=fifelse(firstc[1]<29, firstc[1], pmin(surv_days[1],28)),
               event=fifelse(firstc[1]<29,0L, as.integer(death28[1]==1 & surv_days[1]<=28)),
               z=z[1], censored=as.integer(firstc[1]<29)), by=.(stay_id,z)]
  cl
}
cl <- build_ccw(cc)
cl[, sw:=pmin(pmax(sw,quantile(sw,.01)),quantile(sw,.99))]
ccw <- coxph(Surv(time,event)~z,data=cl,weights=sw,robust=TRUE)
hr_ccw <- c(HR=unname(exp(coef(ccw))),exp(confint(ccw)))
cat("CCW clones=",nrow(cl),"censored clones=",cl[,sum(censored)],
    " HR=",round(hr_ccw[1],3),"\n")

# ---- IPTW Cox (bootstrap CI) ----
iptw <- coxph(Surv(surv_days,death28)~A,data=cc,weights=wt,robust=TRUE)
hr_iptw <- c(HR=unname(exp(coef(iptw))),as.numeric(exp(confint(iptw))))
adj <- coxph(as.formula(paste("Surv(surv_days,death28)~A+",paste(DAG14,collapse="+"))),
             data=cc,weights=wt,robust=TRUE)
hr_adj <- c(HR=unname(exp(coef(adj)["A"])),as.numeric(exp(confint(adj)["A",])))
# PH diagnostic
zph <- cox.zph(coxph(Surv(surv_days,death28)~A,data=cc,weights=wt))
ph_p <- zph$table["A","p"]
# early (within-48h) death: treatment-withdrawal selection signal
early_death <- cc[,.(rate=mean(death28==1&surv_days<=2)),by=A]
# LANDMARK 48h: remove first-48h deaths, shift origin -> removes early
# treatment-withdrawal selection; the bias-robust sensitivity estimate
lmk <- cc[surv_days>2]; lmk[,`:=`(t2=surv_days-2,y=death28)]
cml <- coxph(Surv(t2,y)~A,data=lmk,weights=wt,robust=TRUE)
hr_lmk <- c(HR=unname(exp(coef(cml))),as.numeric(exp(confint(cml))))

# ---- AIPW ----
m1 <- glm(as.formula(paste("death28~",paste(DAG14,collapse="+"))),data=cc[A==1],family=binomial)
m0 <- glm(as.formula(paste("death28~",paste(DAG14,collapse="+"))),data=cc[A==0],family=binomial)
cc$q1<-predict(m1,cc,type="response");cc$q0<-predict(m0,cc,type="response")
cc$ip1<-cc$q1+cc$A/cc$e*(cc$death28-cc$q1)
cc$ip0<-cc$q0+(1-cc$A)/(1-cc$e)*(cc$death28-cc$q0)
r1<-weighted.mean(cc$ip1,cc$wt);r0<-weighted.mean(cc$ip0,cc$wt)
aip<-c(risk1=r1,risk0=r0,RD=r1-r0,RR=r1/r0)

# ---- RMST28 IPTW ----
rmst<-function(d,g){dd<-d[A==g];sf<-survfit(Surv(surv_days,death28)~1,data=dd,weights=wt)
  t<-c(0,sf$time);S<-c(1,sf$surv);k<-t<=28;sum(diff(c(t[k],28))*S[k])}
rm1<-rmst(cc,1);rm0<-rmst(cc,0)

ev<-function(x){x<-ifelse(x<1,1/x,x);1.41*x+sqrt(1.41^2*x*(x-1))-0.41}
evalue<-c(point=ev(hr_iptw[1]),bound=ev(ifelse(hr_iptw[1]<1,hr_iptw[3],hr_iptw[2])))

# ---- bootstrap HR + RD ----
n<-nrow(cc);R<-matrix(NA_real_,B,2)
for(b in seq_len(B)){
  db<-cc[sample(seq_len(n),n,replace=TRUE)]
  fb<-try(glm(as.formula(paste("A~",paste(DAG14,collapse="+"))),data=db,family=binomial),silent=TRUE)
  if(inherits(fb,"try-error"))next
  eb<-predict(fb,type="response");pb<-mean(db$A)
  db$wb<-pmin(pmax(ifelse(db$A==1,pb/eb,(1-pb)/(1-eb)),wq[1]),wq[2])
  cb<-try(coxph(Surv(surv_days,death28)~A,data=db,weights=wb,robust=TRUE),silent=TRUE)
  g1<-try(glm(as.formula(paste("death28~",paste(DAG14,collapse="+"))),data=db[A==1],family=binomial),silent=TRUE)
  g0<-try(glm(as.formula(paste("death28~",paste(DAG14,collapse="+"))),data=db[A==0],family=binomial),silent=TRUE)
  if(!inherits(cb,"try-error"))R[b,1]<-exp(coef(cb))
  if(!inherits(g1,"try-error")&&!inherits(g0,"try-error")){
    q1<-predict(g1,db,type="response");q0<-predict(g0,db,type="response")
    ip1<-q1+db$A/eb*(db$death28-q1);ip0<-q0+(1-db$A)/(1-eb)*(db$death28-q0)
    R[b,2]<-weighted.mean(ip1,db$wb)-weighted.mean(ip0,db$wb)}
}
hrb<-quantile(R[,1],c(.025,.975),na.rm=TRUE);rdb<-quantile(R[,2],c(.025,.975),na.rm=TRUE)

# ---- secondary ----
bin<-function(col){y<-cc[[col]];x1<-weighted.mean(y[cc$A==1],cc$wt[cc$A==1]);x0<-weighted.mean(y[cc$A==0],cc$wt[cc$A==0])
  c(trt_pct=100*x1,ctrl_pct=100*x0,RD_pct=100*(x1-x0),RR=x1/x0)}
sec<-bind_rows(as.data.frame(t(bin("prbc48")))%>%mutate(outcome="prbc48"),
               as.data.frame(t(bin("aki7")))%>%mutate(outcome="aki7"))

# ---- nutrition=none ----
nn<-cc[nutrition_support_type=="none"]
DAGnn <- setdiff(DAG14,"nutrition_support_type")  # constant within subset
fn<-glm(as.formula(paste("A~",paste(DAGnn,collapse="+"))),data=nn,family=binomial)
nn$en<-predict(fn,type="response");pn<-mean(nn$A)
nn$wn<-pmin(pmax(ifelse(nn$A==1,pn/nn$en,(1-pn)/(1-nn$en)),wq[1]),wq[2])
cn<-coxph(Surv(surv_days,death28)~A,data=nn,weights=wn,robust=TRUE)
hrn<-c(n=nrow(nn),HR=unname(exp(coef(cn))),as.numeric(exp(confint(cn))))

# E-value for landmark (bias-robust) HR too
ev_lmk <- c(point=ev(hr_lmk[1]),bound=ev(ifelse(hr_lmk[1]<1,hr_lmk[3],hr_lmk[2])))
main<-data.frame(method=c("CCW_Cox","IPTW_Cox_full0_28","Adjusted_Cox","IPTW_Cox_landmark48h",
  "PH_test_p","early48h_death_IV","early48h_death_none",
  "AIPW_risk1","AIPW_risk0","AIPW_RD","AIPW_RR","RMST28_IV","RMST28_none","RMST_diff",
  "Evalue_full_point","Evalue_full_bound","Evalue_lmk_point","Evalue_lmk_bound"),
  est=c(hr_ccw[1],hr_iptw[1],hr_adj[1],hr_lmk[1],ph_p,
        early_death[A==1]$rate,early_death[A==0]$rate,
        r1,r0,aip["RD"],aip["RR"],rm1,rm0,rm1-rm0,
        evalue[1],evalue[2],ev_lmk[1],ev_lmk[2]),
  lcl=c(hr_ccw[2],hrb[1],hr_adj[2],hr_lmk[2],NA,NA,NA,NA,NA,rdb[1],NA,NA,NA,NA,NA,NA,NA,NA),
  ucl=c(hr_ccw[3],hrb[2],hr_adj[3],hr_lmk[3],NA,NA,NA,NA,NA,rdb[2],NA,NA,NA,NA,NA,NA,NA,NA))
write.csv(main,file.path(OUT,"pc_clinical_main.csv"),row.names=FALSE)
write.csv(sec,file.path(OUT,"pc_clinical_secondary.csv"),row.names=FALSE)
write.csv(data.frame(subset="nutrition_none",n=hrn["n"],HR=hrn["HR"],lcl=unname(hrn[3]),ucl=unname(hrn[4])),
          file.path(OUT,"pc_clinical_nutrition_none.csv"),row.names=FALSE)
cat("\n=== MAIN ===\n");write.table(as.data.frame(main),quote=FALSE,row.names=FALSE,sep="\t")
cat("\n=== SECONDARY ===\n");write.table(as.data.frame(sec),quote=FALSE,row.names=FALSE,sep="\t")
cat("\n=== NUTRITION NONE ===\n");print(hrn)
