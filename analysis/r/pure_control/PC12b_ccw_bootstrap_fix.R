# =====================================================================
# PC12b: corrected CCW-only subject-level bootstrap (review follow-up)
# Bug fixed vs PC12: resampling rows produced duplicated stay_id keys;
#   build_ccw() aggregates by stay_id and merges with allow.cartesian,
#   so duplicated keys created k^2 panel rows and cumulative products
#   over repeated blocks -> crushed weights -> downward-shifted CCW
#   bootstrap distribution (HR 0.51-0.58 vs point 0.654; RD CI wholly
#   below point). Fix: assign a fresh unique id to each resampled row.
# Full-data point estimates unchanged (unique ids in original data).
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
an <- co %>% filter(exposure_w48 %in% c("A_IV","D_none")) %>%
  inner_join(mor,by="stay_id") %>%
  mutate(A=as.integer(exposure_w48=="A_IV"),
         iv_day=as.numeric(difftime(iv_first_w48,t0,units="days")),
         admission_type=ifelse(grepl("ELECTIVE|SURGICAL|AMBULATORY",admission_type),
                               "elective","emergency"),
         nutrition_support_type=ifelse(nutrition_support_type=="none","none","any"))
cc <- as.data.table(an %>% filter(complete.cases(across(all_of(DAG14)))))

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
risk28 <- function(cl0,wcol="swb"){
  s1<-survfit(Surv(time,event)~1,data=cl0[z==1],weights=cl0[z==1][[wcol]])
  s0<-survfit(Surv(time,event)~1,data=cl0[z==0],weights=cl0[z==0][[wcol]])
  c(r1=1-summary(s1,times=28,extend=TRUE)$surv, r0=1-summary(s0,times=28,extend=TRUE)$surv)
}

cl <- build_ccw(cc)
cwq <- quantile(cl$sw,c(.01,.99)); cl[,swb:=pmin(pmax(sw,cwq[1]),cwq[2])]
ccw_hr_p <- unname(exp(coef(coxph(Surv(time,event)~z,data=cl,weights=swb,robust=TRUE)))["z"])
ccw_rd_p <- risk28(cl); rd_p <- ccw_rd_p["r1"]-ccw_rd_p["r0"]

n<-nrow(cc); R<-matrix(NA_real_,B,2); colnames(R)<-c("CCW_HR","CCW_RD28")
t0<-Sys.time()
for(b in seq_len(B)){
  db<-cc[sample(seq_len(n),n,replace=TRUE)]
  db[, stay_id := seq_len(.N)]          # fresh unique ids per resample
  clb<-try(build_ccw(db),silent=TRUE)
  if(inherits(clb,"try-error"))next
  clb[,swb:=pmin(pmax(sw,cwq[1]),cwq[2])]
  hb<-try(coxph(Surv(time,event)~z,data=clb,weights=swb,robust=TRUE),silent=TRUE)
  if(!inherits(hb,"try-error")) R[b,"CCW_HR"]<-unname(exp(coef(hb))["z"])
  rr<-try(risk28(clb),silent=TRUE)
  if(!inherits(rr,"try-error")) R[b,"CCW_RD28"]<-rr["r1"]-rr["r0"]
  if(b%%20==0) cat("b=",b," elapsed=",round(as.numeric(Sys.time()-t0,units="mins"),1),"min\n")
}

qfun <- function(v,nm){
  v<-v[!is.na(v)]; q<-quantile(v,c(.025,.975),na.rm=TRUE)
  data.frame(metric=nm,point=if(nm=="CCW_HR")ccw_hr_p else rd_p,
             boot_mean=mean(v),boot_lcl=as.numeric(q[1]),boot_ucl=as.numeric(q[2]),
             basic_lcl=NA_real_,basic_ucl=NA_real_,
             boot_na=sum(is.na(R[,nm])),row.names=NULL)
}
out <- rbind(qfun(R[,"CCW_HR"],"CCW_HR"),qfun(R[,"CCW_RD28"],"CCW_RD28"))
# basic (reverse-percentile) interval: 2*theta - q on the HR log scale, linear on RD
qh <- quantile(R[,"CCW_HR"],c(.025,.975),na.rm=TRUE)
qr <- quantile(R[,"CCW_RD28"],c(.025,.975),na.rm=TRUE)
out[out$metric=="CCW_HR",c("basic_lcl","basic_ucl")] <-
  exp(2*log(ccw_hr_p)-log(as.numeric(c(qh[2],qh[1]))))
out[out$metric=="CCW_RD28",c("basic_lcl","basic_ucl")] <- 2*rd_p-as.numeric(c(qr[2],qr[1]))
write_csv(out,file.path(OUT,"pc_ccw_bootstrap_fix.csv"))
print(out,digits=4)
cat("\nfull-data point: HR=",round(ccw_hr_p,4)," RD=",round(rd_p,4),
    " clones=",nrow(cl)," censored=",cl[,sum(censored)],"\n")
