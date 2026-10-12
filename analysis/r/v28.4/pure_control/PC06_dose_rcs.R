# =====================================================================
# PC06: Dose-response (first-24h IV mmol): quartiles + restricted cubic
# splines (3-4 knots) adjusted for DAG14, within A_IV; timing x dose
# =====================================================================
suppressMessages({library(DBI);library(RPostgres);library(dplyr);library(tidyr);library(rms);library(readr)})
set.seed(20260629L)
OUT <- "/data/小米云盘/文档/论文/血磷/TTE/v28.4/output/pure_control"
con <- dbConnect(Postgres(),host="localhost",port=5432,dbname="mimiciv",user="postgres",password=Sys.getenv("MIMIC_PG_PASSWORD"))
co <- dbGetQuery(con,"SELECT * FROM public.pure_control_cohort")
dose <- dbGetQuery(con,"SELECT stay_id,dose24,dose48,n_doses FROM public.pure_control_dose")
win <- dbGetQuery(con,"SELECT stay_id,win,phos FROM public.pure_control_phos_windows")
dbDisconnect(con)
co <- co %>% mutate(sexM=as.integer(gender=="M"),sepsis3=as.integer(sepsis3),
  vent_dim3=as.integer(vent_dim3),vasoactive_at_t0=as.integer(vasoactive_at_t0),
  refeeding_risk=as.integer(refeeding_risk))
DAG14 <- c("age","sexM","baseline_sofa","sepsis3","charlson","egfr","vent_dim3",
           "vasoactive_at_t0","admission_type","baseline_phosphate","lactate_t0",
           "refeeding_risk","corrected_calcium_mmol_min","nutrition_support_type")
w12 <- win%>%filter(win=="w12")%>%select(stay_id,p12=phos)
w48 <- win%>%filter(win=="w48")%>%select(stay_id,p48=phos)
d <- co%>%filter(exposure_w48=="A_IV")%>%inner_join(dose,by="stay_id")%>%
  left_join(w12,by="stay_id")%>%left_join(w48,by="stay_id")%>%
  mutate(d12=p12-baseline_phosphate,d48=p48-baseline_phosphate,
         iv_h=as.numeric(difftime(iv_first_w48,t0,units="hours")),
         early=as.integer(iv_h<=10))
cap <- quantile(d$dose24,0.99,na.rm=TRUE)          # trim extreme
d <- d %>% mutate(dose=pmin(dose24,cap))

# ---------- quartile table ----------
qt <- d %>% filter(!is.na(d48)) %>% mutate(q=ntile(dose24,4)) %>%
  group_by(q) %>% summarise(n=n(),dose_med=median(dose24),
    d12=mean(d12,na.rm=TRUE),d48=mean(d48,na.rm=TRUE),
    norm48=100*mean(p48>=0.80,na.rm=TRUE),.groups="drop")
write.csv(qt,file.path(OUT,"pc_dose_quartiles.csv"),row.names=FALSE)

# ---------- RCS dose-response (adjusted for DAG14) ----------
dd <<- datadist(d); options(datadist="dd")
fit_rcs <- function(out){
  da <- d %>% filter(!is.na(.data[[out]]),complete.cases(across(all_of(DAG14))))
  f <- ols(as.formula(paste(out,"~ rcs(dose, c(7.5,15,30,60)) +",paste(DAG14,collapse="+"))),data=da)
  g <- data.frame(dose=0:90)
  for(v in setdiff(DAG14,"baseline_phosphate")){
    g[[v]] <- if(is.numeric(da[[v]])) median(da[[v]],na.rm=TRUE) else names(which.max(table(da[[v]])))
  }
  g$baseline_phosphate <- median(da$baseline_phosphate)
  pred <- predict(f,newdata=g,se.fit=TRUE)
  g$fit <- as.numeric(pred$linear.predictors); g$se<-as.numeric(pred$se.fit)
  g %>% mutate(lcl=fit-1.96*se,ucl=fit+1.96*se,outcome=out)
}
rcs12 <- fit_rcs("d12"); rcs48 <- fit_rcs("d48")
rcs <- bind_rows(rcs12,rcs48)
write.csv(rcs,file.path(OUT,"pc_dose_rcs.csv"),row.names=FALSE)

# ---------- timing x dose 2x2 (48h delta) ----------
tg <- d %>% filter(!is.na(d48)) %>%
  mutate(timing=ifelse(early==1,"early<=10h","late10-48h"),
         dose_grp=ifelse(dose24>=30,"dose>=30","dose<30")) %>%
  group_by(timing,dose_grp) %>% summarise(n=n(),d48=mean(d48),p48=mean(p48),.groups="drop")
write.csv(tg,file.path(OUT,"pc_timing_dose.csv"),row.names=FALSE)

# early vs late vs none (none = D_none, crude delta)
none <- co%>%filter(exposure_w48=="D_none")%>%inner_join(w48,by="stay_id")%>%
  summarise(n=n(),d48=mean(p48-baseline_phosphate),p48=mean(p48))
eln <- bind_rows(
  d%>%filter(!is.na(d48),early==1)%>%summarise(n=n(),d48=mean(d48),p48=mean(p48))%>%mutate(group="early<=10h"),
  d%>%filter(!is.na(d48),early==0)%>%summarise(n=n(),d48=mean(d48),p48=mean(p48))%>%mutate(group="late10-48h"),
  none%>%mutate(group="none"))
write.csv(eln,file.path(OUT,"pc_early_late_none.csv"),row.names=FALSE)

cat("=== DOSE QUARTILES (first-24h mmol) ===\n");write.table(as.data.frame(qt),quote=FALSE,row.names=FALSE,sep="\t")
cat("\n=== RCS predicted delta at dose 0/15/30/60 ===\n")
write.table(as.data.frame(rcs%>%filter(dose%%15==0)%>%select(outcome,dose,fit,lcl,ucl)),quote=FALSE,row.names=FALSE,sep="\t")
cat("\n=== TIMING x DOSE ===\n");write.table(as.data.frame(tg),quote=FALSE,row.names=FALSE,sep="\t")
cat("\n=== EARLY/LATE/NONE ===\n");write.table(as.data.frame(eln),quote=FALSE,row.names=FALSE,sep="\t")
