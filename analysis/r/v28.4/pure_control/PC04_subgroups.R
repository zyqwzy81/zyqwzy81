# =====================================================================
# PC04: Six-layer subgroup analysis (A_IV vs D_none)
# Early pharmacodynamic effect: W12 class, w12 window delta
# Sustained correction:       W48 class, w48 window delta
# Overall stabilized IPTW; strata n<30/arm suppressed
# =====================================================================
suppressMessages({library(DBI);library(RPostgres);library(dplyr);library(tidyr);library(readr)})
set.seed(20260629L)
OUT <- "/data/小米云盘/文档/论文/血磷/TTE/v28.4/output/pure_control"
con <- dbConnect(Postgres(),host="localhost",port=5432,dbname="mimiciv",user="postgres",password=Sys.getenv("MIMIC_PG_PASSWORD"))
co <- dbGetQuery(con,"SELECT * FROM public.pure_control_cohort")
win<- dbGetQuery(con,"SELECT stay_id,win,phos FROM public.pure_control_phos_windows")
dbDisconnect(con)
co <- co %>% mutate(sexM=as.integer(gender=="M"),sepsis3=as.integer(sepsis3),
  vent_dim3=as.integer(vent_dim3),vasoactive_at_t0=as.integer(vasoactive_at_t0),
  refeeding_risk=as.integer(refeeding_risk),rrt_at_t0=as.integer(rrt_at_t0))
DAG14 <- c("age","sexM","baseline_sofa","sepsis3","charlson","egfr","vent_dim3",
           "vasoactive_at_t0","admission_type","baseline_phosphate","lactate_t0",
           "refeeding_risk","corrected_calcium_mmol_min","nutrition_support_type")

psw <- function(d){
  d <- d %>% filter(complete.cases(across(all_of(DAG14)))) %>% mutate(A=as.integer(A))
  f<-glm(as.formula(paste("A~",paste(DAG14,collapse="+"))),data=d,family=binomial)
  e<-predict(f,type="response");p<-mean(d$A)
  d %>% mutate(w=ifelse(A==1,p/e,(1-p)/(1-e))) %>% select(stay_id,w)
}

# assemble analysis frame for a given class column + window
frame <- function(classcol, winname){
  w <- win %>% filter(win==winname) %>% select(stay_id,phos)
  co %>% mutate(cls=.data[[classcol]]) %>%
    filter(cls %in% c("A_IV","D_none")) %>%
    inner_join(w,by="stay_id") %>% mutate(delta=phos-baseline_phosphate,
           A=as.integer(cls=="A_IV"))
}
f12 <- frame("exposure_w12","w12") %>% inner_join(psw(.),by="stay_id")
f48 <- frame("exposure_w48","w48") %>% inner_join(psw(.),by="stay_id")

# stratum labels
add_strata <- function(d){
  d %>% mutate(
    s_baseline_p = case_when(baseline_phosphate<0.50~"<0.50",
      baseline_phosphate>=0.50&baseline_phosphate<0.65~"0.50-0.65",
      baseline_phosphate>=0.65&baseline_phosphate<0.80~"0.65-0.80", TRUE~NA_character_),
    s_renal = case_when(rrt_at_t0==1~"RRT", is.na(egfr)~NA_character_,
      egfr<30~"eGFR<30", egfr<60~"eGFR30-60", TRUE~"eGFR>=60"),
    s_sofa = case_when(baseline_sofa<6~"SOFA<6", baseline_sofa<10~"SOFA6-9", TRUE~"SOFA>=10"),
    s_sepsis = ifelse(sepsis3==1,"sepsis","no_sepsis"),
    s_nutrition = ifelse(nutrition_support_type=="none","none","any_nutrition"),
    s_vent = ifelse(vent_dim3==1,"invasive_vent","no_vent"))
}
f12 <- add_strata(f12); f48 <- add_strata(f48)

wmean_diff <- function(d, svar){
  do.call(bind_rows,lapply(na.omit(unique(d[[svar]])),function(lev){
    ss <- d %>% filter(.data[[svar]]==lev)
    nt <- sum(ss$A==1); nc <- sum(ss$A==0)
    if(nt<30||nc<30)
      return(data.frame(stratum_var=svar,level=lev,n_trt=nt,n_ctrl=nc,
        delta_trt=NA,delta_ctrl=NA,diff=NA,lcl=NA,ucl=NA,note="NA-insufficient sample"))
    g <- ss %>% group_by(A) %>% reframe(wm=weighted.mean(delta,w))
    # weighted variance (inverse-probability, naive independent-group SE)
    se_one <- function(dd){s<-dd$w; sum(s^2*(dd$delta-weighted.mean(dd$delta,s))^2)/(sum(s)^2)}
    vt<-se_one(ss%>%filter(A==1));vc<-se_one(ss%>%filter(A==0));se<-sqrt(vt+vc)
    mt<-g$wm[g$A==1];mc<-g$wm[g$A==0];df<-mt-mc
    data.frame(stratum_var=svar,level=lev,n_trt=nt,n_ctrl=nc,
      delta_trt=round(mt,3),delta_ctrl=round(mc,3),diff=round(df,3),
      lcl=round(df-1.96*se,3),ucl=round(df+1.96*se,3),note="")
  }))
}
strata <- c("s_baseline_p","s_renal","s_sofa","s_sepsis","s_nutrition","s_vent")
r12 <- bind_rows(lapply(strata,function(s) wmean_diff(f12,s))) %>% mutate(analysis="W12_early_12h_delta")
r48 <- bind_rows(lapply(strata,function(s) wmean_diff(f48,s))) %>% mutate(analysis="W48_48h_delta")
res <- bind_rows(r12,r48) %>% select(analysis,stratum_var,level,n_trt,n_ctrl,delta_trt,delta_ctrl,diff,lcl,ucl,note)
write.csv(res,file.path(OUT,"pc_subgroups.csv"),row.names=FALSE)

# IV timing (within A_IV, W48 class) — pharmacodynamic descriptive
tim <- f48 %>% filter(A==1) %>%
  mutate(timing=case_when(is.na(iv_first_w48)~NA_character_,
    difftime(iv_first_w48,t0,units="hours")<=10~"early<=10h",TRUE~"late10-48h")) %>%
  filter(!is.na(timing)) %>% group_by(timing) %>%
  summarise(n=n(),m48=weighted.mean(phos,w),d48=weighted.mean(delta,w),.groups="drop")
write.csv(tim,file.path(OUT,"pc_iv_timing_48h.csv"),row.names=FALSE)
cat("=== SUBGROUP weighted delta difference (IV - none) ===\n")
write.table(as.data.frame(res),quote=FALSE,row.names=FALSE,sep="\t")
cat("\n=== IV TIMING 48h (within IV) ===\n")
write.table(as.data.frame(tim),quote=FALSE,row.names=FALSE,sep="\t")
