# =====================================================================
# PC03: AUC0-48, 48h normalization, safety events (IPTW-weighted)
# =====================================================================
suppressMessages({library(DBI); library(RPostgres); library(dplyr); library(tidyr); library(readr)})
set.seed(20260629L)
OUT <- "/data/小米云盘/文档/论文/血磷/TTE/v28.4/output/pure_control"
con <- dbConnect(Postgres(),host="localhost",port=5432,dbname="mimiciv",user="postgres",password=Sys.getenv("MIMIC_PG_PASSWORD"))
co <- dbGetQuery(con,"SELECT * FROM public.pure_control_cohort")
lng <- dbGetQuery(con,"SELECT stay_id,hph,phos FROM public.pure_control_phos_long WHERE hph<=48")
saf <- dbGetQuery(con,"SELECT * FROM public.pure_control_safety")
w48 <- dbGetQuery(con,"SELECT stay_id,phos FROM public.pure_control_phos_windows WHERE win='w48'")
dbDisconnect(con)

# reuse stabilized IPTW (re-fit quickly, identical DAG14 specification)
co <- co %>% mutate(sexM=as.integer(gender=="M"),sepsis3=as.integer(sepsis3),
  vent_dim3=as.integer(vent_dim3),vasoactive_at_t0=as.integer(vasoactive_at_t0),
  refeeding_risk=as.integer(refeeding_risk))
DAG14 <- c("age","sexM","baseline_sofa","sepsis3","charlson","egfr","vent_dim3",
           "vasoactive_at_t0","admission_type","baseline_phosphate","lactate_t0",
           "refeeding_risk","corrected_calcium_mmol_min","nutrition_support_type")

ps_weights <- function(dat, label){
  d <- dat %>% filter(complete.cases(across(all_of(DAG14)))) %>% mutate(A=as.integer(exposure_w48==label))
  fit <- glm(as.formula(paste("A~",paste(DAG14,collapse="+"))),data=d,family=binomial)
  e <- predict(fit,type="response"); pA<-mean(d$A)
  d %>% mutate(w=ifelse(A==1,pA/e,(1-pA)/(1-e))) %>% select(stay_id,A,w)
}
wC1 <- ps_weights(co %>% filter(exposure_w48 %in% c("A_IV","D_none")),"A_IV")
wC3 <- ps_weights(co %>% filter(exposure_w48 %in% c("B_PO_only","D_none")),"B_PO_only")

# ---------- AUC0-48 trapezoid: >=3 points, first<=12h, last>=30h ----------
auc <- lng %>% arrange(stay_id,hph) %>% group_by(stay_id) %>%
  summarise(npts=n(), t0h=min(hph), tnh=max(hph),
            auc=sum(diff(hph)* (head(phos,-1)+tail(phos,-1))/2)/max(hph),
            .groups="drop") %>%
  filter(npts>=3, t0h<=12, tnh>=30)
auc <- auc %>% inner_join(co %>% select(stay_id,exposure_w48),by="stay_id") %>%
  filter(exposure_w48 %in% c("A_IV","B_PO_only","D_none"))
auc_tab <- auc %>% group_by(exposure_w48) %>%
  summarise(n=n(), mean_auc=mean(auc), sd=sd(auc), median=median(auc),.groups="drop")
write.csv(auc_tab,file.path(OUT,"pc_auc_summary.csv"),row.names=FALSE)
auc_w <- auc %>% inner_join(wC1 %>% select(stay_id,A,w),by="stay_id")
auc_c1 <- auc_w %>% group_by(A) %>% summarise(wmean=weighted.mean(auc,w),.groups="drop")

# ---------- weighted binary outcomes: safety events ----------
mk <- function(w,label){
  w %>% left_join(saf,by="stay_id") %>% mutate(A=ifelse(exposure_w48==label,1,0))
}
c1 <- wC1 %>% inner_join(co %>% select(stay_id,exposure_w48),by="stay_id") %>% mk("A_IV")
c3 <- wC3 %>% inner_join(co %>% select(stay_id,exposure_w48),by="stay_id") %>% mk("B_PO_only")

rate_contrast <- function(d,col,den_col=NULL){
  dd <- d %>% mutate(ev=.data[[col]])
  if(!is.null(den_col)) dd <- dd %>% filter(.data[[den_col]]>0)
  g <- dd %>% group_by(A) %>% summarise(r=weighted.mean(ev,w),.groups="drop")
  data.frame(outcome=col, trt_rate=g$r[g$A==1], ctrl_rate=g$r[g$A==0],
             RD=g$r[g$A==1]-g$r[g$A==0],
             n_trt=sum(dd$A==1),n_ctrl=sum(dd$A==0))
}
events <- c("hyperphos","hypok","hypomg","hypoca","hypoca_severe")
denmap <- c(hyperphos="phos_n48",hypok="k_n",hypomg="mg_n",hypoca="ca_n",hypoca_severe="ca_n")
saf_c1 <- bind_rows(lapply(events,function(v) rate_contrast(c1,v,denmap[[v]]))) %>% mutate(comparison="C1_IV_vs_none")
saf_c3 <- bind_rows(lapply(events,function(v) rate_contrast(c3,v,denmap[[v]]))) %>% mutate(comparison="C3_PO_vs_none")
# 48h normalization among measured
norm48_fn <- function(d){
  dd <- d %>% left_join(w48,by="stay_id") %>% filter(!is.na(phos)) %>%
    mutate(norm=as.integer(phos>=0.80))
  g<-dd%>%group_by(A)%>%summarise(r=weighted.mean(norm,w),.groups="drop")
  data.frame(outcome="norm48_among_measured",trt_rate=g$r[g$A==1],ctrl_rate=g$r[g$A==0],
             RD=g$r[g$A==1]-g$r[g$A==0],n_trt=sum(dd$A==1),n_ctrl=sum(dd$A==0))
}
allres <- bind_rows(saf_c1,saf_c3) %>%
  bind_rows(norm48_fn(c1) %>% mutate(comparison="C1_IV_vs_none")) %>%
  bind_rows(norm48_fn(c3) %>% mutate(comparison="C3_PO_vs_none")) %>%
  mutate(across(c(trt_rate,ctrl_rate,RD),~round(100*.x,1)))
write.csv(allres,file.path(OUT,"pc_safety_weighted.csv"),row.names=FALSE)

cat("=== AUC0-48 (mean normalized mmol/L over 48h) ===\n");print(as.data.frame(auc_tab))
cat("\nC1 weighted mean AUC IV:",round(filter(auc_c1,A==1)$wmean,3),
    " none:",round(filter(auc_c1,A==0)$wmean,3),"\n")
cat("\n=== WEIGHTED RATES (%) ===\n");print(as.data.frame(allres))
