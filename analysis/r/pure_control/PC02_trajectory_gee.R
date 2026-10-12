# =====================================================================
# PC02: Phosphate trajectory — crude summaries + IPCW-weighted GEE
# A_IV vs D_none (C1), B_PO vs D (C3); post-dose alignment for IV
# =====================================================================
suppressMessages({library(DBI); library(RPostgres); library(dplyr); library(tidyr)
  library(geepack); library(readr)}
)
set.seed(20260629L)
OUT <- "/data/小米云盘/文档/论文/血磷/TTE/v28.4/output/pure_control"

con <- dbConnect(Postgres(), host="localhost", port=5432, dbname="mimiciv",
                 user="postgres", password=Sys.getenv("MIMIC_PG_PASSWORD"))
win <- dbGetQuery(con, "SELECT * FROM public.pure_control_phos_windows")
co  <- dbGetQuery(con, "SELECT * FROM public.pure_control_cohort")
lng  <- dbGetQuery(con, "SELECT stay_id, hph, phos FROM public.pure_control_phos_long")
dbDisconnect(con)

co <- co %>% mutate(sexM=as.integer(gender=="M"),
  sepsis3=as.integer(sepsis3), vent_dim3=as.integer(vent_dim3),
  vasoactive_at_t0=as.integer(vasoactive_at_t0), refeeding_risk=as.integer(refeeding_risk))
DAG14 <- c("age","sexM","baseline_sofa","sepsis3","charlson","egfr","vent_dim3",
           "vasoactive_at_t0","admission_type","baseline_phosphate","lactate_t0",
           "refeeding_risk","corrected_calcium_mmol_min","nutrition_support_type")

wcent <- c(w6=6,w12=12,w24=24,w48=48,w72=72)
wlong <- win %>% mutate(time=wcent[win]) %>% select(stay_id, exposure_w48, time, actual_hph, phos)
base <- co %>% select(stay_id, baseline_phosphate, all_of(DAG14))
wlong <- wlong %>% left_join(base, by="stay_id")

# ---------- 1. CRUDE summary: 3 groups x 5 windows ----------
crude <- wlong %>% group_by(exposure_w48, time) %>%
  summarise(n=n(), mean=mean(phos), sd=sd(phos),
            median=median(phos), q1=quantile(phos,.25), q3=quantile(phos,.75),
            mean_delta=mean(phos-baseline_phosphate, na.rm=TRUE), .groups="drop")
gn <- co %>% filter(exposure_w48 %in% c("A_IV","B_PO_only","D_none")) %>%
  count(exposure_w48, name="group_n")
crude <- crude %>% left_join(gn, by="exposure_w48") %>%
  mutate(coverage=round(100*n/group_n,1))
write.csv(crude, file.path(OUT,"pc_crude_trajectory.csv"), row.names=FALSE)

# group contrasts (mean differences vs D_none)
ctr <- crude %>% select(exposure_w48,time,mean,mean_delta) %>%
  pivot_wider(names_from=exposure_w48, values_from=c(mean,mean_delta)) %>%
  mutate(IV_minus_none=mean_A_IV-mean_D_none,
         PO_minus_none=mean_B_PO_only-mean_D_none,
         dIV_vs_dNone=mean_delta_A_IV-mean_delta_D_none)
write.csv(ctr, file.path(OUT,"pc_crude_contrasts.csv"), row.names=FALSE)

# ---------- 2. IPCW measurement probability + IPTW x IPCW GEE ----------
run_gee <- function(dat, trt_label) {
  # dat: A/B vs D window-long with DAG14; complete covariates for weight models
  d <- dat %>% filter(complete.cases(across(all_of(DAG14))))
  # balanced panel grid: each stay x 5 windows; measured flag
  grid <- expand.grid(stay_id=unique(d$stay_id), time=c(6,12,24,48,72),
                      stringsAsFactors=FALSE)
  d0 <- co %>% filter(stay_id %in% d$stay_id) %>%
    select(stay_id, exposure_w48, baseline_phosphate, all_of(DAG14)) %>%
    filter(complete.cases(across(all_of(DAG14))))
  grid <- grid %>% left_join(d0, by="stay_id") %>%
    left_join(wlong %>% select(stay_id,time,phos), by=c("stay_id","time")) %>%
    mutate(measured=as.integer(!is.na(phos)),
           A=as.integer(exposure_w48==trt_label),
           timef=factor(time))
  # PS (one per stay) — reuse stabilized IPTW
  ps <- glm(as.formula(paste("A ~",paste(DAG14,collapse="+"))),
            data=grid %>% distinct(stay_id,.keep_all=TRUE), family=binomial)
  ed <- grid %>% distinct(stay_id,.keep_all=TRUE)
  ed$e <- predict(ps, newdata=ed, type="response")
  pA <- mean(ed$A); ed$iptw <- ifelse(ed$A==1,pA/ed$e,(1-pA)/(1-ed$e))
  # measurement model (pooled over windows)
  mm <- glm(measured ~ A*timef + age + baseline_sofa + sepsis3 + egfr +
              vent_dim3 + vasoactive_at_t0 + baseline_phosphate + charlson,
            data=grid, family=binomial)
  grid$pm <- predict(mm, type="response")
  # stabilized IPCW: numerator P(meas|A,time), denominator full model
  mm0 <- glm(measured ~ A*timef, data=grid, family=binomial)
  grid$pm0 <- predict(mm0, type="response")
  grid <- grid %>% left_join(ed %>% select(stay_id,iptw), by="stay_id") %>%
    mutate(ipcw=pm0/pm, w=iptw*ipcw)
  wq <- quantile(grid$w, c(.01,.99), na.rm=TRUE)
  grid$wt <- pmin(pmax(grid$w,wq[1]),wq[2])
  obs <- grid %>% filter(measured==1)
  fit <- geeglm(phos ~ factor(A)*timef + baseline_phosphate, id=stay_id,
                weights=wt, data=obs, corstr="exchangeable")
  # marginal predicted means via least-squares means style manual contrast
  le <- expand.grid(A=c(1,0), timef=factor(c(6,12,24,48,72)),
                    baseline_phosphate=mean(obs$baseline_phosphate))
  X <- model.matrix(~ factor(A)*timef + baseline_phosphate, data=le)
  eta <- as.numeric(X %*% coef(fit))
  se  <- sqrt(diag(X %*% vcov(fit) %*% t(X)))
  marg <- data.frame(comparison=ifelse(trt_label=="A_IV","C1_IV_vs_none","C3_PO_vs_none"),
                     time=as.numeric(as.character(le$timef)), A=le$A,
                     est=eta, se=se) %>%
    pivot_wider(names_from=A, values_from=c(est,se), names_sep="")
  wdiag <- grid %>% distinct(stay_id,time,.keep_all=TRUE) %>%
    summarise(meas_cov=mean(measured), ipcw_p1=quantile(ipcw,.01,na.rm=TRUE),
              ipcw_p99=quantile(ipcw,.99,na.rm=TRUE),
              totalw_p99=quantile(w,.99,na.rm=TRUE), pct_trunc=100*mean(w<wq[1]|w>wq[2],na.rm=TRUE))
  list(marg=marg, wdiag=wdiag)
}

g1 <- run_gee(wlong %>% filter(exposure_w48 %in% c("A_IV","D_none")), "A_IV")
g3 <- run_gee(wlong %>% filter(exposure_w48 %in% c("B_PO_only","D_none")), "B_PO_only")
gee <- bind_rows(mutate(g1$marg), mutate(g3$marg)) %>%
  mutate(weighted_diff=est1-est0, diff_se=sqrt(se1^2+se0^2),
         lcl=weighted_diff-1.96*diff_se, ucl=weighted_diff+1.96*diff_se)
write.csv(gee, file.path(OUT,"pc_gee_contrasts.csv"), row.names=FALSE)
wdiag <- bind_rows(mutate(g1$wdiag,comparison="C1"),
                   mutate(g3$wdiag,comparison="C3"))
write.csv(wdiag, file.path(OUT,"pc_ipcw_diagnostics.csv"), row.names=FALSE)

# ---------- 3. post-dose alignment (IV only) ----------
iv <- lng %>% inner_join(co %>% filter(exposure_w48=="A_IV") %>%
      select(stay_id, t0, iv_first_w48, baseline_phosphate), by="stay_id") %>%
  mutate(iv_h=as.numeric(difftime(iv_first_w48, t0, units="hours")),
         h_post=hph-iv_h) %>% filter(h_post>=0, h_post<=48)
pd <- iv %>% mutate(pwin=cut(h_post, c(0,4,12,24,48),
                             labels=c("0-4h","4-12h","12-24h","24-48h"), include.lowest=TRUE)) %>%
  group_by(pwin) %>% summarise(n=n(), n_stays=n_distinct(stay_id),
    mean=mean(phos), mean_delta=mean(phos-baseline_phosphate), .groups="drop")
write.csv(pd, file.path(OUT,"pc_postdose_iv.csv"), row.names=FALSE)

# ---------- 4. 1h plot data ----------
plot1h <- lng %>% inner_join(co %>% select(stay_id,exposure_w48), by="stay_id") %>%
  filter(exposure_w48 %in% c("A_IV","B_PO_only","D_none")) %>%
  mutate(h=round(hph)) %>% filter(h>=0,h<=72) %>%
  group_by(exposure_w48,h) %>% summarise(n=n(),mean=mean(phos),.groups="drop") %>%
  filter(n>=20)
write.csv(plot1h, file.path(OUT,"pc_plot_1h.csv"), row.names=FALSE)

cat("=== CRUDE CONTRASTS (IV/PO minus none) ===\n"); print(as.data.frame(ctr))
cat("\n=== WEIGHTED GEE CONTRASTS (W48 strategy, absolute level) ===\n"); print(as.data.frame(gee %>% select(comparison,time,weighted_diff,lcl,ucl)))
cat("\n=== IPCW DIAGNOSTICS ===\n"); print(as.data.frame(wdiag))
cat("\n=== POST-DOSE IV ===\n"); print(as.data.frame(pd))

# =====================================================================
# 5. DELTA-from-baseline GEE (W48)  +  EARLY (W12) cohort GEE
# =====================================================================
run_gee2 <- function(d_class, trt_label, outcol, comp_name) {
  # rebuild window-long with class assignment passed explicitly
  wl <- win %>% mutate(time=wcent[win]) %>%
    select(stay_id,time,actual_hph,phos) %>%
    inner_join(co %>% select(stay_id,cls=all_of(d_class),baseline_phosphate,all_of(DAG14)),
               by="stay_id") %>% filter(cls %in% c(trt_label,"D_none"))
  d <- wl %>% filter(complete.cases(across(all_of(DAG14))))
  d0 <- co %>% select(stay_id,cls=all_of(d_class),baseline_phosphate,all_of(DAG14)) %>%
    filter(stay_id %in% d$stay_id, complete.cases(across(all_of(DAG14))))
  grid <- expand.grid(stay_id=unique(d0$stay_id), time=c(6,12,24,48,72), stringsAsFactors=FALSE) %>%
    left_join(d0, by="stay_id") %>%
    left_join(wl %>% select(stay_id,time,phos), by=c("stay_id","time")) %>%
    mutate(measured=as.integer(!is.na(phos)), A=as.integer(cls==trt_label),
           timef=factor(time), delta=phos-baseline_phosphate)
  ed0 <- grid %>% distinct(stay_id,.keep_all=TRUE)
  ps <- glm(as.formula(paste("A ~",paste(DAG14,collapse="+"))), data=ed0, family=binomial)
  ed0$e <- predict(ps,type="response"); pA<-mean(ed0$A)
  ed0$iptw <- ifelse(ed0$A==1,pA/ed0$e,(1-pA)/(1-ed0$e))
  mm  <- glm(measured~A*timef+age+baseline_sofa+sepsis3+egfr+vent_dim3+vasoactive_at_t0+baseline_phosphate+charlson, data=grid, family=binomial)
  mm0 <- glm(measured~A*timef, data=grid, family=binomial)
  grid$pm<-predict(mm,type="response"); grid$pm0<-predict(mm0,type="response")
  grid <- grid %>% left_join(ed0%>%select(stay_id,iptw),by="stay_id") %>%
    mutate(w=iptw*pm0/pm)
  wq<-quantile(grid$w,c(.01,.99),na.rm=TRUE); grid$wt<-pmin(pmax(grid$w,wq[1]),wq[2])
  obs <- grid %>% filter(measured==1)
  yv <- obs[[outcol]]
  fit <- geeglm(yv ~ factor(A)*timef + baseline_phosphate, id=stay_id, weights=wt,
                data=obs, corstr="exchangeable")
  le <- expand.grid(A=c(1,0),timef=factor(c(6,12,24,48,72)),baseline_phosphate=mean(obs$baseline_phosphate))
  X<-model.matrix(~factor(A)*timef+baseline_phosphate,data=le); eta<-as.numeric(X%*%coef(fit)); se<-sqrt(diag(X%*%vcov(fit)%*%t(X)))
  data.frame(comparison=comp_name,time=as.numeric(as.character(le$timef)),A=le$A,est=eta,se=se) %>%
    pivot_wider(names_from=A,values_from=c(est,se),names_sep="") %>%
    mutate(diff=est1-est0, diff_se=sqrt(se1^2+se0^2),
           lcl=diff-1.96*diff_se,ucl=diff+1.96*diff_se) %>%
    select(comparison,time,diff,lcl,ucl)
}

gee_delta_w48 <- run_gee2("exposure_w48","A_IV","delta","C1_W48_delta")
gee_early_abs <- run_gee2("exposure_w12","A_IV","phos","C1_W12_absolute")
gee_early_delta <- run_gee2("exposure_w12","A_IV","delta","C1_W12_delta")
gee_extra <- bind_rows(gee_delta_w48,gee_early_abs,gee_early_delta)
write.csv(gee_extra, file.path(OUT,"pc_gee_extra.csv"), row.names=FALSE)
cat("\n=== GEE EXTRA: W48 delta / W12 early-treatment absolute + delta ===\n")
print(as.data.frame(gee_extra))
