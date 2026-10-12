# =====================================================================
# PC01: Cohort flowchart + DAG14 balance + stabilized IPTW weights
# Comparisons: C1 = A_IV vs D_none ; C3 = B_PO_only vs D_none
# =====================================================================
suppressMessages({library(DBI); library(RPostgres); library(dplyr)})
set.seed(20260629L)
OUT <- "/data/小米云盘/文档/论文/血磷/TTE/v28.4/output/pure_control"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

con <- dbConnect(Postgres(), host="localhost", port=5432, dbname="mimiciv",
                 user="postgres", password=Sys.getenv("MIMIC_PG_PASSWORD"))
co <- dbGetQuery(con, "SELECT * FROM public.pure_control_cohort")
dbDisconnect(con)

DAG14 <- c("age","gender","baseline_sofa","sepsis3","charlson","egfr",
           "vent_dim3","vasoactive_at_t0","admission_type","baseline_phosphate",
           "lactate_t0","refeeding_risk","corrected_calcium_mmol_min","nutrition_support_type")

co <- co %>% mutate(
  sexM = as.integer(gender=="M"),
  sepsis3 = as.integer(sepsis3),
  vent_dim3 = as.integer(vent_dim3),
  vasoactive_at_t0 = as.integer(vasoactive_at_t0),
  refeeding_risk = as.integer(refeeding_risk),
  rrt_at_t0 = as.integer(rrt_at_t0),
  admission_type = as.character(admission_type),
  nutrition_support_type = as.character(nutrition_support_type)
)

# ---- SMD function (continuous + categorical multi-level) ----
smd_one <- function(d, var, grp) {
  x1 <- d[[var]][grp==1]; x0 <- d[[var]][grp==0]
  if (is.numeric(x1)) {
    if (sd(c(x1,x0), na.rm=TRUE)==0) return(0)
    (mean(x1,na.rm=TRUE)-mean(x0,na.rm=TRUE)) /
      sqrt((var(x1,na.rm=TRUE)+var(x0,na.rm=TRUE))/2)
  } else {
    lv <- na.omit(unique(c(x1,x0)))
    s1 <- sapply(lv, function(l) mean(x1==l,na.rm=TRUE))
    s0 <- sapply(lv, function(l) mean(x0==l,na.rm=TRUE))
    # pooled covariance proxy via overall proportions
    pp <- (s1*length(na.omit(x1)) + s0*length(na.omit(x0)))/(length(na.omit(x1))+length(na.omit(x0)))
    sqrt(sum((s1-s0)^2 / pmax(pp*(1-pp), 1e-6)) / max(1,length(lv)-1))
  }
}

fit_comp <- function(d, label) {
  cc <- d %>% filter(complete.cases(across(all_of(DAG14))))
  cc$A <- as.integer(cc$grp==1)
  # PS
  fml <- as.formula(paste("A ~", paste(DAG14, collapse="+")))
  ps <- glm(fml, data=cc, family=binomial())
  e <- predict(ps, type="response")
  pA <- mean(cc$A)
  w <- ifelse(cc$A==1, pA/e, (1-pA)/(1-e))
  qs <- quantile(w, c(.01,.99))
  wt <- pmin(pmax(w, qs[1]), qs[2])
  cc$ps <- e; cc$w <- wt

  bal <- lapply(DAG14, function(v)
    data.frame(comparison=label, variable=v,
               smd_unweighted=smd_one(cc, v, cc$A),
               smd_weighted=smd_one_w(cc, v, wt))) %>% bind_rows()

  list(data=cc, bal=bal,
       wsum=data.frame(comparison=label, n=nrow(cc), n_trt=sum(cc$A), n_ctrl=sum(cc$A==0),
                       events_trt=NA, events_ctrl=NA,
                       p1=qs[1], p99=qs[2], w_median=median(w), w_max=max(w),
                       pct_truncated=100*mean(w!=wt), max_abs_smd_unw=max(abs(bal$smd_unweighted)),
                       max_abs_smd_w=max(abs(bal$smd_weighted))))
}
# weighted SMD via weighted means (continuous) / props (categorical)
smd_one_w <- function(d, var, w) {
  x <- d[[var]]; A <- d$A
  if (is.numeric(x)) {
    m <- function(ii) weighted.mean(x[ii], w[ii], na.rm=TRUE)
    v <- function(ii) cov.wt(data.frame(x[ii]), w[ii])$cov[1]
    (m(A==1)-m(A==0))/sqrt((v(A==1)+v(A==0))/2)
  } else {
    lv <- na.omit(unique(x))
    wp <- function(ii) sapply(lv, function(l) weighted.mean(x[ii]==l, w[ii], na.rm=TRUE))
    s1 <- wp(A==1); s0 <- wp(A==0)
    pp <- (s1+s0)/2
    sqrt(sum((s1-s0)^2/pmax(pp*(1-pp),1e-6))/max(1,length(lv)-1))
  }
}

c1 <- co %>% filter(exposure_w48 %in% c("A_IV","D_none")) %>% mutate(grp=ifelse(exposure_w48=="A_IV",1,0))
c3 <- co %>% filter(exposure_w48 %in% c("B_PO_only","D_none")) %>% mutate(grp=ifelse(exposure_w48=="B_PO_only",1,0))

r1 <- fit_comp(c1, "C1_IV_vs_none")
r3 <- fit_comp(c3, "C3_PO_vs_none")

bal <- bind_rows(r1$bal, r3$bal)
wsum <- bind_rows(r1$wsum, r3$wsum)
write.csv(bal, file.path(OUT,"pc_balance_smd.csv"), row.names=FALSE)
write.csv(wsum, file.path(OUT,"pc_weights_summary.csv"), row.names=FALSE)
saveRDS(list(C1=r1$data, C3=r3$data), file.path(OUT,"pc_weighted_cohorts.rds"))

# ---- flowchart ----
fc <- data.frame(stage=c("Eligible ICU stays (study_cohort)",
                         "  A_IV (actual IV infusion W48)",
                         "  B_PO_only (oral W48)",
                         "  C_outside_window",
                         "  C_rx_only_uninfused",
                         "  D_none (no supplement ever)",
                         "C1 analysis: A_IV + D_none",
                         "  DAG14 complete cases",
                         "C3 analysis: B_PO + D_none",
                         "  DAG14 complete cases"),
                 n=c(nrow(co), sum(co$exposure_w48=="A_IV"), sum(co$exposure_w48=="B_PO_only"),
                     sum(co$exposure_w48=="C_outside_window"), sum(co$exposure_w48=="C_rx_only_uninfused"),
                     sum(co$exposure_w48=="D_none"),
                     nrow(c1), nrow(r1$data), nrow(c3), nrow(r3$data)))
write.csv(fc, file.path(OUT,"pc_flowchart.csv"), row.names=FALSE)

cat("=== FLOWCHART ===\n"); print(fc)
cat("\n=== WEIGHTS ===\n"); print(wsum)
cat("\n=== BALANCE (|SMD|>0.1 unweighted / weighted) ===\n")
print(bal %>% group_by(comparison) %>% summarise(n_unw_gt01=sum(abs(smd_unweighted)>0.1),
                                                  n_w_gt01=sum(abs(smd_weighted)>0.1),
                                                  max_unw=round(max(abs(smd_unweighted)),3),
                                                  max_w=round(max(abs(smd_weighted)),3)))
cat("\nSaved:", OUT, "\n")
