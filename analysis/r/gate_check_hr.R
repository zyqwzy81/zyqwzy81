#!/usr/bin/env Rscript
# =============================================================================
# GATE-CHECK: Trial-period level IPTW-Cox HR with subject-clustered robust SE
# Route A 的关键前置验证：如果 trial-period HR 与 stay-unique HR=0.982 差太大，
# 说明 estimand 切换导致效应方向变化，需要重新评估路线
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(arrow); library(survival); library(sandwich)})
BASE <- "/data/小米云盘/文档/论文/血磷"
AUG_F <- file.path(BASE,"TTE/v26/output/sequential_trials_augmented.parquet")
MICE_F <- file.path(BASE,"sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet")
OUT <- file.path(BASE,"TTE/v28.4/output/tables/main")
TAU <- 28L; SEED <- 20260629L
BIN <- c("sepsis3","invasive_vent_at_t0","vasoactive_at_t0","refeeding_risk")
CONT<- c("age","baseline_sofa","charlson","egfr","baseline_phosphate",
         "lactate_t0","corrected_calcium_mmol_min")
FACT<- c("gender","admission_type","nutrition_support_type")
LEV <- list(gender=c("F","M"), admission_type=c("Observation","Emergency","Surgical_Elective"),
            nutrition_support_type=c("none","enteral","parenteral","mixed"))

fmt <- function(x) formatC(x, format="d", big.mark=",")
ci_hr <- function(coef, se) exp(c(coef + c(-1.96,1.96)*se))

cat("===== Gate-check: Trial-period vs stay-unique IPTW-Cox =====\n\n")
set.seed(SEED)
aug <- as.data.table(read_parquet(AUG_F))
imp <- as.data.table(read_parquet(sprintf(MICE_F, 1)))

build_person <- function(data_aug, data_imp, fold=FALSE, first_only=FALSE){
  # fold: TRUE → stay-unique (.SD[1] by stay_id)
  p <- data_aug[, .(stay_id, subject_id, trial_num, arm = as.integer(as.logical(iv_within_12h)),
                    death = as.integer(death_28day_num), death_t = as.numeric(death_time_from_t0))]
  for(v in c(CONT,BIN)) p[,(v):=as.numeric(data_imp[[v]])]
  for(v in FACT) p[,(v):=factor(data_aug[[v]], levels=LEV[[v]])]
  p[,(BIN):=lapply(.SD, as.integer), .SDcols=BIN]
  p[, surv_time := fifelse(is.na(death_t), TAU, pmin(TAU, death_t))]
  if(first_only) p <- p[order(stay_id, trial_num), .SD[1], by=stay_id]
  p[, arm_f := factor(arm, labels=c("Ctrl","IV"))]
  p
}

# ------ trial-period level ------
ptp <- build_person(aug, imp, first_only=FALSE)
cat(sprintf("Trial-period level: n = %s  IV = %s (%.1f%%)  deaths = %s  subjects = %s\n",
            fmt(nrow(ptp)), fmt(ptp[arm==1,.N]), ptp[arm==1,.N]/nrow(ptp)*100,
            fmt(ptp[death==1,.N]), fmt(uniqueN(ptp$subject_id))))

# ------ stay-unique level (对照) ------
psu <- build_person(aug, imp, first_only=TRUE)
cat(sprintf("Stay-unique level:  n = %s  IV = %s (%.1f%%)  deaths = %s\n\n",
            fmt(nrow(psu)), fmt(psu[arm==1,.N]), psu[arm==1,.N]/nrow(psu)*100, fmt(psu[death==1,.N])))

run_iptw_cox <- function(pdat, tag){
  ps_fml <- as.formula(paste("arm ~", paste(c(CONT,BIN,FACT), collapse=" + ")))
  ps <- suppressWarnings(predict(glm(ps_fml, data=pdat, family=binomial()), type="response"))
  pa <- mean(pdat$arm)
  w <- fifelse(pdat$arm==1, pa/ps, (1-pa)/(1-ps))
  w <- pmin(pmax(w, quantile(w, 0.01, na.rm=TRUE)), quantile(w, 0.99, na.rm=TRUE))
  cat(sprintf("  [%s] PS=[%.3f,%.3f]  W=[%.2f,%.2f]\n",
              tag, range(ps,na.rm=TRUE)[1], range(ps,na.rm=TRUE)[2],
              quantile(w,0.01,na.rm=TRUE), quantile(w,0.99,na.rm=TRUE)))
  pdat[, w := w]

  # Unadjusted weighted Cox
  cox <- coxph(Surv(surv_time, death) ~ arm_f, weights=w, cluster=subject_id, robust=TRUE, data=pdat)
  sm <- summary(cox)
  coef_row <- sm$coefficients[1,]
  hr_unadj <- exp(coef_row["coef"])
  se_robust_unadj <- coef_row["robust se"]
  ci_unadj <- exp(coef_row["coef"] + c(-1.96,1.96)*se_robust_unadj)
  p_unadj <- coef_row["Pr(>|z|)"]
  cat(sprintf("  [%s] IPTW-Cox unadj HR = %.3f  (%.3f, %.3f)  p=%.3f\n",
              tag, hr_unadj, ci_unadj[1], ci_unadj[2], p_unadj))

  # Fully-adjusted IPTW-Cox (DAG14 as covariates + weights + cluster SE)
  fml_a <- as.formula(paste("Surv(surv_time, death) ~ arm_f +", paste(c(CONT,BIN,FACT), collapse=" + ")))
  cox_a <- coxph(fml_a, weights=w, cluster=subject_id, robust=TRUE, data=pdat)
  sm_a <- summary(cox_a)
  coef_a <- sm_a$coefficients["arm_fIV", ]
  hr_adj <- exp(coef_a["coef"])
  se_robust_adj <- coef_a["robust se"]
  ci_adj <- exp(coef_a["coef"] + c(-1.96,1.96)*se_robust_adj)
  p_adj <- coef_a["Pr(>|z|)"]
  cat(sprintf("  [%s] Fully-adj IPTW-Cox HR = %.3f  (%.3f, %.3f)  p=%.3f\n\n",
              tag, hr_adj, ci_adj[1], ci_adj[2], p_adj))

  data.table(tag, n=nrow(pdat),
             hr_unadj, ci_lo_unadj=ci_unadj[1], ci_hi_unadj=ci_unadj[2], p_unadj,
             hr_adj,   ci_lo_adj=ci_adj[1],   ci_hi_adj=ci_adj[2],   p_adj)
}

res_tp <- run_iptw_cox(ptp, "trial-period")
res_su <- run_iptw_cox(psu, "stay-unique")

res <- rbind(res_tp, res_su)
fwrite(res, file.path(OUT, "gate_check_tp_vs_su_hr.csv"), row.names=FALSE)
cat(res)
cat(sprintf("\n===== Gate check done =====\n"))
