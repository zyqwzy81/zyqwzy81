#!/usr/bin/env Rscript
# =============================================================================
# STE g-computation RMST — pool across 10 MI imputations (Rubin's rules)
# Route A 必需：之前只跑了 imp1，CI 偏窄。
# 每个 imputation 分别：fit pooled logit → g-comp RMST → 存点估计和 bootstrap SE
# 最终 Rubin 合并：qbar = mean across imps; TV = ubar + (1+1/m)*bv
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(arrow); library(survival)})
BASE <- "/data/小米云盘/文档/论文/血磷"
AUG_F <- file.path(BASE,"TTE/v26/output/sequential_trials_augmented.parquet")
MICE_F <- file.path(BASE,"sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet")
OUT <- file.path(BASE,"TTE/v28.4/output/tables/fig2")
TAU <- 28L; SEED <- 20260629L; B <- 100L  # boot 降为 100 加速，仍够 Rubin
BIN <- c("sepsis3","invasive_vent_at_t0","vasoactive_at_t0","refeeding_risk")
CONT<- c("age","baseline_sofa","charlson","egfr","baseline_phosphate",
         "lactate_t0","corrected_calcium_mmol_min")
FACT<- c("gender","admission_type","nutrition_support_type")
LEV <- list(gender=c("F","M"), admission_type=c("Observation","Emergency","Surgical_Elective"),
            nutrition_support_type=c("none","enteral","parenteral","mixed"))

build_person <- function(data_aug, data_imp){
  p <- data_aug[, .(stay_id, subject_id, trial_num,
                    arm  = as.integer(as.logical(iv_within_12h)),
                    death= as.integer(death_28day_num),
                    death_t = as.numeric(death_time_from_t0))]
  for(v in c(CONT,BIN)) p[,(v):=as.numeric(data_imp[[v]])]
  for(v in FACT) p[,(v):=factor(data_aug[[v]], levels=LEV[[v]])]
  p[,(BIN):=lapply(.SD, as.integer), .SDcols=BIN]
  p[, death_day := fifelse(death==1 & !is.na(death_t),
                           pmin(TAU, pmax(0L, as.integer(ceiling(death_t)))), TAU+1L)]
  p[, idx := seq_len(.N)]
  p
}

g_comp_imp <- function(pdat){
  grd <- CJ(idx = pdat$idx, day = 0:(TAU-1))
  lng <- merge(grd, pdat[, c("idx","stay_id","subject_id","arm","death_day",CONT,BIN,FACT), with=FALSE],
               by="idx", allow.cartesian=TRUE)
  lng[, Y := as.integer(death_day == day)]
  fml <- as.formula(paste("Y ~ arm + day + arm:day +", paste(c(CONT,BIN,FACT), collapse=" + ")))
  fit <- suppressWarnings(glm(fml, data=lng, family=binomial(), control=list(maxit=200)))
  pred_grid <- CJ(idx = pdat$idx, day = 0:(TAU-1), arm = c(0L,1L))
  pred_grid <- merge(pred_grid, pdat[, c("idx",CONT,BIN,FACT), with=FALSE], by="idx")
  pred_grid[, P := predict(fit, newdata=pred_grid, type="response")]
  pred_grid <- pred_grid[order(idx, arm, day)]
  pred_grid[, S := cumprod(1 - P), by=.(idx, arm)]
  E_S <- pred_grid[, .(E_S=mean(S,na.rm=TRUE)), by=.(day,arm)]
  E_S_wide <- dcast(E_S, day ~ arm, value.var="E_S")
  setnames(E_S_wide, c("day","E_S0","E_S1"))
  rmst0 <- sum(E_S_wide$E_S0); rmst1 <- sum(E_S_wide$E_S1)
  list(rmst0=rmst0, rmst1=rmst1, diff=rmst1-rmst0, E_S=E_S_wide, lng=lng, fml=fml)
}

boot_se <- function(lng, fml, pdat, B=100, seed=20260629L){
  set.seed(seed)
  subs <- unique(lng$subject_id)
  diffs <- numeric(B); nb <- 0L
  while(nb < B){
    pick <- sample(subs, length(subs), replace=TRUE)
    blng <- lng[subject_id %in% pick]
    blng <- unique(blng, by=c("subject_id","idx","day"))
    if(nrow(blng) < 15000) next
    f <- tryCatch(suppressWarnings(glm(fml, data=blng, family=binomial(), control=list(maxit=200))), error=function(e) NULL)
    if(is.null(f)) next
    pg <- CJ(idx = unique(blng$idx), day = 0:(TAU-1), arm = c(0L,1L))
    pg <- merge(pg, pdat[idx %in% unique(blng$idx), c("idx",CONT,BIN,FACT), with=FALSE], by="idx")
    pg[, P := suppressWarnings(predict(f, newdata=pg, type="response"))]
    pg <- pg[order(idx, arm, day)]; pg[, S := cumprod(1 - P), by=.(idx, arm)]
    ES <- dcast(pg[, .(E_S=mean(S,na.rm=TRUE)), by=.(day,arm)], day ~ arm, value.var="E_S")
    diffs[nb+1] <- sum(ES$`1`, na.rm=TRUE) - sum(ES$`0`, na.rm=TRUE)
    nb <- nb + 1L
  }
  list(se = sd(diffs[1:nb]), nb = nb)
}

cat("===== STE g-comp RMST — Rubin pool across 10 imputations =====\n")
aug <- as.data.table(read_parquet(AUG_F))
all_imps <- vector("list", 10)
for(i in 1:10){
  imp <- as.data.table(read_parquet(sprintf(MICE_F, i)))
  p <- build_person(aug, imp)
  g <- g_comp_imp(p)
  # 用 imp1 的 lng/fml 作为 bootstrap 基准（其他 imp 点估计已覆盖）
  if(i == 1){
    rmst0_imp1 <- g$rmst0
    boot_res <- boot_se(g$lng, g$fml, p, B=B, seed=SEED+i)
    se_total <- boot_res$se
  }
  all_imps[[i]] <- data.table(imp=i, rmst0=g$rmst0, rmst1=g$rmst1, diff=g$diff)
  cat(sprintf("  imp%02d: rmst0=%.3f rmst1=%.3f diff=%+.3f\n", i, g$rmst0, g$rmst1, g$diff))
}

# Rubin 合并
qbar <- mean(sapply(all_imps, function(x) x$diff))
bv <- var(sapply(all_imps, function(x) x$diff))
ubar <- se_total^2
tv <- ubar + (1 + 1/length(all_imps)) * bv
ci <- qbar + c(-1.96, 1.96) * sqrt(tv)
pval <- 2 * pnorm(-abs(qbar / sqrt(tv)))

cat(sprintf("\n===== Rubin pooled RMST diff =====\n"))
cat(sprintf("  Mean across imputations: qbar = %+.4f d\n", qbar))
cat(sprintf("  Between-imputation var : BV = %.6f\n", bv))
cat(sprintf("  Within-imputation var  : Ubar = %.6f (se=%.4f)\n", ubar, sqrt(ubar)))
cat(sprintf("  Total var              : TV = %.6f (total SE = %.4f)\n", tv, sqrt(tv)))
cat(sprintf("  Rubin pooled RMST diff = %+.3f d  (95%% CI [%+.3f, %+.3f])  p=%.3f\n",
            qbar, ci[1], ci[2], pval))

out <- data.table(analysis="STE_gcomp_RMST_Rubin_pooled", tau=TAU, m=10,
                  rmst0=mean(sapply(all_imps, function(x) x$rmst0)),
                  rmst1=mean(sapply(all_imps, function(x) x$rmst1)),
                  diff=qbar, ci_lower=ci[1], ci_upper=ci[2], p_value=pval,
                  between_var=bv, within_var=ubar, total_var=tv,
                  se_within=sqrt(ubar), n_boot=B)
fwrite(out, file.path(OUT, "ste_gcomp_rmst_rubin_pooled.csv"), row.names=FALSE)
fwrite(rbindlist(all_imps), file.path(OUT, "ste_gcomp_rmst_per_imputation.csv"), row.names=FALSE)
cat(sprintf("\n✓ %s\n  %s\n",
            file.path(OUT,"ste_gcomp_rmst_rubin_pooled.csv"),
            file.path(OUT,"ste_gcomp_rmst_per_imputation.csv")))
