#!/usr/bin/env Rscript
# =============================================================================
# STE-RMST (v28.4) — Limozin, Seaman, Su (2025)
#   * trial-period 级事件 (n=32,330)，不折叠
#   * day-specific pooled logistic regression + arm×day + DAG14 协变量
#   * g-computation → E[S(t|do(A=a))] 人群平均 → RMST_28
#   * subject 级 bootstrap B=200（同一 subject 所有 trial-period 一起重抽）
#   * 同时输出 raw unweighted + stay-unique 对照
# 须 LD_LIBRARY_PATH=/opt/miniforge3/lib
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(survival)
})
BASE <- "/data/小米云盘/文档/论文/血磷"
AUG_F <- file.path(BASE,"TTE/v26/output/sequential_trials_augmented.parquet")
MICE_F <- file.path(BASE,"sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet")
OUT <- file.path(BASE,"TTE/v28.4/output/tables/fig2")
TAU <- 28L; SEED <- 20260629L
BIN <- c("sepsis3","invasive_vent_at_t0","vasoactive_at_t0","refeeding_risk")
CONT<- c("age","baseline_sofa","charlson","egfr","baseline_phosphate",
         "lactate_t0","corrected_calcium_mmol_min")
FACT<- c("gender","admission_type","nutrition_support_type")
LEV <- list(gender=c("F","M"), admission_type=c("Observation","Emergency","Surgical_Elective"),
            nutrition_support_type=c("none","enteral","parenteral","mixed"))
fmt <- function(x) formatC(x, format="d", big.mark=",")

# ========== 1. 数据准备 ==========
cat("===== 1. Trial-period 级数据 =====\n")
aug <- as.data.table(read_parquet(AUG_F))
imp <- as.data.table(read_parquet(sprintf(MICE_F, 1)))
p <- aug[, .(stay_id, subject_id, trial_num,
             arm  = as.integer(as.logical(iv_within_12h)),
             death= as.integer(death_28day_num),
             death_t = as.numeric(death_time_from_t0))]
for(v in c(CONT,BIN)) p[,(v):=as.numeric(imp[[v]])]
for(v in FACT) p[,(v):=factor(aug[[v]], levels=LEV[[v]])]
p[,(BIN):=lapply(.SD, as.integer), .SDcols=BIN]
p[, death_day := fifelse(death==1 & !is.na(death_t),
                         pmin(TAU, pmax(0L, as.integer(ceiling(death_t)))), TAU+1L)]
cat(sprintf("  trial-periods n = %s  stays = %s  subjects = %s\n",
            fmt(nrow(p)), fmt(uniqueN(p$stay_id)), fmt(uniqueN(p$subject_id))))
cat(sprintf("  IV=%s  Ctrl=%s  deaths=%s  IV%%=%.1f\n",
            fmt(p[arm==1,.N]), fmt(p[arm==0,.N]), fmt(p[death==1,.N]),
            p[arm==1,.N]/nrow(p)*100))

# stay-unique 对照（取每个 stay 的 first trial）
pu <- p[order(stay_id,trial_num), .SD[1], by=stay_id]

# ========== 2. Raw unweighted RMST (trial-period + stay-unique) ==========
cat("\n===== 2. Raw unweighted RMST =====\n")
raw_rmst <- function(pdat, tag){
  risk0 <- sapply(0:TAU, function(t) mean(pdat[arm==0, death_day <= t], na.rm=TRUE))
  risk1 <- sapply(0:TAU, function(t) mean(pdat[arm==1, death_day <= t], na.rm=TRUE))
  S0 <- 1 - risk0; S1 <- 1 - risk1
  rmst0 <- sum(S0[1:TAU]); rmst1 <- sum(S1[1:TAU])
  cat(sprintf("  %s: RMST Ctrl=%.3f IV=%.3f diff=%+.3f | risk@28 Ctrl=%.3f IV=%.3f RD=%+.3f\n",
              tag, rmst0, rmst1, rmst1-rmst0, risk0[TAU+1], risk1[TAU+1], risk1[TAU+1]-risk0[TAU+1]))
  data.table(tag, rmst0, rmst1, diff=rmst1-rmst0, risk0=risk0[TAU+1], risk1=risk1[TAU+1], rd=risk1[TAU+1]-risk0[TAU+1])
}
res_raw_tp  <- raw_rmst(p,   "trial-period raw")
res_raw_su  <- raw_rmst(pu,  "stay-unique raw")

# ========== 3. IPTW-stabilized weighted RMST (day-level, 避免过度重复计数) ==========
cat("\n===== 3. IPTW stabilized weighted RMST =====\n")
ps_fml <- as.formula(paste("arm ~", paste(c(CONT,BIN,FACT), collapse=" + ")))
ps_fit <- suppressWarnings(glm(ps_fml, data=p, family=binomial(), control=list(maxit=200)))
p[, ps := predict(ps_fit, type="response")]
p[, pa := mean(arm)]
p[, w := fifelse(arm==1, pa/ps, (1-pa)/(1-ps))]
p[, w := pmin(pmax(w, quantile(w,0.01,na.rm=TRUE)), quantile(w,0.99,na.rm=TRUE))]
cat(sprintf("  PS range: [%.3f, %.3f]  W range: [%.2f, %.2f]\n",
            range(p$ps,na.rm=TRUE)[1], range(p$ps,na.rm=TRUE)[2],
            quantile(p$w,0.01,na.rm=TRUE), quantile(p$w,0.99,na.rm=TRUE)))
# day-level weighted cumulative risk
for(t in 0:TAU) p[, sprintf("d%d",t) := as.integer(death_day <= t)]
wrmst <- function(pp, tag){
  tot0 <- sum(pp[arm==0, w]); tot1 <- sum(pp[arm==1, w])
  w0 <- sapply(0:TAU, function(t) sum(pp[arm==0][[sprintf("d%d",t)]] * pp[arm==0, w]) / tot0)
  w1 <- sapply(0:TAU, function(t) sum(pp[arm==1][[sprintf("d%d",t)]] * pp[arm==1, w]) / tot1)
  S0 <- 1-w0; S1 <- 1-w1
  rmst0 <- sum(S0[1:TAU]); rmst1 <- sum(S1[1:TAU])
  cat(sprintf("  %s: wRMST Ctrl=%.3f IV=%.3f diff=%+.3f\n", tag, rmst0, rmst1, rmst1-rmst0))
  data.table(tag, rmst0, rmst1, diff=rmst1-rmst1)
}
res_iptw_tp <- wrmst(p,  "trial-period IPTW")
# stay-unique 的 IPTW (权重只在 stay 级做一次)
pu[, ps_u := predict(ps_fit, newdata=pu, type="response")]
pu[, w_u  := fifelse(arm==1, mean(pu$arm)/ps_u, (1-mean(pu$arm))/(1-ps_u))]
pu[, w_u  := pmin(pmax(w_u, quantile(w_u,0.01,na.rm=TRUE)), quantile(w_u,0.99,na.rm=TRUE))]
# 手动算 stay-unique 加权
tot0u <- sum(pu[arm==0, w_u]); tot1u <- sum(pu[arm==1, w_u])
for(t in 0:TAU) pu[, sprintf("d%d",t) := as.integer(death_day <= t)]
w0u <- sapply(0:TAU, function(t) sum(pu[arm==0][[sprintf("d%d",t)]] * pu[arm==0, w_u]) / tot0u)
w1u <- sapply(0:TAU, function(t) sum(pu[arm==1][[sprintf("d%d",t)]] * pu[arm==1, w_u]) / tot1u)
res_iptw_su <- data.table(tag="stay-unique IPTW", rmst0=sum(1-w0u[1:TAU]),
                           rmst1=sum(1-w1u[1:TAU]), diff=sum(1-w1u[1:TAU])-sum(1-w0u[1:TAU]))
cat(sprintf("  stay-unique IPTW: RMST Ctrl=%.3f IV=%.3f diff=%+.3f\n",
            res_iptw_su$rmst0, res_iptw_su$rmst1, res_iptw_su$diff))

# ========== 4. STE g-computation (trial-period 级) ==========
cat("\n===== 4. STE g-computation (Limozin 2025) =====\n")

# 加权 long format (每个 trial-period × 28 天)
p[, idx := seq_len(.N)]
lng <- merge(CJ(idx = p$idx, day = 0:(TAU-1)),
             p[, c("idx","stay_id","subject_id","arm","death_day",CONT,BIN,FACT), with=FALSE],
             by="idx", allow.cartesian=TRUE)
lng[, Y := as.integer(death_day == day)]
cat(sprintf("  person-day rows: %s\n", fmt(nrow(lng))))

# 拟合 day-specific pooled logistic regression with arm×day
fml <- as.formula(paste("Y ~ arm + day + arm:day +",
                        paste(c(CONT,BIN,FACT), collapse=" + ")))
set.seed(SEED)
fit <- suppressWarnings(glm(fml, data=lng, family=binomial(), control=list(maxit=200)))
cat("  Model fit OK\n")
coef(fit)[1:8] |> print()

# g-computation
pred_grid <- CJ(idx = seq_len(nrow(p)), day = 0:(TAU-1), arm = c(0L,1L))
pred_grid <- merge(pred_grid,
                   p[, c("idx", CONT, BIN, FACT), with=FALSE],
                   by="idx")
pred_grid[, P := predict(fit, newdata=pred_grid, type="response")]
pred_grid <- pred_grid[order(idx, arm, day)]
pred_grid[, S := cumprod(1 - P), by=.(idx, arm)]
E_S <- pred_grid[, .(E_S = mean(S, na.rm=TRUE)), by=.(day, arm)]
E_S_wide <- dcast(E_S, day ~ arm, value.var="E_S")
setnames(E_S_wide, c("day","E_S0","E_S1"))
rmst0 <- sum(E_S_wide$E_S0); rmst1 <- sum(E_S_wide$E_S1)
gcomp_diff <- rmst1 - rmst0
cat(sprintf("  STE g-comp RMST Ctrl = %.3f d\n", rmst0))
cat(sprintf("  STE g-comp RMST IV   = %.3f d\n", rmst1))
cat(sprintf("  STE g-comp diff      = %+.3f d\n", gcomp_diff))
cat(sprintf("  STE g-comp risk@28   Ctrl=%.3f IV=%.3f\n",
            1 - E_S_wide[day==TAU-1]$E_S0, 1 - E_S_wide[day==TAU-1]$E_S1))

# ========== 5. Subject-level bootstrap CI (B=200) ==========
cat("\n===== 5. Subject-level bootstrap CI (B=200) =====\n")
B <- 200L
set.seed(SEED)
subs <- unique(lng$subject_id)
bm <- matrix(NA_real_, B, 3); nb <- 0L
colnames(bm) <- c("rmst0","rmst1","diff")

while(nb < B){
  pick <- sample(subs, length(subs), replace=TRUE)
  blng <- lng[subject_id %in% pick]
  # 同一 subject 可能被多次抽入 → 去重 trial-period
  blng <- unique(blng, by=c("subject_id","idx","day"))
  if(nrow(blng) < 15000) next
  f <- tryCatch(suppressWarnings(
    glm(fml, data=blng, family=binomial(), control=list(maxit=200))),
    error=function(e) NULL)
  if(is.null(f)) next
  pg <- CJ(idx = unique(blng$idx), day = 0:(TAU-1), arm = c(0L,1L))
  pg <- merge(pg, p[idx %in% unique(blng$idx), c("idx",CONT,BIN,FACT), with=FALSE], by="idx")
  pg[, P := suppressWarnings(predict(f, newdata=pg, type="response"))]
  pg <- pg[order(idx, arm, day)]
  pg[, S := cumprod(1 - P), by=.(idx, arm)]
  ES <- pg[, .(E_S=mean(S,na.rm=TRUE)), by=.(day,arm)]
  ESw <- dcast(ES, day ~ arm, value.var="E_S")
  r0 <- sum(ESw$`0`, na.rm=TRUE); r1 <- sum(ESw$`1`, na.rm=TRUE)
  if(any(!is.finite(c(r0,r1)))) next
  nb <- nb + 1L; bm[nb,] <- c(r0, r1, r1-r0)
  if(nb %% 50 == 0) cat(sprintf("  boot %d/%d (rmst diff range %+.2f)\n",
                                nb, B, bm[1:nb,3] |> range()), "\n")
}
ci <- quantile(bm[1:nb,3], c(0.025, 0.975), na.rm=TRUE)
se <- sd(bm[1:nb,3]); p_val <- 2*pnorm(-abs(gcomp_diff/se))
cat(sprintf("  Bootstrap CI: [%+.3f, %+.3f]  SE=%.3f  p=%.3f  (n_boot=%d)\n",
            ci[1], ci[2], se, p_val, nb))

# ========== 6. 汇总输出 ==========
cat("\n===== 6. 汇总结果 =====\n")
result <- data.table(
  method = c("trial-period raw unweighted",
             "stay-unique raw unweighted",
             "trial-period IPTW stabilized",
             "stay-unique IPTW stabilized",
             "** STE g-computation (Limozin 2025) **"),
  event_unit = c("trial-period","stay","trial-period","stay","trial-period"),
  n_unit = c(nrow(p), nrow(pu), nrow(p), nrow(pu), nrow(p)),
  rmst0 = c(res_raw_tp$rmst0, res_raw_su$rmst0,
            res_iptw_tp$rmst0, res_iptw_su$rmst0, rmst0),
  rmst1 = c(res_raw_tp$rmst1, res_raw_su$rmst1,
            res_iptw_tp$rmst1, res_iptw_su$rmst1, rmst1),
  diff = c(res_raw_tp$diff, res_raw_su$diff,
           res_iptw_tp$diff, res_iptw_su$diff, gcomp_diff),
  ci_lower = c(NA, NA, NA, NA, ci[1]),
  ci_upper = c(NA, NA, NA, NA, ci[2]),
  p_value  = c(NA, NA, NA, NA, p_val)
)
print(result, row.names=FALSE)

fwrite(result, file.path(OUT, "ste_rmst_comparison.csv"), row.names=FALSE)
fwrite(E_S_wide, file.path(OUT, "ste_rmst_survival_curve.csv"), row.names=FALSE)
fwrite(data.table(bm[1:nb,]), file.path(OUT, "ste_rmst_boot.csv"), row.names=FALSE)

sink(file.path(OUT,"ste_rmst_model.txt"))
cat("=== STE g-computation pooled logistic regression ===\n")
cat(sprintf("Model formula: %s\n", deparse(fml)))
print(summary(fit))
cat("\n\n=== ANOVA (LRT vs null) ===\n")
print(anova(fit, test="LRT"))
sink()

cat(sprintf("\n✓ 输出:\n  %s\n  %s\n  %s\n  %s\n",
            file.path(OUT,"ste_rmst_comparison.csv"),
            file.path(OUT,"ste_rmst_survival_curve.csv"),
            file.path(OUT,"ste_rmst_boot.csv"),
            file.path(OUT,"ste_rmst_model.txt")))
cat("\n===== 完成 =====\n")
