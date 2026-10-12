#!/usr/bin/env Rscript
# =============================================================================
# STE g-computation survival curve with bootstrap CI band
#   * 读取 v28.4_ste_rmst_gcomp.R 已经拟合的模型
#   * 提取每天 E[S(t|do(A=0))] / E[S(t|do(A=1))]
#   * Bootstrap 200 次得到每天 CI band
#   * ggplot2 geom_step + geom_ribbon 绘制
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(ggplot2)
})
BASE <- "/data/小米云盘/文档/论文/血磷"
AUG_F <- file.path(BASE,"TTE/v26/output/sequential_trials_augmented.parquet")
MICE_F <- file.path(BASE,"sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet")
OUT <- file.path(BASE,"TTE/v28.4/output/figures")
TAU <- 28L; B <- 200L; SEED <- 20260629L
BIN <- c("sepsis3","invasive_vent_at_t0","vasoactive_at_t0","refeeding_risk")
CONT<- c("age","baseline_sofa","charlson","egfr","baseline_phosphate",
         "lactate_t0","corrected_calcium_mmol_min")
FACT<- c("gender","admission_type","nutrition_support_type")
LEV <- list(gender=c("F","M"), admission_type=c("Observation","Emergency","Surgical_Elective"),
            nutrition_support_type=c("none","enteral","parenteral","mixed"))
fmt <- function(x) formatC(x, format="d", big.mark=",")

# ========== 1. 加载数据 ==========
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
p[, idx := seq_len(.N)]
cat(sprintf("trial-periods n = %s  IV=%s  deaths=%s\n",
            fmt(nrow(p)), fmt(p[arm==1,.N]), fmt(p[death==1,.N])))

# ========== 2. 拟合 point-estimate pooled logit ==========
cat("\nFitting pooled logistic regression (day-specific)...\n")
fml <- as.formula(paste("Y ~ arm + day + arm:day +",
                        paste(c(CONT,BIN,FACT), collapse=" + ")))
grd <- CJ(idx = p$idx, day = 0:(TAU-1))
lng <- merge(grd,
             p[, c("idx","stay_id","subject_id","arm","death_day",CONT,BIN,FACT), with=FALSE],
             by="idx", allow.cartesian=TRUE)
lng[, Y := as.integer(death_day == day)]
fit <- suppressWarnings(glm(fml, data=lng, family=binomial(), control=list(maxit=200)))
cat("  model fit ok\n")

# g-computation: 两场景 counterfactual 生存
pred_grid <- CJ(idx = p$idx, day = 0:(TAU-1), arm = c(0L,1L))
pred_grid <- merge(pred_grid, p[, c("idx",CONT,BIN,FACT), with=FALSE], by="idx")
cat(sprintf("  pred_grid cols: %s\n", paste(colnames(pred_grid)[1:12], collapse=",")))
cat(sprintf("  arm type=%s levels=%s\n", class(pred_grid$arm), paste(unique(pred_grid$arm), collapse=",")))

pred_grid[, P := predict(fit, newdata=pred_grid, type="response")]
cat("  predict ok\n")
pred_grid <- pred_grid[order(idx, arm, day)]
pred_grid[, S := cumprod(1 - P), by=.(idx, arm)]
E_S_pt <- pred_grid[, .(E_S = mean(S, na.rm=TRUE)), by=.(day, arm)]
setorder(E_S_pt, arm, day)

# ========== 3. Bootstrap CI band (subject-level) ==========
cat(sprintf("Bootstrap CI band (B=%d, subject-level)...\n", B))
set.seed(SEED)
subs <- unique(lng$subject_id)
# B × 2 × 28 矩阵: [boot_id, arm, day]
boot_E_S <- array(NA_real_, dim=c(B, 2, TAU),
                  dimnames=list(boot=1:B, arm=c("0","1"), day=0:(TAU-1)))
nb <- 0L
while(nb < B){
  pick <- sample(subs, length(subs), replace=TRUE)
  blng <- lng[subject_id %in% pick]
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
  if(nrow(ES) != 56) next  # 28×2
  nb <- nb + 1L
  boot_E_S[nb, "0", ] <- ES[arm==0][order(day)]$E_S
  boot_E_S[nb, "1", ] <- ES[arm==1][order(day)]$E_S
  if(nb %% 50 == 0) cat(sprintf("  boot %d/%d\n", nb, B))
}
cat(sprintf("  有效 bootstrap = %d / %d\n", nb, B))

# 每天 CI 分位数
ci_band <- rbind(
  data.table(arm=0L, day=0:(TAU-1),
             E_S_pt = E_S_pt[arm==0][order(day)]$E_S,
             E_S_lo = apply(boot_E_S[1:nb, "0", ], 2, quantile, 0.025, na.rm=TRUE),
             E_S_hi = apply(boot_E_S[1:nb, "0", ], 2, quantile, 0.975, na.rm=TRUE)),
  data.table(arm=1L, day=0:(TAU-1),
             E_S_pt = E_S_pt[arm==1][order(day)]$E_S,
             E_S_lo = apply(boot_E_S[1:nb, "1", ], 2, quantile, 0.025, na.rm=TRUE),
             E_S_hi = apply(boot_E_S[1:nb, "1", ], 2, quantile, 0.975, na.rm=TRUE))
)

# ========== 4. Counterfactual survival curve ==========
cat("\n4. Counterfactual survival curve plotting...\n")
COL_IV   <- "#C0392B"
COL_CTRL <- "#2874A6"

p_curve <- ggplot(ci_band, aes(x=day, group=factor(arm, labels=c("Conservative management", "IV phosphate therapy")))) +
  # CI band
  geom_ribbon(aes(ymin=E_S_lo, ymax=E_S_hi,
                  fill=factor(arm, labels=c("Conservative management", "IV phosphate therapy"))),
              alpha=0.25, color=NA) +
  # Step lines (counterfactual survival)
  geom_step(aes(y=E_S_pt, color=factor(arm, labels=c("Conservative management", "IV phosphate therapy"))),
            linewidth=1.05, direction="hv") +
  # 28d mark
  geom_hline(yintercept=0, linetype="dashed", color="gray60", linewidth=0.4) +
  # RMST annotation
  annotate("text", x=2, y=0.04, size=3.6, color=COL_CTRL,
           label=sprintf("RMST_28 (Ctrl) = %.2f d",
                         sum(E_S_pt[arm==0][order(day)]$E_S))) +
  annotate("text", x=2, y=0.01, size=3.6, color=COL_IV,
           label=sprintf("RMST_28 (IV)   = %.2f d  diff = %+.2f d",
                         sum(E_S_pt[arm==1][order(day)]$E_S),
                         sum(E_S_pt[arm==1][order(day)]$E_S) - sum(E_S_pt[arm==0][order(day)]$E_S))) +
  scale_color_manual(values=c("Conservative management"=COL_CTRL, "IV phosphate therapy"=COL_IV),
                     name="Counterfactual\nstrategy") +
  scale_fill_manual(values=c("Conservative management"=COL_CTRL, "IV phosphate therapy"=COL_IV),
                    name="Counterfactual\nstrategy") +
  scale_x_continuous(limits=c(0, TAU), breaks=0:4*7, expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 1.005), breaks=seq(0, 1, 0.1), expand=c(0,0)) +
  labs(
    title="Marginal counterfactual survival probability",
    subtitle="STE g-computation (Limozin 2025) — trial-period level, n=32,330\n95% bootstrap CI band (subject-level, B=200)",
    x="Days since target trial origin",
    y=expression(paste("E[", italic(S), "(", italic(t), " | do(", italic(A), "=a))]"))
  ) +
  theme_bw(base_size=11) +
  theme(
    legend.position = c(0.62, 0.98),
    legend.justification = c(0.5, 1),
    legend.box.background = element_rect(color="black"),
    legend.title = element_text(size=10, face="bold"),
    plot.title = element_text(size=13, face="bold"),
    plot.subtitle = element_text(size=9.5, color="gray30"),
    axis.line = element_line(color="black"),
    panel.grid.minor = element_blank(),
    text = element_text(family="Helvetica")
  )

# ========== 5. Risk at t curve (complementary of survival) ==========
cat("5. Risk-at-t (cumulative mortality) plotting...\n")
ci_band[, risk_pt := 1 - E_S_pt]
ci_band[, risk_lo := 1 - E_S_hi]   # survival lo → risk hi
ci_band[, risk_hi := 1 - E_S_lo]   # survival hi → risk lo

p_risk <- ggplot(ci_band, aes(x=day, group=factor(arm, labels=c("Conservative management", "IV phosphate therapy")))) +
  geom_ribbon(aes(ymin=risk_lo, ymax=risk_hi,
                  fill=factor(arm, labels=c("Conservative management", "IV phosphate therapy"))),
              alpha=0.22, color=NA) +
  geom_step(aes(y=risk_pt, color=factor(arm, labels=c("Conservative management", "IV phosphate therapy"))),
            linewidth=1.05, direction="hv") +
  scale_color_manual(values=c("Conservative management"=COL_CTRL, "IV phosphate therapy"=COL_IV),
                     name="Strategy") +
  scale_fill_manual(values=c("Conservative management"=COL_CTRL, "IV phosphate therapy"=COL_IV),
                    name="Strategy") +
  scale_x_continuous(limits=c(0, TAU), breaks=0:4*7, expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 0.25), breaks=seq(0, 0.25, 0.05), expand=c(0,0),
                     labels=scales::percent_format(accuracy=1)) +
  labs(
    title="Counterfactual cumulative mortality risk",
    subtitle="STE g-computation — day-specific pooled logistic regression\nwith arm×time interaction + DAG14 baseline covariates",
    x="Days since target trial origin",
    y=expression(paste("Cumulative risk ", italic(F), "(", italic(t), ") = 1 - E[", italic(S), "(", italic(t), " | do(", italic(A), "=a))]"))
  ) +
  theme_bw(base_size=11) +
  theme(
    legend.position = c(0.62, 0.95),
    legend.justification = c(0.5, 1),
    legend.box.background = element_rect(color="black"),
    legend.title = element_text(size=10, face="bold"),
    plot.title = element_text(size=13, face="bold"),
    plot.subtitle = element_text(size=9.5, color="gray30"),
    axis.line = element_line(color="black"),
    panel.grid.minor = element_blank(),
    text = element_text(family="Helvetica")
  )

# ========== 6. Save ==========
pdf_w_mm <- 176; pdf_h_mm <- 140
ggsave(file.path(OUT, "Fig2_STE_counterfactual_survival.pdf"),
       p_curve, width=pdf_w_mm, height=pdf_h_mm, units="mm", dpi=600)
ggsave(file.path(OUT, "Fig2_STE_counterfactual_survival.png"),
       p_curve, width=pdf_w_mm, height=pdf_h_mm, units="mm", dpi=600)
ggsave(file.path(OUT, "Fig2_STE_cumulative_risk.pdf"),
       p_risk, width=pdf_w_mm, height=pdf_h_mm, units="mm", dpi=600)
ggsave(file.path(OUT, "Fig2_STE_cumulative_risk.png"),
       p_risk, width=pdf_w_mm, height=pdf_h_mm, units="mm", dpi=600)

# 保存完整 band 数据
fwrite(ci_band, file.path(BASE,"TTE/v28.4/output/tables/fig2/v28.4_ste_ci_band.csv"),
       row.names=FALSE)

cat(sprintf("\n✓ 输出:\n  %s\n  %s\n  %s\n  %s\n  %s\n",
            file.path(OUT,"Fig2_STE_counterfactual_survival.pdf"),
            file.path(OUT,"Fig2_STE_counterfactual_survival.png"),
            file.path(OUT,"Fig2_STE_cumulative_risk.pdf"),
            file.path(OUT,"Fig2_STE_cumulative_risk.png"),
            file.path(BASE,"TTE/v28.4/output/tables/fig2/v28.4_ste_ci_band.csv")))
cat("\n===== 完成 =====\n")
