# ============================================================================
# v28.4_stay_common.R
# stay 唯一队列全量重做共享库（v28.3 JIC 投稿版的非主分析）
#
# 统一口径：
#   * ID = stay_id；每 stay 取 order(stay_id,trial_num) 首 trial（stay 唯一）
#   * DAG14 基线；m=20 MICE；B=200 bootstrap；SEED=20261008
#   * 估计方法：
#       - 二分类结局：稳定化 IPTW logistic（cluster 无，单次横截面；MI+bootstrap）
#       - 连续结局：稳定化 IPTW 线性回归
#       - 时间-事件（亚组/交互/逐日/剂型/剂量）：稳定化 IPTW coxph
#   * 主分析及可日级 g-computation 项另用 v28.4_common.R + SEQTaRget
#
# 版本: V28.4   日期: 2026-10-10
# ============================================================================

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(survival)
})

V28_BASE <- "/data/小米云盘/文档/论文/血磷"
AUG_FILE <- file.path(V28_BASE, "TTE/v26/output/sequential_trials_augmented.parquet")
MICE_FMT <- file.path(V28_BASE,
  "sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet")
OC_FILE  <- file.path(V28_BASE, "sql/result/M0_07_outcome_extraction.parquet")
SEED     <- 20260629L  # 与 v28.3 投稿版 bootstrap seed 一致
M_IMP     <- 20L
B_BOOT    <- 200L

DAG_FROM_IMP <- c("age", "baseline_sofa", "sepsis3", "charlson", "egfr",
                  "invasive_vent_at_t0", "vasoactive_at_t0", "baseline_phosphate",
                  "lactate_t0", "refeeding_risk", "corrected_calcium_mmol_min")
DAG_FROM_AUG <- c("gender", "admission_type", "nutrition_support_type")
DAG14        <- c(DAG_FROM_IMP, DAG_FROM_AUG)
GENDER_LEV <- c("F", "M")
ADMIT_LEV  <- c("Observation", "Emergency", "Surgical_Elective")
NUTR_LEV   <- c("none", "enteral", "parenteral", "mixed")
BIN_COLS   <- c("sepsis3", "invasive_vent_at_t0", "vasoactive_at_t0",
                "refeeding_risk")
TRIM <- c(0.01, 0.99)

to01 <- function(x) if (is.numeric(x)) as.integer(x) else as.integer(as.character(x))

# ---- 1. stay 唯一人员表（含 DAG14、结局标志、时间）----
# imp_id: 1..20；cca=TRUE 时仅保留 DAG14 完整病例
load_stay <- function(imp_id) {
  aug <- as.data.table(read_parquet(AUG_FILE))
  imp <- as.data.table(read_parquet(sprintf(MICE_FMT, as.integer(imp_id))))
  oc  <- as.data.table(read_parquet(OC_FILE))
  setorder(aug, stay_id, trial_num)
  setorder(imp, stay_id, trial_num)

  a <- aug[, .SD[1], by = stay_id]
  i <- imp[, .SD[1], by = stay_id]
  stopifnot(identical(a$stay_id, i$stay_id))

  d <- a[, .(stay_id, subject_id, time_zero, itt_group, surv_time, death_28day_num,
             aki_kdigo2_plus, ever_crrt, transfusion_at_t0,
             cardiovascular_event_2comp, cardiovascular_event_3comp,
             new_arrhythmia_rhythm, myocardial_injury)]
  d[, arm := to01(i$treatment)]
  for (v in DAG_FROM_IMP) d[, (v) := i[[v]]]
  for (v in DAG_FROM_AUG) d[, (v) := a[[v]]]
  for (v in BIN_COLS)    d[, (v) := to01(d[[v]])]
  d[, gender := factor(gender, levels = GENDER_LEV)]
  d[, admission_type := factor(admission_type, levels = ADMIT_LEV)]
  d[, nutrition_support_type := factor(nutrition_support_type, levels = NUTR_LEV)]

  # 结局表（stay 唯一）
  setorder(oc, stay_id, time_zero)
  o <- oc[, .SD[1], by = stay_id]
  d <- merge(d,
    o[, .(stay_id, death_icu, death_hospital, death_90day,
          vfd_28, total_mv_days, phosphate_normalization_days, phosphate_normalized,
          icu_los_days, hospital_los_days,
          icu_los_from_t0_days, hospital_los_from_t0_days)],
    by = "stay_id", all.x = TRUE)
  d
}

# ---- 2. 设计矩阵（DAG14）----
dag_matrix <- function(d, covs = DAG14) {
  f <- as.formula(paste("~", paste(covs, collapse = "+")))
  X <- model.matrix(f, data = d)
  X[, !colnames(X) %in% "(Intercept)", drop = FALSE]
}

# ---- 3. 稳定化 IPTW（可换协变量集；返回权重向量）----
stab_weights <- function(d, covs = DAG14) {
  X <- dag_matrix(d, covs)
  ps <- glm.fit(X, d$arm, family = binomial())$fitted
  ps <- pmin(TRIM[2], pmax(TRIM[1], ps))
  p1 <- mean(d$arm)
  ifelse(d$arm == 1, p1 / ps, (1 - p1) / (1 - ps))
}

# 加权风险（用于 RD/RR）
weighted_risk <- function(d, ycol, w) {
  list(risk1 = weighted.mean(d[arm == 1][[ycol]], w[d$arm == 1]),
       risk0 = weighted.mean(d[arm == 0][[ycol]], w[d$arm == 0]))
}

# ---- 4. 二分类：加权 logistic，返回 logOR/SE ----
fit_bin <- function(d, ycol, w = NULL, covs = character(0)) {
  dd <- data.frame(arm = d$arm)
  if (length(covs)) dd <- cbind(dd, as.data.frame(dag_matrix(d, covs)))
  dd$.y <- d[[ycol]]
  if (is.null(w)) w <- rep(1, nrow(d))
  f <- as.formula(paste(".y~", paste(colnames(dd)[colnames(dd) != ".y"],
                                     collapse = "+")))
  fit <- glm(f, data = dd, weights = w, family = binomial())
  list(logor = unname(coef(fit)["arm"]),
       se = unname(sqrt(diag(vcov(fit)))["arm"]))
}

# ---- 5. 连续：加权 lm，返回 beta/SE ----
fit_lin <- function(d, ycol, w = NULL) {
  dd <- data.frame(arm = d$arm, .y = d[[ycol]])
  if (is.null(w)) w <- rep(1, nrow(d))
  fit <- lm(.y ~ arm, data = dd, weights = w)
  list(beta = unname(coef(fit)["arm"]),
       se = unname(sqrt(diag(vcov(fit)))["arm"]))
}

# ---- 6. 时间-事件：加权 coxph，返回系数表行 ----
fit_cox <- function(d, time, event, w = NULL, covs = DAG14,
                    extra = NULL, cluster = FALSE) {
  X <- dag_matrix(d, covs)
  if (!is.null(extra)) X <- cbind(X, extra)
  if (is.null(w)) w <- rep(1, nrow(d))
  df <- as.data.frame(X)
  df$arm <- d$arm
  df$.t <- time; df$.e <- event; df$.w <- w; df$.cl <- d$stay_id
  covn <- setdiff(colnames(X), character(0))
  rhs <- paste(c("arm", covn), collapse = "+")
  f <- as.formula(paste("Surv(.t,.e)~", rhs))
  cx <- coxph(f, data = df, weights = w,
              ties = "efron",
              cluster = if (cluster) .cl else NULL)
  s <- summary(cx)$coefficients
  rn <- if (is.matrix(s)) s["arm", ] else s
  list(loghr = unname(rn["coef"]), se = unname(rn["se(coef)"]),
       p = unname(rn[grep("Pr", names(rn))][1]))
}
