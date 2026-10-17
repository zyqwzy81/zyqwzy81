#!/usr/bin/env Rscript
# ============================================================================
# v28.4_cca_secondary.R
# SA7 扩展：完整病例（DAG14 全部观测、未插补、stay-unique）子集上的
#           全部 28 天二分类重要终点
#
# 背景（应审稿人 R3-M8）：
#   "Please accompany the complete-case analysis results for the primary
#    endpoint and important endpoints in the main body."
#   主终点 CCA 已由 v28.4_sa3_sa7.R 给出 HR=0.962 (0.857-1.078),
#   n=14,293 stays / 2,601 deaths。本脚本补齐同口径重要终点
#   （AKI>=2、PRBC 输血、心律失常、心肌损伤、ICU/院内死亡、CV 复合），
#   方法与 TP 级 Table 3 完全一致：SEQTaRget ITT g-formula（km.curves，
#   风险标准化），区别仅在数据=未插补完整病例、单位=stay-unique。
#
# 推断：无 MI → 单一数据集；B=200 stay_id cluster bootstrap（bootstrap.
#       sample=1，与主 CCA 一致），"se" 法给 log(RR)/RD 的 SE，
#       正态 z CI 与 p（单数据集无 Rubin，B=0 退化口径，同 sa7）。
#
# 用法：
#   验证: LD_LIBRARY_PATH=/opt/miniforge3/lib Rscript v28.4_cca_secondary.R --nboot 0
#   全量: LD_LIBRARY_PATH=/opt/miniforge3/lib Rscript v28.4_cca_secondary.R --nboot 200
# 版本 V28.4  日期 2026-10-14  seed=20260629
# ============================================================================
suppressPackageStartupMessages(library(SEQTaRget))
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/v28.4_multi_common.R")

argv <- commandArgs(trailingOnly = TRUE)
gv <- function(k, d) { i <- match(k, argv); if (is.na(i)) d else argv[i + 1] }
NBOOT <- as.integer(gv("--nboot", "200"))
TAG   <- gv("--tag", "cca")
OCS <- if ("--outcomes" %in% argv)
  strsplit(gv("--outcomes", paste(BIN_OUTCOMES, collapse = ",")), ",")[[1]] else BIN_OUTCOMES

OUT_T <- file.path(V28_BASE, "TTE/v28.4/output/tables/nonprimary_redo")
dir.create(OUT_T, recursive = TRUE, showWarnings = FALSE)

# ---- 未插补 stay-unique 试验行（与 v28.4_sa3_sa7.R 同一 loader 口径）----
load_raw_trials <- function() {
  aug <- as.data.table(read_parquet(AUG_FILE))
  setorder(aug, stay_id, trial_num)
  dt <- aug[, .SD[1], by = stay_id]
  out <- dt[, .(stay_id, trial_num, time_zero, itt_group,
                death_28day_num, surv_time)]
  out[, origin_day := as.integer(as.Date(time_zero) -
                                   min(as.Date(time_zero))), by = stay_id]
  out[, arm := fifelse(itt_group == "treatment", 1L, 0L)]
  for (v in DAG_FROM_IMP) out[, (v) := dt[[v]]]
  for (v in DAG_FROM_AUG) out[, (v) := dt[[v]]]
  for (v in BIN_COLS)    out[, (v) := to01(out[[v]])]
  out[, gender := factor(gender, levels = GENDER_LEV)]
  out[, admission_type := factor(admission_type, levels = ADMIT_LEV)]
  out[, nutrition_support_type := factor(nutrition_support_type,
                                         levels = NUTR_LEV)]
  out
}

# ---- τ=28 处风险/RR/RD/SE 抽取（同 v28.4_secondary_multi.R）----
# bootstrap=TRUE 时 SEQTaRget 在同一 model 内同时返回点估计与 SE/LCI/UCI，
# 因此只需跑一次（B=0 时 bootstrap=FALSE，SE 为 NA）。
extract_at_tau <- function(mod) {
  rd_all <- as.data.table(mod@risk.data[[1]])
  rc_all <- as.data.table(mod@risk.comparison[[1]])
  fu <- TAU
  if (!fu %in% rd_all$Followup) fu <- max(rd_all$Followup)
  rd <- rd_all[Followup == fu]
  r0 <- rd[A == "0"]; r1 <- rd[A == "1"]
  rc <- rc_all[Followup == fu & A_x == "risk_0" & A_y == "risk_1"]
  if (!nrow(rc))  # 健壮兜底：取首个 risk_0→risk_1 行
    rc <- rc_all[A_x == "risk_0" & A_y == "risk_1"][1]
  safe <- function(df, nm) {
    if (nrow(df) && nm %in% names(df)) as.numeric(df[[nm]])[1] else NA_real_
  }
  list(risk0 = safe(r0, "Risk"), risk1 = safe(r1, "Risk"),
       risk0_se = safe(r0, "SE"), risk1_se = safe(r1, "SE"),
       rr = safe(rc, "Risk Ratio"),
       logrr_se = safe(rc, "log(RR) SE"),
       rd = safe(rc, "Risk Difference"),
       rd_se = safe(rc, "RD SE"),
       fu = fu)
}

# ---- 构建 CCA 子集（硬断言锁定与主 CCA 同一人群）----
raw <- load_raw_trials()
comp <- raw[complete.cases(raw[, ..DAG14])]
n_cc <- nrow(comp); n_dth <- sum(comp$death_28day_num)
cat(sprintf("\nCCA n=%d deaths=%d (expect 14293/2601)\n", n_cc, n_dth))
stopifnot(n_cc == 14293L, n_dth == 2601L)
cat("arm table:\n"); print(table(comp$arm))

cat("build event offsets (once)...\n")
off <- build_event_offsets()
daily <- build_multi_daily(comp, off)
check_multi_daily(daily, "cca")

# 每结局 CCA 事件数（stay-unique：daily 列和即 stay 级事件数）
ev_counts <- daily[, lapply(.SD, sum), .SDcols = BIN_OUTCOMES]
print(ev_counts)
stopifnot(ev_counts$oc_death == n_dth)  # 死亡事件须与主 CCA 一致

run_one <- function(oc) {
  do_boot <- NBOOT > 0
  cat(sprintf("\n-- outcome %s (bootstrap=%s, B=%d) --\n", oc, do_boot, NBOOT))
  opts <- SEQopts(km.curves = TRUE,
                  bootstrap = do_boot,
                  bootstrap.nboot = max(1L, NBOOT),
                  bootstrap.sample = 1,
                  followup.max = TAU,
                  seed = SEED,
                  bootstrap.CI_method = "se",
                  risk.times = TAU)
  mod <- SEQuential(daily, id.col = "stay_id", time.col = "day",
                    eligible.col = "eligible", treatment.col = "tx_init",
                    outcome.col = oc, time_varying.cols = DAG14,
                    fixed.cols = list(), method = "ITT",
                    options = opts, verbose = FALSE)
  extract_at_tau(mod)
}

res <- list()
for (oc in OCS) {
  pt <- run_one(oc)   # 单 model：B>0 时点估计与 bootstrap SE 同源
  z <- qnorm(0.975)
  logrr <- log(pt$rr)
  se_rr <- pt$logrr_se
  se_rd <- pt$rd_se
  p_rr  <- if (is.finite(se_rr)) 2 * pnorm(-abs(logrr / se_rr)) else NA_real_
  cat(sprintf("   RR=%.3f logRR_SE=%s RD=%.4f RD_SE=%s\n",
              pt$rr, ifelse(is.finite(se_rr), round(se_rr, 4), "NA"),
              pt$rd, ifelse(is.finite(se_rd), round(se_rd, 4), "NA")))
  res[[oc]] <- data.table(
    outcome = oc, n = n_cc, n_event = as.integer(ev_counts[[oc]]),
    risk_control = pt$risk0, risk_treat = pt$risk1,
    RR = pt$rr,
    RR_lci = if (is.finite(se_rr)) exp(logrr - z * se_rr) else NA_real_,
    RR_uci = if (is.finite(se_rr)) exp(logrr + z * se_rr) else NA_real_,
    RD = pt$rd,
    RD_lci = if (is.finite(se_rd)) pt$rd - z * se_rd else NA_real_,
    RD_uci = if (is.finite(se_rd)) pt$rd + z * se_rd else NA_real_,
    p_RR = p_rr, n_boot = NBOOT)
}
final <- rbindlist(res)
cat("\n========== CCA (DAG14 complete-case, stay-unique, n=", n_cc,
    ") B=", NBOOT, " ==========\n", sep = "")
print(final)

fwrite(final, file.path(OUT_T, paste0("v28.4_secondary_", TAG,
                                     "_b", NBOOT, ".csv")))
cat("\nsaved:", file.path(OUT_T, paste0("v28.4_secondary_", TAG,
                                        "_b", NBOOT, ".csv")), "\n")
