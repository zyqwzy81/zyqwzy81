# ============================================================================
# v28.4_main_multi_dagset.R  — 重跑为 DAG14 m=20
# 多次入组 TP 主结局 g-formula：DAG14 协变量集，m=20 Rubin + B=200 *每个 imp 全 boot (过夜全跑)
# 协议（与所有其他分析一致）：点估计跨 imp1..20，within 方差取 imp1 的 B=200
# stay 聚类 bootstrap，Rubin/Barnard-Rubin t 合并
# 版本 V28.4-revise  日期 2026-10-16
# ============================================================================
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/v28.4_multi_common.R")
suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)

argv <- commandArgs(trailingOnly = TRUE)
gv <- function(k, d) { i <- match(k, argv); if (is.na(i)) d else argv[i + 1] }
IMPS    <- 1:20
NBOOT   <- as.integer(gv("--nboot", "200"))
HORIZON <- 28L
TAG     <- "mimic_dag14_m20b200_fullboot"

# DAG14 (正式手稿协变量集，与 v28.4_common.R / 所有其他分析一致)
DAG14 <- c("age", "gender", "baseline_sofa", "sepsis3", "charlson",
           "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
           "admission_type", "baseline_phosphate", "lactate_t0", "refeeding_risk",
           "corrected_calcium_mmol_min", "nutrition_support_type")
COVS <- DAG14

OUT_T <- file.path(V28_BASE, "TTE/v28.4/output/tables/sensitivity")
dir.create(OUT_T, recursive = TRUE, showWarnings = FALSE)

cat("covset: DAG14 (", length(COVS), " vars)\n", sep = "")
cat("build event offsets once...\n")
off <- build_event_offsets()

run_imp <- function(imp_id) {
  cat(sprintf("\n================ imp %d ================\n", imp_id))
  trials <- load_trials(imp_id)
  daily  <- build_multi_daily(trials, off, H = HORIZON)
  do_boot <- TRUE                                   # 每个 imp 都 boot=200 (过夜全跑)
  opts <- SEQopts(hazard = TRUE,
                  bootstrap = do_boot,
                  bootstrap.nboot = max(1L, NBOOT),
                  bootstrap.sample = 1,
                  followup.max = HORIZON,
                  seed = SEED + if (do_boot) imp_id else 0L,
                  bootstrap.CI_method = "se")
  mod <- SEQuential(daily, id.col = "stay_id", time.col = "day",
                    eligible.col = "eligible", treatment.col = "tx_init",
                    outcome.col = "oc_death", time_varying.cols = COVS,
                    fixed.cols = list(), method = "ITT",
                    options = opts, verbose = FALSE)
  h <- hazard_ratio(mod)
  loghr <- unname(log(h["Hazard ratio"]))
  se <- NA_real_
  if (do_boot && is.finite(h["LCI"]) && is.finite(h["UCI"]))
    se <- (as.numeric(log(h["UCI"])) - as.numeric(log(h["LCI"]))) /
          (2 * qnorm(0.975))
  list(loghr = loghr, se = se)
}

per_imp <- lapply(IMPS, run_imp)

q <- vapply(per_imp, function(r) r$loghr, numeric(1))
within <- vapply(per_imp, function(r) r$se, numeric(1))
m <- length(q); qbar <- mean(q); bv <- if (m > 1) var(q) else 0
se_v <- within[is.finite(within)]; ubar <- mean(se_v^2)
tv <- ubar + (1 + 1 / m) * bv
dfm <- if (bv > 0) (m - 1) * (1 + ubar / ((1 + 1 / m) * bv))^2 else Inf
zt <- if (is.finite(dfm)) qt(0.975, dfm) else qnorm(0.975)
pval <- 2 * (if (is.finite(dfm)) pt(-abs(qbar / sqrt(tv)), dfm)
             else pnorm(-abs(qbar / sqrt(tv))))

res <- data.table(cov_set = "DAG14", n_cov = length(COVS),
                  HR = exp(qbar),
                  HR_lci = exp(qbar - zt * sqrt(tv)),
                  HR_uci = exp(qbar + zt * sqrt(tv)),
                  p = pval, n_imp = m, n_boot = NBOOT,
                  perimp_hr_min = exp(min(q)), perimp_hr_max = exp(max(q)))
cat("\n========== RESULT DAG14 m20/B200 FULLBOOT ==========\n"); print(res)
fp <- file.path(OUT_T, paste0("v28.4_", TAG, ".csv"))
fwrite(res, fp); cat("saved:", fp, "\n")
