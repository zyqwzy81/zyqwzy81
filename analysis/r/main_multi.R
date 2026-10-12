# ============================================================================
# main_multi.R
# 多次入组（multi-enrollment）序贯目标试验：主结局全因死亡的 g-formula HR
#
# 功能开关：
#   --followup 28     随访天数（28 默认；90 跑 90 天死亡）
#   --firstonly       仅首次入组（selection.first_trial，单次入组敏感性）
#   --outcome oc_death 结局列（默认死亡；可切换其他二分类列）
#   --nboot 200       bootstrap 次数（0 仅点估计）
#   --bootimps 1      做 bootstrap 的 imp（默认全部；省时：点估计跨全部 imp，
#                     抽样方差取指定 imp 的 bootstrap + Rubin imp 间方差）
# 版本 V28.4  日期 2026-10-10
# ============================================================================
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/multi_common.R")

argv <- commandArgs(trailingOnly = TRUE)
gv <- function(k, d) { i <- match(k, argv); if (is.na(i)) d else argv[i + 1] }
IMPS    <- as.integer(strsplit(gv("--imps", "1"), ",")[[1]])
NBOOT   <- as.integer(gv("--nboot", "200"))
HORIZON <- as.integer(gv("--followup", "28"))
DEF_BOOT <- if (NBOOT > 0) paste(IMPS, collapse = ",") else ""
BOOTIMPS <- as.integer(strsplit(gv("--bootimps", DEF_BOOT), ",")[[1]])
OUTCOME <- gv("--outcome", "oc_death")
TAG     <- gv("--tag", "main")
FIRST   <- "--firstonly" %in% argv

OUT_T <- file.path(V28_BASE, "TTE/v28.4/output/tables/main_multi")
dir.create(OUT_T, recursive = TRUE, showWarnings = FALSE)

cat("build event offsets once...\n")
off <- build_event_offsets()

# 逐 imp：返回 logHR 以及（若 bootstrap）由 CI 反推的 logHR SE
run_imp <- function(imp_id) {
  cat(sprintf("\n================ imp %d ================\n", imp_id))
  trials <- load_trials(imp_id)
  daily  <- build_multi_daily(trials, off, H = HORIZON)
  do_boot <- NBOOT > 0 && imp_id %in% BOOTIMPS
  opts <- SEQopts(hazard = TRUE,
                  bootstrap = do_boot,
                  bootstrap.nboot = max(1L, NBOOT),
                  bootstrap.sample = 1,
                  followup.max = HORIZON,
                  seed = SEED + if (do_boot) imp_id else 0L,
                  bootstrap.CI_method = "se",
                  selection.first_trial = FIRST)
  mod <- SEQuential(daily, id.col = "stay_id", time.col = "day",
                    eligible.col = "eligible", treatment.col = "tx_init",
                    outcome.col = OUTCOME, time_varying.cols = DAG14,
                    fixed.cols = list(), method = "ITT",
                    options = opts, verbose = FALSE)
  h <- hazard_ratio(mod)            # c("Hazard ratio","LCI","UCI")
  loghr <- unname(log(h["Hazard ratio"]))
  se <- NA_real_
  if (do_boot && is.finite(h["LCI"]) && is.finite(h["UCI"]))
    se <- (as.numeric(log(h["UCI"])) - as.numeric(log(h["LCI"]))) /
          (2 * qnorm(0.975))
  list(loghr = loghr, se = se)
}

per_imp <- lapply(IMPS, run_imp)

# ---- Rubin 合并（点估计跨全部 imp；within 方差取 bootstrap imp 均值）----
g <- function(k) vapply(per_imp, function(r) r[[k]], numeric(1))
q <- g("loghr"); within <- g("se")
m <- length(q)
qbar <- mean(q)
bv <- if (m > 1) var(q) else 0
se_v <- within[is.finite(within)]
has_se <- length(se_v) > 0
if (has_se) {
  ubar <- mean(se_v^2)
  tv <- ubar + (1 + 1 / m) * bv
  dfm <- if (bv > 0) (m - 1) * (1 + ubar / ((1 + 1 / m) * bv))^2 else Inf
} else { tv <- NA_real_; dfm <- Inf }
zt <- if (is.finite(dfm)) qt(0.975, dfm) else qnorm(0.975)
pval <- if (has_se)
  2 * (if (is.finite(dfm)) pt(-abs(qbar / sqrt(tv)), dfm)
       else pnorm(-abs(qbar / sqrt(tv)))) else NA_real_

res <- data.table(outcome = OUTCOME, horizon = HORIZON,
                  HR = exp(qbar),
                  HR_lci = if (has_se) exp(qbar - zt * sqrt(tv)) else NA_real_,
                  HR_uci = if (has_se) exp(qbar + zt * sqrt(tv)) else NA_real_,
                  p = pval, n_imp = m, first_only = FIRST, n_boot = NBOOT)
cat("\n========== MAIN multi-enrollment ==========\n")
print(res)

suffix <- paste0(if (FIRST) "firstonly_" else "", "h", HORIZON, "_", TAG)
fwrite(res, file.path(OUT_T, paste0("main_", suffix, ".csv")))
cat("saved:", file.path(OUT_T, paste0("main_", suffix, ".csv")), "\n")
