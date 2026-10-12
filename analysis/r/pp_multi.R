# ============================================================================
# pp_multi.R
# 多次入组（multi-enrollment）框架下：Per-Protocol（PP）g-formula HR
#
# 机制（包 SEQexpand method="censoring"）：
#   * treatment.col 提供【逐时点实际 IV 状态】轨迹 tx_act（0/1）
#   * 包在克隆窗口内逐日比较实际状态 vs 起点策略：
#       - 治疗策略：起点即 IV（state=1），此后保持 1 → 不删失
#       - 对照策略：保持 0；对照者交叉到 IV 当日 state 0→1 → 首次改变即删失
#   * 包内部对删失建 numerator/denominator IPCW，p99 截断
#
# 开关：--imps/--nboot/--bootimps/--firstonly/--tag
# 版本 V28.4  日期 2026-10-10
# ============================================================================
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/multi_common.R")

argv <- commandArgs(trailingOnly = TRUE)
gv <- function(k, d) { i <- match(k, argv); if (is.na(i)) d else argv[i + 1] }
IMPS   <- as.integer(strsplit(gv("--imps", "1"), ",")[[1]])
NBOOT  <- as.integer(gv("--nboot", "200"))
DEF_BOOT <- if (NBOOT > 0) paste(IMPS, collapse = ",") else ""
BOOTIMPS <- as.integer(strsplit(gv("--bootimps", DEF_BOOT), ",")[[1]])
TAG    <- gv("--tag", "pp")
FIRST  <- "--firstonly" %in% argv

OUT_T <- file.path(V28_BASE, "TTE/v28.4/output/tables/pp_multi")
dir.create(OUT_T, recursive = TRUE, showWarnings = FALSE)

AUG_F <- file.path(V28_BASE, "TTE/v26/output/sequential_trials_augmented.parquet")

cat("build event offsets once...\n")
off <- build_event_offsets()

# 逐 stay 的最早绝对 IV 日（治疗起点日 或 对照交叉日）
build_iv_day <- function(trials) {
  aug <- as.data.table(read_parquet(AUG_F))
  aug <- aug[, .(stay_id, time_zero, itt_group,
                 fivh = first_iv_hours_from_tz)]
  # origin_day 映射
  omap <- unique(trials[, .(stay_id, time_zero, origin_day, arm)])
  z <- merge(aug, omap, by = c("stay_id", "time_zero"), all.x = TRUE)
  # 治疗：IV 日 = 起点日；对照交叉：起点日 + ceil(fivh/24)
  z[, ivd := fifelse(itt_group == "treatment" | arm == 1L, origin_day,
                     fifelse(is.finite(fivh),
                             origin_day + as.integer(ceiling(fivh / 24)),
                             NA_integer_))]
  z[, .(iv_abs_day = min(ivd, na.rm = TRUE)), by = stay_id]
}

run_imp <- function(imp_id) {
  cat(sprintf("\n================ imp %d ================\n", imp_id))
  trials <- load_trials(imp_id)
  daily  <- build_multi_daily(trials, off)
  # 实际 IV 状态轨迹：day>=iv_abs_day → 1
  ivmap <- build_iv_day(trials)
  daily[ivmap, on = "stay_id", iv_abs_day := i.iv_abs_day]
  daily[, tx_act := as.integer(is.finite(iv_abs_day) & day >= iv_abs_day)]
  daily[, iv_abs_day := NULL]

  do_boot <- NBOOT > 0 && imp_id %in% BOOTIMPS
  opts <- SEQopts(hazard = TRUE,
                  weighted = TRUE,
                  weight.preexpansion = TRUE,
                  weight.p99 = TRUE,
                  bootstrap = do_boot,
                  bootstrap.nboot = max(1L, NBOOT),
                  bootstrap.sample = 1,
                  followup.max = TAU,
                  seed = SEED + if (do_boot) imp_id else 0L,
                  bootstrap.CI_method = "se",
                  selection.first_trial = FIRST)
  mod <- SEQuential(daily, id.col = "stay_id", time.col = "day",
                    eligible.col = "eligible", treatment.col = "tx_act",
                    outcome.col = "oc_death", time_varying.cols = DAG14,
                    fixed.cols = list(), method = "censoring",
                    options = opts, verbose = TRUE)
  h <- hazard_ratio(mod)
  loghr <- unname(log(h["Hazard ratio"]))
  se <- NA_real_
  if (do_boot && is.finite(h["LCI"]) && is.finite(h["UCI"]))
    se <- (as.numeric(log(h["UCI"])) - as.numeric(log(h["LCI"]))) /
          (2 * qnorm(0.975))
  list(loghr = loghr, se = se)
}

per_imp <- lapply(IMPS, run_imp)

g <- function(k) vapply(per_imp, function(r) r[[k]], numeric(1))
q <- g("loghr"); within <- g("se"); m <- length(q)
qbar <- mean(q); bv <- if (m > 1) var(q) else 0
se_v <- within[is.finite(within)]; has_se <- length(se_v) > 0
if (has_se) {
  ubar <- mean(se_v^2); tv <- ubar + (1 + 1 / m) * bv
  dfm <- if (bv > 0) (m - 1) * (1 + ubar / ((1 + 1 / m) * bv))^2 else Inf
} else { tv <- NA_real_; dfm <- Inf }
zt <- if (is.finite(dfm)) qt(.975, dfm) else qnorm(.975)
pval <- if (has_se)
  2 * (if (is.finite(dfm)) pt(-abs(qbar / sqrt(tv)), dfm)
       else pnorm(-abs(qbar / sqrt(tv)))) else NA_real_

res <- data.table(analysis = "Per-protocol (censoring IPCW, multi-enrollment)",
                  HR = exp(qbar),
                  HR_lci = if (has_se) exp(qbar - zt * sqrt(tv)) else NA_real_,
                  HR_uci = if (has_se) exp(qbar + zt * sqrt(tv)) else NA_real_,
                  p = pval, n_imp = m, first_only = FIRST, n_boot = NBOOT)
cat("\n========== PP multi-enrollment ==========\n"); print(res)
suffix <- paste0(if (FIRST) "firstonly_" else "", TAG)
fwrite(res, file.path(OUT_T, paste0("pp_", suffix, ".csv")))
cat("saved:", file.path(OUT_T, paste0("pp_", suffix, ".csv")), "\n")
