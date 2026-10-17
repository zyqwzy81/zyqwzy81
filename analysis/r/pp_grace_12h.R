# ============================================================================
# v28.4_pp_grace_12h.R
# Grace-Window PP: g = 12h (0.5 天)
# m = 20, boot = 200 (每个 imp 全 boot)
#
# 机制 (SEQTaRget method="censoring" + excused.cols):
#   grace_h = 12h = 0.5d
#   arm=0 者 cross IV 在 grace 窗口内 → 不删失 (gr0=TRUE)
#   arm=1 者无删失 → gr1=TRUE
#   SEQopts(excused=TRUE, excused.cols=c("gr0","gr1"))
# ============================================================================
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/v28.4_multi_common.R")
suppressPackageStartupMessages({library(data.table); library(SEQTaRget)})
setDTthreads(2L); Sys.setenv(OMP_NUM_THREADS = 2L)
SEED <- 20260629L; TAU <- 28L; GRACE_H <- 12L
AUG_F <- file.path(V28_BASE, "TTE/v26/output/sequential_trials_augmented.parquet")

IMPS  <- 1:20
NBOOT <- 200L

OUT_T  <- file.path(V28_BASE, "TTE/v28.4/output/tables/sensitivity")
dir.create(OUT_T, recursive = TRUE, showWarnings = FALSE)

cat("════════════ Grace-Window PP", GRACE_H, "h ════════════\n")
cat("m=", length(IMPS), ", boot=", NBOOT, "\n", sep="")
cat("开始:", format(Sys.time(), "%H:%M:%S"), "\n")

# ---- IV 日 ----
build_iv_day <- function(trials) {
  aug <- as.data.table(read_parquet(AUG_F))
  aug <- aug[, .(stay_id, time_zero, itt_group, fivh = first_iv_hours_from_tz)]
  omap <- unique(trials[, .(stay_id, time_zero, origin_day, arm)])
  z <- merge(aug, omap, by = c("stay_id", "time_zero"), all.x = TRUE)
  z[, ivd := fifelse(itt_group == "treatment" | arm == 1L, origin_day,
                     fifelse(is.finite(fivh), origin_day + as.integer(ceiling(fivh / 24)),
                             NA_integer_))]
  z[, .(iv_abs_day = min(ivd, na.rm = TRUE)), by = stay_id]
}

# ---- 事件偏移 ----
off <- build_event_offsets()
grace_days <- GRACE_H / 24

# ---- 单 imp ----
run_one <- function(imp_id) {
  trials <- load_trials(imp_id)
  daily  <- build_multi_daily(trials, off, H = TAU)
  ivmap  <- build_iv_day(trials)
  daily[ivmap, on = "stay_id", iv_abs_day := i.iv_abs_day]
  daily[, tx_act := as.integer(is.finite(iv_abs_day) & day >= iv_abs_day)]

  # stay 级 arm
  daily[, stay_arm := fifelse(any(tx_init == 1L), 1L, 0L), by = stay_id]

  # grace excused 列: arm=0 cross IV 在 grace_days 内 → 全程 gr0=TRUE
  cross_dt <- daily[stay_arm == 0L,
                    .(cross_in_grace = any(is.finite(iv_abs_day) & iv_abs_day <= grace_days)),
                    by = stay_id]
  daily[, `:=`(gr0 = FALSE, gr1 = FALSE)]
  daily[cross_dt, on = "stay_id", gr0 := i.cross_in_grace & stay_arm == 0L]
  daily[stay_arm == 1L, gr1 := TRUE]

  # 清理临时列
  daily[, c("stay_arm", "iv_abs_day") := NULL]

  opts <- SEQopts(hazard = TRUE,
                  weighted = TRUE, weight.preexpansion = TRUE, weight.p99 = TRUE,
                  bootstrap = TRUE, bootstrap.nboot = NBOOT,
                  bootstrap.sample = 1, followup.max = TAU,
                  seed = SEED + imp_id,
                  bootstrap.CI_method = "se",
                  excused = TRUE, excused.cols = c("gr0", "gr1"))
  mod <- SEQuential(daily, id.col = "stay_id", time.col = "day",
                    eligible.col = "eligible", treatment.col = "tx_act",
                    outcome.col = "oc_death", time_varying.cols = DAG14,
                    fixed.cols = list(), method = "censoring",
                    options = opts, verbose = FALSE)
  h <- hazard_ratio(mod)
  loghr <- unname(log(h["Hazard ratio"]))
  se <- NA_real_
  if (is.finite(h["LCI"]) && is.finite(h["UCI"]))
    se <- (as.numeric(log(h["UCI"])) - as.numeric(log(h["LCI"]))) /
          (2 * qnorm(0.975))
  rm(mod); invisible(gc())
  data.table(imp = imp_id, loghr = loghr, se = se,
             hr = exp(loghr), lci = h["LCI"], uci = h["UCI"])
}

# ---- 主循环 ----
t0_all <- Sys.time()
per_imp <- vector("list", length(IMPS))
for (k in seq_along(IMPS)) {
  im <- IMPS[k]
  t1 <- Sys.time()
  r <- run_one(im)
  per_imp[[k]] <- r
  elapsed <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  cat(sprintf("  imp %2d/%2d · HR=%.4f · %.0fs\n",
              im, length(IMPS), r$hr, elapsed))
}
per_dt <- rbindlist(per_imp)

# ---- Rubin 合并 ----
q <- per_dt$loghr; within <- per_dt$se; m <- length(q)
qbar <- mean(q); bv <- if (m > 1) var(q) else 0
ubar <- mean(within[is.finite(within)]^2)
tv <- ubar + (1 + 1 / m) * bv
dfm <- if (bv > 0) (m - 1) * (1 + ubar / ((1 + 1 / m) * bv))^2 else Inf
zt <- if (is.finite(dfm)) qt(.975, dfm) else qnorm(.975)
pval <- 2 * (if (is.finite(dfm)) pt(-abs(qbar / sqrt(tv)), dfm)
             else pnorm(-abs(qbar / sqrt(tv))))

final <- data.table(analysis = "Grace-Window PP (censoring + 12h grace)",
                    grace_h = GRACE_H, grace_d = grace_days,
                    m = m, n_boot = NBOOT,
                    HR = exp(qbar),
                    HR_lci = exp(qbar - zt * sqrt(tv)),
                    HR_uci = exp(qbar + zt * sqrt(tv)),
                    p = pval,
                    bv = bv, ubar = ubar, tv = tv)

elapsed_total <- as.numeric(difftime(Sys.time(), t0_all, units = "secs"))
cat("\n════════════ Grace-Window PP 12h 结果 ════════════\n")
print(final)
cat(sprintf("总耗时 %.1f min\n", elapsed_total / 60))

fwrite(final, file.path(OUT_T, "v28.4_pp_grace_12h_pooled.csv"))
fwrite(per_dt, file.path(OUT_T, "v28.4_pp_grace_12h_per_imp.csv"))
cat("\nsaved.\n")
cat("结束:", format(Sys.time(), "%H:%M:%S"), "\n")
