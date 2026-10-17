# ============================================================================
# v28.4_pp_grace.R
# Grace-Window 阶梯敏感性: g ∈ {0h, 6h, 12h, 24h, 48h}
# 每个 g 值 = 独立跑一遍 SEQTaRget censoring + excused.cols
# m=20, boot=200, 全跑
#
# 机制:
#   grace_h 小时 = grace_h/24 天
#   arm=0 (对照): grace 窗口内交叉 IV → 不删失 (gr0=TRUE)
#   arm=1 (治疗): grace 窗口内延迟 → 不删失 (gr1=TRUE)
#   SEQopts(excused=TRUE, excused.cols=c("gr0","gr1")) 传入
# ============================================================================
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/v28.4_multi_common.R")
suppressPackageStartupMessages({library(data.table); library(SEQTaRget); library(parallel)})
setDTthreads(2L); Sys.setenv(OMP_NUM_THREADS = 2L)
SEED <- 20260629L; TAU <- 28L
AUG_F <- file.path(V28_BASE, "TTE/v26/output/sequential_trials_augmented.parquet")

# ---- 实际 IV 日 (同 pp_multi.R) ----
build_iv_day <- function(trials) {
  aug <- as.data.table(read_parquet(AUG_F))
  aug <- aug[, .(stay_id, time_zero, itt_group,
                 fivh = first_iv_hours_from_tz)]
  omap <- unique(trials[, .(stay_id, time_zero, origin_day, arm)])
  z <- merge(aug, omap, by = c("stay_id", "time_zero"), all.x = TRUE)
  z[, ivd := fifelse(itt_group == "treatment" | arm == 1L, origin_day,
                     fifelse(is.finite(fivh),
                             origin_day + as.integer(ceiling(fivh / 24)),
                             NA_integer_))]
  z[, .(iv_abs_day = min(ivd, na.rm = TRUE)), by = stay_id]
}

IMPS   <- 1
NBOOT  <- 200L
NCORES <- 2L

GRACE_HOURS <- c(0, 24, 48)

OUT_T  <- file.path(V28_BASE, "TTE/v28.4/output/tables/sensitivity")
OUT_LOG <- file.path(V28_BASE, "TTE/v28.4/output/logs")
dir.create(OUT_T, recursive = TRUE, showWarnings = FALSE)

cat("Grace-Window PP:", length(GRACE_HOURS), "grace values × m=20 × boot=200\n")
cat("开始:", format(Sys.time(), "%H:%M:%S"), "\n")

# ---- 事件偏移 (一次) ----
off <- build_event_offsets()

# ---- 构造 grace 列 ----
# daily 已经有 iv_abs_day (绝对 IV 日) 和 tx_act (实际 IV 轨迹)
# stay_arm = arm 1=治疗 0=对照 (stay 级)
#
# grace 语义:
#   arm=0 者 cross-IV 在 grace 窗口内 → 全程不删失 (gr0=TRUE)
#   即: grace 够长 → 大部分 cross 都被宽恕 → HR 向 ITT 靠拢
#   grace=0h → strict censoring (cross 日即删)
#   grace=48h → "晚一点也没关系" → 接近 CCW 不删失
#
#   arm=1 起点即 IV, censoring 本来就不删 arm=1 → gr1 无影响
add_grace_cols <- function(daily, grace_days) {
  # stay 级 arm
  daily[, stay_arm := fifelse(any(tx_init == 1L), 1L, 0L), by = stay_id]
  daily[, `:=`(gr0 = FALSE, gr1 = FALSE)]

  # arm=0: cross 日在 grace_days 内 → 从起点到交叉日之后 grace 期都不删失
  # 更简单: 只要 arm=0 stay 有 cross 且 cross 日 <= grace_days → 全程 gr0=TRUE
  # (实际 SEQTaRget censoring 是逐日比较, 我们让 cross 后 grace 内的删失都被 excused)
  cross_in_grace <- daily[stay_arm == 0L,
                          .(cross_in_grace =
                              any(is.finite(iv_abs_day) & iv_abs_day <= grace_days)),
                          by = stay_id]
  daily[cross_in_grace, on = "stay_id",
        gr0 := i.cross_in_grace & stay_arm == 0L]

  # arm=1: censoring 本来就不删 arm=1 → gr1=TRUE 没影响
  daily[stay_arm == 1L, gr1 := TRUE]

  daily[, stay_arm := NULL]
  daily
}

# ---- 单 imp 单 grace 执行 ----
run_one <- function(imp_id, grace_days) {
  trials <- load_trials(imp_id)
  daily  <- build_multi_daily(trials, off, H = TAU)
  ivmap  <- build_iv_day(trials)
  daily[ivmap, on = "stay_id", iv_abs_day := i.iv_abs_day]
  daily[, tx_act := as.integer(is.finite(iv_abs_day) & day >= iv_abs_day)]

  # 加 grace excused 列 (需要 iv_abs_day)
  daily <- add_grace_cols(daily, grace_days)
  daily[, iv_abs_day := NULL]  # 不要传到 SEQuential

  opts <- SEQopts(hazard = TRUE,
                  weighted = TRUE, weight.preexpansion = TRUE, weight.p99 = TRUE,
                  bootstrap = TRUE, bootstrap.nboot = NBOOT,
                  bootstrap.sample = 1,
                  followup.max = TAU,
                  seed = SEED + imp_id * 100L + as.integer(grace_days * 24),
                  bootstrap.CI_method = "se",
                  excused = TRUE, excused.cols = c("gr0", "gr1"),
                  selection.first_trial = FALSE)
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
  list(loghr = loghr, se = se)
}

# ---- 主循环 ----
res_l <- list()
for (gi in seq_along(GRACE_HOURS)) {
  gh <- GRACE_HOURS[gi]
  gd <- gh / 24
  cat(sprintf("\n===== grace = %dh (%.2fd) =====\n", gh, gd))
  t0 <- Sys.time()

  per_imp <- vector("list", length(IMPS))
  for (k in seq_along(IMPS)) {
    im <- IMPS[k]
    t1 <- Sys.time()
    r <- run_one(im, gd)
    per_imp[[k]] <- data.table(imp = im, loghr = r$loghr, se = r$se)
    elapsed <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    cat(sprintf("  imp %2d/%2d · HR=%.4f · %.0fs\n", im, length(IMPS), exp(r$loghr), elapsed))
  }
  per_dt <- rbindlist(per_imp)

  # Rubin
  q <- per_dt$loghr; within <- per_dt$se; m <- length(q)
  qbar <- mean(q); bv <- if (m > 1) var(q) else 0
  se_v <- within[is.finite(within)]
  ubar <- mean(se_v^2); tv <- ubar + (1 + 1 / m) * bv
  dfm <- if (bv > 0) (m - 1) * (1 + ubar / ((1 + 1 / m) * bv))^2 else Inf
  zt <- if (is.finite(dfm)) qt(.975, dfm) else qnorm(.975)
  pval <- 2 * (if (is.finite(dfm)) pt(-abs(qbar / sqrt(tv)), dfm)
               else pnorm(-abs(qbar / sqrt(tv))))

  row <- data.table(grace_h = gh, grace_d = gd,
                    HR = exp(qbar),
                    HR_lci = exp(qbar - zt * sqrt(tv)),
                    HR_uci = exp(qbar + zt * sqrt(tv)),
                    p = pval, n_imp = m, n_boot = NBOOT,
                    hr_sd_per_imp = sd(exp(q)))
  res_l[[gi]] <- row
  elapsed_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  → HR=%.4f (%.4f–%.4f), p=%.3f, %.1f min total\n",
              row$HR, row$HR_lci, row$HR_uci, row$p, elapsed_total/60))
}

# ---- 输出 ----
final <- rbindlist(res_l)
cat("\n========== Grace-Window PP 最终结果 ==========\n")
print(final)

fp <- file.path(OUT_T, "v28.4_pp_grace_sensitivity.csv")
fwrite(final, fp)
cat("\nsaved:", fp, "\n")
cat("结束:", format(Sys.time(), "%H:%M:%S"), "\n")
