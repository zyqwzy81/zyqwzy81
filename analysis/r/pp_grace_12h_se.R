# ============================================================================
# v28.4_pp_grace_12h_se.R
# 目标: 跑 imp1 boot=200 拿 Grace-Window PP 12h 的精确 within-se
#       （旧跑已得 20 imps loghr，但缺 boot=200 的 se）
#
# 机制: 手动 excusal (绕开 SEQTaRget excused=TRUE bug)
#   grace_days = 1 (day 0-1 ≈ 12h，因 iv_abs_day 是 integer)
#   arm=0 cross stay 且 cross day <= 1 → 强制 tx_act := tx_init
#   → SEQTaRget method="censoring" + weighted=TRUE + excused=FALSE
# ============================================================================
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/v28.4_multi_common.R")
suppressPackageStartupMessages({library(data.table); library(SEQTaRget)})
setDTthreads(2L); Sys.setenv(OMP_NUM_THREADS = 2L)
SEED <- 20260629L; TAU <- 28L; GRACE_DAYS <- 1L
AUG_F <- file.path(V28_BASE, "TTE/v26/output/sequential_trials_augmented.parquet")
IM_ID <- 1L
NBOOT <- 200L

OUT_T  <- file.path(V28_BASE, "TTE/v28.4/output/tables/sensitivity")
OUT_LOG <- file.path(V28_BASE, "TTE/v28.4/output/logs")
dir.create(OUT_T, recursive = TRUE, showWarnings = FALSE)

cat("════════ Grace-Window PP 12h se-probe ════════\n")
cat("手动 excusal, grace_days =", GRACE_DAYS, "(day 0-1 ≈ 12h)\n")
cat("imp =", IM_ID, ", boot =", NBOOT, "\n")
cat("开始:", format(Sys.time(), "%H:%M:%S"), "\n")

# ---- IV 日 ----
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

# ---- 事件偏移 ----
cat("[offsets] building event offsets...\n")
off <- build_event_offsets()

# ---- imp1 数据 ----
cat("[data] loading imp", IM_ID, "...\n")
t0 <- Sys.time()
trials <- load_trials(IM_ID)
daily  <- build_multi_daily(trials, off, H = TAU)
ivmap  <- build_iv_day(trials)
daily[ivmap, on = "stay_id", iv_abs_day := i.iv_abs_day]
# ★ 必须先算真实 tx_act 轨迹
daily[, tx_act := as.integer(is.finite(iv_abs_day) & day >= iv_abs_day)]

# stay 级 arm
daily[, stay_arm := fifelse(any(tx_init == 1L), 1L, 0L), by = stay_id]

# cross stays: arm=0 者有 cross IV 且 IV 日 <= grace_days
cross_ids <- daily[stay_arm == 0L & is.finite(iv_abs_day) & iv_abs_day <= GRACE_DAYS,
                   unique(stay_id)]
cat(sprintf("[grace] cross stays in 12h window: %d\n", length(cross_ids)))
cat(sprintf("[grace] total arm=0 stays: %d\n",
            daily[stay_arm == 0L, uniqueN(stay_id)]))
cat(sprintf("[debug] 列数=%d, 有 tx_act=%s\n", ncol(daily), "tx_act" %in% names(daily)))
cat("[debug] 列名:", paste(names(daily), collapse=","), "\n")

# ★ 手动 excusal: 这些 cross stay 的 tx_act := 0 (假装从未 cross)
cat("[mark] before set\n")
set(daily, which(daily$stay_id %in% cross_ids), "tx_act", 0L)
cat("[mark] after set\n")
cat(sprintf("[tx_act] 0=%d, 1=%d, NA=%d\n",
            sum(daily$tx_act == 0, na.rm=TRUE),
            sum(daily$tx_act == 1, na.rm=TRUE),
            sum(is.na(daily$tx_act))))

# ★ 关键: 手动 prune + 删 NA 列
# SEQTaRget 内部 prune 会漏掉 iv_abs_day (492K 行 NA) → 必须先删
needed <- c("stay_id", "day", "eligible", "tx_act", "tx_init", "oc_death",
            DAG14, "trial", "arm", "origin_day", "followup")
needed <- unique(needed[needed %in% names(daily)])
daily <- daily[, needed, with = FALSE]

# 诊断 NA 列
na_by_col <- sapply(names(daily), function(c) sum(is.na(daily[[c]])))
na_cols <- names(daily)[na_by_col > 0]
cat(sprintf("[prune] 保留 %d 列, NA 行 = %d\n",
            length(needed), sum(!complete.cases(daily))))
if (length(na_cols) > 0) {
  cat("[prune] 有 NA 的列:\n")
  print(na_by_col[na_cols])
}

# 删所有含 NA 的行 (SEQTaRget 要求 complete.cases)
n_before <- nrow(daily)
daily <- daily[complete.cases(daily)]
cat(sprintf("[na.omit] 删 %d 行 NA → 剩 %d 行\n",
            n_before - nrow(daily), nrow(daily)))

t_data <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("[data] prepared in %.1fs · rows=%d stays=%d\n",
            t_data, nrow(daily), daily[, uniqueN(stay_id)]))

# ---- SEQuential: censoring + weighted + boot ----
cat("[SEQTaRget] fitting censoring + weighted + boot=200 ...\n")
t1 <- Sys.time()

opts <- SEQopts(hazard = TRUE,
                weighted = TRUE, weight.preexpansion = TRUE, weight.p99 = TRUE,
                bootstrap = TRUE, bootstrap.nboot = NBOOT,
                bootstrap.sample = 1,
                followup.max = TAU,
                seed = SEED + IM_ID,
                bootstrap.CI_method = "se")  # ★ 注意: excused=FALSE (默认)

mod <- SEQuential(daily, id.col = "stay_id", time.col = "day",
                  eligible.col = "eligible", treatment.col = "tx_act",
                  outcome.col = "oc_death", time_varying.cols = DAG14,
                  fixed.cols = list(), method = "censoring",
                  options = opts, verbose = FALSE)

t2 <- Sys.time()
t_seq <- as.numeric(difftime(t2, t1, units = "secs"))
cat(sprintf("[SEQTaRget] done in %.1fs\n", t_seq))

# ---- hazard ratio ----
h <- hazard_ratio(mod)
cat("\n════════ imp1 Grace-Window PP boot=200 ════════\n")
print(h)

loghr <- unname(log(h["Hazard ratio"]))
hr    <- as.numeric(h["Hazard ratio"])
lci   <- as.numeric(h["LCI"])
uci   <- as.numeric(h["UCI"])

# 从 CI 反推 se
se <- if (is.finite(lci) && is.finite(uci))
        (log(uci) - log(lci)) / (2 * qnorm(0.975)) else NA_real_

cat(sprintf("\n★ HR = %.4f (%.4f – %.4f)\n", hr, lci, uci))
cat(sprintf("★ logHR = %.4f  se = %.4f\n", loghr, se))
cat(sprintf("★ boot = %d\n", NBOOT))

rm(mod); invisible(gc())

# ---- 保存 ----
out <- data.table(imp = IM_ID, loghr = loghr, se = se,
                  hr = hr, lci = lci, uci = uci,
                  boot = NBOOT, grace_days = GRACE_DAYS,
                  n_cross_in_grace = length(cross_ids),
                  elapsed_sec = t_seq)
fp <- file.path(OUT_T, "v28.4_pp_grace12h_imp1_se.csv")
fwrite(out, fp)
cat(sprintf("\nsaved: %s\n", fp))
cat("结束:", format(Sys.time(), "%H:%M:%S"), "\n")
