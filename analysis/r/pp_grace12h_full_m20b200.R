# ============================================================================
# v28.4_pp_grace12h_full_m20b200.R
# 跑全 20 imps × boot=200 Grace-Window PP
# 手动 excusal 绕开 SEQTaRget excused=TRUE bug
# 每个 imp 存 loghr, se (从 boot CI 反推), hr, lci, uci
# ============================================================================
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/v28.4_multi_common.R")
suppressPackageStartupMessages({library(data.table); library(SEQTaRget)})
setDTthreads(2L); Sys.setenv(OMP_NUM_THREADS = 2L)
SEED <- 20260629L; TAU <- 28L; GRACE_DAYS <- 1L; M <- 20L; NBOOT <- 200L
AUG_F <- file.path(V28_BASE, "TTE/v26/output/sequential_trials_augmented.parquet")

OUT_T  <- file.path(V28_BASE, "TTE/v28.4/output/tables/sensitivity")
OUT_LOG <- file.path(V28_BASE, "TTE/v28.4/output/logs")
dir.create(OUT_T, recursive = TRUE, showWarnings = FALSE)

cat("════════ Grace-Window PP 24h FULL m=20 boot=200 ════════\n")
cat("手动 excusal, grace_days =", GRACE_DAYS, "\n")
cat("开始:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

build_iv_day <- function(trials) {
  aug <- as.data.table(read_parquet(AUG_F))
  aug <- aug[, .(stay_id, time_zero, itt_group, fivh = first_iv_hours_from_tz)]
  omap <- unique(trials[, .(stay_id, time_zero, origin_day, arm)])
  z <- merge(aug, omap, by = c("stay_id", "time_zero"), all.x = TRUE)
  z[, ivd := fifelse(itt_group == "treatment" | arm == 1L, origin_day,
                     fifelse(is.finite(fivh),
                             origin_day + as.integer(ceiling(fivh/24)),
                             NA_integer_))]
  z[, .(iv_abs_day = min(ivd, na.rm = TRUE)), by = stay_id]
}

cat("[offsets] building event offsets...\n")
off <- build_event_offsets()

all_res <- vector("list", M)
for (imp in seq_len(M)) {
  cat(sprintf("\n─── imp %2d/%d ───\n", imp, M))
  t0 <- Sys.time()
  trials <- load_trials(imp)
  daily  <- build_multi_daily(trials, off, H = TAU)
  ivmap  <- build_iv_day(trials)
  daily[ivmap, on = "stay_id", iv_abs_day := i.iv_abs_day]
  daily[, tx_act := as.integer(is.finite(iv_abs_day) & day >= iv_abs_day)]
  daily[, stay_arm := fifelse(any(tx_init == 1L), 1L, 0L), by = stay_id]

  cross_ids <- daily[stay_arm == 0L & is.finite(iv_abs_day) & iv_abs_day <= GRACE_DAYS,
                     unique(stay_id)]
  set(daily, which(daily$stay_id %in% cross_ids), "tx_act", 0L)

  needed <- c("stay_id", "day", "eligible", "tx_act", "tx_init", "oc_death",
              DAG14[DAG14 %in% names(daily)])
  daily <- daily[, unique(needed), with = FALSE][complete.cases(daily)]

  opts <- SEQopts(hazard = TRUE, weighted = TRUE, weight.preexpansion = TRUE,
                  weight.p99 = TRUE, bootstrap = TRUE, bootstrap.nboot = NBOOT,
                  bootstrap.sample = 1, followup.max = TAU,
                  seed = SEED + imp)
  t1 <- Sys.time()
  mod <- SEQuential(daily, id.col = "stay_id", time.col = "day",
                    eligible.col = "eligible", treatment.col = "tx_act",
                    outcome.col = "oc_death",
                    time_varying.cols = DAG14[DAG14 %in% names(daily)],
                    fixed.cols = list(), method = "censoring",
                    options = opts, verbose = FALSE)
  t_seq <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

  h <- hazard_ratio(mod)
  loghr <- unname(log(h["Hazard ratio"]))
  hr    <- as.numeric(h["Hazard ratio"])
  lci   <- as.numeric(h["LCI"])
  uci   <- as.numeric(h["UCI"])
  se <- if (is.finite(lci) && is.finite(uci))
          (log(uci) - log(lci)) / (2 * qnorm(0.975)) else NA_real_

  all_res[[imp]] <- data.table(imp=imp, loghr=loghr, se=se, hr=hr,
                                lci=lci, uci=uci, boot=NBOOT,
                                n_cross_in_grace=length(cross_ids),
                                elapsed_sec=t_seq)
  cat(sprintf("  HR=%.4f (%.4f-%.4f)  loghr=%.4f  se=%.4f  [%ds]\n",
              hr, lci, uci, loghr, se, round(t_seq)))
  rm(mod, daily, trials, ivmap); invisible(gc())
}

res <- rbindlist(all_res)
fp <- file.path(OUT_T, "v28.4_pp_grace12h_per_imp_m20b200.csv")
fwrite(res, fp)
cat(sprintf("\n✅ 全 20 imps 完成, 保存: %s\n", fp))
cat("完成:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(sprintf("总耗时: %.0f min\n", sum(res$elapsed_sec) / 60))
print(summary(res$loghr))
