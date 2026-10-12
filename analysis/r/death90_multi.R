# ============================================================================
# death90_multi.R
# 多次入组（multi-enrollment）框架下：90 天全因死亡 g-formula HR
#
# 90 天死亡以 mimiciv_hosp.patients.dod 核验（与 binary_continuous 一致）：
#   day90 = floor(dod - time_zero)，0..90 内为事件；
#   在 H=90 的宽日级时间线上新增 oc_death90 列并在该绝对日铺点。
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
TAG    <- gv("--tag", "d90")
FIRST  <- "--firstonly" %in% argv
H90    <- 90L

OUT_T <- file.path(V28_BASE, "TTE/v28.4/output/tables/death90_multi")
dir.create(OUT_T, recursive = TRUE, showWarnings = FALSE)

# ---- dod（日期）----
dod_tsv <- tempfile(fileext = ".tsv")
qry <- "SELECT subject_id, dod::date FROM mimiciv_hosp.patients WHERE dod IS NOT NULL"
pgpw <- Sys.getenv("MIMIC_PG_PASSWORD", unset = "")
envpre <- if (nzchar(pgpw)) paste0("PGPASSWORD=", shQuote(pgpw), " ") else ""
system(paste(envpre, "psql -h localhost -p 5432 -U postgres -d mimiciv",
             "-At -F '|' -c", shQuote(qry, type = "sh"),
             ">", shQuote(dod_tsv, type = "sh")))
dod <- unique(fread(dod_tsv, sep = "|", col.names = c("subject_id", "dod")))
dod[, subject_id := as.integer(subject_id)]

# stay → subject 映射（trials 不含 subject_id，从 augmented 取）
AUG_F <- file.path(V28_BASE, "TTE/v26/output/sequential_trials_augmented.parquet")
subj_map <- unique(as.data.table(read_parquet(AUG_F))[, .(stay_id, subject_id)])

cat("build event offsets once...\n")
off <- build_event_offsets()

# 逐 origin 的 90 天死亡绝对事件日（不随 imp；origin_day 随几何）
# build 时先取 origin 时间线
base_orig <- function(trials) {
  o <- unique(trials[, .(stay_id, time_zero, origin_day)])
  o <- merge(o, subj_map, by = "stay_id", all.x = TRUE)
  o <- merge(o, dod, by = "subject_id", all.x = TRUE)
  o[, d90 := as.integer(as.Date(dod) - as.Date(time_zero))]
  o[, abs90 := fifelse(is.finite(d90) & d90 >= 0 & d90 <= H90,
                       origin_day + d90, NA_integer_)]
  # 同 stay 多起点分别计算
  o[, .(stay_id, time_zero, abs90)]
}

run_imp <- function(imp_id) {
  cat(sprintf("\n================ imp %d ================\n", imp_id))
  trials <- load_trials(imp_id)
  daily  <- build_multi_daily(trials, off, H = H90)
  # 新增 90d 结局列
  daily[, oc_death90 := 0L]
  omap <- merge(unique(trials[, .(stay_id, time_zero)]),
                base_orig(trials), by = c("stay_id", "time_zero"), all.x = TRUE)
  ev <- omap[is.finite(abs90), .(abs_day = unique(abs90)), by = stay_id]
  daily[ev, on = c("stay_id", day = "abs_day"), oc_death90 := 1L]

  do_boot <- NBOOT > 0 && imp_id %in% BOOTIMPS
  opts <- SEQopts(hazard = TRUE,
                  bootstrap = do_boot,
                  bootstrap.nboot = max(1L, NBOOT),
                  bootstrap.sample = 1,
                  followup.max = H90,
                  seed = SEED + if (do_boot) imp_id else 0L,
                  bootstrap.CI_method = "se",
                  selection.first_trial = FIRST)
  mod <- SEQuential(daily, id.col = "stay_id", time.col = "day",
                    eligible.col = "eligible", treatment.col = "tx_init",
                    outcome.col = "oc_death90", time_varying.cols = DAG14,
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
res <- data.table(outcome = "90-day all-cause death",
                  HR = exp(qbar),
                  HR_lci = if (has_se) exp(qbar - zt * sqrt(tv)) else NA_real_,
                  HR_uci = if (has_se) exp(qbar + zt * sqrt(tv)) else NA_real_,
                  p = pval, n_imp = m, first_only = FIRST, n_boot = NBOOT)
cat("\n========== 90d multi-enrollment ==========\n"); print(res)
suffix <- paste0(if (FIRST) "firstonly_" else "", TAG)
fwrite(res, file.path(OUT_T, paste0("death90_", suffix, ".csv")))
cat("saved:", file.path(OUT_T, paste0("death90_", suffix, ".csv")), "\n")
