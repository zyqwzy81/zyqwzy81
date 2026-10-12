#!/usr/bin/env Rscript
# ============================================================================
# v28.4_secondary_multi.R
# 多次入组（multi-enrollment）序贯框架下的全部 28 天二分类结局
#
# 方法：SEQTaRget ITT，km.curves=TRUE（g-formula 标准化），bootstrap=TRUE
#       每 imp 构建一次宽日级表，循环各结局（切换 outcome.col）。
#       逐 imp 抽取：每臂 28d 风险、RR、RD 及各自 SE；跨 imp 用
#       Rubin 规则合并（logRR 尺度；RD 用方差合成）。
# 单次入组：--firstonly 时启用 selection.first_trial=TRUE（包内置）。
#
# 结局列见 v28.4_multi_common.R 的 BIN_OUTCOMES。
# 用法:
#   验证: Rscript v28.4_secondary_multi.R --imps 1 --nboot 0
#   全量: Rscript v28.4_secondary_multi.R --imps 1..20 --nboot 200
#   单次: Rscript ... --firstonly
# 须 LD_LIBRARY_PATH=/opt/miniforge3/lib   V28.4  2026-10-10
# ============================================================================
suppressPackageStartupMessages(library(SEQTaRget))
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/v28.4_multi_common.R")

# ---- 参数 ----
argv <- commandArgs(trailingOnly = TRUE)
gv <- function(k, d) { i <- match(k, argv); if (is.na(i)) d else argv[i + 1] }
IMPS   <- as.integer(strsplit(gv("--imps", "1"), ",")[[1]])
NBOOT  <- as.integer(gv("--nboot", "200"))
# 需 bootstrap 的 imp（默认全部）；省时方案：点估计跨全部 imp，
# bootstrap 仅在少数 imp（如 --bootimps 1），抽样方差 + Rubin imp间方差
DEF_BOOT <- if (NBOOT > 0) paste(IMPS, collapse = ",") else ""
BOOTIMPS <- as.integer(strsplit(gv("--bootimps", DEF_BOOT), ",")[[1]])
TAG    <- gv("--tag", "multi")
FIRST  <- "--firstonly" %in% argv
OCS <- BIN_OUTCOMES
if ("--outcomes" %in% argv)
  OCS <- strsplit(gv("--outcomes", paste(BIN_OUTCOMES, collapse = ",")), ",")[[1]]

OUT_T <- file.path(V28_BASE, "TTE/v28.4/output/tables/secondary_multi")
dir.create(OUT_T, recursive = TRUE, showWarnings = FALSE)

# ---- 逐 imp × 结局 抽取 ----
cat("build event offsets once...\n")
off <- build_event_offsets()

extract_at_tau <- function(mod) {
  rd_all <- as.data.table(mod@risk.data[[1]])
  rc_all <- as.data.table(mod@risk.comparison[[1]])
  # Followup 取 ==TAU；若无则取最小值（不同 risk.times 解析）
  fu <- TAU
  if (!fu %in% rd_all$Followup) fu <- min(rd_all$Followup)
  rd <- rd_all[Followup == fu]
  r0 <- rd[A == "0"]; r1 <- rd[A == "1"]
  rc <- rc_all[Followup == fu & A_x == "risk_0" & A_y == "risk_1"]
  # 无 bootstrap 时 SE/CI 列不存在，安全取 NA
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

run_imp <- function(imp_id) {
  cat(sprintf("\n================ imp %d ================\n", imp_id))
  trials <- load_trials(imp_id)
  daily <- build_multi_daily(trials, off)
  check_multi_daily(daily, paste0("imp", imp_id))

  do_boot <- NBOOT > 0 && imp_id %in% BOOTIMPS
  opts <- SEQopts(km.curves = TRUE,
                  bootstrap = do_boot,
                  bootstrap.nboot = max(1L, NBOOT),
                  bootstrap.sample = 1,
                  followup.max = TAU,
                  seed = SEED + if (do_boot) imp_id else 0L,
                  bootstrap.CI_method = "se",
                  risk.times = TAU,
                  selection.first_trial = FIRST)
  res <- list()
  for (oc in OCS) {
    cat(sprintf("  -- outcome %s\n", oc))
    mod <- SEQuential(daily, id.col = "stay_id", time.col = "day",
                      eligible.col = "eligible", treatment.col = "tx_init",
                      outcome.col = oc, time_varying.cols = DAG14,
                      fixed.cols = list(), method = "ITT",
                      options = opts, verbose = FALSE)
    res[[oc]] <- tryCatch(extract_at_tau(mod),
                          error = function(e) { cat("    ERR:", conditionMessage(e), "\n"); NULL })
  }
  res
}

per_imp <- lapply(IMPS, run_imp)
names(per_imp) <- IMPS

# ---- Rubin 合并 ----
m <- length(IMPS)
pool_one <- function(oc) {
  rows <- lapply(per_imp, function(x) x[[oc]])
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  g <- function(k) vapply(rows, function(r) r[[k]], numeric(1))
  logrr <- log(g("rr"))
  rr_within <- g("logrr_se"); rd_within <- g("rd_se")
  has_se <- any(is.finite(rr_within))
  rr_b <- if (length(logrr) > 1) var(logrr) else 0
  rd_b <- if (length(g("rd")) > 1) var(g("rd")) else 0
  if (has_se) {
    rr_tv <- mean(rr_within^2, na.rm = TRUE) + (1 + 1 / length(rows)) * rr_b
    rd_tv <- mean(rd_within^2, na.rm = TRUE) + (1 + 1 / length(rows)) * rd_b
    df_rr <- if (rr_b > 0) (length(rows) - 1) *
        (1 + mean(rr_within^2, na.rm = TRUE) / ((1 + 1 / length(rows)) * rr_b))^2 else Inf
  } else { rr_tv <- NA_real_; rd_tv <- NA_real_; df_rr <- Inf }
  zt <- if (is.finite(df_rr)) qt(0.975, df_rr) else qnorm(0.975)
  qbar <- mean(logrr); rdbar <- mean(g("rd"))
  risk0 <- mean(g("risk0")); risk1 <- mean(g("risk1"))
  pval <- if (has_se)
    2 * (if (is.finite(df_rr)) pt(-abs(qbar / sqrt(rr_tv)), df_rr)
         else pnorm(-abs(qbar / sqrt(rr_tv)))) else NA_real_
  data.table(outcome = oc,
             risk_control = risk0, risk_treat = risk1,
             RR = exp(qbar),
             RR_lci = if (has_se) exp(qbar - zt * sqrt(rr_tv)) else NA_real_,
             RR_uci = if (has_se) exp(qbar + zt * sqrt(rr_tv)) else NA_real_,
             RD = rdbar,
             RD_lci = if (has_se) rdbar - zt * sqrt(rd_tv) else NA_real_,
             RD_uci = if (has_se) rdbar + zt * sqrt(rd_tv) else NA_real_,
             p_RR = pval,
             n_imp = length(rows), first_only = FIRST)
}
pooled <- rbindlist(lapply(OCS, pool_one))
cat("\n========== POOLED multi-enrollment (m=", m,
    if (FIRST) ", FIRST TRIAL ONLY", ") ==========\n", sep = "")
print(pooled)

suffix <- paste0(if (FIRST) "firstonly_" else "", TAG)
fwrite(pooled, file.path(OUT_T, paste0("v28.4_secondary_", suffix, ".csv")))
cat("\nsaved:", file.path(OUT_T, paste0("v28.4_secondary_", suffix, ".csv")), "\n")
