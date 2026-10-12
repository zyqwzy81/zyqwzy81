# ============================================================================
# continuous_multi.R
# 多次入组（multi-enrollment）框架下：连续有界结局的加权 fractional g-computation
#
# 结局（均 0..28，死亡记 0；Montez-Rath et al., Ann Intern Med 2024 口径）：
#   * VFD-28          呼吸机无事件日（M0_07，逐起点）
#   * ICU-free days   28 天 ICU 无事件日 = 死亡记0；否则 28-min(ICU LOS,28)
#
# 方法：在全部入组（每起点一行，n=32,328）上对 y/28 做 quasibinomial logit
#   fractional 回归（arm + DAG14），分别置 arm=0/1 预测标准化均值；
#   报告 mean difference（天）与 mean ratio。stay 级 bootstrap（重抽 stay，
#   天然处理同一 stay 多起点的聚集）。MI：点估计跨 m=20，bootstrap 在 imp1
#   提供抽样方差 + Rubin imp 间方差。
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
TAG    <- gv("--tag", "cont")
FIRST  <- "--firstonly" %in% argv

OUT_T <- file.path(V28_BASE, "TTE/v28.4/output/tables/continuous_multi")
dir.create(OUT_T, recursive = TRUE, showWarnings = FALSE)

# ---- 逐起点连续结局（M0_07，不随 imp 变）----
m07 <- as.data.table(read_parquet(file.path(RRES, "M0_07_outcome_extraction.parquet")))
m07 <- m07[, .(stay_id, time_zero, death28 = as.integer(death_28day),
               vfd = as.numeric(vfd_28), iculos = as.numeric(icu_los_from_t0_days))]
m07[, icufree := fifelse(death28 == 1L, 0,
                         28 - pmin(pmax(round(iculos), 0), 28))]
YCOLS <- c(vfd = "vfd", icufree = "icufree")

# 组装某次 imp 的逐起点分析表
build_enr <- function(imp_id) {
  trials <- load_trials(imp_id)
  d <- merge(trials, m07, by = c("stay_id", "time_zero"), all.x = TRUE)
  if (FIRST) d <- d[, .SD[1], by = stay_id]
  d
}

# fractional g-computation：返回 diff/ratio/m0/m1
gcomp <- function(d, ycol, covs = DAG14) {
  df <- as.data.frame(d)
  df$s <- df[[ycol]] / 28
  fit <- suppressWarnings(glm(
    as.formula(paste("s ~", paste(c("arm", covs), collapse = "+"))),
    data = df, family = quasibinomial()))
  d0 <- df; d0$arm <- 0L
  d1 <- df; d1$arm <- 1L
  q0 <- mean(predict(fit, d0, type = "response"))
  q1 <- mean(predict(fit, d1, type = "response"))
  c(diff = 28 * (q1 - q0), logratio = log(q1 / q0),
    mean0 = 28 * q0, mean1 = 28 * q1)
}

# stay 级 bootstrap 一次重抽后的估计
boot_once <- function(d, ycol, covs = DAG14) {
  stays <- unique(d$stay_id)
  pick <- sample(stays, length(stays), replace = TRUE)
  b <- data.table(stay_id = pick, bid = seq_along(pick))
  dd <- merge(b, d, by = "stay_id", allow.cartesian = TRUE)
  gcomp(dd, ycol, covs)
}

m <- length(IMPS)
out_l <- list()
for (yn in names(YCOLS)) {
  ycol <- YCOLS[[yn]]
  cat(sprintf("\n========== outcome %s ==========\n", yn))
  pts <- matrix(NA_real_, m, 4)
  boot_sd <- NULL
  for (k in seq_len(m)) {
    imp_id <- IMPS[k]
    d <- build_enr(imp_id)
    pts[k, ] <- gcomp(d, ycol)
    if (NBOOT > 0 && imp_id %in% BOOTIMPS) {
      set.seed(SEED + imp_id)
      br <- t(vapply(seq_len(NBOOT), function(b) boot_once(d, ycol),
                     numeric(4)))
      boot_sd <- apply(br, 2, sd, na.rm = TRUE)
    }
  }
  colnames(pts) <- c("diff", "logratio", "mean0", "mean1")

  # Rubin 合并 diff（与 common.R 完全一致：Barnard-Rubin df + t）
  pool_metric <- function(colname, logscale = FALSE) {
    qbar <- mean(pts[, colname])
    bv <- if (m > 1) var(pts[, colname]) else 0
    se <- if (!is.null(boot_sd)) boot_sd[colname] else NA_real_
    if (is.finite(se)) {
      ubar <- se^2
      tv <- ubar + (1 + 1 / m) * bv
      dfm <- if (bv > 0) (m - 1) * (1 + ubar / ((1 + 1 / m) * bv))^2 else Inf
      zt <- if (is.finite(dfm)) qt(.975, dfm) else qnorm(.975)
      lci <- qbar - zt * sqrt(tv)
      uci <- qbar + zt * sqrt(tv)
      p <- 2 * (if (is.finite(dfm)) pt(-abs(qbar / sqrt(tv)), dfm)
                else pnorm(-abs(qbar / sqrt(tv))))
    } else { lci <- uci <- p <- NA_real_ }
    rv <- c(qbar, lci, uci, p)
    names(rv) <- c("est", "lci", "uci", "p")
    rv
  }
  md <- pool_metric("diff")
  mr <- pool_metric("logratio")
  out_l[[yn]] <- data.table(
    outcome = yn,
    mean_control = mean(pts[, "mean0"]), mean_treat = mean(pts[, "mean1"]),
    diff_days = md["est"], diff_lci = md["lci"], diff_uci = md["uci"], p_diff = md["p"],
    mean_ratio = exp(mr["est"]), ratio_lci = exp(mr["lci"]),
    ratio_uci = exp(mr["uci"]), p_ratio = mr["p"],
    n_imp = m, first_only = FIRST, n_boot = NBOOT)
}
res <- rbindlist(out_l)
cat("\n========== CONTINUOUS multi-enrollment ==========\n"); print(res)
suffix <- paste0(if (FIRST) "firstonly_" else "", TAG)
fwrite(res, file.path(OUT_T, paste0("continuous_", suffix, ".csv")))
cat("saved:", file.path(OUT_T, paste0("continuous_", suffix, ".csv")), "\n")
