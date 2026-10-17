# ============================================================================
# v28.4_continuous_multi.R
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
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/v28.4_multi_common.R")

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

suffix <- paste0(if (FIRST) "firstonly_" else "", TAG)
PER_F <- file.path(OUT_T, paste0("v28.4_continuous_", suffix, "_per_imp.csv"))

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

  # checkpoint per outcome
  done_imps_yn <- if (file.exists(PER_F)) {
    sub <- fread(PER_F, colClasses = list(imp = "integer"))
    sort(unique(sub[outcome == yn]$imp))
  } else integer(0)
  todo_yn <- setdiff(IMPS, done_imps_yn)
  cat(sprintf("[checkpoint %s] todo=%s done=%s\n", yn,
              if (length(todo_yn)) paste(todo_yn, collapse = ",") else "(none)",
              if (length(done_imps_yn)) paste(done_imps_yn, collapse = ",") else "(none)"))

  # ---- imp loop: 只跑 todo，O_APPEND ----
  for (im in todo_yn) {
    d <- build_enr(im)
    pt_row <- gcomp(d, ycol)   # named num: diff, logratio, mean0, mean1
    boot_sd_row <- c(boot_sd_diff = NA_real_, boot_sd_logratio = NA_real_,
                     boot_sd_mean0 = NA_real_, boot_sd_mean1 = NA_real_)
    if (NBOOT > 0 && im %in% BOOTIMPS) {
      set.seed(SEED + im)
      br <- t(vapply(seq_len(NBOOT), function(b) boot_once(d, ycol), numeric(4)))
      boot_sd_row <- setNames(apply(br, 2, sd, na.rm = TRUE),
                              c("boot_sd_diff", "boot_sd_logratio",
                                "boot_sd_mean0", "boot_sd_mean1"))
    }
    per_row <- data.table(imp = im, outcome = yn,
                          diff = pt_row["diff"], logratio = pt_row["logratio"],
                          mean0 = pt_row["mean0"], mean1 = pt_row["mean1"],
                          boot_sd_diff = boot_sd_row["boot_sd_diff"],
                          boot_sd_logratio = boot_sd_row["boot_sd_logratio"],
                          boot_sd_mean0 = boot_sd_row["boot_sd_mean0"],
                          boot_sd_mean1 = boot_sd_row["boot_sd_mean1"])
    fwrite(per_row, PER_F, append = file.exists(PER_F))
    rm(d, per_row); invisible(gc())
  }

  # ---- Rubin 合并（从 PER_F 读全量）----
  sub <- fread(PER_F, colClasses = list(imp = "integer"))[outcome == yn]
  stopifnot(setequal(sort(sub$imp), IMPS))
  setorder(sub, imp)
  pts <- as.matrix(sub[, .(diff, logratio, mean0, mean1)])

  # Rubin 合并 diff（Barnard-Rubin df + t）
  pool_metric <- function(colname, sd_colname, logscale = FALSE) {
    qbar <- mean(pts[, colname])
    bv <- if (m > 1) var(pts[, colname]) else 0
    se <- if (colname %in% names(sub) && any(is.finite(sub[[sd_colname]])))
      sub[[sd_colname]][which.min(abs(IMPS - BOOTIMPS[1]))]  # 取任意一个 bootstrap imp 的 sd
    else NA_real_
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
  # 取第一个有 bootstrap SE 的 imp（对齐旧版 bidx 逻辑）
  boot_idx <- which(apply(sub[, .(boot_sd_diff, boot_sd_logratio, boot_sd_mean0, boot_sd_mean1)],
                          1, function(r) all(is.finite(r))))[1]
  se_boot <- if (length(boot_idx)) {
    c(diff = sub$boot_sd_diff[boot_idx],
      logratio = sub$boot_sd_logratio[boot_idx])
  } else c(diff = NA_real_, logratio = NA_real_)

  # 直接用 Barnard-Rubin（与 pp/subgroup 同口径）
  barnard_pool <- function(vals, se_within) {
    qbar <- mean(vals); bv <- if (length(vals) > 1) var(vals) else 0
    if (is.finite(se_within)) {
      ubar <- se_within^2
      tv <- ubar + (1 + 1 / length(vals)) * bv
      dfm <- if (bv > 0) (length(vals) - 1) * (1 + ubar / ((1 + 1 / length(vals)) * bv))^2 else Inf
      zt <- if (is.finite(dfm)) qt(.975, dfm) else qnorm(.975)
      c(est = qbar,
        lci = if (is.finite(tv)) qbar - zt * sqrt(tv) else NA_real_,
        uci = if (is.finite(tv)) qbar + zt * sqrt(tv) else NA_real_,
        p = if (is.finite(tv)) 2 * (if (is.finite(dfm)) pt(-abs(qbar / sqrt(tv)), dfm)
                                      else pnorm(-abs(qbar / sqrt(tv)))) else NA_real_)
    } else c(est = qbar, lci = NA_real_, uci = NA_real_, p = NA_real_)
  }
  md <- barnard_pool(sub$diff, se_boot["diff"])
  mr <- barnard_pool(sub$logratio, se_boot["logratio"])
  out_l[[yn]] <- data.table(
    outcome = yn,
    mean_control = mean(sub$mean0), mean_treat = mean(sub$mean1),
    diff_days = md["est"], diff_lci = md["lci"], diff_uci = md["uci"], p_diff = md["p"],
    mean_ratio = exp(mr["est"]), ratio_lci = exp(mr["lci"]),
    ratio_uci = exp(mr["uci"]), p_ratio = mr["p"],
    n_imp = m, first_only = FIRST, n_boot = NBOOT)
}
res <- rbindlist(out_l)
cat("\n========== CONTINUOUS multi-enrollment ==========\n"); print(res)
POOL_F <- file.path(OUT_T, paste0("v28.4_continuous_", suffix, ".csv"))
fwrite(res, POOL_F)
cat("saved:", POOL_F, "\n")
cat("per_imp already appended into:", PER_F, "\n")
