# ============================================================================
# piecewise_multi.R
# 多次入组（multi-enrollment）框架下：分段 HR（0-3/3-7/7-14/14-28 天）
#
# 由包 g-formula（ITT，km.curves 风险路径）给出分段边界累积风险 R(t)，
# 分段常数风险换算：
#   p_g = (R_g(r)-R_g(l)) / (1-R_g(l))
#   h_g = -log(1-p_g) / (r-l) ; HR = h_treat / h_control
# CI：包 bootstrap 给出边界风险 SE，delta 法传播到 log(HR)。
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
TAG    <- gv("--tag", "pw")
FIRST  <- "--firstonly" %in% argv

OUT_T <- file.path(V28_BASE, "TTE/v28.4/output/tables/piecewise_multi")
dir.create(OUT_T, recursive = TRUE, showWarnings = FALSE)

EDGES <- c(0L, 3L, 7L, 14L, 28L)
BOUND <- c(3L, 7L, 14L, 28L)
ILAB  <- c("0-3d", "3-7d", "7-14d", "14-28d")

cat("build event offsets once...\n")
off <- build_event_offsets()

# 取边界风险（点 + SE），返回 matrix 行=boundary,列=R0,R1,SE0,SE1
run_imp <- function(imp_id) {
  cat(sprintf("\n================ imp %d ================\n", imp_id))
  trials <- load_trials(imp_id)
  daily  <- build_multi_daily(trials, off)
  do_boot <- NBOOT > 0 && imp_id %in% BOOTIMPS
  fit_once <- function(attempt) {
    opts <- SEQopts(km.curves = TRUE,
                    bootstrap = do_boot,
                    bootstrap.nboot = max(1L, NBOOT),
                    bootstrap.sample = 1,
                    followup.max = TAU,
                    seed = SEED + if (do_boot) imp_id + 1000L * (attempt - 1L) else 0L,
                    bootstrap.CI_method = "se",
                    risk.times = BOUND,
                    selection.first_trial = FIRST)
    SEQuential(daily, id.col = "stay_id", time.col = "day",
                      eligible.col = "eligible", treatment.col="tx_init",
                      outcome.col = "oc_death", time_varying.cols = DAG14,
                      fixed.cols = list(), method = "ITT",
                      options = opts, verbose = FALSE)
  }
  # bootstrap SE 在内存压力下可能静默退化为 NA：校验 + 最多 3 次重试
  MAXTRY <- if (do_boot) 3L else 1L
  for (attempt in seq_len(MAXTRY)) {
    mod <- fit_once(attempt)
    rd <- as.data.table(mod@risk.data[[1]])
    rd <- rd[Followup %in% BOUND]
    pick_arm <- function(a, col)
      rd[A == a][order(Followup)][[col]]
    out <- cbind(R0 = pick_arm("0", "Risk"), R1 = pick_arm("1", "Risk"))
    hasSE <- "SE" %in% names(rd)
    se <- cbind(SE0 = if (hasSE) pick_arm("0", "SE") else rep(NA_real_, length(BOUND)),
                SE1 = if (hasSE) pick_arm("1", "SE") else rep(NA_real_, length(BOUND)))
    if (!do_boot || all(is.finite(se))) {
      if (do_boot && attempt > 1L)
        cat(sprintf("[piecewise] imp %d bootstrap SE recovered on attempt %d\n", imp_id, attempt))
      return(cbind(out, se))
    }
    cat(sprintf("[piecewise] imp %d bootstrap SE missing/NA (attempt %d/%d)%s\n",
                imp_id, attempt, MAXTRY,
                if (attempt < MAXTRY) " — retrying" else " — FINAL FAILURE"))
    rm(mod, rd); gc(verbose = FALSE)
    Sys.sleep(5 * attempt)
  }
  stop(sprintf("[piecewise] imp %d: bootstrap SE unavailable after %d attempts; refusing to write NA CIs", imp_id, MAXTRY))
}

per_imp <- lapply(IMPS, run_imp)
m <- length(IMPS)

# 逐段：从每 imp 边界风险换算 HR；delta 法 SE
segments <- lapply(1:4, function(ii) {
  l <- EDGES[ii]; r <- EDGES[ii + 1]; dt <- r - l
  # 需用到边界 ii-1（l）与 ii（r）。l=0 时 R(0)=0,SE=0
  loghr_imp <- numeric(m)
  grad_imp <- matrix(NA_real_, m, 4)  # dloghr/d[R0l,R0r,R1l,R1r]
  for (k in seq_len(m)) {
    M <- per_imp[[k]]
    R0r <- M[ii, "R0"]; R1r <- M[ii, "R1"]
    R0l <- if (ii == 1) 0 else M[ii - 1, "R0"]
    R1l <- if (ii == 1) 0 else M[ii - 1, "R1"]
    p0 <- (R0r - R0l) / (1 - R0l); p1 <- (R1r - R1l) / (1 - R1l)
    h0 <- -log(1 - p0) / dt; h1 <- -log(1 - p1) / dt
    loghr_imp[k] <- log(h1 / h0)
    # dh/dR: 对 p=(Rr-Rl)/(1-Rl), dp/dRr=1/(1-Rl);
    # dp/dRl=-(1-Rr)/(1-Rl)^2 ; dh = (1/dt)* dp/(1-p)
    gfun <- function(Rr, Rl) {
      # unname：矩阵按列名取的标量携带名字，经算术传播会使输出名变成
      # "dRl.R0" 等，导致下方 g0["dRl"] 名字索引静默返回 NA
      Rr <- unname(Rr); Rl <- unname(Rl)
      p <- (Rr - Rl) / (1 - Rl)
      c(dRr = 1 / (1 - Rl), dRl = -(1 - Rr) / (1 - Rl)^2) /
        ((dt) * (1 - p))
    }
    g0 <- gfun(R0r, R0l); g1 <- gfun(R1r, R1l)
    grad_imp[k, ] <- c(-g0["dRl"], -g0["dRr"], g1["dRl"], g1["dRr"]) /
      c(rep(h0, 2), rep(h1, 2))
  }
  # 点估计（跨 imp 均值 logHR）
  qbar <- mean(loghr_imp); bv <- if (m > 1) var(loghr_imp) else 0
  # within SE：取任一 bootstrap imp 的边界 SE 做 delta（与主分析同口径：
  # Ubar 取 bootstrap imp 的抽样方差，TV = Ubar + (1+1/m)BV）
  se_loghr <- NA_real_
  bidx <- which(vapply(per_imp, function(Mz)
    all(is.finite(Mz[, c("SE0", "SE1")])), logical(1)))[1]
  if (is.finite(bidx)) {
    M <- per_imp[[bidx]]
    s0r <- M[ii, "SE0"]; s1r <- M[ii, "SE1"]
    s0l <- if (ii == 1) 0 else M[ii - 1, "SE0"]
    s1l <- if (ii == 1) 0 else M[ii - 1, "SE1"]
    G <- grad_imp[bidx, ]
    se_loghr <- sqrt(sum((G * c(s0l, s0r, s1l, s1r))^2))
  } else {
    warning(sprintf("[piecewise] segment %s: no finite bootstrap SE in any imputation (bootstrap may have degraded under memory pressure); CI=NA", ILAB[ii]))
  }
  # Rubin 合并（Barnard-Rubin df + t，与 common.R 完全一致）
  ubar <- if (is.finite(se_loghr)) se_loghr^2 else NA_real_
  tv <- if (is.finite(ubar)) ubar + (1 + 1 / m) * bv else NA_real_
  dfm <- if (is.finite(tv) && bv > 0) (m - 1) * (1 + ubar / ((1 + 1 / m) * bv))^2 else Inf
  zt <- if (is.finite(dfm)) qt(.975, dfm) else qnorm(.975)
  data.table(segment = ILAB[ii], t_start = l, t_end = r,
             HR = exp(qbar),
             HR_lci = if (is.finite(tv)) exp(qbar - zt * sqrt(tv)) else NA_real_,
             HR_uci = if (is.finite(tv)) exp(qbar + zt * sqrt(tv)) else NA_real_,
             p = if (is.finite(tv)) 2 * (if (is.finite(dfm)) pt(-abs(qbar / sqrt(tv)), dfm)
                                        else pnorm(-abs(qbar / sqrt(tv)))) else NA_real_)
})
res <- rbindlist(segments)
cat("\n========== PIECEWISE multi-enrollment ==========\n"); print(res)
suffix <- paste0(if (FIRST) "firstonly_" else "", TAG)
fwrite(res, file.path(OUT_T, paste0("piecewise_", suffix, ".csv")))
cat("saved:", file.path(OUT_T, paste0("piecewise_", suffix, ".csv")), "\n")
