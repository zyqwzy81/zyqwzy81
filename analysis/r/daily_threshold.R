# ============================================================================
# daily_threshold.R
# stay 唯一队列：逐日 HR (Day1-14，正确条件风险集) + 血磷阈值剂量反应 (6 切点)
# 稳定化 IPTW + coxph 稳健 SE + Rubin m=20（探索性，无 bootstrap）
# 版本: V28.4  日期: 2026-10-10
# ============================================================================
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/stay_common.R")
OUTDIR <- file.path(V28_BASE, "TTE/v28.4/output/tables/nonprimary_redo")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("loading", M_IMP, "imputations...\n")
DS <- lapply(seq_len(M_IMP), load_stay)
WS <- lapply(DS, stab_weights)

clean_X <- function(X, event) {
  keep <- rep(TRUE, ncol(X))
  for (j in seq_len(ncol(X))) {
    x <- X[, j]
    if (is.na(var(x)) || var(x) == 0) { keep[j] <- FALSE; next }
    if (all(x %in% c(0, 1))) {
      if (sum(event[x == 1]) == 0 || sum(event[x == 0]) == 0) keep[j] <- FALSE
    }
  }
  X[, keep, drop = FALSE]
}
rubin <- function(q, se) {
  m <- length(q); qbar <- mean(q); ubar <- mean(se^2)
  bv <- if (m > 1) var(q) else 0
  tv <- ubar + bv * (1 + 1 / m)
  dfm <- if (m > 1 && bv > 0) (m - 1) * (1 + ubar / (bv * (1 + 1 / m)))^2 else Inf
  zt <- if (is.finite(dfm)) qt(.975, dfm) else qnorm(.975)
  p <- 2 * (if (is.finite(dfm)) pt(-abs(qbar / sqrt(tv)), dfm)
            else pnorm(-abs(qbar / sqrt(tv))))
  list(q = qbar, L = qbar - zt * sqrt(tv), U = qbar + zt * sqrt(tv), p = p)
}

# ---- 逐日 HR：区间 (d-1,d] 的条件风险集 ----
daily_one <- function(d, w, day) {
  st <- d$surv_time; ev <- d$death_28day_num
  risk <- st > (day - 1)
  tc <- pmin(st[risk], day)
  ec <- as.integer(ev[risk] == 1 & st[risk] <= day)
  X <- clean_X(dag_matrix(d[risk], DAG14), ec)
  df <- as.data.frame(X); df$arm <- d$arm[risk]
  df$.t <- tc; df$.e <- ec
  covn <- setdiff(colnames(df), c("arm", ".t", ".e"))
  cx <- coxph(as.formula(paste("Surv(.t,.e)~ arm +",
                               paste(covn, collapse = "+"))),
              data = df, weights = w[risk], ties = "efron", robust = TRUE,
              control = coxph.control(iter.max = 50))
  s <- summary(cx)$coefficients
  rn <- if (is.matrix(s)) s["arm", ] else s
  list(loghr = unname(rn["coef"]), se = unname(rn["se(coef)"]),
       n_risk = sum(risk), n_event = sum(ec))
}

res_daily <- rbindlist(lapply(1:14, function(day) {
  q <- numeric(M_IMP); se <- numeric(M_IMP); nr <- 0L; ne <- 0L
  for (k in seq_len(M_IMP)) {
    r <- daily_one(DS[[k]], WS[[k]], day)
    q[k] <- r$loghr; se[k] <- r$se
    if (k == 1L) { nr <- r$n_risk; ne <- r$n_event }
  }
  rp <- rubin(q, se)
  data.table(day = day, n_risk = nr, n_event = ne,
             hr = exp(rp$q), ci_lower = exp(rp$L), ci_upper = exp(rp$U),
             p_value = rp$p,
             bonf_sig = ifelse(rp$p < 0.05 / 14, "p<0.0036", "NS"))
}))

# ---- 阈值剂量反应：6 切点 ----
cutoffs <- c(0.40, 0.50, 0.60, 0.65, 0.70, 0.80)
cov_no_p <- setdiff(DAG14, "baseline_phosphate")

thr_one <- function(d, w, cut) {
  ix <- d$baseline_phosphate < cut
  dd <- d[ix]; ww <- w[ix]
  X <- clean_X(dag_matrix(dd, cov_no_p), dd$death_28day_num)
  df <- as.data.frame(X); df$arm <- dd$arm
  df$.t <- dd$surv_time; df$.e <- dd$death_28day_num
  covn <- setdiff(colnames(df), c("arm", ".t", ".e"))
  cx <- coxph(as.formula(paste("Surv(.t,.e)~ arm +",
                               paste(covn, collapse = "+"))),
              data = df, weights = ww, ties = "efron", robust = TRUE,
              control = coxph.control(iter.max = 50))
  s <- summary(cx)$coefficients
  rn <- if (is.matrix(s)) s["arm", ] else s
  list(loghr = unname(rn["coef"]), se = unname(rn["se(coef)"]),
       n = sum(ix), n_event = sum(dd$death_28day_num))
}

res_thr <- rbindlist(lapply(cutoffs, function(cut) {
  q <- numeric(M_IMP); se <- numeric(M_IMP); nn <- 0L; ne <- 0L
  for (k in seq_len(M_IMP)) {
    r <- thr_one(DS[[k]], WS[[k]], cut)
    q[k] <- r$loghr; se[k] <- r$se
    if (k == 1L) { nn <- r$n; ne <- r$n_event }
  }
  rp <- rubin(q, se)
  data.table(cutoff = cut, n = nn, n_event = ne,
             hr = exp(rp$q), ci_lower = exp(rp$L), ci_upper = exp(rp$U),
             p_value = rp$p,
             significant = ifelse(rp$p < 0.05, "p<0.05", "NS"))
}))

fwrite(res_daily, file.path(OUTDIR, "daily_hr.csv"))
fwrite(res_thr,   file.path(OUTDIR, "threshold.csv"))
cat("\n=== Daily HR ===\n"); print(res_daily)
cat("\n=== Threshold ===\n"); print(res_thr)
