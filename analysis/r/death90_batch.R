#!/usr/bin/env Rscript
# ============================================================================
# v28.4_death90_batch.R
# 90 天死亡分析 —— 分批并行 bootstrap 方案
#
# 方法: stay-level, DAG14 估 PS → 稳定化 IPTW (1/99 截断) →
#       加权 Cox (cluster=subject_id) → HR
#       m=20 Rubin 合并点估计 + 分批并行 B=200 bootstrap SE
#
# 分批架构 (与主分析一致, 已验证 imp01):
#   Step 1: 点估计 HR (单次加权 Cox, 无 bootstrap)
#   Step 2: n_batch × batch_size 次 bootstrap (mclapply 并行)
#           每批 seed 递增, 保证重采样不重叠
#   Step 3: SE_combined = sqrt(mean(SE_batch²))
#   Step 4: 用点估计 HR + 合并 SE 重建 CI
#   Step 5: 跨 imp 用 Rubin 规则合并
#
# 用法:
#   imp01 pilot: Rscript v28.4_death90_batch.R --imps 1 --nboot 20 --batch-size 10
#   全量 20 imp:  Rscript v28.4_death90_batch.R --imps 1,2,...,20 --nboot 200
#
# 须以 LD_LIBRARY_PATH=/opt/miniforge3/lib 运行
# V28.4-batch  2026-10-13
# ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(arrow); library(survival); library(parallel)
})

BASE   <- "/data/小米云盘/文档/论文/血磷"
AUG_F  <- file.path(BASE, "TTE/v26/output/sequential_trials_augmented.parquet")
MICE_F <- file.path(BASE, "sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet")
SRC_F  <- file.path(BASE, "TTE/v28.4/output/tables/death90/death90_time_source.csv")
OUT_T  <- file.path(BASE, "TTE/v28.4/output/tables/death90")
dir.create(OUT_T, recursive = TRUE, showWarnings = FALSE)
TAU  <- 90L
SEED <- 20260629L

# ---- DAG14 ----
BIN  <- c("sepsis3","invasive_vent_at_t0","vasoactive_at_t0","refeeding_risk")
CONT <- c("age","baseline_sofa","charlson","egfr","baseline_phosphate",
          "lactate_t0","corrected_calcium_mmol_min")
FACT <- c("gender","admission_type","nutrition_support_type")
LEV <- list(gender=c("F","M"),
            admission_type=c("Observation","Emergency","Surgical_Elective"),
            nutrition_support_type=c("none","enteral","parenteral","mixed"))
DAG <- c(CONT, BIN, FACT)

# ---- CLI ----
argv <- commandArgs(trailingOnly = TRUE)
gv <- function(k,d){i<-match(k,argv);if(is.na(i))d else argv[i+1]}
IMPS <- as.integer(strsplit(gv("--imps","1"),",")[[1]])
NBOOT     <- as.integer(gv("--nboot","200"))
BATCH_SZ  <- as.integer(gv("--batch-size", "50"))
NCORES    <- as.integer(gv("--ncores", "4"))
TAG       <- gv("--tag", "batch")
B_NBATCH  <- NBOOT / BATCH_SZ
stopifnot(NBOOT %% BATCH_SZ == 0L)

cat(sprintf("===== 90d death batch === n_boot=%d, batch_size=%d, n_batch=%d, ncores=%d =====\n",
            NBOOT, BATCH_SZ, B_NBATCH, NCORES))

# ---- 核心: 加载 stay 唯一队列 ----
src <- fread(SRC_F)[,.(stay_id, d90day = day_from_pdod)]

load_persons <- function(imp_id) {
  aug <- as.data.table(read_parquet(AUG_F))
  imp <- as.data.table(read_parquet(sprintf(MICE_F, imp_id)))
  a <- aug[order(stay_id, trial_num), .SD[1], by = stay_id]
  m <- imp[order(stay_id, trial_num), .SD[1], by = stay_id]
  stopifnot(identical(a$stay_id, m$stay_id))
  p <- a[, .(stay_id, subject_id,
             arm = as.integer(a$itt_group == "treatment"))]
  p <- merge(p, src, by = "stay_id", sort = FALSE)
  p[, evt  := as.integer(!is.na(d90day) & d90day >= 0 & d90day <= TAU)]
  p[, time := fifelse(evt == 1L, d90day, TAU)]
  for (v in c(CONT, BIN)) p[, (v) := m[[v]]]
  for (v in BIN)          p[, (v) := as.integer(as.character(get(v)))]
  for (v in FACT)         p[, (v) := factor(a[[v]], levels = LEV[[v]])]
  p
}

# ---- IPTW + Cox ----
add_weight <- function(p) {
  fml <- as.formula(paste("arm ~", paste(DAG, collapse = " + ")))
  ps  <- predict(suppressWarnings(glm(fml, data = p, family = binomial(),
                                       control = list(maxit = 100))),
                  type = "response")
  ps[is.na(ps)] <- mean(p$arm)
  pa <- mean(p$arm)
  w  <- fifelse(p$arm == 1, pa / ps, (1 - pa) / (1 - ps))
  p[, wgt := pmin(pmax(w, quantile(w, .01)), quantile(w, .99))]
  p
}

cox_hr <- function(p) {
  fit <- coxph(Surv(time, evt) ~ arm, data = p, weights = wgt,
               cluster = subject_id, robust = TRUE,
               control = coxph.control(iter.max = 50))
  s <- summary(fit)
  c(beta = s$coefficients[1, "coef"],
    se   = s$coefficients[1, "se(coef)"])
}

risk_contrast <- function(p) {
  d <- p[, .(risk = sum(wgt * evt) / sum(wgt)), by = arm][order(-arm)]
  rt <- d[arm == 1]$risk; rc <- d[arm == 0]$risk
  c(risk_t = rt, risk_c = rc, rd = rt - rc, rr = rt / rc,
    or = (rt / (1 - rt)) / (rc / (1 - rc)))
}

# ---- 分批 bootstrap: 单次 imp ----
run_one_imp_batch <- function(imp_id) {
  cat(sprintf("\n================ imp %d ================\n", imp_id))
  t0 <- Sys.time()
  p <- load_persons(imp_id)
  cat(sprintf("  stays=%d, treated=%d, death90=%d (%.1f%%)\n",
              nrow(p), sum(p$arm), sum(p$evt), mean(p$evt) * 100))

  # Step 1: 点估计 (IPTW + Cox, 无 bootstrap)
  cat("  [Step 1] point estimate (IPTW + Cox)... ")
  p_w <- add_weight(p)
  cx  <- cox_hr(p_w)
  pt  <- risk_contrast(p_w)
  pt_hr <- exp(cx["beta"])
  cat(sprintf("HR=%.4f, risk_t=%.4f, risk_c=%.4f, 耗时%.1fs\n",
              pt_hr, pt["risk_t"], pt["risk_c"],
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  # Step 2: 分批 bootstrap (mclapply 并行)
  cat(sprintf("  [Step 2] %d × B=%d bootstrap (并行 ncores=%d)...\n",
              B_NBATCH, BATCH_SZ, NCORES))

  run_one_batch <- function(b_idx) {
    batch_seed <- SEED + (imp_id - 1L) * NBOOT + (b_idx - 1L) * BATCH_SZ
    set.seed(batch_seed)
    nb <- 0L
    Bm <- matrix(NA_real_, BATCH_SZ, 6)
    colnames(Bm) <- c("risk_t", "risk_c", "rd", "rr", "or", "beta")
    while (nb < BATCH_SZ) {
      pick <- sample(unique(p$subject_id), uniqueN(p$subject_id), replace = TRUE)
      pb <- p[data.table(subject_id = pick), on = "subject_id",
              allow.cartesian = TRUE]
      r <- tryCatch({
        pw <- add_weight(pb)
        list(rk = risk_contrast(pw), cx = cox_hr(pw))
      }, error = function(e) NULL)
      if (is.null(r) || any(!is.finite(r$rk), !is.finite(r$cx[1]))) next
      nb <- nb + 1L
      Bm[nb, ] <- c(r$rk, r$cx["beta"])
    }
    betas_batch <- Bm[, "beta"]
    # 该批 beta 均值 + SE
    beta_mean <- mean(betas_batch, na.rm = TRUE)
    se_batch  <- sd(betas_batch, na.rm = TRUE)
    list(
      batch_idx = b_idx,
      batch_seed = batch_seed,
      beta_mean = beta_mean,
      se_batch  = se_batch,
      # 各 risk metric 的 mean + quantiles（可选，诊断用）
      risk_t_mean = mean(Bm[, "risk_t"], na.rm = TRUE),
      risk_c_mean = mean(Bm[, "risk_c"], na.rm = TRUE),
      rd_mean     = mean(Bm[, "rd"],     na.rm = TRUE),
      rr_mean     = mean(Bm[, "rr"],     na.rm = TRUE),
      or_mean     = mean(Bm[, "or"],     na.rm = TRUE)
    )
  }

  t_boot0 <- Sys.time()
  batch_res <- mclapply(seq_len(B_NBATCH), run_one_batch,
                        mc.cores = NCORES)
  t_boot1 <- Sys.time()
  cat(sprintf("    %d 批 bootstrap 完成, 耗时 %.1fs\n", B_NBATCH,
              as.numeric(difftime(t_boot1, t_boot0, units = "secs"))))

  # Step 3: 合并 SE (方差平均)
  ses   <- vapply(batch_res, function(x) x$se_batch,    numeric(1))
  betas <- vapply(batch_res, function(x) x$beta_mean, numeric(1))
  se_combined <- sqrt(mean(ses^2))
  beta_batch_mean <- mean(betas)

  # Step 4: 点估计 HR + 合并 SE 重建 CI
  z <- qnorm(0.975)
  hr_final  <- pt_hr
  lci_final <- exp(cx["beta"] - z * se_combined)
  uci_final <- exp(cx["beta"] + z * se_combined)

  cat(sprintf("  [Step 3+4] final: HR=%.4f, SE_combined=%.4f, CI=(%.4f, %.4f)\n",
              hr_final, se_combined, lci_final, uci_final))

  list(
    loghr = cx["beta"],       # 点估计 logHR
    se    = se_combined,      # 合并 SE
    hr    = hr_final,
    lci   = lci_final,
    uci   = uci_final,
    # risk metrics (点估计, 可选合并 bootstrap)
    risk_t = pt["risk_t"], risk_c = pt["risk_c"],
    rd = pt["rd"], rr = pt["rr"], or = pt["or"],
    # 诊断
    info = list(
      n_boot = NBOOT, batch_size = BATCH_SZ, n_batch = B_NBATCH, ncores = NCORES,
      batch_seeds = vapply(batch_res, function(x) x$batch_seed, integer(1)),
      batch_SEs   = ses,
      beta_pt     = cx["beta"],
      beta_batch_mean = beta_batch_mean,
      duration    = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )
  )
}

# ---- Rubin 合并 ----
rubin_pool <- function(loghr, se, m = length(loghr)) {
  qbar <- mean(loghr); ubar <- mean(se^2)
  bv <- if (m > 1) var(loghr) else 0
  tv <- ubar + bv * (1 + 1 / m)
  dfm <- if (m > 1 && bv > 0)
    (m - 1) * (1 + ubar / (bv * (1 + 1 / m)))^2 else Inf
  zt <- if (is.finite(dfm)) qt(0.975, dfm) else qnorm(0.975)
  list(loghr = qbar, hr = exp(qbar),
       lci = exp(qbar - zt * sqrt(tv)),
       uci = exp(qbar + zt * sqrt(tv)),
       p = 2 * (if (is.finite(dfm)) pt(-abs(qbar / sqrt(tv)), dfm)
                else pnorm(-abs(qbar / sqrt(tv)))),
       T = tv, B = bv, df = dfm)
}

# ---- 主循环 ----
set.seed(SEED)
per_imp <- lapply(IMPS, run_one_imp_batch)
loghr <- vapply(per_imp, function(x) x$loghr, numeric(1))
sev   <- vapply(per_imp, function(x) x$se,    numeric(1))

pooled <- rubin_pool(loghr, sev, m = length(IMPS))
cat("\n========== RUBIN POOLED (m =", length(IMPS), ") ==========\n")
cat(sprintf("90-day death HR = %.4f (%.4f - %.4f), p = %.4f\n",
            pooled$hr, pooled$lci, pooled$uci, pooled$p))

# ---- 落盘 ----
suffix <- paste0("_batch", BATCH_SZ, "_m", length(IMPS))

per_dt <- data.frame(
  imp = IMPS,
  loghr = loghr, se = sev,
  hr    = vapply(per_imp, function(x) x$hr,  numeric(1)),
  lci   = vapply(per_imp, function(x) x$lci, numeric(1)),
  uci   = vapply(per_imp, function(x) x$uci, numeric(1)),
  risk_t = vapply(per_imp, function(x) x$risk_t, numeric(1)),
  risk_c = vapply(per_imp, function(x) x$risk_c, numeric(1)),
  rd     = vapply(per_imp, function(x) x$rd,     numeric(1)),
  rr     = vapply(per_imp, function(x) x$rr,     numeric(1)),
  or_val = vapply(per_imp, function(x) x$or,     numeric(1))
)
write.csv(per_dt,
          file.path(OUT_T, paste0("v28.4_death90", suffix, ".csv")),
          row.names = FALSE)

pool_dt <- data.frame(
  analysis = "90-day death (DAG14 stabilized IPTW, batch bootstrap)",
  version = "V28.4-batch", date = format(Sys.Date()),
  m = length(IMPS), n_boot = NBOOT, batch_size = BATCH_SZ, ncores = NCORES,
  hr = pooled$hr, lci = pooled$lci, uci = pooled$uci,
  p = pooled$p, T = pooled$T, B = pooled$B, df = pooled$df
)
write.csv(pool_dt,
          file.path(OUT_T, paste0("v28.4_death90_pooled", suffix, ".csv")),
          row.names = FALSE)

cat("\nresults saved to:\n")
cat(" ", file.path(OUT_T, paste0("v28.4_death90", suffix, ".csv")), "\n")
cat(" ", file.path(OUT_T, paste0("v28.4_death90_pooled", suffix, ".csv")), "\n")
