# ============================================================================
# binary_continuous.R
# stay 唯一队列：二分类次要结局（死亡/CV/输血）+ 连续结局（VFD/LOS）
# 稳定化 IPTW + m=20 Rubin + B=200 bootstrap（百分位 CI）
# 每个 (imp, boot) 只拟合一次 PS/权重，所有结局共用。
# AKI 已在 common.R 重算；90 天以 patients.dod 核验（4,799）。
# 版本: V28.4  日期: 2026-10-10
# ============================================================================
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/stay_common.R")

OUTDIR <- file.path(V28_BASE, "TTE/v28.4/output/tables/nonprimary_redo")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---- 90 天死亡：patients.dod ----
dod_tsv <- tempfile(fileext = ".tsv")
qry <- "SELECT subject_id, dod::date FROM mimiciv_hosp.patients WHERE dod IS NOT NULL"
# 口令不写入脚本，从环境变量 MIMIC_PG_PASSWORD 读取（未设置则交给 .pgpass/信任连接）
pgpw <- Sys.getenv("MIMIC_PG_PASSWORD", unset = "")
envpre <- if (nzchar(pgpw)) paste0("PGPASSWORD=", shQuote(pgpw), " ") else ""
system(paste(envpre, "psql -h localhost -p 5432 -U postgres -d mimiciv",
             "-At -F '|' -c", shQuote(qry, type = "sh"),
             ">", shQuote(dod_tsv, type = "sh")))
dod <- unique(fread(dod_tsv, sep = "|", col.names = c("subject_id", "dod")))
dod[, subject_id := as.integer(subject_id)]

cat("loading", M_IMP, "imputations...\n")
DS <- lapply(seq_len(M_IMP), function(k) {
  d <- load_stay(k)
  d <- merge(d, dod, by = "subject_id", all.x = TRUE)
  d[, d90 := as.integer(as.Date(dod) - as.Date(time_zero))]
  d[, death90_verified := as.integer(!is.na(d90) & d90 >= 0 & d90 <= 90)]
  d[, `:=`(dod = NULL, d90 = NULL)]
  d
})

bin_specs <- list(
  c("death_icu", "ICU death"),
  c("death_hospital", "Hospital death"),
  c("death90_verified", "90-day death (verified)"),
  c("cardiovascular_event_3comp", "CV composite (3-comp)"),
  c("myocardial_injury", "Myocardial injury"),
  c("new_arrhythmia_rhythm", "New arrhythmia"),
  c("transfusion_at_t0", "Transfusion"))
lin_specs <- list(
  c("vfd_28", "VFD-28 (days)"),
  c("icu_los_from_t0_days", "ICU LOS from t0 (days)"))

# 对给定 d,w，一次性返回全部结局点估计（logor / beta）+ 风险
fit_all <- function(d, w) {
  b <- vapply(bin_specs, function(sp) fit_bin(d, sp[1], w)$logor, numeric(1))
  l <- vapply(lin_specs, function(sp) fit_lin(d, sp[1], w)$beta, numeric(1))
  r <- lapply(bin_specs, function(sp) {
    rr <- weighted_risk(d, sp[1], w)
    c(rr$risk1, rr$risk0)
  })
  c(b, l, unlist(r))
}

NB <- length(bin_specs); NL <- length(lin_specs)
ncol_tot <- NB + NL + 2 * NB

# ---- 点估计：每 imp 一次 ----
pt_mat <- matrix(NA_real_, M_IMP, ncol_tot)
for (k in seq_len(M_IMP)) {
  w <- stab_weights(DS[[k]])
  pt_mat[k, ] <- fit_all(DS[[k]], w)
}

# ---- bootstrap（共享权重）----
set.seed(SEED)
n <- nrow(DS[[1]])
boot_arr <- matrix(NA_real_, B_BOOT * M_IMP, ncol_tot)
ridx <- 1L
for (b in seq_len(B_BOOT)) {
  ix <- sample.int(n, n, replace = TRUE)
  for (k in seq_len(M_IMP)) {
    dd <- DS[[k]][ix]
    boot_arr[ridx, ] <- fit_all(dd, stab_weights(dd)); ridx <- ridx + 1L
  }
  if (b %% 25 == 0) cat(sprintf("boot %d/%d\n", b, B_BOOT))
}

# ---- 汇总二分类：用 bootstrap 估计 SE 做 Rubin p，CI 用百分位 ----
assemble_bin <- function(j, sp) {
  est <- mean(pt_mat[, j])
  se_boot <- sd(boot_arr[, j])
  z <- est / se_boot
  ci <- quantile(boot_arr[, j], c(.025, .975), na.rm = TRUE)
  # 风险（点估计，跨 imp 均值）
  r1 <- mean(pt_mat[, NB + NL + 2 * (j - 1) + 1])
  r0 <- mean(pt_mat[, NB + NL + 2 * (j - 1) + 2])
  data.table(label = sp[2], OR = exp(est), LCI = exp(ci[1]), UCI = exp(ci[2]),
             p = 2 * pnorm(-abs(z)),
             risk_treat = r1, risk_control = r0,
             RD_pct = 100 * (r1 - r0), RR = r1 / r0)
}
res_bin <- rbindlist(Map(assemble_bin, seq_len(NB), bin_specs))

assemble_lin <- function(jj, sp) {
  j <- NB + jj
  est <- mean(pt_mat[, j]); se_boot <- sd(boot_arr[, j])
  ci <- quantile(boot_arr[, j], c(.025, .975), na.rm = TRUE)
  data.table(label = sp[2], beta = est, LCI = ci[1], UCI = ci[2],
             p = 2 * pnorm(-abs(est / se_boot)))
}
res_lin <- rbindlist(Map(assemble_lin, seq_len(NL), lin_specs))

fwrite(res_bin, file.path(OUTDIR, "binary_OR.csv"))
fwrite(res_lin, file.path(OUTDIR, "continuous_beta.csv"))
cat("\n=== Binary ===\n"); print(res_bin)
cat("\n=== Continuous ===\n"); print(res_lin)
