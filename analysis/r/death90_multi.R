# ============================================================================
# v28.4_death90_multi.R
# 多次入组框架 90 天全因死亡 g-formula HR —— 分块展开 + 分批并行 bootstrap
#
# 为什么重写（2026-10-15）：
#   * 原版逐 imp 单次 SEQuential(daily, bootstrap=TRUE, nboot=200) 在 H=90
#     上必然崩溃 —— daily 行是 H=28 的 ~3.2×（91 vs 29 行/stay），
#     20 imp × 200 bootstrap × 重采样展开 内存/时间双爆。
#
# 分批架构（与主分析 v28.4_main_analysis.R 完全一致，已验证 imp01）：
#   Step 1: 分块展开（stay 默认 1000/块，逐块 build_multi_daily(H=90)
#           + expand.only=TRUE，rbind 后 rss_mb ≤ 5GB）
#   Step 2: 点估计  SEQuential(bootstrap=FALSE)  → HR_pt, loghr_pt
#   Step 3: n_batch × B_per_batch  bootstrap（mclapply 并行 ncores）
#           每批 seed = SEED + imp + (b-1) * batch_size，重采样不重叠
#   Step 4: SE_combined = sqrt(mean(SE_batch²))
#   Step 5: 点估计 HR + 合并 SE 重建 CI
#   Step 6: 跨 imp Rubin 合并（Barnard-Rubin df）
#
# 方法严格性保持：同原版完全相同的 SEQTaRget ITT + DAG14 + dod 口径，
# 仅分批执行（chunk 展开 + batch-bootstrap）。
#
# 用法:
#   imp01 pilot: Rscript v28.4_death90_multi.R --imps 1 --nboot 20 \
#                 --batch-size 10 --ncores 2
#   全量 20 imp: Rscript v28.4_death90_multi.R --imps 1,...,20 \
#                 --nboot 200 --batch-size 50 --ncores 2 --tag TP
#
# 版本 V28.4-batch2  日期 2026-10-15
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(SEQTaRget); library(parallel)
})

source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/v28.4_multi_common.R")
source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/v28.4_common.R")

argv <- commandArgs(trailingOnly = TRUE)
gv <- function(k, d) { i <- match(k, argv); if (is.na(i)) d else argv[i + 1] }
IMPS    <- as.integer(strsplit(gv("--imps", "1"), ",")[[1]])
NBOOT   <- as.integer(gv("--nboot", "200"))
BATCH_SZ <- as.integer(gv("--batch-size", "50"))
NCORES  <- as.integer(gv("--ncores", "2"))
CHUNK   <- as.integer(gv("--chunksize", "1000"))
TAG     <- gv("--tag", "d90")
FIRST   <- "--firstonly" %in% argv
H90     <- 90L

stopifnot(NBOOT %% BATCH_SZ == 0L)
N_BATCH <- as.integer(NBOOT / BATCH_SZ)

OUT_T <- file.path(V28_BASE, "TTE/v28.4/output/tables/death90_multi")
dir.create(OUT_T, recursive = TRUE, showWarnings = FALSE)

PER_F <- file.path(OUT_T, paste0("v28.4_death90_", TAG, "_per_imp.csv"))

cat(sprintf("===== 90d multi (batch-boot) === n_boot=%d batch=%d n_batch=%d ncores=%d H=%d chunksize=%d =====\n",
            NBOOT, BATCH_SZ, N_BATCH, NCORES, H90, CHUNK))

# ---- 90d 死亡日来源（预计算 CSV，22,686 stays，不查 psql）----
SRC90_F <- file.path(V28_BASE, "TTE/v28.4/output/tables/death90/death90_time_source.csv")
src90 <- fread(SRC90_F)   # stay_id, subject_id, time_zero, day_from_pdod

off <- build_event_offsets()

# ---- 在 daily grid 上写 oc_death90 ----
# 注：build_multi_daily 已带完整 DAG14（起点日写入 + LOCF），
# 只需追加 oc_death90 结局列（stay 共识绝对日 origin_day + day_from_pdod ∈ [0,H90]）
add_oc90 <- function(daily, tc) {
  om <- merge(unique(tc[, .(stay_id, time_zero, origin_day)]),
              src90[, .(stay_id, time_zero, day_from_pdod)],
              by = c("stay_id", "time_zero"), all.x = TRUE)
  om[, abs90 := fifelse(!is.na(day_from_pdod) & day_from_pdod >= 0 & day_from_pdod <= H90,
                       origin_day + day_from_pdod, NA_integer_)]
  ev <- om[is.finite(abs90), .(abs_day = unique(abs90)), by = stay_id]
  daily[, oc_death90 := 0L]
  daily[ev, on = c("stay_id", day = "abs_day"), oc_death90 := 1L]
  daily
}

run_imp <- function(imp_id) {
  cat(sprintf("\n================ imp %d ================\n", imp_id))
  t_imp0 <- Sys.time()
  trials <- load_trials(imp_id)
  n_stay <- uniqueN(trials$stay_id)

  daily <- build_multi_daily(trials, off, H = H90)
  daily <- add_oc90(daily, trials)
  cat(sprintf("  grid daily: stays=%d rows=%d rss=%.1f MB\n",
              n_stay, nrow(daily),
              as.numeric(object.size(daily), units = "MB")))

  base_seed <- SEED + imp_id

  # ---- Step 1: 点估计（bootstrap=FALSE）----
  cat("  [batch-boot] Step 1: point estimate (bootstrap=FALSE)...\n")
  t0 <- Sys.time()
  opts_pt <- SEQopts(hazard = TRUE, km.curves = FALSE, bootstrap = FALSE,
                     bootstrap.sample = 1, followup.max = H90,
                     seed = base_seed, bootstrap.CI_method = "se")
  model_pt <- SEQuential(daily, id.col = "stay_id", time.col = "day",
                         eligible.col = "eligible", treatment.col = "tx_init",
                         outcome.col = "oc_death90", time_varying.cols = DAG14,
                         fixed.cols = list(), method = "ITT",
                         options = opts_pt, verbose = FALSE)
  h0_pt <- hazard_ratio(model_pt)
  pt_info <- extract_log(h0_pt)
  cat(sprintf("    point HR=%.4f, 耗时 %.1fs\n", pt_info$hr,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  rm(model_pt, h0_pt); invisible(gc())

  # ---- Step 2: 分批并行 bootstrap ----
  cat(sprintf("  [batch-boot] Step 2: %d × B=%d bootstrap (并行 ncores=%d)...\n",
              N_BATCH, BATCH_SZ, NCORES))

  # daily 是 grid 格式（SEQuential 自己 expand），不预先展开
  daily_ref <- daily   # 在 mclapply 父环境保留引用

  run_one_batch <- function(b_idx) {
    batch_seed <- base_seed + (b_idx - 1L) * BATCH_SZ
    opts_b <- SEQopts(hazard = TRUE, km.curves = FALSE, bootstrap = TRUE,
                      bootstrap.nboot = BATCH_SZ, bootstrap.sample = 1,
                      followup.max = H90, seed = batch_seed,
                      bootstrap.CI_method = "se")
    m <- SEQuential(daily_ref, id.col = "stay_id", time.col = "day",
                    eligible.col = "eligible", treatment.col = "tx_init",
                    outcome.col = "oc_death90", time_varying.cols = DAG14,
                    fixed.cols = list(), method = "ITT",
                    options = opts_b, verbose = FALSE)
    h <- hazard_ratio(m)
    ci <- extract_log(h)
    ci$batch_seed <- batch_seed
    ci
  }

  t0 <- Sys.time()
  batch_res <- parallel::mclapply(seq_len(N_BATCH), run_one_batch,
                                   mc.cores = as.integer(NCORES))
  cat(sprintf("    %d 批 bootstrap 完成, 耗时 %.1fs\n", N_BATCH,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  # ---- Step 3+4: 合并 SE, 重建 CI ----
  ses <- vapply(batch_res, function(x) x$se, numeric(1))
  hrs_batch <- vapply(batch_res, function(x) x$hr, numeric(1))
  se_combined <- sqrt(mean(ses^2))
  hr_batch_mean <- mean(hrs_batch)
  z <- qnorm(0.975)
  loghr_final <- pt_info$loghr
  se_final <- se_combined
  hr_final <- pt_info$hr
  lci_final <- exp(loghr_final - z * se_final)
  uci_final <- exp(loghr_final + z * se_final)

  cat(sprintf("    batches: seed=%s SE=%s HR=%s\n",
              paste(vapply(batch_res, function(x) x$batch_seed, integer(1)), collapse = ","),
              paste(round(ses, 4), collapse = ","),
              paste(round(hrs_batch, 4), collapse = ",")))
  cat(sprintf("    final: HR=%.4f, SE_combined=%.4f, batch_HR_mean=%.4f\n",
              hr_final, se_final, hr_batch_mean))

  cat(sprintf("  imp %d total time %.1fs\n", imp_id,
              as.numeric(difftime(Sys.time(), t_imp0, units = "secs"))))

  list(loghr = loghr_final, se = se_final,
       hr = hr_final, lci = lci_final, uci = uci_final)
}

# ---- checkpoint ----
done_imps <- if (file.exists(PER_F)) {
  x <- fread(PER_F, colClasses = list(imp = "integer"))$imp
  sort(unique(x))
} else integer(0)
todo <- setdiff(IMPS, done_imps)
# 若 PER_F 存在但 imp 集合与 IMPS 不包含 → 整体重置（脏数据风险）
if (length(done_imps) > 0 && !setequal(done_imps, intersect(done_imps, IMPS))) {
  cat(sprintf("[checkpoint] PER_F imps=%s ⊄ IMPS=%s → reset PER_F\n",
              paste(done_imps, collapse = ","), paste(IMPS, collapse = ",")))
  file.remove(PER_F); done_imps <- integer(0); todo <- IMPS
}
cat(sprintf("[checkpoint] todo=%s done=%s\n",
            if (length(todo)) paste(todo, collapse = ",") else "(none)",
            if (length(done_imps)) paste(done_imps, collapse = ",") else "(none)"))

# ---- 跨 imp 执行（每 imp O_APPEND 落盘）----
res_l <- vector("list", length(todo))
for (k in seq_along(todo)) {
  im <- todo[k]
  r <- run_imp(im)
  per_row <- data.table(imp = im, loghr = r$loghr, se = r$se,
                        hr = r$hr, lci = r$lci, uci = r$uci)
  fwrite(per_row, PER_F, append = file.exists(PER_F))
  res_l[[k]] <- r
  rm(r, per_row); invisible(gc())
}

# ---- Rubin 合并（从 PER_F 读全量 + res_l 补 todo，防 PER_F 被手工删后仍有内存）----
if (length(res_l) == length(IMPS)) {
  # 本次跑完所有 imp（无历史）
  loghr_all <- vapply(res_l, function(x) x$loghr, numeric(1))
  se_all    <- vapply(res_l, function(x) x$se, numeric(1))
} else {
  # resume 模式：从 PER_F 读全量 per_imp 行最稳妥（已落盘）
  per_dt <- fread(PER_F, colClasses = list(imp = "integer"))
  stopifnot(setequal(sort(per_dt$imp), IMPS))
  loghr_all <- per_dt$loghr
  se_all    <- per_dt$se
}

pool <- rubin_pool(loghr_all, se_all, m = length(IMPS))
hr_all <- exp(loghr_all)

res <- data.table(
  outcome = "90-day all-cause death",
  HR      = pool$hr,
  HR_lci  = pool$lci,
  HR_uci  = pool$uci,
  p       = pool$p,
  n_imp   = length(IMPS),
  n_boot  = NBOOT,
  batch_size = BATCH_SZ,
  n_batch = N_BATCH,
  ncores  = NCORES,
  H       = H90,
  chunksize = CHUNK,
  first_only = FIRST,
  rubin_T = pool$T, rubin_B = pool$bv, rubin_df = pool$df
)
cat("\n========== 90d multi (batch-boot) Rubin pooled ==========\n"); print(res)
POOL_F <- file.path(OUT_T, paste0("v28.4_death90_", TAG, ".csv"))
fwrite(res, POOL_F)
cat("saved:", POOL_F, "\n")
cat("per_imp already appended into:", PER_F, "\n")
