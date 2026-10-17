#!/usr/bin/env Rscript
# ============================================================================
# v28.4_death90_seq20.R
# 90d Death · SEQTaRget 序贯多次入组 · 8块路线 A · imp1-20 · 点估计
#
# 方法学 (28d 已验证):
#   8 块 stay-level expand.only → 合并 pp → 手动 params → internal.analysis + internal.hazard
#   与全量 SEQuential 金标准 HR 差异 0.000% (28d, 不含 trial/period 局部编号问题)
#   90d daily 2.08M 行 → 分块峰值 ~4.5GB (不 OOM)
#   路线 A imp01 HR=1.016 vs stay-level IPTW Cox HR=1.003, 差异 1.3% ✅
#
# SE 策略: 点估计无 bootstrap, Rubin 合并时用 stay-level Cox 的 SE≈0.044
#         (已验证两种方法 SE 相近, null 区域内一致)
#
# 用法: LD_LIBRARY_PATH=/opt/miniforge3/lib Rscript v28.4_death90_seq20.R 1,2,...,20
# 须先 source v28.4_common.R (含 load_trials, DAG14, SEED, na_like, zoo_locf)
# ============================================================================
suppressPackageStartupMessages({library(data.table); library(SEQTaRget); library(zoo)})

CODE_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r"
source(file.path(CODE_DIR, "v28.4_common.R"))
ns_env <- asNamespace("SEQTaRget")
factorize_fn <- get("factorize", envir=ns_env)
ia <- get("internal.analysis", envir=ns_env)
ih <- get("internal.hazard", envir=ns_env)
init_cache <- get("init_formula_cache", envir=ns_env)
clean_fastglm <- get("clean_fastglm", envir=ns_env)

argv <- commandArgs(trailingOnly = TRUE)
IMPS <- as.integer(strsplit(argv[1], ",")[[1]])
BASE <- "/data/小米云盘/文档/论文/血磷"
OUT_DIR <- file.path(BASE, "TTE/v28.4/output/tables/death90")
OUT_CSV <- file.path(OUT_DIR, "v28.4_death90_seq20.csv")
H90 <- 90L
CHUNK_SIZE <- 2835   # 8 块, 已验证

# ---- src90: 90d 绝对死亡日 (从 CSV, 不查 psql) ----
src90 <- fread(file.path(BASE, "TTE/v28.4/output/tables/death90/death90_time_source.csv"))
rss_mb <- function() as.numeric(system("ps -o rss= -p $PPID", intern=TRUE))/1024

# ---- 已完成的 imp 跳过 ----
done_ids <- integer(0)
if (file.exists(OUT_CSV)) {
  done_ids <- fread(OUT_CSV)$imp
}
to_run <- setdiff(IMPS, done_ids)
cat(sprintf("===== 90d Death · SEQTaRget 路线 A · 待跑 %s =====\n",
            if (length(to_run)) paste(to_run, collapse=",") else "(all done)"))

for (imp_id in to_run) {
  cat(sprintf("\n===== imp %02d · start =====\n", imp_id))
  t_imp0 <- Sys.time()
  trials <- load_trials(imp_id)

  # --- daily90 (照抄 build_daily_input 逻辑) ---
  TAU_OLD <- TAU; TAU <- H90

  origin_map <- unique(trials[, .(stay_id, origin_day)])
  src90_i <- merge(copy(src90), origin_map, by="stay_id", all.x=TRUE)
  src90_i[, abs_death_day := fifelse(
    !is.na(day_from_pdod) & day_from_pdod >= 0 & day_from_pdod <= H90,
    origin_day + day_from_pdod, NA_integer_)]
  death90_map <- src90_i[!is.na(abs_death_day), .(abs_death_day = min(abs_death_day)), by = stay_id]

  grid <- trials[, .(max_origin = max(origin_day)), by = stay_id]
  grid[, last_day := max_origin + H90]
  daily90 <- grid[, .(day = seq(0L, last_day[1])), by = stay_id]
  setorder(daily90, stay_id, day)

  daily90[, outcome := 0L]
  daily90[death90_map, on = c("stay_id", day = "abs_death_day"), outcome := 1L]

  orig <- unique(trials[, .(stay_id, origin_day, arm)])
  daily90[, eligible := 0L]; daily90[, tx_init := 0L]
  daily90[orig, on = c(stay_id = "stay_id", day = "origin_day"),
          `:=`(eligible = 1L, tx_init = i.arm)]

  for (v in DAG14) {
    daily90[, (v) := na_like(trials[[v]])]
    daily90[trials, on = c(stay_id = "stay_id", day = "origin_day"),
            (v) := get(paste0("i.", v))]
    setorder(daily90, stay_id, day)
    daily90[, (v) := zoo_locf(get(v)), by = stay_id]
  }
  setorder(daily90, stay_id, day)
  TAU <- TAU_OLD

  n_stay <- uniqueN(daily90$stay_id)
  n_dead <- daily90[outcome==1L, uniqueN(stay_id)]
  dup <- daily90[, .N, by = .(stay_id, day)][N > 1]
  cat(sprintf("  daily90: %d rows, %d stays, 90d死亡=%d, 重复=%d\n",
              nrow(daily90), n_stay, n_dead, nrow(dup)))

  # --- 8 块 expand.only ---
  opts_expand <- SEQopts(hazard=TRUE, bootstrap=FALSE, bootstrap.sample=1,
                          followup.max=H90, seed=SEED, expand.only=TRUE)
  stay_ids <- unique(daily90$stay_id)
  n_chunk <- length(stay_ids) %/% CHUNK_SIZE + 1

  t_exp0 <- Sys.time()
  chunks <- vector("list", n_chunk)
  for (ci in seq_len(n_chunk)) {
    sv <- stay_ids[((ci-1)*CHUNK_SIZE+1):min(ci*CHUNK_SIZE, length(stay_ids))]
    chunk <- daily90[stay_id %in% sv]
    chunks[[ci]] <- SEQuential(as.data.frame(chunk), id.col="stay_id", time.col="day",
                                eligible.col="eligible", treatment.col="tx_init",
                                outcome.col="outcome", time_varying.cols=DAG14,
                                fixed.cols=list(), method="ITT",
                                options=opts_expand, verbose=FALSE)
    rm(chunk); gc(verbose=FALSE)
  }
  pp90 <- as.data.table(rbindlist(chunks)); rm(chunks); gc()
  t_exp1 <- Sys.time()
  cat(sprintf("  8块 expand: %d rows, %.1fs\n",
              nrow(pp90), as.numeric(difftime(t_exp1,t_exp0,units="secs"))))

  # --- params 模板 (小样本 daily90, 可跨 imp 复用! 只在第一次构建) ---
  if (!exists("params_template", inherits = FALSE)) {
    cat(sprintf("  构建 params 模板 (首次)...\n"))
    TAU <- H90
    dt_tpl <- copy(trials[stay_id %in% stay_ids[1:2000]])
    grid_tpl <- dt_tpl[, .(max_origin = max(origin_day)), by = stay_id]
    grid_tpl[, last_day := max_origin + H90]
    daily_tpl <- grid_tpl[, .(day = seq(0L, last_day[1])), by = stay_id]
    setorder(daily_tpl, stay_id, day)
    daily_tpl[, outcome := 0L]
    daily_tpl[death90_map, on = c("stay_id", day = "abs_death_day"), outcome := 1L]
    orig_tpl <- unique(dt_tpl[, .(stay_id, origin_day, arm)])
    daily_tpl[, eligible := 0L]; daily_tpl[, tx_init := 0L]
    daily_tpl[orig_tpl, on = c(stay_id = "stay_id", day = "origin_day"),
              `:=`(eligible = 1L, tx_init = i.arm)]
    for (v in DAG14) {
      daily_tpl[, (v) := na_like(dt_tpl[[v]])]
      daily_tpl[dt_tpl, on = c(stay_id = "stay_id", day = "origin_day"),
                (v) := get(paste0("i.", v))]
      setorder(daily_tpl, stay_id, day)
      daily_tpl[, (v) := zoo_locf(get(v)), by = stay_id]
    }
    setorder(daily_tpl, stay_id, day)
    opts_tpl <- SEQopts(hazard=TRUE, bootstrap=FALSE, followup.max=H90, seed=SEED)
    m_tpl <- SEQuential(as.data.frame(daily_tpl), id.col="stay_id", time.col="day",
                         eligible.col="eligible", treatment.col="tx_init", outcome.col="outcome",
                         time_varying.cols=DAG14, fixed.cols=list(),
                         method="ITT", options=opts_tpl, verbose=FALSE)
    assign("params_template", slot(m_tpl, "params"), envir = .GlobalEnv)
    rm(m_tpl, dt_tpl, grid_tpl, daily_tpl, orig_tpl); gc()
    TAU <- TAU_OLD
    cat(sprintf("    模板 OK\n"))
  }

  # --- 手动 pipeline ---
  params_90 <- get("params_template", envir = .GlobalEnv)
  slot(params_90, "expand.only") <- FALSE
  slot(params_90, "DT") <- factorize_fn(pp90, params_90)
  slot(params_90, "data") <- factorize_fn(as.data.table(daily90), params_90)
  slot(params_90, "followup.max") <- H90
  slot(params_90, "treat.level") <- list(0, 1)
  slot(params_90, "hazard") <- TRUE
  slot(params_90, "bootstrap") <- FALSE
  slot(params_90, "weighted") <- FALSE

  analytic <- ia(params_90)
  clean_model <- lapply(analytic[[1]]$model, function(sg) {
    sg$model <- clean_fastglm(sg$model); sg
  })
  hz <- ih(clean_model, params_90, init_cache(params_90))
  hr <- hz[["Hazard ratio"]]

  # --- 写盘 (append) ---
  row <- data.frame(imp = imp_id, hr = hr)
  if (!file.exists(OUT_CSV)) {
    write.csv(row, OUT_CSV, row.names = FALSE)
  } else {
    write.table(row, OUT_CSV, append = TRUE, sep = ",",
                row.names = FALSE, col.names = FALSE)
  }

  elapsed <- as.numeric(difftime(Sys.time(), t_imp0, units = "secs"))
  cat(sprintf("  ✅ imp %02d · HR=%.4f · %.1fs · RSS=%.0f MB\n",
              imp_id, hr, elapsed, rss_mb()))
  gc(verbose = FALSE)
}

# ---- 全部 imp 完成 → Rubin pooled ----
all <- fread(OUT_CSV)
cat(sprintf("\n===== RUBIN POOLED (m=%d, SEQTaRget 路线 A) =====\n", nrow(all)))
# SE 从 stay-level IPTW Cox 借 (已验证两种方法 null 区域一致)
iptw <- fread(file.path(OUT_DIR, "v28.4_death90_batch50_m20.csv"))
se_ref <- mean(iptw$se)   # within SE ≈ 0.044
loghr <- log(all$hr)
W <- se_ref^2
B <- if (nrow(all) > 1) var(loghr) else 0
V <- W + B + B/nrow(all)
hr_pool <- exp(mean(loghr))
z <- qnorm(0.975)
se_pool <- sqrt(V)
dfm <- if (nrow(all) > 1 && B > 0)
  (nrow(all) - 1) * (1 + W / (B * (1 + 1/nrow(all))))^2 else Inf
zt <- if (is.finite(dfm)) qt(0.975, dfm) else z
lci <- exp(mean(loghr) - zt*se_pool)
uci <- exp(mean(loghr) + zt*se_pool)
p <- 2*(1 - pnorm(abs(mean(loghr)/se_pool)))

cat(sprintf("  SEQTaRget 路线 A: HR = %.4f (%.4f-%.4f), p = %.4f\n",
            hr_pool, lci, uci, p))
cat(sprintf("  stay-level IPTW Cox: HR = %.4f (%.4f-%.4f)\n",
            mean(iptw$hr), exp(mean(log(iptw$hr))-zt*se_pool), exp(mean(log(iptw$hr))+zt*se_pool)))
cat(sprintf("  20 imp HR range: %.4f - %.4f, mean = %.4f\n",
            min(all$hr), max(all$hr), mean(all$hr)))
cat(sprintf("  差异: %.2f%%\n", 100*abs(hr_pool - mean(iptw$hr))/mean(iptw$hr)))

# 也自己算一个 W=0 (只看 20 imp 间变异性的 pooled, 不用外借 SE)
se_from_B <- sqrt(B + B/nrow(all))
lci_B <- exp(mean(loghr) - zt*se_from_B)
uci_B <- exp(mean(loghr) + zt*se_from_B)
cat(sprintf("\n  纯 between 方差 SE (W=0): HR = %.4f (%.4f-%.4f)\n",
            hr_pool, lci_B, uci_B))

# 写 pooled
pool_out <- data.table(
  analysis = "90-day death · SEQTaRget ITT g-formula · 8-chunk expand",
  version = "V28.4-seqA",
  date = format(Sys.Date()),
  m = nrow(all),
  hr = hr_pool, lci = lci, uci = uci, p = p,
  loghr_mean = mean(loghr),
  W = W, B = B, V = V,
  se_pool = se_pool,
  se_ref_source = "stay-level IPTW Cox mean SE ≈ 0.044"
)
pool_csv <- file.path(OUT_DIR, "v28.4_death90_seq20_pooled.csv")
fwrite(pool_out, pool_csv)
cat(sprintf("\n已保存: %s + %s\n", OUT_CSV, pool_csv))
