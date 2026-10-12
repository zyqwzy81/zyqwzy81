# ============================================================================
# death90_chunk.R
# 多次入组框架下 90 天全因死亡 g-formula HR —— 低内存"分块展开→合并→一次分析"
#
# 方法严格性（28 天实证，见会话记录）：
#   * SEQTaRget ITT 未加权路径的 internal.analysis 仅使用展开表 params@DT；
#   * SEQexpand 是纯数据变换，按 stay 分块展开后 rbind 与全量展开逐值相同；
#   * 合并表只拟合一次 outcome 模型并做 hazard，HR 与官方一次跑逐位相等。
# 因此在不触发全量展开爆内存的前提下，得到与全量完全一致的点估计。
#
# 开关：--imps / --chunksize / --tag
# 版本 V28.4  日期 2026-10-11
# ============================================================================
suppressMessages(source("/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/multi_common.R"))
ns <- asNamespace("SEQTaRget")
internal.analysis <- ns$internal.analysis
internal.hazard   <- ns$internal.hazard
init_formula_cache <- ns$init_formula_cache

argv <- commandArgs(trailingOnly = TRUE)
gv <- function(k, d) { i <- match(k, argv); if (is.na(i)) d else argv[i + 1] }
IMPS <- as.integer(strsplit(gv("--imps", "1"), ",")[[1]])
CHUNK_SIZE <- as.integer(gv("--chunksize", "1000"))
TAG  <- gv("--tag", "d90chunk")
H90  <- 90L

OUT_T <- file.path(V28_BASE, "TTE/v28.4/output/tables/death90_multi")
dir.create(OUT_T, recursive = TRUE, showWarnings = FALSE)

# ---- dod（日期，与正式脚本一致）----
dod_tsv <- tempfile(fileext = ".tsv")
qry <- "SELECT subject_id, dod::date FROM mimiciv_hosp.patients WHERE dod IS NOT NULL"
pgpw <- Sys.getenv("MIMIC_PG_PASSWORD", unset = "")
envpre <- if (nzchar(pgpw)) paste0("PGPASSWORD=", shQuote(pgpw), " ") else ""
system(paste(envpre, "psql -h localhost -p 5432 -U postgres -d mimiciv",
             "-At -F '|' -c", shQuote(qry, type = "sh"),
             ">", shQuote(dod_tsv, type = "sh")))
dod <- unique(fread(dod_tsv, sep = "|", col.names = c("subject_id", "dod")))
dod[, subject_id := as.integer(subject_id)]
AUG_F <- file.path(V28_BASE, "TTE/v26/output/sequential_trials_augmented.parquet")
subj_map <- unique(as.data.table(read_parquet(AUG_F))[, .(stay_id, subject_id)])

off <- build_event_offsets()

# 在块 daily 上写 oc_death90（复刻 v28.4_death90_multi.R base_orig 口径）
add_oc90 <- function(daily, tc) {
  o <- unique(tc[, .(stay_id, time_zero, origin_day)])
  o <- merge(o, subj_map, by = "stay_id", all.x = TRUE)
  o <- merge(o, dod, by = "subject_id", all.x = TRUE)
  o[, d90 := as.integer(as.Date(dod) - as.Date(time_zero))]
  o[, abs90 := fifelse(is.finite(d90) & d90 >= 0 & d90 <= H90,
                       origin_day + d90, NA_integer_)]
  om <- merge(unique(tc[, .(stay_id, time_zero)]),
              o[, .(stay_id, time_zero, abs90)],
              by = c("stay_id", "time_zero"), all.x = TRUE)
  ev <- om[is.finite(abs90), .(abs_day = unique(abs90)), by = stay_id]
  daily[, oc_death90 := 0L]
  daily[ev, on = c("stay_id", day = "abs_day"), oc_death90 := 1L]
  list(daily = daily, ev = ev)
}

expand_chunk <- function(tc) {
  daily <- build_multi_daily(tc, off, H = H90)
  z <- add_oc90(daily, tc); daily <- z$daily
  opts <- SEQopts(hazard = TRUE, bootstrap = FALSE,
                  followup.max = H90, seed = SEED,
                  expand.only = TRUE,
                  selection.first_trial = FALSE)
  DT <- SEQuential(daily, id.col = "stay_id", time.col = "day",
                   eligible.col = "eligible", treatment.col = "tx_init",
                   outcome.col = "oc_death90", time_varying.cols = DAG14,
                   fixed.cols = list(), method = "ITT",
                   options = opts, verbose = FALSE)
  list(DT = DT, ev = z$ev)
}

run_imp <- function(imp_id) {
  cat(sprintf("\n================ imp %d ================\n", imp_id))
  trials <- load_trials(imp_id)
  all_stays <- sort(unique(trials$stay_id))
  n_stay <- length(all_stays)
  n_chunk <- ceiling(n_stay / CHUNK_SIZE)
  cat(sprintf("stays=%d chunks=%d (chunk size=%d)\n",
              n_stay, n_chunk, CHUNK_SIZE))

  # ---- 分块展开并累积 ----
  chunks <- vector("list", n_chunk)
  ev_list <- vector("list", n_chunk)
  for (ci in seq_len(n_chunk)) {
    sv <- all_stays[(ci - 1) * CHUNK_SIZE + seq_len(min(CHUNK_SIZE,
                     n_stay - (ci - 1) * CHUNK_SIZE))]
    tc <- trials[stay_id %in% sv]
    t0 <- Sys.time()
    z <- expand_chunk(tc)
    chunks[[ci]] <- z$DT
    ev_list[[ci]] <- z$ev
    cat(sprintf("  chunk %2d/%d rows=%8d time=%5.1fs\n",
                ci, n_chunk, nrow(z$DT),
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    rm(tc, z); invisible(gc())
  }

  full_DT <- rbindlist(chunks)
  rm(chunks); invisible(gc())
  cat(sprintf("merged expanded table: rows=%d cols=%d\n",
              nrow(full_DT), ncol(full_DT)))

  # ---- 事件数校验 ----
  ev_all <- rbindlist(ev_list)
  unique_dead <- uniqueN(ev_all$stay_id)
  cat(sprintf("[check] 90d unique dead stays=%d (expect ~4799)\n", unique_dead))

  # ---- 用一小块跑完整 SEQuential 捕获 params 模板（列规格与全量一致）----
  template_stays <- all_stays[seq_len(min(200L, n_stay))]
  tc_t <- trials[stay_id %in% template_stays]
  daily_t <- build_multi_daily(tc_t, off, H = H90)
  zt <- add_oc90(daily_t, tc_t); daily_t <- zt$daily
  rm(CAP_PARAMS, envir = globalenv())
  trace(internal.analysis,
        tracer = quote(assign("CAP_PARAMS", params, envir = globalenv())),
        print = FALSE, where = ns)
  opts_t <- SEQopts(hazard = TRUE, bootstrap = FALSE,
                    followup.max = H90, seed = SEED,
                    expand.only = FALSE,
                    selection.first_trial = FALSE)
  invisible(SEQuential(daily_t, id.col = "stay_id", time.col = "day",
                       eligible.col = "eligible", treatment.col = "tx_init",
                       outcome.col = "oc_death90", time_varying.cols = DAG14,
                       fixed.cols = list(), method = "ITT",
                       options = opts_t, verbose = FALSE))
  suppressMessages(untrace(internal.analysis, where = ns))
  params <- get("CAP_PARAMS", envir = globalenv())
  rm(CAP_PARAMS, envir = globalenv())

  # ---- 注入全量合并表，一次分析 ----
  p2 <- params
  p2@DT <- full_DT
  t0 <- Sys.time()
  an <- internal.analysis(p2)
  cat(sprintf("internal.analysis: time=%.1fs\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  model_list <- list(an[[1]]$model[[1]])
  cache <- init_formula_cache(p2)
  h <- internal.hazard(model_list, p2, cache)

  data.table(imp = imp_id, hr = unname(h["Hazard ratio"]),
             n_stay = n_stay, n_expanded = nrow(full_DT),
             unique_dead90 = unique_dead)
}

res <- rbindlist(lapply(IMPS, run_imp))
out_f <- file.path(OUT_T, sprintf("death90_%s.csv", TAG))
fwrite(res, out_f)

cat("\n================ 结果 ================\n")
print(res)
cat(sprintf("\n已保存: %s\n", out_f))
