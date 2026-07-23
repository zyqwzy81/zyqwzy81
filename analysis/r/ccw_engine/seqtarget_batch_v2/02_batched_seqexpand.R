# ============================================================================
# 02_batched_seqexpand.R — Phase 2: Batched Daily SEQexpand
#
# 核心改进: 将 patient-period 记录展开为日粒度（与 SEQTaRget 的 SEQexpand 对齐）
#   - 每条 patient-period → 从 day=0 到 floor(surv_time) 的日记录
#   - 死亡日 outcome=1，其余日 outcome=0
#   - 应用 Phase 1 计算的 pooled PS + IPTW 权重
#   - 添加 followup, followup², trial, trial² 列（为 pooled logistic 做准备）
#
# 批次策略: 按 trial_period 分批展开，每批次独立保存为 parquet 分区
# 验证: 每批次保存展开摘要 + 与 Phase 1 权重一致性检查
#
# 关键参数（来自 config_v2.R）:
#   SEQEXPAND_MAX_DAY = 28   — 最大展开天数
#   PERIOD_FILTER_MIN_EVENTS = 5  — 最小事件数
# ============================================================================

source("config_v2.R")

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(jsonlite)
})

module_log_section("Phase 2: Batched SEQexpand（日粒度展开）")

t_start_total <- Sys.time()

# ---- 1. 发现 Phase 1 权重分区 ----
period_dirs <- list.dirs(WEIGHT_DIR, recursive = FALSE)
periods <- sort(as.integer(gsub("trial_period=", "", basename(period_dirs))))
if (length(periods) == 0) {
  stop(sprintf("未在 %s 发现权重分区，请先运行 01_pooled_ps.R", WEIGHT_DIR))
}
module_log_info(sprintf("发现 %d 个权重分区: %s", length(periods),
                         paste(periods, collapse = ", ")))

# ---- 2. SEQexpand 函数 ----
# 实现与 SEQTaRget 内部 SEQexpand 一致的日粒度展开逻辑:
#
# 输入: 一条 patient-period 记录（含 surv_time, outcome, tx_init, iptw_trunc, covariates）
# 输出: floor(surv_time) + 1 条日记录:
#   - day: 0, 1, 2, ..., floor(surv_time)
#   - y_daily: 死亡日=1, 其余=0
#   - 各协变量继承展开前的值
#   - iptw_trunc 权重继承展开前的值（pre-expansion weights）
seqexpand_record <- function(dt) {
  # dt: data.table with columns [stay_id, trial_period, tx_init, outcome,
  #                               surv_time, iptw_trunc, <covariates>]
  # 计算展开天数
  dt[, day_max := as.integer(floor(pmax(surv_time, 0)))]
  dt[day_max > SEQEXPAND_MAX_DAY, day_max := SEQEXPAND_MAX_DAY]
  dt[day_max < 0, day_max := 0L]
  
  # 对每条记录按 (day_max + 1) 行展开
  expanded <- dt[rep(seq_len(.N), day_max + 1L)]
  
  # 添加 day 列和日常结局
  expanded[, day := seq_len(.N) - 1L, by = .(stay_id, trial_period)]
  expanded[, y_daily := 0L]
  expanded[day == day_max & outcome == 1L, y_daily := 1L]
  
  # 添加 followup / trial 时间列（为 pooled logistic 做准备）
  expanded[, followup    := as.numeric(day)]
  expanded[, followup_sq := as.numeric(day^2)]
  expanded[, trial       := as.numeric(trial_period)]
  expanded[, trial_sq    := as.numeric(trial_period^2)]
  
  # 移除中间列
  expanded[, c("day_max", "outcome", "surv_time") := NULL]
  
  # 列顺序: key → weights → temporal → covariates
  setcolorder(expanded, intersect(
    c("stay_id", "trial_period", "day", "tx_init", "y_daily",
      "followup", "followup_sq", "trial", "trial_sq",
      "iptw_trunc", "ps", "iptw"),
    names(expanded)
  ))
  
  expanded[]
}

# ---- 3. 逐期展开 ----
expansion_summaries <- list()
expand_counts <- list()

for (p in periods) {
  t_start <- Sys.time()
  
  # 3a: 读取 Phase 1 权重数据
  in_path <- file.path(WEIGHT_DIR, sprintf("trial_period=%d", p), "part-0.parquet")
  if (!file.exists(in_path)) {
    module_log_warn(sprintf("period=%d: 权重文件不存在，跳过", p))
    next
  }
  
  batch_weights <- as.data.table(arrow::read_parquet(in_path))
  
  # 3b: 基本校验
  if (nrow(batch_weights) < PERIOD_FILTER_MIN_N) {
    module_log_warn(sprintf("period=%d: 样本数 %d < 最小阈值 %d，跳过",
                            p, nrow(batch_weights), PERIOD_FILTER_MIN_N))
    rm(batch_weights); gc(verbose = FALSE)
    next
  }
  
  if (!"iptw_trunc" %in% names(batch_weights)) {
    module_log_warn(sprintf("period=%d: 缺少 iptw_trunc 列，跳过", p))
    rm(batch_weights); gc(verbose = FALSE)
    next
  }
  
  # 3c: 执行日粒度展开
  expanded <- seqexpand_record(batch_weights)
  
  # 3d: 保存展开后的分区
  out_path <- file.path(EXPAND_DIR, sprintf("trial_period=%d", p))
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  arrow::write_parquet(expanded, file.path(out_path, "part-0.parquet"))
  
  # 3e: 批次摘要
  n_orig <- nrow(batch_weights)
  n_expanded <- nrow(expanded)
  n_events_daily <- sum(expanded$y_daily, na.rm = TRUE)
  
  # 展开倍数 = 展开行 / 原始行
  expansion_ratio <- round(n_expanded / n_orig, 2)
  
  # 每日事件分布（检查展开后事件是否集中在死亡日）
  event_by_day <- expanded[y_daily == 1, .N, by = day][order(day)]
  
  summ <- data.table(
    trial_period = p,
    n_original = n_orig,
    n_expanded = n_expanded,
    expansion_ratio = expansion_ratio,
    n_events_original = sum(batch_weights$outcome, na.rm = TRUE),
    n_events_expanded = n_events_daily,
    mean_followup_days = round(mean(batch_weights$surv_time, na.rm = TRUE), 2),
    iptw_mean = round(mean(expanded$iptw_trunc, na.rm = TRUE), 4),
    iptw_sd = round(sd(expanded$iptw_trunc, na.rm = TRUE), 4),
    runtime_sec = round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 2),
    peak_mem_mb = get_peak_mem_mb()
  )
  expansion_summaries[[as.character(p)]] <- summ
  
  # 事件按日分布摘要
  if (nrow(event_by_day) > 0) {
    event_summary <- data.table(
      trial_period = p,
      n_event_days = nrow(event_by_day),
      first_event_day = min(event_by_day$day),
      last_event_day = max(event_by_day$day),
      peak_event_day = event_by_day[which.max(N), day],
      peak_event_count = max(event_by_day$N)
    )
  } else {
    event_summary <- data.table(
      trial_period = p, n_event_days = 0L,
      first_event_day = NA_integer_, last_event_day = NA_integer_,
      peak_event_day = NA_integer_, peak_event_count = 0L
    )
  }
  expand_counts[[as.character(p)]] <- event_summary
  
  module_log_info(sprintf(
    "  period=%d: %d×%.1f→%d行, 事件=%d, 展开倍率=%.1f, 权重均值=%.4f, 时间=%.2fs",
    p, n_orig, summ$mean_followup_days, n_expanded, n_events_daily,
    expansion_ratio, summ$iptw_mean, summ$runtime_sec
  ))
  
  rm(batch_weights, expanded, event_by_day); gc(verbose = FALSE)
}

# ---- 4. 保存展开摘要 ----
expansion_summary_dt <- rbindlist(expansion_summaries, use.names = TRUE, fill = TRUE)
expansion_summary_csv <- file.path(RESULT_DIR, "seqexpand_summary.csv")
fwrite(expansion_summary_dt, expansion_summary_csv)

expand_counts_dt <- rbindlist(expand_counts, use.names = TRUE, fill = TRUE)
expand_counts_csv <- file.path(RESULT_DIR, "seqexpand_event_by_day.csv")
fwrite(expand_counts_dt, expand_counts_csv)

# ---- 5. 全局摘要 ----
total_original <- sum(expansion_summary_dt$n_original)
total_expanded <- sum(expansion_summary_dt$n_expanded)
total_events <- sum(expansion_summary_dt$n_events_expanded)
overall_ratio <- round(total_expanded / total_original, 2)

module_log_section("SEQexpand 汇总")
module_log_info(sprintf("全局: %d 原始记录 → %d 日记录 (展开倍率 %.1f)",
                        total_original, total_expanded, overall_ratio))
module_log_info(sprintf("总事件数: %d (展开后)", total_events))
module_log_info(sprintf("有效 period 数: %d", nrow(expansion_summary_dt)))

# ---- 6. 保存元数据 ----
metadata <- list(
  phase = "02_batched_seqexpand",
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  weight_dir = WEIGHT_DIR,
  expand_dir = EXPAND_DIR,
  seqexpand_max_day = SEQEXPAND_MAX_DAY,
  n_periods_total = length(periods),
  n_periods_expanded = nrow(expansion_summary_dt),
  total_records_original = total_original,
  total_records_expanded = total_expanded,
  overall_expansion_ratio = overall_ratio,
  total_events_expanded = total_events,
  runtime_sec = round(as.numeric(difftime(Sys.time(), t_start_total, units = "secs")), 2),
  peak_mem_mb = get_peak_mem_mb(),
  expansion_summary = expansion_summary_csv,
  event_distribution = expand_counts_csv
)
metadata_json <- file.path(RESULT_DIR, "02_seqexpand_metadata.json")
writeLines(toJSON(metadata, pretty = TRUE, auto_unbox = TRUE), metadata_json)

module_log_info(sprintf("总运行时间: %.2f 秒", metadata$runtime_sec))
module_log_info(sprintf("峰值内存: %s MB", ifelse(is.na(metadata$peak_mem_mb), "NA", as.character(metadata$peak_mem_mb))))

# 清理
rm(expansion_summaries, expand_counts); gc(verbose = FALSE)
