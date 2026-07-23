# ============================================================================
# 01_expand_only.R
# SEQTaRget 分批引擎 —— 步骤一：仅扩展，不落盘建模
# 输出：按 trial_period 分区的 Parquet 文件 + 扩展摘要 JSON
# ============================================================================

# ---- 加载配置 ----
source("config.R")

# ---- 加载依赖 ----
suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(jsonlite)
})

# ---- 峰值内存监控工具 ----
get_peak_mem_mb <- function() {
  # Linux 下读取 /proc/self/status 的 VmPeak
  if (file.exists("/proc/self/status")) {
    lines <- readLines("/proc/self/status", warn = FALSE)
    peak_line <- grep("^VmPeak:", lines, value = TRUE)
    if (length(peak_line) > 0) {
      kb <- as.numeric(gsub("[^0-9]", "", peak_line))
      return(round(kb / 1024, 2))
    }
  }
  NA_real_
}

module_log_section("步骤一：CCW 扩展")
module_log_info(sprintf("读取输入: %s", INPUT_PARQUET))

t_start <- Sys.time()

# ---- 读取原始数据 ----
dt <- as.data.table(arrow::read_parquet(INPUT_PARQUET))
module_log_info(sprintf("原始数据: n=%d 行, %d 列, %d 例患者",
                        nrow(dt), ncol(dt), length(unique(dt[[ID_COL]]))))

# ---- 基本校验 ----
req_cols <- c(ID_COL, PERIOD_COL, TREATMENT_SRC_COL, OUTCOME_COL, SURV_TIME_COL)
missing_cols <- setdiff(req_cols, names(dt))
if (length(missing_cols) > 0) {
  stop(sprintf("输入数据缺少必要列: %s", paste(missing_cols, collapse = ", ")))
}

# ---- 生成扩展数据 ----
# 每条原始 patient-period 记录转换为一个 clone
# clone_id: 唯一克隆标识
# tx_init: 该患者在 trial entry 时的治疗状态（ITT）
# outcome: 28 天死亡结局
# followup_time: 从该期 time_zero 到事件/删失的时间（天）
# eligible: 该期是否合格（默认原始数据已筛选为合格期）
set.seed(RANDOM_SEED)
setorder(dt, stay_id, trial_num)

expanded <- dt[, .(
  stay_id = get(ID_COL),
  subject_id = if ("subject_id" %in% names(dt)) subject_id else NA_integer_,
  trial_period = get(PERIOD_COL),
  time_zero = if ("time_zero" %in% names(dt)) time_zero else NA_character_,
  baseline_phosphate = if ("baseline_phosphate" %in% names(dt)) baseline_phosphate else NA_real_,
  tx_init = as.integer(get(TREATMENT_SRC_COL) == TREATMENT_LABEL),
  outcome = as.integer(get(OUTCOME_COL)),
  followup_time = as.numeric(get(SURV_TIME_COL)),
  eligible = 1L,
  age, gender, baseline_sofa, egfr, sepsis3,
  invasive_vent_at_t0, vasoactive_at_t0, charlson, admission_type
)]

# ---- clone_id 唯一编号 ----
expanded[, clone_id := .I]
setcolorder(expanded, c("clone_id", "stay_id", "trial_period", "tx_init",
                        "outcome", "followup_time", "eligible"))

# ---- 合格性修正（可选） ----
# 若 followup_time <= 0 且未发生结局，则视为不合格（避免零生存时间）
expanded[followup_time <= 0 & outcome == 0, eligible := 0L]

# ---- 按 trial_period 分区保存 ----
periods <- sort(unique(expanded$trial_period))
module_log_info(sprintf("共 %d 个 trial_period: %s",
                        length(periods), paste(periods, collapse = ", ")))

for (p in periods) {
  sub <- expanded[trial_period == p]
  out_path <- file.path(EXPANDED_DIR, sprintf("trial_period=%d", p))
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  arrow::write_parquet(sub, file.path(out_path, "part-0.parquet"))
  module_log_info(sprintf("  period=%d: n_clones=%d, events=%d",
                          p, nrow(sub), sum(sub$outcome, na.rm = TRUE)))
  rm(sub); gc(verbose = FALSE)
}

# ---- 扩展摘要 ----
summary_dt <- expanded[, .(
  n_clones = .N,
  n_patients = length(unique(stay_id)),
  n_events = sum(outcome, na.rm = TRUE),
  event_rate = round(mean(outcome, na.rm = TRUE), 4),
  mean_followup = round(mean(followup_time, na.rm = TRUE), 2),
  mean_baseline_phosphate = round(mean(baseline_phosphate, na.rm = TRUE), 3),
  treatment_rate = round(mean(tx_init, na.rm = TRUE), 4)
), by = trial_period]

overall <- data.table(
  trial_period = 999L,
  n_clones = nrow(expanded),
  n_patients = length(unique(expanded$stay_id)),
  n_events = sum(expanded$outcome, na.rm = TRUE),
  event_rate = round(mean(expanded$outcome, na.rm = TRUE), 4),
  mean_followup = round(mean(expanded$followup_time, na.rm = TRUE), 2),
  mean_baseline_phosphate = round(mean(expanded$baseline_phosphate, na.rm = TRUE), 3),
  treatment_rate = round(mean(expanded$tx_init, na.rm = TRUE), 4)
)
summary_dt <- rbind(summary_dt, overall, use.names = TRUE, fill = TRUE)

summary_csv <- file.path(RESULTS_DIR, "expansion_summary.csv")
fwrite(summary_dt, summary_csv)

summary_json <- list(
  n_patients = length(unique(expanded$stay_id)),
  n_clones = nrow(expanded),
  n_periods = length(periods),
  total_events = sum(expanded$outcome, na.rm = TRUE),
  overall_event_rate = round(mean(expanded$outcome, na.rm = TRUE), 4),
  treatment_rate = round(mean(expanded$tx_init, na.rm = TRUE), 4),
  periods = as.list(periods),
  peak_mem_mb = get_peak_mem_mb(),
  runtime_sec = round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 2),
  output_dir = EXPANDED_DIR,
  summary_csv = summary_csv
)
json_path <- file.path(RESULTS_DIR, "expansion_summary.json")
writeLines(jsonlite::toJSON(summary_json, pretty = TRUE, auto_unbox = TRUE), json_path)

module_log_section("扩展步骤完成")
module_log_info(sprintf("总克隆数: %d, 总患者数: %d, 总事件数: %d",
                        nrow(expanded), length(unique(expanded$stay_id)),
                        sum(expanded$outcome, na.rm = TRUE)))
module_log_info(sprintf("峰值内存: %s MB, 运行时间: %.2f 秒",
                        ifelse(is.na(get_peak_mem_mb()), "NA", as.character(get_peak_mem_mb())),
                        as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
module_log_info(sprintf("扩展摘要 CSV: %s", summary_csv))
module_log_info(sprintf("扩展摘要 JSON: %s", json_path))

# ---- 清理 ----
rm(expanded, dt, summary_dt, overall); gc(verbose = FALSE)
