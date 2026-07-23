# ============================================================================
# 05_verify_outputs.R
# SEQTaRget 分批引擎 —— 步骤五：输出文件格式与字段完整性验证
# 输出：verification_report.json
# ============================================================================

source("config.R")

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(jsonlite)
})

module_log_section("步骤五：验证输出文件格式与字段完整性")

# ---- 验证结果容器 ----
# format_checks: 文件格式与字段完整性
# substantive_checks: 统计/方法学一致性（不影响格式验证 overall_pass）
format_checks <- list()
substantive_checks <- list()

add_format_check <- function(name, status, details = NULL) {
  format_checks[[name]] <<- list(status = status, details = details)
  badge <- ifelse(status, "PASS", "FAIL")
  module_log_info(sprintf("[%s] %s", badge, name))
  if (!is.null(details)) module_log_info(sprintf("  -> %s", details))
}

add_substantive_check <- function(name, status, details = NULL) {
  substantive_checks[[name]] <<- list(status = status, details = details)
  badge <- ifelse(status, "PASS", "WARN")
  module_log_info(sprintf("[%s] %s", badge, name))
  if (!is.null(details)) module_log_info(sprintf("  -> %s", details))
}

# ---- 1. 扩展后 Parquet 批次 ----
module_log_info("检查扩展后 Parquet 批次 ...")
input_df <- arrow::read_parquet(INPUT_PARQUET)
expected_periods <- sort(unique(input_df[[PERIOD_COL]]))
rm(input_df); gc(verbose = FALSE)

expanded_periods <- integer(0)
expanded_cols_ok <- TRUE
expanded_required_cols <- c(
  "clone_id", ID_COL, "trial_period", "tx_init", "outcome",
  "followup_time", "eligible", PS_COVARIATES
)

for (p in expected_periods) {
  pfile <- file.path(EXPANDED_DIR, sprintf("trial_period=%d", p), "part-0.parquet")
  if (!file.exists(pfile)) next
  expanded_periods <- c(expanded_periods, p)
  df <- arrow::read_parquet(pfile)
  missing <- setdiff(expanded_required_cols, names(df))
  if (length(missing) > 0) {
    expanded_cols_ok <- FALSE
    add_format_check(
      sprintf("expanded_period_%d_cols", p), FALSE,
      paste("缺少列:", paste(missing, collapse = ", "))
    )
  }
}
add_format_check(
  "expanded_parquet_exists",
  length(expanded_periods) > 0,
  sprintf("发现 %d 个 period 的扩展 Parquet", length(expanded_periods))
)
add_format_check("expanded_required_columns", expanded_cols_ok)

# ---- 2. 加权后 Parquet 批次 ----
module_log_info("检查加权后 Parquet 批次 ...")
weighted_periods <- integer(0)
weighted_cols_ok <- TRUE
weighted_required_cols <- c(
  "clone_id", ID_COL, "trial_period", "tx_init", "outcome",
  "followup_time", "eligible", "iptw", "iptw_trunc"
)

for (p in expanded_periods) {
  pfile <- file.path(WEIGHTED_DIR, sprintf("trial_period=%d", p), "part-0.parquet")
  if (!file.exists(pfile)) next
  weighted_periods <- c(weighted_periods, p)
  df <- arrow::read_parquet(pfile)
  missing <- setdiff(weighted_required_cols, names(df))
  if (length(missing) > 0) {
    weighted_cols_ok <- FALSE
    add_format_check(
      sprintf("weighted_period_%d_cols", p), FALSE,
      paste("缺少列:", paste(missing, collapse = ", "))
    )
  }
}
add_format_check(
  "weighted_parquet_exists",
  length(weighted_periods) > 0,
  sprintf("发现 %d 个 period 的加权 Parquet", length(weighted_periods))
)
add_format_check("weighted_required_columns", weighted_cols_ok)

# ---- 3. PS/IPTW 摘要 CSV ----
module_log_info("检查 PS/IPTW 摘要 CSV ...")
ps_csv <- file.path(RESULTS_DIR, "batch_ps_iptw_summary.csv")
ps_ok <- file.exists(ps_csv)
ps_cols_ok <- FALSE
ps_sensible <- FALSE
if (ps_ok) {
  ps_dt <- fread(ps_csv)
  required <- c(
    "trial_period", "n_clones", "n_patients", "n_events",
    "treatment_rate", "ps_mean", "ps_sd", "iptw_mean", "iptw_sd",
    "iptw_trunc_mean", "iptw_trunc_sd", "ess", "ps_model_converged"
  )
  ps_cols_ok <- all(required %in% names(ps_dt))
  ps_sensible <- all(ps_dt$n_clones > 0) &&
    all(ps_dt$ps_model_converged %in% c(TRUE, FALSE))
}
add_format_check("ps_iptw_summary_exists", ps_ok)
add_format_check("ps_iptw_summary_columns", ps_cols_ok)
add_format_check("ps_iptw_summary_sensible", ps_sensible)

# ---- 4. Cox HR 摘要 CSV ----
module_log_info("检查 Cox HR 摘要 CSV ...")
hr_summary_csv <- file.path(RESULTS_DIR, "seqtarget_batch_hr_summary.csv")
hr_summary_ok <- file.exists(hr_summary_csv)
hr_summary_cols_ok <- FALSE
hr_summary_sensible <- FALSE
if (hr_summary_ok) {
  hr_dt <- fread(hr_summary_csv)
  required <- c(
    "method", "n_periods", "log_hr", "se", "hr",
    "hr_lower", "hr_upper", "note"
  )
  hr_summary_cols_ok <- all(required %in% names(hr_dt))
  hr_summary_sensible <- all(is.finite(hr_dt$hr)) && all(hr_dt$hr > 0)
}
add_format_check("hr_summary_exists", hr_summary_ok)
add_format_check("hr_summary_columns", hr_summary_cols_ok)
add_format_check("hr_summary_sensible", hr_summary_sensible)

# ---- 5. 每期 HR CSV ----
module_log_info("检查每期 HR CSV ...")
hr_period_csv <- file.path(RESULTS_DIR, "seqtarget_batch_hr_period_specific.csv")
hr_period_ok <- file.exists(hr_period_csv)
hr_period_cols_ok <- FALSE
hr_period_sensible <- FALSE
if (hr_period_ok) {
  hp_dt <- fread(hr_period_csv)
  required <- c(
    "trial_period", "n_clones", "n_patients", "n_events",
    "log_hr", "se", "hr", "hr_lower", "hr_upper", "p_value"
  )
  hr_period_cols_ok <- all(required %in% names(hp_dt))
  hr_period_sensible <- all(is.finite(hp_dt$hr)) && all(hp_dt$hr > 0)
}
add_format_check("hr_period_specific_exists", hr_period_ok)
add_format_check("hr_period_specific_columns", hr_period_cols_ok)
add_format_check("hr_period_specific_sensible", hr_period_sensible)

# ---- 6. SEQTaRget 对比 CSV ----
module_log_info("检查 SEQTaRget 对比 CSV ...")
cmp_csv <- file.path(RESULTS_DIR, "seqtarget_vs_batch_engine_comparison.csv")
cmp_ok <- file.exists(cmp_csv)
cmp_cols_ok <- FALSE
cmp_batch_has_ci <- FALSE
cmp_direction_consistent <- FALSE
seqtarget_na_note <- NULL

if (cmp_ok) {
  cmp_dt <- fread(cmp_csv)
  required <- c(
    "method", "tag", "hr", "hr_lower", "hr_upper", "se",
    "mean_iptw", "ps_mean", "runtime_sec", "peak_mem_mb", "n_periods"
  )
  cmp_cols_ok <- all(required %in% names(cmp_dt))

  batch_rows <- cmp_dt[method == "batch_engine"]
  cmp_batch_has_ci <- all(
    is.finite(batch_rows$hr_lower) &
      is.finite(batch_rows$hr_upper) &
      is.finite(batch_rows$se)
  )

  seq_rows <- cmp_dt[method == "SEQTaRget"]
  seq_has_ci <- all(is.finite(seq_rows$hr_lower) & is.finite(seq_rows$hr_upper))
  if (!seq_has_ci) {
    seqtarget_na_note <- paste(
      "SEQTaRget 行 hr_lower/hr_upper/se 为 NA，",
      "符合无 bootstrap 回退的预期"
    )
    module_log_warn(seqtarget_na_note)
  }

  if (nrow(seq_rows) > 0 && all(is.finite(seq_rows$hr)) &&
      nrow(batch_rows) > 0 && all(is.finite(batch_rows$hr))) {
    cmp_direction_consistent <- all((seq_rows$hr > 1) == (batch_rows$hr > 1))
  }
}
add_format_check("comparison_csv_exists", cmp_ok)
add_format_check("comparison_csv_columns", cmp_cols_ok)
add_format_check("comparison_batch_has_ci", cmp_batch_has_ci)
add_substantive_check(
  "comparison_direction_consistent",
  cmp_direction_consistent,
  ifelse(cmp_direction_consistent, "两方法点估计方向一致", "方向不一致")
)

# ---- 7. 扩展摘要 JSON ----
module_log_info("检查扩展摘要 JSON ...")
exp_json <- file.path(RESULTS_DIR, "expansion_summary.json")
exp_json_ok <- file.exists(exp_json)
exp_json_sensible <- FALSE
if (exp_json_ok) {
  ej <- jsonlite::fromJSON(exp_json)
  required_keys <- c("n_patients", "n_clones", "n_periods", "total_events")
  exp_json_sensible <- all(required_keys %in% names(ej)) &&
    ej$n_patients > 0 && ej$n_clones >= ej$n_patients
}
add_format_check("expansion_summary_json_exists", exp_json_ok)
add_format_check("expansion_summary_json_sensible", exp_json_sensible)

# ---- 8. 可视化输出 ----
module_log_info("检查可视化输出 ...")
fig_pdf <- file.path(RESULTS_DIR, "FigS_batch_vs_seqtarget_hr.pdf")
fig_png <- file.path(RESULTS_DIR, "FigS_batch_vs_seqtarget_hr.png")
add_format_check("comparison_figure_pdf", file.exists(fig_pdf))
add_format_check("comparison_figure_png", file.exists(fig_png))

# ---- 保存验证报告 ----
format_pass <- all(vapply(format_checks, function(x) x$status, logical(1)))
report <- list(
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  output_dir = OUTPUT_DIR,
  format_checks = format_checks,
  substantive_checks = substantive_checks,
  seqtarget_na_note = seqtarget_na_note,
  format_checks_pass = format_pass,
  overall_pass = format_pass
)
report_json <- file.path(RESULTS_DIR, "verification_report.json")
jsonlite::write_json(report, report_json, auto_unbox = TRUE, pretty = TRUE)
module_log_info(sprintf("验证报告已保存: %s", report_json))

status_text <- ifelse(report$overall_pass, "格式验证全部通过", "格式验证存在失败项")
module_log_section(sprintf("验证完成: %s", status_text))
