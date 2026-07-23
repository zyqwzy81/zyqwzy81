# ============================================================================
# V26.9 次要结局时间分辨率分析
#
# 审稿要求 M3: 次要结局时间分辨率
#   对每个次要结局计算效应量，有生存时间信息的做 Cox 分段 HR，
#   仅有事件状态的做 logistic 回归。
#
# 数据: sequential_trials_imputed_v26.parquet
# 输出: TTE/v26/output/v26_9/
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(survival)
  library(ggplot2)
})

# ============================================================================
# 1. 参数配置
# ============================================================================

OUT_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9"
TBL_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")
LOG_DIR <- file.path(OUT_DIR, "logs")
for (d in c(TBL_DIR, FIG_DIR, LOG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

DATA_FILE <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_imputed_v26.parquet"

# DAG14 协变量集
DAG14 <- c(
  "age", "gender", "baseline_sofa", "sepsis3", "charlson",
  "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
  "admission_type", "baseline_phosphate", "lactate_t0", "refeeding_risk",
  "corrected_calcium_mmol_min", "nutrition_support_type"
)

# 次要结局定义
# type: "cox" = Cox 分段HR（有生存时间）, "logistic" = logistic 回归（仅事件状态）
SECONDARY_OUTCOMES <- list(
  death_28day = list(
    var = "death_28day_num", type = "cox", label = "28天死亡",
    time = "surv_time", events_expected = 5610
  ),
  aki = list(
    var = "aki_kdigo2_plus", type = "logistic", label = "AKI KDIGO≥2期",
    events_expected = 5209
  ),
  cv_comprehensive = list(
    var = "cardiovascular_event_comprehensive", type = "logistic",
    label = "CV复合结局（综合）", events_expected = 19238
  ),
  cv_2comp = list(
    var = "cardiovascular_event_2comp", type = "logistic",
    label = "CV复合结局（2成分）", events_expected = 6360
  ),
  cv_3comp = list(
    var = "cardiovascular_event_3comp", type = "logistic",
    label = "CV复合结局（3成分）", events_expected = 7692
  ),
  arrhythmia = list(
    var = "new_arrhythmia_rhythm", type = "logistic",
    label = "新发心律失常", events_expected = 3967
  ),
  myocardial = list(
    var = "myocardial_injury", type = "logistic",
    label = "心肌损伤", events_expected = 1979
  ),
  transfusion = list(
    var = "transfusion_at_t0", type = "logistic",
    label = "输血", events_expected = 3115
  )
)

# 分段时间窗口（天）
TIME_CUTS <- c(0, 3, 7, 14, 28)
SEGMENTS <- list(
  "0-3d"   = c(0, 3),
  "3-7d"   = c(3, 7),
  "7-14d"  = c(7, 14),
  "14-28d" = c(14, 28)
)

SEED <- 20260629L
set.seed(SEED)

# ============================================================================
# 2. 日志系统
# ============================================================================
log_file <- file.path(LOG_DIR, sprintf("v26_9_secondary_time_%s.log",
                                        format(Sys.time(), "%Y%m%d_%H%M%S")))
log_con  <- file(log_file, open = "wt")

log_info <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg),
      file = log_con, append = TRUE)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
  flush(stdout())
}

log_section <- function(title) {
  cat(sprintf("\n========== %s ==========\n", title),
      file = log_con, append = TRUE)
  cat(sprintf("\n========== %s ==========\n", title))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

log_info("V26.9 次要结局时间分辨率分析启动")
log_info(sprintf("种子: %d, 输出: %s", SEED, OUT_DIR))

# ============================================================================
# 3. 数据加载
# ============================================================================
log_section("数据加载")

dt_full <- as.data.table(read_parquet(DATA_FILE))
dt_full[, `:=`(
  tx_int = as.integer(itt_group == "treatment"),
  event  = death_28day_num
)]
setorder(dt_full, stay_id, trial_num)

log_info(sprintf("Parquet 数据: %d 行, %d 事件(28天死亡), %d 列",
                 nrow(dt_full), sum(dt_full$event), ncol(dt_full)))

# DAG14 可用性检查
dag14_avail <- intersect(DAG14, names(dt_full))
log_info(sprintf("DAG14 可用: %d/%d", length(dag14_avail), length(DAG14)))

# 字符型转 factor
for (v in dag14_avail) {
  if (is.character(dt_full[[v]])) {
    dt_full[[v]] <- as.factor(dt_full[[v]])
  }
}

# ============================================================================
# 4. 次要结局分析
# ============================================================================
log_section("次要结局分析")

results_list <- list()
cox_segment_list <- list()  # death_28day 分段 HR

for (oname in names(SECONDARY_OUTCOMES)) {
  oc <- SECONDARY_OUTCOMES[[oname]]
  outcome_var <- oc$var
  outcome_lbl <- oc$label
  oc_type     <- oc$type

  log_info(sprintf("---- %s (%s) ----", outcome_lbl, outcome_var))

  # ---- 4a. 变量存在性 ----
  if (!outcome_var %in% names(dt_full)) {
    log_info(sprintf("  [跳过] 变量不存在: %s", outcome_var))
    next
  }

  # ---- 4b. DAG14 协变量准备 ----
  dag14_active <- setdiff(dag14_avail, outcome_var)
  dag14_active <- dag14_active[sapply(dag14_active, function(v) {
    x <- dt_full[[v]]
    if (is.factor(x))   return(nlevels(x) > 1)
    if (is.character(x)) return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]

  dag14_str <- paste(dag14_active, collapse = " + ")

  # ====================================================================
  # Part A: Cox 分段 HR（有生存时间的结局）
  # ====================================================================
  if (oc_type == "cox") {
    time_var <- oc$time
    if (!time_var %in% names(dt_full)) {
      log_info(sprintf("  [跳过] 生存时间变量不存在: %s", time_var))
      next
    }

    # 整体 Cox
    form_global <- sprintf("Surv(%s, %s) ~ tx_int + %s", time_var, outcome_var, dag14_str)
    fit_global <- tryCatch(
      coxph(as.formula(form_global), data = dt_full, cluster = stay_id),
      error = function(e) NULL
    )

    if (!is.null(fit_global)) {
      s <- summary(fit_global)
      hr <- exp(coef(fit_global)["tx_int"])
      ci <- exp(confint(fit_global)["tx_int", ])
      p  <- s$coefficients["tx_int", "Pr(>|z|)"]
      n_event <- sum(dt_full[[outcome_var]] == 1, na.rm = TRUE)

      log_info(sprintf("  [整体] HR=%.4f (%.4f–%.4f), p=%.6f, events=%d",
                       hr, ci[1], ci[2], p, n_event))

      results_list[[paste0(oname, "_overall")]] <- data.table(
        outcome     = outcome_lbl,
        outcome_var = outcome_var,
        method      = "Cox 整体",
        estimate    = hr,
        ci_lower    = ci[1],
        ci_upper    = ci[2],
        p_value     = p,
        n_events    = n_event
      )
    }

    # 分段 Cox（各时间窗口）
    for (seg_name in names(SEGMENTS)) {
      t_start <- SEGMENTS[[seg_name]][1]
      t_end   <- SEGMENTS[[seg_name]][2]

      dt_seg <- dt_full[get(time_var) > t_start & get(time_var) <= t_end]
      dt_seg[, seg_event := fifelse(get(outcome_var) == 1 & get(time_var) <= t_end, 1L, 0L)]

      if (sum(dt_seg$seg_event) < 10) {
        log_info(sprintf("  [分段 %s] 事件数不足 (%d), 跳过", seg_name, sum(dt_seg$seg_event)))
        next
      }

      form_seg <- sprintf("Surv(%s, seg_event) ~ tx_int + %s", time_var, dag14_str)
      fit_seg <- tryCatch(
        coxph(as.formula(form_seg), data = dt_seg, cluster = stay_id),
        error = function(e) NULL
      )

      if (!is.null(fit_seg)) {
        s  <- summary(fit_seg)
        hr <- exp(coef(fit_seg)["tx_int"])
        ci <- exp(confint(fit_seg)["tx_int", ])
        p  <- s$coefficients["tx_int", "Pr(>|z|)"]

        log_info(sprintf("  [分段 %s] HR=%.4f (%.4f–%.4f), p=%.6f, events=%d, n=%d",
                         seg_name, hr, ci[1], ci[2], p,
                         sum(dt_seg$seg_event), nrow(dt_seg)))

        cox_segment_list[[paste0(oname, "_", seg_name)]] <- data.table(
          outcome     = outcome_lbl,
          outcome_var = outcome_var,
          segment     = seg_name,
          t_start     = t_start,
          t_end       = t_end,
          estimate    = hr,
          ci_lower    = ci[1],
          ci_upper    = ci[2],
          p_value     = p,
          n_events    = sum(dt_seg$seg_event),
          n_total     = nrow(dt_seg)
        )
      }
    }
  }

  # ====================================================================
  # Part B: Logistic 回归（仅有事件状态的结局）
  # ====================================================================
  if (oc_type == "logistic") {
    n_event <- sum(dt_full[[outcome_var]] == 1, na.rm = TRUE)
    n_total <- sum(!is.na(dt_full[[outcome_var]]))

    if (n_event < 10) {
      log_info(sprintf("  [跳过] 事件数不足: events=%d", n_event))
      next
    }

    form_log <- sprintf("%s ~ tx_int + %s", outcome_var, dag14_str)
    fit_log <- tryCatch(
      glm(as.formula(form_log), data = dt_full, family = binomial),
      error = function(e) NULL
    )

    if (!is.null(fit_log)) {
      s <- summary(fit_log)
      if ("tx_int" %in% rownames(s$coefficients)) {
        or   <- exp(coef(fit_log)["tx_int"])
        ci   <- tryCatch(exp(confint(fit_log, "tx_int")), error = function(e) c(NA, NA))
        p    <- s$coefficients["tx_int", "Pr(>|z|)"]

        log_info(sprintf("  OR=%.4f (%.4f–%.4f), p=%.6f, events=%d/%d",
                         or, ci[1], ci[2], p, n_event, n_total))

        results_list[[oname]] <- data.table(
          outcome     = outcome_lbl,
          outcome_var = outcome_var,
          method      = "Logistic",
          estimate    = or,
          ci_lower    = ci[1],
          ci_upper    = ci[2],
          p_value     = p,
          n_events    = n_event
        )
      }
    }
  }
}

# ============================================================================
# 5. 结果汇总输出
# ============================================================================
log_section("结果汇总")

# ---- 5a. 整体结果（Cox 整体 + Logistic） ----
dt_results <- rbindlist(results_list, fill = TRUE)
setorder(dt_results, p_value)

log_info("次要结局整体结果:")
cat(sprintf("\n%-30s %-10s %10s %10s %10s %10s %8s\n",
            "结局", "方法", "估计值", "CI下限", "CI上限", "p值", "事件数"),
    file = log_con, append = TRUE)
cat(sprintf("%-30s %-10s %10s %10s %10s %10s %8s\n",
            "结局", "方法", "估计值", "CI下限", "CI上限", "p值", "事件数"))
for (i in seq_len(nrow(dt_results))) {
  r <- dt_results[i]
  line <- sprintf("  %-28s %-10s %10.4f %10.4f %10.4f %10.6f %8d",
                  r$outcome, r$method, r$estimate, r$ci_lower, r$ci_upper,
                  r$p_value, r$n_events)
  cat(line, "\n", file = log_con, append = TRUE)
  cat(line, "\n")
}

fwrite(dt_results, file.path(TBL_DIR, "v26_9_secondary_time_results.csv"))
log_info(sprintf("结果已保存: %s",
                 file.path(TBL_DIR, "v26_9_secondary_time_results.csv")))

# ---- 5b. Cox 分段结果 ----
if (length(cox_segment_list) > 0) {
  dt_segments <- rbindlist(cox_segment_list)
  fwrite(dt_segments, file.path(TBL_DIR, "v26_9_secondary_time_segments.csv"))
  log_info(sprintf("分段HR已保存: %s",
                   file.path(TBL_DIR, "v26_9_secondary_time_segments.csv")))

  cat(sprintf("\n%-25s %10s %10s %10s %10s %8s %8s\n",
              "分段", "HR", "CI下限", "CI上限", "p值", "事件", "总数"),
      file = log_con, append = TRUE)
  cat(sprintf("%-25s %10s %10s %10s %10s %8s %8s\n",
              "分段", "HR", "CI下限", "CI上限", "p值", "事件", "总数"))
  for (i in seq_len(nrow(dt_segments))) {
    r <- dt_segments[i]
    line <- sprintf("  %-23s %10.4f %10.4f %10.4f %10.6f %8d %8d",
                    sprintf("%s %s", r$outcome, r$segment),
                    r$estimate, r$ci_lower, r$ci_upper, r$p_value,
                    r$n_events, r$n_total)
    cat(line, "\n", file = log_con, append = TRUE)
    cat(line, "\n")
  }
}

# ============================================================================
# 6. 森林图（整体结果）
# ============================================================================
log_section("森林图生成")

if (nrow(dt_results) >= 2) {
  dt_results[, `:=`(
    y_pos = seq_len(.N),
    sig   = fifelse(p_value < 0.05, "p<0.05", "NS")
  )]

  p_forest <- ggplot(dt_results, aes(x = estimate, y = y_pos)) +
    geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = sig),
                   height = 0.2, linewidth = 0.8) +
    geom_point(aes(color = sig, fill = sig), shape = 22, size = 3.5, stroke = 0.8) +
    scale_y_continuous(breaks = dt_results$y_pos,
                       labels = dt_results$outcome) +
    scale_color_manual(values = c("p<0.05" = "#BC3C29", "NS" = "#0072B5"),
                       name = NULL) +
    scale_fill_manual(values = c("p<0.05" = "#BC3C29", "NS" = "#0072B5"),
                      name = NULL) +
    labs(
      title    = "次要结局分析 — DAG14 校正",
      subtitle = sprintf("Logistic 回归 OR / Cox HR · n=%d · %d 项结局",
                         nrow(dt_full), nrow(dt_results)),
      x        = "OR / HR (95% CI)",
      y        = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0, size = 13),
      plot.subtitle    = element_text(color = "grey40", hjust = 0, size = 9),
      legend.position  = "bottom",
      axis.text.y      = element_text(size = 10),
      plot.margin      = margin(12, 16, 10, 10)
    )

  ggsave(file.path(FIG_DIR, "v26_9_secondary_forest.pdf"),
         plot = p_forest, width = 8, height = max(4, nrow(dt_results) * 0.5),
         device = "pdf")
  ggsave(file.path(FIG_DIR, "v26_9_secondary_forest.png"),
         plot = p_forest, width = 8, height = max(4, nrow(dt_results) * 0.5),
         dpi = 300, device = "png")
  log_info(sprintf("森林图已保存: %s",
                   file.path(FIG_DIR, "v26_9_secondary_forest.pdf")))

  # ---- 6b. 分段 HR 森林图 ----
  if (length(cox_segment_list) > 0 && nrow(dt_segments) >= 2) {
    dt_segments[, `:=`(
      y_pos = seq_len(.N),
      sig   = fifelse(p_value < 0.05, "p<0.05", "NS"),
      seg_label = sprintf("%s\n%s", outcome, segment)
    )]

    p_seg <- ggplot(dt_segments, aes(x = estimate, y = reorder(seg_label, y_pos))) +
      geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
      geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = sig),
                     height = 0.2, linewidth = 0.8) +
      geom_point(aes(color = sig, fill = sig), shape = 22, size = 3.5, stroke = 0.8) +
      scale_color_manual(values = c("p<0.05" = "#BC3C29", "NS" = "#0072B5"),
                         name = NULL) +
      scale_fill_manual(values = c("p<0.05" = "#BC3C29", "NS" = "#0072B5"),
                        name = NULL) +
      labs(
        title    = "死亡分段HR分析 — 各时间窗口",
        subtitle = sprintf("Cox 分段回归 · DAG14 校正 · 时间窗口: 0-3/3-7/7-14/14-28天"),
        x        = "HR (95% CI)",
        y        = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title      = element_text(face = "bold", hjust = 0, size = 13),
        plot.subtitle   = element_text(color = "grey40", hjust = 0, size = 9),
        legend.position = "bottom",
        axis.text.y     = element_text(size = 9),
        plot.margin     = margin(12, 16, 10, 10)
      )

    ggsave(file.path(FIG_DIR, "v26_9_death_segmented_forest.pdf"),
           plot = p_seg, width = 8, height = max(4, nrow(dt_segments) * 0.45),
           device = "pdf")
    ggsave(file.path(FIG_DIR, "v26_9_death_segmented_forest.png"),
           plot = p_seg, width = 8, height = max(4, nrow(dt_segments) * 0.45),
           dpi = 300, device = "png")
    log_info(sprintf("分段HR森林图已保存"))
  }
} else {
  log_info("森林图跳过: 不足 2 个成功分析")
}

# ============================================================================
# 7. 保存完整 RDS
# ============================================================================
log_section("结果保存")

rds_out <- list(
  version          = "V26.9",
  script           = "v26_9_secondary_time_analysis.R",
  timestamp        = as.character(Sys.time()),
  seed             = SEED,
  dag14            = DAG14,
  dag14_available  = dag14_avail,
  results_overall  = dt_results,
  results_segments = if (length(cox_segment_list) > 0) dt_segments else data.table()
)
saveRDS(rds_out, file.path(TBL_DIR, "v26_9_secondary_time_analysis.rds"))
log_info(sprintf("RDS 已保存: %s",
                 file.path(TBL_DIR, "v26_9_secondary_time_analysis.rds")))

# ============================================================================
# 8. 流程自审核
# ============================================================================
log_section("流程自审核")

audit <- function(check, pass, detail = "") {
  mark <- if (pass) "PASS" else "FAIL"
  log_info(sprintf("[%s] %s: %s", mark, check, detail))
  return(pass)
}

all_pass <- TRUE
all_pass <- all_pass && audit("数据加载",
  nrow(dt_full) > 10000, sprintf("n=%d", nrow(dt_full)))
all_pass <- all_pass && audit("DAG14 可用性",
  length(dag14_avail) >= 12, sprintf("%d/14", length(dag14_avail)))
all_pass <- all_pass && audit("次要结局分析",
  nrow(dt_results) >= 3, sprintf("%d 结局成功", nrow(dt_results)))
all_pass <- all_pass && audit("森林图",
  file.exists(file.path(FIG_DIR, "v26_9_secondary_forest.pdf")), "")
all_pass <- all_pass && audit("结果 CSV",
  file.exists(file.path(TBL_DIR, "v26_9_secondary_time_results.csv")), "")
all_pass <- all_pass && audit("RDS 输出",
  file.exists(file.path(TBL_DIR, "v26_9_secondary_time_analysis.rds")), "")

log_info(sprintf("\n审核结果: %s", if (all_pass) "全部通过 ✓" else "存在未通过项 ✗"))

# ============================================================================
# 9. 收尾
# ============================================================================
close(log_con)

cat("\n")
cat("################################################################################\n")
cat("      V26.9 次要结局时间分辨率分析完成\n")
cat("################################################################################\n\n")

cat(sprintf("成功分析: %d 项结局\n", nrow(dt_results)))
cat(sprintf("  Cox 整体: %d (28天死亡)\n",
            sum(dt_results$method == "Cox 整体")))
cat(sprintf("  Logistic: %d (AKI/CV/心律失常/心肌损伤/输血)\n",
            sum(dt_results$method == "Logistic")))

if (length(cox_segment_list) > 0) {
  cat(sprintf("  分段 HR: %d 个时间窗口\n", nrow(dt_segments)))
}
cat(sprintf("\n日志: %s\n", log_file))
cat(sprintf("结果: %s\n", file.path(TBL_DIR, "v26_9_secondary_time_results.csv")))
cat(sprintf("森林图: %s\n", file.path(FIG_DIR, "v26_9_secondary_forest.pdf")))
cat("========== V26.9 次要结局时间分析完成 ==========\n")
