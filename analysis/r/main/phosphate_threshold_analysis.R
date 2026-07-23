# ============================================================================
# V28 磷酸保护性阈值分析
#
# 目的：基于 Phosphate < 0.50 mmol/L 亚组的显著保护性效应（HR=0.803, p=0.029），
#       进一步分析保护性效应的阈值边界，识别最佳治疗受益人群。
#
# 分析内容：
#   1. 多切点亚组分析：在 0.32/0.40/0.50/0.60/0.65/0.70 mmol/L 切点分层
#   2. 限制性立方样条（RCS）交互分析：连续磷酸盐 × 治疗交互
#   3. 阈值识别：寻找 HR 越过 1.0 的磷酸盐水平
#
# 方法：Cox + DAG14 + cluster-robust SE (cluster=stay_id) + m=20 Rubin 合并
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(survival)
  library(rms)
})

# ============================================================================
# 1. 参数配置
# ============================================================================
OUT_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v28/output/tables"
LOG_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v28/output/logs"
FIG_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v28/output/figures"
for (d in c(OUT_DIR, LOG_DIR, FIG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

DATA_FILE         <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/sequential_trials_imputed_v26.parquet"
MICE_PARQUET_TEMPL <- "/data/小米云盘/文档/论文/血磷/sql/result/mice_imputed/v26_mice_imputed_imp%02d.parquet"

# DAG14 协变量集（分析时排除 baseline_phosphate，因其作为效应修饰因子）
DAG14_BASE <- c(
  "age", "gender", "baseline_sofa", "sepsis3", "charlson",
  "egfr", "invasive_vent_at_t0", "vasoactive_at_t0",
  "admission_type", "lactate_t0", "refeeding_risk",
  "corrected_calcium_mmol_min", "nutrition_support_type"
)

# 多切点定义（基于临床意义的磷酸盐分层）
THRESHOLDS <- c(0.32, 0.40, 0.50, 0.60, 0.65, 0.70, 0.80)

N_IMPUTE <- 20L
SEED     <- 20260629L
set.seed(SEED)

# ============================================================================
# 2. 日志系统
# ============================================================================
log_file <- file.path(LOG_DIR, sprintf("v28_phosphate_threshold_%s.log",
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

log_info("V28 磷酸保护性阈值分析启动")
log_info(sprintf("切点数: %d, MICE 插补集数: %d", length(THRESHOLDS), N_IMPUTE))

# ============================================================================
# 3. 数据加载
# ============================================================================
log_section("数据加载")

dt_full <- as.data.table(read_parquet(DATA_FILE))
dt_full[, `:=`(
  tx_int   = as.integer(itt_group == "treatment"),
  event    = death_28day_num,
  id_int   = as.integer(factor(stay_id))
)]
setorder(dt_full, stay_id, trial_num)

log_info(sprintf("主数据: %d 行, %d 事件, %d 列",
                 nrow(dt_full), sum(dt_full$event), ncol(dt_full)))

# 磷酸盐分布
log_info("磷酸盐分布:")
log_info(sprintf("  均值: %.4f, 中位数: %.4f, SD: %.4f",
                 mean(dt_full$baseline_phosphate, na.rm=TRUE),
                 median(dt_full$baseline_phosphate, na.rm=TRUE),
                 sd(dt_full$baseline_phosphate, na.rm=TRUE)))
log_info(sprintf("  范围: %.4f - %.4f",
                 min(dt_full$baseline_phosphate, na.rm=TRUE),
                 max(dt_full$baseline_phosphate, na.rm=TRUE)))
for (q in c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)) {
  log_info(sprintf("  P%02d: %.4f", as.integer(q*100),
                   quantile(dt_full$baseline_phosphate, q, na.rm=TRUE)))
}

# 各切点样本量
log_info("各切点样本量分布:")
for (th in THRESHOLDS) {
  n_below <- sum(dt_full$baseline_phosphate < th, na.rm=TRUE)
  n_above <- sum(dt_full$baseline_phosphate >= th, na.rm=TRUE)
  ev_below <- sum(dt_full$death_28day_num[dt_full$baseline_phosphate < th] == 1, na.rm=TRUE)
  log_info(sprintf("  < %.2f: n=%d (%.1f%%), events=%d | >= %.2f: n=%d (%.1f%%)",
                   th, n_below, 100*n_below/nrow(dt_full), ev_below,
                   th, n_above, 100*n_above/nrow(dt_full)))
}

# key_cols
key_cols <- c("stay_id", "subject_id", "trial_num", "time_zero",
              "itt_group", "tx_int", "gender", "admission_type",
              "surv_time", "death_time_from_t0", "dischtime_from_t0",
              "race_group", "icu_type",
              "death_28day_num", "nutrition_support_type", "baseline_phosphate")
key_cols <- intersect(key_cols, names(dt_full))
dt_key <- dt_full[, ..key_cols]
rm(dt_full); gc(TRUE)

# ============================================================================
# 4. 辅助函数
# ============================================================================

# 单切点 Cox 分析
run_threshold_cox <- function(dt_sub, dag_covs) {
  if (sum(dt_sub$tx_int == 1) < 10 || sum(dt_sub$tx_int == 0) < 10) return(NULL)
  n_events <- sum(dt_sub$death_28day_num == 1, na.rm = TRUE)
  if (n_events < 20) return(NULL)

  # 字符型变量转因子
  for (v in dag_covs) {
    if (v %in% names(dt_sub) && is.character(dt_sub[[v]])) {
      dt_sub[, (v) := as.factor(get(v))]
    }
  }

  # 可用协变量
  covs <- intersect(dag_covs, names(dt_sub))
  covs <- covs[sapply(covs, function(v) {
    x <- dt_sub[[v]]
    if (is.factor(x))    return(nlevels(x) > 1)
    if (is.character(x)) return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]

  if (length(covs) < 2) return(NULL)

  formula_str <- sprintf("Surv(surv_time, death_28day_num) ~ tx_int + %s",
                         paste(covs, collapse = " + "))
  fit <- tryCatch(
    coxph(as.formula(formula_str), data = dt_sub, cluster = stay_id),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  if (!"tx_int" %in% names(coef(fit))) return(NULL)

  s <- summary(fit)
  beta <- coef(fit)["tx_int"]
  robust_se <- tryCatch(
    s$coefficients["tx_int", "robust se"],
    error = function(e) s$coefficients["tx_int", "se(coef)"]
  )

  return(list(
    beta     = beta,
    se       = robust_se,
    hr       = exp(beta),
    ci_lower = exp(beta - 1.96 * robust_se),
    ci_upper = exp(beta + 1.96 * robust_se),
    p_value  = 2 * pnorm(abs(beta / robust_se), lower.tail = FALSE),
    n        = nrow(dt_sub),
    n_events = n_events
  ))
}

# Rubin 规则合并
rubin_pool <- function(imp_results) {
  valid <- Filter(function(x) !is.null(x), imp_results)
  if (length(valid) < 3) return(NULL)

  betas <- sapply(valid, function(x) x$beta)
  ses   <- sapply(valid, function(x) x$se)
  m <- length(valid)

  # Rubin 规则
  q_bar <- mean(betas)
  u_bar <- mean(ses^2)
  b_var <- var(betas)
  total_var <- u_bar + (1 + 1/m) * b_var
  se_pooled <- sqrt(total_var)

  # df
  df_old <- (m - 1) * (1 + u_bar / ((1 + 1/m) * b_var))^2
  df <- df_old * (1 + m * u_var_adj(b_var, u_bar, m))
  df <- ifelse(is.na(df) | df <= 0, m - 1, df)

  return(list(
    hr       = exp(q_bar),
    ci_lower = exp(q_bar - 1.96 * se_pooled),
    ci_upper = exp(q_bar + 1.96 * se_pooled),
    p_value  = 2 * pt(abs(q_bar / se_pooled), df = df, lower.tail = FALSE),
    beta     = q_bar,
    se       = se_pooled,
    W        = u_bar,
    B        = b_var,
    df       = df
  ))
}

u_var_adj <- function(b_var, u_bar, m) {
  if (b_var < 1e-10) return(0)
  return(u_bar / ((1 + 1/m) * b_var))
}

# ============================================================================
# 5. Part A: 多切点亚组分析
# ============================================================================
log_section("Part A: 多切点亚组分析")

threshold_results <- list()

for (th in THRESHOLDS) {
  log_info(sprintf("--- 切点: Phosphate < %.2f mmol/L ---", th))

  imp_results <- list()

  for (imp in 1:N_IMPUTE) {
    mice_file <- sprintf(MICE_PARQUET_TEMPL, imp)
    if (!file.exists(mice_file)) {
      log_info(sprintf("  [imp%02d] 文件不存在: %s", imp, mice_file))
      next
    }

    dt_imp <- as.data.table(read_parquet(mice_file))
    imp_covs <- setdiff(intersect(DAG14_BASE, names(dt_imp)), "baseline_phosphate")

    # 合并 key + imp
    dt_comb <- cbind(dt_key, dt_imp[, ..imp_covs])

    # 筛选子组
    dt_sub <- dt_comb[baseline_phosphate < th]

    if (nrow(dt_sub) < 100) {
      log_info(sprintf("  [imp%02d] 样本不足: %d", imp, nrow(dt_sub)))
      next
    }

    res <- run_threshold_cox(dt_sub, DAG14_BASE)
    if (!is.null(res)) {
      imp_results[[length(imp_results) + 1]] <- res
    }

    rm(dt_imp, dt_comb, dt_sub); gc(FALSE)
  }

  pooled <- rubin_pool(imp_results)
  if (is.null(pooled)) {
    log_info(sprintf("  切点 %.2f: Rubin 合并失败", th))
    next
  }

  n_vec <- sapply(imp_results, function(x) x$n)
  ev_vec <- sapply(imp_results, function(x) x$n_events)

  threshold_results[[as.character(th)]] <- list(
    threshold = th,
    hr        = pooled$hr,
    ci_lower  = pooled$ci_lower,
    ci_upper  = pooled$ci_upper,
    p_value   = pooled$p_value,
    n         = as.integer(median(n_vec)),
    n_events  = as.integer(median(ev_vec)),
    n_imp     = length(imp_results)
  )

  sig <- ifelse(pooled$p_value < 0.05, "*", "")
  log_info(sprintf("  Phosphate < %.2f: HR=%.4f (%.4f-%.4f), p=%.4f %s, n=%d, events=%d",
                   th, pooled$hr, pooled$ci_lower, pooled$ci_upper,
                   pooled$p_value, sig,
                   as.integer(median(n_vec)), as.integer(median(ev_vec))))
}

# 输出多切点结果
dt_threshold <- rbindlist(lapply(threshold_results, function(x) {
  data.table(
    threshold = sprintf("< %.2f", x$threshold),
    hr        = round(x$hr, 4),
    ci_lower  = round(x$ci_lower, 4),
    ci_upper  = round(x$ci_upper, 4),
    p_value   = signif(x$p_value, 4),
    n         = x$n,
    n_events  = x$n_events,
    n_imp     = x$n_imp,
    significant = ifelse(x$p_value < 0.05, "p<0.05", "NS")
  )
}))

threshold_csv <- file.path(OUT_DIR, "v28_phosphate_threshold_results.csv")
fwrite(dt_threshold, threshold_csv)
log_info(sprintf("多切点结果已保存: %s", threshold_csv))

# ============================================================================
# 6. Part B: RCS 交互分析（连续磷酸盐效应修饰）
# ============================================================================
log_section("Part B: RCS 交互分析")

# 使用第一个插补集进行 RCS 分析（探索性）
dt_imp1 <- as.data.table(read_parquet(sprintf(MICE_PARQUET_TEMPL, 1)))
imp_covs <- setdiff(intersect(DAG14_BASE, names(dt_imp1)), "baseline_phosphate")
dt_rcs <- cbind(dt_key, dt_imp1[, ..imp_covs])

# 字符变量转因子
for (v in DAG14_BASE) {
  if (v %in% names(dt_rcs) && is.character(dt_rcs[[v]])) {
    dt_rcs[, (v) := as.factor(get(v))]
  }
}

# 排除缺失磷酸盐
dt_rcs <- dt_rcs[!is.na(baseline_phosphate)]
log_info(sprintf("RCS 分析数据: %d 行", nrow(dt_rcs)))

# RCS 节点（4节点，位于 5%/35%/65%/95% 分位数）
phos_knots <- quantile(dt_rcs$baseline_phosphate, c(0.05, 0.35, 0.65, 0.95), na.rm=TRUE)
log_info(sprintf("RCS 节点: %.4f, %.4f, %.4f, %.4f",
                 phos_knots[1], phos_knots[2], phos_knots[3], phos_knots[4]))

# 创建 RCS 变量
ddist <- datadist(dt_rcs)
options(datadist = 'ddist')

# 拟合 RCS 交互模型
# 模型：Surv(surv_time, death) ~ tx_int * rcs(baseline_phosphate, 4) + DAG13
rcs_formula <- as.formula(sprintf(
  "Surv(surv_time, death_28day_num) ~ tx_int * rcs(baseline_phosphate, 4) + %s",
  paste(intersect(DAG14_BASE, names(dt_rcs)), collapse = " + ")
))

fit_rcs <- tryCatch({
  cph(rcs_formula, data = dt_rcs, x = TRUE, y = TRUE, surv = TRUE)
}, error = function(e) {
  log_info(sprintf("RCS cph 拟合失败: %s", e$message))
  # 回退到标准 coxph
  coxph(rcs_formula, data = dt_rcs, cluster = stay_id)
})

# 交互项 p 值
if (inherits(fit_rcs, "cph")) {
  anova_result <- anova(fit_rcs)
  # 查找 tx_int × baseline_phosphate 交互项
  log_info("RCS ANOVA:")
  print(anova_result)
} else {
  s <- summary(fit_rcs)
  # 提取交互项
  int_terms <- grep("tx_int:baseline_phosphate", names(coef(fit_rcs)), value = TRUE)
  log_info(sprintf("交互项: %s", paste(int_terms, collapse = ", ")))
  for (term in int_terms) {
    beta <- coef(fit_rcs)[term]
    se <- summary(fit_rcs)$coefficients[term, "robust se"]
    p <- 2 * pnorm(abs(beta/se), lower.tail = FALSE)
    log_info(sprintf("  %s: beta=%.4f, se=%.4f, p=%.4f", term, beta, se, p))
  }
}

# ============================================================================
# 7. Part C: 治疗效应随磷酸盐变化曲线（基于离散切点结果）
# ============================================================================
log_section("Part C: 治疗效应随磷酸盐变化曲线")

# 直接基于 Part A 的多切点结果绘制阶梯图，避免 Predict 兼容性问题
if (length(threshold_results) > 0) {
  # 构造 HR 网格（基于离散切点的阶梯曲线）
  dt_hr_grid <- rbindlist(lapply(threshold_results, function(x) {
    data.table(
      threshold_cut = x$threshold,
      hr        = x$hr,
      hr_lower  = x$ci_lower,
      hr_upper  = x$ci_upper,
      p_value   = x$p_value,
      n         = x$n,
      n_events  = x$n_events
    )
  }))
  setorder(dt_hr_grid, threshold_cut)

  # 识别 HR 越过 1.0 的磷酸盐水平（线性插值）
  if (nrow(dt_hr_grid) >= 2) {
    cross_points <- c()
    for (i in 1:(nrow(dt_hr_grid) - 1)) {
      y1 <- dt_hr_grid$hr[i]; y2 <- dt_hr_grid$hr[i + 1]
      x1 <- dt_hr_grid$threshold_cut[i]; x2 <- dt_hr_grid$threshold_cut[i + 1]
      if ((y1 - 1) * (y2 - 1) < 0) {
        cross_x <- x1 + (1 - y1) / (y2 - y1) * (x2 - x1)
        cross_points <- c(cross_points, cross_x)
        log_info(sprintf("  HR 越过 1.0 的磷酸盐水平: %.4f mmol/L (线性插值)", cross_x))
      }
    }
    if (length(cross_points) == 0) {
      if (all(dt_hr_grid$hr < 1)) {
        log_info("  所有切点 HR 均 < 1.0（保护性趋势）")
      } else if (all(dt_hr_grid$hr > 1)) {
        log_info("  所有切点 HR 均 > 1.0（有害趋势）")
      }
    }
  }

  # 保存 HR 网格
  hr_grid_csv <- file.path(OUT_DIR, "v28_phosphate_threshold_hr_grid.csv")
  fwrite(dt_hr_grid, hr_grid_csv)
  log_info(sprintf("HR 网格已保存: %s", hr_grid_csv))

  # 绘制效应曲线（阶梯图 + 误差棒）
  png(file.path(FIG_DIR, "FigS_Phosphate_Threshold_HR_Curve.png"),
      width = 10, height = 7, units = "in", res = 300)
  par(mar = c(5.5, 5, 4, 2), family = "sans")

  # X 轴范围：0.32 到 0.85
  xlim <- c(0.32, 0.85)
  ylim <- c(0.65, 1.20)

  plot(NA, xlim = xlim, ylim = ylim,
       xlab = "Baseline Phosphate Threshold (mmol/L)",
       ylab = "Hazard Ratio (IV vs No IV, 28-day mortality)",
       main = "Treatment Effect by Baseline Phosphate Threshold",
       xaxt = "n", las = 1, bty = "l")
  axis(1, at = c(0.32, 0.40, 0.50, 0.60, 0.65, 0.70, 0.80))

  # 灰色背景带（保护性区域 HR<1）
  rect(xlim[1], 0.65, xlim[2], 1.0, col = adjustcolor("#2E8B57", 0.06), border = NA)
  # 灰色背景带（有害区域 HR>1）
  rect(xlim[1], 1.0, xlim[2], 1.20, col = adjustcolor("#BC3C29", 0.06), border = NA)

  # 阶梯连线
  lines(dt_hr_grid$threshold_cut, dt_hr_grid$hr, type = "S",
        lwd = 2, col = "#0072B5")

  # 误差棒
  arrows(dt_hr_grid$threshold_cut, dt_hr_grid$hr_lower,
         dt_hr_grid$threshold_cut, dt_hr_grid$hr_upper,
         angle = 90, code = 3, length = 0.04, col = "#0072B5", lwd = 1.2)

  # 点
  points(dt_hr_grid$threshold_cut, dt_hr_grid$hr, pch = 21,
         bg = ifelse(dt_hr_grid$p_value < 0.05, "#BC3C29", "#0072B5"),
         col = "#0072B5", cex = 1.5, lwd = 1.5)

  # 显著性标注
  for (i in 1:nrow(dt_hr_grid)) {
    if (dt_hr_grid$p_value[i] < 0.05) {
      text(dt_hr_grid$threshold_cut[i], dt_hr_grid$hr_upper[i] + 0.03,
           sprintf("p=%.3f", dt_hr_grid$p_value[i]),
           cex = 0.75, col = "#BC3C29", font = 2)
    }
  }

  # 参考线
  abline(h = 1, col = "#E74C3C", lty = 2, lwd = 1.5)

  # 文字标注：左侧"Protective"，右侧"Null effect"
  text(0.38, 0.72, "Protective", col = "#2E8B57", cex = 0.9, font = 2)
  text(0.78, 1.10, "Null effect", col = "#666666", cex = 0.9, font = 2)

  # 图例
  legend("topright",
         legend = c("HR (95% CI)", "p < 0.05", "Null (HR=1)"),
         pch = c(21, 21, NA),
         pt.bg = c("#0072B5", "#BC3C29", NA),
         lty = c(NA, NA, 2),
         lwd = c(NA, NA, 1.5),
         col = c("#0072B5", "#BC3C29", "#E74C3C"),
         bty = "n", cex = 0.9)

  # 副标题（关键统计量）
  mtext("SEQTaRget CCW + DAG14 + cluster-robust SE + m=20 Rubin's rules",
        side = 3, line = 0.3, cex = 0.7, col = "#666666", adj = 0)

  dev.off()
  log_info("效应曲线图已保存")

  # 同时输出 PDF 矢量版本
  pdf(file.path(FIG_DIR, "FigS_Phosphate_Threshold_HR_Curve.pdf"),
      width = 10, height = 7)
  par(mar = c(5.5, 5, 4, 2), family = "sans")
  plot(NA, xlim = xlim, ylim = ylim,
       xlab = "Baseline Phosphate Threshold (mmol/L)",
       ylab = "Hazard Ratio (IV vs No IV, 28-day mortality)",
       main = "Treatment Effect by Baseline Phosphate Threshold",
       xaxt = "n", las = 1, bty = "l")
  axis(1, at = c(0.32, 0.40, 0.50, 0.60, 0.65, 0.70, 0.80))
  rect(xlim[1], 0.65, xlim[2], 1.0, col = adjustcolor("#2E8B57", 0.06), border = NA)
  rect(xlim[1], 1.0, xlim[2], 1.20, col = adjustcolor("#BC3C29", 0.06), border = NA)
  lines(dt_hr_grid$threshold_cut, dt_hr_grid$hr, type = "S", lwd = 2, col = "#0072B5")
  arrows(dt_hr_grid$threshold_cut, dt_hr_grid$hr_lower,
         dt_hr_grid$threshold_cut, dt_hr_grid$hr_upper,
         angle = 90, code = 3, length = 0.04, col = "#0072B5", lwd = 1.2)
  points(dt_hr_grid$threshold_cut, dt_hr_grid$hr, pch = 21,
         bg = ifelse(dt_hr_grid$p_value < 0.05, "#BC3C29", "#0072B5"),
         col = "#0072B5", cex = 1.5, lwd = 1.5)
  for (i in 1:nrow(dt_hr_grid)) {
    if (dt_hr_grid$p_value[i] < 0.05) {
      text(dt_hr_grid$threshold_cut[i], dt_hr_grid$hr_upper[i] + 0.03,
           sprintf("p=%.3f", dt_hr_grid$p_value[i]),
           cex = 0.75, col = "#BC3C29", font = 2)
    }
  }
  abline(h = 1, col = "#E74C3C", lty = 2, lwd = 1.5)
  text(0.38, 0.72, "Protective", col = "#2E8B57", cex = 0.9, font = 2)
  text(0.78, 1.10, "Null effect", col = "#666666", cex = 0.9, font = 2)
  legend("topright",
         legend = c("HR (95% CI)", "p < 0.05", "Null (HR=1)"),
         pch = c(21, 21, NA),
         pt.bg = c("#0072B5", "#BC3C29", NA),
         lty = c(NA, NA, 2),
         lwd = c(NA, NA, 1.5),
         col = c("#0072B5", "#BC3C29", "#E74C3C"),
         bty = "n", cex = 0.9)
  mtext("SEQTaRget CCW + DAG14 + cluster-robust SE + m=20 Rubin's rules",
        side = 3, line = 0.3, cex = 0.7, col = "#666666", adj = 0)
  dev.off()
  log_info("PDF 矢量图已保存")
} else {
  log_info("[警告] Part A 无有效结果，跳过 Part C 可视化")
}

# ============================================================================
# 8. 汇总输出
# ============================================================================
log_section("分析汇总")

cat("\n============================================================\n")
cat("  V28 磷酸保护性阈值分析完成\n")
cat("============================================================\n")
cat(sprintf("  多切点结果: %s\n", threshold_csv))
cat(sprintf("  HR 网格:    %s\n", file.path(OUT_DIR, "v28_phosphate_threshold_hr_grid.csv")))
cat(sprintf("  效应曲线:   %s\n", file.path(FIG_DIR, "FigS_Phosphate_Threshold_HR_Curve.png")))
cat(sprintf("  日志:       %s\n", log_file))
cat("\n")

# 打印多切点结果表
cat("\n多切点亚组分析结果:\n")
cat(sprintf("%-15s %8s %18s %10s %8s %10s\n",
            "Threshold", "HR", "95% CI", "p", "n", "events"))
cat(paste(rep("-", 75), collapse = ""), "\n")
for (th in THRESHOLDS) {
  key <- as.character(th)
  if (key %in% names(threshold_results)) {
    r <- threshold_results[[key]]
    sig <- ifelse(r$p_value < 0.05, " *", "")
    cat(sprintf("%-15s %8.4f %7.4f-%.4f %10.4f %8d %10d%s\n",
                sprintf("< %.2f", th), r$hr, r$ci_lower, r$ci_upper,
                r$p_value, r$n, r$n_events, sig))
  }
}

# 关闭日志
close(log_con)
