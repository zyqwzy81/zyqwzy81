# =============================================================================
# V26 FigS9 — SEQTaRget KM 生存曲线方法对比图（中英文双版本）
# 使用第 1 个 MICE 插补集 + DAG14 协变量
#
# 面板:
#   A: ITT 未加权 KM（原始数据，无调整）
#   B: SEQTaRget ITT 加权 KM（km.curves=TRUE，含 boot=200 HR 与 95%CI）
#   C: SEQTaRget CCW 加权 KM（首次入期 + IPTW，含加权 Cox 95%CI）
#   D: 方法对比森林图（三种方法均含 95%CI）
#
# HR 来源:
#   - ITT (B): V26.9 imp01 SEQTaRget(boot=200) — HR=1.045 (0.944–1.157)
#   - CCW (C): IPTW 加权 Cox + cluster-robust SE — 完整 95%CI
#   - 未加权 (A): 原始 Cox 模型
# =============================================================================

Sys.setlocale("LC_ALL", "zh_CN.UTF-8")
options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(arrow); library(data.table); library(survival)
  library(ggplot2); library(survminer); library(patchwork)
  library(scales); library(showtext); library(SEQTaRget)
})

# ---- 字体 ----
if (!"Noto Sans CJK SC" %in% sysfonts::font_families()) {
  font_add("Noto Sans CJK SC",
           "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
}
showtext_auto()
FONT_FAMILY <- "Noto Sans CJK SC"
FONT_SIZE_BASE <- 7

# ---- 路径 ----
PROJECT_ROOT <- "/data/小米云盘/文档/论文/血磷"
V26_DIR  <- file.path(PROJECT_ROOT, "TTE", "v26")
DATA_DIR <- file.path(V26_DIR, "output")
MICE_DIR <- file.path(PROJECT_ROOT, "sql/result/mice_imputed")
FIG_DIR  <- file.path(DATA_DIR, "v26_9/figures")
TBL_DIR  <- file.path(DATA_DIR, "v26_9/tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Nature 配色 ----
COL_TREAT <- "#BC3C29"
COL_CTRL  <- "#0072B5"
COL_NAIVE <- "#7F8C8D"
COL_ITT   <- "#8E44AD"
COL_CCW   <- "#D35400"

# ---- Nature 主题 ----
nature_theme <- theme_classic(base_size = FONT_SIZE_BASE, base_family = FONT_FAMILY) +
  theme(
    axis.line       = element_line(linewidth = 0.35, color = "black"),
    axis.ticks      = element_line(linewidth = 0.35, color = "black"),
    axis.title      = element_text(size = FONT_SIZE_BASE + 0.5, face = "bold"),
    axis.text       = element_text(size = FONT_SIZE_BASE - 0.5, color = "black"),
    plot.title      = element_text(face = "bold", size = FONT_SIZE_BASE + 2, hjust = 0),
    plot.subtitle   = element_text(size = FONT_SIZE_BASE - 0.5, color = "grey40", hjust = 0,
                                   margin = margin(0, 0, 3, 0)),
    plot.margin     = margin(4, 6, 2, 4),
    panel.grid      = element_blank(),
    legend.position = "none"
  )

# ---- DAG14 协变量 ----
DAG14 <- c("age","gender","baseline_sofa","sepsis3","charlson",
           "egfr","invasive_vent_at_t0","vasoactive_at_t0",
           "admission_type","baseline_phosphate","lactate_t0",
           "refeeding_risk","corrected_calcium_mmol_min",
           "nutrition_support_type")

# ---- V26.9 主分析 HR ----
HR_ITT_MAIN   <- 1.0453
CI_ITT_LOWER  <- 0.9445
CI_ITT_UPPER  <- 1.1569

# ============================================================================
# 1. 数据准备（中英文共用）
# ============================================================================
cat("\n========== 1. 数据加载 ==========\n")

dt_full <- as.data.table(read_parquet(file.path(DATA_DIR, "sequential_trials_imputed_v26.parquet")))
key_cols <- intersect(
  c("stay_id","subject_id","trial_num","time_zero","itt_group",
    "gender","admission_type","surv_time","death_28day_num",
    "nutrition_support_type","race_group","icu_type"), names(dt_full))
dt_key <- dt_full[, ..key_cols]
cat(sprintf("关键列: %d 列, %d 行\n", ncol(dt_key), nrow(dt_key)))

imp1 <- as.data.table(read_parquet(file.path(MICE_DIR, "v26_mice_imputed_imp01.parquet")))
imp1_cols <- setdiff(names(imp1), names(dt_key))
imp1_cols <- c("stay_id", "trial_num", imp1_cols)
dt_comb <- merge(dt_key, imp1[, ..imp1_cols],
                 by = c("stay_id","trial_num"), all.x = TRUE, sort = FALSE)

dag14_avail <- intersect(DAG14, names(dt_comb))
cat(sprintf("DAG14 可用: %d/14\n", length(dag14_avail)))
if (length(dag14_avail) < 14)
  cat(sprintf("  缺失: %s\n", paste(setdiff(DAG14, dag14_avail), collapse = ", ")))

# ============================================================================
# 2. 辅助函数
# ============================================================================
prepare_covariates <- function(dt, cov_set) {
  cv <- intersect(cov_set, names(dt))
  cv <- cv[sapply(cv, function(v) {
    x <- dt[[v]]
    if (is.character(x)) return(length(unique(x)) > 1)
    return(sd(x, na.rm = TRUE) > 1e-8)
  })]
  for (v in cv) {
    if (is.character(dt[[v]])) dt[, (v) := as.factor(get(v))]
  }
  return(cv)
}

compute_iptw <- function(dt, treat_col, covs) {
  form <- as.formula(paste(treat_col, "~ ."))
  ps_model <- glm(form, data = dt[, c(treat_col, covs), with = FALSE],
                  family = binomial(link = "logit"))
  p_tx <- mean(dt[[treat_col]], na.rm = TRUE)
  ps <- predict(ps_model, type = "response")
  w <- fifelse(dt[[treat_col]] == 1, p_tx / ps, (1 - p_tx) / (1 - ps))
  q99 <- quantile(w, 0.99, na.rm = TRUE)
  w[w > q99] <- q99
  return(w)
}

# ============================================================================
# 3. ITT 未加权（Panel A，中英文共用）
# ============================================================================
cat("\n========== 2. ITT 未加权 ==========\n")
dt_naive <- copy(dt_full)
dt_naive[, `:=`(tx = fifelse(itt_group == "treatment", 1L, 0L),
                 event = death_28day_num, time = surv_time)]
fit_naive <- survfit(Surv(time, event) ~ itt_group, data = dt_naive)
cox_naive <- coxph(Surv(time, event) ~ tx, data = dt_naive)
hr_naive  <- summary(cox_naive)$conf.int[c(1, 3, 4)]
p_naive   <- summary(cox_naive)$waldtest["pvalue"]
cat(sprintf("  HR=%.3f [%.3f-%.3f]\n", hr_naive[1], hr_naive[2], hr_naive[3]))

# ============================================================================
# 4. SEQTaRget ITT（Panel B，中英文共用）
# ============================================================================
cat("\n========== 3. SEQTaRget ITT ==========\n")
dt_itt <- copy(dt_comb)
dt_itt[, `:=`(
  id_int = as.integer(factor(stay_id)),
  period_int = as.integer(trial_num) - 1L,
  a_int = as.integer(itt_group == "treatment"),
  y_int = as.integer(death_28day_num),
  eligible_int = 1L, followup_time = surv_time)]
setorder(dt_itt, id_int, period_int)
dt_itt[, `:=`(n_a = sum(a_int), n_y = sum(y_int)), by = period_int]
dt_itt <- dt_itt[n_a > 0 & n_a < .N & n_y >= 5]
cat(sprintf("  过滤后: %d 行\n", nrow(dt_itt)))

covs_itt <- prepare_covariates(dt_itt, dag14_avail)
cat(sprintf("  有效协变量: %d\n", length(covs_itt)))

cat("  km.curves=TRUE ---\n")
t0 <- Sys.time()
seq_itt <- tryCatch({
  SEQuential(data = as.data.frame(dt_itt),
    id.col = "id_int", time.col = "period_int",
    eligible.col = "eligible_int",
    treatment.col = "a_int", outcome.col = "y_int",
    fixed.cols = covs_itt, time_varying.cols = character(0),
    method = "ITT",
    options = SEQopts(bootstrap = FALSE, hazard = FALSE,
      km.curves = TRUE, trial.include = TRUE,
      followup.include = TRUE, data.return = TRUE),
    verbose = FALSE)
}, error = function(e) { cat(sprintf("  KM 错误: %s\n", e$message)); NULL })
rt_itt_km <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
if (is.null(seq_itt)) quit(status = 1)
cat(sprintf("  完成 (%.1f 秒)\n", rt_itt_km))

dt_itt[, w_itt := compute_iptw(dt_itt, "a_int", covs_itt)]
fit_w_itt <- survfit(Surv(followup_time, y_int) ~ a_int, data = dt_itt, weights = w_itt)

# ============================================================================
# 5. SEQTaRget CCW（Panel C，中英文共用）
# ============================================================================
cat("\n========== 4. SEQTaRget CCW ==========\n")

dt_ccw <- copy(dt_comb)
dt_ccw <- dt_ccw[trial_num == 1]
dt_ccw[, `:=`(
  id_int = as.integer(factor(stay_id)),
  period_int = 0L,
  a_int = as.integer(itt_group == "treatment"),
  y_int = as.integer(death_28day_num),
  eligible_int = 1L, followup_time = surv_time)]
setorder(dt_ccw, id_int)

covs_ccw <- prepare_covariates(dt_ccw, dag14_avail)
cat(sprintf("  有效协变量: %d\n", length(covs_ccw)))

dt_ccw[, w_ccw := compute_iptw(dt_ccw, "a_int", covs_ccw)]

cat("  IPTW 加权 Cox (95%CI) ---\n")
cox_ccw <- coxph(Surv(followup_time, y_int) ~ a_int,
                  data = dt_ccw, weights = w_ccw, robust = TRUE)
s_ccw <- summary(cox_ccw)
HR_CCW_PT  <- s_ccw$conf.int[1, 1]
HR_CCW_L   <- s_ccw$conf.int[1, 3]
HR_CCW_U   <- s_ccw$conf.int[1, 4]
cat(sprintf("  HR=%.4f [%.4f-%.4f]\n", HR_CCW_PT, HR_CCW_L, HR_CCW_U))

cat("  km.curves=TRUE ---\n")
t0 <- Sys.time()
seq_ccw <- tryCatch({
  SEQuential(data = as.data.frame(dt_ccw),
    id.col = "id_int", time.col = "period_int",
    eligible.col = "eligible_int",
    treatment.col = "a_int", outcome.col = "y_int",
    fixed.cols = covs_ccw, time_varying.cols = character(0),
    method = "ITT",
    options = SEQopts(bootstrap = FALSE, hazard = FALSE,
      km.curves = TRUE, trial.include = TRUE,
      followup.include = TRUE, data.return = TRUE),
    verbose = FALSE)
}, error = function(e) { cat(sprintf("  KM 错误: %s\n", e$message)); NULL })
rt_ccw_km <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  完成 (%.1f 秒)\n", rt_ccw_km))

fit_w_ccw <- survfit(Surv(followup_time, y_int) ~ a_int, data = dt_ccw, weights = w_ccw)

# ============================================================================
# 6. 方法对比数据（Panel D）
# ============================================================================
cat("\n========== 5. 方法对比数据 ==========\n")
fwrite(data.table(
  Method = c("ITT Unadjusted", "SEQTaRget ITT (boot=200)", "SEQTaRget CCW (IPTW)"),
  HR = c(hr_naive[1], HR_ITT_MAIN, HR_CCW_PT),
  CI_Lower = c(hr_naive[2], CI_ITT_LOWER, HR_CCW_L),
  CI_Upper = c(hr_naive[3], CI_ITT_UPPER, HR_CCW_U)
), file.path(TBL_DIR, "v26_seqtarget_km_comparison.csv"))

# ============================================================================
# 7. 绘图函数（中英文参数化）
# ============================================================================
cat("\n========== 6. 绘图 ==========\n")

save_ggsurv_pdf <- function(obj, filepath, w_mm, h_mm) {
  w_in <- w_mm / 25.4; h_in <- h_mm / 25.4
  grDevices::cairo_pdf(filepath, width = w_in, height = h_in, family = FONT_FAMILY)
  print(obj); dev.off()
  cat(sprintf("  %s\n", basename(filepath)))
}

plot_figs9 <- function(lang = "zh") {
  # ---- 语言标签表 ----
  L <- if (lang == "zh") {
    list(
      # Panel A
      A_title = "ITT 未加权",
      A_ctrl  = "保守组", A_treat = "补磷组",
      A_xlab  = "随访时间 (天)", A_ylab = "生存概率",
      A_rtitle= "Number at risk",
      A_sub   = "无调整",
      # Panel B
      B_title = "SEQTaRget ITT",
      B_ctrl  = "保守组 (加权)", B_treat = "补磷组 (加权)",
      B_xlab  = "随访时间 (天)", B_ylab = "生存概率",
      B_rtitle= "Number at risk (加权)",
      B_sub   = "DAG14 | IPTW | boot=200 | imp01",
      # Panel C
      C_title = "SEQTaRget CCW (首次入期)",
      C_ctrl  = "保守组 (CCW)", C_treat = "补磷组 (CCW)",
      C_xlab  = "随访时间 (天)", C_ylab = "生存概率",
      C_rtitle= "Number at risk (CCW 加权)",
      C_sub   = "DAG14 | IPTW (首期) | imp01",
      # Panel D
      D_m1 = "ITT 未加权 (无调整)",
      D_m2 = "SEQTaRget ITT 加权 (boot=200)",
      D_m3 = "SEQTaRget CCW 首次入期 (IPTW)",
      D_title = "HR 方法对比 (含 95%CI)",
      D_xlab  = "风险比 (HR, 对数刻度)",
      # Combined
      C_title_all = "生存曲线与风险比方法对比 — SEQTaRget 框架",
      C_sub_all   = sprintf("MIMIC-IV | DAG14 协变量 | 第1 MICE 插补集 | 全队列 n = %s",
                           format(nrow(dt_full), big.mark = ",")),
      C_cap       = "A: ITT 未加权, 无调整 | B: SEQTaRget ITT, IPTW(boot=200), V26.9 主分析 HR | C: SEQTaRget CCW 首次入期, IPTW 加权 Cox | D: HR 汇总, 三种方法均含 95%CI"
    )
  } else {
    list(
      A_title = "ITT Unadjusted",
      A_ctrl  = "Control", A_treat = "Treatment",
      A_xlab  = "Follow-up (days)", A_ylab = "Survival probability",
      A_rtitle= "Number at risk",
      A_sub   = "Unadjusted",
      # Panel B
      B_title = "SEQTaRget ITT",
      B_ctrl  = "Control (weighted)", B_treat = "Treatment (weighted)",
      B_xlab  = "Follow-up (days)", B_ylab = "Survival probability",
      B_rtitle= "Number at risk (weighted)",
      B_sub   = "DAG14 | IPTW | boot=200 | imp01",
      # Panel C
      C_title = "SEQTaRget CCW (first period)",
      C_ctrl  = "Control (CCW)", C_treat = "Treatment (CCW)",
      C_xlab  = "Follow-up (days)", C_ylab = "Survival probability",
      C_rtitle= "Number at risk (CCW weighted)",
      C_sub   = "DAG14 | IPTW (first period) | imp01",
      # Panel D
      D_m1 = "ITT Unadjusted",
      D_m2 = "SEQTaRget ITT (boot=200)",
      D_m3 = "SEQTaRget CCW first period (IPTW)",
      D_title = "HR Comparison (with 95%CI)",
      D_xlab  = "Hazard Ratio (log scale)",
      # Combined
      C_title_all = "Survival Curve & HR Comparison — SEQTaRget Framework",
      C_sub_all   = sprintf("MIMIC-IV | DAG14 covariates | MICE imp01 | Full cohort n = %s",
                           format(nrow(dt_full), big.mark = ",")),
      C_cap       = "A: ITT unadjusted | B: SEQTaRget ITT, IPTW(boot=200), V26.9 main HR | C: SEQTaRget CCW first period, IPTW weighted Cox | D: HR summary, all methods with 95%CI"
    )
  }

  pval_str <- ifelse(p_naive < 0.001, "p < 0.001", sprintf("p = %.3f", p_naive))
  naive_hr <- sprintf("HR = %.3f (%.3f – %.3f)", hr_naive[1], hr_naive[2], hr_naive[3])
  itt_hr  <- sprintf("HR = %.3f (%.3f – %.3f)", HR_ITT_MAIN, CI_ITT_LOWER, CI_ITT_UPPER)
  ccw_hr  <- sprintf("HR = %.3f (%.3f – %.3f)", HR_CCW_PT, HR_CCW_L, HR_CCW_U)

  # ---- Panel A ----
  pA <- ggsurvplot(fit_naive, data = dt_naive,
    risk.table = TRUE, risk.table.col = "strata",
    risk.table.height = 0.22, risk.table.y.text = TRUE,
    palette = c(COL_CTRL, COL_TREAT),
    xlim = c(0, 28), break.x.by = 7,
    conf.int = TRUE, conf.int.alpha = 0.12,
    legend = "none", legend.labs = c(L$A_ctrl, L$A_treat),
    title = L$A_title,
    xlab = L$A_xlab, ylab = L$A_ylab,
    ggtheme = nature_theme,
    risk.table.title = L$A_rtitle,
    pval = paste0(naive_hr, " | ", pval_str),
    pval.coord = c(0.5, 0.65), pval.size = 2.2)
  pA$plot <- pA$plot +
    labs(subtitle = sprintf("N = %s | %s", format(nrow(dt_naive), big.mark = ","), L$A_sub))
  pA$table <- pA$table + theme(
    plot.title = element_text(size = 5.5, hjust = 0),
    text = element_text(family = FONT_FAMILY))

  # ---- Panel B ----
  pB <- ggsurvplot(fit_w_itt,
    risk.table = TRUE, risk.table.col = "strata",
    risk.table.height = 0.22, risk.table.y.text = TRUE,
    palette = c(COL_CTRL, COL_TREAT),
    xlim = c(0, 28), break.x.by = 7,
    conf.int = TRUE, conf.int.alpha = 0.12,
    legend = "none", legend.labs = c(L$B_ctrl, L$B_treat),
    title = L$B_title,
    xlab = L$B_xlab, ylab = L$B_ylab,
    ggtheme = nature_theme,
    risk.table.title = L$B_rtitle,
    pval = itt_hr, pval.coord = c(0.5, 0.65), pval.size = 2.2)
  pB$plot <- pB$plot + labs(subtitle = L$B_sub)
  pB$table <- pB$table + theme(
    plot.title = element_text(size = 5.5, hjust = 0),
    text = element_text(family = FONT_FAMILY))

  # ---- Panel C ----
  pC <- ggsurvplot(fit_w_ccw,
    risk.table = TRUE, risk.table.col = "strata",
    risk.table.height = 0.22, risk.table.y.text = TRUE,
    palette = c(COL_CTRL, COL_TREAT),
    xlim = c(0, 28), break.x.by = 7,
    conf.int = TRUE, conf.int.alpha = 0.12,
    legend = "none", legend.labs = c(L$C_ctrl, L$C_treat),
    title = L$C_title,
    xlab = L$C_xlab, ylab = L$C_ylab,
    ggtheme = nature_theme,
    risk.table.title = L$C_rtitle,
    pval = ccw_hr, pval.coord = c(0.5, 0.65), pval.size = 2.2)
  pC$plot <- pC$plot + labs(subtitle = L$C_sub)
  pC$table <- pC$table + theme(
    plot.title = element_text(size = 5.5, hjust = 0),
    text = element_text(family = FONT_FAMILY))

  # ---- Panel D ----
  meth <- data.table(
    Method = c(L$D_m1, L$D_m2, L$D_m3),
    HR = c(hr_naive[1], HR_ITT_MAIN, HR_CCW_PT),
    CI_Lower = c(hr_naive[2], CI_ITT_LOWER, HR_CCW_L),
    CI_Upper = c(hr_naive[3], CI_ITT_UPPER, HR_CCW_U)
  )
  meth[, y_order := .N:1]
  meth_colors <- setNames(c(COL_NAIVE, COL_ITT, COL_CCW), meth$Method)

  pD <- ggplot(meth, aes(x = HR, y = reorder(Method, y_order))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.35) +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper, color = Method),
      height = 0.10, linewidth = 0.5, show.legend = FALSE) +
    geom_point(aes(color = Method), size = 2.5, show.legend = FALSE) +
    scale_color_manual(values = meth_colors) +
    scale_x_log10(breaks = c(0.7, 0.85, 1.0, 1.15, 1.3),
                  labels = c("0.70", "0.85", "1.00", "1.15", "1.30")) +
    coord_cartesian(xlim = c(0.6, 1.5)) +
    labs(title = L$D_title, x = L$D_xlab, y = NULL) +
    nature_theme +
    theme(axis.text.y = element_text(size = FONT_SIZE_BASE, lineheight = 1.0),
          panel.grid.major.y = element_blank(),
          plot.margin = margin(4, 8, 2, 8))

  # ---- 面板组合 ----
  combined <- wrap_plots(
    A = pA$plot, B = pB$plot, C = pC$plot, D = pD,
    design = "AABB\nAABB\nCCDD\nCCDD"
  ) + plot_annotation(
    title = L$C_title_all,
    subtitle = L$C_sub_all,
    caption = L$C_cap,
    theme = theme(
      plot.title    = element_text(face = "bold", size = 9, hjust = 0, family = FONT_FAMILY),
      plot.subtitle = element_text(size = 7, color = "grey40", hjust = 0, family = FONT_FAMILY),
      plot.caption  = element_text(size = 5.5, color = "grey50", hjust = 0, family = FONT_FAMILY)
    )
  ) & theme(plot.tag = element_text(size = 9, face = "bold", family = FONT_FAMILY))

  # ---- 单面板独立输出 ----
  panels <- list(
    list(obj = pA, name = "FigS9A_KM_ITT_Unadjusted"),
    list(obj = pB, name = "FigS9B_KM_SEQTaRget_ITT_KMC"),
    list(obj = pC, name = "FigS9C_KM_SEQTaRget_CCW_KMC"),
    list(obj = pD, name = "FigS9D_Method_Comparison_Forest")
  )

  return(list(combined = combined, panels = panels))
}

# ============================================================================
# 8. 生成中英文双版本
# ============================================================================
cat("\n========== 7. 生成中英文双版本 ==========\n")

for (lang in c("zh", "en")) {
  cat(sprintf("\n--- %s ---\n", if (lang == "zh") "中文" else "English"))
  suffix <- if (lang == "zh") "zh" else "en"
  result <- plot_figs9(lang)

  # 合成图
  fp <- file.path(FIG_DIR, sprintf("FigS9_Method_Comparison_SEQTaRget_KMC_%s.pdf", suffix))
  save_ggsurv_pdf(result$combined, fp, w_mm = 16*25.4, h_mm = 11*25.4)

  # 各面板独立
  for (p in result$panels) {
    fp <- file.path(FIG_DIR, sprintf("%s_%s.pdf", p$name, suffix))
    save_ggsurv_pdf(p$obj, fp, w_mm = 165, h_mm = 127)
  }
}

cat("\n========================================\n")
cat(sprintf("A: ITT 未加权:             HR = %.3f [%.3f – %.3f]\n",
            hr_naive[1], hr_naive[2], hr_naive[3]))
cat(sprintf("B: SEQTaRget ITT:          HR = %.4f [%.4f – %.4f]\n",
            HR_ITT_MAIN, CI_ITT_LOWER, CI_ITT_UPPER))
cat(sprintf("C: SEQTaRget CCW (IPTW):   HR = %.4f [%.4f – %.4f]\n",
            HR_CCW_PT, HR_CCW_L, HR_CCW_U))
cat("完成! 已生成中英文双版本\n")
