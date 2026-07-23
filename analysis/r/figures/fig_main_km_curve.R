# =============================================================================
# V26 主图 — SEQTaRget ITT 加权 KM 生存曲线（中英文双版本）
# 使用第 1 个 MICE 插补集 + DAG14 协变量
#
# 用途: 插入论文主结果部分，展示 SEQTaRget ITT 的加权生存曲线与完整 CI
#
# HR 来源: V26.9 imp01 SEQTaRget(boot=200) — HR=1.045 (0.945–1.157)
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
FONT_BASE   <- 8

# ---- 路径 ----
PROJECT_ROOT <- "/data/小米云盘/文档/论文/血磷"
V26_DIR  <- file.path(PROJECT_ROOT, "TTE", "v26")
DATA_DIR <- file.path(V26_DIR, "output")
MICE_DIR <- file.path(PROJECT_ROOT, "sql/result/mice_imputed")
FIG_DIR  <- file.path(DATA_DIR, "v26_9/figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Nature 配色 ----
COL_TREAT <- "#BC3C29"
COL_CTRL  <- "#0072B5"

# ---- DAG14 协变量 ----
DAG14 <- c("age","gender","baseline_sofa","sepsis3","charlson",
           "egfr","invasive_vent_at_t0","vasoactive_at_t0",
           "admission_type","baseline_phosphate","lactate_t0",
           "refeeding_risk","corrected_calcium_mmol_min",
           "nutrition_support_type")

# ---- V26.9 主分析 HR ----
HR_MAIN  <- 1.0453
CI_LOWER <- 0.9445
CI_UPPER <- 1.1569

# ============================================================================
# 1. 数据准备
# ============================================================================
cat("\n========== 1. 数据加载 ==========\n")

dt_full <- as.data.table(read_parquet(file.path(DATA_DIR, "sequential_trials_imputed_v26.parquet")))
key_cols <- intersect(
  c("stay_id","subject_id","trial_num","time_zero","itt_group",
    "gender","admission_type","surv_time","death_28day_num",
    "nutrition_support_type"), names(dt_full))
dt_key <- dt_full[, ..key_cols]

imp1 <- as.data.table(read_parquet(file.path(MICE_DIR, "v26_mice_imputed_imp01.parquet")))
imp1_cols <- setdiff(names(imp1), names(dt_key))
imp1_cols <- c("stay_id", "trial_num", imp1_cols)
dt <- merge(dt_key, imp1[, ..imp1_cols],
            by = c("stay_id","trial_num"), all.x = TRUE, sort = FALSE)

dag14_avail <- intersect(DAG14, names(dt))
cat(sprintf("DAG14 可用: %d/14\n", length(dag14_avail)))

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
# 3. SEQTaRget ITT 数据准备（仅一次）
# ============================================================================
cat("\n========== 2. SEQTaRget ITT 数据 ==========\n")
dt_itt <- copy(dt)
dt_itt[, `:=`(
  id_int = as.integer(factor(stay_id)),
  period_int = as.integer(trial_num) - 1L,
  a_int = as.integer(itt_group == "treatment"),
  y_int = as.integer(death_28day_num),
  eligible_int = 1L, followup_time = surv_time)]
setorder(dt_itt, id_int, period_int)
dt_itt[, `:=`(n_a = sum(a_int), n_y = sum(y_int)), by = period_int]
dt_itt <- dt_itt[n_a > 0 & n_a < .N & n_y >= 5]
cat(sprintf("  过滤后: %d 行, %d 个患者\n",
            nrow(dt_itt), uniqueN(dt_itt$stay_id)))

covs <- prepare_covariates(dt_itt, dag14_avail)
cat(sprintf("  有效协变量: %d\n", length(covs)))

# IPTW 权重
dt_itt[, w_itt := compute_iptw(dt_itt, "a_int", covs)]

# SEQTaRget km.curves=TRUE
cat("  SEQTaRget km.curves=TRUE ...\n")
t0 <- Sys.time()
seq_itt <- tryCatch({
  SEQuential(data = as.data.frame(dt_itt),
    id.col = "id_int", time.col = "period_int",
    eligible.col = "eligible_int",
    treatment.col = "a_int", outcome.col = "y_int",
    fixed.cols = covs, time_varying.cols = character(0),
    method = "ITT",
    options = SEQopts(bootstrap = FALSE, hazard = FALSE,
      km.curves = TRUE, trial.include = TRUE,
      followup.include = TRUE, data.return = TRUE),
    verbose = FALSE)
}, error = function(e) { cat(sprintf("  错误: %s\n", e$message)); NULL })
rt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
if (is.null(seq_itt)) quit(status = 1)
cat(sprintf("  完成 (%.1f 秒)\n", rt))

fit_w <- survfit(Surv(followup_time, y_int) ~ a_int, data = dt_itt, weights = w_itt)

# ============================================================================
# 4. 绘图函数（中英文参数化）
# ============================================================================
nature_theme <- theme_classic(base_size = FONT_BASE, base_family = FONT_FAMILY) +
  theme(
    axis.line  = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    axis.title = element_text(size = FONT_BASE + 1, face = "bold"),
    axis.text  = element_text(size = FONT_BASE, color = "black"),
    plot.title = element_text(face = "bold", size = FONT_BASE + 3, hjust = 0),
    plot.subtitle = element_text(size = FONT_BASE, color = "grey40", hjust = 0,
                                 margin = margin(0, 0, 6, 0)),
    plot.caption = element_text(size = FONT_BASE - 1.5, color = "grey45", hjust = 0),
    plot.margin = margin(6, 8, 4, 6),
    panel.grid = element_blank()
  )

plot_km_curve <- function(lang = "zh", n_total) {
  # ---- 语言标签表 ----
  labels <- if (lang == "zh") {
    list(
      title      = "静脉补磷 vs 保守管理 — 28天生存曲线",
      treat_lab  = "静脉补磷组",
      ctrl_lab   = "保守管理组",
      xlab       = "随访时间 (天)",
      ylab       = "生存概率",
      rt_title   = "Number at risk (IPTW 加权)",
      method_note = "CCW 框架 | DAG14 | IPTW 稳定权重 | 200次 Bootstrap",
      caption    = sprintf("数据来源: MIMIC-IV v3.1 (2008–2022) | n = %s trial-periods",
                          format(n_total, big.mark = ",")),
      hr_prefix  = "SEQTaRget ITT"
    )
  } else {
    list(
      title      = "IV Phosphate vs Conservative Management — 28-Day Survival",
      treat_lab  = "IV phosphate",
      ctrl_lab   = "Conservative management",
      xlab       = "Follow-up (days)",
      ylab       = "Survival probability",
      rt_title   = "Number at risk (IPTW weighted)",
      method_note = "CCW framework | DAG14 | IPTW stabilized weights | 200 bootstrap",
      caption    = sprintf("Source: MIMIC-IV v3.1 (2008–2022) | n = %s trial-periods",
                          format(n_total, big.mark = ",")),
      hr_prefix  = "SEQTaRget ITT"
    )
  }

  hr_str <- sprintf("%s: HR = %.3f (95%% CI %.3f – %.3f)",
                    labels$hr_prefix, HR_MAIN, CI_LOWER, CI_UPPER)

  p <- ggsurvplot(fit_w,
    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.height = 0.20,
    risk.table.y.text = TRUE,
    palette = c(COL_CTRL, COL_TREAT),
    xlim = c(0, 28),
    break.x.by = 7,
    conf.int = TRUE,
    conf.int.alpha = 0.18,
    legend.title = "",
    legend.labs = c(labels$ctrl_lab, labels$treat_lab),
    legend = c(0.85, 0.95),
    title = labels$title,
    subtitle = hr_str,
    xlab = labels$xlab,
    ylab = labels$ylab,
    ggtheme = nature_theme,
    risk.table.title = labels$rt_title,
    risk.table.fontsize = 3.2,
    surv.median.line = "none",
    pval = FALSE)

  p$plot <- p$plot +
    annotate("text", x = 14, y = 0.25,
             label = labels$method_note,
             size = 2.8, color = "grey45", hjust = 0.5) +
    labs(caption = labels$caption)

  p$table <- p$table + theme(
    plot.title = element_text(size = 7, hjust = 0, family = FONT_FAMILY),
    text = element_text(family = FONT_FAMILY))
  return(p)
}

# ============================================================================
# 5. 保存函数
# ============================================================================
save_ggsurv <- function(obj, filename_pdf, filename_png, w_mm, h_mm) {
  w_in <- w_mm / 25.4; h_in <- h_mm / 25.4
  grDevices::cairo_pdf(filename_pdf, width = w_in, height = h_in, family = FONT_FAMILY)
  print(obj); dev.off()
  png(filename_png, width = w_in, height = h_in, units = "in", res = 300)
  print(obj); dev.off()
  cat(sprintf("  %s (%.0f × %.0f mm)\n", basename(filename_pdf), w_mm, h_mm))
}

# ============================================================================
# 6. 生成中英文双版本
# ============================================================================
cat("\n========== 3. 生成中英文双版本 ==========\n")

for (lang in c("zh", "en")) {
  cat(sprintf("\n--- %s ---\n", if (lang == "zh") "中文" else "English"))
  suffix <- if (lang == "zh") "zh" else "en"
  p <- plot_km_curve(lang = lang, n_total = nrow(dt_itt))

  # 单栏 89mm
  save_ggsurv(p,
    file.path(FIG_DIR, sprintf("Fig_KM_ITT_Primary_%s.pdf", suffix)),
    file.path(FIG_DIR, sprintf("Fig_KM_ITT_Primary_%s.png", suffix)),
    w_mm = 89, h_mm = 120)

  # 全宽 183mm
  save_ggsurv(p,
    file.path(FIG_DIR, sprintf("Fig_KM_ITT_Primary_FullWidth_%s.pdf", suffix)),
    file.path(FIG_DIR, sprintf("Fig_KM_ITT_Primary_FullWidth_%s.png", suffix)),
    w_mm = 183, h_mm = 140)
}

cat(sprintf("\nHR = %.4f [%.4f – %.4f]\n", HR_MAIN, CI_LOWER, CI_UPPER))
cat("完成! 已生成中英文双版本\n")
