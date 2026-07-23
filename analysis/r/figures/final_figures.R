# ============================================================================
# V26.9 发表级图表统一生成脚本
# 后端: R (ggplot2 + patchwork + survminer)
# 输出: 全部主图 + 补充图, 600 DPI PNG + 矢量 PDF 双格式
# 语言: 纯英文
# 规范: nature-figure + CCM 杂志标准
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(scales)
  library(patchwork)
  library(survival)
  library(survminer)
})

# ---- 路径配置 ----
OUT_DIR  <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9"
FIG_DIR  <- file.path(OUT_DIR, "figures")
TBL_DIR  <- file.path(OUT_DIR, "tables")
TBL_DIR_V26 <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/tables"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 字体 ----
FONT_FAMILY <- "sans"

# ---- JAMA/CCM 专业配色 ----
JAMA_RED     <- "#BC3C29"
JAMA_BLUE    <- "#0072B5"
FOREST_GREEN <- "#2E8B57"
DEEP_NAVY    <- "#1B4F72"
ORANGE_DOT   <- "#D55E00"
EVAL_WEAK    <- "#F39C12"
EVAL_MOD     <- "#E67E22"
EVAL_STRONG  <- "#27AE60"
GRAY_REF     <- "#BDC3C7"
GRAY_50      <- "#7F7F7F"

# ---- 统一 ggplot2 主题 ----
theme_jama <- theme_bw(base_size = 7, base_family = FONT_FAMILY) +
  theme(
    axis.line        = element_line(linewidth = 0.35, color = "black"),
    axis.ticks       = element_line(linewidth = 0.35, color = "black"),
    axis.title       = element_text(size = 9, face = "bold"),
    axis.text        = element_text(size = 7.5, color = "black"),
    plot.title       = element_text(face = "bold", size = 11, hjust = 0,
                                    margin = margin(0, 0, 2, 0)),
    plot.subtitle    = element_text(size = 8, color = "grey30", hjust = 0,
                                    margin = margin(0, 0, 6, 0)),
    plot.caption     = element_text(size = 6.5, color = "grey45", hjust = 0),
    plot.margin      = margin(8, 10, 6, 6),
    panel.grid       = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "grey92"),
    legend.position  = "bottom",
    legend.title     = element_text(size = 7),
    legend.text      = element_text(size = 6.5),
    legend.margin    = margin(2, 0, 2, 0),
    strip.text       = element_text(size = 8, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "grey80")
  )
theme_set(theme_jama)

# ---- 高 DPI 保存函数 ----
save_fig <- function(p, name, w = 183, h = 120, dpi = 600) {
  w_in <- w / 25.4
  h_in <- h / 25.4
  # PDF (矢量)
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p,
         width = w_in, height = h_in, device = grDevices::cairo_pdf)
  # PNG (600 DPI)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p,
         width = w_in, height = h_in, dpi = dpi, device = "png")
  cat(sprintf("  [OK] %s\n", name))
}

# ---- 安全读取 CSV ----
safe_fread <- function(path, default = NULL) {
  if (file.exists(path)) {
    fread(path)
  } else {
    cat(sprintf("  [WARN] Missing: %s\n", path))
    default
  }
}

# ============================================================================
# Fig 2: Covariate Balance Love Plot (数据驱动)
# ============================================================================
cat("\n===== Fig 2: Covariate Balance Love Plot =====\n")

smd <- safe_fread(file.path(TBL_DIR_V26, "v26_balance_smd.csv"))
if (!is.null(smd)) {
  var_name_map <- c(
    "基线血磷" = "Baseline phosphate",
    "基线 SOFA" = "Baseline SOFA",
    "血管活性药物" = "Vasoactive agents",
    "Sepsis-3" = "Sepsis-3",
    "有创通气" = "Invasive ventilation",
    "24h净液体平衡" = "24h net fluid balance",
    "S/F比值" = "S/F ratio",
    "心脏手术" = "Cardiac surgery",
    "FiO2默认值标记" = "FiO2 defaulted",
    "年龄" = "Age",
    "再喂养风险" = "Refeeding risk",
    "TPN热卡(T0前24h)" = "TPN calories (24h pre-T0)",
    "输血" = "Transfusion at T0",
    "最差HCO3" = "Worst HCO3",
    "碳酸氢钠输注" = "NaHCO3 infusion",
    "GNRI" = "GNRI",
    "EN热卡(T0前24h)" = "EN calories (24h pre-T0)",
    "Charlson合并症指数" = "Charlson index",
    "胰岛素" = "Insulin at T0",
    "利尿剂" = "Diuretic at T0",
    "类固醇" = "Steroid at T0",
    "eGFR" = "eGFR",
    "估算P/F" = "Estimated P/F",
    "最差pH" = "Worst pH",
    "校正钙(mmol/L)" = "Corrected calcium",
    "β2激动剂" = "Beta2-agonist"
  )
  smd[, var_label_en := var_name_map[var_label]]
  smd[is.na(var_label_en), var_label_en := var_label]
  smd <- smd[order(-abs(smd_before))]
  smd <- smd[1:min(14, .N)]
  smd[, var_label_en := factor(var_label_en, levels = rev(var_label_en))]

  smd_long <- melt(smd, id.vars = "var_label_en",
                   measure.vars = c("smd_before", "smd_after"),
                   variable.name = "stage", value.name = "smd")
  smd_long[, stage := fcase(
    stage == "smd_before", "Before IPTW",
    stage == "smd_after",  "After IPTW"
  )]
  smd_long[, abs_smd := abs(smd)]

  p_fig2 <- ggplot(smd_long, aes(x = abs_smd, y = var_label_en, color = stage, shape = stage)) +
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
    geom_vline(xintercept = 0.2, linetype = "dotted", color = JAMA_RED, linewidth = 0.4) +
    geom_point(size = 3.5) +
    scale_color_manual(values = c("Before IPTW" = ORANGE_DOT, "After IPTW" = JAMA_BLUE)) +
    scale_shape_manual(values = c("Before IPTW" = 16, "After IPTW" = 17)) +
    scale_x_continuous(limits = c(0, max(smd_long$abs_smd) * 1.05),
                       breaks = seq(0, 0.6, 0.1)) +
    labs(
      title = "Figure 2. Covariate Balance Before and After IPTW Weighting",
      subtitle = "CCCW framework | Top 14 covariates by pre-weighting |SMD|",
      x = "Absolute Standardized Mean Difference",
      y = NULL, color = NULL, shape = NULL,
      caption = "Dashed line: |SMD| = 0.1 (ideal). Dotted line: |SMD| = 0.2 (acceptable)."
    ) +
    theme(panel.grid.major.y = element_blank())
  save_fig(p_fig2, "Fig2_Love_Plot", 170, 120)
}

# ============================================================================
# Fig 3B: Primary & Sensitivity Analysis Forest Plot (数据驱动)
# ============================================================================
cat("\n===== Fig 3B: Primary & Sensitivity Forest Plot =====\n")

rubin <- safe_fread(file.path(TBL_DIR, "v26_9_seqtarget_rubin_results.csv"))
if (!is.null(rubin) && nrow(rubin) > 0) {
  primary_hr <- rubin$hr[1]
  primary_lo <- rubin$ci_lower[1]
  primary_hi <- rubin$ci_upper[1]

  fig3b_data <- data.table(
    label = c(
      "Conventional Cox (no CCW)",
      "SEQTaRget ITT — Primary (DAG14, m=20, boot=200)",
      "DAG13 Sensitivity (excl. sepsis3)",
      "IPTW-Cox (Rubin, m=5)",
      "Weibull AFT (Time Ratio)",
      "Per-Protocol CCW+IPTW (boot=50)",
      "First-Enrollment-Only (n=22,685)",
      "eICU External Validation (11 cov)",
      "Computational Robustness (m=5, boot=50)"
    ),
    hr = c(0.899, primary_hr, 0.969, 1.054, 0.989, 0.979, 1.047, 1.087, 1.057),
    ci_lower = c(0.826, primary_lo, 0.875, 0.948, 0.940, 0.885, 0.951, 0.913, 0.963),
    ci_upper = c(0.980, primary_hi, 1.072, 1.171, 1.041, 1.083, 1.154, 1.293, 1.159),
    type = c("conventional", "primary", "sensitivity", "sensitivity", "sensitivity",
             "sensitivity", "sensitivity", "external", "robustness")
  )
  fig3b_data$y <- nrow(fig3b_data):1

  p_fig3b <- ggplot(fig3b_data, aes(x = hr, y = y)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                   height = 0.15, linewidth = 0.9, color = DEEP_NAVY) +
    geom_point(aes(fill = type), size = 4.5, shape = 22, color = "black") +
    scale_fill_manual(
      values = c("primary" = JAMA_RED, "conventional" = ORANGE_DOT,
                 "sensitivity" = JAMA_BLUE, "external" = FOREST_GREEN,
                 "robustness" = "grey50"),
      guide = "none"
    ) +
    scale_y_continuous(breaks = fig3b_data$y, labels = rev(fig3b_data$label)) +
    scale_x_log10(breaks = c(0.80, 0.90, 1.00, 1.10, 1.20),
                  labels = c("0.80", "0.90", "1.00", "1.10", "1.20")) +
    labs(
      title = "Figure 3B. Primary and Sensitivity Analyses",
      subtitle = "SEQTaRget ITT framework | DAG14 covariates | Cloning-Censoring-Weighting",
      x = "Hazard Ratio (95% CI)", y = NULL,
      caption = "Primary ITT: Rubin's rules across m=20 MICE imputed datasets, 200 bootstrap iterations"
    ) +
    coord_cartesian(xlim = c(0.75, 1.35))

  ann_x <- 1.38
  for (i in seq_len(nrow(fig3b_data))) {
    p_fig3b <- p_fig3b + annotate("text", x = ann_x, y = i,
      label = sprintf("%.3f (%.3f\u2013%.3f)", fig3b_data$hr[i],
                      fig3b_data$ci_lower[i], fig3b_data$ci_upper[i]),
      hjust = 0, size = 3, color = if (i == 2) JAMA_RED else "black")
  }
  save_fig(p_fig3b, "Fig3B_Main_Forest", 183, 140)
}

# ============================================================================
# Fig 4: Subgroup Analysis Forest Plot (数据驱动)
# ============================================================================
cat("\n===== Fig 4: Subgroup Forest Plot =====\n")

subgroup <- safe_fread(file.path(TBL_DIR_V26, "v26_subgroup_results.csv"))
if (!is.null(subgroup) && nrow(subgroup) > 0) {
  subgroup[, y := .N:1]
  subgroup[, is_protective := HR < 1]
  subgroup[, label_en := Subgroup]

  p_fig4 <- ggplot(subgroup, aes(x = HR, y = y)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper),
                   height = 0.15, linewidth = 0.8) +
    geom_point(aes(color = is_protective), size = 3.5, shape = 22,
               fill = ifelse(subgroup$is_protective, FOREST_GREEN, JAMA_RED)) +
    scale_color_manual(values = c("TRUE" = FOREST_GREEN, "FALSE" = JAMA_RED), guide = "none") +
    scale_y_continuous(breaks = subgroup$y, labels = rev(subgroup$label_en)) +
    scale_x_log10(breaks = c(0.60, 0.70, 0.85, 1.00, 1.15, 1.30, 1.50),
                  labels = c("0.60", "0.70", "0.85", "1.00", "1.15", "1.30", "1.50")) +
    labs(
      title = "Figure 4. Subgroup Analyses (SEQTaRget ITT, bootstrap = 50)",
      subtitle = "18 pre-specified subgroups | 4 of 10 interaction tests FDR-significant (q < 0.10)",
      x = "Hazard Ratio (95% CI)", y = NULL,
      caption = "Higher severity subgroups consistently show protective or neutral estimates"
    ) +
    coord_cartesian(xlim = c(0.55, 1.65))
  save_fig(p_fig4, "Fig4_Subgroup_Forest", 183, 160)
}

# ============================================================================
# Fig 5: PO vs Non-PO Stratification Forest Plot
# ============================================================================
cat("\n===== Fig 5: PO vs Non-PO Stratification =====\n")

po_data <- data.table(
  stratum = c("Oral Phosphate (PO)", "Non-PO (IV or Conservative)", "Primary ITT (reference)"),
  hr = c(1.075, 1.053, 1.054),
  ci_lower = c(0.873, 0.931, 0.954),
  ci_upper = c(1.325, 1.190, 1.165),
  n = c("8,848", "23,445", "32,325"),
  events = c("1,329", "4,275", "5,610")
)
po_data$y <- 3:1
po_data$is_ref <- c(FALSE, FALSE, TRUE)

p_fig5 <- ggplot(po_data, aes(x = hr, y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.15, linewidth = 0.9, color = DEEP_NAVY) +
  geom_point(aes(color = is_ref), size = 4.5, shape = 15) +
  scale_color_manual(values = c("TRUE" = FOREST_GREEN, "FALSE" = JAMA_BLUE), guide = "none") +
  scale_y_continuous(breaks = 3:1, labels = po_data$stratum) +
  scale_x_log10(breaks = c(0.80, 0.90, 1.00, 1.10, 1.20, 1.40),
                labels = c("0.80", "0.90", "1.00", "1.10", "1.20", "1.40")) +
  labs(
    title = "Figure 5. PO vs Non-PO Phosphate Stratification",
    subtitle = "SEQTaRget ITT | bootstrap = 50 | DAG14 covariates",
    x = "Hazard Ratio (95% CI)", y = NULL
  ) +
  coord_cartesian(xlim = c(0.75, 1.50))
save_fig(p_fig5, "Fig5_PO_Stratification", 170, 90)

# ============================================================================
# Fig 6: Daily-Resolution HR (Days 1-14) (数据驱动)
# ============================================================================
cat("\n===== Fig 6: Daily-Resolution Hazard Ratios =====\n")

daily <- safe_fread(file.path(TBL_DIR, "v26_9_daily_hr_day1_14.csv"))
if (!is.null(daily) && nrow(daily) > 0) {
  daily[, sig := fcase(
    ci_lower > 1, "Significant harm",
    ci_upper < 1, "Significant benefit",
    default = "Non-significant"
  )]
  daily[, sig := factor(sig, levels = c("Significant harm", "Significant benefit", "Non-significant"))]

  p_fig6 <- ggplot(daily, aes(x = day, y = hr)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "grey85", alpha = 0.6) +
    geom_line(color = DEEP_NAVY, linewidth = 1.2) +
    geom_point(aes(color = sig), size = 3.5) +
    scale_color_manual(
      values = c("Significant harm" = JAMA_RED, "Significant benefit" = FOREST_GREEN,
                 "Non-significant" = "grey50"),
      name = NULL
    ) +
    scale_x_continuous(breaks = 1:14) +
    scale_y_log10(breaks = c(0.40, 0.50, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00),
                  labels = c("0.40", "0.50", "0.60", "0.80", "1.00", "1.25", "1.50", "1.75", "2.00")) +
    labs(
      title = "Figure 6. Daily-Resolution Hazard Ratios (Days 1\u201314)",
      subtitle = "Cox regression | DAG14 covariates | cluster-robust standard errors",
      x = "Follow-up Day", y = "Hazard Ratio (95% CI)",
      caption = "Exploratory analysis: 14 daily comparisons without multiplicity correction"
    )
  save_fig(p_fig6, "Fig6_Daily_HR", 183, 120)
}

# ============================================================================
# Fig 7: Piecewise Hazard Ratios (数据驱动)
# ============================================================================
cat("\n===== Fig 7: Piecewise Hazard Ratios =====\n")

seg <- safe_fread(file.path(TBL_DIR, "v26_9_secondary_time_segments.csv"))
if (!is.null(seg) && nrow(seg) > 0) {
  seg[, label := sprintf("Days %d\u2013%d", t_start, t_end)]
  seg[, y := .N:1]

  p_fig7 <- ggplot(seg, aes(x = estimate, y = y)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                   height = 0.15, linewidth = 1, color = DEEP_NAVY) +
    geom_point(size = 4.5, shape = 18, color = JAMA_BLUE) +
    scale_y_continuous(breaks = seg$y, labels = seg$label) +
    scale_x_log10(breaks = c(0.70, 0.85, 1.00, 1.15, 1.30, 1.50),
                  labels = c("0.70", "0.85", "1.00", "1.15", "1.30", "1.50")) +
    labs(
      title = "Figure 7. Piecewise Hazard Ratios by Follow-up Interval",
      subtitle = "Cox regression | DAG14 covariates | cluster-robust standard errors",
      x = "Hazard Ratio (95% CI)", y = NULL,
      caption = "3\u20137 day segment: HR = 1.308, p < 0.001 (strongest harmful effect)"
    ) +
    coord_cartesian(xlim = c(0.65, 1.65))

  for (i in seq_len(nrow(seg))) {
    p_fig7 <- p_fig7 + annotate("text", x = 1.70, y = i,
      label = sprintf("HR=%.3f (%.3f\u2013%.3f)", seg$estimate[i], seg$ci_lower[i], seg$ci_upper[i]),
      hjust = 0, size = 3.2)
  }
  save_fig(p_fig7, "Fig7_Piecewise_HR", 170, 110)
}

# ============================================================================
# Fig 8: Interaction Analysis Forest Plot (数据驱动)
# ============================================================================
cat("\n===== Fig 8: Interaction Forest Plot =====\n")

int <- safe_fread(file.path(TBL_DIR, "v26_9_interaction_fdr.csv"))
if (!is.null(int) && nrow(int) > 0) {
  int[, label := tools::toTitleCase(gsub("_", " ", var))]
  int[, label := factor(label, levels = rev(label))]
  int[, sig := fcase(
    significance_label == "FDR-significant", "FDR q < 0.10",
    raw_p < 0.05, "Raw p < 0.05",
    default = "Non-significant"
  )]

  p_fig8 <- ggplot(int, aes(x = hr, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                   height = 0.2, linewidth = 0.9) +
    geom_point(aes(color = sig), size = 4) +
    scale_color_manual(
      values = c("FDR q < 0.10" = JAMA_RED, "Raw p < 0.05" = EVAL_MOD, "Non-significant" = "grey50"),
      name = ""
    ) +
    scale_x_log10() +
    labs(
      title = "Figure 8. Interaction Analysis Forest Plot",
      subtitle = "Cox proportional hazards | DAG14 covariates | cluster-robust SE | BH-FDR correction",
      x = "Interaction Hazard Ratio (95% CI)", y = NULL,
      caption = "Interaction HR from Cox proportional hazards models with DAG14 covariates. BH-FDR correction applied across 10 interaction tests. None reached q < 0.10."
    )

  ann_x <- max(int$ci_upper) * 1.08
  for (i in seq_len(nrow(int))) {
    p_fig8 <- p_fig8 + annotate("text", x = ann_x, y = int$label[i],
      label = sprintf("q = %.3f", int$q_value[i]),
      hjust = 0, size = 3,
      color = if (int$sig[i] == "FDR q < 0.10") JAMA_RED else "grey40")
  }
  save_fig(p_fig8, "Fig8_Interaction_Forest", 183, 130)
}

# ============================================================================
# Fig 9: Secondary Outcomes Forest Plot (数据驱动)
# ============================================================================
cat("\n===== Fig 9: Secondary Outcomes Forest Plot =====\n")

sec_unified <- safe_fread(file.path(TBL_DIR, "v26_9_secondary_unified_results.csv"))
if (!is.null(sec_unified) && nrow(sec_unified) > 0) {
  # Filter to main outcomes (primary group)
  sec_main <- sec_unified[group == "primary"]
  sec_main[, y := .N:1]
  sec_main[, sig := fcase(
    p_value < 0.001, "p < 0.001",
    p_value < 0.01,  "p < 0.01",
    p_value < 0.05,  "p < 0.05",
    default = "Non-significant"
  )]
  # Shorten labels
  sec_main[, label_short := gsub("\\s*\\(.*\\)", "", outcome)]
  sec_main[, label_short := gsub("血磷达标.*", "Phosphate normalization", label_short)]
  sec_main[, label_short := gsub("输血需求.*", "Transfusion need", label_short)]
  sec_main[, label_short := gsub("AKI KDIGO.*", "AKI KDIGO \u22652", label_short)]
  sec_main[, label_short := gsub("心血管复合结局.*", "CV composite", label_short)]
  sec_main[, label_short := gsub("90天死亡率.*", "90-day mortality", label_short)]
  sec_main[, label_short := gsub("住院死亡率.*", "In-hospital mortality", label_short)]
  sec_main[, label_short := gsub("ICU死亡率.*", "ICU mortality", label_short)]

  p_fig9 <- ggplot(sec_main, aes(x = hr, y = y)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                   height = 0.18, linewidth = 0.9, color = DEEP_NAVY) +
    geom_point(aes(color = sig), size = 4, shape = 22) +
    scale_color_manual(
      values = c("p < 0.001" = JAMA_RED, "p < 0.01" = EVAL_MOD,
                 "p < 0.05" = EVAL_WEAK, "Non-significant" = "grey50"),
      name = ""
    ) +
    scale_y_continuous(breaks = sec_main$y, labels = rev(sec_main$label_short)) +
    labs(
      title = "Figure 9. Secondary Outcomes \u2014 SEQTaRget CCW ITT",
      subtitle = "DAG14 covariates | bootstrap = 25 | m = 5 imputations | Bonferroni alpha = 0.010",
      x = "Hazard Ratio (95% CI)", y = NULL,
      caption = "HR shown for binary outcomes. Phosphate normalization is a treatment process measure."
    )
  save_fig(p_fig9, "Fig9_Secondary_Outcomes", 183, 150)
}

# ============================================================================
# 补充图 S5: E-value Sensitivity Analysis
# ============================================================================
cat("\n===== Fig S5: E-value Curve =====\n")

evalue_data <- data.table(
  Analysis = c("Point Estimate\nE-value = 1.26", "CI Upper Limit\nE-value = 1.54"),
  E_value = c(1.26, 1.54),
  tier = c("Weak (E < 1.5)", "Moderate (E \u2265 1.5)")
)

p_s5 <- ggplot(evalue_data, aes(x = Analysis, y = E_value, fill = tier)) +
  geom_col(width = 0.45, alpha = 0.85) +
  geom_text(aes(label = sprintf("E = %.2f", E_value)), vjust = -0.5, size = 5) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = EVAL_MOD, linewidth = 0.7) +
  annotate("text", x = 2.5, y = 1.55, label = "Moderate threshold (E = 1.5)",
           hjust = 0, size = 3.2, color = EVAL_MOD) +
  annotate("text", x = 2.5, y = 1.05, label = "Null line (E = 1.0)",
           hjust = 0, size = 3.2, color = "grey40") +
  scale_fill_manual(
    values = c("Weak (E < 1.5)" = EVAL_WEAK, "Moderate (E \u2265 1.5)" = EVAL_MOD),
    name = "Robustness Tier"
  ) +
  labs(
    title = "Figure S5. E-value Sensitivity Analysis",
    subtitle = "Point estimate E = 1.26 (weak) | CI upper limit E = 1.54 (moderate)\nA moderately strong unmeasured confounder could explain the null result",
    x = NULL, y = "E-value"
  ) +
  coord_cartesian(ylim = c(0, 2.5))
save_fig(p_s5, "FigS5_Evalue_Curve", 170, 130)

# ============================================================================
# 补充图 S6: Sequential Trial Period Overview (数据驱动)
# ============================================================================
cat("\n===== Fig S6: Trial Period Overview =====\n")

per <- safe_fread(file.path(TBL_DIR, "v26_9_ccw_period_characteristics.csv"))
if (!is.null(per) && nrow(per) > 0) {
  per[, period_label := sprintf("Period %d", period_int)]
  per[, treatment_pct := n_tx / n_obs * 100]

  ps6a <- ggplot(per, aes(x = period_label, y = n_obs)) +
    geom_col(fill = JAMA_BLUE, alpha = 0.75, width = 0.6) +
    geom_text(aes(label = n_obs), vjust = -0.5, size = 3) +
    labs(title = "A. Enrollment per Sequential Period", x = NULL, y = "N") +
    theme(panel.grid.major.y = element_blank())

  ps6b <- ggplot(per, aes(x = period_label)) +
    geom_col(aes(y = n_events), fill = JAMA_RED, alpha = 0.75, width = 0.6) +
    geom_line(aes(y = event_rate / 100 * max(n_events), group = 1),
              color = DEEP_NAVY, linewidth = 1) +
    geom_point(aes(y = event_rate / 100 * max(n_events)),
               color = DEEP_NAVY, size = 2.5) +
    scale_y_continuous(
      name = "Events (n)",
      sec.axis = sec_axis(~ . / max(per$n_events) * 100, name = "Mortality Rate (%)")
    ) +
    labs(title = "B. Events and Mortality Rate by Period", x = NULL) +
    theme(panel.grid.major.y = element_blank())

  ps6c <- ggplot(per, aes(x = period_label, y = treatment_pct, group = 1)) +
    geom_line(color = FOREST_GREEN, linewidth = 1) +
    geom_point(color = FOREST_GREEN, size = 2.5) +
    geom_hline(yintercept = mean(per$treatment_pct), linetype = "dashed",
               color = "grey50", linewidth = 0.6) +
    labs(title = "C. Treatment Rate by Period", x = "Sequential Trial Period",
         y = "Treated (%)", caption = "Dashed line: global mean treatment rate") +
    theme(panel.grid.major.y = element_blank())

  ps6 <- ps6a / ps6b / ps6c +
    plot_layout(heights = c(1, 1, 1)) +
    plot_annotation(
      title = "Figure S6. Sequential Trial Period Overview",
      subtitle = sprintf("Total: %d trial-periods across %d periods | 6 analytical periods (periods 0\u20135)",
                         sum(per$n_obs), nrow(per)),
      theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0),
                    plot.subtitle = element_text(size = 9, color = "grey30", hjust = 0))
    )
  save_fig(ps6, "FigS6_Trial_Period_Overview", 183, 200)
}

# ============================================================================
# 补充图 S10: 血磷4层分层森林图 (NEW)
# ============================================================================
cat("\n===== Fig S10: Phosphate 4-Level Stratification (NEW) =====\n")

s10_data <- data.table(
  stratum = c("Severe (<0.32)", "Moderate (0.32\u20130.50)", "Mild (0.50\u20130.65)", "Borderline (0.65\u20130.80)"),
  hr = c(0.948, 0.760, 0.969, 0.968),
  ci_lower = c(0.509, 0.630, 0.879, 0.894),
  ci_upper = c(1.763, 0.916, 1.068, 1.049),
  n = c(200, 2652, 10136, 19342)
)
s10_data$y <- 4:1
s10_data$is_protective <- s10_data$hr < 1

p_s10 <- ggplot(s10_data, aes(x = hr, y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.15, linewidth = 0.9, color = DEEP_NAVY) +
  geom_point(aes(fill = is_protective), size = 4.5, shape = 22, color = "black") +
  scale_fill_manual(values = c("TRUE" = FOREST_GREEN, "FALSE" = JAMA_RED), guide = "none") +
  scale_y_continuous(breaks = 4:1, labels = s10_data$stratum) +
  scale_x_log10(breaks = c(0.50, 0.70, 0.85, 1.00, 1.20, 1.50, 2.00),
                labels = c("0.50", "0.70", "0.85", "1.00", "1.20", "1.50", "2.00")) +
  labs(
    title = "Figure S10. Phosphate Level Stratification (4-Level)",
    subtitle = "Conventional multivariable Cox regression | DAG14 covariates",
    x = "Hazard Ratio (95% CI)", y = NULL,
    caption = "Conventional Cox (not CCW) — estimates likely biased toward protection due to immortal time bias"
  ) +
  coord_cartesian(xlim = c(0.40, 2.00))
save_fig(p_s10, "FigS10_Phosphate_4Level", 170, 110)

# ============================================================================
# 补充图 S11: IV 补磷时机分布 (NEW)
# ============================================================================
cat("\n===== Fig S11: IV Phosphate Administration Timing (NEW) =====\n")

s11_data <- data.table(
  time_bin = c("0\u20133h", "3\u20136h", "6\u20139h", "9\u201312h", "12\u201318h", "18\u201324h", "24\u201348h"),
  n = c(1245, 1832, 1256, 738, 532, 268, 200),
  pct = c(20.5, 30.2, 20.7, 12.2, 8.8, 4.4, 3.3)
)
s11_data[, time_bin := factor(time_bin, levels = time_bin)]

p_s11 <- ggplot(s11_data, aes(x = time_bin, y = n)) +
  geom_col(fill = JAMA_BLUE, alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, pct)), vjust = -0.3, size = 3.2, lineheight = 0.9) +
  labs(
    title = "Figure S11. Distribution of Time to First IV Phosphate Administration",
    subtitle = sprintf("Median: 6.9 hours (IQR: 3.8\u201312.4) | N = %d IV-treated trial-periods", sum(s11_data$n)),
    x = "Time from T0 to First IV Phosphate", y = "Number of Trial-Periods"
  ) +
  theme(panel.grid.major.y = element_blank())
save_fig(p_s11, "FigS11_IV_Timing_Distribution", 183, 120)

# ============================================================================
# 补充图 S12: PO 给药途径分析 (NEW)
# ============================================================================
cat("\n===== Fig S12: PO Route Analysis (NEW) =====\n")

s12_data <- data.table(
  route = c("Oral (PO)", "Nasogastric (NG)", "Other"),
  hr = c(0.475, 1.0, NA),
  ci_lower = c(0.254, NA, NA),
  ci_upper = c(0.886, NA, NA),
  n = c(6651, 1047, 354)
)
s12_data <- s12_data[!is.na(hr)]
s12_data$y <- nrow(s12_data):1

p_s12 <- ggplot(s12_data, aes(x = hr, y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.15, linewidth = 0.9, color = DEEP_NAVY) +
  geom_point(size = 4.5, shape = 22, fill = FOREST_GREEN, color = "black") +
  scale_y_continuous(breaks = s12_data$y, labels = s12_data$route) +
  scale_x_log10(breaks = c(0.20, 0.30, 0.50, 0.75, 1.00, 1.50, 2.00),
                labels = c("0.20", "0.30", "0.50", "0.75", "1.00", "1.50", "2.00")) +
  labs(
    title = "Figure S12. Oral vs NG Phosphate Route Analysis",
    subtitle = "Conventional Cox regression | DAG14 covariates + PO cumulative dose",
    x = "Hazard Ratio (95% CI) vs NG Route", y = NULL,
    caption = "Oral route associated with lower mortality (HR=0.475, 95% CI 0.254\u20130.886). Confounded by illness severity."
  ) +
  coord_cartesian(xlim = c(0.15, 2.00))
save_fig(p_s12, "FigS12_PO_Route_Analysis", 170, 90)

# ============================================================================
# 补充图 S13: IV 磷酸盐剂量-反应 (NEW)
# ============================================================================
cat("\n===== Fig S13: IV Phosphate Dose-Response (NEW) =====\n")

s13_data <- data.table(
  group = c("Conservative\n(no phosphate)", "Oral only\n(PO/NG only)", "Low intensity\n(IV <15 mmol)",
            "Moderate intensity\n(IV 15\u201330 mmol)", "High intensity\n(IV >30 mmol)"),
  hr = c(1.00, 0.827, 0.710, 1.108, 0.802),
  ci_lower = c(NA, 0.769, 0.318, 0.654, 0.742),
  ci_upper = c(NA, 0.889, 1.584, 1.877, 0.867),
  n = c(12719, 11781, 65, 112, 7653)
)
s13_data$y <- 5:1

p_s13 <- ggplot(s13_data, aes(x = hr, y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.15, linewidth = 0.9, color = DEEP_NAVY) +
  geom_point(size = 4.5, shape = 22, fill = JAMA_BLUE, color = "black") +
  scale_y_continuous(breaks = 5:1, labels = s13_data$group) +
  scale_x_log10(breaks = c(0.30, 0.50, 0.70, 0.85, 1.00, 1.20, 1.50, 2.00),
                labels = c("0.30", "0.50", "0.70", "0.85", "1.00", "1.20", "1.50", "2.00")) +
  labs(
    title = "Figure S13. IV Phosphate Dose-Response by Exposure Intensity",
    subtitle = "Conventional multivariable Cox regression | DAG14 covariates",
    x = "Hazard Ratio (95% CI)", y = NULL,
    caption = "Conventional Cox (not CCW) \u2014 estimates likely biased toward protection. HR=1.001 per 10 mmol within IV-treated (p=0.889)."
  ) +
  coord_cartesian(xlim = c(0.25, 2.00))
save_fig(p_s13, "FigS13_Dose_Response", 170, 120)

# ============================================================================
# 补充图 S14: K vs Na 磷酸盐对比 (NEW)
# ============================================================================
cat("\n===== Fig S14: K vs Na Phosphate Comparison (NEW) =====\n")

s14_data <- data.table(
  drug = c("Potassium Phosphate", "Sodium Phosphate"),
  hr = c(1.00, 1.026),
  ci_lower = c(NA, 0.914),
  ci_upper = c(NA, 1.151),
  n = c(4732, 3098)
)
s14_data$y <- 2:1

p_s14 <- ggplot(s14_data, aes(x = hr, y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.15, linewidth = 0.9, color = DEEP_NAVY) +
  geom_point(size = 4.5, shape = 22, fill = JAMA_BLUE, color = "black") +
  scale_y_continuous(breaks = 2:1, labels = s14_data$drug) +
  scale_x_log10(breaks = c(0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20),
                labels = c("0.85", "0.90", "0.95", "1.00", "1.05", "1.10", "1.15", "1.20")) +
  labs(
    title = "Figure S14. Potassium vs Sodium Phosphate Comparison",
    subtitle = "Conventional Cox regression | adjusted for cumulative dose + DAG14 covariates",
    x = "Hazard Ratio (95% CI) vs Potassium Phosphate", y = NULL,
    caption = "No significant difference between formulations (p = 0.666). First-administered drug as proxy."
  ) +
  coord_cartesian(xlim = c(0.83, 1.22))
save_fig(p_s14, "FigS14_K_vs_Na_Phosphate", 170, 90)

# ============================================================================
# 汇总
# ============================================================================
cat("\n")
cat("=================================================================\n")
cat("  V26.9 Final Figure Generation Complete\n")
cat("=================================================================\n")
cat(sprintf("  Output directory: %s\n", FIG_DIR))
cat(sprintf("  Font: %s\n", FONT_FAMILY))
cat(sprintf("  DPI: 600 (PNG)\n"))
cat(sprintf("  Format: PDF (vector, cairo_pdf) + PNG (raster, 600dpi)\n"))
cat(sprintf("  Language: English only\n"))
cat("\n")
cat("  Main Figures (9):\n")
cat("    Fig2_Love_Plot.{pdf,png}\n")
cat("    Fig3B_Main_Forest.{pdf,png}\n")
cat("    Fig4_Subgroup_Forest.{pdf,png}\n")
cat("    Fig5_PO_Stratification.{pdf,png}\n")
cat("    Fig6_Daily_HR.{pdf,png}\n")
cat("    Fig7_Piecewise_HR.{pdf,png}\n")
cat("    Fig8_Interaction_Forest.{pdf,png}\n")
cat("    Fig9_Secondary_Outcomes.{pdf,png}\n")
cat("\n")
cat("  Supplementary Figures (5 NEW):\n")
cat("    FigS5_Evalue_Curve.{pdf,png}\n")
cat("    FigS6_Trial_Period_Overview.{pdf,png}\n")
cat("    FigS10_Phosphate_4Level.{pdf,png}  [NEW]\n")
cat("    FigS11_IV_Timing_Distribution.{pdf,png}  [NEW]\n")
cat("    FigS12_PO_Route_Analysis.{pdf,png}  [NEW]\n")
cat("    FigS13_Dose_Response.{pdf,png}  [NEW]\n")
cat("    FigS14_K_vs_Na_Phosphate.{pdf,png}  [NEW]\n")
cat("=================================================================\n")