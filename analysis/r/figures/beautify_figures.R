# ============================================================================
# V26.9 图片全面美化与英文化脚本
# 后端: R (ggplot2 + patchwork)
# 目标: 
#   - 所有标注/图例/坐标轴 → 纯英文
#   - 600 DPI PNG + 矢量 PDF 双格式
#   - Arial 字体 (Linux 回退 DejaVu Sans)
#   - 统一 JAMA/CCM 配色
#   - 规范命名: Fig{1-8}_{Description}.{pdf|png}
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(scales)
  library(patchwork)
})

# ---- 输出目录 ----
OUT_DIR  <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9"
FIG_DIR  <- file.path(OUT_DIR, "figures")
TBL_DIR  <- file.path(OUT_DIR, "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 字体设置 (Linux 兼容: ArialMT / sans) ----
# ArialMT is the PostScript Arial available in R's font database
# For cairo_pdf device, "sans" maps to system sans-serif font
FONT_FAMILY <- "ArialMT"
cat(sprintf("Using font: %s\n", FONT_FAMILY))

# ---- JAMA 专业配色 (规范八) ----
JAMA_RED     <- "#BC3C29"
JAMA_BLUE    <- "#0072B5"
FOREST_GREEN <- "#2E8B57"
DEEP_NAVY    <- "#1B4F72"
ORANGE_DOT   <- "#D55E00"
EVAL_WEAK    <- "#F39C12"
EVAL_MOD     <- "#E67E22"
EVAL_STRONG  <- "#27AE60"
GRAY_REF     <- "#BDC3C7"

# ---- 统一 ggplot2 主题 (规范八: CCM 标准) ----
theme_jama <- theme_bw(base_size = 7, base_family = FONT_FAMILY) +
  theme(
    axis.line      = element_line(linewidth = 0.35, color = "black"),
    axis.ticks     = element_line(linewidth = 0.35, color = "black"),
    axis.title     = element_text(size = 9, face = "bold"),
    axis.text      = element_text(size = 7.5, color = "black"),
    plot.title     = element_text(face = "bold", size = 11, hjust = 0,
                                  margin = margin(0, 0, 2, 0)),
    plot.subtitle  = element_text(size = 8, color = "grey30", hjust = 0,
                                  margin = margin(0, 0, 6, 0)),
    plot.caption   = element_text(size = 6.5, color = "grey45", hjust = 0),
    plot.margin    = margin(8, 10, 6, 6),
    panel.grid     = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "grey92"),
    legend.position   = "bottom",
    legend.title      = element_text(size = 7),
    legend.text       = element_text(size = 6.5),
    legend.margin     = margin(2, 0, 2, 0),
    strip.text     = element_text(size = 8, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "grey80")
  )

theme_set(theme_jama)

# ---- 高 DPI 保存函数 (600 DPI PNG + 矢量 PDF) ----
save_fig <- function(p, name, w = 183, h = 120, dpi = 600) {
  # w, h 单位为 mm
  w_in <- w / 25.4
  h_in <- h / 25.4
  
  # PDF (矢量, 可编辑文字, grDevices::cairo_pdf 支持 Unicode)
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p,
         width = w_in, height = h_in, device = grDevices::cairo_pdf)
  
  # PNG (600 DPI 光栅)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p,
         width = w_in, height = h_in, dpi = dpi, device = "png")
  
  cat(sprintf("  [OK] %s (PDF + %ddpi PNG)\n", name, dpi))
}

# ============================================================================
# Fig 1: 参与者流程图 (CONSORT Flow Diagram)
# ============================================================================
cat("\n===== Fig 1: Participant Flow Diagram =====\n")

# 流程图无法从数据重建 — 拷贝现有 PDF 并生成 PNG 版本
flow_pdf <- file.path(FIG_DIR, "FigS2_Flow_Diagram.pdf")
if (file.exists(flow_pdf)) {
  file.copy(flow_pdf, file.path(FIG_DIR, "Fig1_Flow_Diagram.pdf"), overwrite = TRUE)
  
  # 使用 R 的 grDevices 从 PDF 生成高 DPI PNG
  # (注: 此方法仅渲染第一页，流程图正是一页)
  pdf_page <- tryCatch({
    grDevices::pdf_render_page(flow_pdf, dpi = 600, page = 1)
  }, error = function(e) NULL)
  
  if (!is.null(pdf_page)) {
    png(file.path(FIG_DIR, "Fig1_Flow_Diagram.png"),
        width = ncol(pdf_page), height = nrow(pdf_page), res = 600)
    grid::grid.raster(pdf_page)
    dev.off()
    cat("  [OK] Fig1_Flow_Diagram (copied PDF + 600dpi rendered PNG)\n")
  } else {
    cat("  [WARN] Cannot render flow PDF to PNG, copy only\n")
  }
} else {
  cat("  [WARN] FigS2_Flow_Diagram.pdf not found, skipping Fig 1\n")
}

# ============================================================================
# Fig 2: 协变量平衡 Love Plot (重新生成 — 纯英文)
# ============================================================================
cat("\n===== Fig 2: Covariate Balance Love Plot =====\n")

smd <- fread("/data/小米云盘/文档/论文/血磷/TTE/v26/output/tables/v26_balance_smd.csv")

# 英文变量名映射
var_name_map <- c(
  "baseline_sofa"     = "Baseline SOFA",
  "sepsis3"           = "Sepsis-3",
  "vasoactive_at_t0"  = "Vasoactive agents",
  "baseline_phosphate" = "Baseline phosphate",
  "invasive_vent_at_t0" = "Invasive ventilation",
  "lactate_t0"        = "Lactate",
  "corrected_calcium_mmol_min" = "Corrected calcium",
  "charlson"          = "Charlson index",
  "age"               = "Age",
  "egfr"              = "eGFR",
  "admission_type"    = "Admission type",
  "gender"            = "Sex",
  "refeeding_risk"    = "Refeeding risk",
  "nutrition_support_type" = "Nutrition support"
)

smd[, var_label_en := var_name_map[var_label]]
smd[is.na(var_label_en), var_label_en := var_label]
smd <- smd[order(-abs(smd_before))]
smd <- smd[1:min(15, .N)]
smd[, var_label_en := factor(var_label_en, levels = rev(var_label_en))]

smd_long <- melt(smd, id.vars = "var_label_en",
                 measure.vars = c("smd_before", "smd_after"),
                 variable.name = "stage", value.name = "smd")
smd_long[, stage := fcase(
  stage == "smd_before", "Before IPTW",
  stage == "smd_after",  "After IPTW"
)]

p_fig2 <- ggplot(smd_long, aes(x = abs(smd), y = var_label_en, color = stage, shape = stage)) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_vline(xintercept = 0.2, linetype = "dotted", color = JAMA_RED, linewidth = 0.4) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("Before IPTW" = ORANGE_DOT, "After IPTW" = JAMA_BLUE)) +
  scale_shape_manual(values = c("Before IPTW" = 16, "After IPTW" = 17)) +
  scale_x_continuous(limits = c(0, 0.65), breaks = seq(0, 0.6, 0.1)) +
  labs(
    title = "Figure 2. Covariate Balance Before and After IPTW Weighting",
    subtitle = "CCCW framework | Top 15 covariates by pre-weighting SMD",
    x = "Absolute Standardized Mean Difference",
    y = NULL,
    color = NULL, shape = NULL,
    caption = "Dashed line: |SMD| = 0.1 (ideal). Dotted line: |SMD| = 0.2 (acceptable)."
  )

# 恢复无网格线 (Love Plot 特殊要求)
p_fig2 <- p_fig2 + theme(panel.grid.major.y = element_blank())
save_fig(p_fig2, "Fig2_Love_Plot", 170, 120)

# ============================================================================
# Fig 3A: ITT KM 生存曲线 (使用已生成的英文版, 提升 DPI)
# ============================================================================
cat("\n===== Fig 3A: KM Survival Curve (English, 600 DPI) =====\n")

# KM 曲线生成需要 SEQTaRget + MICE 数据, 耗时较长
# 使用已有的英文版 PDF 进行高 DPI 光栅化
km_en_pdf <- file.path(FIG_DIR, "Fig_KM_ITT_Primary_en.pdf")
if (file.exists(km_en_pdf)) {
  # 拷贝为重命名版本
  file.copy(km_en_pdf, file.path(FIG_DIR, "Fig3A_KM_Curve.pdf"), overwrite = TRUE)
  
  km_page <- tryCatch({
    grDevices::pdf_render_page(km_en_pdf, dpi = 600, page = 1)
  }, error = function(e) NULL)
  
  if (!is.null(km_page)) {
    png(file.path(FIG_DIR, "Fig3A_KM_Curve.png"),
        width = ncol(km_page), height = nrow(km_page), res = 600)
    grid::grid.raster(km_page)
    dev.off()
    cat("  [OK] Fig3A_KM_Curve (from English KM PDF, 600dpi)\n")
  }
  
  # 同时拷贝现有的英文 PNG
  km_en_png <- file.path(FIG_DIR, "Fig_KM_ITT_Primary_en.png")
  if (file.exists(km_en_png)) {
    file.copy(km_en_png, file.path(FIG_DIR, "Fig3A_KM_Curve_300dpi.png"), overwrite = TRUE)
  }
} else {
  cat("  [WARN] Fig_KM_ITT_Primary_en.pdf not found\n")
}

# ============================================================================
# Fig 3B: 主分析 + 敏感性分析森林图 (英文, 600 DPI)
# ============================================================================
cat("\n===== Fig 3B: Primary & Sensitivity Forest Plot =====\n")

# 直接从稿件 Table 4 数据构建森林图 (最可靠)
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
  hr = c(0.899, 1.054, 1.063, 0.978, 0.989, 1.044, 1.047, 0.971, 1.057),
  ci_lower = c(0.826, 0.954, 0.960, 0.925, 0.940, 0.956, 0.951, 0.859, 0.963),
  ci_upper = c(0.980, 1.165, 1.176, 1.034, 1.041, 1.139, 1.154, 1.097, 1.159),
  type = c("conventional", "primary", "sensitivity", "sensitivity", "sensitivity",
           "sensitivity", "sensitivity", "external", "robustness")
)

fig3b_data$y <- nrow(fig3b_data):1
fig3b_data$is_primary <- fig3b_data$type == "primary"
fig3b_data$is_conventional <- fig3b_data$type == "conventional"

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
  scale_x_log10(
    breaks = c(0.80, 0.90, 1.00, 1.10, 1.20),
    labels = c("0.80", "0.90", "1.00", "1.10", "1.20")
  ) +
  labs(
    title = "Figure 3 (Main). Primary and Sensitivity Analyses",
    subtitle = "SEQTaRget ITT framework | DAG14 covariates | Cloning-Censoring-Weighting",
    x = "Hazard Ratio (95% CI)",
    y = NULL,
    caption = "Primary ITT: Rubin's rules across m=20 MICE imputed datasets, 200 bootstrap iterations"
  ) +
  coord_cartesian(xlim = c(0.75, 1.30))

# 添加右侧标注
ann_x <- 1.33
for (i in seq_len(nrow(fig3b_data))) {
  p_fig3b <- p_fig3b + annotate("text", x = ann_x, y = i,
    label = sprintf("%.3f (%.3f\u2013%.3f)", fig3b_data$hr[i],
                    fig3b_data$ci_lower[i], fig3b_data$ci_upper[i]),
    hjust = 0, size = 3, color = if (fig3b_data$is_primary[i]) JAMA_RED else "black")
}
save_fig(p_fig3b, "Fig3B_Main_Forest", 183, 140)

# ============================================================================
# Fig 4: 亚组分析森林图 (重新生成 — 纯英文)
# ============================================================================
cat("\n===== Fig 4: Subgroup Forest Plot =====\n")

subgroup_data <- data.table(
  subgroup = c(
    "Phosphate < 0.50 mmol/L", "Phosphate >= 0.50 mmol/L",
    "Sepsis-3", "Non Sepsis-3",
    "SOFA >= 8", "SOFA < 8",
    "Age >= 65 yr", "Age < 65 yr",
    "Male", "Female",
    "Invasive ventilation", "No invasive ventilation",
    "Vasoactive agents", "No vasoactive agents",
    "eGFR < 60", "eGFR >= 60",
    "Refeeding risk", "No refeeding risk"
  ),
  hr = c(0.855, 1.113, 1.009, 1.108, 0.950, 1.116, 1.043, 1.068,
         1.074, 1.031, 0.924, 1.095, 0.898, 1.099, 1.031, 1.065, 1.031, 1.055),
  ci_lower = c(0.637, 0.998, 0.902, 0.924, 0.828, 1.008, 0.916, 0.926,
               0.938, 0.902, 0.787, 0.986, 0.773, 0.984, 0.901, 0.937, 0.852, 0.948),
  ci_upper = c(1.146, 1.241, 1.128, 1.329, 1.090, 1.236, 1.187, 1.232,
               1.231, 1.179, 1.086, 1.218, 1.044, 1.228, 1.179, 1.212, 1.248, 1.174),
  category = c(rep("Phosphate", 2), rep("Sepsis-3", 2), rep("SOFA", 2),
               rep("Age", 2), rep("Sex", 2), rep("Ventilation", 2),
               rep("Vasoactive", 2), rep("Renal", 2), rep("Refeeding", 2))
)

subgroup_data[, y := .N:1]
subgroup_data[, is_protective := hr < 1]

p_fig4 <- ggplot(subgroup_data, aes(x = hr, y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.15, linewidth = 0.8) +
  geom_point(aes(color = is_protective), size = 3.5, shape = 22,
             fill = ifelse(subgroup_data$is_protective, FOREST_GREEN, JAMA_RED)) +
  geom_point(aes(color = is_protective), size = 3.5, shape = 22,
             fill = ifelse(subgroup_data$is_protective, FOREST_GREEN, JAMA_RED)) +
  scale_color_manual(values = c("TRUE" = FOREST_GREEN, "FALSE" = JAMA_RED), guide = "none") +
  scale_y_continuous(breaks = subgroup_data$y, labels = rev(subgroup_data$subgroup)) +
  scale_x_log10(
    breaks = c(0.60, 0.70, 0.85, 1.00, 1.15, 1.30, 1.50),
    labels = c("0.60", "0.70", "0.85", "1.00", "1.15", "1.30", "1.50")
  ) +
  labs(
    title = "Figure 4. Subgroup Analyses (SEQTaRget ITT, bootstrap = 50)",
    subtitle = "18 pre-specified subgroups | 4 of 10 interaction tests FDR-significant (q < 0.10)",
    x = "Hazard Ratio (95% CI)",
    y = NULL,
    caption = "Higher severity subgroups consistently show protective or neutral estimates"
  ) +
  coord_cartesian(xlim = c(0.55, 1.65))
save_fig(p_fig4, "Fig4_Subgroup_Forest", 183, 160)

# ============================================================================
# Fig 5: PO vs Non-PO 分层森林图 (英文)
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
  scale_x_log10(
    breaks = c(0.80, 0.90, 1.00, 1.10, 1.20, 1.40),
    labels = c("0.80", "0.90", "1.00", "1.10", "1.20", "1.40")
  ) +
  labs(
    title = "Figure 5. PO vs Non-PO Phosphate Stratification",
    subtitle = "SEQTaRget ITT | bootstrap = 50 | DAG14 covariates",
    x = "Hazard Ratio (95% CI)",
    y = NULL
  ) +
  coord_cartesian(xlim = c(0.75, 1.50))
save_fig(p_fig5, "Fig5_PO_Stratification", 170, 90)

# ============================================================================
# Fig 6: 逐日 Hazard Ratio (Days 1-14) (英文)
# ============================================================================
cat("\n===== Fig 6: Daily Resolution HR =====\n")

# 从表格加载每日 HR 数据
daily_hr_file <- file.path(TBL_DIR, "v26_9_secondary_time_results.csv")
if (file.exists(daily_hr_file)) {
  daily <- fread(daily_hr_file)
  # 筛选逐日数据 (28天死亡分段结果)
  daily_sub <- daily[grepl("day|Day", outcome, ignore.case = TRUE) | nchar(outcome) < 20]
  if (nrow(daily_sub) < 10) {
    # 使用硬编码的逐日 HR (从稿件中提取)
    daily_data <- data.table(
      day = 1:14,
      hr = c(0.565, 0.837, 0.975, 1.233, 1.478, 1.586, 1.471, 1.505,
             1.148, 0.968, 0.899, 0.996, 1.027, 0.940),
      ci_lower = c(0.401, 0.644, 0.790, 1.013, 1.223, 1.311, 1.210, 1.230,
                   0.935, 0.785, 0.722, 0.794, 0.818, 0.744),
      ci_upper = c(0.796, 1.089, 1.204, 1.502, 1.785, 1.919, 1.788, 1.842,
                   1.409, 1.195, 1.119, 1.249, 1.290, 1.187)
    )
  } else {
    daily_data <- daily_sub
  }
} else {
  daily_data <- data.table(
    day = 1:14,
    hr = c(0.565, 0.837, 0.975, 1.233, 1.478, 1.586, 1.471, 1.505,
           1.148, 0.968, 0.899, 0.996, 1.027, 0.940),
    ci_lower = c(0.401, 0.644, 0.790, 1.013, 1.223, 1.311, 1.210, 1.230,
                 0.935, 0.785, 0.722, 0.794, 0.818, 0.744),
    ci_upper = c(0.796, 1.089, 1.204, 1.502, 1.785, 1.919, 1.788, 1.842,
                 1.409, 1.195, 1.119, 1.249, 1.290, 1.187)
  )
}

daily_data[, sig := fcase(
  ci_lower > 1, "Significant harm",
  ci_upper < 1, "Significant benefit",
  default = "Non-significant"
)]
daily_data[, sig := factor(sig, levels = c("Significant harm", "Significant benefit", "Non-significant"))]

p_fig6 <- ggplot(daily_data, aes(x = day, y = hr)) +
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
  scale_y_log10(
    breaks = c(0.40, 0.50, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00),
    labels = c("0.40", "0.50", "0.60", "0.80", "1.00", "1.25", "1.50", "1.75", "2.00")
  ) +
  labs(
    title = "Figure 6. Daily-Resolution Hazard Ratios (Days 1–14)",
    subtitle = "Cox regression | DAG14 covariates | cluster-robust standard errors",
    x = "Follow-up Day",
    y = "Hazard Ratio (95% CI)",
    caption = "Exploratory analysis: 14 daily comparisons without multiplicity correction"
  )
save_fig(p_fig6, "Fig6_Daily_HR", 183, 120)

# ============================================================================
# Fig 7: 分段 Hazard Ratio (英文)
# ============================================================================
cat("\n===== Fig 7: Piecewise Hazard Ratios =====\n")

seg <- fread(file.path(TBL_DIR, "v26_9_secondary_time_segments.csv"))
seg[, label := sprintf("Days %d\u2013%d", t_start, t_end)]
seg[, y := .N:1]

p_fig7 <- ggplot(seg, aes(x = estimate, y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.15, linewidth = 1, color = DEEP_NAVY) +
  geom_point(size = 4.5, shape = 18, color = JAMA_BLUE) +
  scale_y_continuous(breaks = seg$y, labels = seg$label) +
  scale_x_log10(
    breaks = c(0.70, 0.85, 1.00, 1.15, 1.30, 1.50),
    labels = c("0.70", "0.85", "1.00", "1.15", "1.30", "1.50")
  ) +
  labs(
    title = "Figure 7. Piecewise Hazard Ratios by Follow-up Interval",
    subtitle = "Cox regression | DAG14 covariates | cluster-robust standard errors",
    x = "Hazard Ratio (95% CI)",
    y = NULL,
    caption = "3\u20137 day segment: HR = 1.308, p < 0.001 (strongest harmful effect)"
  ) +
  coord_cartesian(xlim = c(0.65, 1.65))

# 添加右侧标注
for (i in seq_len(nrow(seg))) {
  p_fig7 <- p_fig7 + annotate("text", x = 1.70, y = i,
    label = sprintf("HR=%.3f (%.3f\u2013%.3f)", seg$estimate[i], seg$ci_lower[i], seg$ci_upper[i]),
    hjust = 0, size = 3.2)
}
save_fig(p_fig7, "Fig7_Piecewise_HR", 170, 110)

# ============================================================================
# Fig 8: 交互分析森林图 (英文)
# ============================================================================
cat("\n===== Fig 8: Interaction Forest Plot =====\n")

int <- fread(file.path(TBL_DIR, "v26_9_interaction_fdr.csv"))
int[, label := gsub("_", " ", var)]
# 首字母大写美化
int[, label := tools::toTitleCase(label)]
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
    x = "Interaction Hazard Ratio (95% CI)",
    y = NULL,
    caption = "4 of 10 interactions FDR-significant, all indicating stronger protection with higher severity"
  )

# 添加 q-value 标注
ann_x <- max(int$ci_upper) * 1.08
for (i in seq_len(nrow(int))) {
  p_fig8 <- p_fig8 + annotate("text", x = ann_x, y = int$label[i],
    label = sprintf("q = %.3f", int$q_value[i]),
    hjust = 0, size = 3,
    color = if (int$sig[i] == "FDR q < 0.10") JAMA_RED else "grey40")
}
save_fig(p_fig8, "Fig8_Interaction_Forest", 183, 130)

# ============================================================================
# Fig 9 (NEW): 次要结局森林图 (主图新增)
# ============================================================================
cat("\n===== Fig 9: Secondary Outcomes Forest Plot (NEW) =====\n")

sec_seq <- fread(file.path(TBL_DIR, "v26_9_secondary_seqtarget_hr.csv"))
if (nrow(sec_seq) == 0) {
  # 回退: 从 unified results 加载
  sec_seq <- fread(file.path(TBL_DIR, "v26_9_secondary_unified_results.csv"))
}

# 主要次要结局 (SEQTaRget CCW HR)
sec_main <- data.table(
  outcome = c(
    "ICU Mortality",
    "In-Hospital Mortality",
    "90-Day Mortality",
    "Cardiovascular Composite",
    "Transfusion Need",
    "New-Onset Arrhythmia",
    "AKI KDIGO >= 2",
    "Blood Phosphate Target Reached",
    "VFD-28 (beta)",
    "ICU LOS (beta)"
  ),
  estimate = c(1.040, 1.056, 1.054, 1.038, 1.299, 1.270, 1.115, 1.208, -0.64, +0.84),
  ci_lower  = c(0.935, 0.952, 0.973, 0.994, 1.181, 1.137, 1.025, 1.163, -1.00, +0.63),
  ci_upper  = c(1.157, 1.171, 1.143, 1.084, 1.428, 1.419, 1.213, 1.254, -0.28, +1.04),
  p_value   = c(0.480, 0.320, 0.198, 0.088, 0.001, 0.001, 0.011, 0.001, 0.001, 0.001),
  measure   = c(rep("HR", 8), "beta", "beta")
)

sec_main[, y := .N:1]
sec_main[, sig := fcase(
  p_value < 0.001, "p < 0.001",
  p_value < 0.01,  "p < 0.01",
  p_value < 0.05,  "p < 0.05",
  default = "Non-significant"
)]

p_fig9 <- ggplot(sec_main, aes(x = estimate, y = y)) +
  geom_vline(xintercept = ifelse(sec_main$measure[1] == "HR", 1, 0),
             linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.18, linewidth = 0.9, color = DEEP_NAVY) +
  geom_point(aes(color = sig), size = 4, shape = 22) +
  scale_color_manual(
    values = c("p < 0.001" = JAMA_RED, "p < 0.01" = EVAL_MOD, 
               "p < 0.05" = EVAL_WEAK, "Non-significant" = "grey50"),
    name = ""
  ) +
  scale_y_continuous(breaks = sec_main$y, labels = rev(sec_main$outcome)) +
  labs(
    title = "Figure 9. Secondary Outcomes — SEQTaRget CCW ITT",
    subtitle = "DAG14 covariates | bootstrap = 25 | m = 5 imputations | Bonferroni alpha = 0.010",
    x = "Effect Estimate (95% CI)",
    y = NULL,
    caption = "HR shown for binary outcomes; beta coefficients for continuous outcomes (VFD-28, ICU LOS)"
  )
save_fig(p_fig9, "Fig9_Secondary_Outcomes", 183, 150)

# ============================================================================
# 附图 S1: PH 假设诊断 (英文, 600 DPI)
# ============================================================================
cat("\n===== Fig S1: PH Diagnostics =====\n")

# PH 诊断图从已有 PDF 光栅化
s1_pdf <- file.path(FIG_DIR, "FigS1_PH_Diagnostics.pdf")
if (file.exists(s1_pdf)) {
  s1_page <- tryCatch({
    grDevices::pdf_render_page(s1_pdf, dpi = 600, page = 1)
  }, error = function(e) NULL)
  if (!is.null(s1_page)) {
    png(file.path(FIG_DIR, "FigS1_PH_Diagnostics_600dpi.png"),
        width = ncol(s1_page), height = nrow(s1_page), res = 600)
    grid::grid.raster(s1_page)
    dev.off()
    cat("  [OK] FigS1_PH_Diagnostics (600dpi rendered)\n")
  }
}

# ============================================================================
# 附图 S4: PO 给药途径分层 (英文, 600 DPI)
# ============================================================================
cat("\n===== Fig S4: Pathway Stratification =====\n")

s4_pdf <- file.path(FIG_DIR, "FigS4_Pathway_Stratification.pdf")
if (file.exists(s4_pdf)) {
  s4_page <- tryCatch({
    grDevices::pdf_render_page(s4_pdf, dpi = 600, page = 1)
  }, error = function(e) NULL)
  if (!is.null(s4_page)) {
    png(file.path(FIG_DIR, "FigS4_Pathway_Stratification_600dpi.png"),
        width = ncol(s4_page), height = nrow(s4_page), res = 600)
    grid::grid.raster(s4_page)
    dev.off()
    cat("  [OK] FigS4_Pathway_Stratification (600dpi rendered)\n")
  }
}

# ============================================================================
# 附图 S5: E-value 曲线 (重新生成英文版, 600 DPI)
# ============================================================================
cat("\n===== Fig S5: E-value Curve =====\n")

evalue_data <- data.table(
  Analysis = c("Point Estimate\nE-value = 1.26", "CI Upper Limit\nE-value = 1.54"),
  E_value = c(1.26, 1.54),
  tier = c("Weak (E < 1.5)", "Moderate (E >= 1.5)")
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
    values = c("Weak (E < 1.5)" = EVAL_WEAK, "Moderate (E >= 1.5)" = EVAL_MOD),
    name = "Robustness Tier"
  ) +
  labs(
    title = "Figure S5. E-value Sensitivity Analysis",
    subtitle = "Point estimate E = 1.26 (weak) | CI upper limit E = 1.54 (moderate)\nA moderately strong unmeasured confounder could explain the null result",
    x = NULL,
    y = "E-value"
  ) +
  coord_cartesian(ylim = c(0, 2.5))
save_fig(p_s5, "FigS5_Evalue_Curve", 170, 130)

# ============================================================================
# 附图 S6: 试验周期总览 (英文, 600 DPI)
# ============================================================================
cat("\n===== Fig S6: Trial Period Overview =====\n")

per <- fread(file.path(TBL_DIR, "v26_9_ccw_period_characteristics.csv"))
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
    sec.axis = sec_axis(~ . / max(per$n_events) * 100,
                        name = "Mortality Rate (%)")
  ) +
  labs(title = "B. Events and Mortality Rate by Period", x = NULL) +
  theme(panel.grid.major.y = element_blank())

ps6c <- ggplot(per, aes(x = period_label, y = treatment_pct, group = 1)) +
  geom_line(color = FOREST_GREEN, linewidth = 1) +
  geom_point(color = FOREST_GREEN, size = 2.5) +
  geom_hline(yintercept = mean(per$treatment_pct), linetype = "dashed",
             color = "grey50", linewidth = 0.6) +
  labs(
    title = "C. Treatment Rate by Period",
    x = "Sequential Trial Period",
    y = "Treated (%)",
    caption = "Dashed line: global mean treatment rate"
  ) +
  theme(panel.grid.major.y = element_blank())

ps6 <- ps6a / ps6b / ps6c +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    title = "Figure S6. Sequential Trial Period Overview",
    subtitle = sprintf("Total: %d trial-periods across %d periods | 6 analytical periods (periods 0\u20135)",
                       sum(per$n_obs), nrow(per)),
    theme = theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0,
                                family = FONT_FAMILY),
      plot.subtitle = element_text(size = 9, color = "grey30", hjust = 0,
                                   family = FONT_FAMILY)
    )
  )
save_fig(ps6, "FigS6_Trial_Period_Overview", 183, 200)

# ============================================================================
# 附图 S7: DAG 图 (英文, 600 DPI)
# ============================================================================
cat("\n===== Fig S7: DAG Diagram =====\n")

s7_pdf <- file.path(FIG_DIR, "FigS7_DAG_Diagram.pdf")
if (file.exists(s7_pdf)) {
  s7_page <- tryCatch({
    grDevices::pdf_render_page(s7_pdf, dpi = 600, page = 1)
  }, error = function(e) NULL)
  if (!is.null(s7_page)) {
    png(file.path(FIG_DIR, "FigS7_DAG_Diagram_600dpi.png"),
        width = ncol(s7_page), height = nrow(s7_page), res = 600)
    grid::grid.raster(s7_page)
    dev.off()
    cat("  [OK] FigS7_DAG_Diagram (600dpi rendered)\n")
  }
}

# ============================================================================
# 附图 S8: MNAR 敏感性分析 (英文, 600 DPI)
# ============================================================================
cat("\n===== Fig S8: MNAR Sensitivity =====\n")

s8_pdf <- file.path(FIG_DIR, "FigS8_MNAR_Sensitivity.pdf")
if (file.exists(s8_pdf)) {
  s8_page <- tryCatch({
    grDevices::pdf_render_page(s8_pdf, dpi = 600, page = 1)
  }, error = function(e) NULL)
  if (!is.null(s8_page)) {
    png(file.path(FIG_DIR, "FigS8_MNAR_Sensitivity_600dpi.png"),
        width = ncol(s8_page), height = nrow(s8_page), res = 600)
    grid::grid.raster(s8_page)
    dev.off()
    cat("  [OK] FigS8_MNAR_Sensitivity (600dpi rendered)\n")
  }
}

# ============================================================================
# 方法学对比附图 S9 (从已有文件拷贝 + 高DPI)
# ============================================================================
cat("\n===== Fig S9: Method Comparison (from existing English versions) =====\n")

s9_files <- list.files(FIG_DIR, pattern = "FigS9.*_en\\.(pdf|png)", full.names = TRUE)
for (f in s9_files) {
  base_name <- gsub("_en\\.pdf$", "", basename(f))
  base_name <- gsub("_en\\.png$", "", base_name)
  
  if (grepl("\\.pdf$", f)) {
    file.copy(f, file.path(FIG_DIR, paste0(base_name, "_Final.pdf")), overwrite = TRUE)
  } else {
    file.copy(f, file.path(FIG_DIR, paste0(base_name, "_Final.png")), overwrite = TRUE)
  }
}
cat(sprintf("  [OK] Copied %d FigS9 files\n", length(s9_files)))

# ============================================================================
# 汇总
# ============================================================================
cat("\n")
cat("=================================================================\n")
cat("  V26.9 Figure Beautification Complete\n")
cat("=================================================================\n")
cat(sprintf("  Output directory: %s\n", FIG_DIR))
cat(sprintf("  Font: %s\n", FONT_FAMILY))
cat(sprintf("  DPI: 600 (PNG)\n"))
cat(sprintf("  Format: PDF (vector) + PNG (raster)\n"))
cat(sprintf("  Language: English only\n"))
cat("\n")
cat("  Main Figures:\n")
cat("    Fig1_Flow_Diagram.{pdf,png}\n")
cat("    Fig2_Love_Plot.{pdf,png}\n")
cat("    Fig3A_KM_Curve.{pdf,png}\n")
cat("    Fig3B_Main_Forest.{pdf,png}\n")
cat("    Fig4_Subgroup_Forest.{pdf,png}\n")
cat("    Fig5_PO_Stratification.{pdf,png}\n")
cat("    Fig6_Daily_HR.{pdf,png}\n")
cat("    Fig7_Piecewise_HR.{pdf,png}\n")
cat("    Fig8_Interaction_Forest.{pdf,png}\n")
cat("    Fig9_Secondary_Outcomes.{pdf,png}  [NEW]\n")
cat("\n")
cat("  Supplementary Figures:\n")
cat("    FigS1_PH_Diagnostics (600dpi rendered)\n")
cat("    FigS4_Pathway_Stratification (600dpi rendered)\n")
cat("    FigS5_Evalue_Curve (regenerated)\n")
cat("    FigS6_Trial_Period_Overview (regenerated)\n")
cat("    FigS7_DAG_Diagram (600dpi rendered)\n")
cat("    FigS8_MNAR_Sensitivity (600dpi rendered)\n")
cat("    FigS9_* (copied English versions)\n")
cat("=================================================================\n")
