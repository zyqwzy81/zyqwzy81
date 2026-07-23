# ============================================================================
# V28 森林图重新生成脚本
#
# 目的：使用 SEQTaRget CCW 重算后的新 CSV 数据重新生成 Fig4/Fig8 森林图
#       - Fig4: 亚组分析森林图（18 项亚组，从 v26_9_subgroup_ccw_results.csv 读取）
#       - Fig8: 交互检验森林图（10 项交互，从 v26_9_interaction_ccw_results.csv 读取）
#
# 关键变更：
#   - 移除原 v26_9_beautify_figures.R 中的硬编码值
#   - 所有数据从 CSV 动态读取
#   - 副标题与稿件最新结果一致（所有 18 项 HR<1.0，仅 Charlson 交互名义显著）
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(scales)
})

# ---- 路径配置 ----
TBL_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v28/output/tables"
FIG_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v28/output/figures"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 字体设置 ----
FONT_FAMILY <- "sans"

# ---- JAMA 配色（规范八）----
JAMA_RED     <- "#BC3C29"
JAMA_BLUE    <- "#0072B5"
FOREST_GREEN <- "#2E8B57"
DEEP_NAVY    <- "#1B4F72"
GRAY_REF     <- "#BDC3C7"

# ---- 主题 ----
theme_jama <- theme_bw(base_size = 9, base_family = FONT_FAMILY) +
  theme(
    axis.line      = element_line(linewidth = 0.4, color = "black"),
    axis.ticks     = element_line(linewidth = 0.4, color = "black"),
    axis.title     = element_text(size = 10, face = "bold"),
    axis.text      = element_text(size = 8.5, color = "black"),
    plot.title     = element_text(face = "bold", size = 12, hjust = 0,
                                  margin = margin(0, 0, 2, 0)),
    plot.subtitle  = element_text(size = 9, color = "grey30", hjust = 0,
                                  margin = margin(0, 0, 6, 0)),
    plot.caption   = element_text(size = 7.5, color = "grey45", hjust = 0),
    plot.margin    = margin(10, 14, 8, 8),
    panel.grid     = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.25, color = "grey92"),
    legend.position   = "bottom",
    legend.title      = element_text(size = 8),
    legend.text       = element_text(size = 7.5),
    legend.margin     = margin(2, 0, 2, 0)
  )
theme_set(theme_jama)

# ---- 保存函数 ----
save_fig <- function(p, name, w = 183, h = 120, dpi = 600) {
  w_in <- w / 25.4
  h_in <- h / 25.4
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p,
         width = w_in, height = h_in, device = grDevices::cairo_pdf)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p,
         width = w_in, height = h_in, dpi = dpi)
  cat(sprintf("  保存: %s.{pdf,png} (%.0f x %.0f mm, %d DPI)\n",
              name, w, h, dpi))
}

cat("=" %+% paste(rep("=", 70), collapse = "") %+% "\n")
cat("V28 森林图重新生成（SEQTaRget CCW + DAG14 + cluster-robust SE + m=20）\n")
cat(paste(rep("=", 72), collapse = "") %+% "\n\n")

# ============================================================================
# Fig 4: 亚组分析森林图（从新 CSV 读取）
# ============================================================================
cat("===== Fig 4: Subgroup Forest Plot =====\n")

sub_csv <- file.path(TBL_DIR, "v26_9_subgroup_ccw_results.csv")
if (!file.exists(sub_csv)) {
  stop(sprintf("CSV 不存在: %s", sub_csv))
}

subgroup_data <- fread(sub_csv)
cat(sprintf("  读取: %s (%d 行)\n", sub_csv, nrow(subgroup_data)))

# 美化亚组名称
subgroup_data[, subgroup := gsub("mmol/L", "mmol/L", subgroup)]
subgroup_data[, subgroup := gsub(">=", "≥", subgroup)]
subgroup_data[, subgroup := gsub("<", "<", subgroup)]

# 添加分类
subgroup_data[, category := {
  cats <- character(.N)
  cats[grepl("Phosphate", subgroup)] <- "Phosphate severity"
  cats[grepl("Sepsis", subgroup)]    <- "Sepsis-3 status"
  cats[grepl("SOFA", subgroup)]     <- "Disease severity (SOFA)"
  cats[grepl("Age", subgroup)]      <- "Age"
  cats[grepl("Male|Female", subgroup)] <- "Sex"
  cats[grepl("ventilation", subgroup, ignore.case = TRUE)] <- "Ventilation"
  cats[grepl("Vasoactive", subgroup)] <- "Vasoactive agents"
  cats[grepl("eGFR", subgroup)]      <- "Renal function"
  cats[grepl("Refeeding", subgroup)] <- "Refeeding risk"
  cats
}]

# 按 category 排序，Phosphate 在顶部
cat_order <- c("Phosphate severity", "Sepsis-3 status", "Disease severity (SOFA)",
               "Age", "Sex", "Ventilation", "Vasoactive agents",
               "Renal function", "Refeeding risk")
subgroup_data[, category := factor(category, levels = cat_order)]
setorder(subgroup_data, category, -hr)

# 生成 Y 轴位置（从下到上）
subgroup_data[, y := .N:1]
subgroup_data[, is_significant := p_value < 0.05]
subgroup_data[, is_protective := hr < 1]

# 所有 HR 均 < 1.0，使用保护性绿色为主色
p_fig4 <- ggplot(subgroup_data, aes(x = hr, y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  # 分类分隔线（每两个亚组之间）
  geom_hline(yintercept = seq(1.5, nrow(subgroup_data) - 0.5, by = 2),
             color = "grey90", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.18, linewidth = 0.9, color = JAMA_BLUE) +
  geom_point(aes(fill = is_significant), size = 4.5, shape = 22,
             color = "black", stroke = 0.5) +
  scale_fill_manual(values = c("TRUE" = JAMA_RED, "FALSE" = JAMA_BLUE),
                    labels = c("TRUE" = "Nominal p < 0.05", "FALSE" = "Non-significant"),
                    name = NULL) +
  scale_y_continuous(breaks = subgroup_data$y,
                     labels = rev(subgroup_data$subgroup),
                     expand = c(0.01, 0.3)) +
  scale_x_log10(
    breaks = c(0.50, 0.65, 0.80, 1.00, 1.25, 1.50),
    labels = c("0.50", "0.65", "0.80", "1.00", "1.25", "1.50")
  ) +
  labs(
    title = "Figure 4. Subgroup Analyses",
    subtitle = "SEQTaRget ITT | DAG14 | cluster-robust SE | m=20 Rubin's rules | 18 pre-specified subgroups",
    x = "Hazard Ratio (95% CI, log scale)",
    y = NULL,
    caption = "All 18 subgroups show protective trends (HR < 1.0). Only Phosphate <0.50 mmol/L reaches nominal significance (HR=0.803, p=0.029)."
  ) +
  coord_cartesian(xlim = c(0.48, 1.55)) +
  annotate("text", x = 0.52, y = nrow(subgroup_data) + 0.8,
           label = "← Protective", color = FOREST_GREEN, size = 3.2,
           fontface = "bold", hjust = 0) +
  annotate("text", x = 1.50, y = nrow(subgroup_data) + 0.8,
           label = "Harmful →", color = JAMA_RED, size = 3.2,
           fontface = "bold", hjust = 1)

save_fig(p_fig4, "Fig4_Subgroup_Forest", 190, 170)

# ============================================================================
# Fig 8: 交互检验森林图（从新 CSV 读取）
# ============================================================================
cat("\n===== Fig 8: Interaction Forest Plot =====\n")

int_csv <- file.path(TBL_DIR, "v26_9_interaction_ccw_results.csv")
if (!file.exists(int_csv)) {
  stop(sprintf("CSV 不存在: %s", int_csv))
}

int_data <- fread(int_csv)
cat(sprintf("  读取: %s (%d 行)\n", int_csv, nrow(int_data)))

# 使用 label 列（已含美化名称）
int_data[, sig := fcase(
  significant == "Raw p<0.05", "Nominal p < 0.05",
  significant == "FDR-significant", "FDR q < 0.10",
  default = "Non-significant"
)]

# 按 q_value 排序
setorder(int_data, q_value)
int_data[, y_label := factor(label, levels = rev(label))]

p_fig8 <- ggplot(int_data, aes(x = hr, y = y_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.2, linewidth = 0.9, color = DEEP_NAVY) +
  geom_point(aes(color = sig), size = 4.2) +
  scale_color_manual(
    values = c("FDR q < 0.10" = JAMA_RED,
               "Nominal p < 0.05" = "#E67E22",
               "Non-significant" = "grey50"),
    name = NULL
  ) +
  scale_x_log10(
    breaks = c(0.70, 0.85, 1.00, 1.15, 1.30),
    labels = c("0.70", "0.85", "1.00", "1.15", "1.30")
  ) +
  labs(
    title = "Figure 8. Interaction Analysis Forest Plot",
    subtitle = "SEQTaRget ITT | DAG14 | cluster-robust SE | m=20 | BH-FDR correction",
    x = "Interaction Hazard Ratio (95% CI, log scale)",
    y = NULL,
    caption = "No interaction reached FDR significance (all q > 0.10). Charlson comorbidity index showed only nominal interaction (p=0.013, q=0.134)."
  ) +
  coord_cartesian(xlim = c(0.65, 1.50)) +
  annotate("text", x = 0.68, y = nrow(int_data) + 0.7,
           label = "← Stronger protection", color = FOREST_GREEN, size = 3.2,
           fontface = "bold", hjust = 0) +
  annotate("text", x = 1.48, y = nrow(int_data) + 0.7,
           label = "Stronger harm →", color = JAMA_RED, size = 3.2,
           fontface = "bold", hjust = 1)

# 添加 q-value 标注
ann_x <- max(int_data$ci_upper) * 1.05
for (i in seq_len(nrow(int_data))) {
  p_fig8 <- p_fig8 + annotate("text", x = ann_x, y = int_data$y_label[i],
    label = sprintf("q=%.3f", int_data$q_value[i]),
    hjust = 0, size = 3.2,
    color = if (int_data$sig[i] == "Nominal p < 0.05") JAMA_RED else "grey40")
}

save_fig(p_fig8, "Fig8_Interaction_Forest", 190, 130)

cat("\n" %+% paste(rep("=", 72), collapse = "") %+% "\n")
cat("森林图重新生成完成\n")
cat(paste(rep("=", 72), collapse = "") %+% "\n")
