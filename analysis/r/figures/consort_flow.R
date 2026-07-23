# =============================================================================
# V28 CONSORT 参与者流程图 (Fig 1A) — 美化版
# 修复点：
#   1. 统一所有主流程方框字号（消除 text_size 4 种不一致值）
#   2. 方框高度按文字行数自适应（避免文字溢出/重叠）
#   3. 排除框使用数据框 + 逐行 geom_text（左对齐文字 + 右对齐数字）
#      替代原多空格对齐方案，彻底解决数字不对齐问题
#   4. 箭头起止点 = 方框边缘精确对齐
#   5. 画布扩大至 22 × 26 cm，留足留白
# 输出：Fig1A_Flow_Diagram.pdf / .png / .tiff
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

OUTPUT_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v28/output/figures"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. 数据 —— V28 版本（与原脚本一致）
# ============================================================================
icu_admissions  <- 76540L   # MIMIC-IV ICU admissions
excl_demo       <- 23000L   # Age < 18 / non-first ICU / < 24h
excl_phos       <- 6611L    # No phosphate / Phosphate > 0.80
excl_washout    <- 14599L   # Washout failure / ≥1 exclusion
excl_total      <- excl_demo + excl_phos + excl_washout  # 44,210
n_periods       <- 32330L   # Sequential trial periods
n_patients      <- 22685L   # Unique patients
n_tx            <- 5147L    # IV phosphate (15.9%)
n_ctrl          <- 27183L   # Conservative (84.1%)
n_death         <- 5610L    # 28-day death
n_alive         <- 26720L   # Alive / censored

fmt <- function(x) formatC(x, format = "d", big.mark = ",")

# ============================================================================
# 2. 配色 —— 期刊级灰度
# ============================================================================
COL_BLACK       <- "#000000"
COL_DARK_GRAY   <- "#2C3E50"
COL_MID_GRAY    <- "#5D6D7E"
COL_LIGHT_GRAY  <- "#95A5A6"
COL_BG_GRAY     <- "#F2F3F4"
COL_WHITE       <- "#FFFFFF"
COL_EXCL_BORDER <- "#7F8C8D"
FONT_SANS <- "sans"
FONT_MONO <- "mono"   # 等宽字体，用于排除框数字对齐

# ============================================================================
# 3. 主流程方框定义（统一字号、高度自适应）
# ============================================================================
# 统一字号
SIZE_MAIN   <- 3.05   # 主流程方框文字
SIZE_BRANCH <- 3.00   # 分支方框（分组/结局）
SIZE_STAGE  <- 3.00   # 阶段标签

# 主流程方框（垂直居中列 x=0.32；右侧留给排除框）
# 字段: id, x, y, w, h, label, size
main_boxes <- data.frame(
  id = 1:7,
  x  = c(0.32, 0.32, 0.32,
         0.18, 0.46,
         0.18, 0.46),
  y  = c(0.920, 0.800, 0.660,
         0.470, 0.470,
         0.285, 0.285),
  w  = c(0.50, 0.50, 0.50,
         0.28, 0.28,
         0.28, 0.28),
  # 高度按行数自适应（每行约 0.022 npc）
  h  = c(0.066, 0.050, 0.090,
         0.072, 0.072,
         0.072, 0.072),
  size = c(SIZE_MAIN, SIZE_MAIN, SIZE_MAIN,
           SIZE_BRANCH, SIZE_BRANCH,
           SIZE_BRANCH, SIZE_BRANCH),
  stringsAsFactors = FALSE
)

main_labels <- c(
  sprintf("ICU Admissions\nAssessed for Eligibility\nn = %s", fmt(icu_admissions)),
  sprintf("Inclusion: Phosphate \u2264 0.80 mmol/L\nn = %s eligible", fmt(46929)),
  sprintf("Sequential Trial Expansion\n\u2500\u2500\u2500\u2500\u25b8 24 nested trials\n\u2500\u2500\u2500\u2500\u25b8 %s trial-periods\n              (%s unique patients)",
          fmt(n_periods), fmt(n_patients)),
  sprintf("IV Phosphate\nSupplementation\nn = %s (15.9%%)", fmt(n_tx)),
  sprintf("Conservative\nManagement\nn = %s (84.1%%)", fmt(n_ctrl)),
  sprintf("28-Day Death\nn = %s (17.4%%)", fmt(n_death)),
  sprintf("Alive / Censored\nat Day 28\nn = %s (82.6%%)", fmt(n_alive))
)

# ============================================================================
# 4. 排除框 —— 数据框形式，逐行绘制实现严格对齐
# ============================================================================
# 排除框整体位置（右移至 x=0.80，与3号方框右边 0.57 留 0.04 间隙避免视觉重叠）
EXCL_X <- 0.80
EXCL_Y <- 0.47
EXCL_W <- 0.38
EXCL_H <- 0.52

EXCL_TITLE <- sprintf("Excluded (n = %s)", fmt(excl_total))

# 排除条目：label 左对齐，count 右对齐
excl_items <- data.frame(
  label = c(
    "\u2460 Demographic exclusion",
    "   Age < 18 / non-first ICU / < 24 h",
    "\u2461 Phosphate-related exclusion",
    "   No phosphate / > 0.80 mmol/L",
    "\u2462 Washout < 48 h",
    "\u2463 Medical exclusion (8 criteria)",
    "   eGFR < 15",
    "   Short stay (elective surgery)",
    "   Maintenance dialysis",
    "   Malignancy hypercalcemia",
    "   Hyperparathyroidism",
    "   Pregnancy",
    "   Phosphate binders > 30 d",
    "   Multiple ICU admissions"
  ),
  count = c(
    fmt(excl_demo), "", fmt(excl_phos), "", "11,637", "2,962",
    "1,374", "796", "269", "215", "212", "92", "4", "(incl.)"
  ),
  # 行类型：header（粗体）/ item（普通）/ sub（缩进）
  type = c(
    "header", "sub",
    "header", "sub",
    "header",
    "header", "sub", "sub", "sub", "sub", "sub", "sub", "sub", "sub"
  ),
  stringsAsFactors = FALSE
)

# 计算每行 y 坐标（从框顶部向下排列）
n_excl_lines <- nrow(excl_items)
line_height  <- 0.030  # 每行间距
excl_title_y <- EXCL_Y + EXCL_H/2 - 0.025  # 标题靠顶
excl_first_y  <- EXCL_Y + EXCL_H/2 - 0.055  # 第一行内容
excl_items$y  <- excl_first_y - (seq_len(n_excl_lines) - 1) * line_height

# 文字左右对齐位置
excl_label_x <- EXCL_X - EXCL_W/2 + 0.012  # 左对齐
excl_count_x <- EXCL_X + EXCL_W/2 - 0.012  # 右对齐

# ============================================================================
# 5. 箭头 —— 起止点严格对齐方框边缘
# ============================================================================
# 主流程箭头：从方框底边 (y - h/2) 到下一方框顶边 (y + h/2)
main_arrows <- data.frame(
  x0 = c(0.32, 0.32, 0.32,
         0.18, 0.46,
         0.18, 0.46),
  y0 = c(main_boxes$y[1] - main_boxes$h[1]/2,
         main_boxes$y[2] - main_boxes$h[2]/2,
         main_boxes$y[3] - main_boxes$h[3]/2,
         main_boxes$y[4] + main_boxes$h[4]/2,
         main_boxes$y[5] + main_boxes$h[5]/2,
         main_boxes$y[6] + main_boxes$h[6]/2,
         main_boxes$y[7] + main_boxes$h[7]/2),
  x1 = c(0.32, 0.32, 0.32,
         0.18, 0.46,
         0.18, 0.46),
  y1 = c(main_boxes$y[2] + main_boxes$h[2]/2,
         main_boxes$y[3] + main_boxes$h[3]/2,
         main_boxes$y[4] + main_boxes$h[4]/2,  # 主→分支1（左）
         main_boxes$y[6] + main_boxes$h[6]/2,
         main_boxes$y[7] + main_boxes$h[7]/2,
         main_boxes$y[6] - main_boxes$h[6]/2,  # 占位，实际不用
         main_boxes$y[7] - main_boxes$h[7]/2),
  type = c("main", "main", "split",
           "branch", "branch",
           "placeholder", "placeholder"),
  stringsAsFactors = FALSE)

# 第3方框 → 两个分支方框的"Y型"分叉箭头
# 主干：3号方框底 → 分叉点
# 分叉点 → 4号、5号方框顶
split_arrows <- data.frame(
  x0 = c(0.32, 0.18, 0.46),
  y0 = c(main_boxes$y[3] - main_boxes$h[3]/2,
         main_boxes$y[4] + main_boxes$h[4]/2 + 0.040,
         main_boxes$y[5] + main_boxes$h[5]/2 + 0.040),
  x1 = c(0.32, 0.18, 0.46),
  y1 = c(main_boxes$y[4] + main_boxes$h[4]/2 + 0.040,
         main_boxes$y[4] + main_boxes$h[4]/2,
         main_boxes$y[5] + main_boxes$h[5]/2),
  stringsAsFactors = FALSE
)

# 横向连接线（分叉点之间的水平线）
split_horiz <- data.frame(
  x0 = 0.18, y0 = main_boxes$y[4] + main_boxes$h[4]/2 + 0.040,
  x1 = 0.46, y1 = main_boxes$y[5] + main_boxes$h[5]/2 + 0.040
)

# 排除框分支箭头（从方框2右侧 → 排除框左侧）
excl_arrows <- data.frame(
  x0 = main_boxes$x[2] + main_boxes$w[2]/2,  # 0.57
  y0 = main_boxes$y[2],
  x1 = EXCL_X - EXCL_W/2,                    # 0.57
  y1 = EXCL_Y + EXCL_H/2 - 0.04
)

# ============================================================================
# 6. 阶段标签
# ============================================================================
stage_labels <- data.frame(
  x     = c(0.015, 0.015, 0.015),
  y     = c(0.960, 0.560, 0.365),
  label = c("ENROLLMENT", "ALLOCATION", "FOLLOW-UP"),
  stringsAsFactors = FALSE
)

# ============================================================================
# 7. 绘图
# ============================================================================
# 准备排除框逐行数据：header 行用粗体，sub 行用普通体
excl_items$fontface <- ifelse(excl_items$type == "header", "bold", "plain")
excl_count_df <- subset(excl_items, count != "")

p <- ggplot() +

  # ---- 阶段标签 ----
  geom_text(data = stage_labels,
            aes(x = x, y = y, label = label),
            hjust = 0, vjust = 0.5,
            size = SIZE_STAGE, fontface = "bold",
            color = COL_LIGHT_GRAY, family = FONT_SANS) +

  # ---- 主流程箭头（垂直主干：1→2, 2→3） ----
  geom_segment(data = main_arrows[1:2, ],
               aes(x = x0, y = y0, xend = x1, yend = y1),
               arrow = arrow(length = unit(2.0, "mm"), type = "closed"),
               color = COL_MID_GRAY, linewidth = 0.55) +

  # ---- Y型分叉箭头（3→分支） ----
  # 主干（3号底 → 水平横线中点上方）
  geom_segment(data = data.frame(
                x0 = 0.32, y0 = main_boxes$y[3] - main_boxes$h[3]/2,
                x1 = 0.32, y1 = main_boxes$y[4] + main_boxes$h[4]/2 + 0.040),
               aes(x = x0, y = y0, xend = x1, yend = y1),
               color = COL_MID_GRAY, linewidth = 0.55) +
  # 水平横线
  geom_segment(data = split_horiz,
               aes(x = x0, y = y0, xend = x1, yend = y1),
               color = COL_MID_GRAY, linewidth = 0.55) +
  # 左分支（横线 → 4号顶）
  geom_segment(data = data.frame(
                x0 = 0.18, y0 = main_boxes$y[4] + main_boxes$h[4]/2 + 0.040,
                x1 = 0.18, y1 = main_boxes$y[4] + main_boxes$h[4]/2),
               aes(x = x0, y = y0, xend = x1, yend = y1),
               arrow = arrow(length = unit(2.0, "mm"), type = "closed"),
               color = COL_MID_GRAY, linewidth = 0.55) +
  # 右分支（横线 → 5号顶）
  geom_segment(data = data.frame(
                x0 = 0.46, y0 = main_boxes$y[5] + main_boxes$h[5]/2 + 0.040,
                x1 = 0.46, y1 = main_boxes$y[5] + main_boxes$h[5]/2),
               aes(x = x0, y = y0, xend = x1, yend = y1),
               arrow = arrow(length = unit(2.0, "mm"), type = "closed"),
               color = COL_MID_GRAY, linewidth = 0.55) +

  # ---- 分支 → 结局 箭头（4→6, 5→7） ----
  geom_segment(data = data.frame(
                x0 = 0.18, y0 = main_boxes$y[4] - main_boxes$h[4]/2,
                x1 = 0.18, y1 = main_boxes$y[6] + main_boxes$h[6]/2),
               aes(x = x0, y = y0, xend = x1, yend = y1),
               arrow = arrow(length = unit(2.0, "mm"), type = "closed"),
               color = COL_MID_GRAY, linewidth = 0.55) +
  geom_segment(data = data.frame(
                x0 = 0.46, y0 = main_boxes$y[5] - main_boxes$h[5]/2,
                x1 = 0.46, y1 = main_boxes$y[7] + main_boxes$h[7]/2),
               aes(x = x0, y = y0, xend = x1, yend = y1),
               arrow = arrow(length = unit(2.0, "mm"), type = "closed"),
               color = COL_MID_GRAY, linewidth = 0.55) +

  # ---- 排除分支箭头 ----
  geom_segment(data = excl_arrows,
               aes(x = x0, y = y0, xend = x1, yend = y1),
               arrow = arrow(length = unit(1.8, "mm"), type = "closed"),
               color = COL_LIGHT_GRAY, linewidth = 0.50,
               linetype = "dashed") +

  # ---- 主流程方框 ----
  geom_rect(data = main_boxes,
            aes(xmin = x - w/2, xmax = x + w/2,
                ymin = y - h/2, ymax = y + h/2),
            fill = COL_WHITE, color = COL_BLACK,
            linewidth = 0.60) +

  # ---- 主流程文字（统一字号，垂直居中） ----
  geom_text(data = data.frame(
              x = main_boxes$x,
              y = main_boxes$y,
              label = main_labels,
              size = main_boxes$size,
              stringsAsFactors = FALSE
            ),
            aes(x = x, y = y, label = label, size = size),
            color = COL_DARK_GRAY, lineheight = 1.15,
            family = FONT_SANS, show.legend = FALSE) +
  scale_size_identity() +

  # ---- 排除大框 ----
  geom_rect(data = data.frame(x = EXCL_X, y = EXCL_Y, w = EXCL_W, h = EXCL_H),
            aes(xmin = x - w/2, xmax = x + w/2,
                ymin = y - h/2, ymax = y + h/2),
            fill = COL_BG_GRAY, color = COL_EXCL_BORDER,
            linewidth = 0.55, linetype = "dashed") +

  # ---- 排除框标题 ----
  geom_text(data = data.frame(x = EXCL_X, y = excl_title_y, label = EXCL_TITLE),
            aes(x = x, y = y, label = label),
            size = SIZE_MAIN, color = COL_DARK_GRAY,
            family = FONT_SANS, fontface = "bold") +

  # ---- 排除框详情：逐行绘制，左对齐文字 + 右对齐数字 ----
  # 文字部分（左对齐）
  geom_text(data = excl_items,
            aes(x = excl_label_x, y = y, label = label, fontface = I(fontface)),
            hjust = 0, vjust = 0.5,
            size = 2.35,
            color = COL_DARK_GRAY,
            family = FONT_MONO) +
  # 数字部分（右对齐）
  geom_text(data = excl_count_df,
            aes(x = excl_count_x, y = y, label = count),
            hjust = 1, vjust = 0.5,
            size = 2.35,
            color = COL_MID_GRAY,
            family = FONT_MONO,
            fontface = "bold") +

  # ---- 标题与注释 ----
  labs(
    title    = "Figure 1A. Participant Flow Diagram",
    subtitle = paste0(
      "Cloning, Censoring, Weighting (CCW) Framework  |  ",
      "IV Phosphate Replacement for ICU Hypophosphatemia",
      "  |  MIMIC-IV v3.1 (2008\u20132022)"
    ),
    caption  = paste0(
      "Exclusion criteria: eGFR < 15, maintenance dialysis, hyperparathyroidism, ",
      "malignancy-related hypercalcemia, phosphate binders > 30 d, pregnancy, ",
      "elective short-stay surgery, multiple ICU admissions per same hospitalization. ",
      "Exclusion counts from study_cohort table (MIMIC-IV v3.1, 2008-2022)."
    )
  ) +

  coord_cartesian(xlim = c(0, 1), ylim = c(-0.02, 1.0), expand = FALSE) +
  theme_void() +
  theme(
    plot.title    = element_text(size = 13, face = "bold", hjust = 0,
                                family = FONT_SANS, margin = margin(b = 1)),
    plot.subtitle = element_text(size = 9.0, color = COL_MID_GRAY, hjust = 0,
                                 family = FONT_SANS, lineheight = 1.3,
                                 margin = margin(b = 4)),
    plot.caption  = element_text(size = 6.0, color = COL_LIGHT_GRAY, hjust = 0,
                                 family = FONT_SANS, lineheight = 1.3,
                                 margin = margin(t = 5)),
    plot.margin   = margin(8, 12, 6, 8),
    panel.background = element_rect(fill = COL_WHITE, color = NA),
    plot.background  = element_rect(fill = COL_WHITE, color = NA)
  )

# ============================================================================
# 8. 保存
# ============================================================================
w_cm <- 22.0; h_cm <- 26.0
w_in <- w_cm / 2.54; h_in <- h_cm / 2.54

pdf_file <- file.path(OUTPUT_DIR, "Fig1A_Flow_Diagram.pdf")
ggsave(pdf_file, p, width = w_in, height = h_in,
       device = cairo_pdf, family = FONT_SANS)
cat(sprintf("PDF: %s (%.1f KB)\n", pdf_file, file.info(pdf_file)$size / 1024))

png_file <- file.path(OUTPUT_DIR, "Fig1A_Flow_Diagram.png")
ggsave(png_file, p, width = w_in, height = h_in, dpi = 600,
       device = "png", type = "cairo")
cat(sprintf("PNG: %s (%.1f KB)\n", png_file, file.info(png_file)$size / 1024))

tiff_file <- file.path(OUTPUT_DIR, "Fig1A_Flow_Diagram.tiff")
grDevices::tiff(tiff_file, width = w_in, height = h_in,
                units = "in", res = 600, compression = "lzw")
print(p)
dev.off()
cat(sprintf("TIFF: %s (%.1f KB)\n", tiff_file, file.info(tiff_file)$size / 1024))

# ============================================================================
# 9. 验证
# ============================================================================
cat("\n===== 验证 =====\n")
cat(sprintf("主方框: %d | 排除条目: %d 行\n", nrow(main_boxes), nrow(excl_items)))
cat(sprintf("总排除: %d = %d (demo) + %d (phos) + %d (washout)\n",
            excl_total, excl_demo, excl_phos, excl_washout))
cat(sprintf("Periods: TX %d + CTRL %d = %d (expected %d)\n",
            n_tx, n_ctrl, n_tx + n_ctrl, n_periods))
cat(sprintf("Outcome: Death %d + Alive %d = %d (expected %d)\n",
            n_death, n_alive, n_death + n_alive, n_periods))
cat(sprintf("字号统一: 主流程 %.2f / 分支 %.2f / 阶段 %.2f\n",
            SIZE_MAIN, SIZE_BRANCH, SIZE_STAGE))

cat("\n输出:\n")
for (f in c(png_file, pdf_file, tiff_file)) {
  if (file.exists(f))
    cat(sprintf("  \u2713 %s (%.1f KB)\n", basename(f), file.info(f)$size / 1024))
}
cat("\n===== 完成 =====\n")
