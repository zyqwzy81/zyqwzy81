# =============================================================================
# V26.9 CONSORT Participant Flow Diagram
# 结构：
#   主列 (x=0.40)
#     [ICU Admissions]              y=0.84
#         ↓
#     [Inclusion: PO4 ≤ 0.80]      y=0.64  ──→ 右侧下方排除大框
#         ↓
#     [Sequential Trial]            y=0.42
#         ↓
#     [TX] x=0.23  [CTRL] x=0.57   y=0.32
#         ↓           ↓
#     [Death]       [Alive]         y=0.22
#
#   排除列 (x=0.80) — Inclusion 下方单一大框
#     [Excluded (n = 44,210)]
#       • Demo exclusion:  23,000
#       • Phosphate excl:   6,611
#       • Washout/other:   14,599
# 输出：Fig1_Participant_Flow.pdf / .png / .tiff
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

OUTPUT_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9/figures"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. 数据
# ============================================================================
icu_admissions  <- 76540L   # MIMIC-IV ICU admissions
excl_demo       <- 23000L   # Age < 18 / non-first ICU / < 24h
excl_phos       <- 6611L    # No phosphate / Phosphate > 0.80
excl_washout    <- 14599L   # Washout failure / ≥1 exclusion
excl_total      <- excl_demo + excl_phos + excl_washout  # 44,210
after_demo      <- icu_admissions - excl_demo             # 53,540
core_cohort     <- after_demo - excl_phos                 # 46,929
n_periods       <- 32330L   # Sequential trial periods
n_patients      <- 22685L   # Unique patients
n_tx            <- 8923L
n_ctrl          <- 23407L
n_death         <- 5610L
n_alive         <- 26720L

# ============================================================================
# 2. 配色 — 期刊级灰度
# ============================================================================
COL_BLACK       <- "#000000"
COL_DARK_GRAY   <- "#2C3E50"
COL_MID_GRAY    <- "#5D6D7E"
COL_LIGHT_GRAY  <- "#95A5A6"
COL_BG_GRAY     <- "#F2F3F4"
COL_WHITE       <- "#FFFFFF"
COL_EXCL_BORDER <- "#7F8C8D"
FONT_FAMILY <- "sans"

# ============================================================================
# 3. 方框定义
# ============================================================================
# 布局结构：
#   主列 (x=0.40)
#     [ICU Admissions]                  y=0.84
#     [Inclusion: Phosphate ≤ 0.80]     y=0.64
#     [Sequential Trial Expansion]      y=0.42
#     [TX] x=0.23  [CTRL] x=0.57       y=0.32
#     [Death] x=0.23  [Alive] x=0.57    y=0.22
#
#   排除列 (x=0.80) — Inclusion 下方
#     [Excluded (n = 44,210)]           y=0.54

# ---- 3a. 主流程（7个框） ----
main_boxes <- data.frame(
  id = 1:7,
  x = c(0.40, 0.40, 0.40,
        0.23, 0.57,
        0.23, 0.57),
  y = c(0.84, 0.64, 0.42,
        0.32, 0.32,
        0.22, 0.22),
  w = c(0.42, 0.42, 0.44,
        0.24, 0.24,
        0.24, 0.24),
  h = c(0.048, 0.048, 0.058,
        0.045, 0.045,
        0.045, 0.045),
  text_size = c(2.9, 2.8, 2.5,
                2.6, 2.6,
                2.6, 2.6),
  stringsAsFactors = FALSE
)

main_labels <- c(
  "ICU Admissions\nAssessed for Eligibility\nn = 76,540",

  "Inclusion: Phosphate \u2264 0.80 mmol/L\nn = 46,929 eligible",

  paste0("Sequential Trial Expansion\n",
         "\u2500\u2500\u2500\u2500\u25b8 24 nested trials\n",
         "\u2500\u2500\u2500\u2500\u25b8 32,330 trial-periods\n",
         "              (22,685 unique patients)"),

  "IV Phosphate\nSupplementation\nn = 8,923 (27.6%)",
  "Conservative\nManagement\nn = 23,407 (72.4%)",

  "28-Day Death\nn = 5,610 (17.4%)",
  "Alive / Censored\nat Day 28\nn = 26,720 (82.6%)"
)

# ---- 3b. 排除方框（单一大框，Inclusion 下方） ----
excl_box <- data.frame(
  id = 1,
  x = 0.80,
  y = 0.54,
  w = 0.32,
  h = 0.170,
  stringsAsFactors = FALSE
)

# 排除框文字
excl_title <- "Excluded (n = 44,210)"

excl_detail <- paste0(
  "\u2460 Demographic exclusion:\n",
  "       Age < 18 / non-first ICU / < 24 h                  23,000\n\n",
  "\u2461 Phosphate-related exclusion:\n",
  "       No phosphate measured / > 0.80 mmol/L         6,611\n\n",
  "\u2462 Washout / other exclusion criteria:\n",
  "       Washout < 48 h / met \u2265 1 exclusion                14,599"
)

excl_footnote <- paste0(
  "Other exclusion criteria (8 total): eGFR < 15, maintenance dialysis, ",
  "hyperparathyroidism, malignancy-related hypercalcemia,\n",
  "phosphate binders > 30 d, pregnancy, elective short-stay surgery, ",
  "multiple ICU admissions per same hospitalization."
)

# ============================================================================
# 4. 箭头
# ============================================================================
# 主流程箭头
main_arrows <- data.frame(
  x0 = c(0.40, 0.40, 0.40,    # Enrollment main + sTTE expansion
         0.23, 0.57,           # TX/CTRL → Death/Alive
         0.23, 0.57),          # Death/Alive → converge
  y0 = c(0.816, 0.616, 0.391,
         0.297, 0.297,
         0.197, 0.197),
  x1 = c(0.40, 0.40, 0.40,
         0.23, 0.57,
         0.40, 0.40),
  y1 = c(0.664, 0.449, 0.374,
         0.242, 0.242,
         0.160, 0.160),
  stringsAsFactors = FALSE
)

# 排除分支箭头（从 Inclusion 框右侧 → 排除框左上方）
excl_arrows <- data.frame(
  x0 = c(0.61),       # Inclusion 框右边缘
  y0 = c(0.64),       # Inclusion 框中心高度
  x1 = c(0.64),       # 排除框左边缘
  y1 = c(0.60),       # 排除框上部区域
  stringsAsFactors = FALSE
)

# ============================================================================
# 5. 阶段标签
# ============================================================================
stage_labels <- data.frame(
  x     = c(0.02, 0.02, 0.02),
  y     = c(0.88, 0.36, 0.26),
  label = c("ENROLLMENT", "ALLOCATION", "FOLLOW-UP"),
  stringsAsFactors = FALSE
)

# ============================================================================
# 6. 绘制
# ============================================================================
p <- ggplot() +

  # ---- 阶段标签 ----
  geom_text(data = stage_labels,
            aes(x = x, y = y, label = label),
            hjust = 0, vjust = 1,
            size = 2.8, fontface = "bold",
            color = COL_LIGHT_GRAY, family = FONT_FAMILY) +

  # ---- 主流程箭头 ----
  geom_segment(data = main_arrows,
               aes(x = x0, y = y0, xend = x1, yend = y1),
               arrow = arrow(length = unit(1.8, "mm"), type = "closed"),
               color = COL_MID_GRAY, linewidth = 0.50) +

  # ---- 排除分支箭头 ----
  geom_segment(data = excl_arrows,
               aes(x = x0, y = y0, xend = x1, yend = y1),
               arrow = arrow(length = unit(1.6, "mm"), type = "closed"),
               color = COL_LIGHT_GRAY, linewidth = 0.45) +

  # ---- 主流程方框 ----
  geom_rect(data = main_boxes,
            aes(xmin = x - w/2, xmax = x + w/2,
                ymin = y - h/2, ymax = y + h/2),
            fill = COL_WHITE, color = COL_BLACK,
            linewidth = 0.55) +

  # ---- 主流程文字 ----
  geom_text(data = data.frame(
              x = main_boxes$x,
              y = main_boxes$y,
              label = main_labels,
              text_size = main_boxes$text_size,
              stringsAsFactors = FALSE
            ),
            aes(x = x, y = y, label = label, size = text_size),
            color = COL_DARK_GRAY, lineheight = 1.20,
            family = FONT_FAMILY, show.legend = FALSE) +
  scale_size_identity() +

  # ---- 排除大框（虚线灰框） ----
  geom_rect(data = excl_box,
            aes(xmin = x - w/2, xmax = x + w/2,
                ymin = y - h/2, ymax = y + h/2),
            fill = COL_BG_GRAY, color = COL_EXCL_BORDER,
            linewidth = 0.50, linetype = "dashed") +

  # ---- 排除框标题（顶部） ----
  geom_text(data = excl_box,
            aes(x = x, y = y + h/2 - 0.014,
                label = excl_title),
            size = 2.4, color = COL_DARK_GRAY,
            family = FONT_FAMILY, fontface = "bold") +

  # ---- 排除框详情（中部） ----
  geom_text(data = excl_box,
            aes(x = x, y = y - 0.010,
                label = excl_detail),
            size = 2.0, color = COL_MID_GRAY,
            family = FONT_FAMILY, lineheight = 1.25,
            hjust = 0.5, vjust = 0.5) +

  # ---- 标题与注释 ----
  labs(
    title    = "Figure 1. Participant Flow Diagram",
    subtitle = paste0(
      "Cloning, Censoring, Weighting (CCW) Framework  |  ",
      "IV Phosphate Replacement for ICU Hypophosphatemia",
      "  |  MIMIC-IV v3.1 (2008\u20132022)"
    ),
    caption  = excl_footnote
  ) +

  coord_cartesian(xlim = c(0, 1), ylim = c(-0.02, 1.0), expand = FALSE) +
  theme_void() +
  theme(
    plot.title    = element_text(size = 13, face = "bold", hjust = 0,
                                  family = FONT_FAMILY, margin = margin(b = 1)),
    plot.subtitle = element_text(size = 8.5, color = COL_MID_GRAY, hjust = 0,
                                  family = FONT_FAMILY, lineheight = 1.3,
                                  margin = margin(b = 3)),
    plot.caption  = element_text(size = 6.0, color = COL_LIGHT_GRAY, hjust = 0,
                                  family = FONT_FAMILY, lineheight = 1.3,
                                  margin = margin(t = 4)),
    plot.margin   = margin(8, 10, 5, 8),
    panel.background = element_rect(fill = COL_WHITE, color = NA),
    plot.background  = element_rect(fill = COL_WHITE, color = NA)
  )

# ============================================================================
# 7. 保存
# ============================================================================
w_cm <- 19.0; h_cm <- 22.0
w_in <- w_cm / 2.54; h_in <- h_cm / 2.54

pdf_file <- file.path(OUTPUT_DIR, "Fig1_Participant_Flow.pdf")
ggsave(pdf_file, p, width = w_in, height = h_in,
       device = cairo_pdf, family = FONT_FAMILY)
cat(sprintf("PDF: %s (%.1f KB)\n", pdf_file, file.info(pdf_file)$size / 1024))

png_file <- file.path(OUTPUT_DIR, "Fig1_Participant_Flow.png")
ggsave(png_file, p, width = w_in, height = h_in, dpi = 300,
       device = "png", family = FONT_FAMILY)
cat(sprintf("PNG: %s (%.1f KB)\n", png_file, file.info(png_file)$size / 1024))

tiff_file <- file.path(OUTPUT_DIR, "Fig1_Participant_Flow.tiff")
grDevices::tiff(tiff_file, width = w_in, height = h_in,
                units = "in", res = 1200, compression = "lzw",
                family = FONT_FAMILY)
print(p)
dev.off()
cat(sprintf("TIFF: %s (%.1f KB)\n", tiff_file, file.info(tiff_file)$size / 1024))

# ============================================================================
# 8. 验证
# ============================================================================
cat("\n===== 验证 =====\n")
cat(sprintf("主方框: %d | 排除大框: 1 | 箭头: %d + %d\n",
            nrow(main_boxes), nrow(main_arrows), nrow(excl_arrows)))
cat(sprintf("总排除: %d = %d (demo) + %d (phos) + %d (washout)\n",
            excl_total, excl_demo, excl_phos, excl_washout))
cat(sprintf("Periods: TX %d + CTRL %d = %d (expected %d)\n",
            n_tx, n_ctrl, n_tx + n_ctrl, n_periods))
cat(sprintf("Outcome: Death %d + Alive %d = %d (expected %d)\n",
            n_death, n_alive, n_death + n_alive, n_periods))

cat("\n输出:\n")
for (f in c(png_file, pdf_file, tiff_file)) {
  if (file.exists(f))
    cat(sprintf("  \u2713 %s (%.1f KB)\n", basename(f), file.info(f)$size / 1024))
}
cat("\n===== 完成 =====\n")
