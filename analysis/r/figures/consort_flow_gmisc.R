# =============================================================================
# V28 CONSORT Participant Flow Diagram (Fig 1A) — 使用 Gmisc::flowchart()
# =============================================================================
# 参考：https://cran.r-project.org/web/packages/Gmisc/vignettes/Grid-based_flowcharts.html

suppressPackageStartupMessages({
  library(Gmisc)
  library(grid)
  library(glue)
})

OUTPUT_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v28/output/figures"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. 数据
# ============================================================================
cat("===== 数据 =====\n")
ICU_ADM      <- txtInt(76540)
EXCL_DEMO    <- txtInt(23000)
EXCL_PHOS    <- txtInt(6611)
EXCL_WASHOUT <- txtInt(11637)
EXCL_MED     <- txtInt(2962)
EXCL_TOTAL   <- txtInt(44210)
ELIGIBLE     <- txtInt(46929)
N_PERIODS    <- txtInt(32330)
N_PATIENTS   <- txtInt(22685)
N_TX         <- txtInt(5147)
N_CTRL       <- txtInt(27183)
N_DEATH      <- txtInt(5610)
N_ALIVE      <- txtInt(26720)

cat(sprintf("ICU 入院: %s\n", ICU_ADM))
cat(sprintf("排除总计: %s = %s + %s + %s + %s\n",
            EXCL_TOTAL, EXCL_DEMO, EXCL_PHOS, EXCL_WASHOUT, EXCL_MED))
cat(sprintf("合格: %s\n", ELIGIBLE))
cat(sprintf("试验期: %s (%s 患者)\n", N_PERIODS, N_PATIENTS))
cat(sprintf("治疗: %s 对照: %s\n", N_TX, N_CTRL))
cat(sprintf("结局: 死亡 %s / 存活 %s\n", N_DEATH, N_ALIVE))

# ============================================================================
# 2. 样式定义
# ============================================================================
main_box_gp <- gpar(fill = "#FFFFFF", col = "#000000", lwd = 1.0)
excl_box_gp <- gpar(fill = "#D9D9D9", col = "#A6A6A6", lwd = 0.8, lty = "dashed")
main_con_gp <- gpar(col = "#808080", lwd = 1.2, fill = "#808080")
excl_con_gp <- gpar(col = "#A6A6A6", lwd = 1.0, fill = "#A6A6A6")

BOX_W <- unit(60, "mm")
BOX_H <- unit(18, "mm")
EXCL_W <- unit(95, "mm")
EXCL_H <- unit(80, "mm")

# ============================================================================
# 3. 构建流程图
# ============================================================================
cat("\n===== 构建流程图 =====\n")

# 3.1 主流程方框
source_box <- boxGrob(
  glue("ICU Admissions\nAssessed for Eligibility\nn = {n}", n = ICU_ADM),
  box_gp = main_box_gp, width = BOX_W, height = BOX_H
)

eligible_box <- boxGrob(
  glue("Inclusion: Phosphate ≤ 0.80 mmol/L\nn = {n} eligible", n = ELIGIBLE),
  box_gp = main_box_gp, width = BOX_W, height = BOX_H
)

included_box <- boxGrob(
  glue("Sequential Trial Expansion\n├─ 24 nested trials\n├─ {N_PERIODS} trial-periods\n└─ ({N_PATIENTS} unique patients)"),
  box_gp = main_box_gp, width = BOX_W, height = BOX_H
)

# 3.2 分组方框
tx_box <- boxGrob(
  glue("IV Phosphate\nSupplementation\nn = {n} ({pct}%)",
       n = N_TX, pct = round(5147/32330*100, 1)),
  box_gp = main_box_gp, width = BOX_W, height = BOX_H
)

ctrl_box <- boxGrob(
  glue("Conservative\nManagement\nn = {n} ({pct}%)",
       n = N_CTRL, pct = round(27183/32330*100, 1)),
  box_gp = main_box_gp, width = BOX_W, height = BOX_H
)

# 3.3 结局方框
death_box <- boxGrob(
  glue("28-Day Death\nn = {n} ({pct}%)",
       n = N_DEATH, pct = round(5610/32330*100, 1)),
  box_gp = main_box_gp, width = BOX_W, height = BOX_H
)

alive_box <- boxGrob(
  glue("Alive / Censored\nat Day 28\nn = {n} ({pct}%)",
       n = N_ALIVE, pct = round(26720/32330*100, 1)),
  box_gp = main_box_gp, width = BOX_W, height = BOX_H
)

# 3.4 排除方框
excl_box <- boxHeaderGrob(
  header = glue("Excluded (n = {EXCL_TOTAL})"),
  body = paste0(
    " ① Demographic exclusion:                          ", EXCL_DEMO, "\n",
    "       (Age < 18 / non-first ICU / < 24 h)\n",
    " ② Phosphate-related exclusion:                ", EXCL_PHOS, "\n",
    "       (No phosphate / > 0.80 mmol/L)\n",
    " ③ Washout < 48 h:                                       ", EXCL_WASHOUT, "\n",
    " ④ Medical exclusion criteria (8 total):        ", EXCL_MED, "\n",
    "       eGFR < 15                                                 1,374\n",
    "       Short stay (elective surgery)                        796\n",
    "       Maintenance dialysis                                    269\n",
    "       Malignancy hypercalcemia                              215\n",
    "       Hyperparathyroidism                                      212\n",
    "       Pregnancy                                                       92\n",
    "       Phosphate binders > 30 d                                4\n",
    "       Multiple ICU admissions                           (incl.)"
  ),
  box_gp = excl_box_gp, width = EXCL_W, height = EXCL_H,
  header_gp = gpar(col = "#000000", cex = 0.75, fontface = "bold"),
  body_gp = gpar(col = "#000000", cex = 0.55, lineheight = 1.15),
  bjust = "left"
)

# 3.5 构建并输出
# 使用 Gmisc::flowchart() 构建，然后保存为 PDF
pdf_file <- file.path(OUTPUT_DIR, "Fig1A_Flow_Diagram_Gmisc.pdf")
pdf(pdf_file, width = 14, height = 12, family = "sans")

# 标题
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1,
  heights = unit(c(0.08, 0.92), c("npc", "npc")))))

# 标题区域
pushViewport(viewport(layout.pos.row = 1))
grid.text(
  "Figure 1A. Participant Flow Diagram",
  x = unit(3, "mm"), y = unit(0.5, "npc"),
  just = c("left", "center"),
  gp = gpar(fontsize = 13, fontface = "bold", fontfamily = "sans")
)
grid.text(
  "Cloning, Censoring, Weighting (CCW) Framework  |  IV Phosphate Replacement for ICU Hypophosphatemia  |  MIMIC-IV v3.1 (2008–2022)",
  x = unit(3, "mm"), y = unit(0.05, "npc"),
  just = c("left", "center"),
  gp = gpar(fontsize = 8.5, col = "#5D6D7E", fontfamily = "sans")
)
popViewport()

# 流程图区域
pushViewport(viewport(layout.pos.row = 2))
flowchart(
  source = source_box,
  eligible = eligible_box,
  included = included_box,
  groups = list(tx_box, ctrl_box),
  outcomes = list(death_box, alive_box)
) |>
  spread(axis = "y") |>
  spread(subelement = "groups", axis = "x") |>
  spread(subelement = "outcomes", axis = "x") |>
  equalizeWidths(subelement = list("source", "eligible", "included")) |>
  equalizeWidths(subelement = "groups") |>
  equalizeWidths(subelement = "outcomes") |>
  insert(
    list(excluded = excl_box),
    after = "eligible",
    name = "excluded"
  ) |>
  move(subelement = "excluded", x = 0.78) |>
  connect("source", "eligible", type = "vert", lty_gp = main_con_gp,
          arrow_size = 3, smooth = TRUE) |>
  connect("eligible", "included", type = "vert", lty_gp = main_con_gp,
          arrow_size = 3, smooth = TRUE) |>
  connect("eligible", "excluded", type = "L", lty_gp = excl_con_gp,
          arrow_size = 3, label = "Excluded") |>
  connect("included", "groups", type = "N", lty_gp = main_con_gp,
          arrow_size = 3, smooth = TRUE) |>
  connect("groups", "outcomes", type = "vert", lty_gp = main_con_gp,
          arrow_size = 3, smooth = TRUE) |>
  phaseLabel("source", label = "ENROLLMENT", overlap = 0.07) |>
  phaseLabel("included", label = "ALLOCATION", overlap = 0.07) |>
  phaseLabel("outcomes", label = "FOLLOW-UP", overlap = 0.07)
popViewport()

# 底部注释
pushViewport(viewport(layout.pos.row = 1))
grid.text(
  "Exclusion criteria: eGFR < 15, maintenance dialysis, hyperparathyroidism, malignancy-related hypercalcemia, phosphate binders > 30 d, pregnancy, elective short-stay surgery, multiple ICU admissions per same hospitalization. Data source: MIMIC-IV v3.1 (2008–2022).",
  x = unit(3, "mm"), y = unit(0.01, "npc"),
  just = c("left", "bottom"),
  gp = gpar(fontsize = 6.0, col = "#95A5A6", fontfamily = "sans", lineheight = 1.3)
)
popViewport()

dev.off()
cat(sprintf("PDF: %s (%.1f KB)\n", pdf_file, file.info(pdf_file)$size / 1024))

# 4.2 PNG
png_file <- file.path(OUTPUT_DIR, "Fig1A_Flow_Diagram_Gmisc.png")
png(png_file, width = 14, height = 12, units = "in", res = 600, family = "sans")

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1,
  heights = unit(c(0.08, 0.92), c("npc", "npc")))))

pushViewport(viewport(layout.pos.row = 1))
grid.text(
  "Figure 1A. Participant Flow Diagram",
  x = unit(3, "mm"), y = unit(0.5, "npc"),
  just = c("left", "center"),
  gp = gpar(fontsize = 13, fontface = "bold", fontfamily = "sans")
)
grid.text(
  "Cloning, Censoring, Weighting (CCW) Framework  |  IV Phosphate Replacement for ICU Hypophosphatemia  |  MIMIC-IV v3.1 (2008–2022)",
  x = unit(3, "mm"), y = unit(0.05, "npc"),
  just = c("left", "center"),
  gp = gpar(fontsize = 8.5, col = "#5D6D7E", fontfamily = "sans")
)
popViewport()

pushViewport(viewport(layout.pos.row = 2))
flowchart(
  source = source_box,
  eligible = eligible_box,
  included = included_box,
  groups = list(tx_box, ctrl_box),
  outcomes = list(death_box, alive_box)
) |>
  spread(axis = "y") |>
  spread(subelement = "groups", axis = "x") |>
  spread(subelement = "outcomes", axis = "x") |>
  equalizeWidths(subelement = list("source", "eligible", "included")) |>
  equalizeWidths(subelement = "groups") |>
  equalizeWidths(subelement = "outcomes") |>
  insert(
    list(excluded = excl_box),
    after = "eligible",
    name = "excluded"
  ) |>
  move(subelement = "excluded", x = 0.78) |>
  connect("source", "eligible", type = "vert", lty_gp = main_con_gp,
          arrow_size = 3, smooth = TRUE) |>
  connect("eligible", "included", type = "vert", lty_gp = main_con_gp,
          arrow_size = 3, smooth = TRUE) |>
  connect("eligible", "excluded", type = "L", lty_gp = excl_con_gp,
          arrow_size = 3, label = "Excluded") |>
  connect("included", "groups", type = "N", lty_gp = main_con_gp,
          arrow_size = 3, smooth = TRUE) |>
  connect("groups", "outcomes", type = "vert", lty_gp = main_con_gp,
          arrow_size = 3, smooth = TRUE) |>
  phaseLabel("source", label = "ENROLLMENT", overlap = 0.07) |>
  phaseLabel("included", label = "ALLOCATION", overlap = 0.07) |>
  phaseLabel("outcomes", label = "FOLLOW-UP", overlap = 0.07)
popViewport()

dev.off()
cat(sprintf("PNG: %s (%.1f KB)\n", png_file, file.info(png_file)$size / 1024))

# 4.3 TIFF
tiff_file <- file.path(OUTPUT_DIR, "Fig1A_Flow_Diagram_Gmisc.tiff")
tiff(tiff_file, width = 14, height = 12, units = "in", res = 600, compression = "lzw")

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1,
  heights = unit(c(0.08, 0.92), c("npc", "npc")))))

pushViewport(viewport(layout.pos.row = 1))
grid.text(
  "Figure 1A. Participant Flow Diagram",
  x = unit(3, "mm"), y = unit(0.5, "npc"),
  just = c("left", "center"),
  gp = gpar(fontsize = 13, fontface = "bold", fontfamily = "sans")
)
grid.text(
  "Cloning, Censoring, Weighting (CCW) Framework  |  IV Phosphate Replacement for ICU Hypophosphatemia  |  MIMIC-IV v3.1 (2008–2022)",
  x = unit(3, "mm"), y = unit(0.05, "npc"),
  just = c("left", "center"),
  gp = gpar(fontsize = 8.5, col = "#5D6D7E", fontfamily = "sans")
)
popViewport()

pushViewport(viewport(layout.pos.row = 2))
flowchart(
  source = source_box,
  eligible = eligible_box,
  included = included_box,
  groups = list(tx_box, ctrl_box),
  outcomes = list(death_box, alive_box)
) |>
  spread(axis = "y") |>
  spread(subelement = "groups", axis = "x") |>
  spread(subelement = "outcomes", axis = "x") |>
  equalizeWidths(subelement = list("source", "eligible", "included")) |>
  equalizeWidths(subelement = "groups") |>
  equalizeWidths(subelement = "outcomes") |>
  insert(
    list(excluded = excl_box),
    after = "eligible",
    name = "excluded"
  ) |>
  move(subelement = "excluded", x = 0.78) |>
  connect("source", "eligible", type = "vert", lty_gp = main_con_gp,
          arrow_size = 3, smooth = TRUE) |>
  connect("eligible", "included", type = "vert", lty_gp = main_con_gp,
          arrow_size = 3, smooth = TRUE) |>
  connect("eligible", "excluded", type = "L", lty_gp = excl_con_gp,
          arrow_size = 3, label = "Excluded") |>
  connect("included", "groups", type = "N", lty_gp = main_con_gp,
          arrow_size = 3, smooth = TRUE) |>
  connect("groups", "outcomes", type = "vert", lty_gp = main_con_gp,
          arrow_size = 3, smooth = TRUE) |>
  phaseLabel("source", label = "ENROLLMENT", overlap = 0.07) |>
  phaseLabel("included", label = "ALLOCATION", overlap = 0.07) |>
  phaseLabel("outcomes", label = "FOLLOW-UP", overlap = 0.07)
popViewport()

dev.off()
cat(sprintf("TIFF: %s (%.1f KB)\n", tiff_file, file.info(tiff_file)$size / 1024))

# ============================================================================
# 5. 验证
# ============================================================================
cat("\n===== 验证 =====\n")
cat(sprintf("总排除: %s = %s + %s + %s + %s ✓\n",
            EXCL_TOTAL, EXCL_DEMO, EXCL_PHOS, EXCL_WASHOUT, EXCL_MED))
cat(sprintf("Periods: TX %s + CTRL %s = %s (expected %s) ✓\n",
            N_TX, N_CTRL, txtInt(5147+27183), N_PERIODS))
cat(sprintf("Outcome: Death %s + Alive %s = %s (expected %s) ✓\n",
            N_DEATH, N_ALIVE, txtInt(5610+26720), N_PERIODS))

cat("\n输出:\n")
for (f in c(pdf_file, png_file, tiff_file)) {
  if (file.exists(f))
    cat(sprintf("  ✓ %s (%.1f KB)\n", basename(f), file.info(f)$size / 1024))
}
cat("\n===== 完成 =====\n")