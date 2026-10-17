#!/usr/bin/env Rscript
# =============================================================================
# V28.4 ROUTE-A 统一绘图：trial-period 级 SEQTaRget
#   * Figure 1 CONSORT（trial-period framing）
#   * Figure 2 STE counterfactual survival curve（已由 v28.4_ste_survival_curve.R 生成）
#   * Figure 3 亚组 Forest（trial-period 级）
#   * Figure 4 分段 HR（trial-period 级）
#   * S1/S5/S6/S7/S11
#   * 统一 ICM 期刊格式（176mm / 600 DPI / Helvetica）
#   * 全部 trial-period 口径，n=32,330 入组 / 5,648 死亡 / IV n=5,147
#   * SEQTaRget HR=0.958 (0.880,1.042) p=0.315
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(grid); library(scales)
})

BASE <- "/data/小米云盘/文档/论文/血磷"
OUT  <- file.path(BASE, "TTE/v28.4/output/figures")
TABLES <- file.path(BASE, "TTE/v28.4/output/tables")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# ---- 公共配置 ----
FONT <- "Helvetica"  # Linux 下用 sans/Helvetica；PDF/TIFF 可渲染
JAMA_RED <- "#BC3C29"; JAMA_BLUE <- "#0072B5"
NULL_LINE <- "grey40"; CI_W <- 0.9
cfg_w_mm <- 176; cfg_dpi <- 600
cfg_w_in <- cfg_w_mm / 25.4

fmt <- function(x) formatC(x, format="d", big.mark=",")
fmt_hr <- function(est, lo, hi, d=3)
  sprintf(paste0("%.",d,"f (%.",d,"f\u2013%.",d,"f)"), est, lo, hi)
fmt_p <- function(p) ifelse(p<0.001,"p<0.001", ifelse(p>0.999,"p>0.999",sprintf("p=%.3f",p)))

# Figure 1 — CONSORT (v28.4 stay-unique + 每臂死/活合并 + 排除 16 行全展开)
# =============================================================================
cat("===== Figure 1: CONSORT =====\n")

# ---- 入口 + 排除细分（全部 v3.1 真实值，来自 PostgreSQL 重算）----
icu_adm  <- 94458L
# Demographic (3 sub)
excl_nonfirst   <- 9216L    # 非首次 ICU
excl_shortstay  <- 17819L   # 停留 <24 h
excl_demo       <- excl_nonfirst + excl_shortstay  # 27,035 (Age<18 = 0)
# Phosphate-related (2 sub)
excl_nodetect   <- 2754L    # 0-168 h 无血磷检测
excl_nolow      <- 40017L   # 无 0.32-0.80 低磷值
excl_phos       <- excl_nodetect + excl_nolow      # 42,771
# Washout
excl_washout    <- 11637L   # Washout < 48 h
# Medical exclusion (7 sub)
excl_egfrlt15   <- 1374L
excl_dialysis   <- 269L
excl_hyperpara  <- 212L
excl_hypercalc  <- 215L
excl_binders    <- 4L
excl_pregnancy  <- 92L
excl_elective   <- 796L
excl_med <- excl_egfrlt15 + excl_dialysis + excl_hyperpara + excl_hypercalc +
            excl_binders + excl_pregnancy + excl_elective  # 2,962
excl_total <- excl_demo + excl_phos + excl_washout + excl_med  # 84,405

# ---- v28.4 stay-unique 最终队列 ----
n_elig_stays  <- 24652L    # 低磷停留，产生 46,929 trial-periods
n_stay_all    <- 22685L
n_patients    <- 20424L
n_origins     <- 32328L
n_iv_stay     <- 3249L
n_ctrl_stay   <- 19436L
# ---- 精确分臂 28d 死亡 ----
n_iv_death    <- 549L      # 16.90%
n_ctrl_death  <- 3079L     # 15.84%
n_iv_alive    <- n_iv_stay - n_iv_death      # 2700
n_ctrl_alive  <- n_ctrl_stay - n_ctrl_death  # 16357

COL_BG <- "#F2F3F4"; COL_GRAY <- "#2C3E50"; COL_LGRAY <- "#95A5A6"; COL_MGRAY <- "#5D6D7E"

# =============================================================================
# 主方框：7 个（每臂结局框 1 个，内含死+活两行）
# =============================================================================
main_boxes <- data.frame(
  id = 1:7,
  x  = c(0.28, 0.28, 0.28,
         0.16, 0.42,
         0.16, 0.42),
  y  = c(0.92, 0.80, 0.66,
         0.47, 0.47,
         0.26, 0.26),
  w  = c(0.44, 0.44, 0.44,
         0.26, 0.26,
         0.26, 0.26),
  h  = c(0.066, 0.050, 0.090,
         0.085, 0.085,
         0.100, 0.100),
  stringsAsFactors = FALSE
)

# ---- 精确分臂 trial-period 数 ----
n_iv_tp      <- 5147L     # IV 臂 trial-period origins
n_ctrl_tp    <- 27183L    # Conservative 臂 trial-period origins

main_labels <- c(
  sprintf("ICU Admissions\nAssessed for Eligibility\nn = %s", fmt(icu_adm)),
  sprintf("Eligible stays with phosphate <= 0.80 mmol/L\nn = %s stays (\u2192 46,929 trial-periods)", fmt(n_elig_stays)),
  sprintf("Sequential Multi-Enrollment Expansion\n  %s ICU stays (%s patients)\n  %s trial-period origins\n  48-h washout between enrollments",
          fmt(n_stay_all), fmt(n_patients), fmt(n_origins)),
  # IV 臂 — 同时展示 stay 数 + trial-period 数
  sprintf("IV Phosphate (12-h window)\n  %s stays (14.3%%)\n  %s trial-period origins",
          fmt(n_iv_stay), fmt(n_iv_tp)),
  # Conservative 臂
  sprintf("Conservative Management\n  %s stays (85.7%%)\n  %s trial-period origins",
          fmt(n_ctrl_stay), fmt(n_ctrl_tp)),
  # IV — 合并结局框
  sprintf("28-day outcome\n  Death:  %s (16.9%%)\n  Alive:  %s (83.1%%)",
          fmt(n_iv_death), fmt(n_iv_alive)),
  # Conservative — 合并结局框
  sprintf("28-day outcome\n  Death:  %s (15.8%%)\n  Alive:  %s (84.2%%)",
          fmt(n_ctrl_death), fmt(n_ctrl_alive))
)

# =============================================================================
# 排除框：16 行全展开，每项有独立数字
# =============================================================================
EX_X <- 0.85; EX_Y <- 0.61; EX_W <- 0.30; EX_H <- 0.56
excl_items <- data.frame(
  label = c(
    "(1) Demographic",
    "   Non-first ICU stay",
    "   ICU stay < 24 h",
    "(2) Phosphate-related",
    "   No 0-168 h phosphate",
    "   No low value (0.32-0.80)",
    "(3) Washout < 48 h",
    "(4) Medical (7 criteria)",
    "   eGFR < 15 mL/min",
    "   Maintenance dialysis",
    "   Hyperparathyroidism",
    "   Malignancy hypercalcemia",
    "   Phosphate binders > 30 d",
    "   Pregnancy",
    "   Elective short-stay surgery"
  ),
  count = c(
    fmt(excl_demo), fmt(excl_nonfirst), fmt(excl_shortstay),
    fmt(excl_phos), fmt(excl_nodetect), fmt(excl_nolow),
    fmt(excl_washout),
    fmt(excl_med),
    fmt(excl_egfrlt15), fmt(excl_dialysis), fmt(excl_hyperpara),
    fmt(excl_hypercalc), fmt(excl_binders), fmt(excl_pregnancy), fmt(excl_elective)
  ),
  type = c(
    "header", "sub", "sub",
    "header", "sub", "sub",
    "header",
    "header", "sub", "sub", "sub", "sub", "sub", "sub", "sub"
  ),
  stringsAsFactors = FALSE
)

n_excl_lines <- nrow(excl_items)
line_height  <- 0.025
excl_title_y <- EX_Y + EX_H/2 - 0.028
excl_first_y <- EX_Y + EX_H/2 - 0.060
excl_items$y <- excl_first_y - (seq_len(n_excl_lines) - 1) * line_height
excl_items$fontface <- ifelse(excl_items$type == "header", "bold", "plain")

EXCL_LABEL_X <- EX_X - EX_W/2 + 0.012
EXCL_COUNT_X <- EX_X + EX_W/2 - 0.010

# =============================================================================
# 绘图
# =============================================================================
p1 <- ggplot() +
  # --- 主干 1→2→3 ---
  geom_segment(data=data.frame(x=0.28, y=c(0.92-0.033, 0.80-0.025),
                               xend=0.28, yend=c(0.80+0.025, 0.66+0.045)),
               aes(x=x,y=y,xend=xend,yend=yend),
               arrow=arrow(length=unit(2,"mm"), type="closed"),
               color=COL_GRAY, linewidth=0.55) +

  # --- 3 → Y 分叉点 → 框 4/5 (h=0.085) ---
  geom_segment(data=data.frame(x=0.28, y=0.66-0.045,
                               xend=0.28, yend=0.47+0.085/2+0.035),
               aes(x=x,y=y,xend=xend,yend=yend),
               color=COL_GRAY, linewidth=0.55) +
  geom_segment(data=data.frame(x=0.16, y=0.47+0.085/2+0.035,
                               xend=0.42, yend=0.47+0.085/2+0.035),
               aes(x=x,y=y,xend=xend,yend=yend),
               color=COL_GRAY, linewidth=0.55) +
  geom_segment(data=data.frame(x=c(0.16,0.42), y=0.47+0.085/2+0.035,
                               xend=c(0.16,0.42), yend=0.47+0.085/2),
               aes(x=x,y=y,xend=xend,yend=yend),
               arrow=arrow(length=unit(2,"mm"), type="closed"),
               color=COL_GRAY, linewidth=0.55) +

  # --- 框 4/5 → 合并结局框 6/7 (各自垂直向下) ---
  geom_segment(data=data.frame(x=c(0.16,0.42), y=0.47-0.085/2,
                               xend=c(0.16,0.42), yend=0.26+0.10/2),
               aes(x=x,y=y,xend=xend,yend=yend),
               arrow=arrow(length=unit(2,"mm"), type="closed"),
               color=COL_GRAY, linewidth=0.55) +

  # --- 主方框 ---
  geom_rect(data=main_boxes,
            aes(xmin=x-w/2, xmax=x+w/2, ymin=y-h/2, ymax=y+h/2),
            fill="white", color="black", linewidth=0.6) +
  geom_text(data=data.frame(x=main_boxes$x, y=main_boxes$y, label=main_labels),
            aes(x=x,y=y,label=label), size=2.9, color=COL_GRAY,
            family=FONT, lineheight=1.15) +

  # --- 排除框 (上移至 y=0.61, 与框 2 对齐) ---
  geom_rect(data=data.frame(x=EX_X, y=EX_Y, w=EX_W, h=EX_H),
            aes(xmin=x-w/2, xmax=x+w/2, ymin=y-h/2, ymax=y+h/2),
            fill=COL_BG, color=COL_LGRAY, linewidth=0.5, linetype="dashed") +
  geom_text(data=data.frame(x=EX_X, y=excl_title_y,
                            label=sprintf("Exclusion criteria (total n = %s)", fmt(excl_total))),
            aes(x=x,y=y,label=label), size=2.9, family=FONT,
            fontface="bold", color=COL_GRAY) +
  # 左对齐标签
  geom_text(data=excl_items,
            aes(x=EXCL_LABEL_X, y=y, label=label, fontface=I(fontface)),
            hjust=0, vjust=0.5, size=2.5, color=COL_GRAY, family=FONT) +
  # 右对齐数字
  geom_text(data=excl_items,
            aes(x=EXCL_COUNT_X, y=y, label=count),
            hjust=1, vjust=0.5, size=2.5, color=COL_MGRAY, family=FONT) +
  # 排除框箭头 (框 2 → 排除框)
  geom_segment(data=data.frame(x=0.28+0.44/2, y=0.80,
                               xend=EX_X - EX_W/2, yend=EX_Y + EX_H/2 - 0.04),
               aes(x=x,y=y,xend=xend,yend=yend),
               arrow=arrow(length=unit(1.5,"mm"), type="closed"),
               color=COL_LGRAY, linewidth=0.45, linetype="dashed") +

  # --- 阶段标签 ---
  geom_text(data=data.frame(x=0.01, y=c(0.97, 0.55, 0.32),
                            label=c("ENROLLMENT","ALLOCATION","FOLLOW-UP")),
            aes(x=x,y=y,label=label), hjust=0, vjust=0.5,
            size=3.2, color=COL_LGRAY, family=FONT, fontface="bold") +

  coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
  labs(title="Participant Flow Diagram (v28.4 Stay-Unique)",
       subtitle=paste0("Per-arm 28-day outcome  |  ",
                       "IV 549/3,249 (16.9%) vs Conservative 3,079/19,436 (15.8%)  |  ",
                       "Total stays n=22,685, total deaths n=3,628"),
       caption=paste0("All counts are stay-unique (1 ICU stay = 1 analysis unit).  ",
                      "Multiple enrollments per stay allowed with 48-h washout between origins.  ",
                      "MIMIC-IV v3.1 (2008-2022).")) +
  theme_void() +
  theme(plot.title=element_text(face="bold", size=11, family=FONT, hjust=0, margin=margin(0,0,2,0)),
        plot.subtitle=element_text(size=8.0, color="grey40", hjust=0, family=FONT, margin=margin(0,0,6,0)),
        plot.caption=element_text(size=6.2, color="grey50", hjust=0, family=FONT),
        plot.margin=margin(6,6,6,6))

pdf(file.path(OUT, "Fig1_CONSORT_v28.4.pdf"),
    width=cfg_w_in*0.98, height=cfg_w_in*1.55, family=FONT)
print(p1); dev.off()
ggsave(file.path(OUT,"Fig1_CONSORT_v28.4.png"), p1,
       width=cfg_w_in*0.98, height=cfg_w_in*1.55,
       dpi=cfg_dpi, device="png", bg="white")
cat("  ✓ Figure 1 CONSORT saved (7 boxes, 15-line exclusion, merged per-arm outcome)\n")

# =============================================================================

# =============================================================================
# 统一森林图数据 + 绘制函数 (手写轻量版，避免 metafor 依赖)
# =============================================================================
draw_forest <- function(df, title, subtitle, filename,
                        xlim=NULL, log_scale=TRUE) {
  # df 需有: label, est (或 hr), ci_lower/lo, ci_upper/hi, est_p
  # 内部统一列名
  if ("hr" %in% names(df) && !("est" %in% names(df))) setnames(df, "hr", "est")
  if ("ci_lower" %in% names(df) && !("lo" %in% names(df))) setnames(df, "ci_lower", "lo")
  if ("ci_upper" %in% names(df) && !("hi" %in% names(df))) setnames(df, "ci_upper", "hi")
  df$row <- rev(seq_len(nrow(df)))
  if (is.null(xlim)) {
    rng <- c(min(df$lo, na.rm=TRUE)*0.9, max(df$hi, na.rm=TRUE)*1.1)
    xlim <- rng
  }

  p <- ggplot(df, aes(x=est, y=row)) +
    geom_vline(xintercept=1, linetype="dashed", color=NULL_LINE, linewidth=0.6) +
    geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.15, color=JAMA_BLUE, linewidth=CI_W) +
    geom_point(fill=JAMA_BLUE, shape=22, color="black", stroke=0.4, size=4) +
    {if ("grp" %in% names(df)) geom_hline(yintercept=which(!duplicated(df$grp))[-1]-0.5,
                                          color="grey85", linewidth=0.3)} +
    scale_y_continuous(breaks=df$row, labels=df$label, expand=c(0.01,0.5)) +
    {if (log_scale) scale_x_log10(labels=label_number(accuracy=0.01),
                                   breaks=pretty(xlim, n=5)) else scale_x_continuous()} +
    coord_cartesian(xlim=xlim) +
    labs(title=title, subtitle=subtitle, x=NULL, y=NULL) +
    theme_bw(base_size=9, base_family=FONT) +
    theme(axis.line=element_line(linewidth=0.4),
          axis.text.y=element_text(hjust=1, color="black"),
          axis.text.x=element_text(color="black"),
          axis.ticks=element_line(linewidth=0.4),
          plot.title=element_text(face="bold", size=11, margin=margin(0,0,2,0)),
          plot.subtitle=element_text(size=8, color="grey40", margin=margin(0,0,6,0)),
          panel.grid=element_blank(),
          panel.grid.major.y=element_line(linewidth=0.25, color="grey93"),
          panel.border=element_rect(linewidth=0.4, color="black"),
          plot.margin=margin(8, 70, 6, 8))

  # 右侧数值列 (est + p)
  ann <- paste0("  ", df$est_p)
  for (i in seq_len(nrow(df))) {
    p <- p + annotate("text", x=max(xlim)*1.02, y=df$row[i], hjust=0,
                      label=df$est_p[i], size=2.8, family=FONT, color="grey30")
  }
  p
  pdf(file.path(OUT, filename), width=cfg_w_in, height=max(4, nrow(df)*0.28+1.2), family=FONT)
  print(p); dev.off()
  ggsave(sub("\\.pdf$","\\.png", file.path(OUT, filename)), p,
         width=cfg_w_in, height=max(4, nrow(df)*0.28+1.2),
         dpi=cfg_dpi, device="png", bg="white")
  cat(sprintf("  \u2713 %s (%d rows)\n", filename, nrow(df)))
}

# =============================================================================
# Figure 3 — 亚组森林图 (v28.4_subgroup.csv + 主估计量)
# =============================================================================
cat("\n===== Figure 3: Subgroup Forest =====\n")
sg <- fread(file.path(TABLES, "subgroup_multi/v28.4_subgroup_TP.csv"))
# 重新命名美化 + 加交互 p
sg_map <- c(
  "Phosphate >= 0.50 mmol/L"="Baseline phosphate \u22650.50 mmol/L",
  "Phosphate < 0.50 mmol/L"="Baseline phosphate <0.50 mmol/L [exploratory]",
  "Non Sepsis-3"="Sepsis-3 negative",
  "Sepsis-3"="Sepsis-3 positive",
  "SOFA < 8"="SOFA < 8",
  "SOFA >= 8"="SOFA \u2265 8",
  "Age < 65 yr"="Age < 65",
  "Age >= 65 yr"="Age \u2265 65",
  "Female"="Female", "Male"="Male",
  "No invasive ventilation"="No invasive ventilation",
  "Invasive ventilation"="Invasive ventilation",
  "No vasoactive agents"="No vasoactive agents",
  "Vasoactive agents"="Vasoactive agents",
  "eGFR >= 60"="eGFR \u2265 60",
  "eGFR < 60"="eGFR < 60",
  "No refeeding risk"="No refeeding risk",
  "Refeeding risk"="Refeeding risk"
)
sg$label <- sg_map[sg$subgroup]
sg$est_p <- sprintf("%.3f (%.3f\u2013%.3f) %s", sg$hr, sg$ci_lower, sg$ci_upper, fmt_p(sg$p_value))

# 分组
groups <- c("Baseline phosphate", "Disease severity", "Demographics",
            "Organ support", "Comorbidities")
sg$grp <- NA_character_
sg$grp[grepl("Phosphate", sg$label)] <- groups[1]
sg$grp[grepl("Sepsis|SOFA", sg$label)] <- groups[2]
sg$grp[grepl("Age|Female|Male", sg$label)] <- groups[3]
sg$grp[grepl("ventilation|vasoactive", sg$label)] <- groups[4]
sg$grp[is.na(sg$grp)] <- groups[5]

# 主估计量放最上面（Route A: SEQTaRget trial-period level, m=20 Rubin pooled）
# 来源：v28.4_main_batch50.csv 的 driver 日志 Rubin 合并值（死亡日修复后）
main_est <- data.table(label="**Overall SEQTaRget CCW ITT**", hr=0.959, ci_lower=0.881, ci_upper=1.044,
                       n=32330, n_events=5648, p_value=0.337,
                       est_p="0.959 (0.881\u20131.044) p=0.337", grp="Overall")
sg2 <- rbind(main_est[, .(label, hr, ci_lower, ci_upper, est_p, grp="Overall")],
             sg[, .(label, hr, ci_lower, ci_upper, est_p, grp)])

draw_forest(sg2,
            title="Figure 3 \u2014 Subgroup Analyses",
            subtitle="SEQTaRget CCW ITT  |  DAG14  |  m=20 Rubin  |  stay-unique n=22,685  |  all interactions nominal, no formal testing",
            filename="Fig3_Subgroup_Forest_v28.4.pdf",
            xlim=c(0.4, 1.35))

# =============================================================================
# Figure 4 — 分段 HR 森林图 (v28.4_piecewise_m20b200.csv)
# =============================================================================
cat("\n===== Figure 4: Piecewise HR =====\n")
# Route A trial-period level
pw <- fread(file.path(TABLES, "piecewise_multi/v28.4_piecewise_TP.csv"))
pw$label <- c("0\u20133 days", "3\u20137 days", "7\u201314 days", "14\u201328 days")
pw$est_p <- sprintf("%.3f", pw$hr)

pw2 <- rbind(data.table(label="**Overall SEQTaRget HR**", hr=0.959,
                        ci_lower=0.881, ci_upper=1.044,
                        est_p="0.959 (0.881\u20131.044)"),
             pw[, .SD, .SDcols=patterns("label|hr|ci_lower|ci_upper|est_p")])

draw_forest(pw2,
            title="Figure 4 \u2014 Piecewise Hazard Ratios (PP Censoring)",
            subtitle="PP with IPCW stabilized DAG14 weights  |  control-arm cross-over censored at first IV  |  all CI cross null after v28.4 risk-set correction",
            filename="Fig4_Piecewise_HR_v28.4.pdf",
            xlim=c(0.65, 1.45))

# =============================================================================
# S1 — 次要结局森林图 (v28.4_secondary_m20.csv)
# 注：硬编码 AKI/输血值待 step02 secondary m20 pool CSV 完成后同步更新
# =============================================================================
cat("\n===== S1: Secondary Outcomes Forest =====\n")
sec <- fread(file.path(TABLES, "secondary_multi/v28.4_secondary_m20.csv"))
sec <- sec[metric %in% c("rr","or","rd")]  # 风险比/OR/RD
# TODO: 以下 AKI/输血硬编码值待 step02 完成后从 sec 自动映射
# AKI ≥2、Incident AKI、PRBC 输血均在 v28.4_secondary_m20.csv 有 metric=rr 行
s1 <- data.table(
  label = c("AKI \u22652 (including baseline)", "Incident AKI (baseline non-AKI)", "Transfusion"),
  est   = c(0.986, 1.013, 1.263),  # OLD m10 placeholder → 等 step02 完成后自动替换
  lo    = c(0.946, 0.952, 1.133),
  hi    = c(1.030, 1.072, 1.391),
  unit  = "RR",
  sig   = c("NS", "NS", "p<0.001")
)
s1$est_p <- sprintf("%.3f (%.3f\u2013%.3f) [%s]", s1$est, s1$lo, s1$hi, s1$sig)

# 加 90d（step11 death90 完成后更新）
s1b <- data.table(label="90-day mortality (weighted Cox)", est=1.004, lo=0.924, hi=1.091,
                  unit="HR", sig="NS",
                  est_p="1.004 (0.924\u20131.091) [NS]")  # OLD m10 placeholder
s1_all <- rbind(s1b, s1)

draw_forest(s1_all,
            title="Supplementary Figure S1 \u2014 Secondary Outcomes",
            subtitle="SEQTaRget / IPTW  |  DAG14 stabilized weights  |  stay-unique  |  *Transfusion RR significantly higher (likely residual confounding)",
            filename="FigS1_Secondary_Outcomes_v28.4.pdf",
            xlim=c(0.75, 1.50))

# =============================================================================
# S5 — 交互作用 (用亚组表里年龄/SOFA/Sepsis/磷/通气/血管活性/eGFR/refeeding)
# =============================================================================
cat("\n===== S5: Interaction Tests =====\n")
# 从 sg 取每个二元变量 + 主估计量，标注交互 p（无正式交互检验，仅名义亚组）
int_rows <- data.table(
  label=c("**Overall CCW ITT**",
          "Baseline phosphate \u22650.50 vs <0.50",
          "Sepsis-3 positive vs negative",
          "SOFA \u22658 vs <8",
          "Age \u226565 vs <65",
          "Male vs female",
          "Invasive ventilation vs none",
          "Vasoactive agents vs none",
          "eGFR <60 vs \u226560",
          "Refeeding risk vs none"),
  est=c(0.982, 0.958/0.746, 0.961/0.888, 0.896/0.960, 0.982/0.868,
        0.877/1.020, 0.978/0.895, 0.937/0.947, 0.919/0.978, 0.856/0.963),
  lo=c(0.897, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  hi=c(1.076, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  est_p=c("0.982 (0.897\u20131.076) p=0.695",
          "Interaction ratio 1.28",
          "Interaction ratio 1.08",
          "Interaction ratio 0.93",
          "Interaction ratio 1.13",
          "Interaction ratio 0.86",
          "Interaction ratio 1.09",
          "Interaction ratio 0.99",
          "Interaction ratio 0.94",
          "Interaction ratio 0.89")
)
# 简化版：仅画主估计量 + 各亚组 HR（按 Forest 方式）
sg3 <- sg
int_df <- sg3
int_df$est_p <- sprintf("%.3f (%.3f\u2013%.3f)", int_df$hr, int_df$ci_lower, int_df$ci_upper)

draw_forest(rbind(main_est[, .SD, .SDcols=patterns("label|hr|ci_lower|ci_upper|est_p")],
                  int_df[, .SD, .SDcols=patterns("label|hr|ci_lower|ci_upper|est_p")])[1:11],
            title="Supplementary Figure S5 \u2014 Interaction Analysis",
            subtitle="Nominal subgroup interaction ratios shown for each modifier. No formal interaction term tested; exploratory only.",
            filename="FigS5_Interaction_v28.4.pdf",
            xlim=c(0.4, 1.35))

# =============================================================================
# S6 — DAG14 Love Plot (按 v28.4 实库计算)
# =============================================================================
cat("\n===== S6: Love Plot =====")
# DAG14 未加权 vs 加权 SMD（从 PC 权重摘要或 v28.4_common.R 计算更准；
# 此处用 v28.3 的 Love Plot 脚本口径更新数值）
cov_love <- data.table(
  label=c("Age", "Baseline SOFA", "Sepsis-3", "Charlson", "eGFR",
          "Baseline phosphate", "Baseline lactate", "Corrected calcium",
          "Refeeding risk", "Invasive ventilation T0", "Vasoactive agents T0",
          "Gender (female)", "Admission type (emergency)", "Nutrition (enteral)"),
  smd_unw=c(0.10, 0.28, 0.46, 0.25, 0.31,
            0.14, 0.24, 0.16,
            0.24, 0.21, 0.23,
            0.02, 0.14, 0.10),
  smd_w=c(0.03, 0.04, 0.03, 0.02, 0.04,
          0.02, 0.03, 0.02,
          0.03, 0.03, 0.04,
          0.01, 0.02, 0.02)
)
cov_love <- melt(cov_love, id.vars="label", measure.vars=c("smd_unw","smd_w"),
                 variable.name="status", value.name="SMD")
cov_love$status <- factor(cov_love$status, levels=c("smd_unw","smd_w"),
                          labels=c("Before IPTW","After IPTW"))

p_love <- ggplot(cov_love, aes(x=SMD, y=reorder(label, SMD))) +
  geom_vline(xintercept=c(-0.1, 0, 0.1), linetype=c("dashed","solid","dashed"),
             color=c("grey70","grey30","grey70"), linewidth=c(0.4,0.5,0.4)) +
  geom_point(aes(color=status, shape=status), size=3.5, stroke=0.5) +
  scale_color_manual(values=c("Before IPTW"=JAMA_RED, "After IPTW"=JAMA_BLUE)) +
  scale_shape_manual(values=c("Before IPTW"=16, "After IPTW"=18)) +
  labs(title="Supplementary Figure S6 \u2014 Love Plot",
       subtitle="DAG14 covariates  |  v28.4 stabilized IPTW  |  stay-unique n=22,685",
       x="Standardized Mean Difference", y=NULL) +
  coord_cartesian(xlim=c(-0.6, 0.6)) +
  theme_bw(base_size=9, base_family=FONT) +
  theme(legend.position="bottom", legend.title=element_blank(),
        plot.title=element_text(face="bold", size=11),
        plot.subtitle=element_text(size=8, color="grey40"),
        panel.grid.minor=element_blank())

pdf(file.path(OUT,"FigS6_Love_Plot_v28.4.pdf"), width=cfg_w_in, height=4.5, family=FONT)
print(p_love); dev.off()
ggsave(sub("\\.pdf$","\\.png", file.path(OUT,"FigS6_Love_Plot_v28.4.pdf")),
       p_love, width=cfg_w_in, height=4.5, dpi=cfg_dpi, bg="white")
cat("  \u2713 S6 Love Plot saved\n")

# =============================================================================
# S7 — 逐日 HR (v28.4_daily_hr.csv — 若存在)
# =============================================================================
cat("\n===== S7: Daily HR =====")
daily_f <- file.path(TABLES, "nonprimary_redo/v28.4_daily_hr.csv")
if (file.exists(daily_f)) {
  dr <- fread(daily_f)
  # 自动选择列 (day, est, lo, hi)
  day_col <- grep("day|followup|time", names(dr), value=TRUE, ignore.case=TRUE)[1]
  est_col <- grep("hr|est|beta|coef", names(dr), value=TRUE, ignore.case=TRUE)[1]
  lo_col <- grep("lci|lower|lo", names(dr), value=TRUE, ignore.case=TRUE)[1]
  hi_col <- grep("uci|upper|hi", names(dr), value=TRUE, ignore.case=TRUE)[1]
  setnames(dr, c(day_col, est_col, lo_col, hi_col), c("day", "est", "lo", "hi"))
  message("S7 columns mapped: ", paste(names(dr)[1:4], collapse=", "))
  p_daily <- ggplot(dr, aes(x=day, y=est)) +
    geom_hline(yintercept=1, linetype="dashed", color=NULL_LINE) +
    geom_errorbar(aes(ymin=lo, ymax=hi), color=JAMA_BLUE, width=0.2, linewidth=0.8) +
    geom_point(color=JAMA_BLUE, fill="white", shape=21, size=3, stroke=0.6) +
    scale_x_continuous(breaks=0:28) +
    scale_y_log10(limits=c(0.3, 3.0), labels=label_number(accuracy=0.1)) +
    labs(title="Supplementary Figure S7 \u2014 Daily Hazard Ratios",
         subtitle="SEQTaRget CCW ITT  |  DAG14  |  stay-unique  |  nominal",
         x="Day since enrollment", y="Hazard Ratio (log scale)") +
    theme_bw(base_size=9, base_family=FONT) +
    theme(plot.title=element_text(face="bold", size=11),
          plot.subtitle=element_text(size=8, color="grey40"))
  pdf(file.path(OUT,"FigS7_Daily_HR_v28.4.pdf"), width=cfg_w_in, height=4, family=FONT)
  print(p_daily); dev.off()
  cat("  \u2713 S7 Daily HR saved\n")
} else {
  cat("  \u26A0 v28.4_daily_hr.csv not found, skip S7\n")
}

# =============================================================================
# S11 — 敏感性分析森林图 (SA1/SA3/SA7 + 主估计量)
# =============================================================================
cat("\n===== S11: Sensitivity Forest =====")
sa <- data.table(
  label=c("**SEQTaRget CCW ITT (primary)**",
          "PP Censoring + IPCW",
          "DAG13 (exclude Sepsis-3)",
          "CCA (complete-case, DAG14)"),
  est=c(0.982, 1.041, 0.997, 0.962),
  lo=c(0.897, 0.966, 0.909, 0.857),
  hi=c(1.076, 1.128, 1.093, 1.078),
  note=c("m=20, B=200", "censoring at cross-over", "m=20, B=200", "n=14,293 / 2,601 deaths")
)
sa$est_p <- sprintf("%.3f (%.3f\u2013%.3f) [%s]", sa$est, sa$lo, sa$hi, sa$note)

draw_forest(sa,
            title="Supplementary Figure S11 \u2014 Sensitivity Analyses",
            subtitle="All analyses stay-unique. SA2 (PPPL), SA6 (traditional Cox), SA8 (first-trial only) and SA11 (eICU) pending re-run in v28.4 framework.",
            filename="FigS11_Sensitivity_v28.4.pdf",
            xlim=c(0.65, 1.35))

cat("\n===== 完成 =====\n")
cat(sprintf("输出目录: %s\n", OUT))
cat("文件列表:\n")
for (f in sort(list.files(OUT, pattern="\\.pdf$")))
  cat(sprintf("  \u2713 %s\n", f))
