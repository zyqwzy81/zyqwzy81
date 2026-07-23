# ============================================================================
# V26.9 发表级图表生成脚本
# 遵循规范八：Critical Care Medicine 杂志标准
# 生成：Fig 2-8 + Fig S1-S6（从已有表格数据）
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(scales)
})

# ---- 输出目录 ----
OUT_DIR  <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9"
FIG_DIR  <- file.path(OUT_DIR, "figures")
TBL_DIR  <- file.path(OUT_DIR, "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- JAMA 配色 ----
JAMA_RED    <- "#BC3C29"
JAMA_BLUE   <- "#0072B5"
FOREST_GREEN <- "#2E8B57"
DEEP_NAVY   <- "#1B4F72"
EVAL_WEAK   <- "#F39C12"
EVAL_MOD    <- "#E67E22"
EVAL_STRONG <- "#27AE60"

save_fig <- function(p, name, w = 10, h = 6) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf")
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300, device = "png")
  cat(sprintf("  ✓ %s\n", name))
}

# ============================================================================
# Fig 2: Primary and Sensitivity Analysis Forest Plot
# ============================================================================
cat("Fig 2: Main Forest Plot...\n")

rubin <- fread(file.path(TBL_DIR, "v26_9_seqtarget_rubin_results.csv"))
fig2_data <- data.table(
  label = c("Primary ITT (DAG14, m=20, boot=200)",
            "DAG13 Sensitivity (excl. sepsis3)",
            "IPTW-Cox (Rubin m=5)",
            "Weibull AFT (Time Ratio)",
            "Per-Protocol CCW+IPTW (boot=50)",
            "First-Enrollment-Only (n=22,685)",
            "eICU External Validation (11 cov)"),
  hr = c(1.054, 1.063, 0.978, 0.989, 1.044, 1.047, 0.971),
  ci_lower = c(0.954, 0.960, 0.925, 0.940, 0.956, 0.951, 0.859),
  ci_upper = c(1.165, 1.176, 1.034, 1.041, 1.139, 1.154, 1.097),
  p = c("0.304", "0.260", "0.438", "0.672", "0.340", "0.351", "0.628"),
  n = c("32,325", "32,325", "32,325", "32,325", "32,325", "22,685", "34,882"),
  events = c("5,610", "5,610", "5,610", "5,610", "5,610", "5,610", "5,508"),
  stringsAsFactors = FALSE
)
fig2_data$y <- 7:1
fig2_data$is_primary <- c(TRUE, rep(FALSE, 6))

p2 <- ggplot(fig2_data, aes(x = hr, y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.15, linewidth = 0.9, color = DEEP_NAVY) +
  geom_point(aes(fill = is_primary), size = 4.5, shape = 22, color = "black") +
  scale_fill_manual(values = c("TRUE" = JAMA_RED, "FALSE" = JAMA_BLUE), guide = "none") +
  scale_y_continuous(breaks = 7:1, labels = rev(fig2_data$label)) +
  scale_x_log10(breaks = c(0.80, 0.90, 1.00, 1.10, 1.20),
                labels = c("0.80", "0.90", "1.00", "1.10", "1.20")) +
  labs(title = "Figure 2. Primary and Sensitivity Analysis Forest Plot",
       subtitle = "SEQTaRget ITT framework · DAG14 covariates · Cloning-Censoring-Weighting",
       x = "Hazard Ratio (95% CI)", y = "",
       caption = "Primary ITT: Rubin's rules across m=20 MICE imputed datasets, bootstrap=200") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(size = 9, color = "gray30", hjust = 0),
        axis.text.y = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(10, 10, 10, 10)) +
  coord_cartesian(xlim = c(0.75, 1.30))

# Add annotation columns
ann_x <- 1.35
for (i in 1:7) {
  p2 <- p2 + annotate("text", x = ann_x, y = i, 
                       label = sprintf("%.3f (%.3f-%.3f)", fig2_data$hr[i], fig2_data$ci_lower[i], fig2_data$ci_upper[i]),
                       hjust = 0, size = 3.2, color = ifelse(fig2_data$is_primary[i], JAMA_RED, "black"))
}
save_fig(p2, "Fig2_Main_Forest", 14, 6)

# ============================================================================
# Fig 5: Piecewise Hazard Ratios
# ============================================================================
cat("Fig 5: Piecewise HR...\n")

seg <- fread(file.path(TBL_DIR, "v26_9_secondary_time_segments.csv"))
seg[, label := sprintf("%d-%d days", t_start, t_end)]
seg[, y := .N:1]

p5 <- ggplot(seg, aes(x = estimate, y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.15, linewidth = 1.0, color = DEEP_NAVY) +
  geom_point(size = 5, shape = 18, color = JAMA_BLUE) +
  scale_y_continuous(breaks = seg$y, labels = seg$label) +
  scale_x_log10(breaks = c(0.70, 0.85, 1.00, 1.15, 1.30, 1.50),
                labels = c("0.70", "0.85", "1.00", "1.15", "1.30", "1.50")) +
  labs(title = "Figure 5. Piecewise Hazard Ratios by Follow-up Interval",
       subtitle = "Cox regression with DAG14 covariates · cluster-robust standard errors",
       x = "Hazard Ratio (95% CI)", y = "") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(size = 9, color = "gray30", hjust = 0),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(10, 10, 10, 10)) +
  coord_cartesian(xlim = c(0.65, 1.65))

# Add annotation
for (i in 1:nrow(seg)) {
  p5 <- p5 + annotate("text", x = 1.70, y = i,
                       label = sprintf("HR=%.3f (%.3f-%.3f)", seg$estimate[i], seg$ci_lower[i], seg$ci_upper[i]),
                       hjust = 0, size = 3.5)
}
save_fig(p5, "Fig5_Piecewise_HR")

# ============================================================================
# Fig 6: CCW Period Characteristics
# ============================================================================
cat("Fig 6: CCW Period Characteristics...\n")

per <- fread(file.path(TBL_DIR, "v26_9_ccw_period_characteristics.csv"))
per[, period_label := sprintf("Period %d", period_int)]

p6a <- ggplot(per, aes(x = period_label, y = n_obs)) +
  geom_col(fill = JAMA_BLUE, alpha = 0.8, width = 0.6) +
  geom_text(aes(label = n_obs), vjust = -0.5, size = 3.5) +
  labs(title = "Figure 6. CCW Trial Period Characteristics",
       subtitle = sprintf("Total: %d trial-periods · %d unique patients",
                          sum(per$n_obs), sum(per$n_unique)),
       x = "Sequential Trial Period", y = "Number of Trial-Periods") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(size = 9, color = "gray30", hjust = 0))

p6b <- ggplot(per, aes(x = period_label, y = event_rate, group = 1)) +
  geom_line(color = JAMA_RED, linewidth = 1) +
  geom_point(color = JAMA_RED, size = 3) +
  geom_text(aes(label = sprintf("%.1f%%", event_rate)), vjust = -1, size = 3.2) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "Sequential Trial Period", y = "Event Rate (%)") +
  theme_bw() +
  theme(plot.title = element_blank())

library(patchwork)
p6 <- p6a / p6b + plot_layout(heights = c(2, 1))
save_fig(p6, "Fig6_CCW_Diagnostics", 10, 8)

# ============================================================================
# Fig 7: E-value Sensitivity Plot
# ============================================================================
cat("Fig 7: E-value...\n")

evalue_data <- data.table(
  Analysis = c("Point Estimate\nE-value = 1.26",
               "CI Upper Limit\nE-value = 1.54"),
  E_value = c(1.26, 1.54),
  tier = c("Weak", "Moderate")
)

p7 <- ggplot(evalue_data, aes(x = Analysis, y = E_value, fill = tier)) +
  geom_col(width = 0.5, alpha = 0.85) +
  geom_text(aes(label = sprintf("E = %.2f", E_value)), vjust = -0.5, size = 5) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = EVAL_MOD, linewidth = 0.8) +
  annotate("text", x = 2.5, y = 1.5, label = "Moderate threshold (E=1.5)", 
           hjust = 0, size = 3.5, color = EVAL_MOD) +
  annotate("text", x = 2.5, y = 1.0, label = "Null line (E=1.0)", 
           hjust = 0, size = 3.5, color = "gray40") +
  scale_fill_manual(values = c("Weak" = EVAL_WEAK, "Moderate" = EVAL_MOD), name = "Robustness") +
  labs(title = "Figure 7. E-value Sensitivity Analysis",
       subtitle = "Point estimate E=1.26 (weak) · UCL E=1.54 (moderate)\nA moderately strong unmeasured confounder could explain the null result",
       x = "", y = "E-value") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(size = 9, color = "gray30", hjust = 0),
        legend.position = "bottom") +
  ylim(0, 2.5)
save_fig(p7, "Fig7_Evalue_Curve")

# ============================================================================
# Fig 8: Interaction Forest Plot (with FDR q-values)
# ============================================================================
cat("Fig 8: Interaction Forest...\n")

int <- fread(file.path(TBL_DIR, "v26_9_interaction_fdr.csv"))
int[, label := gsub("_", " ", var)]
int[, label := factor(label, levels = rev(label))]
int[, sig := fifelse(significance_label == "FDR-significant", "FDR q<0.10", "NS")]

p8 <- ggplot(int, aes(x = hr, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.2, linewidth = 0.9) +
  geom_point(aes(color = sig), size = 4) +
  scale_color_manual(values = c("FDR q<0.10" = JAMA_RED, "NS" = "gray50"), name = "") +
  scale_x_log10() +
  labs(title = "Figure 8. Interaction Test Forest Plot",
       subtitle = "Cox proportional hazards · DAG14 covariates · cluster-robust SE · Benjamini-Hochberg FDR",
       x = "Interaction Hazard Ratio (95% CI)", y = "") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(size = 9, color = "gray30", hjust = 0),
        legend.position = "bottom",
        axis.text.y = element_text(size = 10))

# Add q-value labels
ann_x <- max(int$ci_upper) * 1.1
for (i in 1:nrow(int)) {
  p8 <- p8 + annotate("text", x = ann_x, y = int$label[i],
                       label = sprintf("q=%.3f", int$q_value[i]),
                       hjust = 0, size = 3.2,
                       color = ifelse(int$sig[i] == "FDR q<0.10", JAMA_RED, "gray40"))
}
save_fig(p8, "Fig8_Interaction_Forest", 12, 6)

# ============================================================================
# Fig S3: PO vs Non-PO Stratification Forest
# ============================================================================
cat("Fig S3: PO Stratification Forest...\n")

po_data <- data.table(
  stratum = c("Oral Phosphate", "Non-PO (IV or none)", "Primary ITT (reference)"),
  hr = c(1.075, 1.053, 1.054),
  ci_lower = c(0.873, 0.931, 0.954),
  ci_upper = c(1.325, 1.190, 1.165),
  n = c("8,848", "23,445", "32,325"),
  events = c("1,329", "4,275", "5,610")
)
po_data$y <- 3:1
po_data$is_ref <- c(FALSE, FALSE, TRUE)

ps3 <- ggplot(po_data, aes(x = hr, y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.15, linewidth = 0.9, color = DEEP_NAVY) +
  geom_point(aes(color = is_ref), size = 4.5, shape = 15) +
  scale_color_manual(values = c("TRUE" = FOREST_GREEN, "FALSE" = JAMA_BLUE), guide = "none") +
  scale_y_continuous(breaks = 3:1, labels = po_data$stratum) +
  scale_x_log10(breaks = c(0.80, 0.90, 1.00, 1.10, 1.20, 1.40),
                labels = c("0.80", "0.90", "1.00", "1.10", "1.20", "1.40")) +
  labs(title = "Figure S3. PO vs. Non-PO Phosphate Stratification",
       subtitle = "SEQTaRget ITT, bootstrap = 50, DAG14 covariates",
       x = "Hazard Ratio (95% CI)", y = "") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(size = 9, color = "gray30", hjust = 0),
        axis.text.y = element_text(size = 10),
        panel.grid.major.y = element_blank()) +
  coord_cartesian(xlim = c(0.75, 1.50))
save_fig(ps3, "FigS3_PO_Stratification_Forest", 10, 4.5)

# ============================================================================
# Fig S5: E-value Curve (supplementary)
# ============================================================================
cat("Fig S5: E-value Curve (supplementary)...\n")

ps5 <- p7 + labs(title = "Figure S5. E-value Sensitivity Analysis (Supplementary)")
save_fig(ps5, "FigS5_Evalue_Curve")

# ============================================================================
# Fig S6: Trial Period Overview
# ============================================================================
cat("Fig S6: Trial Period Overview...\n")

per[, treatment_rate := n_tx / n_obs * 100]

ps6a <- ggplot(per, aes(x = period_label, y = n_obs)) +
  geom_col(fill = JAMA_BLUE, alpha = 0.7, width = 0.6) +
  geom_text(aes(label = n_obs), vjust = -0.5, size = 3.2) +
  labs(title = "A. Enrollment per Period", x = "", y = "N") +
  theme_bw()

ps6b <- ggplot(per, aes(x = period_label, y = n_events, group = 1)) +
  geom_col(fill = JAMA_RED, alpha = 0.7, width = 0.6) +
  geom_line(aes(y = event_rate * 100), color = DEEP_NAVY, linewidth = 0.8) +
  geom_point(aes(y = event_rate * 100), color = DEEP_NAVY, size = 2) +
  scale_y_continuous(sec.axis = sec_axis(~./100, name = "Mortality Rate (%)")) +
  labs(title = "B. Events and Mortality Rate", x = "", y = "Events") +
  theme_bw()

ps6c <- ggplot(per, aes(x = period_label, y = treatment_rate, group = 1)) +
  geom_line(color = FOREST_GREEN, linewidth = 0.8) +
  geom_point(color = FOREST_GREEN, size = 2.5) +
  geom_hline(yintercept = mean(per$treatment_rate), linetype = "dashed", color = "gray50") +
  labs(title = "C. Treatment Rate", x = "Sequential Trial Period", y = "Treatment (%)",
       caption = "Dashed line: global mean") +
  theme_bw()

library(patchwork)
ps6 <- ps6a / ps6b / ps6c + plot_layout(heights = c(1, 1, 1))
ps6 <- ps6 + plot_annotation(
  title = "Figure S6. Sequential Trial Period Overview",
  subtitle = sprintf("N = %d trial-periods across %d periods", sum(per$n_obs), nrow(per)),
  theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
                plot.subtitle = element_text(size = 9, color = "gray30", hjust = 0))
)
save_fig(ps6, "FigS6_Trial_Period_Overview", 10, 9)

# ============================================================================
# Secondary Outcomes Forest Plot
# ============================================================================
cat("Secondary outcomes forest plot...\n")

sec <- fread(file.path(TBL_DIR, "v26_9_secondary_time_results.csv"))
sec[, label := outcome]
sec <- sec[!grepl("28天死亡", label)]
sec[, y := .N:1]
sec[, sig := fifelse(p_value < 0.01, "Bonferroni sig.", "NS")]

p_sec <- ggplot(sec, aes(x = estimate, y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.15, linewidth = 0.9) +
  geom_point(aes(color = sig), size = 4) +
  scale_color_manual(values = c("Bonferroni sig." = JAMA_RED, "NS" = "gray50"), name = "") +
  scale_y_continuous(breaks = sec$y, labels = sec$label) +
  scale_x_log10() +
  labs(title = "Secondary Outcomes Forest Plot",
       subtitle = "Logistic regression · DAG14 covariates · Bonferroni α=0.01",
       x = "Odds Ratio (95% CI)", y = "") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(size = 9, color = "gray30", hjust = 0),
        legend.position = "bottom")
save_fig(p_sec, "v26_9_secondary_forest")

# ============================================================================
# Death Segmented Forest (supplementary)
# ============================================================================
cat("Death segmented forest plot...\n")

seg[, label2 := sprintf("Days %d-%d", t_start, t_end)]
seg[, sig_lab := fifelse(p_value < 0.001, "p<0.001", sprintf("p=%.3f", p_value))]

p_death <- ggplot(seg, aes(x = estimate, y = label2)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.15, linewidth = 1.0, color = DEEP_NAVY) +
  geom_point(size = 5, shape = 18, color = JAMA_BLUE) +
  scale_x_log10() +
  labs(title = "Death Segmented Hazard Ratios",
       subtitle = "Piecewise Cox · DAG14 covariates · cluster-robust SE",
       x = "Hazard Ratio (95% CI)", y = "") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(size = 9, color = "gray30", hjust = 0))
save_fig(p_death, "v26_9_death_segmented_forest")

cat("\n=== All figures generated successfully ===\n")
cat(sprintf("Output: %s\n", FIG_DIR))
