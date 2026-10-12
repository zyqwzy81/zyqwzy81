# =====================================================================
# PC10_figures.R — Supplementary Figures S-PC1..S-PC4
# Post-hoc IV phosphate vs completely-untreated mechanistic analysis
# Publication style: Okabe-Ito palette, 300 dpi PNG + vector PDF
# =====================================================================
suppressMessages({library(ggplot2);library(dplyr);library(tidyr);library(patchwork)})

OUT <- "/data/小米云盘/文档/论文/血磷/TTE/v28.4/output/pure_control"
FIG <- file.path(OUT, "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) read.csv(file.path(OUT, f), check.names = FALSE)

# Okabe-Ito
OI <- c(A_IV = "#D55E00", B_PO_only = "#009E73", C_outside_window = "#E69F00",
        C_rx_only_uninfused = "#CC79A7", D_none = "#0072B2")
GL <- c(IV = "#D55E00", none = "#0072B2")
base <- 9
tt <- theme_classic(base_size = base) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text = element_text(size = base), axis.title = element_text(size = base),
        plot.title = element_text(size = base + 1, face = "bold"),
        strip.text = element_text(size = base), panel.spacing = unit(0.4, "cm"))

savefig <- function(p, name, w = 180, h = 95) {
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300)
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf)
}

# =====================================================================
# S-PC1 : three-group 72-h trajectory, absolute + delta-from-baseline
# =====================================================================
three <- c("A_IV","B_PO_only","D_none")
traj <- rd("pc_crude_trajectory.csv") %>% filter(exposure_w48 %in% three)
ph1  <- rd("pc_plot_1h.csv") %>% filter(exposure_w48 %in% three, h <= 72)
gee  <- rd("pc_gee_contrasts.csv") %>% filter(comparison == "C1_IV_vs_none")
ex   <- rd("pc_gee_extra.csv")
w12d <- ex %>% filter(comparison == "C1_W12_delta")
lvl <- c(A_IV = "IV replacement", B_PO_only = "Oral/enteral only",
         C_outside_window = "Late IV (>48 h)", C_rx_only_uninfused = "Prescribed, not infused",
         D_none = "No replacement")

cap1 <- "Points = window means (coverage at 6/12/24/48/72 h, IV: 38.5/48.6/82.2/70.5/60.8%; no replacement: 18.3/26.3/60.6/46.2/36.2%).\nCurves = loess-smoothed hourly means; diamonds = IPTW-IPCW weighted GEE marginal means with 95% CI."
pA <- ggplot() +
  geom_hline(yintercept = 0.80, lty = 2, colour = "grey50", linewidth = .4) +
  annotate("text", 70, 0.775, label = "reference 0.80 mmol/L", hjust = 1, size = 2.4, colour = "grey30") +
  geom_smooth(data = ph1, aes(h, mean, colour = exposure_w48), se = FALSE,
              method = "loess", span = 0.35, alpha = .55, linewidth = .45, show.legend = FALSE) +
  geom_line(data = traj, aes(time, mean, colour = exposure_w48), linewidth = .7) +
  geom_point(data = traj, aes(time, mean, colour = exposure_w48), size = 1.6) +
  # IPTW-IPCW weighted GEE estimates, IV vs none
  geom_errorbar(data = gee, aes(time, est1, ymin = est1 - 1.96*se1, ymax = est1 + 1.96*se1),
                colour = OI["A_IV"], width = 2.2, linewidth = .5) +
  geom_point(data = gee, aes(time, est1), colour = OI["A_IV"], shape = 23, size = 2, fill = "white") +
  geom_errorbar(data = gee, aes(time, est0, ymin = est0 - 1.96*se0, ymax = est0 + 1.96*se0),
                colour = OI["D_none"], width = 2.2, linewidth = .5) +
  geom_point(data = gee, aes(time, est0), colour = OI["D_none"], shape = 23, size = 2, fill = "white") +
  annotate("text", 12, 1.13, label = "post-dose peak 4-12 h: +0.35", size = 2.5, colour = OI["A_IV"]) +
  scale_colour_manual(values = OI, labels = lvl) +
  scale_x_continuous(breaks = c(0,6,12,24,48,72)) +
  coord_cartesian(ylim = c(0.4, 1.16), xlim = c(0,72)) +
  labs(x = "Hours after eligibility (T0)", y = "Serum phosphate (mmol/L)",
       title = "A. Absolute concentration",
       subtitle = "diamonds = IPTW-IPCW weighted GEE (IV vs no replacement)") + tt +
  theme(plot.subtitle = element_text(size = base - 1.5, face = "italic", colour = "grey30"))

pB <- ggplot() +
  geom_hline(yintercept = 0, lty = 1, colour = "grey60", linewidth = .4) +
  geom_line(data = traj, aes(time, mean_delta, colour = exposure_w48), linewidth = .7, show.legend = FALSE) +
  geom_point(data = traj, aes(time, mean_delta, colour = exposure_w48), size = 1.6, show.legend = FALSE) +
  geom_errorbar(data = w12d, aes(time, diff, ymin = lcl, ymax = ucl),
                colour = "black", width = 2.2, linewidth = .5, lty = 1) +
  geom_point(data = w12d, aes(time, diff), shape = 4, size = 1.7, stroke = .7) +
  annotate("text", 14.5, 0.105, label = "weighted IV-none\n+0.054 (0.022 to 0.086) at 12 h, early treated",
           hjust = 0, size = 2.5, colour = "black") +
  scale_colour_manual(values = OI, labels = lvl) +
  scale_x_continuous(breaks = c(0,6,12,24,48,72)) +
  coord_cartesian(ylim = c(-0.16, 0.36), xlim = c(0,72)) +
  labs(x = "Hours after eligibility (T0)", y = "Change from baseline (mmol/L)",
       title = "B. Change from baseline",
       subtitle = "x = weighted difference (early-treatment population)") + tt +
  theme(plot.subtitle = element_text(size = base - 1.5, face = "italic", colour = "grey30"),
        legend.position = "none")

savefig(wrap_plots(pA, pB, ncol = 1) + plot_layout(guides = "collect") +
          plot_annotation(caption = cap1, theme = theme(
            legend.position = "bottom",
            plot.caption = element_text(size = base - 2.5, face = "italic",
                                        colour = "grey30", hjust = 0),
            plot.margin = margin(t = 2, r = 6, b = 8, l = 6))),
        "FigS_PC1_trajectory", 180, 182)

# =====================================================================
# S-PC2 : six-domain subgroup forest (early pharmacologic vs W48 net)
# =====================================================================
sg <- rd("pc_subgroups.csv")
varlab <- c(s_baseline_p = "Baseline phosphate", s_renal = "Renal function",
            s_sofa = "SOFA", s_sepsis = "Sepsis-3", s_nutrition = "Nutrition support",
            s_vent = "Mechanical ventilation")
levlab <- c("0.50-0.65" = "0.50-0.65", "0.65-0.80" = "0.65-0.80", "<0.50" = "<0.50",
            "any_nutrition" = "any", "none" = "none",
            "RRT" = "RRT", "eGFR30-60" = "eGFR 30-60", "eGFR<30" = "eGFR <30",
            "eGFR>=60" = "eGFR >=60", "no_sepsis" = "no sepsis", "sepsis" = "sepsis",
            "SOFA6-9" = "SOFA 6-9", "SOFA<6" = "SOFA <6", "SOFA>=10" = "SOFA >=10",
            "invasive_vent" = "invasive ventilation", "no_vent" = "not ventilated")
sg <- sg %>%
  mutate(domain = varlab[stratum_var],
         lev2 = levlab[level],
         analysis = factor(analysis, levels = c("W12_early_12h_delta","W48_48h_delta"),
                           labels = c("Early treated: 12-h signal",
                                      "W48 strategy: 48-h net difference")),
         sig = (!is.na(lcl) & ((lcl > 0) | (ucl < 0)))) %>%
  arrange(analysis, stratum_var, level) %>%
  group_by(analysis) %>% mutate(ord = row_number()) %>% ungroup() %>%
  # stagger adjacent n labels vertically to prevent overlap
  mutate(stagger = ifelse(ord %% 2 == 0, -0.16, 0.16))
lv <- sg %>% filter(analysis == levels(sg$analysis)[1]) %>% arrange(ord) %>%
  {paste0(.$domain, ": ", .$lev2)}
sg$yl <- factor(paste0(sg$domain, ": ", sg$lev2), levels = rev(lv))

p2 <- ggplot(sg, aes(diff, yl)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = .4) +
  geom_linerange(aes(xmin = lcl, xmax = ucl, colour = sig), linewidth = .5,
                 na.rm = TRUE) +
  geom_point(aes(colour = sig), size = 1.8, na.rm = TRUE) +
  geom_text(aes(x = ifelse(diff >= 0, 0.355, -0.615),
                y = as.numeric(yl) + stagger,
                label = paste0("n=", n_trt, "/", n_ctrl)),
            hjust = ifelse(sg$diff >= 0, 1, 0), size = 1.7, colour = "grey30",
            na.rm = TRUE, check_overlap = FALSE) +
  facet_wrap(~analysis) +
  scale_colour_manual(values = c("FALSE" = "grey55", "TRUE" = "#D55E00"), guide = "none") +
  coord_cartesian(xlim = c(-0.66, 0.44), clip = "off") +
  labs(x = "Weighted IV - no-replacement difference in change (mmol/L)",
       y = NULL, title = "Subgroup effects: pharmacodynamic signal confined to less severe patients",
       caption = "Orange = 95% CI excludes 0. n = IV / no-replacement with complete data in stratum") +
  tt + theme(strip.text = element_text(size = base, face = "bold"),
             plot.caption = element_text(size = base - 2, face = "italic", colour = "grey30"),
             plot.margin = margin(t = 4, r = 12, b = 4, l = 12))
savefig(p2, "FigS_PC2_subgroups", 200, 150)

# =====================================================================
# S-PC3 : dose-response RCS + timing x dose (2 x 2)
# =====================================================================
rcs <- rd("pc_dose_rcs.csv") %>% filter(dose <= 90) %>%
  mutate(outcome = factor(outcome, levels = c("d12","d48"),
                          labels = c("12-h change","48-h change")))
q <- rd("pc_dose_quartiles.csv")
td <- rd("pc_timing_dose.csv")
eln <- rd("pc_early_late_none.csv")
none_d48 <- eln$d48[eln$group == "none"][1]
if (is.na(none_d48)) none_d48 <- 0.323

extrap <- annotate("rect", xmin = 60, xmax = 90, ymin = -Inf, ymax = Inf,
                   fill = "grey70", alpha = .25)
extrap_txt <- annotate("text", 75, -0.155, hjust = 0.5, vjust = 0, size = 2.2,
                       colour = "grey25", label = "no data >60 mmol\n(extrapolation)")
p3a <- ggplot(rcs %>% filter(outcome == "12-h change"), aes(dose, fit)) +
  extrap + extrap_txt +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = "#D55E00", alpha = .15) +
  geom_line(colour = "#D55E00", linewidth = .7) +
  geom_hline(yintercept = 0.035, colour = "grey60", lty = 2, linewidth = .4) +
  annotate("text", 58, 0.05, hjust = 1, size = 2.4, colour = "grey30", label = "zero-dose ~0.035*") +
  annotate("text", 30, 0.215, size = 2.6, colour = "#D55E00", label = "peak ~30 mmol: +0.21") +
  coord_cartesian(ylim = c(-0.18, 0.26), xlim = c(0, 90)) +
  labs(x = "Cumulative IV phosphate in first 24 h (mmol)", y = "Adjusted change (mmol/L)",
       title = "A. Dose-response at 12 h (RCS, nonlinearity p<0.001)") + tt

p3b <- ggplot(rcs %>% filter(outcome == "48-h change"), aes(dose, fit)) +
  extrap + extrap_txt +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = "#0072B2", alpha = .15) +
  geom_line(colour = "#0072B2", linewidth = .7) +
  coord_cartesian(ylim = c(-0.18, 0.26), xlim = c(0, 90)) +
  labs(x = "Cumulative IV phosphate in first 24 h (mmol)", y = "Adjusted change (mmol/L)",
       title = "B. Dose-response at 48 h (flat; p=0.18)") + tt

td2 <- td %>% mutate(timing = factor(timing, levels = c("early<=10h","late10-48h"),
                                     labels = c("Early (<=10 h)","Late (10-48 h)")),
                     dose_grp = factor(dose_grp, levels = c("dose<30","dose>=30"),
                                       labels = c("<30 mmol",">=30 mmol")))
p3c <- ggplot(td2, aes(timing, d48, fill = dose_grp)) +
  geom_hline(yintercept = none_d48, colour = "#0072B2", lty = 2, linewidth = .5) +
  geom_col(position = position_dodge(.8), width = .7) +
  geom_text(aes(label = sprintf("%.2f", d48)), position = position_dodge(.8),
            vjust = -0.5, size = 2.6) +
  annotate("text", 1.5, 0.352, colour = "#0072B2", size = 2.4,
           label = sprintf("untreated spontaneous recovery %.3f", none_d48), hjust = 0.5) +
  coord_cartesian(ylim = c(0, 0.40)) +
  labs(x = NULL, y = "48-h change (mmol/L)", fill = "24-h dose",
       title = "C. Timing x dose: 48-h change") + tt + scale_fill_manual(values = c("#E69F00","#D55E00"))

p3d <- ggplot(td2, aes(timing, p48, fill = dose_grp)) +
  geom_col(position = position_dodge(.8), width = .7) +
  geom_text(aes(label = sprintf("%.2f", p48)), position = position_dodge(.8),
            vjust = -0.5, size = 2.6) +
  coord_cartesian(ylim = c(0, 1.0)) +
  labs(x = NULL, y = "Proportion reaching >=0.80 by 48 h", fill = "24-h dose",
       title = "D. Timing x dose: target attainment") + tt +
  scale_fill_manual(values = c("#E69F00","#D55E00"))

cap3 <- "*dose 0 = 1,106 late-treated stays (infusion >24 h), not untreated patients; dose winsorized at p99 = 60 mmol; shaded region has no data support."
savefig((p3a | p3b) / (p3c | p3d) + plot_layout(guides = "collect") +
          plot_annotation(caption = cap3, theme = theme(
            legend.position = "bottom", legend.direction = "horizontal",
            plot.caption = element_text(size = base - 2.5, face = "italic",
                                        colour = "grey30", hjust = 0),
            plot.margin = margin(t = 2, r = 6, b = 8, l = 6))),
        "FigS_PC3_dose", 190, 166)

# =====================================================================
# S-PC4 : negative controls (sodium trajectory contrast + pre-T0 trend)
# =====================================================================
nc1 <- rd("pc_nc1_sodium.csv")
ph_diff <- rd("pc_gee_contrasts.csv") %>% filter(comparison == "C1_IV_vs_none") %>%
  transmute(time, diff = weighted_diff, lcl, ucl)

soda <- nc1 %>% transmute(time, diff, lcl, ucl, series = "Sodium (negative control)")
phos <- ph_diff %>% transmute(time, diff, lcl, ucl, series = "Phosphate (mechanistic contrast)")
sercol <- c("Sodium (negative control)" = "#0072B2",
            "Phosphate (mechanistic contrast)" = "#D55E00")
p4a <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = .4) +
  geom_line(data = phos, aes(time, diff, colour = series), linewidth = .7) +
  geom_pointrange(data = phos, aes(time, diff, ymin = lcl, ymax = ucl, colour = series),
                  size = .35, linewidth = .5) +
  geom_pointrange(data = soda, aes(time, diff, ymin = lcl, ymax = ucl, colour = series),
                  shape = 15, size = .4, linewidth = .5) +
  annotate("text", 6, 1.08, hjust = 0, vjust = 1, size = 2.4, colour = "#0072B2",
           label = "sodium gap +0.29 to +0.89 mmol/L:\nno peak-trough pharmacologic shape; offset > early phosphate signal (+0.054)") +
  scale_colour_manual(values = sercol) +
  scale_x_continuous(breaks = c(6,12,24,48,72)) +
  coord_cartesian(ylim = c(-0.25, 1.22), clip = "off") +
  labs(x = "Hours after T0", y = "Weighted IV - no-replacement difference (mmol/L)",
       title = "A. Negative-control outcome (sodium) vs phosphate", colour = NULL) + tt +
  theme(plot.margin = margin(t = 4, r = 12, b = 4, l = 10))

nc2 <- rd("pc_nc2_pretrend.csv")
pre <- data.frame(
  grp = rep(c("IV","none"), each = 2),
  val = c(nc2$A_IV[nc2$metric=="level_pre24"], nc2$A_IV[nc2$metric=="level_pre6"],
          nc2$D_none[nc2$metric=="level_pre24"], nc2$D_none[nc2$metric=="level_pre6"]),
  x = rep(c(1, 2), 2))
p4b <- ggplot(pre, aes(x, val, colour = grp)) +
  geom_line(linewidth = .9) + geom_point(size = 2.2) +
  annotate("text", 1.5, 1.095, size = 2.5, colour = GL["IV"],
           label = "IV change -0.216") +
  annotate("text", 1.5, 0.895, size = 2.5, colour = GL["none"],
           label = "none change -0.179") +
  scale_colour_manual(values = GL, labels = c(IV = "IV replacement", none = "No replacement")) +
  scale_x_continuous(breaks = c(1,2), labels = c("pre 24 h","pre 6 h")) +
  coord_cartesian(ylim = c(0.85, 1.2), xlim = c(0.8, 2.2)) +
  labs(x = "Pre-treatment window", y = "Serum phosphate (mmol/L)", colour = NULL,
       title = "B. Pre-treatment trend: group x time interaction p=0.62 (NS)",
       caption = "Point estimates suggest steeper IV decline but slope-difference CI crosses zero (-0.0088 to +0.0023); n = 858/505") + tt +
  theme(plot.caption = element_text(size = base - 2.5, face = "italic", colour = "grey30"))

savefig(p4a / p4b + plot_layout(guides = "collect") +
          plot_annotation(theme = theme(
            legend.position = "bottom",
            plot.margin = margin(t = 2, r = 6, b = 6, l = 8))),
        "FigS_PC4_negative_controls", 170, 164)

cat("Figures written to", FIG, "\n"); print(list.files(FIG))
