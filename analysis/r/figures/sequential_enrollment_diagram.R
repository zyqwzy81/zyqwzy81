# =============================================================================
# V28 Sequential Trial Emulation Enrollment Diagram
# Output: Fig1_SequentialTrial_Enrollment.pdf + .png (600 DPI)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(grid)
  library(patchwork)
})

# ---- Font configuration ----
FONT_FAMILY <- "DejaVu Sans"
# Register DejaVu Sans for PDF
grDevices::pdfFonts("DejaVu Sans" = grDevices::pdfFonts("sans")[[1]])

# ---- Paths ----
BASE_DIR  <- "/data/小米云盘/文档/论文/血磷"
OUT_DIR   <- file.path(BASE_DIR, "TTE", "v28", "output", "figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Colors ----
COL_TREATMENT  <- "#BC3C29"
COL_CONTROL    <- "#0072B5"
COL_THRESHOLD  <- "#E74C3C"
COL_WASHOUT    <- "#F39C12"
COL_PHOSPHATE  <- "#2E8B57"
COL_ICU        <- "#F0F0F0"

# ============================================================
# Build illustrative data
# ============================================================

phosphate_data <- data.table(
  day = c(0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.5, 7.0, 8.0, 9.0,
          10.0, 11.0, 12.5, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
          20.0, 21.0, 22.0, 23.0, 24.0, 24.5),
  phosphate = c(0.65, 0.72, 0.55, 0.48, 0.50, 0.62, 0.70, 0.58, 0.45,
                0.38, 0.42, 0.55, 0.65, 0.60, 0.50, 0.44, 0.48, 0.58,
                0.62, 0.55, 0.48, 0.42, 0.50, 0.58, 0.62, 0.55, 0.50)
)

threshold <- 0.80

# ============================================================
# Panel A: Phosphate trajectory + enrollment timeline
# ============================================================

pA <- ggplot() +
  # ICU stay background
  annotate("rect", xmin = 0, xmax = 25, ymin = -Inf, ymax = Inf,
           fill = COL_ICU, alpha = 0.25) +

  # Washout shaded areas
  annotate("rect", xmin = 2.0, xmax = 4.0, ymin = -Inf, ymax = Inf,
           fill = COL_WASHOUT, alpha = 0.10) +
  annotate("rect", xmin = 6.0, xmax = 10.0, ymin = -Inf, ymax = Inf,
           fill = COL_WASHOUT, alpha = 0.10) +

  # Threshold line
  geom_hline(yintercept = threshold, linetype = "dashed", color = COL_THRESHOLD, linewidth = 0.6) +

  # Phosphate line + points
  geom_line(data = phosphate_data, aes(x = day, y = phosphate),
            color = COL_PHOSPHATE, linewidth = 0.8) +
  geom_point(data = phosphate_data, aes(x = day, y = phosphate),
             color = COL_PHOSPHATE, size = 1.5) +

  # Enrollment #1 vertical line + marker
  annotate("segment", x = 0, xend = 0, y = 0.35, yend = 0.06, color = "black", linewidth = 0.6) +
  annotate("point", x = 0, y = 0.06, color = COL_TREATMENT, size = 3.5, shape = 18) +

  # Enrollment #2 vertical line + marker
  annotate("segment", x = 4, xend = 4, y = 0.35, yend = 0.06, color = "black", linewidth = 0.6) +
  annotate("point", x = 4, y = 0.06, color = COL_CONTROL, size = 3.5, shape = 18) +

  # Enrollment #3 vertical line + marker
  annotate("segment", x = 10, xend = 10, y = 0.35, yend = 0.06, color = "black", linewidth = 0.6) +
  annotate("point", x = 10, y = 0.06, color = COL_TREATMENT, size = 3.5, shape = 18) +

  # IV Phosphate markers (x marks)
  annotate("point", x = 0.2, y = 0.52, color = COL_TREATMENT, size = 2.5, shape = 4, stroke = 0.8) +
  annotate("point", x = 0.8, y = 0.52, color = COL_TREATMENT, size = 2.5, shape = 4, stroke = 0.8) +

  # ---- Labels with white background (above data, below threshold) ----
  # Washout label 1
  annotate("label", x = 3.0, y = 0.88, label = "Washout (>=48h)",
           size = 2.8, color = COL_WASHOUT, fill = alpha("white", 0.85), label.padding = unit(2, "pt")) +

  # Washout label 2
  annotate("label", x = 8.0, y = 0.88, label = "Washout (>=48h)",
           size = 2.8, color = COL_WASHOUT, fill = alpha("white", 0.85), label.padding = unit(2, "pt")) +

  # Enrollment #1 label
  annotate("label", x = 0, y = 0.30, label = "Enrollment #1\n(IV Phosphate)",
           size = 2.8, color = "black", fill = alpha("white", 0.85), label.padding = unit(2, "pt"), lineheight = 0.85) +

  # Enrollment #2 label
  annotate("label", x = 4, y = 0.30, label = "Enrollment #2\n(Control)",
           size = 2.8, color = "black", fill = alpha("white", 0.85), label.padding = unit(2, "pt"), lineheight = 0.85) +

  # Enrollment #3 label
  annotate("label", x = 10, y = 0.30, label = "Enrollment #3\n(IV Phosphate)",
           size = 2.8, color = "black", fill = alpha("white", 0.85), label.padding = unit(2, "pt"), lineheight = 0.85) +

  # IV P labels
  annotate("label", x = 0.2, y = 0.46, label = "IV P",
           size = 2.5, color = COL_TREATMENT, fill = alpha("white", 0.85), label.padding = unit(2, "pt")) +
  annotate("label", x = 0.8, y = 0.46, label = "IV P",
           size = 2.5, color = COL_TREATMENT, fill = alpha("white", 0.85), label.padding = unit(2, "pt")) +

  # Threshold annotation
  annotate("label", x = 0.3, y = threshold + 0.025, label = "Eligibility threshold: 0.80 mmol/L",
           hjust = 0, size = 3, color = COL_THRESHOLD, fontface = "italic",
           fill = alpha("white", 0.85), label.padding = unit(2, "pt")) +

  # ICU Stay annotation
  annotate("label", x = 24.5, y = 0.92, label = "ICU Stay",
           hjust = 1, size = 3.5, color = "grey50", fontface = "italic",
           fill = alpha("white", 0.85), label.padding = unit(2, "pt")) +

  # Axes - clip off to allow labels outside
  scale_x_continuous(breaks = seq(0, 25, by = 5), limits = c(-0.8, 26.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.02, 0.95), expand = c(0, 0),
                     breaks = seq(0.2, 0.8, by = 0.2),
                     labels = scales::label_number(accuracy = 0.01)) +

  coord_cartesian(clip = "off") +

  labs(x = "ICU Length of Stay (days)", y = "Serum Phosphate (mmol/L)",
       title = "A  Sequential Trial Enrollment") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0, margin = margin(b = 4), family = FONT_FAMILY),
    axis.title = element_text(size = 9, margin = margin(t = 6), family = FONT_FAMILY),
    axis.text = element_text(size = 8, family = FONT_FAMILY),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 14, 6, 12)
  )

# ============================================================
# Panel B: Trial period structure
# ============================================================

pB <- ggplot() +
  # Trial period bars
  annotate("segment", x = 0,  xend = 4,  y = 1, yend = 1, color = COL_TREATMENT, linewidth = 6, alpha = 0.75) +
  annotate("segment", x = 4,  xend = 10, y = 2, yend = 2, color = COL_CONTROL,   linewidth = 6, alpha = 0.75) +
  annotate("segment", x = 10, xend = 25, y = 3, yend = 3, color = COL_TREATMENT, linewidth = 6, alpha = 0.75) +

  # Labels on bars
  annotate("text", x = 2,  y = 1, label = "IV Phosphate",  size = 3.5, color = "white", fontface = "bold") +
  annotate("text", x = 7,  y = 2, label = "Conservative",  size = 3.5, color = "white", fontface = "bold") +
  annotate("text", x = 17.5, y = 3, label = "IV Phosphate", size = 3.5, color = "white", fontface = "bold") +

  # Follow-up arrow timeline
  annotate("segment", x = 0, xend = 25, y = 0.30, yend = 0.30,
           arrow = arrow(length = unit(0.06, "inches"), ends = "last", type = "closed"),
           linewidth = 0.5, color = "grey50") +

  annotate("label", x = 25.5, y = 0.50, label = "Follow-up",
           hjust = 0.5, size = 3, color = "grey40", fontface = "italic",
           fill = alpha("white", 0.85), label.padding = unit(2, "pt")) +

  # Enrollment markers on timeline
  annotate("point", x = 0,  y = 0.30, size = 2.5, color = "black") +
  annotate("point", x = 4,  y = 0.30, size = 2.5, color = "black") +
  annotate("point", x = 10, y = 0.30, size = 2.5, color = "black") +
  annotate("point", x = 25, y = 0.30, size = 2.5, color = "black") +

  # Vertical connector lines
  annotate("segment", x = 0,  xend = 0,  y = 0.30, yend = 0.08, linewidth = 0.4, color = "grey50") +
  annotate("segment", x = 4,  xend = 4,  y = 0.30, yend = 0.08, linewidth = 0.4, color = "grey50") +
  annotate("segment", x = 10, xend = 10, y = 0.30, yend = 0.08, linewidth = 0.4, color = "grey50") +
  annotate("segment", x = 25, xend = 25, y = 0.30, yend = 0.08, linewidth = 0.4, color = "grey50") +

  # Labels below timeline
  annotate("label", x = 0,  y = 0.04, label = "Enrollment #1",
           size = 3, color = "black", fill = alpha("white", 0.85), label.padding = unit(2, "pt")) +
  annotate("label", x = 4,  y = 0.04, label = "Enrollment #2",
           size = 3, color = "black", fill = alpha("white", 0.85), label.padding = unit(2, "pt")) +
  annotate("label", x = 10, y = 0.04, label = "Enrollment #3",
           size = 3, color = "black", fill = alpha("white", 0.85), label.padding = unit(2, "pt")) +
  annotate("label", x = 25, y = 0.02, label = "ICU Discharge\nor Day 28",
           size = 3, color = "black", fill = alpha("white", 0.85), label.padding = unit(2, "pt"), lineheight = 0.85) +

  # Y-axis labels
  annotate("text", x = -0.3, y = 1, label = "Period 1", size = 3.5, color = "black", hjust = 1) +
  annotate("text", x = -0.3, y = 2, label = "Period 2", size = 3.5, color = "black", hjust = 1) +
  annotate("text", x = -0.3, y = 3, label = "Period 3", size = 3.5, color = "black", hjust = 1) +

  # Axes
  scale_x_continuous(breaks = seq(0, 25, by = 5), limits = c(-0.8, 26.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.0, 3.6), expand = c(0, 0)) +

  coord_cartesian(clip = "off") +

  labs(x = "ICU Length of Stay (days)", y = "",
       title = "B  Trial Period Structure") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0, margin = margin(b = 4), family = FONT_FAMILY),
    axis.title.x = element_text(size = 9, margin = margin(t = 6), family = FONT_FAMILY),
    axis.text.x = element_text(size = 8, family = FONT_FAMILY),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(4, 14, 14, 12)
  )

# ============================================================
# Combine panels
# ============================================================

fig <- pA / pB +
  plot_layout(heights = c(1.6, 0.8)) +
  plot_annotation(
    title = "Sequential Target Trial Emulation: Enrollment and Trial Period Structure",
    caption = "Panel A: A patient's serum phosphate trajectory during a 25-day ICU stay. Each time phosphate drops below 0.80 mmol/L and >= 48 hours\nhave elapsed since prior enrollment, a new trial period begins. Panel B: Each period assigns the patient to IV phosphate supplementation or\nconservative management, lasting until the next enrollment, ICU discharge, or Day 28.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 8), family = FONT_FAMILY),
      plot.caption = element_text(size = 7.5, color = "grey40", hjust = 0.5, lineheight = 1.3, family = FONT_FAMILY)
    )
  )

# ============================================================
# Output
# ============================================================

# PDF (ggsave for consistent rendering with showtext)
pdf_path <- file.path(OUT_DIR, "Fig1_SequentialTrial_Enrollment.pdf")
ggsave(pdf_path, fig, width = 8.5, height = 6.5, device = cairo_pdf, bg = "white")
cat(sprintf("PDF: %s (%.0f KB)\n", pdf_path, file.size(pdf_path) / 1024))

# PNG (600 DPI)
png_path <- file.path(OUT_DIR, "Fig1_SequentialTrial_Enrollment.png")
ggsave(png_path, fig, width = 8.5, height = 6.5, dpi = 600, device = "png", bg = "white")
cat(sprintf("PNG: %s (%.0f KB)\n", png_path, file.size(png_path) / 1024))

cat("\nDone!\n")