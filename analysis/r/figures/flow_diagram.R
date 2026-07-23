# ============================================================================
# V26.9 CONSORT-style Participant Flow Diagram (Fig S2, Fig 1)
# DAG14 版
# ============================================================================

library(grid)
library(gridExtra)
library(ggplot2)

OUTPUT_DIR <- "/data/小米云盘/文档/论文/血磷/TTE/v26/output/v26_9/figures"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

draw_flow_box <- function(label, x, y, width = 0.3, height = 0.08, 
                          font_size = 10, fill = "white") {
  grid.roundrect(
    x = x, y = y, 
    width = width, height = height,
    gp = gpar(fill = fill, col = "black", lwd = 1.5)
  )
  grid.text(label, x = x, y = y, 
            gp = gpar(fontsize = font_size, fontfamily = "sans"))
}

draw_arrow <- function(x0, y0, x1, y1, lwd = 1) {
  grid.segments(x0, y0, x1, y1, gp = gpar(lwd = lwd))
  dx <- x1 - x0
  dy <- y1 - y0
  len <- sqrt(dx^2 + dy^2)
  ux <- dx / len
  uy <- dy / len
  grid.polygon(
    x = c(x1, x1 - 0.02 * ux - 0.01 * uy, x1 - 0.02 * ux + 0.01 * uy),
    y = c(y1, y1 - 0.02 * uy - 0.01 * ux, y1 - 0.02 * uy + 0.01 * ux),
    gp = gpar(fill = "black", col = "black")
  )
}

render_flow <- function(file_path, is_pdf = TRUE) {
  if (is_pdf) {
    pdf(file = file_path, width = 10, height = 14)
  } else {
    png(file_path, width = 1000, height = 1400, res = 150)
  }
  
  pushViewport(viewport(layout = grid.layout(1, 1)))
  vp <- viewport(x = 0.5, y = 0.5, width = 0.9, height = 0.95)
  pushViewport(vp)
  
  col_main <- "#2C3E50"
  col_excl <- "#E74C3C"
  col_step <- "#3498DB"
  col_result <- "#27AE60"
  
  grid.text("Figure S2. Participant Flow Diagram",
            x = 0.5, y = 0.97, 
            gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "sans"))
  
  # Row 1: Screening
  y_start <- 0.90
  draw_flow_box("MIMIC-IV v3.1 (2008-2022)\nICU admissions screened", 
                0.5, y_start, 0.40, 0.06, 11, col_main)
  
  # Row 2: Exclusions
  y_excl <- 0.78
  draw_arrow(0.5, y_start - 0.035, 0.5, y_excl + 0.035)
  draw_flow_box("Excluded:\n• Age < 18 years\n• ICU stay < 24 hours\n• Serum phosphate > 0.80 mmol/L\n• ESRD or chronic dialysis", 
                0.5, y_excl, 0.40, 0.10, 9, col_excl)
  
  # Row 3: Eligible
  y_elig <- 0.62
  draw_arrow(0.5, y_excl - 0.055, 0.5, y_elig + 0.035)
  draw_flow_box("Eligible patients\nN = 22,685", 
                0.5, y_elig, 0.35, 0.05, 12, col_main)
  
  # Row 4: Sequential trials
  y_seq <- 0.52
  draw_arrow(0.5, y_elig - 0.03, 0.5, y_seq + 0.03)
  draw_flow_box("Sequential trial periods\nN = 32,330\n(Median 1 enrollment, IQR 1-1, range 1-8)", 
                0.5, y_seq, 0.45, 0.06, 11, col_step)
  
  # Row 5: Treatment groups
  y_tx <- 0.40
  draw_arrow(0.5, y_seq - 0.035, 0.35, y_tx + 0.035)
  draw_arrow(0.5, y_seq - 0.035, 0.65, y_tx + 0.035)
  draw_flow_box("IV Phosphate\nn = 5,147 (15.9%)", 
                0.32, y_tx, 0.28, 0.05, 11, col_step)
  draw_flow_box("Conservative Management\nn = 27,183 (84.1%)", 
                0.68, y_tx, 0.28, 0.05, 11, col_step)
  
  # Row 6: Follow-up
  y_fu <- 0.30
  draw_arrow(0.32, y_tx - 0.03, 0.32, y_fu + 0.03)
  draw_arrow(0.68, y_tx - 0.03, 0.68, y_fu + 0.03)
  draw_flow_box("28-day follow-up\n(All trial periods)", 
                0.5, y_fu, 0.40, 0.04, 10, col_step)
  
  # Row 7: Analysis (DAG14 — V26.9 update)
  y_an <- 0.22
  draw_arrow(0.5, y_fu - 0.025, 0.5, y_an + 0.025)
  draw_flow_box("Primary analysis: SEQTaRget ITT (bootstrap = 200, DAG14)\nTotal events: 5,610 (17.4%)", 
                0.5, y_an, 0.55, 0.05, 10, col_result)
  
  # eICU validation
  y_eicu <- 0.12
  draw_arrow(0.5, y_an - 0.03, 0.5, y_eicu + 0.03)
  draw_flow_box("External validation: eICU-CRD v2.0\nn = 34,882 trial periods, 5,508 events\nSEQTaRget ITT (bootstrap = 100, 11 shared covariates)",
                0.5, y_eicu, 0.55, 0.06, 9, "#8E44AD")
  
  popViewport()
  popViewport()
  dev.off()
  cat(sprintf("  ✓ %s\n", basename(file_path)))
}

# Render both PDF and PNG
cat("Updating Fig S2 Flow Diagram (DAG14)...\n")
render_flow(file.path(OUTPUT_DIR, "FigS2_Flow_Diagram.pdf"), is_pdf = TRUE)
render_flow(file.path(OUTPUT_DIR, "FigS2_Flow_Diagram.png"), is_pdf = FALSE)
cat("Done!\n")
