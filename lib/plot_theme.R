# Shared figure styling: colour constants, band palette, theme, and the
# ggsave wrapper. Display-layer only.
#
# Submission typography (Springer Nature artwork guide linked from the BJC
# Guide to Authors, checked 2026-08-28): sans-serif text (Helvetica/Arial),
# 5-7 pt at final size; figures supplied at final width (BJC: 85 mm single /
# 175 mm double column); line art as vector (PDF/EPS). Figures are therefore
# drawn at final width with base_size 7 (axis/strip/legend text 5.6-7 pt),
# in-panel labels 5.5 pt, and no in-figure title/subtitle (the legend carries
# them). Researcher Go 2026-08-28.

# Significance colours (gene-level figures). Legend labels are built by the
# caller from FDR_CUT so the label always states the actual threshold.
COL_UP <- "#eb6834" # up in High
COL_DOWN <- "#2a78d6" # up in Sporadic
COL_NS <- "grey75"

# Exposure-band palette: colorblind-safe, ordered by assigned share
# (blue -> aqua -> yellow -> orange).
PAL_BANDS <- c(
  R_Sporadic = "#2a78d6", R_Low = "#1baf7a", R_Mid = "#eda100",
  R_High = "#eb6834"
)

FONT_FAMILY <- "Liberation Sans" # Arial-metric sans; installed in the canonical image
LABEL_PT <- 5.5 # in-panel text (gene labels, annotations), pt at final size
label_size <- function(pt = LABEL_PT) pt / ggplot2::.pt # geom_text size is in mm

theme_thyr <- function(base_size = 7, legend_position = "top") {
  ggplot2::theme_bw(base_size = base_size, base_family = FONT_FAMILY) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = legend_position,
      legend.key.size = grid::unit(base_size, "pt"),
      panel.spacing = grid::unit(3, "mm"), # keeps facet-edge tick labels apart
      plot.margin = ggplot2::margin(2, 4, 2, 2)
    )
}

# Figures land in output/figures/ (reorg plan v2 s2.6: figure scripts read
# processed/ only and write output/figures/ only). width/height are in mm at
# final display size; three copies are written from the same object:
#   .png 300 dpi  -- embedded in the review/submission docx (Word -> PDF keeps 300 dpi)
#   .tif 600 dpi  -- bitmap submission copy (BJC GTA: >= 300 dpi)
#   .pdf vector   -- line-art submission copy (artwork guide), fonts embedded
save_figure <- function(plot, filename, width, height) {
  out_dir <- file.path(paths$output, "figures")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  out_png <- file.path(out_dir, filename)
  ggplot2::ggsave(out_png, plot,
    width = width, height = height, units = "mm", dpi = 300, type = "cairo"
  )
  message("Saved: ", out_png)
  out_tif <- file.path(out_dir, sub("\\.png$", ".tif", filename))
  ggplot2::ggsave(out_tif, plot,
    width = width, height = height, units = "mm", dpi = 600, device = "tiff",
    type = "cairo", compression = "lzw"
  )
  message("Saved: ", out_tif)
  out_pdf <- file.path(out_dir, sub("\\.png$", ".pdf", filename))
  ggplot2::ggsave(out_pdf, plot,
    width = width, height = height, units = "mm", device = grDevices::cairo_pdf
  )
  message("Saved: ", out_pdf)
  invisible(out_png)
}
