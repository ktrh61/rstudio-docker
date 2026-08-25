# Shared figure styling: colour constants, band palette, theme, and the
# ggsave wrapper. Display-layer only.

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

theme_thyr <- function(base_size = 11, subtitle_size = 9,
                       legend_position = "top") {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(
        size = subtitle_size, colour = "grey30"
      ),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = legend_position
    )
}

# Figures land in output/figures/ (reorg plan v2 s2.6: figure scripts read
# processed/ only and write output/figures/ only).
save_figure <- function(plot, filename, width, height) {
  out_dir <- file.path(paths$output, "figures")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  out_png <- file.path(out_dir, filename)
  ggplot2::ggsave(out_png, plot,
    width = width, height = height, dpi = 160, type = "cairo"
  )
  message("Saved: ", out_png)
  # Submission copy (BJC artwork: minimum 300 dpi) -- 600 dpi LZW TIFF,
  # same physical size as the PNG preview.
  out_tif <- file.path(out_dir, sub("\\.png$", ".tif", filename))
  ggplot2::ggsave(out_tif, plot,
    width = width, height = height, dpi = 600, device = "tiff",
    type = "cairo", compression = "lzw"
  )
  message("Saved: ", out_tif)
  invisible(out_png)
}
