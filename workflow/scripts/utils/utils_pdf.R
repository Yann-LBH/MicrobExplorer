# ==============================================================================
# PROJECT   : MicrobExplorer
# SCRIPT    : utils_pdf.R
# PURPOSE   : Centralize pdf generation
# AUTHOR    : Yann Le Bihan
# DATE      : 2026-09-03
# LINK      : https://github.com/Yann-LBH/MicrobExplorer
# ==============================================================================

# Graphics utilities supporting configuration-driven page dimensions
with_pdf <- function(OUT_PDF, PARAM_PDF_SIZE, code_block) {
  pdf(OUT_PDF, width = PARAM_PDF_SIZE[1], height = PARAM_PDF_SIZE[2])
  
  # Guarantee PDF device cleanup regardless of success or runtime errors
  on.exit({
    if (names(dev.cur()) != "null device") {
      dev.off()
    }
  })
  
  # Execute the plotting logic enclosed in the code block
  force(code_block)
}

render_page <- function(plot_obj) {
  # ggplot compatibility
  if (inherits(plot_obj, "ggplot")) {
    print(plot_obj)
  # ComplexHeatmap compatibility
  } else if (inherits(plot_obj, c("Heatmap", "HeatmapList"))) {
    ComplexHeatmap::draw(plot_obj)
  # baseR compatibility
  } else if (inherits(plot_obj, "recordedplot")) {
    replayPlot(plot_obj)
  } else if (is.function(plot_obj)) {
    plot_obj()
  } else {
    stop("Unsupported plot format provided to render_page().")
  }
}

render_fallback <- function(message_text) {
  grid::grid.newpage()
  grid::grid.text(
    label = message_text,
    gp = grid::gpar(fontsize = 14, col = "red")
  )
}