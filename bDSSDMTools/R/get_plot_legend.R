#' Extract Legend Grob from ggplot Object
#'
#' A helper function to extract the legend graphical object (grob) from a
#' ggplot object. This is useful for arranging plots and legends separately,
#' for example using 'cowplot' or 'gridExtra'.
#'
#' @param plot A ggplot object or a pre-computed grob.
#' @param legend Character string, optional. If the plot has multiple legends
#'   (e.g., from `color` and `size`), specify the name associated with the
#'   desired legend (e.g., "colour", "size"). Usually NULL is sufficient if
#'   there is only one main legend guide box.
#'
#' @return A grid grob (`gtable`) object representing the legend, or NULL if
#'   no legend guide-box is found or the identified legend is empty.
#'
#' @keywords internal
#' @importFrom ggplot2 is.ggplot ggplotGrob
#' @importFrom grid is.grob
#'
get_plot_legend <- function(plot, legend = NULL) {
  # Check if input is ggplot or grob
  if (ggplot2::is.ggplot(plot)) {
    gt <- ggplot2::ggplotGrob(plot)
  } else {
    if (grid::is.grob(plot)) {
      gt <- plot
    } else {
      stop("Plot object is neither a ggplot nor a grob.")
    }
  }

  # Construct pattern to find legend grob name
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }

  # Find grobs matching the pattern
  indices <- grep(pattern, gt$layout$name)

  # Filter out empty grobs (zeroGrobs)
  if (length(indices) > 0) {
    not_empty <- !vapply(
      gt$grobs[indices],
      inherits, what = "zeroGrob", # Check class using base::inherits
      FUN.VALUE = logical(1)
    )
    indices <- indices[not_empty]
  }

  # Return the first non-empty matching grob, or NULL
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }

  return(NULL)
}

# Note: Renamed function from get_legend2 to get_plot_legend for clarity.
# If you keep the old name, adjust the file name and roxygen comments accordingly.
