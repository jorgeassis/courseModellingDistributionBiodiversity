#' Plot Posterior Niche Ellipses
#'
#' Draws the two-dimensional niche regions of each species with
#' `nicheROVER::niche.plot`. Ellipses show the 95% niche region for posterior
#' draws of the niche parameters. At least two environmental variables are
#' required, because each panel is a pair of variables.
#'
#' Species and colours follow the level order of `data$species` (alphabetical
#' when that column is character). `data` is the table returned by
#' [prepareDataNicheOverlap()].
#'
#' @param data A data frame whose first column is named `species`. The remaining
#'   columns are numeric environmental variables. At least two environmental
#'   columns are required.
#' @param nsamples Number of posterior draws. Each draw contributes one ellipse.
#' @param colors Colours, one per species, in species order.
#'
#' @return Called for its plot. Returns the value of `nicheROVER::niche.plot`,
#'   which is `NULL`.
#' @export
#' @importFrom nicheROVER niw.post niche.plot
#' @importFrom graphics par
#' @seealso [prepareDataNicheOverlap()], [meanNicheOverlap()], [nicheSize()]
#'
#' @examples
#' \dontrun{
#' overlapData <- prepareDataNicheOverlap(environment, rasterSPA, rasterSPB)
#' plotNicheOverlap(overlapData, nsamples = 20, colors = c("#1b9e77", "#d95f02"))
#' }
plotNicheOverlap <- function(data, nsamples, colors) {

  species.par <- .nichePosterior(data, nsamples)
  if (ncol(data) < 3) {
    stop("Error :: plotNicheOverlap needs at least two environmental variables")
  }
  colors <- .checkSpeciesColors(colors, species.par)

  species.par.data <- tapply(seq_len(nrow(data)), .speciesLabels(data), function(ii) {
    data[ii, -1, drop = FALSE]
  })

  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)
  nicheROVER::niche.plot(
    niche.par = species.par,
    niche.data = species.par.data,
    ndens = 1000,
    pfrac = 0,
    col = colors,
    xlab = expression("Niche overlap")
  )

}
