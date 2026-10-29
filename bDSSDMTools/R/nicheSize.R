#' Posterior Niche Size
#'
#' Estimates the size of each species' 95% niche region from the posterior
#' draws produced by nicheROVER. Niche size is the hypervolume of that region.
#' The same posterior draws are also plotted with `nicheROVER::niche.plot`, so
#' at least two environmental variables are required.
#'
#' Species and colours follow the level order of `data$species` (alphabetical
#' when that column is character). `data` is the table returned by
#' [prepareDataNicheOverlap()].
#'
#' @param data A data frame whose first column is named `species`. The remaining
#'   columns are numeric environmental variables. At least two environmental
#'   columns are required.
#' @param nsamples Number of posterior draws.
#' @param colors Colours used in the niche plot, one per species, in species
#'   order.
#'
#' @return A matrix with one column per species. Row `est` is the mean niche
#'   size across posterior draws. Row `se` is the standard deviation of those
#'   draws. A niche plot is also drawn.
#' @export
#' @importFrom nicheROVER niw.post niche.plot niche.size
#' @importFrom graphics par
#' @importFrom stats sd
#' @seealso [prepareDataNicheOverlap()], [meanNicheOverlap()], [plotNicheOverlap()]
#'
#' @examples
#' \dontrun{
#' overlapData <- prepareDataNicheOverlap(environment, rasterSPA, rasterSPB)
#' nicheSize(overlapData, nsamples = 1000, colors = c("#1b9e77", "#d95f02"))
#' }
nicheSize <- function(data, nsamples, colors) {

  species.par <- .nichePosterior(data, nsamples)
  if (ncol(data) < 3) {
    stop("Error :: nicheSize needs at least two environmental variables")
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

  niche.size <- sapply(species.par, function(spec) {
    apply(spec$Sigma, 3, nicheROVER::niche.size, alpha = 0.95)
  })
  rbind(est = colMeans(niche.size), se = apply(niche.size, 2, stats::sd))

}
