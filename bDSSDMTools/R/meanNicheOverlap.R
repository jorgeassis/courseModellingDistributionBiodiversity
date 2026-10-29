#' Mean Probabilistic Niche Overlap
#'
#' Estimates pairwise niche overlap with nicheROVER. For every species, niche
#' location and shape are drawn from a normal-inverse-Wishart posterior. Overlap
#' is the probability that an individual drawn from one species falls inside the
#' niche region of another (Swanson et al. 2015).
#'
#' The printed table is the posterior mean overlap, in percent, for a niche
#' region of size 95%. The returned array also contains the 99% niche region.
#' Overlap is not symmetric: rows are the species from which the individual is
#' drawn, and columns are the species whose niche region is tested. Species
#' follow the level order of `data$species` (alphabetical when that column is
#' character). `colors.sp` must follow that same order.
#'
#' `data` is the table returned by [prepareDataNicheOverlap()].
#'
#' @param data A data frame whose first column is named `species`. The remaining
#'   columns are numeric environmental variables. Each species needs more rows
#'   than there are environmental variables.
#' @param nsamples Number of posterior draws, and the number of overlap
#'   replicates. The Monte Carlo sample used inside each replicate is 10,000.
#' @param colors.sp Colours used in the overlap plot, one per species, in
#'   species order.
#'
#' @return Invisibly, a three-dimensional array of mean overlap in percent.
#'   The first two dimensions are the species pair and the third is the niche
#'   region size (`0.95`, then `0.99`). The 95% table is printed, and the
#'   posterior distribution of overlap for the 95% niche region is plotted.
#' @export
#' @importFrom nicheROVER niw.post overlap overlap.plot
#' @importFrom graphics par
#' @seealso [prepareDataNicheOverlap()], [plotNicheOverlap()], [nicheSize()]
#'
#' @examples
#' \dontrun{
#' overlapData <- prepareDataNicheOverlap(environment, rasterSPA, rasterSPB)
#' overlap <- meanNicheOverlap(overlapData, nsamples = 1000, colors.sp = c("#1b9e77", "#d95f02"))
#' overlap[, , 1]
#' }
meanNicheOverlap <- function(data, nsamples, colors.sp) {

  species.par <- .nichePosterior(data, nsamples)
  colors.sp <- .checkSpeciesColors(colors.sp, species.par)
  if (length(species.par) < 2) {
    stop("Error :: data must contain at least two species")
  }

  over.stat <- nicheROVER::overlap(
    species.par,
    nreps = nsamples,
    nprob = 10000,
    alpha = c(0.95, 0.99)
  )

  over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
  print(round(over.mean[, , 1], 2))

  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)
  nicheROVER::overlap.plot(
    over.stat[, , , 1],
    col = colors.sp,
    mean.cred.col = "turquoise",
    equal.axis = TRUE,
    xlab = "Overlap Probability (%) -- Niche Region Size: 95%"
  )

  invisible(over.mean)

}
